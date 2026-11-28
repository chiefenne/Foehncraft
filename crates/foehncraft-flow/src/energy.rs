//! The energy equation of the fluid (stage 7, PR 7a1; D150, D151, D164):
//! the temperature T in the form weighted with ρ c_p,
//!
//! ρ c_p ∂T/∂t + ∇·(Ċ T) − T ∇·Ċ − ∇·(k ∇T) = q''' + Φ,   Ċ_f = ṁ_f c_p,f,
//!
//! solved on the generic scalar transport of D116 ([`crate::scalar`])
//! after the pressure correction of every outer iteration. The convection
//! term is that of the bounded form with the heat capacity flux Ċ instead
//! of the mass flux, so that it equals ρ c_p u·∇T also where c_p varies;
//! the diffusion is the conduction with k at the faces (the series
//! conductance between materials in 7d1 is then a face value as well).
//! c_p(T), k(T) and μ(T) of the material are refreshed per cell and face
//! from the current temperature before each energy solve (μ: the face
//! viscosity of the next momentum solve); ρ stays the one number of the
//! case in 7a1 (a density that follows T is 7c).
//!
//! Case file:
//!
//! ```toml
//! [fluid]
//! material = "air"              # needs specific_heat and conductivity
//!
//! [energy]
//! initial = 300.0               # K: number or expression of x, y, z
//! source = 0.0                  # q''' (W/m³): number or expression of x, y, z (t)
//! viscous_heating = true        # Φ = τ:∇u (default false)
//!
//! [boundary.inlet.energy]
//! temperature = 300.0           # K
//! [boundary.wall.energy]
//! heat_flux = 500.0             # W/m² into the fluid
//! # heat_flow = 2.0             # W into the fluid through the whole boundary
//! # heat_transfer_coefficient = 10.0   # W/(m² K), with
//! # ambient_temperature = 293.15       # K: q = h (T_amb − T_wall)
//! # temperature = 350.0
//!
//! [solver.energy]               # as [solver.scalar]
//! scheme = "linear-upwind"
//! ```
//!
//! The condition follows from the keys of the table (no `type` key, D156).
//! Boundaries without an `energy` table: walls and slip walls are
//! adiabatic (THM005 lists them), symmetry planes and wedges pass no heat
//! and take no table (THM003), outlets have a zero normal gradient (the
//! backflow carries the owner value), periodic pairs carry T across;
//! inlets need a `temperature` (THM011).

use crate::bc::BoundaryCondition;
use crate::diag::{self, Diagnostic};
use crate::expr::{self, Env, Scalar, Vars};
use crate::materials::{Conditions, Conductivity, Curve, MaterialSection, Property, curve, lookup};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The name of the temperature field (output, restart, log columns).
pub const FIELD: &str = "temperature";

/// `[energy]`: the energy equation of the fluid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EnergySection {
    /// Start temperature (K): a number or an expression of x, y, z (D53).
    pub initial: Scalar,
    /// Volumetric heat source q''' (W/m³): a number or an expression of
    /// x, y, z and, in a transient case, t (D53).
    #[serde(default)]
    pub source: Scalar,
    /// The viscous heating Φ = τ:∇u = μ_eff (∇u + ∇uᵀ):∇u as a source
    /// (default off; it matters where the Brinkman number μU²/(kΔT) is not
    /// small: oils, lubrication gaps, high-speed gas flows).
    #[serde(default)]
    pub viscous_heating: bool,
    /// Refresh c_p, k and μ from the temperature every outer iteration
    /// (default true). Off: the properties stay those at `[fluid]
    /// temperature`.
    #[serde(default = "yes")]
    pub temperature_dependent: bool,
    /// The exact temperature, for the verification (an expression of x,
    /// y, z and, in a transient case, t), as `exact` of a scalar.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exact: Option<Scalar>,
    /// Turbulent Prandtl number Pr_t of the `constant` model (D135); the
    /// turbulent conductivity is c_p μ_t/Pr_t.
    #[serde(default = "default_prt")]
    pub turbulent_prandtl: f64,
    /// Model of Pr_t: `constant` or `kays-crawford`.
    #[serde(default)]
    pub turbulent_prandtl_model: crate::scalar::PrandtlModel,
    /// Thermal wall function at the walls of a turbulent flow (D135).
    #[serde(default)]
    pub wall_function: crate::scalar::ThermalWallFunction,
    /// Heat sources of single regions (stage 7d2), in addition to
    /// `source`: `[[energy.sources]]` with a `region` and `power` (W over
    /// the region, spread by volume) or `power_density` (W/m³).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sources: Vec<RegionSource>,
    /// Thermal acceleration of the solid regions in a transient conjugate
    /// run (D213, [`crate::acceleration`]): ρ c_p of the solids scaled by
    /// a factor until their temperatures settle, then restored.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub acceleration: Option<Acceleration>,
    /// `constant` (default): one density, that of `[fluid] temperature`;
    /// `variable`: ρ(T, p₀) from the fluid's density form in every cell, the
    /// low-Mach variable density of stage 7c (D151).
    #[serde(
        default,
        skip_serializing_if = "crate::density::DensityModel::is_default"
    )]
    pub density: crate::density::DensityModel,
    /// The thermodynamic pressure p₀ of the variable density: `constant`
    /// (default; `[fluid] operating_pressure`, a domain open to an ambient)
    /// or `closed-volume` (p₀(t) from the fixed mass of a closed domain).
    #[serde(
        default,
        skip_serializing_if = "crate::density::ThermodynamicPressure::is_default"
    )]
    pub thermodynamic_pressure: crate::density::ThermodynamicPressure,
    /// The pressure work of the energy equation: `none` or `thermodynamic`
    /// (dp₀/dt per unit volume). Default: `thermodynamic` in a closed volume
    /// (otherwise a heated gas would warm with c_p instead of c_v), `none`
    /// in an open domain, where p₀ is constant.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pressure_work: Option<crate::density::PressureWork>,
    /// Cycle-averaged coupling of a periodic flow with its solids (D204,
    /// [`crate::cycle`]): the flow runs over cycles with the solids held,
    /// the cycle-mean interface heat transfer goes to a steady solve of the
    /// solids, and the two alternate until the wall temperatures settle.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cycle: Option<CycleCoupling>,
    /// The frozen-flow march of a transient conjugate run (D118,
    /// [`crate::frozen`]): once the flow has settled it is held, and the
    /// temperatures of fluid and solids march alone with large steps; the
    /// flow is refreshed by coupled steps when the temperatures it depends
    /// on have moved.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub frozen_flow: Option<FrozenFlow>,
}

/// `[energy.frozen_flow]`: the frozen-flow march of a transient conjugate
/// run (D118, §9f of the stage 7 specification). A steady cooling flow
/// settles in a few flow-through times, the solids it cools warm up over
/// minutes or hours: the run steps fully coupled with the `[time]` step
/// until the flow has settled, then holds the velocity, pressure, mass
/// flux and turbulence and marches only the energy equation of fluid and
/// solids (the same matrix, D158) with steps of `time_step`. The history
/// stays physical: the flow is the settled one of the current
/// temperatures. When the flow depends on the temperature (μ(T), the
/// Boussinesq buoyancy) it is refreshed by coupled steps, until it has
/// settled again, whenever a fluid temperature has moved by more than
/// `refresh_change` since the flow was held (or every `refresh_every`
/// frozen steps); a flow that does not depend on the temperature is held
/// once and exactly.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FrozenFlow {
    /// `false` runs the case fully coupled, as if the section were not
    /// there: the reference of the march (default `true`).
    #[serde(default = "default_enabled")]
    pub enabled: bool,
    /// The step of the energy march with the flow held (s), at least the
    /// `[time]` step of the flow.
    pub time_step: f64,
    /// When the flow is first held: `"settled"` (default: by the window
    /// test below) or a time (s).
    #[serde(default)]
    pub start: FrozenStart,
    /// Flow steps per window of the settling test (default 10).
    #[serde(default = "default_window_steps")]
    pub window_steps: u32,
    /// The flow has settled when the relative rms change of the velocity
    /// over the last window and the change still to come (extrapolated
    /// geometrically from the last two windows) are below this (default
    /// 1e-3); for a flow that follows the temperature, which drifts with
    /// it, the change of the last window alone.
    #[serde(default = "default_settle_tolerance")]
    pub settle_tolerance: f64,
    /// K: a flow that depends on the temperature is refreshed when a fluid
    /// temperature has moved by more than this since the flow was held
    /// (default 2).
    #[serde(default = "default_refresh_change")]
    pub refresh_change: f64,
    /// Frozen steps between two refreshes at most (default 0: only by
    /// `refresh_change`); counts also for a flow that does not depend on
    /// the temperature (a boundary value that changes slowly).
    #[serde(default)]
    pub refresh_every: u32,
    /// Solves of the energy equation per frozen step at most (default 20;
    /// a step stops when the largest change of a temperature between two
    /// solves is below 10⁻⁸ of the temperature span).
    #[serde(default = "default_frozen_iterations")]
    pub iterations: u32,
}

/// `[energy.frozen_flow] start`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum FrozenStart {
    /// The time (s) at which the flow is first held.
    Time(f64),
    /// `"settled"`.
    Word(String),
}

impl Default for FrozenStart {
    fn default() -> Self {
        FrozenStart::Word("settled".into())
    }
}

impl FrozenFlow {
    /// The start time, `None` for `"settled"`.
    pub fn start_time(&self) -> Option<f64> {
        match self.start {
            FrozenStart::Time(t) => Some(t),
            FrozenStart::Word(_) => None,
        }
    }
}

fn default_refresh_change() -> f64 {
    2.0
}

fn default_frozen_iterations() -> u32 {
    20
}

/// `[energy.cycle]`: cycle-averaged conjugate coupling of a periodic flow
/// (D204, §9d of the stage 7 specification). Fluid phase: the flow and the
/// fluid energy run `warmup_cycles` + `average_cycles` periods with the
/// solid temperatures held; per interface face the run accumulates the
/// conductance G (fluid side, contact and solid half cell in series) and
/// the heat flow, which give the Robin pair (Ḡ, T̄_ref) that reproduces the
/// cycle-mean heat flow of the face exactly. Solid phase: a steady solve of
/// the solids with these Robin conditions. The new solid temperatures
/// (under-relaxed) are held in the next fluid phase, which continues from
/// the last flow state. The run stops when the wall temperatures and the
/// interface heat flows settle.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CycleCoupling {
    /// `false` runs the case fully coupled, as if the section were not
    /// there: the reference of the loop (default `true`).
    #[serde(default = "default_enabled")]
    pub enabled: bool,
    /// The period of the flow (s); a whole number of time steps.
    pub period: f64,
    /// Cycles run before the averaging, after every update (default 1).
    #[serde(default = "default_one")]
    pub warmup_cycles: u32,
    /// Cycles averaged (default 1).
    #[serde(default = "default_one")]
    pub average_cycles: u32,
    /// Cycles of one fluid phase at most, when the flow is not yet periodic
    /// (default 6).
    #[serde(default = "default_max_cycles")]
    pub max_cycles: u32,
    /// The flow counts as periodic when the mean interface heat flow of
    /// the last cycle differs from the cycle before by less than this
    /// fraction (default 0.005).
    #[serde(default = "default_cycle_tolerance")]
    pub cycle_tolerance: f64,
    /// Updates of the solid temperatures at most (default 5).
    #[serde(default = "default_max_updates")]
    pub max_updates: u32,
    /// K: the loop has converged when the largest change of a wall
    /// temperature in an update is below this (default 0.5) ...
    #[serde(default = "default_cycle_tolerance_k")]
    pub tolerance: f64,
    /// ... and the interface heat flow changes by less than this fraction
    /// (default 0.005).
    #[serde(default = "default_cycle_tolerance")]
    pub heat_tolerance: f64,
    /// Under-relaxation of the solid temperatures in the first update, and
    /// in every update with `update = "relaxation"`; in (0, 1] (default 1).
    #[serde(default = "default_relaxation_one")]
    pub relaxation: f64,
    /// The solid phase: `"steady"` (default: the steady solids of the
    /// cycle-mean heat transfer, the end state) or `"transient"` (a march
    /// of fluid and solids under the cycle-mean flow with large steps over
    /// `march_cycles` periods per update: the warm-up history, D118).
    #[serde(default = "default_solid_phase")]
    pub solid: String,
    /// `solid = "transient"`: the periods of physical time per update, the
    /// fluid phase included (the solids are held during it and the march
    /// covers that time too); the run's clock jumps over the rest, by
    /// whole periods.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub march_cycles: Option<u32>,
    /// `solid = "transient"`: the steps of the march per update (default 4).
    #[serde(default = "default_march_steps")]
    pub march_steps: u32,
    /// How the solid temperatures are updated from the second update on:
    /// `"quasi-newton"` (default: a least-squares model of the loop's
    /// response from all earlier updates, IQN-ILS), `"aitken"` (a dynamic
    /// relaxation factor) or `"relaxation"` (the fixed `relaxation`).
    #[serde(default = "default_cycle_update")]
    pub update: String,
    /// The model of the solid phase: `"mean-flow"` (default: the steady
    /// energy equation of fluid and solids with the cycle-mean mass flux,
    /// corrected per face to the measured cycle-mean heat flow) or
    /// `"robin"` (the solids alone with the Robin pair (Ḡ, T̄_ref) of the
    /// first fluid cell).
    #[serde(default = "default_solid_model")]
    pub solid_model: String,
    /// The swing of the wall temperature over a cycle: `"none"` (default:
    /// the wall is held at its cycle-mean temperature in the fluid phases)
    /// or `"effusivity"` (the periodic solution of a semi-infinite solid
    /// under the measured harmonics of the face heat flux, T̂_w = q̂ /
    /// (e √(iω)), e = √(k ρ c_p): the next fluid phase sees it, and
    /// `cycle_faces.csv` gives its amplitude, minimum and maximum per face).
    #[serde(default = "default_skin")]
    pub skin: String,
    /// The harmonics of the period that the skin swing takes (default 10;
    /// at most half the time steps of a cycle).
    #[serde(default = "default_skin_harmonics")]
    pub skin_harmonics: u32,
    /// Fully coupled cycles after the loop has ended, started from its
    /// solid temperatures: the exact swing, also at thin walls and edges
    /// (default 0).
    #[serde(default)]
    pub resolved_cycles: u32,
}

fn default_enabled() -> bool {
    true
}
fn default_one() -> u32 {
    1
}
fn default_max_cycles() -> u32 {
    6
}
fn default_cycle_tolerance() -> f64 {
    0.005
}
fn default_max_updates() -> u32 {
    5
}
fn default_cycle_tolerance_k() -> f64 {
    0.5
}
fn default_relaxation_one() -> f64 {
    1.0
}
fn default_skin() -> String {
    "none".into()
}
fn default_skin_harmonics() -> u32 {
    10
}
fn default_solid_model() -> String {
    "mean-flow".into()
}
fn default_cycle_update() -> String {
    "quasi-newton".into()
}
fn default_solid_phase() -> String {
    "steady".into()
}
fn default_march_steps() -> u32 {
    4
}

/// `[energy.acceleration]`: thermal acceleration of the solids (D213). The
/// heat capacity ρ c_p of the solid regions is multiplied by
/// `inertia_scale` (0 < f ≤ 1), so that they reach their steady or
/// periodic-mean temperatures f times sooner; the steady temperatures do
/// not depend on ρ c_p of a solid. When the solids' mean temperatures have
/// settled (or at `restore` seconds), the true ρ c_p comes back, at once
/// or raised geometrically over `restore_ramp` windows. The history is not
/// physical before the restore; the result says from when it is.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Acceleration {
    /// The factor f on ρ c_p of the accelerated solids, 0 < f ≤ 1, or
    /// `"auto"`: f from the time constant of each region, measured after
    /// the first window, so that the scaled time constant is
    /// `target_windows` windows.
    pub inertia_scale: InertiaScale,
    /// The accelerated solid regions (default: every solid region).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub regions: Vec<String>,
    /// When the true ρ c_p comes back: `"settled"` (default) or a time (s)
    /// at which the restore starts.
    #[serde(default)]
    pub restore: Restore,
    /// The period of a periodic flow (s): a window is one period, and the
    /// restore is raised geometrically over `restore_ramp` periods. Without
    /// it a window is `window_steps` time steps.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub period: Option<f64>,
    /// Time steps per window without a `period` (default 10).
    #[serde(default = "default_window_steps")]
    pub window_steps: u32,
    /// Settled when, for every accelerated region, the change of its
    /// window-mean temperature still to come (extrapolated geometrically
    /// from the changes of the last two windows) and the change of the
    /// last window are below this fraction of its change since the start
    /// (default 1e-3).
    #[serde(default = "default_settle_tolerance")]
    pub settle_tolerance: f64,
    /// Windows over which f rises geometrically to 1 (default 3 with a
    /// `period`, 0 without: a jump is exact where the flow is not
    /// periodic; in a periodic flow a jump leaves the part of the swing
    /// that the accelerated solid had at that moment in its deep
    /// temperature, which then decays only with the true time constant).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub restore_ramp: Option<u32>,
    /// `"auto"`: the scaled time constant in windows (default 2 with a
    /// `period`, 1 without).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_windows: Option<f64>,
    /// Windows after the restore before the history counts as physical
    /// (default 2; reported, nothing is cut).
    #[serde(default = "default_settle_windows")]
    pub settle_windows: u32,
}

/// The factor of `[energy.acceleration] inertia_scale`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum InertiaScale {
    /// A fixed factor, 0 < f ≤ 1.
    Factor(f64),
    /// `"auto"`.
    Word(String),
}

/// `[energy.acceleration] restore`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum Restore {
    /// The time (s) at which the restore starts.
    Time(f64),
    /// `"settled"`.
    Word(String),
}

impl Default for Restore {
    fn default() -> Self {
        Restore::Word("settled".into())
    }
}

impl Acceleration {
    /// The restore ramp in windows.
    pub fn ramp(&self) -> u32 {
        self.restore_ramp
            .unwrap_or(if self.period.is_some() { 3 } else { 0 })
    }

    /// The target of `"auto"` in windows.
    pub fn target(&self) -> f64 {
        self.target_windows
            .unwrap_or(if self.period.is_some() { 2.0 } else { 1.0 })
    }

    /// The fixed factor, `None` for `"auto"`.
    pub fn factor(&self) -> Option<f64> {
        match self.inertia_scale {
            InertiaScale::Factor(f) => Some(f),
            InertiaScale::Word(_) => None,
        }
    }

    /// The restore time, `None` for `"settled"`.
    pub fn restore_time(&self) -> Option<f64> {
        match self.restore {
            Restore::Time(t) => Some(t),
            Restore::Word(_) => None,
        }
    }
}

fn default_window_steps() -> u32 {
    10
}

fn default_settle_tolerance() -> f64 {
    1e-3
}

fn default_settle_windows() -> u32 {
    2
}

/// Checks `[energy.cycle]` (CYC001).
fn validate_cycle(case: &crate::case::FlowCase, c: &CycleCoupling, accel: bool) -> Result<()> {
    let key = "energy.cycle";
    let Some(t) = case.time.as_ref() else {
        crate::need!(
            false,
            CYC_SETTING,
            "[{key}] in a steady case: the loop couples a periodic flow with its solids; a steady \
             conjugate case is solved monolithically and needs no loop"
        );
        unreachable!()
    };
    crate::need!(
        crate::cht::is_cht(case),
        CYC_SETTING,
        "[{key}] without solid regions: the loop alternates between the flow and the solids"
    );
    crate::need!(
        !accel,
        CYC_SETTING,
        "[{key}] with [energy.acceleration]: the steady solid phase of the loop is already the \
         limit of a vanishing ρ c_p of the solids; use one of the two"
    );
    crate::need!(
        c.period > 0.0 && c.period.is_finite(),
        CYC_SETTING,
        "[{key}] period must be > 0 s, is {}",
        c.period
    );
    crate::need!(
        !t.adaptive,
        CYC_SETTING,
        "[{key}] with an adaptive time step: the cycles need a fixed time step that divides the \
         period"
    );
    let n = c.period / t.time_step;
    crate::need!(
        n >= 2.0 && (n - n.round()).abs() <= 1e-6 * n,
        CYC_SETTING,
        "[{key}] period {} s is not a whole number (>= 2) of time steps of {} s",
        c.period,
        t.time_step
    );
    crate::need!(
        c.average_cycles >= 1,
        CYC_SETTING,
        "[{key}] average_cycles must be >= 1"
    );
    crate::need!(
        c.max_cycles >= (c.warmup_cycles + c.average_cycles).max(2),
        CYC_SETTING,
        "[{key}] max_cycles ({}) must be at least warmup_cycles + average_cycles and at least 2",
        c.max_cycles
    );
    crate::need!(
        c.max_updates >= 1,
        CYC_SETTING,
        "[{key}] max_updates must be >= 1"
    );
    for (name, v) in [
        ("cycle_tolerance", c.cycle_tolerance),
        ("heat_tolerance", c.heat_tolerance),
    ] {
        crate::need!(
            v > 0.0 && v < 1.0,
            CYC_SETTING,
            "[{key}] {name} must be in (0, 1), is {v}"
        );
    }
    crate::need!(
        c.tolerance > 0.0 && c.tolerance.is_finite(),
        CYC_SETTING,
        "[{key}] tolerance must be > 0 K, is {}",
        c.tolerance
    );
    crate::need!(
        c.relaxation > 0.0 && c.relaxation <= 1.0,
        CYC_SETTING,
        "[{key}] relaxation must be in (0, 1], is {}",
        c.relaxation
    );
    crate::need!(
        ["quasi-newton", "aitken", "relaxation"].contains(&c.update.as_str()),
        CYC_SETTING,
        "[{key}] update must be \"quasi-newton\", \"aitken\" or \"relaxation\", is \"{}\"",
        c.update
    );
    crate::need!(
        ["mean-flow", "robin"].contains(&c.solid_model.as_str()),
        CYC_SETTING,
        "[{key}] solid_model must be \"mean-flow\" or \"robin\", is \"{}\"",
        c.solid_model
    );
    crate::need!(
        ["none", "effusivity"].contains(&c.skin.as_str()),
        CYC_SETTING,
        "[{key}] skin must be \"none\" or \"effusivity\", is \"{}\"",
        c.skin
    );
    crate::need!(
        c.skin_harmonics >= 1,
        CYC_SETTING,
        "[{key}] skin_harmonics must be >= 1"
    );
    crate::need!(
        ["steady", "transient"].contains(&c.solid.as_str()),
        CYC_SETTING,
        "[{key}] solid must be \"steady\" or \"transient\", is \"{}\"",
        c.solid
    );
    if c.solid == "transient" {
        let phase = c.warmup_cycles + c.average_cycles;
        crate::need!(
            c.march_cycles.is_some_and(|m| m >= phase),
            CYC_SETTING,
            "[{key}] solid = \"transient\" needs march_cycles, the periods of physical time per \
             update (at least warmup_cycles + average_cycles = {phase}: the march covers the fluid \
             phase, in which the solids were held), is {:?}",
            c.march_cycles
        );
        crate::need!(
            c.march_steps >= 1,
            CYC_SETTING,
            "[{key}] march_steps must be >= 1"
        );
        crate::need!(
            c.relaxation == 1.0,
            CYC_SETTING,
            "[{key}] solid = \"transient\" with relaxation {}: the march gives the new solid \
             temperatures of a history, which an under-relaxation would slow down; leave it at 1",
            c.relaxation
        );
    } else {
        crate::need!(
            c.march_cycles.is_none(),
            CYC_SETTING,
            "[{key}] march_cycles belongs to solid = \"transient\""
        );
    }
    Ok(())
}

/// Checks `[energy.frozen_flow]` (FRZ001).
fn validate_frozen(case: &crate::case::FlowCase, f: &FrozenFlow, e: &EnergySection) -> Result<()> {
    let key = "energy.frozen_flow";
    let Some(t) = case.time.as_ref() else {
        crate::need!(
            false,
            FRZ_SETTING,
            "[{key}] in a steady case: a steady conjugate case is solved monolithically and has no \
             history to march"
        );
        unreachable!()
    };
    crate::need!(
        crate::cht::is_cht(case),
        FRZ_SETTING,
        "[{key}] without solid regions: the fluid temperatures follow the flow within a \
         flow-through time, so there is no slow part to march with large steps"
    );
    crate::need!(
        e.acceleration.is_none(),
        FRZ_SETTING,
        "[{key}] with [energy.acceleration]: the frozen-flow march keeps the true history, the \
         acceleration gives it up; use one of the two"
    );
    crate::need!(
        e.cycle.is_none(),
        FRZ_SETTING,
        "[{key}] with [energy.cycle]: a periodic flow does not settle; the loop's solid = \
         \"transient\" marches the solids under the cycle-mean flow"
    );
    crate::need!(
        e.density != crate::density::DensityModel::Variable
            && case.buoyancy.model != crate::buoyancy::Model::VariableDensity,
        FRZ_SETTING,
        "[{key}] with the variable density: the mass flux follows the temperature through the \
         density, which a held flux cannot; the frozen flow needs a constant density (with the \
         Boussinesq buoyancy if needed)"
    );
    crate::need!(
        !t.adaptive,
        FRZ_SETTING,
        "[{key}] with an adaptive time step: the coupled steps take the fixed [time] step, the \
         march its own"
    );
    crate::need!(
        f.time_step.is_finite() && f.time_step >= t.time_step,
        FRZ_SETTING,
        "[{key}] time_step must be at least the [time] time_step ({} s), is {}",
        t.time_step,
        f.time_step
    );
    if let FrozenStart::Word(w) = &f.start {
        crate::need!(
            w == "settled",
            FRZ_SETTING,
            "[{key}] start must be \"settled\" or a time in s, is \"{w}\""
        );
    }
    if let Some(s) = f.start_time() {
        crate::need!(
            s >= 0.0 && s.is_finite(),
            FRZ_SETTING,
            "[{key}] start must be \"settled\" or a time >= 0 s, is {s}"
        );
    }
    crate::need!(
        f.window_steps >= 2,
        FRZ_SETTING,
        "[{key}] window_steps must be >= 2, is {}",
        f.window_steps
    );
    crate::need!(
        f.settle_tolerance > 0.0 && f.settle_tolerance < 1.0,
        FRZ_SETTING,
        "[{key}] settle_tolerance must be in (0, 1), is {}",
        f.settle_tolerance
    );
    crate::need!(
        f.refresh_change > 0.0 && f.refresh_change.is_finite(),
        FRZ_SETTING,
        "[{key}] refresh_change must be > 0 K, is {}",
        f.refresh_change
    );
    crate::need!(
        f.iterations >= 1,
        FRZ_SETTING,
        "[{key}] iterations must be >= 1"
    );
    Ok(())
}

/// Checks `[energy.acceleration]` (ACC001).
fn validate_acceleration(case: &crate::case::FlowCase, a: &Acceleration) -> Result<()> {
    let key = "energy.acceleration";
    crate::need!(
        case.time.is_some(),
        ACC_SETTING,
        "[{key}] in a steady case: the steady temperatures do not depend on ρ c_p of a solid, \
         there is nothing to accelerate; it is for transient conjugate runs"
    );
    crate::need!(
        crate::cht::is_cht(case),
        ACC_SETTING,
        "[{key}] without solid regions: only the heat capacity of solids is scaled (in a fluid \
         it also multiplies the convection, which would change the answer)"
    );
    match &a.inertia_scale {
        InertiaScale::Factor(f) => crate::need!(
            *f > 0.0 && *f <= 1.0,
            ACC_SETTING,
            "[{key}] inertia_scale must be in (0, 1] or \"auto\", is {f}"
        ),
        InertiaScale::Word(w) => crate::need!(
            w == "auto",
            ACC_SETTING,
            "[{key}] inertia_scale must be a number in (0, 1] or \"auto\", is \"{w}\""
        ),
    }
    match &a.restore {
        Restore::Time(t) => crate::need!(
            *t >= 0.0 && t.is_finite(),
            ACC_SETTING,
            "[{key}] restore must be \"settled\" or a time >= 0 s, is {t}"
        ),
        Restore::Word(w) => crate::need!(
            w == "settled",
            ACC_SETTING,
            "[{key}] restore must be \"settled\" or a time in s, is \"{w}\""
        ),
    }
    for r in &a.regions {
        crate::need!(
            case.regions
                .get(r)
                .is_some_and(|x| x.kind == crate::materials::RegionKind::Solid),
            ACC_SETTING,
            "[{key}] regions: `{r}` is not a solid region of the case"
        );
    }
    if let Some(p) = a.period {
        crate::need!(
            p > 0.0 && p.is_finite(),
            ACC_SETTING,
            "[{key}] period must be > 0 s, is {p}"
        );
    }
    crate::need!(
        a.window_steps >= 2,
        ACC_SETTING,
        "[{key}] window_steps must be >= 2, is {}",
        a.window_steps
    );
    crate::need!(
        a.settle_tolerance > 0.0 && a.settle_tolerance < 1.0,
        ACC_SETTING,
        "[{key}] settle_tolerance must be in (0, 1), is {}",
        a.settle_tolerance
    );
    crate::need!(
        a.target() > 0.0 && a.target().is_finite(),
        ACC_SETTING,
        "[{key}] target_windows must be > 0, is {}",
        a.target()
    );
    Ok(())
}

/// `[[energy.sources]]`: the heat source of one region (stage 7d2, §15 of
/// the stage 7 specification): an electronic component, a heater, the
/// losses of a winding.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RegionSource {
    /// The region: a solid or fluid region of `[regions]`, or `fluid` for
    /// all fluid cells.
    pub region: String,
    /// The heat released in the region (W), spread over its cells by
    /// volume; a number or, in a transient case, an expression of t.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub power: Option<Scalar>,
    /// The power density (W/m³): a number or an expression of x, y, z and,
    /// in a transient case, t.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub power_density: Option<Scalar>,
}

fn yes() -> bool {
    true
}

fn default_prt() -> f64 {
    0.85
}

/// `[boundary.<name>.energy]`: the thermal condition of a boundary. One
/// of `temperature`, `heat_flux`, `heat_flow` and the pair
/// `heat_transfer_coefficient` + `ambient_temperature`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EnergyBoundary {
    /// Fixed temperature (K): a number or an expression of x, y, z and, in
    /// a transient case, t. Inlets and walls.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature: Option<Scalar>,
    /// Heat flux into the fluid (W/m²): a number or an expression. Walls,
    /// slip walls, symmetry planes, wedges.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_flux: Option<Scalar>,
    /// Heat flow into the fluid through the whole boundary (W), spread as
    /// a uniform flux over its area.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_flow: Option<f64>,
    /// Heat transfer coefficient h (W/(m² K)) to an ambient at
    /// `ambient_temperature`: the flux into the fluid is
    /// h (T_amb − T_wall).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_transfer_coefficient: Option<Scalar>,
    /// The ambient temperature of `heat_transfer_coefficient` (K).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ambient_temperature: Option<Scalar>,
}

/// Moves the per-physics sub-tables `energy` out of the
/// `[boundary.<name>]` tables of a parsed case and returns them.
pub fn extract_boundary_tables(v: &mut toml::Table) -> Result<BTreeMap<String, EnergyBoundary>> {
    let mut out = BTreeMap::new();
    let Some(toml::Value::Table(bs)) = v.get_mut("boundary") else {
        return Ok(out);
    };
    for (name, b) in bs.iter_mut() {
        let toml::Value::Table(t) = b else { continue };
        let Some(s) = t.remove("energy") else {
            continue;
        };
        let e: EnergyBoundary = s.try_into().map_err(|e: toml::de::Error| {
            diag::err(
                Diagnostic::new(
                    diag::ENE_BOUNDARY,
                    format!("[boundary.{name}.energy]: {}", e.message()),
                )
                .key(format!("boundary.{name}.energy")),
            )
        })?;
        out.insert(name.clone(), e);
    }
    Ok(out)
}

/// What a boundary of the flow is for the energy equation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Side {
    /// Walls and slip walls: every thermal condition.
    Wall,
    /// Symmetry planes and wedges: zero normal flux only.
    Symmetry,
    Inlet,
    Outlet,
    Periodic,
}

pub(crate) fn side(b: &BoundaryCondition) -> Side {
    use BoundaryCondition as B;
    match b {
        B::Wall { .. } | B::SlipWall => Side::Wall,
        B::Symmetry | B::Wedge => Side::Symmetry,
        B::VelocityInlet { .. }
        | B::MassFlowInlet { .. }
        | B::StaticPressureInlet { .. }
        | B::TotalPressureInlet { .. } => Side::Inlet,
        B::StaticPressureOutlet { .. }
        | B::TotalPressureOutlet { .. }
        | B::ConvectiveOutlet { .. } => Side::Outlet,
        B::Periodic { .. } => Side::Periodic,
    }
}

/// The condition of one boundary, compiled.
#[derive(Debug, Clone)]
pub(crate) enum Condition {
    /// Zero normal gradient (adiabatic walls, outlets).
    ZeroGradient,
    /// Fixed temperature.
    Temperature(expr::Value),
    /// Heat flux into the fluid (W/m²).
    HeatFlux(expr::Value),
    /// h and T_amb.
    Convection(expr::Value, expr::Value),
}

impl EnergyBoundary {
    fn count(&self) -> usize {
        usize::from(self.temperature.is_some())
            + usize::from(self.heat_flux.is_some())
            + usize::from(self.heat_flow.is_some())
            + usize::from(self.heat_transfer_coefficient.is_some())
    }
}

/// The compiled condition of every boundary of the case (by name), with
/// `area` the area of each boundary (for `heat_flow`).
pub(crate) fn conditions(
    case: &crate::case::FlowCase,
    env: &Env,
    area: &BTreeMap<String, f64>,
) -> Result<BTreeMap<String, Condition>> {
    let vars = Vars::SPACE.with(if case.time.is_some() {
        Vars::T
    } else {
        Vars::NONE
    });
    let mut out = BTreeMap::new();
    let coupled = crate::water_walls::boundaries(case);
    for (name, bc) in &case.boundary {
        let key = format!("boundary.{name}.energy");
        if coupled.contains(&name.as_str()) {
            // the water side of a circuit (SF7c): h and T_amb per face from
            // the first circuit solve on; adiabatic until then
            out.insert(
                name.clone(),
                Condition::Convection(expr::Value::Number(0.0), expr::Value::Number(0.0)),
            );
            continue;
        }
        let c = match case.boundary_energy.get(name) {
            None => Condition::ZeroGradient,
            Some(t) => {
                if let Some(v) = &t.temperature {
                    Condition::Temperature(expr::value(
                        v,
                        env,
                        vars,
                        &format!("{key}.temperature"),
                    )?)
                } else if let Some(v) = &t.heat_flux {
                    Condition::HeatFlux(expr::value(v, env, vars, &format!("{key}.heat_flux"))?)
                } else if let Some(q) = t.heat_flow {
                    let a = area.get(name).copied().unwrap_or(0.0);
                    crate::need!(
                        a > 0.0,
                        ENE_BOUNDARY,
                        "[{key}] heat_flow: the boundary has no area"
                    );
                    Condition::HeatFlux(expr::Value::Number(q / a))
                } else if let (Some(h), Some(ta)) =
                    (&t.heat_transfer_coefficient, &t.ambient_temperature)
                {
                    Condition::Convection(
                        expr::value(h, env, vars, &format!("{key}.heat_transfer_coefficient"))?,
                        expr::value(ta, env, vars, &format!("{key}.ambient_temperature"))?,
                    )
                } else {
                    Condition::ZeroGradient
                }
            }
        };
        let _ = bc;
        out.insert(name.clone(), c);
    }
    Ok(out)
}

/// Checks `[energy]` and the `energy` sub-tables of the boundaries.
pub fn validate(case: &crate::case::FlowCase, env: &Env) -> Result<()> {
    let Some(e) = &case.energy else {
        crate::need!(
            case.boundary_energy.is_empty(),
            THM_NO_ENERGY,
            "[boundary.{}.energy] without [energy]: the energy equation is solved only with an \
             [energy] table",
            case.boundary_energy
                .keys()
                .next()
                .map_or("", String::as_str)
        );
        return Ok(());
    };
    let transient = case.time.is_some();
    let time = if transient { Vars::T } else { Vars::NONE };
    expr::value(&e.initial, env, Vars::SPACE, "energy.initial")?;
    expr::value(&e.source, env, Vars::SPACE.with(time), "energy.source")?;
    if let Some(x) = &e.exact {
        expr::value(x, env, Vars::SPACE.with(time), "energy.exact")?;
    }
    crate::need!(
        e.turbulent_prandtl > 0.0 && e.turbulent_prandtl.is_finite(),
        CASE_RANGE,
        "[energy] turbulent_prandtl must be > 0, is {}",
        e.turbulent_prandtl
    );
    crate::need!(
        !case.scalars.contains_key(FIELD),
        ENE_UNSUPPORTED,
        "[scalars.{FIELD}] and [energy]: the energy equation solves the field `{FIELD}`; give \
         the passive scalar another name"
    );
    for (name, r) in &case.regions {
        match r.kind {
            crate::materials::RegionKind::Fluid => {
                // a fluid of its own (D194): its c_p and k
                if own_fluid(case, name) {
                    fluid_of(case, env, Some(name))?;
                }
            }
            crate::materials::RegionKind::Solid => {
                solid(case, env, name)?;
            }
            crate::materials::RegionKind::Porous => {}
        }
    }
    if crate::cht::is_cht(case) {
        crate::need!(
            case.sliding.is_empty(),
            ENE_UNSUPPORTED,
            "solid regions with a sliding mesh ([[sliding]]) are not in this version"
        );
        crate::need!(
            case.run.restart.is_none()
                || case
                    .energy
                    .as_ref()
                    .is_none_or(|e| e.acceleration.is_none()),
            ENE_UNSUPPORTED,
            "a restart with [energy.acceleration] is not in this version (the scaled ρ c_p of \
             the solids and the schedule of the scaling are not in the restart file)"
        );
        crate::need!(
            case.scalars.is_empty(),
            ENE_UNSUPPORTED,
            "passive scalars with solid regions are not in this version (they would diffuse \
             into the solids)"
        );
    }
    case.solver.energy.validate(transient)?;
    if let Some(a) = &e.acceleration {
        validate_acceleration(case, a)?;
    }
    if let Some(c) = &e.cycle {
        validate_cycle(case, c, e.acceleration.is_some())?;
    }
    if let Some(f) = &e.frozen_flow {
        validate_frozen(case, f, e)?;
    }
    for (b, t) in &case.boundary_energy {
        let key = format!("boundary.{b}.energy");
        let Some(bc) = case.boundary.get(b) else {
            return Err(diag::err(
                Diagnostic::new(
                    diag::ENE_BOUNDARY,
                    format!("[{key}] belongs to no boundary condition `[boundary.{b}]`"),
                )
                .key(key),
            ));
        };
        crate::need!(
            t.count() == 1,
            ENE_BOUNDARY,
            "[{key}] needs exactly one of `temperature`, `heat_flux`, `heat_flow` and \
             `heat_transfer_coefficient` (with `ambient_temperature`)"
        );
        crate::need!(
            t.heat_transfer_coefficient.is_some() == t.ambient_temperature.is_some(),
            ENE_BOUNDARY,
            "[{key}]: `heat_transfer_coefficient` and `ambient_temperature` go together"
        );
        let sd = side(bc);
        let ok = match sd {
            Side::Wall => true,
            Side::Inlet => t.temperature.is_some(),
            Side::Symmetry | Side::Outlet | Side::Periodic => false,
        };
        crate::need!(
            ok,
            THM_SIDE,
            "[{key}]: a {} takes {}",
            bc.type_name(),
            match sd {
                Side::Inlet => "a `temperature`, not a flux or a convection condition",
                Side::Symmetry => "no energy condition (no heat crosses it)",
                Side::Outlet => "no energy condition (outlets keep a zero normal gradient)",
                _ =>
                    "no energy condition (T crosses the periodic pair; the periodic heat-flux \
                      mode of the energy equation is not in this version)",
            }
        );
        let vars = Vars::SPACE.with(time);
        for (k, v) in [
            ("temperature", &t.temperature),
            ("heat_flux", &t.heat_flux),
            ("heat_transfer_coefficient", &t.heat_transfer_coefficient),
            ("ambient_temperature", &t.ambient_temperature),
        ] {
            if let Some(v) = v {
                let x = expr::value(v, env, vars, &format!("{key}.{k}"))?;
                let Some(c) = x.constant() else { continue };
                match k {
                    "heat_transfer_coefficient" => crate::need!(
                        c > 0.0 && c.is_finite(),
                        THM_RANGE,
                        "[{key}] heat_transfer_coefficient must be > 0, is {c}"
                    ),
                    "temperature" | "ambient_temperature" => crate::need!(
                        c > 0.0 && c.is_finite(),
                        THM_RANGE,
                        "[{key}] {k} must be > 0 K, is {c}"
                    ),
                    _ => crate::need!(
                        c.is_finite(),
                        THM_RANGE,
                        "[{key}] {k} must be finite, is {c}"
                    ),
                }
            }
        }
        if let Some(q) = t.heat_flow {
            crate::need!(
                q.is_finite(),
                THM_RANGE,
                "[{key}] heat_flow must be a finite number"
            );
        }
    }
    for (b, bc) in &case.boundary {
        crate::need!(
            side(bc) != Side::Inlet
                || case
                    .boundary_energy
                    .get(b)
                    .is_some_and(|t| t.temperature.is_some()),
            ENE_INLET,
            "[boundary.{b}] is an inlet: with [energy] it needs `[boundary.{b}.energy] \
             temperature`"
        );
    }
    if let Some(c) = expr::value(&e.initial, env, Vars::SPACE, "energy.initial")?.constant() {
        crate::need!(
            c > 0.0 && c.is_finite(),
            THM_RANGE,
            "[energy] initial must be > 0 K, is {c}"
        );
    }
    // a steady temperature needs a level: an inlet, a temperature or a
    // convection boundary (THM001); without one, a net heat input has no
    // way out (THM002)
    if !transient {
        let anchored = case.boundary.iter().any(|(b, bc)| {
            side(bc) == Side::Inlet
                || case.boundary_energy.get(b).is_some_and(|t| {
                    t.temperature.is_some() || t.heat_transfer_coefficient.is_some()
                })
        });
        if !anchored {
            let source = expr::value(&e.source, env, Vars::SPACE, "energy.source")?;
            let input = e.viscous_heating
                || source.constant() != Some(0.0)
                || case.boundary_energy.values().any(|t| {
                    t.heat_flow.is_some_and(|q| q != 0.0)
                        || t.heat_flux.as_ref().is_some_and(|q| {
                            expr::value(q, env, Vars::SPACE, "heat_flux")
                                .map_or(true, |v| v.constant() != Some(0.0))
                        })
                });
            if input {
                crate::need!(
                    false,
                    THM_DRIFT,
                    "a steady run with a heat input and no way for the heat to leave: give a \
                     boundary a `temperature` or a `heat_transfer_coefficient`, add an inflow, \
                     or run transient ([time])"
                );
            }
            crate::need!(
                false,
                THM_LEVEL,
                "the boundaries leave the steady temperature undetermined (only adiabatic and \
                 heat-flux walls, no inflow): give a boundary a `temperature` or a \
                 `heat_transfer_coefficient`, add an inflow, or run transient ([time])"
            );
        }
    }
    fluid(case, env)?;
    Ok(())
}

/// The findings of the thermal set-up for the boundary-condition check
/// (D30): the walls left adiabatic (THM005), listed once.
pub fn setup_findings(case: &crate::case::FlowCase) -> Vec<crate::bc::Finding> {
    use crate::bc::{Finding, Level};
    // the fluid view of a conjugate case: the walls of the case as written
    // (the interfaces are no walls of the energy equation)
    let case = case.fluid_view.as_ref().map_or(case, |v| &v.full_case);
    if case.energy.is_none() {
        return Vec::new();
    }
    let coupled = crate::water_walls::boundaries(case);
    let walls: Vec<&str> = case
        .boundary
        .iter()
        .filter(|(b, bc)| {
            side(bc) == Side::Wall
                && !case.boundary_energy.contains_key(*b)
                && !coupled.contains(&b.as_str())
        })
        .map(|(b, _)| b.as_str())
        .collect();
    if walls.is_empty() {
        return Vec::new();
    }
    vec![Finding {
        level: Level::Warning,
        code: diag::THM_ADIABATIC,
        boundary: None,
        message: format!(
            "walls without [boundary.<name>.energy] are adiabatic: {}",
            walls.join(", ")
        ),
    }]
}

/// The range of the temperatures that the case prescribes as numbers (the
/// start value, boundary and ambient temperatures), for the range checks
/// of buoyancy (BUO002); `None` without the energy equation.
pub fn temperature_range(case: &crate::case::FlowCase) -> Option<(f64, f64)> {
    let case = case.fluid_view.as_ref().map_or(case, |v| &v.full_case);
    case.energy.as_ref()?;
    let mut v: Vec<f64> = Vec::new();
    if let Some(t) = case.energy.as_ref().and_then(|e| e.initial.number()) {
        v.push(t);
    }
    for b in case.boundary_energy.values() {
        for s in [&b.temperature, &b.ambient_temperature]
            .into_iter()
            .flatten()
        {
            if let Some(t) = s.number() {
                v.push(t);
            }
        }
    }
    let lo = v.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = v.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (lo <= hi).then_some((lo, hi))
}

/// The thermal properties of the fluid as functions of T.
#[derive(Debug, Clone)]
pub struct Fluid {
    /// The material (`None`: the `[fluid]` shorthand).
    pub material: Option<String>,
    /// ρ (kg/m³), constant in 7a1.
    pub rho: f64,
    /// c_p(T) (J/(kg K)).
    pub cp: Curve,
    /// k(T) (W/(m K)).
    pub k: Curve,
    /// μ(T) (Pa s), when the viscosity is a form of T (else the flow keeps
    /// its viscosity).
    pub mu: Option<Curve>,
}

/// The thermal properties of the `[fluid]` (its material, with the
/// `[fluid]` numbers replacing the material's): MAT007 without c_p or k.
pub fn fluid(case: &crate::case::FlowCase, env: &Env) -> Result<Fluid> {
    fluid_of(case, env, None)
}

/// True when fluid region `name` has a material of its own (D194: another
/// fluid than `[fluid]`).
pub fn own_fluid(case: &crate::case::FlowCase, name: &str) -> bool {
    case.regions.get(name).is_some_and(|r| {
        r.kind == crate::materials::RegionKind::Fluid
            && r.material.is_some()
            && r.material != case.fluid.material
    })
}

/// The thermal properties of a fluid region with a material of its own
/// (D194), or of the `[fluid]` (`region` `None`): MAT007 without c_p or k.
/// A region's own material takes none of the `[fluid]` numbers; the
/// density is that of the property program (at the region's temperature).
pub fn fluid_of(case: &crate::case::FlowCase, env: &Env, region: Option<&str>) -> Result<Fluid> {
    let mut cond: Conditions = case.fluid.conditions();
    let own = region.filter(|r| own_fluid(case, r));
    let sec = region.and_then(|r| case.regions.get(r));
    if let Some(t) = sec.and_then(|r| r.temperature) {
        cond.temperature = t;
    }
    if let Some(p) = sec.and_then(|r| r.operating_pressure) {
        cond.pressure = p;
    }
    let material = match own {
        Some(_) => sec.and_then(|r| r.material.clone()),
        None => case.fluid.material.clone(),
    };
    let (m, owner) = match &material {
        Some(n) => (lookup(n, &case.materials)?, format!("materials.{n}")),
        None => (MaterialSection::default(), "fluid".to_string()),
    };
    let (given_cp, given_k, given_mu) = if own.is_some() {
        (None, None, None)
    } else {
        (
            case.fluid.specific_heat.clone(),
            case.fluid.conductivity.clone(),
            case.fluid.viscosity.clone(),
        )
    };
    let rho = match region {
        Some(r) => {
            let p = case.properties()?;
            p.regions
                .iter()
                .find(|x| x.0 == r)
                .map_or(p.density(), |x| {
                    x.2.density.value.constant().unwrap_or(f64::NAN)
                })
        }
        None => case.rho(),
    };
    let pick = |given: &Option<Scalar>, own: &Option<Property>| -> Option<(Property, bool)> {
        match (given, own) {
            (Some(s), _) => Some((Property::from(s.clone()), true)),
            (None, Some(p)) => Some((p.clone(), false)),
            (None, None) => None,
        }
    };
    let own_k = match &m.conductivity {
        Some(Conductivity::Isotropic(p)) => Some(p.clone()),
        Some(Conductivity::Orthotropic(_)) => {
            return Err(diag::err(
                Diagnostic::new(
                    diag::ENE_UNSUPPORTED,
                    format!(
                        "{owner}: an orthotropic conductivity belongs to a solid; a fluid takes \
                         one value"
                    ),
                )
                .key(format!("{owner}.conductivity")),
            ));
        }
        None => None,
    };
    let cp = pick(&given_cp, &m.specific_heat);
    let k = pick(&given_k, &own_k);
    let missing: Vec<&str> = [
        ("specific_heat", cp.is_none()),
        ("conductivity", k.is_none()),
    ]
    .iter()
    .filter(|x| x.1)
    .map(|x| x.0)
    .collect();
    if !missing.is_empty() {
        let what = match &material {
            Some(n) if own.is_some() => format!(
                "material `{n}` of fluid region `{}` has no {}; the energy equation needs c_p \
                 and k of every fluid (give them in [materials.{n}] with `base = \"{n}\"`)",
                own.unwrap_or(""),
                missing.join(" and no ")
            ),
            Some(n) => format!(
                "material `{n}` has no {}; the energy equation needs c_p and k of the fluid \
                 (give them in [materials.{n}] with `base = \"{n}\"`, or in [fluid])",
                missing.join(" and no ")
            ),
            None => format!(
                "[fluid] has no {}; the energy equation needs c_p and k of the fluid (give \
                 `specific_heat` and `conductivity` in [fluid], or a `material`)",
                missing.join(" and no ")
            ),
        };
        return Err(diag::err(
            Diagnostic::new(diag::MAT_THERMAL, what).key("fluid"),
        ));
    }
    let (cp, cp_given) = cp.unwrap_or((Property::Number(f64::NAN), false));
    let (k, k_given) = k.unwrap_or((Property::Number(f64::NAN), false));
    let label = |given: bool| {
        if given {
            "fluid".to_string()
        } else {
            owner.clone()
        }
    };
    let cp = curve(&cp, "specific_heat", &label(cp_given), &m, env, cond)?;
    let k = curve(&k, "conductivity", &label(k_given), &m, env, cond)?;
    // μ(T): only a form of T of the material's own viscosity (a [fluid]
    // viscosity, a number or an expression of x, y, z stays as it is)
    let mu = match (&given_mu, &m.viscosity) {
        (None, Some(p @ Property::Form(f)))
            if !matches!(f, crate::materials::Form::Expression { .. }) =>
        {
            let c = curve(p, "viscosity", &owner, &m, env, cond)?;
            (!c.is_constant()).then_some(c)
        }
        _ => None,
    };
    Ok(Fluid {
        material,
        rho,
        cp,
        k,
        mu,
    })
}

/// The thermal properties of a solid region (stage 7d1): ρ at the case
/// temperature, c_p(T) and k(T) of its material, an orthotropic k with the
/// axes of the region or of the material (7d2); MAT006 without one of
/// them, MAT008 for an orthotropic k without axes.
pub fn solid(case: &crate::case::FlowCase, env: &Env, region: &str) -> Result<Material> {
    let key = format!("regions.{region}.material");
    let Some(name) = case.regions.get(region).and_then(|r| r.material.clone()) else {
        return Err(diag::err(
            Diagnostic::new(
                diag::MAT_SOLID,
                format!(
                    "[regions.{region}] is a solid without a material: give `material`, with \
                     density, specific_heat and conductivity"
                ),
            )
            .key(key),
        ));
    };
    solid_material(case, env, region, &name, key)
}

/// The material `name` as the solid of `region` (a solid region, or the
/// matrix of a porous region, D271), `key` the key that names it.
pub fn solid_material(
    case: &crate::case::FlowCase,
    env: &Env,
    region: &str,
    name: &str,
    key: String,
) -> Result<Material> {
    let cond: Conditions = case.fluid.conditions();
    let name = name.to_string();
    let m = lookup(&name, &case.materials)?;
    let owner = format!("materials.{name}");
    let mut aniso = None;
    let k = match &m.conductivity {
        Some(Conductivity::Isotropic(p)) => Some(p.clone()),
        Some(Conductivity::Orthotropic(o)) => {
            // the axes of the region (an oriented part of one material),
            // else those of the material (MAT008 without either)
            let region_axes = case.regions.get(region).and_then(|r| r.axes);
            let mut oo = (**o).clone();
            if region_axes.is_some() {
                oo.axes = region_axes;
            }
            let key = if region_axes.is_some() {
                format!("regions.{region}.axes")
            } else {
                format!("{owner}.conductivity")
            };
            let axes = crate::materials::axes_of(&oo, &key)?;
            let mut ks = Vec::with_capacity(3);
            for (i, p) in o.values.iter().enumerate() {
                ks.push(curve(
                    p,
                    &format!("conductivity.values[{i}]"),
                    &owner,
                    &m,
                    env,
                    cond,
                )?);
            }
            aniso = Some(Aniso {
                k: [ks[0].clone(), ks[1].clone(), ks[2].clone()],
                axes,
            });
            // (the first principal value stands for k where one number is
            // asked for; the energy equation takes the tensor)
            o.values.first().cloned()
        }
        None => None,
    };
    crate::need!(
        // (the axes of a porous region orient its resistance and the
        // conductivity of its structure, so a one-valued material is fine)
        aniso.is_some()
            || case.regions.get(region).is_none_or(|r| {
                r.axes.is_none() || r.kind == crate::materials::RegionKind::Porous
            }),
        MAT_AXES,
        "[regions.{region}] axes orient an orthotropic conductivity; material `{name}` has one value"
    );
    let missing: Vec<&str> = [
        ("density", m.density.is_none()),
        ("specific_heat", m.specific_heat.is_none()),
        ("conductivity", k.is_none()),
    ]
    .iter()
    .filter(|x| x.1)
    .map(|x| x.0)
    .collect();
    if !missing.is_empty() {
        return Err(diag::err(
            Diagnostic::new(
                diag::MAT_SOLID,
                format!(
                    "material `{name}` of the solid region `{region}` has no {}: a solid needs \
                     density, specific_heat and conductivity",
                    missing.join(", ")
                ),
            )
            .key(key),
        ));
    }
    let rho = curve(
        m.density.as_ref().unwrap(),
        "density",
        &owner,
        &m,
        env,
        cond,
    )?
    .at(cond.temperature);
    crate::need!(
        rho > 0.0 && rho.is_finite(),
        MAT_SOLID,
        "material `{name}`: the density must be > 0, is {rho}"
    );
    Ok(Material {
        name: name.clone(),
        rho,
        cp: curve(
            m.specific_heat.as_ref().unwrap(),
            "specific_heat",
            &owner,
            &m,
            env,
            cond,
        )?,
        k: curve(k.as_ref().unwrap(), "conductivity", &owner, &m, env, cond)?,
        aniso,
        mu: None,
        fluid: false,
        t_ref: cond.temperature,
    })
}

/// The material of the energy equation of a fluid (D194: one per fluid of
/// the case), `name` its label in messages (`fluid` for the `[fluid]`).
pub fn fluid_material(f: &Fluid, name: &str, t_ref: f64) -> Material {
    Material {
        name: name.to_string(),
        rho: f.rho,
        cp: f.cp.clone(),
        k: f.k.clone(),
        aniso: None,
        mu: f.mu.clone(),
        fluid: true,
        t_ref,
    }
}

/// The thermal properties of one material of the energy equation.
#[derive(Debug, Clone)]
pub struct Material {
    /// The material (or `fluid`).
    pub name: String,
    /// ρ (kg/m³), constant.
    pub rho: f64,
    /// c_p(T) (J/(kg K)).
    pub cp: Curve,
    /// k(T) (W/(m K)); of an orthotropic solid the first principal value
    /// (the equation takes [`Material::tensor`]).
    pub k: Curve,
    /// An orthotropic conductivity (stage 7d2).
    pub aniso: Option<Aniso>,
    /// μ(T) (Pa s) of a fluid whose viscosity follows T.
    pub mu: Option<Curve>,
    /// A fluid (D194: the energy equation has one material per fluid).
    pub fluid: bool,
    /// The temperature of the properties when they do not follow T (K).
    pub t_ref: f64,
}

/// An orthotropic conductivity K = Σᵢ kᵢ(T) eᵢ eᵢᵀ (D116, L48): the three
/// principal values and the orthonormal axes.
#[derive(Debug, Clone)]
pub struct Aniso {
    /// k₁(T), k₂(T), k₃(T) (W/(m K)).
    pub k: [Curve; 3],
    /// e₁, e₂, e₃.
    pub axes: [[f64; 3]; 3],
}

impl Material {
    /// The conductivity tensor at `t` (K): K = Σᵢ kᵢ eᵢ eᵢᵀ, or k I.
    pub fn tensor(&self, t: f64) -> [[f64; 3]; 3] {
        match &self.aniso {
            None => {
                let k = self.k.at(t);
                [[k, 0.0, 0.0], [0.0, k, 0.0], [0.0, 0.0, k]]
            }
            Some(a) => {
                let mut m = [[0.0; 3]; 3];
                for (c, e) in a.k.iter().zip(&a.axes) {
                    let ki = c.at(t);
                    for i in 0..3 {
                        for j in 0..3 {
                            m[i][j] += ki * e[i] * e[j];
                        }
                    }
                }
                m
            }
        }
    }

    /// The conductivity along the unit normal `n` at `t`: n·K·n (k of an
    /// isotropic material).
    pub fn k_normal(&self, t: f64, n: [f64; 3]) -> f64 {
        if self.aniso.is_none() {
            return self.k.at(t);
        }
        let m = self.tensor(t);
        (0..3)
            .map(|i| n[i] * (0..3).map(|j| m[i][j] * n[j]).sum::<f64>())
            .sum()
    }
}

/// The temperature of zero sensible enthalpy (K).
pub const H_REF_T: f64 = 298.15;

/// The sensible enthalpy h_s(T) = ∫ c_p dT from [`H_REF_T`] of a c_p
/// curve: exact for a constant c_p, else the trapezoidal integral on a
/// grid of [`Enthalpy::STEP`] K, linear between the grid points (and with
/// the end values of c_p outside the grid). The energy equation is
/// conservative in this h_s (D150): the convection and storage terms of the
/// matrix in T are corrected with the difference of the h-form at the
/// current temperature, so the residual at convergence is that of h_s.
#[derive(Debug, Clone, PartialEq)]
pub struct Enthalpy {
    cp: Option<f64>,
    t0: f64,
    h: Vec<f64>,
    cp_lo: f64,
    cp_hi: f64,
}

impl Enthalpy {
    /// The grid step (K).
    pub const STEP: f64 = 0.25;

    /// h_s of a constant c_p.
    pub fn constant(cp: f64) -> Self {
        Self {
            cp: Some(cp),
            t0: 0.0,
            h: Vec::new(),
            cp_lo: cp,
            cp_hi: cp,
        }
    }

    /// h_s of the curve `cp`.
    pub fn of(cp: &Curve) -> Self {
        if cp.is_constant() {
            return Self::constant(cp.at(H_REF_T));
        }
        let [lo, hi] = Curve::SAMPLE_RANGE;
        let n = ((hi - lo) / Self::STEP).round() as usize + 1;
        let mut h = Vec::with_capacity(n);
        let mut acc = 0.0;
        let mut prev = cp.at(lo);
        h.push(0.0);
        for i in 1..n {
            let c = cp.at(lo + i as f64 * Self::STEP);
            acc += 0.5 * (prev + c) * Self::STEP;
            h.push(acc);
            prev = c;
        }
        let mut e = Self {
            cp: None,
            t0: lo,
            h,
            cp_lo: cp.at(lo),
            cp_hi: prev,
        };
        let shift = e.at(H_REF_T);
        for v in &mut e.h {
            *v -= shift;
        }
        e
    }

    /// True for a constant c_p (no correction needed).
    pub fn is_constant(&self) -> bool {
        self.cp.is_some()
    }

    /// h_s (J/kg) at `t` (K).
    pub fn at(&self, t: f64) -> f64 {
        if let Some(cp) = self.cp {
            return cp * (t - H_REF_T);
        }
        let n = self.h.len();
        let x = (t - self.t0) / Self::STEP;
        if x <= 0.0 {
            return self.h[0] + self.cp_lo * (t - self.t0);
        }
        let last = (n - 1) as f64;
        if x >= last {
            return self.h[n - 1] + self.cp_hi * (t - (self.t0 + last * Self::STEP));
        }
        let i = (x.floor() as usize).min(n - 2);
        let w = x - i as f64;
        self.h[i] + w * (self.h[i + 1] - self.h[i])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Env;
    use crate::materials::{CaseMaterials, Conditions, Property, curve};

    #[test]
    fn enthalpy_integrates_the_specific_heat() {
        let env = Env::default();
        let m = MaterialSection::default();
        let cond = Conditions::default();
        // c_p = 1000 + 0.5 T: h = 1000 (T − T_r) + 0.25 (T² − T_r²)
        let cp = curve(
            &Property::Expr("1000 + 0.5 * T".into()),
            "specific_heat",
            "test",
            &m,
            &env,
            cond,
        )
        .unwrap();
        let h = Enthalpy::of(&cp);
        assert!(!h.is_constant());
        for t in [250.0, 298.15, 300.0, 412.3, 1500.0] {
            let exact = 1000.0 * (t - H_REF_T) + 0.25 * (t * t - H_REF_T * H_REF_T);
            // linear interpolation of the 0.25 K table: |err| ≤ c_p' Δ²/8 = 0.004 J/kg
            assert!((h.at(t) - exact).abs() < 5e-3, "{t}: {} {exact}", h.at(t));
        }
        let c = Enthalpy::constant(4184.0);
        assert_eq!(c.at(H_REF_T + 1.0), 4184.0);
        let _ = CaseMaterials::default();
    }
}
