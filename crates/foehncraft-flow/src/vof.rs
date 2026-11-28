//! Two immiscible fluids by the volume-of-fluid method (L3, stage MP1).
//!
//! One set of momentum equations carries the mixture; the volume fraction
//! α of phase 1 (α = 1 in phase 1, 0 in phase 2) moves with the flow,
//!
//! ∂α/∂t + ∇·(α u) = 0,
//!
//! and sets the density ρ = α ρ₁ + (1 − α) ρ₂ and the viscosity
//! μ = α μ₁ + (1 − α) μ₂ per cell and face.
//!
//! **A time step.** At the start of each step α is advanced explicitly
//! with the volume fluxes F_f = ṁ_f/ρ_f of the step before, in sub-steps
//! that keep the outflow Courant number of every cell below
//! `max_courant`; the face values come from THINC (default), CICSAM or
//! upwind ([`foehncraft_kernels::vof`]). The mean face value of the
//! sub-steps gives the face density, so that the mass flux ṁ_f = ρ_f F_f
//! that convects the momentum is the one that moved α. The new densities
//! enter the time term of the momentum rows as level n + 1 (level n is the
//! density before the advection). The continuity of the two incompressible
//! phases is that of the volume, Σ_f F_f = 0: the pressure-correction rows
//! take the face coefficients and the fluxes divided by the face density
//! ([`FlowDevice::quotient`]); the corrections of the mass fluxes are the
//! same as for one fluid.
//!
//! **Surface tension** by the continuum surface force (Brackbill, Kothe &
//! Zemach 1992) in balanced-force form (Francois et al. 2006): the face
//! jump σ κ_f (α_N − α_P) is an offset of the reduced pressure jump, beside
//! the gravity offset ½ (b_P + b_N)·d of the balanced-force interpolation
//! (D154), so a static interface of constant curvature is in discrete
//! equilibrium with p′ = σ κ α. The curvature κ = −∇·(∇α̃/|∇α̃|) is taken
//! from α smoothed `smoothing` times; walls and symmetry planes have a
//! contact angle of 90°.
//!
//! **Ledger.** The update is conservative: the volume of phase 1 changes
//! by the boundary fluxes alone, except for the part clipped to keep
//! 0 ≤ α ≤ 1, which is booked and reported (`clipped_volume`).
//!
//! Case file:
//!
//! ```toml
//! [vof]
//! surface_tension = 0.072      # σ (N/m), default 0
//! scheme = "thinc"             # "thinc" (default) | "cicsam" | "upwind"
//! thinc_beta = 3.5             # steepness of THINC
//! max_courant = 0.25           # outflow Courant number of the α sub-steps
//! smoothing = 2                # passes of α̃ for the curvature
//! initial = 0.0                # α of the cells outside every region
//!
//! [[vof.phase]]                # phase 1 (α = 1), then phase 2
//! name = "water"
//! density = 998.2
//! viscosity = 1.0e-3
//!
//! [[vof.phase]]
//! name = "air"
//! density = 1.205
//! viscosity = 1.82e-5
//!
//! [[vof.region]]               # start regions, in order
//! shape = "box"                # "box" | "sphere" | "cylinder" | "field"
//! lower = [0.0, 0.0, 0.0]
//! upper = [0.146, 0.292, 1.0]
//! alpha = 1.0
//!
//! [vof.inflow]                 # α of the inflow per boundary (default 0)
//! inlet = 1.0
//! ```
//!
//! **Energy and phase change (stage MP2, SF14).** With `[energy]` each
//! phase gives `specific_heat` and `conductivity`; the energy equation of
//! the run (in T, D150) then takes the mixture ρ c_p = α ρ₁c₁ + (1 − α)
//! ρ₂c₂ per cell, (ρ c_p)_f of the advected face α in the heat-capacity
//! flux Ċ_f = F_f (ρ c_p)_f, and the conductivity in series, 1/λ = α/λ₁ +
//! (1 − α)/λ₂ per cell and the harmonic mean of the two cells per face
//! (the resistance of a layered interface).
//!
//! `[vof.phase_change]` adds evaporation and condensation by the energy
//! of the interface (after Rattner & Garimella 2014): the interface cells
//! (0 < α < 1, and the cell on the side of a sharp face that can change
//! phase: liquid when it is above T_sat, vapour below) are held at the
//! saturation temperature by the implicit sink S = −G (T − T_sat),
//! G = `relaxation` ρ c_p V/Δt, and the heat G (T − T_sat) of a step
//! turns into the mass rate ṁ = G (T − T_sat)/h_lv of the cell (positive:
//! evaporation), limited to the liquid (vapour) the cell holds; the heat
//! that the limit leaves is given back to the cell in the next step. The
//! rate ṁ then (one step later) adds the volume ṁ (1/ρ_v − 1/ρ_l) to the
//! continuity rows of the volume and takes ṁ/ρ_l of liquid out of α in
//! the next advance (or adds ṁ/ρ_v of vapour, if phase 1 is the vapour),
//! so that the volume source the fluxes carry is the one α loses: the
//! phase masses change by the evaporated mass alone.
//!
//! ```toml
//! [[vof.phase]]
//! name = "water"
//! density = 958.4
//! viscosity = 2.82e-4
//! specific_heat = 4216.0
//! conductivity = 0.679
//!
//! [vof.phase_change]
//! liquid = "water"                # the liquid phase (default: the denser)
//! saturation_temperature = 373.15 # T_sat (K)
//! latent_heat = 2.257e6           # h_lv (J/kg)
//! # pressure = 7.0e6              # or water: T_sat and h_lv of IAPWS-IF97
//! relaxation = 100.0              # G in units of ρ c_p V/Δt
//! ```

use crate::device::{Buf, Buf3, FlowDevice, VofBufs, VofOp};
use crate::diag::{self, Diagnostic};
use anyhow::Result;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The face scheme of α.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Scheme {
    /// THINC with the time-integrated face value (Xiao et al. 2005, 2011).
    #[default]
    Thinc,
    /// CICSAM (Ubbink & Issa 1999).
    Cicsam,
    /// First-order upwind (smears the interface; for checks).
    Upwind,
}

impl Scheme {
    fn code(self) -> u32 {
        use foehncraft_kernels::vof as k;
        match self {
            Self::Thinc => k::SCHEME_THINC,
            Self::Cicsam => k::SCHEME_CICSAM,
            Self::Upwind => k::SCHEME_UPWIND,
        }
    }

    /// The name of the case file.
    pub fn name(self) -> &'static str {
        match self {
            Self::Thinc => "thinc",
            Self::Cicsam => "cicsam",
            Self::Upwind => "upwind",
        }
    }
}

/// One `[[vof.phase]]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PhaseSection {
    /// Name of the phase (output, report).
    pub name: String,
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Dynamic viscosity μ (Pa s).
    pub viscosity: f64,
    /// Specific heat c_p (J/(kg K)); needed with `[energy]` (MP2).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specific_heat: Option<f64>,
    /// Thermal conductivity λ (W/(m K)); needed with `[energy]` (MP2).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conductivity: Option<f64>,
}

/// `[vof.phase_change]`: evaporation and condensation at the interface
/// (stage MP2, SF14).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PhaseChangeSection {
    /// The name of the liquid phase (default: the denser phase).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub liquid: Option<String>,
    /// Saturation temperature T_sat (K).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub saturation_temperature: Option<f64>,
    /// Latent heat h_lv (J/kg).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latent_heat: Option<f64>,
    /// Water at this pressure (Pa): T_sat and h_lv of IAPWS-IF97 (stage
    /// SF7a), where they are not given; also σ(T_sat) where `[vof]
    /// surface_tension` is not given.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pressure: Option<f64>,
    /// The sink G that holds the interface cells at T_sat, in units of
    /// ρ c_p V/Δt (default 100: T − T_sat ends at about 1 % of what one
    /// step of heating would give).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relaxation: Option<f64>,
}

impl PhaseChangeSection {
    /// T_sat and h_lv: given, or of IAPWS-IF97 water at `pressure`.
    pub fn saturation(&self) -> Option<(f64, f64)> {
        let sat = self
            .pressure
            .and_then(foehncraft_kernels::steam::saturation_p);
        let t = self.saturation_temperature.or(sat.as_ref().map(|s| s.t))?;
        let h = self.latent_heat.or(sat.as_ref().map(|s| s.h_fg()))?;
        Some((t, h))
    }

    /// The relaxation of the sink.
    pub fn relax(&self) -> f64 {
        self.relaxation.unwrap_or(100.0)
    }
}

/// One `[[vof.region]]`: α = `alpha` inside the shape (cell fractions by
/// the quadrature of the immersed bodies, D226), applied in order over the
/// cells' current α.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RegionSection {
    /// Shape.
    pub shape: crate::immersed::Shape,
    /// Cylinder: a point on the axis; sphere: the centre (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub centre: Option<[f64; 3]>,
    /// Radius of a cylinder or sphere (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub radius: Option<f64>,
    /// Axis of a cylinder (default z).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub axis: Option<[f64; 3]>,
    /// Length of a cylinder (m); default endless.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub length: Option<f64>,
    /// Lower corner of a box (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lower: Option<[f64; 3]>,
    /// Upper corner of a box (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper: Option<[f64; 3]>,
    /// The fraction (0 to 1) of a `field` region, an expression of x, y, z.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fraction: Option<crate::expr::Scalar>,
    /// α inside (default 1).
    #[serde(default = "one")]
    pub alpha: f64,
}

fn one() -> f64 {
    1.0
}

/// `[vof]`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VofSection {
    /// The two phases: phase 1 (α = 1) first.
    #[serde(default)]
    pub phase: Vec<PhaseSection>,
    /// Surface tension σ (N/m), default 0.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub surface_tension: Option<f64>,
    /// Face scheme of α.
    #[serde(default)]
    pub scheme: Scheme,
    /// Steepness β of THINC (default 3.5).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thinc_beta: Option<f64>,
    /// The largest outflow Courant number of a cell in one α sub-step
    /// (default 0.25, at most 1).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_courant: Option<f64>,
    /// Smoothing passes of α for the curvature (default 2).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub smoothing: Option<u32>,
    /// α outside every region at the start (default 0).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial: Option<f64>,
    /// Start regions.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub region: Vec<RegionSection>,
    /// α of the inflow per boundary (default 0).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub inflow: BTreeMap<String, f64>,
    /// Evaporation and condensation (MP2).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub phase_change: Option<PhaseChangeSection>,
}

impl VofSection {
    /// True without `[vof]`.
    pub fn is_empty(&self) -> bool {
        self.phase.is_empty()
    }

    /// σ (N/m): given, or of IAPWS-IF97 water at T_sat when the phase
    /// change names a pressure.
    pub fn sigma(&self) -> f64 {
        if let Some(s) = self.surface_tension {
            return s;
        }
        match self.phase_change.as_ref().filter(|p| p.pressure.is_some()) {
            Some(p) => p
                .saturation()
                .map_or(0.0, |(t, _)| foehncraft_kernels::steam::surface_tension(t)),
            None => 0.0,
        }
    }

    /// The index of the liquid phase of the phase change.
    pub fn liquid(&self) -> usize {
        let named = self
            .phase_change
            .as_ref()
            .and_then(|p| p.liquid.as_ref())
            .and_then(|n| self.phase.iter().position(|p| &p.name == n));
        named.unwrap_or(
            if self.phase.len() == 2 && self.phase[1].density > self.phase[0].density {
                1
            } else {
                0
            },
        )
    }

    /// β of THINC.
    pub fn beta(&self) -> f64 {
        self.thinc_beta.unwrap_or(3.5)
    }

    /// The Courant number of a sub-step.
    pub fn courant(&self) -> f64 {
        self.max_courant.unwrap_or(0.25)
    }

    /// Smoothing passes.
    pub fn passes(&self) -> u32 {
        self.smoothing.unwrap_or(2)
    }
}

fn bad(key: &str, msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::VOF_CONFIG, msg).key(key.to_string()))
}

fn unsupported(msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::VOF_UNSUPPORTED, msg))
}

/// Ranges of `[vof]` and the settings it does not work with yet (VOF001,
/// VOF002).
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let v = &case.vof;
    if v.is_empty() {
        return Ok(());
    }
    if v.phase.len() != 2 {
        return Err(bad(
            "vof.phase",
            format!(
                "[vof] needs exactly two [[vof.phase]] tables, not {}",
                v.phase.len()
            ),
        ));
    }
    for (i, p) in v.phase.iter().enumerate() {
        if !(p.density > 0.0
            && p.density.is_finite()
            && p.viscosity > 0.0
            && p.viscosity.is_finite())
        {
            return Err(bad(
                "vof.phase",
                format!(
                    "[[vof.phase]] {} (`{}`): density and viscosity must be > 0",
                    i + 1,
                    p.name
                ),
            ));
        }
    }
    let s = v.sigma();
    if !(s >= 0.0 && s.is_finite()) {
        return Err(bad(
            "vof.surface_tension",
            format!("[vof] surface_tension {s} must be >= 0"),
        ));
    }
    let b = v.beta();
    if !(b > 0.0 && b <= 20.0) {
        return Err(bad(
            "vof.thinc_beta",
            format!("[vof] thinc_beta {b} must be in (0, 20]"),
        ));
    }
    let c = v.courant();
    if !(c > 0.0 && c <= 1.0) {
        return Err(bad(
            "vof.max_courant",
            format!("[vof] max_courant {c} must be in (0, 1]"),
        ));
    }
    if v.passes() > 10 {
        return Err(bad(
            "vof.smoothing",
            "[vof] smoothing must be at most 10".into(),
        ));
    }
    let a0 = v.initial.unwrap_or(0.0);
    if !(0.0..=1.0).contains(&a0) {
        return Err(bad(
            "vof.initial",
            format!("[vof] initial {a0} must be in [0, 1]"),
        ));
    }
    for (i, r) in v.region.iter().enumerate() {
        if !(0.0..=1.0).contains(&r.alpha) {
            return Err(bad(
                "vof.region",
                format!(
                    "[[vof.region]] {}: alpha {} must be in [0, 1]",
                    i + 1,
                    r.alpha
                ),
            ));
        }
    }
    for (name, a) in &v.inflow {
        if !(0.0..=1.0).contains(a) {
            return Err(bad(
                "vof.inflow",
                format!("[vof.inflow] {name} = {a} must be in [0, 1]"),
            ));
        }
        if !case.boundary.contains_key(name) {
            return Err(bad(
                "vof.inflow",
                format!("[vof.inflow] names the boundary `{name}`, which the case does not have"),
            ));
        }
    }
    // what the first version does not have
    let mut no = Vec::new();
    if case.time.is_none() {
        no.push("a steady case (VOF needs [time])");
    }
    if case.turbulence.model != crate::turbulence::TurbulenceModel::Laminar {
        no.push("a turbulence model");
    }
    if let Some(e) = &case.energy {
        for (i, p) in v.phase.iter().enumerate() {
            let ok = |x: Option<f64>| x.is_some_and(|x| x > 0.0 && x.is_finite());
            if !(ok(p.specific_heat) && ok(p.conductivity)) {
                return Err(bad(
                    "vof.phase",
                    format!(
                        "[[vof.phase]] {} (`{}`): with [energy] each phase needs \
                         specific_heat and conductivity > 0",
                        i + 1,
                        p.name
                    ),
                ));
            }
        }
        if e.density != crate::density::DensityModel::Constant {
            no.push("[energy] density = \"variable\" (the phases have their own densities)");
        }
    }
    if let Some(pc) = &v.phase_change {
        if case.energy.is_none() {
            return Err(bad(
                "vof.phase_change",
                "[vof.phase_change] needs the energy equation ([energy])".into(),
            ));
        }
        if let Some(n) = &pc.liquid
            && !v.phase.iter().any(|p| &p.name == n)
        {
            return Err(bad(
                "vof.phase_change.liquid",
                format!("[vof.phase_change] liquid = `{n}` is not the name of a [[vof.phase]]"),
            ));
        }
        if let Some(p) = pc.pressure
            && foehncraft_kernels::steam::saturation_p(p).is_none()
        {
            return Err(bad(
                "vof.phase_change.pressure",
                format!(
                    "[vof.phase_change] pressure {p} Pa is outside the saturation line of \
                     IAPWS-IF97 (611.2 Pa to 22.064 MPa)"
                ),
            ));
        }
        match pc.saturation() {
            Some((t, h)) if t > 0.0 && t.is_finite() && h > 0.0 && h.is_finite() => {}
            Some(_) => {
                return Err(bad(
                    "vof.phase_change",
                    "[vof.phase_change] saturation_temperature and latent_heat must be > 0".into(),
                ));
            }
            None => {
                return Err(bad(
                    "vof.phase_change",
                    "[vof.phase_change] needs saturation_temperature and latent_heat, or the \
                     pressure of water"
                        .into(),
                ));
            }
        }
        let r = pc.relax();
        if !(r > 0.0 && r.is_finite()) {
            return Err(bad(
                "vof.phase_change.relaxation",
                format!("[vof.phase_change] relaxation {r} must be > 0"),
            ));
        }
    }
    if !case.scalars.is_empty() {
        no.push("passive scalars or species");
    }
    if case.buoyancy.model != crate::buoyancy::Model::None {
        no.push("a [buoyancy] model (gravity acts through [gravity] alone)");
    }
    if !case.regions.is_empty() {
        no.push("regions or several fluids");
    }
    if !case.spray.is_empty() {
        no.push("sprays");
    }
    if !case.immersed.is_empty() {
        no.push("immersed bodies");
    }
    if case.solver.algorithm == Some(crate::simplec::Algorithm::Coupled) {
        no.push("the coupled algorithm (PISO and SIMPLEC only)");
    }
    if !case.mrf.is_empty() || !case.sliding.is_empty() {
        no.push("rotating zones or sliding interfaces");
    }
    if case.mesh_motion.is_some() {
        no.push("a moving or deforming mesh");
    }
    if case.run.restart.is_some() {
        no.push("a restart");
    }
    if case.boundary.values().any(|b| b.is_periodic()) {
        no.push("periodic boundaries");
    }
    if !no.is_empty() {
        return Err(unsupported(format!(
            "[vof] does not work with {} in this version",
            no.join(", ")
        )));
    }
    Ok(())
}

/// The volume-of-fluid state of a flow.
#[derive(Debug)]
pub struct Vof {
    /// The settings.
    pub section: VofSection,
    rho: [f64; 2],
    mu: [f64; 2],
    /// α per cell.
    pub alpha: Buf,
    alpha_start: Buf,
    alpha_new: Buf,
    alpha_f: Buf,
    alpha_f_mean: Buf,
    flux: Buf,
    inv_rho_f: Buf,
    /// ρ per face on the device (the divisor of the volume fluxes).
    pub rho_f: Buf,
    grad: Buf3,
    gnorm: Buf,
    clipped: Buf,
    work: Buf,
    work2: Buf,
    kappa: Buf,
    tension: Buf,
    alpha_b: Buf,
    /// The raw face offsets o_f of gravity and surface tension (device).
    pub offset: Buf,
    /// ρ per cell (device).
    pub rho_c: Buf,
    /// The balanced body force B_P per unit volume (device).
    pub body: Buf3,
    /// Scratch of the volume-based pressure rows: c_f/ρ_f.
    pub cf_vol: Buf,
    /// Scratch: ṁ_f/ρ_f.
    pub m_vol: Buf,
    /// ρ per cell of the last mix (host).
    pub rho_cell: Vec<f64>,
    /// ρ per face of the last mix (host).
    pub rho_face: Vec<f64>,
    /// α per cell of the last mix (host).
    pub alpha_host: Vec<f64>,
    /// α per face of the last mix (host).
    alpha_face: Vec<f64>,
    /// ρ at the time levels n and n − 1.
    old: Option<[Vec<f64>; 2]>,
    saved: Option<Saved>,
    /// σ κ_f (α_N − α_P) per face of the last refresh.
    pub tension_offset: Vec<f64>,
    /// The totals of the run.
    pub totals: Totals,
    volume: Vec<f64>,
    /// c_p and λ of the two phases, with `[energy]` (MP2).
    thermal: Option<([f64; 2], [f64; 2])>,
    /// The phase change (MP2).
    pub pc: Option<PhaseChange>,
    /// The sink rate of α per cell (1/s) of the phase change (device).
    sink: Buf,
    /// Host copies of the connectivity: owner, neighbour, weight.
    owner: Vec<u32>,
    neighbour: Vec<u32>,
    weight: Vec<f64>,
}

/// The state of the phase change (MP2).
#[derive(Debug, Clone)]
pub struct PhaseChange {
    /// T_sat (K).
    pub t_sat: f64,
    /// h_lv (J/kg).
    pub h_lv: f64,
    /// The index of the liquid phase.
    pub liquid: usize,
    relax: f64,
    /// The interface cells of the current step and their G (W/K).
    rows: Vec<(u32, f64)>,
    /// The heat per cell (W) the limit of the last conversion left, given
    /// back in the current step.
    carry: Vec<f64>,
    /// The mass rate per cell (kg/s, evaporation positive) of the
    /// continuity rows of the current step; it leaves α in the next
    /// advance.
    pub used: Vec<f64>,
    /// `used` moved into the pure vapour cells beside the interface (kg/s
    /// per cell): where the vapour volume appears.
    vap: Vec<f64>,
    /// `used` moved into the pure liquid cells beside the interface: where
    /// the liquid volume vanishes.
    liq: Vec<f64>,
    /// Δt of the current step.
    dt: f64,
}

/// What [`Vof::step_back`] restores.
#[derive(Debug)]
struct Saved {
    old: Option<[Vec<f64>; 2]>,
    rho_cell: Vec<f64>,
    rho_face: Vec<f64>,
    alpha: Vec<f64>,
    alpha_face: Vec<f64>,
    totals: Totals,
    pc: Option<PhaseChange>,
}

/// The totals of a run (for `result.json` and `vof.csv`).
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema, PartialEq)]
pub struct Totals {
    /// Volume of phase 1 at the start (m³).
    pub volume_start: f64,
    /// Volume of phase 1 now (m³).
    pub volume: f64,
    /// Volume of phase 1 that left through the boundaries (m³, net).
    pub volume_out: f64,
    /// Volume clipped to keep 0 ≤ α ≤ 1 (m³, signed sum).
    pub clipped_volume: f64,
    /// Sub-steps of the last step.
    pub substeps: usize,
    /// Most sub-steps in a step.
    pub max_substeps: usize,
    /// Smallest and largest α after the last step.
    pub alpha_range: [f64; 2],
    /// Volume of phase 1 the phase change took out (m³, net; MP2).
    #[serde(default)]
    pub volume_phase_change: f64,
    /// The mass evaporated over the run (kg, net; condensation negative).
    #[serde(default)]
    pub evaporated: f64,
    /// The evaporation rate of the current step (kg/s, net).
    #[serde(default)]
    pub phase_rate: f64,
    /// The interface cells held at T_sat in the current step.
    #[serde(default)]
    pub interface_cells: usize,
}

/// The start field of α from `[vof] initial` and the regions.
pub fn start_alpha(case: &crate::case::FlowCase, mesh: &Mesh, geo: &Geometry) -> Result<Vec<f64>> {
    let v = &case.vof;
    let nc = mesh.n_cells;
    let mut a = vec![v.initial.unwrap_or(0.0); nc];
    if v.region.is_empty() {
        return Ok(a);
    }
    let sec = crate::immersed::ImmersedSection {
        // (the permeability is not used: only the cell fractions)
        permeability: Some(1.0),
        samples: Some(6),
        body: v
            .region
            .iter()
            .enumerate()
            .map(|(i, r)| crate::immersed::BodySection {
                name: format!("vof-region-{}", i + 1),
                shape: r.shape,
                centre: r.centre,
                radius: r.radius,
                axis: r.axis,
                normal: None,
                length: r.length,
                lower: r.lower,
                upper: r.upper,
                fraction: r.fraction.clone(),
                file: None,
                scale: None,
                invert: None,
                omega: None,
                patches: None,
            })
            .collect(),
        topology: None,
        ..Default::default()
    };
    let bodies = crate::immersed::bodies_with(&sec, &case.env()?)
        .map_err(|e| bad("vof.region", format!("[[vof.region]]: {e:#}")))?;
    let fr = crate::immersed::fractions(&bodies, 6, mesh, geo);
    for (r, f) in v.region.iter().zip(&fr) {
        for c in 0..nc {
            a[c] = a[c] * (1.0 - f[c]) + r.alpha * f[c];
        }
    }
    Ok(a)
}

impl Vof {
    /// The VOF state of `case` with the start field `alpha` (per cell) on a
    /// flow device; `open` marks the boundary faces whose curvature takes
    /// the owner's normal, `alpha_b` is the inflow α per boundary face.
    pub fn new(
        section: VofSection,
        dev: &mut dyn FlowDevice,
        disc: &foehncraft_fv::disc::Discretisation,
        alpha: &[f64],
        open: &[u32],
        alpha_b: &[f64],
    ) -> Result<Self> {
        let (nc, nf) = (disc.n_cells, disc.n_faces);
        let cell = |dev: &mut dyn FlowDevice| dev.alloc(nc);
        let face = |dev: &mut dyn FlowDevice| dev.alloc(nf);
        let a = cell(dev)?;
        dev.upload(a, alpha)?;
        let grad = Buf3 {
            x: cell(dev)?,
            y: cell(dev)?,
            z: cell(dev)?,
        };
        let ab = dev.alloc(nf - disc.n_interior)?;
        dev.upload(ab, alpha_b)?;
        dev.vof_open(open)?;
        let p = section.phase.clone();
        let rho = [p[0].density, p[1].density];
        let mu = [p[0].viscosity, p[1].viscosity];
        let mut v = Self {
            rho,
            mu,
            alpha: a,
            alpha_start: cell(dev)?,
            alpha_new: cell(dev)?,
            alpha_f: face(dev)?,
            alpha_f_mean: face(dev)?,
            flux: face(dev)?,
            inv_rho_f: face(dev)?,
            rho_f: face(dev)?,
            grad,
            gnorm: cell(dev)?,
            clipped: cell(dev)?,
            work: cell(dev)?,
            work2: cell(dev)?,
            kappa: cell(dev)?,
            tension: face(dev)?,
            offset: face(dev)?,
            rho_c: cell(dev)?,
            body: Buf3 {
                x: cell(dev)?,
                y: cell(dev)?,
                z: cell(dev)?,
            },
            alpha_b: ab,
            cf_vol: face(dev)?,
            m_vol: face(dev)?,
            rho_cell: Vec::new(),
            rho_face: Vec::new(),
            alpha_host: alpha.to_vec(),
            alpha_face: Vec::new(),
            old: None,
            saved: None,
            tension_offset: vec![0.0; nf],
            totals: Totals::default(),
            volume: disc.cell_volume[..nc].to_vec(),
            thermal: None,
            pc: None,
            sink: cell(dev)?,
            owner: disc.owner[..nf].to_vec(),
            neighbour: disc.neighbour[..disc.n_interior].to_vec(),
            weight: disc.weight[..nf].to_vec(),
            section,
        };
        if p.iter()
            .all(|x| x.specific_heat.is_some() && x.conductivity.is_some())
        {
            let c = |i: usize| p[i].specific_heat.unwrap_or(1.0);
            let k = |i: usize| p[i].conductivity.unwrap_or(1.0);
            v.thermal = Some(([c(0), c(1)], [k(0), k(1)]));
        }
        if let Some((t_sat, h_lv)) = v.section.phase_change.as_ref().and_then(|q| q.saturation()) {
            let q = v.section.phase_change.as_ref().expect("phase change");
            v.pc = Some(PhaseChange {
                t_sat,
                h_lv,
                liquid: v.section.liquid(),
                relax: q.relax(),
                rows: Vec::new(),
                carry: vec![0.0; nc],
                used: vec![0.0; nc],
                vap: vec![0.0; nc],
                liq: vec![0.0; nc],
                dt: 0.0,
            });
        }
        // the start face values: linear between the cells
        let af = linear_faces(disc, alpha);
        dev.upload(v.alpha_f_mean, &af)?;
        v.mix(alpha, &af);
        v.totals.volume_start = v.phase_volume(alpha);
        v.totals.volume = v.totals.volume_start;
        v.totals.alpha_range = range(alpha);
        Ok(v)
    }

    /// Sets the volume fluxes F_f (m³/s) of the next advance (a prescribed
    /// velocity: the advection checks).
    pub fn set_flux(&mut self, dev: &mut dyn FlowDevice, flux: &[f64]) -> Result<()> {
        dev.upload(self.flux, flux)
    }

    /// A face buffer for the density ratio of a step start (free outside
    /// the pressure rows).
    pub fn ratio_buf(&self) -> Buf {
        self.cf_vol
    }

    /// The buffer of the volume fluxes.
    pub fn flux_buf(&self) -> Buf {
        self.flux
    }

    /// α per face of the last mix (host).
    pub fn face_alpha(&self) -> Vec<f64> {
        self.alpha_face.clone()
    }

    /// The densities and viscosities of phase 1 and 2.
    pub fn properties(&self) -> ([f64; 2], [f64; 2]) {
        (self.rho, self.mu)
    }

    fn phase_volume(&self, a: &[f64]) -> f64 {
        a.iter().zip(&self.volume).map(|(a, v)| a * v).exact_sum()
    }

    /// ρ and μ from α per cell and per face (host).
    fn mix(&mut self, a: &[f64], af: &[f64]) {
        let [r1, r2] = self.rho;
        self.rho_cell = a.iter().map(|&x| x * r1 + (1.0 - x) * r2).collect();
        self.rho_face = af.iter().map(|&x| x * r1 + (1.0 - x) * r2).collect();
        self.alpha_host = a.to_vec();
        self.alpha_face = af.to_vec();
    }

    /// μ per face of the last mix.
    pub fn mu_face(&self, af: &[f64]) -> Vec<f64> {
        let [m1, m2] = self.mu;
        af.iter().map(|&x| x * m1 + (1.0 - x) * m2).collect()
    }

    /// Advances α over a step `dt` with the volume fluxes of the current
    /// mass flux `mass_flux` (divided by the face densities of the last
    /// mix): sub-steps, clipping, the mean face values. Returns the face
    /// values α_f (mean over the sub-steps).
    pub fn advance(
        &mut self,
        dev: &mut dyn FlowDevice,
        mass_flux: Buf,
        dt: f64,
        boundary_start: usize,
    ) -> Result<Vec<f64>> {
        dev.quotient(self.flux, 1.0, mass_flux, self.rho_f)?;
        self.advance_flux(dev, dt, boundary_start)
    }

    /// [`Vof::advance`] with the volume fluxes already in their buffer
    /// ([`Vof::set_flux`]).
    pub fn advance_flux(
        &mut self,
        dev: &mut dyn FlowDevice,
        dt: f64,
        boundary_start: usize,
    ) -> Result<Vec<f64>> {
        dev.copy(self.alpha, self.alpha_start)?;
        // the phase change of the step before (MP2): ṁ/ρ_l of liquid leaves
        // α (or ṁ/ρ_v of vapour joins it, if phase 1 is the vapour)
        if let Some(pc) = &self.pc {
            // (the liquid leaves in the pure liquid cells beside the
            // interface, the vapour joins in the pure vapour cells: the
            // rates moved off the interface, `spread`)
            let l = pc.liquid;
            let (r, m) = if l == 0 {
                (self.rho[0], &pc.liq)
            } else {
                (-self.rho[1], &pc.vap)
            };
            let s: Vec<f64> = m
                .iter()
                .zip(&self.volume)
                .map(|(m, v)| m / (r * v))
                .collect();
            dev.upload(self.sink, &s)?;
            let removed: f64 = s.iter().zip(&self.volume).map(|(s, v)| s * v).exact_sum();
            self.totals.volume_phase_change += removed * dt;
            let net: f64 = pc.used.iter().copied().exact_sum();
            self.totals.evaporated += net * dt;
            self.totals.phase_rate = net;
        }
        dev.vof(
            VofOp::OutflowRate,
            &VofBufs {
                flux: Some(self.flux),
                out: Some(self.work),
                ..Default::default()
            },
        )?;
        let rate = dev.min_max(self.work)?.1;
        let co = rate * dt;
        let n = if co.is_finite() && co > 0.0 {
            ((co / self.section.courant()).ceil() as usize).max(1)
        } else {
            1
        };
        anyhow::ensure!(
            n <= 100_000,
            "[vof]: {n} sub-steps of α in one time step (outflow Courant number {co:.3e}); \
             the time step is far too long for the velocity"
        );
        let h = dt / n as f64;
        dev.fill(self.alpha_f_mean, 0.0)?;
        let mut clip = 0.0;
        for _ in 0..n {
            dev.vof(
                VofOp::Gradient,
                &VofBufs {
                    alpha: Some(self.alpha),
                    out: Some(self.gnorm),
                    out3: Some(self.grad),
                    ..Default::default()
                },
            )?;
            dev.vof(
                VofOp::Faces,
                &VofBufs {
                    flux: Some(self.flux),
                    alpha: Some(self.alpha),
                    grad: Some(self.grad),
                    alpha_b: Some(self.alpha_b),
                    scheme: self.section.scheme.code(),
                    dt: h,
                    beta: self.section.beta(),
                    out: Some(self.alpha_f),
                    ..Default::default()
                },
            )?;
            dev.vof(
                VofOp::Update,
                &VofBufs {
                    flux: Some(self.flux),
                    alpha_f: Some(self.alpha_f),
                    alpha: Some(self.alpha),
                    kappa: self.pc.as_ref().map(|_| self.sink),
                    dt: h,
                    out: Some(self.alpha_new),
                    out2: Some(self.clipped),
                    ..Default::default()
                },
            )?;
            dev.copy(self.alpha_new, self.alpha)?;
            dev.lincomb(
                self.alpha_f_mean,
                1.0,
                self.alpha_f_mean,
                1.0 / n as f64,
                self.alpha_f,
            )?;
            let c = dev.download(self.clipped)?;
            clip += c.iter().copied().exact_sum();
            // the volume of phase 1 through the boundary faces
            let af = dev.download(self.alpha_f)?;
            let fl = dev.download(self.flux)?;
            self.totals.volume_out += (boundary_start..fl.len())
                .map(|f| fl[f] * af[f] * h)
                .exact_sum();
        }
        self.totals.clipped_volume += clip;
        self.totals.substeps = n;
        self.totals.max_substeps = self.totals.max_substeps.max(n);
        let a = dev.download(self.alpha)?;
        let af = dev.download(self.alpha_f_mean)?;
        self.totals.volume = self.phase_volume(&a);
        self.totals.alpha_range = range(&a);
        self.mix(&a, &af);
        Ok(af)
    }

    /// Starts a step: the density of now becomes level n (level n − 1 the
    /// one before); keeps what [`Vof::step_back`] restores.
    pub fn start_step(&mut self) {
        let prev = self.old.take();
        self.saved = Some(Saved {
            old: prev.clone(),
            rho_cell: self.rho_cell.clone(),
            rho_face: self.rho_face.clone(),
            alpha: self.alpha_host.clone(),
            alpha_face: self.alpha_face.clone(),
            totals: self.totals.clone(),
            pc: self.pc.clone(),
        });
        let n1 = match prev {
            Some([n, _]) => n,
            None => self.rho_cell.clone(),
        };
        self.old = Some([self.rho_cell.clone(), n1]);
    }

    /// The old density levels of the time term (n, n − 1).
    pub fn old_levels(&self) -> Option<&[Vec<f64>; 2]> {
        self.old.as_ref()
    }

    /// Takes the last step back (a rejected step): α, the densities and
    /// the levels as before [`Vof::start_step`].
    pub fn step_back(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        if let Some(sv) = self.saved.take() {
            self.old = sv.old;
            self.rho_cell = sv.rho_cell;
            self.rho_face = sv.rho_face;
            self.alpha_host = sv.alpha;
            self.alpha_face = sv.alpha_face;
            self.totals = sv.totals;
            self.pc = sv.pc;
            dev.copy(self.alpha_start, self.alpha)?;
            // the volume fluxes of the start of the step (ṁ/ρ_f of the
            // restored densities follow in `upload_faces`)
        }
        Ok(())
    }

    /// Uploads 1/ρ_f and ρ_f of the last mix and rescales the mass flux to
    /// the new face densities: ṁ_f = ρ_f F_f with the volume flux of the
    /// advection.
    pub fn upload_faces(&mut self, dev: &mut dyn FlowDevice, mass_flux: Option<Buf>) -> Result<()> {
        let inv: Vec<f64> = self.rho_face.iter().map(|r| 1.0 / r).collect();
        dev.upload(self.inv_rho_f, &inv)?;
        dev.upload(self.rho_f, &self.rho_face)?;
        if let Some(m) = mass_flux {
            dev.quotient(m, 1.0, self.flux, self.inv_rho_f)?;
        }
        Ok(())
    }

    /// The surface-tension offsets σ κ_f (α_N − α_P) per face from the
    /// current α (zero for σ = 0), into [`Vof::tension_offset`].
    pub fn refresh_tension(&mut self, dev: &mut dyn FlowDevice, delta: f64) -> Result<()> {
        let sigma = self.section.sigma();
        if sigma == 0.0 {
            return Ok(());
        }
        // α̃: `smoothing` passes, alternating between two buffers
        let mut src = self.alpha;
        let bufs = [self.work, self.work2];
        for i in 0..self.section.passes() as usize {
            let dst = bufs[i % 2];
            dev.vof(
                VofOp::Smooth,
                &VofBufs {
                    alpha: Some(src),
                    out: Some(dst),
                    ..Default::default()
                },
            )?;
            src = dst;
        }
        dev.vof(
            VofOp::Gradient,
            &VofBufs {
                alpha: Some(src),
                out: Some(self.gnorm),
                out3: Some(self.grad),
                ..Default::default()
            },
        )?;
        dev.vof(
            VofOp::Curvature,
            &VofBufs {
                grad: Some(self.grad),
                s: delta,
                out: Some(self.kappa),
                ..Default::default()
            },
        )?;
        dev.vof(
            VofOp::Tension,
            &VofBufs {
                alpha: Some(self.alpha),
                kappa: Some(self.kappa),
                gnorm: Some(self.gnorm),
                s: sigma,
                out: Some(self.tension),
                ..Default::default()
            },
        )?;
        self.tension_offset = dev.download(self.tension)?;
        Ok(())
    }

    /// The curvature per cell of the last [`Vof::refresh_tension`].
    pub fn curvature(&self, dev: &mut dyn FlowDevice) -> Result<Vec<f64>> {
        dev.download(self.kappa)
    }
}

impl Vof {
    /// The thermal mixture of the last mix (MP2), with `[energy]`: ρ per
    /// cell and the mixture of the energy rows: c_p per cell
    /// ((ρ c_p)/ρ) and per face (of the advected face α, so that ṁ_f c_p,f
    /// = F_f (ρ c_p)_f), λ per face (the harmonic mean of the cells' series
    /// conductivities) and μ per face.
    pub fn thermal_mixture(&self) -> Option<(Vec<f64>, crate::scalar::MixtureProps)> {
        let (cp, k) = self.thermal?;
        let [r1, r2] = self.rho;
        let rc = |a: f64| a * r1 * cp[0] + (1.0 - a) * r2 * cp[1];
        let cp_cell: Vec<f64> = self
            .alpha_host
            .iter()
            .zip(&self.rho_cell)
            .map(|(&a, r)| rc(a) / r)
            .collect();
        let cp_face: Vec<f64> = self
            .alpha_face
            .iter()
            .zip(&self.rho_face)
            .map(|(&a, r)| rc(a) / r)
            .collect();
        let kc: Vec<f64> = self
            .alpha_host
            .iter()
            .map(|&a| 1.0 / (a / k[0] + (1.0 - a) / k[1]))
            .collect();
        let ni = self.neighbour.len();
        let k_face: Vec<f64> = (0..self.owner.len())
            .map(|f| {
                let o = self.owner[f] as usize;
                if f < ni {
                    let n = self.neighbour[f] as usize;
                    let w = self.weight[f];
                    1.0 / (w / kc[o] + (1.0 - w) / kc[n])
                } else {
                    kc[o]
                }
            })
            .collect();
        let mu_face = self.mu_face(&self.alpha_face);
        Some((
            self.rho_cell.clone(),
            crate::scalar::MixtureProps {
                cp_cell,
                cp_face,
                k_face,
                mu_face,
            },
        ))
    }

    /// ρ c_p per cell of the last mix (MP2).
    fn capacity(&self) -> Option<Vec<f64>> {
        let (cp, _) = self.thermal?;
        let [r1, r2] = self.rho;
        Some(
            self.alpha_host
                .iter()
                .map(|&a| a * r1 * cp[0] + (1.0 - a) * r2 * cp[1])
                .collect(),
        )
    }

    /// The liquid fraction per cell of the last mix.
    fn liquid_fraction(&self, liquid: usize) -> Vec<f64> {
        self.alpha_host
            .iter()
            .map(|&a| if liquid == 0 { a } else { 1.0 - a })
            .collect()
    }

    /// The start of a step `dt` with the phase change (MP2), before the
    /// advance: the heat G (T − T_sat) the interface cells gave up in the
    /// last step at the temperature `t` (per cell) becomes the mass rate
    /// ṁ = G (T − T_sat)/h_lv, limited to the liquid (evaporation) or
    /// vapour (condensation) the cell holds over `dt`; the heat the limit
    /// leaves is carried into this step. Returns the rates (kg/s per cell)
    /// of this step's continuity rows.
    pub fn phase_convert(&mut self, t: &[f64], dt: f64) -> Option<Vec<f64>> {
        let liquid = self.pc.as_ref()?.liquid;
        let al = self.liquid_fraction(liquid);
        let (rl, rv) = (self.rho[liquid], self.rho[1 - liquid]);
        let pc = self.pc.as_mut()?;
        let nc = self.volume.len();
        let mut next = vec![0.0; nc];
        let mut carry = vec![0.0; nc];
        let scale = if pc.dt > 0.0 { pc.dt / dt } else { 1.0 };
        for &(c, g) in &pc.rows {
            let c = c as usize;
            let heat = g * (t[c] - pc.t_sat);
            let m = heat / pc.h_lv;
            let v = self.volume[c];
            let hi = al[c] * rl * v / dt;
            let lo = -(1.0 - al[c]) * rv * v / dt;
            let ml = m.clamp(lo, hi);
            next[c] = ml;
            // the heat of the last step over this one's duration
            carry[c] = (heat - ml * pc.h_lv) * scale;
        }
        pc.carry = carry;
        Some(next)
    }

    /// After the advance of a step `dt` (MP2): the rates `next` become those
    /// of the continuity rows, and the interface cells of the new α at the
    /// temperature `t` (per cell) take G = `relaxation` ρ c_p V/Δt.
    pub fn phase_rows(&mut self, t: &[f64], next: Vec<f64>, dt: f64) {
        let Some(liquid) = self.pc.as_ref().map(|p| p.liquid) else {
            return;
        };
        let al = self.liquid_fraction(liquid);
        let cap = self.capacity().unwrap_or_else(|| vec![0.0; al.len()]);
        let Some(pc) = self.pc.as_mut() else {
            return;
        };
        const EPS: f64 = 1e-3;
        let nc = al.len();
        let mut on = vec![false; nc];
        for (c, &a) in al.iter().enumerate() {
            on[c] = a > EPS && a < 1.0 - EPS;
        }
        // a sharp face: the liquid side when it is above T_sat, the vapour
        // side below (the side that can change phase)
        for f in 0..self.neighbour.len() {
            let (o, n) = (self.owner[f] as usize, self.neighbour[f] as usize);
            let (lo, vn) = (al[o] >= 0.5, al[n] >= 0.5);
            if lo == vn {
                continue;
            }
            let (l, v) = if lo { (o, n) } else { (n, o) };
            if t[l] >= pc.t_sat {
                on[l] = true;
            }
            if t[v] < pc.t_sat {
                on[v] = true;
            }
        }
        pc.rows = (0..nc)
            .filter(|&c| on[c])
            .map(|c| (c as u32, pc.relax * cap[c] * self.volume[c] / dt))
            .collect();
        // carried heat of cells that are no interface cells any more
        for (c, q) in pc.carry.iter().enumerate() {
            if *q != 0.0 && !on[c] {
                pc.rows.push((c as u32, 0.0));
            }
        }
        pc.rows.sort_unstable_by_key(|r| r.0);
        pc.rows.dedup_by_key(|r| r.0);
        pc.vap = spread(&next, &al, false, &self.owner, &self.neighbour);
        pc.liq = spread(&next, &al, true, &self.owner, &self.neighbour);
        pc.used = next;
        pc.dt = dt;
        self.totals.interface_cells = on.iter().filter(|&&x| x).count();
    }

    /// The energy rows of the phase change (MP2): per interface cell (cell,
    /// Q + G T_sat, G) in W and W/K, Q the carried heat.
    pub fn energy_rows(&self) -> Option<Vec<(u32, f64, f64)>> {
        let pc = self.pc.as_ref()?;
        Some(
            pc.rows
                .iter()
                .map(|&(c, g)| (c, pc.carry[c as usize] + g * pc.t_sat, g))
                .collect(),
        )
    }

    /// The source of the continuity rows of the volume (m³/s per cell, in
    /// the form V ∂ρ/∂t − S of the variable density): −(ṁ_v/ρ_v − ṁ_l/ρ_l)
    /// with the rates moved off the interface ([`spread`]).
    pub fn volume_source(&self) -> Option<Vec<f64>> {
        let pc = self.pc.as_ref()?;
        let (rl, rv) = (self.rho[pc.liquid], self.rho[1 - pc.liquid]);
        Some(
            pc.vap
                .iter()
                .zip(&pc.liq)
                .map(|(v, l)| -(v / rv - l / rl))
                .collect(),
        )
    }

    /// The cap of the pressure reconstruction ([`crate::device::VofOp::
    /// Reconstruct`]): 1 with a phase change (MP2), 0 without (MP1 rows
    /// unchanged).
    pub fn reconstruct_cap(&self) -> f64 {
        if self.pc.is_some() { 1.0 } else { 0.0 }
    }

    /// The phase-change rate per cell of the current step (kg/(m³ s),
    /// evaporation positive), for the field output.
    pub fn phase_rate_density(&self) -> Option<Vec<f64>> {
        let pc = self.pc.as_ref()?;
        Some(
            pc.used
                .iter()
                .zip(&self.volume)
                .map(|(m, v)| m / v)
                .collect(),
        )
    }
}

/// Moves the phase-change rates `m` (per cell) off the interface into the
/// pure cells of one phase beside it (MP2): the liquid (`liquid`, α_l >
/// 0.99) or the vapour (α_l < 0.01), with the liquid fraction `al` per cell.
/// A cell that holds a rate and is not pure passes it on to its face
/// neighbours in proportion to their fraction of the phase, up to eight
/// times; what finds no pure cell stays where it is. The sum is kept to
/// round-off. With the volume sources in pure cells, the pressure rows
/// see no source inside the jump of the density (a source in the
/// interface cell of a water–steam interface drives a checkerboard of the
/// velocity), and α stays at 0 and 1 there; the interface moves by the
/// flux of the vapour into it. The shift is that of the smeared sources
/// of Hardt & Wondra (2008, J. Comput. Phys. 227, 5871–5895), local here
/// instead of normalised over the domain.
pub fn spread(m: &[f64], al: &[f64], liquid: bool, owner: &[u32], neighbour: &[u32]) -> Vec<f64> {
    let nc = m.len();
    let frac = |c: usize| if liquid { al[c] } else { 1.0 - al[c] };
    let pure = |c: usize| frac(c) > 0.99;
    let mut a = m.to_vec();
    if a.iter().all(|&x| x == 0.0) {
        return a;
    }
    // the face neighbours of each cell
    let mut ptr = vec![0usize; nc + 1];
    for f in 0..neighbour.len() {
        ptr[owner[f] as usize + 1] += 1;
        ptr[neighbour[f] as usize + 1] += 1;
    }
    for c in 0..nc {
        ptr[c + 1] += ptr[c];
    }
    let mut nb = vec![0u32; ptr[nc]];
    let mut fill = ptr.clone();
    for f in 0..neighbour.len() {
        let (o, n) = (owner[f] as usize, neighbour[f] as usize);
        nb[fill[o]] = n as u32;
        fill[o] += 1;
        nb[fill[n]] = o as u32;
        fill[n] += 1;
    }
    for _ in 0..8 {
        let mut b = vec![0.0; nc];
        let mut moved = false;
        for c in 0..nc {
            if a[c] == 0.0 {
                continue;
            }
            if pure(c) {
                b[c] += a[c];
                continue;
            }
            let w: f64 = nb[ptr[c]..ptr[c + 1]]
                .iter()
                .map(|&n| frac(n as usize))
                .sum();
            if w <= 0.0 {
                b[c] += a[c];
                continue;
            }
            for &n in &nb[ptr[c]..ptr[c + 1]] {
                b[n as usize] += a[c] * frac(n as usize) / w;
            }
            moved = true;
        }
        a = b;
        if !moved {
            break;
        }
    }
    a
}

fn range(a: &[f64]) -> [f64; 2] {
    let lo = a.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = a.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    [lo, hi]
}

/// The linear face values of a cell field (the owner value on boundary
/// faces).
pub fn linear_faces(disc: &foehncraft_fv::disc::Discretisation, a: &[f64]) -> Vec<f64> {
    (0..disc.n_faces)
        .map(|f| {
            let p = disc.owner[f] as usize;
            if f < disc.n_interior {
                let n = disc.neighbour[f] as usize;
                let w = disc.weight[f];
                w * a[p] + (1.0 - w) * a[n]
            } else {
                a[p]
            }
        })
        .collect()
}

/// The open-boundary flags (1: inlets and outlets, 0: walls and symmetry)
/// and the inflow α per boundary face of a case.
pub fn boundary_arrays(
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    n_interior: usize,
    n_faces: usize,
) -> (Vec<u32>, Vec<f64>) {
    let ni = n_interior;
    let nb = n_faces - ni;
    let mut open = vec![0u32; nb];
    let mut ab = vec![0.0; nb];
    for b in &mesh.boundaries {
        let a = case.vof.inflow.get(&b.name).copied().unwrap_or(0.0);
        let o = case.boundary.get(&b.name).is_some_and(|c| c.is_open());
        for f in b.range() {
            if f >= ni && f - ni < nb {
                ab[f - ni] = a;
                open[f - ni] = u32::from(o);
            }
        }
    }
    (open, ab)
}

/// Attaches `[vof]` of `case` to the flow (after the start fields; nothing
/// without `[vof]`).
pub fn attach(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
) -> Result<()> {
    if case.vof.is_empty() {
        return Ok(());
    }
    let alpha = start_alpha(case, mesh, geo)?;
    let (open, ab) = boundary_arrays(case, mesh, fl.disc.n_interior, fl.disc.n_faces);
    let gravity = crate::buoyancy::hydrostatic(case).map(|h| (h.g, h.rho_ref));
    fl.attach_vof(case.vof.clone(), &alpha, &open, &ab, gravity)
}

/// The result of [`advect`].
#[derive(Debug, Clone)]
pub struct Advected {
    /// The cell volumes.
    pub volume: Vec<f64>,
    /// α at the start.
    pub start: Vec<f64>,
    /// α at the end.
    pub end: Vec<f64>,
    /// The totals (volume, clipped volume, sub-steps).
    pub totals: Totals,
}

/// The advection of the start α of `case` alone, in the prescribed
/// velocity `u(x, t)` (m/s) for `steps` steps of `dt`, on the flow device
/// of the case (the kernels of the run; verification: Zalesak's disk, the
/// shear reversal). The face fluxes are u(x_f, t + Δt/2)·S_f.
pub fn advect(
    case: &crate::case::FlowCase,
    u: &dyn Fn([f64; 3], f64) -> [f64; 3],
    dt: f64,
    steps: usize,
    threads: Option<usize>,
) -> Result<Advected> {
    let (mesh, geo, _) = crate::case::mesh_for(case, None)?;
    let (mesh, geo) = (&mesh, &geo);
    let phys = crate::simplec::Physics::constant(case.rho(), case.mu(), [0.0; 3], true);
    let arrays = crate::bc::arrays(
        mesh,
        geo,
        &case.boundary,
        phys.rho,
        crate::bc::At::numbers(None),
    )?;
    let disc = crate::simplec::discretisation(mesh, geo, &arrays)?;
    let mut fl = crate::simplec::Flow::new(
        disc,
        arrays,
        phys,
        case.solver.clone(),
        [0.0; 3],
        0.0,
        Some(0),
        &crate::simplec::Target {
            device: "cpu".into(),
            threads,
            ptx: Default::default(),
            ranks: None,
        },
    )?;
    let alpha = start_alpha(case, mesh, geo)?;
    let (open, ab) = boundary_arrays(case, mesh, fl.disc.n_interior, fl.disc.n_faces);
    let mut v = Vof::new(
        case.vof.clone(),
        fl.dev.as_mut(),
        &fl.disc,
        &alpha,
        &open,
        &ab,
    )?;
    let nf = fl.disc.n_faces;
    let ni = fl.disc.n_interior;
    let mut t = 0.0;
    for _ in 0..steps {
        let tm = t + 0.5 * dt;
        let flux: Vec<f64> = (0..nf)
            .map(|f| {
                let x = geo.face_centroid[f];
                let s = geo.face_area[f];
                let w = u(x, tm);
                w[0] * s[0] + w[1] * s[1] + w[2] * s[2]
            })
            .collect();
        v.set_flux(fl.dev.as_mut(), &flux)?;
        v.advance_flux(fl.dev.as_mut(), dt, ni)?;
        t += dt;
    }
    Ok(Advected {
        volume: geo.cell_volume[..mesh.n_cells].to_vec(),
        start: alpha,
        end: v.alpha_host.clone(),
        totals: v.totals.clone(),
    })
}

/// Units of `vof.csv`.
pub const CSV_UNITS: &str = "step -, time s, volume m3 (phase 1), volume_out m3 (net through the \
    boundaries), clipped m3, substeps -, alpha_min alpha_max -, x_min..z_max m (extent of the \
    cells with alpha > 0.5), evaporated kg (net, MP2), phase_rate kg/s, interface_cells -";

/// The header of `vof.csv`.
pub const CSV_HEAD: &str = "step,time,volume,volume_out,clipped,substeps,alpha_min,alpha_max,\
    x_min,x_max,y_min,y_max,z_min,z_max,evaporated,phase_rate,interface_cells";

impl Vof {
    /// One row of `vof.csv` after step `n` at `time`, with the cell
    /// centroids `x`.
    pub fn csv_row(&self, n: usize, time: f64, x: &[[f64; 3]]) -> String {
        let t = &self.totals;
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for (c, &a) in self.alpha_host.iter().enumerate() {
            if a > 0.5 {
                for i in 0..3 {
                    lo[i] = lo[i].min(x[c][i]);
                    hi[i] = hi[i].max(x[c][i]);
                }
            }
        }
        format!(
            "{n},{time:e},{:e},{:e},{:e},{},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{}",
            t.volume,
            t.volume_out,
            t.clipped_volume,
            t.substeps,
            t.alpha_range[0],
            t.alpha_range[1],
            lo[0],
            hi[0],
            lo[1],
            hi[1],
            lo[2],
            hi[2],
            t.evaporated,
            t.phase_rate,
            t.interface_cells
        )
    }
}
