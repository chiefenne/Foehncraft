//! Passive scalars (D51) on the generic transport equation of D116
//! (stage 6, PR 6a):
//!
//! ∂(ρφ)/∂t + ∇·(ṁ φ) − φ ∇·ṁ − ∇·(ρΓ ∇φ) = ρ S,
//!
//! carried by the face mass fluxes of the flow solver and solved after its
//! pressure correction: once per outer iteration of a steady run, once per
//! outer iteration of a time step (with BDF1 or BDF2; a Crank–Nicolson run
//! takes BDF2 for the scalars, as the stage 6 specification sets for every
//! bounded variable). A temperature with constant properties is such a
//! scalar, with Γ = ν/Pr = λ/(ρ c_p).
//!
//! Case file:
//!
//! ```toml
//! [scalars.temperature]
//! diffusivity = 2.1e-5          # Γ (m²/s): number or expression of x, y, z
//! initial = 300.0               # start value (number or expression)
//! source = 0.0                  # S ([φ]/s): number or expression of x, y, z, t
//! bounds = [0.0, 1.0]           # optional: bounded BDF2 (D28 (1))
//!
//! [boundary.wall.scalar.temperature]   # per-physics sub-table (D116)
//! flux = 1.0e-3                 # diffusive flux into the domain ([φ] m/s)
//! # value = 310.0               # or a fixed value
//!
//! [solver.scalar]               # one block for every scalar
//! scheme = "linear-upwind"
//! limiter = "venkatakrishnan"
//! ```
//!
//! Boundaries without a sub-table: walls, slip walls, symmetry planes and
//! wedges have zero flux, inlets the value 0, outlets a zero normal
//! gradient (the backflow carries the owner value), periodic pairs carry
//! the scalar across unchanged. An inlet takes `value`, a wall-like
//! boundary `value` or `flux`, an outlet neither, one boundary of a
//! periodic pair `periodic = "heat-flux"` (the jump Δ = Q/ṁ of D109 (1))
//! with its `bulk` value.
//!
//! The device work lives in [`Scalars`]: one set of buffers per scalar
//! (value, gradient, limiter, the boundary coefficients φ_f = a φ_P + b of
//! `foehncraft_kernels::fv`, the least-squares weights of its own boundary
//! kinds, the face diffusivity, the source) and one set of work buffers
//! shared by all of them.

use crate::bc::BoundaryCondition;
use crate::device::{
    Buf, Buf3, Eps2, FlowDevice, MlpOps, ScalarFaceArgs, ScalarRowArgs, ScalarTime, ScalarTurbOps,
    ScalarWallOps,
};
use crate::diag::{self, Diagnostic};
use crate::expr::{self, Env, Scalar, Value, Vars};
use crate::limiter::{EpsilonScaling, LimiterFunction};
use anyhow::Result;
use foehncraft_fv::disc::Discretisation;
use foehncraft_fv::problem::BoundaryKind;
use foehncraft_fv::solver::{LimiterKind, Preconditioner, Scheme};
use foehncraft_kernels::turb::thermal;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// `[scalars.<name>]`: a passive scalar.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ScalarSection {
    /// Molecular diffusivity Γ (m²/s): ν/Pr for a temperature, ν/Sc for a
    /// concentration; a number or an expression of x, y, z (D53), evaluated
    /// at the face centroids.
    pub diffusivity: Scalar,
    /// Start value: a number or an expression of x, y, z (D53).
    #[serde(default)]
    pub initial: Scalar,
    /// Source per unit mass S (the unit of the scalar per second): a number
    /// or an expression of x, y, z and, in a transient case, t (D53);
    /// explicit in the equation.
    #[serde(default)]
    pub source: Scalar,
    /// Lower and upper bound of the scalar (a mixture fraction: [0, 1]).
    /// The second-order time integration falls back to first order in the
    /// cells where it would leave them (bounded BDF2, D28 (1)), and the log
    /// reports the extremes against them.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bounds: Option<[f64; 2]>,
    /// The exact solution, for the verification (an expression of x, y, z
    /// and, in a transient case, t): the report gives the L2 and largest
    /// error at the cell centroids, a mesh sweep the observed order
    /// (`[verification] scalar_order`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exact: Option<Scalar>,
    /// Turbulent Prandtl (Schmidt) number Pr_t of the `constant` model; the
    /// turbulent diffusivity is ν_t/Pr_t (D135). Used with a turbulence
    /// model only.
    #[serde(default = "default_prt")]
    pub turbulent_prandtl: f64,
    /// Model of the turbulent Prandtl number: `constant` (the number
    /// above) or `kays-crawford` (Pr_t from ν_t/ν and Pr, Pr_t∞ = 0.85).
    #[serde(default)]
    pub turbulent_prandtl_model: PrandtlModel,
    /// Thermal wall function at the walls of a turbulent flow: `kader`
    /// (all y⁺ and any Pr; it reduces to the conduction law Pr y⁺ where
    /// the wall is resolved), `jayatilleke` (log region, blended) or `none`
    /// (the molecular diffusivity at the wall).
    #[serde(default)]
    pub wall_function: ThermalWallFunction,
}

fn default_prt() -> f64 {
    0.85
}

/// The model of the turbulent Prandtl number (D135).
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum PrandtlModel {
    /// A constant `turbulent_prandtl`.
    #[default]
    Constant,
    /// Kays & Crawford (1993): Pr_t falls from 2 Pr_t∞ at ν_t → 0 to
    /// Pr_t∞ in the turbulent core.
    KaysCrawford,
}

/// The thermal wall function of a scalar (D135).
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum ThermalWallFunction {
    /// Kader (1981), all y⁺, any Pr.
    #[default]
    Kader,
    /// Jayatilleke (1969), log region, blended to the conduction limit.
    Jayatilleke,
    /// No wall function: the molecular diffusivity at the wall.
    None,
}

/// The faces and buffers of the flow a turbulent scalar couples to.
#[derive(Debug, Clone, Copy)]
pub struct TurbCoupling {
    /// Laminar viscosity per face.
    pub mu: Buf,
    /// The face viscosity of the flow (μ + ρ ν_t).
    pub mu_eff: Buf,
    /// 1 on the wall faces (per boundary face).
    pub wall: Buf,
    /// The y⁺ of the thermal law per wall face (y* = C_μ^¼ k_P^½ y_P/ν of
    /// the k-ε family and ζ-f/EB, y⁺ of SST), if the flow has walls.
    pub y_plus: Option<Buf>,
    /// The thermal law is taken at ỹ⁺ = max(y⁺, `y_floor`): 11.06 with the
    /// scalable momentum treatment, 0 otherwise.
    pub y_floor: f64,
}

/// The turbulent part of one scalar on the device (D135).
#[derive(Debug, Clone)]
struct TurbField {
    /// ρ Γ per face, molecular (the working `gamma` follows ν_t).
    gamma0: Buf,
    /// q d_n per boundary face.
    fq: Buf,
    law: u32,
    prt_model: u32,
    prt: f64,
    /// d_n per boundary face.
    dn: Vec<f64>,
}

/// `[boundary.<name>.scalar.<scalar>]`: the condition of a scalar on a
/// boundary (D116). One of the two keys.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ScalarBoundary {
    /// Fixed value of the scalar on the boundary: a number or an expression
    /// of x, y, z and, in a transient case, t (D53). Inlets and walls.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub value: Option<Scalar>,
    /// Diffusive flux into the domain per unit area, Γ ∂φ/∂n (the unit of
    /// the scalar times m/s; a temperature: q_w/(ρ c_p) in K m/s): a number
    /// or an expression (D53). Walls, slip walls, symmetry planes, wedges.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flux: Option<Scalar>,
    /// Periodic pairs (D109 (1)), on one boundary of the pair: `"periodic"`
    /// (the scalar crosses the pair unchanged, the default without a
    /// sub-table) or `"heat-flux"` (it rises by Δ = Q/ṁ from the boundary
    /// the flow enters through to its partner: Q the scalar the walls and
    /// sources put in per second, ρ times their flux and source integrals,
    /// ṁ the mass flow through the pair; the linear trend of a thermally
    /// developed flow with a given wall heat flux).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[schemars(with = "PeriodicMode")]
    pub periodic: Option<PeriodicMode>,
    /// `periodic = "heat-flux"`: the bulk (mixing-cup) value of the scalar
    /// on this boundary, which fixes its level in a steady run (default
    /// 0; a transient run keeps the level of the start value).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bulk: Option<f64>,
}

/// How a scalar crosses a periodic pair (D109 (1), (2)).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum PeriodicMode {
    /// Unchanged (a = 1, b = 0).
    Periodic,
    /// With the jump Δ = Q/ṁ of a uniform heat input (a = 1, b = ±Δ).
    HeatFlux,
}

/// `[solver.scalar]`: the transport equations of the passive scalars.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct ScalarSettings {
    /// Convection scheme of the scalars.
    pub scheme: Scheme,
    /// `scheme = "central"`: fraction of the bounded central scheme.
    pub central_blending: f64,
    /// `scheme = "central"`: β of the normalised-variable bounding.
    pub central_bounding: f64,
    /// Limiter of the scalar gradients.
    pub limiter: LimiterKind,
    /// Constant K of the limiter.
    pub limiter_k: f64,
    /// `limiter = "mlp"`: the limiter function at the vertices (default
    /// Nishikawa's R3).
    pub limiter_function: LimiterFunction,
    /// `limiter = "mlp"`: the threshold ε of the limiter function as a
    /// fraction of the range of the field over the domain, ε = c (φ_max −
    /// φ_min) in the units of the field (Wang 2000): below ε the limiter
    /// lets smooth extrema pass, so a new extremum stays within about ε
    /// (0: strictly bounded).
    pub limiter_epsilon: f64,
    /// `limiter = "mlp"`: ε for the whole domain (`domain`, default) or
    /// per cell, shrinking with the cell size (`cell`: ε_c = c (φ_max −
    /// φ_min) √(V_c / V_domain)). The shrinking ε tightens the bounds as
    /// the mesh is refined and, like ε = 0, can stall a steady residual
    /// on fine meshes (the medium T-junction mesh of gallery 196 with
    /// c = 0.5).
    pub limiter_epsilon_scaling: EpsilonScaling,
    /// `limiter = "mlp"`: relaxation β of the limiter value between outer
    /// iterations, ψ ← ψ + β (ψ_new − ψ). A steady deferred correction can
    /// lock into a two-iteration cycle of ψ in the cells of a mixing layer
    /// (the residual then stalls); β < 1 damps it. The converged ψ is the
    /// unrelaxed one, so the converged solution keeps its bounds (1: no
    /// relaxation).
    pub limiter_relaxation: f64,
    /// Linear-upwind face values clipped to the interval of the two cell
    /// values of the face (the convection boundedness criterion): with a
    /// limiter the converged cell values make no new extrema on any mesh
    /// (default false: the plain linear-upwind face value).
    pub bounded: bool,
    /// Under-relaxation of the scalars (default 0.9 steady, 1 transient).
    pub relaxation: Option<f64>,
    /// Relative tolerance of each scalar solve (default 1e-3 steady, 1e-8
    /// transient).
    pub tolerance: Option<f64>,
    /// Largest number of BiCGStab iterations per scalar solve.
    pub max_iter: usize,
    /// Solves of each scalar per outer iteration (each with new deferred
    /// corrections).
    pub sweeps: usize,
}

impl Default for ScalarSettings {
    fn default() -> Self {
        Self {
            scheme: Scheme::LinearUpwind,
            central_blending: 1.0,
            central_bounding: 0.1,
            limiter: LimiterKind::Venkatakrishnan,
            limiter_k: 5.0,
            limiter_function: LimiterFunction::R3,
            limiter_epsilon: 0.01,
            limiter_epsilon_scaling: EpsilonScaling::Domain,
            limiter_relaxation: 0.5,
            bounded: false,
            relaxation: None,
            tolerance: None,
            max_iter: 50,
            sweeps: 1,
        }
    }
}

impl ScalarSettings {
    /// α of the scalars.
    pub fn alpha(&self, transient: bool) -> f64 {
        self.relaxation.unwrap_or(if transient { 1.0 } else { 0.9 })
    }

    /// Relative tolerance of a solve.
    pub fn tol(&self, transient: bool) -> f64 {
        self.tolerance
            .unwrap_or(if transient { 1e-8 } else { 1e-3 })
    }

    /// Checks the ranges.
    pub fn validate(&self, transient: bool) -> Result<()> {
        let a = self.alpha(transient);
        crate::need!(
            a > 0.0 && a <= 1.0,
            CASE_RANGE,
            "[solver.scalar] relaxation must be in (0, 1], is {a}"
        );
        crate::need!(
            self.tol(transient) > 0.0,
            CASE_RANGE,
            "[solver.scalar] tolerance must be > 0"
        );
        crate::need!(
            self.max_iter >= 1 && self.sweeps >= 1,
            CASE_RANGE,
            "[solver.scalar] max_iter and sweeps must be >= 1"
        );
        crate::need!(
            self.limiter_k > 0.0,
            CASE_RANGE,
            "[solver.scalar] limiter_k must be > 0"
        );
        let top = match self.limiter_epsilon_scaling {
            EpsilonScaling::Domain => 1.0,
            EpsilonScaling::Cell => 100.0,
        };
        crate::need!(
            (0.0..top).contains(&self.limiter_epsilon),
            CASE_RANGE,
            "[solver.scalar] limiter_epsilon must be in [0, {top}) with limiter_epsilon_scaling = {:?}, is {}",
            self.limiter_epsilon_scaling,
            self.limiter_epsilon
        );
        crate::need!(
            self.limiter_relaxation > 0.0 && self.limiter_relaxation <= 1.0,
            CASE_RANGE,
            "[solver.scalar] limiter_relaxation must be in (0, 1], is {}",
            self.limiter_relaxation
        );
        crate::need!(
            (0.0..=1.0).contains(&self.central_blending)
                && (0.0..1.0).contains(&self.central_bounding),
            CASE_RANGE,
            "[solver.scalar] central_blending must be in [0, 1] and central_bounding in [0, 1)"
        );
        Ok(())
    }
}

/// Names a scalar cannot take: the fields of the output and the variables
/// of the expressions.
const RESERVED: &[&str] = &[
    "velocity",
    "velocity_x",
    "velocity_y",
    "velocity_z",
    "speed",
    "pressure",
    "static_pressure",
    "density",
    "velocity_error",
    "wall_distance",
    "x",
    "y",
    "z",
    "t",
];

/// Moves the per-physics sub-tables `scalar` out of the `[boundary.<name>]`
/// tables of a parsed case (the flow conditions keep only their own keys)
/// and returns them per boundary.
pub fn extract_boundary_tables(
    v: &mut toml::Table,
) -> Result<BTreeMap<String, BTreeMap<String, ScalarBoundary>>> {
    let mut out = BTreeMap::new();
    let Some(toml::Value::Table(bs)) = v.get_mut("boundary") else {
        return Ok(out);
    };
    for (name, b) in bs.iter_mut() {
        let toml::Value::Table(t) = b else { continue };
        let Some(s) = t.remove("scalar") else {
            continue;
        };
        let m: BTreeMap<String, ScalarBoundary> = s.try_into().map_err(|e: toml::de::Error| {
            diag::err(
                Diagnostic::new(
                    diag::CASE_BAD_VALUE,
                    format!("[boundary.{name}.scalar]: {}", e.message()),
                )
                .key(format!("boundary.{name}.scalar")),
            )
        })?;
        out.insert(name.clone(), m);
    }
    Ok(out)
}

/// What a boundary of the flow is for a scalar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    /// Walls, slip walls, symmetry planes, wedges: value or flux, default
    /// zero flux.
    WallLike,
    /// Inlets: value, default 0.
    Inlet,
    /// Outlets: zero normal gradient.
    Outlet,
    /// Periodic pairs: across the pair (D109 (1)).
    Periodic,
}

fn side(b: &BoundaryCondition) -> Option<Side> {
    use BoundaryCondition as B;
    Some(match b {
        B::Wall { .. } | B::SlipWall | B::Symmetry | B::Wedge => Side::WallLike,
        B::VelocityInlet { .. }
        | B::MassFlowInlet { .. }
        | B::StaticPressureInlet { .. }
        | B::TotalPressureInlet { .. } => Side::Inlet,
        B::StaticPressureOutlet { .. }
        | B::TotalPressureOutlet { .. }
        | B::ConvectiveOutlet { .. } => Side::Outlet,
        B::Periodic { .. } => Side::Periodic,
    })
}

/// Checks the `[scalars]` tables and the boundary sub-tables of a case.
pub fn validate(
    scalars: &BTreeMap<String, ScalarSection>,
    tables: &BTreeMap<String, BTreeMap<String, ScalarBoundary>>,
    boundary: &BTreeMap<String, BoundaryCondition>,
    settings: &ScalarSettings,
    env: &Env,
    transient: bool,
) -> Result<()> {
    let time = if transient { Vars::T } else { Vars::NONE };
    for (name, s) in scalars {
        let ok = name
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
            && name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_');
        crate::need!(
            ok && !RESERVED.contains(&name.as_str())
                && !expr::reserved(name)
                && !crate::derived::Derived::ALL
                    .iter()
                    .flat_map(|d| d.datasets())
                    .any(|(n, _)| n == name)
                && !name.ends_with("_error"),
            CASE_BAD_VALUE,
            "[scalars.{name}]: a scalar name is letters, digits and _ (not starting with a digit) \
             and not one of the output fields (also the derived ones and names ending in \
             `_error`) or expression names ({})",
            RESERVED.join(", ")
        );
        let d = expr::value(
            &s.diffusivity,
            env,
            Vars::SPACE,
            &format!("scalars.{name}.diffusivity"),
        )?;
        if let Some(v) = d.constant() {
            crate::need!(
                v > 0.0 && v.is_finite(),
                CASE_RANGE,
                "[scalars.{name}] diffusivity must be > 0, is {v}"
            );
        }
        expr::value(
            &s.initial,
            env,
            Vars::SPACE,
            &format!("scalars.{name}.initial"),
        )?;
        if let Some(e) = &s.exact {
            expr::value(
                e,
                env,
                Vars::SPACE.with(time),
                &format!("scalars.{name}.exact"),
            )?;
        }
        expr::value(
            &s.source,
            env,
            Vars::SPACE.with(time),
            &format!("scalars.{name}.source"),
        )?;
        if let Some([lo, hi]) = s.bounds {
            crate::need!(
                lo < hi,
                CASE_RANGE,
                "[scalars.{name}] bounds must be [lower, upper] with lower < upper"
            );
        }
    }
    for (b, m) in tables {
        let Some(bc) = boundary.get(b) else {
            return Err(diag::err(
                Diagnostic::new(
                    diag::SCALAR_BOUNDARY,
                    format!(
                        "[boundary.{b}.scalar] belongs to no boundary condition `[boundary.{b}]`"
                    ),
                )
                .key(format!("boundary.{b}.scalar")),
            ));
        };
        for (n, t) in m {
            let key = format!("boundary.{b}.scalar.{n}");
            if !scalars.contains_key(n) {
                let names: Vec<&str> = scalars.keys().map(String::as_str).collect();
                return Err(diag::err(
                    Diagnostic::new(
                        diag::SCALAR_BOUNDARY,
                        format!(
                            "[{key}]: no scalar `{n}` in [scalars] (the case has {})",
                            if names.is_empty() {
                                "none".into()
                            } else {
                                names.join(", ")
                            }
                        ),
                    )
                    .key(key),
                ));
            }
            let sd = side(bc);
            if sd == Some(Side::Periodic) {
                crate::need!(
                    t.value.is_none() && t.flux.is_none(),
                    SCALAR_BOUNDARY,
                    "[{key}]: a periodic boundary takes `periodic` (\"periodic\" or \"heat-flux\"), \
                     not a value or a flux: the scalar comes across the pair"
                );
                crate::need!(
                    t.periodic.is_some(),
                    SCALAR_BOUNDARY,
                    "[{key}] needs `periodic` (\"periodic\" or \"heat-flux\")"
                );
                crate::need!(
                    t.bulk.is_none() || t.periodic == Some(PeriodicMode::HeatFlux),
                    SCALAR_BOUNDARY,
                    "[{key}]: `bulk` belongs to `periodic = \"heat-flux\"`"
                );
                if let Some(v) = t.bulk {
                    crate::need!(
                        v.is_finite(),
                        CASE_RANGE,
                        "[{key}] bulk must be a finite number"
                    );
                }
                continue;
            }
            crate::need!(
                t.periodic.is_none() && t.bulk.is_none(),
                SCALAR_BOUNDARY,
                "[{key}]: `periodic` and `bulk` belong to periodic boundaries"
            );
            crate::need!(
                t.value.is_some() != t.flux.is_some(),
                SCALAR_BOUNDARY,
                "[{key}] needs one of `value` and `flux`"
            );
            let allowed = match sd {
                Some(Side::WallLike) => true,
                Some(Side::Inlet) => t.value.is_some(),
                _ => false,
            };
            crate::need!(
                allowed,
                SCALAR_BOUNDARY,
                "[{key}]: a {} takes {} (outlets keep a zero normal gradient)",
                bc.type_name(),
                match sd {
                    Some(Side::Inlet) => "a `value`, not a flux",
                    Some(Side::Outlet) => "no scalar condition",
                    _ => "no scalar condition here",
                }
            );
            let vars = Vars::SPACE.with(time);
            if let Some(v) = &t.value {
                expr::value(v, env, vars, &format!("{key}.value"))?;
            }
            if let Some(v) = &t.flux {
                expr::value(v, env, vars, &format!("{key}.flux"))?;
            }
        }
    }
    if !scalars.is_empty() {
        settings.validate(transient)?;
        // the periodic modes (D109 (4)): one sub-table per pair; the
        // heat-flux mode needs a temperature level that no fixed value sets
        let pairs = crate::bc::periodic_pairs(boundary)?;
        for name in scalars.keys() {
            let table = |b: &str| tables.get(b).and_then(|m| m.get(name));
            for p in &pairs {
                let (a, b) = (&p.spec.a, &p.spec.b);
                crate::need!(
                    table(a).is_none() || table(b).is_none(),
                    SCALAR_BOUNDARY,
                    "[boundary.{a}.scalar.{name}] and [boundary.{b}.scalar.{name}]: give the \
                     periodic mode of the pair on one of the two"
                );
                let heat_flux = [a, b]
                    .iter()
                    .any(|x| table(x).is_some_and(|t| t.periodic == Some(PeriodicMode::HeatFlux)));
                if !heat_flux {
                    continue;
                }
                for (bn, bc) in boundary {
                    let t = table(bn);
                    let fixed = match side(bc) {
                        Some(Side::Inlet) => true,
                        Some(Side::WallLike) => t.is_some_and(|t| t.value.is_some()),
                        _ => false,
                    };
                    crate::need!(
                        !fixed,
                        SCALAR_BOUNDARY,
                        "[boundary.{a}.scalar.{name}]: `periodic = \"heat-flux\"` with the fixed \
                         value of `{name}` on `{bn}`: the heat-flux mode is the trend of a flow \
                         heated by given fluxes only; with fixed values the developed field is \
                         plain periodic (use `periodic = \"periodic\"`) or decays exponentially \
                         (the wall-temperature mode of the energy stage)"
                    );
                }
            }
        }
    }
    Ok(())
}

/// How the boundary faces of one scalar are set.
#[derive(Debug, Clone)]
enum FaceRule {
    /// Zero normal gradient (a = 1, b = 0).
    ZeroGradient,
    /// Fixed value.
    Value(Value),
    /// Fixed diffusive flux into the domain; with d_n/Γ_f of the face.
    Flux(Value, f64),
}

/// One scalar on the device.
#[derive(Debug, Clone)]
struct Field {
    name: String,
    phi: Buf,
    grad: Buf3,
    psi: Buf,
    bc_a: Buf,
    bc_b: Buf,
    /// The jump per periodic face of the heat-flux mode (D109 (1)).
    jump: Option<Buf>,
    /// The heat-flux mode, if a pair has it.
    heat_flux: Option<HeatFlux>,
    lsq: Buf3,
    gamma: Buf,
    turb: TurbField,
    su: Option<Buf>,
    eps2: Option<Buf>,
    /// The vertex stencil of `limiter = "mlp"` on the device.
    mlp: Option<usize>,
    /// φⁿ, φⁿ⁻¹ and the φⁿ⁻¹ of the step before (taking a step back).
    old: Option<(Buf, Buf, Buf)>,
    /// The sliding mesh (D182): φⁿ and φⁿ⁻¹ with the data of the moved
    /// cells shifted, which the time term sees instead of `old`.
    content: Option<(Buf, Buf)>,
    /// The last valid value (D63).
    valid: Option<Buf>,
    bounds: (f64, f64),
    rules: Vec<FaceRule>,
    /// The source S, when it is not zero.
    source: Option<Value>,
    time_dependent: bool,
    /// The settings of its solves (`[solver.scalar]`, the energy equation
    /// `[solver.energy]`).
    settings: ScalarSettings,
    /// The energy equation (stage 7a1): the temperature with its property
    /// refresh.
    energy: Option<Box<EnergyField>>,
}

/// A thermal condition per boundary face of the energy equation.
#[derive(Debug, Clone)]
enum ECond {
    /// Zero normal gradient.
    Zero,
    /// Fixed temperature.
    Temperature(Value),
    /// Heat flux into the fluid (W/m²).
    Flux(Value),
    /// h and T_amb.
    Convection(Value, Value),
}

/// The energy equation on the device and its host data (stage 7a1).
#[derive(Debug, Clone)]
struct EnergyField {
    /// The materials: the fluid first, then the solid regions (7d1).
    mats: Vec<crate::energy::Material>,
    /// The material of each cell (index into `mats`; 32 bits, D194: no
    /// fixed limit on the regions).
    mat_cell: Vec<u32>,
    /// The region of each material (`fluid` first).
    regions: Vec<String>,
    /// Refresh c_p, k and μ every outer iteration.
    dependent: bool,
    viscous_heating: bool,
    /// q''' (W/m³), if not zero.
    source: Option<Value>,
    /// ρ c_p per cell.
    capacity: Buf,
    /// Ċ_f = ṁ_f c_p,f per face.
    cflux: Buf,
    /// μ c_p and μ_eff c_p per face (the turbulent conductivity, D135).
    mu_cp: Buf,
    mu_eff_cp: Buf,
    /// Per boundary face.
    conds: Vec<ECond>,
    /// Host copies of the connectivity and geometry.
    owner: Vec<u32>,
    neighbour: Vec<u32>,
    weight: Vec<f64>,
    /// d_n per boundary face.
    dn: Vec<f64>,
    /// Cell volumes.
    volume: Vec<f64>,
    /// c_p per cell and per face, k per face, of the last refresh.
    cp_cell: Vec<f64>,
    cp_face: Vec<f64>,
    k_face: Vec<f64>,
    /// The heat source per cell of the last refresh (W/m³), q''' + Φ.
    heat: Vec<f64>,
    /// Cell centroids (the source expression).
    cells: Vec<[f64; 3]>,
    /// μ(T) per face of the last refresh, if the viscosity follows T.
    mu_face: Option<Vec<f64>>,
    /// h_s(T) per material, in which the equation is conservative (D150).
    h: Vec<crate::energy::Enthalpy>,
    /// The thermal wall function on the fluid side of the interface faces
    /// of a turbulent conjugate case (stage 7d5), with the slot of each
    /// face of the mesh (`u32::MAX`: none).
    iface: Option<(InterfaceWall, Vec<u32>)>,
    /// The thermal law, its Pr_t (of the log layer) and the floor of ỹ⁺ at
    /// the interfaces.
    iface_law: (u32, f64),
    /// Per interface face of the last refresh: the conductance of the fluid
    /// side times d_n (the `k_P/(1 − w)` of the series), for the face
    /// temperature of the report.
    iface_g: Vec<f64>,
    /// The conduction across region interfaces and in orthotropic solids
    /// (stage 7d2): contact resistances, one-sided gradients, cross terms.
    cond: Option<crate::conduction::Conduction>,
    /// The heat sources of single regions (`[[energy.sources]]`): the
    /// cells and the power density (W/m³) as a value of x, y, z, t.
    region_src: Vec<(Vec<u32>, RegionSrc)>,
    /// Unit normals of the faces (the normal conductivity of an
    /// orthotropic solid), when the case has one.
    normals: Vec<[f64; 3]>,
    /// The face area vectors S, when the case has an orthotropic solid.
    areas: Vec<[f64; 3]>,
    /// The factor on ρ c_p per material (thermal acceleration of the
    /// solids, D213; 1 for the fluid and every solid that is not
    /// accelerated).
    inertia: Vec<f64>,
    /// Fixed temperatures (mask, values) of the energy rows: the solids
    /// held while the conductance of the fluid side is measured (D213).
    frozen: Option<(Buf, Buf)>,
    /// Host copy of the mask of `frozen` (the ledger books the heat into
    /// the held cells, D204).
    held: Option<Vec<bool>>,
    /// S_p of the Robin conditions of the solid phase of the cycle loop
    /// (D204), per cell; allocated with the first solid phase.
    robin_sp: Option<Buf>,
    /// Heat added to cells (W per cell): the defect of the mean-flow model
    /// of the cycle loop's solid phase (D204).
    extra_heat: Option<Vec<f64>>,
    /// The variable density of the fluid (stage 7c), set by the flow
    /// after each refresh of ρ(T, p₀).
    var_rho: Option<VarRho>,
}

/// The variable density of the fluid in the energy equation (stage 7c):
/// ρ per cell of the energy mesh (the solids keep their own), ∂ρ/∂t per
/// cell (the mass source of the continuity, for the ledger) and the
/// pressure work dp₀/dt of a closed volume (W/m³ in the fluid cells).
#[derive(Debug, Clone, Default)]
pub struct VarRho {
    /// ρ per cell (kg/m³); NaN where a solid keeps its own density.
    pub rho: Vec<f64>,
    /// ∂ρ/∂t per cell (kg/(m³ s)); empty outside a time step.
    pub drho_dt: Vec<f64>,
    /// dp₀/dt (Pa/s), the heat source of the pressure work in the fluid.
    pub pressure_work: f64,
}

impl EnergyField {
    /// The true density of cell `c` (the variable one of the fluid in stage
    /// 7c, else that of its material; no thermal acceleration).
    fn cell_rho(&self, c: usize) -> f64 {
        match self.var_rho.as_ref().map(|v| v.rho[c]) {
            Some(r) if !r.is_nan() => r,
            _ => self.mats[self.mat_cell[c] as usize].rho,
        }
    }
}

/// The power density of a region source.
#[derive(Debug, Clone)]
enum RegionSrc {
    /// W/m³ as a value of x, y, z (t).
    Density(Value),
    /// W over the region (a value of t) and the region's volume.
    Power(Value, f64),
}

/// The fluid side of the interface faces of a turbulent conjugate case
/// (stage 7d5, stage 7 specification §9.4): the y⁺ of the fluid view's
/// wall treatment per face, with which the energy view puts the wall
/// conductance h_w = ρ c_p u*/T⁺ of the thermal wall function (D135) in
/// series with the solid half cell, q = (T_P − T_N)/(1/h_w + δ_N/k_N).
#[derive(Debug, Clone, Default)]
pub struct InterfaceWall {
    /// The interface faces (interior faces of the full mesh).
    pub faces: Vec<u32>,
    /// True where the fluid cell is the owner of the face.
    pub fluid_owner: Vec<bool>,
    /// The wall distance y_P of the fluid cell (that of the wall
    /// treatment).
    pub y: Vec<f64>,
    /// The distance of the two centroids normal to the face.
    pub dn: Vec<f64>,
    /// y* of the thermal law per face (the fluid view's wall data).
    pub y_star: Vec<f64>,
    /// The laminar viscosity per face (the fluid view's).
    pub mu: Vec<f64>,
    /// The thermal law is taken at ỹ* = max(y*, `y_floor`) (the scalable
    /// treatment).
    pub y_floor: f64,
}

/// The buffers of the flow the energy equation reads (stage 7a1).
#[derive(Debug, Clone, Copy)]
pub struct EnergyLink {
    /// The velocity gradients ∇u_i (viscous heating).
    pub grad_u: [Buf3; 3],
    /// The face viscosity of the momentum equations (μ, or μ + μ_t with a
    /// turbulence model).
    pub mu_face: Buf,
}

/// The heat-flux mode of one scalar on one periodic pair (D109 (1), (2)).
#[derive(Debug, Clone)]
struct HeatFlux {
    /// The pair.
    pair: usize,
    /// Per periodic face: −1 on side `a`, +1 on side `b` of the pair, 0
    /// on the faces of other pairs; jump = Δ · sign.
    sign: Buf,
    /// The input Q/ρ per second of the flux boundaries and the source
    /// (the unit of the scalar times m³/s), at the time of the values.
    q_in: f64,
    /// Δ of the last solve.
    delta: f64,
    /// Bulk value and the side it is given on (−1 `a`, +1 `b`), for the
    /// level of a steady run.
    bulk: (f64, f64),
    /// The faces of side `b`, their owners, the cells across (on side
    /// `a`) and the interpolation weights, for the bulk value.
    faces: Vec<u32>,
    owners: Vec<u32>,
    across: Vec<u32>,
    /// x_f − x_P and x_f − (x_N + T) per face of side `b`, for the face
    /// values from both sides with their gradients.
    to_face: Vec<([f64; 3], [f64; 3])>,
}

/// Device copy of one scalar on the host (restart, output).
#[derive(Debug, Clone, PartialEq)]
pub struct ScalarState {
    /// Name.
    pub name: String,
    /// Cell values.
    pub phi: Vec<f64>,
    /// Transient: the cell values of the level before (φⁿ⁻¹ of the next
    /// step).
    pub old: Option<Vec<f64>>,
    /// The limiter value ψ per cell of `limiter = "mlp"` (D237): its
    /// relaxation (`limiter_relaxation`) carries ψ from one outer iteration
    /// to the next, so a restart needs it to continue bit-identically.
    /// `None` with the other limiters, which keep no memory.
    pub psi: Option<Vec<f64>>,
    /// The boundary coefficients of the last solve (output only; not kept
    /// in restart files).
    pub wall: Option<ScalarWall>,
}

/// The boundary values of a scalar per boundary face: φ_b = a φ_P + b and
/// the face diffusivity Γ_f = ρ (Γ + ν_t/Pr_t) (the wall function included
/// on the walls of a turbulent flow).
#[derive(Debug, Clone, PartialEq)]
pub struct ScalarWall {
    /// a per boundary face (0: fixed value).
    pub a: Vec<f64>,
    /// b per boundary face.
    pub b: Vec<f64>,
    /// ρ Γ per boundary face.
    pub gamma: Vec<f64>,
}

/// Results of the scalar solves of one outer iteration.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ScalarResiduals {
    /// Scaled residual per scalar (before the solve, as the momentum's).
    pub residual: Vec<f64>,
    /// BiCGStab iterations per scalar (all sweeps).
    pub iterations: Vec<usize>,
    /// Cells where the bounded BDF2 fell back to first order, per scalar.
    pub first_order: Vec<usize>,
}

impl ScalarResiduals {
    /// No scalars.
    pub fn is_empty(&self) -> bool {
        self.residual.is_empty()
    }
}

/// Time terms of the scalars in a step.
#[derive(Debug, Clone, Copy)]
pub struct StepTime {
    /// Time of the new level (s).
    pub time: f64,
    /// Δt (s).
    pub dt: f64,
    /// a₀, a₁, a₂.
    pub a: [f64; 3],
}

/// The passive scalars of a flow on its device.
pub struct Scalars {
    fields: Vec<Field>,
    /// Settings.
    pub settings: ScalarSettings,
    nc: usize,
    cp: Buf,
    cn: Buf,
    bf: Buf,
    val: Buf,
    rhs: Buf,
    diag: Buf,
    order: Buf,
    /// Boundary face centroids (for the expressions).
    xf: Vec<[f64; 3]>,
    /// ρ |S| per boundary face and ρ V per cell (the heat input of the
    /// heat-flux mode).
    rho_area: Vec<f64>,
    rho_volume: Vec<f64>,
    /// The last boundary time the time-dependent values were set for.
    t_set: f64,
    /// The time terms of the last step solved (the heat storage of the
    /// balance, stage 7d6).
    last_step: Option<StepTime>,
    /// The vertex stencil on the device, once a field uses
    /// `limiter = "mlp"` (one per mesh, shared by the fields).
    vert: Option<usize>,
    /// The factors of `limiter_epsilon_scaling = "cell"` on the device.
    eps_cell: Option<Buf>,
}

impl Scalars {
    /// The scalars of a case on the device `dev` of a flow with the
    /// discretisation `disc`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        dev: &mut dyn FlowDevice,
        mesh: &Mesh,
        geo: &Geometry,
        disc: &Discretisation,
        rho: f64,
        scalars: &BTreeMap<String, ScalarSection>,
        tables: &BTreeMap<String, BTreeMap<String, ScalarBoundary>>,
        boundary: &BTreeMap<String, BoundaryCondition>,
        settings: &ScalarSettings,
        env: &Env,
        transient: bool,
    ) -> Result<Self> {
        let (nc, nf, ni) = (disc.n_cells, disc.n_faces, disc.n_interior);
        let nb = nf - ni;
        let stored = disc.sell.stored();
        let time = if transient { Vars::T } else { Vars::NONE };
        // boundary of each boundary face
        let mut fb_boundary = vec![0usize; nb];
        for (bi, p) in mesh.boundaries.iter().enumerate() {
            for f in p.range() {
                fb_boundary[f - ni] = bi;
            }
        }
        let xf: Vec<[f64; 3]> = (ni..nf).map(|f| geo.face_centroid[f]).collect();
        // d_n = (x_f − x_P)·n̂ per boundary face
        let dn: Vec<f64> = (ni..nf)
            .map(|f| {
                let (s, _) = geo.normal_area(f);
                let o = mesh.owner[f] as usize;
                let d = [
                    geo.face_centroid[f][0] - geo.cell_centroid[o][0],
                    geo.face_centroid[f][1] - geo.cell_centroid[o][1],
                    geo.face_centroid[f][2] - geo.cell_centroid[o][2],
                ];
                d[0] * s[0] + d[1] * s[1] + d[2] * s[2]
            })
            .collect();
        let mut fields = Vec::new();
        let mut vert = None;
        let eps_cell = if settings.limiter == LimiterKind::Mlp
            && settings.limiter_epsilon_scaling == EpsilonScaling::Cell
            && !scalars.is_empty()
        {
            let b = dev.alloc(nc)?;
            dev.upload(b, &crate::limiter::cell_factors(&disc.cell_volume[..nc]))?;
            Some(b)
        } else {
            None
        };
        for (name, s) in scalars {
            let dif = expr::value(
                &s.diffusivity,
                env,
                Vars::SPACE,
                &format!("scalars.{name}.diffusivity"),
            )?;
            let gam_face: Vec<f64> = (0..nf).map(|f| dif.at(geo.face_centroid[f], 0.0)).collect();
            if let Some(f) = gam_face.iter().position(|g| !(*g > 0.0 && g.is_finite())) {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::EXPR_VALUE,
                        format!(
                            "[scalars.{name}] diffusivity is {} at the face centroid {:?}; it must be > 0",
                            gam_face[f], geo.face_centroid[f]
                        ),
                    )
                    .key(format!("scalars.{name}.diffusivity")),
                ));
            }
            // the face rules and the stencil points of the gradient
            let mut points = BTreeMap::new();
            let mut rules = Vec::with_capacity(nb);
            let mut time_dependent = false;
            let mut per_boundary: Vec<FaceRule> = Vec::new();
            for p in &mesh.boundaries {
                let bc = boundary.get(&p.name);
                let t = tables.get(&p.name).and_then(|m| m.get(name));
                let key = format!("boundary.{}.scalar.{name}", p.name);
                let rule = match (bc.and_then(side), t) {
                    (_, Some(ScalarBoundary { value: Some(v), .. })) => FaceRule::Value(
                        expr::value(v, env, Vars::SPACE.with(time), &format!("{key}.value"))?,
                    ),
                    (_, Some(ScalarBoundary { flux: Some(v), .. })) => FaceRule::Flux(
                        expr::value(v, env, Vars::SPACE.with(time), &format!("{key}.flux"))?,
                        0.0,
                    ),
                    (Some(Side::Inlet), _) => FaceRule::Value(Value::Number(0.0)),
                    _ => FaceRule::ZeroGradient,
                };
                let kind = match rule {
                    FaceRule::Value(_) => BoundaryKind::FixedValue,
                    _ => BoundaryKind::ZeroGradient,
                };
                match &rule {
                    FaceRule::Value(v) | FaceRule::Flux(v, _) => {
                        time_dependent |= v.uses().time();
                    }
                    FaceRule::ZeroGradient => {}
                }
                points.insert(p.name.clone(), kind);
                per_boundary.push(rule);
            }
            for fb in 0..nb {
                let r = per_boundary[fb_boundary[fb]].clone();
                rules.push(match r {
                    FaceRule::Flux(v, _) => {
                        // b = g d_n with g = q/Γ_f
                        FaceRule::Flux(v, dn[fb] / gam_face[ni + fb])
                    }
                    r => r,
                });
            }
            let lsq = crate::simplec::lsq_with(mesh, geo, &points)?;
            let mut put = |v: &[f64]| -> Result<Buf> {
                let b = dev.alloc(v.len())?;
                dev.upload(b, v)?;
                Ok(b)
            };
            let lsq_b = Buf3 {
                x: put(&lsq.x)?,
                y: put(&lsq.y)?,
                z: put(&lsq.z)?,
            };
            let gamma = put(&gam_face.iter().map(|g| rho * g).collect::<Vec<_>>())?;
            let turb = {
                anyhow::ensure!(
                    s.turbulent_prandtl > 0.0 && s.turbulent_prandtl.is_finite(),
                    "[scalars.{name}] turbulent_prandtl is {}; it must be > 0",
                    s.turbulent_prandtl
                );
                use foehncraft_kernels::turb::thermal as th;
                TurbField {
                    gamma0: put(&gam_face.iter().map(|g| rho * g).collect::<Vec<_>>())?,
                    fq: put(&vec![0.0; nb])?,
                    law: match s.wall_function {
                        ThermalWallFunction::Kader => th::KADER,
                        ThermalWallFunction::Jayatilleke => th::JAYATILLEKE,
                        ThermalWallFunction::None => th::NONE,
                    },
                    prt_model: match s.turbulent_prandtl_model {
                        PrandtlModel::Constant => th::PRT_CONSTANT,
                        PrandtlModel::KaysCrawford => th::PRT_KAYS_CRAWFORD,
                    },
                    prt: match s.turbulent_prandtl_model {
                        PrandtlModel::Constant => s.turbulent_prandtl,
                        PrandtlModel::KaysCrawford => th::KC_PRT_INF,
                    },
                    dn: dn.clone(),
                }
            };
            let a: Vec<f64> = rules
                .iter()
                .map(|r| {
                    if matches!(r, FaceRule::Value(_)) {
                        0.0
                    } else {
                        1.0
                    }
                })
                .collect();
            let bc_a = put(&a)?;
            let ini = expr::value(
                &s.initial,
                env,
                Vars::SPACE,
                &format!("scalars.{name}.initial"),
            )?;
            let phi0: Vec<f64> = (0..nc).map(|c| ini.at(geo.cell_centroid[c], 0.0)).collect();
            let phi = put(&phi0)?;
            let bc_b = dev.alloc(nb)?;
            let src = expr::value(
                &s.source,
                env,
                Vars::SPACE.with(time),
                &format!("scalars.{name}.source"),
            )?;
            let source = (src.constant() != Some(0.0)).then_some(src);
            time_dependent |= source.as_ref().is_some_and(|v| v.uses().time());
            let su = match &source {
                Some(_) => Some(dev.alloc(nc)?),
                None => None,
            };
            let eps2 = if settings.limiter == LimiterKind::Venkatakrishnan {
                let e = disc.venkat_eps2(settings.limiter_k);
                let b = dev.alloc(nc)?;
                dev.upload(b, &e)?;
                Some(b)
            } else {
                None
            };
            if settings.limiter == LimiterKind::Mlp && vert.is_none() {
                vert = Some(dev.vertex_setup(&crate::limiter::vertex_stencil(mesh))?);
            }
            let mlp = vert.filter(|_| settings.limiter == LimiterKind::Mlp);
            let psi = dev.alloc(nc)?;
            dev.fill(psi, 1.0)?;
            let old = if transient {
                Some((dev.alloc(nc)?, dev.alloc(nc)?, dev.alloc(nc)?))
            } else {
                None
            };
            let bounds = s
                .bounds
                .map_or((f64::NEG_INFINITY, f64::INFINITY), |[l, h]| (l, h));
            // the heat-flux mode of a periodic pair (D109 (1))
            let mut heat_flux = None;
            let mut jump = None;
            if let Some(per) = &mesh.periodic {
                for (k, pair) in per.pairs.iter().enumerate() {
                    let given =
                        [(&pair.a, -1.0), (&pair.b, 1.0)]
                            .into_iter()
                            .find_map(|(b, sd)| {
                                tables
                                    .get(b)
                                    .and_then(|m| m.get(name))
                                    .filter(|t| t.periodic == Some(PeriodicMode::HeatFlux))
                                    .map(|t| (t.bulk.unwrap_or(0.0), sd))
                            });
                    let Some(bulk) = given else { continue };
                    anyhow::ensure!(
                        heat_flux.is_none(),
                        "[scalars.{name}]: the heat-flux mode on more than one periodic pair"
                    );
                    let np = per.len();
                    let mut sign = vec![0.0; np];
                    for f in pair.a_faces.clone() {
                        sign[f - per.start] = -1.0;
                    }
                    for f in pair.b_faces.clone() {
                        sign[f - per.start] = 1.0;
                    }
                    let sign_b = dev.alloc(np)?;
                    dev.upload(sign_b, &sign)?;
                    jump = Some(dev.alloc(np)?);
                    let faces: Vec<u32> = pair.b_faces.clone().map(|f| f as u32).collect();
                    heat_flux = Some(HeatFlux {
                        pair: k,
                        sign: sign_b,
                        q_in: 0.0,
                        delta: 0.0,
                        bulk,
                        owners: faces.iter().map(|&f| mesh.owner[f as usize]).collect(),
                        across: faces.iter().map(|&f| mesh.neighbour[f as usize]).collect(),
                        to_face: faces
                            .iter()
                            .map(|&f| {
                                let f = f as usize;
                                let xf = geo.face_centroid[f];
                                let xo = geo.cell_centroid[mesh.owner[f] as usize];
                                let xn = geo.cell_centroid[mesh.neighbour[f] as usize];
                                let t = per.shift[f - per.start];
                                (
                                    [xf[0] - xo[0], xf[1] - xo[1], xf[2] - xo[2]],
                                    [
                                        xf[0] - xn[0] - t[0],
                                        xf[1] - xn[1] - t[1],
                                        xf[2] - xn[2] - t[2],
                                    ],
                                )
                            })
                            .collect(),
                        faces,
                    });
                }
            }
            fields.push(Field {
                name: name.clone(),
                phi,
                grad: Buf3 {
                    x: dev.alloc(nc)?,
                    y: dev.alloc(nc)?,
                    z: dev.alloc(nc)?,
                },
                psi,
                bc_a,
                bc_b,
                jump,
                heat_flux,
                lsq: lsq_b,
                gamma,
                turb,
                su,
                eps2,
                mlp,
                old,
                content: None,
                valid: None,
                bounds,
                rules,
                source,
                time_dependent,
                settings: settings.clone(),
                energy: None,
            });
        }
        let mut s = Self {
            fields,
            settings: settings.clone(),
            nc,
            cp: dev.alloc(nf)?,
            cn: dev.alloc(nf)?,
            bf: dev.alloc(nf)?,
            val: dev.alloc(stored)?,
            rhs: dev.alloc(nc)?,
            diag: dev.alloc(nc)?,
            order: dev.alloc(nc)?,
            xf,
            rho_area: (ni..nf).map(|f| rho * geo.normal_area(f).1).collect(),
            rho_volume: geo.cell_volume[..nc].iter().map(|v| rho * v).collect(),
            t_set: f64::NAN,
            last_step: None,
            vert,
            eps_cell,
        };
        let cells: Vec<[f64; 3]> = geo.cell_centroid[..nc].to_vec();
        s.set_values(dev, &cells, 0.0, true)?;
        Ok(s)
    }

    /// Number of scalars.
    pub fn len(&self) -> usize {
        self.fields.len()
    }

    /// No scalars.
    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    /// The names, in the order of the residuals.
    pub fn names(&self) -> Vec<String> {
        self.fields.iter().map(|f| f.name.clone()).collect()
    }

    /// Boundary values and sources at time `t` (only those that depend on
    /// t, unless `all`).
    fn set_values(
        &mut self,
        dev: &mut dyn FlowDevice,
        cells: &[[f64; 3]],
        t: f64,
        all: bool,
    ) -> Result<()> {
        for f in &mut self.fields {
            if !all && !f.time_dependent {
                continue;
            }
            let b: Vec<f64> = f
                .rules
                .iter()
                .zip(&self.xf)
                .map(|(r, &x)| match r {
                    FaceRule::ZeroGradient => 0.0,
                    FaceRule::Value(v) => v.at(x, t),
                    FaceRule::Flux(v, dn_g) => v.at(x, t) * dn_g,
                })
                .collect();
            dev.upload(f.bc_b, &b)?;
            let fq: Vec<f64> = f
                .rules
                .iter()
                .zip(&self.xf)
                .zip(&f.turb.dn)
                .map(|((r, &x), &d)| match r {
                    FaceRule::Flux(v, _) => v.at(x, t) * d,
                    _ => 0.0,
                })
                .collect();
            dev.upload(f.turb.fq, &fq)?;
            let mut src_v = Vec::new();
            if let (Some(src), Some(su)) = (&f.source, f.su) {
                src_v = cells.iter().map(|&x| src.at(x, t)).collect();
                dev.upload(su, &src_v)?;
            }
            if let Some(hf) = f.heat_flux.as_mut() {
                // Q/ρ... times ρ: the flux faces and the source, in face
                // and cell order
                let mut q = 0.0;
                for ((r, &x), &ra) in f.rules.iter().zip(&self.xf).zip(&self.rho_area) {
                    if let FaceRule::Flux(v, _) = r {
                        q += v.at(x, t) * ra;
                    }
                }
                for (s, &rv) in src_v.iter().zip(&self.rho_volume) {
                    q += s * rv;
                }
                hf.q_in = q;
            }
        }
        self.t_set = t;
        Ok(())
    }

    /// The scalar solves of one outer iteration with the fluxes
    /// `mass_flux`; `step` carries the time terms of a time step (and the
    /// new time, for time-dependent values; `cells` are the cell centroids).
    #[allow(clippy::too_many_arguments)]
    pub fn solve(
        &mut self,
        dev: &mut dyn FlowDevice,
        mass_flux: Buf,
        periodic_flow: &[f64],
        step: Option<StepTime>,
        cells: &[[f64; 3]],
        turb: Option<TurbCoupling>,
        link: Option<EnergyLink>,
    ) -> Result<ScalarResiduals> {
        if let Some(st) = step
            && st.time != self.t_set
        {
            self.set_values(dev, cells, st.time, false)?;
        }
        if step.is_some() {
            self.last_step = step;
        }
        let tr = step.is_some();
        let time_now = step.map_or(0.0, |st| st.time);
        let mut out = ScalarResiduals::default();
        // the jumps of the heat-flux mode: Δ = Q/ṁ with the mass flow of
        // the current fluxes, so that the scalar balance of the module
        // closes (D109 (1))
        for f in &mut self.fields {
            if let (Some(hf), Some(jump)) = (f.heat_flux.as_mut(), f.jump) {
                let m = periodic_flow[hf.pair];
                anyhow::ensure!(
                    m.abs() > 0.0,
                    "scalar `{}`: the heat-flux mode needs a mass flow through the periodic pair",
                    f.name
                );
                hf.delta = hf.q_in / m;
                dev.lincomb(jump, hf.delta, hf.sign, 0.0, hf.sign)?;
            }
        }
        let (w_cp, w_cn, w_bf, w_val, w_rhs, w_diag, w_order) = (
            self.cp, self.cn, self.bf, self.val, self.rhs, self.diag, self.order,
        );
        let xf = &self.xf;
        let eps_cell = self.eps_cell;
        for f in self.fields.iter_mut() {
            let s = f.settings.clone();
            let alpha = s.alpha(tr);
            // the energy equation: properties from the current T, the heat
            // capacity flux, the boundary coefficients and the heat source
            let mut conv = mass_flux;
            let mut capacity = None;
            let frozen = f.energy.as_deref().and_then(|e| e.frozen);
            let robin_sp = f
                .energy
                .as_deref()
                .filter(|e| e.cond.as_ref().is_some_and(|c| c.robin.is_some()))
                .and_then(|e| e.robin_sp);
            let mut mus = turb.map(|tc| (tc.mu, tc.mu_eff));
            if let Some(e) = f.energy.as_deref_mut() {
                refresh_energy(
                    dev,
                    e,
                    f.phi,
                    f.bc_a,
                    f.bc_b,
                    f.gamma,
                    f.turb.gamma0,
                    f.turb.fq,
                    f.su,
                    mass_flux,
                    xf,
                    time_now,
                    turb,
                    link,
                    match (step, f.old) {
                        (Some(st), Some((o1, o2, _))) => Some((
                            st,
                            f.content.map_or(o1, |c| c.0),
                            f.content.map_or(o2, |c| c.1),
                        )),
                        _ => None,
                    },
                )?;
                conv = e.cflux;
                capacity = Some(e.capacity);
                if let Some(tc) = turb {
                    // the turbulent conductivity c_p μ_t/Pr_t and the
                    // Prandtl number μ c_p/k of the thermal law
                    let mu = dev.download(tc.mu)?;
                    let me = dev.download(tc.mu_eff)?;
                    let a: Vec<f64> = mu.iter().zip(&e.cp_face).map(|(m, c)| m * c).collect();
                    let b: Vec<f64> = me.iter().zip(&e.cp_face).map(|(m, c)| m * c).collect();
                    dev.upload(e.mu_cp, &a)?;
                    dev.upload(e.mu_eff_cp, &b)?;
                    mus = Some((e.mu_cp, e.mu_eff_cp));
                }
            }
            if let (Some(tc), Some((mu_t, mu_eff_t))) = (turb, mus) {
                // the turbulent diffusivity of this iteration (D135)
                dev.scalar_turb_faces(&ScalarTurbOps {
                    mu: mu_t,
                    mu_eff: mu_eff_t,
                    gamma0: f.turb.gamma0,
                    wall: tc.wall,
                    fq: f.turb.fq,
                    bc_a: f.bc_a,
                    prt: (f.turb.prt_model, f.turb.prt),
                    gamma: f.gamma,
                    bc_b: f.bc_b,
                })?;
                if let (Some(y_plus), true) = (tc.y_plus, f.turb.law != thermal::NONE) {
                    dev.scalar_wall_faces(&ScalarWallOps {
                        y_plus,
                        mu: mu_t,
                        gamma0: f.turb.gamma0,
                        fq: f.turb.fq,
                        bc_a: f.bc_a,
                        law: f.turb.law,
                        prt: f.turb.prt,
                        y_floor: tc.y_floor,
                        gamma: f.gamma,
                        bc_b: f.bc_b,
                        t_plus: None,
                    })?;
                }
            }
            if let Some(e) = f.energy.as_deref() {
                convection_faces(dev, e, f.bc_a, f.bc_b, f.gamma, xf, time_now)?;
            }
            let mut res = f64::NAN;
            let mut its = 0;
            let mut first = 0;
            for sweep in 0..s.sweeps {
                dev.gradient_with(f.lsq, f.phi, f.bc_a, f.bc_b, f.jump, f.grad)?;
                // the one-sided gradients of the interface cells (stage 7d2)
                let (phi_b, grad_b, bca, bcb) = (f.phi, f.grad, f.bc_a, f.bc_b);
                if let Some(cd) = f
                    .energy
                    .as_deref_mut()
                    .and_then(|e| e.cond.as_mut())
                    .filter(|c| c.is_active())
                {
                    cd.one_sided(dev, phi_b, grad_b, bca, bcb)?;
                }
                if let (Some(slot), false) = (f.mlp, s.scheme == Scheme::Upwind) {
                    // ε in the units of the field (Wang 2000)
                    let (mn, mx) = dev.min_max(f.phi)?;
                    dev.mlp_limiter(&MlpOps {
                        slot,
                        phi: f.phi,
                        bc_a: f.bc_a,
                        bc_b: f.bc_b,
                        grad: f.grad,
                        func: s.limiter_function.code(),
                        eps: s.limiter_epsilon * (mx - mn),
                        eps_cell: eps_cell
                            .filter(|_| s.limiter_epsilon_scaling == EpsilonScaling::Cell),
                        relax: s.limiter_relaxation,
                        out: f.psi,
                    })?;
                } else if s.scheme != Scheme::Upwind && s.limiter != LimiterKind::None {
                    let eps2 = match (s.limiter, f.eps2) {
                        (LimiterKind::Venkatakrishnan, Some(b)) => Eps2::Cells(b),
                        // a negative single ε² selects Barth–Jespersen
                        (LimiterKind::BarthJespersen, _) => Eps2::Single(-1.0),
                        _ => {
                            let (mn, mx) = dev.min_max(f.phi)?;
                            let e = s.limiter_k * (mx - mn);
                            Eps2::Single(e * e)
                        }
                    };
                    dev.limiter_with(f.phi, f.bc_a, f.bc_b, f.jump, f.grad, eps2, f.psi)?;
                }
                dev.scalar_faces(&ScalarFaceArgs {
                    bc_a: f.bc_a,
                    bc_b: f.bc_b,
                    jump: f.jump,
                    mass_flux: conv,
                    gamma: f.gamma,
                    phi: f.phi,
                    grad: f.grad,
                    psi: f.psi,
                    blend: if s.scheme == Scheme::Upwind { 0.0 } else { 1.0 },
                    central: if s.scheme == Scheme::Central {
                        s.central_blending
                    } else {
                        0.0
                    },
                    beta: s.central_bounding,
                    nonorth: 1.0,
                    bound: if s.bounded { 1.0 } else { 0.0 },
                    cp: w_cp,
                    cn: w_cn,
                    bf: w_bf,
                })?;
                // the exact interface flux and the orthotropic cross terms
                let gamma_b = f.gamma;
                if let Some(cd) = f
                    .energy
                    .as_deref_mut()
                    .and_then(|e| e.cond.as_mut())
                    .filter(|c| c.is_active())
                {
                    cd.correct(dev, gamma_b, w_bf)?;
                }
                let time = match (step, f.old) {
                    (Some(st), Some((o1, o2, _))) => Some(ScalarTime {
                        inv_dt: 1.0 / st.dt,
                        a: st.a,
                        old: f.content.map_or(o1, |c| c.0),
                        old2: f.content.map_or(o2, |c| c.1),
                        bounds: f.bounds,
                    }),
                    _ => None,
                };
                dev.scalar_rows(&ScalarRowArgs {
                    cp: w_cp,
                    cn: w_cn,
                    bf: w_bf,
                    mass_flux: conv,
                    phi: f.phi,
                    su: f.su,
                    sp: robin_sp,
                    alpha,
                    time,
                    fix: frozen,
                    positive: false,
                    pseudo: None,
                    capacity,
                    val: w_val,
                    rhs: w_rhs,
                    diag: w_diag,
                    order: w_order,
                })?;
                let scale = dev.scalar_scale(f.phi, w_rhs, w_diag, alpha, None)?;
                let (it, r0) = dev.solve(
                    w_val,
                    w_rhs,
                    f.phi,
                    Preconditioner::Jacobi,
                    false,
                    s.tol(tr),
                    s.max_iter,
                )?;
                if sweep == 0 {
                    res = r0 / scale;
                    if time.is_some() {
                        first = dev.sum_sq(w_order)?.round() as usize;
                    }
                }
                its += it;
            }
            // the interface report from the solved temperature
            let phi_b = f.phi;
            if let Some(cd) = f
                .energy
                .as_deref_mut()
                .and_then(|e| e.cond.as_mut())
                .filter(|c| c.is_active())
            {
                cd.settle(dev, phi_b)?;
            }
            out.residual.push(res);
            out.iterations.push(its);
            out.first_order.push(first);
            // the level of a steady heat-flux field: its bulk value where
            // it is given (a transient run keeps the level of its start)
            if let (Some(hf), false) = (&f.heat_flux, tr) {
                let (bulk_b, _) = bulk_of(dev, f.phi, f.grad, mass_flux, hf)?;
                let at = if hf.bulk.1 < 0.0 {
                    bulk_b - hf.delta
                } else {
                    bulk_b
                };
                dev.subtract(f.phi, at - hf.bulk.0)?;
            }
        }
        Ok(out)
    }

    /// Per scalar with the heat-flux mode: the pair, Δ and the bulk values
    /// on side `a` and side `b` (for the report).
    pub fn periodic_bulk(
        &self,
        dev: &mut dyn FlowDevice,
        mass_flux: Buf,
    ) -> Result<Vec<Option<PeriodicBulk>>> {
        self.fields
            .iter()
            .map(|f| match &f.heat_flux {
                Some(hf) => {
                    let (b, _) = bulk_of(dev, f.phi, f.grad, mass_flux, hf)?;
                    Ok(Some((hf.pair, hf.delta, b - hf.delta, b)))
                }
                None => Ok(None),
            })
            .collect()
    }

    /// Start of a time step: φⁿ⁻¹ ← φⁿ (when an old level exists), φⁿ ← φ;
    /// the φⁿ⁻¹ before is kept for [`step_back`](Self::step_back).
    ///
    /// With `extrapolate = Some(ω)` (BDF2 from the second step, ω =
    /// Δt/Δt_prev) the start value of the step is the linear extrapolation
    /// φ ← φⁿ + ω (φⁿ − φⁿ⁻¹): the deferred correction of the convection
    /// is evaluated from φ, so with one outer iteration per step (PISO) it
    /// would otherwise lag by a whole step and make the scheme first order
    /// in time; extrapolated, its lag error is O(Δt²).
    pub fn step_start(
        &mut self,
        dev: &mut dyn FlowDevice,
        has_old: bool,
        extrapolate: Option<f64>,
    ) -> Result<()> {
        for f in &self.fields {
            let Some((o1, o2, keep)) = f.old else {
                continue;
            };
            dev.copy(o2, keep)?;
            if has_old {
                dev.copy(o1, o2)?;
            }
            dev.copy(f.phi, o1)?;
            if has_old && let Some(w) = extrapolate {
                dev.lincomb(f.phi, 1.0 + w, o1, -w, o2)?;
            }
        }
        Ok(())
    }

    /// The sliding mesh (D182): the time terms see the old levels with the
    /// data of the cells that move shifted along (entry i of a shifted
    /// array is entry `perm1[i]` of the Eulerian one for the level n and
    /// `perm2[i]` for the level n−1); call after
    /// [`step_start`](Self::step_start).
    pub fn shift_old(
        &mut self,
        dev: &mut dyn FlowDevice,
        perm1: &[usize],
        perm2: &[usize],
    ) -> Result<()> {
        for f in &mut self.fields {
            let Some((o1, o2, _)) = f.old else {
                continue;
            };
            let (c1, c2) = match f.content {
                Some(c) => c,
                None => {
                    let c = (dev.alloc(self.nc)?, dev.alloc(self.nc)?);
                    f.content = Some(c);
                    c
                }
            };
            for (src, dst, perm) in [(o1, c1, perm1), (o2, c2, perm2)] {
                let v = dev.download(src)?;
                let m: Vec<f64> = perm.iter().map(|&i| v[i]).collect();
                dev.upload(dst, &m)?;
            }
        }
        Ok(())
    }

    /// Takes the last step back exactly: φ ← φⁿ, φⁿ ← φⁿ⁻¹, φⁿ⁻¹ ← the
    /// one before.
    pub fn step_back(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for f in &self.fields {
            let Some((o1, o2, keep)) = f.old else {
                continue;
            };
            dev.copy(o1, f.phi)?;
            dev.copy(o2, o1)?;
            dev.copy(keep, o2)?;
        }
        Ok(())
    }

    /// Keeps a device copy of the values as the last valid state (D63).
    pub fn keep_valid(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for f in &mut self.fields {
            let v = match f.valid {
                Some(v) => v,
                None => dev.alloc(self.nc)?,
            };
            dev.copy(f.phi, v)?;
            f.valid = Some(v);
        }
        Ok(())
    }

    /// Puts the last valid values back.
    pub fn rewind(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for f in &self.fields {
            if let Some(v) = f.valid {
                dev.copy(v, f.phi)?;
            }
        }
        Ok(())
    }

    /// After a rewind of a time step: the old level of n−1 back into φⁿ.
    pub fn rewind_levels(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for f in &self.fields {
            if let Some((o1, o2, _)) = f.old {
                dev.copy(o2, o1)?;
            }
        }
        Ok(())
    }

    /// Host copies (restart, output).
    pub fn states(&self, dev: &mut dyn FlowDevice) -> Result<Vec<ScalarState>> {
        self.fields
            .iter()
            .map(|f| {
                Ok(ScalarState {
                    name: f.name.clone(),
                    phi: dev.download(f.phi)?,
                    old: match f.old {
                        Some((o1, _, _)) => Some(dev.download(o1)?),
                        None => None,
                    },
                    psi: match f.mlp {
                        Some(_) => Some(dev.download(f.psi)?),
                        None => None,
                    },
                    wall: {
                        let g = dev.download(f.gamma)?;
                        let ni = g.len() - f.rules.len();
                        Some(ScalarWall {
                            a: dev.download(f.bc_a)?,
                            b: dev.download(f.bc_b)?,
                            gamma: g[ni..].to_vec(),
                        })
                    },
                })
            })
            .collect()
    }

    /// Loads host copies (a restart): by name; a scalar of the case that
    /// the file does not hold keeps its start value.
    pub fn load(
        &mut self,
        dev: &mut dyn FlowDevice,
        states: &[ScalarState],
    ) -> Result<Vec<String>> {
        let mut missing = Vec::new();
        for f in &self.fields {
            match states.iter().find(|s| s.name == f.name) {
                Some(s) => {
                    dev.upload(f.phi, &s.phi)?;
                    if let (Some((o1, _, _)), Some(o)) = (f.old, &s.old) {
                        dev.upload(o1, o)?;
                    }
                    // the memory of the MLP relaxation (a file without it
                    // starts from ψ = 1, as a new run)
                    if let (Some(_), Some(p)) = (f.mlp, &s.psi)
                        && p.len() == s.phi.len()
                    {
                        dev.upload(f.psi, p)?;
                    }
                }
                None => missing.push(f.name.clone()),
            }
        }
        Ok(missing)
    }

    /// Per scalar: the cell values, the boundary coefficients a and b and
    /// the face diffusivity ρΓ.
    #[allow(clippy::type_complexity)]
    fn host_data(
        &self,
        dev: &mut dyn FlowDevice,
    ) -> Result<Vec<(String, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)>> {
        self.fields
            .iter()
            .map(|f| {
                Ok((
                    f.name.clone(),
                    dev.download(f.phi)?,
                    dev.download(f.bc_a)?,
                    dev.download(f.bc_b)?,
                    dev.download(f.gamma)?,
                ))
            })
            .collect()
    }

    /// What the convection scheme of the energy equation carries out of an
    /// outflow boundary face: the gradient and the limiter of T per cell
    /// and the weights of the linear-upwind and of the central parts, so
    /// that T_out = T_P + w_lu ψ_P ∇T_P·(x_f − x_P) + w_c (T_b − T_P), as in
    /// [`foehncraft_kernels::fv::face_coefficients`].
    fn energy_outflow(&self, dev: &mut dyn FlowDevice) -> Result<Option<Outflow>> {
        let Some(f) = self.fields.iter().find(|f| f.energy.is_some()) else {
            return Ok(None);
        };
        let s = &f.settings;
        let blend = if s.scheme == Scheme::Upwind { 0.0 } else { 1.0 };
        let central = if s.scheme == Scheme::Central {
            s.central_blending
        } else {
            0.0
        };
        Ok(Some(Outflow {
            grad: [
                dev.download(f.grad.x)?,
                dev.download(f.grad.y)?,
                dev.download(f.grad.z)?,
            ],
            psi: dev.download(f.psi)?,
            linear: if central > 0.0 {
                (1.0 - central) * blend
            } else {
                blend
            },
            central,
        }))
    }

    /// The energy equation's h_s(T) and the heat sources over the volume
    /// (W) of the last refresh, if the case has it.
    /// (D194: h_s of the material of each cell, for the boundary enthalpy
    /// flow of every fluid).
    #[allow(clippy::type_complexity)]
    fn energy_host(&self) -> Option<((Vec<crate::energy::Enthalpy>, Vec<u32>), f64)> {
        self.fields.iter().find_map(|f| {
            f.energy.as_deref().map(|e| {
                let src: f64 = e.heat.iter().zip(&e.volume).map(|(h, v)| h * v).sum();
                ((e.h.clone(), e.mat_cell.clone()), src)
            })
        })
    }

    /// The heat stored per region in the last time step (W, transient runs,
    /// stage 7d6): Σ ρ V (a₀ h_s(T) + a₁ h_s(Tⁿ) + a₂ h_s(Tⁿ⁻¹))/Δt over
    /// the cells of each region (`fluid` first), the time terms the energy
    /// equation was solved with; `None` in a steady run.
    pub fn energy_storage(&self, dev: &mut dyn FlowDevice) -> Result<Option<Vec<RegionHeat>>> {
        self.energy_storage_of(dev, None)
    }

    /// The heat stored in the held cells of the cycle loop in the last
    /// time step (W, D204): with the skin swing their temperature moves,
    /// and the heat the holding takes is what flows into them less this.
    fn energy_held_storage(&self, dev: &mut dyn FlowDevice) -> Result<f64> {
        let Some(held) = self
            .fields
            .iter()
            .find_map(|f| f.energy.as_deref())
            .and_then(|e| e.held.clone())
        else {
            return Ok(0.0);
        };
        Ok(self
            .energy_storage_of(dev, Some(&held))?
            .map_or(0.0, |r| r.iter().map(|x| x.storage).sum()))
    }

    /// See [`Scalars::energy_storage`]; with `only`, the cells where it is
    /// true.
    fn energy_storage_of(
        &self,
        dev: &mut dyn FlowDevice,
        only: Option<&[bool]>,
    ) -> Result<Option<Vec<RegionHeat>>> {
        let Some(st) = self.last_step else {
            return Ok(None);
        };
        let Some(f) = self.fields.iter().find(|f| f.energy.is_some()) else {
            return Ok(None);
        };
        let (Some(e), Some((o1, o2, _))) = (f.energy.as_deref(), f.old) else {
            return Ok(None);
        };
        let t = dev.download(f.phi)?;
        let t1 = dev.download(o1)?;
        let t2 = if st.a[2] != 0.0 {
            dev.download(o2)?
        } else {
            vec![0.0; t.len()]
        };
        let [a0, a1, a2] = st.a;
        let mut per = vec![0.0; e.regions.len()];
        for c in 0..e.mat_cell.len() {
            if only.is_some_and(|o| !o[c]) {
                continue;
            }
            let m = e.mat_cell[c] as usize;
            let h = &e.h[m];
            let mut l = a0 * h.at(t[c]) + a1 * h.at(t1[c]);
            if a2 != 0.0 {
                l += a2 * h.at(t2[c]);
            }
            match e.var_rho.as_ref().filter(|v| !v.rho[c].is_nan()) {
                // the bounded convection Σ ṁ (h_f − h_P) adds h_P ∂ρ/∂t
                // to the boundary outflow Σ ṁ h_f: the storage of the
                // ledger is ∂(ρ h)/∂t (stage 7c)
                Some(v) => {
                    let d = v.drho_dt.get(c).copied().unwrap_or(0.0);
                    per[m] += e.volume[c] * (v.rho[c] * l / st.dt + h.at(t[c]) * d);
                }
                None => per[m] += e.mats[m].rho * e.volume[c] * l / st.dt,
            }
        }
        // with a thermal acceleration (D213) the equation stores f times
        // the heat of the true ρ c_p; the rest is booked as inertia scaling
        Ok(Some(
            e.regions
                .iter()
                .zip(per)
                .zip(&e.inertia)
                .map(|((r, s), &f)| RegionHeat {
                    region: r.clone(),
                    storage: f * s,
                    inertia_scaling: (1.0 - f) * s,
                })
                .collect(),
        ))
    }

    /// The heat flow through the faces between two regions of the energy
    /// equation (W, the flux of the last assembly, exact with the one-sided
    /// gradients and the contact resistance of stage 7d2), with the area
    /// and the area-mean face temperature: (name, heat flow, area, mean
    /// temperature). Between the fluid and a solid region the name is
    /// `interface-<solid>`, the heat flows into the fluid and the
    /// temperature is the fluid side's; between two solids
    /// `interface-<a>-<b>`, into `a`.
    pub fn energy_interfaces(&self) -> Vec<(String, f64, f64, f64)> {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref())
            .and_then(|e| e.cond.as_ref())
            .map(|c| c.report())
            .unwrap_or_default()
    }

    /// The heat flow into the domain per boundary face that the cross
    /// terms of orthotropic solids add (W, of the last assembly).
    fn energy_boundary_cross(&self) -> Option<&[f64]> {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref())
            .and_then(|e| e.cond.as_ref())
            .filter(|c| !c.cross.is_empty())
            .map(|c| c.boundary_cross.as_slice())
    }

    /// The regions of the energy equation (`fluid` first, then the solid
    /// regions), the order of [`Scalars::energy_storage`]; empty without
    /// the energy equation.
    pub fn energy_regions(&self) -> Vec<String> {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref().map(|e| e.regions.clone()))
            .unwrap_or_default()
    }

    /// The factor on ρ c_p per region of the energy equation (the order of
    /// [`Scalars::energy_regions`]; thermal acceleration, D213).
    pub fn energy_inertia(&self) -> Vec<f64> {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref().map(|e| e.inertia.clone()))
            .unwrap_or_default()
    }

    /// The temperature per boundary face of the energy equation, T_b = a
    /// T_P + b with the boundary coefficients of the last refresh (stage
    /// 7c: the density of the boundary faces).
    pub fn energy_boundary_temperature(
        &self,
        dev: &mut dyn FlowDevice,
    ) -> Result<Option<Vec<f64>>> {
        let Some(f) = self.fields.iter().find(|f| f.energy.is_some()) else {
            return Ok(None);
        };
        let e = f.energy.as_deref().expect("energy");
        let ni = e.neighbour.len();
        let t = dev.download(f.phi)?;
        let a = dev.download(f.bc_a)?;
        let b = dev.download(f.bc_b)?;
        Ok(Some(
            (0..a.len())
                .map(|fb| a[fb] * t[e.owner[ni + fb] as usize] + b[fb])
                .collect(),
        ))
    }

    /// Sets the variable density of the fluid (stage 7c): it acts from the
    /// next energy solve.
    pub fn set_energy_density(&mut self, v: Option<VarRho>) {
        if let Some(e) = self.fields.iter_mut().find_map(|f| f.energy.as_deref_mut()) {
            e.var_rho = v;
        }
    }

    /// Sets the factor on ρ c_p per region (D213); it acts from the next
    /// energy solve. Only between time steps: the storage of a step is
    /// then that of one factor.
    pub fn set_energy_inertia(&mut self, factors: &[f64]) {
        if let Some(e) = self.fields.iter_mut().find_map(|f| f.energy.as_deref_mut()) {
            for (a, &b) in e.inertia.iter_mut().zip(factors) {
                *a = b;
            }
        }
    }

    /// Per region of the energy equation (the order of
    /// [`Scalars::energy_regions`]): the volume (m³), the true heat
    /// capacity Σ ρ c_p V (J/K, c_p of the last refresh) and the
    /// volume-mean temperature (K).
    pub fn energy_region_state(&self, dev: &mut dyn FlowDevice) -> Result<Vec<(f64, f64, f64)>> {
        let Some(f) = self.fields.iter().find(|f| f.energy.is_some()) else {
            return Ok(Vec::new());
        };
        let Some(e) = f.energy.as_deref() else {
            return Ok(Vec::new());
        };
        let t = dev.download(f.phi)?;
        let mut out = vec![(0.0, 0.0, 0.0); e.regions.len()];
        for (c, (((&m, &v), &cp), &tc)) in e
            .mat_cell
            .iter()
            .zip(&e.volume)
            .zip(&e.cp_cell)
            .zip(&t)
            .enumerate()
        {
            let m = m as usize;
            out[m].0 += v;
            out[m].1 += e.cell_rho(c) * cp * v;
            out[m].2 += tc * v;
        }
        for o in &mut out {
            if o.0 > 0.0 {
                o.2 /= o.0;
            }
        }
        Ok(out)
    }

    /// Per region of the energy equation: the volume means of ρ c_p
    /// (J/(m³ K)) and of the conductivity k (W/(m K); the mean of the
    /// diagonal of an orthotropic region) at the current temperature, from
    /// the materials themselves (D204: the effusivity and the skin depth of
    /// the skin swing, also before the first refresh of the face values).
    pub fn energy_region_material(&self, dev: &mut dyn FlowDevice) -> Result<Vec<(f64, f64)>> {
        let Some(f) = self.fields.iter().find(|f| f.energy.is_some()) else {
            return Ok(Vec::new());
        };
        let Some(e) = f.energy.as_deref() else {
            return Ok(Vec::new());
        };
        let t = dev.download(f.phi)?;
        let mut sum = vec![(0.0, 0.0, 0.0); e.regions.len()];
        for (c, (&m, &v)) in e.mat_cell.iter().zip(&e.volume).enumerate() {
            let mat = &e.mats[m as usize];
            let tc = if e.dependent { t[c] } else { mat.t_ref };
            let k = mat.tensor(tc);
            let s = &mut sum[m as usize];
            s.0 += v;
            s.1 += e.cell_rho(c) * mat.cp.at(tc) * v;
            s.2 += (k[0][0] + k[1][1] + k[2][2]) / 3.0 * v;
        }
        Ok(sum
            .iter()
            .map(|&(v, c, k)| if v > 0.0 { (c / v, k / v) } else { (0.0, 0.0) })
            .collect())
    }

    /// The mean conductivity per region of the energy equation (W/(m K),
    /// over the internal faces between two cells of the region, of the
    /// last refresh; the skin depth of D213).
    pub fn energy_region_conductivity(&self) -> Vec<f64> {
        let Some(e) = self.fields.iter().find_map(|f| f.energy.as_deref()) else {
            return Vec::new();
        };
        let mut sum = vec![(0.0, 0.0); e.regions.len()];
        for (f, &n) in e.neighbour.iter().enumerate() {
            let a = e.mat_cell[e.owner[f] as usize];
            if a == e.mat_cell[n as usize] {
                sum[a as usize].0 += e.k_face[f];
                sum[a as usize].1 += 1.0;
            }
        }
        sum.iter()
            .map(|&(k, n)| if n > 0.0 { k / n } else { 0.0 })
            .collect()
    }

    /// The fluid–solid interface faces of the energy equation with their
    /// state of the last solve (D204): (slot, face, whether the solid cell
    /// is the owner, [area, G, G_fluid, heat into the solid, T_s*, T_w])
    /// as `Conduction::face_states`.
    pub fn energy_fluid_solid_faces(&self) -> Vec<(usize, u32, bool, [f64; 6])> {
        let Some(e) = self.fields.iter().find_map(|f| f.energy.as_deref()) else {
            return Vec::new();
        };
        let Some(cd) = e.cond.as_ref() else {
            return Vec::new();
        };
        let mc = &e.mat_cell;
        let fo: Vec<Option<bool>> = cd
            .faces
            .iter()
            .map(
                |fa| match (mc[fa.o as usize] == 0, mc[fa.n as usize] == 0) {
                    (true, false) => Some(true),
                    (false, true) => Some(false),
                    _ => None,
                },
            )
            .collect();
        cd.face_states(&fo)
            .into_iter()
            .enumerate()
            .zip(&cd.faces)
            .filter_map(|((i, st), fa)| st.map(|x| (i, fa.f, !fo[i].unwrap_or(true), x)))
            .collect()
    }

    /// The report pair (index of [`Scalars::energy_interfaces`]) of the
    /// interface face in `slot`.
    pub fn energy_interface_pair(&self, slot: usize) -> usize {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref())
            .and_then(|e| e.cond.as_ref())
            .map_or(0, |c| c.faces[slot].pair)
    }

    /// The owner and neighbour cells of the interface face in `slot`.
    pub fn energy_interface_cells(&self, slot: usize) -> (u32, u32) {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref())
            .and_then(|e| e.cond.as_ref())
            .map_or((0, 0), |c| (c.faces[slot].o, c.faces[slot].n))
    }

    /// The sum of the heat added to cells from outside (W, D204), if any.
    fn energy_extra_total(&self) -> Option<f64> {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref())
            .and_then(|e| e.extra_heat.as_ref())
            .map(|x| x.iter().sum())
    }

    /// Sets or clears the heat added to cells (W per cell, D204).
    pub(crate) fn set_energy_extra_heat(&mut self, x: Option<Vec<f64>>) {
        if let Some(e) = self.fields.iter_mut().find_map(|f| f.energy.as_deref_mut()) {
            e.extra_heat = x;
        }
    }

    /// The temperature of the energy equation and its old time levels
    /// (φ, φⁿ, φⁿ⁻¹; the levels are φ without a time step).
    pub fn energy_levels(
        &self,
        dev: &mut dyn FlowDevice,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        let Some(f) = self.fields.iter().find(|f| f.energy.is_some()) else {
            return Ok((Vec::new(), Vec::new(), Vec::new()));
        };
        let t = dev.download(f.phi)?;
        Ok(match f.old {
            Some((o1, o2, _)) => (t, dev.download(o1)?, dev.download(o2)?),
            None => (t.clone(), t.clone(), t),
        })
    }

    /// Sets the temperature and its old time levels (D204).
    pub fn set_energy_levels3(
        &mut self,
        dev: &mut dyn FlowDevice,
        t: &[f64],
        o1: &[f64],
        o2: &[f64],
    ) -> Result<()> {
        if let Some(f) = self.fields.iter().find(|f| f.energy.is_some()) {
            dev.upload(f.phi, t)?;
            if let Some((b1, b2, _)) = f.old {
                dev.upload(b1, o1)?;
                dev.upload(b2, o2)?;
            }
        }
        Ok(())
    }

    /// The number of interface faces (slots) of the energy equation.
    pub fn energy_interface_slots(&self) -> usize {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref())
            .and_then(|e| e.cond.as_ref())
            .map_or(0, |c| c.faces.len())
    }

    /// Sets (`Some`, per interface slot) or clears the Robin conditions of
    /// the solid phase of the cycle loop (D204).
    pub(crate) fn set_energy_robin(
        &mut self,
        dev: &mut dyn FlowDevice,
        robin: Option<Vec<Option<crate::conduction::RobinFace>>>,
    ) -> Result<()> {
        let nc = self.nc;
        let Some(e) = self.fields.iter_mut().find_map(|f| f.energy.as_deref_mut()) else {
            return Ok(());
        };
        if robin.is_some() && e.robin_sp.is_none() {
            e.robin_sp = Some(dev.alloc(nc)?);
        }
        if let Some(cd) = e.cond.as_mut() {
            cd.robin = robin;
        }
        Ok(())
    }

    /// Sets the temperature and its old time levels (φ, φⁿ, φⁿ⁻¹) of the
    /// energy equation, so that the next time step stores nothing in the
    /// cells whose value did not change since (D204: the solids updated
    /// between two fluid phases).
    pub fn set_energy_levels(&mut self, dev: &mut dyn FlowDevice, t: &[f64]) -> Result<()> {
        if let Some(f) = self.fields.iter().find(|f| f.energy.is_some()) {
            dev.upload(f.phi, t)?;
            if let Some((o1, o2, _)) = f.old {
                dev.upload(o1, t)?;
                dev.upload(o2, t)?;
            }
        }
        Ok(())
    }

    /// The held cells of the energy equation (D204) and, per held set, the
    /// heat into them through the interface faces from the cells that are
    /// solved (W, the faces without a Robin condition), their heat sources
    /// (W) and the heat into the solids through the Robin conditions (W).
    fn energy_held(&self) -> Option<(&[bool], f64, f64, f64)> {
        let e = self.fields.iter().find_map(|f| f.energy.as_deref())?;
        let held = e.held.as_deref()?;
        let mut iface = 0.0;
        let mut robin = 0.0;
        if let Some(cd) = e.cond.as_ref() {
            for (i, fa) in cd.faces.iter().enumerate() {
                let rb = cd.robin.as_ref().and_then(|r| r[i]);
                // heat[i]: from the owner to the neighbour
                if let Some(r) = rb {
                    robin += if r.solid_owner {
                        -cd.heat[i]
                    } else {
                        cd.heat[i]
                    };
                    continue;
                }
                match (held[fa.o as usize], held[fa.n as usize]) {
                    (false, true) => iface += cd.heat[i],
                    (true, false) => iface -= cd.heat[i],
                    _ => {}
                }
            }
        }
        let src: f64 = held
            .iter()
            .zip(&e.heat)
            .zip(&e.volume)
            .filter(|((h, _), _)| **h)
            .map(|((_, q), v)| q * v)
            .sum();
        Some((held, iface, src, robin))
    }

    /// Per material of the energy equation (as [`Scalars::energy_regions`]):
    /// true for a fluid (D194: one material per fluid of the case).
    pub fn energy_fluid_flags(&self) -> Vec<bool> {
        self.fields
            .iter()
            .find_map(|f| {
                f.energy
                    .as_deref()
                    .map(|e| e.mats.iter().map(|m| m.fluid).collect())
            })
            .unwrap_or_default()
    }

    /// The material (region index) of every cell of the energy equation.
    pub fn energy_cell_regions(&self) -> Vec<u32> {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref().map(|e| e.mat_cell.clone()))
            .unwrap_or_default()
    }

    /// The temperature of the energy equation per cell.
    pub fn energy_temperature(&self, dev: &mut dyn FlowDevice) -> Result<Vec<f64>> {
        match self.fields.iter().find(|f| f.energy.is_some()) {
            Some(f) => dev.download(f.phi),
            None => Ok(Vec::new()),
        }
    }

    /// Sets the temperature of the energy equation per cell (restores the
    /// state after a measurement of the conductances, D213).
    pub fn set_energy_temperature(&mut self, dev: &mut dyn FlowDevice, t: &[f64]) -> Result<()> {
        if let Some(f) = self.fields.iter().find(|f| f.energy.is_some()) {
            dev.upload(f.phi, t)?;
        }
        Ok(())
    }

    /// Holds the temperatures of the cells with `Some(T)` in the following
    /// energy solves (fixed rows), or releases all with `None` (D213).
    pub fn set_energy_frozen(
        &mut self,
        dev: &mut dyn FlowDevice,
        fixed: Option<&[Option<f64>]>,
    ) -> Result<()> {
        let Some(e) = self.fields.iter_mut().find_map(|f| f.energy.as_deref_mut()) else {
            return Ok(());
        };
        e.held = fixed.map(|v| v.iter().map(Option::is_some).collect());
        match fixed {
            None => e.frozen = None,
            Some(v) => {
                let mask: Vec<f64> = v.iter().map(|x| f64::from(u8::from(x.is_some()))).collect();
                let val: Vec<f64> = v.iter().map(|x| x.unwrap_or(0.0)).collect();
                let (bm, bv) = match e.frozen {
                    Some(b) => b,
                    None => (dev.alloc(mask.len())?, dev.alloc(mask.len())?),
                };
                dev.upload(bm, &mask)?;
                dev.upload(bv, &val)?;
                e.frozen = Some((bm, bv));
            }
        }
        Ok(())
    }

    /// The conductance (W/K) from each region of the energy equation to
    /// the outside through its boundaries with a fixed temperature or a
    /// convection condition, at time `t`: Σ A k/d_n, and Σ A h k/d_n/(h +
    /// k/d_n) at convection faces, with the face conductivity of the last
    /// refresh; `areas` are the areas of the boundary faces (D213).
    pub fn energy_outer_conductance(&self, areas: &[f64], t: f64) -> Vec<f64> {
        let Some(e) = self.fields.iter().find_map(|f| f.energy.as_deref()) else {
            return Vec::new();
        };
        let ni = e.neighbour.len();
        let mut g = vec![0.0; e.regions.len()];
        for (fb, c) in e.conds.iter().enumerate() {
            let m = e.mat_cell[e.owner[ni + fb] as usize] as usize;
            let kd = e.k_face[ni + fb] / e.dn[fb];
            let a = areas[fb];
            match c {
                ECond::Temperature(_) => g[m] += a * kd,
                ECond::Convection(h, _) => {
                    let h = h.at(self.xf[fb], t);
                    g[m] += a * h * kd / (h + kd);
                }
                _ => {}
            }
        }
        g
    }

    /// The device buffer of the temperature of the energy equation, if
    /// the scalars have one.
    pub fn energy_phi(&self) -> Option<Buf> {
        self.fields
            .iter()
            .find(|f| f.energy.is_some())
            .map(|f| f.phi)
    }

    /// The laminar viscosity μ(T) per face of the energy equation's last
    /// refresh, when the fluid's viscosity follows the temperature.
    pub fn energy_viscosity(&self) -> Option<&[f64]> {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref().and_then(|e| e.mu_face.as_deref()))
    }

    /// True when the energy equation gives the flow a viscosity that
    /// follows the temperature (μ(T) of the fluid, refreshed with T): the
    /// flow then depends on the temperature (the frozen flow of D118).
    pub fn energy_viscosity_follows(&self) -> bool {
        self.fields
            .iter()
            .find_map(|f| f.energy.as_deref())
            .is_some_and(|e| e.dependent && e.mats.iter().any(|m| m.mu.is_some()))
    }

    /// The fluid side of the interface faces for the next energy solves
    /// (a turbulent conjugate case, stage 7d5): the faces and distances
    /// once, y* and μ of the fluid view every outer iteration.
    pub fn set_energy_interfaces(&mut self, w: InterfaceWall) {
        if let Some(e) = self.fields.iter_mut().find_map(|f| f.energy.as_deref_mut()) {
            let slot = match e.iface.take() {
                Some((_, slot)) => slot,
                None => {
                    let mut slot = vec![u32::MAX; e.owner.len()];
                    for (i, &f) in w.faces.iter().enumerate() {
                        slot[f as usize] = i as u32;
                    }
                    slot
                }
            };
            e.iface = Some((w, slot));
        }
    }

    /// The minimum and maximum of each scalar.
    pub fn extremes(&self, dev: &mut dyn FlowDevice) -> Result<Vec<(f64, f64)>> {
        self.fields.iter().map(|f| dev.min_max(f.phi)).collect()
    }
}

impl Scalars {
    /// Adds the energy equation of `case` (stage 7a1): the field
    /// [`crate::energy::FIELD`] with its thermal conditions, solved with the
    /// settings `[solver.energy]`.
    #[allow(clippy::too_many_arguments)]
    pub fn add_energy(
        &mut self,
        dev: &mut dyn FlowDevice,
        mesh: &Mesh,
        geo: &Geometry,
        disc: &Discretisation,
        case: &crate::case::FlowCase,
        env: &Env,
        transient: bool,
    ) -> Result<()> {
        use crate::energy::Condition;
        let Some(sec) = &case.energy else {
            return Ok(());
        };
        let fluid = crate::energy::fluid(case, env)?;
        let (nc, nf, ni) = (disc.n_cells, disc.n_faces, disc.n_interior);
        let nb = nf - ni;
        let mut area = BTreeMap::new();
        for p in &mesh.boundaries {
            let a: f64 = p.range().map(|f| geo.normal_area(f).1).sum();
            area.insert(p.name.clone(), a);
        }
        let by_name = crate::energy::conditions(case, env, &area)?;
        let mut conds = vec![ECond::Zero; nb];
        let mut points = BTreeMap::new();
        for p in &mesh.boundaries {
            let c = match by_name.get(&p.name) {
                Some(Condition::Temperature(v)) => ECond::Temperature(v.clone()),
                Some(Condition::HeatFlux(v)) => ECond::Flux(v.clone()),
                Some(Condition::Convection(h, t)) => ECond::Convection(h.clone(), t.clone()),
                Some(Condition::ZeroGradient) | None => ECond::Zero,
            };
            let kind = match c {
                ECond::Temperature(_) | ECond::Convection(..) => BoundaryKind::FixedValue,
                _ => BoundaryKind::ZeroGradient,
            };
            points.insert(p.name.clone(), kind);
            for f in p.range() {
                conds[f - ni] = c.clone();
            }
        }
        let dn: Vec<f64> = (ni..nf)
            .map(|f| {
                let (s, _) = geo.normal_area(f);
                let o = mesh.owner[f] as usize;
                (0..3)
                    .map(|d| (geo.face_centroid[f][d] - geo.cell_centroid[o][d]) * s[d])
                    .sum()
            })
            .collect();
        let t_ref = case.fluid.conditions().temperature;
        let k0 = fluid.k.at(t_ref);
        let cp0 = fluid.cp.at(t_ref);
        // the materials: the fluid, the fluid regions with a fluid (or a
        // temperature) of their own (D194), then the solid regions by their
        // zones
        let mut mats = vec![crate::energy::fluid_material(&fluid, "fluid", t_ref)];
        let mut zone_mat = vec![0u32; mesh.zones.len()];
        let mut regions = vec!["fluid".to_string()];
        for (name, r) in &case.regions {
            if r.kind != crate::materials::RegionKind::Fluid
                || !(crate::energy::own_fluid(case, name) || r.temperature.is_some())
            {
                continue;
            }
            let f = crate::energy::fluid_of(case, env, Some(name))?;
            let label = f.material.clone().unwrap_or_else(|| "fluid".to_string());
            mats.push(crate::energy::fluid_material(
                &f,
                &label,
                r.temperature.unwrap_or(t_ref),
            ));
            regions.push(name.clone());
            for z in &r.zones {
                if let Some(i) = mesh.zones.iter().position(|y| y == z) {
                    zone_mat[i] = (mats.len() - 1) as u32;
                }
            }
        }
        for (name, r) in &case.regions {
            if r.kind != crate::materials::RegionKind::Solid {
                continue;
            }
            let m = crate::energy::solid(case, env, name)?;
            mats.push(m);
            regions.push(name.clone());
            for z in &r.zones {
                if let Some(i) = mesh.zones.iter().position(|y| y == z) {
                    zone_mat[i] = (mats.len() - 1) as u32;
                }
            }
        }
        let mat_cell: Vec<u32> = mesh.cell_zone[..nc]
            .iter()
            .map(|&z| zone_mat[z as usize])
            .collect();
        let es = &case.solver.energy;
        let mlp = if es.limiter == LimiterKind::Mlp {
            if self.vert.is_none() {
                self.vert = Some(dev.vertex_setup(&crate::limiter::vertex_stencil(mesh))?);
            }
            if es.limiter_epsilon_scaling == EpsilonScaling::Cell && self.eps_cell.is_none() {
                let b = dev.alloc(nc)?;
                dev.upload(b, &crate::limiter::cell_factors(&disc.cell_volume[..nc]))?;
                self.eps_cell = Some(b);
            }
            self.vert
        } else {
            None
        };
        let mut put = |v: &[f64]| -> Result<Buf> {
            let b = dev.alloc(v.len())?;
            dev.upload(b, v)?;
            Ok(b)
        };
        let lsq = crate::simplec::lsq_with(mesh, geo, &points)?;
        let lsq_b = Buf3 {
            x: put(&lsq.x)?,
            y: put(&lsq.y)?,
            z: put(&lsq.z)?,
        };
        let ini = expr::value(&sec.initial, env, Vars::SPACE, "energy.initial")?;
        let phi0: Vec<f64> = (0..nc).map(|c| ini.at(geo.cell_centroid[c], 0.0)).collect();
        if let Some(c) = phi0.iter().position(|t| !(*t > 0.0 && t.is_finite())) {
            return Err(diag::err(
                Diagnostic::new(
                    diag::EXPR_VALUE,
                    format!(
                        "[energy] initial is {} K at the cell centroid {:?}; a temperature must \
                         be > 0 K",
                        phi0[c], geo.cell_centroid[c]
                    ),
                )
                .key("energy.initial"),
            ));
        }
        let phi = put(&phi0)?;
        // the start coefficients: the fixed temperatures, zero gradients
        // elsewhere (the first refresh sets them all)
        let (mut a0, mut b0) = (vec![1.0; nb], vec![0.0; nb]);
        for (fb, c) in conds.iter().enumerate() {
            if let ECond::Temperature(v) = c {
                a0[fb] = 0.0;
                b0[fb] = v.at(self.xf[fb], 0.0);
            }
        }
        let bc_a = put(&a0)?;
        let bc_b = put(&b0)?;
        let gamma = put(&vec![k0; nf])?;
        let time = if transient { Vars::T } else { Vars::NONE };
        let src = expr::value(&sec.source, env, Vars::SPACE.with(time), "energy.source")?;
        let source = (src.constant() != Some(0.0)).then_some(src);
        let enthalpy: Vec<crate::energy::Enthalpy> = mats
            .iter()
            .map(|m| {
                if sec.temperature_dependent {
                    crate::energy::Enthalpy::of(&m.cp)
                } else {
                    crate::energy::Enthalpy::constant(m.cp.at(m.t_ref))
                }
            })
            .collect();
        let su = if source.is_some()
            || sec.viscous_heating
            || !sec.sources.is_empty()
            || enthalpy.iter().any(|h| !h.is_constant())
            || sec.cycle.is_some()
            || crate::density::wanted(case)
        {
            Some(put(&vec![0.0; nc])?)
        } else {
            None
        };
        let settings = case.solver.energy.clone();
        let eps2 = if settings.limiter == LimiterKind::Venkatakrishnan {
            Some(put(&disc.venkat_eps2(settings.limiter_k))?)
        } else {
            None
        };

        let psi = put(&vec![1.0; nc])?;
        use foehncraft_kernels::turb::thermal as th;
        let turb = TurbField {
            gamma0: put(&vec![k0; nf])?,
            fq: put(&vec![0.0; nb])?,
            law: match sec.wall_function {
                ThermalWallFunction::Kader => th::KADER,
                ThermalWallFunction::Jayatilleke => th::JAYATILLEKE,
                ThermalWallFunction::None => th::NONE,
            },
            prt_model: match sec.turbulent_prandtl_model {
                PrandtlModel::Constant => th::PRT_CONSTANT,
                PrandtlModel::KaysCrawford => th::PRT_KAYS_CRAWFORD,
            },
            prt: match sec.turbulent_prandtl_model {
                PrandtlModel::Constant => sec.turbulent_prandtl,
                PrandtlModel::KaysCrawford => th::KC_PRT_INF,
            },
            dn: dn.clone(),
        };
        let cap0: Vec<f64> = mat_cell
            .iter()
            .map(|&m| {
                let m = &mats[m as usize];
                m.rho * m.cp.at(m.t_ref)
            })
            .collect();
        let capacity = put(&cap0)?;
        let old = if transient {
            Some((dev.alloc(nc)?, dev.alloc(nc)?, dev.alloc(nc)?))
        } else {
            None
        };
        // the region interfaces (stage 7d2): contact resistances, thin
        // walls, one-sided gradients
        let cond = if regions.len() > 1 {
            let specs = crate::cht::interface_specs(case, env, mesh, &mat_cell, &regions)?;
            let fixed: Vec<bool> = conds
                .iter()
                .map(|c| matches!(c, ECond::Temperature(_) | ECond::Convection(..)))
                .collect();
            let (res, of_face) = specs;
            Some(crate::conduction::Conduction::new(
                mesh,
                geo,
                disc,
                &mat_cell,
                &regions,
                &|f| of_face.get(&(f as u32)).copied(),
                res,
                &fixed,
            )?)
        } else {
            None
        };
        // the heat sources of regions
        let mut region_src = Vec::new();
        for (i, src) in sec.sources.iter().enumerate() {
            let key = format!("energy.sources[{i}]");
            let cells = crate::cht::region_cells(case, mesh, &mat_cell, &regions, &src.region)
                .ok_or_else(|| {
                    diag::err(
                        Diagnostic::new(
                            diag::ENE_SOURCE,
                            format!(
                                "[[energy.sources]] region `{}` is not a region of the case \
                                 (`fluid` or a key of [regions]: {})",
                                src.region,
                                case.regions.keys().cloned().collect::<Vec<_>>().join(", ")
                            ),
                        )
                        .key(format!("{key}.region")),
                    )
                })?;
            let v = match (&src.power, &src.power_density) {
                (Some(p), None) => {
                    let vol: f64 = cells.iter().map(|&c| geo.cell_volume[c as usize]).sum();
                    crate::need!(
                        vol > 0.0,
                        ENE_SOURCE,
                        "[[energy.sources]] region `{}` has no cells",
                        src.region
                    );
                    RegionSrc::Power(expr::value(p, env, time, &format!("{key}.power"))?, vol)
                }
                (None, Some(d)) => RegionSrc::Density(expr::value(
                    d,
                    env,
                    Vars::SPACE.with(time),
                    &format!("{key}.power_density"),
                )?),
                _ => {
                    return Err(diag::err(
                        Diagnostic::new(
                            diag::ENE_SOURCE,
                            format!(
                                "[[energy.sources]] region `{}`: give one of `power` (W) and \
                                 `power_density` (W/m³)",
                                src.region
                            ),
                        )
                        .key(key),
                    ));
                }
            };
            region_src.push((cells, v));
        }
        let aniso = mats.iter().any(|m| m.aniso.is_some());
        let (normals, areas) = if aniso {
            (
                (0..nf).map(|f| geo.normal_area(f).0).collect(),
                geo.face_area[..nf].to_vec(),
            )
        } else {
            (Vec::new(), Vec::new())
        };
        let energy = EnergyField {
            dependent: sec.temperature_dependent,
            viscous_heating: sec.viscous_heating,
            source,
            capacity,
            cflux: dev.alloc(nf)?,
            mu_cp: dev.alloc(nf)?,
            mu_eff_cp: dev.alloc(nf)?,
            conds,
            owner: mesh.owner[..nf].to_vec(),
            neighbour: mesh.neighbour[..ni].to_vec(),
            weight: disc.weight[..ni].to_vec(),
            dn,
            volume: geo.cell_volume[..nc].to_vec(),
            cells: geo.cell_centroid[..nc].to_vec(),
            cp_cell: vec![cp0; nc],
            cp_face: vec![cp0; nf],
            k_face: vec![k0; nf],
            heat: vec![0.0; nc],
            mu_face: None,
            h: enthalpy,
            iface: None,
            iface_law: (turb.law, turb.prt),
            iface_g: Vec::new(),
            cond,
            region_src,
            normals,
            areas,
            inertia: vec![1.0; mats.len()],
            frozen: None,
            held: None,
            robin_sp: None,
            extra_heat: None,
            var_rho: None,
            mats,
            mat_cell,
            regions,
        };
        let grad = Buf3 {
            x: dev.alloc(nc)?,
            y: dev.alloc(nc)?,
            z: dev.alloc(nc)?,
        };
        self.fields.push(Field {
            name: crate::energy::FIELD.to_string(),
            phi,
            grad,
            psi,
            bc_a,
            bc_b,
            jump: None,
            heat_flux: None,
            lsq: lsq_b,
            gamma,
            turb,
            su,
            eps2,
            mlp,
            old,
            content: None,
            valid: None,
            bounds: (f64::NEG_INFINITY, f64::INFINITY),
            rules: vec![FaceRule::ZeroGradient; nb],
            source: None,
            time_dependent: false,
            settings,
            energy: Some(Box::new(energy)),
        });
        Ok(())
    }
}

/// The energy refresh before a solve (stage 7a1): c_p, k (and μ) from the
/// current temperature per cell and face, the heat capacity ρ c_p per
/// cell, the heat capacity flux ṁ c_p per face, the face conductivity, the
/// boundary coefficients of the thermal conditions and the heat source
/// (q''' + Φ) per unit heat capacity.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
fn refresh_energy(
    dev: &mut dyn FlowDevice,
    e: &mut EnergyField,
    phi: Buf,
    bc_a: Buf,
    bc_b: Buf,
    gamma: Buf,
    gamma0: Buf,
    fq: Buf,
    su: Option<Buf>,
    mass_flux: Buf,
    xf: &[[f64; 3]],
    t: f64,
    turb: Option<TurbCoupling>,
    link: Option<EnergyLink>,
    step: Option<(StepTime, Buf, Buf)>,
) -> Result<()> {
    let tv = dev.download(phi)?;
    let a_old = dev.download(bc_a)?;
    let b_old = dev.download(bc_b)?;
    let nc = tv.len();
    let nf = e.owner.len();
    let ni = e.neighbour.len();
    let nb = nf - ni;
    let tf: Vec<f64> = (0..nf)
        .map(|f| {
            let o = e.owner[f] as usize;
            if f < ni {
                let w = e.weight[f];
                w * tv[o] + (1.0 - w) * tv[e.neighbour[f] as usize]
            } else {
                a_old[f - ni] * tv[o] + b_old[f - ni]
            }
        })
        .collect();
    let dep = e.dependent;
    let mc = &e.mat_cell;
    let mats = &e.mats;
    // a property at T (or at the material's temperature when they do not
    // follow T)
    let at = |m: &crate::energy::Material, c: &crate::materials::Curve, x: f64| {
        c.at(if dep { x } else { m.t_ref })
    };
    e.cp_cell = (0..nc)
        .map(|c| {
            let m = &mats[mc[c] as usize];
            at(m, &m.cp, tv[c])
        })
        .collect();
    let mut cp_face = Vec::with_capacity(nf);
    let mut k_face = Vec::with_capacity(nf);
    let mut iface_g = vec![0.0; e.iface.as_ref().map_or(0, |(w, _)| w.faces.len())];
    let tt = |m: &crate::energy::Material, x: f64| if dep { x } else { m.t_ref };
    let aniso = !e.normals.is_empty();
    let mut cross = Vec::new();
    // the cross vector K n̂ − k_n n̂ of a cell's material at a face
    let side_cross = |m: &crate::energy::Material, t: f64, n: [f64; 3]| -> [f64; 3] {
        if m.aniso.is_none() {
            return [0.0; 3];
        }
        let k = m.tensor(t);
        let kn = m.k_normal(t, n);
        let mut v = [0.0; 3];
        for i in 0..3 {
            v[i] = (0..3).map(|j| k[i][j] * n[j]).sum::<f64>() - kn * n[i];
        }
        v
    };
    for f in 0..nf {
        let o = e.owner[f] as usize;
        let mo = &mats[mc[o] as usize];
        cp_face.push(at(mo, &mo.cp, tf[f]));
        if f < ni && mc[e.neighbour[f] as usize] != mc[o] {
            // two regions: the series of the two half cells and the
            // contact resistance (D157 §9.1), δ_P/(δ_P+δ_N) = 1 − w
            let n = e.neighbour[f] as usize;
            let w = e.weight[f];
            let cd = e.cond.as_mut().expect("the interfaces of a conjugate case");
            let slot = cd.slot[f] as usize;
            let (nrm, dn, spec) = {
                let fa = &cd.faces[slot];
                (fa.normal, fa.area / fa.dgeo, fa.spec)
            };
            let mn = &mats[mc[n] as usize];
            let kp = mo.k_normal(tt(mo, tv[o]), nrm);
            let kn = mn.k_normal(tt(mn, tv[n]), nrm);
            let rc = spec.map_or(0.0, |i| cd.specs[i].at(tt(mo, 0.5 * (tv[o] + tv[n]))));
            let (mut r_p, mut r_n) = ((1.0 - w) * dn / kp, w * dn / kn);
            match e.iface.as_ref().map(|(iw, slot)| (iw, slot[f])) {
                Some((iw, i)) if i != u32::MAX => {
                    // a turbulent fluid side (§9.4): 1/h_w in place of
                    // δ_P/k_P, h_w = Γ_w/y_P with the wall diffusivity
                    // Γ_w = μ c_p y*/T⁺(ỹ*) of the thermal law (D135); the
                    // solid half cell keeps δ/k
                    let i = i as usize;
                    let cf = if iw.fluid_owner[i] { o } else { n };
                    let fl = &mats[mc[cf] as usize];
                    let cpf = at(fl, &fl.cp, tv[cf]);
                    let kf = at(fl, &fl.k, tv[cf]);
                    let muc = iw.mu[i] * cpf;
                    let yp = iw.y_star[i];
                    let gw = thermal::wall_conductance_at(
                        e.iface_law.0,
                        muc,
                        yp,
                        yp.max(iw.y_floor),
                        muc / kf,
                        e.iface_law.1,
                    );
                    let g = gw / iw.y[i] * iw.dn[i];
                    iface_g[i] = g;
                    if iw.fluid_owner[i] {
                        r_p = dn / g;
                    } else {
                        r_n = dn / g;
                    }
                }
                _ => {}
            }
            k_face.push(dn / (r_p + rc + r_n));
            let fa = &mut cd.faces[slot];
            fa.r_p = r_p;
            fa.r_c = rc;
            fa.r_n = r_n;
            fa.a_p = side_cross(mo, tt(mo, tv[o]), nrm);
            fa.a_n = side_cross(mn, tt(mn, tv[n]), nrm);
        } else if aniso && mo.aniso.is_some() {
            // an orthotropic region: k_n of the face normal in the matrix,
            // the tangential rest K S − k_n S as a deferred correction
            let nrm = e.normals[f];
            let t = tt(mo, tf[f]);
            let kn = mo.k_normal(t, nrm);
            k_face.push(kn);
            let k = mo.tensor(t);
            let sv = e.areas[f];
            let mut tv3 = [0.0; 3];
            for i in 0..3 {
                tv3[i] = (0..3).map(|j| k[i][j] * sv[j]).sum::<f64>() - kn * sv[i];
            }
            cross.push(crate::conduction::CrossFace {
                f: f as u32,
                o: o as u32,
                n: if f < ni { e.neighbour[f] } else { u32::MAX },
                w: if f < ni { e.weight[f] } else { 1.0 },
                t: tv3,
            });
        } else {
            k_face.push(at(mo, &mo.k, tf[f]));
        }
    }
    if let Some(cd) = e.cond.as_mut() {
        cd.cross = cross;
    }
    e.iface_g = iface_g;
    e.cp_face = cp_face;
    e.k_face = k_face;
    for (what, vals, ts, cell) in [
        ("specific_heat", &e.cp_cell, &tv, true),
        ("specific_heat", &e.cp_face, &tf, false),
        ("conductivity", &e.k_face, &tf, false),
    ] {
        if let Some(i) = vals.iter().position(|v| !(*v > 0.0 && v.is_finite())) {
            let m = if cell { mc[i] } else { mc[e.owner[i] as usize] };
            let name = &mats[m as usize].name;
            let (whose, key) = if m == 0 && name == "fluid" {
                ("the fluid's".to_string(), format!("fluid.{what}"))
            } else {
                (
                    format!("material `{name}`'s"),
                    format!("materials.{name}.{what}"),
                )
            };
            return Err(diag::err(
                Diagnostic::new(
                    diag::MAT_VALUE,
                    format!(
                        "{whose} {what} is {} at T = {} K (the temperature of the energy \
                         equation reached it); a property must be positive and finite",
                        vals[i], ts[i]
                    ),
                )
                .key(key),
            ));
        }
    }
    // ρ with the factor of a thermal acceleration (D213): the time term and
    // the storage of h_s take it, the convection (Ċ) does not
    let rho: Vec<f64> = (0..nc)
        .map(|c| e.cell_rho(c) * e.inertia[mc[c] as usize])
        .collect();
    let cap: Vec<f64> = e.cp_cell.iter().zip(&rho).map(|(c, r)| r * c).collect();
    dev.upload(e.capacity, &cap)?;
    let m = dev.download(mass_flux)?;
    let cf: Vec<f64> = m.iter().zip(&e.cp_face).map(|(m, c)| m * c).collect();
    dev.upload(e.cflux, &cf)?;
    // the solid phase of the cycle loop (D204): the fluid–solid faces with
    // a Robin condition carry no conduction in the matrix
    let robin = e
        .cond
        .as_ref()
        .and_then(|cd| cd.robin.as_ref().map(|r| (cd, r)))
        .map(|(cd, r)| {
            let mut kf = e.k_face.clone();
            for (fa, x) in cd.faces.iter().zip(r) {
                if x.is_some() {
                    kf[fa.f as usize] = 0.0;
                }
            }
            (kf, cd.robin_sources())
        });
    let kf = robin.as_ref().map_or(&e.k_face, |r| &r.0);
    dev.upload(gamma, kf)?;
    dev.upload(gamma0, kf)?;
    // μ(T) for the next momentum solve (the laminar viscosity of a
    // turbulence model, which adds μ_t to it)
    // (D194: per face the μ(T) of its fluid; a fluid whose viscosity does
    // not follow T keeps its values)
    if dep && e.mats.iter().any(|m| m.mu.is_some()) {
        let target = turb.map(|tc| tc.mu).or(link.map(|l| l.mu_face));
        if let Some(b) = target {
            let single = e.mats.len() == 1;
            let mut v = if single {
                Vec::new()
            } else {
                match e.mu_face.take() {
                    Some(v) => v,
                    None => dev.download(b)?,
                }
            };
            if single {
                let mu = e.mats[0].mu.as_ref().expect("checked");
                v = tf.iter().map(|&x| mu.at(x)).collect();
            } else {
                let mc = &e.mat_cell;
                for f in 0..nf {
                    let o = mc[e.owner[f] as usize] as usize;
                    let m = if e.mats[o].fluid || f >= ni {
                        o
                    } else {
                        mc[e.neighbour[f] as usize] as usize
                    };
                    if let Some(mu) = &e.mats[m].mu {
                        v[f] = mu.at(tf[f]);
                    }
                }
            }
            dev.upload(b, &v)?;
            e.mu_face = Some(v);
        }
    }
    // the thermal conditions
    let (mut a, mut b, mut q) = (vec![1.0; nb], vec![0.0; nb], vec![0.0; nb]);
    for (fb, c) in e.conds.iter().enumerate() {
        match c {
            ECond::Zero => {}
            ECond::Temperature(v) => {
                a[fb] = 0.0;
                b[fb] = v.at(xf[fb], t);
            }
            ECond::Flux(v) => {
                let qv = v.at(xf[fb], t);
                q[fb] = qv * e.dn[fb];
                b[fb] = q[fb] / e.k_face[ni + fb];
            }
            ECond::Convection(..) => {
                // set from the final face conductivity
                a[fb] = a_old[fb];
                b[fb] = b_old[fb];
            }
        }
    }
    dev.upload(bc_a, &a)?;
    dev.upload(bc_b, &b)?;
    dev.upload(fq, &q)?;
    // the heat source per unit heat capacity
    if let Some(su) = su {
        let mut heat = vec![0.0; nc];
        if let Some(src) = &e.source {
            for (c, h) in heat.iter_mut().enumerate() {
                *h = src.at(e.cells[c], t);
            }
        }
        for (cells, src) in &e.region_src {
            match src {
                RegionSrc::Density(v) => {
                    for &c in cells {
                        let c = c as usize;
                        heat[c] += v.at(e.cells[c], t);
                    }
                }
                RegionSrc::Power(p, vol) => {
                    let q = p.at([0.0; 3], t) / vol;
                    for &c in cells {
                        heat[c as usize] += q;
                    }
                }
            }
        }
        // the pressure work dp₀/dt of a closed volume (stage 7c)
        if let Some(v) = e.var_rho.as_ref().filter(|v| v.pressure_work != 0.0) {
            for (h, r) in heat.iter_mut().zip(&v.rho) {
                if !r.is_nan() {
                    *h += v.pressure_work;
                }
            }
        }
        if let (true, Some(l)) = (e.viscous_heating, link) {
            let mut g = Vec::with_capacity(3);
            for b in &l.grad_u {
                g.push([dev.download(b.x)?, dev.download(b.y)?, dev.download(b.z)?]);
            }
            let mu_f = dev.download(l.mu_face)?;
            let (mut sum, mut cnt) = (vec![0.0; nc], vec![0.0; nc]);
            for f in 0..nf {
                let o = e.owner[f] as usize;
                sum[o] += mu_f[f];
                cnt[o] += 1.0;
                if f < ni {
                    let n = e.neighbour[f] as usize;
                    sum[n] += mu_f[f];
                    cnt[n] += 1.0;
                }
            }
            for c in 0..nc {
                let mu = sum[c] / cnt[c];
                let mut p = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        p += (g[i][j][c] + g[j][i][c]) * g[i][j][c];
                    }
                }
                heat[c] += mu * p;
            }
        }
        // the correction to the h-form (D150): the matrix is in T with
        // the heat capacity flux and ρ c_p; the difference of the upwind
        // convection and of the storage in h_s at the current temperature
        // goes to the right-hand side, so that at convergence the residual
        // is that of h_s (zero for a constant c_p)
        let mut corr = vec![0.0; nc];
        if e.h.iter().any(|h| !h.is_constant()) {
            let hc = |c: usize, x: f64| e.h[e.mat_cell[c] as usize].at(x);
            let hv: Vec<f64> = (0..nc).map(|c| hc(c, tv[c])).collect();
            for f in 0..nf {
                let o = e.owner[f] as usize;
                let mf = m[f];
                let (tu, hu) = if f < ni {
                    let n = e.neighbour[f] as usize;
                    if mf >= 0.0 {
                        (tv[o], hv[o])
                    } else {
                        (tv[n], hv[n])
                    }
                } else if mf >= 0.0 {
                    (tv[o], hv[o])
                } else {
                    (tf[f], hc(o, tf[f]))
                };
                corr[o] += cf[f] * (tu - tv[o]) - mf * (hu - hv[o]);
                if f < ni {
                    let n = e.neighbour[f] as usize;
                    corr[n] += -cf[f] * (tu - tv[n]) + mf * (hu - hv[n]);
                }
            }
            if let Some((st, o1, o2)) = step {
                let t1 = dev.download(o1)?;
                let t2 = if st.a[2] != 0.0 {
                    dev.download(o2)?
                } else {
                    vec![0.0; nc]
                };
                let [a0, a1, a2] = st.a;
                for c in 0..nc {
                    let lt = cap[c] * (a0 * tv[c] + a1 * t1[c] + a2 * t2[c]);
                    let lh = rho[c]
                        * (a0 * hv[c]
                            + a1 * hc(c, t1[c])
                            + if a2 != 0.0 { a2 * hc(c, t2[c]) } else { 0.0 });
                    corr[c] += e.volume[c] * (lt - lh) / st.dt;
                }
            }
        }
        let mut v: Vec<f64> = (0..nc)
            .map(|c| (heat[c] + corr[c] / e.volume[c]) / cap[c])
            .collect();
        // the Robin conditions of the cycle loop's solid phase (D204):
        // Ḡ (T̄_ref − δT) explicit, −Ḡ on the diagonal, per unit heat capacity
        if let Some(x) = &e.extra_heat {
            for (c, q) in x.iter().enumerate() {
                v[c] += q / (e.volume[c] * cap[c]);
            }
        }
        if let (Some((_, rs)), Some(spb)) = (&robin, e.robin_sp) {
            let mut sp = vec![0.0; nc];
            for &(c, q, g) in rs {
                let c = c as usize;
                let w = e.volume[c] * cap[c];
                v[c] += q / w;
                sp[c] -= g / w;
            }
            dev.upload(spb, &sp)?;
        }
        dev.upload(su, &v)?;
        e.heat = heat;
    }
    Ok(())
}

/// The convection conditions h (T_amb - T_b) = k (T_b - T_P)/d_n with the
/// final face conductivity of the iteration (after the turbulent one and
/// the wall function): T_b = a T_P + b.
fn convection_faces(
    dev: &mut dyn FlowDevice,
    e: &EnergyField,
    bc_a: Buf,
    bc_b: Buf,
    gamma: Buf,
    xf: &[[f64; 3]],
    t: f64,
) -> Result<()> {
    if !e.conds.iter().any(|c| matches!(c, ECond::Convection(..))) {
        return Ok(());
    }
    let g = dev.download(gamma)?;
    let mut a = dev.download(bc_a)?;
    let mut b = dev.download(bc_b)?;
    let ni = e.neighbour.len();
    for (fb, c) in e.conds.iter().enumerate() {
        if let ECond::Convection(h, ta) = c {
            let kd = g[ni + fb] / e.dn[fb];
            let h = h.at(xf[fb], t);
            a[fb] = kd / (kd + h);
            b[fb] = h * ta.at(xf[fb], t) / (kd + h);
        }
    }
    dev.upload(bc_a, &a)?;
    dev.upload(bc_b, &b)?;
    Ok(())
}

/// A scalar on a heat-flux pair: the pair, Δ and the bulk values on side
/// `a` and side `b`.
pub type PeriodicBulk = (usize, f64, f64, f64);

/// The bulk (mixing-cup) value Σ ṁ_f φ_f / Σ ṁ_f over the faces of side
/// `b` of a heat-flux pair and Σ ṁ_f: φ_f the mean of the values
/// reconstructed with the gradients from the owner and from the cell
/// across (seen with the jump +Δ), second order also on skewed faces.
fn bulk_of(
    dev: &mut dyn FlowDevice,
    phi: Buf,
    grad: Buf3,
    mass_flux: Buf,
    hf: &HeatFlux,
) -> Result<(f64, f64)> {
    let po = dev.gather(phi, &hf.owners)?;
    let pn = dev.gather(phi, &hf.across)?;
    let m = dev.gather(mass_flux, &hf.faces)?;
    let mut go = Vec::new();
    let mut gn = Vec::new();
    for g in grad.all() {
        go.push(dev.gather(g, &hf.owners)?);
        gn.push(dev.gather(g, &hf.across)?);
    }
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..m.len() {
        let (dxo, dxn) = hf.to_face[i];
        let fo = po[i] + (0..3).map(|d| go[d][i] * dxo[d]).sum::<f64>();
        let fnb = pn[i] + hf.delta + (0..3).map(|d| gn[d][i] * dxn[d]).sum::<f64>();
        let pf = 0.5 * (fo + fnb);
        num += m[i] * pf;
        den += m[i];
    }
    Ok((num / den, den))
}

/// BDF coefficients of the scalars for a step: BDF2 with ω = Δt/Δt_prev
/// when `second`, else BDF1.
pub fn bdf(second: bool, omega: f64) -> [f64; 3] {
    if second {
        [
            (1.0 + 2.0 * omega) / (1.0 + omega),
            -(1.0 + omega),
            omega * omega / (1.0 + omega),
        ]
    } else {
        [1.0, -1.0, 0.0]
    }
}

/// Adds the passive scalars of `case` to the flow solver `fl` (no-op
/// without `[scalars]`); call before loading a restart.
pub fn attach(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
) -> Result<()> {
    if case.scalars.is_empty() && case.energy.is_none() {
        return Ok(());
    }
    // conjugate heat transfer: the energy equation on the full mesh
    if let Some(v) = &case.fluid_view {
        return crate::cht::attach(fl, v.clone());
    }
    let env = case.env()?;
    let transient = case.time.is_some();
    let mut s = Scalars::new(
        fl.dev.as_mut(),
        mesh,
        geo,
        &fl.disc,
        fl.phys.rho,
        &case.scalars,
        &case.boundary_scalar,
        &case.boundary,
        &case.solver.scalar,
        &env,
        transient,
    )?;
    s.add_energy(fl.dev.as_mut(), mesh, geo, &fl.disc, case, &env, transient)?;
    fl.add_scalars(s);
    Ok(())
}

/// The outflow face values of the energy equation's convection scheme
/// ([`Scalars::energy_outflow`]).
struct Outflow {
    grad: [Vec<f64>; 3],
    psi: Vec<f64>,
    linear: f64,
    central: f64,
}

impl Outflow {
    /// The value carried out through boundary face `f` of cell `o` with
    /// the boundary value `tb`.
    fn at(&self, geo: &Geometry, f: usize, o: usize, tp: f64, tb: f64) -> f64 {
        let r: f64 = (0..3)
            .map(|d| self.grad[d][o] * (geo.face_centroid[f][d] - geo.cell_centroid[o][d]))
            .sum();
        tp + self.linear * self.psi[o] * r + self.central * (tb - tp)
    }
}

/// A scalar on one boundary, for the report.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct ScalarBoundarySummary {
    /// Boundary name.
    pub boundary: String,
    /// Area-weighted mean of the face values.
    pub mean: f64,
    /// Mean weighted with the mass flux through the faces (the mixed-mean
    /// or bulk value), where mass crosses the boundary.
    pub mixed_mean: Option<f64>,
    /// Area-weighted mean diffusive flux into the domain, Γ ∂φ/∂n (the unit
    /// of the scalar times m/s); of the temperature: the heat flux k ∂T/∂n
    /// into the domain (W/m²).
    pub flux: f64,
    /// The temperature: the heat flow by conduction into the domain (W).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub heat_flow: Option<f64>,
    /// The temperature: the sensible enthalpy the flow carries into the
    /// domain, −Σ ṁ_f h_s(T_f) (W), where mass crosses the boundary.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub enthalpy_flow: Option<f64>,
}

/// The heat balance of the energy equation over the domain (stage 7a1).
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct EnergyBalance {
    /// Conduction into the domain through the boundaries (W).
    pub conduction: f64,
    /// Enthalpy carried into the domain by the flow (W).
    pub convection: f64,
    /// The heat sources q''' + Φ over the volume (W).
    pub source: f64,
    /// The heat stored in the last time step (W), transient runs (stage
    /// 7d6).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub storage: Option<f64>,
    /// The stored heat per region (`fluid` first, then the solid regions).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub storage_regions: Vec<RegionHeat>,
    /// The heat the true ρ c_p of accelerated solids would have stored in
    /// addition (W, D213), when a region is accelerated: conduction +
    /// convection + source − storage − inertia_scaling is the balance of
    /// the true solids, which the accelerated run does not close.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inertia_scaling: Option<f64>,
    /// The heat into the cells held at fixed temperatures (W, D204): the
    /// solids in a fluid phase of the cycle-averaged loop, the fluid in its
    /// solid phase (their heat flows in through faces and boundaries and
    /// their sources, and is not stored).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub held: Option<f64>,
    /// The heat into the solids from outside the model in the solid phase
    /// of the cycle-averaged loop (W, D204): through the Robin conditions,
    /// or the defect of the mean-flow model.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub robin: Option<f64>,
    /// conduction + convection + source − storage (W): zero in a
    /// converged steady run and in every converged time step.
    pub imbalance: f64,
    /// |imbalance| over the largest single heat flow of a boundary, the
    /// source or the storage.
    pub relative: f64,
}

/// The heat stored in a region in the last time step (stage 7d6).
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct RegionHeat {
    /// The region (`fluid` or a solid region).
    pub region: String,
    /// The heat stored (W), with the heat capacity the equation was solved
    /// with.
    pub storage: f64,
    /// The heat the true ρ c_p would have stored in addition (W): (1 − f)
    /// times the true storage, zero without a thermal acceleration (D213).
    /// storage + inertia_scaling is the true storage.
    #[serde(skip_serializing_if = "is_zero")]
    pub inertia_scaling: f64,
}

fn is_zero(x: &f64) -> bool {
    *x == 0.0
}

/// A scalar in the report of a solve.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct ScalarSummary {
    /// Name.
    pub name: String,
    /// Smallest cell value.
    pub min: f64,
    /// Largest cell value.
    pub max: f64,
    /// `bounds` of the scalar, if set.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bounds: Option<[f64; 2]>,
    /// Mesh size h (as the velocity errors).
    pub h: f64,
    /// Volume-weighted L2 error against `exact` (absolute).
    pub l2: Option<f64>,
    /// Largest error at a cell centroid against `exact` (absolute).
    pub linf: Option<f64>,
    /// Per boundary.
    pub boundaries: Vec<ScalarBoundarySummary>,
    /// The temperature of the energy equation: the heat balance.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub energy: Option<EnergyBalance>,
}

/// The scalars of `fl` for the report, at time `t` (for the exact
/// solutions); `dim` is 2 for a mesh one cell layer thick.
pub fn summaries(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    dim: u8,
    t: f64,
) -> Result<Vec<ScalarSummary>> {
    if fl.cht.is_some() {
        return crate::cht::summaries(fl, dim, t);
    }
    let Some(sc) = fl.scalars.as_ref() else {
        return Ok(Vec::new());
    };
    let data = sc.host_data(fl.dev.as_mut())?;
    // the heat-flux pairs: Δ and the bulk values of the level fix
    let deltas = sc.periodic_bulk(fl.dev.as_mut(), fl.f.mass_flux)?;
    let m = fl.dev.download(fl.f.mass_flux)?;
    let env = case.env()?;
    let transient = case.time.is_some();
    let time = if transient { Vars::T } else { Vars::NONE };
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let rho = fl.phys.rho;
    let vol: f64 = geo.cell_volume[..nc].iter().sum();
    let h = if dim == 2 {
        let (lo, hi) = mesh.bounding_box();
        (vol / (hi[2] - lo[2]) / nc as f64).sqrt()
    } else {
        (vol / nc as f64).cbrt()
    };
    let mut out = Vec::new();
    for ((name, phi, a, b, gamma), delta) in data.into_iter().zip(deltas) {
        let (mut mn, mut mx) = (f64::INFINITY, f64::NEG_INFINITY);
        for &v in &phi {
            mn = mn.min(v);
            mx = mx.max(v);
        }
        let (mut l2, mut linf) = (None, None);
        if let Some((e, key)) = exact_of(case, &name) {
            let ex = expr::value(e, &env, Vars::SPACE.with(time), &key)?;
            let (mut e2, mut ei) = (0.0, 0.0f64);
            for ((&v, &x), &vc) in phi[..nc]
                .iter()
                .zip(&geo.cell_centroid[..nc])
                .zip(&geo.cell_volume[..nc])
            {
                let d = v - ex.at(x, t);
                e2 += d * d * vc;
                ei = ei.max(d.abs());
            }
            l2 = Some((e2 / vol).sqrt());
            linf = Some(ei);
        }
        let mut bs = Vec::new();
        // the energy equation: the flux in W/m² (k, not ρΓ), the heat
        // flows and the balance
        let heat = (name == crate::energy::FIELD)
            .then(|| sc.energy_host())
            .flatten();
        let div = if heat.is_some() { 1.0 } else { rho };
        let outflow = if heat.is_some() {
            sc.energy_outflow(fl.dev.as_mut())?
        } else {
            None
        };
        let cross = if heat.is_some() {
            sc.energy_boundary_cross()
        } else {
            None
        };
        let mut bal = EnergyBalance::default();
        let mut largest = 0.0f64;
        let held = if heat.is_some() {
            sc.energy_held()
        } else {
            None
        };
        let mut held_b = 0.0;
        for p in &mesh.boundaries {
            let (mut area, mut sum, mut flux, mut mf, mut msum) = (0.0, 0.0, 0.0, 0.0, 0.0);
            let mut enthalpy = 0.0;
            for f in p.range() {
                let fb = f - ni;
                let o = mesh.owner[f] as usize;
                let (_, ar) = geo.normal_area(f);
                let pf = a[fb] * phi[o] + b[fb];
                let (flux0, enthalpy0) = (flux, enthalpy);
                area += ar;
                sum += ar * pf;
                // the diffusive flux into the domain, orthogonal part (the
                // whole flux at a flux face): ρΓ |S|²/(d·S) (φ_f − φ_P)
                let dgeo = fl.disc.dgeo[f];
                flux += gamma[f] / div * dgeo * (pf - phi[o]);
                if let (Some(w), true) = (&outflow, heat.is_some()) {
                    // the energy equation: the whole flux of the assembly,
                    // with its non-orthogonal part (1 − a) k ∇T_P·k_f and the
                    // cross term of an orthotropic solid (stage 7d2), so that
                    // the balance closes on skewed meshes too
                    let kv = [fl.disc.k.x[f], fl.disc.k.y[f], fl.disc.k.z[f]];
                    let g: f64 = (0..3).map(|d| w.grad[d][o] * kv[d]).sum();
                    flux += (1.0 - a[fb]) * gamma[f] * g;
                    if let Some(x) = cross {
                        flux += x[fb];
                    }
                }
                mf += m[f].abs();
                msum += m[f].abs() * pf;
                if let Some(((hs, mc), _)) = &heat {
                    let h = &hs[mc[o] as usize];
                    // the value the scheme convects: the boundary value at
                    // inflow, the scheme's reconstruction at outflow
                    let tc = match &outflow {
                        Some(w) if m[f] > 0.0 => w.at(geo, f, o, phi[o], pf),
                        _ => pf,
                    };
                    enthalpy -= m[f] * h.at(tc);
                }
                if held.is_some_and(|h| h.0[o]) {
                    held_b += flux - flux0 + enthalpy - enthalpy0;
                }
            }
            if area > 0.0 {
                let (heat_flow, enthalpy_flow) = match heat {
                    Some(_) => {
                        bal.conduction += flux;
                        bal.convection += enthalpy;
                        largest = largest.max(flux.abs()).max(enthalpy.abs());
                        (Some(flux), (mf > 0.0).then_some(enthalpy))
                    }
                    None => (None, None),
                };
                bs.push(ScalarBoundarySummary {
                    boundary: p.name.clone(),
                    mean: sum / area,
                    mixed_mean: (mf > 0.0).then(|| msum / mf),
                    flux: flux / area,
                    heat_flow,
                    enthalpy_flow,
                });
            }
        }
        let stored = match heat {
            Some(_) => sc.energy_storage(fl.dev.as_mut())?,
            None => None,
        };
        let held_stored = match held {
            Some(_) => sc.energy_held_storage(fl.dev.as_mut())?,
            None => 0.0,
        };
        let energy = heat.as_ref().map(|(_, src)| {
            bal.source = *src;
            largest = largest.max(src.abs());
            let total: f64 = stored.iter().flatten().map(|r| r.storage).sum();
            if let Some(rs) = &stored {
                bal.storage = Some(total);
                largest = largest.max(total.abs());
                for r in rs {
                    largest = largest.max(r.storage.abs());
                }
                bal.storage_regions = rs.clone();
                if rs.iter().any(|r| r.inertia_scaling != 0.0) {
                    bal.inertia_scaling = Some(rs.iter().map(|r| r.inertia_scaling).sum());
                }
            }
            // the held cells of the cycle loop (D204): what flows into them
            // is not stored; the Robin heat comes from outside the mesh
            let mut out = 0.0;
            // the defect of the cycle loop's mean-flow model (D204), heat
            // from outside the model like the Robin heat
            if let Some(x) = sc.energy_extra_total() {
                bal.robin = Some(x);
                largest = largest.max(x.abs());
                out -= x;
            }
            if let Some((_, iface, src_held, robin)) = held {
                // less what the held cells store when their temperature
                // moves (the skin swing)
                let h = held_b + iface + src_held - held_stored;
                bal.held = Some(h);
                largest = largest.max(h.abs());
                out += h;
                if robin != 0.0 {
                    bal.robin = Some(bal.robin.unwrap_or(0.0) + robin);
                    largest = largest.max(robin.abs());
                    out -= robin;
                }
            }
            bal.imbalance = bal.conduction + bal.convection + bal.source - total - out;
            bal.relative = if largest > 0.0 {
                bal.imbalance.abs() / largest
            } else {
                0.0
            };
            bal.clone()
        });
        // the two sides of each periodic pair: the face value interpolated
        // with the cell across (with the jump of the heat-flux mode), the
        // diffusive flux into the domain the axial conduction
        if let Some(per) = &mesh.periodic {
            for (k, pair) in per.pairs.iter().enumerate() {
                let hf = delta.filter(|d| d.0 == k);
                let dk = hf.map_or(0.0, |d| d.1);
                for (bname, faces, sd, bulk) in [
                    (&pair.a, pair.a_faces.clone(), -1.0, hf.map(|d| d.2)),
                    (&pair.b, pair.b_faces.clone(), 1.0, hf.map(|d| d.3)),
                ] {
                    let (mut area, mut sum, mut flux, mut mf, mut msum) = (0.0, 0.0, 0.0, 0.0, 0.0);
                    for f in faces {
                        let o = mesh.owner[f] as usize;
                        let n = mesh.neighbour[f] as usize;
                        let (_, ar) = geo.normal_area(f);
                        let pn = phi[n] + sd * dk;
                        let w = fl.disc.weight[f];
                        let pf = w * phi[o] + (1.0 - w) * pn;
                        area += ar;
                        sum += ar * pf;
                        flux += gamma[f] / div * fl.disc.dgeo[f] * (pn - phi[o]);
                        mf += m[f].abs();
                        msum += m[f].abs() * pf;
                    }
                    if area > 0.0 {
                        bs.push(ScalarBoundarySummary {
                            boundary: bname.clone(),
                            mean: sum / area,
                            mixed_mean: bulk.or((mf > 0.0).then(|| msum / mf)),
                            flux: flux / area,
                            heat_flow: None,
                            enthalpy_flow: None,
                        });
                    }
                }
            }
        }
        out.push(ScalarSummary {
            bounds: case.scalars.get(&name).and_then(|s| s.bounds),
            name,
            min: mn,
            max: mx,
            h,
            l2,
            linf,
            boundaries: bs,
            energy,
        });
    }
    Ok(out)
}

/// The error fields `<name>_error` = φ − φ_exact of the scalars with an
/// exact solution, for the field output; `at` as in the fields file (the
/// time of a transient step).
pub fn field_errors(
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &crate::simplec::Solution,
    at: Option<(usize, Option<f64>)>,
) -> Result<Vec<(String, Vec<f64>)>> {
    let mut out = Vec::new();
    if fl.state.scalars.is_empty() {
        return Ok(out);
    }

    let env = case.env()?;
    let time = if case.time.is_some() {
        Vars::T
    } else {
        Vars::NONE
    };
    let t = at.and_then(|(_, t)| t).unwrap_or(0.0);
    for st in &fl.state.scalars {
        if let Some((e, key)) = exact_of(case, &st.name) {
            let ex = expr::value(e, &env, Vars::SPACE.with(time), &key)?;
            let err = (0..mesh.n_cells)
                .map(|c| st.phi[c] - ex.at(geo.cell_centroid[c], t))
                .collect();
            out.push((format!("{}_error", st.name), err));
        }
    }
    Ok(out)
}

/// One log line of a scalar summary.
pub fn summary_line(s: &ScalarSummary) -> String {
    let mut l = format!("scalar {}: min {:.6e}, max {:.6e}", s.name, s.min, s.max);
    if let Some([lo, hi]) = s.bounds {
        let inside = s.min >= lo && s.max <= hi;
        l.push_str(&format!(
            " ({} the bounds [{lo}, {hi}])",
            if inside { "inside" } else { "OUTSIDE" }
        ));
    }
    if let (Some(a), Some(b)) = (s.l2, s.linf) {
        l.push_str(&format!(", error L2 {a:.4e}, max {b:.4e}"));
    }
    if let Some(e) = &s.energy {
        l.push_str(&format!(
            "; heat balance: conduction {:.6e} W, convection {:.6e} W, sources {:.6e} W, ",
            e.conduction, e.convection, e.source
        ));
        if let Some(st) = e.storage {
            let per: Vec<String> = e
                .storage_regions
                .iter()
                .map(|r| format!("{} {:.6e} W", r.region, r.storage))
                .collect();
            l.push_str(&format!(
                "storage in the last step {st:.6e} W ({}), ",
                per.join(", ")
            ));
        }
        l.push_str(&format!(
            "imbalance {:.3e} W ({:.2e} of the largest flow)",
            e.imbalance, e.relative
        ));
    }
    l
}

/// ENE007: the heat balance of a converged steady run that does not close
/// to 1 % of its largest heat flow. The scaled residual of the temperature
/// can be small while the solve is far from converged, e.g. in a conjugate
/// case, whose solid coefficients dominate the scaling. (The reported
/// boundary fluxes take the orthogonal part of the gradient, so a coarse
/// perturbed mesh shows a gap of up to a few 10⁻³; 1 % stays clear of it.)
pub fn heat_balance_finding(s: &ScalarSummary) -> Option<crate::diag::Diagnostic> {
    let e = s.energy.as_ref()?;
    (e.relative > 1e-2).then(|| {
        crate::diag::Diagnostic::new(
            crate::diag::ENE_BALANCE,
            format!(
                "the run converged, but its heat balance closes only to {:.2e} of the largest \
                 heat flow ({:.3e} W): the temperature is not converged although its residual \
                 is; tighten [solver.energy] tolerance and raise max_iter (conjugate cases: \
                 1e-10 and 300)",
                e.relative, e.imbalance
            ),
        )
    })
}

/// The exact solution of a scalar or of the temperature, with its key.
fn exact_of<'a>(case: &'a crate::case::FlowCase, name: &str) -> Option<(&'a Scalar, String)> {
    if let Some(e) = case.scalars.get(name).and_then(|s| s.exact.as_ref()) {
        return Some((e, format!("scalars.{name}.exact")));
    }
    if name == crate::energy::FIELD {
        return case
            .energy
            .as_ref()
            .and_then(|e| e.exact.as_ref())
            .map(|e| (e, "energy.exact".to_string()));
    }
    None
}

/// The checks of the scalars over the solves of a report: the largest
/// error and, for a mesh sweep, the observed orders.
pub fn check_sweep(
    rep: &mut crate::case::FlowReport,
    v: &crate::case::VerificationSection,
    levels: usize,
    log: &mut dyn FnMut(&str),
) {
    let Some(first) = rep.solves.first() else {
        return;
    };
    if let Some([lo, hi]) = v.scalar_range {
        let mut out = Vec::new();
        for (k, s) in rep.solves.iter().enumerate() {
            for x in &s.scalars {
                // a NaN extreme fails too
                if !(x.min >= lo && x.max <= hi) {
                    out.push(format!(
                        "scalar {}: [{:.6}, {:.6}] on mesh {k} leaves the accepted range [{lo}, {hi}]",
                        x.name, x.min, x.max
                    ));
                }
            }
        }
        rep.failures.extend(out);
    }
    let names: Vec<String> = first
        .scalars
        .iter()
        .filter(|s| s.l2.is_some())
        .map(|s| s.name.clone())
        .collect();
    for name in names {
        let e: Vec<(f64, f64)> = rep
            .solves
            .iter()
            .filter_map(|s| s.scalars.iter().find(|x| x.name == name))
            .filter_map(|x| x.l2.map(|l2| (x.h, l2)))
            .collect();
        if let Some(tol) = v.scalar_max_error
            && let Some(&(_, l2)) = e.last()
            // a NaN error fails too
            && !matches!(
                l2.partial_cmp(&tol),
                Some(std::cmp::Ordering::Less | std::cmp::Ordering::Equal)
            )
        {
            rep.failures.push(format!(
                "scalar {name}: L2 error {l2:.4e} on the finest mesh above {tol}"
            ));
        }
        if levels >= 2 && e.len() == rep.solves.len() {
            let orders: Vec<f64> = e
                .windows(2)
                .map(|w| (w[0].1 / w[1].1).ln() / (w[0].0 / w[1].0).ln())
                .collect();
            log(&format!("scalar {name}: observed orders {orders:.3?}"));
            if let (Some([lo, hi]), Some(&o)) = (v.scalar_order, orders.last())
                && !(lo..=hi).contains(&o)
            {
                rep.failures.push(format!(
                    "scalar {name}: observed order {o:.3} on the finest pair, outside [{lo}, {hi}]"
                ));
            }
            rep.scalar_orders.push((name, orders));
        }
    }
}
