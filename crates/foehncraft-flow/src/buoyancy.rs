//! Gravity, the Boussinesq buoyancy and the balanced-force interpolation
//! (stage 7, PR 7b1; D153, D154 of the stage 7 specification).
//!
//! The solved pressure is the reduced pressure p′ = p − ρ_ref g·(x − x_ref),
//! so that a fluid at rest at the reference density needs no pressure
//! gradient at all; the body force is b = (ρ − ρ_ref) g, for the Boussinesq
//! model
//!
//! b_P = [(ρ₀ − ρ_ref) − ρ₀ β (T_P − T_ref)] g.
//!
//! **Balanced force (D154).** The spec writes the momentum force of a cell
//! as −G̃_P, the least-squares operator of D24 applied to the reduced face
//! jumps J_f = (p′_N − p′_P) − o_f with the face offsets o_f = b̄_f·d. The
//! least-squares gradient is linear in the differences, so
//! G̃_P = ∇p′_P − B_P with B_P the same operator applied to the offsets
//! alone. The flow solver therefore keeps its pressure gradient and gets
//!
//! - the cell source B_P V in the momentum rows,
//! - the offset q_f = o_f − B̄_f·d subtracted from the Rhie–Chow jump
//!   (B̄_f interpolated with the weight of ∇p̄_f), so the face flux reads
//!   (J_f − G̃_f·d),
//! - the term e_b = (b_P − B_P)·r in the extrapolated boundary pressure,
//!   p_b = p_P + ∇p_P·r + e_b = p_P + (G̃_P + b_P)·r (a zero normal
//!   gradient of the reduced jump at walls and symmetry planes),
//!
//! all three computed here from the fluid temperature once per outer
//! iteration (and corrector) and uploaded to the device: the CPU and GPU
//! paths share them with no kernel of their own. If the discrete pressure
//! satisfies J_f = 0 on every face (discrete hydrostatics), the cell force
//! and the face flux are zero exactly, on every mesh.
//!
//! **Face force.** The offset takes b̄_f = ½ (b_P + b_N), not the
//! interpolation weight: for a temperature that is linear along g (a
//! stable stratification at rest) the mean of the two centre values is the
//! value at the midpoint of the line P–N projected on g, so the offsets are
//! the exact differences of the quadratic hydrostatic pressure and a
//! discrete hydrostatic state exists on any mesh (with the weight w it
//! would exist only where w = ½).
//!
//! Case file:
//!
//! ```toml
//! [gravity]
//! vector = [0.0, -9.81, 0.0]        # m/s²
//! reference_point = [0.0, 0.0, 0.0] # x_ref of the reduced pressure
//!
//! [buoyancy]
//! model = "boussinesq"              # "none" (default) | "boussinesq"
//! reference_temperature = 300.0     # T_ref (K); default [fluid] temperature
//! expansion_coefficient = 3.3e-3    # β (1/K); default from the material
//! # reference_density = 1.16        # ρ_ref (kg/m³); default ρ₀
//!
//! [boundary.outlet]
//! type = "static-pressure-outlet"
//! pressure = 0.0
//! pressure_kind = "reduced"         # (default) the value is p′; "static": p
//! ```

use crate::diag::{self, Diagnostic};
use anyhow::Result;
use serde::{Deserialize, Serialize};

/// `[gravity]`: the gravitational acceleration, seen by every module.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GravitySection {
    /// The acceleration g (m/s²), for example `[0.0, 0.0, -9.81]`.
    pub vector: [f64; 3],
    /// The reference point x_ref of the reduced pressure
    /// p′ = p − ρ_ref g·(x − x_ref) (m); default the origin.
    #[serde(default)]
    pub reference_point: [f64; 3],
}

/// The buoyancy model (D153).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Model {
    /// No buoyancy force: with a constant density gravity is then only the
    /// hydrostatic part of the static pressure.
    #[default]
    None,
    /// ρ = ρ₀ everywhere except in the body force
    /// −ρ₀ β (T − T_ref) g.
    Boussinesq,
    /// ρ(T, p₀) everywhere (stage 7c, D151): the body force (ρ − ρ_ref) g
    /// with the variable density of `[energy] density = "variable"`, which
    /// this model switches on.
    VariableDensity,
}

/// `[buoyancy]`: the buoyancy model and its constants.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BuoyancySection {
    /// `none` (default), `boussinesq` or `variable-density` (the body force
    /// (ρ − ρ_ref) g with ρ(T, p₀) in every cell, stage 7c).
    #[serde(default)]
    pub model: Model,
    /// T_ref (K) of the body force; default the `[fluid] temperature`, at
    /// which ρ₀ is evaluated.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_temperature: Option<f64>,
    /// β (1/K); default from the fluid's density form at T_ref
    /// (−(1/ρ) dρ/dT; 1/T_ref for an ideal gas).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expansion_coefficient: Option<f64>,
    /// ρ_ref (kg/m³) of the reduced pressure; default ρ₀, the density of
    /// the flow.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_density: Option<f64>,
}

impl BuoyancySection {
    /// True for the default table (nothing written).
    pub fn is_default(&self) -> bool {
        *self == Self::default()
    }
}

/// What the value of a pressure boundary means when gravity is set
/// (`[boundary.<name>] pressure_kind`).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum PressureKind {
    /// The reduced pressure p′ = p − ρ_ref g·(x − x_ref) (default): a value
    /// of 0 is an opening to an ambient at rest at ρ_ref.
    #[default]
    Reduced,
    /// The static pressure p, converted per face to
    /// p′_b = p_b − ρ_ref g·(x_f − x_ref).
    Static,
}

impl PressureKind {
    /// True for the default.
    pub fn is_default(&self) -> bool {
        *self == Self::Reduced
    }
}

/// The hydrostatic part of the reference state: p = p′ + ρ_ref g·(x − x_ref).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hydrostatic {
    /// g (m/s²).
    pub g: [f64; 3],
    /// x_ref (m).
    pub x_ref: [f64; 3],
    /// ρ_ref (kg/m³).
    pub rho_ref: f64,
}

impl Hydrostatic {
    /// ρ_ref g·(x − x_ref): static minus reduced pressure at `x`.
    pub fn head(&self, x: [f64; 3]) -> f64 {
        self.rho_ref
            * (self.g[0] * (x[0] - self.x_ref[0])
                + self.g[1] * (x[1] - self.x_ref[1])
                + self.g[2] * (x[2] - self.x_ref[2]))
    }
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The resolved constants of the Boussinesq model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Constants {
    /// g (m/s²).
    pub g: [f64; 3],
    /// ρ₀ (kg/m³), the density of the flow.
    pub rho0: f64,
    /// ρ_ref (kg/m³).
    pub rho_ref: f64,
    /// β (1/K).
    pub beta: f64,
    /// T_ref (K).
    pub t_ref: f64,
    /// Where β comes from (for the log).
    pub beta_source: &'static str,
    /// The variable-density model (stage 7c): the force per cell is
    /// (ρ_P − ρ_ref) g from the cell density; β only sets the buoyancy
    /// velocity of the pseudo-time step.
    pub variable: bool,
}

impl Constants {
    /// The body force per unit volume of a cell of density `rho` (N/m³),
    /// variable density.
    pub fn force_of_density(&self, rho: f64) -> [f64; 3] {
        let c = rho - self.rho_ref;
        [c * self.g[0], c * self.g[1], c * self.g[2]]
    }

    /// The body force per unit volume at temperature `t` (N/m³), Boussinesq.
    pub fn force(&self, t: f64) -> [f64; 3] {
        let c = (self.rho0 - self.rho_ref) - self.rho0 * self.beta * (t - self.t_ref);
        [c * self.g[0], c * self.g[1], c * self.g[2]]
    }
}

/// The hydrostatic reference of a case: `Some` when `[gravity]` is set.
pub fn hydrostatic(case: &crate::case::FlowCase) -> Option<Hydrostatic> {
    let g = case.gravity.as_ref()?;
    let rho0 = case.rho();
    Some(Hydrostatic {
        g: g.vector,
        x_ref: g.reference_point,
        rho_ref: case.buoyancy.reference_density.unwrap_or(rho0),
    })
}

/// The resolved Boussinesq constants of a case (`None` without the model).
pub fn constants(case: &crate::case::FlowCase) -> Result<Option<Constants>> {
    if case.buoyancy.model == Model::None {
        return Ok(None);
    }
    let g = case.gravity.as_ref().map(|g| g.vector).unwrap_or([0.0; 3]);
    let rho0 = case.rho();
    let t_ref = case
        .buoyancy
        .reference_temperature
        .unwrap_or_else(|| case.fluid.conditions().temperature);
    let (beta, beta_source) = match case.buoyancy.expansion_coefficient {
        Some(b) => (b, "case"),
        None => (material_beta(case, t_ref)?, "material"),
    };
    Ok(Some(Constants {
        g,
        rho0,
        rho_ref: case.buoyancy.reference_density.unwrap_or(rho0),
        beta,
        t_ref,
        beta_source,
        variable: case.buoyancy.model == Model::VariableDensity,
    }))
}

/// β = −(1/ρ) dρ/dT of the fluid's density form at `t` (a central
/// difference over ±0.05 K of the material curve; exactly 1/T for an ideal
/// gas up to the difference); BUO003 for a constant density form.
fn material_beta(case: &crate::case::FlowCase, t: f64) -> Result<f64> {
    use crate::materials::{MaterialSection, Property, curve, lookup};
    let env = case.env()?;
    let cond = case.fluid.conditions();
    let (m, owner) = match &case.fluid.material {
        Some(n) => (lookup(n, &case.materials)?, format!("materials.{n}")),
        None => (MaterialSection::default(), "fluid".to_string()),
    };
    let p = match (&case.fluid.density, &m.density) {
        (Some(s), _) => Property::from(s.clone()),
        (None, Some(p)) => p.clone(),
        (None, None) => Property::Number(f64::NAN),
    };
    let c = curve(&p, "density", &owner, &m, &env, cond)?;
    let h = 0.05;
    let (lo, mid, hi) = (c.at(t - h), c.at(t), c.at(t + h));
    let beta = -(hi - lo) / (2.0 * h * mid);
    if c.is_constant() || !beta.is_finite() || beta == 0.0 {
        return Err(diag::err(
            Diagnostic::new(
                diag::BUO_NO_EXPANSION,
                format!(
                    "the Boussinesq model needs the expansion coefficient β, but the density of \
                     {owner} does not depend on the temperature; give [buoyancy] \
                     expansion_coefficient (1/K) or a density form of T for the fluid"
                ),
            )
            .key("buoyancy.expansion_coefficient"),
        ));
    }
    Ok(beta)
}

/// The checks of `[gravity]` and `[buoyancy]` (BUO001, BUO003, ENE006).
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    if let Some(g) = &case.gravity {
        let n = norm(g.vector);
        if !(n.is_finite() && n > 0.0) || g.reference_point.iter().any(|x| !x.is_finite()) {
            return Err(diag::err(
                Diagnostic::new(
                    diag::BUO_NO_GRAVITY,
                    format!(
                        "[gravity] vector {:?} has no length (or a value is not finite); give \
                         the acceleration in m/s², for example [0.0, 0.0, -9.81], or remove \
                         the table",
                        g.vector
                    ),
                )
                .key("gravity.vector"),
            ));
        }
    }
    match case.buoyancy.model {
        Model::None => Ok(()),
        // the checks of the variable density are those of crate::density
        // (called after this one)
        Model::VariableDensity => Ok(()),
        Model::Boussinesq => {
            if case.gravity.is_none() {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::BUO_NO_GRAVITY,
                        "[buoyancy] model = \"boussinesq\" without [gravity]: add the table \
                         [gravity] with `vector` (m/s²)",
                    )
                    .key("buoyancy.model"),
                ));
            }
            if case.energy.is_none() {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::BUO_NO_EXPANSION,
                        "the Boussinesq model needs the temperature: add the table [energy] \
                         (the energy equation)",
                    )
                    .key("buoyancy.model"),
                ));
            }
            if let Some(r) = case.buoyancy.reference_density
                && !(r.is_finite() && r > 0.0)
            {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::BUO_NO_EXPANSION,
                        format!("[buoyancy] reference_density must be > 0 (kg/m³), is {r}"),
                    )
                    .key("buoyancy.reference_density"),
                ));
            }
            let c = constants(case)?.expect("Boussinesq");
            if !(c.t_ref.is_finite() && c.t_ref > 0.0) {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::BUO_NO_EXPANSION,
                        format!(
                            "[buoyancy] reference_temperature must be > 0 K, is {}",
                            c.t_ref
                        ),
                    )
                    .key("buoyancy.reference_temperature"),
                ));
            }
            if !c.beta.is_finite() {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::BUO_NO_EXPANSION,
                        format!(
                            "[buoyancy] expansion_coefficient must be finite, is {}",
                            c.beta
                        ),
                    )
                    .key("buoyancy.expansion_coefficient"),
                ));
            }
            Ok(())
        }
    }
}

/// The set-up findings of buoyancy: BUO002 (β ΔT above 0.1, the
/// Boussinesq model outside its range), BUO004 (gravity and the energy
/// equation without a buoyancy model), BUO006 (a static-pressure value over
/// a tall opening), and the note that the turbulence models have no
/// buoyancy production yet (7e1).
pub fn setup_findings(
    case: &crate::case::FlowCase,
    mesh: &foehncraft_mesh::Mesh,
    geo: &foehncraft_mesh::Geometry,
) -> Vec<crate::bc::Finding> {
    use crate::bc::{Finding, Level};
    let mut out = Vec::new();
    let Some(gs) = &case.gravity else {
        return out;
    };
    let finding = |level, code: &'static str, boundary: &str, message: String| Finding {
        level,
        code,
        boundary: (!boundary.is_empty()).then(|| boundary.to_string()),
        message,
    };
    if let Ok(Some(c)) = constants(case)
        && !c.variable
        && let Some((lo, hi)) = crate::energy::temperature_range(case)
    {
        let d = c.beta.abs()
            * (hi - lo)
                .max((hi - c.t_ref).abs())
                .max((lo - c.t_ref).abs());
        if d > 0.1 {
            out.push(finding(
                Level::Warning,
                diag::BUO_RANGE,
                "",
                format!(
                    "β ΔT = {d:.3} above 0.1 (β = {:.4e} 1/K, temperatures {lo} K to {hi} K, \
                     T_ref = {} K): the density changes by more than 10 %, outside the range of \
                     the Boussinesq model (Gray & Giorgini 1976); `[buoyancy] model = \
                     \"variable-density\"` (stage 7c) takes ρ(T, p₀) in every cell",
                    c.beta, c.t_ref
                ),
            ));
        }
    }
    if case.buoyancy.model == Model::None && case.energy.is_some() {
        out.push(finding(
            Level::Warning,
            diag::BUO_NONE,
            "",
            "[gravity] and the energy equation without a buoyancy model: the temperature \
             drives no flow (set [buoyancy] model = \"boussinesq\" for natural or mixed \
             convection)"
                .into(),
        ));
    }
    if case.buoyancy.model != Model::None && !case.turbulence.is_laminar() {
        out.push(finding(
            Level::Note,
            diag::BUO_NONE,
            "",
            "the turbulence model has no buoyancy production of k (G_b) yet (stage 7e1): the \
             mean buoyancy force is in the momentum equations, its effect on the turbulence \
             is not"
                .into(),
        ));
    }
    // BUO006: a static-pressure value over an opening with a height span
    let g = gs.vector;
    let gn = norm(g);
    let rho_ref = case.buoyancy.reference_density.unwrap_or(case.rho());
    for p in &mesh.boundaries {
        let Some(bc) = case.boundary.get(&p.name) else {
            continue;
        };
        if bc.pressure_kind() != Some(PressureKind::Static) || !static_number(bc) {
            continue;
        }
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for f in p.range() {
            let h = dot(geo.face_centroid[f], g) / gn;
            lo = lo.min(h);
            hi = hi.max(h);
        }
        let dp = rho_ref * gn * (hi - lo);
        let u = inlet_speed(case);
        let q = 0.5 * case.rho() * u * u;
        if dp > q && dp > 0.0 {
            out.push(finding(
                Level::Warning,
                diag::BUO_TALL_OPENING,
                &p.name,
                format!(
                    "static pressure over a height span of {:.3e} m: the hydrostatic \
                     difference ρ_ref g Δz = {dp:.3e} Pa is above the dynamic pressure \
                     {q:.3e} Pa, so a uniform static value drives a flow through the opening; \
                     a uniform value of the reduced pressure (pressure_kind = \"reduced\", \
                     the default) is an ambient at rest",
                    hi - lo
                ),
            ));
        }
    }
    out
}

/// True for a pressure boundary whose value is one number.
fn static_number(bc: &crate::bc::BoundaryCondition) -> bool {
    use crate::bc::BoundaryCondition as B;
    match bc {
        B::StaticPressureInlet { pressure, .. } | B::StaticPressureOutlet { pressure, .. } => {
            pressure.number().is_some()
        }
        B::TotalPressureInlet { total_pressure, .. }
        | B::TotalPressureOutlet { total_pressure, .. } => total_pressure.number().is_some(),
        _ => false,
    }
}

/// The largest inlet speed given as numbers (m/s; 0 without one).
fn inlet_speed(case: &crate::case::FlowCase) -> f64 {
    use crate::bc::BoundaryCondition as B;
    case.boundary
        .values()
        .filter_map(|bc| match bc {
            B::VelocityInlet { velocity, .. } => {
                velocity.as_ref().and_then(crate::expr::numbers).map(norm)
            }
            _ => None,
        })
        .fold(0.0, f64::max)
}

/// The balanced-force arrays of one flow (host side, D154).
#[derive(Debug, Clone)]
pub struct Balanced {
    /// The model constants.
    pub constants: Constants,
    n_interior: usize,
    owner: Vec<u32>,
    neighbour: Vec<u32>,
    weight: Vec<f64>,
    cell_ptr: Vec<u32>,
    cell_faces: Vec<u32>,
    /// x_N − x_P per interior face (with the periodic shift).
    d: Vec<[f64; 3]>,
    /// x_f − x_P per boundary face.
    r: Vec<[f64; 3]>,
    /// The least-squares weights of the pressure per cell-face entry.
    w: Vec<[f64; 3]>,
    /// 1 on boundary faces that do not enter the pressure gradient with
    /// their value (walls, symmetry planes of the least-squares gradient).
    zero_gradient: Vec<f64>,
    volume: Vec<f64>,
    /// The cell forces of the last refresh (N/m³), for the output.
    pub force: Vec<[f64; 3]>,
    /// The temperature difference of the case (K): the larger of the range
    /// of its prescribed temperatures and of the last cell temperatures.
    pub delta_t: f64,
}

/// The arrays of one refresh: B_P V per cell and component, q_f per face,
/// e_b per boundary face.
pub struct Arrays {
    /// B_P V (N) per component.
    pub source: [Vec<f64>; 3],
    /// q_f per face.
    pub jump: Vec<f64>,
    /// e_b per boundary face.
    pub pb: Vec<f64>,
    /// The plain face offsets o_f (the D-weighted Rhie–Chow form of a
    /// volume-of-fluid case takes them in place of q_f).
    pub offset: Vec<f64>,
}

impl Balanced {
    /// The arrays of the discretisation `d` with the pressure weights `w`
    /// and the zero-gradient mask of the pressure gradient.
    pub fn new(
        constants: Constants,
        disc: &foehncraft_fv::disc::Discretisation,
        w: &foehncraft_fv::disc::Soa3,
        zero_gradient: &[f64],
    ) -> Self {
        let ni = disc.n_interior;
        let nf = disc.n_faces;
        let c3 = |s: &foehncraft_fv::disc::Soa3, i: usize| [s.x[i], s.y[i], s.z[i]];
        let d = (0..ni)
            .map(|f| {
                let p = c3(&disc.cell_centroid, disc.owner[f] as usize);
                let mut n = c3(&disc.cell_centroid, disc.neighbour[f] as usize);
                if f >= disc.periodic {
                    let s = c3(&disc.shift, f - disc.periodic);
                    n = [n[0] + s[0], n[1] + s[1], n[2] + s[2]];
                }
                [n[0] - p[0], n[1] - p[1], n[2] - p[2]]
            })
            .collect();
        let r = (ni..nf)
            .map(|f| {
                let p = c3(&disc.cell_centroid, disc.owner[f] as usize);
                let x = c3(&disc.face_centroid, f);
                [x[0] - p[0], x[1] - p[1], x[2] - p[2]]
            })
            .collect();
        let w = (0..w.x.len()).map(|k| c3(w, k)).collect();
        Self {
            constants,
            n_interior: ni,
            owner: disc.owner.clone(),
            neighbour: disc.neighbour.clone(),
            weight: disc.weight.clone(),
            cell_ptr: disc.cell_ptr.clone(),
            cell_faces: disc.cell_faces.clone(),
            d,
            r,
            w,
            zero_gradient: zero_gradient.to_vec(),
            volume: disc.cell_volume.clone(),
            force: Vec::new(),
            delta_t: 0.0,
        }
    }

    /// The buoyancy velocity √(|g| |β| ΔT L) over the length `l` (m/s; 0
    /// for an isothermal fluid): the speed scale of a flow driven only by
    /// buoyancy, which bounds the automatic pseudo-time step of the coupled
    /// algorithm (an explicit buoyancy force with steps much longer than
    /// the buoyancy time scale does not settle).
    pub fn speed(&self, l: f64) -> f64 {
        let c = &self.constants;
        (norm(c.g) * c.beta.abs() * self.delta_t * l).sqrt()
    }

    /// The arrays for the cell temperatures `t` (Boussinesq) or the cell
    /// densities `rho` (variable density).
    pub fn arrays(&mut self, t: &[f64], rho: Option<&[f64]>) -> Arrays {
        self.arrays_with(t, rho, None)
    }

    /// [`Balanced::arrays`] with an extra offset per face added to the
    /// face offsets o_f (the surface tension σ κ_f (α_N − α_P) of a
    /// volume-of-fluid case, L3; zero on boundary faces).
    pub fn arrays_with(&mut self, t: &[f64], rho: Option<&[f64]>, extra: Option<&[f64]>) -> Arrays {
        let nc = self.volume.len();
        let ni = self.n_interior;
        let nf = self.owner.len();
        let b: Vec<[f64; 3]> = match rho {
            Some(r) => r
                .iter()
                .map(|&x| self.constants.force_of_density(x))
                .collect(),
            None => t.iter().map(|&x| self.constants.force(x)).collect(),
        };
        // face offsets o_f = ½ (b_P + b_N)·d, o_b = b_P·r
        let o: Vec<f64> = (0..nf)
            .map(|f| {
                let p = self.owner[f] as usize;
                if f < ni {
                    let n = self.neighbour[f] as usize;
                    let bf = [
                        0.5 * (b[p][0] + b[n][0]),
                        0.5 * (b[p][1] + b[n][1]),
                        0.5 * (b[p][2] + b[n][2]),
                    ];
                    dot(bf, self.d[f])
                } else {
                    dot(b[p], self.r[f - ni])
                }
            })
            .collect();
        let o: Vec<f64> = match extra {
            Some(e) => o.iter().zip(e).map(|(a, b)| a + b).collect(),
            None => o,
        };
        // B_P: the least-squares operator on the offsets, the faces in the
        // order of the gradient kernel
        let mut bc = vec![[0.0; 3]; nc];
        for (c, bcc) in bc.iter_mut().enumerate() {
            let (lo, hi) = (self.cell_ptr[c] as usize, self.cell_ptr[c + 1] as usize);
            let mut g = [0.0; 3];
            for k in lo..hi {
                let (f, is_nb) = foehncraft_kernels::unpack_face(self.cell_faces[k]);
                let dk = if f < ni {
                    if is_nb { -o[f] } else { o[f] }
                } else if self.zero_gradient[f - ni] != 0.0 {
                    0.0
                } else {
                    o[f]
                };
                let wk = self.w[k];
                g[0] += wk[0] * dk;
                g[1] += wk[1] * dk;
                g[2] += wk[2] * dk;
            }
            *bcc = g;
        }
        let jump: Vec<f64> = (0..nf)
            .map(|f| {
                let p = self.owner[f] as usize;
                if f < ni {
                    let n = self.neighbour[f] as usize;
                    let (w, wn) = (self.weight[f], 1.0 - self.weight[f]);
                    let bf = [
                        w * bc[p][0] + wn * bc[n][0],
                        w * bc[p][1] + wn * bc[n][1],
                        w * bc[p][2] + wn * bc[n][2],
                    ];
                    o[f] - dot(bf, self.d[f])
                } else {
                    o[f] - dot(bc[p], self.r[f - ni])
                }
            })
            .collect();
        let pb: Vec<f64> = (ni..nf)
            .map(|f| {
                let p = self.owner[f] as usize;
                let e = [b[p][0] - bc[p][0], b[p][1] - bc[p][1], b[p][2] - bc[p][2]];
                dot(e, self.r[f - ni])
            })
            .collect();
        let source = std::array::from_fn(|i| (0..nc).map(|c| bc[c][i] * self.volume[c]).collect());
        let lo = t.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = t.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !t.is_empty() && hi >= lo {
            self.delta_t = self.delta_t.max(hi - lo);
        }
        self.force = b;
        Arrays {
            source,
            jump,
            pb,
            offset: o,
        }
    }
}

/// Attaches the Boussinesq buoyancy of `case` to the flow (after the energy
/// equation, so that the start temperature gives the first force).
pub fn attach(fl: &mut crate::simplec::Flow, case: &crate::case::FlowCase) -> Result<()> {
    match constants(case)? {
        Some(c) => {
            fl.attach_buoyancy(c)?;
            if let (Some(b), Some((lo, hi))) =
                (fl.buoyancy.as_mut(), crate::energy::temperature_range(case))
            {
                b.delta_t = b.delta_t.max(hi - lo);
            }
            Ok(())
        }
        None => Ok(()),
    }
}
