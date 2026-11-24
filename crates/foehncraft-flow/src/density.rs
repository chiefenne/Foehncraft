//! Low-Mach variable density (stage 7c, D151).
//!
//! The density follows the temperature and the thermodynamic pressure p₀
//! through the fluid's equation of state, ρ = ρ(T, p₀): an ideal gas
//! ρ = p₀ M/(R T), a liquid by its density form of T (independent of p₀).
//! It does not depend on the dynamic pressure, so the pressure-correction
//! equation keeps its form (Rehm & Baum 1978; Paolucci 1982). After every
//! energy solve the cell and face densities are refreshed and uploaded; in
//! a closed volume p₀ follows from the fixed total mass, an exact reduction
//! (D31). In a time step the continuity rows take the mass source
//! V ∂ρ/∂t of the BDF formula and the momentum time term the density of
//! each time level.
//!
//! Switched on by `[energy] density = "variable"` or by
//! `[buoyancy] model = "variable-density"`; without it the density stays
//! one number and every path is bit for bit as before.

use crate::diag::{self, Diagnostic};
use crate::materials::Curve;
use anyhow::Result;
use foehncraft_core::exact::ExactSumIter as _;
use serde::{Deserialize, Serialize};

/// `[energy] density`: whether the density follows the temperature.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum DensityModel {
    /// One density, that of `[fluid] temperature` (default).
    #[default]
    Constant,
    /// ρ(T, p₀) from the fluid's density form in every cell and face
    /// (stage 7c, D151).
    Variable,
}

impl DensityModel {
    /// True for the default.
    pub fn is_default(&self) -> bool {
        *self == Self::Constant
    }
}

/// `[energy] thermodynamic_pressure`: how p₀ is found.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum ThermodynamicPressure {
    /// p₀ = `[fluid] operating_pressure`, constant (a domain open to an
    /// ambient through a pressure boundary).
    #[default]
    Constant,
    /// p₀(t) from the fixed mass of a closed volume (no inlet, no outlet).
    ClosedVolume,
}

impl ThermodynamicPressure {
    /// True for the default.
    pub fn is_default(&self) -> bool {
        *self == Self::Constant
    }
}

/// `[energy] pressure_work`: the work of the pressure in the energy
/// equation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum PressureWork {
    /// None.
    None,
    /// dp₀/dt per unit volume (W/m³): the heating by compression of a closed
    /// volume. The default with `thermodynamic_pressure = "closed-volume"`,
    /// where leaving it out would heat a gas with c_p instead of c_v.
    Thermodynamic,
}

/// The variable density of a flow: the equation of state, p₀, the cell and
/// face densities and their time levels.
#[derive(Debug, Clone)]
pub struct Variable {
    /// ρ(T) at the operating pressure.
    curve: Curve,
    /// p M/R of an ideal gas at the operating pressure (`None`: a liquid).
    pm_over_r: Option<f64>,
    /// The operating pressure (Pa) the curve was built for.
    p_op: f64,
    /// The thermodynamic pressure p₀ (Pa) of the last refresh.
    pub p0: f64,
    /// p₀ at the old time levels n and n − 1.
    p0_old: [f64; 2],
    /// A closed volume: p₀ from the mass `mass`.
    closed: bool,
    /// The mass of the closed volume (kg).
    pub mass: f64,
    /// The pressure work dp₀/dt enters the energy equation.
    pub pressure_work: bool,
    /// The specific heat c_p(T) of the fluid (the speed of sound).
    cp: Curve,
    /// ρ per cell of the last refresh.
    pub rho: Vec<f64>,
    /// ρ per face of the last refresh.
    pub rho_face: Vec<f64>,
    /// ρ per cell at the old time levels n and n − 1 (in a transient run).
    pub old: Option<[Vec<f64>; 2]>,
    /// The levels before the last `start_step` (taking a step back).
    saved: Option<([Vec<f64>; 2], [f64; 2])>,
    owner: Vec<u32>,
    neighbour: Vec<u32>,
    weight: Vec<f64>,
    volume: Vec<f64>,
    /// The factor M/M_fluid of a species mixture (stage SF1) per cell and
    /// per boundary face: ρ = p₀ M/(R T) with the mixture's molar mass.
    molar: Option<(Vec<f64>, Vec<f64>)>,
}

impl Variable {
    /// The density at temperature `t` and the current p₀.
    pub fn at(&self, t: f64) -> f64 {
        let r = self.curve.at(t);
        match self.pm_over_r {
            Some(_) => r * (self.p0 / self.p_op),
            None => r,
        }
    }

    /// True for an ideal gas (ρ ∝ p₀/T).
    pub fn is_ideal_gas(&self) -> bool {
        self.pm_over_r.is_some()
    }

    /// True in a closed volume (p₀ from the mass).
    pub fn is_closed(&self) -> bool {
        self.closed
    }

    /// The speed of sound √(γ R T / M) of an ideal gas at temperature `t`
    /// (m/s), γ = c_p/(c_p − R/M); `None` for a liquid.
    pub fn speed_of_sound(&self, t: f64) -> Option<f64> {
        let pmr = self.pm_over_r?;
        let r_spec = self.p_op / pmr; // R/M (J/(kg K))
        let cp = self.cp.at(t);
        let gamma = cp / (cp - r_spec);
        (gamma.is_finite() && gamma > 1.0).then(|| (gamma * r_spec * t).sqrt())
    }

    /// Sets the factor M/M_fluid of a species mixture per cell and per
    /// boundary face (stage SF1), for the next refresh.
    pub fn set_molar(&mut self, f: Option<(Vec<f64>, Vec<f64>)>) {
        self.molar = f;
    }

    /// The molar factor of cell `c` (1 without species).
    fn fc(&self, c: usize) -> f64 {
        self.molar
            .as_ref()
            .map_or(1.0, |(m, _)| m.get(c).copied().unwrap_or(1.0))
    }

    /// The molar factor of boundary face `fb` (1 without species).
    fn fb(&self, fb: usize) -> f64 {
        self.molar
            .as_ref()
            .map_or(1.0, |(_, m)| m.get(fb).copied().unwrap_or(1.0))
    }

    /// Recomputes p₀ (closed volume), the cell and the face densities from
    /// the cell temperatures `t` and the boundary-face temperatures `tb`.
    pub fn refresh(&mut self, t: &[f64], tb: &[f64]) {
        if self.closed && self.pm_over_r.is_some() {
            // M = Σ ρ V = (p₀/p_op) Σ ρ_op(T) V, exact
            let s: f64 = t
                .iter()
                .zip(&self.volume)
                .enumerate()
                .map(|(c, (&tc, v))| self.curve.at(tc) * self.fc(c) * v)
                .exact_sum();
            self.p0 = self.mass * self.p_op / s;
        }
        self.rho = t
            .iter()
            .enumerate()
            .map(|(c, &tc)| self.at(tc) * self.fc(c))
            .collect();
        let ni = self.neighbour.len();
        let nf = self.owner.len();
        self.rho_face = (0..nf)
            .map(|f| {
                if f < ni {
                    let w = self.weight[f];
                    w * self.rho[self.owner[f] as usize]
                        + (1.0 - w) * self.rho[self.neighbour[f] as usize]
                } else {
                    self.at(tb[f - ni]) * self.fb(f - ni)
                }
            })
            .collect();
    }

    /// A new time step: the current density and p₀ become level n, level n
    /// becomes n − 1 (the first step starts both from the start state).
    pub fn start_step(&mut self) {
        let prev = self.old.take();
        self.saved = prev.clone().map(|o| (o, self.p0_old));
        let n1 = match prev {
            Some([n, _]) => n,
            None => self.rho.clone(),
        };
        let p1 = if self.saved.is_some() {
            self.p0_old[0]
        } else {
            self.p0
        };
        self.old = Some([self.rho.clone(), n1]);
        self.p0_old = [self.p0, p1];
    }

    /// Takes the last `start_step` back (a rejected step).
    pub fn step_back(&mut self) {
        if let Some(o) = self.old.as_ref() {
            self.rho = o[0].clone();
            self.p0 = self.p0_old[0];
        }
        match self.saved.take() {
            Some((o, p)) => {
                self.old = Some(o);
                self.p0_old = p;
            }
            None => self.old = None,
        }
    }

    /// The mass source V ∂ρ/∂t per cell with the backward-difference
    /// coefficients `a` (a₀, −a₁, −a₂ as the momentum time terms) and Δt.
    pub fn mass_source(&self, a: [f64; 3], dt: f64) -> Option<Vec<f64>> {
        let [n, n1] = self.old.as_ref()?;
        Some(
            (0..self.rho.len())
                .map(|c| {
                    let mut d = a[0] * self.rho[c] - a[1] * n[c];
                    if a[2] != 0.0 {
                        d -= a[2] * n1[c];
                    }
                    d * self.volume[c] / dt
                })
                .collect(),
        )
    }

    /// dp₀/dt (Pa/s) with the coefficients `a` and Δt (0 outside a time
    /// step or with a constant p₀).
    pub fn dp0_dt(&self, a: [f64; 3], dt: f64) -> f64 {
        if !self.closed || self.old.is_none() {
            return 0.0;
        }
        let mut d = a[0] * self.p0 - a[1] * self.p0_old[0];
        if a[2] != 0.0 {
            d -= a[2] * self.p0_old[1];
        }
        d / dt
    }

    /// The mass in the domain, Σ ρ V (kg).
    pub fn total_mass(&self) -> f64 {
        self.rho
            .iter()
            .zip(&self.volume)
            .map(|(r, v)| r * v)
            .exact_sum()
    }
}

/// True when the case asks for the variable density.
pub fn wanted(case: &crate::case::FlowCase) -> bool {
    case.energy
        .as_ref()
        .is_some_and(|e| e.density == DensityModel::Variable)
        || case.buoyancy.model == crate::buoyancy::Model::VariableDensity
}

/// The density curve ρ(T) of the fluid at the operating pressure and its
/// owner label.
pub fn fluid_curve(case: &crate::case::FlowCase) -> Result<(Curve, String)> {
    use crate::materials::{MaterialSection, Property, curve, lookup};
    let env = case.env()?;
    let cond = case.fluid.conditions();
    let (m, owner) = match &case.fluid.material {
        Some(n) => (lookup(n, &case.materials)?, format!("materials.{n}")),
        None => (MaterialSection::default(), "fluid".to_string()),
    };
    let (p, owner) = match (&case.fluid.density, &m.density) {
        (Some(s), _) => (Property::from(s.clone()), "fluid".to_string()),
        (None, Some(p)) => (p.clone(), owner),
        (None, None) => (Property::Number(f64::NAN), owner),
    };
    Ok((curve(&p, "density", &owner, &m, &env, cond)?, owner))
}

/// The checks of the variable density (ENE009, BUO005, ENE006).
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let closed_asked = case
        .energy
        .as_ref()
        .is_some_and(|e| e.thermodynamic_pressure == ThermodynamicPressure::ClosedVolume);
    if !wanted(case) {
        if closed_asked {
            return Err(diag::err(
                Diagnostic::new(
                    diag::VD_SETTINGS,
                    "[energy] thermodynamic_pressure = \"closed-volume\" needs the variable \
                     density: add [energy] density = \"variable\"",
                )
                .key("energy.thermodynamic_pressure"),
            ));
        }
        return Ok(());
    }
    let Some(e) = case.energy.as_ref() else {
        return Err(diag::err(
            Diagnostic::new(
                diag::BUO_NO_EXPANSION,
                "the variable density follows the temperature: add the table [energy] (the \
                 energy equation)",
            )
            .key("buoyancy.model"),
        ));
    };
    if case.buoyancy.model == crate::buoyancy::Model::Boussinesq {
        return Err(diag::err(
            Diagnostic::new(
                diag::VD_SETTINGS,
                "[energy] density = \"variable\" with [buoyancy] model = \"boussinesq\": the \
                 Boussinesq model keeps the density constant; use model = \
                 \"variable-density\" for the buoyancy of a variable density",
            )
            .key("buoyancy.model"),
        ));
    }
    if case.buoyancy.model == crate::buoyancy::Model::VariableDensity && case.gravity.is_none() {
        return Err(diag::err(
            Diagnostic::new(
                diag::BUO_NO_GRAVITY,
                "[buoyancy] model = \"variable-density\" without [gravity]: add the table \
                 [gravity] with `vector` (m/s²), or use [energy] density = \"variable\" for a \
                 variable density without gravity",
            )
            .key("buoyancy.model"),
        ));
    }
    let (c, owner) = fluid_curve(case)?;
    if c.is_constant() {
        return Err(diag::err(
            Diagnostic::new(
                diag::VD_SETTINGS,
                format!(
                    "the variable density needs a density that depends on the temperature, but \
                     the density of {owner} is one number; give the fluid a density form of T \
                     (for a gas `{{ form = \"ideal-gas\" }}`), or a material that has one"
                ),
            )
            .key("fluid.density"),
        ));
    }
    let closed = e.thermodynamic_pressure == ThermodynamicPressure::ClosedVolume;
    let open = case.boundary.values().any(|b| b.is_open());
    if closed && open {
        return Err(diag::err(
            Diagnostic::new(
                diag::VD_SETTINGS,
                "[energy] thermodynamic_pressure = \"closed-volume\" in a domain with an inlet or \
                 an outlet: p₀ of an open domain is the operating pressure; remove the key, or \
                 close the domain",
            )
            .key("energy.thermodynamic_pressure"),
        ));
    }
    if !closed && !open {
        return Err(diag::err(
            Diagnostic::new(
                diag::BUO_CLOSED,
                "a variable density in a closed domain (no inlet, no outlet) with a constant \
                 thermodynamic pressure: heating would change the mass in the domain; set \
                 [energy] thermodynamic_pressure = \"closed-volume\"",
            )
            .key("energy.thermodynamic_pressure"),
        ));
    }
    if closed && c.ideal_gas().is_none() {
        return Err(diag::err(
            Diagnostic::new(
                diag::VD_SETTINGS,
                format!(
                    "a closed volume of a liquid ({owner}): its density does not depend on p₀, \
                     so a heated liquid in a rigid closed volume has no low-Mach solution; open \
                     the domain to an ambient (a pressure boundary)"
                ),
            )
            .key("energy.thermodynamic_pressure"),
        ));
    }
    let unsupported = if !case.scalars.is_empty() {
        Some("passive scalars")
    } else if case.run.restart.is_some() {
        Some("a restart")
    } else if case.boundary.values().any(|b| b.is_periodic()) {
        Some("a periodic mesh")
    } else if case.boundary.values().any(|b| b.is_convective_outlet()) {
        Some("a convective outlet")
    } else if case.boundary.values().any(|b| b.has_wall_mass_flux()) {
        Some("a wall with a mass flux")
    } else if !case.mrf.is_empty() || !case.sliding.is_empty() {
        Some("a rotating zone")
    } else {
        None
    };
    if let Some(what) = unsupported {
        return Err(diag::err(
            Diagnostic::new(
                diag::ENE_UNSUPPORTED,
                format!(
                    "the variable density (stage 7c) does not yet work with {what}; use a \
                     constant density, or leave {what} out"
                ),
            )
            .key("energy.density"),
        ));
    }
    Ok(())
}

/// The variable density of `case` on the discretisation `d` (`None`
/// without it), at the start temperatures `t` (cells) and `tb` (boundary
/// faces).
pub fn setup(
    case: &crate::case::FlowCase,
    d: &foehncraft_fv::disc::Discretisation,
    t: &[f64],
    tb: &[f64],
    molar: Option<(Vec<f64>, Vec<f64>)>,
) -> Result<Option<Variable>> {
    if !wanted(case) {
        return Ok(None);
    }
    let e = case.energy.as_ref().expect("checked by validate");
    let (curve, _) = fluid_curve(case)?;
    let env = case.env()?;
    let fl = crate::energy::fluid(case, &env)?;
    let cond = case.fluid.conditions();
    let closed = e.thermodynamic_pressure == ThermodynamicPressure::ClosedVolume;
    let pressure_work = match e.pressure_work {
        Some(PressureWork::Thermodynamic) => true,
        Some(PressureWork::None) => false,
        None => closed,
    };
    let ni = d.n_interior;
    let mut v = Variable {
        pm_over_r: curve.ideal_gas(),
        curve,
        p_op: cond.pressure,
        p0: cond.pressure,
        p0_old: [cond.pressure; 2],
        closed,
        mass: 0.0,
        pressure_work,
        cp: fl.cp,
        rho: Vec::new(),
        rho_face: Vec::new(),
        old: None,
        saved: None,
        owner: d.owner.clone(),
        neighbour: d.neighbour[..ni].to_vec(),
        weight: d.weight.clone(),
        volume: d.cell_volume[..d.n_cells].to_vec(),
        molar,
    };
    // the start state at the operating pressure: its mass is that of the
    // closed volume
    v.closed = false;
    v.refresh(t, tb);
    v.closed = closed;
    v.mass = v.total_mass();
    Ok(Some(v))
}

/// Attaches the variable density of `case` to the flow (after the energy
/// equation, before the buoyancy, so that the start temperature gives the
/// first density and force).
pub fn attach(fl: &mut crate::simplec::Flow, case: &crate::case::FlowCase) -> Result<()> {
    if !wanted(case) {
        return Ok(());
    }
    let Some(t) = fl.fluid_temperature()? else {
        anyhow::bail!("the variable density needs the temperature of the energy equation");
    };
    let tb = fl.fluid_boundary_temperature(&t)?;
    let molar = fl.scalars().and_then(|s| s.species_molar_factor());
    if let Some(v) = setup(case, &fl.disc, &t, &tb, molar)? {
        fl.density = Some(Box::new(v));
    }
    fl.refresh_density()
}

/// The variable density at the end of a solve (stage 7c), for the report.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Summary {
    /// The thermodynamic pressure p₀ (Pa).
    pub p0: f64,
    /// p₀ is that of a closed volume (from its mass).
    pub closed_volume: bool,
    /// The mass in the domain, Σ ρ V (kg), and at the start.
    pub mass: f64,
    /// See `mass`.
    pub mass_start: f64,
    /// The smallest and the largest cell density (kg/m³).
    pub rho_min: f64,
    /// See `rho_min`.
    pub rho_max: f64,
    /// The largest Mach number |u|/a of the cells (an ideal gas; `None` for
    /// a liquid).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mach_max: Option<f64>,
}

/// The summary of the variable density of `fl` (`None` without it), with
/// ENE004 into `warnings` above Mach 0.3.
pub fn summary(
    fl: &mut crate::simplec::Flow,
    warnings: &mut Vec<Diagnostic>,
) -> Result<Option<Summary>> {
    if fl.density.is_none() {
        return Ok(None);
    }
    let t = fl.fluid_temperature()?;
    let u = fl.f.u;
    let ux = fl.dev.download(u.x)?;
    let uy = fl.dev.download(u.y)?;
    let uz = fl.dev.download(u.z)?;
    let d = fl.density.as_deref().expect("density");
    let (lo, hi) = d
        .rho
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &r| {
            (a.min(r), b.max(r))
        });
    let mach_max = match (&t, d.is_ideal_gas()) {
        (Some(t), true) => Some(
            (0..t.len())
                .filter_map(|c| {
                    let s = (ux[c] * ux[c] + uy[c] * uy[c] + uz[c] * uz[c]).sqrt();
                    d.speed_of_sound(t[c]).map(|a| s / a)
                })
                .fold(0.0, f64::max),
        ),
        _ => None,
    };
    if let Some(m) = mach_max.filter(|&m| m > 0.3) {
        warnings.push(Diagnostic::new(
            diag::VD_MACH,
            format!(
                "the largest Mach number of the variable-density run is {m:.2} (above 0.3): the \
                 low-Mach form, in which the density follows only the temperature and p₀, loses \
                 its validity; a compressible solver is needed for this flow"
            ),
        ));
    }
    Ok(Some(Summary {
        p0: d.p0,
        closed_volume: d.is_closed(),
        mass: d.total_mass(),
        mass_start: d.mass,
        rho_min: lo,
        rho_max: hi,
        mach_max,
    }))
}
