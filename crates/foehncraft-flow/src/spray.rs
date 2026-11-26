//! The sprays of a flow case (L12, SP1): the `[spray]` table checks, the gas
//! seen by the parcels, and the output of a run.
//!
//! One-way in this version (D240 to D246): the parcels feel the gas of the
//! last time step, the gas does not feel them. The parcel layer itself is
//! the crate `foehncraft-particles`.

use crate::bc::BoundaryCondition;
use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use crate::simplec::Solution;
use anyhow::Result;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_mesh::Mesh;
use foehncraft_mesh::geometry::Geometry;
use foehncraft_particles::{BoundaryMode, BoundaryRule, Gas, GasState, SpraySection};

fn err(code: &'static str, key: &str, message: String) -> anyhow::Error {
    diag::err(Diagnostic::new(code, message).key(key))
}

/// Reads the spray file of `[spray] include` and checks the values (SPR001).
pub fn load(s: &mut SpraySection) -> Result<()> {
    s.load()
        .map_err(|e| err(diag::SPRAY_CONFIG, "spray", format!("{e:#}")))
}

/// The checks of the `[spray]` table against the rest of the case.
pub fn validate(case: &FlowCase) -> Result<()> {
    let s = &case.spray;
    if s.is_empty() {
        return Ok(());
    }
    if s.dispersion != foehncraft_particles::Dispersion::None && case.turbulence.is_laminar() {
        return Err(err(
            diag::SPRAY_CONFIG,
            "spray.dispersion",
            "[spray] dispersion needs a turbulence model ([turbulence] model): a laminar gas has \
             no k and ε for the eddies"
                .into(),
        ));
    }
    match (&case.time, &s.steady) {
        (Some(_), Some(_)) => {
            return Err(err(
                diag::SPRAY_STEADY,
                "spray.steady",
                "[spray.steady] is the steady spray mode (D260): a case with [time] runs the \
                 parcels in time; leave out one of the two"
                    .into(),
            ));
        }
        (Some(_), None) => {
            for i in &s.injectors {
                if !(i.parcels_per_second.is_finite() && i.parcels_per_second > 0.0) {
                    return Err(err(
                        diag::SPRAY_CONFIG,
                        "spray.injector",
                        format!(
                            "[[spray.injector]] {}: parcels_per_second must be positive in a \
                             transient case",
                            i.name
                        ),
                    ));
                }
            }
        }
        (None, _) => {
            // the steady mode (D260 (6)): a mass flow that changes in
            // time has no steady state
            for i in &s.injectors {
                if i.start != 0.0 || i.end.is_some() {
                    return Err(err(
                        diag::SPRAY_STEADY,
                        "spray.injector",
                        format!(
                            "[[spray.injector]] {}: start and end have no meaning in a steady \
                             case (D260): the injector's mass flow must not change in time",
                            i.name
                        ),
                    ));
                }
            }
            if case.fluid_view.is_some() {
                return Err(err(
                    diag::SPRAY_UNSUPPORTED,
                    "spray",
                    "the steady spray mode does not work with solid regions (conjugate heat \
                     transfer) in this version"
                        .into(),
                ));
            }
            if s.coupling.two_way()
                && case.solver.algorithm(false) == crate::simplec::Algorithm::Coupled
            {
                return Err(err(
                    diag::SPRAY_UNSUPPORTED,
                    "spray.coupling",
                    "the steady spray mode with two-way coupling needs SIMPLE or SIMPLEC \
                     ([solver] algorithm), not the coupled algorithm, in this version"
                        .into(),
                ));
            }
        }
    }
    // the vapour of the include file is known once it is read: `create`
    // checks it again
    if s.evaporation.on() && (s.liquid.vapour.is_some() || s.include.is_none()) {
        check_evaporation(case)?;
    }
    if !case.mrf.is_empty() || !case.sliding.is_empty() {
        return Err(err(
            diag::SPRAY_UNSUPPORTED,
            "spray",
            "[spray] with [[mrf]] or [[sliding]] is not available in this version (parcels in \
             turning zones follow with SP8)"
                .into(),
        ));
    }
    for (name, b) in &case.boundary {
        if matches!(b, BoundaryCondition::Periodic { .. }) {
            return Err(err(
                diag::SPRAY_UNSUPPORTED,
                &format!("boundary.{name}"),
                format!(
                    "[spray] with the periodic boundary {name} is not available in this version"
                ),
            ));
        }
    }
    if s.coupling.two_way() {
        let refuse = |what: &str| {
            Err(err(
                diag::SPRAY_UNSUPPORTED,
                "spray.coupling",
                format!("two-way spray coupling does not work with {what} in this version"),
            ))
        };
        if case.fluid_view.is_some() {
            return refuse("solid regions (conjugate heat transfer)");
        }
        if case.time.is_some() && case.solver.algorithm(true) != crate::simplec::Algorithm::Piso {
            return refuse("a transient algorithm other than PISO");
        }
    }
    for name in s.boundary.keys() {
        if !case.boundary.contains_key(name) {
            return Err(err(
                diag::SPRAY_CONFIG,
                &format!("spray.boundary.{name}"),
                format!("[spray.boundary.{name}]: the case has no boundary of that name"),
            ));
        }
    }
    Ok(())
}

/// The checks of an evaporating spray (SP4, SP5) against the rest of the
/// case.
fn check_evaporation(case: &FlowCase) -> Result<()> {
    let s = &case.spray;
    let key = "spray.evaporation";
    let Some(v) = s.liquid.vapour.as_ref() else {
        return Err(err(
            diag::SPRAY_CONFIG,
            key,
            "[spray.evaporation] needs [spray.liquid] vapour, the species the liquid \
                 evaporates into (for water \"H2O\")"
                .into(),
        ));
    };
    if case.energy.is_none() {
        return Err(err(
            diag::SPRAY_CONFIG,
            key,
            "[spray.evaporation] needs the energy equation ([energy]): the droplets \
                 evaporate by the heat of the gas"
                .into(),
        ));
    }
    let Some(sec) = case.species.as_ref() else {
        return Err(err(
            diag::SPRAY_CONFIG,
            key,
            format!(
                "[spray.evaporation] needs the vapour as a species of the gas: [species] \
                     names = [\"{v}\", \"O2\", \"N2\"] (the closure species last)"
            ),
        ));
    };
    let pos = sec.names.iter().position(|n| n == v);
    if pos.is_none() || pos == Some(sec.names.len() - 1) {
        return Err(err(
            diag::SPRAY_CONFIG,
            "spray.liquid.vapour",
            format!(
                "[spray.liquid] vapour = \"{v}\" must be a transported species of \
                     [species] names, not the closure species (the last one)"
            ),
        ));
    }
    if case.energy.as_ref().is_some_and(|e| {
        e.thermodynamic_pressure == crate::density::ThermodynamicPressure::ClosedVolume
    }) {
        return Err(err(
            diag::SPRAY_UNSUPPORTED,
            key,
            "evaporation in a closed volume (thermodynamic_pressure = \"closed-volume\") \
                 is not available in this version: p₀ would have to follow the mass the \
                 vapour adds"
                .into(),
        ));
    }
    if case.fluid_view.is_some() {
        return Err(err(
            diag::SPRAY_UNSUPPORTED,
            key,
            "evaporating sprays with solid regions (conjugate heat transfer) are not \
                 available in this version"
                .into(),
        ));
    }
    Ok(())
}

/// What each boundary of the mesh does to a parcel: the choice of
/// `[spray.boundary.<name>]`, else by the boundary condition (walls:
/// `wall_mode`; slip walls, symmetry planes and wedge sides: mirror;
/// everything else lets the parcel leave).
pub fn rules(case: &FlowCase, mesh: &Mesh) -> Vec<BoundaryRule> {
    let s = &case.spray;
    mesh.boundaries
        .iter()
        .map(|b| {
            let own = s.boundary.get(&b.name);
            let default = match case.boundary.get(&b.name) {
                Some(BoundaryCondition::Wall { .. }) => s.wall_mode,
                Some(
                    BoundaryCondition::SlipWall
                    | BoundaryCondition::Symmetry
                    | BoundaryCondition::Wedge,
                ) => BoundaryMode::Mirror,
                _ => BoundaryMode::Escape,
            };
            BoundaryRule {
                mode: own.and_then(|o| o.mode).unwrap_or(default),
                e_normal: own
                    .and_then(|o| o.restitution_normal)
                    .unwrap_or(s.restitution_normal),
                e_tangential: own
                    .and_then(|o| o.restitution_tangential)
                    .unwrap_or(s.restitution_tangential),
            }
        })
        .collect()
}

/// The gravity vector of the case (zero without `[gravity]`).
pub fn gravity(case: &FlowCase) -> [f64; 3] {
    case.gravity.as_ref().map_or([0.0; 3], |g| g.vector)
}

/// The gas of a flow solution as the parcels see it: the velocity of the
/// cell with its least-squares gradient (D117 (2)), the cell density and
/// the viscosity.
pub struct FlowGas<'a> {
    /// The solution.
    pub sol: &'a Solution<'a>,
    /// The geometry of the mesh.
    pub geo: &'a Geometry,
    /// The gas temperature and properties per cell, with the energy
    /// equation (SP2b).
    pub heat: Option<&'a GasHeat>,
    /// The vapour of the liquid per cell (evaporation, SP4, SP5).
    pub vapour: Option<&'a GasVapour>,
}

/// How the gas of a case carries the vapour of the liquid (SP4): the
/// species field, the diffusivity model and the operating pressure.
#[derive(Debug, Clone)]
pub struct VapourLink {
    /// The species (`H2O`).
    pub name: String,
    /// Its index in `[species] names`.
    pub index: usize,
    /// Its scalar field (`Y_H2O`).
    pub field: String,
    /// ρ D = λ/(c_p Le) (`diffusion = "lewis"`), else μ/Sc.
    pub lewis: bool,
    /// Sc or Le of the vapour.
    pub number: f64,
    /// The operating pressure (Pa): p₀ of a constant thermodynamic
    /// pressure.
    pub p_op: f64,
}

impl VapourLink {
    /// The link of `case`, when its droplets evaporate.
    pub fn of(case: &FlowCase) -> Option<Self> {
        let s = &case.spray;
        if s.is_empty() || !s.evaporation.on() {
            return None;
        }
        let name = s.liquid.vapour.clone()?;
        let sec = case.species.as_ref()?;
        let index = sec.names.iter().position(|n| *n == name)?;
        let lewis = sec.diffusion == crate::species::DiffusionModel::Lewis;
        let number = if lewis {
            sec.lewis.of(&name, 1.0)
        } else {
            sec.schmidt.of(&name, 0.7)
        };
        Some(Self {
            field: crate::species::field_name(&name),
            name,
            index,
            lewis,
            number,
            p_op: case.fluid.operating_pressure.unwrap_or(101_325.0),
        })
    }
}

/// The vapour of the gas per cell as the evaporating parcels see it (SP4).
#[derive(Debug, Clone, Default)]
pub struct GasVapour {
    /// Vapour mass fraction.
    pub y: Vec<f64>,
    /// Molar mass of the gas without the vapour (kg/mol).
    pub m_gas: Vec<f64>,
    /// Diffusivity of the vapour (m²/s).
    pub diff: Vec<f64>,
    /// The thermodynamic pressure p₀ (Pa).
    pub p0: f64,
}

impl GasVapour {
    /// The vapour of the flow `fl` by `link`, with the gas properties of
    /// `heat` (`None` without species or energy equation).
    pub fn of(
        fl: &mut crate::simplec::Flow,
        link: &VapourLink,
        heat: Option<&GasHeat>,
    ) -> Result<Option<Self>> {
        let phys = fl.physics();
        let Some(sc) = fl.scalars.as_ref() else {
            return Ok(None);
        };
        Self::from_parts(
            sc,
            fl.dev.as_mut(),
            &phys,
            fl.density.as_deref(),
            link,
            heat,
        )
    }

    /// As [`GasVapour::of`] from the parts of a flow.
    pub fn from_parts(
        sc: &crate::scalar::Scalars,
        dev: &mut dyn crate::device::FlowDevice,
        phys: &crate::simplec::Physics,
        density: Option<&crate::density::Variable>,
        link: &VapourLink,
        heat: Option<&GasHeat>,
    ) -> Result<Option<Self>> {
        let Some(heat) = heat else {
            return Ok(None);
        };
        let Some(sp) = sc.species() else {
            return Ok(None);
        };
        let Some(y) = sc.field_values(dev, &link.field)? else {
            return Ok(None);
        };
        let Some(mv) = sp.res.data.get(link.index).map(|d| d.molar_mass) else {
            return Ok(None);
        };
        let n = y.len().min(sp.cell.len());
        let mut m_gas = Vec::with_capacity(n);
        let mut diff = Vec::with_capacity(n);
        for (c, &yc) in y.iter().enumerate().take(n) {
            let yv = yc.clamp(0.0, 1.0 - 1e-9);
            // 1/M = Y_v/M_v + (1 − Y_v)/M_g
            let inv = 1.0 / sp.cell[c].molar_mass - yv / mv;
            m_gas.push(if inv > 0.0 {
                (1.0 - yv) / inv
            } else {
                sp.cell[c].molar_mass
            });
            let rho = phys.rho_at(c);
            diff.push(if link.lewis {
                heat.k[c] / (rho * heat.cp[c] * link.number)
            } else {
                phys.mu / (rho * link.number)
            });
        }
        let p0 = density.map_or(link.p_op, |d| d.p0);
        Ok(Some(Self { y, m_gas, diff, p0 }))
    }
}

/// The gas per cell for the heat exchange of the parcels (SP2b).
#[derive(Debug, Clone, Default)]
pub struct GasHeat {
    /// Temperature (K).
    pub t: Vec<f64>,
    /// Specific heat (J/(kg K)).
    pub cp: Vec<f64>,
    /// Thermal conductivity (W/(m K)).
    pub k: Vec<f64>,
}

impl GasHeat {
    /// The gas of the energy equation of `fl`, if it has one.
    pub fn of(fl: &mut crate::simplec::Flow) -> Result<Option<Self>> {
        let Some(sc) = fl.scalars.as_ref() else {
            return Ok(None);
        };
        Ok(sc
            .energy_gas(fl.dev.as_mut())?
            .map(|(t, cp, k)| Self { t, cp, k }))
    }
}

impl Gas for FlowGas<'_> {
    fn at(&self, cell: usize, x: [f64; 3]) -> GasState {
        let c = self.geo.cell_centroid[cell];
        let d = [x[0] - c[0], x[1] - c[1], x[2] - c[2]];
        let mut u = self.sol.state.u.get(cell);
        for (i, ui) in u.iter_mut().enumerate() {
            let g = self.sol.state.grad_u[i].get(cell);
            *ui += g[0] * d[0] + g[1] * d[1] + g[2] * d[2];
        }
        let mut g = GasState::isothermal(u, self.sol.phys.rho_at(cell), self.sol.phys.mu);
        if let Some(h) = self.heat {
            g.t = h.t[cell];
            g.cp = h.cp[cell];
            g.k = h.k[cell];
        }
        if let Some(v) = self.vapour {
            g.y_v = v.y[cell];
            g.p0 = v.p0;
            g.m_gas = v.m_gas[cell];
            g.diff = v.diff[cell];
        }
        if let Some(t) = &self.sol.state.turbulence {
            // the turbulence of the cell for the dispersion (D243): ε of the
            // k-ε family, or ε = β* k ω with β* = C_μ = 0.09 (SST)
            g.tke = t.k[cell].max(0.0);
            g.eps = if t.epsilon {
                t.omega[cell]
            } else {
                0.09 * g.tke * t.omega[cell]
            }
            .max(0.0);
        }
        g
    }
}

/// Units of `spray.csv`.
pub const CSV_UNITS: &str = "step -; time s; parcels -; injected_kg kg; in_flow_kg kg; \
     deposited_kg kg; escaped_kg kg; lost_kg kg; lost_parcels -; rebounds -; substeps -; \
     com_x com_y com_z m (centre of mass of the liquid in the flow); d32_m m (Sauter mean \
     diameter of the droplets in the flow); t_mean_k K (mass mean temperature of the droplets \
     in the flow); per injector <name>.parcels -, <name>.in_flow_kg kg, \
     <name>.penetration_m m (distance from the injector, along its axis or straight as [spray] \
     penetration says, within which the share [spray] penetration_fraction of the injector's \
     liquid in the flow lies), <name>.farthest_m m (the farthest parcel, outlier-sensitive); \
     with [spray.ecn] (D315) <name>.liquid_length_<k>_m m (the farthest axial position where \
     the projected liquid volume across the spray falls to the k-th threshold of [spray.ecn] \
     plv) and, evaporating, <name>.vapour_penetration_m m (the farthest axial distance of a cell \
     centre with the fuel mass fraction at or above [spray.ecn] vapour_fraction); \
     two-way coupling (SP2): passes - (parcel passes of the step), coupling_lag - (largest change \
     of the gas velocity of a spray cell against the frozen state of the last pass, over the \
     larger of the largest slip velocity and the largest change of the gas velocity of a spray \
     cell in the step), unconverged - (1: the passes ended at max_passes above the \
     tolerance), gas_force_x/y/z N (momentum the gas received per second), gas_power_w W, \
     coupling_max - (largest Δt/τ of the gas momentum of a cell, β Δt / (ρ V)), coupling_share - \
     (share of the spray cells with Δt/τ above 1), liquid_fraction_max - (largest liquid volume \
     fraction of a cell), parcels_per_cell_max -, gas_heat_w W (heat the gas received per \
     second, SP2b), coupling_heat_max - (largest Δt/τ of the gas temperature of a cell, \
     G Δt / (ρ c_p V)); evaporation (SP5): evaporated_kg kg (liquid evaporated since the \
     start), evaporated_parcels - (parcels that evaporated completely); with two-way coupling \
     also gas_vapour_kg_s kg/s (vapour the gas received per second), vapour_lag_kg_s kg/s (mass \
     source of the continuity less gas_vapour_kg_s), coupling_vapour_max - (largest Δt/τ of the \
     vapour of a cell, β Δt / (ρ V)), limiter_cells - (cells where the saturation limiter of \
     D253 held vapour back in the last pass), limiter_kg_s kg/s (the vapour held back); breakup \
     (SP6): child_parcels - (child parcels the blobs shed since the start), rt_breakups - \
     (Rayleigh-Taylor breakups of drops since the start); KH-ACT: act_children - (of the child \
     parcels, those the nozzle turbulence or the cavitation stripped)";

/// The header of `spray.csv`.
pub fn csv_head(sp: &foehncraft_particles::Spray) -> String {
    let mut h = String::from(
        "step,time,parcels,injected_kg,in_flow_kg,deposited_kg,escaped_kg,\
         lost_kg,lost_parcels,rebounds,substeps,com_x,com_y,com_z,d32_m,t_mean_k",
    );
    for i in &sp.config().injectors {
        let n = &i.name;
        h.push_str(&format!(
            ",{n}.parcels,{n}.in_flow_kg,{n}.penetration_m,{n}.farthest_m"
        ));
        if let Some(e) = &sp.config().ecn {
            for k in 1..=e.plv.len() {
                h.push_str(&format!(",{n}.liquid_length_{k}_m"));
            }
            if sp.evaporates() {
                h.push_str(&format!(",{n}.vapour_penetration_m"));
            }
        }
    }
    if sp.config().coupling.two_way() {
        h.push_str(
            ",passes,coupling_lag,unconverged,gas_force_x,gas_force_y,gas_force_z,\
             gas_power_w,coupling_max,coupling_share,liquid_fraction_max,parcels_per_cell_max,\
             gas_heat_w,coupling_heat_max",
        );
    }
    if sp.evaporates() {
        h.push_str(",evaporated_kg,evaporated_parcels");
        if sp.config().coupling.two_way() {
            h.push_str(
                ",gas_vapour_kg_s,vapour_lag_kg_s,coupling_vapour_max,limiter_cells,limiter_kg_s",
            );
        }
    }
    if sp.config().breakup.on() {
        h.push_str(",child_parcels,rt_breakups");
        if sp.config().breakup.act() {
            h.push_str(",act_children");
        }
    }
    h
}

/// The sprays in `result.json` (D265, D267).
#[derive(Debug, Clone, serde::Serialize)]
pub struct SpraySummary {
    /// The mass book.
    pub tally: foehncraft_particles::Tally,
    /// Liquid in the flow at the end (kg).
    pub in_flow_kg: f64,
    /// Centre of mass of the liquid in the flow (m).
    pub centre_of_mass: [f64; 3],
    /// Per injector: parcels, liquid, penetration, farthest parcel.
    pub injectors: Vec<foehncraft_particles::InjectorStats>,
    /// Droplet statistics of the parcels in the flow at the end.
    pub in_flow: foehncraft_particles::Stats,
    /// How the numbers are defined.
    pub definitions: &'static str,
    /// The two-way coupling over the run (SP2), if any.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coupling: Option<CouplingSummary>,
}

/// The definitions written with the statistics (D265).
pub const DEFINITIONS: &str = "Dpq = (sum N d^p / sum N d^q)^(1/(p-q)) over the droplets N of the \
     parcels; Dv10, Dv50, Dv90: diameters below which 10, 50, 90 % of the liquid mass lies, \
     interpolated in the mass histogram of the datum (log bins in ln d); span = (Dv90 - Dv10) / \
     Dv50; velocities weighted by droplet number (mean, standard deviation) and by mass; at a \
     surface datum every crossing counts (both directions, flux-weighted), in a volume the \
     parcels inside; penetration: distance from the injector within which the share \
     penetration_fraction of the injector's liquid in the flow lies (axial: along the injector \
     axis; radial: straight distance); rates are net (along minus against the datum's \
     direction) over the report interval, the mean rate over the run with the relative \
     standard error of ten batch means of the interval rates";

/// The spray summary at the end of a run.
pub fn summary(sp: &foehncraft_particles::Spray, zone: &[u32]) -> SpraySummary {
    let spec = foehncraft_particles::HistSpec::default();
    SpraySummary {
        tally: sp.tally.clone(),
        in_flow_kg: sp.mass_in_flow(),
        centre_of_mass: sp.centre_of_mass(),
        injectors: sp.injector_stats(),
        in_flow: sp.moments(None, zone).stats(&spec),
        definitions: DEFINITIONS,
        coupling: None,
    }
}

/// The carrier gas of the evaporation film (D321): the species of the
/// start composition without the vapour `v`, by mass. An expression in the
/// start composition leaves the closure species alone as the carrier.
fn carrier(case: &FlowCase, v: &str) -> Result<Vec<(foehncraft_kernels::species::Species, f64)>> {
    let Some(sec) = case.species.as_ref() else {
        return Ok(Vec::new());
    };
    let res = crate::species::Resolved::of(sec)?;
    let ys = res.mass_fractions(&sec.initial, "species.initial")?;
    let index = |n: &str| res.names.iter().position(|m| m == n);
    let closure = res.names[res.names.len() - 1].clone();
    let numbers: Option<Vec<f64>> = ys.iter().map(|y| y.number()).collect();
    let mut out = Vec::new();
    match numbers {
        Some(ys) => {
            let rest: f64 = ys.iter().sum();
            for (n, y) in res.transported().iter().zip(&ys) {
                if n != v
                    && *y > 0.0
                    && let Some(i) = index(n)
                {
                    out.push((res.data[i], *y));
                }
            }
            if let Some(i) = index(&closure)
                && closure != v
                && rest < 1.0
            {
                out.push((res.data[i], 1.0 - rest));
            }
        }
        None => {
            if let Some(i) = index(&closure) {
                out.push((res.data[i], 1.0));
            }
        }
    }
    Ok(out)
}

/// The sprays of a run, or `None` without `[spray]`.
pub fn create(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    loc: &foehncraft_mesh::locate::Locator,
) -> Result<Option<foehncraft_particles::Spray>> {
    if case.spray.is_empty() {
        return Ok(None);
    }
    if case.spray.evaporation.on() {
        check_evaporation(case)?;
    }
    let env = foehncraft_particles::Env { mesh, geo, loc };
    let mut sp =
        foehncraft_particles::Spray::new(&case.spray, &env, rules(case, mesh), gravity(case))
            .map_err(|e| err(diag::SPRAY_CONFIG, "spray", format!("{e:#}")))?;
    // a two-dimensional case lies in the x-y plane: eddies without z
    sp.set_planar(crate::case::dimension(case) == 2);
    // the vapour of an evaporating liquid (SP4): its species data
    if let Some(v) = case.spray.liquid.vapour.as_ref()
        && case.spray.evaporation.on()
    {
        let rec = crate::species::database().get(v).ok_or_else(|| {
            err(
                diag::SPRAY_CONFIG,
                "spray.liquid.vapour",
                format!("[spray.liquid] vapour = \"{v}\" is not a species of the data file"),
            )
        })?;
        sp.set_vapour(Some(foehncraft_particles::Vapour {
            species: rec.kernel(),
        }));
        sp.set_carrier(carrier(case, v)?);
    }
    // the Peng-Robinson surface at the run's p₀ (D311): the operating
    // pressure, as `GasVapour` gives it (a closed volume is refused above)
    if let Some(link) = VapourLink::of(case) {
        sp.set_surface_pressure(link.p_op).map_err(|e| {
            err(
                diag::SPRAY_CONFIG,
                "spray.evaporation.peng_robinson",
                format!("{e:#}"),
            )
        })?;
    }
    Ok(Some(sp))
}

/// The state for a restart file.
pub fn state(spray: Option<&foehncraft_particles::Spray>) -> Result<Option<Vec<u8>>> {
    spray.map(|s| s.state()).transpose()
}

/// One row of `spray.csv`; with two-way coupling the passes, the
/// coupling lag and the health numbers of the step (D254) close the row.
pub fn csv_row(
    step: usize,
    sp: &foehncraft_particles::Spray,
    zone: &[u32],
    c: Option<&StepCoupling>,
) -> String {
    let t = &sp.tally;
    let com = sp.centre_of_mass();
    let d32 = sp
        .moments(None, zone)
        .stats(&foehncraft_particles::HistSpec::default())
        .d32;
    let t_mean = match sp.mass_in_flow() {
        m if m > 0.0 => sp.enthalpy() / (m * sp.config().liquid.specific_heat),
        _ => f64::NAN,
    };
    let mut extra = format!(",{:e},{:e},{:e},{d32:e},{t_mean:e}", com[0], com[1], com[2]);
    for i in sp.injector_stats() {
        extra.push_str(&format!(
            ",{},{:e},{:e},{:e}",
            i.parcels, i.mass, i.penetration, i.farthest
        ));
        if sp.config().ecn.is_some() {
            for l in &i.liquid_length {
                extra.push_str(&format!(",{l:e}"));
            }
            if sp.evaporates() {
                extra.push_str(&format!(",{:e}", i.vapour_penetration.unwrap_or(0.0)));
            }
        }
    }
    if let Some(c) = c {
        extra.push_str(&format!(
            ",{},{:e},{},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{},{:e},{:e}",
            c.passes,
            c.lag,
            u8::from(c.unconverged),
            c.gas_momentum[0],
            c.gas_momentum[1],
            c.gas_momentum[2],
            c.gas_power,
            c.coupling_max,
            c.coupling_share,
            c.liquid_fraction_max,
            c.parcels_per_cell_max,
            c.gas_heat,
            c.coupling_heat_max
        ));
    }
    if sp.evaporates() {
        extra.push_str(&format!(
            ",{:e},{}",
            sp.tally.evaporated_mass, sp.tally.evaporated_parcels
        ));
        if sp.config().coupling.two_way() {
            let z = StepCoupling::default();
            let c = c.unwrap_or(&z);
            extra.push_str(&format!(
                ",{:e},{:e},{:e},{},{:e}",
                c.gas_vapour,
                c.gas_mass - c.gas_vapour,
                c.coupling_vapour_max,
                c.limiter_cells,
                c.limiter_mass
            ));
        }
    }
    if sp.config().breakup.on() {
        extra.push_str(&format!(",{},{}", t.child_parcels, t.rt_breakups));
        if sp.config().breakup.act() {
            extra.push_str(&format!(",{}", t.act_children));
        }
    }
    format!(
        "{step},{:e},{},{:e},{:e},{:e},{:e},{:e},{},{},{}{extra}",
        sp.time,
        sp.parcels.len(),
        t.injected_mass,
        sp.mass_in_flow(),
        t.deposited_mass.iter().sum::<f64>(),
        t.escaped_mass.iter().sum::<f64>(),
        t.lost_mass,
        t.lost_parcels,
        t.rebounds,
        t.substeps
    )
}

/// Advances the parcels over the step that ended at `t_end` (size `h`) in
/// the gas of the solution, books the row of `spray.csv` and writes the
/// parcel file when it is due. Returns the parcels lost in this step and
/// the parcel file written.
#[allow(clippy::too_many_arguments)]
pub fn step(
    case: &FlowCase,
    sp: &mut foehncraft_particles::Spray,
    sol: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    loc: &foehncraft_mesh::locate::Locator,
    n: usize,
    t_end: f64,
    h: f64,
    out: &std::path::Path,
    name: &str,
    csv: &mut crate::monitor::History,
    heat: Option<&GasHeat>,
    vapour: Option<&GasVapour>,
) -> Result<(u64, Option<std::path::PathBuf>)> {
    let env = foehncraft_particles::Env { mesh, geo, loc };
    let lost0 = sp.tally.lost_parcels;
    let gas = FlowGas {
        sol,
        geo,
        heat,
        vapour,
    };
    sp.advance(&env, &gas, t_end - h, h)?;
    if let Some(v) = vapour {
        sp.set_vapour_penetration(&v.y, &geo.cell_centroid);
    }
    csv.row(&csv_row(n, sp, &mesh.cell_zone, None));
    let every = case.spray.write_every;
    let mut file = None;
    if every > 0 && n.is_multiple_of(every) {
        file = write_parcels(sp, out, name, n)?;
    }
    Ok((sp.tally.lost_parcels - lost0, file))
}

/// As [`step`] for a two-way step: the closing step of the passes `done`
/// (D250) instead of a one-way advance. Returns the parcels lost, the
/// parcel file written and the coupling of the step.
#[allow(clippy::too_many_arguments)]
pub fn step_coupled(
    case: &FlowCase,
    sp: &mut foehncraft_particles::Spray,
    fl: &mut crate::simplec::Flow,
    done: Passes,
    mesh: &Mesh,
    geo: &Geometry,
    loc: &foehncraft_mesh::locate::Locator,
    n: usize,
    out: &std::path::Path,
    name: &str,
    csv: &mut crate::monitor::History,
) -> Result<(u64, Option<std::path::PathBuf>, StepCoupling)> {
    let env = foehncraft_particles::Env { mesh, geo, loc };
    let lost0 = sp.tally.lost_parcels;
    let c = close(sp, &env, fl, done)?;
    if sp.config().ecn.is_some()
        && let Some(link) = VapourLink::of(case)
        && let Some(sc) = fl.scalars.as_ref()
        && let Some(y) = sc.field_values(fl.dev.as_mut(), &link.field)?
    {
        sp.set_vapour_penetration(&y, &geo.cell_centroid);
    }
    csv.row(&csv_row(n, sp, &mesh.cell_zone, Some(&c)));
    let every = case.spray.write_every;
    let mut file = None;
    if every > 0 && n.is_multiple_of(every) {
        file = write_parcels(sp, out, name, n)?;
    }
    Ok((sp.tally.lost_parcels - lost0, file, c))
}

/// The coupling over a run (D252, D254): health numbers in `result.json`.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct CouplingSummary {
    /// Coupled steps.
    pub steps: u64,
    /// Parcel passes over the run.
    pub passes: u64,
    /// Steps that ended at max_passes with the lag above the tolerance.
    pub unconverged_steps: u64,
    /// Largest coupling lag at the end of a step.
    pub lag_max: f64,
    /// Largest Δt/τ of the gas momentum of a spray cell.
    pub coupling_max: f64,
    /// Largest liquid volume fraction of a cell.
    pub liquid_fraction_max: f64,
    /// Most parcels in one cell.
    pub parcels_per_cell_max: u32,
    /// Momentum the gas received over the run (kg m/s).
    pub gas_impulse: [f64; 3],
    /// Heat the gas received over the run (J, SP2b).
    pub gas_heat: f64,
    /// Largest Δt/τ of the gas temperature of a spray cell.
    pub coupling_heat_max: f64,
    /// Vapour the gas received over the run (kg, SP5).
    #[serde(skip_serializing_if = "is_zero")]
    pub gas_vapour: f64,
    /// Steps in which the saturation limiter acted (D253).
    #[serde(skip_serializing_if = "is_zero_u64")]
    pub limiter_steps: u64,
    /// The vapour it held back over the run (kg).
    #[serde(skip_serializing_if = "is_zero")]
    pub limiter_mass: f64,
}

fn is_zero(x: &f64) -> bool {
    *x == 0.0
}

fn is_zero_u64(x: &u64) -> bool {
    *x == 0
}

impl CouplingSummary {
    /// Adds a step of size `dt`.
    pub fn add(&mut self, c: &StepCoupling, dt: f64) {
        self.steps += 1;
        self.passes += u64::from(c.passes);
        self.unconverged_steps += u64::from(c.unconverged);
        self.lag_max = self.lag_max.max(c.lag);
        self.coupling_max = self.coupling_max.max(c.coupling_max);
        self.liquid_fraction_max = self.liquid_fraction_max.max(c.liquid_fraction_max);
        self.parcels_per_cell_max = self.parcels_per_cell_max.max(c.parcels_per_cell_max);
        for k in 0..3 {
            self.gas_impulse[k] += dt * c.gas_momentum[k];
        }
        self.gas_heat += dt * c.gas_heat;
        self.coupling_heat_max = self.coupling_heat_max.max(c.coupling_heat_max);
        self.gas_vapour += dt * c.gas_vapour;
        if c.limiter_cells > 0 {
            self.limiter_steps += 1;
            self.limiter_mass += dt * c.limiter_mass;
        }
    }
}

/// Writes `spray/parcels_<step>.vtkhdf`; returns its path (`None` in a
/// build without HDF5).
pub fn write_parcels(
    sp: &foehncraft_particles::Spray,
    out: &std::path::Path,
    name: &str,
    n: usize,
) -> Result<Option<std::path::PathBuf>> {
    #[cfg(feature = "hdf5")]
    {
        let dir = out.join("spray");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join(format!("parcels_{n:06}.vtkhdf"));
        sp.write(&path, name, n as u64)?;
        Ok(Some(path))
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = (sp, out, name, n);
        Ok(None)
    }
}

/// The mass book at the end of the run: released = in the flow + stuck +
/// left + lost + evaporated, to round-off (SPR005).
pub fn finish(sp: &foehncraft_particles::Spray) -> Result<()> {
    let t = &sp.tally;
    let total = sp.mass_in_flow()
        + t.deposited_mass.iter().sum::<f64>()
        + t.escaped_mass.iter().sum::<f64>()
        + t.lost_mass
        + t.evaporated_mass;
    if (t.injected_mass - total).abs() > 1e-9 * t.injected_mass.max(1e-300) {
        return Err(err(
            diag::SPRAY_BOOK,
            "spray",
            format!(
                "the spray mass book does not close: released {:e} kg, accounted for {total:e} kg",
                t.injected_mass
            ),
        ));
    }
    Ok(())
}

/// The statistics of a parcel file (the tools `post.droplets` and
/// `post.penetration`, D267): the same sums as the run.
#[derive(Debug, Clone, serde::Serialize)]
pub struct FileStatistics {
    /// The parcel file.
    pub file: std::path::PathBuf,
    /// Time of the file (s).
    pub time: Option<f64>,
    /// Step of the file.
    pub step: Option<u64>,
    /// Droplet statistics of all parcels of the file (or of those in the
    /// volume datum asked for).
    pub statistics: foehncraft_particles::Stats,
    /// The histogram: lower bin edges (m), droplets and mass (kg) per bin.
    pub histogram: Histogram,
    /// Per injector of the case: parcels, liquid, penetration.
    pub injectors: Vec<foehncraft_particles::InjectorStats>,
    /// How the numbers are defined.
    pub definitions: &'static str,
}

/// A droplet size histogram.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Histogram {
    /// Bin edges (m), bins + 1.
    pub edges: Vec<f64>,
    /// Droplets per bin.
    pub droplets: Vec<f64>,
    /// Mass per bin (kg).
    pub mass: Vec<f64>,
    /// Droplets and mass below the first edge.
    pub under: [f64; 2],
    /// Droplets and mass above the last edge.
    pub over: [f64; 2],
}

/// The statistics of the parcel file `file` of a run of `case` (its
/// liquid, injectors and penetration settings), with the bins `spec`;
/// with `datum` only the parcels inside that volume datum of the case.
/// A field the file lacks is reported with the fields it has.
#[cfg(feature = "hdf5")]
pub fn file_statistics(
    case: &FlowCase,
    file: &std::path::Path,
    spec: foehncraft_particles::HistSpec,
    datum: Option<&str>,
) -> Result<FileStatistics> {
    let pf = foehncraft_mesh::vtkhdf::read_points(file)?;
    let field = |name: &str, comps: usize| -> Result<&Vec<f64>> {
        match pf.fields.get(name) {
            Some((v, c)) if *c == comps => Ok(v),
            _ => {
                let have: Vec<String> = pf
                    .fields
                    .iter()
                    .map(|(k, (_, c))| format!("{k} ({c})"))
                    .collect();
                Err(err(
                    diag::SPRAY_CONFIG,
                    "spray",
                    format!(
                        "{}: no point field `{name}` with {comps} component(s); the file has {}",
                        file.display(),
                        have.join(", ")
                    ),
                ))
            }
        }
    };
    let d = field("diameter", 1)?;
    let count = field("droplets", 1)?;
    let vel = field("velocity", 3)?;
    let id: Vec<u64> = field("parcel_id", 1)?.iter().map(|&x| x as u64).collect();
    let v: Vec<[f64; 3]> = vel.as_chunks::<3>().0.to_vec();
    let rho = pf
        .numbers
        .get("liquid_density")
        .copied()
        .unwrap_or(case.spray.liquid.density);
    let mut keep: Vec<bool> = vec![true; id.len()];
    if let Some(name) = datum {
        let names: Vec<&str> = case
            .datums
            .iter()
            .filter(|x| matches!(x.kind.as_str(), "box" | "cylinder"))
            .map(|x| x.name.as_str())
            .collect();
        let sec = case
            .datums
            .iter()
            .find(|x| x.name == name && matches!(x.kind.as_str(), "box" | "cylinder"))
            .ok_or_else(|| {
                err(
                    diag::DATUM_CONFIG,
                    "datum",
                    format!(
                        "the case has no box or cylinder datum `{name}`; it has {}",
                        if names.is_empty() {
                            "none".into()
                        } else {
                            names.join(", ")
                        }
                    ),
                )
            })?;
        let vol = match sec.kind.as_str() {
            "box" => foehncraft_particles::Volume::Box {
                lower: sec.lower.unwrap_or_default(),
                upper: sec.upper.unwrap_or_default(),
            },
            _ => foehncraft_particles::Volume::Cylinder {
                start: sec.start.unwrap_or_default(),
                end: sec.end.unwrap_or_default(),
                radius: sec.radius.unwrap_or_default(),
            },
        };
        for (i, k) in keep.iter_mut().enumerate() {
            *k = vol.contains(pf.points[i], 0, &[]);
        }
    }
    let pick =
        |a: &[f64]| -> Vec<f64> { a.iter().zip(&keep).filter(|x| *x.1).map(|x| *x.0).collect() };
    let ids: Vec<u64> = id
        .iter()
        .zip(&keep)
        .filter(|x| *x.1)
        .map(|x| *x.0)
        .collect();
    let vs: Vec<[f64; 3]> = v.iter().zip(&keep).filter(|x| *x.1).map(|x| *x.0).collect();
    let xs: Vec<[f64; 3]> = pf
        .points
        .iter()
        .zip(&keep)
        .filter(|x| *x.1)
        .map(|x| *x.0)
        .collect();
    // the mirror planes of the ECN projection the run wrote (D315)
    let mirrors = match pf.text.get("ecn_mirrors") {
        Some(t) => foehncraft_particles::spray::parse_mirrors(t).map_err(|e| {
            err(
                diag::SPRAY_CONFIG,
                "ecn_mirrors",
                format!("{}: {e:#}", file.display()),
            )
        })?,
        None => Vec::new(),
    };
    let (ds, cs) = (pick(d), pick(count));
    let m = foehncraft_particles::moments_of(&spec, &ids, &ds, &cs, &vs, rho);
    let masses: Vec<f64> = ds
        .iter()
        .zip(&cs)
        .map(|(d, c)| c * rho * std::f64::consts::PI / 6.0 * d * d * d)
        .collect();
    Ok(FileStatistics {
        file: file.to_path_buf(),
        time: pf.numbers.get("time").copied(),
        step: pf.numbers.get("step").map(|&s| s as u64),
        statistics: m.stats(&spec),
        histogram: Histogram {
            edges: spec.edges(),
            droplets: m.hist_n.clone(),
            mass: m.hist_m.clone(),
            under: m.under,
            over: m.over,
        },
        injectors: foehncraft_particles::injector_stats_mirrored(
            &case.spray,
            &ids,
            &xs,
            &masses,
            None,
            &mirrors,
        ),
        definitions: DEFINITIONS,
    })
}

/// Without HDF5 there are no parcel files to read.
#[cfg(not(feature = "hdf5"))]
pub fn file_statistics(
    _case: &FlowCase,
    _file: &std::path::Path,
    _spec: foehncraft_particles::HistSpec,
    _datum: Option<&str>,
) -> Result<FileStatistics> {
    anyhow::bail!("parcel files are VTKHDF: this build has no HDF5 (feature `hdf5`)")
}

/// The bins of a droplet size histogram.
pub use foehncraft_particles::HistSpec;

// ---------------------------------------------------------------------------
// Two-way coupling (SP2, D247 to D252)

/// The momentum sources of a spray in value-plus-slope form (D250), per
/// spray cell: S = E − β (u − u*), with E the exchange of the pass (N), β
/// the summed slope (kg/s) and u* the frozen cell velocity of the pass.
#[derive(Debug, Clone, Default)]
pub struct MomentumSource {
    /// Cells, ascending.
    pub cells: Vec<u32>,
    /// E per cell (N).
    pub e: Vec<[f64; 3]>,
    /// β per cell (kg/s).
    pub beta: Vec<f64>,
    /// u* per cell (m/s).
    pub ustar: Vec<[f64; 3]>,
}

impl MomentumSource {
    /// The sources of the per-cell sums `src` of a pass over a step `dt`,
    /// frozen at the cell velocities `u`.
    pub fn of(
        src: &foehncraft_particles::CellSources,
        dt: f64,
        u: &dyn Fn(usize) -> [f64; 3],
    ) -> Self {
        Self {
            cells: src.cells.clone(),
            e: src
                .impulse
                .iter()
                .map(|i| [i[0] / dt, i[1] / dt, i[2] / dt])
                .collect(),
            beta: src.slope.iter().map(|b| b / dt).collect(),
            ustar: src.cells.iter().map(|&c| u(c as usize)).collect(),
        }
    }
}

impl crate::simplec::Flow {
    /// Adds the momentum sources of a two-way coupled spray to the rows
    /// just assembled (relaxation `alpha`), as the porous sink does: β on
    /// the diagonal of each component (Patankar 1980), E + β u* on the
    /// right-hand side, and β to D, D^s and D^C of the pressure correction.
    pub(crate) fn add_spray_momentum(&mut self, alpha: f64, theta: f64) -> Result<()> {
        let Some(src) = self.spray_momentum.as_ref() else {
            return Ok(());
        };
        if src.cells.is_empty() {
            return Ok(());
        }
        let (f, w) = (self.f, self.w);
        let u: Vec<Vec<f64>> =
            f.u.all()
                .into_iter()
                .map(|b| self.dev.download(b))
                .collect::<Result<_>>()?;
        let mut diag: Vec<Vec<f64>> = w
            .diag
            .all()
            .into_iter()
            .map(|b| self.dev.download(b))
            .collect::<Result<_>>()?;
        let mut rhs: Vec<Vec<f64>> = w
            .rhs
            .all()
            .into_iter()
            .map(|b| self.dev.download(b))
            .collect::<Result<_>>()?;
        let mut d = self.dev.download(w.d)?;
        let mut ds = self.dev.download(w.ds)?;
        let mut dc = self.dev.download(w.dc)?;
        let vol = &self.disc.cell_volume;
        for (j, &c) in src.cells.iter().enumerate() {
            let c = c as usize;
            let b = src.beta[j];
            for i in 0..3 {
                diag[i][c] += b / alpha;
                rhs[i][c] += src.e[j][i] + b * src.ustar[j][i] + (b / alpha - b) * u[i][c];
            }
            let add = b / (alpha * vol[c]);
            d[c] = 1.0 / (1.0 / d[c] + add);
            dc[c] = 1.0 / (1.0 / dc[c] + add);
            ds[c] = 1.0 / (1.0 / ds[c] + add / theta);
        }
        for (b, x) in w.diag.all().into_iter().zip(&diag) {
            self.dev.upload(b, x)?;
        }
        for (b, x) in w.rhs.all().into_iter().zip(&rhs) {
            self.dev.upload(b, x)?;
        }
        self.dev.upload(w.d, &d)?;
        self.dev.upload(w.ds, &ds)?;
        self.dev.upload(w.dc, &dc)?;
        Ok(())
    }
}

/// The parcel passes of one time step of a two-way coupled spray (D252):
/// a pass from the old level in every outer iteration until the coupling
/// lag is below the tolerance or the largest number of passes is reached.
pub struct TwoWay<'a> {
    sp: &'a foehncraft_particles::Spray,
    env: foehncraft_particles::Env<'a>,
    geo: &'a Geometry,
    t0: f64,
    dt: f64,
    limit: u32,
    converged: bool,
    tolerance: f64,
    /// Passes run so far in this step.
    pub passes: u32,
    /// The coupling lag after each pass.
    pub lags: Vec<f64>,
    /// The last pass and its frozen cell velocities and temperatures.
    pass: Option<Frozen>,
    /// The cell velocities of the start of the step (by cell).
    u_old: Vec<(u32, [f64; 3])>,
    /// The cell temperatures of the start of the step (by cell; SP2b).
    t_old: Vec<(u32, f64)>,
    /// The vapour of an evaporating spray (SP5).
    vapour: Option<VapourLink>,
    /// The vapour mass fractions of the start of the step (by cell).
    y_old: Vec<(u32, f64)>,
    /// The saturation limiter on the last pass (D253): cells, mass held
    /// back per second.
    limited: (u32, f64),
}

/// A pass with the frozen gas state of its cells.
pub struct Frozen {
    /// The pass.
    pub pass: foehncraft_particles::Pass,
    /// u* per cell of the sources.
    pub ustar: Vec<[f64; 3]>,
    /// T* per cell of the sources (empty without an energy equation).
    pub tstar: Vec<f64>,
    /// Y* of the vapour per cell of the sources (empty without
    /// evaporation, SP5).
    pub ystar: Vec<f64>,
}

impl<'a> TwoWay<'a> {
    /// The passes of the step from `t0` of size `dt`.
    pub fn new(
        sp: &'a foehncraft_particles::Spray,
        env: foehncraft_particles::Env<'a>,
        geo: &'a Geometry,
        t0: f64,
        dt: f64,
        vapour: Option<VapourLink>,
    ) -> Self {
        let c = sp.config().coupling;
        let (limit, converged) = c.pass_limit();
        Self {
            sp,
            env,
            geo,
            t0,
            dt,
            limit,
            converged,
            tolerance: c.tolerance,
            passes: 0,
            lags: Vec::new(),
            pass: None,
            u_old: Vec::new(),
            t_old: Vec::new(),
            vapour: vapour.filter(|_| sp.evaporates()),
            y_old: Vec::new(),
            limited: (0, 0.0),
        }
    }

    /// The passes of the step, for [`close`].
    pub fn finish(mut self) -> Passes {
        Passes {
            passes: self.passes,
            lag: self.lag(),
            unconverged: self.unconverged(),
            last: self.pass.take(),
            limited: self.limited,
        }
    }

    /// The last coupling lag (0 before a pass).
    pub fn lag(&self) -> f64 {
        self.lags.last().copied().unwrap_or(0.0)
    }

    /// Whether the step ended at the largest number of passes with the lag
    /// above the tolerance (a health number, D252).
    pub fn unconverged(&self) -> bool {
        self.converged && self.lag() > self.tolerance
    }
}

/// The passes of one step ([`TwoWay::finish`]).
pub struct Passes {
    /// Passes run.
    pub passes: u32,
    /// Coupling lag of the last.
    pub lag: f64,
    /// Ended at the largest number with the lag above the tolerance.
    pub unconverged: bool,
    /// The last pass and its frozen gas state.
    pub last: Option<Frozen>,
    /// The saturation limiter on the last pass (D253): cells, kg/s.
    pub limited: (u32, f64),
}

/// The coupling lag of D252: the largest change of the gas velocity of a
/// spray cell against the frozen state of the pass, over a velocity scale
/// of the step: the larger of the largest slip velocity of the spray cells
/// (|impulse| / slope of a cell) and the largest change of the gas velocity
/// of a spray cell over the step (against `u_old`, the state the step
/// started from). Both are velocity differences, so the measure does not
/// depend on the frame; the second keeps it finite where the parcels have
/// come to rest in the gas. The first pass of a step has a lag of 1 or
/// less by construction.
pub fn coupling_lag(
    src: &foehncraft_particles::CellSources,
    ustar: &[[f64; 3]],
    u_old: &[[f64; 3]],
    u: &dyn Fn(usize) -> [f64; 3],
) -> f64 {
    let norm = |a: [f64; 3], b: [f64; 3]| {
        let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    };
    let (mut du, mut scale, mut umax) = (0.0f64, 0.0f64, 0.0f64);
    for (j, &c) in src.cells.iter().enumerate() {
        let un = u(c as usize);
        du = du.max(norm(un, ustar[j]));
        scale = scale.max(norm(un, u_old[j]));
        let i = src.impulse[j];
        if src.slope[j] > 0.0 {
            scale = scale.max(norm(i, [0.0; 3]) / src.slope[j]);
        }
        umax = umax.max(norm(un, [0.0; 3]));
    }
    du / scale.max(1e-12 * (1.0 + umax))
}

/// The coupling lag of the temperature (SP2b), as [`coupling_lag`]: the
/// largest change of the gas temperature of a spray cell against the
/// frozen state of the pass, over the larger of the largest temperature
/// difference of parcels and gas (|heat| / heat slope of a cell) and the
/// largest change of the gas temperature of a spray cell over the step.
pub fn coupling_lag_heat(
    src: &foehncraft_particles::CellSources,
    tstar: &[f64],
    t_old: &[f64],
    t: &[f64],
) -> f64 {
    let (mut dt, mut scale, mut tmax) = (0.0f64, 0.0f64, 0.0f64);
    for (j, &c) in src.cells.iter().enumerate() {
        let tn = t[c as usize];
        dt = dt.max((tn - tstar[j]).abs());
        scale = scale.max((tn - t_old[j]).abs());
        if src.heat_slope[j] > 0.0 {
            scale = scale.max(src.heat[j].abs() / src.heat_slope[j]);
        }
        tmax = tmax.max(tn.abs());
    }
    dt / scale.max(1e-12 * (1.0 + tmax))
}

/// The coupling lag of the vapour (SP5), as [`coupling_lag_heat`]: the
/// largest change of the vapour mass fraction of a spray cell against the
/// frozen state of the pass, over the larger of the largest gap to the
/// saturation (|vapour| / vapour slope of a cell) and the largest change of
/// the vapour mass fraction of a spray cell over the step.
pub fn coupling_lag_vapour(
    src: &foehncraft_particles::CellSources,
    ystar: &[f64],
    y_old: &[f64],
    y: &[f64],
) -> f64 {
    let (mut dy, mut scale) = (0.0f64, 0.0f64);
    for (j, &c) in src.cells.iter().enumerate() {
        let yn = y[c as usize];
        dy = dy.max((yn - ystar[j]).abs());
        scale = scale.max((yn - y_old[j]).abs());
        if src.vapour_slope[j] > 0.0 {
            scale = scale.max(src.vapour[j].abs() / src.vapour_slope[j]);
        }
    }
    dy / scale.max(1e-12)
}

/// The vapour of a spray cell at the vapour mass fraction `y` (kg/s, SP5,
/// D250): S_v = (E − β (Y − Y*)) / Δt.
fn vapour_source(
    src: &foehncraft_particles::CellSources,
    j: usize,
    ystar: f64,
    y: f64,
    dt: f64,
) -> f64 {
    (src.vapour[j] - src.vapour_slope[j] * (y - ystar)) / dt
}

/// The gas that flows into each cell (D253): per cell the mass flow in
/// through its faces (kg/s) and the flow-weighted vapour mass fraction and
/// temperature of what comes in (the upwind cell of an interior face; the
/// cell's own state for a boundary face, which keeps the capacity on the
/// small side next to an inlet).
pub struct Inflow {
    m: Vec<f64>,
    y: Vec<f64>,
    t: Vec<f64>,
}

impl Inflow {
    /// The inflow of every cell of the solution `sol` with the gas states of `gas`.
    pub fn of(sol: &Solution, gas: &dyn Gas, geo: &Geometry) -> Self {
        let d = sol.disc;
        let n = d.n_cells;
        let state = |c: usize| {
            let g = gas.at(c, geo.cell_centroid[c]);
            (g.y_v, g.t)
        };
        let (mut m, mut y, mut t) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        let mut add = |c: usize, q: f64, (yv, tv): (f64, f64)| {
            m[c] += q;
            y[c] += q * yv;
            t[c] += q * tv;
        };
        for (f, &q) in sol.state.mass_flux.iter().enumerate().take(d.n_faces) {
            let o = d.owner[f] as usize;
            if o >= n {
                continue;
            }
            if f < d.n_interior {
                let nb = d.neighbour[f] as usize;
                if nb >= n {
                    continue;
                }
                if q > 0.0 {
                    add(nb, q, state(o));
                } else if q < 0.0 {
                    add(o, -q, state(nb));
                }
            } else if q < 0.0 {
                add(o, -q, state(o));
            }
        }
        for c in 0..n {
            if m[c] > 0.0 {
                y[c] /= m[c];
                t[c] /= m[c];
            }
        }
        Inflow { m, y, t }
    }
}

/// The factors of the saturation limiter of D253 per spray cell of
/// `src`: the vapour that brings the gas the cell holds (`store(c)`, kg; 0
/// in the steady mode) at its state and the gas that flows in over `tau`
/// (the time step; 1 in the steady mode, kg/s) at its own state to
/// adiabatic saturation is the capacity of the cell; a cell whose net
/// vapour is above it gets the factor of its evaporating records, else 1.
pub fn limiter_factors(
    sp: &foehncraft_particles::Spray,
    src: &foehncraft_particles::CellSources,
    gas: &dyn Gas,
    geo: &Geometry,
    inflow: &Inflow,
    tau: f64,
    store: &dyn Fn(usize) -> f64,
) -> Vec<f64> {
    src.cells
        .iter()
        .enumerate()
        .map(|(j, &c)| {
            let c = c as usize;
            if src.vapour[j] <= 0.0 {
                return 1.0;
            }
            let gs = gas.at(c, geo.cell_centroid[c]);
            let mut cap = sp.saturation_capacity(&gs, store(c));
            if inflow.m[c] > 0.0 {
                let gi = foehncraft_particles::GasState {
                    y_v: inflow.y[c],
                    t: inflow.t[c],
                    ..gs
                };
                cap += sp.saturation_capacity(&gi, inflow.m[c] * tau);
            }
            foehncraft_particles::limiter_factor(src.evaporated[j], src.vapour[j], cap)
        })
        .collect()
}

impl TwoWay<'_> {
    /// The sources of an evaporating spray at the current iterate (SP5):
    /// the vapour and the dilution in the species equations, the mass in
    /// the continuity (S_c = S_v of the iterate) and the heat of the pass
    /// with the enthalpy the vapour brings, in the T-form of the energy
    /// equation: q = Q − G (T − T*) − (β h)(Y − Y*)/Δt −
    /// S_v [h_v(T*) + c_pv (T − T*)], with h_v the vapour enthalpy of the
    /// parcels' reference (the liquid at 0 K) and S_v c_pv on the diagonal
    /// where S_v > 0.
    fn vapour_rows(&mut self, fl: &mut crate::simplec::Flow) -> Result<()> {
        let (Some(link), Some(fr)) = (self.vapour.as_ref(), self.pass.as_ref()) else {
            return Ok(());
        };
        if fr.ystar.len() != fr.pass.sources.cells.len() {
            return Ok(());
        }
        let Some(sc) = fl.scalars.as_ref() else {
            return Ok(());
        };
        let Some(y) = sc.field_values(fl.dev.as_mut(), &link.field)? else {
            return Ok(());
        };
        let t = sc.energy_temperature(fl.dev.as_mut())?;
        let phys = fl.physics();
        let vol = &self.geo.cell_volume;
        let src = &fr.pass.sources;
        let dt = self.dt;
        let mut mass = Vec::with_capacity(src.cells.len());
        let mut species = Vec::with_capacity(src.cells.len());
        let mut heat = Vec::with_capacity(src.cells.len());
        for (j, &c) in src.cells.iter().enumerate() {
            let cu = c as usize;
            let (ys, beta) = (fr.ystar[j], src.vapour_slope[j] / dt);
            let sv = vapour_source(src, j, ys, y[cu], dt);
            mass.push((c, sv));
            species.push((
                c,
                src.vapour[j] / dt + beta * ys,
                beta,
                sv,
                phys.rho_at(cu) * vol[cu],
            ));
            if fr.tstar.len() == src.cells.len() {
                let ts = fr.tstar[j];
                let g = src.heat_slope[j] / dt;
                let hv = self.sp.vapour_enthalpy(ts);
                let cpv = self.sp.vapour_cp(ts);
                let mut a = src.heat[j] / dt + g * ts
                    - src.vapour_slope_h[j] / dt * (y[cu] - ys)
                    - sv * (hv - cpv * ts);
                let mut gg = g;
                if sv > 0.0 {
                    gg += sv * cpv;
                } else {
                    a -= sv * cpv * t[cu];
                }
                heat.push((c, a, gg));
            }
        }
        fl.set_spray_mass(Some(mass))?;
        if let Some(sc) = fl.scalars.as_mut() {
            sc.set_spray_species(
                fl.dev.as_mut(),
                Some(crate::scalar::SpraySpecies {
                    vapour: link.field.clone(),
                    rows: species,
                }),
            )?;
            if !heat.is_empty() {
                sc.set_energy_spray(fl.dev.as_mut(), Some(heat))?;
            }
        }
        Ok(())
    }
}

fn fl_old<T: Copy>(old: &[(u32, T)], c: u32) -> Option<T> {
    old.binary_search_by_key(&c, |x| x.0).ok().map(|i| old[i].1)
}

/// Inserts the value of the cells of `cells` that `old` does not hold yet
/// (the current iterate: the cell was not reached by an earlier pass).
fn add_missing<T: Copy>(old: &mut Vec<(u32, T)>, cells: &[u32], now: &dyn Fn(u32) -> T) {
    for &c in cells {
        if let Err(k) = old.binary_search_by_key(&c, |x| x.0) {
            old.insert(k, (c, now(c)));
        }
    }
}

fn cell_velocity(st: &crate::simplec::State, c: usize) -> [f64; 3] {
    st.u.get(c)
}

impl crate::transient::OuterHook for TwoWay<'_> {
    fn max_outer(&self) -> usize {
        self.limit as usize
    }

    fn before_momentum(&mut self, fl: &mut crate::simplec::Flow, _k: usize) -> Result<()> {
        if self.passes >= self.limit {
            // the sources of the last pass stay (the flow iterates on); the
            // vapour's follow the gas (SP5)
            return self.vapour_rows(fl);
        }
        let heat = GasHeat::of(fl)?;
        let vap = match &self.vapour {
            Some(l) => GasVapour::of(fl, l, heat.as_ref())?,
            None => None,
        };
        let sol = Solution {
            state: fl.state()?,
            disc: &fl.disc,
            phys: fl.physics(),
        };
        let gas = FlowGas {
            sol: &sol,
            geo: self.geo,
            heat: heat.as_ref(),
            vapour: vap.as_ref(),
        };
        let mut pass = self.sp.pass(&self.env, &gas, self.t0, self.dt)?;
        if vap.is_some() {
            // the saturation limiter (D253): the gas of a cell and what
            // flows in over the step
            let inflow = Inflow::of(&sol, &gas, self.geo);
            let vol = &self.geo.cell_volume;
            let f = limiter_factors(
                self.sp,
                &pass.sources,
                &gas,
                self.geo,
                &inflow,
                self.dt,
                &|c| sol.phys.rho_at(c) * vol[c],
            );
            let cells = pass.sources.cells.clone();
            let (n, back) = self.sp.limit_vapour(&mut pass, &|c| {
                cells.binary_search(&c).map_or(1.0, |j| f[j])
            });
            self.limited = (n, back / self.dt);
        }
        let src = MomentumSource::of(&pass.sources, self.dt, &|c| cell_velocity(&sol.state, c));
        let ustar = src.ustar.clone();
        let cells = src.cells.clone();
        if self.passes == 0 {
            // the state the step starts from (D248: the first pass sees it)
            self.u_old = cells.iter().copied().zip(ustar.iter().copied()).collect();
        }
        add_missing(&mut self.u_old, &cells, &|c| {
            cell_velocity(&sol.state, c as usize)
        });
        drop(sol);
        fl.spray_momentum = Some(src);
        // the heat of the pass (SP2b): S = Q − G (T − T*), Q = heat / Δt,
        // G = heat slope / Δt
        let tstar: Vec<f64> = match &heat {
            Some(h) => {
                let ts: Vec<f64> = cells.iter().map(|&c| h.t[c as usize]).collect();
                if self.passes == 0 {
                    self.t_old = cells.iter().copied().zip(ts.iter().copied()).collect();
                }
                add_missing(&mut self.t_old, &cells, &|c| h.t[c as usize]);
                let s = &pass.sources;
                let rows: Vec<(u32, f64, f64)> = cells
                    .iter()
                    .enumerate()
                    .map(|(j, &c)| {
                        let g = s.heat_slope[j] / self.dt;
                        (c, s.heat[j] / self.dt + g * ts[j], g)
                    })
                    .collect();
                if let Some(sc) = fl.scalars.as_mut() {
                    sc.set_energy_spray(fl.dev.as_mut(), Some(rows))?;
                }
                ts
            }
            None => Vec::new(),
        };
        let ystar: Vec<f64> = match &vap {
            Some(v) => {
                let ys: Vec<f64> = cells.iter().map(|&c| v.y[c as usize]).collect();
                if self.passes == 0 {
                    self.y_old = cells.iter().copied().zip(ys.iter().copied()).collect();
                }
                add_missing(&mut self.y_old, &cells, &|c| v.y[c as usize]);
                ys
            }
            None => Vec::new(),
        };
        self.pass = Some(Frozen {
            pass,
            ustar,
            tstar,
            ystar,
        });
        self.passes += 1;
        self.vapour_rows(fl)
    }

    fn after_solve(&mut self, fl: &mut crate::simplec::Flow, _k: usize) -> Result<bool> {
        if self.pass.is_none() {
            return Ok(false);
        }
        let heat = GasHeat::of(fl)?;
        let st = fl.state()?;
        let Some(fr) = self.pass.as_ref() else {
            return Ok(false);
        };
        let cells = &fr.pass.sources.cells;
        let old: Vec<[f64; 3]> = cells
            .iter()
            .map(|&c| fl_old(&self.u_old, c).unwrap_or([0.0; 3]))
            .collect();
        let mut lag = coupling_lag(&fr.pass.sources, &fr.ustar, &old, &|c| {
            cell_velocity(&st, c)
        });
        if let Some(h) = &heat
            && fr.tstar.len() == cells.len()
        {
            let t_old: Vec<f64> = cells
                .iter()
                .map(|&c| fl_old(&self.t_old, c).unwrap_or(h.t[c as usize]))
                .collect();
            lag = lag.max(coupling_lag_heat(&fr.pass.sources, &fr.tstar, &t_old, &h.t));
        }
        if let Some(l) = &self.vapour
            && fr.ystar.len() == cells.len()
            && let Some(y) = fl
                .scalars
                .as_ref()
                .map(|sc| sc.field_values(fl.dev.as_mut(), &l.field))
                .transpose()?
                .flatten()
        {
            let y_old: Vec<f64> = cells
                .iter()
                .map(|&c| fl_old(&self.y_old, c).unwrap_or(y[c as usize]))
                .collect();
            lag = lag.max(coupling_lag_vapour(&fr.pass.sources, &fr.ystar, &y_old, &y));
        }
        self.lags.push(lag);
        if self.passes >= self.limit {
            return Ok(false);
        }
        Ok(!self.converged || lag > self.tolerance)
    }
}

/// What the closing step of a two-way step leaves for the ledger and the
/// health numbers (D244, D254).
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct StepCoupling {
    /// Passes of the step.
    pub passes: u32,
    /// Coupling lag of the last pass.
    pub lag: f64,
    /// The passes ended at the largest number with the lag above the
    /// tolerance.
    pub unconverged: bool,
    /// Momentum the gas received per second, Σ over the spray cells of
    /// E − β (u_new − u*) (N), summed in cell order.
    pub gas_momentum: [f64; 3],
    /// Power of that source on the gas, Σ u_new · S (W).
    pub gas_power: f64,
    /// Heat the gas received per second, Σ over the spray cells of
    /// Q − G (T_new − T*) (W), summed in cell order (SP2b).
    pub gas_heat: f64,
    /// Largest Δt/τ of the gas temperature in a spray cell,
    /// G Δt / (ρ c_p V) (SP2b, D254 (3)).
    pub coupling_heat_max: f64,
    /// Whether the step exchanged heat (an energy equation).
    #[serde(skip)]
    pub heat: bool,
    /// The momentum books of the parcels (kg m/s over the step).
    #[serde(skip)]
    pub books: foehncraft_particles::StepBooks,
    /// Largest Δt/τ of the gas momentum in a spray cell, β Δt / (ρ V) (D254 (3)).
    pub coupling_max: f64,
    /// Share of the spray cells with Δt/τ above 1.
    pub coupling_share: f64,
    /// Largest liquid volume fraction of a cell (D254 (5)).
    pub liquid_fraction_max: f64,
    /// Most parcels in one cell.
    pub parcels_per_cell_max: u32,
    /// Whether the step evaporated liquid into the gas (SP5).
    #[serde(skip)]
    pub vapour: bool,
    /// Vapour the gas received per second, Σ over the spray cells of
    /// (E − β (Y_new − Y*)) / Δt (kg/s), summed in cell order (SP5).
    pub gas_vapour: f64,
    /// The mass the continuity took per second (kg/s): S_v of the iterate
    /// before the last pressure solve; less `gas_vapour` it is the lag of
    /// the mass source (zero when the passes have converged).
    pub gas_mass: f64,
    /// The energy the gas received per second in the form of its own
    /// energy equation (W): the heat and the vapour of the pass with the
    /// enthalpy the vapour brings, as the balance of the energy equation
    /// sees it; less `gas_heat` it is the term `spray vapour enthalpy` of
    /// the ledger (the change from the liquid's to the gas's reference).
    pub gas_energy: f64,
    /// Largest Δt/τ of the vapour of a spray cell, β Δt / (ρ V) (SP5,
    /// D254 (3)).
    pub coupling_vapour_max: f64,
    /// Cells where the saturation limiter acted on the last pass (D253).
    pub limiter_cells: u32,
    /// The vapour it held back per second (kg/s).
    pub limiter_mass: f64,
}

/// The scalar field of the vapour of an evaporating spray, if set (SP5).
fn vapour_field(fl: &crate::simplec::Flow) -> Option<String> {
    fl.scalars
        .as_ref()?
        .spray_species()
        .map(|s| s.vapour.clone())
}

/// The closing step of a two-way step (D250): the parcels receive what the
/// gas gave up, the pass is committed, and the books of the step are
/// returned.
pub fn close(
    sp: &mut foehncraft_particles::Spray,
    env: &foehncraft_particles::Env,
    fl: &mut crate::simplec::Flow,
    done: Passes,
) -> Result<StepCoupling> {
    fl.spray_momentum = None;
    let Passes {
        passes,
        lag,
        unconverged,
        last,
        limited,
    } = done;
    let Some(Frozen {
        pass,
        ustar,
        tstar,
        ystar,
    }) = last
    else {
        return Ok(StepCoupling::default());
    };
    let heat = GasHeat::of(fl)?.filter(|_| tstar.len() == pass.sources.cells.len());
    // the vapour of the new level (SP5)
    let y_new: Option<Vec<f64>> = match (fl.scalars.as_ref(), vapour_field(fl)) {
        (Some(sc), Some(name)) if ystar.len() == pass.sources.cells.len() => {
            sc.field_values(fl.dev.as_mut(), &name)?
        }
        _ => None,
    };
    let (gas_mass, gas_energy) = match (fl.spray_mass.as_ref(), fl.scalars.as_ref()) {
        (Some(rows), Some(sc)) if y_new.is_some() => {
            let t = sc.energy_temperature(fl.dev.as_mut())?;
            (
                rows.iter().map(|r| r.1).exact_sum(),
                sc.energy_spray_total(&t).unwrap_or(0.0),
            )
        }
        _ => (0.0, 0.0),
    };
    let st = fl.state()?;
    let phys = fl.physics();
    let dt = pass.dt;
    let src = pass.sources.clone();
    let du = |c: u32| -> foehncraft_particles::Delta {
        match src.find(c) {
            Some(j) => {
                let u = st.u.get(c as usize);
                let dtemp = heat.as_ref().map_or(0.0, |h| h.t[c as usize] - tstar[j]);
                foehncraft_particles::Delta {
                    u: [u[0] - ustar[j][0], u[1] - ustar[j][1], u[2] - ustar[j][2]],
                    t: dtemp,
                    y: y_new.as_ref().map_or(0.0, |y| y[c as usize] - ystar[j]),
                }
            }
            None => foehncraft_particles::Delta::default(),
        }
    };
    let mut out = StepCoupling {
        passes,
        lag,
        unconverged,
        limiter_cells: limited.0,
        limiter_mass: limited.1,
        ..Default::default()
    };
    let (mut above, vol) = (0usize, &env.geo.cell_volume);
    out.heat = heat.is_some();
    out.vapour = y_new.is_some();
    let mut vap = Vec::new();
    for (j, &c) in src.cells.iter().enumerate() {
        let dl = du(c);
        let (d, dtemp) = (dl.u, dl.t);
        let u = st.u.get(c as usize);
        if out.vapour {
            vap.push(vapour_source(&src, j, 0.0, dl.y, dt));
            let m = phys.rho_at(c as usize) * vol[c as usize];
            out.coupling_vapour_max = out.coupling_vapour_max.max(src.vapour_slope[j] / m);
        }
        if let Some(h) = &heat {
            out.gas_heat +=
                (src.heat[j] - src.heat_slope[j] * dtemp - src.vapour_slope_h[j] * dl.y) / dt;
            let cap = phys.rho_at(c as usize) * h.cp[c as usize] * vol[c as usize];
            out.coupling_heat_max = out.coupling_heat_max.max(src.heat_slope[j] / cap);
        }
        let mut s = [0.0; 3];
        for k in 0..3 {
            s[k] = (src.impulse[j][k] - src.slope[j] * d[k]) / dt;
            out.gas_momentum[k] += s[k];
        }
        out.gas_power += u[0] * s[0] + u[1] * s[1] + u[2] * s[2];
        let cn = src.slope[j] / (phys.rho_at(c as usize) * vol[c as usize]);
        out.coupling_max = out.coupling_max.max(cn);
        if cn > 1.0 {
            above += 1;
        }
    }
    if !src.cells.is_empty() {
        out.coupling_share = above as f64 / src.cells.len() as f64;
    }
    if out.vapour {
        out.gas_vapour = vap.iter().exact_sum();
        out.gas_mass = gas_mass;
        out.gas_energy = gas_energy;
    }
    out.books = sp.commit(env, pass, Some(&du));
    for (c, n, v) in sp.cell_loading() {
        out.parcels_per_cell_max = out.parcels_per_cell_max.max(n);
        out.liquid_fraction_max = out.liquid_fraction_max.max(v / vol[c as usize]);
    }
    Ok(out)
}

/// Pair ids of the gas-parcel momentum exchange in the ledger (D244).
pub const PAIR_MOMENTUM: [u32; 3] = [2801, 2802, 2803];

/// Pair id of the gas-parcel heat exchange in the ledger (SP2b).
pub const PAIR_HEAT: u32 = 2804;

/// Pair id of the vapour, the gas-parcel mass exchange, in the ledger
/// (SP5).
pub const PAIR_VAPOUR: u32 = 2805;

/// The ledger rows of a two-way spray at a state (D119, D244): the gas side
/// of the momentum exchange on the momentum accounts of the domain (which
/// hold only this term, so their imbalance is the force of the walls and
/// boundaries on the gas), its power on the kinetic energy, and the
/// accounts of the parcel momentum, which close to round-off. `c` is the
/// coupling of the step that ended at this state (`None` at the start).
pub fn book(
    r: &mut crate::balance::Rates,
    sp: &foehncraft_particles::Spray,
    c: Option<&StepCoupling>,
    dt: f64,
    gas_momentum: [f64; 3],
    energy: bool,
) {
    use crate::balance::{Account, StepRule, Term, TermKind, names};
    let z = StepCoupling::default();
    let c = c.unwrap_or(&z);
    let per = |x: f64| if dt > 0.0 { x / dt } else { 0.0 };
    let b = &c.books;
    let content = sp.momentum();
    for k in 0..3 {
        let mut g = Account::new(names::MOMENTUM[k], "kg m/s");
        g.content = Some(gas_momentum[k]);
        g.terms.push(
            Term::volume("spray drag", TermKind::Exchange, c.gas_momentum[k])
                .rule(StepRule::End)
                .paired(PAIR_MOMENTUM[k]),
        );
        r.book(g);
        let mut a = Account::new(&format!("parcel {}", names::MOMENTUM[k]), "kg m/s");
        a.region = "parcels".into();
        a.content = Some(content[k]);
        for (name, v) in [
            ("weight", per(b.gravity[k])),
            ("injected", per(b.injected[k])),
            ("left the flow", per(b.removed[k])),
            ("walls", per(b.walls[k])),
        ] {
            a.terms
                .push(Term::volume(name, TermKind::Source, v).rule(StepRule::End));
        }
        a.terms.push(
            Term::volume("from the gas", TermKind::Exchange, -c.gas_momentum[k])
                .rule(StepRule::End)
                .paired(PAIR_MOMENTUM[k]),
        );
        r.book(a);
    }
    let mut ke = Account::new(names::KINETIC, "J");
    ke.terms
        .push(Term::volume("spray drag power", TermKind::Source, c.gas_power).rule(StepRule::End));
    r.book(ke);
    // the vapour (SP5): the gas mass receives what the parcels lose; the
    // continuity took the source of the iterate before the last pressure
    // solve, the difference is the lag of the mass source
    let evap = sp.evaporates();
    if evap {
        let mut m = Account::new(names::MASS, "kg");
        m.terms.push(
            Term::volume("spray vapour", TermKind::Exchange, c.gas_vapour)
                .rule(StepRule::End)
                .paired(PAIR_VAPOUR),
        );
        m.terms.push(
            Term::volume(
                "spray vapour lag",
                TermKind::Source,
                c.gas_mass - c.gas_vapour,
            )
            .rule(StepRule::End),
        );
        r.book(m);
        let mut a = Account::new("parcel mass", "kg");
        a.region = "parcels".into();
        a.content = Some(sp.mass_in_flow_by_id());
        for (name, v) in [
            ("injected", per(b.mass[3])),
            ("left the flow", per(b.mass[4])),
        ] {
            a.terms
                .push(Term::volume(name, TermKind::Source, v).rule(StepRule::End));
        }
        a.terms.push(
            Term::volume("to the gas", TermKind::Exchange, -c.gas_vapour)
                .rule(StepRule::End)
                .paired(PAIR_VAPOUR),
        );
        r.book(a);
    }
    if !energy {
        return;
    }
    // the heat (SP2b): the gas side on the thermal account of the energy
    // equation, the parcels' sensible enthalpy in an account of its own
    let mut th = Account::new(names::THERMAL, "J");
    th.terms.push(
        Term::volume("spray heat", TermKind::Exchange, c.gas_heat)
            .rule(StepRule::End)
            .paired(PAIR_HEAT),
    );
    if evap {
        // the heat and the vapour enthalpy arrive in the parcels' reference
        // (the liquid at 0 K); the energy equation of the gas has its own:
        // the difference, with the enthalpy the added mass carries out of
        // the boundaries
        th.terms.push(
            Term::volume(
                "spray vapour enthalpy",
                TermKind::Source,
                c.gas_energy - c.gas_heat,
            )
            .rule(StepRule::End),
        );
    }
    r.book(th);
    let [_, _, _, injected, removed] = b.enthalpy;
    let mut a = Account::new("parcel enthalpy", "J");
    a.region = "parcels".into();
    a.content = Some(sp.enthalpy());
    for (name, v) in [("injected", per(injected)), ("left the flow", per(removed))] {
        a.terms
            .push(Term::volume(name, TermKind::Source, v).rule(StepRule::End));
    }
    a.terms.push(
        Term::volume("from the gas", TermKind::Exchange, -c.gas_heat)
            .rule(StepRule::End)
            .paired(PAIR_HEAT),
    );
    r.book(a);
}

/// Σ ρ u V of the gas, in cell order (kg m/s).
pub fn gas_momentum(sol: &Solution, vol: &[f64]) -> [f64; 3] {
    let mut m = [0.0; 3];
    for (c, &v) in vol.iter().enumerate().take(sol.disc.n_cells) {
        let u = sol.state.u.get(c);
        let r = sol.phys.rho_at(c) * v;
        for k in 0..3 {
            m[k] += r * u[k];
        }
    }
    m
}
