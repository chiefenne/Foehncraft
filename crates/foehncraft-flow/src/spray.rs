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
    if case.time.is_none() {
        return Err(err(
            diag::SPRAY_STEADY,
            "spray",
            "[spray] needs a transient case ([time])".into(),
        ));
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
        if case.solver.algorithm(true) != crate::simplec::Algorithm::Piso {
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
     two-way coupling (SP2): passes - (parcel passes of the step), coupling_lag - (largest change \
     of the gas velocity of a spray cell against the frozen state of the last pass, over the \
     larger of the largest slip velocity and the largest change of the gas velocity of a spray \
     cell in the step), unconverged - (1: the passes ended at max_passes above the \
     tolerance), gas_force_x/y/z N (momentum the gas received per second), gas_power_w W, \
     coupling_max - (largest Δt/τ of the gas momentum of a cell, β Δt / (ρ V)), coupling_share - \
     (share of the spray cells with Δt/τ above 1), liquid_fraction_max - (largest liquid volume \
     fraction of a cell), parcels_per_cell_max -, gas_heat_w W (heat the gas received per \
     second, SP2b), coupling_heat_max - (largest Δt/τ of the gas temperature of a cell, \
     G Δt / (ρ c_p V))";

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
    }
    if sp.config().coupling.two_way() {
        h.push_str(
            ",passes,coupling_lag,unconverged,gas_force_x,gas_force_y,gas_force_z,\
             gas_power_w,coupling_max,coupling_share,liquid_fraction_max,parcels_per_cell_max,\
             gas_heat_w,coupling_heat_max",
        );
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
    let env = foehncraft_particles::Env { mesh, geo, loc };
    foehncraft_particles::Spray::new(&case.spray, &env, rules(case, mesh), gravity(case))
        .map(Some)
        .map_err(|e| err(diag::SPRAY_CONFIG, "spray", format!("{e:#}")))
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
) -> Result<(u64, Option<std::path::PathBuf>)> {
    let env = foehncraft_particles::Env { mesh, geo, loc };
    let lost0 = sp.tally.lost_parcels;
    let gas = FlowGas { sol, geo, heat };
    sp.advance(&env, &gas, t_end - h, h)?;
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
/// left + lost, to round-off (SPR005).
pub fn finish(sp: &foehncraft_particles::Spray) -> Result<()> {
    let t = &sp.tally;
    let total = sp.mass_in_flow()
        + t.deposited_mass.iter().sum::<f64>()
        + t.escaped_mass.iter().sum::<f64>()
        + t.lost_mass;
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
        injectors: foehncraft_particles::injector_stats(&case.spray, &ids, &xs, &masses, None),
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
}

/// A pass with the frozen gas state of its cells.
pub struct Frozen {
    /// The pass.
    pub pass: foehncraft_particles::Pass,
    /// u* per cell of the sources.
    pub ustar: Vec<[f64; 3]>,
    /// T* per cell of the sources (empty without an energy equation).
    pub tstar: Vec<f64>,
}

impl<'a> TwoWay<'a> {
    /// The passes of the step from `t0` of size `dt`.
    pub fn new(
        sp: &'a foehncraft_particles::Spray,
        env: foehncraft_particles::Env<'a>,
        geo: &'a Geometry,
        t0: f64,
        dt: f64,
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
        }
    }

    /// The passes of the step, for [`close`].
    pub fn finish(mut self) -> Passes {
        Passes {
            passes: self.passes,
            lag: self.lag(),
            unconverged: self.unconverged(),
            last: self.pass.take(),
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
            // the sources of the last pass stay (the flow iterates on)
            return Ok(());
        }
        let heat = GasHeat::of(fl)?;
        let sol = Solution {
            state: fl.state()?,
            disc: &fl.disc,
            phys: fl.physics(),
        };
        let gas = FlowGas {
            sol: &sol,
            geo: self.geo,
            heat: heat.as_ref(),
        };
        let pass = self.sp.pass(&self.env, &gas, self.t0, self.dt)?;
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
        self.pass = Some(Frozen { pass, ustar, tstar });
        self.passes += 1;
        Ok(())
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
    } = done;
    let Some(Frozen { pass, ustar, tstar }) = last else {
        return Ok(StepCoupling::default());
    };
    let heat = GasHeat::of(fl)?.filter(|_| tstar.len() == pass.sources.cells.len());
    let st = fl.state()?;
    let phys = fl.physics();
    let dt = pass.dt;
    let src = pass.sources.clone();
    let du = |c: u32| -> ([f64; 3], f64) {
        match src.find(c) {
            Some(j) => {
                let u = st.u.get(c as usize);
                let dtemp = heat.as_ref().map_or(0.0, |h| h.t[c as usize] - tstar[j]);
                (
                    [u[0] - ustar[j][0], u[1] - ustar[j][1], u[2] - ustar[j][2]],
                    dtemp,
                )
            }
            None => ([0.0; 3], 0.0),
        }
    };
    let mut out = StepCoupling {
        passes,
        lag,
        unconverged,
        ..Default::default()
    };
    let (mut above, vol) = (0usize, &env.geo.cell_volume);
    out.heat = heat.is_some();
    for (j, &c) in src.cells.iter().enumerate() {
        let (d, dtemp) = du(c);
        let u = st.u.get(c as usize);
        if let Some(h) = &heat {
            out.gas_heat += (src.heat[j] - src.heat_slope[j] * dtemp) / dt;
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
