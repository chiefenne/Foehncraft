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
        GasState {
            u,
            rho: self.sol.phys.rho_at(cell),
            mu: self.sol.phys.mu,
        }
    }
}

/// Units of `spray.csv`.
pub const CSV_UNITS: &str = "step -; time s; parcels -; injected_kg kg; in_flow_kg kg; \
     deposited_kg kg; escaped_kg kg; lost_kg kg; lost_parcels -; rebounds -; substeps -; \
     com_x com_y com_z m (centre of mass of the liquid in the flow); d32_m m (Sauter mean \
     diameter of the droplets in the flow); per injector <name>.parcels -, <name>.in_flow_kg kg, \
     <name>.penetration_m m (distance from the injector, along its axis or straight as [spray] \
     penetration says, within which the share [spray] penetration_fraction of the injector's \
     liquid in the flow lies), <name>.farthest_m m (the farthest parcel, outlier-sensitive)";

/// The header of `spray.csv`.
pub fn csv_head(sp: &foehncraft_particles::Spray) -> String {
    let mut h = String::from(
        "step,time,parcels,injected_kg,in_flow_kg,deposited_kg,escaped_kg,\
         lost_kg,lost_parcels,rebounds,substeps,com_x,com_y,com_z,d32_m",
    );
    for i in &sp.config().injectors {
        let n = &i.name;
        h.push_str(&format!(
            ",{n}.parcels,{n}.in_flow_kg,{n}.penetration_m,{n}.farthest_m"
        ));
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

/// One row of `spray.csv`.
pub fn csv_row(step: usize, sp: &foehncraft_particles::Spray, zone: &[u32]) -> String {
    let t = &sp.tally;
    let com = sp.centre_of_mass();
    let d32 = sp
        .moments(None, zone)
        .stats(&foehncraft_particles::HistSpec::default())
        .d32;
    let mut extra = format!(",{:e},{:e},{:e},{d32:e}", com[0], com[1], com[2]);
    for i in sp.injector_stats() {
        extra.push_str(&format!(
            ",{},{:e},{:e},{:e}",
            i.parcels, i.mass, i.penetration, i.farthest
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
) -> Result<(u64, Option<std::path::PathBuf>)> {
    let env = foehncraft_particles::Env { mesh, geo, loc };
    let lost0 = sp.tally.lost_parcels;
    let gas = FlowGas { sol, geo };
    sp.advance(&env, &gas, t_end - h, h)?;
    csv.row(&csv_row(n, sp, &mesh.cell_zone));
    let every = case.spray.write_every;
    let mut file = None;
    if every > 0 && n.is_multiple_of(every) {
        file = write_parcels(sp, out, name, n)?;
    }
    Ok((sp.tally.lost_parcels - lost0, file))
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
