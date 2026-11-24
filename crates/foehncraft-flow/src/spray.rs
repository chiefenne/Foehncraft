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
                format!("[spray] with the periodic boundary {name} is not available in this version"),
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
     deposited_kg kg; escaped_kg kg; lost_kg kg; lost_parcels -; rebounds -; substeps -";

/// The header of `spray.csv`.
pub const CSV_HEAD: &str = "step,time,parcels,injected_kg,in_flow_kg,deposited_kg,escaped_kg,\
     lost_kg,lost_parcels,rebounds,substeps";

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
pub fn csv_row(step: usize, sp: &foehncraft_particles::Spray) -> String {
    let t = &sp.tally;
    format!(
        "{step},{:e},{},{:e},{:e},{:e},{:e},{:e},{},{},{}",
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
/// parcel file when it is due. Returns the parcels lost in this step.
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
) -> Result<u64> {
    let env = foehncraft_particles::Env { mesh, geo, loc };
    let lost0 = sp.tally.lost_parcels;
    let gas = FlowGas { sol, geo };
    sp.advance(&env, &gas, t_end - h, h)?;
    csv.row(&csv_row(n, sp));
    let every = case.spray.write_every;
    if every > 0 && n.is_multiple_of(every) {
        write_parcels(sp, out, name, n)?;
    }
    Ok(sp.tally.lost_parcels - lost0)
}

/// Writes `spray/parcels_<step>.vtkhdf`.
pub fn write_parcels(
    sp: &foehncraft_particles::Spray,
    out: &std::path::Path,
    name: &str,
    n: usize,
) -> Result<()> {
    #[cfg(feature = "hdf5")]
    {
        let dir = out.join("spray");
        std::fs::create_dir_all(&dir)?;
        sp.write(&dir.join(format!("parcels_{n:06}.vtkhdf")), name, n as u64)?;
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = (sp, out, name, n);
    }
    Ok(())
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
