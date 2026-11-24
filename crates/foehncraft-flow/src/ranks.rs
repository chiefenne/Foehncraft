//! The flow solver on several ranks, checked against one rank (MPI step
//! W5, D93).
//!
//! A rank solves on its part of the mesh with one layer of halo cells
//! ([`foehncraft_mesh::decompose::extract`]): the CPU flow device keeps
//! room for the halo in every cell buffer, fills it from the owners after
//! each operation that writes the own cells, takes only the own cells and
//! the faces of own owner cells into a reduction and sums exactly over all
//! ranks; the linear solvers and the multigrid run on ranks as in steps W3
//! and W4. The algorithm is the one of one rank, unchanged, so every rank
//! computes bitwise the values of its cells in the undecomposed run.
//!
//! [`solve_one`] and [`solve_part`] run a number of outer iterations or
//! time steps of a case (with its scalars, energy equation and turbulence
//! model) on the whole mesh and on one rank's part; [`check`] compares a
//! part with the whole cell by cell (by global cell identity) and the
//! residual histories entry by entry. [`check_cases`] are the cases of the
//! tests and of the example `mpi_flow`.

use crate::case::FlowCase;
use crate::simplec::{Flow, Target};
use anyhow::Result;
use foehncraft_mesh::decompose::{Ranks, SharedComm, extract};
use foehncraft_mesh::{Geometry, Mesh};

/// The fields of a solve on the own cells, in local order, and the
/// residual history.
#[derive(Debug, Clone)]
pub struct FlowSolve {
    /// Global identity of every own cell.
    pub cell_id: Vec<u64>,
    /// Velocity components.
    pub u: [Vec<f64>; 3],
    /// Pressure.
    pub p: Vec<f64>,
    /// Pressure gradient components.
    pub grad_p: [Vec<f64>; 3],
    /// Limiter of the velocity components.
    pub psi: [Vec<f64>; 3],
    /// The passive scalars and the energy (by name) and the turbulence
    /// variables (`k`, `omega`, `nu_t`).
    pub extra: Vec<(String, Vec<f64>)>,
    /// One line per outer iteration or time step: every residual and
    /// iteration count.
    pub history: Vec<String>,
    /// Wall time of the outer iterations or time steps (s), without the
    /// setup (not compared).
    pub seconds: f64,
}

/// The result of comparing a rank's part with the one-rank solve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RankCheck {
    /// The rank.
    pub rank: usize,
    /// Its own cells.
    pub cells: usize,
    /// Cells whose velocity, pressure, scalars or turbulence variables
    /// differ in any bit.
    pub differing_fields: usize,
    /// Cells whose pressure gradient or limiter differs in any bit.
    pub differing_gradients: usize,
    /// The residual histories are equal.
    pub same_history: bool,
}

impl RankCheck {
    /// Everything is bitwise the one-rank solve.
    pub fn bitwise(&self) -> bool {
        self.differing_fields == 0 && self.differing_gradients == 0 && self.same_history
    }
}

/// A case of the checks.
#[derive(Debug, Clone)]
pub struct CheckCase {
    /// What it covers.
    pub name: &'static str,
    /// The case (its `[mesh]` and `[solver]` sections are used).
    pub case: FlowCase,
    /// Outer iterations.
    pub iterations: usize,
}

impl CheckCase {
    /// The mesh of the case.
    pub fn mesh(&self) -> Result<Mesh> {
        Ok(crate::case::mesh_for(&self.case, None)?.0)
    }
}

/// The flow solver of `case` on `mesh` (whole or a rank's part) with its
/// passive scalars, energy equation, turbulence model and time stepping;
/// `wall_distance` of the cells of `mesh`, computed on the whole mesh.
struct Solver {
    fl: Flow,
    tr: Option<crate::transient::Transient>,
    dt: f64,
}

impl Solver {
    fn new(
        case: &FlowCase,
        mesh: &Mesh,
        geo: &Geometry,
        target: &Target,
        wall_distance: Option<Vec<f64>>,
    ) -> Result<Self> {
        let phys = crate::case::physics(case)?;
        let values = case.compiled_boundaries()?;
        let turning = case.turning();
        let t0 = case.time.map_or(0.0, |t| t.start_time);
        let at = crate::bc::At {
            analytic: case.analytic.as_ref(),
            values: &values,
            t: t0,
            turning: &turning,
            hydrostatic: crate::buoyancy::hydrostatic(case),
        };
        let arrays = crate::bc::arrays(mesh, geo, &case.boundary, phys.rho, at)?;
        let disc = crate::simplec::discretisation(mesh, geo, &arrays)?;
        let phys = crate::case::physics_on(case, &phys, mesh, geo)?;
        let mut fl = Flow::on(
            disc,
            arrays,
            phys,
            case.solver.clone(),
            case.time.is_some(),
            None,
            target,
        )?;
        crate::scalar::attach(&mut fl, case, mesh, geo)?;
        crate::density::attach(&mut fl, case)?;
        crate::turbulence::attach_with(&mut fl, case, mesh, geo, wall_distance)?;
        let (tr, dt) = match case.time {
            Some(t) => (
                Some(crate::transient::Transient::new(
                    &mut fl,
                    t.scheme,
                    t.start_time,
                )?),
                t.time_step,
            ),
            None => (None, 0.0),
        };
        Ok(Self { fl, tr, dt })
    }

    /// Start fields of the own cells.
    fn start(&mut self, u: &foehncraft_fv::disc::Soa3, p: &[f64]) -> Result<()> {
        self.fl.initialise(u, p)?;
        if let Some(tr) = self.tr.as_mut() {
            tr.start(&mut self.fl)?;
        }
        Ok(())
    }

    /// One outer iteration or one time step: its record (every residual,
    /// iteration count and Courant number, exactly) and the largest
    /// residual.
    fn advance(&mut self) -> Result<(String, f64)> {
        Ok(match self.tr.as_mut() {
            Some(tr) => {
                let r = tr.step(&mut self.fl, self.dt)?;
                (format!("{r:?}"), r.last.max())
            }
            None => {
                let r = self.fl.iterate()?;
                (format!("{r:?}"), r.max())
            }
        })
    }

    /// The restart of the own part (local order, the identities of
    /// `mesh`): what a run of one rank writes, on a rank its part.
    fn restart(&mut self, mesh: &Mesh, name: &str) -> Result<crate::restart::Restart> {
        let fl = &mut self.fl;
        Ok(crate::restart::Restart {
            case: name.into(),
            cell_id: mesh.cell_id[..mesh.n_cells].to_vec(),
            face_id: mesh.face_id.clone(),
            state: fl.state()?,
            iteration: fl.iteration(),
            amg_setups: fl.amg_setups,
            amg_matrix: fl.amg_matrix.clone(),
            sweep: fl.sweep.clone(),
            transient: match &self.tr {
                Some(tr) => Some((tr.old_levels(fl)?, tr.clock, self.dt, tr.scheme)),
                None => None,
            },
            geometry: None,
            cht: None,
        })
    }

    /// Continues from `r`, a restart in the local order of the part (as
    /// the steady and the transient run of one rank load it).
    fn load(&mut self, case: &FlowCase, r: &crate::restart::Restart) -> Result<()> {
        anyhow::ensure!(
            r.transient.is_some() == self.tr.is_some(),
            "the restart file is {}, the case {}",
            if r.transient.is_some() {
                "transient"
            } else {
                "steady"
            },
            if self.tr.is_some() {
                "transient"
            } else {
                "steady"
            }
        );
        self.fl.load(&r.state)?;
        if let (Some(tr), Some((old, clock, _, _))) = (self.tr.as_mut(), &r.transient) {
            tr.restore(&mut self.fl, old, *clock)?;
        }
        self.fl.iteration = r.iteration;
        if let Some(m) = &r.amg_matrix {
            self.fl.restore_amg(m, r.amg_setups)?;
        }
        if self.tr.is_none()
            && case.solver.algorithm(false) == crate::simplec::Algorithm::Coupled
            && case.solver.coupled.momentum_sweep == crate::coupled::MomentumSweep::Multigrid
        {
            self.fl.restore_sweep(&r.sweep)?;
        }
        Ok(())
    }

    /// The fields of the own cells (`n`, global identities `ids`).
    fn fields(&mut self, n: usize, ids: &[u64], history: Vec<String>) -> Result<FlowSolve> {
        let s = self.fl.state()?;
        let own = |v: &[f64]| v[..n].to_vec();
        let mut extra: Vec<(String, Vec<f64>)> = s
            .scalars
            .iter()
            .map(|sc| (sc.name.clone(), own(&sc.phi)))
            .collect();
        if let Some(t) = &s.turbulence {
            extra.push(("k".into(), own(&t.k)));
            extra.push(("omega".into(), own(&t.omega)));
            extra.push(("nu_t".into(), own(&t.nu_t)));
        }
        Ok(FlowSolve {
            cell_id: ids[..n].to_vec(),
            u: [own(&s.u.x), own(&s.u.y), own(&s.u.z)],
            p: own(&s.p),
            grad_p: [own(&s.grad_p.x), own(&s.grad_p.y), own(&s.grad_p.z)],
            psi: [own(&s.psi.x), own(&s.psi.y), own(&s.psi.z)],
            extra,
            history,
            seconds: 0.0,
        })
    }
}

/// The wall distance of every cell of `global` when `case` has a
/// turbulence model (the exact distance, D48), else `None`.
fn global_wall_distance(case: &FlowCase, global: &Mesh) -> Result<Option<Vec<f64>>> {
    if case.turbulence.is_laminar() {
        return Ok(None);
    }
    let geo = Geometry::compute(global);
    Ok(Some(crate::derived::wall_distance(case, global, &geo)?))
}

/// The start velocity: a smooth field that is not uniform, so that every
/// term of the first iterations is exercised.
fn start(geo: &Geometry, n: usize) -> foehncraft_fv::disc::Soa3 {
    let mut u = foehncraft_fv::disc::Soa3::zeros(n);
    for c in 0..n {
        let x = geo.cell_centroid[c];
        u.x[c] = 0.5 + 0.3 * (2.0 * x[1]).sin();
        u.y[c] = 0.1 * (3.0 * x[0]).cos();
        u.z[c] = 0.0;
    }
    u
}

fn run(sv: &mut Solver, geo: &Geometry, n: usize, ids: &[u64], steps: usize) -> Result<FlowSolve> {
    sv.start(&start(geo, n), &vec![0.0; n])?;
    let mut history = Vec::with_capacity(steps);
    let t = std::time::Instant::now();
    for _ in 0..steps {
        history.push(sv.advance()?.0);
    }
    let seconds = t.elapsed().as_secs_f64();
    let mut s = sv.fields(n, ids, history)?;
    s.seconds = seconds;
    Ok(s)
}

/// The solve on the undecomposed mesh (the one-rank reference): `steps`
/// outer iterations, or time steps of a transient case.
pub fn solve_one(
    case: &FlowCase,
    mesh: &Mesh,
    steps: usize,
    threads: Option<usize>,
) -> Result<FlowSolve> {
    let geo = Geometry::compute(mesh);
    let target = Target {
        threads,
        ..Target::default()
    };
    let mut sv = Solver::new(case, mesh, &geo, &target, global_wall_distance(case, mesh)?)?;
    run(&mut sv, &geo, mesh.n_cells, &mesh.cell_id, steps)
}

/// The solve of `comm`'s rank on its part of `global` (`part`: the rank of
/// every global cell, one halo layer).
pub fn solve_part(
    case: &FlowCase,
    global: &Mesh,
    part: &[u32],
    comm: SharedComm,
    steps: usize,
    threads: Option<usize>,
) -> Result<FlowSolve> {
    let local = extract(global, part, comm.rank(), 1).map_err(anyhow::Error::msg)?;
    local.mesh.validate().map_err(anyhow::Error::msg)?;
    let geo = local.geometry(&*comm);
    let n = local.mesh.n_cells;
    let d = global_wall_distance(case, global)?.map(|d| {
        local.cell_position[..n]
            .iter()
            .map(|&g| d[g as usize])
            .collect()
    });
    let target = Target {
        threads,
        ranks: Some(Ranks::new(comm, &local, part)),
        ..Target::default()
    };
    let mut sv = Solver::new(case, &local.mesh, &geo, &target, d)?;
    run(&mut sv, &geo, n, &local.mesh.cell_id, steps)
}

/// `mine` (a rank's part) against `reference` (one rank), cell by cell.
pub fn check(rank: usize, mine: &FlowSolve, reference: &FlowSolve) -> RankCheck {
    let at: std::collections::HashMap<u64, usize> = reference
        .cell_id
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let same = |a: f64, b: f64| a.to_bits() == b.to_bits();
    let (mut fields, mut grads) = (0, 0);
    for (i, g) in mine.cell_id.iter().enumerate() {
        let r = at[g];
        let f = (0..3).any(|k| !same(mine.u[k][i], reference.u[k][r]))
            || !same(mine.p[i], reference.p[r])
            || mine
                .extra
                .iter()
                .zip(&reference.extra)
                .any(|((na, a), (nb, b))| na != nb || !same(a[i], b[r]));
        let d = (0..3).any(|k| {
            !same(mine.grad_p[k][i], reference.grad_p[k][r])
                || !same(mine.psi[k][i], reference.psi[k][r])
        });
        fields += f as usize;
        grads += d as usize;
    }
    RankCheck {
        rank,
        cells: mine.cell_id.len(),
        differing_fields: fields,
        differing_gradients: grads,
        same_history: mine.history == reference.history
            && mine.extra.len() == reference.extra.len(),
    }
}

/// The cases of the checks, a channel with a velocity inlet, walls,
/// symmetry planes and a pressure outlet on polyhedra, mixed cells and
/// tetrahedra: steady SIMPLEC with linear upwind and the Venkatakrishnan
/// limiter and bounded central, CG with Jacobi and with multigrid for the
/// pressure; transient PISO with BDF2; a passive scalar; the energy
/// equation (transient); the k-ω SST and the k-ε model with scalable wall
/// functions; the coupled solver.
pub fn check_cases() -> Result<Vec<CheckCase>> {
    let with = |mesh: &str, solver: &str, fluid: &str, extra: &str| -> Result<FlowCase> {
        let text = format!(
            r#"
            case_format = "0.6"
            [mesh.generate]
            {mesh}
            [fluid]
            density = 1.0
            viscosity = 0.02
            {fluid}
            [boundary.xmin]
            type = "velocity-inlet"
            velocity = [1.0, 0.1, 0.0]
            [boundary.xmax]
            type = "static-pressure-outlet"
            pressure = 0.0
            [boundary.ymin]
            type = "wall"
            [boundary.ymax]
            type = "wall"
            velocity = [0.5, 0.0, 0.0]
            [boundary.zmin]
            type = "symmetry"
            [boundary.zmax]
            type = "wall"
            [solver]
            max_iterations = 100
            {solver}
            {extra}
            "#
        );
        FlowCase::from_toml(&text, std::path::Path::new("."))
    };
    let case = |mesh: &str, solver: &str| with(mesh, solver, "", "");
    let poly = r#"type = "poly"
                dimension = 3
                upper = [2.0, 1.0, 0.5]
                cells = [8, 5, 3]
                perturb = 0.15"#;
    let mixed = r#"type = "mixed"
                dimension = 3
                upper = [2.0, 1.0, 0.5]
                cells = [8, 6, 3]
                perturb = 0.15"#;
    Ok(vec![
        CheckCase {
            name: "channel, polyhedra, linear upwind, CG + multigrid",
            case: case(
                r#"type = "poly"
                dimension = 3
                upper = [2.0, 1.0, 0.5]
                cells = [8, 5, 3]
                perturb = 0.15"#,
                r#"[solver.momentum]
                limiter = "venkatakrishnan""#,
            )?,
            iterations: 12,
        },
        CheckCase {
            name: "channel, mixed cells, bounded central",
            case: case(
                r#"type = "mixed"
                dimension = 3
                upper = [2.0, 1.0, 0.5]
                cells = [8, 6, 3]
                perturb = 0.15"#,
                r#"[solver.momentum]
                scheme = "central"
                limiter = "venkatakrishnan-wang""#,
            )?,
            iterations: 10,
        },
        CheckCase {
            name: "channel, tetrahedra, Jacobi pressure solver",
            case: case(
                r#"type = "tet"
                dimension = 3
                upper = [2.0, 1.0, 0.5]
                cells = [6, 4, 3]
                perturb = 0.1"#,
                r#"[solver.pressure]
                preconditioner = "jacobi"
                max_iter = 400"#,
            )?,
            iterations: 10,
        },
        CheckCase {
            name: "channel, polyhedra, transient, BDF2",
            case: with(
                poly,
                r#"[solver.momentum]
                limiter = "venkatakrishnan""#,
                "",
                r#"[time]
                scheme = "bdf2"
                time_step = 0.05
                steps = 6"#,
            )?,
            iterations: 6,
        },
        CheckCase {
            name: "channel, polyhedra, passive scalar",
            case: with(
                poly,
                r#"[solver.momentum]
                limiter = "venkatakrishnan""#,
                "",
                r#"[scalars.c]
                diffusivity = 0.01
                initial = "0.5 + 0.4 * sin(3 * x) * y"
                [boundary.xmin.scalar.c]
                value = 1.0
                [boundary.ymin.scalar.c]
                value = 0.0"#,
            )?,
            iterations: 8,
        },
        CheckCase {
            name: "channel, mixed cells, energy, transient",
            case: with(
                mixed,
                "",
                r#"specific_heat = 1000.0
                conductivity = 20.0"#,
                r#"[energy]
                initial = 300.0
                [boundary.xmin.energy]
                temperature = 300.0
                [boundary.ymin.energy]
                temperature = 320.0
                [time]
                scheme = "bdf2"
                time_step = 0.05
                steps = 5"#,
            )?,
            iterations: 5,
        },
        CheckCase {
            name: "channel, polyhedra, SST",
            case: with(
                poly,
                r#"[solver.momentum]
                limiter = "venkatakrishnan""#,
                "",
                r#"[turbulence]
                model = "k-omega-sst"
                [initial]
                turbulence_intensity = 0.05
                viscosity_ratio = 10.0"#,
            )?,
            iterations: 8,
        },
        CheckCase {
            name: "channel, mixed cells, k-epsilon, scalable wall functions",
            case: with(
                mixed,
                "",
                "",
                r#"[turbulence]
                model = "k-epsilon"
                wall_treatment = "scalable"
                [initial]
                turbulence_intensity = 0.05
                viscosity_ratio = 10.0"#,
            )?,
            iterations: 8,
        },
        CheckCase {
            name: "channel, polyhedra, coupled solver",
            case: case(
                poly,
                r#"algorithm = "coupled"
                [solver.momentum]
                limiter = "venkatakrishnan""#,
            )?,
            iterations: 6,
        },
    ])
}

/// What `[run] ranks` found for one number of ranks.
#[derive(Debug, Clone)]
pub struct RanksRun {
    /// Number of ranks.
    pub ranks: usize,
    /// Outer iterations (the same stopping rule as one rank).
    pub iterations: usize,
    /// Cells whose velocity, pressure, pressure gradient or limiter differs
    /// in any bit from one rank, over all ranks.
    pub differing_cells: usize,
    /// Every rank's residual history is that of one rank.
    pub same_history: bool,
    /// The restart gathered from the ranks is byte for byte the restart
    /// of one rank (checked for the largest count, with HDF5).
    pub restart_bytes: Option<bool>,
}

impl RanksRun {
    /// Every rank is bitwise one rank.
    pub fn bitwise(&self) -> bool {
        self.differing_cells == 0 && self.same_history && self.restart_bytes != Some(false)
    }
}

/// Fails unless `[run] ranks` can solve `case`: steady or transient
/// SIMPLEC, PISO or the coupled solver, laminar or with the k-ω SST or a
/// k-ε model, with passive scalars and the energy equation of one fluid;
/// from the start fields or a restart file (W6); without gravity, turning
/// zones, a schedule, levels, regions or the model options that are not on
/// ranks yet.
pub fn supported(case: &FlowCase) -> Result<()> {
    use crate::turbulence::{QcrMode, TurbulenceModel as Model};
    let n = &case.run.ranks;
    anyhow::ensure!(
        n.iter().all(|&r| (1..=256).contains(&r)),
        "[run] ranks must lie in 1..=256"
    );
    let t = &case.turbulence;
    let missing: Vec<&str> = [
        (
            !matches!(
                t.model,
                Model::Laminar | Model::KOmegaSst | Model::KEpsilon | Model::RealizableKEpsilon
            ),
            "a turbulence model other than k-ω SST and k-ε",
        ),
        (
            t.curvature_correction || t.qcr != QcrMode::None,
            "the curvature correction or the QCR",
        ),
        (case.gravity.is_some(), "gravity"),
        (
            !case.mrf.is_empty() || !case.sliding.is_empty(),
            "turning zones",
        ),
        (!case.schedule.is_empty(), "a [[schedule]]"),
        (case.verification.levels.len() > 1, "a sweep over levels"),
        (
            !case.regions.is_empty(),
            "regions (conjugate heat transfer)",
        ),
    ]
    .into_iter()
    .filter_map(|(on, what)| on.then_some(what))
    .collect();
    anyhow::ensure!(
        missing.is_empty(),
        "[run] ranks does not solve this case yet (MPI steps W5, W6): it has {}",
        missing.join(", ")
    );
    Ok(())
}

/// The start fields of `case` on the whole mesh `mesh` (its start values,
/// D238): computed once, replicated, so that every rank starts from the
/// values of its cells in the one-rank run (a part's start pressure comes
/// from boundaries that may lie on other ranks).
fn start_fields(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
) -> Result<(foehncraft_fv::disc::Soa3, Vec<f64>)> {
    let p0 = case.initial.pressure;
    let n = mesh.n_cells;
    Ok(
        match crate::start::initial_flow(case, mesh, geo, case.initial.velocity, p0, false)? {
            Some((u, p, _)) => (u, p),
            None => {
                let v = case.initial.velocity;
                (
                    foehncraft_fv::disc::Soa3 {
                        x: vec![v[0]; n],
                        y: vec![v[1]; n],
                        z: vec![v[2]; n],
                    },
                    vec![p0; n],
                )
            }
        },
    )
}

/// Where a solve of [`solve_case`] begins.
enum Begin {
    /// The start fields `u`, `p` of the own cells (the potential flow
    /// first when the case asks for it).
    Fields((foehncraft_fv::disc::Soa3, Vec<f64>)),
    /// A restart in the local order of the part.
    Restart(Box<crate::restart::Restart>),
}

/// The case's flow from `begin`, iterated until the largest residual is
/// below the tolerance or `max_iterations` (or its time steps); the
/// solver is returned for its restart.
fn solve_case(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    target: &Target,
    begin: Begin,
    wall_distance: Option<Vec<f64>>,
) -> Result<(FlowSolve, Solver)> {
    let mut sv = Solver::new(case, mesh, geo, target, wall_distance)?;
    let n = mesh.n_cells;
    match begin {
        Begin::Fields((u, p)) => {
            sv.start(&u, &p)?;
            if case.initial.potential_flow {
                sv.fl.potential_flow(mesh, geo)?;
                if let Some(tr) = sv.tr.as_mut() {
                    tr.start(&mut sv.fl)?;
                }
            }
        }
        Begin::Restart(r) => sv.load(case, &r)?,
    }
    let mut history = Vec::new();
    match case.time {
        // a fixed number of time steps (D93 compares fixed steps; the
        // adaptive step is the same on every rank, its Courant number a
        // global maximum)
        Some(t) => {
            let steps = t.steps.unwrap_or_else(|| {
                t.end_time.map_or(1, |e| {
                    ((e - t.start_time) / t.time_step).round().max(1.0) as usize
                })
            });
            for _ in 0..steps {
                history.push(sv.advance()?.0);
            }
        }
        None => {
            for _ in 0..case.solver.max_iterations {
                let (line, r) = sv.advance()?;
                history.push(line);
                if r < case.solver.tolerance {
                    break;
                }
            }
        }
    }
    let s = sv.fields(n, &mesh.cell_id, history)?;
    Ok((s, sv))
}

/// `[run] ranks` of a case: the flow on one rank, then on every listed
/// number of ranks (threads of this process, recursive coordinate
/// bisection), each rank compared with one rank. A `[run] restart` file is
/// read on every rank count (each rank takes its part by global identity,
/// MPI step W6). The restart of the largest count is gathered into the
/// global order and compared byte for byte with the restart of one rank.
/// With HDF5, `ranks.vtkhdf` gets the solution of the largest count
/// gathered on the global mesh and the cell field `rank`.
pub fn verify(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    out: &std::path::Path,
) -> Result<Vec<RanksRun>> {
    use foehncraft_mesh::partition::{Method, partition};
    use foehncraft_parallel::ThreadComm;
    supported(case)?;
    let cores = case
        .run
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |c| c.get()));
    let start = start_fields(case, mesh, geo)?;
    let distance = global_wall_distance(case, mesh)?;
    // a restart file in the order of the global mesh
    let restart = match &case.run.restart {
        Some(path) => {
            let n = mesh.n_cells;
            let (r, _) = crate::restart::read(path)?.to_local(
                &mesh.cell_id[..n],
                &mesh.face_id,
                mesh.n_interior_faces(),
            )?;
            local_only(&r)?;
            Some(r)
        }
        None => None,
    };
    let begin = |pick: &dyn Fn(&[f64]) -> Vec<f64>| {
        Begin::Fields((
            foehncraft_fv::disc::Soa3 {
                x: pick(&start.0.x),
                y: pick(&start.0.y),
                z: pick(&start.0.z),
            },
            pick(&start.1),
        ))
    };
    let (one, mut one_solver) = solve_case(
        case,
        mesh,
        geo,
        &Target {
            threads: Some(cores),
            ..Target::default()
        },
        match &restart {
            Some(r) => Begin::Restart(Box::new(r.clone())),
            None => begin(&|v: &[f64]| v.to_vec()),
        },
        distance.clone(),
    )?;
    let one_restart = if crate::restart::AVAILABLE {
        Some(one_solver.restart(mesh, "ranks")?)
    } else {
        None
    };
    drop(one_solver);
    let mut runs = Vec::new();
    let largest = case.run.ranks.iter().copied().max().unwrap_or(1);
    for &n in &case.run.ranks {
        let part = partition(mesh, geo, n, Method::Rcb).map_err(anyhow::Error::msg)?;
        let per_rank = (cores / n).max(1);
        let gather = n == largest && one_restart.is_some();
        let solves = ThreadComm::run(
            n,
            |c| -> Result<(FlowSolve, Option<crate::restart::Restart>)> {
                let comm: SharedComm = c;
                let local = extract(mesh, &part, comm.rank(), 1).map_err(anyhow::Error::msg)?;
                local.mesh.validate().map_err(anyhow::Error::msg)?;
                let lg = local.geometry(&*comm);
                let layout = Layout::new(comm.clone(), mesh, &local)?;
                let target = Target {
                    threads: Some(per_rank),
                    ranks: Some(Ranks::new(comm, &local, &part)),
                    ..Target::default()
                };
                // the own cells of the replicated start fields
                let own = &local.cell_position[..local.mesh.n_cells];
                let pick = |v: &[f64]| own.iter().map(|&g| v[g as usize]).collect::<Vec<f64>>();
                let d = distance.as_ref().map(|d| pick(d));
                let b = match &restart {
                    Some(r) => Begin::Restart(Box::new(layout.scatter(r.clone())?)),
                    None => begin(&pick),
                };
                let (s, mut sv) = solve_case(case, &local.mesh, &lg, &target, b, d)?;
                let g = if gather {
                    let r = sv.restart(&local.mesh, "ranks")?;
                    Some(layout.gather(r)?)
                } else {
                    None
                };
                Ok((s, g))
            },
        );
        let mut run = RanksRun {
            ranks: n,
            iterations: one.history.len(),
            differing_cells: 0,
            same_history: true,
            restart_bytes: None,
        };
        let mut parts = Vec::new();
        let mut gathered = None;
        for (rank, s) in solves.into_iter().enumerate() {
            let (s, g) = s?;
            let c = check(rank, &s, &one);
            run.iterations = s.history.len();
            run.differing_cells += c.differing_fields.max(c.differing_gradients);
            run.same_history &= c.same_history;
            parts.push(s);
            if rank == 0 {
                gathered = g;
            }
        }
        if n == largest {
            write_ranks(out, mesh, &part, &parts)?;
        }
        if let (Some(g), Some(r1)) = (gathered, &one_restart) {
            let a = out.join("ranks.restart.h5");
            let b = out.join("ranks.one.restart.h5");
            crate::restart::write(&a, &g)?;
            crate::restart::write(&b, r1)?;
            let same = std::fs::read(&a)? == std::fs::read(&b)?;
            std::fs::remove_file(&b)?;
            run.restart_bytes = Some(same);
        }
        runs.push(run);
    }
    Ok(runs)
}

/// What the restart check of MPI step W6 found.
#[derive(Debug, Clone)]
pub struct RestartCheck {
    /// Ranks that wrote the restart file.
    pub written_on: usize,
    /// Ranks that continued from it.
    pub read_on: usize,
    /// The file gathered from the ranks is byte for byte the file of one
    /// rank at the same iteration or step.
    pub same_bytes: bool,
    /// Every rank of the continuation against the uninterrupted run of one
    /// rank.
    pub continuation: Vec<RankCheck>,
}

impl RestartCheck {
    /// The file and the continuation are bitwise those of one rank.
    pub fn bitwise(&self) -> bool {
        self.same_bytes && self.continuation.iter().all(RankCheck::bitwise)
    }
}

/// The restart on ranks (MPI step W6, D93): `first` outer iterations or
/// time steps of `check` on `written_on` ranks, the restart gathered into
/// one file in `dir` (`ranks.restart.h5`, compared byte for byte with
/// `one.restart.h5` of one rank at the same point), then `second` more on
/// `read_on` ranks from that file, compared cell by cell and residual by
/// residual with one rank that ran `first + second` without the
/// interruption. Ranks are threads of this process ([`ThreadComm`]).
///
/// [`ThreadComm`]: foehncraft_parallel::ThreadComm
pub fn restart_check(
    check_case: &CheckCase,
    first: usize,
    second: usize,
    written_on: usize,
    read_on: usize,
    dir: &std::path::Path,
) -> Result<RestartCheck> {
    use foehncraft_mesh::partition::{Method, partition};
    use foehncraft_parallel::ThreadComm;
    let case = &check_case.case;
    let mesh = check_case.mesh()?;
    let geo = Geometry::compute(&mesh);
    let start = (start(&geo, mesh.n_cells), vec![0.0; mesh.n_cells]);
    let distance = global_wall_distance(case, &mesh)?;
    std::fs::create_dir_all(dir)?;
    let one_file = dir.join("one.restart.h5");
    let ranks_file = dir.join("ranks.restart.h5");

    // one rank, uninterrupted; its restart after `first`
    let target = Target {
        threads: Some(1),
        ..Target::default()
    };
    let mut sv = Solver::new(case, &mesh, &geo, &target, distance.clone())?;
    sv.start(&start.0, &start.1)?;
    for _ in 0..first {
        sv.advance()?;
    }
    crate::restart::write(&one_file, &sv.restart(&mesh, "restart-check")?)?;
    let mut history = Vec::new();
    for _ in 0..second {
        history.push(sv.advance()?.0);
    }
    let reference = sv.fields(mesh.n_cells, &mesh.cell_id, history)?;
    drop(sv);

    let on_part = |comm: SharedComm,
                   part: &[u32]|
     -> Result<(Solver, foehncraft_mesh::decompose::LocalMesh, Layout)> {
        let local = extract(&mesh, part, comm.rank(), 1).map_err(anyhow::Error::msg)?;
        local.mesh.validate().map_err(anyhow::Error::msg)?;
        let lg = local.geometry(&*comm);
        let layout = Layout::new(comm.clone(), &mesh, &local)?;
        let n = local.mesh.n_cells;
        let d = distance.as_ref().map(|d| {
            local.cell_position[..n]
                .iter()
                .map(|&g| d[g as usize])
                .collect()
        });
        let target = Target {
            threads: Some(1),
            ranks: Some(Ranks::new(comm, &local, part)),
            ..Target::default()
        };
        let sv = Solver::new(case, &local.mesh, &lg, &target, d)?;
        Ok((sv, local, layout))
    };

    // `first` on `written_on` ranks, gathered, written by rank 0
    let part = partition(&mesh, &geo, written_on, Method::Rcb).map_err(anyhow::Error::msg)?;
    for r in ThreadComm::run(written_on, |c| -> Result<()> {
        let rank = foehncraft_parallel::Communicator::rank(&*c);
        let (mut sv, local, layout) = on_part(c, &part)?;
        let own = &local.cell_position[..local.mesh.n_cells];
        let pick = |v: &[f64]| own.iter().map(|&g| v[g as usize]).collect::<Vec<f64>>();
        let u = foehncraft_fv::disc::Soa3 {
            x: pick(&start.0.x),
            y: pick(&start.0.y),
            z: pick(&start.0.z),
        };
        sv.start(&u, &pick(&start.1))?;
        for _ in 0..first {
            sv.advance()?;
        }
        let g = layout.gather(sv.restart(&local.mesh, "restart-check")?)?;
        if rank == 0 {
            crate::restart::write(&ranks_file, &g)?;
        }
        Ok(())
    }) {
        r?;
    }
    let same_bytes = std::fs::read(&one_file)? == std::fs::read(&ranks_file)?;

    // `second` more on `read_on` ranks from the gathered file
    let part = partition(&mesh, &geo, read_on, Method::Rcb).map_err(anyhow::Error::msg)?;
    let continuation = ThreadComm::run(read_on, |c| -> Result<RankCheck> {
        let rank = foehncraft_parallel::Communicator::rank(&*c);
        let (mut sv, local, layout) = on_part(c, &part)?;
        let r = layout.scatter(crate::restart::read(&ranks_file)?)?;
        sv.load(case, &r)?;
        let mut history = Vec::new();
        for _ in 0..second {
            history.push(sv.advance()?.0);
        }
        let mine = sv.fields(local.mesh.n_cells, &local.mesh.cell_id, history)?;
        Ok(check(rank, &mine, &reference))
    })
    .into_iter()
    .collect::<Result<Vec<_>>>()?;
    Ok(RestartCheck {
        written_on,
        read_on,
        same_bytes,
        continuation,
    })
}

/// `[run] ranks` of a case after its run: [`verify`] on its mesh, one
/// line per rank count to `line`; returns the failures.
pub(crate) fn report(
    case: &FlowCase,
    out: &std::path::Path,
    line: &mut dyn FnMut(&str),
) -> Result<Vec<String>> {
    let (mesh, geo, _) = crate::case::mesh_for(case, None)?;
    let what = if case.time.is_some() {
        "time steps"
    } else {
        "outer iterations"
    };
    let mut failures = Vec::new();
    for r in verify(case, &mesh, &geo, out)? {
        let l = format!(
            "solve on {} ranks: {} {what}, {}",
            r.ranks,
            r.iterations,
            if r.differing_cells == 0 && r.same_history {
                "bitwise equal to one rank".to_string()
            } else {
                format!(
                    "{} cells and {} the history differ from one rank",
                    r.differing_cells,
                    if r.same_history { "not" } else { "also" }
                )
            }
        );
        line(&l);
        if !(r.differing_cells == 0 && r.same_history) {
            failures.push(l);
        }
        if let Some(same) = r.restart_bytes {
            let l = format!(
                "restart file gathered from {} ranks (ranks.restart.h5): {}",
                r.ranks,
                if same {
                    "byte for byte the file of one rank"
                } else {
                    "differs from the file of one rank"
                }
            );
            line(&l);
            if !same {
                failures.push(l);
            }
        }
    }
    Ok(failures)
}

/// `ranks.vtkhdf`: the solution of the ranks on the global mesh and the
/// rank of every cell.
#[cfg(feature = "hdf5")]
fn write_ranks(
    out: &std::path::Path,
    mesh: &Mesh,
    part: &[u32],
    parts: &[FlowSolve],
) -> Result<()> {
    use foehncraft_mesh::vtkhdf::{CellField, write_all};
    let at: std::collections::HashMap<u64, usize> = mesh.cell_id[..mesh.n_cells]
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let n = mesh.n_cells;
    let (mut u, mut p) = ([vec![0.0; n], vec![0.0; n], vec![0.0; n]], vec![0.0; n]);
    // the scalars, the temperature and the turbulence variables
    let mut extra: Vec<(String, Vec<f64>)> = parts
        .first()
        .map(|s| {
            s.extra
                .iter()
                .map(|(name, _)| (name.clone(), vec![0.0; n]))
                .collect()
        })
        .unwrap_or_default();
    for s in parts {
        for (i, g) in s.cell_id.iter().enumerate() {
            let c = at[g];
            for (uk, sk) in u.iter_mut().zip(&s.u) {
                uk[c] = sk[i];
            }
            p[c] = s.p[i];
            for ((_, e), (_, se)) in extra.iter_mut().zip(&s.extra) {
                e[c] = se[i];
            }
        }
    }
    let speed: Vec<f64> = (0..n)
        .map(|c| (u[0][c] * u[0][c] + u[1][c] * u[1][c] + u[2][c] * u[2][c]).sqrt())
        .collect();
    let rank: Vec<f64> = part.iter().map(|&r| r as f64).collect();
    let fields = [
        CellField {
            name: "velocity_x",
            values: &u[0],
        },
        CellField {
            name: "velocity_y",
            values: &u[1],
        },
        CellField {
            name: "velocity_z",
            values: &u[2],
        },
        CellField {
            name: "speed",
            values: &speed,
        },
        CellField {
            name: "pressure",
            values: &p,
        },
        CellField {
            name: "rank",
            values: &rank,
        },
    ];
    let mut fields = Vec::from(fields);
    fields.extend(
        extra
            .iter()
            .map(|(name, values)| CellField { name, values }),
    );
    write_all(out.join("ranks.vtkhdf"), mesh, &fields, &[], &[])
}

#[cfg(not(feature = "hdf5"))]
fn write_ranks(_: &std::path::Path, _: &Mesh, _: &[u32], _: &[FlowSolve]) -> Result<()> {
    Ok(())
}

/// A rank's place in the global order of the mesh (MPI step W6): the
/// gather of its part of a restart into the order of the one-rank run, and
/// the scatter of a restart in that order to its part. The gathered file is
/// byte for byte the one-rank file (D84), so a run written on some number
/// of ranks continues on any other number (by global identity).
pub struct Layout {
    comm: SharedComm,
    /// Identities of the global mesh (cells, faces).
    global_cell_id: Vec<u64>,
    global_face_id: Vec<u64>,
    global_interior: usize,
    /// Identities of the local mesh (own cells, faces).
    cell_id: Vec<u64>,
    face_id: Vec<u64>,
    local_interior: usize,
    /// Global position of every own and halo cell.
    position: Vec<u32>,
    /// Global position of every local face.
    faces: Vec<usize>,
    /// Local faces whose owner is an own cell: each global face is
    /// gathered from the rank of its owner.
    own_faces: Vec<usize>,
    /// Every rank's own cells, own-owner faces and boundary faces as
    /// global positions (gathered once).
    all_cells: Vec<Vec<usize>>,
    all_faces: Vec<Vec<usize>>,
    all_boundary: Vec<Vec<usize>>,
    local_sell: foehncraft_fv::disc::SellStructure,
    global_sell: foehncraft_fv::disc::SellStructure,
    /// Matrix entries per row (diagonal and interior faces), local and
    /// global.
    local_row: Vec<usize>,
    global_row: Vec<usize>,
}

fn row_lengths(mesh: &Mesh) -> Vec<usize> {
    let ni = mesh.n_interior_faces();
    (0..mesh.n_cells)
        .map(|c| 1 + mesh.cell_face_iter(c).filter(|&(f, _)| f < ni).count())
        .collect()
}

impl Layout {
    /// The layout of `comm`'s rank whose local mesh `local` is a part of
    /// `global`. Collective.
    pub fn new(
        comm: SharedComm,
        global: &Mesh,
        local: &foehncraft_mesh::decompose::LocalMesh,
    ) -> Result<Self> {
        let lm = &local.mesh;
        let n = lm.n_cells;
        let at: std::collections::HashMap<u64, usize> = global
            .face_id
            .iter()
            .enumerate()
            .map(|(i, &g)| (g, i))
            .collect();
        let faces: Vec<usize> = lm
            .face_id
            .iter()
            .map(|g| {
                at.get(g)
                    .copied()
                    .ok_or_else(|| anyhow::anyhow!("face {g} is not in the global mesh"))
            })
            .collect::<Result<_>>()?;
        let li = lm.n_interior_faces();
        let gi = global.n_interior_faces();
        let own_faces: Vec<usize> = (0..lm.n_faces())
            .filter(|&f| (lm.owner[f] as usize) < n)
            .collect();
        let as_f64 = |v: &mut dyn Iterator<Item = usize>| v.map(|x| x as f64).collect::<Vec<f64>>();
        let back = |v: Vec<Vec<f64>>| -> Vec<Vec<usize>> {
            v.into_iter()
                .map(|r| r.into_iter().map(|x| x as usize).collect())
                .collect()
        };
        let all_cells = back(comm.all_gather_vec_f64(&as_f64(
            &mut local.cell_position[..n].iter().map(|&g| g as usize),
        )));
        let all_faces =
            back(comm.all_gather_vec_f64(&as_f64(&mut own_faces.iter().map(|&f| faces[f]))));
        let all_boundary =
            back(comm.all_gather_vec_f64(&as_f64(&mut faces[li..].iter().map(|&g| g - gi))));
        Ok(Self {
            global_cell_id: global.cell_id[..global.n_cells].to_vec(),
            global_face_id: global.face_id.clone(),
            global_interior: gi,
            cell_id: lm.cell_id[..n].to_vec(),
            face_id: lm.face_id.clone(),
            local_interior: li,
            position: local.cell_position.clone(),
            faces,
            own_faces,
            all_cells,
            all_faces,
            all_boundary,
            local_sell: foehncraft_fv::disc::sell_structure(lm),
            global_sell: foehncraft_fv::disc::sell_structure(global),
            local_row: row_lengths(lm),
            global_row: row_lengths(global),
            comm,
        })
    }

    /// The restart of the whole run from this rank's part `r`, in the
    /// global order (the content of the one-rank file). Collective; every
    /// rank receives it.
    pub fn gather(&self, r: crate::restart::Restart) -> Result<crate::restart::Restart> {
        local_only(&r)?;
        r.map_fields(
            &mut Gather(self),
            &self.global_cell_id,
            &self.global_face_id,
        )
    }

    /// This rank's part of `r`, a restart of the whole mesh (any cell and
    /// face order: it is put into the global order by identity first).
    pub fn scatter(&self, r: crate::restart::Restart) -> Result<crate::restart::Restart> {
        local_only(&r)?;
        let (r, _) = r.to_local(
            &self.global_cell_id,
            &self.global_face_id,
            self.global_interior,
        )?;
        r.map_fields(&mut Scatter(self), &self.cell_id, &self.face_id)
    }

    /// Whether this rank and all others agree that `ok`.
    fn agree(&self, ok: bool) -> bool {
        self.comm.min(if ok { 1.0 } else { 0.0 }) > 0.0
    }

    /// The global value array from every rank's values at `positions`.
    fn assemble(&self, mine: &[f64], len: usize, positions: &[Vec<usize>]) -> Result<Vec<f64>> {
        let parts = self.comm.all_gather_vec_f64(mine);
        let mut out = vec![0.0; len];
        let mut filled = 0;
        for (vals, pos) in parts.iter().zip(positions) {
            anyhow::ensure!(vals.len() == pos.len(), "gather: the parts do not match");
            for (&v, &g) in vals.iter().zip(pos) {
                out[g] = v;
            }
            filled += pos.len();
        }
        anyhow::ensure!(filled == len, "gather: {filled} of {len} values");
        Ok(out)
    }
}

/// The parts of a restart that are not mapped on ranks (a conjugate
/// case's full mesh, the points of a moving mesh) are refused there.
fn local_only(r: &crate::restart::Restart) -> Result<()> {
    anyhow::ensure!(
        r.cht.is_none() && r.geometry.is_none(),
        "restart files of conjugate cases and moving meshes are not on ranks yet"
    );
    Ok(())
}

struct Gather<'a>(&'a Layout);

impl crate::restart::FieldMap for Gather<'_> {
    fn cells(&mut self, v: &[f64]) -> Result<Vec<f64>> {
        let l = self.0;
        let n = l.cell_id.len();
        let ok = l.agree(v.len() >= n);
        anyhow::ensure!(ok, "gather: a cell array is shorter than the part");
        l.assemble(&v[..n], l.global_cell_id.len(), &l.all_cells)
    }
    fn faces(&mut self, v: &[f64]) -> Result<Vec<f64>> {
        let l = self.0;
        let ok = l.agree(v.len() == l.face_id.len());
        anyhow::ensure!(ok, "gather: a face array has another length than the part");
        let mine: Vec<f64> = l.own_faces.iter().map(|&f| v[f]).collect();
        l.assemble(&mine, l.global_face_id.len(), &l.all_faces)
    }
    fn boundary(&mut self, v: &[f64]) -> Result<Vec<f64>> {
        let l = self.0;
        let ok = l.agree(v.len() == l.face_id.len() - l.local_interior);
        anyhow::ensure!(
            ok,
            "gather: a boundary array has another length than the part"
        );
        l.assemble(
            v,
            l.global_face_id.len() - l.global_interior,
            &l.all_boundary,
        )
    }
    fn matrix(&mut self, v: &[f64]) -> Result<Option<Vec<f64>>> {
        let l = self.0;
        // the entries of the own rows in slot order; the slots of a row
        // are the same on ranks and on one rank (diagonal, then the
        // interior faces in global face order, W1)
        let mut mine = Vec::new();
        let mut ok = v.len() == l.local_sell.stored();
        if ok {
            'rows: for (c, &len) in l.local_row.iter().enumerate() {
                let g = l.position[c] as usize;
                ok &= l.global_row[g] == len;
                for j in 0..len {
                    let k = l.local_sell.pos(c, j);
                    let col = l.position[l.local_sell.col[k] as usize];
                    ok &= l.global_sell.col[l.global_sell.pos(g, j)] == col;
                    if !ok {
                        break 'rows;
                    }
                    mine.push(v[k]);
                }
            }
        }
        anyhow::ensure!(
            l.agree(ok),
            "gather: a matrix is not in the storage of the part"
        );
        let parts = l.comm.all_gather_vec_f64(&mine);
        let mut out = vec![0.0; l.global_sell.stored()];
        for (vals, cells) in parts.iter().zip(&l.all_cells) {
            let mut k = 0;
            for &g in cells {
                for j in 0..l.global_row[g] {
                    out[l.global_sell.pos(g, j)] = vals[k];
                    k += 1;
                }
            }
        }
        Ok(Some(out))
    }
}

struct Scatter<'a>(&'a Layout);

impl crate::restart::FieldMap for Scatter<'_> {
    fn cells(&mut self, v: &[f64]) -> Result<Vec<f64>> {
        let l = self.0;
        anyhow::ensure!(v.len() == l.global_cell_id.len(), "scatter: another mesh");
        Ok(l.position[..l.cell_id.len()]
            .iter()
            .map(|&g| v[g as usize])
            .collect())
    }
    fn faces(&mut self, v: &[f64]) -> Result<Vec<f64>> {
        let l = self.0;
        anyhow::ensure!(v.len() == l.global_face_id.len(), "scatter: another mesh");
        Ok(l.faces.iter().map(|&g| v[g]).collect())
    }
    fn boundary(&mut self, v: &[f64]) -> Result<Vec<f64>> {
        let l = self.0;
        let gi = l.global_interior;
        anyhow::ensure!(
            v.len() == l.global_face_id.len() - gi,
            "scatter: another mesh"
        );
        Ok(l.faces[l.local_interior..]
            .iter()
            .map(|&g| v[g - gi])
            .collect())
    }
    fn matrix(&mut self, v: &[f64]) -> Result<Option<Vec<f64>>> {
        let l = self.0;
        anyhow::ensure!(
            v.len() == l.global_sell.stored(),
            "scatter: a matrix of another mesh"
        );
        let mut out = vec![0.0; l.local_sell.stored()];
        for (c, &len) in l.local_row.iter().enumerate() {
            let g = l.position[c] as usize;
            for j in 0..len {
                let k = l.local_sell.pos(c, j);
                let gk = l.global_sell.pos(g, j);
                anyhow::ensure!(
                    l.global_sell.col[gk] == l.position[l.local_sell.col[k] as usize],
                    "scatter: a matrix row of another layout"
                );
                out[k] = v[gk];
            }
        }
        Ok(Some(out))
    }
}
