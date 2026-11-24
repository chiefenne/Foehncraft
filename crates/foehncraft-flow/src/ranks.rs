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
    for _ in 0..steps {
        history.push(sv.advance()?.0);
    }
    sv.fields(n, ids, history)
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
}

impl RanksRun {
    /// Every rank is bitwise one rank.
    pub fn bitwise(&self) -> bool {
        self.differing_cells == 0 && self.same_history
    }
}

/// Fails unless `[run] ranks` can solve `case`: steady or transient
/// SIMPLEC, PISO or the coupled solver, laminar or with the k-ω SST or a
/// k-ε model, with passive scalars and the energy equation of one fluid;
/// without gravity, turning zones, a schedule, levels, a restart, regions
/// or the model options that are not on ranks yet (MPI step W5 so far).
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
        (case.run.restart.is_some(), "a restart"),
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
        "[run] ranks does not solve this case yet (MPI step W5): it has {}",
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

/// The case's flow from the start fields `u`, `p` of the own cells (the
/// potential flow first when the case asks for it), iterated until the
/// largest residual is below the tolerance or `max_iterations`.
fn solve_case(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    target: &Target,
    (u, p): (foehncraft_fv::disc::Soa3, Vec<f64>),
    wall_distance: Option<Vec<f64>>,
) -> Result<FlowSolve> {
    let mut sv = Solver::new(case, mesh, geo, target, wall_distance)?;
    let n = mesh.n_cells;
    sv.start(&u, &p)?;
    if case.initial.potential_flow {
        sv.fl.potential_flow(mesh, geo)?;
        if let Some(tr) = sv.tr.as_mut() {
            tr.start(&mut sv.fl)?;
        }
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
    sv.fields(n, &mesh.cell_id, history)
}

/// `[run] ranks` of a case: the flow on one rank, then on every listed
/// number of ranks (threads of this process, recursive coordinate
/// bisection), each rank compared with one rank. With HDF5,
/// `ranks.vtkhdf` gets the solution of the largest count gathered on the
/// global mesh and the cell field `rank`.
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
    let one = solve_case(
        case,
        mesh,
        geo,
        &Target {
            threads: Some(cores),
            ..Target::default()
        },
        start.clone(),
        distance.clone(),
    )?;
    let mut runs = Vec::new();
    let largest = case.run.ranks.iter().copied().max().unwrap_or(1);
    for &n in &case.run.ranks {
        let part = partition(mesh, geo, n, Method::Rcb).map_err(anyhow::Error::msg)?;
        let per_rank = (cores / n).max(1);
        let solves = ThreadComm::run(n, |c| -> Result<FlowSolve> {
            let comm: SharedComm = c;
            let local = extract(mesh, &part, comm.rank(), 1).map_err(anyhow::Error::msg)?;
            local.mesh.validate().map_err(anyhow::Error::msg)?;
            let lg = local.geometry(&*comm);
            let target = Target {
                threads: Some(per_rank),
                ranks: Some(Ranks::new(comm, &local, &part)),
                ..Target::default()
            };
            // the own cells of the replicated start fields
            let own = &local.cell_position[..local.mesh.n_cells];
            let pick = |v: &[f64]| own.iter().map(|&g| v[g as usize]).collect::<Vec<f64>>();
            let u = foehncraft_fv::disc::Soa3 {
                x: pick(&start.0.x),
                y: pick(&start.0.y),
                z: pick(&start.0.z),
            };
            let d = distance.as_ref().map(|d| pick(d));
            solve_case(case, &local.mesh, &lg, &target, (u, pick(&start.1)), d)
        });
        let mut run = RanksRun {
            ranks: n,
            iterations: one.history.len(),
            differing_cells: 0,
            same_history: true,
        };
        let mut parts = Vec::new();
        for (rank, s) in solves.into_iter().enumerate() {
            let s = s?;
            let c = check(rank, &s, &one);
            run.iterations = s.history.len();
            run.differing_cells += c.differing_fields.max(c.differing_gradients);
            run.same_history &= c.same_history;
            parts.push(s);
        }
        if n == largest {
            write_ranks(out, mesh, &part, &parts)?;
        }
        runs.push(run);
    }
    Ok(runs)
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
            if r.bitwise() {
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
        if !r.bitwise() {
            failures.push(l);
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
