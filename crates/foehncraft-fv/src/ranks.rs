//! Scalar solves on the ranks of a decomposed mesh, and their check against
//! one rank (MPI step W3, D93).
//!
//! [`solve_part`] runs the whole steady scalar solver of [`crate::solver`]
//! (gradient, limiter, face coefficients, assembly, BiCGStab or CG with the
//! Jacobi preconditioner) on one rank's local mesh with
//! [`CpuBackend::on_ranks`]; [`solve_one`] runs it on the undecomposed mesh.
//! [`check`] compares the two bit for bit, keyed by global identities: the
//! solution, its gradient and limiter, every matrix entry and right-hand
//! side, the residual history and the iteration counts. The tests run it
//! with ranks as threads ([`foehncraft_parallel::ThreadComm`]), the example
//! `mpi_scalar` across MPI processes in CI.

use crate::amg::ReductionSettings;
use crate::backend::cpu::CpuBackend;
use crate::backend::{Backend, V};
use crate::disc::Discretisation;
use crate::problem::Problem;
use crate::solver::{self, Outcome, Settings};
use anyhow::Result;
use foehncraft_kernels::unpack_face;
use foehncraft_mesh::decompose::{Ranks, SharedComm, extract};
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// The result of a scalar solve, keyed by global identities.
#[derive(Debug, Clone, Default)]
pub struct ScalarSolve {
    /// Global identity of each own cell.
    pub cell_id: Vec<u64>,
    /// φ per own cell.
    pub phi: Vec<f64>,
    /// Gradient per own cell.
    pub grad: [Vec<f64>; 3],
    /// Limiter value per own cell.
    pub psi: Vec<f64>,
    /// Right-hand side per own cell (of the last assembly).
    pub rhs: Vec<f64>,
    /// Matrix entries (row, column) of the own rows by global identities,
    /// diagonal included, from the last assembly.
    pub matrix: HashMap<(u64, u64), f64>,
    /// The solver's outcome (residual history, iterations).
    pub outcome: Outcome,
    /// φ after two cycles of the reduction multigrid (D179) on the last
    /// matrix, when the case asks for it (empty otherwise).
    pub reduction_phi: Vec<f64>,
}

/// What [`check`] found on one rank.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct RankCheck {
    /// The rank.
    pub rank: usize,
    /// Own cells.
    pub cells: usize,
    /// Own cells whose φ differs in any bit from the one-rank run.
    pub differing_phi: usize,
    /// Own cells whose gradient or limiter value differs in any bit.
    pub differing_gradients: usize,
    /// Matrix entries and right-hand sides of the own rows that differ in
    /// any bit (or are missing).
    pub differing_matrix: usize,
    /// Whether the residual of every outer iteration and the number of
    /// linear iterations are bitwise those of the one-rank run.
    pub same_history: bool,
}

impl RankCheck {
    /// True when nothing differs.
    pub fn bitwise_equal(&self) -> bool {
        self.differing_phi == 0
            && self.differing_gradients == 0
            && self.differing_matrix == 0
            && self.same_history
    }
}

fn collect(
    be: &mut dyn Backend,
    d: &Discretisation,
    ids: &[u64],
    outcome: Outcome,
) -> Result<ScalarSolve> {
    let n = d.n_cells;
    // the gradient and limiter of the final φ, as the next outer iteration
    // would see them
    let phi = be.download(V::Phi)?;
    let grad = [
        be.download(V::Gx)?,
        be.download(V::Gy)?,
        be.download(V::Gz)?,
    ];
    let psi = be.download(V::Psi)?;
    let rhs = be.download(V::Rhs)?;
    let val = be.matrix_values()?;
    let mut matrix = HashMap::new();
    for c in 0..n {
        matrix.insert((ids[c], ids[c]), val[d.sell.pos(c, 0)]);
        let mut j = 1;
        for &e in &d.cell_faces[d.cell_ptr[c] as usize..d.cell_ptr[c + 1] as usize] {
            let (f, is_nb) = unpack_face(e);
            if f < d.n_interior {
                let other = if is_nb { d.owner[f] } else { d.neighbour[f] } as usize;
                matrix.insert((ids[c], ids[other]), val[d.sell.pos(c, j)]);
                j += 1;
            }
        }
    }
    Ok(ScalarSolve {
        cell_id: ids[..n].to_vec(),
        phi,
        grad,
        psi,
        rhs,
        matrix,
        outcome,
        reduction_phi: Vec::new(),
    })
}

fn run(
    mut be: Box<dyn Backend>,
    d: &Discretisation,
    ids: &[u64],
    phi0: f64,
    settings: &Settings,
    reduction: Option<&ReductionSettings>,
) -> Result<ScalarSolve> {
    be.upload(V::Phi, &vec![phi0; d.n_cells])?;
    let outcome = solver::solve(be.as_mut(), settings, |_| {})?;
    let mut out = collect(be.as_mut(), d, ids, outcome)?;
    if let Some(rs) = reduction {
        // the reduction multigrid of the coupled momentum sweep on the last
        // matrix (not symmetric), as a stationary iteration from φ = 0
        solver::setup_reduction(be.as_mut(), rs)?;
        solver::cycles(be.as_mut(), 2)?;
        out.reduction_phi = be.download(V::Phi)?;
    }
    Ok(out)
}

/// The start value of φ: the offset of the manufactured solution, or 0.
fn start(problem: &Problem) -> f64 {
    problem.manufactured.map_or(0.0, |m| m.offset)
}

/// The solve on the undecomposed mesh (the one-rank reference).
pub fn solve_one(
    mesh: &Mesh,
    problem: &Problem,
    settings: &Settings,
    reduction: Option<&ReductionSettings>,
    threads: Option<usize>,
) -> Result<ScalarSolve> {
    let geo = Geometry::compute(mesh);
    let d = Discretisation::new(mesh, &geo, problem)?;
    let be = Box::new(CpuBackend::new(d.clone(), threads)?);
    run(be, &d, &mesh.cell_id, start(problem), settings, reduction)
}

/// The solve of this rank on its part of `global`: the local mesh with
/// `layers` halo layers (1 suffices for the scalar operators), the geometry
/// with the halo from the owners, the discretisation of the local mesh, and
/// the solver on [`CpuBackend::on_ranks`]. Every rank calls it with the
/// same arguments.
#[allow(clippy::too_many_arguments)]
pub fn solve_part(
    global: &Mesh,
    part: &[u32],
    layers: usize,
    comm: SharedComm,
    problem: &Problem,
    settings: &Settings,
    reduction: Option<&ReductionSettings>,
    threads: Option<usize>,
) -> Result<ScalarSolve> {
    let local = extract(global, part, comm.rank(), layers).map_err(anyhow::Error::msg)?;
    local.mesh.validate().map_err(anyhow::Error::msg)?;
    let geo = local.geometry(&*comm);
    let d = Discretisation::new(&local.mesh, &geo, problem)?;
    let ranks = Ranks::new(comm, &local, part);
    let be = Box::new(CpuBackend::on_ranks(d.clone(), threads, ranks)?);
    run(
        be,
        &d,
        &local.mesh.cell_id,
        start(problem),
        settings,
        reduction,
    )
}

/// Compares this rank's solve with the one-rank reference, bit for bit.
pub fn check(rank: usize, reference: &ScalarSolve, mine: &ScalarSolve) -> RankCheck {
    let at: HashMap<u64, usize> = reference
        .cell_id
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let same = |a: f64, b: f64| a.to_bits() == b.to_bits();
    let mut out = RankCheck {
        rank,
        cells: mine.cell_id.len(),
        ..Default::default()
    };
    for (i, id) in mine.cell_id.iter().enumerate() {
        let Some(&g) = at.get(id) else {
            out.differing_phi += 1;
            continue;
        };
        if !same(mine.phi[i], reference.phi[g]) {
            out.differing_phi += 1;
        }
        if (0..3).any(|k| !same(mine.grad[k][i], reference.grad[k][g]))
            || !same(mine.psi[i], reference.psi[g])
        {
            out.differing_gradients += 1;
        }
        if !same(mine.rhs[i], reference.rhs[g]) {
            out.differing_matrix += 1;
        }
        if !mine.reduction_phi.is_empty()
            && !same(mine.reduction_phi[i], reference.reduction_phi[g])
        {
            out.differing_phi += 1;
        }
    }
    for (key, &v) in &mine.matrix {
        if reference.matrix.get(key).is_none_or(|&r| !same(r, v)) {
            out.differing_matrix += 1;
        }
    }
    let (a, b) = (&mine.outcome, &reference.outcome);
    out.same_history = a.linear_iterations == b.linear_iterations
        && a.residuals.len() == b.residuals.len()
        && a.residuals
            .iter()
            .zip(&b.residuals)
            .all(|(x, y)| same(*x, *y));
    out
}

/// One case of the multi-rank check: a generated mesh, a problem and solver
/// settings.
#[derive(Debug, Clone)]
pub struct CheckCase {
    /// Name.
    pub name: &'static str,
    /// The global mesh.
    pub mesh: Mesh,
    /// The problem.
    pub problem: Problem,
    /// Solver settings.
    pub settings: Settings,
    /// Also two cycles of the reduction multigrid on the last matrix.
    pub reduction: Option<ReductionSettings>,
}

/// The cases of the multi-rank check (tests and the example `mpi_scalar`):
/// every scheme, limiter and Krylov method of the scalar solver on
/// perturbed meshes of several cell types, few outer iterations each.
///
/// - convection of an oblique step (linear upwind, Venkatakrishnan limiter,
///   non-orthogonal correction, BiCGStab) on polyhedra;
/// - the same with the bounded central scheme and the Venkatakrishnan–Wang
///   limiter (a global extremum) on a mixed 3D mesh;
/// - diffusion of a manufactured solution (CG) on tetrahedra;
/// - the same with the smoothed-aggregation multigrid (CG + AMG, mixed
///   precision, several levels) on polyhedra (MPI step W4);
/// - the oblique step with BiCGStab + AMG on hexahedra, then two cycles of
///   the reduction multigrid of D179 on its last matrix (W4).
pub fn check_cases() -> Result<Vec<CheckCase>> {
    use foehncraft_mesh::generate::{BoxSpec, Kind};
    let mesh = |kind: Kind, dim: u8, cells: [usize; 3]| -> Result<Mesh> {
        let mut spec = BoxSpec::new(kind, dim, cells);
        spec.perturb = 0.2;
        spec.generate().map_err(anyhow::Error::msg)
    };
    let step = |extra: &str| -> Result<Problem> {
        let mut p: Problem = toml::from_str(&format!(
            r#"
            diffusivity = 1e-3
            velocity = [0.8, 0.5, 0.3]
            {extra}
            [boundary.xmin]
            type = "fixed-value"
            profile = "step"
            axis = 1
            at = 0.4
            below = 1.0
            above = 0.0
            [boundary.ymin]
            type = "fixed-value"
            value = 1.0
            [boundary.zmin]
            type = "fixed-gradient"
            value = 0.5
            [boundary.xmax]
            type = "zero-gradient"
            [boundary.ymax]
            type = "zero-gradient"
            [boundary.zmax]
            type = "zero-gradient"
            "#
        ))?;
        p.link();
        Ok(p)
    };
    let mut mms: Problem = toml::from_str(
        r#"
        diffusivity = 0.2
        velocity = [0.0, 0.0, 0.0]
        manufactured = {}
        [boundary.xmin]
        type = "fixed-value"
        profile = "manufactured"
        [boundary.xmax]
        type = "fixed-gradient"
        profile = "manufactured"
        [boundary.ymin]
        type = "fixed-value"
        profile = "manufactured"
        [boundary.ymax]
        type = "fixed-value"
        profile = "manufactured"
        [boundary.zmin]
        type = "fixed-value"
        profile = "manufactured"
        [boundary.zmax]
        type = "fixed-value"
        profile = "manufactured"
        "#,
    )?;
    mms.link();
    let lu = Settings {
        limiter: crate::solver::LimiterKind::Venkatakrishnan,
        limiter_k: 1.0,
        max_outer: 12,
        linear_max_iter: 40,
        ..Settings::default()
    };
    let central = Settings {
        scheme: crate::solver::Scheme::Central,
        central_bounding: 0.5,
        limiter: crate::solver::LimiterKind::VenkatakrishnanWang,
        limiter_k: 0.5,
        relaxation: 0.9,
        ..lu.clone()
    };
    let cg = Settings {
        linear_solver: crate::solver::LinearSolver::Cg,
        max_outer: 6,
        linear_max_iter: 60,
        ..Settings::default()
    };
    let amg = crate::amg::AmgSettings {
        max_coarse: 24,
        ..Default::default()
    };
    let cg_amg = Settings {
        preconditioner: crate::solver::Preconditioner::Amg,
        amg,
        ..cg.clone()
    };
    let lu_short = Settings {
        max_outer: 6,
        ..lu.clone()
    };
    Ok(vec![
        CheckCase {
            name: "oblique step, polyhedra, linear upwind",
            mesh: mesh(Kind::Poly, 3, [6, 5, 4])?,
            problem: step("")?,
            settings: lu,
            reduction: None,
        },
        CheckCase {
            name: "oblique step, mixed, bounded central",
            mesh: mesh(Kind::Mixed, 3, [7, 6, 5])?,
            problem: step("")?,
            settings: central,
            reduction: None,
        },
        CheckCase {
            name: "manufactured diffusion, tetrahedra, CG",
            mesh: mesh(Kind::Tet, 3, [5, 5, 4])?,
            problem: mms.clone(),
            settings: cg,
            reduction: None,
        },
        CheckCase {
            name: "manufactured diffusion, hexahedra, CG + AMG",
            mesh: mesh(Kind::Hex, 3, [10, 9, 8])?,
            problem: mms,
            settings: cg_amg,
            reduction: None,
        },
        CheckCase {
            name: "oblique step, hexahedra, reduction multigrid",
            mesh: mesh(Kind::Hex, 3, [14, 12, 10])?,
            problem: step("")?,
            settings: lu_short,
            reduction: Some(ReductionSettings {
                max_coarse: 24,
                ..Default::default()
            }),
        },
    ])
}
