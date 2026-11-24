//! The steady scalar solver: outer iterations of deferred correction, each a
//! gradient, a limiter, the face coefficients, the assembly and a linear
//! solve: BiCGStab (D36) or, for symmetric matrices such as the pressure
//! equation, CG (D35), preconditioned by Jacobi or by the algebraic
//! multigrid of [`crate::amg`].
//!
//! The solver drives a [`Backend`] and only moves scalars (dot products,
//! norms) to the host, so the same code runs on the CPU and on the GPU.

use crate::amg::{self, AmgSettings, HierarchyInfo, ReductionSettings};
use crate::backend::{Backend, FaceSettings, Limiter, V};
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Convection scheme (D25).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Scheme {
    /// First-order upwind.
    Upwind,
    /// Second-order linear upwind (with the limiter of `limiter`).
    #[default]
    LinearUpwind,
    /// Bounded central: the linear interpolation on interior faces with
    /// the normalised-variable bounding of `central_bounding`, blended with
    /// linear upwind by `central_blending` (D25).
    Central,
}

/// Limiter of the case file.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum LimiterKind {
    /// No limiter.
    #[default]
    None,
    /// Venkatakrishnan.
    Venkatakrishnan,
    /// Venkatakrishnan–Wang.
    VenkatakrishnanWang,
    /// Barth–Jespersen: ψ = min(1, Δ₊/Δ₋), strictly inside the
    /// neighbourhood bounds (not differentiable).
    BarthJespersen,
    /// The multi-dimensional limiting process (Park & Kim 2012): the
    /// reconstruction at the cell vertices bounded by the cells that share
    /// each vertex; scalars and the energy equation only.
    Mlp,
}

/// Krylov method of the linear solves.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "kebab-case")]
pub enum LinearSolver {
    /// BiCGStab (van der Vorst 1992), any matrix (D36).
    #[default]
    Bicgstab,
    /// Flexible preconditioned conjugate gradients (Notay 2000), symmetric
    /// positive definite matrices (D35).
    Cg,
}

/// Preconditioner of the linear solves.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Preconditioner {
    /// Diagonal (Jacobi).
    #[default]
    Jacobi,
    /// One V-cycle of the smoothed-aggregation multigrid (stage 3).
    Amg,
}

/// `[solver]` settings.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Settings {
    /// Convection scheme.
    pub scheme: Scheme,
    /// `scheme = "central"`: fraction of the central scheme, the rest
    /// linear upwind (0 to 1).
    pub central_blending: f64,
    /// `scheme = "central"`: β of the normalised-variable bounding, the
    /// normalised upwind value below which the scheme moves to upwind
    /// (0: plain central, unbounded).
    pub central_bounding: f64,
    /// Gradient limiter.
    pub limiter: LimiterKind,
    /// Constant K of the limiter (Venkatakrishnan: ε² = (K h)³;
    /// Venkatakrishnan–Wang: ε = K (φ_max − φ_min)).
    pub limiter_k: f64,
    /// Non-orthogonal diffusion correction (D26).
    pub nonorth_correction: bool,
    /// Limit of the correction as a multiple of the orthogonal part
    /// (negative: no limit).
    pub nonorth_limit: f64,
    /// Implicit under-relaxation factor α (0 < α ≤ 1).
    pub relaxation: f64,
    /// Largest number of outer iterations.
    pub max_outer: usize,
    /// Outer tolerance on the normalised residual ‖b − Aφ‖ / ‖b‖.
    pub tolerance: f64,
    /// Largest number of BiCGStab iterations per outer iteration.
    pub linear_max_iter: usize,
    /// Relative tolerance of each linear solve.
    pub linear_tolerance: f64,
    /// Krylov method.
    pub linear_solver: LinearSolver,
    /// Preconditioner.
    pub preconditioner: Preconditioner,
    /// Multigrid settings (`preconditioner = "amg"`).
    pub amg: AmgSettings,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            scheme: Scheme::LinearUpwind,
            central_blending: 1.0,
            central_bounding: 0.1,
            limiter: LimiterKind::None,
            limiter_k: 5.0,
            nonorth_correction: true,
            nonorth_limit: -1.0,
            relaxation: 1.0,
            max_outer: 200,
            tolerance: 1e-10,
            linear_max_iter: 200,
            linear_tolerance: 1e-3,
            linear_solver: LinearSolver::Bicgstab,
            preconditioner: Preconditioner::Jacobi,
            amg: AmgSettings::default(),
        }
    }
}

impl Settings {
    /// Checks the ranges.
    pub fn validate(&self) -> anyhow::Result<()> {
        anyhow::ensure!(
            self.relaxation > 0.0 && self.relaxation <= 1.0,
            "[solver] relaxation must lie in (0, 1]"
        );
        anyhow::ensure!(self.limiter_k >= 0.0, "[solver] limiter_k must be >= 0");
        anyhow::ensure!(
            self.limiter != LimiterKind::Mlp,
            "[solver] limiter = \"mlp\" is for the scalars and the energy equation of the flow solver"
        );
        anyhow::ensure!(
            (0.0..=1.0).contains(&self.central_blending),
            "[solver] central_blending must lie in [0, 1]"
        );
        anyhow::ensure!(
            (0.0..1.0).contains(&self.central_bounding),
            "[solver] central_bounding must lie in [0, 1)"
        );
        anyhow::ensure!(self.tolerance > 0.0, "[solver] tolerance must be > 0");
        anyhow::ensure!(
            self.linear_tolerance > 0.0 && self.linear_tolerance < 1.0,
            "[solver] linear_tolerance must lie in (0, 1)"
        );
        self.amg.validate()?;
        Ok(())
    }

    fn limiter(&self) -> Limiter {
        match self.limiter {
            LimiterKind::None => Limiter::None,
            LimiterKind::Venkatakrishnan => Limiter::Venkatakrishnan { k: self.limiter_k },
            LimiterKind::VenkatakrishnanWang => Limiter::VenkatakrishnanWang { k: self.limiter_k },
            LimiterKind::BarthJespersen => Limiter::BarthJespersen,
            // refused by `validate`
            LimiterKind::Mlp => Limiter::BarthJespersen,
        }
    }
}

/// Wall-clock seconds per phase (R6), measured with device synchronisation.
#[derive(Debug, Clone, Default, Serialize)]
pub struct Timings {
    /// Gradient.
    pub gradient: f64,
    /// Limiter.
    pub limiter: f64,
    /// Face coefficients.
    pub faces: f64,
    /// Assembly.
    pub assemble: f64,
    /// Linear solves.
    pub linear: f64,
    /// Multigrid setup (host), including the transfer of the matrix.
    pub amg_setup: f64,
}

impl Timings {
    /// Total seconds.
    pub fn total(&self) -> f64 {
        self.gradient + self.limiter + self.faces + self.assemble + self.linear + self.amg_setup
    }
}

/// Result of a solve.
#[derive(Debug, Clone, Default, Serialize)]
pub struct Outcome {
    /// Outer iterations done.
    pub outer: usize,
    /// Normalised residual per outer iteration (before its linear solve).
    pub residuals: Vec<f64>,
    /// BiCGStab iterations per outer iteration.
    pub linear_iterations: Vec<usize>,
    /// Converged to the outer tolerance.
    pub converged: bool,
    /// Phase timings.
    pub timings: Timings,
    /// Multigrid hierarchy, when used.
    pub amg: Option<HierarchyInfo>,
}

impl Outcome {
    /// Figure of merit (R6): nanoseconds per cell and outer iteration.
    pub fn ns_per_cell_iteration(&self, cells: usize) -> f64 {
        1e9 * self.timings.total() / (cells as f64 * self.outer.max(1) as f64)
    }
}

fn timed(
    be: &mut dyn Backend,
    acc: &mut f64,
    f: impl FnOnce(&mut dyn Backend) -> anyhow::Result<()>,
) -> anyhow::Result<()> {
    let t = Instant::now();
    f(be)?;
    be.sync()?;
    *acc += t.elapsed().as_secs_f64();
    Ok(())
}

/// Gradient (and limiter) of the current `Phi`.
pub fn gradient(be: &mut dyn Backend, s: &Settings, tm: &mut Timings) -> anyhow::Result<()> {
    timed(be, &mut tm.gradient, |b| b.gradient())?;
    let lim = if s.scheme == Scheme::Upwind {
        Limiter::None
    } else {
        s.limiter()
    };
    timed(be, &mut tm.limiter, |b| b.limiter(lim))
}

/// Solves the steady problem from the `Phi` held by the backend.
pub fn solve(
    be: &mut dyn Backend,
    s: &Settings,
    mut log: impl FnMut(&str),
) -> anyhow::Result<Outcome> {
    s.validate()?;
    let mut out = Outcome::default();
    let fs = FaceSettings {
        blend: if s.scheme == Scheme::Upwind { 0.0 } else { 1.0 },
        central: if s.scheme == Scheme::Central {
            s.central_blending
        } else {
            0.0
        },
        beta: s.central_bounding,
        nonorth: s.nonorth_correction,
        nonorth_limit: s.nonorth_limit,
    };
    for it in 0..s.max_outer {
        let mut tm = std::mem::take(&mut out.timings);
        gradient(be, s, &mut tm)?;
        timed(be, &mut tm.faces, |b| b.faces(fs))?;
        timed(be, &mut tm.assemble, |b| b.assemble(s.relaxation))?;
        // residual of the unrelaxed system at the current φ (the relaxation
        // terms cancel there)
        let t = Instant::now();
        be.residual(V::Phi, V::Rhs, V::R)?;
        let rr = be.dot(V::R, V::R)?.sqrt();
        let bb = be.dot(V::Rhs, V::Rhs)?.sqrt();
        let res = rr / bb.max(f64::MIN_POSITIVE);
        out.residuals.push(res);
        out.outer = it + 1;
        if res < s.tolerance {
            out.converged = true;
            out.linear_iterations.push(0);
            be.sync()?;
            tm.linear += t.elapsed().as_secs_f64();
            out.timings = tm;
            log(&format!(
                "outer {:4}  residual {res:.3e}  converged",
                it + 1
            ));
            break;
        }
        if s.preconditioner == Preconditioner::Amg && out.amg.is_none() {
            // the matrix of the stage-2/3 problems does not change between
            // outer iterations (constant mass flux and diffusivity): one
            // setup serves all of them
            let ts = Instant::now();
            out.amg = Some(setup_amg(be, &s.amg)?);
            tm.amg_setup += ts.elapsed().as_secs_f64();
        }
        let lin = match s.linear_solver {
            LinearSolver::Bicgstab => bicgstab(
                be,
                s.preconditioner,
                s.linear_max_iter,
                s.linear_tolerance,
                rr,
            )?,
            LinearSolver::Cg => {
                cg(
                    be,
                    s.preconditioner,
                    s.linear_max_iter,
                    s.linear_tolerance,
                    rr,
                )?
                .0
            }
        };
        be.sync()?;
        tm.linear += t.elapsed().as_secs_f64();
        out.timings = tm;
        out.linear_iterations.push(lin);
        if it % 10 == 0 || it + 1 == s.max_outer {
            log(&format!(
                "outer {:4}  residual {res:.3e}  bicgstab {lin}",
                it + 1
            ));
        }
        anyhow::ensure!(res.is_finite(), "the solution diverged (residual {res})");
    }
    Ok(out)
}

/// Builds the multigrid hierarchy of the matrix assembled in the backend
/// and loads it there.
pub fn setup_amg(be: &mut dyn Backend, s: &AmgSettings) -> anyhow::Result<HierarchyInfo> {
    setup_amg_with(be, s, false)
}

/// As [`setup_amg`]; with `symmetric_part` the hierarchy is built from
/// ½(A + Aᵀ), the preconditioner of BiCGStab for a transport equation with
/// convection (D36; stage 6 PR 6d).
pub fn setup_amg_with(
    be: &mut dyn Backend,
    s: &AmgSettings,
    symmetric_part: bool,
) -> anyhow::Result<HierarchyInfo> {
    let t = Instant::now();
    let a = be.amg_matrix()?;
    let a = if symmetric_part {
        a.symmetric_part()
    } else {
        a
    };
    let mut h = amg::build(a, s)?;
    be.amg_load(&h)?;
    be.sync()?;
    h.info.setup_seconds = t.elapsed().as_secs_f64();
    Ok(h.info)
}

/// Builds the reduction multigrid (D179, [`amg::reduction`]) of the matrix
/// assembled in the backend, which need not be symmetric, and loads it into
/// the active slot.
pub fn setup_reduction(
    be: &mut dyn Backend,
    s: &ReductionSettings,
) -> anyhow::Result<HierarchyInfo> {
    let t = Instant::now();
    let a = be.amg_matrix()?;
    let mut h = amg::reduction::build(a, s)?;
    be.amg_load(&h)?;
    be.sync()?;
    h.info.setup_seconds = t.elapsed().as_secs_f64();
    Ok(h.info)
}

/// Phi = the result of `cycles` multigrid cycles of the active hierarchy on
/// A Phi = Rhs from Phi = 0: a stationary iteration, a fixed linear
/// operator on Rhs (the first cycle needs no matrix product outside the
/// hierarchy; the later ones take the residual with the matrix of the
/// backend).
pub fn cycles(be: &mut dyn Backend, cycles: usize) -> anyhow::Result<()> {
    amg::cycle::vcycle(be, V::Rhs, V::Phi)?;
    for _ in 1..cycles {
        be.residual(V::Phi, V::Rhs, V::R)?;
        amg::cycle::vcycle(be, V::R, V::Z)?;
        be.axpy(1.0, V::Z, V::Phi)?;
    }
    Ok(())
}

/// out = M⁻¹ x.
pub fn precondition(be: &mut dyn Backend, pc: Preconditioner, x: V, out: V) -> anyhow::Result<()> {
    match pc {
        Preconditioner::Jacobi => be.precondition(x, out),
        Preconditioner::Amg => amg::cycle::vcycle(be, x, out),
    }
}

/// Flexible preconditioned CG (Notay 2000: the Polak–Ribière form of β,
/// robust to a preconditioner that is symmetric only up to the rounding of
/// its float32 hierarchy) on A Phi = Rhs; expects R = Rhs − A Phi with norm
/// `r0`. Stops when ‖r‖ ≤ `rel_tol` r0. Returns the iterations done and the
/// final ‖r‖ / r0 of the recursive residual.
pub fn cg(
    be: &mut dyn Backend,
    pc: Preconditioner,
    max_iter: usize,
    rel_tol: f64,
    r0: f64,
) -> anyhow::Result<(usize, f64)> {
    let target = rel_tol * r0;
    if r0 == 0.0 {
        return Ok((0, 0.0));
    }
    precondition(be, pc, V::R, V::Z)?;
    be.copy(V::Z, V::P)?;
    let mut rho = be.dot(V::R, V::Z)?;
    let mut rn = r0;
    for it in 1..=max_iter {
        if rho == 0.0 {
            // z = M⁻¹ r vanishes although r does not: the residual lies below
            // the range of the float32 multigrid hierarchy (an exactly
            // converged one-dimensional flow, about 1e-40 per cell), so no
            // further step can change Phi. As in BiCGStab, stop.
            return Ok((it - 1, rn / r0));
        }
        be.spmv(V::P, V::Q)?;
        let pq = be.dot(V::P, V::Q)?;
        anyhow::ensure!(
            pq > 0.0 && rho > 0.0,
            "CG breakdown in iteration {it} ((p, A p) = {pq:.3e}, (r, z) = {rho:.3e}): \
             the matrix or the preconditioner is not positive definite"
        );
        let alpha = rho / pq;
        be.axpy(alpha, V::P, V::Phi)?;
        be.copy(V::R, V::Rold)?;
        be.axpy(-alpha, V::Q, V::R)?;
        rn = be.dot(V::R, V::R)?.sqrt();
        if rn <= target {
            return Ok((it, rn / r0));
        }
        precondition(be, pc, V::R, V::Z)?;
        let rho_new = be.dot(V::R, V::Z)?;
        let zr_old = be.dot(V::Z, V::Rold)?;
        let beta = (rho_new - zr_old) / rho;
        rho = rho_new;
        be.xpby(V::Z, beta, V::P)?;
    }
    Ok((max_iter, rn / r0))
}

/// BiCGStab (van der Vorst 1992) with right preconditioning on
/// A Phi = Rhs; expects R = Rhs − A Phi with norm `r0`. Returns the
/// iterations done.
pub fn bicgstab(
    be: &mut dyn Backend,
    pc: Preconditioner,
    max_iter: usize,
    rel_tol: f64,
    r0: f64,
) -> anyhow::Result<usize> {
    let target = rel_tol * r0;
    be.copy(V::R, V::Rh)?;
    be.copy(V::R, V::P)?;
    let mut rho = be.dot(V::Rh, V::R)?;
    for it in 1..=max_iter {
        if rho == 0.0 {
            return Ok(it - 1);
        }
        precondition(be, pc, V::P, V::Ph)?;
        be.spmv(V::Ph, V::Vv)?;
        let rv = be.dot(V::Rh, V::Vv)?;
        if rv == 0.0 {
            return Ok(it - 1);
        }
        let alpha = rho / rv;
        be.sub_scaled(V::R, V::Vv, alpha, V::S)?;
        let ss = be.dot(V::S, V::S)?.sqrt();
        if ss <= target {
            // Phi += α p̂
            be.copy(V::S, V::Sh)?;
            be.update_x(alpha, 0.0)?;
            return Ok(it);
        }
        precondition(be, pc, V::S, V::Sh)?;
        be.spmv(V::Sh, V::T)?;
        let ts = be.dot(V::T, V::S)?;
        let tt = be.dot(V::T, V::T)?;
        let omega = if tt > 0.0 { ts / tt } else { 0.0 };
        be.update_x(alpha, omega)?;
        be.sub_scaled(V::S, V::T, omega, V::R)?;
        let rn = be.dot(V::R, V::R)?.sqrt();
        if rn <= target || omega == 0.0 {
            return Ok(it);
        }
        let rho_new = be.dot(V::Rh, V::R)?;
        let beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;
        be.update_p(beta, omega)?;
    }
    Ok(max_iter)
}
