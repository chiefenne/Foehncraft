//! The incompressible flow solver state and the steady algorithm: SIMPLEC
//! (Van Doormaal & Raithby 1984), or SIMPLE (Patankar & Spalding 1972), on
//! the collocated, cell-centred grid with the Rhie–Chow face mass flux (D19,
//! D27). The transient PISO algorithm on the same state is in
//! [`crate::transient`].
//!
//! One outer iteration:
//!
//! 1. boundary values of velocity and pressure (extrapolated, symmetry,
//!    from the face flux or the total pressure) from the current fields;
//! 2. least-squares gradients of u, v, w and p (D24), limiters (D25);
//! 3. the relaxation term of the Rhie–Chow flux from the current fields,
//!    (1 − α_u)(ṁ_f − ρ ū_f·S);
//! 4. momentum: face coefficients, one matrix for the three components
//!    (relaxed with α_u), BiCGStab with Jacobi (D36) per component;
//! 5. predicted mass flux ṁ*_f by Rhie–Chow (relaxation independent: the
//!    converged result does not depend on α_u, Majumdar 1988);
//! 6. pressure correction: Σ c_f (p'_P − p'_N) = −Σ ±ṁ*_f with
//!    c_f = ρ D̄^C_f |S|²/(d·S), D^C = V / (a_P/α_u − Σ|a_N|) (SIMPLEC) or
//!    V / (a_P/α_u) (SIMPLE); CG with the smoothed-aggregation multigrid of
//!    stage 3 (D35), the hierarchy rebuilt every `amg_rebuild` iterations;
//! 7. corrections: ṁ_f −= c_f (p'_N − p'_P) (so the face fluxes satisfy
//!    continuity to the tolerance of the pressure solve), u −= D^C ∇p',
//!    p += α_p p'.
//!
//! The algorithm is a sequence of operations of a [`FlowDevice`]: the CPU
//! device runs the single-source bodies of `foehncraft_kernels::flow` in
//! threaded loops, the CUDA device the same bodies on the GPU. Each output
//! is written by one index, so the result does not depend on the number of
//! threads (D31).
//!
//! Residuals (dimensionless, per outer iteration, before the solves):
//! momentum ‖b − A u_i‖₂ / ‖a_P u‖₂ per component with ‖a_P u‖ over all
//! three components, continuity ‖Σ_f ±ṁ*_f‖₂ / ‖Σ_f |ṁ_f|‖₂ (net outflow
//! of the predicted fluxes over the cell throughflow).

use crate::bc::BoundaryArrays;
use crate::device::{
    self, BoundaryArgs, Buf, Buf3, Eps2, FlowDevice, FluxArgs, Lsq, Mask, MomentumFaceArgs,
    MomentumRowArgs, Setup, TimeTerms,
};
use anyhow::Result;
use foehncraft_fv::amg::AmgSettings;
use foehncraft_fv::disc::{Discretisation, Soa3};
use foehncraft_fv::problem::{Boundary, BoundaryKind, Problem};
use foehncraft_fv::solver::{LimiterKind, Preconditioner, Scheme};
use foehncraft_kernels::flow as k;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::time::Instant;

pub use crate::device::cpu::SERIAL_BELOW;

/// Pressure–velocity algorithm.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Algorithm {
    /// SIMPLEC (Van Doormaal & Raithby 1984), steady.
    #[default]
    Simplec,
    /// SIMPLE (Patankar & Spalding 1972), steady.
    Simple,
    /// PISO (Issa 1986), transient: per time step a momentum predictor and
    /// `correctors` pressure corrections, optionally repeated in
    /// `outer_iterations` (a PIMPLE-type loop).
    Piso,
    /// Coupled (L11), steady: the momentum and continuity equations of an
    /// outer iteration are solved as one linear system for the increments
    /// of the velocity and the pressure ([`crate::coupled`]).
    Coupled,
}

/// Gradient of the pressure (D24, D50), used in the momentum source
/// −∇p V, the Rhie–Chow flux, the extrapolated boundary pressures and the
/// velocity correction.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum PressureGradient {
    /// Weighted least squares (D24), the stencil points of walls and
    /// symmetry planes projected onto the face normal (zero normal
    /// gradient there).
    #[default]
    LeastSquares,
    /// Gauss: the face sum (1/V) Σ_f p_f S_f, p_f = w p_P + (1 − w) p_N on
    /// interior faces and the boundary pressure on boundary faces. Robust on
    /// tetrahedra with skewed faces, where least squares can make the steady
    /// iteration diverge; exact for a linear pressure only where faces are
    /// not skewed.
    Gauss,
}

impl PressureGradient {
    /// Name used in case files.
    pub fn name(self) -> &'static str {
        match self {
            Self::LeastSquares => "least-squares",
            Self::Gauss => "gauss",
        }
    }
}

/// Face velocity of the Rhie–Chow mass flux (D50).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum FaceVelocity {
    /// Linear interpolation to the point of the line P–N nearest the face
    /// plane, moved to the face centroid with the interpolated cell
    /// gradient: w u_P + (1 − w) u_N + ∇ū_f·(x_f − x_w).
    #[default]
    Interpolated,
    /// The mean of the two limited extrapolations to the face centroid,
    /// w (u_P + ψ_P ∇u_P·r_P) + (1 − w)(u_N + ψ_N ∇u_N·r_N), with the
    /// limiter ψ of the velocity (Venkatakrishnan when none is set). Exact
    /// for a linear velocity like the interpolated form, but the face value
    /// stays between the values the two cells see: on a skewed face
    /// (‖x_f − x_w‖ as long as the centroid distance) the shift with the
    /// unlimited gradient feeds a velocity mode back into continuity with a
    /// gain above one.
    Bounded,
}

impl FaceVelocity {
    /// Name used in case files.
    pub fn name(self) -> &'static str {
        match self {
            Self::Interpolated => "interpolated",
            Self::Bounded => "bounded",
        }
    }
}

/// The transposed part ∇·[μ_eff (∇u)ᵀ] of the viscous stress in turbulent
/// flow (stage 6 specification 2.4), explicit (D50).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum TransposedStress {
    /// −Σ_f μ_f (∇u)ᵀ_f·S_f with the interpolated cell gradients.
    #[default]
    Full,
    /// Without the trace: −Σ_f μ_f ((∇u)ᵀ_f − (∇·u)_f I)·S_f. The same term
    /// where ∇·u = 0; discretely it drops the grad-div part μ ∇(∇·u),
    /// which the cell gradients of tetrahedra do not make zero and which,
    /// explicit and as large as the implicit Laplacian where the eddy
    /// viscosity is large, amplifies a velocity mode and with it k.
    TraceFree,
}

impl TransposedStress {
    /// Name used in case files.
    pub fn name(self) -> &'static str {
        match self {
            Self::Full => "full",
            Self::TraceFree => "trace-free",
        }
    }
}

/// Cell-quality remediation (D50).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Remediation {
    /// On a mesh with faces whose skewness exceeds `skewness_limit`: the
    /// Gauss pressure gradient, for a second-order convection scheme
    /// without a limiter the Venkatakrishnan limiter, and the bounded face
    /// velocity of the mass flux; reported in the log,
    /// `flow.md` and `result.json`. A mesh without such faces runs the
    /// schemes as set.
    #[default]
    Auto,
    /// Run the schemes as set on every mesh.
    Off,
}

/// The pressure corrector of PISO (`solver.corrector`).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum PisoCorrector {
    /// `simplec` with outer iterations (`outer_iterations` > 1), `issa` for
    /// plain PISO.
    #[default]
    Auto,
    /// Issa's corrector: D = V/(a_P/α) and the neighbour correction
    /// H(u')/a_P between two correctors. It damps more in plain PISO
    /// (Crank–Nicolson at Courant number 4 settles with it), but where a_P
    /// is almost all neighbour coupling (thin wedge cells about an axis at
    /// high viscosity) the pressure can grow from one outer iteration to the
    /// next.
    Issa,
    /// The SIMPLEC form: D^C = V/(a_P/α − Σ|a_N|) and the neighbour
    /// correction Σ_N |a_N| (u'_N − u'_P)/a_P; stable on such cells with
    /// outer iterations.
    Simplec,
}

/// The cell-quality remediation of a run (D50): what the mesh check found
/// and what was changed.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct RemediationReport {
    /// The skewness limit of the check.
    pub skewness_limit: f64,
    /// Interior faces whose skewness exceeds the limit.
    pub skewed_faces: usize,
    /// Cells with such a face.
    pub poor_cells: usize,
    /// Largest face skewness of the mesh.
    pub max_skewness: f64,
    /// Changes of the schemes, `key: old -> new`; empty with
    /// `remediation = "off"` or when the case already uses them.
    pub changes: Vec<String>,
}

/// Face skewness (D50, as the mesh quality report): the distance between
/// the face centroid and the point where the line through the two cell
/// centroids crosses the face plane, over that distance; per interior
/// face.
pub fn face_skewness(d: &Discretisation) -> Vec<f64> {
    let c = &d.cell_centroid;
    let fc = &d.face_centroid;
    let a = &d.face_area;
    (0..d.n_interior)
        .map(|f| {
            let p = d.owner[f] as usize;
            let xp = [c.x[p], c.y[p], c.z[p]];
            let xn = d.neighbour_centroid(f);
            let dd = [xn[0] - xp[0], xn[1] - xp[1], xn[2] - xp[2]];
            let r = [fc.x[f] - xp[0], fc.y[f] - xp[1], fc.z[f] - xp[2]];
            let s = [a.x[f], a.y[f], a.z[f]];
            let t = (r[0] * s[0] + r[1] * s[1] + r[2] * s[2])
                / (dd[0] * s[0] + dd[1] * s[1] + dd[2] * s[2]);
            let e = [t * dd[0] - r[0], t * dd[1] - r[1], t * dd[2] - r[2]];
            ((e[0] * e[0] + e[1] * e[1] + e[2] * e[2])
                / (dd[0] * dd[0] + dd[1] * dd[1] + dd[2] * dd[2]))
                .sqrt()
        })
        .collect()
}

/// The weights of the Gauss gradient in the form of the least-squares
/// weights (D24), so the same gradient body computes it:
/// ∇φ_c = Σ_k w_k (φ_k − φ_c) with w_k = (1 − w_f) S_f / V_c for the owner
/// of an interior face, w_f (−S_f) / V_c for its neighbour and S_f / V_c on
/// a boundary face (φ_k the boundary value). With Σ_f S_f = 0 this is
/// (1/V) Σ_f φ_f S_f, φ_f = w_f φ_P + (1 − w_f) φ_N.
pub fn gauss_weights(d: &Discretisation) -> Soa3 {
    let ni = d.n_interior;
    let mut out = Soa3::zeros(d.cell_faces.len());
    for c in 0..d.n_cells {
        let v = d.cell_volume[c];
        for k in d.cell_ptr[c] as usize..d.cell_ptr[c + 1] as usize {
            let (f, is_nb) = foehncraft_kernels::unpack_face(d.cell_faces[k]);
            let g = if f >= ni {
                1.0 / v
            } else if is_nb {
                -d.weight[f] / v
            } else {
                (1.0 - d.weight[f]) / v
            };
            out.set(
                k,
                [
                    g * d.face_area.x[f],
                    g * d.face_area.y[f],
                    g * d.face_area.z[f],
                ],
            );
        }
    }
    out
}

/// Default α_u of the coupled algorithm: a local Courant number of at most
/// α/(1 − α) = 49. (A larger cap saves outer iterations, but the momentum
/// sweeps of the preconditioner, BiCGStab with a diagonal preconditioner,
/// then need more than they save.)
pub const COUPLED_RELAXATION: f64 = 0.98;

/// The same default for the k-ε family when its transport equations take
/// the reduction multigrid (D180): a local Courant number of at most 199.
/// On the backward-facing step the run settles in 65 instead of 90
/// iterations on 32 800 cells (with the turbulence relaxation 0.995 in
/// both); the laminar cavity does not gain (449 iterations at 0.98, 493
/// at 0.995), so it keeps 0.98.
pub const COUPLED_RELAXATION_REDUCTION: f64 = 0.995;

/// `[solver]` settings of a flow case.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Settings {
    /// Pressure–velocity algorithm (default SIMPLEC, PISO with `[time]`).
    pub algorithm: Option<Algorithm>,
    /// The momentum equations (`[solver.momentum]`): convection scheme,
    /// relaxation and linear solver.
    pub momentum: MomentumSettings,
    /// The pressure (correction) equation (`[solver.pressure]`).
    pub pressure: PressureSettings,
    /// The linear solve of the coupled algorithm (`[solver.coupled]`).
    pub coupled: crate::coupled::CoupledSettings,
    /// The pseudo-time steps of a steady run (`[solver.pseudo_time]`).
    pub pseudo_time: crate::coupled::PseudoTimeSettings,
    /// The passive scalars (`[solver.scalar]`, D51).
    pub scalar: crate::scalar::ScalarSettings,
    /// The energy equation (`[solver.energy]`, stage 7a1): the settings of
    /// a scalar.
    pub energy: crate::scalar::ScalarSettings,
    /// The turbulence equations (`[solver.turbulence]`, stage 6).
    pub turbulence: crate::turbulence::TurbulenceSettings,
    /// The elliptic equations of the turbulence models (`[solver.elliptic]`:
    /// f of k-ε-ζ-f, stage 6 PR 6d).
    pub elliptic: crate::turbulence::EllipticSettings,
    /// Gradient of the pressure (D24, D50): `least-squares`, or `gauss`,
    /// the face sum of the face pressures (robust on tetrahedra with skewed
    /// faces).
    pub pressure_gradient: PressureGradient,
    /// Face velocity of the Rhie–Chow mass flux (D50): `interpolated`, or
    /// `bounded`, the mean of the two limited extrapolations to the face
    /// centroid (robust on skewed faces).
    pub face_velocity: FaceVelocity,
    /// Transposed part of the viscous stress in turbulent flow (D50):
    /// `full`, or `trace-free` (the same term where ∇·u = 0, without the
    /// grad-div part that cell gradients on tetrahedra leave in it; a
    /// further safety net, not switched by the remediation).
    pub transposed_stress: TransposedStress,
    /// Cell-quality remediation (D50): `auto` switches to the Gauss pressure
    /// gradient, the Venkatakrishnan limiter and the bounded face velocity
    /// on a mesh with faces whose skewness exceeds `skewness_limit`, and
    /// reports it; `off` keeps the schemes as set.
    pub remediation: Remediation,
    /// Faces with a larger skewness (distance of the face centroid from the
    /// line through the two cell centroids, over the centroid distance) are
    /// poor (D50).
    pub skewness_limit: f64,
    /// The operator check (D236): a cell whose velocity gradient (relative
    /// to h |H|) or velocity Hessian (relative to |H|) has a larger error
    /// for a quadratic field is flagged.
    pub operator_tolerance: f64,
    /// Non-orthogonal correction of the viscous term (D26).
    pub nonorth_correction: bool,
    /// Largest number of outer iterations (steady).
    pub max_iterations: usize,
    /// Convergence: every residual below this value (steady).
    pub tolerance: f64,
    /// Rebuild the multigrid hierarchy every N outer iterations (steady) or
    /// N time steps (transient).
    pub amg_rebuild: usize,
    /// Log a line every N outer iterations or time steps.
    pub log_every: usize,
    /// PISO: pressure corrections per outer iteration.
    pub correctors: usize,
    /// PISO: outer iterations per time step (1: plain PISO).
    pub outer_iterations: usize,
    /// PISO: the pressure corrector (`auto`: `simplec` with outer
    /// iterations, `issa` for plain PISO).
    pub corrector: PisoCorrector,
    /// PISO: an outer iteration ends the time step when the momentum and
    /// continuity residuals of its start are below this value.
    pub outer_tolerance: f64,
    /// PISO: in the first outer iteration of a second-order step the
    /// explicit quantities are extrapolated from the two old levels to the
    /// new time: the convecting mass flux and the velocity gradients and
    /// limiters of the deferred corrections. Without it the lag of these
    /// terms makes plain PISO (one outer iteration) first order in time.
    pub extrapolation: bool,
    /// Set from the whole case, not read: the k-ε family runs the coupled
    /// algorithm with the reduction multigrid in its transport equations,
    /// so the default α_u is [`COUPLED_RELAXATION_REDUCTION`] (D180).
    #[serde(skip)]
    #[schemars(skip)]
    pub reduction_defaults: bool,
}

/// `[solver.momentum]`: the momentum equations (the velocity column of the
/// "solution flags" of the log).
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct MomentumSettings {
    /// Convection scheme of the momentum equations (during a run only
    /// between `upwind` and `linear-upwind`: a first-order start, then
    /// second order, with `[[schedule]]`).
    pub scheme: Scheme,
    /// `scheme = "central"`: fraction of the bounded central scheme on the
    /// interior faces, the rest linear upwind (0 to 1).
    pub central_blending: f64,
    /// `scheme = "central"`: β of the normalised-variable bounding; the
    /// scheme moves to upwind where the normalised upwind value leaves
    /// [β, 1) (0: plain central, unbounded).
    pub central_bounding: f64,
    /// Limiter of the velocity gradients.
    pub limiter: LimiterKind,
    /// Constant K of the limiter.
    pub limiter_k: f64,
    /// Under-relaxation of the momentum equations α_u (default 0.7 steady,
    /// 1 transient; with the coupled algorithm a cap on the local step:
    /// 0.98, or 0.995 for the k-ε family with the reduction multigrid in
    /// its transport equations, D180).
    pub relaxation: Option<f64>,
    /// Relative tolerance of each momentum solve (default 1e-2 steady, 1e-8
    /// transient).
    pub tolerance: Option<f64>,
    /// Largest number of BiCGStab iterations per momentum solve.
    pub max_iter: usize,
}

impl Default for MomentumSettings {
    fn default() -> Self {
        Self {
            scheme: Scheme::LinearUpwind,
            central_blending: 1.0,
            central_bounding: 0.1,
            limiter: LimiterKind::None,
            limiter_k: 5.0,
            relaxation: None,
            tolerance: None,
            max_iter: 50,
        }
    }
}

/// `[solver.pressure]`: the pressure (correction) equation.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct PressureSettings {
    /// Relaxation of the pressure correction α_p (default 1 for SIMPLEC and
    /// PISO, 0.3 for SIMPLE).
    pub relaxation: Option<f64>,
    /// Relative tolerance of each pressure-correction solve (default 1e-2
    /// steady, 1e-8 transient).
    pub tolerance: Option<f64>,
    /// Largest number of CG iterations per pressure solve.
    pub max_iter: usize,
    /// Preconditioner of the pressure solve.
    pub preconditioner: Preconditioner,
    /// Multigrid settings.
    pub amg: AmgSettings,
}

impl Default for PressureSettings {
    fn default() -> Self {
        Self {
            relaxation: None,
            tolerance: None,
            max_iter: 100,
            preconditioner: Preconditioner::Amg,
            amg: AmgSettings::default(),
        }
    }
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            algorithm: None,
            momentum: MomentumSettings::default(),
            pressure: PressureSettings::default(),
            coupled: crate::coupled::CoupledSettings::default(),
            pseudo_time: crate::coupled::PseudoTimeSettings::default(),
            scalar: crate::scalar::ScalarSettings::default(),
            energy: crate::scalar::ScalarSettings::default(),
            turbulence: crate::turbulence::TurbulenceSettings::default(),
            elliptic: crate::turbulence::EllipticSettings::default(),
            pressure_gradient: PressureGradient::LeastSquares,
            face_velocity: FaceVelocity::Interpolated,
            transposed_stress: TransposedStress::Full,
            remediation: Remediation::Auto,
            skewness_limit: 0.9,
            operator_tolerance: 0.01,
            nonorth_correction: true,
            max_iterations: 2000,
            tolerance: 1e-8,
            amg_rebuild: 10,
            log_every: 10,
            correctors: 2,
            corrector: PisoCorrector::Auto,
            outer_iterations: 1,
            outer_tolerance: 1e-6,
            extrapolation: true,
            reduction_defaults: false,
        }
    }
}

impl Settings {
    /// The algorithm in force: the one given, else PISO in a transient case
    /// and SIMPLEC in a steady one.
    pub fn algorithm(&self, transient: bool) -> Algorithm {
        self.algorithm.unwrap_or(if transient {
            Algorithm::Piso
        } else {
            Algorithm::Simplec
        })
    }

    /// Relative tolerance of the momentum solves. A time step needs
    /// accurate solves (their error is not iterated away), a steady outer
    /// iteration does not.
    pub fn momentum_tol(&self, transient: bool) -> f64 {
        self.momentum
            .tolerance
            .unwrap_or(if transient { 1e-8 } else { 1e-2 })
    }

    /// Relative tolerance of the pressure solves (see
    /// [`momentum_tol`](Self::momentum_tol)).
    pub fn pressure_tol(&self, transient: bool) -> f64 {
        self.pressure
            .tolerance
            .unwrap_or(if transient { 1e-8 } else { 1e-2 })
    }

    /// PISO with the SIMPLEC-type corrector ([`PisoCorrector`]).
    pub fn simplec_corrector(&self) -> bool {
        match self.corrector {
            PisoCorrector::Auto => self.outer_iterations > 1,
            PisoCorrector::Issa => false,
            PisoCorrector::Simplec => true,
        }
    }

    /// α_u. With the coupled algorithm it only caps the local Courant
    /// number of the pseudo-time step at α/(1 − α).
    pub fn alpha_u(&self, transient: bool) -> f64 {
        self.momentum.relaxation.unwrap_or(if transient {
            1.0
        } else if self.algorithm(false) == Algorithm::Coupled {
            if self.reduction_defaults {
                COUPLED_RELAXATION_REDUCTION
            } else {
                COUPLED_RELAXATION
            }
        } else {
            0.7
        })
    }

    /// α_p.
    pub fn alpha_p(&self, transient: bool) -> f64 {
        self.pressure
            .relaxation
            .unwrap_or(match self.algorithm(transient) {
                Algorithm::Simplec | Algorithm::Piso | Algorithm::Coupled => 1.0,
                Algorithm::Simple => 0.3,
            })
    }

    /// The settings with the cell-quality remediation (D50) applied to a
    /// mesh whose largest face skewness is `max_skewness`, and the changes
    /// made (`key: old -> new`; none with `remediation = "off"`, on a mesh
    /// without poor faces or where the case already sets them).
    pub fn remedied(&self, max_skewness: f64) -> (Self, Vec<String>) {
        let mut s = self.clone();
        let mut changes = Vec::new();
        if self.remediation == Remediation::Off || max_skewness <= self.skewness_limit {
            return (s, changes);
        }
        if s.pressure_gradient == PressureGradient::LeastSquares {
            s.pressure_gradient = PressureGradient::Gauss;
            changes.push("solver.pressure_gradient: \"least-squares\" -> \"gauss\"".into());
        }
        if s.momentum.scheme != Scheme::Upwind && s.momentum.limiter == LimiterKind::None {
            s.momentum.limiter = LimiterKind::Venkatakrishnan;
            changes.push("solver.momentum.limiter: \"none\" -> \"venkatakrishnan\"".into());
        }
        if s.face_velocity == FaceVelocity::Interpolated {
            s.face_velocity = FaceVelocity::Bounded;
            changes.push("solver.face_velocity: \"interpolated\" -> \"bounded\"".into());
        }
        (s, changes)
    }

    /// Checks the ranges; `transient` when the case has a `[time]` table.
    pub fn validate(&self, transient: bool) -> Result<()> {
        let alg = self.algorithm(transient);
        let au = self.alpha_u(transient);
        crate::need!(
            au > 0.0 && au <= 1.0,
            CASE_RANGE,
            "[solver.momentum] relaxation must lie in (0, 1]"
        );
        let ap = self.alpha_p(transient);
        crate::need!(
            ap > 0.0 && ap <= 1.0,
            CASE_RANGE,
            "[solver.pressure] relaxation must lie in (0, 1]"
        );
        if transient {
            crate::need!(
                alg == Algorithm::Piso,
                CASE_INCONSISTENT,
                "[solver] algorithm = \"{}\" is a steady algorithm; a case with a [time] table \
                 runs PISO (remove `algorithm` or set it to \"piso\")",
                match alg {
                    Algorithm::Simple => "simple",
                    Algorithm::Coupled => "coupled",
                    _ => "simplec",
                }
            );
            crate::need!(
                self.correctors >= 1,
                CASE_RANGE,
                "[solver] correctors must be >= 1"
            );
            crate::need!(
                self.outer_iterations >= 1,
                CASE_RANGE,
                "[solver] outer_iterations must be >= 1"
            );
            crate::need!(
                self.outer_tolerance > 0.0,
                CASE_RANGE,
                "[solver] outer_tolerance must be > 0"
            );
            crate::need!(
                !(self.outer_iterations == 1 && au < 1.0),
                CASE_INCONSISTENT,
                "[solver.momentum] relaxation < 1 needs outer_iterations > 1 in a transient case \
                 (with one outer iteration the relaxation would lag the solution)"
            );
        } else {
            crate::need!(
                alg != Algorithm::Piso,
                CASE_INCONSISTENT,
                "[solver] algorithm = \"piso\" is transient and needs a [time] table"
            );
            self.pseudo_time.validate()?;
            if alg == Algorithm::Coupled {
                self.coupled.validate()?;
                crate::need!(
                    au < 1.0,
                    CASE_RANGE,
                    "[solver.momentum] the coupled algorithm needs relaxation < 1 \
                     (it caps the local Courant number of the pseudo-time step at \
                     α/(1 − α)); try 0.98"
                );
            }
            crate::diag::check_fix(
                !(alg == Algorithm::Simplec && au >= 1.0),
                crate::diag::CASE_INCONSISTENT,
                &["solver.momentum.relaxation=0.7"],
                || {
                    String::from(
                        "[solver.momentum] SIMPLEC needs relaxation < 1 for a steady case \
                 (D^C = V/(a_P/α − Σ|a_N|) is singular at α = 1); try 0.7",
                    )
                },
            )?;
        }
        crate::need!(
            self.tolerance > 0.0,
            CASE_RANGE,
            "[solver] tolerance must be > 0"
        );
        crate::need!(
            self.momentum.limiter_k >= 0.0,
            CASE_RANGE,
            "[solver.momentum] limiter_k must be >= 0"
        );
        crate::need!(
            self.skewness_limit > 0.0,
            CASE_RANGE,
            "[solver] skewness_limit must be > 0"
        );
        crate::need!(
            self.operator_tolerance > 0.0,
            CASE_RANGE,
            "[solver] operator_tolerance must be > 0"
        );
        crate::need!(
            (0.0..=1.0).contains(&self.momentum.central_blending),
            CASE_RANGE,
            "[solver.momentum] central_blending must lie in [0, 1]"
        );
        crate::need!(
            (0.0..1.0).contains(&self.momentum.central_bounding),
            CASE_RANGE,
            "[solver.momentum] central_bounding must lie in [0, 1)"
        );
        crate::need!(
            self.momentum.limiter != LimiterKind::Mlp,
            CASE_RANGE,
            "[solver.momentum] limiter = \"mlp\" is not built for the velocity yet (scalars and \
             energy only)"
        );
        for (name, t) in [
            ("[solver.momentum] tolerance", self.momentum_tol(transient)),
            ("[solver.pressure] tolerance", self.pressure_tol(transient)),
        ] {
            crate::need!(t > 0.0 && t < 1.0, CASE_RANGE, "{name} must lie in (0, 1)");
        }
        crate::need!(
            self.amg_rebuild >= 1,
            CASE_RANGE,
            "[solver] amg_rebuild must be >= 1"
        );
        crate::need!(
            self.log_every >= 1,
            CASE_RANGE,
            "[solver] log_every must be >= 1"
        );
        self.pressure.amg.validate().map_err(|e| {
            let m = e.to_string();
            let mut d = crate::diag::Diagnostic::new(crate::diag::CASE_RANGE, m.clone());
            d.key = crate::diag::key_of(&m).map(|k| format!("solver.{k}"));
            crate::diag::err(d)
        })
    }
}

/// Physical inputs.
#[derive(Debug, Clone)]
pub struct Physics {
    /// Density ρ (kg/m³).
    pub rho: f64,
    /// Dynamic viscosity μ (Pa s); the reference value when it varies.
    pub mu: f64,
    /// Uniform body force per unit volume (N/m³); a body force given by
    /// expressions (D53) is a momentum source of the [`Flow`] instead
    /// ([`Flow::set_momentum_source`]) and this is zero.
    pub body_force: [f64; 3],
    /// Convection on (false: Stokes flow).
    pub convection: bool,
    /// Viscosity per face (D94; `None`: `mu` on every face, the constant
    /// path of the kernels bit for bit).
    pub mu_face: Option<std::sync::Arc<Vec<f64>>>,
    /// The momentum source per unit volume (N/m³) per cell of a body force
    /// given by expressions (D53), as last set ([`Flow::set_momentum_source`]).
    pub source: Option<std::sync::Arc<[Vec<f64>; 3]>>,
    /// The balanced buoyancy force B_P per unit volume (N/m³) per cell of
    /// the momentum rows (stage 7b1, D154), as last refreshed.
    pub buoyancy: Option<std::sync::Arc<[Vec<f64>; 3]>>,
    /// The variable density per cell and per face (stage 7c), as last
    /// refreshed (`None`: `rho` everywhere).
    pub rho_var: Option<std::sync::Arc<(Vec<f64>, Vec<f64>)>>,
    /// The density differs between the zones of several fluids (D194): the
    /// coupled solver takes the pseudo-time term and the scale of the
    /// continuity row per cell (false: `rho`, as before).
    pub zones: bool,
}

impl Physics {
    /// Constant properties (the face viscosity `mu` everywhere).
    pub fn constant(rho: f64, mu: f64, body_force: [f64; 3], convection: bool) -> Self {
        Self {
            rho,
            mu,
            body_force,
            convection,
            mu_face: None,
            source: None,
            buoyancy: None,
            rho_var: None,
            zones: false,
        }
    }

    /// The density of cell `c`.
    #[inline]
    pub fn rho_at(&self, c: usize) -> f64 {
        self.rho_var.as_ref().map_or(self.rho, |r| r.0[c])
    }

    /// The density of face `f`.
    #[inline]
    pub fn rho_face_at(&self, f: usize) -> f64 {
        self.rho_var.as_ref().map_or(self.rho, |r| r.1[f])
    }

    /// The viscosity of face `f`.
    #[inline]
    pub fn mu_at(&self, f: usize) -> f64 {
        self.mu_face.as_ref().map_or(self.mu, |m| m[f])
    }
}

/// Residuals of one outer iteration.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct Residuals {
    /// Outer iteration (from 1).
    pub iteration: usize,
    /// Momentum residuals of u, v, w.
    pub momentum: [f64; 3],
    /// Continuity residual of the predicted fluxes.
    pub continuity: f64,
    /// BiCGStab iterations of the three momentum solves.
    pub momentum_iterations: [usize; 3],
    /// CG iterations of the pressure solve (all correctors).
    pub pressure_iterations: usize,
    /// Krylov iterations of the coupled solve (0 with a segregated
    /// algorithm).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub coupled_iterations: usize,
    /// The pseudo-time step Δτ of the momentum equations in s (coupled
    /// algorithm, D177); none where only the relaxation acts.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pseudo_step: Option<f64>,
    /// The passive scalars (D51): residuals and iterations per scalar, in
    /// name order.
    #[serde(
        default,
        skip_serializing_if = "crate::scalar::ScalarResiduals::is_empty"
    )]
    pub scalars: crate::scalar::ScalarResiduals,
    /// The turbulence model (stage 6): residuals of k and ω.
    #[serde(
        default,
        skip_serializing_if = "crate::turbulence::TurbResiduals::is_empty"
    )]
    pub turbulence: crate::turbulence::TurbResiduals,
}

fn is_zero(n: &usize) -> bool {
    *n == 0
}

impl Residuals {
    /// Every residual is finite (a NaN does not hide behind
    /// [`max`](Self::max), which skips it).
    pub fn finite(&self) -> bool {
        self.continuity.is_finite()
            && self.momentum.iter().all(|r| r.is_finite())
            && self.scalars.residual.iter().all(|r| r.is_finite())
            && self.turbulence.residual.iter().all(|r| r.is_finite())
    }

    /// The largest residual.
    pub fn max(&self) -> f64 {
        self.momentum
            .iter()
            .chain(&self.scalars.residual)
            .chain(&self.turbulence.residual)
            .copied()
            .fold(self.continuity, f64::max)
    }
}

/// Wall-clock seconds per phase (R6).
#[derive(Debug, Clone, Default, Serialize)]
pub struct Timings {
    /// Boundary values, gradients, limiters.
    pub gradients: f64,
    /// Momentum faces and rows.
    pub momentum_assembly: f64,
    /// Momentum solves.
    pub momentum_solve: f64,
    /// Rhie–Chow fluxes, pressure rows, corrections.
    pub pressure_assembly: f64,
    /// Pressure solves.
    pub pressure_solve: f64,
    /// Multigrid setups.
    pub amg_setup: f64,
    /// Passive scalars (all their phases).
    pub scalars: f64,
    /// Turbulence model (all its phases).
    pub turbulence: f64,
}

impl Timings {
    /// Total seconds.
    pub fn total(&self) -> f64 {
        self.gradients
            + self.momentum_assembly
            + self.momentum_solve
            + self.pressure_assembly
            + self.pressure_solve
            + self.amg_setup
            + self.scalars
            + self.turbulence
    }
}

/// Host copy of the solution fields.
#[derive(Debug, Clone)]
pub struct State {
    /// Cell velocity.
    pub u: Soa3,
    /// Cell pressure.
    pub p: Vec<f64>,
    /// Face mass flux (kg/s) out of the owner.
    pub mass_flux: Vec<f64>,
    /// Boundary face velocity.
    pub ub: Soa3,
    /// Boundary face pressure.
    pub pb: Vec<f64>,
    /// Velocity gradients.
    pub grad_u: [Soa3; 3],
    /// Pressure gradient.
    pub grad_p: Soa3,
    /// Limiters of the velocity components.
    pub psi: Soa3,
    /// The passive scalars (D51), in name order.
    pub scalars: Vec<crate::scalar::ScalarState>,
    /// The driving force β (N/m³) of each periodic pair (D29), along its
    /// translation; empty without periodic pairs.
    pub periodic_force: Vec<f64>,
    /// The turbulence model's fields (stage 6); `None` in laminar flow.
    /// With a model, `p` and `pb` are the static pressure and the solver's
    /// modified pressure p̂ = p + ⅔ρk is kept here.
    pub turbulence: Option<crate::turbulence::TurbState>,
}

/// Host copies of the boundary values, gradients and limiters: what an
/// outer iteration takes over from the previous one besides the cell
/// values and the face fluxes ([`Flow::lagged`]).
pub struct Lagged {
    ub: Soa3,
    pb: Vec<f64>,
    grad_u: [Soa3; 3],
    grad_p: Soa3,
    psi: Soa3,
}

/// A host copy of the fields with the discretisation and the physics, for
/// the output.
pub struct Solution<'a> {
    /// The fields.
    pub state: State,
    /// Discretisation.
    pub disc: &'a Discretisation,
    /// Physics.
    pub phys: Physics,
}

impl Solution<'_> {
    /// The face viscosity of the momentum equation: μ + ρν_t with a
    /// turbulence model (μ_w of the wall treatment on the walls), else μ.
    pub fn mu_eff_at(&self, f: usize) -> f64 {
        match &self.state.turbulence {
            Some(t) => t.mu_eff[f],
            None => self.phys.mu_at(f),
        }
    }
}

impl std::ops::Deref for Solution<'_> {
    type Target = State;
    fn deref(&self) -> &State {
        &self.state
    }
}

/// The field buffers that carry the solution from one iteration (or time
/// step) to the next.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Fields {
    pub u: Buf3,
    pub p: Buf,
    pub mass_flux: Buf,
    pub ub: Buf3,
    pub ub_set: Buf3,
    pub pb: Buf,
    pub grad_u: [Buf3; 3],
    pub grad_p: Buf3,
    pub psi: Buf3,
}

/// Work buffers of one iteration.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Work {
    pub zero_face: Buf,
    pub relax: Buf,
    pub cp: Buf,
    pub cn: Buf,
    pub bf: Buf3,
    pub val: Buf,
    pub rhs: Buf3,
    pub diag: Buf3,
    pub d: Buf,
    pub ds: Buf,
    pub dc: Buf,
    pub m: Buf,
    pub defect: Buf,
    pub through: Buf,
    pub cf: Buf,
    pub pval: Buf,
    pub prhs: Buf,
    pub pc: Buf,
    pub mc: Buf,
    pub gpc: Buf3,
    pub zero_boundary: Buf,
}

/// The flow solver state on one partition.
pub struct Flow {
    /// Host copy of the discretisation (geometry for the output).
    pub disc: Discretisation,
    /// Boundary arrays.
    pub bc: BoundaryArrays,
    /// Physics.
    pub phys: Physics,
    /// Settings.
    pub settings: Settings,
    /// Reference cell of the pressure (no pressure boundary).
    pub reference_cell: Option<usize>,
    /// Phase timings.
    pub timings: Timings,
    /// Multigrid setups done.
    pub amg_setups: usize,
    /// The pressure matrix of the last multigrid setup (for a restart that
    /// continues bit-identically).
    pub amg_matrix: Option<Vec<f64>>,
    pub(crate) dev: Box<dyn FlowDevice>,
    pub(crate) f: Fields,
    pub(crate) w: Work,
    pub(crate) transient: bool,
    pub(crate) iteration: usize,
    /// Device copies of the last valid u, p and ṁ and its iteration, for
    /// the divergence report (D63).
    valid: Option<(Buf3, Buf, Buf, usize)>,
    /// The residuals of the outer iteration or corrector that diverged.
    pub(crate) diverged: Option<Residuals>,
    /// The velocity of the central correction when it is not the current
    /// one (the extrapolated velocity of a time step, set by
    /// [`crate::transient`]).
    pub(crate) central_u: Option<Buf3>,
    /// Some boundary face suppresses backflow (D29).
    no_backflow: bool,
    /// The convective outlets (D29), if any.
    pub(crate) conv: Option<Convective>,
    /// The momentum source S·V per cell (N) of a body force given by
    /// expressions (D53), added to the right-hand side of the momentum
    /// rows; `None` for a uniform body force (inside the rows).
    pub(crate) source: Option<Buf3>,
    /// The cell-quality remediation (D50): the mesh check and the changes.
    pub remediation: RemediationReport,
    /// The passive scalars (D51), solved after the pressure correction of
    /// every outer iteration ([`Flow::add_scalars`]).
    pub(crate) scalars: Option<crate::scalar::Scalars>,
    /// The energy equation on the full mesh of a conjugate heat transfer
    /// case (stage 7d1): this flow runs on the fluid view.
    pub(crate) cht: Option<Box<crate::cht::Cht>>,
    /// The turbulence model (stage 6), solved before the scalars.
    pub(crate) turbulence: Option<crate::turbulence::Turbulence>,
    /// The driving of the periodic pairs (D29), one per pair.
    pub(crate) periodic: Vec<PeriodicDriving>,
    /// Work face buffer of the mass-flow driving (the flux of a unit
    /// increment of the driving force).
    unit_flux: Option<Buf>,
    /// Buffers of the coupled algorithm (allocated at its first iteration).
    pub(crate) coupled: Option<crate::coupled::CoupledWork>,
    /// Outer iterations of the coupled algorithm in which the multigrid
    /// setup of the momentum sweep failed and BiCGStab did the sweeps.
    pub sweep_fallbacks: usize,
    /// The message of the first such failure, until the run loop reports
    /// it.
    pub sweep_note: Option<String>,
    /// Multigrid setups of the momentum sweeps of the coupled algorithm.
    pub sweep_setups: usize,
    /// The hierarchies of the momentum sweeps kept between outer
    /// iterations (the content of a restart file).
    pub sweep: crate::coupled::SweepState,
    /// The length of the automatic pseudo-time step (computed at its first
    /// use).
    pub(crate) pseudo_length: Option<f64>,
    /// The largest prescribed boundary speed for that step (computed at
    /// its first use, again after new boundary velocities).
    pub(crate) pseudo_speed: Option<f64>,
    /// The rotating zones of the moving reference frame (D55), if any.
    pub(crate) frames: Option<crate::frame::Frames>,
    /// The host copy of the least-squares weights of the pressure (the
    /// balanced force of gravity applies the same operator, D154).
    pub(crate) lsq_p: foehncraft_fv::disc::Soa3,
    /// The balanced-force arrays of the Boussinesq buoyancy (stage 7b1).
    pub buoyancy: Option<Box<crate::buoyancy::Balanced>>,
    /// The device buffer of the buoyancy source B_P V of the momentum rows.
    buoyancy_source: Option<Buf3>,
    /// The variable density ρ(T, p₀) of stage 7c.
    pub density: Option<Box<crate::density::Variable>>,
    /// The backward-difference coefficients (a₀, −a₁, −a₂, as the
    /// momentum time terms) and Δt of the current time step (`None` in a
    /// steady run): the mass source V ∂ρ/∂t of the continuity.
    pub(crate) density_time: Option<([f64; 3], f64)>,
    /// The mass source V ∂ρ/∂t per cell on the device (the continuity
    /// residual adds it).
    density_source: Option<Buf>,
}

/// The driving of one periodic pair: a uniform body force β t̂ along the
/// translation, fixed (`pressure_gradient`) or adjusted after every
/// pressure correction so that the mass flow through side `b` is the
/// prescribed one (`mass_flow`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct PeriodicDriving {
    pub drive: crate::bc::Drive,
    pub direction: [f64; 3],
    /// 1 on the faces of side `b`, 0 elsewhere (a face field).
    pub plane: Buf,
    /// β (N/m³).
    pub beta: f64,
}

/// The convective outlets of a transient run (D29): the device buffers of
/// the outlet equation and of the global mass balance.
///
/// Per time step [`crate::transient`] sets the old-level part c and β of
/// [`foehncraft_kernels::flow::UB_CONVECTIVE`]. When no boundary fixes the
/// pressure, the outlet velocities are shifted along the face normals after
/// every boundary update, by one common speed, so that the outflow through
/// the convective outlets equals the prescribed inflow exactly (the
/// pressure-correction equation is then compatible; its level is fixed in
/// the reference cell). Between time steps the outlets *hold*: β = 0 and
/// c = u_b, so a boundary update reproduces the boundary values of the last
/// step bitwise (a restart, a rejected step, a rewind).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Convective {
    /// c (per boundary face; the device buffer the boundary values read).
    pub c: Buf3,
    /// β (per boundary face).
    pub beta: Buf,
    /// U_c/Δn on the faces of outlets with a given convection velocity.
    given_dn: Buf,
    /// 1/Δn on the faces of outlets with the default convection velocity.
    default_dn: Buf,
    /// Area vectors of the convective faces (0 elsewhere).
    area: Buf3,
    /// Unit normals of the convective faces (0 elsewhere).
    normal: Buf3,
    /// Total area of the convective faces (m²).
    area_total: f64,
    /// Mass flux out of the domain through the other boundaries that fix
    /// the flux (kg/s; negative for a net inflow).
    pub q_other: f64,
    /// Shift the outlet velocities to the mass balance (no pressure
    /// boundary).
    balance: bool,
    /// Between time steps: β = 0, c = u_b, no balance shift.
    pub hold: bool,
}

impl Convective {
    /// The default convection velocity: the mean outflow speed through the
    /// convective outlets that carries the prescribed net inflow, ≥ 0.
    pub fn default_velocity(&self, rho: f64) -> f64 {
        if self.area_total > 0.0 {
            (-self.q_other / (rho * self.area_total)).max(0.0)
        } else {
            0.0
        }
    }
}

pub(crate) fn disc_with(
    mesh: &Mesh,
    geo: &Geometry,
    points: &std::collections::BTreeMap<String, BoundaryKind>,
) -> Result<Discretisation> {
    Discretisation::new(mesh, geo, &stencil_problem(mesh, points))
}

/// The least-squares weights of [`disc_with`] alone.
pub(crate) fn lsq_with(
    mesh: &Mesh,
    geo: &Geometry,
    points: &std::collections::BTreeMap<String, BoundaryKind>,
) -> Result<Soa3> {
    Discretisation::least_squares(mesh, geo, &stencil_problem(mesh, points))
}

/// The boundary stencil points as a problem: a fixed value where the point
/// is the face centroid, a zero gradient where it is projected.
fn stencil_problem(
    mesh: &Mesh,
    points: &std::collections::BTreeMap<String, BoundaryKind>,
) -> Problem {
    let mut prob = Problem::default();
    for p in &mesh.boundaries {
        let kind = points
            .get(&p.name)
            .copied()
            .unwrap_or(BoundaryKind::FixedValue);
        let b = match kind {
            BoundaryKind::FixedValue => Boundary::value(0.0),
            _ => Boundary::zero_gradient(),
        };
        prob.boundary.insert(p.name.clone(), b);
    }
    prob
}

/// The discretisations of the flow solver: the geometry, the matrix
/// structure and the least-squares weights of the velocity (boundary
/// stencil points at the face centroids, projected onto the face normal at
/// outlets and symmetry planes), and the least-squares weights of the
/// pressure (projected at walls and symmetry planes, where the pressure has
/// a zero normal gradient).
pub fn discretisation(
    mesh: &Mesh,
    geo: &Geometry,
    bc: &BoundaryArrays,
) -> Result<(Discretisation, Soa3)> {
    let d = disc_with(mesh, geo, &bc.u_points)?;
    let lsq_p = lsq_with(mesh, geo, &bc.p_points)?;
    Ok((d, lsq_p))
}

/// Where a flow runs.
#[derive(Debug, Clone)]
pub struct Target {
    /// Device name: `cpu` or `cuda`.
    pub device: String,
    /// CPU threads.
    pub threads: Option<usize>,
    /// PTX module of the CUDA device.
    pub ptx: std::path::PathBuf,
}

impl Default for Target {
    fn default() -> Self {
        Self {
            device: "cpu".into(),
            threads: None,
            ptx: std::path::PathBuf::new(),
        }
    }
}

impl Flow {
    /// A solver with the fields at `u0`, `p0` and the fluxes of `u0`;
    /// `reference_cell` fixes the pressure level when no boundary does.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        disc: (Discretisation, Soa3),
        bc: BoundaryArrays,
        phys: Physics,
        settings: Settings,
        u0: [f64; 3],
        p0: f64,
        reference_cell: Option<usize>,
        target: &Target,
    ) -> Result<Self> {
        Self::on(disc, bc, phys, settings, false, reference_cell, target).and_then(|mut s| {
            let nc = s.disc.n_cells;
            s.initialise(
                &Soa3 {
                    x: vec![u0[0]; nc],
                    y: vec![u0[1]; nc],
                    z: vec![u0[2]; nc],
                },
                &vec![p0; nc],
            )?;
            Ok(s)
        })
    }

    /// A solver on `target` with zero fields; set them with
    /// [`initialise`](Self::initialise) or [`load`](Self::load).
    pub fn on(
        (disc, lsq_p): (Discretisation, Soa3),
        mut bc: BoundaryArrays,
        phys: Physics,
        settings: Settings,
        transient: bool,
        reference_cell: Option<usize>,
        target: &Target,
    ) -> Result<Self> {
        settings.validate(transient)?;
        // cell-quality remediation (D50); with several ranks (stage 7) the
        // largest skewness is to be reduced over the ranks, so that all of
        // them take the same schemes
        let skew = face_skewness(&disc);
        let max_skewness = skew.iter().copied().fold(0.0, f64::max);
        let (settings, changes) = settings.remedied(max_skewness);
        let remediation = {
            let lim = settings.skewness_limit;
            let mut poor = vec![false; disc.n_cells];
            let mut skewed_faces = 0;
            for (f, &v) in skew.iter().enumerate() {
                if v > lim {
                    skewed_faces += 1;
                    poor[disc.owner[f] as usize] = true;
                    let n = disc.neighbour[f] as usize;
                    if n < disc.n_cells {
                        poor[n] = true;
                    }
                }
            }
            RemediationReport {
                skewness_limit: lim,
                skewed_faces,
                poor_cells: poor.iter().filter(|&&p| p).count(),
                max_skewness,
                changes,
            }
        };
        // the Gauss pressure gradient: its weights, every boundary pressure
        // enters
        let lsq_p = if settings.pressure_gradient == PressureGradient::Gauss {
            bc.p_zero_gradient.iter_mut().for_each(|a| *a = 0.0);
            gauss_weights(&disc)
        } else {
            lsq_p
        };
        let nc = disc.n_cells;
        let nf = disc.n_faces;
        let nb = nf - disc.n_interior;
        // the per-cell ε² of the Venkatakrishnan limiter, also when only the
        // bounded face velocity (D50) takes it (`gradients` uses that limiter
        // when the convection scheme has none)
        let venkat = match settings.momentum.limiter {
            LimiterKind::Venkatakrishnan => true,
            LimiterKind::None => settings.face_velocity == FaceVelocity::Bounded,
            LimiterKind::VenkatakrishnanWang | LimiterKind::BarthJespersen | LimiterKind::Mlp => {
                false
            }
        };
        let eps2 = if venkat {
            disc.venkat_eps2(settings.momentum.limiter_k)
        } else {
            Vec::new()
        };
        let mut dev = device::create(
            &target.device,
            &Setup {
                disc: &disc,
                lsq_p: &lsq_p,
                bc: &bc,
                rho: phys.rho,
                mu: phys.mu,
                mu_face: phys.mu_face.as_deref().map(|v| v.as_slice()),
                eps2: &eps2,
                threads: target.threads,
            },
            &target.ptx,
        )?;
        // a density per cell and face (D194: several fluids; 7c sets its
        // own after the energy equation)
        if let Some(r) = &phys.rho_var {
            dev.set_density(&r.0, &r.1)?;
        }
        let stored = disc.sell.stored();
        let cell3 = |dev: &mut Box<dyn FlowDevice>| -> Result<Buf3> {
            Ok(Buf3 {
                x: dev.alloc(nc)?,
                y: dev.alloc(nc)?,
                z: dev.alloc(nc)?,
            })
        };
        let u = cell3(&mut dev)?;
        let grad_u = [cell3(&mut dev)?, cell3(&mut dev)?, cell3(&mut dev)?];
        let grad_p = cell3(&mut dev)?;
        let psi = cell3(&mut dev)?;
        for b in psi.all() {
            dev.fill(b, 1.0)?;
        }
        let ub_set = Buf3 {
            x: dev.alloc(nb)?,
            y: dev.alloc(nb)?,
            z: dev.alloc(nb)?,
        };
        dev.upload(ub_set.x, &bc.ub[0])?;
        dev.upload(ub_set.y, &bc.ub[1])?;
        dev.upload(ub_set.z, &bc.ub[2])?;
        let ub = Buf3 {
            x: dev.alloc(nb)?,
            y: dev.alloc(nb)?,
            z: dev.alloc(nb)?,
        };
        let f = Fields {
            u,
            p: dev.alloc(nc)?,
            mass_flux: dev.alloc(nf)?,
            ub,
            ub_set,
            pb: dev.alloc(nb)?,
            grad_u,
            grad_p,
            psi,
        };
        let w = Work {
            zero_face: dev.alloc(nf)?,
            relax: dev.alloc(nf)?,
            cp: dev.alloc(nf)?,
            cn: dev.alloc(nf)?,
            bf: Buf3 {
                x: dev.alloc(nf)?,
                y: dev.alloc(nf)?,
                z: dev.alloc(nf)?,
            },
            val: dev.alloc(stored)?,
            rhs: cell3(&mut dev)?,
            diag: cell3(&mut dev)?,
            d: dev.alloc(nc)?,
            ds: dev.alloc(nc)?,
            dc: dev.alloc(nc)?,
            m: dev.alloc(nf)?,
            defect: dev.alloc(nc)?,
            through: dev.alloc(nc)?,
            cf: dev.alloc(nf)?,
            pval: dev.alloc(stored)?,
            prhs: dev.alloc(nc)?,
            pc: dev.alloc(nc)?,
            mc: dev.alloc(nf)?,
            gpc: cell3(&mut dev)?,
            zero_boundary: dev.alloc(nb)?,
        };
        let reference_cell = if bc.has_pressure() {
            None
        } else {
            Some(reference_cell.unwrap_or(0))
        };
        let no_backflow = bc.flux_kind.contains(&k::FLUX_PRESSURE_NO_BACKFLOW);
        let conv = if bc.ub_rule.contains(&k::UB_CONVECTIVE) {
            let (c, beta) = dev.convective_buffers();
            let ni = disc.n_interior;
            let mut given = vec![0.0; nb];
            let mut default = vec![0.0; nb];
            let mut area = [vec![0.0; nb], vec![0.0; nb], vec![0.0; nb]];
            let mut normal = [vec![0.0; nb], vec![0.0; nb], vec![0.0; nb]];
            let mut area_total = 0.0;
            for fb in 0..nb {
                if bc.ub_rule[fb] != k::UB_CONVECTIVE {
                    continue;
                }
                let f = fb + ni;
                let sv = [
                    disc.face_area.x[f],
                    disc.face_area.y[f],
                    disc.face_area.z[f],
                ];
                let a = (sv[0] * sv[0] + sv[1] * sv[1] + sv[2] * sv[2]).sqrt();
                let o = disc.owner[f] as usize;
                // the distance from the owner centroid to the face plane
                let dn = (0..3)
                    .map(|i| {
                        let xf = [
                            disc.face_centroid.x[f],
                            disc.face_centroid.y[f],
                            disc.face_centroid.z[f],
                        ];
                        let xc = [
                            disc.cell_centroid.x[o],
                            disc.cell_centroid.y[o],
                            disc.cell_centroid.z[o],
                        ];
                        (xf[i] - xc[i]) * sv[i] / a
                    })
                    .sum::<f64>();
                anyhow::ensure!(
                    dn > 0.0,
                    "a convective-outlet face lies behind the centroid of its cell"
                );
                if bc.convection[fb] > 0.0 {
                    given[fb] = bc.convection[fb] / dn;
                } else {
                    default[fb] = 1.0 / dn;
                }
                for i in 0..3 {
                    area[i][fb] = sv[i];
                    normal[i][fb] = sv[i] / a;
                }
                area_total += a;
            }
            let mut put = |v: &[f64]| -> Result<Buf> {
                let b = dev.alloc(nb)?;
                dev.upload(b, v)?;
                Ok(b)
            };
            let given_dn = put(&given)?;
            let default_dn = put(&default)?;
            let area = Buf3 {
                x: put(&area[0])?,
                y: put(&area[1])?,
                z: put(&area[2])?,
            };
            let normal = Buf3 {
                x: put(&normal[0])?,
                y: put(&normal[1])?,
                z: put(&normal[2])?,
            };
            Some(Convective {
                c,
                beta,
                given_dn,
                default_dn,
                area,
                normal,
                area_total,
                q_other: fixed_outflow(&disc, &bc, &bc.ub, &phys),
                balance: !bc.has_pressure(),
                hold: true,
            })
        } else {
            None
        };
        let mut periodic = Vec::new();
        for d in &bc.drives {
            let mut v = vec![0.0; nf];
            v[d.b_faces.clone()].iter_mut().for_each(|x| *x = 1.0);
            let plane = dev.alloc(nf)?;
            dev.upload(plane, &v)?;
            periodic.push(PeriodicDriving {
                drive: d.drive,
                direction: d.direction,
                plane,
                beta: match d.drive {
                    crate::bc::Drive::PressureGradient(g) => g,
                    _ => 0.0,
                },
            });
        }
        let unit_flux = if periodic
            .iter()
            .any(|d| matches!(d.drive, crate::bc::Drive::MassFlow(_)))
        {
            Some(dev.alloc(nf)?)
        } else {
            None
        };
        Ok(Self {
            disc,
            bc,
            phys,
            settings,
            reference_cell,
            timings: Timings::default(),
            amg_setups: 0,
            amg_matrix: None,
            dev,
            f,
            w,
            transient,
            iteration: 0,
            valid: None,
            diverged: None,
            central_u: None,
            no_backflow,
            conv,
            source: None,
            remediation,
            scalars: None,
            cht: None,
            turbulence: None,
            periodic,
            unit_flux,
            coupled: None,
            sweep_fallbacks: 0,
            sweep_note: None,
            sweep_setups: 0,
            sweep: Default::default(),
            pseudo_length: None,
            pseudo_speed: None,
            frames: None,
            lsq_p,
            buoyancy: None,
            buoyancy_source: None,
            density: None,
            density_time: None,
            density_source: None,
        })
    }

    /// The uniform body force (N/m³): the one of the case plus the driving
    /// of the periodic pairs (D29).
    pub fn body_force(&self) -> [f64; 3] {
        let mut b = self.phys.body_force;
        for d in &self.periodic {
            for (bi, ti) in b.iter_mut().zip(d.direction) {
                *bi += d.beta * ti;
            }
        }
        b
    }

    /// The physics with the driving of the periodic pairs in the body
    /// force (for the output and the balances).
    pub fn physics(&self) -> Physics {
        Physics {
            body_force: self.body_force(),
            ..self.phys.clone()
        }
    }

    /// The driving force β (N/m³) of each periodic pair, along its
    /// translation.
    pub fn periodic_force(&self) -> Vec<f64> {
        self.periodic.iter().map(|d| d.beta).collect()
    }

    /// The mass flow (kg/s) through side `b` of each periodic pair, from
    /// the current fluxes.
    pub fn periodic_mass_flow(&mut self) -> Result<Vec<f64>> {
        let m = self.f.mass_flux;
        let planes: Vec<Buf> = self.periodic.iter().map(|d| d.plane).collect();
        planes.into_iter().map(|p| self.dev.dot(m, p)).collect()
    }

    /// The mass-flow driving of the periodic pairs (D29), after the
    /// corrections of the fluxes `w.mc` and the velocity: β changes by δβ
    /// so that the flux through side `b` is the prescribed one. The
    /// velocity follows a uniform force increment by δβ D^C t̂ (the
    /// response of the SIMPLEC rows), the fluxes by the flux of that
    /// increment.
    pub(crate) fn drive_mass_flow(&mut self) -> Result<()> {
        let Some(unit) = self.unit_flux else {
            return Ok(());
        };
        let (f, w) = (self.f, self.w);
        for i in 0..self.periodic.len() {
            let d = self.periodic[i];
            let crate::bc::Drive::MassFlow(target) = d.drive else {
                continue;
            };
            for (g, t) in w.gpc.all().into_iter().zip(d.direction) {
                self.dev.lincomb(g, t, w.dc, 0.0, w.dc)?;
            }
            self.dev.delta_flux(w.gpc, unit)?;
            let k = self.dev.dot(unit, d.plane)?;
            let q = self.dev.dot(w.mc, d.plane)?;
            anyhow::ensure!(
                k > 0.0,
                "the periodic pair cannot be driven: no flux follows a force along its translation"
            );
            let db = (target - q) / k;
            self.dev.lincomb(w.mc, 1.0, w.mc, db, unit)?;
            for (u, g) in f.u.all().into_iter().zip(w.gpc.all()) {
                self.dev.lincomb(u, 1.0, u, db, g)?;
            }
            self.periodic[i].beta += db;
        }
        Ok(())
    }

    /// Sets the momentum source of a body force given by expressions
    /// (D53): `s` per unit volume (N/m³) per cell; uploaded as S·V.
    pub fn set_momentum_source(&mut self, s: [Vec<f64>; 3]) -> Result<()> {
        let nc = self.disc.n_cells;
        let b = match self.source {
            Some(b) => b,
            None => {
                let b = Buf3 {
                    x: self.dev.alloc(nc)?,
                    y: self.dev.alloc(nc)?,
                    z: self.dev.alloc(nc)?,
                };
                self.source = Some(b);
                b
            }
        };
        for (i, bi) in b.all().into_iter().enumerate() {
            let sv: Vec<f64> = s[i]
                .iter()
                .zip(&self.disc.cell_volume)
                .map(|(a, v)| a * v)
                .collect();
            self.dev.upload(bi, &sv)?;
        }
        self.phys.source = Some(std::sync::Arc::new(s));
        Ok(())
    }

    /// Keeps a device copy of the current u, p and ṁ as the last valid
    /// state (D63): three cell fields and one face field, copied on the
    /// device, against the dozens of passes of an iteration.
    pub fn keep_valid(&mut self) -> Result<()> {
        let (nc, nf) = (self.disc.n_cells, self.disc.n_faces);
        let v = match self.valid {
            Some(v) => v,
            None => {
                let u = Buf3 {
                    x: self.dev.alloc(nc)?,
                    y: self.dev.alloc(nc)?,
                    z: self.dev.alloc(nc)?,
                };
                (u, self.dev.alloc(nc)?, self.dev.alloc(nf)?, 0)
            }
        };
        for (a, b) in self.f.u.all().into_iter().zip(v.0.all()) {
            self.dev.copy(a, b)?;
        }
        self.dev.copy(self.f.p, v.1)?;
        self.dev.copy(self.f.mass_flux, v.2)?;
        self.valid = Some((v.0, v.1, v.2, self.iteration));
        if let Some(sc) = self.scalars.as_mut() {
            sc.keep_valid(self.dev.as_mut())?;
        }
        // the temperature of a conjugate case lives in the energy view
        if let Some((sc, dev)) = crate::cht::energy_fields(self) {
            sc.keep_valid(dev)?;
        }
        if let Some(t) = self.turbulence.as_mut() {
            t.keep_valid(self.dev.as_mut())?;
        }
        Ok(())
    }

    /// Puts the last valid state of [`keep_valid`](Self::keep_valid) back
    /// with its boundary values and gradients; false when none was kept.
    pub fn rewind(&mut self) -> Result<bool> {
        let Some((u, p, m, it)) = self.valid else {
            return Ok(false);
        };
        for (a, b) in u.all().into_iter().zip(self.f.u.all()) {
            self.dev.copy(a, b)?;
        }
        self.dev.copy(p, self.f.p)?;
        self.dev.copy(m, self.f.mass_flux)?;
        if let Some(sc) = self.scalars.as_mut() {
            sc.rewind(self.dev.as_mut())?;
        }
        if let Some((sc, dev)) = crate::cht::energy_fields(self) {
            sc.rewind(dev)?;
        }
        if let Some(t) = self.turbulence.as_mut() {
            t.rewind(self.dev.as_mut())?;
        }
        self.iteration = it;
        // the kept hierarchies of the momentum sweeps belong to later
        // matrices
        self.sweep = Default::default();
        self.refresh_buoyancy()?;
        self.update_boundaries()?;
        self.gradients()?;
        Ok(true)
    }

    /// Sets velocity and pressure and starts from the fluxes of the
    /// interpolated velocity (no pressure terms).
    pub fn initialise(&mut self, u: &Soa3, p: &[f64]) -> Result<()> {
        let started;
        let u = if self
            .periodic
            .iter()
            .any(|d| matches!(d.drive, crate::bc::Drive::MassFlow(_)))
        {
            started = self.with_mass_flow(u);
            &started
        } else {
            u
        };
        self.dev.upload(self.f.u.x, &u.x)?;
        self.dev.upload(self.f.u.y, &u.y)?;
        self.dev.upload(self.f.u.z, &u.z)?;
        self.dev.upload(self.f.p, p)?;
        if let Some(cv) = self.conv {
            // the convective outlets start from the owner velocity
            let (ni, nb) = (
                self.disc.n_interior,
                self.disc.n_faces - self.disc.n_interior,
            );
            let mut c = [vec![0.0; nb], vec![0.0; nb], vec![0.0; nb]];
            for (fb, &rule) in self.bc.ub_rule.iter().enumerate().take(nb) {
                if rule == k::UB_CONVECTIVE {
                    let o = self.disc.owner[fb + ni] as usize;
                    c[0][fb] = u.x[o];
                    c[1][fb] = u.y[o];
                    c[2][fb] = u.z[o];
                }
            }
            for (b, v) in cv.c.all().into_iter().zip(&c) {
                self.dev.upload(b, v)?;
            }
            self.dev.fill(cv.beta, 0.0)?;
            self.conv = Some(Convective { hold: true, ..cv });
        }
        self.update_boundaries()?;
        self.gradients()?;
        let a = self.flux_args(None, self.w.zero_face, None);
        self.dev.mass_flux(&a, self.f.mass_flux)
    }

    /// The start velocity `u` with a uniform velocity along the
    /// translation of each periodic pair driven by its mass flow (D29)
    /// added, so that the flow through side `b` (from the owner
    /// velocities) is the prescribed one.
    fn with_mass_flow(&self, u: &Soa3) -> Soa3 {
        let mut u = u.clone();
        let d = &self.disc;
        for (pd, bd) in self.periodic.iter().zip(&self.bc.drives) {
            let crate::bc::Drive::MassFlow(target) = pd.drive else {
                continue;
            };
            let (mut q, mut q_unit) = (0.0, 0.0);
            for f in bd.b_faces.clone() {
                let o = d.owner[f] as usize;
                let s = [d.face_area.x[f], d.face_area.y[f], d.face_area.z[f]];
                let rho = self.phys.rho_at(o);
                q += rho * (u.x[o] * s[0] + u.y[o] * s[1] + u.z[o] * s[2]);
                q_unit += rho
                    * (pd.direction[0] * s[0] + pd.direction[1] * s[1] + pd.direction[2] * s[2]);
            }
            if q_unit > 0.0 {
                let du = (target - q) / q_unit;
                for c in 0..d.n_cells {
                    u.x[c] += du * pd.direction[0];
                    u.y[c] += du * pd.direction[1];
                    u.z[c] += du * pd.direction[2];
                }
            }
        }
        u
    }

    /// Potential-flow start (D52): the velocity of the potential φ that
    /// carries the prescribed boundary fluxes through the domain, ∇²φ = 0
    /// with ∂φ/∂n from the velocity inlets, ∂φ/∂n = 0 on walls and symmetry
    /// planes and φ = 0 on pressure boundaries. It is one pressure-correction
    /// solve from zero cell velocities with D = 1, so the AMG of the
    /// pressure solves it with the same boundary conditions; the face fluxes
    /// are continuous to the solver tolerance. The cell velocity comes from
    /// the face fluxes, u_P = Σ_f ±ṁ_f (x_f − x_P)/(ρ V_P), which is exact
    /// for a uniform flow and, unlike the gradient of φ, also holds next to
    /// the inlets (where the correction gradient assumes ∂φ/∂n = 0). The
    /// pressure stays at its start value, the multigrid of the run is built
    /// anew. Returns the CG iterations.
    pub fn potential_flow(&mut self, mesh: &Mesh, geo: &Geometry) -> Result<usize> {
        let (f, w) = (self.f, self.w);
        for b in f.u.all() {
            self.dev.fill(b, 0.0)?;
        }
        let kept = self.settings.clone();
        self.settings.pressure.tolerance = Some(1e-10);
        self.settings.pressure.max_iter = self.settings.pressure.max_iter.max(500);
        self.update_boundaries()?;
        let a = self.flux_args(None, w.zero_face, None);
        self.dev.mass_flux(&a, f.mass_flux)?;
        self.dev.fill(w.dc, 1.0)?;
        // the boundary fluxes, and per pass the non-orthogonal part of the
        // face flux −ρ (∇φ)_f·k_f of the potential of the pass before
        // (deferred, Jasak 1996), so that the potential is that of the full
        // Laplacian on non-orthogonal meshes
        let m0 = self.dev.download(f.mass_flux)?;
        let phys = self.phys.clone();
        let mut it = 0;
        let mut r = Ok(());
        for pass in 0..=POTENTIAL_NONORTH_PASSES {
            if pass > 0 {
                self.dev.gradient(
                    Lsq::Pressure,
                    w.pc,
                    Mask::PressureCorrection,
                    w.zero_boundary,
                    w.gpc,
                )?;
                let g = [
                    self.dev.download(w.gpc.x)?,
                    self.dev.download(w.gpc.y)?,
                    self.dev.download(w.gpc.z)?,
                ];
                let d = &self.disc;
                let mut ms = m0.clone();
                for (fi, m) in ms.iter_mut().enumerate().take(d.n_interior) {
                    let (p, n, wt) = (d.owner[fi] as usize, d.neighbour[fi] as usize, d.weight[fi]);
                    let k = [d.k.x[fi], d.k.y[fi], d.k.z[fi]];
                    let gk: f64 = (0..3)
                        .map(|i| (wt * g[i][p] + (1.0 - wt) * g[i][n]) * k[i])
                        .sum();
                    *m -= phys.rho_face_at(fi) * gk;
                }
                self.dev.upload(f.mass_flux, &ms)?;
            }
            match self.pressure_correction(f.mass_flux, pass == 0) {
                Ok(n) => it += n,
                Err(e) => {
                    r = Err(e);
                    break;
                }
            }
        }
        self.settings = kept;
        r?;
        self.dev.correct_fluxes(w.cf, f.mass_flux, w.pc, w.mc)?;
        self.dev.copy(w.mc, f.mass_flux)?;
        let m = self.dev.download(f.mass_flux)?;
        let u = velocity_of_fluxes(mesh, geo, &m, &self.phys, self.disc.n_cells);
        self.dev.upload(f.u.x, &u.x)?;
        self.dev.upload(f.u.y, &u.y)?;
        self.dev.upload(f.u.z, &u.z)?;
        self.amg_setups = 0;
        self.amg_matrix = None;
        self.sweep = Default::default();
        self.update_boundaries()?;
        self.gradients()?;
        self.dev.sync()?;
        Ok(it)
    }

    /// Outer iterations done.
    pub fn iteration(&self) -> usize {
        self.iteration
    }

    /// Device name and description.
    pub fn device(&self) -> String {
        self.dev.device()
    }

    /// Host copy of the fields.
    pub fn state(&mut self) -> Result<State> {
        let f = self.f;
        let mut p = self.dev.download(f.p)?;
        let mut pb = self.dev.download(f.pb)?;
        let turbulence = match &self.turbulence {
            Some(t) => {
                let ts = t.state(self.dev.as_mut(), &p, &pb)?;
                // the static pressure p = p̂ − ⅔ρk
                for (c, (pi, ki)) in p.iter_mut().zip(&ts.k).enumerate() {
                    *pi -= 2.0 / 3.0 * self.phys.rho_at(c) * ki;
                }
                let kf = t.k_face(&ts.k, &ts.k_b, &self.disc.owner);
                let ni = self.disc.n_interior;
                for (fb, (pi, ki)) in pb.iter_mut().zip(&kf).enumerate() {
                    *pi -= 2.0 / 3.0 * self.phys.rho_face_at(ni + fb) * ki;
                }
                Some(ts)
            }
            None => None,
        };
        let mut v3 = |b: Buf3| -> Result<Soa3> {
            Ok(Soa3 {
                x: self.dev.download(b.x)?,
                y: self.dev.download(b.y)?,
                z: self.dev.download(b.z)?,
            })
        };
        Ok(State {
            u: v3(f.u)?,
            ub: v3(f.ub)?,
            grad_u: [v3(f.grad_u[0])?, v3(f.grad_u[1])?, v3(f.grad_u[2])?],
            grad_p: v3(f.grad_p)?,
            psi: v3(f.psi)?,
            p,
            mass_flux: self.dev.download(f.mass_flux)?,
            pb,
            scalars: match &self.scalars {
                Some(sc) => sc.states(self.dev.as_mut())?,
                None => Vec::new(),
            },
            periodic_force: self.periodic_force(),
            turbulence,
        })
    }

    /// Replaces the fields by `s` (a restart).
    pub fn load(&mut self, s: &State) -> Result<()> {
        let f = self.f;
        let mut put = |b: Buf3, v: &Soa3| -> Result<()> {
            self.dev.upload(b.x, &v.x)?;
            self.dev.upload(b.y, &v.y)?;
            self.dev.upload(b.z, &v.z)
        };
        put(f.u, &s.u)?;
        put(f.ub, &s.ub)?;
        for i in 0..3 {
            put(f.grad_u[i], &s.grad_u[i])?;
        }
        put(f.grad_p, &s.grad_p)?;
        put(f.psi, &s.psi)?;
        self.dev.upload(f.mass_flux, &s.mass_flux)?;
        match (&s.turbulence, self.turbulence.as_mut()) {
            (Some(ts), Some(t)) => {
                let (p, pb) = crate::turbulence::Turbulence::pressure_of(ts);
                self.dev.upload(f.p, p)?;
                self.dev.upload(f.pb, pb)?;
                t.load(self.dev.as_mut(), ts)?;
                // p̂ of the pressure boundaries from the loaded k (the
                // start set them from the start fields)
                t.modified_pressure(self.dev.as_mut(), &self.bc)?;
            }
            _ => {
                self.dev.upload(f.p, &s.p)?;
                self.dev.upload(f.pb, &s.pb)?;
            }
        }
        if let Some(sc) = self.scalars.as_mut() {
            sc.load(self.dev.as_mut(), &s.scalars)?;
        }
        if s.periodic_force.len() == self.periodic.len() {
            for (d, &b) in self.periodic.iter_mut().zip(&s.periodic_force) {
                d.beta = b;
            }
        }
        self.hold_convective(f.ub)
    }

    /// The boundary values, gradients and limiters as they are.
    pub fn lagged(&mut self) -> Result<Lagged> {
        let f = self.f;
        let mut v3 = |b: Buf3| -> Result<Soa3> {
            Ok(Soa3 {
                x: self.dev.download(b.x)?,
                y: self.dev.download(b.y)?,
                z: self.dev.download(b.z)?,
            })
        };
        Ok(Lagged {
            ub: v3(f.ub)?,
            grad_u: [v3(f.grad_u[0])?, v3(f.grad_u[1])?, v3(f.grad_u[2])?],
            grad_p: v3(f.grad_p)?,
            psi: v3(f.psi)?,
            pb: self.dev.download(f.pb)?,
        })
    }

    /// Puts the boundary values, gradients and limiters back. An observer
    /// between two outer iterations (an intermediate result, a balance
    /// table) brings them up to date for its output
    /// ([`update_boundaries`](Self::update_boundaries),
    /// [`gradients`](Self::gradients)). Without this the next iteration
    /// would start from values that have seen one more update, and the run
    /// would depend on its output schedule (D31).
    pub fn set_lagged(&mut self, l: &Lagged) -> Result<()> {
        let f = self.f;
        let mut put = |b: Buf3, v: &Soa3| -> Result<()> {
            self.dev.upload(b.x, &v.x)?;
            self.dev.upload(b.y, &v.y)?;
            self.dev.upload(b.z, &v.z)
        };
        put(f.ub, &l.ub)?;
        for (g, v) in f.grad_u.into_iter().zip(&l.grad_u) {
            put(g, v)?;
        }
        put(f.grad_p, &l.grad_p)?;
        put(f.psi, &l.psi)?;
        self.dev.upload(f.pb, &l.pb)
    }

    /// The convective outlets hold the boundary velocity `ub` (β = 0,
    /// c = `ub`), so that boundary updates reproduce it bitwise.
    pub(crate) fn hold_convective(&mut self, ub: Buf3) -> Result<()> {
        let Some(cv) = self.conv else { return Ok(()) };
        for (a, b) in ub.all().into_iter().zip(cv.c.all()) {
            self.dev.copy(a, b)?;
        }
        self.dev.fill(cv.beta, 0.0)?;
        self.conv = Some(Convective { hold: true, ..cv });
        Ok(())
    }

    /// The convective outlets of a time step from the old boundary levels:
    /// c = (−a₁ u_bⁿ − a₂ u_bⁿ⁻¹)/a₀ and β = U_c Δt/(a₀ Δn).
    pub(crate) fn step_convective(
        &mut self,
        ub_old: Buf3,
        ub_old2: Buf3,
        [a0, a1, a2]: [f64; 3],
        dt: f64,
    ) -> Result<()> {
        let Some(cv) = self.conv else { return Ok(()) };
        for ((c, a), b) in cv.c.all().into_iter().zip(ub_old.all()).zip(ub_old2.all()) {
            self.dev.lincomb(c, -a1 / a0, a, -a2 / a0, b)?;
        }
        let uc = cv.default_velocity(self.phys.rho);
        let s = dt / a0;
        self.dev
            .lincomb(cv.beta, s, cv.given_dn, s * uc, cv.default_dn)?;
        self.conv = Some(Convective { hold: false, ..cv });
        Ok(())
    }

    /// New prescribed boundary pressures (time-dependent conditions, D53);
    /// faces that do not fix the pressure ignore them.
    pub fn set_boundary_pressure(&mut self, pset: &[f64]) -> Result<()> {
        self.bc.pset.copy_from_slice(pset);
        self.dev.set_pressure_values(pset)
    }

    /// Backflow suppression after the predicted fluxes `m` (D29).
    pub(crate) fn close_backflow(&mut self, m: Buf) -> Result<()> {
        if self.no_backflow {
            self.dev.close_backflow(m)?;
        }
        Ok(())
    }

    /// New prescribed boundary velocity (time-dependent conditions); faces
    /// with another update rule ignore it.
    pub fn set_boundary_velocity(&mut self, ub: &[Vec<f64>; 3]) -> Result<()> {
        for (own, new) in self.bc.ub.iter_mut().zip(ub) {
            own.copy_from_slice(new);
        }
        self.pseudo_speed = None;
        let s = self.f.ub_set;
        self.dev.upload(s.x, &ub[0])?;
        self.dev.upload(s.y, &ub[1])?;
        self.dev.upload(s.z, &ub[2])?;
        if let Some(cv) = self.conv {
            let q = fixed_outflow(&self.disc, &self.bc, ub, &self.phys);
            self.conv = Some(Convective { q_other: q, ..cv });
        }
        Ok(())
    }

    pub(crate) fn flux_args(
        &self,
        d: Option<Buf>,
        relax: Buf,
        time: Option<(Buf, Buf, f64)>,
    ) -> FluxArgs {
        let f = &self.f;
        FluxArgs {
            u: f.u,
            grad_u: f.grad_u,
            p: f.p,
            grad_p: f.grad_p,
            pb: f.pb,
            ub: f.ub,
            d,
            relax,
            time,
            bounded: (self.settings.face_velocity == FaceVelocity::Bounded).then_some(f.psi),
        }
    }

    /// Step 1: boundary values from the current fields.
    pub fn update_boundaries(&mut self) -> Result<()> {
        let f = &self.f;
        self.dev.boundary_values(&BoundaryArgs {
            u: f.u,
            grad_u: f.grad_u,
            p: f.p,
            grad_p: f.grad_p,
            mass_flux: f.mass_flux,
            ub_set: f.ub_set,
            ub: f.ub,
            pb: f.pb,
        })?;
        // the mass balance of the convective outlets: one normal speed δ on
        // all of them, ρ Σ (u_b + δ n̂)·S = −(the other fixed fluxes)
        if let Some(cv) = self.conv
            && cv.balance
            && !cv.hold
            && cv.area_total > 0.0
        {
            let rho = self.phys.rho;
            let ub = self.f.ub;
            let mut q = 0.0;
            for (u, a) in ub.all().into_iter().zip(cv.area.all()) {
                q += self.dev.dot(u, a)?;
            }
            let delta = -(rho * q + cv.q_other) / (rho * cv.area_total);
            for (u, n) in ub.all().into_iter().zip(cv.normal.all()) {
                self.dev.lincomb(u, 1.0, u, delta, n)?;
            }
        }
        Ok(())
    }

    /// Step 2: gradients and limiters.
    pub fn gradients(&mut self) -> Result<()> {
        let f = self.f;
        for i in 0..3 {
            let (u, ub) = (f.u.all()[i], f.ub.all()[i]);
            self.dev
                .gradient(Lsq::Velocity, u, Mask::Velocity, ub, f.grad_u[i])?;
        }
        self.dev
            .gradient(Lsq::Pressure, f.p, Mask::Pressure, f.pb, f.grad_p)?;
        let s = &self.settings;
        let bounded = s.face_velocity == FaceVelocity::Bounded;
        if !bounded
            && (s.momentum.scheme == Scheme::Upwind || s.momentum.limiter == LimiterKind::None)
        {
            return Ok(());
        }
        // the bounded face velocity (D50) takes the limiter also where the
        // convection scheme has none
        let lim = match s.momentum.limiter {
            LimiterKind::None => LimiterKind::Venkatakrishnan,
            l => l,
        };
        let kk = s.momentum.limiter_k;
        for i in 0..3 {
            let (u, ub, psi) = (f.u.all()[i], f.ub.all()[i], f.psi.all()[i]);
            let eps2 = match lim {
                LimiterKind::VenkatakrishnanWang => {
                    let (mn, mx) = self.dev.min_max(u)?;
                    let e = kk * (mx - mn);
                    Eps2::Single(e * e)
                }
                // a negative single ε² selects Barth–Jespersen
                LimiterKind::BarthJespersen => Eps2::Single(-1.0),
                _ => Eps2::PerCell,
            };
            self.dev.limiter(u, ub, f.grad_u[i], eps2, psi)?;
        }
        Ok(())
    }

    /// Step 3: the relaxation term (1 − α)(ṁ − ρ ū·S) of the current fields
    /// into `w.relax` (and returns nothing else).
    pub(crate) fn relaxation_term(&mut self, alpha: f64) -> Result<()> {
        let a = self.flux_args(None, self.w.zero_face, None);
        self.dev.mass_flux(&a, self.w.relax)?;
        self.dev
            .scaled_difference(self.w.relax, 1.0 - alpha, self.f.mass_flux, self.w.relax)
    }

    /// Step 4 (assembly): momentum faces with the convecting flux `conv`
    /// and the rows; returns the residual scale.
    pub(crate) fn momentum_assembly(
        &mut self,
        conv: Buf,
        alpha: f64,
        theta: f64,
        time: Option<TimeTerms>,
    ) -> Result<f64> {
        self.momentum_assembly_faces(conv)?;
        let (f, w) = (self.f, self.w);
        self.dev.momentum_rows(&MomentumRowArgs {
            cp: w.cp,
            cn: w.cn,
            bf: w.bf,
            grad_p: f.grad_p,
            body_force: self.body_force(),
            u: f.u,
            alpha,
            theta,
            time,
            val: w.val,
            rhs: w.rhs,
            d: w.d,
            ds: w.ds,
            dc: w.dc,
            diag: w.diag,
        })?;
        if let Some(src) = self.source {
            for (r, q) in w.rhs.all().into_iter().zip(src.all()) {
                self.dev.lincomb(r, 1.0, r, 1.0, q)?;
            }
        }
        if let Some(src) = self.buoyancy_source {
            for (r, q) in w.rhs.all().into_iter().zip(src.all()) {
                self.dev.lincomb(r, 1.0, r, 1.0, q)?;
            }
        }
        self.add_coriolis()?;
        let alg = self.settings.algorithm(self.transient);
        if alg == Algorithm::Simple
            || (alg == Algorithm::Piso && !self.settings.simplec_corrector())
        {
            self.dev.copy(w.d, w.dc)?;
        }
        // residual scale ‖a_P u‖ + ‖b‖ over the three components (unrelaxed
        // a_P and b), so that it is finite also from a fluid at rest
        self.dev.residual_scale(f.u, w.rhs, w.diag, alpha)
    }

    /// The momentum face coefficients of the current fields convected by
    /// `conv` into `w.cp`, `w.cn`, `w.bf`.
    pub(crate) fn momentum_assembly_faces(&mut self, conv: Buf) -> Result<()> {
        let s = &self.settings;
        let (f, w) = (self.f, self.w);
        let uc = self.central_u.unwrap_or(f.u);
        self.dev.momentum_faces(&MomentumFaceArgs {
            mass_flux: conv,
            u: f.u,
            uc,
            grad_u: f.grad_u,
            psi: f.psi,
            ub: f.ub,
            blend: if s.momentum.scheme == Scheme::Upwind {
                0.0
            } else {
                1.0
            },
            central: if s.momentum.scheme == Scheme::Central {
                s.momentum.central_blending
            } else {
                0.0
            },
            beta: s.momentum.central_bounding,
            nonorth: if s.nonorth_correction { 1.0 } else { 0.0 },
            convection: if self.phys.convection { 1.0 } else { 0.0 },
            // the variable density (stage 7c): the full stress μ (∇u +
            // ∇uᵀ − ⅔ (∇·u) I), the transposed part with ⅔ of the trace
            transposed: if self.turbulence.is_some() || self.density.is_some() {
                1.0
            } else {
                0.0
            },
            trace_free: if self.density.is_some() {
                2.0 / 3.0
            } else if s.transposed_stress == TransposedStress::TraceFree {
                1.0
            } else {
                0.0
            },
            cp: w.cp,
            cn: w.cn,
            bf: w.bf,
        })?;
        match self.turbulence.as_ref() {
            Some(t) => t.add_explicit_stress(&mut *self.dev, f.grad_u, f.u, f.ub, w.bf),
            None => Ok(()),
        }
    }

    /// Step 4 (solve): the three momentum components; residuals over
    /// `scale` and iterations into `res`.
    pub(crate) fn momentum_solve(&mut self, scale: f64, res: &mut Residuals) -> Result<()> {
        let (f, w) = (self.f, self.w);
        let sym = self.bc.mom_kind.contains(&k::MOM_SYMMETRY);
        let s = &self.settings;
        let (tol, max) = (s.momentum_tol(self.transient), s.momentum.max_iter);
        for i in 0..3 {
            if sym {
                // the diagonal of this component (symmetry planes)
                self.dev.set_diagonal(w.val, w.diag.all()[i])?;
            }
            let (it, r0) = self.dev.solve(
                w.val,
                w.rhs.all()[i],
                f.u.all()[i],
                Preconditioner::Jacobi,
                false,
                tol,
                max,
            )?;
            res.momentum[i] = r0 / scale;
            res.momentum_iterations[i] = it;
        }
        Ok(())
    }

    /// Continuity residual of the fluxes `m`.
    pub(crate) fn continuity(&mut self, m: Buf) -> Result<f64> {
        let w = self.w;
        self.dev.net_outflow(m, w.defect)?;
        if let (Some(src), Some(_)) = (self.density_source, self.density_time) {
            // the mass source V ∂ρ/∂t of a variable density (stage 7c)
            self.dev.lincomb(w.defect, 1.0, w.defect, 1.0, src)?;
        }
        self.dev.throughflow(m, w.through)?;
        let d = self.dev.sum_sq(w.defect)?.sqrt();
        let t = self.dev.sum_sq(w.through)?.sqrt();
        Ok(d / t.max(f64::MIN_POSITIVE))
    }

    /// The multigrid setup on the pressure matrix `w.pval` (the first one,
    /// and again when `rebuild`).
    pub(crate) fn pressure_amg(&mut self, rebuild: bool) -> Result<()> {
        let t = Instant::now();
        let w = self.w;
        let s = self.settings.clone();
        if s.pressure.preconditioner == Preconditioner::Amg && (self.amg_setups == 0 || rebuild) {
            if let Err(e) = self.dev.amg_setup(w.pval, &s.pressure.amg) {
                if self.amg_setups == 0 {
                    return Err(e);
                }
                // the matrix was symmetric positive definite at the earlier
                // setups: the fields that build it (the momentum diagonal
                // of SIMPLEC) have run away, a divergence (D63), not a
                // missing pressure level
                return Err(crate::diag::err(crate::diag::Diagnostic::new(
                    crate::diag::RUN_DIVERGED,
                    format!(
                        "the solution diverged in outer iteration {}: the pressure-correction \
                         matrix is no longer positive definite (the multigrid setup failed, \
                         after {} successful setups)",
                        self.iteration, self.amg_setups
                    ),
                )));
            }
            self.amg_matrix = Some(self.dev.download(w.pval)?);
            self.amg_setups += 1;
            self.timings.amg_setup += t.elapsed().as_secs_f64();
        }
        Ok(())
    }

    /// Step 6: the pressure correction of the fluxes `m` into `w.pc`
    /// (shifted to zero in the reference cell); `rebuild` sets the
    /// multigrid up first. Returns the CG iterations.
    pub(crate) fn pressure_correction(&mut self, m: Buf, rebuild: bool) -> Result<usize> {
        let t = Instant::now();
        let w = self.w;
        self.dev.pressure_faces(w.dc, w.cf)?;
        self.dev.pressure_rows(w.cf, m, w.pval, w.prhs)?;
        if let Some(r) = self.reference_cell {
            let pos = self.disc.sell.pos(r, 0);
            self.dev.scale_entry(w.pval, pos, 2.0)?;
        }
        self.dev.sync()?;
        self.timings.pressure_assembly += t.elapsed().as_secs_f64();
        self.pressure_amg(rebuild)?;
        let s = self.settings.clone();
        let t = Instant::now();
        self.dev.fill(w.pc, 0.0)?;
        let (it, _) = self.dev.solve(
            w.pval,
            w.prhs,
            w.pc,
            s.pressure.preconditioner,
            true,
            s.pressure_tol(self.transient),
            s.pressure.max_iter,
        )?;
        if let Some(r) = self.reference_cell {
            let shift = self.dev.value(w.pc, r)?;
            self.dev.subtract(w.pc, shift)?;
        }
        self.dev.sync()?;
        self.timings.pressure_solve += t.elapsed().as_secs_f64();
        Ok(it)
    }

    /// Step 7: corrected fluxes `m` → `w.mc`, u −= D^C ∇p', p += α_p p'.
    pub(crate) fn corrections(&mut self, m: Buf, alpha_p: f64) -> Result<()> {
        let (f, w) = (self.f, self.w);
        self.dev.correct_fluxes(w.cf, m, w.pc, w.mc)?;
        // gradient of p': zero on pressure boundaries, zero normal gradient
        // elsewhere
        self.dev.gradient(
            Lsq::Pressure,
            w.pc,
            Mask::PressureCorrection,
            w.zero_boundary,
            w.gpc,
        )?;
        self.dev.correct_velocity(w.dc, w.gpc, f.u)?;
        self.drive_mass_flow()?;
        self.dev.lincomb(f.p, 1.0, f.p, alpha_p, w.pc)
    }

    /// Rebuilds the AMG hierarchy of a restart from the stored pressure
    /// matrix, so that the run continues bit-identically.
    pub fn restore_amg(&mut self, m: &[f64], setups: usize) -> Result<()> {
        let w = self.w;
        self.dev.upload(w.pval, m)?;
        let amg = self.settings.pressure.amg.clone();
        self.dev.amg_setup(w.pval, &amg)?;
        self.amg_setups = setups;
        self.amg_matrix = Some(m.to_vec());
        Ok(())
    }

    /// Re-validated settings that may change during a run (D40).
    pub fn set_settings(&mut self, s: Settings) -> Result<()> {
        s.validate(self.transient)?;
        if s.coupled.multigrid != self.settings.coupled.multigrid {
            // hierarchies of other settings are not kept
            self.sweep = Default::default();
        }
        // the remediation of the start stays in force (D50)
        self.settings = s.remedied(self.remediation.max_skewness).0;
        Ok(())
    }

    /// One steady outer iteration; returns the residuals before the solves.
    pub fn iterate(&mut self) -> Result<Residuals> {
        if self.settings.algorithm(false) == Algorithm::Coupled {
            return self.iterate_coupled();
        }
        self.iteration += 1;
        let mut res = Residuals {
            iteration: self.iteration,
            ..Default::default()
        };
        let s = self.settings.clone();
        let alpha = s.alpha_u(false);
        let t = Instant::now();
        self.update_boundaries()?;
        self.gradients()?;
        self.dev.sync()?;
        self.timings.gradients += t.elapsed().as_secs_f64();

        // step 3: relaxation term of the Rhie–Chow flux; with buoyancy also
        // the pseudo-time term of the buoyancy step (stage 7b1), with the
        // current iterate as the old level: Y = ṁ − ρ ū·S into the scratch
        // `w.mc` (free until the corrections)
        let t = Instant::now();
        let step = self.buoyancy_step();
        res.pseudo_step = step;
        let w = self.w;
        let time = match step {
            Some(dt) => {
                let a = self.flux_args(None, w.zero_face, None);
                self.dev.mass_flux(&a, w.mc)?;
                self.dev
                    .scaled_difference(w.mc, 1.0, self.f.mass_flux, w.mc)?;
                self.dev.lincomb(w.relax, 1.0 - alpha, w.mc, 0.0, w.mc)?;
                Some(TimeTerms {
                    a: [1.0, 1.0, 0.0],
                    dt,
                    u_old: self.f.u,
                    u_old2: self.f.u,
                    explicit: (0.0, self.f.u),
                })
            }
            None => {
                self.relaxation_term(alpha)?;
                None
            }
        };
        // step 4: momentum
        let conv = self.convecting_flux()?;
        let scale = self.momentum_assembly(conv, alpha, 1.0, time)?;
        self.dev.sync()?;
        self.timings.momentum_assembly += t.elapsed().as_secs_f64();
        let t = Instant::now();
        self.momentum_solve(scale, &mut res)?;
        self.dev.sync()?;
        self.timings.momentum_solve += t.elapsed().as_secs_f64();

        // step 5: predicted fluxes
        let t = Instant::now();
        let a = self.flux_args(Some(w.d), w.relax, step.map(|_| (w.ds, w.mc, alpha)));
        self.dev.mass_flux(&a, self.w.m)?;
        self.close_backflow(self.w.m)?;
        res.continuity = self.continuity(self.w.m)?;
        self.dev.sync()?;
        self.timings.pressure_assembly += t.elapsed().as_secs_f64();

        // step 6: pressure correction
        let rebuild = (self.iteration - 1).is_multiple_of(s.amg_rebuild);
        res.pressure_iterations = self.pressure_correction(self.w.m, rebuild)?;

        // step 7: corrections
        let t = Instant::now();
        self.corrections(self.w.m, s.alpha_p(false))?;
        self.dev.sync()?;
        self.timings.pressure_assembly += t.elapsed().as_secs_f64();
        if !res.finite() {
            self.diverged = Some(res.clone());
            return Err(crate::diag::err(crate::diag::Diagnostic::new(
                crate::diag::RUN_DIVERGED,
                format!(
                    "the solution diverged in outer iteration {} (residuals {:?}, continuity {:.3e})",
                    self.iteration, res.momentum, res.continuity
                ),
            )));
        }
        std::mem::swap(&mut self.f.mass_flux, &mut self.w.mc);
        if let Some(t) = self.turbulence.as_mut() {
            t.amg_rebuild = rebuild;
        }
        self.turbulence_step(None)?;
        self.solve_scalars(None, &mut res)?;
        Ok(res)
    }

    /// Adds the turbulence model (stage 6): its start viscosity enters the
    /// next momentum assembly, its equations are solved after every
    /// pressure correction.
    pub fn add_turbulence(&mut self, mut t: crate::turbulence::Turbulence) -> Result<()> {
        let f = self.f;
        t.start(self.dev.as_mut(), f.u, f.ub, &self.bc)?;
        self.turbulence = Some(t);
        Ok(())
    }

    /// The turbulence model, if the case has one.
    pub fn turbulence(&self) -> Option<&crate::turbulence::Turbulence> {
        self.turbulence.as_ref()
    }

    /// The turbulence model, mutable (verification set-ups).
    pub fn turbulence_mut(
        &mut self,
    ) -> Option<(&mut crate::turbulence::Turbulence, &mut dyn FlowDevice)> {
        match self.turbulence.as_mut() {
            Some(t) => Some((t, self.dev.as_mut())),
            None => None,
        }
    }

    /// The force per unit volume of the explicit stress of the QCR
    /// (`turbulence.qcr`) on every cell, from the velocity gradient and
    /// the eddy viscosity as they are: the face sums that the momentum
    /// equations take (verification: manufactured stress, stage 6
    /// specification 9.2). Zero without the QCR.
    pub fn explicit_stress_force(&mut self) -> Result<Vec<[f64; 3]>> {
        let nf = self.disc.n_faces;
        let mut bf = Vec::new();
        for _ in 0..3 {
            let b = self.dev.alloc(nf)?;
            self.dev.fill(b, 0.0)?;
            bf.push(b);
        }
        let bf = Buf3 {
            x: bf[0],
            y: bf[1],
            z: bf[2],
        };
        if let Some(t) = self.turbulence.as_ref() {
            t.add_explicit_stress(&mut *self.dev, self.f.grad_u, self.f.u, self.f.ub, bf)?;
        }
        let mut v = Vec::new();
        for b in bf.all() {
            v.push(self.dev.download(b)?);
        }
        let d = &self.disc;
        let mut force = vec![[0.0; 3]; d.n_cells];
        for f in 0..nf {
            // b_f is subtracted from the owner row and added to the
            // neighbour row of the momentum equations: the force on the
            // owner is the flux through the face, ∮ τ·dS = −b_f
            for c in 0..3 {
                force[d.owner[f] as usize][c] -= v[c][f];
                if f < d.n_interior {
                    force[d.neighbour[f] as usize][c] += v[c][f];
                }
            }
        }
        for (c, fc) in force.iter_mut().enumerate() {
            for x in fc.iter_mut() {
                *x /= d.cell_volume[c];
            }
        }
        Ok(force)
    }

    /// One steady outer iteration of the turbulence equations alone, with
    /// the velocity and the fluxes as they are (verification: manufactured
    /// solutions with a prescribed velocity, stage 6 specification 9.2).
    pub fn turbulence_iteration(&mut self) -> Result<crate::turbulence::TurbResiduals> {
        let mut res = Residuals::default();
        self.solve_turbulence(None, &mut res)?;
        self.iteration += 1;
        Ok(res.turbulence)
    }

    /// The turbulence solves of an outer iteration (after the pressure
    /// correction), into `res`.
    pub(crate) fn solve_turbulence(
        &mut self,
        step: Option<crate::scalar::StepTime>,
        res: &mut Residuals,
    ) -> Result<()> {
        if self.turbulence.is_none() {
            return Ok(());
        }
        let conv = self.convecting_flux()?;
        let Some(t) = self.turbulence.as_mut() else {
            return Ok(());
        };
        let tm = Instant::now();
        let f = self.f;
        res.turbulence = t.solve(self.dev.as_mut(), f.u, f.ub, conv, f.grad_p, &self.bc, step)?;
        self.dev.sync()?;
        self.timings.turbulence += tm.elapsed().as_secs_f64();
        if !res.finite() {
            self.diverged = Some(res.clone());
            return Err(crate::diag::err(crate::diag::Diagnostic::new(
                crate::diag::RUN_DIVERGED,
                format!(
                    "the solution diverged in outer iteration {}: turbulence residuals {:?}",
                    self.iteration, res.turbulence.residual
                ),
            )));
        }
        Ok(())
    }

    /// Adds the passive scalars (D51); they are solved from the next outer
    /// iteration on.
    pub fn add_scalars(&mut self, s: crate::scalar::Scalars) {
        self.scalars = Some(s);
    }

    /// The passive scalars, if the case has any.
    pub fn scalars(&self) -> Option<&crate::scalar::Scalars> {
        self.scalars.as_ref()
    }

    /// The scalar solves of an outer iteration, with the current fluxes
    /// (after the pressure correction), into `res`; then the buoyancy of
    /// the new temperature.
    pub(crate) fn solve_scalars(
        &mut self,
        step: Option<crate::scalar::StepTime>,
        res: &mut Residuals,
    ) -> Result<()> {
        self.solve_turbulence(step, res)?;
        if self.cht.is_some() {
            let t = Instant::now();
            crate::cht::solve(self, step, res)?;
            self.timings.scalars += t.elapsed().as_secs_f64();
        } else {
            let coupling = self
                .turbulence
                .as_ref()
                .map(|t| t.scalar_coupling(self.dev.as_ref()));
            self.solve_scalar_fields(step, res, coupling)?;
        }
        self.refresh_density()?;
        self.refresh_buoyancy()
    }

    /// The fluid temperature per boundary face of the flow (`t` per cell):
    /// T_b of the energy equation, the face of the full mesh in a
    /// conjugate case (the fluid cell's temperature at an interface).
    pub(crate) fn fluid_boundary_temperature(&mut self, t: &[f64]) -> Result<Vec<f64>> {
        let ni = self.disc.n_interior;
        let nb = self.disc.n_faces - ni;
        let owner = |fb: usize| self.disc.owner[ni + fb] as usize;
        if let Some(c) = self.cht.as_deref_mut() {
            let full = match c.fe.scalars.as_ref() {
                Some(s) => s.energy_boundary_temperature(c.fe.dev.as_mut())?,
                None => None,
            };
            let nif = c.fe.disc.n_interior;
            return Ok((0..nb)
                .map(|fb| {
                    let g = c.view.face[ni + fb] as usize;
                    match &full {
                        Some(tb) if g >= nif => tb[g - nif],
                        _ => t[owner(fb)],
                    }
                })
                .collect());
        }
        let tb = match self.scalars.as_ref() {
            Some(s) => s.energy_boundary_temperature(self.dev.as_mut())?,
            None => None,
        };
        Ok(tb.unwrap_or_else(|| (0..nb).map(|fb| t[owner(fb)]).collect()))
    }

    /// Recomputes the variable density (stage 7c) from the current
    /// temperature: p₀ of a closed volume, ρ per cell and face on the
    /// device, the mass source and the old levels of a time step, and ρ,
    /// ∂ρ/∂t and the pressure work of the energy equation (nothing with a
    /// constant density).
    pub fn refresh_density(&mut self) -> Result<()> {
        if self.density.is_none() {
            return Ok(());
        }
        let Some(t) = self.fluid_temperature()? else {
            return Ok(());
        };
        let tb = self.fluid_boundary_temperature(&t)?;
        let d = self.density.as_deref_mut().expect("density");
        d.refresh(&t, &tb);
        self.dev.set_density(&d.rho, &d.rho_face)?;
        self.phys.rho_var = Some(std::sync::Arc::new((d.rho.clone(), d.rho_face.clone())));
        let (src, dp) = match self.density_time {
            Some((a, dt)) => (
                d.mass_source(a, dt),
                if d.pressure_work {
                    d.dp0_dt(a, dt)
                } else {
                    0.0
                },
            ),
            None => (None, 0.0),
        };
        let old = d.old.as_ref().filter(|_| self.density_time.is_some());
        self.dev.set_variable_density(
            src.as_deref(),
            old.map(|o| o[0].as_slice()),
            old.map(|o| o[1].as_slice()),
        )?;
        if let Some(s) = &src {
            let b = match self.density_source {
                Some(b) => b,
                None => {
                    let b = self.dev.alloc(self.disc.n_cells)?;
                    self.density_source = Some(b);
                    b
                }
            };
            self.dev.upload(b, s)?;
        }
        let vol = &self.disc.cell_volume;
        let drho_dt: Vec<f64> = src
            .map(|s| s.iter().zip(vol).map(|(m, v)| m / v).collect())
            .unwrap_or_default();
        let rho = d.rho.clone();
        if let Some(c) = self.cht.as_deref_mut() {
            let n = c.fe.disc.n_cells;
            let mut r = vec![f64::NAN; n];
            let mut dr = if drho_dt.is_empty() {
                Vec::new()
            } else {
                vec![0.0; n]
            };
            for (i, &g) in c.view.cell.iter().enumerate() {
                r[g as usize] = rho[i];
                if !drho_dt.is_empty() {
                    dr[g as usize] = drho_dt[i];
                }
            }
            if let Some(s) = c.fe.scalars.as_mut() {
                s.set_energy_density(Some(crate::scalar::VarRho {
                    rho: r,
                    drho_dt: dr,
                    pressure_work: dp,
                }));
            }
        } else if let Some(s) = self.scalars.as_mut() {
            s.set_energy_density(Some(crate::scalar::VarRho {
                rho,
                drho_dt,
                pressure_work: dp,
            }));
        }
        Ok(())
    }

    /// Starts a time step of a variable density (stage 7c) with the
    /// backward-difference coefficients `a` (a₀, −a₁, −a₂) and Δt: the
    /// current density becomes level n.
    pub(crate) fn density_step_start(&mut self, a: [f64; 3], dt: f64) -> Result<()> {
        let Some(d) = self.density.as_deref_mut() else {
            return Ok(());
        };
        d.start_step();
        self.density_time = Some((a, dt));
        self.refresh_density()
    }

    /// Takes the last time step of a variable density back (after the
    /// temperature was restored).
    pub(crate) fn density_step_back(&mut self) -> Result<()> {
        let Some(d) = self.density.as_deref_mut() else {
            return Ok(());
        };
        d.step_back();
        self.density_time = None;
        self.refresh_density()
    }

    /// Sets up the Boussinesq buoyancy with the constants `c` (stage 7b1)
    /// and computes its arrays from the current temperature.
    pub fn attach_buoyancy(&mut self, c: crate::buoyancy::Constants) -> Result<()> {
        let b =
            crate::buoyancy::Balanced::new(c, &self.disc, &self.lsq_p, &self.bc.p_zero_gradient);
        self.buoyancy = Some(Box::new(b));
        self.refresh_buoyancy()
    }

    /// The fluid temperature per cell of the flow's cells: that of its own
    /// energy equation, or in a conjugate case of the energy view's.
    pub fn fluid_temperature(&mut self) -> Result<Option<Vec<f64>>> {
        if let Some(c) = self.cht.as_deref_mut() {
            let Some(phi) = c.fe.scalars.as_ref().and_then(|s| s.energy_phi()) else {
                return Ok(None);
            };
            let full = c.fe.dev.download(phi)?;
            return Ok(Some(
                c.view.cell.iter().map(|&i| full[i as usize]).collect(),
            ));
        }
        match self.scalars.as_ref().and_then(|s| s.energy_phi()) {
            Some(phi) => Ok(Some(self.dev.download(phi)?)),
            None => Ok(None),
        }
    }

    /// The pseudo-time step of the momentum equations of a steady SIMPLEC
    /// or SIMPLE run with buoyancy (stage 7b1): `factor` × L / U_b, with
    /// U_b = √(|g| |β| ΔT L) the buoyancy velocity, which is `factor` / N
    /// with N = √(|g| |β| ΔT / L) the buoyancy frequency. The relaxation
    /// alone takes much longer steps on coarse meshes (α/(1 − α) ρ h²/μ in
    /// a fluid at rest), and the explicit buoyancy force then overshoots
    /// the internal waves of a stable stratification instead of letting
    /// them settle. None without buoyancy or with an isothermal fluid.
    pub(crate) fn buoyancy_step(&mut self) -> Option<f64> {
        let l = self.pseudo_length();
        let speed = self.buoyancy.as_ref()?.speed(l);
        (speed > 0.0 && speed.is_finite()).then(|| self.settings.pseudo_time.factor * l / speed)
    }

    /// Recomputes the balanced-force arrays of the buoyancy from the
    /// current temperature and uploads them (nothing without buoyancy).
    pub fn refresh_buoyancy(&mut self) -> Result<()> {
        if self.buoyancy.is_none() {
            return Ok(());
        }
        let Some(t) = self.fluid_temperature()? else {
            return Ok(());
        };
        let nc = self.disc.n_cells;
        let rho = self.density.as_ref().map(|d| d.rho.clone());
        let a = self
            .buoyancy
            .as_mut()
            .expect("buoyancy")
            .arrays(&t, rho.as_deref());
        let b = match self.buoyancy_source {
            Some(b) => b,
            None => {
                let b = Buf3 {
                    x: self.dev.alloc(nc)?,
                    y: self.dev.alloc(nc)?,
                    z: self.dev.alloc(nc)?,
                };
                self.buoyancy_source = Some(b);
                b
            }
        };
        for (bi, v) in b.all().into_iter().zip(&a.source) {
            self.dev.upload(bi, v)?;
        }
        let vol = &self.disc.cell_volume;
        self.phys.buoyancy = Some(std::sync::Arc::new(std::array::from_fn(|i| {
            a.source[i].iter().zip(vol).map(|(s, v)| s / v).collect()
        })));
        self.dev.set_balanced_force(&a.jump, &a.pb)
    }

    /// The scalar solves with the turbulent coupling `coupling` (that of
    /// the flow's own turbulence model, or, in the energy view of a
    /// conjugate case, the fluid view's mapped onto the full mesh).
    pub(crate) fn solve_scalar_fields(
        &mut self,
        step: Option<crate::scalar::StepTime>,
        res: &mut Residuals,
        coupling: Option<crate::scalar::TurbCoupling>,
    ) -> Result<()> {
        if self.scalars.is_none() {
            return Ok(());
        }
        let t = Instant::now();
        let cells = &self.disc.cell_centroid;
        let xs: Vec<[f64; 3]> = if step.is_some() {
            (0..self.disc.n_cells)
                .map(|c| [cells.x[c], cells.y[c], cells.z[c]])
                .collect()
        } else {
            Vec::new()
        };
        let flows = self.periodic_mass_flow()?;
        let conv = self.convecting_flux()?;
        let Some(sc) = self.scalars.as_mut() else {
            return Ok(());
        };
        let link = crate::scalar::EnergyLink {
            grad_u: self.f.grad_u,
            mu_face: self.dev.face_viscosity(),
        };
        res.scalars = sc.solve(
            self.dev.as_mut(),
            conv,
            &flows,
            step,
            &xs,
            coupling,
            Some(link),
        )?;
        // μ(T) of the energy equation: the host copy of the output (forces,
        // y⁺) follows the device's
        if let Some(mu) = sc.energy_viscosity() {
            self.phys.mu_face = Some(std::sync::Arc::new(mu.to_vec()));
        }
        self.dev.sync()?;
        self.timings.scalars += t.elapsed().as_secs_f64();
        if !res.finite() {
            self.diverged = Some(res.clone());
            return Err(crate::diag::err(crate::diag::Diagnostic::new(
                crate::diag::RUN_DIVERGED,
                format!(
                    "the solution diverged in outer iteration {}: scalar residuals {:?}",
                    self.iteration, res.scalars.residual
                ),
            )));
        }
        Ok(())
    }
}

/// Mass flux (kg/s) out of the domain through the boundary faces that fix
/// the flux with the prescribed velocity `ub`, without the convective
/// outlets (in face order).
fn fixed_outflow(
    disc: &Discretisation,
    bc: &BoundaryArrays,
    ub: &[Vec<f64>; 3],
    phys: &Physics,
) -> f64 {
    let ni = disc.n_interior;
    let mut q = 0.0;
    for (fb, (&kind, &rule)) in bc.flux_kind.iter().zip(&bc.ub_rule).enumerate() {
        if kind == k::FLUX_FIXED && rule == k::UB_FIXED {
            let f = fb + ni;
            q += phys.rho_face_at(f)
                * (ub[0][fb] * disc.face_area.x[f]
                    + ub[1][fb] * disc.face_area.y[f]
                    + ub[2][fb] * disc.face_area.z[f]);
        }
    }
    q
}

/// Cell velocities from face mass fluxes, u_P = Σ_f ±ṁ_f (x_f − x_P)/(ρ V_P):
/// exact for a uniform velocity, since Σ_f (x_f − x_P) S_fᵀ = V_P I for a
/// closed cell (the start of D52).
fn velocity_of_fluxes(mesh: &Mesh, geo: &Geometry, m: &[f64], phys: &Physics, nc: usize) -> Soa3 {
    // a periodic face counts for its owner only: its partner closes the
    // cell across (D29)
    let ni = mesh.periodic_start();
    let mut u = Soa3::zeros(nc);
    for (fi, &mf) in m.iter().enumerate() {
        let xf = geo.face_centroid[fi];
        let sides = [(mesh.owner[fi] as usize, 1.0), (usize::MAX, 0.0)];
        let sides = if fi < ni {
            [sides[0], (mesh.neighbour[fi] as usize, -1.0)]
        } else {
            sides
        };
        for (c, s) in sides {
            if c < nc {
                let xp = geo.cell_centroid[c];
                u.x[c] += s * mf * (xf[0] - xp[0]);
                u.y[c] += s * mf * (xf[1] - xp[1]);
                u.z[c] += s * mf * (xf[2] - xp[2]);
            }
        }
    }
    for c in 0..nc {
        let r = 1.0 / (phys.rho_at(c) * geo.cell_volume[c]);
        u.x[c] *= r;
        u.y[c] *= r;
        u.z[c] *= r;
    }
    u
}

/// Passes of the non-orthogonal correction of the potential-flow start
/// (D52) after the first solve.
const POTENTIAL_NONORTH_PASSES: usize = 4;
