//! Aggregation-based algebraic multigrid (stage 3, D35, D37): the
//! preconditioner of CG for the pressure equation.
//!
//! **Setup** (host, float64, deterministic; [`setup`]): smoothed aggregation
//! after Vaněk, Mandel and Brezina (1996):
//!
//! 1. strength of connection, j strong for i when
//!    |a_ij| ≥ θ_l √(|a_ii a_jj|), θ_l = θ · 2⁻ˡ on level l;
//! 2. three-phase greedy aggregation in the stored row order (the cell
//!    order of the mesh): new aggregates from untouched neighbourhoods,
//!    remaining nodes join the aggregate of their strongest aggregated
//!    neighbour, the rest form aggregates of their own; nodes without strong
//!    connections stay out of every aggregate and are only smoothed;
//! 3. tentative prolongator: the constant vector on each aggregate,
//!    normalised;
//! 4. smoothed prolongator P = (I − ω D_F⁻¹ A_F) P_tent with the filtered
//!    matrix A_F (weak connections lumped onto the diagonal) and
//!    ω = 4 / (3 λ_max(D_F⁻¹ A_F)); plain (unsmoothed) aggregation skips this
//!    step;
//! 5. Galerkin coarse operator A_c = Pᵀ A P, restriction R = Pᵀ stored
//!    explicitly (so the restriction is a gather, no atomics);
//! 6. down to at most `max_coarse` unknowns, where a dense inverse (by
//!    Cholesky factorisation) solves the coarsest system.
//!
//! **Cycle** ([`cycle`]): a symmetric V-cycle with the same polynomial
//! smoother before and after the coarse correction: Chebyshev of degree k
//! on D⁻¹A over [a λ_max, b λ_max] (Adams, Brezina, Hu, Tuminaro 2003) or k
//! sweeps of damped Jacobi. λ_max comes from a fixed power iteration. A
//! symmetric V-cycle with polynomial smoothers is a symmetric positive
//! definite preconditioner, as CG requires.
//!
//! **Precision** (D37): in the default mixed mode the hierarchy (including a
//! float32 copy of the fine matrix) is stored and applied in float32; CG,
//! its residuals and the solution stay in float64. The float64 mode runs the
//! identical hierarchy in float64.
//!
//! [`reduction`] holds a second setup for the non-symmetric matrices of
//! transport equations (D179): a reduction multigrid with a splitting into
//! coarse and fine points. Its hierarchies run through the same cycle and
//! the same level operations.

pub mod csr;
pub mod cycle;
pub mod dist;
pub mod reduction;
pub mod setup;

pub use csr::{Csr, SellMatrix};
pub use reduction::ReductionSettings;
pub use setup::{Hierarchy, HierarchyInfo, build};

use serde::{Deserialize, Serialize};

/// Aggregation flavour.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Aggregation {
    /// Smoothed aggregation (Vaněk, Mandel, Brezina 1996).
    #[default]
    Smoothed,
    /// Plain (unsmoothed) aggregation: P = P_tent.
    Plain,
}

/// Smoother of the V-cycle.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Smoother {
    /// Chebyshev polynomial in D⁻¹A.
    #[default]
    Chebyshev,
    /// Damped Jacobi.
    Jacobi,
}

/// Precision of the hierarchy.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Precision {
    /// Hierarchy in float32, Krylov iteration in float64 (D37).
    #[default]
    Mixed,
    /// Everything in float64.
    Double,
}

impl Precision {
    /// Name for reports.
    pub fn name(self) -> &'static str {
        match self {
            Precision::Mixed => "mixed (f32 hierarchy)",
            Precision::Double => "double",
        }
    }
}

/// `[solver.amg]` settings.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct AmgSettings {
    /// Smoothed or plain aggregation.
    pub aggregation: Aggregation,
    /// Strength threshold θ on the finest level (halved on each coarser
    /// level).
    pub strength: f64,
    /// Smoother.
    pub smoother: Smoother,
    /// Chebyshev degree, or Jacobi sweeps, before and after the coarse
    /// correction.
    pub degree: usize,
    /// Chebyshev interval as fractions [a, b] of the estimated λ_max.
    pub chebyshev_interval: [f64; 2],
    /// Damping of the Jacobi smoother, as a multiple of 1/λ_max.
    pub jacobi_weight: f64,
    /// Damping of the prolongator smoothing, as a multiple of 1/λ_max.
    pub prolongation_weight: f64,
    /// Largest coarsest level (solved by a dense inverse).
    pub max_coarse: usize,
    /// Largest number of levels.
    pub max_levels: usize,
    /// Power iterations for λ_max.
    pub power_iterations: usize,
    /// Precision of the hierarchy.
    pub precision: Precision,
    /// Tests the pivots of the dense coarsest level against their own
    /// diagonals instead of the largest one: for the transport and
    /// elliptic equations of the turbulence models, whose coefficients
    /// span many decades on a fine mesh. Not a case key; the pressure keeps
    /// the strict test, which finds a matrix that stops being positive
    /// definite early.
    #[serde(skip)]
    #[schemars(skip)]
    pub scaled_pivots: bool,
}

impl Default for AmgSettings {
    fn default() -> Self {
        Self {
            aggregation: Aggregation::Smoothed,
            strength: 0.04,
            smoother: Smoother::Chebyshev,
            degree: 2,
            chebyshev_interval: [0.1, 1.1],
            jacobi_weight: 4.0 / 3.0,
            prolongation_weight: 4.0 / 3.0,
            max_coarse: 256,
            max_levels: 20,
            power_iterations: 20,
            precision: Precision::Mixed,
            scaled_pivots: false,
        }
    }
}

impl AmgSettings {
    /// Checks the ranges.
    pub fn validate(&self) -> anyhow::Result<()> {
        anyhow::ensure!(
            (0.0..1.0).contains(&self.strength),
            "[amg] strength must lie in [0, 1)"
        );
        anyhow::ensure!(
            (1..=8).contains(&self.degree),
            "[amg] degree must lie in 1..=8"
        );
        let [a, b] = self.chebyshev_interval;
        anyhow::ensure!(
            a > 0.0 && a < b && b >= 1.0,
            "[amg] chebyshev_interval must be [a, b] with 0 < a < b and b >= 1"
        );
        anyhow::ensure!(
            self.jacobi_weight > 0.0 && self.jacobi_weight < 2.0,
            "[amg] jacobi_weight must lie in (0, 2)"
        );
        anyhow::ensure!(
            self.prolongation_weight > 0.0 && self.prolongation_weight < 2.0,
            "[amg] prolongation_weight must lie in (0, 2)"
        );
        anyhow::ensure!(
            (1..=4000).contains(&self.max_coarse),
            "[amg] max_coarse must lie in 1..=4000 (dense inverse)"
        );
        anyhow::ensure!(self.max_levels >= 1, "[amg] max_levels must be >= 1");
        anyhow::ensure!(
            self.power_iterations >= 1,
            "[amg] power_iterations must be >= 1"
        );
        Ok(())
    }
}
