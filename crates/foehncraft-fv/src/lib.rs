//! Finite-volume operators on a scalar transport equation (stage 2) and
//! the linear solvers of the pressure equation (stage 3).
//!
//! - [`problem`]: the steady convection–diffusion problem, boundary
//!   conditions and manufactured solutions;
//! - [`disc`]: the host-side discretisation (geometry factors,
//!   least-squares weights, sliced-ELL structure);
//! - [`backend`]: where the kernel bodies of `foehncraft_kernels` run: CPU
//!   threads or the NVIDIA GPU (crate `foehncraft-cuda`, feature `cuda`);
//! - [`amg`]: smoothed-aggregation algebraic multigrid, the pressure
//!   preconditioner (stage 3);
//! - [`solver`]: outer deferred-correction iterations, BiCGStab and CG;
//! - [`poisson`]: the `foehn poisson` linear-solver study (stage 3);
//! - `hypre` (feature `hypre`): hypre's BoomerAMG as the reference solver;
//! - [`ranks`]: scalar solves on the ranks of a decomposed mesh and their
//!   bitwise check against one rank (MPI step W3);
//! - [`verify`]: gradient exactness, limiter bounds and order-of-accuracy
//!   sweeps (the acceptance checks of stage 2);
//! - [`case`]: the `foehn scalar` case file and run.

pub mod amg;
pub mod backend;
pub mod case;
pub mod disc;
pub mod exactness;
#[cfg(feature = "hypre")]
pub mod hypre;
pub mod poisson;
pub mod problem;
pub mod ranks;
pub mod solver;
pub mod verify;

pub use disc::Discretisation;
pub use problem::Problem;
