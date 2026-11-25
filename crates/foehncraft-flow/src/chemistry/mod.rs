//! Gas-phase chemistry (steam furnace stage SF2; the gas-phase slice of
//! the chemistry program D120).
//!
//! - [`mechanism`]: the built-in mechanisms (GRI-Mech 3.0, the global
//!   mechanisms of Westbrook & Dryer 1981 and Jones & Lindstedt 1988) and
//!   the importer of the YAML mechanism format of D103;
//! - [`equilibrium`]: chemical equilibrium by Gibbs minimisation at fixed
//!   T and p or h and p;
//! - [`reactor`]: the 0D constant-pressure adiabatic reactor, run by the
//!   tool `chemistry.reactor` (`foehn chemistry reactor`).
//!
//! The rates, their analytic Jacobian and the stiff integrator are the
//! kernels of [`foehncraft_kernels::chemistry`], the same code for the
//! cells of a flow (SF2b) and the GPU.

pub mod equilibrium;
pub mod mechanism;
pub mod reactor;

pub use mechanism::Mechanism;

#[cfg(test)]
mod tests;
