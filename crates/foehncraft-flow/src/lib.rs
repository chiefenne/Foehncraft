//! The incompressible flow solver of Foehncraft: steady laminar flow with the
//! SIMPLEC algorithm (stage 4) and transient flow with a PISO-type
//! algorithm (stage 5), with the Rhie–Chow interpolation on polyhedral
//! meshes, on the CPU or an NVIDIA GPU.
//!
//! - [`bc`]: boundary conditions (D29) and the check of their combination
//!   with guidance (D30);
//! - [`simplec`]: the outer iterations (D19, D27);
//! - [`analytic`]: exact solutions and reference data (Couette, Poiseuille,
//!   Ghia et al.);
//! - [`case`]: the `foehn flow` case file with validation (D39), the run
//!   with the log and CSV output (D43) and the verification checks;
//! - [`device`]: the devices of the flow loop (CPU, CUDA);
//! - [`transient`]: the time step (D28): BDF1, BDF2, Crank–Nicolson,
//!   PISO correctors, outer iterations;
//! - [`unsteady`]: the transient case driver with the adaptive time step;
//! - [`restart`]: restart files (D42);
//! - [`checkpoint`]: binomial checkpointing of transient runs for the
//!   reverse sweep of an adjoint (D59);
//! - [`control`]: control files and re-reading the case file (D40);
//! - [`derived`]: derived fields such as vorticity, Q, λ₂ and the wall
//!   shear stress, at write time and after the run (D72);
//! - [`expr`]: expressions of space and time for boundary values, the body
//!   force and the viscosity, checked and compiled when the case is read
//!   (D53);
//! - [`materials`]: the material database, case materials and regions, and
//!   the resolved property program of density and viscosity (D94–D102,
//!   D112);
//! - [`monitor`]: probes, force coefficients, histories, benchmarks (D43);
//! - [`plot`]: small SVG charts of the reports.

pub mod acceleration;
#[cfg(test)]
mod adjoint_tests;
pub mod analytic;
pub mod balance;
pub mod banner_art;
pub mod bc;
pub mod budget;
pub mod buoyancy;
pub mod case;
pub mod checkpoint;
pub mod cht;
pub mod conduction;
pub mod console;
pub mod control;
pub mod coupled;
pub mod cycle;
pub mod density;
pub mod derived;
pub mod device;
pub mod diag;
pub mod divergence;
pub mod energy;
pub mod engine;
pub mod expr;
pub mod fluids;
pub mod frame;
pub mod frozen;
pub mod limitations;
pub mod limiter;
pub mod live;
pub mod logblocks;
pub mod logcols;
pub mod materials;
pub mod monitor;
pub mod operators;
pub mod plot;
pub mod pressure_values;
pub mod ranks;
pub mod residual;
pub mod restart;
pub mod rotor;
pub mod scalar;
pub mod schedule;
pub mod schema;
pub mod simplec;
pub mod sliding;
pub mod spray;
pub mod start;
pub mod statistics;
pub mod transient;
pub mod turbulence;
pub mod uncertainty;
pub mod unsteady;
