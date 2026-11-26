//! The integration tests of foehncraft-flow, first of two test programs
//! (`flow`, this file, and `flow2`, `tests/main2.rs`).
//!
//! Each file of this folder is a module of one of the two programs, not a
//! program of its own: one program per file compiled and linked the crate's
//! generic code again for each file, which made the test build the longest
//! step of CI. One program for all files compiled for nine minutes, because
//! the compiler front end works on one crate with one thread; two halves
//! compile side by side. Run one file with a filter on its module name, in
//! the program that lists it (`--test flow2` for the modules of
//! `tests/main2.rs`):
//!
//! ```sh
//! cargo test --release -p foehncraft-flow --test flow energy::
//! ```
//!
//! A new test file needs a `mod` line in one of the two files, in the
//! shorter one (`tools/check/test_modules.py` checks it in CI). Every file
//! writes into its own temporary folder, so the files can share the process.

mod adjoint;
mod bcs;
mod broken_inputs;
mod buoyancy;
mod cht;
mod curvature;
mod derived;
mod energy;
mod fluids;
mod health;
mod immersed;
mod inlet_table;
mod mesh_motion;
mod porous_energy;
mod pressure_values;
mod ptx_entries;
mod ranks;
mod scalars;
mod spray_coupling;
mod spray_evaporation;
mod spray_steady;
mod stage4;
mod start_values;
mod symmetry_stress;
mod transient;
mod turbulence;
mod variable_density;
mod vof;
mod zeta_f0;
