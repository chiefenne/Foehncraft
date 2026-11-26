//! The integration tests of foehncraft-flow, built as one test program.
//!
//! Each file of this folder is a module of this program, not a program of
//! its own: one program per file compiled and linked the crate's generic
//! code again for each file, which made the test build the longest step of
//! CI. Run one file with a filter on its module name:
//!
//! ```sh
//! cargo test --release -p foehncraft-flow --test flow energy::
//! ```
//!
//! A new test file needs a `mod` line here (`tools/check/test_modules.py`
//! checks it in CI). Every file writes into its own temporary folder, so the
//! files can share the process.

mod adjoint;
mod bcs;
mod broken_inputs;
mod buoyancy;
mod cht;
mod cht_acceleration;
mod cht_cycle;
mod cht_frozen;
mod coupled;
mod curvature;
mod datums;
mod derived;
mod droplet_validation;
mod eb;
mod energy;
mod expressions;
mod fluids;
mod health;
mod immersed;
mod impinging_jet;
mod inlet_table;
mod keps;
mod log_detail;
mod mesh_motion;
mod mrf;
mod periodic;
mod porous;
mod porous_energy;
mod pressure_values;
mod ptx_entries;
mod quality;
mod ranks;
mod run_ranks;
mod scalars;
mod schedule;
mod sliding;
mod species;
mod spray_coupling;
mod spray_evaporation;
mod spray_steady;
mod stage4;
mod stage5;
mod start_values;
mod statistics;
mod steam;
mod symmetry_stress;
mod thermal;
mod transient;
mod turbulence;
mod turbulence_derived;
mod uncertainty;
mod variable_density;
mod vof;
mod wall_crossing;
mod zeta_f;
mod zeta_f0;
