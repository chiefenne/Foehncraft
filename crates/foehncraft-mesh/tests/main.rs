//! The integration tests of foehncraft-mesh, built as one test program.
//!
//! Each file of this folder is a module of this program, not a program of
//! its own: one program per file compiled and linked the crate's generic
//! code again for each file, which made the test build the longest step of
//! CI. Run one file with a filter on its module name:
//!
//! ```sh
//! cargo test --release -p foehncraft-mesh --test mesh gmsh::
//! ```
//!
//! A new test file needs a `mod` line here (`tools/check/test_modules.py`
//! checks it in CI). Every file writes into its own temporary folder, so the
//! files can share the process.

mod broken_inputs;
mod canonical_local_mesh;
mod cgns;
mod cut;
mod decompose;
mod exact_reductions;
mod fluent_msh;
mod generated;
mod gmsh;
mod gmsh_fingerprint;
mod golden;
mod locate;
mod multizone;
mod openfoam;
mod ray;
mod su2;
mod vtk;
mod vtkhdf;
mod vtkhdf_read;
mod walk;
