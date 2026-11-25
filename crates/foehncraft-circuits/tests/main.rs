//! The integration tests of foehncraft-circuits, built as one test
//! program.
//!
//! Each file of this folder is a module of this program. Run one file
//! with a filter on its module name:
//!
//! ```sh
//! cargo test --release -p foehncraft-circuits --test circuits exact::
//! ```
//!
//! A new test file needs a `mod` line here (`tools/check/test_modules.py`
//! checks it in CI).

mod boiling;
mod exact;
