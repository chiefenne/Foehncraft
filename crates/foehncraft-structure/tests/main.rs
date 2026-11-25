//! The integration tests of foehncraft-structure, built as one test
//! program.
//!
//! Each file of this folder is a module of this program, not a program of
//! its own (one program per file linked the crate again for each file).
//! Run one file with a filter on its module name:
//!
//! ```sh
//! cargo test --release -p foehncraft-structure --test structure solids::
//! ```
//!
//! A new test file needs a `mod` line here (`tools/check/test_modules.py`
//! checks it in CI).

mod beams;
mod contact;
mod dynamics;
mod shells;
mod solids;
mod thermal;
