//! Core types shared by every Foehncraft crate.
//!
//! The crate is `no_std`, so that the same definitions can be compiled for
//! the CPU and for GPU targets (see `gpu/ptx-kernels`).
//!
//! - [`Real`]: the floating-point type is a generic parameter of the code
//!   (requirement D7). The reference path is `f64`.
//! - [`date`]: calendar dates from Unix time (D84).
//! - [`dual`]: forward-mode dual numbers, a [`Real`] whose kernels return
//!   exact derivatives next to their results (requirement D59).
//! - [`exact`]: exact, order-independent sums, dot products and min/max
//!   (D31 (1), D93 (4)): the same bits for any thread and rank count.
//! - [`index`]: local (per rank) and global index types.
//! - [`math`]: elementary functions written in Foehncraft, the same bits on
//!   every CPU and GPU (D93 (6)), for expressions and properties.
//! - [`repro`]: reductions with a fixed summation order (requirement D31).
//! - [`rng`]: counter-based random numbers keyed by identity, step and event,
//!   independent of thread or rank order (requirement D47).
#![no_std]

pub mod date;
pub mod dual;
pub mod exact;
pub mod index;
pub mod math;
pub mod real;
pub mod repro;
pub mod rng;

pub use dual::Dual;
pub use index::{GlobalIndex, LocalIndex};
pub use real::Real;

/// Version of the Foehncraft crates, from the workspace `Cargo.toml`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
