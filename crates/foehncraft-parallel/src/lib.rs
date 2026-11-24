//! Parallel layer of Foehncraft.
//!
//! - [`Communicator`]: the interface between ranks. [`SerialComm`] is the
//!   one-rank implementation; `MpiComm` (feature `mpi`) uses the system MPI
//!   library (D6, D9). Domain decomposition is part of the data structures
//!   from the start (D6), so solver code is written against this trait even
//!   when it runs on one rank.
//! - [`ThreadComm`]: ranks as threads of one process, so multi-rank
//!   algorithms are tested bit for bit by `cargo test` without MPI.
//! - [`repro`]: multi-threaded reductions whose result does not depend on
//!   the number of threads (D31).

pub mod comm;
pub mod repro;
pub mod threads;

#[cfg(feature = "mpi")]
pub use comm::mpi_comm::MpiComm;
pub use comm::{Communicator, SerialComm};
pub use threads::ThreadComm;
