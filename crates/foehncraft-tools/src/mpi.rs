//! The MPI processes of a run (MPI step W7): a build with the `mpi`
//! feature started under an MPI launcher (`mpirun -n 4 foehn flow
//! case.toml`) runs `flow.run` on its processes, rank 0 writing the
//! output; without the feature, or with one process, nothing changes.

use foehncraft_flow::ranks::RunComm;
use foehncraft_mesh::decompose::SharedComm;
use std::sync::Mutex;

static WORLD: Mutex<Option<SharedComm>> = Mutex::new(None);

/// Holds MPI while it lives: dropping it finalises MPI (drop it before
/// `std::process::exit`).
#[must_use]
pub struct Guard(());

impl Drop for Guard {
    fn drop(&mut self) {
        let w = WORLD.lock().map(|mut w| w.take());
        drop(w);
    }
}

/// Whether an MPI launcher started this process (Open MPI, MPICH, PMIx).
pub fn launched() -> bool {
    [
        "OMPI_COMM_WORLD_SIZE",
        "PMI_SIZE",
        "PMIX_RANK",
        "MPI_LOCALNRANKS",
    ]
    .iter()
    .any(|k| std::env::var_os(k).is_some())
}

/// Starts MPI when this process runs under an MPI launcher with more than
/// one process and the build has the `mpi` feature.
pub fn start() -> Guard {
    #[cfg(feature = "mpi")]
    if launched()
        && let Some(c) = foehncraft_parallel::MpiComm::init()
        && foehncraft_parallel::Communicator::size(&c) > 1
        && let Ok(mut w) = WORLD.lock()
    {
        *w = Some(std::sync::Arc::new(c));
    }
    Guard(())
}

/// The communicator of the MPI processes, when there are several.
pub fn world() -> Option<RunComm> {
    WORLD.lock().ok().and_then(|w| w.clone()).map(RunComm)
}

/// The rank of this process (0 without MPI).
pub fn rank() -> usize {
    world().map_or(0, |c| c.0.rank())
}
