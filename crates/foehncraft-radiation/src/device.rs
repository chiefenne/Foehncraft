//! Devices that move the bundles (D259): the GPU backend registered by
//! `foehncraft_cuda::register`, and the same data path on the CPU threads.
//!
//! A device receives the bundles as flat records
//! ([`Bundle::to_words`](foehncraft_kernels::radiation::Bundle::to_words)),
//! runs the kernel body of `foehncraft_kernels::radiation` on each and
//! returns the integer tallies per batch. Integer sums do not depend on
//! the order of the additions, so every device gives bitwise the tallies of
//! the CPU threads.

use foehncraft_kernels::radiation::{BUNDLE_WORDS, Medium, Tally, transport_words};
use foehncraft_kernels::walk::WalkMesh;
use rayon::prelude::*;
use std::path::Path;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicI64, Ordering};

/// The tallies of one launch: quanta per batch and cell, per batch and
/// boundary face (arriving and absorbed), batch-major, and the escaped,
/// roulette and lost quanta.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct DeviceTallies {
    /// Absorbed per batch and cell (`batch × n_cells + cell`).
    pub cell: Vec<i64>,
    /// Arriving per batch and boundary face.
    pub incident: Vec<i64>,
    /// Absorbed per batch and boundary face.
    pub absorbed: Vec<i64>,
    /// Escaped, roulette, lost.
    pub ledger: [i64; 3],
}

/// A device that moves bundles to their end.
pub trait Transport: Send {
    /// Name of the device (`cuda`, `cpu-device-path`).
    fn name(&self) -> &str;

    /// Moves every bundle of `words` ([`BUNDLE_WORDS`] words each) on the
    /// mesh `m` (a GPU holds its own copy, uploaded when it was created)
    /// through the medium `md`, and returns the tallies of `batches`
    /// batches.
    fn run(
        &mut self,
        m: &WalkMesh<'_>,
        md: &Medium<'_>,
        words: &[u64],
        call: u64,
        batches: usize,
    ) -> anyhow::Result<DeviceTallies>;
}

/// Builds a GPU device for a walk mesh on device `ordinal` from the kernel
/// module (for CUDA the PTX file of `gpu/ptx-kernels/build.py`).
pub type GpuFactory = fn(&WalkMesh<'_>, usize, &Path) -> anyhow::Result<Box<dyn Transport>>;

static GPU: OnceLock<(&'static str, GpuFactory)> = OnceLock::new();

/// Registers the GPU backend of this build under `name` (`cuda`); called
/// once by `foehncraft_cuda::register`. Later calls are ignored.
pub fn register_gpu(name: &'static str, f: GpuFactory) {
    let _ = GPU.set((name, f));
}

/// The name of the registered GPU backend, if any.
pub fn gpu_backend() -> Option<&'static str> {
    GPU.get().map(|g| g.0)
}

pub(crate) fn create_gpu(
    m: &WalkMesh<'_>,
    ordinal: usize,
    ptx: &Path,
) -> anyhow::Result<Box<dyn Transport>> {
    match GPU.get() {
        Some((_, f)) => f(m, ordinal, ptx),
        None => anyhow::bail!("this build has no GPU backend; build with `--features cuda`"),
    }
}

/// The device data path on the CPU threads: flat records in, the kernel
/// body of the PTX entry per record, integer tallies out. It checks the
/// path of the GPU without a GPU.
#[derive(Debug, Default)]
pub struct CpuDevicePath;

struct Atomics<'a> {
    cell: &'a [AtomicI64],
    incident: &'a [AtomicI64],
    absorbed: &'a [AtomicI64],
    ledger: &'a [AtomicI64; 3],
    nc: usize,
    nb: usize,
}

impl Tally for Atomics<'_> {
    fn cell(&self, batch: u32, cell: usize, q: i64) {
        self.cell[batch as usize * self.nc + cell].fetch_add(q, Ordering::Relaxed);
    }

    fn wall(&self, batch: u32, face: usize, incident: i64, absorbed: i64) {
        let i = batch as usize * self.nb + face;
        self.incident[i].fetch_add(incident, Ordering::Relaxed);
        self.absorbed[i].fetch_add(absorbed, Ordering::Relaxed);
    }

    fn escaped(&self, q: i64) {
        self.ledger[0].fetch_add(q, Ordering::Relaxed);
    }

    fn roulette(&self, q: i64) {
        self.ledger[1].fetch_add(q, Ordering::Relaxed);
    }

    fn lost(&self, q: i64) {
        self.ledger[2].fetch_add(q, Ordering::Relaxed);
    }
}

impl Transport for CpuDevicePath {
    fn name(&self) -> &str {
        "cpu-device-path"
    }

    fn run(
        &mut self,
        m: &WalkMesh<'_>,
        md: &Medium<'_>,
        words: &[u64],
        call: u64,
        batches: usize,
    ) -> anyhow::Result<DeviceTallies> {
        anyhow::ensure!(
            words.len().is_multiple_of(BUNDLE_WORDS),
            "incomplete bundle records"
        );
        let (nc, nb) = (m.n_cells, md.normal.len());
        let zeros = |n: usize| -> Vec<AtomicI64> { (0..n).map(|_| AtomicI64::new(0)).collect() };
        let (cell, incident, absorbed) = (
            zeros(batches * nc),
            zeros(batches * nb),
            zeros(batches * nb),
        );
        let ledger = [AtomicI64::new(0), AtomicI64::new(0), AtomicI64::new(0)];
        let t = Atomics {
            cell: &cell,
            incident: &incident,
            absorbed: &absorbed,
            ledger: &ledger,
            nc,
            nb,
        };
        (0..words.len() / BUNDLE_WORDS)
            .into_par_iter()
            .for_each(|i| {
                transport_words(m, md, words, i, call, &t);
            });
        let take = |v: Vec<AtomicI64>| v.into_iter().map(AtomicI64::into_inner).collect();
        Ok(DeviceTallies {
            cell: take(cell),
            incident: take(incident),
            absorbed: take(absorbed),
            ledger: ledger.map(AtomicI64::into_inner),
        })
    }
}
