//! Devices that integrate the chemistry of the cells (SF2b): the GPU
//! backend registered by `foehncraft_cuda::register` (PTX entry
//! `chemistry_react`), and the CPU threads. Both run the one kernel body
//! [`react_cell`](foehncraft_kernels::chemistry::react_cell) per cell, with
//! the operations in the same order, so the GPU gives bitwise the states of
//! the CPU (D31).

use anyhow::Result;
use foehncraft_kernels::chemistry::{self as kc, Program, Settings, Stats};
use rayon::prelude::*;
use std::path::Path;
use std::sync::OnceLock;

/// A device that integrates the reactors of many cells over one time step.
pub trait Device: Send {
    /// Name of the device (`cuda`).
    fn name(&self) -> &str;

    /// Integrates every cell of `state` (species + 1 values per cell: the
    /// mass fractions, then the temperature) at the pressure `p` over `dt`
    /// (or over its own interval `dt_cell[c]`, the residence time of the
    /// fine structures of the eddy dissipation concept, SF3) with the
    /// mechanism the device was created with; writes the statistics of
    /// each cell into `stats`.
    #[allow(clippy::too_many_arguments)]
    fn react(
        &mut self,
        p: f64,
        dt: f64,
        dt_cell: Option<&[f64]>,
        min_temperature: f64,
        set: &Settings,
        state: &mut [f64],
        stats: &mut [Stats],
    ) -> Result<()>;
}

/// Builds a GPU device for a mechanism on device `ordinal` from the kernel
/// module (for CUDA the PTX file of `gpu/ptx-kernels/build.py`).
pub type GpuFactory = fn(&Program<'_>, usize, &Path) -> Result<Box<dyn Device>>;

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

/// Creates the GPU device `name` for the mechanism `prog`.
pub fn create_gpu(
    name: &str,
    prog: &Program<'_>,
    ordinal: usize,
    ptx: &Path,
) -> Result<Box<dyn Device>> {
    match GPU.get() {
        Some((n, f)) if *n == name => f(prog, ordinal, ptx),
        _ => anyhow::bail!(
            "this build has no `{name}` backend for the chemistry; build with `--features cuda`"
        ),
    }
}

/// The CPU threads: each cell on its own (any number of threads gives the
/// same bits); `dt_cell` as in [`Device::react`].
#[allow(clippy::too_many_arguments)]
pub fn react_cpu(
    prog: &Program<'_>,
    p: f64,
    dt: f64,
    dt_cell: Option<&[f64]>,
    min_temperature: f64,
    set: &Settings,
    state: &mut [f64],
    stats: &mut [Stats],
) {
    let n = prog.n();
    state
        .par_chunks_mut(n + 1)
        .zip(stats.par_iter_mut())
        .enumerate()
        .for_each_init(
            || vec![0.0; kc::work_len(n)],
            |work, (c, (z, st))| {
                let d = dt_cell.map_or(dt, |d| d[c]);
                *st = kc::react_cell(prog, p, z, d, set, min_temperature, work);
            },
        );
}
