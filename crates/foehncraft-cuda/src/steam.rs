//! Water and steam states of IAPWS-IF97 (stage SF7a) on an NVIDIA GPU: the
//! entry `steam_state` of `gpu/ptx-kernels`, one GPU thread per state, the
//! kernel body that the CPU runs (`foehncraft_kernels::steam::state_entry`).
//! No reductions and no atomics, and the math routines are the in-house
//! ones of `foehncraft_core::math`, so the outputs are bitwise those of
//! `foehncraft_kernels::steam::state_outputs` on the CPU.

use anyhow::{Context, Result, ensure};
use cudarc::driver::{CudaContext, DevicePtr, DeviceRepr, LaunchConfig, PushKernelArg};
use cudarc::nvrtc::Ptx;
use foehncraft_kernels::entry::Raw;
use foehncraft_kernels::steam::{Input, OUT_WORDS, StateArgs};
use std::path::Path;

/// Threads per block (the states iterate and branch: small blocks).
const THREADS: u32 = 128;

/// The argument block as the kernel parameter (passed by value).
#[repr(transparent)]
struct Arg(StateArgs);
// SAFETY: StateArgs is a #[repr(C)] struct of an integer and device
// addresses with lengths; the PTX entry takes it by value with the same
// layout (one definition in foehncraft-kernels).
unsafe impl DeviceRepr for Arg {}

/// The states of the points (a[i], b[i]) on GPU `ordinal`: [`OUT_WORDS`]
/// outputs per point (ρ, h, T, c_p, μ, k, x, region), as
/// `foehncraft_kernels::steam::state_outputs`.
pub fn states(ordinal: usize, ptx: &Path, input: Input, a: &[f64], b: &[f64]) -> Result<Vec<f64>> {
    crate::driver()?;
    ensure!(
        a.len() == b.len(),
        "steam states: inputs of different length"
    );
    ensure!(
        ptx.exists(),
        "PTX module `{}` not found; build it with python3 gpu/ptx-kernels/build.py",
        ptx.display()
    );
    let n = a.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let ctx = CudaContext::new(ordinal).context("CUDA device (is the NVIDIA driver installed?)")?;
    let module = ctx
        .load_module(Ptx::from_src(std::fs::read_to_string(ptx)?))
        .context("loading the PTX module")?;
    let f = module
        .load_function("steam_state")
        .context("entry point steam_state (rebuild the PTX module)")?;
    let s = ctx.default_stream();
    let da = s.clone_htod(a)?;
    let db = s.clone_htod(b)?;
    let out = s.alloc_zeros::<f64>(n * OUT_WORDS)?;
    let args = StateArgs {
        input: input as u32,
        a: Raw::new(da.device_ptr(&s).0, n),
        b: Raw::new(db.device_ptr(&s).0, n),
        out: Raw::new(out.device_ptr(&s).0, n * OUT_WORDS),
    };
    let cfg = LaunchConfig {
        grid_dim: ((n as u32).div_ceil(THREADS), 1, 1),
        block_dim: (THREADS, 1, 1),
        shared_mem_bytes: 0,
    };
    let arg = Arg(args);
    let mut l = s.launch_builder(&f);
    l.arg(&arg);
    // SAFETY: steam_state takes (StateArgs); its addresses are the device
    // buffers above, alive until after the synchronisation, with n points
    // and OUT_WORDS outputs per point.
    unsafe { l.launch(cfg) }?;
    s.synchronize()?;
    Ok(s.clone_dtoh(&out)?)
}
