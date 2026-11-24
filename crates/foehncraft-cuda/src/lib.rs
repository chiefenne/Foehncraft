//! NVIDIA GPU backends of Foehncraft (feature `cuda`). This is the only crate
//! with CUDA types, calls or sources (D60, checked in CI by
//! `tools/check/backend.py`): the solver crates use the vendor-neutral
//! interfaces (`foehncraft_fv::backend::Backend`,
//! `foehncraft_flow::device::kernel::Exec`, `foehncraft_bench::backend::Backend`)
//! and find the GPU backends through the factories that [`register`]
//! installs. A backend for another vendor is another crate of this kind.
//!
//! - `flow`: the executor of the flow loop (entries `flow_*`);
//! - `fv`: the scalar operators, Krylov vector operations and AMG V-cycle
//!   (entries `fv_*`, `blas_*`, `amg_*`);
//! - `radiation`: the transport of the radiation bundles (entry
//!   `radiation_transport`, step 29);
//! - `bench`: the G1 benchmark backends `cuda-reference` (hand-written
//!   CUDA C, `kernels/g1_reference.cu`) and `cuda-rust` (`g1_*`).
//!
//! All PTX entries come from `gpu/ptx-kernels`, built with
//! `python3 gpu/ptx-kernels/build.py` (nightly Rust, sm_70 and newer).
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

#[cfg(feature = "cuda")]
pub mod bench;
#[cfg(feature = "cuda")]
pub mod flow;
#[cfg(feature = "cuda")]
pub mod fv;
#[cfg(feature = "cuda")]
pub mod radiation;

/// Fails with a clear message when the NVIDIA driver library is missing
/// (the bindings would panic on first use).
#[cfg(feature = "cuda")]
pub(crate) fn driver() -> anyhow::Result<()> {
    // SAFETY: only tries to open the shared library.
    let present = unsafe { cudarc::driver::sys::is_culib_present() };
    anyhow::ensure!(
        present,
        "the cuda backend needs an NVIDIA GPU and its driver (libcuda was not found)"
    );
    Ok(())
}

/// Registers the GPU backends of this build with the solver crates (`cuda`
/// for the flow loop, the scalar operators and the radiation, `cuda-reference` and
/// `cuda-rust` for G1). Without the `cuda` feature it does nothing. Call it
/// once at start-up, before any case is run.
pub fn register() {
    #[cfg(feature = "cuda")]
    {
        foehncraft_fv::backend::register_gpu("cuda", |d, ordinal, ptx| {
            driver()?;
            Ok(Box::new(fv::CudaBackend::new(d, ordinal, ptx)?))
        });
        foehncraft_flow::device::register_gpu("cuda", flow::create);
        foehncraft_radiation::register_gpu("cuda", |m, ordinal, ptx| {
            driver()?;
            Ok(Box::new(radiation::CudaRadiation::new(m, ordinal, ptx)?))
        });
        foehncraft_bench::backend::register_gpu(&["cuda-reference", "cuda-rust"], |name, o| {
            driver()?;
            Ok(Box::new(match name {
                "cuda-reference" => bench::CudaBackend::reference(o.device)?,
                _ => bench::CudaBackend::rust(o.device, &o.rust_ptx)?,
            }))
        });
    }
}

/// True when this build has the GPU backends.
pub const AVAILABLE: bool = cfg!(feature = "cuda");
