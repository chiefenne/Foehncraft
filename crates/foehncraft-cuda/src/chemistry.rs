//! Gas-phase chemistry in the cells of a flow (SF2b) on an NVIDIA GPU: the
//! entry `chemistry_react` of `gpu/ptx-kernels`, one GPU thread per cell.
//! The mechanism is uploaded once; each step uploads the states of the
//! cells in batches (the work space of the stiff integrator is about
//! 2 (N + 1)² values per cell, so a batch is limited to [`WORK_BYTES`]),
//! runs the kernel body that the CPU threads run, and reads the states and
//! statistics back. Each cell is one thread with the operations of the CPU
//! in the same order (no fused multiply-add, R10), so the GPU gives the
//! bits of the CPU.

use anyhow::{Context, Result, ensure};
use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, DevicePtr, DeviceRepr, LaunchConfig,
    PushKernelArg, ValidAsZeroBits,
};
use cudarc::nvrtc::Ptx;
use foehncraft_flow::chemistry::device::Device;
use foehncraft_kernels::chemistry::{
    Efficiency, Program, ReactArgs, Reaction, Settings, Stats, Term, work_len,
};
use foehncraft_kernels::entry::Raw;
use foehncraft_kernels::species::Species;
use std::path::Path;
use std::sync::Arc;

/// Threads per block (a thread holds a whole stiff integration: small
/// blocks).
const THREADS: u32 = 64;

/// The largest work space of one launch (bytes).
pub const WORK_BYTES: usize = 1 << 30;

/// The argument block as the kernel parameter (passed by value).
#[repr(transparent)]
struct Arg(ReactArgs);
// SAFETY: ReactArgs is a #[repr(C)] struct of integers and floating-point
// numbers; the PTX entry takes it by value with the same layout (one
// definition in foehncraft-kernels).
unsafe impl DeviceRepr for Arg {}

/// A plain-data record copied to the device as bytes.
#[repr(transparent)]
#[derive(Clone, Copy)]
struct Pod<T>(T);
// SAFETY: only instantiated for the #[repr(C)] records of the mechanism
// and Stats, which hold integers and floating-point numbers only.
unsafe impl<T: Copy> DeviceRepr for Pod<T> {}
// SAFETY: all-zero bits are a valid value of these records.
unsafe impl<T: Copy> ValidAsZeroBits for Pod<T> {}

/// A device copy of a host array and its address for the argument block.
struct Up<T> {
    buf: CudaSlice<Pod<T>>,
    raw: Raw<T>,
}

fn up<T: Copy>(s: &Arc<CudaStream>, v: &[T]) -> Result<Up<T>> {
    let w: Vec<Pod<T>> = v.iter().map(|&x| Pod(x)).collect();
    let buf = if w.is_empty() {
        s.alloc_zeros::<Pod<T>>(1)?
    } else {
        s.clone_htod(&w)?
    };
    let raw = Raw::new(buf.device_ptr(s).0, v.len());
    Ok(Up { buf, raw })
}

/// The chemistry of the cells on one GPU.
pub struct CudaChemistry {
    ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    f: CudaFunction,
    n: usize,
    _species: Up<Species>,
    _reactions: Up<Reaction>,
    _terms: Up<Term>,
    _efficiencies: Up<Efficiency>,
    /// The argument block with the mechanism filled in.
    mech: ReactArgs,
}

impl CudaChemistry {
    /// Uploads the mechanism to GPU `ordinal` and loads the PTX module.
    pub fn new(prog: &Program<'_>, ordinal: usize, ptx: &Path) -> Result<Self> {
        ensure!(
            ptx.exists(),
            "PTX module `{}` not found; build it with python3 gpu/ptx-kernels/build.py",
            ptx.display()
        );
        let ctx =
            CudaContext::new(ordinal).context("CUDA device (is the NVIDIA driver installed?)")?;
        let module = ctx
            .load_module(Ptx::from_src(std::fs::read_to_string(ptx)?))
            .context("loading the PTX module")?;
        let f = module
            .load_function("chemistry_react")
            .context("entry point chemistry_react (rebuild the PTX module)")?;
        let s = ctx.default_stream();
        let species = up(&s, prog.species)?;
        let reactions = up(&s, prog.reactions)?;
        let terms = up(&s, prog.terms)?;
        let efficiencies = up(&s, prog.efficiencies)?;
        let mech = ReactArgs {
            species: species.raw,
            reactions: reactions.raw,
            terms: terms.raw,
            efficiencies: efficiencies.raw,
            settings: Settings::default(),
            pressure: 0.0,
            dt: 0.0,
            min_temperature: 0.0,
            n_cells: 0,
            state: Raw::EMPTY,
            work: Raw::EMPTY,
            stats: Raw::EMPTY,
        };
        Ok(Self {
            ctx,
            stream: s,
            f,
            n: prog.n(),
            _species: species,
            _reactions: reactions,
            _terms: terms,
            _efficiencies: efficiencies,
            mech,
        })
    }

    /// Name and compute capability of the device.
    pub fn device(&self) -> String {
        let name = self.ctx.name().unwrap_or_else(|_| "NVIDIA GPU".into());
        let (ma, mi) = self.ctx.compute_capability().unwrap_or((0, 0));
        format!("{name} (sm_{ma}{mi})")
    }
}

impl Device for CudaChemistry {
    fn name(&self) -> &str {
        "cuda"
    }

    fn react(
        &mut self,
        p: f64,
        dt: f64,
        min_temperature: f64,
        set: &Settings,
        state: &mut [f64],
        stats: &mut [Stats],
    ) -> Result<()> {
        let nz = self.n + 1;
        let wl = work_len(self.n);
        ensure!(
            state.len() == stats.len() * nz,
            "chemistry launch: {} states for {} cells of {nz} values",
            state.len(),
            stats.len()
        );
        let batch = (WORK_BYTES / (8 * wl)).max(1);
        let s = self.stream.clone();
        for (zs, ss) in state.chunks_mut(batch * nz).zip(stats.chunks_mut(batch)) {
            let m = ss.len();
            let z = up(&s, zs)?;
            let work = s.alloc_zeros::<f64>(m * wl)?;
            let st = up(&s, ss)?;
            let a = ReactArgs {
                settings: *set,
                pressure: p,
                dt,
                min_temperature,
                n_cells: m as u64,
                state: z.raw,
                work: Raw::new(work.device_ptr(&s).0, m * wl),
                stats: st.raw,
                ..self.mech
            };
            let cfg = LaunchConfig {
                grid_dim: ((m as u32).div_ceil(THREADS), 1, 1),
                block_dim: (THREADS, 1, 1),
                shared_mem_bytes: 0,
            };
            let arg = Arg(a);
            let mut l = s.launch_builder(&self.f);
            l.arg(&arg);
            // SAFETY: chemistry_react takes (ReactArgs); every address in
            // it is a device buffer of this launch or of the uploaded
            // mechanism, alive until after the synchronisation, sized for
            // m cells.
            unsafe { l.launch(cfg) }?;
            s.synchronize()?;
            let back = s.clone_dtoh(&z.buf)?;
            for (d, b) in zs.iter_mut().zip(back) {
                *d = b.0;
            }
            let back = s.clone_dtoh(&st.buf)?;
            for (d, b) in ss.iter_mut().zip(back) {
                *d = b.0;
            }
            drop(work);
        }
        Ok(())
    }
}
