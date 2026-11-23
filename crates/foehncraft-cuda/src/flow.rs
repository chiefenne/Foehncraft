//! NVIDIA GPU executor of the flow loop: the entries of
//! `foehncraft_kernels::entry` as the PTX entry points `flow_*`
//! (`gpu/ptx-kernels/src/flow.rs`), one thread per index, on buffers that
//! stay in GPU memory for the whole run. The argument block of an entry is
//! passed by value as a kernel parameter; it is the same `#[repr(C)]` struct
//! the CPU executor passes to the same body, so the two cannot disagree on
//! the arguments; every entry gives the CPU's bits (no fused multiply-adds
//! in the PTX, IEEE division). The reductions differ: fixed blocks and tree
//! here, exact sums on the CPU since MPI step W2, so whole runs agree to
//! rounding, not bitwise, until v1.1 (owner decision of 2026-11-23).
//!
//! The linear solves run on the GPU backend of `foehncraft-fv`
//! (`crate::fv::CudaBackend`) with device-to-device copies of the matrix,
//! right-hand side and solution; the multigrid setup is built on the host
//! from a download of the matrix (only when the hierarchy is rebuilt).

use crate::fv::CudaBackend;
use anyhow::{Context as _, Result, ensure};
use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, DevicePtr, DeviceRepr, LaunchConfig,
    PushKernelArg,
};
use cudarc::nvrtc::Ptx;
use foehncraft_flow::device::kernel::{Exec, KernelFlow};
use foehncraft_flow::device::{Buf, FlowDevice, Setup};
use foehncraft_fv::amg::{AmgSettings, ReductionSettings};
use foehncraft_fv::backend::{Backend, V};
use foehncraft_fv::solver::{self as linear, Preconditioner};
use foehncraft_kernels::entry::{self as e, Entry, Mesh, Raw};
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

/// Threads per block.
const THREADS: u32 = 256;

/// An argument block as a kernel parameter (passed by value, like a struct
/// argument of a CUDA C kernel).
#[repr(transparent)]
struct Arg<A>(A);
// SAFETY: the argument blocks are #[repr(C)] structs of integers and
// floating-point numbers; the PTX entries take them by value with the same
// layout (one definition in foehncraft-kernels).
unsafe impl<A: Copy> DeviceRepr for Arg<A> {}

fn cfg(n: usize) -> LaunchConfig {
    LaunchConfig {
        grid_dim: ((n as u32).div_ceil(THREADS).max(1), 1, 1),
        block_dim: (THREADS, 1, 1),
        shared_mem_bytes: 0,
    }
}

/// The GPU executor.
pub struct CudaExec {
    ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    f: HashMap<&'static str, CudaFunction>,
    bufs: Vec<CudaSlice<f64>>,
    /// Address and length of each buffer (allocations hold at least one
    /// entry; the length is the requested one).
    raws: Vec<Raw<f64>>,
    statics: Vec<CudaSlice<u32>>,
    /// The matrix layout for the Jacobi diagonal of a device matrix.
    layout: Mesh,
    be: CudaBackend,
}

impl CudaExec {
    /// An executor on GPU `ordinal` with the PTX module at `ptx`.
    pub fn new(s: &Setup, ordinal: usize, ptx: &Path) -> Result<Self> {
        crate::driver()?;
        ensure!(
            ptx.exists(),
            "PTX module `{}` not found; build it with python3 gpu/ptx-kernels/build.py",
            ptx.display()
        );
        let be = CudaBackend::new(s.disc, ordinal, ptx)?;
        let ctx =
            CudaContext::new(ordinal).context("CUDA device (is the NVIDIA driver installed?)")?;
        let module = ctx
            .load_module(Ptx::from_src(std::fs::read_to_string(ptx)?))
            .context("loading the PTX module")?;
        let mut f = HashMap::new();
        for name in e::NAMES {
            let entry = format!("flow_{name}");
            f.insert(
                name,
                module.load_function(&entry).with_context(|| {
                    format!(
                        "entry point {entry} (rebuild the PTX module with python3 gpu/ptx-kernels/build.py)"
                    )
                })?,
            );
        }
        let stream = be.stream().clone();
        let mut x = Self {
            ctx,
            stream,
            f,
            bufs: Vec::new(),
            raws: Vec::new(),
            statics: Vec::new(),
            layout: Mesh::default(),
            be,
        };
        x.layout.n_cells = s.disc.n_cells as u64;
        x.layout.chunk = s.disc.sell.chunk as u64;
        x.layout.slice_ptr = x.alloc_u32(&s.disc.sell.slice_ptr)?;
        Ok(x)
    }

    fn ptr<T>(&self, s: &CudaSlice<T>) -> u64 {
        s.device_ptr(&self.stream).0
    }

    /// Two buffers at once, the second mutable (distinct).
    fn pair(&mut self, a: Buf, b: Buf) -> (&CudaSlice<f64>, &mut CudaSlice<f64>) {
        let (i, j) = (a.index(), b.index());
        assert_ne!(i, j);
        if i < j {
            let (l, r) = self.bufs.split_at_mut(j);
            (&l[i], &mut r[0])
        } else {
            let (l, r) = self.bufs.split_at_mut(i);
            (&r[0], &mut l[j])
        }
    }

    /// Copies a flow buffer into a vector of the linear solvers.
    fn copy_to_solver(&mut self, b: Buf, v: V) -> Result<()> {
        let s = self.stream.clone();
        let src = &self.bufs[b.index()];
        let n = self.raws[b.index()].len as usize;
        let dst = self.be.vector_mut(v);
        ensure!(
            dst.len() == n,
            "solver vector {v:?}: {} entries, buffer {n}",
            dst.len()
        );
        s.memcpy_dtod(&src.slice(..n), dst)?;
        Ok(())
    }

    /// Puts the matrix `val` into the linear solvers, with its Jacobi
    /// diagonal computed on the device (as the CPU backend does on the
    /// host, `1 / val[slot 0]`).
    fn set_matrix(&mut self, val: Buf) -> Result<()> {
        let s = self.stream.clone();
        let n = self.raws[val.index()].len as usize;
        {
            let src = &self.bufs[val.index()];
            let dst = self.be.matrix_mut();
            ensure!(dst.len() == n, "matrix of {n} values, solver {}", dst.len());
            s.memcpy_dtod(&src.slice(..n), dst)?;
        }
        let dinv = self.be.vector(V::Dinv);
        let a = e::UnaryArgs {
            m: self.layout,
            s: 0.0,
            x: self.raws[val.index()],
            out: Raw::new(self.ptr(dinv), dinv.len()),
        };
        let nc = self.layout.n_cells as usize;
        self.launch::<e::DiagonalInverse>(nc, &a)
    }
}

impl Exec for CudaExec {
    fn name(&self) -> String {
        "cuda".into()
    }

    fn device(&self) -> String {
        let name = self.ctx.name().unwrap_or_else(|_| "NVIDIA GPU".into());
        let (ma, mi) = self.ctx.compute_capability().unwrap_or((0, 0));
        format!("{name} (sm_{ma}{mi})")
    }

    fn alloc(&mut self, n: usize) -> Result<Buf> {
        let b = self.stream.alloc_zeros::<f64>(n.max(1))?;
        let r = Raw::new(self.ptr(&b), n);
        self.bufs.push(b);
        self.raws.push(r);
        Ok(Buf::new(self.bufs.len() - 1))
    }

    fn alloc_u32(&mut self, v: &[u32]) -> Result<Raw<u32>> {
        if v.is_empty() {
            return Ok(Raw::EMPTY);
        }
        let b = self.stream.clone_htod(v)?;
        let r = Raw::new(self.ptr(&b), v.len());
        self.statics.push(b);
        Ok(r)
    }

    fn raw(&self, b: Buf) -> Raw<f64> {
        self.raws[b.index()]
    }

    fn upload(&mut self, b: Buf, v: &[f64]) -> Result<()> {
        ensure!(
            self.raws[b.index()].len as usize == v.len(),
            "upload of {} values into {}",
            v.len(),
            self.raws[b.index()].len
        );
        if !v.is_empty() {
            let s = self.stream.clone();
            s.memcpy_htod(v, &mut self.bufs[b.index()])?;
        }
        Ok(())
    }

    fn download(&mut self, b: Buf) -> Result<Vec<f64>> {
        let n = self.raws[b.index()].len as usize;
        if n == 0 {
            return Ok(Vec::new());
        }
        Ok(self.stream.clone_dtoh(&self.bufs[b.index()].slice(..n))?)
    }

    fn copy(&mut self, src: Buf, dst: Buf) -> Result<()> {
        let n = self.raws[src.index()].len;
        ensure!(
            n == self.raws[dst.index()].len,
            "copy between buffers of different lengths"
        );
        let s = self.stream.clone();
        let (a, b) = self.pair(src, dst);
        s.memcpy_dtod(a, b)?;
        Ok(())
    }

    fn value(&mut self, b: Buf, i: usize) -> Result<f64> {
        ensure!(
            i < self.raws[b.index()].len as usize,
            "value: index {i} out of range"
        );
        Ok(self
            .stream
            .clone_dtoh(&self.bufs[b.index()].slice(i..i + 1))?[0])
    }

    fn launch<K: Entry>(&mut self, n: usize, a: &K::Args) -> Result<()> {
        let f = self
            .f
            .get(K::NAME)
            .with_context(|| format!("no entry point flow_{}", K::NAME))?;
        let n64 = n as u64;
        let arg = Arg(*a);
        let mut l = self.stream.launch_builder(f);
        l.arg(&n64).arg(&arg);
        // SAFETY: flow_<NAME> takes (u64, K::Args); the argument block holds
        // device buffers of this executor that live as long as it does.
        unsafe { l.launch(cfg(n)) }?;
        Ok(())
    }

    fn solve(
        &mut self,
        val: Buf,
        rhs: Buf,
        x: Buf,
        pc: Preconditioner,
        cg: bool,
        tol: f64,
        max_iter: usize,
    ) -> Result<(usize, f64)> {
        self.set_matrix(val)?;
        self.copy_to_solver(rhs, V::Rhs)?;
        self.copy_to_solver(x, V::Phi)?;
        let be: &mut dyn Backend = &mut self.be;
        be.residual(V::Phi, V::Rhs, V::R)?;
        let r0 = be.dot(V::R, V::R)?.sqrt();
        let it = if r0 == 0.0 {
            0
        } else if cg {
            linear::cg(be, pc, max_iter, tol, r0)?.0
        } else {
            linear::bicgstab(be, pc, max_iter, tol, r0)?
        };
        let s = self.stream.clone();
        let n = self.raws[x.index()].len as usize;
        let src = self.be.vector(V::Phi);
        s.memcpy_dtod(src, &mut self.bufs[x.index()].slice_mut(..n))?;
        Ok((it, r0))
    }

    fn amg_setup(&mut self, val: Buf, s: &AmgSettings, symmetric_part: bool) -> Result<()> {
        let host = self.download(val)?;
        self.be.set_matrix(&host)?;
        linear::setup_amg_with(&mut self.be, s, symmetric_part).map(|_| ())
    }

    fn amg_slot(&mut self, slot: usize) -> Result<()> {
        self.be.amg_slot(slot)
    }

    fn reduction_setup(&mut self, val: Buf, s: &ReductionSettings) -> Result<()> {
        let host = self.download(val)?;
        self.be.set_matrix(&host)?;
        linear::setup_reduction(&mut self.be, s).map(|_| ())
    }

    fn cycles(&mut self, val: Buf, rhs: Buf, x: Buf, cycles: usize) -> Result<()> {
        if cycles > 1 {
            // the residual between two cycles takes the matrix of the backend
            self.set_matrix(val)?;
        }
        self.copy_to_solver(rhs, V::Rhs)?;
        let be: &mut dyn Backend = &mut self.be;
        linear::cycles(be, cycles)?;
        let s = self.stream.clone();
        let n = self.raws[x.index()].len as usize;
        let src = self.be.vector(V::Phi);
        s.memcpy_dtod(src, &mut self.bufs[x.index()].slice_mut(..n))?;
        Ok(())
    }

    fn sync(&mut self) -> Result<()> {
        Ok(self.stream.synchronize()?)
    }
}

/// A flow device on GPU 0 (registered as `cuda` by [`crate::register`]).
pub fn create(s: &Setup, ptx: &Path) -> Result<Box<dyn FlowDevice>> {
    Ok(Box::new(KernelFlow::new(CudaExec::new(s, 0, ptx)?, s)?))
}
