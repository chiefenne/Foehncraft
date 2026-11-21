//! CPU executor: the entries of `foehncraft_kernels::entry` in
//! multi-threaded loops (rayon) over host memory, the linear solves on the
//! CPU backend of `foehncraft-fv`. Each output is written by one index, so
//! the result is the same for any number of threads (D31), and the same as
//! on the GPU, which runs the same bodies.
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use super::kernel::{Exec, KernelFlow};
use super::{Buf, Setup};
use anyhow::{Result, ensure};
use foehncraft_fv::amg::{AmgSettings, ReductionSettings};
use foehncraft_fv::backend::{Backend, V};
use foehncraft_fv::solver::{self as linear, Preconditioner};
use foehncraft_kernels::entry::{Entry, Raw};
use rayon::prelude::*;

/// Indices per task of the loops.
const GRAIN: usize = 1 << 12;

/// Meshes with fewer cells run on one thread.
pub const SERIAL_BELOW: usize = 20_000;

/// The CPU flow device.
pub type CpuFlow = KernelFlow<CpuExec>;

impl CpuFlow {
    /// A CPU device for the static data of `s`.
    pub fn cpu(s: &Setup) -> Result<Self> {
        KernelFlow::new(CpuExec::new(s)?, s)
    }
}

/// A heap array addressed only through its raw pointer (the entries write
/// through the addresses of their argument blocks).
struct HostArray<T> {
    ptr: *mut T,
    len: usize,
}

impl<T: Copy + Default> HostArray<T> {
    fn new(v: Vec<T>) -> Self {
        let len = v.len();
        let ptr = Box::into_raw(v.into_boxed_slice()) as *mut T;
        Self { ptr, len }
    }

    fn raw(&self) -> Raw<T> {
        Raw::new(self.ptr as u64, self.len)
    }

    fn slice(&self) -> &[T] {
        // SAFETY: ptr holds len initialised values; no launch runs while the
        // executor is borrowed.
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }

    fn slice_mut(&mut self) -> &mut [T] {
        // SAFETY: as above, exclusive through &mut self.
        unsafe { std::slice::from_raw_parts_mut(self.ptr, self.len) }
    }
}

impl<T> Drop for HostArray<T> {
    fn drop(&mut self) {
        // SAFETY: ptr and len come from Box::into_raw of a boxed slice.
        unsafe {
            drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(
                self.ptr, self.len,
            )))
        }
    }
}

// SAFETY: the arrays are plain heap memory owned by the executor.
unsafe impl<T: Send> Send for HostArray<T> {}

/// The CPU executor.
pub struct CpuExec {
    bufs: Vec<HostArray<f64>>,
    statics: Vec<HostArray<u32>>,
    pool: rayon::ThreadPool,
    /// The pool has one thread.
    serial: bool,
    be: Box<dyn Backend>,
}

impl CpuExec {
    /// An executor with the linear solvers on the discretisation of `s`.
    pub fn new(s: &Setup) -> Result<Self> {
        // small meshes run on one thread: the latency of waking the pool
        // would dominate (the result is the same for any number of threads)
        let threads = if s.disc.n_cells < SERIAL_BELOW {
            Some(1)
        } else {
            s.threads
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.unwrap_or(0))
            .build()?;
        let be = foehncraft_fv::backend::create("cpu", s.disc, threads, std::path::Path::new(""))?;
        Ok(Self {
            bufs: Vec::new(),
            statics: Vec::new(),
            serial: pool.current_num_threads() == 1,
            pool,
            be,
        })
    }

    fn b(&self, h: Buf) -> &[f64] {
        self.bufs[h.0 as usize].slice()
    }
}

impl Exec for CpuExec {
    fn name(&self) -> String {
        "cpu".into()
    }

    fn device(&self) -> String {
        format!("CPU, {} threads", self.pool.current_num_threads())
    }

    fn alloc(&mut self, n: usize) -> Result<Buf> {
        self.bufs.push(HostArray::new(vec![0.0; n]));
        Ok(Buf(self.bufs.len() as u32 - 1))
    }

    fn alloc_u32(&mut self, v: &[u32]) -> Result<Raw<u32>> {
        let a = HostArray::new(v.to_vec());
        let r = a.raw();
        self.statics.push(a);
        Ok(r)
    }

    fn raw(&self, b: Buf) -> Raw<f64> {
        self.bufs[b.0 as usize].raw()
    }

    fn upload(&mut self, b: Buf, v: &[f64]) -> Result<()> {
        let dst = self.bufs[b.0 as usize].slice_mut();
        ensure!(
            dst.len() == v.len(),
            "upload of {} values into {}",
            v.len(),
            dst.len()
        );
        dst.copy_from_slice(v);
        Ok(())
    }

    fn download(&mut self, b: Buf) -> Result<Vec<f64>> {
        Ok(self.b(b).to_vec())
    }

    fn copy(&mut self, src: Buf, dst: Buf) -> Result<()> {
        ensure!(src != dst, "copy onto itself");
        let (s, d) = (&self.bufs[src.0 as usize], &self.bufs[dst.0 as usize]);
        ensure!(s.len == d.len, "copy of {} values into {}", s.len, d.len);
        // SAFETY: two distinct arrays of the same length.
        unsafe { std::ptr::copy_nonoverlapping(s.ptr, d.ptr, s.len) };
        Ok(())
    }

    fn value(&mut self, b: Buf, i: usize) -> Result<f64> {
        Ok(self.b(b)[i])
    }

    fn launch<K: Entry>(&mut self, n: usize, a: &K::Args) -> Result<()> {
        let task = |ch: usize| {
            for i in ch * GRAIN..((ch + 1) * GRAIN).min(n) {
                // SAFETY: i < n; the device built the argument block from
                // live buffers and checked that outputs are distinct.
                unsafe { K::run(i, a) }
            }
        };
        // Handing a loop to the pool costs two thread switches, more than
        // the loop itself on a small mesh: a pool of one thread and a loop
        // of one task run on the calling thread (same values, D31).
        if self.serial || n <= GRAIN {
            (0..n.div_ceil(GRAIN)).for_each(task);
        } else {
            self.pool
                .install(|| (0..n.div_ceil(GRAIN)).into_par_iter().for_each(task));
        }
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
        let be = self.be.as_mut();
        be.set_matrix(self.bufs[val.0 as usize].slice())?;
        be.upload(V::Rhs, self.bufs[rhs.0 as usize].slice())?;
        be.upload(V::Phi, self.bufs[x.0 as usize].slice())?;
        be.residual(V::Phi, V::Rhs, V::R)?;
        let r0 = be.dot(V::R, V::R)?.sqrt();
        let it = if r0 == 0.0 {
            0
        } else if cg {
            linear::cg(be, pc, max_iter, tol, r0)?.0
        } else {
            linear::bicgstab(be, pc, max_iter, tol, r0)?
        };
        let phi = be.download(V::Phi)?;
        self.bufs[x.0 as usize].slice_mut().copy_from_slice(&phi);
        Ok((it, r0))
    }

    fn amg_setup(&mut self, val: Buf, s: &AmgSettings, symmetric_part: bool) -> Result<()> {
        self.be.set_matrix(self.bufs[val.0 as usize].slice())?;
        linear::setup_amg_with(self.be.as_mut(), s, symmetric_part).map(|_| ())
    }

    fn amg_slot(&mut self, slot: usize) -> Result<()> {
        self.be.amg_slot(slot)
    }

    fn reduction_setup(&mut self, val: Buf, s: &ReductionSettings) -> Result<()> {
        self.be.set_matrix(self.bufs[val.0 as usize].slice())?;
        linear::setup_reduction(self.be.as_mut(), s).map(|_| ())
    }

    fn cycles(&mut self, val: Buf, rhs: Buf, x: Buf, cycles: usize) -> Result<()> {
        let be = self.be.as_mut();
        if cycles > 1 {
            // the residual between two cycles takes the matrix of the backend
            be.set_matrix(self.bufs[val.0 as usize].slice())?;
        }
        be.upload(V::Rhs, self.bufs[rhs.0 as usize].slice())?;
        linear::cycles(be, cycles)?;
        let phi = be.download(V::Phi)?;
        self.bufs[x.0 as usize].slice_mut().copy_from_slice(&phi);
        Ok(())
    }

    fn sync(&mut self) -> Result<()> {
        Ok(())
    }
}
