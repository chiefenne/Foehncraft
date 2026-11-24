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
use foehncraft_core::exact::{ExactSum, max_total, min_total};
use foehncraft_fv::amg::{AmgSettings, ReductionSettings};
use foehncraft_fv::backend::{Backend, V};
use foehncraft_fv::solver::{self as linear, Preconditioner};
use foehncraft_kernels::entry::{Entry, Raw};
use foehncraft_mesh::decompose::Ranks;
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

/// What the entries of a buffer are, on ranks (MPI step W5): which of them
/// a reduction takes and whether the halo is exchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    /// One value per own and halo cell (allocated for `n_cells` or
    /// `n_cells + n_halo` values): reductions take the own cells, the halo
    /// is filled from the owners.
    Cell,
    /// One value per face: reductions take the faces whose owner cell is
    /// own (a processor face counts once, on the rank of its owner).
    Face,
    /// Anything else: every entry is the rank's own.
    Plain,
}

/// The rank's link to the other ranks and the layout of its local mesh.
struct OnRanks {
    ranks: Ranks,
    n_cells: usize,
    n_halo: usize,
    n_faces: usize,
    /// The faces whose owner cell is own.
    face_own: Vec<bool>,
}

/// The CPU executor.
pub struct CpuExec {
    bufs: Vec<HostArray<f64>>,
    /// Kind and requested length of every buffer.
    kinds: Vec<(Kind, usize)>,
    statics: Vec<HostArray<u32>>,
    pool: rayon::ThreadPool,
    /// The pool has one thread.
    serial: bool,
    be: Box<dyn Backend>,
    ranks: Option<OnRanks>,
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
        let (be, ranks): (Box<dyn Backend>, _) = match &s.ranks {
            None => (
                foehncraft_fv::backend::create("cpu", s.disc, threads, std::path::Path::new(""))?,
                None,
            ),
            Some(r) => {
                let d = s.disc;
                ensure!(
                    d.n_faces != d.n_cells,
                    "on ranks a part needs more faces than cells"
                );
                let face_own = (0..d.n_faces)
                    .map(|f| (d.owner[f] as usize) < d.n_cells)
                    .collect();
                (
                    Box::new(foehncraft_fv::backend::cpu::CpuBackend::on_ranks(
                        d.clone(),
                        threads,
                        r.clone(),
                    )?),
                    Some(OnRanks {
                        ranks: r.clone(),
                        n_cells: d.n_cells,
                        n_halo: d.n_halo,
                        n_faces: d.n_faces,
                        face_own,
                    }),
                )
            }
        };
        Ok(Self {
            bufs: Vec::new(),
            kinds: Vec::new(),
            statics: Vec::new(),
            serial: pool.current_num_threads() == 1,
            pool,
            be,
            ranks,
        })
    }

    /// The buffer whose values start at `ptr`, if any.
    fn buffer_at(&self, ptr: u64) -> Option<usize> {
        self.bufs.iter().position(|b| b.ptr as u64 == ptr)
    }

    /// The entries of `r` a reduction takes on ranks: the own cells of a
    /// cell buffer, the own faces of a face buffer (as a mask), all others.
    fn own_part(&self, r: Raw<f64>) -> (usize, Option<&[bool]>) {
        let Some(o) = &self.ranks else {
            return (r.len as usize, None);
        };
        match self.buffer_at(r.ptr).map(|i| self.kinds[i].0) {
            Some(Kind::Cell) => (o.n_cells, None),
            Some(Kind::Face) => (o.n_faces, Some(&o.face_own)),
            _ => (r.len as usize, None),
        }
    }

    fn b(&self, h: Buf) -> &[f64] {
        self.bufs[h.0 as usize].slice()
    }
}

/// Values per task of the exact reductions.
const EXACT_GRAIN: usize = 1 << 12;

impl CpuExec {
    /// The values of a raw array of this executor.
    fn host(&self, r: Raw<f64>) -> &[f64] {
        if r.len == 0 {
            return &[];
        }
        // SAFETY: every Raw handed to the executor points into one of its
        // host arrays (buffers or statics), which live as long as `self`;
        // no launch runs while `self` is borrowed.
        unsafe { std::slice::from_raw_parts(r.ptr as *const f64, r.len as usize) }
    }

    fn exact_acc(&self, x: &[f64], y: Option<&[f64]>) -> ExactSum {
        let chunk = |(a, b): (&[f64], Option<&[f64]>)| {
            let mut s = ExactSum::new();
            match b {
                Some(b) => s.add_products(a, b),
                None => s.add_products(a, a),
            }
            s
        };
        if self.serial || x.len() <= EXACT_GRAIN {
            return chunk((x, y));
        }
        let merge = |mut a: ExactSum, b: ExactSum| {
            a.merge(&b);
            a
        };
        self.pool.install(|| match y {
            Some(y) => x
                .par_chunks(EXACT_GRAIN)
                .zip(y.par_chunks(EXACT_GRAIN))
                .map(|(a, b)| chunk((a, Some(b))))
                .reduce(ExactSum::new, merge),
            None => x
                .par_chunks(EXACT_GRAIN)
                .map(|a| chunk((a, None)))
                .reduce(ExactSum::new, merge),
        })
    }
}

impl Exec for CpuExec {
    fn name(&self) -> String {
        "cpu".into()
    }

    fn exact_dot(&mut self, x: Raw<f64>, y: Option<Raw<f64>>) -> Option<f64> {
        let xs = self.host(x);
        let ys = y.map(|y| self.host(y));
        if let Some(ys) = ys {
            assert_eq!(xs.len(), ys.len(), "dot product of unequal lengths");
        }
        let Some(o) = &self.ranks else {
            return Some(self.exact_acc(xs, ys).value());
        };
        // on ranks: the own entries, summed exactly over all ranks
        let (n, mask) = self.own_part(x);
        let acc = match mask {
            None => self.exact_acc(&xs[..n], ys.map(|y| &y[..n])),
            Some(m) => {
                let mut s = ExactSum::new();
                for (i, &own) in m.iter().enumerate() {
                    if own {
                        let a = xs[i];
                        s.add_products(&[a], &[ys.map_or(a, |y| y[i])]);
                    }
                }
                s
            }
        };
        Some(o.ranks.comm.sum_exact(&acc))
    }

    fn exact_min_max(&mut self, x: Raw<f64>) -> Option<(f64, f64)> {
        if let Some(o) = &self.ranks {
            // on ranks: the own entries, the extrema over all ranks
            let xs = self.host(x);
            let (n, mask) = self.own_part(x);
            let (lo, hi) = xs[..n]
                .iter()
                .enumerate()
                .filter(|(i, _)| mask.is_none_or(|m| m[*i]))
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), (_, &v)| {
                    (min_total(a, v), max_total(b, v))
                });
            return Some((o.ranks.comm.min(lo), o.ranks.comm.max(hi)));
        }
        let xs = self.host(x);
        let fold = |c: &[f64]| {
            c.iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
                    (min_total(a, v), max_total(b, v))
                })
        };
        if self.serial || xs.len() <= EXACT_GRAIN {
            return Some(fold(xs));
        }
        Some(self.pool.install(|| {
            xs.par_chunks(EXACT_GRAIN).map(fold).reduce(
                || (f64::INFINITY, f64::NEG_INFINITY),
                |(a, b), (c, d)| (min_total(a, c), max_total(b, d)),
            )
        }))
    }

    fn device(&self) -> String {
        format!("CPU, {} threads", self.pool.current_num_threads())
    }

    fn alloc(&mut self, n: usize) -> Result<Buf> {
        // on ranks a cell buffer has room for the halo; the kind is agreed
        // by all ranks (a boundary buffer may have as many entries as the
        // part has cells on one rank), so that every rank exchanges the
        // same buffers
        let (kind, len) = match &self.ranks {
            Some(o) => {
                let all = |local: bool| o.ranks.comm.min(if local { 1.0 } else { 0.0 }) > 0.0;
                if all(n == o.n_cells || n == o.n_cells + o.n_halo) {
                    (Kind::Cell, o.n_cells + o.n_halo)
                } else if all(n == o.n_faces) {
                    (Kind::Face, n)
                } else {
                    (Kind::Plain, n)
                }
            }
            None => (Kind::Plain, n),
        };
        self.bufs.push(HostArray::new(vec![0.0; len]));
        self.kinds.push((kind, n));
        Ok(Buf(self.bufs.len() as u32 - 1))
    }

    fn alloc_plain(&mut self, n: usize) -> Result<Buf> {
        self.bufs.push(HostArray::new(vec![0.0; n]));
        self.kinds.push((Kind::Plain, n));
        Ok(Buf(self.bufs.len() as u32 - 1))
    }

    fn on_ranks(&self) -> bool {
        self.ranks.is_some()
    }

    fn exchange(&mut self, bufs: &[Buf]) -> Result<()> {
        let Some(o) = &self.ranks else {
            return Ok(());
        };
        let mut idx: Vec<usize> = bufs
            .iter()
            .map(|b| b.index())
            .filter(|&i| self.kinds[i].0 == Kind::Cell)
            .collect();
        idx.sort_unstable();
        idx.dedup();
        if idx.is_empty() {
            return Ok(());
        }
        let mut fields: Vec<&mut [f64]> = idx
            .iter()
            .map(|&i| {
                let a = &self.bufs[i];
                // SAFETY: distinct buffers (deduplicated), each holding len
                // initialised values; no launch runs during the exchange.
                unsafe { std::slice::from_raw_parts_mut(a.ptr, a.len) }
            })
            .collect();
        o.ranks.exchange(&mut fields);
        Ok(())
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
        let (kind, n) = self.kinds[b.index()];
        let dst = self.bufs[b.0 as usize].slice_mut();
        ensure!(
            dst.len() == v.len() || n == v.len(),
            "upload of {} values into {}",
            v.len(),
            n
        );
        dst[..v.len()].copy_from_slice(v);
        // own cells only: the halo from the owners
        if kind == Kind::Cell && v.len() < dst.len() {
            self.exchange(&[b])?;
        }
        Ok(())
    }

    fn download(&mut self, b: Buf) -> Result<Vec<f64>> {
        // the requested length: a cell buffer without its halo
        let n = self.kinds[b.index()].1;
        let n = match &self.ranks {
            Some(o) if self.kinds[b.index()].0 == Kind::Cell => o.n_cells,
            _ => n,
        };
        Ok(self.b(b)[..n].to_vec())
    }

    fn download_with_halo(&mut self, b: Buf) -> Result<Vec<f64>> {
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
        self.bufs[x.0 as usize].slice_mut()[..phi.len()].copy_from_slice(&phi);
        self.exchange(&[x])?;
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
        self.bufs[x.0 as usize].slice_mut()[..phi.len()].copy_from_slice(&phi);
        self.exchange(&[x])?;
        Ok(())
    }

    fn sync(&mut self) -> Result<()> {
        Ok(())
    }
}
