//! CPU backend: the kernel bodies in multi-threaded loops (rayon).
//!
//! Each output value is computed by one iteration from its index, so the
//! results are bitwise the same for any number of threads (D31).
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use super::{Backend, FaceSettings, Limiter, V};
use crate::amg::cycle::LevelOps;
use crate::amg::{Hierarchy, Precision, SellMatrix};
use crate::disc::Discretisation;
use foehncraft_core::Real;
use foehncraft_core::exact::{ExactSum, max_total, min_total};
use foehncraft_kernels::fv::{self, FaceInputs, SellLayout, Vec3};
use foehncraft_kernels::{amg, blas, spmv};
use foehncraft_mesh::decompose::Ranks;
use rayon::prelude::*;

/// Indices per task of the parallel loops.
const GRAIN: usize = 1 << 14;

/// Values per task of the exact reductions.
const EXACT_GRAIN: usize = 1 << 12;

/// Exact dot product of `x` and `y` (D31 (1)): the exact sum of the rounded
/// products, rounded once, so the bits depend neither on the thread count
/// nor on the order or the split of the cells between ranks.
pub fn exact_dot(w: &Workers, x: &[f64], y: &[f64]) -> f64 {
    exact_dot_acc(w, x, y).value()
}

/// The accumulator of [`exact_dot`], for a reduction across ranks.
pub fn exact_dot_acc(w: &Workers, x: &[f64], y: &[f64]) -> ExactSum {
    assert_eq!(x.len(), y.len(), "dot product of unequal lengths");
    let chunk = |(a, b): (&[f64], &[f64])| {
        let mut s = ExactSum::new();
        s.add_products(a, b);
        s
    };
    let merge = |mut a: ExactSum, b: ExactSum| {
        a.merge(&b);
        a
    };
    if w.inline(x.len()) {
        chunk((x, y))
    } else {
        w.install(|| {
            x.par_chunks(EXACT_GRAIN)
                .zip(y.par_chunks(EXACT_GRAIN))
                .map(chunk)
                .reduce(ExactSum::new, merge)
        })
    }
}

/// Minimum and maximum of `x`, independent of order and grouping
/// ([`min_total`], [`max_total`]); (+∞, −∞) for no values.
pub fn min_max(w: &Workers, x: &[f64]) -> (f64, f64) {
    let fold = |c: &[f64]| {
        c.iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
                (min_total(a, v), max_total(b, v))
            })
    };
    if w.inline(x.len()) {
        fold(x)
    } else {
        w.install(|| {
            x.par_chunks(EXACT_GRAIN).map(fold).reduce(
                || (f64::INFINITY, f64::NEG_INFINITY),
                |(a, b), (c, d)| (min_total(a, c), max_total(b, d)),
            )
        })
    }
}

/// The threads of the backend. Handing a loop to the pool costs two thread
/// switches, more than the loop itself on a small mesh or on a coarse level
/// of a multigrid cycle. So a pool of one thread and a loop of one task run
/// on the calling thread. The values do not depend on where a loop runs
/// (D31).
pub struct Workers {
    pool: rayon::ThreadPool,
    serial: bool,
}

impl Workers {
    /// A pool of `threads` threads (`None`: one per core).
    pub fn new(threads: Option<usize>) -> anyhow::Result<Self> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.unwrap_or(0))
            .build()?;
        let serial = pool.current_num_threads() == 1;
        Ok(Self { pool, serial })
    }

    fn threads(&self) -> usize {
        self.pool.current_num_threads()
    }

    /// A loop over `n` indices runs on the calling thread.
    fn inline(&self, n: usize) -> bool {
        self.serial || n <= GRAIN
    }

    fn install<R: Send>(&self, f: impl FnOnce() -> R + Send) -> R {
        self.pool.install(f)
    }
}

/// CPU backend.
pub struct CpuBackend {
    d: Discretisation,
    pool: Workers,
    vecs: Vec<Vec<f64>>,
    scratch: Vec<f64>,
    val: Vec<f64>,
    eps2: Option<(f64, Vec<f64>)>,
    amg: Option<Amg>,
    /// The hierarchies of the inactive slots, and the active slot.
    parked: Vec<Option<Amg>>,
    slot: usize,
    /// The other ranks of a decomposed run (MPI step W3); `None` on one
    /// partition.
    ranks: Option<Ranks>,
}

/// One level of the hierarchy in precision `T`, with its cycle vectors.
struct LevelData<T> {
    a: SellMatrix<T>,
    dinv: Vec<T>,
    p: SellMatrix<T>,
    r: SellMatrix<T>,
    steps: Vec<(f64, f64)>,
}

/// The hierarchy in precision `T`; `x`, `b` have L + 1 entries (the
/// coarsest last), `w`, `d` have L.
struct AmgData<T> {
    levels: Vec<LevelData<T>>,
    coarse_n: usize,
    coarse_inv: Vec<T>,
    x: Vec<Vec<T>>,
    b: Vec<Vec<T>>,
    w: Vec<Vec<T>>,
    d: Vec<Vec<T>>,
}

enum Amg {
    F32(AmgData<f32>),
    F64(AmgData<f64>),
}

/// Fills `out[i] = f(i)` in parallel, one value per index.
fn par_fill<T: Real>(pool: &Workers, out: &mut [T], f: impl Fn(usize, T) -> T + Sync) {
    if pool.inline(out.len()) {
        for (i, o) in out.iter_mut().enumerate() {
            *o = f(i, *o);
        }
        return;
    }
    pool.install(|| {
        out.par_chunks_mut(GRAIN).enumerate().for_each(|(c, ch)| {
            for (j, o) in ch.iter_mut().enumerate() {
                *o = f(c * GRAIN + j, *o);
            }
        })
    });
}

impl<T: Real> AmgData<T> {
    fn new(h: &Hierarchy) -> Self {
        let cv = |v: &[f64]| v.iter().map(|&x| T::from_f64(x)).collect::<Vec<T>>();
        let levels: Vec<LevelData<T>> = h
            .levels
            .iter()
            .map(|l| LevelData {
                a: SellMatrix::from_csr(&l.a),
                dinv: cv(&l.dinv),
                p: SellMatrix::from_csr(&l.p),
                r: SellMatrix::from_csr(&l.r),
                steps: l.steps.clone(),
            })
            .collect();
        let mut sizes: Vec<usize> = h.levels.iter().map(|l| l.a.n_rows).collect();
        sizes.push(h.coarse_n);
        let zeros = |n: usize| vec![T::ZERO; n];
        Self {
            coarse_n: h.coarse_n,
            coarse_inv: cv(&h.coarse_inv),
            x: sizes.iter().map(|&n| zeros(n)).collect(),
            b: sizes.iter().map(|&n| zeros(n)).collect(),
            w: sizes[..levels.len()].iter().map(|&n| zeros(n)).collect(),
            d: sizes[..levels.len()].iter().map(|&n| zeros(n)).collect(),
            levels,
        }
    }

    fn load_rhs(&mut self, pool: &Workers, r: &[f64]) {
        // SAFETY: b_0 and r have one entry per cell.
        par_fill(pool, &mut self.b[0], |i, _| unsafe { amg::convert(i, r) });
    }

    fn store(&self, pool: &Workers, z: &mut [f64]) {
        let x = &self.x[0];
        // SAFETY: as above.
        par_fill(pool, z, |i, _| unsafe { amg::convert(i, x) });
    }

    fn smooth_first(&mut self, pool: &Workers, l: usize, c2: f64) {
        let c2 = T::from_f64(c2);
        let (dinv, b) = (&self.levels[l].dinv, &self.b[l]);
        // SAFETY: level vectors have one entry per row of A_l.
        par_fill(pool, &mut self.d[l], |i, _| unsafe {
            amg::smooth_first(i, dinv, b, c2)
        });
        self.x[l].copy_from_slice(&self.d[l]);
    }

    fn smooth_step(&mut self, pool: &Workers, l: usize, c1: f64, c2: f64) {
        let (c1, c2) = (T::from_f64(c1), T::from_f64(c2));
        let lv = &self.levels[l];
        let (a, dinv, b, x) = (lv.a.view(), &lv.dinv, &self.b[l], &self.x[l]);
        // SAFETY: A_l was built by `SellMatrix::from_csr` for these sizes.
        par_fill(pool, &mut self.d[l], |i, di| unsafe {
            amg::smooth_step(i, &a, dinv, b, x, di, c1, c2)
        });
        let d = &self.d[l];
        par_fill(pool, &mut self.x[l], |i, xi| amg::add(xi, d[i]));
    }

    fn restrict(&mut self, pool: &Workers, l: usize) {
        let lv = &self.levels[l];
        let (a, b, x) = (lv.a.view(), &self.b[l], &self.x[l]);
        // SAFETY: as above.
        par_fill(pool, &mut self.w[l], |i, _| unsafe {
            amg::residual(i, &a, b, x)
        });
        let (r, w) = (lv.r.view(), &self.w[l]);
        // SAFETY: R_l has the rows of level l + 1 and the columns of level l.
        par_fill(pool, &mut self.b[l + 1], |i, _| unsafe { r.row(i, w) });
    }

    fn coarse(&mut self, pool: &Workers) {
        let nl = self.levels.len();
        let (n, m, b) = (self.coarse_n, &self.coarse_inv, &self.b[nl]);
        // SAFETY: the dense inverse is n × n.
        par_fill(pool, &mut self.x[nl], |i, _| unsafe {
            amg::dense_row(i, n, m, b)
        });
    }

    fn prolong(&mut self, pool: &Workers, l: usize) {
        let (lo, hi) = self.x.split_at_mut(l + 1);
        let (p, xc) = (self.levels[l].p.view(), &hi[0]);
        // SAFETY: P_l has the rows of level l and the columns of level l + 1.
        par_fill(pool, &mut lo[l], |i, xi| unsafe {
            amg::prolong_add(i, &p, xc, xi)
        });
    }
}

macro_rules! with_amg {
    ($self:ident, $a:ident => $e:expr) => {
        match $self.amg.as_mut() {
            Some(Amg::F32($a)) => $e,
            Some(Amg::F64($a)) => $e,
            None => anyhow::bail!("no multigrid hierarchy loaded"),
        }
    };
}

/// A raw pointer that may be shared by the threads of a gather loop in
/// which every index writes its own, disjoint outputs.
#[derive(Clone, Copy)]
struct Shared(*mut f64);
// SAFETY: used only where each thread writes disjoint entries.
unsafe impl Send for Shared {}
// SAFETY: as above.
unsafe impl Sync for Shared {}

impl CpuBackend {
    /// New backend on `threads` threads (`None`: all cores).
    pub fn new(d: Discretisation, threads: Option<usize>) -> anyhow::Result<Self> {
        Self::build(d, threads, None)
    }

    /// New backend for one rank of a decomposed run (MPI step W3): `d` is
    /// the discretisation of the rank's local mesh, `ranks` its link to the
    /// others. Cell vectors carry the halo cells after the own ones; an
    /// operation that reads a vector across faces first fills its halo
    /// from the owners (φ before the gradient, the gradient and ψ before
    /// the face coefficients, x before a matrix product), dot products and
    /// extrema are exact over all ranks. So every rank computes bitwise the
    /// values of its cells in the one-rank run.
    pub fn on_ranks(
        d: Discretisation,
        threads: Option<usize>,
        ranks: Ranks,
    ) -> anyhow::Result<Self> {
        let halo: usize = ranks.halo.neighbours.iter().map(|n| n.recv.len()).sum();
        anyhow::ensure!(
            halo == d.n_halo,
            "the halo plan has {halo} halo cells, the discretisation {}",
            d.n_halo
        );
        Self::build(d, threads, Some(ranks))
    }

    fn build(
        d: Discretisation,
        threads: Option<usize>,
        ranks: Option<Ranks>,
    ) -> anyhow::Result<Self> {
        let pool = Workers::new(threads)?;
        let mut vecs = vec![Vec::new(); V::COUNT];
        for v in V::CELL {
            vecs[v.index()] = vec![0.0; d.n_cells + d.n_halo];
        }
        for v in V::FACE {
            vecs[v.index()] = vec![0.0; d.n_faces];
        }
        vecs[V::Psi.index()].fill(1.0);
        let val = vec![0.0; d.sell.stored()];
        Ok(Self {
            d,
            pool,
            vecs,
            scratch: Vec::new(),
            val,
            eps2: None,
            amg: None,
            parked: Vec::new(),
            slot: 0,
            ranks,
        })
    }

    fn v(&self, v: V) -> &[f64] {
        &self.vecs[v.index()]
    }

    /// Fills the halo entries of the cell vectors `vs` from the owning
    /// ranks (nothing on one partition).
    fn exchange(&mut self, vs: &[V]) {
        let Some(r) = &self.ranks else { return };
        let mut taken: Vec<Vec<f64>> = vs
            .iter()
            .map(|v| std::mem::take(&mut self.vecs[v.index()]))
            .collect();
        {
            let mut refs: Vec<&mut [f64]> = taken.iter_mut().map(|x| x.as_mut_slice()).collect();
            r.exchange(&mut refs);
        }
        for (v, x) in vs.iter().zip(taken) {
            self.vecs[v.index()] = x;
        }
    }

    /// Computes `out[i] = body(self, i)` for the own cells into the scratch
    /// buffer and swaps it in, so `out` may also be an input; the halo
    /// entries of `out` are kept.
    fn map_into(&mut self, out: V, body: impl Fn(&Self, usize) -> f64 + Sync) {
        let len = self.vecs[out.index()].len();
        let n = self.d.n_cells;
        let mut s = std::mem::take(&mut self.scratch);
        s.resize(len, 0.0);
        s[n..].copy_from_slice(&self.vecs[out.index()][n..]);
        let s_all = &mut s;
        let s = &mut s_all[..n];
        let this = &*self;
        if this.pool.inline(n) {
            for (i, o) in s.iter_mut().enumerate() {
                *o = body(this, i);
            }
        } else {
            this.pool.install(|| {
                s.par_chunks_mut(GRAIN).enumerate().for_each(|(c, ch)| {
                    for (j, o) in ch.iter_mut().enumerate() {
                        *o = body(this, c * GRAIN + j);
                    }
                })
            });
        }
        let mut s = std::mem::take(s_all);
        std::mem::swap(&mut s, &mut self.vecs[out.index()]);
        self.scratch = s;
    }

    fn grad_view(&self) -> Vec3<'_, f64> {
        Vec3 {
            x: self.v(V::Gx),
            y: self.v(V::Gy),
            z: self.v(V::Gz),
        }
    }

    fn sell_row(&self, r: usize, x: &[f64]) -> f64 {
        let s = &self.d.sell;
        // SAFETY: the structure was built by `sell_structure` for this
        // mesh; every column is a cell index.
        unsafe { spmv::sell_row(r, s.chunk, &s.slice_ptr, &s.slice_len, &s.col, &self.val, x) }
    }
}

impl Backend for CpuBackend {
    fn name(&self) -> String {
        if self.pool.threads() == 1 {
            "cpu-serial".into()
        } else {
            "cpu".into()
        }
    }

    fn device(&self) -> String {
        format!("CPU, {} threads", self.pool.threads())
    }

    fn n_cells(&self) -> usize {
        self.d.n_cells
    }

    fn upload(&mut self, v: V, data: &[f64]) -> anyhow::Result<()> {
        // a cell vector takes the own cells (the halo is filled by the
        // exchange before it is read) or the own and halo cells
        let own = self.d.n_cells;
        let t = &mut self.vecs[v.index()];
        let is_cell = V::CELL.contains(&v);
        anyhow::ensure!(
            t.len() == data.len() || (is_cell && data.len() == own),
            "upload of {v:?}: wrong length"
        );
        t[..data.len()].copy_from_slice(data);
        Ok(())
    }

    fn download(&mut self, v: V) -> anyhow::Result<Vec<f64>> {
        // a cell vector without its halo
        let x = &self.vecs[v.index()];
        Ok(if V::CELL.contains(&v) {
            x[..self.d.n_cells].to_vec()
        } else {
            x.clone()
        })
    }

    fn sync(&mut self) -> anyhow::Result<()> {
        Ok(())
    }

    fn gradient(&mut self) -> anyhow::Result<()> {
        self.exchange(&[V::Phi]);
        let nc = self.d.n_cells;
        let mut gx_all = std::mem::take(&mut self.vecs[V::Gx.index()]);
        let mut gy_all = std::mem::take(&mut self.vecs[V::Gy.index()]);
        let mut gz_all = std::mem::take(&mut self.vecs[V::Gz.index()]);
        let (gx, gy, gz) = (&mut gx_all[..nc], &mut gy_all[..nc], &mut gz_all[..nc]);
        let (d, phi) = (&self.d, &self.vecs[V::Phi.index()]);
        let (topo, bc, w) = (d.topology(), d.boundary(), d.lsq.view());
        let fill = |first: usize, x: &mut [f64], y: &mut [f64], z: &mut [f64]| {
            for j in 0..x.len() {
                // SAFETY: the discretisation was built from a validated
                // mesh; first + j is a local cell.
                let g = unsafe { fv::lsq_gradient(first + j, &topo, &w, phi, &bc) };
                x[j] = g[0];
                y[j] = g[1];
                z[j] = g[2];
            }
        };
        if self.pool.inline(gx.len()) {
            fill(0, gx, gy, gz);
        } else {
            self.pool.install(|| {
                gx.par_chunks_mut(GRAIN)
                    .zip(gy.par_chunks_mut(GRAIN))
                    .zip(gz.par_chunks_mut(GRAIN))
                    .enumerate()
                    .for_each(|(ci, ((x, y), z))| fill(ci * GRAIN, x, y, z))
            });
        }
        self.vecs[V::Gx.index()] = gx_all;
        self.vecs[V::Gy.index()] = gy_all;
        self.vecs[V::Gz.index()] = gz_all;
        // the face coefficients read the gradient of both cells of a face
        self.exchange(&[V::Gx, V::Gy, V::Gz]);
        Ok(())
    }

    fn limiter(&mut self, l: Limiter) -> anyhow::Result<()> {
        let eps2: Vec<f64> = match l {
            Limiter::None => {
                self.vecs[V::Psi.index()].fill(1.0);
                return Ok(());
            }
            Limiter::Venkatakrishnan { k } => match &self.eps2 {
                Some((kk, e)) if *kk == k => e.clone(),
                _ => {
                    let e = self.d.venkat_eps2(k);
                    self.eps2 = Some((k, e.clone()));
                    e
                }
            },
            Limiter::VenkatakrishnanWang { k } => {
                let (mn, mx) = self.min_max(V::Phi)?;
                let e = k * (mx - mn);
                vec![e * e]
            }
            // a negative single ε² selects Barth–Jespersen
            Limiter::BarthJespersen => vec![-1.0],
        };
        self.map_into(V::Psi, |s, c| {
            let d = &s.d;
            // SAFETY: as in `gradient`; eps2 has one entry or one per cell.
            unsafe {
                fv::venkat_limiter(
                    c,
                    &d.topology(),
                    s.v(V::Phi),
                    &d.boundary(),
                    &s.grad_view(),
                    &d.cell_centroid.view(),
                    &d.face_centroid.view(),
                    &eps2,
                )
            }
        });
        // the face coefficients read ψ of the upwind cell
        self.exchange(&[V::Psi]);
        Ok(())
    }

    fn faces(&mut self, st: FaceSettings) -> anyhow::Result<()> {
        let mut cp = std::mem::take(&mut self.vecs[V::Cp.index()]);
        let mut cn = std::mem::take(&mut self.vecs[V::Cn.index()]);
        let mut bf = std::mem::take(&mut self.vecs[V::Bf.index()]);
        let d = &self.d;
        let inp = FaceInputs {
            topo: d.topology(),
            bc: d.boundary(),
            mass_flux: &d.mass_flux,
            gamma: &d.gamma,
            dgeo: &d.dgeo,
            k: d.k.view(),
            weight: &d.weight,
            face_centroid: d.face_centroid.view(),
            cell_centroid: d.cell_centroid.view(),
            shift: d.shift.view(),
            phi: self.v(V::Phi),
            grad: self.grad_view(),
            psi: self.v(V::Psi),
            blend: st.blend,
            central: st.central,
            beta: st.beta,
            nonorth: if st.nonorth { 1.0 } else { 0.0 },
            nonorth_limit: st.nonorth_limit,
            bound: 0.0,
        };
        let fill = |first: usize, a: &mut [f64], b: &mut [f64], c: &mut [f64]| {
            for j in 0..a.len() {
                // SAFETY: first + j is a face of the discretisation.
                let (x, y, z) = unsafe { fv::face_coefficients(first + j, &inp) };
                a[j] = x;
                b[j] = y;
                c[j] = z;
            }
        };
        if self.pool.inline(cp.len()) {
            fill(0, &mut cp, &mut cn, &mut bf);
        } else {
            self.pool.install(|| {
                cp.par_chunks_mut(GRAIN)
                    .zip(cn.par_chunks_mut(GRAIN))
                    .zip(bf.par_chunks_mut(GRAIN))
                    .enumerate()
                    .for_each(|(ci, ((a, b), c))| fill(ci * GRAIN, a, b, c))
            });
        }
        self.vecs[V::Cp.index()] = cp;
        self.vecs[V::Cn.index()] = cn;
        self.vecs[V::Bf.index()] = bf;
        Ok(())
    }

    fn assemble(&mut self, alpha: f64) -> anyhow::Result<()> {
        let mut rhs = std::mem::take(&mut self.vecs[V::Rhs.index()]);
        let mut dinv = std::mem::take(&mut self.vecs[V::Dinv.index()]);
        let (pv, pr, pd) = (
            Shared(self.val.as_mut_ptr()),
            Shared(rhs.as_mut_ptr()),
            Shared(dinv.as_mut_ptr()),
        );
        let d = &self.d;
        let topo = d.topology();
        let layout = SellLayout {
            chunk: d.sell.chunk,
            slice_ptr: &d.sell.slice_ptr,
        };
        let (cp, cn, bf, phi) = (self.v(V::Cp), self.v(V::Cn), self.v(V::Bf), self.v(V::Phi));
        let n = d.n_cells;
        let rows = |ci: usize| {
            let (pv, pr, pd) = (pv, pr, pd);
            for c in ci * GRAIN..((ci + 1) * GRAIN).min(n) {
                // SAFETY: row c and rhs[c], dinv[c] are written by this
                // iteration only; the SELL storage has a slot for every
                // interior face of c.
                unsafe {
                    fv::assemble_row(
                        c, &topo, cp, cn, bf, &d.source, phi, alpha, &layout, pv.0, pr.0, pd.0,
                    )
                }
            }
        };
        if self.pool.inline(n) {
            (0..n.div_ceil(GRAIN)).for_each(rows);
        } else {
            self.pool
                .install(|| (0..n.div_ceil(GRAIN)).into_par_iter().for_each(rows));
        }
        self.vecs[V::Rhs.index()] = rhs;
        self.vecs[V::Dinv.index()] = dinv;
        Ok(())
    }

    fn spmv(&mut self, x: V, out: V) -> anyhow::Result<()> {
        self.exchange(&[x]);
        self.map_into(out, |s, r| s.sell_row(r, s.v(x)));
        Ok(())
    }

    fn residual(&mut self, x: V, b: V, out: V) -> anyhow::Result<()> {
        self.exchange(&[x]);
        self.map_into(out, |s, r| s.v(b)[r] - s.sell_row(r, s.v(x)));
        Ok(())
    }

    fn dot(&mut self, a: V, b: V) -> anyhow::Result<f64> {
        let n = self.d.n_cells;
        let (x, y) = (&self.v(a)[..n], &self.v(b)[..n]);
        Ok(match &self.ranks {
            None => exact_dot(&self.pool, x, y),
            Some(r) => r.comm.sum_exact(&exact_dot_acc(&self.pool, x, y)),
        })
    }

    fn min_max(&mut self, a: V) -> anyhow::Result<(f64, f64)> {
        let (lo, hi) = min_max(&self.pool, &self.v(a)[..self.d.n_cells]);
        Ok(match &self.ranks {
            None => (lo, hi),
            Some(r) => (r.comm.min(lo), r.comm.max(hi)),
        })
    }

    fn update_p(&mut self, beta: f64, omega: f64) -> anyhow::Result<()> {
        // SAFETY: all cell vectors have the same length.
        self.map_into(V::P, |s, i| unsafe {
            blas::update_p(i, s.v(V::R), s.v(V::P), s.v(V::Vv), beta, omega)
        });
        Ok(())
    }

    fn sub_scaled(&mut self, x: V, z: V, alpha: f64, out: V) -> anyhow::Result<()> {
        // SAFETY: as above.
        self.map_into(out, |s, i| unsafe {
            blas::sub_scaled(i, s.v(x), s.v(z), alpha)
        });
        Ok(())
    }

    fn update_x(&mut self, alpha: f64, omega: f64) -> anyhow::Result<()> {
        // SAFETY: as above.
        self.map_into(V::Phi, |s, i| unsafe {
            blas::update_x(i, s.v(V::Phi), s.v(V::Ph), s.v(V::Sh), alpha, omega)
        });
        Ok(())
    }

    fn precondition(&mut self, x: V, out: V) -> anyhow::Result<()> {
        // SAFETY: as above.
        self.map_into(out, |s, i| unsafe { blas::scale(i, s.v(V::Dinv), s.v(x)) });
        Ok(())
    }

    fn copy(&mut self, src: V, dst: V) -> anyhow::Result<()> {
        let s = self.vecs[src.index()].clone();
        self.vecs[dst.index()] = s;
        Ok(())
    }

    fn axpy(&mut self, alpha: f64, x: V, y: V) -> anyhow::Result<()> {
        // SAFETY: all cell vectors have the same length.
        self.map_into(y, |s, i| unsafe { blas::axpy(i, alpha, s.v(x), s.v(y)) });
        Ok(())
    }

    fn xpby(&mut self, x: V, beta: f64, y: V) -> anyhow::Result<()> {
        // SAFETY: as above.
        self.map_into(y, |s, i| unsafe { blas::xpby(i, s.v(x), beta, s.v(y)) });
        Ok(())
    }

    fn matrix_values(&mut self) -> anyhow::Result<Vec<f64>> {
        Ok(self.val.clone())
    }

    fn set_matrix(&mut self, val: &[f64]) -> anyhow::Result<()> {
        anyhow::ensure!(val.len() == self.val.len(), "set_matrix: wrong length");
        self.val.copy_from_slice(val);
        let dinv = diagonal_inverse(&self.d, val);
        self.vecs[V::Dinv.index()][..dinv.len()].copy_from_slice(&dinv);
        Ok(())
    }

    fn discretisation(&self) -> Option<&Discretisation> {
        Some(&self.d)
    }

    fn amg_slot(&mut self, slot: usize) -> anyhow::Result<()> {
        if slot != self.slot {
            let n = slot.max(self.slot) + 1;
            if self.parked.len() < n {
                self.parked.resize_with(n, || None);
            }
            self.parked[self.slot] = self.amg.take();
            self.amg = self.parked[slot].take();
            self.slot = slot;
        }
        Ok(())
    }

    fn amg_load(&mut self, h: &Hierarchy) -> anyhow::Result<()> {
        anyhow::ensure!(
            self.ranks.is_none(),
            "the multigrid preconditioner runs on one rank (the distributed multigrid is MPI step W4); use the Jacobi preconditioner on ranks"
        );
        anyhow::ensure!(
            h.levels.first().map_or(h.coarse_n, |l| l.a.n_rows) == self.d.n_cells,
            "hierarchy size does not match the mesh"
        );
        self.amg = Some(match h.precision {
            Precision::Mixed => Amg::F32(AmgData::new(h)),
            Precision::Double => Amg::F64(AmgData::new(h)),
        });
        Ok(())
    }
}

impl LevelOps for CpuBackend {
    fn amg_levels(&self) -> usize {
        match &self.amg {
            Some(Amg::F32(a)) => a.levels.len(),
            Some(Amg::F64(a)) => a.levels.len(),
            None => 0,
        }
    }

    fn amg_steps(&self, l: usize) -> Vec<(f64, f64)> {
        match &self.amg {
            Some(Amg::F32(a)) => a.levels[l].steps.clone(),
            Some(Amg::F64(a)) => a.levels[l].steps.clone(),
            None => Vec::new(),
        }
    }

    fn amg_load_rhs(&mut self, r: V) -> anyhow::Result<()> {
        let rv = std::mem::take(&mut self.vecs[r.index()]);
        let pool = &self.pool;
        let res = (|| -> anyhow::Result<()> {
            with_amg!(self, a => a.load_rhs(pool, &rv));
            Ok(())
        })();
        self.vecs[r.index()] = rv;
        res
    }

    fn amg_store(&mut self, z: V) -> anyhow::Result<()> {
        let mut zv = std::mem::take(&mut self.vecs[z.index()]);
        let pool = &self.pool;
        let res = (|| -> anyhow::Result<()> {
            with_amg!(self, a => a.store(pool, &mut zv));
            Ok(())
        })();
        self.vecs[z.index()] = zv;
        res
    }

    fn amg_smooth_first(&mut self, l: usize, c2: f64) -> anyhow::Result<()> {
        let pool = &self.pool;
        with_amg!(self, a => a.smooth_first(pool, l, c2));
        Ok(())
    }

    fn amg_smooth_step(&mut self, l: usize, c1: f64, c2: f64) -> anyhow::Result<()> {
        let pool = &self.pool;
        with_amg!(self, a => a.smooth_step(pool, l, c1, c2));
        Ok(())
    }

    fn amg_restrict(&mut self, l: usize) -> anyhow::Result<()> {
        let pool = &self.pool;
        with_amg!(self, a => a.restrict(pool, l));
        Ok(())
    }

    fn amg_coarse(&mut self) -> anyhow::Result<()> {
        let pool = &self.pool;
        with_amg!(self, a => a.coarse(pool));
        Ok(())
    }

    fn amg_prolong(&mut self, l: usize) -> anyhow::Result<()> {
        let pool = &self.pool;
        with_amg!(self, a => a.prolong(pool, l));
        Ok(())
    }
}

/// 1 / diagonal of a matrix in the sliced-ELL slots of `d` (slot 0).
pub fn diagonal_inverse(d: &Discretisation, val: &[f64]) -> Vec<f64> {
    (0..d.n_cells)
        .map(|c| 1.0 / val[d.sell.pos(c, 0)])
        .collect()
}
