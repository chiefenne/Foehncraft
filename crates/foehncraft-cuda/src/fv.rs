//! NVIDIA GPU backend of `foehncraft-fv` (feature `cuda`): the kernel bodies of
//! `foehncraft_kernels` compiled by rustc to PTX (`gpu/ptx-kernels`, entry
//! points `fv_*` and `blas_*`), launched through the CUDA driver. All data
//! stay on the GPU; only dot products (one value per block of 256) and the
//! limiter range come back to the host.
//!
//! The PTX has no fused multiply-adds and IEEE-rounded division. The dot
//! products are block sums in a fixed tree; the CPU backend sums exactly
//! since MPI step W2, so CPU and GPU results agree to rounding, not bitwise,
//! until the GPU sums exactly too (v1.1, owner decision of 2026-11-23;
//! `foehn scalar` with `backends = ["cpu", "cuda"]` reports the difference).

use anyhow::Context as _;
use cudarc::driver::{
    CudaContext, CudaFunction, CudaModule, CudaSlice, CudaStream, DeviceRepr, LaunchConfig,
    PushKernelArg, ValidAsZeroBits,
};
use cudarc::nvrtc::Ptx;
use foehncraft_core::Real;
use foehncraft_core::repro::{BLOCK, combine_block_sums};
use foehncraft_fv::amg::cycle::LevelOps;
use foehncraft_fv::amg::{Hierarchy, Precision, SellMatrix};
use foehncraft_fv::backend::{Backend, FaceSettings, Limiter, V};
use foehncraft_fv::disc::Discretisation;
use std::path::Path;
use std::sync::Arc;

/// Threads per block.
const THREADS: u32 = 256;

const ENTRIES: [&str; 14] = [
    "fv_gradient",
    "fv_limiter",
    "fv_faces",
    "fv_assemble",
    "fv_residual",
    "g1_spmv_sell",
    "blas_block_dot",
    "blas_block_min_max",
    "blas_update_p",
    "blas_sub_scaled",
    "blas_update_x",
    "blas_scale",
    "blas_axpy",
    "blas_xpby",
];

/// Entry points of the V-cycle, per precision (suffix `_f32`, `_f64`).
const AMG_ENTRIES: [&str; 9] = [
    "amg_smooth_first",
    "amg_smooth_step",
    "amg_add",
    "amg_residual",
    "amg_spmv",
    "amg_prolong",
    "amg_dense",
    "amg_to",
    "amg_from",
];
const SMOOTH_FIRST: usize = 0;
const SMOOTH_STEP: usize = 1;
const ADD: usize = 2;
const RESIDUAL: usize = 3;
const SPMV: usize = 4;
const PROLONG: usize = 5;
const DENSE: usize = 6;
const TO: usize = 7;
const FROM: usize = 8;

/// A precision of the hierarchy on the GPU.
trait DevReal: Real + DeviceRepr + ValidAsZeroBits {}
impl DevReal for f32 {}
impl DevReal for f64 {}

/// A sliced-ELL matrix on the GPU.
struct DevSell<T> {
    n_rows: u32,
    n_cols: u32,
    chunk: u32,
    n_slices: u32,
    stored: u32,
    slice_ptr: CudaSlice<u32>,
    slice_len: CudaSlice<u32>,
    col: CudaSlice<u32>,
    val: CudaSlice<T>,
}

impl<T: DevReal> DevSell<T> {
    fn new(s: &Arc<CudaStream>, m: &SellMatrix<T>) -> anyhow::Result<Self> {
        Ok(Self {
            n_rows: m.n_rows as u32,
            n_cols: m.n_cols as u32,
            chunk: m.chunk as u32,
            n_slices: m.slice_ptr.len() as u32,
            stored: m.stored() as u32,
            slice_ptr: s.clone_htod(&m.slice_ptr)?,
            slice_len: s.clone_htod(&m.slice_len)?,
            col: s.clone_htod(&m.col)?,
            val: s.clone_htod(&m.val)?,
        })
    }
}

struct DevLevel<T> {
    a: DevSell<T>,
    p: DevSell<T>,
    r: DevSell<T>,
    dinv: CudaSlice<T>,
    steps: Vec<(f64, f64)>,
}

/// The hierarchy on the GPU with its cycle vectors (`x`, `b`: L + 1 levels,
/// `w`, `d`: L levels) and its entry points.
struct DevAmg<T> {
    levels: Vec<DevLevel<T>>,
    coarse_n: u32,
    coarse_inv: CudaSlice<T>,
    x: Vec<CudaSlice<T>>,
    b: Vec<CudaSlice<T>>,
    w: Vec<CudaSlice<T>>,
    d: Vec<CudaSlice<T>>,
    f: Vec<CudaFunction>,
}

enum AmgGpu {
    F32(DevAmg<f32>),
    F64(DevAmg<f64>),
}

impl<T: DevReal> DevAmg<T> {
    fn new(s: &Arc<CudaStream>, module: &Arc<CudaModule>, h: &Hierarchy) -> anyhow::Result<Self> {
        let mut f = Vec::new();
        for e in AMG_ENTRIES {
            let name = format!("{e}_{}", T::NAME);
            f.push(
                module
                    .load_function(&name)
                    .with_context(|| format!("entry point {name} (rebuild the PTX module)"))?,
            );
        }
        let cv = |v: &[f64]| v.iter().map(|&x| T::from_f64(x)).collect::<Vec<T>>();
        let mut levels = Vec::new();
        for l in &h.levels {
            levels.push(DevLevel {
                a: DevSell::new(s, &SellMatrix::<T>::from_csr(&l.a))?,
                p: DevSell::new(s, &SellMatrix::<T>::from_csr(&l.p))?,
                r: DevSell::new(s, &SellMatrix::<T>::from_csr(&l.r))?,
                dinv: s.clone_htod(&cv(&l.dinv))?,
                steps: l.steps.clone(),
            });
        }
        let mut sizes: Vec<usize> = h.levels.iter().map(|l| l.a.n_rows).collect();
        sizes.push(h.coarse_n);
        let alloc = |range: &[usize]| -> anyhow::Result<Vec<CudaSlice<T>>> {
            range.iter().map(|&n| Ok(s.alloc_zeros::<T>(n)?)).collect()
        };
        let nl = levels.len();
        Ok(Self {
            x: alloc(&sizes)?,
            b: alloc(&sizes)?,
            w: alloc(&sizes[..nl])?,
            d: alloc(&sizes[..nl])?,
            coarse_n: h.coarse_n as u32,
            coarse_inv: s.clone_htod(&cv(&h.coarse_inv))?,
            levels,
            f,
        })
    }

    fn load_rhs(&mut self, s: &Arc<CudaStream>, r: &CudaSlice<f64>) -> anyhow::Result<()> {
        let n = self.b[0].len() as u32;
        let mut l = s.launch_builder(&self.f[TO]);
        l.arg(&n).arg(r).arg(&mut self.b[0]);
        // SAFETY: argument list of amg_to_*; both vectors have n entries.
        unsafe { l.launch(cfg(n as usize)) }?;
        Ok(())
    }

    fn store(&mut self, s: &Arc<CudaStream>, z: &mut CudaSlice<f64>) -> anyhow::Result<()> {
        let n = self.x[0].len() as u32;
        let mut l = s.launch_builder(&self.f[FROM]);
        l.arg(&n).arg(&self.x[0]).arg(z);
        // SAFETY: argument list of amg_from_*.
        unsafe { l.launch(cfg(n as usize)) }?;
        Ok(())
    }

    fn smooth_first(&mut self, s: &Arc<CudaStream>, lv: usize, c2: f64) -> anyhow::Result<()> {
        let n = self.levels[lv].a.n_rows;
        let c2 = T::from_f64(c2);
        let mut l = s.launch_builder(&self.f[SMOOTH_FIRST]);
        l.arg(&n)
            .arg(&c2)
            .arg(&self.levels[lv].dinv)
            .arg(&self.b[lv])
            .arg(&mut self.d[lv])
            .arg(&mut self.x[lv]);
        // SAFETY: argument list of amg_smooth_first_*.
        unsafe { l.launch(cfg(n as usize)) }?;
        Ok(())
    }

    fn smooth_step(
        &mut self,
        s: &Arc<CudaStream>,
        lv: usize,
        c1: f64,
        c2: f64,
    ) -> anyhow::Result<()> {
        let (c1, c2) = (T::from_f64(c1), T::from_f64(c2));
        let a = &self.levels[lv].a;
        let n = a.n_rows;
        {
            let mut l = s.launch_builder(&self.f[SMOOTH_STEP]);
            l.arg(&n).arg(&a.chunk).arg(&a.n_slices).arg(&a.stored);
            l.arg(&a.slice_ptr)
                .arg(&a.slice_len)
                .arg(&a.col)
                .arg(&a.val);
            l.arg(&c1).arg(&c2);
            l.arg(&self.levels[lv].dinv)
                .arg(&self.b[lv])
                .arg(&self.x[lv])
                .arg(&mut self.d[lv]);
            // SAFETY: argument list of amg_smooth_step_*.
            unsafe { l.launch(cfg(n as usize)) }?;
        }
        let mut l = s.launch_builder(&self.f[ADD]);
        l.arg(&n).arg(&self.d[lv]).arg(&mut self.x[lv]);
        // SAFETY: argument list of amg_add_*.
        unsafe { l.launch(cfg(n as usize)) }?;
        Ok(())
    }

    fn restrict(&mut self, s: &Arc<CudaStream>, lv: usize) -> anyhow::Result<()> {
        let a = &self.levels[lv].a;
        {
            let n = a.n_rows;
            let mut l = s.launch_builder(&self.f[RESIDUAL]);
            l.arg(&n).arg(&a.chunk).arg(&a.n_slices).arg(&a.stored);
            l.arg(&a.slice_ptr)
                .arg(&a.slice_len)
                .arg(&a.col)
                .arg(&a.val);
            l.arg(&self.b[lv]).arg(&self.x[lv]).arg(&mut self.w[lv]);
            // SAFETY: argument list of amg_residual_*.
            unsafe { l.launch(cfg(n as usize)) }?;
        }
        let r = &self.levels[lv].r;
        let next = &mut self.b[lv + 1];
        let mut l = s.launch_builder(&self.f[SPMV]);
        l.arg(&r.n_rows).arg(&r.n_cols);
        l.arg(&r.chunk).arg(&r.n_slices).arg(&r.stored);
        l.arg(&r.slice_ptr)
            .arg(&r.slice_len)
            .arg(&r.col)
            .arg(&r.val);
        l.arg(&self.w[lv]).arg(next);
        // SAFETY: argument list of amg_spmv_*; R_l maps level l to l + 1.
        unsafe { l.launch(cfg(r.n_rows as usize)) }?;
        Ok(())
    }

    fn coarse(&mut self, s: &Arc<CudaStream>) -> anyhow::Result<()> {
        let nl = self.levels.len();
        let n = self.coarse_n;
        let mut l = s.launch_builder(&self.f[DENSE]);
        l.arg(&n)
            .arg(&self.coarse_inv)
            .arg(&self.b[nl])
            .arg(&mut self.x[nl]);
        // SAFETY: argument list of amg_dense_*; n × n inverse.
        unsafe { l.launch(cfg(n as usize)) }?;
        Ok(())
    }

    fn prolong(&mut self, s: &Arc<CudaStream>, lv: usize) -> anyhow::Result<()> {
        let p = &self.levels[lv].p;
        let (lo, hi) = self.x.split_at_mut(lv + 1);
        let mut l = s.launch_builder(&self.f[PROLONG]);
        l.arg(&p.n_rows).arg(&p.n_cols);
        l.arg(&p.chunk).arg(&p.n_slices).arg(&p.stored);
        l.arg(&p.slice_ptr)
            .arg(&p.slice_len)
            .arg(&p.col)
            .arg(&p.val);
        l.arg(&hi[0]).arg(&mut lo[lv]);
        // SAFETY: argument list of amg_prolong_*; P_l maps level l + 1 to l.
        unsafe { l.launch(cfg(p.n_rows as usize)) }?;
        Ok(())
    }
}

macro_rules! with_amg {
    ($self:ident, $a:ident => $e:expr) => {
        match $self.amg.as_mut() {
            Some(AmgGpu::F32($a)) => $e,
            Some(AmgGpu::F64($a)) => $e,
            None => anyhow::bail!("no multigrid hierarchy loaded"),
        }
    };
}

fn cfg(n: usize) -> LaunchConfig {
    let n = n as u32;
    LaunchConfig {
        grid_dim: (n.div_ceil(THREADS).max(1), 1, 1),
        block_dim: (THREADS, 1, 1),
        shared_mem_bytes: 0,
    }
}

struct Soa {
    x: CudaSlice<f64>,
    y: CudaSlice<f64>,
    z: CudaSlice<f64>,
}

/// Device copy of the discretisation.
struct Dev {
    owner: CudaSlice<u32>,
    neighbour: CudaSlice<u32>,
    cell_ptr: CudaSlice<u32>,
    cell_faces: CudaSlice<u32>,
    cc: Soa,
    fc: Soa,
    k: Soa,
    lsq: Soa,
    dgeo: CudaSlice<f64>,
    weight: CudaSlice<f64>,
    bc_a: CudaSlice<f64>,
    bc_b: CudaSlice<f64>,
    mass_flux: CudaSlice<f64>,
    gamma: CudaSlice<f64>,
    source: CudaSlice<f64>,
    slice_ptr: CudaSlice<u32>,
    slice_len: CudaSlice<u32>,
    col: CudaSlice<u32>,
    val: CudaSlice<f64>,
    eps2_cell: CudaSlice<f64>,
    eps2_one: CudaSlice<f64>,
    sums: CudaSlice<f64>,
    sums2: CudaSlice<f64>,
}

/// CUDA backend.
pub struct CudaBackend {
    ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    module: Arc<CudaModule>,
    f: Vec<CudaFunction>,
    d: Dev,
    vecs: Vec<CudaSlice<f64>>,
    nc: u32,
    nf: u32,
    ni: u32,
    ne: u32,
    chunk: u32,
    n_slices: u32,
    stored: u32,
    venkat_k: Option<f64>,
    host: Discretisation,
    amg: Option<AmgGpu>,
    /// The hierarchies of the inactive slots, and the active slot.
    parked: Vec<Option<AmgGpu>>,
    slot: usize,
}

impl CudaBackend {
    /// Uploads the discretisation to GPU `ordinal` and loads the PTX module.
    pub fn new(disc: &Discretisation, ordinal: usize, ptx: &Path) -> anyhow::Result<Self> {
        anyhow::ensure!(
            ptx.exists(),
            "PTX module `{}` not found; build it with python3 gpu/ptx-kernels/build.py",
            ptx.display()
        );
        let ctx =
            CudaContext::new(ordinal).context("CUDA device (is the NVIDIA driver installed?)")?;
        let src = std::fs::read_to_string(ptx)?;
        let module = ctx
            .load_module(Ptx::from_src(src))
            .context("loading the PTX module")?;
        let mut f = Vec::new();
        for e in ENTRIES {
            f.push(
                module
                    .load_function(e)
                    .with_context(|| format!("entry point {e} (rebuild the PTX module)"))?,
            );
        }
        let stream = ctx.default_stream();
        let s = &stream;
        let soa = |v: &foehncraft_fv::disc::Soa3| -> anyhow::Result<Soa> {
            Ok(Soa {
                x: s.clone_htod(&v.x)?,
                y: s.clone_htod(&v.y)?,
                z: s.clone_htod(&v.z)?,
            })
        };
        let nc = disc.n_cells;
        let nb = nc.div_ceil(BLOCK);
        let d = Dev {
            owner: s.clone_htod(&disc.owner)?,
            neighbour: s.clone_htod(&nonempty_u32(&disc.neighbour))?,
            cell_ptr: s.clone_htod(&disc.cell_ptr)?,
            cell_faces: s.clone_htod(&disc.cell_faces)?,
            cc: soa(&disc.cell_centroid)?,
            fc: soa(&disc.face_centroid)?,
            k: soa(&disc.k)?,
            lsq: soa(&disc.lsq)?,
            dgeo: s.clone_htod(&disc.dgeo)?,
            weight: s.clone_htod(&disc.weight)?,
            bc_a: s.clone_htod(&nonempty_f64(&disc.bc_a))?,
            bc_b: s.clone_htod(&nonempty_f64(&disc.bc_b))?,
            mass_flux: s.clone_htod(&disc.mass_flux)?,
            gamma: s.clone_htod(&disc.gamma)?,
            source: s.clone_htod(&disc.source)?,
            slice_ptr: s.clone_htod(&disc.sell.slice_ptr)?,
            slice_len: s.clone_htod(&disc.sell.slice_len)?,
            col: s.clone_htod(&disc.sell.col)?,
            val: s.alloc_zeros(disc.sell.stored())?,
            eps2_cell: s.alloc_zeros(nc)?,
            eps2_one: s.alloc_zeros(1)?,
            sums: s.alloc_zeros(nb)?,
            sums2: s.alloc_zeros(nb)?,
        };
        let mut vecs = Vec::with_capacity(V::COUNT);
        for i in 0..V::COUNT {
            let n = if V::FACE.iter().any(|v| v.index() == i) {
                disc.n_faces
            } else {
                nc
            };
            vecs.push(s.alloc_zeros(n)?);
        }
        let mut be = Self {
            ctx,
            stream,
            module,
            f,
            d,
            vecs,
            nc: nc as u32,
            nf: disc.n_faces as u32,
            ni: disc.n_interior as u32,
            ne: disc.cell_faces.len() as u32,
            chunk: disc.sell.chunk as u32,
            n_slices: disc.sell.slice_ptr.len() as u32,
            stored: disc.sell.stored() as u32,
            venkat_k: None,
            host: disc.clone(),
            amg: None,
            parked: Vec::new(),
            slot: 0,
        };
        be.upload(V::Psi, &vec![1.0; nc])?;
        be.stream.synchronize()?;
        Ok(be)
    }

    /// The stream all work of the backend runs on.
    pub(crate) fn stream(&self) -> &Arc<CudaStream> {
        &self.stream
    }

    /// A named vector on the device.
    pub(crate) fn vector(&self, v: V) -> &CudaSlice<f64> {
        &self.vecs[v.index()]
    }

    /// A named vector on the device, for writing.
    pub(crate) fn vector_mut(&mut self, v: V) -> &mut CudaSlice<f64> {
        &mut self.vecs[v.index()]
    }

    /// The matrix values (sliced-ELL slots of the discretisation), for
    /// writing. The caller also sets `Dinv`.
    pub(crate) fn matrix_mut(&mut self) -> &mut CudaSlice<f64> {
        &mut self.d.val
    }

    fn func(&self, name: &str) -> &CudaFunction {
        let i = ENTRIES.iter().position(|e| *e == name).expect("entry");
        &self.f[i]
    }

    /// Two vectors at once, the second mutable (distinct).
    fn pair(&mut self, a: V, out: V) -> (&CudaSlice<f64>, &mut CudaSlice<f64>) {
        assert_ne!(a, out);
        let (i, j) = (a.index(), out.index());
        if i < j {
            let (l, r) = self.vecs.split_at_mut(j);
            (&l[i], &mut r[0])
        } else {
            let (l, r) = self.vecs.split_at_mut(i);
            (&r[0], &mut l[j])
        }
    }

    /// Takes a vector out of the table (put it back with `put`).
    fn take(&mut self, v: V) -> anyhow::Result<CudaSlice<f64>> {
        let empty = self.stream.alloc_zeros(1)?;
        Ok(std::mem::replace(&mut self.vecs[v.index()], empty))
    }

    fn put(&mut self, v: V, s: CudaSlice<f64>) {
        self.vecs[v.index()] = s;
    }

    fn block_dot(&mut self, a: V, b: V) -> anyhow::Result<f64> {
        let nb = (self.nc as usize).div_ceil(BLOCK);
        let nbu = nb as u32;
        let s = self.stream.clone();
        let f = self.func("blas_block_dot").clone();
        let (va, vb) = (&self.vecs[a.index()], &self.vecs[b.index()]);
        let mut l = s.launch_builder(&f);
        l.arg(&self.nc)
            .arg(&nbu)
            .arg(va)
            .arg(vb)
            .arg(&mut self.d.sums);
        // SAFETY: argument list of blas_block_dot; arrays of length n_cells
        // and n_blocks.
        unsafe { l.launch(cfg(nb)) }?;
        let sums = s.clone_dtoh(&self.d.sums)?;
        Ok(combine_block_sums(&sums))
    }
}

fn nonempty_u32(v: &[u32]) -> Vec<u32> {
    if v.is_empty() { vec![0] } else { v.to_vec() }
}

fn nonempty_f64(v: &[f64]) -> Vec<f64> {
    if v.is_empty() { vec![0.0] } else { v.to_vec() }
}

impl Backend for CudaBackend {
    fn name(&self) -> String {
        "cuda".into()
    }

    fn device(&self) -> String {
        let name = self.ctx.name().unwrap_or_else(|_| "NVIDIA GPU".into());
        let (ma, mi) = self.ctx.compute_capability().unwrap_or((0, 0));
        format!("{name} (sm_{ma}{mi})")
    }

    fn n_cells(&self) -> usize {
        self.nc as usize
    }

    fn upload(&mut self, v: V, data: &[f64]) -> anyhow::Result<()> {
        let s = self.stream.clone();
        let t = &mut self.vecs[v.index()];
        anyhow::ensure!(t.len() == data.len(), "upload of {v:?}: wrong length");
        s.memcpy_htod(data, t)?;
        Ok(())
    }

    fn download(&mut self, v: V) -> anyhow::Result<Vec<f64>> {
        Ok(self.stream.clone_dtoh(&self.vecs[v.index()])?)
    }

    fn sync(&mut self) -> anyhow::Result<()> {
        Ok(self.stream.synchronize()?)
    }

    fn gradient(&mut self) -> anyhow::Result<()> {
        let mut gx = self.take(V::Gx)?;
        let mut gy = self.take(V::Gy)?;
        let mut gz = self.take(V::Gz)?;
        {
            let d = &self.d;
            let f = self.func("fv_gradient");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&self.nf).arg(&self.ni).arg(&self.ne);
            l.arg(&d.owner)
                .arg(&d.neighbour)
                .arg(&d.cell_ptr)
                .arg(&d.cell_faces);
            l.arg(&d.lsq.x).arg(&d.lsq.y).arg(&d.lsq.z);
            l.arg(&self.vecs[V::Phi.index()]).arg(&d.bc_a).arg(&d.bc_b);
            l.arg(&mut gx).arg(&mut gy).arg(&mut gz);
            // SAFETY: argument list of fv_gradient (gpu/ptx-kernels/src/fv.rs).
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::Gx, gx);
        self.put(V::Gy, gy);
        self.put(V::Gz, gz);
        Ok(())
    }

    fn limiter(&mut self, lim: Limiter) -> anyhow::Result<()> {
        let n_eps2: u32 = match lim {
            Limiter::None => {
                let ones = vec![1.0; self.nc as usize];
                return self.upload(V::Psi, &ones);
            }
            Limiter::Venkatakrishnan { k } => {
                if self.venkat_k != Some(k) {
                    let e = self.host.venkat_eps2(k);
                    self.stream.memcpy_htod(&e, &mut self.d.eps2_cell)?;
                    self.venkat_k = Some(k);
                }
                self.nc
            }
            Limiter::VenkatakrishnanWang { k } => {
                let (mn, mx) = self.min_max(V::Phi)?;
                let e = k * (mx - mn);
                self.stream.memcpy_htod(&[e * e], &mut self.d.eps2_one)?;
                1
            }
            // a negative single ε² selects Barth–Jespersen
            Limiter::BarthJespersen => {
                self.stream.memcpy_htod(&[-1.0], &mut self.d.eps2_one)?;
                1
            }
        };
        let mut psi = self.take(V::Psi)?;
        {
            let d = &self.d;
            let eps2 = if n_eps2 == 1 {
                &d.eps2_one
            } else {
                &d.eps2_cell
            };
            let f = self.func("fv_limiter");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc)
                .arg(&self.nf)
                .arg(&self.ni)
                .arg(&self.ne)
                .arg(&n_eps2);
            l.arg(&d.owner)
                .arg(&d.neighbour)
                .arg(&d.cell_ptr)
                .arg(&d.cell_faces);
            l.arg(&self.vecs[V::Phi.index()]).arg(&d.bc_a).arg(&d.bc_b);
            l.arg(&self.vecs[V::Gx.index()])
                .arg(&self.vecs[V::Gy.index()])
                .arg(&self.vecs[V::Gz.index()]);
            l.arg(&d.cc.x).arg(&d.cc.y).arg(&d.cc.z);
            l.arg(&d.fc.x).arg(&d.fc.y).arg(&d.fc.z);
            l.arg(eps2).arg(&mut psi);
            // SAFETY: argument list of fv_limiter.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::Psi, psi);
        Ok(())
    }

    fn faces(&mut self, st: FaceSettings) -> anyhow::Result<()> {
        let mut cp = self.take(V::Cp)?;
        let mut cn = self.take(V::Cn)?;
        let mut bf = self.take(V::Bf)?;
        let nonorth = if st.nonorth { 1.0f64 } else { 0.0 };
        {
            let d = &self.d;
            let f = self.func("fv_faces");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&self.nf).arg(&self.ni);
            l.arg(&st.blend).arg(&st.central).arg(&st.beta);
            l.arg(&nonorth).arg(&st.nonorth_limit);
            l.arg(&d.owner).arg(&d.neighbour).arg(&d.bc_a).arg(&d.bc_b);
            l.arg(&d.mass_flux).arg(&d.gamma).arg(&d.dgeo);
            l.arg(&d.k.x).arg(&d.k.y).arg(&d.k.z).arg(&d.weight);
            l.arg(&d.fc.x).arg(&d.fc.y).arg(&d.fc.z);
            l.arg(&d.cc.x).arg(&d.cc.y).arg(&d.cc.z);
            l.arg(&self.vecs[V::Phi.index()]);
            l.arg(&self.vecs[V::Gx.index()])
                .arg(&self.vecs[V::Gy.index()])
                .arg(&self.vecs[V::Gz.index()]);
            l.arg(&self.vecs[V::Psi.index()]);
            l.arg(&mut cp).arg(&mut cn).arg(&mut bf);
            // SAFETY: argument list of fv_faces.
            unsafe { l.launch(cfg(self.nf as usize)) }?;
        }
        self.put(V::Cp, cp);
        self.put(V::Cn, cn);
        self.put(V::Bf, bf);
        Ok(())
    }

    fn assemble(&mut self, alpha: f64) -> anyhow::Result<()> {
        let mut rhs = self.take(V::Rhs)?;
        let mut dinv = self.take(V::Dinv)?;
        {
            let d = &mut self.d;
            let i = ENTRIES.iter().position(|e| *e == "fv_assemble").unwrap();
            let mut l = self.stream.launch_builder(&self.f[i]);
            l.arg(&self.nc).arg(&self.nf).arg(&self.ni).arg(&self.ne);
            l.arg(&self.chunk).arg(&self.n_slices).arg(&alpha);
            l.arg(&d.owner)
                .arg(&d.neighbour)
                .arg(&d.cell_ptr)
                .arg(&d.cell_faces);
            l.arg(&self.vecs[V::Cp.index()])
                .arg(&self.vecs[V::Cn.index()])
                .arg(&self.vecs[V::Bf.index()]);
            l.arg(&d.source).arg(&self.vecs[V::Phi.index()]);
            l.arg(&d.slice_ptr)
                .arg(&mut d.val)
                .arg(&mut rhs)
                .arg(&mut dinv);
            // SAFETY: argument list of fv_assemble.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::Rhs, rhs);
        self.put(V::Dinv, dinv);
        Ok(())
    }

    fn spmv(&mut self, x: V, out: V) -> anyhow::Result<()> {
        let mut y = self.take(out)?;
        {
            let d = &self.d;
            let f = self.func("g1_spmv_sell");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc)
                .arg(&self.chunk)
                .arg(&self.n_slices)
                .arg(&self.stored);
            l.arg(&d.slice_ptr)
                .arg(&d.slice_len)
                .arg(&d.col)
                .arg(&d.val);
            l.arg(&self.vecs[x.index()]).arg(&mut y);
            // SAFETY: argument list of g1_spmv_sell; σ = 1, so the sorted
            // rows are the cell rows.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(out, y);
        Ok(())
    }

    fn residual(&mut self, x: V, b: V, out: V) -> anyhow::Result<()> {
        anyhow::ensure!(out != x && out != b, "residual: output aliases an input");
        let mut y = self.take(out)?;
        {
            let d = &self.d;
            let f = self.func("fv_residual");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc)
                .arg(&self.chunk)
                .arg(&self.n_slices)
                .arg(&self.stored);
            l.arg(&d.slice_ptr)
                .arg(&d.slice_len)
                .arg(&d.col)
                .arg(&d.val);
            l.arg(&self.vecs[x.index()])
                .arg(&self.vecs[b.index()])
                .arg(&mut y);
            // SAFETY: argument list of fv_residual.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(out, y);
        Ok(())
    }

    fn dot(&mut self, a: V, b: V) -> anyhow::Result<f64> {
        self.block_dot(a, b)
    }

    fn min_max(&mut self, a: V) -> anyhow::Result<(f64, f64)> {
        let nb = (self.nc as usize).div_ceil(BLOCK);
        let nbu = nb as u32;
        {
            let d = &mut self.d;
            let i = ENTRIES
                .iter()
                .position(|e| *e == "blas_block_min_max")
                .unwrap();
            let mut l = self.stream.launch_builder(&self.f[i]);
            l.arg(&self.nc)
                .arg(&nbu)
                .arg(&self.vecs[a.index()])
                .arg(&mut d.sums)
                .arg(&mut d.sums2);
            // SAFETY: argument list of blas_block_min_max.
            unsafe { l.launch(cfg(nb)) }?;
        }
        let mins = self.stream.clone_dtoh(&self.d.sums)?;
        let maxs = self.stream.clone_dtoh(&self.d.sums2)?;
        Ok((
            mins.iter().copied().fold(f64::INFINITY, f64::min),
            maxs.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        ))
    }

    fn update_p(&mut self, beta: f64, omega: f64) -> anyhow::Result<()> {
        let mut p = self.take(V::P)?;
        {
            let f = self.func("blas_update_p");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&beta).arg(&omega);
            l.arg(&self.vecs[V::R.index()])
                .arg(&self.vecs[V::Vv.index()])
                .arg(&mut p);
            // SAFETY: argument list of blas_update_p.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::P, p);
        Ok(())
    }

    fn sub_scaled(&mut self, x: V, z: V, alpha: f64, out: V) -> anyhow::Result<()> {
        anyhow::ensure!(out != x && out != z, "sub_scaled: output aliases an input");
        let mut y = self.take(out)?;
        {
            let f = self.func("blas_sub_scaled");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&alpha);
            l.arg(&self.vecs[x.index()])
                .arg(&self.vecs[z.index()])
                .arg(&mut y);
            // SAFETY: argument list of blas_sub_scaled.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(out, y);
        Ok(())
    }

    fn update_x(&mut self, alpha: f64, omega: f64) -> anyhow::Result<()> {
        let mut x = self.take(V::Phi)?;
        {
            let f = self.func("blas_update_x");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc).arg(&alpha).arg(&omega);
            l.arg(&self.vecs[V::Ph.index()])
                .arg(&self.vecs[V::Sh.index()])
                .arg(&mut x);
            // SAFETY: argument list of blas_update_x.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(V::Phi, x);
        Ok(())
    }

    fn precondition(&mut self, x: V, out: V) -> anyhow::Result<()> {
        anyhow::ensure!(
            out != x && out != V::Dinv,
            "precondition: output aliases an input"
        );
        let mut y = self.take(out)?;
        {
            let f = self.func("blas_scale");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc);
            l.arg(&self.vecs[V::Dinv.index()])
                .arg(&self.vecs[x.index()])
                .arg(&mut y);
            // SAFETY: argument list of blas_scale.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(out, y);
        Ok(())
    }

    fn copy(&mut self, src: V, dst: V) -> anyhow::Result<()> {
        let s = self.stream.clone();
        let (a, b) = self.pair(src, dst);
        s.memcpy_dtod(a, b)?;
        Ok(())
    }

    fn axpy(&mut self, alpha: f64, x: V, y: V) -> anyhow::Result<()> {
        anyhow::ensure!(x != y, "axpy: x aliases y");
        let mut yv = self.take(y)?;
        {
            let f = self.func("blas_axpy");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc)
                .arg(&alpha)
                .arg(&self.vecs[x.index()])
                .arg(&mut yv);
            // SAFETY: argument list of blas_axpy.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(y, yv);
        Ok(())
    }

    fn xpby(&mut self, x: V, beta: f64, y: V) -> anyhow::Result<()> {
        anyhow::ensure!(x != y, "xpby: x aliases y");
        let mut yv = self.take(y)?;
        {
            let f = self.func("blas_xpby");
            let mut l = self.stream.launch_builder(f);
            l.arg(&self.nc)
                .arg(&self.vecs[x.index()])
                .arg(&beta)
                .arg(&mut yv);
            // SAFETY: argument list of blas_xpby.
            unsafe { l.launch(cfg(self.nc as usize)) }?;
        }
        self.put(y, yv);
        Ok(())
    }

    fn matrix_values(&mut self) -> anyhow::Result<Vec<f64>> {
        Ok(self.stream.clone_dtoh(&self.d.val)?)
    }

    fn set_matrix(&mut self, val: &[f64]) -> anyhow::Result<()> {
        anyhow::ensure!(
            val.len() == self.host.sell.stored(),
            "set_matrix: wrong length"
        );
        let s = self.stream.clone();
        s.memcpy_htod(val, &mut self.d.val)?;
        let dinv = foehncraft_fv::backend::cpu::diagonal_inverse(&self.host, val);
        self.upload(V::Dinv, &dinv)
    }

    fn discretisation(&self) -> Option<&Discretisation> {
        Some(&self.host)
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
            h.levels.first().map_or(h.coarse_n, |l| l.a.n_rows) == self.nc as usize,
            "hierarchy size does not match the mesh"
        );
        self.amg = None;
        self.amg = Some(match h.precision {
            Precision::Mixed => AmgGpu::F32(DevAmg::new(&self.stream, &self.module, h)?),
            Precision::Double => AmgGpu::F64(DevAmg::new(&self.stream, &self.module, h)?),
        });
        Ok(())
    }
}

impl LevelOps for CudaBackend {
    fn amg_levels(&self) -> usize {
        match &self.amg {
            Some(AmgGpu::F32(a)) => a.levels.len(),
            Some(AmgGpu::F64(a)) => a.levels.len(),
            None => 0,
        }
    }

    fn amg_steps(&self, l: usize) -> Vec<(f64, f64)> {
        match &self.amg {
            Some(AmgGpu::F32(a)) => a.levels[l].steps.clone(),
            Some(AmgGpu::F64(a)) => a.levels[l].steps.clone(),
            None => Vec::new(),
        }
    }

    fn amg_load_rhs(&mut self, r: V) -> anyhow::Result<()> {
        let s = self.stream.clone();
        let rv = &self.vecs[r.index()];
        match self.amg.as_mut() {
            Some(AmgGpu::F32(a)) => a.load_rhs(&s, rv),
            Some(AmgGpu::F64(a)) => a.load_rhs(&s, rv),
            None => anyhow::bail!("no multigrid hierarchy loaded"),
        }
    }

    fn amg_store(&mut self, z: V) -> anyhow::Result<()> {
        let s = self.stream.clone();
        let zv = &mut self.vecs[z.index()];
        match self.amg.as_mut() {
            Some(AmgGpu::F32(a)) => a.store(&s, zv),
            Some(AmgGpu::F64(a)) => a.store(&s, zv),
            None => anyhow::bail!("no multigrid hierarchy loaded"),
        }
    }

    fn amg_smooth_first(&mut self, l: usize, c2: f64) -> anyhow::Result<()> {
        let s = self.stream.clone();
        with_amg!(self, a => a.smooth_first(&s, l, c2))
    }

    fn amg_smooth_step(&mut self, l: usize, c1: f64, c2: f64) -> anyhow::Result<()> {
        let s = self.stream.clone();
        with_amg!(self, a => a.smooth_step(&s, l, c1, c2))
    }

    fn amg_restrict(&mut self, l: usize) -> anyhow::Result<()> {
        let s = self.stream.clone();
        with_amg!(self, a => a.restrict(&s, l))
    }

    fn amg_coarse(&mut self) -> anyhow::Result<()> {
        let s = self.stream.clone();
        with_amg!(self, a => a.coarse(&s))
    }

    fn amg_prolong(&mut self, l: usize) -> anyhow::Result<()> {
        let s = self.stream.clone();
        with_amg!(self, a => a.prolong(&s, l))
    }
}
