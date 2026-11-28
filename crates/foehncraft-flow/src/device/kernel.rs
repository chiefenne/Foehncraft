//! The flow device written once over an executor: [`KernelFlow`] turns the
//! operations of [`FlowDevice`] into launches of the entries of
//! `foehncraft_kernels::entry` (argument blocks of buffer addresses and
//! scalars), and an [`Exec`] runs them: the CPU executor
//! ([`super::cpu::CpuExec`]) in multi-threaded loops over host memory, the
//! GPU executor (`foehncraft-cuda`) as PTX entries on device memory. Both
//! run the same body for an index, and the reductions use the same blocks
//! and the same fixed tree (D31), so CPU and GPU give bitwise the same
//! results. All field data stay where the executor keeps them; only scalars
//! (norms, extrema, single values) and the explicit downloads of the
//! algorithms (output, restart) reach the host.

use super::{
    BoilingBufs, BoilingOp, BoundaryArgs, Buf, Buf3, CurvatureOps, Eps2, FlowDevice, FluxArgs,
    KepsCellOps, Lsq, Mask, MlpOps, MomentumFaceArgs, MomentumRowArgs, ScalarFaceArgs,
    ScalarRowArgs, ScalarTurbOps, ScalarWallOps, Setup, SstCellOps, TurbFaceOps, TwoFluidBufs,
    TwoFluidOp, VertexStencil, VofBufs, VofOp, WallCellOps, WallFaceOps, WallStatic, ZetaFCellOps,
};
use anyhow::{Result, ensure};
use foehncraft_core::repro::{BLOCK, combine_block_sums};
use foehncraft_fv::amg::{AmgSettings, ReductionSettings};
use foehncraft_fv::solver::Preconditioner;
use foehncraft_kernels::entry::{self as e, Entry, Mesh, Raw, Raw3, RawGrad};
use foehncraft_kernels::flow as k;

/// Where the entries run: storage of the buffers, launches, the linear
/// solvers.
pub trait Exec {
    /// Backend name (`cpu`, `cuda`).
    fn name(&self) -> String;
    /// Device description.
    fn device(&self) -> String;
    /// A new buffer of `n` zeros.
    fn alloc(&mut self, n: usize) -> Result<Buf>;
    /// A new buffer of `n` zeros that is neither a cell nor a face field
    /// (on ranks every entry is the rank's own).
    fn alloc_plain(&mut self, n: usize) -> Result<Buf> {
        self.alloc(n)
    }
    /// On ranks (MPI step W5): fills the halo entries of the cell buffers
    /// among `bufs` from their owners; other buffers are left alone.
    fn exchange(&mut self, _bufs: &[Buf]) -> Result<()> {
        Ok(())
    }
    /// The executor runs one rank of several.
    fn on_ranks(&self) -> bool {
        false
    }
    /// The communicator of the ranks.
    fn comm(&self) -> Option<foehncraft_mesh::decompose::SharedComm> {
        None
    }
    /// A static array of the run (mesh data), kept until the executor ends.
    fn alloc_u32(&mut self, v: &[u32]) -> Result<Raw<u32>>;
    /// The address and length of a buffer.
    fn raw(&self, b: Buf) -> Raw<f64>;
    /// Copies host values into a buffer.
    fn upload(&mut self, b: Buf, v: &[f64]) -> Result<()>;
    /// Copies a buffer to the host.
    fn download(&mut self, b: Buf) -> Result<Vec<f64>>;
    /// Copies a buffer to the host with the halo of a cell buffer on ranks.
    fn download_with_halo(&mut self, b: Buf) -> Result<Vec<f64>> {
        self.download(b)
    }
    /// dst = src.
    fn copy(&mut self, src: Buf, dst: Buf) -> Result<()>;
    /// One value of a buffer.
    fn value(&mut self, b: Buf, i: usize) -> Result<f64>;
    /// Exact dot product of two arrays of the executor (`y = None`: Σ x²),
    /// when the executor can compute it exactly (D31 (1)); `None` tells the
    /// caller to use the fixed block sums. The CPU executor is exact; the
    /// GPU executor keeps the block sums until v1.1 (owner decision of
    /// 2026-11-23), so CPU and GPU agree to rounding, not bitwise.
    fn exact_dot(&mut self, _x: Raw<f64>, _y: Option<Raw<f64>>) -> Option<f64> {
        None
    }
    /// Order-independent minimum and maximum of an array, when the executor
    /// computes it on the host; `None`: the block path.
    fn exact_min_max(&mut self, _x: Raw<f64>) -> Option<(f64, f64)> {
        None
    }
    /// Runs the body of `K` for the indices `0..n` (n > 0).
    fn launch<K: Entry>(&mut self, n: usize, a: &K::Args) -> Result<()>;
    /// Solves `val` x = `rhs` (see [`FlowDevice::solve`]).
    #[allow(clippy::too_many_arguments)]
    fn solve(
        &mut self,
        val: Buf,
        rhs: Buf,
        x: Buf,
        pc: Preconditioner,
        cg: bool,
        tol: f64,
        max_iter: usize,
    ) -> Result<(usize, f64)>;
    /// Multigrid setup on the matrix `val` into the active slot (with
    /// `symmetric_part` from ½(A + Aᵀ)).
    fn amg_setup(&mut self, val: Buf, s: &AmgSettings, symmetric_part: bool) -> Result<()>;
    /// Makes multigrid slot `slot` the active one.
    fn amg_slot(&mut self, slot: usize) -> Result<()>;
    /// Setup of the reduction multigrid (D179) on the matrix `val`, which
    /// need not be symmetric, into the active slot.
    fn reduction_setup(&mut self, val: Buf, s: &ReductionSettings) -> Result<()>;
    /// x = `cycles` multigrid cycles of the active slot on `val` x = `rhs`
    /// from x = 0 (see [`FlowDevice::multigrid_cycles`]).
    fn cycles(&mut self, val: Buf, rhs: Buf, x: Buf, cycles: usize) -> Result<()>;
    /// Waits for queued work.
    fn sync(&mut self) -> Result<()>;
}

/// A [`FlowDevice`] on an executor.
pub struct KernelFlow<E: Exec> {
    x: E,
    m: Mesh,
    rho: f64,
    nc: usize,
    nf: usize,
    nb: usize,
    lsq_u: Raw3,
    lsq_p: Raw3,
    /// Gradient masks: velocity, pressure, pressure correction.
    masks: [Buf; 3],
    eps2: Buf,
    eps2_one: Buf,
    scratch: Buf,
    zero_cells: Buf,
    /// Block sums (minima) and maxima of up to max(3 n_cells, n_faces) values.
    blocks: Buf,
    blocks2: Buf,
    /// The parts of the residual scale, 3 n_cells each.
    rs_v: Buf,
    rs_w: Buf,
    /// Convective-outlet buffers c and β (per boundary face), addressed by
    /// the mesh block.
    conv: Buf3,
    beta: Buf,
    /// Prescribed pressures per boundary face (time-dependent values, D53).
    pset: Buf,
    /// The face viscosity (the buffer behind `m.mu_face`).
    mu_buf: Buf,
    /// The density per cell and per face (the buffers behind `m.rho` and
    /// `m.rho_face`; D150).
    rho_buf: Buf,
    rho_face_buf: Buf,
    /// The balanced-force face arrays (jump offset per face, boundary
    /// pressure term per boundary face; D154), allocated by the first
    /// `set_balanced_force`.
    balanced: Option<(Buf, Buf)>,
    /// The variable-density cell arrays (stage 7c): mass source, ρⁿ,
    /// ρⁿ⁻¹, allocated by the first `set_variable_density`.
    variable: Option<[Buf; 3]>,
    /// The wall arrays of the turbulence models: faces, cells, ptr, slots.
    wall: Option<[Raw<u32>; 5]>,
    /// The stencils of the quadratic fits (D236): cells, ptr, idx.
    quad: Vec<[Raw<u32>; 3]>,
    /// The vertex stencils of the MLP limiter.
    vert: Vec<VertSlot>,
    /// On ranks: the shared vertices of each stencil of `vert` (MPI step
    /// W7).
    vert_shared: Vec<Option<VertShared>>,
    /// The open-boundary flags of the volume-of-fluid curvature.
    vof_open: Raw<u32>,
    /// The cell body force of the D-weighted Rhie–Chow flux.
    body: Option<[Buf; 3]>,
    /// The inflow flags of the two-fluid face arrays (MP3).
    two_fluid_inflow: Raw<u32>,
    /// The heated wall slots of the boiling (MP4).
    boil_slots: Option<BoilSlots>,
}

/// The vertices a rank shares with other ranks (on their processor faces)
/// and, per shared vertex, where the other ranks hold it: the bounds of a
/// shared vertex are the minimum and the maximum over the ranks.
#[derive(Debug, Clone)]
struct VertShared {
    /// The local vertex of each shared vertex.
    local: Vec<u32>,
    /// Per shared vertex: (rank, position in that rank's list).
    others: Vec<Vec<(u32, u32)>>,
}

/// A vertex stencil on the device ([`FlowDevice::vertex_setup`]).
#[derive(Debug, Clone, Copy)]
struct VertSlot {
    /// vert_ptr, vert_cells, bface_ptr, bfaces, cell_ptr, cell_verts.
    lists: [Raw<u32>; 6],
    points: Buf3,
    vmin: Buf,
    vmax: Buf,
    n_vert: usize,
}

impl<E: Exec> KernelFlow<E> {
    /// Uploads the static data of `s` to the executor.
    pub fn new(mut x: E, s: &Setup) -> Result<Self> {
        let d = s.disc;
        let (nc, nf, ni) = (d.n_cells, d.n_faces, d.n_interior);
        let nb = nf - ni;
        let mut f64s = |v: &[f64]| -> Result<Raw<f64>> {
            let b = x.alloc(v.len())?;
            x.upload(b, v)?;
            Ok(x.raw(b))
        };
        let mut soa = |v: &foehncraft_fv::disc::Soa3| -> Result<Raw3> {
            Ok(Raw3 {
                x: f64s(&v.x)?,
                y: f64s(&v.y)?,
                z: f64s(&v.z)?,
            })
        };
        let area = soa(&d.face_area)?;
        let face_centroid = soa(&d.face_centroid)?;
        let cell_centroid = soa(&d.cell_centroid)?;
        let kk = soa(&d.k)?;
        let shift = soa(&d.shift)?;
        let lsq_u = soa(&d.lsq)?;
        let lsq_p = soa(s.lsq_p)?;
        let mut f64s = |v: &[f64]| -> Result<Raw<f64>> {
            let b = x.alloc(v.len())?;
            x.upload(b, v)?;
            Ok(x.raw(b))
        };
        let dgeo = f64s(&d.dgeo)?;
        let weight = f64s(&d.weight)?;
        let volume = f64s(&d.cell_volume)?;

        // the viscosity per face: one value, or the property program's
        // values (D94; stage 6 adds the eddy viscosity)
        let mu_buf = x.alloc(nf)?;
        match s.mu_face {
            Some(m) => x.upload(mu_buf, m)?,
            None => x.upload(mu_buf, &vec![s.mu; nf])?,
        }
        let mu_face = x.raw(mu_buf);
        // the density per cell and per face: today's constant density
        // filled in (a variable density uploads its fields with
        // `set_density`)
        let rho_buf = x.alloc(nc)?;
        x.upload(rho_buf, &vec![s.rho; nc])?;
        let rho_face_buf = x.alloc(nf)?;
        x.upload(rho_face_buf, &vec![s.rho; nf])?;
        let (rho_raw, rho_face) = (x.raw(rho_buf), x.raw(rho_face_buf));
        let pset_buf = x.alloc(nb)?;
        x.upload(pset_buf, &s.bc.pset)?;
        let pset = x.raw(pset_buf);
        let closed = x.alloc(nb)?;
        let conv = Buf3 {
            x: x.alloc(nb)?,
            y: x.alloc(nb)?,
            z: x.alloc(nb)?,
        };
        let beta = x.alloc(nb)?;
        let m = Mesh {
            n_cells: nc as u64,
            n_faces: nf as u64,
            n_interior: ni as u64,
            owner: x.alloc_u32(&d.owner)?,
            neighbour: x.alloc_u32(&d.neighbour)?,
            cell_ptr: x.alloc_u32(&d.cell_ptr)?,
            cell_faces: x.alloc_u32(&d.cell_faces)?,
            area,
            face_centroid,
            cell_centroid,
            dgeo,
            k: kk,
            weight,
            volume,
            lsq_u,
            lsq_p,
            chunk: d.sell.chunk as u64,
            slice_ptr: x.alloc_u32(&d.sell.slice_ptr)?,
            transpose: x.alloc_u32(&d.transpose_map())?,
            mom_kind: x.alloc_u32(&s.bc.mom_kind)?,
            flux_kind: x.alloc_u32(&s.bc.flux_kind)?,
            ub_rule: x.alloc_u32(&s.bc.ub_rule)?,
            pb_rule: x.alloc_u32(&s.bc.pb_rule)?,
            pset,
            mu_face,
            rho: rho_raw,
            rho_face,
            closed: x.raw(closed),
            conv: Raw3 {
                x: x.raw(conv.x),
                y: x.raw(conv.y),
                z: x.raw(conv.z),
            },
            beta: x.raw(beta),
            periodic: d.periodic as u64,
            shift,
            jump_offset: Raw::EMPTY,
            pb_offset: Raw::EMPTY,
            mass_source: Raw::EMPTY,
            rho_old: Raw::EMPTY,
            rho_old2: Raw::EMPTY,
            body: Raw3::default(),
        };
        let pca: Vec<f64> =
            s.bc.flux_kind
                .iter()
                .map(|&kind| if kind == k::FLUX_FIXED { 1.0 } else { 0.0 })
                .collect();
        let mut buf = |v: &[f64]| -> Result<Buf> {
            let b = x.alloc(v.len())?;
            x.upload(b, v)?;
            Ok(b)
        };
        let masks = [
            buf(&s.bc.u_zero_gradient)?,
            buf(&s.bc.p_zero_gradient)?,
            buf(&pca)?,
        ];
        let eps2 = buf(s.eps2)?;
        let n_blocks = (3 * nc).max(nf).div_ceil(BLOCK);
        Ok(Self {
            eps2,
            eps2_one: x.alloc(1)?,
            scratch: x.alloc(nb)?,
            zero_cells: x.alloc(nc)?,
            blocks: x.alloc_plain(n_blocks)?,
            blocks2: x.alloc_plain(n_blocks)?,
            rs_v: x.alloc_plain(3 * nc)?,
            rs_w: x.alloc_plain(3 * nc)?,
            masks,
            conv,
            beta,
            pset: pset_buf,
            mu_buf,
            rho_buf,
            rho_face_buf,
            balanced: None,
            variable: None,
            wall: None,
            quad: Vec::new(),
            vert: Vec::new(),
            vert_shared: Vec::new(),
            vof_open: Raw::EMPTY,
            two_fluid_inflow: Raw::EMPTY,
            boil_slots: None,
            body: None,
            lsq_u,
            lsq_p,
            rho: s.rho,
            nc,
            nf,
            nb,
            m,
            x,
        })
    }

    /// The executor.
    pub fn exec(&mut self) -> &mut E {
        &mut self.x
    }

    fn r(&self, b: Buf) -> Raw<f64> {
        self.x.raw(b)
    }

    fn r3(&self, b: Buf3) -> Raw3 {
        Raw3 {
            x: self.r(b.x),
            y: self.r(b.y),
            z: self.r(b.z),
        }
    }

    fn rg(&self, g: &[Buf3; 3]) -> RawGrad {
        RawGrad {
            g: [self.r3(g[0]), self.r3(g[1]), self.r3(g[2])],
        }
    }

    /// [`e::CurvatureArgs`] with the gradient and the output set, the rest
    /// empty.
    fn curvature_args(&self, grad_u: [Buf3; 3], out: Buf) -> e::CurvatureArgs {
        let none = Raw::EMPTY;
        e::CurvatureArgs {
            grad_u: self.rg(&grad_u),
            comp: 0,
            u: Raw3::default(),
            v_frame: Raw3::default(),
            grad: Raw3::default(),
            s: none,
            s_old: none,
            inv_dt: 0.0,
            dsdt: [none; 6],
            w: none,
            omega_rot: Raw3::default(),
            rot: [0.0; 3],
            spin: [0.0; 3],
            c: foehncraft_kernels::turb::sst::SstConstants::default(),
            smooth: 0.0,
            out: self.r(out),
        }
    }

    fn mask(&self, m: Mask) -> Buf {
        self.masks[match m {
            Mask::Velocity => 0,
            Mask::Pressure => 1,
            Mask::PressureCorrection => 2,
        }]
    }

    fn go<K: Entry>(&mut self, n: usize, a: &K::Args) -> Result<()> {
        if n == 0 {
            return Ok(());
        }
        self.x.launch::<K>(n, a)
    }

    /// The ε² array of a limiter.
    fn eps2_raw(&mut self, eps2: Eps2) -> Result<Raw<f64>> {
        Ok(match eps2 {
            Eps2::PerCell => {
                let r = self.r(self.eps2);
                ensure!(
                    r.len as usize >= self.nc,
                    "flow device: the limiter has no per-cell epsilon (it was not set up)"
                );
                r
            }
            Eps2::Single(v) => {
                self.x.upload(self.eps2_one, &[v])?;
                self.r(self.eps2_one)
            }
            Eps2::Cells(b) => self.r(b),
        })
    }

    /// Boundary constants of `value` with mask `m` into the scratch array.
    fn constants(&mut self, m: Mask, value: Buf) -> Result<()> {
        let a = e::ConstArgs {
            mask: self.r(self.mask(m)),
            value: self.r(value),
            out: self.r(self.scratch),
        };
        self.go::<e::BoundaryConstants>(self.nb, &a)
    }

    /// Σ x² of the first `n` entries of the raw array in the fixed blocks.
    fn sum_sq_raw(&mut self, x: Raw<f64>) -> Result<f64> {
        if let Some(v) = self.x.exact_dot(x, None) {
            return Ok(v);
        }
        let n = x.len as usize;
        let nbk = n.div_ceil(BLOCK);
        let mut out = self.r(self.blocks);
        ensure!(
            nbk as u64 <= out.len,
            "sum of squares of {n} values: too long"
        );
        out.len = nbk as u64;
        let a = e::BlockArgs {
            x,
            out,
            out2: Raw::EMPTY,
        };
        self.go::<e::BlockSumSq>(nbk, &a)?;
        let sums = self.x.download(self.blocks)?;
        Ok(combine_block_sums(&sums[..nbk]))
    }

    fn lin<K: Entry<Args = e::LinArgs>>(&mut self, n: usize, a: e::LinArgs) -> Result<()> {
        self.go::<K>(n, &a)
    }

    fn unary<K: Entry<Args = e::UnaryArgs>>(
        &mut self,
        n: usize,
        s: f64,
        x: Buf,
        out: Buf,
    ) -> Result<()> {
        check(&[out], &[x])?;
        let a = e::UnaryArgs {
            m: self.m,
            s,
            x: self.r(x),
            out: self.r(out),
        };
        self.go::<K>(n, &a)?;
        // an own-cell loop: the halo of its output from the owners
        if n == self.nc {
            self.x.exchange(&[out])?;
        }
        Ok(())
    }
}

impl<E: Exec> KernelFlow<E> {
    /// The bounds of the shared vertices of a stencil on ranks (MPI step
    /// W7): the minimum and the maximum of the bounds of every rank that
    /// holds the vertex, which are those of all cells and boundary faces
    /// of the vertex in the whole mesh (the same in any order).
    fn combine_bounds(&mut self, v: &VertSlot, sh: &VertShared) -> Result<()> {
        let Some(c) = self.x.comm() else {
            return Ok(());
        };
        let mut lo = self.x.download(v.vmin)?;
        let mut hi = self.x.download(v.vmax)?;
        let n = sh.local.len();
        let mut mine = Vec::with_capacity(2 * n);
        mine.extend(sh.local.iter().map(|&i| lo[i as usize]));
        mine.extend(sh.local.iter().map(|&i| hi[i as usize]));
        let all = c.all_gather_vec_f64(&mine);
        for (k, &i) in sh.local.iter().enumerate() {
            let (mut a, mut b) = (lo[i as usize], hi[i as usize]);
            for &(r, j) in &sh.others[k] {
                let theirs = &all[r as usize];
                let m = theirs.len() / 2;
                a = a.min(theirs[j as usize]);
                b = b.max(theirs[m + j as usize]);
            }
            lo[i as usize] = a;
            hi[i as usize] = b;
        }
        self.x.upload(v.vmin, &lo)?;
        self.x.upload(v.vmax, &hi)
    }

    /// Fails on ranks: the operation does not run there yet.
    fn not_on_ranks(&self, what: &str) -> Result<()> {
        ensure!(
            !self.x.on_ranks(),
            "{what} does not run on ranks yet (MPI step W5)"
        );
        Ok(())
    }

    fn product<K: Entry<Args = e::ProductArgs>>(
        &mut self,
        val: Buf,
        x: Buf,
        out: Buf,
    ) -> Result<()> {
        check(&[out], &[val, x])?;
        let args = e::ProductArgs {
            m: self.m,
            val: self.r(val),
            x: self.r(x),
            out: self.r(out),
        };
        self.go::<K>(self.nc, &args)?;
        self.x.exchange(&[out])
    }
}

/// Fails when an output is also an input or appears twice.
fn check(outs: &[Buf], inputs: &[Buf]) -> Result<()> {
    for (i, o) in outs.iter().enumerate() {
        ensure!(
            !inputs.contains(o) && !outs[..i].contains(o),
            "flow device: output buffer {o:?} is also an input"
        );
    }
    Ok(())
}

/// The heated wall slots of the boiling (MP4): per slot the near-wall
/// cell, area, wall distance and heat flux; the near-wall cells with the
/// start of their slots.
#[derive(Debug, Clone, Copy)]
struct BoilSlots {
    cell: Raw<u32>,
    area: Buf,
    y: Buf,
    q: Buf,
    wcell: Raw<u32>,
    wptr: Raw<u32>,
    n: usize,
}

impl<E: Exec> FlowDevice for KernelFlow<E> {
    fn name(&self) -> String {
        self.x.name()
    }

    fn device(&self) -> String {
        self.x.device()
    }

    fn alloc(&mut self, n: usize) -> Result<Buf> {
        self.x.alloc(n)
    }

    fn upload(&mut self, b: Buf, v: &[f64]) -> Result<()> {
        // the executor checks the length (on ranks a cell buffer takes the
        // own cells or the own and halo cells)
        self.x.upload(b, v)
    }

    fn download(&mut self, b: Buf) -> Result<Vec<f64>> {
        self.x.download(b)
    }

    fn download_with_halo(&mut self, b: Buf) -> Result<Vec<f64>> {
        self.x.download_with_halo(b)
    }

    fn comm(&self) -> Option<foehncraft_mesh::decompose::SharedComm> {
        self.x.comm()
    }

    fn copy(&mut self, src: Buf, dst: Buf) -> Result<()> {
        if src != dst {
            self.x.copy(src, dst)?;
        }
        Ok(())
    }

    fn fill(&mut self, b: Buf, v: f64) -> Result<()> {
        let r = self.r(b);
        self.lin::<e::Fill>(
            r.len as usize,
            e::LinArgs {
                a: v,
                b: 0.0,
                x: Raw::EMPTY,
                y: Raw::EMPTY,
                out: r,
            },
        )
    }

    fn schur_diagonal(&mut self, out: Buf, a: f64, b: f64, ds: Buf, dc: Buf) -> Result<()> {
        let r = self.r(out);
        let args = e::LinArgs {
            a,
            b,
            x: self.r(ds),
            y: self.r(dc),
            out: r,
        };
        self.lin::<e::SchurDiagonal>(r.len as usize, args)
    }

    fn pseudo_coefficient(
        &mut self,
        out: Buf,
        inv_dt: f64,
        factor: f64,
        num: Option<Buf>,
        den: Option<Buf>,
    ) -> Result<()> {
        let ins: Vec<Buf> = num.into_iter().chain(den).collect();
        check(&[out], &ins)?;
        let opt = |b: Option<Buf>| b.map_or(Raw::EMPTY, |b| self.r(b));
        let args = e::PseudoArgs {
            m: self.m,
            inv_dt,
            factor,
            num: opt(num),
            den: opt(den),
            out: self.r(out),
        };
        self.go::<e::PseudoCoefficient>(self.nc, &args)?;
        self.x.exchange(&[out])
    }

    fn scaled_difference(&mut self, out: Buf, a: f64, x: Buf, y: Buf) -> Result<()> {
        let r = self.r(out);
        let args = e::LinArgs {
            a,
            b: 0.0,
            x: self.r(x),
            y: self.r(y),
            out: r,
        };
        self.lin::<e::ScaledDifference>(r.len as usize, args)
    }

    fn squares_add(&mut self, v: Buf3, out: Buf) -> Result<()> {
        check(&[out], &v.all())?;
        let args = e::VecArgs {
            s: Raw::EMPTY,
            v: self.r3(v),
            out3: Raw3::default(),
            out: self.r(out),
        };
        self.go::<e::SquaresAdd>(self.r(out).len as usize, &args)
    }

    fn scale_vector(&mut self, s: Buf, v: Buf3, out: Buf3) -> Result<()> {
        check(&out.all(), &[s])?;
        let args = e::VecArgs {
            s: self.r(s),
            v: self.r3(v),
            out3: self.r3(out),
            out: Raw::EMPTY,
        };
        self.go::<e::ScaleVector>(self.r(s).len as usize, &args)
    }

    fn product(&mut self, out: Buf, a: f64, x: Buf, y: Buf) -> Result<()> {
        let r = self.r(out);
        let args = e::LinArgs {
            a,
            b: 0.0,
            x: self.r(x),
            y: self.r(y),
            out: r,
        };
        self.lin::<e::Product>(r.len as usize, args)
    }

    fn lincomb(&mut self, out: Buf, a: f64, x: Buf, b: f64, y: Buf) -> Result<()> {
        let r = self.r(out);
        let args = e::LinArgs {
            a,
            b,
            x: self.r(x),
            y: self.r(y),
            out: r,
        };
        self.lin::<e::Lincomb>(r.len as usize, args)
    }

    fn sum_sq(&mut self, x: Buf) -> Result<f64> {
        let r = self.r(x);
        self.sum_sq_raw(r)
    }

    fn dot(&mut self, x: Buf, y: Buf) -> Result<f64> {
        let (rx, ry) = (self.r(x), self.r(y));
        ensure!(
            rx.len == ry.len,
            "dot product of {} and {} values",
            rx.len,
            ry.len
        );
        if let Some(v) = self.x.exact_dot(rx, Some(ry)) {
            return Ok(v);
        }
        let n = rx.len as usize;
        let nbk = n.div_ceil(BLOCK);
        let mut out = self.r(self.blocks);
        ensure!(nbk as u64 <= out.len, "dot product of {n} values: too long");
        out.len = nbk as u64;
        let a = e::DotArgs { x: rx, y: ry, out };
        self.go::<e::BlockDot>(nbk, &a)?;
        let sums = self.x.download(self.blocks)?;
        Ok(combine_block_sums(&sums[..nbk]))
    }

    fn min_max(&mut self, x: Buf) -> Result<(f64, f64)> {
        let r = self.r(x);
        if let Some(v) = self.x.exact_min_max(r) {
            return Ok(v);
        }
        let nbk = (r.len as usize).div_ceil(BLOCK);
        let (mut o1, mut o2) = (self.r(self.blocks), self.r(self.blocks2));
        ensure!(
            nbk as u64 <= o1.len,
            "minimum and maximum of {} values: too long",
            r.len
        );
        o1.len = nbk as u64;
        o2.len = nbk as u64;
        self.go::<e::BlockMinMax>(
            nbk,
            &e::BlockArgs {
                x: r,
                out: o1,
                out2: o2,
            },
        )?;
        let mins = self.x.download(self.blocks)?;
        let maxs = self.x.download(self.blocks2)?;
        Ok((
            mins[..nbk].iter().copied().fold(f64::INFINITY, f64::min),
            maxs[..nbk]
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
        ))
    }

    fn value(&mut self, b: Buf, i: usize) -> Result<f64> {
        self.x.value(b, i)
    }

    fn subtract(&mut self, b: Buf, s: f64) -> Result<()> {
        let r = self.r(b);
        self.lin::<e::Subtract>(
            r.len as usize,
            e::LinArgs {
                a: s,
                b: 0.0,
                x: Raw::EMPTY,
                y: Raw::EMPTY,
                out: r,
            },
        )
    }

    fn scale_entry(&mut self, b: Buf, i: usize, s: f64) -> Result<()> {
        let r = self.r(b);
        ensure!(i < r.len as usize, "scale_entry: index {i} out of range");
        self.lin::<e::ScaleEntry>(
            1,
            e::LinArgs {
                a: s,
                b: i as f64,
                x: Raw::EMPTY,
                y: Raw::EMPTY,
                out: r,
            },
        )
    }

    fn boundary_values(&mut self, a: &BoundaryArgs) -> Result<()> {
        let outs = [a.ub.x, a.ub.y, a.ub.z, a.pb];
        let mut ins = vec![a.p, a.mass_flux];
        ins.extend(a.u.all());
        ins.extend(a.grad_p.all());
        ins.extend(a.ub_set.all());
        for g in &a.grad_u {
            ins.extend(g.all());
        }
        check(&outs, &ins)?;
        let args = e::BoundaryValuesArgs {
            m: self.m,
            mass_flux: self.r(a.mass_flux),
            u: self.r3(a.u),
            grad_u: self.rg(&a.grad_u),
            p: self.r(a.p),
            grad_p: self.r3(a.grad_p),
            ub_set: self.r3(a.ub_set),
            ub: self.r3(a.ub),
            pb: self.r(a.pb),
        };
        self.go::<e::BoundaryValues>(self.nb, &args)
    }

    fn gradient(&mut self, w: Lsq, phi: Buf, mask: Mask, value: Buf, out: Buf3) -> Result<()> {
        check(&out.all(), &[phi, value])?;
        self.constants(mask, value)?;
        let args = e::GradientArgs {
            m: self.m,
            w: match w {
                Lsq::Velocity => self.lsq_u,
                Lsq::Pressure => self.lsq_p,
            },
            phi: self.r(phi),
            bc_a: self.r(self.mask(mask)),
            bc_b: self.r(self.scratch),
            jump: Raw::EMPTY,
            out: self.r3(out),
        };
        self.go::<e::Gradient>(self.nc, &args)?;
        self.x.exchange(&out.all())
    }

    fn limiter(&mut self, phi: Buf, value: Buf, grad: Buf3, eps2: Eps2, out: Buf) -> Result<()> {
        let mut ins = vec![phi, value];
        ins.extend(grad.all());
        check(&[out], &ins)?;
        self.constants(Mask::Velocity, value)?;
        let e2 = self.eps2_raw(eps2)?;
        let args = e::LimiterArgs {
            m: self.m,
            phi: self.r(phi),
            bc_a: self.r(self.mask(Mask::Velocity)),
            bc_b: self.r(self.scratch),
            jump: Raw::EMPTY,
            grad: self.r3(grad),
            eps2: e2,
            out: self.r(out),
        };
        self.go::<e::Limiter>(self.nc, &args)?;
        self.x.exchange(&[out])
    }

    fn gradient_with(
        &mut self,
        w: Buf3,
        phi: Buf,
        bc_a: Buf,
        bc_b: Buf,
        jump: Option<Buf>,
        out: Buf3,
    ) -> Result<()> {
        let mut ins = vec![phi, bc_a, bc_b];
        ins.extend(jump);
        ins.extend(w.all());
        check(&out.all(), &ins)?;
        let args = e::GradientArgs {
            m: self.m,
            w: self.r3(w),
            phi: self.r(phi),
            bc_a: self.r(bc_a),
            bc_b: self.r(bc_b),
            jump: jump.map_or(Raw::EMPTY, |b| self.r(b)),
            out: self.r3(out),
        };
        self.go::<e::Gradient>(self.nc, &args)?;
        self.x.exchange(&out.all())
    }

    fn limiter_with(
        &mut self,
        phi: Buf,
        bc_a: Buf,
        bc_b: Buf,
        jump: Option<Buf>,
        grad: Buf3,
        eps2: Eps2,
        out: Buf,
    ) -> Result<()> {
        let mut ins = vec![phi, bc_a, bc_b];
        ins.extend(jump);
        ins.extend(grad.all());
        check(&[out], &ins)?;
        let e2 = self.eps2_raw(eps2)?;
        let args = e::LimiterArgs {
            m: self.m,
            phi: self.r(phi),
            bc_a: self.r(bc_a),
            bc_b: self.r(bc_b),
            jump: jump.map_or(Raw::EMPTY, |b| self.r(b)),
            grad: self.r3(grad),
            eps2: e2,
            out: self.r(out),
        };
        self.go::<e::Limiter>(self.nc, &args)?;
        self.x.exchange(&[out])
    }

    fn scalar_faces(&mut self, a: &ScalarFaceArgs) -> Result<()> {
        let mut ins = vec![a.bc_a, a.bc_b, a.mass_flux, a.gamma, a.phi, a.psi];
        ins.extend(a.grad.all());
        check(&[a.cp, a.cn, a.bf], &ins)?;
        let args = e::ScalarFaceArgs {
            m: self.m,
            bc_a: self.r(a.bc_a),
            bc_b: self.r(a.bc_b),
            jump: a.jump.map_or(Raw::EMPTY, |b| self.r(b)),
            mass_flux: self.r(a.mass_flux),
            gamma: self.r(a.gamma),
            phi: self.r(a.phi),
            grad: self.r3(a.grad),
            psi: self.r(a.psi),
            blend: a.blend,
            central: a.central,
            beta: a.beta,
            nonorth: a.nonorth,
            nonorth_limit: -1.0,
            bound: a.bound,
            cp: self.r(a.cp),
            cn: self.r(a.cn),
            bf: self.r(a.bf),
        };
        self.go::<e::ScalarFaces>(self.nf, &args)
    }

    fn scalar_rows(&mut self, a: &ScalarRowArgs) -> Result<()> {
        let mut ins = vec![a.cp, a.cn, a.bf, a.mass_flux, a.phi];
        ins.extend(a.su);
        ins.extend(a.sp);
        if let Some((m, v)) = a.fix {
            ins.push(m);
            ins.push(v);
        }
        if let Some(t) = &a.time {
            ins.push(t.old);
            ins.push(t.old2);
        }
        ins.extend(a.pseudo);
        ins.extend(a.capacity);
        check(&[a.val, a.rhs, a.diag, a.order], &ins)?;
        let opt = |b: Option<Buf>| b.map_or(Raw::EMPTY, |b| self.r(b));
        let (transient, inv_dt, coef, old, old2, (lo, hi)) = match &a.time {
            Some(t) => (
                1,
                t.inv_dt,
                t.a,
                self.r(t.old),
                if t.a[2] != 0.0 {
                    self.r(t.old2)
                } else {
                    Raw::EMPTY
                },
                t.bounds,
            ),
            None => (
                0,
                0.0,
                [1.0, 0.0, 0.0],
                Raw::EMPTY,
                Raw::EMPTY,
                (f64::NEG_INFINITY, f64::INFINITY),
            ),
        };
        let args = e::ScalarRowArgs {
            m: self.m,
            cp: self.r(a.cp),
            cn: self.r(a.cn),
            bf: self.r(a.bf),
            mass_flux: self.r(a.mass_flux),
            phi: self.r(a.phi),
            su: opt(a.su),
            sp: opt(a.sp),
            alpha: a.alpha,
            transient,
            inv_dt,
            a: coef,
            old,
            old2,
            lo,
            hi,
            fix_mask: a.fix.map_or(Raw::EMPTY, |(m, _)| self.r(m)),
            fix_val: a.fix.map_or(Raw::EMPTY, |(_, v)| self.r(v)),
            pseudo: opt(a.pseudo),
            capacity: opt(a.capacity),
            positive: u32::from(a.positive),
            val: self.r(a.val),
            rhs: self.r(a.rhs),
            diag: self.r(a.diag),
            order: self.r(a.order),
        };
        self.go::<e::ScalarRows>(self.nc, &args)?;
        self.x.exchange(&[a.rhs, a.diag, a.order])
    }

    fn scalar_scale(
        &mut self,
        phi: Buf,
        rhs: Buf,
        diag: Buf,
        alpha: f64,
        pseudo: Option<Buf>,
    ) -> Result<f64> {
        let nc = self.nc;
        let mut v = self.r(self.rs_v);
        let mut w = self.r(self.rs_w);
        v.len = nc as u64;
        w.len = nc as u64;
        let args = e::ScalarScaleArgs {
            phi: self.r(phi),
            rhs: self.r(rhs),
            diag: self.r(diag),
            alpha,
            pseudo: pseudo.map_or(Raw::EMPTY, |b| self.r(b)),
            v,
            w,
        };
        self.go::<e::ScalarScale>(nc, &args)?;
        let a = self.sum_sq_raw(v)?;
        let b = self.sum_sq_raw(w)?;
        Ok((a.sqrt() + b.sqrt()).max(f64::MIN_POSITIVE))
    }

    fn mass_flux(&mut self, a: &FluxArgs, out: Buf) -> Result<()> {
        let mut ins = vec![a.p, a.pb, a.relax];
        ins.extend(a.u.all());
        ins.extend(a.grad_p.all());
        ins.extend(a.ub.all());
        for g in &a.grad_u {
            ins.extend(g.all());
        }
        ins.extend(a.d);
        if let Some(psi) = a.bounded {
            ins.extend(psi.all());
        }
        if let Some((ds, y, _)) = a.time {
            ins.push(ds);
            ins.push(y);
        }
        check(&[out], &ins)?;
        let (transient, ds, time, alpha) = match a.time {
            Some((ds, y, al)) => (1, self.r(ds), self.r(y), al),
            None => (0, Raw::EMPTY, Raw::EMPTY, 1.0),
        };
        let args = e::MassFluxArgs {
            m: self.m,
            u: self.r3(a.u),
            grad_u: self.rg(&a.grad_u),
            p: self.r(a.p),
            grad_p: self.r3(a.grad_p),
            pb: self.r(a.pb),
            ub: self.r3(a.ub),
            d: self.r(a.d.unwrap_or(self.zero_cells)),
            relax: self.r(a.relax),
            transient,
            ds,
            time,
            alpha,
            psi: a.bounded.map_or(
                Raw3 {
                    x: Raw::EMPTY,
                    y: Raw::EMPTY,
                    z: Raw::EMPTY,
                },
                |b| self.r3(b),
            ),
            bounded: u32::from(a.bounded.is_some()),
            out: self.r(out),
        };
        self.go::<e::MassFlux>(self.nf, &args)
    }

    fn momentum_faces(&mut self, a: &MomentumFaceArgs) -> Result<()> {
        let outs = [a.cp, a.cn, a.bf.x, a.bf.y, a.bf.z];
        let mut ins = vec![a.mass_flux];
        ins.extend(a.u.all());
        ins.extend(a.uc.all());
        ins.extend(a.psi.all());
        ins.extend(a.ub.all());
        for g in &a.grad_u {
            ins.extend(g.all());
        }
        check(&outs, &ins)?;
        let args = e::MomentumFaceArgs {
            m: self.m,
            mass_flux: self.r(a.mass_flux),
            u: self.r3(a.u),
            uc: self.r3(a.uc),
            grad_u: self.rg(&a.grad_u),
            psi: self.r3(a.psi),
            ub: self.r3(a.ub),
            blend: a.blend,
            central: a.central,
            beta: a.beta,
            nonorth: a.nonorth,
            convection: a.convection,
            transposed: a.transposed,
            trace_free: a.trace_free,
            cp: self.r(a.cp),
            cn: self.r(a.cn),
            bf: self.r3(a.bf),
        };
        self.go::<e::MomentumFaces>(self.nf, &args)
    }

    fn momentum_rows(&mut self, a: &MomentumRowArgs) -> Result<()> {
        let outs = [
            a.val, a.rhs.x, a.rhs.y, a.rhs.z, a.d, a.ds, a.dc, a.diag.x, a.diag.y, a.diag.z,
        ];
        let mut ins = vec![a.cp, a.cn];
        ins.extend(a.bf.all());
        ins.extend(a.grad_p.all());
        ins.extend(a.u.all());
        if let Some(t) = &a.time {
            ins.extend(t.u_old.all());
            ins.extend(t.u_old2.all());
            ins.extend(t.explicit.1.all());
        }
        if let Some(s) = a.sink {
            ins.push(s);
        }
        check(&outs, &ins)?;
        let empty = Raw3::default();
        let (transient, time_a, dt, u_old, u_old2, ew, ex) = match &a.time {
            Some(t) => (
                1,
                t.a,
                t.dt,
                self.r3(t.u_old),
                self.r3(t.u_old2),
                t.explicit.0,
                self.r3(t.explicit.1),
            ),
            None => (0, [0.0; 3], 1.0, empty, empty, 0.0, empty),
        };
        let args = e::MomentumRowArgs {
            m: self.m,
            cp: self.r(a.cp),
            cn: self.r(a.cn),
            bf: self.r3(a.bf),
            grad_p: self.r3(a.grad_p),
            body_force: a.body_force,
            u: self.r3(a.u),
            alpha: a.alpha,
            theta: a.theta,
            transient,
            time_a,
            dt,
            u_old,
            u_old2,
            explicit_weight: ew,
            explicit: ex,
            sink: a.sink.map(|s| self.r(s)).unwrap_or_default(),
            val: self.r(a.val),
            rhs: self.r3(a.rhs),
            d: self.r(a.d),
            ds: self.r(a.ds),
            dc: self.r(a.dc),
            diag: self.r3(a.diag),
        };
        self.go::<e::MomentumRows>(self.nc, &args)?;
        self.x.exchange(&[
            a.rhs.x, a.rhs.y, a.rhs.z, a.d, a.ds, a.dc, a.diag.x, a.diag.y, a.diag.z,
        ])
    }

    fn momentum_operator(&mut self, cp: Buf, cn: Buf, bf: Buf3, u: Buf3, out: Buf3) -> Result<()> {
        let mut ins = vec![cp, cn];
        ins.extend(bf.all());
        ins.extend(u.all());
        check(&out.all(), &ins)?;
        let args = e::OperatorArgs {
            m: self.m,
            cp: self.r(cp),
            cn: self.r(cn),
            bf: self.r3(bf),
            u: self.r3(u),
            out: self.r3(out),
        };
        self.go::<e::MomentumOperator>(self.nc, &args)?;
        self.x.exchange(&out.all())
    }

    fn residual_scale(&mut self, u: Buf3, rhs: Buf3, diag: Buf3, alpha: f64) -> Result<f64> {
        let args = e::ResidualScaleArgs {
            u: self.r3(u),
            rhs: self.r3(rhs),
            diag: self.r3(diag),
            alpha,
            v: self.r(self.rs_v),
            w: self.r(self.rs_w),
        };
        self.go::<e::ResidualScale>(self.nc, &args)?;
        let (v, w) = (self.r(self.rs_v), self.r(self.rs_w));
        let a = self.sum_sq_raw(v)?;
        let b = self.sum_sq_raw(w)?;
        Ok((a.sqrt() + b.sqrt()).max(f64::MIN_POSITIVE))
    }

    fn set_diagonal(&mut self, val: Buf, diag: Buf) -> Result<()> {
        self.unary::<e::SetDiagonal>(self.nc, 0.0, diag, val)
    }

    fn net_outflow(&mut self, m: Buf, out: Buf) -> Result<()> {
        self.unary::<e::NetOutflow>(self.nc, 0.0, m, out)
    }

    fn throughflow(&mut self, m: Buf, out: Buf) -> Result<()> {
        self.unary::<e::Throughflow>(self.nc, 0.0, m, out)
    }

    fn courant(&mut self, m: Buf, dt: f64, out: Buf) -> Result<()> {
        let factor = dt / (2.0 * self.rho);
        self.unary::<e::Courant>(self.nc, factor, m, out)
    }

    fn convective_buffers(&self) -> (Buf3, Buf) {
        (self.conv, self.beta)
    }

    fn set_pressure_values(&mut self, v: &[f64]) -> Result<()> {
        let b = self.pset;
        self.x.upload(b, v)
    }

    fn close_backflow(&mut self, m: Buf) -> Result<()> {
        let args = e::CloseArgs {
            m: self.m,
            mass_flux: self.r(m),
        };
        self.go::<e::CloseBackflow>(self.nb, &args)
    }

    fn pressure_faces(&mut self, dc: Buf, out: Buf) -> Result<()> {
        let rho = self.rho;
        self.unary::<e::PressureFaces>(self.nf, rho, dc, out)
    }

    fn pressure_rows(&mut self, cf: Buf, m: Buf, val: Buf, rhs: Buf) -> Result<()> {
        check(&[val, rhs], &[cf, m])?;
        let args = e::PressureRowArgs {
            m: self.m,
            cf: self.r(cf),
            mass_flux: self.r(m),
            val: self.r(val),
            rhs: self.r(rhs),
        };
        self.go::<e::PressureRows>(self.nc, &args)?;
        self.x.exchange(&[rhs])
    }

    fn quotient(&mut self, out: Buf, a: f64, x: Buf, y: Buf) -> Result<()> {
        let r = self.r(out);
        let args = e::LinArgs {
            a,
            b: 0.0,
            x: self.r(x),
            y: self.r(y),
            out: r,
        };
        self.lin::<e::Quotient>(r.len as usize, args)
    }

    fn set_body_force(&mut self, b: &[Vec<f64>; 3]) -> Result<()> {
        self.not_on_ranks("the volume-of-fluid method")?;
        ensure!(
            b.iter().all(|v| v.len() == self.nc),
            "flow device: a body force of {} cells expected",
            self.nc
        );
        let bufs = match self.body {
            Some(x) => x,
            None => {
                let x = [
                    self.x.alloc(self.nc)?,
                    self.x.alloc(self.nc)?,
                    self.x.alloc(self.nc)?,
                ];
                self.body = Some(x);
                x
            }
        };
        for (bb, v) in bufs.iter().zip(b) {
            self.x.upload(*bb, v)?;
        }
        self.m.body = Raw3 {
            x: self.x.raw(bufs[0]),
            y: self.x.raw(bufs[1]),
            z: self.x.raw(bufs[2]),
        };
        Ok(())
    }

    fn vof_open(&mut self, open: &[u32]) -> Result<()> {
        self.not_on_ranks("the volume-of-fluid method")?;
        ensure!(
            open.len() == self.nb,
            "flow device: {} open-boundary flags expected",
            self.nb
        );
        self.vof_open = if open.is_empty() {
            Raw::EMPTY
        } else {
            self.x.alloc_u32(open)?
        };
        Ok(())
    }

    fn vof(&mut self, op: VofOp, b: &VofBufs) -> Result<()> {
        self.not_on_ranks("the volume-of-fluid method")?;
        let ins: Vec<Buf> = [b.flux, b.alpha, b.alpha_f, b.alpha_b, b.kappa, b.gnorm]
            .into_iter()
            .flatten()
            .chain(b.grad.map(|g| g.all()).into_iter().flatten())
            .collect();
        let outs: Vec<Buf> = [b.out, b.out2]
            .into_iter()
            .flatten()
            .chain(b.out3.map(|g| g.all()).into_iter().flatten())
            .collect();
        check(&outs, &ins)?;
        let o = |x: Option<Buf>| x.map_or(Raw::EMPTY, |x| self.r(x));
        let args = e::VofArgs {
            m: self.m,
            flux: o(b.flux),
            alpha: o(b.alpha),
            alpha_f: o(b.alpha_f),
            grad: b.grad.map_or(Raw3::default(), |g| self.r3(g)),
            alpha_b: o(b.alpha_b),
            open: self.vof_open,
            kappa: o(b.kappa),
            gnorm: o(b.gnorm),
            scheme: b.scheme as u64,
            dt: b.dt,
            beta: b.beta,
            s: b.s,
            out: o(b.out),
            out2: o(b.out2),
            out3: b.out3.map_or(Raw3::default(), |g| self.r3(g)),
        };
        let need = |x: bool, what: &str| -> Result<()> {
            ensure!(
                x,
                "flow device: the volume-of-fluid operation {op:?} needs {what}"
            );
            Ok(())
        };
        match op {
            VofOp::Faces => {
                need(
                    b.flux.is_some() && b.alpha.is_some() && b.grad.is_some() && b.out.is_some(),
                    "flux, alpha, grad, out",
                )?;
                need(b.alpha_b.is_some() || self.nb == 0, "alpha_b")?;
                self.go::<e::VofFaces>(self.nf, &args)
            }
            VofOp::Update => {
                need(
                    b.flux.is_some()
                        && b.alpha_f.is_some()
                        && b.alpha.is_some()
                        && b.out.is_some()
                        && b.out2.is_some(),
                    "flux, alpha_f, alpha, out, out2",
                )?;
                self.go::<e::VofUpdate>(self.nc, &args)
            }
            VofOp::OutflowRate => {
                need(b.flux.is_some() && b.out.is_some(), "flux, out")?;
                self.go::<e::VofOutflowRate>(self.nc, &args)
            }
            VofOp::Smooth => {
                need(b.alpha.is_some() && b.out.is_some(), "alpha, out")?;
                self.go::<e::VofSmooth>(self.nc, &args)
            }
            VofOp::Gradient => {
                need(
                    b.alpha.is_some() && b.out.is_some() && b.out3.is_some(),
                    "alpha, out, out3",
                )?;
                self.go::<e::VofGradient>(self.nc, &args)
            }
            VofOp::Curvature => {
                need(b.grad.is_some() && b.out.is_some(), "grad, out")?;
                need(self.vof_open.len as usize == self.nb, "vof_open first")?;
                self.go::<e::VofCurvature>(self.nc, &args)
            }
            VofOp::Reconstruct => {
                need(
                    b.alpha.is_some() && b.flux.is_some() && b.kappa.is_some() && b.out3.is_some(),
                    "alpha, flux (rho_f), kappa (rho), out3",
                )?;
                need(self.vof_open.len as usize == self.nb, "vof_open first")?;
                self.go::<e::VofReconstruct>(self.nc, &args)?;
                self.x.exchange(&b.out3.expect("out3").all())
            }
            VofOp::Tension => {
                need(
                    b.alpha.is_some() && b.kappa.is_some() && b.gnorm.is_some() && b.out.is_some(),
                    "alpha, kappa, gnorm, out",
                )?;
                self.go::<e::VofTension>(self.nf, &args)
            }
        }
    }

    fn two_fluid_inflow(&mut self, inflow: &[u32]) -> Result<()> {
        self.not_on_ranks("the two-fluid model")?;
        ensure!(
            inflow.len() == self.nb,
            "flow device: {} two-fluid inflow flags expected",
            self.nb
        );
        self.two_fluid_inflow = if inflow.is_empty() {
            Raw::EMPTY
        } else {
            self.x.alloc_u32(inflow)?
        };
        Ok(())
    }

    fn two_fluid(
        &mut self,
        op: TwoFluidOp,
        b: &TwoFluidBufs,
        kc: &foehncraft_kernels::twofluid::Constants<f64>,
    ) -> Result<()> {
        self.not_on_ranks("the two-fluid model")?;
        let one = [
            b.alpha, b.wall_y, b.nu_t, b.alpha_b, b.alpha_gf, b.alpha_lf, b.mu_eff, b.cell_h,
            b.diam,
        ];
        let three = [b.ul, b.ug, b.wall_n, b.grad_alpha];
        let mut ins: Vec<Buf> = one
            .into_iter()
            .flatten()
            .chain(three.into_iter().flatten().flat_map(|g| g.all()))
            .chain(b.grad_ul.into_iter().flatten().flat_map(|g| g.all()))
            .collect();
        let mut outs: Vec<Buf> = Vec::new();
        match op {
            TwoFluidOp::Cells => {
                outs.extend([b.kv, b.mu_bit, b.disp].into_iter().flatten());
                outs.extend(
                    [b.fv, b.body_g, b.body_l]
                        .into_iter()
                        .flatten()
                        .flat_map(|g| g.all()),
                );
            }
            TwoFluidOp::Faces => {
                ins.extend([b.mu_bit, b.disp].into_iter().flatten());
                outs.extend(
                    [b.mu_g, b.mu_l, b.off_g, b.off_l, b.td]
                        .into_iter()
                        .flatten(),
                );
            }
        }
        check(&outs, &ins)?;
        let o = |x: Option<Buf>| x.map_or(Raw::EMPTY, |x| self.r(x));
        let o3 = |x: Option<Buf3>| x.map_or(Raw3::default(), |x| self.r3(x));
        let args = e::TwoFluidArgs {
            m: self.m,
            k: *kc,
            alpha: o(b.alpha),
            ul: o3(b.ul),
            ug: o3(b.ug),
            grad_ul: b.grad_ul.map_or(RawGrad::default(), |g| RawGrad {
                g: [self.r3(g[0]), self.r3(g[1]), self.r3(g[2])],
            }),
            wall_y: o(b.wall_y),
            wall_n: o3(b.wall_n),
            nu_t: o(b.nu_t),
            grad_alpha: o3(b.grad_alpha),
            cell_h: o(b.cell_h),
            alpha_b: o(b.alpha_b),
            inflow: self.two_fluid_inflow,
            alpha_gf: o(b.alpha_gf),
            alpha_lf: o(b.alpha_lf),
            mu_eff: o(b.mu_eff),
            mu_bit: o(b.mu_bit),
            disp: o(b.disp),
            kv: o(b.kv),
            fv: o3(b.fv),
            body_g: o3(b.body_g),
            body_l: o3(b.body_l),
            mu_g: o(b.mu_g),
            mu_l: o(b.mu_l),
            off_g: o(b.off_g),
            off_l: o(b.off_l),
            td: o(b.td),
            diam: o(b.diam),
        };
        let need = |x: bool, what: &str| -> Result<()> {
            ensure!(
                x,
                "flow device: the two-fluid operation {op:?} needs {what}"
            );
            Ok(())
        };
        match op {
            TwoFluidOp::Cells => {
                need(
                    b.alpha.is_some()
                        && b.ul.is_some()
                        && b.ug.is_some()
                        && b.grad_ul.is_some()
                        && b.wall_y.is_some()
                        && b.wall_n.is_some()
                        && b.grad_alpha.is_some()
                        && b.cell_h.is_some(),
                    "alpha, ul, ug, grad_ul, wall_y, wall_n, grad_alpha, cell_h",
                )?;
                need(
                    b.kv.is_some()
                        && b.fv.is_some()
                        && b.mu_bit.is_some()
                        && b.disp.is_some()
                        && b.body_g.is_some()
                        && b.body_l.is_some(),
                    "kv, fv, mu_bit, disp, body_g, body_l",
                )?;
                self.go::<e::TwoFluidCells>(self.nc, &args)
            }
            TwoFluidOp::Faces => {
                need(
                    b.alpha.is_some()
                        && b.alpha_gf.is_some()
                        && b.alpha_lf.is_some()
                        && b.mu_eff.is_some()
                        && b.mu_bit.is_some()
                        && b.disp.is_some(),
                    "alpha, alpha_gf, alpha_lf, mu_eff, mu_bit, disp",
                )?;
                need(
                    self.nb == 0
                        || (b.alpha_b.is_some() && self.two_fluid_inflow.len as usize == self.nb),
                    "alpha_b and two_fluid_inflow first",
                )?;
                need(
                    b.mu_g.is_some()
                        && b.mu_l.is_some()
                        && b.off_g.is_some()
                        && b.off_l.is_some()
                        && b.td.is_some(),
                    "mu_g, mu_l, off_g, off_l, td",
                )?;
                self.go::<e::TwoFluidFaces>(self.nf, &args)
            }
        }
    }

    fn boiling_slots(
        &mut self,
        cell: &[u32],
        area: &[f64],
        y: &[f64],
        q: &[f64],
        wcell: &[u32],
        wptr: &[u32],
    ) -> Result<()> {
        self.not_on_ranks("wall boiling")?;
        let n = cell.len();
        ensure!(
            area.len() == n && y.len() == n && q.len() == n && wptr.len() == wcell.len() + 1,
            "flow device: inconsistent boiling slots"
        );
        ensure!(
            cell.iter().chain(wcell).all(|&c| (c as usize) < self.nc)
                && wptr.windows(2).all(|w| w[0] <= w[1])
                && wptr.last().is_some_and(|&l| l as usize == n),
            "flow device: boiling slots outside the mesh"
        );
        let mut put = |v: &[f64]| -> Result<Buf> {
            let b = self.alloc(v.len())?;
            self.upload(b, v)?;
            Ok(b)
        };
        let (area, y, q) = (put(area)?, put(y)?, put(q)?);
        self.boil_slots = Some(BoilSlots {
            cell: self.x.alloc_u32(cell)?,
            area,
            y,
            q,
            wcell: self.x.alloc_u32(wcell)?,
            wptr: self.x.alloc_u32(wptr)?,
            n,
        });
        Ok(())
    }

    fn boiling(
        &mut self,
        op: BoilingOp,
        b: &BoilingBufs,
        kc: &foehncraft_kernels::boiling::Constants<f64>,
        dt: f64,
    ) -> Result<()> {
        self.not_on_ranks("wall boiling")?;
        let (ins, outs): (Vec<Buf>, Vec<Buf>) = match op {
            BoilingOp::Wall => (
                vec![b.t_l, b.t_g, b.tk, b.alpha],
                vec![b.t_w, b.q_l, b.m_e, b.q_e, b.area, b.q_g],
            ),
            BoilingOp::WallCells => (
                vec![b.q_l, b.m_e, b.q_g],
                vec![b.q_wall, b.m_wall, b.q_wall_g],
            ),
            BoilingOp::Cells => (
                [b.alpha, b.t_l, b.t_g, b.q_wall, b.m_wall, b.q_wall_g]
                    .into_iter()
                    .chain(b.ul.all())
                    .chain(b.ug.all())
                    .collect(),
                vec![
                    b.alpha_out,
                    b.diam,
                    b.gamma_i,
                    b.gamma,
                    b.h_l,
                    b.h_g,
                    b.su_l,
                    b.sp_l,
                    b.su_g,
                    b.sp_g,
                    b.cap_l,
                    b.cap_g,
                    b.src,
                ],
            ),
        };
        check(&outs, &ins)?;
        let sl = self.boil_slots;
        let args = e::BoilingArgs {
            m: self.m,
            k: *kc,
            dt,
            slot_cell: sl.map_or(Raw::EMPTY, |s| s.cell),
            slot_area: sl.map_or(Raw::EMPTY, |s| self.r(s.area)),
            slot_y: sl.map_or(Raw::EMPTY, |s| self.r(s.y)),
            slot_q: sl.map_or(Raw::EMPTY, |s| self.r(s.q)),
            t_l: self.r(b.t_l),
            t_g: self.r(b.t_g),
            tk: self.r(b.tk),
            t_w: self.r(b.t_w),
            q_l: self.r(b.q_l),
            m_e: self.r(b.m_e),
            q_e: self.r(b.q_e),
            area: self.r(b.area),
            q_g: self.r(b.q_g),
            wcell: sl.map_or(Raw::EMPTY, |s| s.wcell),
            wptr: sl.map_or(Raw::EMPTY, |s| s.wptr),
            q_wall: self.r(b.q_wall),
            m_wall: self.r(b.m_wall),
            q_wall_g: self.r(b.q_wall_g),
            alpha: self.r(b.alpha),
            ul: self.r3(b.ul),
            ug: self.r3(b.ug),
            alpha_out: self.r(b.alpha_out),
            diam: self.r(b.diam),
            gamma_i: self.r(b.gamma_i),
            gamma: self.r(b.gamma),
            h_l: self.r(b.h_l),
            h_g: self.r(b.h_g),
            su_l: self.r(b.su_l),
            sp_l: self.r(b.sp_l),
            su_g: self.r(b.su_g),
            sp_g: self.r(b.sp_g),
            cap_l: self.r(b.cap_l),
            cap_g: self.r(b.cap_g),
            src: self.r(b.src),
        };
        match op {
            BoilingOp::Wall | BoilingOp::WallCells => {
                let Some(s) = sl else {
                    return Ok(());
                };
                ensure!(
                    args.t_w.len as usize >= s.n,
                    "flow device: the boiling slot outputs need {} entries",
                    s.n
                );
                if op == BoilingOp::Wall {
                    self.go::<e::BoilingWall>(s.n, &args)
                } else {
                    self.go::<e::BoilingWallCells>(s.wcell.len as usize, &args)
                }
            }
            BoilingOp::Cells => self.go::<e::BoilingCells>(self.nc, &args),
        }
    }

    fn set_mass_source(&mut self, b: Option<Buf>) -> Result<()> {
        self.m.mass_source = match b {
            Some(b) => {
                ensure!(
                    self.r(b).len as usize >= self.nc,
                    "flow device: a mass source of {} cells expected",
                    self.nc
                );
                self.r(b)
            }
            None => Raw::EMPTY,
        };
        Ok(())
    }

    fn two_fluid_eliminate(&mut self, b: &super::EliminateBufs) -> Result<()> {
        self.not_on_ranks("the two-fluid model")?;
        let outs: Vec<Buf> =
            b.ug.all()
                .into_iter()
                .chain(b.ul.all())
                .chain(b.d_g)
                .chain(b.d_l)
                .collect();
        let ins: Vec<Buf> = [b.alpha, b.kv]
            .into_iter()
            .chain(b.diag_g.all())
            .chain(b.diag_l.all())
            .chain(b.ul_lag.all())
            .collect();
        check(&outs, &ins)?;
        let args = e::TwoFluidElimArgs {
            m: self.m,
            alpha: self.r(b.alpha),
            kv: self.r(b.kv),
            diag_g: self.r3(b.diag_g),
            diag_l: self.r3(b.diag_l),
            ug: self.r3(b.ug),
            ul: self.r3(b.ul),
            ul_lag: self.r3(b.ul_lag),
            d_g: b.d_g.map(|x| self.r(x)),
            d_l: b.d_l.map(|x| self.r(x)),
        };
        self.go::<e::TwoFluidEliminate>(self.nc, &args)
    }

    fn material_bufs(&self) -> super::MaterialBufs {
        super::MaterialBufs {
            rho: self.rho_buf,
            rho_face: self.rho_face_buf,
            mu: self.mu_buf,
            body: self.body.map(|b| Buf3 {
                x: b[0],
                y: b[1],
                z: b[2],
            }),
            jump: self.balanced.map(|b| b.0),
        }
    }

    fn correct_fluxes(&mut self, cf: Buf, m: Buf, pc: Buf, out: Buf) -> Result<()> {
        check(&[out], &[cf, m, pc])?;
        let args = e::CorrectFluxArgs {
            m: self.m,
            cf: self.r(cf),
            mass_flux: self.r(m),
            pc: self.r(pc),
            out: self.r(out),
        };
        self.go::<e::CorrectFluxes>(self.nf, &args)
    }

    fn correct_velocity(&mut self, dc: Buf, grad: Buf3, u: Buf3) -> Result<()> {
        let mut ins = vec![dc];
        ins.extend(grad.all());
        check(&u.all(), &ins)?;
        let args = e::CorrectVelocityArgs {
            dc: self.r(dc),
            grad: self.r3(grad),
            u: self.r3(u),
        };
        self.go::<e::CorrectVelocity>(self.nc, &args)?;
        self.x.exchange(&u.all())
    }

    fn neighbour_correction(
        &mut self,
        val: Buf,
        x: Buf,
        diag: Buf,
        out: Buf,
        own: f64,
    ) -> Result<()> {
        check(&[out], &[val, x, diag])?;
        let args = e::NeighbourArgs {
            m: self.m,
            val: self.r(val),
            x: self.r(x),
            diag: self.r(diag),
            out: self.r(out),
            own,
        };
        self.go::<e::NeighbourCorrection>(self.nc, &args)?;
        self.x.exchange(&[out])
    }

    fn delta_flux(&mut self, delta: Buf3, out: Buf) -> Result<()> {
        check(&[out], &delta.all())?;
        let args = e::DeltaFluxArgs {
            m: self.m,
            delta: self.r3(delta),
            out: self.r(out),
        };
        self.go::<e::DeltaFlux>(self.nf, &args)
    }

    fn transpose(&mut self, val: Buf, out: Buf) -> Result<()> {
        self.not_on_ranks("transpose")?;
        let n = self.m.transpose.len as usize;
        ensure!(
            self.r(val).len as usize == n && self.r(out).len as usize == n,
            "transpose: matrices of {n} slots expected"
        );
        self.unary::<e::Transpose>(n, 0.0, val, out)
    }

    fn matrix_product(&mut self, val: Buf, x: Buf, out: Buf) -> Result<()> {
        self.product::<e::MatrixProduct>(val, x, out)
    }

    fn transposed_product(&mut self, val: Buf, x: Buf, out: Buf) -> Result<()> {
        self.product::<e::TransposedProduct>(val, x, out)
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
        self.x.solve(val, rhs, x, pc, cg, tol, max_iter)
    }

    fn amg_setup(&mut self, val: Buf, s: &AmgSettings) -> Result<()> {
        self.x.amg_setup(val, s, false)
    }

    fn amg_setup_slot(
        &mut self,
        slot: usize,
        val: Buf,
        s: &AmgSettings,
        symmetric_part: bool,
    ) -> Result<()> {
        self.x.amg_slot(slot)?;
        self.x.amg_setup(val, s, symmetric_part)
    }

    fn amg_select(&mut self, slot: usize) -> Result<()> {
        self.x.amg_slot(slot)
    }

    fn reduction_setup_slot(&mut self, slot: usize, val: Buf, s: &ReductionSettings) -> Result<()> {
        self.x.amg_slot(slot)?;
        self.x.reduction_setup(val, s)
    }

    fn multigrid_cycles(&mut self, val: Buf, rhs: Buf, x: Buf, cycles: usize) -> Result<()> {
        self.x.cycles(val, rhs, x, cycles)
    }

    fn face_viscosity(&self) -> Buf {
        self.mu_buf
    }

    fn set_density(&mut self, cell: &[f64], face: &[f64]) -> Result<()> {
        // on ranks the cells may come with their halo
        ensure!(
            cell.len() >= self.nc && face.len() == self.nf,
            "flow device: density fields of {} cells and {} faces expected",
            self.nc,
            self.nf
        );
        self.x.upload(self.rho_buf, cell)?;
        self.x.upload(self.rho_face_buf, face)
    }

    fn set_balanced_force(&mut self, jump: &[f64], pb: &[f64]) -> Result<()> {
        ensure!(
            jump.len() == self.nf && pb.len() == self.nb,
            "flow device: balanced-force arrays of {} faces and {} boundary faces expected",
            self.nf,
            self.nb
        );
        let (bj, bp) = match self.balanced {
            Some(b) => b,
            None => {
                let b = (self.x.alloc(self.nf)?, self.x.alloc(self.nb)?);
                self.balanced = Some(b);
                self.m.jump_offset = self.x.raw(b.0);
                self.m.pb_offset = self.x.raw(b.1);
                b
            }
        };
        self.x.upload(bj, jump)?;
        self.x.upload(bp, pb)
    }

    fn set_variable_density(
        &mut self,
        source: Option<&[f64]>,
        old: Option<&[f64]>,
        old2: Option<&[f64]>,
    ) -> Result<()> {
        let b = match self.variable {
            Some(b) => b,
            None => {
                let b = [
                    self.x.alloc(self.nc)?,
                    self.x.alloc(self.nc)?,
                    self.x.alloc(self.nc)?,
                ];
                self.variable = Some(b);
                b
            }
        };
        let mut raws = [Raw::EMPTY; 3];
        for (i, v) in [source, old, old2].into_iter().enumerate() {
            if let Some(v) = v {
                // on ranks the levels may come with their halo
                ensure!(
                    v.len() >= self.nc,
                    "flow device: a variable-density array of {} cells expected",
                    self.nc
                );
                self.x.upload(b[i], v)?;
                raws[i] = self.x.raw(b[i]);
            }
        }
        self.m.mass_source = raws[0];
        self.m.rho_old = raws[1];
        self.m.rho_old2 = raws[2];
        Ok(())
    }

    fn wall_setup(&mut self, w: &WallStatic) -> Result<()> {
        ensure!(
            w.ptr.len() == w.cells.len() + 1 && w.slots.len() == w.faces.len(),
            "flow device: inconsistent wall arrays"
        );
        let mut up = |v: &[u32]| -> Result<Raw<u32>> {
            if v.is_empty() {
                Ok(Raw::EMPTY)
            } else {
                self.x.alloc_u32(v)
            }
        };
        ensure!(
            w.sample.is_empty() || w.sample.len() == w.faces.len(),
            "flow device: inconsistent wall sample cells"
        );
        let arrays = [
            up(&w.faces)?,
            up(&w.cells)?,
            up(&w.ptr)?,
            up(&w.slots)?,
            up(&w.sample)?,
        ];
        self.wall = Some(arrays);
        Ok(())
    }

    fn sst_cells(&mut self, a: &SstCellOps) -> Result<()> {
        let mut ins = vec![a.k, a.w, a.d];
        for g in &a.grad_u {
            ins.extend(g.all());
        }
        ins.extend(a.grad_k.all());
        ins.extend(a.grad_w.all());
        if let Some((x, y)) = a.add {
            ins.push(x);
            ins.push(y);
        }
        if let Some((x, y)) = a.wall {
            ins.push(x);
            ins.push(y);
        }
        if let Some(x) = a.fr1 {
            ins.push(x);
        }
        let mut outs = vec![a.nu_t, a.f1, a.pk, a.dk, a.dw, a.capped];
        if let Some(s) = a.sources {
            outs.extend(s);
        }
        check(&outs, &ins)?;
        let src = a.sources.map(|s| s.map(|b| self.r(b)));
        let none = Raw::EMPTY;
        let args = e::SstCellArgs {
            grad_u: self.rg(&a.grad_u),
            k: self.r(a.k),
            w: self.r(a.w),
            grad_k: self.r3(a.grad_k),
            grad_w: self.r3(a.grad_w),
            d: self.r(a.d),
            nu: a.nu,
            nu_cell: a.nu_cell.map_or(Raw::EMPTY, |b| self.r(b)),
            c: a.c,
            o: a.o,
            add_k: a.add.map_or(none, |(x, _)| self.r(x)),
            add_w: a.add.map_or(none, |(_, y)| self.r(y)),
            wall_mask: a.wall.map_or(none, |(x, _)| self.r(x)),
            wall_pk: a.wall.map_or(none, |(_, y)| self.r(y)),
            fr1: a.fr1.map_or(none, |x| self.r(x)),
            nu_t: self.r(a.nu_t),
            f1: self.r(a.f1),
            pk: self.r(a.pk),
            su_k: src.map_or(none, |s| s[0]),
            sp_k: src.map_or(none, |s| s[1]),
            su_w: src.map_or(none, |s| s[2]),
            sp_w: src.map_or(none, |s| s[3]),
            dk: self.r(a.dk),
            dw: self.r(a.dw),
            capped: self.r(a.capped),
        };
        self.go::<e::SstCells>(self.nc, &args)?;
        // (ν_t and F1 are read at faces)
        self.x.exchange(&outs)
    }

    fn strain_component(&mut self, grad_u: [Buf3; 3], comp: u32, out: Buf) -> Result<()> {
        self.not_on_ranks("strain component")?;
        let ins: Vec<Buf> = grad_u.iter().flat_map(|g| g.all()).collect();
        check(&[out], &ins)?;
        let args = self.curvature_args(grad_u, out);
        let args = e::CurvatureArgs { comp, ..args };
        self.go::<e::StrainCells>(self.nc, &args)
    }

    fn convective_derivative(
        &mut self,
        u: Buf3,
        v_frame: Option<Buf3>,
        grad: Buf3,
        old: Option<(Buf, Buf, f64)>,
        out: Buf,
    ) -> Result<()> {
        self.not_on_ranks("convective derivative")?;
        let mut ins = u.all().to_vec();
        ins.extend(grad.all());
        if let Some(v) = v_frame {
            ins.extend(v.all());
        }
        if let Some((a, b, _)) = old {
            ins.push(a);
            ins.push(b);
        }
        check(&[out], &ins)?;
        let none = Raw::EMPTY;
        let base = self.curvature_args([grad, grad, grad], out);
        let args = e::CurvatureArgs {
            u: self.r3(u),
            v_frame: v_frame.map_or(Raw3::default(), |v| self.r3(v)),
            grad: self.r3(grad),
            s: old.map_or(none, |(a, _, _)| self.r(a)),
            s_old: old.map_or(none, |(_, b, _)| self.r(b)),
            inv_dt: old.map_or(0.0, |(_, _, x)| x),
            ..base
        };
        self.go::<e::ConvectiveDerivative>(self.nc, &args)
    }

    fn curvature_cells(&mut self, a: &CurvatureOps) -> Result<()> {
        self.not_on_ranks("curvature cells")?;
        let mut ins: Vec<Buf> = a.grad_u.iter().flat_map(|g| g.all()).collect();
        ins.extend(a.dsdt);
        ins.push(a.w);
        if let Some(o) = a.omega_rot {
            ins.extend(o.all());
        }
        check(&[a.out], &ins)?;
        let base = self.curvature_args(a.grad_u, a.out);
        let args = e::CurvatureArgs {
            dsdt: a.dsdt.map(|b| self.r(b)),
            w: self.r(a.w),
            omega_rot: a.omega_rot.map_or(Raw3::default(), |o| self.r3(o)),
            rot: a.rot,
            spin: a.spin,
            c: a.c,
            smooth: a.smooth,
            ..base
        };
        self.go::<e::CurvatureCells>(self.nc, &args)
    }

    fn qcr_cells(
        &mut self,
        grad_u: [Buf3; 3],
        nu_t: Buf,
        c: &foehncraft_kernels::turb::sst::SstConstants,
        mode: u32,
        tau: [Buf; 6],
    ) -> Result<()> {
        self.not_on_ranks("qcr cells")?;
        let mut ins: Vec<Buf> = grad_u.iter().flat_map(|g| g.all()).collect();
        ins.push(nu_t);
        check(&tau, &ins)?;
        let args = e::QcrArgs {
            m: self.m,
            grad_u: self.rg(&grad_u),
            nu_t: self.r(nu_t),
            c: *c,
            mode,
            tau: tau.map(|b| self.r(b)),
            bf: Raw3::default(),
            gtau: [Raw3::default(); 6],
        };
        self.go::<e::QcrCells>(self.nc, &args)
    }

    fn explicit_stress_faces(
        &mut self,
        tau: [Buf; 6],
        gtau: Option<[Buf3; 6]>,
        bf: Buf3,
    ) -> Result<()> {
        self.not_on_ranks("explicit stress faces")?;
        let mut ins = tau.to_vec();
        if let Some(g) = gtau {
            ins.extend(g.iter().flat_map(|b| b.all()));
        }
        check(&bf.all(), &ins)?;
        let args = e::QcrArgs {
            m: self.m,
            grad_u: RawGrad::default(),
            nu_t: Raw::EMPTY,
            c: foehncraft_kernels::turb::sst::SstConstants::default(),
            mode: 0,
            tau: tau.map(|b| self.r(b)),
            bf: self.r3(bf),
            gtau: gtau.map_or([Raw3::default(); 6], |g| g.map(|b| self.r3(b))),
        };
        self.go::<e::ExplicitStressFaces>(self.nf, &args)
    }

    fn vertex_setup(&mut self, s: &VertexStencil) -> Result<usize> {
        let nv = s.vert_ptr.len().saturating_sub(1);
        ensure!(
            nv > 0
                && s.bface_ptr.len() == nv + 1
                && s.cell_ptr.len() == self.nc + 1
                && s.points.iter().all(|p| p.len() == nv)
                && s.vert_ptr[nv] as usize == s.vert_cells.len()
                && s.bface_ptr[nv] as usize == s.bfaces.len()
                && s.cell_ptr[self.nc] as usize == s.cell_verts.len(),
            "flow device: inconsistent vertex stencil"
        );
        let lists = [
            self.x.alloc_u32(&s.vert_ptr)?,
            self.x.alloc_u32(&s.vert_cells)?,
            self.x.alloc_u32(&s.bface_ptr)?,
            self.x.alloc_u32(&s.bfaces)?,
            self.x.alloc_u32(&s.cell_ptr)?,
            self.x.alloc_u32(&s.cell_verts)?,
        ];
        // vertex arrays (on ranks neither cell nor face fields)
        let points = Buf3 {
            x: self.x.alloc_plain(nv)?,
            y: self.x.alloc_plain(nv)?,
            z: self.x.alloc_plain(nv)?,
        };
        self.x.upload(points.x, &s.points[0])?;
        self.x.upload(points.y, &s.points[1])?;
        self.x.upload(points.z, &s.points[2])?;
        let (vmin, vmax) = (self.x.alloc_plain(nv)?, self.x.alloc_plain(nv)?);
        // on ranks: where the other ranks hold the shared vertices (their
        // global ids, gathered once)
        let shared = match self.x.comm().filter(|_| self.x.on_ranks()) {
            Some(c) => {
                ensure!(
                    s.shared.len() == s.shared_id.len(),
                    "flow device: inconsistent shared vertices"
                );
                // the ids as f64: exact below 2⁵³
                let ids: Vec<f64> = s.shared_id.iter().map(|&i| i as f64).collect();
                let all = c.all_gather_vec_f64(&ids);
                let me = c.rank();
                let mut at: std::collections::HashMap<u64, Vec<(u32, u32)>> =
                    std::collections::HashMap::new();
                for (r, list) in all.iter().enumerate() {
                    if r == me {
                        continue;
                    }
                    for (j, &id) in list.iter().enumerate() {
                        at.entry(id as u64).or_default().push((r as u32, j as u32));
                    }
                }
                Some(VertShared {
                    local: s.shared.clone(),
                    others: s
                        .shared_id
                        .iter()
                        .map(|id| at.remove(id).unwrap_or_default())
                        .collect(),
                })
            }
            None => None,
        };
        self.vert.push(VertSlot {
            lists,
            points,
            vmin,
            vmax,
            n_vert: nv,
        });
        self.vert_shared.push(shared);
        Ok(self.vert.len() - 1)
    }

    fn mlp_limiter(&mut self, a: &MlpOps) -> Result<()> {
        let Some(&v) = self.vert.get(a.slot) else {
            anyhow::bail!("flow device: no vertex stencil {}", a.slot);
        };
        let mut ins = vec![a.phi, a.bc_a, a.bc_b];
        ins.extend(a.grad.all());
        ins.extend(a.eps_cell);
        check(&[a.out], &ins)?;
        let [vp, vc, bp, bf, cp, cv] = v.lists;
        let args = e::MlpArgs {
            m: self.m,
            vert_ptr: vp,
            vert_cells: vc,
            bface_ptr: bp,
            bfaces: bf,
            cell_ptr: cp,
            cell_verts: cv,
            points: self.r3(v.points),
            phi: self.r(a.phi),
            bc_a: self.r(a.bc_a),
            bc_b: self.r(a.bc_b),
            grad: self.r3(a.grad),
            vmin: self.r(v.vmin),
            vmax: self.r(v.vmax),
            func: a.func as u64,
            eps: a.eps,
            eps_cell: a.eps_cell.map_or(Raw::EMPTY, |b| self.r(b)),
            relax: a.relax,
            out: self.r(a.out),
        };
        self.go::<e::VertexBounds>(v.n_vert, &args)?;
        if let Some(sh) = self.vert_shared.get(a.slot).cloned().flatten() {
            self.combine_bounds(&v, &sh)?;
        }
        self.go::<e::MlpLimiter>(self.nc, &args)?;
        // an own-cell loop: the halo of ψ from the owners
        self.x.exchange(&[a.out])
    }

    fn quadratic_setup(&mut self, cells: &[u32], ptr: &[u32], idx: &[u32]) -> Result<usize> {
        self.not_on_ranks("quadratic setup")?;
        ensure!(
            !cells.is_empty()
                && ptr.len() == cells.len() + 1
                && ptr.last().copied().unwrap_or(0) as usize == idx.len(),
            "flow device: empty or inconsistent quadratic stencil"
        );
        let s = [
            self.x.alloc_u32(cells)?,
            self.x.alloc_u32(ptr)?,
            self.x.alloc_u32(idx)?,
        ];
        self.quad.push(s);
        Ok(self.quad.len() - 1)
    }

    fn quadratic_gradient(
        &mut self,
        slot: usize,
        w: Buf3,
        dims: u32,
        u: Buf3,
        ub: Buf3,
        out: [Buf3; 3],
    ) -> Result<()> {
        self.not_on_ranks("quadratic gradient")?;
        let Some(&[cells, ptr, idx]) = self.quad.get(slot) else {
            anyhow::bail!("flow device: no quadratic stencil {slot}");
        };
        let mut ins = w.all().to_vec();
        ins.extend(u.all());
        ins.extend(ub.all());
        let outs: Vec<Buf> = out.iter().flat_map(|b| b.all()).collect();
        check(&outs, &ins)?;
        let args = e::QuadGradArgs {
            cells,
            ptr,
            idx,
            w: self.r3(w),
            n_cells: self.nc as u64,
            dims: dims as u64,
            u: self.r3(u),
            ub: self.r3(ub),
            out: self.rg(&out),
        };
        self.go::<e::QuadraticGradient>(cells.len as usize, &args)
    }

    fn keps_cells(&mut self, a: &KepsCellOps) -> Result<()> {
        let mut ins = vec![a.k, a.e];
        for g in &a.grad_u {
            ins.extend(g.all());
        }
        if let Some((x, y)) = a.add {
            ins.push(x);
            ins.push(y);
        }
        if let Some((x, y)) = a.wall {
            ins.push(x);
            ins.push(y);
        }
        let mut outs = vec![a.nu_t, a.c_mu, a.pk, a.dk, a.de, a.capped];
        if let Some(s) = a.sources {
            outs.extend(s);
        }
        check(&outs, &ins)?;
        let src = a.sources.map(|s| s.map(|b| self.r(b)));
        let none = Raw::EMPTY;
        let args = e::KepsCellArgs {
            grad_u: self.rg(&a.grad_u),
            k: self.r(a.k),
            e: self.r(a.e),
            nu: a.nu,
            nu_cell: a.nu_cell.map_or(Raw::EMPTY, |b| self.r(b)),
            c: a.c,
            o: a.o,
            add_k: a.add.map_or(none, |(x, _)| self.r(x)),
            add_e: a.add.map_or(none, |(_, y)| self.r(y)),
            wall_mask: a.wall.map_or(none, |(x, _)| self.r(x)),
            wall_pk: a.wall.map_or(none, |(_, y)| self.r(y)),
            nu_t: self.r(a.nu_t),
            c_mu: self.r(a.c_mu),
            pk: self.r(a.pk),
            su_k: src.map_or(none, |s| s[0]),
            sp_k: src.map_or(none, |s| s[1]),
            su_e: src.map_or(none, |s| s[2]),
            sp_e: src.map_or(none, |s| s[3]),
            dk: self.r(a.dk),
            de: self.r(a.de),
            capped: self.r(a.capped),
        };
        self.go::<e::KepsCells>(self.nc, &args)?;
        self.x.exchange(&outs)
    }

    fn zeta_f_cells(&mut self, a: &ZetaFCellOps) -> Result<()> {
        self.not_on_ranks("zeta f cells")?;
        let mut ins = vec![a.k, a.e, a.zeta, a.f];
        for g in &a.grad_u {
            ins.extend(g.all());
        }
        if let Some(x) = a.add {
            ins.extend(x);
        }
        if let Some((x, y)) = a.wall {
            ins.push(x);
            ins.push(y);
        }
        if let Some((gk, gz, g2, dkt)) = a.eb {
            ins.extend(gk.all());
            ins.extend(gz.all());
            ins.extend([g2, dkt]);
        }
        let mut outs = vec![a.nu_t, a.t, a.l2, a.pk, a.dk, a.de, a.dz, a.capped];
        if let Some(s) = a.sources {
            outs.extend(s);
        }
        check(&outs, &ins)?;
        let src = a.sources.map(|s| s.map(|b| self.r(b)));
        let add = a.add.map(|s| s.map(|b| self.r(b)));
        let none = Raw::EMPTY;
        let args = e::ZetaFCellArgs {
            grad_u: self.rg(&a.grad_u),
            k: self.r(a.k),
            e: self.r(a.e),
            zeta: self.r(a.zeta),
            f: self.r(a.f),
            nu: a.nu,
            nu_cell: a.nu_cell.map_or(Raw::EMPTY, |b| self.r(b)),
            c: a.c,
            o: a.o,
            add_k: add.map_or(none, |s| s[0]),
            add_e: add.map_or(none, |s| s[1]),
            add_z: add.map_or(none, |s| s[2]),
            add_f: add.map_or(none, |s| s[3]),
            wall_mask: a.wall.map_or(none, |(x, _)| self.r(x)),
            wall_pk: a.wall.map_or(none, |(_, y)| self.r(y)),
            nu_t: self.r(a.nu_t),
            t: self.r(a.t),
            l2: self.r(a.l2),
            pk: self.r(a.pk),
            su_k: src.map_or(none, |s| s[0]),
            sp_k: src.map_or(none, |s| s[1]),
            su_e: src.map_or(none, |s| s[2]),
            sp_e: src.map_or(none, |s| s[3]),
            su_z: src.map_or(none, |s| s[4]),
            sp_z: src.map_or(none, |s| s[5]),
            rhs_f: src.map_or(none, |s| s[6]),
            sp_f: src.map_or(none, |s| s[7]),
            dk: self.r(a.dk),
            de: self.r(a.de),
            dz: self.r(a.dz),
            capped: self.r(a.capped),
            grad_k: a.eb.map_or(Raw3::default(), |b| self.r3(b.0)),
            grad_z: a.eb.map_or(Raw3::default(), |b| self.r3(b.1)),
            g2: a.eb.map_or(none, |b| self.r(b.2)),
            dkt: a.eb.map_or(none, |b| self.r(b.3)),
            dist: self.r(a.dist),
        };
        self.go::<e::ZetaFCells>(self.nc, &args)
    }

    fn turb_faces(&mut self, a: &TurbFaceOps) -> Result<()> {
        let mut outs = vec![a.mu_eff, a.gamma_k, a.gamma_w];
        let mut ins = vec![a.nu_t, a.dk, a.dw, a.mu, a.wall];
        if let Some([dz, l2, gz, gf]) = a.zeta_f {
            ins.extend([dz, l2]);
            outs.extend([gz, gf]);
        }
        check(&outs, &ins)?;
        let none = Raw::EMPTY;
        let args = e::TurbFaceArgs {
            m: self.m,
            nu_t: self.r(a.nu_t),
            dk: self.r(a.dk),
            dw: self.r(a.dw),
            mu: self.r(a.mu),
            wall: self.r(a.wall),
            mol: a.mol,
            mu_eff: self.r(a.mu_eff),
            gamma_k: self.r(a.gamma_k),
            gamma_w: self.r(a.gamma_w),
            dz: a.zeta_f.map_or(none, |z| self.r(z[0])),
            l2: a.zeta_f.map_or(none, |z| self.r(z[1])),
            gamma_z: a.zeta_f.map_or(none, |z| self.r(z[2])),
            gamma_f: a.zeta_f.map_or(none, |z| self.r(z[3])),
        };
        self.go::<e::TurbFaces>(self.nf, &args)
    }

    fn scalar_turb_faces(&mut self, a: &ScalarTurbOps) -> Result<()> {
        check(
            &[a.gamma, a.bc_b],
            &[a.mu, a.mu_eff, a.gamma0, a.wall, a.fq, a.bc_a],
        )?;
        let args = e::ScalarTurbArgs {
            m: self.m,
            mu: self.r(a.mu),
            mu_eff: self.r(a.mu_eff),
            gamma0: self.r(a.gamma0),
            wall: self.r(a.wall),
            fq: self.r(a.fq),
            bc_a: self.r(a.bc_a),
            prt_model: a.prt.0,
            prt: a.prt.1,
            gamma: self.r(a.gamma),
            bc_b: self.r(a.bc_b),
        };
        self.go::<e::ScalarTurbFaces>(self.nf, &args)
    }

    fn scalar_wall_faces(&mut self, a: &ScalarWallOps) -> Result<()> {
        let Some([faces, ..]) = self.wall else {
            anyhow::bail!("flow device: scalar_wall_faces before wall_setup");
        };
        let mut outs = vec![a.gamma, a.bc_b];
        outs.extend(a.t_plus);
        check(&outs, &[a.y_plus, a.mu, a.gamma0, a.fq, a.bc_a])?;
        let args = e::ScalarWallArgs {
            m: self.m,
            faces,
            y_plus: self.r(a.y_plus),
            mu: self.r(a.mu),
            gamma0: self.r(a.gamma0),
            fq: self.r(a.fq),
            bc_a: self.r(a.bc_a),
            law: a.law,
            prt: a.prt,
            y_floor: a.y_floor,
            gamma: self.r(a.gamma),
            bc_b: self.r(a.bc_b),
            t_plus: a.t_plus.map_or(Raw::EMPTY, |b| self.r(b)),
        };
        self.go::<e::ScalarWallFaces>(faces.len as usize, &args)
    }

    fn wall_faces(&mut self, a: &WallFaceOps) -> Result<()> {
        let Some([faces, _, _, _, sample]) = self.wall else {
            anyhow::bail!("flow device: wall_faces before wall_setup");
        };
        let mut ins = vec![a.y, a.ks, a.k, a.cs];
        if let Some(y) = a.y_sample {
            ins.push(y);
        }
        let sample = if a.y_sample.is_some() {
            sample
        } else {
            Raw::EMPTY
        };
        ins.extend(a.normal.all());
        ins.extend(a.u.all());
        ins.extend(a.ub.all());
        let mut outs = vec![
            a.y_plus,
            a.y_star,
            a.omega,
            a.blending,
            a.viscous_weight,
            a.mu_eff,
            a.k_b,
            a.w_b,
            a.u_tau,
            a.pk,
            a.prod_w,
        ];
        if let Some((z, fb)) = a.zeta_f {
            ins.push(z);
            outs.push(fb);
        }
        if let Some((gu, gp)) = &a.non_eq {
            for g in gu {
                ins.extend(g.all());
            }
            ins.extend(gp.all());
        }
        check(&outs, &ins)?;
        let args = e::WallFaceArgs {
            m: self.m,
            faces,
            normal: self.r3(a.normal),
            y: self.r(a.y),
            ks: self.r(a.ks),
            u: self.r3(a.u),
            ub: self.r3(a.ub),
            nu: a.nu,
            nu_cell: a.nu_cell.map_or(Raw::EMPTY, |b| self.r(b)),
            treatment: a.treatment,
            c: a.c,
            ke: a.ke,
            k: self.r(a.k),
            cs: self.r(a.cs),
            u_tau_old: self.r(a.u_tau),
            u_tau: self.r(a.u_tau),
            y_plus: self.r(a.y_plus),
            y_star: self.r(a.y_star),
            omega: self.r(a.omega),
            blending: self.r(a.blending),
            viscous_weight: self.r(a.viscous_weight),
            pk: self.r(a.pk),
            prod_w: self.r(a.prod_w),
            mu_eff: self.r(a.mu_eff),
            k_b: self.r(a.k_b),
            w_b: self.r(a.w_b),
            zeta: a.zeta_f.map_or(Raw::EMPTY, |(z, _)| self.r(z)),
            f_b: a.zeta_f.map_or(Raw::EMPTY, |(_, f)| self.r(f)),
            n_zeta: a.n_zeta,
            grad_u: a.non_eq.map_or(Default::default(), |(gu, _)| self.rg(&gu)),
            grad_p: a.non_eq.map_or(Default::default(), |(_, gp)| self.r3(gp)),
            sample,
            y_sample: a.y_sample.map_or(Raw::EMPTY, |b| self.r(b)),
        };
        self.go::<e::WallFaces>(faces.len as usize, &args)
    }

    fn floor(
        &mut self,
        x: Buf,
        lo: f64,
        hi: f64,
        rel_old: Option<(Buf, f64)>,
        work: Buf,
    ) -> Result<usize> {
        check(&[work], &[x])?;
        if let Some((old, _)) = rel_old {
            check(&[work], &[old])?;
            ensure!(old != x, "flow device: floor with old = x");
        }
        let args = e::FloorArgs {
            x: self.r(x),
            lo,
            hi,
            rel: rel_old.map_or(0.0, |r| r.1),
            old: rel_old.map_or(Raw::EMPTY, |r| self.r(r.0)),
            hit: self.r(work),
        };
        let n = self.r(x).len as usize;
        self.go::<e::Floor>(n, &args)?;
        Ok(self.sum_sq(work)?.round() as usize)
    }

    fn wall_cells(&mut self, a: &WallCellOps) -> Result<()> {
        let Some([_, cells, ptr, slots, _]) = self.wall else {
            anyhow::bail!("flow device: wall_cells before wall_setup");
        };
        check(
            &[a.fix_val, a.fix_mask, a.wall_pk, a.prod_mask],
            &[a.omega, a.weight, a.blending, a.pk, a.prod_w],
        )?;
        let args = e::WallCellArgs {
            cells,
            ptr,
            slots,
            weight: self.r(a.weight),
            omega: self.r(a.omega),
            blending: self.r(a.blending),
            pk: self.r(a.pk),
            prod_w: self.r(a.prod_w),
            wall_pk: self.r(a.wall_pk),
            prod_mask: self.r(a.prod_mask),
            fix_val: self.r(a.fix_val),
            fix_mask: self.r(a.fix_mask),
            phi: a.phi.map_or(Raw::EMPTY, |b| self.r(b)),
        };
        self.go::<e::WallCells>(cells.len as usize, &args)?;
        let mut outs = vec![a.fix_val, a.fix_mask, a.wall_pk, a.prod_mask];
        outs.extend(a.phi);
        self.x.exchange(&outs)
    }

    fn sync(&mut self) -> Result<()> {
        self.x.sync()
    }
}
