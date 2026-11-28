//! Launch entries of the flow solver (stage 5): one argument block and one
//! per-index body for every operation of the flow loop, shared by the CPU
//! loops (`foehncraft-flow`, `device::kernel`) and the PTX entry points
//! (`gpu/ptx-kernels/src/flow.rs`, entry names `flow_*`).
//!
//! An argument block is a `#[repr(C)]` struct of device addresses with
//! their lengths ([`Raw`]) and scalars. The host fills it with the
//! addresses of its buffers (host memory on the CPU, device memory on the
//! GPU) and passes it by value: as a kernel parameter on the GPU (the
//! NVPTX back end passes aggregates like CUDA C, `.param .align 8 .b8`), as
//! a reference into the loop body on the CPU. Both sides compile this one
//! definition, so the argument lists cannot drift apart, and both run the
//! same body for an index, so CPU and GPU agree bitwise (D31).
//!
//! The bodies call the kernel bodies of [`crate::flow`], [`crate::fv`] and
//! [`crate::blas`] and write their outputs; each output is written by one
//! index only. Inputs and outputs of one launch are distinct buffers, except
//! where an entry says that an output may be an input (element-wise
//! updates, which read index i before they write it, through raw pointers).

use crate::adjoint;
use crate::flow::{self as k, FaceGeometry, VelocityGradient};
use crate::fv::{self, Boundary, SellLayout, Topology, Vec3};
use crate::{blas, ld};
use core::marker::PhantomData;

/// A device array: address and length. The address is that of host memory
/// for the CPU loops and of GPU memory for the PTX entries.
#[repr(C)]
pub struct Raw<T> {
    /// Address of the first entry.
    pub ptr: u64,
    /// Number of entries.
    pub len: u64,
    _t: PhantomData<fn() -> T>,
}

impl<T> Clone for Raw<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Raw<T> {}

impl<T> Default for Raw<T> {
    fn default() -> Self {
        Self::EMPTY
    }
}

impl<T> core::fmt::Debug for Raw<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Raw({:#x}, {})", self.ptr, self.len)
    }
}

impl<T> Raw<T> {
    /// No array.
    pub const EMPTY: Self = Self {
        ptr: 0,
        len: 0,
        _t: PhantomData,
    };

    /// An array at `ptr` with `len` entries.
    pub const fn new(ptr: u64, len: usize) -> Self {
        Self {
            ptr,
            len: len as u64,
            _t: PhantomData,
        }
    }

    /// The entries as a slice (empty for [`Raw::EMPTY`]).
    ///
    /// # Safety
    /// `ptr` addresses `len` initialised entries that no thread writes while
    /// the slice is used.
    #[inline(always)]
    pub unsafe fn get<'a>(self) -> &'a [T] {
        if self.len == 0 {
            &[]
        } else {
            // SAFETY: guaranteed by the caller.
            unsafe { core::slice::from_raw_parts(self.ptr as usize as *const T, self.len as usize) }
        }
    }

    /// The address for writing.
    #[inline(always)]
    pub fn out(self) -> *mut T {
        self.ptr as usize as *mut T
    }
}

impl<T: Copy> Raw<T> {
    /// Entry `i` read through the raw address (for arrays that the same
    /// launch also writes, at the same index).
    ///
    /// # Safety
    /// `i < len`.
    #[inline(always)]
    pub unsafe fn read(self, i: usize) -> T {
        // SAFETY: guaranteed by the caller.
        unsafe { *(self.ptr as usize as *const T).add(i) }
    }
}

/// Three arrays, one per component.
#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct Raw3 {
    /// x components.
    pub x: Raw<f64>,
    /// y components.
    pub y: Raw<f64>,
    /// z components.
    pub z: Raw<f64>,
}

impl Raw3 {
    /// As a [`Vec3`].
    ///
    /// # Safety
    /// As [`Raw::get`] for all three.
    #[inline(always)]
    pub unsafe fn v<'a>(self) -> Vec3<'a, f64> {
        // SAFETY: guaranteed by the caller.
        unsafe {
            Vec3 {
                x: self.x.get(),
                y: self.y.get(),
                z: self.z.get(),
            }
        }
    }
}

/// The nine arrays of a velocity gradient.
#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct RawGrad {
    /// ∇u_x, ∇u_y, ∇u_z.
    pub g: [Raw3; 3],
}

impl RawGrad {
    /// As a [`VelocityGradient`].
    ///
    /// # Safety
    /// As [`Raw::get`] for all nine.
    #[inline(always)]
    pub unsafe fn v<'a>(self) -> VelocityGradient<'a, f64> {
        // SAFETY: guaranteed by the caller.
        unsafe {
            VelocityGradient {
                x: self.g[0].v(),
                y: self.g[1].v(),
                z: self.g[2].v(),
            }
        }
    }
}

/// Static data of a mesh: connectivity, geometry, matrix layout and the
/// boundary-condition arrays (uploaded once per run).
#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct Mesh {
    /// Number of local cells.
    pub n_cells: u64,
    /// Number of faces.
    pub n_faces: u64,
    /// Number of interior faces (boundary faces follow).
    pub n_interior: u64,
    /// Owner cell of each face.
    pub owner: Raw<u32>,
    /// Neighbour cell of each interior face.
    pub neighbour: Raw<u32>,
    /// Start of the face list of each cell.
    pub cell_ptr: Raw<u32>,
    /// Packed faces of each cell.
    pub cell_faces: Raw<u32>,
    /// Face area vectors.
    pub area: Raw3,
    /// Face centroids.
    pub face_centroid: Raw3,
    /// Cell centroids.
    pub cell_centroid: Raw3,
    /// |S|²/(d·S).
    pub dgeo: Raw<f64>,
    /// Non-orthogonal part k of the face area vector.
    pub k: Raw3,
    /// Owner weight of the face interpolation.
    pub weight: Raw<f64>,
    /// Cell volumes.
    pub volume: Raw<f64>,
    /// Least-squares weights of the velocity (per cell-face entry).
    pub lsq_u: Raw3,
    /// Least-squares weights of the pressure.
    pub lsq_p: Raw3,
    /// SELL slice height.
    pub chunk: u64,
    /// SELL slice starts.
    pub slice_ptr: Raw<u32>,
    /// Transpose map of the SELL storage (D59; see [`crate::adjoint`]).
    pub transpose: Raw<u32>,
    /// Momentum kind per boundary face.
    pub mom_kind: Raw<u32>,
    /// Flux kind per boundary face.
    pub flux_kind: Raw<u32>,
    /// Boundary velocity rule per boundary face.
    pub ub_rule: Raw<u32>,
    /// Boundary pressure rule per boundary face.
    pub pb_rule: Raw<u32>,
    /// Prescribed (static or total) pressure per boundary face.
    pub pset: Raw<f64>,
    /// Viscosity per face.
    pub mu_face: Raw<f64>,
    /// Density per cell (D150; a constant fill reproduces the constant
    /// density bit for bit).
    pub rho: Raw<f64>,
    /// Density per face, the face value of the mass flux and of the
    /// pressure-correction coefficient.
    pub rho_face: Raw<f64>,
    /// Closed state per boundary face of the backflow suppression (1
    /// closed, 0 open), written by [`CloseBackflow`] in every outer
    /// iteration.
    pub closed: Raw<f64>,
    /// Old-level part c of the convective outlets per boundary face, set
    /// per time step ([`k::UB_CONVECTIVE`]).
    pub conv: Raw3,
    /// β of the convective outlets per boundary face, set per time step.
    pub beta: Raw<f64>,
    /// First periodic face (`n_interior` without periodic faces; D29).
    pub periodic: u64,
    /// Shift of the neighbour centroid per periodic face.
    pub shift: Raw3,
    /// The balanced-force offset q_f of the face pressure jump per face
    /// (D154, stage 7b1): the Rhie–Chow jump becomes (p_N − p_P) −
    /// ∇p̄_f·d − q_f; empty without gravity (the plain jump bit for bit).
    pub jump_offset: Raw<f64>,
    /// The balanced-force term e_b of the extrapolated boundary pressure
    /// per boundary face (D154 (4)): p_b = p_P + ∇p_P·r + e_b; empty
    /// without gravity.
    pub pb_offset: Raw<f64>,
    /// The mass source V (a₀ρⁿ⁺¹ + a₁ρⁿ + a₂ρⁿ⁻¹)/Δt per cell of a
    /// variable density in a time step (stage 7c, D151): the continuity
    /// rows read Σ ṁ + source = 0; empty otherwise.
    pub mass_source: Raw<f64>,
    /// The cell density of the old time levels n and n − 1 (stage 7c);
    /// empty: the current density.
    pub rho_old: Raw<f64>,
    /// See `rho_old`.
    pub rho_old2: Raw<f64>,
    /// The body force per unit volume B_P per cell of the balanced force
    /// (N/m³) of a volume-of-fluid case (L3): with it the Rhie–Chow flux
    /// takes the D-weighted form of [`k::face_mass_flux`] and `jump_offset`
    /// holds the plain face offsets o_f; empty otherwise (bit for bit the
    /// former flux).
    pub body: Raw3,
}

impl Mesh {
    /// Connectivity.
    ///
    /// # Safety
    /// As [`Raw::get`].
    #[inline(always)]
    pub unsafe fn topo<'a>(&self) -> Topology<'a> {
        // SAFETY: guaranteed by the caller.
        unsafe {
            Topology {
                owner: self.owner.get(),
                neighbour: self.neighbour.get(),
                cell_ptr: self.cell_ptr.get(),
                cell_faces: self.cell_faces.get(),
                n_interior: self.n_interior as usize,
                periodic: self.periodic as usize,
            }
        }
    }

    /// Geometry.
    ///
    /// # Safety
    /// As [`Raw::get`].
    #[inline(always)]
    pub unsafe fn geo<'a>(&self) -> FaceGeometry<'a, f64> {
        // SAFETY: guaranteed by the caller.
        unsafe {
            FaceGeometry {
                area: self.area.v(),
                face_centroid: self.face_centroid.v(),
                cell_centroid: self.cell_centroid.v(),
                dgeo: self.dgeo.get(),
                k: self.k.v(),
                weight: self.weight.get(),
                shift: self.shift.v(),
            }
        }
    }

    /// Matrix layout.
    ///
    /// # Safety
    /// As [`Raw::get`].
    #[inline(always)]
    pub unsafe fn layout<'a>(&self) -> SellLayout<'a> {
        SellLayout {
            chunk: self.chunk as usize,
            // SAFETY: guaranteed by the caller.
            slice_ptr: unsafe { self.slice_ptr.get() },
        }
    }

    /// Number of boundary faces.
    #[inline(always)]
    pub fn n_boundary(&self) -> usize {
        (self.n_faces - self.n_interior) as usize
    }
}

/// An operation of the flow loop: its argument block and its body for one
/// index. The PTX entry `flow_<NAME>` runs [`Entry::run`] for every index
/// below the count given at launch.
pub trait Entry {
    /// Argument block.
    type Args: Copy + Send + Sync;
    /// Name of the operation (PTX entry `flow_<NAME>`).
    const NAME: &'static str;
    /// The body for index `i`.
    ///
    /// # Safety
    /// `i` is below the count of the launch, every address in `a` is valid
    /// for its length and the outputs are written by no other launch at the
    /// same time.
    unsafe fn run(i: usize, a: &Self::Args);
}

// SAFETY (of the impls below): the argument blocks hold plain addresses and
// numbers; sharing them between threads is what a launch does.
macro_rules! send_sync {
    ($($t:ty),*) => {$(
        // SAFETY: see above.
        unsafe impl Send for $t {}
        // SAFETY: see above.
        unsafe impl Sync for $t {}
    )*};
}

/// Arguments of the element-wise entries.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct LinArgs {
    /// a.
    pub a: f64,
    /// b.
    pub b: f64,
    /// x.
    pub x: Raw<f64>,
    /// y.
    pub y: Raw<f64>,
    /// Output (may be `x` or `y`).
    pub out: Raw<f64>,
}

/// out = a x + b y.
pub struct Lincomb;
impl Entry for Lincomb {
    type Args = LinArgs;
    const NAME: &'static str = "lincomb";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe { *a.out.out().add(i) = a.a * a.x.read(i) + a.b * a.y.read(i) }
    }
}

/// out = a (x − y).
pub struct ScaledDifference;
impl Entry for ScaledDifference {
    type Args = LinArgs;
    const NAME: &'static str = "scaled_difference";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe { *a.out.out().add(i) = a.a * (a.x.read(i) - a.y.read(i)) }
    }
}

/// out = min(y, 1/(a/x + b)): D̂ of the coupled preconditioner
/// ([`crate::flow::schur_diagonal`], x = D^s, y = D^C).
pub struct SchurDiagonal;
impl Entry for SchurDiagonal {
    type Args = LinArgs;
    const NAME: &'static str = "schur_diagonal";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe {
            *a.out.out().add(i) = crate::flow::schur_diagonal(a.a, a.b, a.x.read(i), a.y.read(i))
        }
    }
}

/// out = a x y per entry (the density per cell of the coupled
/// preconditioner with several fluids, D194).
pub struct Product;
impl Entry for Product {
    type Args = LinArgs;
    const NAME: &'static str = "product";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe { *a.out.out().add(i) = a.a * a.x.read(i) * a.y.read(i) }
    }
}

/// Arguments of the element-wise entries on three-component arrays.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct VecArgs {
    /// A factor per entry ([`ScaleVector`]).
    pub s: Raw<f64>,
    /// The three components.
    pub v: Raw3,
    /// Output of [`ScaleVector`] (may be `v`).
    pub out3: Raw3,
    /// Output of [`SquaresAdd`].
    pub out: Raw<f64>,
}

/// out += v_x² + v_y² + v_z².
pub struct SquaresAdd;
impl Entry for SquaresAdd {
    type Args = VecArgs;
    const NAME: &'static str = "squares_add";
    #[inline(always)]
    unsafe fn run(i: usize, a: &VecArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe {
            let (x, y, z) = (a.v.x.read(i), a.v.y.read(i), a.v.z.read(i));
            *a.out.out().add(i) = a.out.read(i) + (x * x + y * y + z * z);
        }
    }
}

/// out3 = s v.
pub struct ScaleVector;
impl Entry for ScaleVector {
    type Args = VecArgs;
    const NAME: &'static str = "scale_vector";
    #[inline(always)]
    unsafe fn run(i: usize, a: &VecArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe {
            let s = a.s.read(i);
            *a.out3.x.out().add(i) = s * a.v.x.read(i);
            *a.out3.y.out().add(i) = s * a.v.y.read(i);
            *a.out3.z.out().add(i) = s * a.v.z.read(i);
        }
    }
}

/// Arguments of [`PseudoCoefficient`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct PseudoArgs {
    /// Mesh (cell volumes).
    pub m: Mesh,
    /// 1/Δτ_max.
    pub inv_dt: f64,
    /// Factor f of the rate f · num / den.
    pub factor: f64,
    /// Numerator per cell (empty: 1).
    pub num: Raw<f64>,
    /// Denominator per cell (empty: 1).
    pub den: Raw<f64>,
    /// Output: ρ V/Δτ per cell.
    pub out: Raw<f64>,
}

/// Pseudo-time coefficient of a cell
/// ([`crate::scalar::pseudo_coefficient`]).
pub struct PseudoCoefficient;
impl Entry for PseudoCoefficient {
    type Args = PseudoArgs;
    const NAME: &'static str = "pseudo_coefficient";
    #[inline(always)]
    unsafe fn run(c: usize, a: &PseudoArgs) {
        // SAFETY: guaranteed by the caller; index c only.
        unsafe {
            let (num, den) = (a.num.get(), a.den.get());
            *a.out.out().add(c) = crate::scalar::pseudo_coefficient(
                ld(a.m.rho.get(), c) * ld(a.m.volume.get(), c),
                a.inv_dt,
                a.factor,
                if num.is_empty() { 1.0 } else { ld(num, c) },
                if den.is_empty() { 1.0 } else { ld(den, c) },
            );
        }
    }
}

/// out = a (every entry).
pub struct Fill;
impl Entry for Fill {
    type Args = LinArgs;
    const NAME: &'static str = "fill";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(i) = a.a }
    }
}

/// out −= a.
pub struct Subtract;
impl Entry for Subtract {
    type Args = LinArgs;
    const NAME: &'static str = "subtract";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe { *a.out.out().add(i) = a.out.read(i) - a.a }
    }
}

/// out\[b as index\] *= a (launched with count 1).
pub struct ScaleEntry;
impl Entry for ScaleEntry {
    type Args = LinArgs;
    const NAME: &'static str = "scale_entry";
    #[inline(always)]
    unsafe fn run(_: usize, a: &LinArgs) {
        let j = a.b as usize;
        // SAFETY: guaranteed by the caller (j inside out).
        unsafe { *a.out.out().add(j) = a.out.read(j) * a.a }
    }
}

/// Arguments of the block reductions.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct BlockArgs {
    /// The array.
    pub x: Raw<f64>,
    /// Block sums (or minima), one per block of `foehncraft_core::repro::BLOCK`.
    pub out: Raw<f64>,
    /// Block maxima (min_max only).
    pub out2: Raw<f64>,
}

/// Σ x² per block (the leaves of the fixed tree of D31).
pub struct BlockSumSq;
impl Entry for BlockSumSq {
    type Args = BlockArgs;
    const NAME: &'static str = "block_sum_sq";
    #[inline(always)]
    unsafe fn run(b: usize, a: &BlockArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let x = a.x.get();
            *a.out.out().add(b) = blas::block_dot(b, x, x);
        }
    }
}

/// Arguments of [`BlockDot`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct DotArgs {
    /// First array.
    pub x: Raw<f64>,
    /// Second array (as long as `x`).
    pub y: Raw<f64>,
    /// Block sums, one per block of `foehncraft_core::repro::BLOCK`.
    pub out: Raw<f64>,
}

/// Σ x y per block (the leaves of the fixed tree of D31).
pub struct BlockDot;
impl Entry for BlockDot {
    type Args = DotArgs;
    const NAME: &'static str = "block_dot";
    #[inline(always)]
    unsafe fn run(b: usize, a: &DotArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(b) = blas::block_dot(b, a.x.get(), a.y.get());
        }
    }
}

/// Minimum and maximum per block.
pub struct BlockMinMax;
impl Entry for BlockMinMax {
    type Args = BlockArgs;
    const NAME: &'static str = "block_min_max";
    #[inline(always)]
    unsafe fn run(b: usize, a: &BlockArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let (mn, mx) = blas::block_min_max(b, a.x.get());
            *a.out.out().add(b) = mn;
            *a.out2.out().add(b) = mx;
        }
    }
}

/// Arguments of [`BoundaryConstants`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ConstArgs {
    /// Zero-gradient mask per boundary face.
    pub mask: Raw<f64>,
    /// Boundary values.
    pub value: Raw<f64>,
    /// Output per boundary face.
    pub out: Raw<f64>,
}

/// Boundary constants of a gradient ([`k::boundary_constant`]).
pub struct BoundaryConstants;
impl Entry for BoundaryConstants {
    type Args = ConstArgs;
    const NAME: &'static str = "boundary_constants";
    #[inline(always)]
    unsafe fn run(fb: usize, a: &ConstArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(fb) = k::boundary_constant(fb, a.mask.get(), a.value.get()) }
    }
}

/// Arguments of [`BoundaryValues`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct BoundaryValuesArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face mass flux.
    pub mass_flux: Raw<f64>,
    /// Cell velocity.
    pub u: Raw3,
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// Cell pressure.
    pub p: Raw<f64>,
    /// Pressure gradient.
    pub grad_p: Raw3,
    /// Prescribed boundary velocity.
    pub ub_set: Raw3,
    /// Output: boundary velocity.
    pub ub: Raw3,
    /// Output: boundary pressure.
    pub pb: Raw<f64>,
}

/// New boundary velocity and pressure of boundary face `fb`.
pub struct BoundaryValues;
impl Entry for BoundaryValues {
    type Args = BoundaryValuesArgs;
    const NAME: &'static str = "boundary_values";
    #[inline(always)]
    unsafe fn run(fb: usize, a: &BoundaryValuesArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::BoundaryInputs {
                topo: a.m.topo(),
                geo: a.m.geo(),
                rho_face: a.m.rho_face.get(),
                mass_flux: a.mass_flux.get(),
            };
            let (u, g, ub0) = (a.u.v(), a.grad_u.v(), a.ub_set.v());
            let rule = ld(a.m.ub_rule.get(), fb);
            let conv = if rule == k::UB_CONVECTIVE {
                (a.m.conv.v().at(fb), ld(a.m.beta.get(), fb))
            } else {
                ([0.0; 3], 0.0)
            };
            let v = k::boundary_velocity(fb, &inp, rule, ub0.at(fb), conv, &u, &g);
            *a.ub.x.out().add(fb) = v[0];
            *a.ub.y.out().add(fb) = v[1];
            *a.ub.z.out().add(fb) = v[2];
            let (p, gp) = (a.p.get(), a.grad_p.v());
            *a.pb.out().add(fb) = k::boundary_pressure(
                fb,
                &inp,
                ld(a.m.pb_rule.get(), fb),
                ld(a.m.pset.get(), fb),
                p,
                &gp,
                a.m.pb_offset.get(),
            );
        }
    }
}

/// Arguments of [`Gradient`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct GradientArgs {
    /// Mesh.
    pub m: Mesh,
    /// Least-squares weights.
    pub w: Raw3,
    /// Cell values.
    pub phi: Raw<f64>,
    /// Boundary coefficient a (the zero-gradient mask).
    pub bc_a: Raw<f64>,
    /// Boundary constant b.
    pub bc_b: Raw<f64>,
    /// Jump per periodic face (empty: plain periodic).
    pub jump: Raw<f64>,
    /// Output gradient.
    pub out: Raw3,
}

/// Least-squares gradient of cell `c` ([`fv::lsq_gradient`]).
pub struct Gradient;
impl Entry for Gradient {
    type Args = GradientArgs;
    const NAME: &'static str = "gradient";
    #[inline(always)]
    unsafe fn run(c: usize, a: &GradientArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let bc = Boundary {
                a: a.bc_a.get(),
                b: a.bc_b.get(),
                jump: a.jump.get(),
            };
            let g = fv::lsq_gradient(c, &a.m.topo(), &a.w.v(), a.phi.get(), &bc);
            *a.out.x.out().add(c) = g[0];
            *a.out.y.out().add(c) = g[1];
            *a.out.z.out().add(c) = g[2];
        }
    }
}

/// Arguments of [`Limiter`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct LimiterArgs {
    /// Mesh.
    pub m: Mesh,
    /// Cell values.
    pub phi: Raw<f64>,
    /// Boundary coefficient a.
    pub bc_a: Raw<f64>,
    /// Boundary constant b.
    pub bc_b: Raw<f64>,
    /// Jump per periodic face (empty: plain periodic).
    pub jump: Raw<f64>,
    /// Gradient.
    pub grad: Raw3,
    /// ε² per cell or one value.
    pub eps2: Raw<f64>,
    /// Output ψ.
    pub out: Raw<f64>,
}

/// Venkatakrishnan limiter of cell `c` ([`fv::venkat_limiter`]).
pub struct Limiter;
impl Entry for Limiter {
    type Args = LimiterArgs;
    const NAME: &'static str = "limiter";
    #[inline(always)]
    unsafe fn run(c: usize, a: &LimiterArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let bc = Boundary {
                a: a.bc_a.get(),
                b: a.bc_b.get(),
                jump: a.jump.get(),
            };
            *a.out.out().add(c) = fv::venkat_limiter(
                c,
                &a.m.topo(),
                a.phi.get(),
                &bc,
                &a.grad.v(),
                &a.m.cell_centroid.v(),
                &a.m.face_centroid.v(),
                a.eps2.get(),
            );
        }
    }
}

/// Arguments of [`VertexBounds`] and [`MlpLimiter`] (the
/// multi-dimensional limiting process, Park & Kim 2012).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct MlpArgs {
    /// Mesh (topology, cell centroids).
    pub m: Mesh,
    /// Start of the cells of each vertex (one more entry than vertices).
    pub vert_ptr: Raw<u32>,
    /// The cells of the vertices.
    pub vert_cells: Raw<u32>,
    /// Start of the boundary faces of each vertex.
    pub bface_ptr: Raw<u32>,
    /// The boundary faces (boundary face indices) of the vertices.
    pub bfaces: Raw<u32>,
    /// Start of the vertices of each cell.
    pub cell_ptr: Raw<u32>,
    /// The vertices of the cells.
    pub cell_verts: Raw<u32>,
    /// Vertex coordinates.
    pub points: Raw3,
    /// Cell values.
    pub phi: Raw<f64>,
    /// Boundary coefficient a.
    pub bc_a: Raw<f64>,
    /// Boundary constant b.
    pub bc_b: Raw<f64>,
    /// Gradient.
    pub grad: Raw3,
    /// Bounds per vertex (written by [`VertexBounds`]).
    pub vmin: Raw<f64>,
    /// Bounds per vertex.
    pub vmax: Raw<f64>,
    /// Limiter function ([`fv::mlp_limiter`]).
    pub func: u64,
    /// Threshold ε in the units of φ.
    pub eps: f64,
    /// Factor of ε per cell (empty: none).
    pub eps_cell: Raw<f64>,
    /// Relaxation β of ψ between calls, ψ ← ψ + β (ψ_new − ψ) (1: the new
    /// value).
    pub relax: f64,
    /// Output ψ per cell (read as the last ψ when `relax` < 1).
    pub out: Raw<f64>,
}

/// The bounds of every vertex over its MLP stencil ([`fv::vertex_bounds`]).
pub struct VertexBounds;
impl Entry for VertexBounds {
    type Args = MlpArgs;
    const NAME: &'static str = "vertex_bounds";
    #[inline(always)]
    unsafe fn run(v: usize, a: &MlpArgs) {
        // SAFETY: guaranteed by the caller; index v only.
        unsafe {
            let bc = Boundary {
                a: a.bc_a.get(),
                b: a.bc_b.get(),
                jump: &[],
            };
            let (mn, mx) = fv::vertex_bounds(
                v,
                a.vert_ptr.get(),
                a.vert_cells.get(),
                a.bface_ptr.get(),
                a.bfaces.get(),
                &a.m.topo(),
                a.phi.get(),
                &bc,
            );
            *a.vmin.out().add(v) = mn;
            *a.vmax.out().add(v) = mx;
        }
    }
}

/// The MLP limiter of cell `c` ([`fv::mlp_limiter`]).
pub struct MlpLimiter;
impl Entry for MlpLimiter {
    type Args = MlpArgs;
    const NAME: &'static str = "mlp_limiter";
    #[inline(always)]
    unsafe fn run(c: usize, a: &MlpArgs) {
        // SAFETY: guaranteed by the caller; index c only.
        unsafe {
            let psi = fv::mlp_limiter(
                c,
                a.cell_ptr.get(),
                a.cell_verts.get(),
                &a.points.v(),
                a.vmin.get(),
                a.vmax.get(),
                a.phi.get(),
                &a.grad.v(),
                &a.m.cell_centroid.v(),
                a.func as u32,
                a.eps,
                a.eps_cell.get(),
            );
            let o = a.out.out().add(c);
            *o = if a.relax < 1.0 {
                *o + a.relax * (psi - *o)
            } else {
                psi
            };
        }
    }
}

/// Arguments of [`MassFlux`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct MassFluxArgs {
    /// Mesh.
    pub m: Mesh,
    /// Cell velocity.
    pub u: Raw3,
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// Cell pressure.
    pub p: Raw<f64>,
    /// Pressure gradient.
    pub grad_p: Raw3,
    /// Boundary pressure.
    pub pb: Raw<f64>,
    /// Boundary velocity.
    pub ub: Raw3,
    /// D per cell (zeros: the plain interpolated flux).
    pub d: Raw<f64>,
    /// Relaxation term per face.
    pub relax: Raw<f64>,
    /// 1 with the time-step term.
    pub transient: u32,
    /// D^s per cell.
    pub ds: Raw<f64>,
    /// Old-level part Y per face.
    pub time: Raw<f64>,
    /// α_u.
    pub alpha: f64,
    /// Velocity limiter per cell (read when `bounded` is 1).
    pub psi: Raw3,
    /// 1 for the bounded face velocity (D50).
    pub bounded: u32,
    /// Output per face.
    pub out: Raw<f64>,
}

/// Rhie–Chow face mass flux ([`k::face_mass_flux`]).
pub struct MassFlux;
impl Entry for MassFlux {
    type Args = MassFluxArgs;
    const NAME: &'static str = "mass_flux";
    #[inline(always)]
    unsafe fn run(f: usize, a: &MassFluxArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::FluxInputs {
                topo: a.m.topo(),
                geo: a.m.geo(),
                rho_face: a.m.rho_face.get(),
                u: a.u.v(),
                grad: a.grad_u.v(),
                p: a.p.get(),
                grad_p: a.grad_p.v(),
                pb: a.pb.get(),
                pset: a.m.pset.get(),
                ub: a.ub.v(),
                flux_kind: a.m.flux_kind.get(),
                mom_kind: a.m.mom_kind.get(),
                d: a.d.get(),
                relax: a.relax.get(),
                transient: a.transient != 0,
                ds: a.ds.get(),
                time: a.time.get(),
                alpha: a.alpha,
                psi: a.psi.v(),
                bounded: a.bounded != 0,
                jump_offset: a.m.jump_offset.get(),
                body: a.m.body.v(),
            };
            *a.out.out().add(f) = k::face_mass_flux(f, &inp);
        }
    }
}

/// Arguments of [`MomentumFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct MomentumFaceArgs {
    /// Mesh.
    pub m: Mesh,
    /// Convecting mass flux.
    pub mass_flux: Raw<f64>,
    /// Cell velocity.
    pub u: Raw3,
    /// Velocity of the central correction.
    pub uc: Raw3,
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// Limiter per component.
    pub psi: Raw3,
    /// Boundary velocity.
    pub ub: Raw3,
    /// 1 linear upwind, 0 upwind.
    pub blend: f64,
    /// Fraction of the bounded central scheme.
    pub central: f64,
    /// β of its bounding.
    pub beta: f64,
    /// 1 with the non-orthogonal correction.
    pub nonorth: f64,
    /// 1 with convection.
    pub convection: f64,
    /// 1 with the transposed viscous stress (turbulent flow).
    pub transposed: f64,
    /// 1 for its trace-free form (D50).
    pub trace_free: f64,
    /// Output c_P.
    pub cp: Raw<f64>,
    /// Output c_N.
    pub cn: Raw<f64>,
    /// Output b_f.
    pub bf: Raw3,
}

/// Momentum face coefficients ([`k::momentum_face`]).
pub struct MomentumFaces;
impl Entry for MomentumFaces {
    type Args = MomentumFaceArgs;
    const NAME: &'static str = "momentum_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &MomentumFaceArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::MomentumInputs {
                topo: a.m.topo(),
                geo: a.m.geo(),
                mass_flux: a.mass_flux.get(),
                mu: a.m.mu_face.get(),
                u: a.u.v(),
                uc: a.uc.v(),
                grad: a.grad_u.v(),
                psi: a.psi.v(),
                ub: a.ub.v(),
                mom_kind: a.m.mom_kind.get(),
                blend: a.blend,
                central: a.central,
                beta: a.beta,
                nonorth: a.nonorth,
                convection: a.convection,
                transposed: a.transposed,
                trace_free: a.trace_free,
            };
            let (cp, cn, b) = k::momentum_face(f, &inp);
            *a.cp.out().add(f) = cp;
            *a.cn.out().add(f) = cn;
            *a.bf.x.out().add(f) = b[0];
            *a.bf.y.out().add(f) = b[1];
            *a.bf.z.out().add(f) = b[2];
        }
    }
}

/// Arguments of [`MomentumRows`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct MomentumRowArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cp: Raw<f64>,
    /// See `cp`.
    pub cn: Raw<f64>,
    /// Deferred parts.
    pub bf: Raw3,
    /// Pressure gradient.
    pub grad_p: Raw3,
    /// Body force per unit volume.
    pub body_force: [f64; 3],
    /// Current velocity.
    pub u: Raw3,
    /// α_u.
    pub alpha: f64,
    /// θ.
    pub theta: f64,
    /// 1 in a time step.
    pub transient: u32,
    /// a₀, −a₁, −a₂ (the density is `m.rho`).
    pub time_a: [f64; 3],
    /// Δt.
    pub dt: f64,
    /// uⁿ.
    pub u_old: Raw3,
    /// uⁿ⁻¹.
    pub u_old2: Raw3,
    /// 1 − θ of the explicit operator.
    pub explicit_weight: f64,
    /// Explicit operator.
    pub explicit: Raw3,
    /// Brinkman sink s_P per cell (empty: none).
    pub sink: Raw<f64>,
    /// Output: matrix values.
    pub val: Raw<f64>,
    /// Output: right-hand sides.
    pub rhs: Raw3,
    /// Output: D.
    pub d: Raw<f64>,
    /// Output: D^s.
    pub ds: Raw<f64>,
    /// Output: D^C.
    pub dc: Raw<f64>,
    /// Output: relaxed diagonals.
    pub diag: Raw3,
}

/// Momentum rows ([`k::momentum_row`]).
pub struct MomentumRows;
impl Entry for MomentumRows {
    type Args = MomentumRowArgs;
    const NAME: &'static str = "momentum_rows";
    #[inline(always)]
    unsafe fn run(c: usize, a: &MomentumRowArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::MomentumRowInputs {
                topo: a.m.topo(),
                cp: a.cp.get(),
                cn: a.cn.get(),
                bf: a.bf.v(),
                grad_p: a.grad_p.v(),
                body_force: a.body_force,
                volume: a.m.volume.get(),
                u: a.u.v(),
                alpha: a.alpha,
                layout: a.m.layout(),
                mom_kind: a.m.mom_kind.get(),
                face_area: a.m.area.v(),
                theta: a.theta,
                transient: a.transient != 0,
                rho: a.m.rho.get(),
                rho_old: a.m.rho_old.get(),
                rho_old2: a.m.rho_old2.get(),
                time_a: a.time_a,
                dt: a.dt,
                u_old: a.u_old.v(),
                u_old2: a.u_old2.v(),
                explicit_weight: a.explicit_weight,
                explicit: a.explicit.v(),
                sink: a.sink.get(),
            };
            let out = k::MomentumRowOut {
                val: a.val.out(),
                rhs: [a.rhs.x.out(), a.rhs.y.out(), a.rhs.z.out()],
                d: a.d.out(),
                ds: a.ds.out(),
                dc: a.dc.out(),
                diag: [a.diag.x.out(), a.diag.y.out(), a.diag.z.out()],
            };
            k::momentum_row(c, &inp, &out);
        }
    }
}

/// Arguments of [`MomentumOperator`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct OperatorArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cp: Raw<f64>,
    /// See `cp`.
    pub cn: Raw<f64>,
    /// Deferred parts.
    pub bf: Raw3,
    /// Velocity.
    pub u: Raw3,
    /// Output.
    pub out: Raw3,
}

/// Spatial momentum operator ([`k::momentum_operator`]).
pub struct MomentumOperator;
impl Entry for MomentumOperator {
    type Args = OperatorArgs;
    const NAME: &'static str = "momentum_operator";
    #[inline(always)]
    unsafe fn run(c: usize, a: &OperatorArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = k::OperatorInputs {
                topo: a.m.topo(),
                cp: a.cp.get(),
                cn: a.cn.get(),
                bf: a.bf.v(),
                mom_kind: a.m.mom_kind.get(),
                face_area: a.m.area.v(),
                u: a.u.v(),
            };
            let r = k::momentum_operator(c, &inp);
            *a.out.x.out().add(c) = r[0];
            *a.out.y.out().add(c) = r[1];
            *a.out.z.out().add(c) = r[2];
        }
    }
}

/// Arguments of [`ResidualScale`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ResidualScaleArgs {
    /// Velocity.
    pub u: Raw3,
    /// Right-hand sides.
    pub rhs: Raw3,
    /// Relaxed diagonals.
    pub diag: Raw3,
    /// α_u.
    pub alpha: f64,
    /// Output a u (3 per cell, interleaved).
    pub v: Raw<f64>,
    /// Output b (3 per cell, interleaved).
    pub w: Raw<f64>,
}

/// Parts of the momentum residual scale ([`k::residual_scale`]).
pub struct ResidualScale;
impl Entry for ResidualScale {
    type Args = ResidualScaleArgs;
    const NAME: &'static str = "residual_scale";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ResidualScaleArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let (v, w) = k::residual_scale(c, &a.u.v(), &a.rhs.v(), &a.diag.v(), a.alpha);
            let mut i = 0;
            while i < 3 {
                *a.v.out().add(3 * c + i) = v[i];
                *a.w.out().add(3 * c + i) = w[i];
                i += 1;
            }
        }
    }
}

/// Arguments of the cell and face entries with one input and one output.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct UnaryArgs {
    /// Mesh.
    pub m: Mesh,
    /// Scalar parameter.
    pub s: f64,
    /// Input.
    pub x: Raw<f64>,
    /// Output.
    pub out: Raw<f64>,
}

/// Slot 0 of row `c` of `out` (matrix values) = `x[c]`.
pub struct SetDiagonal;
impl Entry for SetDiagonal {
    type Args = UnaryArgs;
    const NAME: &'static str = "set_diagonal";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(a.m.layout().pos(c, 0)) = ld(a.x.get(), c) }
    }
}

/// `out[c]` = 1 / slot 0 of row `c` of `x` (the Jacobi diagonal of a
/// matrix assembled on the device).
pub struct DiagonalInverse;
impl Entry for DiagonalInverse {
    type Args = UnaryArgs;
    const NAME: &'static str = "diagonal_inverse";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(c) = 1.0 / ld(a.x.get(), a.m.layout().pos(c, 0)) }
    }
}

/// Net outflow Σ ±ṁ of cell `c`.
pub struct NetOutflow;
impl Entry for NetOutflow {
    type Args = UnaryArgs;
    const NAME: &'static str = "net_outflow";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(c) = k::net_outflow(c, &a.m.topo(), a.x.get()) }
    }
}

/// Throughflow Σ |ṁ| of cell `c`.
pub struct Throughflow;
impl Entry for Throughflow {
    type Args = UnaryArgs;
    const NAME: &'static str = "throughflow";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe { *a.out.out().add(c) = k::throughflow(c, &a.m.topo(), a.x.get()) }
    }
}

/// Courant number of cell `c` (`s` = Δt / (2ρ)).
pub struct Courant;
impl Entry for Courant {
    type Args = UnaryArgs;
    const NAME: &'static str = "courant";
    #[inline(always)]
    unsafe fn run(c: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) = k::courant(c, &a.m.topo(), a.x.get(), a.m.volume.get(), a.s)
        }
    }
}

/// Pressure-correction coefficient of face `f` from D^C (`s` unused; the
/// density is `m.rho_face`).
pub struct PressureFaces;
impl Entry for PressureFaces {
    type Args = UnaryArgs;
    const NAME: &'static str = "pressure_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(f) = k::pressure_face(
                f,
                &a.m.topo(),
                &a.m.geo(),
                a.m.flux_kind.get(),
                a.m.rho_face.get(),
                a.x.get(),
                a.m.closed.get(),
            )
        }
    }
}

/// Arguments of [`CloseBackflow`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CloseArgs {
    /// Mesh (its `closed` array is written).
    pub m: Mesh,
    /// Predicted face mass flux, updated in place on the boundary faces.
    pub mass_flux: Raw<f64>,
}

/// Backflow suppression of boundary face `fb` ([`k::close_backflow`]):
/// writes the flux of the face (the output is the input, at the same
/// index) and its closed state.
pub struct CloseBackflow;
impl Entry for CloseBackflow {
    type Args = CloseArgs;
    const NAME: &'static str = "close_backflow";
    #[inline(always)]
    unsafe fn run(fb: usize, a: &CloseArgs) {
        // SAFETY: guaranteed by the caller; index fb only.
        unsafe {
            let f = fb + a.m.n_interior as usize;
            let kind = ld(a.m.flux_kind.get(), fb);
            if kind == k::FLUX_PRESSURE_NO_BACKFLOW {
                let (v, c) = k::close_backflow(kind, a.mass_flux.read(f));
                *a.mass_flux.out().add(f) = v;
                *a.m.closed.out().add(fb) = c;
            } else {
                *a.m.closed.out().add(fb) = 0.0;
            }
        }
    }
}

/// Arguments of [`PressureRows`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct PressureRowArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cf: Raw<f64>,
    /// Mass flux.
    pub mass_flux: Raw<f64>,
    /// Output: matrix values.
    pub val: Raw<f64>,
    /// Output: right-hand side.
    pub rhs: Raw<f64>,
}

/// Pressure-correction rows ([`k::pressure_row`]).
pub struct PressureRows;
impl Entry for PressureRows {
    type Args = PressureRowArgs;
    const NAME: &'static str = "pressure_rows";
    #[inline(always)]
    unsafe fn run(c: usize, a: &PressureRowArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            k::pressure_row(
                c,
                &a.m.topo(),
                a.cf.get(),
                a.mass_flux.get(),
                a.m.mass_source.get(),
                &a.m.layout(),
                a.val.out(),
                a.rhs.out(),
            )
        }
    }
}

/// Arguments of [`CorrectFluxes`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CorrectFluxArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cf: Raw<f64>,
    /// Predicted mass flux.
    pub mass_flux: Raw<f64>,
    /// Pressure correction.
    pub pc: Raw<f64>,
    /// Output.
    pub out: Raw<f64>,
}

/// Corrected face flux ([`k::correct_flux`]).
pub struct CorrectFluxes;
impl Entry for CorrectFluxes {
    type Args = CorrectFluxArgs;
    const NAME: &'static str = "correct_fluxes";
    #[inline(always)]
    unsafe fn run(f: usize, a: &CorrectFluxArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(f) =
                k::correct_flux(f, &a.m.topo(), a.cf.get(), a.mass_flux.get(), a.pc.get())
        }
    }
}

/// Arguments of [`CorrectVelocity`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CorrectVelocityArgs {
    /// D^C.
    pub dc: Raw<f64>,
    /// Gradient of the pressure correction.
    pub grad: Raw3,
    /// Velocity, updated in place.
    pub u: Raw3,
}

/// u_i −= D^C ∇p'_i of cell `c`.
pub struct CorrectVelocity;
impl Entry for CorrectVelocity {
    type Args = CorrectVelocityArgs;
    const NAME: &'static str = "correct_velocity";
    #[inline(always)]
    unsafe fn run(c: usize, a: &CorrectVelocityArgs) {
        // SAFETY: guaranteed by the caller; index c only.
        unsafe {
            let d = a.dc.read(c);
            *a.u.x.out().add(c) = a.u.x.read(c) - d * a.grad.x.read(c);
            *a.u.y.out().add(c) = a.u.y.read(c) - d * a.grad.y.read(c);
            *a.u.z.out().add(c) = a.u.z.read(c) - d * a.grad.z.read(c);
        }
    }
}

/// Arguments of [`NeighbourCorrection`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct NeighbourArgs {
    /// Mesh.
    pub m: Mesh,
    /// Matrix values.
    pub val: Raw<f64>,
    /// Cell values.
    pub x: Raw<f64>,
    /// Diagonal.
    pub diag: Raw<f64>,
    /// Output.
    pub out: Raw<f64>,
    /// 1 after a SIMPLEC-type corrector (subtract Σ_N val_PN x_P), 0 after
    /// Issa's.
    pub own: f64,
}

/// −Σ_N val_PN (x_N − own x_P) / diag of cell `c`: the PISO neighbour
/// correction H(x)/a_P after Issa's corrector (`own` = 0), or
/// Σ_N |a_N| (x_N − x_P)/a_P after a corrector with D^C = V/(a_P − Σ|a_N|)
/// (`own` = 1).
pub struct NeighbourCorrection;
impl Entry for NeighbourCorrection {
    type Args = NeighbourArgs;
    const NAME: &'static str = "neighbour_correction";
    #[inline(always)]
    unsafe fn run(c: usize, a: &NeighbourArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) = -k::neighbour_difference(
                c,
                &a.m.topo(),
                &a.m.layout(),
                a.val.get(),
                a.x.get(),
                a.own,
            ) / a.diag.read(c)
        }
    }
}

/// Slot `k` of Aᵀ in the storage of A: `out[k] = x[transpose[k]]`, for
/// every stored slot (D59).
pub struct Transpose;
impl Entry for Transpose {
    type Args = UnaryArgs;
    const NAME: &'static str = "transpose";
    #[inline(always)]
    unsafe fn run(k: usize, a: &UnaryArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(k) = adjoint::transposed_value(k, a.m.transpose.get(), a.x.get())
        }
    }
}

/// Arguments of [`MatrixProduct`] and [`TransposedProduct`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ProductArgs {
    /// Mesh.
    pub m: Mesh,
    /// Matrix values.
    pub val: Raw<f64>,
    /// Cell values.
    pub x: Raw<f64>,
    /// Output.
    pub out: Raw<f64>,
}

/// Row `c` of A x.
pub struct MatrixProduct;
impl Entry for MatrixProduct {
    type Args = ProductArgs;
    const NAME: &'static str = "matrix_product";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ProductArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) =
                adjoint::row_product(c, &a.m.topo(), &a.m.layout(), a.val.get(), a.x.get())
        }
    }
}

/// Row `c` of Aᵀ x, without storing Aᵀ (D59).
pub struct TransposedProduct;
impl Entry for TransposedProduct {
    type Args = ProductArgs;
    const NAME: &'static str = "transposed_product";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ProductArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) = adjoint::transposed_row_product(
                c,
                &a.m.topo(),
                &a.m.layout(),
                a.m.transpose.get(),
                a.val.get(),
                a.x.get(),
            )
        }
    }
}

/// Arguments of [`DeltaFlux`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct DeltaFluxArgs {
    /// Mesh.
    pub m: Mesh,
    /// Velocity increment.
    pub delta: Raw3,
    /// Output per face.
    pub out: Raw<f64>,
}

/// Mass flux of a velocity increment ([`k::delta_flux`]).
pub struct DeltaFlux;
impl Entry for DeltaFlux {
    type Args = DeltaFluxArgs;
    const NAME: &'static str = "delta_flux";
    #[inline(always)]
    unsafe fn run(f: usize, a: &DeltaFluxArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(f) = k::delta_flux(
                f,
                &a.m.topo(),
                &a.m.geo(),
                a.m.flux_kind.get(),
                a.m.rho_face.get(),
                &a.delta.v(),
                a.m.closed.get(),
            )
        }
    }
}

/// Arguments of [`ScalarFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ScalarFaceArgs {
    /// Mesh.
    pub m: Mesh,
    /// Boundary coefficient a of the scalar per boundary face.
    pub bc_a: Raw<f64>,
    /// Boundary constant b.
    pub bc_b: Raw<f64>,
    /// Jump per periodic face (empty: plain periodic).
    pub jump: Raw<f64>,
    /// Convecting face mass flux.
    pub mass_flux: Raw<f64>,
    /// ρΓ per face.
    pub gamma: Raw<f64>,
    /// Cell values.
    pub phi: Raw<f64>,
    /// Cell gradient.
    pub grad: Raw3,
    /// Limiter per cell.
    pub psi: Raw<f64>,
    /// 1 linear upwind, 0 upwind.
    pub blend: f64,
    /// Fraction of the bounded central scheme.
    pub central: f64,
    /// β of its bounding.
    pub beta: f64,
    /// 1 with the non-orthogonal correction.
    pub nonorth: f64,
    /// Limit of the non-orthogonal correction (< 0: none).
    pub nonorth_limit: f64,
    /// 1: linear-upwind face values between the two cell values.
    pub bound: f64,
    /// Output c_P.
    pub cp: Raw<f64>,
    /// Output c_N.
    pub cn: Raw<f64>,
    /// Output b_f.
    pub bf: Raw<f64>,
}

/// Face coefficients of a scalar ([`fv::face_coefficients`]).
pub struct ScalarFaces;
impl Entry for ScalarFaces {
    type Args = ScalarFaceArgs;
    const NAME: &'static str = "scalar_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &ScalarFaceArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = fv::FaceInputs {
                topo: a.m.topo(),
                bc: Boundary {
                    a: a.bc_a.get(),
                    b: a.bc_b.get(),
                    jump: a.jump.get(),
                },
                mass_flux: a.mass_flux.get(),
                gamma: a.gamma.get(),
                dgeo: a.m.dgeo.get(),
                k: a.m.k.v(),
                weight: a.m.weight.get(),
                face_centroid: a.m.face_centroid.v(),
                cell_centroid: a.m.cell_centroid.v(),
                shift: a.m.shift.v(),
                phi: a.phi.get(),
                grad: a.grad.v(),
                psi: a.psi.get(),
                blend: a.blend,
                central: a.central,
                beta: a.beta,
                nonorth: a.nonorth,
                nonorth_limit: a.nonorth_limit,
                bound: a.bound,
            };
            let (cp, cn, b) = fv::face_coefficients(f, &inp);
            *a.cp.out().add(f) = cp;
            *a.cn.out().add(f) = cn;
            *a.bf.out().add(f) = b;
        }
    }
}

/// Arguments of [`ScalarRows`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ScalarRowArgs {
    /// Mesh.
    pub m: Mesh,
    /// Face coefficients.
    pub cp: Raw<f64>,
    /// See `cp`.
    pub cn: Raw<f64>,
    /// Deferred parts.
    pub bf: Raw<f64>,
    /// Convecting face mass flux (bounded form).
    pub mass_flux: Raw<f64>,
    /// Current cell values.
    pub phi: Raw<f64>,
    /// S_u per cell (empty: none).
    pub su: Raw<f64>,
    /// S_p per cell (empty: none).
    pub sp: Raw<f64>,
    /// α.
    pub alpha: f64,
    /// 1 in a time step.
    pub transient: u64,
    /// 1/Δt.
    pub inv_dt: f64,
    /// a₀, a₁, a₂.
    pub a: [f64; 3],
    /// φⁿ.
    pub old: Raw<f64>,
    /// φⁿ⁻¹.
    pub old2: Raw<f64>,
    /// Lower bound of the bounded BDF2.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
    /// 1 in the cells with a fixed value (empty: none).
    pub fix_mask: Raw<f64>,
    /// The fixed values.
    pub fix_val: Raw<f64>,
    /// Pseudo-time coefficient ρ V/Δτ per cell (empty: none).
    pub pseudo: Raw<f64>,
    /// The heat capacity ρ c_p per cell of the energy equation (stage
    /// 7a1), in place of the density of the mesh in the time and source
    /// terms (empty: the density).
    pub capacity: Raw<f64>,
    /// 1: a positive variable (the deferred corrections split by sign).
    pub positive: u32,
    /// Output: matrix values.
    pub val: Raw<f64>,
    /// Output: right-hand side.
    pub rhs: Raw<f64>,
    /// Output: relaxed diagonal.
    pub diag: Raw<f64>,
    /// Output: 1 where the bounded BDF2 fell back to first order.
    pub order: Raw<f64>,
}

/// Row of a scalar ([`crate::scalar::scalar_row`]).
pub struct ScalarRows;
impl Entry for ScalarRows {
    type Args = ScalarRowArgs;
    const NAME: &'static str = "scalar_rows";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ScalarRowArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let inp = crate::scalar::ScalarRowInputs {
                topo: a.m.topo(),
                cp: a.cp.get(),
                cn: a.cn.get(),
                bf: a.bf.get(),
                mass_flux: a.mass_flux.get(),
                volume: a.m.volume.get(),
                phi: a.phi.get(),
                rho: if a.capacity.len > 0 {
                    a.capacity.get()
                } else {
                    a.m.rho.get()
                },
                su: a.su.get(),
                sp: a.sp.get(),
                alpha: a.alpha,
                layout: a.m.layout(),
                transient: a.transient != 0,
                inv_dt: a.inv_dt,
                a: a.a,
                old: a.old.get(),
                old2: a.old2.get(),
                lo: a.lo,
                hi: a.hi,
                fix_mask: a.fix_mask.get(),
                fix_val: a.fix_val.get(),
                positive: a.positive != 0,
                pseudo: a.pseudo.get(),
            };
            let r = crate::scalar::scalar_row(c, &inp, a.val.out());
            *a.rhs.out().add(c) = r.rhs;
            *a.diag.out().add(c) = r.diag;
            *a.order.out().add(c) = r.first_order;
        }
    }
}

/// Arguments of [`ScalarScale`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ScalarScaleArgs {
    /// Cell values.
    pub phi: Raw<f64>,
    /// Right-hand side.
    pub rhs: Raw<f64>,
    /// Relaxed diagonal.
    pub diag: Raw<f64>,
    /// α.
    pub alpha: f64,
    /// Pseudo-time coefficient ρ V/Δτ per cell (empty: none).
    pub pseudo: Raw<f64>,
    /// Output a φ.
    pub v: Raw<f64>,
    /// Output b − (a/α − a) φ.
    pub w: Raw<f64>,
}

/// Parts of the residual scale of a scalar
/// ([`crate::scalar::residual_parts`]).
pub struct ScalarScale;
impl Entry for ScalarScale {
    type Args = ScalarScaleArgs;
    const NAME: &'static str = "scalar_scale";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ScalarScaleArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let ps = a.pseudo.get();
            let (v, w) = crate::scalar::residual_parts(
                ld(a.phi.get(), c),
                ld(a.rhs.get(), c),
                ld(a.diag.get(), c),
                a.alpha,
                if ps.is_empty() { 0.0 } else { ld(ps, c) },
            );
            *a.v.out().add(c) = v;
            *a.w.out().add(c) = w;
        }
    }
}

/// Arguments of [`SstCells`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct SstCellArgs {
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// k.
    pub k: Raw<f64>,
    /// ω.
    pub w: Raw<f64>,
    /// ∇k.
    pub grad_k: Raw3,
    /// ∇ω.
    pub grad_w: Raw3,
    /// Wall distance.
    pub d: Raw<f64>,
    /// ν.
    pub nu: f64,
    /// ν per cell (D194: several fluids), or [`Raw::EMPTY`]: `nu` in
    /// every cell.
    pub nu_cell: Raw<f64>,
    /// Constants.
    pub c: crate::turb::sst::SstConstants,
    /// Options.
    pub o: crate::turb::sst::SstOptions,
    /// Additional explicit sources of k and ω per cell (MMS, D104);
    /// empty: none.
    pub add_k: Raw<f64>,
    /// See `add_k`.
    pub add_w: Raw<f64>,
    /// Weight m of the wall production per cell (the automatic wall
    /// treatment, 0 off the near-wall cells); empty: none.
    pub wall_mask: Raw<f64>,
    /// The weighted wall production per cell, added to (1 − m) P̃_k.
    pub wall_pk: Raw<f64>,
    /// f_r1 of the curvature correction per cell; empty: 1.
    pub fr1: Raw<f64>,
    /// Output ν_t.
    pub nu_t: Raw<f64>,
    /// Output F₁.
    pub f1: Raw<f64>,
    /// Output P̃_k.
    pub pk: Raw<f64>,
    /// Output S_u of k (empty: ν_t only).
    pub su_k: Raw<f64>,
    /// Output S_p of k.
    pub sp_k: Raw<f64>,
    /// Output S_u of ω.
    pub su_w: Raw<f64>,
    /// Output S_p of ω.
    pub sp_w: Raw<f64>,
    /// Output σ_k ν_t.
    pub dk: Raw<f64>,
    /// Output σ_ω ν_t.
    pub dw: Raw<f64>,
    /// Output 1 where the cap of ν_t/ν acted.
    pub capped: Raw<f64>,
}

/// The k-ω SST model per cell ([`crate::turb::sst::sst_cell`]).
pub struct SstCells;
impl Entry for SstCells {
    type Args = SstCellArgs;
    const NAME: &'static str = "sst_cells";
    #[inline(always)]
    unsafe fn run(c: usize, a: &SstCellArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = a.grad_u.v();
            let gu = [g.x.at(c), g.y.at(c), g.z.at(c)];
            let inp = crate::turb::sst::SstCellInputs {
                k: ld(a.k.get(), c),
                w: ld(a.w.get(), c),
                grad_k: a.grad_k.v().at(c),
                grad_w: a.grad_w.v().at(c),
                grad_u: gu,
                d: ld(a.d.get(), c),
                nu: if a.nu_cell.len > 0 {
                    ld(a.nu_cell.get(), c)
                } else {
                    a.nu
                },
                fr1: if a.fr1.len > 0 {
                    ld(a.fr1.get(), c)
                } else {
                    1.0
                },
            };
            let r = crate::turb::sst::sst_cell(&inp, &a.c, &a.o);
            *a.nu_t.out().add(c) = r.nu_t;
            *a.f1.out().add(c) = r.f1;
            *a.pk.out().add(c) = r.pk;
            *a.dk.out().add(c) = r.dk;
            *a.dw.out().add(c) = r.dw;
            *a.capped.out().add(c) = r.capped;
            if a.su_k.len > 0 {
                let (ak, aw) = if a.add_k.len > 0 {
                    (ld(a.add_k.get(), c), ld(a.add_w.get(), c))
                } else {
                    (0.0, 0.0)
                };
                // an added source splits like the others: its positive
                // part explicit, its negative part implicit
                let (k, w) = (inp.k, inp.w);
                // near-wall cells of the automatic treatment: the
                // production (1 − m) P̃_k + m u_τ³/(κ y_P) (the weighted
                // wall part comes from [`WallCells`])
                let su_k = if a.wall_mask.len > 0 {
                    let m = ld(a.wall_mask.get(), c);
                    r.su_k - m * r.pk + ld(a.wall_pk.get(), c)
                } else {
                    r.su_k
                };
                *a.su_k.out().add(c) = su_k + ak.max(0.0);
                *a.sp_k.out().add(c) = r.sp_k + ak.min(0.0) / k.max(1e-300);
                *a.su_w.out().add(c) = r.su_w + aw.max(0.0);
                *a.sp_w.out().add(c) = r.sp_w + aw.min(0.0) / w.max(1e-300);
            }
        }
    }
}

/// Arguments of [`StrainCells`], [`ConvectiveDerivative`] and
/// [`CurvatureCells`] (the curvature correction of SST, D170).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CurvatureArgs {
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// [`StrainCells`]: the component of S_ij (0 … 5: xx, yy, zz, xy, yz,
    /// zx).
    pub comp: u32,
    /// [`ConvectiveDerivative`]: the convecting cell velocity.
    pub u: Raw3,
    /// The frame velocity per cell, subtracted from `u` (empty: 0).
    pub v_frame: Raw3,
    /// The gradient of the component.
    pub grad: Raw3,
    /// The component now and at the last time level (`s_old` empty: no
    /// time part).
    pub s: Raw<f64>,
    /// See `s`.
    pub s_old: Raw<f64>,
    /// 1/Δt of the time part.
    pub inv_dt: f64,
    /// [`CurvatureCells`]: DS_ij/Dt (order of `comp`).
    pub dsdt: [Raw<f64>; 6],
    /// ω.
    pub w: Raw<f64>,
    /// Ω^rot per cell (empty: `rot`).
    pub omega_rot: Raw3,
    /// Ω^rot of every cell when `omega_rot` is empty.
    pub rot: [f64; 3],
    /// The rotation of the frame of a relative velocity (spin).
    pub spin: [f64; 3],
    /// Constants (c_r1 … c_r3).
    pub c: crate::turb::sst::SstConstants,
    /// Width of the smooth limits (0: hard).
    pub smooth: f64,
    /// Output: S_comp, the derivative, or f_r1.
    pub out: Raw<f64>,
}

/// out = S_comp of the velocity gradient.
pub struct StrainCells;
impl Entry for StrainCells {
    type Args = CurvatureArgs;
    const NAME: &'static str = "strain_cells";
    #[inline(always)]
    unsafe fn run(c: usize, a: &CurvatureArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = a.grad_u.v();
            let gu = [g.x.at(c), g.y.at(c), g.z.at(c)];
            let (i, j) = match a.comp {
                0 => (0, 0),
                1 => (1, 1),
                2 => (2, 2),
                3 => (0, 1),
                4 => (1, 2),
                _ => (2, 0),
            };
            *a.out.out().add(c) = 0.5 * (gu[i][j] + gu[j][i]);
        }
    }
}

/// out = (u − v_frame)·grad + (s − s_old)/Δt: the material derivative of
/// one component of S_ij (BDF1 time part).
pub struct ConvectiveDerivative;
impl Entry for ConvectiveDerivative {
    type Args = CurvatureArgs;
    const NAME: &'static str = "convective_derivative";
    #[inline(always)]
    unsafe fn run(c: usize, a: &CurvatureArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let mut u = [a.u.x.read(c), a.u.y.read(c), a.u.z.read(c)];
            if a.v_frame.x.len > 0 {
                u[0] -= a.v_frame.x.read(c);
                u[1] -= a.v_frame.y.read(c);
                u[2] -= a.v_frame.z.read(c);
            }
            let mut d = u[0] * a.grad.x.read(c) + u[1] * a.grad.y.read(c) + u[2] * a.grad.z.read(c);
            if a.s_old.len > 0 {
                d += (a.s.read(c) - a.s_old.read(c)) * a.inv_dt;
            }
            *a.out.out().add(c) = d;
        }
    }
}

/// out = f_r1 ([`crate::turb::sst::rotation_function`]).
pub struct CurvatureCells;
impl Entry for CurvatureCells {
    type Args = CurvatureArgs;
    const NAME: &'static str = "curvature_cells";
    #[inline(always)]
    unsafe fn run(c: usize, a: &CurvatureArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = a.grad_u.v();
            let gu = [g.x.at(c), g.y.at(c), g.z.at(c)];
            let ds = [
                a.dsdt[0].read(c),
                a.dsdt[1].read(c),
                a.dsdt[2].read(c),
                a.dsdt[3].read(c),
                a.dsdt[4].read(c),
                a.dsdt[5].read(c),
            ];
            let rot = if a.omega_rot.x.len > 0 {
                [
                    a.omega_rot.x.read(c),
                    a.omega_rot.y.read(c),
                    a.omega_rot.z.read(c),
                ]
            } else {
                a.rot
            };
            *a.out.out().add(c) = crate::turb::sst::rotation_function(
                &gu,
                &ds,
                rot,
                a.spin,
                a.w.read(c),
                &a.c,
                a.smooth,
            );
        }
    }
}

/// Arguments of [`QcrCells`] and [`ExplicitStressFaces`] (the quadratic
/// constitutive relation of SST, D170).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct QcrArgs {
    /// Mesh ([`ExplicitStressFaces`]).
    pub m: Mesh,
    /// Velocity gradient ([`QcrCells`]).
    pub grad_u: RawGrad,
    /// ν_t per cell.
    pub nu_t: Raw<f64>,
    /// Constants (c_cr1, c_cr2).
    pub c: crate::turb::sst::SstConstants,
    /// [`crate::turb::sst::QCR_2000`] or [`crate::turb::sst::QCR_2013`].
    pub mode: u32,
    /// The cell stress (xx, yy, zz, xy, yz, zx): output of [`QcrCells`],
    /// input of [`ExplicitStressFaces`].
    pub tau: [Raw<f64>; 6],
    /// [`ExplicitStressFaces`]: b_f, updated in place.
    pub bf: Raw3,
    /// [`ExplicitStressFaces`]: the cell gradients of the six stress
    /// components for the skewness correction of the face stress (D236);
    /// empty: the plain interpolation.
    pub gtau: [Raw3; 6],
}

/// tau = the non-linear stress of a cell ([`crate::turb::sst::qcr_stress`]).
pub struct QcrCells;
impl Entry for QcrCells {
    type Args = QcrArgs;
    const NAME: &'static str = "qcr_cells";
    #[inline(always)]
    unsafe fn run(c: usize, a: &QcrArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = a.grad_u.v();
            let gu = [g.x.at(c), g.y.at(c), g.z.at(c)];
            let t = crate::turb::sst::qcr_stress(
                &gu,
                ld(a.m.rho.get(), c) * a.nu_t.read(c),
                &a.c,
                a.mode,
            );
            let mut i = 0;
            while i < 6 {
                *a.tau[i].out().add(c) = t[i];
                i += 1;
            }
        }
    }
}

/// b_f −= τ_f·S with τ_f = w τ_P + (1 − w) τ_N on the interior faces and
/// τ_P on the boundary faces: the divergence of an explicit cell stress
/// in the momentum equations (b_f is subtracted from the owner row).
///
/// With the stress gradients `gtau` (D236) the face stress is corrected
/// for skewness as the face velocity of the mass flux:
/// τ_f = w τ_P + (1 − w) τ_N + ∇τ̄_f·(x_f − x_w), ∇τ̄_f = w ∇τ_P +
/// (1 − w) ∇τ_N, x_w = w x_P + (1 − w) x_N, and on a boundary face
/// τ_f = τ_P + ∇τ_P·(x_f − x_P); the face stress is then second order
/// where the cell stress is (exact for a linear stress).
pub struct ExplicitStressFaces;
impl Entry for ExplicitStressFaces {
    type Args = QcrArgs;
    const NAME: &'static str = "explicit_stress_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &QcrArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let o = a.m.owner.read(f) as usize;
            let mut t = [0.0; 6];
            let mut i = 0;
            while i < 6 {
                t[i] = a.tau[i].read(o);
                i += 1;
            }
            let skew = a.gtau[0].x.len > 0;
            let xf = [
                a.m.face_centroid.x.read(f),
                a.m.face_centroid.y.read(f),
                a.m.face_centroid.z.read(f),
            ];
            let xp = [
                a.m.cell_centroid.x.read(o),
                a.m.cell_centroid.y.read(o),
                a.m.cell_centroid.z.read(o),
            ];
            if f < a.m.n_interior as usize {
                let n = a.m.neighbour.read(f) as usize;
                let w = a.m.weight.read(f);
                let mut i = 0;
                while i < 6 {
                    t[i] = w * t[i] + (1.0 - w) * a.tau[i].read(n);
                    i += 1;
                }
                if skew {
                    let g = a.m.geo();
                    let xn = crate::fv::neighbour_centroid(
                        f,
                        n,
                        &a.m.topo(),
                        &g.cell_centroid,
                        &g.shift,
                    );
                    let r = [
                        xf[0] - (w * xp[0] + (1.0 - w) * xn[0]),
                        xf[1] - (w * xp[1] + (1.0 - w) * xn[1]),
                        xf[2] - (w * xp[2] + (1.0 - w) * xn[2]),
                    ];
                    let mut i = 0;
                    while i < 6 {
                        let gt = a.gtau[i];
                        let gx = w * gt.x.read(o) + (1.0 - w) * gt.x.read(n);
                        let gy = w * gt.y.read(o) + (1.0 - w) * gt.y.read(n);
                        let gz = w * gt.z.read(o) + (1.0 - w) * gt.z.read(n);
                        t[i] += gx * r[0] + gy * r[1] + gz * r[2];
                        i += 1;
                    }
                }
            } else if skew {
                let r = [xf[0] - xp[0], xf[1] - xp[1], xf[2] - xp[2]];
                let mut i = 0;
                while i < 6 {
                    let gt = a.gtau[i];
                    t[i] += gt.x.read(o) * r[0] + gt.y.read(o) * r[1] + gt.z.read(o) * r[2];
                    i += 1;
                }
            }
            let s = [a.m.area.x.read(f), a.m.area.y.read(f), a.m.area.z.read(f)];
            // xx, yy, zz, xy, yz, zx
            let fx = t[0] * s[0] + t[3] * s[1] + t[5] * s[2];
            let fy = t[3] * s[0] + t[1] * s[1] + t[4] * s[2];
            let fz = t[5] * s[0] + t[4] * s[1] + t[2] * s[2];
            *a.bf.x.out().add(f) -= fx;
            *a.bf.y.out().add(f) -= fy;
            *a.bf.z.out().add(f) -= fz;
        }
    }
}

/// Arguments of [`QuadraticGradient`] (D236).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct QuadGradArgs {
    /// The treated cells.
    pub cells: Raw<u32>,
    /// Start of the stencil of each treated cell (one more entry).
    pub ptr: Raw<u32>,
    /// Stencil points: a cell, or `n_cells` + a boundary face index.
    pub idx: Raw<u32>,
    /// Gradient weights per stencil point.
    pub w: Raw3,
    /// Number of cells.
    pub n_cells: u64,
    /// Components of the gradient the fit gives (2: x and y; the others
    /// keep the face-neighbour gradient).
    pub dims: u64,
    /// Cell velocity.
    pub u: Raw3,
    /// Boundary face velocity.
    pub ub: Raw3,
    /// The velocity gradient, overwritten at the treated cells.
    pub out: RawGrad,
}

/// The quadratic least-squares velocity gradient of one treated cell
/// (`foehncraft_fv::exactness::QuadraticGradient`): g = Σ_k w_k (u_k − u_c)
/// per component, written over the face-neighbour gradient.
pub struct QuadraticGradient;
impl Entry for QuadraticGradient {
    type Args = QuadGradArgs;
    const NAME: &'static str = "quadratic_gradient";
    #[inline(always)]
    unsafe fn run(t: usize, a: &QuadGradArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let c = a.cells.read(t) as usize;
            let lo = a.ptr.read(t) as usize;
            let hi = a.ptr.read(t + 1) as usize;
            let nc = a.n_cells as usize;
            let uc = [a.u.x.read(c), a.u.y.read(c), a.u.z.read(c)];
            let mut g = [[0.0f64; 3]; 3];
            let mut k = lo;
            while k < hi {
                let i = a.idx.read(k) as usize;
                let v = if i < nc {
                    [a.u.x.read(i), a.u.y.read(i), a.u.z.read(i)]
                } else {
                    let b = i - nc;
                    [a.ub.x.read(b), a.ub.y.read(b), a.ub.z.read(b)]
                };
                let w = [a.w.x.read(k), a.w.y.read(k), a.w.z.read(k)];
                let mut comp = 0;
                while comp < 3 {
                    let d = v[comp] - uc[comp];
                    g[comp][0] += w[0] * d;
                    g[comp][1] += w[1] * d;
                    g[comp][2] += w[2] * d;
                    comp += 1;
                }
                k += 1;
            }
            let mut comp = 0;
            while comp < 3 {
                let o = a.out.g[comp];
                *o.x.out().add(c) = g[comp][0];
                if a.dims > 1 {
                    *o.y.out().add(c) = g[comp][1];
                }
                if a.dims > 2 {
                    *o.z.out().add(c) = g[comp][2];
                }
                comp += 1;
            }
        }
    }
}

/// Arguments of [`KepsCells`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct KepsCellArgs {
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// k.
    pub k: Raw<f64>,
    /// ε.
    pub e: Raw<f64>,
    /// ν.
    pub nu: f64,
    /// ν per cell (D194: several fluids), or [`Raw::EMPTY`]: `nu` in
    /// every cell.
    pub nu_cell: Raw<f64>,
    /// Constants.
    pub c: crate::turb::keps::KepsConstants,
    /// Options.
    pub o: crate::turb::keps::KepsOptions,
    /// Additional explicit sources of k and ε per cell (MMS, D104);
    /// empty: none.
    pub add_k: Raw<f64>,
    /// See `add_k`.
    pub add_e: Raw<f64>,
    /// Weight m_p of the wall production per cell (0 off the near-wall
    /// cells); empty: none.
    pub wall_mask: Raw<f64>,
    /// The weighted wall production per cell, added to (1 − m_p) P_k.
    pub wall_pk: Raw<f64>,
    /// Output ν_t.
    pub nu_t: Raw<f64>,
    /// Output C_μ.
    pub c_mu: Raw<f64>,
    /// Output P_k.
    pub pk: Raw<f64>,
    /// Output S_u of k (empty: ν_t only).
    pub su_k: Raw<f64>,
    /// Output S_p of k.
    pub sp_k: Raw<f64>,
    /// Output S_u of ε.
    pub su_e: Raw<f64>,
    /// Output S_p of ε.
    pub sp_e: Raw<f64>,
    /// Output ν_t/σ_k.
    pub dk: Raw<f64>,
    /// Output ν_t/σ_ε.
    pub de: Raw<f64>,
    /// Output 1 where the cap of ν_t/ν acted.
    pub capped: Raw<f64>,
}

/// The standard or realizable k-ε model per cell
/// ([`crate::turb::keps::keps_cell`]).
pub struct KepsCells;
impl Entry for KepsCells {
    type Args = KepsCellArgs;
    const NAME: &'static str = "keps_cells";
    #[inline(always)]
    unsafe fn run(c: usize, a: &KepsCellArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = a.grad_u.v();
            let inp = crate::turb::keps::KepsCellInputs {
                k: ld(a.k.get(), c),
                e: ld(a.e.get(), c),
                grad_u: [g.x.at(c), g.y.at(c), g.z.at(c)],
                nu: if a.nu_cell.len > 0 {
                    ld(a.nu_cell.get(), c)
                } else {
                    a.nu
                },
            };
            let r = crate::turb::keps::keps_cell(&inp, &a.c, &a.o);
            *a.nu_t.out().add(c) = r.nu_t;
            *a.c_mu.out().add(c) = r.c_mu;
            *a.pk.out().add(c) = r.pk;
            *a.dk.out().add(c) = r.dk;
            *a.de.out().add(c) = r.de;
            *a.capped.out().add(c) = r.capped;
            if a.su_k.len > 0 {
                let (ak, ae) = if a.add_k.len > 0 {
                    (ld(a.add_k.get(), c), ld(a.add_e.get(), c))
                } else {
                    (0.0, 0.0)
                };
                let (k, e) = (inp.k, inp.e);
                // near-wall cells: (1 − m_p) P_k + the wall production
                let su_k = if a.wall_mask.len > 0 {
                    let m = ld(a.wall_mask.get(), c);
                    r.su_k - m * r.pk + ld(a.wall_pk.get(), c)
                } else {
                    r.su_k
                };
                *a.su_k.out().add(c) = su_k + ak.max(0.0);
                *a.sp_k.out().add(c) = r.sp_k + ak.min(0.0) / k.max(1e-300);
                *a.su_e.out().add(c) = r.su_e + ae.max(0.0);
                *a.sp_e.out().add(c) = r.sp_e + ae.min(0.0) / e.max(1e-300);
            }
        }
    }
}

/// Arguments of [`ZetaFCells`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ZetaFCellArgs {
    /// Velocity gradient.
    pub grad_u: RawGrad,
    /// k.
    pub k: Raw<f64>,
    /// ε.
    pub e: Raw<f64>,
    /// ζ.
    pub zeta: Raw<f64>,
    /// f.
    pub f: Raw<f64>,
    /// ν.
    pub nu: f64,
    /// ν per cell (D194: several fluids), or [`Raw::EMPTY`]: `nu` in
    /// every cell.
    pub nu_cell: Raw<f64>,
    /// Constants.
    pub c: crate::turb::zeta_f::ZetaFConstants,
    /// Options.
    pub o: crate::turb::zeta_f::ZetaFOptions,
    /// Additional explicit sources of k, ε and ζ per cell (MMS, D104);
    /// empty: none.
    pub add_k: Raw<f64>,
    /// See `add_k`.
    pub add_e: Raw<f64>,
    /// See `add_k`.
    pub add_z: Raw<f64>,
    /// Added to the right-hand side of the f (α) equation (with `add_k`).
    pub add_f: Raw<f64>,
    /// Weight m_p of the wall production per cell (0 off the near-wall
    /// cells); empty: none.
    pub wall_mask: Raw<f64>,
    /// The weighted wall production per cell, added to (1 − m_p) P_k.
    pub wall_pk: Raw<f64>,
    /// Output ν_t.
    pub nu_t: Raw<f64>,
    /// Output T.
    pub t: Raw<f64>,
    /// Output L².
    pub l2: Raw<f64>,
    /// Output P_k.
    pub pk: Raw<f64>,
    /// Output S_u of k (empty: ν_t, T, L², P_k only).
    pub su_k: Raw<f64>,
    /// Output S_p of k.
    pub sp_k: Raw<f64>,
    /// Output S_u of ε.
    pub su_e: Raw<f64>,
    /// Output S_p of ε.
    pub sp_e: Raw<f64>,
    /// Output S_u of ζ.
    pub su_z: Raw<f64>,
    /// Output S_p of ζ.
    pub sp_z: Raw<f64>,
    /// Output: the right-hand side of the f equation.
    pub rhs_f: Raw<f64>,
    /// Output: S_p of the f (α) equation.
    pub sp_f: Raw<f64>,
    /// Output ν_t/σ_k.
    pub dk: Raw<f64>,
    /// Output ν_t/σ_ε.
    pub de: Raw<f64>,
    /// Output ν_t/σ_ζ.
    pub dz: Raw<f64>,
    /// Output 1 where the cap of ν_t/ν acted.
    pub capped: Raw<f64>,
    /// ∇k per cell (the elliptic-blending model; empty otherwise).
    pub grad_k: Raw3,
    /// ∇φ per cell (with `grad_k`).
    pub grad_z: Raw3,
    /// Σ_ijk (∂²U_i/∂x_j∂x_k)² per cell (with `grad_k`).
    pub g2: Raw<f64>,
    /// ∇·(ν_t/σ_k ∇k) per cell (with `grad_k`).
    pub dkt: Raw<f64>,
    /// The wall distance per cell.
    pub dist: Raw<f64>,
}

/// The k-ε-ζ-f model per cell ([`crate::turb::zeta_f::zeta_f_cell`]).
pub struct ZetaFCells;
impl Entry for ZetaFCells {
    type Args = ZetaFCellArgs;
    const NAME: &'static str = "zeta_f_cells";
    #[inline(always)]
    unsafe fn run(c: usize, a: &ZetaFCellArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = a.grad_u.v();
            let (wm, wp) = if a.wall_mask.len > 0 {
                (ld(a.wall_mask.get(), c), ld(a.wall_pk.get(), c))
            } else {
                (0.0, 0.0)
            };
            let inp = crate::turb::zeta_f::ZetaFCellInputs {
                k: ld(a.k.get(), c),
                e: ld(a.e.get(), c),
                zeta: ld(a.zeta.get(), c),
                f: ld(a.f.get(), c),
                grad_u: [g.x.at(c), g.y.at(c), g.z.at(c)],
                nu: if a.nu_cell.len > 0 {
                    ld(a.nu_cell.get(), c)
                } else {
                    a.nu
                },
                wall_m: wm,
                wall_pk: wp,
                grad_k: if a.grad_k.x.len > 0 {
                    a.grad_k.v().at(c)
                } else {
                    [0.0; 3]
                },
                grad_z: if a.grad_z.x.len > 0 {
                    a.grad_z.v().at(c)
                } else {
                    [0.0; 3]
                },
                g2: if a.g2.len > 0 { ld(a.g2.get(), c) } else { 0.0 },
                dkt: if a.dkt.len > 0 {
                    ld(a.dkt.get(), c)
                } else {
                    0.0
                },
                d: if a.dist.len > 0 {
                    ld(a.dist.get(), c)
                } else {
                    1e10
                },
            };
            let r = crate::turb::zeta_f::zeta_f_cell(&inp, &a.c, &a.o);
            *a.nu_t.out().add(c) = r.nu_t;
            *a.t.out().add(c) = r.t;
            *a.l2.out().add(c) = r.l2;
            *a.pk.out().add(c) = r.pk;
            *a.dk.out().add(c) = r.dk;
            *a.de.out().add(c) = r.de;
            *a.dz.out().add(c) = r.dz;
            *a.capped.out().add(c) = r.capped;
            if a.su_k.len > 0 {
                let (ak, ae, az, af) = if a.add_k.len > 0 {
                    (
                        ld(a.add_k.get(), c),
                        ld(a.add_e.get(), c),
                        ld(a.add_z.get(), c),
                        ld(a.add_f.get(), c),
                    )
                } else {
                    (0.0, 0.0, 0.0, 0.0)
                };
                let lo = crate::turb::zeta_f::ZETA_MIN;
                *a.su_k.out().add(c) = r.su_k + ak.max(0.0);
                *a.sp_k.out().add(c) = r.sp_k + ak.min(0.0) / inp.k.max(1e-300);
                *a.su_e.out().add(c) = r.su_e + ae.max(0.0);
                *a.sp_e.out().add(c) = r.sp_e + ae.min(0.0) / inp.e.max(1e-300);
                *a.su_z.out().add(c) = r.su_z + az.max(0.0);
                *a.sp_z.out().add(c) = r.sp_z + az.min(0.0) / inp.zeta.max(lo);
                *a.rhs_f.out().add(c) = r.rhs_f + af;
                *a.sp_f.out().add(c) = r.sp_f;
            }
        }
    }
}

/// Arguments of [`TurbFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct TurbFaceArgs {
    /// Mesh.
    pub m: Mesh,
    /// ν_t per cell.
    pub nu_t: Raw<f64>,
    /// σ_k ν_t per cell.
    pub dk: Raw<f64>,
    /// σ_ω ν_t per cell.
    pub dw: Raw<f64>,
    /// Laminar viscosity per face.
    pub mu: Raw<f64>,
    /// 1 on the wall faces (per boundary face): no turbulent diffusion
    /// through a smooth wall.
    pub wall: Raw<f64>,
    /// Factor of the molecular viscosity in the diffusivities of k, ε/ω
    /// and ζ: 1; ½ for the elliptic-blending model.
    pub mol: f64,
    /// Output μ + ρ ν_t per face.
    pub mu_eff: Raw<f64>,
    /// Output ρ(ν + σ_k ν_t) per face.
    pub gamma_k: Raw<f64>,
    /// Output ρ(ν + σ_ω ν_t) per face.
    pub gamma_w: Raw<f64>,
    /// ν_t/σ_ζ per cell (k-ε-ζ-f), or [`Raw::EMPTY`].
    pub dz: Raw<f64>,
    /// L² per cell (k-ε-ζ-f), or [`Raw::EMPTY`].
    pub l2: Raw<f64>,
    /// Output ρ(ν + ν_t/σ_ζ) per face (with `dz`).
    pub gamma_z: Raw<f64>,
    /// Output ρ L² per face (with `l2`): interpolated inside, the owner
    /// value on every boundary face (walls included, where f is fixed).
    pub gamma_f: Raw<f64>,
}

/// Face viscosities of the turbulent flow: linear interpolation of the
/// cell values on interior and periodic faces, the owner value on boundary
/// faces, none on walls.
pub struct TurbFaces;
impl Entry for TurbFaces {
    type Args = TurbFaceArgs;
    const NAME: &'static str = "turb_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &TurbFaceArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let ni = a.m.n_interior as usize;
            let p = ld(a.m.owner.get(), f) as usize;
            let mu = ld(a.mu.get(), f);
            let rho = ld(a.m.rho_face.get(), f);
            let wall = f >= ni && ld(a.wall.get(), f - ni) > 0.0;
            let (nt, gk, gw) = if f < ni {
                let n = ld(a.m.neighbour.get(), f) as usize;
                let w = ld(a.m.weight.get(), f);
                let wn = 1.0 - w;
                let at = |b: Raw<f64>| w * ld(b.get(), p) + wn * ld(b.get(), n);
                (at(a.nu_t), at(a.dk), at(a.dw))
            } else {
                let s = 1.0 - ld(a.wall.get(), f - ni);
                (
                    s * ld(a.nu_t.get(), p),
                    s * ld(a.dk.get(), p),
                    s * ld(a.dw.get(), p),
                )
            };
            // wall faces: μ_w of the wall treatment ([`WallFaces`])
            if !wall {
                *a.mu_eff.out().add(f) = mu + rho * nt;
            }
            let mm = a.mol * mu;
            *a.gamma_k.out().add(f) = mm + rho * gk;
            *a.gamma_w.out().add(f) = mm + rho * gw;
            if a.dz.len > 0 {
                let (gz, gf) = if f < ni {
                    let n = ld(a.m.neighbour.get(), f) as usize;
                    let w = ld(a.m.weight.get(), f);
                    let wn = 1.0 - w;
                    (
                        w * ld(a.dz.get(), p) + wn * ld(a.dz.get(), n),
                        w * ld(a.l2.get(), p) + wn * ld(a.l2.get(), n),
                    )
                } else {
                    let s = 1.0 - ld(a.wall.get(), f - ni);
                    (s * ld(a.dz.get(), p), ld(a.l2.get(), p))
                };
                *a.gamma_z.out().add(f) = mm + rho * gz;
                *a.gamma_f.out().add(f) = rho * gf;
            }
        }
    }
}

/// Arguments of [`ScalarTurbFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ScalarTurbArgs {
    /// Mesh.
    pub m: Mesh,
    /// Laminar viscosity μ per face.
    pub mu: Raw<f64>,
    /// μ + ρ ν_t per face (the face viscosity of the flow).
    pub mu_eff: Raw<f64>,
    /// ρ Γ per face: the molecular diffusivity of the scalar.
    pub gamma0: Raw<f64>,
    /// 1 on the wall faces (per boundary face).
    pub wall: Raw<f64>,
    /// q d_n per boundary face (fixed-flux faces; 0 elsewhere).
    pub fq: Raw<f64>,
    /// 0 on fixed-value faces, 1 otherwise (per boundary face).
    pub bc_a: Raw<f64>,
    /// Model of the turbulent Prandtl number ([`crate::turb::thermal`]).
    pub prt_model: u32,
    /// Its constant value.
    pub prt: f64,
    /// Output ρ(Γ + ν_t/Pr_t) per face.
    pub gamma: Raw<f64>,
    /// Output (per boundary face, flux faces off the walls): q d_n/Γ_f.
    pub bc_b: Raw<f64>,
}

/// The turbulent diffusivity of a scalar per face: ρ Γ + (μ_eff − μ)/Pr_t,
/// the owner value on boundary faces; walls keep ρ Γ ([`ScalarWallFaces`]).
pub struct ScalarTurbFaces;
impl Entry for ScalarTurbFaces {
    type Args = ScalarTurbArgs;
    const NAME: &'static str = "scalar_turb_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &ScalarTurbArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let ni = a.m.n_interior as usize;
            let g0 = ld(a.gamma0.get(), f);
            let wall = f >= ni && ld(a.wall.get(), f - ni) > 0.0;
            let mut g = g0;
            if !wall {
                let mu = ld(a.mu.get(), f);
                let mt = foehncraft_core::math::pmax(ld(a.mu_eff.get(), f) - mu, 0.0);
                let pr = mu / g0;
                let prt = crate::turb::thermal::prandtl_t(a.prt_model, a.prt, mt / mu, pr);
                g = g0 + mt / prt;
            }
            *a.gamma.out().add(f) = g;
            if f >= ni && !wall && ld(a.bc_a.get(), f - ni) > 0.5 {
                *a.bc_b.out().add(f - ni) = ld(a.fq.get(), f - ni) / g;
            }
        }
    }
}

/// Arguments of [`ScalarWallFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ScalarWallArgs {
    /// Mesh.
    pub m: Mesh,
    /// The wall faces (global face index).
    pub faces: Raw<u32>,
    /// y⁺ of the thermal law per wall face (`y_star` of [`WallFaces`]).
    pub y_plus: Raw<f64>,
    /// Laminar viscosity μ per face.
    pub mu: Raw<f64>,
    /// ρ Γ per face.
    pub gamma0: Raw<f64>,
    /// q d_n per boundary face.
    pub fq: Raw<f64>,
    /// 0 on fixed-value faces, 1 otherwise (per boundary face).
    pub bc_a: Raw<f64>,
    /// Thermal wall law ([`crate::turb::thermal::KADER`], …).
    pub law: u32,
    /// Pr_t of the log layer.
    pub prt: f64,
    /// The y⁺ at which T⁺ is taken at least: 11.06 with the scalable
    /// momentum treatment (its log law at ỹ⁺ = max(y⁺, 11.06), section 4.3
    /// of the stage-6 specification), 0 otherwise.
    pub y_floor: f64,
    /// Output Γ_w = h_w y_P per face (the wall faces).
    pub gamma: Raw<f64>,
    /// Output (per boundary face): q d_n/Γ_w on flux faces.
    pub bc_b: Raw<f64>,
    /// Output T⁺ per wall face (or [`Raw::EMPTY`]).
    pub t_plus: Raw<f64>,
}

/// The thermal wall function per wall face: the face diffusivity is
/// Γ_w = μ y*/T⁺(ỹ*), so that its conductance is h_w = ρ u*/T⁺(ỹ*), with
/// ỹ* = max(y*, `y_floor`): with the scalable treatment the thermal law
/// sees the same virtual wall distance as the momentum law. The velocity
/// scale u* = C_μ^¼ k_P^½ of the k-ε family stays finite where the wall
/// shear vanishes (a stagnation point or a reattachment line).
pub struct ScalarWallFaces;
impl Entry for ScalarWallFaces {
    type Args = ScalarWallArgs;
    const NAME: &'static str = "scalar_wall_faces";
    #[inline(always)]
    unsafe fn run(i: usize, a: &ScalarWallArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let f = ld(a.faces.get(), i) as usize;
            let fb = f - a.m.n_interior as usize;
            let mu = ld(a.mu.get(), f);
            let pr = mu / ld(a.gamma0.get(), f);
            let yp = ld(a.y_plus.get(), i);
            let yt = yp.max(a.y_floor);
            let g = crate::turb::thermal::wall_conductance_at(a.law, mu, yp, yt, pr, a.prt);
            *a.gamma.out().add(f) = g;
            if ld(a.bc_a.get(), fb) > 0.5 {
                *a.bc_b.out().add(fb) = ld(a.fq.get(), fb) / g;
            }
            if a.t_plus.len > 0 {
                *a.t_plus.out().add(i) = crate::turb::thermal::t_plus(a.law, yt, pr, a.prt);
            }
        }
    }
}

/// Arguments of [`WallFaces`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct WallFaceArgs {
    /// Mesh.
    pub m: Mesh,
    /// The wall faces (global face index).
    pub faces: Raw<u32>,
    /// Unit normal of each wall face (outward).
    pub normal: Raw3,
    /// y_P of each wall face.
    pub y: Raw<f64>,
    /// k_s of each wall face.
    pub ks: Raw<f64>,
    /// Cell velocity.
    pub u: Raw3,
    /// Boundary face velocity (per boundary face).
    pub ub: Raw3,
    /// ν.
    pub nu: f64,
    /// ν per cell (D194: several fluids), or [`Raw::EMPTY`]: `nu` in
    /// every cell.
    pub nu_cell: Raw<f64>,
    /// Treatment ([`crate::turb::wall::LOW_RE`], `AUTOMATIC`).
    pub treatment: u32,
    /// Constants.
    pub c: crate::turb::sst::SstConstants,
    /// Constants of the k-ε family (treatments `SCALABLE`, `STANDARD`,
    /// `COMPOUND`).
    pub ke: crate::turb::keps::KepsConstants,
    /// k per cell (the k-ε treatments).
    pub k: Raw<f64>,
    /// C_s of Cebeci & Bradshaw per wall face (the k-ε treatments).
    pub cs: Raw<f64>,
    /// u_τ of the last iteration per wall face (rough walls); may be the
    /// output `u_tau` (read before written, same index).
    pub u_tau_old: Raw<f64>,
    /// Output u_τ per wall face.
    pub u_tau: Raw<f64>,
    /// Output y⁺ per wall face.
    pub y_plus: Raw<f64>,
    /// Output per wall face: the y⁺ of the thermal wall law, y* = u* y_P/ν
    /// with u* = C_μ^¼ k_P^½ for the k-ε family and ζ-f/EB, the y⁺ of the
    /// blended u_τ for SST (section 5.2 of the stage-6 specification).
    pub y_star: Raw<f64>,
    /// Output ω per wall face.
    pub omega: Raw<f64>,
    /// Output: the viscous weight of the blending per wall face.
    pub blending: Raw<f64>,
    /// Output (or [`Raw::EMPTY`]): the viscous weight of the wall treatment
    /// per wall face, for the output (the wall field `blending` of the
    /// stage-6 specification, 2.3): 1 for the integration to the wall
    /// (`low-re`), 1 − tanh((y⁺/5)²) of the automatic treatment of SST,
    /// e^{−Γ} of the compound treatment, Γ = dU⁺/dy⁺ of `compound-2023`,
    /// 1 in the sublayer of the standard wall function and 0 above it,
    /// 0 for the scalable wall function. Not used by the solver.
    pub viscous_weight: Raw<f64>,
    /// Output: the weighted wall production per wall face (SST automatic:
    /// (1 − b) u_τ³/(κ y_P)).
    pub pk: Raw<f64>,
    /// Output: the weight m_p of the wall production per wall face (the
    /// cell keeps (1 − m_p) of its own production).
    pub prod_w: Raw<f64>,
    /// Output (per face): μ_w on the wall faces.
    pub mu_eff: Raw<f64>,
    /// Output (per boundary face): the wall value of k (bc_b of k).
    pub k_b: Raw<f64>,
    /// Output (per boundary face): the wall value of ω (bc_b of ω).
    pub w_b: Raw<f64>,
    /// ζ per cell (k-ε-ζ-f), or [`Raw::EMPTY`]: with it, `LOW_RE` is the
    /// integration to the wall of ζ-f and `COMPOUND` also writes f_w.
    pub zeta: Raw<f64>,
    /// Output (per boundary face, with `zeta`): the wall value of f.
    pub f_b: Raw<f64>,
    /// N_ζ of the ζ-f₀ form (1: the 2004 ζ-f, f_w = −2νζ_P/y_P²; 2: f₀_w
    /// = 0), with `zeta`.
    pub n_zeta: f64,
    /// Velocity gradient per cell (treatment `COMPOUND_2023`; empty: C_U
    /// = 0).
    pub grad_u: RawGrad,
    /// Pressure gradient per cell (treatment `COMPOUND_2023`; empty: no
    /// pressure part of C_U).
    pub grad_p: Raw3,
    /// The sample cell of each wall face (cut-cell walls, L13 stage 2c),
    /// `u32::MAX` for none, or [`Raw::EMPTY`]: the wall function takes the
    /// tangential velocity, k and distance of the sample cell for the wall
    /// shear stress, the near-wall cell keeps its own ε or ω.
    pub sample: Raw<u32>,
    /// The exact wall distance of the sample cell's centroid per wall face
    /// (with `sample`).
    pub y_sample: Raw<f64>,
}

/// The tangential speed of cell `c` relative to the wall face `fb` with the
/// unit normal `n`.
#[inline(always)]
unsafe fn wall_tangential(a: &WallFaceArgs, c: usize, fb: usize, n: [f64; 3]) -> f64 {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let up = a.u.v().at(c);
        let uw = a.ub.v().at(fb);
        let du = [up[0] - uw[0], up[1] - uw[1], up[2] - uw[2]];
        let dn = du[0] * n[0] + du[1] * n[1] + du[2] * n[2];
        let t = [du[0] - dn * n[0], du[1] - dn * n[1], du[2] - dn * n[2]];
        foehncraft_core::math::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2])
    }
}

/// The wall treatment per wall face ([`crate::turb::wall::wall_face`]).
pub struct WallFaces;
impl Entry for WallFaces {
    type Args = WallFaceArgs;
    const NAME: &'static str = "wall_faces";
    #[inline(always)]
    unsafe fn run(i: usize, a: &WallFaceArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let f = ld(a.faces.get(), i) as usize;
            let fb = f - a.m.n_interior as usize;
            let p = ld(a.m.owner.get(), f) as usize;
            let nu_p = if a.nu_cell.len > 0 {
                ld(a.nu_cell.get(), p)
            } else {
                a.nu
            };
            let n = a.normal.v().at(i);
            let up = a.u.v().at(p);
            let uw = a.ub.v().at(fb);
            let du = [up[0] - uw[0], up[1] - uw[1], up[2] - uw[2]];
            let dn = du[0] * n[0] + du[1] * n[1] + du[2] * n[2];
            let t = [du[0] - dn * n[0], du[1] - dn * n[1], du[2] - dn * n[2]];
            let ut = foehncraft_core::math::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
            let zf = a.zeta.len > 0;
            if zf && a.treatment == crate::turb::wall::LOW_RE {
                // k-ε-ζ-f integrated to the wall: k_w = 0, ζ_w = 0, ε_w
                // and f_w of the first cell (Dirichlet values)
                let (y, nu) = (ld(a.y.get(), i), nu_p);
                let (e_w, f_w) = crate::turb::wall::zeta_f_wall_values(
                    ld(a.k.get(), p),
                    ld(a.zeta.get(), p),
                    y,
                    nu,
                    1.0,
                    a.n_zeta,
                    a.ke.e_wall,
                );
                let u_tau = foehncraft_core::math::sqrt(nu * ut / y);
                let u_star = foehncraft_core::math::sqrt(foehncraft_core::math::sqrt(a.ke.c_mu))
                    * foehncraft_core::math::sqrt(ld(a.k.get(), p).max(0.0));
                *a.u_tau.out().add(i) = u_tau;
                *a.y_plus.out().add(i) = u_tau * y / nu;
                *a.y_star.out().add(i) = u_star * y / nu;
                *a.omega.out().add(i) = 0.0;
                *a.blending.out().add(i) = 1.0;
                if a.viscous_weight.len > 0 {
                    *a.viscous_weight.out().add(i) = 1.0;
                }
                *a.pk.out().add(i) = 0.0;
                *a.prod_w.out().add(i) = 0.0;
                *a.mu_eff.out().add(f) = ld(a.m.rho_face.get(), f) * nu;
                *a.k_b.out().add(fb) = 0.0;
                *a.w_b.out().add(fb) = e_w;
                *a.f_b.out().add(fb) = f_w;
                return;
            }
            if a.treatment >= crate::turb::wall::SCALABLE {
                // the k-ε family: ε fixed in the cell (weight 1), k and ε
                // with zero gradient at the wall
                // C_U = t·[(U·∇)U + ∇p/ρ] along the tangential direction t
                // of the relative velocity (Popovac 2023, eq. (14); steady
                // part)
                let mut cu = 0.0;
                if a.treatment == crate::turb::wall::COMPOUND_2023 && ut > 0.0 {
                    let inv = 1.0 / ut;
                    let tn = [t[0] * inv, t[1] * inv, t[2] * inv];
                    if a.grad_u.g[0].x.len > 0 {
                        let g = a.grad_u.v();
                        let gu = [g.x.at(p), g.y.at(p), g.z.at(p)];
                        for (ti, gi) in tn.iter().zip(gu.iter()) {
                            cu += ti * (up[0] * gi[0] + up[1] * gi[1] + up[2] * gi[2]);
                        }
                    }
                    if a.grad_p.x.len > 0 {
                        let gp = a.grad_p.v().at(p);
                        cu +=
                            (tn[0] * gp[0] + tn[1] * gp[1] + tn[2] * gp[2]) / ld(a.m.rho.get(), p);
                    }
                }
                let inp = crate::turb::wall::KepsWallInputs {
                    ut,
                    y: ld(a.y.get(), i),
                    nu: nu_p,
                    k: ld(a.k.get(), p),
                    ks: ld(a.ks.get(), i),
                    cs: ld(a.cs.get(), i),
                    cu,
                };
                let mut r = crate::turb::wall::keps_wall_face(&inp, &a.ke, a.treatment);
                if a.sample.len > 0 {
                    let q = ld(a.sample.get(), i);
                    if q != u32::MAX {
                        // the wall shear stress from the sample cell one
                        // layer out (L13 stage 2c); ε, the production and the
                        // wall viscosity stay those of this cell
                        let q = q as usize;
                        let mut si = inp;
                        si.ut = wall_tangential(a, q, fb, n);
                        si.y = ld(a.y_sample.get(), i);
                        si.k = ld(a.k.get(), q);
                        si.cu = 0.0;
                        if a.nu_cell.len > 0 {
                            si.nu = ld(a.nu_cell.get(), q);
                        }
                        let rs = crate::turb::wall::keps_wall_face(&si, &a.ke, a.treatment);
                        let (tau_s, tau_p) = (rs.u_tau * rs.u_tau, r.u_tau * r.u_tau);
                        if tau_p > 1e-300 {
                            r.pk_wall *= tau_s / tau_p;
                        }
                        if ut > 1e-300 {
                            r.nu_wall = tau_s * inp.y / ut;
                        }
                        r.u_tau = rs.u_tau;
                        r.y_plus = rs.u_tau * inp.y / nu_p;
                    }
                }
                *a.u_tau.out().add(i) = r.u_tau;
                *a.y_plus.out().add(i) = r.y_plus;
                *a.y_star.out().add(i) = r.y_star;
                *a.omega.out().add(i) = r.epsilon;
                *a.blending.out().add(i) = 0.0;
                if a.viscous_weight.len > 0 {
                    *a.viscous_weight.out().add(i) = r.viscous;
                }
                *a.pk.out().add(i) = r.pk_wall;
                *a.prod_w.out().add(i) = r.prod_weight;
                *a.mu_eff.out().add(f) = ld(a.m.rho_face.get(), f) * r.nu_wall;
                *a.k_b.out().add(fb) = 0.0;
                *a.w_b.out().add(fb) = 0.0;
                if zf {
                    // the compound f_w: e^{−Γ} (−2ν ζ_P/y_P²) (2023: Γ =
                    // dU⁺/dy⁺; f₀_w = 0)
                    let (_, f_w) = crate::turb::wall::zeta_f_wall_values(
                        inp.k,
                        ld(a.zeta.get(), p),
                        inp.y,
                        nu_p,
                        r.viscous,
                        a.n_zeta,
                        a.ke.e_wall,
                    );
                    *a.f_b.out().add(fb) = f_w;
                }
                return;
            }
            let inp = crate::turb::wall::WallInputs {
                ut,
                y: ld(a.y.get(), i),
                nu: nu_p,
                ks: ld(a.ks.get(), i),
                u_tau_old: a.u_tau_old.read(i),
            };
            let mut r = crate::turb::wall::wall_face(&inp, &a.c, a.treatment);
            if a.sample.len > 0 && a.treatment != crate::turb::wall::LOW_RE {
                let q = ld(a.sample.get(), i);
                if q != u32::MAX {
                    // the wall shear stress from the sample cell (L13
                    // stage 2c); ω and the wall values stay this cell's
                    let q = q as usize;
                    let mut si = inp;
                    si.ut = wall_tangential(a, q, fb, n);
                    si.y = ld(a.y_sample.get(), i);
                    if a.nu_cell.len > 0 {
                        si.nu = ld(a.nu_cell.get(), q);
                    }
                    let rs = crate::turb::wall::wall_face(&si, &a.c, a.treatment);
                    let (tau_s, tau_p) = (rs.u_tau * rs.u_tau, r.u_tau * r.u_tau);
                    if tau_p > 1e-300 {
                        r.pk_wall *= tau_s / tau_p;
                    }
                    if ut > 1e-300 {
                        r.nu_wall = tau_s * inp.y / ut;
                    }
                    r.u_tau = rs.u_tau;
                    r.y_plus = rs.u_tau * inp.y / nu_p;
                }
            }
            let m = r.prod_weight;
            *a.u_tau.out().add(i) = r.u_tau;
            *a.y_plus.out().add(i) = r.y_plus;
            *a.y_star.out().add(i) = r.y_plus;
            *a.omega.out().add(i) = r.omega;
            *a.blending.out().add(i) = r.blending;
            if a.viscous_weight.len > 0 {
                *a.viscous_weight.out().add(i) = r.blending;
            }
            *a.pk.out().add(i) = m * r.pk_wall;
            *a.prod_w.out().add(i) = m;
            *a.mu_eff.out().add(f) = ld(a.m.rho_face.get(), f) * r.nu_wall;
            *a.k_b.out().add(fb) = r.k_wall;
            *a.w_b.out().add(fb) = r.omega_wall;
        }
    }
}

/// Arguments of [`WallCells`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct WallCellArgs {
    /// The near-wall cells.
    pub cells: Raw<u32>,
    /// Start of the wall-face list of each near-wall cell (length + 1).
    pub ptr: Raw<u32>,
    /// Wall-face slots (indices into the wall-face arrays).
    pub slots: Raw<u32>,
    /// Area weight of each slot (the weights of a cell sum to 1).
    pub weight: Raw<f64>,
    /// ω per wall face.
    pub omega: Raw<f64>,
    /// Viscous weight of the blending per wall face (1: low-Re).
    pub blending: Raw<f64>,
    /// Weighted wall production per wall face.
    pub pk: Raw<f64>,
    /// Weight of the wall production per wall face.
    pub prod_w: Raw<f64>,
    /// Output: the area mean of `pk` per cell.
    pub wall_pk: Raw<f64>,
    /// Output: the area mean of `prod_w` per cell.
    pub prod_mask: Raw<f64>,
    /// Output: the fixed value per cell.
    pub fix_val: Raw<f64>,
    /// Output: the weight of the fixed value per cell, the area mean of
    /// 1 − blending over the smooth wall faces.
    pub fix_mask: Raw<f64>,
    /// The fixed variable, or [`Raw::EMPTY`]: if given, its value in the
    /// near-wall cells moves to the fixed one by the weight, the start
    /// value of the solve (the k-ε family: a start far from the fixed ε
    /// left the solve of the first iteration negative there).
    pub phi: Raw<f64>,
}

/// The area-weighted mean of the wall-face values of each near-wall cell,
/// in the stored order of its wall faces (stage 6 specification 4.1).
pub struct WallCells;
impl Entry for WallCells {
    type Args = WallCellArgs;
    const NAME: &'static str = "wall_cells";
    #[inline(always)]
    unsafe fn run(i: usize, a: &WallCellArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let c = ld(a.cells.get(), i) as usize;
            let lo = ld(a.ptr.get(), i) as usize;
            let hi = ld(a.ptr.get(), i + 1) as usize;
            let (mut v, mut m, mut p, mut q) = (0.0, 0.0, 0.0, 0.0);
            let mut j = lo;
            while j < hi {
                let s = ld(a.slots.get(), j) as usize;
                let w = ld(a.weight.get(), s);
                v += w * ld(a.omega.get(), s);
                m += w * (1.0 - ld(a.blending.get(), s));
                p += w * ld(a.pk.get(), s);
                q += w * ld(a.prod_w.get(), s);
                j += 1;
            }
            *a.fix_val.out().add(c) = v;
            *a.fix_mask.out().add(c) = m;
            if a.phi.len > 0 {
                let old = a.phi.read(c);
                *a.phi.out().add(c) = old + m * (v - old);
            }
            *a.wall_pk.out().add(c) = p;
            *a.prod_mask.out().add(c) = q;
        }
    }
}

/// Arguments of [`Floor`].
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct FloorArgs {
    /// Values (input and output).
    pub x: Raw<f64>,
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
    /// Relative lower bound: x ≥ rel · old where rel > 0 (k-ε-ζ-f: a
    /// solve may lower a value by at most this factor; 0: unused).
    pub rel: f64,
    /// The values before the solve (read only where rel > 0).
    pub old: Raw<f64>,
    /// Output: 1 where a bound acted, else 0.
    pub hit: Raw<f64>,
}

/// x ← min(max(x, lo, rel · old), hi), and where a bound acted (the positivity floors
/// of the turbulence variables, stage 6 specification 6.2; counted in the
/// health section).
pub struct Floor;
impl Entry for Floor {
    type Args = FloorArgs;
    const NAME: &'static str = "floor";
    #[inline(always)]
    unsafe fn run(i: usize, a: &FloorArgs) {
        // SAFETY: guaranteed by the caller; x is read before it is written.
        unsafe {
            let v = a.x.read(i);
            let lo = if a.rel > 0.0 {
                a.lo.max(a.rel * a.old.read(i))
            } else {
                a.lo
            };
            let r = v.max(lo).min(a.hi);
            *a.hit.out().add(i) = if r != v { 1.0 } else { 0.0 };
            *a.x.out().add(i) = r;
        }
    }
}

/// out = a x / y per entry (the volume flux ṁ_f/ρ_f of the
/// volume-of-fluid continuity, L3).
pub struct Quotient;
impl Entry for Quotient {
    type Args = LinArgs;
    const NAME: &'static str = "quotient";
    #[inline(always)]
    unsafe fn run(i: usize, a: &LinArgs) {
        // SAFETY: guaranteed by the caller; index i only.
        unsafe { *a.out.out().add(i) = a.a * a.x.read(i) / a.y.read(i) }
    }
}

/// Arguments of the volume-of-fluid entries (L3, [`crate::vof`]).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct VofArgs {
    /// Mesh.
    pub m: Mesh,
    /// Volume flux per face (m³/s).
    pub flux: Raw<f64>,
    /// α per cell (or the field the entry reads).
    pub alpha: Raw<f64>,
    /// α per face of the last [`VofFaces`].
    pub alpha_f: Raw<f64>,
    /// A cell gradient (∇α of the face schemes, ∇α̃ of the curvature).
    pub grad: Raw3,
    /// α of the inflow per boundary face.
    pub alpha_b: Raw<f64>,
    /// 1 where a boundary face takes the interface normal of its owner
    /// (open boundaries), 0 on walls and symmetry planes (90°).
    pub open: Raw<u32>,
    /// Curvature per cell.
    pub kappa: Raw<f64>,
    /// |∇α̃| per cell.
    pub gnorm: Raw<f64>,
    /// Face scheme ([`crate::vof::SCHEME_THINC`] …).
    pub scheme: u64,
    /// Sub-step Δt (s).
    pub dt: f64,
    /// β of THINC.
    pub beta: f64,
    /// σ (N/m) of [`VofTension`], δ of [`VofCurvature`].
    pub s: f64,
    /// Output.
    pub out: Raw<f64>,
    /// Second output (the clipped volume of [`VofUpdate`]).
    pub out2: Raw<f64>,
    /// Vector output ([`VofGradient`]).
    pub out3: Raw3,
}

/// Face value α_f of face `f` ([`crate::vof::face_alpha`]).
pub struct VofFaces;
impl Entry for VofFaces {
    type Args = VofArgs;
    const NAME: &'static str = "vof_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &VofArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let topo = a.m.topo();
            let cc = a.m.cell_centroid.v();
            let g = a.grad.v();
            let inp = crate::vof::FaceIn {
                topo: &topo,
                cell_centroid: &cc,
                volume: a.m.volume.get(),
                flux: a.flux.get(),
                alpha: a.alpha.get(),
                grad: &g,
                alpha_b: a.alpha_b.get(),
                scheme: a.scheme as u32,
                dt: a.dt,
                beta: a.beta,
            };
            *a.out.out().add(f) = crate::vof::face_alpha(f, &inp);
        }
    }
}

/// Explicit update of cell `c` ([`crate::vof::update`]): `out` the new α,
/// `out2` the clipped volume.
pub struct VofUpdate;
impl Entry for VofUpdate {
    type Args = VofArgs;
    const NAME: &'static str = "vof_update";
    #[inline(always)]
    unsafe fn run(c: usize, a: &VofArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let (v, clip) = crate::vof::update(
                c,
                &a.m.topo(),
                a.m.volume.get(),
                a.flux.get(),
                a.alpha_f.get(),
                a.alpha.get(),
                a.kappa.get(),
                a.dt,
            );
            *a.out.out().add(c) = v;
            *a.out2.out().add(c) = clip;
        }
    }
}

/// Outflow rate Σ_out |F|/V of cell `c` ([`crate::vof::outflow_rate`]).
pub struct VofOutflowRate;
impl Entry for VofOutflowRate {
    type Args = VofArgs;
    const NAME: &'static str = "vof_outflow_rate";
    #[inline(always)]
    unsafe fn run(c: usize, a: &VofArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) =
                crate::vof::outflow_rate(c, &a.m.topo(), a.m.volume.get(), a.flux.get())
        }
    }
}

/// Smoothed α̃ of cell `c` ([`crate::vof::smooth`]).
pub struct VofSmooth;
impl Entry for VofSmooth {
    type Args = VofArgs;
    const NAME: &'static str = "vof_smooth";
    #[inline(always)]
    unsafe fn run(c: usize, a: &VofArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) = crate::vof::smooth(
                c,
                &a.m.topo(),
                &a.m.area.v(),
                a.m.weight.get(),
                a.alpha.get(),
            )
        }
    }
}

/// Green–Gauss gradient of `alpha` in cell `c` into `out3`, its norm into
/// `out` ([`crate::vof::gauss_gradient`]).
pub struct VofGradient;
impl Entry for VofGradient {
    type Args = VofArgs;
    const NAME: &'static str = "vof_gradient";
    #[inline(always)]
    unsafe fn run(c: usize, a: &VofArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = crate::vof::gauss_gradient(
                c,
                &a.m.topo(),
                &a.m.area.v(),
                a.m.weight.get(),
                a.m.volume.get(),
                a.alpha.get(),
            );
            *a.out3.x.out().add(c) = g[0];
            *a.out3.y.out().add(c) = g[1];
            *a.out3.z.out().add(c) = g[2];
            *a.out.out().add(c) =
                foehncraft_core::Real::sqrt_of(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
        }
    }
}

/// Curvature of cell `c` ([`crate::vof::curvature`], δ = `s`).
pub struct VofCurvature;
impl Entry for VofCurvature {
    type Args = VofArgs;
    const NAME: &'static str = "vof_curvature";
    #[inline(always)]
    unsafe fn run(c: usize, a: &VofArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(c) = crate::vof::curvature(
                c,
                &a.m.topo(),
                &a.m.area.v(),
                a.m.weight.get(),
                a.m.volume.get(),
                &a.grad.v(),
                a.open.get(),
                a.s,
            )
        }
    }
}

/// The balanced pressure-force reconstruction of cell `c`
/// ([`crate::vof::reconstruct`]): φ = `alpha`, offsets = `alpha_f`, ρ_f =
/// `flux`, ρ = `kappa`, boundary values = `alpha_b`, added vector = `grad`,
/// cap = `s`, into `out3`.
pub struct VofReconstruct;
impl Entry for VofReconstruct {
    type Args = VofArgs;
    const NAME: &'static str = "vof_reconstruct";
    #[inline(always)]
    unsafe fn run(c: usize, a: &VofArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let v = crate::vof::reconstruct(
                c,
                &a.m.topo(),
                &a.m.area.v(),
                &a.m.face_centroid.v(),
                &a.m.cell_centroid.v(),
                a.alpha.get(),
                a.alpha_f.get(),
                a.flux.get(),
                a.kappa.get(),
                a.open.get(),
                a.alpha_b.get(),
                &a.grad.v(),
                a.s,
            );
            let o = &a.out3;
            *o.x.out().add(c) = v[0];
            *o.y.out().add(c) = v[1];
            *o.z.out().add(c) = v[2];
        }
    }
}

/// Surface-tension face offset of face `f` ([`crate::vof::tension`],
/// σ = `s`).
pub struct VofTension;
impl Entry for VofTension {
    type Args = VofArgs;
    const NAME: &'static str = "vof_tension";
    #[inline(always)]
    unsafe fn run(f: usize, a: &VofArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            *a.out.out().add(f) = crate::vof::tension(
                f,
                &a.m.topo(),
                a.alpha.get(),
                a.kappa.get(),
                a.gnorm.get(),
                a.s,
            )
        }
    }
}

/// Arguments of the two-fluid entries (L3 MP3, [`crate::twofluid`]).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct TwoFluidArgs {
    /// Mesh.
    pub m: Mesh,
    /// The constants of the mixture.
    pub k: crate::twofluid::Constants<f64>,
    /// α_g per cell.
    pub alpha: Raw<f64>,
    /// Liquid velocity.
    pub ul: Raw3,
    /// Gas velocity.
    pub ug: Raw3,
    /// The liquid velocity gradient.
    pub grad_ul: RawGrad,
    /// Wall distance per cell.
    pub wall_y: Raw<f64>,
    /// Unit wall normal per cell.
    pub wall_n: Raw3,
    /// ν_t of the liquid per cell (empty: laminar).
    pub nu_t: Raw<f64>,
    /// ∇α_g per cell.
    pub grad_alpha: Raw3,
    /// The cell size h per cell (the dispersion limiter).
    pub cell_h: Raw<f64>,
    /// α_g of the inflow per boundary face.
    pub alpha_b: Raw<f64>,
    /// 1 per inflow boundary face.
    pub inflow: Raw<u32>,
    /// α_g per face.
    pub alpha_gf: Raw<f64>,
    /// α_l per face.
    pub alpha_lf: Raw<f64>,
    /// μ_l + μ_t per face.
    pub mu_eff: Raw<f64>,
    /// μ_BIT per cell (output of [`TwoFluidCells`], input of
    /// [`TwoFluidFaces`]).
    pub mu_bit: Raw<f64>,
    /// D per cell (as `mu_bit`).
    pub disp: Raw<f64>,
    /// K V per cell (output).
    pub kv: Raw<f64>,
    /// (lift + wall) V per cell (output).
    pub fv: Raw3,
    /// b_g per cell (output).
    pub body_g: Raw3,
    /// b_l per cell (output).
    pub body_l: Raw3,
    /// The face viscosity of the gas (output).
    pub mu_g: Raw<f64>,
    /// The face viscosity of the liquid (output).
    pub mu_l: Raw<f64>,
    /// The offsets of the gas (output).
    pub off_g: Raw<f64>,
    /// The offsets of the liquid (output).
    pub off_l: Raw<f64>,
}

/// The interphase forces of cell `c` ([`crate::twofluid::cell_forces`])
/// and the cell forces b_g, b_l ([`crate::twofluid::cell_body`]).
pub struct TwoFluidCells;
impl Entry for TwoFluidCells {
    type Args = TwoFluidArgs;
    const NAME: &'static str = "two_fluid_cells";
    #[inline(always)]
    unsafe fn run(c: usize, a: &TwoFluidArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let g = a.grad_ul.g;
            let inp = crate::twofluid::CellIn {
                volume: a.m.volume.get(),
                alpha: a.alpha.get(),
                ul: a.ul.v(),
                ug: a.ug.v(),
                grad_ul: [g[0].v(), g[1].v(), g[2].v()],
                wall_y: a.wall_y.get(),
                wall_n: a.wall_n.v(),
                nu_t: a.nu_t.get(),
                h: a.cell_h.get(),
                grad_alpha: a.grad_alpha.v(),
            };
            let r = crate::twofluid::cell_forces(c, &inp, &a.k);
            *a.kv.out().add(c) = r.kv;
            *a.fv.x.out().add(c) = r.fv[0];
            *a.fv.y.out().add(c) = r.fv[1];
            *a.fv.z.out().add(c) = r.fv[2];
            *a.mu_bit.out().add(c) = r.mu_bit;
            *a.disp.out().add(c) = r.disp;
            let al = a.alpha.read(c);
            let ga = [
                a.grad_alpha.x.read(c),
                a.grad_alpha.y.read(c),
                a.grad_alpha.z.read(c),
            ];
            let bg = crate::twofluid::cell_body(&a.k, true, r.alpha_td, r.disp, ga);
            let bl = crate::twofluid::cell_body(&a.k, false, al, r.disp, ga);
            *a.body_g.x.out().add(c) = bg[0];
            *a.body_g.y.out().add(c) = bg[1];
            *a.body_g.z.out().add(c) = bg[2];
            *a.body_l.x.out().add(c) = bl[0];
            *a.body_l.y.out().add(c) = bl[1];
            *a.body_l.z.out().add(c) = bl[2];
        }
    }
}

/// The face viscosities and offsets of face `f`
/// ([`crate::twofluid::face_arrays`]).
pub struct TwoFluidFaces;
impl Entry for TwoFluidFaces {
    type Args = TwoFluidArgs;
    const NAME: &'static str = "two_fluid_faces";
    #[inline(always)]
    unsafe fn run(f: usize, a: &TwoFluidArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let topo = a.m.topo();
            let cc = a.m.cell_centroid.v();
            let fc = a.m.face_centroid.v();
            let sh = a.m.shift.v();
            let inp = crate::twofluid::FaceIn {
                topo: &topo,
                weight: a.m.weight.get(),
                cell_centroid: &cc,
                face_centroid: &fc,
                shift: &sh,
                alpha: a.alpha.get(),
                alpha_b: a.alpha_b.get(),
                alpha_gf: a.alpha_gf.get(),
                alpha_lf: a.alpha_lf.get(),
                mu_eff: a.mu_eff.get(),
                mu_bit: a.mu_bit.get(),
                disp: a.disp.get(),
                inflow: a.inflow.get(),
            };
            let r = crate::twofluid::face_arrays(f, &inp, &a.k);
            *a.mu_g.out().add(f) = r.mu_g;
            *a.mu_l.out().add(f) = r.mu_l;
            *a.off_g.out().add(f) = r.off_g;
            *a.off_l.out().add(f) = r.off_l;
        }
    }
}

/// Arguments of [`TwoFluidEliminate`].
#[repr(C)]
#[derive(Clone, Copy)]
pub struct TwoFluidElimArgs {
    /// Mesh.
    pub m: Mesh,
    /// α_g per cell.
    pub alpha: Raw<f64>,
    /// K V per cell.
    pub kv: Raw<f64>,
    /// The diagonals of the gas rows.
    pub diag_g: Raw3,
    /// The diagonals of the liquid rows.
    pub diag_l: Raw3,
    /// The solved gas velocity (in and out).
    pub ug: Raw3,
    /// The solved liquid velocity (in and out).
    pub ul: Raw3,
    /// The liquid velocity of the gas rows.
    pub ul_lag: Raw3,
    /// d, d^s, d^C of the gas (V/a without α in, D out).
    pub d_g: [Raw<f64>; 3],
    /// d, d^s, d^C of the liquid (as `d_g`).
    pub d_l: [Raw<f64>; 3],
}

/// The partial elimination of the drag of cell `c`
/// ([`crate::twofluid::eliminate`], [`crate::twofluid::eliminate_d`]).
pub struct TwoFluidEliminate;
impl Entry for TwoFluidEliminate {
    type Args = TwoFluidElimArgs;
    const NAME: &'static str = "two_fluid_eliminate";
    #[inline(always)]
    unsafe fn run(c: usize, a: &TwoFluidElimArgs) {
        // SAFETY: guaranteed by the caller.
        unsafe {
            let kv = a.kv.read(c);
            let alpha = a.alpha.read(c);
            let v = a.m.volume.read(c);
            let comps = [
                (a.diag_g.x, a.diag_l.x, a.ug.x, a.ul.x, a.ul_lag.x),
                (a.diag_g.y, a.diag_l.y, a.ug.y, a.ul.y, a.ul_lag.y),
                (a.diag_g.z, a.diag_l.z, a.ug.z, a.ul.z, a.ul_lag.z),
            ];
            for (dg, dl, ug, ul, lag) in comps {
                let (g, l) = crate::twofluid::eliminate(
                    kv,
                    dg.read(c),
                    dl.read(c),
                    ug.read(c),
                    ul.read(c),
                    lag.read(c),
                );
                *ug.out().add(c) = g;
                *ul.out().add(c) = l;
            }
            for (g, l) in a.d_g.into_iter().zip(a.d_l) {
                let (x, y) = crate::twofluid::eliminate_d(v, kv, alpha, g.read(c), l.read(c));
                *g.out().add(c) = x;
                *l.out().add(c) = y;
            }
        }
    }
}

send_sync!(
    LinArgs,
    VofArgs,
    TwoFluidArgs,
    TwoFluidElimArgs,
    VecArgs,
    PseudoArgs,
    BlockArgs,
    DotArgs,
    CloseArgs,
    ConstArgs,
    BoundaryValuesArgs,
    GradientArgs,
    LimiterArgs,
    MlpArgs,
    MassFluxArgs,
    MomentumFaceArgs,
    MomentumRowArgs,
    OperatorArgs,
    ResidualScaleArgs,
    UnaryArgs,
    PressureRowArgs,
    CorrectFluxArgs,
    CorrectVelocityArgs,
    NeighbourArgs,
    DeltaFluxArgs,
    ProductArgs,
    ScalarFaceArgs,
    ScalarRowArgs,
    ScalarScaleArgs,
    SstCellArgs,
    CurvatureArgs,
    QcrArgs,
    QuadGradArgs,
    KepsCellArgs,
    ZetaFCellArgs,
    TurbFaceArgs,
    WallFaceArgs,
    WallCellArgs,
    FloorArgs
);

/// Names of all entries (the PTX module has `flow_<name>` for each).
pub const NAMES: [&str; 69] = [
    Lincomb::NAME,
    ScaledDifference::NAME,
    SquaresAdd::NAME,
    ScaleVector::NAME,
    SchurDiagonal::NAME,
    Product::NAME,
    VertexBounds::NAME,
    MlpLimiter::NAME,
    PseudoCoefficient::NAME,
    Fill::NAME,
    Subtract::NAME,
    ScaleEntry::NAME,
    BlockSumSq::NAME,
    BlockDot::NAME,
    BlockMinMax::NAME,
    BoundaryConstants::NAME,
    BoundaryValues::NAME,
    Gradient::NAME,
    Limiter::NAME,
    MassFlux::NAME,
    MomentumFaces::NAME,
    MomentumRows::NAME,
    MomentumOperator::NAME,
    ResidualScale::NAME,
    SetDiagonal::NAME,
    DiagonalInverse::NAME,
    NetOutflow::NAME,
    Throughflow::NAME,
    Courant::NAME,
    PressureFaces::NAME,
    CloseBackflow::NAME,
    PressureRows::NAME,
    CorrectFluxes::NAME,
    CorrectVelocity::NAME,
    NeighbourCorrection::NAME,
    DeltaFlux::NAME,
    Transpose::NAME,
    MatrixProduct::NAME,
    TransposedProduct::NAME,
    ScalarFaces::NAME,
    ScalarRows::NAME,
    ScalarScale::NAME,
    SstCells::NAME,
    StrainCells::NAME,
    ConvectiveDerivative::NAME,
    CurvatureCells::NAME,
    QcrCells::NAME,
    ExplicitStressFaces::NAME,
    QuadraticGradient::NAME,
    KepsCells::NAME,
    ZetaFCells::NAME,
    TurbFaces::NAME,
    ScalarTurbFaces::NAME,
    ScalarWallFaces::NAME,
    WallFaces::NAME,
    WallCells::NAME,
    Floor::NAME,
    Quotient::NAME,
    VofFaces::NAME,
    VofUpdate::NAME,
    VofOutflowRate::NAME,
    VofSmooth::NAME,
    VofGradient::NAME,
    VofCurvature::NAME,
    VofTension::NAME,
    VofReconstruct::NAME,
    TwoFluidCells::NAME,
    TwoFluidFaces::NAME,
    TwoFluidEliminate::NAME,
];

#[cfg(test)]
mod tests {
    use super::*;
    use core::mem::{align_of, size_of};

    /// CUDA limits the parameters of a kernel to 4096 bytes (Volta and
    /// newer accept more only with CUDA 12.1 and later); the count `n` takes
    /// another 8.
    #[test]
    fn argument_blocks_fit_the_kernel_parameter_limit() {
        let sizes = [
            size_of::<LinArgs>(),
            size_of::<PseudoArgs>(),
            size_of::<BlockArgs>(),
            size_of::<DotArgs>(),
            size_of::<CloseArgs>(),
            size_of::<ConstArgs>(),
            size_of::<BoundaryValuesArgs>(),
            size_of::<GradientArgs>(),
            size_of::<LimiterArgs>(),
            size_of::<MlpArgs>(),
            size_of::<MassFluxArgs>(),
            size_of::<MomentumFaceArgs>(),
            size_of::<MomentumRowArgs>(),
            size_of::<OperatorArgs>(),
            size_of::<ResidualScaleArgs>(),
            size_of::<UnaryArgs>(),
            size_of::<PressureRowArgs>(),
            size_of::<CorrectFluxArgs>(),
            size_of::<CorrectVelocityArgs>(),
            size_of::<NeighbourArgs>(),
            size_of::<DeltaFluxArgs>(),
            size_of::<ProductArgs>(),
            size_of::<ScalarFaceArgs>(),
            size_of::<ScalarRowArgs>(),
            size_of::<ScalarScaleArgs>(),
            size_of::<SstCellArgs>(),
            size_of::<KepsCellArgs>(),
            size_of::<ZetaFCellArgs>(),
            size_of::<TurbFaceArgs>(),
            size_of::<WallFaceArgs>(),
            size_of::<WallCellArgs>(),
            size_of::<FloorArgs>(),
            size_of::<QcrArgs>(),
            size_of::<QuadGradArgs>(),
            size_of::<TwoFluidArgs>(),
            size_of::<TwoFluidElimArgs>(),
            size_of::<VofArgs>(),
        ];
        for s in sizes {
            assert!(s + 8 <= 4096, "argument block of {s} bytes");
        }
    }

    /// The layout the PTX entries see: 64-bit addresses and lengths, 8-byte
    /// alignment, no padding inside an array descriptor.
    #[test]
    fn array_descriptors_are_two_words() {
        assert_eq!(size_of::<Raw<f64>>(), 16);
        assert_eq!(size_of::<Raw<u32>>(), 16);
        assert_eq!(align_of::<Raw<f64>>(), 8);
        assert_eq!(size_of::<Raw3>(), 48);
        assert_eq!(
            size_of::<Mesh>(),
            3 * 8
                + 4 * 16
                + 3 * 48
                + 16
                + 48
                + 16
                + 16
                + 2 * 48
                + 8
                + 2 * 16
                + 4 * 16
                + 16
                + 16
                + 2 * 16
                + 16
                + 48
                + 16
                + 8
                + 48
                + 2 * 16
                + 3 * 16
                // body (VOF, D-weighted Rhie–Chow)
                + 48
        );
    }

    /// Every entry body of the element-wise operations computes what the
    /// documentation says, on host memory.
    #[test]
    fn element_wise_entries() {
        let x = [1.0, 2.0, 3.0];
        let mut y = [10.0, 20.0, 30.0];
        let r = |v: &[f64]| Raw::<f64>::new(v.as_ptr() as u64, v.len());
        let a = LinArgs {
            a: 2.0,
            b: -1.0,
            x: r(&x),
            y: r(&y),
            out: Raw::new(y.as_mut_ptr() as u64, 3),
        };
        for i in 0..3 {
            // SAFETY: index i of arrays of length 3; out is y (allowed).
            unsafe { Lincomb::run(i, &a) };
        }
        assert_eq!(y, [-8.0, -16.0, -24.0]);
        // out = a x y (D194)
        let mut z = [0.0; 3];
        let a = LinArgs {
            a: 0.5,
            b: 0.0,
            x: r(&x),
            y: r(&y),
            out: Raw::new(z.as_mut_ptr() as u64, 3),
        };
        for i in 0..3 {
            // SAFETY: index i of arrays of length 3.
            unsafe { Product::run(i, &a) };
        }
        assert_eq!(z, [-4.0, -16.0, -36.0]);
        assert_eq!(NAMES.len(), 69);
        // the elliptic-blending helpers: out += |v|², out3 = s v
        let (vx, vy, vz) = ([1.0, 0.0], [2.0, -3.0], [2.0, 4.0]);
        let s = [0.5, -2.0];
        let mut acc = [1.0, 0.0];
        let (mut ox, mut oy, mut oz) = ([0.0; 2], [0.0; 2], [0.0; 2]);
        let w = |v: &mut [f64]| Raw::<f64>::new(v.as_mut_ptr() as u64, v.len());
        let a = VecArgs {
            s: r(&s),
            v: Raw3 {
                x: r(&vx),
                y: r(&vy),
                z: r(&vz),
            },
            out3: Raw3 {
                x: w(&mut ox),
                y: w(&mut oy),
                z: w(&mut oz),
            },
            out: w(&mut acc),
        };
        for i in 0..2 {
            // SAFETY: index i of arrays of length 2.
            unsafe {
                SquaresAdd::run(i, &a);
                ScaleVector::run(i, &a);
            }
        }
        assert_eq!(acc, [10.0, 25.0]);
        assert_eq!((ox, oy, oz), ([0.5, 0.0], [1.0, 6.0], [1.0, -8.0]));
    }

    /// The pseudo-time kernels of the coupled algorithm (D176 to D178).
    #[test]
    fn pseudo_time_kernels() {
        use crate::scalar::{pseudo_coefficient, residual_parts};
        // D̂ = min(D^C, 1/(a/D^s + b)): the bound holds where D^C is large
        // (fluxes that do not satisfy continuity) and D^C where it is small
        let ds = [2.0, 2.0, 2.0];
        let dc = [1.0, 1.0e30, 5.0];
        let mut out = [0.0; 3];
        let r = |v: &[f64]| Raw::<f64>::new(v.as_ptr() as u64, v.len());
        let args = LinArgs {
            a: 0.02,
            b: 0.5,
            x: r(&ds),
            y: r(&dc),
            out: Raw::new(out.as_mut_ptr() as u64, 3),
        };
        for i in 0..3 {
            // SAFETY: index i of arrays of length 3.
            unsafe { SchurDiagonal::run(i, &args) };
        }
        let bound = 1.0 / (0.02 / 2.0 + 0.5);
        assert_eq!(out, [1.0, bound, bound]);
        assert_eq!(crate::flow::schur_diagonal(0.0, 0.5, 2.0, 5.0), 2.0);

        // ρV max(1/Δτ, f·num/den): the smaller of the two steps
        assert_eq!(pseudo_coefficient(3.0, 0.2, 2.0, 1.0, 4.0), 3.0 * 0.5);
        assert_eq!(pseudo_coefficient(3.0, 0.2, 2.0, 0.1, 4.0), 3.0 * 0.2);
        // no time scale (k = 0): the global step; no global step: the scale
        assert_eq!(pseudo_coefficient(3.0, 0.2, 2.0, 1.0, 0.0), 3.0 * 0.2);
        assert_eq!(pseudo_coefficient(3.0, 0.0, 2.0, 1.0, 4.0), 3.0 * 0.5);

        // the residual scale of a scalar row does not contain the
        // pseudo-time term: a row with the coefficient τ on the diagonal
        // and τ φ on the right-hand side has the parts of the row without
        let (phi, alpha) = (1.7_f64, 0.9_f64);
        let (ap, b, tau) = (4.0_f64, 3.0_f64, 2.5_f64);
        let diag = |t: f64| (ap + t) / alpha;
        let rhs = |t: f64| b + t * phi + (diag(t) - (ap + t)) * phi;
        let plain = residual_parts(phi, rhs(0.0), diag(0.0), alpha, 0.0);
        let with = residual_parts(phi, rhs(tau), diag(tau), alpha, tau);
        assert!((plain.0 - with.0).abs() < 1e-12 && (plain.1 - with.1).abs() < 1e-12);
        assert!((plain.0 - ap * phi).abs() < 1e-12 && (plain.1 - b).abs() < 1e-12);
    }
}
