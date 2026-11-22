//! Where the loops of the flow solver run.
//!
//! The algorithms ([`crate::simplec`], [`crate::transient`]) are written once,
//! as a sequence of operations on device buffers ([`Buf`]); a
//! [`FlowDevice`] carries them out. The CPU device runs the single-source
//! bodies of `foehncraft_kernels::flow` in multi-threaded loops; a GPU device
//! (next step of stage 5) launches the same bodies compiled to PTX
//! (`gpu/ptx-kernels`) on device-resident data behind the same trait. Every
//! operation writes each output from one index only, and the reductions use
//! the fixed blocks of `foehncraft_core::repro`, so the result does not depend on
//! the number of threads, and CPU and GPU agree bitwise (D31).
//!
//! A device holds the static data of the run (mesh, geometry, least-squares
//! weights, boundary-condition arrays, matrix structure) and the linear
//! solvers (a backend of `foehncraft-fv`); the field buffers are allocated by
//! the algorithm and addressed by handle.

pub mod cpu;
pub mod kernel;

use crate::bc::BoundaryArrays;
use anyhow::Result;
use foehncraft_fv::amg::{AmgSettings, ReductionSettings};
use foehncraft_fv::disc::{Discretisation, Soa3};
use foehncraft_fv::solver::Preconditioner;

/// Handle of a device buffer of f64 values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Buf(pub(crate) u32);

impl Buf {
    /// The handle of the `i`-th buffer of an executor (executors number
    /// their buffers from 0 in the order of allocation).
    pub fn new(i: usize) -> Self {
        Self(i as u32)
    }

    /// Index of the buffer in its executor.
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Three buffers, one per component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Buf3 {
    /// x components.
    pub x: Buf,
    /// y components.
    pub y: Buf,
    /// z components.
    pub z: Buf,
}

impl Buf3 {
    /// The three handles in order.
    pub fn all(self) -> [Buf; 3] {
        [self.x, self.y, self.z]
    }
}

/// Which least-squares weights a gradient uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Lsq {
    /// Weights of the velocity (outlet and symmetry stencil points
    /// projected).
    Velocity,
    /// Weights of the pressure (wall and symmetry points projected).
    Pressure,
}

/// ε² of the limiter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Eps2 {
    /// Per cell, Venkatakrishnan (stored on the device at creation).
    PerCell,
    /// One value, Venkatakrishnan–Wang.
    Single(f64),
    /// Per cell, from a buffer of the caller (a limiter constant K of its
    /// own, e.g. a scalar's).
    Cells(Buf),
}

/// Inputs and outputs of [`FlowDevice::boundary_values`].
#[derive(Debug, Clone, Copy)]
pub struct BoundaryArgs {
    /// Cell velocity.
    pub u: Buf3,
    /// Velocity gradients.
    pub grad_u: [Buf3; 3],
    /// Cell pressure.
    pub p: Buf,
    /// Pressure gradient.
    pub grad_p: Buf3,
    /// Face mass flux.
    pub mass_flux: Buf,
    /// Prescribed boundary velocity (faces with the fixed rule).
    pub ub_set: Buf3,
    /// New boundary velocity (output).
    pub ub: Buf3,
    /// New boundary pressure (output).
    pub pb: Buf,
}

/// Inputs of [`FlowDevice::momentum_faces`].
#[derive(Debug, Clone, Copy)]
pub struct MomentumFaceArgs {
    /// Convecting face mass flux.
    pub mass_flux: Buf,
    /// Cell velocity.
    pub u: Buf3,
    /// Velocity of the central correction (`u`, or extrapolated in time).
    pub uc: Buf3,
    /// Velocity gradients.
    pub grad_u: [Buf3; 3],
    /// Limiter per component.
    pub psi: Buf3,
    /// Boundary velocity.
    pub ub: Buf3,
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
    /// 1 for its trace-free form (cell-quality remediation, D50).
    pub trace_free: f64,
    /// Output c_P.
    pub cp: Buf,
    /// Output c_N.
    pub cn: Buf,
    /// Output b_f per component.
    pub bf: Buf3,
}

/// Time terms of a momentum row.
#[derive(Debug, Clone, Copy)]
pub struct TimeTerms {
    /// a₀, −a₁ and −a₂ of the backward difference; the rows take
    /// ρ_P a / Δt with the cell density.
    pub a: [f64; 3],
    /// Δt.
    pub dt: f64,
    /// uⁿ.
    pub u_old: Buf3,
    /// uⁿ⁻¹ (used when `a[2]` is not zero).
    pub u_old2: Buf3,
    /// 1 − θ and the old spatial operator (Crank–Nicolson); weight 0: none.
    pub explicit: (f64, Buf3),
}

/// Inputs and outputs of [`FlowDevice::momentum_rows`].
#[derive(Debug, Clone, Copy)]
pub struct MomentumRowArgs {
    /// Face coefficients.
    pub cp: Buf,
    /// See `cp`.
    pub cn: Buf,
    /// Deferred parts.
    pub bf: Buf3,
    /// Pressure gradient.
    pub grad_p: Buf3,
    /// Body force per unit volume.
    pub body_force: [f64; 3],
    /// Current velocity (relaxation).
    pub u: Buf3,
    /// α_u.
    pub alpha: f64,
    /// θ.
    pub theta: f64,
    /// Time terms (None: steady).
    pub time: Option<TimeTerms>,
    /// Output: matrix values.
    pub val: Buf,
    /// Output: right-hand sides.
    pub rhs: Buf3,
    /// Output: D.
    pub d: Buf,
    /// Output: D^s.
    pub ds: Buf,
    /// Output: D^C.
    pub dc: Buf,
    /// Output: relaxed diagonal per component.
    pub diag: Buf3,
}

/// Inputs of [`FlowDevice::mass_flux`].
#[derive(Debug, Clone, Copy)]
pub struct FluxArgs {
    /// Cell velocity.
    pub u: Buf3,
    /// Velocity gradients.
    pub grad_u: [Buf3; 3],
    /// Cell pressure.
    pub p: Buf,
    /// Pressure gradient.
    pub grad_p: Buf3,
    /// Boundary pressure.
    pub pb: Buf,
    /// Boundary velocity.
    pub ub: Buf3,
    /// D per cell; None: the plain interpolated flux (no pressure terms).
    pub d: Option<Buf>,
    /// Relaxation term per face.
    pub relax: Buf,
    /// Time-step term: (D^s, Y, α_u).
    pub time: Option<(Buf, Buf, f64)>,
    /// The velocity limiter of the bounded face velocity (cell-quality
    /// remediation, D50); None: the face velocity with the interpolated
    /// gradient.
    pub bounded: Option<Buf3>,
}

/// Inputs and outputs of [`FlowDevice::scalar_faces`].
#[derive(Debug, Clone, Copy)]
pub struct ScalarFaceArgs {
    /// Boundary coefficient a per boundary face.
    pub bc_a: Buf,
    /// Boundary constant b per boundary face.
    pub bc_b: Buf,
    /// Jump per periodic face (`None`: plain periodic).
    pub jump: Option<Buf>,
    /// Convecting face mass flux.
    pub mass_flux: Buf,
    /// ρΓ per face (kg/(m s)).
    pub gamma: Buf,
    /// Cell values.
    pub phi: Buf,
    /// Cell gradient.
    pub grad: Buf3,
    /// Limiter per cell.
    pub psi: Buf,
    /// 1 linear upwind, 0 upwind.
    pub blend: f64,
    /// Fraction of the bounded central scheme.
    pub central: f64,
    /// β of its bounding.
    pub beta: f64,
    /// 1 with the non-orthogonal correction.
    pub nonorth: f64,
    /// Output c_P.
    pub cp: Buf,
    /// Output c_N.
    pub cn: Buf,
    /// Output b_f.
    pub bf: Buf,
}

/// Time terms of a scalar row.
#[derive(Debug, Clone, Copy)]
pub struct ScalarTime {
    /// 1/Δt.
    pub inv_dt: f64,
    /// a₀, a₁, a₂ of the backward difference.
    pub a: [f64; 3],
    /// φⁿ.
    pub old: Buf,
    /// φⁿ⁻¹ (read when a₂ ≠ 0).
    pub old2: Buf,
    /// Bounds of the bounded BDF2 (±∞: none).
    pub bounds: (f64, f64),
}

/// Inputs and outputs of [`FlowDevice::scalar_rows`].
#[derive(Debug, Clone, Copy)]
pub struct ScalarRowArgs {
    /// Face coefficients.
    pub cp: Buf,
    /// See `cp`.
    pub cn: Buf,
    /// Deferred parts.
    pub bf: Buf,
    /// Convecting face mass flux.
    pub mass_flux: Buf,
    /// Current cell values.
    pub phi: Buf,
    /// S_u per cell (the unit of φ per second), if any.
    pub su: Option<Buf>,
    /// S_p ≤ 0 per cell (1/s), if any.
    pub sp: Option<Buf>,
    /// α.
    pub alpha: f64,
    /// Time terms (None: steady).
    pub time: Option<ScalarTime>,
    /// Fixed values: (mask, values) per cell (the ω of SST in the near-wall
    /// cells), if any.
    pub fix: Option<(Buf, Buf)>,
    /// A positive variable: the deferred corrections split by sign (k, ω).
    pub positive: bool,
    /// Pseudo-time coefficient ρ V/Δτ per cell of a steady iteration
    /// ([`FlowDevice::pseudo_coefficient`]), if any.
    pub pseudo: Option<Buf>,
    /// The heat capacity ρ c_p per cell of the energy equation (stage
    /// 7a1), in place of the density in the time and source terms; `None`:
    /// the density.
    pub capacity: Option<Buf>,
    /// Output: matrix values.
    pub val: Buf,
    /// Output: right-hand side.
    pub rhs: Buf,
    /// Output: relaxed diagonal.
    pub diag: Buf,
    /// Output: 1 where the bounded BDF2 fell back to first order.
    pub order: Buf,
}

/// The static arrays of the wall treatment (stage 6, D134), uploaded once
/// by [`FlowDevice::wall_setup`].
#[derive(Debug, Clone, Default)]
pub struct WallStatic {
    /// The wall faces (global face index), in face order.
    pub faces: Vec<u32>,
    /// The near-wall cells (owners of the wall faces), ascending.
    pub cells: Vec<u32>,
    /// Start of each near-wall cell's slots (length cells + 1).
    pub ptr: Vec<u32>,
    /// Wall-face slots of the cells (indices into `faces`), in face order.
    pub slots: Vec<u32>,
}

/// Inputs and outputs of [`FlowDevice::sst_cells`].
#[derive(Debug, Clone, Copy)]
pub struct SstCellOps {
    /// Velocity gradient.
    pub grad_u: [Buf3; 3],
    /// k.
    pub k: Buf,
    /// ω.
    pub w: Buf,
    /// ∇k.
    pub grad_k: Buf3,
    /// ∇ω.
    pub grad_w: Buf3,
    /// Wall distance.
    pub d: Buf,
    /// ν.
    pub nu: f64,
    /// ν per cell (D194: several fluids); `None`: `nu` everywhere.
    pub nu_cell: Option<Buf>,
    /// Constants.
    pub c: foehncraft_kernels::turb::sst::SstConstants,
    /// Options.
    pub o: foehncraft_kernels::turb::sst::SstOptions,
    /// Added sources of k and ω, if any.
    pub add: Option<(Buf, Buf)>,
    /// The weight and the value of the wall production per cell (the
    /// automatic wall treatment), if any.
    pub wall: Option<(Buf, Buf)>,
    /// f_r1 of the curvature correction per cell (D170), if any.
    pub fr1: Option<Buf>,
    /// Output ν_t.
    pub nu_t: Buf,
    /// Output F₁.
    pub f1: Buf,
    /// Output P̃_k.
    pub pk: Buf,
    /// Outputs S_u, S_p of k and of ω (None: ν_t, F₁, P̃_k only).
    pub sources: Option<[Buf; 4]>,
    /// Output σ_k ν_t.
    pub dk: Buf,
    /// Output σ_ω ν_t.
    pub dw: Buf,
    /// Output: 1 where the ν_t cap acted.
    pub capped: Buf,
}

/// Inputs and outputs of [`FlowDevice::curvature_cells`] (D170).
#[derive(Debug, Clone, Copy)]
pub struct CurvatureOps {
    /// Velocity gradient.
    pub grad_u: [Buf3; 3],
    /// DS_ij/Dt (xx, yy, zz, xy, yz, zx).
    pub dsdt: [Buf; 6],
    /// ω.
    pub w: Buf,
    /// Ω^rot per cell (None: `rot` in every cell).
    pub omega_rot: Option<Buf3>,
    /// See `omega_rot`.
    pub rot: [f64; 3],
    /// Rotation of the frame of a relative velocity.
    pub spin: [f64; 3],
    /// Constants.
    pub c: foehncraft_kernels::turb::sst::SstConstants,
    /// Width of the smooth limits.
    pub smooth: f64,
    /// Output f_r1.
    pub out: Buf,
}

/// Inputs and outputs of [`FlowDevice::keps_cells`].
#[derive(Debug, Clone, Copy)]
pub struct KepsCellOps {
    /// Velocity gradient.
    pub grad_u: [Buf3; 3],
    /// k.
    pub k: Buf,
    /// ε.
    pub e: Buf,
    /// ν.
    pub nu: f64,
    /// ν per cell (D194: several fluids); `None`: `nu` everywhere.
    pub nu_cell: Option<Buf>,
    /// Constants.
    pub c: foehncraft_kernels::turb::keps::KepsConstants,
    /// Options.
    pub o: foehncraft_kernels::turb::keps::KepsOptions,
    /// Added sources of k and ε, if any.
    pub add: Option<(Buf, Buf)>,
    /// The weight and the weighted value of the wall production per cell,
    /// if any.
    pub wall: Option<(Buf, Buf)>,
    /// Output ν_t.
    pub nu_t: Buf,
    /// Output C_μ.
    pub c_mu: Buf,
    /// Output P_k.
    pub pk: Buf,
    /// Outputs S_u, S_p of k and of ε (None: ν_t, C_μ, P_k only).
    pub sources: Option<[Buf; 4]>,
    /// Output ν_t/σ_k.
    pub dk: Buf,
    /// Output ν_t/σ_ε.
    pub de: Buf,
    /// Output: 1 where the ν_t cap acted.
    pub capped: Buf,
}

/// Inputs and outputs of [`FlowDevice::zeta_f_cells`].
#[derive(Debug, Clone, Copy)]
pub struct ZetaFCellOps {
    /// Velocity gradient.
    pub grad_u: [Buf3; 3],
    /// k.
    pub k: Buf,
    /// ε.
    pub e: Buf,
    /// ζ.
    pub zeta: Buf,
    /// f.
    pub f: Buf,
    /// ν.
    pub nu: f64,
    /// ν per cell (D194: several fluids); `None`: `nu` everywhere.
    pub nu_cell: Option<Buf>,
    /// Constants.
    pub c: foehncraft_kernels::turb::zeta_f::ZetaFConstants,
    /// Options.
    pub o: foehncraft_kernels::turb::zeta_f::ZetaFOptions,
    /// Added sources of k, ε and ζ and the added right-hand side of f, if
    /// any.
    pub add: Option<[Buf; 4]>,
    /// The weight and the weighted value of the wall production per cell,
    /// if any.
    pub wall: Option<(Buf, Buf)>,
    /// The wall distance per cell.
    pub dist: Buf,
    /// Output ν_t.
    pub nu_t: Buf,
    /// Output T.
    pub t: Buf,
    /// Output L².
    pub l2: Buf,
    /// Output P_k.
    pub pk: Buf,
    /// Outputs S_u, S_p of k, of ε and of ζ, and S_u (the right-hand
    /// side) and S_p of f (None: ν_t, T, L², P_k only).
    pub sources: Option<[Buf; 8]>,
    /// Output ν_t/σ_k.
    pub dk: Buf,
    /// Output ν_t/σ_ε.
    pub de: Buf,
    /// Output ν_t/σ_ζ.
    pub dz: Buf,
    /// Output: 1 where the ν_t cap acted.
    pub capped: Buf,
    /// The elliptic-blending model: ∇k, ∇φ, Σ (∂_l ∂_j U_i)² and
    /// ∇·(ν_t/σ_k ∇k) per cell.
    pub eb: Option<(Buf3, Buf3, Buf, Buf)>,
}

/// Inputs and outputs of [`FlowDevice::turb_faces`].
#[derive(Debug, Clone, Copy)]
pub struct TurbFaceOps {
    /// ν_t per cell.
    pub nu_t: Buf,
    /// σ_k ν_t per cell.
    pub dk: Buf,
    /// σ_ω ν_t per cell.
    pub dw: Buf,
    /// Laminar viscosity per face.
    pub mu: Buf,
    /// 1 on the wall faces (per boundary face).
    pub wall: Buf,
    /// Output μ + ρ ν_t per face (the device's face viscosity).
    pub mu_eff: Buf,
    /// Output ρ(ν + σ_k ν_t) per face.
    pub gamma_k: Buf,
    /// Output ρ(ν + σ_ω ν_t) per face.
    pub gamma_w: Buf,
    /// k-ε-ζ-f: ν_t/σ_ζ and L² per cell, and the outputs ρ(ν + ν_t/σ_ζ)
    /// and ρ L² per face.
    pub zeta_f: Option<[Buf; 4]>,
    /// Factor of the molecular viscosity in the diffusivities of k, ε/ω
    /// and ζ: 1; ½ for the elliptic-blending model.
    pub mol: f64,
}

/// Inputs and outputs of [`FlowDevice::scalar_turb_faces`].
#[derive(Debug, Clone, Copy)]
pub struct ScalarTurbOps {
    /// Laminar viscosity per face.
    pub mu: Buf,
    /// The face viscosity of the flow (μ + ρ ν_t).
    pub mu_eff: Buf,
    /// ρ Γ per face (molecular).
    pub gamma0: Buf,
    /// 1 on the wall faces (per boundary face).
    pub wall: Buf,
    /// q d_n per boundary face.
    pub fq: Buf,
    /// 0 on fixed-value faces, 1 otherwise (per boundary face).
    pub bc_a: Buf,
    /// Model of Pr_t (`turb::thermal::PRT_*`) and its constant.
    pub prt: (u32, f64),
    /// Output: the face diffusivity.
    pub gamma: Buf,
    /// Output: b of the flux faces (per boundary face).
    pub bc_b: Buf,
}

/// Inputs and outputs of [`FlowDevice::scalar_wall_faces`].
#[derive(Debug, Clone, Copy)]
pub struct ScalarWallOps {
    /// y⁺ per wall face.
    pub y_plus: Buf,
    /// Laminar viscosity per face.
    pub mu: Buf,
    /// ρ Γ per face (molecular).
    pub gamma0: Buf,
    /// q d_n per boundary face.
    pub fq: Buf,
    /// 0 on fixed-value faces, 1 otherwise (per boundary face).
    pub bc_a: Buf,
    /// Thermal wall law (`turb::thermal::KADER`, …).
    pub law: u32,
    /// Pr_t of the log layer.
    pub prt: f64,
    /// T⁺ is taken at ỹ⁺ = max(y⁺, `y_floor`) (11.06 with the scalable
    /// momentum treatment, else 0).
    pub y_floor: f64,
    /// Output: the wall-face diffusivity (per face).
    pub gamma: Buf,
    /// Output: b of the flux faces (per boundary face).
    pub bc_b: Buf,
    /// Output T⁺ per wall face, if wanted.
    pub t_plus: Option<Buf>,
}

/// Inputs and outputs of [`FlowDevice::wall_faces`].
#[derive(Debug, Clone, Copy)]
pub struct WallFaceOps {
    /// Unit normals per wall face.
    pub normal: Buf3,
    /// y_P per wall face.
    pub y: Buf,
    /// k_s per wall face.
    pub ks: Buf,
    /// Cell velocity.
    pub u: Buf3,
    /// Boundary velocity.
    pub ub: Buf3,
    /// ν.
    pub nu: f64,
    /// ν per cell (D194: several fluids); `None`: `nu` everywhere.
    pub nu_cell: Option<Buf>,
    /// Treatment.
    pub treatment: u32,
    /// Constants.
    pub c: foehncraft_kernels::turb::sst::SstConstants,
    /// Constants of the k-ε family.
    pub ke: foehncraft_kernels::turb::keps::KepsConstants,
    /// k per cell.
    pub k: Buf,
    /// C_s per wall face.
    pub cs: Buf,
    /// u_τ per wall face (input of the last iteration and output).
    pub u_tau: Buf,
    /// Output y⁺ per wall face.
    pub y_plus: Buf,
    /// Output per wall face: the y⁺ of the thermal law (y* = C_μ^¼ k_P^½
    /// y_P/ν of the k-ε family and ζ-f/EB, y⁺ of SST).
    pub y_star: Buf,
    /// Output ω per wall face.
    pub omega: Buf,
    /// Output viscous weight per wall face.
    pub blending: Buf,
    /// Output: the viscous weight of the wall treatment per wall face for
    /// the output (the wall field `blending`; not used by the solver).
    pub viscous_weight: Buf,
    /// Output: the weighted wall production per wall face.
    pub pk: Buf,
    /// Output: the weight of the wall production per wall face.
    pub prod_w: Buf,
    /// Output μ_w into the face viscosity.
    pub mu_eff: Buf,
    /// Output wall values of k and ω (bc_b, per boundary face).
    pub k_b: Buf,
    /// See `k_b`.
    pub w_b: Buf,
    /// k-ε-ζ-f: ζ per cell and the output wall value of f per boundary
    /// face.
    pub zeta_f: Option<(Buf, Buf)>,
    /// N_ζ of the f wall value (1: ζ-f, 2: ζ-f₀).
    pub n_zeta: f64,
    /// The compound treatment of 2023: the velocity gradient and the
    /// pressure gradient per cell (C_U of the non-equilibrium factor).
    pub non_eq: Option<([Buf3; 3], Buf3)>,
}

/// Inputs and outputs of [`FlowDevice::wall_cells`].
#[derive(Debug, Clone, Copy)]
pub struct WallCellOps {
    /// ω per wall face.
    pub omega: Buf,
    /// Area weight per wall face (0 on rough faces).
    pub weight: Buf,
    /// Viscous weight of the blending per wall face.
    pub blending: Buf,
    /// Weighted wall production per wall face.
    pub pk: Buf,
    /// Weight of the wall production per wall face.
    pub prod_w: Buf,
    /// Output: the area mean of `pk` per cell.
    pub wall_pk: Buf,
    /// Output: the area mean of `prod_w` per cell.
    pub prod_mask: Buf,
    /// Output: the fixed value of ω per cell.
    pub fix_val: Buf,
    /// Output: its weight per cell (0 off the near-wall cells).
    pub fix_mask: Buf,
    /// The fixed variable, moved to the fixed value by the weight in the
    /// near-wall cells (the k-ε family), or none.
    pub phi: Option<Buf>,
}

/// Static data of a device.
pub struct Setup<'a> {
    /// Discretisation (geometry, velocity weights, matrix structure).
    pub disc: &'a Discretisation,
    /// Least-squares weights of the pressure.
    pub lsq_p: &'a Soa3,
    /// Boundary arrays.
    pub bc: &'a BoundaryArrays,
    /// Density.
    pub rho: f64,
    /// Dynamic viscosity (when it is the same on every face).
    pub mu: f64,
    /// Viscosity per face (D94), or `None` for `mu` everywhere.
    pub mu_face: Option<&'a [f64]>,
    /// ε² of the Venkatakrishnan limiter per cell (empty: unused).
    pub eps2: &'a [f64],
    /// CPU threads.
    pub threads: Option<usize>,
}

/// The operations of the flow algorithms. Buffers are created by
/// [`alloc`](FlowDevice::alloc); an output buffer is never one of the
/// inputs of the same operation unless the operation says so.
pub trait FlowDevice {
    /// Backend name (`cpu`, `cuda`).
    fn name(&self) -> String;
    /// Device description.
    fn device(&self) -> String;

    /// A new buffer of `n` zeros.
    fn alloc(&mut self, n: usize) -> Result<Buf>;
    /// Copies host values into a buffer.
    fn upload(&mut self, b: Buf, v: &[f64]) -> Result<()>;
    /// Copies a buffer to the host.
    fn download(&mut self, b: Buf) -> Result<Vec<f64>>;
    /// dst = src.
    fn copy(&mut self, src: Buf, dst: Buf) -> Result<()>;
    /// b = v everywhere.
    fn fill(&mut self, b: Buf, v: f64) -> Result<()>;
    /// out = a (x − y); `out` may be `x` or `y`.
    fn scaled_difference(&mut self, out: Buf, a: f64, x: Buf, y: Buf) -> Result<()>;
    /// out = min(dc, 1/(a/ds + b)) per cell: D̂ of the coupled
    /// preconditioner ([`foehncraft_kernels::flow::schur_diagonal`]).
    fn schur_diagonal(&mut self, out: Buf, a: f64, b: f64, ds: Buf, dc: Buf) -> Result<()>;
    /// out = ρ V max(inv_dt, factor · num/den) per cell: the pseudo-time
    /// coefficient of a scalar row with the step min(Δτ_max, T/factor)
    /// ([`foehncraft_kernels::scalar::pseudo_coefficient`]; `num`, `den`
    /// absent: 1).
    fn pseudo_coefficient(
        &mut self,
        out: Buf,
        inv_dt: f64,
        factor: f64,
        num: Option<Buf>,
        den: Option<Buf>,
    ) -> Result<()>;
    /// out = a x + b y; `out` may be `x` or `y`.
    fn lincomb(&mut self, out: Buf, a: f64, x: Buf, b: f64, y: Buf) -> Result<()>;
    /// out += v_x² + v_y² + v_z² per entry.
    fn squares_add(&mut self, v: Buf3, out: Buf) -> Result<()>;
    /// out = s v per entry (`out` may be `v`).
    fn scale_vector(&mut self, s: Buf, v: Buf3, out: Buf3) -> Result<()>;
    /// Σ x² in the fixed block order.
    fn sum_sq(&mut self, x: Buf) -> Result<f64>;
    /// Σ x y in the fixed block order (x and y of equal length).
    fn dot(&mut self, x: Buf, y: Buf) -> Result<f64>;
    /// Minimum and maximum.
    fn min_max(&mut self, x: Buf) -> Result<(f64, f64)>;
    /// One value.
    fn value(&mut self, b: Buf, i: usize) -> Result<f64>;
    /// The values at the indices `idx` (a few, such as the cells next to
    /// a periodic plane).
    fn gather(&mut self, b: Buf, idx: &[u32]) -> Result<Vec<f64>> {
        let v = self.download(b)?;
        Ok(idx.iter().map(|&i| v[i as usize]).collect())
    }
    /// b −= s.
    fn subtract(&mut self, b: Buf, s: f64) -> Result<()>;
    /// `b[i] *= s`.
    fn scale_entry(&mut self, b: Buf, i: usize, s: f64) -> Result<()>;

    /// New boundary velocity and pressure from the current fields.
    fn boundary_values(&mut self, a: &BoundaryArgs) -> Result<()>;
    /// The per-boundary-face buffers of the convective outlets that the
    /// boundary values read ([`foehncraft_kernels::flow::UB_CONVECTIVE`]):
    /// the old-level part c (three components) and β, zero at creation.
    fn convective_buffers(&self) -> (Buf3, Buf);
    /// Backflow suppression after the predicted fluxes `m`: closes the
    /// pressure-outlet faces whose predicted flux enters the domain (their
    /// flux in `m` becomes 0, the closed state is kept on the device for
    /// the pressure correction of this iteration).
    fn close_backflow(&mut self, m: Buf) -> Result<()>;
    /// New prescribed pressures (or total pressures) of the boundary faces
    /// (time-dependent boundary values, D53).
    fn set_pressure_values(&mut self, v: &[f64]) -> Result<()>;
    /// Least-squares gradient of `phi` with boundary values `value`, where
    /// `zero_gradient` (one of [`Lsq`]'s masks, by `mask`) is 1 the owner
    /// value is used instead.
    fn gradient(&mut self, w: Lsq, phi: Buf, mask: Mask, value: Buf, out: Buf3) -> Result<()>;
    /// Venkatakrishnan limiter of the velocity component `phi` with its
    /// gradient and boundary values.
    fn limiter(&mut self, phi: Buf, value: Buf, grad: Buf3, eps2: Eps2, out: Buf) -> Result<()>;
    /// Least-squares gradient of a scalar with its own weights `w` (one
    /// entry per cell-face entry) and boundary coefficients φ_f = a φ_P + b.
    /// `jump`: the jump per periodic face (D109 (1); `None`: plain
    /// periodic).
    #[allow(clippy::too_many_arguments)]
    fn gradient_with(
        &mut self,
        w: Buf3,
        phi: Buf,
        bc_a: Buf,
        bc_b: Buf,
        jump: Option<Buf>,
        out: Buf3,
    ) -> Result<()>;
    /// Venkatakrishnan limiter of a scalar with its boundary coefficients.
    #[allow(clippy::too_many_arguments)]
    fn limiter_with(
        &mut self,
        phi: Buf,
        bc_a: Buf,
        bc_b: Buf,
        jump: Option<Buf>,
        grad: Buf3,
        eps2: Eps2,
        out: Buf,
    ) -> Result<()>;
    /// Face coefficients of a scalar.
    fn scalar_faces(&mut self, a: &ScalarFaceArgs) -> Result<()>;
    /// Rows of a scalar.
    fn scalar_rows(&mut self, a: &ScalarRowArgs) -> Result<()>;
    /// Residual scale ‖a φ‖ + ‖b − (a/α − a) φ‖ of a scalar row (unrelaxed
    /// a and b, without the pseudo-time coefficient `pseudo` of the row).
    fn scalar_scale(
        &mut self,
        phi: Buf,
        rhs: Buf,
        diag: Buf,
        alpha: f64,
        pseudo: Option<Buf>,
    ) -> Result<f64>;
    /// Face mass fluxes (Rhie–Chow).
    fn mass_flux(&mut self, a: &FluxArgs, out: Buf) -> Result<()>;
    /// Momentum face coefficients.
    fn momentum_faces(&mut self, a: &MomentumFaceArgs) -> Result<()>;
    /// Momentum rows.
    fn momentum_rows(&mut self, a: &MomentumRowArgs) -> Result<()>;
    /// Spatial operator of the momentum equations with the coefficients
    /// `cp`, `cn`, `bf` and the velocity `u` (Crank–Nicolson).
    fn momentum_operator(&mut self, cp: Buf, cn: Buf, bf: Buf3, u: Buf3, out: Buf3) -> Result<()>;
    /// Momentum residual scale ‖a u‖ + ‖b − (a/α − a) u‖ over the three
    /// components (unrelaxed a and b).
    fn residual_scale(&mut self, u: Buf3, rhs: Buf3, diag: Buf3, alpha: f64) -> Result<f64>;
    /// Slot 0 of every row of `val` = `diag`.
    fn set_diagonal(&mut self, val: Buf, diag: Buf) -> Result<()>;
    /// Net outflow Σ ±ṁ per cell.
    fn net_outflow(&mut self, m: Buf, out: Buf) -> Result<()>;
    /// Throughflow Σ |ṁ| per cell.
    fn throughflow(&mut self, m: Buf, out: Buf) -> Result<()>;
    /// Courant number per cell, Δt Σ|ṁ| / (2 ρ V).
    fn courant(&mut self, m: Buf, dt: f64, out: Buf) -> Result<()>;
    /// Pressure-correction face coefficients from D^C.
    fn pressure_faces(&mut self, dc: Buf, out: Buf) -> Result<()>;
    /// Pressure-correction rows.
    fn pressure_rows(&mut self, cf: Buf, m: Buf, val: Buf, rhs: Buf) -> Result<()>;
    /// Corrected fluxes ṁ − c_f (p'_N − p'_P).
    fn correct_fluxes(&mut self, cf: Buf, m: Buf, pc: Buf, out: Buf) -> Result<()>;
    /// u_i −= D^C ∇p'_i.
    fn correct_velocity(&mut self, dc: Buf, grad: Buf3, u: Buf3) -> Result<()>;
    /// out = −Σ_N val_PN (x_N − own x_P) / diag: the neighbour correction
    /// H(x)/A of PISO after Issa's corrector (`own` = 0), or after a
    /// SIMPLEC-type corrector, whose D^C already holds Σ_N val_PN x_P
    /// (`own` = 1).
    fn neighbour_correction(
        &mut self,
        val: Buf,
        x: Buf,
        diag: Buf,
        out: Buf,
        own: f64,
    ) -> Result<()>;
    /// Mass flux of a cell velocity increment.
    fn delta_flux(&mut self, delta: Buf3, out: Buf) -> Result<()>;
    /// out = Aᵀ in the storage of A = `val` (D59: the same solvers and the
    /// multigrid then work on Aᵀ).
    fn transpose(&mut self, val: Buf, out: Buf) -> Result<()>;
    /// out = A x for the matrix `val`.
    fn matrix_product(&mut self, val: Buf, x: Buf, out: Buf) -> Result<()>;
    /// out = Aᵀ x for the matrix `val`, without storing Aᵀ (D59).
    fn transposed_product(&mut self, val: Buf, x: Buf, out: Buf) -> Result<()>;

    /// Solves `val` x = `rhs` (x holds the start value): CG with `pc`
    /// when `cg`, else BiCGStab. Returns the iterations and ‖b − A x₀‖.
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
    /// Multigrid setup on the matrix `val` (the pressure: slot 0, which
    /// must be the active slot).
    fn amg_setup(&mut self, val: Buf, s: &AmgSettings) -> Result<()>;
    /// Multigrid setup on the matrix `val` in slot `slot`, which becomes the
    /// active one (the turbulence equations, stage 6 PR 6d: slot 1, so that
    /// the pressure hierarchy of slot 0 is kept); with `symmetric_part` the
    /// hierarchy is built from ½(A + Aᵀ) (a transport matrix with
    /// convection, preconditioning BiCGStab).
    fn amg_setup_slot(
        &mut self,
        slot: usize,
        val: Buf,
        s: &AmgSettings,
        symmetric_part: bool,
    ) -> Result<()>;
    /// Makes multigrid slot `slot` the active one (back to 0 after a
    /// turbulence solve).
    fn amg_select(&mut self, slot: usize) -> Result<()>;
    /// Setup of the reduction multigrid (D179) on the matrix `val`, which
    /// need not be symmetric (a momentum matrix: slots 5 to 7), in slot
    /// `slot`, which becomes the active one.
    fn reduction_setup_slot(&mut self, slot: usize, val: Buf, s: &ReductionSettings) -> Result<()>;
    /// x = the result of `cycles` multigrid cycles of the active slot on
    /// `val` x = `rhs` from x = 0: a stationary iteration, so a fixed
    /// linear operator on `rhs` (the sweep of a preconditioner).
    fn multigrid_cycles(&mut self, val: Buf, rhs: Buf, x: Buf, cycles: usize) -> Result<()>;
    /// The face viscosity of the momentum equation (a face buffer that a
    /// turbulence model overwrites with μ + ρ ν_t, stage 6).
    fn face_viscosity(&self) -> Buf;
    /// Replaces the density per cell (`n_cells` values) and per face
    /// (`n_faces` values) that the flux, momentum, scalar and turbulence
    /// operations read (D150). The device starts with the constant
    /// density of its setup on every cell and face.
    fn set_density(&mut self, cell: &[f64], face: &[f64]) -> Result<()>;
    /// Sets the balanced-force arrays of gravity (D154, stage 7b1): the
    /// offset q_f subtracted from the Rhie–Chow pressure jump of every face
    /// (`n_faces` values) and the term e_b added to the extrapolated
    /// boundary pressure (`n_boundary` values). Without a call both are
    /// absent and the operations are those without gravity, bit for bit.
    fn set_balanced_force(&mut self, jump: &[f64], pb: &[f64]) -> Result<()>;
    /// The variable-density cell arrays of stage 7c (D151): the mass
    /// source V ∂ρ/∂t of the continuity rows and the densities ρⁿ, ρⁿ⁻¹ of
    /// the momentum time term; `None` clears an array (no source, the
    /// current density for the old levels).
    fn set_variable_density(
        &mut self,
        source: Option<&[f64]>,
        old: Option<&[f64]>,
        old2: Option<&[f64]>,
    ) -> Result<()>;
    /// Uploads the static arrays of the wall treatment.
    fn wall_setup(&mut self, s: &WallStatic) -> Result<()>;
    /// The k-ω SST model per cell.
    fn sst_cells(&mut self, a: &SstCellOps) -> Result<()>;
    /// S_comp (xx, yy, zz, xy, yz, zx = 0 … 5) of the velocity gradient.
    fn strain_component(&mut self, grad_u: [Buf3; 3], comp: u32, out: Buf) -> Result<()>;
    /// out = (u − v_frame)·grad + (s − s_old)/Δt (`old` = (s, s_old, 1/Δt)).
    fn convective_derivative(
        &mut self,
        u: Buf3,
        v_frame: Option<Buf3>,
        grad: Buf3,
        old: Option<(Buf, Buf, f64)>,
        out: Buf,
    ) -> Result<()>;
    /// f_r1 of the curvature correction of SST per cell (D170).
    fn curvature_cells(&mut self, a: &CurvatureOps) -> Result<()>;
    /// The non-linear stress of the quadratic constitutive relation of SST
    /// per cell (D170; `mode` QCR_2000 or QCR_2013, μ_t = ρ ν_t).
    #[allow(clippy::too_many_arguments)]
    fn qcr_cells(
        &mut self,
        grad_u: [Buf3; 3],
        nu_t: Buf,
        c: &foehncraft_kernels::turb::sst::SstConstants,
        mode: u32,
        tau: [Buf; 6],
    ) -> Result<()>;
    /// b_f −= τ_f·S: an explicit cell stress `tau` (xx, yy, zz, xy, yz, zx)
    /// in the momentum face parts `bf`; with the cell gradients `gtau` of
    /// the six components, skewness-corrected face stresses (D236).
    fn explicit_stress_faces(
        &mut self,
        tau: [Buf; 6],
        gtau: Option<[Buf3; 6]>,
        bf: Buf3,
    ) -> Result<()>;
    /// Uploads the stencil of the quadratic least-squares gradient (D236):
    /// the treated cells, the start of each stencil and its points (a cell,
    /// or the number of cells plus a boundary face index); returns the
    /// slot of the stencil.
    fn quadratic_setup(&mut self, cells: &[u32], ptr: &[u32], idx: &[u32]) -> Result<usize>;
    /// Overwrites the velocity gradient `out` at the treated cells of the
    /// stencil `slot` of [`quadratic_setup`](Self::quadratic_setup) with
    /// Σ_k w_k (u_k − u_c), components 0 … `dims` − 1.
    #[allow(clippy::too_many_arguments)]
    fn quadratic_gradient(
        &mut self,
        slot: usize,
        w: Buf3,
        dims: u32,
        u: Buf3,
        ub: Buf3,
        out: [Buf3; 3],
    ) -> Result<()>;
    /// The standard or realizable k-ε model per cell.
    fn keps_cells(&mut self, a: &KepsCellOps) -> Result<()>;
    /// The k-ε-ζ-f model per cell.
    fn zeta_f_cells(&mut self, a: &ZetaFCellOps) -> Result<()>;
    /// Face viscosities of a turbulent flow.
    fn turb_faces(&mut self, a: &TurbFaceOps) -> Result<()>;
    /// The turbulent diffusivity of a passive scalar per face (D135).
    fn scalar_turb_faces(&mut self, a: &ScalarTurbOps) -> Result<()>;
    /// The thermal wall function per wall face (D135, after
    /// [`scalar_turb_faces`](Self::scalar_turb_faces)).
    fn scalar_wall_faces(&mut self, a: &ScalarWallOps) -> Result<()>;
    /// The wall treatment per wall face (after [`wall_setup`](Self::wall_setup)).
    fn wall_faces(&mut self, a: &WallFaceOps) -> Result<()>;
    /// The area-weighted wall-face values `omega` of the near-wall cells
    /// into `fix_val`.
    fn wall_cells(&mut self, a: &WallCellOps) -> Result<()>;
    /// x ← min(max(x, lo, rel · old), hi) with `rel_old` = (old, rel)
    /// (else without the relative bound); returns the number of entries a
    /// bound changed (`work` is a scratch buffer of the length of x).
    fn floor(
        &mut self,
        x: Buf,
        lo: f64,
        hi: f64,
        rel_old: Option<(Buf, f64)>,
        work: Buf,
    ) -> Result<usize>;
    /// Waits for queued work (timers).
    fn sync(&mut self) -> Result<()>;
}

/// Boundary masks of the gradients: where the owner value replaces the face
/// value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mask {
    /// Zero normal gradient of the velocity (outlets).
    Velocity,
    /// Zero normal gradient of the pressure (walls, symmetry planes).
    Pressure,
    /// The pressure correction: zero on pressure boundaries, zero normal
    /// gradient elsewhere.
    PressureCorrection,
}

/// Builds a GPU flow device from the static data and the kernel module
/// (for CUDA the PTX file of `gpu/ptx-kernels/build.py`).
pub type GpuFactory = fn(&Setup, &std::path::Path) -> Result<Box<dyn FlowDevice>>;

static GPU: std::sync::OnceLock<(&'static str, GpuFactory)> = std::sync::OnceLock::new();

/// Registers the GPU flow device of this build under `name` (`cuda`);
/// called once by the front end (`foehncraft_cuda::register`, D60: no
/// vendor code in this crate).
pub fn register_gpu(name: &'static str, f: GpuFactory) {
    let _ = GPU.set((name, f));
}

/// The flow backends of this build: `cpu`, and the GPU backend when one is
/// registered.
pub fn available() -> Vec<&'static str> {
    let mut v = vec!["cpu"];
    if let Some((n, _)) = GPU.get() {
        v.push(n);
    }
    v
}

/// Creates a device: `cpu`, or `cuda` when the build registered it (`ptx`
/// names its kernel module, from `gpu/ptx-kernels`).
pub fn create(name: &str, s: &Setup, ptx: &std::path::Path) -> Result<Box<dyn FlowDevice>> {
    match name {
        "cpu" => Ok(Box::new(cpu::CpuFlow::cpu(s)?)),
        "cuda" => match GPU.get() {
            Some((n, f)) if *n == name => f(s, ptx),
            _ => anyhow::bail!(
                "this build has no GPU backend; build with `--features cuda` \
                 (cargo build --release -p foehncraft-cli --features cuda)"
            ),
        },
        other => anyhow::bail!("unknown backend `{other}` (cpu, cuda)"),
    }
}
