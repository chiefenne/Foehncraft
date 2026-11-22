//! Kernel bodies of the incompressible flow solver (stage 4): momentum face
//! coefficients, the momentum and pressure-correction rows, the Rhie–Chow
//! face mass flux and the boundary values, for the SIMPLEC algorithm of
//! `foehncraft-flow` (Van Doormaal & Raithby 1984) on the collocated,
//! cell-centred grid (D19, D23).
//!
//! Like [`crate::fv`], every body computes the outputs of one face, one
//! boundary face or one cell from its inputs (a gather: no atomics, faces in
//! stored order), so a CPU loop and one GPU thread per index give bitwise
//! the same results (D31).
//!
//! **Boundary faces.** Each boundary face carries
//!
//! - a momentum kind ([`MOM_FIXED`]: the face velocity `ub` is known, walls
//!   and inlets; [`MOM_OUTFLOW`]: the face velocity is extrapolated from the
//!   owner, outlets and openings; [`MOM_SYMMETRY`]: no shear stress, the
//!   normal component vanishes);
//! - a flux kind ([`FLUX_FIXED`]: ṁ_f = ρ u_b·S; [`FLUX_PRESSURE`]: the
//!   face pressure is known and ṁ_f follows from the Rhie–Chow formula);
//! - the boundary values of the velocity `ub` and the pressure `pb`, updated
//!   once per outer iteration by [`boundary_velocity`] and
//!   [`boundary_pressure`].
//!
//! The gradients of all variables use the least-squares body
//! [`crate::fv::lsq_gradient`] with `a = 0` and `b` the boundary value at the
//! face centroid, so a linear field is reproduced exactly at every face.

use crate::fv::{SellLayout, Topology, Vec3};
use crate::{ld, unpack_face};
use foehncraft_core::Real;

/// Momentum kind of a boundary face: the face velocity is prescribed.
pub const MOM_FIXED: u32 = 0;
/// Momentum kind of a boundary face: the face velocity is extrapolated.
pub const MOM_OUTFLOW: u32 = 1;
/// Momentum kind of a boundary face: symmetry plane. The diffusive flux is
/// μ D (u_P·n̂) n̂ (no shear stress, zero normal velocity on the face); its
/// part μ D n̂_i² u_P,i is implicit in the diagonal of component i (so the
/// three components get different diagonals), the cross-component part
/// explicit. Treating the whole flux with the lagged face value instead
/// would act as a strong extra under-relaxation of the tangential
/// components.
pub const MOM_SYMMETRY: u32 = 2;
/// Flux kind of a boundary face: ṁ_f = ρ u_b·S.
pub const FLUX_FIXED: u32 = 0;
/// Flux kind of a boundary face: pressure prescribed, Rhie–Chow flux.
pub const FLUX_PRESSURE: u32 = 1;
/// Flux kind of a boundary face: as [`FLUX_PRESSURE`], but a face whose
/// predicted flux enters the domain is closed for that outer iteration
/// (ṁ_f = 0 and no pressure-correction coefficient; backflow suppression
/// at pressure outlets, see [`close_backflow`]).
pub const FLUX_PRESSURE_NO_BACKFLOW: u32 = 2;

/// Update rule of the boundary velocity: keep the prescribed value.
pub const UB_FIXED: u32 = 0;
/// Update rule: symmetry, u_b = u_P − (u_P·n̂) n̂.
pub const UB_SYMMETRY: u32 = 1;
/// Update rule: zero normal gradient with linear extrapolation along the
/// face, u_b = u_P + ∇u_P·r_t, r_t the part of x_f − x_P tangential to the
/// face (outlets).
pub const UB_EXTRAPOLATE: u32 = 2;
/// Update rule: normal to the face from the face mass flux,
/// u_b = ṁ_f n̂ / (ρ |S|) (pressure inlets).
pub const UB_FROM_FLUX: u32 = 3;
/// Update rule of a pressure outlet with backflow treatment: the
/// tangential extrapolation of [`UB_EXTRAPOLATE`] where the flux leaves
/// the domain (ṁ_f ≥ 0), its normal part (u_e·n̂) n̂ where it enters (the
/// backflow enters normal to the face).
pub const UB_OUTLET: u32 = 4;
/// Update rule of a convective outlet (transient runs): the discrete
/// convection equation ∂u/∂t + U_c ∂u/∂n = 0 between the owner and the
/// face, u_b = (c + β u_e)/(1 + β), with u_e the owner velocity
/// extrapolated along the face (as [`UB_EXTRAPOLATE`]), c the old-level
/// part of the backward difference over a₀ and β = U_c Δt /(a₀ Δn) per
/// face (both set per time step by the solver; β = 0 keeps u_b = c).
pub const UB_CONVECTIVE: u32 = 5;
/// Update rule of a pressure outlet that suppresses backflow: the
/// extrapolation of [`UB_EXTRAPOLATE`] where the flux leaves the domain
/// (ṁ_f > 0), its tangential part u_e − (u_e·n̂) n̂ where the face is closed
/// (ṁ_f = 0 after [`close_backflow`]): the closed face is a slip wall.
pub const UB_NO_BACKFLOW: u32 = 6;

/// Update rule of the boundary pressure: keep the prescribed value.
pub const PB_FIXED: u32 = 0;
/// Update rule: linear extrapolation p_b = p_P + ∇p_P·(x_f − x_P).
pub const PB_EXTRAPOLATE: u32 = 1;
/// Update rule: total pressure, p_b = p₀ − ½ ρ (ṁ_f / (ρ |S|))².
pub const PB_TOTAL: u32 = 2;
/// Update rule of a static-pressure outlet with backflow treatment: the
/// prescribed static pressure where the flux leaves the domain, and where
/// it enters the prescribed value taken as the total pressure of the
/// entering fluid, p_b = p_set − ½ ρ (ṁ_f / (ρ |S|))².
pub const PB_OUTLET: u32 = 3;
/// Update rule of a static-pressure outlet that suppresses backflow: the
/// prescribed value where the flux leaves the domain (ṁ_f > 0), the
/// extrapolation of [`PB_EXTRAPOLATE`] where the face is closed (ṁ_f = 0).
/// The rules read only the mass flux, so a restart or a rejected step
/// needs no state of the closed faces.
pub const PB_NO_BACKFLOW: u32 = 4;
/// As [`PB_NO_BACKFLOW`] for a total-pressure outlet: [`PB_TOTAL`] where
/// the flux leaves the domain, the extrapolation where the face is closed.
pub const PB_TOTAL_NO_BACKFLOW: u32 = 5;

#[inline(always)]
fn dot<T: Real>(a: [T; 3], b: [T; 3]) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
fn sub<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Gradients of the three velocity components.
#[derive(Clone, Copy)]
pub struct VelocityGradient<'a, T> {
    /// ∇u_x.
    pub x: Vec3<'a, T>,
    /// ∇u_y.
    pub y: Vec3<'a, T>,
    /// ∇u_z.
    pub z: Vec3<'a, T>,
}

impl<T: Real> VelocityGradient<'_, T> {
    /// ∇u_x, ∇u_y, ∇u_z of cell `c`.
    ///
    /// # Safety
    /// `c` is inside all nine slices.
    #[inline(always)]
    pub unsafe fn at(&self, c: usize) -> [[T; 3]; 3] {
        // SAFETY: guaranteed by the caller.
        unsafe { [self.x.at(c), self.y.at(c), self.z.at(c)] }
    }
}

/// Geometry of the faces and cells (see `foehncraft_fv::Discretisation`).
#[derive(Clone, Copy)]
pub struct FaceGeometry<'a, T> {
    /// Face area vectors S.
    pub area: Vec3<'a, T>,
    /// Face centroids.
    pub face_centroid: Vec3<'a, T>,
    /// Cell centroids.
    pub cell_centroid: Vec3<'a, T>,
    /// |S|²/(d·S), d from owner to neighbour or to the boundary face centroid.
    pub dgeo: &'a [T],
    /// k = S − d |S|²/(d·S).
    pub k: Vec3<'a, T>,
    /// Owner weight of the face interpolation.
    pub weight: &'a [T],
    /// Shift of the neighbour centroid per periodic face
    /// ([`crate::fv::neighbour_centroid`]).
    pub shift: Vec3<'a, T>,
}

/// Inputs of [`momentum_face`].
#[derive(Clone, Copy)]
pub struct MomentumInputs<'a, T> {
    /// Connectivity.
    pub topo: Topology<'a>,
    /// Geometry.
    pub geo: FaceGeometry<'a, T>,
    /// Face mass flux ṁ_f (kg/s), out of the owner.
    pub mass_flux: &'a [T],
    /// Dynamic viscosity on the face (Pa s).
    pub mu: &'a [T],
    /// Cell velocity.
    pub u: Vec3<'a, T>,
    /// Cell velocity of the central correction: `u`, or in the first outer
    /// iteration of a time step the velocity extrapolated to tⁿ⁺¹ (as the
    /// gradients of the linear-upwind correction are).
    pub uc: Vec3<'a, T>,
    /// Cell velocity gradient.
    pub grad: VelocityGradient<'a, T>,
    /// Limiter of each velocity component per cell (1: none).
    pub psi: Vec3<'a, T>,
    /// Boundary face velocity (one entry per boundary face).
    pub ub: Vec3<'a, T>,
    /// Momentum kind per boundary face ([`MOM_FIXED`], [`MOM_OUTFLOW`]).
    pub mom_kind: &'a [u32],
    /// 1 for linear upwind, 0 for first-order upwind.
    pub blend: T,
    /// Fraction of the bounded central scheme on interior faces (0: none;
    /// the rest takes the linear-upwind correction of `blend`).
    pub central: T,
    /// β of the normalised-variable bounding of the central value (0:
    /// plain central), see [`crate::fv::central_weight`].
    pub beta: T,
    /// 1 with the non-orthogonal diffusion correction, 0 without.
    pub nonorth: T,
    /// 1 with convection, 0 for Stokes flow.
    pub convection: T,
    /// 1 with the transposed part ∇·[μ (∇u)ᵀ] of the viscous stress as a
    /// deferred term (variable eddy viscosity, stage 6 specification 2.4),
    /// 0 without (laminar flow: it vanishes for constant μ).
    pub transposed: T,
    /// 1 to take the transposed part without its trace,
    /// ∇·[μ ((∇u)ᵀ − (∇·u) I)] (the cell-quality remediation, D50): the
    /// same term where ∇·u = 0, without the grad-div part μ ∇(∇·u) that
    /// cell gradients on tetrahedra leave in it; 0 for the full part.
    pub trace_free: T,
}

/// Momentum flux out of the owner through face `f` for the three velocity
/// components, F_f,i = c_P u_P,i + c_N u_N,i + b_f,i (the same c_P, c_N for
/// all components). Returns `(c_P, c_N, [b_x, b_y, b_z])`.
///
/// - Interior faces: first-order upwind convection and the orthogonal part of
///   the diffusion implicit; the linear-upwind correction
///   ṁ ψ_U ∇u_U·(x_f − x_U) (or, for the fraction `central`, the bounded
///   central correction ṁ γ (u_f^C − u_U) per component, with u_f^C
///   interpolated as in [`face_mass_flux`]) and the non-orthogonal part
///   −μ k·∇u_f deferred (as [`crate::fv::face_coefficients`]).
/// - [`MOM_FIXED`] faces: convection ṁ⁺ u_P (with the linear-upwind
///   correction) + ṁ⁻ u_b, diffusion −μ [D (u_b − u_P) + k·∇u_P].
/// - [`MOM_OUTFLOW`] faces: convection ṁ⁺ u_P (with the correction) + ṁ⁻ u_b
///   (back flow takes the boundary value), no diffusive flux (zero normal
///   gradient, a developed outflow).
/// - [`MOM_SYMMETRY`] faces (no mass flux): c_P = μ D, which
///   [`momentum_row`] puts into the diagonal of component i weighted with
///   n̂_i², and b_i = μ D n̂_i Σ_{j≠i} n̂_j u_P,j; with `transposed` 2μ D
///   in both (the normal stress 2μ ∂u_n/∂n of the full stress).
///
/// # Safety
/// `f` is a face of the topology and every array of `inp` is long enough
/// for the indices it is read at.
#[inline(always)]
pub unsafe fn momentum_face<T: Real>(f: usize, inp: &MomentumInputs<T>) -> (T, T, [T; 3]) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let t = &inp.topo;
        let g = &inp.geo;
        let p = ld(t.owner, f) as usize;
        let m = inp.convection * ld(inp.mass_flux, f);
        let mu = ld(inp.mu, f);
        let dd = mu * ld(g.dgeo, f);
        let mp = m.max_of(T::ZERO);
        let mm = m.min_of(T::ZERO);
        let xf = g.face_centroid.at(f);
        let kf = g.k.at(f);
        let gp = inp.grad.at(p);
        let rp = sub(xf, g.cell_centroid.at(p));
        let psi_p = inp.psi.at(p);
        let mut b = [T::ZERO; 3];
        if f < t.n_interior {
            let n = ld(t.neighbour, f) as usize;
            let gn = inp.grad.at(n);
            let xc_n = crate::fv::neighbour_centroid(f, n, t, &g.cell_centroid, &g.shift);
            let rn = sub(xf, xc_n);
            let psi_n = inp.psi.at(n);
            let w = ld(g.weight, f);
            let wn = T::ONE - w;
            // skewness: x_f − x_m, x_m the midpoint of P–N, where the
            // difference u_N − u_P gives the exact derivative of a quadratic
            let half = T::from_f64(0.5);
            let rm = [
                xf[0] - half * (xc_n[0] + xf[0] - rp[0]),
                xf[1] - half * (xc_n[1] + xf[1] - rp[1]),
                xf[2] - half * (xc_n[2] + xf[2] - rp[2]),
            ];
            let xc_p = g.cell_centroid.at(p);
            let xw = [
                w * xc_p[0] + wn * xc_n[0],
                w * xc_p[1] + wn * xc_n[1],
                w * xc_p[2] + wn * xc_n[2],
            ];
            let rw = sub(xf, xw);
            let (up, un) = (inp.uc.at(p), inp.uc.at(n));
            let mut i = 0;
            while i < 3 {
                let ho = if m >= T::ZERO {
                    psi_p[i] * dot(gp[i], rp)
                } else {
                    psi_n[i] * dot(gn[i], rn)
                };
                let gf = [
                    w * gp[i][0] + wn * gn[i][0],
                    w * gp[i][1] + wn * gn[i][1],
                    w * gp[i][2] + wn * gn[i][2],
                ];
                let conv = if inp.central > T::ZERO {
                    let uc = w * up[i] + wn * un[i] + dot(gf, rw);
                    let (uu, ud, gu, d, psi_u) = if m >= T::ZERO {
                        (up[i], un[i], gp[i], sub(xc_n, xc_p), psi_p[i])
                    } else {
                        (un[i], up[i], gn[i], sub(xc_p, xc_n), psi_n[i])
                    };
                    let uc = crate::fv::bounded_between(uc, uu, ud, inp.beta);
                    let gam = psi_u * crate::fv::central_weight(ud - uu, dot(gu, d), inp.beta);
                    m * ((T::ONE - inp.central) * inp.blend * ho + inp.central * gam * (uc - uu))
                } else {
                    inp.blend * m * ho
                };
                let skew = dot(sub(gn[i], gp[i]), rm);
                b[i] = conv - inp.nonorth * (mu * dot(kf, gf) + dd * skew);
                i += 1;
            }
            if inp.transposed > T::ZERO {
                // −μ Σ_j (∂u_j/∂x_i)_f S_j with the interpolated gradients,
                // less −μ (∇·u)_f S_i in the trace-free form
                let s = g.area.at(f);
                let div = inp.trace_free
                    * (w * (gp[0][0] + gp[1][1] + gp[2][2])
                        + wn * (gn[0][0] + gn[1][1] + gn[2][2]));
                let mut i = 0;
                while i < 3 {
                    let mut t = -div * s[i];
                    let mut jj = 0;
                    while jj < 3 {
                        t += (w * gp[jj][i] + wn * gn[jj][i]) * s[jj];
                        jj += 1;
                    }
                    b[i] -= mu * t;
                    i += 1;
                }
            }
            (mp + dd, mm - dd, b)
        } else {
            let fb = f - t.n_interior;
            let ub = inp.ub.at(fb);
            let mut i = 0;
            let kind = ld(inp.mom_kind, fb);
            // the central fraction takes the boundary value as the face
            // value of the outflow, ṁ⁺ (u_b − u_P)
            let lud = if inp.central > T::ZERO {
                (T::ONE - inp.central) * inp.blend
            } else {
                inp.blend
            };
            let up = inp.uc.at(p);
            if kind == MOM_OUTFLOW {
                while i < 3 {
                    b[i] = lud * mp * psi_p[i] * dot(gp[i], rp) + mm * ub[i];
                    if inp.central > T::ZERO {
                        b[i] += inp.central * mp * (ub[i] - up[i]);
                    }
                    i += 1;
                }
                (mp, T::ZERO, b)
            } else if kind == MOM_SYMMETRY {
                // n̂_i n̂_j = S_i S_j / |S|²; with the transposed stress the
                // normal stress is 2μ ∂u_n/∂n: on the plane u_n = 0, so the
                // transposed traction μ (∇u)ᵀ·S = μ |S| ∇u_n has only the
                // normal part μ |S| ∂u_n/∂n n̂, the part the Laplacian
                // takes once more. On the planes of a wedge it is half of
                // the hoop stress, −2μ v/r² (stage 6 part 6g)
                let dd = dd * (T::ONE + inp.transposed);
                let s = g.area.at(f);
                let a2 = dot(s, s);
                let up = inp.u.at(p);
                let su = dot(s, up);
                while i < 3 {
                    b[i] = dd * s[i] * (su - s[i] * up[i]) / a2;
                    i += 1;
                }
                (dd, T::ZERO, b)
            } else {
                let s = g.area.at(f);
                let div = inp.trace_free * (gp[0][0] + gp[1][1] + gp[2][2]);
                while i < 3 {
                    b[i] = lud * mp * psi_p[i] * dot(gp[i], rp) + mm * ub[i]
                        - dd * ub[i]
                        - inp.nonorth * mu * dot(kf, gp[i]);
                    if inp.transposed > T::ZERO {
                        b[i] -=
                            mu * (gp[0][i] * s[0] + gp[1][i] * s[1] + gp[2][i] * s[2] - div * s[i]);
                    }
                    if inp.central > T::ZERO {
                        b[i] += inp.central * mp * (ub[i] - up[i]);
                    }
                    i += 1;
                }
                (mp + dd, T::ZERO, b)
            }
        }
    }
}

/// Outputs of [`momentum_row`].
#[derive(Clone, Copy)]
pub struct MomentumRowOut<T> {
    /// Matrix values (SELL slots).
    pub val: *mut T,
    /// Right-hand sides of the three components.
    pub rhs: [*mut T; 3],
    /// D_P = V_P / a_P (relaxed diagonal with the time term), for the
    /// Rhie–Chow flux.
    pub d: *mut T,
    /// D^s_P = V_P / a_P (relaxed diagonal of the full operator, θ = 1,
    /// without the time term), for the time-step term of the Rhie–Chow flux
    /// (equal to `d` in steady flow).
    pub ds: *mut T,
    /// D_P^C = V_P / (a_P − Σ |a_N|) (SIMPLEC), for the pressure correction.
    pub dc: *mut T,
    /// Relaxed diagonal of each component (with the symmetry parts); the
    /// matrix itself holds the diagonal without them.
    pub diag: [*mut T; 3],
}

/// Inputs of [`momentum_row`].
#[derive(Clone, Copy)]
pub struct MomentumRowInputs<'a, T> {
    /// Connectivity.
    pub topo: Topology<'a>,
    /// Face coefficients c_P, c_N of [`momentum_face`].
    pub cp: &'a [T],
    /// See `cp`.
    pub cn: &'a [T],
    /// Deferred parts b_f of the three components.
    pub bf: Vec3<'a, T>,
    /// Pressure gradient per cell.
    pub grad_p: Vec3<'a, T>,
    /// Body force per unit volume (N/m³), the same in every cell.
    pub body_force: [T; 3],
    /// Cell volumes.
    pub volume: &'a [T],
    /// Current velocity (for the implicit under-relaxation).
    pub u: Vec3<'a, T>,
    /// Under-relaxation factor α_u.
    pub alpha: T,
    /// Matrix layout.
    pub layout: SellLayout<'a>,
    /// Momentum kinds of the boundary faces.
    pub mom_kind: &'a [u32],
    /// Face area vectors.
    pub face_area: Vec3<'a, T>,
    /// Weight θ of the implicit spatial terms (1 for steady flow, backward
    /// Euler and BDF2, ½ for Crank–Nicolson).
    pub theta: T,
    /// True in a time step: the time-derivative and old-time terms below
    /// enter; false (steady) skips them and their arrays may be empty.
    pub transient: bool,
    /// Cell density ρ_P, one entry per cell (read in a time step).
    pub rho: &'a [T],
    /// Coefficients a₀, −a₁, −a₂ of the backward differences
    /// du/dt ≈ (a₀ uⁿ⁺¹ + a₁ uⁿ + a₂ uⁿ⁻¹) / Δt: the diagonal gets
    /// τ_P = ρ_P a₀ / Δt V_P, the right-hand side the old levels uⁿ and
    /// uⁿ⁻¹ with ρ_P (−a₁) / Δt and ρ_P (−a₂) / Δt.
    pub time_a: [T; 3],
    /// Δt (or the pseudo-time step).
    pub dt: T,
    /// Velocity at the old time level n.
    pub u_old: Vec3<'a, T>,
    /// Velocity at the old time level n − 1 (read only when `c_old[1]` is
    /// not zero).
    pub u_old2: Vec3<'a, T>,
    /// Weight 1 − θ of the explicit spatial operator (Crank–Nicolson); 0
    /// skips `explicit`.
    pub explicit_weight: T,
    /// The spatial operator Σ_f ±F_f(uⁿ) of the old level per component
    /// ([`momentum_operator`]).
    pub explicit: Vec3<'a, T>,
}

/// Momentum row of cell `c`: Σ_f ±F_f = (−∇p_P + f_b) V_P per component,
/// relaxed implicitly (a_P/α, Patankar 1980), as [`crate::fv::assemble_row`]
/// with three right-hand sides. Also writes D_P = V_P / (a_P/α) and the
/// SIMPLEC factor D^C_P = V_P / (a_P/α + Σ a_N) (off-diagonals a_N ≤ 0),
/// bounded below by 10⁻¹² a_P/α.
///
/// In a time step (`transient`) the face terms are weighted with θ, the
/// diagonal gets the time term τ_P = ρ a₀ V_P / Δt and the right-hand side
/// the old levels ρ V_P (−a₁ uⁿ − a₂ uⁿ⁻¹) / Δt and, for Crank–Nicolson,
/// −(1 − θ) Σ_f ±F_f(uⁿ):
///
/// (θ a_P + τ_P)/α u_P + θ Σ a_N u_N = θ b + ρV(−a₁uⁿ − a₂uⁿ⁻¹)/Δt
///   − (1 − θ) R(uⁿ) + (1 − α)/α (θ a_P + τ_P) u_P^prev.
///
/// D_P then contains τ_P; D^s_P = V_P / (θ a_P/α) does not. With θ = 1 and
/// no time term the arithmetic is exactly that of the steady row.
///
/// # Safety
/// As [`crate::fv::assemble_row`]; the output arrays have an entry for `c`
/// and no other thread writes row `c`.
#[inline(always)]
pub unsafe fn momentum_row<T: Real>(c: usize, inp: &MomentumRowInputs<T>, out: &MomentumRowOut<T>) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let t = &inp.topo;
        let lo = ld(t.cell_ptr, c) as usize;
        let hi = ld(t.cell_ptr, c + 1) as usize;
        let vol = ld(inp.volume, c);
        let gp = inp.grad_p.at(c);
        let mut b = [
            (inp.body_force[0] - gp[0]) * vol,
            (inp.body_force[1] - gp[1]) * vol,
            (inp.body_force[2] - gp[2]) * vol,
        ];
        let mut diag = T::ZERO;
        let mut sym = [T::ZERO; 3];
        let mut off = T::ZERO;
        let mut j = 1;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(t.cell_faces, k));
            let bf = inp.bf.at(f);
            let th = inp.theta;
            if is_nb {
                diag -= th * ld(inp.cn, f);
                b[0] += th * bf[0];
                b[1] += th * bf[1];
                b[2] += th * bf[2];
                let a = -(th * ld(inp.cp, f));
                off += a;
                *out.val.add(inp.layout.pos(c, j)) = a;
                j += 1;
            } else {
                if f >= t.n_interior && ld(inp.mom_kind, f - t.n_interior) == MOM_SYMMETRY {
                    let s = inp.face_area.at(f);
                    let a2 = dot(s, s);
                    let c = ld(inp.cp, f);
                    sym[0] += th * (c * s[0] * s[0] / a2);
                    sym[1] += th * (c * s[1] * s[1] / a2);
                    sym[2] += th * (c * s[2] * s[2] / a2);
                } else {
                    diag += th * ld(inp.cp, f);
                }
                b[0] -= th * bf[0];
                b[1] -= th * bf[1];
                b[2] -= th * bf[2];
                if f < t.n_interior {
                    let a = th * ld(inp.cn, f);
                    off += a;
                    *out.val.add(inp.layout.pos(c, j)) = a;
                    j += 1;
                }
            }
            k += 1;
        }
        let u = inp.u.at(c);
        // the face terms of the diagonal carry θ; D^s is that of the full
        // operator (Crank–Nicolson: the steady state is then the SIMPLEC
        // one, and the face flux follows the same θ-scheme as the cells);
        // × θ is exact for θ = 1
        let ds = vol / (diag / inp.alpha) * inp.theta;
        let mut src = [T::ZERO; 3];
        if inp.transient {
            // ρ_P a₀ / Δt in this order: a constant density gives the
            // former single coefficient bit for bit
            let r = ld(inp.rho, c);
            diag += r * inp.time_a[0] / inp.dt * vol;
            let uo = inp.u_old.at(c);
            let c1 = r * inp.time_a[1] / inp.dt;
            let mut i = 0;
            while i < 3 {
                src[i] = c1 * vol * uo[i];
                i += 1;
            }
            if inp.time_a[2] != T::ZERO {
                let c2 = r * inp.time_a[2] / inp.dt;
                let uo2 = inp.u_old2.at(c);
                i = 0;
                while i < 3 {
                    src[i] += c2 * vol * uo2[i];
                    i += 1;
                }
            }
            if inp.explicit_weight != T::ZERO {
                let r = inp.explicit.at(c);
                i = 0;
                while i < 3 {
                    src[i] -= inp.explicit_weight * r[i];
                    i += 1;
                }
            }
        }
        let dr = diag / inp.alpha;
        let mut i = 0;
        while i < 3 {
            let a = diag + sym[i];
            let ar = a / inp.alpha;
            let rhs = if inp.transient { b[i] + src[i] } else { b[i] };
            *out.rhs[i].add(c) = rhs + (ar - a) * u[i];
            *out.diag[i].add(c) = ar;
            i += 1;
        }
        *out.val.add(inp.layout.pos(c, 0)) = dr;
        *out.d.add(c) = vol / dr;
        *out.ds.add(c) = ds;
        let dcd = (dr + off).max_of(T::from_f64(1e-12) * dr);
        *out.dc.add(c) = vol / dcd;
    }
}

/// Inputs of [`face_mass_flux`].
#[derive(Clone, Copy)]
pub struct FluxInputs<'a, T> {
    /// Connectivity.
    pub topo: Topology<'a>,
    /// Geometry.
    pub geo: FaceGeometry<'a, T>,
    /// Face density ρ_f (kg/m³), one entry per face.
    pub rho_face: &'a [T],
    /// Cell velocity.
    pub u: Vec3<'a, T>,
    /// Gradient of the velocity used for the skewness correction of the
    /// face interpolation.
    pub grad: VelocityGradient<'a, T>,
    /// Cell pressure.
    pub p: &'a [T],
    /// Cell pressure gradient.
    pub grad_p: Vec3<'a, T>,
    /// Boundary face pressure.
    pub pb: &'a [T],
    /// Prescribed boundary pressure: the trial flux of a
    /// [`FLUX_PRESSURE_NO_BACKFLOW`] face takes it in place of `pb` (which
    /// is the extrapolated pressure while the face is closed), so a closed
    /// face opens only when the prescribed pressure would drive the flow out.
    pub pset: &'a [T],
    /// Boundary face velocity.
    pub ub: Vec3<'a, T>,
    /// Flux kind per boundary face ([`FLUX_FIXED`], [`FLUX_PRESSURE`],
    /// [`FLUX_PRESSURE_NO_BACKFLOW`]).
    pub flux_kind: &'a [u32],
    /// Momentum kind per boundary face: on [`MOM_OUTFLOW`] faces the
    /// boundary velocity `ub` (zero normal gradient) enters the Rhie–Chow
    /// flux, on pressure inlets the owner velocity extrapolated to the face.
    pub mom_kind: &'a [u32],
    /// D_P = V_P / a_P per cell (0 for the plain interpolated flux).
    pub d: &'a [T],
    /// Relaxation term per face, (1 − α_u)(ṁ_f − ρ ū_f·S)ⁿ⁻¹, added as is.
    pub relax: &'a [T],
    /// True to add the time-step term (a time step with the pressure
    /// terms); false skips it, `ds` and `time` may then be empty.
    pub transient: bool,
    /// D^s_P per cell, D_P without the time term ([`MomentumRowOut::ds`]).
    pub ds: &'a [T],
    /// Old-level part per face, Y_f = Σ_k (−a_k/a₀)(ṁ_f − ρ ū_f·S)ᵏ over
    /// the old levels k = n, n − 1.
    pub time: &'a [T],
    /// Under-relaxation factor α_u of the momentum equations.
    pub alpha: T,
    /// Limiter of the velocity gradient per cell and component (read when
    /// `bounded`).
    pub psi: Vec3<'a, T>,
    /// True for the bounded face velocity of the cell-quality remediation
    /// (D50): the mean of the two limited extrapolations to the face
    /// centroid, w (u_P + ψ_P ∇u_P·r_P) + (1 − w)(u_N + ψ_N ∇u_N·r_N),
    /// in place of the interpolation with the interpolated gradient.
    pub bounded: bool,
    /// Balanced-force offset q_f per face (D154): subtracted from the
    /// pressure jump of interior and pressure-boundary faces; empty
    /// without gravity.
    pub jump_offset: &'a [T],
}

/// Face mass flux (kg/s) out of the owner, with the Rhie–Chow interpolation
/// (Rhie & Chow 1983) in the form whose converged result does not depend on
/// the under-relaxation (Majumdar 1988):
///
/// ṁ_f = ρ ū_f·S − ρ D̄_f |S|²/(d·S) [(p_N − p_P) − ∇p̄_f·d] + relax_f,
///
/// with ū_f the velocity interpolated linearly to the point x_w = w x_P +
/// (1 − w) x_N of the line P–N and corrected for skewness with the
/// interpolated gradient, ū_f = w u_P + (1 − w) u_N + ∇ū_f·(x_f − x_w), and
/// D̄_f, ∇p̄_f interpolated with the same weight w. (The interpolated
/// gradient, not the two cell gradients, carries the correction: the
/// first-order errors of the least-squares gradients of neighbouring cells
/// largely cancel in it.) On [`FLUX_FIXED`] boundary faces
/// ṁ_f = ρ u_b·S; on [`FLUX_PRESSURE`] faces the face velocity (the
/// boundary value at outlets, the owner value extrapolated to the face at
/// pressure inlets) and the boundary pressure enter the same formula. A linear
/// velocity and a linear pressure give the exact flux on every mesh.
///
/// **Time-step term (D27).** In a time step the flux gets
/// γ_f Y_f with Y_f the old-level differences (ṁ − ρ ū·S) of
/// [`FluxInputs::time`] and γ_f = α (1 − D̄_f / D̄^s_f), D̄^s the
/// interpolated D without the time term. At a steady state
/// (ṁ − ρ ū·S)(α − γ_f) = −ρ D̄_f |S|²/(d·S) \[jump\], so the face flux equals
/// the one of the steady solver with D̄^s/α: the converged result depends
/// neither on Δt nor on α (the analogue of Majumdar's relaxation term; the
/// cell-centred form ρ D a₀/Δt is only approximately independent on
/// non-uniform meshes).
///
/// # Safety
/// As [`momentum_face`].
#[inline(always)]
pub unsafe fn face_mass_flux<T: Real>(f: usize, inp: &FluxInputs<T>) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let t = &inp.topo;
        let g = &inp.geo;
        let p = ld(t.owner, f) as usize;
        let s = g.area.at(f);
        let xf = g.face_centroid.at(f);
        let xp = g.cell_centroid.at(p);
        let rp = sub(xf, xp);
        let up = inp.u.at(p);
        let gp = inp.grad.at(p);
        let gpp = inp.grad_p.at(p);
        let dgeo = ld(g.dgeo, f);
        if f < t.n_interior {
            let n = ld(t.neighbour, f) as usize;
            let xn = crate::fv::neighbour_centroid(f, n, t, &g.cell_centroid, &g.shift);
            let un = inp.u.at(n);
            let gn = inp.grad.at(n);
            let w = ld(g.weight, f);
            let wn = T::ONE - w;
            // linear interpolation to the point x_w of the line P–N, then
            // along the face with the interpolated gradient
            let xw = [
                w * xp[0] + wn * xn[0],
                w * xp[1] + wn * xn[1],
                w * xp[2] + wn * xn[2],
            ];
            let rw = sub(xf, xw);
            let mut uf = [T::ZERO; 3];
            let mut i = 0;
            while i < 3 {
                let gf = [
                    w * gp[i][0] + wn * gn[i][0],
                    w * gp[i][1] + wn * gn[i][1],
                    w * gp[i][2] + wn * gn[i][2],
                ];
                uf[i] = if inp.bounded {
                    let (sp, sn) = (inp.psi.at(p), inp.psi.at(n));
                    let rn = sub(xf, xn);
                    w * (up[i] + sp[i] * dot(gp[i], rp)) + wn * (un[i] + sn[i] * dot(gn[i], rn))
                } else {
                    w * up[i] + wn * un[i] + dot(gf, rw)
                };
                i += 1;
            }
            let gpn = inp.grad_p.at(n);
            let gpf = [
                w * gpp[0] + wn * gpn[0],
                w * gpp[1] + wn * gpn[1],
                w * gpp[2] + wn * gpn[2],
            ];
            let d = sub(xn, xp);
            let df = w * ld(inp.d, p) + wn * ld(inp.d, n);
            let mut jump = (ld(inp.p, n) - ld(inp.p, p)) - dot(gpf, d);
            if !inp.jump_offset.is_empty() {
                jump -= ld(inp.jump_offset, f);
            }
            let m = ld(inp.rho_face, f) * (dot(uf, s) - df * dgeo * jump) + ld(inp.relax, f);
            if inp.transient {
                let dsf = w * ld(inp.ds, p) + wn * ld(inp.ds, n);
                m + inp.alpha * (T::ONE - df / dsf) * ld(inp.time, f)
            } else {
                m
            }
        } else {
            let fb = f - t.n_interior;
            let kind = ld(inp.flux_kind, fb);
            if kind != FLUX_FIXED {
                let uf = if ld(inp.mom_kind, fb) == MOM_OUTFLOW {
                    inp.ub.at(fb)
                } else {
                    [
                        up[0] + dot(gp[0], rp),
                        up[1] + dot(gp[1], rp),
                        up[2] + dot(gp[2], rp),
                    ]
                };
                let pb = if kind == FLUX_PRESSURE_NO_BACKFLOW {
                    ld(inp.pset, fb)
                } else {
                    ld(inp.pb, fb)
                };
                let mut jump = (pb - ld(inp.p, p)) - dot(gpp, rp);
                if !inp.jump_offset.is_empty() {
                    jump -= ld(inp.jump_offset, f);
                }
                let dp = ld(inp.d, p);
                let m = ld(inp.rho_face, f) * (dot(uf, s) - dp * dgeo * jump) + ld(inp.relax, f);
                if inp.transient {
                    m + inp.alpha * (T::ONE - dp / ld(inp.ds, p)) * ld(inp.time, f)
                } else {
                    m
                }
            } else {
                ld(inp.rho_face, f) * dot(inp.ub.at(fb), s)
            }
        }
    }
}

/// Coefficient c_f of the pressure correction on face `f`: the mass-flux
/// correction is ṁ'_f = −c_f (p'_N − p'_P), c_f = ρ D̄^C_f |S|²/(d·S)
/// (SIMPLEC, orthogonal part); on pressure boundary faces
/// ([`FLUX_PRESSURE`], open [`FLUX_PRESSURE_NO_BACKFLOW`]) p'_b = 0 and
/// c_f = ρ D^C_P |S|²/(d·S); on [`FLUX_FIXED`] faces and closed faces
/// (`closed` = 1 per boundary face, see [`close_backflow`]) c_f = 0.
///
/// # Safety
/// As [`momentum_face`]; `dc` has an entry per cell, `closed` one per
/// boundary face.
#[inline(always)]
pub unsafe fn pressure_face<T: Real>(
    f: usize,
    topo: &Topology,
    geo: &FaceGeometry<T>,
    flux_kind: &[u32],
    rho_face: &[T],
    dc: &[T],
    closed: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let p = ld(topo.owner, f) as usize;
        let dgeo = ld(geo.dgeo, f);
        if f < topo.n_interior {
            let n = ld(topo.neighbour, f) as usize;
            let w = ld(geo.weight, f);
            ld(rho_face, f) * (w * ld(dc, p) + (T::ONE - w) * ld(dc, n)) * dgeo
        } else {
            let fb = f - topo.n_interior;
            let kind = ld(flux_kind, fb);
            if kind == FLUX_PRESSURE
                || (kind == FLUX_PRESSURE_NO_BACKFLOW && ld(closed, fb) == T::ZERO)
            {
                ld(rho_face, f) * ld(dc, p) * dgeo
            } else {
                T::ZERO
            }
        }
    }
}

/// Backflow suppression of a boundary face of flux kind `kind` with the
/// predicted flux `m` (after the predicted face fluxes): a
/// [`FLUX_PRESSURE_NO_BACKFLOW`] face whose predicted flux enters the
/// domain (ṁ*_f < 0) is closed for this outer iteration, ṁ*_f = 0 and
/// closed = 1, so that [`pressure_face`] and [`delta_flux`] give it no
/// correction; every other face keeps its flux and gets closed = 0. The
/// face opens again in an iteration whose predicted flux leaves the domain.
/// Returns `(ṁ_f, closed)`.
#[inline(always)]
pub fn close_backflow<T: Real>(kind: u32, m: T) -> (T, T) {
    if kind == FLUX_PRESSURE_NO_BACKFLOW && m < T::ZERO {
        (T::ZERO, T::ONE)
    } else {
        (m, T::ZERO)
    }
}

/// Pressure-correction row of cell `c`: Σ_f c_f (p'_P − p'_N) = −Σ_f ±ṁ*_f
/// (the continuity defect of the predicted fluxes). Writes the SELL slots
/// and `rhs[c]`; returns nothing else. The matrix is symmetric and, with at
/// least one pressure boundary or a reference cell, positive definite.
///
/// # Safety
/// As [`crate::fv::assemble_row`].
#[inline(always)]
pub unsafe fn pressure_row<T: Real>(
    c: usize,
    topo: &Topology,
    cf: &[T],
    mass_flux: &[T],
    layout: &SellLayout,
    val: *mut T,
    rhs: *mut T,
) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut diag = T::ZERO;
        let mut b = T::ZERO;
        let mut j = 1;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let a = ld(cf, f);
            diag += a;
            if is_nb {
                b += ld(mass_flux, f);
            } else {
                b -= ld(mass_flux, f);
            }
            if f < topo.n_interior {
                *val.add(layout.pos(c, j)) = -a;
                j += 1;
            }
            k += 1;
        }
        *val.add(layout.pos(c, 0)) = diag;
        *rhs.add(c) = b;
    }
}

/// Net mass flux out of cell `c`, Σ_f ±ṁ_f (kg/s): the continuity defect.
///
/// # Safety
/// As [`crate::fv::lsq_gradient`]; `mass_flux` has an entry per face.
#[inline(always)]
pub unsafe fn net_outflow<T: Real>(c: usize, topo: &Topology, mass_flux: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut s = T::ZERO;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            if is_nb {
                s -= ld(mass_flux, f);
            } else {
                s += ld(mass_flux, f);
            }
            k += 1;
        }
        s
    }
}

/// Throughflow of cell `c`, Σ_f |ṁ_f| (kg/s), the scale of the continuity
/// residual.
///
/// # Safety
/// As [`net_outflow`].
#[inline(always)]
pub unsafe fn throughflow<T: Real>(c: usize, topo: &Topology, mass_flux: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut s = T::ZERO;
        let mut k = lo;
        while k < hi {
            let (f, _) = unpack_face(ld(topo.cell_faces, k));
            s += ld(mass_flux, f).abs();
            k += 1;
        }
        s
    }
}

/// Corrected mass flux of face `f`, ṁ_f = ṁ*_f − c_f (p'_N − p'_P), with
/// p'_b = 0 on pressure boundary faces (where c_f is non-zero).
///
/// # Safety
/// As [`pressure_face`]; `pc` has an entry per cell.
#[inline(always)]
pub unsafe fn correct_flux<T: Real>(
    f: usize,
    topo: &Topology,
    cf: &[T],
    mass_flux: &[T],
    pc: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let p = ld(topo.owner, f) as usize;
        let pn = if f < topo.n_interior {
            ld(pc, ld(topo.neighbour, f) as usize)
        } else {
            T::ZERO
        };
        ld(mass_flux, f) - ld(cf, f) * (pn - ld(pc, p))
    }
}

/// Inputs of the boundary updates.
#[derive(Clone, Copy)]
pub struct BoundaryInputs<'a, T> {
    /// Connectivity.
    pub topo: Topology<'a>,
    /// Geometry.
    pub geo: FaceGeometry<'a, T>,
    /// Face density ρ_f, one entry per face.
    pub rho_face: &'a [T],
    /// Face mass flux.
    pub mass_flux: &'a [T],
}

/// New velocity of boundary face `fb` by its update rule (`UB_*`): the
/// prescribed value `ub`, the symmetry value u_P − (u_P·n̂) n̂, the
/// tangential extrapolation u_e = u_P + ∇u_P·r_t, the velocity normal to
/// the face from the mass flux, the outlet value with backflow treatment
/// ([`UB_OUTLET`]), the tangential part on a closed face ([`UB_NO_BACKFLOW`])
/// or the convective outlet value (c + β u_e)/(1 + β)
/// with `conv` = (c, β) of the face ([`UB_CONVECTIVE`]).
///
/// # Safety
/// As [`momentum_face`]; `fb` is a boundary face index.
#[inline(always)]
pub unsafe fn boundary_velocity<T: Real>(
    fb: usize,
    inp: &BoundaryInputs<T>,
    rule: u32,
    ub: [T; 3],
    conv: ([T; 3], T),
    u: &Vec3<T>,
    grad: &VelocityGradient<T>,
) -> [T; 3] {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let f = fb + inp.topo.n_interior;
        let p = ld(inp.topo.owner, f) as usize;
        let s = inp.geo.area.at(f);
        let up = u.at(p);
        if rule == UB_SYMMETRY {
            let ss = dot(s, s);
            let un = dot(up, s) / ss;
            [up[0] - un * s[0], up[1] - un * s[1], up[2] - un * s[2]]
        } else if rule == UB_EXTRAPOLATE
            || rule == UB_OUTLET
            || rule == UB_CONVECTIVE
            || rule == UB_NO_BACKFLOW
        {
            if rule == UB_CONVECTIVE && conv.1 == T::ZERO {
                return conv.0;
            }
            let r = sub(inp.geo.face_centroid.at(f), inp.geo.cell_centroid.at(p));
            let rn = dot(r, s) / dot(s, s);
            let r = [r[0] - rn * s[0], r[1] - rn * s[1], r[2] - rn * s[2]];
            let g = grad.at(p);
            let ue = [
                up[0] + dot(g[0], r),
                up[1] + dot(g[1], r),
                up[2] + dot(g[2], r),
            ];
            if rule == UB_OUTLET && ld(inp.mass_flux, f) < T::ZERO {
                let q = dot(ue, s) / dot(s, s);
                [q * s[0], q * s[1], q * s[2]]
            } else if rule == UB_NO_BACKFLOW && ld(inp.mass_flux, f) <= T::ZERO {
                let q = dot(ue, s) / dot(s, s);
                [ue[0] - q * s[0], ue[1] - q * s[1], ue[2] - q * s[2]]
            } else if rule == UB_CONVECTIVE {
                let (c, beta) = conv;
                let d = T::ONE + beta;
                [
                    (c[0] + beta * ue[0]) / d,
                    (c[1] + beta * ue[1]) / d,
                    (c[2] + beta * ue[2]) / d,
                ]
            } else {
                ue
            }
        } else if rule == UB_FROM_FLUX {
            let ss = dot(s, s);
            let q = ld(inp.mass_flux, f) / (ld(inp.rho_face, f) * ss);
            [q * s[0], q * s[1], q * s[2]]
        } else {
            ub
        }
    }
}

/// New pressure of boundary face `fb` by its update rule (`PB_*`): the
/// prescribed value `pset`, the extrapolation p_P + ∇p_P·(x_f − x_P), the
/// static pressure of the total pressure `pset`,
/// p₀ − ½ ρ (ṁ_f / (ρ |S|))², or at an outlet with backflow treatment the
/// static pressure `pset` on outflow and the total pressure `pset` on
/// backflow ([`PB_OUTLET`]); with backflow suppressed the extrapolation on
/// the closed faces ([`PB_NO_BACKFLOW`], [`PB_TOTAL_NO_BACKFLOW`]). With
/// gravity the extrapolation adds the balanced-force term `offset[fb]`
/// (D154 (4)), so that the reduced jump has a zero normal gradient;
/// `offset` is empty without gravity.
///
/// # Safety
/// As [`boundary_velocity`]; `offset` is empty or has an entry per
/// boundary face.
#[inline(always)]
pub unsafe fn boundary_pressure<T: Real>(
    fb: usize,
    inp: &BoundaryInputs<T>,
    rule: u32,
    pset: T,
    p: &[T],
    grad_p: &Vec3<T>,
    offset: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let f = fb + inp.topo.n_interior;
        let c = ld(inp.topo.owner, f) as usize;
        let m = ld(inp.mass_flux, f);
        let closed = (rule == PB_NO_BACKFLOW || rule == PB_TOTAL_NO_BACKFLOW) && m <= T::ZERO;
        if rule == PB_EXTRAPOLATE || closed {
            let r = sub(inp.geo.face_centroid.at(f), inp.geo.cell_centroid.at(c));
            let pb = ld(p, c) + dot(grad_p.at(c), r);
            if offset.is_empty() {
                pb
            } else {
                pb + ld(offset, fb)
            }
        } else if rule == PB_TOTAL
            || rule == PB_TOTAL_NO_BACKFLOW
            || (rule == PB_OUTLET && m < T::ZERO)
        {
            let s = inp.geo.area.at(f);
            let rho = ld(inp.rho_face, f);
            let un = ld(inp.mass_flux, f) / rho;
            pset - T::from_f64(0.5) * rho * un * un / dot(s, s)
        } else {
            pset
        }
    }
}

/// Inputs of [`momentum_operator`]: the face coefficients of
/// [`momentum_face`] for the old level.
#[derive(Clone, Copy)]
pub struct OperatorInputs<'a, T> {
    /// Connectivity.
    pub topo: Topology<'a>,
    /// Face coefficients c_P, c_N.
    pub cp: &'a [T],
    /// See `cp`.
    pub cn: &'a [T],
    /// Deferred parts b_f of the three components.
    pub bf: Vec3<'a, T>,
    /// Momentum kinds of the boundary faces.
    pub mom_kind: &'a [u32],
    /// Face area vectors.
    pub face_area: Vec3<'a, T>,
    /// The velocity the coefficients act on.
    pub u: Vec3<'a, T>,
}

/// The spatial operator of the momentum equations of cell `c`,
/// R_i = Σ_f ±F_f,i (out of the cell positive) with
/// F_f,i = c_P u_owner,i + c_N u_neighbour,i + b_f,i, and on symmetry faces
/// c_P n̂_i² u_P,i + b_f,i: the same terms as the row of [`momentum_row`],
/// evaluated with a known velocity. Crank–Nicolson uses it for the explicit
/// half of the old level.
///
/// # Safety
/// As [`momentum_row`].
#[inline(always)]
pub unsafe fn momentum_operator<T: Real>(c: usize, inp: &OperatorInputs<T>) -> [T; 3] {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let t = &inp.topo;
        let lo = ld(t.cell_ptr, c) as usize;
        let hi = ld(t.cell_ptr, c + 1) as usize;
        let uc = inp.u.at(c);
        let mut r = [T::ZERO; 3];
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(t.cell_faces, k));
            let cp = ld(inp.cp, f);
            let cn = ld(inp.cn, f);
            let bf = inp.bf.at(f);
            let mut i = 0;
            if f < t.n_interior {
                let o = if is_nb {
                    ld(t.owner, f) as usize
                } else {
                    ld(t.neighbour, f) as usize
                };
                let uo = inp.u.at(o);
                while i < 3 {
                    if is_nb {
                        r[i] -= cp * uo[i] + cn * uc[i] + bf[i];
                    } else {
                        r[i] += cp * uc[i] + cn * uo[i] + bf[i];
                    }
                    i += 1;
                }
            } else if ld(inp.mom_kind, f - t.n_interior) == MOM_SYMMETRY {
                let s = inp.face_area.at(f);
                let a2 = dot(s, s);
                while i < 3 {
                    r[i] += (cp * s[i] * s[i] / a2) * uc[i] + bf[i];
                    i += 1;
                }
            } else {
                while i < 3 {
                    r[i] += cp * uc[i] + bf[i];
                    i += 1;
                }
            }
            k += 1;
        }
        r
    }
}

/// Σ_N val_PN x_N over the off-diagonal slots of row `c` (the interior
/// faces of the cell in stored order, slot 1, 2, …), for the neighbour
/// correction of the PISO correctors: H(x)_P = −Σ_N val_PN x_N.
///
/// # Safety
/// As [`momentum_row`]; `val` is the SELL storage of `layout`, `x` has an
/// entry per cell.
#[inline(always)]
pub unsafe fn neighbour_sum<T: Real>(
    c: usize,
    topo: &Topology,
    layout: &SellLayout,
    val: &[T],
    x: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut s = T::ZERO;
        let mut j = 1;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            if f < topo.n_interior {
                let o = if is_nb {
                    ld(topo.owner, f)
                } else {
                    ld(topo.neighbour, f)
                };
                s += ld(val, layout.pos(c, j)) * ld(x, o as usize);
                j += 1;
            }
            k += 1;
        }
        s
    }
}

/// Σ_N val_PN (x_N − own x_P) of cell `c`: with `own` = 0
/// [`neighbour_sum`], with `own` = 1 the neighbour sum less the
/// off-diagonal row sum times x_P (the neighbour correction of PISO after a
/// SIMPLEC-type corrector, whose D^C already counts Σ_N a_PN x_P).
///
/// # Safety
/// As [`neighbour_sum`].
#[inline(always)]
pub unsafe fn neighbour_difference<T: Real>(
    c: usize,
    topo: &Topology,
    layout: &SellLayout,
    val: &[T],
    x: &[T],
    own: T,
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let xp = own * ld(x, c);
        let mut s = T::ZERO;
        let mut j = 1;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            if f < topo.n_interior {
                let o = if is_nb {
                    ld(topo.owner, f)
                } else {
                    ld(topo.neighbour, f)
                };
                s += ld(val, layout.pos(c, j)) * (ld(x, o as usize) - xp);
                j += 1;
            }
            k += 1;
        }
        s
    }
}

/// Mass flux of a cell velocity increment δ through face `f`:
/// ρ (w δ_P + (1 − w) δ_N)·S on interior faces, ρ δ_P·S on open pressure
/// boundaries and 0 where the boundary fixes the flux or a closed face
/// (`closed`, [`close_backflow`]) has none.
///
/// # Safety
/// As [`pressure_face`]; `delta` has an entry per cell.
#[inline(always)]
pub unsafe fn delta_flux<T: Real>(
    f: usize,
    topo: &Topology,
    geo: &FaceGeometry<T>,
    flux_kind: &[u32],
    rho_face: &[T],
    delta: &Vec3<T>,
    closed: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let p = ld(topo.owner, f) as usize;
        let s = geo.area.at(f);
        let dp = delta.at(p);
        if f < topo.n_interior {
            let n = ld(topo.neighbour, f) as usize;
            let dn = delta.at(n);
            let w = ld(geo.weight, f);
            let wn = T::ONE - w;
            let df = [
                w * dp[0] + wn * dn[0],
                w * dp[1] + wn * dn[1],
                w * dp[2] + wn * dn[2],
            ];
            ld(rho_face, f) * dot(df, s)
        } else {
            let fb = f - topo.n_interior;
            let kind = ld(flux_kind, fb);
            if kind == FLUX_PRESSURE
                || (kind == FLUX_PRESSURE_NO_BACKFLOW && ld(closed, fb) == T::ZERO)
            {
                ld(rho_face, f) * dot(dp, s)
            } else {
                T::ZERO
            }
        }
    }
}

/// Courant number of cell `c`, `factor` Σ_f |ṁ_f| / V_c with
/// `factor` = Δt / (2ρ): the outflow volume of one time step over the cell
/// volume (half the throughflow is the outflow in a divergence-free field).
///
/// # Safety
/// As [`throughflow`]; `volume` has an entry for `c`.
#[inline(always)]
pub unsafe fn courant<T: Real>(
    c: usize,
    topo: &Topology,
    mass_flux: &[T],
    volume: &[T],
    factor: T,
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { factor * throughflow(c, topo, mass_flux) / ld(volume, c) }
}

/// Parts of the momentum residual scale of cell `c`: per component the
/// unrelaxed a_P u_P,i and b_i − (a_P/α − a_P) u_P,i, from the relaxed
/// diagonal `diag` (a_P/α) and the right-hand side `rhs` of
/// [`momentum_row`]. Returns `(a u, b)` per component.
///
/// # Safety
/// Every array has an entry for `c`.
#[inline(always)]
pub unsafe fn residual_scale<T: Real>(
    c: usize,
    u: &Vec3<T>,
    rhs: &Vec3<T>,
    diag: &Vec3<T>,
    alpha: T,
) -> ([T; 3], [T; 3]) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let uc = u.at(c);
        let b = rhs.at(c);
        let dg = diag.at(c);
        let mut v = [T::ZERO; 3];
        let mut w = [T::ZERO; 3];
        let mut i = 0;
        while i < 3 {
            let dr = dg[i];
            let a = dr * alpha;
            v[i] = a * uc[i];
            w[i] = b[i] - (dr - a) * uc[i];
            i += 1;
        }
        (v, w)
    }
}

/// D̂ of the coupled preconditioner (the momentum inverse of its Schur
/// complement): the SIMPLEC factor D^C = V/(a_P/α + Σ a_N) bounded by its
/// value for face fluxes that satisfy continuity, 1/(a/D^s + b) with
/// a = 1 − α and b = τ/α (τ = ρ/Δτ of the pseudo-time term, D^s the factor
/// without it). With fluxes that do not satisfy continuity yet (a start, a
/// loose solve) a_P + Σ a_N can be negative and D^C without the bound
/// arbitrarily large.
#[inline(always)]
pub fn schur_diagonal<T: Real>(a: T, b: T, ds: T, dc: T) -> T {
    dc.min_of(T::ONE / (a / ds + b))
}

/// a (x_i − y_i).
///
/// # Safety
/// `i` is inside both slices.
#[inline(always)]
pub unsafe fn scaled_difference<T: Real>(i: usize, a: T, x: &[T], y: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { a * (ld(x, i) - ld(y, i)) }
}

/// a x_i + b y_i.
///
/// # Safety
/// `i` is inside both slices.
#[inline(always)]
pub unsafe fn lincomb<T: Real>(i: usize, a: T, x: &[T], b: T, y: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe { a * ld(x, i) + b * ld(y, i) }
}

/// Boundary value b of the least-squares gradient of boundary face `fb`:
/// 0 where the variable has a zero normal gradient (`zero_gradient` = 1,
/// the face value then follows the owner), else the face value.
///
/// # Safety
/// `fb` is inside both slices.
#[inline(always)]
pub unsafe fn boundary_constant<T: Real>(fb: usize, zero_gradient: &[T], value: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        if ld(zero_gradient, fb) == T::ONE {
            T::ZERO
        } else {
            ld(value, fb)
        }
    }
}
