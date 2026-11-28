//! Kernel bodies of the two-fluid (Euler–Euler) model of dispersed bubbles
//! in a liquid (L3, stage MP3, D352, D353): the interphase forces per cell
//! and the face arrays of the two momentum equations.
//!
//! Like [`crate::flow`], every body computes the outputs of one cell or one
//! face from its inputs (a gather, no atomics), so a CPU loop and one GPU
//! thread per index give bitwise the same results (D31).
//!
//! **Forces** (per unit volume, on the gas; the liquid takes the opposite):
//!
//! - drag K (u_l − u_g), K = ¾ C_D ρ_l α_g |u_r| / d with C_D of Ishii &
//!   Zuber (1979): the largest of the viscous sphere 24/Re (1 + 0.1
//!   Re^0.75) and the smaller of the distorted bubble ⅔ √Eo and the cap
//!   8/3, Re with the mixture viscosity μ_l/(1 − α_g) of the bubbly
//!   regime. K is written as ¾ (C_D Re) μ_m α_g / d², which stays finite
//!   for u_r → 0 (Stokes drag);
//! - lift −C_L ρ_l α_g (u_g − u_l) × (∇ × u_l) with C_L of Tomiyama et al.
//!   (2002): min(0.288 tanh(0.121 Re), f(Eo_d)) for Eo_d < 4, f(Eo_d) =
//!   0.00105 Eo_d³ − 0.0159 Eo_d² − 0.0204 Eo_d + 0.474 for 4 ≤ Eo_d ≤ 10,
//!   −0.27 above, Eo_d of the horizontal diameter d_h = d (1 + 0.163
//!   Eo^0.757)^⅓ (Wellek et al. 1966); or a constant C_L;
//! - wall lubrication of Hosokawa et al. (2002) in the form of Rzehak et al.
//!   (2012): (2/d) C_W α_g ρ_l |u_r,∥|² n̂_w with C_W = f(Eo) (d/(2y))²,
//!   f = max(7/Re^1.9, 0.0217 Eo), u_r,∥ the slip along the wall, y at
//!   least d/2;
//! - turbulent dispersion of Burns et al. (2004, the Favre-averaged drag):
//!   −C_TD K (ν_t/σ_TD) (1/α_g + 1/α_l) ∇α_g, σ_TD = 0.9, a gradient force
//!   that the face offsets carry (below), like gravity;
//! - bubble-induced turbulence of Sato et al. (1981): μ_BIT = C_μb ρ_l
//!   α_g d |u_r| on the liquid's viscosity (C_μb = 0.6).
//!
//! **Pressure and gravity** (D352). The pressure is p_rgh = p − ρ_l g·x:
//! the force α_k (−∇p + ρ_k g) of phase k is −α_k ∇p_rgh + α_k (ρ_k − ρ_l)
//! g. The rows take α_k (∇p − b_k) with the cell force per unit phase
//! volume b_k (gravity and dispersion over α_k); the Rhie–Chow flux of
//! phase k, whose D is α_k V/a_P, takes the same b_k as its cell body
//! force and the face offsets o_k,f/α_k,f, so that the face flux reacts
//! to the pressure and to these forces as the cell velocities do.

use crate::fv::{Topology, Vec3};
use crate::ld;
use foehncraft_core::Real;

/// No lift.
pub const LIFT_NONE: u32 = 0;
/// Lift of Tomiyama et al. (2002).
pub const LIFT_TOMIYAMA: u32 = 1;
/// Lift with a constant coefficient.
pub const LIFT_CONSTANT: u32 = 2;

/// No wall lubrication.
pub const WALL_NONE: u32 = 0;
/// Wall lubrication of Hosokawa et al. (2002).
pub const WALL_HOSOKAWA: u32 = 1;

/// σ_TD of the turbulent dispersion of Burns et al. (2004).
pub const SIGMA_TD: f64 = 0.9;

/// The constants of the bubbly mixture.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct Constants<T> {
    /// Liquid density ρ_l (kg/m³).
    pub rho_l: T,
    /// Gas density ρ_g (kg/m³).
    pub rho_g: T,
    /// Liquid viscosity μ_l (Pa s).
    pub mu_l: T,
    /// Gas viscosity μ_g (Pa s).
    pub mu_g: T,
    /// Surface tension σ (N/m).
    pub sigma: T,
    /// Bubble diameter d (m).
    pub d: T,
    /// Gravity g (m/s²).
    pub g: [T; 3],
    /// Lift model ([`LIFT_NONE`] …).
    pub lift: u32,
    /// C_L of [`LIFT_CONSTANT`].
    pub c_lift: T,
    /// 1: the lift is damped near walls ([`lift_wall_damping`]), 0: not.
    pub lift_damping: u32,
    /// Wall lubrication ([`WALL_NONE`] …).
    pub wall: u32,
    /// C_TD of the turbulent dispersion (0: none).
    pub c_td: T,
    /// C_μb of the bubble-induced viscosity (0: none).
    pub sato: T,
    /// The smallest phase fraction of the divisions (the residual
    /// fraction, 10⁻⁶).
    pub alpha_min: T,
}

#[inline(always)]
fn dot<T: Real>(a: [T; 3], b: [T; 3]) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
fn cross<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline(always)]
fn sub<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// x^p for x > 0 (0 for x ≤ 0).
#[inline(always)]
fn pow<T: Real>(x: T, p: f64) -> T {
    if x > T::ZERO {
        (x.ln_of() * T::from_f64(p)).exp_of()
    } else {
        T::ZERO
    }
}

/// |g| (m/s²).
#[inline(always)]
fn g_abs<T: Real>(c: &Constants<T>) -> T {
    dot(c.g, c.g).sqrt_of()
}

/// The Eötvös number g Δρ d²/σ of the diameter `d` (0 without σ).
#[inline(always)]
pub fn eotvos<T: Real>(c: &Constants<T>, d: T) -> T {
    if c.sigma > T::ZERO {
        g_abs(c) * (c.rho_l - c.rho_g) * d * d / c.sigma
    } else {
        T::ZERO
    }
}

/// C_D Re of Ishii & Zuber (1979) at the bubble Reynolds number `re` and
/// the Eötvös number `eo`, gas fraction `alpha` (the viscous sphere at
/// small Re, else the smaller of the distorted bubble and the cap).
#[inline(always)]
pub fn drag_cd_re<T: Real>(re: T, eo: T, alpha: T) -> T {
    let sphere = T::from_f64(24.0) * (T::ONE + T::from_f64(0.1) * pow(re, 0.75));
    let al = (T::ONE - alpha).max_of(T::from_f64(1e-3));
    let ellipse = T::from_f64(2.0 / 3.0) * eo.sqrt_of() * re;
    let cap = T::from_f64(8.0 / 3.0) * al * al * re;
    sphere.max_of(ellipse.min_of(cap))
}

/// The drag coefficient K (kg/(m³ s)) of gas fraction `alpha` (floored at
/// the residual fraction) and slip `ur` (m/s).
#[inline(always)]
pub fn drag_k<T: Real>(c: &Constants<T>, alpha: T, ur: T) -> T {
    let a = alpha.max_of(c.alpha_min).min_of(T::ONE - c.alpha_min);
    let mu_m = c.mu_l / (T::ONE - a);
    let re = c.rho_l * ur * c.d / mu_m;
    let eo = eotvos(c, c.d);
    T::from_f64(0.75) * drag_cd_re(re, eo, a) * mu_m * a / (c.d * c.d)
}

/// C_L of Tomiyama et al. (2002) at the bubble Reynolds number `re`.
#[inline(always)]
pub fn lift_tomiyama<T: Real>(c: &Constants<T>, re: T) -> T {
    let eo = eotvos(c, c.d);
    let dh = c.d * pow(T::ONE + T::from_f64(0.163) * pow(eo, 0.757), 1.0 / 3.0);
    let e = eotvos(c, dh);
    let f =
        T::from_f64(0.00105) * e * e * e - T::from_f64(0.0159) * e * e - T::from_f64(0.0204) * e
            + T::from_f64(0.474);
    if e < T::from_f64(4.0) {
        (T::from_f64(0.288) * (T::from_f64(0.121) * re).tanh_of()).min_of(f)
    } else if e <= T::from_f64(10.0) {
        f
    } else {
        T::from_f64(-0.27)
    }
}

/// The near-wall damping of the lift of Shaver & Podowski (2015) at the
/// wall distance `s` = y/d in bubble diameters: 0 for s < ½ (a bubble
/// whose centre is closer than its radius touches the wall, the shear
/// lift of the free stream does not apply), 3(2s − 1)² − 2(2s − 1)³
/// between ½ and 1 and 1 beyond. On wall-resolving grids the first cells
/// lie well inside a bubble radius, where the undamped lift with the wall
/// shear drives the bubbles into the wall at metres per second.
#[inline(always)]
pub fn lift_wall_damping<T: Real>(s: T) -> T {
    if !(s.to_f64().is_finite()) || s >= T::ONE {
        T::ONE
    } else if s <= T::from_f64(0.5) {
        T::ZERO
    } else {
        let x = T::from_f64(2.0) * s - T::ONE;
        T::from_f64(3.0) * x * x - T::from_f64(2.0) * x * x * x
    }
}

/// What [`cell_forces`] gives per cell.
#[derive(Clone, Copy, Debug, Default)]
pub struct CellForces<T> {
    /// K V (kg/s): the implicit drag of the rows.
    pub kv: T,
    /// (lift + wall lubrication) V on the gas (N).
    pub fv: [T; 3],
    /// μ_BIT of Sato (Pa s).
    pub mu_bit: T,
    /// The dispersion coefficient D = C_TD K (ν_t/σ_TD)(1/α_g + 1/α_l)
    /// (kg/(m s²)... per unit gradient of α: N/m³).
    pub disp: T,
    /// α* of the dispersion ([`alpha_dispersion`]), the α_g of the gas's
    /// [`cell_body`].
    pub alpha_td: T,
}

/// Inputs of [`cell_forces`].
#[derive(Clone, Copy)]
pub struct CellIn<'a, T> {
    /// Cell volumes.
    pub volume: &'a [T],
    /// α_g per cell.
    pub alpha: &'a [T],
    /// Liquid velocity.
    pub ul: Vec3<'a, T>,
    /// Gas velocity.
    pub ug: Vec3<'a, T>,
    /// The gradients of the three liquid velocity components (∇u_x, ∇u_y,
    /// ∇u_z).
    pub grad_ul: [Vec3<'a, T>; 3],
    /// Wall distance per cell (m).
    pub wall_y: &'a [T],
    /// Unit wall normal into the fluid per cell.
    pub wall_n: Vec3<'a, T>,
    /// ν_t of the liquid per cell (m²/s); empty: laminar.
    pub nu_t: &'a [T],
    /// The cell size h (m): the shortest distance to a neighbour centre.
    pub h: &'a [T],
    /// ∇α_g per cell.
    pub grad_alpha: Vec3<'a, T>,
}

/// The gas fraction α* of the turbulent dispersion of a cell: α_g, but at
/// least h |∇α_g| (and the residual fraction). The dispersion of Burns et
/// al. drifts the gas at (ν_t/σ_TD) ∇α_g/α_g, a diffusion velocity that
/// grows without bound where a cell is nearly empty beside a full one;
/// with α* it is at most about (ν_t/σ_TD)/h, the diffusion velocity of
/// the grid, and where α_g varies over several cells (α_g > h |∇α_g|) the
/// model is unchanged (D354).
#[inline(always)]
pub fn alpha_dispersion<T: Real>(alpha: T, h: T, grad: [T; 3], alpha_min: T) -> T {
    let g = dot(grad, grad).sqrt_of();
    alpha.max_of(h * g).max_of(alpha_min)
}

/// The interphase forces of cell `c`.
///
/// # Safety
/// The arrays have an entry for `c`.
#[inline(always)]
pub unsafe fn cell_forces<T: Real>(c: usize, inp: &CellIn<T>, k: &Constants<T>) -> CellForces<T> {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let vol = ld(inp.volume, c);
        let a = ld(inp.alpha, c).max_of(T::ZERO).min_of(T::ONE);
        let ul = inp.ul.at(c);
        let ug = inp.ug.at(c);
        let ur = sub(ug, ul);
        let urm = dot(ur, ur).sqrt_of();
        let kk = drag_k(k, a, urm);
        let re = k.rho_l * urm * k.d / k.mu_l;
        let mut f = [T::ZERO; 3];
        // lift
        let cl = match k.lift {
            LIFT_TOMIYAMA => lift_tomiyama(k, re),
            LIFT_CONSTANT => k.c_lift,
            _ => T::ZERO,
        };
        let cl = if k.lift_damping != 0 {
            cl * lift_wall_damping(ld(inp.wall_y, c) / k.d)
        } else {
            cl
        };
        if cl != T::ZERO {
            let g = [
                inp.grad_ul[0].at(c),
                inp.grad_ul[1].at(c),
                inp.grad_ul[2].at(c),
            ];
            // ω = ∇ × u (g[i][j] = ∂u_i/∂x_j)
            let w = [g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1]];
            let x = cross(ur, w);
            let s = -cl * k.rho_l * a;
            f = [s * x[0], s * x[1], s * x[2]];
        }
        // wall lubrication
        if k.wall == WALL_HOSOKAWA {
            // (a bubble centre is at least d/2 from the wall: closer cell
            // centres take the force at y = d/2, where (d/2y)² = 1; the
            // force of the formula there grows as 1/y² and drove the gas
            // of the wall cells at tens of m/s at a pipe inlet, D353)
            let y = ld(inp.wall_y, c).max_of(T::from_f64(0.5) * k.d);
            if y > T::ZERO && y.to_f64().is_finite() {
                let n = inp.wall_n.at(c);
                let un = dot(ur, n);
                let par = [ur[0] - un * n[0], ur[1] - un * n[1], ur[2] - un * n[2]];
                let u2 = dot(par, par);
                let eo = eotvos(k, k.d);
                let fw = if re > T::ZERO {
                    (T::from_f64(7.0) / pow(re, 1.9)).max_of(T::from_f64(0.0217) * eo)
                } else {
                    T::from_f64(0.0217) * eo
                };
                // (2/d) f (d/2y)² = f d/(2 y²)
                let s = fw * k.d / (T::from_f64(2.0) * y * y) * a * k.rho_l * u2;
                f = [f[0] + s * n[0], f[1] + s * n[1], f[2] + s * n[2]];
            }
        }
        let mu_bit = k.sato * k.rho_l * a * k.d * urm;
        let a_td = alpha_dispersion(a, ld(inp.h, c), inp.grad_alpha.at(c), k.alpha_min);
        let disp = if k.c_td > T::ZERO && !inp.nu_t.is_empty() {
            let ag = a_td;
            let al = (T::ONE - a).max_of(k.alpha_min);
            k.c_td * kk * ld(inp.nu_t, c) / T::from_f64(SIGMA_TD) * (T::ONE / ag + T::ONE / al)
        } else {
            T::ZERO
        };
        CellForces {
            kv: kk * vol,
            fv: [f[0] * vol, f[1] * vol, f[2] * vol],
            mu_bit,
            disp,
            alpha_td: a_td,
        }
    }
}

/// Partial elimination of the drag (Spalding 1980, D353) of one velocity
/// component of a cell. The rows were solved with the drag K V on the
/// diagonals a_g, a_l and the other phase lagged in the source: the gas
/// with the liquid `ul_lag`, the liquid with the solved gas `ug`. With
/// H_g = a_g u_g − K V u_l,lag and H_l = a_l u_l − K V u_g (everything of
/// the rows but the drag) the pair
///   a_g u_g − K V u_l = H_g,   a_l u_l − K V u_g = H_l
/// is solved exactly: the drag is implicit between the phases, and a
/// light phase follows the heavy one at once in place of over the outer
/// iterations. Returns (u_g, u_l); a singular pair (no drag) keeps the
/// solved values.
#[inline(always)]
pub fn eliminate<T: Real>(kv: T, ag: T, al: T, ug: T, ul: T, ul_lag: T) -> (T, T) {
    let det = ag * al - kv * kv;
    if det.partial_cmp(&T::ZERO) != Some(core::cmp::Ordering::Greater) || kv == T::ZERO {
        return (ug, ul);
    }
    let hg = ag * ug - kv * ul_lag;
    let hl = al * ul - kv * ug;
    ((al * hg + kv * hl) / det, (ag * hl + kv * hg) / det)
}

/// The phase coefficients D_k = ∂u_k/∂(−∇p) V of the eliminated pair
/// (D353) from the coefficients of the separate rows d_k = V/a_k (without
/// the phase fraction; `v` the cell volume): with a_k = V/d_k and
/// det = a_g a_l − (K V)²,
///   D_g = V (a_l α_g + K V α_l)/det,  D_l = V (a_g α_l + K V α_g)/det.
/// Without drag that is α_k d_k; a gas without volume (α_g = 0) gets the
/// liquid's response through the drag. Returns (D_g, D_l).
#[inline(always)]
pub fn eliminate_d<T: Real>(v: T, kv: T, alpha_g: T, dg: T, dl: T) -> (T, T) {
    let al_ = T::ONE - alpha_g;
    let (ag, al) = (v / dg, v / dl);
    let det = ag * al - kv * kv;
    if det.partial_cmp(&T::ZERO) != Some(core::cmp::Ordering::Greater) || kv == T::ZERO {
        return (alpha_g * dg, al_ * dl);
    }
    (
        v * (al * alpha_g + kv * al_) / det,
        v * (ag * al_ + kv * alpha_g) / det,
    )
}

/// The cell force per unit phase volume b_k of the gas (`gas`) or the
/// liquid: gravity (ρ_k − ρ_l) g and the dispersion ∓(D/α_k) ∇α_g
/// (`alpha`: α* of [`alpha_dispersion`] for the gas, α_g for the liquid).
#[inline(always)]
pub fn cell_body<T: Real>(k: &Constants<T>, gas: bool, alpha: T, disp: T, grad: [T; 3]) -> [T; 3] {
    if gas {
        let a = alpha.max_of(k.alpha_min);
        let s = disp / a;
        let r = k.rho_g - k.rho_l;
        [
            r * k.g[0] - s * grad[0],
            r * k.g[1] - s * grad[1],
            r * k.g[2] - s * grad[2],
        ]
    } else {
        let a = (T::ONE - alpha).max_of(k.alpha_min);
        let s = disp / a;
        [s * grad[0], s * grad[1], s * grad[2]]
    }
}

/// What [`face_arrays`] gives per face.
#[derive(Clone, Copy, Debug, Default)]
pub struct FaceArrays<T> {
    /// The face viscosity α_g,f (μ_g + μ_t,f ρ_g/ρ_l) of the gas.
    pub mu_g: T,
    /// α_l,f (μ_eff,f + μ_BIT,f) of the liquid.
    pub mu_l: T,
    /// The offset o_g,f/α_g,f of the gas.
    pub off_g: T,
    /// The offset o_l,f/α_l,f of the liquid.
    pub off_l: T,
    /// The conductance C_TD (ν_t/σ_TD)/α_l,f |S|²/(d·S) of the dispersion
    /// in the α equation (interior faces; 0 on the boundary and without
    /// dispersion): the drift −C_TD (ν_t/σ_TD)(1/α_g + 1/α_l) ∇α_g of the
    /// gas carries the volume C_TD (ν_t/σ_TD) ∇α_g/α_l (D354).
    pub td: T,
}

/// Inputs of [`face_arrays`].
#[derive(Clone, Copy)]
pub struct FaceIn<'a, T> {
    /// Connectivity.
    pub topo: &'a Topology<'a>,
    /// Interpolation weights of the owner.
    pub weight: &'a [T],
    /// Cell centroids.
    pub cell_centroid: &'a Vec3<'a, T>,
    /// Face centroids.
    pub face_centroid: &'a Vec3<'a, T>,
    /// The shift of the periodic faces (the neighbour centroid); empty.
    pub shift: &'a Vec3<'a, T>,
    /// α_g per cell.
    pub alpha: &'a [T],
    /// α_g of the inflow per boundary face.
    pub alpha_b: &'a [T],
    /// α_g per face (the face value of the gas flux).
    pub alpha_gf: &'a [T],
    /// α_l per face (the face value of the liquid flux).
    pub alpha_lf: &'a [T],
    /// μ_l + μ_t of the liquid per face (the turbulence model's face
    /// viscosity).
    pub mu_eff: &'a [T],
    /// μ_BIT per cell.
    pub mu_bit: &'a [T],
    /// The dispersion coefficient D per cell.
    pub disp: &'a [T],
    /// 1 per boundary face of an inflow (α_b given), 0 elsewhere (α_b =
    /// α_P).
    pub inflow: &'a [u32],
    /// |S|²/(d·S) per face.
    pub dgeo: &'a [T],
}

/// The face arrays of face `f`.
///
/// # Safety
/// The arrays have an entry for `f` and its cells.
#[inline(always)]
pub unsafe fn face_arrays<T: Real>(f: usize, inp: &FaceIn<T>, k: &Constants<T>) -> FaceArrays<T> {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let t = inp.topo;
        let p = ld(t.owner, f) as usize;
        let xp = inp.cell_centroid.at(p);
        let ap = ld(inp.alpha, p);
        let (d, an, w, nb) = if f < t.n_interior {
            let n = ld(t.neighbour, f) as usize;
            let xn = crate::fv::neighbour_centroid(f, n, t, inp.cell_centroid, inp.shift);
            (sub(xn, xp), ld(inp.alpha, n), ld(inp.weight, f), Some(n))
        } else {
            let fb = f - t.n_interior;
            let ab = if ld(inp.inflow, fb) != 0 {
                ld(inp.alpha_b, fb)
            } else {
                ap
            };
            (sub(inp.face_centroid.at(f), xp), ab, T::ONE, None)
        };
        let wn = T::ONE - w;
        let (bit, disp) = match nb {
            Some(n) => (
                w * ld(inp.mu_bit, p) + wn * ld(inp.mu_bit, n),
                w * ld(inp.disp, p) + wn * ld(inp.disp, n),
            ),
            None => (ld(inp.mu_bit, p), ld(inp.disp, p)),
        };
        let agf = ld(inp.alpha_gf, f);
        let alf = ld(inp.alpha_lf, f);
        let mu_eff = ld(inp.mu_eff, f);
        let mu_t = (mu_eff - k.mu_l).max_of(T::ZERO);
        let jump = an - ap;
        // (the face α* of the dispersion, as alpha_dispersion: at least the
        // jump of α across the face)
        let ag = agf.max_of(k.alpha_min);
        let ag_td = ag.max_of(jump.abs());
        let al = alf.max_of(k.alpha_min);
        let r = k.rho_g - k.rho_l;
        let td = if nb.is_some() && k.c_td > T::ZERO {
            // (1/α_l at most 10: the conductance only damps the step)
            let al = al.max_of(T::from_f64(0.1));
            k.c_td * mu_t / (k.rho_l * T::from_f64(SIGMA_TD) * al) * ld(inp.dgeo, f)
        } else {
            T::ZERO
        };
        FaceArrays {
            mu_g: ag * (k.mu_g + mu_t * k.rho_g / k.rho_l),
            mu_l: al * (mu_eff + bit),
            off_g: r * dot(k.g, d) - disp / ag_td * jump,
            off_l: disp / al * jump,
            td,
        }
    }
}
