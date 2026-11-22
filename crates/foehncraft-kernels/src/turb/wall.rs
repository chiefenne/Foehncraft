//! Wall treatments of SST per wall face (D134; stage 6 specification 4.1,
//! 4.3, 4.5).
//!
//! - **`low-re`** (integration to the wall): τ_w = μ U_t/y_P; the first
//!   cell takes ω_P = 6ν/(β₁ y_P²), the near-wall solution of the ω
//!   equation (Menter 1994), and k_w = 0.
//! - **`automatic`** (Menter, Kuntz & Langtry 2003): u_τ = (u_τ,vis⁴ +
//!   u_τ,log⁴)^¼ with u_τ,vis = (ν U_t/y_P)^½ and u_τ,log from the log law
//!   U_t = (u_τ/κ) ln(E y_P u_τ/ν) (fixed number of Newton steps), and
//!   ω_P = (ω_vis^n + ω_log^n)^{1/n}, ω_vis = 6ν/(β₁ y_P²),
//!   ω_log = u_τ/(√β* κ y_P), n = 1.2 (Menter et al.: n = 2); the wall
//!   shear τ_w = ρ u_τ² enters the momentum equation as the effective wall
//!   viscosity μ_w = τ_w y_P/U_t, k has zero flux through the wall. The
//!   near-wall cell takes ω_P, and the equilibrium production of the
//!   constant-stress layer P = (u_τ⁴/ν)(1 − g) g with g = du⁺/dy⁺ of
//!   Reichardt's law in place of its own, with the weight tanh((y⁺/5)²);
//!   below y⁺ ≈ 2 the ω equation with the wall-face value and the cell's
//!   production hold (integration to the wall). Stage 6g: the quartic
//!   weight u_τ,vis⁴/u_τ⁴ of the first version and the log-layer
//!   production u_τ³/(κ y_P) overpredicted C_f of the flat plate by up to
//!   17 % in the buffer layer; with these weights C_f lies within 2.6 % of
//!   our own `low-re` result on the finest plate mesh for first-cell y⁺
//!   from 0.24 to 85 on two meshes (D201 rework: the constants were first
//!   fitted to the value of another solver and are checked again against
//!   our own, `docs/verification/stage6.md`).
//! - **Rough walls** (Aupoix 2015, the set fitted to Colebrook's
//!   correlation): k_w = u_τ² k_w⁺ and ω_w = (u_τ²/ν) ω_w⁺ as wall values
//!   with
//!   k_w⁺ = max(0, (1/√β*) tanh[(ln(k_s⁺/30)/ln 10 + 1 − tanh(k_s⁺/125))
//!   tanh(k_s⁺/125)]),
//!   ω_w⁺ = 300/k_s⁺² [tanh(15/(4 k_s⁺))]⁻¹ + 191/k_s⁺ [1 − exp(−k_s⁺/250)],
//!   k_s⁺ = u_τ k_s/ν (B. Aupoix, *Roughness corrections for the k–ω shear
//!   stress transport model: status and proposals*, J. Fluids Eng. 137,
//!   021202 (2015); formulas checked on 2026-11-17 against their
//!   reproduction as eq. (11) of EUCASS2017-372, which names the source).
//!   With k_w > 0 the eddy viscosity at the wall is k_w/ω_w, so the wall
//!   shear is (μ + ρ ν_t,w) U_t/y_P.
//!
//! The k-ε family (stage 6 PR 6c; specification 4.3) has its own kernel,
//! [`keps_wall_face`], with the velocity scale u* = C_μ^¼ k_P^½ of the
//! near-wall cell in all three treatments (the equilibrium scale of
//! Launder & Spalding 1974): y* = u* y_P/ν, and the wall shear stress
//! τ_w/ρ = u* U_t/u⁺(y*), i.e. the effective wall viscosity
//! ν_w = ν y*/u⁺(y*). The viscous limit u⁺ = y* gives ν_w = ν exactly,
//! for any k; the log limit is the wall function of Launder & Spalding.
//! No iteration is needed, and u* stays finite where U_t vanishes
//! (separation and reattachment), unlike u_τ.
//!
//! - **`scalable`** (G. Grotjans, F. R. Menter, *Wall functions for general
//!   application CFD codes*, ECCOMAS 98, 1112–1117 (1998)): the log law at
//!   ỹ* = max(y*, 11.06), ε_P = u*³/(κ ỹ), P_k,P = (τ_w/ρ) u*/(κ ỹ) with
//!   ỹ = ỹ* ν/u*;
//! - **`standard`** (Launder & Spalding 1974): the log law for
//!   y* ≥ 11.06 and the viscous sublayer u⁺ = y* below (with P_k,P = 0
//!   there), ε_P = u*³/(κ y_P);
//! - **`compound`** (M. Popovac, K. Hanjalić, *Compound wall treatment for
//!   RANS computation of complex turbulent flows and heat transfer*, Flow
//!   Turbul. Combust. 78, 177–202 (2007)): u⁺ = e^{−Γ} y* + e^{−1/Γ} u⁺_log,
//!   Γ = 0.01 y*⁴/(1 + 5 y*) (the blending of Kader 1981 with Pr = 1),
//!   ε_P = e^{−Γ} 2νk_P/y_P² + e^{−1/Γ} u*³/(κ y_P), and the production
//!   of the cell e^{−Γ} P_k + e^{−1/Γ} (τ_w/ρ) u*/(κ y_P). Γ is evaluated
//!   with y* (the specification writes y⁺; the two agree in the log layer,
//!   and in the sublayer, where y* < y⁺, the viscous branch is exact).
//! - **Rough walls** (all three): u⁺_log = ln(E y*)/κ − ΔB(k_s⁺),
//!   k_s⁺ = u* k_s/ν, with ΔB of T. Cebeci, P. Bradshaw, *Momentum
//!   transfer in boundary layers*, Hemisphere (1977): 0 for k_s⁺ < 2.25,
//!   (1/κ) ln[(k_s⁺ − 2.25)/87.75 + C_s k_s⁺] sin[0.4258 (ln k_s⁺ − 0.811)]
//!   up to 90, (1/κ) ln(1 + C_s k_s⁺) above (the pieces join
//!   continuously). The log-law u⁺ is kept at least 1 (a roughness height
//!   above the first-cell centre would make it vanish; BC027 warns).

use super::sst::SstConstants;
use foehncraft_core::Real;

/// E of the log law u⁺ = (1/κ) ln(E y⁺) (B = 5.2 with κ = 0.41:
/// E = e^{κB} = 8.43; E = 9.793 is the value of the specification 4.2
/// with B ≈ 5.56). The specification fixes E = 9.793.
pub const LOG_LAW_E: f64 = 9.793;

/// Newton steps of the log law (fixed: no data-dependent loop on the GPU).
pub const NEWTON_STEPS: usize = 12;

/// Wall treatment: integration to the wall.
pub const LOW_RE: u32 = 0;
/// Wall treatment: the automatic blending of Menter et al. (2003).
pub const AUTOMATIC: u32 = 1;
/// Wall treatment (k-ε family): scalable wall functions (Grotjans & Menter
/// 1998).
pub const SCALABLE: u32 = 2;
/// Wall treatment (k-ε family): standard wall functions (Launder &
/// Spalding 1974).
pub const STANDARD: u32 = 3;
/// Wall treatment (k-ε family): the compound treatment (Popovac &
/// Hanjalić 2007).
pub const COMPOUND: u32 = 4;
/// Wall treatment (k-ε family): the compound treatment of Popovac (2023),
/// the harmonic blend of the viscous and the logarithmic layer with the
/// non-equilibrium factor ψ.
pub const COMPOUND_2023: u32 = 5;

/// E of the log law of the 2023 compound treatment (Popovac 2023,
/// eq. (11): κ = 0.41, E = 8.9).
pub const E_2023: f64 = 8.9;
/// The exponent n_b = 10/3 of the harmonic blend (Popovac 2023, eq. (17)).
pub const N_BLEND: f64 = 10.0 / 3.0;
/// The viscous sublayer thickness y_v⁺ = 11 of the non-equilibrium factor
/// (Popovac 2023, eq. (15)).
pub const Y_V_PLUS: f64 = 11.0;
/// Bounds of the non-equilibrium factor ψ (Popovac 2023 bounds ψ but
/// gives no values; ours, 6d2).
pub const PSI_MIN: f64 = 0.5;
/// See [`PSI_MIN`].
pub const PSI_MAX: f64 = 2.0;

/// The lower end of the log layer, y⁺ = 11.06 (where u⁺ = y⁺ meets the log
/// law with κ = 0.41, E = 9.793; the limit of the scalable wall functions
/// of Grotjans & Menter 1998).
pub const Y_PLUS_LOG: f64 = 11.06;

/// u_τ of the log law U = (u_τ/κ) ln(E max(y u_τ/ν, 11.06)), by Newton's
/// method from u₀ = κ U.
///
/// The floor of y⁺ at the lower end of the log layer makes u_τ,log vanish
/// in the viscous sublayer (u_τ,log = κU/ln(11.06 E) ≈ 0.088 U there,
/// against u_τ,vis = (νU/y)^½ ≫ U for y⁺ < 1), so that the quartic blend
/// of the automatic treatment tends to the viscous value; without it the
/// log law has a spurious root at E y⁺ ≈ 1 whenever the cell lies in the
/// sublayer. The function g(u) = (u/κ) ln(E max(y u/ν, 11.06)) − U is
/// convex and increasing, and g(κU) ≥ U (ln(11.06 E) − 1) > 0, so the
/// iterates fall monotonically to the root.
#[inline(always)]
pub fn log_law_u_tau<T: Real>(u: T, y: T, nu: T, kappa: T) -> T {
    let e = T::from_f64(LOG_LAW_E);
    let y_nu = y / nu;
    let floor = T::from_f64(Y_PLUS_LOG);
    let mut ut = kappa * u;
    let mut i = 0;
    while i < NEWTON_STEPS {
        let yp = y_nu * ut;
        let lin = yp < floor;
        let l = (e * yp.max_of(floor)).ln_of();
        let g = ut * l / kappa - u;
        // d/du [u ln(E y u/ν)] = ln + 1 in the log layer, ln(E 11.06) below
        let dg = if lin { l / kappa } else { (l + T::ONE) / kappa };
        ut -= g / dg;
        i += 1;
    }
    ut
}

/// Inputs of [`wall_face`] for one wall face.
#[derive(Clone, Copy, Debug)]
pub struct WallInputs<T> {
    /// Tangential speed of the near-wall cell relative to the wall (m/s).
    pub ut: T,
    /// Normal distance y_P of the cell centroid from the face (m).
    pub y: T,
    /// ν (m²/s).
    pub nu: T,
    /// Equivalent sand-grain roughness k_s (m), 0 for a smooth wall.
    pub ks: T,
    /// u_τ of the last outer iteration (rough walls: k_s⁺).
    pub u_tau_old: T,
}

/// Results of [`wall_face`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallFace<T> {
    /// u_τ (m/s).
    pub u_tau: T,
    /// y⁺ = u_τ y_P/ν.
    pub y_plus: T,
    /// ω of the near-wall cell (automatic treatment: weighted by
    /// 1 − `blending` in its row) (1/s).
    pub omega: T,
    /// ω at the wall face (Dirichlet value): smooth walls
    /// 10·6ν/(β₁ d₁²) with d₁ = 2 y_P, the wall-cell height (Menter 1994);
    /// rough walls the value of Aupoix (1/s).
    pub omega_wall: T,
    /// k at the wall (rough walls; 0 smooth).
    pub k_wall: T,
    /// Effective kinematic wall viscosity ν_w = τ_w y_P/(ρ U_t) (m²/s).
    pub nu_wall: T,
    /// The viscous weight of the near-wall ω (1: the ω equation with the
    /// wall-face value, integration to the wall; automatic treatment
    /// 1 − tanh((y⁺/5)²), the cell takes ω_P with the weight 1 − `blending`).
    pub blending: T,
    /// The equilibrium production of k at the near-wall cell (automatic
    /// treatment), τ_w/ρ · dU/dy of the constant-stress layer with the
    /// slope g = du⁺/dy⁺ of Reichardt's law: (u_τ⁴/ν)(1 − g) g, which tends
    /// to 0 in the viscous sublayer and to u_τ³/(κ y_P) in the log layer
    /// (m²/s³).
    pub pk_wall: T,
    /// The weight of `pk_wall` in the production of the near-wall cell
    /// (automatic treatment tanh((y⁺/5)²), the weight of ω_P: above the
    /// sublayer it replaces the cell's own production, whose velocity
    /// gradient is too large there; 0 otherwise).
    pub prod_weight: T,
}

/// k_w⁺ of Aupoix (2015), Colebrook set.
#[inline(always)]
pub fn aupoix_k_plus<T: Real>(ks_plus: T, beta_star: T) -> T {
    let f = T::from_f64;
    let t = (ks_plus / f(125.0)).tanh_of();
    let a = (ks_plus / f(30.0)).ln_of() / f(core::f64::consts::LN_10) + T::ONE - t;
    ((a * t).tanh_of() / beta_star.sqrt_of()).max_of(T::ZERO)
}

/// ω_w⁺ of Aupoix (2015), Colebrook set.
#[inline(always)]
pub fn aupoix_omega_plus<T: Real>(ks_plus: T) -> T {
    let f = T::from_f64;
    let ks2 = ks_plus * ks_plus;
    f(300.0) / ks2 / (f(15.0) / (f(4.0) * ks_plus)).tanh_of()
        + f(191.0) / ks_plus * (T::ONE - (-ks_plus / f(250.0)).exp_of())
}

/// The wall treatment `treatment` ([`LOW_RE`] or [`AUTOMATIC`]) for one
/// wall face.
#[inline(always)]
pub fn wall_face<T: Real>(i: &WallInputs<T>, c: &SstConstants, treatment: u32) -> WallFace<T> {
    let f = T::from_f64;
    let (ut, y, nu) = (i.ut, i.y, i.nu);
    let kappa = f(c.kappa);
    let w_vis = f(6.0) * nu / (f(c.beta1) * y * y);
    // 10·6ν/(β₁ (2 y_P)²)
    let omega_wall = f(2.5) * w_vis;
    let tiny = f(1e-300);
    if i.ks > T::ZERO {
        // rough wall: wall values of Aupoix with u_τ of the last iteration
        // (at least that of the viscous wall shear)
        let u_vis = (nu * ut / y).sqrt_of();
        let u_tau = i.u_tau_old.max_of(u_vis);
        let ksp = (u_tau * i.ks / nu).max_of(tiny);
        let u2 = u_tau * u_tau;
        let k_wall = u2 * aupoix_k_plus(ksp, f(c.beta_star));
        let omega = u2 / nu * aupoix_omega_plus(ksp);
        return WallFace {
            u_tau,
            y_plus: u_tau * y / nu,
            omega,
            omega_wall: omega,
            k_wall,
            nu_wall: nu + k_wall / omega,
            blending: T::ONE,
            pk_wall: T::ZERO,
            prod_weight: T::ZERO,
        };
    }
    if treatment == LOW_RE {
        let u_tau = (nu * ut / y).sqrt_of();
        return WallFace {
            u_tau,
            y_plus: u_tau * y / nu,
            omega: w_vis,
            omega_wall,
            k_wall: T::ZERO,
            nu_wall: nu,
            blending: T::ONE,
            pk_wall: T::ZERO,
            prod_weight: T::ZERO,
        };
    }
    let u_vis2 = nu * ut / y;
    let u_log = log_law_u_tau(ut, y, nu, kappa);
    let v4 = u_vis2 * u_vis2;
    let l2 = u_log * u_log;
    let l4 = l2 * l2;
    let u4 = v4 + l4;
    let u_tau = u4.sqrt_of().sqrt_of();
    let w_log = u_tau / (f(c.beta_star).sqrt_of() * kappa * y);
    // (ω_vis^n + ω_log^n)^{1/n}, n = 1.2 (ω > 0: as exp/ln)
    let n = f(OMEGA_BLEND_N);
    let omega =
        (((w_vis.ln_of() * n).exp_of() + (w_log.ln_of() * n).exp_of()).ln_of() / n).exp_of();
    let y_plus = u_tau * y / nu;
    // the weight of the cell values: ω_P and the equilibrium production
    let q = y_plus / f(OMEGA_FIX_Y_PLUS);
    let fix = (q * q).tanh_of();
    // ν_w = u_τ² y/U_t (U_t > 0; ν in a fluid at rest)
    let nu_wall = if ut > T::ZERO {
        u4.sqrt_of() * y / ut
    } else {
        nu
    };
    WallFace {
        u_tau,
        y_plus,
        omega,
        omega_wall,
        k_wall: T::ZERO,
        nu_wall,
        blending: T::ONE - fix,
        pk_wall: {
            let g = reichardt_slope(y_plus, kappa);
            u4 / nu * (T::ONE - g).max_of(T::ZERO) * g
        },
        prod_weight: fix,
    }
}

/// The exponent n of the near-wall ω of the automatic treatment,
/// (ω_vis^n + ω_log^n)^{1/n} (6g: 2 in Menter et al. 2003; 1.2 brings
/// C_f of the flat plate within 2.5 % of our own `low-re` result for
/// first-cell y⁺ from 0.8 to 85; rechecked against our own result, the
/// optimum of n lies between 1.2 and 1.3).
pub const OMEGA_BLEND_N: f64 = 1.2;

/// The y⁺ scale of the weight tanh((y⁺/y₀)²) with which the near-wall cell
/// takes ω_P in the automatic treatment (6g; below y⁺ ≈ 2 the ω equation
/// with the wall-face value, integration to the wall).
pub const OMEGA_FIX_Y_PLUS: f64 = 5.0;

/// du⁺/dy⁺ of Reichardt's law u⁺ = ln(1 + κy⁺)/κ + 7.8 (1 − e^{−y⁺/11} −
/// (y⁺/11) e^{−y⁺/3}).
#[inline(always)]
pub fn reichardt_slope<T: Real>(yp: T, kappa: T) -> T {
    let f = T::from_f64;
    let e11 = (-yp / f(11.0)).exp_of();
    let e3 = (-yp / f(3.0)).exp_of();
    T::ONE / (T::ONE + kappa * yp) + f(7.8) * (e11 / f(11.0) - e3 / f(11.0) + yp / f(33.0) * e3)
}

/// ΔB(k_s⁺) of Cebeci & Bradshaw (1977) with the roughness constant
/// C_s (0 below k_s⁺ = 2.25, transitional up to 90, fully rough above).
#[inline(always)]
pub fn roughness_shift<T: Real>(ks_plus: T, cs: T, kappa: T) -> T {
    let f = T::from_f64;
    let lo = f(2.25);
    let hi = f(90.0);
    let k = ks_plus.max_of(lo);
    // transitional
    let tr = ((k - lo) / f(87.75) + cs * k).ln_of() * (f(0.4258) * (k.ln_of() - f(0.811))).sin_of();
    // fully rough
    let fr = (T::ONE + cs * k).ln_of();
    let v = if ks_plus < lo {
        T::ZERO
    } else if ks_plus < hi {
        tr
    } else {
        fr
    };
    v / kappa
}

/// The wall values of the k-ε-ζ-f model (Hanjalić, Popovac & Hadžiabdić
/// 2004) at a wall face from its near-wall cell: ε_w = 2ν k_P/y_P² and
/// f_w = −2ν ζ_P/y_P² (the limits for y → 0 with the first-cell values,
/// k ~ y² and ζ ~ y² near the wall). With the compound treatment
/// (Popovac & Hanjalić 2007) f_w is blended towards 0 with the viscous
/// weight e^{−Γ} (`viscous`, 1 for the integration to the wall). The
/// ζ-f₀ form (Popovac 2023, eq. (21)) generalises f_w to
/// f₀_w = (N_ζ − 2) 2ν ζ_P/y_P², 0 with its N_ζ = 2; `n_zeta` = 1 is the
/// 2004 model. The elliptic-blending model (Billard & Laurence 2012)
/// solves ε with half the molecular diffusion and takes ε_w = ν k_P/y_P²
/// (`e_wall` = 1 instead of 2); its α_w = 0 is the f₀_w with N_ζ = 2.
/// Returns (ε_w, f_w).
#[inline(always)]
pub fn zeta_f_wall_values<T: Real>(
    k: T,
    zeta: T,
    y: T,
    nu: T,
    viscous: T,
    n_zeta: T,
    e_wall: T,
) -> (T, T) {
    let nu_y2 = nu / (y * y);
    (
        e_wall * nu_y2 * k.max_of(T::ZERO),
        (n_zeta - T::from_f64(2.0)) * viscous * T::from_f64(2.0) * nu_y2 * zeta.max_of(T::ZERO),
    )
}

/// Γ = 0.01 y⁴/(1 + 5 y) of the compound blending (Kader 1981, Pr = 1).
#[inline(always)]
pub fn compound_gamma<T: Real>(y: T) -> T {
    let y2 = y * y;
    T::from_f64(0.01) * y2 * y2 / (T::ONE + T::from_f64(5.0) * y)
}

/// Inputs of [`keps_wall_face`] for one wall face.
#[derive(Clone, Copy, Debug)]
pub struct KepsWallInputs<T> {
    /// Tangential speed of the near-wall cell relative to the wall (m/s).
    pub ut: T,
    /// Normal distance y_P of the cell centroid from the face (m).
    pub y: T,
    /// ν (m²/s).
    pub nu: T,
    /// k of the near-wall cell (m²/s²).
    pub k: T,
    /// Equivalent sand-grain roughness k_s (m), 0 for a smooth wall.
    pub ks: T,
    /// C_s of Cebeci & Bradshaw.
    pub cs: T,
    /// The non-equilibrium term C_U = t·[(U·∇)U + ∇p/ρ] of the near-wall
    /// cell along the wall-tangential flow direction t (m/s²), used by
    /// [`COMPOUND_2023`] only (0: equilibrium).
    pub cu: T,
}

/// Results of [`keps_wall_face`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KepsWallFace<T> {
    /// u_τ = (τ_w/ρ)^½ (m/s).
    pub u_tau: T,
    /// y⁺ = u_τ y_P/ν.
    pub y_plus: T,
    /// y* = C_μ^¼ k_P^½ y_P/ν.
    pub y_star: T,
    /// ε of the near-wall cell (fixed value) (m²/s³).
    pub epsilon: T,
    /// Effective kinematic wall viscosity ν_w = τ_w y_P/(ρ U_t) (m²/s).
    pub nu_wall: T,
    /// Weight m_p of the wall production: the cell takes
    /// (1 − m_p) P_k + `pk_wall`.
    pub prod_weight: T,
    /// The wall part of the production (already weighted) (m²/s³).
    pub pk_wall: T,
    /// The viscous weight e^{−Γ} of the compound blending (1 in the
    /// sublayer of the standard wall function, else 0; the blended
    /// dU⁺/dy⁺ with [`COMPOUND_2023`]).
    pub viscous: T,
    /// The non-equilibrium factor ψ ([`COMPOUND_2023`], else 1).
    pub psi: T,
}

/// The wall treatment of Popovac (2023) for one wall face, in the
/// velocity scale u* = C_μ^¼ k_P^½ of the k-ε family (y* = u* y_P/ν).
///
/// U⁺ = [(U⁺_vis)^{−n} + (U⁺_log)^{−n}]^{−1/n}, n = 10/3 (eq. (17)), with
/// U⁺_vis = y*, U⁺_log = ln(E y*)/(κψ) (eq. (15)), E = 8.9. The log value
/// is taken at max(y*, 1): below y* ≈ 0.12 ln(E y*) would vanish and the
/// harmonic blend would follow it to 0; with the guard it is ≥ 5.3 there,
/// and U⁺ = y* to 0.1 % for y* ≤ 1.
///
/// ψ = 1 − C_U (y_P − y_v/2)/(U_t κ u*), y_v = 11 ν/u*, bounded to
/// [ψ_min, ψ_max] (eq. (15)). The paper applies it for y_P > y_v (eq. (14));
/// it is applied here from y_v/2 on, where it is 1, so that ψ is
/// continuous in y_P.
///
/// Returns (y*/U⁺, Γ = dU⁺/dy⁺, ψ): the effective wall viscosity is
/// ν_w = ν y*/U⁺ (eq. (16)); Γ (eq. (19), derived from eq. (17)) =
/// (1 + q^n)^{−(n+1)/n} [1 + q^{n+1}/(κψ y*)], q = y*/U⁺_log, tends to 1
/// in the sublayer and to 1/(κψ y*) in the log layer.
#[inline(always)]
pub fn compound_2023<T: Real>(ys: T, y: T, nu: T, u_star: T, ut: T, cu: T, kappa: T) -> (T, T, T) {
    let f = T::from_f64;
    let tiny = f(1e-300);
    let n = f(N_BLEND);
    // ψ (1 without a velocity scale or a tangential velocity)
    let yv = f(Y_V_PLUS) * nu / u_star.max_of(tiny);
    let den = ut * kappa * u_star;
    let psi = if den > tiny {
        (T::ONE - cu * (y - f(0.5) * yv).max_of(T::ZERO) / den)
            .max_of(f(PSI_MIN))
            .min_of(f(PSI_MAX))
    } else {
        T::ONE
    };
    let yl = ys.max_of(T::ONE);
    let u_log = (f(E_2023) * yl).ln_of() / (kappa * psi);
    let q = (ys / u_log).max_of(tiny);
    let lq = q.ln_of();
    // 1 + q^n
    let s = T::ONE + (n * lq).exp_of();
    let ls = s.ln_of();
    // y*/U⁺ = (1 + q^n)^{1/n}
    let r = (ls / n).exp_of();
    // Γ = (1 + q^n)^{−(n+1)/n} [1 + q^{n+1}/(κψ y*)]
    let tail = ((n + T::ONE) * lq).exp_of() / (kappa * psi * yl);
    let g = (-(n + T::ONE) / n * ls).exp_of() * (T::ONE + tail);
    (r, g.min_of(T::ONE), psi)
}

/// The wall treatment `treatment` ([`SCALABLE`], [`STANDARD`] or
/// [`COMPOUND`]) of the k-ε family for one wall face.
#[inline(always)]
pub fn keps_wall_face<T: Real>(
    i: &KepsWallInputs<T>,
    c: &super::keps::KepsConstants,
    treatment: u32,
) -> KepsWallFace<T> {
    let f = T::from_f64;
    let (ut, y, nu) = (i.ut, i.y, i.nu);
    let kappa = f(c.kappa);
    let tiny = f(1e-300);
    let k = i.k.max_of(T::ZERO);
    let u_star = f(foehncraft_core::math::sqrt(foehncraft_core::math::sqrt(
        c.c_mu,
    ))) * k.sqrt_of();
    let ys = u_star * y / nu;
    let shift = if i.ks > T::ZERO {
        roughness_shift(u_star * i.ks / nu, i.cs, kappa)
    } else {
        T::ZERO
    };
    let e = f(LOG_LAW_E);
    let floor = f(Y_PLUS_LOG);
    let log_plus = |yy: T| ((e * yy).ln_of() / kappa - shift).max_of(T::ONE);
    let u3 = u_star * u_star * u_star;
    // ratio r = y*/u⁺ (ν_w = ν r), the length of the log-layer values and
    // the weights
    let mut psi = T::ONE;
    let (r, y_len, viscous, w_log, eps_vis) = if treatment == COMPOUND_2023 {
        // Popovac (2023): ε_P = Γ 2νk/y² + (1 − Γ) u*³/(κ y) (eq. (20))
        let (r, g, ps) = compound_2023(ys, y, nu, u_star, ut, i.cu, kappa);
        psi = ps;
        (r, y, g, T::ONE - g, g * f(c.e_wall) * nu * k / (y * y))
    } else if treatment == COMPOUND {
        let g = compound_gamma(ys).max_of(tiny);
        let a = (-g).exp_of();
        let b = (-(T::ONE / g)).exp_of();
        let yl = ys.max_of(T::ONE);
        // u⁺/y* = a + b u⁺_log/y* (finite as y* → 0, where b = 0)
        let r = T::ONE / (a + b * log_plus(yl) / ys.max_of(tiny));
        (r, y, a, b, a * f(c.e_wall) * nu * k / (y * y))
    } else if treatment == STANDARD {
        let log = ys >= floor;
        if log {
            (ys / log_plus(ys), y, T::ZERO, T::ONE, T::ZERO)
        } else {
            (T::ONE, y, T::ONE, T::ZERO, T::ZERO)
        }
    } else {
        // scalable: the log law at ỹ* = max(y*, 11.06), ỹ = ỹ* ν/u*
        let yt = ys.max_of(floor);
        let yl = y.max_of(floor * nu / u_star.max_of(tiny));
        (ys / log_plus(yt), yl, T::ZERO, T::ONE, T::ZERO)
    };
    let nu_wall = nu * r;
    // τ_w/ρ = ν_w U_t/y_P
    let tau = nu_wall * ut / y;
    let u_tau = tau.sqrt_of();
    let eps = eps_vis + w_log * u3 / (kappa * y_len);
    let pk_log = tau * u_star / (kappa * y_len);
    // standard: no production in the sublayer (w_log = 0 there)
    let (prod_weight, pk_wall) = if treatment == COMPOUND_2023 {
        // Popovac (2023): P = |uv| dU/dy fixed in the cell, |uv| = u*²,
        // dU/dy = (τ_w/ρ)/ν Γ (Γ = dU⁺/dy⁺; the log limit is
        // τ_w u*/(ρ κ y), the production of the other treatments)
        (T::ONE, u_star * u_star * tau / nu * viscous)
    } else if treatment == COMPOUND {
        (T::ONE - viscous, w_log * pk_log)
    } else {
        (T::ONE, w_log * pk_log)
    };
    KepsWallFace {
        u_tau,
        y_plus: u_tau * y / nu,
        y_star: ys,
        epsilon: eps.max_of(tiny),
        nu_wall,
        prod_weight,
        pk_wall,
        viscous,
        psi,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_law_newton_converges_over_the_y_plus_range() {
        let (nu, kappa) = (1e-5, 0.41);
        for &yp in &[11.06, 11.5, 30.0, 100.0, 1e3, 1e4, 1e5] {
            // build U from a chosen u_τ = 1 and y = y⁺ ν
            let y = yp * nu;
            let u = (1.0 / kappa) * foehncraft_core::math::ln(LOG_LAW_E * yp);
            let ut = log_law_u_tau(u, y, nu, kappa);
            assert!((ut - 1.0).abs() < 1e-12, "y+ {yp}: {ut}");
        }
    }

    #[test]
    fn automatic_treatment_tends_to_its_limits() {
        let c = SstConstants::default();
        let nu = 1e-5;
        // viscous sublayer: y⁺ = 0.1 with u_τ = 1: U = y⁺ u_τ
        let y = 0.1 * nu;
        let i = WallInputs {
            ut: 0.1,
            y,
            nu,
            ks: 0.0,
            u_tau_old: 0.0,
        };
        let r = wall_face(&i, &c, AUTOMATIC);
        assert!((r.u_tau - 1.0).abs() < 1e-6, "{}", r.u_tau);
        assert!((r.nu_wall - nu).abs() < 1e-5 * nu);
        let w_vis = 6.0 * nu / (0.075 * y * y);
        // the log part adds (ω_log/ω_vis)^1.2/1.2 relative, ω_log/ω_vis =
        // 0.075 y⁺/(6 · 0.3 · 0.41) = 0.010 here
        assert!(
            (r.omega - w_vis).abs() < 5e-3 * w_vis,
            "{} {w_vis}",
            r.omega
        );
        assert!(r.blending > 0.999);
        // the equilibrium production vanishes in the sublayer, (1 − g) g
        // with g = 1 − O(y⁺²)
        assert!(
            r.pk_wall < 0.02 * (r.u_tau * r.u_tau * r.u_tau * r.u_tau) / nu,
            "{}",
            r.pk_wall
        );
        assert!(r.prod_weight < 1e-3);
        // log layer: y⁺ = 300 on the log law
        let y = 300.0 * nu;
        let u = foehncraft_core::math::ln(LOG_LAW_E * 300.0) / 0.41;
        let i = WallInputs {
            ut: u,
            y,
            nu,
            ks: 0.0,
            u_tau_old: 0.0,
        };
        let r = wall_face(&i, &c, AUTOMATIC);
        // the viscous part adds (νU/y)²/4 relative: U/300 = 0.048, +0.06 %
        let v = u / 300.0;
        let blended = foehncraft_core::math::sqrt(foehncraft_core::math::sqrt(1.0 + v * v));
        assert!((r.u_tau - blended).abs() < 1e-12, "{}", r.u_tau);
        let w_log = r.u_tau / (0.3 * 0.41 * y);
        let w_vis = 6.0 * nu / (0.075 * y * y);
        let n = OMEGA_BLEND_N;
        let w = foehncraft_core::math::pow(
            foehncraft_core::math::pow(w_vis, n) + foehncraft_core::math::pow(w_log, n),
            1.0 / n,
        );
        assert!((r.omega - w).abs() < 1e-12 * w, "{} {w}", r.omega);
        // ω_vis/ω_log = 0.033 at y⁺ 300: within 2 % of the log value
        assert!((r.omega - w_log).abs() < 2e-2 * w_log);
        // the cell takes ω_P fully, and the production tends to the log
        // layer's u_τ³/(κ y): g = 1/(κ y⁺) (1 + O(1/(κy⁺))) there
        assert!(r.blending < 1e-12);
        let pk_log = (r.u_tau * r.u_tau * r.u_tau) / (0.41 * y);
        assert!(
            (r.pk_wall - pk_log).abs() < 2e-2 * pk_log,
            "{} {pk_log}",
            r.pk_wall
        );
        assert!((r.y_plus - 300.0 * r.u_tau).abs() < 1e-9);
        // low-re: the viscous values
        let r = wall_face(&i, &c, LOW_RE);
        assert_eq!(r.nu_wall, nu);
    }

    /// The near-wall ω of the automatic treatment moves from the ω
    /// equation (y⁺ ≲ 2) to the cell value (y⁺ ≳ 10) with tanh((y⁺/5)²), and
    /// the equilibrium production of the constant-stress layer has its
    /// maximum u_τ⁴/(4ν) where du⁺/dy⁺ = ½, near y⁺ 11 (Reichardt).
    #[test]
    fn automatic_buffer_layer_weights() {
        let c = SstConstants::default();
        let nu = 1e-5;
        let (mut best, mut at) = (0.0, 0.0);
        for j in 1..400 {
            let yp = 0.1 * j as f64;
            let g = reichardt_slope(yp, 0.41);
            let p = (1.0 - g).max(0.0) * g;
            if p > best {
                (best, at) = (p, yp);
            }
        }
        assert!((best - 0.25).abs() < 1e-3, "{best}");
        assert!((9.0..14.0).contains(&at), "{at}");
        // the slope of Reichardt's law is 1 at the wall, 1/(κ y⁺) far out
        assert!((reichardt_slope(0.0, 0.41) - 1.0).abs() < 1e-14);
        let g = reichardt_slope(1000.0, 0.41);
        assert!((g * 0.41 * 1000.0 - 1.0).abs() < 3e-3, "{g}");
        // the weight of ω_P against y⁺ of the face
        for (yp, lo, hi) in [(1.0, 0.95, 1.0), (5.0, 0.2, 0.3), (10.0, 0.0, 1e-3)] {
            let y = yp * nu; // u_τ ≈ 1 on Reichardt's law
            let u = (1.0 / 0.41) * foehncraft_core::math::ln(1.0 + 0.41 * yp)
                + 7.8
                    * (1.0
                        - foehncraft_core::math::exp(-yp / 11.0)
                        - yp / 11.0 * foehncraft_core::math::exp(-yp / 3.0));
            let r = wall_face(
                &WallInputs {
                    ut: u,
                    y,
                    nu,
                    ks: 0.0,
                    u_tau_old: 0.0,
                },
                &c,
                AUTOMATIC,
            );
            let q = r.y_plus / OMEGA_FIX_Y_PLUS;
            assert!((r.blending - (1.0 - foehncraft_core::math::tanh(q * q))).abs() < 1e-15);
            assert!(
                r.blending >= lo && r.blending <= hi,
                "y+ {yp}: {}",
                r.blending
            );
        }
    }

    #[test]
    fn aupoix_values_of_fully_rough_and_smooth_limits() {
        // ω_w⁺ ~ 300/k_s⁺² · 4k_s⁺/15 ... for small k_s⁺ the smooth limit
        // grows like k_s⁺⁻²; fully rough: ω_w⁺ → 191/k_s⁺ + 300/k_s⁺²·
        let w = aupoix_omega_plus(1000.0_f64);
        let exact = 300.0 / 1e6 / foehncraft_core::math::tanh(15.0 / 4000.0)
            + 0.191 * (1.0 - foehncraft_core::math::exp(-4.0));
        assert!((w - exact).abs() < 1e-14);
        // k_w⁺ = 0 below the onset of roughness effects (k_s⁺ small)
        assert_eq!(aupoix_k_plus(1.0_f64, 0.09), 0.0);
        // fully rough: tanh(k_s⁺/125) = 1, k_w⁺ = tanh(log₁₀(k_s⁺/30))/√β*
        let k = aupoix_k_plus(1e4_f64, 0.09);
        use foehncraft_core::math as m;
        let exact = m::tanh(m::log10(1e4 / 30.0)) / 0.3;
        assert!((k - exact).abs() < 1e-12, "{k} {exact}");
    }

    #[test]
    fn keps_wall_functions_in_the_log_layer() {
        // equilibrium log layer: k = u_τ²/√C_μ, U = u_τ u⁺_log(y⁺): every
        // treatment returns u_τ, ε = u_τ³/(κ y), P = u_τ³/(κ y)
        let c = crate::turb::keps::KepsConstants::standard();
        let nu = 1e-5;
        for &yp in &[30.0, 100.0, 300.0] {
            let y = yp * nu;
            let u = foehncraft_core::math::ln(LOG_LAW_E * yp) / 0.41;
            let i = KepsWallInputs {
                ut: u,
                y,
                nu,
                k: 1.0 / 0.3,
                ks: 0.0,
                cs: 0.5,
                cu: 0.0,
            };
            for t in [SCALABLE, STANDARD, COMPOUND] {
                let r = keps_wall_face(&i, &c, t);
                // the compound blend reaches the log law slowly:
                // e^{−1/Γ} = 0.98 at y⁺ = 30
                let tol = if t == COMPOUND { 2.5e-2 } else { 1e-12 };
                assert!((r.u_tau - 1.0).abs() < tol, "{t} {yp}: {}", r.u_tau);
                assert!((r.y_star - yp).abs() < 1e-9 * yp);
                let e = 1.0 / (0.41 * y);
                assert!((r.epsilon - e).abs() < tol * e, "{t}: {} {e}", r.epsilon);
                assert!((r.pk_wall - e).abs() < 2.0 * tol * e, "{t}: {}", r.pk_wall);
                assert!((r.prod_weight - 1.0).abs() < tol);
            }
        }
    }

    #[test]
    fn keps_wall_treatments_in_the_viscous_sublayer() {
        let c = crate::turb::keps::KepsConstants::standard();
        let nu = 1e-5;
        // y⁺ = 0.5 with u_τ = 1, k⁺ = 0.1 y⁺²
        let y = 0.5 * nu;
        let k: f64 = 0.1 * 0.25;
        let i = KepsWallInputs {
            ut: 0.5,
            y,
            nu,
            k,
            ks: 0.0,
            cs: 0.5,
            cu: 0.0,
        };
        for t in [STANDARD, COMPOUND] {
            let r = keps_wall_face(&i, &c, t);
            // the viscous wall shear exactly
            assert!((r.nu_wall - nu).abs() < 1e-6 * nu, "{t}: {}", r.nu_wall);
            assert!((r.u_tau - 1.0).abs() < 1e-6);
        }
        let r = keps_wall_face(&i, &c, COMPOUND);
        assert!((r.epsilon - 2.0 * nu * k / (y * y)).abs() < 1e-6 * r.epsilon);
        assert!(r.prod_weight.abs() < 1e-6 && r.pk_wall == 0.0);
        // scalable takes the log law at ỹ* = max(y*, 11.06) (it moves the
        // wall; accurate for y⁺ ≥ 30 only)
        let r = keps_wall_face(&i, &c, SCALABLE);
        let m = foehncraft_core::math::sqrt;
        let ys = m(m(c.c_mu)) * m(k) * y / nu;
        let up = foehncraft_core::math::ln(LOG_LAW_E * Y_PLUS_LOG) / 0.41;
        assert!((r.nu_wall - nu * ys / up).abs() < 1e-12 * nu);
        // no k: still the viscous shear (compound, standard) and finite
        let i0 = KepsWallInputs::<f64> { k: 0.0, ..i };
        for t in [SCALABLE, STANDARD, COMPOUND] {
            let r = keps_wall_face(&i0, &c, t);
            assert!(r.nu_wall.is_finite() && r.epsilon.is_finite() && r.pk_wall.is_finite());
        }
        assert_eq!(keps_wall_face(&i0, &c, COMPOUND).nu_wall, nu);
    }

    #[test]
    fn cebeci_bradshaw_shift_is_continuous() {
        let (cs, kappa) = (0.5, 0.41);
        assert_eq!(roughness_shift(1.0_f64, cs, kappa), 0.0);
        let at = |x: f64| roughness_shift(x, cs, kappa);
        // (to the four digits of the constants 0.811 = ln 2.25)
        assert!(at(2.25 + 1e-9).abs() < 1e-4);
        assert!((at(90.0 - 1e-9) - at(90.0)).abs() < 1e-6);
        // fully rough: u⁺ = (1/κ) ln(E y⁺/(1 + C_s k_s⁺))
        let v = at(1000.0);
        assert!((v - foehncraft_core::math::ln(501.0) / 0.41).abs() < 1e-12);
    }

    #[test]
    fn keps_wall_face_derivatives_by_dual_numbers() {
        use foehncraft_core::Dual;
        let c = crate::turb::keps::KepsConstants::standard();
        let base = KepsWallInputs {
            ut: 5.0,
            y: 2e-4,
            nu: 1.5e-5,
            k: 0.08,
            ks: 1e-4,
            cs: 0.5,
            cu: 0.0,
        };
        let lift = |x: f64| Dual::<1>::constant(x);
        for t in [SCALABLE, STANDARD, COMPOUND, COMPOUND_2023] {
            let di = KepsWallInputs {
                ut: lift(base.ut),
                y: lift(base.y),
                nu: lift(base.nu),
                k: Dual::<1>::variable(base.k, 0),
                ks: lift(base.ks),
                cs: lift(base.cs),
                cu: lift(0.0),
            };
            let r = keps_wall_face(&di, &c, t);
            let h = 1e-7;
            let at = |k: f64| keps_wall_face(&KepsWallInputs { k, ..base }, &c, t);
            let (p, m) = (at(base.k + h), at(base.k - h));
            for (dv, fd) in [
                (r.epsilon.d[0], (p.epsilon - m.epsilon) / (2.0 * h)),
                (r.nu_wall.d[0], (p.nu_wall - m.nu_wall) / (2.0 * h)),
                (r.pk_wall.d[0], (p.pk_wall - m.pk_wall) / (2.0 * h)),
            ] {
                assert!(
                    (dv - fd).abs() <= 1e-5 * fd.abs().max(1e-10),
                    "{t}: {dv} {fd}"
                );
            }
        }
    }

    #[test]
    fn compound_2023_tends_to_its_limits() {
        let c = crate::turb::keps::KepsConstants::standard();
        let nu = 1e-5;
        let m = foehncraft_core::math::sqrt;
        // log layer, k in equilibrium with u_τ = 1 (u* = 1): u⁺ = ln(8.9 y⁺)/κ
        for &(yp, tol, tg) in &[
            (30.0, 2.5e-2, 0.3),
            (100.0, 1e-3, 3e-2),
            (300.0, 1e-4, 5e-3),
        ] {
            let y = yp * nu;
            let u = foehncraft_core::math::ln(E_2023 * yp) / 0.41;
            let i = KepsWallInputs {
                ut: u,
                y,
                nu,
                k: 1.0 / m(c.c_mu),
                ks: 0.0,
                cs: 0.5,
                cu: 0.0,
            };
            let r = keps_wall_face(&i, &c, COMPOUND_2023);
            assert!((r.u_tau - 1.0).abs() < tol, "{yp}: {}", r.u_tau);
            let e = 1.0 / (0.41 * y);
            // Γ = dU⁺/dy⁺ → 1/(κ y⁺) decays slowly: ε keeps the weight Γ
            // of 2νk/y² (eq. (20)), about 9 % at y⁺ = 30
            let g = 1.0 / (0.41 * yp);
            assert!((r.viscous - g).abs() < tg * g, "{yp}: {} {g}", r.viscous);
            let ew = 2.0 * nu * i.k / (y * y);
            let eb = r.viscous * ew + (1.0 - r.viscous) * e;
            assert!(
                (r.epsilon - eb).abs() < 1e-12 * eb,
                "{yp}: {} {eb}",
                r.epsilon
            );
            // P = u*² (τ_w/ρν) Γ
            let pk = r.u_tau * r.u_tau / nu * r.viscous;
            assert!((r.pk_wall - pk).abs() < 1e-12 * pk, "{yp}: {}", r.pk_wall);
            assert!(
                (r.pk_wall - e).abs() < (tg + 3.0 * tol) * e,
                "{yp}: {}",
                r.pk_wall
            );
            assert_eq!(r.prod_weight, 1.0);
            assert_eq!(r.psi, 1.0);
        }
        // viscous sublayer: ν_w = ν, ε = 2νk/y², Γ = 1 (to 0.1 %)
        for yp in [1e-3, 0.1, 0.5, 1.0] {
            let y = yp * nu;
            let k: f64 = 0.1 * yp * yp;
            let i = KepsWallInputs {
                ut: yp,
                y,
                nu,
                k,
                ks: 0.0,
                cs: 0.5,
                cu: 0.0,
            };
            let r = keps_wall_face(&i, &c, COMPOUND_2023);
            assert!((r.nu_wall - nu).abs() < 1e-3 * nu, "{yp}: {}", r.nu_wall);
            assert!(r.viscous > 0.99, "{yp}: {}", r.viscous);
            let ew = 2.0 * nu * k / (y * y);
            assert!(
                (r.epsilon - ew).abs() < 1e-2 * ew,
                "{yp}: {} {ew}",
                r.epsilon
            );
        }
        // no k: the viscous shear, finite values
        let i0 = KepsWallInputs::<f64> {
            ut: 1.0,
            y: 1e-4,
            nu,
            k: 0.0,
            ks: 0.0,
            cs: 0.5,
            cu: 3.0,
        };
        let r = keps_wall_face(&i0, &c, COMPOUND_2023);
        assert_eq!(r.nu_wall, nu);
        assert!(r.epsilon.is_finite() && r.pk_wall.is_finite() && r.psi == 1.0);
    }

    #[test]
    fn compound_2023_gamma_is_the_derivative_of_the_blend() {
        // U⁺(y*) = y*/(y*/U⁺) from the kernel; Γ against its derivative
        let (nu, kappa) = (1.0, 0.41);
        let u_plus = |ys: f64| ys / compound_2023(ys, ys, nu, 1.0, 10.0, 0.0, kappa).0;
        for ys in [1.5, 3.0, 8.0, 15.0, 40.0, 200.0] {
            let h = 1e-6 * ys;
            let d = (u_plus(ys + h) - u_plus(ys - h)) / (2.0 * h);
            let g = compound_2023(ys, ys, nu, 1.0, 10.0, 0.0, kappa).1;
            assert!((g - d).abs() < 1e-6 * d.max(1e-3), "{ys}: {g} {d}");
        }
    }

    #[test]
    fn compound_2023_non_equilibrium_factor() {
        let (nu, kappa, us, ut) = (1e-5, 0.41, 1.0, 15.0);
        let at = |yp: f64, cu: f64| compound_2023(yp, yp * nu, nu, us, ut, cu, kappa);
        // 1 below y_v/2, continuous there, < 1 in an adverse gradient
        assert_eq!(at(5.0, 1e3).2, 1.0);
        assert!((at(5.5 + 1e-9, 1e3).2 - 1.0).abs() < 1e-9);
        let (yp, cu) = (50.0, 1e3);
        let psi = 1.0 - cu * (yp - 5.5) * nu / (ut * kappa * us);
        assert!((at(yp, cu).2 - psi).abs() < 1e-12, "{} {psi}", at(yp, cu).2);
        assert!(
            at(yp, cu).0 < at(yp, 0.0).0,
            "smaller wall shear in an adverse gradient"
        );
        // bounded
        assert_eq!(at(yp, 1e9).2, PSI_MIN);
        assert_eq!(at(yp, -1e9).2, PSI_MAX);
    }
}
