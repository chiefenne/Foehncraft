//! The k-ε-ζ-f model (D132; stage 6 specification 3.3, 6.2) per cell:
//! the time and length scales, the eddy viscosity, the point-implicit
//! source split of the k, ε and ζ equations and the right-hand side of the
//! elliptic f equation.
//!
//! K. Hanjalić, M. Popovac, M. Hadžiabdić, *A robust near-wall
//! elliptic-relaxation eddy-viscosity turbulence model for CFD*, Int. J.
//! Heat Fluid Flow 25, 1047–1051 (2004). Per unit mass, constant density:
//!
//!   ν_t = C_μ ζ k T
//!
//!   Dk/Dt = P_k − ε + ∇·[(ν + ν_t/σ_k)∇k]
//!
//!   Dε/Dt = (C_ε1 P_k − C_ε2 ε)/T + ∇·[(ν + ν_t/σ_ε)∇ε]
//!
//!   Dζ/Dt = f − (ζ/k) P_k + ∇·[(ν + ν_t/σ_ζ)∇ζ]
//!
//!   L²∇²f − f = (1/T)(C₁ − 1 + C₂′ P_k/ε)(ζ − 2/3)
//!
//! T = max[min(k/ε, a/(√6 C_μ |S| ζ)), C_T (ν/ε)^½],
//! L = C_L max[min(k^{3/2}/ε, k^½/(√6 C_μ |S| ζ)), C_η ν^¾/ε^¼],
//! P_k = ν_t S², S = (2 S_ij S_ij)^½, |S| = S (6d amendment of the
//! specification 2.1, which read (S_ij S_ij)^½: the authors'
//! implementation and the TMBWG page use S), C_ε1 = 1.4 (1 + 0.012/ζ).
//!
//! Constants: C_μ = 0.22, C_ε1 = 1.4 (1 + 0.012/ζ), C_ε2 = 1.9,
//! C₁ = 1.4 (C₁ − 1 = 0.4), C₂′ = 0.65, σ_k = 1, σ_ε = 1.3, σ_ζ = 1.2,
//! a = 0.6, C_T = 6, C_L = 0.36, C_η = 85. Checked on 2026-11-17 (the ◇
//! of the specification) against the model page "k-e-zeta-f2004" of the
//! Turbulence Modeling Benchmarking Working Group
//! (tmbwg.github.io/turbmodels), which reproduces the paper's equations and
//! table and cites the 2004 paper; it notes that the 2007 paper on the
//! compound wall treatment has a typo in T, so T is taken from the 2004
//! form above.
//!
//! The wall values of the integration to the wall are k_w = 0, ζ_w = 0,
//! ε_w = 2ν k_P/y_P² and f_w = −2ν ζ_P/y_P² (first-cell values), see
//! [`super::wall`]; the compound treatment blends ε_P, P_k,P and f_w.
//!
//! The source split (Patankar 1980; specification 6.2) keeps S_u ≥ 0 and
//! S_p ≤ 0 in the k, ε and ζ equations, with select operations only; the f
//! equation is linear (S_u = its right-hand side of either sign, S_p = −1).
//!
//! C_ε1 P_k is evaluated as 1.4 (P_k + 0.012 P_k/ζ) with P_k/ζ = C_μ k T S²
//! (ν_t/ζ), which stays finite as ζ → 0 at the wall, where the
//! coefficient 1 + 0.012/ζ alone does not.
//!
//! **The ζ-f₀ form** (`form` = [`FORM_F0`]; M. Popovac, *Robust eddy
//! viscosity turbulence modeling with elliptic relaxation and compound
//! wall treatment*, Energies 16, 3685 (2023), CC BY) solves
//! f₀ = f + (N_ζ − 1) ζ/T instead of f (the idea of Lien, Kalitzin &
//! Durbin 1998), which is zero at the wall with N_ζ = 2 (eq. (21)):
//!
//!   ν_t = (c_μ/ζ_eq) ζ k T, c_μ = 0.09, ζ_eq = 0.41 (eq. (4))
//!
//!   Dζ/Dt = f₀ − ζ [P_k/ε + (N_ζ − 1)]/T + ∇·[(ν + ν_t/σ_ζ)∇ζ] (eq. (8))
//!
//!   L²∇²f₀ − f₀ = (1/T)[(C_f1 + C_f2 P_k/ε)(ζ − 2/3) − (N_ζ − 1) ζ]
//!
//! with C_ε1 = 1.4 (1 + 0.012/√ζ) (our default: 1 + 0.009/ζ, below) and
//! the other constants of the 2004
//! model (the paper's table 1: C_f1 = C₁ − 1 = 0.4, C_f2 = C₂′ = 0.65,
//! C_T = 0.6 = a, C_τ = 6 = C_T here, C_L = 0.36, C_η = 85). The sign of
//! the (N_ζ − 1) ζ term is the one that follows from substituting f₀
//! into the 2004 f equation (the gradient of ζ/T neglected); eq. (7)
//! prints + (N_ζ − 1) ζ, which would count the term twice in the ζ
//! equation and pull the log-layer ζ from ≈ 0.4 to ≈ 0.17. The
//! realizability bounds of T and L take the coefficient of ν_t,
//! c_μ/ζ_eq = 0.2195 (the 0.22 of 2004), for which they are derived;
//! eq. (9) writes c_μ.
//!
//! Choices of ours (6d2, #141) make the form the default `k-epsilon-zeta-f0`:
//!
//! - The (N_ζ − 1) terms take ε/k instead of 1/T. Off the wall the two
//!   agree (T = k/ε); at the wall ε/k = 2ν/y², which is what makes
//!   f₀_w = (N_ζ − 2) 2νζ/y² = 0 (eq. (21)) consistent, whereas the
//!   Kolmogorov-bounded 1/T vanishes there.
//! - The sink of k ([`k_sink_rate`], both ζ-f models, #141) is
//!   max(k/T_K, 2νk/d²), T_K = max(k/ε, C_T (ν/ε)^½), d the wall
//!   distance, instead of ε (eq. (5)). With −ε the laminar part of a
//!   boundary layer has no positive steady k, the positivity floor holds
//!   it and the k residual stalls (the 6d finding on the flat plate).
//!   k/T_K alone (6d2 to #140) removes that, but its rate stays finite at
//!   the wall, so k falls like y instead of y² and ε_w = 2νk_P/y_P² grows
//!   with every refinement (the channel's U_b⁺ drifted from 17.79 to 18.39
//!   over four meshes). 2νk/d² is ε at the wall and negligible away from
//!   it, and it has no feedback from ε (blends with ε drove a limit cycle
//!   in the wall jet of the impinging jet). The realizability bound of T
//!   is left out of T_K: it would raise the sink in shear layers and
//!   relaminarise them (U_b⁺ +30 % in the channel). `k_sink` keeps −ε
//!   ([`K_SINK_EPSILON`]) and k/T_K ([`K_SINK_KOLMOGOROV`]) as options.
//! - C_ε1 = 1.4 (1 + 0.009/ζ) instead of the printed 1.4 (1 + 0.012/√ζ)
//!   (`c_e1_power` 1 instead of ½): with a sink that is ε at the wall the
//!   printed form gives U_b⁺ −11 % in the channel, and a larger √ζ
//!   coefficient spoils the step (√ζ is not small in its shear layer);
//!   with 0.009/ζ the channel is mesh-converged at +0.2 % from the DNS.
//! - The realizability bound of T takes a = 0.6 √2 = 0.849 instead
//!   of 0.6, which with S = (2 S_ij S_ij)^½ is the bound of Durbin (1996)
//!   written with (S_ij S_ij)^½. In the separated shear layer of the
//!   backward-facing step the bound with 0.6 is active and halves ν_t:
//!   the reattachment length came out at 7.4 to 7.7 step heights against
//!   6.26 measured; with 0.849 it was 6.38 in 6d2, 6.44 since #141 (the
//!   bound no longer acts there:
//!   a = 1.2 gives the same). The channel moves by 0.1 %, C_f of the flat
//!   plate by +1.5 % (within 3 % of White's correlation from y⁺ 1 to 30).
//!
//! The elliptic-blending k-ε (D133) shares the fields, the solvers and
//! these constants ([`FORM_EB`], [`ZetaFConstants::billard_laurence_2012`]);
//! its cell body is [`super::eb::eb_cell`].

use super::invariants::{Tensor, s2_omega2};
use super::{smax, smin};
use foehncraft_core::Real;

/// Lower bound of ζ in the divisions of the source split (the floor of the
/// solved ζ is the same, specification 6.2).
pub const ZETA_MIN: f64 = 1e-10;
/// Upper bound of ζ = v²/k (v² ≤ 2k).
pub const ZETA_MAX: f64 = 2.0;

/// The 2004 ζ-f (Hanjalić, Popovac & Hadžiabdić).
pub const FORM_ZETA_F: u32 = 0;
/// The ζ-f₀ of Popovac (2023).
pub const FORM_F0: u32 = 1;
/// The elliptic-blending k-ε of Billard & Laurence (2012): φ in the slot of
/// ζ, α in the slot of f.
pub const FORM_EB: u32 = 2;

/// The sink of k max(k/T_K, 2νk/d²), T_K = max(k/ε, C_T (ν/ε)^½), d the
/// wall distance (the default of ζ-f and ζ-f₀; module documentation).
pub const K_SINK_WALL: u32 = 0;
/// The sink of k ε of the papers (Hanjalić et al. 2004, eq. (5) of
/// Popovac 2023).
pub const K_SINK_EPSILON: u32 = 1;
/// The sink of k k/T_K (the ζ-f₀ default of 6d2, not grid-convergent at
/// the wall).
pub const K_SINK_KOLMOGOROV: u32 = 2;

/// The sink of k divided by k (1/s, the S_p of the k equation is its
/// negative) at k, ε, the viscosity ν and the wall distance d.
#[inline(always)]
pub fn k_sink_rate<T: Real>(k: T, e: T, nu: T, d: T, c: &ZetaFConstants) -> T {
    let f = T::from_f64;
    if c.k_sink == K_SINK_EPSILON {
        return e / k;
    }
    let t_kol = f(c.c_t) * (nu / e).sqrt_of();
    let kol = T::ONE / (k / e).max_of(t_kol);
    if c.k_sink == K_SINK_KOLMOGOROV {
        kol
    } else {
        kol.max_of(f(2.0) * nu / (d * d).max_of(f(1e-300)))
    }
}

/// The constants of the model (the paper's by default; a case may change
/// them, TURB007).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZetaFConstants {
    /// C_μ = 0.22 (in ν_t = C_μ ζ k T).
    pub c_mu: f64,
    /// The base of C_ε1 = 1.4 (1 + 0.012/ζ): 1.4.
    pub c_e1: f64,
    /// The ζ term of C_ε1: 0.012.
    pub c_e1_zeta: f64,
    /// C_ε2 = 1.9.
    pub c_e2: f64,
    /// C₁ = 1.4 (the f equation takes C₁ − 1).
    pub c1: f64,
    /// C₂′ = 0.65.
    pub c2: f64,
    /// σ_k = 1.0.
    pub sigma_k: f64,
    /// σ_ε = 1.3.
    pub sigma_e: f64,
    /// σ_ζ = 1.2.
    pub sigma_z: f64,
    /// a = 0.6 of the realizability bound of T and L.
    pub a: f64,
    /// C_T = 6 (Kolmogorov time-scale bound).
    pub c_t: f64,
    /// C_L = 0.36.
    pub c_l: f64,
    /// C_η = 85 (Kolmogorov length-scale bound).
    pub c_eta: f64,
    /// The equilibrium C_μ = 0.09 of the wall treatment (u* = C_μ^¼ k^½)
    /// and of the inlet and start values (ν_t ≈ 0.09 k²/ε where ζ ≈ 0.41
    /// in the log layer).
    pub c_mu_eq: f64,
    /// κ = 0.41 (wall treatment).
    pub kappa: f64,
    /// N_ζ of the ζ-f₀ form (2; 1 in the 2004 model, where f₀ = f; 2 for
    /// EB, where the wall kernel then writes α_w = 0).
    pub n_zeta: f64,
    /// [`FORM_ZETA_F`], [`FORM_F0`] or [`FORM_EB`].
    pub form: u32,
    /// C_ε3 = 2.3 of the E term of the k equation (EB).
    pub c_e3: f64,
    /// C_ε4 = 0.4 of C*_ε2 (EB).
    pub c_e4: f64,
    /// Factor of the molecular diffusion of k, ε and ζ (φ): 1; ½ for EB.
    pub mol: f64,
    /// The wall value ε_w = e_wall ν k_P/y_P²: 2; 1 for EB.
    pub e_wall: f64,
    /// The sink of k: [`K_SINK_WALL`] (the default of ζ-f and ζ-f₀),
    /// [`K_SINK_EPSILON`] or [`K_SINK_KOLMOGOROV`]; not used by EB.
    pub k_sink: u32,
    /// The power of ζ in C_ε1 = C_ε1′ (1 + C_ε1″/ζ^p): 1 (2004 and our
    /// ζ-f₀ default) or ½ (Popovac 2023 as printed).
    pub c_e1_power: f64,
}

impl ZetaFConstants {
    /// The constants of Hanjalić, Popovac & Hadžiabdić (2004).
    pub const fn paper() -> Self {
        Self {
            c_mu: 0.22,
            c_e1: 1.4,
            c_e1_zeta: 0.012,
            c_e2: 1.9,
            c1: 1.4,
            c2: 0.65,
            sigma_k: 1.0,
            sigma_e: 1.3,
            sigma_z: 1.2,
            a: 0.6,
            c_t: 6.0,
            c_l: 0.36,
            c_eta: 85.0,
            c_mu_eq: 0.09,
            kappa: 0.41,
            n_zeta: 1.0,
            form: FORM_ZETA_F,
            c_e3: 0.0,
            c_e4: 0.0,
            mol: 1.0,
            e_wall: 2.0,
            k_sink: K_SINK_WALL,
            c_e1_power: 1.0,
        }
    }

    /// The constants of the elliptic-blending k-ε (Billard & Laurence 2012;
    /// module documentation). `a` = 0.6 is the coefficient of T_lim, which
    /// takes (S_ij S_ij)^½.
    pub const fn billard_laurence_2012() -> Self {
        Self {
            c_mu: 0.22,
            c_e1: 1.44,
            c_e1_zeta: 0.0,
            c_e2: 1.83,
            c1: 1.7,
            c2: 0.9,
            sigma_k: 1.0,
            sigma_e: 1.5,
            sigma_z: 1.0,
            a: 0.6,
            c_t: 4.0,
            c_l: 0.164,
            c_eta: 75.0,
            c_mu_eq: 0.09,
            kappa: 0.41,
            n_zeta: 2.0,
            form: FORM_EB,
            c_e3: 2.3,
            c_e4: 0.4,
            mol: 0.5,
            e_wall: 1.0,
            k_sink: K_SINK_EPSILON,
            c_e1_power: 1.0,
        }
    }

    /// The constants of the ζ-f₀ model of Popovac (2023), table 1
    /// (C_μ = c_μ/ζ_eq = 0.09/0.41), with the realizability coefficient
    /// a = 0.6 √2 of our default (module documentation).
    pub const fn popovac_2023() -> Self {
        let mut c = Self::paper();
        c.c_mu = 0.09 / 0.41;
        c.a = 0.6 * core::f64::consts::SQRT_2;
        c.n_zeta = 2.0;
        c.form = FORM_F0;
        c.c_e1_zeta = 0.009;
        c
    }

    /// The constants the wall kernel of the k-ε family takes for ζ-f
    /// (C_μ the equilibrium value, κ).
    pub const fn wall(&self) -> super::keps::KepsConstants {
        let mut c = super::keps::KepsConstants::standard();
        c.c_mu = self.c_mu_eq;
        c.kappa = self.kappa;
        c.e_wall = self.e_wall;
        c
    }
}

impl Default for ZetaFConstants {
    fn default() -> Self {
        Self::paper()
    }
}

/// Run-time options of the model.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZetaFOptions {
    /// 1: Newton linearisation of the destruction of ε
    /// (S_p = −2 C_ε2/T), 0: Picard (−C_ε2/T, the default).
    pub newton: u32,
    /// Relative width of the smooth switches (0: the published hard form).
    pub smooth: f64,
    /// Ambient k of the sustaining terms (Spalart & Rumsey 2007); 0: off.
    pub k_amb: f64,
    /// Ambient ε of the sustaining terms; 0: off.
    pub e_amb: f64,
    /// Largest ν_t/ν (`turbulence.max_viscosity_ratio`).
    pub max_ratio: f64,
}

impl Default for ZetaFOptions {
    fn default() -> Self {
        Self {
            newton: 0,
            smooth: 0.0,
            k_amb: 0.0,
            e_amb: 0.0,
            max_ratio: 1e5,
        }
    }
}

/// Inputs of [`zeta_f_cell`] for one cell.
#[derive(Clone, Copy, Debug)]
pub struct ZetaFCellInputs<T> {
    /// k (m²/s²).
    pub k: T,
    /// ε (m²/s³).
    pub e: T,
    /// ζ = v²/k.
    pub zeta: T,
    /// f (1/s), the solution of the elliptic equation.
    pub f: T,
    /// Velocity gradient g_ij = ∂u_i/∂x_j.
    pub grad_u: Tensor<T>,
    /// ν (m²/s).
    pub nu: T,
    /// Weight m_p of the wall production (0 off the near-wall cells): the
    /// cell takes (1 − m_p) P_k + `wall_pk`.
    pub wall_m: T,
    /// The weighted wall production (m²/s³).
    pub wall_pk: T,
    /// ∇k (EB: the cross diffusion of φ; zero otherwise).
    pub grad_k: [T; 3],
    /// ∇ζ = ∇φ (EB).
    pub grad_z: [T; 3],
    /// (∂_l ∂_j U_i)(∂_l ∂_j U_i) (1/(m² s²); EB: the E term).
    pub g2: T,
    /// ∇·(ν_t/σ_k ∇k) (m²/s³; EB: C*_ε2).
    pub dkt: T,
    /// The wall distance (m; 1e10 without walls).
    pub d: T,
}

/// Results of [`zeta_f_cell`] for one cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZetaFCell<T> {
    /// ν_t (m²/s).
    pub nu_t: T,
    /// The time scale T (s).
    pub t: T,
    /// L² (m²), the diffusivity of the f equation (1 for the
    /// elliptic-blending model, whose α equation divides by L²).
    pub l2: T,
    /// The production entering the equations (with the wall part) (m²/s³).
    pub pk: T,
    /// S_u, S_p of k.
    pub su_k: T,
    /// See `su_k`.
    pub sp_k: T,
    /// S_u, S_p of ε.
    pub su_e: T,
    /// See `su_e`.
    pub sp_e: T,
    /// S_u, S_p of ζ.
    pub su_z: T,
    /// See `su_z`.
    pub sp_z: T,
    /// The right-hand side of −∇·(L²∇f) + f = R (1/s); S_u of the α
    /// equation (EB).
    pub rhs_f: T,
    /// S_p of the elliptic equation: −1; −1/L² for EB.
    pub sp_f: T,
    /// ν_t/σ_k.
    pub dk: T,
    /// ν_t/σ_ε.
    pub de: T,
    /// ν_t/σ_ζ.
    pub dz: T,
    /// 1 where the cap `max_ratio` of ν_t/ν acted, else 0.
    pub capped: T,
}

/// T, L and ν_t/ζ = C_μ k T of one cell (shared by [`zeta_f_cell`] and the
/// output of the time scale).
#[inline(always)]
pub fn scales<T: Real>(k: T, e: T, zeta: T, s_abs: T, nu: T, c: &ZetaFConstants, w: T) -> (T, T) {
    let f = T::from_f64;
    let tiny = f(1e-300);
    let den = f(foehncraft_core::math::sqrt(6.0) * c.c_mu) * s_abs * zeta;
    let sk = k.sqrt_of();
    let t_k = k / e;
    let l_k = k * sk / e;
    let (t_r, l_r) = if den > tiny {
        (smin(t_k, f(c.a) / den, w), smin(l_k, sk / den, w))
    } else {
        (t_k, l_k)
    };
    let t = smax(t_r, f(c.c_t) * (nu / e).sqrt_of(), w);
    let l = f(c.c_l) * smax(l_r, f(c.c_eta) * (nu * nu * nu / e).sqrt_of().sqrt_of(), w);
    (t, l)
}

/// ν_t, the scales and the source split of the k, ε and ζ equations and
/// the right-hand side of the f equation in one cell, with the values
/// k*, ε*, ζ*, f of the previous outer iteration (f of the current one:
/// the f equation is solved before ζ).
#[inline(always)]
pub fn zeta_f_cell<T: Real>(
    i: &ZetaFCellInputs<T>,
    c: &ZetaFConstants,
    o: &ZetaFOptions,
) -> ZetaFCell<T> {
    if c.form == FORM_EB {
        return super::eb::eb_cell(i, c, o);
    }
    let f = T::from_f64;
    let w = f(o.smooth);
    let tiny = f(1e-300);
    let (k, e, nu) = (i.k.max_of(T::ZERO), i.e.max_of(tiny), i.nu);
    let z = i.zeta.max_of(T::ZERO);
    let (s2, _) = s2_omega2(&i.grad_u);
    // |S| = S = (2 S_ij S_ij)^½ in the realizability bounds of T and L,
    // as in the authors' implementation (Popovac 2012) and the TMBWG
    // page; (S_ij S_ij)^½ = S/√2 would loosen both bounds by √2 (it gave
    // U_b⁺ +2.9 % instead of +1.9 % in the channel at Re_τ 395)
    let s_abs = s2.sqrt_of();
    let (t, l) = scales(k, e, z, s_abs, nu, c, w);
    // ν_t = a ζ, a = C_μ k T, capped at max_ratio ν
    let a = f(c.c_mu) * k * t;
    let raw = a * z;
    let cap = f(o.max_ratio) * nu;
    let over = raw > cap;
    let capped = if over { T::ONE } else { T::ZERO };
    let nu_t = raw.min_of(cap);
    // ν_t/ζ (a unless capped, where ζ > 0)
    let nt_z = if over { cap / z.max_of(tiny) } else { a };
    let m = i.wall_m;
    let keep = T::ONE - m;
    let pk = keep * nu_t * s2 + i.wall_pk;
    let f0 = c.form == FORM_F0;
    // P_k/ζ (2004) or P_k/√ζ (2023), finite as ζ → 0
    let p_over_z = if c.c_e1_power < 0.75 {
        keep * nt_z * z.sqrt_of() * s2 + i.wall_pk / z.max_of(f(ZETA_MIN)).sqrt_of()
    } else {
        keep * nt_z * s2 + i.wall_pk / z.max_of(f(ZETA_MIN))
    };
    let ce1p = f(c.c_e1) * (pk + f(c.c_e1_zeta) * p_over_z);
    let inv_t = T::ONE / t;
    let k_div = k.max_of(tiny);
    // k: S_u = P_k (+ the sink at the ambient state, which balances the
    // decay of the free stream), S_p = −(the sink)/k*
    let amb_k = if o.k_amb > 0.0 && o.e_amb > 0.0 {
        let (ka, ea) = (f(o.k_amb), f(o.e_amb));
        ka * k_sink_rate(ka, ea, nu, i.d, c)
    } else {
        T::ZERO
    };
    let su_k = pk + amb_k;
    let sp_k = -k_sink_rate(k_div, e, nu, i.d, c);
    // ε: (C_ε1 P_k − C_ε2 ε)/T (+ the ambient destruction)
    let amb_e = if o.k_amb > 0.0 && o.e_amb > 0.0 {
        let nu0 = i.nu.to_f64();
        let ta = (o.k_amb / o.e_amb).max(c.c_t * foehncraft_core::math::sqrt(nu0 / o.e_amb));
        f(c.c_e2 * o.e_amb / ta)
    } else {
        T::ZERO
    };
    let r = f(c.c_e2) * inv_t;
    let (su_e, sp_e) = if o.newton != 0 {
        (ce1p * inv_t + r * e + amb_e, -f(2.0) * r)
    } else {
        (ce1p * inv_t + amb_e, -r)
    };
    // ζ: f − (ζ/k) P_k (2004); f₀ − ζ [P_k/ε + (N_ζ − 1)]/T (2023)
    let n1 = f(c.n_zeta - 1.0);
    let su_z = i.f.max_of(T::ZERO);
    // the (N_ζ − 1) terms take ε/k (= 1/T off the wall): at the wall
    // ε/k = 2ν/y², which makes f₀_w = (N_ζ − 2) 2νζ/y² (eq. (21))
    // consistent; the Kolmogorov-bounded 1/T vanishes there and leaves ζ
    // without its wall sink (6d2 finding: U_b⁺ −12 % in the channel)
    let e_k = e / k_div;
    let sink = if f0 {
        pk / e * inv_t + n1 * e_k
    } else {
        pk / k_div
    };
    // f: −∇·(L²∇f) + f = (1/T)(C₁ − 1 + C₂′ P_k/ε)(2/3 − ζ)
    // [+ (N_ζ − 1) ζ/T for f₀]
    let g = inv_t * (f(c.c1 - 1.0) + f(c.c2) * pk / e);
    let mut rhs_f = g * (f(2.0 / 3.0) - z);
    if f0 {
        rhs_f += e_k * n1 * z;
    }
    // The ζ equation takes f from the previous outer iteration, and f
    // falls with ζ at the local rate b = −∂(rhs_f)/∂ζ = G − (N_ζ − 1) ε/k
    // (G = (C₁ − 1 + C₂′ P_k/ε)/T). Where P_k/ε is large (a jet started
    // into gas at rest), b Δt ≫ 1 and the lagged pair overshoots: ζ jumps
    // to its cap, f turns large and negative, ζ to its floor, and the run
    // diverges within a few steps. Linearising the source about ζ*,
    // f(ζ) ≈ f* − b (ζ − ζ*), puts b on the diagonal (b ζ* into S_u):
    // implicit in the stiff part and identical at convergence (ζ = ζ*).
    let b = (g - if f0 { e_k * n1 } else { T::ZERO }).max_of(T::ZERO);
    let su_z = su_z + b * z;
    let sp_z = i.f.min_of(T::ZERO) / z.max_of(f(ZETA_MIN)) - sink - b;
    ZetaFCell {
        nu_t,
        t,
        l2: l * l,
        pk,
        su_k,
        sp_k,
        su_e,
        sp_e,
        su_z,
        sp_z,
        rhs_f,
        sp_f: -T::ONE,
        dk: nu_t / f(c.sigma_k),
        de: nu_t / f(c.sigma_e),
        dz: nu_t / f(c.sigma_z),
        capped,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use foehncraft_core::Dual;

    fn shear(k: f64, e: f64, zeta: f64, f: f64, dudy: f64) -> ZetaFCellInputs<f64> {
        let mut g = [[0.0; 3]; 3];
        g[0][1] = dudy;
        ZetaFCellInputs {
            k,
            e,
            zeta,
            f,
            grad_u: g,
            nu: 1.5e-5,
            wall_m: 0.0,
            wall_pk: 0.0,
            grad_k: [0.0; 3],
            grad_z: [0.0; 3],
            g2: 0.0,
            dkt: 0.0,
            d: 1e10,
        }
    }

    #[test]
    fn scales_and_eddy_viscosity_in_simple_shear() {
        let c = ZetaFConstants::paper();
        let (k, e, z, s) = (0.01, 0.1, 0.4, 50.0);
        let r = zeta_f_cell(&shear(k, e, z, 0.0, s), &c, &ZetaFOptions::default());
        // |S| = S = |dU/dy|
        let s_abs = s;
        let t = (k / e)
            .min(0.6 / (6f64.sqrt() * 0.22 * s_abs * z))
            .max(6.0 * (1.5e-5 / e).sqrt());
        assert!((r.t - t).abs() < 1e-15 * t, "{} {t}", r.t);
        let l = 0.36
            * (k.powf(1.5) / e)
                .min(k.sqrt() / (6f64.sqrt() * 0.22 * s_abs * z))
                .max(85.0 * (1.5e-5f64.powi(3) / e).powf(0.25));
        assert!((r.l2 - l * l).abs() < 1e-12 * l * l);
        let nu_t = 0.22 * z * k * t;
        assert!((r.nu_t - nu_t).abs() < 1e-15 * nu_t);
        assert!((r.pk - nu_t * s * s).abs() < 1e-12 * r.pk);
        // C_ε1 P_k = 1.4 (1 + 0.012/ζ) P_k, Picard
        let ce1 = 1.4 * (1.0 + 0.012 / z);
        assert!((r.su_e - ce1 * r.pk / t).abs() < 1e-12 * r.su_e);
        assert!((r.sp_e + 1.9 / t).abs() < 1e-12 / t);
        // f: R = (1/T)(0.4 + 0.65 P/ε)(2/3 − ζ)
        let rf = (0.4 + 0.65 * r.pk / e) * (2.0 / 3.0 - z) / t;
        assert!((r.rhs_f - rf).abs() < 1e-12 * rf.abs());
        assert!((r.dz - nu_t / 1.2).abs() < 1e-15 * nu_t);
    }

    #[test]
    fn kolmogorov_bounds_act_at_small_k() {
        // k → 0 at fixed ε: T = C_T (ν/ε)^½, L = C_L C_η ν^¾/ε^¼
        let c = ZetaFConstants::paper();
        let r = zeta_f_cell(
            &shear(1e-12, 1.0, 0.01, 0.0, 0.0),
            &c,
            &ZetaFOptions::default(),
        );
        let nu: f64 = 1.5e-5;
        assert!((r.t - 6.0 * nu.sqrt()).abs() < 1e-15);
        let l = 0.36 * 85.0 * nu.powf(0.75);
        assert!((r.l2.sqrt() - l).abs() < 1e-12 * l);
    }

    #[test]
    fn production_term_of_epsilon_stays_finite_as_zeta_vanishes() {
        let c = ZetaFConstants::paper();
        let o = ZetaFOptions::default();
        let r0 = zeta_f_cell(&shear(1e-3, 1e-2, 0.0, 0.0, 100.0), &c, &o);
        assert_eq!(r0.nu_t, 0.0);
        assert!(r0.su_e.is_finite() && r0.su_e > 0.0);
        // the limit of 1.4 (P + 0.012 P/ζ) for ζ → 0 is 1.4·0.012 C_μ k T S²
        let r1 = zeta_f_cell(&shear(1e-3, 1e-2, 1e-9, 0.0, 100.0), &c, &o);
        assert!((r1.su_e - r0.su_e).abs() < 1e-6 * r0.su_e);
        assert!(r0.sp_z.is_finite() && r0.su_z == 0.0);
    }

    #[test]
    fn sources_keep_their_signs() {
        let c = ZetaFConstants::paper();
        for newton in [0, 1] {
            for f in [-50.0, 0.0, 30.0] {
                for z in [0.0, 1e-3, 0.4, 1.9] {
                    let o = ZetaFOptions {
                        newton,
                        ..Default::default()
                    };
                    let mut i = shear(1e-3, 5e-2, z, f, 300.0);
                    i.wall_m = 0.5;
                    i.wall_pk = 1.0;
                    let r = zeta_f_cell(&i, &c, &o);
                    assert!(r.su_k >= 0.0 && r.sp_k <= 0.0, "{r:?}");
                    assert!(r.su_e >= 0.0 && r.sp_e <= 0.0, "{r:?}");
                    assert!(r.su_z >= 0.0 && r.sp_z <= 0.0, "{r:?}");
                    assert!(r.rhs_f.is_finite() && r.l2 > 0.0);
                }
            }
        }
    }

    #[test]
    fn homogeneous_decay_relaxes_zeta_to_two_thirds() {
        // no strain: f = R = (0.4/T)(2/3 − ζ) drives ζ towards 2/3
        let c = ZetaFConstants::paper();
        let o = ZetaFOptions::default();
        for z in [0.2, 2.0 / 3.0, 1.2] {
            let r = zeta_f_cell(&shear(1.0, 1.0, z, 0.0, 0.0), &c, &o);
            assert!((r.rhs_f - 0.4 * (2.0 / 3.0 - z)).abs() < 1e-15);
            assert_eq!(r.pk, 0.0);
            // k, ε: dk/dt = −ε, dε/dt = −C_ε2 ε²/k with T = k/ε
            assert!((r.sp_e + 1.9).abs() < 1e-15 && r.su_e == 0.0);
        }
    }

    #[test]
    fn sustaining_terms_balance_the_ambient_decay() {
        let c = ZetaFConstants::paper();
        let (ka, ea) = (1e-4, 2e-3);
        let o = ZetaFOptions {
            k_amb: ka,
            e_amb: ea,
            ..Default::default()
        };
        let r = zeta_f_cell(&shear(ka, ea, 2.0 / 3.0, 0.0, 0.0), &c, &o);
        assert!((r.su_k + r.sp_k * ka).abs() < 1e-18);
        assert!((r.su_e + r.sp_e * ea).abs() < 1e-15);
    }

    #[test]
    fn dual_numbers_give_the_derivative_of_the_sources() {
        let mut i = shear(0.02, 0.3, 0.35, -4.0, 40.0);
        i.grad_u[1][0] = 5.0;
        i.grad_u[0][0] = 3.0;
        i.grad_u[1][1] = -2.0;
        i.grad_u[2][2] = -1.0;
        let c = ZetaFConstants::paper();
        let o = ZetaFOptions {
            smooth: 1e-3,
            ..Default::default()
        };
        let lift = |x: f64| Dual::<1>::constant(x);
        for which in 0..3 {
            let pick = |j: usize, x: f64| {
                if j == which {
                    Dual::<1>::variable(x, 0)
                } else {
                    lift(x)
                }
            };
            let di = ZetaFCellInputs {
                k: pick(0, i.k),
                e: pick(1, i.e),
                zeta: pick(2, i.zeta),
                f: lift(i.f),
                grad_u: i.grad_u.map(|r| r.map(lift)),
                nu: lift(i.nu),
                wall_m: lift(0.0),
                wall_pk: lift(0.0),
                grad_k: [lift(0.0); 3],
                grad_z: [lift(0.0); 3],
                g2: lift(0.0),
                dkt: lift(0.0),
                d: lift(i.d),
            };
            let r = zeta_f_cell(&di, &c, &o);
            let x0 = [i.k, i.e, i.zeta][which];
            let h = 1e-5 * x0;
            let at = |x: f64| {
                let mut j = i;
                match which {
                    0 => j.k = x,
                    1 => j.e = x,
                    _ => j.zeta = x,
                }
                zeta_f_cell(&j, &c, &o)
            };
            let (p, m) = (at(x0 + h), at(x0 - h));
            for (dv, fd) in [
                (r.su_e.d[0], (p.su_e - m.su_e) / (2.0 * h)),
                (r.nu_t.d[0], (p.nu_t - m.nu_t) / (2.0 * h)),
                (r.rhs_f.d[0], (p.rhs_f - m.rhs_f) / (2.0 * h)),
                (r.sp_z.d[0], (p.sp_z - m.sp_z) / (2.0 * h)),
                (r.l2.d[0], (p.l2 - m.l2) / (2.0 * h)),
            ] {
                assert!(
                    (dv - fd).abs() <= 1e-5 * fd.abs().max(1e-8),
                    "{which}: {dv} {fd}"
                );
            }
        }
    }

    #[test]
    fn f0_form_matches_the_2004_zeta_source() {
        // with f = f₀ − (N_ζ − 1) ζ/T and T = k/ε the ζ equation of the
        // f₀ form is the 2004 one; C_ε1 takes √ζ; ν_t takes 0.09/0.41
        let c4 = ZetaFConstants::paper();
        let c0 = ZetaFConstants::popovac_2023();
        let o = ZetaFOptions::default();
        let (k, e, z, s) = (0.01, 0.1, 0.4, 3.0);
        let base = zeta_f_cell(&shear(k, e, z, 0.0, s), &c0, &o);
        assert!((base.t - k / e).abs() < 1e-15, "T = k/ε here: {}", base.t);
        let t = base.t;
        for f0 in [-3.0, 0.0, 0.7, 5.0] {
            let r0 = zeta_f_cell(&shear(k, e, z, f0, s), &c0, &o);
            let f = f0 - z / t;
            let mut c4b = c4;
            c4b.c_mu = c0.c_mu;
            let r4 = zeta_f_cell(&shear(k, e, z, f, s), &c4b, &o);
            let net = |r: &ZetaFCell<f64>| r.su_z + r.sp_z * z;
            assert!(
                (net(&r0) - net(&r4)).abs() < 1e-12,
                "{f0}: {} {}",
                net(&r0),
                net(&r4)
            );
            // f₀ = f + ζ/T in the right-hand side too
            assert!((r0.rhs_f - (r4.rhs_f + z / t)).abs() < 1e-12);
        }
        let nu_t = 0.09 / 0.41 * z * k * t;
        assert!((base.nu_t - nu_t).abs() < 1e-15 * nu_t);
        // the sink of k: k/T_K, ε where T_K = k/ε
        assert!((base.sp_k + e / k).abs() < 1e-12 * e / k);
        // and below ε where the Kolmogorov bound acts (k → 0)
        let (k1, e1) = (1e-8, 1.0);
        let r = zeta_f_cell(&shear(k1, e1, z, 0.0, s), &c0, &o);
        let tk = 6.0 * (1.5e-5f64 / e1).sqrt();
        assert!(
            (r.sp_k + 1.0 / tk).abs() < 1e-12 / tk,
            "{} {}",
            r.sp_k,
            -1.0 / tk
        );
        assert!(-r.sp_k * k1 < e1);
        let ce1 = 1.4 * (1.0 + 0.009 / z);
        assert!((base.su_e - ce1 * base.pk / t).abs() < 1e-12 * base.su_e);
        // the printed form C_ε1″/√ζ with `c_epsilon1_power` = ½
        let mut c5 = c0;
        c5.c_e1_power = 0.5;
        c5.c_e1_zeta = 0.012;
        let r5 = zeta_f_cell(&shear(k, e, z, 0.0, s), &c5, &o);
        let ce1 = 1.4 * (1.0 + 0.012 / z.sqrt());
        assert!((r5.su_e - ce1 * r5.pk / t).abs() < 1e-12 * r5.su_e);
        // homogeneous decay: the net ζ source is (0.4/T)(2/3 − ζ)
        for z in [0.2, 1.2] {
            let r = zeta_f_cell(&shear(1.0, 1.0, z, 0.0, 0.0), &c0, &o);
            let rf = r.rhs_f;
            let net = rf.max(0.0) + (rf.min(0.0) / z - 1.0) * z;
            assert!((net - 0.4 * (2.0 / 3.0 - z)).abs() < 1e-14, "{z}: {net}");
        }
        // finite at ζ = 0 with signs kept
        let r = zeta_f_cell(&shear(1e-3, 1e-2, 0.0, -2.0, 100.0), &c0, &o);
        assert!(r.su_e.is_finite() && r.sp_z.is_finite() && r.sp_z <= 0.0 && r.su_z == 0.0);
    }
    #[test]
    fn wall_sink_is_epsilon_at_the_wall_and_k_over_t_away_from_it() {
        let nu = 1.5e-5;
        let o = ZetaFOptions::default();
        for c in [ZetaFConstants::paper(), ZetaFConstants::popovac_2023()] {
            assert_eq!(c.k_sink, K_SINK_WALL);
            // the sublayer: k = ε d²/(2ν), the Kolmogorov bound acts; the
            // sink is 2νk/d² = ε exactly
            let (e, d) = (0.3, 1e-4);
            let k = e * d * d / (2.0 * nu);
            let mut i = shear(k, e, 0.01, 0.0, 10.0);
            i.d = d;
            let r = zeta_f_cell(&i, &c, &o);
            assert!((-r.sp_k * k - e).abs() < 1e-12 * e, "{}", -r.sp_k * k);
            // the same cell away from walls: k/T_K, far below ε
            i.d = 1e10;
            let r = zeta_f_cell(&i, &c, &o);
            let tk = 6.0 * (nu / e).sqrt();
            assert!((r.sp_k + 1.0 / tk).abs() < 1e-12 / tk);
            assert!(-r.sp_k * k < 0.01 * e);
            // the other sinks
            for (code, want) in [(K_SINK_EPSILON, e / k), (K_SINK_KOLMOGOROV, 1.0 / tk)] {
                let mut c2 = c;
                c2.k_sink = code;
                i.d = d;
                let r = zeta_f_cell(&i, &c2, &o);
                assert!((r.sp_k + want).abs() < 1e-12 * want, "{code}");
            }
            // the ambient source balances the sink of the ambient state
            let (ka, ea) = (1e-8, 1e-4);
            let oa = ZetaFOptions {
                k_amb: ka,
                e_amb: ea,
                ..Default::default()
            };
            let mut j = shear(ka, ea, 2.0 / 3.0, 0.0, 0.0);
            j.d = 0.01;
            let r = zeta_f_cell(&j, &c, &oa);
            assert!((r.su_k + r.sp_k * ka).abs() < 1e-12 * r.su_k);
        }
    }
}
