//! Kernel body of the droplet evaporation of the Lagrangian sprays (L12,
//! SP5, D251): the diameter, the temperature and the evaporated mass of one
//! droplet over an interval in which the gas state at the droplet and the
//! transfer coefficients are frozen.
//!
//! **Surface state.** The vapour at the surface is saturated at the droplet
//! temperature: mole fraction X_s = p_sat(T_p)/p₀, mass fraction
//! Y_s = X_s M_v / (X_s M_v + (1 − X_s) M_g), with M_g the molar mass of the
//! gas without the vapour. Spalding's mass transfer number is
//! B_M = (Y_s − Y_∞)/(1 − Y_s).
//!
//! **Rates** (a uniform droplet temperature, Biot number small). With the
//! Sherwood and Nusselt numbers of a sphere in a stream, Sh₀ = 2 + 0.6
//! Re^½ Sc^⅓ and Nu₀ = 2 + 0.6 Re^½ Pr^⅓ (Ranz and Marshall 1952):
//!
//! - `ranz-marshall` (the classical model, Spalding 1953; Faeth 1977):
//!   ṁ = π d ρ D Sh₀ ln(1 + B_M), Nu = Nu₀;
//! - `abramzon-sirignano` (Abramzon and Sirignano 1989): the film
//!   thickening by the blowing, Sh* = 2 + (Sh₀ − 2)/F(B_M),
//!   Nu* = 2 + (Nu₀ − 2)/F(B_T), F(B) = (1 + B)^0.7 ln(1 + B)/B, with
//!   B_T = (1 + B_M)^φ − 1, φ = (c_p,v/c_p,g)(Sh*/Nu*)/Le, found by a fixed
//!   number of substitutions (the same bits on every target).
//!
//! The heat to the droplet is Q = π d k Nu (T_∞ − T_p) ln(1 + B_T)/B_T (the
//! correction for the outflow of the vapour; it tends to 1 as B → 0), and
//! the droplet balance is m c_l dT_p/dt = Q − ṁ L.
//!
//! **Exact integration with frozen coefficients** (D251). ṁ and Q are both
//! proportional to d, so with Sh, Nu, B_M and the gas frozen the diameter
//! follows the d² law d² = d₀² − K t, K = 4 ρ D Sh ln(1 + B_M)/ρ_l (Godsave
//! 1953; Spalding 1953), the droplet's life ends at t = d₀²/K inside the
//! interval when K Δt ≥ d₀², and the temperature relaxes towards the
//! equilibrium T_eq = T_∞ − ṁ L/(Q/(T_∞ − T_p)) (the wet-bulb state of the
//! frozen coefficients) with the response time τ(t) = m c_l/(π d k Nu
//! ln(1+B_T)/B_T) ∝ d², which falls with d²: the decay over the interval is
//! exactly (d₁²/d₀²)^(d₀²/(K τ₀)) (exp(−Δt/τ₀) as K → 0). Mass is never
//! negative and the temperature never passes T_eq. A negative B_M (the gas
//! above saturation at the droplet temperature) gives condensation by the
//! same law.
//!
//! The body is a gather over one droplet with the in-house exp and ln of
//! `foehncraft_core::math`, so CPU and GPU agree bit for bit (D93 (6)).

use foehncraft_core::Real;

/// The evaporation model (`[spray.evaporation] model`).
pub mod model {
    /// Ranz-Marshall correlations, the classical film theory.
    pub const RANZ_MARSHALL: u32 = 1;
    /// Abramzon and Sirignano (1989): film thickening by the blowing.
    pub const ABRAMZON_SIRIGNANO: u32 = 2;
}

/// Water: the saturation pressure of the auxiliary equation of the IAPWS-95
/// formulation (Wagner and Pruss 1993, eq. 2.5 of Wagner and Pruss 2002).
pub mod water {
    use foehncraft_core::Real;

    /// Critical temperature (K).
    pub const T_CRITICAL: f64 = 647.096;
    /// Critical pressure (Pa).
    pub const P_CRITICAL: f64 = 22.064e6;
    /// Molar mass (kg/mol).
    pub const MOLAR_MASS: f64 = 0.018_015_268;
    /// The coefficients a₁ … a₆ of the vapour-pressure equation.
    const A: [f64; 6] = [
        -7.859_517_83,
        1.844_082_59,
        -11.786_649_7,
        22.680_741_1,
        -15.961_871_9,
        1.801_225_02,
    ];

    /// Saturation pressure of water (Pa) at `t` (K):
    /// ln(p/p_c) = (T_c/T)(a₁ϑ + a₂ϑ^1.5 + a₃ϑ³ + a₄ϑ^3.5 + a₅ϑ⁴ + a₆ϑ^7.5),
    /// ϑ = 1 − T/T_c. Valid from the triple point (273.16 K) to the
    /// critical point; below 273.16 K the same formula is used for the
    /// supercooled liquid, at and above T_c it returns p_c.
    #[inline]
    pub fn saturation_pressure(t: f64) -> f64 {
        if t >= T_CRITICAL {
            return P_CRITICAL;
        }
        let th = 1.0 - t / T_CRITICAL;
        let s = th.sqrt_of();
        let th3 = th * th * th;
        let sum = A[0] * th
            + A[1] * th * s
            + A[2] * th3
            + A[3] * th3 * s
            + A[4] * th3 * th
            + A[5] * th3 * th3 * th * s;
        P_CRITICAL * (T_CRITICAL / t * sum).exp_of()
    }
}

/// The state of one droplet and the gas for the evaporation step.
#[derive(Clone, Copy, Debug)]
pub struct EvapInput {
    /// Droplet temperature (K).
    pub tp: f64,
    /// Gas temperature at the droplet (K).
    pub tg: f64,
    /// Vapour mass fraction of the gas at the droplet.
    pub yg: f64,
    /// Thermodynamic pressure p₀ (Pa).
    pub p0: f64,
    /// Molar mass of the vapour (kg/mol).
    pub m_vapour: f64,
    /// Molar mass of the gas without the vapour (kg/mol).
    pub m_gas: f64,
    /// Saturation pressure of the liquid at `tp` (Pa).
    pub psat: f64,
    /// The vapour mass fraction at the surface from a phase equilibrium
    /// (a high-pressure model, D310 to D319); negative: from `psat` by
    /// Raoult's law ([`surface_fraction`]), the default.
    pub ys_surface: f64,
    /// Latent heat at `tp` (J/kg).
    pub latent: f64,
    /// Particle Reynolds number (of the drag step).
    pub re: f64,
    /// Gas Prandtl number μ c_p/k.
    pub pr: f64,
    /// Schmidt number of the vapour in the gas, μ/(ρ D).
    pub sc: f64,
    /// Gas density (kg/m³).
    pub rho_g: f64,
    /// Diffusivity of the vapour in the gas (m²/s).
    pub diff: f64,
    /// Gas thermal conductivity (W/(m K)).
    pub k_g: f64,
    /// Gas specific heat (J/(kg K)).
    pub cp_g: f64,
    /// Specific heat of the vapour (J/(kg K)).
    pub cp_v: f64,
    /// Droplet diameter (m).
    pub d: f64,
    /// Liquid density (kg/m³).
    pub rho_l: f64,
    /// Liquid specific heat (J/(kg K)).
    pub cp_l: f64,
    /// Interval (s).
    pub dt: f64,
    /// [`model`].
    pub model: u32,
}

/// The result of [`evaporation_step`].
#[derive(Clone, Copy, Debug, Default)]
pub struct EvapStep {
    /// Diameter at the end of the interval (m); 0 when the droplet has
    /// evaporated in it.
    pub d: f64,
    /// Temperature at the end (K).
    pub tp: f64,
    /// Time until the droplet has evaporated, from the start of the
    /// interval (s); infinite when it lasts.
    pub life: f64,
    /// Rate of the d² law K (m²/s): d² = d₀² − K t (negative: condensation).
    pub k: f64,
    /// Vapour mass fraction at the surface.
    pub ys: f64,
    /// Spalding mass transfer number B_M.
    pub bm: f64,
    /// Equilibrium temperature of the frozen coefficients (K).
    pub t_eq: f64,
    /// Thermal response time at the start (s).
    pub tau: f64,
    /// Decay of the temperature difference over the interval.
    pub decay: f64,
    /// Sherwood number used.
    pub sh: f64,
    /// Nusselt number used.
    pub nu: f64,
}

/// ln(1 + x) for x > −1, accurate near 0.
#[inline(always)]
fn ln1p(x: f64) -> f64 {
    if x.abs() < 1e-4 {
        x * (1.0 - x * (0.5 - x * (1.0 / 3.0 - 0.25 * x)))
    } else {
        (1.0 + x).ln_of()
    }
}

/// ln(1 + x)/x, 1 at x = 0.
#[inline(always)]
fn ln1p_over(x: f64) -> f64 {
    if x.abs() < 1e-4 {
        1.0 - x * (0.5 - x * (1.0 / 3.0 - 0.25 * x))
    } else {
        ln1p(x) / x
    }
}

/// The film-thickening function F(B) = (1 + B)^0.7 ln(1 + B)/B of
/// Abramzon and Sirignano (1989); 1 for B ≤ 0 (no blowing).
#[inline(always)]
pub fn film_factor(b: f64) -> f64 {
    if b <= 0.0 {
        return 1.0;
    }
    (0.7 * ln1p(b)).exp_of() * ln1p_over(b)
}

/// x^(1/3) for x > 0.
#[inline(always)]
fn cbrt(x: f64) -> f64 {
    (x.ln_of() / 3.0).exp_of()
}

/// Vapour mass fraction at the surface of a droplet with the saturation
/// pressure `psat` at the pressure `p0` (X_s held below 0.999: a droplet
/// at its boiling point is not resolved by this model).
#[inline]
pub fn surface_fraction(psat: f64, p0: f64, m_vapour: f64, m_gas: f64) -> f64 {
    let xs = (psat / p0).min_of(0.999).max_of(0.0);
    xs * m_vapour / (xs * m_vapour + (1.0 - xs) * m_gas)
}

/// Diameter, temperature and life of one droplet over `dt` with the
/// coefficients frozen at the start (module documentation).
#[inline]
pub fn evaporation_step(i: &EvapInput) -> EvapStep {
    let ys = if i.ys_surface >= 0.0 {
        i.ys_surface.min_of(0.999)
    } else {
        surface_fraction(i.psat, i.p0, i.m_vapour, i.m_gas)
    };
    // B_M > −1 for any Y_∞ < 1
    let bm = (ys - i.yg) / (1.0 - ys);
    let lm = ln1p(bm);
    let root = if i.re > 1e-30 { i.re.sqrt_of() } else { 0.0 };
    let sh0 = 2.0 + 0.6 * root * cbrt(i.sc.max_of(1e-30));
    let nu0 = 2.0 + 0.6 * root * cbrt(i.pr.max_of(1e-30));
    let le = i.k_g / (i.rho_g * i.cp_g * i.diff);
    let (sh, nu, bt) = if i.model == model::ABRAMZON_SIRIGNANO {
        let sh = 2.0 + (sh0 - 2.0) / film_factor(bm);
        let mut nu = nu0;
        let mut bt = 0.0;
        // φ and B_T by a fixed number of substitutions (converged to
        // round-off for the B of sprays; a fixed count keeps the bits)
        let mut k = 0;
        while k < 8 {
            let phi = i.cp_v / i.cp_g * sh / nu / le;
            bt = (phi * lm).exp_of() - 1.0;
            nu = 2.0 + (nu0 - 2.0) / film_factor(bt);
            k += 1;
        }
        (sh, nu, bt)
    } else {
        let phi = i.cp_v / i.cp_g * sh0 / nu0 / le;
        (sh0, nu0, (phi * lm).exp_of() - 1.0)
    };
    let pi = core::f64::consts::PI;
    // per unit diameter: ṁ/d and the heat conductance h/d (W/(m K))
    let mdot_d = pi * i.rho_g * i.diff * sh * lm;
    let cond_d = pi * i.k_g * nu * ln1p_over(bt);
    let k_rate = 4.0 * i.rho_g * i.diff * sh * lm / i.rho_l;
    let t_eq = i.tg - mdot_d * i.latent / cond_d;
    let d0 = i.d;
    let d02 = d0 * d0;
    // τ₀ = m c_l / (h) = ρ_l π d³/6 c_l / (cond_d d) = ρ_l c_l d²/(6/π cond_d)
    let tau = i.rho_l * i.cp_l * pi * d02 / (6.0 * cond_d);
    let life = if k_rate > 0.0 {
        d02 / k_rate
    } else {
        f64::INFINITY
    };
    if i.dt >= life {
        return EvapStep {
            d: 0.0,
            tp: i.tp,
            life,
            k: k_rate,
            ys,
            bm,
            t_eq,
            tau,
            decay: 0.0,
            sh,
            nu,
        };
    }
    let s = k_rate * i.dt / d02;
    let d12 = d02 * (1.0 - s);
    let decay = if s.abs() < 1e-8 {
        (-(i.dt / tau)).exp_of()
    } else {
        // exp(∫ −dt/τ(t)) with τ(t) = τ₀ (1 − K t/d₀²)
        (d02 / (k_rate * tau) * ln1p(-s)).exp_of()
    };
    EvapStep {
        d: d12.sqrt_of(),
        tp: t_eq + (i.tp - t_eq) * decay,
        life,
        k: k_rate,
        ys,
        bm,
        t_eq,
        tau,
        decay,
        sh,
        nu,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input(tp: f64, tg: f64, yg: f64, d: f64, dt: f64, model: u32) -> EvapInput {
        // dry air at 1 atm; vapour diffusivity 2.6e-5, Sc = 0.6
        EvapInput {
            tp,
            tg,
            yg,
            p0: 101_325.0,
            m_vapour: water::MOLAR_MASS,
            m_gas: 0.028_97,
            psat: water::saturation_pressure(tp),
            ys_surface: -1.0,
            latent: 2.44e6,
            re: 0.0,
            pr: 0.7,
            sc: 0.6,
            rho_g: 1.0,
            diff: 2.6e-5,
            k_g: 0.03,
            cp_g: 1010.0,
            cp_v: 1880.0,
            d,
            rho_l: 998.0,
            cp_l: 4181.0,
            dt,
            model,
        }
    }

    #[test]
    fn saturation_pressure_of_water_at_its_anchor_points() {
        // triple point 611.657 Pa, 25 °C 3169.9 Pa, normal boiling point
        // (373.124 K on ITS-90) 101.325 kPa, all from the IAPWS-95 tables
        for (t, p, tol) in [
            (273.16, 611.657, 2e-4),
            (298.15, 3169.9, 1e-3),
            (373.124, 101_325.0, 2e-4),
        ] {
            let ps = water::saturation_pressure(t);
            assert!((ps / p - 1.0).abs() < tol, "{t} K: {ps} Pa, expected {p}");
        }
    }

    #[test]
    fn the_d2_law_is_exact_with_frozen_coefficients_at_any_interval() {
        // at the equilibrium temperature K is constant: one interval and
        // 100 give the same diameter
        for m in [model::RANZ_MARSHALL, model::ABRAMZON_SIRIGNANO] {
            // the wet-bulb temperature of the frozen coefficients, by
            // bisection of T_eq(T) − T
            let (mut lo, mut hi) = (280.0f64, 399.0f64);
            for _ in 0..80 {
                let t = 0.5 * (lo + hi);
                let e = evaporation_step(&EvapInput {
                    psat: water::saturation_pressure(t),
                    ..input(t, 400.0, 0.0, 50e-6, 0.0, m)
                });
                if e.t_eq > t {
                    lo = t;
                } else {
                    hi = t;
                }
            }
            let teq = 0.5 * (lo + hi);
            let one = evaporation_step(&EvapInput {
                psat: water::saturation_pressure(teq),
                ..input(teq, 400.0, 0.0, 50e-6, 0.05, m)
            });
            let exact = (50e-6f64 * 50e-6 - one.k * 0.05).sqrt();
            assert!((one.d - exact).abs() < 1e-15, "{} {exact}", one.d);
            let mut d = 50e-6;
            let mut t = teq;
            for _ in 0..100 {
                let s = evaporation_step(&EvapInput {
                    psat: water::saturation_pressure(t),
                    ..input(t, 400.0, 0.0, d, 5e-4, m)
                });
                d = s.d;
                t = s.tp;
            }
            assert!((d / one.d - 1.0).abs() < 1e-6, "{d} {}", one.d);
            // the life ends inside an interval: diameter zero, mass not
            // negative
            let long = evaporation_step(&EvapInput {
                psat: water::saturation_pressure(teq),
                ..input(teq, 400.0, 0.0, 50e-6, 10.0, m)
            });
            assert_eq!(long.d, 0.0);
            assert!(long.life < 10.0);
        }
    }

    #[test]
    fn the_temperature_relaxes_to_the_wet_bulb_state_without_overshoot() {
        // Ranz-Marshall at Re = 0 (Nu = Sh = 2): the wet-bulb temperature
        // solves k (T_g − T) ln(1+B_T)/B_T = ρ D ln(1 + B_M) L
        let mut t = 293.15;
        let mut d = 100e-6;
        let mut last = t;
        for _ in 0..2000 {
            let s = evaporation_step(&EvapInput {
                psat: water::saturation_pressure(t),
                ..input(t, 373.15, 0.0, d, 1e-3, model::RANZ_MARSHALL)
            });
            assert!(s.tp >= last - 1e-12, "monotone");
            last = s.tp;
            t = s.tp;
            d = s.d;
            if d < 60e-6 {
                break;
            }
        }
        // the wet bulb of dry air at 373 K, 1 atm (own fixed point of the
        // same balance with the latent heat frozen): between 300 and 315 K
        assert!(t > 300.0 && t < 315.0, "wet bulb {t}");
        let s = evaporation_step(&EvapInput {
            psat: water::saturation_pressure(t),
            ..input(t, 373.15, 0.0, d, 0.0, model::RANZ_MARSHALL)
        });
        assert!((s.t_eq - t).abs() < 1e-3, "{} {t}", s.t_eq);
    }

    #[test]
    fn saturated_gas_neither_evaporates_nor_condenses_at_the_droplet_temperature() {
        let t = 310.0;
        let ys = surface_fraction(
            water::saturation_pressure(t),
            101_325.0,
            water::MOLAR_MASS,
            0.028_97,
        );
        let s = evaporation_step(&input(t, t, ys, 40e-6, 0.1, model::ABRAMZON_SIRIGNANO));
        assert!(s.k.abs() < 1e-20, "{}", s.k);
        assert!((s.d - 40e-6).abs() < 1e-18);
        // above saturation the droplet grows
        let c = evaporation_step(&input(t, t, ys * 1.1, 40e-6, 0.1, model::RANZ_MARSHALL));
        assert!(c.k < 0.0 && c.d > 40e-6);
    }

    #[test]
    fn abramzon_sirignano_evaporates_more_slowly_than_the_classical_model_in_a_stream() {
        let mut a = input(330.0, 800.0, 0.0, 50e-6, 1e-4, model::ABRAMZON_SIRIGNANO);
        a.re = 50.0;
        let mut r = a;
        r.model = model::RANZ_MARSHALL;
        let (sa, sr) = (evaporation_step(&a), evaporation_step(&r));
        assert!(sa.k < sr.k && sa.k > 0.5 * sr.k, "{} {}", sa.k, sr.k);
        assert!(sa.sh < sr.sh);
    }
}
