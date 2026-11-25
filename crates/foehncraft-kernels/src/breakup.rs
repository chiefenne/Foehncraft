//! Kernel bodies of the droplet breakup of the sprays (step 28, stage SP6,
//! D300 to D305): the fastest growing Kelvin-Helmholtz and Rayleigh-Taylor
//! waves on a drop, the exact relaxation of the Kelvin-Helmholtz radius, the
//! exact solution of the TAB distortion oscillator, and the drag factor of a
//! distorted drop.
//!
//! **Kelvin-Helmholtz (KH).** A liquid column (an injected blob, Reitz and
//! Diwakar 1987) or a drop of radius a in a gas stream of relative speed U
//! grows surface waves; Reitz (1987) fitted the wavelength Λ and growth rate
//! Ω of the fastest wave of the linear dispersion relation of a viscous jet
//! in an inviscid gas:
//!
//! Λ/a = 9.02 (1 + 0.45 Z^½)(1 + 0.4 T^0.7) / (1 + 0.87 We_g^1.67)^0.6,
//! Ω (ρ_l a³/σ)^½ = (0.34 + 0.38 We_g^1.5) / ((1 + Z)(1 + 1.4 T^0.6)),
//!
//! with We_g = ρ_g U² a/σ, the Ohnesorge number Z = μ_l/(ρ_l a σ)^½ (Reitz
//! writes it We_l^½/Re_l, which is the same and stays finite at U = 0) and the
//! Taylor number T = Z We_g^½. The drop sheds child drops of radius
//! r_c = B₀ Λ (B₀ = 0.61) and its radius relaxes to r_c on the breakup time
//! τ = 3.726 B₁ a/(Λ Ω): da/dt = −(a − r_c)/τ. With the gas frozen over a
//! sub-step this is linear and solved exactly, a(Δt) = r_c + (a − r_c)
//! e^(−Δt/τ) ([`kh_radius`]).
//!
//! **Rayleigh-Taylor (RT).** A drop decelerated by the drag at
//! g_t = 3 C_D ρ_g U²/(4 ρ_l d) is unstable on its front (Beale and Reitz
//! 1999; Patterson and Reitz 1998): the fastest wave has
//! Ω_RT = (2/(3√3 σ))^½ (g_t (ρ_l − ρ_g))^¾ / (ρ_l + ρ_g)^½ and wave number
//! K_RT = (g_t (ρ_l − ρ_g)/(3σ))^½; the drop breaks into drops of the size
//! Λ_RT = 2π C_RT/K_RT once the wave had the time C_τ/Ω_RT to grow and fits
//! on the drop.
//!
//! **TAB distortion.** The Taylor analogy (O'Rourke and Amsden 1987) treats
//! the distortion y of the drop's equator (y = x/(C_b r)) as a damped
//! oscillator driven by the gas,
//!
//! ÿ = C_F/C_b ρ_g U²/(ρ_l r²) − C_k σ/(ρ_l r³) y − C_d μ_l/(ρ_l r²) ẏ,
//!
//! C_F = 1/3, C_k = 8, C_d = 5, C_b = 1/2. With the coefficients frozen over
//! a sub-step the equation is linear with constant forcing and solved exactly
//! ([`tab_step`]), under-, critically and over-damped. Here y only shapes the
//! drag: Liu, Mather and Reitz (1993) let the drag coefficient of a distorted
//! drop grow from the sphere's to that of a disc, C_D = C_D,sphere (1 + 2.632
//! y), 0 ≤ y ≤ 1 ([`dynamic_drag`]).
//!
//! The bodies are plain f64 and use the in-house elementary functions of
//! `foehncraft_core::math`, so a value is the same on every platform (D93
//! (6)); each is a gather over one drop.

use foehncraft_core::math;

/// TAB constant C_F (O'Rourke and Amsden 1987).
pub const TAB_CF: f64 = 1.0 / 3.0;
/// TAB constant C_k.
pub const TAB_CK: f64 = 8.0;
/// TAB constant C_d.
pub const TAB_CD: f64 = 5.0;
/// TAB constant C_b.
pub const TAB_CB: f64 = 0.5;
/// Slope of the drag of a distorted drop (Liu, Mather and Reitz 1993):
/// C_D = C_D,sphere (1 + 2.632 y).
pub const DYNAMIC_DRAG_SLOPE: f64 = 2.632;

/// The fastest growing wave on a drop.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Wave {
    /// Wavelength (m).
    pub lambda: f64,
    /// Growth rate (1/s).
    pub omega: f64,
}

/// A drop in the gas, for [`kh_wave`] and [`tab_step`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DropInput {
    /// Drop radius (m).
    pub a: f64,
    /// Speed of the gas relative to the drop (m/s).
    pub u: f64,
    /// Liquid density (kg/m³).
    pub rho_l: f64,
    /// Gas density (kg/m³).
    pub rho_g: f64,
    /// Surface tension (N/m).
    pub sigma: f64,
    /// Liquid dynamic viscosity (Pa s).
    pub mu_l: f64,
}

impl DropInput {
    /// Gas Weber number ρ_g U² a/σ (on the radius).
    pub fn weber_gas(&self) -> f64 {
        self.rho_g * self.u * self.u * self.a / self.sigma
    }

    /// Ohnesorge number Z = μ_l/(ρ_l a σ)^½ (on the radius).
    pub fn ohnesorge(&self) -> f64 {
        self.mu_l / math::sqrt(self.rho_l * self.a * self.sigma)
    }
}

/// The fastest Kelvin-Helmholtz wave of Reitz's (1987) curve fits (module
/// documentation).
pub fn kh_wave(i: &DropInput) -> Wave {
    let we = i.weber_gas();
    let z = i.ohnesorge();
    let t = z * math::sqrt(we);
    let lambda = 9.02 * i.a * (1.0 + 0.45 * math::sqrt(z)) * (1.0 + 0.4 * math::pow(t, 0.7))
        / math::pow(1.0 + 0.87 * math::pow(we, 1.67), 0.6);
    let omega = (0.34 + 0.38 * math::pow(we, 1.5)) / ((1.0 + z) * (1.0 + 1.4 * math::pow(t, 0.6)))
        * math::sqrt(i.sigma / (i.rho_l * i.a * i.a * i.a));
    Wave { lambda, omega }
}

/// The Kelvin-Helmholtz breakup time τ = 3.726 B₁ a/(Λ Ω) (s).
pub fn kh_time(a: f64, w: &Wave, b1: f64) -> f64 {
    3.726 * b1 * a / (w.lambda * w.omega)
}

/// The radius after `dt` of a drop of radius `a` that relaxes to the child
/// radius `rc` on the time `tau`: the exact solution rc + (a − rc) e^(−dt/τ)
/// of da/dt = −(a − rc)/τ (only for rc < a: the drop sheds, it does not
/// grow).
pub fn kh_radius(a: f64, rc: f64, tau: f64, dt: f64) -> f64 {
    if rc >= a || tau.is_nan() || tau <= 0.0 {
        return a;
    }
    rc + (a - rc) * math::exp(-dt / tau)
}

/// The deceleration of a drop of diameter `d` by the drag (m/s²):
/// g_t = 3 C_D ρ_g U²/(4 ρ_l d).
pub fn drag_deceleration(cd: f64, rho_g: f64, u: f64, rho_l: f64, d: f64) -> f64 {
    0.75 * cd * rho_g * u * u / (rho_l * d)
}

/// The fastest Rayleigh-Taylor wave on the front of a drop decelerated at
/// `g_t` (m/s², module documentation): wavelength 2π/K_RT (before the size
/// constant C_RT) and growth rate Ω_RT. `None` without deceleration.
pub fn rt_wave(g_t: f64, rho_l: f64, rho_g: f64, sigma: f64) -> Option<Wave> {
    let drive = g_t * (rho_l - rho_g);
    if !(drive > 0.0 && sigma > 0.0) {
        return None;
    }
    let k = math::sqrt(drive / (3.0 * sigma));
    let omega = math::sqrt(2.0 / (3.0 * math::sqrt(3.0) * sigma)) * math::pow(drive, 0.75)
        / math::sqrt(rho_l + rho_g);
    Some(Wave {
        lambda: 2.0 * core::f64::consts::PI / k,
        omega,
    })
}

/// Breakup length of the liquid core (Levich 1962; Patterson and Reitz
/// 1998): L_b = C_bl (ρ_l/ρ_g)^½ d₀ (m).
pub fn breakup_length(c_bl: f64, rho_l: f64, rho_g: f64, d0: f64) -> f64 {
    c_bl * math::sqrt(rho_l / rho_g) * d0
}

/// The TAB distortion after a sub-step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TabStep {
    /// Distortion y at the end.
    pub y: f64,
    /// Its rate ẏ at the end (1/s).
    pub ydot: f64,
    /// The distortion the drop relaxes to, y_s = C_F/(C_k C_b) We_g
    /// (We_g on the radius).
    pub y_steady: f64,
}

/// The exact solution of the TAB oscillator over `dt` from (y₀, ẏ₀) with the
/// gas frozen (module documentation).
pub fn tab_step(i: &DropInput, y0: f64, ydot0: f64, dt: f64) -> TabStep {
    let r = i.a;
    // ÿ + 2β ẏ + k y = F
    let k = TAB_CK * i.sigma / (i.rho_l * r * r * r);
    let beta = 0.5 * TAB_CD * i.mu_l / (i.rho_l * r * r);
    let f = TAB_CF / TAB_CB * i.rho_g * i.u * i.u / (i.rho_l * r * r);
    let ys = f / k;
    let a0 = y0 - ys;
    let disc = k - beta * beta;
    let e = math::exp(-beta * dt);
    let (y, ydot) = if disc > 1e-12 * k {
        // under-damped
        let w = math::sqrt(disc);
        let b = (ydot0 + beta * a0) / w;
        let (s, c) = (math::sin(w * dt), math::cos(w * dt));
        (
            ys + e * (a0 * c + b * s),
            e * ((b * w - beta * a0) * c - (beta * b + w * a0) * s),
        )
    } else if disc < -1e-12 * k {
        // over-damped: roots s₁,₂ = −β ± √(β² − k)
        let q = math::sqrt(-disc);
        let s1 = -beta + q;
        let s2 = -beta - q;
        let c1 = (ydot0 - s2 * a0) / (s1 - s2);
        let c2 = a0 - c1;
        let (e1, e2) = (math::exp(s1 * dt), math::exp(s2 * dt));
        (ys + c1 * e1 + c2 * e2, s1 * c1 * e1 + s2 * c2 * e2)
    } else {
        // critically damped
        let b = ydot0 + beta * a0;
        (ys + (a0 + b * dt) * e, (b - beta * (a0 + b * dt)) * e)
    };
    TabStep {
        y,
        ydot,
        y_steady: ys,
    }
}

/// The drag factor C_D/C_D,sphere = 1 + 2.632 y of a drop of distortion y,
/// with y taken in [0, 1] (Liu, Mather and Reitz 1993).
pub fn dynamic_drag(y: f64) -> f64 {
    let y = if y > 1.0 {
        1.0
    } else if y > 0.0 {
        y
    } else {
        0.0
    };
    1.0 + DYNAMIC_DRAG_SLOPE * y
}

#[cfg(test)]
mod tests {
    extern crate std;
    use super::*;

    fn drop(a: f64, u: f64) -> DropInput {
        // a light hydrocarbon at about 360 K in a gas of 3.5 kg/m³
        DropInput {
            a,
            u,
            rho_l: 690.0,
            rho_g: 3.5,
            sigma: 0.0145,
            mu_l: 2.6e-4,
        }
    }

    /// The modified Bessel functions I₀, I₁, K₀, K₁ (Abramowitz and Stegun
    /// 9.8.1 to 9.8.8, polynomial approximations, |ε| < 2e-7).
    fn bessel(x: f64) -> (f64, f64, f64, f64) {
        let (i0, i1) = if x <= 3.75 {
            let t = (x / 3.75).powi(2);
            let i0 = 1.0
                + t * (3.5156229
                    + t * (3.0899424
                        + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
            let i1 = x
                * (0.5
                    + t * (0.87890594
                        + t * (0.51498869
                            + t * (0.15084934
                                + t * (0.02658733 + t * (0.00301532 + t * 0.00032411))))));
            (i0, i1)
        } else {
            let t = 3.75 / x;
            let f = x.exp() / x.sqrt();
            let i0 = f
                * (0.39894228
                    + t * (0.01328592
                        + t * (0.00225319
                            + t * (-0.00157565
                                + t * (0.00916281
                                    + t * (-0.02057706
                                        + t * (0.02635537
                                            + t * (-0.01647633 + t * 0.00392377))))))));
            let i1 = f
                * (0.39894228
                    + t * (-0.03988024
                        + t * (-0.00362018
                            + t * (0.00163801
                                + t * (-0.01031555
                                    + t * (0.02282967
                                        + t * (-0.02895312
                                            + t * (0.01787654 - t * 0.00420059))))))));
            (i0, i1)
        };
        let (k0, k1) = if x <= 2.0 {
            let t = x * x / 4.0;
            let l = (x / 2.0).ln();
            let k0 = -l * i0
                + (-0.57721566
                    + t * (0.42278420
                        + t * (0.23069756
                            + t * (0.03488590
                                + t * (0.00262698 + t * (0.00010750 + t * 0.0000074))))));
            let k1 = l * i1
                + (1.0
                    + t * (0.15443144
                        + t * (-0.67278579
                            + t * (-0.18156897
                                + t * (-0.01919402 + t * (-0.00110404 - t * 0.00004686))))))
                    / x;
            (k0, k1)
        } else {
            let t = 2.0 / x;
            let f = (-x).exp() / x.sqrt();
            let k0 = f
                * (1.25331414
                    + t * (-0.07832358
                        + t * (0.02189568
                            + t * (-0.01062446
                                + t * (0.00587872 + t * (-0.00251540 + t * 0.00053208))))));
            let k1 = f
                * (1.25331414
                    + t * (0.23498619
                        + t * (-0.03655620
                            + t * (0.01504268
                                + t * (-0.00780353 + t * (0.00325614 - t * 0.00068245))))));
            (k0, k1)
        };
        (i0, i1, k0, k1)
    }

    /// The largest growth rate of the inviscid dispersion relation of a
    /// liquid jet in a gas stream (Rayleigh's capillary term plus the
    /// aerodynamic term of the gas, the Z = 0 limit of Reitz 1987 eq. 1):
    /// ω² = σk/(ρ_l a²)(1 − k²a²) I₁/I₀ + ρ_g U² k²/ρ_l (I₁ K₀)/(I₀ K₁),
    /// with the wavelength 2π/k of the maximum.
    fn inviscid_max(i: &DropInput) -> (f64, f64) {
        let a = i.a;
        let mut best = (0.0, 0.0);
        let n = 20000;
        for j in 1..n {
            // ka from 0 to a few hundred on a log scale
            let ka = 1e-3 * (1e6f64).powf(j as f64 / n as f64);
            let k = ka / a;
            let (i0, i1, k0, k1) = bessel(ka);
            let w2 = i.sigma * k / (i.rho_l * a * a) * (1.0 - ka * ka) * i1 / i0
                + i.rho_g * i.u * i.u * k * k / i.rho_l * (i1 * k0) / (i0 * k1);
            if w2 > best.0 {
                best = (w2, k);
            }
        }
        (best.0.sqrt(), 2.0 * std::f64::consts::PI / best.1)
    }

    /// The fits against the inviscid dispersion relation at Z → 0. The
    /// wavelength fit holds within 5 % at every gas Weber number; the growth
    /// rate fit is the sum of the two limits (Rayleigh's 0.343 and the
    /// aerodynamic (4/27)^½ We^1.5 = 0.385 We^1.5) and matches within 3 %
    /// below We = 0.1 and above We = 50, but overshoots by up to 46 % around
    /// We = 2 (own derivation, D301): the Spray G blobs (We ≈ 1000) lie in
    /// the accurate range.
    #[test]
    fn kh_fits_match_the_dispersion_relation() {
        let mut worst_mid = 0.0f64;
        for we in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 200.0, 1000.0] {
            let mut i = drop(80e-6, 0.0);
            i.mu_l = 1e-12; // Z → 0
            i.u = (we * i.sigma / (i.rho_g * i.a)).sqrt();
            let fit = kh_wave(&i);
            let (omega, lambda) = inviscid_max(&i);
            let eo = (fit.omega - omega).abs() / omega;
            let el = (fit.lambda - lambda).abs() / lambda;
            std::println!(
                "We {we}: Omega fit {:.4e} exact {:.4e} ({:.1}%), Lambda fit {:.4e} exact {:.4e} ({:.1}%)",
                fit.omega,
                omega,
                100.0 * eo,
                fit.lambda,
                lambda,
                100.0 * el
            );
            assert!(el < 0.05, "wavelength off by {el} at We {we}");
            if !(0.1..50.0).contains(&we) {
                assert!(eo < 0.03, "growth rate off by {eo} at We {we}");
            } else {
                worst_mid = worst_mid.max(eo);
            }
        }
        assert!(
            worst_mid < 0.47,
            "growth rate off by {worst_mid} between the limits"
        );
    }

    /// At U = 0 the fit is Rayleigh's capillary instability of a column:
    /// Ω = 0.34 (σ/(ρ_l a³))^½ and Λ = 9.02 a (Rayleigh's maximum: ka = 0.697,
    /// Λ = 9.01 a, Ω = 0.343).
    #[test]
    fn kh_at_rest_is_rayleigh() {
        let mut i = drop(50e-6, 0.0);
        i.mu_l = 0.0;
        let w = kh_wave(&i);
        let s = (i.sigma / (i.rho_l * i.a.powi(3))).sqrt();
        assert!((w.lambda / i.a - 9.02).abs() < 1e-12);
        assert!((w.omega / s - 0.34).abs() < 1e-12);
    }

    #[test]
    fn kh_radius_relaxes_exactly() {
        let (a, rc, tau) = (80e-6, 5e-6, 2e-5);
        let one = kh_radius(a, rc, tau, 3e-5);
        // two half steps give the same as one step (linear, exact)
        let half = kh_radius(kh_radius(a, rc, tau, 1.5e-5), rc, tau, 1.5e-5);
        assert!((one - half).abs() < 1e-18);
        assert!((one - (rc + (a - rc) * (-1.5f64).exp())).abs() < 1e-17);
        // no growth
        assert_eq!(kh_radius(4e-6, 5e-6, tau, 1.0), 4e-6);
    }

    #[test]
    fn rt_wave_matches_the_formulae() {
        let (g, rl, rg, s) = (1e7, 690.0, 3.5, 0.0145);
        let w = rt_wave(g, rl, rg, s).unwrap();
        let k = (g * (rl - rg) / (3.0 * s)).sqrt();
        let om = (2.0 / (3.0 * 3f64.sqrt() * s) * (g * (rl - rg)).powf(1.5) / (rl + rg)).sqrt();
        assert!((w.lambda - 2.0 * std::f64::consts::PI / k).abs() / w.lambda < 1e-13);
        assert!((w.omega - om).abs() / om < 1e-13);
        assert!(rt_wave(-1.0, rl, rg, s).is_none());
    }

    /// The exact TAB solution against a fine Runge-Kutta integration, in the
    /// under-, critically and over-damped regimes, and two half steps against
    /// one step.
    #[test]
    fn tab_step_is_exact() {
        let cases = [
            drop(40e-6, 30.0), // under-damped
            DropInput {
                mu_l: 0.5,
                ..drop(5e-6, 30.0)
            }, // over-damped
        ];
        for i in cases {
            let k = TAB_CK * i.sigma / (i.rho_l * i.a.powi(3));
            let beta = 0.5 * TAB_CD * i.mu_l / (i.rho_l * i.a * i.a);
            let f = TAB_CF / TAB_CB * i.rho_g * i.u * i.u / (i.rho_l * i.a * i.a);
            let dt = 3.0 / k.sqrt();
            let (y0, v0) = (0.1, -2.0 * k.sqrt());
            let s = tab_step(&i, y0, v0, dt);
            // RK4
            let n = 200_000;
            let h = dt / n as f64;
            let (mut y, mut v) = (y0, v0);
            let acc = |y: f64, v: f64| f - k * y - 2.0 * beta * v;
            for _ in 0..n {
                let (k1y, k1v) = (v, acc(y, v));
                let (k2y, k2v) = (v + 0.5 * h * k1v, acc(y + 0.5 * h * k1y, v + 0.5 * h * k1v));
                let (k3y, k3v) = (v + 0.5 * h * k2v, acc(y + 0.5 * h * k2y, v + 0.5 * h * k2v));
                let (k4y, k4v) = (v + h * k3v, acc(y + h * k3y, v + h * k3v));
                y += h / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
                v += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
            }
            let scale = 1.0 + s.y_steady.abs();
            assert!((s.y - y).abs() < 1e-8 * scale, "y {} vs {}", s.y, y);
            assert!(
                (s.ydot - v).abs() < 1e-8 * scale * k.sqrt(),
                "ydot {} vs {}",
                s.ydot,
                v
            );
            let h1 = tab_step(&i, y0, v0, 0.5 * dt);
            let h2 = tab_step(&i, h1.y, h1.ydot, 0.5 * dt);
            assert!((h2.y - s.y).abs() < 1e-10 * scale);
        }
        // critically damped: the limit of both neighbours
        let mut i = drop(10e-6, 20.0);
        let k = TAB_CK * i.sigma / (i.rho_l * i.a.powi(3));
        // β² = k
        i.mu_l = 2.0 * k.sqrt() * i.rho_l * i.a * i.a / TAB_CD;
        let dt = 1.0 / k.sqrt();
        let c = tab_step(&i, 0.2, 0.0, dt);
        let mut lo = i;
        lo.mu_l *= 1.0 - 1e-6;
        let mut hi = i;
        hi.mu_l *= 1.0 + 1e-6;
        let (a, b) = (tab_step(&lo, 0.2, 0.0, dt), tab_step(&hi, 0.2, 0.0, dt));
        assert!((c.y - a.y).abs() < 1e-5 && (c.y - b.y).abs() < 1e-5);
    }

    /// The steady distortion is We_g/12 (C_F/(C_k C_b) = 1/12, We_g on the
    /// radius): y = 1 at We_g = 12 in a steady stream, at 6 with the
    /// overshoot of a drop that starts at rest.
    #[test]
    fn tab_steady_distortion() {
        let i = drop(40e-6, 30.0);
        let s = tab_step(&i, 0.0, 0.0, 1.0);
        assert!((s.y_steady - i.weber_gas() / 12.0).abs() < 1e-12 * s.y_steady);
        assert!((s.y - s.y_steady).abs() < 1e-9);
    }

    #[test]
    fn dynamic_drag_is_bounded() {
        assert_eq!(dynamic_drag(-0.3), 1.0);
        assert_eq!(dynamic_drag(0.0), 1.0);
        assert!((dynamic_drag(0.5) - 2.316).abs() < 1e-12);
        assert!((dynamic_drag(7.0) - 3.632).abs() < 1e-12);
    }
}
