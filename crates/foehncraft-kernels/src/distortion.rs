//! The drag, heat and mass transfer of a distorted, viscous drop: the
//! Corrected Distortion model of Nguyen, Dahms, Pickett and Tagliante
//! (2022), Int. J. Multiphase Flow 148, 103927 (doi
//! 10.1016/j.ijmultiphaseflow.2021.103927; equation numbers below are
//! theirs), for ECN Spray A (SA4, D316).
//!
//! A drop of equivalent radius r and TAB distortion y (O'Rourke and Amsden
//! 1987; y > 0 oblate, y < 0 prolate) is a spheroid with the axis
//! L = r (1 − C_b y) along the flow and the crosswise half-axis
//! H = r/√(1 − C_b y), so that L H² = r³ (eqs. 24, 28, 29; C_b = 1/2). Its
//! sphericity φ = 4πr²/S and crosswise sphericity φ⊥ = (r/H)² (eqs. 26, 27)
//! follow from the spheroid surface S (eqs. 30, 31).
//!
//! - **Drag** (eq. 17, 18, 23): the drag coefficient of a viscous sphere of
//!   Feng and Michaelides (2001a) for the viscosity ratio λ = μ_l/μ_g, times
//!   the shape ratio of Richter and Nikrityuk (2012) in l/d = 1 − C_b y.
//! - **Heat and mass transfer** (eqs. 19, 20, 25, 32): the Sherwood (or
//!   Nusselt) number of a viscous sphere of Feng and Michaelides (2001b),
//!   times the shape ratio of Richter and Nikrityuk in φ and φ⊥, times the
//!   surface ratio S/(π d²) = 1/φ.
//!
//! The functions here return these corrections as ratios to the solid
//! sphere of the same correlations (eqs. 18c, 20c), so that a code applies
//! them as factors on its own sphere laws (Schiller-Naumann drag, the
//! Ranz-Marshall or Abramzon-Sirignano Sherwood and Nusselt numbers with
//! their blowing corrections). At y = 0 and λ → ∞ every factor is 1. The
//! correlations of Feng and Michaelides hold for 5 < Re < 1000; Re is held
//! in that range for their ratio (below 5 the bubble term turns negative).
//! The shape ratios take Re ≤ 1000 (sub-critical, as fitted).

use foehncraft_core::math;

/// TAB constant C_b, the axis ratio of a distorted drop (O'Rourke and
/// Amsden 1987; eq. 24).
pub const C_B: f64 = 0.5;
/// The range of Re of the correlations of Feng and Michaelides (2001).
pub const RE_MIN: f64 = 5.0;
/// Upper end of the Reynolds number of the correlations.
pub const RE_MAX: f64 = 1000.0;

/// The geometry of a drop of distortion y (eqs. 24 to 31).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Shape {
    /// The length along the flow over the equivalent diameter, l/d = 1 − C_b y.
    pub l_over_d: f64,
    /// Sphericity φ = 4πr²/S.
    pub phi: f64,
    /// Crosswise sphericity φ⊥ = (r/H)² = 1 − C_b y.
    pub phi_perp: f64,
    /// Surface over the surface of the equivalent sphere, S/(π d²) = 1/φ.
    pub area: f64,
}

/// atanh(e)/e, 1 at e = 0.
fn atanh_over(e: f64) -> f64 {
    if e < 1e-4 {
        let e2 = e * e;
        1.0 + e2 * (1.0 / 3.0 + e2 / 5.0)
    } else {
        0.5 * math::ln((1.0 + e) / (1.0 - e)) / e
    }
}

/// asin(e)/e, 1 at e = 0.
fn asin_over(e: f64) -> f64 {
    if e < 1e-4 {
        let e2 = e * e;
        1.0 + e2 * (1.0 / 6.0 + 3.0 / 40.0 * e2)
    } else {
        math::asin(e) / e
    }
}

/// The spheroid of a drop of distortion `y`, taken in [−1, 1] (the TAB
/// drop breaks at y = 1), with the axis constant `cb`.
pub fn shape(y: f64, cb: f64) -> Shape {
    let y = y.clamp(-1.0, 1.0);
    let a = 1.0 - cb * y;
    // per r: L = a, H = 1/√a
    let l = a;
    let h2 = 1.0 / a;
    let h = math::sqrt(h2);
    let (lo, hi) = if h < l { (h, l) } else { (l, h) };
    let q = lo / hi;
    let e = math::sqrt((1.0 - q * q).max(0.0));
    // eq. 30, per r²: oblate (H the equatorial radius) or prolate
    let s = if y >= 0.0 {
        2.0 * core::f64::consts::PI * h2 * (1.0 + (1.0 - e * e) * atanh_over(e))
    } else {
        2.0 * core::f64::consts::PI * h2 * (1.0 + l / h * asin_over(e))
    };
    let phi = 4.0 * core::f64::consts::PI / s;
    Shape {
        l_over_d: a,
        phi,
        phi_perp: a,
        area: 1.0 / phi,
    }
}

fn clamp_re(re: f64) -> f64 {
    re.clamp(RE_MIN, RE_MAX)
}

/// x^p for x > 0.
fn pow(x: f64, p: f64) -> f64 {
    math::exp(p * math::ln(x))
}

/// The drag coefficient of a solid sphere of the model, eq. 18c,
/// 24/Re (1 + Re^⅔/6).
pub fn drag_solid(re: f64) -> f64 {
    24.0 / re * (1.0 + pow(re, 2.0 / 3.0) / 6.0)
}

/// The drag coefficient of a viscous sphere of viscosity ratio `lambda`
/// (eqs. 17, 18; Feng and Michaelides 2001a) at `re` in [5, 1000].
pub fn drag_viscous(re: f64, lambda: f64) -> f64 {
    let re = clamp_re(re);
    let bubble = 48.0 / re * (1.0 + 2.21 / math::sqrt(re) - 2.14 / re);
    let two = 17.0 * pow(re, -2.0 / 3.0);
    if lambda <= 2.0 {
        let l = lambda.max(0.0);
        (2.0 - l) / 2.0 * bubble + 4.0 * l / (6.0 + l) * two
    } else {
        4.0 / (2.0 + lambda) * two + (lambda - 2.0) / (lambda + 2.0) * drag_solid(re)
    }
}

/// The drag of a viscous sphere over that of a solid sphere at `re`.
pub fn drag_viscous_ratio(re: f64, lambda: f64) -> f64 {
    drag_viscous(re, lambda) / drag_solid(clamp_re(re))
}

/// The shape ratio of the drag C_D/C_D,vis,sph of Richter and Nikrityuk
/// (2012), eq. 23: (0.21 + 20/Re (l/d)^0.58 + 6.9/√Re (l/d)^−1.4) /
/// (0.21 + 20/Re + 6.9/√Re).
pub fn drag_shape_ratio(re: f64, l_over_d: f64) -> f64 {
    let re = re.clamp(1e-12, RE_MAX);
    let a = 20.0 / re;
    let b = 6.9 / math::sqrt(re);
    (0.21 + a * pow(l_over_d, 0.58) + b * pow(l_over_d, -1.4)) / (0.21 + a + b)
}

/// The three Sherwood (or Nusselt) numbers of Feng and Michaelides (2001b),
/// eq. 20, of the bubble, of λ = 2 and of the solid sphere, at `re` and
/// the Péclet number `pe` (Sc Re or Pr Re).
pub fn transfer_parts(re: f64, pe: f64) -> [f64; 3] {
    let re = clamp_re(re);
    let pe = pe.max(1e-30);
    let f = 0.61 * re / (re + 21.0);
    let bubble = 0.651 * math::sqrt(pe) * (1.032 + f) + (1.60 - f);
    let r287 = pow(re, 0.287);
    let two = 0.64 * pow(pe, 0.43) * (1.0 + 0.233 * r287) + 1.41 - 0.15 * r287;
    let solid = 0.852 * pow(pe, 1.0 / 3.0) * (1.0 + 0.233 * r287) + 1.3 - 0.182 * pow(re, 0.355);
    [bubble, two, solid]
}

/// The Sherwood (or Nusselt) number of a viscous sphere over that of a
/// solid sphere (eq. 19 over eq. 20c) at `re`, the Schmidt (or Prandtl)
/// number `x` and the viscosity ratio `lambda`.
pub fn transfer_viscous_ratio(re: f64, x: f64, lambda: f64) -> f64 {
    let rc = clamp_re(re);
    let [bubble, two, solid] = transfer_parts(rc, x * rc);
    let vis = if lambda <= 2.0 {
        let l = lambda.max(0.0);
        (2.0 - l) / 2.0 * bubble + 4.0 * l / (6.0 + l) * two
    } else {
        4.0 / (2.0 + lambda) * two + (lambda - 2.0) / (lambda + 2.0) * solid
    };
    vis / solid
}

/// The shape ratio Sh/Sh_vis,sph of Richter and Nikrityuk (2012), eq. 25,
/// at `re`, the Schmidt (or Prandtl) number `x`, φ and φ⊥.
pub fn transfer_shape_ratio(re: f64, x: f64, phi: f64, phi_perp: f64) -> f64 {
    let re = re.clamp(0.0, RE_MAX);
    let x3 = pow(x.max(1e-30), 1.0 / 3.0);
    let a = 0.55 * x3 * math::sqrt(re);
    let b = if re > 0.0 {
        0.014 * x3 * pow(re, 2.0 / 3.0)
    } else {
        0.0
    };
    (1.76 + a * phi * pow(phi_perp, 0.075) + b * pow(phi / phi_perp, 7.2)) / (1.76 + a + b)
}

/// The drag factor C_D/C_D,solid sphere of a drop of distortion `y` and
/// viscosity ratio `lambda` at `re` (eq. 17 times eq. 23).
pub fn drag_factor(re: f64, lambda: f64, y: f64) -> f64 {
    let s = shape(y, C_B);
    drag_viscous_ratio(re, lambda) * drag_shape_ratio(re, s.l_over_d)
}

/// The factor on the transfer rate (Sherwood with `x` = Sc, Nusselt with
/// `x` = Pr) of a drop of distortion `y` and viscosity ratio `lambda` at
/// `re`: eq. 19 over 20c, times eq. 25, times the surface ratio of eq. 32.
pub fn transfer_factor(re: f64, x: f64, lambda: f64, y: f64) -> f64 {
    let s = shape(y, C_B);
    transfer_viscous_ratio(re, x, lambda) * transfer_shape_ratio(re, x, s.phi, s.phi_perp) * s.area
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The surface of a spheroid of revolution by quadrature of the
    /// generating ellipse (polar half-axis c, equatorial a): an independent
    /// check of eq. 30.
    fn quad_surface(a: f64, c: f64) -> f64 {
        let n = 200_000;
        let mut s = 0.0;
        for k in 0..n {
            let t = core::f64::consts::PI * (k as f64 + 0.5) / n as f64;
            let (x, z) = (a * t.sin(), c * t.cos());
            let (dx, dz) = (a * t.cos(), -c * t.sin());
            s += 2.0 * core::f64::consts::PI * x * (dx * dx + dz * dz).sqrt();
            let _ = z;
        }
        s * core::f64::consts::PI / n as f64
    }

    #[test]
    fn the_sphere_is_the_sphere() {
        let s = shape(0.0, C_B);
        assert_eq!((s.l_over_d, s.phi_perp), (1.0, 1.0));
        assert!((s.phi - 1.0).abs() < 1e-15 && (s.area - 1.0).abs() < 1e-15);
        assert!((drag_shape_ratio(100.0, 1.0) - 1.0).abs() < 1e-15);
        assert!((transfer_shape_ratio(100.0, 2.0, 1.0, 1.0) - 1.0).abs() < 1e-15);
        // a solid drop (λ → ∞) is the solid sphere
        assert!((drag_viscous_ratio(100.0, 1e12) - 1.0).abs() < 1e-10);
        assert!((transfer_viscous_ratio(100.0, 2.0, 1e12) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn the_spheroid_surface_and_volume() {
        // eq. 30 against quadrature of the surface of revolution, both
        // shapes; the volume L H² = r³ by construction
        for y in [-1.0, -0.4, 0.3, 1.0] {
            let s = shape(y, C_B);
            let l = 1.0 - C_B * y;
            let h = 1.0 / l.sqrt();
            assert!((l * h * h - 1.0).abs() < 1e-14);
            let q = quad_surface(h, l) / (4.0 * core::f64::consts::PI);
            assert!((s.area / q - 1.0).abs() < 1e-9, "{y}: {} {q}", s.area);
        }
        // the most distorted drop has about 20 % more surface (the paper,
        // section 4.3): 22.7 %
        assert!((shape(1.0, C_B).area - 1.2272).abs() < 1e-4);
    }

    #[test]
    fn the_viscous_sphere_branches_meet_at_lambda_two() {
        for re in [5.0, 50.0, 500.0] {
            let lo = drag_viscous(re, 2.0);
            let hi = drag_viscous(re, 2.0 + 1e-12);
            assert!((lo / hi - 1.0).abs() < 1e-10);
            let lo = transfer_viscous_ratio(re, 2.0, 2.0);
            let hi = transfer_viscous_ratio(re, 2.0, 2.0 + 1e-12);
            assert!((lo / hi - 1.0).abs() < 1e-10);
            // from λ = 2 on (Spray A: λ of 2 to 30) the internal
            // circulation lowers the drag and raises the transfer, less so
            // the more viscous the drop (below λ = 2 the bubble term of the
            // correlation is not monotone at low Re)
            let mut last = (0.0, f64::INFINITY);
            for lambda in [2.0, 10.0, 100.0] {
                let d = drag_viscous_ratio(re, lambda);
                let t = transfer_viscous_ratio(re, 2.0, lambda);
                assert!(d < 1.0 && d > last.0, "{re} {lambda}: {d}");
                assert!(t > 1.0 && t < last.1, "{re} {lambda}: {t}");
                last = (d, t);
            }
        }
    }

    #[test]
    fn a_flattened_drop_drags_more_and_evaporates_faster() {
        // the trends of figs. 18 and 19 of the paper (λ from 2 to 10,
        // Spray A drops): at y = 1 the drag 1.1 to 2 times and the mass
        // transfer about 1.5 to 4 times that of a sphere (fig. 18 is
        // normalised by the Frossling number, so only roughly the same
        // ratio); a prolate drop less drag
        for lambda in [2.5, 10.0] {
            for re in [50.0, 200.0] {
                let d1 = drag_factor(re, lambda, 1.0);
                let m1 = transfer_factor(re, 2.0, lambda, 1.0);
                assert!(d1 > 1.1 && d1 < 2.0, "{lambda} {re}: {d1}");
                assert!(m1 > 1.5 && m1 < 4.5, "{lambda} {re}: {m1}");
                assert!(drag_factor(re, lambda, -0.5) < drag_factor(re, lambda, 0.0));
                let mut last = 0.0;
                for y in [0.0, 0.25, 0.5, 0.75, 1.0] {
                    let m = transfer_factor(re, 2.0, lambda, y);
                    assert!(m > last);
                    last = m;
                }
            }
        }
    }
}
