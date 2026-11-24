//! Exact solutions of the verification cases (step 29, section 3).
//!
//! - the exponential integrals `E_n(x)` (Abramowitz & Stegun, *Handbook of
//!   Mathematical Functions*, 1964: the power series 5.1.12 for x ≤ 1, the
//!   continued fraction 5.1.22 for x > 1);
//! - the isothermal grey slab between cold black plates: the flux leaving
//!   each face, `σT⁴ (1 − 2 E₃(τ))`, and the net source in a layer of the
//!   slab, from the incident radiation `G(x) = 2σT⁴/π · π (2 − E₂(κx) −
//!   E₂(κ(L − x)))` (Modest & Mazumder, *Radiative Heat Transfer*, 4th ed.,
//!   2021, the one-dimensional slab);
//! - the configuration factor between two directly opposed parallel
//!   rectangles (Howell, Mengüç, Daun & Siegel, *Thermal Radiation Heat
//!   Transfer*, 7th ed., 2021, catalogue of configuration factors).

use std::f64::consts::PI;

const EULER: f64 = 0.577_215_664_901_532_9;

/// The exponential integral `E_n(x) = ∫₁^∞ e^{−xt} / tⁿ dt`, for n ≥ 1 and
/// x ≥ 0 (x > 0 for n = 1), to about 1e-15 relative.
pub fn expint(n: u32, x: f64) -> f64 {
    assert!(n >= 1 && x >= 0.0, "expint: n >= 1, x >= 0");
    let nm1 = f64::from(n - 1);
    if x == 0.0 {
        assert!(n > 1, "E_1(0) is infinite");
        return 1.0 / nm1;
    }
    if x > 1.0 {
        // continued fraction (modified Lentz)
        let tiny = 1e-300;
        let mut b = x + f64::from(n);
        let mut c = 1.0 / tiny;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..1000 {
            let i = f64::from(i);
            let an = -i * (nm1 + i);
            b += 2.0;
            d = 1.0 / (an * d + b);
            c = b + an / c;
            let del = c * d;
            h *= del;
            if (del - 1.0).abs() < 1e-16 {
                break;
            }
        }
        return h * (-x).exp();
    }
    // power series
    let mut ans = if n > 1 { 1.0 / nm1 } else { -x.ln() - EULER };
    let mut fact = 1.0;
    for i in 1..1000u32 {
        fact *= -x / f64::from(i);
        let del = if i != n - 1 {
            -fact / (f64::from(i) - nm1)
        } else {
            let psi = -EULER + (1..n).map(|k| 1.0 / f64::from(k)).sum::<f64>();
            fact * (-x.ln() + psi)
        };
        ans += del;
        if del.abs() < ans.abs() * 1e-17 {
            break;
        }
    }
    ans
}

/// Flux leaving each face of an isothermal grey slab of optical thickness
/// `tau` between cold black plates, divided by σT⁴.
pub fn slab_wall_flux(tau: f64) -> f64 {
    1.0 - 2.0 * expint(3, tau)
}

/// Net source (absorbed − emitted) averaged over the layer `x0 < x < x1`
/// of an isothermal grey slab of thickness `l` and absorption coefficient
/// `kappa` between cold black plates, divided by σT⁴ (unit 1/m).
pub fn slab_source(kappa: f64, l: f64, x0: f64, x1: f64) -> f64 {
    // mean of E₂(κx) over the layer: (E₃(κx0) − E₃(κx1)) / (κ (x1 − x0))
    let mean = |a: f64, b: f64| (expint(3, kappa * a) - expint(3, kappa * b)) / (kappa * (b - a));
    -2.0 * kappa * (mean(x0, x1) + mean(l - x1, l - x0))
}

/// Configuration factor from a rectangle `a × b` to the directly opposed
/// parallel rectangle of the same size at distance `c`.
pub fn parallel_rectangles(a: f64, b: f64, c: f64) -> f64 {
    let x = a / c;
    let y = b / c;
    let (x2, y2) = (x * x, y * y);
    let s1 = ((1.0 + x2) * (1.0 + y2) / (1.0 + x2 + y2)).sqrt().ln();
    let s2 = x * (1.0 + y2).sqrt() * (x / (1.0 + y2).sqrt()).atan();
    let s3 = y * (1.0 + x2).sqrt() * (y / (1.0 + x2).sqrt()).atan();
    let s4 = x * x.atan() + y * y.atan();
    2.0 / (PI * x * y) * (s1 + s2 + s3 - s4)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exponential_integrals() {
        // tabulated values (A&S table 5.1 and the closed forms)
        assert!((expint(1, 1.0) - 0.219_383_934_395_520_3).abs() < 1e-15);
        assert!((expint(2, 1.0) - 0.148_495_506_775_922).abs() < 1e-15);
        assert!((expint(3, 0.0) - 0.5).abs() < 1e-16);
        assert!((expint(1, 0.5) - 0.559_773_594_776_160_8).abs() < 1e-15);
        // recurrence n E_{n+1}(x) = e^{-x} - x E_n(x), across both branches
        for &x in &[0.01, 0.3, 0.99, 1.0, 1.01, 2.5, 7.0, 20.0] {
            for n in 1..4 {
                let lhs = f64::from(n) * expint(n + 1, x);
                let rhs = (-x).exp() - x * expint(n, x);
                assert!(
                    (lhs - rhs).abs() < 1e-14 * lhs.abs().max(1e-300) + 1e-16,
                    "{n} {x}"
                );
            }
        }
    }

    #[test]
    fn slab_limits() {
        // thin: flux ~ 2τ (optically thin emission); thick: flux -> 1
        assert!((slab_wall_flux(1e-6) / 2e-6 - 1.0).abs() < 1e-4);
        assert!((slab_wall_flux(30.0) - 1.0).abs() < 1e-12);
        // the sources of all layers balance the two wall fluxes
        let (k, l, n) = (2.0, 1.0, 50);
        let total: f64 = (0..n)
            .map(|i| {
                let (a, b) = (l * i as f64 / n as f64, l * (i + 1) as f64 / n as f64);
                slab_source(k, l, a, b) * (b - a)
            })
            .sum();
        assert!(
            (total + 2.0 * slab_wall_flux(k * l)).abs() < 1e-12,
            "{total}"
        );
    }

    #[test]
    fn opposed_squares() {
        // unit squares at distance 1: F = 0.19982 (catalogue value)
        assert!((parallel_rectangles(1.0, 1.0, 1.0) - 0.199_825).abs() < 1e-5);
        // close plates see each other entirely
        assert!(parallel_rectangles(1.0, 1.0, 1e-4) > 0.999);
    }
}
