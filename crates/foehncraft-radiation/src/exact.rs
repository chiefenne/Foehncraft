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
//!   Transfer*, 7th ed., 2021, catalogue of configuration factors);
//! - the flux into a cold black wall of a cube filled with an isothermal
//!   grey gas (case 5): the integral of the gas emission over the other
//!   five walls, `q(p) = σT⁴ Σ ∫ (1 − e^{−κr}) cos θ_p cos θ_q / (π r²) dA`,
//!   by adaptive Gauss–Legendre quadrature to about 1e-10, and its mean
//!   over a part of the wall. The exact reference of the cube case, so no
//!   published numbers are needed.

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

/// Nodes and weights of the 6-point Gauss–Legendre rule on [0, 1].
const GL6: [(f64, f64); 6] = [
    (0.033_765_242_898_423_975, 0.085_662_246_189_585_17),
    (0.169_395_306_766_867_74, 0.180_380_786_524_069_3),
    (0.380_690_406_958_401_5, 0.233_956_967_286_345_53),
    (0.619_309_593_041_598_5, 0.233_956_967_286_345_53),
    (0.830_604_693_233_132_3, 0.180_380_786_524_069_3),
    (0.966_234_757_101_576, 0.085_662_246_189_585_17),
];

/// The 6 × 6 Gauss rule of `f` on the rectangle [a0, a1] × [b0, b1].
fn gauss2(f: &impl Fn(f64, f64) -> f64, a0: f64, a1: f64, b0: f64, b1: f64) -> f64 {
    let (da, db) = (a1 - a0, b1 - b0);
    let mut s = 0.0;
    for &(u, wu) in &GL6 {
        for &(v, wv) in &GL6 {
            s += wu * wv * f(a0 + u * da, b0 + v * db);
        }
    }
    s * da * db
}

/// Adaptive 2D quadrature: splits a rectangle into four until the rule
/// on the four parts agrees with the rule on the whole to `tol`.
fn adapt(f: &impl Fn(f64, f64) -> f64, r: [f64; 4], whole: f64, tol: f64, depth: u32) -> f64 {
    let [a0, a1, b0, b1] = r;
    let (am, bm) = (0.5 * (a0 + a1), 0.5 * (b0 + b1));
    let parts = [
        [a0, am, b0, bm],
        [am, a1, b0, bm],
        [a0, am, bm, b1],
        [am, a1, bm, b1],
    ];
    let q: Vec<f64> = parts
        .iter()
        .map(|p| gauss2(f, p[0], p[1], p[2], p[3]))
        .collect();
    let sum: f64 = q.iter().sum();
    if (sum - whole).abs() <= tol || depth >= 12 {
        return sum;
    }
    parts
        .iter()
        .zip(&q)
        .map(|(p, &w)| adapt(f, *p, w, 0.25 * tol, depth + 1))
        .sum()
}

fn integrate(f: impl Fn(f64, f64) -> f64, tol: f64) -> f64 {
    let whole = gauss2(&f, 0.0, 1.0, 0.0, 1.0);
    adapt(&f, [0.0, 1.0, 0.0, 1.0], whole, tol, 0)
}

/// Flux into a cold black wall of the unit cube at the point `(x, y)` of
/// the wall (0 < x, y < 1), with the cube filled with an isothermal grey
/// gas of optical side `kl` (κ times the side), divided by σT⁴. For a cube
/// of side `L`, take `kl = κ L`; the flux does not depend on L otherwise.
pub fn cube_wall_flux(kl: f64, x: f64, y: f64) -> f64 {
    let emit = |r: f64| -> f64 {
        let t = kl * r;
        if t < 1e-5 {
            t * (1.0 - t * (0.5 - t / 6.0))
        } else {
            1.0 - (-t).exp()
        }
    };
    let tol = 1e-12;
    // opposite wall z = 1: cos θ_p = cos θ_q = 1/r
    let top = integrate(
        |u, v| {
            let r2 = (u - x).powi(2) + (v - y).powi(2) + 1.0;
            emit(r2.sqrt()) / (PI * r2 * r2)
        },
        tol,
    );
    // a side wall at distance h from the point, coordinates (s, w): s
    // along the wall, w the height; cos θ_p = w / r, cos θ_q = h / r
    let side = |h: f64, c: f64| {
        integrate(
            |s, w| {
                let r2 = h * h + (s - c).powi(2) + w * w;
                emit(r2.sqrt()) * h * w / (PI * r2 * r2)
            },
            tol,
        )
    };
    top + side(x, y) + side(1.0 - x, y) + side(y, x) + side(1.0 - y, x)
}

/// Mean of [`cube_wall_flux`] over the rectangle `x0 < x < x1`,
/// `y0 < y < y1` of the wall: the 6 × 6 Gauss rule on `n × n` parts. At the
/// edges of the cube the flux goes as `x ln x`, so the error of a part that
/// touches an edge falls as the square of its size: on the whole wall of
/// `kl = 1`, about 1e-5 with n = 8.
pub fn cube_wall_flux_mean(kl: f64, x: [f64; 2], y: [f64; 2], n: usize) -> f64 {
    let (dx, dy) = ((x[1] - x[0]) / n as f64, (y[1] - y[0]) / n as f64);
    let mut s = 0.0;
    for i in 0..n {
        for j in 0..n {
            let (a, b) = (x[0] + i as f64 * dx, y[0] + j as f64 * dy);
            for &(u, wu) in &GL6 {
                for &(v, wv) in &GL6 {
                    s += wu * wv * cube_wall_flux(kl, a + u * dx, b + v * dy);
                }
            }
        }
    }
    s / (n * n) as f64
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
    fn cube_flux_limits() {
        // optically thick: the wall sees black gas everywhere
        assert!((cube_wall_flux(200.0, 0.5, 0.5) - 1.0).abs() < 1e-6);
        assert!((cube_wall_flux(200.0, 0.1, 0.3) - 1.0).abs() < 1e-4);
        // transparent walls: the view factors of the five walls sum to 1
        // (the gas emission → 0, so q/κ ≈ mean chord · cos weight finite)
        let thin = cube_wall_flux(1e-6, 0.5, 0.5) / 1e-6;
        assert!(thin > 0.5 && thin < 2.0, "{thin}");
        // symmetry of the wall
        let a = cube_wall_flux(1.0, 0.2, 0.7);
        for (x, y) in [(0.8, 0.7), (0.2, 0.3), (0.7, 0.2), (0.3, 0.8)] {
            assert!((cube_wall_flux(1.0, x, y) - a).abs() < 1e-11);
        }
        // the flux falls towards the edges
        assert!(cube_wall_flux(1.0, 0.5, 0.5) > cube_wall_flux(1.0, 0.1, 0.5));
        // the mean over the wall converges as the square of the part size
        // (the flux goes as x ln x at the edges)
        let m: Vec<f64> = [2, 4, 8]
            .map(|n| cube_wall_flux_mean(1.0, [0.0, 1.0], [0.0, 1.0], n))
            .to_vec();
        let ratio = (m[0] - m[1]) / (m[1] - m[2]);
        assert!((ratio - 4.0).abs() < 0.2, "{m:?}");
        assert!((m[2] - 0.446_027_5).abs() < 1e-5, "{m:?}");
    }

    #[test]
    fn opposed_squares() {
        // unit squares at distance 1: F = 0.19982 (catalogue value)
        assert!((parallel_rectangles(1.0, 1.0, 1.0) - 0.199_825).abs() < 1e-5);
        // close plates see each other entirely
        assert!(parallel_rectangles(1.0, 1.0, 1e-4) > 0.999);
    }
}
