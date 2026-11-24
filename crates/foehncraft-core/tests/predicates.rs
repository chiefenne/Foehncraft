//! Exact orientation predicates (D225) against an independent oracle: the
//! same determinants in 128-bit integer arithmetic on points whose
//! coordinates are integers times a power of two (scaling by a power of two
//! is exact in binary64 and does not change a sign).

use foehncraft_core::predicates::{orient2d, orient2d_exact, orient3d, orient3d_exact};
use foehncraft_core::rng;

fn sign_f(x: f64) -> i32 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

fn sign_i(x: i128) -> i32 {
    x.signum() as i32
}

fn det3_i(a: [i64; 3], b: [i64; 3], c: [i64; 3], d: [i64; 3]) -> i128 {
    let s = |p: [i64; 3]| [0, 1, 2].map(|k| (p[k] - d[k]) as i128);
    let (a, b, c) = (s(a), s(b), s(c));
    a[2] * (b[0] * c[1] - c[0] * b[1])
        + b[2] * (c[0] * a[1] - a[0] * c[1])
        + c[2] * (a[0] * b[1] - b[0] * a[1])
}

fn det2_i(a: [i64; 2], b: [i64; 2], c: [i64; 2]) -> i128 {
    (a[0] - c[0]) as i128 * (b[1] - c[1]) as i128 - (a[1] - c[1]) as i128 * (b[0] - c[0]) as i128
}

/// Uniform integer in [-r, r].
fn int(key: [u64; 4], r: i64) -> i64 {
    (rng::keyed(key) % (2 * r as u64 + 1)) as i64 - r
}

fn to_f3(p: [i64; 3], scale: f64) -> [f64; 3] {
    p.map(|v| v as f64 * scale)
}

#[test]
fn orient3d_random_and_nearly_coplanar() {
    let mut naive_wrong = 0;
    let mut checked = 0;
    for i in 0..20_000u64 {
        // exactly coplanar integer points, then one coordinate moved by
        // -1, 0 or +1: the determinant is tiny or zero
        let a = [0, 1, 2].map(|k| int([1, i, k, 0], 1 << 18));
        let b = [0, 1, 2].map(|k| int([1, i, k, 1], 1 << 18));
        let c = [0, 1, 2].map(|k| int([1, i, k, 2], 1 << 18));
        let (s, t) = (int([2, i, 0, 0], 1 << 10), int([2, i, 1, 0], 1 << 10));
        let mut d = [0, 1, 2].map(|k| a[k] + s * (b[k] - a[k]) + t * (c[k] - a[k]));
        d[(i % 3) as usize] += int([3, i, 0, 0], 1);
        let want = sign_i(det3_i(a, b, c, d));
        // power-of-two scales, exact
        for scale in [1.0, 2f64.powi(-70), 2f64.powi(40)] {
            let (fa, fb, fc, fd) = (
                to_f3(a, scale),
                to_f3(b, scale),
                to_f3(c, scale),
                to_f3(d, scale),
            );
            assert_eq!(
                sign_f(orient3d(fa, fb, fc, fd)),
                want,
                "case {i} scale {scale}"
            );
            assert_eq!(
                sign_f(orient3d_exact(fa, fb, fc, fd)),
                want,
                "exact case {i}"
            );
            // every permutation changes the sign by its parity
            assert_eq!(sign_f(orient3d(fb, fa, fc, fd)), -want);
            assert_eq!(sign_f(orient3d(fa, fb, fd, fc)), -want);
            checked += 1;
            let naive = {
                let s = |p: [f64; 3]| [0, 1, 2].map(|k| p[k] - fd[k]);
                let (x, y, z) = (s(fa), s(fb), s(fc));
                x[2] * (y[0] * z[1] - z[0] * y[1])
                    + y[2] * (z[0] * x[1] - x[0] * z[1])
                    + z[2] * (x[0] * y[1] - y[0] * x[1])
            };
            if sign_f(naive) != want {
                naive_wrong += 1;
            }
        }
    }
    // random points far from degenerate
    for i in 0..20_000u64 {
        let p = |j: u64| [0, 1, 2].map(|k| int([4, i, k, j], 1 << 40));
        let (a, b, c, d) = (p(0), p(1), p(2), p(3));
        let want = sign_i(det3_i(a, b, c, d));
        let s = 2f64.powi(-20);
        assert_eq!(
            sign_f(orient3d(to_f3(a, s), to_f3(b, s), to_f3(c, s), to_f3(d, s))),
            want
        );
    }
    eprintln!(
        "orient3d: {checked} nearly coplanar cases, plain floating point wrong in {naive_wrong}"
    );
    assert!(naive_wrong > 0, "the test cases are not hard enough");
}

/// The classic near-collinear grid (points 0.5 + i ulp against a line
/// through (12, 12) and (24, 24)), where plain floating point gives a
/// scattered, wrong pattern of signs.
#[test]
fn orient2d_near_collinear_grid() {
    let u = 2f64.powi(-53);
    let mut naive_wrong = 0;
    for i in 0..256i64 {
        for j in 0..256i64 {
            let p = [0.5 + i as f64 * u, 0.5 + j as f64 * u];
            let (q, r) = ([12.0, 12.0], [24.0, 24.0]);
            // the same points as integers times 2⁻⁵³
            let pi = [(1i64 << 52) + i, (1i64 << 52) + j];
            let qi = [12i64 << 53, 12i64 << 53];
            let ri = [24i64 << 53, 24i64 << 53];
            let want = sign_i(det2_i(pi, qi, ri));
            assert_eq!(sign_f(orient2d(p, q, r)), want, "{i} {j}");
            assert_eq!(sign_f(orient2d_exact(p, q, r)), want, "exact {i} {j}");
            let naive = (p[0] - r[0]) * (q[1] - r[1]) - (p[1] - r[1]) * (q[0] - r[0]);
            if sign_f(naive) != want {
                naive_wrong += 1;
            }
        }
    }
    eprintln!("orient2d: plain floating point wrong in {naive_wrong} of 65536 points");
    assert!(naive_wrong > 0);
}

#[test]
fn orient2d_random_integers() {
    for i in 0..50_000u64 {
        let p = |j: u64| [0, 1].map(|k| int([5, i, k, j], 1 << 50));
        let (a, b, c) = (p(0), p(1), p(2));
        // nearly collinear in half of the cases
        let c = if i % 2 == 0 {
            let t = int([6, i, 0, 0], 8);
            [
                a[0] + t * (b[0] - a[0]) / 8,
                a[1] + t * (b[1] - a[1]) / 8 + int([7, i, 0, 0], 1),
            ]
        } else {
            c
        };
        let ok = c.iter().all(|v| v.abs() < (1 << 52));
        if !ok {
            continue;
        }
        let want = sign_i(det2_i(a, b, c));
        let f = |q: [i64; 2]| q.map(|v| v as f64);
        assert_eq!(sign_f(orient2d(f(a), f(b), f(c))), want, "{i}");
    }
}
