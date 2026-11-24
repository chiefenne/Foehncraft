//! Exact geometric predicates: the signs of the 2D and 3D orientation
//! determinants of binary64 points, always right (D225, D108).
//!
//! Written from the paper of J. R. Shewchuk, *Adaptive precision
//! floating-point arithmetic and fast robust geometric predicates*,
//! Discrete & Computational Geometry 18 (1997) 305–363: the error-free
//! transformations (Two-Sum, Two-Diff, Split, Two-Product, after Knuth,
//! Dekker and Priest), expansions (sums of non-overlapping binary64
//! components in increasing magnitude), and the first-stage error bounds of
//! the orientation tests. No code was copied (D201).
//!
//! Each predicate first evaluates the determinant in floating point and
//! compares it with the paper's a-priori bound of its rounding error
//! (stage A). Only when the sign is uncertain (points coplanar or nearly so)
//! does it evaluate the determinant **exactly** as an expansion, from the
//! exact coordinate differences (Two-Diff) and exact products
//! (Two-Product); the paper's intermediate stages B and C are left out,
//! because the exact path is rare and its cost does not matter. The sign
//! of the result is therefore exact for every finite input whose products
//! neither overflow nor underflow (the paper's assumption).
//!
//! No fused multiply-add is used (D93 (6)): Rust never contracts a product
//! and a sum into one, and Two-Product splits its factors after Dekker, so
//! the result is the same on every machine.

/// Machine epsilon of the paper for round to nearest: 2⁻⁵³.
const EPS: f64 = 1.1102230246251565e-16;
/// 2²⁷ + 1: splits a binary64 into two halves of 26 bits.
const SPLITTER: f64 = 134_217_729.0;
/// Stage-A error bound of orient2d (paper, section 4.3): (3 + 16ε)ε.
const O2D_BOUND: f64 = (3.0 + 16.0 * EPS) * EPS;
/// Stage-A error bound of orient3d: (7 + 56ε)ε.
const O3D_BOUND: f64 = (7.0 + 56.0 * EPS) * EPS;

/// Largest expansion an exact orientation test needs (orient3d: three
/// terms of 64 components).
const CAP: usize = 192;

/// `x + y = a + b` exactly, `x = fl(a + b)` (Knuth's Two-Sum).
#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bv = x - a;
    let av = x - bv;
    let br = b - bv;
    let ar = a - av;
    (x, ar + br)
}

/// `x + y = a - b` exactly, `x = fl(a - b)` (Two-Diff).
#[inline]
fn two_diff(a: f64, b: f64) -> (f64, f64) {
    let x = a - b;
    let bv = a - x;
    let av = x + bv;
    let br = bv - b;
    let ar = a - av;
    (x, ar + br)
}

/// `a = hi + lo` with both halves of at most 26 significant bits (Dekker's
/// split).
#[inline]
fn split(a: f64) -> (f64, f64) {
    let c = SPLITTER * a;
    let big = c - a;
    let hi = c - big;
    (hi, a - hi)
}

/// `x + y = a b` exactly, `x = fl(a b)` (Dekker's Two-Product).
#[inline]
fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;
    let (ah, al) = split(a);
    let (bh, bl) = split(b);
    let e1 = x - ah * bh;
    let e2 = e1 - al * bh;
    let e3 = e2 - ah * bl;
    (x, al * bl - e3)
}

/// An expansion: the exact sum of its components, which do not overlap
/// and grow in magnitude; zero components are dropped.
#[derive(Clone, Copy)]
struct Expansion {
    c: [f64; CAP],
    n: usize,
}

impl Expansion {
    const ZERO: Expansion = Expansion {
        c: [0.0; CAP],
        n: 0,
    };

    /// The exact value of `a - b`.
    fn diff(a: f64, b: f64) -> Self {
        let (x, y) = two_diff(a, b);
        let mut e = Self::ZERO;
        e.push(y);
        e.push(x);
        e
    }

    fn push(&mut self, v: f64) {
        if v != 0.0 {
            assert!(self.n < CAP, "expansion capacity");
            self.c[self.n] = v;
            self.n += 1;
        }
    }

    /// Adds `b` exactly (Grow-Expansion, paper Theorem 10, with the zero
    /// components eliminated).
    fn grow(&mut self, b: f64) {
        let mut q = b;
        let mut out = Self::ZERO;
        for i in 0..self.n {
            let (s, h) = two_sum(q, self.c[i]);
            out.push(h);
            q = s;
        }
        out.push(q);
        *self = out;
    }

    /// `self + other`, exactly (repeated Grow-Expansion).
    fn add(&self, other: &Expansion) -> Expansion {
        let mut r = *self;
        for i in 0..other.n {
            r.grow(other.c[i]);
        }
        r
    }

    /// `self - other`, exactly.
    fn sub(&self, other: &Expansion) -> Expansion {
        let mut r = *self;
        for i in 0..other.n {
            r.grow(-other.c[i]);
        }
        r
    }

    /// `self * b`, exactly: every component times `b` by Two-Product, both
    /// parts added by Grow-Expansion.
    fn scale(&self, b: f64) -> Expansion {
        let mut r = Self::ZERO;
        for i in 0..self.n {
            let (x, y) = two_product(self.c[i], b);
            r.grow(y);
            r.grow(x);
        }
        r
    }

    /// `self * other`, exactly.
    fn mul(&self, other: &Expansion) -> Expansion {
        let mut r = Self::ZERO;
        for j in 0..other.n {
            r = r.add(&self.scale(other.c[j]));
        }
        r
    }

    /// The component of largest magnitude: it has the sign of the exact
    /// value (0 for a zero expansion).
    fn most_significant(&self) -> f64 {
        if self.n == 0 { 0.0 } else { self.c[self.n - 1] }
    }
}

/// Twice the signed area of the triangle (a, b, c): positive when a, b, c
/// turn counter-clockwise, negative when clockwise, zero when collinear.
/// The sign is exact; the magnitude is an approximation.
pub fn orient2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let left = (a[0] - c[0]) * (b[1] - c[1]);
    let right = (a[1] - c[1]) * (b[0] - c[0]);
    let det = left - right;
    let sum = if left > 0.0 {
        if right <= 0.0 {
            return det;
        }
        left + right
    } else if left < 0.0 {
        if right >= 0.0 {
            return det;
        }
        -left - right
    } else {
        return det;
    };
    let bound = O2D_BOUND * sum;
    if det >= bound || -det >= bound {
        return det;
    }
    orient2d_exact(a, b, c)
}

/// [`orient2d`] by exact expansion arithmetic only (for the tests).
pub fn orient2d_exact(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let acx = Expansion::diff(a[0], c[0]);
    let acy = Expansion::diff(a[1], c[1]);
    let bcx = Expansion::diff(b[0], c[0]);
    let bcy = Expansion::diff(b[1], c[1]);
    acx.mul(&bcy).sub(&acy.mul(&bcx)).most_significant()
}

/// Six times the signed volume of the tetrahedron (a, b, c, d), the
/// determinant (a − d) · ((b − d) × (c − d)): positive when `d` lies below
/// the plane of a, b, c, where "below" means a, b, c turn counter-clockwise
/// seen from above; negative above; zero when the four points are coplanar.
/// The sign is exact; the magnitude is an approximation.
pub fn orient3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let (adx, ady, adz) = (a[0] - d[0], a[1] - d[1], a[2] - d[2]);
    let (bdx, bdy, bdz) = (b[0] - d[0], b[1] - d[1], b[2] - d[2]);
    let (cdx, cdy, cdz) = (c[0] - d[0], c[1] - d[1], c[2] - d[2]);
    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * adz.abs()
        + (cdxady.abs() + adxcdy.abs()) * bdz.abs()
        + (adxbdy.abs() + bdxady.abs()) * cdz.abs();
    let bound = O3D_BOUND * permanent;
    if det > bound || -det > bound {
        return det;
    }
    orient3d_exact(a, b, c, d)
}

/// [`orient3d`] by exact expansion arithmetic only (for the tests).
pub fn orient3d_exact(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let ad = [0, 1, 2].map(|k| Expansion::diff(a[k], d[k]));
    let bd = [0, 1, 2].map(|k| Expansion::diff(b[k], d[k]));
    let cd = [0, 1, 2].map(|k| Expansion::diff(c[k], d[k]));
    // the same expansion along the third column as the filter above:
    // adz (bdx cdy − cdx bdy) + bdz (cdx ady − adx cdy) + cdz (adx bdy − bdx ady)
    let m_a = bd[0].mul(&cd[1]).sub(&cd[0].mul(&bd[1]));
    let m_b = cd[0].mul(&ad[1]).sub(&ad[0].mul(&cd[1]));
    let m_c = ad[0].mul(&bd[1]).sub(&bd[0].mul(&ad[1]));
    ad[2]
        .mul(&m_a)
        .add(&bd[2].mul(&m_b))
        .add(&cd[2].mul(&m_c))
        .most_significant()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_free_transformations() {
        let (x, y) = two_sum(1.0, 1e-30);
        assert_eq!((x, y), (1.0, 1e-30));
        let (x, y) = two_product(1.0 + EPS * 2.0, 1.0 + EPS * 2.0);
        // (1 + 2ε)² = 1 + 4ε + 4ε²
        assert_eq!(x, 1.0 + 4.0 * EPS);
        assert_eq!(y, 4.0 * EPS * EPS);
    }

    #[test]
    fn simple_orientations() {
        assert!(orient2d([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]) > 0.0);
        assert!(orient2d([0.0, 0.0], [0.0, 1.0], [1.0, 0.0]) < 0.0);
        assert_eq!(orient2d([0.0, 0.0], [1.0, 1.0], [2.0, 2.0]), 0.0);
        let (a, b, c) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(orient3d(a, b, c, [0.0, 0.0, -1.0]) > 0.0);
        assert!(orient3d(a, b, c, [0.0, 0.0, 1.0]) < 0.0);
        assert_eq!(orient3d(a, b, c, [0.3, 0.7, 0.0]), 0.0);
    }
}
