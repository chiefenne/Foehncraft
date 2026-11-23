//! Exact, order-independent sums of binary64 values (D31 (1), D93 (4)).
//!
//! [`ExactSum`] is a long accumulator in the sense of Kulisch (1971) and
//! Collange, Defour, Graillat & Iakymchuk (2015): the whole range of
//! binary64, from 2⁻¹⁰⁷⁴ (the smallest subnormal) up past 2¹⁰²⁴, held as a
//! fixed-point integer in 68 limbs of 32 bits. Every finite value is added
//! without rounding, so the accumulator after any set of additions is the
//! same integer whatever the order of the terms and however they were split
//! between threads or ranks; only the final conversion to binary64 rounds,
//! once, to nearest with ties to even. A sum, a dot product or a norm is
//! therefore bitwise the same on one thread, many threads and any number of
//! ranks.
//!
//! A dot product adds the rounded products `a[i] * b[i]`: each product is a
//! fixed function of its two factors (no fused multiply-add, D93 (6)), and
//! their sum is exact, so the result does not depend on the order either.
//!
//! Non-finite terms are counted apart: a NaN, or +∞ together with −∞, gives
//! NaN; otherwise an infinite term gives that infinity. A finite exact sum
//! beyond the binary64 range rounds to ±∞.
//!
//! Cost (spike of 2026-11-22, `foehncraft-parallel/examples/exact_sum_spike.rs`):
//! about 5 ns per value against 0.8 ns for the pairwise sum, data in cache;
//! one reduction across ranks sends [`WORDS`] integers.

/// Number of 32-bit limbs: 68 × 32 = 2176 bits ≥ 2098 bits of the binary64
/// range plus headroom for 2³¹ terms of the largest magnitude.
pub const LIMBS: usize = 68;

/// Length of the integer image of an accumulator ([`ExactSum::to_words`]):
/// the limbs and three counters (NaN, +∞, −∞).
pub const WORDS: usize = LIMBS + 3;

/// Additions (or merged additions) after which the carries are propagated:
/// each limb gains less than 2³² per addition, so 2³⁰ additions keep every
/// limb, also after the merge of two accumulators, below 2⁶³.
const PENDING_LIMIT: u64 = 1 << 30;

const MASK: i64 = 0xffff_ffff;

/// Exact accumulator of binary64 values. See the [module](self) text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExactSum {
    limb: [i64; LIMBS],
    pending: u64,
    nan: u64,
    pos_inf: u64,
    neg_inf: u64,
}

impl Default for ExactSum {
    fn default() -> Self {
        Self::new()
    }
}

impl ExactSum {
    /// The empty sum (+0).
    pub const fn new() -> Self {
        Self {
            limb: [0; LIMBS],
            pending: 0,
            nan: 0,
            pos_inf: 0,
            neg_inf: 0,
        }
    }

    /// Adds `x` without rounding.
    #[inline]
    pub fn add(&mut self, x: f64) {
        let b = x.to_bits();
        let e = ((b >> 52) & 0x7ff) as usize;
        let frac = b & ((1u64 << 52) - 1);
        if e == 0x7ff {
            if frac != 0 {
                self.nan += 1;
            } else if b >> 63 == 0 {
                self.pos_inf += 1;
            } else {
                self.neg_inf += 1;
            }
            return;
        }
        // x = m · 2^(pos − 1074) with the lowest accumulator bit = 2^−1074
        let (m, pos) = if e == 0 {
            (frac, 0)
        } else {
            (frac | (1u64 << 52), e - 1)
        };
        if m == 0 {
            return;
        }
        let (l, s) = (pos / 32, pos % 32);
        let v = (m as u128) << s; // at most 53 + 31 bits: three limbs
        let (v0, v1, v2) = (
            (v & 0xffff_ffff) as i64,
            ((v >> 32) & 0xffff_ffff) as i64,
            (v >> 64) as i64,
        );
        if b >> 63 == 0 {
            self.limb[l] += v0;
            self.limb[l + 1] += v1;
            self.limb[l + 2] += v2;
        } else {
            self.limb[l] -= v0;
            self.limb[l + 1] -= v1;
            self.limb[l + 2] -= v2;
        }
        self.pending += 1;
        if self.pending >= PENDING_LIMIT {
            self.normalise();
        }
    }

    /// Adds every value of `x` (batched: one integer bin per binary exponent).
    pub fn add_slice(&mut self, x: &[f64]) {
        if x.len() < SHORT {
            for &v in x {
                self.add(v);
            }
            return;
        }
        let mut b = Bins::new();
        for &v in x {
            b.push(self, v);
        }
        b.flush(self);
    }

    /// Adds the rounded products `a[i] * b[i]`.
    ///
    /// # Panics
    /// When the lengths differ.
    pub fn add_products(&mut self, a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len(), "dot product of unequal lengths");
        if a.len() < SHORT {
            for (&u, &v) in a.iter().zip(b) {
                self.add(u * v);
            }
            return;
        }
        let mut bins = Bins::new();
        for (&u, &v) in a.iter().zip(b) {
            bins.push(self, u * v);
        }
        bins.flush(self);
    }

    /// Adds the signed integer `m` times 2^(pos − 1074), |m| < 2⁶³.
    #[inline]
    fn add_scaled(&mut self, m: i64, pos: usize) {
        if m == 0 {
            return;
        }
        let (l, s) = (pos / 32, pos % 32);
        let v = (m.unsigned_abs() as u128) << s; // at most 63 + 31 bits
        let (v0, v1, v2) = (
            (v & 0xffff_ffff) as i64,
            ((v >> 32) & 0xffff_ffff) as i64,
            (v >> 64) as i64,
        );
        if m > 0 {
            self.limb[l] += v0;
            self.limb[l + 1] += v1;
            self.limb[l + 2] += v2;
        } else {
            self.limb[l] -= v0;
            self.limb[l + 1] -= v1;
            self.limb[l + 2] -= v2;
        }
        self.pending += 1;
        if self.pending >= PENDING_LIMIT {
            self.normalise();
        }
    }

    /// Adds the accumulated terms of `o` (exact, so merging in any order and
    /// any grouping gives the same accumulator).
    pub fn merge(&mut self, o: &ExactSum) {
        for (a, b) in self.limb.iter_mut().zip(&o.limb) {
            *a += b;
        }
        self.pending += o.pending;
        self.nan += o.nan;
        self.pos_inf += o.pos_inf;
        self.neg_inf += o.neg_inf;
        if self.pending >= PENDING_LIMIT {
            self.normalise();
        }
    }

    /// Propagates the carries: limbs `0..LIMBS-1` in `0..2³²`, the last
    /// limb signed. Does not change the value.
    fn normalise(&mut self) {
        let mut carry = 0i64;
        for l in &mut self.limb[..LIMBS - 1] {
            let t = *l + carry;
            *l = t & MASK;
            carry = t >> 32; // arithmetic shift: floor division by 2³²
        }
        self.limb[LIMBS - 1] += carry;
        self.pending = 0;
    }

    /// The integer image for a reduction across ranks: normalised limbs and
    /// the three counters. The element-wise integer sum of the images of any
    /// set of accumulators is the image of their merge (for fewer than 2³⁰
    /// accumulators), so an `MPI_Allreduce` with `MPI_SUM` on `i64` is exact
    /// and independent of the reduction tree.
    pub fn to_words(&self) -> [i64; WORDS] {
        let mut c = self.clone();
        c.normalise();
        let mut w = [0i64; WORDS];
        w[..LIMBS].copy_from_slice(&c.limb);
        w[LIMBS] = c.nan as i64;
        w[LIMBS + 1] = c.pos_inf as i64;
        w[LIMBS + 2] = c.neg_inf as i64;
        w
    }

    /// The accumulator of an integer image (see [`Self::to_words`]).
    pub fn from_words(w: &[i64; WORDS]) -> Self {
        let mut limb = [0i64; LIMBS];
        limb.copy_from_slice(&w[..LIMBS]);
        let mut s = Self {
            limb,
            pending: 0,
            nan: w[LIMBS] as u64,
            pos_inf: w[LIMBS + 1] as u64,
            neg_inf: w[LIMBS + 2] as u64,
        };
        s.normalise();
        s
    }

    /// The sum rounded once to binary64, to nearest with ties to even.
    /// An exact zero is +0.
    pub fn value(&self) -> f64 {
        if self.nan > 0 || (self.pos_inf > 0 && self.neg_inf > 0) {
            return f64::NAN;
        }
        if self.pos_inf > 0 {
            return f64::INFINITY;
        }
        if self.neg_inf > 0 {
            return f64::NEG_INFINITY;
        }
        let mut c = self.clone();
        c.normalise();
        let neg = c.limb[LIMBS - 1] < 0;
        let mut mag = c.limb;
        if neg {
            // two's complement negation of the multi-limb integer
            let mut carry = 1i64;
            for l in &mut mag[..LIMBS - 1] {
                let t = (!*l & MASK) + carry;
                *l = t & MASK;
                carry = t >> 32;
            }
            mag[LIMBS - 1] = -mag[LIMBS - 1] - 1 + carry;
        }
        let bits = round_magnitude(&mag);
        f64::from_bits(bits | if neg { 1u64 << 63 } else { 0 })
    }
}

/// Below this length a slice is added value by value: setting up the
/// 16 kB of bins would cost more than it saves.
const SHORT: usize = 256;

/// Values added to one bin before the bins are flushed: 2⁹ significands
/// below 2⁵³ keep a bin below 2⁶².
const BIN_BATCH: u32 = 512;

/// The fast path of [`ExactSum::add_slice`] and
/// [`ExactSum::add_products`]: one signed 64-bit integer per binary
/// exponent collects the significands of that exponent without rounding
/// (one integer addition per value), and the bins in the range touched are
/// flushed into the long accumulator every [`BIN_BATCH`] values. Exact, so
/// it gives the accumulator of [`ExactSum::add`] value by value.
struct Bins {
    bin: [i64; 2047],
    count: u32,
    lo: usize,
    hi: usize,
}

impl Bins {
    fn new() -> Self {
        Self {
            bin: [0; 2047],
            count: 0,
            lo: usize::MAX,
            hi: 0,
        }
    }

    #[inline(always)]
    fn push(&mut self, acc: &mut ExactSum, x: f64) {
        let b = x.to_bits();
        let e = ((b >> 52) & 0x7ff) as usize;
        if e == 0x7ff {
            acc.add(x); // NaN and infinities: counted by the slow path
            return;
        }
        let frac = (b & ((1u64 << 52) - 1)) as i64;
        // subnormals (e = 0) share the scale of e = 1, without the hidden bit
        let (m, k) = if e == 0 {
            (frac, 1)
        } else {
            (frac | (1i64 << 52), e)
        };
        let m = if b >> 63 == 0 { m } else { -m };
        self.bin[k] += m;
        self.lo = self.lo.min(k);
        self.hi = self.hi.max(k);
        self.count += 1;
        if self.count == BIN_BATCH {
            self.flush(acc);
        }
    }

    fn flush(&mut self, acc: &mut ExactSum) {
        if self.lo <= self.hi {
            for k in self.lo..=self.hi {
                let m = core::mem::take(&mut self.bin[k]);
                acc.add_scaled(m, k - 1);
            }
        }
        self.count = 0;
        self.lo = usize::MAX;
        self.hi = 0;
    }
}

/// Bit `i` of the non-negative multi-limb integer `m`.
#[inline]
fn bit(m: &[i64; LIMBS], i: usize) -> bool {
    (m[i / 32] >> (i % 32)) & 1 == 1
}

/// The binary64 bit pattern (sign clear) of the non-negative integer `m`
/// times 2⁻¹⁰⁷⁴, rounded to nearest with ties to even.
fn round_magnitude(m: &[i64; LIMBS]) -> u64 {
    // the top limb may exceed 32 bits only for sums far beyond the range
    if m[LIMBS - 1] >= 1 << 32 {
        return f64::INFINITY.to_bits();
    }
    let Some(top) = (0..LIMBS).rev().find(|&i| m[i] != 0) else {
        return 0;
    };
    let p = top * 32 + (63 - (m[top] as u64).leading_zeros() as usize); // index of the leading bit
    if p <= 52 {
        // fewer than 54 bits: exact; the integer is the bit pattern (the
        // subnormals and the first binade share the ulp 2⁻¹⁰⁷⁴)
        return (m[0] as u64) | ((m[1] as u64) << 32);
    }
    let shift = p - 52; // lowest bit of the 53-bit significand
    let mut mant: u64 = 0;
    for i in (shift..=p).rev() {
        mant = (mant << 1) | bit(m, i) as u64;
    }
    let round = bit(m, shift - 1);
    let sticky = round && sticky_below(m, shift - 1);
    let mut shift = shift as u64;
    if round && (sticky || mant & 1 == 1) {
        mant += 1;
        if mant == 1 << 53 {
            mant >>= 1;
            shift += 1;
        }
    }
    // value = mant · 2^(shift − 1074), mant in [2⁵², 2⁵³): biased exponent
    // shift + 1, so the pattern is (shift + 1) << 52 | (mant − 2⁵²)
    if shift + 1 >= 0x7ff {
        return f64::INFINITY.to_bits();
    }
    (shift << 52) + mant
}

/// Whether any bit below index `k` of `m` is set.
fn sticky_below(m: &[i64; LIMBS], k: usize) -> bool {
    let (l, s) = (k / 32, k % 32);
    m[..l].iter().any(|&x| x != 0) || (m[l] & ((1i64 << s) - 1)) != 0
}

/// Exact sum of `x`, rounded once (see [`ExactSum`]).
pub fn exact_sum(x: &[f64]) -> f64 {
    let mut s = ExactSum::new();
    s.add_slice(x);
    s.value()
}

/// Exact sum of the rounded products `a[i] * b[i]`, rounded once.
pub fn exact_dot(a: &[f64], b: &[f64]) -> f64 {
    let mut s = ExactSum::new();
    s.add_products(a, b);
    s.value()
}

/// `.exact_sum()` on an iterator of `f64` (or `&f64`): the exact sum of
/// [`ExactSum`], rounded once, in place of `.sum::<f64>()` wherever a sum
/// over cells, faces or boundary patches feeds the solution or a written
/// number (D93 (4)).
pub trait ExactSumIter {
    /// Exact sum of the items, rounded once.
    fn exact_sum(self) -> f64;
}

impl<I> ExactSumIter for I
where
    I: Iterator,
    I::Item: core::borrow::Borrow<f64>,
{
    fn exact_sum(self) -> f64 {
        let mut s = ExactSum::new();
        for v in self {
            s.add(*core::borrow::Borrow::<f64>::borrow(&v));
        }
        s.value()
    }
}

/// Smaller of two values, commutative and associative on every input, so a
/// minimum over any order or grouping has the same bits: a NaN operand
/// gives NaN, and −0 is smaller than +0.
#[inline]
pub fn min_total(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a < b || (a == b && a.is_sign_negative()) {
        a
    } else {
        b
    }
}

/// Larger of two values, the counterpart of [`min_total`]: a NaN operand
/// gives NaN, and +0 is larger than −0.
#[inline]
pub fn max_total(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a > b || (a == b && a.is_sign_positive()) {
        a
    } else {
        b
    }
}

#[cfg(test)]
mod tests {
    extern crate std;
    use super::*;
    use crate::repro::pairwise_sum;
    use crate::rng::uniform;
    use std::vec::Vec;

    fn wide(n: usize, seed: u64) -> Vec<f64> {
        (0..n)
            .map(|i| {
                let u = uniform([i as u64, seed, 3]) - 0.5;
                let k = (i % 23) as i32 - 11;
                u * 10f64.powi(k)
            })
            .collect()
    }

    #[test]
    fn simple_values_are_exact() {
        assert_eq!(exact_sum(&[]).to_bits(), 0.0f64.to_bits());
        assert_eq!(exact_sum(&[1.0, 2.0, 3.5]), 6.5);
        assert_eq!(exact_sum(&[1e300, 1.0, -1e300]), 1.0);
        assert_eq!(exact_sum(&[-2.5]), -2.5);
        assert_eq!(exact_sum(&[-0.0]).to_bits(), 0.0f64.to_bits());
        assert_eq!(exact_sum(&[1.0, -1.0]).to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn subnormals_and_extremes() {
        let tiny = f64::from_bits(1);
        assert_eq!(exact_sum(&[tiny, tiny, tiny]), f64::from_bits(3));
        assert_eq!(
            exact_sum(&[f64::MIN_POSITIVE, -tiny]).to_bits(),
            f64::MIN_POSITIVE.to_bits() - 1
        );
        assert_eq!(exact_sum(&[f64::MAX, -f64::MAX, f64::MAX]), f64::MAX);
        assert_eq!(exact_sum(&[f64::MAX, f64::MAX]), f64::INFINITY);
        assert_eq!(exact_sum(&[-f64::MAX, -f64::MAX]), f64::NEG_INFINITY);
        assert_eq!(exact_sum(&[f64::MAX, f64::MAX, -f64::MAX]), f64::MAX);
    }

    #[test]
    fn non_finite_terms() {
        assert!(exact_sum(&[1.0, f64::NAN]).is_nan());
        assert!(exact_sum(&[f64::INFINITY, f64::NEG_INFINITY]).is_nan());
        assert_eq!(exact_sum(&[f64::INFINITY, -1e308]), f64::INFINITY);
        assert_eq!(exact_sum(&[f64::NEG_INFINITY, 1.0]), f64::NEG_INFINITY);
    }

    #[test]
    fn rounding_is_to_nearest_even() {
        // 1 + 2^-53 is a tie: rounds to 1 (even); 1 + 2^-53 + 2^-80 rounds up
        let h = 2f64.powi(-53);
        assert_eq!(exact_sum(&[1.0, h]), 1.0);
        assert_eq!(exact_sum(&[1.0, h, 2f64.powi(-80)]), 1.0 + 2.0 * h);
        // (1 + 2^-52) + 2^-53 is a tie with an odd significand: rounds up
        let one_ulp = 1.0 + 2.0 * h;
        assert_eq!(exact_sum(&[one_ulp, h]), 1.0 + 4.0 * h);
        assert_eq!(exact_sum(&[-1.0, -h]), -1.0);
        // carry out of the significand
        let below_two = 2.0 - 2.0 * h;
        assert_eq!(exact_sum(&[below_two, h]), 2.0);
    }

    #[test]
    fn agrees_with_a_compensated_sum_where_that_is_exact() {
        // two-sum chains of a few terms are exact: compare the rounding
        for i in 0..2000u64 {
            let a = (uniform([i, 1, 1]) - 0.5) * 1e8;
            let b = (uniform([i, 2, 1]) - 0.5) * 1e-3;
            let s = a + b;
            assert_eq!(exact_sum(&[a, b]), s, "case {i}");
        }
    }

    #[test]
    fn correctly_rounded_against_an_integer_reference() {
        // values k · 2^-60 with |k| < 2^62: their exact sum is an i128, and
        // `i128 as f64` rounds to nearest even (Rust reference), then the
        // scaling by 2^-60 is exact
        let scale = 2f64.powi(-60);
        for case in 0..400u64 {
            let n = 2 + (case % 37) as usize;
            let ks: Vec<i64> = (0..n)
                .map(|i| {
                    let u = uniform([case, i as u64, 5]);
                    let bits = 10 + (case as usize * 7 + i * 13) % 52;
                    ((u - 0.5) * 2f64.powi(bits as i32)) as i64
                })
                .collect();
            let xs: Vec<f64> = ks.iter().map(|&k| k as f64 * scale).collect();
            // the f64 values are exact multiples of 2^-60 (|k| < 2^53 after
            // rounding to f64)
            let total: i128 = xs.iter().map(|&x| (x / scale) as i128).sum();
            let want = total as f64 * scale;
            assert_eq!(exact_sum(&xs).to_bits(), want.to_bits(), "case {case}");
        }
    }

    #[test]
    fn independent_of_order_and_grouping() {
        let v = wide(100_003, 7);
        let reference = exact_sum(&v).to_bits();
        let mut rev = v.clone();
        rev.reverse();
        assert_eq!(exact_sum(&rev).to_bits(), reference);
        for parts in [2usize, 3, 4, 7, 16] {
            let cuts: Vec<usize> = (0..=parts).map(|k| k * v.len() / parts).collect();
            let mut total = ExactSum::new();
            for k in (0..parts).rev() {
                let mut local = ExactSum::new();
                local.add_slice(&v[cuts[k]..cuts[k + 1]]);
                // through the integer image, as across ranks
                let w = local.to_words();
                total.merge(&ExactSum::from_words(&w));
            }
            assert_eq!(total.value().to_bits(), reference, "{parts} parts");
        }
        // the word-wise integer sum of the images is the image of the merge
        let (a, b) = (&v[..500], &v[500..]);
        let (mut ea, mut eb) = (ExactSum::new(), ExactSum::new());
        ea.add_slice(a);
        eb.add_slice(b);
        let (wa, wb) = (ea.to_words(), eb.to_words());
        let mut w = [0i64; WORDS];
        for i in 0..WORDS {
            w[i] = wa[i] + wb[i];
        }
        assert_eq!(ExactSum::from_words(&w).value().to_bits(), reference);
    }

    #[test]
    fn more_accurate_than_pairwise() {
        let v = wide(200_000, 11);
        let e = exact_sum(&v);
        let p = pairwise_sum(&v);
        let scale: f64 = v.iter().map(|x| x.abs()).sum();
        assert!((e - p).abs() <= 1e-14 * scale);
    }

    #[test]
    fn dot_is_order_independent() {
        let a = wide(10_001, 1);
        let b = wide(10_001, 2);
        let d = exact_dot(&a, &b);
        let (mut ra, mut rb) = (a.clone(), b.clone());
        ra.reverse();
        rb.reverse();
        assert_eq!(exact_dot(&ra, &rb).to_bits(), d.to_bits());
        let prods: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x * y).collect();
        assert_eq!(exact_sum(&prods).to_bits(), d.to_bits());
    }

    #[test]
    fn many_terms_force_normalisation() {
        let mut s = ExactSum::new();
        s.pending = PENDING_LIMIT - 3;
        for _ in 0..10 {
            s.add(f64::MAX);
            s.add(-f64::MAX);
        }
        s.add(0.5);
        assert_eq!(s.value(), 0.5);
    }

    #[test]
    fn batched_path_equals_value_by_value() {
        let mut v = wide(5000, 3);
        v.extend([
            f64::from_bits(1),
            -f64::from_bits(7),
            f64::MAX,
            -f64::MAX,
            0.0,
            -0.0,
        ]);
        v.extend((0..3000).map(|i| 2f64.powi(-1060 + (i % 40))));
        let mut slow = ExactSum::new();
        for &x in &v {
            slow.add(x);
        }
        let mut fast = ExactSum::new();
        fast.add_slice(&v);
        assert_eq!(fast.to_words(), slow.to_words());
        assert_eq!(fast.value().to_bits(), slow.value().to_bits());
        let mut nan = ExactSum::new();
        nan.add_slice(&[1.0, f64::NAN, 2.0]);
        assert!(nan.value().is_nan());
    }

    #[test]
    fn total_min_and_max() {
        assert_eq!(min_total(0.0, -0.0).to_bits(), (-0.0f64).to_bits());
        assert_eq!(min_total(-0.0, 0.0).to_bits(), (-0.0f64).to_bits());
        assert_eq!(max_total(-0.0, 0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(max_total(0.0, -0.0).to_bits(), 0.0f64.to_bits());
        assert!(min_total(f64::NAN, 1.0).is_nan());
        assert!(max_total(1.0, f64::NAN).is_nan());
        assert_eq!(min_total(1.0, 2.0), 1.0);
        assert_eq!(max_total(1.0, 2.0), 2.0);
    }
}
