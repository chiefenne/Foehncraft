//! Multi-threaded reductions with a result independent of the thread count.
//!
//! The data are cut into blocks of [`foehncraft_core::repro::BLOCK`] values,
//! each block is summed sequentially (in parallel over blocks), and the block
//! sums are combined in the fixed tree order of
//! [`foehncraft_core::repro::pairwise_sum`]. The result is therefore bitwise
//! equal to the sequential `pairwise_sum` for any number of threads (D31).

use foehncraft_core::Real;
use foehncraft_core::exact::{ExactSum, max_total, min_total};
use foehncraft_core::repro::{BLOCK, block_sum, combine_block_sums};
use rayon::prelude::*;

/// Values per task of the exact reductions: large enough that the merge of
/// two accumulators (68 limbs) is negligible against the additions.
pub const EXACT_CHUNK: usize = 4096;

/// Exact accumulator of `x` (D31 (1)): the same integer for any number of
/// threads, because the chunk accumulators merge exactly.
pub fn par_exact_sum_acc(x: &[f64]) -> ExactSum {
    x.par_chunks(EXACT_CHUNK)
        .map(|c| {
            let mut s = ExactSum::new();
            s.add_slice(c);
            s
        })
        .reduce(ExactSum::new, |mut a, b| {
            a.merge(&b);
            a
        })
}

/// Exact accumulator of the rounded products `a[i] * b[i]`.
pub fn par_exact_dot_acc(a: &[f64], b: &[f64]) -> ExactSum {
    assert_eq!(a.len(), b.len(), "dot product of unequal lengths");
    a.par_chunks(EXACT_CHUNK)
        .zip(b.par_chunks(EXACT_CHUNK))
        .map(|(x, y)| {
            let mut s = ExactSum::new();
            s.add_products(x, y);
            s
        })
        .reduce(ExactSum::new, |mut a, b| {
            a.merge(&b);
            a
        })
}

/// Exact sum of `x`, rounded once: bitwise the same for any thread count,
/// any order of the values and any split between ranks.
pub fn par_exact_sum(x: &[f64]) -> f64 {
    par_exact_sum_acc(x).value()
}

/// Exact dot product (sum of the rounded products), rounded once.
pub fn par_exact_dot(a: &[f64], b: &[f64]) -> f64 {
    par_exact_dot_acc(a, b).value()
}

/// Minimum and maximum of `x` with [`min_total`] and [`max_total`]:
/// independent of order and grouping, a NaN gives NaN. The empty slice
/// gives (+∞, −∞).
pub fn par_min_max(x: &[f64]) -> (f64, f64) {
    x.par_chunks(EXACT_CHUNK)
        .map(|c| {
            c.iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
                    (min_total(a, v), max_total(b, v))
                })
        })
        .reduce(
            || (f64::INFINITY, f64::NEG_INFINITY),
            |(a, b), (c, d)| (min_total(a, c), max_total(b, d)),
        )
}

/// Sum of `x`, bitwise equal to `foehncraft_core::repro::pairwise_sum(x)` for
/// any number of threads.
pub fn par_sum<T: Real>(x: &[T]) -> T {
    let sums: Vec<T> = x.par_chunks(BLOCK).map(block_sum).collect();
    combine_block_sums(&sums)
}

/// Dot product of `a` and `b` in the same fixed order as [`par_sum`].
pub fn par_dot<T: Real>(a: &[T], b: &[T]) -> T {
    assert_eq!(a.len(), b.len());
    let sums: Vec<T> = a
        .par_chunks(BLOCK)
        .zip(b.par_chunks(BLOCK))
        .map(|(x, y)| {
            let mut s = T::ZERO;
            for (&u, &v) in x.iter().zip(y) {
                s += u * v;
            }
            s
        })
        .collect();
    combine_block_sums(&sums)
}

#[cfg(test)]
mod tests {
    use super::*;
    use foehncraft_core::repro::pairwise_sum;

    fn data(n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| foehncraft_core::rng::uniform([i as u64, 1, 2]) - 0.3)
            .collect()
    }

    #[test]
    fn independent_of_thread_count() {
        let v = data(1_234_567);
        let reference = pairwise_sum(&v).to_bits();
        for threads in [1, 2, 3, 4, 7] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let s = pool.install(|| par_sum(&v));
            assert_eq!(s.to_bits(), reference, "{threads} threads");
        }
    }

    #[test]
    fn exact_reductions_independent_of_thread_count() {
        let v = data(1_234_567);
        let w: Vec<f64> = v.iter().map(|x| x * 1e6 - 3.0).collect();
        let s = foehncraft_core::exact::exact_sum(&v).to_bits();
        let d = foehncraft_core::exact::exact_dot(&v, &w).to_bits();
        for threads in [1, 2, 3, 4, 7] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                assert_eq!(par_exact_sum(&v).to_bits(), s, "{threads} threads");
                assert_eq!(par_exact_dot(&v, &w).to_bits(), d, "{threads} threads");
                let (mn, mx) = par_min_max(&w);
                assert_eq!(mn, w.iter().copied().fold(f64::INFINITY, f64::min));
                assert_eq!(mx, w.iter().copied().fold(f64::NEG_INFINITY, f64::max));
            });
        }
    }

    #[test]
    fn dot_matches_sum_of_products() {
        let a = data(100_001);
        let b: Vec<f64> = data(100_001).iter().map(|x| x * 0.5 + 1.0).collect();
        let prod: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x * y).collect();
        assert_eq!(par_dot(&a, &b).to_bits(), pairwise_sum(&prod).to_bits());
    }
}
