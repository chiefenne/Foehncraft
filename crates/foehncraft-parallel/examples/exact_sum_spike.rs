//! Sizing spike for the reproducible mode of D93 (4): the cost of an exact,
//! order-independent sum against the fixed-order pairwise sum, and the proof
//! that it does not depend on the partition.
//!
//! ```text
//! cargo run --release -p foehncraft-parallel --example exact_sum_spike
//! ```
//!
//! The accumulator is a long accumulator in the sense of Kulisch (1971) and
//! Collange et al. (2015): the 2098 bit range of binary64 in 32 bit limbs held
//! in 64 bit words, so that carries can wait until the end. Adding a value
//! touches at most three limbs.

use foehncraft_core::repro::pairwise_sum;
use std::time::Instant;

const LIMBS: usize = 68; // 68 * 32 = 2176 bits > 2098 + headroom

#[derive(Clone)]
struct Exact {
    limb: [i64; LIMBS],
    adds: u64,
}

impl Exact {
    fn new() -> Self {
        Self {
            limb: [0; LIMBS],
            adds: 0,
        }
    }

    /// Adds `x` (finite) exactly.
    #[inline]
    fn add(&mut self, x: f64) {
        let b = x.to_bits();
        let e = ((b >> 52) & 0x7ff) as i32;
        let m = b & ((1u64 << 52) - 1);
        let (m, e) = if e == 0 {
            (m, 1)
        } else {
            (m | (1u64 << 52), e)
        };
        if m == 0 {
            return;
        }
        // value = m * 2^(e - 1075); bit position of the lowest bit: e - 1 + 1074 - 1074
        let pos = (e - 1) as usize; // bit index in the accumulator (0 = 2^-1074)
        let (l, s) = (pos / 32, pos % 32);
        let v = (m as u128) << s; // up to 53 + 31 bits
        let sign = if (b >> 63) != 0 { -1i64 } else { 1 };
        self.limb[l] += sign * (v & 0xffff_ffff) as i64;
        self.limb[l + 1] += sign * ((v >> 32) & 0xffff_ffff) as i64;
        self.limb[l + 2] += sign * (v >> 64) as i64;
        self.adds += 1;
        if self.adds & 0x3fff_ffff == 0 {
            self.normalise();
        }
    }

    fn merge(&mut self, o: &Exact) {
        for i in 0..LIMBS {
            self.limb[i] += o.limb[i];
        }
        self.adds += o.adds;
    }

    fn normalise(&mut self) {
        let mut carry = 0i64;
        for l in &mut self.limb {
            let t = *l + carry;
            *l = t & 0xffff_ffff;
            carry = t >> 32; // arithmetic shift: floor
        }
        // the last limb keeps the sign in `carry`
        self.limb[LIMBS - 1] += carry << 32;
        self.adds = 0;
    }

    /// The sum rounded to binary64 (round to nearest through the top 128 bits
    /// plus a sticky bit is not needed here: any fixed function of the exact
    /// value is reproducible; this one is correctly rounded except in
    /// ties-to-even corner cases).
    fn value(&mut self) -> f64 {
        self.normalise();
        // sign from the top limb
        let neg = self.limb[LIMBS - 1] < 0;
        let mut mag = self.limb;
        if neg {
            // two's complement negate the multi-limb number
            let mut carry = 1i64;
            for l in mag.iter_mut() {
                let t = (!(*l) & 0xffff_ffff) + carry;
                *l = t & 0xffff_ffff;
                carry = t >> 32;
            }
        }
        let top = (0..LIMBS).rev().find(|&i| mag[i] != 0);
        let Some(top) = top else { return 0.0 };
        // gather the top 96 bits as f64 pieces, most significant first
        let mut r = 0.0f64;
        let mut scale_bits = (top as i32) * 32 - 1074;
        let mut pieces = 0;
        let mut i = top as i64;
        while i >= 0 && pieces < 4 {
            r += (mag[i as usize] as f64) * 2f64.powi(scale_bits);
            scale_bits -= 32;
            i -= 1;
            pieces += 1;
        }
        if neg { -r } else { r }
    }
}

fn data(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| {
            (foehncraft_core::rng::uniform([i as u64, 7, 9]) - 0.5) * 10f64.powi((i % 9) as i32 - 4)
        })
        .collect()
}

fn main() {
    let n = 20_000_000;
    let v = data(n);
    // partition independence: 1, 2, 3, 4, 7 "ranks" cutting the array at
    // different places and in a shuffled order
    let mut bits = Vec::new();
    for parts in [1usize, 2, 3, 4, 7] {
        let mut total = Exact::new();
        let cuts: Vec<usize> = (0..=parts).map(|k| k * n / parts).collect();
        for k in (0..parts).rev() {
            let mut local = Exact::new();
            for &x in &v[cuts[k]..cuts[k + 1]] {
                local.add(x);
            }
            total.merge(&local);
        }
        bits.push(total.value().to_bits());
    }
    println!("exact sum bits for 1,2,3,4,7 parts: {bits:x?}");
    assert!(bits.iter().all(|&b| b == bits[0]), "partition dependent");
    // the fixed-order pairwise sum follows the partition: per-part sums
    // combined in rank order (what `sum_reproducible` does today)
    let ref_pair = pairwise_sum(&v);
    for parts in [2usize, 3, 4, 7] {
        let cuts: Vec<usize> = (0..=parts).map(|k| k * n / parts).collect();
        let per: Vec<f64> = (0..parts)
            .map(|k| pairwise_sum(&v[cuts[k]..cuts[k + 1]]))
            .collect();
        let c = pairwise_sum(&per);
        println!(
            "pairwise, {parts} parts: {:x} (equal to 1 part: {})",
            c.to_bits(),
            c.to_bits() == ref_pair.to_bits()
        );
    }

    // cost
    let reps = 5;
    let t = Instant::now();
    let mut s = 0.0;
    for _ in 0..reps {
        s += pairwise_sum(std::hint::black_box(&v));
    }
    let t_pair = t.elapsed().as_secs_f64() / reps as f64;
    let t = Instant::now();
    let mut s2 = 0.0;
    for _ in 0..reps {
        let mut e = Exact::new();
        for &x in std::hint::black_box(&v) {
            e.add(x);
        }
        s2 += e.value();
    }
    let t_exact = t.elapsed().as_secs_f64() / reps as f64;
    println!(
        "{n} values: pairwise {:.1} ms ({:.2} ns/value), exact {:.1} ms ({:.2} ns/value), ratio {:.1} ({s:.3e} {s2:.3e})",
        t_pair * 1e3,
        t_pair * 1e9 / n as f64,
        t_exact * 1e3,
        t_exact * 1e9 / n as f64,
        t_exact / t_pair
    );
}
