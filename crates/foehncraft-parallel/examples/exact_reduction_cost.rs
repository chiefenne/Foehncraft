//! Cost of the exact dot product of step W2 against the fixed-order block
//! sum it replaces (D93 asks for the measured cost).
//!
//! ```text
//! cargo run --release -p foehncraft-parallel --example exact_reduction_cost
//! ```
//!
//! For vectors of 10⁴ to 10⁶ values and 1 to all threads: nanoseconds per
//! value of `par_dot` (blocks of 256 in a fixed tree, the old path) and of
//! `par_exact_dot` (long accumulator), and their ratio. Output as a table and
//! as JSON lines (`{"n":…,"threads":…,"block_ns":…,"exact_ns":…}`).

use foehncraft_parallel::repro::{par_dot, par_exact_dot};
use std::hint::black_box;
use std::time::Instant;

fn data(n: usize, seed: u64) -> Vec<f64> {
    (0..n)
        .map(|i| foehncraft_core::rng::uniform([i as u64, seed, 1]) - 0.4)
        .collect()
}

fn time(reps: usize, mut f: impl FnMut() -> f64) -> f64 {
    let mut acc = 0.0;
    f(); // warm-up
    let t = Instant::now();
    for _ in 0..reps {
        acc += f();
    }
    black_box(acc);
    t.elapsed().as_secs_f64() / reps as f64
}

fn main() {
    let max_threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let mut threads = vec![1];
    if max_threads > 1 {
        threads.push(max_threads);
    }
    println!(
        "{:>9} {:>7} {:>10} {:>10} {:>6}",
        "n", "threads", "block ns", "exact ns", "ratio"
    );
    for n in [10_000usize, 100_000, 1_000_000] {
        let (a, b) = (data(n, 1), data(n, 2));
        let reps = (20_000_000 / n).max(5);
        for &t in &threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .unwrap();
            let (tb, te) = pool.install(|| {
                (
                    time(reps, || par_dot(black_box(&a), black_box(&b))),
                    time(reps, || par_exact_dot(black_box(&a), black_box(&b))),
                )
            });
            let (nb, ne) = (tb * 1e9 / n as f64, te * 1e9 / n as f64);
            println!("{n:>9} {t:>7} {nb:>10.3} {ne:>10.3} {:>6.1}", ne / nb);
            println!("{{\"n\":{n},\"threads\":{t},\"block_ns\":{nb:.4},\"exact_ns\":{ne:.4}}}");
        }
    }
}
