//! The parcel step for a batch of parcels (structure of arrays): the CPU
//! backends of the kernel that the GPU entry `particle_step`
//! (`gpu/ptx-kernels/src/particle.rs`) runs, one thread per parcel. All three
//! call the same body, `foehncraft_kernels::particle::parcel_step_at`, with
//! the in-house exponential and logarithm, so they agree bit for bit (D93
//! (6)); the test `backends_agree_bit_for_bit` checks the two CPU loops
//! against the body, and CI compiles the PTX entry and checks it for atomics
//! and fused arithmetic (R10).

use foehncraft_kernels::particle::{ParcelArrays, parcel_step_checked};
use rayon::prelude::*;

/// The results of a batch, one entry per parcel.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Batch {
    /// Velocity after the step.
    pub v: [Vec<f64>; 3],
    /// Displacement.
    pub dx: [Vec<f64>; 3],
    /// Response time at the start (s).
    pub tau: Vec<f64>,
    /// Particle Reynolds number at the start.
    pub re: Vec<f64>,
}

impl Batch {
    fn zeros(n: usize) -> Self {
        let z = || vec![0.0; n];
        Self {
            v: [z(), z(), z()],
            dx: [z(), z(), z()],
            tau: z(),
            re: z(),
        }
    }

    fn put(&mut self, i: usize, s: &foehncraft_kernels::particle::ParcelStep<f64>) {
        for k in 0..3 {
            self.v[k][i] = s.v[k];
            self.dx[k][i] = s.dx[k];
        }
        self.tau[i] = s.tau;
        self.re[i] = s.re;
    }
}

fn check(a: &ParcelArrays<f64>) -> usize {
    let n = a.dt.len();
    let ok = a.v.iter().chain(a.u.iter()).all(|s| s.len() == n)
        && [a.rho_g, a.mu_g, a.d, a.rho_p].iter().all(|s| s.len() == n);
    assert!(ok, "the arrays of a parcel batch must have one length");
    n
}

/// One loop over the parcels.
pub fn step_sequential(a: &ParcelArrays<f64>) -> Batch {
    let n = check(a);
    let mut out = Batch::zeros(n);
    for i in 0..n {
        let s = parcel_step_checked(i, a);
        out.put(i, &s);
    }
    out
}

/// The loop split over the threads of the pool.
pub fn step_parallel(a: &ParcelArrays<f64>) -> Batch {
    let n = check(a);
    let steps: Vec<_> = (0..n)
        .into_par_iter()
        .map(|i| parcel_step_checked(i, a))
        .collect();
    let mut out = Batch::zeros(n);
    for (i, s) in steps.iter().enumerate() {
        out.put(i, s);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use foehncraft_core::rng;
    use foehncraft_kernels::particle::{ParcelInput, parcel_step};

    #[test]
    fn backends_agree_bit_for_bit() {
        let n = 5000;
        let r = |k: u64, i: usize, lo: f64, hi: f64| lo + (hi - lo) * rng::uniform([77, k, i as u64]);
        let col = |k: u64, lo: f64, hi: f64| (0..n).map(|i| r(k, i, lo, hi)).collect::<Vec<_>>();
        let v = [col(1, -5.0, 5.0), col(2, -5.0, 5.0), col(3, -5.0, 5.0)];
        let u = [col(4, -5.0, 5.0), col(5, -5.0, 5.0), col(6, -5.0, 5.0)];
        let rho_g = col(7, 0.5, 2.0);
        let mu_g = col(8, 1e-5, 4e-5);
        // diameters over four decades: Stokes to Newton regimes
        let d: Vec<f64> = (0..n).map(|i| 10f64.powf(r(9, i, -6.0, -2.5))).collect();
        let rho_p = col(10, 500.0, 1500.0);
        let dt = col(11, 1e-6, 1e-1);
        let g = [0.0, -9.81, 0.0];
        let a = ParcelArrays {
            v: [&v[0], &v[1], &v[2]],
            u: [&u[0], &u[1], &u[2]],
            rho_g: &rho_g,
            mu_g: &mu_g,
            d: &d,
            rho_p: &rho_p,
            dt: &dt,
            g,
        };
        let seq = step_sequential(&a);
        for threads in [1, 2, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let par = pool.install(|| step_parallel(&a));
            assert!(seq == par, "{threads} threads differ from the sequential loop");
        }
        // and the per-parcel body
        for i in (0..n).step_by(97) {
            let s = parcel_step(&ParcelInput {
                v: [v[0][i], v[1][i], v[2][i]],
                u: [u[0][i], u[1][i], u[2][i]],
                rho_g: rho_g[i],
                mu_g: mu_g[i],
                d: d[i],
                rho_p: rho_p[i],
                g,
                dt: dt[i],
            });
            assert_eq!(s.v[0].to_bits(), seq.v[0][i].to_bits());
            assert_eq!(s.dx[2].to_bits(), seq.dx[2][i].to_bits());
            assert!(s.v.iter().chain(s.dx.iter()).all(|x| x.is_finite()));
        }
    }
}
