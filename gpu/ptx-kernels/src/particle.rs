//! The parcel step of the Lagrangian sprays (L12, SP1, D245):
//! `foehncraft_kernels::particle::parcel_step` as a PTX entry point, one GPU
//! thread per parcel. The host side (arrays of the parcels and of the gas at
//! the parcels, structure of arrays) is the batch backend of
//! `crates/foehncraft-particles/src/batch.rs`; the argument list there must
//! match this signature.

use crate::tid;
use core::slice::{from_raw_parts as sl, from_raw_parts_mut as sl_mut};
use foehncraft_kernels::particle::{self, ParcelArrays};

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn particle_step(
    n: u32,
    vx: *const f64,
    vy: *const f64,
    vz: *const f64,
    ux: *const f64,
    uy: *const f64,
    uz: *const f64,
    rho_g: *const f64,
    mu_g: *const f64,
    d: *const f64,
    rho_p: *const f64,
    dt: *const f64,
    gx: f64,
    gy: f64,
    gz: f64,
    out_vx: *mut f64,
    out_vy: *mut f64,
    out_vz: *mut f64,
    out_dx: *mut f64,
    out_dy: *mut f64,
    out_dz: *mut f64,
    out_tau: *mut f64,
    out_re: *mut f64,
) {
    let i = tid();
    let n = n as usize;
    if i >= n {
        return;
    }
    // SAFETY: the host passes device arrays of n entries.
    unsafe {
        let a = ParcelArrays {
            v: [sl(vx, n), sl(vy, n), sl(vz, n)],
            u: [sl(ux, n), sl(uy, n), sl(uz, n)],
            rho_g: sl(rho_g, n),
            mu_g: sl(mu_g, n),
            d: sl(d, n),
            rho_p: sl(rho_p, n),
            dt: sl(dt, n),
            g: [gx, gy, gz],
        };
        let s = particle::parcel_step_at(i, &a);
        sl_mut(out_vx, n)[i] = s.v[0];
        sl_mut(out_vy, n)[i] = s.v[1];
        sl_mut(out_vz, n)[i] = s.v[2];
        sl_mut(out_dx, n)[i] = s.dx[0];
        sl_mut(out_dy, n)[i] = s.dx[1];
        sl_mut(out_dz, n)[i] = s.dx[2];
        sl_mut(out_tau, n)[i] = s.tau;
        sl_mut(out_re, n)[i] = s.re;
    }
}
