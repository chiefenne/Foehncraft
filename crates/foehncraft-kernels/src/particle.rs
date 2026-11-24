//! Kernel body of the parcel step of the Lagrangian sprays (L12, SP1, D245):
//! the velocity and the displacement of one parcel over a time interval in
//! which the gas state at the parcel is constant.
//!
//! The equation of motion of a small sphere in a gas (Maxey and Riley 1983,
//! reduced to drag, gravity and buoyancy; Crowe, Schwarzkopf, Sommerfeld and
//! Tsuji 2011) is
//!
//! dv/dt = (u − v)/τ + g (1 − ρ_g/ρ_p),   τ = ρ_p d² / (18 μ f(Re)),
//!
//! with the particle Reynolds number Re = ρ_g d |u − v| / μ and the drag
//! correction f = C_D Re/24 of Schiller and Naumann (1935): f = 1 + 0.15
//! Re^0.687 up to Re = 1000 and f = 0.0183 Re above (C_D = 0.44, Newton's
//! regime). With τ and the gas velocity frozen over the interval the
//! equation is linear and solved **exactly**:
//!
//! v(Δt) = a + (v₀ − a) e^(−Δt/τ),
//! x(Δt) − x₀ = a Δt + (v₀ − a) τ (1 − e^(−Δt/τ)),
//!
//! with a = u + g_eff τ, the velocity the parcel relaxes to. The step is
//! therefore not limited by the response time τ (a parcel with τ → 0 follows
//! the gas). The only error is the freezing of f(Re) and of the gas state,
//! which the tracker keeps small by cutting the interval at every cell
//! crossing and at a fraction of the relaxation time.
//!
//! The body is a gather over one parcel: no atomics, no reduction, and the
//! exponential and logarithm are the in-house routines of
//! `foehncraft_core::math`, so CPU and GPU agree bit for bit (D93 (6)).

use crate::ld;
use foehncraft_core::Real;

/// Drag correction f(Re) = C_D Re / 24 (Schiller and Naumann 1935; above
/// Re = 1000 the Newton regime C_D = 0.44).
#[inline(always)]
pub fn drag_correction<T: Real>(re: T) -> T {
    let one = T::ONE;
    let newton = T::from_f64(0.44 / 24.0) * re;
    if re > T::from_f64(1000.0) {
        newton
    } else if re > T::from_f64(1e-30) {
        // Re^0.687 = exp(0.687 ln Re)
        one + T::from_f64(0.15) * (T::from_f64(0.687) * re.ln_of()).exp_of()
    } else {
        one
    }
}

/// The state of one parcel and the gas at the parcel.
#[derive(Clone, Copy, Debug)]
pub struct ParcelInput<T> {
    /// Parcel velocity (m/s).
    pub v: [T; 3],
    /// Gas velocity at the parcel (m/s).
    pub u: [T; 3],
    /// Gas density (kg/m³).
    pub rho_g: T,
    /// Gas dynamic viscosity (Pa s).
    pub mu_g: T,
    /// Droplet diameter (m).
    pub d: T,
    /// Droplet density (kg/m³).
    pub rho_p: T,
    /// Gravity vector (m/s²).
    pub g: [T; 3],
    /// Interval (s).
    pub dt: T,
}

/// The result of [`parcel_step`].
#[derive(Clone, Copy, Debug)]
pub struct ParcelStep<T> {
    /// Velocity at the end of the interval (m/s).
    pub v: [T; 3],
    /// Displacement over the interval (m).
    pub dx: [T; 3],
    /// Response time τ at the start of the interval (s).
    pub tau: T,
    /// Particle Reynolds number at the start of the interval.
    pub re: T,
}

/// Velocity and displacement of one parcel over `dt`, by the exact solution
/// of the linear relaxation equation (module documentation).
#[inline]
pub fn parcel_step<T: Real>(i: &ParcelInput<T>) -> ParcelStep<T> {
    let mut rel2 = T::ZERO;
    for k in 0..3 {
        let w = i.u[k] - i.v[k];
        rel2 += w * w;
    }
    let re = i.rho_g * i.d * rel2.sqrt_of() / i.mu_g;
    let f = drag_correction(re);
    let tau = i.rho_p * i.d * i.d / (T::from_f64(18.0) * i.mu_g * f);
    let g_eff = T::ONE - i.rho_g / i.rho_p;
    let e = (-(i.dt / tau)).exp_of();
    let one_minus_e = T::ONE - e;
    let mut v = [T::ZERO; 3];
    let mut dx = [T::ZERO; 3];
    for k in 0..3 {
        let a = i.u[k] + i.g[k] * g_eff * tau;
        let w0 = i.v[k] - a;
        v[k] = a + w0 * e;
        dx[k] = a * i.dt + w0 * tau * one_minus_e;
    }
    ParcelStep { v, dx, tau, re }
}

/// Arrays of the parcel step, one entry per parcel (structure of arrays,
/// D47).
#[derive(Clone, Copy)]
pub struct ParcelArrays<'a, T> {
    /// Parcel velocity components.
    pub v: [&'a [T]; 3],
    /// Gas velocity components at the parcels.
    pub u: [&'a [T]; 3],
    /// Gas density per parcel.
    pub rho_g: &'a [T],
    /// Gas viscosity per parcel.
    pub mu_g: &'a [T],
    /// Droplet diameter per parcel.
    pub d: &'a [T],
    /// Droplet density per parcel.
    pub rho_p: &'a [T],
    /// Interval per parcel (s).
    pub dt: &'a [T],
    /// Gravity (m/s²).
    pub g: [T; 3],
}

/// Body of the parcel-step kernel for parcel `i`.
///
/// # Safety
/// `i` is smaller than the length of every array.
#[inline]
pub unsafe fn parcel_step_at<T: Real>(i: usize, a: &ParcelArrays<'_, T>) -> ParcelStep<T> {
    // SAFETY: guaranteed by the caller.
    unsafe {
        parcel_step(&ParcelInput {
            v: [ld(a.v[0], i), ld(a.v[1], i), ld(a.v[2], i)],
            u: [ld(a.u[0], i), ld(a.u[1], i), ld(a.u[2], i)],
            rho_g: ld(a.rho_g, i),
            mu_g: ld(a.mu_g, i),
            d: ld(a.d, i),
            rho_p: ld(a.rho_p, i),
            g: a.g,
            dt: ld(a.dt, i),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input(v: [f64; 3], u: [f64; 3], d: f64, dt: f64) -> ParcelInput<f64> {
        ParcelInput {
            v,
            u,
            rho_g: 1.2,
            mu_g: 1.8e-5,
            d,
            rho_p: 1000.0,
            g: [0.0, 0.0, -9.81],
            dt,
        }
    }

    #[test]
    fn stokes_relaxation_and_terminal_velocity() {
        // 10 µm droplet: Re ≪ 1, f ≈ 1, τ = ρ_p d²/(18 μ)
        let d = 10e-6;
        let tau = 1000.0 * d * d / (18.0 * 1.8e-5);
        let g_eff = 9.81 * (1.0 - 1.2 / 1000.0);
        // released at rest in still air: v_z(t) = −g_eff τ (1 − e^(−t/τ))
        let dt = 3.0 * tau;
        let s = parcel_step(&input([0.0; 3], [0.0; 3], d, dt));
        let exact = -g_eff * tau * (1.0 - (-3.0f64).exp());
        assert!((s.v[2] - exact).abs() < 1e-9 * g_eff * tau, "{} {exact}", s.v[2]);
        assert!((s.tau / tau - 1.0).abs() < 1e-4, "{} {tau}", s.tau);
        // a long step reaches the terminal velocity, however large τ/Δt
        let s = parcel_step(&input([0.0; 3], [0.0; 3], d, 1e3 * tau));
        assert!((s.v[2] + g_eff * tau).abs() < 1e-12);
    }

    #[test]
    fn relaxes_to_the_gas_velocity_without_gravity() {
        let mut i = input([0.0; 3], [2.0, 0.0, 0.0], 20e-6, 1e-3);
        i.g = [0.0; 3];
        let s = parcel_step(&i);
        let tau = s.tau;
        let exact = 2.0 * (1.0 - (-(1e-3) / tau).exp());
        assert!((s.v[0] - exact).abs() < 1e-12);
        // the displacement is the integral of the velocity
        let dx = 2.0 * 1e-3 - 2.0 * tau * (1.0 - (-(1e-3) / tau).exp());
        assert!((s.dx[0] - dx).abs() < 1e-15);
    }

    #[test]
    fn a_tiny_droplet_follows_the_gas() {
        let mut i = input([0.0; 3], [1.0, 2.0, 3.0], 1e-7, 1e-3);
        i.g = [0.0; 3];
        let s = parcel_step(&i);
        for k in 0..3 {
            assert!((s.v[k] - i.u[k]).abs() < 1e-9);
        }
    }

    #[test]
    fn drag_correction_is_continuous_at_the_newton_limit() {
        let lo = drag_correction(999.999_f64);
        let hi = drag_correction(1000.001_f64);
        // the two correlations meet within 0.5 % (Schiller and Naumann 1935)
        assert!((lo - hi).abs() / hi < 1e-2, "{lo} {hi}");
        assert!((drag_correction(1e-12_f64) - 1.0).abs() < 1e-8);
    }
}
