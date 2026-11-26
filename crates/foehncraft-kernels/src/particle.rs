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
    parcel_step_drag(i, T::ONE)
}

/// As [`parcel_step`] with the drag coefficient multiplied by `factor` (the
/// drag of a distorted drop, D303: C_D = C_D,sphere (1 + 2.632 y)); a factor
/// of one gives the bits of [`parcel_step`].
#[inline]
pub fn parcel_step_drag<T: Real>(i: &ParcelInput<T>, factor: T) -> ParcelStep<T> {
    let mut rel2 = T::ZERO;
    for k in 0..3 {
        let w = i.u[k] - i.v[k];
        rel2 += w * w;
    }
    let re = i.rho_g * i.d * rel2.sqrt_of() / i.mu_g;
    let f = drag_correction(re) * factor;
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

/// Nusselt number of a sphere in a gas stream, Nu = 2 + 0.6 Re^½ Pr^⅓
/// (Ranz and Marshall 1952).
#[inline(always)]
pub fn nusselt<T: Real>(re: T, pr: T) -> T {
    let two = T::from_f64(2.0);
    if re > T::from_f64(1e-30) && pr > T::from_f64(1e-30) {
        two + T::from_f64(0.6) * re.sqrt_of() * (pr.ln_of() / T::from_f64(3.0)).exp_of()
    } else {
        two
    }
}

/// The state of one parcel and the gas for the heat exchange (SP2b).
#[derive(Clone, Copy, Debug)]
pub struct HeatInput<T> {
    /// Parcel temperature (K).
    pub tp: T,
    /// Gas temperature at the parcel (K).
    pub tg: T,
    /// Particle Reynolds number (of the drag step).
    pub re: T,
    /// Gas Prandtl number μ c_p / k.
    pub pr: T,
    /// Gas thermal conductivity (W/(m K)).
    pub k_g: T,
    /// Droplet diameter (m).
    pub d: T,
    /// Droplet density (kg/m³).
    pub rho_p: T,
    /// Droplet specific heat (J/(kg K)).
    pub cp_p: T,
    /// Interval (s).
    pub dt: T,
}

/// The result of [`heat_step`].
#[derive(Clone, Copy, Debug)]
pub struct HeatStep<T> {
    /// Parcel temperature at the end of the interval (K).
    pub tp: T,
    /// Thermal response time τ_h = ρ_p c_p d² / (6 Nu k_g) (s).
    pub tau: T,
    /// exp(−Δt/τ_h).
    pub decay: T,
}

/// Temperature of one parcel over `dt` with the gas temperature and the
/// heat transfer coefficient frozen: the exact solution of
/// dT_p/dt = (T_g − T_p)/τ_h, T_p(Δt) = T_g + (T_p − T_g) e^(−Δt/τ_h)
/// (sensible heat only, no evaporation; a uniform droplet temperature,
/// Biot number small).
#[inline]
pub fn heat_step<T: Real>(i: &HeatInput<T>) -> HeatStep<T> {
    let nu = nusselt(i.re, i.pr);
    let tau = i.rho_p * i.cp_p * i.d * i.d / (T::from_f64(6.0) * nu * i.k_g);
    let decay = (-(i.dt / tau)).exp_of();
    HeatStep {
        tp: i.tg + (i.tp - i.tg) * decay,
        tau,
        decay,
    }
}

/// The eddy of the eddy-lifetime model (Gosman and Ioannides 1983).
#[derive(Clone, Copy, Debug)]
pub struct EddyTimes<T> {
    /// Eddy length scale l_e = C_μ^¾ k^{3/2} / ε (m).
    pub length: T,
    /// Eddy lifetime τ_e = l_e / (2k/3)^½ (s).
    pub lifetime: T,
    /// Interaction time: the shorter of the lifetime and the time the
    /// parcel needs to cross the eddy, t_cr = −τ_p ln(1 − l_e / (τ_p |w|))
    /// (the eddy is never crossed when l_e ≥ τ_p |w|) (s).
    pub interaction: T,
}

/// The times of an eddy with the turbulence `k`, `eps` at the parcel, the
/// parcel's momentum response time `tau_p` and its slip `slip` = |u + u' −
/// v| at the start of the eddy (Gosman and Ioannides 1983; C_μ = 0.09).
#[inline]
pub fn eddy_times<T: Real>(k: T, eps: T, tau_p: T, slip: T) -> EddyTimes<T> {
    let zero = T::from_f64(0.0);
    let c = T::from_f64(0.164_316_767_251_549_8); // 0.09^¾
    let length = c * k * k.sqrt_of() / eps;
    let lifetime = length / (T::from_f64(2.0 / 3.0) * k).sqrt_of();
    let reach = tau_p * slip;
    let interaction = if reach > length && reach > zero {
        let cross = -(tau_p * (T::from_f64(1.0) - length / reach).ln_of());
        if cross < lifetime { cross } else { lifetime }
    } else {
        lifetime
    };
    EddyTimes {
        length,
        lifetime,
        interaction,
    }
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

/// As [`parcel_step_at`] with the index checked: the entry for the CPU
/// backends outside this crate.
///
/// # Panics
/// When `i` is not smaller than the length of every array.
#[inline]
pub fn parcel_step_checked<T: Real>(i: usize, a: &ParcelArrays<'_, T>) -> ParcelStep<T> {
    let n = a.dt.len();
    assert!(
        i < n
            && a.v.iter().chain(a.u.iter()).all(|s| s.len() == n)
            && [a.rho_g, a.mu_g, a.d, a.rho_p].iter().all(|s| s.len() == n),
        "parcel {i} outside the arrays of the batch"
    );
    // SAFETY: i < n and every array has n entries (checked above).
    unsafe { parcel_step_at(i, a) }
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
        assert!(
            (s.v[2] - exact).abs() < 1e-9 * g_eff * tau,
            "{} {exact}",
            s.v[2]
        );
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

#[cfg(test)]
mod heat_tests {
    use super::*;

    #[test]
    fn eddy_times_follow_gosman_and_ioannides() {
        let (k, eps) = (0.6f64, 2.0f64);
        let l = 0.09f64.powf(0.75) * k.powf(1.5) / eps;
        let te = l / (2.0 * k / 3.0).sqrt();
        // a tracer never crosses its eddy: the lifetime
        let e = eddy_times(k, eps, 1e-9, 1.0);
        assert!((e.length - l).abs() < 1e-15 * l && (e.lifetime - te).abs() < 1e-14 * te);
        assert_eq!(e.interaction, e.lifetime);
        // a heavy parcel with a large slip crosses it first
        let (tau, w) = (1.0, 10.0 * l / te);
        let e = eddy_times(k, eps, tau, w);
        let cross = -tau * (1.0 - l / (tau * w)).ln();
        assert!(cross < te && (e.interaction - cross).abs() < 1e-14 * cross);
        // the crossing time tends to l/|w| for a parcel that is not
        // slowed in the eddy
        let e = eddy_times(k, eps, 1e6, w);
        assert!((e.interaction - l / w).abs() < 1e-6 * l / w);
    }

    #[test]
    fn heat_step_relaxes_exactly_and_halves_compose() {
        let base = HeatInput {
            tp: 293.15,
            tg: 573.15,
            re: 12.0,
            pr: 0.7,
            k_g: 0.044,
            d: 100e-6,
            rho_p: 998.0,
            cp_p: 4181.0,
            dt: 0.02,
        };
        let one = heat_step(&base);
        let a = heat_step(&HeatInput { dt: 0.01, ..base });
        let b = heat_step(&HeatInput {
            tp: a.tp,
            dt: 0.01,
            ..base
        });
        assert!((one.tp - b.tp).abs() < 1e-9 * one.tp);
        let nu = 2.0 + 0.6 * 12.0f64.sqrt() * 0.7f64.powf(1.0 / 3.0);
        assert!((nusselt(12.0, 0.7) - nu).abs() < 1e-12);
        let tau = 998.0 * 4181.0 * 1e-8 / (6.0 * nu * 0.044);
        assert!((one.tau - tau).abs() < 1e-12 * tau);
        // far beyond τ_h the parcel takes the gas temperature
        let long = heat_step(&HeatInput {
            dt: 100.0 * tau,
            ..base
        });
        assert!((long.tp - 573.15).abs() < 1e-9);
    }
}
