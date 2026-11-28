//! Kernel bodies of the wall boiling of the two-fluid model (L3, stage MP4,
//! D357, D358): the heat flux partitioning of Kurul & Podowski at a heated
//! wall, the heat and mass transfer at the bubble interfaces in the bulk,
//! and the sources of the energy equations of both phases.
//!
//! As the other kernel modules, every body computes the outputs of one
//! index from its inputs (a gather), so CPU loops and GPU threads give
//! bitwise the same values (D31).
//!
//! **Wall** (Kurul & Podowski 1991; the closures as summarised by Končar et
//! al. 2025, arXiv:2509.17637, and Fu & Anglart 2017, arXiv:1709.01783).
//! The wall heat flux q_w of a heated face is split into single-phase
//! convection, quenching and evaporation,
//!
//!   q_w = (1 − A) h_c (T_w − T_P) + A h_q (T_w − T_ref) + q_e,
//!
//! and the wall temperature T_w is the root of that sum (bracketed, then
//! bisected; a fixed count, no data-dependent loop length beyond the
//! bracket). With ΔT_sup = T_w − T_sat and ΔT_sub = T_sat − T_ref:
//!
//! - h_c = ρ_l c_p,l u_τ / T⁺(y⁺_P) of the near-wall cell P, T⁺ of Kader
//!   (1981), u_τ = C_μ^¼ √k, at least the conduction k_l/y_P;
//! - T_ref, the liquid temperature of the quenching and of the subcooling,
//!   at a fixed y⁺ (250 by default, Frank et al. 2010) by the ratio of the
//!   Kader profiles: T_ref = T_w − T⁺(y⁺_ref)/T⁺(y⁺_P) (T_w − T_P); a
//!   reference y⁺ of 0 takes T_P;
//! - nucleation site density n = (m ΔT_sup)^p, m = 185, p = 1.805
//!   (Lemmert & Chawla 1977);
//! - departure diameter d_w = min(d_ref exp(−ΔT_sub/ΔT_ref), d_max), 0.6
//!   mm, 45 K, 1.4 mm (Tolubinsky & Kostanchuk 1970);
//! - departure frequency f = √(4 g (ρ_l − ρ_g)/(3 C_D d_w ρ_l)), C_D = 1
//!   (Cole 1960); waiting time t_w = 0.8/f;
//! - bubble influence area A = min(K π d_w² n/4, A_max), K = 4 (an
//!   influence diameter of two departure diameters, Del Valle & Kenning
//!   1985);
//! - quenching h_q = 2 k_l f √(t_w/(π a_l)), a_l = k_l/(ρ_l c_p,l) (Del
//!   Valle & Kenning 1985);
//! - evaporation q_e = (π/6) d_w³ ρ_g h_fg f n; the vapour made per unit
//!   wall area is q_e/(h_fg + c_p,l (T_sat − T_P)), so that the
//!   liquid leaving the near-wall cell at its own enthalpy and the vapour
//!   made at saturation carry exactly q_e.
//!
//! **Dryout** (the vapour energy equation, D357). As the vapour fraction α
//! of the near-wall cell passes the break point α_dry (0.9), a dry part
//! K_dry = β²(3 − 2β), β = (α − α_dry)/(1 − α_dry), of the wall heats the
//! vapour, K_dry h_g (T_w − T_g) with h_g = ρ_g c_p,g u_τ/T⁺ of the vapour,
//! and the three liquid terms above take the wet part 1 − K_dry (Yoon,
//! Agostinelli & Baglietto 2017, INL/CON-16-40542, eqs. 9, 11, 12). A wall
//! cell full of vapour then stops making vapour and superheats it.
//!
//! **Interface** (two resistances, an interface at T_sat). The heat to the
//! interface from phase k is H_k (T_k − T_sat), H_k = Nu_k (k_k/d) (6 α_g
//! /d), with Nu_l of Ranz & Marshall (1952), 2 + 0.6 Re^½ Pr^⅓, and Nu_g = 2
//! inside the bubble (no internal convection, as Shaver et al., OSTI
//! 1607420). The net heat makes vapour, Γ_i = (H_l (T_l − T_sat) + H_g (T_g
//! − T_sat))/L (condensation where negative, L below). Γ_i is proportional to α_g;
//! a condensing cell takes it with the α_g of the new level, so α_g stays
//! positive for any step.
//!
//! **Bubble size** in the bulk: linear in the local subcooling between
//! d_cold at ΔT_cold and d_hot at ΔT_hot, constant outside (Kurul &
//! Podowski 1991, Anglart et al. 1997; 0.1 mm at 13.5 K and 2 mm at −5 K,
//! as given by Frank et al. 2010).
//!
//! **Energy sources** (the temperature form of each phase, D357). With the
//! mass leaving a phase at the interface enthalpy and arriving in the
//! other at saturation, the non-conservative energy equations take
//!
//!   liquid: H_l (T_sat − T_l) + Γ_→l c_p,l (T_sat − T_l) + q‴_wall,
//!   vapour: H_g (T_sat − T_g) + Γ_→g c_p,g (T_sat − T_g),
//!
//! Γ_→l the condensate, Γ_→g the vapour made (interface and wall); mass
//! leaves a phase at its own enthalpy (no term) and arrives at saturation.
//! Γ_i therefore divides the net interface heat by the enthalpy jump of
//! the mass crossing it, h_fg + c_p,l (T_sat − T_l) evaporating and h_fg +
//! c_p,g (T_g − T_sat) condensing, and every term damps (implicit). Liquid
//! and vapour then exchange exactly the interface heat.

use crate::turb::thermal;
use foehncraft_core::Real;

/// C_μ of u_τ = C_μ^¼ √k.
pub const C_MU: f64 = 0.09;

/// The largest α_g after the phase change of a step.
pub const ALPHA_MAX: f64 = 0.999;

/// The constants of the boiling (saturation state and closures).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct Constants<T> {
    /// Saturation temperature T_sat (K).
    pub t_sat: T,
    /// Enthalpy of evaporation h_fg (J/kg).
    pub h_fg: T,
    /// Liquid density, specific heat, conductivity and viscosity.
    pub rho_l: T,
    /// See `rho_l`.
    pub cp_l: T,
    /// See `rho_l`.
    pub k_l: T,
    /// See `rho_l`.
    pub mu_l: T,
    /// Vapour density, specific heat, conductivity and viscosity.
    pub rho_g: T,
    /// See `rho_g`.
    pub cp_g: T,
    /// See `rho_g`.
    pub k_g: T,
    /// See `rho_g`.
    pub mu_g: T,
    /// |g| (m/s²).
    pub g: T,
    /// m and p of the site density (m ΔT)^p (1/m², ΔT in K).
    pub site_m: T,
    /// See `site_m`.
    pub site_p: T,
    /// d_ref, ΔT_ref and d_max of the departure diameter (m, K, m).
    pub dw_ref: T,
    /// See `dw_ref`.
    pub dw_dt: T,
    /// See `dw_ref`.
    pub dw_max: T,
    /// C_D of the departure frequency.
    pub cole_cd: T,
    /// The influence factor K of the bubble area and its largest value.
    pub influence: T,
    /// See `influence`.
    pub area_max: T,
    /// The waiting time as a fraction of 1/f.
    pub wait: T,
    /// y⁺ of the reference liquid temperature (0: the near-wall cell).
    pub y_ref: T,
    /// The bulk bubble size: d_cold at the subcooling dt_cold, d_hot at
    /// dt_hot (dt_cold > dt_hot; equal diameters: a constant size).
    pub d_cold: T,
    /// See `d_cold`.
    pub dt_cold: T,
    /// See `d_cold`.
    pub d_hot: T,
    /// See `d_cold`.
    pub dt_hot: T,
    /// Nu of the vapour side of the interface.
    pub nu_g: T,
    /// The dryout break point α_dry of the near-wall vapour fraction
    /// (1: no dryout).
    pub alpha_dry: T,
    /// The residual phase fraction.
    pub alpha_min: T,
}

/// x^p for x > 0 (0 for x ≤ 0).
#[inline(always)]
fn pow<T: Real>(x: T, p: T) -> T {
    if x > T::ZERO {
        (x.ln_of() * p).exp_of()
    } else {
        T::ZERO
    }
}

/// The bulk bubble diameter at the liquid temperature `t_l`.
#[inline(always)]
pub fn bubble_diameter<T: Real>(c: &Constants<T>, t_l: T) -> T {
    let sub = c.t_sat - t_l;
    if c.dt_cold <= c.dt_hot || sub >= c.dt_cold {
        c.d_cold
    } else if sub <= c.dt_hot {
        c.d_hot
    } else {
        (c.d_cold * (sub - c.dt_hot) + c.d_hot * (c.dt_cold - sub)) / (c.dt_cold - c.dt_hot)
    }
}

/// The departure diameter d_w (m) at the subcooling `sub` (K).
#[inline(always)]
pub fn departure_diameter<T: Real>(c: &Constants<T>, sub: T) -> T {
    (c.dw_ref * (-(sub / c.dw_dt)).exp_of()).min_of(c.dw_max)
}

/// The departure frequency f (1/s) of the departure diameter `dw`.
#[inline(always)]
pub fn departure_frequency<T: Real>(c: &Constants<T>, dw: T) -> T {
    (T::from_f64(4.0) * c.g * (c.rho_l - c.rho_g) / (T::from_f64(3.0) * c.cole_cd * dw * c.rho_l))
        .sqrt_of()
}

/// The heat to make 1 kg of saturated vapour from liquid at `t_l`:
/// h_g,sat − h_l = h_fg + c_p,l (T_sat − T_l), at least h_fg/10.
#[inline(always)]
pub fn latent_liquid<T: Real>(c: &Constants<T>, t_l: T) -> T {
    (c.h_fg + c.cp_l * (c.t_sat - t_l)).max_of(T::from_f64(0.1) * c.h_fg)
}

/// The heat released when 1 kg of vapour at `t_g` condenses to saturated
/// liquid: h_g − h_l,sat = h_fg + c_p,g (T_g − T_sat), at least h_fg/10.
#[inline(always)]
pub fn latent_vapour<T: Real>(c: &Constants<T>, t_g: T) -> T {
    (c.h_fg + c.cp_g * (t_g - c.t_sat)).max_of(T::from_f64(0.1) * c.h_fg)
}

/// The near-wall state of a heated face.
#[derive(Clone, Copy, Debug)]
pub struct WallIn<T> {
    /// Wall heat flux q_w (W/m², into the fluid).
    pub q_w: T,
    /// Liquid temperature of the near-wall cell T_P (K).
    pub t_p: T,
    /// Turbulent kinetic energy of the near-wall cell (m²/s²).
    pub k: T,
    /// Distance of the cell centre from the wall y_P (m).
    pub y: T,
    /// Vapour fraction and vapour temperature of the near-wall cell.
    pub alpha: T,
    /// See `alpha`.
    pub t_g: T,
}

/// The split of the wall heat flux.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WallSplit<T> {
    /// Wall temperature T_w (K).
    pub t_w: T,
    /// Convection, quenching and evaporation (W/m²).
    pub q_c: T,
    /// See `q_c`.
    pub q_q: T,
    /// See `q_c`.
    pub q_e: T,
    /// The bubble influence area fraction A.
    pub area: T,
    /// The vapour made per unit wall area (kg/(m² s)).
    pub m_e: T,
    /// Convection to the vapour K_dry h_g (T_w − T_g) (W/m²).
    pub q_g: T,
    /// The dry fraction K_dry of the wall.
    pub dry: T,
}

/// The dry fraction of the wall K_dry = β² (3 − 2β), β = (α − α_dry)/(1 −
/// α_dry) above the break point α_dry, 0 below (Yoon, Agostinelli &
/// Baglietto 2017, INL/CON-16-40542, eq. 12, α_dry = 0.9; the α of the
/// near-wall cell, their "wall cell" layer).
#[inline(always)]
pub fn dry_fraction<T: Real>(c: &Constants<T>, alpha: T) -> T {
    if c.alpha_dry >= T::ONE || alpha <= c.alpha_dry {
        return T::ZERO;
    }
    let b = ((alpha - c.alpha_dry) / (T::ONE - c.alpha_dry)).min_of(T::ONE);
    b * b * (T::from_f64(3.0) - T::from_f64(2.0) * b)
}

/// The wall-function conductance ρ c_p u_τ/T⁺(y⁺) of a phase at the wall
/// distance `y`, at least the conduction k/y.
#[inline(always)]
fn wall_conductance<T: Real>(rho: T, cp: T, k: T, mu: T, ut: T, y: T) -> T {
    let yp = rho * ut * y / mu;
    let tp = thermal::kader(yp, mu * cp / k);
    if tp > T::ZERO {
        (rho * cp * ut / tp).max_of(k / y)
    } else {
        k / y
    }
}

/// The parts of the wall heat flux at the wall temperature `t_w`, with
/// h_c and the ratio r = T⁺(y⁺_ref)/T⁺(y⁺_P) of the face.
#[inline(always)]
fn parts<T: Real>(
    c: &Constants<T>,
    w: &WallIn<T>,
    h_c: T,
    h_g: T,
    dry: T,
    r: T,
    t_w: T,
) -> WallSplit<T> {
    let t_ref = t_w - r * (t_w - w.t_p);
    let sup = t_w - c.t_sat;
    let sub = c.t_sat - t_ref;
    let dw = departure_diameter(c, sub);
    let f = departure_frequency(c, dw);
    let n = pow(c.site_m * sup, c.site_p);
    let pi = T::from_f64(core::f64::consts::PI);
    let area = (c.influence * pi * dw * dw * n / T::from_f64(4.0)).min_of(c.area_max);
    let tw = c.wait / f;
    let a_l = c.k_l / (c.rho_l * c.cp_l);
    let h_q = T::from_f64(2.0) * c.k_l * f * (tw / (pi * a_l)).sqrt_of();
    let q_e = pi / T::from_f64(6.0) * dw * dw * dw * c.rho_g * c.h_fg * f * n;
    let wet = T::ONE - dry;
    let q_e = wet * q_e;
    let q_c = wet * (T::ONE - area) * h_c * (t_w - w.t_p);
    let q_q = wet * area * h_q * (t_w - t_ref);
    let m_e = q_e / latent_liquid(c, w.t_p);
    WallSplit {
        t_w,
        q_c,
        q_q,
        q_e,
        area,
        m_e,
        q_g: dry * h_g * (t_w - w.t_g),
        dry,
    }
}

/// The heat flux partitioning of Kurul & Podowski at one heated face: the
/// wall temperature where convection, quenching and evaporation of the wet
/// part (1 − K_dry) and convection to the vapour of the dry part K_dry
/// carry q_w,
///
///   q_w = (1 − K_dry)(q_c + q_q + q_e) + K_dry h_g (T_w − T_g)
///
/// (Yoon et al. 2017 eq. 9, h_g the wall function of the vapour, their eq.
/// 11).
pub fn wall_split<T: Real>(c: &Constants<T>, w: &WallIn<T>) -> WallSplit<T> {
    let pr = c.mu_l * c.cp_l / c.k_l;
    let y = w.y.max_of(T::from_f64(1e-12));
    let ut = T::from_f64(C_MU).sqrt_of().sqrt_of() * w.k.max_of(T::ZERO).sqrt_of();
    let yp = c.rho_l * ut * y / c.mu_l;
    let tp_p = thermal::kader(yp, pr);
    let h_c = wall_conductance(c.rho_l, c.cp_l, c.k_l, c.mu_l, ut, y);
    let dry = dry_fraction(c, w.alpha);
    let h_g = if dry > T::ZERO {
        wall_conductance(c.rho_g, c.cp_g, c.k_g, c.mu_g, ut, y)
    } else {
        T::ZERO
    };
    let r = if c.y_ref > T::ZERO && tp_p > T::ZERO {
        thermal::kader(c.y_ref, pr) / tp_p
    } else {
        T::ONE
    };
    let total = |t: T| {
        let s = parts(c, w, h_c, h_g, dry, r, t);
        s.q_c + s.q_q + s.q_e + s.q_g - w.q_w
    };
    // the bracket: no evaporation at or below T_sat, no convection at T_P
    let wet = T::ONE - dry;
    let h = wet * h_c + dry * h_g;
    let lo0 = w.t_p.min_of(c.t_sat);
    let lo0 = if dry > T::ZERO {
        lo0.min_of(w.t_g)
    } else {
        lo0
    };
    let mut lo = lo0;
    let mut step = w.q_w.abs() / h + T::ONE;
    let hi0 = w.t_p.max_of(c.t_sat);
    let hi0 = if dry > T::ZERO {
        hi0.max_of(w.t_g)
    } else {
        hi0
    };
    let mut hi = hi0 + step;
    let mut i = 0;
    while total(hi) < T::ZERO && i < 60 {
        lo = hi;
        step = step + step;
        hi = hi + step;
        i += 1;
    }
    if total(lo) > T::ZERO {
        // (a cooled wall, q_w < 0, below both phases: convection alone)
        let t_w = (w.q_w + wet * h_c * w.t_p + dry * h_g * w.t_g) / h;
        return WallSplit {
            t_w,
            q_c: wet * h_c * (t_w - w.t_p),
            q_g: dry * h_g * (t_w - w.t_g),
            dry,
            ..Default::default()
        };
    }
    let mut k = 0;
    while k < 64 {
        let mid = T::from_f64(0.5) * (lo + hi);
        if total(mid) < T::ZERO {
            lo = mid;
        } else {
            hi = mid;
        }
        k += 1;
    }
    parts(c, w, h_c, h_g, dry, r, T::from_f64(0.5) * (lo + hi))
}

/// What [`phase_change`] gives per cell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CellOut<T> {
    /// α_g of the new level.
    pub alpha: T,
    /// The bulk bubble diameter (m).
    pub d: T,
    /// Γ_i and Γ = Γ_i + Γ_w (kg/(m³ s), vapour made).
    pub gamma_i: T,
    /// See `gamma_i`.
    pub gamma: T,
    /// H_l, H_g (W/(m³ K)).
    pub h_l: T,
    /// See `h_l`.
    pub h_g: T,
}

/// The phase change of one cell over the step `dt`: from α* (after the
/// advection), the temperatures and the slip of the step start and the
/// wall vapour `m_w` (kg/(m³ s)), the new α_g and the rates.
pub fn phase_change<T: Real>(
    c: &Constants<T>,
    alpha: T,
    t_l: T,
    t_g: T,
    slip: T,
    m_w: T,
    dt: T,
) -> CellOut<T> {
    let a = alpha.max_of(T::ZERO);
    let d = bubble_diameter(c, t_l);
    let re = c.rho_l * slip * d / c.mu_l;
    let pr = c.mu_l * c.cp_l / c.k_l;
    let nu_l = T::from_f64(2.0)
        + T::from_f64(0.6) * re.max_of(T::ZERO).sqrt_of() * pow(pr, T::from_f64(1.0 / 3.0));
    let six = T::from_f64(6.0) / d;
    let hl = nu_l * c.k_l / d * six;
    let hg = c.nu_g * c.k_g / d * six;
    // Γ_i = α r: the net heat to the interface over the enthalpy jump of
    // the mass that crosses it (leaving its phase at that phase's own
    // enthalpy, arriving at saturation)
    let q = hl * (t_l - c.t_sat) + hg * (t_g - c.t_sat);
    let r = if q > T::ZERO {
        q / latent_liquid(c, t_l)
    } else {
        q / latent_vapour(c, t_g)
    };
    let (an, ai) = if r < T::ZERO {
        let n = (a + dt * m_w / c.rho_g) / (T::ONE - dt * r / c.rho_g);
        (n, n)
    } else {
        (a + dt * (a * r + m_w) / c.rho_g, a)
    };
    let an = an.min_of(T::from_f64(ALPHA_MAX));
    let gamma = c.rho_g * (an - a) / dt;
    CellOut {
        alpha: an,
        d,
        gamma_i: gamma - m_w,
        gamma,
        h_l: hl * ai,
        h_g: hg * ai,
    }
}

/// The sources of the two temperature equations of one cell, per unit of
/// the heat capacity α_k ρ_k c_p,k (the S_u, S_p of the scalar rows):
/// (S_u,l, S_p,l, S_u,g, S_p,g, capacity_l, capacity_g). `q_wall` and
/// `q_wall_g` are the wall heat of the liquid and of the vapour (W/m³),
/// `m_w` the wall vapour (kg/(m³ s)). Mass
/// leaves a phase at its own enthalpy (no term in the temperature form)
/// and arrives at saturation, which pulls the receiving phase towards
/// T_sat: every exchange term damps, all are implicit.
pub fn energy_sources<T: Real>(
    c: &Constants<T>,
    o: &CellOut<T>,
    q_wall: T,
    q_wall_g: T,
    m_w: T,
) -> [T; 6] {
    let cap_l = (T::ONE - o.alpha).max_of(c.alpha_min) * c.rho_l * c.cp_l;
    let cap_g = o.alpha.max_of(c.alpha_min) * c.rho_g * c.cp_g;
    let ts = c.t_sat;
    // the mass into the liquid (condensation) and into the vapour
    let to_l = (-o.gamma_i).max_of(T::ZERO);
    let to_g = o.gamma_i.max_of(T::ZERO) + m_w.max_of(T::ZERO);
    let al = o.h_l + to_l * c.cp_l;
    let ag = o.h_g + to_g * c.cp_g;
    [
        (al * ts + q_wall) / cap_l,
        -al / cap_l,
        (ag * ts + q_wall_g) / cap_g,
        -ag / cap_g,
        cap_l,
        cap_g,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Water at 6.84 MPa (IAPWS-IF97 saturation, rounded).
    fn water() -> Constants<f64> {
        Constants {
            t_sat: 557.0,
            h_fg: 1.51e6,
            rho_l: 740.0,
            cp_l: 5400.0,
            k_l: 0.56,
            mu_l: 9.2e-5,
            rho_g: 35.6,
            cp_g: 5900.0,
            k_g: 0.06,
            mu_g: 1.9e-5,
            g: 9.81,
            site_m: 185.0,
            site_p: 1.805,
            dw_ref: 0.6e-3,
            dw_dt: 45.0,
            dw_max: 1.4e-3,
            cole_cd: 1.0,
            influence: 4.0,
            area_max: 1.0,
            wait: 0.8,
            y_ref: 250.0,
            d_cold: 0.1e-3,
            dt_cold: 13.5,
            d_hot: 2.0e-3,
            dt_hot: -5.0,
            nu_g: 2.0,
            alpha_dry: 0.9,
            alpha_min: 1e-6,
        }
    }

    fn wall(t_p: f64, q_w: f64) -> WallIn<f64> {
        WallIn {
            q_w,
            t_p,
            k: 0.01,
            y: 1e-4,
            alpha: 0.0,
            t_g: 557.0,
        }
    }

    /// A wall cell full of vapour: the wet part of the split vanishes, the
    /// wall heats the vapour, and the split still closes q_w.
    #[test]
    fn dry_wall_heats_the_vapour() {
        let c = water();
        assert_eq!(dry_fraction(&c, 0.5), 0.0);
        assert_eq!(dry_fraction(&c, 0.9), 0.0);
        assert!((dry_fraction(&c, 0.95) - 0.5).abs() < 1e-12);
        assert!((dry_fraction(&c, 1.0) - 1.0).abs() < 1e-12);
        let q = 1.13e6;
        let mut prev = 0.0;
        for a in [0.91, 0.95, 0.99, 0.999] {
            let w = WallIn {
                alpha: a,
                ..wall(556.0, q)
            };
            let s = wall_split(&c, &w);
            let sum = s.q_c + s.q_q + s.q_e + s.q_g;
            assert!((sum - q).abs() < 1e-8 * q, "{a}: {sum} vs {q}");
            assert!(s.q_g > prev, "{a}: vapour share {} after {prev}", s.q_g);
            assert!((s.m_e * latent_liquid(&c, 556.0) - s.q_e).abs() < 1e-6 * q);
            prev = s.q_g;
        }
        // a hotter wall below a vapour blanket than below wet boiling
        let wet = wall_split(&c, &wall(556.0, q));
        let dry = wall_split(
            &c,
            &WallIn {
                alpha: 0.999,
                ..wall(556.0, q)
            },
        );
        assert!(dry.t_w > wet.t_w + 5.0, "{} vs {}", dry.t_w, wet.t_w);
    }

    #[test]
    fn split_closes_the_wall_flux() {
        let c = water();
        for (t_p, q) in [
            (500.0, 1.13e6),
            (540.0, 1.13e6),
            (555.0, 0.44e6),
            (450.0, 2e5),
        ] {
            let s = wall_split(&c, &wall(t_p, q));
            let sum = s.q_c + s.q_q + s.q_e;
            assert!((sum - q).abs() < 1e-8 * q, "{t_p}: {sum} vs {q}");
            assert!(s.t_w > t_p);
        }
    }

    #[test]
    fn cold_liquid_does_not_boil() {
        let c = water();
        let s = wall_split(&c, &wall(400.0, 2e5));
        assert!(s.t_w < c.t_sat, "{}", s.t_w);
        assert_eq!(s.q_e, 0.0);
        assert_eq!(s.area, 0.0);
    }

    #[test]
    fn boiling_holds_the_wall_near_saturation() {
        // a heated wall over nearly saturated liquid: a superheat of a few
        // to ~20 K, most of the flux into quenching and evaporation
        let c = water();
        let s = wall_split(&c, &wall(553.0, 1.13e6));
        let sup = s.t_w - c.t_sat;
        assert!(sup > 1.0 && sup < 30.0, "superheat {sup}");
        assert!(s.q_e > 0.0 && s.m_e > 0.0);
        // the vapour made carries q_e with the liquid's subcooling
        let back = s.m_e * latent_liquid(&c, 553.0);
        assert!((back - s.q_e).abs() < 1e-9 * s.q_e);
    }

    #[test]
    fn bubble_size_follows_the_subcooling() {
        let c = water();
        assert_eq!(bubble_diameter(&c, c.t_sat - 20.0), 0.1e-3);
        assert_eq!(bubble_diameter(&c, c.t_sat + 10.0), 2.0e-3);
        let mid = bubble_diameter(&c, c.t_sat - 4.25);
        assert!((mid - 1.05e-3).abs() < 1e-12, "{mid}");
    }

    #[test]
    fn condensation_keeps_alpha_positive() {
        // very cold liquid and a long step: the implicit condensation
        // shrinks α but never below zero
        let c = water();
        let o = phase_change(&c, 0.2, 450.0, c.t_sat, 0.1, 0.0, 10.0);
        assert!(o.alpha > 0.0 && o.alpha < 0.2, "{}", o.alpha);
        assert!(o.gamma < 0.0);
        assert!((o.gamma - c.rho_g * (o.alpha - 0.2) / 10.0).abs() < 1e-12);
    }

    #[test]
    fn interface_exchanges_the_latent_heat() {
        // the heat both phases give the interface is Γ_i h_fg (no wall)
        let c = water();
        let (tl, tg) = (540.0, 560.0);
        let o = phase_change(&c, 0.1, tl, tg, 0.2, 0.0, 1e-5);
        let q = o.h_l * (tl - c.t_sat) + o.h_g * (tg - c.t_sat);
        // (condensing: the vapour leaves at its own enthalpy)
        assert!(q < 0.0);
        assert!((q - o.gamma_i * latent_vapour(&c, tg)).abs() < 1e-6 * q.abs());
    }
}
