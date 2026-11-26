//! Friction, void, heat transfer and subcooled boiling correlations of the
//! tube circuits. All in SI units; every function names its source.

use crate::G;

/// The Darcy friction factor by Churchill (1977), valid in laminar,
/// transitional and turbulent flow of smooth and rough tubes; 64/Re in
/// laminar flow to round-off.
///
/// Churchill, S. W. (1977). Friction-factor equation spans all
/// fluid-flow regimes. Chemical Engineering 84(24), 91–92.
pub fn darcy_churchill(re: f64, rel_roughness: f64) -> f64 {
    let re = re.max(1e-12);
    let a = (2.457 * (1.0 / ((7.0 / re).powf(0.9) + 0.27 * rel_roughness)).ln()).powi(16);
    let b = (37530.0 / re).powi(16);
    8.0 * ((8.0 / re).powi(12) + (a + b).powf(-1.5)).powf(1.0 / 12.0)
}

/// The two-phase multiplier φ²_lo of Friedel (1979) on the pressure
/// gradient of the whole flow as liquid.
///
/// Friedel, L. (1979). Improved friction pressure drop correlations for
/// horizontal and vertical two-phase pipe flow. European Two-Phase Flow
/// Group Meeting, Ispra, paper E2. In the form of Whalley (1987) and the
/// VDI Heat Atlas (2010), section L2.1.
#[allow(clippy::too_many_arguments)]
pub fn friedel_multiplier(
    x: f64,
    g_flux: f64,
    d: f64,
    rho_l: f64,
    rho_g: f64,
    mu_l: f64,
    mu_g: f64,
    sigma: f64,
    f_lo: f64,
    f_go: f64,
) -> f64 {
    let x = x.clamp(0.0, 1.0);
    let e = (1.0 - x).powi(2) + x * x * rho_l * f_go / (rho_g * f_lo);
    let f = x.powf(0.78) * (1.0 - x).powf(0.224);
    let h = (rho_l / rho_g).powf(0.91) * (mu_g / mu_l).powf(0.19) * (1.0 - mu_g / mu_l).powf(0.7);
    let rho_h = 1.0 / (x / rho_g + (1.0 - x) / rho_l);
    let fr = g_flux * g_flux / (G * d * rho_h * rho_h);
    let we = g_flux * g_flux * d / (sigma * rho_h);
    e + 3.24 * f * h / (fr.powf(0.045) * we.powf(0.035))
}

/// The void fraction of the drift-flux model of Zuber & Findlay (1965)
/// with the distribution parameter and drift velocity of Rouhani &
/// Axelsson (1970): C₀ = 1 + 0.12 (1 − x), V_gj = 1.18 (1 − x) [σ g (ρ_l −
/// ρ_g)/ρ_l²]^¼. `sin_up` is the sine of the angle of the flow direction
/// above the horizontal; the drift velocity takes its sign and size (no
/// drift in a horizontal tube). With `drift = false` the homogeneous void
/// (C₀ = 1, V_gj = 0).
///
/// Rouhani, S. Z., Axelsson, E. (1970). Calculation of void volume
/// fraction in the subcooled and quality boiling regions. Int. J. Heat
/// Mass Transfer 13, 383–393.
#[allow(clippy::too_many_arguments)]
pub fn void_fraction(
    x: f64,
    g_flux: f64,
    rho_l: f64,
    rho_g: f64,
    sigma: f64,
    sin_up: f64,
    drift: bool,
) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    if !drift {
        return x / (x + rho_g / rho_l * (1.0 - x));
    }
    let c0 = 1.0 + 0.12 * (1.0 - x);
    let vgj =
        1.18 * (1.0 - x) * (sigma * G * (rho_l - rho_g) / (rho_l * rho_l)).powf(0.25) * sin_up;
    let a = x / (c0 * (x + rho_g / rho_l * (1.0 - x)) + rho_g * vgj / g_flux.max(1e-12));
    a.clamp(0.0, 1.0)
}

/// The momentum density ρ′ of a separated flow: 1/ρ′ = x²/(α ρ_g) +
/// (1 − x)²/((1 − α) ρ_l); the flux of momentum is G²/ρ′.
pub fn momentum_density(x: f64, alpha: f64, rho_l: f64, rho_g: f64) -> f64 {
    let x = x.clamp(0.0, 1.0);
    if x <= 0.0 || alpha <= 0.0 {
        return rho_l;
    }
    if x >= 1.0 || alpha >= 1.0 {
        return rho_g;
    }
    1.0 / (x * x / (alpha * rho_g) + (1.0 - x).powi(2) / ((1.0 - alpha) * rho_l))
}

/// The Nusselt number of a fully developed tube flow with a uniform heat
/// flux: Gnielinski (1976) in turbulent flow with the friction factor of
/// Konakov–Filonenko, 48/11 in laminar flow, and the linear blend of
/// Gnielinski (2013) between Re = 2300 and 10⁴.
///
/// Gnielinski, V. (1976). New equations for heat and mass transfer in
/// turbulent pipe and channel flow. Int. Chem. Eng. 16, 359–368.
/// Gnielinski, V. (2013). On heat transfer in tubes. Int. J. Heat Mass
/// Transfer 63, 134–140.
pub fn nusselt_tube(re: f64, pr: f64) -> f64 {
    const NU_LAM: f64 = 48.0 / 11.0;
    let turb = |re: f64| {
        let f = (1.8 * re.log10() - 1.5).powi(-2);
        (f / 8.0) * (re - 1000.0) * pr
            / (1.0 + 12.7 * (f / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0))
    };
    if re <= 2300.0 {
        NU_LAM
    } else if re >= 1e4 {
        turb(re)
    } else {
        let g = (re - 2300.0) / (1e4 - 2300.0);
        (1.0 - g) * NU_LAM + g * turb(1e4)
    }
}

/// The heat transfer coefficient of Dittus & Boelter (1930) in the form
/// of McAdams for heating, 0.023 Re^0.8 Pr^0.4 k/D (the liquid term of
/// Chen and of Gungor–Winterton).
pub fn dittus_boelter(re: f64, pr: f64, k: f64, d: f64) -> f64 {
    0.023 * re.max(1.0).powf(0.8) * pr.powf(0.4) * k / d
}

/// The Martinelli parameter of turbulent liquid and vapour, X_tt.
pub fn x_tt(x: f64, rho_l: f64, rho_g: f64, mu_l: f64, mu_g: f64) -> f64 {
    ((1.0 - x) / x).powf(0.9) * (rho_g / rho_l).sqrt() * (mu_l / mu_g).powf(0.1)
}

/// The inputs of the flow boiling correlations at one segment.
#[derive(Clone, Copy, Debug)]
pub struct BoilingInput {
    /// Wall heat flux into the fluid (W/m²).
    pub q: f64,
    /// Mass flux (kg/(m² s)).
    pub g_flux: f64,
    /// Inner diameter (m).
    pub d: f64,
    /// Equilibrium quality (below 0 in subcooled boiling).
    pub x: f64,
    /// Bulk temperature of the liquid (K), T_sat in saturated boiling.
    pub t_bulk: f64,
    /// Pressure (Pa).
    pub p: f64,
    /// Saturated liquid.
    pub rho_l: f64,
    /// Saturated vapour.
    pub rho_g: f64,
    /// Liquid viscosity (Pa s).
    pub mu_l: f64,
    /// Vapour viscosity (Pa s).
    pub mu_g: f64,
    /// Liquid conductivity (W/(m K)).
    pub k_l: f64,
    /// Liquid isobaric specific heat (J/(kg K)).
    pub cp_l: f64,
    /// Surface tension (N/m).
    pub sigma: f64,
    /// Enthalpy of evaporation (J/kg).
    pub h_fg: f64,
    /// Saturation temperature (K).
    pub t_sat: f64,
}

/// The wall temperature (K) of flow boiling by Chen (1966) with the
/// subcooled extension of Butterworth (1977): q = h_mac (T_w − T_b) +
/// h_mic (T_w − T_sat), solved for T_w by bisection. `p_sat` gives the
/// saturation pressure at a temperature.
///
/// Chen, J. C. (1966). Correlation for boiling heat transfer to saturated
/// fluids in convective flow. Ind. Eng. Chem. Process Des. Dev. 5(3),
/// 322–329. Butterworth, D. (1977), in Two-phase flow and heat transfer,
/// Oxford University Press, chapter 11.
pub fn chen_wall_temperature(b: &BoilingInput, p_sat: impl Fn(f64) -> f64) -> f64 {
    let xs = b.x.clamp(0.0, 1.0 - 1e-9);
    let re_l = b.g_flux * (1.0 - xs) * b.d / b.mu_l;
    let pr_l = b.mu_l * b.cp_l / b.k_l;
    let f = if xs <= 0.0 {
        1.0
    } else {
        let inv = 1.0 / x_tt(xs, b.rho_l, b.rho_g, b.mu_l, b.mu_g);
        if inv <= 0.1 {
            1.0
        } else {
            2.35 * (inv + 0.213).powf(0.736)
        }
    };
    let h_mac = dittus_boelter(re_l, pr_l, b.k_l, b.d) * f;
    let s = 1.0 / (1.0 + 2.53e-6 * (re_l * f.powf(1.25)).powf(1.17));
    let fz = 0.00122 * b.k_l.powf(0.79) * b.cp_l.powf(0.45) * b.rho_l.powf(0.49)
        / (b.sigma.powf(0.5) * b.mu_l.powf(0.29) * b.h_fg.powf(0.24) * b.rho_g.powf(0.24));
    let flux = |tw: f64| {
        let dts = (tw - b.t_sat).max(0.0);
        let dps = (p_sat(tw) - b.p).max(0.0);
        let h_mic = fz * dts.powf(0.24) * dps.powf(0.75) * s;
        h_mac * (tw - b.t_bulk) + h_mic * dts
    };
    // the flux rises with T_w; bracket from the bulk upward
    let mut lo = b.t_bulk.max(b.t_sat.min(b.t_bulk));
    let mut hi = b.t_sat.max(b.t_bulk) + 1.0;
    while flux(hi) < b.q {
        hi += 2.0 * (hi - lo).max(1.0);
        if hi > 2000.0 {
            return hi;
        }
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if flux(mid) < b.q {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-9 * hi {
            break;
        }
    }
    0.5 * (lo + hi)
}

/// The pool boiling coefficient of Cooper (1984) for water (M = 18.015):
/// h = 55 p_r^0.12 (−log₁₀ p_r)^−0.55 M^−0.5 q^0.67 (smooth surface).
///
/// Cooper, M. G. (1984). Heat flow rates in saturated nucleate pool
/// boiling — a wide-ranging examination using reduced properties.
/// Advances in Heat Transfer 16, 157–239.
pub fn cooper(p_reduced: f64, q: f64) -> f64 {
    55.0 * p_reduced.powf(0.12)
        * (-p_reduced.log10()).powf(-0.55)
        * 18.015_f64.powf(-0.5)
        * q.powf(0.67)
}

/// The wall temperature (K) of flow boiling by Gungor & Winterton (1986):
/// saturated h = E h_l + S h_pool; subcooled q = h_l (T_w − T_b) +
/// S h_pool (T_w − T_sat) with E = 1. Horizontal tubes take the Froude
/// corrections of the paper when `horizontal` is set.
///
/// Gungor, K. E., Winterton, R. H. S. (1986). A general correlation for
/// flow boiling in tubes and annuli. Int. J. Heat Mass Transfer 29(3),
/// 351–358.
pub fn gungor_winterton_wall_temperature(b: &BoilingInput, p_crit: f64, horizontal: bool) -> f64 {
    let xs = b.x.clamp(0.0, 1.0 - 1e-9);
    let re_l = b.g_flux * (1.0 - xs) * b.d / b.mu_l;
    let pr_l = b.mu_l * b.cp_l / b.k_l;
    let h_l = dittus_boelter(re_l, pr_l, b.k_l, b.d);
    let h_pool = cooper(b.p / p_crit, b.q);
    let bo = b.q / (b.g_flux * b.h_fg);
    let (mut e, mut s_fac) = (1.0, 1.0);
    if horizontal {
        let fr = b.g_flux * b.g_flux / (b.rho_l * b.rho_l * G * b.d);
        if fr < 0.05 {
            e = fr.powf(0.1 - 2.0 * fr);
            s_fac = fr.sqrt();
        }
    }
    if b.x <= 0.0 {
        let s = s_fac / (1.0 + 1.15e-6 * re_l.powf(1.17));
        (b.q + h_l * b.t_bulk + s * h_pool * b.t_sat) / (h_l + s * h_pool)
    } else {
        let inv = 1.0 / x_tt(xs, b.rho_l, b.rho_g, b.mu_l, b.mu_g);
        let e_tp = e * (1.0 + 24000.0 * bo.powf(1.16) + 1.37 * inv.powf(0.86));
        let s = s_fac / (1.0 + 1.15e-6 * e_tp * e_tp * re_l.powf(1.17));
        b.t_sat + b.q / (e_tp * h_l + s * h_pool)
    }
}

/// The equilibrium quality of the onset of significant void by Saha &
/// Zuber (1974): ΔT_sub = 0.0022 q D/k for Pe ≤ 70 000, 153.8 q/(G c_p)
/// above; x_d = −c_p ΔT_sub/h_fg.
///
/// Saha, P., Zuber, N. (1974). Point of net vapor generation and vapor
/// void fraction in subcooled boiling. Proc. 5th Int. Heat Transfer
/// Conference, Tokyo, B4.7, 175–179.
pub fn saha_zuber_x_d(q: f64, g_flux: f64, d: f64, k_l: f64, cp_l: f64, h_fg: f64) -> f64 {
    let pe = g_flux * d * cp_l / k_l;
    let dt = if pe <= 70_000.0 {
        0.0022 * q * d / k_l
    } else {
        153.8 * q / (g_flux * cp_l)
    };
    -cp_l * dt / h_fg
}

/// The true (flowing) quality of subcooled and low-quality boiling by the
/// profile of Levy (1967): x = x_eq − x_d exp(x_eq/x_d − 1) for x_eq ≥
/// x_d (x_d < 0), 0 below.
///
/// Levy, S. (1967). Forced convection subcooled boiling — prediction of
/// vapor volumetric fraction. Int. J. Heat Mass Transfer 10, 951–965.
pub fn levy_quality(x_eq: f64, x_d: f64) -> f64 {
    if x_d >= 0.0 {
        return x_eq.max(0.0);
    }
    if x_eq <= x_d {
        return 0.0;
    }
    (x_eq - x_d * (x_eq / x_d - 1.0).exp()).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn churchill_is_laminar_and_colebrook_like() {
        for re in [10.0, 100.0, 1000.0] {
            assert!((darcy_churchill(re, 0.0) * re / 64.0 - 1.0).abs() < 1e-12);
        }
        // smooth tube, Re = 1e5: Colebrook 0.01799
        let f = darcy_churchill(1e5, 0.0);
        assert!((f - 0.01799).abs() / 0.01799 < 0.02, "f = {f}");
        // rough: e/D = 1e-3, Re = 1e6, Colebrook 0.01993
        let f = darcy_churchill(1e6, 1e-3);
        assert!((f - 0.01993).abs() / 0.01993 < 0.02, "f = {f}");
    }

    #[test]
    fn gnielinski_is_close_to_dittus_boelter_and_exact_at_pr_one() {
        // at Pr = 1 the denominator is one: Nu = (f/8)(Re - 1000)
        let f = (1.8 * 5.0_f64 - 1.5).powi(-2);
        assert!((nusselt_tube(1e5, 1.0) - f / 8.0 * 99_000.0).abs() < 1e-9);
        for (re, pr) in [(2e4, 0.9), (5e4, 2.0), (1e5, 5.0)] {
            let db = 0.023 * f64::powf(re, 0.8) * f64::powf(pr, 0.4);
            let nu = nusselt_tube(re, pr);
            assert!(
                (nu / db - 1.0).abs() < 0.2, // Gnielinski is up to 17 % above at Pr = 5
                "Re {re} Pr {pr}: {nu} vs {db}"
            );
        }
        assert_eq!(nusselt_tube(1000.0, 5.0), 48.0 / 11.0);
    }

    #[test]
    fn levy_profile_meets_its_ends() {
        let xd = -0.05;
        assert_eq!(levy_quality(-0.06, xd), 0.0);
        assert!(levy_quality(xd, xd).abs() < 1e-15);
        assert!((levy_quality(0.5, xd) - 0.5).abs() < 1e-4);
        // monotone
        let mut last = 0.0;
        for i in 0..100 {
            let x = levy_quality(xd + i as f64 * 0.002, xd);
            assert!(x >= last);
            last = x;
        }
    }

    #[test]
    fn homogeneous_void_and_momentum_density() {
        let a = void_fraction(0.1, 1000.0, 740.0, 36.5, 0.018, 1.0, false);
        let expect = 0.1 / (0.1 + 36.5 / 740.0 * 0.9);
        assert!((a - expect).abs() < 1e-15);
        // with no slip the momentum density is the homogeneous density
        let rho = momentum_density(0.1, a, 740.0, 36.5);
        let rho_h = 1.0 / (0.1 / 36.5 + 0.9 / 740.0);
        assert!((rho - rho_h).abs() / rho_h < 1e-12);
        // drift flux gives less void in upflow
        assert!(void_fraction(0.1, 1000.0, 740.0, 36.5, 0.018, 1.0, true) < a);
    }
}
