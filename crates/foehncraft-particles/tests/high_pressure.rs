//! The high-pressure phase equilibrium of n-dodecane and nitrogen at the
//! droplet surface (ECN Spray A, stage SA1, checks VA.1 to VA.3 of
//! `docs/requirements/step28-spray-a.md`). Every measured number is quoted
//! here with its source (`docs/verification/DATA_SOURCES.md`, section 3e).

use foehncraft_kernels::realfluid::{self as rf, MAX, Mixture, Root, fuels::N_DODECANE, gases};
use foehncraft_particles::{SurfaceTable, TableEnd, fit_mathias_copeman};

/// The vapour pressure of n-dodecane (Pa) by eq. (3) and table 1 of
/// Ambrose, D., Walton, J. (1989), Vapour pressures up to their critical
/// temperatures of normal alkanes and 1-alkanols, Pure Appl. Chem. 61,
/// 1395-1403, doi 10.1351/pac198961081395: ln(p/p_c) = (a τ + b τ^1.5 +
/// c τ^2.5 + d τ^5)/T_r, τ = 1 − T_r, T_c = 658.65 K, p_c = 1.83 MPa.
fn ambrose_walton(t: f64) -> f64 {
    let (tc, pc) = (658.65, 1.83e6);
    let (a, b, c, d) = (-9.08593, 2.77846, -5.1985, -4.173);
    let tr = t / tc;
    let s = 1.0 - tr;
    pc * ((a * s + b * s.powf(1.5) + c * s.powf(2.5) + d * s.powi(5)) / tr).exp()
}

fn dodecane_n2() -> Mixture {
    let mut m = Mixture::new(&[N_DODECANE.component, gases::N2]);
    m.set_k(0, 1, N_DODECANE.k_with("N2"));
    m
}

fn nitrogen_ratios() -> [f64; MAX] {
    let mut r = [0.0; MAX];
    r[1] = 1.0;
    r
}

#[test]
fn the_wagner_equation_gives_the_normal_boiling_point() {
    // the normal boiling point of n-dodecane, 489.44 K (Lemmon, E. W.,
    // Goodwin, A. R. H. (2000), J. Phys. Chem. Ref. Data 29, 1-39, doi
    // 10.1063/1.556054, table 2): a check of the transcribed coefficients
    let p = ambrose_walton(489.439);
    assert!((p / 101_325.0 - 1.0).abs() < 2e-3, "{p}");
}

#[test]
fn va1_the_peng_robinson_vapour_pressure_follows_the_measured_equation() {
    let fuel = Mixture::new(&[N_DODECANE.component]);
    let mut plain = N_DODECANE.component;
    plain.mc = [0.0; 3];
    let plain = Mixture::new(&[plain]);
    let mut worst: f64 = 0.0;
    let mut t = 300.0;
    while t <= 655.0 {
        let p = fuel.saturation_pressure(0, t).expect("two phases below Tc");
        worst = worst.max((p / ambrose_walton(t) - 1.0).abs());
        t += 5.0;
    }
    assert!(worst < 0.01, "Mathias-Copeman: {worst}");
    // without the fitted α the cold end is far off: 16 % at 300 K, where
    // the injected fuel (363 K in Spray A) starts
    let p = plain.saturation_pressure(0, 300.0).unwrap();
    assert!(
        p / ambrose_walton(300.0) - 1.0 > 0.1,
        "{}",
        p / ambrose_walton(300.0)
    );
}

#[test]
fn fit_mathias_copeman_reproduces_the_dodecane_coefficients() {
    let ts: Vec<f64> = (0..72).map(|i| 300.0 + 5.0 * i as f64).collect();
    let mut c = N_DODECANE.component;
    c.mc = [0.0; 3];
    let (mc, worst) = fit_mathias_copeman(c, &ambrose_walton, &ts).unwrap();
    for (got, want) in mc.iter().zip(N_DODECANE.component.mc.iter()) {
        assert!((got - want).abs() < 2e-3, "{mc:?}");
    }
    assert!(worst < 0.01, "{worst}");
}

/// Bubble points of N₂ + n-dodecane, a subset (three per isotherm, near 2,
/// 5 and 12 MPa) of table 2 of García-Córdova, T., Justo-García, D. N.,
/// García-Flores, B. E., García-Sánchez, F. (2011), Vapor-liquid
/// equilibrium data for the nitrogen + dodecane system at temperatures
/// from (344 to 593) K and at pressures up to 60 MPa, J. Chem. Eng. Data 56,
/// 1555-1564, doi 10.1021/je1012372: (T in K, p in kPa, x_N2).
const BUBBLE: [(f64, f64, f64); 21] = [
    (344.4, 1430.0, 0.0204),
    (344.4, 5030.0, 0.0669),
    (344.4, 12380.0, 0.1435),
    (410.7, 2170.0, 0.0339),
    (410.7, 5020.0, 0.0754),
    (410.7, 12980.0, 0.1698),
    (463.9, 2050.0, 0.0372),
    (463.9, 5040.0, 0.0875),
    (463.9, 12500.0, 0.1945),
    (503.4, 2160.0, 0.0436),
    (503.4, 5310.0, 0.1041),
    (503.4, 12440.0, 0.2213),
    (532.9, 2030.0, 0.0423),
    (532.9, 6010.0, 0.1271),
    (532.9, 12450.0, 0.2462),
    (562.1, 2130.0, 0.0435),
    (562.1, 5960.0, 0.1315),
    (562.1, 11970.0, 0.2582),
    (593.5, 2270.0, 0.0504),
    (593.5, 5070.0, 0.1308),
    (593.5, 10960.0, 0.2837),
];

#[test]
fn va2_bubble_pressures_of_nitrogen_in_dodecane() {
    let mix = dodecane_n2();
    let mut sum = 0.0;
    let mut worst: f64 = 0.0;
    let mut hot = 0.0;
    for &(t, p_kpa, x_n2) in &BUBBLE {
        let mut x = [0.0; MAX];
        x[0] = 1.0 - x_n2;
        x[1] = x_n2;
        let (p, y) = rf::bubble_pressure(&mix, t, &x, p_kpa * 1e3).expect("bubble point");
        let e = p / (p_kpa * 1e3) - 1.0;
        sum += e.abs();
        worst = worst.max(e.abs());
        if t > 550.0 {
            hot += e / 6.0;
        }
        // the vapour is mostly nitrogen (84 % at 533 K and 2 MPa)
        assert!(y[1] > 0.5 && y[1] < 1.0, "{t} {p_kpa} {y:?}");
    }
    let mean = sum / BUBBLE.len() as f64;
    // k_ij = 0.205 fitted to all 75 points up to 20 MPa: mean 2.7 %, at
    // most 9.1 %; above 550 K Peng-Robinson stays 3 to 9 % low at one
    // k_ij (the measured Henry constant falls faster with T than the
    // equation's), recorded as the model limit
    assert!(mean < 0.04, "mean {mean}");
    assert!(worst < 0.10, "worst {worst}");
    assert!(hot < 0.0 && hot > -0.09, "bias above 550 K {hot}");
}

#[test]
fn va3_low_pressure_limit_is_raoult() {
    // the surface vapour fraction over Raoult's law (the Peng-Robinson
    // saturation pressure over p) tends to one as p → 0, nearly linearly
    // (the Poynting factor and the second virial coefficients of the
    // vapour): 0.3 % at 10 kPa, 2 to 3 % at 0.1 MPa, growing with p
    let mix = dodecane_n2();
    let fuel = Mixture::new(&[N_DODECANE.component]);
    let excess = |t: f64, p: f64| {
        let rf::SurfaceState::TwoPhase(s) = rf::surface(&mix, t, p, &nitrogen_ratios(), None)
        else {
            panic!("two phases at {t} K and {p} Pa");
        };
        s.y[0] / (fuel.saturation_pressure(0, t).unwrap() / p) - 1.0
    };
    for t in [300.0, 350.0, 400.0] {
        let (e1, e3, e10) = (excess(t, 1e4), excess(t, 3e4), excess(t, 1e5));
        assert!(e1 > 0.0 && e1 < 0.005, "{t} K: {e1}");
        assert!(e10 < 0.03, "{t} K: {e10}");
        assert!(e1 < e3 && e3 < e10, "{t} K: {e1} {e3} {e10}");
    }
}

#[test]
fn the_surface_at_six_megapascal_is_richer_than_raoult_and_ends_critical() {
    let mix = dodecane_n2();
    let tab = SurfaceTable::build(&mix, &nitrogen_ratios(), 6e6, 363.0, 700.0, 1.0).unwrap();
    assert_eq!(
        tab.end,
        TableEnd::Critical,
        "{:?} at {}",
        tab.end,
        tab.t_end
    );
    // above the fuel's critical pressure the two-phase surface ends at the
    // mixture critical temperature, between 640 K and the fuel's T_c
    assert!(
        tab.t_end > 645.0 && tab.t_end < N_DODECANE.component.tc,
        "{}",
        tab.t_end
    );
    for t in [450.0, 500.0, 550.0, 600.0] {
        let s = tab.at(t);
        let raoult = ambrose_walton(t) / 6e6;
        let f = s.y_fuel / raoult;
        assert!(s.two_phase && f > 1.5 && f < 1.9, "{t} K: {f}");
        // 10 to 20 mol % nitrogen dissolves in the surface liquid
        assert!(
            s.dissolved > 0.08 && s.dissolved < 0.2,
            "{t} K: {}",
            s.dissolved
        );
        // the latent heat is 10 to 15 % below its saturation value (most
        // of it the departure enthalpy of the fuel in the dense vapour)
        assert!(
            s.latent_ratio > 0.8 && s.latent_ratio < 0.92,
            "{t} K: {}",
            s.latent_ratio
        );
    }
    // and towards the critical end it goes to zero, monotonically
    let mut last = tab.at(550.0).latent_ratio;
    let mut t = 560.0;
    while t < tab.t_end {
        let l = tab.at(t).latent_ratio;
        assert!(l < last, "{t} K: {l} after {last}");
        last = l;
        t += 10.0;
    }
    let near = tab.at(tab.t_end - 0.1);
    assert!(near.latent_ratio < 0.15, "{}", near.latent_ratio);
    let past = tab.at(tab.t_end + 1.0);
    assert!(!past.two_phase && past.latent_ratio == 0.0);
}

#[test]
fn the_critical_mixing_temperature_falls_with_pressure() {
    // above the fuel's critical pressure the end of the two-phase surface
    // follows the critical locus of the binary: 656.6 K at 3 MPa down to
    // 647.8 K at 8 MPa with these constants
    let mix = dodecane_n2();
    let mut last = N_DODECANE.component.tc;
    for p in [3e6, 4e6, 6e6, 8e6] {
        let tab = SurfaceTable::build(&mix, &nitrogen_ratios(), p, 363.0, 700.0, 1.0).unwrap();
        assert_eq!(tab.end, TableEnd::Critical);
        assert!(tab.t_end < last, "{p} Pa: {} after {last}", tab.t_end);
        last = tab.t_end;
    }
}

#[test]
fn below_the_critical_pressure_the_table_ends_boiling() {
    let mix = dodecane_n2();
    let tab = SurfaceTable::build(&mix, &nitrogen_ratios(), 1e6, 363.0, 700.0, 2.0).unwrap();
    assert_eq!(tab.end, TableEnd::Boiling, "{:?} at {}", tab.end, tab.t_end);
    // the boiling point at 1 MPa is above the pure fuel's (the dissolved
    // nitrogen is the lighter component, so it lowers it a little: the
    // two are within a few kelvin)
    let fuel = Mixture::new(&[N_DODECANE.component]);
    let (mut lo, mut hi) = (450.0, 650.0);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if fuel.saturation_pressure(0, mid).unwrap() < 1e6 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    assert!(
        (tab.t_end - lo).abs() < 5.0,
        "table end {} pure boiling {lo}",
        tab.t_end
    );
}

#[test]
fn the_table_is_a_pure_function_of_its_inputs() {
    let mix = dodecane_n2();
    let a = SurfaceTable::build(&mix, &nitrogen_ratios(), 6e6, 363.0, 700.0, 1.0).unwrap();
    let b = SurfaceTable::build(&mix, &nitrogen_ratios(), 6e6, 363.0, 700.0, 1.0).unwrap();
    assert_eq!(a, b);
    for i in 0..a.t.len() {
        assert_eq!(a.y_fuel[i].to_bits(), b.y_fuel[i].to_bits());
        assert_eq!(a.latent_ratio[i].to_bits(), b.latent_ratio[i].to_bits());
    }
}

#[test]
fn the_liquid_is_denser_than_the_vapour_at_the_surface() {
    let mix = dodecane_n2();
    let rf::SurfaceState::TwoPhase(s) = rf::surface(&mix, 550.0, 6e6, &nitrogen_ratios(), None)
    else {
        panic!("two phases");
    };
    let l = mix.phase(&s.x, 550.0, 6e6, Root::Liquid);
    let v = mix.phase(&s.y, 550.0, 6e6, Root::Vapour);
    assert!(l.z < 0.5 * v.z, "{} {}", l.z, v.z);
}
