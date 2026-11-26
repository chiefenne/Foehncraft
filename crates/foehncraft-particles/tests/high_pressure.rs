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
fn the_table_reaches_its_end_from_any_cold_start() {
    // the spray starts the table 60 K below the injected liquid, at most at
    // 250 K: from every start the table runs to the critical mixing state
    // (a stalled iteration near 267 K at 6 MPa once ended it there, as if
    // critical)
    let mix = dodecane_n2();
    let r = nitrogen_ratios();
    for p in [4e6, 6e6, 8e6] {
        let mut end: Option<f64> = None;
        for lo in [150.0, 200.0, 240.0, 250.0] {
            let tab = SurfaceTable::build(&mix, &r, p, lo, 700.0, 1.0).unwrap();
            assert_eq!(tab.end, TableEnd::Critical, "{p} Pa from {lo} K");
            assert!(tab.t_end > 600.0, "{p} Pa from {lo} K: {}", tab.t_end);
            if let Some(e) = end {
                assert!((tab.t_end - e).abs() < 1e-6, "{p} Pa from {lo} K");
            }
            end = Some(tab.t_end);
        }
    }
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

// ---- the surface in the spray (`[spray.evaporation] equilibrium`) ----

mod in_the_spray {
    use foehncraft_mesh::Geometry;
    use foehncraft_mesh::generate::{BoxSpec, Kind};
    use foehncraft_mesh::locate::Locator;
    use foehncraft_particles::{
        BoundaryMode, BoundaryRule, Env, Equilibrium, EvaporationModel, Gas, GasState, Injector,
        Spray, SpraySection, Vapour,
    };

    struct Cell {
        t: f64,
        p0: f64,
    }

    impl Gas for Cell {
        fn at(&self, _: usize, _: [f64; 3]) -> GasState {
            GasState {
                u: [0.0; 3],
                rho: self.p0 * 0.028 / (8.314_462_618 * self.t),
                mu: 1.85e-5,
                t: self.t,
                cp: 1040.0,
                k: 0.03,
                tke: 0.0,
                eps: 0.0,
                y_v: 0.0,
                p0: self.p0,
                m_gas: 0.028_013_4,
                diff: 2.5e-5,
            }
        }
    }

    fn one_cell() -> (foehncraft_mesh::Mesh, Geometry) {
        let mut s = BoxSpec::new(Kind::Hex, 3, [1, 1, 1]);
        s.lower = [0.0; 3];
        s.upper = [0.01; 3];
        let m = s.generate().unwrap();
        let g = Geometry::compute(&m);
        (m, g)
    }

    fn section(name: &str, temp: f64) -> SpraySection {
        let mut cfg = SpraySection {
            injectors: vec![Injector {
                name: "n".into(),
                position: [0.005; 3],
                direction: [1.0, 0.0, 0.0],
                speed: 1e-4,
                half_angle: 0.0,
                mass_flow: 1e-6,
                diameter: 20e-6,
                size: None,
                parcels_per_second: 1e4,
                start: 0.0,
                end: Some(1e-3),
                temperature: Some(temp),
                ..Injector::default()
            }],
            ..SpraySection::default()
        };
        cfg.liquid.name = name.into();
        cfg.liquid.vapour = Some(if name == "water" { "H2O" } else { "NC12H26" }.into());
        cfg.evaporation.model = EvaporationModel::RanzMarshall;
        cfg.evaporation.equilibrium = Equilibrium::PengRobinson;
        cfg
    }

    fn rules(env: &Env) -> Vec<BoundaryRule> {
        vec![
            BoundaryRule {
                mode: BoundaryMode::Rebound,
                e_normal: 0.5,
                e_tangential: 0.5,
            };
            env.mesh.boundaries.len()
        ]
    }

    #[test]
    fn the_spray_takes_the_surface_of_the_table_and_closes_its_books() {
        let (m, geo) = one_cell();
        let loc = Locator::new(&m, &geo);
        let env = Env {
            mesh: &m,
            geo: &geo,
            loc: &loc,
        };
        // water with its critical point (Wagner and Pruss 2002) and the
        // acentric factor of DATA_SOURCES 3e, in nitrogen at 0.1 MPa
        let mut cfg = section("water", 300.0);
        let pr = &mut cfg.evaporation.peng_robinson;
        pr.critical_temperature = Some(647.096);
        pr.critical_pressure = Some(22.064e6);
        pr.acentric_factor = Some(0.3443);
        pr.molar_mass = Some(0.018_015_28);
        let mut sp = Spray::new(&cfg, &env, rules(&env), [0.0; 3]).unwrap();
        sp.set_vapour(Some(Vapour::water()));
        // without the table it is Raoult's law and the library latent heat
        let raoult = sp.saturation_pressure(320.0);
        let low = sp.latent_heat(320.0);
        sp.set_surface_pressure(101_325.0).unwrap();
        let tab = sp.surface_table().expect("table");
        for t in [290.0, 320.0, 350.0] {
            let s = sp.surface_state(t).unwrap();
            assert!(s.two_phase);
            assert_eq!(
                sp.saturation_pressure(t).to_bits(),
                (s.y_fuel * tab.p0).to_bits()
            );
            // the gas receives the vapour at c_l T + L: the books close by
            // construction (the same function on both sides)
            let c_l = sp.config().liquid.specific_heat;
            assert_eq!(
                sp.vapour_enthalpy(t).to_bits(),
                (c_l * t + sp.latent_heat(t)).to_bits()
            );
            // at 0.1 MPa the latent heat stays within 2 % of the library's
            assert!(
                s.latent_ratio > 0.98 && s.latent_ratio <= 1.0,
                "{t} K: {}",
                s.latent_ratio
            );
        }
        assert!((sp.latent_heat(320.0) / low - 1.0).abs() < 0.02);
        // the Peng-Robinson vapour pressure of water with the plain κ(ω)
        // (no fitted α for this polar liquid) is 11 % below the IAPWS
        // equation at 320 K: why the built-in fuels carry a fitted α (VA.1)
        let e = sp.saturation_pressure(320.0) / raoult - 1.0;
        assert!(e < -0.05 && e > -0.15, "{e}");
        // a droplet evaporates with it
        let gas = Cell {
            t: 400.0,
            p0: 101_325.0,
        };
        sp.advance(&env, &gas, 0.0, 1e-3).unwrap();
        let m0 = sp.tally.injected_mass;
        sp.advance(&env, &gas, 1e-3, 1e-3).unwrap();
        assert!(sp.tally.evaporated_mass > 0.0 && sp.tally.evaporated_mass < m0);
        assert_eq!(sp.tally.critical_parcels, 0);
    }

    #[test]
    fn parcels_above_the_critical_mixing_temperature_are_counted() {
        let (m, geo) = one_cell();
        let loc = Locator::new(&m, &geo);
        let env = Env {
            mesh: &m,
            geo: &geo,
            loc: &loc,
        };
        // n-dodecane injected at 655 K into nitrogen at 6 MPa: above the
        // critical mixing temperature 651.3 K of its surface (D314); no
        // vapour species, so the droplets keep their mass and temperature
        let cfg = section("n-dodecane", 655.0);
        let mut sp = Spray::new(&cfg, &env, rules(&env), [0.0; 3]).unwrap();
        sp.set_surface_pressure(6e6).unwrap();
        assert!(sp.surface_table().unwrap().t_end < 655.0);
        let gas = Cell { t: 655.0, p0: 6e6 };
        sp.advance(&env, &gas, 0.0, 1e-3).unwrap();
        let n = sp.parcels.len() as u64;
        assert!(n > 0);
        assert_eq!(sp.tally.critical_parcels, n);
        assert!((sp.tally.critical_mass / sp.mass_in_flow() - 1.0).abs() < 1e-12);
        assert_eq!(sp.latent_heat(655.0), 0.0);
    }

    #[test]
    fn the_keys_are_checked() {
        let (m, geo) = one_cell();
        let loc = Locator::new(&m, &geo);
        let env = Env {
            mesh: &m,
            geo: &geo,
            loc: &loc,
        };
        let fails = |cfg: SpraySection, what: &str| {
            let mut sp = Spray::new(&cfg, &env, rules(&env), [0.0; 3]).unwrap();
            let e = format!("{:#}", sp.set_surface_pressure(6e6).unwrap_err());
            assert!(e.contains(what), "{e}");
        };
        fails(section("kerosene", 363.0), "no high-pressure data");
        let mut c = section("water", 300.0);
        c.evaporation.peng_robinson.critical_temperature = Some(647.0);
        fails(c, "all four");
        let mut c = section("n-dodecane", 363.0);
        c.evaporation.peng_robinson.ambient = [("HE".to_string(), 1.0)].into();
        fails(c, "unknown");
        let mut c = section("n-dodecane", 363.0);
        c.evaporation.peng_robinson.k_ij = [("O2".to_string(), 0.1)].into();
        fails(c, "names no gas");
        // Raoult: no table
        let mut c = section("n-dodecane", 363.0);
        c.evaporation.equilibrium = Equilibrium::Raoult;
        let mut sp = Spray::new(&c, &env, rules(&env), [0.0; 3]).unwrap();
        sp.set_surface_pressure(6e6).unwrap();
        assert!(sp.surface_table().is_none());
    }
}
