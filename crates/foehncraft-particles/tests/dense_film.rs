//! The real-fluid gas film of a droplet at high pressure (Spray A, SA2,
//! D313): the dense-fluid transport of Chung et al. (1988) against the
//! NIST reference tables, the Peng-Robinson heat capacity, and the film of
//! an n-dodecane droplet in nitrogen at 6 MPa against the ideal-gas film.

use foehncraft_kernels::realfluid::{MAX, Mixture, Root, fuels::N_DODECANE, gases};
use foehncraft_kernels::species::{R, Species};
use foehncraft_particles::ideal_gas_species;

fn nitrogen() -> Species {
    ideal_gas_species(
        include_str!("../../foehncraft-flow/data/species-gri30.toml"),
        "N2",
    )
    .unwrap()
}

fn pure(c: foehncraft_kernels::realfluid::Component) -> (Mixture, [f64; MAX]) {
    let mut x = [0.0; MAX];
    x[0] = 1.0;
    (Mixture::new(&[c]), x)
}

/// Nitrogen, NIST Chemistry WebBook fluid tables (equation of state of Span
/// et al. 2000, transport of Lemmon and Jacobsen 2004): T (K), p (MPa),
/// ρ (kg/m³), η (Pa s), λ (W/(m K)), c_p (J/(g K)).
const N2: [(f64, f64, f64, f64, f64, f64); 11] = [
    (300.0, 0.1, 1.12328, 1.78899e-05, 0.0259682, 1.04133),
    (600.0, 0.1, 0.561301, 2.95770e-05, 0.0448405, 1.07513),
    (900.0, 0.1, 0.374220, 3.87800e-05, 0.0605199, 1.14573),
    (1200.0, 0.1, 0.280684, 4.67665e-05, 0.0746785, 1.20404),
    (300.0, 6.0, 67.5629, 1.89442e-05, 0.0287804, 1.13602),
    (600.0, 6.0, 32.8266, 2.99556e-05, 0.0460242, 1.09303),
    (900.0, 6.0, 21.9700, 3.89937e-05, 0.0612592, 1.15276),
    (1200.0, 6.0, 16.5453, 4.69097e-05, 0.0752132, 1.20758),
    (300.0, 15.0, 164.239, 2.14870e-05, 0.0342397, 1.25646),
    (600.0, 15.0, 78.7843, 3.07601e-05, 0.0482899, 1.11640),
    (900.0, 15.0, 53.1469, 3.94277e-05, 0.0626566, 1.16229),
];

/// The NIST row at 0.1 MPa of the same temperature.
fn low(t: f64) -> (f64, f64, f64, f64, f64, f64) {
    *N2.iter().find(|r| r.0 == t && r.1 == 0.1).unwrap()
}

#[test]
fn chung_follows_the_nist_nitrogen_tables() {
    // with the NIST density as the input (the method alone): the viscosity
    // within 1.5 % and the conductivity within 7 % from 0.1 to 15 MPa; the
    // dense-fluid factors η/η° and λ/λ° that the film takes within 1.5 %
    // and 5 % of the NIST ratios to 0.1 MPa
    let (mix, x) = pure(gases::N2);
    let fluid = mix.chung(&x);
    let n2 = nitrogen();
    let m = gases::N2.molar_mass;
    let mut worst = [0.0f64; 4];
    for (t, p, rho, eta, lambda, _) in N2 {
        let rm = rho / m;
        let cv = n2.cp(t) * m - R;
        let e = fluid.viscosity(t, rm) / eta - 1.0;
        let l = fluid.conductivity(t, rm, cv) / lambda - 1.0;
        worst[0] = worst[0].max(e.abs());
        worst[1] = worst[1].max(l.abs());
        if p > 0.1 {
            let lo = low(t);
            let f_eta = fluid.viscosity_factor(t, rm) / fluid.viscosity_factor(t, lo.2 / m);
            let f_lambda =
                fluid.conductivity_factor(t, rm, cv) / fluid.conductivity_factor(t, lo.2 / m, cv);
            let fe = f_eta / (eta / lo.3) - 1.0;
            let fl = f_lambda / (lambda / lo.4) - 1.0;
            worst[2] = worst[2].max(fe.abs());
            worst[3] = worst[3].max(fl.abs());
            // the dense gas is more viscous and conducts better
            assert!(f_eta > 1.0 && f_lambda > 1.0, "{t} K {p} MPa");
        }
    }
    // own results: 1.3 %, 6.8 %, 1.3 %, 4.8 % (recorded so a change shows)
    assert!(worst[0] < 0.015, "{worst:?}");
    assert!(worst[1] < 0.07, "{worst:?}");
    assert!(worst[2] < 0.015, "{worst:?}");
    assert!(worst[3] < 0.05, "{worst:?}");
}

#[test]
fn the_factors_vanish_at_low_density_and_riazi_whitson_slows_dense_diffusion() {
    let (mix, x) = pure(gases::N2);
    let fluid = mix.chung(&x);
    let cv = nitrogen().cp(600.0) * gases::N2.molar_mass - R;
    // a gas at 1 Pa: the dilute-gas values
    let rm = 1.0 / (R * 600.0);
    assert!((fluid.viscosity_factor(600.0, rm) - 1.0).abs() < 1e-6);
    assert!((fluid.conductivity_factor(600.0, rm, cv) - 1.0).abs() < 1e-6);
    assert!((mix.riazi_whitson(&x, 1.0, 1.0) - 1.0).abs() < 1e-15);
    // a denser gas has the larger viscosity factor and so the smaller ρD
    let mut last = 1.0;
    for p in [2e6, 6e6, 15e6] {
        let ph = mix.phase(&x, 600.0, p, Root::Vapour);
        let f = mix.riazi_whitson(&x, p, fluid.viscosity_factor(600.0, 1.0 / ph.v));
        assert!(f < last && f > 0.8, "{p}: {f}");
        last = f;
    }
}

#[test]
fn the_peng_robinson_heat_capacity_and_density_of_nitrogen() {
    // c_p = c_p^ig + ∂(h − h^ig)/∂T against NIST within 1 %, the density
    // within 2.5 % (own results: 0.8 % and 2.2 %, the largest at 300 K and
    // 15 MPa)
    let (mix, x) = pure(gases::N2);
    let n2 = nitrogen();
    let m = gases::N2.molar_mass;
    for (t, p, rho, _, _, cp) in N2 {
        let ph = mix.phase(&x, t, p * 1e6, Root::Vapour);
        let c = n2.cp(t) + mix.cp_departure(&x, t, p * 1e6, Root::Vapour) / m;
        assert!((c / (cp * 1e3) - 1.0).abs() < 0.01, "{t} K {p} MPa: {c}");
        assert!((m / ph.v / rho - 1.0).abs() < 0.025, "{t} K {p} MPa");
    }
}

/// n-dodecane, NIST WebBook (Lemmon and Huber 2004; Huber, Laesecke and
/// Perkins 2004): T (K), η and λ at 0.1 MPa (vapour) and ρ (kg/m³), η, λ
/// at 6 MPa.
const DODECANE: [(f64, f64, f64, f64, f64, f64, f64); 3] = [
    (
        500.0,
        4.38359,
        7.59345e-06,
        0.0256701,
        600.339,
        0.000205248,
        0.0979494,
    ),
    (
        600.0,
        3.52778,
        9.27663e-06,
        0.0355249,
        508.188,
        0.000104313,
        0.0824061,
    ),
    (
        700.0,
        2.97980,
        1.09291e-05,
        0.0460392,
        387.424,
        4.83236e-05,
        0.0723930,
    ),
];

#[test]
fn chung_overestimates_the_liquid_like_dodecane() {
    // the gap of the method (D313): for n-dodecane at 6 MPa and liquid-like
    // densities the factor on the viscosity is 9 to 17 % and the one on the
    // conductivity 40 to 55 % above NIST (own results); the film of a
    // Spray A droplet is mostly nitrogen gas
    let (mix, x) = pure(N_DODECANE.component);
    let fluid = mix.chung(&x);
    let m = N_DODECANE.component.molar_mass;
    let ig = N_DODECANE.ideal_gas;
    for (t, rho0, eta0, l0, rho, eta, lambda) in DODECANE {
        let cv = ig.cp(t) * m - R;
        let fe = fluid.viscosity_factor(t, rho / m) / fluid.viscosity_factor(t, rho0 / m);
        let fl =
            fluid.conductivity_factor(t, rho / m, cv) / fluid.conductivity_factor(t, rho0 / m, cv);
        let ee = fe / (eta / eta0) - 1.0;
        let el = fl / (lambda / l0) - 1.0;
        assert!(ee > 0.0 && ee < 0.2, "{t} K: {ee}");
        assert!(el > 0.3 && el < 0.6, "{t} K: {el}");
    }
}

mod film {
    use super::*;
    use foehncraft_mesh::Geometry;
    use foehncraft_mesh::generate::{BoxSpec, Kind};
    use foehncraft_mesh::locate::Locator;
    use foehncraft_particles::{
        BoundaryMode, BoundaryRule, Env, Equilibrium, EvaporationModel, FilmFluid, GasState,
        Injector, Spray, SpraySection, Vapour,
    };

    fn gas(t: f64, p0: f64) -> GasState {
        GasState {
            u: [0.0; 3],
            rho: p0 * 0.028 / (R * t),
            mu: 1.85e-5,
            t,
            cp: 1040.0,
            k: 0.03,
            tke: 0.0,
            eps: 0.0,
            y_v: 0.0,
            p0,
            m_gas: 0.028_013_4,
            diff: 2.5e-5,
        }
    }

    /// A spray of 363 K n-dodecane with the Peng-Robinson surface in
    /// nitrogen at `p0` and the film model `film`.
    fn spray(p0: f64, film: FilmFluid) -> Spray {
        let mut s = BoxSpec::new(Kind::Hex, 3, [1, 1, 1]);
        s.lower = [0.0; 3];
        s.upper = [0.01; 3];
        let m = s.generate().unwrap();
        let geo = Geometry::compute(&m);
        let loc = Locator::new(&m, &geo);
        let env = Env {
            mesh: &m,
            geo: &geo,
            loc: &loc,
        };
        let mut cfg = SpraySection {
            injectors: vec![Injector {
                name: "n".into(),
                position: [0.005; 3],
                direction: [1.0, 0.0, 0.0],
                speed: 1e-4,
                mass_flow: 1e-6,
                diameter: 20e-6,
                parcels_per_second: 1e4,
                end: Some(1e-3),
                temperature: Some(363.0),
                ..Injector::default()
            }],
            ..SpraySection::default()
        };
        cfg.liquid.library = Some("n-dodecane".into());
        cfg.liquid.temperature = 363.0;
        cfg.liquid.vapour = Some("NC12H26".into());
        cfg.liquid.apply_library().unwrap();
        cfg.evaporation.model = EvaporationModel::AbramzonSirignano;
        cfg.evaporation.equilibrium = Equilibrium::PengRobinson;
        cfg.evaporation.peng_robinson.film = film;
        let rules = vec![
            BoundaryRule {
                mode: BoundaryMode::Rebound,
                e_normal: 0.5,
                e_tangential: 0.5,
            };
            m.boundaries.len()
        ];
        let mut sp = Spray::new(&cfg, &env, rules, [0.0; 3]).unwrap();
        sp.set_vapour(Some(Vapour {
            species: N_DODECANE.ideal_gas,
        }));
        sp.set_carrier(vec![(nitrogen(), 1.0)]);
        sp.set_surface_pressure(p0).unwrap();
        sp
    }

    #[test]
    fn a_gas_above_the_critical_mixing_state_takes_vapour() {
        // the saturation limiter of D253 with the Peng-Robinson surface:
        // above the critical mixing temperature (651 K at 6 MPa) the
        // surface takes no latent heat, so the gas takes vapour without
        // cooling, up to the fuel fraction of the critical surface. The
        // capacity was zero there, and in a Spray A run the limiter held
        // back all the evaporation in the hot gas (35 times what still
        // evaporated, droplets of 1 µm heating to 900 K)
        let sp = spray(6e6, FilmFluid::RealFluid);
        let g = gas(900.0, 6e6);
        let s = sp.surface_state(g.t).unwrap();
        assert!(!s.two_phase);
        assert_eq!(sp.latent_heat(g.t), 0.0);
        let ys = foehncraft_kernels::evaporation::surface_fraction(
            sp.saturation_pressure(g.t),
            g.p0,
            N_DODECANE.ideal_gas.molar_mass,
            g.m_gas,
        );
        let m = 1e-9;
        let cap = sp.saturation_capacity(&g, m);
        assert!(ys > 0.3 && ys < 0.999, "{ys}");
        assert!((cap / (m * ys / (1.0 - ys)) - 1.0).abs() < 1e-12, "{cap}");
        // below it the gas cools as it takes the vapour: less, but some
        let cold = sp.saturation_capacity(&gas(600.0, 6e6), m);
        assert!(cold > 0.0 && cold < cap, "{cold} {cap}");
        // a gas at the critical fraction takes none
        let full = foehncraft_particles::GasState { y_v: ys, ..g };
        assert_eq!(sp.saturation_capacity(&full, m), 0.0);
    }

    #[test]
    fn at_low_pressure_the_real_fluid_film_is_the_ideal_film() {
        let real = spray(1e5, FilmFluid::RealFluid);
        let ideal = spray(1e5, FilmFluid::IdealGas);
        let g = gas(600.0, 1e5);
        for (tp, ys) in [(400.0, 0.05), (450.0, 0.3)] {
            let a = real.film_properties(&g, tp, ys);
            let b = ideal.film_properties(&g, tp, ys);
            for k in 0..6 {
                assert!((a[k] / b[k] - 1.0).abs() < 0.01, "{tp} K, {k}: {a:?} {b:?}");
            }
        }
    }

    #[test]
    fn the_film_of_a_spray_a_droplet() {
        // 900 K nitrogen at 6 MPa, droplets at 450 to 600 K with the fuel
        // fraction of the Peng-Robinson surface: at the film temperature
        // (600 to 700 K) the film is a dense gas only 2 % above the ideal
        // (Z ≈ 1.02), so the real fluid moves its properties by 1 to 4 %
        // (own results: ρ 0.98, μ 1.01 to 1.03, k 1.02 to 1.04, c_p 1.02,
        // D 1.01 to 1.02 times the ideal film)
        let real = spray(6e6, FilmFluid::RealFluid);
        let ideal = spray(6e6, FilmFluid::IdealGas);
        let g = gas(900.0, 6e6);
        let mf = N_DODECANE.component.molar_mass;
        let bounds = [
            (0.97, 0.99),
            (1.005, 1.04),
            (1.01, 1.05),
            (1.01, 1.03),
            (1.0, 1.03),
        ];
        for tp in [450.0, 550.0, 600.0] {
            let s = real.surface_state(tp).unwrap();
            assert!(s.two_phase, "{tp} K");
            // the mass fraction of the fuel at the surface
            let ys = s.y_fuel * mf / (s.y_fuel * mf + (1.0 - s.y_fuel) * 0.028_013_4);
            let a = real.film_properties(&g, tp, ys);
            let b = ideal.film_properties(&g, tp, ys);
            for (k, (lo, hi)) in bounds.into_iter().enumerate() {
                let r = a[k] / b[k];
                assert!(r > lo && r < hi, "{tp} K, Y_s {ys:.3}, property {k}: {r}");
            }
            // the vapour c_p stays the ideal gas's (the enthalpy books)
            assert_eq!(a[5].to_bits(), b[5].to_bits());
        }
    }
}
