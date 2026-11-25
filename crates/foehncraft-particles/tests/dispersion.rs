//! Turbulent dispersion of the parcels by the eddy-lifetime model (SP2c,
//! D243; Gosman and Ioannides 1983): tracers in homogeneous isotropic
//! turbulence at rest spread as the random walk the model is. Each eddy
//! gives a tracer the displacement u' τ_e with u' Gaussian of variance 2k/3
//! per component, the eddies are independent, so after a time T = n τ_e + r
//! the expected squared displacement per component is (2k/3)(n τ_e² + r²)
//! (Taylor 1921 for a velocity that is constant over τ_e and then
//! forgotten). The mean displacement is zero.

use foehncraft_kernels::particle::eddy_times;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, Dispersion, Env, Gas, GasState, Injector, Spray, SpraySection,
};

const K: f64 = 1.0;
const EPS: f64 = 1.0;

/// Gas at rest with uniform turbulence k, ε.
struct Still(f64, f64);

impl Gas for Still {
    fn at(&self, _: usize, _: [f64; 3]) -> GasState {
        let mut g = GasState::isothermal([0.0; 3], 1.2, 1.8e-5);
        g.tke = self.0;
        g.eps = self.1;
        g
    }
}

#[test]
fn tracers_spread_as_the_random_walk_of_the_eddies() {
    let mut s = BoxSpec::new(Kind::Hex, 3, [10, 10, 10]);
    s.lower = [-10.0; 3];
    s.upper = [10.0; 3];
    let m = s.generate().unwrap();
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    // 1 µm water droplets: τ_p ≈ 3 µs against τ_e ≈ 0.2 s, tracers
    let dt = 0.01;
    let cfg = SpraySection {
        injectors: vec![Injector {
            name: "puff".into(),
            position: [0.1, -0.2, 0.3],
            direction: [1.0, 0.0, 0.0],
            speed: 0.0,
            half_angle: 0.0,
            mass_flow: 1e-9,
            diameter: 1e-6,
            size: None,
            parcels_per_second: 2e6,
            start: 0.0,
            end: Some(dt),
            temperature: None,
        }],
        dispersion: Dispersion::EddyLifetime,
        ..SpraySection::default()
    };
    let rules = vec![
        BoundaryRule {
            mode: BoundaryMode::Rebound,
            e_normal: 1.0,
            e_tangential: 1.0,
        };
        m.boundaries.len()
    ];
    let mut sp = Spray::new(&cfg, &env, rules, [0.0; 3]).unwrap();
    let steps = 200;
    for n in 0..steps {
        sp.advance(&env, &Still(K, EPS), n as f64 * dt, dt).unwrap();
    }
    let t_end = steps as f64 * dt;
    let n = sp.parcels.len();
    assert!(n >= 19_000, "{n} parcels");
    let te = eddy_times(K, EPS, 1e-9, 0.0).lifetime;
    let var = 2.0 / 3.0 * K;
    let (mut sum, mut sum2, mut expect) = ([0.0; 3], [0.0; 3], 0.0);
    for i in 0..n {
        let age = t_end - sp.parcels.born[i];
        let eddies = (age / te).floor();
        let r = age - eddies * te;
        expect += var * (eddies * te * te + r * r);
        for k in 0..3 {
            let dx = sp.parcels.x[i][k] - cfg.injectors[0].position[k];
            sum[k] += dx;
            sum2[k] += dx * dx;
        }
        // every parcel is in its cell, far from the walls
        assert!(sp.parcels.x[i].iter().all(|x| x.abs() < 9.0));
    }
    let rms = (expect / n as f64).sqrt();
    for k in 0..3 {
        // the mean: zero within four standard errors
        let mean = sum[k] / n as f64;
        assert!(
            mean.abs() < 4.0 * rms / (n as f64).sqrt(),
            "component {k}: mean {mean}, rms {rms}"
        );
        // the mean square: within four standard errors of a Gaussian
        // sample, (2/n)^½
        let ratio = sum2[k] / expect;
        eprintln!("component {k}: <x²> / theory = {ratio:.4}");
        assert!(
            (ratio - 1.0).abs() < 4.0 * (2.0 / n as f64).sqrt(),
            "component {k}: <x²> / theory = {ratio}"
        );
    }
    assert_eq!(sp.tally.lost_parcels, 0);
}

/// A heavy droplet that falls through the eddies crosses them before they
/// die (crossing trajectories): it changes eddies more often and spreads
/// less sideways than a tracer in the same time.
#[test]
fn falling_droplets_cross_the_eddies_and_spread_less() {
    let mut s = BoxSpec::new(Kind::Hex, 3, [10, 10, 10]);
    s.lower = [-10.0; 3];
    s.upper = [10.0; 3];
    let m = s.generate().unwrap();
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let dt = 0.01;
    let spread = |d: f64| {
        let cfg = SpraySection {
            injectors: vec![Injector {
                name: "puff".into(),
                position: [0.0, 0.0, 8.0],
                direction: [0.0, 0.0, -1.0],
                speed: 0.0,
                half_angle: 0.0,
                mass_flow: 1e-9,
                diameter: d,
                size: None,
                parcels_per_second: 5e5,
                start: 0.0,
                end: Some(dt),
                temperature: None,
            }],
            dispersion: Dispersion::EddyLifetime,
            ..SpraySection::default()
        };
        let rules = vec![
            BoundaryRule {
                mode: BoundaryMode::Escape,
                e_normal: 1.0,
                e_tangential: 1.0,
            };
            m.boundaries.len()
        ];
        let mut sp = Spray::new(&cfg, &env, rules, [0.0, 0.0, -9.81]).unwrap();
        for n in 0..100 {
            // weak turbulence, u' ≈ 0.26 m/s, eddies of 5 cm living 0.2 s
            sp.advance(&env, &Still(0.1, 0.1), n as f64 * dt, dt)
                .unwrap();
        }
        let n = sp.parcels.len() as f64;
        let x2: f64 = sp.parcels.x.iter().map(|x| x[0] * x[0] + x[1] * x[1]).sum();
        let eddies: f64 = sp.parcels.eddy.iter().map(|e| e.n as f64).sum();
        (x2 / n, eddies / n)
    };
    // 1 µm: a tracer; 300 µm: falls at about 1.1 m/s, crosses an eddy in
    // about 50 ms
    let (tracer, e_tracer) = spread(1e-6);
    let (heavy, e_heavy) = spread(300e-6);
    eprintln!(
        "sideways <r²>: tracer {tracer:.4}, heavy {heavy:.4}; eddies {e_tracer:.1} {e_heavy:.1}"
    );
    assert!(e_heavy > 1.5 * e_tracer, "eddies {e_tracer} {e_heavy}");
    assert!(heavy < 0.7 * tracer, "<r²> {tracer} {heavy}");
}
