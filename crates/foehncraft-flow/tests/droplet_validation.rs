//! Single-droplet evaporation against measurements (D320, D321): n-heptane
//! droplets in still nitrogen at 0.1 MPa, Chauveau, Birouk, Halter and
//! Gökalp (2019), the microgravity rates K0(T) of their eq. (4). The droplet
//! runs through `Spray::evaporation`, the path of the spray solver, with the
//! liquid library, the species data of the flow solver and both film states.
//!
//! The model over-predicts these data (a recorded gap, not a pass): with the
//! far-field properties by 80-90 %, with the 1/3 rule by 40-60 %. The bands
//! below hold the present numbers so that a change shows; the 1/3 rule must
//! stay the closer one.

use foehncraft_flow::species::database;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, Env, EvaporationModel, Film, GasState, Spray, SpraySection, Vapour,
};

const P0: f64 = 101_325.0;

/// The spray of one still droplet of n-heptane in N2 with the film `film`.
fn spray(env: &Env, film: Film) -> Spray {
    let mut cfg = SpraySection::default();
    cfg.liquid.vapour = Some("NC7H16".into());
    cfg.liquid.library = Some("n-heptane".into());
    cfg.liquid.apply_library().unwrap();
    cfg.evaporation.model = EvaporationModel::AbramzonSirignano;
    cfg.evaporation.film = film;
    let rules = vec![
        BoundaryRule {
            mode: BoundaryMode::Rebound,
            e_normal: 0.5,
            e_tangential: 0.5,
        };
        env.mesh.boundaries.len()
    ];
    let mut sp = Spray::new(&cfg, env, rules, [0.0; 3]).unwrap();
    let db = database();
    sp.set_vapour(Some(Vapour {
        species: db["NC7H16"].kernel(),
    }));
    sp.set_carrier(vec![(db["N2"].kernel(), 1.0)]);
    sp
}

/// K (mm²/s, the slope of d²/d0² between 0.8 and 0.3) and the wet-bulb
/// temperature of a droplet of d0 at 300 K in still N2 at `tg`.
fn rate(sp: &Spray, tg: f64, d0: f64) -> (f64, f64) {
    let db = database();
    let n2 = db["N2"].kernel();
    let mu = n2.viscosity(tg);
    // the far-field gas: N2 at tg (the film = "far-field" takes these)
    let gs = GasState {
        u: [0.0; 3],
        rho: P0 * n2.molar_mass / (8.314_462_618 * tg),
        mu,
        t: tg,
        cp: n2.cp(tg),
        k: n2.conductivity(tg, mu),
        tke: 0.0,
        eps: 0.0,
        y_v: 0.0,
        p0: P0,
        m_gas: n2.molar_mass,
        diff: foehncraft_kernels::species::binary_diffusion(&db["NC7H16"].kernel(), &n2, tg, P0),
    };
    let (mut d, mut tp, mut t) = (d0, 300.0, 0.0);
    let mut hist = Vec::new();
    while d > 0.05 * d0 {
        let dt = 1e-4 * (d / d0).powi(2) + 1e-7;
        let s = sp.evaporation(&gs, d, tp, 0.0, dt).unwrap();
        d = s.d;
        tp = s.tp;
        t += dt;
        hist.push((t, (d / d0).powi(2), tp));
        if d == 0.0 {
            break;
        }
    }
    let a = hist.iter().find(|h| h.1 < 0.8).unwrap();
    let b = hist.iter().find(|h| h.1 < 0.3).unwrap();
    ((a.1 - b.1) / (b.0 - a.0) * d0 * d0 * 1e6, b.2)
}

#[test]
fn n_heptane_droplets_against_chauveau_2019() {
    let mut s = BoxSpec::new(Kind::Hex, 3, [1, 1, 1]);
    s.lower = [0.0; 3];
    s.upper = [1.0; 3];
    let mesh = s.generate().unwrap();
    let geo = Geometry::compute(&mesh);
    let loc = Locator::new(&mesh, &geo);
    let env = Env {
        mesh: &mesh,
        geo: &geo,
        loc: &loc,
    };
    let third = spray(&env, Film::OneThird);
    let far = spray(&env, Film::FarField);
    assert!(third.film_rule() && !far.film_rule());
    for tg in [473.0, 548.0, 623.0, 698.0, 748.0] {
        // mm²/s, eq. (4): microgravity, 14 µm cross fibres, d0 ≈ 0.37 mm
        let k0 = 3.6552e-4 * tg - 0.1078;
        let (kt, tw) = rate(&third, tg, 370e-6);
        let (kf, _) = rate(&far, tg, 370e-6);
        let (et, ef) = (kt / k0 - 1.0, kf / k0 - 1.0);
        eprintln!(
            "T {tg} K: measured {k0:.4} mm2/s; 1/3 rule {kt:.4} ({:+.1} %), far field \
             {kf:.4} ({:+.1} %), wet bulb {tw:.1} K",
            100.0 * et,
            100.0 * ef
        );
        assert!(et.abs() < ef.abs(), "the 1/3 rule must be closer at {tg} K");
        assert!((0.0..0.75).contains(&et), "1/3 rule at {tg} K: {et:+.3}");
        assert!((0.6..1.05).contains(&ef), "far field at {tg} K: {ef:+.3}");
        assert!((320.0..350.0).contains(&tw), "wet bulb at {tg} K: {tw}");
    }
}
