//! The Corrected Distortion model of the drops (Nguyen, Dahms, Pickett and
//! Tagliante 2022; ECN Spray A, SA4, D316) in the spray:
//!
//! - the evaporation rate of a drop of distortion y is that of the sphere
//!   times the Sherwood factor of the kernel, exactly (the d² rate is
//!   proportional to Sh, whatever the film thickening of Abramzon and
//!   Sirignano does to it), and the other drag laws leave it alone;
//! - a flattened drop evaporates faster and a prolate one about as fast as
//!   a sphere (fig. 18 of the paper);
//! - the viscosity ratio follows the library liquid's viscosity at the
//!   drop temperature.

use foehncraft_kernels::distortion as cdm;
use foehncraft_kernels::realfluid::fuels::N_DODECANE;
use foehncraft_kernels::species::R;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, DragLaw, Env, EvaporationModel, GasState, Injector, Spray,
    SpraySection, Vapour, ideal_gas_species,
};

fn gas(t: f64, p0: f64, slip: f64) -> GasState {
    GasState {
        u: [slip, 0.0, 0.0],
        rho: p0 * 0.028 / (R * t),
        mu: 3.9e-5,
        t,
        cp: 1150.0,
        k: 0.061,
        tke: 0.0,
        eps: 0.0,
        y_v: 0.0,
        p0,
        m_gas: 0.028_013_4,
        diff: 2.5e-6,
    }
}

/// A spray of the library n-dodecane at 363 K (Raoult surface) in
/// nitrogen with the drag law `drag`.
fn spray(drag: DragLaw) -> Spray {
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
    cfg.breakup.drag = drag;
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
    sp.set_carrier(vec![(
        ideal_gas_species(
            include_str!("../../foehncraft-flow/data/species-gri30.toml"),
            "N2",
        )
        .unwrap(),
        1.0,
    )]);
    sp
}

#[test]
fn the_rate_of_a_distorted_drop_is_the_sherwood_factor() {
    let cd = spray(DragLaw::CorrectedDistortion);
    let dynamic = spray(DragLaw::Dynamic);
    let g = gas(900.0, 2e6, 0.0);
    let (d, tp, re) = (20e-6, 450.0, 60.0);
    let sphere = dynamic.evaporation_of(&g, d, tp, re, 1e-7, 0.0).unwrap();
    // the dynamic drag does not touch the evaporation, whatever y
    let flat = dynamic.evaporation_of(&g, d, tp, re, 1e-7, 1.0).unwrap();
    assert_eq!(flat.k.to_bits(), sphere.k.to_bits());
    let mut last = 0.0;
    for y in [-1.0, -0.5, 0.0, 0.25, 0.5, 0.75, 1.0] {
        let e = cd.evaporation_of(&g, d, tp, re, 1e-7, y).unwrap();
        // the Sherwood factor the spray gave the kernel, read back from
        // the Sherwood numbers: the rates differ by exactly that factor
        let f = e.sh / sphere.sh;
        assert!((e.k / sphere.k / f - 1.0).abs() < 1e-12, "{y}");
        if y < 0.0 {
            // a prolate drop: about the sphere (fig. 18: 1 to 1.2)
            assert!(f > 0.95 && f < 1.25, "{y}: {f}");
        } else {
            // a flattened drop: faster the flatter
            assert!(f > last, "{y}: {f}");
            last = f;
        }
    }
    // at y = 0 only the circulation inside the viscous drop is left (λ
    // about 9 here): 11 %; the most flattened drop 1.5 to 4 times
    let f0 = cd.evaporation_of(&g, d, tp, re, 1e-7, 0.0).unwrap().sh / sphere.sh;
    let f1 = cd.evaporation_of(&g, d, tp, re, 1e-7, 1.0).unwrap().sh / sphere.sh;
    assert!(f0 > 1.0 && f0 < 1.2, "{f0}");
    assert!(f1 > 1.5 && f1 < 4.0, "{f1}");
}

#[test]
fn the_factor_is_the_kernel_factor_at_the_film_state() {
    // the film rule rescales Re to the film viscosity; the factor the spray
    // applies is the kernel's at that Re, Sc and λ = μ_l(T_p)/μ_film
    let cd = spray(DragLaw::CorrectedDistortion);
    let dynamic = spray(DragLaw::Dynamic);
    let g = gas(900.0, 2e6, 0.0);
    let (d, tp, re, y) = (20e-6, 450.0, 60.0, 0.7);
    let a = cd.evaporation_of(&g, d, tp, re, 1e-7, y).unwrap();
    let b = dynamic.evaporation_of(&g, d, tp, re, 1e-7, y).unwrap();
    let ys = foehncraft_kernels::evaporation::surface_fraction(
        cd.saturation_pressure(tp),
        g.p0,
        N_DODECANE.ideal_gas.molar_mass,
        g.m_gas,
    );
    let [rho, mu, _, _, diff, _] = cd.film_properties(&g, tp, ys);
    let mu_l = foehncraft_kernels::liquids::by_name("n-dodecane")
        .unwrap()
        .viscosity(tp);
    let f = cdm::transfer_factor(re * g.mu / mu, mu / (rho * diff), mu_l / mu, y);
    assert!((a.sh / b.sh / f - 1.0).abs() < 1e-12, "{} {f}", a.sh / b.sh);
}
