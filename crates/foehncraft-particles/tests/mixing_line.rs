//! The adiabatic mixing line of n-dodecane and nitrogen at Spray A
//! conditions and the mixing-limited liquid length from it (D319, VA.4,
//! VA.5): the thermodynamics checked without CFD.

use foehncraft_kernels::realfluid::{
    Flash, MAX, Mixing, Mixture, Streams, fuels::N_DODECANE, gases,
};
use foehncraft_kernels::species::{R, Species};
use foehncraft_particles::{ideal_gas_species, siebers_liquid_length};

fn nitrogen() -> Species {
    ideal_gas_species(
        include_str!("../../foehncraft-flow/data/species-gri30.toml"),
        "N2",
    )
    .unwrap()
}

fn mixture() -> Mixture {
    let mut m = Mixture::new(&[N_DODECANE.component, gases::N2]);
    m.set_k(0, 1, N_DODECANE.k_with("N2"));
    m
}

/// Runs `f` with the Spray A streams: 363 K n-dodecane into nitrogen at
/// `t_gas` and `p`.
fn with_line<T>(t_gas: f64, p: f64, f: impl FnOnce(&Streams) -> T) -> T {
    let mix = mixture();
    let n2 = nitrogen();
    let fuel = N_DODECANE.ideal_gas;
    let ideal = |i: usize, t: f64| {
        R * t
            * if i == 0 {
                fuel.h_over_rt(t)
            } else {
                n2.h_over_rt(t)
            }
    };
    let mut gas = [0.0; MAX];
    gas[1] = 1.0;
    let s = Streams {
        mix: &mix,
        ideal: &ideal,
        p,
        t_fuel: 363.0,
        t_gas,
        gas,
    };
    f(&s)
}

#[test]
fn the_ideal_gas_of_dodecane_follows_its_heat_capacity_table() {
    // Scott (1974) by the CH₂ increment (tools/data/dodecane_ideal_gas.py):
    // the NASA-7 fit within 0.3 %, and h(900) − h(363) the integral of c_p
    let s = N_DODECANE.ideal_gas;
    for (t, cp) in [
        (300.0, 279.70),
        (500.0, 426.21),
        (700.0, 532.61),
        (900.0, 605.44),
        (1200.0, 675.70),
    ] {
        let got = R * s.cp_over_r(t);
        assert!((got / cp - 1.0).abs() < 3e-3, "{t} K: {got} {cp}");
    }
    let h = |t: f64| R * t * s.h_over_rt(t);
    let n = 2000;
    let mut int = 0.0;
    for k in 0..n {
        let t = 363.0 + (f64::from(k) + 0.5) * 537.0 / f64::from(n);
        int += R * s.cp_over_r(t) * 537.0 / f64::from(n);
    }
    assert!(((h(900.0) - h(363.0)) / int - 1.0).abs() < 1e-6);
    // the enthalpy of formation at 298.15 K (NIST WebBook)
    assert!((h(298.15) + 290.9e3).abs() < 1.0, "{}", h(298.15));
}

#[test]
fn va4_the_mixing_line_closes_the_energy_balance() {
    // 363 K n-dodecane into 900 K nitrogen at 6 MPa: at every mixture
    // fraction the state found has the enthalpy of the mixed streams to
    // round-off, frozen, in equilibrium and ideal, and the end points are
    // the streams
    with_line(900.0, 6e6, |line| {
        for model in [Mixing::Frozen, Mixing::Equilibrium, Mixing::Ideal] {
            let (hf, hg) = line.stream_enthalpies(model);
            let scale = (hg - hf).abs();
            let mut last = f64::INFINITY;
            for k in 0..=40 {
                let z = f64::from(k) / 40.0;
                let pt = line.point(z, model);
                assert!(
                    pt.residual.abs() < 1e-13 * scale,
                    "{model:?} Z {z}: {:e}",
                    pt.residual
                );
                if model == Mixing::Equilibrium {
                    // the equilibrium line cools monotonically
                    assert!(pt.t < last, "Z {z}: {} after {last}", pt.t);
                }
                last = pt.t;
            }
            assert_eq!(line.point(0.0, model).t, 900.0, "{model:?}");
            assert!((line.point(1.0, model).t - 363.0).abs() < 1e-9, "{model:?}");
        }
    });
}

#[test]
fn va4_the_two_phase_region_of_the_mixing_line() {
    with_line(900.0, 6e6, |line| {
        let eq = line.saturation(Mixing::Equilibrium).unwrap();
        let frozen = line.saturation(Mixing::Frozen).unwrap();
        let ideal = line.saturation(Mixing::Ideal).unwrap();
        // the frozen line is the equilibrium line up to Z_sat, and its one
        // phase fails the stability test exactly where the flash splits
        assert!((eq - frozen).abs() < 1e-6, "{eq} {frozen}");
        for z in [0.1, 0.2, 0.3] {
            let a = line.point(z, Mixing::Frozen);
            let b = line.point(z, Mixing::Equilibrium);
            assert_eq!(a.t, b.t, "Z {z}");
        }
        // at 6 MPa the vapour holds more fuel than Raoult's law says, so the
        // last liquid evaporates later on the line (own result: 0.357 and
        // 0.320)
        assert!((0.34..0.37).contains(&eq), "{eq}");
        assert!((0.30..0.34).contains(&ideal), "{ideal}");
        assert!(eq > ideal + 0.02);
        // inside the region: two phases of equal fugacity, the liquid
        // carrying nitrogen, and more liquid the richer the mixture
        let mut last = 0.0;
        for z in [0.4, 0.6, 0.8] {
            let pt = line.point(z, Mixing::Equilibrium);
            assert!(pt.split && pt.liquid > last, "Z {z}: {pt:?}");
            last = pt.liquid;
            let x = line.composition(z);
            let Flash::Two(f) = line.mix.flash(&x, pt.t, line.p, None) else {
                panic!("Z {z}: no split")
            };
            for i in 0..2 {
                let fl = f.l.ln_phi[i] + f.x[i].ln();
                let fv = f.v.ln_phi[i] + f.y[i].ln();
                assert!((fl - fv).abs() < 1e-9, "Z {z}, {i}: {fl} {fv}");
            }
            assert!(f.x[1] > 0.05, "Z {z}: dissolved N2 {}", f.x[1]);
        }
    });
}

/// Measured liquid lengths of Spray A (diffused back-illumination,
/// injector 210675, 150 MPa, 0 % O₂; ECN, ecn.sandia.gov/data/dbi675/):
/// ambient temperature (K), density (kg/m³), liquid length (m).
const MEASURED: [(f64, f64, f64); 4] = [
    (900.0, 22.8, 11.7e-3),
    (700.0, 22.8, 16.3e-3),
    (1200.0, 22.8, 7.4e-3),
    (900.0, 7.6, 20.6e-3),
];

#[test]
fn va5_the_siebers_liquid_length_with_our_thermodynamics() {
    // Siebers (1999) with B from our mixing line, the Spray A nozzle
    // (d = 90 µm, C_a = 0.95, full angle 21.5°: Pickett et al. 2011) and
    // the fuel density of Lemmon and Huber (2004) at 363 K, 6 MPa
    let len = |t_gas: f64, rho: f64, model: Mixing| {
        let mix = mixture();
        let mut x = [0.0; MAX];
        x[1] = 1.0;
        let p = mix.pressure(&x, t_gas, gases::N2.molar_mass / rho);
        let z = with_line(t_gas, p, |l| l.saturation(model).unwrap());
        siebers_liquid_length(z, 703.71, rho, 0.95, 0.090e-3, 21.5)
    };
    // own results (millimetres; recorded here so a change shows):
    // equilibrium 7.16, 11.33, 4.86, 10.67; ideal 8.21, 12.95, 5.67, 11.27
    let want_eq = [7.163, 11.333, 4.858, 10.668];
    let want_id = [8.214, 12.952, 5.669, 11.267];
    let mut eq = [0.0; 4];
    let mut id = [0.0; 4];
    for (k, (t, rho, measured)) in MEASURED.into_iter().enumerate() {
        eq[k] = len(t, rho, Mixing::Equilibrium);
        id[k] = len(t, rho, Mixing::Ideal);
        assert!(
            (eq[k] * 1e3 / want_eq[k] - 1.0).abs() < 2e-3,
            "{t} K {rho}: {}",
            eq[k]
        );
        assert!(
            (id[k] * 1e3 / want_id[k] - 1.0).abs() < 2e-3,
            "{t} K {rho}: {}",
            id[k]
        );
        // the real fluid shortens the liquid length by 5 to 14 %; with
        // C_L = 0.62 the equilibrium length is 30 to 48 % below the
        // measured one and the ideal 20 to 45 % (the gap of VA.5)
        let r = eq[k] / id[k];
        assert!((0.85..0.96).contains(&r), "{t} K {rho}: {r}");
        assert!(
            eq[k] < 0.72 * measured && id[k] < 0.82 * measured,
            "{t} K {rho}"
        );
    }
    // the trends: shorter in hotter and in denser gas, as measured
    for l in [eq, id] {
        assert!(l[1] > l[0] && l[0] > l[2] && l[3] > l[0]);
    }
}
