//! Checks of the chemistry of SF2: the analytic Jacobian against finite
//! differences, the detailed balance of the equilibrium constants, the
//! Gibbs equilibrium against the adiabatic flame temperature, the reactor
//! against the equilibrium, element conservation.

use super::equilibrium;
use super::mechanism::Mechanism;
use super::reactor::{self, ReactorCase, ReactorSection};
use foehncraft_kernels::chemistry as kc;
use std::collections::BTreeMap;

fn map(v: &[(&str, f64)]) -> BTreeMap<String, f64> {
    v.iter().map(|(k, x)| (k.to_string(), *x)).collect()
}

fn air_methane(mech: &str, phi: f64, t: f64, end: f64) -> ReactorCase {
    ReactorCase {
        reactor: ReactorSection {
            mechanism: mech.into(),
            pressure: 101_325.0,
            temperature: t,
            end_time: end,
            mole_fractions: None,
            mass_fractions: None,
            equivalence_ratio: Some(phi),
            fuel: Some(map(&[("CH4", 1.0)])),
            oxidiser: Some(map(&[("O2", 0.21), ("N2", 0.79)])),
            samples: 200,
            spacing: reactor::Spacing::Linear,
            integrator: Default::default(),
            check: None,
        },
    }
}

/// A partly burnt state of the mechanism (the equilibrium at 2000 K mixed
/// with the fresh gas), so that radicals and every reaction are active.
fn partly_burnt(m: &Mechanism, t: f64) -> Vec<f64> {
    let c = air_methane(&m.name, 1.0, 300.0, 1.0);
    let x0 = reactor::mixture(m, &c.reactor).unwrap();
    let y0 = reactor::mass_of_mole(m, &x0);
    let eq = equilibrium::tp(m, &y0, 2000.0, 101_325.0, None).unwrap();
    let mut z: Vec<f64> = y0
        .iter()
        .zip(&eq.mass_fractions)
        .map(|(a, b)| 0.5 * a + 0.5 * b)
        .collect();
    z.push(t);
    z
}

#[test]
fn the_analytic_jacobian_matches_finite_differences() {
    for (mech, t) in [
        ("gri30", 1600.0),
        ("gri30", 900.0),
        ("jones-lindstedt-4", 1700.0),
        ("westbrook-dryer-2", 1700.0),
    ] {
        let m = Mechanism::load(mech, None).unwrap();
        let prog = m.program();
        let n = m.len();
        let nz = n + 1;
        let p = 101_325.0;
        let z = partly_burnt(&m, t);
        let mut f = vec![0.0; nz];
        let mut jac = vec![0.0; nz * nz];
        let mut w = vec![0.0; kc::rhs_work_len(n)];
        assert!(kc::reactor_rhs(
            &prog,
            p,
            &z,
            &mut f,
            Some(&mut jac),
            &mut w
        ));
        let mut f2 = vec![0.0; nz];
        assert!(kc::reactor_rhs(&prog, p, &z, &mut f2, None, &mut w));
        assert_eq!(f, f2, "{mech}: f with and without the Jacobian");
        let (mut fp, mut fm) = (vec![0.0; nz], vec![0.0; nz]);
        // the scale of each row
        let mut worst: f64 = 0.0;
        for c in 0..nz {
            let h = if c == n {
                1e-4 * t
            } else {
                1e-6 * z[c].abs() + 1e-9
            };
            let mut zp = z.clone();
            let mut zm = z.clone();
            zp[c] += h;
            zm[c] -= h;
            kc::reactor_rhs(&prog, p, &zp, &mut fp, None, &mut w);
            kc::reactor_rhs(&prog, p, &zm, &mut fm, None, &mut w);
            for r in 0..nz {
                let fd = (fp[r] - fm[r]) / (2.0 * h);
                let an = jac[r * nz + c];
                // the row scale: the largest entry of the row
                let scale = (0..nz)
                    .map(|k| jac[r * nz + k].abs() * if k == n { t } else { 1.0 })
                    .fold(0.0, f64::max)
                    .max(1e-30);
                let colw = if c == n { t } else { 1.0 };
                let e = (fd - an).abs() * colw / scale;
                worst = worst.max(e);
                assert!(
                    e < 1e-4,
                    "{mech} T = {t}: J[{r}][{c}] analytic {an:e}, differences {fd:e}"
                );
            }
        }
        eprintln!("{mech} at {t} K: largest scaled Jacobian error {worst:.2e}");
    }
}

#[test]
fn reversible_reactions_stop_at_the_gibbs_equilibrium() {
    // at the Gibbs equilibrium the net rate of every reversible reaction
    // of GRI-Mech vanishes: the equilibrium constants of the kinetics and
    // the minimisation agree
    let m = Mechanism::load("gri30", None).unwrap();
    let prog = m.program();
    let n = m.len();
    let c = air_methane("gri30", 1.0, 300.0, 1.0);
    let x0 = reactor::mixture(&m, &c.reactor).unwrap();
    let y0 = reactor::mass_of_mole(&m, &x0);
    for t in [1500.0, 2200.0] {
        let eq = equilibrium::tp(&m, &y0, t, 101_325.0, None).unwrap();
        let mut conc = vec![0.0; n];
        kc::concentrations(&m.data, 101_325.0, t, &eq.mass_fractions, &mut conc);
        // net and forward rates per reaction, one at a time
        for (j, rx) in m.reactions.iter().enumerate() {
            if rx.reversible == 0 {
                continue;
            }
            let one = kc::Program {
                reactions: std::slice::from_ref(rx),
                ..prog
            };
            let mut wdot = vec![0.0; n];
            let mut w = vec![0.0; kc::rates_work_len(n)];
            kc::rates(&one, t, &conc, &mut wdot, None, None, &mut w);
            let fwd = kc::Program {
                reactions: &[kc::Reaction {
                    reversible: 0,
                    ..*rx
                }],
                ..prog
            };
            let mut wf = vec![0.0; n];
            kc::rates(&fwd, t, &conc, &mut wf, None, None, &mut w);
            let k = prog.terms[rx.terms[0] as usize].species as usize;
            let rel = wdot[k].abs() / wf[k].abs().max(1e-300);
            assert!(
                rel < 1e-6,
                "T = {t}: reaction {} `{}` net/forward = {rel:e}",
                j + 1,
                m.equations[j]
            );
        }
    }
}

/// The adiabatic flame temperature of stoichiometric methane–air from 298
/// K at 1 atm with dissociation: 2226 K (Turns 2012, "An introduction to
/// combustion", 3rd ed., chapter 2, with the JANAF data); within 10 K.
#[test]
fn adiabatic_flame_temperature_of_methane_and_air() {
    let m = Mechanism::load("gri30", None).unwrap();
    let c = air_methane("gri30", 1.0, 298.15, 1.0);
    let x0 = reactor::mixture(&m, &c.reactor).unwrap();
    let i = m.index("CH4").unwrap();
    assert!((x0[i] - 1.0 / 10.52).abs() < 1e-3, "{}", x0[i]);
    let y0 = reactor::mass_of_mole(&m, &x0);
    let eq = equilibrium::hp(&m, &y0, 298.15, 101_325.0).unwrap();
    assert!(
        (eq.temperature - 2226.0).abs() < 10.0,
        "T_ad = {}",
        eq.temperature
    );
    // the enthalpy is that of the fresh mixture
    let h0 = equilibrium::enthalpy(&m, &y0, 298.15);
    let h1 = equilibrium::enthalpy(&m, &eq.mass_fractions, eq.temperature);
    assert!((h1 - h0).abs() < 1e-6 * h0.abs().max(1e3));
    // CO and NO of the equilibrium at T_ad, of the order of 1e-3 (CO
    // 0.009, NO 0.002 by mole, Turns 2012, figure 2.14 shape): positive
    // and below 2 %
    for s in ["CO", "NO", "OH"] {
        let x = eq.mole_fractions[m.index(s).unwrap()];
        assert!(x > 1e-4 && x < 0.02, "{s}: {x}");
    }
}

#[test]
fn the_reactor_ends_at_the_equilibrium_and_conserves_the_elements() {
    // GRI-Mech from 1400 K: ignition, then the approach to equilibrium
    let c = air_methane("gri30", 1.0, 1400.0, 0.2);
    let r = reactor::run(&c, None).unwrap();
    assert!(r.ignition_delay.is_some());
    assert!(
        r.temperature_difference.abs() < 0.5,
        "T_end − T_eq = {}",
        r.temperature_difference
    );
    assert!(
        r.largest_mole_fraction_difference < 2e-4,
        "{}",
        r.largest_mole_fraction_difference
    );
    assert!(r.element_drift < 1e-9, "{}", r.element_drift);
    // the four-step mechanism ends at the equilibrium of its species
    let c = air_methane("jones-lindstedt-4", 1.0, 1400.0, 0.2);
    let r = reactor::run(&c, None).unwrap();
    assert!(
        r.temperature_difference.abs() < 1.0,
        "JL4 {}",
        r.temperature_difference
    );
    assert!(r.element_drift < 1e-9);
}

#[test]
fn reactor_inputs_are_checked() {
    let mut c = air_methane("gri30", 1.0, 1400.0, 0.01);
    c.reactor.mole_fractions = Some(map(&[("CH4", 1.0)]));
    assert!(reactor::run(&c, None).is_err());
    let mut c = air_methane("gri30", 1.0, 1400.0, 0.01);
    c.reactor.fuel = Some(map(&[("XYZ", 1.0)]));
    assert!(reactor::run(&c, None).is_err());
    let mut c = air_methane("gri30", 1.0, -3.0, 0.01);
    c.reactor.samples = 3;
    assert!(reactor::run(&c, None).is_err());
}
