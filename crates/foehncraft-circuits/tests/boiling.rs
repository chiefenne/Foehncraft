//! Flow boiling in the circuits (SF7b): the regimes along an evaporator
//! tube, the two boiling correlations side by side, the friction and void
//! models, the CHF warning, and the refusals of the case reader.

use crate::exact::{cases, solve};
use foehncraft_circuits::case::Case;
use foehncraft_circuits::{codes, output};
use foehncraft_kernels::steam;

fn variant(name: &str, edit: impl Fn(&str) -> String) -> Case {
    let text = std::fs::read_to_string(cases().join(name)).unwrap();
    Case::from_toml(&edit(&text)).unwrap()
}

#[test]
fn evaporator_passes_through_the_regimes_in_order() {
    let (_, sol) = solve("evaporator_tube.toml");
    let t = sol.tube("evaporator").unwrap();
    let t_sat = steam::t_sat(7e6);
    let order = ["liquid", "subcooled boiling", "saturated boiling"];
    let mut last = 0;
    for s in &t.segments {
        let k = order
            .iter()
            .position(|r| *r == s.regime)
            .unwrap_or_else(|| panic!("{}", s.regime));
        assert!(k >= last, "regime {} after {}", s.regime, order[last]);
        last = k;
        // the flowing quality lies between zero and the equilibrium quality
        assert!(s.x >= 0.0 && s.x >= s.x_eq - 1e-12);
        // nucleate boiling: the wall a few kelvin above saturation
        if s.regime != "liquid" {
            assert!(
                s.t_wall > t_sat && s.t_wall < t_sat + 30.0,
                "T_w {}",
                s.t_wall
            );
        }
    }
    assert_eq!(last, 2);
    // vapour appears before the bulk reaches saturation (subcooled boiling)
    let first = t.segments.iter().find(|s| s.x > 0.0).unwrap();
    assert!(first.x_eq < -0.01, "first vapour at x_eq {}", first.x_eq);
    // energy: the outlet quality of the heat balance
    let h_in = steam::state_pt(7e6, 530.0).unwrap().h;
    let q = 4e5 * std::f64::consts::PI * 0.02 * 12.0;
    assert!((t.h_end - (h_in + q / 0.35)).abs() < 1e-6);
    assert!(sol.warnings.is_empty(), "{:?}", sol.warnings);
}

#[test]
fn chen_and_gungor_winterton_agree_on_the_wall_superheat() {
    let (_, chen) = solve("evaporator_tube.toml");
    let gw = output::run(&variant("evaporator_tube.toml", |t| {
        t.replace("boiling = \"chen\"", "boiling = \"gungor-winterton\"")
    }))
    .unwrap();
    let t_sat = steam::t_sat(7e6);
    let (a, b) = (
        chen.tube("evaporator").unwrap(),
        gw.tube("evaporator").unwrap(),
    );
    for (sa, sb) in a.segments.iter().zip(&b.segments) {
        if sa.regime == "saturated boiling" && sb.regime == "saturated boiling" {
            let (da, db) = (sa.t_wall - t_sat, sb.t_wall - t_sat);
            // the two correlations, each fitted to some thousand points with
            // about 20 % scatter, give superheats within a factor of two
            assert!(
                da / db > 0.5 && da / db < 2.0,
                "s {} Chen {da} K, GW {db} K",
                sa.s
            );
        }
    }
}

#[test]
fn friedel_and_drift_flux_differ_from_homogeneous_as_expected() {
    let (_, base) = solve("evaporator_tube.toml");
    let homo = output::run(&variant("evaporator_tube.toml", |t| {
        t.replace("void = \"drift-flux\"", "void = \"homogeneous\"")
    }))
    .unwrap();
    let (a, b) = (
        base.tube("evaporator").unwrap(),
        homo.tube("evaporator").unwrap(),
    );
    // slip in upflow: less void, so a heavier column, with drift flux
    assert!(a.dp_grav > b.dp_grav);
    for (sa, sb) in a.segments.iter().zip(&b.segments) {
        assert!(sa.alpha <= sb.alpha + 1e-12);
    }
    let homo_f = output::run(&variant("evaporator_tube.toml", |t| {
        t.replace("friction = \"friedel\"", "friction = \"homogeneous\"")
    }))
    .unwrap();
    // the homogeneous model is known to under-predict the friction of
    // boiling water at moderate pressure against Friedel
    let (fa, fb) = (a.dp_fric, homo_f.tube("evaporator").unwrap().dp_fric);
    assert!(
        fa > fb && fa < 3.0 * fb,
        "Friedel {fa} Pa, homogeneous {fb} Pa"
    );
}

#[test]
fn a_high_heat_flux_reaches_the_chf_and_warns() {
    // 2.5 MW/m2 at a low flow: Bowring's CHF is exceeded near the outlet
    let case = variant("evaporator_tube.toml", |t| {
        t.replace("heat_flux = 4.0e5", "heat_flux = 1.0e6")
            .replace("mass_flow = 0.35", "mass_flow = 0.25")
    });
    let sol = output::run(&case).unwrap();
    let t = sol.tube("evaporator").unwrap();
    let (_, r) = output::min_chf(t).unwrap();
    assert!(r < 1.0, "CHF ratio {r}");
    assert!(sol.warnings.iter().any(|w| w.code == codes::CHF));
    let post: Vec<_> = t
        .segments
        .iter()
        .filter(|s| s.regime == "post-CHF")
        .collect();
    assert!(!post.is_empty());
    // beyond the CHF the wall runs far hotter than in nucleate boiling
    let t_sat = steam::t_sat(7e6);
    assert!(post.iter().all(|s| s.t_wall > t_sat + 50.0));
}

#[test]
fn broken_cases_are_refused() {
    let base = std::fs::read_to_string(cases().join("drum_natural_circulation.toml")).unwrap();
    for (what, edit) in [
        (
            "a tube to a missing node",
            base.replace("to = \"lower_header\"", "to = \"nowhere\""),
        ),
        (
            "a rise above the length",
            base.replace("rise = -22.0", "rise = -30.0"),
        ),
        (
            "no pressure level",
            base.replace("kind = \"drum\"", "kind = \"junction\""),
        ),
        (
            "a drum without feedwater",
            base.replace("feed_temperature = 490.0", ""),
        ),
        ("an unknown key", base.replace("loss = 1.5", "los = 1.5")),
        (
            "a supercritical drum",
            base.replace("pressure = 7.0e6", "pressure = 23.0e6"),
        ),
    ] {
        let e = Case::from_toml(&edit).expect_err(what);
        assert_eq!(e.code, codes::CASE, "{what}: {e}");
    }
}
