//! Verification of the thermal strain loads (MS5 part 3): a temperature
//! change and a difference across the thickness on beams and plates,
//! against the stress-free thermal deformation of a statically determinate
//! structure and the force in a fully held bar. And the transient case
//! keys: histories, damping, the direct and the modal method.
//!
//! References: Timoshenko & Goodier, *Theory of Elasticity* (1970), the
//! free thermal expansion and the thermal curvature α ΔT / h of beams and
//! plates (a linear temperature across the thickness of a free plate bends
//! it into a sphere without stress); Boley & Weiner, *Theory of Thermal
//! Stresses* (1960).

#![allow(clippy::needless_range_loop)]

use foehncraft_structure::case::{self, Case};
use foehncraft_structure::codes;

const E: f64 = 206e9;
const ALPHA: f64 = 1.2e-5;
const L: f64 = 0.03;
const H: f64 = 0.2e-3;

fn rel(a: f64, b: f64) -> f64 {
    ((a - b) / b).abs()
}

fn beam_case(theory: &str, supports: &str, temps: &str) -> Case {
    let text = format!(
        r#"
[analysis]
static = true
[materials.steel]
young = {E}
poisson = 0.3
density = 7850.0
expansion = {ALPHA}
[sections.reed]
shape = "rectangle"
width = 6e-3
thickness = {H}
[[beams]]
name = "reed"
from = [0.0, 0.0, 0.0]
to = [{L}, 0.0, 0.0]
elements = 10
section = "reed"
material = "steel"
theory = "{theory}"
[points]
tip = [{L}, 0.0, 0.0]
{supports}
{temps}
"#
    );
    Case::from_toml(&text).unwrap()
}

const CLAMP_ROOT: &str = r#"
[[supports]]
at = [0.0, 0.0, 0.0]
fix = ["all"]
"#;

fn tip(out: &case::Outcome) -> [f64; 6] {
    let s = out.statics.as_ref().unwrap();
    let n = out.built.points["tip"];
    let mut v = [0.0; 6];
    v.copy_from_slice(&s.u[6 * n..6 * n + 6]);
    v
}

#[test]
fn free_beam_expands_and_bends_without_stress() {
    for theory in ["euler-bernoulli", "timoshenko"] {
        // ΔT = 80 K: the tip moves α ΔT L, nothing else
        let out = case::run(&beam_case(
            theory,
            CLAMP_ROOT,
            "[[temperatures]]\npart = \"reed\"\nchange = 80.0\n",
        ))
        .unwrap();
        let u = tip(&out);
        assert!(rel(u[0], ALPHA * 80.0 * L) < 1e-10, "{theory}: {u:?}");
        assert!(u[2].abs() < 1e-12 * L, "{theory}: {u:?}");
        // a difference of 5 K across the thickness, +z hotter: the
        // curvature α ΔT / h bends the tip towards −z by κ L²/2
        let out = case::run(&beam_case(
            theory,
            CLAMP_ROOT,
            "[[temperatures]]\npart = \"reed\"\ndifference = 5.0\n",
        ))
        .unwrap();
        let u = tip(&out);
        let kappa = ALPHA * 5.0 / H;
        assert!(rel(u[2], -kappa * L * L / 2.0) < 1e-9, "{theory}: {u:?}");
        // the slope θy = −dw/dx = κ L
        assert!(rel(u[4], kappa * L) < 1e-9, "{theory}: {u:?}");
        // a free (determinate) structure carries no load: no reaction
        let rf = out.to_json("t")["static"]["reaction_force"].clone();
        for d in 0..3 {
            assert!(rf[d].as_f64().unwrap().abs() < 1e-9, "{rf}");
        }
    }
}

#[test]
fn held_bar_carries_the_thermal_force() {
    // both ends held along x: the bar is pressed by E A α ΔT
    let both = r#"
[[supports]]
at = [0.0, 0.0, 0.0]
fix = ["all"]
[[supports]]
at = [0.03, 0.0, 0.0]
fix = ["ux"]
"#;
    let out = case::run(&beam_case(
        "timoshenko",
        both,
        "[[temperatures]]\npart = \"reed\"\nchange = 50.0\n",
    ))
    .unwrap();
    let s = out.statics.as_ref().unwrap();
    let n = out.built.points["tip"];
    let force = E * 6e-3 * H * ALPHA * 50.0;
    // the support at the tip pushes back along −x
    assert!(
        rel(s.reactions[6 * n], -force) < 1e-9,
        "{} {force}",
        s.reactions[6 * n]
    );
    assert!(tip(&out)[0].abs() < 1e-15);
}

fn plate_case(element: &str, temps: &str) -> Case {
    let text = format!(
        r#"
[analysis]
static = true
[materials.steel]
young = {E}
poisson = 0.3
density = 7850.0
expansion = {ALPHA}
[[plates]]
name = "p"
origin = [0.0, 0.0, 0.0]
edge_a = [0.02, 0.0, 0.0]
edge_b = [0.0, 0.01, 0.0]
divisions = [8, 4]
thickness = {H}
material = "steel"
element = "{element}"
[[supports]]
at = [0.0, 0.0, 0.0]
fix = ["all"]
[[supports]]   # the in-plane rotation (the drilling springs do not hold it)
at = [0.02, 0.0, 0.0]
fix = ["uy"]
[points]
far = [0.02, 0.01, 0.0]
mid = [0.01, 0.005, 0.0]
{temps}
"#
    );
    Case::from_toml(&text).unwrap()
}

#[test]
fn free_plate_expands_and_bends_into_a_sphere() {
    for element in ["mitc4", "mitc3+"] {
        // held at one corner and against in-plane rotation: free to deform
        let out = case::run(&plate_case(
            element,
            "[[temperatures]]\npart = \"p\"\nchange = 60.0\n",
        ))
        .unwrap();
        let s = out.statics.as_ref().unwrap();
        for (name, x) in [("far", [0.02, 0.01]), ("mid", [0.01, 0.005])] {
            let n = out.built.points[name];
            for d in 0..2 {
                let want = ALPHA * 60.0 * x[d];
                assert!(
                    rel(s.u[6 * n + d], want) < 1e-8,
                    "{element} {name}: {} {want}",
                    s.u[6 * n + d]
                );
            }
            assert!(s.u[6 * n + 2].abs() < 1e-12, "{element}");
        }
        // top hotter by 4 K: w = −κ (x² + y²) / 2 with κ = α ΔT / h
        let out = case::run(&plate_case(
            element,
            "[[temperatures]]\npart = \"p\"\ndifference = 4.0\n",
        ))
        .unwrap();
        let s = out.statics.as_ref().unwrap();
        let kappa = ALPHA * 4.0 / H;
        for (name, x) in [("far", [0.02, 0.01]), ("mid", [0.01, 0.005])] {
            let n = out.built.points[name];
            let want = -kappa * (x[0] * x[0] + x[1] * x[1]) / 2.0;
            let w = s.u[6 * n + 2];
            eprintln!("{element} {name}: w = {w:.6e}, sphere {want:.6e}");
            assert!(rel(w, want) < 1e-6, "{element} {name}: {w} {want}");
        }
    }
}

#[test]
fn thermal_case_errors() {
    // no expansion coefficient
    let text = r#"
[analysis]
static = true
[materials.steel]
young = 206e9
poisson = 0.3
density = 7850.0
[sections.round]
shape = "circle"
diameter = 1e-3
[[beams]]
name = "b"
from = [0.0, 0.0, 0.0]
to = [0.01, 0.0, 0.0]
elements = 2
section = "round"
material = "steel"
[[supports]]
at = [0.0, 0.0, 0.0]
fix = ["all"]
[[temperatures]]
part = "b"
change = 10.0
"#;
    let e = case::run(&Case::from_toml(text).unwrap()).err().unwrap();
    assert_eq!(e.code, codes::CASE);
    assert!(e.message.contains("expansion"), "{e}");
    // a part that does not exist, a history that does not exist
    for extra in [
        "[[temperatures]]\npart = \"nope\"\nchange = 1.0\n",
        "[transient]\nend_time = 0.01\ntime_step = 1e-3\n[[loads]]\nat = [0.03, 0.0, 0.0]\nforce = [0.0, 0.0, 1.0]\nhistory = \"x\"\n",
        "[transient]\nend_time = 0.01\ntime_step = 1e-3\n[histories.x]\ntable = [[0.0, 0.0], [0.0, 1.0]]\n[[loads]]\nat = [0.03, 0.0, 0.0]\nforce = [0.0, 0.0, 1.0]\nhistory = \"x\"\n",
        "[transient]\nend_time = 0.01\ntime_step = 1e-3\nmethod = \"modal\"\n",
        "[transient]\nend_time = 0.01\ntime_step = 1e-3\n[damping]\nrayleigh = [1.0, 0.0]\nratios = [[1.0, 0.01], [2.0, 0.01]]\n",
    ] {
        let Err(e) = case::run(&beam_case("timoshenko", CLAMP_ROOT, extra)) else {
            panic!("accepted: {extra}");
        };
        assert_eq!(e.code, codes::CASE, "{extra}: {e}");
    }
}

#[test]
fn transient_case_direct_and_modal_agree_on_a_ramped_load() {
    // a tip force ramped up over 5 ms and held, with 2 % damping at 100 Hz
    // and 2 kHz: both methods give the same tip history
    let keys = |method: &str, dt: f64| {
        format!(
            r#"
[transient]
end_time = 0.02
time_step = {dt}
method = "{method}"
modes = 8
rho_infinity = 1.0
[damping]
ratios = [[100.0, 0.02], [2000.0, 0.02]]
[histories.ramp]
table = [[0.0, 0.0], [0.005, 1.0]]
[[loads]]
at = [0.03, 0.0, 0.0]
force = [0.0, 0.0, 0.01]
history = "ramp"
"#
        )
    };
    let direct = case::run(&beam_case(
        "euler-bernoulli",
        CLAMP_ROOT,
        &keys("direct", 2e-6),
    ))
    .unwrap();
    let modal = case::run(&beam_case(
        "euler-bernoulli",
        CLAMP_ROOT,
        &keys("modal", 1e-4),
    ))
    .unwrap();
    let jd = direct.to_json("d");
    let jm = modal.to_json("m");
    let last = |j: &serde_json::Value| {
        let h = j["transient"]["points"]["tip"]["displacement"]
            .as_array()
            .unwrap()
            .clone();
        h[h.len() - 1][2].as_f64().unwrap()
    };
    let (wd, wm) = (last(&jd), last(&jm));
    eprintln!("tip at 20 ms: direct {wd:.6e}, modal {wm:.6e}");
    assert!(rel(wd, wm) < 2e-3, "{wd} {wm}");
    // the static deflection under the full force, which the damped
    // response approaches
    let ei = E * 6e-3 * H.powi(3) / 12.0;
    let st = 0.01 * L.powi(3) / (3.0 * ei);
    assert!(rel(wd, st) < 0.2, "{wd} {st}");
    // outputs: at most 1001 kept times, energies, a CSV with one row each
    let times = jd["transient"]["times"].as_array().unwrap().len();
    assert!(times <= 1001 && times > 100, "{times}");
    let csv = direct.to_csv().unwrap();
    assert_eq!(csv.lines().count(), times + 1);
    assert!(csv.starts_with("time,kinetic_energy,strain_energy,tip_ux"));
    assert_eq!(jm["transient"]["modes"].as_array().unwrap().len(), 8);
}
