//! The deforming mesh (D91, stage MS0 of the moving solids plan): checks
//! MS0.1 (a uniform flow on a deforming mesh stays uniform to round-off,
//! the space conservation law) and MS0.2 (the L22 piston: a slow adiabatic
//! compression keeps p₀ V^γ, the error falls as Δt²), and the refusals
//! MOV001 to MOV004.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

const R_GAS: f64 = 8.314462618 / 0.0289647;
const CP: f64 = 1005.0;
const P_OP: f64 = 101_325.0;
const V0: f64 = 0.02 * 0.04 * 0.002;

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = std::env::temp_dir().join(format!("foehncraft_motion_{name}_{}", std::process::id()));
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

/// The code and message a case text is refused with.
fn refused(text: &str, name: &str) -> (String, String) {
    let dir = std::env::temp_dir().join(format!("foehncraft_motion_{name}_{}", std::process::id()));
    let e = match FlowCase::from_toml(text, &cases()) {
        Err(e) => e,
        Ok(c) => match run(&c, name, &dir, |_| {}) {
            Err(e) => e,
            Ok(_) => panic!("the case was accepted"),
        },
    };
    let _ = std::fs::remove_dir_all(&dir);
    (diag::of_error(&e).code.to_string(), format!("{e:#}"))
}

/// MS0.1: the uniform flow of motion_free_stream.toml on a mesh whose
/// inner nodes swing by up to 3 cm stays uniform to round-off, for the
/// three cell types of the test mesher and the three time schemes.
#[test]
fn uniform_flow_on_a_deforming_mesh_stays_uniform() {
    let base = text("motion_free_stream.toml").replace("end_time = 0.5", "end_time = 0.2");
    for (name, set) in [
        ("hex_bdf2", ("", "")),
        ("poly", ("type = \"hex\"", "type = \"poly\"")),
        ("tet", ("type = \"hex\"", "type = \"tet\"")),
        ("bdf1", ("scheme = \"bdf2\"", "scheme = \"bdf1\"")),
        ("cn", ("scheme = \"bdf2\"", "scheme = \"crank-nicolson\"")),
    ] {
        let t = if set.0.is_empty() {
            base.clone()
        } else {
            base.replace(set.0, set.1)
        };
        let rep = run_case(&parse(&t), name);
        let s = &rep.solves[0];
        let m = s.transient.as_ref().unwrap().mesh_motion.as_ref().unwrap();
        assert!(m.closure_max < 1e-13, "{name}: closure {}", m.closure_max);
        assert!(
            m.volume_ratio_min < 0.9 && m.volume_ratio_max > 1.1,
            "{name}: the mesh did not deform ({m:?})"
        );
        let b = |n: &str| s.boundaries.iter().find(|b| b.name == n).unwrap();
        let (inl, outl) = (b("xmin"), b("xmax"));
        assert!(
            (inl.mass_flow + outl.mass_flow).abs() < 1e-12 * outl.mass_flow,
            "{name}: mass flows {} and {}",
            inl.mass_flow,
            outl.mass_flow
        );
        // p = 0 everywhere: ½ρU² = 0.5 Pa is the scale
        assert!(
            inl.mean_pressure.abs() < 1e-10,
            "{name}: inlet pressure {}",
            inl.mean_pressure
        );
        let r = &s.residuals;
        assert!(
            r.momentum.iter().all(|x| *x < 1e-12) && r.continuity < 1e-12,
            "{name}: residuals {r:?}"
        );
    }
}

/// p₀ at the end of the adiabatic piston run with the step `dt` until
/// `end`, and the exact p₀(0) (V₀/V)^γ of the volume at that time.
fn piston(dt: f64, end: f64) -> (f64, f64, FlowReport) {
    let t = text("motion_piston_adiabatic.toml")
        .replace("time_step = 0.01", &format!("time_step = {dt}"))
        .replace("end_time = 1.0", &format!("end_time = {end}"));
    let rep = run_case(&parse(&t), &format!("piston_{dt}_{end}"));
    let s = &rep.solves[0];
    let v = s
        .transient
        .as_ref()
        .unwrap()
        .mesh_motion
        .as_ref()
        .unwrap()
        .volume;
    let gamma = CP / (CP - R_GAS);
    let exact = P_OP * (V0 / v).powf(gamma);
    (s.density.as_ref().unwrap().p0, exact, rep)
}

/// MS0.2, the L22 piston: compressed to a quarter of its volume, the gas
/// is on the isentrope to the time error of BDF2, which falls as Δt²; the
/// mass stays to round-off and the heat ledger of every step closes.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn piston_compression_keeps_p_v_gamma() {
    let (p1, e1, rep) = piston(0.01, 0.5);
    let (p2, e2, _) = piston(0.005, 0.5);
    let (err1, err2) = ((p1 / e1 - 1.0).abs(), (p2 / e2 - 1.0).abs());
    // the top dead centre: V = V₀/4, p₀ = 4^γ p₀(0) = 705.5 kPa
    assert!((e1 / 705_500.75 - 1.0).abs() < 1e-6, "exact p₀ {e1}");
    assert!(err1 < 1e-3, "p₀ {p1} against {e1}");
    let order = (err1 / err2).log2();
    assert!(
        (1.8..2.2).contains(&order),
        "order {order} ({err1:e}, {err2:e})"
    );
    let s = &rep.solves[0];
    let d = s.density.as_ref().unwrap();
    assert!((d.mass - d.mass_start).abs() < 1e-13 * d.mass_start);
    let t = s.scalars.iter().find(|x| x.name == "temperature").unwrap();
    // uniform, on T = T₀ (V₀/V)^(γ−1) = 522.2 K to the time error
    let te = 300.0 * (e1 / P_OP).powf(R_GAS / CP);
    assert!(
        (t.max - t.min) < 1e-8 * t.max,
        "T from {} to {}",
        t.min,
        t.max
    );
    assert!((t.max / te - 1.0).abs() < 1e-3, "T {} against {te}", t.max);
    let tr = s.transient.as_ref().unwrap();
    assert!(
        tr.heat_relative_max.unwrap() < 1e-8,
        "heat ledger {:?}",
        tr.heat_relative_max
    );
    // the whole period: back at the start pressure
    let (p, e, _) = piston(0.01, 1.0);
    assert!((e / P_OP - 1.0).abs() < 1e-12);
    assert!((p / e - 1.0).abs() < 1e-4, "p₀ {p} after the period");
}

/// MOV001 to MOV004: a steady case, a folding displacement, an inlet that
/// moves out of its plane and a closed volume of constant density that is
/// compressed are refused with their codes.
#[test]
fn deforming_mesh_settings_are_checked() {
    let free = text("motion_free_stream.toml");
    let steady = free.replace(
        "[time]\nscheme = \"bdf2\"\ntime_step = 0.01\nend_time = 0.5\n",
        "",
    );
    let (code, msg) = refused(&steady, "steady");
    assert_eq!(code, "MOV001", "{msg}");
    let fold = free
        .replace("amp = 0.03", "amp = 0.5")
        .replace("end_time = 0.5", "end_time = 0.2");
    let (code, msg) = refused(&fold, "fold");
    assert_eq!(code, "MOV002", "{msg}");
    let inlet = free.replace(
        "\"amp * sin(pi * x) * sin(4 * pi * y) * sin(freq * t)\"",
        "\"0.01 * sin(freq * t)\"",
    );
    let (code, msg) = refused(&inlet, "inlet");
    assert_eq!(code, "MOV003", "{msg}");
    assert!(msg.contains("xmin"), "{msg}");
    let incompressible = text("motion_piston_adiabatic.toml")
        .replace("density = \"variable\"\n", "")
        .replace("thermodynamic_pressure = \"closed-volume\"", "#")
        .replace("end_time = 1.0", "end_time = 0.05");
    let (code, msg) = refused(&incompressible, "closed");
    assert_eq!(code, "MOV004", "{msg}");
}
