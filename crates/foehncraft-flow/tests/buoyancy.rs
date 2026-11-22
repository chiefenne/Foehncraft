//! Stage 7, PR 7b1: gravity, the Boussinesq buoyancy, the reduced
//! pressure and the balanced-force interpolation (D153, D154),
//! docs/requirements/stage7-energy-buoyancy-cht.md, checks V7.7 (a) and the
//! buoyant slot.

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_buoyancy_{name}_{}", std::process::id()))
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

/// Runs a case and returns its report and the largest cell speed of the
/// final state (from the restart file, which needs the `hdf5` feature).
#[cfg(feature = "hdf5")]
fn run_speed(case: &FlowCase, name: &str) -> (foehncraft_flow::case::FlowReport, f64) {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let r = foehncraft_flow::restart::read(&dir.join("restart.h5")).unwrap();
    let u = &r.state.u;
    let s = (0..u.x.len())
        .map(|c| (u.x[c] * u.x[c] + u.y[c] * u.y[c] + u.z[c] * u.z[c]).sqrt())
        .fold(0.0, f64::max);
    let _ = std::fs::remove_dir_all(&dir);
    (rep, s)
}

/// The code and message a case text is refused with.
fn refused(text: &str) -> (String, String) {
    let e = match FlowCase::from_toml(text, &cases()) {
        Err(e) => e,
        Ok(c) => match run(&c, "refused", &out("refused"), |_| {}) {
            Err(e) => e,
            Ok(_) => panic!("the case was accepted"),
        },
    };
    (diag::of_error(&e).code.to_string(), format!("{e:#}"))
}

/// The stratified box of `buoyancy_box_at_rest.toml` on mesh type `kind`
/// (with `dimension`, `cells`) and the algorithm `alg`.
#[cfg(feature = "hdf5")]
fn box_case(kind: &str, dimension: u32, cells: &str, alg: &str, iterations: usize) -> FlowCase {
    let t = text("buoyancy_box_at_rest.toml")
        .replace("type = \"hex\"", &format!("type = \"{kind}\""))
        .replace("dimension = 2", &format!("dimension = {dimension}"))
        .replace("cells = [16, 16, 1]", &format!("cells = {cells}"))
        .replace(
            "max_iterations = 1000",
            &format!("max_iterations = {iterations}\nalgorithm = \"{alg}\""),
        );
    let t = if dimension == 3 {
        let wall = "type = \"wall\"\n[boundary.zmin.energy]\ntemperature = \"T\"\n";
        t.replace(
            "[boundary.zmin]\ntype = \"symmetry\"\n[boundary.zmax]\ntype = \"symmetry\"\n",
            &format!(
                "[boundary.zmin]\n{wall}[boundary.zmax]\n{}",
                wall.replace("zmin", "zmax")
            ),
        )
        .replace(
            "temperature = \"T\"",
            "temperature = \"300.0 + 10.0*(0.6*x + 0.8*y)\"",
        )
        .replace("upper = [1.0, 1.0, 0.1]", "upper = [1.0, 1.0, 1.0]")
        .replace("[-5.886, -7.848, 0.0]", "[-3.924, -5.886, -6.867]")
        .replace(
            "300.0 + 10.0*(0.6*x + 0.8*y)",
            "300.0 + 10.0*(0.4*x + 0.6*y + 0.7*z)",
        )
    } else {
        t
    };
    parse(&t)
}

/// V7.7 (a): a closed box at rest with a stable linear stratification
/// along a gravity oblique to the mesh axes stays at rest to round-off on
/// perturbed meshes of every type, with SIMPLEC and with the coupled
/// algorithm: the balanced force of D154 has a discrete hydrostatic state.
#[test]
#[cfg(feature = "hdf5")]
fn stratified_box_stays_at_rest() {
    let g = 9.81f64;
    let bound = 1e-10 * (g * 1.0f64).sqrt();
    for (kind, dim, cells, alg, its) in [
        ("hex", 2, "[16, 16, 1]", "simplec", 1000),
        ("tet", 2, "[12, 12, 1]", "simplec", 1000),
        ("poly", 2, "[12, 12, 1]", "simplec", 1000),
        ("mixed", 3, "[5, 5, 5]", "simplec", 1000),
        ("hex", 2, "[16, 16, 1]", "coupled", 300),
        ("poly", 2, "[12, 12, 1]", "coupled", 300),
    ] {
        let (_, s) = run_speed(
            &box_case(kind, dim, cells, alg, its),
            &format!("box_{kind}_{alg}"),
        );
        assert!(
            s < bound,
            "{kind} {alg}: largest speed {s:e} m/s, bound {bound:e}"
        );
    }
}

/// The buoyant slot between a hot and a cold vertical wall: the cubic
/// velocity profile v = (gβΔT D²/(48ν))(ξ³ − ξ) to second order on
/// perturbed meshes.
#[test]
fn buoyant_slot_has_second_order() {
    let r = run(
        &parse(&text("buoyancy_slot.toml")),
        "slot",
        &out("slot"),
        |_| {},
    )
    .unwrap();
    let _ = std::fs::remove_dir_all(out("slot"));
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    assert!(r.orders.len() == 2, "{:?}", r.orders);
    for o in &r.orders {
        assert!((1.8..2.3).contains(o), "observed orders {:?}", r.orders);
    }
}

/// The reduced and the static outlet value describe the same flow: a tall
/// channel run once with `pressure_kind = "reduced"` (0) and once with
/// `"static"` (the hydrostatic profile) gives the same velocity.
#[test]
#[cfg(feature = "hdf5")]
fn reduced_and_static_outlet_give_the_same_flow() {
    let t = text("buoyancy_channel_reduced.toml");
    let (a, sa) = run_speed(&parse(&t), "chan_reduced");
    let (b, sb) = run_speed(&parse(&text("buoyancy_channel_static.toml")), "chan_static");
    assert!(a.failures.is_empty() && b.failures.is_empty());
    assert!((sa - sb).abs() < 1e-9 * sa, "{sa} against {sb}");
    let ma = a.solves[0].mass_imbalance;
    assert!(ma.abs() < 1e-10, "{ma}");
}

/// Without the hydrostatic part a uniform static outlet value of 0 Pa on a
/// tall opening drives a flow through a box that should be at rest
/// (BUO006 warns), while the default reduced value keeps it at rest.
#[test]
#[cfg(feature = "hdf5")]
fn static_outlet_without_hydrostatics_drives_a_false_flow() {
    let t = text("buoyancy_open_box.toml");
    let (r, s_reduced) = run_speed(&parse(&t), "open_reduced");
    assert!(s_reduced < 1e-10, "{s_reduced}");
    assert!(!r.findings.iter().any(|f| f.code == "BUO006"));
    let t2 = t.replace(
        "pressure = 0.0",
        "pressure = 0.0\npressure_kind = \"static\"",
    );
    let (r2, s_static) = run_speed(&parse(&t2), "open_static");
    assert!(s_static > 1e-3, "{s_static}");
    assert!(r2.findings.iter().any(|f| f.code == "BUO006"));
}

#[test]
fn buoyancy_needs_gravity_and_energy() {
    let base = text("buoyancy_box_at_rest.toml");
    let (code, msg) =
        refused(&base.replace("vector = [-5.886, -7.848, 0.0]", "vector = [0.0, 0.0, 0.0]"));
    assert_eq!(code, "BUO001", "{msg}");
    let no_g: String = base
        .lines()
        .filter(|l| !l.starts_with("vector") && *l != "[gravity]")
        .collect::<Vec<_>>()
        .join("\n");
    let (code, msg) = refused(&no_g);
    assert_eq!(code, "BUO001", "{msg}");
    let (code, msg) = refused(&base.replace("expansion_coefficient = 1e-3\n", ""));
    assert_eq!(code, "BUO003", "{msg}");
    let (code, msg) =
        refused(&base.replace("model = \"boussinesq\"", "model = \"variable-density\""));
    assert_eq!(code, "ENE006", "{msg}");
}
