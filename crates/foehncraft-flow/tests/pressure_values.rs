//! Pressure values absolute by default, gauge optional (D239),
//! docs/requirements.md D239 and the codes PRS001–PRS005: an absolute case
//! and its gauge twin run the same bits, older files keep meaning gauge, and
//! every fluid part takes its operating pressure from its own outlet.

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::diag;
use foehncraft_flow::pressure_values::PressureValues;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "foehncraft_pressure_values_{name}_{}",
        std::process::id()
    ))
}

/// A laminar channel: velocity inlet, static-pressure outlet at `outlet`.
fn channel(head: &str, outlet: &str) -> String {
    format!(
        r#"{head}
[mesh.generate]
type = "hex"
dimension = 2
upper = [4.0, 1.0, 0.1]
cells = [24, 8, 1]

[fluid]
density = 1.0
viscosity = 0.05

[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
pressure = {outlet}
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"

[solver]
max_iterations = 60
tolerance = 1e-12
"#
    )
}

fn read(dir: &Path, f: &str) -> String {
    std::fs::read_to_string(dir.join(f)).unwrap()
}

/// The columns of `boundaries.csv` after its comment lines, per boundary.
fn boundary_rows(dir: &Path) -> Vec<Vec<String>> {
    read(dir, "boundaries.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .skip(1)
        .map(|l| l.split(',').map(str::to_string).collect())
        .collect()
}

/// The residual history without its comment lines (and without times).
fn residuals(dir: &Path) -> String {
    read(dir, "residuals.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect::<Vec<_>>()
        .join("\n")
}

/// An absolute case at 1 bar and its gauge twin (operating pressure 1 bar,
/// outlet 0 Pa) give the same bits: the residual history and the mass flows
/// are equal, and the boundary pressures differ by the operating pressure.
#[test]
fn an_absolute_case_and_its_gauge_twin_run_the_same_bits() {
    let abs = FlowCase::from_toml(&channel("case_format = \"0.7\"", "1.0e5"), &cases()).unwrap();
    let gauge = FlowCase::from_toml_with(
        &channel("case_format = \"0.7\"", "0.0"),
        &cases(),
        &[
            "fluid.pressure_values=\"gauge\"".into(),
            "fluid.operating_pressure=1.0e5".into(),
        ],
    )
    .unwrap();
    let fa = abs.pressure_frame.as_deref().unwrap();
    assert_eq!(fa.pressure_values, PressureValues::Absolute);
    assert_eq!(fa.parts.len(), 1);
    assert_eq!(fa.parts[0].operating_pressure, 1.0e5);
    assert_eq!(fa.parts[0].source, "boundary xmax");
    // the reader made the gauge twin
    assert_eq!(abs.fluid.pressure_values, Some(PressureValues::Gauge));
    assert_eq!(abs.fluid.operating_pressure, Some(1.0e5));
    let (da, dg) = (out("abs"), out("gauge"));
    run(&abs, "abs", &da, |_| {}).unwrap();
    run(&gauge, "gauge", &dg, |_| {}).unwrap();
    assert_eq!(residuals(&da), residuals(&dg), "bitwise the same history");
    let (ra, rg) = (boundary_rows(&da), boundary_rows(&dg));
    for (a, g) in ra.iter().zip(&rg) {
        assert_eq!(a[0], g[0]);
        assert_eq!(a[3], g[3], "mass flow of {}", a[0]);
        let (pa, pg): (f64, f64) = (a[4].parse().unwrap(), g[4].parse().unwrap());
        assert!(
            (pa - 1.0e5 - pg).abs() <= 1e-9 * 1.0e5,
            "{}: {pa} vs {pg}",
            a[0]
        );
    }
    // the report states the convention and the operating pressure
    let ja: serde_json::Value = serde_json::from_str(&read(&da, "flow.json")).unwrap();
    assert_eq!(ja["pressure_values"]["pressure_values"], "absolute");
    assert_eq!(
        ja["pressure_values"]["parts"][0]["operating_pressure"],
        1.0e5
    );
    let log = read(&da, "flow.log");
    assert!(
        log.contains("[PRS005] pressure values are absolute"),
        "{log}"
    );
    assert!(
        log.contains("operating pressure 100000 Pa (boundary xmax)"),
        "{log}"
    );
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&dg);
}

/// A file of format 0.6 (or without `case_format`) means gauge values: the
/// migration puts `pressure_values = "gauge"` in, with a note, and the
/// operating pressure keeps its default.
#[test]
fn older_files_keep_meaning_gauge() {
    for head in ["", "case_format = \"0.6\""] {
        let c = FlowCase::from_toml(&channel(head, "0.0"), &cases()).unwrap();
        let f = c.pressure_frame.as_deref().unwrap();
        assert_eq!(f.pressure_values, PressureValues::Gauge, "{head}");
        assert_eq!(f.parts[0].operating_pressure, 101_325.0);
        assert_eq!(c.boundary_pressure("xmax"), Some(0.0));
        assert!(
            c.notes.iter().any(|n| n.contains("pressure_values")),
            "{:?}",
            c.notes
        );
    }
}

/// PRS001: an absolute pressure of 0 Pa (a gauge value in a new file).
#[test]
fn an_absolute_zero_is_refused() {
    let e = FlowCase::from_toml(&channel("case_format = \"0.7\"", "0.0"), &cases()).unwrap_err();
    assert_eq!(diag::of_error(&e).code, "PRS001", "{e:#}");
    let m = format!("{e:#}");
    assert!(m.contains("pressure_values = \"gauge\""), "{m}");
}

/// Expressions and a total-pressure inlet: the operating pressure is the
/// outlet's, the inlet becomes its gauge difference exactly (Sterbenz).
#[test]
fn values_become_gauge_relative_to_the_outlet() {
    let text = channel("case_format = \"0.7\"", "\"1.5e5 + 0*y\"").replace(
        "[boundary.xmin]\ntype = \"velocity-inlet\"\nvelocity = [1.0, 0.0, 0.0]",
        "[boundary.xmin]\ntype = \"total-pressure-inlet\"\ntotal_pressure = 150010.0",
    );
    let c = FlowCase::from_toml(&text, &cases()).unwrap();
    let f = c.pressure_frame.as_deref().unwrap();
    // the outlet is an expression: the inlet's number gives p_op
    assert_eq!(f.parts[0].operating_pressure, 150_010.0);
    assert_eq!(c.boundary_pressure("xmin"), Some(0.0));
    // a number outlet ranks first
    let text = channel("case_format = \"0.7\"", "1.5e5").replace(
        "[boundary.xmin]\ntype = \"velocity-inlet\"\nvelocity = [1.0, 0.0, 0.0]",
        "[boundary.xmin]\ntype = \"total-pressure-inlet\"\ntotal_pressure = 150010.0",
    );
    let c = FlowCase::from_toml(&text, &cases()).unwrap();
    assert_eq!(c.boundary_pressure("xmin"), Some(10.0));
    assert_eq!(c.boundary_pressure("xmax"), Some(0.0));
    // only expressions: PRS003
    let text = channel("case_format = \"0.7\"", "\"1.5e5 + 0*y\"");
    let e = FlowCase::from_toml(&text, &cases()).unwrap_err();
    assert_eq!(diag::of_error(&e).code, "PRS003", "{e:#}");
    // with an explicit operating pressure the expression is fine
    let c = FlowCase::from_toml_with(&text, &cases(), &["fluid.operating_pressure=1.5e5".into()])
        .unwrap();
    assert_eq!(
        c.pressure_frame.as_deref().unwrap().parts[0].source,
        "[fluid] operating_pressure"
    );
}

/// V7.45 (absolute and gauge parts): two fluid parts across a solid plate
/// with outlets at 1.0·10⁵ Pa and 1.5·10⁵ Pa. Each part takes its own
/// outlet as its operating pressure and starts there (gauge 0 in every
/// cell); the gauge twin with the region operating pressure gives the same
/// bits.
#[test]
fn two_parts_take_their_own_outlet_pressures() {
    let base = std::fs::read_to_string(cases().join("cht_counterflow_two_fluids.toml"))
        .unwrap()
        .replace("\r\n", "\n")
        .replace("case_format = \"0.6\"", "case_format = \"0.7\"")
        .replace("max_iterations = 3000", "max_iterations = 40");
    let abs = base
        .replace(
            "[boundary.xmax]\ntype = \"static-pressure-outlet\"\npressure = 0.0",
            "[boundary.xmax]\ntype = \"static-pressure-outlet\"\npressure = 1.0e5",
        )
        .replace(
            "[boundary.oil_out]\ntype = \"static-pressure-outlet\"\npressure = 0.0",
            "[boundary.oil_out]\ntype = \"static-pressure-outlet\"\npressure = 1.5e5",
        );
    let a = FlowCase::from_toml(&abs, &cases()).unwrap();
    let f = a.pressure_frame.as_deref().unwrap();
    assert_eq!(f.parts.len(), 2, "{f:?}");
    let op: Vec<(Vec<String>, f64)> = f
        .parts
        .iter()
        .map(|p| (p.regions.clone(), p.operating_pressure))
        .collect();
    assert!(op.contains(&(vec!["water".into()], 1.0e5)), "{op:?}");
    assert!(op.contains(&(vec!["oil".into()], 1.5e5)), "{op:?}");
    assert_eq!(a.regions["oil"].operating_pressure, Some(1.5e5));
    assert_eq!(a.boundary_pressure("oil_out"), Some(0.0));
    assert_eq!(a.boundary_pressure("xmax"), Some(0.0));
    let gauge = base.replace(
        "[fluid]\nmaterial = \"water\"",
        "[fluid]\npressure_values = \"gauge\"\noperating_pressure = 1.0e5\nmaterial = \"water\"",
    ).replace(
        "[regions.oil]\nkind = \"fluid\"",
        "[regions.oil]\nkind = \"fluid\"\noperating_pressure = 1.5e5",
    );
    let g = FlowCase::from_toml(&gauge, &cases()).unwrap();
    let (da, dg) = (out("parts_abs"), out("parts_gauge"));
    run(&a, "abs", &da, |_| {}).unwrap();
    run(&g, "gauge", &dg, |_| {}).unwrap();
    assert_eq!(residuals(&da), residuals(&dg), "bitwise the same history");
    // the oil outlet reports 1.5 bar, the water outlet 1 bar
    let rows = boundary_rows(&da);
    let p = |n: &str| -> f64 { rows.iter().find(|r| r[0] == n).unwrap()[4].parse().unwrap() };
    assert!((p("oil_out") - 1.5e5).abs() < 1.0, "{}", p("oil_out"));
    assert!((p("xmax") - 1.0e5).abs() < 1.0, "{}", p("xmax"));
    // a case-wide absolute start pressure for both parts: PRS004
    let e = FlowCase::from_toml(
        &abs.replace("[energy]", "[initial]\npressure = 1.2e5\n\n[energy]"),
        &cases(),
    )
    .unwrap_err();
    assert_eq!(diag::of_error(&e).code, "PRS004", "{e:#}");
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&dg);
}
