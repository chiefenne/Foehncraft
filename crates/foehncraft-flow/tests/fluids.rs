//! Several fluids and materials in one case (D194, stage MF-A),
//! docs/requirements/stage7-energy-buoyancy-cht.md §9c, check V7.33 and the
//! codes REG003, REG004 and REG006.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_fluids_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn text() -> String {
    std::fs::read_to_string(cases().join("cht_counterflow_two_fluids.toml"))
        .unwrap()
        .replace("\r\n", "\n")
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

/// V7.33 (planar form): water and a light oil in counterflow on the two
/// sides of a steel plate, C_r = 1. The effectiveness from the water
/// outlet against ε–NTU with the fully developed Nu = 5.385 of each
/// channel; the heat the water gives off is the heat the oil takes up, and
/// both equal the heat through the plate.
#[test]
fn counterflow_of_water_and_oil_follows_eps_ntu() {
    let case = FlowCase::from_toml(&text(), &cases()).unwrap();
    let r = run_case(&case, "counterflow");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0];
    assert!(s.converged);
    let t = s.scalars.iter().find(|x| x.name == "temperature").unwrap();
    let b = |n: &str| t.boundaries.iter().find(|b| b.boundary == n).unwrap();
    let tw_out = b("xmax").mixed_mean.unwrap();
    let to_out = b("oil_out").mixed_mean.unwrap();
    let eps = (340.0 - tw_out) / 40.0;
    // the analytic effectiveness (see the case file)
    let (kw, ko, h, l, t_p, k_p) = (0.6094956547006251, 0.13, 0.002, 0.5, 0.0005, 15.0);
    let (hw, ho) = (5.385 * kw / (2.0 * h), 5.385 * ko / (2.0 * h));
    let u = 1.0 / (1.0 / hw + t_p / k_p + 1.0 / ho);
    let c = 996.5570590375 * 0.0025 * h * 4180.413119939875;
    let ntu = u * l / c;
    let eps_ref = ntu / (1.0 + ntu);
    assert!(
        (eps / eps_ref - 1.0).abs() < 0.02,
        "ε = {eps}, ε–NTU {eps_ref}"
    );
    // C_r = 1: the oil warms as much as the water cools
    let eps_oil = (to_out - 300.0) / 40.0;
    assert!((eps_oil / eps - 1.0).abs() < 2e-3, "{eps_oil} vs {eps}");
    // the ledger: water out = through the plate = oil in
    let q_water = b("xmin").enthalpy_flow.unwrap() + b("xmax").enthalpy_flow.unwrap();
    let q_oil = b("oil_in").enthalpy_flow.unwrap() + b("oil_out").enthalpy_flow.unwrap();
    let q_plate = b("interface-plate").heat_flow.unwrap();
    assert!(q_water > 0.0 && q_oil < 0.0, "{q_water} {q_oil}");
    assert!(
        (q_water + q_oil).abs() < 1e-3 * q_water,
        "{q_water} {q_oil}"
    );
    assert!((q_plate + q_water).abs() < 1e-3 * q_water, "{q_plate} {q_water}");
}

/// REG004: the oil directly on the water (the plate a fluid region of the
/// oil) is a two-phase interface.
#[test]
fn two_fluids_touching_are_refused() {
    let t = text().replace(
        "[regions.plate]\nkind = \"solid\"\nzones = [\"plate\"]\nmaterial = \"plate-steel\"",
        "[regions.plate]\nkind = \"fluid\"\nzones = [\"plate\"]\nmaterial = \"light-oil\"",
    );
    let (code, msg) = refused(&t);
    assert_eq!(code, "REG004", "{msg}");
    assert!(msg.contains("two-phase"), "{msg}");
}

/// REG003: a fluid part without a pressure boundary in a case of two.
#[test]
fn a_part_without_pressure_level_is_refused() {
    let t = text().replace(
        "[boundary.oil_out]\ntype = \"static-pressure-outlet\"\npressure = 0.0",
        "[boundary.oil_out]\ntype = \"wall\"",
    );
    let (code, msg) = refused(&t);
    assert_eq!(code, "REG003", "{msg}");
    assert!(msg.contains("regions oil"), "{msg}");
}

/// REG006: a setting that several fluids do not yet support.
#[test]
fn gravity_with_several_fluids_is_refused() {
    let t = format!("{}\n[gravity]\nvector = [0.0, -9.81, 0.0]\n", text());
    let (code, msg) = refused(&t);
    assert_eq!(code, "REG006", "{msg}");
}
