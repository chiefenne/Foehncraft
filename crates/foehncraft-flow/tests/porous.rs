//! Porous regions (D270, porous plan stage PM1): the Darcy–Forchheimer
//! sink in the axes of the region, against exact pressure drops (check V1
//! of the porous plan), and the checks of the keys (POR001, POR002).

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_porous_{name}_{}", std::process::id()))
}

fn text() -> String {
    std::fs::read_to_string(cases().join("porous_darcy_block.toml"))
        .unwrap()
        .replace("\r\n", "\n")
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn run_case(text: &str, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(&parse(text), name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

/// The inlet pressure over the outlet pressure (0).
fn inlet_pressure(rep: &FlowReport) -> f64 {
    let s = &rep.solves[0];
    assert!(s.converged);
    s.boundaries
        .iter()
        .find(|b| b.name == "xmin")
        .unwrap()
        .mean_pressure
}

/// The code and message a case text is refused with.
fn refused(text: &str) -> (String, String) {
    let e = FlowCase::from_toml(text, &cases()).unwrap_err();
    (diag::of_error(&e).code.to_string(), format!("{e:#}"))
}

const MU: f64 = 1.8e-5;
const RHO: f64 = 1.2;

/// The whole duct is porous: no jump of the resistance, so the pressure is
/// linear and the drop is exact, (μ K⁻¹ U + ½ ρ F U²) over the length
/// 1 m, to the solver tolerance.
#[test]
fn a_duct_filled_with_a_porous_matrix_has_the_exact_pressure_drop() {
    let t = text().replace("upper = [0.6, 0.1, 0.01] }", "upper = [1.0, 0.1, 0.01] }");
    let t = t.replace("lower = [0.4, 0.0, 0.0]", "lower = [0.0, 0.0, 0.0]");
    let exact = MU * 1e6 + 0.5 * RHO * 10.0;
    let p = inlet_pressure(&run_case(&t, "filled"));
    assert!((p - exact).abs() < 1e-8 * exact, "{p} against {exact}");
}

/// Channels across the duct (the first axis along y): the flow along x
/// meets the second principal value; Darcy only, exact.
#[test]
fn channels_across_the_flow_give_the_across_resistance() {
    let t = text()
        .replace("upper = [0.6, 0.1, 0.01] }", "upper = [1.0, 0.1, 0.01] }")
        .replace("lower = [0.4, 0.0, 0.0]", "lower = [0.0, 0.0, 0.0]")
        .replace(
            "viscous = 1.0e6\ninertial = 10.0",
            "viscous = [1.0e4, 1.0e6, 1.0e6]",
        )
        .replace(
            "zones = [\"block\"]",
            "zones = [\"block\"]\naxes = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]",
        );
    let exact = MU * 1e6;
    let p = inlet_pressure(&run_case(&t, "across"));
    assert!((p - exact).abs() < 1e-8 * exact, "{p} against {exact}");
    // and along the channels (the first axis along x) the first value
    let t = t.replace(
        "axes = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]",
        "axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]",
    );
    let exact = MU * 1e4;
    let p = inlet_pressure(&run_case(&t, "along"));
    assert!((p - exact).abs() < 1e-8 * exact, "{p} against {exact}");
}

/// V1: a block in the duct, clear gas before and after it. The drop
/// (μ K⁻¹ U + ½ ρ F U²) L = 4.8 Pa is reached as the mesh is refined: the
/// resistance jumps at the block faces, which the pressure gradient of the
/// neighbouring cells sees over one cell, an error that falls with the
/// cell size (observed order about 2).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn a_block_in_a_duct_converges_to_the_exact_pressure_drop() {
    let exact = (MU * 1e6 + 0.5 * RHO * 10.0) * 0.2;
    let errs: Vec<f64> = [50, 100, 200]
        .iter()
        .map(|n| {
            let t = text().replace("cells = [100, 4, 1]", &format!("cells = [{n}, 4, 1]"));
            (inlet_pressure(&run_case(&t, &format!("block{n}"))) - exact).abs()
        })
        .collect();
    assert!(errs[1] < 5e-3 * exact, "{errs:?}");
    assert!(errs[2] < 1.5e-3 * exact, "{errs:?}");
    for w in errs.windows(2) {
        let order = (w[0] / w[1]).log2();
        assert!(order > 1.5, "observed order {order} ({errs:?})");
    }
}

/// The checks of the keys.
#[test]
fn porous_keys_are_checked() {
    let t = text();
    // no matrix
    let (c, m) = refused(&t.replace(
        "[regions.block.porous]\nviscous = 1.0e6\ninertial = 10.0\nsource = \"chosen values for the exact check\"\n",
        "",
    ));
    assert_eq!(c, "POR001", "{m}");
    // direction-dependent values without axes
    let (c, m) = refused(&t.replace("viscous = 1.0e6", "viscous = [1.0e6, 1.0e8, 1.0e8]"));
    assert_eq!(c, "POR001", "{m}");
    // both honeycomb and resistances
    let (c, m) = refused(&t.replace(
        "inertial = 10.0",
        "inertial = 10.0\nhoneycomb = { channel = 6.9e-3, pitch = 13.6e-3 }",
    ));
    assert_eq!(c, "POR001", "{m}");
    // a matrix in a fluid region
    let (c, m) = refused(&t.replace("kind = \"porous\"", "kind = \"fluid\""));
    assert_eq!(c, "POR001", "{m}");
    // with the energy equation: stage PM2
    let (c, m) = refused(&t.replace(
        "[boundary.xmax]",
        "[boundary.xmin.energy]\ntemperature = 300.0\n[boundary.xmax]",
    ).replace(
        "[fluid]\ndensity = 1.2\nviscosity = 1.8e-5",
        "[fluid]\ndensity = 1.2\nviscosity = 1.8e-5\nspecific_heat = 1005.0\nconductivity = 0.026\n\n[energy]\ninitial = 300.0",
    ));
    assert_eq!(c, "POR002", "{m}");
    // direction-dependent values with axes turned against the mesh
    let (c, m) = refused(
        &t.replace("viscous = 1.0e6", "viscous = [1.0e6, 1.0e8, 1.0e8]")
            .replace(
                "zones = [\"block\"]",
                "zones = [\"block\"]\naxes = [[0.866, 0.5, 0.0], [-0.5, 0.866, 0.0]]",
            ),
    );
    assert_eq!(c, "POR002", "{m}");
    // a honeycomb is accepted with its axes
    parse(
        &t.replace(
            "viscous = 1.0e6\ninertial = 10.0",
            "honeycomb = { channel = 6.9e-3, pitch = 13.6e-3, length = 0.15 }",
        )
        .replace(
            "zones = [\"block\"]",
            "zones = [\"block\"]\naxes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]",
        ),
    );
}
