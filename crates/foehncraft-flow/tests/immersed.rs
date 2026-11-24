//! L13 stage 1: immersed bodies by Brinkman penalisation (D226).
//!
//! Plane Poiseuille flow between two penalised walls that do not lie on
//! cell faces, driven by a pressure gradient through a periodic pair: the
//! mass flow converges to the exact ρ G H³ / (12 μ) per depth as the mesh
//! is refined (about first order: the wall is smeared over one cell); the
//! SIMPLEC and the coupled solver give the same mass flow; the settings that
//! do not go with a penalised body are refused.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_immersed_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

/// The channel between y = 0.13 and y = 0.87 (H = 0.74 m) in a box mesh of
/// height 1 with `ny` cells; the walls are two penalised boxes.
fn channel(ny: usize, k: f64, extra: &str) -> String {
    format!(
        r#"
case_format = "0.7"
[mesh.generate]
type = "hex"
dimension = 2
upper = [0.5, 1.0, 0.1]
cells = [4, {ny}, 1]
[fluid]
density = 1.0
viscosity = 0.1
pressure_values = "gauge"
[immersed]
permeability = {k:e}
[[immersed.body]]
name = "lower"
shape = "box"
lower = [-1.0, -1.0, -1.0]
upper = [2.0, 0.13, 1.0]
[[immersed.body]]
name = "upper"
shape = "box"
lower = [-1.0, 0.87, -1.0]
upper = [2.0, 2.0, 1.0]
[boundary.xmin]
type = "periodic"
partner = "xmax"
pressure_gradient = 1.2
[boundary.xmax]
type = "periodic"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
tolerance = 1e-9
max_iterations = 3000
{extra}
[output]
vtkhdf = false
wall_csv = false
"#
    )
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

/// The exact mass flow ρ G H³ / (12 μ) × depth (kg/s).
const EXACT: f64 = 1.2 * 0.74 * 0.74 * 0.74 / (12.0 * 0.1) * 0.1;

fn mass_flow(r: &FlowReport) -> f64 {
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    r.solves[0].periodic[0].mass_flow
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn penalised_poiseuille_converges_to_the_exact_mass_flow() {
    let mut err = Vec::new();
    for ny in [16, 32, 64] {
        let r = run_case(&parse(&channel(ny, 1e-8, "")), &format!("pois{ny}"));
        let m = mass_flow(&r);
        err.push((m - EXACT).abs() / EXACT);
        // the walls carry the driving force: Σ s u over both bodies equals
        // G times the volume of the box (fluid and solid, at the steady
        // state the solid cells take their share of G V too)
        let f: f64 = r.solves[0].immersed.iter().map(|b| b.force[0]).sum();
        let gv = 1.2 * 0.5 * 1.0 * 0.1;
        assert!((f / gv - 1.0).abs() < 1e-4, "{ny}: force {f}, G V {gv}");
    }
    // measured 20.5, 8.5, 2.3 %: first order or better
    assert!(err[0] < 0.25 && err[1] < 0.10 && err[2] < 0.03, "{err:?}");
    assert!(err[1] < 0.6 * err[0] && err[2] < 0.6 * err[1], "{err:?}");
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn coupled_and_simplec_give_the_same_penalised_flow() {
    // a sink about 10⁵ times the fluid diagonal: the coupled solver
    // converges from rest (with 10⁷, K = 1e-8, its linear solver breaks
    // down; the gap is written in features.toml)
    let s = mass_flow(&run_case(&parse(&channel(32, 1e-6, "")), "simplec"));
    let c = mass_flow(&run_case(
        &parse(&channel(32, 1e-6, "algorithm = \"coupled\"")),
        "coupled",
    ));
    assert!((s / c - 1.0).abs() < 1e-5, "SIMPLEC {s}, coupled {c}");
}

fn refused(text: &str) -> String {
    let e = match FlowCase::from_toml(text, &cases()) {
        Err(e) => e,
        Ok(c) => match run(&c, "refused", &out("refused"), |_| {}) {
            Err(e) => e,
            Ok(_) => panic!("the case was accepted"),
        },
    };
    diag::of_error(&e).code.to_string()
}

#[test]
fn unsupported_settings_are_refused() {
    let base = channel(8, 1e-8, "");
    assert_eq!(
        refused(&base.replace("permeability = 1e-8\n", "")),
        "IMM101"
    );
    assert_eq!(
        refused(&base.replace("upper = [2.0, 0.13, 1.0]", "upper = [2.0, -2.0, 1.0]")),
        "IMM101"
    );
    let turb = base.replace(
        "[immersed]",
        "[turbulence]\nmodel = \"k-omega-sst\"\n\n[immersed]",
    );
    assert_eq!(refused(&turb), "IMM102");
}
