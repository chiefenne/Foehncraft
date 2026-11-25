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

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file)).unwrap()
}

/// The manufactured solution of the penalised equations (a smooth solid
/// fraction, the sink everywhere): second order on hexahedra and polygons,
/// as without the sink. On the triangular prisms of `type = "tet"` the
/// same solution without the sink does not converge either (a gap of the
/// flow discretisation on that mesh, not of the penalisation; features.toml).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn manufactured_penalised_solution_converges_at_second_order() {
    let t = text("immersed_mms.toml");
    for (kind, band) in [("hex", [1.9, 2.1]), ("poly", [1.9, 2.2])] {
        let rep = run_case(
            &parse(&t.replace("type = \"hex\"", &format!("type = \"{kind}\""))),
            &format!("mms_{kind}"),
        );
        let o = *rep.orders.last().unwrap();
        // measured 2.00 (hex) and 2.05 (poly) on the finest pair
        assert!(
            o >= band[0] && o <= band[1],
            "{kind}: order {o} outside {band:?} ({:?})",
            rep.orders
        );
    }
}

/// A channel with a grey plug (the design region, γ = `start`) and one
/// design update after the flow has converged: Φ and Σ ∂Φ/∂γ of the
/// report.
fn plug(start: f64, stokes: bool) -> (f64, f64) {
    let c = format!(
        r#"
case_format = "0.7"
[mesh.generate]
type = "hex"
dimension = 2
upper = [2.0, 1.0, 0.05]
cells = [40, 20, 1]
[fluid]
density = 1.0
viscosity = {mu}
pressure_values = "gauge"
[flow]
convection = {conv}
[immersed]
permeability = 1e-3
[immersed.topology]
lower = [0.7, 0.2, -1.0]
upper = [1.2, 0.8, 1.0]
fluid_fraction = 0.5
start = {start:?}
every = 1000000
designs = 1
[boundary.xmin]
type = "velocity-inlet"
profile = "parabolic"
max_velocity = 1.0
profile_axes = ["y"]
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
tolerance = 1e-12
max_iterations = 20000
[output]
vtkhdf = false
wall_csv = false
"#,
        mu = if stokes { 1.0 } else { 0.05 },
        conv = !stokes,
    );
    let r = run_case(&parse(&c), &format!("plug{start}{stokes}"));
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let t = r.solves[0].topology.as_ref().expect("topology report");
    (t.first_objective, t.first_gradient)
}

/// The gradient of the topology search (self-adjoint, no adjoint solve)
/// against central differences of converged runs. For Stokes flow it is
/// the exact derivative of the continuous problem; the discrete operator is
/// self-adjoint up to the pressure-stabilising Rhie–Chow term, which
/// vanishes as the mesh is refined: 8e-5 relative on these 40 × 20 cells,
/// 3e-6 on 80 × 40. With convection it is an approximation (Re = 20: 0.5 %).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn topology_gradient_matches_central_differences() {
    let d = 1e-3;
    for (stokes, tol) in [(true, 2e-4), (false, 2e-2)] {
        let (_, g) = plug(0.5, stokes);
        let (fp, _) = plug(0.5 + d, stokes);
        let (fm, _) = plug(0.5 - d, stokes);
        let fd = (fp - fm) / (2.0 * d);
        assert!(g < 0.0, "more fluid must lower the power: {g}");
        assert!(
            (g / fd - 1.0).abs() < tol,
            "stokes {stokes}: gradient {g}, central difference {fd}"
        );
    }
}

/// The pipe bend of Borrvall & Petersson (2003): the search converges to a
/// black-and-white design with the fluid volume held, the power falls by
/// an order of magnitude, the surface is written.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn topology_search_of_a_pipe_bend() {
    let t =
        text("immersed_topology_bend.toml").replace("cells = [60, 60, 1]", "cells = [30, 30, 1]");
    let dir = out("bend");
    let r = run(&parse(&t), "bend", &dir, |_| {}).unwrap();
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0];
    assert!(s.converged);
    let tp = s.topology.as_ref().expect("topology report");
    assert!(tp.converged, "{tp:?}");
    assert!((tp.fluid_fraction - 0.25).abs() < 1e-3, "{tp:?}");
    assert!(tp.grey < 0.12, "{tp:?}");
    assert!(tp.objective < 0.15 * tp.first_objective, "{tp:?}");
    assert!(tp.triangles > 0 && tp.surface_area > 0.0, "{tp:?}");
    let stl = std::fs::read_to_string(dir.join("topology.stl")).unwrap();
    assert!(stl.starts_with("solid topology"));
    assert!(dir.join("topology.csv").exists());
    let _ = std::fs::remove_dir_all(&dir);
}

/// Stokes flow through a simple cubic array of penalised spheres: the drag
/// coefficient F/(6π μ a U) against the series of Hasimoto (1959) and
/// Sangani & Acrivos (1982), 2.152 at c = 0.0335; with K = 0.01 h² 6.2 %
/// high on 16³ cells and 5.5 % on 24³ (3.9 % on the 32³ of the case file).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn stokes_sphere_array_drag_approaches_the_series() {
    let c = 4.0 / 3.0 * std::f64::consts::PI * 0.2f64.powi(3);
    let x = c.cbrt();
    let k =
        1.0 / (1.0 - 1.7601 * x + c - 1.5593 * c * c + 3.9799 * x.powi(8) - 3.0734 * x.powi(10));
    assert!((k - 2.1518).abs() < 1e-3, "{k}");
    let t = text("immersed_stokes_sphere.toml");
    let mut errs = Vec::new();
    for n in [16usize, 24] {
        let h = 1.0 / n as f64;
        let case = t
            .replace("cells = [32, 32, 32]", &format!("cells = [{n}, {n}, {n}]"))
            .replace(
                "permeability = 9.765625e-6",
                &format!("permeability = {:e}", 0.01 * h * h),
            )
            .replace("derived = [\"solid-fraction\"]", "vtkhdf = false");
        let r = run_case(&parse(&case), &format!("sphere{n}"));
        assert!(r.failures.is_empty(), "{:?}", r.failures);
        let s = &r.solves[0];
        assert!(s.converged);
        let u = s.periodic[0].mass_flow; // ρ = 1, L = 1
        let f = s.immersed[0].force[0];
        assert!((f - 1.0).abs() < 1e-5, "the body carries G L³: {f}");
        let kn = f / (6.0 * std::f64::consts::PI * 0.2 * u);
        errs.push(kn / k - 1.0);
    }
    assert!(
        errs[0] > 0.0 && errs[0] < 0.07 && errs[1] < errs[0],
        "{errs:?}"
    );
}
