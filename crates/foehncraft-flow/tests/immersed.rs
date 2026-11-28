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

/// A case file of the catalogue with `\n` line ends (a Windows checkout
/// has `\r\n`, and the tests edit the text line by line).
fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
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

/// The channel of [`channel`] with the two walls cut out of the mesh (L13
/// stage 2, D227): the walls become wall boundaries `lower` and `upper`.
fn cut_channel(ny: usize) -> String {
    channel(ny, 1e-8, "").replace(
        "[immersed]\npermeability = 1e-8\n",
        "[immersed]\nmethod = \"cut\"\n",
    ) + "[boundary.lower]\ntype = \"wall\"\n[boundary.upper]\ntype = \"wall\"\n"
}

/// Plane Poiseuille flow between two cut walls off the cell faces: the
/// cut cells give a sharp wall, the mass flow converges at second order
/// (the penalised walls of the same channel: first order), and the wall
/// forces carry the driving force of the fluid volume exactly.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn cut_poiseuille_converges_at_second_order() {
    let mut err = Vec::new();
    for ny in [16, 32, 64] {
        let r = run_case(&parse(&cut_channel(ny)), &format!("cutpois{ny}"));
        err.push((mass_flow(&r) - EXACT).abs() / EXACT);
    }
    let order = (err[1] / err[2]).log2();
    eprintln!("cut Poiseuille mass flow errors {err:?}, order {order:.2}");
    assert!(err[2] < 0.01, "{err:?}");
    assert!(order > 1.5, "order {order}: {err:?}");
}

#[test]
fn cut_settings_are_refused() {
    let base = cut_channel(8);
    // a cut is not a penalisation: no permeability needed, but the key
    // table needs the method
    assert_eq!(
        refused(&base.replace(
            "method = \"cut\"\n",
            "[immersed.cut]\nmin_volume_fraction = 0.2\n"
        )),
        "IMM101"
    );
    assert_eq!(
        refused(&base.replace(
            "method = \"cut\"\n",
            "method = \"cut\"\n[immersed.cut]\nmin_volume_fraction = 0.95\n"
        )),
        "IMM201"
    );
    let regions = base.replace(
        "[immersed]",
        "[regions.solid]\nmaterial = \"steel\"\n\n[immersed]",
    );
    assert_eq!(refused(&regions), "IMM202");
}

/// K = F/(6π μ a U) of the sphere array (a = 0.2 m, μ = 1 Pa s, ρ = 1,
/// L = 1 m, F = G L³ = 1 N): U is the periodic mass flow.
fn array_drag(r: &FlowReport) -> f64 {
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0];
    assert!(s.converged);
    1.0 / (6.0 * std::f64::consts::PI * 0.2 * s.periodic[0].mass_flow)
}

/// The series of Hasimoto (1959) and Sangani & Acrivos (1982) at c = 0.0335.
const ARRAY_K: f64 = 2.1518;

/// The STL sphere of the array cut out of the mesh (L13 stage 2, D227):
/// the drag coefficient approaches the series from below at second order
/// (measured −9.0, −3.6, −1.9 and −0.9 % on 16³, 24³, 32³ and 48³; the
/// penalised sphere: +6.2, +5.5, +3.9 %), and the wall carries the
/// driving force G L³ to the error of the cut.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn stl_sphere_array_cut_converges_to_the_series() {
    let t = text("immersed_cut_stl_sphere_array.toml");
    let mut errs = Vec::new();
    for n in [16usize, 24] {
        let case = t
            .replace("cells = [32, 32, 32]", &format!("cells = [{n}, {n}, {n}]"))
            .replace("wall_csv = false", "wall_csv = false\nvtkhdf = false");
        let r = run_case(&parse(&case), &format!("stlsphere{n}"));
        errs.push(array_drag(&r) / ARRAY_K - 1.0);
        let wall = r.solves[0]
            .boundaries
            .iter()
            .find(|b| b.name == "sphere")
            .expect("the wall of the cut body");
        assert!((wall.force[0] - 1.0).abs() < 0.06, "{n}: {:?}", wall.force);
    }
    assert!(errs[0] < 0.0 && errs[1] < 0.0, "{errs:?}");
    assert!(errs[1].abs() < 0.045, "{errs:?}");
    assert!(errs[1].abs() < 0.5 * errs[0].abs(), "{errs:?}");
}

/// The same STL sphere as a penalised body: its solid fractions are the
/// exact cut volumes of the triangles, so the drag is that of the
/// primitive sphere (whose fractions are sampled) to well within the
/// penalisation error.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn penalised_stl_sphere_matches_the_primitive() {
    let t = text("immersed_stokes_sphere.toml")
        .replace("cells = [32, 32, 32]", "cells = [16, 16, 16]")
        .replace("permeability = 9.765625e-6", "permeability = 3.90625e-5")
        .replace("derived = [\"solid-fraction\"]", "vtkhdf = false");
    let stl = t.replace(
        "shape = \"sphere\"\ncentre = [0.5, 0.5, 0.5]\nradius = 0.2\n",
        "shape = \"stl\"\nfile = \"immersed_sphere.stl\"\n",
    );
    assert_ne!(stl, t);
    let kp = array_drag(&run_case(&parse(&t), "penprim"));
    let ks = array_drag(&run_case(&parse(&stl), "penstl"));
    assert!((ks / kp - 1.0).abs() < 0.01, "STL {ks}, primitive {kp}");
    assert!(ks > ARRAY_K, "the penalised wall is too thick: {ks}");
}

#[test]
fn stl_body_settings_are_refused() {
    let base = text("immersed_cut_stl_sphere_array.toml");
    let file = "file = \"immersed_sphere.stl\"\n";
    assert_eq!(refused(&base.replace(file, "")), "IMM101");
    assert_eq!(
        refused(&base.replace(file, &format!("{file}scale = 0.0\n"))),
        "IMM101"
    );
    assert_eq!(
        refused(&base.replace("immersed_sphere.stl", "no_such_body.stl")),
        "CASE006"
    );
    assert_eq!(
        refused(&base.replace(
            "shape = \"stl\"\n",
            "shape = \"sphere\"\ncentre = [0.5, 0.5, 0.5]\nradius = 0.2\n"
        )),
        "IMM101"
    );
    // a file without triangles, or with a corner that is not a number,
    // cannot be read (the reader refuses it before the surface check)
    let dir = out("badstl");
    std::fs::create_dir_all(&dir).unwrap();
    let with = |name: &str, text: &str| {
        let path = dir.join(name);
        std::fs::write(&path, text).unwrap();
        refused(&base.replace(
            "immersed_sphere.stl",
            &path.display().to_string().replace('\\', "/"),
        ))
    };
    let empty = with("empty.stl", "solid empty\nendsolid empty\n");
    let nan = with(
        "nan.stl",
        "solid nan\nfacet normal 0 0 1\nouter loop\nvertex nan 0 0\nvertex 1 0 0\n\
         vertex 0 1 0\nendloop\nendfacet\nendsolid nan\n",
    );
    let _ = std::fs::remove_dir_all(&dir);
    assert_eq!(empty, "CASE006");
    assert_eq!(nan, "CASE006");
}

/// The tilted turbulent plate of `immersed_cut_plate_keps.toml` on cells
/// twice as large (2 mm at the wall, 75 × 65), with and without the sample
/// cell of the wall function (L13 stage 2c): the cut cells along the plate
/// take every fraction, so their own velocities scatter the wall shear
/// stress; the sampled wall function stays close to White's correlation
/// and scatters much less. Returns the mean of C_f/C_f,White − 1 over
/// 0.8 < s < 1.15 and the relative standard deviation of C_f there.
fn cut_plate(sample: bool) -> (f64, f64) {
    let text = text("immersed_cut_plate_keps.toml")
        .replace("cells = [150, 130, 1]", "cells = [75, 115, 1]")
        .replace(
            "{ to = 1.0, cells = 30, ratio = 142.27 }",
            "{ to = 1.0, cells = 15, ratio = 302.14 }",
        )
        .replace("max_iterations = 2000", "max_iterations = 1500")
        .replace("vtkhdf = true", "vtkhdf = false")
        .replace(
            "method = \"cut\"\n",
            &format!("method = \"cut\"\n[immersed.cut]\nwall_sample = {sample}\n"),
        );
    let case = parse(&text);
    let name = format!("cutplate{sample}");
    let dir = out(&name);
    let rep = run(&case, &name, &dir, |_| {}).unwrap();
    // the k residual stalls near 1e-3 at the leading edge (the inlet), the
    // skin friction downstream is settled: judged by C_f, not by residuals
    assert!(
        rep.failures.iter().all(|f| f.contains("not converged")),
        "{:?}",
        rep.failures
    );
    let csv = std::fs::read_to_string(dir.join("wall_plate.csv")).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let mut lines = csv.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let (cx, cu) = (col("x"), col("u_tau"));
    let cos = 2f64.to_radians().cos();
    let (mut rel, mut cf) = (Vec::new(), Vec::new());
    for l in lines {
        let v: Vec<f64> = l.split(',').map(|x| x.parse().unwrap()).collect();
        let s = v[cx] / cos;
        if !(0.8..1.15).contains(&s) {
            continue;
        }
        let c = 2.0 * v[cu] * v[cu];
        let white = 0.455 / (0.06 * s / 2e-7f64).ln().powi(2);
        rel.push(c / white - 1.0);
        cf.push(c);
    }
    let n = rel.len() as f64;
    let mean = rel.iter().sum::<f64>() / n;
    let m = cf.iter().sum::<f64>() / n;
    let sd = (cf.iter().map(|c| (c - m) * (c - m)).sum::<f64>() / n).sqrt() / m;
    (mean, sd)
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn cut_plate_wall_function_samples_one_layer_out() {
    let (mean_s, sd_s) = cut_plate(true);
    let (mean_c, sd_c) = cut_plate(false);
    eprintln!(
        "tilted cut plate: sampled C_f {:+.1} % (scatter {:.1} %), cut cell itself {:+.1} % \
         (scatter {:.1} %)",
        100.0 * mean_s,
        100.0 * sd_s,
        100.0 * mean_c,
        100.0 * sd_c
    );
    assert!(mean_s.abs() < 0.08, "{mean_s}");
    assert!(sd_s < 0.6 * sd_c, "{sd_s} {sd_c}");
}
