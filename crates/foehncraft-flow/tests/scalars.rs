//! Stage 6, part 6a: passive scalars (D51) with the generic scalar
//! transport equation (D116) and the bounded BDF2 (D28 (1)),
//! docs/verification/stage6.md.
//!
//! The case files of the gallery (`channel_heated.toml` on three meshes)
//! also run in CI; the tests here use small meshes or two levels.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_scalar_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

/// A case file with `\n` line ends (a Windows checkout has `\r\n`).
fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
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

/// The laminar heated channel (D51 acceptance): the Nusselt number of the
/// thermally developed flow with uniform wall flux, Nu = 140/17, within
/// 0.5 % on the second mesh, and second order of the temperature error.
/// The case file runs the three meshes in CI.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn heated_channel_nusselt_number() {
    let t = text("channel_heated.toml").replace("levels = [0, 1, 2]", "levels = [0, 1]");
    let r = run_case(&parse(&t), "heated");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let (g, l, q, gamma) = (2e-3, 2.0, 1e-3, 1e-3);
    let nu: Vec<f64> = r
        .solves
        .iter()
        .map(|s| {
            let sc = &s.scalars[0];
            let b = |n: &str| sc.boundaries.iter().find(|b| b.boundary == n).unwrap();
            let tw = 0.5 * (b("ymin").mean + b("ymax").mean);
            let tb = b("xmin").mixed_mean.unwrap() + g * l / 2.0;
            // the wall flux is the prescribed one
            assert!((b("ymin").flux - q).abs() < 1e-12, "{}", b("ymin").flux);
            2.0 * q / (gamma * (tw - tb))
        })
        .collect();
    let exact = 140.0 / 17.0;
    assert!(
        (nu[1] / exact - 1.0).abs() < 5e-3,
        "Nu {nu:?}, exact {exact}"
    );
    assert!((nu[1] - exact).abs() < (nu[0] - exact).abs(), "{nu:?}");
    let o = &r.scalar_orders[0].1;
    assert!(o[0] > 1.8 && o[0] < 2.2, "orders {o:?}");
}

/// A manufactured solution φ = 1 + cos(πx/2) sin(πy) in the Poiseuille
/// flow of the channel, u = 6 y (1 − y), with the source
/// S = u ∂φ/∂x − Γ ∇²φ; its x derivative vanishes at the outlet, so the
/// zero gradient there is exact.
fn mms(mesh: &str, levels: &str) -> FlowCase {
    parse(&format!(
        r#"
[mesh.generate]
type = "{mesh}"
dimension = 2
upper = [2.0, 1.0, 0.1]
cells = [16, 8, 1]
perturb = 0.15
[fluid]
viscosity = 0.1
[constants]
gamma = 0.02
[scalars.c]
diffusivity = "gamma"
exact = "1 + cos(pi * x / 2) * sin(pi * y)"
source = "-3 * pi * y * (1 - y) * sin(pi * x / 2) * sin(pi * y) + gamma * 1.25 * pi^2 * cos(pi * x / 2) * sin(pi * y)"
[boundary.xmin]
type = "velocity-inlet"
profile = "analytic"
[boundary.xmin.scalar.c]
value = "1 + cos(pi * x / 2) * sin(pi * y)"
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "wall"
[boundary.ymin.scalar.c]
value = 1.0
[boundary.ymax]
type = "wall"
[boundary.ymax.scalar.c]
value = 1.0
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
tolerance = 1e-8
max_iterations = 5000
[solver.momentum]
relaxation = 0.9
[analytic]
type = "poiseuille"
mean_velocity = 1.0
x_outlet = 2.0
[verification]
levels = {levels}
require_convergence = true
[output]
vtkhdf = false
wall_csv = false
"#
    ))
}

/// The manufactured solution converges at second order on hex cells and
/// within the band of the velocity on polyhedra and triangular prisms
/// (the non-orthogonal correction and the least-squares gradients of the
/// scalar are those of the momentum).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn manufactured_solution_orders() {
    for (mesh, lo, hi) in [("hex", 1.8, 2.3), ("poly", 1.6, 2.6), ("tet", 1.6, 2.6)] {
        let r = run_case(&mms(mesh, "[0, 1, 2]"), &format!("mms_{mesh}"));
        assert!(r.failures.is_empty(), "{mesh}: {:?}", r.failures);
        let o = &r.scalar_orders[0].1;
        let last = *o.last().unwrap();
        eprintln!("{mesh}: scalar orders {o:?}");
        assert!(last > lo && last < hi, "{mesh}: orders {o:?}");
    }
}

/// One mesh of the manufactured solution in a debug build: the error is
/// small and the scalar is reported.
#[test]
fn manufactured_solution_on_one_mesh() {
    let r = run_case(&mms("hex", "[0]"), "mms_one");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0].scalars[0];
    assert_eq!(s.name, "c");
    let l2 = s.l2.unwrap();
    assert!(l2 > 0.0 && l2 < 2e-2, "L2 {l2:e}");
    assert!(s.min > 0.0 && s.max < 2.05, "{} {}", s.min, s.max);
}

/// A front of a scalar in [0, 1] entering a uniform flow, upwind in space,
/// BDF2 in time.
fn front(bounds: &str) -> FlowCase {
    parse(&format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [2.0, 0.1, 0.1]
cells = [40, 1, 1]
[fluid]
viscosity = 0.01
[flow]
convection = true
[scalars.c]
diffusivity = 1e-6
initial = 0.0
{bounds}
[boundary.xmin]
type = "velocity-inlet"
profile = "analytic"
[boundary.xmin.scalar.c]
value = 1.0
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "slip-wall"
[boundary.ymax]
type = "slip-wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver.scalar]
scheme = "upwind"
[analytic]
type = "uniform"
velocity = [1.0, 0.0, 0.0]
length = 1.0
[time]
scheme = "bdf2"
time_step = 0.04
steps = 25
[output]
vtkhdf = false
wall_csv = false
"#
    ))
}

/// The bounded BDF2 (D28 (1)): with `bounds = [0, 1]` the front stays in
/// its bounds, and the cells where the second-order backward difference
/// would leave them fall back to first order; without bounds BDF2
/// overshoots.
#[test]
fn bounded_bdf2_keeps_the_front_in_its_bounds() {
    let dir = out("front");
    let r = run(&front("bounds = [0.0, 1.0]"), "front", &dir, |_| {}).unwrap();
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0].scalars[0];
    assert!(
        s.min >= -1e-10 && s.max <= 1.0 + 1e-10,
        "{} {}",
        s.min,
        s.max
    );
    // the history counts the cells at first order
    let h = std::fs::read_to_string(dir.join("history.csv")).unwrap();
    let mut lines = h.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = head.iter().position(|c| *c == "first_order_c").unwrap();
    let fallen: usize = lines
        .map(|l| l.split(',').nth(col).unwrap().parse::<usize>().unwrap())
        .sum();
    assert!(fallen > 0, "no cell fell back to first order");
    let _ = std::fs::remove_dir_all(&dir);
    let free = run_case(&front(""), "front_free");
    let m = free.solves[0].scalars[0].max;
    assert!(m > 1.0 + 1e-6, "BDF2 without bounds: max {m}");
}

/// The Taylor–Green vortex carrying a smooth scalar: second order in time
/// with BDF2 and with Crank–Nicolson (whose scalars take BDF2), from the
/// differences of successive solutions of a time-step sweep.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn scalar_time_order() {
    for scheme in ["bdf2", "crank-nicolson"] {
        let c = parse(&format!(
            r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.03125]
cells = [32, 32, 1]
[fluid]
viscosity = 0.01
[scalars.c]
diffusivity = 0.01
initial = "2 + cos(2 * pi * x) * cos(2 * pi * y)"
[boundary.xmin]
type = "velocity-inlet"
profile = "analytic"
[boundary.xmin.scalar.c]
value = "2 + cos(2 * pi * x) * cos(2 * pi * y)"
[boundary.xmax]
type = "velocity-inlet"
profile = "analytic"
[boundary.xmax.scalar.c]
value = "2 + cos(2 * pi * x) * cos(2 * pi * y)"
[boundary.ymin]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymin.scalar.c]
value = "2 + cos(2 * pi * x) * cos(2 * pi * y)"
[boundary.ymax]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymax.scalar.c]
value = "2 + cos(2 * pi * x) * cos(2 * pi * y)"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[analytic]
type = "taylor-green"
[time]
scheme = "{scheme}"
time_step = 0.04
end_time = 0.2
[verification]
time_steps = [0.04, 0.02, 0.01, 0.005]
scalar_order = [1.8, 2.6]
[output]
vtkhdf = false
wall_csv = false
"#
        ));
        let r = run_case(&c, &format!("tg_{scheme}"));
        assert!(r.failures.is_empty(), "{scheme}: {:?}", r.failures);
        let o = &r
            .scalar_orders
            .iter()
            .find(|(n, _)| n == "c (time)")
            .unwrap()
            .1;
        eprintln!("{scheme}: scalar time orders {o:?}");
        assert!(!o.is_empty(), "{scheme}");
    }
}

/// Wrong scalar tables are refused with their codes.
#[test]
fn bad_scalar_tables_are_refused() {
    let t = text("channel_heated.toml");
    // a boundary table of an unknown scalar
    let (code, msg) =
        refused(&t.replace("[boundary.ymin.scalar.theta]", "[boundary.ymin.scalar.phi]"));
    assert_eq!(code, "SCAL001", "{msg}");
    // both a value and a flux
    let (code, msg) = refused(&t.replace(
        "flux = \"q\"\n[boundary.ymax]",
        "flux = \"q\"\nvalue = 0.0\n[boundary.ymax]",
    ));
    assert_eq!(code, "SCAL001", "{msg}");
    // a flux at an inlet
    let (code, msg) = refused(&t.replace(
        "value = \"h * (2 * (y / h)^3 - (y / h)^4 - y / h)\"",
        "flux = 1.0",
    ));
    assert_eq!(code, "SCAL001", "{msg}");
    // anything at an outlet
    let (code, msg) = refused(&t.replace(
        "[boundary.xmax]\ntype = \"static-pressure-outlet\"\npressure = 0.0\n",
        "[boundary.xmax]\ntype = \"static-pressure-outlet\"\npressure = 0.0\n[boundary.xmax.scalar.theta]\nvalue = 0.0\n",
    ));
    assert_eq!(code, "SCAL001", "{msg}");
    // a reserved name (the temperature of the energy equation)
    let (code, msg) = refused(&t.replace("[scalars.theta]", "[scalars.T]"));
    assert_eq!(code, "CASE003", "{msg}");
    // a negative diffusivity
    let (code, msg) = refused(&t.replace("diffusivity = 1e-3", "diffusivity = -1e-3"));
    assert_eq!(code, "CASE004", "{msg}");
}

/// The resolved case (`case.resolved.toml`) keeps the boundary tables of
/// the scalars, so it reads back to the same case.
#[test]
fn resolved_case_keeps_the_scalar_tables() {
    let c = parse(&text("channel_heated.toml"));
    let r = c.resolved_toml().unwrap();
    assert!(r.contains("[boundary.ymin.scalar.theta]"), "{r}");
    let back = FlowCase::from_toml(&r, &cases()).unwrap();
    assert_eq!(back.boundary_scalar, c.boundary_scalar);
}

/// A transient run with a scalar continued from a restart file ends
/// bit-identically to the run without the interruption (D42): the scalar
/// and its old level are part of the restart file.
#[cfg(feature = "hdf5")]
#[test]
fn scalar_restart_continues_bit_identically() {
    use foehncraft_flow::restart;
    let case = |extra: &str| {
        parse(&format!(
            r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [12, 12, 1]
[fluid]
viscosity = 0.01
[scalars.c]
diffusivity = 1e-3
initial = "x"
bounds = [0.0, 1.0]
[boundary.xmin]
type = "wall"
[boundary.xmin.scalar.c]
value = 0.0
[boundary.xmax]
type = "wall"
[boundary.xmax.scalar.c]
value = 1.0
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[output]
vtkhdf = false
wall_csv = false
[time]
scheme = "bdf2"
time_step = 0.02
steps = 16
adaptive = true
max_cfl = 2.0
{extra}
"#
        ))
    };
    let da = out("restart_full");
    let ra = run(&case("[run]\nrestart_at = [8]"), "full", &da, |_| {}).unwrap();
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let from = da.join("restart_000008.h5");
    let db = out("restart_cont");
    let rb = run(
        &case(&format!("[run]\nrestart = '{}'", from.display())),
        "cont",
        &db,
        |_| {},
    )
    .unwrap();
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert_eq!(a.state.scalars.len(), 1);
    assert!(a.state.scalars[0].old.is_some());
    assert_eq!(
        a.state.scalars, b.state.scalars,
        "the continued scalar differs"
    );
    assert_eq!(a.state.p, b.state.p);
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// The bounded convection of the scalars (D237) on the oblique step of
/// `scalar_oblique_step_mlp.toml` (perturbed polyhedra): the face-neighbour
/// limiter leaves [0, 1] by almost 2 %, the MLP limiter keeps the
/// overshoot within its threshold ε = 0.01 (c_max − c_min) and with ε = 0
/// the step is strictly bounded; all converge.
#[test]
fn mlp_limiter_bounds_the_oblique_step() {
    let base = text("scalar_oblique_step_mlp.toml").replace("vtkhdf = true", "vtkhdf = false");
    let range = |t: &str, name: &str| {
        let r = run_case(&parse(t), name);
        assert!(r.failures.is_empty(), "{name}: {:?}", r.failures);
        let s = &r.solves[0];
        assert!(s.converged, "{name}: {} iterations", s.iterations);
        (s.scalars[0].min, s.scalars[0].max)
    };
    let (lo, hi) = range(&base, "mlp");
    assert!(lo > -0.005 && hi < 1.005, "MLP R3: [{lo}, {hi}]");
    let strict = base.replace("limiter_epsilon = 0.01", "limiter_epsilon = 0.0");
    let (lo, hi) = range(&strict, "mlp_strict");
    assert!(lo > -1e-9 && hi < 1.0 + 1e-9, "MLP R3, ε = 0: [{lo}, {hi}]");
    let face = base.replace("limiter = \"mlp\"", "limiter = \"venkatakrishnan\"");
    let (lo, hi) = range(&face, "venkatakrishnan");
    assert!(hi > 1.01 && lo < -0.01, "Venkatakrishnan: [{lo}, {hi}]");
}

/// `limiter = "mlp"` is for the scalars and the energy equation; the
/// momentum refuses it with its code.
#[test]
fn mlp_limiter_is_refused_for_the_momentum() {
    let t = text("scalar_oblique_step_mlp.toml").replace(
        "[solver.scalar]",
        "[solver.momentum]\nlimiter = \"mlp\"\n\n[solver.scalar]",
    );
    let (code, msg) = refused(&t);
    assert_eq!(code, "CASE004", "{msg}");
    assert!(msg.contains("scalars and energy only"), "{msg}");
}
