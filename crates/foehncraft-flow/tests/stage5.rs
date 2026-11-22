//! Stage 5 acceptance: transient flow with PISO (docs/verification/stage5.md).
//!
//! The full Taylor–Green sweeps and the cylinder benchmark also run from the
//! case files in CI; the tests here use small meshes.
#![cfg(feature = "hdf5")]

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::restart;
use foehncraft_fv::disc::Soa3;
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_stage5_{name}_{}", std::process::id()))
}

fn run_in(case: &FlowCase, name: &str, dir: &Path) -> FlowReport {
    run(case, name, dir, |_| {}).unwrap()
}

/// The lid-driven cavity on an `n` x `n` hex mesh with a `[time]` table.
fn cavity(n: usize, time: &str, extra: &str) -> FlowCase {
    let text = format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [{n}, {n}, 1]
[fluid]
viscosity = 0.01
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
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
{time}
{extra}
"#
    );
    FlowCase::from_toml(&text, Path::new(".")).unwrap()
}

/// The Taylor–Green vortex on the unit square with the exact velocity on
/// the four sides.
fn taylor_green(n: usize, time: &str, extra: &str) -> FlowCase {
    let text = format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [{n}, {n}, 1]
[fluid]
viscosity = 0.01
[flow]
convection = true
[boundary.xmin]
type = "velocity-inlet"
profile = "analytic"
[boundary.xmax]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymin]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymax]
type = "velocity-inlet"
profile = "analytic"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[analytic]
type = "taylor-green"
[output]
vtkhdf = false
wall_csv = false
[time]
{time}
{extra}
"#
    );
    FlowCase::from_toml(&text, Path::new(".")).unwrap()
}

fn same3(a: &Soa3, b: &Soa3) -> bool {
    a.x == b.x && a.y == b.y && a.z == b.z
}

fn same_restart(a: &restart::Restart, b: &restart::Restart) -> bool {
    let (s, t) = (&a.state, &b.state);
    let fields = same3(&s.u, &t.u)
        && s.p == t.p
        && s.mass_flux == t.mass_flux
        && same3(&s.ub, &t.ub)
        && s.pb == t.pb
        && (0..3).all(|i| same3(&s.grad_u[i], &t.grad_u[i]))
        && same3(&s.grad_p, &t.grad_p)
        && same3(&s.psi, &t.psi);
    let levels = match (&a.transient, &b.transient) {
        (Some((o, c, dt, _)), Some((p, d, du, _))) => {
            o.u_old == p.u_old
                && o.m_old == p.m_old
                && o.x_cur == p.x_cur
                && o.x_old == p.x_old
                && o.g_old == p.g_old
                && o.u_old2 == p.u_old2
                && c == d
                && dt == du
        }
        _ => false,
    };
    fields && levels && a.iteration == b.iteration && a.amg_matrix == b.amg_matrix
}

/// A run continued from a restart file ends bit-identically to the run
/// without the interruption (D42): BDF2 with the adaptive time step, so the
/// clock, the old levels, the extrapolation and the multigrid all carry
/// over.
#[test]
fn restart_continues_bit_identically() {
    let time = "scheme = \"bdf2\"\ntime_step = 0.01\nsteps = 24\nadaptive = true\nmax_cfl = 2.0";
    let full = cavity(
        16,
        time,
        "[run]\nrestart_at = [10]\n[solver]\namg_rebuild = 5",
    );
    let da = out("restart_full");
    let ra = run_in(&full, "full", &da);
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let from = da.join("restart_000010.h5");
    let cont = cavity(
        16,
        time,
        &format!(
            "[run]\nrestart = '{}'\n[solver]\namg_rebuild = 5",
            from.display()
        ),
    );
    let db = out("restart_cont");
    let rb = run_in(&cont, "cont", &db);
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert_eq!(
        a.transient.as_ref().unwrap().1.step,
        24,
        "the full run ends at step 24"
    );
    assert!(same_restart(&a, &b), "the continued run differs");
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// The time step from the error estimate (D65) with rejected steps: the
/// run steps back exactly and continues from a restart file bit-identically
/// (the level uⁿ⁻² of the predictor and the step before the last are part
/// of the restart).
#[test]
fn error_controlled_steps_restart_bit_identically() {
    let time = "scheme = \"bdf2\"\ntime_step = 0.02\nsteps = 30\nadaptive = true\nmax_cfl = 2.0\n\
                error_tolerance = 1e-3";
    let full = cavity(16, time, "[run]\nrestart_at = [12]");
    let da = out("restart_err_full");
    let ra = run_in(&full, "full", &da);
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let log = std::fs::read_to_string(da.join("flow.log")).unwrap();
    assert!(
        log.contains(" rejected: error "),
        "no step was rejected:\n{log}"
    );
    assert!(
        log.contains(" Error "),
        "no error estimate in the log:\n{log}"
    );
    let from = da.join("restart_000012.h5");
    let r = restart::read(&from).unwrap();
    let (o, c, ..) = r.transient.as_ref().unwrap();
    assert!(o.u_old2.is_some() && c.dt_prev2 > 0.0);
    let cont = cavity(16, time, &format!("[run]\nrestart = '{}'", from.display()));
    let db = out("restart_err_cont");
    let rb = run_in(&cont, "cont", &db);
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert!(same_restart(&a, &b), "the continued run differs");
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// Output in simulated time (D66): the adaptive steps land on every
/// scheduled time, so the monitor rows are equally spaced in time, and the
/// numbered restart file sits at its scheduled time; no step is shorter
/// than half the step before it.
#[test]
fn output_lands_on_scheduled_times() {
    let time =
        "scheme = \"bdf2\"\ntime_step = 0.013\nend_time = 0.5\nadaptive = true\nmax_cfl = 1.5";
    let mut case = cavity(
        16,
        time,
        "[run]\nrestart_times = [0.3337]\nrestart_interval = 0.2",
    );
    case.output.monitor_interval = Some(0.05);
    let d = out("scheduled");
    let r = run_in(&case, "scheduled", &d);
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let forces = std::fs::read_to_string(d.join("forces.csv")).unwrap();
    let times: Vec<f64> = forces
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("step"))
        .map(|l| l.split(',').nth(1).unwrap().parse().unwrap())
        .collect();
    assert_eq!(times.len(), 10, "{times:?}");
    for (k, t) in times.iter().enumerate() {
        let want = 0.05 * (k + 1) as f64;
        assert!((t - want).abs() < 1e-12, "row {k} at t = {t}, not {want}");
    }
    let numbered: Vec<_> = std::fs::read_dir(&d)
        .unwrap()
        .filter_map(|e| {
            let n = e.unwrap().file_name().into_string().unwrap();
            n.starts_with("restart_").then_some(n)
        })
        .collect();
    assert_eq!(numbered.len(), 1, "{numbered:?}");
    let rs = restart::read(&d.join(&numbered[0])).unwrap();
    let t = rs.transient.unwrap().1.time;
    assert!((t - 0.3337).abs() < 1e-12, "restart at t = {t}");
    let hist = std::fs::read_to_string(d.join("history.csv")).unwrap();
    let dts: Vec<f64> = hist
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("step"))
        .map(|l| l.split(',').nth(2).unwrap().parse().unwrap())
        .collect();
    for w in dts.windows(2).take(dts.len() - 2) {
        assert!(
            w[1] >= 0.5 * w[0] * (1.0 - 1e-9),
            "step {} after {}",
            w[1],
            w[0]
        );
    }
    let _ = std::fs::remove_dir_all(&d);
}

/// Crank–Nicolson with two outer iterations continues bit-identically too
/// (the explicit operator of the old level and the outer loop).
#[test]
fn restart_crank_nicolson_outer() {
    let time = "scheme = \"crank-nicolson\"\ntime_step = 0.01\nsteps = 12";
    let extra = "[solver]\nouter_iterations = 2";
    let full = taylor_green(16, time, &format!("[run]\nrestart_at = [5]\n{extra}"));
    let da = out("restart_cn_full");
    run_in(&full, "full", &da);
    let from = da.join("restart_000005.h5");
    let cont = taylor_green(
        16,
        time,
        &format!("[run]\nrestart = '{}'\n{extra}", from.display()),
    );
    let db = out("restart_cn_cont");
    run_in(&cont, "cont", &db);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert!(same_restart(&a, &b), "the continued run differs");
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// `.stop` ends the run after the current step with the restart file and
/// the results (D40).
#[test]
fn stop_file_ends_the_run() {
    let case = cavity(
        12,
        "time_step = 0.01\nsteps = 50",
        "[solver]\nlog_every = 1",
    );
    let dir = out("stop");
    std::fs::create_dir_all(&dir).unwrap();
    let stop = dir.join(".stop");
    let rep = run(&case, "stop", &dir, |l: &str| {
        if l.starts_with(" ") && l.split_whitespace().next() == Some("5") && l.contains("E-") {
            std::fs::write(&stop, "").unwrap();
        }
    })
    .unwrap();
    assert!(rep.stopped);
    let r = restart::read(&dir.join("restart.h5")).unwrap();
    assert_eq!(r.transient.unwrap().1.step, 5);
    // a start with an existing .stop fails at once
    let e = run(&case, "stop", &dir, |_| {}).unwrap_err();
    assert!(format!("{e:#}").contains(".stop"), "{e:#}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// `.result_N` writes the fields of step N; `.read_toml_N` re-reads the
/// case file at step N and applies only what may change during a run.
#[test]
fn result_and_reread_files() {
    let dir = out("reread");
    std::fs::create_dir_all(&dir).unwrap();
    let text = |correctors: usize, viscosity: f64, steps: usize| {
        format!(
            r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [12, 12, 1]
[fluid]
viscosity = {viscosity}
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
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
wall_csv = false
[solver]
log_every = 1
correctors = {correctors}
[time]
time_step = 0.01
steps = {steps}
"#
        )
    };
    let file = dir.join("case.toml");
    std::fs::write(&file, text(2, 0.01, 20)).unwrap();
    std::fs::write(dir.join(".result_3"), "").unwrap();
    std::fs::write(dir.join(".read_toml_4"), "").unwrap();
    let case = FlowCase::from_file(&file).unwrap();
    let rep = run(&case, "reread", &dir, |l: &str| {
        if l.starts_with(" ") && l.split_whitespace().next() == Some("2") && l.contains("E-") {
            std::fs::write(&file, text(3, 0.02, 8)).unwrap();
        }
    })
    .unwrap();
    assert!(!rep.stopped);
    let log = std::fs::read_to_string(dir.join("flow.log")).unwrap();
    assert!(log.contains("changed solver.correctors: 2 -> 3"), "{log}");
    assert!(log.contains("changed time.steps: 20 -> 8"), "{log}");
    assert!(log.contains("not applied"), "{log}");
    assert!(log.contains("fluid.viscosity"), "{log}");
    assert!(dir.join("fields_000003.vtkhdf").exists());
    assert!(dir.join("case.read_toml_4.toml").exists());
    let r = restart::read(&dir.join("restart.h5")).unwrap();
    assert_eq!(r.transient.unwrap().1.step, 8);
    let _ = std::fs::remove_dir_all(&dir);
}

/// A channel with a parabolic inflow between slip walls (the profile
/// flattens downstream) and a convective outlet, error-controlled steps.
fn convective_channel(time: &str, extra: &str) -> FlowCase {
    let text = format!(
        r#"
[mesh.generate]
type = "mixed"
dimension = 2
upper = [2.0, 1.0, 0.1]
cells = [24, 12, 1]
perturb = 0.2
[fluid]
viscosity = 0.01
[flow]
convection = true
[boundary.xmin]
type = "velocity-inlet"
profile = "parabolic"
max_velocity = 1.5
profile_axes = ["y"]
[boundary.xmax]
type = "convective-outlet"
[boundary.ymin]
type = "slip-wall"
[boundary.ymax]
type = "slip-wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[output]
vtkhdf = false
wall_csv = false
[time]
{time}
{extra}
"#
    );
    FlowCase::from_toml(&text, Path::new(".")).unwrap()
}

/// The convective outlet keeps its old boundary levels: a run with
/// rejected steps continued from a restart file ends bit-identically to the
/// run without the interruption (D42), and the outlet carries the inflow.
#[test]
fn convective_outlet_restarts_bit_identically() {
    let time = "scheme = \"bdf2\"\ntime_step = 0.1\nsteps = 30\nadaptive = true\nmax_cfl = 4.0\n\
                error_tolerance = 1e-3";
    let full = convective_channel(time, "[run]\nrestart_at = [12]");
    let da = out("restart_conv_full");
    let ra = run_in(&full, "full", &da);
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    assert!(
        ra.solves[0].mass_imbalance < 1e-12,
        "{}",
        ra.solves[0].mass_imbalance
    );
    let log = std::fs::read_to_string(da.join("flow.log")).unwrap();
    assert!(
        log.contains(" rejected: error "),
        "no step was rejected:\n{log}"
    );
    let from = da.join("restart_000012.h5");
    let r = restart::read(&from).unwrap();
    assert!(r.transient.as_ref().unwrap().0.ub_old.is_some());
    let cont = convective_channel(time, &format!("[run]\nrestart = '{}'", from.display()));
    let db = out("restart_conv_cont");
    let rb = run_in(&cont, "cont", &db);
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    assert!(same_restart(&a, &b), "the continued run differs");
    let (oa, ob) = (
        &a.transient.as_ref().unwrap().0,
        &b.transient.as_ref().unwrap().0,
    );
    assert!(oa.ub_old.is_some() && oa.ub_old == ob.ub_old);
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// A transient run may start from the restart file of a steady run: its
/// fields are the start fields at `start_time`, with no older levels. A
/// converged steady cavity stays where it is over a few time steps.
#[test]
fn transient_run_starts_from_a_steady_restart() {
    let steady_text = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [16, 16, 1]
[fluid]
viscosity = 0.01
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
tolerance = 1e-10
max_iterations = 2000
[output]
vtkhdf = false
wall_csv = false
"#;
    let steady = FlowCase::from_toml(steady_text, Path::new(".")).unwrap();
    let ds = out("steady_start_steady");
    let rs = run_in(&steady, "steady", &ds);
    assert!(
        rs.solves.iter().all(|s| s.converged),
        "the steady cavity converges"
    );
    let from = ds.join("restart.h5");
    let tr = cavity(
        16,
        "scheme = \"bdf2\"\ntime_step = 0.05\nsteps = 6",
        &format!("[run]\nrestart = '{}'", from.display()),
    );
    let dt = out("steady_start_transient");
    let rt = run_in(&tr, "transient", &dt);
    assert!(rt.failures.is_empty(), "{:?}", rt.failures);
    let a = restart::read(&from).unwrap();
    let b = restart::read(&dt.join("restart.h5")).unwrap();
    assert!(a.transient.is_none(), "the steady file has no time levels");
    assert_eq!(
        b.transient.as_ref().unwrap().1.step,
        6,
        "six steps from t = 0"
    );
    let du = a
        .state
        .u
        .x
        .iter()
        .zip(&b.state.u.x)
        .chain(a.state.u.y.iter().zip(&b.state.u.y))
        .map(|(p, q)| (p - q).abs())
        .fold(0.0, f64::max);
    assert!(du < 1e-5, "the steady state moved by {du:e}");
    let log = std::fs::read_to_string(dt.join("flow.log")).unwrap();
    assert!(
        log.contains("start from the steady solution"),
        "the log says so"
    );
    let _ = std::fs::remove_dir_all(&ds);
    let _ = std::fs::remove_dir_all(&dt);
}
