//! `[[schedule]]` (D128): settings that change at prepared iterations,
//! steps and times; docs/recipes/runs.md.

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::control;
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("foehncraft_sched_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

const CAVITY: &str = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [12, 12, 1]
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
tolerance = 1e-5
max_iterations = 400
[output]
vtkhdf = false
wall_csv = false
"#;

const TRANSIENT: &str = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [8, 8, 1]
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
scheme = "bdf2"
time_step = 0.01
end_time = 0.1
"#;

fn parse(text: &str) -> anyhow::Result<FlowCase> {
    FlowCase::from_toml(text, Path::new("."))
}

fn message(text: &str) -> String {
    format!("{:#}", parse(text).unwrap_err())
}

#[test]
fn a_valid_schedule_is_read_and_listed() {
    let c = parse(&format!(
        "{CAVITY}\n[[schedule]]\nat_iteration = 50\nset = {{ solver.momentum.relaxation = 0.5 }}\n\
         [[schedule]]\nat_iteration = 100\nset = {{ \"solver.pressure.relaxation\" = 0.3, \"solver.log_every\" = 10 }}\n"
    ))
    .unwrap();
    assert_eq!(c.schedule.len(), 2);
    assert_eq!(c.schedule[1].set.len(), 2);
    // the resolved case keeps it (the run is repeatable from its case file)
    let t = toml::to_string_pretty(&c).unwrap();
    assert!(t.contains("[[schedule]]"), "{t}");
    let again = parse(&t).unwrap();
    assert_eq!(again.schedule, c.schedule);
    // no schedule: nothing in the resolved case
    assert!(
        !toml::to_string_pretty(&parse(CAVITY).unwrap())
            .unwrap()
            .contains("schedule")
    );
}

#[test]
fn invalid_entries_are_refused_when_the_case_is_read() {
    let sched = |e: &str| format!("{CAVITY}\n[[schedule]]\n{e}\n");
    // no trigger, two triggers, the wrong kind for the case
    assert!(message(&sched("set = { solver.tolerance = 1e-6 }")).contains("exactly one"));
    assert!(
        message(&sched(
            "at_iteration = 5\nat_step = 5\nset = { solver.tolerance = 1e-6 }"
        ))
        .contains("exactly one")
    );
    assert!(
        message(&sched("at_step = 5\nset = { solver.tolerance = 1e-6 }"))
            .contains("does not fit a steady case")
    );
    assert!(
        message(&sched(
            "at_iteration = 0\nset = { solver.tolerance = 1e-6 }"
        ))
        .contains("at least 1")
    );
    // a key that cannot change during a run, with the reason
    let m = message(&sched("at_iteration = 5\nset = { fluid.viscosity = 0.1 }"));
    assert!(m.contains("cannot change during a run"), "{m}");
    let m = message(&sched("at_iteration = 5\nset = { run.reread_every = 10 }"));
    assert!(m.contains("cannot be scheduled"), "{m}");
    // a misspelled key gets the nearest one
    let m = message(&sched(
        "at_iteration = 5\nset = { solver.momentum.relaxatoin = 0.5 }",
    ));
    assert!(
        m.contains("did you mean `solver.momentum.relaxation`"),
        "{m}"
    );
    // a transient key in a steady case
    let m = message(&sched("at_iteration = 5\nset = { time.time_step = 0.1 }"));
    assert!(m.contains("belongs to a transient case"), "{m}");
    // a value that does not validate names the entry
    let m = message(&sched(
        "at_iteration = 5\nset = { solver.momentum.relaxation = 1.5 }",
    ));
    assert!(m.contains("entry 1"), "{m}");
    // out of order
    let m = message(&format!(
        "{CAVITY}\n[[schedule]]\nat_iteration = 50\nset = {{ solver.tolerance = 1e-6 }}\n\
         [[schedule]]\nat_iteration = 20\nset = {{ solver.tolerance = 1e-7 }}\n"
    ));
    assert!(m.contains("ascending order"), "{m}");
    // the code of the diagnostic
    let e = parse(&sched("at_iteration = 5\nset = { fluid.viscosity = 0.1 }")).unwrap_err();
    assert_eq!(diag::of_error(&e).code, "SCHED001");
}

#[test]
fn transient_triggers_are_checked() {
    let sched = |e: &str| format!("{TRANSIENT}\n[[schedule]]\n{e}\n");
    assert!(
        message(&sched(
            "at_iteration = 5\nset = { solver.tolerance = 1e-6 }"
        ))
        .contains("does not fit a transient case")
    );
    assert!(
        message(&sched("at_time = 0.0\nset = { time.time_step = 0.005 }"))
            .contains("later than the start time")
    );
    let m = message(&format!(
        "{TRANSIENT}\n[[schedule]]\nat_step = 3\nset = {{ time.time_step = 0.005 }}\n\
         [[schedule]]\nat_time = 0.05\nset = {{ time.time_step = 0.01 }}\n"
    ));
    assert!(m.contains("one kind of trigger"), "{m}");
    parse(&sched("at_step = 3\nset = { time.time_step = 0.005 }")).unwrap();
    parse(&sched("at_time = 0.03\nset = { time.time_step = 0.005 }")).unwrap();
}

/// A schedule that sets what is already in force changes nothing: the run
/// is bitwise the one without a schedule; one that changes a setting logs
/// the change and the run differs.
#[test]
fn steady_schedule_logs_changes_and_noop_changes_nothing() {
    let base = |extra: &str, name: &str| {
        let d = out(name);
        let c = parse(&format!("{CAVITY}{extra}")).unwrap();
        let r = run(&c, "cavity", &d, |_| {}).unwrap();
        let log = std::fs::read_to_string(d.join("flow.log")).unwrap();
        let hist = std::fs::read_to_string(d.join("residuals.csv")).unwrap();
        let _ = std::fs::remove_dir_all(&d);
        (r, log, hist)
    };
    let (_, log0, hist0) = base("", "none");
    assert!(!log0.contains("schedule"));
    // sets what is in force already: identical history, no change logged
    let (_, log1, hist1) = base(
        "[[schedule]]\nat_iteration = 10\nset = { solver.tolerance = 1e-5 }\n",
        "noop",
    );
    assert_eq!(hist0, hist1);
    assert!(!log1.contains("schedule (iteration 10)"), "{log1}");
    // a real change
    let (_, log2, hist2) = base(
        "[[schedule]]\nat_iteration = 10\nset = { solver.momentum.relaxation = 0.4, solver.log_every = 7 }\n",
        "change",
    );
    assert!(
        log2.contains("iteration 10: schedule (iteration 10): solver.momentum.relaxation"),
        "{log2}"
    );
    assert!(log2.contains("solver.log_every"), "{log2}");
    assert!(log2.contains("Scheduled: iteration 10"), "{log2}");
    // before iteration 10 the history is the one of the plain run
    let head = |h: &str| {
        h.lines()
            .filter(|l| !l.starts_with('#'))
            .take(10)
            .map(String::from)
            .collect::<Vec<_>>()
    };
    assert_eq!(head(&hist0), head(&hist2));
    assert_ne!(hist0, hist2);
}

#[test]
fn transient_schedule_lands_on_the_time_and_changes_the_step() {
    let d = out("transient");
    let c = parse(&format!(
        "{TRANSIENT}\n[[schedule]]\nat_time = 0.0355\nset = {{ time.time_step = 0.005 }}\n"
    ))
    .unwrap();
    run(&c, "cavity", &d, |_| {}).unwrap();
    let log = std::fs::read_to_string(d.join("flow.log")).unwrap();
    assert!(log.contains("time.time_step"), "{log}");
    assert!(log.contains("schedule (t = 0.0355 s)"), "{log}");
    let _ = std::fs::remove_dir_all(&d);
}

/// A re-read keeps the entries that have fired in force over the values
/// of the file, and reports a changed schedule as not applied.
#[test]
fn reread_keeps_fired_entries() {
    let d = out("reread");
    std::fs::create_dir_all(&d).unwrap();
    let f = d.join("case.toml");
    let text = format!(
        "{CAVITY}[[schedule]]\nat_iteration = 10\nset = {{ solver.momentum.relaxation = 0.4 }}\n"
    );
    std::fs::write(&f, &text).unwrap();
    let c = FlowCase::from_file(&f).unwrap();
    let mut r = foehncraft_flow::schedule::Runner::default();
    let (now, ch) = r
        .fire(
            &c,
            foehncraft_flow::schedule::Position {
                n: 10,
                time: 0.0,
                tol: 0.0,
            },
        )
        .unwrap()
        .unwrap();
    assert_eq!(ch.len(), 1);
    assert_eq!(now.solver.momentum.relaxation, Some(0.4));
    // nothing is due before iteration 10
    let mut r0 = foehncraft_flow::schedule::Runner::default();
    assert!(
        r0.fire(
            &c,
            foehncraft_flow::schedule::Position {
                n: 9,
                time: 0.0,
                tol: 0.0,
            },
        )
        .unwrap()
        .is_none()
    );
    // re-reading the unchanged file does not undo the entry
    let rr = control::reread(&f, &now, r.fired).unwrap();
    assert!(rr.applied.is_empty(), "{:?}", rr.applied);
    assert_eq!(rr.case.solver.momentum.relaxation, Some(0.4));
    // editing the schedule in the file is not applied
    std::fs::write(&f, text.replace("0.4", "0.2")).unwrap();
    let rr = control::reread(&f, &now, r.fired).unwrap();
    assert!(rr.not_applied.iter().any(|l| l.starts_with("schedule")));
    assert_eq!(rr.case.solver.momentum.relaxation, Some(0.4));
    let _ = std::fs::remove_dir_all(&d);
}

/// The scalar and energy sub-tables of the boundaries are kept by the
/// reader outside the TOML of the case; a schedule entry, which re-reads
/// the case with its settings, keeps them (an inlet of an energy case
/// needs its temperature).
#[test]
fn a_schedule_keeps_the_energy_conditions_of_the_boundaries() {
    let text = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [2.0, 1.0, 0.05]
cells = [8, 4, 1]
[fluid]
viscosity = 0.01
specific_heat = 1000.0
conductivity = 0.1
[energy]
initial = 300.0
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmin.energy]
temperature = 300.0
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "wall"
[boundary.ymin.energy]
heat_flux = 100.0
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[[schedule]]
at_iteration = 10
set = { solver.momentum.relaxation = 0.8 }
"#;
    let c = parse(text).unwrap();
    assert_eq!(c.boundary_energy.len(), 2);
    let e = foehncraft_flow::schedule::apply_entry(&c, &c.schedule[0], None).unwrap();
    assert_eq!(e.boundary_energy.len(), 2);
    assert_eq!(e.solver.momentum.relaxation, Some(0.8));
}

/// A first-order start: the momentum scheme switches from `upwind` to
/// `linear-upwind` at a prepared iteration. Before it the history is the
/// one of an upwind run, the change is logged, and the run converges to the
/// solution of the plain linear-upwind run. A switch to or from the bounded
/// central scheme is refused when the case is read.
#[test]
fn the_momentum_scheme_switches_from_upwind_to_linear_upwind() {
    let go = |extra: &str, name: &str| {
        let d = out(name);
        let text = CAVITY
            .replace("max_iterations = 400", "max_iterations = 1500")
            .replace("tolerance = 1e-5", "tolerance = 1e-9");
        let line = "lines = [{ name = \"mid\", from = [0.5, 0.0, 0.025], to = [0.5, 1.0, 0.025], points = 21 }]\n";
        let c = parse(&format!("{text}{line}{extra}")).unwrap();
        let r = run(&c, "cavity", &d, |_| {}).unwrap();
        let log = std::fs::read_to_string(d.join("flow.log")).unwrap();
        let hist = std::fs::read_to_string(d.join("residuals.csv")).unwrap();
        let u: Vec<f64> = std::fs::read_to_string(d.join("line_mid.csv"))
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#') && !l.starts_with('s'))
            .map(|l| l.split(',').nth(4).unwrap().parse().unwrap())
            .collect();
        let _ = std::fs::remove_dir_all(&d);
        (r, log, hist, u)
    };
    let (lu, _, _, u_lu) = go("", "lu");
    let (up, _, hist_up, u_up) = go("[solver.momentum]\nscheme = \"upwind\"\n", "up");
    let (sw, log_sw, hist_sw, u_sw) = go(
        "[solver.momentum]\nscheme = \"upwind\"\n[[schedule]]\nat_iteration = 20\nset = { solver.momentum.scheme = \"linear-upwind\" }\n",
        "switch",
    );
    assert!(lu.solves[0].converged && up.solves[0].converged && sw.solves[0].converged);
    assert!(
        log_sw.contains("schedule (iteration 20): solver.momentum.scheme"),
        "{log_sw}"
    );
    let head = |h: &str| {
        h.lines()
            .filter(|l| !l.starts_with('#'))
            .take(20)
            .map(String::from)
            .collect::<Vec<_>>()
    };
    assert_eq!(head(&hist_up), head(&hist_sw), "upwind before the switch");
    // the converged state is the one of linear upwind, not of upwind
    let diff = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max)
    };
    assert!(
        diff(&u_lu, &u_sw) < 1e-6,
        "switched {:e}",
        diff(&u_lu, &u_sw)
    );
    assert!(
        diff(&u_lu, &u_up) > 1e-3,
        "upwind differs {:e}",
        diff(&u_lu, &u_up)
    );
    let m = message(&format!(
        "{CAVITY}\n[solver.momentum]\nscheme = \"central\"\n[[schedule]]\nat_iteration = 5\nset = {{ solver.momentum.scheme = \"upwind\" }}\n"
    ));
    assert!(m.contains("only between"), "{m}");
}
