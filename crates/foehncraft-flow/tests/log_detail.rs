//! The log additions of D214: `[log] detail`, `log.columns` during a run,
//! the Limitations block before the run and the Events block before the
//! closing frame (docs/requirements/log-redesign-review.md).

use foehncraft_flow::case::FlowCase;
use foehncraft_flow::engine::{RunOptions, run};
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("foehncraft_logd_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

const CAVITY: &str = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [10, 10, 1]
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
max_iterations = 300
log_every = 10
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
steps = 6
"#;

/// Runs `text` through the driver; the log and the run record.
fn go(text: &str, name: &str) -> (String, String) {
    let d = out(name);
    let c = FlowCase::from_toml(text, Path::new(".")).unwrap();
    run(&c, "cavity", &d, &RunOptions::default(), |_| {}).unwrap();
    let log = std::fs::read_to_string(d.join("flow.log")).unwrap();
    let rec = std::fs::read_to_string(d.join("run.jsonl")).unwrap();
    let _ = std::fs::remove_dir_all(&d);
    (log, rec)
}

fn rows(rec: &str, kind: &str) -> usize {
    rec.lines()
        .filter(|l| l.contains(&format!("\"type\":\"{kind}\"")))
        .count()
}

#[test]
fn summary_detail_leaves_out_the_table_but_not_the_record() {
    let (full, rec_full) = go(CAVITY, "full");
    let (sum, rec_sum) = go(&format!("{CAVITY}[log]\ndetail = \"summary\"\n"), "summary");
    assert!(full.contains("ResU"), "{full}");
    assert!(!sum.contains("ResU"), "{sum}");
    assert!(sum.contains("[log] detail = \"summary\""), "{sum}");
    // the header blocks, the summaries and the closing frame stay
    for s in ["Boundary conditions", "Convergence", "Exit code"] {
        assert!(sum.contains(s), "{s}: {sum}");
    }
    // the run record keeps every row
    assert!(rows(&rec_full, "table_row") > 0);
    assert_eq!(rows(&rec_full, "table_row"), rows(&rec_sum, "table_row"));
}

#[test]
fn steps_detail_leaves_out_the_outer_iterations() {
    let (full, _) = go(TRANSIENT, "tfull");
    let (steps, _) = go(&format!("{TRANSIENT}[log]\ndetail = \"steps\"\n"), "tsteps");
    let (sum, _) = go(&format!("{TRANSIENT}[log]\ndetail = \"summary\"\n"), "tsum");
    // the header of the table of the outer iterations of a step
    let outer = "Outer      ResU";
    assert!(full.contains(outer), "{full}");
    assert!(full.contains("TIME STEP"));
    assert!(!steps.contains(outer), "{steps}");
    assert!(steps.contains("TIME STEP"), "{steps}");
    assert!(steps.contains("Step "), "{steps}");
    assert!(!sum.contains("TIME STEP"), "{sum}");
    assert!(sum.contains("Exit code"), "{sum}");
}

#[test]
fn limitations_are_listed_before_the_run() {
    let (log, _) = go(CAVITY, "lim");
    let at = log
        .find("Limitations known before the run")
        .expect("the Limitations block");
    let bc = log.find("Boundary conditions").unwrap();
    let iters = log.find("Iterations (").unwrap();
    assert!(bc < at && at < iters, "{log}");
    assert!(log[at..iters].contains("LIM004"), "{log}");
}

#[test]
fn scheduled_columns_reprint_the_header_and_land_in_the_events_block() {
    let (plain, _) = go(CAVITY, "noevents");
    assert!(!plain.contains("Events:"), "{plain}");
    let (log, rec) = go(
        &format!(
            "{CAVITY}[[schedule]]\nat_iteration = 21\nset = {{ log.columns = [\"residual.continuity\", \"max.velocity_magnitude\"] }}\n"
        ),
        "cols",
    );
    // two `columns` rows: the start and the change
    assert_eq!(rows(&rec, "columns"), 2, "{rec}");
    // the new header after the change
    let change = log.find("schedule (iteration 21)").expect("the change");
    assert!(log[change..].contains("ResM"), "{log}");
    assert!(!log[change..].contains("ResU"), "{log}");
    // the Events block before the closing frame
    let ev = log.find("Events:").expect("the Events block");
    let close = log.rfind("Exit code").unwrap();
    assert!(ev < close);
    assert!(log[ev..close].contains("RUN015"), "{log}");
    assert!(log[ev..close].contains("log.columns"), "{log}");
}

#[test]
fn a_column_change_through_the_settings_validates() {
    // `log.columns` is a run-time key; a wrong name is refused as before
    let c = FlowCase::from_toml(CAVITY, Path::new(".")).unwrap();
    let ok = foehncraft_flow::control::with_setting(
        &c,
        "log.columns",
        &toml::Value::Array(vec![toml::Value::String("residual.momentum".into())]),
    );
    assert!(ok.is_ok(), "{ok:?}");
    let bad = foehncraft_flow::control::with_setting(
        &c,
        "log.columns",
        &toml::Value::Array(vec![toml::Value::String("residual.nothing".into())]),
    );
    assert!(bad.is_err());
    assert!(foehncraft_flow::live::settable("log.columns"));
    assert!(foehncraft_flow::schedule::schedulable("log.columns"));
}
