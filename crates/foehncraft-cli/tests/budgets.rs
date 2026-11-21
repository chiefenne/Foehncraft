//! Iteration and study budgets of `foehn flow` (D198 (c), (d)): the budget
//! in `--check`, the diagnosis of a steady run in the run record, the stop
//! by the diagnosis, the study ledger, and that a run without the new keys
//! computes the same fields.

use serde_json::Value;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const FOEHN: &str = env!("CARGO_BIN_EXE_foehn");

fn dir(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("budgets")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A steady cavity on 16 x 16 cells.
fn cavity(d: &Path) -> PathBuf {
    let p = d.join("cavity.toml");
    std::fs::write(
        &p,
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.0625]
cells = [16, 16, 1]

[fluid]
density = 1.0
viscosity = 0.01

[boundary.ymax]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
[boundary.ymin]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"

[solver]
tolerance = 1e-8
max_iterations = 2000
[solver.momentum]
relaxation = 0.9
"#,
    )
    .unwrap();
    p
}

fn foehn(args: &[&str]) -> Output {
    Command::new(FOEHN).args(args).output().expect("foehn runs")
}

fn json(o: &Output) -> Value {
    serde_json::from_slice(&o.stdout).unwrap_or_else(|e| {
        panic!(
            "stdout is not one JSON object ({e}):\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        )
    })
}

fn rows(out: &Path) -> Vec<Value> {
    std::fs::read_to_string(out.join("run.jsonl"))
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

fn run(case: &Path, out: &Path, sets: &[&str]) -> Output {
    let mut a = vec!["--json", "--quiet", "flow", s(case), "--out", s(out)];
    for x in sets {
        a.push("--set");
        a.push(x);
    }
    foehn(&a)
}

#[test]
fn check_shows_the_iteration_budget() {
    let d = dir("check");
    let case = cavity(&d);
    let o = foehn(&[
        "--json",
        "flow",
        "--check",
        s(&case),
        "--set",
        "run.limits.max_iterations=500",
    ]);
    let v = json(&o);
    let b = &v["estimate"]["iteration_budget"];
    assert_eq!(b["iterations"], 500);
    assert_eq!(b["source"], "run.limits.max_iterations");
    // the default diagnosis lies after this budget
    assert_eq!(b["diagnose_at"], Value::Null);
    assert_eq!(b["diagnosis"], "record");
    // a budget above 10 000 draws RUN021
    let o = foehn(&[
        "--json",
        "flow",
        "--check",
        s(&case),
        "--set",
        "solver.max_iterations=20000",
    ]);
    let v = json(&o);
    let codes: Vec<&str> = v["diagnostics"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| d["code"].as_str().unwrap())
        .collect();
    assert!(codes.contains(&"RUN021"), "{codes:?}");
    assert_eq!(v["ok"], true);
}

#[test]
fn diagnosis_is_recorded_and_changes_no_result() {
    let d = dir("record");
    let case = cavity(&d);
    let (a, b) = (d.join("a"), d.join("b"));
    let oa = run(&case, &a, &["run.limits.diagnose_at=30"]);
    assert_eq!(
        oa.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&oa.stderr)
    );
    let ob = run(&case, &b, &["run.limits.diagnose_at=0"]);
    assert_eq!(ob.status.code(), Some(0));
    let ra = rows(&a);
    let budget = ra
        .iter()
        .find(|r| r["type"] == "budget")
        .expect("budget row");
    assert_eq!(budget["iterations"]["iterations"], 2000);
    assert_eq!(budget["iterations"]["diagnose_at"], 30);
    let dg = ra
        .iter()
        .find(|r| r["type"] == "diagnosis")
        .expect("diagnosis row");
    assert_eq!(dg["iteration"], 30);
    assert_eq!(dg["verdict"], "on-track");
    assert!(
        dg["equations"]
            .as_array()
            .unwrap()
            .iter()
            .all(|e| e["equation"] != "momentum_z")
    );
    assert!(
        ra.iter()
            .any(|r| r["type"] == "event" && r["code"] == "RUN016")
    );
    assert!(!rows(&b).iter().any(|r| r["type"] == "diagnosis"));
    // the diagnosis only reads: the same iterations and fields
    for f in ["residuals.csv", "fields.vtkhdf"] {
        let (x, y) = (a.join(f), b.join(f));
        if x.exists() || y.exists() {
            assert!(
                std::fs::read(&x).unwrap() == std::fs::read(&y).unwrap(),
                "{f} differs"
            );
        }
    }
}

#[test]
fn diagnosis_stops_a_run_off_track() {
    let d = dir("stop");
    let case = cavity(&d);
    let out = d.join("out");
    // at iteration 30 the residual falls, but reaches 1e-8 only long after
    // a budget of 35 iterations
    let o = run(
        &case,
        &out,
        &[
            "solver.max_iterations=35",
            "run.limits.diagnose_at=30",
            "run.limits.diagnosis=\"stop\"",
        ],
    );
    assert_eq!(
        o.status.code(),
        Some(3),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let v = json(&o);
    assert_eq!(v["stop"]["code"], "RUN017");
    assert_eq!(v["stop"]["kind"], "diagnosis");
    let dg = rows(&out)
        .into_iter()
        .find(|r| r["type"] == "diagnosis")
        .unwrap();
    assert_eq!(dg["verdict"], "slow");
}

#[test]
fn study_budget_stops_and_refuses() {
    let d = dir("study");
    let case = cavity(&d);
    let study = d.join("study");
    std::fs::create_dir_all(&study).unwrap();
    let led = study.join("study_budget.jsonl");
    // a budget smaller than any run: the first run starts, stops at it
    let o = run(
        &case,
        &study.join("v0000"),
        &["run.limits.study_cell_hours=1e-9"],
    );
    assert_eq!(
        o.status.code(),
        Some(3),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert_eq!(json(&o)["stop"]["code"], "RUN018");
    let r = rows(&study.join("v0000"));
    let b = r.iter().find(|r| r["type"] == "budget").unwrap();
    assert_eq!(b["study"]["cell_hours"]["budget"], 1e-9);
    assert_eq!(b["study"]["cell_hours"]["spent"], 0.0);
    assert_eq!(b["estimate"]["cells"], 256);
    let spend = r.iter().find(|r| r["type"] == "spend").unwrap();
    assert!(spend["cell_hours"].as_f64().unwrap() > 0.0);
    assert_eq!(std::fs::read_to_string(&led).unwrap().lines().count(), 1);
    // the next run of the study does not start
    let o = run(
        &case,
        &study.join("v0001"),
        &["run.limits.study_cell_hours=1e-9"],
    );
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stdout).contains("RUN019"));
    let o = foehn(&[
        "--json",
        "flow",
        "--check",
        s(&case),
        "--out",
        s(&study.join("v0002")),
        "--set",
        "run.limits.study_cell_hours=1e-9",
    ]);
    let v = json(&o);
    assert_eq!(v["ok"], false);
    assert!(
        v["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|d| d["code"] == "RUN019")
    );
    // a generous budget: the run finishes and adds its line
    let o = run(
        &case,
        &study.join("v0003"),
        &[
            "run.limits.study_cell_hours=100.0",
            "run.limits.study_output=1e9",
        ],
    );
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stdout)
    );
    assert_eq!(std::fs::read_to_string(&led).unwrap().lines().count(), 2);
}
