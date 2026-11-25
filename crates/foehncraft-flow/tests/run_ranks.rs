//! MPI step W7 (D93, D84): a steady run on 2 to 4 ranks (threads of this
//! process) writes the output files of the run on one rank byte for byte:
//! the residuals, the ledger, the probes, the log columns of the run
//! record, the fields and the restart file; the result summary differs
//! only in the rank count of its provenance.

use foehncraft_flow::engine::{self, RunOptions};
use std::path::Path;

const CASE: &str = r#"
case_format = "0.7"
[mesh.generate]
type = "poly"
dimension = 2
upper = [3.0, 1.0, 0.1]
cells = [30, 10, 1]
perturb = 0.3
[mesh]
renumber = "rcm"
[fluid]
density = 1.0
viscosity = 0.01
[initial]
potential_flow = true
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 101325.0
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
velocity = [-1.0, 0.0, 0.0]
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
tolerance = 1e-30
max_iterations = 30
log_every = 5
[verification]
require_convergence = false
[solver.momentum]
limiter = "venkatakrishnan"
[balance]
every = 10
[[output.probes]]
name = "eddy"
point = [2.5, 0.8, 0.05]
[log]
columns = ["residual.momentum", "probe.eddy.velocity_x", "max.velocity_magnitude", "boundary.xmax.mass_flow"]
[run]
restart_every = 20
"#;

fn run(dir: &Path, ranks: usize) -> serde_json::Value {
    let mut case = foehncraft_flow::case::FlowCase::from_toml(CASE, &std::env::temp_dir()).unwrap();
    case.run.rank_count = ranks;
    case.run.threads = Some(ranks);
    let out = engine::run(&case, "w7", dir, &RunOptions::default(), |_| {}).unwrap();
    assert_eq!(
        out.exit_code(),
        0,
        "{ranks} ranks: {:?}",
        out.summary.diagnostics
    );
    let mut v = serde_json::to_value(&out.summary).unwrap();
    assert_eq!(v["provenance"]["ranks"], ranks);
    v["provenance"] = serde_json::Value::Null;
    // the hashes of the files without timings or settings
    if let Some(files) = v["files"].as_array_mut() {
        files.retain(|f| {
            !matches!(
                f["name"].as_str().unwrap_or(""),
                "flow.log"
                    | "flow.md"
                    | "flow.json"
                    | "run.jsonl"
                    | "case.resolved.toml"
                    | "ranks.vtkhdf"
            )
        });
    }
    v
}

/// The rows of the run record that do not carry timings.
fn record(dir: &Path) -> Vec<String> {
    std::fs::read_to_string(dir.join("run.jsonl"))
        .unwrap()
        .lines()
        .filter_map(|l| {
            let v: serde_json::Value = serde_json::from_str(l).ok()?;
            let kind = v["type"].as_str()?.to_string();
            matches!(kind.as_str(), "iteration" | "balance" | "table_row").then(|| l.to_string())
        })
        .collect()
}

#[test]
fn steady_run_on_ranks_writes_the_files_of_one_rank() {
    let base = std::env::temp_dir().join(format!("foehn-w7-run-{}", std::process::id()));
    let one = base.join("1");
    let reference = run(&one, 1);
    let files = [
        "residuals.csv",
        "balance.csv",
        "probes.csv",
        "boundaries.csv",
        #[cfg(feature = "hdf5")]
        "fields.vtkhdf",
        #[cfg(feature = "hdf5")]
        "restart.h5",
    ];
    for n in 2..=4 {
        let dir = base.join(n.to_string());
        let summary = run(&dir, n);
        for f in files {
            let (a, b) = (one.join(f), dir.join(f));
            if !a.exists() {
                assert!(!b.exists(), "{f} only on {n} ranks");
                continue;
            }
            assert!(
                std::fs::read(&a).unwrap() == std::fs::read(&b).unwrap(),
                "{f} on {n} ranks differs from one rank"
            );
        }
        assert_eq!(record(&one), record(&dir), "run.jsonl rows on {n} ranks");
        assert_eq!(
            summary["solves"], reference["solves"],
            "result.json solves on {n} ranks"
        );
        assert_eq!(
            summary["files"].as_array().map(Vec::len),
            reference["files"].as_array().map(Vec::len)
        );
    }
    let _ = std::fs::remove_dir_all(&base);
}

#[test]
fn refuses_what_does_not_run_on_ranks_yet() {
    let mut case = foehncraft_flow::case::FlowCase::from_toml(
        &format!("{CASE}\n[time]\nscheme = \"bdf2\"\ntime_step = 0.05\nsteps = 2\n"),
        &std::env::temp_dir(),
    )
    .unwrap();
    case.run.rank_count = 2;
    let err = foehncraft_flow::ranks::supported_run(&case).unwrap_err();
    assert!(err.to_string().contains("[time]"), "{err}");
}
