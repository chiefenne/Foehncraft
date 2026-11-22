//! The run contract of `foehn flow` (D62, D64, D78, D81, D84): JSON on
//! standard output, `result.json`, exit codes, `--check --json`, limits,
//! `--restart latest` and, on unix, a clean stop on SIGTERM.

use serde_json::Value;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const FOEHN: &str = env!("CARGO_BIN_EXE_foehn");

/// A fresh directory for one test.
fn dir(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("run_contract")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

const MESH_AND_WALLS: &str = r#"
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
"#;

/// A steady cavity on 16 x 16 cells.
fn steady(d: &Path) -> PathBuf {
    let p = d.join("cavity.toml");
    let text = format!(
        "{MESH_AND_WALLS}\n[solver]\ntolerance = 1e-8\nmax_iterations = 2000\n[solver.momentum]\nrelaxation = 0.9\n"
    );
    std::fs::write(&p, text).unwrap();
    p
}

/// A transient cavity on 16 x 16 cells.
#[cfg(unix)]
fn transient(d: &Path, steps: usize) -> PathBuf {
    let p = d.join("cavity_transient.toml");
    let text = format!(
        "{MESH_AND_WALLS}\n[solver]\nlog_every = 1000\n\n[time]\nscheme = \"bdf2\"\ntime_step = 0.01\nsteps = {steps}\n"
    );
    std::fs::write(&p, text).unwrap();
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

fn bytes(p: &Path) -> Vec<u8> {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

#[test]
fn json_result_and_byte_identical_files() {
    let d = dir("identical");
    let case = steady(&d);
    let (a, b) = (d.join("a"), d.join("b"));
    for out in [&a, &b] {
        let o = foehn(&["--json", "--quiet", "flow", s(&case), "--out", s(out)]);
        assert_eq!(
            o.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&o.stderr)
        );
        let v = json(&o);
        assert_eq!(v["format"], "foehncraft.result");
        assert_eq!(v["status"], "finished");
        assert_eq!(v["exit_code"], 0);
        assert_eq!(v["solves"][0]["cells"], 256);
        // one mesh: a limitation (D199) that changes neither status nor exit code
        assert!(
            v["limitations"]
                .as_array()
                .unwrap()
                .iter()
                .any(|l| l["code"] == "LIM004"),
            "{}",
            v["limitations"]
        );
        let log = std::fs::read_to_string(out.join("flow.log")).unwrap();
        assert!(log.contains("Run finished with"), "{log}");
        let md = std::fs::read_to_string(out.join("flow.md")).unwrap();
        assert!(md.contains("## Limitations ("), "{md}");
        assert!(v["provenance"]["case_hash"].as_str().unwrap().len() == 64);
        // result.json is the object on standard output
        let file: Value = serde_json::from_slice(&bytes(&out.join("result.json"))).unwrap();
        assert_eq!(file, v);
        // the run record starts with the header and ends with the end row
        let rec = std::fs::read_to_string(out.join("run.jsonl")).unwrap();
        let rows: Vec<Value> = rec
            .lines()
            .map(|l| serde_json::from_str(l).unwrap())
            .collect();
        assert_eq!(rows[0]["type"], "header");
        assert_eq!(rows[0]["format"], "foehncraft.run_record");
        assert!(rows.iter().any(|r| r["type"] == "iteration"));
        assert_eq!(rows.last().unwrap()["type"], "end");
    }
    for f in ["result.json", "fields.vtkhdf", "residuals.csv"] {
        assert!(
            bytes(&a.join(f)) == bytes(&b.join(f)),
            "{f} differs between two runs"
        );
    }
    // result.json lists the hash of every output except the measured files
    let v: Value = serde_json::from_slice(&bytes(&a.join("result.json"))).unwrap();
    let names: Vec<&str> = v["files"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f["name"].as_str().unwrap())
        .collect();
    assert!(names.contains(&"fields.vtkhdf"), "{names:?}");
    assert!(
        !names.contains(&"flow.log") && !names.contains(&"run.jsonl"),
        "{names:?}"
    );
}

/// D199 evidence (schema): every row of `run.jsonl` is one JSON object with a
/// known `type` and the keys of that type; the header names the format.
#[test]
fn run_record_rows_follow_their_schema() {
    let d = dir("record_schema");
    let case = steady(&d);
    let out = d.join("out");
    let o = foehn(&["--json", "--quiet", "flow", s(&case), "--out", s(&out)]);
    assert_eq!(o.status.code(), Some(0));
    let rec = std::fs::read_to_string(out.join("run.jsonl")).unwrap();
    let required: &[(&str, &[&str])] = &[
        ("header", &["format", "format_version", "case", "provenance"]),
        ("columns", &["columns"]),
        ("iteration", &["iteration"]),
        ("timings", &["seconds"]),
        ("result", &["result"]),
        ("end", &["status", "exit_code"]),
    ];
    let optional = [
        "section", "table_row", "balance", "event", "diagnosis", "budget", "spend", "call",
    ];
    let mut seen = std::collections::BTreeSet::new();
    for (n, line) in rec.lines().enumerate() {
        let row: Value = serde_json::from_str(line).unwrap_or_else(|e| panic!("row {n}: {e}"));
        let t = row["type"].as_str().unwrap_or_else(|| panic!("row {n} has no type"));
        match required.iter().find(|(name, _)| *name == t) {
            Some((_, keys)) => {
                for k in *keys {
                    assert!(row.get(*k).is_some(), "row {n} ({t}) lacks `{k}`");
                }
            }
            None => assert!(optional.contains(&t), "row {n}: unknown type `{t}`"),
        }
        seen.insert(t.to_string());
    }
    for (name, _) in required.iter().filter(|(n, _)| *n != "columns" && *n != "timings") {
        assert!(seen.contains(*name), "no `{name}` row in {seen:?}");
    }
}

#[test]
fn exit_codes() {
    let d = dir("exit_codes");
    let case = steady(&d);
    // 1: a value out of range, as a JSON error object with code and key
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&case),
        "--set",
        "fluid.viscosity=-1",
        "--out",
        s(&d.join("e1")),
    ]);
    assert_eq!(o.status.code(), Some(1));
    let v = json(&o);
    assert_eq!(v["format"], "foehncraft.error");
    assert_eq!(v["diagnostics"][0]["code"], "CASE004");
    assert_eq!(v["diagnostics"][0]["key"], "fluid.viscosity");
    // 1: an unknown key, with its location and a hint
    let bad = d.join("typo.toml");
    std::fs::write(
        &bad,
        std::fs::read_to_string(&case)
            .unwrap()
            .replace("cells = [16", "cels = [16"),
    )
    .unwrap();
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&bad),
        "--out",
        s(&d.join("e1b")),
    ]);
    assert_eq!(o.status.code(), Some(1));
    let v = json(&o);
    assert_eq!(v["diagnostics"][0]["code"], "CASE002");
    assert_eq!(v["diagnostics"][0]["key"], "mesh.generate.cels");
    // 2: the run ends, but it did not converge
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&case),
        "--set",
        "solver.max_iterations=5",
        "--set",
        "verification.require_convergence=true",
        "--out",
        s(&d.join("e2")),
    ]);
    assert_eq!(
        o.status.code(),
        Some(2),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let v = json(&o);
    assert_eq!(v["status"], "failed-checks");
    assert!(
        v["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|d| d["code"] == "RUN006")
    );
    // 3: stopped at a limit
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&case),
        "--set",
        "run.limits.max_iterations=10",
        "--out",
        s(&d.join("e3")),
    ]);
    assert_eq!(o.status.code(), Some(3));
    let v = json(&o);
    assert_eq!(v["status"], "stopped");
    assert_eq!(v["stop"]["code"], "RUN004");
    // 3: the .stop control file of an earlier run blocks the start
    std::fs::write(d.join("e3").join(".stop"), "").unwrap();
    let o = foehn(&[
        "--json",
        "--quiet",
        "flow",
        s(&case),
        "--out",
        s(&d.join("e3")),
    ]);
    assert_eq!(o.status.code(), Some(1));
}

#[test]
fn check_json_with_estimate() {
    let d = dir("check");
    let case = steady(&d);
    let o = foehn(&["--json", "flow", s(&case), "--check"]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let v = json(&o);
    assert_eq!(v["format"], "foehncraft.check");
    assert_eq!(v["ok"], true);
    // the pressure level is held in a reference cell: a note with its code
    assert!(
        v["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|d| d["code"] == "BC015")
    );
    let e = &v["estimate"];
    assert_eq!(e["cells"], 256);
    assert!(e["host_memory"].as_f64().unwrap() > 1e6);
    assert!(e["output_bytes"].as_f64().unwrap() > 0.0);
    // a broken case is a report with ok = false, not a crash
    let o = foehn(&[
        "--json",
        "flow",
        s(&case),
        "--check",
        "--set",
        "boundary.xmin.type=\"nozzle\"",
    ]);
    assert_ne!(o.status.code(), Some(0));
    let v = json(&o);
    assert_eq!(v["ok"], false);
}

#[test]
fn schema_is_json_schema() {
    let o = foehn(&["schema"]);
    assert_eq!(o.status.code(), Some(0));
    let v: Value = serde_json::from_slice(&o.stdout).unwrap();
    assert!(v["$schema"].as_str().unwrap().contains("2020-12"));
    assert!(v["properties"]["fluid"].is_object());
}

#[test]
fn steady_restart_latest_is_bitwise() {
    let d = dir("steady_restart");
    let case = steady(&d);
    let (straight, split) = (d.join("straight"), d.join("split"));
    let o = foehn(&["--quiet", "flow", s(&case), "--out", s(&straight)]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let o = foehn(&[
        "--quiet",
        "flow",
        s(&case),
        "--set",
        "run.limits.max_iterations=40",
        "--out",
        s(&split),
    ]);
    assert_eq!(o.status.code(), Some(3));
    assert!(split.join("restart.h5").exists());
    let o = foehn(&[
        "--quiet",
        "flow",
        s(&case),
        "--restart",
        "latest",
        "--out",
        s(&split),
    ]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    for f in ["fields.vtkhdf", "restart.h5"] {
        assert!(
            bytes(&straight.join(f)) == bytes(&split.join(f)),
            "{f}: stop and restart differs from the straight run"
        );
    }
}

#[cfg(unix)]
#[test]
fn sigterm_stops_cleanly_and_restart_latest_continues_bitwise() {
    let d = dir("sigterm");
    let case = transient(&d, 600);
    let (straight, split) = (d.join("straight"), d.join("split"));
    let o = foehn(&["--quiet", "flow", s(&case), "--out", s(&straight)]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );

    let mut child = Command::new(FOEHN)
        .args(["--quiet", "flow", s(&case), "--out", s(&split)])
        .spawn()
        .unwrap();
    // wait for the first step rows of the run record
    let rec = split.join("run.jsonl");
    let t0 = std::time::Instant::now();
    loop {
        let n = std::fs::read_to_string(&rec).map_or(0, |t| t.matches("\"step\"").count());
        if n > 0 {
            break;
        }
        assert!(t0.elapsed().as_secs() < 120, "no step rows in run.jsonl");
        std::thread::sleep(std::time::Duration::from_millis(20));
    }
    let k = Command::new("kill")
        .args(["-TERM", &child.id().to_string()])
        .status()
        .unwrap();
    assert!(k.success());
    let st = child.wait().unwrap();
    assert_eq!(st.code(), Some(3));
    let v: Value = serde_json::from_slice(&bytes(&split.join("result.json"))).unwrap();
    assert_eq!(v["stop"]["code"], "RUN002");
    assert_eq!(v["stop"]["signal"], "SIGTERM");

    let o = foehn(&[
        "--quiet",
        "flow",
        s(&case),
        "--restart",
        "latest",
        "--out",
        s(&split),
    ]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    for f in ["fields.vtkhdf", "restart.h5"] {
        assert!(
            bytes(&straight.join(f)) == bytes(&split.join(f)),
            "{f}: stop and restart differs from the straight run"
        );
    }
}

/// `foehn post derive` (D72): the fields of a run with `[output] derived`
/// again from its restart file, byte for byte, and its errors.
#[test]
fn post_derive_from_the_restart_file() {
    let d = dir("post_derive");
    let case = steady(&d);
    let out = d.join("run");
    let derived = r#"output.derived=["vorticity","q-criterion","lambda2","wall-shear-stress"]"#;
    let o = foehn(&[
        "--quiet",
        "flow",
        s(&case),
        "--set",
        derived,
        "--set",
        "run.restart_final=true",
        "--out",
        s(&out),
    ]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let o = foehn(&[
        "--json",
        "--quiet",
        "post",
        "derive",
        s(&case),
        "--set",
        derived,
        "--out",
        s(&out),
    ]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let v = json(&o);
    assert_eq!(v["format"], "foehncraft.command");
    assert_eq!(v["command"], "post");
    assert_eq!(v["cells"], 256);
    let names: Vec<&str> = v["datasets"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap())
        .collect();
    assert!(names.contains(&"q_criterion") && names.contains(&"wall_shear_stress_magnitude"));
    assert!(
        bytes(&out.join("fields.vtkhdf")) == bytes(&out.join("derived.vtkhdf")),
        "post derive differs from the fields of the run"
    );
    // another list into another file
    let file = d.join("vortices.vtkhdf");
    let o = foehn(&[
        "--quiet",
        "post",
        "derive",
        s(&case),
        "--fields",
        "enstrophy,strain-rate",
        "--file",
        s(&file),
        "--out",
        s(&out),
    ]);
    assert_eq!(o.status.code(), Some(0));
    assert!(file.exists());
    // errors with their codes
    for (args, code) in [
        (vec!["--fields", "helicity"], "CASE003"),
        (vec!["--fields", "pressure-coefficient"], "CASE005"),
        (vec!["--out", s(&d.join("empty"))], "RUN011"),
    ] {
        let mut a = vec!["--json", "--quiet", "post", "derive", s(&case)];
        if !args.contains(&"--out") {
            a.extend(["--out", s(&out)]);
        }
        a.extend(args.iter().copied());
        let o = foehn(&a);
        assert_eq!(o.status.code(), Some(1), "{a:?}");
        assert_eq!(json(&o)["diagnostics"][0]["code"], code, "{a:?}");
    }
}

/// `[log] columns` (D125): selectable columns with extremes and cells, the
/// same content in `flow.log` and `run.jsonl`, and no change of any result.
#[test]
fn configurable_log_columns_change_no_result() {
    let d = dir("log_columns");
    let plain = steady(&d);
    let custom = d.join("columns").join("cavity.toml");
    let mut text = std::fs::read_to_string(&plain).unwrap();
    // the case name is part of the files, so both cases have one name; the
    // probe is in both, because the columns refer to it
    text.push_str("\n[[output.probes]]\nname = \"mid\"\npoint = [0.5, 0.5, 0.03]\n");
    std::fs::write(&plain, &text).unwrap();
    let d2 = d.join("columns");
    std::fs::create_dir_all(&d2).unwrap();
    text.push_str(
        "\n[log]\ndigits = 4\ncolumns = [\n  \"residual.momentum\",\n  \"residual.continuity\",\n  \"probe.mid.velocity_x\",\n  \"cell.100.pressure\",\n  \"max.velocity_magnitude\",\n  \"min.pressure\",\n  \"boundary.ymax.mass_flow\",\n  \"time.per_iteration\",\n]\n",
    );
    std::fs::write(&custom, &text).unwrap();
    let (a, b) = (d.join("a"), d.join("b"));
    for (case, out) in [(&plain, &a), (&custom, &b)] {
        let o = foehn(&["--json", "--quiet", "flow", s(case), "--out", s(out)]);
        assert_eq!(
            o.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&o.stderr)
        );
    }
    // the results are the same bytes with and without the columns
    for f in ["fields.vtkhdf", "residuals.csv", "balance.csv"] {
        if a.join(f).exists() {
            assert!(bytes(&a.join(f)) == bytes(&b.join(f)), "{f} differs");
        }
    }
    let log = std::fs::read_to_string(b.join("flow.log")).unwrap();
    assert!(log.contains("probe mid"), "{log}");
    assert!(log.contains("extremes"), "{log}");
    assert!(log.contains("Umax") && log.contains("Pmin"), "{log}");
    let rec = std::fs::read_to_string(b.join("run.jsonl")).unwrap();
    let rows: Vec<Value> = rec
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(rows[0]["format_version"], "1.1");
    let cols = rows
        .iter()
        .find(|r| r["type"] == "columns")
        .expect("a columns row");
    let ids: Vec<&str> = cols["columns"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c["id"].as_str().unwrap())
        .collect();
    assert!(ids.contains(&"max.velocity_magnitude.cell"), "{ids:?}");
    let row = rows
        .iter()
        .find(|r| r["type"] == "table_row")
        .expect("a table_row row");
    assert!(row["values"]["max.velocity_magnitude.cell"]["cell"].is_u64());
    let sections: Vec<&str> = rows
        .iter()
        .filter(|r| r["type"] == "section")
        .map(|r| r["title"].as_str().unwrap())
        .collect();
    assert!(sections.contains(&"Fluid properties"), "{sections:?}");
    // the default table gets its columns row as well
    let rec = std::fs::read_to_string(a.join("run.jsonl")).unwrap();
    assert!(rec.contains("\"type\":\"columns\""));
    // a wrong name is refused before the run
    let bad = d.join("bad.toml");
    std::fs::write(
        &bad,
        format!("{text}\n").replace("min.pressure", "min.speed"),
    )
    .unwrap();
    let o = foehn(&["flow", s(&bad), "--check"]);
    assert_ne!(o.status.code(), Some(0));
    let err = format!(
        "{}{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(err.contains("`min.speed`: unknown quantity"), "{err}");
}

/// `[[schedule]]` (D128): every change of a scheduled setting is an `event`
/// row of the run record with its key, old and new value, and a case
/// without a schedule has none.
#[test]
fn schedule_changes_are_events_of_the_run_record() {
    let d = dir("schedule");
    let case = steady(&d);
    let mut text = std::fs::read_to_string(&case).unwrap();
    text.push_str(
        "[[schedule]]\nat_iteration = 20\nset = { solver.momentum.relaxation = 0.6, solver.log_every = 5 }\n",
    );
    std::fs::write(&case, text).unwrap();
    let out = d.join("out");
    let o = foehn(&["--json", "--quiet", "flow", s(&case), "--out", s(&out)]);
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let rec = std::fs::read_to_string(out.join("run.jsonl")).unwrap();
    let rows: Vec<Value> = rec
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let ev: Vec<&Value> = rows.iter().filter(|r| r["code"] == "RUN015").collect();
    assert_eq!(ev.len(), 2, "{rec}");
    assert!(ev.iter().all(|e| e["type"] == "event" && e["at"] == 20));
    let k = ev
        .iter()
        .find(|e| e["key"] == "solver.momentum.relaxation")
        .unwrap();
    assert_eq!(k["old"], "0.9");
    assert_eq!(k["new"], "0.6");
    let log = std::fs::read_to_string(out.join("flow.log")).unwrap();
    assert!(
        log.contains("iteration 20: schedule (iteration 20)"),
        "{log}"
    );
    // the plain run has no such events
    let plain = d.join("plain");
    let c0 = d.join("plain.toml");
    std::fs::write(&c0, std::fs::read_to_string(steady(&d)).unwrap()).unwrap();
    let o = foehn(&["--json", "--quiet", "flow", s(&c0), "--out", s(&plain)]);
    assert_eq!(o.status.code(), Some(0));
    let rec = std::fs::read_to_string(plain.join("run.jsonl")).unwrap();
    assert!(!rec.contains("RUN015"));
}
