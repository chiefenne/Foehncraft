//! The tool registry from the command line (D196): `foehn commands` lists
//! every tool, and `foehn call <id>` gives the same object as the clap
//! command of the same tool with `--json`.

use serde_json::{Value, json};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const FOEHN: &str = env!("CARGO_BIN_EXE_foehn");

fn dir(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("tools_contract")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

const CAVITY: &str = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.0625]
cells = [8, 8, 1]

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

fn foehn(args: &[&str]) -> Output {
    Command::new(FOEHN).args(args).output().expect("foehn runs")
}

fn json_of(o: &Output) -> Value {
    serde_json::from_slice(&o.stdout).unwrap_or_else(|e| {
        panic!(
            "stdout is not JSON ({e}):\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        )
    })
}

/// `foehn call id` with `input` on standard input: (exit code, envelope).
fn call(id: &str, input: &Value) -> (i32, Value) {
    let mut c = Command::new(FOEHN)
        .args(["--quiet", "call", id, "--input", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("foehn runs");
    c.stdin
        .take()
        .unwrap()
        .write_all(input.to_string().as_bytes())
        .unwrap();
    let o = c.wait_with_output().unwrap();
    let v = json_of(&o);
    let code = o.status.code().unwrap();
    assert_eq!(v["format"], "foehn.tool-result");
    assert_eq!(v["tool"], id);
    assert_eq!(v["exit_code"], code, "{v}");
    (code, v)
}

/// The `--json` object of a clap command: (exit code, object).
fn cli(args: &[&str]) -> (i32, Value) {
    let mut a = vec!["--json", "--quiet"];
    a.extend_from_slice(args);
    let o = foehn(&a);
    (o.status.code().unwrap(), json_of(&o))
}

#[test]
fn commands_lists_every_tool_with_schemas() {
    let o = foehn(&["commands"]);
    assert!(o.status.success());
    let v = json_of(&o);
    assert_eq!(v["format"], "foehn.tool-catalogue");
    let tools = v["tools"].as_array().unwrap();
    for id in [
        "system.info",
        "system.commands",
        "case.check",
        "case.resolve",
        "diag.codes",
        "materials.list",
        "expression.eval",
    ] {
        let t = tools
            .iter()
            .find(|t| t["id"] == id)
            .unwrap_or_else(|| panic!("{id} missing"));
        assert_eq!(t["input_schema"]["type"], "object", "{id}");
        assert!(t["output_schema"].is_object(), "{id}");
        assert!(t["action"].is_string(), "{id}");
    }
    let (code, e) = call("system.commands", &Value::Null);
    assert_eq!(code, 0);
    assert_eq!(e["result"], v);
}

#[test]
fn call_matches_the_clap_commands() {
    let (c1, a) = cli(&["info"]);
    let (c2, b) = call("system.info", &json!({}));
    assert_eq!((c1, &a), (c2, &b["result"]));

    let (_, a) = cli(&["bench", "backends"]);
    let (_, b) = call("system.backends", &json!({}));
    assert_eq!(a["backends"], b["result"]["backends"]);

    let (_, a) = cli(&["schema"]);
    let (_, b) = call("case.schema", &json!({}));
    assert_eq!(a, b["result"]);

    let (_, a) = cli(&["materials", "--temperature", "350"]);
    let (_, b) = call("materials.list", &json!({"temperature": 350.0}));
    assert_eq!(a, b["result"]);

    let (_, a) = cli(&["materials", "water"]);
    let (_, b) = call("materials.show", &json!({"name": "water"}));
    assert_eq!(a, b["result"]);

    let (_, a) = cli(&["diag", "codes"]);
    let (_, b) = call("diag.codes", &json!({}));
    assert_eq!(a["codes"], b["result"]["codes"]);

    let (_, a) = cli(&["expression", "x*x + t", "--point", "2,0,0,1"]);
    let (_, b) = call(
        "expression.eval",
        &json!({"source": "x*x + t", "points": [[2.0, 0.0, 0.0, 1.0]]}),
    );
    assert_eq!(a["values"], json!([5.0]));
    assert_eq!(a["values"], b["result"]["values"]);
    assert_eq!(a["uses"], b["result"]["uses"]);
}

#[test]
fn check_and_resolve_match() {
    let d = dir("check");
    let case = d.join("cavity.toml");
    std::fs::write(&case, CAVITY).unwrap();
    let c = case.to_str().unwrap();
    let (c1, a) = cli(&["flow", c, "--check", "--set", "fluid.viscosity=0.02"]);
    let (c2, b) = call(
        "case.check",
        &json!({"case": c, "set": ["fluid.viscosity=0.02"]}),
    );
    assert_eq!(c1, 0);
    assert_eq!((c1, &a), (c2, &b["result"]));
    assert_eq!(a["ok"], true);

    // a case that cannot run: exit code 1 from both
    let (c1, a) = cli(&["flow", c, "--check", "--set", "fluid.density=-1"]);
    let (c2, b) = call(
        "case.check",
        &json!({"case": c, "set": ["fluid.density=-1"]}),
    );
    assert_eq!(c1, 1);
    assert_eq!((c1, &a), (c2, &b["result"]));
    assert_eq!(b["status"], "error");

    let (c1, a) = cli(&["case", "resolve", c, "--set", "fluid.viscosity=0.02"]);
    let (c2, b) = call(
        "case.resolve",
        &json!({"case": c, "set": ["fluid.viscosity=0.02"]}),
    );
    assert_eq!(c1, 0);
    assert_eq!(a["toml"], b["result"]["toml"]);
    assert_eq!(c2, 0);
    assert!(a["toml"].as_str().unwrap().contains("viscosity = 0.02"));
}

#[test]
fn bad_calls_are_diagnostics() {
    let (code, e) = call("case.chek", &json!({}));
    assert_eq!(code, 1);
    assert_eq!(e["diagnostics"][0]["code"], "TOOL001");
    let (code, e) = call("case.check", &json!({"cas": "x.toml"}));
    assert_eq!(code, 1);
    assert_eq!(e["diagnostics"][0]["code"], "TOOL002");
    let (code, e) = call("case.check", &json!({"case": "missing.toml"}));
    // the check never fails: the missing file is a diagnostic of the report
    assert_eq!(code, 1);
    assert_eq!(e["result"]["ok"], false);
}

#[test]
fn markdown_is_docs_tools_md() {
    let o = foehn(&["commands", "--markdown"]);
    let doc = include_str!("../../../docs/tools.md").replace("\r\n", "\n");
    assert_eq!(
        String::from_utf8_lossy(&o.stdout).replace("\r\n", "\n"),
        doc,
        "regenerate docs/tools.md with foehn commands --markdown"
    );
}

#[test]
fn runs_and_files_match() {
    let d = dir("runs");
    let case = d.join("cavity.toml");
    std::fs::write(&case, CAVITY).unwrap();
    let c = case.to_str().unwrap();
    let a_out = d.join("cli");
    let b_out = d.join("call");
    let (c1, a) = cli(&["flow", c, "--out", a_out.to_str().unwrap()]);
    let (c2, b) = call("flow.run", &json!({"case": c, "out": b_out}));
    assert_eq!(c1, c2);
    let b = &b["result"];
    for k in [
        "format",
        "status",
        "exit_code",
        "solves",
        "orders",
        "diagnostics",
    ] {
        assert_eq!(a[k], b[k], "{k}");
    }
    let bytes = |p: &Path| std::fs::read(p).unwrap();
    assert_eq!(
        bytes(&a_out.join("restart.h5")),
        bytes(&b_out.join("restart.h5")),
        "the same solution, bitwise"
    );

    // the envelope carries the limitations of the run (D199): one mesh
    let (_, env) = call("flow.run", &json!({"case": c, "out": b_out}));
    assert_eq!(env["limitations"], env["result"]["limitations"]);
    assert!(
        env["limitations"]
            .as_array()
            .unwrap()
            .iter()
            .any(|l| l["code"] == "LIM004"),
        "{}",
        env["limitations"]
    );
    let (_, r) = call("result.read", &json!({"dir": b_out}));
    assert_eq!(r["result"]["solves"], b["solves"]);
    assert_eq!(r["limitations"], env["limitations"]);
    let (_, i) = call("system.info", &json!({}));
    assert_eq!(i["limitations"], json!([]));

    let (_, e) = call("case.estimate", &json!({"case": c}));
    let (_, k) = cli(&["flow", c, "--check"]);
    assert_eq!(e["result"], k["estimate"]);

    let (c1, a) = cli(&[
        "post",
        "derive",
        c,
        "--out",
        a_out.to_str().unwrap(),
        "--fields",
        "vorticity",
        "--file",
        d.join("a.vtkhdf").to_str().unwrap(),
    ]);
    let (c2, b) = call(
        "post.derive",
        &json!({"case": c, "out": a_out, "fields": ["vorticity"], "file": d.join("b.vtkhdf")}),
    );
    assert_eq!((c1, c2), (0, 0));
    assert_eq!(a["derived"], b["result"]["derived"]);
    assert_eq!(bytes(&d.join("a.vtkhdf")), bytes(&d.join("b.vtkhdf")));

    let (c1, a) = cli(&["case", "upgrade", c]);
    let (c2, b) = call("case.upgrade", &json!({"files": [c]}));
    assert_eq!(c1, c2);
    assert_eq!(a["files"], b["result"]["files"]);

    let mesh = d.join("mesh.toml");
    std::fs::write(
        &mesh,
        "[mesh.generate]\ntype = \"hex\"\ndimension = 2\ncells = [4, 4, 1]\n",
    )
    .unwrap();
    let m = mesh.to_str().unwrap();
    let (c1, a) = cli(&["mesh", m, "--out", d.join("ma").to_str().unwrap()]);
    let (c2, b) = call("mesh.process", &json!({"input": m, "out": d.join("mb")}));
    assert_eq!((c1, c2), (0, 0));
    assert_eq!(a["valid"], b["result"]["valid"]);
}

#[test]
fn studies_match() {
    let d = dir("studies");
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let scalar = root.join("examples/cases/scalar/mms_poly_2d.toml");
    let s = scalar.to_str().unwrap();
    let (c1, a) = cli(&["scalar", s, "--out", d.join("a").to_str().unwrap()]);
    let (c2, b) = call("study.scalar", &json!({"case": s, "out": d.join("b")}));
    assert_eq!(c1, c2);
    assert_eq!(
        a["orders"], b["result"]["orders"],
        "bitwise the same orders"
    );
    assert_eq!(a["failures"], b["result"]["failures"]);

    let g1 = root.join("examples/cases/bench/g1_ci_small.toml");
    let g = g1.to_str().unwrap();
    let (c1, a) = cli(&["bench", "g1", g, "--out", d.join("ga").to_str().unwrap()]);
    let (c2, b) = call("study.g1", &json!({"case": g, "out": d.join("gb")}));
    assert_eq!(c1, c2);
    let pass = |v: &Value| -> Vec<Value> {
        v["verdicts"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| json!([x["candidate"], x["reference"]]))
            .collect()
    };
    assert_eq!(pass(&a), pass(&b["result"]));
}

/// D197: `--journal` records every call with its typed input, a run writes
/// its call into `run.jsonl`, and `foehn run --script` replays both to the
/// same files.
#[test]
fn journal_and_replay() {
    let d = dir("journal");
    let case = d.join("cavity.toml");
    std::fs::write(&case, CAVITY).unwrap();
    let c = case.to_str().unwrap();
    let j = d.join("session.jsonl");
    let js = j.to_str().unwrap();
    let out = d.join("run");
    let o = out.to_str().unwrap();
    assert!(
        foehn(&["--quiet", "--journal", js, "flow", c, "--check"])
            .status
            .success()
    );
    assert!(
        foehn(&["--quiet", "--journal", js, "flow", c, "--out", o])
            .status
            .success()
    );
    let first = std::fs::read(out.join("restart.h5")).unwrap();
    let rows: Vec<Value> = std::fs::read_to_string(&j)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(rows[0]["format"], "foehn.journal");
    assert_eq!(rows[1]["tool"], "case.check");
    assert_eq!(rows[2]["tool"], "flow.run");
    assert_eq!(rows[2]["via"], "cli");
    assert_eq!(rows[2]["input"]["out"], o);
    assert_eq!(
        rows[2]["input"]["no_vtkhdf"], false,
        "typed: defaults filled in"
    );
    assert_eq!(rows[2]["exit_code"], 0);

    // the run record holds the call that started the run
    let record = std::fs::read_to_string(out.join("run.jsonl")).unwrap();
    let call_row = record
        .lines()
        .map(|l| serde_json::from_str::<Value>(l).unwrap())
        .find(|v| v["type"] == "call")
        .expect("a call row in run.jsonl");
    assert_eq!(call_row["tool"], "flow.run");
    assert_eq!(call_row["input"], rows[2]["input"]);

    std::fs::remove_dir_all(&out).unwrap();
    let (code, r) = cli(&["run", "--script", js]);
    assert_eq!(code, 0, "{r}");
    assert_eq!(r["format"], "foehn.replay");
    assert_eq!(
        (r["replayed"].clone(), r["different"].clone()),
        (json!(2), json!(0))
    );
    assert_eq!(
        std::fs::read(out.join("restart.h5")).unwrap(),
        first,
        "the replayed run, bitwise"
    );
    // a replay without --journal adds nothing to the journal
    let rows_after = std::fs::read_to_string(&j).unwrap().lines().count();
    assert_eq!(
        rows_after,
        rows.len(),
        "the replay without --journal adds nothing"
    );

    let rec = out.join("run.jsonl");
    let (code, r) = cli(&["run", "--script", rec.to_str().unwrap()]);
    assert_eq!(code, 0, "{r}");
    assert_eq!(r["calls"][0]["tool"], "flow.run");
    assert_eq!(std::fs::read(out.join("restart.h5")).unwrap(), first);

    let (code, r) = call(
        "journal.replay",
        &json!({"script": d.join("missing.jsonl")}),
    );
    assert_eq!(code, 1);
    assert_eq!(r["diagnostics"][0]["code"], "TOOL004");
}
