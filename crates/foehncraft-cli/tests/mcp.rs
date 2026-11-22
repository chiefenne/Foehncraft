//! The MCP server `foehn mcp` (L18): the protocol on standard input and
//! output, the roots it needs, the action table (D198) with its questions
//! (MCP elicitation), cancellation, paging and the journal a person replays
//! with `foehn run --script` (D197).

use serde_json::{Value, json};
use std::io::{BufRead as _, BufReader, Write as _};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

const FOEHN: &str = env!("CARGO_BIN_EXE_foehn");

fn dir(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("mcp")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("runs")).unwrap();
    d
}

/// A cavity of 8 × 8 cells; `limits` adds `[run.limits]`.
fn cavity(limits: bool, iterations: usize) -> String {
    let mut s = format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.0625]
cells = [8, 8, 1]

[fluid]
density = 1.0
viscosity = 0.01

[solver]
max_iterations = {iterations}
tolerance = 1e-14

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
"#
    );
    if limits {
        s.push_str("\n[run.limits]\nmax_wall_time = 600.0\nmax_output = 1e9\n");
    }
    s
}

struct Client {
    child: Child,
    input: Option<ChildStdin>,
    output: BufReader<ChildStdout>,
    next: u64,
    /// Answers to the questions of the server, in order.
    answers: Vec<&'static str>,
    /// The questions the server asked.
    asked: Vec<String>,
    /// Progress notifications seen.
    progress: usize,
}

impl Client {
    fn start(project: &Path, elicitation: bool) -> Client {
        let mut child = Command::new(FOEHN)
            .args(["--quiet", "mcp", "--read-root", ".", "--write-root", "runs"])
            .current_dir(project)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .expect("foehn mcp starts");
        let input = child.stdin.take();
        let output = BufReader::new(child.stdout.take().unwrap());
        let mut c = Client {
            child,
            input,
            output,
            next: 0,
            answers: Vec::new(),
            asked: Vec::new(),
            progress: 0,
        };
        let caps = if elicitation {
            json!({"elicitation": {}})
        } else {
            json!({})
        };
        let r = c.request(
            "initialize",
            json!({"protocolVersion": "2025-06-18", "capabilities": caps,
                   "clientInfo": {"name": "test", "version": "0"}}),
        );
        assert_eq!(r["result"]["protocolVersion"], "2025-06-18");
        assert_eq!(r["result"]["serverInfo"]["name"], "foehn");
        c.send(&json!({"jsonrpc": "2.0", "method": "notifications/initialized"}));
        c
    }

    fn send(&mut self, v: &Value) {
        let i = self.input.as_mut().unwrap();
        writeln!(i, "{v}").unwrap();
        i.flush().unwrap();
    }

    fn read(&mut self) -> Value {
        let mut l = String::new();
        assert!(
            self.output.read_line(&mut l).unwrap() > 0,
            "the server ended"
        );
        serde_json::from_str(&l).unwrap()
    }

    fn post(&mut self, method: &str, params: Value) -> u64 {
        self.next += 1;
        let id = self.next;
        self.send(&json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params}));
        id
    }

    /// Waits for the response `id`, answering questions on the way; `on`
    /// sees every notification.
    fn wait(&mut self, id: u64, mut on: impl FnMut(&mut Client, &Value)) -> Value {
        loop {
            let m = self.read();
            if m["method"] == "elicitation/create" {
                self.asked
                    .push(m["params"]["message"].as_str().unwrap().to_string());
                let a = if self.answers.is_empty() {
                    "refuse-once"
                } else {
                    self.answers.remove(0)
                };
                let r = json!({"jsonrpc": "2.0", "id": m["id"],
                               "result": {"action": "accept", "content": {"answer": a}}});
                self.send(&r);
            } else if m["method"] == "notifications/progress" {
                self.progress += 1;
                on(self, &m);
            } else if m["id"] == id {
                return m;
            }
        }
    }

    fn request(&mut self, method: &str, params: Value) -> Value {
        let id = self.post(method, params);
        self.wait(id, |_, _| {})
    }

    /// Calls a tool and returns the envelope (`structuredContent`).
    fn call(&mut self, name: &str, args: Value) -> Value {
        let r = self.request(
            "tools/call",
            json!({"name": name, "arguments": args, "_meta": {"progressToken": 1}}),
        );
        let res = &r["result"];
        let text = res["content"][0]["text"].as_str().unwrap();
        let env = res["structuredContent"].clone();
        assert_eq!(env["tool"], name, "{text}");
        assert!(text.starts_with(&format!("{name}: {}", env["status"].as_str().unwrap())));
        assert_eq!(res["isError"], env["status"] == "error");
        env
    }

    fn end(mut self) -> i32 {
        drop(self.input.take());
        self.child.wait().unwrap().code().unwrap_or(-1)
    }
}

/// The run ran (the cavity does not reach its tolerance of 1e-14: failed
/// checks, exit code 2).
fn ran(env: &Value) {
    assert_eq!(env["status"], "failed-checks", "{env}");
    assert_eq!(env["exit_code"], 2);
}

fn code(env: &Value) -> &str {
    env["diagnostics"][0]["code"].as_str().unwrap_or("")
}

fn journal(project: &Path) -> Vec<Value> {
    std::fs::read_to_string(project.join("foehn-mcp-journal.jsonl"))
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .filter(|v: &Value| v["type"] == "call")
        .collect()
}

#[test]
fn refuses_to_start_without_both_roots() {
    let d = dir("roots");
    for args in [
        vec!["mcp"],
        vec!["mcp", "--read-root", "."],
        vec!["mcp", "--write-root", "runs"],
        vec![
            "mcp",
            "--read-root",
            ".",
            "--write-root",
            ".",
            "--policy",
            "foehn-mcp.toml",
        ],
    ] {
        let o = Command::new(FOEHN)
            .args(&args)
            .current_dir(&d)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(o.status.code(), Some(1), "{args:?}");
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(
            err.contains("--read-root and --write-root") || err.contains("inside a write root"),
            "{err}"
        );
    }
}

#[test]
fn the_tools_are_the_registry_ids() {
    let d = dir("list");
    let mut c = Client::start(&d, false);
    let l = c.request("tools/list", json!({}));
    let names: Vec<String> = l["result"]["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t["name"].as_str().unwrap().to_string())
        .collect();
    let cat = Command::new(FOEHN).arg("commands").output().unwrap();
    let cat: Value = serde_json::from_slice(&cat.stdout).unwrap();
    let ids: Vec<String> = cat["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t["id"].as_str().unwrap().to_string())
        .collect();
    assert_eq!(names, ids);
    let check = &l["result"]["tools"][4];
    assert_eq!(check["name"], "case.check");
    assert_eq!(check["inputSchema"]["type"], "object");
    assert_eq!(check["annotations"]["readOnlyHint"], true);
    // an unknown tool is a diagnostic, journaled
    let env = c.call("case.chek", json!({}));
    assert_eq!(code(&env), "TOOL001");
    assert_eq!(c.end(), 0);
    assert_eq!(journal(&d)[0]["via"], "mcp");
}

#[test]
fn a_session_with_questions_refusals_paging_and_its_replay() {
    let d = dir("session");
    std::fs::write(d.join("cavity.toml"), cavity(false, 30)).unwrap();
    let mut c = Client::start(&d, true);

    let env = c.call(
        "case.check",
        json!({"case": "cavity.toml", "out": "runs/a"}),
    );
    assert_eq!(env["status"], "finished", "{env}");

    // no [run.limits]: refused, with the fix
    let env = c.call("flow.run", json!({"case": "cavity.toml", "out": "runs/a"}));
    assert_eq!(code(&env), "MCP002");
    let fix = env["diagnostics"][0]["fix"]["set"].as_array().unwrap();
    assert_eq!(fix.len(), 2);
    assert!(c.asked.is_empty(), "no question for a refused call");

    // with limits: the server asks; "allow always" goes into the table
    let limits = json!(["run.limits.max_wall_time=600", "run.limits.max_output=1e9"]);
    c.answers = vec!["allow-always"];
    let env = c.call(
        "flow.run",
        json!({"case": "cavity.toml", "out": "runs/a", "set": limits}),
    );
    ran(&env);
    assert_eq!(c.asked.len(), 1);
    assert!(c.asked[0].contains("start-run"), "{}", c.asked[0]);
    assert!(c.asked[0].contains("Estimate: 64 cells"), "{}", c.asked[0]);
    assert!(env["result"]["limitations"].is_array());
    let table = std::fs::read_to_string(d.join("foehn-mcp.toml")).unwrap();
    assert!(table.contains("start-run = \"allow\""), "{table}");

    // a new directory: no question now
    let env = c.call(
        "flow.run",
        json!({"case": "cavity.toml", "out": "runs/b", "set": limits}),
    );
    ran(&env);
    assert_eq!(c.asked.len(), 1);

    // over earlier results: asked, refused once
    let env = c.call(
        "flow.run",
        json!({"case": "cavity.toml", "out": "runs/a", "set": limits}),
    );
    assert_eq!(code(&env), "MCP004");
    assert!(c.asked[1].contains("overwrite-results"), "{}", c.asked[1]);

    // outside the write roots: refused by the table
    let env = c.call(
        "flow.run",
        json!({"case": "cavity.toml", "out": "elsewhere", "set": limits}),
    );
    assert_eq!(code(&env), "MCP003");
    assert!(!d.join("elsewhere").exists());

    // paging: the series, then the residual history by cursor
    let env = c.call("result.page", json!({"dir": "runs/a"}));
    let names: Vec<&str> = env["result"]["available"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s["name"].as_str().unwrap())
        .collect();
    assert!(
        names.contains(&"residuals") && names.contains(&"record"),
        "{names:?}"
    );
    let p1 = c.call(
        "result.page",
        json!({"dir": "runs/a", "series": "residuals", "limit": 20}),
    );
    assert_eq!(p1["result"]["total"], 30);
    assert_eq!(p1["result"]["rows"].as_array().unwrap().len(), 20);
    let cursor = p1["result"]["next_cursor"].as_str().unwrap().to_string();
    let p2 = c.call(
        "result.page",
        json!({"dir": "runs/a", "series": "residuals", "limit": 20, "cursor": cursor}),
    );
    assert_eq!(p2["result"]["first"], 20);
    assert_eq!(p2["result"]["rows"].as_array().unwrap().len(), 10);
    assert!(p2["result"]["next_cursor"].is_null());
    let bad = c.call("result.page", json!({"dir": "runs/a", "series": "nothing"}));
    assert_eq!(code(&bad), "TOOL006");
    assert_eq!(c.end(), 0);

    // the journal: every call, the answers, the refusals
    let j = journal(&d);
    assert_eq!(j.len(), 10);
    assert!(j.iter().all(|r| r["via"] == "mcp"));
    assert_eq!(j[2]["approval"]["answer"], "allow-always");
    assert_eq!(j[3]["approval"]["by"], "action table");
    assert_eq!(j[4]["diagnostic"]["code"], "MCP004");

    // a person replays the session without a window
    let o = Command::new(FOEHN)
        .args([
            "--json",
            "run",
            "--script",
            "foehn-mcp-journal.jsonl",
            "--read-root",
            ".",
            "--write-root",
            "runs",
        ])
        .current_dir(&d)
        .output()
        .unwrap();
    let r: Value = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(o.status.code(), Some(0), "{r}");
    assert_eq!(r["skipped"], 3, "{r}");
    assert_eq!(r["different"], 0, "{r}");

    // results through MCP are the results of the command line
    let o = Command::new(FOEHN)
        .args([
            "--quiet",
            "flow",
            "cavity.toml",
            "--out",
            "cli",
            "--set",
            "run.limits.max_wall_time=600",
            "--set",
            "run.limits.max_output=1e9",
        ])
        .current_dir(&d)
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert_eq!(
        std::fs::read(d.join("cli/residuals.csv")).unwrap(),
        std::fs::read(d.join("runs/b/residuals.csv")).unwrap()
    );
}

#[test]
fn a_host_that_cannot_ask_gets_a_refusal() {
    let d = dir("no_ask");
    std::fs::write(d.join("cavity.toml"), cavity(true, 5)).unwrap();
    let mut c = Client::start(&d, false);
    let env = c.call("flow.run", json!({"case": "cavity.toml", "out": "runs/a"}));
    assert_eq!(code(&env), "MCP005");
    assert!(!d.join("runs/a").exists());
    assert_eq!(c.end(), 0);
    // the person allows runs in the table
    std::fs::write(
        d.join("foehn-mcp.toml"),
        "[actions]\nstart-run = \"allow\"\n",
    )
    .unwrap();
    let mut c = Client::start(&d, false);
    let env = c.call("flow.run", json!({"case": "cavity.toml", "out": "runs/a"}));
    ran(&env);
    assert_eq!(c.end(), 0);
}

#[test]
fn a_cancelled_run_stops_with_its_output() {
    let d = dir("cancel");
    std::fs::write(d.join("cavity.toml"), cavity(true, 2000)).unwrap();
    std::fs::write(
        d.join("foehn-mcp.toml"),
        "[actions]\nstart-run = \"allow\"\n",
    )
    .unwrap();
    let mut c = Client::start(&d, false);
    let id = c.post(
        "tools/call",
        json!({"name": "flow.run", "arguments": {"case": "cavity.toml", "out": "runs/a"},
               "_meta": {"progressToken": "run"}}),
    );
    let mut sent = false;
    let r = c.wait(id, |c, _| {
        if !sent {
            sent = true;
            c.send(
                &json!({"jsonrpc": "2.0", "method": "notifications/cancelled",
                           "params": {"requestId": id, "reason": "test"}}),
            );
        }
    });
    let env = &r["result"]["structuredContent"];
    assert!(c.progress >= 1);
    assert_eq!(env["status"], "stopped", "{env}");
    assert_eq!(env["exit_code"], 3);
    assert!(d.join("runs/a/result.json").exists());
    assert_eq!(c.end(), 0);
}

/// L19: while `flow.run` is going, `run.status` answers at once with the
/// live residuals, `run.set` changes a run-time setting (and cannot touch the
/// limits), `run.stop` ends the run with its output; the journal holds the
/// live calls and a replay skips them.
#[test]
fn a_running_case_is_steered_and_stopped_through_its_control_port() {
    let d = dir("control");
    std::fs::write(d.join("cavity.toml"), cavity(true, 2_000_000)).unwrap();
    std::fs::write(
        d.join("foehn-mcp.toml"),
        "[actions]\nstart-run = \"allow\"\n",
    )
    .unwrap();
    let mut c = Client::start(&d, false);
    let run = c.post(
        "tools/call",
        json!({"name": "flow.run", "arguments": {"case": "cavity.toml", "out": "runs/live"}}),
    );
    // the port opens with the run; wait for the first iterations
    let mut status = Value::Null;
    for _ in 0..600 {
        if d.join("runs/live/control.json").exists() {
            let env = c.call("run.status", json!({"dir": "runs/live"}));
            if env["status"] == "finished" && env["result"]["last"]["iteration"].as_u64() > Some(5)
            {
                status = env["result"].clone();
                break;
            }
        }
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
    assert_eq!(status["format"], "foehn.run-status", "{status}");
    assert_eq!(status["state"], "running");
    assert!(status["last"]["continuity"].is_number());

    // a limit is the person's: refused; a run-time setting is queued
    let env = c.call(
        "run.set",
        json!({"dir": "runs/live", "set": {"run.limits.max_wall_time": 1e9}}),
    );
    assert_eq!(code(&env), "TOOL007", "{env}");
    let env = c.call(
        "run.set",
        json!({"dir": "runs/live", "set": {"solver.momentum.relaxation": 0.5}}),
    );
    assert_eq!(env["status"], "finished", "{env}");
    // the run reports what it did with it
    let mut applied = false;
    for _ in 0..200 {
        let env = c.call("run.status", json!({"dir": "runs/live"}));
        applied = env["result"]["control"]
            .as_array()
            .is_some_and(|l| l.iter().any(|e| e["status"] == "applied"));
        if applied {
            let e = &env["result"]["control"][0];
            assert_eq!(e["key"], "solver.momentum.relaxation");
            assert!(e["schedule"].as_str().unwrap().contains("at_iteration"));
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
    assert!(applied, "the run did not report the change");
    // an invalid value is reported at the boundary, not applied
    c.call(
        "run.set",
        json!({"dir": "runs/live", "set": {"solver.momentum.relaxation": "high"}}),
    );
    let mut rejected = false;
    for _ in 0..200 {
        let env = c.call("run.status", json!({"dir": "runs/live"}));
        rejected = env["result"]["control"]
            .as_array()
            .is_some_and(|l| l.iter().any(|e| e["status"] == "rejected"));
        if rejected {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
    assert!(rejected);

    let env = c.call("run.stop", json!({"dir": "runs/live"}));
    assert_eq!(env["status"], "finished", "{env}");
    let r = c.wait(run, |_, _| {});
    let env = &r["result"]["structuredContent"];
    assert_eq!(env["exit_code"], 3, "{env}");
    assert_eq!(env["result"]["stop"]["code"], "RUN022", "{env}");
    // the port is closed and its files are gone
    assert!(!d.join("runs/live/control.json").exists());
    assert!(!d.join("runs/live/control.token").exists());
    let env = c.call("run.status", json!({"dir": "runs/live"}));
    assert_eq!(code(&env), "TOOL007");
    // the run record says what happened, with the equivalent schedule entry
    let rec = std::fs::read_to_string(d.join("runs/live/run.jsonl")).unwrap();
    for code in ["RUN022", "RUN023", "RUN024"] {
        assert!(rec.contains(code), "{code} missing in the run record");
    }
    assert_eq!(c.end(), 0);

    // the journal holds the live calls; a replay skips them
    let j = journal(&d);
    assert!(j.iter().any(|v| v["tool"] == "run.stop"));
    let script = d.join("live-script.jsonl");
    let rows: Vec<String> = j
        .iter()
        .filter(|v| v["tool"].as_str().is_some_and(|t| t.starts_with("run.")))
        .map(|v| v.to_string())
        .collect();
    std::fs::write(&script, rows.join("\n")).unwrap();
    let o = Command::new(FOEHN)
        .args(["--json", "run", "--script", "live-script.jsonl"])
        .args(["--read-root", ".", "--write-root", "runs"])
        .current_dir(&d)
        .output()
        .unwrap();
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(
        out.contains("\"skipped\": ") || out.contains("\"skipped\":"),
        "{out}"
    );
    assert!(
        out.contains("\"replayed\": 0") || out.contains("\"replayed\":0"),
        "{out}"
    );
}
