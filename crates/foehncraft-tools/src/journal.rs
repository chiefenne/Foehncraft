//! The journal of tool calls and its replay (D197).
//!
//! With a journal file set ([`Ctx::journal`], `foehn --journal FILE`, the
//! environment variable `FOEHN_JOURNAL`), every top-level tool call of a
//! session, by a person or an agent, is appended to it as one JSON line
//! with its typed input (the input after parsing, defaults filled in), its
//! status, exit code and duration. A run started by a call also writes that
//! call into its run record (`run.jsonl`, D69). `journal.replay`
//! (`foehn run --script FILE`) runs the calls of a journal, or of a run
//! record, again in order, without a window, and reports each call whose
//! status or exit code differs from its record.
//!
//! Journal file (`foehn.journal`, JSON lines):
//!
//! ```text
//! {"type":"header","format":"foehn.journal","format_version":1,"version":"0.6.0","commit":"…"}
//! {"type":"call","seq":1,"date":"…","via":"cli","tool":"case.check","tool_version":1,"input":{…},"status":"finished","exit_code":0,"seconds":0.01}
//! ```
//!
//! A call refused before it ran (an unknown id, `TOOL001`, an input that
//! does not match the schema, `TOOL002`, or a call that the MCP server
//! refused or the person declined, `MCP002` to `MCP005`) is journaled with
//! the input as it came and its code under `diagnostic`; replay skips such
//! calls. A call that the person approved carries the answer under
//! `approval`.

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::diag::{self, Diagnostic};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::io::Write as _;
use std::path::{Path, PathBuf};

/// Format of the journal file.
pub const JOURNAL_FORMAT: &str = "foehn.journal";
/// Format of the output of `journal.replay`.
pub const REPLAY_FORMAT: &str = "foehn.replay";

/// The environment variable that sets the journal file of a front end.
pub const JOURNAL_ENV: &str = "FOEHN_JOURNAL";

/// An open journal: the file and the number of calls in it.
pub struct Journal {
    path: PathBuf,
    via: &'static str,
    seq: Option<u64>,
}

impl Journal {
    /// The journal at `path`, written by the front end `via` (`cli`,
    /// `python`, …). The file is created at the first call; calls are
    /// appended to an existing journal.
    pub fn new(path: impl Into<PathBuf>, via: &'static str) -> Self {
        Self {
            path: path.into(),
            via,
            seq: None,
        }
    }

    /// The journal of the environment variable [`JOURNAL_ENV`], if set.
    pub fn from_env(via: &'static str) -> Option<Self> {
        std::env::var_os(JOURNAL_ENV)
            .filter(|v| !v.is_empty())
            .map(|p| Self::new(PathBuf::from(p), via))
    }

    /// The path of the file.
    pub fn path(&self) -> &Path {
        &self.path
    }

    fn append(&mut self, build: &foehncraft_flow::engine::BuildInfo, row: Value) -> Result<()> {
        let seq = match self.seq {
            Some(n) => n,
            None => {
                if let Some(d) = self.path.parent().filter(|d| !d.as_os_str().is_empty()) {
                    std::fs::create_dir_all(d)?;
                }
                let text = std::fs::read_to_string(&self.path).unwrap_or_default();
                let n = text
                    .lines()
                    .filter_map(|l| serde_json::from_str::<Value>(l).ok())
                    .filter(|v| v["type"] == "call")
                    .count() as u64;
                if text.trim().is_empty() {
                    let h = json!({
                        "type": "header",
                        "format": JOURNAL_FORMAT,
                        "format_version": crate::FORMAT_VERSION,
                        "version": foehncraft_core::VERSION,
                        "commit": build.commit,
                    });
                    self.write(&h)?;
                }
                n
            }
        };
        // fixed key order: type, seq, date, via, then the call
        let mut o = serde_json::Map::new();
        o.insert("type".into(), json!("call"));
        o.insert("seq".into(), json!(seq + 1));
        o.insert("date".into(), json!(foehncraft_flow::engine::run_date()));
        o.insert("via".into(), json!(self.via));
        if let Value::Object(m) = row {
            for (k, v) in m {
                o.entry(k).or_insert(v);
            }
        }
        self.write(&Value::Object(o))?;
        self.seq = Some(seq + 1);
        Ok(())
    }

    fn write(&self, v: &Value) -> Result<()> {
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.path)?;
        writeln!(f, "{}", serde_json::to_string(v)?)?;
        Ok(())
    }
}

/// The call of tool `T` with input `input`, as the journal and the run
/// record hold it.
pub(crate) fn call_row(tool: &str, version: u32, input: Value) -> Value {
    json!({"tool": tool, "tool_version": version, "input": input})
}

/// Appends a finished or refused call to the journal of `ctx`, if any. A
/// journal that cannot be written is a warning (`TOOL003`) in the text of
/// the call; the call itself stands.
pub(crate) fn record(
    ctx: &mut Ctx<'_>,
    call: &Value,
    code: i32,
    diagnostic: Option<&Diagnostic>,
    seconds: f64,
) {
    let Some(mut j) = ctx.journal.take() else {
        return;
    };
    let mut row = call.clone();
    row["status"] = json!(crate::status_of(code));
    row["exit_code"] = json!(code);
    row["seconds"] = json!((seconds * 1e4).round() / 1e4);
    if let Some(d) = diagnostic {
        row["diagnostic"] = json!({"code": d.code, "message": d.message});
    }
    if let Some(a) = &ctx.approval {
        row["approval"] = a.clone();
    }
    if let Err(e) = j.append(&ctx.build, row) {
        let msg = format!(
            "{}: the journal {} could not be written: {e}",
            diag::TOOL_JOURNAL,
            j.path.display()
        );
        ctx.line(&foehncraft_flow::console::tagged("warning", &msg));
    }
    ctx.journal = Some(j);
}

/// Journals a call of the tool `tool` that a host refused before it ran
/// (the action table of `foehn mcp`, a declined approval: `MCP002` to
/// `MCP005`), with its input as it came; a replay skips it.
pub fn record_refused(ctx: &mut Ctx<'_>, tool: &str, input: Value, d: &Diagnostic) {
    let version = crate::spec(tool).map_or(0, |s| s.version);
    record(ctx, &call_row(tool, version, input), 1, Some(d), 0.0);
}

/// Input of `journal.replay`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReplayInput {
    /// A journal (`foehn.journal`) or a run record (`run.jsonl`): JSON
    /// lines; the rows with `"type": "call"` are run in order.
    pub script: PathBuf,
    /// Go on after a call that ends differently from its record.
    #[serde(default)]
    pub keep_going: bool,
}

/// A call of the script and how its replay ended.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct ReplayedCall {
    /// Number of the call in the script (its `seq`, else its position).
    pub seq: u64,
    /// Tool id.
    pub tool: String,
    /// Status in the record (`finished` when the script has none).
    pub recorded_status: String,
    /// Exit code in the record (0 when the script has none).
    pub recorded_exit_code: i32,
    /// Status of the replay (`skipped` for a call refused when it was
    /// recorded, `not-run` after the replay stopped).
    pub status: String,
    /// Exit code of the replay.
    pub exit_code: Option<i32>,
    /// Duration of the replayed call.
    pub seconds: Option<f64>,
    /// The replay ended like the record.
    pub same: bool,
    /// The diagnostic of a replayed call that failed.
    pub diagnostic: Option<Diagnostic>,
}

/// The output of `journal.replay`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct ReplayOutput {
    /// `foehn.replay`.
    pub format: &'static str,
    /// The script.
    pub script: String,
    /// Every call of the script.
    pub calls: Vec<ReplayedCall>,
    /// Calls run.
    pub replayed: usize,
    /// Calls refused when recorded, not run.
    pub skipped: usize,
    /// Calls that ended differently from their record.
    pub different: usize,
}

fn script_error(msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::TOOL_SCRIPT, msg))
}

/// The calls of a script.
fn calls_of(path: &Path) -> Result<Vec<Value>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| script_error(format!("reading {}: {e}", path.display())))?;
    let mut calls = Vec::new();
    for (i, l) in text.lines().enumerate() {
        if l.trim().is_empty() {
            continue;
        }
        let v: Value = serde_json::from_str(l).map_err(|e| {
            script_error(format!("{}, line {}: not JSON: {e}", path.display(), i + 1))
        })?;
        if v["type"] == "call" {
            if !v["tool"].is_string() {
                return Err(script_error(format!(
                    "{}, line {}: a call without `tool`",
                    path.display(),
                    i + 1
                )));
            }
            calls.push(v);
        }
    }
    if calls.is_empty() {
        return Err(script_error(format!(
            "{} has no calls (rows with \"type\": \"call\")",
            path.display()
        )));
    }
    Ok(calls)
}

/// `journal.replay`.
pub struct Replay;

impl Tool for Replay {
    const ID: &'static str = "journal.replay";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Replay the tool calls of a journal or a run record in order, without a window, and report each call that ends differently from its record (D197).";
    const CATEGORY: Category = Category::Journal;
    const ACTION: Action = Action::StartRun;
    type Input = ReplayInput;
    type Output = ReplayOutput;

    fn call(a: ReplayInput, ctx: &mut Ctx<'_>) -> Result<ReplayOutput> {
        ctx.roots.check_read(&a.script, "the script")?;
        let calls = calls_of(&a.script)?;
        let mut out = ReplayOutput {
            format: REPLAY_FORMAT,
            script: a.script.display().to_string(),
            calls: Vec::new(),
            replayed: 0,
            skipped: 0,
            different: 0,
        };
        let mut stopped = false;
        for (i, c) in calls.iter().enumerate() {
            let tool = c["tool"].as_str().unwrap_or_default().to_string();
            let seq = c["seq"].as_u64().unwrap_or(i as u64 + 1);
            let recorded_exit_code = c["exit_code"].as_i64().map_or(0, |v| v as i32);
            let recorded_status = c["status"]
                .as_str()
                .map_or_else(|| crate::status_of(recorded_exit_code).into(), String::from);
            let mut r = ReplayedCall {
                seq,
                tool: tool.clone(),
                recorded_status,
                recorded_exit_code,
                status: "not-run".into(),
                exit_code: None,
                seconds: None,
                same: false,
                diagnostic: None,
            };
            let refused = c["diagnostic"]["code"]
                .as_str()
                .is_some_and(diag::refused_before_run);
            if refused {
                r.status = "skipped".into();
                r.same = true;
                out.skipped += 1;
                ctx.line(&format!(
                    "{seq:>4}  {tool:<18} skipped (refused when recorded)"
                ));
            } else if !stopped {
                ctx.line(&format!("{seq:>4}  {tool}"));
                let t0 = std::time::Instant::now();
                let res = crate::call(&tool, c["input"].clone(), ctx);
                r.seconds = Some((t0.elapsed().as_secs_f64() * 1e4).round() / 1e4);
                r.status = res.status.into();
                r.exit_code = Some(res.exit_code);
                r.same = res.exit_code == r.recorded_exit_code && r.status == r.recorded_status;
                r.diagnostic = res.diagnostics.into_iter().next();
                out.replayed += 1;
                if r.same {
                    ctx.line(&format!(
                        "      {} (exit {}), as recorded",
                        r.status, res.exit_code
                    ));
                } else {
                    out.different += 1;
                    let msg = format!(
                        "{}: call {seq} ({tool}) ended {} (exit {}), recorded {} (exit {})",
                        diag::TOOL_REPLAY_DIFFERS,
                        r.status,
                        res.exit_code,
                        r.recorded_status,
                        r.recorded_exit_code
                    );
                    ctx.line(&foehncraft_flow::console::tagged("warning", &msg));
                    if let Some(d) = &r.diagnostic {
                        ctx.line(&format!("      {}: {}", d.code, d.message));
                    }
                    stopped = !a.keep_going;
                }
            }
            out.calls.push(r);
        }
        ctx.line(&format!(
            "replayed {} of {} calls, {} skipped, {} different",
            out.replayed,
            calls.len(),
            out.skipped,
            out.different
        ));
        Ok(out)
    }

    /// 0 when every call ended as recorded, else 1.
    fn exit_code(o: &ReplayOutput) -> i32 {
        i32::from(o.different > 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use foehncraft_flow::engine::BuildInfo;

    fn dir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("foehn-journal-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn rows(p: &Path) -> Vec<Value> {
        std::fs::read_to_string(p)
            .unwrap()
            .lines()
            .map(|l| serde_json::from_str(l).unwrap())
            .collect()
    }

    fn ctx(j: &Path) -> Ctx<'static> {
        Ctx::new(BuildInfo::default()).journal(Some(Journal::new(j, "test")))
    }

    #[test]
    fn calls_are_journaled_with_their_typed_input() {
        let d = dir("typed");
        let j = d.join("j.jsonl");
        let mut c = ctx(&j);
        let r = crate::call(
            "expression.eval",
            json!({"source": "2*x", "points": [[1.0, 0.0, 0.0, 0.0]]}),
            &mut c,
        );
        assert_eq!(r.exit_code, 0);
        crate::call("case.chek", Value::Null, &mut c);
        crate::call("diag.codes", json!({"bogus": 1}), &mut c);
        // a second context appends and counts on
        crate::call("system.backends", Value::Null, &mut ctx(&j));
        let r = rows(&j);
        assert_eq!(r[0]["format"], JOURNAL_FORMAT);
        let calls = &r[1..];
        assert_eq!(calls.len(), 4);
        let seq: Vec<u64> = calls.iter().map(|c| c["seq"].as_u64().unwrap()).collect();
        assert_eq!(seq, [1, 2, 3, 4]);
        assert_eq!(calls[0]["tool"], "expression.eval");
        assert_eq!(calls[0]["via"], "test");
        assert_eq!(calls[0]["input"]["allowed"], "xyzt", "defaults filled in");
        assert_eq!(calls[0]["status"], "finished");
        assert_eq!(calls[1]["diagnostic"]["code"], diag::TOOL_UNKNOWN);
        assert_eq!(calls[2]["diagnostic"]["code"], diag::TOOL_INPUT);
        assert_eq!(calls[2]["input"], json!({"bogus": 1}), "as it came");
        assert_eq!(calls[3]["tool"], "system.backends");
    }

    #[test]
    fn replay_runs_the_calls_and_skips_the_refused_ones() {
        let d = dir("replay");
        let j = d.join("j.jsonl");
        let mut c = ctx(&j);
        crate::call("materials.show", json!({"name": "water"}), &mut c);
        crate::call("case.chek", Value::Null, &mut c);
        let mut plain = Ctx::new(BuildInfo::default());
        let r = crate::call("journal.replay", json!({"script": j}), &mut plain);
        assert_eq!(r.exit_code, 0, "{:?}", r.diagnostics);
        assert_eq!(r.result["replayed"], 1);
        assert_eq!(r.result["skipped"], 1);
        assert_eq!(r.result["calls"][0]["same"], true);
        // a replay in a journaled context is one call: the inner calls are
        // not journaled again
        let j2 = d.join("j2.jsonl");
        crate::call("journal.replay", json!({"script": j}), &mut ctx(&j2));
        let r2 = rows(&j2);
        assert_eq!(r2.len(), 2);
        assert_eq!(r2[1]["tool"], "journal.replay");
    }

    #[test]
    fn replay_stops_where_it_differs_unless_keep_going() {
        let d = dir("differs");
        let s = d.join("script.jsonl");
        let lines = [
            json!({"type": "call", "tool": "diag.codes", "input": {}, "status": "failed-checks", "exit_code": 2}),
            json!({"type": "call", "tool": "system.backends", "input": {}}),
        ];
        let text: Vec<String> = lines.iter().map(Value::to_string).collect();
        std::fs::write(&s, text.join("\n")).unwrap();
        let mut c = Ctx::new(BuildInfo::default());
        let r = crate::call("journal.replay", json!({"script": s}), &mut c);
        assert_eq!(r.exit_code, 1);
        assert_eq!(r.result["different"], 1);
        assert_eq!(r.result["calls"][1]["status"], "not-run");
        let r = crate::call(
            "journal.replay",
            json!({"script": s, "keep_going": true}),
            &mut c,
        );
        assert_eq!(r.result["replayed"], 2);
        assert_eq!(
            r.result["calls"][1]["same"], true,
            "no record: expected to finish"
        );
    }

    #[test]
    fn a_script_without_calls_is_refused() {
        let d = dir("empty");
        let s = d.join("s.jsonl");
        std::fs::write(&s, "{\"type\": \"header\"}\n").unwrap();
        let mut c = Ctx::new(BuildInfo::default());
        let r = crate::call("journal.replay", json!({"script": s}), &mut c);
        assert_eq!(r.diagnostics[0].code, diag::TOOL_SCRIPT);
        std::fs::write(&s, "not json\n").unwrap();
        let r = crate::call("journal.replay", json!({"script": s}), &mut c);
        assert_eq!(r.diagnostics[0].code, diag::TOOL_SCRIPT);
    }
}
