//! `run.status`, `run.stop` and `run.set`: the live control of a running
//! case through its control port (L19, D198 (f)).
//!
//! The run must have been started with the control port (`foehn mcp` does it
//! for every `flow.run`; `foehn flow --control-port` for a person). The tools
//! find the port through `control.json` and the token through
//! `control.token` in the run's output directory ([`foehncraft_flow::live`]).
//! They are *live* calls: they act on a process that is going now, so a
//! journal replay skips them ([`LIVE_TOOLS`]); what they did is in the run
//! record of the run (`event` rows `RUN022` to `RUN024`, with the equivalent
//! `[[schedule]]` entry).

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::diag::{self, Diagnostic};
use foehncraft_flow::live::{self, CONTROL_FILE};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// The tools that act on a run that is going; a replay skips them.
pub const LIVE_TOOLS: &[&str] = &["run.status", "run.stop", "run.set"];

/// Input of `run.status` and `run.stop`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DirInput {
    /// The output directory of the running case (the one `flow.run` writes).
    pub dir: PathBuf,
}

/// Input of `run.set`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SetInput {
    /// The output directory of the running case.
    pub dir: PathBuf,
    /// Run-time settings to change, `{"solver.momentum.relaxation": 0.3}`:
    /// the keys that a re-read of the case file may change (D40), never
    /// `run.limits.*`. Applied at the end of the current iteration or step;
    /// the status shows what the run did with each one.
    pub set: BTreeMap<String, Value>,
}

/// What `run.stop` and `run.set` answer.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct Ack {
    /// `foehn.run-control`.
    pub format: &'static str,
    /// The run directory.
    pub dir: String,
    /// `stop` or `set`.
    pub op: &'static str,
    /// Settings queued (`run.set`).
    pub queued: Vec<String>,
    /// When it takes effect.
    pub applied_at: &'static str,
}

fn control_error(dir: &Path, e: &anyhow::Error) -> anyhow::Error {
    diag::err(Diagnostic::new(
        diag::TOOL_CONTROL,
        format!("{}: {e:#}", dir.display()),
    ))
}

fn ask(dir: &Path, req: Value, write: bool, ctx: &mut Ctx<'_>) -> Result<Value> {
    let file = dir.join(CONTROL_FILE);
    if write {
        ctx.roots.check_write(&file)?;
    } else {
        ctx.roots.check_read(&file, "the control file")?;
    }
    let answer = live::request(dir, req).map_err(|e| control_error(dir, &e))?;
    if answer["ok"] != true {
        return Err(diag::err(Diagnostic::new(
            diag::TOOL_CONTROL,
            format!(
                "{}: {}",
                dir.display(),
                answer["error"]
                    .as_str()
                    .unwrap_or("the run refused the request")
            ),
        )));
    }
    Ok(answer)
}

/// `run.status`.
pub struct Status;

impl Tool for Status {
    const ID: &'static str = "run.status";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "The live state of a running case through its control port: the last iteration or step with its residuals, the wall time, a pending stop and what the run did with earlier `run.set` calls.";
    const CATEGORY: Category = Category::Run;
    const ACTION: Action = Action::Read;
    type Input = DirInput;
    type Output = Value;

    fn call(a: DirInput, ctx: &mut Ctx<'_>) -> Result<Value> {
        let mut s = ask(&a.dir, json!({"op": "status"}), false, ctx)?;
        if let Some(o) = s.as_object_mut() {
            o.remove("ok");
        }
        let it = &s["last"];
        let at = it["iteration"]
            .as_u64()
            .map(|n| format!("iteration {n}"))
            .or_else(|| it["step"].as_u64().map(|n| format!("step {n}")))
            .unwrap_or_else(|| "no iteration yet".into());
        ctx.line(&format!(
            "{}: {} after {:.1} s, {at}{}",
            a.dir.display(),
            s["state"].as_str().unwrap_or("?"),
            s["wall_time"].as_f64().unwrap_or(0.0),
            if s["stop_requested"] == true {
                ", stop requested"
            } else {
                ""
            }
        ));
        Ok(s)
    }
}

/// `run.stop`.
pub struct Stop;

impl Tool for Stop {
    const ID: &'static str = "run.stop";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Stop a running case through its control port at the end of the current iteration or step: restart file, final output and the result summary as for a stop by signal (exit code 3, RUN022).";
    const CATEGORY: Category = Category::Run;
    const ACTION: Action = Action::Control;
    type Input = DirInput;
    type Output = Ack;

    fn call(a: DirInput, ctx: &mut Ctx<'_>) -> Result<Ack> {
        ask(&a.dir, json!({"op": "stop"}), true, ctx)?;
        ctx.line(&format!("{}: stop requested", a.dir.display()));
        Ok(Ack {
            format: "foehn.run-control",
            dir: a.dir.display().to_string(),
            op: "stop",
            queued: Vec::new(),
            applied_at: "the end of the current iteration or step",
        })
    }
}

/// `run.set`.
pub struct Set;

impl Tool for Set {
    const ID: &'static str = "run.set";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Change run-time settings of a running case through its control port (the keys a re-read of the case file may change, D40; never `run.limits.*`); applied at the end of the current iteration or step, reported in the run record with the equivalent `[[schedule]]` entry.";
    const CATEGORY: Category = Category::Run;
    const ACTION: Action = Action::Control;
    type Input = SetInput;
    type Output = Ack;

    fn call(a: SetInput, ctx: &mut Ctx<'_>) -> Result<Ack> {
        let r = ask(&a.dir, json!({"op": "set", "set": a.set}), true, ctx)?;
        let queued: Vec<String> = r["queued"]
            .as_array()
            .map(|q| {
                q.iter()
                    .filter_map(|k| k.as_str().map(String::from))
                    .collect()
            })
            .unwrap_or_default();
        ctx.line(&format!(
            "{}: queued {}",
            a.dir.display(),
            queued.join(", ")
        ));
        Ok(Ack {
            format: "foehn.run-control",
            dir: a.dir.display().to_string(),
            op: "set",
            queued,
            applied_at: "the end of the current iteration or step",
        })
    }
}
