//! `flow.run` and `result.read`: runs and their results.

use crate::case::FlowInput;
use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::case::output_dir;
use foehncraft_flow::engine::{self, Outcome, ResultSummary, RunOptions};
use schemars::JsonSchema;
use serde::Deserialize;
use std::path::PathBuf;

/// Runs a flow case: the case read inside the roots with the overrides,
/// the interrupts, observer and build of `ctx`, the log lines to its text.
/// The outcome keeps the report for the front ends that print its failed
/// checks.
pub fn run_outcome(a: &FlowInput, ctx: &mut Ctx<'_>) -> Result<Outcome> {
    let case = a.load(ctx)?;
    let (name, dir) = output_dir(&case, &a.case, a.out.clone());
    let opts = RunOptions {
        interrupts: ctx.interrupts.clone(),
        build: ctx.build.clone(),
        roots: ctx.roots.clone(),
        observer: ctx.observer.clone(),
    };
    engine::run(&case, &name, &dir, &opts, |l| ctx.line(l))
}

/// `flow.run`.
pub struct Run;

impl Tool for Run {
    const ID: &'static str = "flow.run";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Run a flow case (steady, transient, sweeps) and return the result summary of `result.json`; streams the rows of the run record and stops on request.";
    const CATEGORY: Category = Category::Flow;
    const ACTION: Action = Action::StartRun;
    const STREAMS: bool = true;
    type Input = FlowInput;
    type Output = ResultSummary;

    fn call(a: FlowInput, ctx: &mut Ctx<'_>) -> Result<ResultSummary> {
        Ok(run_outcome(&a, ctx)?.summary)
    }

    fn exit_code(s: &ResultSummary) -> i32 {
        s.exit_code
    }
}

/// Input of `result.read`.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReadInput {
    /// Output directory of a finished or stopped run.
    pub dir: PathBuf,
}

/// `result.read`.
pub struct Read;

impl Tool for Read {
    const ID: &'static str = "result.read";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str =
        "The result summary (`result.json`, `foehncraft.result`) of an earlier run.";
    const CATEGORY: Category = Category::Result;
    const ACTION: Action = Action::Read;
    type Input = ReadInput;
    type Output = serde_json::Value;

    fn call(a: ReadInput, ctx: &mut Ctx<'_>) -> Result<serde_json::Value> {
        ctx.roots
            .check_read(&a.dir.join("result.json"), "the result file")?;
        engine::read_result(&a.dir)
    }
}
