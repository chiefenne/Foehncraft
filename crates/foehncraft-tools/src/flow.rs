//! `flow.run`, `result.read` and `result.page`: runs and their results.

use crate::case::FlowInput;
use crate::{Action, Category, Ctx, Effects, Tool};
use anyhow::Result;
use foehncraft_flow::case::output_dir;
use foehncraft_flow::engine::{self, Outcome, ResultSummary, RunOptions};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Runs a flow case as the call `flow.run` (journaled, D197): the case
/// read inside the roots with the overrides, the interrupts, observer and
/// build of `ctx`, the log lines to its text. The outcome keeps the report
/// for the front ends that print its failed checks.
pub fn run_outcome(a: FlowInput, ctx: &mut Ctx<'_>) -> Result<Outcome> {
    crate::run_as::<Run, _>(a, ctx, |a, c| outcome(&a, c), |o| o.summary.exit_code)
}

fn outcome(a: &FlowInput, ctx: &mut Ctx<'_>) -> Result<Outcome> {
    let case = a.load(ctx)?;
    let (name, dir) = output_dir(&case, &a.case, a.out.clone());
    let opts = RunOptions {
        interrupts: ctx.interrupts.clone(),
        build: ctx.build.clone(),
        roots: ctx.roots.clone(),
        observer: ctx.observer.clone(),
        call: ctx.call.clone(),
        // only rank 0 of the MPI processes opens the control port
        control: ctx.control && crate::mpi::rank() == 0,
        comm: crate::mpi::world(),
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
        Ok(outcome(&a, ctx)?.summary)
    }

    /// The output directory, the estimate (D81) and `[run.limits]`.
    fn effects(a: &FlowInput, ctx: &mut Ctx<'_>) -> Result<Effects> {
        let case = a.load(ctx)?;
        let (_, dir) = output_dir(&case, &a.case, a.out.clone());
        let estimate = engine::estimate(&case, &dir, &ctx.build)?;
        Ok(Effects {
            writes: vec![dir],
            estimate: Some(estimate),
            limits: Some(case.run.limits),
            publishes: false,
        })
    }

    fn exit_code(s: &ResultSummary) -> i32 {
        s.exit_code
    }
}

/// Input of `result.read`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
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

/// Input of `result.page`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PageInput {
    /// Output directory of a finished, stopped or running run.
    pub dir: PathBuf,
    /// The series: the name of a CSV file of the run without `.csv`
    /// (`residuals`, `history`, `probes`, `statistics`, `forces`, `heat`,
    /// `balance`, …) or `record` for the rows of the run record
    /// (`run.jsonl`). Without it, the call lists the series of the run.
    #[serde(default)]
    pub series: Option<String>,
    /// Where to go on: the `next_cursor` of the page before (none: from the
    /// first row).
    #[serde(default)]
    pub cursor: Option<String>,
    /// Rows per page (default 200, at most 2000).
    #[serde(default)]
    pub limit: Option<usize>,
    /// For `record`: only rows of these types (`iteration`, `table_row`,
    /// `section`, `budget`, …).
    #[serde(default)]
    pub types: Vec<String>,
}

/// One series of a run.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct SeriesInfo {
    /// Name, as `series` takes it.
    pub name: String,
    /// File in the run directory.
    pub file: String,
    /// Rows (without the header).
    pub rows: usize,
    /// Column names (empty for `record`).
    pub columns: Vec<String>,
}

/// The output of `result.page`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct PageOutput {
    /// `foehn.result-page`.
    pub format: &'static str,
    /// The run directory.
    pub dir: String,
    /// The series of the run (when the call named none).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub available: Vec<SeriesInfo>,
    /// The series of this page.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub series: Option<String>,
    /// Column names (CSV series).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub columns: Vec<String>,
    /// The rows of this page: numbers where a cell is a number (CSV), the
    /// JSON rows (`record`).
    pub rows: Vec<serde_json::Value>,
    /// Number of the first row of the page (from 0).
    pub first: usize,
    /// Rows of the series in all (after the `types` filter).
    pub total: usize,
    /// The cursor of the next page; none on the last page.
    pub next_cursor: Option<String>,
}

fn series_error(msg: String) -> anyhow::Error {
    foehncraft_flow::diag::err(foehncraft_flow::diag::Diagnostic::new(
        foehncraft_flow::diag::TOOL_SERIES,
        msg,
    ))
}

/// Lines of a CSV file without comments and blank lines.
fn csv_lines(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .map(str::trim_end)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
}

fn cell(c: &str) -> serde_json::Value {
    let c = c.trim();
    match c.parse::<f64>() {
        Ok(x) if x.is_finite() => serde_json::json!(x),
        _ => serde_json::json!(c),
    }
}

/// The series of a run directory, sorted by name.
fn series_of(dir: &std::path::Path) -> Vec<SeriesInfo> {
    let mut out = Vec::new();
    let Ok(rd) = std::fs::read_dir(dir) else {
        return out;
    };
    let mut files: Vec<PathBuf> = rd.flatten().map(|e| e.path()).collect();
    files.sort();
    for f in files {
        let Some(name) = f.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if let Some(stem) = name.strip_suffix(".csv") {
            let text = std::fs::read_to_string(&f).unwrap_or_default();
            let mut lines = csv_lines(&text);
            let columns = lines
                .next()
                .map(|h| h.split(',').map(|c| c.trim().to_string()).collect())
                .unwrap_or_default();
            out.push(SeriesInfo {
                name: stem.into(),
                file: name.into(),
                rows: lines.count(),
                columns,
            });
        } else if name == "run.jsonl" {
            let text = std::fs::read_to_string(&f).unwrap_or_default();
            out.push(SeriesInfo {
                name: "record".into(),
                file: name.into(),
                rows: text.lines().filter(|l| !l.trim().is_empty()).count(),
                columns: Vec::new(),
            });
        }
    }
    out
}

/// `result.page`.
pub struct Page;

impl Tool for Page {
    const ID: &'static str = "result.page";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "One page of a series of a run (residual history, monitors, probes, statistics, the rows of the run record), paged by cursor; without a series, the series of the run.";
    const CATEGORY: Category = Category::Result;
    const ACTION: Action = Action::Read;
    type Input = PageInput;
    type Output = PageOutput;

    fn call(a: PageInput, ctx: &mut Ctx<'_>) -> Result<PageOutput> {
        let limit = a.limit.unwrap_or(200).clamp(1, 2000);
        let mut out = PageOutput {
            format: "foehn.result-page",
            dir: a.dir.display().to_string(),
            available: Vec::new(),
            series: a.series.clone(),
            columns: Vec::new(),
            rows: Vec::new(),
            first: 0,
            total: 0,
            next_cursor: None,
        };
        ctx.roots
            .check_read(&a.dir.join("run.jsonl"), "the run directory")?;
        let available = series_of(&a.dir);
        let Some(name) = a.series.as_deref() else {
            for s in &available {
                ctx.line(&format!("{:<14} {:>7} rows  {}", s.name, s.rows, s.file));
            }
            out.total = available.len();
            out.available = available;
            return Ok(out);
        };
        let Some(info) = available.iter().find(|s| s.name == name) else {
            let names: Vec<&str> = available.iter().map(|s| s.name.as_str()).collect();
            return Err(series_error(format!(
                "{} has no series `{name}` (it has: {})",
                a.dir.display(),
                names.join(", ")
            )));
        };
        let first = match a.cursor.as_deref() {
            None => 0,
            Some(c) => c
                .strip_prefix(&format!("{name}:"))
                .and_then(|n| n.parse::<usize>().ok())
                .ok_or_else(|| {
                    series_error(format!(
                        "the cursor `{c}` does not belong to the series `{name}`"
                    ))
                })?,
        };
        let text = std::fs::read_to_string(a.dir.join(&info.file))?;
        let rows: Vec<serde_json::Value> = if name == "record" {
            text.lines()
                .filter(|l| !l.trim().is_empty())
                .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
                .filter(|v| {
                    a.types.is_empty()
                        || v["type"]
                            .as_str()
                            .is_some_and(|t| a.types.iter().any(|x| x == t))
                })
                .collect()
        } else {
            out.columns = info.columns.clone();
            csv_lines(&text)
                .skip(1)
                .map(|l| serde_json::Value::Array(l.split(',').map(cell).collect()))
                .collect()
        };
        out.total = rows.len();
        out.first = first.min(out.total);
        out.rows = rows.into_iter().skip(out.first).take(limit).collect();
        let end = out.first + out.rows.len();
        if end < out.total {
            out.next_cursor = Some(format!("{name}:{end}"));
        }
        ctx.line(&format!(
            "{name}: rows {} to {} of {}",
            out.first, end, out.total
        ));
        Ok(out)
    }
}
