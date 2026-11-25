//! `case.*`: the case-file schema, the check and the resolved case.

use crate::{Action, Category, Ctx, NoInput, Tool};
use anyhow::Result;
use foehncraft_flow::case::FlowCase;
use foehncraft_flow::diag;
use foehncraft_flow::engine::{self, CheckReport, Estimate};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// A flow case file with its overrides: the input of `case.check`,
/// `case.resolve` and `flow.run`, the same options as `foehn flow`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FlowInput {
    /// Flow case file (TOML).
    pub case: PathBuf,
    /// Case keys to override, `key=value` with a TOML value, in order (D62).
    #[serde(default)]
    pub set: Vec<String>,
    /// Output directory (overrides `output.directory`).
    #[serde(default)]
    pub out: Option<PathBuf>,
    /// CPU threads (overrides `run.threads`).
    #[serde(default)]
    pub threads: Option<usize>,
    /// Refinement levels of `[mesh.generate]` (overrides `verification.levels`).
    #[serde(default)]
    pub levels: Option<Vec<u32>>,
    /// Largest number of outer iterations (overrides `solver.max_iterations`).
    #[serde(default)]
    pub max_iterations: Option<usize>,
    /// Do not write the VTKHDF file.
    #[serde(default)]
    pub no_vtkhdf: bool,
    /// Restart file, or `latest` (overrides `run.restart`).
    #[serde(default)]
    pub restart: Option<PathBuf>,
    /// End time of a transient case in s (replaces `time.end_time` and
    /// `time.steps`).
    #[serde(default)]
    pub end_time: Option<f64>,
    /// Where the flow loop runs: cpu or cuda (overrides `run.backend`).
    #[serde(default)]
    pub backend: Option<String>,
    /// PTX module of the cuda backend (overrides `run.ptx`).
    #[serde(default)]
    pub ptx: Option<PathBuf>,
    /// Solve on this number of ranks (overrides `run.rank_count`, MPI step
    /// W7): threads of this process, or the processes of an MPI launcher.
    #[serde(default)]
    pub ranks: Option<usize>,
}

impl FlowInput {
    /// The input for a case file without overrides.
    pub fn new(case: impl Into<PathBuf>) -> Self {
        Self {
            case: case.into(),
            ..Default::default()
        }
    }

    /// The options that override case keys, applied to `case`, then the
    /// validation.
    pub fn apply(&self, case: &mut FlowCase) -> Result<()> {
        if self.threads.is_some() {
            case.run.threads = self.threads;
        }
        if let Some(l) = &self.levels {
            case.verification.levels = l.clone();
        }
        if let Some(n) = self.max_iterations {
            case.solver.max_iterations = n;
        }
        if self.no_vtkhdf {
            case.output.vtkhdf = false;
        }
        if let Some(b) = &self.backend {
            case.run.backend = b.clone();
        }
        if let Some(p) = &self.ptx {
            case.run.ptx = p.clone();
        }
        if let Some(r) = &self.restart {
            case.run.restart = Some(r.clone());
        }
        if let Some(n) = self.ranks {
            case.run.rank_count = n;
        }
        if let Some(e) = self.end_time {
            let Some(t) = case.time.as_mut() else {
                return Err(diag::err(diag::Diagnostic::new(
                    diag::CASE_INCONSISTENT,
                    "--end-time needs a transient case (a [time] table)",
                )));
            };
            t.end_time = Some(e);
            t.steps = None;
            t.validate()?;
        }
        case.validate()
    }

    /// Reads the case inside the roots with the overrides applied.
    pub fn load(&self, ctx: &Ctx<'_>) -> Result<FlowCase> {
        let mut c = engine::load(&self.case, &self.set, &ctx.roots)?;
        self.apply(&mut c)?;
        Ok(c)
    }
}

/// `case.schema`.
pub struct Schema;

impl Tool for Schema {
    const ID: &'static str = "case.schema";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str =
        "The JSON Schema of the flow case file (D77) in its format envelope.";
    const CATEGORY: Category = Category::Case;
    const ACTION: Action = Action::Read;
    type Input = NoInput;
    type Output = serde_json::Value;

    fn call(_: NoInput, _: &mut Ctx<'_>) -> Result<serde_json::Value> {
        use foehncraft_flow::schema;
        Ok(schema::envelope(schema::json_schema()))
    }
}

/// `case.check`.
pub struct Check;

impl Tool for Check {
    const ID: &'static str = "case.check";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Check a flow case: validation, boundary conditions on the first mesh and the cost estimate (D78, D81).";
    const CATEGORY: Category = Category::Case;
    const ACTION: Action = Action::Read;
    type Input = FlowInput;
    type Output = CheckReport;

    fn call(a: FlowInput, ctx: &mut Ctx<'_>) -> Result<CheckReport> {
        let mut rep = CheckReport {
            format: "foehncraft.check",
            format_version: engine::FORMAT_VERSION,
            case_file: a.case.display().to_string(),
            ok: false,
            diagnostics: Vec::new(),
            estimate: None,
            operators: None,
        };
        match a.load(ctx) {
            Ok(case) => engine::check_case(
                &case,
                &a.case,
                a.out.as_deref(),
                &ctx.roots,
                &ctx.build,
                &mut rep,
            ),
            Err(e) => rep.diagnostics.push(diag::of_error(&e)),
        }
        ctx.line(&format!("{}:", a.case.display()));
        for d in &rep.diagnostics {
            ctx.line(&format!(
                "{} {}: {}{}",
                d.severity.name(),
                d.code,
                d.boundary
                    .as_ref()
                    .map_or(String::new(), |b| format!("boundary `{b}`: ")),
                d.message
            ));
            if let Some(f) = &d.fix {
                ctx.line(&format!("  fix: --set {}", f.set.join(" --set ")));
            }
        }
        if let Some(e) = &rep.estimate {
            ctx.line(&estimate_text(e));
        }
        ctx.line(if rep.ok {
            "the case can be run"
        } else {
            "the case cannot be run"
        });
        Ok(rep)
    }

    fn exit_code(rep: &CheckReport) -> i32 {
        if rep.ok { 0 } else { 1 }
    }
}

/// The estimate of a check as text.
pub fn estimate_text(e: &Estimate) -> String {
    let mb = |b: f64| {
        if b >= 1e9 {
            format!("{:.1} GB", b / 1e9)
        } else {
            format!("{:.1} MB", b / 1e6)
        }
    };
    let mut s = format!(
        "estimate: {} cells{}, {} rank(s), {} GPU(s); {} {}; host memory about {}{}; output about {}",
        e.cells,
        if e.level_cells.len() > 1 {
            format!(" (levels {:?})", e.level_cells)
        } else {
            String::new()
        },
        e.ranks,
        e.gpus,
        if e.steps_approximate { "up to" } else { "" },
        e.steps
            .iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(" + ")
            + " iterations or steps",
        mb(e.host_memory),
        if e.device_memory > 0.0 {
            format!(", device memory about {}", mb(e.device_memory))
        } else {
            String::new()
        },
        mb(e.output_bytes),
    );
    for o in &e.outputs {
        s.push_str(&format!(
            "\n  {}{}: {}",
            o.name,
            if o.count > 1 {
                format!(" ×{}", o.count)
            } else {
                String::new()
            },
            mb(o.bytes)
        ));
    }
    match e.wall_time {
        Some(t) => s.push_str(&format!(
            "\n  wall time about {t:.0} s ({})",
            e.wall_time_basis
        )),
        None => s.push_str(&format!("\n  wall time: {}", e.wall_time_basis)),
    }
    if let Some(c) = e.cell_hours {
        s.push_str(&format!(
            "\n  cell-hours about {c:.3}{}",
            e.gpu_hours
                .filter(|g| *g > 0.0)
                .map_or(String::new(), |g| format!(", GPU-hours about {g:.3}"))
        ));
    }
    if let Some(b) = &e.iteration_budget {
        s.push_str(&format!(
            "\n  iteration budget {} ({}); {}",
            b.iterations,
            b.source,
            b.diagnose_at
                .map_or("no diagnosis".to_string(), |n| format!(
                    "diagnosis at iteration {n} ({})",
                    b.diagnosis
                ))
        ));
    }
    s
}

/// The resolved case.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct ResolvedOutput {
    /// `foehncraft.resolved`.
    pub format: &'static str,
    /// Format version of the engine.
    pub format_version: &'static str,
    /// Case file.
    pub case_file: String,
    /// The case with the defaults filled in and the overrides applied, as
    /// TOML (the `case.resolved.toml` of a run).
    pub toml: String,
}

/// `case.resolve`.
pub struct Resolve;

impl Tool for Resolve {
    const ID: &'static str = "case.resolve";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Read a flow case the way a run does and return it with the defaults filled in and the overrides applied, as TOML.";
    const CATEGORY: Category = Category::Case;
    const ACTION: Action = Action::Read;
    type Input = FlowInput;
    type Output = ResolvedOutput;

    fn call(a: FlowInput, ctx: &mut Ctx<'_>) -> Result<ResolvedOutput> {
        let toml = a.load(ctx)?.resolved_toml()?;
        ctx.line(toml.trim_end_matches('\n'));
        Ok(ResolvedOutput {
            format: "foehncraft.resolved",
            format_version: engine::FORMAT_VERSION,
            case_file: a.case.display().to_string(),
            toml,
        })
    }
}

/// `case.estimate`.
pub struct EstimateTool;

impl Tool for EstimateTool {
    const ID: &'static str = "case.estimate";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "What a run of a flow case will cost before it starts (D81): cells, steps, memory, output volume, wall time.";
    const CATEGORY: Category = Category::Case;
    const ACTION: Action = Action::Read;
    type Input = FlowInput;
    type Output = Estimate;

    fn call(a: FlowInput, ctx: &mut Ctx<'_>) -> Result<Estimate> {
        let case = a.load(ctx)?;
        let (_, dir) = foehncraft_flow::case::output_dir(&case, &a.case, a.out.clone());
        let e = engine::estimate(&case, &dir, &ctx.build)?;
        ctx.line(&estimate_text(&e));
        Ok(e)
    }
}

/// Input of `case.upgrade`, the options of `foehn case upgrade`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct UpgradeInput {
    /// Case files (TOML) of any kind: flow, mesh, scalar or Poisson.
    pub files: Vec<PathBuf>,
    /// Write the files (default: list the changes only).
    #[serde(default)]
    pub write: bool,
    /// Exit code 1 when a file is not in the current format.
    #[serde(default)]
    pub check: bool,
}

/// The upgrade of case files.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct UpgradeOutput {
    /// `foehncraft.case_upgrade`.
    pub format: &'static str,
    /// Format version of the engine.
    pub format_version: &'static str,
    /// The current case format.
    pub current: &'static str,
    /// One row per file: `file`, and `changed`, `written`, `notes`, or
    /// `error`.
    pub files: Vec<serde_json::Value>,
    /// Files that need an upgrade.
    #[serde(skip)]
    pub pending: usize,
    /// Files that cannot be upgraded.
    #[serde(skip)]
    pub failed: usize,
    #[serde(skip)]
    check: bool,
}

/// `case.upgrade`.
pub struct Upgrade;

impl Tool for Upgrade {
    const ID: &'static str = "case.upgrade";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Bring case files of an older case format to the current one (renamed keys, `case_format`), keeping their comments; writes them with `write`.";
    const CATEGORY: Category = Category::Case;
    const ACTION: Action = Action::WriteFiles;
    type Input = UpgradeInput;
    type Output = UpgradeOutput;

    /// The files, when `write` rewrites them.
    fn effects(a: &UpgradeInput, _: &mut Ctx<'_>) -> Result<crate::Effects> {
        Ok(crate::Effects {
            writes: if a.write { a.files.clone() } else { Vec::new() },
            ..Default::default()
        })
    }

    fn call(a: UpgradeInput, ctx: &mut Ctx<'_>) -> Result<UpgradeOutput> {
        use anyhow::Context as _;
        use foehncraft_mesh::case_format as cf;
        use serde_json::json;
        let mut pending = 0usize;
        let mut failed = 0usize;
        let mut rows = Vec::new();
        for f in &a.files {
            ctx.roots.check_read(f, "the case file")?;
            let text =
                std::fs::read_to_string(f).with_context(|| format!("reading {}", f.display()))?;
            match cf::upgrade_text(&text) {
                Err(e) => {
                    failed += 1;
                    ctx.line(&format!("{}: {e}", f.display()));
                    rows.push(json!({ "file": f.display().to_string(), "error": e }));
                }
                Ok((new, notes)) => {
                    let changed = new != text;
                    if changed {
                        pending += 1;
                    }
                    if !changed {
                        ctx.line(&format!(
                            "{}: format {}, nothing to do",
                            f.display(),
                            cf::CURRENT
                        ));
                    } else {
                        ctx.line(&format!("{}:", f.display()));
                        for n in &notes {
                            ctx.line(&format!("  {n}"));
                        }
                        if notes.is_empty() {
                            ctx.line(&format!(
                                "  the line case_format = \"{}\" is added",
                                cf::CURRENT
                            ));
                        }
                        if a.write {
                            if let Some(d) = f.parent().filter(|d| !d.as_os_str().is_empty()) {
                                ctx.roots.check_write(d)?;
                            }
                            std::fs::write(f, &new)
                                .with_context(|| format!("writing {}", f.display()))?;
                            ctx.line("  written");
                        }
                    }
                    rows.push(json!({
                        "file": f.display().to_string(),
                        "changed": changed,
                        "written": changed && a.write,
                        "notes": notes,
                    }));
                }
            }
        }
        if pending > 0 && !a.write {
            ctx.line("run again with --write to rewrite the files");
        }
        Ok(UpgradeOutput {
            format: "foehncraft.case_upgrade",
            format_version: engine::FORMAT_VERSION,
            current: cf::CURRENT,
            files: rows,
            pending,
            failed,
            check: a.check,
        })
    }

    fn exit_code(o: &UpgradeOutput) -> i32 {
        i32::from(o.failed > 0 || (o.check && o.pending > 0))
    }
}
