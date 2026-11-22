//! `study.*`: the verification and solver studies of the stages before the
//! flow solver: scalar transport (stage 2), the pressure solver (stage 3)
//! and the kernel gate G1 (stage 0).

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::diag;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

fn code_of(failures: &[String]) -> i32 {
    if failures.is_empty() { 0 } else { 2 }
}

/// Input of `study.scalar`, the options of `foehn scalar`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ScalarInput {
    /// Scalar case file (TOML).
    pub case: PathBuf,
    /// Backends: `cpu`, `cpu-serial`, `cuda` (overrides `run.backends`).
    #[serde(default)]
    pub backends: Option<Vec<String>>,
    /// CPU threads (overrides `run.threads`).
    #[serde(default)]
    pub threads: Option<usize>,
    /// PTX module of the GPU backend (overrides `run.rust_ptx`).
    #[serde(default)]
    pub rust_ptx: Option<PathBuf>,
    /// Refinement levels of `[mesh.generate]` (overrides `verification.levels`).
    #[serde(default)]
    pub levels: Option<Vec<u32>>,
    /// Output directory (overrides `output.directory`).
    #[serde(default)]
    pub out: Option<PathBuf>,
}

/// Output of `study.scalar`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct ScalarOutput {
    /// Observed orders of accuracy of a sweep.
    pub orders: Vec<f64>,
    /// Failed checks.
    pub failures: Vec<String>,
    /// The report (`scalar.md`).
    pub report: PathBuf,
}

/// `study.scalar`.
pub struct Scalar;

impl Tool for Scalar {
    const ID: &'static str = "study.scalar";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Solve a steady scalar convection-diffusion case (stage 2): one mesh, an order-of-accuracy sweep or solver variants, on the CPU and the GPU.";
    const CATEGORY: Category = Category::Study;
    const ACTION: Action = Action::StartRun;
    type Input = ScalarInput;
    type Output = ScalarOutput;

    /// The output directory.
    fn effects(a: &ScalarInput, ctx: &mut Ctx<'_>) -> Result<crate::Effects> {
        use foehncraft_fv::case::{ScalarCase, output_dir};
        ctx.roots.check_read(&a.case, "the case file")?;
        let case = ScalarCase::from_file(&a.case).map_err(diag::case_error)?;
        let (_, dir) = output_dir(&case, &a.case, a.out.clone());
        Ok(crate::Effects {
            writes: vec![dir],
            ..Default::default()
        })
    }

    fn call(a: ScalarInput, ctx: &mut Ctx<'_>) -> Result<ScalarOutput> {
        use foehncraft_fv::case::{ScalarCase, output_dir, run};
        ctx.roots.check_read(&a.case, "the case file")?;
        let mut case = ScalarCase::from_file(&a.case).map_err(diag::case_error)?;
        if let Some(b) = a.backends {
            case.run.backends = b;
        }
        if a.threads.is_some() {
            case.run.threads = a.threads;
        }
        if let Some(p) = a.rust_ptx {
            case.run.rust_ptx = p;
        }
        if let Some(l) = a.levels {
            case.verification.levels = l;
        }
        let (name, dir) = output_dir(&case, &a.case, a.out);
        ctx.roots.check_write(&dir)?;
        let rep = run(&case, &name, &dir, |l| ctx.line(l))?;
        if !rep.orders.is_empty() {
            ctx.line(&format!(
                "observed order: {}",
                rep.orders
                    .iter()
                    .map(|o| format!("{o:.3}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
        ctx.line(&format!("wrote {}", dir.join("scalar.md").display()));
        Ok(ScalarOutput {
            orders: rep.orders,
            failures: rep.failures,
            report: dir.join("scalar.md"),
        })
    }

    fn exit_code(o: &ScalarOutput) -> i32 {
        code_of(&o.failures)
    }
}

/// Input of `study.poisson`, the options of `foehn poisson`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PoissonInput {
    /// Poisson case file (TOML).
    pub case: PathBuf,
    /// Backends: `cpu`, `cpu-serial`, `cuda` (overrides `study.backends`).
    #[serde(default)]
    pub backends: Option<Vec<String>>,
    /// Solvers: `amg`, `amg-double`, `amg-plain`, `jacobi`, `hypre`
    /// (overrides `study.solvers`).
    #[serde(default)]
    pub solvers: Option<Vec<String>>,
    /// CPU threads (overrides `study.threads`).
    #[serde(default)]
    pub threads: Option<usize>,
    /// PTX module of the GPU backend (overrides `study.rust_ptx`).
    #[serde(default)]
    pub rust_ptx: Option<PathBuf>,
    /// Refinement levels of `[mesh.generate]` (overrides `study.levels`).
    #[serde(default)]
    pub levels: Option<Vec<u32>>,
    /// Repetitions of each solve (overrides `study.repeat`).
    #[serde(default)]
    pub repeat: Option<usize>,
    /// Output directory (overrides `output.directory`).
    #[serde(default)]
    pub out: Option<PathBuf>,
}

/// Output of `study.poisson`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct PoissonOutput {
    /// Iterations on the finest over the coarsest level, per solver.
    pub iteration_ratios: Vec<(String, f64)>,
    /// Failed checks.
    pub failures: Vec<String>,
    /// The report (`poisson.md`).
    pub report: PathBuf,
}

/// `study.poisson`.
pub struct Poisson;

impl Tool for Poisson {
    const ID: &'static str = "study.poisson";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Study the pressure solver on a Poisson problem (stage 3): CG with multigrid, Jacobi and hypre over refinement levels, on the CPU and the GPU.";
    const CATEGORY: Category = Category::Study;
    const ACTION: Action = Action::StartRun;
    type Input = PoissonInput;
    type Output = PoissonOutput;

    /// The output directory.
    fn effects(a: &PoissonInput, ctx: &mut Ctx<'_>) -> Result<crate::Effects> {
        use foehncraft_fv::poisson::{PoissonCase, output_dir};
        ctx.roots.check_read(&a.case, "the case file")?;
        let case = PoissonCase::from_file(&a.case).map_err(diag::case_error)?;
        let (_, dir) = output_dir(&case, &a.case, a.out.clone());
        Ok(crate::Effects {
            writes: vec![dir],
            ..Default::default()
        })
    }

    fn call(a: PoissonInput, ctx: &mut Ctx<'_>) -> Result<PoissonOutput> {
        use foehncraft_fv::poisson::{PoissonCase, StudySolver, output_dir, run};
        ctx.roots.check_read(&a.case, "the case file")?;
        let mut case = PoissonCase::from_file(&a.case).map_err(diag::case_error)?;
        if let Some(b) = a.backends {
            case.study.backends = b;
        }
        if let Some(s) = a.solvers {
            case.study.solvers = s
                .iter()
                .map(|n| n.parse::<StudySolver>())
                .collect::<Result<_>>()?;
        }
        if a.threads.is_some() {
            case.study.threads = a.threads;
        }
        if let Some(p) = a.rust_ptx {
            case.study.rust_ptx = p;
        }
        if let Some(l) = a.levels {
            case.study.levels = l;
        }
        if let Some(r) = a.repeat {
            case.study.repeat = r;
        }
        case.validate()?;
        let (name, dir) = output_dir(&case, &a.case, a.out);
        ctx.roots.check_write(&dir)?;
        let rep = run(&case, &name, &dir, |l| ctx.line(l))?;
        for (n, r) in &rep.iteration_ratios {
            ctx.line(&format!("iterations finest / coarsest level, {n}: {r:.2}"));
        }
        ctx.line(&format!("wrote {}", dir.join("poisson.md").display()));
        Ok(PoissonOutput {
            iteration_ratios: rep.iteration_ratios,
            failures: rep.failures,
            report: dir.join("poisson.md"),
        })
    }

    fn exit_code(o: &PoissonOutput) -> i32 {
        code_of(&o.failures)
    }
}

/// Input of `study.g1`, the options of `foehn bench g1`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct G1Input {
    /// Case file (TOML); every key is optional.
    #[serde(default)]
    pub case: Option<PathBuf>,
    /// Approximate number of cells (overrides `mesh.cells`).
    #[serde(default)]
    pub cells: Option<usize>,
    /// Cell numbering: `natural` or `shuffled` (overrides `mesh.ordering`).
    #[serde(default)]
    pub ordering: Option<String>,
    /// Backends (overrides `run.backends`).
    #[serde(default)]
    pub backends: Option<Vec<String>>,
    /// Timed repetitions (overrides `run.repeat`).
    #[serde(default)]
    pub repeat: Option<usize>,
    /// CPU threads (overrides `run.threads`).
    #[serde(default)]
    pub threads: Option<usize>,
    /// PTX module of the single-source GPU route (overrides `run.rust_ptx`).
    #[serde(default)]
    pub rust_ptx: Option<PathBuf>,
    /// Output directory (overrides `output.directory`).
    #[serde(default)]
    pub out: Option<PathBuf>,
}

/// One verdict of G1.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct G1Verdict {
    /// The route measured.
    pub candidate: String,
    /// The route it is compared with.
    pub reference: String,
    /// It passes the gate.
    pub pass: bool,
}

/// Output of `study.g1`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct G1Output {
    /// The verdicts of the gate.
    pub verdicts: Vec<G1Verdict>,
    /// The files written.
    pub files: Vec<PathBuf>,
    /// Failed checks (a backend that did not reproduce the reference).
    #[serde(skip)]
    pub failures: Vec<String>,
}

/// `study.g1`.
pub struct G1;

impl Tool for G1 {
    const ID: &'static str = "study.g1";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Gate G1: the kernel benchmark, the single-source route against hand-written references (docs/benchmarks/g1.md).";
    const CATEGORY: Category = Category::Study;
    const ACTION: Action = Action::StartRun;
    type Input = G1Input;
    type Output = G1Output;

    /// The output directory.
    fn effects(a: &G1Input, ctx: &mut Ctx<'_>) -> Result<crate::Effects> {
        use foehncraft_bench::case::G1Case;
        let dir = match (&a.out, &a.case) {
            (Some(o), _) => o.clone(),
            (None, Some(p)) => {
                ctx.roots.check_read(p, "the case file")?;
                G1Case::from_file(p)
                    .map_err(diag::case_error)?
                    .output
                    .directory
            }
            (None, None) => G1Case::default().output.directory,
        };
        Ok(crate::Effects {
            writes: vec![dir],
            ..Default::default()
        })
    }

    fn call(a: G1Input, ctx: &mut Ctx<'_>) -> Result<G1Output> {
        use foehncraft_bench::G1Case;
        let mut case = match &a.case {
            Some(p) => {
                ctx.roots.check_read(p, "the case file")?;
                G1Case::from_file(p).map_err(diag::case_error)?
            }
            None => G1Case::default(),
        };
        if let Some(n) = a.cells {
            case.mesh.cells = n;
            case.mesh.nx = None;
        }
        if let Some(o) = &a.ordering {
            case.mesh.ordering = match o.as_str() {
                "natural" => foehncraft_bench::mesh::Ordering::Natural,
                "shuffled" => foehncraft_bench::mesh::Ordering::Shuffled,
                _ => anyhow::bail!("--ordering: natural or shuffled"),
            };
        }
        if let Some(b) = a.backends {
            case.run.backends = b;
        }
        if let Some(r) = a.repeat {
            case.run.repeat = r;
        }
        if a.threads.is_some() {
            case.run.threads = a.threads;
        }
        if let Some(p) = a.rust_ptx {
            case.run.rust_ptx = p;
        }
        if let Some(o) = a.out {
            case.output.directory = o;
        }
        ctx.roots.check_write(&case.output.directory)?;
        let report = foehncraft_bench::run_g1(&case, &mut |line| ctx.line(line))?;
        for v in &report.verdicts {
            ctx.line(&format!(
                "G1 {} against {}: {}",
                v.candidate,
                v.reference,
                if v.pass { "PASS" } else { "FAIL" }
            ));
        }
        let files = report.write(&case.output.directory, &case.output.name)?;
        for f in &files {
            ctx.line(&format!("wrote {}", f.display()));
        }
        let all_correct = report
            .results
            .iter()
            .all(|b| b.kernels.iter().all(|k| k.within_tolerance));
        Ok(G1Output {
            verdicts: report
                .verdicts
                .iter()
                .map(|v| G1Verdict {
                    candidate: v.candidate.to_string(),
                    reference: v.reference.to_string(),
                    pass: v.pass,
                })
                .collect(),
            files,
            failures: if all_correct {
                Vec::new()
            } else {
                vec!["a backend did not reproduce the reference results".into()]
            },
        })
    }

    fn exit_code(o: &G1Output) -> i32 {
        code_of(&o.failures)
    }
}
