//! `foehn`: the command line of Foehncraft (D4).
//!
//! Commands: `info` (build and platform), `mesh` (build or read a mesh,
//! check, partition and write it; stage 1), `scalar` (steady scalar
//! transport, verification sweeps, CPU and GPU; stage 2), `poisson` (the
//! pressure solver study: CG with algebraic multigrid against hypre; stage
//! 3), `flow` (steady incompressible laminar flow with SIMPLEC; stage 4)
//! and `bench g1` (the GPU kernel gate of stage 0), `schema` (the JSON
//! Schema of the flow case file, D77), `post derive` (derived fields from
//! a restart file, D72), `materials` (the material database, D94),
//! `case resolve`, `diag codes`, `expression`, and `commands` and `call`
//! (the tool registry, D196: the catalogue of every tool, and any tool run
//! from JSON). The commands build the typed input of their tool of
//! `foehncraft-tools` and call it.
//!
//! Output streams (D79): human text on standard output, or with `--json`
//! on standard error and one JSON object on standard output (the result
//! summary of a flow run, the check report of `--check`, a command summary
//! for the others, an error object with a diagnostic code on failure).
//! `--quiet` silences the human text. No command asks questions. Exit
//! codes (D64): 0 finished with every check passed, 1 error, 2 finished
//! with failed checks, 3 stopped on request (signal, `.stop`, a limit).

mod mcp;

use anyhow::Result;
use clap::{Parser, Subcommand};
use foehncraft_flow::diag;
use foehncraft_flow::engine::{self, BuildInfo, Interrupts, Roots};
use foehncraft_tools::{self as tools, Ctx, NoInput};
use serde_json::json;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "foehn",
    version,
    about = "Foehncraft – general purpose CFD solver"
)]
struct Cli {
    /// Machine output: one JSON object on standard output, human text on
    /// standard error (D79).
    #[arg(long, global = true)]
    json: bool,
    /// No human text.
    #[arg(long, global = true)]
    quiet: bool,
    /// Append every tool call to this journal (D197; default: the
    /// environment variable `FOEHN_JOURNAL`); `foehn run --script`
    /// replays it.
    #[arg(long, global = true, value_name = "FILE")]
    journal: Option<PathBuf>,
    #[command(subcommand)]
    command: Command,
}

/// Where human text goes.
#[derive(Clone, Copy)]
struct Out {
    json: bool,
    quiet: bool,
}

impl Out {
    fn line(self, l: &str) {
        if self.quiet {
        } else if self.json {
            eprintln!("{l}");
        } else {
            stdout(l);
        }
    }
}

/// One line on standard output. A closed pipe (`foehn info | head -1`)
/// ends the output instead of a panic.
fn stdout(l: &str) {
    use std::io::Write;
    let _ = writeln!(std::io::stdout().lock(), "{l}");
}

/// A finished command: exit code and the JSON summary of `--json`.
struct Done {
    code: i32,
    json: serde_json::Value,
}

#[derive(Subcommand)]
enum Command {
    /// Print the version, the enabled features and the platform.
    Info,
    /// Build or read a mesh, check it, partition it and write the quality
    /// report and a VTKHDF file for ParaView.
    Mesh(MeshArgs),
    /// Solve a steady scalar convection–diffusion case (stage 2 operators):
    /// one mesh, an order-of-accuracy sweep, or several solver variants, on
    /// the CPU and the GPU.
    Scalar(ScalarArgs),
    /// Study the pressure solver on a Poisson problem (stage 3): CG with the
    /// smoothed-aggregation multigrid, Jacobi and hypre, over refinement
    /// levels, on the CPU and the GPU.
    Poisson(PoissonArgs),
    /// Solve an incompressible laminar flow case: steady with SIMPLEC or
    /// SIMPLE (stage 4), transient with PISO and BDF2 or Crank–Nicolson
    /// when the case has a `[time]` table (stage 5); on one mesh, an
    /// order-of-accuracy sweep in space or a sweep over time steps.
    /// `--check` only checks the case file and the boundary conditions.
    /// Control files in the output directory steer a running case: `.stop`,
    /// `.result_N`, `.read_toml_N`.
    Flow(FlowArgs),
    /// Benchmarks.
    Bench {
        #[command(subcommand)]
        bench: Bench,
    },
    /// Print the JSON Schema of the flow case file (D77): type, range,
    /// default, unit and description of every key, and whether it may
    /// change during a run.
    Schema(SchemaArgs),
    /// List the built-in material database (D94) with the density and
    /// viscosity of each material at a temperature, or show one material
    /// with its property forms, data range and sources.
    Materials(MaterialsArgs),
    /// Case files: bring a file of an older case format to the current
    /// one (`foehn case upgrade`).
    Case {
        #[command(subcommand)]
        case: CaseCmd,
    },
    /// Post-processing of a finished or stopped flow run.
    Post {
        #[command(subcommand)]
        post: Post,
    },
    /// Diagnostic codes (D78).
    Diag {
        #[command(subcommand)]
        diag: DiagCmd,
    },
    /// Compile an expression (D53) with the solver's own compiler and
    /// evaluate it at points.
    Expression(ExpressionArgs),
    /// Print the catalogue of every tool (D196) as JSON: ids, versions,
    /// summaries, action kinds, input and output schemas.
    Commands(CommandsArgs),
    /// Run any tool of the catalogue from JSON input and print the result
    /// envelope (`foehn.tool-result`); human text goes to standard error.
    Call(CallArgs),
    /// Replay the tool calls of a journal (`--journal`) or of a run record
    /// (`run.jsonl`) in order, without a window (D197); exit 1 when a call
    /// ends differently from its record.
    Run(RunArgs),
    /// Serve the tools to an assistant as an MCP server on standard input
    /// and output (L18). Needs `--read-root` and `--write-root`; the action
    /// table `foehn-mcp.toml` says which calls run, which ask the person and
    /// which are refused (D198); every call is journaled (default
    /// `foehn-mcp-journal.jsonl` beside the table). docs/agents.md says how
    /// to connect an assistant.
    Mcp(McpArgs),
}

#[derive(clap::Args)]
struct McpArgs {
    /// Allow reading only below this directory (repeatable; required, D198 (a)).
    #[arg(long, value_name = "DIR")]
    read_root: Vec<PathBuf>,
    /// Allow writing only below this directory (repeatable; required).
    #[arg(long, value_name = "DIR")]
    write_root: Vec<PathBuf>,
    /// The action table of the project (default: foehn-mcp.toml in the
    /// working directory; must lie outside the write roots).
    #[arg(long, value_name = "FILE")]
    policy: Option<PathBuf>,
    /// How tool ids appear as MCP tool names: `dotted` (the registry ids)
    /// or `underscore` (for hosts that refuse dots in names).
    #[arg(long, value_enum, default_value = "dotted")]
    names: mcp::Names,
}

#[derive(clap::Args)]
struct RunArgs {
    /// The journal or run record to replay.
    #[arg(long, value_name = "FILE")]
    script: PathBuf,
    /// Go on after a call that ends differently from its record.
    #[arg(long)]
    keep_going: bool,
    /// Allow reading only below this directory (repeatable, D82).
    #[arg(long, value_name = "DIR")]
    read_root: Vec<PathBuf>,
    /// Allow writing only below this directory (repeatable, D82).
    #[arg(long, value_name = "DIR")]
    write_root: Vec<PathBuf>,
}

#[derive(Subcommand)]
enum DiagCmd {
    /// List the diagnostic codes with their usual severity and meaning.
    Codes,
}

#[derive(clap::Args)]
struct ExpressionArgs {
    /// The expression, as in a case file.
    source: String,
    /// A TOML file with `[constants]` and `[tables]` the expression may
    /// refer to (relative table files from its directory).
    #[arg(long, value_name = "FILE")]
    env: Option<PathBuf>,
    /// The variables it may use: a string of x, y, z, t.
    #[arg(long, default_value = "xyzt")]
    allowed: String,
    /// A point `x,y,z,t` at which it is evaluated (repeatable).
    #[arg(long = "point", value_name = "X,Y,Z,T")]
    points: Vec<String>,
}

#[derive(clap::Args)]
struct CommandsArgs {
    /// The catalogue as Markdown (docs/tools.md) instead of JSON.
    #[arg(long)]
    markdown: bool,
}

#[derive(clap::Args)]
struct CallArgs {
    /// Tool id (`foehn commands` lists them), e.g. `case.check`.
    id: String,
    /// JSON file with the input of the tool, `-` for standard input
    /// (default: no fields).
    #[arg(long, value_name = "FILE")]
    input: Option<PathBuf>,
    /// Allow reading only below this directory (repeatable, D82).
    #[arg(long, value_name = "DIR")]
    read_root: Vec<PathBuf>,
    /// Allow writing only below this directory (repeatable, D82).
    #[arg(long, value_name = "DIR")]
    write_root: Vec<PathBuf>,
}

#[derive(Subcommand)]
enum CaseCmd {
    /// Rewrite case files in the current case format (renamed keys, the
    /// line `case_format`), keeping their comments. Without `--write` only
    /// the changes are listed; `--check` fails when a file needs one.
    Upgrade(UpgradeArgs),
    /// Print a flow case the way a run reads it: defaults filled in,
    /// overrides applied (the `case.resolved.toml` of a run).
    Resolve(ResolveArgs),
}

#[derive(clap::Args)]
struct ResolveArgs {
    /// Flow case file (TOML).
    case: PathBuf,
    /// Override a case key, `--set fluid.viscosity=1e-3` (repeatable, D62).
    #[arg(long = "set", value_name = "KEY=VALUE")]
    set: Vec<String>,
    /// Allow reading only below this directory (repeatable, D82).
    #[arg(long, value_name = "DIR")]
    read_root: Vec<PathBuf>,
}

#[derive(clap::Args)]
struct UpgradeArgs {
    /// Case files (TOML) of any kind: flow, mesh, scalar or Poisson.
    #[arg(required = true)]
    files: Vec<PathBuf>,
    /// Write the files (default: list the changes only).
    #[arg(long)]
    write: bool,
    /// Exit with code 1 when a file is not in the current format.
    #[arg(long)]
    check: bool,
}

#[derive(Subcommand)]
enum Post {
    /// Write the fields with derived fields (D72: vorticity, q-criterion,
    /// lambda2, enstrophy, strain-rate, pressure-coefficient,
    /// total-pressure, wall-shear-stress, y-plus) from a restart file of
    /// the case, the same numbers as `[output] derived` at write time.
    Derive(DeriveArgs),
}

#[derive(clap::Args)]
struct DeriveArgs {
    /// Flow case file (TOML) of the run.
    case: PathBuf,
    /// Restart file, or `latest` for the newest one of the output
    /// directory.
    #[arg(long, default_value = "latest")]
    restart: PathBuf,
    /// Derived fields, comma-separated (default: `output.derived` of the
    /// case, or every field when that is empty).
    #[arg(long, value_delimiter = ',')]
    fields: Option<Vec<String>>,
    /// File to write (default: `derived.vtkhdf` for `restart.h5`,
    /// `derived_<N>.vtkhdf` for `restart_<N>.h5`, in the directory of the
    /// restart file).
    #[arg(long)]
    file: Option<PathBuf>,
    /// Override a case key, `--set output.reference.velocity=2` (repeatable,
    /// D62; the value is TOML, a bare word is a string).
    #[arg(long = "set", value_name = "KEY=VALUE")]
    set: Vec<String>,
    /// Allow reading only below this directory (repeatable, D82).
    #[arg(long, value_name = "DIR")]
    read_root: Vec<PathBuf>,
    /// Allow writing only below this directory (repeatable, D82).
    #[arg(long, value_name = "DIR")]
    write_root: Vec<PathBuf>,
    /// Output directory of the run (overrides `output.directory`), where
    /// `latest` looks.
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(clap::Args)]
struct MaterialsArgs {
    /// One material (default: all).
    name: Option<String>,
    /// Temperature (K) at which the properties are evaluated.
    #[arg(long, default_value_t = foehncraft_flow::materials::DEFAULT_TEMPERATURE)]
    temperature: f64,
    /// Operating pressure (Pa) of the ideal-gas densities.
    #[arg(long, default_value_t = foehncraft_flow::materials::DEFAULT_PRESSURE)]
    pressure: f64,
    /// A user material database file (repeatable), searched before the
    /// built-in database like `[materials] include` of a case; the files
    /// of `FOEHN_MATERIALS` are read as well.
    #[arg(long = "include", value_name = "FILE")]
    include: Vec<PathBuf>,
}

#[derive(clap::Args)]
struct SchemaArgs {
    /// The key reference as Markdown (docs/keys.md) instead of JSON.
    #[arg(long)]
    markdown: bool,
    /// Only the JSON Schema document, without the format envelope (for
    /// editors that validate TOML against a schema).
    #[arg(long)]
    schema_only: bool,
}

#[derive(Subcommand)]
enum Bench {
    /// Gate G1: kernel benchmark, single-source route against hand-written
    /// references (docs/benchmarks/g1.md).
    G1(G1Args),
    /// List the backends compiled into this build.
    Backends,
}

#[derive(clap::Args)]
struct MeshArgs {
    /// Mesh case file (TOML) or a mesh file (Gmsh or Fluent `.msh`, SU2
    /// `.su2`, CGNS `.cgns`, VTK `.vtk`/`.vtu`, OpenFOAM case or polyMesh
    /// directory).
    input: PathBuf,
    /// Number of parts (overrides `partition.parts`).
    #[arg(long)]
    parts: Option<usize>,
    /// Partitioning method: metis or rcb (overrides `partition.method`).
    #[arg(long)]
    method: Option<String>,
    /// Cell order: none, rcm or shuffle (overrides `mesh.renumber`).
    #[arg(long)]
    renumber: Option<String>,
    /// Random points of the point-location test (overrides
    /// `check.locate_points`).
    #[arg(long)]
    locate: Option<usize>,
    /// Do not write the VTKHDF file.
    #[arg(long)]
    no_vtkhdf: bool,
    /// Also write the mesh as `mesh.vtu` (VTK XML).
    #[arg(long)]
    vtu: bool,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(clap::Args)]
struct ScalarArgs {
    /// Scalar case file (TOML).
    case: PathBuf,
    /// Comma-separated backends: cpu, cpu-serial, cuda (overrides
    /// `run.backends`).
    #[arg(long, value_delimiter = ',')]
    backends: Option<Vec<String>>,
    /// CPU threads (overrides `run.threads`).
    #[arg(long)]
    threads: Option<usize>,
    /// PTX module of the GPU backend (overrides `run.rust_ptx`).
    #[arg(long)]
    rust_ptx: Option<PathBuf>,
    /// Refinement levels of `[mesh.generate]`, comma-separated (overrides
    /// `verification.levels`).
    #[arg(long, value_delimiter = ',')]
    levels: Option<Vec<u32>>,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(clap::Args)]
struct PoissonArgs {
    /// Poisson case file (TOML).
    case: PathBuf,
    /// Comma-separated backends: cpu, cpu-serial, cuda (overrides
    /// `study.backends`).
    #[arg(long, value_delimiter = ',')]
    backends: Option<Vec<String>>,
    /// Comma-separated solvers: amg, amg-double, amg-plain, jacobi, hypre
    /// (overrides `study.solvers`).
    #[arg(long, value_delimiter = ',')]
    solvers: Option<Vec<String>>,
    /// CPU threads (overrides `study.threads`).
    #[arg(long)]
    threads: Option<usize>,
    /// PTX module of the GPU backend (overrides `study.rust_ptx`).
    #[arg(long)]
    rust_ptx: Option<PathBuf>,
    /// Refinement levels of `[mesh.generate]`, comma-separated (overrides
    /// `study.levels`).
    #[arg(long, value_delimiter = ',')]
    levels: Option<Vec<u32>>,
    /// Repetitions of each solve (overrides `study.repeat`).
    #[arg(long)]
    repeat: Option<usize>,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(clap::Args)]
struct FlowArgs {
    /// Flow case file (TOML).
    case: PathBuf,
    /// Check the case file and the boundary conditions, estimate the cost
    /// of the run (D81), then stop.
    #[arg(long)]
    check: bool,
    /// Override a case key, `--set fluid.viscosity=1e-3` (repeatable, D62;
    /// the value is TOML, a bare word is a string).
    #[arg(long = "set", value_name = "KEY=VALUE")]
    set: Vec<String>,
    /// Allow reading only below this directory (repeatable, D82).
    #[arg(long, value_name = "DIR")]
    read_root: Vec<PathBuf>,
    /// Allow writing only below this directory (repeatable, D82).
    #[arg(long, value_name = "DIR")]
    write_root: Vec<PathBuf>,
    /// Open the live control port of the run (L19): `run.status`,
    /// `run.stop` and `run.set` (`foehn call`) reach the run through
    /// `control.json` and `control.token` in its output directory.
    #[arg(long)]
    control_port: bool,
    /// CPU threads (overrides `run.threads`).
    #[arg(long)]
    threads: Option<usize>,
    /// Refinement levels of `[mesh.generate]`, comma-separated (overrides
    /// `verification.levels`).
    #[arg(long, value_delimiter = ',')]
    levels: Option<Vec<u32>>,
    /// Largest number of outer iterations (overrides
    /// `solver.max_iterations`).
    #[arg(long)]
    max_iterations: Option<usize>,
    /// Do not write the VTKHDF file.
    #[arg(long)]
    no_vtkhdf: bool,
    /// Continue from this restart file, or `latest` for the newest one of
    /// the output directory (overrides `run.restart`).
    #[arg(long)]
    restart: Option<PathBuf>,
    /// End time of a transient case in s (replaces `time.end_time` and
    /// `time.steps`).
    #[arg(long)]
    end_time: Option<f64>,
    /// Where the flow loop runs: cpu or cuda (overrides `run.backend`).
    #[arg(long)]
    backend: Option<String>,
    /// PTX module of the cuda backend (overrides `run.ptx`).
    #[arg(long)]
    ptx: Option<PathBuf>,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(clap::Args)]
struct G1Args {
    /// Case file (TOML); every key is optional.
    case: Option<PathBuf>,
    /// Approximate number of cells (overrides `mesh.cells`).
    #[arg(long)]
    cells: Option<usize>,
    /// Cell numbering: natural or shuffled (overrides `mesh.ordering`).
    #[arg(long)]
    ordering: Option<String>,
    /// Comma-separated backends (overrides `run.backends`).
    #[arg(long, value_delimiter = ',')]
    backends: Option<Vec<String>>,
    /// Timed repetitions (overrides `run.repeat`).
    #[arg(long)]
    repeat: Option<usize>,
    /// CPU threads (overrides `run.threads`).
    #[arg(long)]
    threads: Option<usize>,
    /// PTX module of the single-source GPU route (overrides `run.rust_ptx`).
    #[arg(long)]
    rust_ptx: Option<PathBuf>,
    /// Output directory (overrides `output.directory`).
    #[arg(long)]
    out: Option<PathBuf>,
}

fn main() {
    // GPU backends (feature `cuda`; nothing without it)
    foehncraft_cuda::register();
    let cli = Cli::parse();
    let _ = JOURNAL.set(cli.journal.clone());
    let out = Out {
        json: cli.json,
        quiet: cli.quiet,
    };
    let name = match &cli.command {
        Command::Info => "info",
        Command::Mesh(_) => "mesh",
        Command::Scalar(_) => "scalar",
        Command::Poisson(_) => "poisson",
        Command::Flow(_) => "flow",
        Command::Bench { .. } => "bench",
        Command::Schema(_) => "schema",
        Command::Materials(_) => "materials",
        Command::Case { .. } => "case",
        Command::Post { .. } => "post",
        Command::Diag { .. } => "diag",
        Command::Expression(_) => "expression",
        Command::Commands(_) => "commands",
        Command::Call(_) => "call",
        Command::Run(_) => "run",
        Command::Mcp(_) => "mcp",
    };
    let res = match cli.command {
        Command::Info => info(out),
        Command::Mesh(a) => mesh(a, out),
        Command::Scalar(a) => scalar(a, out),
        Command::Poisson(a) => poisson(a, out),
        Command::Flow(a) => flow(a, out),
        Command::Schema(a) => schema(a, out),
        Command::Materials(a) => materials(a, out),
        Command::Case {
            case: CaseCmd::Upgrade(a),
        } => upgrade(a, out),
        Command::Case {
            case: CaseCmd::Resolve(a),
        } => resolve(a, out),
        Command::Diag {
            diag: DiagCmd::Codes,
        } => typed::<tools::diag_codes::Codes>(NoInput {}, out),
        Command::Expression(a) => expression(a, out),
        Command::Commands(a) => commands(a),
        Command::Call(a) => call(a, out),
        Command::Run(a) => replay(a, out),
        Command::Mcp(a) => mcp_serve(a, out),
        Command::Post {
            post: Post::Derive(a),
        } => derive(a, out),
        Command::Bench {
            bench: Bench::Backends,
        } => typed::<tools::system::Backends>(NoInput {}, out),
        Command::Bench {
            bench: Bench::G1(a),
        } => g1(a, out),
    };
    let code = match res {
        Ok(d) => {
            if out.json && !d.json.is_null() {
                let v = if d.json.get("format").is_some() {
                    d.json
                } else {
                    let mut v = json!({
                        "format": "foehncraft.command",
                        "format_version": engine::FORMAT_VERSION,
                        "command": name,
                        "exit_code": d.code,
                    });
                    if let (Some(m), serde_json::Value::Object(x)) = (v.as_object_mut(), d.json) {
                        m.extend(x);
                    }
                    v
                };
                stdout(&serde_json::to_string_pretty(&v).unwrap_or_default());
            }
            d.code
        }
        Err(e) => {
            let d = diag::of_error(&e);
            if !out.quiet || !out.json {
                eprintln!("Error: {e:#}");
            }
            if out.json {
                let v = json!({
                    "format": "foehncraft.error",
                    "format_version": engine::FORMAT_VERSION,
                    "command": name,
                    "exit_code": 1,
                    "diagnostics": [d],
                });
                stdout(&serde_json::to_string_pretty(&v).unwrap_or_default());
            }
            1
        }
    };
    std::process::exit(code);
}

/// `foehn mcp`.
fn mcp_serve(a: McpArgs, out: Out) -> Result<Done> {
    let journal = JOURNAL
        .get()
        .cloned()
        .flatten()
        .or_else(|| tools::journal::Journal::from_env("mcp").map(|j| j.path().to_path_buf()));
    let code = mcp::serve(mcp::Options {
        read_roots: a.read_root,
        write_roots: a.write_root,
        policy: a.policy,
        journal,
        names: a.names,
        quiet: out.quiet,
        build: build_info(),
    })?;
    Ok(Done {
        code,
        json: serde_json::Value::Null,
    })
}

/// `foehn post derive`.
fn derive(a: DeriveArgs, out: Out) -> Result<Done> {
    let roots = Roots {
        read: a.read_root,
        write: a.write_root,
    };
    let input = tools::post::DeriveInput {
        case: a.case,
        restart: a.restart,
        fields: a.fields,
        file: a.file,
        set: a.set,
        out: a.out,
    };
    let (o, code) = tools::run::<tools::post::Derive>(input, &mut ctx(out).roots(roots))?;
    Ok(Done { code, json: o })
}

/// The features of this build.
fn features() -> Vec<String> {
    [
        ("cuda", cfg!(feature = "cuda")),
        ("mpi", cfg!(feature = "mpi")),
        ("checked-indexing", cfg!(feature = "checked-indexing")),
        ("hdf5", cfg!(feature = "hdf5")),
        ("metis", cfg!(feature = "metis")),
        ("hypre", cfg!(feature = "hypre")),
    ]
    .into_iter()
    .filter(|&(_, on)| on)
    .map(|(n, _)| n.to_string())
    .collect()
}

/// Build and host for the provenance (D85).
fn build_info() -> BuildInfo {
    BuildInfo {
        commit: env!("FOEHNCRAFT_GIT_COMMIT").into(),
        features: features(),
        platform: tools::system::platform(),
        cpu: tools::system::cpu_model(),
        gpu: None,
    }
}

/// Connects the signals of D64 to the requests of the run: SIGUSR1 writes a
/// restart file and goes on; SIGTERM, SIGINT and SIGUSR2 stop with a
/// restart file and the final output; a second SIGINT ends at once (exit
/// code 3).
fn signals() -> Interrupts {
    let i = Interrupts::new();
    #[cfg(unix)]
    {
        use engine::StopSignal;
        use signal_hook::consts::{SIGINT, SIGTERM, SIGUSR1, SIGUSR2};
        use signal_hook::flag;
        use std::sync::Arc;
        use std::sync::atomic::AtomicBool;
        let seen = Arc::new(AtomicBool::new(false));
        let mut ok = flag::register_conditional_shutdown(SIGINT, 3, seen.clone()).is_ok();
        ok &= flag::register(SIGINT, seen).is_ok();
        for (sig, s) in [
            (SIGINT, StopSignal::Int),
            (SIGTERM, StopSignal::Term),
            (SIGUSR2, StopSignal::Usr2),
        ] {
            ok &= flag::register_usize(sig, i.stop_flag(), s.value()).is_ok();
        }
        ok &= flag::register(SIGUSR1, i.restart_flag()).is_ok();
        if !ok {
            eprintln!("warning: the signal handlers could not all be installed");
        }
    }
    i
}

fn materials(a: MaterialsArgs, out: Out) -> Result<Done> {
    use tools::materials as m;
    match a.name {
        Some(name) => typed::<m::Show>(
            m::ShowInput {
                name,
                temperature: a.temperature,
                pressure: a.pressure,
                include: a.include,
            },
            out,
        ),
        None => typed::<m::List>(
            m::ListInput {
                temperature: a.temperature,
                pressure: a.pressure,
                include: a.include,
            },
            out,
        ),
    }
}

/// The context of a tool call: the build of this binary, the human text
/// to `out`.
/// The journal of `--journal` (D197).
static JOURNAL: std::sync::OnceLock<Option<PathBuf>> = std::sync::OnceLock::new();

/// The journal of this command: `--journal`, else `FOEHN_JOURNAL`.
fn journal() -> Option<tools::journal::Journal> {
    match JOURNAL.get().cloned().flatten() {
        Some(p) => Some(tools::journal::Journal::new(p, "cli")),
        None => tools::journal::Journal::from_env("cli"),
    }
}

fn ctx(out: Out) -> Ctx<'static> {
    Ctx::new(build_info())
        .journal(journal())
        .text(move |l| out.line(l))
}

/// `foehn run --script`.
fn replay(a: RunArgs, out: Out) -> Result<Done> {
    let mut c = ctx(out).roots(Roots {
        read: a.read_root,
        write: a.write_root,
    });
    c.interrupts = signals();
    let input = tools::journal::ReplayInput {
        script: a.script,
        keep_going: a.keep_going,
    };
    let (o, code) = tools::run::<tools::journal::Replay>(input, &mut c)?;
    Ok(Done {
        code,
        json: serde_json::to_value(o)?,
    })
}

/// Runs tool `T` with a typed input; its output is the `--json` object.
fn typed<T: tools::Tool>(input: T::Input, out: Out) -> Result<Done> {
    let (o, code) = tools::run::<T>(input, &mut ctx(out))?;
    Ok(Done {
        code,
        json: serde_json::to_value(o)?,
    })
}

/// `foehn case resolve`.
fn resolve(a: ResolveArgs, out: Out) -> Result<Done> {
    let input = tools::case::FlowInput {
        set: a.set,
        ..tools::case::FlowInput::new(a.case)
    };
    let roots = Roots {
        read: a.read_root,
        write: Vec::new(),
    };
    let (o, code) = tools::run::<tools::case::Resolve>(input, &mut ctx(out).roots(roots))?;
    Ok(Done {
        code,
        json: serde_json::to_value(o)?,
    })
}

/// `foehn expression`.
fn expression(a: ExpressionArgs, out: Out) -> Result<Done> {
    use anyhow::Context as _;
    let mut input = tools::expression::EvalInput::new(a.source);
    input.allowed = a.allowed;
    if let Some(f) = &a.env {
        input.env =
            std::fs::read_to_string(f).with_context(|| format!("reading {}", f.display()))?;
        input.base = f.parent().map(std::path::Path::to_path_buf);
    }
    for p in &a.points {
        let v: Vec<f64> = p
            .split(',')
            .map(|x| x.trim().parse::<f64>())
            .collect::<std::result::Result<_, _>>()
            .with_context(|| format!("--point {p}: four numbers x,y,z,t"))?;
        anyhow::ensure!(v.len() == 4, "--point {p}: four numbers x,y,z,t");
        input.points.push([v[0], v[1], v[2], v[3]]);
    }
    typed::<tools::expression::Eval>(input, out)
}

/// `foehn commands`.
fn commands(a: CommandsArgs) -> Result<Done> {
    if a.markdown {
        stdout(tools::markdown::catalogue().trim_end_matches('\n'));
    } else {
        stdout(&serde_json::to_string_pretty(&tools::catalogue())?);
    }
    Ok(Done {
        code: 0,
        json: serde_json::Value::Null,
    })
}

/// `foehn call`: the envelope on standard output, human text on standard
/// error.
fn call(a: CallArgs, out: Out) -> Result<Done> {
    use anyhow::Context as _;
    use std::io::Read as _;
    let text = match &a.input {
        None => String::new(),
        Some(p) if p.as_os_str() == "-" => {
            let mut t = String::new();
            std::io::stdin()
                .read_to_string(&mut t)
                .context("reading the input from standard input")?;
            t
        }
        Some(p) => {
            std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?
        }
    };
    let input: serde_json::Value = if text.trim().is_empty() {
        serde_json::Value::Null
    } else {
        serde_json::from_str(&text).map_err(|e| {
            diag::err(diag::Diagnostic::new(
                diag::TOOL_INPUT,
                format!("the input is not JSON: {e}"),
            ))
        })?
    };
    let quiet = out.quiet;
    let mut c = Ctx::new(build_info())
        .journal(journal())
        .roots(Roots {
            read: a.read_root,
            write: a.write_root,
        })
        .text(move |l| {
            if !quiet {
                eprintln!("{l}");
            }
        });
    c.interrupts = signals();
    let r = tools::call(&a.id, input, &mut c);
    for d in &r.diagnostics {
        if !quiet {
            eprintln!("Error {}: {}", d.code, d.message);
        }
    }
    stdout(&serde_json::to_string_pretty(&r)?);
    Ok(Done {
        code: r.exit_code,
        json: serde_json::Value::Null,
    })
}

/// `foehn case upgrade`.
fn upgrade(a: UpgradeArgs, out: Out) -> Result<Done> {
    typed::<tools::case::Upgrade>(
        tools::case::UpgradeInput {
            files: a.files,
            write: a.write,
            check: a.check,
        },
        out,
    )
}

fn schema(a: SchemaArgs, out: Out) -> Result<Done> {
    use foehncraft_flow::schema;
    if a.markdown {
        stdout(schema::markdown().trim_end_matches('\n'));
        return Ok(Done {
            code: 0,
            json: serde_json::Value::Null,
        });
    }
    if a.schema_only || !out.json {
        let (env, _) = tools::run::<tools::case::Schema>(NoInput {}, &mut ctx(out))?;
        stdout(&serde_json::to_string_pretty(&env["schema"])?);
        return Ok(Done {
            code: 0,
            json: serde_json::Value::Null,
        });
    }
    typed::<tools::case::Schema>(NoInput {}, out)
}

fn info(out: Out) -> Result<Done> {
    typed::<tools::system::Info>(NoInput {}, out)
}

fn mesh(a: MeshArgs, out: Out) -> Result<Done> {
    let input = tools::mesh::MeshInput {
        input: a.input,
        parts: a.parts,
        method: a.method,
        renumber: a.renumber,
        locate: a.locate,
        no_vtkhdf: a.no_vtkhdf,
        vtu: a.vtu,
        out: a.out,
    };
    let (o, code) = tools::run::<tools::mesh::Process>(input, &mut ctx(out))?;
    if !o.valid {
        eprintln!("the mesh is invalid (see the report)");
    }
    Ok(Done {
        code,
        json: serde_json::to_value(o)?,
    })
}

fn scalar(a: ScalarArgs, out: Out) -> Result<Done> {
    let input = tools::study::ScalarInput {
        case: a.case,
        backends: a.backends,
        threads: a.threads,
        rust_ptx: a.rust_ptx,
        levels: a.levels,
        out: a.out,
    };
    let (o, _) = tools::run::<tools::study::Scalar>(input, &mut ctx(out))?;
    Ok(failures(&o.failures, serde_json::to_value(&o)?))
}

/// Exit code 2 with the failed checks on standard error, else 0.
fn failures(f: &[String], json: serde_json::Value) -> Done {
    if !f.is_empty() {
        eprintln!("{} check(s) failed:\n  {}", f.len(), f.join("\n  "));
    }
    Done {
        code: if f.is_empty() { 0 } else { 2 },
        json,
    }
}

/// The typed input of `case.check` and `flow.run` from the options.
fn flow_input(a: &FlowArgs) -> tools::case::FlowInput {
    tools::case::FlowInput {
        case: a.case.clone(),
        set: a.set.clone(),
        out: a.out.clone(),
        threads: a.threads,
        levels: a.levels.clone(),
        max_iterations: a.max_iterations,
        no_vtkhdf: a.no_vtkhdf,
        restart: a.restart.clone(),
        end_time: a.end_time,
        backend: a.backend.clone(),
        ptx: a.ptx.clone(),
    }
}

fn flow(a: FlowArgs, out: Out) -> Result<Done> {
    let roots = Roots {
        read: a.read_root.clone(),
        write: a.write_root.clone(),
    };
    let input = flow_input(&a);
    if a.check {
        let (rep, code) =
            tools::run::<tools::case::Check>(input, &mut ctx(out).roots(roots.clone()))?;
        return Ok(Done {
            code,
            json: serde_json::to_value(&rep)?,
        });
    }
    let mut c = ctx(out).roots(roots);
    c.interrupts = signals();
    c.control = a.control_port;
    let o = tools::flow::run_outcome(input, &mut c)?;
    let s = &o.summary;
    for d in s
        .diagnostics
        .iter()
        .filter(|d| o.report.is_none() && d.severity == diag::Severity::Error)
    {
        eprintln!("Error {}: {}", d.code, d.message);
    }
    if let Some(r) = &o.report
        && !r.failures.is_empty()
    {
        eprintln!(
            "{} check(s) failed:\n  {}",
            r.failures.len(),
            r.failures.join("\n  ")
        );
    }
    Ok(Done {
        code: s.exit_code,
        json: serde_json::to_value(s)?,
    })
}

fn poisson(a: PoissonArgs, out: Out) -> Result<Done> {
    let input = tools::study::PoissonInput {
        case: a.case,
        backends: a.backends,
        solvers: a.solvers,
        threads: a.threads,
        rust_ptx: a.rust_ptx,
        levels: a.levels,
        repeat: a.repeat,
        out: a.out,
    };
    let (o, _) = tools::run::<tools::study::Poisson>(input, &mut ctx(out))?;
    Ok(failures(&o.failures, serde_json::to_value(&o)?))
}

fn g1(a: G1Args, out: Out) -> Result<Done> {
    let input = tools::study::G1Input {
        case: a.case,
        cells: a.cells,
        ordering: a.ordering,
        backends: a.backends,
        repeat: a.repeat,
        threads: a.threads,
        rust_ptx: a.rust_ptx,
        out: a.out,
    };
    let (o, _) = tools::run::<tools::study::G1>(input, &mut ctx(out))?;
    Ok(failures(&o.failures, serde_json::to_value(&o)?))
}
