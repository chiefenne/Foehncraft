//! `foehn`: the command line of Foehncraft (D4).
//!
//! Commands: `info` (build and platform), `mesh` (build or read a mesh,
//! check, partition and write it; stage 1), `scalar` (steady scalar
//! transport, verification sweeps, CPU and GPU; stage 2), `poisson` (the
//! pressure solver study: CG with algebraic multigrid against hypre; stage
//! 3), `flow` (steady incompressible laminar flow with SIMPLEC; stage 4)
//! and `bench g1` (the GPU kernel gate of stage 0), `schema` (the JSON
//! Schema of the flow case file, D77), `post derive` (derived fields from
//! a restart file, D72), `materials` (the material database, D94).
//!
//! Output streams (D79): human text on standard output, or with `--json`
//! on standard error and one JSON object on standard output (the result
//! summary of a flow run, the check report of `--check`, a command summary
//! for the others, an error object with a diagnostic code on failure).
//! `--quiet` silences the human text. No command asks questions. Exit
//! codes (D64): 0 finished with every check passed, 1 error, 2 finished
//! with failed checks, 3 stopped on request (signal, `.stop`, a limit).

use anyhow::Result;
use clap::{Parser, Subcommand};
use foehncraft_bench::G1Case;
use foehncraft_flow::diag;
use foehncraft_flow::engine::{self, BuildInfo, Interrupts, Roots, RunOptions};
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

impl Done {
    fn ok(json: serde_json::Value) -> Self {
        Self { code: 0, json }
    }
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
}

#[derive(Subcommand)]
enum CaseCmd {
    /// Rewrite case files in the current case format (renamed keys, the
    /// line `case_format`), keeping their comments. Without `--write` only
    /// the changes are listed; `--check` fails when a file needs one.
    Upgrade(UpgradeArgs),
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
    /// Mesh case file (TOML) or a mesh file (Gmsh or industry `.msh`, SU2
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
        Command::Post {
            post: Post::Derive(a),
        } => derive(a, out),
        Command::Bench {
            bench: Bench::Backends,
        } => {
            let b = foehncraft_bench::backend::available();
            for x in &b {
                out.line(x);
            }
            Ok(Done::ok(json!({ "backends": b })))
        }
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

/// `foehn post derive`.
fn derive(a: DeriveArgs, out: Out) -> Result<Done> {
    use foehncraft_flow::case::output_dir;
    use foehncraft_flow::derived::{self, Derived};
    let roots = Roots {
        read: a.read_root.clone(),
        write: a.write_root.clone(),
    };
    let case = engine::load(&a.case, &a.set, &roots)?;
    case.validate()?;
    let (_, dir) = output_dir(&case, &a.case, a.out.clone());
    let restart = if a.restart == std::path::Path::new("latest") {
        engine::latest_restart(&dir)?.ok_or_else(|| {
            diag::err(diag::Diagnostic::new(
                diag::RUN_NO_RESTART,
                format!(
                    "--restart latest: no restart file (restart.h5, restart_<N>.h5) in {}",
                    dir.display()
                ),
            ))
        })?
    } else {
        a.restart.clone()
    };
    roots.check_read(&restart, "the restart file")?;
    let list: Vec<Derived> = match &a.fields {
        Some(names) => names
            .iter()
            .map(|n| {
                Derived::parse(n.trim()).ok_or_else(|| {
                    let all: Vec<&str> = Derived::ALL.iter().map(|d| d.name()).collect();
                    diag::err(diag::Diagnostic::new(
                        diag::CASE_BAD_VALUE,
                        format!("unknown derived field `{n}`; known: {}", all.join(", ")),
                    ))
                })
            })
            .collect::<Result<_>>()?,
        None if !case.output.derived.is_empty() => case.output.derived.clone(),
        // every field; the pressure coefficient only with a reference
        None => Derived::ALL
            .into_iter()
            .filter(|d| {
                *d != Derived::PressureCoefficient || derived::validate_list(&case, &[*d]).is_ok()
            })
            .collect(),
    };
    let file = a.file.clone().unwrap_or_else(|| {
        let stem = restart
            .file_stem()
            .map_or("restart".into(), |s| s.to_string_lossy().into_owned());
        let name = match stem.strip_prefix("restart") {
            Some(rest) => format!("derived{rest}.vtkhdf"),
            None => format!("{stem}_derived.vtkhdf"),
        };
        restart.with_file_name(name)
    });
    if let Some(d) = file.parent().filter(|d| !d.as_os_str().is_empty()) {
        roots.check_write(d)?;
    }
    let rep = derived::from_restart(&case, &restart, &list, &file)?;
    let at = match (rep.step, rep.time) {
        (Some(n), Some(t)) => format!("step {n}, t = {t} s"),
        _ => format!("iteration {}", rep.iteration),
    };
    out.line(&format!(
        "{}: {} cells at {at} from {}",
        rep.file.display(),
        rep.cells,
        rep.restart.display()
    ));
    out.line(&format!(
        "  derived: {}",
        rep.derived
            .iter()
            .map(|d| d.name())
            .collect::<Vec<_>>()
            .join(", ")
    ));
    Ok(Done::ok(serde_json::to_value(&rep)?))
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
        platform: foehncraft_bench::sysinfo::platform(),
        cpu: foehncraft_bench::sysinfo::cpu_model(),
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
    use foehncraft_flow::materials as m;
    anyhow::ensure!(
        a.temperature > 0.0 && a.pressure > 0.0,
        "--temperature and --pressure must be > 0"
    );
    let cond = m::Conditions {
        temperature: a.temperature,
        pressure: a.pressure,
    };
    let mut lib = m::CaseMaterials {
        include: a.include.clone(),
        ..Default::default()
    };
    lib.load()?;
    let names: Vec<String> = match &a.name {
        Some(n) => vec![n.clone()],
        None => lib.names(),
    };
    let mut list = Vec::new();
    for n in &names {
        list.push(m::describe_in(n, &lib, cond)?);
    }
    let num = |d: &serde_json::Value, k: &str| -> Option<f64> {
        d[k]["value"]
            .as_f64()
            .or_else(|| d[k]["tensor"][0][0].as_f64())
    };
    let cell = |v: Option<f64>, w: usize, p: usize| match v {
        Some(x) => format!("{x:>w$.p$e}"),
        None => format!("{:>w$}", "–"),
    };
    if a.name.is_none() {
        out.line(&format!(
            "material database ({} materials) at T = {} K, p = {} Pa",
            list.len(),
            a.temperature,
            a.pressure
        ));
        out.line(&format!(
            "  {:<18} {:<7} {:>11} {:>11} {:>11} {:>11}  {:<16} description",
            "name", "kind", "ρ (kg/m³)", "μ (Pa s)", "c_p (J/kgK)", "k (W/mK)", "range (K)"
        ));
        for d in &list {
            let range = d["temperature_range"]
                .as_array()
                .map(|r| format!("{}–{}", r[0], r[1]))
                .unwrap_or_default();
            out.line(&format!(
                "  {:<18} {:<7} {} {} {} {}  {:<16} {}",
                d["name"].as_str().unwrap_or(""),
                d["kind"].as_str().unwrap_or(""),
                cell(num(d, "density"), 11, 4),
                cell(num(d, "viscosity"), 11, 3),
                cell(num(d, "specific_heat"), 11, 4),
                cell(num(d, "conductivity"), 11, 3),
                range,
                d["description"].as_str().unwrap_or("")
            ));
        }
    } else if let Some(d) = list.first() {
        out.line(&format!(
            "{}: {}",
            d["name"].as_str().unwrap_or(""),
            d["description"].as_str().unwrap_or("")
        ));
        for what in ["density", "viscosity", "specific_heat", "conductivity"] {
            if d[what].is_null() {
                continue;
            }
            let v = match num(d, what) {
                Some(x) => format!("{x}"),
                None => "–".into(),
            };
            out.line(&format!(
                "  {what:<14} {v} ({})",
                d[what]["origin"].as_str().unwrap_or("")
            ));
        }
        for w in d["warnings"].as_array().into_iter().flatten() {
            out.line(&format!("  warning: {}", w.as_str().unwrap_or("")));
        }
    }
    Ok(Done::ok(json!({
        "format": "foehncraft.materials",
        "format_version": engine::FORMAT_VERSION,
        "database_sha256": m::database_hash(),
        "materials": list,
    })))
}

/// `foehn case upgrade`.
fn upgrade(a: UpgradeArgs, out: Out) -> Result<Done> {
    use anyhow::Context as _;
    use foehncraft_mesh::case_format as cf;
    let mut pending = 0usize;
    let mut failed = 0usize;
    let mut rows = Vec::new();
    for f in &a.files {
        let text =
            std::fs::read_to_string(f).with_context(|| format!("reading {}", f.display()))?;
        match cf::upgrade_text(&text) {
            Err(e) => {
                failed += 1;
                out.line(&format!("{}: {e}", f.display()));
                rows.push(json!({ "file": f.display().to_string(), "error": e }));
            }
            Ok((new, notes)) => {
                let changed = new != text;
                if changed {
                    pending += 1;
                }
                if !changed {
                    out.line(&format!(
                        "{}: format {}, nothing to do",
                        f.display(),
                        cf::CURRENT
                    ));
                } else {
                    out.line(&format!("{}:", f.display()));
                    for n in &notes {
                        out.line(&format!("  {n}"));
                    }
                    if notes.is_empty() {
                        out.line(&format!(
                            "  the line case_format = \"{}\" is added",
                            cf::CURRENT
                        ));
                    }
                    if a.write {
                        std::fs::write(f, &new)
                            .with_context(|| format!("writing {}", f.display()))?;
                        out.line("  written");
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
    let code = i32::from(failed > 0 || (a.check && pending > 0));
    if pending > 0 && !a.write {
        out.line("run again with --write to rewrite the files");
    }
    Ok(Done {
        code,
        json: json!({
            "format": "foehncraft.case_upgrade",
            "format_version": engine::FORMAT_VERSION,
            "current": cf::CURRENT,
            "files": rows,
        }),
    })
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
    let s = schema::json_schema();
    if a.schema_only || !out.json {
        stdout(&serde_json::to_string_pretty(&s)?);
        return Ok(Done {
            code: 0,
            json: serde_json::Value::Null,
        });
    }
    Ok(Done::ok(schema::envelope(s)))
}

fn info(out: Out) -> Result<Done> {
    let b = build_info();
    out.line(&format!("Foehncraft {}", foehncraft_core::VERSION));
    out.line(&format!("commit       {}", b.commit));
    out.line(&format!("platform     {}", b.platform));
    out.line(&format!("cpu          {}", b.cpu));
    out.line(&format!("threads      {}", rayon::current_num_threads()));
    out.line("precision    f64 (reference path, D7)");
    let ranks = {
        use foehncraft_parallel::Communicator as _;
        let c = foehncraft_parallel::SerialComm;
        out.line(&format!(
            "ranks        {} ({} communicator; MPI runs start with mpirun)",
            c.size(),
            c.name()
        ));
        c.size()
    };
    out.line(&format!(
        "features     {}",
        if b.features.is_empty() {
            "none".to_string()
        } else {
            b.features.join(", ")
        }
    ));
    let backends = foehncraft_bench::backend::available();
    out.line(&format!("backends     {}", backends.join(", ")));
    Ok(Done::ok(json!({
        "format": "foehncraft.info",
        "format_version": engine::FORMAT_VERSION,
        "version": foehncraft_core::VERSION,
        "commit": b.commit,
        "platform": b.platform,
        "cpu": b.cpu,
        "threads": rayon::current_num_threads(),
        "ranks": ranks,
        "features": b.features,
        "backends": backends,
        "flow_backends": foehncraft_flow::device::available(),
    })))
}

fn mesh(a: MeshArgs, out: Out) -> Result<Done> {
    use foehncraft_mesh::case::MeshCase;
    let mut case = if foehncraft_mesh::format::is_mesh_file(&a.input) {
        MeshCase::for_file(&a.input)
    } else {
        MeshCase::from_file(&a.input).map_err(diag::case_error)?
    };
    let name = a
        .input
        .file_stem()
        .map_or("mesh".into(), |s| s.to_string_lossy().into_owned());
    if let Some(p) = a.parts {
        case.partition.parts = p;
    }
    if let Some(m) = a.method {
        case.partition.method = match m.as_str() {
            "metis" => foehncraft_mesh::partition::Method::Metis,
            "rcb" => foehncraft_mesh::partition::Method::Rcb,
            _ => anyhow::bail!("--method: metis or rcb"),
        };
    }
    if let Some(r) = a.renumber {
        use foehncraft_mesh::renumber::Ordering;
        case.mesh.renumber = match r.as_str() {
            "none" => Ordering::None,
            "rcm" => Ordering::Rcm,
            "shuffle" => Ordering::Shuffle,
            _ => anyhow::bail!("--renumber: none, rcm or shuffle"),
        };
    }
    if let Some(n) = a.locate {
        case.check.locate_points = n;
    }
    if a.no_vtkhdf {
        case.output.vtkhdf = false;
    }
    if a.vtu {
        case.output.vtu = true;
    }
    if a.out.is_some() {
        case.output.directory = a.out;
    }
    let (report, files) = foehncraft_mesh::case::run(&case, &name, &mut |l| out.line(l))?;
    for f in &files {
        out.line(&format!("wrote {}", f.display()));
    }
    let ok = report.quality.ok();
    if !ok {
        eprintln!("the mesh is invalid (see the report)");
    }
    Ok(Done {
        code: if ok { 0 } else { 2 },
        json: json!({ "valid": ok, "files": files }),
    })
}

fn scalar(a: ScalarArgs, out: Out) -> Result<Done> {
    use foehncraft_fv::case::{ScalarCase, output_dir, run};
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
    let rep = run(&case, &name, &dir, |l| out.line(l))?;
    if !rep.orders.is_empty() {
        out.line(&format!(
            "observed order: {}",
            rep.orders
                .iter()
                .map(|o| format!("{o:.3}"))
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    out.line(&format!("wrote {}", dir.join("scalar.md").display()));
    Ok(failures(
        &rep.failures,
        json!({
            "orders": rep.orders,
            "failures": rep.failures,
            "report": dir.join("scalar.md"),
        }),
    ))
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

/// The CLI options that override case keys.
fn flow_overrides(case: &mut foehncraft_flow::case::FlowCase, a: &FlowArgs) -> Result<()> {
    if a.threads.is_some() {
        case.run.threads = a.threads;
    }
    if let Some(l) = &a.levels {
        case.verification.levels = l.clone();
    }
    if let Some(n) = a.max_iterations {
        case.solver.max_iterations = n;
    }
    if a.no_vtkhdf {
        case.output.vtkhdf = false;
    }
    if let Some(b) = &a.backend {
        case.run.backend = b.clone();
    }
    if let Some(p) = &a.ptx {
        case.run.ptx = p.clone();
    }
    if let Some(r) = &a.restart {
        case.run.restart = Some(r.clone());
    }
    if let Some(e) = a.end_time {
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

fn flow(a: FlowArgs, out: Out) -> Result<Done> {
    use foehncraft_flow::case::output_dir;
    let roots = Roots {
        read: a.read_root.clone(),
        write: a.write_root.clone(),
    };
    let build = build_info();
    let loaded = engine::load(&a.case, &a.set, &roots).and_then(|mut c| {
        flow_overrides(&mut c, &a)?;
        Ok(c)
    });
    if a.check {
        let mut rep = engine::CheckReport {
            format: "foehncraft.check",
            format_version: engine::FORMAT_VERSION,
            case_file: a.case.display().to_string(),
            ok: false,
            diagnostics: Vec::new(),
            estimate: None,
        };
        match &loaded {
            Ok(case) => {
                engine::check_case(case, &a.case, a.out.as_deref(), &roots, &build, &mut rep)
            }
            Err(e) => rep.diagnostics.push(diag::of_error(e)),
        }
        out.line(&format!("{}:", a.case.display()));
        for d in &rep.diagnostics {
            out.line(&format!(
                "{} {}: {}{}",
                d.severity.name(),
                d.code,
                d.boundary
                    .as_ref()
                    .map_or(String::new(), |b| format!("boundary `{b}`: ")),
                d.message
            ));
            if let Some(f) = &d.fix {
                out.line(&format!("  fix: --set {}", f.set.join(" --set ")));
            }
        }
        if let Some(e) = &rep.estimate {
            out.line(&estimate_text(e));
        }
        out.line(if rep.ok {
            "the case can be run"
        } else {
            "the case cannot be run"
        });
        return Ok(Done {
            code: if rep.ok { 0 } else { 1 },
            json: serde_json::to_value(&rep)?,
        });
    }
    let case = loaded?;
    let (name, dir) = output_dir(&case, &a.case, a.out.clone());
    let opts = RunOptions {
        interrupts: signals(),
        build,
        roots,
        observer: None,
    };
    let o = engine::run(&case, &name, &dir, &opts, |l| out.line(l))?;
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

/// The estimate of `--check` as text.
fn estimate_text(e: &engine::Estimate) -> String {
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
    s
}

fn poisson(a: PoissonArgs, out: Out) -> Result<Done> {
    use foehncraft_fv::poisson::{PoissonCase, StudySolver, output_dir, run};
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
    let rep = run(&case, &name, &dir, |l| out.line(l))?;
    for (n, r) in &rep.iteration_ratios {
        out.line(&format!("iterations finest / coarsest level, {n}: {r:.2}"));
    }
    out.line(&format!("wrote {}", dir.join("poisson.md").display()));
    Ok(failures(
        &rep.failures,
        json!({
            "iteration_ratios": rep.iteration_ratios,
            "failures": rep.failures,
            "report": dir.join("poisson.md"),
        }),
    ))
}

fn g1(a: G1Args, out: Out) -> Result<Done> {
    let mut case = match &a.case {
        Some(p) => G1Case::from_file(p).map_err(diag::case_error)?,
        None => G1Case::default(),
    };
    if let Some(n) = a.cells {
        case.mesh.cells = n;
        case.mesh.nx = None;
    }
    if let Some(o) = a.ordering {
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
    let report = foehncraft_bench::run_g1(&case, &mut |line| out.line(line))?;
    for v in &report.verdicts {
        out.line(&format!(
            "G1 {} against {}: {}",
            v.candidate,
            v.reference,
            if v.pass { "PASS" } else { "FAIL" }
        ));
    }
    let files = report.write(&case.output.directory, &case.output.name)?;
    for f in &files {
        out.line(&format!("wrote {}", f.display()));
    }
    let all_correct = report
        .results
        .iter()
        .all(|b| b.kernels.iter().all(|k| k.within_tolerance));
    let f: Vec<String> = if all_correct {
        Vec::new()
    } else {
        vec!["a backend did not reproduce the reference results".into()]
    };
    Ok(failures(
        &f,
        json!({
            "verdicts": report.verdicts.iter().map(|v| json!({
                "candidate": v.candidate, "reference": v.reference, "pass": v.pass,
            })).collect::<Vec<_>>(),
            "files": files,
        }),
    ))
}
