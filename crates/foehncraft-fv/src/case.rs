//! The `foehn scalar` case: a steady scalar transport problem on one
//! mesh, or an order-of-accuracy sweep over refinement levels, on one or
//! more backends.
//!
//! Case file (TOML, every key optional unless noted, unknown keys are
//! errors):
//!
//! ```toml
//! [mesh]                       # as in `foehn mesh` (file or generate)
//! [mesh.generate]
//! type = "poly"
//! cells = [8, 8, 8]
//! [problem]                    # crate::problem::Problem
//! diffusivity = 0.2
//! velocity = [1.0, 0.6, 0.3]
//! manufactured = {}            # source and boundary values of the MMS
//! [problem.boundary.xmin]      # one table per boundary (required)
//! type = "fixed-value"         # fixed-value, fixed-gradient, zero-gradient
//! profile = "manufactured"     # constant (value), step, linear, manufactured
//! [solver]                     # crate::solver::Settings
//! limiter = "venkatakrishnan"
//! [run]
//! backends = ["cpu", "cuda"]   # results are compared with the first one
//! [verification]
//! levels = [1, 2, 3]           # refine levels of [mesh.generate]: a sweep
//! [[variants]]                 # optional: the same case with other solver settings
//! name = "upwind"
//! solver = { scheme = "upwind" }
//! [output]
//! vtkhdf = true
//! profile = { from = [0.75, 0.0, 0.05], to = [0.75, 1.0, 0.05], points = 200 }
//! ```

use crate::backend::{self, V};
use crate::disc::Discretisation;
use crate::problem::Problem;
use crate::solver::{self, Outcome, Settings};
use crate::verify::{MmsResult, observed_orders};
use anyhow::{Context, Result, bail, ensure};
use foehncraft_mesh::case::{MeshSection, build_mesh};
use foehncraft_mesh::locate::Locator;
use foehncraft_mesh::renumber::renumber;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// A scalar case.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ScalarCase {
    /// Format of the case file (`case_format`).
    pub case_format: Option<String>,
    /// The mesh.
    pub mesh: MeshSection,
    /// The problem.
    pub problem: Problem,
    /// Solver settings.
    pub solver: Settings,
    /// Backends.
    pub run: RunSection,
    /// Order-of-accuracy sweep.
    pub verification: VerificationSection,
    /// Variants of the solver settings.
    pub variants: Vec<Variant>,
    /// Output.
    pub output: OutputSection,
}

/// `[run]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RunSection {
    /// Backends: `cpu`, `cpu-serial`, `cuda`.
    pub backends: Vec<String>,
    /// CPU threads (default: all cores).
    pub threads: Option<usize>,
    /// PTX module of the GPU backend.
    pub rust_ptx: PathBuf,
    /// Also solve every variant on these numbers of ranks (MPI step W3):
    /// the mesh is cut by recursive coordinate bisection and every rank
    /// runs as a thread of this process with its own communicator, the
    /// threads of `threads` shared among them. Each solve must give
    /// bitwise the φ of the first backend, or the run fails.
    pub ranks: Vec<usize>,
}

impl Default for RunSection {
    fn default() -> Self {
        Self {
            backends: vec!["cpu".into()],
            ranks: Vec::new(),
            threads: None,
            rust_ptx: "target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx".into(),
        }
    }
}

/// `[verification]`
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct VerificationSection {
    /// Refinement levels of `[mesh.generate]`; two or more give an
    /// order-of-accuracy sweep (needs a manufactured solution).
    pub levels: Vec<u32>,
    /// Accepted range of the observed order on the finest pair; the run
    /// fails outside it.
    pub expected_order: Option<[f64; 2]>,
    /// Accepted range of the solution, for bounded problems; the run fails
    /// outside it (plus 1e-6).
    pub bounds: Option<[f64; 2]>,
}

/// `[[variants]]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Variant {
    /// Name (output sub-directory and legend).
    pub name: String,
    /// Settings that replace `[solver]` (unset keys take their defaults).
    #[serde(default)]
    pub solver: Settings,
    /// Accepted range of φ for this variant (default:
    /// `verification.bounds`).
    #[serde(default)]
    pub bounds: Option<[f64; 2]>,
}

/// A line along which φ is sampled.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Line {
    /// Start point.
    pub from: [f64; 3],
    /// End point.
    pub to: [f64; 3],
    /// Number of points.
    #[serde(default = "default_points")]
    pub points: usize,
}

fn default_points() -> usize {
    200
}

/// `[output]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct OutputSection {
    /// Output directory (default `output/<case name>`).
    pub directory: Option<PathBuf>,
    /// Write the fields (φ, error, limiter) as VTKHDF.
    pub vtkhdf: bool,
    /// Sample φ along a line (CSV and SVG).
    pub profile: Option<Line>,
}

impl Default for OutputSection {
    fn default() -> Self {
        Self {
            directory: None,
            vtkhdf: true,
            profile: None,
        }
    }
}

impl ScalarCase {
    /// Reads a case file; relative paths are taken from its directory.
    pub fn from_file(path: &Path) -> Result<Self> {
        let text =
            std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let (mut c, _notes): (Self, _) = foehncraft_mesh::case_format::read(&text)
            .with_context(|| format!("in case file {}", path.display()))?;
        c.case_format = Some(foehncraft_mesh::case_format::CURRENT.into());
        let dir = path.parent().unwrap_or(Path::new("."));
        if let Some(f) = &c.mesh.file {
            c.mesh.file = Some(dir.join(f));
        }
        c.problem.link();
        c.solver.validate()?;
        for v in &c.variants {
            v.solver
                .validate()
                .with_context(|| format!("variant {}", v.name))?;
        }
        ensure!(!c.run.backends.is_empty(), "[run] backends is empty");
        ensure!(
            c.run.ranks.iter().all(|&n| (1..=256).contains(&n)),
            "[run] ranks must lie in 1..=256"
        );
        if c.verification.levels.len() >= 2 {
            ensure!(
                c.mesh.generate.is_some(),
                "[verification] levels need a [mesh.generate] table"
            );
            ensure!(
                c.problem.manufactured.is_some(),
                "[verification] levels need [problem] manufactured"
            );
        }
        Ok(c)
    }
}

/// One solve on one backend.
#[derive(Debug, Clone, Serialize)]
pub struct Solve {
    /// Variant name (empty: the main settings).
    pub variant: String,
    /// Backend.
    pub backend: String,
    /// Device.
    pub device: String,
    /// Cells.
    pub cells: usize,
    /// Solver outcome.
    pub outcome: Outcome,
    /// ns per cell and outer iteration (R6).
    pub ns_per_cell_iteration: f64,
    /// Error norms against the manufactured solution.
    pub mms: Option<MmsResult>,
    /// Range of φ.
    pub range: [f64; 2],
    /// Largest difference to the first backend, relative to max |φ|.
    pub diff_to_first: Option<f64>,
    /// Bitwise equal to the first backend.
    pub bitwise_equal: Option<bool>,
}

/// Report of a scalar case.
#[derive(Debug, Clone, Serialize)]
pub struct ScalarReport {
    /// Case name.
    pub name: String,
    /// Mesh description per level.
    pub meshes: Vec<String>,
    /// All solves.
    pub solves: Vec<Solve>,
    /// Observed orders of the sweep (first backend).
    pub orders: Vec<f64>,
    /// Problems found by the checks.
    pub failures: Vec<String>,
    /// Seconds for everything.
    pub seconds: f64,
}

fn mesh_for(case: &ScalarCase, level: Option<u32>) -> Result<(Mesh, Geometry, String)> {
    let mut sec = case.mesh.clone();
    if let (Some(l), Some(g)) = (level, sec.generate.as_mut()) {
        g.refine = l;
    }
    let (mut mesh, source, _) = build_mesh(&sec)?;
    renumber(&mut mesh, sec.renumber, sec.seed).map_err(anyhow::Error::msg)?;
    let geo = Geometry::compute(&mesh);
    let desc = format!("{source}: {} cells", mesh.n_cells);
    Ok((mesh, geo, desc))
}

fn dimension(case: &ScalarCase) -> u8 {
    case.mesh.generate.as_ref().map_or(3, |g| g.dimension)
}

/// Runs a case; `log` receives progress lines; files go to `out`.
pub fn run(
    case: &ScalarCase,
    name: &str,
    out: &Path,
    mut log: impl FnMut(&str),
) -> Result<ScalarReport> {
    let t0 = Instant::now();
    std::fs::create_dir_all(out)?;
    let levels: Vec<Option<u32>> = if case.verification.levels.is_empty() {
        vec![None]
    } else {
        case.verification.levels.iter().map(|&l| Some(l)).collect()
    };
    let variants: Vec<(String, Settings, Option<[f64; 2]>)> = if case.variants.is_empty() {
        vec![(String::new(), case.solver.clone(), case.verification.bounds)]
    } else {
        case.variants
            .iter()
            .map(|v| {
                (
                    v.name.clone(),
                    v.solver.clone(),
                    v.bounds.or(case.verification.bounds),
                )
            })
            .collect()
    };
    let mut rep = ScalarReport {
        name: name.into(),
        meshes: Vec::new(),
        solves: Vec::new(),
        orders: Vec::new(),
        failures: Vec::new(),
        seconds: 0.0,
    };
    let mut profiles: Vec<(String, Vec<(f64, f64)>)> = Vec::new();
    let last_level = levels.len() - 1;
    for (li, level) in levels.iter().enumerate() {
        let (mesh, geo, desc) = mesh_for(case, *level)?;
        log(&format!("mesh  {desc}"));
        rep.meshes.push(desc);
        let disc = Discretisation::new(&mesh, &geo, &case.problem)?;
        for (vname, settings, bounds) in &variants {
            let mut first: Option<Vec<f64>> = None;
            for be_name in &case.run.backends {
                let tb = Instant::now();
                let mut be = backend::create(be_name, &disc, case.run.threads, &case.run.rust_ptx)?;
                let device = be.device();
                log(&format!(
                    "solve {}{} on {be_name} ({device})",
                    if vname.is_empty() { "" } else { "variant " },
                    vname
                ));
                let phi0 = vec![case.problem.manufactured.map_or(0.0, |m| m.offset); mesh.n_cells];
                be.upload(V::Phi, &phi0)?;
                let outcome = solver::solve(be.as_mut(), settings, |l| log(&format!("  {l}")))?;
                let phi = be.download(V::Phi)?;
                let psi = be.download(V::Psi)?;
                if !outcome.converged {
                    rep.failures.push(format!(
                        "{be_name} {vname}: not converged in {} outer iterations (residual {:.2e})",
                        outcome.outer,
                        outcome.residuals.last().copied().unwrap_or(f64::NAN)
                    ));
                }
                let range = [
                    phi.iter().copied().fold(f64::INFINITY, f64::min),
                    phi.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                ];
                if let Some([lo, hi]) = *bounds
                    && (range[0] < lo - 1e-6 || range[1] > hi + 1e-6)
                {
                    rep.failures.push(format!(
                        "{be_name} {vname}: φ in [{:.6}, {:.6}], outside [{lo}, {hi}]",
                        range[0], range[1]
                    ));
                }
                let mms = case.problem.manufactured.map(|m| {
                    let (mut e2, mut einf, mut vol) = (0.0, 0.0f64, 0.0);
                    for (c, &v) in phi.iter().enumerate() {
                        let e = v - m.phi(geo.cell_centroid[c]);
                        e2 += e * e * geo.cell_volume[c];
                        einf = einf.max(e.abs());
                        vol += geo.cell_volume[c];
                    }
                    let n = mesh.n_cells as f64;
                    let h = if dimension(case) == 2 {
                        let (lo, hi) = mesh.bounding_box();
                        (vol / (hi[2] - lo[2]) / n).sqrt()
                    } else {
                        (vol / n).cbrt()
                    };
                    MmsResult {
                        kind: case
                            .mesh
                            .generate
                            .as_ref()
                            .map_or("file".into(), |g| g.kind.name().into()),
                        cells: mesh.n_cells,
                        h,
                        l2: (e2 / vol).sqrt(),
                        linf: einf,
                        max_nonorth_deg: disc.max_nonorth_deg,
                        outer: outcome.outer,
                        converged: outcome.converged,
                        timings: outcome.timings.clone(),
                        ns_per_cell_iteration: outcome.ns_per_cell_iteration(mesh.n_cells),
                    }
                });
                let (diff, bitwise) = match &first {
                    None => (None, None),
                    Some(f) => {
                        let scale = f.iter().fold(0.0f64, |a, &b| a.max(b.abs())).max(1e-300);
                        let d = f
                            .iter()
                            .zip(&phi)
                            .fold(0.0f64, |a, (x, y)| a.max((x - y).abs()))
                            / scale;
                        let eq = f.iter().zip(&phi).all(|(x, y)| x.to_bits() == y.to_bits());
                        (Some(d), Some(eq))
                    }
                };
                if let Some(d) = diff
                    && d > 1e-10
                {
                    rep.failures.push(format!(
                        "{be_name} {vname}: differs from {} by {d:.2e} (relative)",
                        case.run.backends[0]
                    ));
                }
                log(&format!(
                    "      {} outer iterations, {:.2} s, {:.0} ns/cell/iteration{}{}",
                    outcome.outer,
                    tb.elapsed().as_secs_f64(),
                    outcome.ns_per_cell_iteration(mesh.n_cells),
                    mms.as_ref()
                        .map_or(String::new(), |m| format!(", L2 error {:.3e}", m.l2)),
                    match bitwise {
                        Some(true) => ", bitwise equal to the first backend".to_string(),
                        Some(false) => format!(", max rel. difference {:.2e}", diff.unwrap()),
                        None => String::new(),
                    }
                ));
                rep.solves.push(Solve {
                    variant: vname.clone(),
                    backend: be_name.clone(),
                    device,
                    cells: mesh.n_cells,
                    ns_per_cell_iteration: outcome.ns_per_cell_iteration(mesh.n_cells),
                    outcome,
                    mms,
                    range,
                    diff_to_first: diff,
                    bitwise_equal: bitwise,
                });
                if first.is_none() {
                    if li == last_level {
                        let dir = if vname.is_empty() {
                            out.to_path_buf()
                        } else {
                            out.join(vname)
                        };
                        std::fs::create_dir_all(&dir)?;
                        let rank = case.run.ranks.iter().max().map(|&n| {
                            rcb(&mesh, &geo, n)
                                .iter()
                                .map(|&r| r as f64)
                                .collect::<Vec<_>>()
                        });
                        write_fields(case, &dir, &mesh, &geo, &phi, &psi, rank.as_deref())?;
                        if let Some(line) = &case.output.profile {
                            let p = sample(&mesh, &geo, &phi, line);
                            let mut csv = String::from("s,phi\n");
                            for (s, v) in &p {
                                writeln!(csv, "{s},{v}")?;
                            }
                            std::fs::write(dir.join("profile.csv"), csv)?;
                            profiles.push((
                                if vname.is_empty() {
                                    "φ".into()
                                } else {
                                    vname.clone()
                                },
                                p,
                            ));
                        }
                    }
                    first = Some(phi);
                }
            }
            for &n in &case.run.ranks {
                let solve = solve_on_ranks(case, &mesh, &geo, settings, n, first.as_deref())?;
                log(&format!(
                    "solve {}{} on {n} ranks: {} outer iterations, {}",
                    if vname.is_empty() { "" } else { "variant " },
                    vname,
                    solve.outcome.outer,
                    if solve.bitwise_equal == Some(true) {
                        "bitwise equal to the first backend".to_string()
                    } else {
                        format!(
                            "DIFFERS from the first backend ({:.2e} relative)",
                            solve.diff_to_first.unwrap_or(f64::NAN)
                        )
                    }
                ));
                if solve.bitwise_equal != Some(true) {
                    rep.failures.push(format!(
                        "{vname} on {n} ranks: not bitwise equal to {}",
                        case.run.backends[0]
                    ));
                }
                rep.solves.push(Solve {
                    variant: vname.clone(),
                    ..solve
                });
            }
        }
    }
    // order of accuracy of the first backend and variant
    if levels.len() >= 2 {
        let r: Vec<MmsResult> = rep
            .solves
            .iter()
            .filter(|s| s.backend == case.run.backends[0] && s.variant == variants[0].0)
            .filter_map(|s| s.mms.clone())
            .collect();
        rep.orders = observed_orders(&r, dimension(case));
        if let (Some([lo, hi]), Some(&last)) = (case.verification.expected_order, rep.orders.last())
            && !(lo..=hi).contains(&last)
        {
            rep.failures.push(format!(
                "observed order {last:.3} on the finest pair, outside [{lo}, {hi}]"
            ));
        }
        std::fs::write(
            out.join("convergence.svg"),
            convergence_svg(name, &r, &rep.orders),
        )?;
    }
    if !profiles.is_empty() {
        std::fs::write(out.join("profile.svg"), profile_svg(name, &profiles))?;
    }
    rep.seconds = t0.elapsed().as_secs_f64();
    std::fs::write(out.join("scalar.json"), serde_json::to_string_pretty(&rep)?)?;
    std::fs::write(out.join("scalar.md"), rep.to_markdown())?;
    Ok(rep)
}

#[cfg(feature = "hdf5")]
fn write_fields(
    case: &ScalarCase,
    dir: &Path,
    mesh: &Mesh,
    geo: &Geometry,
    phi: &[f64],
    psi: &[f64],
    rank: Option<&[f64]>,
) -> Result<()> {
    use foehncraft_mesh::vtkhdf::{CellField, write};
    if !case.output.vtkhdf {
        return Ok(());
    }
    let mut fields = vec![
        CellField {
            name: "phi",
            values: phi,
        },
        CellField {
            name: "limiter",
            values: psi,
        },
    ];
    let err: Vec<f64>;
    if let Some(m) = case.problem.manufactured {
        err = phi
            .iter()
            .enumerate()
            .map(|(c, v)| v - m.phi(geo.cell_centroid[c]))
            .collect();
        fields.push(CellField {
            name: "error",
            values: &err,
        });
    }
    if let Some(r) = rank {
        fields.push(CellField {
            name: "rank",
            values: r,
        });
    }
    write(dir.join("fields.vtkhdf"), mesh, &fields)
}

#[cfg(not(feature = "hdf5"))]
fn write_fields(
    _: &ScalarCase,
    _: &Path,
    _: &Mesh,
    _: &Geometry,
    _: &[f64],
    _: &[f64],
    _: Option<&[f64]>,
) -> Result<()> {
    Ok(())
}

/// The partition of `mesh` into `n` ranks of the `[run] ranks` solves:
/// recursive coordinate bisection (deterministic, no C library).
fn rcb(mesh: &Mesh, geo: &Geometry, n: usize) -> Vec<u32> {
    use foehncraft_mesh::partition::{Method, partition};
    partition(mesh, geo, n, Method::Rcb).unwrap_or_else(|_| vec![0; mesh.n_cells])
}

/// One variant solved on `n` ranks, threads of this process (MPI step
/// W3), φ gathered by global cell identity and compared bit for bit with
/// `first` (the first backend's φ).
fn solve_on_ranks(
    case: &ScalarCase,
    mesh: &Mesh,
    geo: &Geometry,
    settings: &Settings,
    n: usize,
    first: Option<&[f64]>,
) -> Result<Solve> {
    use foehncraft_mesh::decompose::SharedComm;
    use foehncraft_parallel::ThreadComm;
    let t = Instant::now();
    let part = rcb(mesh, geo, n);
    let cores = case
        .run
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |c| c.get()));
    let per_rank = (cores / n).max(1);
    let solves = ThreadComm::run(n, |c| {
        let comm: SharedComm = c;
        crate::ranks::solve_part(
            mesh,
            &part,
            1,
            comm,
            &case.problem,
            settings,
            None,
            Some(per_rank),
        )
    });
    let at: std::collections::HashMap<u64, usize> = mesh
        .cell_id
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let mut phi = vec![f64::NAN; mesh.n_cells];
    let mut outcome = None;
    for s in solves {
        let s = s?;
        for (id, v) in s.cell_id.iter().zip(&s.phi) {
            phi[at[id]] = *v;
        }
        outcome.get_or_insert(s.outcome);
    }
    let outcome = outcome.unwrap_or_default();
    let range = [
        phi.iter().copied().fold(f64::INFINITY, f64::min),
        phi.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    ];
    let (diff, bitwise) = match first {
        None => (None, None),
        Some(f) => {
            let scale = f.iter().fold(0.0f64, |a, &b| a.max(b.abs())).max(1e-300);
            let d = f
                .iter()
                .zip(&phi)
                .fold(0.0f64, |a, (x, y)| a.max((x - y).abs()))
                / scale;
            let eq = f.iter().zip(&phi).all(|(x, y)| x.to_bits() == y.to_bits());
            (Some(d), Some(eq))
        }
    };
    let seconds = t.elapsed().as_secs_f64();
    Ok(Solve {
        variant: String::new(),
        backend: format!("cpu, {n} ranks"),
        device: format!("{n} ranks as threads of one process, {per_rank} threads each, RCB"),
        cells: mesh.n_cells,
        ns_per_cell_iteration: 1e9 * seconds / (mesh.n_cells as f64 * outcome.outer.max(1) as f64),
        outcome,
        mms: None,
        range,
        diff_to_first: diff,
        bitwise_equal: bitwise,
    })
}

/// φ of the cells that contain the points of `line` (cell values, no
/// interpolation): pairs (distance along the line, φ).
fn sample(mesh: &Mesh, geo: &Geometry, phi: &[f64], line: &Line) -> Vec<(f64, f64)> {
    let loc = Locator::new(mesh, geo);
    let n = line.points.max(2);
    let d: [f64; 3] = std::array::from_fn(|i| line.to[i] - line.from[i]);
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let mut out = Vec::new();
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let x = std::array::from_fn(|k| line.from[k] + t * d[k]);
        if let Some(c) = loc.locate_search(x) {
            out.push((t * len, phi[c]));
        }
    }
    out
}

impl ScalarReport {
    /// Markdown report.
    pub fn to_markdown(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "# Scalar case `{}`\n", self.name);
        for m in &self.meshes {
            let _ = writeln!(s, "- mesh: {m}");
        }
        let _ = writeln!(s, "- total time: {:.2} s\n", self.seconds);
        let _ = writeln!(
            s,
            "| Variant | Backend | Device | Cells | Outer it. | Linear it. | Final residual | L2 error | max error | φ range | ns/cell/it. | vs first backend |"
        );
        let _ = writeln!(
            s,
            "|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---|"
        );
        for v in &self.solves {
            let _ = writeln!(
                s,
                "| {} | {} | {} | {} | {} | {} | {:.2e} | {} | {} | [{:.4}, {:.4}] | {:.0} | {} |",
                if v.variant.is_empty() {
                    "–"
                } else {
                    &v.variant
                },
                v.backend,
                v.device,
                v.cells,
                v.outcome.outer,
                v.outcome.linear_iterations.iter().sum::<usize>(),
                v.outcome.residuals.last().copied().unwrap_or(f64::NAN),
                v.mms
                    .as_ref()
                    .map_or("–".into(), |m| format!("{:.3e}", m.l2)),
                v.mms
                    .as_ref()
                    .map_or("–".into(), |m| format!("{:.3e}", m.linf)),
                v.range[0],
                v.range[1],
                v.ns_per_cell_iteration,
                match (v.bitwise_equal, v.diff_to_first) {
                    (Some(true), _) => "bitwise equal".into(),
                    (Some(false), Some(d)) => format!("{d:.2e}"),
                    _ => "–".into(),
                }
            );
        }
        let _ = writeln!(
            s,
            "\nTime per phase (seconds; R6), per solve: gradient / limiter / faces / assembly / linear solver / multigrid setup\n"
        );
        for v in &self.solves {
            let t = &v.outcome.timings;
            let _ = writeln!(
                s,
                "- {} {} {} cells: {:.3} / {:.3} / {:.3} / {:.3} / {:.3} / {:.3}",
                v.variant,
                v.backend,
                v.cells,
                t.gradient,
                t.limiter,
                t.faces,
                t.assemble,
                t.linear,
                t.amg_setup
            );
        }
        if !self.orders.is_empty() {
            let _ = writeln!(
                s,
                "\nObserved order of accuracy (L2, successive levels): {}",
                self.orders
                    .iter()
                    .map(|o| format!("{o:.3}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
        if self.failures.is_empty() {
            let _ = writeln!(s, "\nAll checks passed.");
        } else {
            let _ = writeln!(s, "\n**Checks failed:**\n");
            for f in &self.failures {
                let _ = writeln!(s, "- {f}");
            }
        }
        s
    }
}

const SVG_HEAD: &str = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 440" width="640" height="440" font-family="Helvetica, Arial, sans-serif" font-size="12"><rect width="640" height="440" fill="#ffffff"/>"##;
const COLOURS: [&str; 5] = ["#2a62c9", "#d1495b", "#2e933c", "#e09f3e", "#6c4ab6"];

/// Log-log plot of the L2 error over h, with a slope-2 reference line.
fn convergence_svg(name: &str, r: &[MmsResult], orders: &[f64]) -> String {
    let (x0, x1, y0, y1) = (70.0, 600.0, 390.0, 50.0);
    let lh: Vec<f64> = r.iter().map(|m| m.h.log10()).collect();
    let le: Vec<f64> = r.iter().map(|m| m.l2.log10()).collect();
    let (hmin, hmax) = (
        lh.iter().copied().fold(f64::INFINITY, f64::min).floor(),
        lh.iter().copied().fold(f64::NEG_INFINITY, f64::max).ceil(),
    );
    let (emin, emax) = (
        le.iter().copied().fold(f64::INFINITY, f64::min).floor(),
        le.iter().copied().fold(f64::NEG_INFINITY, f64::max).ceil(),
    );
    let px = |v: f64| x0 + (v - hmin) / (hmax - hmin).max(1e-9) * (x1 - x0);
    let py = |v: f64| y0 - (v - emin) / (emax - emin).max(1e-9) * (y0 - y1);
    let mut s = String::from(SVG_HEAD);
    let _ = write!(
        s,
        r##"<text x="320" y="28" text-anchor="middle" font-size="15">{name}: L2 error of φ over mesh size h</text><rect x="{x0}" y="{y1}" width="{}" height="{}" fill="none" stroke="#444"/>"##,
        x1 - x0,
        y0 - y1
    );
    let mut d = hmin;
    while d <= hmax + 1e-9 {
        let _ = write!(
            s,
            r##"<text x="{:.1}" y="{}" text-anchor="middle">1e{}</text>"##,
            px(d),
            y0 + 18.0,
            d as i32
        );
        d += 1.0;
    }
    let mut d = emin;
    while d <= emax + 1e-9 {
        let _ = write!(
            s,
            r##"<line x1="{x0}" x2="{x1}" y1="{y:.1}" y2="{y:.1}" stroke="#ddd"/><text x="{}" y="{:.1}" text-anchor="end">1e{}</text>"##,
            x0 - 6.0,
            py(d) + 4.0,
            d as i32,
            y = py(d)
        );
        d += 1.0;
    }
    // slope-2 reference through the finest point
    if let (Some(&h), Some(&e)) = (lh.last(), le.last()) {
        let ha = lh[0];
        let ea = e + 2.0 * (ha - h);
        let _ = write!(
            s,
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#999" stroke-dasharray="6 4"/><text x="{:.1}" y="{:.1}" fill="#777">slope 2</text>"##,
            px(ha),
            py(ea) - 14.0,
            px(h),
            py(e) - 14.0,
            px(ha) + 8.0,
            py(ea) - 20.0
        );
    }
    let pts: Vec<String> = lh
        .iter()
        .zip(&le)
        .map(|(h, e)| format!("{:.1},{:.1}", px(*h), py(*e)))
        .collect();
    let _ = write!(
        s,
        r##"<polyline points="{}" fill="none" stroke="{c}" stroke-width="2"/>"##,
        pts.join(" "),
        c = COLOURS[0]
    );
    for (i, (h, e)) in lh.iter().zip(&le).enumerate() {
        let _ = write!(
            s,
            r##"<circle cx="{:.1}" cy="{:.1}" r="4" fill="{}"/>"##,
            px(*h),
            py(*e),
            COLOURS[0]
        );
        if i > 0 {
            let _ = write!(
                s,
                r##"<text x="{:.1}" y="{:.1}" fill="#333">{:.2}</text>"##,
                (px(*h) + px(lh[i - 1])) / 2.0,
                (py(*e) + py(le[i - 1])) / 2.0 + 18.0,
                orders[i - 1]
            );
        }
    }
    let _ = write!(
        s,
        r##"<text x="320" y="428" text-anchor="middle">h (m); numbers: observed order between levels; {} to {} cells</text></svg>"##,
        r.first().map_or(0, |m| m.cells),
        r.last().map_or(0, |m| m.cells)
    );
    s
}

/// φ along the sampling line, one curve per variant.
fn profile_svg(name: &str, p: &[(String, Vec<(f64, f64)>)]) -> String {
    let (x0, x1, y0, y1) = (70.0, 600.0, 390.0, 50.0);
    let smax = p
        .iter()
        .flat_map(|(_, v)| v.iter().map(|x| x.0))
        .fold(0.0f64, f64::max)
        .max(1e-12);
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for (_, v) in p {
        for &(_, y) in v {
            lo = lo.min(y);
            hi = hi.max(y);
        }
    }
    let pad = 0.05 * (hi - lo).max(1e-12);
    let (lo, hi) = (lo - pad, hi + pad);
    let px = |s: f64| x0 + s / smax * (x1 - x0);
    let py = |v: f64| y0 - (v - lo) / (hi - lo) * (y0 - y1);
    let mut s = String::from(SVG_HEAD);
    let _ = write!(
        s,
        r##"<text x="320" y="28" text-anchor="middle" font-size="15">{name}: φ along the sampling line</text><rect x="{x0}" y="{y1}" width="{}" height="{}" fill="none" stroke="#444"/>"##,
        x1 - x0,
        y0 - y1
    );
    for v in [0.0, 1.0] {
        if v > lo && v < hi {
            let _ = write!(
                s,
                r##"<line x1="{x0}" x2="{x1}" y1="{y:.1}" y2="{y:.1}" stroke="#bbb" stroke-dasharray="4 3"/><text x="{}" y="{:.1}" text-anchor="end">{v}</text>"##,
                x0 - 6.0,
                py(v) + 4.0,
                y = py(v)
            );
        }
    }
    for (i, (label, v)) in p.iter().enumerate() {
        let c = COLOURS[i % COLOURS.len()];
        let pts: Vec<String> = v
            .iter()
            .map(|(a, b)| format!("{:.1},{:.1}", px(*a), py(*b)))
            .collect();
        let _ = write!(
            s,
            r##"<polyline points="{}" fill="none" stroke="{c}" stroke-width="2"/><rect x="{}" y="{}" width="14" height="4" fill="{c}"/><text x="{}" y="{}">{label}</text>"##,
            pts.join(" "),
            x0 + 12.0,
            y1 + 14.0 + 18.0 * i as f64,
            x0 + 32.0,
            y1 + 19.0 + 18.0 * i as f64
        );
    }
    let _ = write!(
        s,
        r##"<text x="320" y="420" text-anchor="middle">distance along the line (m)</text></svg>"##
    );
    s
}

/// Case name and output directory of a case file.
pub fn output_dir(case: &ScalarCase, path: &Path, out: Option<PathBuf>) -> (String, PathBuf) {
    let name = path
        .file_stem()
        .map_or("scalar".into(), |s| s.to_string_lossy().into_owned());
    let dir = out
        .or_else(|| case.output.directory.clone())
        .unwrap_or_else(|| PathBuf::from("output").join(&name));
    (name, dir)
}

/// Fails with the list of failed checks.
pub fn check(rep: &ScalarReport) -> Result<()> {
    if rep.failures.is_empty() {
        Ok(())
    } else {
        bail!("checks failed:\n  {}", rep.failures.join("\n  "))
    }
}
