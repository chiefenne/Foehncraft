//! The engine API under every front end (D86): load and validate a case,
//! check it and estimate its cost, run it, stop it, read its results. The
//! CLI (`foehn flow`), and later the Python package (D61), the MCP server
//! (L18) and the desktop app (L7), are thin layers on these functions;
//! nothing is reachable only from the CLI's `main.rs`.
//!
//! The run contract (D62, D64, D69, D79–D85):
//!
//! - [`load`] reads a case with `--set key=value` overrides and the file
//!   confinement of [`Roots`] (D82);
//! - [`check_file`] gives the diagnostics of the case and its boundary
//!   conditions (codes of [`crate::diag`], D78) and the estimate of
//!   [`Estimate`] (D81);
//! - [`run`] runs it and writes, next to the outputs of the solver,
//!   `result.json` (the result summary, D62) and `run.jsonl` (the run
//!   record, D69), both with the provenance of [`Provenance`] (D85); it
//!   ends with a [`Status`] whose exit code is that of D64;
//! - [`Interrupts`] carries the requests of signals or of another thread:
//!   write a restart file and go on, or stop (D64); `[run.limits]` stops a
//!   run at a wall time, step count or output volume (D83).
//!
//! `result.json` holds no timings: two runs of the same case with the same
//! binary on the same hardware give the same bytes, and its `files` list
//! holds the size and SHA-256 hash of every other output file except the
//! four that record wall-clock times (`flow.log`, `flow.md`, `flow.json`,
//! `run.jsonl`), so comparing two `result.json` files compares every
//! result of two runs (D84).

use crate::case::{FlowCase, FlowReport, LATEST, Limits};
use crate::diag::{self, Diagnostic, Severity};
use anyhow::{Context, Result};
use foehncraft_mesh::Mesh;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Instant;

/// Version of the machine formats of this module (D80): `result.json`,
/// `run.jsonl`, the check report and the error object. Within a major
/// version only additive changes.
pub const FORMAT_VERSION: &str = "1.1";

/// Output files that record wall-clock times; every other file of a run is
/// byte-identical between two runs of the same case (D84).
pub const MEASURED_FILES: &[&str] = &[
    "flow.log",
    "flow.md",
    "flow.json",
    "run.jsonl",
    "result.json",
];

// ---------------------------------------------------------------------------
// interrupts and stops (D64, D83)

/// A request to stop, by signal or from the API.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, schemars::JsonSchema)]
pub enum StopSignal {
    /// SIGTERM (a scheduler ends the job).
    #[serde(rename = "SIGTERM")]
    Term = 1,
    /// SIGINT (Ctrl-C).
    #[serde(rename = "SIGINT")]
    Int = 2,
    /// SIGUSR2.
    #[serde(rename = "SIGUSR2")]
    Usr2 = 3,
    /// [`Interrupts::request_stop`] from another thread.
    #[serde(rename = "api")]
    Api = 4,
}

impl StopSignal {
    /// The value a signal handler stores ([`Interrupts::stop_flag`]).
    pub const fn value(self) -> usize {
        self as usize
    }

    fn from_value(v: usize) -> Option<Self> {
        match v {
            1 => Some(Self::Term),
            2 => Some(Self::Int),
            3 => Some(Self::Usr2),
            4 => Some(Self::Api),
            _ => None,
        }
    }

    /// Name (`SIGTERM` …).
    pub fn name(self) -> &'static str {
        match self {
            Self::Term => "SIGTERM",
            Self::Int => "SIGINT",
            Self::Usr2 => "SIGUSR2",
            Self::Api => "a stop request",
        }
    }
}

/// Requests to a running case. The front end connects signals to the
/// flags (`foehn` with `signal-hook`: SIGUSR1 sets [`restart_flag`], SIGTERM,
/// SIGINT and SIGUSR2 store their [`StopSignal::value`] in [`stop_flag`]);
/// the run looks at them after every outer iteration or time step.
///
/// [`restart_flag`]: Self::restart_flag
/// [`stop_flag`]: Self::stop_flag
#[derive(Debug, Clone, Default)]
pub struct Interrupts {
    restart: Arc<AtomicBool>,
    stop: Arc<AtomicUsize>,
}

impl Interrupts {
    /// No requests.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set to write a restart file at the end of the current iteration or
    /// step and go on (SIGUSR1).
    pub fn restart_flag(&self) -> Arc<AtomicBool> {
        self.restart.clone()
    }

    /// Holds the [`StopSignal::value`] of a stop request (0: none).
    pub fn stop_flag(&self) -> Arc<AtomicUsize> {
        self.stop.clone()
    }

    /// Asks the run to stop after the current iteration or step, with a
    /// restart file and the final output.
    pub fn request_stop(&self, s: StopSignal) {
        self.stop.store(s.value(), Ordering::SeqCst);
    }

    /// Asks the run for a restart file after the current iteration or step.
    pub fn request_restart(&self) {
        self.restart.store(true, Ordering::SeqCst);
    }

    fn take_restart(&self) -> bool {
        self.restart.swap(false, Ordering::SeqCst)
    }

    fn stop(&self) -> Option<StopSignal> {
        StopSignal::from_value(self.stop.load(Ordering::SeqCst))
    }
}

/// Why a run stopped before its end (exit code 3, D64).
#[derive(Debug, Clone, PartialEq, Serialize, schemars::JsonSchema)]
pub struct Stop {
    /// Diagnostic code (`RUN001` … `RUN005`, or `BAL001` … `BAL003`).
    pub code: &'static str,
    /// `control-file`, `signal`, `wall-time`, `steps`, `iterations`,
    /// `output`, `balance`, `diagnosis`, `study-budget`.
    pub kind: &'static str,
    /// The signal, for `signal`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub signal: Option<StopSignal>,
    /// For the log.
    pub message: String,
}

/// What the loops of a run need from the engine: the requests, the limits
/// and the run record.
pub struct Steer {
    interrupts: Interrupts,
    start: Instant,
    out: PathBuf,
    pub(crate) record: Record,
    /// The stop decided by [`poll`](Self::poll).
    pub(crate) stop: Option<Stop>,
    /// A stop that the run itself asked for (`[balance] action = "stop"`),
    /// returned by the next [`poll`](Self::poll).
    pending: Option<Stop>,
    /// The build and host, for the banner of the log.
    pub(crate) build: Option<BuildInfo>,
    /// The spend of this run against the study budgets (D198 (c)).
    spend: SpendMeter,
    /// Lines for the log after its banner (the budgets, D198).
    pub(crate) notes: Vec<String>,
    /// The live control port of the run (L19), when it has one.
    pub(crate) live: Option<crate::live::Live>,
    /// The changes made during the run, for the Events block of the log
    /// (D214).
    pub(crate) events: Vec<EventLine>,
}

/// One line of the Events block (D214): a change made during the run.
#[derive(Debug, Clone, PartialEq)]
pub struct EventLine {
    /// The code (`RUN015`, `RUN023`, ...).
    pub code: &'static str,
    /// The iteration or step.
    pub at: usize,
    /// Wall-clock seconds since the start of the run.
    pub elapsed: f64,
    /// The message.
    pub message: String,
}

/// The codes whose events go into the Events block (D214): every change of
/// a setting (schedule, control port, re-read), every stop, the diagnosis
/// and the steps of the thermal acceleration, the cycle-averaged loop and
/// the frozen-flow march.
pub const EVENT_CODES: &[&str] = &[
    diag::RUN_STOP_FILE,
    diag::RUN_SIGNAL,
    diag::RUN_LIMIT_WALL_TIME,
    diag::RUN_LIMIT_STEPS,
    diag::RUN_LIMIT_OUTPUT,
    diag::RUN_NOT_APPLIED,
    diag::RUN_REREAD_INVALID,
    diag::RUN_SCHEDULED,
    diag::RUN_DIAGNOSIS,
    diag::RUN_DIAGNOSIS_STOP,
    diag::RUN_STUDY_BUDGET,
    diag::RUN_CONTROL_STOP,
    diag::RUN_CONTROL_SET,
    diag::RUN_CONTROL_REJECTED,
    diag::RUN_REREAD,
    diag::ACC_PHASE,
    diag::CYC_PHASE,
    diag::FRZ_PHASE,
];

/// What a run has spent so far (cell-seconds, GPUs), and what the earlier
/// runs of its study spent (D198 (c)).
#[derive(Debug, Clone)]
struct SpendMeter {
    cells: usize,
    gpus: usize,
    cell_seconds: f64,
    last: Instant,
    /// The spend of the earlier runs of the study (only with a budget).
    before: Option<crate::budget::Spend>,
}

impl SpendMeter {
    fn new() -> Self {
        Self {
            cells: 0,
            gpus: 0,
            cell_seconds: 0.0,
            last: Instant::now(),
            before: None,
        }
    }

    fn tick(&mut self) {
        let now = Instant::now();
        self.cell_seconds += self.cells as f64 * (now - self.last).as_secs_f64();
        self.last = now;
    }
}

impl Steer {
    /// Requests from `interrupts`; the run record goes to `out/run.jsonl`
    /// when `record` is set.
    pub fn new(interrupts: Interrupts, out: &Path, record: bool) -> Result<Self> {
        Ok(Self {
            interrupts,
            start: Instant::now(),
            out: out.to_path_buf(),
            record: if record {
                Record::create(&out.join("run.jsonl"))?
            } else {
                Record::none()
            },
            stop: None,
            pending: None,
            build: None,
            spend: SpendMeter::new(),
            notes: Vec::new(),
            live: None,
            events: Vec::new(),
        })
    }

    /// Without requests, limits or record (tests, the plain [`crate::case::run`]).
    pub fn quiet(out: &Path) -> Self {
        Self {
            interrupts: Interrupts::new(),
            start: Instant::now(),
            out: out.to_path_buf(),
            record: Record::none(),
            stop: None,
            pending: None,
            build: None,
            spend: SpendMeter::new(),
            notes: Vec::new(),
            live: None,
            events: Vec::new(),
        }
    }

    /// The solve now runs on `cells` cells and `gpus` GPUs (the spend of a
    /// study budget counts from here, D198 (c)).
    pub(crate) fn set_cells(&mut self, cells: usize, gpus: usize) {
        self.spend.tick();
        self.spend.cells = cells;
        self.spend.gpus = gpus;
    }

    /// What this run has spent so far.
    fn spent(&mut self) -> crate::budget::Spend {
        self.spend.tick();
        crate::budget::Spend {
            cell_hours: self.spend.cell_seconds / 3600.0,
            gpu_hours: self.spend.gpus as f64 * self.start.elapsed().as_secs_f64() / 3600.0,
            output_bytes: dir_bytes(&self.out) as f64,
        }
    }

    /// The diagnosis of a steady run (D198 (d)) into the run record; with
    /// `stop` and a verdict off track the run stops at the next
    /// [`poll`](Self::poll).
    pub(crate) fn diagnosis(&mut self, d: &crate::budget::Diagnosis, stop: bool) {
        let mut row = serde_json::to_value(d).unwrap_or_default();
        row["type"] = "diagnosis".into();
        self.record.row(&row);
        self.event(
            diag::RUN_DIAGNOSIS,
            d.iteration,
            &format!("diagnosis at iteration {}: {}", d.iteration, d.advice),
        );
        if stop && d.verdict.off_track() && self.pending.is_none() {
            self.pending = Some(Stop {
                code: diag::RUN_DIAGNOSIS_STOP,
                kind: "diagnosis",
                signal: None,
                message: format!(
                    "stopped by the diagnosis at iteration {} ({:?}; [run.limits] diagnosis = \"stop\")",
                    d.iteration, d.verdict
                )
                .to_lowercase(),
            });
        }
    }

    /// Stops the run at the next [`poll`](Self::poll) as at a limit
    /// (restart file, final output, exit code 3) with `code` (kind
    /// `balance`).
    pub(crate) fn stop_for(&mut self, code: &'static str, message: String) {
        if self.pending.is_none() {
            self.pending = Some(Stop {
                code,
                kind: "balance",
                signal: None,
                message,
            });
        }
    }

    /// After an iteration or step: stop now? `done` counts the iterations
    /// or steps of this run, `stop_file` is the `.stop` control file.
    pub(crate) fn poll(
        &mut self,
        limits: &Limits,
        done: usize,
        transient: bool,
        stop_file: bool,
    ) -> Option<Stop> {
        let stop = |code, kind, message: String| Stop {
            code,
            kind,
            signal: None,
            message,
        };
        let s = if let Some(p) = self.pending.take() {
            Some(p)
        } else if stop_file {
            Some(stop(
                diag::RUN_STOP_FILE,
                "control-file",
                "stopped by .stop".into(),
            ))
        } else if self.live.as_ref().is_some_and(|l| l.stop_requested()) {
            Some(stop(
                diag::RUN_CONTROL_STOP,
                "control-port",
                "stopped through the control port (run.stop)".into(),
            ))
        } else if let Some(sig) = self.interrupts.stop() {
            Some(Stop {
                code: diag::RUN_SIGNAL,
                kind: "signal",
                signal: Some(sig),
                message: format!("stopped by {}", sig.name()),
            })
        } else if let Some(t) = limits.max_wall_time
            && self.start.elapsed().as_secs_f64() >= t
        {
            Some(stop(
                diag::RUN_LIMIT_WALL_TIME,
                "wall-time",
                format!("stopped at the wall-time limit of {t} s ([run.limits] max_wall_time)"),
            ))
        } else if transient && limits.max_steps.is_some_and(|n| done >= n) {
            Some(stop(
                diag::RUN_LIMIT_STEPS,
                "steps",
                format!("stopped at the limit of {done} time steps ([run.limits] max_steps)"),
            ))
        } else if !transient && limits.max_iterations.is_some_and(|n| done >= n) {
            Some(stop(
                diag::RUN_LIMIT_STEPS,
                "iterations",
                format!(
                    "stopped at the limit of {done} outer iterations ([run.limits] max_iterations)"
                ),
            ))
        } else if let Some(b) = limits.max_output
            && dir_bytes(&self.out) as f64 >= b
        {
            Some(stop(
                diag::RUN_LIMIT_OUTPUT,
                "output",
                format!("stopped at the output limit of {b:e} bytes ([run.limits] max_output)"),
            ))
        } else if let Some(before) = self.spend.before.clone()
            && limits.has_study_budget()
            && let Some((k, b, used)) = {
                let now = self.spent();
                [
                    (
                        "study_cell_hours",
                        limits.study_cell_hours,
                        before.cell_hours + now.cell_hours,
                    ),
                    (
                        "study_gpu_hours",
                        limits.study_gpu_hours,
                        before.gpu_hours + now.gpu_hours,
                    ),
                    (
                        "study_output",
                        limits.study_output,
                        before.output_bytes + now.output_bytes,
                    ),
                ]
                .into_iter()
                .find_map(|(k, b, u)| b.filter(|&b| u >= b).map(|b| (k, b, u)))
            }
        {
            Some(stop(
                diag::RUN_STUDY_BUDGET,
                "study-budget",
                format!(
                    "stopped at the study budget ([run.limits] {k} = {b}; spent {used:.4} with this run)"
                ),
            ))
        } else {
            None
        };
        if let Some(s) = &s {
            self.record.row(&serde_json::json!({
                "type": "event",
                "code": s.code,
                "severity": "note",
                "message": s.message,
                "at": done,
            }));
            self.remember(s.code, done, &s.message);
            self.stop = Some(s.clone());
        }
        s
    }

    /// The settings asked for through the control port since the last call
    /// (L19).
    pub(crate) fn take_settings(&self) -> Vec<crate::live::Pending> {
        self.live
            .as_ref()
            .map(|l| l.take_settings())
            .unwrap_or_default()
    }

    /// What the run did with a request of the control port, for its status.
    pub(crate) fn live_note(&self, entry: serde_json::Value) {
        if let Some(l) = &self.live {
            l.note(entry);
        }
    }

    /// A restart file was asked for (SIGUSR1) since the last call.
    pub(crate) fn restart_requested(&self) -> bool {
        self.interrupts.take_restart()
    }

    /// A warning or note with a code, into the run record (and the Events
    /// block of the log when it is a change, D214).
    pub(crate) fn event(&mut self, code: &'static str, at: usize, message: &str) {
        let sev = diag::info(code).map_or(Severity::Note, |i| i.severity);
        self.record.row(&serde_json::json!({
            "type": "event",
            "code": code,
            "severity": sev,
            "message": message,
            "at": at,
        }));
        self.remember(code, at, message);
    }

    /// Keeps an event for the Events block of the log (D214) when its code
    /// is one of [`EVENT_CODES`].
    pub(crate) fn remember(&mut self, code: &'static str, at: usize, message: &str) {
        if EVENT_CODES.contains(&code) {
            self.events.push(EventLine {
                code,
                at,
                elapsed: self.start.elapsed().as_secs_f64(),
                message: message.to_string(),
            });
        }
    }
}

/// Bytes of the regular files in `dir` (not recursive).
fn dir_bytes(dir: &Path) -> u64 {
    std::fs::read_dir(dir)
        .map(|rd| {
            rd.flatten()
                .filter_map(|e| e.metadata().ok())
                .filter(|m| m.is_file())
                .map(|m| m.len())
                .sum()
        })
        .unwrap_or(0)
}

// ---------------------------------------------------------------------------
// run record (D69)

/// A function that sees every row of the run record as it is written
/// (header, iterations, steps, balances, events, end): the hook of the
/// Python package for monitors and stop criteria at iteration or step
/// level (D61). It runs on the thread of the run, between iterations or
/// steps; it reads, it never changes the state, so results do not depend
/// on it (D31).
#[derive(Clone)]
pub struct Observer(Arc<std::sync::Mutex<ObserverFn>>);

/// The function of an [`Observer`].
pub type ObserverFn = dyn FnMut(&serde_json::Value) + Send;

impl Observer {
    /// Wraps `f`.
    pub fn new(f: impl FnMut(&serde_json::Value) + Send + 'static) -> Self {
        Self(Arc::new(std::sync::Mutex::new(f)))
    }

    fn call(&self, v: &serde_json::Value) {
        if let Ok(mut f) = self.0.lock() {
            f(v);
        }
    }
}

impl std::fmt::Debug for Observer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Observer")
    }
}

/// The run record: JSON lines, buffered and flushed with the log lines.
pub(crate) struct Record {
    file: Option<std::io::BufWriter<std::fs::File>>,
    observer: Option<Observer>,
    /// The state for the live control port (L19).
    live: Option<crate::live::Live>,
}

impl Record {
    fn create(path: &Path) -> Result<Self> {
        let f =
            std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
        Ok(Self {
            file: Some(std::io::BufWriter::new(f)),
            observer: None,
            live: None,
        })
    }

    fn none() -> Self {
        Self {
            file: None,
            observer: None,
            live: None,
        }
    }

    /// Appends one row.
    pub(crate) fn row(&mut self, v: &impl Serialize) {
        if (self.observer.is_some() || self.live.is_some())
            && let Ok(j) = serde_json::to_value(v)
        {
            if let Some(o) = &self.observer {
                o.call(&j);
            }
            if let Some(l) = &self.live {
                l.observe(&j);
            }
        }
        if let Some(f) = self.file.as_mut()
            && let Ok(s) = serde_json::to_string(v)
        {
            let _ = writeln!(f, "{s}");
        }
    }

    /// Writes the buffered rows.
    pub(crate) fn flush(&mut self) {
        if let Some(f) = self.file.as_mut() {
            let _ = f.flush();
        }
    }

    /// Hands the rows to the live control port (L19).
    pub(crate) fn live(&mut self, l: crate::live::Live) {
        self.live = Some(l);
    }

    /// True when rows are written or observed.
    pub(crate) fn on(&self) -> bool {
        self.file.is_some() || self.observer.is_some()
    }
}

// ---------------------------------------------------------------------------
// provenance (D85)

/// What the front end knows about the build and the host.
#[derive(Debug, Clone, Default, Serialize)]
pub struct BuildInfo {
    /// Git commit of the build (`unknown` outside a Git checkout; `-dirty`
    /// with uncommitted changes).
    pub commit: String,
    /// Cargo features of the build.
    pub features: Vec<String>,
    /// Operating system and architecture.
    pub platform: String,
    /// CPU model.
    pub cpu: String,
    /// GPU of the `cuda` backend, when the front end knows it.
    pub gpu: Option<String>,
}

/// Where a result comes from (D85): enough to tell whether two results
/// are comparable.
#[derive(Debug, Clone, Serialize, schemars::JsonSchema)]
pub struct Provenance {
    /// Foehncraft version.
    pub version: String,
    /// Git commit.
    pub commit: String,
    /// Build features.
    pub features: Vec<String>,
    /// Flow backends of the build (`cpu`, `cuda`).
    pub backends: Vec<String>,
    /// Backend of the run.
    pub backend: String,
    /// Operating system and architecture.
    pub platform: String,
    /// CPU model.
    pub cpu: String,
    /// GPU of the run, if one.
    pub gpu: Option<String>,
    /// MPI ranks.
    pub ranks: usize,
    /// GPUs.
    pub gpus: usize,
    /// The strict reproducible mode of D31 (2) and D93 (partition and rank
    /// independence); not available yet. The always-on part of D31 holds
    /// for every run: bitwise the same results for any thread count and on
    /// CPU and GPU.
    pub reproducible_mode: bool,
    /// SHA-256 of the resolved case without the keys that do not change
    /// the results (`run.threads`, `run.backend`, `run.ptx`,
    /// `output.directory`, the `[log]` table).
    pub case_hash: String,
    /// SHA-256 of the resolved property program (D94 (3)).
    pub properties_hash: String,
    /// SHA-256 of the built-in material database.
    pub materials_database: String,
}

fn provenance(case: &FlowCase, build: &BuildInfo, gpu: Option<String>) -> Provenance {
    Provenance {
        version: foehncraft_core::VERSION.into(),
        commit: build.commit.clone(),
        features: build.features.clone(),
        backends: crate::device::available()
            .iter()
            .map(|s| s.to_string())
            .collect(),
        backend: case.run.backend.clone(),
        platform: build.platform.clone(),
        cpu: build.cpu.clone(),
        gpus: usize::from(gpu.is_some()),
        gpu,
        ranks: 1,
        reproducible_mode: false,
        case_hash: case_hash(case),
        properties_hash: properties_hash(case),
        materials_database: crate::materials::database_hash(),
    }
}

/// SHA-256 (hex) of `bytes`.
pub fn sha256(bytes: &[u8]) -> String {
    hex(&Sha256::digest(bytes))
}

fn hex(b: &[u8]) -> String {
    use std::fmt::Write as _;
    b.iter()
        .fold(String::with_capacity(2 * b.len()), |mut s, x| {
            let _ = write!(s, "{x:02x}");
            s
        })
}

/// SHA-256 of the resolved case, without the keys that only say where it
/// runs (D31: the results do not depend on them).
pub fn case_hash(case: &FlowCase) -> String {
    let mut c = case.clone();
    c.run.threads = None;
    c.run.backend = "cpu".into();
    c.run.ptx = crate::case::RunSection::default().ptx;
    c.output.directory = None;
    // how the log looks does not change a result
    c.log = crate::console::LogSection::default();
    let mut text = toml::to_string_pretty(&c).unwrap_or_default();
    // the contents of the profile tables (the case names only their paths)
    for (name, b) in &case.boundary {
        if let Some(t) = b.table() {
            text.push_str(&format!("\n# table of {name}: {}", t.sha256));
        }
    }
    // the scalar sub-tables of the boundaries (not serialised with the
    // boundary conditions)
    if !case.boundary_scalar.is_empty() {
        text.push_str(&format!("\n# boundary scalars: {:?}", case.boundary_scalar));
    }
    for (name, t) in &case.tables {
        if let Some(h) = t.sha256() {
            text.push_str(&format!("\n# table {name}: {h}"));
        }
    }
    if let Ok(p) = case.properties() {
        text.push('\n');
        text.push_str(&p.text());
    }
    sha256(text.as_bytes())
}

/// SHA-256 of the resolved property program of a case (D94 (3)).
pub fn properties_hash(case: &FlowCase) -> String {
    case.properties()
        .map(|p| sha256(p.text().as_bytes()))
        .unwrap_or_default()
}

/// SHA-256 of a mesh: points, faces, cells, boundaries and the global
/// identities, as stored (D85).
pub fn mesh_hash(m: &Mesh) -> String {
    let mut h = Sha256::new();
    let put_u32 = |h: &mut Sha256, v: &[u32]| {
        h.update((v.len() as u64).to_le_bytes());
        for x in v {
            h.update(x.to_le_bytes());
        }
    };
    h.update((m.points.len() as u64).to_le_bytes());
    for p in &m.points {
        for x in p {
            h.update(x.to_bits().to_le_bytes());
        }
    }
    put_u32(&mut h, &m.face_ptr);
    put_u32(&mut h, &m.face_nodes);
    put_u32(&mut h, &m.owner);
    put_u32(&mut h, &m.neighbour);
    put_u32(&mut h, &m.cell_zone);
    for ids in [&m.cell_id, &m.face_id, &m.point_id] {
        h.update((ids.len() as u64).to_le_bytes());
        for x in ids.iter() {
            h.update(x.to_le_bytes());
        }
    }
    for p in &m.boundaries {
        h.update(p.name.as_bytes());
        h.update([0]);
        h.update((p.start as u64).to_le_bytes());
        h.update((p.len as u64).to_le_bytes());
    }
    hex(&h.finalize())
}

/// One output file of `result.json`.
#[derive(Debug, Clone, Serialize, PartialEq, Eq, schemars::JsonSchema)]
pub struct FileHash {
    /// Name in the output directory.
    pub name: String,
    /// Size.
    pub bytes: u64,
    /// SHA-256 (hex).
    pub sha256: String,
}

/// Size and hash of every output file of `dir` except the measured ones
/// and the control files (D40, L19), by name.
pub fn file_hashes(dir: &Path) -> Result<Vec<FileHash>> {
    let mut out = Vec::new();
    let Ok(rd) = std::fs::read_dir(dir) else {
        return Ok(out);
    };
    for e in rd.flatten() {
        let name = e.file_name().to_string_lossy().into_owned();
        // the files of the live control port exist only while the run goes
        // (L19) and hold a token: never part of the result
        if name.starts_with('.')
            || name.ends_with(".tmp")
            || MEASURED_FILES.contains(&name.as_str())
            || name == crate::live::CONTROL_FILE
            || name == crate::live::TOKEN_FILE
        {
            continue;
        }
        if !e.file_type().is_ok_and(|t| t.is_file()) {
            continue;
        }
        let bytes = std::fs::read(e.path())?;
        out.push(FileHash {
            name,
            bytes: bytes.len() as u64,
            sha256: sha256(&bytes),
        });
    }
    out.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(out)
}

/// The date of the run: `SOURCE_DATE_EPOCH` when set (D84), else now
/// (ISO 8601, UTC).
pub fn run_date() -> String {
    let secs = std::env::var("SOURCE_DATE_EPOCH")
        .ok()
        .and_then(|v| v.trim().parse::<i64>().ok())
        .unwrap_or_else(|| {
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs() as i64)
                .unwrap_or(0)
        });
    let [y, mo, d, h, mi, s] = foehncraft_core::date::civil(secs);
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{mi:02}:{s:02}Z")
}

// ---------------------------------------------------------------------------
// file confinement (D82)

/// Read and write roots (`--read-root`, `--write-root`): with any given,
/// a path outside them is refused. Relative paths of a case file resolve
/// against its directory; writes go to the output directory only.
#[derive(Debug, Clone, Default)]
pub struct Roots {
    /// Directories that may be read.
    pub read: Vec<PathBuf>,
    /// Directories that may be written.
    pub write: Vec<PathBuf>,
}

/// The absolute form of `p`, resolving symbolic links of the part that
/// exists (an output directory may not exist yet).
fn absolute(p: &Path) -> PathBuf {
    let p = if p.is_absolute() {
        p.to_path_buf()
    } else {
        std::env::current_dir()
            .map(|d| d.join(p))
            .unwrap_or_else(|_| p.to_path_buf())
    };
    let mut existing = p.as_path();
    let mut rest = Vec::new();
    while !existing.exists() {
        let (Some(parent), Some(name)) = (existing.parent(), existing.file_name()) else {
            break;
        };
        rest.push(name.to_os_string());
        existing = parent;
    }
    let mut out = existing
        .canonicalize()
        .unwrap_or_else(|_| existing.to_path_buf());
    for r in rest.iter().rev() {
        out.push(r);
    }
    // lexical clean-up of `..` in the part that does not exist
    let mut clean = PathBuf::new();
    for c in out.components() {
        match c {
            std::path::Component::ParentDir => {
                clean.pop();
            }
            std::path::Component::CurDir => {}
            other => clean.push(other),
        }
    }
    clean
}

impl Roots {
    fn inside(roots: &[PathBuf], p: &Path) -> bool {
        let p = absolute(p);
        roots.iter().any(|r| p.starts_with(absolute(r)))
    }

    /// Fails unless `p` may be read.
    pub fn check_read(&self, p: &Path, what: &str) -> Result<()> {
        if self.read.is_empty() || Self::inside(&self.read, p) {
            return Ok(());
        }
        Err(diag::err(Diagnostic::new(
            diag::CASE_OUTSIDE_ROOT,
            format!(
                "{what} {} is outside the read roots ({})",
                p.display(),
                show_roots(&self.read)
            ),
        )))
    }

    /// Fails unless the directory `p` may be written.
    pub fn check_write(&self, p: &Path) -> Result<()> {
        if self.write.is_empty() || Self::inside(&self.write, p) {
            return Ok(());
        }
        Err(diag::err(Diagnostic::new(
            diag::CASE_OUTSIDE_ROOT,
            format!(
                "the output directory {} is outside the write roots ({})",
                p.display(),
                show_roots(&self.write)
            ),
        )))
    }
}

fn show_roots(r: &[PathBuf]) -> String {
    r.iter()
        .map(|p| p.display().to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

// ---------------------------------------------------------------------------
// load, restart latest

/// Reads a case file with `--set` overrides; the case file, the mesh file
/// and the restart file must lie inside the read roots.
pub fn load(path: &Path, sets: &[String], roots: &Roots) -> Result<FlowCase> {
    roots.check_read(path, "the case file")?;
    let mut case = FlowCase::from_file_unloaded(path, sets)?;
    check_paths(&case, roots)?;
    case.load_tables()?;
    Ok(case)
}

/// The input files of a case against the read roots.
pub fn check_paths(case: &FlowCase, roots: &Roots) -> Result<()> {
    if let Some(f) = &case.mesh.file {
        roots.check_read(f, "the mesh file")?;
    }
    if let Some(r) = &case.run.restart
        && r.as_os_str() != LATEST
    {
        roots.check_read(r, "the restart file")?;
    }
    for b in case.boundary.values() {
        if let Some(t) = b.table_path() {
            roots.check_read(t, "the profile table")?;
        }
    }
    Ok(())
}

/// The newest restart file of `dir`: the one with the largest step or
/// iteration (read from the files, D75); `restart.h5` wins a tie.
pub fn latest_restart(dir: &Path) -> Result<Option<PathBuf>> {
    let mut best: Option<(u64, bool, PathBuf)> = None;
    let Ok(rd) = std::fs::read_dir(dir) else {
        return Ok(None);
    };
    let mut names: Vec<PathBuf> = rd
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.file_name().and_then(|n| n.to_str()).is_some_and(|n| {
                n == "restart.h5" || (n.starts_with("restart_") && n.ends_with(".h5"))
            })
        })
        .collect();
    names.sort();
    for p in names {
        let Ok((_, n)) = crate::restart::progress(&p) else {
            continue;
        };
        let main = p.file_name().is_some_and(|f| f == "restart.h5");
        let better = match &best {
            None => true,
            Some((m, was_main, _)) => n > *m || (n == *m && main && !was_main),
        };
        if better {
            best = Some((n, main, p));
        }
    }
    Ok(best.map(|b| b.2))
}

/// Replaces `run.restart = "latest"` by the newest restart file of `dir`.
pub fn resolve_latest(case: &mut FlowCase, dir: &Path) -> Result<()> {
    if case.run.restart.as_deref() != Some(Path::new(LATEST)) {
        return Ok(());
    }
    match latest_restart(dir)? {
        Some(p) => {
            case.run.restart = Some(p);
            Ok(())
        }
        None => Err(diag::err(Diagnostic::new(
            diag::RUN_NO_RESTART,
            format!(
                "--restart latest: no restart file (restart.h5, restart_<N>.h5) in {}",
                dir.display()
            ),
        ))),
    }
}

// ---------------------------------------------------------------------------
// check and estimate (D78, D81)

/// One output file of the estimate.
#[derive(Debug, Clone, Serialize, schemars::JsonSchema)]
pub struct OutputEstimate {
    /// File name (`fields_<step>.vtkhdf` for a series).
    pub name: String,
    /// Number of files.
    pub count: usize,
    /// Bytes of each.
    pub bytes: f64,
}

/// What a run will cost (D81), before it starts.
#[derive(Debug, Clone, Serialize, schemars::JsonSchema)]
pub struct Estimate {
    /// Cells of the largest level.
    pub cells: usize,
    /// Faces of the largest level.
    pub faces: usize,
    /// Interior faces of the largest level.
    pub interior_faces: usize,
    /// Cells of every level of a sweep.
    pub level_cells: Vec<usize>,
    /// MPI ranks.
    pub ranks: usize,
    /// GPUs.
    pub gpus: usize,
    /// Outer iterations (steady, an upper bound: `solver.max_iterations`)
    /// or time steps (transient; with an adaptive step from the first one)
    /// of every level.
    pub steps: Vec<usize>,
    /// The step counts are upper bounds or guesses (steady convergence,
    /// adaptive time steps).
    pub steps_approximate: bool,
    /// Host memory (bytes).
    pub host_memory: f64,
    /// Device memory (bytes; 0 on the CPU).
    pub device_memory: f64,
    /// Output volume (bytes).
    pub output_bytes: f64,
    /// Output files.
    pub outputs: Vec<OutputEstimate>,
    /// Wall time (s), when a benchmark of this hardware exists.
    pub wall_time: Option<f64>,
    /// Where the wall time comes from, or why there is none.
    pub wall_time_basis: String,
    /// Cell-hours (cells times wall hours, summed over the levels), with
    /// the wall time (D198 (c)).
    pub cell_hours: Option<f64>,
    /// GPU-hours, with the wall time (0 on the CPU).
    pub gpu_hours: Option<f64>,
    /// Steady: the iteration budget, its source and the diagnosis
    /// (D198 (d)).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub iteration_budget: Option<crate::budget::IterationBudget>,
}

/// Host memory of a run: a base and bytes per interior and per boundary
/// face, fitted to the peak resident memory of eleven cases of the
/// catalogue (2D and 3D, hexahedra to polyhedra, steady and transient, 2
/// threads; 128 to 28 227 cells): within 6 % on each
/// (`tools/plot/estimate_chart.py`, gallery 18). An interior face carries two
/// off-diagonal entries of each matrix (SELL storage), its fluxes and its
/// share of the multigrid hierarchy; a boundary face only its geometry and
/// boundary values, which matters for the one-cell-thick 2D meshes, whose
/// front and back faces are half of all faces.
const HOST_BYTES_BASE: f64 = 14e6;
const HOST_BYTES_PER_INTERIOR_FACE: f64 = 1960.0;
const HOST_BYTES_PER_BOUNDARY_FACE: f64 = 710.0;
/// Device memory of the GPU flow loop: the host model without the base and
/// the host copies of the mesh, about three quarters of it (an assumption,
/// not measured on a GPU yet).
const DEVICE_SHARE: f64 = 0.75;

/// Measured time per cell and iteration or step (ns): hardware (a part of
/// the device or CPU name), backend, transient, value. From
/// `docs/benchmarks/stage5_gpu_t4.md`.
const BENCHMARKS: &[(&str, &str, bool, f64)] = &[
    ("Tesla T4", "cuda", false, 722.0),
    ("Tesla T4", "cuda", true, 1049.0),
];

/// Output files of a case on `mesh` and their sizes, for `steps`
/// iterations or steps.
pub(crate) fn outputs(case: &FlowCase, mesh: &Mesh, steps: usize) -> Vec<OutputEstimate> {
    use crate::bc::BoundaryCondition;
    let nc = mesh.n_cells;
    let o = &case.output;
    let r = &case.run;
    let mut v = Vec::new();
    let mut add = |name: &str, count: usize, bytes: f64| {
        if count > 0 {
            v.push(OutputEstimate {
                name: name.into(),
                count,
                bytes,
            });
        }
    };
    let transient = case.time.is_some();
    let rows = steps / o.history_every + 1;
    if transient {
        add("history.csv", 1, rows as f64 * 120.0);
        add("balance.csv", 1, rows as f64 * 240.0);
    } else {
        add("residuals.csv", 1, rows as f64 * 90.0);
        let brows = steps.checked_div(case.balance.every).unwrap_or(0) + 1;
        add("balance.csv", 1, 400.0 + brows as f64 * 230.0);
    }
    let fields = fields_bytes(case, mesh);
    if o.vtkhdf {
        add("fields.vtkhdf", 1, fields);
        if transient {
            let dt = case.time.map_or(0.0, |t| t.time_step);
            let series = steps.checked_div(o.fields_every).unwrap_or(0)
                + o.fields_at.len()
                + crate::unsteady::scheduled_outputs(case, steps, dt).0;
            add("fields_<step>.vtkhdf", series, fields);
        }
    }
    // a transient sweep writes the restart file of its last run, a steady
    // one none
    let single = transient
        || (case.verification.levels.len() <= 1 && case.verification.time_steps.is_empty());
    if crate::restart::AVAILABLE
        && single
        && (r.restart_final
            || r.restart_every > 0
            || !r.restart_at.is_empty()
            || r.restart_interval.is_some()
            || !r.restart_times.is_empty())
    {
        // per cell: u, p, the velocity and pressure gradients, the limiters
        // (19 values) and for a transient run the old levels (3 + 12); per
        // face the mass flux (and two old fluxes); per boundary face 4; the
        // multigrid matrix (diagonal and two slots per interior face)
        let per_cell = if transient { 34 } else { 19 };
        let per_face = if transient { 3 } else { 1 };
        let restart = 8.0
            * (nc * (per_cell + 2)
                + mesh.n_faces() * (per_face + 1)
                + mesh.n_boundary_faces() * 4
                + nc
                + 2 * mesh.n_interior_faces()) as f64
            + 16384.0;
        add("restart.h5", 1, restart);
        let dt = case.time.map_or(0.0, |t| t.time_step);
        add(
            "restart_<step>.h5",
            r.restart_at.len() + crate::unsteady::scheduled_outputs(case, steps, dt).1,
            restart,
        );
    }
    let walls: Vec<&String> = case
        .boundary
        .iter()
        .filter(|(_, b)| matches!(b, BoundaryCondition::Wall { .. }))
        .map(|(n, _)| n)
        .collect();
    if transient && (!walls.is_empty() || o.coefficients.is_some()) {
        let monitor_rows = steps / o.monitor_every + 1;
        add(
            "forces.csv",
            1,
            monitor_rows as f64 * (30.0 + 40.0 * walls.len() as f64),
        );
        if !o.probes.is_empty() {
            add(
                "probes.csv",
                1,
                monitor_rows as f64 * (30.0 + 56.0 * o.probes.len() as f64),
            );
        }
    }
    if o.wall_csv {
        for w in &walls {
            if let Some(p) = mesh.boundary(w) {
                add(&format!("wall_{w}.csv"), 1, p.len as f64 * 100.0);
            }
        }
    }
    for l in &o.lines {
        add(&format!("line_{}.csv", l.name), 1, l.points as f64 * 110.0);
    }
    add(
        "flow.log, flow.json, flow.md, result.json, run.jsonl, boundaries.csv, case.resolved.toml, SVG charts",
        1,
        2e4 + rows as f64 * 200.0,
    );
    v
}

/// The size of one `fields*.vtkhdf` of `case` on `mesh`.
fn fields_bytes(case: &FlowCase, mesh: &Mesh) -> f64 {
    #[cfg(feature = "hdf5")]
    {
        let mut n = if case.analytic.is_some() { 6 } else { 5 };
        let mut faces = 0;
        for d in crate::derived::sorted(&case.output.derived) {
            if d.on_faces() {
                faces += d.datasets().len();
            } else {
                n += d.datasets().len();
            }
        }
        foehncraft_mesh::vtkhdf::estimate_bytes(mesh, n)
            + (faces * 8 * mesh.n_boundary_faces()) as f64
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = case;
        (mesh.n_cells * 8 * 8 + mesh.face_nodes.len() * 16 + mesh.points.len() * 24) as f64
    }
}

/// The number of iterations or steps of a run and whether it is a bound.
fn step_count(case: &FlowCase) -> (usize, bool) {
    match case.time {
        Some(t) => {
            let dt = t.time_step;
            let n = t.steps.unwrap_or_else(|| {
                t.end_time
                    .map_or(0, |e| ((e - t.start_time) / dt).ceil() as usize)
            });
            let n = case.run.limits.max_steps.map_or(n, |m| m.min(n));
            (n, t.adaptive)
        }
        None => {
            let n = case.solver.max_iterations;
            (case.run.limits.max_iterations.map_or(n, |m| m.min(n)), true)
        }
    }
}

/// The measured time per cell and iteration or step of an earlier run of
/// this case in `dir` on the same CPU and backend (its `run.jsonl`).
fn earlier_run(dir: &Path, build: &BuildInfo, backend: &str) -> Option<f64> {
    let text = std::fs::read_to_string(dir.join("run.jsonl")).ok()?;
    let mut same_host = false;
    let mut ns = None;
    for l in text.lines() {
        let v: serde_json::Value = serde_json::from_str(l).ok()?;
        match v["type"].as_str() {
            Some("header") => {
                let p = &v["provenance"];
                same_host = p["cpu"] == build.cpu.as_str() && p["backend"] == backend;
            }
            Some("timings") => ns = v["ns_per_cell_step"].as_f64(),
            _ => {}
        }
    }
    ns.filter(|_| same_host)
}

/// The estimate of `case` (meshes built, nothing solved).
pub fn estimate(case: &FlowCase, dir: &Path, build: &BuildInfo) -> Result<Estimate> {
    let levels: Vec<Option<u32>> = if case.verification.levels.is_empty() {
        vec![None]
    } else {
        case.verification.levels.iter().map(|&l| Some(l)).collect()
    };
    let (n, approx) = step_count(case);
    let runs = if case.verification.time_steps.is_empty() {
        levels.len()
    } else {
        case.verification.time_steps.len()
    };
    let mut level_cells = Vec::new();
    let (mut faces, mut interior_faces) = (0usize, 0usize);
    let (mut host, mut dev, mut work) = (0.0f64, 0.0f64, 0.0f64);
    let mut cell_work = 0.0f64;
    let mut outs = Vec::new();
    let gpu = case.run.backend != "cpu";
    let mut steps = Vec::new();
    for (i, l) in levels.iter().enumerate() {
        let (mesh, _, _) = crate::case::mesh_for(case, *l)?;
        let (nc, nf) = (mesh.n_cells as f64, mesh.n_faces() as f64);
        let ni = mesh.n_interior_faces() as f64;
        let per_faces =
            HOST_BYTES_PER_INTERIOR_FACE * ni + HOST_BYTES_PER_BOUNDARY_FACE * (nf - ni);
        level_cells.push(mesh.n_cells);
        faces = faces.max(mesh.n_faces());
        interior_faces = interior_faces.max(mesh.n_interior_faces());
        host = host.max(HOST_BYTES_BASE + per_faces);
        if gpu {
            dev = dev.max(DEVICE_SHARE * per_faces);
        }
        // a time-step sweep runs the same mesh with smaller steps
        let level_steps: Vec<usize> = if case.verification.time_steps.is_empty() {
            vec![n]
        } else {
            let t = case.time.unwrap_or_default();
            case.verification
                .time_steps
                .iter()
                .map(|&dt| {
                    t.end_time
                        .map_or(n, |e| ((e - t.start_time) / dt).ceil() as usize)
                })
                .collect()
        };
        for s in &level_steps {
            work += nc * *s as f64;
            cell_work += nc * nc * *s as f64;
            steps.push(*s);
        }
        if i == levels.len() - 1 {
            outs = outputs(case, &mesh, *level_steps.last().unwrap_or(&n));
        }
    }
    let _ = runs;
    let output_bytes = outs.iter().map(|o| o.count as f64 * o.bytes).sum();
    let transient = case.time.is_some();
    let (wall_time, basis) = if let Some(ns) = earlier_run(dir, build, &case.run.backend) {
        (
            Some(1e-9 * ns * work),
            format!(
                "an earlier run of this case in {} on this CPU ({ns:.0} ns per cell and step)",
                dir.display()
            ),
        )
    } else if let Some(b) = BENCHMARKS.iter().find(|b| {
        let hw = if case.run.backend == "cpu" {
            build.cpu.as_str()
        } else {
            build.gpu.as_deref().unwrap_or("")
        };
        b.1 == case.run.backend && b.2 == transient && hw.contains(b.0)
    }) {
        (
            Some(1e-9 * b.3 * work),
            format!("the benchmark of {} ({} ns per cell and step)", b.0, b.3),
        )
    } else {
        (
            None,
            format!(
                "{}: no benchmark of this hardware ({}, backend {}); a first run writes one into run.jsonl",
                diag::RUN_NO_BENCHMARK,
                if build.cpu.is_empty() {
                    "unknown CPU"
                } else {
                    &build.cpu
                },
                case.run.backend
            ),
        )
    };
    // cell-hours: each level's cells for its share of the wall time
    // (the time is proportional to cells times steps, so the sum of
    // cells times its time is the work-weighted mean)
    let cell_hours = wall_time.map(|w| {
        if work > 0.0 {
            w / work * cell_work / 3600.0
        } else {
            0.0
        }
    });
    let gpu_hours = wall_time.map(|w| f64::from(u8::from(gpu)) * w / 3600.0);
    let iteration_budget = case
        .time
        .is_none()
        .then(|| crate::budget::IterationBudget::of(case.solver.max_iterations, &case.run.limits));
    Ok(Estimate {
        cell_hours,
        gpu_hours,
        iteration_budget,
        cells: level_cells.iter().copied().max().unwrap_or(0),
        faces,
        interior_faces,
        level_cells,
        ranks: 1,
        gpus: usize::from(gpu),
        steps,
        steps_approximate: approx,
        host_memory: host,
        device_memory: dev,
        output_bytes,
        outputs: outs,
        wall_time,
        wall_time_basis: basis,
    })
}

/// The findings of the iteration and study budgets (D198) for a run of
/// `case` into `out` with the estimate `e`.
fn budget_findings(case: &FlowCase, e: &Estimate, out: &Path) -> Vec<Diagnostic> {
    let mut v = Vec::new();
    if let Some(b) = &e.iteration_budget
        && b.iterations > crate::budget::LARGE_ITERATION_BUDGET
    {
        v.push(
            Diagnostic::new(
                diag::RUN_LARGE_BUDGET,
                format!(
                    "the steady iteration budget is {} iterations ({}); a steady case needs hundreds, at most a few thousand: \
                     set it to about 3000 and look at the diagnosis at [run.limits] diagnose_at",
                    b.iterations, b.source
                ),
            )
            .key(b.source),
        );
    }
    if case.run.limits.has_study_budget() {
        let sb = study_budget(case, e, out);
        if let Some((k, b, used)) = sb.exhausted() {
            v.push(
                Diagnostic::new(
                    diag::RUN_STUDY_SPENT,
                    format!(
                        "[run.limits] {k} = {b} is spent: the earlier runs of the study spent {used:.4} ({})",
                        sb.ledger
                    ),
                )
                .key(format!("run.limits.{k}")),
            );
        } else {
            for k in sb.over_estimate() {
                v.push(
                    Diagnostic::new(
                        diag::RUN_STUDY_ESTIMATE,
                        format!("the estimate of this run does not fit into what is left of [run.limits] {k}"),
                    )
                    .key(format!("run.limits.{k}")),
                );
            }
        }
    }
    v
}

/// The study budgets of `case` against the ledger of the study of `out`.
fn study_budget(case: &FlowCase, e: &Estimate, out: &Path) -> crate::budget::StudyBudget {
    crate::budget::StudyBudget::new(
        &case.run.limits,
        &crate::budget::ledger_path(out),
        (e.cell_hours, e.gpu_hours, e.output_bytes),
    )
}

/// The `"budget"` row of the run record: the iteration budget, the D81
/// estimate and the study budgets (D198).
fn budget_row(
    case: &FlowCase,
    e: Option<&Estimate>,
    study: Option<&crate::budget::StudyBudget>,
) -> serde_json::Value {
    let mut row = serde_json::json!({"type": "budget"});
    if case.time.is_none() {
        row["iterations"] = serde_json::to_value(crate::budget::IterationBudget::of(
            case.solver.max_iterations,
            &case.run.limits,
        ))
        .unwrap_or_default();
    }
    if let Some(e) = e {
        row["estimate"] = serde_json::json!({
            "cells": e.cells,
            "steps": e.steps,
            "wall_time": e.wall_time,
            "wall_time_basis": e.wall_time_basis,
            "cell_hours": e.cell_hours,
            "gpu_hours": e.gpu_hours,
            "output_bytes": e.output_bytes,
        });
    }
    if let Some(s) = study {
        row["study"] = serde_json::to_value(s).unwrap_or_default();
    }
    row
}

/// The budgets at the start of a run (D198): the `"budget"` row of the run
/// record, the log lines, and with a study budget its ledger; a spent
/// study budget stops the run before it starts (`RUN019`).
fn start_budgets(
    case: &FlowCase,
    dir: &Path,
    build: &BuildInfo,
    roots: &Roots,
    steer: &mut Steer,
) -> Result<()> {
    let limits = &case.run.limits;
    if let Some(b) = case
        .time
        .is_none()
        .then(|| crate::budget::IterationBudget::of(case.solver.max_iterations, limits))
    {
        steer.notes.push(format!(
            "iteration budget {} ({}); {}",
            b.iterations,
            b.source,
            match b.diagnose_at {
                Some(n) => format!(
                    "diagnosis at iteration {n} ([run.limits] diagnose_at, {})",
                    if b.diagnosis == "stop" {
                        "stops a run off track"
                    } else {
                        "recorded"
                    }
                ),
                None => "no diagnosis".into(),
            }
        ));
    }
    if !limits.has_study_budget() {
        if case.time.is_none() {
            steer.record.row(&budget_row(case, None, None));
        }
        return Ok(());
    }
    // the ledger is in the parent of the output directory (D82)
    roots.check_write(&crate::budget::ledger_path(dir))?;
    let e = estimate(case, dir, build)?;
    let sb = study_budget(case, &e, dir);
    steer.record.row(&budget_row(case, Some(&e), Some(&sb)));
    steer.notes.extend(sb.lines());
    if let Some((k, b, used)) = sb.exhausted() {
        return Err(diag::err(
            Diagnostic::new(
                diag::RUN_STUDY_SPENT,
                format!(
                    "[run.limits] {k} = {b} is spent: the earlier runs of the study spent {used:.4} ({}); the run does not start",
                    sb.ledger
                ),
            )
            .key(format!("run.limits.{k}")),
        ));
    }
    for k in sb.over_estimate() {
        let m =
            format!("the estimate of this run does not fit into what is left of [run.limits] {k}");
        steer.event(diag::RUN_STUDY_ESTIMATE, 0, &m);
        steer.notes.push(m);
    }
    steer.spend.before = Some(crate::budget::spent(&crate::budget::ledger_path(dir)));
    Ok(())
}

/// The report of `foehn flow --check` (D78, D81).
#[derive(Debug, Clone, Serialize, schemars::JsonSchema)]
pub struct CheckReport {
    /// `foehncraft.check`.
    pub format: &'static str,
    /// [`FORMAT_VERSION`].
    pub format_version: &'static str,
    /// Case file.
    pub case_file: String,
    /// The case can be run (no error).
    pub ok: bool,
    /// Findings of the validation and the boundary-condition check.
    pub diagnostics: Vec<Diagnostic>,
    /// The estimate (when the case could be read).
    pub estimate: Option<Estimate>,
    /// The operator exactness check of the first mesh (D236): the cells
    /// where the velocity gradient or Hessian is not exact, worst first.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub operators: Option<crate::operators::OperatorReport>,
}

/// Checks a case file: reads it (with overrides and roots), checks the
/// boundary conditions on the mesh of the first level and estimates the
/// run. Never fails: every problem is a diagnostic.
pub fn check_file(
    path: &Path,
    sets: &[String],
    roots: &Roots,
    dir: Option<&Path>,
    build: &BuildInfo,
) -> CheckReport {
    let mut rep = CheckReport {
        format: "foehncraft.check",
        format_version: FORMAT_VERSION,
        case_file: path.display().to_string(),
        ok: false,
        diagnostics: Vec::new(),
        estimate: None,
        operators: None,
    };
    let case = match load(path, sets, roots) {
        Ok(c) => c,
        Err(e) => {
            rep.diagnostics.push(diag::of_error(&e));
            return rep;
        }
    };
    check_case(&case, path, dir, roots, build, &mut rep);
    rep
}

/// As [`check_file`] for a case already read.
pub fn check_case(
    case: &FlowCase,
    path: &Path,
    dir: Option<&Path>,
    roots: &Roots,
    build: &BuildInfo,
    rep: &mut CheckReport,
) {
    let (_, out) = crate::case::output_dir(case, path, dir.map(Path::to_path_buf));
    if let Err(e) = roots.check_write(&out) {
        rep.diagnostics.push(diag::of_error(&e));
    }
    let level = case.verification.levels.first().copied();
    match crate::case::mesh_for(case, level) {
        Ok((mesh, geo, _)) => {
            // the property program on this mesh: regions name existing zones (D112)
            let phys_on = match crate::case::physics(case)
                .and_then(|p| crate::case::physics_on(case, &p, &mesh, &geo))
            {
                Ok(p) => Some(p),
                Err(e) => {
                    rep.diagnostics.push(diag::of_error(&e));
                    None
                }
            };
            let rho: crate::bc::Density = match &phys_on {
                Some(p) => p.into(),
                None => case.rho().into(),
            };
            let t0 = case.time.map_or(0.0, |t| t.start_time);
            let values = case.compiled_boundaries().unwrap_or_default();
            let mut findings = crate::bc::check(
                &mesh,
                &geo,
                &case.boundary,
                rho,
                crate::bc::At {
                    analytic: case.analytic_at(t0).as_ref(),
                    values: &values,
                    t: t0,
                    turning: &case.turning(),
                    hydrostatic: crate::buoyancy::hydrostatic(case),
                },
                case.check_context(),
            );
            findings.extend(crate::turbulence::setup_findings(case));
            findings.extend(crate::energy::setup_findings(case));
            findings.extend(crate::buoyancy::setup_findings(case, &mesh, &geo));
            rep.diagnostics
                .extend(findings.iter().map(crate::bc::Finding::diagnostic));
            // the operator exactness check (D236)
            if let Ok(mut r) = crate::operators::check_case(case, &mesh, &geo) {
                let ds = crate::operators::diagnostics(&r, case);
                r.findings = ds.iter().map(|d| d.message.clone()).collect();
                rep.diagnostics.extend(ds);
                rep.operators = Some(r);
            }
        }
        Err(e) => rep.diagnostics.push(diag::of_error(&e)),
    }
    match estimate(case, &out, build) {
        Ok(e) => {
            rep.diagnostics.extend(budget_findings(case, &e, &out));
            rep.estimate = Some(e);
        }
        Err(e) => rep.diagnostics.push(diag::of_error(&e)),
    }
    rep.ok = !rep
        .diagnostics
        .iter()
        .any(|d| d.severity == Severity::Error);
}

// ---------------------------------------------------------------------------
// run and result (D62, D64)

/// How a run ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Status {
    /// Finished, every check passed (exit code 0).
    Finished,
    /// An error, including divergence (exit code 1).
    Error,
    /// Finished with failed checks (exit code 2).
    FailedChecks,
    /// Stopped on request: signal, `.stop` or a limit (exit code 3).
    Stopped,
}

impl Status {
    /// Exit code of D64.
    pub fn exit_code(self) -> i32 {
        match self {
            Self::Finished => 0,
            Self::Error => 1,
            Self::FailedChecks => 2,
            Self::Stopped => 3,
        }
    }
}

/// Options of [`run`].
#[derive(Debug, Clone, Default)]
pub struct RunOptions {
    /// Requests (signals).
    pub interrupts: Interrupts,
    /// Build and host.
    pub build: BuildInfo,
    /// File confinement.
    pub roots: Roots,
    /// Sees every row of the run record ([`Observer`]).
    pub observer: Option<Observer>,
    /// The tool call that started the run (tool, version, typed input,
    /// D197): written into the run record after its header, so that
    /// `foehn run --script run.jsonl` runs it again.
    pub call: Option<serde_json::Value>,
    /// Listen on the live control port (L19, D198 (f)): `run.status`,
    /// `run.stop` and `run.set` reach the run through `control.json` and
    /// `control.token` in its output directory.
    pub control: bool,
}

/// The result summary (`result.json`, D62): what a study or an agent needs
/// from a run, without timings.
#[derive(Debug, Clone, Serialize, schemars::JsonSchema)]
pub struct ResultSummary {
    /// `foehncraft.result`.
    pub format: &'static str,
    /// [`FORMAT_VERSION`].
    pub format_version: &'static str,
    /// How the run ended.
    pub status: Status,
    /// Exit code of D64.
    pub exit_code: i32,
    /// Why it stopped early.
    pub stop: Option<Stop>,
    /// Case name.
    pub case: String,
    /// Overrides applied (`--set`).
    pub overrides: Vec<String>,
    /// Where it comes from.
    pub provenance: Option<Provenance>,
    /// Errors, failed checks, warnings and notes, with codes.
    pub diagnostics: Vec<Diagnostic>,
    /// Reynolds number of the case, if it defines one.
    pub reynolds: Option<f64>,
    /// One entry per solve (level or time step of a sweep): mesh, cells,
    /// mesh hash, iterations, convergence, residuals, mass imbalance,
    /// boundaries (mass flow, mean pressure, force), coefficients, errors,
    /// benchmark, transient summary.
    pub solves: Vec<serde_json::Value>,
    /// Observed orders in space.
    pub orders: Vec<f64>,
    /// Observed orders in time.
    pub time_orders: Vec<f64>,
    /// Observed orders of the passive scalars with an exact solution (D51):
    /// name and orders; `<name> (time)` for a time-step sweep.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub scalar_orders: Vec<(String, Vec<f64>)>,
    /// Discretisation uncertainty of a sweep (D67): per monitored quantity
    /// the observed order, the extrapolated value and the GCI.
    pub uncertainty: Vec<crate::uncertainty::Estimate>,
    /// Output files with sizes and SHA-256 hashes (D84).
    pub files: Vec<FileHash>,
    /// What the user must know before trusting the result (D199): codes
    /// LIM001 to LIM006, the feature and a message. They change neither
    /// `status` nor the exit code.
    pub limitations: Vec<crate::limitations::Limitation>,
}

/// The outcome of [`run`].
#[derive(Debug, Clone)]
pub struct Outcome {
    /// The report of the solver (`None` after an error).
    pub report: Option<FlowReport>,
    /// The result summary, also in `result.json`.
    pub summary: ResultSummary,
}

impl Outcome {
    /// Exit code of D64.
    pub fn exit_code(&self) -> i32 {
        self.summary.exit_code
    }
}

/// Keys of a solve that hold timings, left out of `result.json`.
const TIMING_KEYS: &[&str] = &[
    "timings",
    "ns_per_cell_iteration",
    "ns_per_cell_step",
    "device",
];

fn strip_timings(v: &mut serde_json::Value) {
    match v {
        serde_json::Value::Object(m) => {
            m.retain(|k, _| !TIMING_KEYS.contains(&k.as_str()));
            for x in m.values_mut() {
                strip_timings(x);
            }
        }
        serde_json::Value::Array(a) => a.iter_mut().for_each(strip_timings),
        _ => {}
    }
}

/// An error result (no report).
pub fn error_summary(
    case: Option<&FlowCase>,
    name: &str,
    e: &anyhow::Error,
    build: &BuildInfo,
) -> ResultSummary {
    ResultSummary {
        format: "foehncraft.result",
        format_version: FORMAT_VERSION,
        status: Status::Error,
        exit_code: Status::Error.exit_code(),
        stop: None,
        case: name.into(),
        overrides: case.map(|c| c.overrides.clone()).unwrap_or_default(),
        provenance: case.map(|c| provenance(c, build, None)),
        diagnostics: vec![diag::of_error(e)],
        reynolds: case.and_then(FlowCase::reynolds),
        solves: Vec::new(),
        orders: Vec::new(),
        time_orders: Vec::new(),
        scalar_orders: Vec::new(),
        uncertainty: Vec::new(),
        files: Vec::new(),
        limitations: Vec::new(),
    }
}

/// Runs a case (`name` is the case name, `dir` the output directory) and
/// writes `result.json` and `run.jsonl` next to the solver's output. An
/// error of the solver is an [`Outcome`] with [`Status::Error`]; `Err` only
/// when the output directory cannot be written.
pub fn run(
    case: &FlowCase,
    name: &str,
    dir: &Path,
    opts: &RunOptions,
    log: impl FnMut(&str),
) -> Result<Outcome> {
    let mut log = log;
    opts.roots.check_write(dir)?;
    std::fs::create_dir_all(dir).map_err(|e| {
        diag::err(Diagnostic::new(
            diag::RUN_IO,
            format!("creating the output directory {}: {e}", dir.display()),
        ))
    })?;
    let mut case = case.clone();
    let mut steer = Steer::new(opts.interrupts.clone(), dir, true)?;
    // box lines need a UTF-8 terminal; the asterisks are plain ASCII (D214)
    if case.log.style == crate::console::Style::Box && !crate::console::utf8_locale() {
        case.log.style = crate::console::Style::Asterisk;
        steer.notes.push(
            "the terminal is not UTF-8 (LC_ALL, LC_CTYPE, LANG): log.style \"box\" falls back to \"asterisk\"".into(),
        );
    }
    steer.record.observer = opts.observer.clone();
    steer.build = Some(opts.build.clone());
    // the live control port lives as long as the run (L19)
    let _port = if opts.control {
        let live = crate::live::Live::default();
        steer.record.live(live.clone());
        steer.live = Some(live.clone());
        Some(crate::live::Port::open(dir, &live, &opts.interrupts)?)
    } else {
        None
    };
    let res = resolve_latest(&mut case, dir)
        .and_then(|()| check_paths(&case, &opts.roots))
        .and_then(|()| {
            let prov = provenance(&case, &opts.build, None);
            steer.record.row(&serde_json::json!({
                "type": "header",
                "format": "foehncraft.run_record",
                "format_version": FORMAT_VERSION,
                "case": name,
                "date": run_date(),
                "provenance": prov,
                "threads": case.run.threads,
                "overrides": case.overrides,
                "case_file": case.source.as_ref().map(|p| p.display().to_string()),
                "settings": case,
            }));
            if let Some(c) = &opts.call {
                let mut row = serde_json::json!({"type": "call"});
                if let (Some(r), Some(c)) = (row.as_object_mut(), c.as_object()) {
                    r.extend(c.clone());
                }
                steer.record.row(&row);
            }
            start_budgets(&case, dir, &opts.build, &opts.roots, &mut steer)?;
            steer.record.flush();
            crate::case::run_steered(&case, name, dir, &mut steer, &mut log)
        });
    // the spend of this run into the study ledger (D198 (c))
    if case.run.limits.has_study_budget() && steer.spend.before.is_some() {
        let s = steer.spent();
        steer.record.row(&serde_json::json!({
            "type": "spend",
            "cell_hours": s.cell_hours,
            "gpu_hours": s.gpu_hours,
            "output_bytes": s.output_bytes,
        }));
        if let Err(e) = crate::budget::append(&crate::budget::ledger_path(dir), name, dir, &s) {
            log(&crate::console::tagged(
                "warning",
                &format!("the study ledger could not be written: {e}"),
            ));
        }
    }
    let summary = match &res {
        Ok(rep) => summary_of(&case, rep, dir, &opts.build)?,
        Err(e) => {
            let mut s = error_summary(Some(&case), name, e, &opts.build);
            s.files = file_hashes(dir)?;
            s
        }
    };
    if let Ok(rep) = &res {
        let cells: usize = rep.solves.iter().map(|s| s.cells).sum();
        let steps: usize = rep
            .solves
            .iter()
            .map(|s| s.transient.as_ref().map_or(s.iterations, |t| t.steps))
            .sum();
        let work: f64 = rep
            .solves
            .iter()
            .map(|s| s.cells as f64 * s.transient.as_ref().map_or(s.iterations, |t| t.steps) as f64)
            .sum();
        steer.record.row(&serde_json::json!({
            "type": "timings",
            "seconds": rep.seconds,
            "cells": cells,
            "steps": steps,
            "ns_per_cell_step": if work > 0.0 { 1e9 * rep.seconds / work } else { 0.0 },
            "solves": rep.solves.iter().map(|s| serde_json::json!({
                "mesh": s.mesh,
                "device": s.device,
                "timings": s.timings,
                "ns_per_cell_iteration": s.ns_per_cell_iteration,
                "ns_per_cell_step": s.transient.as_ref().map(|t| t.ns_per_cell_step),
            })).collect::<Vec<_>>(),
        }));
    }
    let text = serde_json::to_string_pretty(&summary)?;
    std::fs::write(dir.join("result.json"), format!("{text}\n"))?;
    events_frame(&case, &steer.events, dir, &mut log);
    closing_frame(
        &case,
        &summary,
        dir,
        steer.start.elapsed().as_secs_f64(),
        &mut log,
    );
    steer
        .record
        .row(&serde_json::json!({"type": "result", "result": summary}));
    steer.record.row(&serde_json::json!({"type": "end", "status": summary.status, "exit_code": summary.exit_code}));
    steer.record.flush();
    Ok(Outcome {
        report: res.ok(),
        summary,
    })
}

/// The Events block (D214): every change made during the run, with its
/// iteration or step, its wall-clock time and its code, before the closing
/// frame. Nothing when nothing changed. A steered run can be reproduced
/// from it; the rows stay in `run.jsonl`.
fn events_frame(case: &FlowCase, events: &[EventLine], dir: &Path, log: &mut impl FnMut(&str)) {
    use std::io::Write as _;
    let Some(f) = crate::logblocks::events(events, case.time.is_some()) else {
        return;
    };
    let mut lines = vec![String::new()];
    lines.extend(f.render(case.log.style));
    for l in &lines {
        log(l);
    }
    if let Ok(mut file) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(dir.join("flow.log"))
    {
        let _ = writeln!(file, "{}", lines.join("\n"));
    }
}

/// The last frame of the log: the verdict and the exit code. It goes to the
/// terminal and to the end of `flow.log`.
fn closing_frame(
    case: &FlowCase,
    s: &ResultSummary,
    dir: &Path,
    secs: f64,
    log: &mut impl FnMut(&str),
) {
    use std::io::Write as _;
    let first_error = s
        .diagnostics
        .iter()
        .find(|d| d.severity == diag::Severity::Error);
    let n = s.limitations.len();
    let title = match s.status {
        Status::Finished if n > 0 => format!(
            "Run finished with {n} limitation{}: every check passed",
            if n == 1 { "" } else { "s" }
        ),
        Status::Finished => "Run finished: every check passed".to_string(),
        Status::FailedChecks => "Run finished with failed checks".to_string(),
        Status::Stopped => format!(
            "Run stopped: {}",
            s.stop.as_ref().map_or("", |x| x.message.as_str())
        ),
        Status::Error => format!(
            "Run FAILED: {}",
            first_error.map_or(String::new(), |d| d.code.to_string())
        ),
    };
    let mut f = crate::console::Frame::new(title).summary();
    if let (Status::Error, Some(d)) = (s.status, first_error) {
        f.push_kv("Reason", &d.message);
    }
    let mut f = f
        .kv("Exit code", s.exit_code)
        .kv("Elapsed time", crate::console::clock(secs))
        .kv("Results in", dir.display());
    if dir.join("flow.md").exists() {
        f.push_kv("Report", "flow.md, result.json");
    }
    for l in &s.limitations {
        f.push_kv(l.code, &l.message);
    }
    let mut lines = vec![String::new()];
    lines.extend(f.render(case.log.style));
    for l in &lines {
        log(l);
    }
    if let Ok(mut file) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(dir.join("flow.log"))
    {
        let _ = writeln!(file, "{}", lines.join("\n"));
    }
}

/// The result summary of a finished report.
pub fn summary_of(
    case: &FlowCase,
    rep: &FlowReport,
    dir: &Path,
    build: &BuildInfo,
) -> Result<ResultSummary> {
    let mut diagnostics: Vec<Diagnostic> = rep
        .findings
        .iter()
        .map(crate::bc::Finding::diagnostic)
        .collect();
    diagnostics.extend(rep.warnings.iter().cloned());
    for f in &rep.failures {
        let code = if f.contains("not converged") {
            diag::RUN_NOT_CONVERGED
        } else {
            diag::RUN_CHECK_FAILED
        };
        diagnostics.push(Diagnostic::new(code, f.clone()));
    }
    if let Some(s) = &rep.stop {
        diagnostics.push(Diagnostic::new(s.code, s.message.clone()));
    }
    let status = if rep.stop.is_some() {
        Status::Stopped
    } else if !rep.failures.is_empty() {
        Status::FailedChecks
    } else {
        Status::Finished
    };
    let limitations = crate::limitations::of_run(case, rep);
    if !limitations.is_empty() {
        // the list goes to the top of the report (D199)
        let md = dir.join("flow.md");
        if let Ok(text) = std::fs::read_to_string(&md) {
            let block = crate::limitations::markdown(&limitations);
            let at = text.find("\n\n").map_or(0, |i| i + 2);
            let _ = std::fs::write(&md, format!("{}{block}{}", &text[..at], &text[at..]));
        }
    }
    let gpu = rep
        .solves
        .iter()
        .map(|s| s.device.clone())
        .find(|d| !d.starts_with("CPU"));
    let solves = rep
        .solves
        .iter()
        .map(|s| {
            let mut v = serde_json::to_value(s).unwrap_or_default();
            strip_timings(&mut v);
            v
        })
        .collect();
    Ok(ResultSummary {
        format: "foehncraft.result",
        format_version: FORMAT_VERSION,
        status,
        exit_code: status.exit_code(),
        stop: rep.stop.clone(),
        case: rep.name.clone(),
        overrides: case.overrides.clone(),
        provenance: Some(provenance(case, build, gpu)),
        diagnostics,
        reynolds: rep.reynolds,
        solves,
        orders: rep.orders.clone(),
        time_orders: rep.time_orders.clone(),
        scalar_orders: rep.scalar_orders.clone(),
        uncertainty: rep.uncertainty.clone(),
        files: file_hashes(dir)?,
        limitations,
    })
}

/// Reads the result summary of a finished run in `dir`.
pub fn read_result(dir: &Path) -> Result<serde_json::Value> {
    let p = dir.join("result.json");
    let text = std::fs::read_to_string(&p).with_context(|| format!("reading {}", p.display()))?;
    Ok(serde_json::from_str(&text)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exit_codes() {
        assert_eq!(Status::Finished.exit_code(), 0);
        assert_eq!(Status::Error.exit_code(), 1);
        assert_eq!(Status::FailedChecks.exit_code(), 2);
        assert_eq!(Status::Stopped.exit_code(), 3);
    }

    #[test]
    fn interrupts_round_trip() {
        let i = Interrupts::new();
        assert_eq!(i.stop(), None);
        assert!(!i.take_restart());
        i.request_restart();
        assert!(i.take_restart());
        assert!(!i.take_restart());
        i.stop_flag()
            .store(StopSignal::Term.value(), Ordering::SeqCst);
        assert_eq!(i.stop(), Some(StopSignal::Term));
    }

    #[test]
    fn roots_confine_paths() {
        let base = std::env::temp_dir().join("foehncraft_roots_test");
        let _ = std::fs::create_dir_all(base.join("cases"));
        let r = Roots {
            read: vec![base.join("cases")],
            write: vec![base.join("out")],
        };
        assert!(r.check_read(&base.join("cases/a.toml"), "case").is_ok());
        assert!(
            r.check_read(&base.join("cases/../secret.toml"), "case")
                .is_err()
        );
        assert!(r.check_write(&base.join("out/run1")).is_ok());
        let e = r.check_write(&base.join("cases")).unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::CASE_OUTSIDE_ROOT);
        assert!(Roots::default().check_write(Path::new("/anywhere")).is_ok());
    }

    #[test]
    fn dates_follow_source_date_epoch_format() {
        let d = run_date();
        assert_eq!(d.len(), 20, "{d}");
        assert!(d.ends_with('Z'));
    }

    #[test]
    fn sha256_of_known_input() {
        assert_eq!(
            sha256(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }
}
