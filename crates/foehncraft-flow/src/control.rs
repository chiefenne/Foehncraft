//! Steering a running case (D40): control files in the output directory and
//! re-reading the case file.
//!
//! | file | effect |
//! |---|---|
//! | `.stop` | the run writes the restart file and the results of the current iteration or time step and ends ("stopped") |
//! | `.result_N` | at iteration or step N the run writes `fields_N.vtkhdf` and goes on |
//! | `.read_toml_N` | at iteration or step N the run reads the case file again and applies what may change during a run |
//!
//! The files are checked after every outer iteration (steady) or time step
//! (transient) and never deleted by the solver. A start with an existing
//! `.stop` fails at once; a `.result_N` or `.read_toml_N` whose N lies in
//! the past is ignored. `[run] reread_every = N` re-reads the case file
//! every N iterations or steps without a control file.
//!
//! What a re-read applies ([`RUNTIME_KEYS`]): the output settings, the
//! restart and re-read schedule, the stop criteria, the relaxation factors,
//! the momentum scheme between `upwind` and `linear-upwind` (a first-order
//! start, then second order), the linear-solver tolerances, the time-step settings and the balance
//! report and tolerances. Every other
//! difference is reported as not applied (the state of the run belongs to
//! the case it started with; stop the run and restart from the restart file
//! to change those). A case file that no longer validates leaves everything
//! as it was. The log lists every change (`key: old -> new`) and
//! `case.read_toml_N.toml` holds the settings in force afterwards. The
//! entries of `[[schedule]]` that have fired stay in force over the file's
//! values; the schedule itself cannot be re-read (D128).

use crate::case::FlowCase;
use anyhow::{Result, bail};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// Keys (`section.key`) that a re-read may change.
pub const RUNTIME_KEYS: &[&str] = &[
    "output.fields_every",
    "output.fields_at",
    "output.fields_interval",
    "output.fields_times",
    "output.history_every",
    "output.monitor_every",
    "output.monitor_interval",
    "output.vtkhdf",
    "output.wall_csv",
    "output.derived",
    "output.reference.velocity",
    "output.reference.pressure",
    "output.reference.density",
    "run.restart_every",
    "run.restart_at",
    "run.restart_interval",
    "run.restart_times",
    "run.restart_final",
    "run.reread_every",
    "run.limits.max_wall_time",
    "run.limits.max_steps",
    "run.limits.max_iterations",
    "run.limits.max_output",
    "run.limits.study_cell_hours",
    "run.limits.study_gpu_hours",
    "run.limits.study_output",
    "solver.max_iterations",
    "solver.tolerance",
    "solver.momentum.relaxation",
    "solver.momentum.scheme",
    "solver.pressure.relaxation",
    "solver.momentum.tolerance",
    "solver.momentum.max_iter",
    "solver.pressure.tolerance",
    "solver.pressure.max_iter",
    "solver.coupled.tolerance",
    "solver.coupled.max_iter",
    "solver.coupled.momentum_sweep",
    "solver.coupled.momentum_cycles",
    "solver.coupled.momentum_reuse",
    "solver.coupled.momentum_tolerance",
    "solver.coupled.momentum_max_iter",
    "solver.coupled.pressure_tolerance",
    "solver.coupled.pressure_max_iter",
    "solver.coupled.pressure_reuse",
    "solver.pseudo_time.step",
    "solver.pseudo_time.factor",
    "solver.pseudo_time.momentum",
    "solver.pseudo_time.turbulence",
    "solver.correctors",
    "solver.outer_iterations",
    "solver.outer_tolerance",
    "solver.amg_rebuild",
    "solver.log_every",
    "time.end_time",
    "time.steps",
    "time.time_step",
    "time.adaptive",
    "time.max_cfl",
    "time.max_growth",
    "time.min_time_step",
    "time.max_time_step",
    "time.error_tolerance",
    "time.error_absolute",
    "balance.every",
    "balance.tolerance_step",
    "balance.tolerance_run",
    "balance.tolerance_energy",
    "balance.action",
];

/// What the control files ask for at one iteration or step.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Requests {
    /// `.stop` exists.
    pub stop: bool,
    /// A `.result_N` for this step.
    pub result: bool,
    /// A `.read_toml_N` for this step (or the re-read schedule).
    pub reread: bool,
}

/// The control files of one output directory.
pub struct Control {
    dir: PathBuf,
}

impl Control {
    /// Control files in `dir`; fails if `.stop` exists already.
    pub fn new(dir: &Path) -> Result<Self> {
        if dir.join(".stop").exists() {
            bail!(
                "{} exists: remove it before starting the run (`rm {}`)",
                dir.join(".stop").display(),
                dir.join(".stop").display()
            );
        }
        Ok(Self {
            dir: dir.to_path_buf(),
        })
    }

    /// The requests at iteration or step `n`.
    pub fn check(&self, n: usize, reread_every: usize) -> Requests {
        let mut r = Requests {
            stop: self.dir.join(".stop").exists(),
            reread: reread_every > 0 && n.is_multiple_of(reread_every),
            ..Default::default()
        };
        if let Ok(rd) = std::fs::read_dir(&self.dir) {
            for e in rd.flatten() {
                let name = e.file_name();
                let Some(name) = name.to_str() else { continue };
                if let Some(k) = name.strip_prefix(".result_") {
                    r.result |= k.parse::<usize>() == Ok(n);
                } else if let Some(k) = name.strip_prefix(".read_toml_") {
                    r.reread |= k.parse::<usize>() == Ok(n);
                }
            }
        }
        r
    }
}

pub(crate) fn flatten(prefix: &str, v: &toml::Value, out: &mut BTreeMap<String, toml::Value>) {
    match v {
        toml::Value::Table(t) => {
            for (k, x) in t {
                let key = if prefix.is_empty() {
                    k.clone()
                } else {
                    format!("{prefix}.{k}")
                };
                flatten(&key, x, out);
            }
        }
        other => {
            out.insert(prefix.to_string(), other.clone());
        }
    }
}

pub(crate) fn set(v: &mut toml::Value, key: &str, x: Option<&toml::Value>) {
    let mut parts: Vec<&str> = key.split('.').collect();
    let last = parts.pop().unwrap_or_default();
    let mut cur = v;
    for p in parts {
        let toml::Value::Table(t) = cur else { return };
        cur = t
            .entry(p.to_string())
            .or_insert_with(|| toml::Value::Table(Default::default()));
    }
    if let toml::Value::Table(t) = cur {
        match x {
            Some(x) => {
                t.insert(last.to_string(), x.clone());
            }
            None => {
                t.remove(last);
            }
        }
    }
}

/// A case with one setting changed (L19: `run.set` through the live control
/// port): the new case, the old and the new value as the case holds them.
/// An error when the value does not fit the key's type or the case no
/// longer validates.
pub fn with_setting(
    current: &FlowCase,
    key: &str,
    value: &toml::Value,
) -> Result<(FlowCase, String, String)> {
    let mut v = toml::Value::try_from(current)?;
    let mut flat = BTreeMap::new();
    flatten("", &v, &mut flat);
    let old = show(flat.get(key));
    set(&mut v, key, Some(value));
    let mut case: FlowCase = v.try_into()?;
    case.source = current.source.clone();
    case.overrides = current.overrides.clone();
    case.pressure_frame = current.pressure_frame.clone();
    case.validate()?;
    let mut flat = BTreeMap::new();
    flatten("", &toml::Value::try_from(&case)?, &mut flat);
    let new = show(flat.get(key));
    Ok((case, old, new))
}

/// The result of a re-read.
#[derive(Debug, Clone)]
pub struct Reread {
    /// The case in force afterwards.
    pub case: FlowCase,
    /// Applied changes, `key: old -> new`.
    pub applied: Vec<String>,
    /// Differences that cannot change during a run.
    pub not_applied: Vec<String>,
}

pub(crate) fn show(v: Option<&toml::Value>) -> String {
    v.map_or("(default)".into(), |v| v.to_string())
}

/// Reads the case file `path` again and applies its run-time settings to
/// `current`. An unreadable or invalid file is an error (the caller keeps
/// `current`).
pub fn reread(path: &Path, current: &FlowCase, fired: usize) -> Result<Reread> {
    let mut new = FlowCase::from_file_with(path, &current.overrides)?;
    // the entries of `[[schedule]]` that have fired stay in force over the
    // values of the file (D128)
    if fired > 0 {
        new = crate::schedule::apply_fired(&new, current, fired)?;
    }
    let old_v = toml::Value::try_from(current)?;
    let new_v = toml::Value::try_from(&new)?;
    let (mut a, mut b) = (BTreeMap::new(), BTreeMap::new());
    flatten("", &old_v, &mut a);
    flatten("", &new_v, &mut b);
    let keys: std::collections::BTreeSet<&String> = a.keys().chain(b.keys()).collect();
    let mut merged = old_v.clone();
    let (mut applied, mut not_applied) = (Vec::new(), Vec::new());
    for k in keys {
        let (x, y) = (a.get(k), b.get(k));
        if x == y {
            continue;
        }
        let line = format!("{k}: {} -> {}", show(x), show(y));
        // the scheme switches only between upwind and linear upwind
        let central = k == "solver.momentum.scheme"
            && [x, y]
                .iter()
                .any(|v| v.and_then(|v| v.as_str()) == Some("central"));
        if RUNTIME_KEYS.contains(&k.as_str()) && !central {
            set(&mut merged, k, y);
            applied.push(line);
        } else {
            not_applied.push(line);
        }
    }
    // the scalar and energy sub-tables of the boundaries are kept outside
    // the TOML of the case: a change of them is reported, not applied (the
    // boundary conditions belong to the case the run started with)
    for (name, changed) in [
        (
            "boundary energy conditions",
            format!("{:?}", current.boundary_energy) != format!("{:?}", new.boundary_energy),
        ),
        (
            "boundary scalar conditions",
            format!("{:?}", current.boundary_scalar) != format!("{:?}", new.boundary_scalar),
        ),
    ] {
        if changed {
            not_applied.push(format!("{name}: changed in the file"));
        }
    }
    let mut case: FlowCase = merged.try_into()?;
    case.source = current.source.clone();
    case.overrides = current.overrides.clone();
    case.notes = current.notes.clone();
    case.pressure_frame = current.pressure_frame.clone();
    // the fields the round trip through TOML drops (as in
    // `schedule::apply_entry`)
    case.boundary_scalar = current.boundary_scalar.clone();
    case.boundary_energy = current.boundary_energy.clone();
    case.fluid_view = current.fluid_view.clone();
    case.validate()?;
    Ok(Reread {
        case,
        applied,
        not_applied,
    })
}
