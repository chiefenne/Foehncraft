//! Settings that change at prepared points of a run (`[[schedule]]`, D128).
//!
//! ```toml
//! [[schedule]]
//! at_iteration = 200                       # steady; transient: at_step or at_time
//! set = { solver.momentum.relaxation = 0.9, solver.pressure.relaxation = 0.4 }
//! ```
//!
//! A value holds *from* the named point on: before the start of iteration
//! 200 (step 500, time 2 ms) the run uses the value of the case file, from
//! there on the value of the entry, until a later entry sets the key again.
//! The whole run is set up before it starts and the case file alone says
//! what was run.
//!
//! **Schedulable keys** are the keys a re-read may change during a run
//! ([`crate::control::RUNTIME_KEYS`], D40) without those that steer the run
//! itself ([`EXCLUDED`]). Everything else belongs to the state the run
//! starts from (mesh, boundary values, schemes, the case's physics) and
//! cannot change without a new start; asking for it is an error when the
//! case is read, with the reason.
//!
//! Rules checked when the case is read ([`validate`]): one trigger per
//! entry (`at_iteration` in a steady case, `at_step` or `at_time` in a
//! transient one, not both kinds in one case), the entries in ascending
//! order, the keys schedulable, and the values valid: the entries are
//! applied one after the other to a copy of the case, which is validated
//! after each of them.
//!
//! In the run ([`Runner`]): before iteration or step N the entries due are
//! applied in order, every change is a ` [I]` line in the log and an
//! `event` row (`RUN015`) in `run.jsonl` with `key`, `old` and `new`. A
//! restart applies the entries that lie before its start at once. A
//! re-read of the case file (D40) keeps the entries that have fired in
//! force over the file's values; the schedule itself cannot be re-read.

use crate::case::FlowCase;
use crate::control::{RUNTIME_KEYS, flatten, set, show};
use crate::diag;
use anyhow::{Context, Result};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Keys that a re-read may change but a schedule may not, with the reason.
pub const EXCLUDED: &[(&str, &str)] = &[
    (
        "run.reread_every",
        "it steers the re-reading itself; the schedule is the prepared form of it",
    ),
    (
        "run.limits.max_wall_time",
        "the limits are guard rails of the machine and the session, not settings of the numerics",
    ),
    (
        "run.limits.max_steps",
        "the limits are guard rails of the machine and the session, not settings of the numerics",
    ),
    (
        "run.limits.max_iterations",
        "the limits are guard rails of the machine and the session, not settings of the numerics",
    ),
    (
        "run.limits.max_output",
        "the limits are guard rails of the machine and the session, not settings of the numerics",
    ),
];

/// One entry of `[[schedule]]`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct ScheduleEntry {
    /// Steady case: the iteration from whose start the values hold (≥ 1).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub at_iteration: Option<usize>,
    /// Transient case: the time step from whose start the values hold (≥ 1).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub at_step: Option<usize>,
    /// Transient case: the time (s) from which the values hold; the run
    /// lands a step on it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub at_time: Option<f64>,
    /// The settings that change: `key = value`, with the dotted key path
    /// of the case file (`solver.momentum.relaxation = 0.9`).
    pub set: BTreeMap<String, serde_json::Value>,
}

impl ScheduleEntry {
    /// The trigger as text (`iteration 200`, `step 500`, `t = 0.002 s`).
    pub fn when(&self) -> String {
        match (self.at_iteration, self.at_step, self.at_time) {
            (Some(i), _, _) => format!("iteration {i}"),
            (_, Some(s), _) => format!("step {s}"),
            (_, _, Some(t)) => format!("t = {t} s"),
            _ => "?".into(),
        }
    }

    fn position(&self) -> f64 {
        match (self.at_iteration, self.at_step, self.at_time) {
            (Some(i), _, _) => i as f64,
            (_, Some(s), _) => s as f64,
            (_, _, Some(t)) => t,
            _ => f64::NAN,
        }
    }

    /// The settings flattened to dotted keys.
    fn flat(&self) -> Result<BTreeMap<String, toml::Value>> {
        let v = toml::Value::try_from(&self.set)?;
        let mut out = BTreeMap::new();
        flatten("", &v, &mut out);
        Ok(out)
    }
}

/// Whether a schedule may set `key`.
pub fn schedulable(key: &str) -> bool {
    RUNTIME_KEYS.contains(&key) && !EXCLUDED.iter().any(|(k, _)| *k == key)
}

/// One change made by an entry.
#[derive(Debug, Clone, PartialEq)]
pub struct Change {
    /// The entry's trigger.
    pub when: String,
    /// The key.
    pub key: String,
    /// The value before (`(default)` when the key was not set).
    pub old: String,
    /// The value after.
    pub new: String,
}

fn bad(msg: String, i: usize) -> anyhow::Error {
    diag::err(diag::Diagnostic::new(diag::SCHED_ENTRY, msg).key(format!("schedule[{i}]")))
}

/// Checks the schedule of `case` when the case is read.
pub fn validate(case: &FlowCase) -> Result<()> {
    if case.schedule.is_empty() {
        return Ok(());
    }
    let transient = case.time.is_some();
    let mut kind: Option<&str> = None;
    let mut last = f64::NEG_INFINITY;
    for (i, e) in case.schedule.iter().enumerate() {
        let n = i + 1;
        let given = usize::from(e.at_iteration.is_some())
            + usize::from(e.at_step.is_some())
            + usize::from(e.at_time.is_some());
        if given != 1 {
            return Err(bad(
                format!(
                    "[[schedule]] entry {n}: give exactly one of `at_iteration` (steady), \
                     `at_step` or `at_time` (transient)"
                ),
                n,
            ));
        }
        let k = if e.at_iteration.is_some() {
            "at_iteration"
        } else if e.at_step.is_some() {
            "at_step"
        } else {
            "at_time"
        };
        if transient == (k == "at_iteration") {
            return Err(bad(
                format!(
                    "[[schedule]] entry {n}: `{k}` does not fit a {} case; use {}",
                    if transient { "transient" } else { "steady" },
                    if transient {
                        "`at_step` or `at_time`"
                    } else {
                        "`at_iteration`"
                    }
                ),
                n,
            ));
        }
        if let Some(prev) = kind
            && prev != k
        {
            return Err(bad(
                format!(
                    "[[schedule]] entry {n}: `{k}` after `{prev}`; use one kind of trigger in a case \
                     (steps and times cannot be put in order before the run)"
                ),
                n,
            ));
        }
        kind = Some(k);
        let p = e.position();
        let start = case.time.map_or(0.0, |t| t.start_time);
        let ok = match k {
            "at_time" => p.is_finite() && p > start,
            _ => p >= 1.0,
        };
        if !ok {
            return Err(bad(
                format!(
                    "[[schedule]] entry {n}: `{k}` must be {} (the value of the case file holds \
                     before it)",
                    if k == "at_time" {
                        format!("later than the start time {start} s")
                    } else {
                        "at least 1".into()
                    }
                ),
                n,
            ));
        }
        if p < last {
            return Err(bad(
                format!(
                    "[[schedule]] entry {n}: `{k}` = {p} comes before the entry above; list the \
                     entries in ascending order"
                ),
                n,
            ));
        }
        last = p;
        if e.set.is_empty() {
            return Err(bad(
                format!("[[schedule]] entry {n}: `set` names no setting"),
                n,
            ));
        }
        for key in e.flat()?.keys() {
            check_key(key, n, transient)?;
        }
    }
    // the values: apply the entries one after the other, validating each
    let mut c = case.clone();
    for (i, e) in case.schedule.iter().enumerate() {
        c = apply_entry(&c, e, None)
            .with_context(|| format!("[[schedule]] entry {} ({})", i + 1, e.when()))?;
        scheme_switch(case, &c, i + 1)?;
    }
    Ok(())
}

/// The momentum scheme may switch between `upwind` and `linear-upwind`
/// during a run (a first-order start, then second order); the bounded
/// central scheme keeps state of its own and is chosen at the start.
pub(crate) fn scheme_switch(old: &FlowCase, new: &FlowCase, n: usize) -> Result<()> {
    use foehncraft_fv::solver::Scheme;
    let (a, b) = (old.solver.momentum.scheme, new.solver.momentum.scheme);
    if a != b && (a == Scheme::Central || b == Scheme::Central) {
        return Err(bad(
            format!(
                "[[schedule]] entry {n}: `solver.momentum.scheme` may switch only between \
                 \"upwind\" and \"linear-upwind\" during a run, not to or from \"central\""
            ),
            n,
        ));
    }
    Ok(())
}

fn check_key(key: &str, n: usize, transient: bool) -> Result<()> {
    if let Some((_, why)) = EXCLUDED.iter().find(|(k, _)| *k == key) {
        return Err(bad(
            format!("[[schedule]] entry {n}: `{key}` cannot be scheduled: {why}"),
            n,
        ));
    }
    if !RUNTIME_KEYS.contains(&key) {
        let hint = crate::bc::did_you_mean(key, RUNTIME_KEYS.iter().copied())
            .map(|c| format!("; did you mean `{c}`?"))
            .unwrap_or_default();
        return Err(bad(
            format!(
                "[[schedule]] entry {n}: `{key}` cannot change during a run (it belongs to the \
                 state the run starts from: mesh, boundaries, schemes, physics); stop the run and \
                 restart from the restart file to change it{hint}"
            ),
            n,
        ));
    }
    if key.starts_with("time.") && !transient {
        return Err(bad(
            format!("[[schedule]] entry {n}: `{key}` belongs to a transient case"),
            n,
        ));
    }
    Ok(())
}

/// The case with the settings of `e` applied and validated; `changes`
/// receives what changed.
pub fn apply_entry(
    case: &FlowCase,
    e: &ScheduleEntry,
    changes: Option<&mut Vec<Change>>,
) -> Result<FlowCase> {
    let old_v = toml::Value::try_from(case)?;
    let mut old_flat = BTreeMap::new();
    flatten("", &old_v, &mut old_flat);
    let mut merged = old_v;
    let mut list = Vec::new();
    for (k, v) in e.flat()? {
        let before = old_flat.get(&k);
        if before != Some(&v) {
            list.push(Change {
                when: e.when(),
                key: k.clone(),
                old: show(before),
                new: show(Some(&v)),
            });
        }
        // a pressure value in the file's convention (absolute, D239) is
        // set as the solver's gauge value
        let v = case
            .pressure_frame
            .as_ref()
            .and_then(|f| f.gauge_value(&k, &v))
            .unwrap_or(v);
        set(&mut merged, &k, Some(&v));
    }
    let mut c: FlowCase = merged
        .try_into()
        .map_err(|err: toml::de::Error| bad(format!("a value does not fit its key: {err}"), 0))?;
    // the fields the round trip through TOML drops: the boundary
    // sub-tables, the fluid view, the pressure convention, ...
    c.keep_reader_fields(case);
    // the schedule was checked with the case; validating the setting values
    // must not check it again
    let sched = std::mem::take(&mut c.schedule);
    c.validate()?;
    c.schedule = sched;
    if let Some(ch) = changes {
        ch.extend(list);
    }
    Ok(c)
}

/// `case` with the first `fired` entries of its own schedule applied: the
/// settings in force after a re-read of the case file.
pub fn apply_fired(case: &FlowCase, base: &FlowCase, fired: usize) -> Result<FlowCase> {
    let mut c = case.clone();
    for e in &base.schedule[..fired.min(base.schedule.len())] {
        c = apply_entry(&c, e, None)?;
    }
    Ok(c)
}

/// Where the run is: the iteration or step about to start, and its time.
#[derive(Debug, Clone, Copy)]
pub struct Position {
    /// The iteration or step about to start (1-based).
    pub n: usize,
    /// Its start time (transient; 0 in a steady case).
    pub time: f64,
    /// Tolerance on times (s).
    pub tol: f64,
}

/// The progress of a run through its schedule.
#[derive(Debug, Default, Clone)]
pub struct Runner {
    /// Entries applied so far.
    pub fired: usize,
}

impl Runner {
    fn due(e: &ScheduleEntry, p: Position) -> bool {
        match (e.at_iteration, e.at_step, e.at_time) {
            (Some(i), _, _) => i <= p.n,
            (_, Some(s), _) => s <= p.n,
            (_, _, Some(t)) => t <= p.time + p.tol,
            _ => false,
        }
    }

    /// The case after the entries that are due at `p`, with what changed;
    /// `None` when none is due.
    pub fn fire(
        &mut self,
        case: &FlowCase,
        p: Position,
    ) -> Result<Option<(FlowCase, Vec<Change>)>> {
        let mut c: Option<FlowCase> = None;
        let mut changes = Vec::new();
        while let Some(e) = case.schedule.get(self.fired)
            && Self::due(e, p)
        {
            let base = c.as_ref().unwrap_or(case);
            c = Some(apply_entry(base, e, Some(&mut changes))?);
            self.fired += 1;
        }
        Ok(c.map(|c| (c, changes)))
    }

    /// The next time (s) a not yet fired `at_time` entry lies at or after
    /// `t`: the run lands a step on it.
    pub fn next_time(&self, case: &FlowCase, t: f64, tol: f64) -> Option<f64> {
        case.schedule[self.fired.min(case.schedule.len())..]
            .iter()
            .find_map(|e| e.at_time)
            .filter(|&x| x > t + tol)
    }
}

/// Rows of the "Run control" frame: `(when, settings)` per entry.
pub fn rows(case: &FlowCase) -> Vec<(String, String)> {
    case.schedule
        .iter()
        .map(|e| {
            let list = e
                .flat()
                .map(|f| {
                    f.iter()
                        .map(|(k, v)| format!("{k} = {v}"))
                        .collect::<Vec<_>>()
                        .join("\n")
                })
                .unwrap_or_default();
            (e.when(), list)
        })
        .collect()
}
