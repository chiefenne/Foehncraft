//! Iteration and study budgets of a run (D198 (c), (d)).
//!
//! - **Iteration budget** (D198 (d)): the number of outer iterations a
//!   steady run may take, shown in the estimate of `--check` and in the run
//!   record with where it comes from (`solver.max_iterations` or
//!   `[run.limits] max_iterations`). At `[run.limits] diagnose_at`
//!   (default [`DEFAULT_DIAGNOSE_AT`]) a steady run that has not converged
//!   is diagnosed from its residual history: per equation the trend over
//!   the last iterations (falling, stalled, rising, oscillating), the
//!   iteration at which the tolerance would be reached at that rate, and
//!   advice. The diagnosis goes into the run record (`"type":
//!   "diagnosis"`) and the log; with `diagnosis = "stop"` a run that is not
//!   on its way to the tolerance within its budget stops there (exit code
//!   3, restart and final output written). The diagnosis only reads the
//!   residuals: the run itself is unchanged.
//! - **Study budgets** (D198 (c)): `[run.limits] study_cell_hours`,
//!   `study_gpu_hours` and `study_output` bound what all runs of one study
//!   may spend together. The runs of a study are the output directories
//!   next to each other (the variants of `foehncraft.study.Study` under its
//!   root); their spend is kept in the ledger `study_budget.jsonl` in the
//!   parent of the output directory, one line per finished run. The run
//!   record shows the budget, the spend so far and this run's D81 estimate
//!   (`"type": "budget"`); a run whose study budget is spent does not
//!   start, and a run that spends the rest stops as at a limit.
//!
//! A cell-hour is one cell computed for one hour of wall time (cells of the
//! mesh times hours), independent of the threads; a GPU-hour is one GPU for
//! one hour. Without any of these keys nothing changes.

use crate::case::Limits;
use anyhow::Result;
use serde::Serialize;
use std::path::{Path, PathBuf};

/// Iteration of the diagnosis of a steady run when `diagnose_at` is not
/// set (the project rule: look at a steady run after 2000 to 3000
/// iterations).
pub const DEFAULT_DIAGNOSE_AT: usize = 2000;

/// A steady iteration budget above this draws a warning in the check: a
/// steady case here needs hundreds, at most a few thousand, iterations.
pub const LARGE_ITERATION_BUDGET: usize = 10_000;

/// The file of the study ledger, in the parent of the output directory.
pub const LEDGER: &str = "study_budget.jsonl";

/// The iteration budget of a steady run, as shown in the estimate and the
/// run record.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct IterationBudget {
    /// Outer iterations the run may take.
    pub iterations: usize,
    /// The key that sets it: `solver.max_iterations` or
    /// `run.limits.max_iterations`.
    pub source: &'static str,
    /// Iteration of the diagnosis (none: off, or after the budget).
    pub diagnose_at: Option<usize>,
    /// `record` or `stop`.
    pub diagnosis: &'static str,
}

impl IterationBudget {
    /// The budget of a steady case with `max_iterations` and `limits`.
    pub fn of(max_iterations: usize, limits: &Limits) -> Self {
        let (iterations, source) = match limits.max_iterations {
            Some(m) if m < max_iterations => (m, "run.limits.max_iterations"),
            _ => (max_iterations, "solver.max_iterations"),
        };
        Self {
            iterations,
            source,
            diagnose_at: diagnose_at(limits).filter(|&n| n <= iterations),
            diagnosis: if limits.stops_at_diagnosis() {
                "stop"
            } else {
                "record"
            },
        }
    }
}

/// The iteration of the diagnosis (none when switched off by 0).
pub fn diagnose_at(limits: &Limits) -> Option<usize> {
    match limits.diagnose_at {
        Some(0) => None,
        Some(n) => Some(n),
        None => Some(DEFAULT_DIAGNOSE_AT),
    }
}

/// Iterations of the window the trend is fitted over.
fn window(at: usize) -> usize {
    (at / 4).clamp(10, 500)
}

/// The residual history the diagnosis reads: log10 of each equation's
/// residual for the last iterations before the diagnosis.
#[derive(Debug, Default)]
pub struct History {
    names: Vec<String>,
    from: usize,
    at: Option<usize>,
    rows: Vec<(usize, Vec<f64>)>,
}

impl History {
    /// Collects from the iterations before `at` (none: nothing collected).
    pub fn new(names: Vec<String>, at: Option<usize>) -> Self {
        let from = at.map_or(usize::MAX, |a| a.saturating_sub(window(a)) + 1);
        Self {
            names,
            from,
            at,
            rows: Vec::new(),
        }
    }

    /// Records the residuals of iteration `it` (in the order of the names).
    pub fn push(&mut self, it: usize, residuals: &[f64]) {
        if it >= self.from && self.at.is_some_and(|a| it <= a) {
            self.rows.push((it, residuals.to_vec()));
        }
    }

    /// The iteration is the one of the diagnosis.
    pub fn due(&self, it: usize) -> bool {
        self.at == Some(it)
    }
}

/// The trend of one equation's residual.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Trend {
    /// Falls by more than half a decade per 1000 iterations.
    Falling,
    /// Neither falls nor rises.
    Stalled,
    /// Rises by more than 0.2 decades per 1000 iterations.
    Rising,
    /// Scatters by more than 0.3 decades about its trend.
    Oscillating,
}

/// One equation of the diagnosis.
#[derive(Debug, Clone, Serialize)]
pub struct EquationTrend {
    /// `momentum_x`, `continuity`, `k`, a scalar …
    pub equation: String,
    /// The residual at the diagnosis.
    pub residual: f64,
    /// Decades per 1000 iterations over the window (negative: falling).
    pub decades_per_1000: f64,
    /// Root-mean-square scatter about the trend (decades).
    pub scatter: f64,
    /// The trend.
    pub trend: Trend,
    /// Iteration at which the residual reaches the tolerance at this rate
    /// (falling only).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tolerance_at: Option<usize>,
}

/// The verdict of a diagnosis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Verdict {
    /// Every equation falls fast enough to reach the tolerance within the
    /// budget.
    OnTrack,
    /// Falls, but reaches the tolerance only after the budget.
    Slow,
    /// The controlling equation stalls.
    Stalled,
    /// An equation rises.
    Rising,
    /// The controlling equation oscillates.
    Oscillating,
}

impl Verdict {
    /// The run is not on its way to the tolerance within its budget.
    pub fn off_track(self) -> bool {
        self != Self::OnTrack
    }
}

/// The diagnosis of a steady run at 2000 to 3000 iterations.
#[derive(Debug, Clone, Serialize)]
pub struct Diagnosis {
    /// The iteration.
    pub iteration: usize,
    /// The window of the trend (iterations).
    pub window: usize,
    /// The iteration budget.
    pub budget: usize,
    /// The tolerance of the run (`solver.tolerance`).
    pub tolerance: f64,
    /// The equation with the largest residual.
    pub controlling: String,
    /// The verdict.
    pub verdict: Verdict,
    /// Every equation.
    pub equations: Vec<EquationTrend>,
    /// What to look at or change.
    pub advice: String,
}

/// Least-squares slope and scatter of `y(x)`.
fn fit(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for (a, b) in x.iter().zip(y) {
        sxy += (a - mx) * (b - my);
        sxx += (a - mx) * (a - mx);
    }
    let s = if sxx > 0.0 { sxy / sxx } else { 0.0 };
    let var = x
        .iter()
        .zip(y)
        .map(|(a, b)| {
            let r = b - (my + s * (a - mx));
            r * r
        })
        .sum::<f64>()
        / n;
    (s, var.sqrt())
}

/// Diagnoses `h` at its iteration with the tolerance and budget of the run;
/// none before enough iterations are recorded.
pub fn diagnose(h: &History, tolerance: f64, budget: usize) -> Option<Diagnosis> {
    let at = h.at?;
    if h.rows.len() < 5 {
        return None;
    }
    let x: Vec<f64> = h.rows.iter().map(|(i, _)| *i as f64).collect();
    let mut eqs = Vec::new();
    for (k, name) in h.names.iter().enumerate() {
        let vals: Vec<(f64, f64)> = h
            .rows
            .iter()
            .zip(&x)
            .filter_map(|((_, r), xi)| {
                let v = *r.get(k)?;
                (v.is_finite() && v > 0.0).then(|| (*xi, v.log10()))
            })
            .collect();
        // an equation not solved (w in 2D) has no residual
        if vals.len() < 5 {
            continue;
        }
        let (xs, ys): (Vec<f64>, Vec<f64>) = vals.into_iter().unzip();
        let (s, scatter) = fit(&xs, &ys);
        let last = h.rows.last().and_then(|(_, r)| r.get(k).copied()).unwrap_or(0.0);
        let per_1000 = 1000.0 * s;
        let trend = if per_1000 < -0.5 {
            Trend::Falling
        } else if per_1000 > 0.2 {
            Trend::Rising
        } else if scatter > 0.3 {
            Trend::Oscillating
        } else {
            Trend::Stalled
        };
        let tolerance_at = (trend == Trend::Falling && last > tolerance && s < 0.0).then(|| {
            let more = (tolerance.log10() - last.log10()) / s;
            at + more.ceil().max(0.0) as usize
        });
        eqs.push(EquationTrend {
            equation: name.clone(),
            residual: last,
            decades_per_1000: per_1000,
            scatter,
            trend,
            tolerance_at,
        });
    }
    let ctl = eqs
        .iter()
        .max_by(|a, b| a.residual.total_cmp(&b.residual))?
        .clone();
    let rising = eqs.iter().find(|e| e.trend == Trend::Rising);
    let verdict = if rising.is_some() {
        Verdict::Rising
    } else {
        match ctl.trend {
            Trend::Stalled => Verdict::Stalled,
            Trend::Oscillating => Verdict::Oscillating,
            Trend::Rising => Verdict::Rising,
            Trend::Falling => {
                let worst = eqs
                    .iter()
                    .filter(|e| e.residual > tolerance)
                    .map(|e| e.tolerance_at.unwrap_or(usize::MAX))
                    .max()
                    .unwrap_or(at);
                if worst <= budget {
                    Verdict::OnTrack
                } else {
                    Verdict::Slow
                }
            }
        }
    };
    let sci = |v: f64| format!("{v:.1e}");
    let advice = match verdict {
        Verdict::OnTrack => format!(
            "on track: every residual falls and reaches the tolerance {} within the budget of {budget} iterations",
            sci(tolerance)
        ),
        Verdict::Slow => format!(
            "falls, but at {:.2} decades per 1000 iterations {} reaches the tolerance {} only at about iteration {}, after the budget of {budget}. \
             The residual is no stop measure: if the monitored quantities have settled, the run is done; otherwise raise the momentum \
             relaxation after the start (a [[schedule]] ramp) or use the coupled algorithm",
            -ctl.decades_per_1000,
            ctl.equation,
            sci(tolerance),
            ctl.tolerance_at.map_or("?".into(), |n| n.to_string())
        ),
        Verdict::Stalled => format!(
            "{} stalls at {}: look at the monitored quantities first (a run whose quantities have settled is done). A stall at a level comes \
             from a relaxation set too high, from limiters switching, or from mesh quality; lower the momentum relaxation, check the mesh \
             (non-orthogonality, skewness) and the scheme before running longer",
            ctl.equation,
            sci(ctl.residual)
        ),
        Verdict::Rising => {
            let e = rising.unwrap_or(&ctl);
            format!(
                "{} rises ({:+.2} decades per 1000 iterations, now {}): lower the relaxation, check the boundary conditions and the mesh; \
                 stop before it diverges",
                e.equation,
                e.decades_per_1000,
                sci(e.residual)
            )
        }
        Verdict::Oscillating => format!(
            "{} oscillates by {:.2} decades about {}: lower the relaxation; a flow that keeps cycling may be unsteady (run it transient)",
            ctl.equation,
            ctl.scatter,
            sci(ctl.residual)
        ),
    };
    Some(Diagnosis {
        iteration: at,
        window: h.rows.len(),
        budget,
        tolerance,
        controlling: ctl.equation,
        verdict,
        equations: eqs,
        advice,
    })
}

// study budgets

/// One finished run in the study ledger.
#[derive(Debug, Clone, Default, Serialize, serde::Deserialize, PartialEq)]
pub struct Spend {
    /// Cell-hours.
    pub cell_hours: f64,
    /// GPU-hours.
    pub gpu_hours: f64,
    /// Bytes written into the output directory.
    pub output_bytes: f64,
}

impl Spend {
    fn add(&mut self, o: &Self) {
        self.cell_hours += o.cell_hours;
        self.gpu_hours += o.gpu_hours;
        self.output_bytes += o.output_bytes;
    }
}

/// The ledger of the study of the output directory `out`.
pub fn ledger_path(out: &Path) -> PathBuf {
    let abs = if out.is_absolute() {
        out.to_path_buf()
    } else {
        std::env::current_dir().unwrap_or_default().join(out)
    };
    abs.parent()
        .map_or_else(|| PathBuf::from(LEDGER), |p| p.join(LEDGER))
}

/// The spend of the earlier runs of the study, without the runs of `out`
/// itself when `skip` (a run that writes its own line later).
pub fn spent(ledger: &Path) -> Spend {
    let mut s = Spend::default();
    if let Ok(text) = std::fs::read_to_string(ledger) {
        for l in text.lines() {
            if let Ok(v) = serde_json::from_str::<serde_json::Value>(l)
                && let Ok(x) = serde_json::from_value::<Spend>(v["spend"].clone())
            {
                s.add(&x);
            }
        }
    }
    s
}

/// Appends the spend of a finished run to the ledger.
pub fn append(ledger: &Path, run: &str, dir: &Path, spend: &Spend) -> Result<()> {
    use std::io::Write as _;
    let line = serde_json::to_string(&serde_json::json!({
        "run": run,
        "directory": dir.file_name().map(|n| n.to_string_lossy().into_owned()),
        "date": crate::engine::run_date(),
        "spend": spend,
    }))?;
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(ledger)?;
    // one write per line: lines of runs in parallel do not interleave
    f.write_all(format!("{line}\n").as_bytes())?;
    Ok(())
}

/// One quantity of a study budget in the run record.
#[derive(Debug, Clone, Serialize)]
pub struct BudgetLine {
    /// The budget.
    pub budget: f64,
    /// Spent by the earlier runs of the study.
    pub spent: f64,
    /// Left before this run.
    pub remaining: f64,
    /// This run's D81 estimate (none without a benchmark of the hardware).
    pub estimate: Option<f64>,
    /// The estimate fits into what is left (none without an estimate).
    pub fits: Option<bool>,
}

impl BudgetLine {
    fn new(budget: f64, spent: f64, estimate: Option<f64>) -> Self {
        let remaining = (budget - spent).max(0.0);
        Self {
            budget,
            spent,
            remaining,
            estimate,
            fits: estimate.map(|e| e <= remaining),
        }
    }
}

/// The study budgets of a run, its estimate and the spend so far: the
/// `"budget"` row of the run record.
#[derive(Debug, Clone, Serialize)]
pub struct StudyBudget {
    /// The ledger.
    pub ledger: String,
    /// Earlier runs in the ledger.
    pub runs: usize,
    /// Cell-hours.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_hours: Option<BudgetLine>,
    /// GPU-hours.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gpu_hours: Option<BudgetLine>,
    /// Output volume (bytes).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<BudgetLine>,
}

impl StudyBudget {
    /// The budgets of `limits` against the ledger `ledger` and the estimate
    /// (cell-hours, GPU-hours, output bytes) of this run.
    pub fn new(limits: &Limits, ledger: &Path, est: (Option<f64>, Option<f64>, f64)) -> Self {
        let s = spent(ledger);
        let runs = std::fs::read_to_string(ledger)
            .map(|t| t.lines().filter(|l| !l.trim().is_empty()).count())
            .unwrap_or(0);
        Self {
            ledger: ledger.display().to_string(),
            runs,
            cell_hours: limits
                .study_cell_hours
                .map(|b| BudgetLine::new(b, s.cell_hours, est.0)),
            gpu_hours: limits
                .study_gpu_hours
                .map(|b| BudgetLine::new(b, s.gpu_hours, est.1)),
            output: limits
                .study_output
                .map(|b| BudgetLine::new(b, s.output_bytes, Some(est.2))),
        }
    }

    /// The first budget that is already spent: (key, budget, spent).
    pub fn exhausted(&self) -> Option<(&'static str, f64, f64)> {
        [
            ("study_cell_hours", &self.cell_hours),
            ("study_gpu_hours", &self.gpu_hours),
            ("study_output", &self.output),
        ]
        .into_iter()
        .find_map(|(k, l)| {
            l.as_ref()
                .filter(|l| l.spent >= l.budget)
                .map(|l| (k, l.budget, l.spent))
        })
    }

    /// The budgets whose estimate does not fit into what is left.
    pub fn over_estimate(&self) -> Vec<&'static str> {
        [
            ("study_cell_hours", &self.cell_hours),
            ("study_gpu_hours", &self.gpu_hours),
            ("study_output", &self.output),
        ]
        .into_iter()
        .filter(|(_, l)| l.as_ref().and_then(|l| l.fits) == Some(false))
        .map(|(k, _)| k)
        .collect()
    }

    /// Text lines for the log.
    pub fn lines(&self) -> Vec<String> {
        let mut v = Vec::new();
        let fmt = |what: &str, l: &BudgetLine, unit: &dyn Fn(f64) -> String| {
            format!(
                "study budget {what}: {} spent of {} by {} earlier run(s), {} left; this run about {}",
                unit(l.spent),
                unit(l.budget),
                self.runs,
                unit(l.remaining),
                l.estimate.map_or("? (no benchmark of this hardware)".into(), unit)
            )
        };
        let hours = |x: f64| format!("{x:.3} h");
        let bytes = |x: f64| {
            if x >= 1e9 {
                format!("{:.2} GB", x / 1e9)
            } else {
                format!("{:.1} MB", x / 1e6)
            }
        };
        if let Some(l) = &self.cell_hours {
            v.push(fmt("cell-hours", l, &|x| format!("{x:.3}")));
        }
        if let Some(l) = &self.gpu_hours {
            v.push(fmt("GPU-hours", l, &hours));
        }
        if let Some(l) = &self.output {
            v.push(fmt("output", l, &bytes));
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn history(f: impl Fn(usize) -> f64, at: usize) -> History {
        let mut h = History::new(vec!["continuity".into()], Some(at));
        for it in 1..=at {
            h.push(it, &[f(it)]);
        }
        h
    }

    #[test]
    fn trends_of_model_histories() {
        // 2 decades per 1000 iterations from 1e-2: 1e-6 at iteration 2000 + 0
        let d = diagnose(&history(|i| 1e-2 * 10f64.powf(-2e-3 * i as f64), 2000), 1e-8, 3000)
            .unwrap();
        assert_eq!(d.equations[0].trend, Trend::Falling);
        assert!((d.equations[0].decades_per_1000 + 2.0).abs() < 1e-6);
        // 1e-6 at 2000, 1e-8 two decades later: iteration 3000
        let t = d.equations[0].tolerance_at.unwrap();
        assert!((2999..=3001).contains(&t), "{t}");
        assert_eq!(d.verdict, Verdict::OnTrack);
        let slow = diagnose(&history(|i| 1e-2 * 10f64.powf(-2e-3 * i as f64), 2000), 1e-8, 2500)
            .unwrap();
        assert_eq!(slow.verdict, Verdict::Slow);
        let flat = diagnose(&history(|_| 1e-4, 2000), 1e-8, 3000).unwrap();
        assert_eq!(flat.verdict, Verdict::Stalled);
        let up = diagnose(&history(|i| 1e-4 * 10f64.powf(1e-3 * i as f64), 2000), 1e-8, 3000)
            .unwrap();
        assert_eq!(up.verdict, Verdict::Rising);
        let osc = diagnose(
            &history(|i| if i % 2 == 0 { 1e-3 } else { 1e-4 }, 2000),
            1e-8,
            3000,
        )
        .unwrap();
        assert_eq!(osc.verdict, Verdict::Oscillating);
        assert_eq!(osc.window, 500);
    }

    #[test]
    fn unsolved_equations_are_left_out() {
        let mut h = History::new(vec!["momentum_z".into(), "continuity".into()], Some(100));
        for it in 1..=100 {
            h.push(it, &[0.0, 1e-3]);
        }
        let d = diagnose(&h, 1e-8, 200).unwrap();
        assert_eq!(d.equations.len(), 1);
        assert_eq!(d.controlling, "continuity");
    }

    #[test]
    fn iteration_budget_source() {
        let mut l = Limits::default();
        let b = IterationBudget::of(3000, &l);
        assert_eq!(b.iterations, 3000);
        assert_eq!(b.source, "solver.max_iterations");
        assert_eq!(b.diagnose_at, Some(DEFAULT_DIAGNOSE_AT));
        l.max_iterations = Some(1500);
        let b = IterationBudget::of(3000, &l);
        assert_eq!((b.iterations, b.source), (1500, "run.limits.max_iterations"));
        // the diagnosis after the budget does not happen
        assert_eq!(b.diagnose_at, None);
        l.diagnose_at = Some(0);
        assert_eq!(diagnose_at(&l), None);
    }

    #[test]
    fn ledger_sums_runs() {
        let d = std::env::temp_dir().join(format!("foehn_ledger_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let led = d.join(LEDGER);
        let _ = std::fs::remove_file(&led);
        let one = Spend {
            cell_hours: 1.5,
            gpu_hours: 0.25,
            output_bytes: 1e6,
        };
        append(&led, "a", &d.join("v0000"), &one).unwrap();
        append(&led, "b", &d.join("v0001"), &one).unwrap();
        let s = spent(&led);
        assert!((s.cell_hours - 3.0).abs() < 1e-12);
        let l = Limits {
            study_cell_hours: Some(3.0),
            study_output: Some(1e7),
            ..Limits::default()
        };
        let b = StudyBudget::new(&l, &led, (Some(0.1), None, 1e6));
        assert_eq!(b.runs, 2);
        assert_eq!(b.exhausted().map(|e| e.0), Some("study_cell_hours"));
        assert!(b.gpu_hours.is_none());
        assert_eq!(b.output.as_ref().unwrap().fits, Some(true));
        assert_eq!(ledger_path(&d.join("v0002")), led);
        std::fs::remove_dir_all(&d).unwrap();
    }
}
