//! The frozen-flow march of a transient conjugate run (D118, stage 7
//! specification §9f).
//!
//! A steady cooling flow settles within a few flow-through times, the
//! solids it cools warm up over minutes or hours: fully coupled, almost all
//! of the run's steps go to a flow that no longer changes. The run steps
//! fully coupled with the `[time]` step until the flow has settled, then
//! holds the velocity, pressure, mass flux and turbulence and marches only
//! the energy equation of fluid and solids, on the same matrix (D158), with
//! steps of `time_step` ([`crate::transient::Transient::energy_step`]). The
//! history stays physical: the fluid temperatures follow the solids at
//! every step (their own time scale is that of the flow, so a long step
//! makes them quasi-steady, which they are), and the flow is the settled
//! flow of the current temperatures.
//!
//! - **Settling** (a coupled phase). Every `window_steps` steps the
//!   velocity is compared with that of the window before: c = |u − u₋|/|u|
//!   (rms over the cells). The flow has settled when c and the change still
//!   to come, extrapolated geometrically from the last two windows (Aitken's
//!   Δ², c r/(1 − r) with r = c/c₋), are both below `settle_tolerance`. A
//!   flow that follows the temperature drifts with it while the solids
//!   warm (r stays near 1 and the extrapolation sees a change without end):
//!   for it c alone counts, and the refreshes follow the drift. `start` = a
//!   time holds the flow first at that time instead.
//! - **Refresh.** A flow that depends on the temperature (the viscosity
//!   μ(T) of its material, the Boussinesq buoyancy) is released for a new
//!   coupled phase when a fluid temperature has moved by more than
//!   `refresh_change` since the flow was held; it is held again when it
//!   has settled again. `refresh_every` frozen steps force a refresh as well
//!   (also for a flow that does not depend on T: a boundary value that
//!   changes slowly, FRZ004). A flow that does not depend on T is held once,
//!   and exactly.
//! - **Time levels.** The first frozen step after a coupled phase is BDF1
//!   (the ratio of the long step to the short one is beyond what BDF2
//!   takes); later ones are BDF2 as the case asks. The first coupled step
//!   of a refresh sees a step ratio far below one, which makes it
//!   BDF1-like.
//!
//! Ledger: every frozen step closes its heat balance as a coupled step does
//! (`heat.csv`), with the held mass flux in the convection. Outputs: the
//! `frozen_flow` section of the result (phases, steps, refreshes), the log
//! lines and run-record events (FRZ003).

use crate::diag::Diagnostic;
use crate::energy::FrozenFlow;
use crate::simplec::Flow;
use crate::transient::StepReport;
use anyhow::Result;
use serde::{Deserialize, Serialize};

/// Coupled windows without settling after which the end of the run warns
/// (FRZ002).
const UNSETTLED_WINDOWS: u32 = 20;

/// A phase of the run in the result.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct FrozenPhase {
    /// `"coupled"` or `"frozen"`.
    pub kind: String,
    /// The time at its start (s).
    pub start: f64,
    /// The time at its end (s).
    pub end: f64,
    /// Its time steps.
    pub steps: u32,
    /// Why it started: `"start"`, `"settled"`, `"start time"`,
    /// `"temperature"` (a fluid temperature moved by more than
    /// `refresh_change`) or `"refresh_every"`.
    pub reason: String,
}

/// The `frozen_flow` section of the transient summary (D118).
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct FrozenSummary {
    /// The step of the march (s).
    pub time_step: f64,
    /// True when the flow depends on the temperature (μ(T), the Boussinesq
    /// buoyancy): it is refreshed by `refresh_change`.
    pub flow_follows_temperature: bool,
    /// The time at which the flow was first held (s).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub frozen_at: Option<f64>,
    /// Steps of the march with the flow held.
    pub frozen_steps: u32,
    /// Time of the march with the flow held (s).
    pub frozen_time: f64,
    /// Steps fully coupled (the settling and the refreshes).
    pub coupled_steps: u32,
    /// Time fully coupled (s).
    pub coupled_time: f64,
    /// Refreshes of the flow.
    pub refreshes: u32,
    /// Mean solves of the energy equation per frozen step.
    pub solves_mean: f64,
    /// The largest change of a fluid temperature since the flow was last
    /// held, at the end of the run (K).
    pub fluid_change: f64,
    /// The phases of the run in order.
    pub phases: Vec<FrozenPhase>,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
enum Phase {
    Coupled,
    Frozen,
}

/// The state for a restart file (JSON).
#[derive(Debug, Clone, Serialize, Deserialize)]
struct State {
    phase: Phase,
    w_steps: u32,
    w_ref: Option<Vec<f64>>,
    w_change: Option<f64>,
    windows: u32,
    f_steps: u32,
    t_ref: Vec<f64>,
    solves: u64,
    sum: FrozenSummary,
}

/// The frozen-flow march during a run.
pub struct Controller {
    cfg: FrozenFlow,
    follows: bool,
    phase: Phase,
    /// The settling test: steps of the current window, the velocity at the
    /// end of the window before (x, y, z per cell), the change of the
    /// window before, windows of the current coupled phase.
    w_steps: u32,
    w_ref: Option<Vec<f64>>,
    w_change: Option<f64>,
    windows: u32,
    /// Frozen steps since the flow was held, and the fluid temperatures
    /// then (the fluid cells of the energy equation in order).
    f_steps: u32,
    t_ref: Vec<f64>,
    solves: u64,
    sum: FrozenSummary,
    /// Diagnostics for the run's warnings (FRZ002, FRZ004).
    pub diagnostics: Vec<Diagnostic>,
}

impl Controller {
    /// The march of a conjugate case from time `t0`; `time_dependent`:
    /// the flow's boundary values or body force depend on time (FRZ004).
    pub fn new(
        fl: &mut Flow,
        cfg: &FrozenFlow,
        buoyancy: bool,
        time_dependent: bool,
        t0: f64,
    ) -> Result<Self> {
        let (sc, _) = crate::cht::energy_fields(fl).ok_or_else(|| {
            anyhow::anyhow!("[energy.frozen_flow] needs the energy equation of a conjugate case")
        })?;
        let follows = buoyancy || sc.energy_viscosity_follows();
        let mut diagnostics = Vec::new();
        if time_dependent {
            diagnostics.push(Diagnostic::new(
                crate::diag::FRZ_TIME,
                "[energy.frozen_flow]: a boundary value or the body force of the flow depends on \
                 time; the held flow does not follow it between refreshes (set refresh_every, or \
                 let the change end before the flow settles)"
                    .to_string(),
            ));
        }
        Ok(Self {
            cfg: cfg.clone(),
            follows,
            phase: Phase::Coupled,
            w_steps: 0,
            w_ref: None,
            w_change: None,
            windows: 0,
            f_steps: 0,
            t_ref: Vec::new(),
            solves: 0,
            sum: FrozenSummary {
                time_step: cfg.time_step,
                flow_follows_temperature: follows,
                phases: vec![FrozenPhase {
                    kind: "coupled".into(),
                    start: t0,
                    end: t0,
                    steps: 0,
                    reason: "start".into(),
                }],
                ..Default::default()
            },
            diagnostics,
        })
    }

    /// True while the flow is held: the next step is an energy-only step
    /// of [`time_step`](Self::time_step).
    pub fn frozen(&self) -> bool {
        self.phase == Phase::Frozen
    }

    /// The step of the march (s).
    pub fn time_step(&self) -> f64 {
        self.cfg.time_step
    }

    /// The solves of the energy equation per frozen step at most.
    pub fn iterations(&self) -> u32 {
        self.cfg.iterations
    }

    /// True for the first frozen step after a coupled phase (BDF1).
    pub fn restart_levels(&self) -> bool {
        self.phase == Phase::Frozen && self.f_steps == 0
    }

    /// The lines of the start for the log.
    pub fn start_lines(&self) -> Vec<String> {
        let mut v = vec![format!(
            "frozen flow: coupled steps until the flow has settled ({}), then the energy marches \
             alone with steps of {} s",
            match self.cfg.start_time() {
                Some(t) => format!("held at t = {t} s"),
                None => format!(
                    "windows of {} steps, settle_tolerance {:e}",
                    self.cfg.window_steps, self.cfg.settle_tolerance
                ),
            },
            self.cfg.time_step
        )];
        v.push(if self.follows {
            format!(
                "frozen flow: the flow depends on the temperature; refreshed when a fluid \
                 temperature has moved by more than {} K{}",
                self.cfg.refresh_change,
                if self.cfg.refresh_every > 0 {
                    format!(" or every {} frozen steps", self.cfg.refresh_every)
                } else {
                    String::new()
                }
            )
        } else if self.cfg.refresh_every > 0 {
            format!(
                "frozen flow: the flow does not depend on the temperature; refreshed every {} \
                 frozen steps",
                self.cfg.refresh_every
            )
        } else {
            "frozen flow: the flow does not depend on the temperature; it is held once (exact)"
                .to_string()
        });
        v
    }

    /// After an accepted step `rp` (`frozen`: an energy-only step):
    /// counts it, runs the settling test of a coupled phase or the refresh
    /// test of a frozen one. Returns the lines for the log.
    pub fn after_step(
        &mut self,
        fl: &mut Flow,
        rp: &StepReport,
        frozen: bool,
    ) -> Result<Vec<String>> {
        let mut lines = Vec::new();
        let (t, dt) = (rp.time, rp.dt);
        if let Some(p) = self.sum.phases.last_mut() {
            p.end = t;
            p.steps += 1;
        }
        if frozen {
            self.sum.frozen_steps += 1;
            self.sum.frozen_time += dt;
            self.f_steps += 1;
            self.solves += rp.outer as u64;
            let change = self.fluid_change(fl)?;
            self.sum.fluid_change = change;
            if self.follows && change > self.cfg.refresh_change {
                self.release(
                    fl,
                    t,
                    "temperature",
                    format!(
                        "a fluid temperature has moved by {change:.4} K (> refresh_change {} K) \
                         since the flow was held",
                        self.cfg.refresh_change
                    ),
                    &mut lines,
                )?;
            } else if self.cfg.refresh_every > 0 && self.f_steps >= self.cfg.refresh_every {
                self.release(
                    fl,
                    t,
                    "refresh_every",
                    format!("{} frozen steps (refresh_every)", self.f_steps),
                    &mut lines,
                )?;
            }
            return Ok(lines);
        }
        self.sum.coupled_steps += 1;
        self.sum.coupled_time += dt;
        let first = self.sum.frozen_at.is_none();
        if first && let Some(ts) = self.cfg.start_time() {
            if t >= ts - 1e-6 * dt {
                self.hold(
                    fl,
                    t,
                    "start time",
                    format!("the start time {ts} s"),
                    &mut lines,
                )?;
            }
            return Ok(lines);
        }
        self.w_steps += 1;
        if self.w_steps < self.cfg.window_steps {
            return Ok(lines);
        }
        self.w_steps = 0;
        self.windows += 1;
        let u = velocity(fl)?;
        let Some(r) = self.w_ref.replace(u.clone()) else {
            return Ok(lines);
        };
        let (mut dd, mut uu) = (0.0, 0.0);
        for (a, b) in u.iter().zip(&r) {
            dd += (a - b) * (a - b);
            uu += a * a;
        }
        let c = if uu > 0.0 {
            (dd / uu).sqrt()
        } else {
            f64::INFINITY
        };
        let prev = self.w_change.replace(c);
        let Some(c0) = prev else {
            return Ok(lines);
        };
        let ratio = if c0 > 0.0 { c / c0 } else { 0.0 };
        let rest = if c == 0.0 {
            0.0
        } else if (0.0..1.0).contains(&ratio) {
            c * ratio / (1.0 - ratio)
        } else {
            f64::INFINITY
        };
        let tol = self.cfg.settle_tolerance;
        // a flow that follows the temperature drifts with it (the window
        // changes then hardly fall, and the extrapolation sees a change
        // without end): for it the change of the last window counts, the
        // drift is what the refreshes follow
        if c <= tol && (rest <= tol || self.follows) {
            self.hold(
                fl,
                t,
                "settled",
                format!(
                    "the flow has settled (velocity change {c:.2e} over the last window, \
                     {rest:.2e} still to come; settle_tolerance {tol:e})"
                ),
                &mut lines,
            )?;
        }
        Ok(lines)
    }

    /// Holds the flow at time `t`.
    fn hold(
        &mut self,
        fl: &mut Flow,
        t: f64,
        reason: &str,
        why: String,
        lines: &mut Vec<String>,
    ) -> Result<()> {
        self.phase = Phase::Frozen;
        self.f_steps = 0;
        self.t_ref = fluid_temperatures(fl)?;
        self.sum.fluid_change = 0.0;
        if self.sum.frozen_at.is_none() {
            self.sum.frozen_at = Some(t);
        }
        self.push_phase("frozen", t, reason);
        lines.push(format!(
            "frozen flow: t = {t:.6e} s, {why}: the flow is held, the energy marches with steps \
             of {} s",
            self.cfg.time_step
        ));
        Ok(())
    }

    /// Releases the flow at time `t` for a refresh.
    fn release(
        &mut self,
        fl: &mut Flow,
        t: f64,
        reason: &str,
        why: String,
        lines: &mut Vec<String>,
    ) -> Result<()> {
        self.phase = Phase::Coupled;
        self.w_steps = 0;
        self.w_ref = Some(velocity(fl)?);
        self.w_change = None;
        self.windows = 0;
        self.sum.refreshes += 1;
        self.push_phase("coupled", t, reason);
        lines.push(format!(
            "frozen flow: t = {t:.6e} s, {why}: refresh {} of the flow, coupled steps until it \
             has settled again",
            self.sum.refreshes
        ));
        Ok(())
    }

    fn push_phase(&mut self, kind: &str, t: f64, reason: &str) {
        self.sum.phases.push(FrozenPhase {
            kind: kind.into(),
            start: t,
            end: t,
            steps: 0,
            reason: reason.into(),
        });
    }

    /// The largest change of a fluid temperature since the flow was held.
    fn fluid_change(&self, fl: &mut Flow) -> Result<f64> {
        let now = fluid_temperatures(fl)?;
        Ok(now
            .iter()
            .zip(&self.t_ref)
            .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs())))
    }

    /// The state for a restart file (JSON).
    pub fn save(&self) -> Result<Vec<u8>> {
        let st = State {
            phase: self.phase,
            w_steps: self.w_steps,
            w_ref: self.w_ref.clone(),
            w_change: self.w_change,
            windows: self.windows,
            f_steps: self.f_steps,
            t_ref: self.t_ref.clone(),
            solves: self.solves,
            sum: self.sum.clone(),
        };
        Ok(serde_json::to_vec(&st)?)
    }

    /// Continues from the state of a restart file.
    pub fn restore(&mut self, bytes: &[u8]) -> Result<()> {
        let st: State = serde_json::from_slice(bytes)?;
        self.phase = st.phase;
        self.w_steps = st.w_steps;
        self.w_ref = st.w_ref;
        self.w_change = st.w_change;
        self.windows = st.windows;
        self.f_steps = st.f_steps;
        self.t_ref = st.t_ref;
        self.solves = st.solves;
        let (step, follows) = (self.sum.time_step, self.sum.flow_follows_temperature);
        self.sum = st.sum;
        self.sum.time_step = step;
        self.sum.flow_follows_temperature = follows;
        Ok(())
    }

    /// The summary at the end of the run; FRZ002 when the flow never
    /// settled, or a refresh did not settle within twenty windows.
    pub fn summary(&mut self) -> FrozenSummary {
        let mut s = self.sum.clone();
        s.solves_mean = if s.frozen_steps > 0 {
            self.solves as f64 / f64::from(s.frozen_steps)
        } else {
            0.0
        };
        let unsettled = self.phase == Phase::Coupled
            && (s.frozen_at.is_none() && self.cfg.start_time().is_none()
                || self.windows >= UNSETTLED_WINDOWS);
        if unsettled {
            self.diagnostics.push(Diagnostic::new(
                crate::diag::FRZ_NOT_SETTLED,
                format!(
                    "[energy.frozen_flow]: {} did not settle by the window test ({} windows of {} \
                     steps, settle_tolerance {:e}); the run stayed fully coupled from t = {:.6e} \
                     s on (an unsteady flow, or a tolerance too tight; `start` = a time holds it \
                     anyway)",
                    if s.frozen_at.is_none() {
                        "the flow"
                    } else {
                        "the flow of the last refresh"
                    },
                    self.windows,
                    self.cfg.window_steps,
                    self.cfg.settle_tolerance,
                    s.phases.last().map_or(0.0, |p| p.start)
                ),
            ));
        }
        s
    }
}

/// The velocity of the fluid view (x, y, z of every cell).
fn velocity(fl: &mut Flow) -> Result<Vec<f64>> {
    let u = fl.f.u;
    let mut v = fl.dev.download(u.x)?;
    v.extend(fl.dev.download(u.y)?);
    v.extend(fl.dev.download(u.z)?);
    Ok(v)
}

/// The temperatures of the fluid cells of the energy equation, in the
/// order of the cells.
fn fluid_temperatures(fl: &mut Flow) -> Result<Vec<f64>> {
    let (sc, dev) = crate::cht::energy_fields(fl).ok_or_else(|| {
        anyhow::anyhow!("[energy.frozen_flow] needs the energy equation of a conjugate case")
    })?;
    let regions = sc.energy_cell_regions();
    let t = sc.energy_temperature(dev)?;
    Ok(t.iter()
        .zip(&regions)
        .filter(|(_, m)| **m == 0)
        .map(|(&x, _)| x)
        .collect())
}
