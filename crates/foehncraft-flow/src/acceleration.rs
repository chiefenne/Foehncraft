//! Thermal acceleration of solid regions in transient conjugate heat
//! transfer (D213, stage 7 specification §9e).
//!
//! The heat capacity ρ c_p of the accelerated solids is multiplied by a
//! factor 0 < f ≤ 1 in the time term and the heat storage of the energy
//! equation (never in the fluid, where ρ c_p also multiplies the
//! convection). The steady energy equation of a solid has no ρ c_p in it,
//! so the steady temperatures and the cycle mean of a periodic flow are
//! those of the true solid; only the way there is f times shorter. When
//! the solids have settled (or at a set time) the true ρ c_p comes back,
//! at once or raised geometrically over a few windows (a periodic flow:
//! the solid that the jump would catch in the middle of its swing keeps
//! that part in its deep temperature, which then decays only with the true
//! time constant; a ramp lets the swing die out on the way).
//!
//! With `inertia_scale = "auto"` the factor is set after the first window
//! from the time constant τ = C/G of each accelerated region: C = Σ ρ c_p V
//! its true heat capacity, G its conductance to the rest: the change of
//! the heat flow out of the region through its interfaces for a rise of
//! its temperature by 1 K, with the fluid quasi-steady (two steady solves
//! of the energy equation with the solids held), plus the conductance of
//! its outer boundaries with a temperature or a convection condition. The
//! internal conduction of the solid is not in τ (a small Biot number);
//! where it matters, τ is underestimated and the acceleration is milder
//! than asked for (never stronger).
//!
//! Ledger: the equation stores f times the heat of the true solid; the
//! rest, (1 − f) of the true storage, is booked per region as
//! `inertia_scaling` (W per step, J over the run), so that the balance of
//! the solved equation closes as before and the heat the acceleration
//! skipped is in the record.

use crate::energy::Acceleration;
use crate::simplec::Flow;
use anyhow::Result;
use serde::{Deserialize, Serialize};

/// The smallest factor of `"auto"`.
pub const FACTOR_MIN: f64 = 1e-4;

/// Steady solves per state of the conductance measurement (at most; it
/// stops when the heat flows have settled to 10⁻⁶).
const PROBE_ITERATIONS: usize = 40;

/// One accelerated region in the result.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct RegionAcceleration {
    /// The solid region.
    pub region: String,
    /// The factor f on ρ c_p while accelerated.
    pub factor: f64,
    /// The true heat capacity Σ ρ c_p V (J/K).
    pub heat_capacity: f64,
    /// `"auto"`: the measured conductance G to the rest (W/K).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conductance: Option<f64>,
    /// `"auto"`: the true time constant C/G (s).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub time_constant: Option<f64>,
    /// The heat the true ρ c_p would have stored in addition over the run
    /// (J): Σ (1 − f) × true storage × Δt.
    pub inertia_scaling_heat: f64,
    /// With a `period`: the skin depth √(2k/(f ρ c_p ω)) while accelerated
    /// (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skin_depth: Option<f64>,
    /// The thickness V/A of the region (m), A its interface area.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thickness: Option<f64>,
}

/// The `acceleration` section of the transient summary (D213).
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct AccelerationSummary {
    /// `"factor"` or `"auto"`.
    pub mode: String,
    /// The length of a window (s; one period, or `window_steps` steps of
    /// the last window).
    pub window: f64,
    /// Windows completed.
    pub windows: u32,
    /// The accelerated regions.
    pub regions: Vec<RegionAcceleration>,
    /// The time from which the solids were accelerated (s).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub accelerated_from: Option<f64>,
    /// The time at which the solids were found settled (s).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub settled_at: Option<f64>,
    /// The time at which the restore started (s).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub restore_started: Option<f64>,
    /// The time from which the solids have the true ρ c_p (s).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub restored_at: Option<f64>,
    /// The time from which the history counts as physical (s):
    /// `restored_at` plus `settle_windows` windows.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub physical_from: Option<f64>,
    /// True when the run went past `physical_from`.
    pub history_physical: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Phase {
    /// `"auto"`: the first window with the true ρ c_p, at whose end the
    /// conductances are measured.
    Measure,
    Accelerated,
    /// The factor rises geometrically from the time and step it started.
    Ramp {
        time: f64,
        steps: u32,
    },
    Restored,
}

/// The thermal acceleration of a transient conjugate run.
#[derive(Debug, Clone)]
pub struct Controller {
    cfg: Acceleration,
    /// The regions of the energy equation (`fluid` first).
    names: Vec<String>,
    /// The accelerated regions (indices into `names`).
    idx: Vec<usize>,
    /// The factor of each accelerated region while accelerated.
    f0: Vec<f64>,
    phase: Phase,
    /// The current window: steps, time, Σ T̄ Δt per accelerated region.
    w_steps: u32,
    w_time: f64,
    w_sum: Vec<f64>,
    /// The end of the current window (a period).
    w_end: f64,
    /// Accelerated windows completed.
    acc_windows: u32,
    /// The mean temperatures at the start and of the last window.
    start: Vec<f64>,
    prev: Option<Vec<f64>>,
    /// The last window length (s).
    window: f64,
    sum: AccelerationSummary,
}

impl Controller {
    /// The controller of a case with `[energy.acceleration]`, `None`
    /// without; a fixed factor acts from the first step.
    pub fn new(fl: &mut Flow, cfg: &Acceleration, t0: f64) -> Result<Option<Self>> {
        let Some((sc, dev)) = crate::cht::energy_fields(fl) else {
            return Ok(None);
        };
        let names = sc.energy_regions();
        let idx: Vec<usize> = if cfg.regions.is_empty() {
            (1..names.len()).collect()
        } else {
            cfg.regions
                .iter()
                .filter_map(|r| names.iter().position(|n| n == r))
                .collect()
        };
        let state = sc.energy_region_state(dev)?;
        let start: Vec<f64> = idx.iter().map(|&i| state[i].2).collect();
        let fixed = cfg.factor();
        let mut c = Controller {
            cfg: cfg.clone(),
            f0: vec![fixed.unwrap_or(1.0); idx.len()],
            phase: if fixed.is_some() {
                Phase::Accelerated
            } else {
                Phase::Measure
            },
            w_steps: 0,
            w_time: 0.0,
            w_sum: vec![0.0; idx.len()],
            w_end: t0 + cfg.period.unwrap_or(0.0),
            acc_windows: 0,
            start,
            prev: None,
            window: cfg.period.unwrap_or(0.0),
            sum: AccelerationSummary {
                mode: if fixed.is_some() { "factor" } else { "auto" }.into(),
                regions: idx
                    .iter()
                    .map(|&i| RegionAcceleration {
                        region: names[i].clone(),
                        factor: fixed.unwrap_or(1.0),
                        heat_capacity: state[i].1,
                        ..Default::default()
                    })
                    .collect(),
                ..Default::default()
            },
            names,
            idx,
        };
        if fixed.is_some() {
            c.sum.accelerated_from = Some(t0);
            c.apply(fl, &c.f0.clone());
        }
        Ok(Some(c))
    }

    /// The factor of each region of the energy equation for factors `f`
    /// of the accelerated ones.
    fn apply(&self, fl: &mut Flow, f: &[f64]) {
        if let Some((sc, _)) = crate::cht::energy_fields(fl) {
            let mut all = vec![1.0; self.names.len()];
            for (&i, &x) in self.idx.iter().zip(f) {
                all[i] = x;
            }
            sc.set_energy_inertia(&all);
        }
    }

    /// True while a region has a factor below 1.
    pub fn active(&self) -> bool {
        !matches!(self.phase, Phase::Measure | Phase::Restored)
    }

    /// True once the true ρ c_p is back.
    pub fn restored(&self) -> bool {
        self.phase == Phase::Restored
    }

    /// After every accepted time step (at the new time `t`, the step `dt`):
    /// books the inertia scaling, closes windows, measures, restores.
    /// Returns the lines for the log (and the run record).
    pub fn after_step(&mut self, fl: &mut Flow, t: f64, dt: f64) -> Result<Vec<String>> {
        let mut lines = Vec::new();
        let Some((sc, dev)) = crate::cht::energy_fields(fl) else {
            return Ok(lines);
        };
        if let Some(st) = sc.energy_storage(dev)? {
            for (r, &i) in self.sum.regions.iter_mut().zip(&self.idx) {
                r.inertia_scaling_heat += st[i].inertia_scaling * dt;
            }
        }
        let state = sc.energy_region_state(dev)?;
        for (s, &i) in self.w_sum.iter_mut().zip(&self.idx) {
            *s += state[i].2 * dt;
        }
        self.w_steps += 1;
        self.w_time += dt;
        // the restore at a set time
        if self.phase == Phase::Accelerated
            && let Some(tr) = self.cfg.restore_time()
            && t >= tr - 1e-9 * dt
        {
            self.start_restore(t, &mut lines);
        }
        // the ramp
        if let Phase::Ramp { time, steps } = self.phase {
            let n = f64::from(self.cfg.ramp());
            let x = match self.cfg.period {
                Some(p) => (t - time) / (n * p),
                None => f64::from(steps + 1) / (n * f64::from(self.cfg.window_steps)),
            };
            if x >= 1.0 - 1e-9 {
                self.finish_restore(fl, t, dt, &mut lines);
            } else {
                let f: Vec<f64> = self.f0.iter().map(|&f| f.powf(1.0 - x)).collect();
                self.apply(fl, &f);
                self.phase = Phase::Ramp {
                    time,
                    steps: steps + 1,
                };
            }
        }
        // the end of a window
        let end = match self.cfg.period {
            Some(p) => t >= self.w_end - 1e-9 * p,
            None => self.w_steps >= self.cfg.window_steps,
        };
        if end {
            self.close_window(fl, t, &mut lines)?;
        }
        Ok(lines)
    }

    fn close_window(&mut self, fl: &mut Flow, t: f64, lines: &mut Vec<String>) -> Result<()> {
        let mean: Vec<f64> = self.w_sum.iter().map(|s| s / self.w_time).collect();
        self.window = match self.cfg.period {
            Some(p) => p,
            None => self.w_time,
        };
        self.sum.window = self.window;
        self.sum.windows += 1;
        if let Some(p) = self.cfg.period {
            self.w_end += p;
        }
        self.w_steps = 0;
        self.w_time = 0.0;
        self.w_sum.iter_mut().for_each(|s| *s = 0.0);
        match self.phase {
            Phase::Measure => {
                self.measure(fl, t, lines)?;
            }
            Phase::Accelerated => {
                self.acc_windows += 1;
                if self.cfg.restore_time().is_none()
                    && self.acc_windows >= 2
                    && let Some(prev) = &self.prev
                {
                    let tol = self.cfg.settle_tolerance;
                    let settled = mean.iter().zip(prev).zip(&self.start).all(|((m, p), s)| {
                        let total = (m - s).abs().max(1e-12 * m.abs());
                        (m - p).abs() <= tol * total
                    });
                    if settled {
                        self.sum.settled_at = Some(t);
                        lines.push(format!(
                            "thermal acceleration: the solids have settled at t = {t:.6e} s \
                             (window-mean temperatures change by less than {:e} of their rise)",
                            tol
                        ));
                        self.start_restore(t, lines);
                        if self.cfg.ramp() == 0 {
                            self.finish_restore(fl, t, self.window, lines);
                        }
                    }
                }
            }
            _ => {}
        }
        self.prev = Some(mean);
        Ok(())
    }

    fn start_restore(&mut self, t: f64, lines: &mut Vec<String>) {
        self.sum.restore_started = Some(t);
        self.phase = Phase::Ramp { time: t, steps: 0 };
        if self.cfg.ramp() > 0 {
            lines.push(format!(
                "thermal acceleration: restoring the true ρ c_p over {} windows from t = {t:.6e} s",
                self.cfg.ramp()
            ));
        }
    }

    fn finish_restore(&mut self, fl: &mut Flow, t: f64, dt: f64, lines: &mut Vec<String>) {
        self.apply(fl, &vec![1.0; self.idx.len()]);
        self.phase = Phase::Restored;
        self.sum.restored_at = Some(t);
        let w = if self.window > 0.0 {
            self.window
        } else {
            f64::from(self.cfg.window_steps) * dt
        };
        self.sum.physical_from = Some(t + f64::from(self.cfg.settle_windows) * w);
        lines.push(format!(
            "thermal acceleration: true ρ c_p from t = {t:.6e} s; the history counts as \
             physical from t = {:.6e} s",
            t + f64::from(self.cfg.settle_windows) * w
        ));
    }

    /// `"auto"`: the conductances, time constants and factors.
    fn measure(&mut self, fl: &mut Flow, t: f64, lines: &mut Vec<String>) -> Result<()> {
        let g = conductances(fl, &self.idx, t)?;
        let Some((sc, dev)) = crate::cht::energy_fields(fl) else {
            return Ok(());
        };
        let state = sc.energy_region_state(dev)?;
        let target = self.cfg.target() * self.window;
        for (k, &i) in self.idx.iter().enumerate() {
            let c = state[i].1;
            let tau = if g[k] > 0.0 { c / g[k] } else { f64::INFINITY };
            let f = (target / tau).clamp(FACTOR_MIN, 1.0);
            self.f0[k] = f;
            let r = &mut self.sum.regions[k];
            r.heat_capacity = c;
            r.conductance = Some(g[k]);
            r.time_constant = Some(tau);
            r.factor = f;
            lines.push(format!(
                "thermal acceleration: region {} G = {:.4e} W/K, C = {c:.4e} J/K, τ = {tau:.4e} \
                 s, factor {f:.4e} on ρ c_p",
                self.names[i], g[k]
            ));
        }
        self.phase = Phase::Accelerated;
        self.sum.accelerated_from = Some(t);
        let f = self.f0.clone();
        self.apply(fl, &f);
        Ok(())
    }

    /// The skin depths of the accelerated regions against their thickness
    /// (with a `period`): the regions where the skin is thicker (ACC003).
    pub fn skin_check(&mut self, fl: &mut Flow) -> Result<Vec<String>> {
        let Some(p) = self.cfg.period else {
            return Ok(Vec::new());
        };
        let Some((sc, dev)) = crate::cht::energy_fields(fl) else {
            return Ok(Vec::new());
        };
        let state = sc.energy_region_state(dev)?;
        let k = sc.energy_region_conductivity();
        let faces = sc.energy_interfaces();
        let omega = 2.0 * std::f64::consts::PI / p;
        let mut thick = Vec::new();
        for (r, (&i, &f)) in self
            .sum
            .regions
            .iter_mut()
            .zip(self.idx.iter().zip(&self.f0))
        {
            let (v, c) = (state[i].0, state[i].1);
            let rho_cp = if v > 0.0 { c / v } else { 0.0 };
            let name = &self.names[i];
            let area: f64 = faces
                .iter()
                .filter(|(n, ..)| interface_touches(n, name, &self.names))
                .map(|x| x.2)
                .sum();
            if rho_cp > 0.0 && k[i] > 0.0 {
                let d = (2.0 * k[i] / (f * rho_cp * omega)).sqrt();
                r.skin_depth = Some(d);
                if area > 0.0 {
                    let l = v / area;
                    r.thickness = Some(l);
                    if d > l {
                        thick.push(name.clone());
                    }
                }
            }
        }
        Ok(thick)
    }

    /// The summary at the end of the run at time `t`.
    pub fn summary(&self, t: f64) -> AccelerationSummary {
        let mut s = self.sum.clone();
        s.history_physical = s
            .physical_from
            .is_some_and(|p| t >= p - 1e-9 * p.abs().max(1.0));
        s
    }
}

/// True when interface `name` (as [`crate::scalar::Scalars::energy_interfaces`]
/// names them) touches region `r`.
fn interface_touches(name: &str, r: &str, names: &[String]) -> bool {
    side_of(name, names).is_some_and(|(a, b)| a == r || b == r)
}

/// The two regions of an interface name: `interface-<solid>` (the fluid
/// and the solid) or `interface-<a>-<b>` (two solids).
fn side_of<'a>(name: &str, names: &'a [String]) -> Option<(&'a str, &'a str)> {
    let rest = name.strip_prefix(crate::cht::INTERFACE)?;
    if let Some(s) = names.iter().skip(1).find(|n| n.as_str() == rest) {
        return Some((names[0].as_str(), s.as_str()));
    }
    for a in names.iter().skip(1) {
        for b in names.iter().skip(1) {
            if a != b && rest == format!("{a}-{b}") {
                return Some((a.as_str(), b.as_str()));
            }
        }
    }
    None
}

/// The heat flow out of each region through its interfaces (W), from the
/// interface fluxes of the last assembly.
fn heat_out(sc: &crate::scalar::Scalars, names: &[String]) -> Vec<f64> {
    let mut out = vec![0.0; names.len()];
    for (n, q, ..) in sc.energy_interfaces() {
        let Some((a, b)) = side_of(&n, names) else {
            continue;
        };
        let ia = names.iter().position(|x| x == a);
        let ib = names.iter().position(|x| x == b);
        if let (Some(ia), Some(ib)) = (ia, ib) {
            if ia == 0 {
                // `interface-<solid>`: q into the fluid, out of the solid
                out[ib] += q;
                out[0] -= q;
            } else {
                // `interface-<a>-<b>`: q into a
                out[ia] -= q;
                out[ib] += q;
            }
        }
    }
    out
}

/// The conductance of each region `idx` to the rest (W/K): the change of
/// its interface heat flow when it is 1 K warmer, the fluid quasi-steady
/// (steady energy solves with every solid held), plus the conductance of
/// its outer boundaries. The state of the energy equation is put back.
fn conductances(fl: &mut Flow, idx: &[usize], t: f64) -> Result<Vec<f64>> {
    let areas: Vec<f64> = match fl.cht.as_deref() {
        Some(c) => {
            let ni = c.view.mesh.neighbour.len();
            let nf = c.view.mesh.owner.len();
            (ni..nf).map(|f| c.view.geo.normal_area(f).1).collect()
        }
        None => return Ok(vec![0.0; idx.len()]),
    };
    let Some((sc, dev)) = crate::cht::energy_fields(fl) else {
        return Ok(vec![0.0; idx.len()]);
    };
    let names = sc.energy_regions();
    let t0 = sc.energy_temperature(dev)?;
    let cell = sc.energy_cell_regions();
    let outer = sc.energy_outer_conductance(&areas, t);
    let held = |shift: Option<usize>| -> Vec<Option<f64>> {
        t0.iter()
            .zip(&cell)
            .map(|(&x, &m)| match m {
                0 => None,
                m if Some(m as usize) == shift => Some(x + 1.0),
                _ => Some(x),
            })
            .collect()
    };
    let base = steady_heat_out(fl, &held(None), &names)?;
    let mut g = Vec::with_capacity(idx.len());
    for &i in idx {
        let q = steady_heat_out(fl, &held(Some(i)), &names)?;
        g.push((q[i] - base[i]).max(0.0) + outer[i]);
    }
    if let Some((sc, dev)) = crate::cht::energy_fields(fl) {
        sc.set_energy_frozen(dev, None)?;
        sc.set_energy_temperature(dev, &t0)?;
    }
    Ok(g)
}

/// The steady heat flow out of each region with the cells of `fixed` held.
fn steady_heat_out(fl: &mut Flow, fixed: &[Option<f64>], names: &[String]) -> Result<Vec<f64>> {
    if let Some((sc, dev)) = crate::cht::energy_fields(fl) {
        sc.set_energy_frozen(dev, Some(fixed))?;
        let t: Vec<f64> = fixed.iter().map(|x| x.unwrap_or(f64::NAN)).collect();
        let cur = sc.energy_temperature(dev)?;
        let start: Vec<f64> = t
            .iter()
            .zip(&cur)
            .map(|(&a, &b)| if a.is_nan() { b } else { a })
            .collect();
        sc.set_energy_temperature(dev, &start)?;
    }
    let mut last: Option<Vec<f64>> = None;
    for _ in 0..PROBE_ITERATIONS {
        let mut res = crate::simplec::Residuals::default();
        crate::cht::solve_with(fl, None, &mut res, false)?;
        let Some((sc, _)) = crate::cht::energy_fields(fl) else {
            break;
        };
        let q = heat_out(sc, names);
        let done = last.as_ref().is_some_and(|l| {
            let scale = q.iter().fold(0.0_f64, |a, b| a.max(b.abs())).max(1e-300);
            q.iter().zip(l).all(|(a, b)| (a - b).abs() <= 1e-6 * scale)
        });
        last = Some(q);
        if done {
            break;
        }
    }
    Ok(last.unwrap_or_else(|| vec![0.0; names.len()]))
}
