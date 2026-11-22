//! Cycle-averaged conjugate coupling of a periodic flow with its solids
//! (D204, stage 7 specification §9d).
//!
//! A periodic flow (an engine exhaust manifold over the engine cycle) and
//! solids that need thousands of cycles to warm up: the fully coupled
//! transient run is far too costly, and below a thin skin the solids carry
//! an almost steady temperature field. The loop alternates two phases on
//! the one mesh and the one energy matrix of the conjugate solve (D158):
//!
//! - **Fluid phase.** The flow and the fluid energy run transient over
//!   `warmup_cycles` + `average_cycles` periods with every solid cell held
//!   at its temperature (the fixed rows of the energy matrix). Per
//!   fluid–solid interface face the run accumulates over the averaging
//!   cycles ∫ G dt, ∫ G T_s* dt and ∫ Q dt, with G = |S|/(R_fluid + R″_c +
//!   R_solid) the conductance of the face (the fluid side, a contact
//!   resistance and the solid half cell in series), T_s* the solid cell's
//!   temperature at the face and Q the heat flow into the solid (the exact
//!   interface flux of stage 7d2). Ḡ = ⟨G⟩ and T̄_ref = (⟨G T_s*⟩ + ⟨Q⟩)/⟨G⟩
//!   give Ḡ (T̄_ref − T_s*) = ⟨Q⟩: the cycle-mean heat flow of every face
//!   exactly, with its first-order response to a change of the solid
//!   temperature. Each phase also checks that the flow is periodic (the
//!   mean interface heat flow of the last two cycles) and adds cycles up to
//!   `max_cycles` when it is not.
//! - **Solid phase** (`solid_model`), steady, from the state at the end
//!   of the fluid phase:
//!   - `"mean-flow"` (default): the steady energy equation of fluid and
//!     solids with the cycle-mean mass flux of the averaging cycles,
//!     corrected per face to the measured cycle mean. A first solve with
//!     the solids held gives the model's heat flow Q_m per face at the
//!     held wall; the second, with everything free, adds the defect
//!     Q̄ − Q_m to the solid cell of every face. At the held wall the
//!     solids thus receive the measured cycle-mean heat flow exactly, and
//!     their response to a change of the wall is that of the mean flow:
//!     the warming of the fluid along the wall included, which a local
//!     Robin pair cannot carry (defect correction of a cheaper model, as
//!     in multigrid and multi-fidelity iterations).
//!   - `"robin"`: the solids alone, the fluid held, the Robin pair
//!     (Ḡ, T̄_ref) on the solid side of every fluid–solid face (the face
//!     itself carries no conduction in the matrix; −Ḡ on the diagonal of
//!     the solid cell, Ḡ T̄_ref on its right-hand side). This is the
//!     classic practice. Its flaw: Ḡ of the first fluid cell is far
//!     stiffer than the response of the whole fluid, whose near-wall
//!     temperature follows the wall and whose bulk warms with the heat it
//!     takes up. The updates then contract slowly: about 4 % per update
//!     in the gallery channel, where the mean-flow model needs 3 or 4
//!     updates in all.
//!
//!   Interfaces between two solids, contact resistances, sources and the
//!   outer boundaries of the solids stay as they are.
//! - **Update** (`update`). The new solid temperatures are held in the
//!   next fluid phase, which continues from the last flow state and its
//!   time levels. The first update takes the solved field with
//!   `relaxation`. Later ones take the interface quasi-Newton update
//!   (`"quasi-newton"`, the default, IQN-ILS), the Aitken factor
//!   (`"aitken"`) or the fixed factor again (`"relaxation"`). The loop stops when the largest change of a wall
//!   temperature is below `tolerance` and the interface heat flow changes
//!   by less than `heat_tolerance` between two fluid phases, or after
//!   `max_updates` (CYC004).
//!
//! The steady solid phase is the limit of a vanishing ρ c_p of the solids
//! (the thermal acceleration of D213 with f → 0), with the cycle averaging
//! taking out the swing that a small f would blow up. The swing of the wall
//! temperature over a cycle (`skin = "effusivity"`), fully coupled final
//! cycles and a restart in any phase follow in the next step of D204.
//!
//! Ledger: in a fluid phase the heat into the held solids is booked as
//! `held` (it is not stored), so the balance of every step closes. In a
//! solid phase the balance of the solids closes to the solver tolerance,
//! with the Robin heat or the defect of the mean-flow model booked as heat
//! from outside the model (`cycle_coupling.csv`, `result.json`).

use crate::conduction::RobinFace;
use crate::diag::Diagnostic;
use crate::energy::CycleCoupling;
use crate::simplec::Flow;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;

/// Steady solves of a solid phase at most (it stops when the largest
/// change of a solid temperature is below 10⁻⁹ of the temperature span).
const SOLID_ITERATIONS: usize = 200;

/// The bounds of the Aitken factor of the updates.
const OMEGA_MIN: f64 = 0.05;
const OMEGA_MAX: f64 = 50.0;

/// One update of the loop in the result and in `cycle_coupling.csv`.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct CycleUpdate {
    /// 1, 2, … .
    pub update: u32,
    /// The time at the end of the fluid phase (s).
    pub time: f64,
    /// The cycles of the fluid phase.
    pub cycles: u32,
    /// The fluid phase was periodic within `cycle_tolerance`.
    pub periodic: bool,
    /// The relative change of the mean interface heat flow between the
    /// last two cycles of the fluid phase.
    pub cycle_change: f64,
    /// Per interface (`interface-<solid>`): the cycle-mean heat flow into
    /// the fluid over the averaging cycles (W).
    pub heat_flow: Vec<(String, f64)>,
    /// The relative change of the total interface heat flow against the
    /// fluid phase before (none at the first update).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_change: Option<f64>,
    /// The largest and the area-mean change of the wall temperature
    /// (solid side of the interface faces) of the solid phase (K).
    pub wall_change_max: f64,
    /// See `wall_change_max`.
    pub wall_change_mean: f64,
    /// The area-mean wall temperature after the solid phase (K).
    pub wall_mean: f64,
    /// The steady solves of the solid phase and the balance of the solids
    /// (W, and relative to the largest heat flow).
    pub solid_iterations: usize,
    /// The relaxation factor of this update (`relaxation` at the first
    /// update and with `update = "relaxation"`, the Aitken factor with
    /// `"aitken"`, 1 for the quasi-Newton update).
    pub relaxation: f64,
    /// See `solid_iterations`.
    pub solid_imbalance: f64,
    /// See `solid_iterations`.
    pub solid_relative: f64,
}

/// The loop in the result (`transient.cycle_coupling`).
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct CycleSummary {
    /// The period (s) and the time steps per cycle.
    pub period: f64,
    /// See `period`.
    pub steps_per_cycle: usize,
    /// The fluid–solid interface faces.
    pub faces: usize,
    /// The updates.
    pub updates: Vec<CycleUpdate>,
    /// The loop converged within `tolerance` and `heat_tolerance`.
    pub converged: bool,
    /// The cycles of all fluid phases.
    pub cycles: u32,
    /// The time when the loop ended (s), none when the run ended first.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ended_at: Option<f64>,
}

/// Per interface face: the integrals over a cycle.
#[derive(Debug, Clone, Default)]
struct Acc {
    time: f64,
    /// ∫ G, ∫ G T_s*, ∫ Q, ∫ G_fluid, ∫ T_w (per face, times dt).
    g: Vec<f64>,
    gts: Vec<f64>,
    q: Vec<f64>,
    gf: Vec<f64>,
    tw: Vec<f64>,
    /// ∫ ṁ dt per face of the fluid view (the mean-flow model).
    flux: Vec<f64>,
}

impl Acc {
    fn new(n: usize) -> Self {
        Self {
            flux: Vec::new(),
            time: 0.0,
            g: vec![0.0; n],
            gts: vec![0.0; n],
            q: vec![0.0; n],
            gf: vec![0.0; n],
            tw: vec![0.0; n],
        }
    }

    fn add(&mut self, o: &Acc) {
        self.time += o.time;
        for (a, b) in [
            (&mut self.g, &o.g),
            (&mut self.gts, &o.gts),
            (&mut self.q, &o.q),
            (&mut self.gf, &o.gf),
            (&mut self.tw, &o.tw),
        ] {
            for (x, y) in a.iter_mut().zip(b) {
                *x += y;
            }
        }
        if self.flux.is_empty() {
            self.flux = o.flux.clone();
        } else {
            for (x, y) in self.flux.iter_mut().zip(&o.flux) {
                *x += y;
            }
        }
    }
}

/// One fluid–solid interface face.
#[derive(Debug, Clone)]
struct Face {
    slot: usize,
    face: u32,
    solid_owner: bool,
    pair: usize,
    area: f64,
    solid_cell: u32,
}

/// The averages of a fluid phase per face, for `cycle_faces.csv`.
#[derive(Debug, Clone, Default)]
struct FaceMeans {
    /// Ḡ, Ḡ_fluid (W/K), T̄_ref (K), Q̄ into the solid (W), T_w after the
    /// solid phase (K).
    g: Vec<f64>,
    gf: Vec<f64>,
    t_ref: Vec<f64>,
    q: Vec<f64>,
    tw: Vec<f64>,
}

/// The loop during a run.
pub struct Controller {
    cfg: CycleCoupling,
    dim: u8,
    faces: Vec<Face>,
    pairs: Vec<String>,
    solid: Vec<bool>,
    cycle_start: f64,
    cycles_in_phase: u32,
    cycles_total: u32,
    cur: Acc,
    hist: VecDeque<Acc>,
    cycle_q: Vec<f64>,
    tw_prev: Option<Vec<f64>>,
    q_prev: Option<f64>,
    cycle_abs: f64,
    omega: f64,
    /// The residuals x̃ − x and the solved fields x̃ of the solid cells of
    /// the earlier updates (the quasi-Newton and Aitken updates).
    hist_r: Vec<Vec<f64>>,
    hist_xt: Vec<Vec<f64>>,
    updates: Vec<CycleUpdate>,
    means: FaceMeans,
    done: bool,
    converged: bool,
    ended_at: Option<f64>,
    /// Diagnostics of the loop (CYC002, CYC004) for the run's warnings.
    pub diagnostics: Vec<Diagnostic>,
    steps_per_cycle: usize,
}

impl Controller {
    /// The loop of a conjugate case from time `t0`: the solids are held at
    /// their current temperatures from the first step on.
    pub fn new(fl: &mut Flow, cfg: &CycleCoupling, t0: f64, dt: f64, dim: u8) -> Result<Self> {
        let (faces, pairs, solid) = {
            let (sc, dev) = crate::cht::energy_fields(fl).ok_or_else(|| {
                anyhow::anyhow!("[energy.cycle] needs the energy equation of a conjugate case")
            })?;
            let faces: Vec<Face> = sc
                .energy_fluid_solid_faces()
                .into_iter()
                .map(|(slot, face, solid_owner, st)| {
                    let (o, nb) = sc.energy_interface_cells(slot);
                    let solid_cell = if solid_owner { o } else { nb };
                    Face {
                        slot,
                        face,
                        solid_owner,
                        pair: sc.energy_interface_pair(slot),
                        area: st[0],
                        solid_cell,
                    }
                })
                .collect();
            let pairs: Vec<String> = sc.energy_interfaces().into_iter().map(|x| x.0).collect();
            let solid: Vec<bool> = sc.energy_cell_regions().iter().map(|&m| m != 0).collect();
            let t = sc.energy_temperature(dev)?;
            hold(sc, dev, &solid, &t)?;
            (faces, pairs, solid)
        };
        anyhow::ensure!(
            !faces.is_empty(),
            "[energy.cycle]: the case has no interface faces between the fluid and a solid"
        );
        let n = faces.len();
        Ok(Self {
            steps_per_cycle: (cfg.period / dt).round() as usize,
            cfg: cfg.clone(),
            dim,
            faces,
            pairs,
            solid,
            cycle_start: t0,
            cycles_in_phase: 0,
            cycles_total: 0,
            cur: Acc::new(n),
            hist: VecDeque::new(),
            cycle_q: Vec::new(),
            tw_prev: None,
            q_prev: None,
            cycle_abs: 0.0,
            omega: cfg.relaxation,
            hist_r: Vec::new(),
            hist_xt: Vec::new(),

            updates: Vec::new(),
            means: FaceMeans::default(),
            done: false,
            converged: false,
            ended_at: None,
            diagnostics: Vec::new(),
        })
    }

    /// The scale of a change of the interface heat flow: the net flow, or
    /// the sum of the magnitudes over the faces of the last cycle where the
    /// net flow is small (heat that goes in at one end and out at the
    /// other).
    fn heat_scale(&self, net: f64) -> f64 {
        net.max(self.cycle_abs).max(1e-300)
    }

    /// True when the loop has ended (converged or at `max_updates`): the
    /// run stops.
    pub fn done(&self) -> bool {
        self.done
    }

    /// After a time step that ended at `t` and took `dt`: accumulates the
    /// interface faces, closes a cycle, runs a solid phase and an update
    /// at the end of a fluid phase. Returns the lines for the log.
    pub fn after_step(&mut self, fl: &mut Flow, t: f64, dt: f64) -> Result<Vec<String>> {
        let mut lines = Vec::new();
        if self.done {
            return Ok(lines);
        }
        {
            let (sc, _) = fields(fl)?;
            let st = sc.energy_fluid_solid_faces();
            let c = &mut self.cur;
            c.time += dt;
            for (i, (_, _, _, x)) in st.iter().enumerate() {
                let [_, g, gf, q, ts, tw] = *x;
                c.g[i] += g * dt;
                c.gts[i] += g * ts * dt;
                c.q[i] += q * dt;
                c.gf[i] += gf * dt;
                c.tw[i] += tw * dt;
            }
        }
        if self.cfg.solid_model != "robin" {
            let b = fl.convecting_flux()?;
            let m = fl.dev.download(b)?;
            let c = &mut self.cur;
            if c.flux.is_empty() {
                c.flux = vec![0.0; m.len()];
            }
            for (x, y) in c.flux.iter_mut().zip(&m) {
                *x += y * dt;
            }
        }
        if t - self.cycle_start < self.cfg.period - 1e-6 * dt {
            return Ok(lines);
        }
        // a cycle is complete
        self.cycle_start = t;
        self.cycles_in_phase += 1;
        self.cycles_total += 1;
        let done = std::mem::replace(&mut self.cur, Acc::new(self.faces.len()));
        self.cycle_q.push(done.q.iter().sum::<f64>() / done.time);
        self.cycle_abs = done.q.iter().map(|x| x.abs()).sum::<f64>() / done.time;
        self.hist.push_back(done);
        while self.hist.len() > self.cfg.average_cycles as usize {
            self.hist.pop_front();
        }
        let need = self.cfg.warmup_cycles + self.cfg.average_cycles;
        if self.cycles_in_phase < need.max(2) {
            return Ok(lines);
        }
        let k = self.cycle_q.len();
        let (a, b) = (self.cycle_q[k - 1], self.cycle_q[k - 2]);
        let change = (a - b).abs() / self.heat_scale(a.abs().max(b.abs()));
        let periodic = change <= self.cfg.cycle_tolerance;
        if !periodic && self.cycles_in_phase < self.cfg.max_cycles {
            lines.push(format!(
                "[energy.cycle] t = {t} s: cycle {} of the fluid phase, the mean interface heat \
                 flow still changes by {:.2e} (> {:.2e}): one more cycle",
                self.cycles_in_phase, change, self.cfg.cycle_tolerance
            ));
            return Ok(lines);
        }
        if !periodic {
            let m = format!(
                "[energy.cycle] update {}: the fluid phase ended at max_cycles = {} with the mean \
                 interface heat flow still changing by {:.2e} between its last two cycles (> \
                 cycle_tolerance = {:.2e}); the averages are not those of a periodic flow",
                self.updates.len() + 1,
                self.cfg.max_cycles,
                change,
                self.cfg.cycle_tolerance
            );
            self.diagnostics
                .push(Diagnostic::new(crate::diag::CYC_NOT_PERIODIC, m.clone()));
            lines.push(m);
        }
        lines.extend(self.update(fl, t, periodic, change)?);
        Ok(lines)
    }

    /// The solid phase and the update at the end of a fluid phase.
    fn update(
        &mut self,
        fl: &mut Flow,
        t: f64,
        periodic: bool,
        change: f64,
    ) -> Result<Vec<String>> {
        let mut lines = Vec::new();
        let n = self.faces.len();
        let mut sum = Acc::new(n);
        for h in &self.hist {
            sum.add(h);
        }
        let tm = sum.time;
        let g: Vec<f64> = sum.g.iter().map(|x| x / tm).collect();
        let gf: Vec<f64> = sum.gf.iter().map(|x| x / tm).collect();
        let q: Vec<f64> = sum.q.iter().map(|x| x / tm).collect();
        // T_s* of the phase (the held solid cell at the face) and the Robin
        // pair Ḡ (T̄_ref − T_s*) = Q̄ of every face
        let ts: Vec<f64> = (0..n).map(|i| sum.gts[i] / sum.g[i]).collect();
        let t_ref: Vec<f64> = (0..n).map(|i| ts[i] + q[i] / g[i]).collect();
        let tw_fluid: Vec<f64> = sum.tw.iter().map(|x| x / tm).collect();
        // the cycle-mean heat flow per interface, into the fluid
        let mut per = vec![0.0; self.pairs.len()];
        for (f, qi) in self.faces.iter().zip(&q) {
            per[f.pair] -= qi;
        }
        let q_total: f64 = per.iter().sum();
        let heat_change = self
            .q_prev
            .map(|p| (q_total - p).abs() / self.heat_scale(q_total.abs().max(p.abs())));
        // the solid phase, from the state at the end of the fluid phase
        // (its time levels are put back afterwards)
        let saved = {
            let (sc, dev) = fields(fl)?;
            sc.energy_levels(dev)?
        };
        let t_before = saved.0.clone();
        let span = {
            let (lo, hi) = t_ref
                .iter()
                .chain(t_before.iter())
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| {
                    (a.min(x), b.max(x))
                });
            (hi - lo).max(1.0)
        };
        let (last, its) = if self.cfg.solid_model == "robin" {
            // the solids alone with the Robin pairs, the fluid held
            let (sc, dev) = fields(fl)?;
            let fluid: Vec<bool> = self.solid.iter().map(|s| !s).collect();
            hold(sc, dev, &fluid, &t_before)?;
            let mut robin = vec![None; sc.energy_interface_slots()];
            for (i, f) in self.faces.iter().enumerate() {
                robin[f.slot] = Some(RobinFace {
                    g: g[i],
                    t_ref: t_ref[i],
                    solid_owner: f.solid_owner,
                });
            }
            sc.set_energy_robin(dev, Some(robin))?;
            let r = steady(fl, None, span)?;
            let (sc, dev) = fields(fl)?;
            sc.set_energy_robin(dev, None)?;
            r
        } else {
            // the mean-flow model with its defect: the steady energy
            // equation of fluid and solids with the cycle-mean mass flux;
            // first with the solids held, which gives the model's heat flow
            // Q_m per face at the held wall, then with everything free and
            // the defect Q̄ − Q_m added to the solid cell of every face, so
            // that at the held wall the solids receive the measured
            // cycle-mean heat flow exactly, and their response to a change
            // of the wall (the warming of the fluid along the wall included)
            // is that of the mean flow. The defect goes into the solids
            // only (the model's fluid stays that of the held wall); the
            // ledger books it with the Robin heat as heat from outside the
            // model.
            let mf: Vec<f64> = sum.flux.iter().map(|x| x / tm).collect();
            {
                let (sc, dev) = fields(fl)?;
                hold(sc, dev, &self.solid, &t_before)?;
            }
            let (_, i1) = steady(fl, Some(&mf), span)?;
            let (sc, dev) = fields(fl)?;
            let qm: Vec<f64> = sc
                .energy_fluid_solid_faces()
                .iter()
                .map(|x| x.3[3])
                .collect();
            let mut extra = vec![0.0; t_before.len()];
            for (i, f) in self.faces.iter().enumerate() {
                let d = q[i] - qm[i];
                extra[f.solid_cell as usize] += d;
            }
            sc.set_energy_frozen(dev, None)?;
            sc.set_energy_extra_heat(Some(extra));
            let (last, i2) = steady(fl, Some(&mf), span)?;
            (last, i1 + i2)
        };
        // the balance of the solids in the solid phase: the fluid counted
        // as held (what flows into it, convection included, is booked as
        // `held`, so the balance is that of the solids with the Robin heat
        // or the defect), the old levels at the solved field (a steady
        // solve stores nothing)
        {
            let (sc, dev) = fields(fl)?;
            sc.set_energy_levels(dev, &last)?;
            let fluid: Vec<bool> = self.solid.iter().map(|s| !s).collect();
            hold(sc, dev, &fluid, &last)?;
        }
        let sums = crate::cht::summaries(fl, self.dim, t)?;
        let bal = sums
            .iter()
            .find(|s| s.name == crate::energy::FIELD)
            .and_then(|s| s.energy.clone())
            .unwrap_or_default();
        let (sc, dev) = fields(fl)?;
        sc.set_energy_extra_heat(None);
        let st = sc.energy_fluid_solid_faces();
        let tw: Vec<f64> = st.iter().map(|x| x.3[5]).collect();
        // the update: the solids from the solved field, held with their old
        // time levels in the next fluid phase. The plain update x ← x +
        // ω (x̃ − x) contracts slowly: the Robin pair of the first fluid cell
        // is stiffer than the true response of the fluid (its near-wall
        // temperature follows the wall). "quasi-newton" (default): the
        // interface quasi-Newton method with a least-squares model of the
        // inverse Jacobian of the residual r = x̃ − x from all earlier
        // updates (IQN-ILS, Degroote, Bathe & Vierendeels 2009): with V the
        // differences of the residuals and W those of the solved fields,
        // c = argmin |V c + r_k| and x_{k+1} = x̃_k + W c; "aitken": the
        // dynamic factor ω_k = −ω_{k−1} r_{k−1}·(r_k − r_{k−1})/|r_k −
        // r_{k−1}|² (Irons & Tuck 1969; Küttler & Wall 2008); "relaxation":
        // the fixed factor. The first update is the plain one with
        // `relaxation`.
        let idx: Vec<usize> = (0..self.solid.len()).filter(|&c| self.solid[c]).collect();
        let x: Vec<f64> = idx.iter().map(|&c| t_before[c]).collect();
        let xt: Vec<f64> = idx.iter().map(|&c| last[c]).collect();
        let r: Vec<f64> = xt.iter().zip(&x).map(|(a, b)| a - b).collect();
        let (x_new, w) = match (self.cfg.update.as_str(), self.hist_r.is_empty()) {
            (_, true) | ("relaxation", _) => {
                let w = self.cfg.relaxation;
                (x.iter().zip(&r).map(|(a, d)| a + w * d).collect(), w)
            }
            ("aitken", false) => {
                let rp = &self.hist_r[self.hist_r.len() - 1];
                let (mut num, mut den) = (0.0, 0.0);
                for (a, b) in rp.iter().zip(&r) {
                    let d = b - a;
                    num += a * d;
                    den += d * d;
                }
                let w = if den > 0.0 {
                    (-self.omega * num / den).clamp(OMEGA_MIN, OMEGA_MAX)
                } else {
                    self.omega
                };
                (x.iter().zip(&r).map(|(a, d)| a + w * d).collect(), w)
            }
            _ => (quasi_newton(&self.hist_r, &self.hist_xt, &r, &xt), 1.0),
        };
        self.omega = w;
        self.hist_r.push(r);
        self.hist_xt.push(xt);
        // the fluid back at its state and time levels of the end of the
        // fluid phase, the solids at the new temperatures on every level
        let (mut t_new, mut o1, mut o2) = saved;
        for (&c, v) in idx.iter().zip(x_new) {
            t_new[c] = v;
            o1[c] = v;
            o2[c] = v;
        }
        sc.set_energy_levels3(dev, &t_new, &o1, &o2)?;
        hold(sc, dev, &self.solid, &t_new)?;
        let reference = self.tw_prev.as_ref().unwrap_or(&tw_fluid);
        let area: f64 = self.faces.iter().map(|f| f.area).sum();
        let (mut dmax, mut dmean, mut wmean) = (0.0_f64, 0.0, 0.0);
        for ((f, a), b) in self.faces.iter().zip(&tw).zip(reference) {
            dmax = dmax.max((a - b).abs());
            dmean += f.area * (a - b).abs();
            wmean += f.area * a;
        }
        dmean /= area;
        wmean /= area;
        let k = self.updates.len() as u32 + 1;
        let row = CycleUpdate {
            update: k,
            time: t,
            cycles: self.cycles_in_phase,
            periodic,
            cycle_change: change,
            heat_flow: self
                .pairs
                .iter()
                .cloned()
                .zip(per.iter().copied())
                .collect(),
            heat_change,
            wall_change_max: dmax,
            wall_change_mean: dmean,
            wall_mean: wmean,
            solid_iterations: its,
            relaxation: w,
            solid_imbalance: bal.imbalance,
            solid_relative: bal.relative,
        };
        lines.push(format!(
            "[energy.cycle] update {k} at t = {t} s after {} cycles: interface heat flow {:.6e} W \
             into the fluid{}, wall temperature change max {:.4} K, mean {:.4} K, mean wall {:.4} \
             K; solid phase {its} solves, balance {:.2e}, relaxation {:.3}",
            self.cycles_in_phase,
            q_total,
            heat_change.map_or(String::new(), |c| format!(" (change {c:.2e})")),
            dmax,
            dmean,
            wmean,
            bal.relative,
            w
        ));
        self.updates.push(row);
        self.means = FaceMeans {
            g,
            gf,
            t_ref,
            q,
            tw: tw.clone(),
        };
        let converged =
            heat_change.is_some_and(|c| c <= self.cfg.heat_tolerance) && dmax <= self.cfg.tolerance;
        self.tw_prev = Some(tw);
        self.q_prev = Some(q_total);
        if converged {
            self.done = true;
            self.converged = true;
            self.ended_at = Some(t);
            lines.push(format!(
                "[energy.cycle] converged after {k} updates and {} cycles: the run stops",
                self.cycles_total
            ));
        } else if k >= self.cfg.max_updates {
            self.done = true;
            self.ended_at = Some(t);
            let m = format!(
                "[energy.cycle] stopped at max_updates = {k} before the tolerances: the last \
                 wall temperature change {dmax:.4} K (tolerance {} K), interface heat flow change \
                 {} (heat_tolerance {})",
                self.cfg.tolerance,
                heat_change.map_or("-".to_string(), |c| format!("{c:.2e}")),
                self.cfg.heat_tolerance
            );
            self.diagnostics
                .push(Diagnostic::new(crate::diag::CYC_MAX_UPDATES, m.clone()));
            lines.push(m);
        }
        // the next fluid phase
        self.cycles_in_phase = 0;
        self.hist.clear();
        self.cycle_q.clear();
        self.cur = Acc::new(n);
        Ok(lines)
    }

    /// The loop for the result at the end of the run; a run that ended
    /// before the loop gets CYC004.
    pub fn summary(&mut self) -> CycleSummary {
        if !self.done {
            let m = format!(
                "[energy.cycle]: the run ended (end_time or steps) after {} updates, before the \
                 loop converged or reached max_updates",
                self.updates.len()
            );
            self.diagnostics
                .push(Diagnostic::new(crate::diag::CYC_MAX_UPDATES, m));
            self.done = true;
        }
        CycleSummary {
            period: self.cfg.period,
            steps_per_cycle: self.steps_per_cycle,
            faces: self.faces.len(),
            updates: self.updates.clone(),
            converged: self.converged,
            cycles: self.cycles_total,
            ended_at: self.ended_at,
        }
    }

    /// The header and rows of `cycle_coupling.csv`.
    pub fn csv(&self) -> (String, Vec<String>) {
        let mut head = String::from("update,time,cycles,periodic,cycle_change");
        for p in &self.pairs {
            head.push_str(&format!(",{p}_heat"));
        }
        head.push_str(
            ",heat_change,wall_change_max,wall_change_mean,wall_mean,solid_iterations,\
             relaxation,solid_imbalance,solid_relative",
        );
        let rows = self
            .updates
            .iter()
            .map(|u| {
                let mut r = format!(
                    "{},{:e},{},{},{:e}",
                    u.update,
                    u.time,
                    u.cycles,
                    u8::from(u.periodic),
                    u.cycle_change
                );
                for (_, q) in &u.heat_flow {
                    r.push_str(&format!(",{q:e}"));
                }
                r.push_str(&format!(
                    ",{},{:e},{:e},{:e},{},{:e},{:e},{:e}",
                    u.heat_change.map_or(String::new(), |c| format!("{c:e}")),
                    u.wall_change_max,
                    u.wall_change_mean,
                    u.wall_mean,
                    u.solid_iterations,
                    u.relaxation,
                    u.solid_imbalance,
                    u.solid_relative
                ));
                r
            })
            .collect();
        (head, rows)
    }

    /// The header and rows of `cycle_faces.csv`: per fluid–solid interface
    /// face of the last update, with the face centroids `centroid`.
    pub fn faces_csv(&self, centroid: &[[f64; 3]]) -> (String, Vec<String>) {
        let head = "face,interface,x,y,z,area,h_fluid,conductance,t_ref,heat_flux,t_wall".into();
        let m = &self.means;
        let rows = if m.g.is_empty() {
            Vec::new()
        } else {
            self.faces
                .iter()
                .enumerate()
                .map(|(i, f)| {
                    let c = centroid[f.face as usize];
                    format!(
                        "{},{},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e}",
                        f.face,
                        self.pairs[f.pair],
                        c[0],
                        c[1],
                        c[2],
                        f.area,
                        m.gf[i] / f.area,
                        m.g[i] / f.area,
                        m.t_ref[i],
                        m.q[i] / f.area,
                        m.tw[i]
                    )
                })
                .collect()
        };
        (head, rows)
    }
}

/// Holds the cells with `which` at the temperatures `t`, releases the rest.
fn hold(
    sc: &mut crate::scalar::Scalars,
    dev: &mut dyn crate::device::FlowDevice,
    which: &[bool],
    t: &[f64],
) -> Result<()> {
    let fixed: Vec<Option<f64>> = which.iter().zip(t).map(|(&w, &x)| w.then_some(x)).collect();
    sc.set_energy_frozen(dev, Some(&fixed))
}

/// The quasi-Newton update of IQN-ILS (Degroote, Bathe & Vierendeels 2009)
/// from the residuals `hr` and solved fields `hx` of the earlier updates and
/// the current residual `r` and solved field `xt`: x̃ + W c with c the
/// least-squares solution of V c = −r, V and W the differences of the
/// current against each earlier residual and solved field. The normal
/// equations are solved with a small ridge (the columns can be nearly
/// dependent once the loop has converged); columns that add nothing are
/// left out.
fn quasi_newton(hr: &[Vec<f64>], hx: &[Vec<f64>], r: &[f64], xt: &[f64]) -> Vec<f64> {
    let m = hr.len();
    let v: Vec<Vec<f64>> = hr
        .iter()
        .map(|p| r.iter().zip(p).map(|(a, b)| a - b).collect())
        .collect();
    let w: Vec<Vec<f64>> = hx
        .iter()
        .map(|p| xt.iter().zip(p).map(|(a, b)| a - b).collect())
        .collect();
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let mut a = vec![vec![0.0; m]; m];
    let mut b = vec![0.0; m];
    for i in 0..m {
        for j in 0..m {
            a[i][j] = dot(&v[i], &v[j]);
        }
        b[i] = -dot(&v[i], r);
    }
    let scale = (0..m).map(|i| a[i][i]).fold(0.0_f64, f64::max);
    if scale <= 0.0 {
        return xt.to_vec();
    }
    for (i, row) in a.iter_mut().enumerate() {
        row[i] += 1e-10 * scale;
    }
    let c = solve_small(a, b);
    let mut x = xt.to_vec();
    for (ci, wi) in c.iter().zip(&w) {
        for (xv, d) in x.iter_mut().zip(wi) {
            *xv += ci * d;
        }
    }
    x
}

/// Gaussian elimination with partial pivoting of a small dense system.
fn solve_small(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .unwrap_or(k);
        a.swap(k, p);
        b.swap(k, p);
        let d = a[k][k];
        if d == 0.0 {
            continue;
        }
        for i in k + 1..n {
            let f = a[i][k] / d;
            for j in k..n {
                a[i][j] -= f * a[k][j];
            }
            b[i] -= f * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for k in (0..n).rev() {
        let s: f64 = (k + 1..n).map(|j| a[k][j] * x[j]).sum();
        x[k] = if a[k][k] == 0.0 {
            0.0
        } else {
            (b[k] - s) / a[k][k]
        };
    }
    x
}

/// Steady energy solves (with the full-mesh face mass flux `flux` of the
/// fluid view in place of the current one, if given) until the largest
/// change of a temperature is below 10⁻⁹ of `span`, at most
/// `SOLID_ITERATIONS`; the solved temperatures and the number of solves.
fn steady(fl: &mut Flow, flux: Option<&[f64]>, span: f64) -> Result<(Vec<f64>, usize)> {
    let mut last = {
        let (sc, dev) = fields(fl)?;
        sc.energy_temperature(dev)?
    };
    let mut its = 0;
    for _ in 0..SOLID_ITERATIONS {
        let mut res = crate::simplec::Residuals::default();
        crate::cht::solve_with_flux(fl, None, &mut res, false, flux)?;
        its += 1;
        let (sc, dev) = fields(fl)?;
        let now = sc.energy_temperature(dev)?;
        let d = now
            .iter()
            .zip(&last)
            .fold(0.0_f64, |a, (x, y)| a.max((x - y).abs()));
        last = now;
        if d <= 1e-9 * span {
            break;
        }
    }
    Ok((last, its))
}

/// The energy equation of the conjugate case and its device.
fn fields(
    fl: &mut Flow,
) -> Result<(
    &mut crate::scalar::Scalars,
    &mut dyn crate::device::FlowDevice,
)> {
    crate::cht::energy_fields(fl).ok_or_else(|| {
        anyhow::anyhow!("[energy.cycle] needs the energy equation of a conjugate case")
    })
}
