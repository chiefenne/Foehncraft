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
//! taking out the swing that a small f would blow up.
//!
//! - **Skin swing** (`skin = "effusivity"`). Below a skin of depth
//!   δ = √(2a/ω) the solid does not see the swing of the gas; at the
//!   surface it does. The harmonics q̂_n of the heat flux into the solid
//!   of every face over the averaging cycles give the periodic surface
//!   temperature of a semi-infinite solid, T̂_w,n = q̂_n/(e √(iω_n)) with
//!   e = √(k ρ c_p) (Carslaw & Jaeger 1959, §2.6). In the next fluid phase
//!   the held solid cell of the face swings by (Z_n − r_s) q̂_n, Z_n =
//!   1/(e √(iω_n)) and r_s the resistance of the solid half cell, so that
//!   the face temperature T_s + r_s q swings like the surface of the
//!   semi-infinite solid (with the flux of the phase before); the fluid sees
//!   the swing, and the correlation of its heat transfer with the wall
//!   temperature enters the cycle means. The mean of every held cell is
//!   that of the solid phase. CYC003 where δ is not small against the
//!   thickness V/A of a region.
//! - **Resolved cycles** (`resolved_cycles`). After the loop has ended, the
//!   solids are released and the run goes on fully coupled for that many
//!   cycles from the loop's solid temperatures: the exact swing (also at
//!   thin walls and edges) and the heat flow of a fully coupled cycle.
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
use foehncraft_core::exact::ExactSumIter as _;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;

/// Steady solves of a solid phase at most (it stops when the largest
/// change of a solid temperature is below 10⁻⁹ of the temperature span).
const SOLID_ITERATIONS: usize = 200;

/// Solves per step of a transient solid phase at most (the same stop).
const MARCH_ITERATIONS: usize = 40;

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
    /// With `skin = "effusivity"`: the largest amplitude of the surface
    /// swing of a face (half its range over a cycle, K) that the next fluid
    /// phase sees.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub swing_max: Option<f64>,
    /// With `skin = "effusivity"`: the largest change of a harmonic of the
    /// held-cell swing (modulus of the complex difference, K) against the
    /// swing of the fluid phase; it counts against `tolerance` too.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub swing_change: Option<f64>,
    /// With `solid = "transient"`: the physical time the march covered (s,
    /// the fluid phase included, in which the solids were held) ...
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub march_time: Option<f64>,
    /// ... and the jump of the run's clock after the fluid phase (s, whole
    /// periods): the solids are at `time + march_jump` after the update.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub march_jump: Option<f64>,
}

/// The skin swing in the result.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SkinSummary {
    /// The harmonics of the period taken.
    pub harmonics: usize,
    /// Per solid region at the interfaces: its name, the skin depth of the
    /// period √(2k/(ρ c_p ω)) (m) and its thickness V/A (m).
    pub regions: Vec<(String, f64, f64)>,
    /// The interface faces in regions whose skin is not small against
    /// their thickness (CYC003).
    pub thick_faces: usize,
    /// The largest amplitude of the surface swing of a face (K) after the
    /// last update.
    pub swing_max: f64,
}

/// The fully coupled cycles after the loop in the result.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct ResolvedSummary {
    /// The fully coupled cycles run.
    pub cycles: u32,
    /// Per interface: the mean heat flow into the fluid over the last of
    /// them (W).
    pub heat_flow: Vec<(String, f64)>,
    /// The largest amplitude of the wall temperature of a face over the last
    /// cycle (half its range, K).
    pub swing_max: f64,
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
    /// The solid phase: `"steady"` or `"transient"` (a history).
    #[serde(default)]
    pub solid: String,
    /// The skin swing (`skin = "effusivity"`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skin: Option<SkinSummary>,
    /// The fully coupled cycles after the loop (`resolved_cycles`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resolved: Option<ResolvedSummary>,
}

/// Per interface face: the integrals over a cycle.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
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
    /// With the skin swing: ∫ Q² dt and ∫ (T_w − T_s*) Q dt per face (the
    /// resistance of the solid half cell, by least squares), and the
    /// harmonics ∫ q e^(−iω_n t) dt of the flux into the solid (W/m², face
    /// i, harmonic n at i·H + n − 1).
    sq: Vec<f64>,
    stq: Vec<f64>,
    fr: Vec<f64>,
    fi: Vec<f64>,
}

impl Acc {
    fn new(n: usize, h: usize) -> Self {
        Self {
            flux: Vec::new(),
            sq: vec![0.0; if h > 0 { n } else { 0 }],
            stq: vec![0.0; if h > 0 { n } else { 0 }],
            fr: vec![0.0; n * h],
            fi: vec![0.0; n * h],
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
            (&mut self.sq, &o.sq),
            (&mut self.stq, &o.stq),
            (&mut self.fr, &o.fr),
            (&mut self.fi, &o.fi),
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

/// The state of the loop in a restart file (`cht_cycle_state`, JSON): every
/// counter, accumulator and update history, so that a run continued from it
/// is the run without the interruption.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct LoopState {
    faces: usize,
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
    hist_r: Vec<Vec<f64>>,
    hist_xt: Vec<Vec<f64>>,
    updates: Vec<CycleUpdate>,
    means: FaceMeans,
    done: bool,
    converged: bool,
    ended_at: Option<f64>,
    diagnostics: Vec<(String, String)>,
    held_mean: Option<Vec<f64>>,
    swing: Swing,
    skin: Option<SkinSummary>,
    resolved_left: Option<u32>,
    res_q: Vec<f64>,
    res_time: f64,
    res_min: Vec<f64>,
    res_max: Vec<f64>,
    resolved: Option<ResolvedSummary>,
}

/// Where a solid cell of the skin layer gets its swing from: the held
/// cell, the centroid and unit normal of its interface face, and the depth
/// of the held cell below that face (m).
type Source = (u32, [f64; 3], [f64; 3], f64);

/// The swing of the held cells: per cell the complex amplitudes of the
/// harmonics (K).
type Swing = Vec<(u32, Vec<(f64, f64)>)>;

/// One fluid–solid interface face.
#[derive(Debug, Clone)]
struct Face {
    slot: usize,
    face: u32,
    solid_owner: bool,
    pair: usize,
    area: f64,
    solid_cell: u32,
    /// The solid region of the face and its effusivity √(k ρ c_p).
    region: usize,
    effusivity: f64,
}

/// The averages of a fluid phase per face, for `cycle_faces.csv`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct FaceMeans {
    /// Ḡ, Ḡ_fluid (W/K), T̄_ref (K), Q̄ into the solid (W), T_w after the
    /// solid phase (K).
    g: Vec<f64>,
    gf: Vec<f64>,
    t_ref: Vec<f64>,
    q: Vec<f64>,
    tw: Vec<f64>,
    /// With the skin swing: the amplitude, minimum and maximum of the
    /// surface temperature over a cycle (K).
    amp: Vec<f64>,
    tmin: Vec<f64>,
    tmax: Vec<f64>,
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
    /// The harmonics of the skin swing (0 without it).
    harmonics: usize,
    /// The skin depth of the period √(2k/(ρ c_p ω)) per region (m; 0 for
    /// the fluid).
    skin_depth: Vec<f64>,
    /// The held temperatures of the solids (the mean, without the swing).
    held_mean: Option<Vec<f64>>,
    /// The swing of the held cells: per cell the complex amplitudes of the
    /// harmonics.
    swing: Swing,
    skin: Option<SkinSummary>,
    /// The fully coupled cycles after the loop: those still to run, the
    /// heat flow and the range of the wall temperature per face of the
    /// running cycle, and the summary.
    resolved_left: Option<u32>,
    res_q: Vec<f64>,
    res_time: f64,
    res_min: Vec<f64>,
    res_max: Vec<f64>,
    resolved: Option<ResolvedSummary>,
    /// The end time of the run (a transient solid phase does not march
    /// past it), and the jump of the clock that the last update asks for.
    end_time: Option<f64>,
    jump: Option<f64>,
}

impl Controller {
    /// The loop of a conjugate case from time `t0`: the solids are held at
    /// their current temperatures from the first step on.
    pub fn new(
        fl: &mut Flow,
        cfg: &CycleCoupling,
        t0: f64,
        dt: f64,
        dim: u8,
        end_time: Option<f64>,
    ) -> Result<Self> {
        let mut skin_notes = Vec::new();
        let mut skin = None;
        let (faces, pairs, solid, skin_depth) = {
            let (sc, dev) = crate::cht::energy_fields(fl).ok_or_else(|| {
                anyhow::anyhow!("[energy.cycle] needs the energy equation of a conjugate case")
            })?;
            let regions = sc.energy_cell_regions();
            let state = sc.energy_region_state(dev)?;
            // ρ c_p and k from the materials (the face values of the energy
            // equation are those of the fluid until its first refresh)
            let mat = sc.energy_region_material(dev)?;
            let k: Vec<f64> = mat.iter().map(|m| m.1).collect();
            let rho_cp: Vec<f64> = mat.iter().map(|m| m.0).collect();
            let effusivity: Vec<f64> = rho_cp
                .iter()
                .zip(&k)
                .map(|(&c, &k)| (k * c).sqrt())
                .collect();
            let faces: Vec<Face> = sc
                .energy_fluid_solid_faces()
                .into_iter()
                .map(|(slot, face, solid_owner, st)| {
                    let (o, nb) = sc.energy_interface_cells(slot);
                    let solid_cell = if solid_owner { o } else { nb };
                    let region = regions[solid_cell as usize] as usize;
                    Face {
                        slot,
                        face,
                        solid_owner,
                        pair: sc.energy_interface_pair(slot),
                        area: st[0],
                        solid_cell,
                        region,
                        effusivity: effusivity[region],
                    }
                })
                .collect();
            let pairs: Vec<String> = sc.energy_interfaces().into_iter().map(|x| x.0).collect();
            let names = sc.energy_regions();
            // the skin depth of the period against the thickness V/A of
            // every solid region at the interfaces (CYC003)
            if cfg.skin == "effusivity" {
                let omega = 2.0 * std::f64::consts::PI / cfg.period;
                let mut area = vec![0.0; state.len()];
                for f in &faces {
                    area[f.region] += f.area;
                }
                let mut regions = Vec::new();
                let mut thick_faces = 0;
                for (r, a) in area.iter().enumerate() {
                    if *a <= 0.0 || state[r].0 <= 0.0 {
                        continue;
                    }
                    let d = (2.0 * k[r] / (rho_cp[r] * omega)).sqrt();
                    let l = state[r].0 / a;
                    let name = names.get(r).cloned().unwrap_or_default();
                    if d > l / 3.0 {
                        let n = faces.iter().filter(|f| f.region == r).count();
                        thick_faces += n;
                        skin_notes.push(Diagnostic::new(
                            crate::diag::CYC_SKIN,
                            format!(
                                "[energy.cycle] skin = \"effusivity\": the skin depth of the \
                                 period in region `{name}` is {d:.3e} m, more than a third of its \
                                 thickness V/A = {l:.3e} m ({n} interface faces): the swing \
                                 reaches through the wall and the semi-infinite solution does not \
                                 hold there; resolved_cycles gives the exact swing"
                            ),
                        ));
                    }
                    regions.push((name, d, l));
                }
                skin = Some(SkinSummary {
                    harmonics: 0,
                    regions,
                    thick_faces,
                    swing_max: 0.0,
                });
            }
            let fluid = sc.energy_fluid_flags();
            let solid: Vec<bool> = sc
                .energy_cell_regions()
                .iter()
                .map(|&m| !fluid[m as usize])
                .collect();
            let t = sc.energy_temperature(dev)?;
            hold(sc, dev, &solid, &t)?;
            let omega = 2.0 * std::f64::consts::PI / cfg.period;
            let skin_depth: Vec<f64> = rho_cp
                .iter()
                .zip(&k)
                .map(|(&c, &k)| {
                    if c > 0.0 {
                        (2.0 * k / (c * omega)).sqrt()
                    } else {
                        0.0
                    }
                })
                .collect();
            (faces, pairs, solid, skin_depth)
        };
        anyhow::ensure!(
            !faces.is_empty(),
            "[energy.cycle]: the case has no interface faces between the fluid and a solid"
        );
        let n = faces.len();
        let steps_per_cycle = (cfg.period / dt).round() as usize;
        let harmonics = if cfg.skin == "effusivity" {
            (cfg.skin_harmonics as usize)
                .min(steps_per_cycle.saturating_sub(1) / 2)
                .max(1)
        } else {
            0
        };
        Ok(Self {
            steps_per_cycle,
            harmonics,
            skin_depth,
            held_mean: None,
            swing: Vec::new(),
            skin: skin.map(|s| SkinSummary { harmonics, ..s }),
            resolved_left: None,
            res_q: vec![0.0; n],
            res_time: 0.0,
            res_min: vec![f64::MAX; n],
            res_max: vec![f64::MIN; n],
            resolved: None,
            end_time,
            jump: None,
            cfg: cfg.clone(),
            dim,
            faces,
            pairs,
            solid,
            cycle_start: t0,
            cycles_in_phase: 0,
            cycles_total: 0,
            cur: Acc::new(n, harmonics),
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
            diagnostics: skin_notes,
        })
    }

    /// The state of the loop for a restart file (JSON).
    pub fn save(&self) -> Result<Vec<u8>> {
        let st = LoopState {
            faces: self.faces.len(),
            cycle_start: self.cycle_start,
            cycles_in_phase: self.cycles_in_phase,
            cycles_total: self.cycles_total,
            cur: self.cur.clone(),
            hist: self.hist.clone(),
            cycle_q: self.cycle_q.clone(),
            tw_prev: self.tw_prev.clone(),
            q_prev: self.q_prev,
            cycle_abs: self.cycle_abs,
            omega: self.omega,
            hist_r: self.hist_r.clone(),
            hist_xt: self.hist_xt.clone(),
            updates: self.updates.clone(),
            means: self.means.clone(),
            done: self.done,
            converged: self.converged,
            ended_at: self.ended_at,
            diagnostics: self
                .diagnostics
                .iter()
                .map(|d| (d.code.to_string(), d.message.clone()))
                .collect(),
            held_mean: self.held_mean.clone(),
            swing: self.swing.clone(),
            skin: self.skin.clone(),
            resolved_left: self.resolved_left,
            res_q: self.res_q.clone(),
            res_time: self.res_time,
            res_min: self.res_min.clone(),
            res_max: self.res_max.clone(),
            resolved: self.resolved.clone(),
        };
        Ok(serde_json::to_vec(&st)?)
    }

    /// Continues the loop from the state of a restart file (`save`): the
    /// solids are held again as they were (with the swing at time `t`), or
    /// left free in the fully coupled cycles at the end.
    pub fn restore(&mut self, fl: &mut Flow, bytes: &[u8], t: f64) -> Result<()> {
        let st: LoopState = serde_json::from_slice(bytes)
            .map_err(|e| anyhow::anyhow!("the loop state of the restart file: {e}"))?;
        anyhow::ensure!(
            st.faces == self.faces.len(),
            "the loop state of the restart file has {} interface faces, the case {}",
            st.faces,
            self.faces.len()
        );
        self.cycle_start = st.cycle_start;
        self.cycles_in_phase = st.cycles_in_phase;
        self.cycles_total = st.cycles_total;
        self.cur = st.cur;
        self.hist = st.hist;
        self.cycle_q = st.cycle_q;
        self.tw_prev = st.tw_prev;
        self.q_prev = st.q_prev;
        self.cycle_abs = st.cycle_abs;
        self.omega = st.omega;
        self.hist_r = st.hist_r;
        self.hist_xt = st.hist_xt;
        self.updates = st.updates;
        self.means = st.means;
        self.done = st.done;
        self.converged = st.converged;
        self.ended_at = st.ended_at;
        // the notes of the start (CYC003) come again from the case
        for (c, m) in st.diagnostics {
            if let Some(info) = crate::diag::info(&c)
                && !self
                    .diagnostics
                    .iter()
                    .any(|d| d.code == info.code && d.message == m)
            {
                self.diagnostics.push(Diagnostic::new(info.code, m));
            }
        }
        self.held_mean = st.held_mean;
        self.swing = st.swing;
        if st.skin.is_some() {
            self.skin = st.skin;
        }
        self.resolved_left = st.resolved_left;
        self.res_q = st.res_q;
        self.res_time = st.res_time;
        self.res_min = st.res_min;
        self.res_max = st.res_max;
        self.resolved = st.resolved;
        if self.resolved_left.is_some() || self.done {
            let (sc, dev) = fields(fl)?;
            sc.set_energy_frozen(dev, None)?;
        } else if !self.swing.is_empty() {
            self.apply_swing(fl, t)?;
        } else if let Some(h) = self.held_mean.clone() {
            let (sc, dev) = fields(fl)?;
            hold(sc, dev, &self.solid, &h)?;
        }
        Ok(())
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

    /// The jump of the run's clock after a transient solid phase (s, whole
    /// periods), once: the solids have marched past the end of the fluid
    /// phase.
    pub fn take_jump(&mut self) -> Option<f64> {
        self.jump.take()
    }

    /// After a time step that ended at `t` and took `dt`: accumulates the
    /// interface faces, closes a cycle, runs a solid phase and an update
    /// at the end of a fluid phase. Returns the lines for the log.
    pub fn after_step(&mut self, fl: &mut Flow, t: f64, dt: f64) -> Result<Vec<String>> {
        if self.done {
            return Ok(Vec::new());
        }
        if self.resolved_left.is_some() {
            return self.resolved_step(fl, t, dt);
        }
        let lines = self.advance(fl, t, dt)?;
        if !self.done && self.resolved_left.is_none() && !self.swing.is_empty() {
            self.apply_swing(fl, t + dt)?;
        }
        Ok(lines)
    }

    /// Releases the solids for the fully coupled cycles at the end of the
    /// loop, at time `t` with steps of `dt`. With the skin swing the solids
    /// start from the periodic state of a semi-infinite solid: every solid
    /// cell within six skin depths of a held cell takes the swing of that
    /// cell damped and delayed by its depth beyond it, s_n e^{−(1+i) d/δ_n}
    /// with δ_n = δ/√n, on the temperature and on its old time levels (at
    /// t − dt and t − 2 dt). Released from the mean alone, the solids would
    /// start with the energy of the missing swing in the skin, which leaves
    /// with the time constant of the whole wall, not of the skin.
    fn release(&mut self, fl: &mut Flow, t: f64, dt: f64) -> Result<()> {
        let profile = match (&self.held_mean, self.swing.is_empty()) {
            (Some(mean), false) => Some((mean.clone(), self.depth_layer(fl)?)),
            _ => None,
        };
        let (sc, dev) = fields(fl)?;
        sc.set_energy_frozen(dev, None)?;
        if let Some((mean, layer)) = profile {
            let (mut p, mut o1, mut o2) = sc.energy_levels(dev)?;
            let omega = 2.0 * std::f64::consts::PI / self.cfg.period;
            let regions = sc.energy_cell_regions();
            let at: std::collections::HashMap<u32, &Vec<(f64, f64)>> =
                self.swing.iter().map(|(c, a)| (*c, a)).collect();
            for (cell, held, d) in layer {
                let Some(amp) = at.get(&held) else { continue };
                let delta = self.skin_depth[regions[cell as usize] as usize];
                let mut v = [0.0; 3];
                for (n, &(re, im)) in amp.iter().enumerate() {
                    let dn = if delta > 0.0 {
                        d * ((n + 1) as f64).sqrt() / delta
                    } else {
                        0.0
                    };
                    // s_n e^{−(1+i) d/δ_n}
                    let (fr, fi) = ((-dn).exp() * dn.cos(), -(-dn).exp() * dn.sin());
                    let (ar, ai) = (re * fr - im * fi, re * fi + im * fr);
                    for (j, x) in v.iter_mut().enumerate() {
                        let w = (n + 1) as f64 * omega * (t - j as f64 * dt);
                        *x += ar * w.cos() - ai * w.sin();
                    }
                }
                let c = cell as usize;
                p[c] = mean[c] + v[0];
                o1[c] = mean[c] + v[1];
                o2[c] = mean[c] + v[2];
            }
            sc.set_energy_levels3(dev, &p, &o1, &o2)?;
        }
        self.swing.clear();
        Ok(())
    }

    /// The solid cells within six skin depths of the held cells: (cell,
    /// held cell, depth beyond the held cell along the normal of its
    /// interface face, m), by a breadth-first walk through the solid cells
    /// from the held cells (the held cells themselves with depth 0).
    fn depth_layer(&self, fl: &Flow) -> Result<Vec<(u32, u32, f64)>> {
        let c = fl
            .cht
            .as_deref()
            .ok_or_else(|| anyhow::anyhow!("[energy.cycle] needs a conjugate case"))?;
        let (mesh, geo) = (&c.view.mesh, &c.view.geo);
        let delta_max = self.skin_depth.iter().fold(0.0_f64, |a, &b| a.max(b));
        let reach = 6.0 * delta_max;
        let nc = mesh.n_cells;
        // per cell: the source (held cell, face centroid, unit normal, depth
        // of the held cell)
        let mut src: Vec<Option<Source>> = vec![None; nc];
        let mut queue = VecDeque::new();
        let mut out = Vec::new();
        let depth = |x: &[f64; 3], xf: &[f64; 3], nf: &[f64; 3]| -> f64 {
            ((x[0] - xf[0]) * nf[0] + (x[1] - xf[1]) * nf[1] + (x[2] - xf[2]) * nf[2]).abs()
        };
        for f in &self.faces {
            let h = f.solid_cell as usize;
            if src[h].is_some() || h >= nc {
                continue;
            }
            let fi = f.face as usize;
            let a = geo.face_area[fi];
            let m = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt().max(1e-300);
            let nf = [a[0] / m, a[1] / m, a[2] / m];
            let xf = geo.face_centroid[fi];
            let d0 = depth(&geo.cell_centroid[h], &xf, &nf);
            src[h] = Some((f.solid_cell, xf, nf, d0));
            out.push((f.solid_cell, f.solid_cell, 0.0));
            queue.push_back(h);
        }
        while let Some(cell) = queue.pop_front() {
            let (held, xf, nf, d0) = src[cell].expect("visited");
            let (a, b) = (
                mesh.cell_ptr[cell] as usize,
                mesh.cell_ptr[cell + 1] as usize,
            );
            for &face in &mesh.cell_faces[a..b] {
                let face = face as usize;
                if face >= mesh.n_inner_faces {
                    continue;
                }
                let (o, n) = (mesh.owner[face] as usize, mesh.neighbour[face] as usize);
                let other = if o == cell { n } else { o };
                if other >= nc || !self.solid[other] || src[other].is_some() {
                    continue;
                }
                let d = (depth(&geo.cell_centroid[other], &xf, &nf) - d0).max(0.0);
                if d > reach {
                    continue;
                }
                src[other] = Some((held, xf, nf, d0));
                out.push((other as u32, held, d));
                queue.push_back(other);
            }
        }
        Ok(out)
    }

    /// Holds the solids at their mean temperatures plus the skin swing at
    /// time `t` (the held cells at the interface faces).
    fn apply_swing(&mut self, fl: &mut Flow, t: f64) -> Result<()> {
        let Some(mean) = &self.held_mean else {
            return Ok(());
        };
        let omega = 2.0 * std::f64::consts::PI / self.cfg.period;
        let mut v = mean.clone();
        for (c, amp) in &self.swing {
            let mut d = 0.0;
            for (n, &(re, im)) in amp.iter().enumerate() {
                let w = (n + 1) as f64 * omega * t;
                d += re * w.cos() - im * w.sin();
            }
            v[*c as usize] += d;
        }
        let (sc, dev) = fields(fl)?;
        hold(sc, dev, &self.solid, &v)
    }

    /// A step of the fully coupled cycles after the loop.
    fn resolved_step(&mut self, fl: &mut Flow, t: f64, dt: f64) -> Result<Vec<String>> {
        let mut lines = Vec::new();
        {
            let (sc, _) = fields(fl)?;
            for (i, (_, _, _, x)) in sc.energy_fluid_solid_faces().iter().enumerate() {
                self.res_q[i] += x[3] * dt;
                self.res_min[i] = self.res_min[i].min(x[5]);
                self.res_max[i] = self.res_max[i].max(x[5]);
            }
        }
        self.res_time += dt;
        if t - self.cycle_start < self.cfg.period - 1e-6 * dt {
            return Ok(lines);
        }
        self.cycle_start = t;
        self.cycles_total += 1;
        let left = self.resolved_left.unwrap_or(1).saturating_sub(1);
        let mut per = vec![0.0; self.pairs.len()];
        for (f, q) in self.faces.iter().zip(&self.res_q) {
            per[f.pair] -= q / self.res_time;
        }
        let swing = self
            .res_min
            .iter()
            .zip(&self.res_max)
            .fold(0.0_f64, |a, (lo, hi)| a.max(0.5 * (hi - lo)));
        let done = self.cfg.resolved_cycles - left;
        lines.push(format!(
            "[energy.cycle] fully coupled cycle {done} of {} at t = {t} s: interface heat flow \
             {:.6e} W into the fluid, largest wall swing {swing:.4} K",
            self.cfg.resolved_cycles,
            per.iter().sum::<f64>()
        ));
        if left == 0 {
            self.resolved = Some(ResolvedSummary {
                cycles: done,
                heat_flow: self.pairs.iter().cloned().zip(per).collect(),
                swing_max: swing,
            });
            self.resolved_left = None;
            self.done = true;
            self.ended_at = Some(t);
            lines.push("[energy.cycle] the fully coupled cycles are done: the run stops".into());
        } else {
            self.resolved_left = Some(left);
            self.res_q.iter_mut().for_each(|x| *x = 0.0);
            self.res_min.iter_mut().for_each(|x| *x = f64::MAX);
            self.res_max.iter_mut().for_each(|x| *x = f64::MIN);
            self.res_time = 0.0;
        }
        Ok(lines)
    }

    /// A step of a fluid phase: accumulates the interface faces, closes a
    /// cycle, runs a solid phase and an update at the end of the phase.
    fn advance(&mut self, fl: &mut Flow, t: f64, dt: f64) -> Result<Vec<String>> {
        let mut lines = Vec::new();
        {
            let (sc, _) = fields(fl)?;
            let st = sc.energy_fluid_solid_faces();
            let h = self.harmonics;
            let omega = 2.0 * std::f64::consts::PI / self.cfg.period;
            let c = &mut self.cur;
            c.time += dt;
            for (i, (_, _, _, x)) in st.iter().enumerate() {
                let [area, g, gf, q, ts, tw] = *x;
                c.g[i] += g * dt;
                c.gts[i] += g * ts * dt;
                c.q[i] += q * dt;
                c.gf[i] += gf * dt;
                c.tw[i] += tw * dt;
                if h > 0 {
                    c.sq[i] += q * q * dt;
                    c.stq[i] += (tw - ts) * q * dt;
                    let qa = if area > 0.0 { q / area } else { 0.0 };
                    for n in 0..h {
                        let w = (n + 1) as f64 * omega * t;
                        c.fr[i * h + n] += qa * w.cos() * dt;
                        c.fi[i * h + n] -= qa * w.sin() * dt;
                    }
                }
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
        let done = std::mem::replace(&mut self.cur, Acc::new(self.faces.len(), self.harmonics));
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
        lines.extend(self.update(fl, t, dt, periodic, change)?);
        Ok(lines)
    }

    /// The solid phase and the update at the end of a fluid phase.
    fn update(
        &mut self,
        fl: &mut Flow,
        t: f64,
        dt: f64,
        periodic: bool,
        change: f64,
    ) -> Result<Vec<String>> {
        let mut lines = Vec::new();
        let n = self.faces.len();
        let mut sum = Acc::new(n, self.harmonics);
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
        // a transient solid phase: the march covers the time since the
        // solids were last updated (the fluid phase, in which they were
        // held) and the rest of `march_cycles` periods, over which the
        // clock then jumps (whole periods, so the flow continues in phase),
        // not past the end of the run
        let transient = self.cfg.solid == "transient";
        let span_phase = f64::from(self.cycles_in_phase) * self.cfg.period;
        let (march, jump) = if transient {
            let want = f64::from(self.cfg.march_cycles.unwrap_or(0)) * self.cfg.period;
            let mut m = want.max(span_phase);
            if let Some(e) = self.end_time {
                m = m.min(e - (t - span_phase)).max(span_phase);
            }
            (Some(m), Some((m - span_phase).max(0.0)))
        } else {
            (None, None)
        };
        let march_dt = march.map(|m| m / f64::from(self.cfg.march_steps.max(1)));
        // the solid phase, from the state at the end of the fluid phase
        // (its time levels are put back afterwards)
        let saved = {
            let (sc, dev) = fields(fl)?;
            sc.energy_levels(dev)?
        };
        // the solids at their held mean (without the skin swing)
        let mut t_before = saved.0.clone();
        if let Some(h) = &self.held_mean {
            for (c, &s) in self.solid.iter().enumerate() {
                if s {
                    t_before[c] = h[c];
                }
            }
        }
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
            let r = match march_dt {
                Some(h) => {
                    sc.set_energy_temperature(dev, &t_before)?;
                    march_steps(fl, None, t, h, self.cfg.march_steps, span)?
                }
                None => steady(fl, None, span)?,
            };
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
            let (last, i2) = match march_dt {
                Some(h) => march_steps(fl, Some(&mf), t, h, self.cfg.march_steps, span)?,
                None => steady(fl, Some(&mf), span)?,
            };
            (last, i1 + i2)
        };
        // the balance of the solids in the solid phase: the fluid counted
        // as held (what flows into it, convection included, is booked as
        // `held`, so the balance is that of the solids with the Robin heat
        // or the defect), the old levels at the solved field (a steady
        // solve stores nothing; a march keeps those of its last step, whose
        // balance it is)
        {
            let (sc, dev) = fields(fl)?;
            if !transient {
                sc.set_energy_levels(dev, &last)?;
            }
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
        //
        // With the skin swing the unknowns also hold the swing of the held
        // cells (real and imaginary part of every harmonic): the swing the
        // next fluid phase sees comes from the flux of this one, which the
        // swing itself changes, so the same update takes both.
        let (amp, tmin, tmax, swing_max, s_new) = if self.harmonics > 0 {
            let (a, lo, hi, sw) = self.skin_swing(&sum, tm, &tw);
            let m = a.iter().fold(0.0_f64, |x, &y| x.max(y));
            if let Some(sk) = &mut self.skin {
                sk.swing_max = m;
            }
            (a, lo, hi, Some(m), sw)
        } else {
            (Vec::new(), Vec::new(), Vec::new(), None, Vec::new())
        };
        let flat = |sw: &Swing| -> Vec<f64> {
            sw.iter()
                .flat_map(|(_, v)| v.iter().flat_map(|&(a, b)| [a, b]))
                .collect()
        };
        let s_old: Vec<f64> = if self.swing.is_empty() {
            vec![0.0; flat(&s_new).len()]
        } else {
            flat(&self.swing)
        };
        let s_next = flat(&s_new);
        let swing_change = s_next
            .chunks(2)
            .zip(s_old.chunks(2))
            .fold(0.0_f64, |a, (p, q)| a.max((p[0] - q[0]).hypot(p[1] - q[1])));
        let idx: Vec<usize> = (0..self.solid.len()).filter(|&c| self.solid[c]).collect();
        let nx = idx.len();
        let x: Vec<f64> = idx
            .iter()
            .map(|&c| t_before[c])
            .chain(s_old.iter().copied())
            .collect();
        let xt: Vec<f64> = idx
            .iter()
            .map(|&c| last[c])
            .chain(s_next.iter().copied())
            .collect();
        let r: Vec<f64> = xt.iter().zip(&x).map(|(a, b)| a - b).collect();
        // (a transient solid phase: the march is a history, not a fixed
        // point; its result is the new state)
        let (x_new, w) = match (self.cfg.update.as_str(), self.hist_r.is_empty()) {
            _ if transient => (xt.clone(), 1.0),
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
        if !transient {
            self.hist_r.push(r);
            self.hist_xt.push(xt);
        }
        // the fluid back at its state and time levels of the end of the
        // fluid phase, the solids at the new temperatures on every level
        if !s_new.is_empty() {
            let mut k = nx;
            let mut sw = s_new;
            for (_, v) in sw.iter_mut() {
                for z in v.iter_mut() {
                    *z = (x_new[k], x_new[k + 1]);
                    k += 2;
                }
            }
            self.swing = sw;
        }
        let (mut t_new, mut o1, mut o2) = saved;
        for (&c, &v) in idx.iter().zip(&x_new[..nx]) {
            t_new[c] = v;
            o1[c] = v;
            o2[c] = v;
        }
        sc.set_energy_levels3(dev, &t_new, &o1, &o2)?;
        hold(sc, dev, &self.solid, &t_new)?;
        self.held_mean = Some(t_new);
        let reference = self.tw_prev.as_ref().unwrap_or(&tw_fluid);
        let area: f64 = self.faces.iter().map(|f| f.area).exact_sum();
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
            swing_max,
            swing_change: swing_max.map(|_| swing_change),
            march_time: march,
            march_jump: jump,
        };
        lines.push(format!(
            "[energy.cycle] update {k} at t = {t} s after {} cycles: interface heat flow {:.6e} W \
             into the fluid{}, wall temperature change max {:.4} K, mean {:.4} K, mean wall {:.4} \
             K; solid phase {its} solves, balance {:.2e}, relaxation {:.3}{}",
            self.cycles_in_phase,
            q_total,
            heat_change.map_or(String::new(), |c| format!(" (change {c:.2e})")),
            dmax,
            dmean,
            wmean,
            bal.relative,
            w,
            match (march, jump) {
                (Some(m), Some(j)) => format!(
                    "; the solids marched {m} s in {} steps, the clock jumps {j} s to t = {} s",
                    self.cfg.march_steps,
                    t + j
                ),
                _ => String::new(),
            }
        ));
        self.updates.push(row);
        self.means = FaceMeans {
            g,
            gf,
            t_ref,
            q,
            tw: tw.clone(),
            amp,
            tmin,
            tmax,
        };
        let converged = heat_change.is_some_and(|c| c <= self.cfg.heat_tolerance)
            && dmax <= self.cfg.tolerance
            && swing_change <= self.cfg.tolerance;
        self.tw_prev = Some(tw);
        self.q_prev = Some(q_total);
        let resolved = self.cfg.resolved_cycles;
        let end = |c: &mut Self, lines: &mut Vec<String>, fl: &mut Flow| -> Result<()> {
            if resolved > 0 {
                // release the solids: fully coupled cycles from here
                c.release(fl, t, dt)?;
                c.resolved_left = Some(resolved);
                lines.push(format!(
                    "[energy.cycle] {resolved} fully coupled cycles follow from the loop's solid \
                     temperatures"
                ));
            } else {
                c.done = true;
                c.ended_at = Some(t);
            }
            Ok(())
        };
        if converged {
            self.converged = true;
            lines.push(format!(
                "[energy.cycle] converged after {k} updates and {} cycles{}",
                self.cycles_total,
                if resolved > 0 { "" } else { ": the run stops" }
            ));
            end(self, &mut lines, fl)?;
        } else if k >= self.cfg.max_updates {
            let m = format!(
                "[energy.cycle] stopped at max_updates = {k} before the tolerances: the last \
                 wall temperature change {dmax:.4} K (tolerance {} K), interface heat flow change \
                 {} (heat_tolerance {}){}",
                self.cfg.tolerance,
                heat_change.map_or("-".to_string(), |c| format!("{c:.2e}")),
                self.cfg.heat_tolerance,
                if self.harmonics > 0 {
                    format!(", skin swing change {swing_change:.4} K (tolerance)")
                } else {
                    String::new()
                }
            );
            self.diagnostics
                .push(Diagnostic::new(crate::diag::CYC_MAX_UPDATES, m.clone()));
            lines.push(m);
            end(self, &mut lines, fl)?;
        }
        // the next fluid phase (after the jump of a transient solid phase)
        if let Some(j) = jump.filter(|&j| j > 0.0) {
            self.jump = Some(j);
            self.cycle_start += j;
            if let Some(e) = self.ended_at.as_mut() {
                *e += j;
            }
        }
        self.cycles_in_phase = 0;
        self.hist.clear();
        self.cycle_q.clear();
        self.cur = Acc::new(n, self.harmonics);
        Ok(lines)
    }

    /// The skin swing from the harmonics of the face flux of the averaging
    /// cycles `sum` (over `tm` s), with the mean wall temperatures `tw`:
    /// sets the swing of the held cells for the next fluid phase and gives
    /// per face the amplitude, minimum and maximum of the surface
    /// temperature over a cycle.
    #[allow(clippy::type_complexity)]
    fn skin_swing(&self, sum: &Acc, tm: f64, tw: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Swing) {
        const SAMPLES: usize = 128;
        let h = self.harmonics;
        let omega = 2.0 * std::f64::consts::PI / self.cfg.period;
        let n = self.faces.len();
        let (mut amp, mut lo, mut hi) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        // per held cell: Σ area · s_n and Σ area
        let mut cells: std::collections::BTreeMap<u32, (f64, Vec<(f64, f64)>)> =
            std::collections::BTreeMap::new();
        for (i, f) in self.faces.iter().enumerate() {
            // the resistance of the solid half cell (K m²/W) from
            // T_w − T_s* = r_s q/A by least squares over the cycles
            let rs = if sum.sq[i] > 0.0 {
                f.area * sum.stq[i] / sum.sq[i]
            } else {
                0.0
            };
            let mut theta = Vec::with_capacity(h);
            let entry = cells
                .entry(f.solid_cell)
                .or_insert_with(|| (0.0, vec![(0.0, 0.0); h]));
            entry.0 += f.area;
            for k in 0..h {
                let w = (k + 1) as f64 * omega;
                let (cr, ci) = (2.0 * sum.fr[i * h + k] / tm, 2.0 * sum.fi[i * h + k] / tm);
                // Z = (1 − i)/(√2 e √ω), the surface response of a
                // semi-infinite solid to a periodic flux (K m²/W)
                let m = if f.effusivity > 0.0 {
                    1.0 / (std::f64::consts::SQRT_2 * f.effusivity * w.sqrt())
                } else {
                    0.0
                };
                let (zr, zi) = (m, -m);
                theta.push((zr * cr - zi * ci, zr * ci + zi * cr));
                // the held cell: (Z − r_s) q̂
                let (yr, yi) = (zr - rs, zi);
                entry.1[k].0 += f.area * (yr * cr - yi * ci);
                entry.1[k].1 += f.area * (yr * ci + yi * cr);
            }
            let (mut a, mut b) = (f64::INFINITY, f64::NEG_INFINITY);
            for j in 0..SAMPLES {
                let tau = j as f64 / SAMPLES as f64 * self.cfg.period;
                let v: f64 = theta
                    .iter()
                    .enumerate()
                    .map(|(k, &(re, im))| {
                        let x = (k + 1) as f64 * omega * tau;
                        re * x.cos() - im * x.sin()
                    })
                    .sum();
                a = a.min(v);
                b = b.max(v);
            }
            amp[i] = 0.5 * (b - a);
            lo[i] = tw[i] + a;
            hi[i] = tw[i] + b;
        }
        let swing = cells
            .into_iter()
            .map(|(c, (a, v))| (c, v.into_iter().map(|(x, y)| (x / a, y / a)).collect()))
            .collect();
        (amp, lo, hi, swing)
    }

    /// The loop for the result at the end of the run; a run that ended
    /// before the loop gets CYC004 (not with a transient solid phase,
    /// whose history ends with the run).
    pub fn summary(&mut self) -> CycleSummary {
        if !self.done && self.cfg.solid != "transient" {
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
            solid: self.cfg.solid.clone(),
            skin: self.skin.clone(),
            resolved: self.resolved.clone(),
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
             relaxation,solid_imbalance,solid_relative,swing_max,swing_change,march_time,\
             march_jump",
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
                    ",{},{:e},{:e},{:e},{},{:e},{:e},{:e},{},{},{},{}",
                    u.heat_change.map_or(String::new(), |c| format!("{c:e}")),
                    u.wall_change_max,
                    u.wall_change_mean,
                    u.wall_mean,
                    u.solid_iterations,
                    u.relaxation,
                    u.solid_imbalance,
                    u.solid_relative,
                    u.swing_max.map_or(String::new(), |x| format!("{x:e}")),
                    u.swing_change.map_or(String::new(), |x| format!("{x:e}")),
                    u.march_time.map_or(String::new(), |x| format!("{x:e}")),
                    u.march_jump.map_or(String::new(), |x| format!("{x:e}"))
                ));
                r
            })
            .collect();
        (head, rows)
    }

    /// The header and rows of `cycle_faces.csv`: per fluid–solid interface
    /// face of the last update, with the face centroids `centroid`.
    pub fn faces_csv(&self, centroid: &[[f64; 3]]) -> (String, Vec<String>) {
        let head = "face,interface,x,y,z,area,h_fluid,conductance,t_ref,heat_flux,t_wall,\
                    t_wall_amplitude,t_wall_min,t_wall_max,resolved_t_wall_min,resolved_t_wall_max"
            .into();
        let opt = |v: &[f64], i: usize| v.get(i).map_or(String::new(), |x| format!("{x:e}"));
        let resolved = self.resolved.is_some();
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
                        "{},{},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{},{},{},{},{}",
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
                        m.tw[i],
                        opt(&m.amp, i),
                        opt(&m.tmin, i),
                        opt(&m.tmax, i),
                        if resolved {
                            opt(&self.res_min, i)
                        } else {
                            String::new()
                        },
                        if resolved {
                            opt(&self.res_max, i)
                        } else {
                            String::new()
                        },
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
        let (top, rest) = a.split_at_mut(k + 1);
        let pivot = &top[k];
        for (i, row) in rest.iter_mut().enumerate() {
            let f = row[k] / d;
            for (x, p) in row[k..].iter_mut().zip(&pivot[k..]) {
                *x -= f * p;
            }
            b[k + 1 + i] -= f * b[k];
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

/// The march of a transient solid phase: `steps` steps of `dt` of the
/// energy equation (with the full-mesh face mass flux `flux` of the fluid
/// view in place of the current one, if given) from the current
/// temperatures, BDF1 first and BDF2 after, each solved until the largest
/// change of a temperature is below 10⁻⁹ of `span` (at most
/// `MARCH_ITERATIONS`). The boundary values and sources stay at those of
/// the time `t` of the update (their swing over a cycle is in the defect
/// and in the Robin pairs). The temperatures at the end and the solves.
fn march_steps(
    fl: &mut Flow,
    flux: Option<&[f64]>,
    t: f64,
    dt: f64,
    steps: u32,
    span: f64,
) -> Result<(Vec<f64>, usize)> {
    let mut last = {
        let (sc, dev) = fields(fl)?;
        let t0 = sc.energy_temperature(dev)?;
        sc.set_energy_levels3(dev, &t0, &t0, &t0)?;
        t0
    };
    let mut its = 0;
    for k in 0..steps.max(1) {
        let second = k > 0;
        {
            let (sc, dev) = fields(fl)?;
            sc.step_start(dev, true, second.then_some(1.0))?;
        }
        let st = crate::scalar::StepTime {
            time: t,
            dt,
            a: crate::scalar::bdf(second, 1.0),
        };
        for _ in 0..MARCH_ITERATIONS {
            let mut res = crate::simplec::Residuals::default();
            crate::cht::solve_with_flux(fl, Some(st), &mut res, false, flux)?;
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
