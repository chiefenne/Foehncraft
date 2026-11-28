//! The water and steam circuits behind the walls of a steady flow (steam
//! furnace stage SF7c, D350): the 1D tube circuits of SF7b
//! (`foehncraft-circuits`) coupled to the wall faces of the 3D energy
//! equation.
//!
//! Each `[[circuits.walls]]` entry maps the faces of one wall boundary to
//! the segments of one tube of the circuit: the face centroid is projected
//! on the tube axis from `start` (the tube's `from` end) to `end`, and the
//! fraction s ∈ [0, 1] along it picks the segment. A tube may take the
//! faces of several boundaries; a bank of `count` parallel tubes shares the
//! faces of its height.
//!
//! **Implicit coupling.** The water side enters the 3D energy equation as
//! a Robin condition on every mapped face,
//!
//!   q″ = U (T_water − T_b),   1/U = R_wall + φ (R_fouling + 1/h_in),
//!
//! with T_water and h_in the bulk temperature and the inner heat transfer
//! coefficient of the face's segment (Gnielinski, Chen or Gungor–Winterton,
//! post-CHF: the correlations of SF7b), φ the ratio of the 3D face area of
//! the segment to the inner wall area of its tubes, `wall_resistance` the
//! metal between the 3D face and the inner wall (per 3D area; zero when the
//! tube metal is a solid region of the mesh), `fouling` the scale on the
//! inner wall (per inner area). The condition is part of the energy matrix
//! (the D160 face form, as `heat_transfer_coefficient` with
//! `ambient_temperature`), so the wall temperature and heat follow the gas
//! within each outer iteration; only U and T_water lag. Every `every`
//! outer iterations the heat the 3D faces send into each segment is summed
//! (exact sums in face order), relaxed and handed to the circuit, whose
//! Newton solve starts from its last solution. Explicit, lagged heat-flux
//! or temperature exchange, which oscillates or needs heavy relaxation, is
//! the flaw this avoids (development path, section 8).
//!
//! **Convergence.** The run converges when the flow residuals are below
//! `[solver] tolerance` and the coupling residual, Σ|Q_3D − Q_circuit| /
//! Σ|Q_3D| over the segments (the heat the faces send now against the
//! heat of the last circuit solve), is below `[circuits] tolerance`; when
//! the flow converges first, the circuit is solved at once with the new
//! heats and the iterations go on.
//!
//! **Ledger.** The circuit books its own thermal account (region
//! `circuit`): the heat from the walls (one side of the exchange pair
//! [`PAIR_WALLS`]), the heat of tubes that are not coupled, and the
//! enthalpy flows of feedwater, steam and the pressure nodes. The other side
//! of the pair is the conduction of the coupled boundaries in the thermal
//! account of the 3D domain, so the pair residual is the coupling residual
//! in watts.

use crate::balance::{Account, Term, TermKind, names};
use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use anyhow::Result;
use foehncraft_circuits::coupling::{Coupled, WaterSide};
use foehncraft_core::exact::exact_sum;
use foehncraft_mesh::Mesh;
use foehncraft_mesh::geometry::Geometry;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

/// The exchange pair between the coupled 3D walls and the circuit.
pub const PAIR_WALLS: u32 = 2901;

/// `[circuits]`: the water and steam circuit behind the walls (SF7c).
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CircuitsSection {
    /// The circuit case file (SF7b, the format of the tool `circuits.run`), relative
    /// to the flow case. The heat of a coupled tube in it is the start
    /// value of the coupling; the walls replace it at the first update.
    pub file: PathBuf,
    /// Outer iterations between two circuit solves (default 10).
    #[serde(default = "default_every")]
    pub every: usize,
    /// The coupling residual Σ|Q_3D − Q_circuit| / Σ|Q_3D| at which the
    /// coupling counts as converged (default 1e-3).
    #[serde(default = "default_tolerance")]
    pub tolerance: f64,
    /// Relaxation α of the segment heats handed to the circuit, Q ← Q +
    /// α (Q_3D − Q), in (0, 1] (default 1: the Robin condition is
    /// implicit, the heat needs no damping in most cases).
    #[serde(default = "default_relaxation")]
    pub relaxation: f64,
    /// The walls and their tubes.
    pub walls: Vec<WallMap>,
}

fn default_every() -> usize {
    10
}
fn default_tolerance() -> f64 {
    1e-3
}
fn default_relaxation() -> f64 {
    1.0
}

/// `[[circuits.walls]]`: one wall boundary on one tube of the circuit.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WallMap {
    /// The wall boundary of the mesh.
    pub boundary: String,
    /// The tube of the circuit (a bank of `count` parallel tubes).
    pub tube: String,
    /// The point of the tube's `from` end (m), in the case coordinates.
    pub start: [f64; 3],
    /// The point of the tube's `to` end (m).
    pub end: [f64; 3],
    /// Fouling resistance of the inner wall (m² K/W, per inner wall area;
    /// default 0).
    #[serde(default)]
    pub fouling: f64,
    /// Thermal resistance of the metal between the 3D face and the inner
    /// wall (m² K/W, per 3D face area; default 0: the face is the inner
    /// wall, or the tube metal is a solid region of the mesh).
    #[serde(default)]
    pub wall_resistance: f64,
}

impl CircuitsSection {
    /// Makes the circuit file path relative to the case folder `dir`.
    pub fn resolve_paths(&mut self, dir: &Path) {
        if self.file.is_relative() {
            self.file = dir.join(&self.file);
        }
    }
}

fn bad(key: &str, msg: impl Into<String>) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::CIR_COUPLING, msg.into()).key(key.to_string()))
}

/// The names of the boundaries the circuits drive (empty without
/// `[circuits]`).
pub fn boundaries(case: &FlowCase) -> Vec<&str> {
    case.circuits
        .as_ref()
        .map(|c| c.walls.iter().map(|w| w.boundary.as_str()).collect())
        .unwrap_or_default()
}

/// The checks of `[circuits]` (CIR005); reads the circuit file.
pub fn validate(case: &FlowCase) -> Result<()> {
    let Some(c) = &case.circuits else {
        return Ok(());
    };
    if case.energy.is_none() {
        return Err(bad(
            "circuits",
            "[circuits] needs the energy equation ([energy]): the walls exchange heat with the \
             water and steam",
        ));
    }
    if case.time.is_some() {
        return Err(bad(
            "circuits",
            "[circuits] couples steady runs only: transient circuits (drum level, swell and \
             shrink) with a transient gas side are stage SF18",
        ));
    }
    if c.every == 0 {
        return Err(bad("circuits.every", "must be ≥ 1"));
    }
    if !(c.tolerance > 0.0 && c.tolerance.is_finite()) {
        return Err(bad("circuits.tolerance", "must be a number > 0"));
    }
    if !(c.relaxation > 0.0 && c.relaxation <= 1.0) {
        return Err(bad("circuits.relaxation", "must lie in (0, 1]"));
    }
    if c.walls.is_empty() {
        return Err(bad(
            "circuits.walls",
            "[circuits] needs at least one [[circuits.walls]] entry",
        ));
    }
    let circuit = foehncraft_circuits::case::Case::from_file(&c.file)
        .map_err(|e| bad("circuits.file", format!("{}: {e}", c.file.display())))?;
    let mut seen = std::collections::BTreeSet::new();
    for (i, w) in c.walls.iter().enumerate() {
        let key = format!("circuits.walls[{i}]");
        if !seen.insert(w.boundary.as_str()) {
            return Err(bad(
                &format!("{key}.boundary"),
                format!(
                    "boundary `{}` is on two tubes: give each wall boundary one tube (split the \
                     boundary in the mesh for several tube groups)",
                    w.boundary
                ),
            ));
        }
        match case.boundary.get(&w.boundary) {
            None => {
                return Err(bad(
                    &format!("{key}.boundary"),
                    format!("`{}` is not a boundary of the case", w.boundary),
                ));
            }
            Some(b) if crate::energy::side(b) != crate::energy::Side::Wall => {
                return Err(bad(
                    &format!("{key}.boundary"),
                    format!("`{}` is not a wall", w.boundary),
                ));
            }
            Some(_) => {}
        }
        if case.boundary_energy.contains_key(&w.boundary) {
            return Err(bad(
                &format!("{key}.boundary"),
                format!(
                    "[boundary.{}.energy] and the circuit both set the thermal condition of the \
                     wall: remove the energy table, the circuit gives it",
                    w.boundary
                ),
            ));
        }
        if !circuit.tubes.iter().any(|t| t.name == w.tube) {
            return Err(bad(
                &format!("{key}.tube"),
                format!("no tube `{}` in {}", w.tube, c.file.display()),
            ));
        }
        let d: [f64; 3] = std::array::from_fn(|k| w.end[k] - w.start[k]);
        if d.iter().map(|v| v * v).sum::<f64>() <= 0.0
            || w.start.iter().chain(&w.end).any(|v| !v.is_finite())
        {
            return Err(bad(
                &format!("{key}.end"),
                "the tube axis needs two different finite points `start` and `end`",
            ));
        }
        for (name, v) in [
            ("fouling", w.fouling),
            ("wall_resistance", w.wall_resistance),
        ] {
            if !(v >= 0.0 && v.is_finite()) {
                return Err(bad(&format!("{key}.{name}"), "must be a number ≥ 0"));
            }
        }
    }
    Ok(())
}

/// One mapped boundary face.
#[derive(Debug, Clone, Copy)]
struct Face {
    /// Boundary face index of the energy equation.
    fb: usize,
    /// The `[[circuits.walls]]` entry.
    wall: usize,
    /// Tube of the circuit.
    tube: usize,
    /// Segment of the tube.
    seg: usize,
    /// Position along the tube in segments, s n (0 at `from`, n at `to`).
    pos: f64,
    /// Area (m²).
    area: f64,
}

/// One circuit solve of the run (a row of `circuit.csv`).
#[derive(Debug, Clone, Copy, Serialize, JsonSchema)]
pub struct Update {
    /// The outer iteration before which it ran (0: the start).
    pub iteration: usize,
    /// Coupling residual of the heats it took (NaN for the start).
    pub residual: f64,
    /// Heat of the coupled tubes handed to the circuit (W).
    pub heat: f64,
    /// Newton iterations.
    pub newton: usize,
    /// The solve failed and the last water side was kept.
    pub failed: bool,
}

/// A drum of the circuit at the end.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct DrumReport {
    /// Node name.
    pub name: String,
    /// Steam flow (kg/s).
    pub steam_flow: f64,
    /// Circulation ratio.
    pub circulation_ratio: f64,
}

/// A coupled tube at the end.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct TubeReport {
    /// Tube name.
    pub name: String,
    /// Heat from the walls (W, all parallel tubes).
    pub heat: f64,
    /// Flow (kg/s, all parallel tubes).
    pub flow: f64,
    /// Outlet temperature (K).
    pub outlet_temperature: f64,
    /// Equilibrium quality of the last segment.
    pub outlet_quality: f64,
    /// Highest inner wall temperature (K).
    pub wall_temperature_max: f64,
    /// Smallest CHF ratio (NaN where not computed).
    pub chf_ratio_min: f64,
}

/// The coupled circuit in `result.json` (SF7c).
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct CircuitSummary {
    /// The circuit file.
    pub file: String,
    /// Circuit solves (the start included).
    pub solves: usize,
    /// Solves that failed (the water side of the solve before kept).
    pub failed: usize,
    /// Coupling residual at the end: Σ|Q_3D − Q_circuit| / Σ|Q_3D|.
    pub residual: f64,
    /// The coupling converged (`residual` ≤ `[circuits] tolerance`).
    pub converged: bool,
    /// Heat the 3D walls send into the water at the end (W).
    pub wall_heat: f64,
    /// Heat of the last circuit solve from the walls (W).
    pub circuit_heat: f64,
    /// Energy closure of the circuit: (in + heat − out) / heat.
    pub circuit_closure: f64,
    /// The drums.
    pub drums: Vec<DrumReport>,
    /// The coupled tubes.
    pub tubes: Vec<TubeReport>,
    /// Every circuit solve.
    pub updates: Vec<Update>,
}

/// The coupled circuit of a run.
#[derive(Debug)]
pub struct Walls {
    sec: CircuitsSection,
    circuit: Coupled,
    faces: Vec<Face>,
    /// 3D face area per tube segment (m²).
    area3: Vec<Vec<f64>>,
    /// The heats of the last circuit solve, per tube segment (W).
    applied: Vec<Vec<f64>>,
    /// The water side of the last successful solve.
    water: Vec<Vec<WaterSide>>,
    /// U (W/(m² K)) and T_water (K) per mapped face, as applied.
    robin: Vec<(f64, f64)>,
    updates: Vec<Update>,
    first: usize,
    /// The coupling residual of the last test.
    residual: f64,
    /// Failure messages not yet logged.
    notes: Vec<String>,
}

impl Walls {
    /// The coupling of `case` on the energy mesh `mesh` (the full mesh of a
    /// conjugate case), or `None` without `[circuits]`; `first` is the
    /// outer iteration the run starts from.
    pub fn create(
        case: &FlowCase,
        mesh: &Mesh,
        geo: &Geometry,
        first: usize,
    ) -> Result<Option<Self>> {
        let Some(sec) = case.circuits.clone() else {
            return Ok(None);
        };
        let ccase = foehncraft_circuits::case::Case::from_file(&sec.file)
            .map_err(|e| bad("circuits.file", format!("{}: {e}", sec.file.display())))?;
        let tubes: Vec<&str> = sec.walls.iter().map(|w| w.tube.as_str()).collect();
        let circuit = Coupled::new(ccase, &tubes).map_err(|e| bad("circuits.walls", e.message))?;
        let ni = mesh.n_interior_faces();
        let nt = circuit.model.tubes.len();
        let mut area3: Vec<Vec<f64>> = (0..nt).map(|i| vec![0.0; circuit.segments(i)]).collect();
        let mut faces = Vec::new();
        for (wi, w) in sec.walls.iter().enumerate() {
            let key = format!("circuits.walls[{wi}].boundary");
            let p = mesh
                .boundaries
                .iter()
                .find(|b| b.name == w.boundary)
                .ok_or_else(|| bad(&key, format!("no boundary `{}` in the mesh", w.boundary)))?;
            let ti = circuit.tube(&w.tube).unwrap_or(0);
            let n = circuit.segments(ti);
            let d: [f64; 3] = std::array::from_fn(|k| w.end[k] - w.start[k]);
            let dd: f64 = d.iter().map(|v| v * v).sum();
            for f in p.range() {
                let c = geo.face_centroid[f];
                let s = (0..3).map(|k| (c[k] - w.start[k]) * d[k]).sum::<f64>() / dd;
                if !(-1e-6..=1.0 + 1e-6).contains(&s) {
                    return Err(bad(
                        &key,
                        format!(
                            "face {f} of `{}` at ({:.4}, {:.4}, {:.4}) m lies outside the tube \
                             axis from start to end (s = {s:.4} of 0 to 1): `start` and `end` \
                             must span the whole boundary",
                            w.boundary, c[0], c[1], c[2]
                        ),
                    ));
                }
                let seg = ((s.clamp(0.0, 1.0) * n as f64).floor() as usize).min(n - 1);
                let area = geo.normal_area(f).1;
                area3[ti][seg] += area;
                faces.push(Face {
                    fb: f - ni,
                    wall: wi,
                    tube: ti,
                    seg,
                    pos: s.clamp(0.0, 1.0) * n as f64,
                    area,
                });
            }
        }
        faces.sort_by_key(|f| f.fb);
        let applied = (0..nt).map(|i| circuit.heat(i)).collect();
        Ok(Some(Self {
            robin: vec![(0.0, 0.0); faces.len()],
            sec,
            circuit,
            faces,
            area3,
            applied,
            water: Vec::new(),
            updates: Vec::new(),
            first,
            residual: f64::NAN,
            notes: Vec::new(),
        }))
    }

    /// The flow that holds the energy equation (the energy view of a
    /// conjugate case).
    fn energy(fl: &mut crate::simplec::Flow) -> &mut crate::simplec::Flow {
        if fl.cht.is_some() {
            return &mut fl.cht.as_deref_mut().expect("conjugate case").fe;
        }
        fl
    }

    /// The first circuit solve, with the heat of the circuit file, and its
    /// Robin conditions on the walls.
    pub fn start(&mut self, fl: &mut crate::simplec::Flow) -> Result<()> {
        let heat = self.applied.clone();
        let water = self.circuit.solve(&heat).map_err(|e| {
            bad(
                "circuits.file",
                format!("the start solve of the circuit: {e}"),
            )
        })?;
        self.updates.push(Update {
            iteration: self.first,
            residual: f64::NAN,
            heat: self.coupled_heat(&heat),
            newton: self.circuit.iterations,
            failed: false,
        });
        self.water = water;
        self.apply(fl)
    }

    /// Σ of the heats of the coupled tubes (W).
    fn coupled_heat(&self, heat: &[Vec<f64>]) -> f64 {
        let v: Vec<f64> = heat
            .iter()
            .enumerate()
            .filter(|(i, _)| self.circuit.coupled[*i])
            .flat_map(|(_, r)| r.iter().copied())
            .collect();
        exact_sum(&v)
    }

    /// The conductance U of segment `i` of a tube for the faces of wall
    /// entry `wall`.
    fn conductance(&self, wall: usize, tube: usize, i: usize) -> f64 {
        let w = &self.sec.walls[wall];
        let ws = self.water[tube][i];
        let phi = self.area3[tube][i] / self.circuit.inner_area(tube);
        1.0 / (w.wall_resistance + phi * (w.fouling + 1.0 / ws.htc))
    }

    /// The Robin conditions of the current water side on the walls: U and
    /// T_water of the segments at their midpoints, linear in between along
    /// the tube and extended linearly over the half segments at the two
    /// ends, so that the wall sees a smooth water side and not steps at the
    /// segment ends (the heat of each face still goes to its own segment).
    fn apply(&mut self, fl: &mut crate::simplec::Flow) -> Result<()> {
        let mut set = Vec::with_capacity(self.faces.len());
        for k in 0..self.faces.len() {
            let f = self.faces[k];
            let n = self.water[f.tube].len();
            let z = f.pos - 0.5;
            let i0 = (z.max(0.0).floor() as usize).min(n.saturating_sub(2));
            let (i1, a) = if n > 1 {
                (i0 + 1, z - i0 as f64)
            } else {
                (i0, 0.0)
            };
            let (u0, u1) = (
                self.conductance(f.wall, f.tube, i0),
                self.conductance(f.wall, f.tube, i1),
            );
            let (t0, t1) = (self.water[f.tube][i0].t_bulk, self.water[f.tube][i1].t_bulk);
            let u = u0 + a * (u1 - u0);
            let t = t0 + a * (t1 - t0);
            self.robin[k] = (u, t);
            set.push((f.fb, u, t));
        }
        let e = Self::energy(fl);
        if let Some(s) = e.scalars.as_mut() {
            s.set_energy_convection(&set);
        }
        Ok(())
    }

    /// The heat (W) each tube segment receives from the 3D faces now, with
    /// the Robin conditions as applied, and the heat flux into the water
    /// per mapped face (W/m²).
    fn heats(&self, fl: &mut crate::simplec::Flow) -> Result<(Vec<Vec<f64>>, Vec<f64>)> {
        let e = Self::energy(fl);
        let tb = match e.scalars.as_ref() {
            Some(s) => s.energy_boundary_temperature(e.dev.as_mut())?,
            None => None,
        };
        let tb =
            tb.ok_or_else(|| anyhow::anyhow!("[circuits]: the case has no energy equation"))?;
        let mut parts: Vec<Vec<Vec<f64>>> = self
            .applied
            .iter()
            .map(|r| vec![Vec::new(); r.len()])
            .collect();
        let mut flux = Vec::with_capacity(self.faces.len());
        for (k, f) in self.faces.iter().enumerate() {
            let (u, tw) = self.robin[k];
            let q = u * (tb[f.fb] - tw);
            flux.push(q);
            parts[f.tube][f.seg].push(q * f.area);
        }
        let heat = parts
            .iter()
            .enumerate()
            .map(|(i, r)| {
                if self.circuit.coupled[i] {
                    r.iter().map(|v| exact_sum(v)).collect()
                } else {
                    self.applied[i].clone()
                }
            })
            .collect();
        Ok((heat, flux))
    }

    /// Σ|Q_3D − Q_applied| / Σ|Q_3D| over the coupled segments.
    fn residual_of(&self, heat: &[Vec<f64>]) -> f64 {
        let (mut num, mut den) = (Vec::new(), Vec::new());
        for (i, (a, b)) in heat.iter().zip(&self.applied).enumerate() {
            if !self.circuit.coupled[i] {
                continue;
            }
            for (x, y) in a.iter().zip(b) {
                num.push((x - y).abs());
                den.push(x.abs());
            }
        }
        let den = exact_sum(&den);
        if den > 0.0 {
            exact_sum(&num) / den
        } else {
            exact_sum(&num)
        }
    }

    /// A circuit solve with the heats of the walls now (relaxed), at outer
    /// iteration `it`; a failed solve keeps the last water side, unless
    /// `strict`.
    fn update(
        &mut self,
        fl: &mut crate::simplec::Flow,
        it: usize,
        heat: Vec<Vec<f64>>,
        strict: bool,
    ) -> Result<()> {
        let residual = self.residual_of(&heat);
        let a = self.sec.relaxation;
        let next: Vec<Vec<f64>> = heat
            .iter()
            .zip(&self.applied)
            .map(|(n, o)| n.iter().zip(o).map(|(n, o)| o + a * (n - o)).collect())
            .collect();
        match self.circuit.solve(&next) {
            Ok(water) => {
                self.updates.push(Update {
                    iteration: it,
                    residual,
                    heat: self.coupled_heat(&next),
                    newton: self.circuit.iterations,
                    failed: false,
                });
                self.water = water;
                self.applied = next;
                self.apply(fl)
            }
            Err(e) => {
                self.updates.push(Update {
                    iteration: it,
                    residual,
                    heat: self.coupled_heat(&next),
                    newton: 0,
                    failed: true,
                });
                if strict {
                    return Err(diag::err(Diagnostic::new(
                        diag::CIR_CONVERGENCE,
                        format!("the circuit solve at iteration {it}: {e}"),
                    )));
                }
                self.notes.push(format!(
                    "{}: the circuit solve at iteration {it} failed ({e}); the water side of the \
                     solve before is kept",
                    diag::CIR_NOT_COUPLED
                ));
                Ok(())
            }
        }
    }

    /// Before outer iteration `next`: a circuit solve every `every`
    /// iterations.
    pub fn before_iteration(&mut self, fl: &mut crate::simplec::Flow, next: usize) -> Result<()> {
        let done = next.saturating_sub(self.first + 1);
        if done > 0 && done.is_multiple_of(self.sec.every) {
            let (heat, _) = self.heats(fl)?;
            self.update(fl, next, heat, false)?;
        }
        Ok(())
    }

    /// The last solution of the circuit.
    pub fn solution(&self) -> Option<&foehncraft_circuits::output::Solution> {
        self.circuit.solution.as_ref()
    }

    /// Messages for the log (failed solves).
    pub fn take_notes(&mut self) -> Vec<String> {
        std::mem::take(&mut self.notes)
    }

    /// Whether the coupling is converged after outer iteration `it`, given
    /// that the flow is (`gas_done`); if the flow is converged and the
    /// coupling not, the circuit is solved at once with the new heats.
    pub fn settled(
        &mut self,
        fl: &mut crate::simplec::Flow,
        it: usize,
        gas_done: bool,
    ) -> Result<bool> {
        if !gas_done {
            return Ok(false);
        }
        let (heat, _) = self.heats(fl)?;
        self.residual = self.residual_of(&heat);
        if self.residual <= self.sec.tolerance {
            return Ok(true);
        }
        self.update(fl, it + 1, heat, false)?;
        Ok(false)
    }

    /// The summary of the end of the run, its files (`circuit.csv`,
    /// `circuit/result.json`, `circuit/profiles.csv`, `circuit/walls.csv`
    /// with every mapped face) when `write`, and the circuit's account.
    pub fn finish(
        &mut self,
        fl: &mut crate::simplec::Flow,
        mesh: &Mesh,
        geo: &Geometry,
        out: &Path,
        name: &str,
        write: bool,
    ) -> Result<(CircuitSummary, Account)> {
        let (heat, flux) = self.heats(fl)?;
        self.residual = self.residual_of(&heat);
        let sol = self
            .circuit
            .solution
            .clone()
            .ok_or_else(|| anyhow::anyhow!("[circuits]: the circuit was not solved"))?;
        let wall_heat = self.coupled_heat(&heat);
        let circuit_heat = self.coupled_heat(&self.applied);
        let m = &sol.model;
        let other: Vec<f64> = m
            .tubes
            .iter()
            .enumerate()
            .filter(|(i, _)| !self.circuit.coupled[*i])
            .map(|(_, t)| t.heat() * t.count)
            .collect();
        let other = exact_sum(&other);
        let (inflow, outflow) = sol.enthalpy_flows();
        let total = circuit_heat + other;
        let closure = (inflow + total - outflow) / total.abs().max(1e-300);
        let drums = m
            .nodes
            .iter()
            .filter_map(|n| {
                sol.drum(&n.name).map(|d| DrumReport {
                    name: n.name.clone(),
                    steam_flow: d.steam,
                    circulation_ratio: d.circulation_ratio,
                })
            })
            .collect();
        let tubes = m
            .tubes
            .iter()
            .zip(&sol.marches)
            .enumerate()
            .filter(|(i, _)| self.circuit.coupled[*i])
            .map(|(i, (t, mr))| TubeReport {
                name: t.name.clone(),
                heat: exact_sum(&self.applied[i]),
                flow: sol.x[t.im] * t.count,
                outlet_temperature: foehncraft_kernels::steam::state_ph(mr.p_end, mr.h_end)
                    .map_or(f64::NAN, |s| s.t),
                outlet_quality: mr.segments.last().map_or(f64::NAN, |s| s.x_eq),
                wall_temperature_max: mr
                    .segments
                    .iter()
                    .map(|s| s.t_wall)
                    .fold(f64::NAN, f64::max),
                chf_ratio_min: foehncraft_circuits::output::min_chf(mr).map_or(f64::NAN, |c| c.1),
            })
            .collect();
        let summary = CircuitSummary {
            file: self.sec.file.display().to_string(),
            solves: self.updates.len(),
            failed: self.updates.iter().filter(|u| u.failed).count(),
            residual: self.residual,
            converged: self.residual <= self.sec.tolerance,
            wall_heat,
            circuit_heat,
            circuit_closure: closure,
            drums,
            tubes,
            updates: self.updates.clone(),
        };
        if write {
            let dir = out.join("circuit");
            std::fs::create_dir_all(&dir)?;
            let json = sol.to_json(name, &self.circuit.case);
            std::fs::write(
                dir.join("result.json"),
                serde_json::to_string_pretty(&json)?,
            )?;
            std::fs::write(dir.join("profiles.csv"), sol.to_csv())?;
            std::fs::write(out.join("circuit.csv"), self.csv(name))?;
            std::fs::write(dir.join("walls.csv"), self.walls_csv(mesh, geo, &flux))?;
        }
        let mut a = Account::new(names::THERMAL, "J");
        a.region = "circuit".into();
        a.terms = vec![
            Term::volume("from the walls", TermKind::Exchange, circuit_heat).paired(PAIR_WALLS),
            Term::volume("other tubes", TermKind::Source, other),
            Term::volume("water in", TermKind::Source, inflow),
            Term::volume("water and steam out", TermKind::Source, -outflow),
        ];
        a.storage = Some(0.0);
        Ok((summary, a))
    }

    /// Marks the conduction of the coupled boundaries in the thermal
    /// account of the domain as the other side of [`PAIR_WALLS`].
    pub fn pair_domain(case: &FlowCase, a: &mut Account) {
        let b = boundaries(case);
        for t in &mut a.terms {
            if t.name == "conduction" && t.boundary.as_deref().is_some_and(|n| b.contains(&n)) {
                t.pair = Some(PAIR_WALLS);
            }
        }
    }

    /// `circuit.csv`: one row per circuit solve.
    fn csv(&self, name: &str) -> String {
        let mut s = format!(
            "# circuit: the circuit solves of the coupled walls (SF7c) of `{name}`\n# units: \
             iteration - (the outer iteration it ran before); residual - (Σ|Q_3D − Q_circuit| / \
             Σ|Q_3D| of the heats it took); heat_w W (heat of the coupled tubes it was given); \
             newton - (iterations); failed - (1: the water side before was kept)\n\
             iteration,residual,heat_w,newton,failed\n"
        );
        for u in &self.updates {
            let _ = writeln!(
                s,
                "{},{:e},{:e},{},{}",
                u.iteration,
                u.residual,
                u.heat,
                u.newton,
                u8::from(u.failed)
            );
        }
        s
    }

    /// `circuit/walls.csv`: every mapped face with its segment, the Robin
    /// condition and the heat flux into the water.
    fn walls_csv(&self, mesh: &Mesh, geo: &Geometry, flux: &[f64]) -> String {
        let ni = mesh.n_interior_faces();
        let mut s = String::from(
            "boundary,tube,segment,x_m,y_m,z_m,area_m2,U_W_m2K,T_water_K,q_W_m2,T_wall_K\n",
        );
        for (k, f) in self.faces.iter().enumerate() {
            let c = geo.face_centroid[ni + f.fb];
            let (u, tw) = self.robin[k];
            let _ = writeln!(
                s,
                "{},{},{},{:.6},{:.6},{:.6},{:.6e},{:.6e},{:.6},{:.6e},{:.6}",
                self.sec.walls[f.wall].boundary,
                self.circuit.model.tubes[f.tube].name,
                f.seg,
                c[0],
                c[1],
                c[2],
                f.area,
                u,
                tw,
                flux[k],
                tw + flux[k] / u
            );
        }
        s
    }
}

/// The log frame of the coupled circuit.
pub fn frame(s: &CircuitSummary, tolerance: f64) -> crate::console::Frame {
    use crate::console::sci3;
    let mut f = crate::console::Frame::new("Water and steam circuit (SF7c)").summary();
    f.push_kv(
        "Coupling",
        format!(
            "{} circuit solves ({} failed), residual {} (tolerance {}){}",
            s.solves,
            s.failed,
            sci3(s.residual),
            sci3(tolerance),
            if s.converged { "" } else { ": not converged" }
        ),
    );
    f.push_kv(
        "Heat",
        format!(
            "walls {} W, circuit {} W, circuit closure {}",
            sci3(s.wall_heat),
            sci3(s.circuit_heat),
            sci3(s.circuit_closure)
        ),
    );
    for d in &s.drums {
        f.push_kv(
            &format!("Drum {}", d.name),
            format!(
                "steam {} kg/s, circulation ratio {:.2}",
                sci3(d.steam_flow),
                d.circulation_ratio
            ),
        );
    }
    for t in &s.tubes {
        f.push_kv(
            &format!("Tube {}", t.name),
            format!(
                "heat {} W, flow {} kg/s, outlet {:.2} K (x_eq {:.4}), wall up to {:.2} K, CHF \
                 ratio ≥ {}",
                sci3(t.heat),
                sci3(t.flow),
                t.outlet_temperature,
                t.outlet_quality,
                t.wall_temperature_max,
                if t.chf_ratio_min.is_finite() {
                    format!("{:.2}", t.chf_ratio_min)
                } else {
                    "–".into()
                }
            ),
        );
    }
    f
}
