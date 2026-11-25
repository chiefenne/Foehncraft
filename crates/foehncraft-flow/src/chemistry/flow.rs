//! Chemistry in the cells of a flow (steam furnace stage SF2b): laminar
//! reacting flows with the mechanisms of [`super::mechanism`].
//!
//! `[chemistry]` of a transient flow case with `[species]`, `[energy]` and
//! the variable density. The species of the case are the species of the
//! mechanism. Each time step is split (Lie splitting, chemistry first): at
//! the first outer iteration of the step every cell is integrated as a
//! constant-pressure adiabatic reactor over the step Δt from its state of
//! the old level, (Y, T) → (Y*, T*), with the stiff integrator of
//! [`foehncraft_kernels::chemistry`] (the same code as the 0D reactor and
//! the GPU). The change is the source of the step:
//!
//! - species i: S_i = (Y*_i − Y_i)/Δt (1/s, the scalar source, times ρ V
//!   in the rows);
//! - energy: q = −ρ Σ_i h_i(T) (Y*_i − Y_i)/Δt (W/m³), h_i the absolute
//!   enthalpy (formation included) per kg, so that the heat release is the
//!   enthalpy the reactor turned into sensible heat.
//!
//! The transport of the step then solves with these sources held, which is
//! first order in Δt in the coupling (Lie–Trotter splitting; the symmetric
//! splitting of Strang 1968 is second order, not done). The sources conserve every element exactly as the
//! reactor does (the rows of the species and the closure carry them), and
//! the ledger books them as the term `reaction` of each species account
//! (zero in each element account). Cells below `min_temperature` are not
//! integrated (no source).

use super::mechanism::{Mechanism, bad};
use super::reactor::IntegratorSection;
use crate::diag;
use crate::transient::OuterHook;
use anyhow::{Context, Result};
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_kernels::chemistry as kc;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// `[chemistry]` of a flow case (SF2b).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ChemistrySection {
    /// The mechanism: a built-in one (`gri30`, `westbrook-dryer-2`,
    /// `jones-lindstedt-4`) or a YAML mechanism file (D103), relative to
    /// the case file. Its species must be the species of `[species]`.
    pub mechanism: String,
    /// Cells colder than this are not integrated (K); 0 integrates all.
    #[serde(default)]
    pub min_temperature: f64,
    /// The stiff integrator in every cell.
    #[serde(default)]
    pub integrator: IntegratorSection,
    /// The speed of a premixed flame from the fuel consumption, for the
    /// report and its check.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flame_speed: Option<FlameSpeedSection>,
}

/// `[chemistry.flame_speed]`: the burning speed of a premixed flame from
/// the fuel the reactions consume, S_c = −ṁ_F/(ρ_u Y_F,u A), with ṁ_F the
/// net production of the fuel (kg/s), ρ_u and Y_F,u the density and fuel
/// mass fraction of the fresh mixture at the inlet and A the inlet area;
/// averaged over the end of the run. For a flat flame it is the laminar
/// burning velocity, whatever the flame does in the frame of the mesh.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FlameSpeedSection {
    /// The fuel species.
    pub fuel: String,
    /// The inlet of the fresh mixture.
    pub inlet: String,
    /// The average is taken over the part of the run after this fraction
    /// of its time.
    #[serde(default = "default_average_from")]
    pub average_from: f64,
    /// The band the speed must lie in (m/s): the run fails outside it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub range: Option<[f64; 2]>,
}

fn default_average_from() -> f64 {
    0.5
}

/// The chemistry of a run for the report.
#[derive(Debug, Clone, Serialize)]
pub struct Summary {
    /// The mechanism.
    pub mechanism: String,
    /// The heat release of the last step (W).
    pub heat_release: f64,
    /// Integrator steps of all cells over the run, accepted and rejected.
    pub steps: [u64; 2],
    /// The flame speed S_c (m/s), with `[chemistry.flame_speed]`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flame_speed: Option<f64>,
    /// Its band.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flame_speed_range: Option<[f64; 2]>,
    /// The time window of the average (s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub flame_speed_window: Option<[f64; 2]>,
}

impl Summary {
    /// The failure of the flame-speed band, if any.
    pub fn failure(&self) -> Option<String> {
        let (s, [lo, hi]) = (self.flame_speed?, self.flame_speed_range?);
        (!(s >= lo && s <= hi)).then(|| {
            format!("flame speed {s:.4} m/s outside [chemistry.flame_speed] range [{lo}, {hi}]")
        })
    }
}

impl ChemistrySection {
    /// Resolves a mechanism file against the folder of the case file.
    pub fn resolve_paths(&mut self, dir: &std::path::Path) {
        let m = &self.mechanism;
        let file = m.ends_with(".yaml") || m.ends_with(".yml") || m.contains('/');
        if file && std::path::Path::new(m).is_relative() {
            self.mechanism = dir.join(m).to_string_lossy().into_owned();
        }
    }
}

/// What a reacting flow does not have yet (CHM003).
fn unsupported(key: &str, what: &str) -> anyhow::Error {
    diag::err(
        diag::Diagnostic::new(diag::CHEMISTRY_UNSUPPORTED, format!("[chemistry] {what}"))
            .key(key.to_string()),
    )
}

/// The checks of `[chemistry]` (CHM001, CHM003).
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let Some(ch) = &case.chemistry else {
        return Ok(());
    };
    let Some(sp) = &case.species else {
        return Err(bad(
            "chemistry",
            "[chemistry] needs the species of its mechanism in [species]",
        ));
    };
    if case.energy.is_none() {
        return Err(bad(
            "chemistry",
            "[chemistry] needs the energy equation ([energy]) for the temperature and the heat \
             release",
        ));
    }
    if !crate::density::wanted(case) {
        return Err(bad(
            "chemistry",
            "[chemistry] needs the variable density ([energy] density = \"variable\"): the \
             density follows the heat release",
        ));
    }
    if case.time.is_none() {
        return Err(unsupported(
            "chemistry",
            "runs in transient cases only (a time step splits transport and chemistry); \
             steady reacting flows come with the turbulence–chemistry interaction (SF3)",
        ));
    }
    if case.turbulence.model != crate::turbulence::TurbulenceModel::Laminar {
        return Err(unsupported(
            "chemistry",
            "with a turbulence model: the mean reaction rates of turbulent flames come with \
             SF3 (eddy dissipation concept)",
        ));
    }
    if !case.spray.is_empty() {
        return Err(unsupported(
            "chemistry",
            "with sprays: evaporating and reacting sprays come later (SP5, SF15)",
        ));
    }
    if !(ch.min_temperature >= 0.0 && ch.min_temperature.is_finite()) {
        return Err(bad("chemistry.min_temperature", "must be a number ≥ 0 K"));
    }
    ch.integrator.check()?;
    if let Some(f) = &ch.flame_speed {
        if !case.boundary.contains_key(&f.inlet) || !case.boundary_species.contains_key(&f.inlet) {
            return Err(bad(
                "chemistry.flame_speed.inlet",
                format!(
                    "`{}` is not an inlet with [boundary.{}.species]",
                    f.inlet, f.inlet
                ),
            ));
        }
        if !(f.average_from >= 0.0 && f.average_from < 1.0) {
            return Err(bad(
                "chemistry.flame_speed.average_from",
                "must lie in [0, 1)",
            ));
        }
        if let Some([lo, hi]) = f.range
            && !(lo < hi)
        {
            return Err(bad(
                "chemistry.flame_speed.range",
                "needs [low, high] with low < high",
            ));
        }
    }
    let m = Mechanism::load(&ch.mechanism, None)?;
    let names = crate::species::Resolved::of(sp)?.names;
    let missing: Vec<&str> = m
        .species
        .iter()
        .filter(|s| !names.contains(s))
        .map(String::as_str)
        .collect();
    let extra: Vec<&str> = names
        .iter()
        .filter(|s| !m.species.contains(s))
        .map(String::as_str)
        .collect();
    if let Some(f) = &ch.flame_speed
        && !names.contains(&f.fuel)
    {
        return Err(bad(
            "chemistry.flame_speed.fuel",
            format!("`{}` is not a species of [species]", f.fuel),
        ));
    }
    if !missing.is_empty() || !extra.is_empty() {
        return Err(bad(
            "species.names",
            format!(
                "[species] names must be the species of the mechanism `{}`{}{}",
                m.name,
                if missing.is_empty() {
                    String::new()
                } else {
                    format!("; missing: {}", missing.join(", "))
                },
                if extra.is_empty() {
                    String::new()
                } else {
                    format!("; not in the mechanism: {}", extra.join(", "))
                }
            ),
        ));
    }
    Ok(())
}

/// The reaction state of a run: the mechanism in the order of the
/// species of the case, and the rows of the reaction log.
pub struct Chemistry {
    mech: Mechanism,
    /// For each species of the mechanism its index among the species of
    /// the case (closure last).
    to_flow: Vec<usize>,
    settings: kc::Settings,
    min_temperature: f64,
    pressure: f64,
    volume: Vec<f64>,
    /// The reaction of the last step (a rejected step is overwritten by
    /// its repetition).
    pub last: Option<Row>,
    /// Per step: (end time, Δt, net production of the flame-speed fuel).
    fuel: Vec<(f64, f64, f64)>,
    /// The index of the flame-speed fuel among the species of the case.
    fuel_index: Option<usize>,
    /// Integrator steps of the run.
    steps: [u64; 2],
    /// The GPU that integrates the cells (`[run] backend` other than
    /// `cpu`); none: the CPU threads.
    gpu: Option<Box<dyn super::device::Device>>,
}

/// One row of the reaction log.
#[derive(Debug, Clone, Serialize)]
pub struct Row {
    /// The time at the end of the step (s).
    pub time: f64,
    /// The heat release Σ q V (W).
    pub heat_release: f64,
    /// The net production Σ ρ V S per species of the case (kg/s).
    pub production: Vec<f64>,
    /// Accepted and rejected integrator steps over all cells.
    pub steps: [u64; 2],
    /// Cells integrated.
    pub cells: usize,
}

impl Chemistry {
    /// The chemistry of `case`, if it has `[chemistry]`; `volume` are the
    /// cell volumes.
    pub fn of(case: &crate::case::FlowCase, volume: &[f64]) -> Result<Option<Self>> {
        let Some(ch) = &case.chemistry else {
            return Ok(None);
        };
        validate(case)?;
        let mech = Mechanism::load(&ch.mechanism, None)?;
        let gpu = if case.run.backend == "cpu" {
            None
        } else {
            Some(super::device::create_gpu(
                &case.run.backend,
                &mech.program(),
                0,
                &case.run.ptx,
            )?)
        };
        let sp = case.species.as_ref().context("species")?;
        let names = crate::species::Resolved::of(sp)?.names;
        let to_flow = mech
            .species
            .iter()
            .map(|s| names.iter().position(|n| n == s).context("species"))
            .collect::<Result<Vec<_>>>()?;
        Ok(Some(Self {
            mech,
            to_flow,
            settings: ch.integrator.settings(),
            min_temperature: ch.min_temperature,
            pressure: case.fluid.conditions().pressure,
            volume: volume.to_vec(),
            last: None,
            fuel: Vec::new(),
            fuel_index: ch
                .flame_speed
                .as_ref()
                .and_then(|f| names.iter().position(|n| *n == f.fuel)),
            steps: [0, 0],
            gpu,
        }))
    }

    /// The mechanism.
    pub fn mechanism(&self) -> &Mechanism {
        &self.mech
    }

    /// Integrates every cell over `dt` from the current fields and sets
    /// the sources of the step; `time` is the time at the end of the step.
    pub fn react(&mut self, fl: &mut crate::simplec::Flow, dt: f64, time: f64) -> Result<()> {
        let p = fl.density.as_ref().map_or(self.pressure, |d| d.p0);
        let sc = fl.scalars.as_ref().context("scalars")?;
        let link = sc.species().context("species")?;
        let t = sc.temperature(fl.dev.as_mut())?.context("energy")?;
        let n = self.mech.len();
        let nflow = link.res.names.len();
        let nc = self.volume.len();
        let y = &link.y_cell;
        let h_flow: Vec<&foehncraft_kernels::species::Species> = link.res.data.iter().collect();
        let set = self.settings;
        let tmin = self.min_temperature;
        let map = &self.to_flow;
        // the states (Y ≥ 0, T) of all cells, cell-major
        let nz = n + 1;
        let mut state = vec![0.0; nc * nz];
        for (c, z) in state.chunks_mut(nz).enumerate() {
            for (k, &j) in map.iter().enumerate() {
                z[k] = y[j][c].max(0.0);
            }
            z[n] = t[c];
        }
        let mut stats = vec![kc::Stats::default(); nc];
        match self.gpu.as_mut() {
            Some(g) => g.react(p, dt, tmin, &set, &mut state, &mut stats)?,
            None => super::device::react_cpu(
                &self.mech.program(),
                p,
                dt,
                tmin,
                &set,
                &mut state,
                &mut stats,
            ),
        }
        // per cell: (ΔY per flow species, heat release W/m³ per unit ρ,
        // steps, rejected, status); a cell that did not react has no ΔY
        let out: Vec<(Vec<f64>, f64, u32, u32, u32)> = (0..nc)
            .into_par_iter()
            .map(|c| {
                let tc = t[c];
                if !(tc >= tmin) || !(tc > 0.0) {
                    return (Vec::new(), 0.0, 0, 0, 0);
                }
                let z = &state[c * nz..(c + 1) * nz];
                let mut dy = vec![0.0; nflow];
                for (k, &j) in map.iter().enumerate() {
                    dy[j] = z[k] - y[j][c].max(0.0);
                }
                // −Σ h_i(T) ΔY_i, h_i = (h/RT)_i R T / W_i (J/kg)
                let dh: f64 = (0..nflow)
                    .map(|j| {
                        let sp = h_flow[j];
                        -sp.h_over_rt(tc) * foehncraft_kernels::species::R * tc / sp.molar_mass
                            * dy[j]
                    })
                    .exact_sum();
                let st = stats[c];
                (dy, dh, st.steps, st.rejected, st.status)
            })
            .collect();
        if let Some(c) = out.iter().position(|o| o.4 != 0) {
            return Err(diag::err(
                diag::Diagnostic::new(
                    diag::CHEMISTRY_INTEGRATOR,
                    format!(
                        "the chemistry of cell {c} (T = {:.1} K) did not reach the end of the \
                         time step {dt:.3e} s ({})",
                        t[c],
                        if out[c].4 == 1 {
                            "too many steps"
                        } else {
                            "the step became too small"
                        }
                    ),
                )
                .key("chemistry.integrator"),
            ));
        }
        let rho: Vec<f64> = (0..nc).map(|c| fl.phys.rho_at(c)).collect();
        let transported = nflow - 1;
        let mut src = vec![vec![0.0; nc]; transported];
        let mut heat = vec![0.0; nc];
        let (mut steps, mut rejected, mut cells) = (0u64, 0u64, 0usize);
        for (c, (dy, dh, s, r, _)) in out.iter().enumerate() {
            if dy.is_empty() {
                continue;
            }
            cells += 1;
            steps += u64::from(*s);
            rejected += u64::from(*r);
            for (j, row) in src.iter_mut().enumerate() {
                row[c] = dy[j] / dt;
            }
            heat[c] = rho[c] * dh / dt;
        }
        let production: Vec<f64> = (0..nflow)
            .map(|j| {
                (0..nc)
                    .map(|c| match out[c].0.get(j) {
                        Some(d) => rho[c] * self.volume[c] * d / dt,
                        None => 0.0,
                    })
                    .exact_sum()
            })
            .collect();
        let heat_release = (0..nc).map(|c| heat[c] * self.volume[c]).exact_sum();
        // a step taken back and repeated replaces the rows from its start
        let start = time - dt;
        while self.fuel.last().is_some_and(|r| r.0 > start + 1e-9 * dt) {
            self.fuel.pop();
        }
        if let Some(i) = self.fuel_index {
            self.fuel.push((time, dt, production[i]));
        }
        self.steps[0] += steps;
        self.steps[1] += rejected;
        self.last = Some(Row {
            time,
            heat_release,
            production,
            steps: [steps, rejected],
            cells,
        });
        let sc = fl.scalars.as_mut().context("scalars")?;
        sc.set_reaction(fl.dev.as_mut(), src, heat)
    }

    /// The report of the run; `boundaries` are the boundary summaries of
    /// the end (the inlet area).
    pub fn summary(
        &self,
        case: &crate::case::FlowCase,
        boundaries: &[crate::case::BoundarySummary],
    ) -> Result<Summary> {
        let fs = case.chemistry.as_ref().and_then(|c| c.flame_speed.as_ref());
        let (mut speed, mut window) = (None, None);
        if let (Some(f), Some(sp), Some((t_end, _, _))) =
            (fs, case.species.as_ref(), self.fuel.last())
        {
            let res = crate::species::Resolved::of(sp)?;
            let comp = case.boundary_species.get(&f.inlet).context("inlet")?;
            // the transported species, then the closure
            let mut y: Vec<f64> = res
                .mass_fractions(comp, &format!("boundary.{}.species", f.inlet))?
                .iter()
                .map(|v| match v {
                    crate::expr::Scalar::Number(x) => Ok(*x),
                    _ => Err(bad(
                        "chemistry.flame_speed.inlet",
                        "the fresh mixture of the inlet must be numbers",
                    )),
                })
                .collect::<Result<_>>()?;
            y.push(1.0 - y.iter().exact_sum());
            let t_u = case
                .boundary_energy
                .get(&f.inlet)
                .and_then(|e| match &e.temperature {
                    Some(crate::expr::Scalar::Number(t)) => Some(*t),
                    _ => None,
                })
                .context("the inlet needs a fixed temperature")?;
            let w: f64 = y
                .iter()
                .zip(&res.data)
                .map(|(y, s)| y / s.molar_mass)
                .exact_sum();
            let rho_u = self.pressure / (foehncraft_kernels::species::R * t_u * w);
            let area = boundaries
                .iter()
                .find(|b| b.name == f.inlet)
                .map(|b| b.area)
                .context("inlet area")?;
            let k = res
                .names
                .iter()
                .position(|n| *n == f.fuel)
                .context("fuel")?;
            let t0 = self.fuel.first().map_or(0.0, |r| r.0 - r.1);
            let from = t0 + f.average_from * (t_end - t0);
            let rows: Vec<&(f64, f64, f64)> = self
                .fuel
                .iter()
                .filter(|r| r.0 - r.1 >= from - 1e-12)
                .collect();
            let dt: f64 = rows.iter().map(|r| r.1).exact_sum();
            let m: f64 = rows.iter().map(|r| r.1 * r.2).exact_sum();
            if dt > 0.0 && y[k] > 0.0 && area > 0.0 {
                speed = Some(-m / dt / (rho_u * y[k] * area));
                window = Some([from, *t_end]);
            }
        }
        Ok(Summary {
            mechanism: self.mech.name.clone(),
            heat_release: self.last.as_ref().map_or(0.0, |r| r.heat_release),
            steps: self.steps,
            flame_speed: speed,
            flame_speed_range: fs.and_then(|f| f.range),
            flame_speed_window: window,
        })
    }

    /// The units line of `chemistry.csv`.
    pub const CSV_UNITS: &'static str = "step -; time s; heat_release W; production_<species> kg/s \
         (net production by the reactions, ρ V dY/dt summed over the cells); steps, rejected \
         (integrator steps over all cells); cells (cells integrated)";

    /// The header of `chemistry.csv`.
    pub fn csv_head(names: &[String]) -> String {
        let mut s = String::from("step,time,heat_release");
        for n in names {
            s.push_str(&format!(",production_{n}"));
        }
        s.push_str(",steps,rejected,cells");
        s
    }

    /// The row of `chemistry.csv` of step `step`, if a step reacted.
    pub fn csv_row(&self, step: usize) -> Option<String> {
        let r = self.last.as_ref()?;
        let mut s = format!("{step},{:.9e},{:.9e}", r.time, r.heat_release);
        for v in &r.production {
            s.push_str(&format!(",{v:.9e}"));
        }
        s.push_str(&format!(",{},{},{}", r.steps[0], r.steps[1], r.cells));
        Some(s)
    }
}

/// The reaction step of one time step as an outer-iteration hook: the
/// cells are integrated at the first outer iteration.
pub struct Step<'a> {
    chem: &'a mut Chemistry,
    dt: f64,
    time: f64,
}

impl<'a> Step<'a> {
    /// The step of size `dt` that ends at `time`.
    pub fn new(chem: &'a mut Chemistry, dt: f64, time: f64) -> Self {
        Self { chem, dt, time }
    }
}

impl OuterHook for Step<'_> {
    fn max_outer(&self) -> usize {
        0
    }

    fn before_momentum(&mut self, fl: &mut crate::simplec::Flow, k: usize) -> Result<()> {
        if k == 1 {
            self.chem.react(fl, self.dt, self.time)?;
        }
        Ok(())
    }

    fn after_solve(&mut self, _fl: &mut crate::simplec::Flow, _k: usize) -> Result<bool> {
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use crate::case::FlowCase;
    use crate::diag;
    use std::path::Path;

    const CASE: &str = include_str!("../../../../examples/cases/flow/flame_methane_air_jl4.toml");

    fn code(text: &str) -> Option<String> {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow");
        FlowCase::from_toml(text, &dir)
            .err()
            .map(|e| diag::of_error(&e).code.to_string())
    }

    /// The flame case reads; steady and turbulent runs, a constant density,
    /// species other than the mechanism's and a broken flame-speed table
    /// are refused with their codes.
    #[test]
    fn chemistry_cases_are_checked() {
        assert_eq!(code(CASE), None);
        // the case without its [time] table
        let (head, rest) = CASE.split_once("[time]\n").unwrap();
        let steady = format!("{head}{}", &rest[rest.find("\n\n").unwrap()..]);
        assert_eq!(code(&steady).as_deref(), Some("CHM003"));
        let turbulent = format!("{CASE}\n[turbulence]\nmodel = \"k-omega-sst\"\n");
        assert_eq!(code(&turbulent).as_deref(), Some("CHM003"));
        let constant = CASE.replace("density = \"variable\"", "density = \"constant\"");
        assert_ne!(constant, CASE);
        assert!(code(&constant).is_some());
        let fewer = CASE.replace("\"H2\", ", "");
        assert_ne!(fewer, CASE);
        assert_eq!(code(&fewer).as_deref(), Some("CHM001"));
        let fuel = CASE.replace("fuel = \"CH4\"", "fuel = \"C3H8\"");
        assert_ne!(fuel, CASE);
        assert_eq!(code(&fuel).as_deref(), Some("CHM001"));
        let range = CASE.replace("range = [0.38, 0.48]", "range = [0.48, 0.38]");
        assert_ne!(range, CASE);
        assert_eq!(code(&range).as_deref(), Some("CHM001"));
        let inlet = CASE.replace("inlet = \"xmin\"", "inlet = \"xmax\"");
        assert_ne!(inlet, CASE);
        assert_eq!(code(&inlet).as_deref(), Some("CHM001"));
    }
}
