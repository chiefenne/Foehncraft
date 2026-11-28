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
//!
//! # Turbulent flames (SF3, D349)
//!
//! With a turbulence model `interaction` is `"edc"` or `"edm"` (see
//! [`super::edc`]) and the mean rate of a cell is R_i = ρ κ (Y*_i − Y_i)
//! (EDC, the fine structures integrated over their residence time τ* from
//! the mean state) or the rates of the eddy dissipation model. Steady and
//! transient cases:
//!
//! - **transient**: at the first outer iteration of each step the rates of
//!   the old level are held over the step, as for a laminar flame. The EDC
//!   increment is the exact solution of dY/dt = κ (Y* − Y) with Y* held,
//!   ΔY = (1 − e^(−κ Δt)) (Y* − Y): bounded by Y* for any Δt, one factor for
//!   all species of a cell, so every element is conserved exactly. The EDM
//!   rates of a cell are scaled down together when a species would become
//!   negative within the step (also exact in the elements);
//! - **steady**: every `update_every` outer iterations the sources are
//!   renewed in value-plus-slope form, S_i = S_u,i + S_p,i Y_i with
//!   S_p ≤ 0 on the diagonal of the species rows (EDC: S_u = κ Y*,
//!   S_p = −κ; EDM: a net sink as S_p = S/Y at the iterate), relaxed with
//!   `relaxation` between renewals; the heat release is held from the
//!   renewal. The rows conserve the elements at the converged state (where
//!   Y is the state of the renewal), and the ledger books S_u + S_p Y.
//!
//! Cells where the fine-structure reactor fails (too many steps, a step too
//! small) get no source in that renewal and are counted (CHM004), as are the
//! cells where γ* is held at its guard.

use super::edc::{EdcSection, EdmSection, Interaction};
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
    /// The turbulence–chemistry interaction (SF3): `"laminar"` (the rates
    /// of the mean state; laminar flows only), `"edc"` (eddy dissipation
    /// concept) or `"edm"` (eddy dissipation model, global mechanisms).
    #[serde(default)]
    pub interaction: Interaction,
    /// The constants of the eddy dissipation concept.
    #[serde(default)]
    pub edc: EdcSection,
    /// The constants of the eddy dissipation model.
    #[serde(default)]
    pub edm: EdmSection,
    /// Steady cases: the outer iterations between renewals of the
    /// reaction sources.
    #[serde(default = "one")]
    pub update_every: usize,
    /// Steady cases: the relaxation of the reaction sources between
    /// renewals, S ← S + α (S_new − S).
    #[serde(default = "unit")]
    pub relaxation: f64,
}

fn one() -> usize {
    1
}
fn unit() -> f64 {
    1.0
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
    /// The turbulence–chemistry interaction.
    pub interaction: &'static str,
    /// The fine structures of the last renewal (EDC).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fine_structures: Option<FineStats>,
    /// Cells whose fine-structure reactor failed, summed over the renewals.
    pub failed_cells: u64,
}

/// The fine structures of a renewal (EDC).
#[derive(Debug, Clone, Copy, Default, Serialize, PartialEq)]
pub struct FineStats {
    /// Cells with fine structures (turbulent and hot enough).
    pub cells: usize,
    /// Cells with γ* held at `max_fine_fraction`.
    pub clipped: usize,
    /// The smallest and the largest residence time τ* (s).
    pub tau: [f64; 2],
    /// The largest fine-structure fraction γ*.
    pub fraction_max: f64,
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
fn unsupported(key: &str, what: impl std::fmt::Display) -> anyhow::Error {
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
    let laminar = case.turbulence.model == crate::turbulence::TurbulenceModel::Laminar;
    match (ch.interaction, laminar) {
        (Interaction::Laminar, false) => {
            return Err(unsupported(
                "chemistry.interaction",
                "with a turbulence model needs a turbulence–chemistry interaction: \
                 interaction = \"edc\" (eddy dissipation concept) or \"edm\" (eddy \
                 dissipation model)",
            ));
        }
        (Interaction::Edc | Interaction::Edm, true) => {
            return Err(unsupported(
                "chemistry.interaction",
                format!(
                    "\"{}\" needs a turbulence model ([turbulence] model): laminar flames take \
                     interaction = \"laminar\"",
                    ch.interaction.name()
                ),
            ));
        }
        _ => {}
    }
    if case.time.is_none() && laminar {
        return Err(unsupported(
            "chemistry",
            "in a steady laminar case: a laminar flame runs in a transient case (a time step \
             splits transport and chemistry); steady reacting flows take a turbulence model \
             and interaction = \"edc\" or \"edm\"",
        ));
    }
    ch.edc.check()?;
    if ch.update_every == 0 {
        return Err(bad("chemistry.update_every", "must be ≥ 1"));
    }
    if !(ch.relaxation > 0.0 && ch.relaxation <= 1.0) {
        return Err(bad("chemistry.relaxation", "must lie in (0, 1]"));
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
    if ch.interaction == Interaction::Edm {
        ch.edm.check(&m.program(), &m.name)?;
    }
    if ch.flame_speed.is_some() && case.time.is_none() {
        return Err(bad(
            "chemistry.flame_speed",
            "the flame speed is an average over the time of a transient run",
        ));
    }
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
    /// The turbulence–chemistry interaction (SF3).
    interaction: Interaction,
    edc: EdcSection,
    edm: EdmSection,
    update_every: usize,
    relaxation: f64,
    /// Steady: the sources held between renewals (S_u, S_p per
    /// transported species, the heat release).
    held: Option<Held>,
    /// The reaction of the last step (a rejected step is overwritten by
    /// its repetition).
    pub last: Option<Row>,
    /// Per step: (end time, Δt, net production of the flame-speed fuel).
    fuel: Vec<(f64, f64, f64)>,
    /// The index of the flame-speed fuel among the species of the case.
    fuel_index: Option<usize>,
    /// Integrator steps of the run.
    steps: [u64; 2],
    /// Cells whose fine-structure reactor failed, over the run.
    failed: u64,
    /// The GPU that integrates the cells (`[run] backend` other than
    /// `cpu`); none: the CPU threads.
    gpu: Option<Box<dyn super::device::Device>>,
}

/// The sources of a steady run held between renewals.
#[derive(Debug, Clone)]
struct Held {
    su: Vec<Vec<f64>>,
    sp: Vec<Vec<f64>>,
    heat: Vec<f64>,
}

/// One row of the reaction log.
#[derive(Debug, Clone, Serialize)]
pub struct Row {
    /// The time at the end of the step (s); 0 in a steady run.
    pub time: f64,
    /// The heat release Σ q V (W).
    pub heat_release: f64,
    /// The net production Σ ρ V S per species of the case (kg/s).
    pub production: Vec<f64>,
    /// Accepted and rejected integrator steps over all cells.
    pub steps: [u64; 2],
    /// Cells integrated.
    pub cells: usize,
    /// Cells whose fine-structure reactor failed (EDC).
    pub failed: usize,
    /// The fine structures (EDC).
    pub fine: Option<FineStats>,
}

/// The result of the reactions of one cell.
#[derive(Debug, Clone, Default)]
struct CellOut {
    /// S_u per species of the case (1/s); empty: no reaction.
    su: Vec<f64>,
    /// S_p (≤ 0) per species of the case (1/s), steady runs.
    sp: Vec<f64>,
    /// The heat release per unit density (W/kg).
    heat: f64,
    steps: u32,
    rejected: u32,
    failed: bool,
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
        // the eddy dissipation model runs on the CPU threads (no stiff
        // integration)
        let gpu = if case.run.backend == "cpu" || ch.interaction == Interaction::Edm {
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
            interaction: ch.interaction,
            edc: ch.edc.clone(),
            edm: ch.edm.clone(),
            update_every: ch.update_every,
            relaxation: ch.relaxation,
            held: None,
            last: None,
            fuel: Vec::new(),
            fuel_index: ch
                .flame_speed
                .as_ref()
                .and_then(|f| names.iter().position(|n| *n == f.fuel)),
            steps: [0, 0],
            failed: 0,
            gpu,
        }))
    }

    /// The mechanism.
    pub fn mechanism(&self) -> &Mechanism {
        &self.mech
    }

    /// The interaction.
    pub fn interaction(&self) -> Interaction {
        self.interaction
    }

    /// Steady runs: renews the reaction sources before outer iteration
    /// `it` (every `update_every` iterations after `first`, the iteration
    /// the run starts from); true when it renewed them.
    pub fn before_iteration(
        &mut self,
        fl: &mut crate::simplec::Flow,
        it: usize,
        first: usize,
    ) -> Result<bool> {
        if self.held.is_some() && (it.saturating_sub(first + 1)) % self.update_every != 0 {
            return Ok(false);
        }
        self.react_turbulent(fl, None, 0.0)?;
        Ok(true)
    }

    /// Integrates every cell over `dt` from the current fields and sets
    /// the sources of the step; `time` is the time at the end of the step.
    pub fn react(&mut self, fl: &mut crate::simplec::Flow, dt: f64, time: f64) -> Result<()> {
        if self.interaction != Interaction::Laminar {
            return self.react_turbulent(fl, Some(dt), time);
        }
        self.react_laminar(fl, dt, time)
    }

    /// The finding of the cells whose fine-structure reactor failed.
    pub fn finding(&self) -> Option<diag::Diagnostic> {
        (self.failed > 0).then(|| {
            diag::Diagnostic::new(
                diag::CHEMISTRY_FINE_STRUCTURES,
                format!(
                    "the fine-structure reactor failed {} times over the run (cells × renewals); \
                     those cells got no reaction source in that renewal",
                    self.failed
                ),
            )
            .key("chemistry.integrator")
        })
    }

    /// The reactions of a turbulent flame (SF3) from the current fields:
    /// over the step `dt` of a transient run, or the renewal of the held
    /// sources of a steady one (`dt` none, `time` the outer iteration).
    fn react_turbulent(
        &mut self,
        fl: &mut crate::simplec::Flow,
        dt: Option<f64>,
        time: f64,
    ) -> Result<()> {
        let p = fl.density.as_ref().map_or(self.pressure, |d| d.p0);
        let sc = fl.scalars.as_ref().context("scalars")?;
        let link = sc.species().context("species")?;
        let t = sc.temperature(fl.dev.as_mut())?.context("energy")?;
        let (k, eps) = fl
            .turbulence
            .as_ref()
            .context("turbulence")?
            .k_and_epsilon(fl.dev.as_mut())?;
        let n = self.mech.len();
        let nflow = link.res.names.len();
        let transported = nflow - 1;
        let nc = self.volume.len();
        let y = &link.y_cell;
        let h_flow: Vec<&foehncraft_kernels::species::Species> = link.res.data.iter().collect();
        let rho: Vec<f64> = (0..nc).map(|c| fl.phys.rho_at(c)).collect();
        // the viscosity of a cell: the mean of its faces
        let mu: Vec<f64> = match fl.phys.mu_face.as_deref() {
            Some(mf) => {
                let d = &fl.disc;
                let (mut sum, mut cnt) = (vec![0.0; nc], vec![0u32; nc]);
                for (f, &m) in mf.iter().enumerate().take(d.n_faces) {
                    let o = d.owner[f] as usize;
                    if o < nc {
                        sum[o] += m;
                        cnt[o] += 1;
                    }
                    if f < d.n_interior {
                        let nb = d.neighbour[f] as usize;
                        if nb < nc {
                            sum[nb] += m;
                            cnt[nb] += 1;
                        }
                    }
                }
                sum.iter()
                    .zip(&cnt)
                    .map(|(s, &c)| s / f64::from(c.max(1)))
                    .collect()
            }
            None => vec![fl.phys.mu; nc],
        };
        let tmin = self.min_temperature;
        let map = &self.to_flow;
        let nz = n + 1;
        let mut state = vec![0.0; nc * nz];
        for (c, z) in state.chunks_mut(nz).enumerate() {
            for (kk, &j) in map.iter().enumerate() {
                z[kk] = y[j][c].max(0.0);
            }
            z[n] = t[c];
        }
        let hot = |c: usize| t[c] >= tmin && t[c] > 0.0;
        // h_i(T) per kg of species j of the case at the temperature of c
        let h_of = |j: usize, c: usize| {
            let sp = h_flow[j];
            sp.h_over_rt(t[c]) * foehncraft_kernels::species::R * t[c] / sp.molar_mass
        };
        let mut fine = None;
        let mut fields: Vec<(String, Vec<f64>)> = Vec::new();
        let out: Vec<CellOut> = match self.interaction {
            Interaction::Edc => {
                let fs: Vec<Option<super::edc::FineStructure>> = (0..nc)
                    .map(|c| {
                        hot(c)
                            .then(|| {
                                super::edc::fine_structure(&self.edc, mu[c] / rho[c], k[c], eps[c])
                            })
                            .flatten()
                    })
                    .collect();
                let taus: Vec<f64> = fs.iter().map(|f| f.map_or(0.0, |f| f.tau)).collect();
                let mut stats = vec![kc::Stats::default(); nc];
                let set = self.settings;
                match self.gpu.as_mut() {
                    Some(g) => g.react(p, 0.0, Some(&taus), tmin, &set, &mut state, &mut stats)?,
                    None => super::device::react_cpu(
                        &self.mech.program(),
                        p,
                        0.0,
                        Some(&taus),
                        tmin,
                        &set,
                        &mut state,
                        &mut stats,
                    ),
                }
                let mut st = FineStats {
                    tau: [f64::INFINITY, 0.0],
                    ..FineStats::default()
                };
                for f in fs.iter().flatten() {
                    st.cells += 1;
                    st.clipped += usize::from(f.clipped);
                    st.tau = [st.tau[0].min(f.tau), st.tau[1].max(f.tau)];
                    st.fraction_max = st.fraction_max.max(f.fraction);
                }
                if st.cells == 0 {
                    st.tau = [0.0, 0.0];
                }
                fine = Some(st);
                fields.push((
                    "fine_structure_fraction".into(),
                    fs.iter().map(|f| f.map_or(0.0, |f| f.fraction)).collect(),
                ));
                fields.push(("fine_structure_time".into(), taus.clone()));
                (0..nc)
                    .into_par_iter()
                    .map(|c| {
                        let Some(f) = fs[c] else {
                            return CellOut::default();
                        };
                        let s = stats[c];
                        if s.status != 0 {
                            return CellOut {
                                steps: s.steps,
                                rejected: s.rejected,
                                failed: true,
                                ..CellOut::default()
                            };
                        }
                        let z = &state[c * nz..(c + 1) * nz];
                        // Y* per species of the case, and Y* − Y
                        let mut ys = vec![0.0; nflow];
                        let mut dy = vec![0.0; nflow];
                        for (kk, &j) in map.iter().enumerate() {
                            ys[j] = z[kk];
                            dy[j] = z[kk] - y[j][c].max(0.0);
                        }
                        let (su, sp, heat) = match dt {
                            // the exact increment of dY/dt = κ (Y* − Y) over
                            // the step with Y* held
                            Some(dt) => {
                                let g = -(-f.kappa * dt).exp_m1();
                                let inc: Vec<f64> = dy.iter().map(|d| g * d).collect();
                                let dh: f64 = (0..nflow).map(|j| -h_of(j, c) * inc[j]).exact_sum();
                                (inc.iter().map(|d| d / dt).collect(), Vec::new(), dh / dt)
                            }
                            None => {
                                let dh: f64 = (0..nflow).map(|j| -h_of(j, c) * dy[j]).exact_sum();
                                (
                                    ys.iter().map(|v| f.kappa * v).collect(),
                                    vec![-f.kappa; nflow],
                                    f.kappa * dh,
                                )
                            }
                        };
                        CellOut {
                            su,
                            sp,
                            heat,
                            steps: s.steps,
                            rejected: s.rejected,
                            failed: false,
                        }
                    })
                    .collect()
            }
            Interaction::Edm => {
                let prog = self.mech.program();
                let edm = &self.edm;
                let wk: Vec<f64> = self
                    .mech
                    .program()
                    .species
                    .iter()
                    .map(|s| s.molar_mass)
                    .collect();
                (0..nc)
                    .into_par_iter()
                    .map_init(
                        || {
                            (
                                vec![0.0; n],
                                vec![0.0; super::edc::edm_work_len(n)],
                                vec![0.0; n],
                            )
                        },
                        |(w, work, conc), c| {
                            if !hot(c) || !(k[c] > 0.0 && eps[c] > 0.0) {
                                return CellOut::default();
                            }
                            let z = &state[c * nz..(c + 1) * nz];
                            super::edc::edm_rates(
                                &prog,
                                edm,
                                p,
                                t[c],
                                &z[..n],
                                eps[c] / k[c],
                                w,
                                work,
                            );
                            let (rk, _) = kc::concentrations(prog.species, p, t[c], &z[..n], conc);
                            // dY/dt per species of the case
                            let mut s = vec![0.0; nflow];
                            for (kk, &j) in map.iter().enumerate() {
                                s[j] = wk[kk] * w[kk] / rk;
                            }
                            match dt {
                                Some(dt) => {
                                    // one factor for the cell: no species
                                    // below zero within the step
                                    let mut f: f64 = 1.0;
                                    for kk in 0..n {
                                        if w[kk] < 0.0 {
                                            f = f.min(conc[kk] / (-w[kk] * dt));
                                        }
                                    }
                                    let f = f.max(0.0);
                                    let su: Vec<f64> = s.iter().map(|v| f * v).collect();
                                    let dh: f64 =
                                        (0..nflow).map(|j| -h_of(j, c) * su[j]).exact_sum();
                                    CellOut {
                                        su,
                                        heat: dh,
                                        ..CellOut::default()
                                    }
                                }
                                None => {
                                    let dh: f64 =
                                        (0..nflow).map(|j| -h_of(j, c) * s[j]).exact_sum();
                                    let mut su = vec![0.0; nflow];
                                    let mut sp = vec![0.0; nflow];
                                    for j in 0..nflow {
                                        if s[j] >= 0.0 {
                                            su[j] = s[j];
                                        } else {
                                            sp[j] = s[j] / y[j][c].max(1e-12);
                                        }
                                    }
                                    CellOut {
                                        su,
                                        sp,
                                        heat: dh,
                                        ..CellOut::default()
                                    }
                                }
                            }
                        },
                    )
                    .collect()
            }
            Interaction::Laminar => unreachable!("laminar flames take react_laminar"),
        };
        // the sources of the transported species and the heat (W/m³)
        let mut su = vec![vec![0.0; nc]; transported];
        let mut sp = vec![vec![0.0; nc]; transported];
        let mut heat = vec![0.0; nc];
        let (mut steps, mut rejected, mut cells, mut failed) = (0u64, 0u64, 0usize, 0usize);
        for (c, o) in out.iter().enumerate() {
            steps += u64::from(o.steps);
            rejected += u64::from(o.rejected);
            if o.failed {
                failed += 1;
            }
            if o.su.is_empty() {
                continue;
            }
            cells += 1;
            for j in 0..transported {
                su[j][c] = o.su[j];
                if !o.sp.is_empty() {
                    sp[j][c] = o.sp[j];
                }
            }
            heat[c] = rho[c] * o.heat;
        }
        // steady: relaxed against the sources held
        if dt.is_none() {
            if let Some(h) = self.held.as_ref() {
                let a = self.relaxation;
                let mix = |old: &[f64], new: &mut [f64]| {
                    for (n, o) in new.iter_mut().zip(old) {
                        *n = o + a * (*n - o);
                    }
                };
                for j in 0..transported {
                    mix(&h.su[j], &mut su[j]);
                    mix(&h.sp[j], &mut sp[j]);
                }
                mix(&h.heat, &mut heat);
            }
            self.held = Some(Held {
                su: su.clone(),
                sp: sp.clone(),
                heat: heat.clone(),
            });
        }
        // the net production at the state of the renewal, S_u + S_p Y; the
        // closure species takes the rest
        let mut production: Vec<f64> = (0..transported)
            .map(|j| {
                (0..nc)
                    .map(|c| rho[c] * self.volume[c] * (su[j][c] + sp[j][c] * y[j][c]))
                    .exact_sum()
            })
            .collect();
        production.push(-production.iter().exact_sum());
        let heat_release = (0..nc).map(|c| heat[c] * self.volume[c]).exact_sum();
        if let Some(t) = dt {
            let start = time - t;
            while self.fuel.last().is_some_and(|r| r.0 > start + 1e-9 * t) {
                self.fuel.pop();
            }
            if let Some(i) = self.fuel_index {
                self.fuel.push((time, t, production[i]));
            }
        }
        self.steps[0] += steps;
        self.steps[1] += rejected;
        self.failed += failed as u64;
        self.last = Some(Row {
            time,
            heat_release,
            production,
            steps: [steps, rejected],
            cells,
            failed,
            fine,
        });
        fields.insert(0, ("heat_release_rate".into(), heat.clone()));
        fl.phys.reaction_fields = Some(std::sync::Arc::new(fields));
        let sc = fl.scalars.as_mut().context("scalars")?;
        match dt {
            Some(_) => sc.set_reaction(fl.dev.as_mut(), su, heat),
            None => sc.set_reaction_with(fl.dev.as_mut(), su, Some(sp), heat),
        }
    }

    /// A laminar flame (SF2b): every cell integrated over `dt` from the
    /// current fields, the sources of the step.
    fn react_laminar(&mut self, fl: &mut crate::simplec::Flow, dt: f64, time: f64) -> Result<()> {
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
            Some(g) => g.react(p, dt, None, tmin, &set, &mut state, &mut stats)?,
            None => super::device::react_cpu(
                &self.mech.program(),
                p,
                dt,
                None,
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
            failed: 0,
            fine: None,
        });
        fl.phys.reaction_fields = Some(std::sync::Arc::new(vec![(
            "heat_release_rate".into(),
            heat.clone(),
        )]));
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
            interaction: self.interaction.name(),
            fine_structures: self.last.as_ref().and_then(|r| r.fine),
            failed_cells: self.failed,
        })
    }

    /// The units line of `chemistry.csv`.
    pub const CSV_UNITS: &'static str = "step - (the outer iteration in a steady run); time s; \
         heat_release W; production_<species> kg/s (net production by the reactions, ρ V dY/dt \
         summed over the cells); steps, rejected (integrator steps over all cells); cells (cells \
         that reacted); failed (cells whose fine-structure reactor failed, EDC); clipped (cells \
         with γ* at max_fine_fraction, EDC); tau_min s, tau_max s (the range of the \
         fine-structure residence time τ*, EDC); fraction_max - (the largest γ*, EDC)";

    /// The header of `chemistry.csv`.
    pub fn csv_head(names: &[String]) -> String {
        let mut s = String::from("step,time,heat_release");
        for n in names {
            s.push_str(&format!(",production_{n}"));
        }
        s.push_str(",steps,rejected,cells,failed,clipped,tau_min,tau_max,fraction_max");
        s
    }

    /// The row of `chemistry.csv` of step `step`, if a step reacted.
    pub fn csv_row(&self, step: usize) -> Option<String> {
        let r = self.last.as_ref()?;
        let mut s = format!("{step},{:.9e},{:.9e}", r.time, r.heat_release);
        for v in &r.production {
            s.push_str(&format!(",{v:.9e}"));
        }
        let f = r.fine.unwrap_or_default();
        s.push_str(&format!(
            ",{},{},{},{},{},{:.6e},{:.6e},{:.6e}",
            r.steps[0],
            r.steps[1],
            r.cells,
            r.failed,
            f.clipped,
            f.tau[0],
            f.tau[1],
            f.fraction_max
        ));
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
        // the checkout may carry CRLF line ends (Windows)
        let case = CASE.replace("\r\n", "\n");
        let case = case.as_str();
        assert_eq!(code(case), None);
        // the case without its [time] table
        let (head, rest) = case.split_once("[time]\n").unwrap();
        let steady = format!("{head}{}", &rest[rest.find("\n\n").unwrap()..]);
        assert_eq!(code(&steady).as_deref(), Some("CHM003"));
        let turbulent = format!("{case}\n[turbulence]\nmodel = \"k-omega-sst\"\n");
        assert_eq!(code(&turbulent).as_deref(), Some("CHM003"));
        // SF3: an interaction needs a turbulence model, the eddy
        // dissipation model a mechanism in one direction, and the
        // constants are checked
        let with = |extra: &str| {
            case.replace(
                "min_temperature = 400.0",
                &format!("min_temperature = 400.0\n{extra}"),
            )
        };
        assert_eq!(
            code(&with("interaction = \"edc\"")).as_deref(),
            Some("CHM003")
        );
        let edc = with("interaction = \"edc\"");
        let edc_turb = format!("{edc}\n[turbulence]\nmodel = \"k-omega-sst\"\n");
        assert_eq!(code(&edc_turb), None);
        let edm = format!(
            "{}\n[turbulence]\nmodel = \"k-epsilon\"\n",
            with("interaction = \"edm\"")
        );
        assert_eq!(code(&edm).as_deref(), Some("CHM001"));
        let guard = format!(
            "{}\n[turbulence]\nmodel = \"k-epsilon\"\n",
            with("interaction = \"edc\"\nedc = { max_fine_fraction = 1.0 }")
        );
        assert_eq!(code(&guard).as_deref(), Some("CHM001"));
        let every = format!(
            "{}\n[turbulence]\nmodel = \"k-epsilon\"\n",
            with("interaction = \"edc\"\nupdate_every = 0")
        );
        assert_eq!(code(&every).as_deref(), Some("CHM001"));
        let constant = case.replace("density = \"variable\"", "density = \"constant\"");
        assert_ne!(constant, case);
        assert!(code(&constant).is_some());
        let fewer = case.replace("\"H2\", ", "");
        assert_ne!(fewer, case);
        assert_eq!(code(&fewer).as_deref(), Some("CHM001"));
        let fuel = case.replace("fuel = \"CH4\"", "fuel = \"C3H8\"");
        assert_ne!(fuel, case);
        assert_eq!(code(&fuel).as_deref(), Some("CHM001"));
        let range = case.replace("range = [0.38, 0.48]", "range = [0.48, 0.38]");
        assert_ne!(range, case);
        assert_eq!(code(&range).as_deref(), Some("CHM001"));
        let inlet = case.replace("inlet = \"xmin\"", "inlet = \"xmax\"");
        assert_ne!(inlet, case);
        assert_eq!(code(&inlet).as_deref(), Some("CHM001"));
    }
}
