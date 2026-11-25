//! `chemistry.*`: gas-phase chemistry (stage SF2): the 0D reactor and the
//! chemical equilibrium of a mechanism.

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::chemistry::{equilibrium, mechanism::Mechanism, reactor};
use foehncraft_flow::diag;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// The reactor of a call: a case file or the table inline, and the folder
/// relative mechanism files are read from.
fn section(
    case: &Option<PathBuf>,
    inline: &Option<reactor::ReactorSection>,
    ctx: &mut Ctx<'_>,
) -> Result<(reactor::ReactorSection, Option<PathBuf>)> {
    match (case, inline) {
        (Some(p), None) => {
            ctx.roots.check_read(p, "the reactor case file")?;
            let (c, base) = reactor::read(p)?;
            Ok((c.reactor, base))
        }
        (None, Some(r)) => Ok((r.clone(), None)),
        _ => Err(diag::err(
            diag::Diagnostic::new(diag::CHEMISTRY, "give one of `case` and `reactor`")
                .key("reactor"),
        )),
    }
}

/// Input of `chemistry.reactor`, the options of `foehn chemistry reactor`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReactorInput {
    /// A reactor case file (TOML with `[reactor]`).
    #[serde(default)]
    pub case: Option<PathBuf>,
    /// Or the `[reactor]` table inline.
    #[serde(default)]
    pub reactor: Option<reactor::ReactorSection>,
    /// The CSV file of the trace (time, temperature, mole fractions);
    /// default `output/<case name>/reactor.csv` (as the output directory
    /// of a flow case), none for an inline reactor.
    #[serde(default)]
    pub csv: Option<PathBuf>,
}

/// Output of `chemistry.reactor`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct ReactorOutput {
    /// `foehncraft.reactor`.
    pub format: &'static str,
    /// The mechanism.
    pub mechanism: String,
    /// Number of species and reactions.
    pub species: usize,
    /// Number of reactions.
    pub reactions: usize,
    /// The ignition delay (s), the time of the largest dT/dt.
    pub ignition_delay: Option<f64>,
    /// The start and end temperature (K).
    pub temperature: [f64; 2],
    /// The equilibrium temperature at the start enthalpy (K).
    pub equilibrium_temperature: f64,
    /// T_end − T_eq (K).
    pub temperature_difference: f64,
    /// The largest |x_end − x_eq|.
    pub largest_mole_fraction_difference: f64,
    /// The largest relative change of an element mass.
    pub element_drift: f64,
    /// Accepted and rejected integrator steps.
    pub steps: [u32; 2],
    /// Wall time of the integration (s).
    pub seconds: f64,
    /// The end mole fractions above 1e-6, by species.
    pub final_mole_fractions: Vec<(String, f64)>,
    /// The CSV file written, if one.
    pub csv: Option<PathBuf>,
    /// The bands of `[reactor.check]` not met.
    pub failures: Vec<String>,
}

fn csv_path(a: &ReactorInput) -> Option<PathBuf> {
    a.csv.clone().or_else(|| {
        a.case.as_ref().and_then(|c| {
            c.file_stem()
                .map(|s| Path::new("output").join(s).join("reactor.csv"))
        })
    })
}

/// `chemistry.reactor`.
pub struct Reactor;

impl Tool for Reactor {
    const ID: &'static str = "chemistry.reactor";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Integrate the 0D constant-pressure adiabatic reactor of a mechanism (SF2) and compare its end state with the chemical equilibrium.";
    const CATEGORY: Category = Category::Chemistry;
    const ACTION: Action = Action::WriteFiles;
    type Input = ReactorInput;
    type Output = ReactorOutput;

    fn effects(a: &ReactorInput, _ctx: &mut Ctx<'_>) -> Result<crate::Effects> {
        Ok(crate::Effects {
            writes: csv_path(a).into_iter().collect(),
            ..Default::default()
        })
    }

    fn call(a: ReactorInput, ctx: &mut Ctx<'_>) -> Result<ReactorOutput> {
        let (r, base) = section(&a.case, &a.reactor, ctx)?;
        let csv = csv_path(&a);
        if let Some(p) = &csv {
            ctx.roots.check_write(p)?;
        }
        let res = reactor::run(
            &reactor::ReactorCase { reactor: r.clone() },
            base.as_deref(),
        )?;
        let n = res.species.len();
        let last = res.mole_fractions.last().cloned().unwrap_or_default();
        let t1 = *res.temperature.last().unwrap_or(&r.temperature);
        ctx.line(&format!(
            "reactor: {} ({} species, {} reactions), p = {} Pa, T₀ = {} K, {} s",
            res.mechanism, n, res.reactions, r.pressure, r.temperature, r.end_time
        ));
        match res.ignition_delay {
            Some(t) => ctx.line(&format!("  ignition delay {t:.6e} s")),
            None => ctx.line("  no ignition (the temperature rose by less than 400 K)"),
        }
        ctx.line(&format!(
            "  T end {t1:.3} K, equilibrium {:.3} K, difference {:.3e} K",
            res.equilibrium.temperature, res.temperature_difference
        ));
        ctx.line(&format!(
            "  largest |x_end − x_eq| {:.3e}, element drift {:.3e}",
            res.largest_mole_fraction_difference, res.element_drift
        ));
        ctx.line(&format!(
            "  {} steps, {} rejected, {:.3} s",
            res.steps, res.rejected, res.seconds
        ));
        if let Some(p) = &csv {
            if let Some(d) = p.parent().filter(|d| !d.as_os_str().is_empty()) {
                std::fs::create_dir_all(d).map_err(|e| anyhow::anyhow!("{}: {e}", d.display()))?;
            }
            std::fs::write(p, reactor::csv(&res))
                .map_err(|e| anyhow::anyhow!("{}: {e}", p.display()))?;
            ctx.line(&format!("wrote {}", p.display()));
        }
        let mut fin: Vec<(String, f64)> = res
            .species
            .iter()
            .cloned()
            .zip(last)
            .filter(|(_, x)| *x > 1e-6)
            .collect();
        fin.sort_by(|a, b| b.1.total_cmp(&a.1));
        let failures = r
            .check
            .as_ref()
            .map(|c| c.failures(&res))
            .unwrap_or_default();
        for f in &failures {
            ctx.line(&format!("check failed: {f}"));
        }
        if r.check.is_some() && failures.is_empty() {
            ctx.line("checks passed");
        }
        Ok(ReactorOutput {
            format: "foehncraft.reactor",
            mechanism: res.mechanism,
            species: n,
            reactions: res.reactions,
            ignition_delay: res.ignition_delay,
            temperature: [r.temperature, t1],
            equilibrium_temperature: res.equilibrium.temperature,
            temperature_difference: res.temperature_difference,
            largest_mole_fraction_difference: res.largest_mole_fraction_difference,
            element_drift: res.element_drift,
            steps: [res.steps, res.rejected],
            seconds: res.seconds,
            final_mole_fractions: fin,
            csv,
            failures,
        })
    }

    fn exit_code(o: &ReactorOutput) -> i32 {
        if o.failures.is_empty() { 0 } else { 2 }
    }
}

/// What the equilibrium holds fixed.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Deserialize, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Constant {
    /// Enthalpy and pressure: the adiabatic flame state.
    #[default]
    Enthalpy,
    /// Temperature and pressure.
    Temperature,
}

/// Input of `chemistry.equilibrium`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EquilibriumInput {
    /// A reactor case file (TOML with `[reactor]`; `end_time` is not used).
    #[serde(default)]
    pub case: Option<PathBuf>,
    /// Or the `[reactor]` table inline.
    #[serde(default)]
    pub reactor: Option<reactor::ReactorSection>,
    /// What is held fixed with the pressure.
    #[serde(default)]
    pub constant: Constant,
}

/// Output of `chemistry.equilibrium`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct EquilibriumOutput {
    /// `foehncraft.equilibrium`.
    pub format: &'static str,
    /// The mechanism (its species are the species of the equilibrium).
    pub mechanism: String,
    /// The species.
    pub species: Vec<String>,
    /// The state.
    pub state: equilibrium::State,
}

/// `chemistry.equilibrium`.
pub struct Equilibrium;

impl Tool for Equilibrium {
    const ID: &'static str = "chemistry.equilibrium";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "The chemical equilibrium of a mixture by Gibbs minimisation (SF2), at fixed enthalpy and pressure (adiabatic flame) or temperature and pressure.";
    const CATEGORY: Category = Category::Chemistry;
    const ACTION: Action = Action::Read;
    type Input = EquilibriumInput;
    type Output = EquilibriumOutput;

    fn call(a: EquilibriumInput, ctx: &mut Ctx<'_>) -> Result<EquilibriumOutput> {
        let (r, base) = section(&a.case, &a.reactor, ctx)?;
        if !(r.pressure > 0.0 && r.temperature > 0.0) {
            return Err(diag::err(
                diag::Diagnostic::new(diag::CHEMISTRY, "pressure and temperature must be > 0")
                    .key("reactor"),
            ));
        }
        let m = Mechanism::load(&r.mechanism, base.as_deref())?;
        let x = reactor::mixture(&m, &r)?;
        let y = reactor::mass_of_mole(&m, &x);
        let s = match a.constant {
            Constant::Enthalpy => equilibrium::hp(&m, &y, r.temperature, r.pressure)?,
            Constant::Temperature => equilibrium::tp(&m, &y, r.temperature, r.pressure, None)?,
        };
        ctx.line(&format!(
            "equilibrium ({}) of {}: T = {:.3} K, p = {} Pa",
            match a.constant {
                Constant::Enthalpy => "h, p",
                Constant::Temperature => "T, p",
            },
            m.name,
            s.temperature,
            s.pressure
        ));
        let mut rows: Vec<(usize, f64)> = s.mole_fractions.iter().copied().enumerate().collect();
        rows.sort_by(|a, b| b.1.total_cmp(&a.1));
        for (k, v) in rows.iter().take_while(|(_, v)| *v > 1e-6) {
            ctx.line(&format!("  {:<8} {v:.6e}", m.species[*k]));
        }
        Ok(EquilibriumOutput {
            format: "foehncraft.equilibrium",
            mechanism: m.name.clone(),
            species: m.species.clone(),
            state: s,
        })
    }
}
