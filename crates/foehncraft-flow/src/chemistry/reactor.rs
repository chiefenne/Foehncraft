//! The 0D constant-pressure adiabatic reactor (SF2): the check of a
//! mechanism and of the integrator of the cells, against the equilibrium
//! of [`super::equilibrium`] and against measured ignition delays.
//!
//! Case file:
//!
//! ```toml
//! [reactor]
//! mechanism = "gri30"            # or "westbrook-dryer-2", "jones-lindstedt-4", a YAML file
//! pressure = 101325.0            # Pa
//! temperature = 1500.0           # K, the start
//! end_time = 0.05                # s
//! equivalence_ratio = 1.0        # with fuel and oxidiser (mole fractions) ...
//! fuel = { CH4 = 1.0 }
//! oxidiser = { O2 = 0.21, N2 = 0.79 }
//! # mole_fractions = { CH4 = 0.095, O2 = 0.19, N2 = 0.715 }   # ... or the mixture
//! samples = 2000                 # output times, evenly spaced
//!
//! [reactor.integrator]           # the settings of the cells as well (SF2b)
//! rtol = 1e-6
//! atol = 1e-12
//! ```

use super::equilibrium;
use super::mechanism::{Mechanism, bad};
use anyhow::Result;
use foehncraft_kernels::chemistry as kc;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

fn default_samples() -> usize {
    2000
}

fn default_pressure() -> f64 {
    101_325.0
}

/// `[reactor.integrator]`: the settings of the stiff integrator.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct IntegratorSection {
    /// Relative tolerance.
    #[serde(default = "rtol")]
    pub rtol: f64,
    /// Absolute tolerance of the mass fractions.
    #[serde(default = "atol")]
    pub atol: f64,
    /// Absolute tolerance of the temperature (K).
    #[serde(default = "atol_t")]
    pub atol_temperature: f64,
    /// The largest number of steps of one interval.
    #[serde(default = "max_steps")]
    pub max_steps: u32,
}

fn rtol() -> f64 {
    kc::Settings::default().rtol
}
fn atol() -> f64 {
    kc::Settings::default().atol
}
fn atol_t() -> f64 {
    kc::Settings::default().atol_t
}
fn max_steps() -> u32 {
    kc::Settings::default().max_steps
}

impl Default for IntegratorSection {
    fn default() -> Self {
        Self {
            rtol: rtol(),
            atol: atol(),
            atol_temperature: atol_t(),
            max_steps: max_steps(),
        }
    }
}

impl IntegratorSection {
    /// The kernel settings.
    pub fn settings(&self) -> kc::Settings {
        kc::Settings {
            rtol: self.rtol,
            atol: self.atol,
            atol_t: self.atol_temperature,
            h0: 0.0,
            max_steps: self.max_steps,
        }
    }

    pub(crate) fn check(&self) -> Result<()> {
        let finite = self.atol.is_finite() && self.atol_temperature.is_finite();
        if !(finite
            && self.rtol > 0.0
            && self.rtol < 1.0
            && self.atol > 0.0
            && self.atol_temperature > 0.0)
        {
            return Err(bad(
                "reactor.integrator",
                "rtol must lie in (0, 1) and atol, atol_temperature be > 0",
            ));
        }
        if self.max_steps == 0 {
            return Err(bad("reactor.integrator.max_steps", "max_steps must be > 0"));
        }
        Ok(())
    }
}

/// The spacing of the output times of a reactor.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Spacing {
    /// Evenly in t.
    #[default]
    Linear,
    /// Evenly in log t, from `end_time`/10⁶.
    Log,
}

/// `[reactor]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReactorSection {
    /// A built-in mechanism (`gri30`, `westbrook-dryer-2`,
    /// `jones-lindstedt-4`) or a YAML mechanism file (D103).
    pub mechanism: String,
    /// Pressure (Pa).
    #[serde(default = "default_pressure")]
    pub pressure: f64,
    /// Start temperature (K).
    pub temperature: f64,
    /// The time to integrate (s).
    pub end_time: f64,
    /// The mixture as mole fractions (a species left out has 0).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mole_fractions: Option<BTreeMap<String, f64>>,
    /// The mixture as mass fractions.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mass_fractions: Option<BTreeMap<String, f64>>,
    /// The fuel/oxidiser ratio over the stoichiometric one, with `fuel` and
    /// `oxidiser`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub equivalence_ratio: Option<f64>,
    /// The fuel as mole fractions.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fuel: Option<BTreeMap<String, f64>>,
    /// The oxidiser as mole fractions.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub oxidiser: Option<BTreeMap<String, f64>>,
    /// The number of output times after the start.
    #[serde(default = "default_samples")]
    pub samples: usize,
    /// How the output times are spaced: evenly over `end_time`
    /// (`linear`), or evenly in log t from `end_time`/10⁶ (`log`).
    #[serde(default)]
    pub spacing: Spacing,
    /// The integrator.
    #[serde(default)]
    pub integrator: IntegratorSection,
    /// Bands the result must meet (R5: a case checks itself; a band not
    /// met is a failure and a non-zero exit code).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub check: Option<CheckSection>,
}

/// `[reactor.check]`: the bands of a reactor case.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckSection {
    /// The largest |T_end − T_eq| (K).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature_difference: Option<f64>,
    /// The largest |x_end − x_eq| of any species.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mole_fraction_difference: Option<f64>,
    /// The range of the equilibrium temperature (K).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub equilibrium_temperature: Option<[f64; 2]>,
    /// The range of the ignition delay (s).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ignition_delay: Option<[f64; 2]>,
    /// The largest relative drift of an element mass.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub element_drift: Option<f64>,
}

impl CheckSection {
    /// The bands not met by `r`, one line each.
    pub fn failures(&self, r: &ReactorResult) -> Vec<String> {
        let mut out = Vec::new();
        if let Some(b) = self.temperature_difference
            && !(r.temperature_difference.abs() <= b)
        {
            out.push(format!(
                "|T_end − T_eq| = {:.3e} K above {b} K",
                r.temperature_difference.abs()
            ));
        }
        if let Some(b) = self.mole_fraction_difference
            && !(r.largest_mole_fraction_difference <= b)
        {
            out.push(format!(
                "|x_end − x_eq| = {:.3e} above {b}",
                r.largest_mole_fraction_difference
            ));
        }
        if let Some([lo, hi]) = self.equilibrium_temperature
            && !(r.equilibrium.temperature >= lo && r.equilibrium.temperature <= hi)
        {
            out.push(format!(
                "T_eq = {:.3} K outside [{lo}, {hi}] K",
                r.equilibrium.temperature
            ));
        }
        if let Some([lo, hi]) = self.ignition_delay {
            match r.ignition_delay {
                Some(t) if t >= lo && t <= hi => {}
                Some(t) => out.push(format!("ignition delay {t:.4e} s outside [{lo}, {hi}] s")),
                None => out.push("no ignition".into()),
            }
        }
        if let Some(b) = self.element_drift
            && !(r.element_drift <= b)
        {
            out.push(format!("element drift {:.3e} above {b}", r.element_drift));
        }
        out
    }
}

/// A reactor case file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReactorCase {
    /// The reactor.
    pub reactor: ReactorSection,
}

/// The mole fractions of a mixture, in the mechanism's order.
pub fn mixture(m: &Mechanism, r: &ReactorSection) -> Result<Vec<f64>> {
    let given = [
        r.mole_fractions.is_some(),
        r.mass_fractions.is_some(),
        r.equivalence_ratio.is_some(),
    ]
    .iter()
    .filter(|v| **v)
    .count();
    if given != 1 {
        return Err(bad(
            "reactor",
            "give one of `mole_fractions`, `mass_fractions` and `equivalence_ratio` (with `fuel` \
             and `oxidiser`)",
        ));
    }
    let vec_of = |map: &BTreeMap<String, f64>, key: &str| -> Result<Vec<f64>> {
        let mut v = vec![0.0; m.len()];
        for (s, x) in map {
            let k = m.index(s).ok_or_else(|| {
                bad(
                    key,
                    format!("`{s}` is not a species of the mechanism {}", m.name),
                )
            })?;
            if !(*x >= 0.0 && x.is_finite()) {
                return Err(bad(key, format!("the fraction of `{s}` is {x}")));
            }
            v[k] = *x;
        }
        let s: f64 = v.iter().sum();
        if !(s > 0.0) {
            return Err(bad(key, "the fractions sum to zero"));
        }
        Ok(v.iter().map(|x| x / s).collect())
    };
    if let Some(x) = &r.mole_fractions {
        return vec_of(x, "reactor.mole_fractions");
    }
    if let Some(y) = &r.mass_fractions {
        let y = vec_of(y, "reactor.mass_fractions")?;
        let inv: Vec<f64> = y
            .iter()
            .zip(&m.data)
            .map(|(y, s)| y / s.molar_mass)
            .collect();
        let t: f64 = inv.iter().sum();
        return Ok(inv.iter().map(|v| v / t).collect());
    }
    let phi = r.equivalence_ratio.unwrap_or(1.0);
    let (Some(fu), Some(ox)) = (&r.fuel, &r.oxidiser) else {
        return Err(bad(
            "reactor.equivalence_ratio",
            "`equivalence_ratio` needs `fuel` and `oxidiser`",
        ));
    };
    if !(phi > 0.0 && phi.is_finite()) {
        return Err(bad("reactor.equivalence_ratio", "must be > 0"));
    }
    let xf = vec_of(fu, "reactor.fuel")?;
    let xo = vec_of(ox, "reactor.oxidiser")?;
    // the oxygen atoms a mole needs for CO2 and H2O, less what it has
    let atoms = |e: &str, x: &[f64]| -> f64 {
        m.elements
            .iter()
            .position(|n| n == e)
            .map_or(0.0, |i| m.atoms[i].iter().zip(x).map(|(a, x)| a * x).sum())
    };
    let need = |x: &[f64]| 2.0 * atoms("C", x) + 0.5 * atoms("H", x) - atoms("O", x);
    let (nf, no) = (need(&xf), need(&xo));
    if !(nf > 0.0 && no < 0.0) {
        return Err(bad(
            "reactor.fuel",
            "the fuel must need oxygen and the oxidiser supply it",
        ));
    }
    // moles of oxidiser per mole of fuel at stoichiometry: nf/(−no)
    let a = nf / (-no) / phi;
    let tot = 1.0 + a;
    Ok(xf.iter().zip(&xo).map(|(f, o)| (f + a * o) / tot).collect())
}

/// The result of a reactor run.
#[derive(Debug, Clone, PartialEq, Serialize, schemars::JsonSchema)]
pub struct ReactorResult {
    /// The mechanism.
    pub mechanism: String,
    /// Its species.
    pub species: Vec<String>,
    /// Number of reactions.
    pub reactions: usize,
    /// The start mole fractions.
    pub initial_mole_fractions: Vec<f64>,
    /// Output times (s).
    pub time: Vec<f64>,
    /// Temperature at each time (K).
    pub temperature: Vec<f64>,
    /// Mole fractions at each time (`[time][species]`).
    pub mole_fractions: Vec<Vec<f64>>,
    /// The ignition delay (s): the time of the largest dT/dt between two
    /// samples, if the temperature rose by more than 400 K.
    pub ignition_delay: Option<f64>,
    /// The equilibrium at the start enthalpy and pressure.
    pub equilibrium: equilibrium::State,
    /// T at the end minus the equilibrium temperature (K).
    pub temperature_difference: f64,
    /// The largest |x_end − x_eq| over the species.
    pub largest_mole_fraction_difference: f64,
    /// The largest relative change of an element mass over the run.
    pub element_drift: f64,
    /// Accepted and rejected integrator steps.
    pub steps: u32,
    /// Rejected steps.
    pub rejected: u32,
    /// Wall time of the integration (s).
    pub seconds: f64,
}

/// Mass fractions of mole fractions.
pub fn mass_of_mole(m: &Mechanism, x: &[f64]) -> Vec<f64> {
    let w: f64 = x.iter().zip(&m.data).map(|(x, s)| x * s.molar_mass).sum();
    x.iter()
        .zip(&m.data)
        .map(|(x, s)| x * s.molar_mass / w)
        .collect()
}

/// Mole fractions of mass fractions.
pub fn mole_of_mass(m: &Mechanism, y: &[f64]) -> Vec<f64> {
    let inv: Vec<f64> = y
        .iter()
        .zip(&m.data)
        .map(|(y, s)| y / s.molar_mass)
        .collect();
    let t: f64 = inv.iter().sum();
    inv.iter().map(|v| v / t).collect()
}

/// Reads a reactor case file (TOML with `[reactor]`); returns the case and
/// the folder that relative mechanism files are read from. Every error
/// carries `CHM001`.
pub fn read(path: &std::path::Path) -> Result<(ReactorCase, Option<std::path::PathBuf>)> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| bad("reactor", format!("{}: {e}", path.display())))?;
    let c: ReactorCase = toml::from_str(&text)
        .map_err(|e| bad("reactor", format!("{}: {}", path.display(), e.message())))?;
    Ok((c, path.parent().map(std::path::Path::to_path_buf)))
}

/// Checks a reactor case before it runs: the state, the samples, the
/// integrator settings, the mechanism and the start mixture. Returns the
/// mechanism and the start mole fractions.
pub fn prepare(
    case: &ReactorCase,
    base: Option<&std::path::Path>,
) -> Result<(Mechanism, Vec<f64>)> {
    let r = &case.reactor;
    let finite = r.pressure.is_finite() && r.temperature.is_finite() && r.end_time.is_finite();
    if !(finite && r.pressure > 0.0 && r.temperature > 0.0 && r.end_time > 0.0) {
        return Err(bad(
            "reactor",
            "pressure, temperature and end_time must be finite and > 0",
        ));
    }
    if r.samples == 0 || r.samples > 1_000_000 {
        return Err(bad("reactor.samples", "samples must lie in 1 … 1 000 000"));
    }
    r.integrator.check()?;
    let m = Mechanism::load(&r.mechanism, base)?;
    let x0 = mixture(&m, r)?;
    Ok((m, x0))
}

/// Runs the reactor of `case`; `base` is the folder of relative mechanism
/// files.
pub fn run(case: &ReactorCase, base: Option<&std::path::Path>) -> Result<ReactorResult> {
    let r = &case.reactor;
    let (m, x0) = prepare(case, base)?;
    let y0 = mass_of_mole(&m, &x0);
    let n = m.len();
    let prog = m.program();
    let set = r.integrator.settings();
    let mut z: Vec<f64> = y0.clone();
    z.push(r.temperature);
    let mut work = vec![0.0; kc::work_len(n)];
    let at = |i: usize| match r.spacing {
        Spacing::Linear => r.end_time * i as f64 / r.samples as f64,
        Spacing::Log => {
            r.end_time * foehncraft_core::math::pow(1e-6, 1.0 - i as f64 / r.samples as f64)
        }
    };
    let mut time = vec![0.0];
    let mut temp = vec![r.temperature];
    let mut xs = vec![x0.clone()];
    let (mut steps, mut rejected) = (0, 0);
    let mut h = 0.0;
    let started = std::time::Instant::now();
    for i in 1..=r.samples {
        let (t0, t1) = (if i == 1 { 0.0 } else { at(i - 1) }, at(i));
        let dt = t1 - t0;
        let s = kc::Settings { h0: h, ..set };
        let st = kc::integrate(&prog, r.pressure, &mut z, dt, &s, &mut work);
        steps += st.steps;
        rejected += st.rejected;
        if st.status != 0 {
            return Err(crate::diag::err(
                crate::diag::Diagnostic::new(
                    crate::diag::CHEMISTRY_INTEGRATOR,
                    format!(
                        "the integrator stopped at t = {:.6e} s (status {}): {}",
                        t0,
                        st.status,
                        if st.status == 1 {
                            "too many steps"
                        } else {
                            "the step became too small"
                        }
                    ),
                )
                .key("reactor.integrator"),
            ));
        }
        h = st.h_last;
        time.push(t1);
        temp.push(z[n]);
        xs.push(mole_of_mass(&m, &z[..n]));
    }
    let seconds = started.elapsed().as_secs_f64();
    let eq = equilibrium::hp(&m, &y0, r.temperature, r.pressure)?;
    let last = xs.last().cloned().unwrap_or_default();
    let tend = *temp.last().unwrap_or(&r.temperature);
    let dx = last
        .iter()
        .zip(&eq.mole_fractions)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let b0 = equilibrium::elements_of(&m, &y0);
    let b1 = equilibrium::elements_of(&m, &z[..n]);
    let drift = b0
        .iter()
        .zip(&b1)
        .filter(|(a, _)| **a > 0.0)
        .map(|(a, b)| ((b - a) / a).abs())
        .fold(0.0, f64::max);
    let mut ignition = None;
    if tend - r.temperature > 400.0 {
        let mut best = 0.0;
        for i in 1..temp.len() {
            let d = (temp[i] - temp[i - 1]) / (time[i] - time[i - 1]);
            if d > best {
                best = d;
                ignition = Some(0.5 * (time[i] + time[i - 1]));
            }
        }
    }
    Ok(ReactorResult {
        mechanism: m.name.clone(),
        species: m.species.clone(),
        reactions: m.reactions.len(),
        initial_mole_fractions: x0,
        time,
        temperature: temp,
        mole_fractions: xs,
        ignition_delay: ignition,
        temperature_difference: tend - eq.temperature,
        largest_mole_fraction_difference: dx,
        equilibrium: eq,
        element_drift: drift,
        steps,
        rejected,
        seconds,
    })
}

/// The CSV text of a result: time, temperature and the mole fractions.
pub fn csv(r: &ReactorResult) -> String {
    let mut s = String::from("time,temperature");
    for n in &r.species {
        s.push_str(&format!(",X_{n}"));
    }
    s.push('\n');
    for i in 0..r.time.len() {
        s.push_str(&format!("{:.9e},{:.9e}", r.time[i], r.temperature[i]));
        for x in &r.mole_fractions[i] {
            s.push_str(&format!(",{x:.9e}"));
        }
        s.push('\n');
    }
    s
}
