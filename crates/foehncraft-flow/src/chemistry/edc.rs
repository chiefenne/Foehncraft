//! Turbulence–chemistry interaction (steam furnace stage SF3, D349): the
//! mean reaction rates of turbulent flames from the eddy dissipation
//! concept (EDC) with finite-rate chemistry in the fine structures, and the
//! eddy dissipation model (EDM) for global mechanisms.
//!
//! # Eddy dissipation concept
//!
//! Magnussen (1981), "On the structure of turbulence and a generalized eddy
//! dissipation concept for chemical reaction in turbulent flow", 19th AIAA
//! Aerospace Sciences Meeting, AIAA 81-0042; Gran & Magnussen (1996), "A
//! numerical study of a bluff-body stabilized diffusion flame. Part 2.
//! Influence of combustion modeling and finite-rate chemistry", Combust.
//! Sci. Technol. 119, 191–217; Magnussen (2005), "The eddy dissipation
//! concept — a bridge between science and technology", ECCOMAS Thematic
//! Conference on Computational Combustion, Lisbon. The reactions take
//! place in the fine structures, the regions of the dissipative eddies.
//! From the cascade model of the turbulence, with the kinematic viscosity
//! ν, k and ε of the mean flow:
//!
//! - the length scale of the fine structures relative to the flow,
//!   γ_λ = C_γ (ν ε/k²)^(1/4), C_γ = 2.1377;
//! - their mass fraction γ* = γ_λ² (version 2005; γ_λ³ in version 1981);
//! - their residence time τ* = C_τ (ν/ε)^(1/2), C_τ = 0.4083.
//!
//! The fine structures are a constant-pressure adiabatic reactor that
//! starts from the mean state (Y, T) and reacts over τ* (the plug-flow
//! reactor of Gran & Magnussen 1996), with the stiff integrator of SF2:
//! (Y, T) → (Y*, T*). The exchange of mass between the fine structures and
//! their surroundings at the rate 1/τ* per unit mass of fine structures
//! gives the mean rate (per unit volume, χ = 1: every fine structure
//! reacts)
//!
//! R_i = ρ κ (Y*_i − Y_i),  κ = γ*/(τ* (1 − γ*)),
//!
//! from ρ γ* (Y*_i − Y°_i)/τ* with the surroundings Y° and the mean
//! Y = γ* Y* + (1 − γ*) Y°. The heat release is
//! q = −ρ κ Σ_i h_i(T) (Y*_i − Y_i), h_i the absolute enthalpy per kg.
//! γ* is held below `max_fine_fraction` (a numerical guard, not a part of
//! the model: κ grows without bound as γ* → 1 where the turbulence is weak,
//! ν ε/k² large, as near walls); the cells held are counted.
//!
//! # Eddy dissipation model
//!
//! Magnussen & Hjertager (1977), "On mathematical modeling of turbulent
//! combustion with special emphasis on soot formation and combustion",
//! 16th Symposium (International) on Combustion, 719–729: the reactions
//! are as fast as the turbulence mixes their reactants. Reaction j of a
//! global mechanism (forward only) has the mixing rate (mol/(m³ s))
//!
//! q_mix,j = A (ε/k) min( min_R C_R/ν'_R , B Σ_P ρ Y_P / Σ_P ν''_P W_P ),
//!
//! C_R = ρ Y_R/W_R the concentrations of its reactants and P its products,
//! A = 4, B = 0.5. With `kinetic_limit` the rate is the smaller of q_mix,j
//! and the rate of the mechanism q_kin,j (the kinetics hold the reactions
//! back where the gas is cold, as at the nozzle); without it, a reaction
//! runs wherever its reactants and products meet. A reaction and its exact
//! reverse (the CO oxidation and the CO₂ dissociation of Westbrook & Dryer)
//! are taken as a pair: the mixing rate caps their net kinetic rate in its
//! direction. Limiting each of the two alone lets a fast dissociation undo
//! a mixing-limited oxidation, and a jet flame then kept CO at the level of
//! CO₂ and released heat with the wrong sign where its hot pilot gas met
//! slow mixing.

use super::mechanism::bad;
use anyhow::Result;
use foehncraft_kernels::chemistry as kc;
use foehncraft_kernels::species::Species;
use serde::{Deserialize, Serialize};

/// `[chemistry] interaction`: how the mean reaction rate of a cell
/// follows from the chemistry.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Interaction {
    /// The rates of the mean state (laminar flames, SF2b).
    #[default]
    Laminar,
    /// The eddy dissipation concept with finite-rate chemistry (SF3).
    Edc,
    /// The eddy dissipation model of a global mechanism (SF3).
    Edm,
}

impl Interaction {
    /// The name in the case file.
    pub fn name(self) -> &'static str {
        match self {
            Self::Laminar => "laminar",
            Self::Edc => "edc",
            Self::Edm => "edm",
        }
    }
}

/// `[chemistry.edc] version`.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
pub enum EdcVersion {
    /// Magnussen (2005): γ* = γ_λ².
    #[default]
    #[serde(rename = "2005")]
    V2005,
    /// Magnussen (1981): γ* = γ_λ³.
    #[serde(rename = "1981")]
    V1981,
}

/// `[chemistry.edc]`: the constants of the eddy dissipation concept.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EdcSection {
    /// The version: `"2005"` (γ* = γ_λ², default) or `"1981"` (γ_λ³).
    #[serde(default)]
    pub version: EdcVersion,
    /// C_γ of the length scale γ_λ = C_γ (ν ε/k²)^(1/4).
    #[serde(default = "c_gamma")]
    pub c_gamma: f64,
    /// C_τ of the residence time τ* = C_τ (ν/ε)^(1/2).
    #[serde(default = "c_tau")]
    pub c_tau: f64,
    /// The largest mass fraction of fine structures γ* (a numerical guard).
    #[serde(default = "max_fine_fraction")]
    pub max_fine_fraction: f64,
}

fn c_gamma() -> f64 {
    2.1377
}
fn c_tau() -> f64 {
    0.4083
}
fn max_fine_fraction() -> f64 {
    0.75
}

impl Default for EdcSection {
    fn default() -> Self {
        Self {
            version: EdcVersion::default(),
            c_gamma: c_gamma(),
            c_tau: c_tau(),
            max_fine_fraction: max_fine_fraction(),
        }
    }
}

impl EdcSection {
    /// The checks of the constants (CHM001).
    pub fn check(&self) -> Result<()> {
        if !(self.c_gamma > 0.0 && self.c_gamma.is_finite()) {
            return Err(bad("chemistry.edc.c_gamma", "must be a number > 0"));
        }
        if !(self.c_tau > 0.0 && self.c_tau.is_finite()) {
            return Err(bad("chemistry.edc.c_tau", "must be a number > 0"));
        }
        if !(self.max_fine_fraction > 0.0 && self.max_fine_fraction < 1.0) {
            return Err(bad("chemistry.edc.max_fine_fraction", "must lie in (0, 1)"));
        }
        Ok(())
    }
}

/// The fine structures of a cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FineStructure {
    /// The length scale γ_λ (unclipped).
    pub gamma_lambda: f64,
    /// The mass fraction γ* (held below the guard).
    pub fraction: f64,
    /// The residence time τ* (s).
    pub tau: f64,
    /// κ = γ*/(τ* (1 − γ*)) (1/s).
    pub kappa: f64,
    /// γ* was held at `max_fine_fraction`.
    pub clipped: bool,
}

/// The fine structures of a cell with the kinematic viscosity `nu` (m²/s),
/// k (m²/s²) and ε (m²/s³); `None` without turbulence (k or ε not > 0).
pub fn fine_structure(e: &EdcSection, nu: f64, k: f64, eps: f64) -> Option<FineStructure> {
    if !(k > 0.0 && eps > 0.0 && nu > 0.0) {
        return None;
    }
    let gl = e.c_gamma * (nu * eps / (k * k)).powf(0.25);
    let g = match e.version {
        EdcVersion::V2005 => gl * gl,
        EdcVersion::V1981 => gl * gl * gl,
    };
    let clipped = !(g <= e.max_fine_fraction);
    let g = if clipped { e.max_fine_fraction } else { g };
    let tau = e.c_tau * (nu / eps).sqrt();
    let v = FineStructure {
        gamma_lambda: gl,
        fraction: g,
        tau,
        kappa: g / (tau * (1.0 - g)),
        clipped,
    };
    (v.kappa.is_finite() && v.tau.is_finite()).then_some(v)
}

/// `[chemistry.edm]`: the constants of the eddy dissipation model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EdmSection {
    /// A of the mixing rate A ε/k.
    #[serde(default = "edm_a")]
    pub a: f64,
    /// B of the product term.
    #[serde(default = "edm_b")]
    pub b: f64,
    /// The rate of a reaction is at most that of its kinetics.
    #[serde(default = "yes")]
    pub kinetic_limit: bool,
}

fn edm_a() -> f64 {
    4.0
}
fn edm_b() -> f64 {
    0.5
}
fn yes() -> bool {
    true
}

impl Default for EdmSection {
    fn default() -> Self {
        Self {
            a: edm_a(),
            b: edm_b(),
            kinetic_limit: true,
        }
    }
}

impl EdmSection {
    /// The checks of the constants and of the mechanism (CHM001): the
    /// model needs reactions in one direction.
    pub fn check(&self, prog: &kc::Program<'_>, mech: &str) -> Result<()> {
        if !(self.a > 0.0 && self.a.is_finite()) {
            return Err(bad("chemistry.edm.a", "must be a number > 0"));
        }
        if !(self.b > 0.0 && self.b.is_finite()) {
            return Err(bad("chemistry.edm.b", "must be a number > 0"));
        }
        if let Some(j) = prog.reactions.iter().position(|r| r.reversible != 0) {
            return Err(bad(
                "chemistry.interaction",
                format!(
                    "the eddy dissipation model needs reactions in one direction (`=>`); reaction \
                     {} of `{mech}` is reversible: take interaction = \"edc\", or a global \
                     mechanism such as `westbrook-dryer-2`",
                    j + 1
                ),
            ));
        }
        Ok(())
    }
}

/// The length of the work array of [`edm_rates`] for `n` species.
pub const fn edm_work_len(n: usize) -> usize {
    2 * n + kc::rates_work_len(n)
}

/// The molar production rates ω̇_k (mol/(m³ s)) of the eddy dissipation
/// model in a cell of the mixture (Y, T) at the pressure `p` with the
/// turbulence rate ε/k (1/s), written to `wdot`; `work` holds
/// [`edm_work_len`] values. The density is that of the ideal gas.
#[allow(clippy::too_many_arguments)]
pub fn edm_rates(
    prog: &kc::Program<'_>,
    e: &EdmSection,
    p: f64,
    t: f64,
    y: &[f64],
    eps_over_k: f64,
    wdot: &mut [f64],
    work: &mut [f64],
) {
    let n = prog.n();
    let (c, rest) = work.split_at_mut(n);
    let (w1, rw) = rest.split_at_mut(n);
    let (rho, _) = kc::concentrations(prog.species, p, t, y, c);
    wdot[..n].iter_mut().for_each(|w| *w = 0.0);
    if !(eps_over_k > 0.0) || !(t > 0.0) {
        return;
    }
    let nr = prog.reactions.len();
    let terms_of = |j: usize| {
        let r = &prog.reactions[j];
        &prog.terms[r.terms[0] as usize..r.terms[1] as usize]
    };
    // the mixing rate of reaction j in its own direction (reactants, then
    // the product term), or with `back` in the opposite one
    let mixing = |j: usize, back: bool| {
        let mut lim = f64::INFINITY;
        let (mut yp, mut wp) = (0.0, 0.0);
        for s in terms_of(j) {
            let k = s.species as usize;
            let (nr_, np_) = if back {
                (s.product, s.reactant)
            } else {
                (s.reactant, s.product)
            };
            if nr_ > 0.0 {
                lim = lim.min(c[k].max(0.0) / nr_);
            }
            if np_ > 0.0 {
                yp += y[k].max(0.0);
                wp += np_ * prog.species[k].molar_mass;
            }
        }
        if wp > 0.0 {
            lim = lim.min(e.b * rho * yp / wp);
        }
        if lim.is_finite() {
            e.a * eps_over_k * lim
        } else {
            0.0
        }
    };
    // the rate of the mechanism of reaction j alone (mol/(m³ s)): its
    // production rates, divided by the net coefficient of one of its species
    let mut kinetic = |j: usize| {
        let one = kc::Program {
            reactions: &prog.reactions[j..j + 1],
            ..*prog
        };
        kc::rates(&one, t, c, w1, None, None, rw);
        terms_of(j)
            .iter()
            .find(|s| s.product != s.reactant)
            .map_or(0.0, |s| w1[s.species as usize] / (s.product - s.reactant))
            .max(0.0)
    };
    // a reaction and its exact reverse (the CO₂ dissociation of
    // Westbrook & Dryer) form a pair: with the kinetic limit the mixing
    // rate caps their net rate, in the direction of the kinetics, so that
    // the pair keeps the equilibrium of its kinetics where mixing is fast
    // and a mixing-limited forward step is not undone by an unlimited
    // reverse one
    let reverse_of = |j: usize| {
        (0..nr).find(|&i| {
            i != j && {
                let (a, b) = (terms_of(j), terms_of(i));
                let net = |t: &[kc::Term]| -> Vec<(u32, f64, f64)> {
                    let mut v: Vec<(u32, f64, f64)> = t
                        .iter()
                        .filter(|s| s.reactant != 0.0 || s.product != 0.0)
                        .map(|s| (s.species, s.reactant, s.product))
                        .collect();
                    v.sort_by_key(|x| x.0);
                    v
                };
                let (na, nb) = (net(a), net(b));
                na.len() == nb.len()
                    && na
                        .iter()
                        .zip(&nb)
                        .all(|(x, z)| x.0 == z.0 && x.1 == z.2 && x.2 == z.1)
            }
        })
    };
    for j in 0..nr {
        let pair = if e.kinetic_limit { reverse_of(j) } else { None };
        let q = match pair {
            // the pair is done with its first reaction
            Some(i) if i < j => continue,
            Some(i) => {
                let net = kinetic(j) - kinetic(i);
                if net >= 0.0 {
                    net.min(mixing(j, false))
                } else {
                    -(-net).min(mixing(j, true))
                }
            }
            None => {
                let q = mixing(j, false);
                if e.kinetic_limit {
                    q.min(kinetic(j))
                } else {
                    q
                }
            }
        };
        for s in terms_of(j) {
            wdot[s.species as usize] += (s.product - s.reactant) * q;
        }
    }
}

/// The absolute enthalpy h_i(T) per kg (J/kg) of each species.
pub fn enthalpies(sp: &[Species], t: f64) -> impl Iterator<Item = f64> + '_ {
    sp.iter()
        .map(move |s| s.h_over_rt(t) * foehncraft_kernels::species::R * t / s.molar_mass)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::Mechanism;

    /// The formulas of the fine structures, the versions and the guard.
    #[test]
    fn fine_structures_follow_the_cascade() {
        let e = EdcSection::default();
        let (nu, k, eps) = (1.5e-5, 2.0, 400.0);
        let f = fine_structure(&e, nu, k, eps).unwrap();
        let gl = 2.1377 * (nu * eps / (k * k)).powf(0.25);
        assert!((f.gamma_lambda - gl).abs() < 1e-15);
        assert!((f.fraction - gl * gl).abs() < 1e-15);
        assert!((f.tau - 0.4083 * (nu / eps).sqrt()).abs() < 1e-18);
        assert!((f.kappa - gl * gl / (f.tau * (1.0 - gl * gl))).abs() < 1e-9 * f.kappa);
        assert!(!f.clipped);
        let old = EdcSection {
            version: EdcVersion::V1981,
            ..e.clone()
        };
        let g = fine_structure(&old, nu, k, eps).unwrap();
        assert!((g.fraction - gl * gl * gl).abs() < 1e-15);
        // weak turbulence: γ* is held, κ = 3/τ* at the default guard 0.75
        let w = fine_structure(&e, 1e-3, 1e-4, 1e-2).unwrap();
        assert!(w.clipped && w.fraction == 0.75);
        assert!((w.kappa * w.tau - 3.0).abs() < 1e-12);
        assert!(fine_structure(&e, nu, 0.0, eps).is_none());
        assert!(fine_structure(&e, nu, k, -1.0).is_none());
    }

    /// The eddy dissipation model of the two-step mechanism: the fuel
    /// step is limited by the scarcer of fuel and oxygen, the products
    /// hold it back where there are none, the CO/CO₂ pair is capped by
    /// mixing as a pair, and every rate conserves the mass.
    #[test]
    fn edm_rates_are_mixing_limited_and_conserve_mass() {
        let m = Mechanism::load("westbrook-dryer-2", None).unwrap();
        let prog = m.program();
        let n = m.len();
        let mut y = vec![0.0; n];
        let ix = |s: &str| m.index(s).unwrap();
        y[ix("CH4")] = 0.03;
        y[ix("O2")] = 0.20;
        y[ix("CO2")] = 0.05;
        y[ix("H2O")] = 0.04;
        y[ix("N2")] = 1.0 - 0.32;
        let e = EdmSection {
            kinetic_limit: false,
            ..EdmSection::default()
        };
        let mut w = vec![0.0; n];
        let mut work = vec![0.0; edm_work_len(n)];
        let p = 101_325.0;
        edm_rates(&prog, &e, p, 1500.0, &y, 50.0, &mut w, &mut work);
        let mut c = vec![0.0; n];
        let (rho, _) = kc::concentrations(prog.species, p, 1500.0, &y, &mut c);
        // fuel step: min(C_CH4/1, C_O2/1.5, B ρ Y_P/(W_CO + 2 W_H2O)) with
        // no CO present, Y_P = Y_H2O
        let wco = prog.species[ix("CO")].molar_mass;
        let wh2o = prog.species[ix("H2O")].molar_mass;
        let q1 = 4.0
            * 50.0
            * (c[ix("CH4")])
                .min(c[ix("O2")] / 1.5)
                .min(0.5 * rho * y[ix("H2O")] / (wco + 2.0 * wh2o));
        assert!(
            (-w[ix("CH4")] - q1).abs() < 1e-9 * q1,
            "{} {q1}",
            w[ix("CH4")]
        );
        // mass: Σ W_k ω̇_k = 0
        let mass: f64 = (0..n).map(|k| prog.species[k].molar_mass * w[k]).sum();
        let scale: f64 = (0..n)
            .map(|k| (prog.species[k].molar_mass * w[k]).abs())
            .sum();
        assert!(mass.abs() < 1e-12 * scale);
        // no products: nothing burns
        let mut fresh = y.clone();
        fresh[ix("CO2")] = 0.0;
        fresh[ix("H2O")] = 0.0;
        fresh[ix("N2")] = 1.0 - 0.23;
        edm_rates(&prog, &e, p, 1500.0, &fresh, 50.0, &mut w, &mut work);
        assert!(w[ix("CH4")] == 0.0);
        // the kinetic limit: cold gas does not react
        let kin = EdmSection::default();
        edm_rates(&prog, &kin, p, 300.0, &y, 50.0, &mut w, &mut work);
        assert!(-w[ix("CH4")] < 1e-12 * q1);
        // the CO oxidation and the CO₂ dissociation are a pair: in burnt
        // gas without CO and with slow mixing the net dissociation is held
        // to the mixing rate of CO₂ towards CO and O₂ (its products, O₂,
        // are there), not to its Arrhenius rate alone
        let mut burnt = vec![0.0; n];
        burnt[ix("O2")] = 0.05;
        burnt[ix("CO2")] = 0.12;
        burnt[ix("H2O")] = 0.10;
        burnt[ix("N2")] = 1.0 - 0.27;
        let slow = 1.0;
        edm_rates(&prog, &kin, p, 2000.0, &burnt, slow, &mut w, &mut work);
        let (rho_b, _) = kc::concentrations(prog.species, p, 2000.0, &burnt, &mut c);
        let wo2 = prog.species[ix("O2")].molar_mass;
        let cap = 4.0 * slow * c[ix("CO2")].min(0.5 * rho_b * burnt[ix("O2")] / (wco + 0.5 * wo2));
        assert!(w[ix("CO")] > 0.0 && w[ix("CO")] <= cap * (1.0 + 1e-12));
        let mass: f64 = (0..n).map(|k| prog.species[k].molar_mass * w[k]).sum();
        assert!(mass.abs() < 1e-12 * cap * wco);
        // reversible reactions are refused
        let jl = Mechanism::load("jones-lindstedt-4", None).unwrap();
        assert!(kin.check(&jl.program(), "jones-lindstedt-4").is_err());
        assert!(kin.check(&prog, "westbrook-dryer-2").is_ok());
    }
}
