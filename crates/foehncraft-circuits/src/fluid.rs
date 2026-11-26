//! Water and steam states for the circuits, on IAPWS-IF97 (stage SF7a,
//! `foehncraft_kernels::steam`).

use crate::{Error, codes};
use foehncraft_kernels::steam::{self, MixtureViscosity, Phase, Saturation, State};

/// The lowest specific enthalpy of a circuit state (J/kg): liquid at
/// 273.16 K.
pub const H_MIN: f64 = 1.0e3;

/// A state of the flowing water or steam with the properties the
/// correlations need. In a two-phase state the liquid and vapour
/// properties are those of the saturation line.
#[derive(Clone, Copy, Debug)]
pub struct Point {
    /// Pressure (Pa).
    pub p: f64,
    /// Specific enthalpy of the mixture (J/kg).
    pub h: f64,
    /// The state at (p, h) (homogeneous density in a two-phase state).
    pub state: State,
    /// The saturation line at p, `None` at and above the critical
    /// pressure.
    pub sat: Option<Saturation>,
}

impl Point {
    /// The state at (p, h).
    pub fn new(p: f64, h: f64) -> Result<Self, Error> {
        let h = h.max(H_MIN);
        let state = steam::state_ph(p, h).ok_or_else(|| {
            Error::new(
                codes::RANGE,
                format!(
                    "no IAPWS-IF97 state at p = {:.6e} Pa, h = {:.6e} J/kg (pressure 611 Pa to 100 MPa, temperature 273.15 K to 1073.15 K)",
                    p, h
                ),
            )
        })?;
        let sat = steam::saturation_p(p);
        Ok(Self { p, h, state, sat })
    }

    /// The equilibrium quality (h − h′)/(h″ − h′); NaN at and above p_c.
    pub fn x_eq(&self) -> f64 {
        match &self.sat {
            Some(s) => (self.h - s.liquid.h) / s.h_fg(),
            None => f64::NAN,
        }
    }

    /// Temperature (K).
    pub fn t(&self) -> f64 {
        self.state.t
    }

    /// Saturation temperature (K), NaN at and above p_c.
    pub fn t_sat(&self) -> f64 {
        self.sat.as_ref().map_or(f64::NAN, |s| s.t)
    }

    /// Whether the state is single phase (liquid, vapour or supercritical).
    pub fn single_phase(&self) -> bool {
        self.state.phase != Phase::TwoPhase
    }
}

/// Transport properties of one phase.
#[derive(Clone, Copy, Debug)]
pub struct Phase1 {
    /// Density (kg/m³).
    pub rho: f64,
    /// Dynamic viscosity (Pa s).
    pub mu: f64,
    /// Thermal conductivity (W/(m K)).
    pub k: f64,
    /// Isobaric specific heat (J/(kg K)).
    pub cp: f64,
}

impl Phase1 {
    /// The properties of a single-phase state.
    pub fn of(s: &State) -> Self {
        Self {
            rho: s.rho,
            mu: steam::viscosity(s.rho, s.t),
            k: steam::conductivity(s),
            cp: s.cp,
        }
    }

    /// The Prandtl number.
    pub fn pr(&self) -> f64 {
        self.mu * self.cp / self.k
    }
}

/// The properties of the saturated phases.
#[derive(Clone, Copy, Debug)]
pub struct SatProps {
    /// Saturated liquid.
    pub l: Phase1,
    /// Saturated vapour.
    pub g: Phase1,
    /// Enthalpy of evaporation (J/kg).
    pub h_fg: f64,
    /// Surface tension (N/m).
    pub sigma: f64,
    /// Saturation temperature (K).
    pub t: f64,
    /// Saturated liquid enthalpy (J/kg).
    pub h_l: f64,
}

impl SatProps {
    /// The saturated phase properties of a saturation line.
    pub fn of(sat: &Saturation) -> Self {
        Self {
            l: Phase1::of(&sat.liquid),
            g: Phase1::of(&sat.vapour),
            h_fg: sat.h_fg(),
            sigma: steam::surface_tension(sat.t),
            t: sat.t,
            h_l: sat.liquid.h,
        }
    }
}

/// The homogeneous mixture viscosity (McAdams) at the quality x.
pub fn mixture_viscosity(sat: &Saturation, x: f64) -> f64 {
    steam::mixture_viscosity(sat, x.clamp(0.0, 1.0), MixtureViscosity::McAdams)
}

/// The saturation pressure at the temperature t (K), clamped to the
/// saturation line (273.15 K to T_c).
pub fn p_sat(t: f64) -> f64 {
    steam::p_sat(t.clamp(steam::T_MIN, steam::T_CRIT - 1e-9))
}

/// The critical pressure (Pa).
pub const P_CRIT: f64 = steam::P_CRIT;
