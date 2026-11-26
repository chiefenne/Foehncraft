//! Water and steam states of IAPWS-IF97 (stage SF7a) for the tools: one
//! state from any pair of (p, T), (p, h), (p, s), (p, x), (T, x), or the
//! saturation line at p or T, with the transport properties. The
//! formulation itself is the single-source kernel
//! [`foehncraft_kernels::steam`], re-exported here.

pub use foehncraft_kernels::steam::*;

use crate::diag::{self, Diagnostic};
use anyhow::Result;
use serde_json::{Value, json};

/// The inputs of a state query (SI units); exactly one of the pairs (p, T),
/// (p, h), (p, s), (p, x), (T, x), or p or T alone (the saturation line).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Query {
    /// Pressure (Pa).
    pub pressure: Option<f64>,
    /// Temperature (K).
    pub temperature: Option<f64>,
    /// Specific enthalpy (J/kg).
    pub enthalpy: Option<f64>,
    /// Specific entropy (J/(kg K)).
    pub entropy: Option<f64>,
    /// Quality (0 to 1).
    pub quality: Option<f64>,
}

fn err(m: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::MAT_STEAM, m))
}

/// The JSON of one state: the thermodynamic properties, the transport
/// properties of a single-phase state, and for a two-phase state the
/// saturated phases and the mixture viscosities.
pub fn state_json(s: &State) -> Value {
    let num = |v: f64| if v.is_finite() { json!(v) } else { Value::Null };
    let phase = match s.phase {
        Phase::Liquid => "liquid",
        Phase::Vapour => "vapour",
        Phase::TwoPhase => "two-phase",
        Phase::Supercritical => "supercritical",
    };
    let single = s.phase != Phase::TwoPhase;
    let (mu, k) = if single {
        (viscosity(s.rho, s.t), conductivity(s))
    } else {
        (f64::NAN, f64::NAN)
    };
    let mut v = json!({
        "region": s.region,
        "phase": phase,
        "pressure": num(s.p),
        "temperature": num(s.t),
        "density": num(s.rho),
        "specific_volume": num(1.0 / s.rho),
        "enthalpy": num(s.h),
        "internal_energy": num(s.u),
        "entropy": num(s.s),
        "specific_heat": num(s.cp),
        "specific_heat_cv": num(s.cv),
        "speed_of_sound": num(s.w),
        "quality": num(s.x),
        "void_fraction": num(s.alpha),
        "viscosity": num(mu),
        "conductivity": num(k),
        "prandtl": num(mu * s.cp / k),
    });
    if s.t < T_CRIT {
        v["surface_tension"] = num(surface_tension(s.t));
    }
    if !single && let Some(sat) = saturation_p(s.p) {
        v["liquid"] = state_json(&sat.liquid);
        v["vapour"] = state_json(&sat.vapour);
        v["mixture_viscosity"] = json!({
            "mcadams": num(mixture_viscosity(&sat, s.x, MixtureViscosity::McAdams)),
            "cicchitti": num(mixture_viscosity(&sat, s.x, MixtureViscosity::Cicchitti)),
            "dukler": num(mixture_viscosity(&sat, s.x, MixtureViscosity::Dukler)),
        });
    }
    v
}

/// The saturation line at the pressure or the temperature of `sat`.
pub fn saturation_json(sat: &Saturation) -> Value {
    json!({
        "pressure": sat.p,
        "temperature": sat.t,
        "enthalpy_of_evaporation": sat.h_fg(),
        "surface_tension": surface_tension(sat.t),
        "liquid": state_json(&sat.liquid),
        "vapour": state_json(&sat.vapour),
    })
}

/// Evaluates a query: `{"input": …, "state": …}` for a state, `{"input":
/// …, "saturation": …}` for p or T alone.
pub fn query(q: Query) -> Result<Value> {
    let given: Vec<&str> = [
        ("pressure", q.pressure.is_some()),
        ("temperature", q.temperature.is_some()),
        ("enthalpy", q.enthalpy.is_some()),
        ("entropy", q.entropy.is_some()),
        ("quality", q.quality.is_some()),
    ]
    .iter()
    .filter(|g| g.1)
    .map(|g| g.0)
    .collect();
    let outside = |what: String| {
        err(format!(
            "{what} lies outside IAPWS-IF97 (273.15 K to 1073.15 K up to 100 MPa, to 2273.15 K up \
             to 50 MPa; the saturation line from 611.2 Pa to 22.064 MPa)"
        ))
    };
    let input = given.join(", ");
    let state = |s: Option<State>, what: String| -> Result<Value> {
        let s = s.ok_or_else(|| outside(what))?;
        Ok(json!({ "standard": "IAPWS-IF97", "input": input, "state": state_json(&s) }))
    };
    match (q.pressure, q.temperature, q.enthalpy, q.entropy, q.quality) {
        (Some(p), Some(t), None, None, None) => {
            state(state_pt(p, t), format!("p = {p} Pa, T = {t} K"))
        }
        (Some(p), None, Some(h), None, None) => {
            state(state_ph(p, h), format!("p = {p} Pa, h = {h} J/kg"))
        }
        (Some(p), None, None, Some(s), None) => {
            state(state_ps(p, s), format!("p = {p} Pa, s = {s} J/(kg K)"))
        }
        (Some(p), None, None, None, Some(x)) => {
            state(state_px(p, x), format!("p = {p} Pa, x = {x}"))
        }
        (None, Some(t), None, None, Some(x)) => {
            state(state_tx(t, x), format!("T = {t} K, x = {x}"))
        }
        (Some(p), None, None, None, None) => {
            let sat = saturation_p(p).ok_or_else(|| outside(format!("p = {p} Pa")))?;
            Ok(
                json!({ "standard": "IAPWS-IF97", "input": input, "saturation": saturation_json(&sat) }),
            )
        }
        (None, Some(t), None, None, None) => {
            let sat = saturation_t(t).ok_or_else(|| outside(format!("T = {t} K")))?;
            Ok(
                json!({ "standard": "IAPWS-IF97", "input": input, "saturation": saturation_json(&sat) }),
            )
        }
        _ => Err(err(format!(
            "a water or steam state needs (p, T), (p, h), (p, s), (p, x), (T, x), or p or T alone \
             for the saturation line; given: {}",
            if given.is_empty() {
                "nothing".to_string()
            } else {
                input.clone()
            }
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn queries_give_states_and_saturation_lines() {
        let v = query(Query {
            pressure: Some(3e6),
            temperature: Some(300.0),
            ..Query::default()
        })
        .unwrap();
        let h = v["state"]["enthalpy"].as_f64().unwrap();
        assert!((h - 115.331273e3).abs() < 1e-3, "{h}");
        assert_eq!(v["state"]["region"], 1);
        let s = query(Query {
            pressure: Some(101_325.0),
            ..Query::default()
        })
        .unwrap();
        let t = s["saturation"]["temperature"].as_f64().unwrap();
        assert!((t - 373.124).abs() < 2e-3, "{t}");
        let m = query(Query {
            pressure: Some(7e6),
            quality: Some(0.2),
            ..Query::default()
        })
        .unwrap();
        assert_eq!(m["state"]["phase"], "two-phase");
        assert!(m["state"]["mixture_viscosity"]["mcadams"].as_f64().unwrap() > 0.0);
        let e = query(Query {
            pressure: Some(1e6),
            ..Query {
                temperature: Some(400.0),
                enthalpy: Some(1e6),
                ..Query::default()
            }
        })
        .unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MAT_STEAM, "{e:#}");
        let e = query(Query {
            pressure: Some(1e6),
            temperature: Some(3000.0),
            ..Query::default()
        })
        .unwrap_err();
        assert!(format!("{e:#}").contains("outside IAPWS-IF97"), "{e:#}");
    }
}
