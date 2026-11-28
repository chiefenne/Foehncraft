//! The solution of a circuit: per node, per tube and per segment, as JSON
//! (`result.json`) and CSV (`profiles.csv`).

use crate::case::{Case, NodeKind};
use crate::network::{March, Model};
use crate::{Error, codes};
use serde_json::{Value, json};

/// A warning of a solved circuit.
#[derive(Clone, Debug)]
pub struct Warning {
    /// Diagnostic code.
    pub code: &'static str,
    /// Message.
    pub message: String,
}

/// The solved circuit.
#[derive(Clone, Debug)]
pub struct Solution {
    /// The model.
    pub model: Model,
    /// The unknowns.
    pub x: Vec<f64>,
    /// The recorded marches, one per tube.
    pub marches: Vec<March>,
    /// Largest scaled residual per Newton iteration.
    pub history: Vec<f64>,
    /// Warnings (CHF margin, enthalpy limit).
    pub warnings: Vec<Warning>,
}

/// A drum of the solution.
#[derive(Clone, Copy, Debug)]
pub struct DrumBalance {
    /// Steam leaving the drum (kg/s), equal to the feedwater.
    pub steam: f64,
    /// Flow of the risers into the drum (kg/s).
    pub riser_flow: f64,
    /// Circulation ratio: riser flow over steam flow.
    pub circulation_ratio: f64,
    /// Enthalpy of the water to the downcomers (J/kg).
    pub h_downcomer: f64,
}

/// Solve a case.
pub fn run(case: &Case) -> Result<Solution, Error> {
    let model = Model::new(case)?;
    let (x, marches, history) = model.solve()?;
    Ok(Solution::of(model, x, marches, history))
}

impl Solution {
    /// The solution of a solved model, with its warnings.
    pub fn of(model: Model, x: Vec<f64>, marches: Vec<March>, history: Vec<f64>) -> Self {
        let mut warnings = Vec::new();
        for (t, m) in model.tubes.iter().zip(&marches) {
            if m.clipped {
                warnings.push(Warning {
                code: codes::RANGE,
                message: format!(
                    "tube `{}`: the enthalpy reached its limit of 4.1 MJ/kg: the flow is too small for the heat",
                    t.name
                ),
            });
            }
            if let Some((s, r)) = min_chf(m)
                && r < model.chf_margin
            {
                {
                    warnings.push(Warning {
                        code: codes::CHF,
                        message: format!(
                            "tube `{}`: CHF ratio {r:.3} at s = {s:.3} m, below the margin {:.2}{}",
                            t.name,
                            model.chf_margin,
                            if r < 1.0 {
                                " (the wall is beyond the critical heat flux there)"
                            } else {
                                ""
                            }
                        ),
                    });
                }
            }
        }
        Solution {
            model,
            x,
            marches,
            history,
            warnings,
        }
    }
}

/// The smallest CHF ratio of a tube and its position.
pub fn min_chf(m: &March) -> Option<(f64, f64)> {
    m.segments
        .iter()
        .filter(|s| s.chf_ratio.is_finite())
        .map(|s| (s.s, s.chf_ratio))
        .min_by(|a, b| a.1.total_cmp(&b.1))
}

impl Solution {
    /// Flow of one tube (kg/s).
    pub fn flow(&self, tube: &str) -> Option<f64> {
        self.model
            .tubes
            .iter()
            .find(|t| t.name == tube)
            .map(|t| self.x[t.im])
    }

    /// The march of a tube.
    pub fn tube(&self, tube: &str) -> Option<&March> {
        self.model
            .tubes
            .iter()
            .position(|t| t.name == tube)
            .map(|i| &self.marches[i])
    }

    /// Pressure (Pa) of a node.
    pub fn pressure(&self, node: &str) -> Option<f64> {
        self.model
            .nodes
            .iter()
            .position(|n| n.name == node)
            .map(|i| self.model.p(&self.x, i))
    }

    /// Enthalpy (J/kg) of a node.
    pub fn enthalpy(&self, node: &str) -> Option<f64> {
        self.model
            .nodes
            .iter()
            .position(|n| n.name == node)
            .map(|i| self.model.h(&self.x, i))
    }

    /// The balance of a drum.
    pub fn drum(&self, node: &str) -> Option<DrumBalance> {
        let k = self
            .model
            .nodes
            .iter()
            .position(|n| n.name == node && n.kind == NodeKind::Drum)?;
        let p = self.model.p(&self.x, k);
        let sat = foehncraft_kernels::steam::saturation_p(p)?;
        let (hf, hg) = (sat.liquid.h, sat.vapour.h);
        let (mut steam, mut rise) = (0.0, 0.0);
        for (t, mr) in self.model.tubes.iter().zip(&self.marches) {
            let m = self.x[t.im] * t.count;
            let into = (t.to == k && m > 0.0) || (t.from == k && m < 0.0);
            if into {
                rise += m.abs();
                steam += m.abs() * ((mr.h_end - hf) / (hg - hf)).clamp(0.0, 1.0);
            }
        }
        Some(DrumBalance {
            steam,
            riser_flow: rise,
            circulation_ratio: if steam > 0.0 {
                rise / steam
            } else {
                f64::INFINITY
            },
            h_downcomer: self.model.h(&self.x, k),
        })
    }

    /// The enthalpy flows (W) that enter and leave the circuit: the water
    /// of the pressure nodes and of the junction sources, the feedwater
    /// and the steam of the drums. With the heat of the tubes they close
    /// the energy balance of the circuit: in + heat − out = 0 to the
    /// Newton tolerance (the ledger account of SF7c).
    pub fn enthalpy_flows(&self) -> (f64, f64) {
        let m = &self.model;
        let (mut inflow, mut outflow) = (0.0, 0.0);
        for (t, mr) in m.tubes.iter().zip(&self.marches) {
            let flow = self.x[t.im] * t.count;
            let (up, down) = if flow >= 0.0 {
                (t.from, t.to)
            } else {
                (t.to, t.from)
            };
            if m.nodes[up].kind == NodeKind::Pressure {
                inflow += flow.abs() * m.h(&self.x, up);
            }
            if m.nodes[down].kind == NodeKind::Pressure {
                outflow += flow.abs() * mr.h_end;
            }
        }
        for (k, n) in m.nodes.iter().enumerate() {
            if n.m_source > 0.0 {
                inflow += n.m_source * n.h_given.unwrap_or(0.0);
            } else if n.m_source < 0.0 {
                outflow -= n.m_source * m.h(&self.x, k);
            }
            if n.kind == NodeKind::Drum
                && let Some(d) = self.drum(&n.name)
                && let Some(sat) = foehncraft_kernels::steam::saturation_p(m.p(&self.x, k))
            {
                inflow += d.steam * n.h_feed;
                outflow += d.steam * sat.vapour.h;
            }
        }
        (inflow, outflow)
    }

    /// The result as JSON (format `foehncraft.circuit-result` 1.0).
    pub fn to_json(&self, name: &str, case: &Case) -> Value {
        let m = &self.model;
        let nodes: Vec<Value> = m
            .nodes
            .iter()
            .enumerate()
            .map(|(k, n)| {
                let p = m.p(&self.x, k);
                let h = m.h(&self.x, k);
                let t = foehncraft_kernels::steam::state_ph(p, h).map_or(f64::NAN, |s| s.t);
                let mut v = json!({
                    "name": n.name,
                    "kind": format!("{:?}", n.kind).to_lowercase(),
                    "pressure": {"value": p, "unit": "Pa"},
                    "enthalpy": {"value": h, "unit": "J/kg"},
                    "temperature": {"value": t, "unit": "K"},
                });
                if let Some(d) = self.drum(&n.name) {
                    v["drum"] = json!({
                        "steam_flow": {"value": d.steam, "unit": "kg/s"},
                        "riser_flow": {"value": d.riser_flow, "unit": "kg/s"},
                        "circulation_ratio": d.circulation_ratio,
                        "downcomer_enthalpy": {"value": d.h_downcomer, "unit": "J/kg"},
                    });
                }
                v
            })
            .collect();
        let tubes: Vec<Value> = m
            .tubes
            .iter()
            .zip(&self.marches)
            .map(|(t, mr)| {
                let flow = self.x[t.im];
                let out = foehncraft_kernels::steam::state_ph(mr.p_end, mr.h_end);
                let t_max = mr.segments.iter().map(|s| s.t_wall).fold(f64::NAN, f64::max);
                let chf = min_chf(mr);
                let x_out = mr.segments.last().map_or(f64::NAN, |s| s.x_eq);
                json!({
                    "name": t.name,
                    "count": t.count,
                    "flow": {"value": flow, "unit": "kg/s", "of": "one tube"},
                    "flow_total": {"value": flow * t.count, "unit": "kg/s"},
                    "mass_flux": {"value": flow / t.area, "unit": "kg/(m2 s)"},
                    "heat": {"value": t.heat(), "unit": "W", "of": "one tube"},
                    "dp_friction": {"value": mr.dp_fric, "unit": "Pa"},
                    "dp_gravity": {"value": mr.dp_grav, "unit": "Pa"},
                    "dp_acceleration": {"value": mr.dp_acc, "unit": "Pa"},
                    "outlet": {
                        "pressure": {"value": mr.p_end, "unit": "Pa"},
                        "enthalpy": {"value": mr.h_end, "unit": "J/kg"},
                        "temperature": {"value": out.map_or(f64::NAN, |s| s.t), "unit": "K"},
                        "quality_last_segment": x_out,
                    },
                    "wall_temperature_max": {"value": t_max, "unit": "K"},
                    "chf_ratio_min": chf.map(|c| json!({"value": c.1, "at": {"value": c.0, "unit": "m"}})),
                })
            })
            .collect();
        json!({
            "format": "foehncraft.circuit-result",
            "version": "1.0",
            "name": name,
            "models": {
                "friction": case.circuit.friction,
                "void": case.circuit.void,
                "boiling": case.circuit.boiling,
                "subcooled": case.circuit.subcooled,
                "chf": case.circuit.chf,
            },
            "converged": true,
            "newton": {"iterations": self.history.len() - 1, "residuals": self.history},
            "nodes": nodes,
            "tubes": tubes,
            "warnings": self.warnings.iter().map(|w| json!({"code": w.code, "message": w.message})).collect::<Vec<_>>(),
        })
    }

    /// Every segment of every tube as CSV.
    pub fn to_csv(&self) -> String {
        let mut s = String::from(
            "tube,s_m,p_Pa,h_J_kg,T_K,x_eq,x_flow,void,q_W_m2,T_wall_K,htc_W_m2K,regime,q_chf_W_m2,chf_ratio\n",
        );
        for (t, mr) in self.model.tubes.iter().zip(&self.marches) {
            for g in &mr.segments {
                s += &format!(
                    "{},{:.6},{:.6e},{:.6e},{:.6},{:.6},{:.6},{:.6},{:.6e},{:.6},{:.6e},{},{:.6e},{:.6}\n",
                    t.name,
                    g.s,
                    g.p,
                    g.h,
                    g.t,
                    g.x_eq,
                    g.x,
                    g.alpha,
                    g.q,
                    g.t_wall,
                    g.htc,
                    g.regime,
                    g.q_chf,
                    g.chf_ratio
                );
            }
        }
        s
    }
}
