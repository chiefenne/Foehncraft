//! Discretisation uncertainty of level and time-step sweeps (D67).
//!
//! For every monitored quantity of a sweep with three or more solutions
//! (boundary mass flows, mean pressures and forces, force coefficients,
//! wall probes and wall crossings, benchmark values, viscous and
//! numerical dissipation) the report gives
//! the observed order, the extrapolated value, the uncertainty band and
//! the grid convergence index of the three finest solutions, after the procedure of
//! Celik, Ghia, Roache, Freitas, Coleman and Raad (2008), with the verdict
//! of the ratio of the successive differences:
//!
//! - `monotonic`: 0 < R < 1, the differences fall and keep their sign;
//! - `oscillatory`: −1 < R < 0, they fall and change sign (the order is
//!   that of |ε|, the GCI is an estimate);
//! - `divergent`: |R| ≥ 1, the differences do not fall (no order);
//! - `unchanged`: the solutions agree to round-off (nothing to estimate).
//!
//! The sizes are h ∝ N^(−1/d) for a level sweep (N cells, d the dimension
//! of the mesh) and the time step for a time-step sweep; only their ratios
//! enter. The error norms against an exact solution have their observed
//! orders of their own and are not repeated here.

use serde::Serialize;
use std::fmt::Write as _;

/// The verdict of the ratio of the successive differences.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Verdict {
    /// 0 < R < 1.
    Monotonic,
    /// −1 < R < 0.
    Oscillatory,
    /// |R| ≥ 1.
    Divergent,
    /// The three values agree to round-off.
    Unchanged,
}

impl Verdict {
    /// Name.
    pub fn name(self) -> &'static str {
        match self {
            Self::Monotonic => "monotonic",
            Self::Oscillatory => "oscillatory",
            Self::Divergent => "divergent",
            Self::Unchanged => "unchanged",
        }
    }
}

/// The estimate for one quantity.
#[derive(Debug, Clone, PartialEq, Serialize, schemars::JsonSchema)]
pub struct Estimate {
    /// The quantity (`boundaries.<name>.force_x`, `benchmark.C_D` …).
    pub quantity: String,
    /// Its values, finest solution first (φ1, φ2, φ3).
    pub values: [f64; 3],
    /// The refinement ratios r21 = h2/h1 and r32 = h3/h2.
    pub ratios: [f64; 2],
    /// Observed order p (None when divergent or unchanged).
    pub order: Option<f64>,
    /// The extrapolated value φ_ext (Richardson, from the two finest).
    pub extrapolated: Option<f64>,
    /// Approximate relative error |(φ1 − φ2)/φ1|.
    pub relative_error: f64,
    /// Extrapolated relative error |(φ_ext − φ1)/φ_ext|.
    pub extrapolated_error: Option<f64>,
    /// Grid convergence index of the finest solution,
    /// 1.25 e_a / (r21^p − 1), relative to |φ1|.
    pub gci: Option<f64>,
    /// The same band in the units of the quantity,
    /// 1.25 |φ1 − φ2| / (r21^p − 1): φ1 ± band. For a quantity whose
    /// limit is near zero (a mean pressure at an arbitrary level, a
    /// vanishing force component) only this absolute band is meaningful.
    pub band: Option<f64>,
    /// Verdict.
    pub verdict: Verdict,
}

/// The Celik et al. (2008) estimate from the sizes `h` (finest first,
/// increasing) and the values `phi` of a quantity.
pub fn celik(quantity: &str, h: [f64; 3], phi: [f64; 3]) -> Estimate {
    let (r21, r32) = (h[1] / h[0], h[2] / h[1]);
    let e21 = phi[1] - phi[0];
    let e32 = phi[2] - phi[1];
    let scale = phi.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let relative_error = if phi[0] != 0.0 {
        (e21 / phi[0]).abs()
    } else {
        e21.abs()
    };
    let mut est = Estimate {
        quantity: quantity.to_string(),
        values: phi,
        ratios: [r21, r32],
        order: None,
        extrapolated: None,
        relative_error,
        extrapolated_error: None,
        gci: None,
        band: None,
        verdict: Verdict::Unchanged,
    };
    let tiny = 1e-12 * scale.max(f64::MIN_POSITIVE);
    if e21.abs() <= tiny && e32.abs() <= tiny {
        return est;
    }
    if e32 == 0.0 || e21 == 0.0 {
        // one pair agrees and the other not: no ratio to judge
        est.verdict = Verdict::Divergent;
        return est;
    }
    let ratio = e21 / e32;
    est.verdict = if ratio >= 1.0 || ratio <= -1.0 {
        Verdict::Divergent
    } else if ratio > 0.0 {
        Verdict::Monotonic
    } else {
        Verdict::Oscillatory
    };
    if est.verdict == Verdict::Divergent {
        return est;
    }
    // the order by fixed-point iteration (equal ratios: q = 0, one step)
    let s = (e32 / e21).signum();
    let lr = (e32 / e21).abs().ln();
    let mut p = lr / r21.ln();
    for _ in 0..200 {
        let q = ((r21.powf(p) - s) / (r32.powf(p) - s)).ln();
        let next = ((lr + q) / r21.ln()).abs();
        if !next.is_finite() {
            break;
        }
        let done = (next - p).abs() <= 1e-12 * next.max(1.0);
        p = next;
        if done {
            break;
        }
    }
    if !(p.is_finite() && p > 0.0) {
        return est;
    }
    let rp = r21.powf(p);
    let ext = (rp * phi[0] - phi[1]) / (rp - 1.0);
    est.order = Some(p);
    est.extrapolated = Some(ext);
    est.extrapolated_error = Some(if ext != 0.0 {
        ((ext - phi[0]) / ext).abs()
    } else {
        (ext - phi[0]).abs()
    });
    est.gci = Some(1.25 * relative_error / (rp - 1.0));
    est.band = Some(1.25 * e21.abs() / (rp - 1.0));
    est
}

/// The monitored quantities of one solve (its JSON form), by name. The
/// error norms against an exact solution are not among them: their limit
/// is zero and their observed order is reported on its own.
pub fn quantities(solve: &serde_json::Value) -> Vec<(String, f64)> {
    let mut out = Vec::new();
    let num = |v: &serde_json::Value| v.as_f64().filter(|x| x.is_finite());
    if let Some(ps) = solve.get("boundaries").and_then(|p| p.as_array()) {
        for p in ps {
            let name = p.get("name").and_then(|n| n.as_str()).unwrap_or("?");
            for k in ["mass_flow", "mean_pressure"] {
                if let Some(v) = p.get(k).and_then(num) {
                    out.push((format!("boundaries.{name}.{k}"), v));
                }
            }
            if let Some(f) = p.get("force").and_then(|f| f.as_array()) {
                for (c, v) in ["x", "y", "z"].iter().zip(f) {
                    if let Some(v) = num(v) {
                        out.push((format!("boundaries.{name}.force_{c}"), v));
                    }
                }
            }
        }
    }
    if let Some(c) = solve.get("coefficients").and_then(|c| c.as_array()) {
        for (k, v) in ["drag", "lift"].iter().zip(c) {
            if let Some(v) = num(v) {
                out.push((format!("coefficients.{k}"), v));
            }
        }
    }
    if let Some(ps) = solve.get("wall_probes").and_then(|p| p.as_array()) {
        for p in ps {
            let name = p.get("name").and_then(|n| n.as_str()).unwrap_or("?");
            for k in ["magnitude", "skin_friction"] {
                if let Some(v) = p.get(k).and_then(num) {
                    out.push((format!("wall_probes.{name}.{k}"), v));
                }
            }
            // the passive scalars at the probe: wall value and flux
            for sc in p
                .get("scalars")
                .and_then(|s| s.as_array())
                .into_iter()
                .flatten()
            {
                let s = sc.get("name").and_then(|n| n.as_str()).unwrap_or("?");
                for k in ["wall", "flux"] {
                    if let Some(v) = sc.get(k).and_then(num) {
                        out.push((format!("wall_probes.{name}.{s}_{k}"), v));
                    }
                }
            }
        }
    }
    if let Some(ps) = solve.get("wall_crossings").and_then(|p| p.as_array()) {
        for p in ps {
            let name = p.get("name").and_then(|n| n.as_str()).unwrap_or("?");
            if let Some(v) = p.get("position").and_then(num) {
                out.push((format!("wall_crossings.{name}.position"), v));
            }
        }
    }
    if let Some(ps) = solve.get("periodic").and_then(|p| p.as_array()) {
        for p in ps {
            let name = p.get("from").and_then(|n| n.as_str()).unwrap_or("?");
            for k in ["mass_flow", "driving_force"] {
                if let Some(v) = p.get(k).and_then(num) {
                    out.push((format!("periodic.{name}.{k}"), v));
                }
            }
        }
    }
    if let Some(vals) = solve
        .get("benchmark")
        .and_then(|b| b.get("values"))
        .and_then(|v| v.as_array())
    {
        for b in vals {
            if let (Some(n), Some(v)) = (
                b.get("range")
                    .and_then(|r| r.get("name"))
                    .and_then(|n| n.as_str()),
                b.get("value").and_then(num),
            ) {
                out.push((format!("benchmark.{n}"), v));
            }
        }
    }
    if let Some(k) = solve
        .get("balance")
        .and_then(|b| b.get("last"))
        .and_then(|l| l.get("kinetic_energy"))
    {
        for q in ["viscous_dissipation", "numerical_dissipation"] {
            if let Some(v) = k.get(q).and_then(num) {
                out.push((format!("balance.{q}"), v));
            }
        }
    }
    out
}

/// The estimates of a sweep: `solves` in the order of the run with their
/// sizes (`h ∝ N^(−1/d)` or the time step); the three finest solutions
/// are used. Quantities that are zero to round-off in all three are left
/// out.
pub fn sweep(solves: &[serde_json::Value], sizes: &[f64]) -> Vec<Estimate> {
    if solves.len() < 3 || solves.len() != sizes.len() {
        return Vec::new();
    }
    let mut idx: Vec<usize> = (0..sizes.len()).collect();
    idx.sort_by(|&a, &b| sizes[a].total_cmp(&sizes[b]));
    let pick = [idx[0], idx[1], idx[2]];
    let h = [sizes[pick[0]], sizes[pick[1]], sizes[pick[2]]];
    if !(h[0] < h[1] && h[1] < h[2]) {
        return Vec::new();
    }
    let q: Vec<Vec<(String, f64)>> = pick.iter().map(|&i| quantities(&solves[i])).collect();
    let mut out = Vec::new();
    for (name, v1) in &q[0] {
        let find = |l: &[(String, f64)]| l.iter().find(|(n, _)| n == name).map(|(_, v)| *v);
        let (Some(v2), Some(v3)) = (find(&q[1]), find(&q[2])) else {
            continue;
        };
        let phi = [*v1, v2, v3];
        if phi.iter().all(|v| v.abs() <= 1e-12) {
            continue;
        }
        out.push(celik(name, h, phi));
    }
    out
}

/// The estimates as a Markdown table.
pub fn markdown(est: &[Estimate], what: &str) -> String {
    let mut s = String::new();
    if est.is_empty() {
        return s;
    }
    let _ = writeln!(
        s,
        "\nDiscretisation uncertainty of the three finest {what} (Celik et al. 2008, D67):\n"
    );
    let _ = writeln!(
        s,
        "| Quantity | φ1 (finest) | φ2 | φ3 | Order | Extrapolated | ± band | GCI fine | Verdict |"
    );
    let _ = writeln!(s, "|---|---|---|---|---|---|---|---|---|");
    let f = |v: Option<f64>, fmt: &dyn Fn(f64) -> String| v.map_or("-".to_string(), fmt);
    for e in est {
        let _ = writeln!(
            s,
            "| {} | {:.6e} | {:.6e} | {:.6e} | {} | {} | {} | {} | {} |",
            e.quantity,
            e.values[0],
            e.values[1],
            e.values[2],
            f(e.order, &|p| format!("{p:.3}")),
            f(e.extrapolated, &|x| format!("{x:.6e}")),
            f(e.band, &|b| format!("{b:.3e}")),
            // relative to φ1 only where φ1 is larger than its band
            f(
                e.gci
                    .filter(|_| e.band.is_some_and(|b| b < e.values[0].abs())),
                &|g| format!("{:.3} %", 100.0 * g)
            ),
            e.verdict.name()
        );
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn second_order_data_give_order_two_and_the_exact_value() {
        // φ(h) = 1 + 0.3 h², h = 1, 2, 4
        let f = |h: f64| 1.0 + 0.3 * h * h;
        let e = celik("q", [1.0, 2.0, 4.0], [f(1.0), f(2.0), f(4.0)]);
        assert_eq!(e.verdict, Verdict::Monotonic);
        assert!((e.order.unwrap() - 2.0).abs() < 1e-12);
        assert!((e.extrapolated.unwrap() - 1.0).abs() < 1e-12);
        // GCI = 1.25 |φ1 − φ2|/|φ1| / (r^p − 1)
        let ea = (f(2.0) - f(1.0)) / f(1.0);
        assert!((e.gci.unwrap() - 1.25 * ea / 3.0).abs() < 1e-14);
    }

    #[test]
    fn unequal_ratios_need_the_fixed_point() {
        // φ(h) = 2 − 0.1 h^1.5 at h = 1, 1.5, 3.3 (ratios 1.5 and 2.2)
        let f = |h: f64| 2.0 - 0.1 * h.powf(1.5);
        let e = celik("q", [1.0, 1.5, 3.3], [f(1.0), f(1.5), f(3.3)]);
        assert_eq!(e.verdict, Verdict::Monotonic);
        assert!((e.order.unwrap() - 1.5).abs() < 1e-9, "{:?}", e.order);
        assert!((e.extrapolated.unwrap() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn verdicts() {
        let osc = celik("q", [1.0, 2.0, 4.0], [1.0, 1.1, 0.7]);
        assert_eq!(osc.verdict, Verdict::Oscillatory);
        assert!(osc.order.is_some());
        let div = celik("q", [1.0, 2.0, 4.0], [1.0, 1.5, 1.6]);
        assert_eq!(div.verdict, Verdict::Divergent);
        assert!(div.order.is_none() && div.gci.is_none());
        let same = celik("q", [1.0, 2.0, 4.0], [3.0, 3.0, 3.0]);
        assert_eq!(same.verdict, Verdict::Unchanged);
    }

    #[test]
    fn quantities_of_a_solve() {
        let s = serde_json::json!({
            "errors": {"l2": 1e-3, "linf": 2e-3, "pressure_linf": 0.1, "h": 0.5, "cells": 10},
            "boundaries": [{"name": "wall", "area": 1.0, "mass_flow": 0.0, "mean_pressure": 2.0,
                         "force": [0.5, 0.0, 0.0]}],
            "benchmark": {"values": [{"range": {"name": "C_D", "lower": 5.5, "upper": 5.6},
                                      "value": 5.58, "deviation": 0.0}]},
        });
        let q: Vec<String> = quantities(&s).into_iter().map(|(n, _)| n).collect();
        assert_eq!(
            q,
            [
                "boundaries.wall.mass_flow",
                "boundaries.wall.mean_pressure",
                "boundaries.wall.force_x",
                "boundaries.wall.force_y",
                "boundaries.wall.force_z",
                "benchmark.C_D"
            ]
        );
        let sw = sweep(&[s.clone(), s.clone(), s], &[4.0, 2.0, 1.0]);
        // zero force components and the zero mass flow are left out
        assert!(sw.iter().all(|e| e.verdict == Verdict::Unchanged));
        assert_eq!(sw.len(), 3);
    }

    #[test]
    fn band_in_the_units_of_the_quantity() {
        // φ = 1 + h², h = 1, 2, 4: p = 2, φ_ext = 1, band = 1.25 · 3 / 3
        let e = celik("q", [1.0, 2.0, 4.0], [2.0, 5.0, 17.0]);
        assert!((e.order.unwrap() - 2.0).abs() < 1e-10);
        assert!((e.extrapolated.unwrap() - 1.0).abs() < 1e-12);
        assert!((e.band.unwrap() - 1.25).abs() < 1e-12);
        assert!((e.gci.unwrap() - 1.25 * 1.5 / 3.0).abs() < 1e-12);
        let md = markdown(&[e], "meshes");
        assert!(md.contains("| q |") && md.contains("1.250e0"));
    }
}
