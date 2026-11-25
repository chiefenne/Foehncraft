//! The critical heat flux of water in round tubes.
//!
//! Two sources: the correlation of Bowring (1972), built in, and a look-up
//! table in the layout of Groeneveld et al. (2007), "The 2006 CHF look-up
//! table" (pressure, mass flux, equilibrium quality → CHF of an 8 mm
//! tube), read from a CSV file the user supplies, with the diameter
//! correction of the paper, (0.008/D)^½ for 3 mm ≤ D ≤ 25 mm.
//!
//! Bowring, R. W. (1972). A simple but accurate round tube uniform heat
//! flux, dryout correlation over the pressure range 0.7–17 MN/m².
//! AEEW-R 789, UKAEA Winfrith. In the local-conditions form of Todreas &
//! Kazimi (2012), Nuclear Systems I, section 13.

use crate::{Error, codes};

/// The critical heat flux (W/m²) of Bowring (1972) at the local
/// equilibrium quality `x` (local-conditions form, q = (A − B h_fg x)/C)
/// for pressure `p` (Pa), mass flux `g` (kg/(m² s)), diameter `d` (m)
/// and enthalpy of evaporation `h_fg` (J/kg). Range of the data: 0.2 to
/// 19 MPa, 136 to 18 600 kg/(m² s), 2 to 45 mm.
pub fn bowring(p: f64, g: f64, d: f64, h_fg: f64, x: f64) -> f64 {
    let pr = 0.145e-6 * p;
    let n = 2.0 - 0.5 * pr;
    let (f1, f2, f3, f4) = if pr <= 1.0 {
        let f1 = (pr.powf(18.942) * (20.89 * (1.0 - pr)).exp() + 0.917) / 1.917;
        let f1_f2 = (pr.powf(1.316) * (2.444 * (1.0 - pr)).exp() + 0.309) / 1.309;
        let f3 = (pr.powf(17.023) * (16.658 * (1.0 - pr)).exp() + 0.667) / 1.667;
        (f1, f1 / f1_f2, f3, f3 * pr.powf(1.649))
    } else {
        let f1 = pr.powf(-0.368) * (0.648 * (1.0 - pr)).exp();
        let f1_f2 = pr.powf(-0.448) * (0.245 * (1.0 - pr)).exp();
        let f3 = pr.powf(0.219);
        (f1, f1 / f1_f2, f3, f3 * pr.powf(1.649))
    };
    let a = 2.317 * (h_fg * d * g / 4.0) * f1 / (1.0 + 0.0143 * f2 * d.sqrt() * g);
    let b = d * g / 4.0;
    let c = 0.077 * f3 * d * g / (1.0 + 0.347 * f4 * (g / 1356.0).powf(n));
    ((a - b * h_fg * x) / c).max(0.0)
}

/// A CHF look-up table on a grid of pressure, mass flux and equilibrium
/// quality for an 8 mm tube (the layout of the 2006 table).
#[derive(Clone, Debug)]
pub struct LookUpTable {
    /// Pressures (Pa), increasing.
    pub p: Vec<f64>,
    /// Mass fluxes (kg/(m² s)), increasing.
    pub g: Vec<f64>,
    /// Equilibrium qualities, increasing.
    pub x: Vec<f64>,
    /// CHF (W/m²), index [ip][ig][ix]; NaN where the table has no value.
    pub chf: Vec<f64>,
}

fn grid(v: &mut Vec<f64>) {
    v.sort_by(f64::total_cmp);
    v.dedup();
}

impl LookUpTable {
    /// Read a table from CSV text with the header `p_kpa,g,x,chf_kw`:
    /// pressure in kPa, mass flux in kg/(m² s), equilibrium quality and
    /// CHF in kW/m² (the units of the 2006 table). Every grid point must
    /// be given once; empty CHF cells are NaN.
    pub fn from_csv(text: &str) -> Result<Self, Error> {
        let bad = |m: String| Error::new(codes::CASE, format!("CHF table: {m}"));
        let mut rows = Vec::new();
        let mut lines = text
            .lines()
            .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
        let head = lines.next().ok_or_else(|| bad("empty".into()))?;
        let cols: Vec<&str> = head.split(',').map(str::trim).collect();
        if cols != ["p_kpa", "g", "x", "chf_kw"] {
            return Err(bad(format!(
                "header must be p_kpa,g,x,chf_kw, found {head}"
            )));
        }
        for (i, l) in lines.enumerate() {
            let f: Vec<&str> = l.split(',').map(str::trim).collect();
            if f.len() != 4 {
                return Err(bad(format!("line {}: four values expected", i + 2)));
            }
            let num = |s: &str| -> Result<f64, Error> {
                if s.is_empty() {
                    return Ok(f64::NAN);
                }
                s.parse::<f64>()
                    .map_err(|_| bad(format!("line {}: not a number: {s}", i + 2)))
            };
            rows.push([num(f[0])? * 1e3, num(f[1])?, num(f[2])?, num(f[3])? * 1e3]);
        }
        let (mut p, mut g, mut x) = (Vec::new(), Vec::new(), Vec::new());
        for r in &rows {
            p.push(r[0]);
            g.push(r[1]);
            x.push(r[2]);
        }
        grid(&mut p);
        grid(&mut g);
        grid(&mut x);
        let (np, ng, nx) = (p.len(), g.len(), x.len());
        if np < 2 || ng < 2 || nx < 2 {
            return Err(bad("at least two values of p, G and x are needed".into()));
        }
        let mut chf = vec![f64::NAN; np * ng * nx];
        let mut seen = vec![false; np * ng * nx];
        let find = |v: &[f64], a: f64| v.iter().position(|&b| b == a).unwrap_or(0);
        for r in &rows {
            let k = (find(&p, r[0]) * ng + find(&g, r[1])) * nx + find(&x, r[2]);
            if seen[k] {
                return Err(bad(format!(
                    "the point p = {} kPa, G = {}, x = {} twice",
                    r[0] * 1e-3,
                    r[1],
                    r[2]
                )));
            }
            seen[k] = true;
            chf[k] = r[3];
        }
        if let Some(k) = seen.iter().position(|s| !s) {
            let (ip, ig, ix) = (k / (ng * nx), (k / nx) % ng, k % nx);
            return Err(bad(format!(
                "the point p = {} kPa, G = {}, x = {} is missing",
                p[ip] * 1e-3,
                g[ig],
                x[ix]
            )));
        }
        Ok(Self { p, g, x, chf })
    }

    /// The CHF (W/m²) of a tube of diameter `d` (m) by trilinear
    /// interpolation, the arguments clamped to the table, with the
    /// diameter factor (0.008/D)^½ (D clamped to 3 to 25 mm); NaN where a
    /// corner of the cell has no value.
    pub fn value(&self, p: f64, g: f64, x: f64, d: f64) -> f64 {
        let at = |v: &[f64], a: f64| {
            let a = a.clamp(v[0], v[v.len() - 1]);
            let i = v.partition_point(|&b| b <= a).clamp(1, v.len() - 1) - 1;
            (i, (a - v[i]) / (v[i + 1] - v[i]))
        };
        let (ip, tp) = at(&self.p, p);
        let (ig, tg) = at(&self.g, g);
        let (ix, tx) = at(&self.x, x);
        let (ng, nx) = (self.g.len(), self.x.len());
        let mut s = 0.0;
        for (a, wa) in [(0, 1.0 - tp), (1, tp)] {
            for (b, wb) in [(0, 1.0 - tg), (1, tg)] {
                for (c, wc) in [(0, 1.0 - tx), (1, tx)] {
                    let w = wa * wb * wc;
                    if w != 0.0 {
                        s += w * self.chf[((ip + a) * ng + ig + b) * nx + ix + c];
                    }
                }
            }
        }
        s * (0.008 / d.clamp(0.003, 0.025)).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bowring_is_continuous_at_its_reference_pressure_and_falls_with_quality() {
        // p_R = 1 at 6.895 MPa: both branches give F = 1
        let p1 = 1.0 / 0.145e-6;
        let below = bowring(p1 * (1.0 - 1e-9), 2000.0, 0.01, 1.5e6, 0.1);
        let above = bowring(p1 * (1.0 + 1e-9), 2000.0, 0.01, 1.5e6, 0.1);
        assert!((below / above - 1.0).abs() < 1e-6);
        // at p_R = 1: A = 2.317 h_fg D G/4/(1 + 0.0143 D^0.5 G), C = 0.077 D G/(1 + 0.347 (G/1356)^1.5)
        let (d, g, hfg) = (0.01_f64, 2000.0, 1.5e6);
        let a = 2.317 * hfg * d * g / 4.0 / (1.0 + 0.0143 * d.sqrt() * g);
        let c = 0.077 * d * g / (1.0 + 0.347 * (g / 1356.0_f64).powf(1.5));
        assert!((bowring(p1, g, d, hfg, 0.0) / (a / c) - 1.0).abs() < 1e-9);
        let mut last = f64::INFINITY;
        for i in 0..10 {
            let q = bowring(7e6, 1500.0, 0.012, 1.5e6, -0.2 + 0.1 * i as f64);
            assert!(q < last);
            last = q;
        }
    }

    #[test]
    fn table_reads_and_interpolates() {
        let mut csv = String::from("p_kpa,g,x,chf_kw\n");
        for p in [1000.0, 5000.0] {
            for g in [0.0, 1000.0] {
                for x in [-0.5, 0.5] {
                    // linear in each argument: trilinear interpolation is exact
                    let q = 1000.0 + 0.1 * p + g - 2000.0 * x;
                    csv += &format!("{p},{g},{x},{q}\n");
                }
            }
        }
        let t = LookUpTable::from_csv(&csv).unwrap();
        let q = t.value(3e6, 400.0, 0.1, 0.008);
        let exact = (1000.0 + 300.0 + 400.0 - 200.0) * 1e3;
        assert!((q - exact).abs() < 1e-6);
        // diameter factor
        assert!((t.value(3e6, 400.0, 0.1, 0.016) / q - (0.5_f64).sqrt()).abs() < 1e-12);
        // a missing point is refused
        let short: String = csv.lines().take(8).collect::<Vec<_>>().join("\n");
        assert!(LookUpTable::from_csv(&short).is_err());
    }
}
