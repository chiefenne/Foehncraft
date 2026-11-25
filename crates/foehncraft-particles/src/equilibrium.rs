//! The phase equilibrium at the droplet surface at high pressure, tabulated
//! over the droplet temperature (ECN Spray A, stage SA1, D311, D312, D314).
//!
//! The surface state of a fuel droplet depends only on its temperature T,
//! the thermodynamic pressure p₀ (one number per domain, D248) and the
//! ratios of the ambient gases in the vapour (D311), not on the far-field
//! fuel fraction. So it is solved once per p₀ on a temperature grid with
//! the Peng-Robinson kernel of `foehncraft_kernels::realfluid`, by
//! continuation in T from the cold end (each node starts from the K values
//! of the node below), up to the temperature where the two-phase solution
//! ends: the critical mixing state (above the fuel's critical pressure) or
//! the boiling point (below it). A parcel then reads the table with linear
//! interpolation; the table is a pure function of (mixture, ratios, p₀), so
//! every pass, thread and restart sees the same bits (D247).
//!
//! Each node holds the fuel mole fraction of the vapour at the surface, the
//! dissolved gas in the surface liquid and the ratio of the latent heat at
//! p₀ to the latent heat at saturation of the pure fuel, both from the
//! partial molar departure enthalpies of the same equation of state:
//! L_PR(T, p₀) = h̄_f,v − h̄_f,l (surface vapour minus surface liquid),
//! L_PR,sat(T) = h_v − h_l of the pure fuel at its Peng-Robinson saturation
//! pressure. The evaporation uses L(T, p₀) = L_lib(T) L_PR/L_PR,sat, so the
//! measured latent heat of the liquid library is kept at low pressure and
//! the heat falls to zero at the critical mixing state (own derivation,
//! D312).

use anyhow::{Result, bail};
use foehncraft_kernels::realfluid::{self as rf, MAX, Mixture, Root, SurfaceState};

/// Where the two-phase solution of the table ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TableEnd {
    /// The K values collapse: the critical mixing state (D314).
    Critical,
    /// The fuel's vapour fraction reaches one: the boiling point at p₀.
    Boiling,
    /// The upper end of the requested range.
    Range,
}

/// The surface state at one droplet temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfacePoint {
    /// False at or above the end of the two-phase solution.
    pub two_phase: bool,
    /// Fuel mole fraction of the vapour at the surface.
    pub y_fuel: f64,
    /// Mole fraction of dissolved gas in the surface liquid.
    pub dissolved: f64,
    /// L(T, p₀)/L_sat(T) of the equation of state (D312).
    pub latent_ratio: f64,
}

/// The surface states over the droplet temperature at one p₀.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceTable {
    /// The thermodynamic pressure (Pa).
    pub p0: f64,
    /// Temperature of the first node (K).
    pub t0: f64,
    /// Node spacing (K); the last node sits at [`SurfaceTable::t_end`].
    pub dt: f64,
    /// Node temperatures (K).
    pub t: Vec<f64>,
    /// Fuel mole fraction of the vapour per node.
    pub y_fuel: Vec<f64>,
    /// Dissolved gas mole fraction per node.
    pub dissolved: Vec<f64>,
    /// Latent heat ratio per node.
    pub latent_ratio: Vec<f64>,
    /// The last two-phase temperature (K).
    pub t_end: f64,
    /// Why the table ends there.
    pub end: TableEnd,
}

/// The latent heat ratio L_PR(T, p)/L_PR,sat(T) at a two-phase surface
/// state (module documentation); 1 when the pure fuel has no saturation
/// state at `t` (at or above its critical temperature).
fn latent_ratio(mix: &Mixture, s: &rf::Surface, t: f64, p: f64) -> f64 {
    let hv = mix.partial_h_dep(&s.y, t, p, Root::Vapour)[0];
    let hl = mix.partial_h_dep(&s.x, t, p, Root::Liquid)[0];
    let pure = Mixture::new(&[mix.c[0]]);
    let Some(ps) = pure.saturation_pressure(0, t) else {
        return 1.0;
    };
    let mut x = [0.0; MAX];
    x[0] = 1.0;
    let sat = pure.phase(&x, t, ps, Root::Vapour).h_dep - pure.phase(&x, t, ps, Root::Liquid).h_dep;
    if sat <= 0.0 {
        return 1.0;
    }
    ((hv - hl) / sat).max(0.0)
}

impl SurfaceTable {
    /// Builds the table of `mix` (component 0 the fuel) with the gas ratios
    /// `r` in the vapour at `p0` (Pa), from `t_lo` up to at most `t_hi`
    /// (K) with the node spacing `dt` (K).
    pub fn build(
        mix: &Mixture,
        r: &[f64; MAX],
        p0: f64,
        t_lo: f64,
        t_hi: f64,
        dt: f64,
    ) -> Result<Self> {
        if !(p0 > 0.0 && t_lo > 0.0 && t_hi > t_lo && dt > 0.0) {
            bail!("surface table: p0 {p0} Pa, range {t_lo} to {t_hi} K, step {dt} K");
        }
        let mut tab = Self {
            p0,
            t0: t_lo,
            dt,
            t: Vec::new(),
            y_fuel: Vec::new(),
            dissolved: Vec::new(),
            latent_ratio: Vec::new(),
            t_end: t_hi,
            end: TableEnd::Range,
        };
        let mut start: Option<[f64; MAX]> = None;
        let mut k = 0usize;
        loop {
            let t = t_lo + k as f64 * dt;
            if t > t_hi {
                break;
            }
            match rf::surface(mix, t, p0, r, start.as_ref()) {
                SurfaceState::TwoPhase(s) => {
                    tab.push(mix, &s, t, p0);
                    start = Some(s.ln_k);
                    k += 1;
                }
                other => {
                    if k == 0 {
                        bail!(
                            "surface table: no two-phase state at the coldest droplet temperature \
                             {t_lo} K and {p0} Pa ({other:?})"
                        );
                    }
                    // the end lies between the last node and t: bisect,
                    // each try from the last two-phase K values
                    let (mut lo, mut hi) = (t - dt, t);
                    let mut last = start.unwrap_or_default();
                    let mut kind = other;
                    for _ in 0..40 {
                        let mid = 0.5 * (lo + hi);
                        match rf::surface(mix, mid, p0, r, Some(&last)) {
                            SurfaceState::TwoPhase(s) => {
                                lo = mid;
                                last = s.ln_k;
                            }
                            o => {
                                hi = mid;
                                kind = o;
                            }
                        }
                    }
                    if lo > *tab.t.last().unwrap_or(&t_lo)
                        && let SurfaceState::TwoPhase(s) = rf::surface(mix, lo, p0, r, Some(&last))
                    {
                        tab.push(mix, &s, lo, p0);
                    }
                    tab.t_end = *tab.t.last().unwrap_or(&lo);
                    tab.end = match kind {
                        SurfaceState::Boiling => TableEnd::Boiling,
                        _ => TableEnd::Critical,
                    };
                    return Ok(tab);
                }
            }
        }
        tab.t_end = *tab.t.last().unwrap_or(&t_hi);
        Ok(tab)
    }

    fn push(&mut self, mix: &Mixture, s: &rf::Surface, t: f64, p: f64) {
        let mut dis = 0.0;
        for v in s.x.iter().take(mix.n).skip(1) {
            dis += v;
        }
        self.t.push(t);
        self.y_fuel.push(s.y[0]);
        self.dissolved.push(dis);
        self.latent_ratio.push(latent_ratio(mix, s, t, p));
    }

    /// The surface state at the droplet temperature `t` (K): linear
    /// between the nodes, the first node below the table, and the end state
    /// (with `two_phase` false and no latent heat) at and above `t_end`.
    pub fn at(&self, t: f64) -> SurfacePoint {
        let n = self.t.len();
        let last = n - 1;
        if t >= self.t_end && self.end != TableEnd::Range {
            return SurfacePoint {
                two_phase: false,
                y_fuel: self.y_fuel[last],
                dissolved: self.dissolved[last],
                latent_ratio: 0.0,
            };
        }
        if t <= self.t[0] || n == 1 {
            return self.node(0);
        }
        if t >= self.t[last] {
            return self.node(last);
        }
        // uniform nodes up to the last, which may sit closer
        let mut i = ((t - self.t0) / self.dt) as usize;
        if i >= last {
            i = last - 1;
        }
        while i + 1 < last && self.t[i + 1] < t {
            i += 1;
        }
        let w = (t - self.t[i]) / (self.t[i + 1] - self.t[i]);
        let lerp = |v: &[f64]| v[i] + w * (v[i + 1] - v[i]);
        SurfacePoint {
            two_phase: true,
            y_fuel: lerp(&self.y_fuel),
            dissolved: lerp(&self.dissolved),
            latent_ratio: lerp(&self.latent_ratio),
        }
    }

    fn node(&self, i: usize) -> SurfacePoint {
        SurfacePoint {
            two_phase: true,
            y_fuel: self.y_fuel[i],
            dissolved: self.dissolved[i],
            latent_ratio: self.latent_ratio[i],
        }
    }
}

/// Fits the Mathias-Copeman coefficients of `c` so that the Peng-Robinson
/// saturation pressure follows `psat` (Pa, a function of T in K, from
/// measured data) at the temperatures `ts` (K, below T_c), by Gauss-Newton
/// on ln p with difference derivatives, from the κ(ω) of Peng and Robinson.
/// Returns the coefficients and the largest relative deviation left.
pub fn fit_mathias_copeman(
    c: rf::Component,
    psat: &dyn Fn(f64) -> f64,
    ts: &[f64],
) -> Result<([f64; 3], f64)> {
    if ts.len() < 3 {
        bail!("Mathias-Copeman fit: at least three temperatures");
    }
    let mut mc = [c.kappa(), 0.0, 0.0];
    let resid = |mc: [f64; 3]| -> Option<Vec<f64>> {
        let mut cc = c;
        cc.mc = mc;
        let m = Mixture::new(&[cc]);
        ts.iter()
            .map(|&t| m.saturation_pressure(0, t).map(|p| (p / psat(t)).ln()))
            .collect()
    };
    for _ in 0..30 {
        let Some(f) = resid(mc) else {
            bail!("Mathias-Copeman fit: no saturation state at the coefficients {mc:?}");
        };
        // J (n × 3)
        let mut jac = vec![[0.0; 3]; f.len()];
        for k in 0..3 {
            let h = 1e-6;
            let mut m2 = mc;
            m2[k] += h;
            let Some(fp) = resid(m2) else {
                bail!("Mathias-Copeman fit: derivative failed");
            };
            for (row, (a, b)) in jac.iter_mut().zip(fp.iter().zip(f.iter())) {
                row[k] = (a - b) / h;
            }
        }
        // normal equations (3 × 3)
        let mut a = [[0.0; 3]; 3];
        let mut g = [0.0; 3];
        for (row, fi) in jac.iter().zip(f.iter()) {
            for i in 0..3 {
                g[i] -= row[i] * fi;
                for j in 0..3 {
                    a[i][j] += row[i] * row[j];
                }
            }
        }
        let d = solve3(a, g).ok_or_else(|| anyhow::anyhow!("Mathias-Copeman fit: singular"))?;
        let mut step: f64 = 0.0;
        for k in 0..3 {
            mc[k] += d[k];
            step = step.max(d[k].abs());
        }
        if step < 1e-12 {
            break;
        }
    }
    let f = resid(mc).ok_or_else(|| anyhow::anyhow!("Mathias-Copeman fit: final state"))?;
    let worst = f.iter().fold(0.0f64, |m, v| m.max(v.exp_m1().abs()));
    Ok((mc, worst))
}

#[allow(clippy::needless_range_loop)]
fn solve3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> Option<[f64; 3]> {
    for c in 0..3 {
        let mut p = c;
        for r in c + 1..3 {
            if a[r][c].abs() > a[p][c].abs() {
                p = r;
            }
        }
        if a[p][c] == 0.0 {
            return None;
        }
        a.swap(c, p);
        b.swap(c, p);
        for r in c + 1..3 {
            let m = a[r][c] / a[c][c];
            for k in c..3 {
                a[r][k] -= m * a[c][k];
            }
            b[r] -= m * b[c];
        }
    }
    let mut x = [0.0; 3];
    for c in (0..3).rev() {
        let mut s = b[c];
        for k in c + 1..3 {
            s -= a[c][k] * x[k];
        }
        x[c] = s / a[c][c];
    }
    Some(x)
}
