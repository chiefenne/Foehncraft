//! Density-field topology search with the Brinkman penalisation (L13
//! stage 1, D226).
//!
//! The cells of a design region carry a fluid fraction γ ∈ [γ_min, 1]
//! (1 fluid, 0 solid), the design variable of Borrvall & Petersson (2003).
//! Its solid fraction in the momentum rows is their interpolation
//!
//! χ(γ) = q (1 − γ) / (q + γ),
//!
//! so the Brinkman sink is μ χ(γ̃) V / K with the permeability K of
//! `[immersed]` (the sink of a fully solid cell); q sets how convex the
//! interpolation is (small q: grey cells are penalised hard). γ̃ is γ after
//! the density filter of Bruns & Tortorelli (2001) and Bourdin (2001) with a
//! linear hat of radius `filter_radius` (no filter with 0, the default).
//!
//! **Objective.** The power the flow takes from its boundaries,
//!
//! Φ = −Σ_b (p_b / ρ + ½ |u_b|²) ṁ_b,
//!
//! the total-pressure flux into the domain; at the steady state it is the
//! dissipation (viscous and Brinkman). The search minimises Φ with the fluid
//! volume Σ γ̃ V of the design region held at `fluid_fraction` of the
//! region's volume.
//!
//! **Gradient.** For Stokes flow the minimum-dissipation principle makes
//! the problem self-adjoint: the velocity minimises ½ a_γ(v, v) over the
//! divergence-free fields with the boundary values, Φ = a_γ(u, u), and the
//! derivative is that of the bilinear form at the solution (Borrvall &
//! Petersson 2003, §4):
//!
//! ∂Φ/∂γ̃_c = (μ / K) χ'(γ̃_c) |u_c|² V_c ≤ 0.
//!
//! No adjoint solve is needed (the adjoint solver L17 is not built). With
//! convection the formula is an approximation that is exact as Re → 0; the
//! gradient check (`tests/immersed.rs`) measures both against central
//! differences of converged runs.
//!
//! **Update.** Optimality criteria (Bendsøe & Sigmund 2003, §1.2): γ ←
//! clamp(γ (−∂Φ/∂γ / (λ ∂V/∂γ))^½) within the move limit, λ by bisection so
//! that the fluid volume is met. The design is updated during the flow
//! iterations ("one-shot", D206): every `every` iterations or when the flow
//! has converged for the current design; after the last update (the largest
//! change of γ below `change`, or `designs` updates) the flow converges on
//! the final design.
//!
//! The run writes `topology.csv` (one row per design update) and
//! `topology.stl`, the surface γ̃ = ½ extracted with marching tetrahedra
//! (`foehncraft_mesh::isosurface`).
//!
//! ```toml
//! [immersed]
//! permeability = 1e-5          # K of a fully solid cell (m²)
//!
//! [immersed.topology]
//! lower = [0.0, 0.0, -1.0]     # the design region (a box)
//! upper = [1.0, 1.0, 1.0]
//! fluid_fraction = 0.25        # the fluid volume allowed in it
//! # interpolation = 0.1        # q
//! # filter_radius = 0.0        # m (0: no filter)
//! # move_limit = 0.2
//! # every = 20                 # flow iterations between design updates
//! # designs = 100              # most design updates
//! # change = 0.01              # stop when no γ changes by more
//! ```

use crate::diag::{self, Diagnostic};
use anyhow::Result;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;

/// The lower bound of the fluid fraction (the multiplicative update cannot
/// leave 0).
pub const GAMMA_MIN: f64 = 1e-3;

/// `[immersed.topology]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopologySection {
    /// Lower corner of the design region (m): cells whose centroid lies in
    /// the box are design cells, the others stay fluid (or the primitive
    /// bodies they are in).
    pub lower: [f64; 3],
    /// Upper corner of the design region (m).
    pub upper: [f64; 3],
    /// The fluid volume allowed, as a fraction of the design region's
    /// volume (0 to 1).
    pub fluid_fraction: f64,
    /// q of the interpolation χ(γ) = q (1 − γ)/(q + γ) (default 0.1).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interpolation: Option<f64>,
    /// Radius of the density filter (m); 0 or absent: no filter.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub filter_radius: Option<f64>,
    /// Largest change of γ per update (default 0.2).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub move_limit: Option<f64>,
    /// Flow iterations between two design updates (default 20).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub every: Option<usize>,
    /// The most design updates (default 100).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub designs: Option<usize>,
    /// The search stops when no γ changes by more than this (default 0.01).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub change: Option<f64>,
    /// The start value of γ in the design region (default
    /// `fluid_fraction`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub start: Option<f64>,
}

fn bad(key: &str, msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::IMM_CONFIG, msg).key(format!("immersed.topology.{key}")))
}

impl TopologySection {
    /// q.
    pub fn q(&self) -> f64 {
        self.interpolation.unwrap_or(0.1)
    }
    fn move_limit(&self) -> f64 {
        self.move_limit.unwrap_or(0.2)
    }
    fn every(&self) -> usize {
        self.every.unwrap_or(20)
    }
    fn designs(&self) -> usize {
        self.designs.unwrap_or(100)
    }
    fn change(&self) -> f64 {
        self.change.unwrap_or(0.01)
    }

    /// Checks the keys (IMM101).
    pub fn check(&self) -> Result<()> {
        if (0..3).any(|i| self.upper[i] <= self.lower[i]) {
            return Err(bad("upper", "`upper` must exceed `lower`".into()));
        }
        let unit = |v: f64| v > 0.0 && v <= 1.0;
        if !unit(self.fluid_fraction) {
            return Err(bad(
                "fluid_fraction",
                format!("fluid_fraction = {}: give 0 < f ≤ 1", self.fluid_fraction),
            ));
        }
        if let Some(s) = self.start
            && !(GAMMA_MIN..=1.0).contains(&s)
        {
            return Err(bad("start", format!("start = {s}: give {GAMMA_MIN} to 1")));
        }
        if !(self.q() > 0.0 && self.q().is_finite()) {
            return Err(bad("interpolation", "interpolation (q) must be > 0".into()));
        }
        if let Some(r) = self.filter_radius
            && !(r >= 0.0 && r.is_finite())
        {
            return Err(bad("filter_radius", "filter_radius must be ≥ 0".into()));
        }
        if !unit(self.move_limit()) {
            return Err(bad("move_limit", "move_limit: give 0 < m ≤ 1".into()));
        }
        if self.every() == 0 || self.designs() == 0 {
            return Err(bad("every", "every and designs must be ≥ 1".into()));
        }
        if !(self.change() > 0.0) {
            return Err(bad("change", "change must be > 0".into()));
        }
        Ok(())
    }
}

/// χ(γ) = q (1 − γ)/(q + γ).
#[inline]
pub fn chi(q: f64, g: f64) -> f64 {
    q * (1.0 - g) / (q + g)
}

/// dχ/dγ = −q (1 + q)/(q + γ)².
#[inline]
pub fn chi_prime(q: f64, g: f64) -> f64 {
    -q * (1.0 + q) / ((q + g) * (q + g))
}

/// The objective Φ = −Σ_b (p_b/ρ + ½|u_b|²) ṁ_b (W) of a state.
pub fn power(st: &crate::simplec::State, mesh: &Mesh, rho: f64) -> f64 {
    let ni = mesh.n_interior_faces();
    let mut s = 0.0;
    for (j, &pb) in st.pb.iter().enumerate() {
        let m = st.mass_flux[ni + j];
        if m == 0.0 {
            continue;
        }
        let u = st.ub.get(j);
        s -= (pb / rho + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2])) * m;
    }
    s
}

/// The density filter: for each design cell its neighbours within the
/// radius and their weights (r − d) V, normalised.
#[derive(Debug, Clone)]
struct Filter {
    start: Vec<usize>,
    index: Vec<usize>,
    weight: Vec<f64>,
}

impl Filter {
    fn identity(n: usize) -> Self {
        Self {
            start: (0..=n).collect(),
            index: (0..n).collect(),
            weight: vec![1.0; n],
        }
    }

    /// The hat filter over the design cells `x` (centroids) with volumes
    /// `v`: neighbours found on a bucket grid of width r, in a fixed order.
    fn hat(x: &[[f64; 3]], v: &[f64], r: f64) -> Self {
        use std::collections::BTreeMap;
        let key = |p: [f64; 3]| -> [i64; 3] { std::array::from_fn(|i| (p[i] / r).floor() as i64) };
        let mut grid: BTreeMap<[i64; 3], Vec<usize>> = BTreeMap::new();
        for (i, &p) in x.iter().enumerate() {
            grid.entry(key(p)).or_default().push(i);
        }
        let mut start = vec![0];
        let mut index = Vec::new();
        let mut weight = Vec::new();
        for &p in x {
            let k = key(p);
            let first = index.len();
            let mut sum = 0.0;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = grid.get(&[k[0] + dx, k[1] + dy, k[2] + dz]) {
                            for &j in list {
                                let d = (0..3).map(|i| (x[j][i] - p[i]).powi(2)).sum::<f64>().sqrt();
                                if d < r {
                                    let w = (r - d) * v[j];
                                    index.push(j);
                                    weight.push(w);
                                    sum += w;
                                }
                            }
                        }
                    }
                }
            }
            for w in &mut weight[first..] {
                *w /= sum;
            }
            start.push(index.len());
        }
        Self {
            start,
            index,
            weight,
        }
    }

    /// γ̃ = W γ.
    fn apply(&self, g: &[f64]) -> Vec<f64> {
        (0..self.start.len() - 1)
            .map(|c| {
                (self.start[c]..self.start[c + 1])
                    .map(|k| self.weight[k] * g[self.index[k]])
                    .sum()
            })
            .collect()
    }

    /// Wᵀ y (the chain rule of the filter), deterministic order.
    fn transpose(&self, y: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; y.len()];
        for (c, &yc) in y.iter().enumerate() {
            for k in self.start[c]..self.start[c + 1] {
                out[self.index[k]] += self.weight[k] * yc;
            }
        }
        out
    }
}

/// One row of `topology.csv`.
#[derive(Debug, Clone, Serialize)]
pub struct DesignRow {
    /// Design update (1, 2, …).
    pub design: usize,
    /// Flow iteration of the update.
    pub iteration: usize,
    /// Φ at the update (W).
    pub objective: f64,
    /// Fluid volume of the design region over its volume (after the
    /// update).
    pub fluid_fraction: f64,
    /// Largest change of γ.
    pub change: f64,
    /// Grey level 4 γ̃ (1 − γ̃) averaged over the design volume (0: black
    /// and white).
    pub grey: f64,
}

/// The summary in `result.json`.
#[derive(Debug, Clone, Serialize)]
pub struct TopologyReport {
    /// Design updates done.
    pub designs: usize,
    /// The search stopped on `change` (not on `designs`).
    pub converged: bool,
    /// Φ of the final state (W).
    pub objective: f64,
    /// Φ at the first update (W).
    pub first_objective: f64,
    /// Fluid volume fraction of the design region.
    pub fluid_fraction: f64,
    /// Grey level of the final design.
    pub grey: f64,
    /// Σ_c ∂Φ/∂γ_c at the first update (W): the derivative of Φ when every
    /// γ of the design region moves together (the gradient check).
    pub first_gradient: f64,
    /// Design cells.
    pub cells: usize,
    /// Triangles of `topology.stl`.
    pub triangles: usize,
    /// Area of `topology.stl` (m²).
    pub surface_area: f64,
}

/// A running search.
#[derive(Debug, Clone)]
pub struct Search {
    sec: TopologySection,
    /// The design cells (local indices).
    pub cells: Vec<usize>,
    gamma: Vec<f64>,
    filter: Filter,
    vol: Vec<f64>,
    region: f64,
    /// χ of the primitive bodies per cell.
    base: Vec<f64>,
    mu_k: f64,
    rho: f64,
    last: usize,
    /// Design updates done.
    pub done: usize,
    /// No more updates.
    pub finished: bool,
    converged: bool,
    first_gradient: f64,
    /// The rows of `topology.csv`.
    pub history: Vec<DesignRow>,
}

impl Search {
    /// The search of `case` on the mesh; `base` is χ of the primitive
    /// bodies per cell. None without `[immersed.topology]`.
    pub fn new(
        case: &crate::case::FlowCase,
        mesh: &Mesh,
        geo: &Geometry,
        base: Vec<f64>,
    ) -> Result<Option<Self>> {
        let Some(sec) = case.immersed.topology.clone() else {
            return Ok(None);
        };
        sec.check()?;
        let inside = |x: [f64; 3]| (0..3).all(|i| x[i] >= sec.lower[i] && x[i] <= sec.upper[i]);
        let cells: Vec<usize> = (0..mesh.n_cells)
            .filter(|&c| inside(geo.cell_centroid[c]) && base[c] < 1.0)
            .collect();
        if cells.is_empty() {
            return Err(bad(
                "lower",
                "the design region holds no cell centroid (outside the mesh?)".into(),
            ));
        }
        let vol: Vec<f64> = cells.iter().map(|&c| geo.cell_volume[c]).collect();
        let region: f64 = vol.iter().sum();
        let x: Vec<[f64; 3]> = cells.iter().map(|&c| geo.cell_centroid[c]).collect();
        let filter = match sec.filter_radius {
            Some(r) if r > 0.0 => Filter::hat(&x, &vol, r),
            _ => Filter::identity(cells.len()),
        };
        let g0 = sec.start.unwrap_or(sec.fluid_fraction);
        let k = case.immersed.permeability.unwrap_or(f64::INFINITY);
        Ok(Some(Self {
            gamma: vec![g0; cells.len()],
            cells,
            filter,
            vol,
            region,
            base,
            mu_k: case.mu() / k,
            rho: case.rho(),
            last: 0,
            done: 0,
            finished: false,
            converged: false,
            first_gradient: 0.0,
            history: Vec::new(),
            sec,
        }))
    }

    /// γ̃ of the design cells.
    pub fn filtered(&self) -> Vec<f64> {
        self.filter.apply(&self.gamma)
    }

    /// χ per cell in the momentum rows: the bodies plus χ(γ̃) in the
    /// design region (capped at 1).
    pub fn chi(&self) -> Vec<f64> {
        let mut out = self.base.clone();
        let q = self.sec.q();
        for (&c, &g) in self.cells.iter().zip(&self.filtered()) {
            out[c] = (out[c] + chi(q, g)).min(1.0);
        }
        out
    }

    /// The solid fraction for the output: the bodies plus 1 − γ̃.
    pub fn solid(&self) -> Vec<f64> {
        let mut out = self.base.clone();
        for (&c, &g) in self.cells.iter().zip(&self.filtered()) {
            out[c] = (out[c] + 1.0 - g).min(1.0);
        }
        out
    }

    /// The sink μ χ V / K per cell.
    pub fn sink(&self, geo: &Geometry) -> Vec<f64> {
        self.chi()
            .iter()
            .zip(&geo.cell_volume)
            .map(|(x, v)| self.mu_k * x * v)
            .collect()
    }

    /// The design part of χ per cell (0 outside the design region), for
    /// the force on the design.
    pub fn design_chi(&self, nc: usize) -> Vec<f64> {
        let mut out = vec![0.0; nc];
        let q = self.sec.q();
        for (&c, &g) in self.cells.iter().zip(&self.filtered()) {
            out[c] = chi(q, g);
        }
        out
    }

    /// True when a design update is due at flow iteration `it` (`flow_done`:
    /// the flow has converged for the current design).
    pub fn due(&self, it: usize, flow_done: bool) -> bool {
        !self.finished && (flow_done || it >= self.last + self.sec.every())
    }

    /// ∂Φ/∂γ of the design variables at the velocity `u` (self-adjoint
    /// form, exact for Stokes flow), through the filter.
    pub fn gradient(&self, u: &foehncraft_fv::disc::Soa3) -> Vec<f64> {
        let q = self.sec.q();
        let gt = self.filtered();
        let dt: Vec<f64> = self
            .cells
            .iter()
            .zip(&gt)
            .zip(&self.vol)
            .map(|((&c, &g), &v)| {
                let w = u.get(c);
                self.mu_k * chi_prime(q, g) * (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]) * v
            })
            .collect();
        self.filter.transpose(&dt)
    }

    /// Fluid volume of the design region over its volume.
    pub fn fluid_fraction(&self) -> f64 {
        self.filtered()
            .iter()
            .zip(&self.vol)
            .map(|(g, v)| g * v)
            .sum::<f64>()
            / self.region
    }

    fn grey(&self) -> f64 {
        self.filtered()
            .iter()
            .zip(&self.vol)
            .map(|(g, v)| 4.0 * g * (1.0 - g) * v)
            .sum::<f64>()
            / self.region
    }

    /// The design update at flow iteration `it` from the state `st`;
    /// returns the log line.
    pub fn update(&mut self, st: &crate::simplec::State, mesh: &Mesh, it: usize) -> String {
        let phi = power(st, mesh, self.rho);
        let df = self.gradient(&st.u);
        if self.done == 0 {
            self.first_gradient = df.iter().sum();
        }
        // ∂V/∂γ through the filter
        let dv = self.filter.transpose(&self.vol);
        let m = self.sec.move_limit();
        let target = self.sec.fluid_fraction * self.region;
        let step = |lambda: f64| -> Vec<f64> {
            self.gamma
                .iter()
                .zip(&df)
                .zip(&dv)
                .map(|((&g, &d), &v)| {
                    let b = (-d).max(0.0) / (lambda * v);
                    let g1 = g * b.sqrt();
                    g1.clamp((g - m).max(GAMMA_MIN), (g + m).min(1.0))
                })
                .collect()
        };
        let volume = |g: &[f64]| -> f64 {
            self.filter
                .apply(g)
                .iter()
                .zip(&self.vol)
                .map(|(a, v)| a * v)
                .sum()
        };
        // bisection on λ in log space: the volume falls as λ grows
        let scale = df.iter().map(|d| d.abs()).fold(0.0f64, f64::max).max(f64::MIN_POSITIVE);
        let (mut lo, mut hi) = (scale * 1e-12, scale * 1e12);
        for _ in 0..200 {
            let mid = (lo * hi).sqrt();
            if volume(&step(mid)) > target {
                lo = mid;
            } else {
                hi = mid;
            }
            if hi / lo < 1.0 + 1e-12 {
                break;
            }
        }
        let new = step((lo * hi).sqrt());
        let change = new
            .iter()
            .zip(&self.gamma)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        self.gamma = new;
        self.done += 1;
        self.last = it;
        self.converged = change < self.sec.change();
        self.finished = self.converged || self.done >= self.sec.designs();
        let row = DesignRow {
            design: self.done,
            iteration: it,
            objective: phi,
            fluid_fraction: self.fluid_fraction(),
            change,
            grey: self.grey(),
        };
        let line = format!(
            "topology design {}: Φ = {:.6e} W, fluid {:.4}, change {:.3e}, grey {:.3}{}",
            row.design,
            row.objective,
            row.fluid_fraction,
            row.change,
            row.grey,
            if self.finished {
                if self.converged {
                    " (converged; the flow now converges on this design)"
                } else {
                    " (the most design updates; the flow now converges on this design)"
                }
            } else {
                ""
            }
        );
        self.history.push(row);
        line
    }

    /// `topology.csv`.
    pub fn csv(&self) -> String {
        let mut s = String::from(
            "# topology search (L13 stage 1, D226): one row per design update\n\
             # units: design -, iteration -, objective W, fluid_fraction -, change -, grey -\n\
             design,iteration,objective,fluid_fraction,change,grey\n",
        );
        for r in &self.history {
            let _ = writeln!(
                s,
                "{},{},{:e},{:e},{:e},{:e}",
                r.design, r.iteration, r.objective, r.fluid_fraction, r.change, r.grey
            );
        }
        s
    }

    /// The summary, with Φ of the final state and the surface.
    pub fn report(&self, objective: f64, surface: &foehncraft_mesh::isosurface::Surface) -> TopologyReport {
        TopologyReport {
            designs: self.done,
            converged: self.converged,
            objective,
            first_objective: self.history.first().map_or(objective, |r| r.objective),
            fluid_fraction: self.fluid_fraction(),
            grey: self.grey(),
            first_gradient: self.first_gradient,
            cells: self.cells.len(),
            triangles: surface.triangles.len(),
            surface_area: surface.area(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolation_and_its_derivative() {
        let q = 0.1;
        assert!((chi(q, 1.0)).abs() < 1e-15);
        assert!((chi(q, 0.0) - 1.0).abs() < 1e-15);
        for g in [0.01, 0.3, 0.7, 0.99] {
            let h = 1e-6;
            let fd = (chi(q, g + h) - chi(q, g - h)) / (2.0 * h);
            assert!((fd - chi_prime(q, g)).abs() < 1e-6 * fd.abs(), "{g}");
        }
    }

    #[test]
    fn filter_keeps_constants_and_its_transpose_is_exact() {
        let x: Vec<[f64; 3]> = (0..30)
            .map(|i| [(i % 6) as f64 * 0.1, (i / 6) as f64 * 0.1, 0.0])
            .collect();
        let v: Vec<f64> = (0..30).map(|i| 1.0 + 0.01 * i as f64).collect();
        let f = Filter::hat(&x, &v, 0.25);
        let one = f.apply(&vec![0.7; 30]);
        assert!(one.iter().all(|a| (a - 0.7).abs() < 1e-14));
        // ⟨W a, b⟩ = ⟨a, Wᵀ b⟩
        let a: Vec<f64> = (0..30).map(|i| ((i * 7) % 11) as f64).collect();
        let b: Vec<f64> = (0..30).map(|i| ((i * 5) % 13) as f64).collect();
        let l: f64 = f.apply(&a).iter().zip(&b).map(|(p, q)| p * q).sum();
        let r: f64 = a.iter().zip(&f.transpose(&b)).map(|(p, q)| p * q).sum();
        assert!((l - r).abs() < 1e-12 * l.abs());
    }
}
