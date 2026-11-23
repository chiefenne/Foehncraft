//! Monitors of a run (D43): point probes, force coefficients, CSV histories
//! that survive a restart, and the evaluation of the cylinder benchmark of
//! Schäfer & Turek (1996).

use crate::analytic::{Benchmark, Range};
use crate::case::{BoundarySummary, FlowCase, locate, value_at};
use crate::simplec::State;
use anyhow::Result;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::geometry::dot;
use foehncraft_mesh::locate::Locator;
use serde::Serialize;
use std::path::{Path, PathBuf};

/// Point probes with their cells.
pub struct Probes {
    /// Names.
    pub names: Vec<String>,
    /// Points.
    pub points: Vec<[f64; 3]>,
    cells: Vec<usize>,
    /// The operating pressure added to the pressure of each probe: absolute
    /// pressure values (D239); empty with gauge values.
    offsets: Vec<f64>,
}

impl Probes {
    /// The probes of the case plus `extra` (name, point).
    pub fn new(
        case: &FlowCase,
        extra: &[(String, [f64; 3])],
        loc: &Locator,
        geo: &Geometry,
    ) -> Self {
        let mut names: Vec<String> = case.output.probes.iter().map(|p| p.name.clone()).collect();
        let mut points: Vec<[f64; 3]> = case.output.probes.iter().map(|p| p.point).collect();
        for (n, x) in extra {
            names.push(n.clone());
            points.push(*x);
        }
        let cells: Vec<usize> = points.iter().map(|&x| locate(loc, geo, x)).collect();
        let offsets = case
            .pressure_frame
            .as_ref()
            .and_then(|f| f.cell_offsets(case, loc.mesh()))
            .map(|o| {
                cells
                    .iter()
                    .map(|&c| o.get(c).copied().unwrap_or(0.0))
                    .collect()
            })
            .unwrap_or_default();
        Self {
            names,
            points,
            cells,
            offsets,
        }
    }

    /// The probes found again in the cells of `mesh` (a sliding mesh by
    /// imprinting, D182: the cells of the rotor are somewhere else in every
    /// step).
    pub fn relocated(&self, mesh: &foehncraft_mesh::Mesh, geo: &Geometry) -> Self {
        let loc = Locator::new(mesh, geo);
        Self {
            names: self.names.clone(),
            points: self.points.clone(),
            cells: self.points.iter().map(|&x| locate(&loc, geo, x)).collect(),
            offsets: self.offsets.clone(),
        }
    }

    /// Velocity and pressure at every probe (reconstructed with the cell
    /// gradients).
    pub fn values(&self, sol: &State, geo: &Geometry) -> Vec<([f64; 3], f64)> {
        self.cells
            .iter()
            .zip(&self.points)
            .enumerate()
            .map(|(i, (&c, &x))| {
                let (u, p) = value_at(sol, geo, c, x);
                match self.offsets.get(i) {
                    Some(o) => (u, p + o),
                    None => (u, p),
                }
            })
            .collect()
    }

    /// The temperature of the energy equation at every probe (the value of
    /// its cell, stage 7d6); NaN without the energy equation.
    pub fn temperatures(&self, sol: &State) -> Vec<f64> {
        let t = sol.scalars.iter().find(|s| s.name == crate::energy::FIELD);
        self.cells
            .iter()
            .map(|&c| t.map_or(f64::NAN, |t| t.phi[c]))
            .collect()
    }

    /// The velocity relative to the turning zone at every probe (D183),
    /// from the values `vals` of [`Probes::values`].
    pub fn relative(
        &self,
        frames: &crate::derived::Frames,
        mesh: &foehncraft_mesh::Mesh,
        vals: &[([f64; 3], f64)],
    ) -> Vec<[f64; 3]> {
        self.cells
            .iter()
            .zip(&self.points)
            .zip(vals)
            .map(|((&c, &x), (u, _))| frames.at(mesh, c, x, *u).w)
            .collect()
    }

    /// Value of the probe `name`.
    pub fn index(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }
}

/// Drag and lift coefficients of the case's `[output.coefficients]` from
/// the boundary forces.
pub fn coefficients(case: &FlowCase, boundaries: &[BoundarySummary]) -> Option<[f64; 2]> {
    let c = case.output.coefficients.as_ref()?;
    let mut f = [0.0; 3];
    for p in boundaries.iter().filter(|p| c.boundaries.contains(&p.name)) {
        for (fi, pi) in f.iter_mut().zip(p.force) {
            *fi += pi;
        }
    }
    let q = 0.5 * case.rho() * c.velocity * c.velocity * c.area;
    Some([dot(f, c.drag_direction) / q, dot(f, c.lift_direction) / q])
}

/// A CSV history that continues after a restart: rows of later steps than
/// the restart are dropped, new rows are appended.
pub struct History {
    path: PathBuf,
    buf: String,
}

impl History {
    /// A history `name` in `dir` with the comment lines `preamble` (D75)
    /// and the column `header`; with `continue_after` the existing rows up
    /// to that step (first column) are kept.
    pub fn new(
        dir: &Path,
        name: &str,
        preamble: &str,
        header: &str,
        continue_after: Option<usize>,
    ) -> Result<Self> {
        let path = dir.join(name);
        let mut text = format!("{preamble}{header}\n");
        if let Some(n) = continue_after
            && let Ok(old) = std::fs::read_to_string(&path)
        {
            let mut lines = old.lines().skip_while(|l| l.starts_with('#'));
            if lines.next() == Some(header) {
                for l in lines {
                    let step = l.split(',').next().and_then(|v| v.parse::<usize>().ok());
                    if step.is_some_and(|s| s <= n) {
                        text.push_str(l);
                        text.push('\n');
                    }
                }
            }
        }
        std::fs::write(&path, &text)?;
        Ok(Self {
            path,
            buf: String::new(),
        })
    }

    /// Adds a row (written on [`flush`](Self::flush)).
    pub fn row(&mut self, r: &str) {
        self.buf.push_str(r);
        self.buf.push('\n');
    }

    /// Appends the buffered rows to the file.
    pub fn flush(&mut self) -> Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        use std::io::Write as _;
        let mut f = std::fs::OpenOptions::new().append(true).open(&self.path)?;
        f.write_all(self.buf.as_bytes())?;
        self.buf.clear();
        Ok(())
    }
}

/// One quantity of the benchmark.
#[derive(Debug, Clone, Serialize)]
pub struct BenchmarkValue {
    /// The reference range.
    pub range: Range,
    /// The computed value.
    pub value: f64,
    /// Distance from the range relative to its midpoint (0 inside).
    pub deviation: f64,
}

/// Comparison with the benchmark.
#[derive(Debug, Clone, Serialize)]
pub struct BenchmarkResult {
    /// The benchmark.
    pub benchmark: String,
    /// Computed values against the reference ranges.
    pub values: Vec<BenchmarkValue>,
    /// Periods of C_L used for the Strouhal number (periodic case).
    pub periods: usize,
}

impl BenchmarkResult {
    /// Largest relative deviation.
    pub fn max_deviation(&self) -> f64 {
        self.values.iter().map(|v| v.deviation).fold(0.0, f64::max)
    }
}

fn result(b: Benchmark, values: &[f64], periods: usize) -> BenchmarkResult {
    BenchmarkResult {
        benchmark: b.label().into(),
        values: b
            .ranges()
            .into_iter()
            .zip(values)
            .map(|(range, &value)| BenchmarkValue {
                deviation: range.deviation(value),
                range,
                value,
            })
            .collect(),
        periods,
    }
}

/// The steady benchmark (2D-1) from the final C_D, C_L and Δp.
pub fn steady_benchmark(b: Benchmark, cd: f64, cl: f64, dp: f64) -> BenchmarkResult {
    result(b, &[cd, cl, dp], 0)
}

/// Times of the upward zero crossings of `y − mean` (linear
/// interpolation).
fn up_crossings(t: &[f64], y: &[f64]) -> Vec<f64> {
    if y.is_empty() {
        return Vec::new();
    }
    let mean = y.iter().sum::<f64>() / y.len() as f64;
    let mut out = Vec::new();
    for i in 1..y.len() {
        let (a, b) = (y[i - 1] - mean, y[i] - mean);
        if a < 0.0 && b >= 0.0 {
            out.push(t[i - 1] + (t[i] - t[i - 1]) * (-a) / (b - a));
        }
    }
    out
}

fn interp(t: &[f64], y: &[f64], x: f64) -> f64 {
    match t.iter().position(|&ti| ti >= x) {
        Some(0) => y[0],
        Some(i) => y[i - 1] + (y[i] - y[i - 1]) * (x - t[i - 1]) / (t[i] - t[i - 1]),
        None => *y.last().unwrap_or(&0.0),
    }
}

/// The periodic benchmark (2D-2) from the histories of C_D, C_L and Δp in
/// the window t ≥ `start`: C_D,max, C_L,max, St = D f / Ū from the mean
/// period of the last (up to five) C_L periods, and Δp at t₀ + 1/(2f) with
/// t₀ the time of the last C_L maximum that leaves room for half a period.
/// None while the window has less than two full periods.
pub fn periodic_benchmark(
    b: Benchmark,
    t: &[f64],
    cd: &[f64],
    cl: &[f64],
    dp: &[f64],
    start: f64,
) -> Option<BenchmarkResult> {
    let i0 = t.iter().position(|&x| x >= start)?;
    let (t, cd, cl, dp) = (&t[i0..], &cd[i0..], &cl[i0..], &dp[i0..]);
    let ups = up_crossings(t, cl);
    if ups.len() < 3 {
        return None;
    }
    let k = (ups.len() - 1).min(5);
    let period = (ups[ups.len() - 1] - ups[ups.len() - 1 - k]) / k as f64;
    let f = 1.0 / period;
    let st = Benchmark::DIAMETER * f / b.mean_velocity();
    let cd_max = cd.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let cl_max = cl.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    // last C_L maximum with half a period after it inside the window
    let t_end = *t.last()?;
    let mut t0 = None;
    for i in 1..cl.len().saturating_sub(1) {
        if cl[i] >= cl[i - 1] && cl[i] > cl[i + 1] && t[i] + 0.5 * period <= t_end {
            // parabola through the three samples for the time of the maximum
            let (y0, y1, y2) = (cl[i - 1], cl[i], cl[i + 1]);
            let den = y0 - 2.0 * y1 + y2;
            let h = 0.5 * (t[i + 1] - t[i - 1]);
            let off = if den != 0.0 {
                0.5 * (y0 - y2) / den * h
            } else {
                0.0
            };
            t0 = Some(t[i] + off);
        }
    }
    let dp_half = interp(t, dp, t0? + 0.5 * period);
    Some(result(b, &[cd_max, cl_max, st, dp_half], k))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthetic histories with known values: C_L = sin 2πft with St = 0.3
    /// (f = 3 Hz for D = 0.1 m, Ū = 1 m/s), C_D = 3.23 + 0.005 sin 4πft,
    /// Δp = 2.48 + 0.02 sin 2πft (so Δp is 2.46 half a period after the
    /// lift maximum).
    #[test]
    fn periodic_benchmark_of_known_signals() {
        let b = Benchmark::SchaferTurek2d2;
        let f = 0.3 * b.mean_velocity() / Benchmark::DIAMETER;
        let w = std::f64::consts::TAU * f;
        let dt = 1e-4;
        let t: Vec<f64> = (0..20_000).map(|i| 0.0123 + i as f64 * dt).collect();
        let cl: Vec<f64> = t.iter().map(|&t| (w * t).sin()).collect();
        let cd: Vec<f64> = t
            .iter()
            .map(|&t| 3.23 + 0.005 * (2.0 * w * t).sin())
            .collect();
        let dp: Vec<f64> = t.iter().map(|&t| 2.48 + 0.02 * (w * t).sin()).collect();
        assert!(
            periodic_benchmark(b, &t, &cd, &cl, &dp, 1.9).is_none(),
            "no periods left"
        );
        let r = periodic_benchmark(b, &t, &cd, &cl, &dp, 0.2).unwrap();
        let v: Vec<f64> = r.values.iter().map(|v| v.value).collect();
        assert!((v[0] - 3.235).abs() < 1e-8, "C_D,max {}", v[0]);
        assert!((v[1] - 1.0).abs() < 1e-6, "C_L,max {}", v[1]);
        assert!((v[2] - 0.3).abs() < 1e-6, "St {}", v[2]);
        assert!((v[3] - 2.46).abs() < 1e-6, "Δp {}", v[3]);
        assert_eq!(r.periods, 5);
        assert_eq!(r.max_deviation(), 0.0);
    }
}
