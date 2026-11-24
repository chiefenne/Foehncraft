//! Datums of the sprays: the liquid through surfaces and in regions, and
//! the droplet statistics (D262 to D265).
//!
//! A datum seen by the parcels is either a **surface** that parcels cross
//! or a **volume** that holds parcels:
//!
//! - [`Surface::Plane`]: a plane through a point with a unit normal,
//!   optionally bounded by a disc or a rectangle. A sub-step from x₀ to x₁
//!   crosses it when the signed distances s(x) = (x − p)·n at both ends
//!   have different signs (s < 0 is one side, s ≥ 0 the other), at the
//!   chord parameter s₀/(s₀ − s₁), inside the bound. Each end point of a
//!   sub-step is the start of the next, so a parcel that ends a sub-step on
//!   the plane is counted once: watertight, no double counts.
//! - [`Surface::Cells`]: a set of mesh faces given as a side per cell
//!   (+1, −1, 0 for cells not next to it): the faces between two cell zones,
//!   or a plane moved to the faces of the mesh (`snap = "faces"`). A
//!   sub-step crosses it when the parcel's cell changes side, which is
//!   exactly when the tracker took the parcel through one of the faces.
//! - [`Surface::Boundary`]: a boundary of the mesh; a parcel that reaches
//!   it is booked outwards, a parcel that bounces off it is booked back in.
//! - [`Surface::Triangles`]: a triangulated surface (STL), with the exact
//!   segment-triangle test of the surface hierarchy (D108 part 2).
//! - [`Volume`]: a box, a cylinder or a cell zone, for the parcels inside.
//!
//! Every crossing is a [`Crossing`] record. The records of a time step are
//! sorted by (datum, parcel id, sub-step, order in the sub-step) and summed
//! in that order (D249, D263): no atomics, the same bits for any number of
//! threads and any storage order of the parcels. The primary result is the
//! **cumulative** mass, droplet number and crossing count through a datum,
//! per direction; rates are differences over a report interval.
//!
//! The droplet statistics are sums of moments (D264): Σ N dᵏ for k = 0 to
//! 4, the mass, the velocity sums and a histogram of number and mass, so
//! that D10, D20, D30, D32, D43, the volume diameters Dv10, Dv50, Dv90 and
//! the span, the mean velocity and its standard deviation follow from
//! sums, and sums of several ranks merge (SP3).

use foehncraft_mesh::query::Bvh;
use serde::{Deserialize, Serialize};

/// The bound of a datum plane, in the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Bound {
    /// The whole plane.
    None,
    /// A disc about the plane's point.
    Disc {
        /// Radius (m).
        radius: f64,
    },
    /// A rectangle centred on the plane's point.
    Rectangle {
        /// Unit vector of the width, in the plane.
        e1: [f64; 3],
        /// Unit vector of the height, in the plane.
        e2: [f64; 3],
        /// Half width (m).
        half_width: f64,
        /// Half height (m).
        half_height: f64,
    },
}

impl Bound {
    /// Whether `x`, a point of the plane through `origin`, lies inside.
    pub fn contains(&self, origin: [f64; 3], x: [f64; 3]) -> bool {
        let r = sub(x, origin);
        match *self {
            Bound::None => true,
            Bound::Disc { radius } => dot(r, r) <= radius * radius,
            Bound::Rectangle {
                e1,
                e2,
                half_width,
                half_height,
            } => dot(r, e1).abs() <= half_width && dot(r, e2).abs() <= half_height,
        }
    }
}

/// A surface the parcels cross.
#[derive(Debug, Clone)]
pub enum Surface {
    /// A plane: point, unit normal (the positive direction), bound.
    Plane {
        /// A point of the plane.
        point: [f64; 3],
        /// Unit normal.
        normal: [f64; 3],
        /// The bound.
        bound: Bound,
    },
    /// Mesh faces between cells of different sides.
    Cells {
        /// Side of each local cell: −1, +1, or 0 for no side.
        side: Vec<i8>,
        /// The plane of a snapped plane: the crossing point is where the
        /// chord meets it (inside its bound); `None` for the faces between
        /// two zones (the crossing point is where the chord meets the
        /// faces, from `faces`).
        plane: Option<([f64; 3], [f64; 3], Bound)>,
        /// The faces, for the crossing point.
        faces: Option<Box<Bvh>>,
    },
    /// A boundary of the mesh (index into `mesh.boundaries`).
    Boundary {
        /// The boundary.
        boundary: usize,
    },
    /// A triangulated surface; the positive direction is that of the
    /// triangle normals.
    Triangles {
        /// The surface hierarchy.
        bvh: Box<Bvh>,
    },
}

/// A volume that holds parcels.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Volume {
    /// Axis-aligned box.
    Box {
        /// Lower corner.
        lower: [f64; 3],
        /// Upper corner.
        upper: [f64; 3],
    },
    /// Cylinder from `start` to `end` with `radius`.
    Cylinder {
        /// Centre of one end.
        start: [f64; 3],
        /// Centre of the other end.
        end: [f64; 3],
        /// Radius.
        radius: f64,
    },
    /// The cells of a cell zone (index into `mesh.zones`).
    Zone {
        /// The zone.
        zone: u32,
    },
}

impl Volume {
    /// Whether the point `x` of cell `cell` lies inside (`zone` is the
    /// cell zone of each cell).
    pub fn contains(&self, x: [f64; 3], cell: usize, zone: &[u32]) -> bool {
        match *self {
            Volume::Box { lower, upper } => (0..3).all(|k| x[k] >= lower[k] && x[k] <= upper[k]),
            Volume::Cylinder { start, end, radius } => {
                let a = sub(end, start);
                let l2 = dot(a, a);
                let r = sub(x, start);
                let s = dot(r, a);
                if s < 0.0 || s > l2 {
                    return false;
                }
                let along = s / l2;
                let q = [
                    r[0] - along * a[0],
                    r[1] - along * a[1],
                    r[2] - along * a[2],
                ];
                dot(q, q) <= radius * radius
            }
            Volume::Zone { zone: z } => zone.get(cell).is_some_and(|&c| c == z),
        }
    }
}

/// What a datum is for the parcels.
#[derive(Debug, Clone)]
pub enum Shape {
    /// A surface: crossings.
    Surface(Surface),
    /// A volume: the parcels inside.
    Volume(Volume),
}

/// The bins of a droplet size histogram.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HistSpec {
    /// Number of bins.
    pub bins: usize,
    /// Lower edge of the first bin (m).
    pub min: f64,
    /// Upper edge of the last bin (m).
    pub max: f64,
    /// Logarithmic bins (equal ratio) instead of linear ones.
    pub log: bool,
}

impl Default for HistSpec {
    fn default() -> Self {
        Self {
            bins: 40,
            min: 1e-6,
            max: 1e-3,
            log: true,
        }
    }
}

impl HistSpec {
    /// The bin edges (`bins + 1`).
    pub fn edges(&self) -> Vec<f64> {
        (0..=self.bins)
            .map(|i| {
                let f = i as f64 / self.bins as f64;
                if self.log {
                    self.min * (self.max / self.min).powf(f)
                } else {
                    self.min + f * (self.max - self.min)
                }
            })
            .collect()
    }

    /// The bin of diameter `d`: `Err(false)` below the first, `Err(true)`
    /// above the last.
    fn bin(&self, d: f64) -> Result<usize, bool> {
        if d < self.min {
            return Err(false);
        }
        if d >= self.max {
            return Err(true);
        }
        let f = if self.log {
            (d / self.min).ln() / (self.max / self.min).ln()
        } else {
            (d - self.min) / (self.max - self.min)
        };
        Ok(((f * self.bins as f64) as usize).min(self.bins - 1))
    }
}

/// A datum of the parcels.
#[derive(Debug, Clone)]
pub struct Datum {
    /// Name (column prefix of `datums.csv`).
    pub name: String,
    /// Surface or volume.
    pub shape: Shape,
    /// Histogram bins.
    pub hist: HistSpec,
}

/// One crossing of a datum surface by a parcel (D263).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Crossing {
    /// The datum (index).
    pub datum: u32,
    /// Parcel identity.
    pub id: u64,
    /// Sub-step of the parcel in this time step.
    pub sub: u64,
    /// Order within the sub-step.
    pub k: u32,
    /// +1 along the datum's direction, −1 against it.
    pub dir: i8,
    /// Time of the crossing (s).
    pub t: f64,
    /// Where (m).
    pub x: [f64; 3],
    /// Velocity (m/s).
    pub v: [f64; 3],
    /// Droplet diameter (m).
    pub d: f64,
    /// Droplets of the parcel.
    pub count: f64,
}

impl Crossing {
    /// The sort key of D263.
    pub fn key(&self) -> (u32, u64, u64, u32) {
        (self.datum, self.id, self.sub, self.k)
    }
}

/// Moment sums of droplets (D264): the basis of every statistic.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Moments {
    /// Parcels (or crossings) counted.
    pub parcels: u64,
    /// Σ N dᵏ for k = 0 to 4 (Σ N is the number of droplets).
    pub s: [f64; 5],
    /// Liquid mass (kg).
    pub mass: f64,
    /// Σ N v (droplet-number-weighted velocity sums).
    pub nv: [f64; 3],
    /// Σ N v² per component.
    pub nv2: [f64; 3],
    /// Σ m v (mass-weighted velocity sums).
    pub mv: [f64; 3],
    /// Droplets per bin.
    pub hist_n: Vec<f64>,
    /// Mass per bin (kg).
    pub hist_m: Vec<f64>,
    /// Droplets and mass below the first bin.
    pub under: [f64; 2],
    /// Droplets and mass above the last bin.
    pub over: [f64; 2],
}

impl Moments {
    /// Empty sums for `spec`.
    pub fn new(spec: &HistSpec) -> Self {
        Self {
            hist_n: vec![0.0; spec.bins],
            hist_m: vec![0.0; spec.bins],
            ..Default::default()
        }
    }

    /// Adds a parcel of `count` droplets of diameter `d`, mass `mass`,
    /// velocity `v`.
    pub fn add(&mut self, spec: &HistSpec, d: f64, count: f64, mass: f64, v: [f64; 3]) {
        self.parcels += 1;
        let mut dk = 1.0;
        for k in 0..5 {
            self.s[k] += count * dk;
            dk *= d;
        }
        self.mass += mass;
        for k in 0..3 {
            self.nv[k] += count * v[k];
            self.nv2[k] += count * v[k] * v[k];
            self.mv[k] += mass * v[k];
        }
        match spec.bin(d) {
            Ok(b) => {
                self.hist_n[b] += count;
                self.hist_m[b] += mass;
            }
            Err(false) => {
                self.under[0] += count;
                self.under[1] += mass;
            }
            Err(true) => {
                self.over[0] += count;
                self.over[1] += mass;
            }
        }
    }

    /// The statistics of the sums.
    pub fn stats(&self, spec: &HistSpec) -> Stats {
        let s = &self.s;
        let q = |a: f64, b: f64| if b > 0.0 { a / b } else { 0.0 };
        let n = s[0];
        let mut vm = [0.0; 3];
        let mut vs = [0.0; 3];
        let mut vmass = [0.0; 3];
        for k in 0..3 {
            vm[k] = q(self.nv[k], n);
            vs[k] = (q(self.nv2[k], n) - vm[k] * vm[k]).max(0.0).sqrt();
            vmass[k] = q(self.mv[k], self.mass);
        }
        let dv = |f: f64| self.volume_diameter(spec, f);
        let (dv10, dv50, dv90) = (dv(0.1), dv(0.5), dv(0.9));
        Stats {
            parcels: self.parcels,
            droplets: n,
            mass: self.mass,
            d10: q(s[1], s[0]),
            d20: q(s[2], s[0]).sqrt(),
            d30: q(s[3], s[0]).cbrt(),
            d32: q(s[3], s[2]),
            d43: q(s[4], s[3]),
            dv10,
            dv50,
            dv90,
            span: q(dv90 - dv10, dv50),
            velocity_mean: vm,
            velocity_std: vs,
            velocity_mass_mean: vmass,
            rosin_rammler: self.rosin_rammler(spec),
        }
    }

    /// The diameter below which the share `f` of the mass lies, from the
    /// mass histogram (interpolated in the bin, in ln d for log bins); the
    /// edge of the histogram when it lies in the under- or overflow.
    pub fn volume_diameter(&self, spec: &HistSpec, f: f64) -> f64 {
        let total: f64 = self.under[1] + self.hist_m.iter().sum::<f64>() + self.over[1];
        if total <= 0.0 {
            return 0.0;
        }
        let target = f * total;
        let edges = spec.edges();
        let mut cum = self.under[1];
        if cum >= target {
            return spec.min;
        }
        for (b, &m) in self.hist_m.iter().enumerate() {
            if m > 0.0 && cum + m >= target {
                let w = (target - cum) / m;
                let (a, c) = (edges[b], edges[b + 1]);
                return if spec.log {
                    a * (c / a).powf(w)
                } else {
                    a + w * (c - a)
                };
            }
            cum += m;
        }
        spec.max
    }

    /// The Rosin-Rammler fit of the mass histogram: least squares of
    /// ln(−ln Y) over ln d at the upper bin edges, Y the mass share above
    /// the edge, for 0.005 < Y < 0.995; `None` with fewer than 3 points.
    pub fn rosin_rammler(&self, spec: &HistSpec) -> Option<RosinRammlerFit> {
        let total: f64 = self.under[1] + self.hist_m.iter().sum::<f64>() + self.over[1];
        if total <= 0.0 {
            return None;
        }
        let edges = spec.edges();
        let mut cum = self.under[1];
        let mut pts = Vec::new();
        for (b, &m) in self.hist_m.iter().enumerate() {
            cum += m;
            let y = 1.0 - cum / total;
            if y > 0.005 && y < 0.995 {
                pts.push((edges[b + 1].ln(), (-y.ln()).ln()));
            }
        }
        if pts.len() < 3 {
            return None;
        }
        let k = pts.len() as f64;
        let mx = pts.iter().map(|p| p.0).sum::<f64>() / k;
        let my = pts.iter().map(|p| p.1).sum::<f64>() / k;
        let sxx: f64 = pts.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum();
        let sxy: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
        if sxx <= 0.0 {
            return None;
        }
        let n = sxy / sxx;
        let c = my - n * mx;
        let res = (pts
            .iter()
            .map(|p| (p.1 - (c + n * p.0)).powi(2))
            .sum::<f64>()
            / k)
            .sqrt();
        Some(RosinRammlerFit {
            mean: (-c / n).exp(),
            spread: n,
            residual: res,
            points: pts.len(),
        })
    }
}

/// The Rosin-Rammler fit of a histogram.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct RosinRammlerFit {
    /// Size parameter X (m).
    pub mean: f64,
    /// Spread n.
    pub spread: f64,
    /// Root mean square residual of ln(−ln Y).
    pub residual: f64,
    /// Points of the fit.
    pub points: usize,
}

/// The droplet statistics of moment sums (D264). Definitions: Dpq =
/// (Σ N dᵖ / Σ N d^q)^(1/(p−q)); Dv_f the diameter below which the share f
/// of the liquid mass lies (from the mass histogram); span = (Dv90 −
/// Dv10)/Dv50; velocities weighted by droplet number, and by mass.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Stats {
    /// Parcels (or crossings).
    pub parcels: u64,
    /// Droplets.
    pub droplets: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Arithmetic mean diameter (m).
    pub d10: f64,
    /// Surface mean diameter (m).
    pub d20: f64,
    /// Volume mean diameter (m).
    pub d30: f64,
    /// Sauter mean diameter (m).
    pub d32: f64,
    /// De Brouckere mean diameter (m).
    pub d43: f64,
    /// Dv10 (m).
    pub dv10: f64,
    /// Dv50, the mass median diameter (m).
    pub dv50: f64,
    /// Dv90 (m).
    pub dv90: f64,
    /// Relative span.
    pub span: f64,
    /// Mean droplet velocity (m/s).
    pub velocity_mean: [f64; 3],
    /// Standard deviation of the droplet velocity (m/s).
    pub velocity_std: [f64; 3],
    /// Mass-weighted mean velocity (m/s).
    pub velocity_mass_mean: [f64; 3],
    /// Rosin-Rammler fit of the mass histogram.
    pub rosin_rammler: Option<RosinRammlerFit>,
}

/// Liquid through a datum surface in one direction and both together.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Flux {
    /// Mass along the datum's direction (kg).
    pub mass_pos: f64,
    /// Mass against it (kg).
    pub mass_neg: f64,
    /// Droplets along the direction.
    pub droplets_pos: f64,
    /// Droplets against it.
    pub droplets_neg: f64,
    /// Crossings along the direction.
    pub crossings_pos: u64,
    /// Crossings against it.
    pub crossings_neg: u64,
    /// Moment sums of all crossings, both directions (flux-weighted
    /// statistics, D264).
    pub all: Moments,
}

impl Flux {
    fn new(spec: &HistSpec) -> Self {
        Self {
            all: Moments::new(spec),
            ..Default::default()
        }
    }

    /// Net mass along the direction (kg).
    pub fn net(&self) -> f64 {
        self.mass_pos - self.mass_neg
    }

    /// Crossings in both directions.
    pub fn crossings(&self) -> u64 {
        self.crossings_pos + self.crossings_neg
    }
}

/// The running sums of a datum surface: since the start and in the report
/// interval that is open (kept in the restart file).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DatumTally {
    /// Name of the datum.
    pub name: String,
    /// Since the start of the run (cumulative, the primary result).
    pub total: Flux,
    /// In the open report interval.
    pub interval: Flux,
    /// Start of the open interval (s).
    pub interval_start: f64,
    /// Start of the sums (s).
    pub start: f64,
    /// Net mass rate of every closed interval (kg/s), for the batch means.
    pub rates: Vec<f64>,
}

impl DatumTally {
    /// Empty sums from time `t0`.
    pub fn new(name: &str, spec: &HistSpec, t0: f64) -> Self {
        Self {
            name: name.into(),
            total: Flux::new(spec),
            interval: Flux::new(spec),
            interval_start: t0,
            start: t0,
            rates: Vec::new(),
        }
    }

    /// Books one crossing of mass `mass`.
    pub fn book(&mut self, spec: &HistSpec, c: &Crossing, mass: f64) {
        for f in [&mut self.total, &mut self.interval] {
            if c.dir > 0 {
                f.mass_pos += mass;
                f.droplets_pos += c.count;
                f.crossings_pos += 1;
            } else {
                f.mass_neg += mass;
                f.droplets_neg += c.count;
                f.crossings_neg += 1;
            }
            f.all.add(spec, c.d, c.count, mass, c.v);
        }
    }

    /// Closes the open interval at time `t`: its sums, its length; opens
    /// the next.
    pub fn close(&mut self, spec: &HistSpec, t: f64) -> (Flux, f64) {
        let dt = t - self.interval_start;
        let f = std::mem::replace(&mut self.interval, Flux::new(spec));
        if dt > 0.0 {
            self.rates.push(f.net() / dt);
        }
        self.interval_start = t;
        (f, dt)
    }

    /// The mean net mass rate since the start (kg/s) at time `t`.
    pub fn mean_rate(&self, t: f64) -> f64 {
        let dt = t - self.start;
        if dt > 0.0 { self.total.net() / dt } else { 0.0 }
    }

    /// Relative standard error of the mean rate from ten batch means of
    /// the interval rates (D71); `None` with fewer than ten intervals or a
    /// zero mean.
    pub fn rate_relative_error(&self) -> Option<f64> {
        batch_relative_error(&self.rates)
    }
}

/// Relative standard error of the mean of `x` from ten batch means.
pub fn batch_relative_error(x: &[f64]) -> Option<f64> {
    let b = 10usize;
    let size = x.len() / b;
    if size == 0 {
        return None;
    }
    let means: Vec<f64> = (0..b)
        .map(|i| x[i * size..(i + 1) * size].iter().sum::<f64>() / size as f64)
        .collect();
    let m = means.iter().sum::<f64>() / b as f64;
    if m == 0.0 {
        return None;
    }
    let var = means.iter().map(|y| (y - m) * (y - m)).sum::<f64>() / (b - 1) as f64;
    Some((var / b as f64).sqrt() / m.abs())
}

pub(crate) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Where the chord from `x0` to `x1` crosses the plane, as the chord
/// parameter and the direction, when the end points are on different sides
/// (s < 0 against s ≥ 0) and the crossing point lies in the bound.
pub fn plane_crossing(
    point: [f64; 3],
    normal: [f64; 3],
    bound: &Bound,
    x0: [f64; 3],
    x1: [f64; 3],
) -> Option<(f64, i8)> {
    let s0 = dot(sub(x0, point), normal);
    let s1 = dot(sub(x1, point), normal);
    if (s0 < 0.0) == (s1 < 0.0) {
        return None;
    }
    let f = if s0 == s1 {
        1.0
    } else {
        (s0 / (s0 - s1)).clamp(0.0, 1.0)
    };
    let x = [
        x0[0] + f * (x1[0] - x0[0]),
        x0[1] + f * (x1[1] - x0[1]),
        x0[2] + f * (x1[2] - x0[2]),
    ];
    if !bound.contains(point, x) {
        return None;
    }
    Some((f, if s1 >= 0.0 { 1 } else { -1 }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diameters_of_one_size_are_that_size() {
        let spec = HistSpec::default();
        let mut m = Moments::new(&spec);
        for _ in 0..10 {
            m.add(&spec, 50e-6, 3.0, 1e-9, [1.0, 0.0, 0.0]);
        }
        let s = m.stats(&spec);
        for d in [s.d10, s.d20, s.d30, s.d32, s.d43] {
            assert!((d / 50e-6 - 1.0).abs() < 1e-12, "{d}");
        }
        // Dv50 lies in the bin of 50 µm
        let e = spec.edges();
        let b = e.iter().position(|&x| x > 50e-6).unwrap();
        assert!(s.dv50 >= e[b - 1] && s.dv50 <= e[b]);
        assert_eq!(s.velocity_std, [0.0; 3]);
    }

    #[test]
    fn plane_crossings_are_watertight() {
        let p = [0.5, 0.0, 0.0];
        let n = [1.0, 0.0, 0.0];
        // a path in steps that ends exactly on the plane, then goes on
        let pts = [
            [0.1, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.9, 0.0, 0.0],
            [0.2, 0.0, 0.0],
        ];
        let mut net = 0;
        let mut count = 0;
        for w in pts.windows(2) {
            if let Some((_, d)) = plane_crossing(p, n, &Bound::None, w[0], w[1]) {
                net += d as i32;
                count += 1;
            }
        }
        assert_eq!((count, net), (2, 0));
        // outside a disc of radius 0.1
        let disc = Bound::Disc { radius: 0.1 };
        assert!(plane_crossing(p, n, &disc, [0.0, 0.2, 0.0], [1.0, 0.2, 0.0]).is_none());
        assert!(plane_crossing(p, n, &disc, [0.0, 0.05, 0.0], [1.0, 0.05, 0.0]).is_some());
    }

    #[test]
    fn rosin_rammler_fit_recovers_its_parameters() {
        let spec = HistSpec {
            bins: 200,
            min: 1e-6,
            max: 1e-3,
            log: true,
        };
        let mut m = Moments::new(&spec);
        // mass in each bin from the exact distribution X = 80 µm, n = 2.5
        let e = spec.edges();
        let y = |d: f64| (-(d / 80e-6f64).powf(2.5)).exp();
        for b in 0..spec.bins {
            m.hist_m[b] = y(e[b]) - y(e[b + 1]);
        }
        m.under[1] = 1.0 - y(spec.min);
        m.over[1] = y(spec.max);
        let f = m.rosin_rammler(&spec).unwrap();
        assert!((f.mean / 80e-6 - 1.0).abs() < 1e-6, "{f:?}");
        assert!((f.spread / 2.5 - 1.0).abs() < 1e-6, "{f:?}");
    }

    #[test]
    fn batch_error_needs_ten_intervals() {
        assert!(batch_relative_error(&[1.0; 9]).is_none());
        assert_eq!(batch_relative_error(&[2.0; 20]), Some(0.0));
    }
}
