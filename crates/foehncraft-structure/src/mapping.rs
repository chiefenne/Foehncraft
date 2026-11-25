//! Temperatures of a solid region of a flow case on the structural nodes
//! (L53, stage MS6 of the moving solids plan).
//!
//! The flow solver keeps one temperature per cell, at the cell centre; the
//! structure needs it at its nodes (and, for beams and shells, its slope
//! across the section). At each point the [`Mapper`] fits a linear field
//! T(y) = a + g·(y − x) to the `k` nearest cell centres of the source by
//! weighted least squares and gives a (the temperature) and g (its
//! gradient):
//!
//! - exact for a linear temperature field on any mesh, second order for a
//!   smooth one; no search for the cell that contains the point, so a
//!   structural node a little outside the cells of the solid (a curved
//!   surface, a different mesh) takes the field extrapolated along its
//!   slope;
//! - directions the centres do not span (one layer of cells across a 2D
//!   case or a thin wall) take no slope (the minimum-norm solution of the
//!   fit): the field is constant across them;
//! - the weights 1/(d² + ε²) favour the nearest centres; they do not
//!   change the exactness for linear fields.
//!
//! A point farther from the nearest centre than the size of that cell
//! (the diagonal of its bounding box) or the given distance is refused: the
//! structural mesh does not lie in the solid region (units, position).

use std::path::Path;

/// The centres and values of the cells of a source field.
#[derive(Debug, Clone, Default)]
pub struct Source {
    /// Cell centres.
    pub centres: Vec<[f64; 3]>,
    /// Values at the centres.
    pub values: Vec<f64>,
    /// Size of each cell: the diagonal of its bounding box.
    pub sizes: Vec<f64>,
}

/// What a source read from a file holds, for the report.
#[derive(Debug, Clone, Default)]
pub struct SourceInfo {
    /// Cells of the whole file.
    pub cells: usize,
    /// Zones of the file.
    pub zones: Vec<String>,
}

impl Source {
    /// The cell field `field` of the cells in `zones` (every cell if
    /// `zones` is empty) of a VTKHDF result file of the flow solver.
    pub fn from_vtkhdf(
        path: &Path,
        field: &str,
        zones: &[String],
    ) -> Result<(Self, SourceInfo), String> {
        #[cfg(not(feature = "hdf5"))]
        {
            let _ = (path, field, zones);
            Err("this build has no HDF5: VTKHDF result files cannot be read".into())
        }
        #[cfg(feature = "hdf5")]
        Self::read_vtkhdf(path, field, zones)
    }

    #[cfg(feature = "hdf5")]
    fn read_vtkhdf(
        path: &Path,
        field: &str,
        zones: &[String],
    ) -> Result<(Self, SourceInfo), String> {
        use foehncraft_mesh::format::{ReadOptions, vtkhdf};
        use foehncraft_mesh::geometry::Geometry;
        let mf = vtkhdf::read_file(path, &ReadOptions::default())?;
        let mesh = mf.mesh;
        let values = vtkhdf::read_cell_fields(path, &[field])?
            .pop()
            .unwrap_or_default();
        if values.len() != mesh.n_cells {
            return Err(format!(
                "the field `{field}` has {} values for {} cells",
                values.len(),
                mesh.n_cells
            ));
        }
        let info = SourceInfo {
            cells: mesh.n_cells,
            zones: mesh.zones.clone(),
        };
        let mut want = vec![zones.is_empty(); mesh.zones.len()];
        for z in zones {
            match mesh.zones.iter().position(|m| m == z) {
                Some(i) => want[i] = true,
                None => {
                    return Err(format!(
                        "no zone `{z}` in the file; it has {}",
                        mesh.zones
                            .iter()
                            .map(|m| format!("`{m}`"))
                            .collect::<Vec<_>>()
                            .join(", ")
                    ));
                }
            }
        }
        let geo = Geometry::compute(&mesh);
        let mut src = Source::default();
        for c in 0..mesh.n_cells {
            if !want[mesh.cell_zone[c] as usize] {
                continue;
            }
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for n in mesh.cell_nodes(c) {
                let p = mesh.points[n as usize];
                for d in 0..3 {
                    lo[d] = lo[d].min(p[d]);
                    hi[d] = hi[d].max(p[d]);
                }
            }
            src.centres.push(geo.cell_centroid[c]);
            src.values.push(values[c]);
            src.sizes
                .push((0..3).map(|d| (hi[d] - lo[d]).powi(2)).sum::<f64>().sqrt());
        }
        if src.values.is_empty() {
            return Err("the zones have no cells".into());
        }
        if let Some(i) = src.values.iter().position(|v| !v.is_finite()) {
            return Err(format!(
                "the field `{field}` is not finite at the cell centred at {:?}",
                src.centres[i]
            ));
        }
        Ok((src, info))
    }
}

/// A linear fit at a point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fit {
    /// The value at the point.
    pub value: f64,
    /// The gradient.
    pub gradient: [f64; 3],
    /// Distance to the nearest centre.
    pub distance: f64,
    /// Size of the nearest cell.
    pub size: f64,
}

/// Nearest-centre search and the linear least-squares fit.
#[derive(Debug, Clone)]
pub struct Mapper {
    src: Source,
    k: usize,
    lo: [f64; 3],
    h: f64,
    dims: [usize; 3],
    /// Start of each bucket in `order` (buckets + 1).
    start: Vec<usize>,
    /// Centres sorted by bucket.
    order: Vec<usize>,
}

impl Mapper {
    /// A mapper over `src` that fits to the `k` nearest centres (at least
    /// 4, at most the number of centres).
    pub fn new(src: Source, k: usize) -> Self {
        let n = src.centres.len();
        let k = k.max(4).min(n.max(1));
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &src.centres {
            for d in 0..3 {
                lo[d] = lo[d].min(p[d]);
                hi[d] = hi[d].max(p[d]);
            }
        }
        if n == 0 {
            lo = [0.0; 3];
            hi = [0.0; 3];
        }
        // buckets of about two centres: the side from the extents that the
        // centres span (a flat set of centres spans two)
        let ext: Vec<f64> = (0..3).map(|d| hi[d] - lo[d]).collect();
        let big = ext.iter().cloned().fold(0.0, f64::max);
        let spanned: Vec<f64> = ext.iter().copied().filter(|&e| e > 1e-9 * big).collect();
        let h = if spanned.is_empty() {
            1.0
        } else {
            let vol: f64 = spanned.iter().product();
            (2.0 * vol / n as f64).powf(1.0 / spanned.len() as f64)
        }
        .max(big * 1e-6)
        .max(f64::MIN_POSITIVE);
        let dims = [0, 1, 2].map(|d| ((ext[d] / h).floor() as usize + 1).min(1 << 10));
        let nb = dims[0] * dims[1] * dims[2];
        let mut count = vec![0usize; nb + 1];
        let bucket: Vec<usize> = src
            .centres
            .iter()
            .map(|&p| {
                let b = Self::bucket_of(lo, h, dims, p);
                b[0] + dims[0] * (b[1] + dims[1] * b[2])
            })
            .collect();
        for &b in &bucket {
            count[b + 1] += 1;
        }
        for i in 0..nb {
            count[i + 1] += count[i];
        }
        let mut fill = count.clone();
        let mut order = vec![0; n];
        for (i, &b) in bucket.iter().enumerate() {
            order[fill[b]] = i;
            fill[b] += 1;
        }
        Self {
            src,
            k,
            lo,
            h,
            dims,
            start: count,
            order,
        }
    }

    fn bucket_of(lo: [f64; 3], h: f64, dims: [usize; 3], p: [f64; 3]) -> [usize; 3] {
        [0, 1, 2].map(|d| {
            let i = ((p[d] - lo[d]) / h).floor();
            if i.is_nan() || i < 0.0 {
                0
            } else {
                (i as usize).min(dims[d] - 1)
            }
        })
    }

    /// The source.
    pub fn source(&self) -> &Source {
        &self.src
    }

    /// The `k` nearest centres of `x`: (squared distance, centre), nearest
    /// first.
    pub fn nearest(&self, x: [f64; 3]) -> Vec<(f64, usize)> {
        let n = self.src.centres.len();
        let mut best: Vec<(f64, usize)> = Vec::with_capacity(self.k + 1);
        if n == 0 {
            return best;
        }
        let b = Self::bucket_of(self.lo, self.h, self.dims, x);
        // distance from x to the faces of its bucket box: what ring r
        // covers for certain
        let inside = (0..3)
            .map(|d| {
                let a = x[d] - (self.lo[d] + b[d] as f64 * self.h);
                let c = self.lo[d] + (b[d] + 1) as f64 * self.h - x[d];
                a.min(c).max(0.0)
            })
            .fold(f64::INFINITY, f64::min);
        let rmax = *self.dims.iter().max().unwrap_or(&1);
        for r in 0..=rmax {
            let lo = [0, 1, 2].map(|d| b[d].saturating_sub(r));
            let hi = [0, 1, 2].map(|d| (b[d] + r).min(self.dims[d] - 1));
            for k in lo[2]..=hi[2] {
                for j in lo[1]..=hi[1] {
                    for i in lo[0]..=hi[0] {
                        let ring = [i.abs_diff(b[0]), j.abs_diff(b[1]), k.abs_diff(b[2])]
                            .into_iter()
                            .max()
                            .unwrap_or(0);
                        if ring != r {
                            continue;
                        }
                        let id = i + self.dims[0] * (j + self.dims[1] * k);
                        for &c in &self.order[self.start[id]..self.start[id + 1]] {
                            let p = self.src.centres[c];
                            let d2: f64 = (0..3).map(|d| (p[d] - x[d]).powi(2)).sum();
                            if best.len() < self.k || d2 < best[best.len() - 1].0 {
                                let at = best.partition_point(|e| (e.0, e.1) < (d2, c));
                                best.insert(at, (d2, c));
                                best.truncate(self.k);
                            }
                        }
                    }
                }
            }
            // after ring r every centre within (r h + inside) has been seen
            let covered = r as f64 * self.h + inside;
            if best.len() == self.k.min(n) && best[best.len() - 1].0 <= covered * covered {
                break;
            }
        }
        best
    }

    /// The linear least-squares fit at `x`.
    pub fn fit(&self, x: [f64; 3]) -> Fit {
        let near = self.nearest(x);
        let Some(&(d0, c0)) = near.first() else {
            return Fit {
                value: f64::NAN,
                gradient: [0.0; 3],
                distance: f64::INFINITY,
                size: 0.0,
            };
        };
        let dk = near[near.len() - 1].0.sqrt();
        let eps2 = (1e-3 * dk).powi(2) + f64::MIN_POSITIVE;
        let mut wsum = 0.0;
        let mut xm = [0.0; 3];
        let mut tm = 0.0;
        for &(d2, c) in &near {
            let w = 1.0 / (d2 + eps2);
            let p = self.src.centres[c];
            wsum += w;
            for d in 0..3 {
                xm[d] += w * (p[d] - x[d]);
            }
            tm += w * self.src.values[c];
        }
        for v in &mut xm {
            *v /= wsum;
        }
        tm /= wsum;
        let mut cov = [[0.0; 3]; 3];
        let mut rhs = [0.0; 3];
        for &(d2, c) in &near {
            let w = 1.0 / (d2 + eps2);
            let p = self.src.centres[c];
            let y = [0, 1, 2].map(|d| p[d] - x[d] - xm[d]);
            let dt = self.src.values[c] - tm;
            for a in 0..3 {
                rhs[a] += w * y[a] * dt;
                for b in 0..3 {
                    cov[a][b] += w * y[a] * y[b];
                }
            }
        }
        let g = pseudo_solve(cov, rhs);
        let value = tm - (0..3).map(|d| g[d] * xm[d]).sum::<f64>();
        Fit {
            value,
            gradient: g,
            distance: d0.sqrt(),
            size: self.src.sizes.get(c0).copied().unwrap_or(0.0),
        }
    }
}

/// The minimum-norm solution of the symmetric positive semi-definite
/// system `a x = b`: eigen-directions with eigenvalues below 10⁻⁸ of the
/// largest are left out.
fn pseudo_solve(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let (lam, v) = jacobi3(a);
    let top = lam.iter().cloned().fold(0.0, f64::max);
    let mut x = [0.0; 3];
    if top <= 0.0 {
        return x;
    }
    for k in 0..3 {
        if lam[k] > 1e-8 * top {
            let col = [v[0][k], v[1][k], v[2][k]];
            let c = (0..3).map(|d| col[d] * b[d]).sum::<f64>() / lam[k];
            for d in 0..3 {
                x[d] += c * col[d];
            }
        }
    }
    x
}

/// Eigenvalues and eigenvectors (columns) of a symmetric 3 × 3 matrix by
/// cyclic Jacobi rotations.
fn jacobi3(mut a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let off = a[0][1].powi(2) + a[0][2].powi(2) + a[1][2].powi(2);
        let diag = a[0][0].powi(2) + a[1][1].powi(2) + a[2][2].powi(2);
        if off <= 1e-30 * diag || off == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q] == 0.0 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for k in 0..3 {
                let (akp, akq) = (a[k][p], a[k][q]);
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for k in 0..3 {
                let (apk, aqk) = (a[p][k], a[q][k]);
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for k in 0..3 {
                let (vkp, vkq) = (v[k][p], v[k][q]);
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jacobi_and_pseudo_inverse() {
        let a = [[4.0, 1.0, 0.5], [1.0, 3.0, 0.2], [0.5, 0.2, 2.0]];
        let x = [0.3, -1.2, 2.5];
        let b = [0, 1, 2].map(|i| (0..3).map(|j| a[i][j] * x[j]).sum::<f64>());
        let y = pseudo_solve(a, b);
        for d in 0..3 {
            assert!((x[d] - y[d]).abs() < 1e-12, "{y:?}");
        }
        // singular along z: the z component is left out
        let a = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]];
        let y = pseudo_solve(a, [2.0, 3.0, 0.0]);
        assert_eq!(y, [1.0, 3.0, 0.0]);
    }

    #[test]
    fn nearest_matches_brute_force() {
        let mut centres = Vec::new();
        let mut s = 12345u64;
        let mut rnd = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..500 {
            centres.push([rnd(), 0.3 * rnd(), 0.01]);
        }
        let n = centres.len();
        let src = Source {
            centres: centres.clone(),
            values: vec![0.0; n],
            sizes: vec![0.01; n],
        };
        let m = Mapper::new(src, 12);
        for _ in 0..50 {
            let x = [1.4 * rnd() - 0.2, 0.5 * rnd() - 0.1, 0.05 * rnd()];
            let got = m.nearest(x);
            let mut all: Vec<(f64, usize)> = centres
                .iter()
                .enumerate()
                .map(|(i, p)| ((0..3).map(|d| (p[d] - x[d]).powi(2)).sum(), i))
                .collect();
            all.sort_by(|a, b| a.partial_cmp(b).unwrap());
            assert_eq!(got, all[..12].to_vec());
        }
    }
}
