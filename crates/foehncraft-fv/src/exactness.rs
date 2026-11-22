//! The operator exactness check of a mesh (D236): the least-squares
//! gradient and the gradient of the gradient (the velocity Hessian of the
//! curvature correction of SST) applied to known linear and quadratic
//! fields on the real mesh before a solve, with the error per cell.
//!
//! A cell-centred operator is consistent when it reproduces the fields its
//! order promises: the least-squares gradient (D24) is exact for a linear
//! field on every mesh; it is second order where it is also exact for a
//! quadratic field, which needs a stencil symmetric about the cell
//! centroid (a regular hexahedral or polyhedral mesh). On irregular cells
//! (tetrahedra, perturbed or stretched cells, the cells next to a
//! boundary) the gradient of a quadratic has an error of order h, and an
//! operator that differentiates it once more (the stress divergence of the
//! QCR, the convective derivative of the strain of the curvature
//! correction) has an error of order one that does not fall with the mesh.
//!
//! Three dimensionless errors per cell, the largest over the test fields
//! (`dims` coordinates; a two-dimensional case leaves z out):
//!
//! - `linear`: |∇ₕφ − ∇φ| L for φ = (x_i − x̄_i)/L; round-off on a sound
//!   stencil, larger on a degenerate or badly conditioned one;
//! - `gradient`: |∇ₕq − ∇q| / (h_c |H|) for the quadratic fields
//!   q = (x_i − x̄_i)(x_j − x̄_j)/L², with H their Hessian (Frobenius norm)
//!   and h_c the longest stencil distance of the cell: zero where the
//!   gradient is second order, of order one where it is first order;
//! - `hessian`: |Hₕ − H| / |H| with Hₕ the least-squares gradient of the
//!   least-squares gradient, the second pass with a zero normal gradient
//!   on every boundary (as the curvature correction takes it): zero where
//!   that operator is consistent.
//!
//! The boundary stencil points take the exact values of the test field
//! (the face centroid of a fixed value, the projected point of the other
//! kinds), so the check measures the geometry of the stencils, not the
//! boundary data of a case. x̄ is the centre of the bounding box and L its
//! largest side.

use crate::disc::lsq_weights;
use crate::problem::{BoundaryKind, Problem};
use anyhow::{Result, bail};
use foehncraft_kernels::unpack_face;
use foehncraft_mesh::geometry::{add, dot, norm, scale, sub};
use foehncraft_mesh::{Geometry, Mesh};
use rayon::prelude::*;

/// The errors of the operators per cell (see the module documentation).
#[derive(Debug, Clone, Default)]
pub struct OperatorErrors {
    /// The gradient of linear fields.
    pub linear: Vec<f64>,
    /// The gradient of quadratic fields, relative to h_c |H|.
    pub gradient: Vec<f64>,
    /// The gradient of the gradient of quadratic fields, relative to |H|.
    pub hessian: Vec<f64>,
}

/// The kind of every boundary face from the conditions of `problem`
/// (the stencil points of [`crate::Discretisation::least_squares`]).
pub fn boundary_kinds(mesh: &Mesh, problem: &Problem) -> Result<Vec<BoundaryKind>> {
    let ni = mesh.n_interior_faces();
    let mut kind = vec![BoundaryKind::ZeroGradient; mesh.n_faces() - ni];
    for p in &mesh.boundaries {
        let Some(bc) = problem.boundary.get(&p.name) else {
            bail!("no boundary condition for boundary `{}`", p.name);
        };
        for f in p.range() {
            kind[f - ni] = bc.kind;
        }
    }
    Ok(kind)
}

/// The stencil of a cell: per entry of `cell_faces` the point (the other
/// centroid, shifted across a periodic face, or the boundary point) and
/// whether it is a boundary point.
fn stencil(
    mesh: &Mesh,
    geo: &Geometry,
    kind: &[BoundaryKind],
    c: usize,
) -> Vec<(usize, [f64; 3], Option<usize>)> {
    let ni = mesh.n_interior_faces();
    let periodic = mesh.periodic_start();
    let pshift: &[[f64; 3]] = mesh.periodic.as_ref().map_or(&[], |p| &p.shift);
    let xc = &geo.cell_centroid;
    let lo = mesh.cell_ptr[c] as usize;
    let hi = mesh.cell_ptr[c + 1] as usize;
    mesh.cell_faces[lo..hi]
        .iter()
        .map(|&e| {
            let (f, is_nb) = unpack_face(e);
            if f < ni {
                if is_nb {
                    let o = mesh.owner[f] as usize;
                    // the owner seen from the neighbour of a periodic face
                    let x = if f >= periodic {
                        sub(xc[o], pshift[f - periodic])
                    } else {
                        xc[o]
                    };
                    (f, x, None)
                } else {
                    let n = mesh.neighbour[f] as usize;
                    let x = if f >= periodic {
                        add(xc[n], pshift[f - periodic])
                    } else {
                        xc[n]
                    };
                    (f, x, None)
                }
            } else {
                let d = sub(geo.face_centroid[f], xc[c]);
                let d = match kind[f - ni] {
                    BoundaryKind::FixedValue => d,
                    _ => {
                        let s = geo.face_area[f];
                        let n = scale(s, 1.0 / norm(s));
                        scale(n, dot(d, n))
                    }
                };
                (f, add(xc[c], d), Some(f - ni))
            }
        })
        .collect()
}

/// The errors of the least-squares gradient with the boundary kinds `kind`
/// (the velocity of a flow) and of its gradient (the second pass with a zero
/// normal gradient on every boundary), on `dims` coordinates (2 or 3).
pub fn operator_errors(
    mesh: &Mesh,
    geo: &Geometry,
    kind: &[BoundaryKind],
    dims: usize,
) -> Result<OperatorErrors> {
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let dims = dims.clamp(1, 3);
    let w1 = lsq_weights(mesh, geo, kind)?;
    let zero = vec![BoundaryKind::ZeroGradient; kind.len()];
    let w2 = lsq_weights(mesh, geo, &zero)?;
    let (lo, hi) = mesh.bounding_box();
    let centre: [f64; 3] = std::array::from_fn(|i| 0.5 * (lo[i] + hi[i]));
    let len = (0..3)
        .map(|i| hi[i] - lo[i])
        .fold(0.0, f64::max)
        .max(1e-300);
    let xc = &geo.cell_centroid;
    let rel = |x: [f64; 3], i: usize| (x[i] - centre[i]) / len;
    // the quadratic fields (i, j), i ≤ j
    let pairs: Vec<(usize, usize)> = (0..dims)
        .flat_map(|i| (i..dims).map(move |j| (i, j)))
        .collect();
    let q = |x: [f64; 3], (i, j): (usize, usize)| rel(x, i) * rel(x, j);
    let dq = |x: [f64; 3], (i, j): (usize, usize)| {
        let mut g = [0.0; 3];
        g[i] += rel(x, j) / len;
        g[j] += rel(x, i) / len;
        g
    };
    let hess = |(i, j): (usize, usize)| {
        let mut h = [[0.0; 3]; 3];
        h[i][j] += 1.0 / (len * len);
        h[j][i] += 1.0 / (len * len);
        h
    };
    let frob = |h: &[[f64; 3]; 3]| h.iter().flatten().map(|v| v * v).sum::<f64>().sqrt();
    let lsq = |w: &crate::disc::Soa3, c: usize, vals: &[f64], vc: f64| {
        let lo = mesh.cell_ptr[c] as usize;
        let mut g = [0.0; 3];
        for (k, v) in vals.iter().enumerate() {
            let wk = [w.x[lo + k], w.y[lo + k], w.z[lo + k]];
            for r in 0..3 {
                g[r] += wk[r] * (v - vc);
            }
        }
        g
    };
    let stencils = |c: usize| stencil(mesh, geo, kind, c);
    let h_c: Vec<f64> = (0..nc)
        .into_par_iter()
        .map(|c| {
            stencils(c)
                .iter()
                .map(|(_, x, _)| norm(sub(*x, xc[c])))
                .fold(0.0, f64::max)
                .max(1e-300)
        })
        .collect();
    // the linear fields
    let linear: Vec<f64> = (0..nc)
        .into_par_iter()
        .map(|c| {
            let st = stencils(c);
            let mut worst = 0.0f64;
            for i in 0..dims {
                let vals: Vec<f64> = st.iter().map(|(_, x, _)| rel(*x, i)).collect();
                let mut e = lsq(&w1, c, &vals, rel(xc[c], i));
                e[i] -= 1.0 / len;
                worst = worst.max(norm(e) * len);
            }
            worst
        })
        .collect();
    let mut gradient = vec![0.0f64; nc];
    let mut hessian = vec![0.0f64; nc];
    let periodic = mesh.periodic_start();
    let pshift: &[[f64; 3]] = mesh.periodic.as_ref().map_or(&[], |p| &p.shift);
    for &p in &pairs {
        let h = hess(p);
        let hn = frob(&h);
        // pass 1: the gradient of the quadratic field
        let gq: Vec<[f64; 3]> = (0..nc)
            .into_par_iter()
            .map(|c| {
                let vals: Vec<f64> = stencils(c).iter().map(|(_, x, _)| q(*x, p)).collect();
                lsq(&w1, c, &vals, q(xc[c], p))
            })
            .collect();
        // pass 2: the gradient of each of its components, zero normal
        // gradient on the boundary; across a periodic face the neighbour's
        // gradient continued with the exact Hessian (the test fields are
        // not periodic, the fields of a run are)
        let hc: Vec<f64> = (0..nc)
            .into_par_iter()
            .map(|c| {
                let lo = mesh.cell_ptr[c] as usize;
                let hi = mesh.cell_ptr[c + 1] as usize;
                let gc = gq[c];
                let mut hh = [[0.0; 3]; 3];
                for b in 0..3 {
                    let vals: Vec<f64> = mesh.cell_faces[lo..hi]
                        .iter()
                        .map(|&e| {
                            let (f, is_nb) = unpack_face(e);
                            if f >= ni {
                                return gc[b];
                            }
                            let (other, sign) = if is_nb {
                                (mesh.owner[f] as usize, -1.0)
                            } else {
                                (mesh.neighbour[f] as usize, 1.0)
                            };
                            let mut v = gq[other][b];
                            if f >= periodic {
                                let s = pshift[f - periodic];
                                v += sign * (h[b][0] * s[0] + h[b][1] * s[1] + h[b][2] * s[2]);
                            }
                            v
                        })
                        .collect();
                    // H_ab = ∂_a (∂_b q)
                    let g = lsq(&w2, c, &vals, gc[b]);
                    for a in 0..3 {
                        hh[a][b] = g[a] - h[a][b];
                    }
                }
                frob(&hh) / hn
            })
            .collect();
        for c in 0..nc {
            let e = sub(gq[c], dq(xc[c], p));
            gradient[c] = gradient[c].max(norm(e) / (h_c[c] * hn));
            hessian[c] = hessian[c].max(hc[c]);
        }
    }
    Ok(OperatorErrors {
        linear,
        gradient,
        hessian,
    })
}

/// The quadratic least-squares gradient of the cells `cells` (D236, step
/// 2): a weighted fit of a linear and a quadratic part,
/// φ_k − φ_c = ∇φ·d_k + ½ d_kᵀ H d_k, over the face neighbours, their face
/// neighbours and the fixed-value boundary faces of both (the boundary
/// value at the face centroid), with the weights 1/|d_k|² of the
/// face-neighbour gradient (D24). The gradient is exact for a quadratic
/// field wherever the fit has full rank. Per treated cell the rows of the
/// gradient: g = Σ_k w_k (φ_k − φ_c), with `idx` a cell (< n_cells) or a
/// boundary face (n_cells + boundary face index).
#[derive(Debug, Clone, Default)]
pub struct QuadraticGradient {
    /// The treated cells.
    pub cells: Vec<u32>,
    /// Start of the stencil of each treated cell in `idx` and `w`
    /// (`cells.len() + 1` entries).
    pub ptr: Vec<u32>,
    /// Stencil points: a cell, or `n_cells` + a boundary face index.
    pub idx: Vec<u32>,
    /// Gradient weights per stencil point.
    pub w: crate::disc::Soa3,
    /// Hessian weights per stencil point, H_mn = Σ_k h[s][k] (φ_k − φ_c)
    /// with s = 0 … 5 for xx, yy, zz, xy, yz, zx (the order of the strain
    /// components of the curvature correction); zero beyond `dims`.
    pub h: [Vec<f64>; 6],
    /// Cells asked for whose fit is rank deficient or badly conditioned
    /// (left to the face-neighbour gradient).
    pub failed: Vec<usize>,
    /// Coordinates of the fit (2: x and y).
    pub dims: usize,
}

/// Smallest pivot of the scaled normal equations, relative to the
/// largest, for a fit to count as full rank.
pub const QUADRATIC_PIVOT: f64 = 1e-9;

/// The points of a fitted stencil with their gradient and Hessian weights.
type Stencil = Vec<(u32, [f64; 3], [f64; 6])>;

/// The index of the Hessian component (i, j), i ≤ j, in the order xx, yy,
/// zz, xy, yz, zx.
const fn hessian_slot(i: usize, j: usize) -> usize {
    match (i, j) {
        (0, 0) => 0,
        (1, 1) => 1,
        (2, 2) => 2,
        (0, 1) => 3,
        (1, 2) => 4,
        _ => 5,
    }
}

/// Cholesky solve of the symmetric positive definite `n`×`n` matrix `a`
/// (row major) for the columns of `b` (`n` × `m`), in place; `None` when a
/// pivot falls below `QUADRATIC_PIVOT` times the largest diagonal entry.
fn cholesky_solve(a: &mut [f64], n: usize, b: &mut [f64], m: usize) -> Option<()> {
    let big = (0..n).map(|i| a[i * n + i]).fold(0.0, f64::max);
    if big.is_nan() || big <= 0.0 {
        return None;
    }
    for j in 0..n {
        let mut d = a[j * n + j];
        for k in 0..j {
            d -= a[j * n + k] * a[j * n + k];
        }
        if d.is_nan() || d <= QUADRATIC_PIVOT * big {
            return None;
        }
        let d = d.sqrt();
        a[j * n + j] = d;
        for i in j + 1..n {
            let mut v = a[i * n + j];
            for k in 0..j {
                v -= a[i * n + k] * a[j * n + k];
            }
            a[i * n + j] = v / d;
        }
    }
    for c in 0..m {
        // L y = b, Lᵀ x = y
        for i in 0..n {
            let mut v = b[i * m + c];
            for k in 0..i {
                v -= a[i * n + k] * b[k * m + c];
            }
            b[i * m + c] = v / a[i * n + i];
        }
        for i in (0..n).rev() {
            let mut v = b[i * m + c];
            for k in i + 1..n {
                v -= a[k * n + i] * b[k * m + c];
            }
            b[i * m + c] = v / a[i * n + i];
        }
    }
    Some(())
}

/// The stencil points of the quadratic fit of cell `c`: (index, position).
fn quadratic_points(
    mesh: &Mesh,
    geo: &Geometry,
    kind: &[BoundaryKind],
    c: usize,
    min: usize,
) -> Vec<(u32, [f64; 3])> {
    let ni = mesh.n_interior_faces();
    let nc = mesh.n_cells;
    let periodic = mesh.periodic_start();
    let pshift: &[[f64; 3]] = mesh.periodic.as_ref().map_or(&[], |p| &p.shift);
    let xc = &geo.cell_centroid;
    // the face neighbours of cell `a` seen with the shift `s` of `a`
    let ring = |a: usize, s: [f64; 3]| -> Vec<(u32, [f64; 3], [f64; 3])> {
        let mut out = Vec::new();
        for (f, is_nb) in mesh.cell_face_iter(a) {
            if f < ni {
                let (o, sh) = if is_nb {
                    let sh = if f >= periodic {
                        scale(pshift[f - periodic], -1.0)
                    } else {
                        [0.0; 3]
                    };
                    (mesh.owner[f] as usize, sh)
                } else {
                    let sh = if f >= periodic {
                        pshift[f - periodic]
                    } else {
                        [0.0; 3]
                    };
                    (mesh.neighbour[f] as usize, sh)
                };
                let s2 = add(s, sh);
                out.push((o as u32, add(xc[o], s2), s2));
            } else if kind[f - ni] == BoundaryKind::FixedValue {
                out.push(((nc + f - ni) as u32, add(geo.face_centroid[f], s), s));
            }
        }
        out
    };
    let same = |a: &[f64; 3], b: &[f64; 3]| {
        (0..3).all(|i| (a[i] - b[i]).abs() < 1e-12 * (1.0 + a[i].abs()))
    };
    let mut pts: Vec<(u32, [f64; 3])> = Vec::new();
    let seen = |pts: &mut Vec<(u32, [f64; 3])>, i: u32, x: [f64; 3]| {
        if (i as usize == c && same(&x, &xc[c])) || pts.iter().any(|(j, y)| *j == i && same(y, &x))
        {
            return;
        }
        pts.push((i, x));
    };
    // the neighbours of the neighbours; a third ring where that leaves
    // fewer than `min` points (corner cells)
    let mut frontier: Vec<(u32, [f64; 3])> = vec![(c as u32, [0.0; 3])];
    for layer in 0..3 {
        if layer == 2 && pts.len() >= min {
            break;
        }
        let mut next = Vec::new();
        for &(a, s) in &frontier {
            if (a as usize) >= nc {
                continue;
            }
            for (j, y, s2) in ring(a as usize, s) {
                let before = pts.len();
                seen(&mut pts, j, y);
                if pts.len() > before {
                    next.push((j, s2));
                }
            }
        }
        frontier = next;
    }
    pts
}

/// The quadratic least-squares gradient of `cells` (see
/// [`QuadraticGradient`]) on `dims` coordinates, with the boundary kinds
/// `kind` of the field (fixed-value faces enter the fit).
pub fn quadratic_gradient(
    mesh: &Mesh,
    geo: &Geometry,
    kind: &[BoundaryKind],
    dims: usize,
    cells: &[usize],
) -> QuadraticGradient {
    let dims = dims.clamp(1, 3);
    let nq = dims * (dims + 1) / 2;
    let m = dims + nq;
    let xc = &geo.cell_centroid;
    let fits: Vec<Option<Stencil>> = cells
        .par_iter()
        .map(|&c| {
            let pts = quadratic_points(mesh, geo, kind, c, 2 * m);
            if pts.len() < m + 1 {
                return None;
            }
            let ds: Vec<[f64; 3]> = pts.iter().map(|(_, x)| sub(*x, xc[c])).collect();
            let h = ds.iter().map(|d| norm(*d)).fold(0.0, f64::max);
            if h.is_nan() || h <= 0.0 {
                return None;
            }
            let n = pts.len();
            // rows a_k of the scaled design matrix and the weights
            let mut rows = vec![0.0; n * m];
            let mut wt = vec![0.0; n];
            for (k, d) in ds.iter().enumerate() {
                let e = scale(*d, 1.0 / h);
                let r = &mut rows[k * m..(k + 1) * m];
                r[..dims].copy_from_slice(&e[..dims]);
                let mut q = dims;
                for i in 0..dims {
                    for j in i..dims {
                        r[q] = if i == j {
                            0.5 * e[i] * e[i]
                        } else {
                            e[i] * e[j]
                        };
                        q += 1;
                    }
                }
                wt[k] = 1.0 / dot(e, e).max(1e-300);
            }
            let mut a = vec![0.0; m * m];
            let mut b = vec![0.0; m * n];
            for k in 0..n {
                let r = &rows[k * m..(k + 1) * m];
                for i in 0..m {
                    for j in 0..m {
                        a[i * m + j] += wt[k] * r[i] * r[j];
                    }
                    b[i * n + k] = wt[k] * r[i];
                }
            }
            cholesky_solve(&mut a, m, &mut b, n)?;
            Some(
                pts.iter()
                    .enumerate()
                    .map(|(k, (i, _))| {
                        let mut g = [0.0; 3];
                        for (r, gr) in g.iter_mut().enumerate().take(dims) {
                            *gr = b[r * n + k] / h;
                        }
                        // the quadratic rows hold ½H_ii and H_ij (i < j)
                        let mut hw = [0.0; 6];
                        let mut q = dims;
                        for a in 0..dims {
                            for c in a..dims {
                                hw[hessian_slot(a, c)] = b[q * n + k] / (h * h);
                                q += 1;
                            }
                        }
                        (*i, g, hw)
                    })
                    .collect(),
            )
        })
        .collect();
    let mut out = QuadraticGradient {
        ptr: vec![0],
        dims,
        ..Default::default()
    };
    for (&c, fit) in cells.iter().zip(fits) {
        match fit {
            Some(v) => {
                out.cells.push(c as u32);
                for (i, g, hw) in v {
                    out.idx.push(i);
                    out.w.x.push(g[0]);
                    out.w.y.push(g[1]);
                    out.w.z.push(g[2]);
                    for (o, x) in out.h.iter_mut().zip(hw) {
                        o.push(x);
                    }
                }
                out.ptr.push(out.idx.len() as u32);
            }
            None => out.failed.push(c),
        }
    }
    out
}

/// The gradient error of the treated cells of `q` (built with the
/// boundary kinds `kind`) for the quadratic test fields, as `gradient` of
/// [`OperatorErrors`] (relative to h_c |H| with h_c the longest distance of
/// the cell's face-neighbour stencil): zero up to round-off where the fit
/// has full rank. Per treated cell, in the order of `q.cells`.
pub fn quadratic_gradient_errors(
    mesh: &Mesh,
    geo: &Geometry,
    kind: &[BoundaryKind],
    q: &QuadraticGradient,
) -> Vec<f64> {
    let ni = mesh.n_interior_faces();
    let dims = q.dims;
    let (lo, hi) = mesh.bounding_box();
    let centre: [f64; 3] = std::array::from_fn(|i| 0.5 * (lo[i] + hi[i]));
    let len = (0..3)
        .map(|i| hi[i] - lo[i])
        .fold(0.0, f64::max)
        .max(1e-300);
    let xc = &geo.cell_centroid;
    (0..q.cells.len())
        .into_par_iter()
        .map(|t| {
            let c = q.cells[t] as usize;
            let pts = quadratic_points(mesh, geo, kind, c, 2 * (dims + dims * (dims + 1) / 2));
            let hc = mesh
                .cell_face_iter(c)
                .map(|(f, nb)| {
                    let x = if f < ni {
                        let o = if nb { mesh.owner[f] } else { mesh.neighbour[f] } as usize;
                        xc[o]
                    } else {
                        geo.face_centroid[f]
                    };
                    norm(sub(x, xc[c]))
                })
                .fold(0.0, f64::max)
                .max(1e-300);
            let a = q.ptr[t] as usize;
            let mut worst = 0.0f64;
            for i in 0..dims {
                for j in i..dims {
                    let f = |x: [f64; 3]| (x[i] - centre[i]) * (x[j] - centre[j]) / (len * len);
                    let mut g = [0.0; 3];
                    for (k, (_, x)) in pts.iter().enumerate() {
                        let v = f(*x) - f(xc[c]);
                        g[0] += q.w.x[a + k] * v;
                        g[1] += q.w.y[a + k] * v;
                        g[2] += q.w.z[a + k] * v;
                    }
                    let mut e = [0.0; 3];
                    e[i] += (xc[c][j] - centre[j]) / (len * len);
                    e[j] += (xc[c][i] - centre[i]) / (len * len);
                    let hn = if i == j { 2.0 } else { 2f64.sqrt() } / (len * len);
                    worst = worst.max(norm(sub(g, e)) / (hc * hn));
                }
            }
            worst
        })
        .collect()
}

/// The Hessian error of the treated cells of `q` for the quadratic test
/// fields, as `hessian` of [`OperatorErrors`] (|Hₕ − H|/|H|): zero up to
/// round-off where the fit has full rank. Per treated cell, in the order
/// of `q.cells`.
pub fn quadratic_hessian_errors(
    mesh: &Mesh,
    geo: &Geometry,
    kind: &[BoundaryKind],
    q: &QuadraticGradient,
) -> Vec<f64> {
    let dims = q.dims;
    let (lo, hi) = mesh.bounding_box();
    let centre: [f64; 3] = std::array::from_fn(|i| 0.5 * (lo[i] + hi[i]));
    let len = (0..3)
        .map(|i| hi[i] - lo[i])
        .fold(0.0, f64::max)
        .max(1e-300);
    let xc = &geo.cell_centroid;
    (0..q.cells.len())
        .into_par_iter()
        .map(|t| {
            let c = q.cells[t] as usize;
            let pts = quadratic_points(mesh, geo, kind, c, 2 * (dims + dims * (dims + 1) / 2));
            let a = q.ptr[t] as usize;
            let mut worst = 0.0f64;
            for i in 0..dims {
                for j in i..dims {
                    let f = |x: [f64; 3]| (x[i] - centre[i]) * (x[j] - centre[j]) / (len * len);
                    let mut diff = 0.0;
                    for m in 0..dims {
                        for n in m..dims {
                            let s = hessian_slot(m, n);
                            let hm: f64 = pts
                                .iter()
                                .enumerate()
                                .map(|(k, (_, x))| q.h[s][a + k] * (f(*x) - f(xc[c])))
                                .sum();
                            let exact = if (m, n) == (i, j) {
                                if i == j { 2.0 } else { 1.0 }
                            } else {
                                0.0
                            } / (len * len);
                            let mult = if m == n { 1.0 } else { 2.0 };
                            diff += mult * (hm - exact).powi(2);
                        }
                    }
                    let hn = if i == j { 2.0 } else { 2f64.sqrt() } / (len * len);
                    worst = worst.max(diff.sqrt() / hn);
                }
            }
            worst
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verify::test_mesh;
    use foehncraft_mesh::generate::Kind;

    fn errors(kind: Kind, perturb: f64) -> (Mesh, OperatorErrors) {
        let (mesh, geo) = test_mesh(kind, 3, 6, 0, perturb).unwrap();
        let kinds = vec![BoundaryKind::FixedValue; mesh.n_faces() - mesh.n_interior_faces()];
        let e = operator_errors(&mesh, &geo, &kinds, 3).unwrap();
        (mesh, e)
    }

    /// The largest value over the cells two layers away from the boundary.
    fn inner_max(mesh: &Mesh, vals: &[f64]) -> f64 {
        let ni = mesh.n_interior_faces();
        let at_boundary: Vec<bool> = (0..mesh.n_cells)
            .map(|c| mesh.cell_face_iter(c).any(|(f, _)| f >= ni))
            .collect();
        (0..mesh.n_cells)
            .filter(|&c| {
                !at_boundary[c]
                    && mesh.cell_face_iter(c).all(|(f, nb)| {
                        let o = if nb { mesh.owner[f] } else { mesh.neighbour[f] };
                        !at_boundary[o as usize]
                    })
            })
            .map(|c| vals[c])
            .fold(0.0, f64::max)
    }

    #[test]
    fn linear_fields_are_exact_on_every_mesh() {
        for (m, p) in [
            (Kind::Hex, 0.0),
            (Kind::Hex, 0.3),
            (Kind::Tet, 0.0),
            (Kind::Tet, 0.2),
        ] {
            let (_, e) = errors(m, p);
            let worst = e.linear.iter().cloned().fold(0.0, f64::max);
            assert!(worst < 1e-10, "{m:?} {p}: {worst:e}");
        }
    }

    #[test]
    fn regular_hexahedra_are_second_order_inside() {
        let (mesh, e) = errors(Kind::Hex, 0.0);
        let g = inner_max(&mesh, &e.gradient);
        let h = inner_max(&mesh, &e.hessian);
        assert!(g < 1e-10 && h < 1e-8, "gradient {g:e}, hessian {h:e}");
        // the cells at the boundary are not (one-sided stencils)
        assert!(e.gradient.iter().cloned().fold(0.0, f64::max) > 0.01);
    }

    #[test]
    fn irregular_cells_are_flagged() {
        for (m, p) in [(Kind::Hex, 0.3), (Kind::Tet, 0.0)] {
            let (mesh, e) = errors(m, p);
            let g = inner_max(&mesh, &e.gradient);
            let h = inner_max(&mesh, &e.hessian);
            assert!(
                g > 0.01 && h > 0.05,
                "{m:?} {p}: gradient {g:e}, hessian {h:e}"
            );
        }
    }

    #[test]
    fn the_quadratic_gradient_is_exact_for_quadratic_fields() {
        for (m, p) in [
            (Kind::Hex, 0.3),
            (Kind::Tet, 0.0),
            (Kind::Tet, 0.2),
            (Kind::Poly, 0.2),
        ] {
            let (mesh, geo) = test_mesh(m, 3, 4, 0, p).unwrap();
            let kinds = vec![BoundaryKind::FixedValue; mesh.n_faces() - mesh.n_interior_faces()];
            let cells: Vec<usize> = (0..mesh.n_cells).collect();
            let q = quadratic_gradient(&mesh, &geo, &kinds, 3, &cells);
            assert!(q.failed.is_empty(), "{m:?} {p}: {} failed", q.failed.len());
            let e = quadratic_gradient_errors(&mesh, &geo, &kinds, &q);
            let worst = e.iter().cloned().fold(0.0, f64::max);
            assert!(worst < 1e-8, "{m:?} {p}: {worst:e}");
            let e = quadratic_hessian_errors(&mesh, &geo, &kinds, &q);
            let worst = e.iter().cloned().fold(0.0, f64::max);
            assert!(worst < 1e-7, "{m:?} {p}: Hessian {worst:e}");
        }
    }

    #[test]
    fn the_quadratic_gradient_in_two_dimensions() {
        let (mesh, geo) = test_mesh(Kind::Tet, 2, 6, 0, 0.2).unwrap();
        let ni = mesh.n_interior_faces();
        // walls fixed, the z sides symmetry (not in the fit)
        let mut kinds = vec![BoundaryKind::FixedValue; mesh.n_faces() - ni];
        for p in &mesh.boundaries {
            if p.name.starts_with('z') {
                for f in p.range() {
                    kinds[f - ni] = BoundaryKind::ZeroGradient;
                }
            }
        }
        let cells: Vec<usize> = (0..mesh.n_cells).collect();
        let q = quadratic_gradient(&mesh, &geo, &kinds, 2, &cells);
        let names: Vec<_> = mesh.boundaries.iter().map(|p| p.name.clone()).collect();
        assert!(
            q.failed.is_empty(),
            "{} of {} failed; boundaries {names:?}; points of the first {}",
            q.failed.len(),
            mesh.n_cells,
            quadratic_points(&mesh, &geo, &kinds, q.failed[0], 10).len()
        );
        let e = quadratic_gradient_errors(&mesh, &geo, &kinds, &q);
        assert!(e.iter().cloned().fold(0.0, f64::max) < 1e-8);
        let e = quadratic_hessian_errors(&mesh, &geo, &kinds, &q);
        assert!(e.iter().cloned().fold(0.0, f64::max) < 1e-7);
    }
}
