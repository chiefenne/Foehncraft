//! The operator exactness check of a mesh (D235): the least-squares
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
}
