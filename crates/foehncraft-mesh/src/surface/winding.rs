//! Generalised winding numbers (Jacobson, Kavan & Sorkine-Hornung 2013)
//! with the fast hierarchical evaluation of Barill, Dickson, Schmidt, Levin
//! & Jacobson 2018 ("Fast winding numbers for soups and clouds", ACM TOG
//! 37(4)), on the surface hierarchy of D108. Written from the papers, no
//! code read (D201).
//!
//! The winding number of a point `q` is the solid angle the surface
//! subtends at `q`, divided by 4π: 1 inside a closed outward surface, 0
//! outside, and a smooth "confidence" in between for a surface with holes
//! or overlaps, which is why it is the inside test for dirty STL.
//!
//! **Exact part.** The solid angle of one triangle is the formula of
//! Van Oosterom & Strackee 1983 (IEEE Trans. Biomed. Eng. 30(2)):
//! `Ω = 2 atan2(a·(b×c), |a||b||c| + (a·b)|c| + (a·c)|b| + (b·c)|a|)` with
//! `a, b, c` the corners relative to `q`.
//!
//! **Far field.** For a node of the hierarchy whose centre `p̃` (the
//! area-weighted centroid of its triangles) is farther from `q` than `β`
//! times the node radius `r` (the largest distance from `p̃` to a corner of
//! its box), the triangles are replaced by the Taylor expansion of the
//! dipole kernel `K(x) = (x − q)/(4π|x − q|³)` about `p̃` with three terms
//! (Barill et al., section 3.2): `w ≈ K·N + Σ_jk J_jk M_jk + ½ Σ_jkl H_jkl
//! T_jkl`, with `N = Σ A_t` the sum of the triangle area vectors,
//! `J = ∂K/∂x`, `H = ∂²K/∂x²`, `M = Σ (c_t − p̃) A_tᵀ` and
//! `T_jkl = Σ n_tj ∫_t d_k d_l dA` (exact moments over each triangle,
//! `d = x − p̃`). The error falls like `β⁻⁴`; [`BETA`] = 2 gives at most
//! 3e-3 next to the surface in the tests and far less away from it.
//!
//! The traversal is depth-first in node order, so the sum, and the result,
//! is the same on every machine and thread count (D93).

use super::Surface;
use crate::geometry::{add, cross, dot, norm, scale, sub};
use crate::query::Bvh;
use rayon::prelude::*;

/// Default ratio of distance to node radius above which the far-field
/// expansion is used.
pub const BETA: f64 = 2.0;

const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;

/// The far-field data of every node of a surface hierarchy.
#[derive(Debug, Clone, Default)]
pub struct Winding {
    /// Area-weighted centroid of the node.
    centre: Vec<[f64; 3]>,
    /// Radius: largest distance from the centre to a corner of the box.
    radius: Vec<f64>,
    /// Sum of the area vectors (half the cross products).
    normal: Vec<[f64; 3]>,
    /// `M = Σ (c_t − p̃) A_tᵀ`, row major.
    moment: Vec<[f64; 9]>,
    /// `T_jkl = Σ_t n_tj ∫_t d_k d_l dA`, `d = x − p̃`, index `9 j + 3 k + l`.
    quad: Vec<[f64; 27]>,
}

/// The signed solid angle of the triangle `t` seen from `q` (steradian):
/// positive when `q` lies behind the triangle (on the side its corners
/// turn clockwise), in (−2π, 2π).
pub fn solid_angle(t: &[[f64; 3]; 3], q: [f64; 3]) -> f64 {
    let a = sub(t[0], q);
    let b = sub(t[1], q);
    let c = sub(t[2], q);
    let (la, lb, lc) = (norm(a), norm(b), norm(c));
    let det = dot(a, cross(b, c));
    let den = la * lb * lc + dot(a, b) * lc + dot(a, c) * lb + dot(b, c) * la;
    2.0 * det.atan2(den)
}

impl Winding {
    /// The far-field data of `bvh` (children after their parent, so one
    /// pass from the last node to the root).
    pub fn build(bvh: &Bvh) -> Self {
        let n = bvh.nodes.len();
        let mut area = vec![0.0f64; n];
        let mut first = vec![[0.0f64; 3]; n]; // Σ a_t c_t
        let mut normal = vec![[0.0f64; 3]; n];
        let mut second = vec![[0.0f64; 9]; n]; // Σ c_t A_tᵀ
        for i in (0..n).rev() {
            let nd = bvh.nodes[i];
            if nd.count > 0 {
                let s = nd.first as usize;
                for &t in &bvh.order[s..s + nd.count as usize] {
                    let tri = &bvh.faces.tris[t as usize];
                    let av = scale(cross(sub(tri[1], tri[0]), sub(tri[2], tri[0])), 0.5);
                    let a = norm(av);
                    let c = scale(add(add(tri[0], tri[1]), tri[2]), 1.0 / 3.0);
                    area[i] += a;
                    first[i] = add(first[i], scale(c, a));
                    normal[i] = add(normal[i], av);
                    for j in 0..3 {
                        for k in 0..3 {
                            second[i][3 * j + k] += c[j] * av[k];
                        }
                    }
                }
            } else {
                for c in [nd.first as usize, nd.first as usize + 1] {
                    area[i] += area[c];
                    first[i] = add(first[i], first[c]);
                    normal[i] = add(normal[i], normal[c]);
                    for m in 0..9 {
                        second[i][m] += second[c][m];
                    }
                }
            }
        }
        let mut centre = vec![[0.0; 3]; n];
        let mut radius = vec![0.0; n];
        let mut moment = vec![[0.0; 9]; n];
        for i in 0..n {
            let bx = bvh.nodes[i].bx;
            let p = if area[i] > 0.0 {
                scale(first[i], 1.0 / area[i])
            } else {
                scale(add(bx.lo, bx.hi), 0.5)
            };
            let mut r2 = 0.0f64;
            for k in 0..8 {
                let corner: [f64; 3] =
                    std::array::from_fn(|d| if k >> d & 1 == 0 { bx.lo[d] } else { bx.hi[d] });
                let e = sub(corner, p);
                r2 = r2.max(dot(e, e));
            }
            centre[i] = p;
            radius[i] = r2.sqrt();
            for j in 0..3 {
                for k in 0..3 {
                    moment[i][3 * j + k] = second[i][3 * j + k] - p[j] * normal[i][k];
                }
            }
        }
        // node ranges in the leaf order (a subtree is a contiguous range)
        let mut range = vec![(usize::MAX, 0usize); n];
        for i in (0..n).rev() {
            let nd = bvh.nodes[i];
            range[i] = if nd.count > 0 {
                (nd.first as usize, nd.first as usize + nd.count as usize)
            } else {
                let (l, r) = (range[nd.first as usize], range[nd.first as usize + 1]);
                (l.0.min(r.0), l.1.max(r.1))
            };
        }
        let quad: Vec<[f64; 27]> = (0..n)
            .into_par_iter()
            .map(|i| {
                let p = centre[i];
                let mut t = [0.0f64; 27];
                let (s, e) = range[i];
                for &k in &bvh.order[s..e] {
                    let tri = &bvh.faces.tris[k as usize];
                    let av = scale(cross(sub(tri[1], tri[0]), sub(tri[2], tri[0])), 0.5);
                    let v = tri.map(|c| sub(c, p));
                    let sum = add(add(v[0], v[1]), v[2]);
                    for kk in 0..3 {
                        for l in 0..3 {
                            let m = (v[0][kk] * v[0][l]
                                + v[1][kk] * v[1][l]
                                + v[2][kk] * v[2][l]
                                + sum[kk] * sum[l])
                                / 12.0;
                            for j in 0..3 {
                                t[9 * j + 3 * kk + l] += av[j] * m;
                            }
                        }
                    }
                }
                t
            })
            .collect();
        Winding {
            centre,
            radius,
            normal,
            moment,
            quad,
        }
    }

    /// The far-field value of node `i` at `q`.
    fn far(&self, i: usize, q: [f64; 3]) -> f64 {
        let r = sub(self.centre[i], q);
        let d2 = dot(r, r);
        let d = d2.sqrt();
        let inv3 = 1.0 / (d2 * d);
        let inv5 = inv3 / d2;
        let mut w = dot(r, self.normal[i]) * inv3;
        let m = &self.moment[i];
        for j in 0..3 {
            for k in 0..3 {
                let jac = if j == k { inv3 } else { 0.0 } - 3.0 * r[j] * r[k] * inv5;
                w += jac * m[3 * j + k];
            }
        }
        // third term: ½ Σ T_jkl ∂²K_j/∂x_k∂x_l
        let t = &self.quad[i];
        let inv7 = inv5 / d2;
        let mut trace = 0.0; // Σ_jk T_jkj r_k (= Σ_jl T_jjl r_l)
        let mut diag = 0.0; // Σ_jk T_jkk r_j
        let mut cubic = 0.0;
        for j in 0..3 {
            for k in 0..3 {
                trace += t[9 * j + 3 * k + j] * r[k];
                diag += t[9 * j + 3 * k + k] * r[j];
                for l in 0..3 {
                    cubic += t[9 * j + 3 * k + l] * r[j] * r[k] * r[l];
                }
            }
        }
        w += 0.5 * (-3.0 * inv5 * (2.0 * trace + diag) + 15.0 * inv7 * cubic);
        w / FOUR_PI
    }
}

impl Surface {
    /// The generalised winding number at `q` by the hierarchy with the
    /// far-field ratio `beta` ([`BETA`] by default; `f64::INFINITY` sums
    /// every triangle exactly).
    pub fn winding_number_with(&self, q: [f64; 3], beta: f64) -> f64 {
        let bvh = &self.bvh;
        if bvh.nodes.is_empty() {
            return 0.0;
        }
        let mut w = 0.0;
        let mut stack: Vec<u32> = Vec::with_capacity(64);
        stack.push(0);
        while let Some(i) = stack.pop() {
            let i = i as usize;
            let nd = bvh.nodes[i];
            let r = sub(self.winding.centre[i], q);
            if norm(r) > beta * self.winding.radius[i] {
                w += self.winding.far(i, q);
                continue;
            }
            if nd.count > 0 {
                let s = nd.first as usize;
                for &t in &bvh.order[s..s + nd.count as usize] {
                    w += solid_angle(&bvh.faces.tris[t as usize], q) / FOUR_PI;
                }
            } else {
                // right child pushed first: the left subtree is summed first
                stack.push(nd.first + 1);
                stack.push(nd.first);
            }
        }
        w
    }

    /// The generalised winding number at `q` ([`BETA`]).
    pub fn winding_number(&self, q: [f64; 3]) -> f64 {
        self.winding_number_with(q, BETA)
    }

    /// The winding number by summing every triangle (the reference of the
    /// tests and the benchmark).
    pub fn winding_number_exact(&self, q: [f64; 3]) -> f64 {
        (0..self.len())
            .map(|t| solid_angle(&self.corners(t), q))
            .sum::<f64>()
            / FOUR_PI
    }

    /// The winding numbers of many points, in parallel (each point summed
    /// in the same order whatever the threads).
    pub fn winding_numbers(&self, points: &[[f64; 3]]) -> Vec<f64> {
        points
            .par_iter()
            .map(|&q| self.winding_number(q))
            .collect()
    }

    /// Whether `q` lies inside the body: winding number above one half.
    /// Within about 1e-3 of the local triangle size of the surface the far
    /// field could decide it wrongly, so there the exact sum over the near
    /// triangles is used (the hierarchy with a larger ratio).
    pub fn inside(&self, q: [f64; 3]) -> bool {
        let w = self.winding_number(q);
        if (w - 0.5).abs() > 0.05 {
            return w > 0.5;
        }
        self.winding_number_with(q, 8.0) > 0.5
    }
}
