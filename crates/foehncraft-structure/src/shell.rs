//! Four-node shells with assumed transverse shear strains: MITC4
//! (Dvorkin & Bathe 1984, *Engineering Computations* 1; Bathe, *Finite
//! Element Procedures*, 2nd ed. 2014, section 5.4).
//!
//! The element is the degenerated continuum: the position is the
//! bilinear interpolation of the nodes plus the thickness coordinate t
//! along the nodal directors V_k, and the displacement adds the turn of
//! each director by the nodal rotation θ_k (global axes, so the element
//! fits the six unknowns per node of the beams):
//!
//!   x = Σ N_k (x_k + t h_k V_k),   u = Σ N_k (u_k + t h_k θ_k × V_k),
//!
//! with h_k half the thickness. The in-plane covariant strains come from
//! these fields at the integration points; the transverse shear strains
//! e_rt and e_st are interpolated from their values at the midpoints of
//! the edges (the tying points A, C at r = 0, s = ±1 and B, D at r = ±1,
//! s = 0), which removes shear locking. The normal strain e_tt is not
//! used and the stress normal to the shell is zero (plane stress), with a
//! shear correction factor κ on the transverse shear. The covariant
//! strains go to an orthonormal frame on the shell (e₁ along g_r, e₃
//! normal) by the contravariant base vectors.
//!
//! The rotation about the director (drilling) has no stiffness in this
//! formulation; a small spring on it (10⁻⁶ of the largest rotational
//! diagonal of the element at that node) keeps K regular without changing
//! the response of a smooth shell to the digits that matter. Integration:
//! 2 × 2 Gauss points on the mid-surface and 2 through the thickness.

/// Material and thickness of a shell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellProps {
    /// Young's modulus (Pa).
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Thickness (m).
    pub thickness: f64,
    /// Shear correction factor (5/6).
    pub kappa: f64,
}

const G2: f64 = 0.577_350_269_189_625_8; // 1/√3
/// Natural coordinates of the four nodes.
const NODES: [(f64, f64); 4] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];

fn shape(r: f64, s: f64) -> ([f64; 4], [f64; 4], [f64; 4]) {
    let mut n = [0.0; 4];
    let mut nr = [0.0; 4];
    let mut ns = [0.0; 4];
    for (k, (ri, si)) in NODES.iter().enumerate() {
        n[k] = 0.25 * (1.0 + ri * r) * (1.0 + si * s);
        nr[k] = 0.25 * ri * (1.0 + si * s);
        ns[k] = 0.25 * si * (1.0 + ri * r);
    }
    (n, nr, ns)
}

use crate::beam::{cross, dot3, norm, scale};

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// −skew(V): the matrix with (θ × V) = −skew(V) θ.
fn turn(v: [f64; 3]) -> [[f64; 3]; 3] {
    [[0.0, v[2], -v[1]], [-v[2], 0.0, v[0]], [v[1], -v[0], 0.0]]
}

/// The geometry of an element: nodes, unit directors, half thickness.
pub struct Geometry {
    /// Node coordinates.
    pub x: [[f64; 3]; 4],
    /// Unit directors (normals) at the nodes.
    pub v: [[f64; 3]; 4],
    /// Half the thickness.
    pub h: f64,
}

/// A 3 × 24 matrix: a vector field at a point in terms of the 24 unknowns.
type Field = [[f64; 24]; 3];

impl Geometry {
    /// Base vectors g_r, g_s, g_t at (r, s, t).
    fn base(&self, r: f64, s: f64, t: f64) -> [[f64; 3]; 3] {
        let (n, nr, ns) = shape(r, s);
        let mut g = [[0.0; 3]; 3];
        for k in 0..4 {
            let p = add(self.x[k], scale(self.v[k], t * self.h));
            for d in 0..3 {
                g[0][d] += nr[k] * p[d];
                g[1][d] += ns[k] * p[d];
                g[2][d] += n[k] * self.h * self.v[k][d];
            }
        }
        g
    }

    /// ∂u/∂r, ∂u/∂s, ∂u/∂t and u at (r, s, t) as fields of the unknowns
    /// (per node u, v, w, θx, θy, θz).
    fn fields(&self, r: f64, s: f64, t: f64) -> [Field; 4] {
        let (n, nr, ns) = shape(r, s);
        let mut f = [[[0.0; 24]; 3]; 4];
        for k in 0..4 {
            let m = turn(self.v[k]);
            // (translation weight, rotation weight) of ∂/∂r, ∂/∂s, ∂/∂t and u
            let w = [
                (nr[k], nr[k] * t * self.h),
                (ns[k], ns[k] * t * self.h),
                (0.0, n[k] * self.h),
                (n[k], n[k] * t * self.h),
            ];
            for (which, &(tr, rot)) in w.iter().enumerate() {
                for d in 0..3 {
                    f[which][d][6 * k + d] += tr;
                    for e in 0..3 {
                        f[which][d][6 * k + 3 + e] += rot * m[d][e];
                    }
                }
            }
        }
        f
    }

    /// The covariant strain e_ij = ½(g_i·u,j + g_j·u,i) as a row of the
    /// unknowns.
    fn covariant(&self, g: &[[f64; 3]; 3], du: &[Field; 4], i: usize, j: usize) -> [f64; 24] {
        let mut row = [0.0; 24];
        for c in 0..24 {
            let mut v = 0.0;
            for d in 0..3 {
                v += g[i][d] * du[j][d][c] + g[j][d] * du[i][d][c];
            }
            row[c] = 0.5 * v;
        }
        row
    }
}

fn inverse(m: [[f64; 3]; 3]) -> ([[f64; 3]; 3], f64) {
    // m has the base vectors as rows: J = mᵀ; returns the contravariant
    // vectors as rows (g^i · g_j = δ_ij) and det
    let a = m;
    let det = dot3(a[0], cross(a[1], a[2]));
    let c0 = cross(a[1], a[2]);
    let c1 = cross(a[2], a[0]);
    let c2 = cross(a[0], a[1]);
    (
        [scale(c0, 1.0 / det), scale(c1, 1.0 / det), scale(c2, 1.0 / det)],
        det,
    )
}

/// Stiffness and consistent mass (24 × 24, row major) of an MITC4 element.
pub fn mitc4(geo: &Geometry, p: &ShellProps) -> (Vec<f64>, Vec<f64>) {
    let e = p.young;
    let nu = p.poisson;
    let c = e / (1.0 - nu * nu);
    let gk = p.kappa * e / (2.0 * (1.0 + nu));
    let gm = e / (2.0 * (1.0 + nu));
    let mut k = vec![0.0; 576];
    let mut m = vec![0.0; 576];
    let pts = [-G2, G2];
    for &t in &pts {
        // the shear strains at the tying points, at this t
        let tie = |r: f64, s: f64, i: usize| {
            let g = geo.base(r, s, t);
            let du = geo.fields(r, s, t);
            geo.covariant(&g, &du, i, 2)
        };
        let ert_a = tie(0.0, 1.0, 0);
        let ert_c = tie(0.0, -1.0, 0);
        let est_b = tie(1.0, 0.0, 1);
        let est_d = tie(-1.0, 0.0, 1);
        for &r in &pts {
            for &s in &pts {
                let g = geo.base(r, s, t);
                let du = geo.fields(r, s, t);
                let (gc, det) = inverse(g);
                let mut cov = [[0.0f64; 24]; 6]; // rr, ss, rs, rt, st
                cov[0] = geo.covariant(&g, &du, 0, 0);
                cov[1] = geo.covariant(&g, &du, 1, 1);
                cov[2] = geo.covariant(&g, &du, 0, 1);
                for q in 0..24 {
                    cov[3][q] = 0.5 * (1.0 + s) * ert_a[q] + 0.5 * (1.0 - s) * ert_c[q];
                    cov[4][q] = 0.5 * (1.0 + r) * est_b[q] + 0.5 * (1.0 - r) * est_d[q];
                }
                // local orthonormal frame
                let e3 = {
                    let n = cross(g[0], g[1]);
                    scale(n, 1.0 / norm(n))
                };
                let e1 = scale(g[0], 1.0 / norm(g[0]));
                let e2 = cross(e3, e1);
                let ea = [e1, e2, e3];
                // Cartesian strain ε_ab = Σ (e_a·g^i)(e_b·g^j) e_ij, e_tt = 0
                let pairs = [(0, 0, 0), (1, 1, 1), (0, 1, 2), (0, 2, 3), (1, 2, 4)];
                let mut eps = [[0.0f64; 24]; 5]; // ε11, ε22, γ12, γ13, γ23
                let outs = [(0, 0, 1.0), (1, 1, 1.0), (0, 1, 2.0), (0, 2, 2.0), (1, 2, 2.0)];
                for (o, &(a, b, f)) in outs.iter().enumerate() {
                    for &(i, j, idx) in &pairs {
                        let w_ij = dot3(ea[a], gc[i]) * dot3(ea[b], gc[j]);
                        let w_ji = dot3(ea[a], gc[j]) * dot3(ea[b], gc[i]);
                        let w = if i == j { w_ij } else { w_ij + w_ji };
                        if w != 0.0 {
                            for q in 0..24 {
                                eps[o][q] += f * w * cov[idx][q];
                            }
                        }
                    }
                }
                let wdet = det.abs(); // weights are 1
                // D: plane stress in-plane, κG transverse shear
                let dmat = |o1: usize, o2: usize| -> f64 {
                    match (o1, o2) {
                        (0, 0) | (1, 1) => c,
                        (0, 1) | (1, 0) => c * nu,
                        (2, 2) => gm,
                        (3, 3) | (4, 4) => gk,
                        _ => 0.0,
                    }
                };
                for a in 0..24 {
                    for b in a..24 {
                        let mut v = 0.0;
                        for o1 in 0..5 {
                            if eps[o1][a] == 0.0 {
                                continue;
                            }
                            for o2 in 0..5 {
                                let dd = dmat(o1, o2);
                                if dd != 0.0 {
                                    v += eps[o1][a] * dd * eps[o2][b];
                                }
                            }
                        }
                        k[a * 24 + b] += wdet * v;
                    }
                }
                let u = &du[3];
                for a in 0..24 {
                    for b in a..24 {
                        let v: f64 = (0..3).map(|d| u[d][a] * u[d][b]).sum();
                        m[a * 24 + b] += wdet * p.density * v;
                    }
                }
            }
        }
    }
    for a in 0..24 {
        for b in 0..a {
            k[a * 24 + b] = k[b * 24 + a];
            m[a * 24 + b] = m[b * 24 + a];
        }
    }
    // drilling springs
    for n in 0..4 {
        let mut big: f64 = 0.0;
        for d in 3..6 {
            big = big.max(k[(6 * n + d) * 25]);
        }
        let kd = 1e-6 * big;
        let v = geo.v[n];
        for a in 0..3 {
            for b in 0..3 {
                k[(6 * n + 3 + a) * 24 + 6 * n + 3 + b] += kd * v[a] * v[b];
            }
        }
    }
    (k, m)
}

/// Work-equivalent nodal forces (24) of a load per area `q` (global axes)
/// on the mid-surface, and the area.
pub fn area_load(geo: &Geometry, q: [f64; 3]) -> ([f64; 24], f64) {
    let mut f = [0.0; 24];
    let mut area = 0.0;
    for &r in &[-G2, G2] {
        for &s in &[-G2, G2] {
            let g = geo.base(r, s, 0.0);
            let da = norm(cross(g[0], g[1]));
            area += da;
            let (n, _, _) = shape(r, s);
            for k in 0..4 {
                for d in 0..3 {
                    f[6 * k + d] += da * n[k] * q[d];
                }
            }
        }
    }
    (f, area)
}

/// The unit normal of the mid-surface at the centre.
pub fn centre_normal(x: &[[f64; 3]; 4]) -> [f64; 3] {
    let d1 = crate::beam::sub(x[2], x[0]);
    let d2 = crate::beam::sub(x[3], x[1]);
    let n = cross(d1, d2);
    scale(n, 1.0 / norm(n))
}

// ------------------------------------------------------------------ MITC3+

/// Seven-point Gauss rule on the triangle (r, s ≥ 0, r + s ≤ 1), exact
/// for polynomials of degree 5 (Dunavant 1985), weights summing to ½.
const TRI7: [(f64, f64, f64); 7] = [
    (1.0 / 3.0, 1.0 / 3.0, 0.225 / 2.0),
    (0.059_715_871_789_769_82, 0.470_142_064_105_115_1, 0.132_394_152_788_506_2 / 2.0),
    (0.470_142_064_105_115_1, 0.059_715_871_789_769_82, 0.132_394_152_788_506_2 / 2.0),
    (0.470_142_064_105_115_1, 0.470_142_064_105_115_1, 0.132_394_152_788_506_2 / 2.0),
    (0.797_426_985_353_087_3, 0.101_286_507_323_456_3, 0.125_939_180_544_827_1 / 2.0),
    (0.101_286_507_323_456_3, 0.797_426_985_353_087_3, 0.125_939_180_544_827_1 / 2.0),
    (0.101_286_507_323_456_3, 0.101_286_507_323_456_3, 0.125_939_180_544_827_1 / 2.0),
];

/// The geometry of a three-node shell element.
pub struct TriGeometry {
    /// Node coordinates.
    pub x: [[f64; 3]; 3],
    /// Unit directors at the nodes.
    pub v: [[f64; 3]; 3],
    /// Half the thickness.
    pub h: f64,
}

/// Unknowns of the element before condensation: 6 per node and the two
/// rotations of the bubble (about V₁ and V₂ of the centroid director).
const NT: usize = 20;
type TriField = [[f64; NT]; 3];

impl TriGeometry {
    /// Centroid director and two unit vectors normal to it.
    fn bubble_frame(&self) -> [[f64; 3]; 3] {
        let vn = {
            let s = add(add(self.v[0], self.v[1]), self.v[2]);
            scale(s, 1.0 / norm(s))
        };
        let pick = if vn[1].abs() < 0.9 { [0.0, 1.0, 0.0] } else { [1.0, 0.0, 0.0] };
        let v1 = {
            let c = cross(pick, vn);
            scale(c, 1.0 / norm(c))
        };
        let v2 = cross(vn, v1);
        [v1, v2, vn]
    }

    /// Base vectors g_r, g_s, g_t (the geometry is linear: the bubble of
    /// the directors cancels with the mean director, Lee & Bathe 2014).
    fn base(&self, t: f64) -> [[f64; 3]; 3] {
        let mut g = [[0.0; 3]; 3];
        for d in 0..3 {
            let p = |k: usize| self.x[k][d] + t * self.h * self.v[k][d];
            g[0][d] = p(1) - p(0);
            g[1][d] = p(2) - p(0);
        }
        let h = [1.0 / 3.0; 3];
        let _ = h;
        g
    }

    fn base_at(&self, r: f64, s: f64, t: f64) -> [[f64; 3]; 3] {
        let mut g = self.base(t);
        let hh = [1.0 - r - s, r, s];
        for d in 0..3 {
            g[2][d] = (0..3).map(|k| hh[k] * self.h * self.v[k][d]).sum();
        }
        g
    }

    /// ∂u/∂r, ∂u/∂s, ∂u/∂t and u at (r, s, t).
    fn fields(&self, r: f64, s: f64, t: f64) -> [TriField; 4] {
        let h = [1.0 - r - s, r, s];
        let hr = [-1.0, 1.0, 0.0];
        let hs = [-1.0, 0.0, 1.0];
        let b = 27.0 * r * s * (1.0 - r - s);
        let br = 27.0 * s * (1.0 - 2.0 * r - s);
        let bs = 27.0 * r * (1.0 - r - 2.0 * s);
        // rotation functions f_i = h_i − f₄/3, f₄ the bubble
        let f = [h[0] - b / 3.0, h[1] - b / 3.0, h[2] - b / 3.0];
        let fr = [hr[0] - br / 3.0, hr[1] - br / 3.0, hr[2] - br / 3.0];
        let fs = [hs[0] - bs / 3.0, hs[1] - bs / 3.0, hs[2] - bs / 3.0];
        let mut out = [[[0.0; NT]; 3]; 4];
        for k in 0..3 {
            let m = turn(self.v[k]);
            let w = [
                (hr[k], fr[k] * t * self.h),
                (hs[k], fs[k] * t * self.h),
                (0.0, f[k] * self.h),
                (h[k], f[k] * t * self.h),
            ];
            for (which, &(tr, rot)) in w.iter().enumerate() {
                for d in 0..3 {
                    out[which][d][6 * k + d] += tr;
                    for e in 0..3 {
                        out[which][d][6 * k + 3 + e] += rot * m[d][e];
                    }
                }
            }
        }
        // the bubble: θ₄ = α V₁ + β V₂ turns the centroid director V₄
        let fr4 = self.bubble_frame();
        let m4 = turn(fr4[2]);
        let w = [(br * t * self.h), (bs * t * self.h), (b * self.h), (b * t * self.h)];
        for (which, &rot) in w.iter().enumerate() {
            for d in 0..3 {
                for (j, axis) in [fr4[0], fr4[1]].iter().enumerate() {
                    // (θ × V₄)_d = Σ_e m4[d][e] θ_e, θ = axis · dof
                    let c: f64 = (0..3).map(|e| m4[d][e] * axis[e]).sum();
                    out[which][d][18 + j] += rot * c;
                }
            }
        }
        out
    }

    fn covariant(g: &[[f64; 3]; 3], du: &[TriField; 4], i: usize, j: usize) -> [f64; NT] {
        let mut row = [0.0; NT];
        for (c, rc) in row.iter_mut().enumerate() {
            let mut v = 0.0;
            for d in 0..3 {
                v += g[i][d] * du[j][d][c] + g[j][d] * du[i][d][c];
            }
            *rc = 0.5 * v;
        }
        row
    }
}

/// Stiffness and consistent mass (18 × 18, row major) of an MITC3+
/// element (Lee, Lee & Bathe 2014, *Computers & Structures* 138): the
/// rotations carry a cubic bubble at the centroid, condensed out here,
/// and the transverse shear strains are tied at six points inside the
/// element:
///
///   ẽ_rt = ⅔(e_rt^B − ½e_st^B) + ⅓(e_rt^C + e_st^C) + ⅓ĉ(3s − 1),
///   ẽ_st = ⅔(e_st^A − ½e_rt^A) + ⅓(e_rt^C + e_st^C) + ⅓ĉ(1 − 3r),
///   ĉ = e_rt^F − e_rt^D − e_st^F + e_st^E,
///
/// A = (1/6, 2/3), B = (2/3, 1/6), C = (1/6, 1/6), D = (1/3 + d, 1/3 − 2d),
/// E = (1/3 − 2d, 1/3 + d), F = (1/3 + d, 1/3 + d), d = 10⁻⁴.
pub fn mitc3plus(geo: &TriGeometry, p: &ShellProps) -> (Vec<f64>, Vec<f64>) {
    let e = p.young;
    let nu = p.poisson;
    let c = e / (1.0 - nu * nu);
    let gk = p.kappa * e / (2.0 * (1.0 + nu));
    let gm = e / (2.0 * (1.0 + nu));
    let mut k = vec![0.0; NT * NT];
    let mut m = vec![0.0; NT * NT];
    let dd = 1e-4;
    let tying = [
        (1.0 / 6.0, 2.0 / 3.0),
        (2.0 / 3.0, 1.0 / 6.0),
        (1.0 / 6.0, 1.0 / 6.0),
        (1.0 / 3.0 + dd, 1.0 / 3.0 - 2.0 * dd),
        (1.0 / 3.0 - 2.0 * dd, 1.0 / 3.0 + dd),
        (1.0 / 3.0 + dd, 1.0 / 3.0 + dd),
    ];
    for &t in &[-G2, G2] {
        let at = |r: f64, s: f64| {
            let g = geo.base_at(r, s, t);
            let du = geo.fields(r, s, t);
            (
                TriGeometry::covariant(&g, &du, 0, 2),
                TriGeometry::covariant(&g, &du, 1, 2),
            )
        };
        let tp: Vec<([f64; NT], [f64; NT])> = tying.iter().map(|&(r, s)| at(r, s)).collect();
        let (a, b, cc, d, ee, f) = (&tp[0], &tp[1], &tp[2], &tp[3], &tp[4], &tp[5]);
        for &(r, s, w) in &TRI7 {
            let g = geo.base_at(r, s, t);
            let du = geo.fields(r, s, t);
            let (gc, det) = inverse(g);
            let mut cov = [[0.0f64; NT]; 5];
            cov[0] = TriGeometry::covariant(&g, &du, 0, 0);
            cov[1] = TriGeometry::covariant(&g, &du, 1, 1);
            cov[2] = TriGeometry::covariant(&g, &du, 0, 1);
            for q in 0..NT {
                let chat = f.0[q] - d.0[q] - f.1[q] + ee.1[q];
                let common = (cc.0[q] + cc.1[q]) / 3.0;
                cov[3][q] = 2.0 / 3.0 * (b.0[q] - 0.5 * b.1[q]) + common + chat / 3.0 * (3.0 * s - 1.0);
                cov[4][q] = 2.0 / 3.0 * (a.1[q] - 0.5 * a.0[q]) + common + chat / 3.0 * (1.0 - 3.0 * r);
            }
            let e3 = {
                let n = cross(g[0], g[1]);
                scale(n, 1.0 / norm(n))
            };
            let e1 = scale(g[0], 1.0 / norm(g[0]));
            let e2 = cross(e3, e1);
            let ea = [e1, e2, e3];
            let pairs = [(0, 0, 0), (1, 1, 1), (0, 1, 2), (0, 2, 3), (1, 2, 4)];
            let outs = [(0, 0, 1.0), (1, 1, 1.0), (0, 1, 2.0), (0, 2, 2.0), (1, 2, 2.0)];
            let mut eps = [[0.0f64; NT]; 5];
            for (o, &(a1, b1, fac)) in outs.iter().enumerate() {
                for &(i, j, idx) in &pairs {
                    let wij = dot3(ea[a1], gc[i]) * dot3(ea[b1], gc[j]);
                    let wji = dot3(ea[a1], gc[j]) * dot3(ea[b1], gc[i]);
                    let wt = if i == j { wij } else { wij + wji };
                    if wt != 0.0 {
                        for q in 0..NT {
                            eps[o][q] += fac * wt * cov[idx][q];
                        }
                    }
                }
            }
            let wdet = w * det.abs();
            let dmat = [[c, c * nu, 0.0, 0.0, 0.0], [c * nu, c, 0.0, 0.0, 0.0], [0.0, 0.0, gm, 0.0, 0.0], [0.0, 0.0, 0.0, gk, 0.0], [0.0, 0.0, 0.0, 0.0, gk]];
            for a1 in 0..NT {
                for b1 in a1..NT {
                    let mut v = 0.0;
                    for o1 in 0..5 {
                        if eps[o1][a1] == 0.0 {
                            continue;
                        }
                        for o2 in 0..5 {
                            v += eps[o1][a1] * dmat[o1][o2] * eps[o2][b1];
                        }
                    }
                    k[a1 * NT + b1] += wdet * v;
                    let u = &du[3];
                    let mv: f64 = (0..3).map(|d| u[d][a1] * u[d][b1]).sum();
                    m[a1 * NT + b1] += wdet * p.density * mv;
                }
            }
        }
    }
    for a1 in 0..NT {
        for b1 in 0..a1 {
            k[a1 * NT + b1] = k[b1 * NT + a1];
            m[a1 * NT + b1] = m[b1 * NT + a1];
        }
    }
    // condense the bubble: T = [I; −K_bb⁻¹ K_ba], K* = Tᵀ K T, M* = Tᵀ M T
    let (kb00, kb01, kb11) = (k[18 * NT + 18], k[18 * NT + 19], k[19 * NT + 19]);
    let det = kb00 * kb11 - kb01 * kb01;
    let inv = [[kb11 / det, -kb01 / det], [-kb01 / det, kb00 / det]];
    let mut tmat = vec![0.0; NT * 18];
    for i in 0..18 {
        tmat[i * 18 + i] = 1.0;
        for bi in 0..2 {
            let v: f64 = (0..2).map(|bj| -inv[bi][bj] * k[(18 + bj) * NT + i]).sum();
            tmat[(18 + bi) * 18 + i] = v;
        }
    }
    let reduce = |a: &[f64]| -> Vec<f64> {
        let mut at = vec![0.0; NT * 18];
        for i in 0..NT {
            for j in 0..18 {
                at[i * 18 + j] = (0..NT).map(|l| a[i * NT + l] * tmat[l * 18 + j]).sum();
            }
        }
        let mut out = vec![0.0; 18 * 18];
        for i in 0..18 {
            for j in 0..18 {
                out[i * 18 + j] = (0..NT).map(|l| tmat[l * 18 + i] * at[l * 18 + j]).sum();
            }
        }
        out
    };
    let mut kr = reduce(&k);
    let mr = reduce(&m);
    for n in 0..3 {
        let mut big: f64 = 0.0;
        for d in 3..6 {
            big = big.max(kr[(6 * n + d) * 19]);
        }
        let kd = 1e-6 * big;
        let v = geo.v[n];
        for a1 in 0..3 {
            for b1 in 0..3 {
                kr[(6 * n + 3 + a1) * 18 + 6 * n + 3 + b1] += kd * v[a1] * v[b1];
            }
        }
    }
    (kr, mr)
}

/// Work-equivalent nodal forces (18) of a load per area `q` on a
/// three-node element, and its area.
pub fn tri_area_load(x: &[[f64; 3]; 3], q: [f64; 3]) -> ([f64; 18], f64) {
    let area = 0.5 * norm(cross(crate::beam::sub(x[1], x[0]), crate::beam::sub(x[2], x[0])));
    let mut f = [0.0; 18];
    for k in 0..3 {
        for d in 0..3 {
            f[6 * k + d] = area / 3.0 * q[d];
        }
    }
    (f, area)
}
