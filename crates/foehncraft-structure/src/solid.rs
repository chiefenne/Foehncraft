//! Quadratic solid elements: the ten-node tetrahedron and the twenty-node
//! hexahedron (serendipity), isoparametric, three translations per node.
//! Zienkiewicz, Taylor & Zhu, *The Finite Element Method: Its Basis and
//! Fundamentals*, 7th ed. (2013), chapter 6; Bathe, *Finite Element
//! Procedures*, 2nd ed. (2014), section 5.3.
//!
//! - Stiffness K = ∫ Bᵀ D B dV and consistent mass M = ∫ ρ Nᵀ N dV with
//!   the isotropic elasticity matrix D (strains in Voigt order xx, yy, zz,
//!   and the engineering shears xy, yz, zx).
//! - Thermal strain ε_th = α ΔT on the normal strains; its load is
//!   f = ∫ Bᵀ D ε_th dV.
//! - Integration: the hexahedron with 3 × 3 × 3 Gauss points (full: no
//!   spurious zero-energy modes); the tetrahedron with the 14-point rule of
//!   degree 5 (Walkington 2000, *Quadrature on simplices of arbitrary
//!   dimension*, also in Jaśkowiec & Sukumar 2020), exact for the mass of a
//!   straight-sided element.
//! - Node order: the VTK order (corners, then the edges); [`GMSH_TET10`]
//!   and [`GMSH_HEX20`] map Gmsh's order to it.
//! - Faces carry pressures and forces per area: the six-node triangle with
//!   the 7-point rule of degree 5 (Dunavant 1985) and the eight-node
//!   quadrilateral with 3 × 3 Gauss points.

#![allow(clippy::needless_range_loop)] // index loops over tensor components read as the formulas

/// The kind of a solid element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// Ten-node tetrahedron.
    Tet10,
    /// Twenty-node hexahedron.
    Hex20,
}

/// Material of a solid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidProps {
    /// Young's modulus (Pa).
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Density (kg/m³).
    pub density: f64,
}

/// Most nodes of an element.
pub const MAX_NODES: usize = 20;

/// Corners and edges (as corner pairs) of the ten-node tetrahedron.
pub const TET_EDGES: [[usize; 2]; 6] = [[0, 1], [1, 2], [0, 2], [0, 3], [1, 3], [2, 3]];

/// Edges of the twenty-node hexahedron (VTK order).
pub const HEX_EDGES: [[usize; 2]; 12] = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 4],
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7],
];

/// Gmsh node of each node of the ten-node tetrahedron (Gmsh type 11).
/// Gmsh numbers the edges 01, 12, 20, 30, 32, 31.
pub const GMSH_TET10: [usize; 10] = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8];

/// Gmsh node of each node of the twenty-node hexahedron (Gmsh type 17).
/// Gmsh numbers the edges 01, 03, 04, 12, 15, 23, 26, 37, 45, 47, 56, 67.
pub const GMSH_HEX20: [usize; 20] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15,
];

/// The faces of the elements, counter-clockwise seen from outside: the
/// corners, then the edge nodes (a six-node triangle or an eight-node
/// quadrilateral).
pub const TET_FACES: [[usize; 6]; 4] = [
    [0, 2, 1, 6, 5, 4],
    [0, 1, 3, 4, 8, 7],
    [1, 2, 3, 5, 9, 8],
    [0, 3, 2, 7, 9, 6],
];

/// The faces of the hexahedron (see [`TET_FACES`]).
pub const HEX_FACES: [[usize; 8]; 6] = [
    [0, 3, 2, 1, 11, 10, 9, 8],
    [4, 5, 6, 7, 12, 13, 14, 15],
    [0, 1, 5, 4, 8, 17, 12, 16],
    [2, 3, 7, 6, 10, 19, 14, 18],
    [0, 4, 7, 3, 16, 15, 19, 11],
    [1, 2, 6, 5, 9, 18, 13, 17],
];

impl Kind {
    /// Nodes per element.
    pub fn nodes(self) -> usize {
        match self {
            Kind::Tet10 => 10,
            Kind::Hex20 => 20,
        }
    }

    /// Corners per element.
    pub fn corners(self) -> usize {
        match self {
            Kind::Tet10 => 4,
            Kind::Hex20 => 8,
        }
    }

    /// The faces (node lists, see [`TET_FACES`]).
    pub fn faces(self) -> Vec<&'static [usize]> {
        match self {
            Kind::Tet10 => TET_FACES.iter().map(|f| &f[..]).collect(),
            Kind::Hex20 => HEX_FACES.iter().map(|f| &f[..]).collect(),
        }
    }

    /// The VTK cell type (24 quadratic tetrahedron, 25 quadratic
    /// hexahedron).
    pub fn vtk_type(self) -> u8 {
        match self {
            Kind::Tet10 => 24,
            Kind::Hex20 => 25,
        }
    }

    /// The natural coordinates of the nodes.
    pub fn natural_nodes(self) -> Vec<[f64; 3]> {
        match self {
            Kind::Tet10 => {
                let c = [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ];
                let mut v = c.to_vec();
                for [a, b] in TET_EDGES {
                    v.push(mid(c[a], c[b]));
                }
                v
            }
            Kind::Hex20 => {
                let c = HEX_CORNERS;
                let mut v = c.to_vec();
                for [a, b] in HEX_EDGES {
                    v.push(mid(c[a], c[b]));
                }
                v
            }
        }
    }
}

const HEX_CORNERS: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
];

fn mid(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        0.5 * (a[0] + b[0]),
        0.5 * (a[1] + b[1]),
        0.5 * (a[2] + b[2]),
    ]
}

/// Shape functions and their derivatives with respect to the natural
/// coordinates at `p`.
pub fn shape(kind: Kind, p: [f64; 3]) -> ([f64; MAX_NODES], [[f64; 3]; MAX_NODES]) {
    let mut n = [0.0; MAX_NODES];
    let mut d = [[0.0; 3]; MAX_NODES];
    match kind {
        Kind::Tet10 => {
            let l = [1.0 - p[0] - p[1] - p[2], p[0], p[1], p[2]];
            let dl = [
                [-1.0, -1.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ];
            for i in 0..4 {
                n[i] = l[i] * (2.0 * l[i] - 1.0);
                for k in 0..3 {
                    d[i][k] = (4.0 * l[i] - 1.0) * dl[i][k];
                }
            }
            for (e, [a, b]) in TET_EDGES.iter().enumerate() {
                n[4 + e] = 4.0 * l[*a] * l[*b];
                for k in 0..3 {
                    d[4 + e][k] = 4.0 * (l[*b] * dl[*a][k] + l[*a] * dl[*b][k]);
                }
            }
        }
        Kind::Hex20 => {
            let nodes = Kind::Hex20.natural_nodes();
            for (i, c) in nodes.iter().enumerate() {
                let [a, b, cc] = *c;
                let (x, y, z) = (p[0], p[1], p[2]);
                if i < 8 {
                    let (pp, q, r) = (1.0 + a * x, 1.0 + b * y, 1.0 + cc * z);
                    let s = a * x + b * y + cc * z - 2.0;
                    n[i] = 0.125 * pp * q * r * s;
                    d[i] = [
                        0.125 * a * q * r * (s + pp),
                        0.125 * b * pp * r * (s + q),
                        0.125 * cc * pp * q * (s + r),
                    ];
                } else if a == 0.0 {
                    let (q, r) = (1.0 + b * y, 1.0 + cc * z);
                    n[i] = 0.25 * (1.0 - x * x) * q * r;
                    d[i] = [
                        -0.5 * x * q * r,
                        0.25 * (1.0 - x * x) * b * r,
                        0.25 * (1.0 - x * x) * q * cc,
                    ];
                } else if b == 0.0 {
                    let (pp, r) = (1.0 + a * x, 1.0 + cc * z);
                    n[i] = 0.25 * (1.0 - y * y) * pp * r;
                    d[i] = [
                        0.25 * (1.0 - y * y) * a * r,
                        -0.5 * y * pp * r,
                        0.25 * (1.0 - y * y) * pp * cc,
                    ];
                } else {
                    let (pp, q) = (1.0 + a * x, 1.0 + b * y);
                    n[i] = 0.25 * (1.0 - z * z) * pp * q;
                    d[i] = [
                        0.25 * (1.0 - z * z) * a * q,
                        0.25 * (1.0 - z * z) * pp * b,
                        -0.5 * z * pp * q,
                    ];
                }
            }
        }
    }
    (n, d)
}

/// The integration points and weights of an element.
pub fn rule(kind: Kind) -> Vec<([f64; 3], f64)> {
    match kind {
        Kind::Tet10 => {
            let mut v = Vec::with_capacity(14);
            for (a, w) in [
                (0.092_735_250_310_891_2, 0.012_248_840_519_393_66),
                (0.310_885_919_263_300_6, 0.018_781_320_953_002_64),
            ] {
                let b = 1.0 - 3.0 * a;
                for p in [[a, a, a], [b, a, a], [a, b, a], [a, a, b]] {
                    v.push((p, w));
                }
            }
            let (a, b, w) = (
                0.454_496_295_874_350_4,
                0.045_503_704_125_649_6,
                0.007_091_003_462_846_911,
            );
            for p in [
                [a, a, b],
                [a, b, a],
                [b, a, a],
                [a, b, b],
                [b, a, b],
                [b, b, a],
            ] {
                v.push((p, w));
            }
            v
        }
        Kind::Hex20 => {
            let g = [-(0.6f64.sqrt()), 0.0, 0.6f64.sqrt()];
            let w = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
            let mut v = Vec::with_capacity(27);
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        v.push(([g[i], g[j], g[k]], w[i] * w[j] * w[k]));
                    }
                }
            }
            v
        }
    }
}

/// The elasticity matrix (Voigt order xx, yy, zz, xy, yz, zx with
/// engineering shears).
pub fn elasticity(p: &SolidProps) -> [[f64; 6]; 6] {
    let (e, nu) = (p.young, p.poisson);
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let mut d = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            d[i][j] = lam;
        }
        d[i][i] = lam + 2.0 * mu;
        d[3 + i][3 + i] = mu;
    }
    d
}

/// Shape functions, their gradients in global axes and the Jacobian
/// determinant at a point.
type Gradients = ([f64; MAX_NODES], [[f64; 3]; MAX_NODES], f64);

/// Stiffness, mass (row-major) and thermal load of an element.
pub type ElementMatrices = (Vec<f64>, Vec<f64>, Vec<f64>);

/// Gradients of the shape functions in global axes and the Jacobian
/// determinant at a point.
fn gradients(kind: Kind, x: &[[f64; 3]], p: [f64; 3]) -> Result<Gradients, String> {
    let nn = kind.nodes();
    let (n, dn) = shape(kind, p);
    // J[i][k] = ∂x_k/∂ξ_i
    let mut j = [[0.0; 3]; 3];
    for a in 0..nn {
        for i in 0..3 {
            for k in 0..3 {
                j[i][k] += dn[a][i] * x[a][k];
            }
        }
    }
    let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
        - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
        + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
    if det.is_nan() || det <= 0.0 {
        return Err(format!(
            "a {} has a Jacobian determinant ≤ 0 (inverted or badly distorted)",
            match kind {
                Kind::Tet10 => "ten-node tetrahedron",
                Kind::Hex20 => "twenty-node hexahedron",
            }
        ));
    }
    let inv = [
        [
            (j[1][1] * j[2][2] - j[1][2] * j[2][1]) / det,
            (j[0][2] * j[2][1] - j[0][1] * j[2][2]) / det,
            (j[0][1] * j[1][2] - j[0][2] * j[1][1]) / det,
        ],
        [
            (j[1][2] * j[2][0] - j[1][0] * j[2][2]) / det,
            (j[0][0] * j[2][2] - j[0][2] * j[2][0]) / det,
            (j[0][2] * j[1][0] - j[0][0] * j[1][2]) / det,
        ],
        [
            (j[1][0] * j[2][1] - j[1][1] * j[2][0]) / det,
            (j[0][1] * j[2][0] - j[0][0] * j[2][1]) / det,
            (j[0][0] * j[1][1] - j[0][1] * j[1][0]) / det,
        ],
    ];
    // ∂N/∂x_k = Σ_i (J⁻¹)[k][i] ∂N/∂ξ_i
    let mut g = [[0.0; 3]; MAX_NODES];
    for a in 0..nn {
        for k in 0..3 {
            g[a][k] = (0..3).map(|i| inv[k][i] * dn[a][i]).sum();
        }
    }
    Ok((n, g, det))
}

/// The strain-displacement rows of node `a` (6 × 3).
fn b_node(g: [f64; 3]) -> [[f64; 3]; 6] {
    [
        [g[0], 0.0, 0.0],
        [0.0, g[1], 0.0],
        [0.0, 0.0, g[2]],
        [g[1], g[0], 0.0],
        [0.0, g[2], g[1]],
        [g[2], 0.0, g[0]],
    ]
}

/// Stiffness, consistent mass and thermal load of an element with the
/// thermal strain α ΔT at its nodes, `strain` (one per node, interpolated
/// by the shape functions; empty: none), row-major, 3 unknowns per node
/// in node order.
pub fn matrices(
    kind: Kind,
    x: &[[f64; 3]],
    p: &SolidProps,
    strain: &[f64],
) -> Result<ElementMatrices, String> {
    let nn = kind.nodes();
    let nd = 3 * nn;
    let d = elasticity(p);
    // D ε_th per unit strain: the thermal stress of a held element
    let mut sth1 = [0.0; 6];
    for i in 0..6 {
        sth1[i] = (0..3).map(|j| d[i][j]).sum();
    }
    let mut k = vec![0.0; nd * nd];
    let mut m = vec![0.0; nd * nd];
    let mut f = vec![0.0; nd];
    for (pt, w) in rule(kind) {
        let (n, g, det) = gradients(kind, x, pt)?;
        let wd = w * det;
        let eth: f64 = strain.iter().zip(&n[..nn]).map(|(a, b)| a * b).sum();
        let sth = sth1.map(|v| v * eth);
        let bs: Vec<[[f64; 3]; 6]> = (0..nn).map(|a| b_node(g[a])).collect();
        // D B per node
        let dbs: Vec<[[f64; 3]; 6]> = bs
            .iter()
            .map(|b| {
                let mut r = [[0.0; 3]; 6];
                for i in 0..6 {
                    for c in 0..3 {
                        r[i][c] = (0..6).map(|l| d[i][l] * b[l][c]).sum();
                    }
                }
                r
            })
            .collect();
        for a in 0..nn {
            for c in 0..3 {
                let row = 3 * a + c;
                f[row] += wd * (0..6).map(|i| bs[a][i][c] * sth[i]).sum::<f64>();
                for b in 0..nn {
                    let mab = wd * p.density * n[a] * n[b];
                    m[row * nd + 3 * b + c] += mab;
                    for e in 0..3 {
                        let v: f64 = (0..6).map(|i| bs[a][i][c] * dbs[b][i][e]).sum();
                        k[row * nd + 3 * b + e] += wd * v;
                    }
                }
            }
        }
    }
    Ok((k, m, f))
}

/// The volume of an element.
pub fn volume(kind: Kind, x: &[[f64; 3]]) -> Result<f64, String> {
    let mut v = 0.0;
    for (pt, w) in rule(kind) {
        v += w * gradients(kind, x, pt)?.2;
    }
    Ok(v)
}

/// The stresses (Voigt order) at the nodes of an element with the nodal
/// displacements `u` (3 per node) and the thermal strain at the nodes
/// `strain` (one per node; empty: none).
pub fn nodal_stresses(
    kind: Kind,
    x: &[[f64; 3]],
    p: &SolidProps,
    u: &[f64],
    strain: &[f64],
) -> Result<Vec<[f64; 6]>, String> {
    let nn = kind.nodes();
    let d = elasticity(p);
    let mut out = Vec::with_capacity(nn);
    for (node, pt) in kind.natural_nodes().into_iter().enumerate() {
        let (_, g, _) = gradients(kind, x, pt)?;
        let eth = strain.get(node).copied().unwrap_or(0.0);
        let mut eps = [0.0; 6];
        for a in 0..nn {
            let b = b_node(g[a]);
            for i in 0..6 {
                for c in 0..3 {
                    eps[i] += b[i][c] * u[3 * a + c];
                }
            }
        }
        for i in 0..3 {
            eps[i] -= eth;
        }
        let mut s = [0.0; 6];
        for i in 0..6 {
            s[i] = (0..6).map(|l| d[i][l] * eps[l]).sum();
        }
        out.push(s);
    }
    Ok(out)
}

/// The von Mises stress of a stress (Voigt order).
pub fn von_mises(s: &[f64; 6]) -> f64 {
    let (a, b, c) = (s[0] - s[1], s[1] - s[2], s[2] - s[0]);
    (0.5 * (a * a + b * b + c * c) + 3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])).sqrt()
}

/// Shape functions and natural derivatives of a face: the six-node
/// triangle (corners (0,0), (1,0), (0,1)) or the eight-node quadrilateral
/// (corners (±1, ±1)).
fn face_shape(quad: bool, r: f64, s: f64) -> ([f64; 8], [[f64; 2]; 8]) {
    let mut n = [0.0; 8];
    let mut d = [[0.0; 2]; 8];
    if quad {
        let c = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
        for i in 0..4 {
            let (a, b) = (c[i][0], c[i][1]);
            let (p, q) = (1.0 + a * r, 1.0 + b * s);
            let t = a * r + b * s - 1.0;
            n[i] = 0.25 * p * q * t;
            d[i] = [0.25 * a * q * (t + p), 0.25 * b * p * (t + q)];
        }
        for (e, [i, j]) in [[0, 1], [1, 2], [2, 3], [3, 0]].iter().enumerate() {
            let a = 0.5 * (c[*i][0] + c[*j][0]);
            let b = 0.5 * (c[*i][1] + c[*j][1]);
            if a == 0.0 {
                n[4 + e] = 0.5 * (1.0 - r * r) * (1.0 + b * s);
                d[4 + e] = [-r * (1.0 + b * s), 0.5 * (1.0 - r * r) * b];
            } else {
                n[4 + e] = 0.5 * (1.0 - s * s) * (1.0 + a * r);
                d[4 + e] = [0.5 * (1.0 - s * s) * a, -s * (1.0 + a * r)];
            }
        }
    } else {
        let l = [1.0 - r - s, r, s];
        let dl = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
        for i in 0..3 {
            n[i] = l[i] * (2.0 * l[i] - 1.0);
            d[i] = [(4.0 * l[i] - 1.0) * dl[i][0], (4.0 * l[i] - 1.0) * dl[i][1]];
        }
        for (e, [a, b]) in [[0, 1], [1, 2], [2, 0]].iter().enumerate() {
            n[3 + e] = 4.0 * l[*a] * l[*b];
            for k in 0..2 {
                d[3 + e][k] = 4.0 * (l[*b] * dl[*a][k] + l[*a] * dl[*b][k]);
            }
        }
    }
    (n, d)
}

/// The integration points of a face.
fn face_rule(quad: bool) -> Vec<([f64; 2], f64)> {
    if quad {
        let g = [-(0.6f64.sqrt()), 0.0, 0.6f64.sqrt()];
        let w = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
        let mut v = Vec::with_capacity(9);
        for i in 0..3 {
            for j in 0..3 {
                v.push(([g[i], g[j]], w[i] * w[j]));
            }
        }
        v
    } else {
        // Dunavant 1985, degree 5, weights summing to the area 1/2
        let mut v = vec![([1.0 / 3.0, 1.0 / 3.0], 0.5 * 0.225)];
        for (a, b, w) in [
            (
                0.059_715_871_789_770,
                0.470_142_064_105_115,
                0.132_394_152_788_506,
            ),
            (
                0.797_426_985_353_087,
                0.101_286_507_323_456,
                0.125_939_180_544_827,
            ),
        ] {
            for p in [[b, b], [a, b], [b, a]] {
                v.push((p, 0.5 * w));
            }
        }
        v
    }
}

/// The nodal forces of a face (6 or 8 nodes, counter-clockwise seen from
/// outside) under a force per area `q` (global axes) and a pressure `p`
/// pushing into the solid (positive against the outward normal); and the
/// face area.
pub fn face_load(x: &[[f64; 3]], q: [f64; 3], p: f64) -> (Vec<f64>, f64) {
    let quad = x.len() == 8;
    let nn = x.len();
    let mut f = vec![0.0; 3 * nn];
    let mut area = 0.0;
    for ([r, s], w) in face_rule(quad) {
        let (n, d) = face_shape(quad, r, s);
        let mut tr = [0.0; 3];
        let mut ts = [0.0; 3];
        for a in 0..nn {
            for k in 0..3 {
                tr[k] += d[a][0] * x[a][k];
                ts[k] += d[a][1] * x[a][k];
            }
        }
        let c = crate::beam::cross(tr, ts);
        let j = crate::beam::norm(c);
        area += w * j;
        for a in 0..nn {
            for k in 0..3 {
                f[3 * a + k] += w * n[a] * (q[k] * j - p * c[k]);
            }
        }
    }
    (f, area)
}

/// The outward normal (not unit) of a face at its centre.
pub fn face_normal(x: &[[f64; 3]]) -> [f64; 3] {
    let quad = x.len() == 8;
    let (r, s) = if quad {
        (0.0, 0.0)
    } else {
        (1.0 / 3.0, 1.0 / 3.0)
    };
    let (_, d) = face_shape(quad, r, s);
    let mut tr = [0.0; 3];
    let mut ts = [0.0; 3];
    for a in 0..x.len() {
        for k in 0..3 {
            tr[k] += d[a][0] * x[a][k];
            ts[k] += d[a][1] * x[a][k];
        }
    }
    crate::beam::cross(tr, ts)
}

/// A face with its node order turned over (the normal flips).
pub fn flip_face(f: &[usize]) -> Vec<usize> {
    if f.len() == 8 {
        vec![f[0], f[3], f[2], f[1], f[7], f[6], f[5], f[4]]
    } else {
        vec![f[0], f[2], f[1], f[5], f[4], f[3]]
    }
}
