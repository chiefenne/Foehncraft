//! Two-node beams in space: Euler–Bernoulli and Timoshenko.
//!
//! Each plane of bending uses the interdependent interpolation of the
//! shear-deformable beam (Friedman & Kosmatka 1993, *Computers &
//! Structures* 47; Reddy 1997, *Comput. Methods Appl. Mech. Engrg.* 149):
//! a cubic for the deflection and a quadratic for the rotation, both
//! depending on Φ = 12EI/(κGAL²), so the shear strain is constant along
//! the element and the element has no shear locking. The stiffness so
//! obtained is the exact stiffness of the Timoshenko beam (nodal values
//! exact for end loads), and with Φ = 0 the functions are the cubic
//! Hermite polynomials of the Euler–Bernoulli beam. Stiffness, consistent
//! mass and the work-equivalent loads are integrated with four Gauss
//! points, exact for these polynomials. Axial and torsional terms are
//! linear.
//!
//! Local axes: x along the beam from node 1 to node 2, y and z the
//! principal axes of the section. Unknowns per node: u, v, w, θx, θy, θz.
//! Bending in the x–y plane uses I_z and θz = dv/dx; in the x–z plane I_y
//! and θy = −dw/dx.

/// Beam theory.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Theory {
    /// Euler–Bernoulli: no shear strain, no rotary inertia.
    EulerBernoulli,
    /// Timoshenko: shear strain and rotary inertia.
    Timoshenko,
}

/// Section and material of a beam, in its local axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamProps {
    /// Young's modulus E (Pa).
    pub young: f64,
    /// Shear modulus G (Pa).
    pub shear: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Area A (m²).
    pub area: f64,
    /// Second moment about local y, I_y (bending in x–z) (m⁴).
    pub iy: f64,
    /// Second moment about local z, I_z (bending in x–y) (m⁴).
    pub iz: f64,
    /// Torsion constant J (m⁴).
    pub torsion: f64,
    /// Shear coefficient for shear along local y (κ_y).
    pub kappa_y: f64,
    /// Shear coefficient for shear along local z (κ_z).
    pub kappa_z: f64,
    /// Beam theory.
    pub theory: Theory,
}

/// Gauss–Legendre points and weights on [0, 1] (four points).
const GAUSS4: [(f64, f64); 4] = [
    (0.069_431_844_202_973_71, 0.173_927_422_568_726_93),
    (0.330_009_478_207_571_87, 0.326_072_577_431_273_07),
    (0.669_990_521_792_428_1, 0.326_072_577_431_273_07),
    (0.930_568_155_797_026_3, 0.173_927_422_568_726_93),
];

/// Shape functions of one bending plane at ξ ∈ [0, 1]: deflection N_w,
/// its derivative dN_w/dx, rotation N_θ and dN_θ/dx, for the unknowns
/// (w₁, θ₁, w₂, θ₂) with θ = dw/dx when Φ = 0.
fn shapes(xi: f64, len: f64, phi: f64) -> [[f64; 4]; 4] {
    let mu = 1.0 / (1.0 + phi);
    let (x2, x3) = (xi * xi, xi * xi * xi);
    let nw = [
        mu * (2.0 * x3 - 3.0 * x2 - phi * xi + 1.0 + phi),
        mu * len * (x3 - (2.0 + 0.5 * phi) * x2 + (1.0 + 0.5 * phi) * xi),
        -mu * (2.0 * x3 - 3.0 * x2 - phi * xi),
        mu * len * (x3 - (1.0 - 0.5 * phi) * x2 - 0.5 * phi * xi),
    ];
    let dnw = [
        mu * (6.0 * x2 - 6.0 * xi - phi) / len,
        mu * (3.0 * x2 - (4.0 + phi) * xi + 1.0 + 0.5 * phi),
        -mu * (6.0 * x2 - 6.0 * xi - phi) / len,
        mu * (3.0 * x2 - (2.0 - phi) * xi - 0.5 * phi),
    ];
    let nt = [
        6.0 * mu / len * (x2 - xi),
        mu * (3.0 * x2 - (4.0 + phi) * xi + 1.0 + phi),
        -6.0 * mu / len * (x2 - xi),
        mu * (3.0 * x2 - (2.0 - phi) * xi),
    ];
    let dnt = [
        6.0 * mu / (len * len) * (2.0 * xi - 1.0),
        mu * (6.0 * xi - (4.0 + phi)) / len,
        -6.0 * mu / (len * len) * (2.0 * xi - 1.0),
        mu * (6.0 * xi - (2.0 - phi)) / len,
    ];
    [nw, dnw, nt, dnt]
}

/// Stiffness and mass (4 × 4, row major) of one bending plane, and the
/// load vector of a unit distributed transverse load.
fn plane(
    len: f64,
    ei: f64,
    kga: Option<f64>,
    rho_a: f64,
    rho_i: f64,
) -> ([f64; 16], [f64; 16], [f64; 4]) {
    let phi = kga.map_or(0.0, |s| 12.0 * ei / (s * len * len));
    let mut k = [0.0; 16];
    let mut m = [0.0; 16];
    let mut f = [0.0; 4];
    for (xi, w) in GAUSS4 {
        let [nw, dnw, nt, dnt] = shapes(xi, len, phi);
        let wl = w * len;
        for a in 0..4 {
            f[a] += wl * nw[a];
            for b in 0..4 {
                let mut kab = ei * dnt[a] * dnt[b];
                if let Some(s) = kga {
                    kab += s * (dnw[a] - nt[a]) * (dnw[b] - nt[b]);
                }
                k[a * 4 + b] += wl * kab;
                m[a * 4 + b] += wl * (rho_a * nw[a] * nw[b] + rho_i * nt[a] * nt[b]);
            }
        }
    }
    (k, m, f)
}

/// The local unknowns of the two bending planes: x–y (v, θz) with signs
/// (+, +), x–z (w, θy) with signs (+, −) because θy = −dw/dx.
const PLANE_XY: ([usize; 4], [f64; 4]) = ([1, 5, 7, 11], [1.0, 1.0, 1.0, 1.0]);
const PLANE_XZ: ([usize; 4], [f64; 4]) = ([2, 4, 8, 10], [1.0, -1.0, 1.0, -1.0]);

/// Local stiffness and consistent mass (12 × 12, row major) of a beam of
/// length `len`.
pub fn local_matrices(p: &BeamProps, len: f64) -> ([f64; 144], [f64; 144]) {
    let mut k = [0.0; 144];
    let mut m = [0.0; 144];
    let timo = p.theory == Theory::Timoshenko;
    let rho = p.density;
    // axial and torsion
    let ea = p.young * p.area / len;
    let gj = p.shear * p.torsion / len;
    let ma = rho * p.area * len / 6.0;
    let mt = rho * (p.iy + p.iz) * len / 6.0;
    for (i, j, kk, mm) in [(0, 6, ea, ma), (3, 9, gj, mt)] {
        k[i * 12 + i] += kk;
        k[j * 12 + j] += kk;
        k[i * 12 + j] -= kk;
        k[j * 12 + i] -= kk;
        m[i * 12 + i] += 2.0 * mm;
        m[j * 12 + j] += 2.0 * mm;
        m[i * 12 + j] += mm;
        m[j * 12 + i] += mm;
    }
    // bending: x–y uses I_z and shear along y; x–z uses I_y and shear along z
    for ((idx, sg), inertia, kappa) in [(PLANE_XY, p.iz, p.kappa_y), (PLANE_XZ, p.iy, p.kappa_z)] {
        let kga = timo.then_some(kappa * p.shear * p.area);
        let rho_i = if timo { rho * inertia } else { 0.0 };
        let (kb, mb, _) = plane(len, p.young * inertia, kga, rho * p.area, rho_i);
        for a in 0..4 {
            for b in 0..4 {
                let s = sg[a] * sg[b];
                k[idx[a] * 12 + idx[b]] += s * kb[a * 4 + b];
                m[idx[a] * 12 + idx[b]] += s * mb[a * 4 + b];
            }
        }
    }
    (k, m)
}

/// Work-equivalent local loads (12) of a uniform load per length `q`
/// given in local axes.
pub fn local_line_load(p: &BeamProps, len: f64, q: [f64; 3]) -> [f64; 12] {
    let mut f = [0.0; 12];
    f[0] = 0.5 * q[0] * len;
    f[6] = 0.5 * q[0] * len;
    let timo = p.theory == Theory::Timoshenko;
    for ((idx, sg), inertia, kappa, qq) in [
        (PLANE_XY, p.iz, p.kappa_y, q[1]),
        (PLANE_XZ, p.iy, p.kappa_z, q[2]),
    ] {
        let kga = timo.then_some(kappa * p.shear * p.area);
        let (_, _, fb) = plane(len, p.young * inertia, kga, 0.0, 0.0);
        for a in 0..4 {
            f[idx[a]] += sg[a] * qq * fb[a];
        }
    }
    f
}

/// The rotation from global to local axes (rows: local x, y, z) of a beam
/// from `a` to `b`; `y_hint` is a direction in the local x–y plane
/// (the local y axis after orthogonalisation). Returns the length too.
pub fn frame(a: [f64; 3], b: [f64; 3], y_hint: [f64; 3]) -> Option<([[f64; 3]; 3], f64)> {
    let d = sub(b, a);
    let len = norm(d);
    if len.is_nan() || len <= 0.0 {
        return None;
    }
    let ex = scale(d, 1.0 / len);
    let y = sub(y_hint, scale(ex, dot3(y_hint, ex)));
    let ny = norm(y);
    if ny.is_nan() || ny <= 1e-9 * norm(y_hint).max(f64::MIN_POSITIVE) {
        return None;
    }
    let ey = scale(y, 1.0 / ny);
    let ez = cross(ex, ey);
    Some(([ex, ey, ez], len))
}

/// The default local y axis of a beam along `axis`: the section's
/// "thickness" direction z is the global z where possible (a reed lying
/// in the x–y plane bends about local y), otherwise global x.
pub fn default_y(axis: [f64; 3]) -> [f64; 3] {
    let n = norm(axis);
    let ex = scale(axis, 1.0 / n);
    let z = if dot3(ex, [0.0, 0.0, 1.0]).abs() < 0.99 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    cross(z, ex)
}

/// Global = Tᵀ local T for a 12 × 12 matrix with the rotation `r`.
pub fn to_global(r: &[[f64; 3]; 3], kl: &[f64; 144]) -> [f64; 144] {
    // T is block diagonal with four copies of r
    let mut t = [0.0; 144];
    for blk in 0..4 {
        for i in 0..3 {
            for j in 0..3 {
                t[(3 * blk + i) * 12 + 3 * blk + j] = r[i][j];
            }
        }
    }
    let mut kt = [0.0; 144];
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for l in 0..12 {
                s += kl[i * 12 + l] * t[l * 12 + j];
            }
            kt[i * 12 + j] = s;
        }
    }
    let mut g = [0.0; 144];
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for l in 0..12 {
                s += t[l * 12 + i] * kt[l * 12 + j];
            }
            g[i * 12 + j] = s;
        }
    }
    g
}

/// Tᵀ f for a local vector of 12.
pub fn vec_to_global(r: &[[f64; 3]; 3], fl: &[f64; 12]) -> [f64; 12] {
    let mut g = [0.0; 12];
    for blk in 0..4 {
        for j in 0..3 {
            g[3 * blk + j] = (0..3).map(|i| r[i][j] * fl[3 * blk + i]).sum();
        }
    }
    g
}

/// R v.
pub fn rotate(r: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [dot3(r[0], v), dot3(r[1], v), dot3(r[2], v)]
}

pub(crate) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(crate) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(crate) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(crate) fn norm(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn props(theory: Theory) -> BeamProps {
        BeamProps {
            young: 200e9,
            shear: 80e9,
            density: 7800.0,
            area: 1e-4,
            iy: 2e-9,
            iz: 3e-9,
            torsion: 1e-9,
            kappa_y: 5.0 / 6.0,
            kappa_z: 5.0 / 6.0,
            theory,
        }
    }

    #[test]
    fn closed_form_stiffness_terms() {
        // the classical entries: 12EI/(L³(1+Φ)), (4+Φ)EI/(L(1+Φ)),
        // (2−Φ)EI/(L(1+Φ)) (Przemieniecki 1968, ch. 5)
        let len = 0.05;
        for th in [Theory::EulerBernoulli, Theory::Timoshenko] {
            let p = props(th);
            let (k, m) = local_matrices(&p, len);
            let ei = p.young * p.iy;
            let phi = if th == Theory::Timoshenko {
                12.0 * ei / (p.kappa_z * p.shear * p.area * len * len)
            } else {
                0.0
            };
            let c = ei / (1.0 + phi);
            let rel = |a: f64, b: f64| ((a - b) / b).abs();
            assert!(rel(k[2 * 12 + 2], 12.0 * c / len.powi(3)) < 1e-12);
            assert!(rel(k[4 * 12 + 4], (4.0 + phi) * c / len) < 1e-12);
            assert!(rel(k[4 * 12 + 10], (2.0 - phi) * c / len) < 1e-12);
            assert!(rel(k[2 * 12 + 4], -6.0 * c / len.powi(2)) < 1e-12);
            // symmetric
            for i in 0..12 {
                for j in 0..12 {
                    assert!((k[i * 12 + j] - k[j * 12 + i]).abs() <= 1e-9 * k[i * 12 + i].abs());
                    assert!((m[i * 12 + j] - m[j * 12 + i]).abs() <= 1e-12 * m[i * 12 + i].abs());
                }
            }
            // translational mass sums to ρAL in each direction
            let mass: f64 = (0..12)
                .map(|j| m[2 * 12 + j] * [0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0.][j])
                .sum::<f64>()
                + (0..12)
                    .map(|j| m[8 * 12 + j] * [0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0.][j])
                    .sum::<f64>();
            assert!(rel(mass, p.density * p.area * len) < 1e-12);
        }
    }

    #[test]
    fn rigid_motions_have_no_strain() {
        let p = props(Theory::Timoshenko);
        let a = [0.1, 0.2, -0.3];
        let b = [0.4, -0.1, 0.5];
        let (r, len) = frame(a, b, default_y(sub(b, a))).unwrap();
        let (kl, _) = local_matrices(&p, len);
        let kg = to_global(&r, &kl);
        // a small rigid rotation ω about a point plus a translation
        let om = [0.3, -0.2, 0.7];
        let t = [1.0, 2.0, 3.0];
        let mut u = [0.0; 12];
        for (n, x) in [a, b].iter().enumerate() {
            let rot = cross(om, *x);
            for d in 0..3 {
                u[6 * n + d] = t[d] + rot[d];
                u[6 * n + 3 + d] = om[d];
            }
        }
        // global dof order per node is (u, v, w, θx, θy, θz): matches the
        // block T, so K u must vanish
        for i in 0..12 {
            let f: f64 = (0..12).map(|j| kg[i * 12 + j] * u[j]).sum();
            assert!(f.abs() < 1e-6 * kg[i * 12 + i].abs(), "{i}: {f}");
        }
    }
}
