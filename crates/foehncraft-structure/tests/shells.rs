//! Verification of the MITC4 shell (MS5 part 2; plan section 6.4, rows
//! "cantilever plate" and "shell patch tests and locking tests").
//!
//! References: MacNeal & Harder (1985), *Finite Elements in Analysis and
//! Design* 1, the standard set (patch tests on the distorted five-element
//! mesh, Scordelis–Lo roof); Belytschko et al. (1985), *Comput. Methods
//! Appl. Mech. Engrg.* 51 (the roof's 0.3024); Timoshenko &
//! Woinowsky-Krieger, *Theory of Plates and Shells* (1959), the Navier
//! series of the simply supported plate; Leissa (1969), *Vibration of
//! Plates*, NASA SP-160, the cantilever square plate.

#![allow(clippy::needless_range_loop)]

use foehncraft_structure::analysis;
use foehncraft_structure::model::{AreaLoad, DOFS, Model, Shell};
use foehncraft_structure::shell::ShellProps;
use std::f64::consts::PI;

fn props(e: f64, nu: f64, rho: f64, t: f64) -> ShellProps {
    ShellProps {
        young: e,
        poisson: nu,
        density: rho,
        thickness: t,
        kappa: 5.0 / 6.0,
    }
}

fn rel(a: f64, b: f64) -> f64 {
    ((a - b) / b).abs()
}

/// The MacNeal–Harder patch: a 0.24 × 0.12 rectangle cut into five
/// distorted quadrilaterals by four interior nodes.
fn patch(t: f64, tri: bool) -> (Model, Vec<usize>) {
    let mut m = Model::default();
    let pts = [
        [0.0, 0.0],
        [0.24, 0.0],
        [0.24, 0.12],
        [0.0, 0.12],
        [0.04, 0.02],
        [0.18, 0.03],
        [0.16, 0.08],
        [0.08, 0.08],
    ];
    for p in pts {
        m.add_node([p[0], p[1], 0.0]);
    }
    m.shell_props.push(props(1.0e6, 0.25, 1.0, t));
    for q in [[0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7], [4, 5, 6, 7]] {
        if tri {
            for t in [[q[0], q[1], q[2]], [q[0], q[2], q[3]]] {
                m.shells.push(Shell { nodes: [t[0], t[1], t[2], usize::MAX], tri: true, props: 0 });
            }
        } else {
            m.shells.push(Shell { nodes: q, tri: false, props: 0 });
        }
    }
    (m, vec![0, 1, 2, 3])
}

#[test]
fn membrane_patch_test() {
    for tri in [false, true] {
        membrane_patch(tri);
    }
}

fn membrane_patch(tri: bool) {
    // u = 1e-3 (x + y/2), v = 1e-3 (y + x/2) on the boundary: constant
    // strain ε_x = ε_y = γ_xy = 1e-3 everywhere
    let (mut m, boundary) = patch(0.001, tri);
    let field = |x: [f64; 3]| [1e-3 * (x[0] + 0.5 * x[1]), 1e-3 * (x[1] + 0.5 * x[0])];
    for n in 0..m.nodes.len() {
        // out of plane held everywhere: the test is of the membrane
        for d in [2, 3, 4, 5] {
            m.fixed[DOFS * n + d] = true;
        }
    }
    for &n in &boundary {
        let v = field(m.nodes[n]);
        for d in 0..2 {
            m.fixed[DOFS * n + d] = true;
            m.prescribed.push((DOFS * n + d, v[d]));
        }
    }
    let s = analysis::solve_static(&m).unwrap();
    for n in 4..8 {
        let v = field(m.nodes[n]);
        for d in 0..2 {
            assert!((s.u[DOFS * n + d] - v[d]).abs() < 1e-12, "tri {tri} node {n} {d}: {} {}", s.u[DOFS * n + d], v[d]);
        }
    }
}

#[test]
fn bending_patch_test() {
    for tri in [false, true] {
        bending_patch(tri);
    }
}

fn bending_patch(tri: bool) {
    // w = 1e-3 (x² + xy + y²)/2, θx = ∂w/∂y, θy = −∂w/∂x on the boundary:
    // constant curvatures, the interior nodes follow exactly
    let (mut m, boundary) = patch(0.001, tri);
    let w = |x: [f64; 3]| 1e-3 * (x[0] * x[0] + x[0] * x[1] + x[1] * x[1]) / 2.0;
    let rx = |x: [f64; 3]| 1e-3 * (x[0] + 2.0 * x[1]) / 2.0;
    let ry = |x: [f64; 3]| -1e-3 * (2.0 * x[0] + x[1]) / 2.0;
    for n in 0..m.nodes.len() {
        for d in [0, 1, 5] {
            m.fixed[DOFS * n + d] = true;
        }
    }
    for &n in &boundary {
        let x = m.nodes[n];
        for (d, v) in [(2, w(x)), (3, rx(x)), (4, ry(x))] {
            m.fixed[DOFS * n + d] = true;
            m.prescribed.push((DOFS * n + d, v));
        }
    }
    let s = analysis::solve_static(&m).unwrap();
    for n in 4..8 {
        let x = m.nodes[n];
        for (d, v) in [(2, w(x)), (3, rx(x)), (4, ry(x))] {
            assert!((s.u[DOFS * n + d] - v).abs() < 1e-9 * 1e-4, "tri {tri} node {n} {d}: {} {v}", s.u[DOFS * n + d]);
        }
    }
}

/// w at the centre of a simply supported square plate under a uniform
/// load: the Navier series (Timoshenko & Woinowsky-Krieger, §30), with the
/// shear deflection of the Reissner–Mindlin plate, which is exact for the
/// hard simple support (Wang, Reddy & Lee, *Shear Deformable Beams and
/// Plates*, 2000, ch. 7): each term grows by 1 + D k²/(κGh),
/// k² = π²(m² + n²)/a².
fn navier_centre(q: f64, a: f64, d: f64, kgh: f64) -> f64 {
    let mut s = 0.0;
    for m in (1..400).step_by(2) {
        for n in (1..400).step_by(2) {
            let (mf, nf) = (m as f64, n as f64);
            let sign = ((-1i64).pow(((m + n) / 2 - 1) as u32)) as f64;
            let k2 = PI * PI * (mf * mf + nf * nf) / (a * a);
            s += sign / (mf * nf * (mf * mf + nf * nf).powi(2)) * (1.0 + d * k2 / kgh);
        }
    }
    16.0 * q * a.powi(4) / (PI.powi(6) * d) * s
}

#[test]
fn simply_supported_plate_converges_to_navier() {
    for tri in [false, true] {
        navier(tri);
    }
}

fn grid(m: &mut Model, tri: bool, at: impl Fn(f64, f64) -> [f64; 3], n: usize) -> Vec<Vec<usize>> {
    if tri { m.add_tri_grid(at, n, n, 0) } else { m.add_shell_grid(at, n, n, 0) }
}

fn navier(tri: bool) {
    for t in [0.01f64, 0.0001] {
        navier_t(tri, t);
    }
}

fn navier_t(tri: bool, t: f64) {
    let (a, e, nu, q): (f64, f64, f64, f64) = (1.0, 2.0e11, 0.3, 1.0e3);
    let d = e * t.powi(3) / (12.0 * (1.0 - nu * nu));
    let kgh = 5.0 / 6.0 * e / (2.0 * (1.0 + nu)) * t;
    let exact = navier_centre(q, a, d, kgh);
    let mut errs = Vec::new();
    for n in [8, 16, 32] {
        let mut m = Model::default();
        m.shell_props.push(props(e, nu, 7800.0, t));
        let g = grid(&mut m, tri, |u, v| [a * u, a * v, 0.0], n);
        for i in 0..=n {
            for j in 0..=n {
                let k = g[i][j];
                // membrane and drilling held: a plate test
                for dd in [0, 1, 5] {
                    m.fixed[DOFS * k + dd] = true;
                }
                // hard simple support: w = 0 and the rotation along the edge
                if i == 0 || i == n {
                    m.fixed[DOFS * k + 2] = true;
                    m.fixed[DOFS * k + 3] = true;
                }
                if j == 0 || j == n {
                    m.fixed[DOFS * k + 2] = true;
                    m.fixed[DOFS * k + 4] = true;
                }
            }
        }
        for el in 0..m.shells.len() {
            m.area_loads.push(AreaLoad { element: el, q: [0.0; 3], pressure: -q });
        }
        let s = analysis::solve_static(&m).unwrap();
        let w = s.u[DOFS * g[n / 2][n / 2] + 2];
        errs.push(rel(-w, exact));
    }
    eprintln!("Navier, tri {tri}, t {t}: {errs:?}");
    // MITC3+ on the same grid (twice the elements) is about nine times
    // less accurate here than MITC4, but converges at the same order and
    // does not lock as the plate thins
    let band = if tri { 4e-3 } else { 4e-4 };
    assert!(errs[2] < band, "{errs:?}");
    let order = (errs[1] / errs[2]).log2();
    assert!((1.8..2.2).contains(&order), "order {order}: {errs:?}");
}

#[test]
fn cantilever_square_plate_modes_against_leissa() {
    for tri in [false, true] {
        leissa(tri);
    }
}

fn leissa(tri: bool) {
    // λ = ω a² √(ρh/D) of the square plate clamped on one edge, ν = 0.3
    // (Leissa 1969, table 4.30, after Young 1950 and Barton 1951)
    const LEISSA: [f64; 5] = [3.492, 8.525, 21.43, 27.33, 31.11];
    let (a, t, e, nu, rho): (f64, f64, f64, f64, f64) = (0.1, 0.5e-3, 2.0e11, 0.3, 7800.0);
    let d = e * t.powi(3) / (12.0 * (1.0 - nu * nu));
    let n = 24;
    let mut m = Model::default();
    m.shell_props.push(props(e, nu, rho, t));
    let g = grid(&mut m, tri, |u, v| [a * u, a * v, 0.0], n);
    for i in 0..=n {
        for j in 0..=n {
            let k = g[i][j];
            for dd in [0, 1, 5] {
                m.fixed[DOFS * k + dd] = true;
            }
            if i == 0 {
                for dd in 0..DOFS {
                    m.fixed[DOFS * k + dd] = true;
                }
            }
        }
    }
    let md = analysis::solve_modal(&m, 5, 0.0).unwrap();
    assert!(md.sturm.ok(), "{:?}", md.sturm);
    for (i, mo) in md.modes.iter().enumerate() {
        let lam = mo.lambda.sqrt() * a * a * (rho * t / d).sqrt();
        eprintln!("tri {tri} mode {}: λ = {lam:.4} (Leissa {})", i + 1, LEISSA[i]);
        assert!(rel(lam, LEISSA[i]) < 1e-2, "mode {}: {lam} {}", i + 1, LEISSA[i]);
    }
}

/// The Scordelis–Lo roof (MacNeal & Harder 1985): a cylindrical panel of
/// radius 25, length 50, 80° opening, thickness 0.25, E = 4.32e8, ν = 0,
/// its own weight 90 per area, on rigid diaphragms at the ends, the long
/// edges free. One quarter by symmetry; the vertical displacement at the
/// middle of the free edge is 0.3024.
fn scordelis_lo(n: usize, tri: bool) -> f64 {
    let (r, len, t) = (25.0, 50.0, 0.25);
    let phi = 40.0f64.to_radians();
    let mut m = Model::default();
    m.shell_props.push(props(4.32e8, 0.0, 1.0, t));
    // a along x from the diaphragm (0) to the middle (L/2), b along the arc
    // from the crown (0) to the free edge (40°)
    let g = grid(
        &mut m,
        tri,
        |u, v| {
            let a = v * phi;
            [u * len / 2.0, r * a.sin(), r * a.cos()]
        },
        n,
    );
    for i in 0..=n {
        for j in 0..=n {
            let k = g[i][j];
            let mut fix = |d: usize| m.fixed[DOFS * k + d] = true;
            if i == 0 {
                fix(1); // diaphragm: no motion in its plane
                fix(2);
            }
            if i == n {
                fix(0); // symmetry about x = L/2
                fix(4);
                fix(5);
            }
            if j == 0 {
                fix(1); // symmetry about y = 0
                fix(3);
                fix(5);
            }
        }
    }
    for el in 0..m.shells.len() {
        m.area_loads.push(AreaLoad { element: el, q: [0.0, 0.0, -90.0], pressure: 0.0 });
    }
    let s = analysis::solve_static(&m).unwrap();
    -s.u[DOFS * g[n][n] + 2]
}

#[test]
fn scordelis_lo_roof_converges() {
    // MITC4: within 2 % on 16 × 16. MITC3+ converges more slowly from
    // below; Ko, Lee, Lee & Bathe (2017, Comput. Struct. 193, table 10)
    // report w/w_ref = 0.6695, 0.8577, 0.9550, 0.9851 on the regular
    // N × N meshes N = 4, 8, 16, 32 (their mesh pattern differs from our
    // alternating diagonals, hence 3 %).
    const KO_2017: [f64; 4] = [0.6695, 0.8577, 0.9550, 0.9851];
    for tri in [false, true] {
        let w: Vec<f64> = [4, 8, 16, 32].iter().map(|&n| scordelis_lo(n, tri)).collect();
        eprintln!("Scordelis–Lo, tri {tri}: {w:?}");
        if tri {
            for k in 2..4 {
                let ratio = w[k] / 0.3024;
                assert!((ratio - KO_2017[k]).abs() < 3e-2, "N index {k}: {ratio} {w:?}");
            }
            assert!(w.windows(2).all(|p| p[1] > p[0]), "not monotone: {w:?}");
        } else {
            assert!(rel(w[2], 0.3024) < 2e-2, "{w:?}");
        }
    }
}
