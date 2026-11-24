//! L17 A0 option B spike: compiler reverse mode (Enzyme, rustc nightly
//! `#[autodiff_reverse]`) on the real momentum face kernel, interior face,
//! compared with option A (local Dual<26> Jacobian, then Jᵀ ȳ).
#![feature(autodiff)]
use foehncraft_kernels::flow::{self as k, FaceGeometry, MomentumInputs, VelocityGradient};
use foehncraft_kernels::fv::{Topology, Vec3};
use foehncraft_kernels::reverse as rv;
use std::autodiff::autodiff_reverse;
use std::hint::black_box;

pub struct Geo {
    pub area: [f64; 3],
    pub fc: [f64; 3],
    pub k: [f64; 3],
    pub dgeo: f64,
    pub w: f64,
    pub cc: [[f64; 2]; 3],
    pub psi: [[f64; 2]; 3],
}

fn v<'a>(a: &'a [f64], b: &'a [f64], c: &'a [f64]) -> Vec3<'a, f64> {
    Vec3 { x: a, y: b, z: c }
}

/// The local face problem on f64: x in the layout of `reverse` (26 inputs).
#[inline(never)]
#[autodiff_reverse(d_face, Const, Duplicated, Duplicated)]
pub fn face(g: &Geo, x: &[f64; 26], y: &mut [f64; 5]) {
    let m = [x[0]];
    let mu = [x[1]];
    let u = [[x[2], x[5]], [x[3], x[6]], [x[4], x[7]]];
    let mut gr = [[[0.0; 2]; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            gr[i][j] = [x[8 + 3 * i + j], x[17 + 3 * i + j]];
        }
    }
    let owner = [0u32];
    let nb = [1u32];
    let (a0, a1, a2) = ([g.area[0]], [g.area[1]], [g.area[2]]);
    let (f0, f1, f2) = ([g.fc[0]], [g.fc[1]], [g.fc[2]]);
    let (k0, k1, k2) = ([g.k[0]], [g.k[1]], [g.k[2]]);
    let dgeo = [g.dgeo];
    let w = [g.w];
    let e: &[f64] = &[];
    let inp = MomentumInputs {
        topo: Topology { owner: &owner, neighbour: &nb, cell_ptr: &[], cell_faces: &[], n_interior: 1, periodic: 1 },
        geo: FaceGeometry {
            area: v(&a0, &a1, &a2),
            face_centroid: v(&f0, &f1, &f2),
            cell_centroid: v(&g.cc[0], &g.cc[1], &g.cc[2]),
            dgeo: &dgeo,
            k: v(&k0, &k1, &k2),
            weight: &w,
            shift: v(e, e, e),
        },
        mass_flux: &m,
        mu: &mu,
        u: v(&u[0], &u[1], &u[2]),
        uc: v(&u[0], &u[1], &u[2]),
        grad: VelocityGradient {
            x: v(&gr[0][0], &gr[0][1], &gr[0][2]),
            y: v(&gr[1][0], &gr[1][1], &gr[1][2]),
            z: v(&gr[2][0], &gr[2][1], &gr[2][2]),
        },
        psi: v(&g.psi[0], &g.psi[1], &g.psi[2]),
        ub: v(e, e, e),
        mom_kind: &[],
        blend: 1.0,
        central: 0.0,
        beta: 0.0,
        nonorth: 1.0,
        convection: 1.0,
        transposed: 0.0,
        trace_free: 0.0,
    };
    // SAFETY: one interior face, two cells, every array long enough.
    let (cp, cn, b) = unsafe { k::momentum_face(0, &inp) };
    *y = [cp, cn, b[0], b[1], b[2]];
}

fn main() {
    let g = Geo {
        area: [0.31, 0.02, -0.05],
        fc: [0.5, 0.21, 0.13],
        k: [0.01, 0.02, -0.03],
        dgeo: 3.3,
        w: 0.46,
        cc: [[0.33, 0.68], [0.19, 0.24], [0.11, 0.16]],
        psi: [[0.8, 0.6], [1.0, 0.9], [0.7, 0.5]],
    };
    let mut x = [0.0; 26];
    for (i, a) in x.iter_mut().enumerate() {
        *a = 0.3 + 0.07 * (i as f64 * 1.37).sin();
    }
    let ybar = [0.4, -0.7, 1.1, 0.25, -0.6];
    // option B
    let mut y = [0.0; 5];
    let mut dx = [0.0; 26];
    let mut dy = ybar;
    d_face(&g, &x, &mut dx, &mut y, &mut dy);
    // option A on the same face (copy into a two-cell mesh)
    let owner = [0u32];
    let nb = [1u32];
    let cols = |o: usize| -> [f64; 2] { [x[o], x[o + 3]] };
    let ux = cols(2); let uy = cols(3); let uz = cols(4);
    let gcol = |i: usize, j: usize| -> [f64; 2] { [x[8 + 3 * i + j], x[17 + 3 * i + j]] };
    let gr: Vec<[f64; 2]> = (0..9).map(|n| gcol(n / 3, n % 3)).collect();
    let (a0, a1, a2) = ([g.area[0]], [g.area[1]], [g.area[2]]);
    let (f0, f1, f2) = ([g.fc[0]], [g.fc[1]], [g.fc[2]]);
    let (k0, k1, k2) = ([g.k[0]], [g.k[1]], [g.k[2]]);
    let (m, mu, dgeo, w) = ([x[0]], [x[1]], [g.dgeo], [g.w]);
    let e: &[f64] = &[];
    let inp = MomentumInputs {
        topo: Topology { owner: &owner, neighbour: &nb, cell_ptr: &[], cell_faces: &[], n_interior: 1, periodic: 1 },
        geo: FaceGeometry { area: v(&a0, &a1, &a2), face_centroid: v(&f0, &f1, &f2), cell_centroid: v(&g.cc[0], &g.cc[1], &g.cc[2]), dgeo: &dgeo, k: v(&k0, &k1, &k2), weight: &w, shift: v(e, e, e) },
        mass_flux: &m, mu: &mu, u: v(&ux, &uy, &uz), uc: v(&ux, &uy, &uz),
        grad: VelocityGradient { x: v(&gr[0], &gr[1], &gr[2]), y: v(&gr[3], &gr[4], &gr[5]), z: v(&gr[6], &gr[7], &gr[8]) },
        psi: v(&g.psi[0], &g.psi[1], &g.psi[2]), ub: v(e, e, e), mom_kind: &[],
        blend: 1.0, central: 0.0, beta: 0.0, nonorth: 1.0, convection: 1.0, transposed: 0.0, trace_free: 0.0,
    };
    let (jac, val) = rv::momentum_face_jacobian(0, &inp).unwrap();
    let xa = rv::pullback(&jac, &ybar);
    let mut worst: f64 = 0.0;
    for i in 0..26 {
        let s = xa[i].abs().max(1e-300);
        worst = worst.max((xa[i] - dx[i]).abs() / s.max(1e-12));
    }
    println!("values: enzyme {y:?}\n        dual   {val:?}");
    println!("largest relative difference Jᵀȳ (option B vs option A): {worst:.2e}");
    // timing
    let n = 2_000_000;
    let t = std::time::Instant::now();
    let mut s = 0.0;
    for _ in 0..n { let mut yy = [0.0; 5]; face(black_box(&g), black_box(&x), &mut yy); s += yy[2]; }
    let tp = t.elapsed().as_secs_f64() / n as f64;
    let t = std::time::Instant::now();
    for _ in 0..n { let mut yy = [0.0; 5]; let mut d = [0.0; 26]; let mut dyy = ybar; d_face(black_box(&g), black_box(&x), &mut d, &mut yy, &mut dyy); s += d[3]; }
    let tb = t.elapsed().as_secs_f64() / n as f64;
    let n2 = n / 10;
    let t = std::time::Instant::now();
    for _ in 0..n2 { let (j, _) = rv::momentum_face_jacobian(0, black_box(&inp)).unwrap(); s += rv::pullback(&j, &ybar)[3]; }
    let ta = t.elapsed().as_secs_f64() / n2 as f64;
    black_box(s);
    println!("per face: primal (local copy) {:.1} ns, option B reverse {:.1} ns ({:.1} x), option A recompute {:.1} ns ({:.1} x)", tp*1e9, tb*1e9, tb/tp, ta*1e9, ta/tp);
}
