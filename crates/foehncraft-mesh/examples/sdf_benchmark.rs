//! Signed distance on a grid by three methods (L13 stage 0, D225 item 4;
//! owner rule: benchmark the computer-graphics methods; decision Q11: the
//! default is decided by this benchmark).
//!
//! ```sh
//! cargo run --release -p foehncraft-mesh --example sdf_benchmark -- output/sdf 48 96
//! python3 -m tools.bench.sdf output/sdf
//! ```
//!
//! Methods, each on the nodes of an n³ grid over [-1.5, 1.5]³:
//!
//! * `ours`: exact closest point on the surface hierarchy, sign from the
//!   generalised winding number (`foehncraft_mesh::surface`).
//! * `sweep`: fast sweeping for the eikonal equation (Zhao 2005, Math.
//!   Comp. 74; Godunov upwind update, 8 orderings until no node changes),
//!   started from the exact distance at the nodes within 1.5 cells of the
//!   surface; the sign by scan conversion: parity of the crossings along
//!   each grid line in x (the classical choice, as in Roosing et al. 2019;
//!   it needs a closed, correctly oriented surface).
//! * `heat`: the signed heat method (Feng & Crane 2024, ACM TOG 43(4)):
//!   the area-weighted normals of the surface, spread to the nodes, are
//!   diffused for one step `(I − tΔ) Y = N δ_Γ` with `t = h²`, normalised to
//!   `X = Y/|Y|`, and `φ` is the least-squares integral of `X`
//!   (`Δφ = ∇·X`), shifted to zero mean on the surface. Written from the
//!   paper's description on a grid (7-point Laplacian, zero-flux sides,
//!   conjugate gradients), no code read.
//!
//! Written per surface and grid: `<shape>_<n>_<method>.f64` (raw
//! little-endian binary64, x fastest) and `bench.json` with the times.

use foehncraft_mesh::format::stl::{Triangle, to_binary};
use foehncraft_mesh::surface::Surface;
use rayon::prelude::*;
use std::time::Instant;

const LO: f64 = -1.5;
const HI: f64 = 1.5;

fn icosphere(c: [f64; 3], r: f64, levels: usize) -> Vec<Triangle> {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut v: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut f: Vec<[usize; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    let unit = |p: [f64; 3]| {
        let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        p.map(|x| x / n)
    };
    for p in v.iter_mut() {
        *p = unit(*p);
    }
    for _ in 0..levels {
        let mut mid = std::collections::HashMap::new();
        let mut nf = Vec::with_capacity(4 * f.len());
        for tri in &f {
            let mut m = [0usize; 3];
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                m[k] = *mid.entry((a.min(b), a.max(b))).or_insert_with(|| {
                    v.push(unit([0, 1, 2].map(|d| 0.5 * (v[a][d] + v[b][d]))));
                    v.len() - 1
                });
            }
            nf.push([tri[0], m[0], m[2]]);
            nf.push([tri[1], m[1], m[0]]);
            nf.push([tri[2], m[2], m[1]]);
            nf.push([m[0], m[1], m[2]]);
        }
        f = nf;
    }
    f.iter()
        .map(|tri| tri.map(|i| [0, 1, 2].map(|d| c[d] + r * v[i][d])))
        .collect()
}

fn torus(big: f64, small: f64, n: usize, m: usize) -> Vec<Triangle> {
    let p = |i: usize, j: usize| {
        let u = std::f64::consts::TAU * (i % n) as f64 / n as f64;
        let v = std::f64::consts::TAU * (j % m) as f64 / m as f64;
        [
            (big + small * v.cos()) * u.cos(),
            (big + small * v.cos()) * u.sin(),
            small * v.sin(),
        ]
    };
    let mut out = Vec::new();
    for i in 0..n {
        for j in 0..m {
            let (a, b, c, d) = (p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1));
            out.push([a, b, c]);
            out.push([a, c, d]);
        }
    }
    out
}

/// Box with half sides `h`, rotated by `ax` about x, then `az` about z
/// (the same rotation as `tools/bench/isosurface.py`).
fn rotated_box(h: [f64; 3], ax: f64, az: f64) -> Vec<Triangle> {
    let (cx, sx, cz, sz) = (ax.cos(), ax.sin(), az.cos(), az.sin());
    let rot = |p: [f64; 3]| {
        let q = [p[0], cx * p[1] - sx * p[2], sx * p[1] + cx * p[2]];
        [cz * q[0] - sz * q[1], sz * q[0] + cz * q[1], q[2]]
    };
    let v = |i: usize| {
        rot([
            if i & 1 == 0 { -h[0] } else { h[0] },
            if i >> 1 & 1 == 0 { -h[1] } else { h[1] },
            if i >> 2 & 1 == 0 { -h[2] } else { h[2] },
        ])
    };
    let quads = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let mut out = Vec::new();
    for q in quads {
        out.push([v(q[0]), v(q[1]), v(q[2])]);
        out.push([v(q[0]), v(q[2]), v(q[3])]);
    }
    out
}

fn centroid(t: &Triangle) -> [f64; 3] {
    [0, 1, 2].map(|d| (t[0][d] + t[1][d] + t[2][d]) / 3.0)
}

/// The surfaces: clean ones with an exact distance, and dirty ones whose
/// reference is the clean sphere.
fn shapes() -> Vec<(&'static str, Vec<Triangle>)> {
    let sphere = icosphere([0.0; 3], 1.0, 5);
    let holed: Vec<Triangle> = sphere
        .iter()
        .copied()
        .filter(|t| {
            let c = centroid(t);
            let n = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
            (0.3 * c[0] + 0.954 * c[2]) / n < 0.97
        })
        .collect();
    let mut flipped = sphere.clone();
    for t in flipped.iter_mut() {
        let c = centroid(t);
        if c[0] < -0.5 && c[2] < -0.5 {
            *t = [t[0], t[2], t[1]];
        }
    }
    // duplicated half and an inner shell overlapping the outer surface
    let mut soup = sphere.clone();
    soup.extend(sphere.iter().filter(|t| centroid(t)[1] > 0.2).copied());
    soup.extend(icosphere([0.55, 0.1, 0.0], 0.5, 3));
    vec![
        ("sphere", sphere),
        ("torus", torus(0.9, 0.35, 128, 48)),
        ("box", rotated_box([0.55, 0.38, 0.27], 0.4, 0.3)),
        ("holed", holed),
        ("flipped", flipped),
        ("soup", soup),
    ]
}

/// Node coordinate of index `i` on an n-cell axis.
fn coord(i: usize, n: usize) -> f64 {
    LO + (HI - LO) * i as f64 / n as f64
}

fn nodes(n: usize) -> Vec<[f64; 3]> {
    let m = n + 1;
    (0..m * m * m)
        .map(|k| [coord(k % m, n), coord(k / m % m, n), coord(k / (m * m), n)])
        .collect()
}

/// Fast sweeping (Zhao 2005) with scan-conversion sign.
fn sweep(s: &Surface, n: usize) -> Vec<f64> {
    let m = n + 1;
    let h = (HI - LO) / n as f64;
    let pts = nodes(n);
    let band = 1.5 * h;
    let mut u: Vec<f64> = pts
        .par_iter()
        .map(|&x| {
            s.bvh()
                .closest_within(x, band)
                .map_or(f64::INFINITY, |c| c.distance)
        })
        .collect();
    let fixed: Vec<bool> = u.iter().map(|v| v.is_finite()).collect();
    let idx = |i: usize, j: usize, k: usize| i + m * (j + m * k);
    loop {
        let mut changed = false;
        for dir in 0..8 {
            let rev = |a: usize, bit: usize| if dir >> bit & 1 == 1 { m - 1 - a } else { a };
            for kk in 0..m {
                let k = rev(kk, 2);
                for jj in 0..m {
                    let j = rev(jj, 1);
                    for ii in 0..m {
                        let i = rev(ii, 0);
                        let p = idx(i, j, k);
                        if fixed[p] {
                            continue;
                        }
                        let nb = |a: usize, b: usize, c: usize| u[idx(a, b, c)];
                        let mut a = [
                            f64::min(
                                if i > 0 {
                                    nb(i - 1, j, k)
                                } else {
                                    f64::INFINITY
                                },
                                if i + 1 < m {
                                    nb(i + 1, j, k)
                                } else {
                                    f64::INFINITY
                                },
                            ),
                            f64::min(
                                if j > 0 {
                                    nb(i, j - 1, k)
                                } else {
                                    f64::INFINITY
                                },
                                if j + 1 < m {
                                    nb(i, j + 1, k)
                                } else {
                                    f64::INFINITY
                                },
                            ),
                            f64::min(
                                if k > 0 {
                                    nb(i, j, k - 1)
                                } else {
                                    f64::INFINITY
                                },
                                if k + 1 < m {
                                    nb(i, j, k + 1)
                                } else {
                                    f64::INFINITY
                                },
                            ),
                        ];
                        a.sort_by(f64::total_cmp);
                        if !a[0].is_finite() {
                            continue;
                        }
                        // Godunov: the largest root using 1, 2 or 3 directions
                        let mut x = a[0] + h;
                        if x > a[1] {
                            let s2 = a[0] + a[1];
                            let d = s2 * s2 - 2.0 * (a[0] * a[0] + a[1] * a[1] - h * h);
                            x = 0.5 * (s2 + d.max(0.0).sqrt());
                            if x > a[2] {
                                let s3 = a[0] + a[1] + a[2];
                                let q = a[0] * a[0] + a[1] * a[1] + a[2] * a[2] - h * h;
                                let d = s3 * s3 - 3.0 * q;
                                x = (s3 + d.max(0.0).sqrt()) / 3.0;
                            }
                        }
                        if x < u[p] {
                            u[p] = x;
                            changed = true;
                        }
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
    // scan conversion: parity of crossings along each x line
    let signs: Vec<Vec<bool>> = (0..m * m)
        .into_par_iter()
        .map(|jk| {
            let (j, k) = (jk % m, jk / m);
            let (y, z) = (coord(j, n), coord(k, n));
            let p = [LO - 0.1, y, z];
            let q = [HI + 0.1, y, z];
            let len = q[0] - p[0];
            let cross = s.crossings(p, q);
            let mut inside = vec![false; m];
            let mut c = 0;
            let mut state = false;
            for (i, flag) in inside.iter_mut().enumerate() {
                let t = (coord(i, n) - p[0]) / len;
                while c < cross.len() && cross[c].t <= t {
                    state = !state;
                    c += 1;
                }
                *flag = state;
            }
            inside
        })
        .collect();
    for k in 0..m {
        for j in 0..m {
            for i in 0..m {
                if signs[j + m * k][i] {
                    u[idx(i, j, k)] = -u[idx(i, j, k)];
                }
            }
        }
    }
    u
}

/// Matrix-free conjugate gradients for `A x = b` (A symmetric positive
/// (semi-)definite), Jacobi preconditioned; returns the iterations.
fn cg(
    apply: &(dyn Fn(&[f64], &mut [f64]) + Sync),
    diag: &[f64],
    b: &[f64],
    x: &mut [f64],
    tol: f64,
) -> usize {
    let n = b.len();
    let mut r = vec![0.0; n];
    apply(x, &mut r);
    r.par_iter_mut().zip(b).for_each(|(r, b)| *r = b - *r);
    let mut z: Vec<f64> = r.iter().zip(diag).map(|(r, d)| r / d).collect();
    let mut p = z.clone();
    let mut rz: f64 = r.par_iter().zip(&z).map(|(a, b)| a * b).sum();
    let bn: f64 = b.par_iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
    let mut ap = vec![0.0; n];
    for it in 0..20000 {
        let rn: f64 = r.par_iter().map(|v| v * v).sum::<f64>().sqrt();
        if rn <= tol * bn {
            return it;
        }
        apply(&p, &mut ap);
        let pap: f64 = p.par_iter().zip(&ap).map(|(a, b)| a * b).sum();
        let alpha = rz / pap;
        x.par_iter_mut().zip(&p).for_each(|(x, p)| *x += alpha * p);
        r.par_iter_mut().zip(&ap).for_each(|(r, a)| *r -= alpha * a);
        z.par_iter_mut()
            .zip(&r)
            .zip(diag)
            .for_each(|((z, r), d)| *z = r / d);
        let rz2: f64 = r.par_iter().zip(&z).map(|(a, b)| a * b).sum();
        let beta = rz2 / rz;
        rz = rz2;
        p.par_iter_mut()
            .zip(&z)
            .for_each(|(p, z)| *p = z + beta * *p);
    }
    20000
}

/// Graph Laplacian of the grid (zero-flux sides): `(L x)_i = Σ (x_i − x_j)`.
fn laplacian(m: usize, x: &[f64], y: &mut [f64]) {
    y.par_iter_mut().enumerate().for_each(|(p, y)| {
        let (i, j, k) = (p % m, p / m % m, p / (m * m));
        let mut s = 0.0;
        let c = x[p];
        if i > 0 {
            s += c - x[p - 1];
        }
        if i + 1 < m {
            s += c - x[p + 1];
        }
        if j > 0 {
            s += c - x[p - m];
        }
        if j + 1 < m {
            s += c - x[p + m];
        }
        if k > 0 {
            s += c - x[p - m * m];
        }
        if k + 1 < m {
            s += c - x[p + m * m];
        }
        *y = s;
    });
}

fn degree(m: usize) -> Vec<f64> {
    (0..m * m * m)
        .map(|p| {
            let (i, j, k) = (p % m, p / m % m, p / (m * m));
            [i, j, k]
                .iter()
                .map(|&a| (a > 0) as usize + (a + 1 < m) as usize)
                .sum::<usize>() as f64
        })
        .collect()
}

/// Points on a triangle with spacing below `h` (barycentric grid), each
/// with its share of the area.
fn samples(t: &Triangle, h: f64) -> Vec<([f64; 3], f64)> {
    let e = (0..3)
        .map(|k| {
            let d: [f64; 3] = [0, 1, 2].map(|c| t[(k + 1) % 3][c] - t[k][c]);
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .fold(0.0, f64::max);
    let s = ((e / (0.5 * h)).ceil() as usize).max(1);
    let a = {
        let u: [f64; 3] = [0, 1, 2].map(|c| t[1][c] - t[0][c]);
        let v: [f64; 3] = [0, 1, 2].map(|c| t[2][c] - t[0][c]);
        let c = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
    };
    let mut out = Vec::with_capacity(s * s);
    let w = a / (s * s) as f64;
    for i in 0..s {
        for j in 0..s - i {
            // the two sub-triangles of the cell (i, j): centroids
            for (di, dj, up) in [(1.0 / 3.0, 1.0 / 3.0, true), (2.0 / 3.0, 2.0 / 3.0, false)] {
                if !up && i + j + 1 >= s {
                    continue;
                }
                let b1 = (i as f64 + di) / s as f64;
                let b2 = (j as f64 + dj) / s as f64;
                let b0 = 1.0 - b1 - b2;
                out.push((
                    [0, 1, 2].map(|c| b0 * t[0][c] + b1 * t[1][c] + b2 * t[2][c]),
                    w,
                ));
            }
        }
    }
    out
}

/// Trilinear weights of `x`: (node, weight) for the 8 corners of its cell.
fn trilinear(x: [f64; 3], n: usize) -> [(usize, f64); 8] {
    let m = n + 1;
    let h = (HI - LO) / n as f64;
    let mut base = [0usize; 3];
    let mut f = [0.0; 3];
    for d in 0..3 {
        let g = ((x[d] - LO) / h).clamp(0.0, n as f64 - 1e-9);
        base[d] = g.floor() as usize;
        f[d] = g - base[d] as f64;
    }
    std::array::from_fn(|c| {
        let (a, b, cc) = (c & 1, c >> 1 & 1, c >> 2 & 1);
        let w = (if a == 1 { f[0] } else { 1.0 - f[0] })
            * (if b == 1 { f[1] } else { 1.0 - f[1] })
            * (if cc == 1 { f[2] } else { 1.0 - f[2] });
        ((base[0] + a) + m * ((base[1] + b) + m * (base[2] + cc)), w)
    })
}

/// The signed heat method (Feng & Crane 2024) on the grid.
fn heat(tris: &[Triangle], n: usize) -> (Vec<f64>, usize) {
    let m = n + 1;
    let h = (HI - LO) / n as f64;
    let nn = m * m * m;
    // source: area-weighted normals spread trilinearly, per unit volume
    let mut src = vec![[0.0f64; 3]; nn];
    let mut surface_pts: Vec<([f64; 3], f64)> = Vec::new();
    for t in tris {
        let u: [f64; 3] = [0, 1, 2].map(|c| t[1][c] - t[0][c]);
        let v: [f64; 3] = [0, 1, 2].map(|c| t[2][c] - t[0][c]);
        let c = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let a = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        if a == 0.0 {
            continue;
        }
        let nrm = c.map(|x| x / a);
        for (p, w) in samples(t, h) {
            for (node, tw) in trilinear(p, n) {
                for d in 0..3 {
                    src[node][d] += nrm[d] * w * tw / (h * h * h);
                }
            }
            surface_pts.push((p, w));
        }
    }
    // (I − tΔ) Y = S with t = h²: in grid-Laplacian form (I + L) Y = S
    let diag: Vec<f64> = degree(m).iter().map(|d| 1.0 + d).collect();
    let apply = |x: &[f64], y: &mut [f64]| {
        laplacian(m, x, y);
        y.par_iter_mut().zip(x).for_each(|(y, x)| *y += x);
    };
    let mut iters = 0;
    let mut y = vec![[0.0f64; 3]; nn];
    for d in 0..3 {
        let b: Vec<f64> = src.iter().map(|s| s[d]).collect();
        let mut x = vec![0.0; nn];
        iters += cg(&apply, &diag, &b, &mut x, 1e-10);
        for p in 0..nn {
            y[p][d] = x[p];
        }
    }
    let x: Vec<[f64; 3]> = y
        .iter()
        .map(|v| {
            let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if l > 0.0 { v.map(|c| c / l) } else { [0.0; 3] }
        })
        .collect();
    // least squares: Σ_edges (φ_j − φ_i − h X_e)², X_e the mean of the ends
    let mut b = vec![0.0; nn];
    for p in 0..nn {
        let (i, j, k) = (p % m, p / m % m, p / (m * m));
        for (d, (ok, q)) in [
            (i + 1 < m, p + 1),
            (j + 1 < m, p + m),
            (k + 1 < m, p + m * m),
        ]
        .into_iter()
        .enumerate()
        {
            if ok {
                let g = h * 0.5 * (x[p][d] + x[q][d]);
                b[q] += g;
                b[p] -= g;
            }
        }
    }
    let lap = |x: &[f64], y: &mut [f64]| laplacian(m, x, y);
    let mut phi = vec![0.0; nn];
    iters += cg(&lap, &degree(m), &b, &mut phi, 1e-10);
    // zero mean on the surface
    let (mut num, mut den) = (0.0, 0.0);
    for (p, w) in &surface_pts {
        let v: f64 = trilinear(*p, n).iter().map(|(q, tw)| tw * phi[*q]).sum();
        num += w * v;
        den += w;
    }
    let shift = num / den;
    phi.iter_mut().for_each(|v| *v -= shift);
    (phi, iters)
}

fn write(path: &std::path::Path, v: &[f64]) {
    let bytes: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).collect();
    std::fs::write(path, bytes).expect("write field");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out = std::path::PathBuf::from(args.get(1).cloned().unwrap_or("output/sdf".into()));
    let grids: Vec<usize> = if args.len() > 2 {
        args[2..]
            .iter()
            .map(|a| a.parse().expect("grid size"))
            .collect()
    } else {
        vec![48, 96]
    };
    std::fs::create_dir_all(&out).expect("output directory");
    let mut runs = Vec::new();
    for (name, tris) in shapes() {
        std::fs::write(out.join(format!("{name}.stl")), to_binary(&tris)).unwrap();
        let t0 = Instant::now();
        let s = Surface::new(&tris);
        let build = t0.elapsed().as_secs_f64();
        for &n in &grids {
            let pts = nodes(n);
            let t = Instant::now();
            let ours = s.signed_distances(&pts, None);
            let t_ours = t.elapsed().as_secs_f64();
            let t = Instant::now();
            let band = s.signed_distances(&pts, Some(3.0 * (HI - LO) / n as f64));
            let t_band = t.elapsed().as_secs_f64();
            drop(band);
            let t = Instant::now();
            let sw = sweep(&s, n);
            let t_sweep = t.elapsed().as_secs_f64();
            let t = Instant::now();
            let (ht, iters) = heat(&tris, n);
            let t_heat = t.elapsed().as_secs_f64();
            for (m, f) in [("ours", &ours), ("sweep", &sw), ("heat", &ht)] {
                write(&out.join(format!("{name}_{n}_{m}.f64")), f);
            }
            eprintln!(
                "{name} n={n}: ours {t_ours:.3} s (band {t_band:.3}), sweep {t_sweep:.3} s, heat {t_heat:.3} s ({iters} CG its)"
            );
            runs.push(serde_json::json!({
                "shape": name, "n": n, "nodes": pts.len(), "triangles": tris.len(),
                "build_s": build, "ours_s": t_ours, "ours_band_s": t_band,
                "sweep_s": t_sweep, "heat_s": t_heat, "heat_cg_iterations": iters,
            }));
        }
    }
    let json = serde_json::json!({
        "lower": LO, "upper": HI, "threads": rayon::current_num_threads(), "runs": runs,
    });
    std::fs::write(
        out.join("bench.json"),
        serde_json::to_string_pretty(&json).unwrap(),
    )
    .unwrap();
}
