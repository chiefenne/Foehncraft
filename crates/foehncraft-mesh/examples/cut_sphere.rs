//! A sphere cut out of box meshes (L13 stage 0, D225 item 5; gallery entry
//! of `docs/gallery/immersed.md`).
//!
//! ```sh
//! cargo run --release -p foehncraft-mesh --example cut_sphere -- output/cut
//! python3 -m tools.plot.cut_sphere output/cut docs/figures/gallery/immersed_cut_sphere.png
//! ```
//!
//! The sphere (radius 0.3 about (0.51, 0.47, 0.5) in the unit cube) is an
//! icosphere of 7 subdivisions (327 680 triangles, facet error 1.4e-5 of
//! the radius). On hexahedral meshes of 8 to 64 cells per side the exact cut
//! (`Surface::cut_moments`) is compared with the cut a level set gives: the
//! exact signed distance of the sphere at the mesh nodes, interpolated
//! linearly on the six tetrahedra of each cube (marching tetrahedra volume).
//! The reference is the analytic sphere through a sub-grid integration
//! of each cut cell (linear level set on 16³ and 32³ sub-cubes,
//! Richardson-extrapolated); the exact cut of the icosphere differs from
//! the sphere by the facet error only. Written: `convergence.json` and, for the 24³ mesh,
//! `slice.csv` (cells of the slice through the centre: fraction, face
//! aperture in x) and `walls.csv` (wall piece centroids near the slice).

use foehncraft_mesh::Geometry;
use foehncraft_mesh::format::stl::Triangle;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::surface::Surface;
use rayon::prelude::*;
use std::io::Write;
use std::time::Instant;

const C: [f64; 3] = [0.51, 0.47, 0.5];
const R: f64 = 0.3;

fn icosphere(levels: usize) -> Vec<Triangle> {
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
        .map(|tri| tri.map(|i| [0, 1, 2].map(|d| C[d] + R * v[i][d])))
        .collect()
}

fn phi(x: [f64; 3]) -> f64 {
    ((x[0] - C[0]).powi(2) + (x[1] - C[1]).powi(2) + (x[2] - C[2]).powi(2)).sqrt() - R
}

/// Volume of `{φ < 0}` in a tetrahedron with φ linear between its corners.
fn tet_inside(p: [[f64; 3]; 4], f: [f64; 4]) -> f64 {
    let vol = |a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]| {
        let u = [0, 1, 2].map(|k| b[k] - a[k]);
        let v = [0, 1, 2].map(|k| c[k] - a[k]);
        let w = [0, 1, 2].map(|k| d[k] - a[k]);
        (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0]))
            .abs()
            / 6.0
    };
    let full = vol(p[0], p[1], p[2], p[3]);
    let neg: Vec<usize> = (0..4).filter(|&i| f[i] < 0.0).collect();
    let pos: Vec<usize> = (0..4).filter(|&i| f[i] >= 0.0).collect();
    let cut = |a: usize, b: usize| {
        let t = f[a] / (f[a] - f[b]);
        [0, 1, 2].map(|k| p[a][k] + t * (p[b][k] - p[a][k]))
    };
    match neg.len() {
        0 => 0.0,
        4 => full,
        1 => {
            let a = neg[0];
            vol(p[a], cut(a, pos[0]), cut(a, pos[1]), cut(a, pos[2]))
        }
        3 => {
            let b = pos[0];
            full - vol(p[b], cut(b, neg[0]), cut(b, neg[1]), cut(b, neg[2]))
        }
        _ => {
            // prism between the two negative corners: split into 3 tets
            let (a, b) = (neg[0], neg[1]);
            let (c, d) = (pos[0], pos[1]);
            let (ac, ad, bc, bd) = (cut(a, c), cut(a, d), cut(b, c), cut(b, d));
            vol(p[a], p[b], ac, ad) + vol(p[b], ac, ad, bd) + vol(p[b], ac, bc, bd)
        }
    }
}

/// The analytic sphere's volume in the box `lo`–`lo + h` by a sub-grid of
/// `s³` boxes, each with the linear level-set volume of its six tets (the
/// error of the sub-grid falls like (h/s)²).
fn reference(lo: [f64; 3], h: f64, s: usize) -> f64 {
    let g = h / s as f64;
    let kuhn = [
        [0, 1, 3, 7],
        [0, 3, 2, 7],
        [0, 2, 6, 7],
        [0, 6, 4, 7],
        [0, 4, 5, 7],
        [0, 5, 1, 7],
    ];
    let mut v = 0.0;
    for i in 0..s {
        for j in 0..s {
            for k in 0..s {
                let o = [
                    lo[0] + i as f64 * g,
                    lo[1] + j as f64 * g,
                    lo[2] + k as f64 * g,
                ];
                let corner = |c: usize| {
                    [
                        o[0] + (c & 1) as f64 * g,
                        o[1] + (c >> 1 & 1) as f64 * g,
                        o[2] + (c >> 2 & 1) as f64 * g,
                    ]
                };
                for t in kuhn {
                    let p = t.map(corner);
                    v += tet_inside(p, p.map(phi));
                }
            }
        }
    }
    v
}

fn main() {
    let out = std::path::PathBuf::from(std::env::args().nth(1).unwrap_or("output/cut".into()));
    std::fs::create_dir_all(&out).unwrap();
    let tris = icosphere(7);
    let s = Surface::new(&tris);
    let mut rows = Vec::new();
    for n in [8usize, 16, 32, 64] {
        let mesh = BoxSpec::new(Kind::Hex, 3, [n, n, n]).generate().unwrap();
        let geo = Geometry::compute(&mesh);
        let t0 = Instant::now();
        let cut = s.cut_moments(&mesh, &geo);
        let t_cut = t0.elapsed().as_secs_f64();
        let h = 1.0 / n as f64;
        // the cut cells: neither empty nor full by the exact cut
        let cells: Vec<usize> = (0..mesh.n_cells)
            .filter(|&c| {
                let f = cut.cell_inside_volume[c] / geo.cell_volume[c];
                f > 1e-12 && f < 1.0 - 1e-12
            })
            .collect();
        let errs: Vec<(f64, f64)> = cells
            .par_iter()
            .map(|&c| {
                let nodes = mesh.cell_nodes(c);
                let mut lo = [f64::INFINITY; 3];
                for &v in &nodes {
                    for k in 0..3 {
                        lo[k] = lo[k].min(mesh.points[v as usize][k]);
                    }
                }
                // Richardson: the sub-grid error falls like (h/s)²
                let refv = (4.0 * reference(lo, h, 32) - reference(lo, h, 16)) / 3.0;
                let lin = reference(lo, h, 1);
                let vol = h * h * h;
                (
                    (cut.cell_inside_volume[c] - refv).abs() / vol,
                    (lin - refv).abs() / vol,
                )
            })
            .collect();
        let max_cut = errs.iter().map(|e| e.0).fold(0.0, f64::max);
        let max_lin = errs.iter().map(|e| e.1).fold(0.0, f64::max);
        let l1_cut = errs.iter().map(|e| e.0).sum::<f64>() / errs.len() as f64;
        let l1_lin = errs.iter().map(|e| e.1).sum::<f64>() / errs.len() as f64;
        let total: f64 = cut.cell_inside_volume.iter().sum();
        let wall: f64 = cut
            .walls
            .iter()
            .map(|w| w.area.iter().map(|x| x * x).sum::<f64>().sqrt())
            .sum();
        eprintln!(
            "n={n}: {} cut cells, fraction error max {max_cut:.2e} / {max_lin:.2e} (exact cut / linear φ), mean {l1_cut:.2e} / {l1_lin:.2e}, {t_cut:.2} s",
            cells.len()
        );
        rows.push(serde_json::json!({
            "n": n, "cut_cells": cells.len(), "cut_seconds": t_cut,
            "max_exact": max_cut, "max_linear": max_lin, "mean_exact": l1_cut, "mean_linear": l1_lin,
            "total_volume": total, "sphere_volume": 4.0 / 3.0 * std::f64::consts::PI * R.powi(3),
            "wall_area": wall, "sphere_area": 4.0 * std::f64::consts::PI * R * R,
        }));
    }
    std::fs::write(
        out.join("convergence.json"),
        serde_json::to_string_pretty(&serde_json::json!({"triangles": tris.len(), "rows": rows}))
            .unwrap(),
    )
    .unwrap();

    // slice data on 24³
    let n = 24;
    let mesh = BoxSpec::new(Kind::Hex, 3, [n, n, n]).generate().unwrap();
    let geo = Geometry::compute(&mesh);
    let cut = s.cut_moments(&mesh, &geo);
    let k = ((C[2] * n as f64) as usize).min(n - 1);
    let mut f = std::io::BufWriter::new(std::fs::File::create(out.join("slice.csv")).unwrap());
    writeln!(f, "x0,y0,h,fraction,xface_aperture").unwrap();
    let h = 1.0 / n as f64;
    for c in 0..mesh.n_cells {
        let cc = geo.cell_centroid[c];
        if ((cc[2] * n as f64) as usize) != k {
            continue;
        }
        // aperture of the face on the low x side
        let mut ap = 0.0;
        for (fc, _) in mesh.cell_face_iter(c) {
            let fcen = geo.face_centroid[fc];
            if (fcen[0] - (cc[0] - 0.5 * h)).abs() < 1e-9 {
                let a = geo.face_area[fc].iter().map(|x| x * x).sum::<f64>().sqrt();
                ap = cut.face_inside_area[fc] / a;
            }
        }
        writeln!(
            f,
            "{},{},{},{},{}",
            cc[0] - 0.5 * h,
            cc[1] - 0.5 * h,
            h,
            cut.cell_inside_volume[c] / geo.cell_volume[c],
            ap
        )
        .unwrap();
    }
    let mut w = std::io::BufWriter::new(std::fs::File::create(out.join("walls.csv")).unwrap());
    writeln!(w, "x,y,z,ax,ay,az").unwrap();
    for p in &cut.walls {
        let cc = geo.cell_centroid[p.cell];
        if ((cc[2] * n as f64) as usize) == k {
            writeln!(
                w,
                "{},{},{},{},{},{}",
                p.centroid[0], p.centroid[1], p.centroid[2], p.area[0], p.area[1], p.area[2]
            )
            .unwrap();
        }
    }
    println!("slice z cells {k}, wrote {}", out.display());
}
