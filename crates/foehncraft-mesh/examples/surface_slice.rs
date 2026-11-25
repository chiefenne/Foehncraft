//! A dirty STL body and its winding number and signed distance on a slice
//! (L13 stage 0, D225; gallery entry of `docs/gallery/immersed.md`).
//!
//! ```sh
//! cargo run --release -p foehncraft-mesh --example surface_slice -- output/surface
//! python3 -m tools.plot.surface_slice output/surface docs/figures/gallery/immersed_surface_winding.png
//! ```
//!
//! The body is a lumpy ball (an icosphere with a smooth radial bump field)
//! with three defects: a hole in its top, a patch of flipped triangles and a
//! second, smaller lump that overlaps it (self-intersections). Written: the
//! surface (`body.stl`), the surface check (`check.json`), the slice values
//! (`slice.csv`: x, z, fast and exact winding number, signed distance) and
//! the trace of the surface in the slice plane (`trace.csv`).

use foehncraft_mesh::format::stl::{Triangle, to_binary};
use foehncraft_mesh::surface::{Surface, check};
use std::io::Write;
use std::time::Instant;

fn icosphere(levels: usize) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
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
    (v, f)
}

/// A lumpy ball: unit directions scaled by a smooth bump field.
fn lump(c: [f64; 3], r: f64, levels: usize, phase: f64) -> Vec<Triangle> {
    let (v, f) = icosphere(levels);
    let p: Vec<[f64; 3]> = v
        .iter()
        .map(|u| {
            let s = 1.0
                + 0.16 * (2.1 * u[0] + 0.4 + phase).sin() * (1.7 * u[1] - 0.3).cos()
                + 0.10 * (3.3 * u[2] + 1.1 - phase).sin()
                + 0.06 * (4.1 * u[0] * u[1] + 0.8).cos();
            [0, 1, 2].map(|d| c[d] + r * s * u[d])
        })
        .collect();
    f.iter().map(|t| t.map(|i| p[i])).collect()
}

fn main() {
    let out = std::path::PathBuf::from(
        std::env::args()
            .nth(1)
            .unwrap_or_else(|| "output/surface".into()),
    );
    std::fs::create_dir_all(&out).expect("output directory");

    let main_body = lump([0.0; 3], 1.0, 5, 0.0);
    let centroid = |t: &Triangle| [0, 1, 2].map(|d| (t[0][d] + t[1][d] + t[2][d]) / 3.0);
    // hole: a cap around the direction h
    let h = [0.35, 0.0, 0.94];
    let mut tris: Vec<Triangle> = main_body
        .into_iter()
        .filter(|t| {
            let c = centroid(t);
            let n = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
            (c[0] * h[0] + c[1] * h[1] + c[2] * h[2]) / n < 0.985
        })
        .collect();
    // flipped patch on the lower left
    for t in tris.iter_mut() {
        let c = centroid(t);
        if c[0] < -0.55 && c[2] < -0.45 {
            *t = [t[0], t[2], t[1]];
        }
    }
    // an overlapping second lump on the right
    tris.extend(lump([0.95, 0.05, -0.35], 0.48, 4, 1.3));
    std::fs::write(out.join("body.stl"), to_binary(&tris)).expect("body.stl");

    let t0 = Instant::now();
    let s = Surface::new(&tris);
    let build = t0.elapsed().as_secs_f64();
    let report = check(&s);
    std::fs::write(
        out.join("check.json"),
        serde_json::to_string_pretty(&report).unwrap(),
    )
    .expect("check.json");
    for f in &report.findings {
        eprintln!("{} {}", f.code, f.message);
    }

    // slice y = 0, x in [-1.8, 1.9], z in [-1.6, 1.6]
    let (nx, nz) = (370usize, 320usize);
    let pts: Vec<[f64; 3]> = (0..nz)
        .flat_map(|j| {
            (0..nx).map(move |i| {
                [
                    -1.8 + 3.7 * (i as f64 + 0.5) / nx as f64,
                    0.0,
                    -1.6 + 3.2 * (j as f64 + 0.5) / nz as f64,
                ]
            })
        })
        .collect();
    let t1 = Instant::now();
    let fast = s.winding_numbers(&pts);
    let t_fast = t1.elapsed().as_secs_f64();
    let t2 = Instant::now();
    let exact: Vec<f64> = {
        use rayon::prelude::*;
        pts.par_iter().map(|&q| s.winding_number_exact(q)).collect()
    };
    let t_exact = t2.elapsed().as_secs_f64();
    let t3 = Instant::now();
    let phi = s.signed_distances(&pts, None);
    let t_phi = t3.elapsed().as_secs_f64();
    let mut f = std::io::BufWriter::new(std::fs::File::create(out.join("slice.csv")).unwrap());
    writeln!(f, "x,z,winding,winding_exact,phi").unwrap();
    for k in 0..pts.len() {
        writeln!(
            f,
            "{:.6},{:.6},{:.9e},{:.9e},{:.9e}",
            pts[k][0], pts[k][2], fast[k], exact[k], phi[k]
        )
        .unwrap();
    }
    // trace of the triangles in the plane y = 0
    let mut g = std::io::BufWriter::new(std::fs::File::create(out.join("trace.csv")).unwrap());
    writeln!(g, "x0,z0,x1,z1,flipped").unwrap();
    for (ti, t) in tris.iter().enumerate() {
        let mut ends = Vec::new();
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if (a[1] > 0.0) != (b[1] > 0.0) {
                let w = a[1] / (a[1] - b[1]);
                ends.push([a[0] + w * (b[0] - a[0]), a[2] + w * (b[2] - a[2])]);
            }
        }
        if ends.len() == 2 {
            let c = centroid(t);
            let flipped = ti < tris.len() && c[0] < -0.55 && c[2] < -0.45 && c[0] > -1.6;
            writeln!(
                g,
                "{},{},{},{},{}",
                ends[0][0], ends[0][1], ends[1][0], ends[1][1], flipped as u8
            )
            .unwrap();
        }
    }
    let worst = fast
        .iter()
        .zip(&exact)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let summary = serde_json::json!({
        "triangles": s.len(),
        "points": pts.len(),
        "build_s": build,
        "fast_winding_s": t_fast,
        "exact_winding_s": t_exact,
        "signed_distance_s": t_phi,
        "fast_vs_exact_max": worst,
        "threads": rayon::current_num_threads(),
    });
    std::fs::write(
        out.join("summary.json"),
        serde_json::to_string_pretty(&summary).unwrap(),
    )
    .unwrap();
    println!("{summary}");
}
