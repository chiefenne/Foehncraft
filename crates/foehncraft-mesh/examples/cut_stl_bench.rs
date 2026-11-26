//! Benchmark of cut cells for STL bodies (L13 stage 2, maintainer rule:
//! benchmark the computer-graphics algorithms too), against the exact cut
//! of the same triangles (`Surface::cut_moments`, Gauss–Green by dimension,
//! stage 0).
//!
//! Three surfaces cut out of hexahedral, tetrahedral and polyhedral box
//! meshes of 8³ to 32³ cells on [−1, 1]³:
//!
//! - `sphere`: an icosphere of 20 480 triangles (smooth, convex);
//! - `box`: a rotated box of 12 triangles (flat sides, sharp edges and
//!   corners);
//! - `torus`: 24 576 triangles (convex and concave curvature).
//!
//! Six ways to place the cut surface:
//!
//! - `MC`: crossing points at the zero of the linear interpolant of the node
//!   distances, one wall face per loop (marching cubes / tetrahedra);
//! - `exact`: crossing points at the exact segment–triangle crossings;
//! - `closest`: as `exact`, wall faces fanned around the closest point of
//!   the surface to the loop centroid;
//! - `DC`: as `closest`, but where the normals of the loop open more than
//!   15° the fan point is the vertex of dual contouring (QEF of the loop
//!   points and normals, Ju et al. 2002), kept when it lies on the surface;
//! - `closest+face`, `DC+face`: every chord of a cut face bent through a
//!   point of the surface in the face plane (the bisector crossing, or with
//!   `DC` the point where a sharp edge pierces the face).
//!
//! Errors against the exact cut of the same mesh: the fluid volume and the
//! wall area (relative to the body), the cell volumes (L1 over the groups
//! of input cells that make up the cut cells, relative to the body volume;
//! max relative to the group's uncut volume) and the face apertures (fluid
//! area of every input face: mean and max of the error over the cut faces,
//! relative to the face area). Printed as a Markdown table; with an argument
//! also written as JSON for `tools/plot/immersed_stl.py`.
//!
//! ```sh
//! cargo run --release -p foehncraft-mesh --example cut_stl_bench -- output/cut_stl_bench.json
//! ```

use foehncraft_mesh::Geometry;
use foehncraft_mesh::cut::{CutOptions, cut};
use foehncraft_mesh::format::stl::Triangle;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::surface::{Surface, SurfaceBody};
use std::io::Write;
use std::time::Instant;

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

fn torus(c: [f64; 3], big: f64, small: f64, n: usize, m: usize) -> Vec<Triangle> {
    // tilted by 0.4 rad about x so that no mesh plane is a symmetry plane
    let (ca, sa) = (0.4f64.cos(), 0.4f64.sin());
    let p = |i: usize, j: usize| {
        let u = std::f64::consts::TAU * (i % n) as f64 / n as f64;
        let v = std::f64::consts::TAU * (j % m) as f64 / m as f64;
        let q = [
            (big + small * v.cos()) * u.cos(),
            (big + small * v.cos()) * u.sin(),
            small * v.sin(),
        ];
        [
            c[0] + q[0],
            c[1] + ca * q[1] - sa * q[2],
            c[2] + sa * q[1] + ca * q[2],
        ]
    };
    let mut out = Vec::new();
    for i in 0..n {
        for j in 0..m {
            let (a, b, cc, d) = (p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1));
            out.push([a, b, cc]);
            out.push([a, cc, d]);
        }
    }
    out
}

/// A box with half sides `h` about `c`, rotated by `ax` about x, then `az`
/// about z.
fn rotated_box(c: [f64; 3], h: [f64; 3], ax: f64, az: f64) -> Vec<Triangle> {
    let (cx, sx, cz, sz) = (ax.cos(), ax.sin(), az.cos(), az.sin());
    let rot = |p: [f64; 3]| {
        let q = [p[0], cx * p[1] - sx * p[2], sx * p[1] + cx * p[2]];
        [
            c[0] + cz * q[0] - sz * q[1],
            c[1] + sz * q[0] + cz * q[1],
            c[2] + q[2],
        ]
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

struct Uf(Vec<usize>);

impl Uf {
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
}

fn main() {
    let out = std::env::args().nth(1);
    let shapes: Vec<(&str, Vec<Triangle>)> = vec![
        ("sphere", icosphere([0.031, -0.017, 0.012], 0.55, 5)),
        (
            "box",
            rotated_box([0.021, -0.013, 0.017], [0.55, 0.4, 0.3], 0.35, 0.6),
        ),
        ("torus", torus([0.013, 0.021, -0.011], 0.55, 0.22, 192, 64)),
    ];
    let methods: [(&str, bool, bool, bool, bool); 6] = [
        ("MC", false, false, false, false),
        ("exact", true, false, false, false),
        ("closest", true, true, false, false),
        ("DC", true, true, true, false),
        ("closest+face", true, true, false, true),
        ("DC+face", true, true, true, true),
    ];
    println!(
        "| body | mesh | n | method | volume | order | area | order | cells L1 | cells max | apertures mean | apertures max | merged | bad faces | simplified | cut (ms) |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|");
    let mut rows = Vec::new();
    for (sname, tris) in &shapes {
        let surface = Surface::new(tris);
        for kind in [Kind::Hex, Kind::Tet, Kind::Poly] {
            let meshes: Vec<_> = [8usize, 16, 32]
                .iter()
                .map(|&n| {
                    let mut s = BoxSpec::new(kind, 3, [n, n, n]);
                    s.lower = [-1.0; 3];
                    s.upper = [1.0; 3];
                    let m = s.generate().expect("mesh");
                    let g = Geometry::compute(&m);
                    let t = Instant::now();
                    let mo = surface.cut_moments(&m, &g);
                    let ms = t.elapsed().as_secs_f64() * 1e3;
                    (n, m, g, mo, ms)
                })
                .collect();
            for (mname, exact, wall_point, features, face_point) in methods {
                let mut prev: Option<(f64, f64)> = None;
                for (n, m, g, mo, mo_ms) in &meshes {
                    let body = SurfaceBody {
                        surface: &surface,
                        exact,
                        features,
                    };
                    let t = Instant::now();
                    let c = cut(
                        m,
                        &body,
                        &["body".to_string()],
                        &CutOptions {
                            min_fraction: 0.3,
                            wall_point,
                            face_point,
                        },
                    )
                    .expect("cut");
                    let ms = t.elapsed().as_secs_f64() * 1e3;
                    let gc = Geometry::compute(&c.mesh);
                    let v_body: f64 = mo.cell_inside_volume.iter().sum();
                    let v_ex: f64 = g.cell_volume.iter().sum::<f64>() - v_body;
                    let a_ex: f64 = mo
                        .walls
                        .iter()
                        .map(|w| {
                            (w.area[0] * w.area[0] + w.area[1] * w.area[1] + w.area[2] * w.area[2])
                                .sqrt()
                        })
                        .sum();
                    let ev = (c.report.fluid_volume - v_ex).abs() / v_body;
                    let ea = (c.report.wall_area[0] - a_ex).abs() / a_ex;
                    // groups of input cells joined by the cut cells
                    let mut uf = Uf((0..m.n_cells).collect());
                    for cov in &c.covers {
                        for w in cov.windows(2) {
                            let (a, b) = (uf.find(w[0] as usize), uf.find(w[1] as usize));
                            uf.0[a.max(b)] = a.min(b);
                        }
                    }
                    let mut vcut = vec![0.0; m.n_cells];
                    let mut vex = vec![0.0; m.n_cells];
                    let mut vun = vec![0.0; m.n_cells];
                    let mut touched = vec![false; m.n_cells];
                    for (k, cov) in c.covers.iter().enumerate() {
                        let r = uf.find(cov[0] as usize);
                        vcut[r] += gc.cell_volume[k];
                    }
                    for i in 0..m.n_cells {
                        let r = uf.find(i);
                        vex[r] += g.cell_volume[i] - mo.cell_inside_volume[i];
                        vun[r] += g.cell_volume[i];
                        if mo.cell_inside_volume[i] > 0.0 {
                            touched[r] = true;
                        }
                    }
                    let (mut l1, mut cmax) = (0.0f64, 0.0f64);
                    for i in 0..m.n_cells {
                        if uf.find(i) == i && (touched[i] || (vcut[i] - vex[i]).abs() > 0.0) {
                            let d = (vcut[i] - vex[i]).abs();
                            l1 += d;
                            cmax = cmax.max(d / vun[i]);
                        }
                    }
                    let l1 = l1 / v_body;
                    let (mut asum, mut amax, mut acount) = (0.0f64, 0.0f64, 0usize);
                    for f in 0..m.n_faces() {
                        if mo.face_inside_area[f] > 0.0 || c.face_fluid_area[f] == 0.0 {
                            let af = g.face_area[f];
                            let af = (af[0] * af[0] + af[1] * af[1] + af[2] * af[2]).sqrt();
                            let ex = af - mo.face_inside_area[f];
                            if ex >= af {
                                continue;
                            }
                            let d = (c.face_fluid_area[f] - ex).abs() / af;
                            asum += d;
                            amax = amax.max(d);
                            acount += 1;
                        }
                    }
                    let amean = asum / acount.max(1) as f64;
                    let (ov, oa) = match prev {
                        Some((pv, pa)) => (
                            format!("{:.2}", (pv / ev).log2()),
                            format!("{:.2}", (pa / ea).log2()),
                        ),
                        None => (String::new(), String::new()),
                    };
                    prev = Some((ev, ea));
                    println!(
                        "| {sname} | {} | {n} | {mname} | {ev:.2e} | {ov} | {ea:.2e} | {oa} | {l1:.2e} | {cmax:.2e} | {amean:.2e} | {amax:.2e} | {} | {} | {} | {ms:.0} |",
                        kind.name(),
                        c.report.cells_merged,
                        c.report.bad_faces,
                        c.report.simplified
                    );
                    rows.push(format!(
                        "{{\"body\": \"{sname}\", \"mesh\": \"{}\", \"n\": {n}, \"method\": \"{mname}\", \"volume\": {ev:e}, \"area\": {ea:e}, \"cells_l1\": {l1:e}, \"cells_max\": {cmax:e}, \"apertures_mean\": {amean:e}, \"apertures_max\": {amax:e}, \"cells\": {}, \"merged\": {}, \"failed\": {}, \"bad_faces\": {}, \"simplified\": {}, \"rounds\": {}, \"cut_ms\": {ms:.1}, \"moments_ms\": {mo_ms:.1}}}",
                        kind.name(),
                        c.mesh.n_cells,
                        c.report.cells_merged,
                        c.report.merge_failed,
                        c.report.bad_faces,
                        c.report.simplified,
                        c.report.rounds
                    ));
                }
            }
        }
    }
    if let Some(path) = out {
        let mut f = std::fs::File::create(&path).expect("output file");
        writeln!(f, "[\n{}\n]", rows.join(",\n")).expect("write");
    }
}
