//! The wall faces of a box turned against the mesh and cut out of it (L13
//! stage 2, gallery 320): with the closest point of the surface as the fan
//! point (the edges and corners are rounded off inside the cut cells), and
//! with the dual-contouring wall point and the face points (kept sharp).
//! Writes the section of the wall faces with the plane z = `Z0` (and of the
//! box triangles, the exact section) as JSON for
//! `tools/plot/immersed_stl.py`.
//!
//! ```sh
//! cargo run --release -p foehncraft-mesh --example cut_stl_walls -- output/cut_stl_walls.json
//! ```

use foehncraft_mesh::Geometry;
use foehncraft_mesh::cut::{CutOptions, cut};
use foehncraft_mesh::format::stl::Triangle;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::surface::{Surface, SurfaceBody};
use std::fmt::Write as _;

/// The box of the benchmark (`cut_stl_bench`).
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

/// The plane of the section, between two node planes of the mesh.
const Z0: f64 = 0.05;

/// The pieces of polygon `ring` in the plane z = `Z0`: its crossing
/// points in ring order, paired along the line they lie on.
fn section(ring: &[[f64; 3]]) -> Vec<[[f64; 2]; 2]> {
    let k = ring.len();
    let mut x: Vec<[f64; 2]> = Vec::new();
    for i in 0..k {
        let (a, b) = (ring[i], ring[(i + 1) % k]);
        let (da, db) = (a[2] - Z0, b[2] - Z0);
        if (da < 0.0) != (db < 0.0) {
            let t = da / (da - db);
            x.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
        }
    }
    if x.len() < 2 {
        return Vec::new();
    }
    let d = [x[x.len() - 1][0] - x[0][0], x[x.len() - 1][1] - x[0][1]];
    x.sort_by(|p, q| (p[0] * d[0] + p[1] * d[1]).total_cmp(&(q[0] * d[0] + q[1] * d[1])));
    x.chunks_exact(2).map(|c| [c[0], c[1]]).collect()
}

fn segments(s: &[[[f64; 2]; 2]]) -> String {
    let v: Vec<String> = s
        .iter()
        .map(|[a, b]| format!("[[{:e},{:e}],[{:e},{:e}]]", a[0], a[1], b[0], b[1]))
        .collect();
    format!("[{}]", v.join(","))
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "output/cut_stl_walls.json".into());
    let n = 16;
    let tris = rotated_box([0.021, -0.013, 0.017], [0.55, 0.4, 0.3], 0.35, 0.6);
    let surface = Surface::new(&tris);
    let mut s = BoxSpec::new(Kind::Hex, 3, [n, n, n]);
    s.lower = [-1.0; 3];
    s.upper = [1.0; 3];
    let mesh = s.generate().expect("mesh");
    let h = 2.0 / n as f64;
    let mut json = String::from("{");
    let exact: Vec<_> = tris.iter().flat_map(|t| section(t)).collect();
    writeln!(
        json,
        "\"h\": {h:e}, \"z0\": {Z0:e}, \"exact\": {},",
        segments(&exact)
    )
    .unwrap();
    let cases = [("closest", false, false), ("sharp", true, true)];
    for (k, (name, features, face_point)) in cases.into_iter().enumerate() {
        let body = SurfaceBody {
            surface: &surface,
            exact: true,
            features,
        };
        let c = cut(
            &mesh,
            &body,
            &["box".to_string()],
            &CutOptions {
                min_fraction: 0.3,
                wall_point: true,
                face_point,
            },
        )
        .expect("cut");
        let geo = Geometry::compute(&c.mesh);
        let b = c.mesh.boundary("box").expect("wall boundary");
        let mut walls = Vec::new();
        let mut dist = Vec::new();
        for f in b.range() {
            let ring: Vec<[f64; 3]> = c
                .mesh
                .face(f)
                .iter()
                .map(|&i| c.mesh.points[i as usize])
                .collect();
            for seg in section(&ring) {
                walls.push(seg);
                dist.push(format!("{:e}", surface.distance(geo.face_centroid[f]) / h));
            }
        }
        // the other faces of the cut mesh: the cells in the plane
        let mut cells = Vec::new();
        for f in 0..c.mesh.n_faces() {
            if b.range().contains(&f) {
                continue;
            }
            let ring: Vec<[f64; 3]> = c
                .mesh
                .face(f)
                .iter()
                .map(|&i| c.mesh.points[i as usize])
                .collect();
            cells.extend(section(&ring));
        }
        eprintln!(
            "{name}: {} wall faces, wall area {:.6}, bad faces {}, simplified {}",
            b.len, c.report.wall_area[0], c.report.bad_faces, c.report.simplified
        );
        let sep = if k + 1 < cases.len() { "," } else { "" };
        writeln!(
            json,
            "\"{name}\": {{\"area\": {:e}, \"walls\": {}, \"distance\": [{}], \"cells\": {}}}{sep}",
            c.report.wall_area[0],
            segments(&walls),
            dist.join(","),
            segments(&cells)
        )
        .unwrap();
    }
    json.push('}');
    std::fs::write(&path, json).expect("write");
    eprintln!("wrote {path}");
}
