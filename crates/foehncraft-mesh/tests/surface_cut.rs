//! L13 stage 0 (D225 item 5): cut-cell moments of a mesh cut by a
//! triangle surface (`Surface::cut_moments`).

#![allow(clippy::needless_range_loop)]

use crate::surface::{icosphere, torus};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::format::stl::Triangle;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::surface::{Surface, check};

fn volume_area_centroid(tris: &[Triangle]) -> (f64, f64, [f64; 3]) {
    let (mut v, mut a, mut m) = (0.0, 0.0, [0.0; 3]);
    for t in tris {
        let [p, q, r] = *t;
        let c = [
            q[1] * r[2] - q[2] * r[1],
            q[2] * r[0] - q[0] * r[2],
            q[0] * r[1] - q[1] * r[0],
        ];
        let tv = (p[0] * c[0] + p[1] * c[1] + p[2] * c[2]) / 6.0;
        v += tv;
        for k in 0..3 {
            m[k] += tv * (p[k] + q[k] + r[k]) / 4.0;
        }
        let u = [0, 1, 2].map(|k| q[k] - p[k]);
        let w = [0, 1, 2].map(|k| r[k] - p[k]);
        let x = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        a += 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    }
    (v, a, m.map(|x| x / v))
}

fn check_cut(kind: Kind, cells: [usize; 3], tris: &[Triangle]) {
    let mesh = BoxSpec::new(kind, 3, cells).generate().unwrap();
    let geo = Geometry::compute(&mesh);
    let s = Surface::new(tris);
    assert!(check(&s).clean());
    let cut = s.cut_moments(&mesh, &geo);
    let (v, a, cen) = volume_area_centroid(tris);
    let total: f64 = cut.cell_inside_volume.iter().sum();
    assert!(
        (total - v).abs() < 1e-12 * v,
        "{kind:?}: volume {total} against {v}"
    );
    let mut first = [0.0; 3];
    for c in 0..mesh.n_cells {
        let vc = cut.cell_inside_volume[c];
        assert!(
            vc >= -1e-14 && vc <= geo.cell_volume[c] * (1.0 + 1e-12),
            "cell {c}: {vc}"
        );
        for k in 0..3 {
            first[k] += vc * cut.cell_inside_centroid[c][k];
        }
    }
    for k in 0..3 {
        assert!(
            (first[k] / total - cen[k]).abs() < 1e-11,
            "{kind:?}: centroid"
        );
    }
    let wall: f64 = cut
        .walls
        .iter()
        .map(|w| w.area.iter().map(|x| x * x).sum::<f64>().sqrt())
        .sum();
    assert!(
        (wall - a).abs() < 1e-11 * a,
        "{kind:?}: wall {wall} against {a}"
    );
    // closure of every cell: Σ inside face areas · n_out + Σ wall area = 0
    // (checked where the faces are planar; a warped polygon is split into
    // fan triangles with their own normals)
    if !matches!(kind, Kind::Hex | Kind::Tet) {
        return;
    }
    let mut worst: f64 = 0.0;
    let mut closure = vec![[0.0f64; 3]; mesh.n_cells];
    for w in &cut.walls {
        for k in 0..3 {
            closure[w.cell][k] += w.area[k];
        }
    }
    for c in 0..mesh.n_cells {
        for (f, nb) in mesh.cell_face_iter(c) {
            let s = geo.face_area[f];
            let full = s.iter().map(|x| x * x).sum::<f64>().sqrt();
            let frac = cut.face_inside_area[f] / full;
            let sign = if nb { -1.0 } else { 1.0 };
            for k in 0..3 {
                closure[c][k] += sign * frac * s[k];
            }
        }
        let scale = geo.cell_volume[c].powf(2.0 / 3.0);
        worst = worst.max(closure[c].iter().map(|x| x.abs()).fold(0.0, f64::max) / scale);
    }
    assert!(worst < 1e-12, "{kind:?}: closure {worst}");
}

#[test]
fn cut_sphere_and_torus_on_every_cell_type() {
    let sphere: Vec<Triangle> = icosphere([0.47, 0.52, 0.49], 0.31, 3);
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        check_cut(kind, [7, 6, 5], &sphere);
    }
    let t: Vec<Triangle> = torus(0.28, 0.1, 48, 24)
        .iter()
        .map(|t| t.map(|p| [p[0] + 0.51, p[1] + 0.48, p[2] + 0.513]))
        .collect::<Vec<Triangle>>()
        .into_iter()
        .collect();
    check_cut(Kind::Hex, [9, 9, 4], &t);
    // a ring of surface edges exactly in the mesh plane z = 0.5
    let t: Vec<Triangle> = torus(0.28, 0.1, 48, 24)
        .iter()
        .map(|t| t.map(|p| [p[0] + 0.51, p[1] + 0.48, p[2] + 0.5]))
        .collect();
    check_cut(Kind::Hex, [9, 9, 4], &t);
}

#[test]
fn cut_moments_are_the_same_on_any_thread_count() {
    let mesh = BoxSpec::new(Kind::Poly, 3, [6, 6, 6]).generate().unwrap();
    let geo = Geometry::compute(&mesh);
    let s = Surface::new(&icosphere([0.5, 0.5, 0.5], 0.33, 3));
    let run = |n: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| s.cut_moments(&mesh, &geo))
    };
    let a = run(1);
    for n in [2, 4] {
        let b = run(n);
        assert!(
            a.cell_inside_volume
                .iter()
                .zip(&b.cell_inside_volume)
                .all(|(x, y)| x.to_bits() == y.to_bits())
        );
        assert!(
            a.face_inside_area
                .iter()
                .zip(&b.face_inside_area)
                .all(|(x, y)| x.to_bits() == y.to_bits())
        );
        assert_eq!(a.walls, b.walls);
    }
}
