//! Cut cells (L13 stage 2a/2b): closure, volume and area of a cut sphere
//! against the analytic values, exact volumes for flat bodies, split cells,
//! small-cell merging and reproducibility.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::cut::{CutOptions, ImplicitBody, Sphere, cut, edge_root};
use foehncraft_mesh::generate::{BoxSpec, Kind};
use std::f64::consts::PI;

fn box_mesh(kind: Kind, n: usize) -> foehncraft_mesh::Mesh {
    let mut s = BoxSpec::new(kind, 3, [n, n, n]);
    s.lower = [-1.0; 3];
    s.upper = [1.0; 3];
    s.generate().unwrap()
}

fn names() -> Vec<String> {
    vec!["body".to_string()]
}

/// A slab |x − x0| ≤ w/2: flat faces, so linear crossings are exact.
struct Slab {
    x0: f64,
    w: f64,
}

impl ImplicitBody for Slab {
    fn value(&self, x: [f64; 3]) -> f64 {
        (x[0] - self.x0).abs() - 0.5 * self.w
    }
    fn edge_point(&self, a: [f64; 3], b: [f64; 3], fa: f64, fb: f64) -> [f64; 3] {
        edge_root(|x| self.value(x), a, b, fa, fb)
    }
}

/// Fluid only near two opposite corners of a cell: its fluid falls into
/// two pieces.
struct Corners;

impl ImplicitBody for Corners {
    fn value(&self, x: [f64; 3]) -> f64 {
        let d = |c: f64| ((x[0] - c).powi(2) + (x[1] - c).powi(2) + (x[2] - c).powi(2)).sqrt();
        0.8 - d(-1.0).min(d(1.0))
    }
}

#[test]
fn sphere_closure_volume_and_area_converge() {
    let r = 0.55;
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly] {
        let mut ev = Vec::new();
        let mut ea = Vec::new();
        for n in [8, 16, 32] {
            let m = box_mesh(kind, n);
            let body = Sphere {
                centre: [0.03, -0.02, 0.01],
                radius: r,
                exact: true,
            };
            let c = cut(&m, &body, &names(), &CutOptions::default()).unwrap();
            c.mesh.validate().unwrap();
            let rep = &c.report;
            assert!(
                rep.max_closure < 1e-12,
                "{kind:?} n={n}: closure {}",
                rep.max_closure
            );
            assert_eq!(rep.bad_faces, 0, "{kind:?} n={n}: {rep:?}");
            let v = 8.0 - 4.0 / 3.0 * PI * r * r * r;
            ev.push((rep.fluid_volume - v).abs() / (4.0 / 3.0 * PI * r * r * r));
            ea.push((rep.wall_area[0] - 4.0 * PI * r * r).abs() / (4.0 * PI * r * r));
            // every cell has a positive volume
            let g = Geometry::compute(&c.mesh);
            assert!(g.cell_volume.iter().all(|&v| v > 0.0));
        }
        for i in 1..ev.len() {
            let pv = (ev[i - 1] / ev[i]).log2();
            let pa = (ea[i - 1] / ea[i]).log2();
            assert!(pv > 1.6, "{kind:?}: volume order {pv} ({ev:?})");
            assert!(pa > 1.6, "{kind:?}: area order {pa} ({ea:?})");
        }
    }
}

#[test]
fn flat_body_is_exact_and_splits_cells() {
    let m = box_mesh(Kind::Hex, 5);
    let body = Slab { x0: 0.1, w: 0.5 };
    let c = cut(
        &m,
        &body,
        &names(),
        &CutOptions {
            min_fraction: 0.0,
            wall_point: false,
            face_point: false,
        },
    )
    .unwrap();
    let rep = &c.report;
    assert!(
        (rep.fluid_volume - (8.0 - 0.5 * 4.0)).abs() < 1e-12,
        "{rep:?}"
    );
    assert!((rep.wall_area[0] - 8.0).abs() < 1e-12, "{rep:?}");
    assert!(rep.max_closure < 1e-12);
    // one cell whose fluid is two corner pieces becomes two cells
    let one = box_mesh(Kind::Hex, 1);
    let c = cut(
        &one,
        &Corners,
        &names(),
        &CutOptions {
            min_fraction: 0.0,
            wall_point: false,
            face_point: false,
        },
    )
    .unwrap();
    assert_eq!(c.report.cells_split, 1, "{:?}", c.report);
    assert_eq!(c.mesh.n_cells, 2);
    assert!(c.report.max_closure < 1e-12);
}

#[test]
fn small_cells_merge() {
    let m = box_mesh(Kind::Hex, 10);
    let body = Sphere {
        centre: [0.013, 0.021, -0.007],
        radius: 0.43,
        exact: true,
    };
    let none = cut(
        &m,
        &body,
        &names(),
        &CutOptions {
            min_fraction: 0.0,
            wall_point: false,
            face_point: false,
        },
    )
    .unwrap();
    let merged = cut(
        &m,
        &body,
        &names(),
        &CutOptions {
            min_fraction: 0.3,
            wall_point: false,
            face_point: false,
        },
    )
    .unwrap();
    assert!(merged.report.cells_merged > 0);
    assert_eq!(
        merged.report.cells_out + merged.report.cells_merged,
        none.report.cells_out
    );
    assert!(merged.report.min_fraction >= 0.3 || merged.report.merge_failed > 0);
    // merging moves no volume and no wall
    assert!((merged.report.fluid_volume - none.report.fluid_volume).abs() < 1e-12);
    assert!((merged.report.wall_area[0] - none.report.wall_area[0]).abs() < 1e-12);
    assert!(merged.report.max_closure < 1e-12);
    merged.mesh.validate().unwrap();
}

#[test]
fn wall_point_follows_the_sphere() {
    let m = box_mesh(Kind::Tet, 12);
    let r = 0.5;
    let err = |wall_point: bool, exact: bool| {
        let body = Sphere {
            centre: [0.0; 3],
            radius: r,
            exact,
        };
        let c = cut(
            &m,
            &body,
            &names(),
            &CutOptions {
                min_fraction: 0.2,
                wall_point,
                face_point: false,
            },
        )
        .unwrap();
        assert!(c.report.max_closure < 1e-12);
        (c.report.fluid_volume - (8.0 - 4.0 / 3.0 * PI * r * r * r)).abs()
    };
    let lin = err(false, false);
    let ex = err(false, true);
    let wp = err(true, true);
    assert!(ex < lin, "exact crossings {ex} vs linear {lin}");
    assert!(wp < ex, "wall point {wp} vs plain {ex}");
}

#[test]
fn cut_is_independent_of_threads() {
    let m = box_mesh(Kind::Poly, 10);
    let body = Sphere {
        centre: [0.1, 0.0, 0.0],
        radius: 0.6,
        exact: true,
    };
    let run = |t: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build()
            .unwrap()
            .install(|| cut(&m, &body, &names(), &CutOptions::default()).unwrap())
    };
    let (a, b) = (run(1), run(4));
    assert_eq!(a.mesh.points, b.mesh.points);
    assert_eq!(a.mesh.face_nodes, b.mesh.face_nodes);
    assert_eq!(a.mesh.owner, b.mesh.owner);
}

struct Cyl;
impl ImplicitBody for Cyl {
    fn value(&self, x: [f64; 3]) -> f64 {
        ((x[0] - 0.2).powi(2) + (x[1] - 0.2).powi(2)).sqrt() - 0.05
    }
    fn edge_point(&self, a: [f64; 3], b: [f64; 3], fa: f64, fb: f64) -> [f64; 3] {
        edge_root(|x| self.value(x), a, b, fa, fb)
    }
    fn project(&self, x: [f64; 3]) -> Option<[f64; 3]> {
        foehncraft_mesh::cut::project_newton(|y| self.value(y), x, 1e-7)
    }
}

/// The cylinder of Schäfer & Turek on the 2D box mesh of h = D/20: node
/// (0.16, 0.17) lies on the surface to round-off, which once gave crossing
/// points 1e-17 m from a node and a wall triangle the run refused.
#[test]
fn node_on_the_surface_is_snapped() {
    let mut s = BoxSpec::new(Kind::Hex, 2, [440, 82, 1]);
    s.upper = [2.2, 0.41, 0.01];
    let m = s.generate().unwrap();
    let c = cut(
        &m,
        &Cyl,
        &names(),
        &CutOptions {
            min_fraction: 0.3,
            wall_point: true,
            face_point: false,
        },
    )
    .unwrap();
    assert_eq!(c.report.bad_faces, 0, "{:?}", c.report);
    assert_eq!(c.report.merge_failed, 0, "{:?}", c.report);
    assert!(c.report.max_closure < 1e-12);
    let a = 2.0 * PI * 0.05 * 0.01;
    assert!(
        (c.report.wall_area[0] - a).abs() < 1e-3 * a,
        "{:?}",
        c.report
    );
}
