//! L13 stage 2 (D227): a triangle surface as the body of a cut. The cut
//! cells against the exact moments of the same surface
//! (`Surface::cut_moments`): fluid volume, wall area, the fluid volume of
//! every group of input cells a merged cell covers, and the aperture of
//! every input face; the sharp-feature wall points on a box; the fallback
//! that keeps every face usable; the same mesh on any thread count.

#![allow(clippy::needless_range_loop)]

use crate::surface::icosphere;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::cut::{CutMesh, CutOptions, cut};
use foehncraft_mesh::format::stl::Triangle;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::surface::{Surface, SurfaceBody};

fn box_mesh(kind: Kind, n: usize) -> foehncraft_mesh::Mesh {
    let mut s = BoxSpec::new(kind, 3, [n, n, n]);
    s.lower = [-1.0; 3];
    s.upper = [1.0; 3];
    s.generate().unwrap()
}

fn names() -> Vec<String> {
    vec!["body".to_string()]
}

/// The box of half-sides `h` about `c`, turned by `ax` about x, then by
/// `az` about z: sharp edges and corners at no angle to the mesh.
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

/// The errors of a cut against the exact moments, over the body volume
/// (volumes), the exact wall area (area) and the face area (apertures).
#[derive(Debug)]
struct Errors {
    volume: f64,
    area: f64,
    /// Σ |fluid volume − exact| over the groups of input cells joined by
    /// merged cells.
    cells: f64,
    /// Mean |fluid area − exact| over the cut faces.
    apertures: f64,
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

fn errors(m: &foehncraft_mesh::Mesh, s: &Surface, c: &CutMesh) -> Errors {
    let g = Geometry::compute(m);
    let mo = s.cut_moments(m, &g);
    let gc = Geometry::compute(&c.mesh);
    let mag = |a: [f64; 3]| (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    let v_body: f64 = mo.cell_inside_volume.iter().sum();
    let v_ex = g.cell_volume.iter().sum::<f64>() - v_body;
    let a_ex: f64 = mo.walls.iter().map(|w| mag(w.area)).sum();
    // groups of input cells joined by the merged cells
    let mut uf = Uf((0..m.n_cells).collect());
    for cov in &c.covers {
        for w in cov.windows(2) {
            let (a, b) = (uf.find(w[0] as usize), uf.find(w[1] as usize));
            uf.0[a.max(b)] = a.min(b);
        }
    }
    let mut vcut = vec![0.0; m.n_cells];
    let mut vex = vec![0.0; m.n_cells];
    for (k, cov) in c.covers.iter().enumerate() {
        let r = uf.find(cov[0] as usize);
        vcut[r] += gc.cell_volume[k];
    }
    for i in 0..m.n_cells {
        let r = uf.find(i);
        vex[r] += g.cell_volume[i] - mo.cell_inside_volume[i];
    }
    let cells: f64 = (0..m.n_cells).map(|i| (vcut[i] - vex[i]).abs()).sum();
    let (mut sum, mut count) = (0.0, 0usize);
    for f in 0..m.n_faces() {
        let af = mag(g.face_area[f]);
        let ex = af - mo.face_inside_area[f];
        if mo.face_inside_area[f] > 0.0 && ex > 0.0 {
            sum += (c.face_fluid_area[f] - ex).abs() / af;
            count += 1;
        }
    }
    Errors {
        volume: (c.report.fluid_volume - v_ex).abs() / v_body,
        area: (c.report.wall_area[0] - a_ex).abs() / a_ex,
        cells: cells / v_body,
        apertures: sum / count.max(1) as f64,
    }
}

fn opts(wall_point: bool, face_point: bool) -> CutOptions {
    CutOptions {
        min_fraction: 0.3,
        wall_point,
        face_point,
    }
}

fn usable(c: &CutMesh) {
    assert!(
        c.report.max_closure < 1e-12,
        "closure {}",
        c.report.max_closure
    );
    assert_eq!(c.report.bad_faces, 0, "{:?}", c.report);
    c.mesh.validate().unwrap();
}

/// The polyhedral sphere on every cell type: the cut converges at second
/// order to the exact moments of the same triangles, and the face points
/// make it more accurate.
#[test]
fn stl_sphere_cut_matches_the_exact_moments() {
    let s = Surface::new(&icosphere([0.031, -0.017, 0.012], 0.55, 4));
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly] {
        let mut prev: Option<f64> = None;
        for n in [8, 16] {
            let m = box_mesh(kind, n);
            let c = cut(&m, &s.body(), &names(), &opts(true, true)).unwrap();
            usable(&c);
            let e = errors(&m, &s, &c);
            if let Some(p) = prev {
                let order = (p / e.volume).log2();
                assert!(order > 1.6, "{kind:?}: volume order {order}, {e:?}");
            }
            prev = Some(e.volume);
            if n == 16 {
                assert!(e.volume < 1e-2, "{kind:?}: {e:?}");
                assert!(e.area < 1e-2, "{kind:?}: {e:?}");
                // cell by cell: what one cell lacks shows up in no other
                assert!(e.cells < 5e-2, "{kind:?}: {e:?}");
                if kind != Kind::Poly {
                    // the faces of the poly mesh are not planar: the fan of
                    // a clipped face is not that of the input face
                    assert!(e.apertures < 6e-3, "{kind:?}: {e:?}");
                }
                let plain = cut(&m, &s.body(), &names(), &opts(true, false)).unwrap();
                let ep = errors(&m, &s, &plain);
                assert!(e.volume < 0.8 * ep.volume, "{kind:?}: {e:?} vs {ep:?}");
                assert!(e.apertures <= ep.apertures, "{kind:?}: {e:?} vs {ep:?}");
            }
        }
    }
}

/// A box turned against the mesh: the dual-contouring wall point puts the
/// sharp edges and corners into the cut cells, the wall area converges at
/// second order instead of first.
#[test]
fn stl_box_edges_and_corners_are_kept() {
    let tris = rotated_box([0.021, -0.013, 0.017], [0.55, 0.4, 0.3], 0.35, 0.6);
    let s = Surface::new(&tris);
    let mut prev: Option<f64> = None;
    for n in [16, 32] {
        let m = box_mesh(Kind::Hex, n);
        let sharp = cut(&m, &s.body(), &names(), &opts(true, true)).unwrap();
        usable(&sharp);
        let smooth = SurfaceBody {
            surface: &s,
            exact: true,
            features: false,
        };
        let round = cut(&m, &smooth, &names(), &opts(true, false)).unwrap();
        usable(&round);
        let (e, er) = (errors(&m, &s, &sharp), errors(&m, &s, &round));
        assert!(e.area < 0.5 * er.area, "n = {n}: {e:?} vs {er:?}");
        assert!(e.volume < 0.5 * er.volume, "n = {n}: {e:?} vs {er:?}");
        if let Some(p) = prev {
            let order = (p / e.area).log2();
            assert!(order > 1.5, "area order {order}: {e:?}");
        }
        prev = Some(e.area);
    }
}

/// Where a sharp-edge point or a face point gives a face the
/// discretisation refuses, the cells next to it are cut with fewer
/// points: every face is usable and the rest keeps the finer points.
#[test]
fn fallback_keeps_every_face_usable() {
    let tris = rotated_box([0.021, -0.013, 0.017], [0.55, 0.4, 0.3], 0.35, 0.6);
    let s = Surface::new(&tris);
    for (surface, kind) in [(&s, Kind::Hex), (&s, Kind::Tet), (&s, Kind::Poly)] {
        let m = box_mesh(kind, 16);
        let c = cut(&m, &surface.body(), &names(), &opts(true, true)).unwrap();
        usable(&c);
        assert!(
            c.report.simplified > 0 && c.report.rounds > 1,
            "{:?}",
            c.report
        );
        assert!(
            c.report.simplified < c.report.cells_cut / 4,
            "{:?}",
            c.report
        );
        let whole = cut(&m, &surface.body(), &names(), &opts(false, false)).unwrap();
        usable(&whole);
        let (e, ew) = (errors(&m, surface, &c), errors(&m, surface, &whole));
        assert!(e.volume < ew.volume, "{kind:?}: {e:?} vs {ew:?}");
    }
}

#[test]
fn stl_cut_is_independent_of_threads() {
    let s = Surface::new(&rotated_box(
        [0.021, -0.013, 0.017],
        [0.55, 0.4, 0.3],
        0.35,
        0.6,
    ));
    let m = box_mesh(Kind::Tet, 10);
    let run = |t: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build()
            .unwrap()
            .install(|| cut(&m, &s.body(), &names(), &opts(true, true)).unwrap())
    };
    let (a, b) = (run(1), run(4));
    assert_eq!(a.mesh.points, b.mesh.points);
    assert_eq!(a.mesh.face_nodes, b.mesh.face_nodes);
    assert_eq!(a.mesh.owner, b.mesh.owner);
    assert_eq!(a.report, b.report);
}
