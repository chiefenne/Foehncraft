//! Near-body refinement with hanging-node faces (L13 stage 2f, D57):
//! valid closed cells, volume and boundary area conserved, 2:1 balance
//! across faces, edges and corners, planar meshes refined in the plane,
//! other cell types kept, the same result on every run.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::Mesh;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::quality::QualityReport;
use foehncraft_mesh::refine::{Candidate, RefineOptions, Refined, refine};

fn box_mesh(kind: Kind, dim: u8, n: usize, perturb: f64) -> Mesh {
    let mut s = BoxSpec::new(kind, dim, [n, n, n]);
    s.lower = [-1.0, -1.0, if dim == 2 { 0.0 } else { -1.0 }];
    s.upper = [1.0, 1.0, if dim == 2 { 0.1 } else { 1.0 }];
    s.perturb = perturb;
    s.generate().unwrap()
}

/// A band around the sphere of radius 0.55.
fn near_sphere(c: &Candidate) -> bool {
    let x = c.centre;
    let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    (r - 0.55).abs() - c.radius < 2.0 * c.size
}

fn near_circle(c: &Candidate) -> bool {
    let x = c.centre;
    ((x[0] * x[0] + x[1] * x[1]).sqrt() - 0.55).abs() - c.radius < 2.0 * c.size
}

/// The boundary area vector and area of every boundary.
fn boundary_areas(m: &Mesh, g: &Geometry) -> Vec<(String, [f64; 3], f64)> {
    m.boundaries
        .iter()
        .map(|b| {
            let mut s = [0.0; 3];
            let mut a = 0.0;
            for f in b.range() {
                let v = g.face_area[f];
                s = [s[0] + v[0], s[1] + v[1], s[2] + v[2]];
                a += (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            }
            (b.name.clone(), s, a)
        })
        .collect()
}

/// The checks every refined mesh must pass.
fn check(input: &Mesh, r: &Refined) {
    let m = &r.mesh;
    m.validate().unwrap();
    let g0 = Geometry::compute(input);
    let g = Geometry::compute(m);
    let (q, _) = QualityReport::compute(m, &g);
    assert_eq!(q.negative_volumes, 0);
    assert_eq!(q.inverted_faces, 0);
    assert!(q.closure.max < 1e-13, "closure {}", q.closure.max);
    // watertight: the volume of the cells is the volume the boundary encloses
    assert!(
        (q.total_volume - q.boundary_volume).abs() < 1e-12 * q.total_volume,
        "{} {}",
        q.total_volume,
        q.boundary_volume
    );
    // the boundary faces of each input face add up to its area vector
    for ((n0, s0, a0), (n1, s1, a1)) in boundary_areas(input, &g0).iter().zip(boundary_areas(m, &g))
    {
        assert_eq!(n0, &n1);
        for d in 0..3 {
            assert!(
                (s0[d] - s1[d]).abs() < 1e-12 * a0.max(1e-300),
                "{n0}: {s0:?} {s1:?}"
            );
        }
        assert!((a0 - a1).abs() < 1e-12 * a0, "{n0}: {a0} {a1}");
    }
    // 2:1 across faces
    for f in 0..m.n_interior_faces() {
        let (a, b) = (m.owner[f] as usize, m.neighbour[f] as usize);
        assert!(
            r.level[a].abs_diff(r.level[b]) <= 1,
            "face {f}: {} {}",
            r.level[a],
            r.level[b]
        );
    }
    // 2:1 across edges and corners: around every node
    let mut lo = vec![u8::MAX; m.n_points()];
    let mut hi = vec![0u8; m.n_points()];
    for (c, &l) in r.level.iter().enumerate() {
        for n in m.cell_nodes(c) {
            lo[n as usize] = lo[n as usize].min(l);
            hi[n as usize] = hi[n as usize].max(l);
        }
    }
    for n in 0..m.n_points() {
        assert!(
            hi[n] - lo[n] <= 1,
            "node {n}: levels {} to {}",
            lo[n],
            hi[n]
        );
    }
    // the input nodes keep their place
    assert_eq!(&m.points[..input.n_points()], &input.points[..]);
    assert_eq!(r.parent.len(), m.n_cells);
    assert_eq!(r.report.cells_out, m.n_cells);
    let _ = g0;
}

/// The volume of the cells of each input cell, against its volume.
fn parent_volumes(input: &Mesh, r: &Refined) -> f64 {
    let g0 = Geometry::compute(input);
    let g = Geometry::compute(&r.mesh);
    let mut v = vec![0.0; input.n_cells];
    for c in 0..r.mesh.n_cells {
        v[r.parent[c] as usize] += g.cell_volume[c];
    }
    (0..input.n_cells)
        .map(|c| (v[c] / g0.cell_volume[c] - 1.0).abs())
        .fold(0.0, f64::max)
}

#[test]
fn a_band_around_a_sphere_keeps_volume_and_area() {
    let input = box_mesh(Kind::Hex, 3, 8, 0.0);
    let opt = RefineOptions {
        levels: 2,
        keep_direction: None,
    };
    let r = refine(&input, &opt, near_sphere).unwrap();
    check(&input, &r);
    let rep = &r.report;
    eprintln!("sphere band: {rep:?}");
    assert_eq!(rep.per_level.len(), 3);
    assert!(rep.per_level[2] > 0 && rep.hanging_faces > 0 && rep.faces_with_hanging_nodes > 0);
    // box cells: every refined cell keeps its volume to round-off
    let e = parent_volumes(&input, &r);
    assert!(e < 1e-13, "volume of the parts {e}");
    // the same result on every run (no hash order in the output)
    let again = refine(&input, &opt, near_sphere).unwrap();
    assert_eq!(again.mesh.points, r.mesh.points);
    assert_eq!(again.mesh.face_nodes, r.mesh.face_nodes);
    assert_eq!(again.mesh.owner, r.mesh.owner);
}

#[test]
fn perturbed_hexahedra_are_refined_into_closed_cells() {
    let input = box_mesh(Kind::Hex, 3, 8, 0.2);
    let opt = RefineOptions {
        levels: 3,
        keep_direction: None,
    };
    let r = refine(&input, &opt, near_sphere).unwrap();
    check(&input, &r);
    // the faces of a warped hexahedron are not planes, yet its parts fill
    // it (measured 2e-15): the trilinear parts of a face pass through the
    // centre the fan of the input face uses
    let e = parent_volumes(&input, &r);
    eprintln!(
        "perturbed: {:?}, volume of the parts within {e:.2e}",
        r.report
    );
    assert!(e < 1e-12, "volume of the parts {e}");
    assert_eq!(r.report.per_level.len(), 4);
}

#[test]
fn a_planar_mesh_is_refined_in_its_plane() {
    let input = box_mesh(Kind::Hex, 2, 16, 0.15);
    let opt = RefineOptions {
        levels: 3,
        keep_direction: Some([0.0, 0.0, 1.0]),
    };
    let r = refine(&input, &opt, near_circle).unwrap();
    check(&input, &r);
    let m = &r.mesh;
    assert!(m.points.iter().all(|p| p[2] == 0.0 || p[2] == 0.1));
    // every cell still spans the layer: two faces on the z planes
    let zmin = m.boundary("zmin").unwrap();
    let zmax = m.boundary("zmax").unwrap();
    assert_eq!(zmin.len, m.n_cells);
    assert_eq!(zmax.len, m.n_cells);
    let e = parent_volumes(&input, &r);
    assert!(e < 1e-13, "volume of the parts {e}");
    // a coarse cell meets at most two finer cells on each of its four sides
    let mut count = std::collections::HashMap::new();
    for f in 0..m.n_interior_faces() {
        let (a, b) = (m.owner[f] as usize, m.neighbour[f] as usize);
        if r.level[a] != r.level[b] {
            let coarse = if r.level[a] < r.level[b] { a } else { b };
            *count.entry(coarse).or_insert(0usize) += 1;
        }
    }
    eprintln!("planar: {:?}", r.report);
    assert!(!count.is_empty());
    assert!(
        count.values().all(|&n| n <= 8),
        "{:?}",
        count.values().max()
    );
}

#[test]
fn one_point_refined_three_times_is_balanced() {
    let input = box_mesh(Kind::Hex, 3, 6, 0.0);
    let opt = RefineOptions {
        levels: 3,
        keep_direction: None,
    };
    // only the cells around a point just inside an input cell, next to its
    // face x = 0: the cell across that face must follow
    let p = [0.001, 0.17, 0.17];
    let r = refine(&input, &opt, |c| {
        let d = [c.centre[0] - p[0], c.centre[1] - p[1], c.centre[2] - p[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() < c.radius
    })
    .unwrap();
    check(&input, &r);
    eprintln!("point: {:?}", r.report);
    assert!(r.report.balance > 0);
    assert_eq!(r.report.per_level.len(), 4);
}

#[test]
fn other_cells_are_kept_and_bad_settings_refused() {
    let input = box_mesh(Kind::Tet, 3, 3, 0.0);
    let opt = RefineOptions {
        levels: 1,
        keep_direction: None,
    };
    let r = refine(&input, &opt, |_| true).unwrap();
    assert_eq!(r.mesh.n_cells, input.n_cells);
    assert_eq!(r.report.not_hexahedra, input.n_cells);
    assert_eq!(r.report.refined, 0);
    let hex = box_mesh(Kind::Hex, 3, 2, 0.0);
    for levels in [0, 7] {
        assert!(
            refine(
                &hex,
                &RefineOptions {
                    levels,
                    keep_direction: None
                },
                |_| true
            )
            .is_err()
        );
    }
    // a 3D mesh has no direction across a layer
    assert!(
        refine(
            &hex,
            &RefineOptions {
                levels: 1,
                keep_direction: Some([1.0, 1.0, 1.0])
            },
            |_| true
        )
        .is_err()
    );
}
