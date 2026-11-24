//! Ray casting on the surface hierarchy (D108 part 2): watertight on closed
//! surfaces of every mesh type, the same answer as testing every triangle,
//! visibility between surface points, no self-hits from a wall.

use foehncraft_core::rng;
use foehncraft_mesh::generate::{BoxSpec, Kind, RemoveBox};
use foehncraft_mesh::query::{Bvh, SurfaceFaces};
use foehncraft_mesh::{Geometry, Mesh};

fn point(seed: u64, i: u64, lo: [f64; 3], hi: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|d| lo[d] + (hi[d] - lo[d]) * rng::uniform([seed, i, d as u64]))
}

/// Turns every mesh point about an oblique axis, so that no wall is
/// aligned with a coordinate plane and the corners are not representable
/// sums of the cell sizes.
fn rotate(m: &mut Mesh) {
    let (s, c) = 0.7f64.sin_cos();
    let (s2, c2) = 0.4f64.sin_cos();
    for p in &mut m.points {
        let [x, y, z] = *p;
        let (x, y) = (c * x - s * y, s * x + c * y);
        let (y, z) = (c2 * y - s2 * z, s2 * y + c2 * z);
        *p = [x, y, z];
    }
}

fn surface(m: &Mesh, g: &Geometry) -> Bvh {
    let all: Vec<usize> = (0..m.boundaries.len()).collect();
    Bvh::build(SurfaceFaces::of_boundaries(m, g, &all))
}

#[test]
fn rays_from_inside_never_leak_through_a_closed_surface() {
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        let mut s = BoxSpec::new(kind, 3, [6, 5, 4]);
        s.perturb = 0.25;
        s.upper = [2.0, 1.0, 0.5];
        let mut m = s.generate().unwrap();
        rotate(&mut m);
        let g = Geometry::compute(&m);
        let bvh = surface(&m, &g);
        let tris = &bvh.faces().tris;
        // aim at every triangle corner: mesh nodes on the surface and the
        // centres of the fans of polygon faces, the hardest directions
        let mut targets: Vec<[f64; 3]> = tris.iter().flat_map(|t| t.iter().copied()).collect();
        targets.sort_by(|a, b| a.partial_cmp(b).unwrap());
        targets.dedup();
        let corners = targets.len();
        // and at edge midpoints (rounded: on the surface or a little off)
        targets.extend(
            tris.iter()
                .map(|t| std::array::from_fn(|k| 0.5 * (t[0][k] + t[1][k]))),
        );
        let name = kind.name();
        let mut hits = 0usize;
        for c in [0, m.n_cells / 3, m.n_cells / 2, m.n_cells - 1] {
            let o = g.cell_centroid[c];
            for (j, &x) in targets.iter().enumerate() {
                let d: [f64; 3] = std::array::from_fn(|k| x[k] - o[k]);
                let h = bvh
                    .ray_first_hit(o, d)
                    .unwrap_or_else(|| panic!("{name}: ray from cell {c} to {x:?} leaked"));
                // the domain is convex: the first hit is the aimed point
                let e: f64 = (0..3).map(|k| (h.point[k] - x[k]).powi(2)).sum();
                assert!(e.sqrt() < 1e-9, "{name}: hit {:?} for {x:?}", h.point);
                assert!((h.t - 1.0).abs() < 1e-9, "{name}: t {}", h.t);
                assert!(h.entering, "{name}: leaves the domain into the wall");
                // the segment to an aimed corner ends on the surface
                if j < corners {
                    let hs = bvh.first_hit(o, x).expect("closed end meets the wall");
                    assert!((hs.t - 1.0).abs() < 1e-9);
                    assert!(!bvh.occluded(o, x), "{name}: between {o:?} and {x:?}");
                }
                hits += 1;
            }
        }
        assert!(hits > 100, "{name}: {hits}");
    }
}

#[test]
fn the_hierarchy_agrees_with_testing_every_triangle() {
    for kind in [Kind::Hex, Kind::Tet] {
        // the removed block on cell faces (unperturbed), then all turned
        let mut s = BoxSpec::new(kind, 3, [16, 8, 8]);
        s.remove.push(RemoveBox {
            lower: [0.5, 0.25, 0.25],
            upper: [1.25, 0.75, 0.75],
            name: "block".into(),
        });
        s.upper = [2.0, 1.0, 1.0];
        let mut m = s.generate().unwrap();
        rotate(&mut m);
        let g = Geometry::compute(&m);
        let bvh = surface(&m, &g);
        let (lo, hi) = m.bounding_box();
        let lo = lo.map(|v| v - 0.3);
        let hi = hi.map(|v| v + 0.3);
        let mut found = 0;
        for i in 0..3000u64 {
            let p = point(31, 2 * i, lo, hi);
            let q = point(31, 2 * i + 1, lo, hi);
            let a = bvh.first_hit(p, q);
            let b = bvh.first_hit_brute_force(p, q);
            assert_eq!(a, b, "{}: {p:?} {q:?}", kind.name());
            assert_eq!(bvh.occluded(p, q), a.is_some_and(|h| h.t < 1.0));
            found += usize::from(a.is_some());
        }
        assert!(found > 1000, "{found}");
    }
}

#[test]
fn visibility_between_wall_points_and_no_self_hit() {
    let mut s = BoxSpec::new(Kind::Hex, 3, [20, 10, 10]);
    s.upper = [2.0, 1.0, 1.0];
    s.remove.push(RemoveBox {
        lower: [0.8, 0.3, 0.3],
        upper: [1.2, 0.7, 0.7],
        name: "block".into(),
    });
    let m = s.generate().unwrap();
    let g = Geometry::compute(&m);
    let bvh = surface(&m, &g);
    let block = m.boundaries.iter().position(|b| b.name == "block").unwrap();
    let on_block = |h: &foehncraft_mesh::query::Hit| {
        m.boundaries[block]
            .range()
            .contains(&bvh.faces().mesh_face[h.face])
    };
    // wall to wall past the block: visible; through it: hidden
    assert!(!bvh.occluded([0.0, 0.5, 0.1], [2.0, 0.5, 0.1]));
    assert!(bvh.occluded([0.0, 0.5, 0.5], [2.0, 0.5, 0.5]));
    // started on the wall x = 0 (on its plane, on a face edge, at a
    // corner): never the start wall, the far wall at t = 1
    for p in [[0.0, 0.5, 0.15], [0.0, 0.5, 0.2], [0.0, 0.0, 0.0]] {
        let q = [2.0, p[1], p[2]];
        let h = bvh.first_hit(p, q).unwrap();
        assert_eq!(h.t, 1.0, "{p:?}");
        assert!((h.point[0] - 2.0).abs() < 1e-12);
        assert!(h.entering);
    }
    // the block is met at x = 0.8 from the left, entering it
    let h = bvh.first_hit([0.1, 0.5, 0.5], [1.9, 0.5, 0.5]).unwrap();
    assert!(on_block(&h));
    assert!((h.point[0] - 0.8).abs() < 1e-12);
    assert!(h.entering);
    assert!((h.normal[0] + 1.0).abs() < 1e-12);
    // from inside the block (a segment that leaves a wall) it is not
    // entering
    let h = bvh.first_hit([1.0, 0.5, 0.5], [1.9, 0.5, 0.5]).unwrap();
    assert!(on_block(&h));
    assert!(!h.entering);
    // through an edge and a corner of the block: the smallest global id of
    // the faces there wins, the same in both directions of traversal order
    for x in [[0.8, 0.5, 0.3], [0.8, 0.3, 0.3]] {
        let p = [0.1, x[1] + 0.05, x[2] + 0.05];
        let q: [f64; 3] = std::array::from_fn(|k| p[k] + 2.0 * (x[k] - p[k]));
        let a = bvh.first_hit(p, q).unwrap();
        let b = bvh.first_hit_brute_force(p, q).unwrap();
        assert_eq!(a, b);
        assert!(on_block(&a));
    }
}
