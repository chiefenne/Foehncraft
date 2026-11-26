//! L13 stage 0 (D225): the surface check, winding numbers, signed distance
//! and crossings of `foehncraft_mesh::surface`.

use foehncraft_mesh::format::stl::Triangle;
use foehncraft_mesh::surface::{Surface, check};
use std::f64::consts::PI;

/// Icosphere of radius `r` about `c`, `levels` subdivisions, outward.
pub fn icosphere(c: [f64; 3], r: f64, levels: usize) -> Vec<Triangle> {
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
                let key = (a.min(b), a.max(b));
                m[k] = *mid.entry(key).or_insert_with(|| {
                    let p = [0, 1, 2].map(|d| 0.5 * (v[a][d] + v[b][d]));
                    v.push(unit(p));
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

/// Torus about the z axis, radii `big` and `small`, `n` by `m` quads.
pub fn torus(big: f64, small: f64, n: usize, m: usize) -> Vec<Triangle> {
    let p = |i: usize, j: usize| {
        let u = 2.0 * PI * (i % n) as f64 / n as f64;
        let v = 2.0 * PI * (j % m) as f64 / m as f64;
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

/// Unit cube [0, 1]³, two triangles per side, outward.
pub fn cube() -> Vec<Triangle> {
    let v = |i: usize| [(i & 1) as f64, (i >> 1 & 1) as f64, (i >> 2 & 1) as f64];
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

/// Deterministic points in a box (a simple linear congruential sequence).
pub fn points(n: usize, lo: f64, hi: f64) -> Vec<[f64; 3]> {
    let mut s: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut next = || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        lo + (hi - lo) * ((s >> 11) as f64 / (1u64 << 53) as f64)
    };
    (0..n).map(|_| [next(), next(), next()]).collect()
}

#[test]
fn clean_sphere_passes_the_check() {
    let s = Surface::new(&icosphere([0.0; 3], 1.0, 3));
    let r = check(&s);
    assert!(r.clean(), "{:?}", r.findings);
    assert!(r.findings.is_empty(), "{:?}", r.findings);
    assert_eq!(r.components, 1);
    assert_eq!(r.triangles, 1280);
    assert_eq!(r.vertices, 642);
    // inscribed polyhedron: slightly less than the sphere
    let exact = 4.0 / 3.0 * PI;
    assert!(r.volume < exact && r.volume > 0.98 * exact, "{}", r.volume);
    let t = Surface::new(&torus(1.0, 0.3, 48, 24));
    let rt = check(&t);
    assert!(rt.clean(), "{:?}", rt.findings);
}

#[test]
fn defects_are_reported_not_repaired() {
    let tris = icosphere([0.0; 3], 1.0, 2);
    // hole: three triangles removed
    let holed: Vec<Triangle> = tris[3..].to_vec();
    let r = check(&Surface::new(&holed));
    assert!(!r.closed);
    assert!(r.boundary_edges > 0);
    assert!(r.findings.iter().any(|f| f.code == "IMM002"));
    // inward normals
    let flipped: Vec<Triangle> = tris.iter().map(|t| [t[0], t[2], t[1]]).collect();
    let r = check(&Surface::new(&flipped));
    assert!(r.closed && r.volume < 0.0);
    assert!(r.findings.iter().any(|f| f.code == "IMM008"));
    // one flipped triangle
    let mut one = tris.clone();
    one[5] = [one[5][0], one[5][2], one[5][1]];
    let r = check(&Surface::new(&one));
    assert_eq!(r.flipped_edges, 3);
    assert!(r.findings.iter().any(|f| f.code == "IMM004"));
    // duplicate and degenerate triangles
    let mut dup = tris.clone();
    dup.push(tris[7]);
    dup.push([[0.0; 3], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]);
    let r = check(&Surface::new(&dup));
    assert_eq!(r.n_duplicates, 1);
    assert_eq!(r.n_degenerate, 1);
    assert!(r.findings.iter().any(|f| f.code == "IMM006"));
    assert!(r.findings.iter().any(|f| f.code == "IMM007"));
    assert!(r.findings.iter().any(|f| f.code == "IMM003"));
    // two overlapping spheres intersect
    let mut two = icosphere([0.0; 3], 1.0, 2);
    two.extend(icosphere([1.2, 0.1, 0.05], 1.0, 2));
    let r = check(&Surface::new(&two));
    assert!(r.n_intersecting > 0);
    assert_eq!(r.components, 2);
    assert!(r.findings.iter().any(|f| f.code == "IMM005"));
    // millimetres
    let mm: Vec<Triangle> = tris
        .iter()
        .map(|t| t.map(|p| p.map(|x| 50.0 * x)))
        .collect();
    let r = check(&Surface::new(&mm));
    assert!(r.findings.iter().any(|f| f.code == "IMM009"));
}

#[test]
fn degenerate_contacts_are_exact() {
    let base: Triangle = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let pair = |b: Triangle| check(&Surface::new(&[base, b])).n_intersecting;
    // a corner exactly on the face
    assert_eq!(
        pair([[0.25, 0.25, 0.0], [0.3, 0.2, 1.0], [0.2, 0.3, 1.0]]),
        1
    );
    // a corner a hair above the face (1e-300 is exact in binary64)
    assert_eq!(
        pair([[0.25, 0.25, 1e-300], [0.3, 0.2, 1.0], [0.2, 0.3, 1.0]]),
        0
    );
    // coplanar overlap
    assert_eq!(pair([[0.2, 0.2, 0.0], [2.0, 0.2, 0.0], [0.2, 2.0, 0.0]]), 1);
    // coplanar, apart
    assert_eq!(pair([[1.0, 1.0, 0.0], [2.0, 1.0, 0.0], [1.0, 2.0, 0.0]]), 0);
    // an edge exactly through the edge of the other
    assert_eq!(
        pair([[0.5, -1.0, 0.0], [0.5, 1.0, 0.0], [0.5, 0.0, 1.0]]),
        1
    );
    // shared edge, folded onto each other
    assert_eq!(pair([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.1, 0.1, 0.0]]), 1);
    // shared edge, a proper neighbour
    assert_eq!(pair([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]), 0);
    // shared corner, the other triangle pierces this one
    assert_eq!(
        pair([[0.0, 0.0, 0.0], [0.5, 0.3, -1.0], [0.3, 0.5, 1.0]]),
        1
    );
    // shared corner only
    assert_eq!(
        pair([[0.0, 0.0, 0.0], [-1.0, 0.0, 1.0], [0.0, -1.0, 1.0]]),
        0
    );
}

#[test]
fn fast_winding_numbers_match_the_exact_sum() {
    let s = Surface::new(&icosphere([0.1, -0.2, 0.05], 1.0, 4));
    let q = points(2000, -2.0, 2.0);
    let mut worst: f64 = 0.0;
    for &x in &q {
        let fast = s.winding_number(x);
        let exact = s.winding_number_exact(x);
        worst = worst.max((fast - exact).abs());
        // infinite ratio: the hierarchy sums every triangle exactly
        let all = s.winding_number_with(x, f64::INFINITY);
        assert!((all - exact).abs() < 1e-12);
        let r = ((x[0] - 0.1).powi(2) + (x[1] + 0.2).powi(2) + (x[2] - 0.05).powi(2)).sqrt();
        if r < 0.98 {
            assert!((exact - 1.0).abs() < 1e-9, "{exact} at {x:?}");
        } else if r > 1.02 {
            assert!(exact.abs() < 1e-9, "{exact} at {x:?}");
        }
    }
    assert!(worst < 5e-3, "fast against exact: {worst}");
    // torus: the hole is outside
    let t = Surface::new(&torus(1.0, 0.3, 64, 32));
    assert!(t.winding_number([0.0; 3]).abs() < 5e-3);
    assert!((t.winding_number([1.0, 0.0, 0.0]) - 1.0).abs() < 5e-3);
    assert!(t.inside([0.0, -1.05, 0.1]));
    assert!(!t.inside([0.0, -1.35, 0.0]));
}

#[test]
fn winding_numbers_tolerate_dirty_surfaces() {
    let tris = icosphere([0.0; 3], 1.0, 3);
    // a hole of about 2 % of the area
    let holed: Vec<Triangle> = tris
        .iter()
        .copied()
        .filter(|t| t.iter().map(|p| p[2]).sum::<f64>() / 3.0 < 0.97)
        .collect();
    assert!(holed.len() < tris.len());
    let s = Surface::new(&holed);
    let w = s.winding_number([0.0; 3]);
    assert!(w > 0.9 && w < 1.0, "{w}");
    assert!(s.inside([0.0, 0.0, 0.5]));
    assert!(!s.inside([0.0, 0.0, 1.5]));
    // duplicated triangles and an inner overlapping shell: the half-way
    // test still separates inside from outside
    let mut over = tris.clone();
    over.extend(icosphere([0.3, 0.0, 0.0], 0.5, 2));
    let s = Surface::new(&over);
    assert!(s.inside([0.0, 0.0, 0.0]));
    assert!(s.inside([-0.7, 0.0, 0.0]));
    assert!(!s.inside([1.5, 0.0, 0.0]));
}

#[test]
fn signed_distance_converges_to_the_sphere() {
    let q = points(500, -1.6, 1.6);
    let mut last = f64::INFINITY;
    for level in 2..6 {
        let s = Surface::new(&icosphere([0.0; 3], 1.0, level));
        let mut worst: f64 = 0.0;
        for &x in &q {
            let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            let d = s.signed_distance(x);
            // exact closest point of the hierarchy against brute force
            let b = s.bvh().brute_force(x).unwrap().distance;
            assert_eq!(d.abs(), b);
            worst = worst.max((d - (r - 1.0)).abs());
        }
        // inscribed polyhedron: the error is the sagitta, h² / 8
        assert!(worst < 0.4 * last, "level {level}: {worst} after {last}");
        last = worst;
    }
    assert!(last < 2e-3, "{last}");
}

#[test]
fn narrow_band_and_threads() {
    let s = Surface::new(&icosphere([0.0; 3], 1.0, 4));
    let q = points(4000, -1.5, 1.5);
    let band = 0.1;
    let full = s.signed_distances(&q, None);
    let banded = s.signed_distances(&q, Some(band));
    for (f, b) in full.iter().zip(&banded) {
        assert_eq!(*b, f.clamp(-band, band));
    }
    let run = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| (s.signed_distances(&q, Some(band)), s.winding_numbers(&q)))
    };
    let one = run(1);
    for t in [2, 4] {
        let other = run(t);
        assert!(
            one.0
                .iter()
                .zip(&other.0)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
        assert!(
            one.1
                .iter()
                .zip(&other.1)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
}

#[test]
fn crossings_find_the_exact_roots() {
    let s = Surface::new(&icosphere([0.0; 3], 1.0, 3));
    let c = s.crossings([-2.0, 0.013, 0.021], [2.0, 0.013, 0.021]);
    assert_eq!(c.len(), 2);
    assert!(c[0].entering && !c[1].entering);
    assert!(c[0].t < c[1].t);
    // the cube: a segment through an edge of two triangles (the side
    // diagonal) counts once, and the root is exact
    let s = Surface::new(&cube());
    let c = s.crossings([0.5, 0.5, -1.0], [0.5, 0.5, 0.5]);
    assert_eq!(c.len(), 1, "{c:?}");
    assert!(c[0].entering);
    assert_eq!(c[0].point[2], 0.0);
    let f = s
        .first_crossing([0.25, 0.5, 2.0], [0.25, 0.5, 0.5])
        .unwrap();
    assert!(f.entering);
    assert_eq!(f.point, [0.25, 0.5, 1.0]);
    assert_eq!(s.signed_distance([0.5, 0.5, 0.25]), -0.25);
    assert_eq!(s.signed_distance([0.5, 0.5, 1.5]), 0.5);
}
