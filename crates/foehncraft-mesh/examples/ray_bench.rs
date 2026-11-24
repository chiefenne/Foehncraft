//! Ray casting (D108 part 2): the exact segment–triangle test from
//! orientation predicates against the floating-point Möller–Trumbore test
//! (1997), the common choice of renderers.
//!
//! ```text
//! cargo run --release -p foehncraft-mesh --example ray_bench
//! ```
//!
//! For the closed wall surface of generated meshes, turned about an
//! oblique axis so that no wall lies in a coordinate plane:
//!
//! * **Leaks.** Rays from interior points aimed exactly at every surface
//!   corner (mesh nodes and the centres of polygon fans) and at random
//!   points on every surface edge. The domain is convex, so the ray leaves
//!   it exactly there and only the triangles at that corner or edge can
//!   catch it; a leak is a ray that none of them catches. Counted for the
//!   exact test, for Möller–Trumbore without a tolerance, and with a
//!   tolerance of 1e-9 on the barycentric weights.
//! * **Cost.** Nanoseconds per ray–triangle test on random pairs (mostly
//!   misses, as in a traversal), and microseconds per first-hit query on
//!   the hierarchy for random rays from inside.
//!
//! Prints markdown tables.

use foehncraft_core::rng;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::query::{Bvh, End, SurfaceFaces, segment_triangle};
use foehncraft_mesh::{Geometry, Mesh};
use std::collections::HashMap;
use std::hint::black_box;
use std::time::Instant;

type V = [f64; 3];

fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Möller & Trumbore 1997: the ray o + t d against the triangle, with a
/// tolerance `eps` on the barycentric weights (0: none).
fn moller_trumbore(o: V, d: V, tri: &[V; 3], eps: f64) -> Option<f64> {
    let e1 = sub(tri[1], tri[0]);
    let e2 = sub(tri[2], tri[0]);
    let pv = cross(d, e2);
    let det = dot(e1, pv);
    if det == 0.0 {
        return None;
    }
    let inv = 1.0 / det;
    let tv = sub(o, tri[0]);
    let u = dot(tv, pv) * inv;
    if u < -eps || u > 1.0 + eps {
        return None;
    }
    let qv = cross(tv, e1);
    let v = dot(d, qv) * inv;
    if v < -eps || u + v > 1.0 + eps {
        return None;
    }
    let t = dot(e2, qv) * inv;
    (t > 0.0).then_some(t)
}

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

fn key(p: V) -> [u64; 3] {
    p.map(f64::to_bits)
}

fn main() {
    let make = |kind: Kind, n: [usize; 3]| {
        let mut s = BoxSpec::new(kind, 3, n);
        s.perturb = 0.2;
        let mut m = s.generate().expect("mesh");
        rotate(&mut m);
        m
    };
    let meshes = [
        ("hex 64³", make(Kind::Hex, [64, 64, 64])),
        ("tet 24³", make(Kind::Tet, [24, 24, 24])),
        ("poly 20³", make(Kind::Poly, [20, 20, 20])),
    ];
    println!(
        "| Mesh | Triangles | Rays | Leaks exact | Leaks Möller–Trumbore | Leaks Möller–Trumbore 1e-9 | Exact test (ns) | Möller–Trumbore (ns) | First hit on the hierarchy (µs/ray) | Brute force (µs/ray) |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|");
    for (name, m) in &meshes {
        let g = Geometry::compute(m);
        let all: Vec<usize> = (0..m.boundaries.len()).collect();
        let bvh = Bvh::build(SurfaceFaces::of_boundaries(m, &g, &all));
        let tris = &bvh.faces().tris;
        // corners and edges with their triangles
        let mut at_corner: HashMap<[u64; 3], Vec<usize>> = HashMap::new();
        let mut at_edge: HashMap<[[u64; 3]; 2], Vec<usize>> = HashMap::new();
        for (t, tri) in tris.iter().enumerate() {
            for k in 0..3 {
                at_corner.entry(key(tri[k])).or_default().push(t);
                let (a, b) = (key(tri[k]), key(tri[(k + 1) % 3]));
                at_edge
                    .entry(if a < b { [a, b] } else { [b, a] })
                    .or_default()
                    .push(t);
            }
        }
        let mut targets: Vec<(V, Vec<usize>)> = Vec::new();
        let mut corners: Vec<_> = at_corner.into_iter().collect();
        corners.sort();
        for (k, ts) in corners {
            targets.push((k.map(f64::from_bits), ts));
        }
        let mut edges: Vec<_> = at_edge.into_iter().collect();
        edges.sort();
        for (i, ([a, b], ts)) in edges.into_iter().enumerate() {
            let (a, b) = (a.map(f64::from_bits), b.map(f64::from_bits));
            for j in 0..2u64 {
                let s = rng::uniform([5, i as u64, j]);
                targets.push((
                    std::array::from_fn(|k| a[k] + s * (b[k] - a[k])),
                    ts.clone(),
                ));
            }
        }
        let origins: Vec<V> = (0..4)
            .map(|i| g.cell_centroid[i * m.n_cells / 4 + m.n_cells / 8])
            .collect();
        let (mut rays, mut leak_exact, mut leak_mt, mut leak_mt_eps) = (0usize, 0, 0, 0);
        for &o in &origins {
            for (x, ts) in &targets {
                let d = sub(*x, o);
                let q: V = std::array::from_fn(|k| o[k] + 2.0 * d[k]);
                rays += 1;
                if !ts
                    .iter()
                    .any(|&t| segment_triangle(o, q, &tris[t], End::Closed).is_some())
                {
                    leak_exact += 1;
                }
                if !ts
                    .iter()
                    .any(|&t| moller_trumbore(o, d, &tris[t], 0.0).is_some())
                {
                    leak_mt += 1;
                }
                if !ts
                    .iter()
                    .any(|&t| moller_trumbore(o, d, &tris[t], 1e-9).is_some())
                {
                    leak_mt_eps += 1;
                }
            }
        }
        // cost per test on random ray–triangle pairs
        let (lo, hi) = m.bounding_box();
        let pt = |s: u64, i: u64| -> V {
            std::array::from_fn(|k| lo[k] + (hi[k] - lo[k]) * rng::uniform([s, i, k as u64]))
        };
        let n = 2_000_000u64;
        let pairs: Vec<(V, V, usize)> = (0..n)
            .map(|i| {
                let t = (rng::uniform([9, i, 0]) * tris.len() as f64) as usize % tris.len();
                (pt(11, i), pt(12, i), t)
            })
            .collect();
        let start = Instant::now();
        let mut c = 0usize;
        for &(p, q, t) in &pairs {
            c += usize::from(segment_triangle(p, q, black_box(&tris[t]), End::Closed).is_some());
        }
        let exact_ns = start.elapsed().as_secs_f64() * 1e9 / n as f64;
        let start = Instant::now();
        let mut c2 = 0usize;
        for &(p, q, t) in &pairs {
            let d = sub(q, p);
            c2 += usize::from(
                moller_trumbore(p, d, black_box(&tris[t]), 0.0).is_some_and(|s| s <= 1.0),
            );
        }
        let mt_ns = start.elapsed().as_secs_f64() * 1e9 / n as f64;
        black_box((c, c2));
        // first hit for random rays from inside
        let nr = 100_000u64;
        let rays_in: Vec<(V, V)> = (0..nr)
            .map(|i| {
                let o = g.cell_centroid
                    [(rng::uniform([13, i, 0]) * m.n_cells as f64) as usize % m.n_cells];
                let d: V = std::array::from_fn(|k| rng::uniform([14, i, k as u64]) - 0.5);
                (o, d)
            })
            .collect();
        let start = Instant::now();
        let mut missed = 0usize;
        for &(o, d) in &rays_in {
            missed += usize::from(bvh.ray_first_hit(o, d).is_none());
        }
        let tree_us = start.elapsed().as_secs_f64() * 1e6 / nr as f64;
        assert_eq!(
            missed, 0,
            "{name}: a ray from inside missed the closed surface"
        );
        // brute force on a few rays, the same answers
        let nb = 200usize;
        let start = Instant::now();
        let far: Vec<(V, V)> = rays_in[..nb]
            .iter()
            .map(|&(o, d)| (o, std::array::from_fn(|k| o[k] + 8.0 * d[k])))
            .collect();
        let b: Vec<_> = far
            .iter()
            .map(|&(o, q)| bvh.first_hit_brute_force(o, q))
            .collect();
        let brute_us = start.elapsed().as_secs_f64() * 1e6 / nb as f64;
        for (&(o, q), b) in far.iter().zip(&b) {
            assert_eq!(bvh.first_hit(o, q), *b, "{name}: brute force disagrees");
        }
        println!(
            "| {name} | {} | {rays} | {leak_exact} | {leak_mt} | {leak_mt_eps} | {exact_ns:.1} | {mt_ns:.1} | {tree_us:.2} | {brute_us:.0} |",
            tris.len()
        );
    }
}
