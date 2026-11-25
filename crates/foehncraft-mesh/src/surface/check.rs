//! The surface check (L13 stage 0, D225 item 1): what is wrong with a
//! triangle surface before it is used as an immersed body. A report, no
//! silent repair (owner decision Q2, 2026-11-22).
//!
//! Every decision is exact: corners are the same vertex only when their
//! coordinates are bitwise equal; a triangle is degenerate when its corners
//! are exactly collinear, and two triangles intersect when the exact
//! orientation predicates of D225 (Shewchuk 1997) say so. Self-intersection
//! candidates come from the surface hierarchy (box overlap), so the check
//! costs O(n log n).
//!
//! Two triangles that share an edge intersect only when they fold onto
//! each other (coplanar, the third corners on the same side); two that
//! share a corner intersect when they meet anywhere else. Touching counts
//! as intersecting (a corner on another triangle is a defect of a solid).

use super::Surface;
use crate::geometry::{cross, dot, sub};
use foehncraft_core::predicates::{orient2d, orient3d};
use rayon::prelude::*;
use schemars::JsonSchema;
use serde::Serialize;
use std::collections::HashMap;

/// Most items of a list kept in the report.
const LISTED: usize = 20;

/// One finding of the check.
#[derive(Debug, Clone, Serialize, PartialEq, JsonSchema)]
pub struct Finding {
    /// Diagnostic code (`IMM001` to `IMM009`).
    pub code: &'static str,
    /// What was found.
    pub message: String,
}

/// The report of [`check`].
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct SurfaceReport {
    /// Triangles.
    pub triangles: usize,
    /// Distinct corners.
    pub vertices: usize,
    /// Lower corner of the bounding box.
    pub lower: [f64; 3],
    /// Upper corner of the bounding box.
    pub upper: [f64; 3],
    /// Length of the box diagonal.
    pub diagonal: f64,
    /// Total area.
    pub area: f64,
    /// Enclosed volume, Σ a·(b×c)/6: positive for an outward closed
    /// surface (meaningful only when closed).
    pub volume: f64,
    /// Connected pieces (triangles joined by edges).
    pub components: usize,
    /// Every edge has exactly two triangles that run it in opposite
    /// directions: a closed, consistently oriented 2-manifold.
    pub closed: bool,
    /// Edges with only one triangle (holes, open borders).
    pub boundary_edges: usize,
    /// Edges with more than two triangles.
    pub non_manifold_edges: usize,
    /// Edges whose two triangles run them in the same direction.
    pub flipped_edges: usize,
    /// Triangles with collinear or repeated corners (first ones).
    pub degenerate: Vec<usize>,
    /// Number of degenerate triangles.
    pub n_degenerate: usize,
    /// Pairs of triangles with the same three corners (first ones).
    pub duplicates: Vec<[usize; 2]>,
    /// Number of duplicate pairs.
    pub n_duplicates: usize,
    /// Pairs of triangles that intersect (first ones).
    pub intersecting: Vec<[usize; 2]>,
    /// Number of intersecting pairs.
    pub n_intersecting: usize,
    /// Distinct corners closer than 1e-9 of the box diagonal (welding
    /// candidates; not welded).
    pub near_duplicate_vertices: usize,
    /// The findings, worst first.
    pub findings: Vec<Finding>,
}

impl SurfaceReport {
    /// The surface can be used as a solid body as it is: closed, oriented
    /// outward, no self-intersections, no degenerate or duplicate
    /// triangles.
    pub fn clean(&self) -> bool {
        self.closed
            && self.volume > 0.0
            && self.n_intersecting == 0
            && self.n_degenerate == 0
            && self.n_duplicates == 0
    }
}

/// Checks the surface.
pub fn check(s: &Surface) -> SurfaceReport {
    let n = s.len();
    let (lower, upper) = s.bounding_box();
    let diagonal = if n == 0 {
        0.0
    } else {
        let d = sub(upper, lower);
        dot(d, d).sqrt()
    };
    let mut area = 0.0;
    let mut volume = 0.0;
    for t in 0..n {
        let [a, b, c] = s.corners(t);
        area += 0.5 * dot(cross(sub(b, a), sub(c, a)), cross(sub(b, a), sub(c, a))).sqrt();
        volume += dot(a, cross(b, c)) / 6.0;
    }

    // degenerate triangles: repeated corner or exactly collinear corners
    let degenerate_all: Vec<usize> = (0..n).filter(|&t| degenerate(s, t)).collect();

    // duplicate triangles: the same set of corners
    let mut by_set: HashMap<[u32; 3], usize> = HashMap::new();
    let mut dup_all: Vec<[usize; 2]> = Vec::new();
    for (t, tri) in s.triangles.iter().enumerate() {
        let mut k = *tri;
        k.sort_unstable();
        if let Some(&first) = by_set.get(&k) {
            dup_all.push([first, t]);
        } else {
            by_set.insert(k, t);
        }
    }

    // edges: (low, high) -> (triangles, runs low -> high)
    let mut edges: HashMap<(u32, u32), (u32, u32)> = HashMap::new();
    let mut tri_edges: Vec<[(u32, u32); 3]> = Vec::with_capacity(n);
    for tri in &s.triangles {
        let mut te = [(0, 0); 3];
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            let key = (a.min(b), a.max(b));
            let e = edges.entry(key).or_insert((0, 0));
            e.0 += 1;
            if a < b {
                e.1 += 1;
            }
            te[k] = key;
        }
        tri_edges.push(te);
    }
    let (mut boundary, mut non_manifold, mut flipped) = (0, 0, 0);
    for (&(a, b), &(count, fwd)) in &edges {
        if a == b {
            continue;
        }
        match count {
            1 => boundary += 1,
            2 => {
                if fwd != 1 {
                    flipped += 1;
                }
            }
            _ => non_manifold += 1,
        }
    }

    // components by union–find over shared edges
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut owner: HashMap<(u32, u32), usize> = HashMap::new();
    for (t, te) in tri_edges.iter().enumerate() {
        for e in te {
            if let Some(&o) = owner.get(e) {
                let (ra, rb) = (find(&mut parent, o), find(&mut parent, t));
                if ra != rb {
                    parent[ra.max(rb)] = ra.min(rb);
                }
            } else {
                owner.insert(*e, t);
            }
        }
    }
    let components = (0..n).filter(|&t| find(&mut parent, t) == t).count();

    let closed = n > 0 && boundary == 0 && non_manifold == 0 && flipped == 0;
    let intersecting_all = self_intersections(s, &degenerate_all, &dup_all);
    let near = near_duplicates(s, 1e-9 * diagonal);

    let mut findings = Vec::new();
    if n == 0 {
        findings.push(Finding {
            code: "IMM001",
            message: "the surface has no triangles".into(),
        });
    }
    if s.vertices.iter().any(|v| v.iter().any(|x| !x.is_finite())) {
        findings.push(Finding {
            code: "IMM001",
            message: "corners with non-finite coordinates".into(),
        });
    }
    if !intersecting_all.is_empty() {
        findings.push(Finding {
            code: "IMM005",
            message: format!(
                "{} pairs of triangles intersect (first: {:?}); inside/outside comes from the winding number",
                intersecting_all.len(),
                &intersecting_all[..intersecting_all.len().min(3)]
            ),
        });
    }
    if boundary > 0 {
        findings.push(Finding {
            code: "IMM002",
            message: format!(
                "open surface: {boundary} edges with one triangle (holes); the winding number closes them smoothly"
            ),
        });
    }
    if non_manifold > 0 {
        findings.push(Finding {
            code: "IMM003",
            message: format!("{non_manifold} edges shared by more than two triangles"),
        });
    }
    if flipped > 0 {
        findings.push(Finding {
            code: "IMM004",
            message: format!(
                "{flipped} edges whose two triangles run them in the same direction (inconsistent orientation)"
            ),
        });
    }
    if !degenerate_all.is_empty() {
        findings.push(Finding {
            code: "IMM006",
            message: format!(
                "{} degenerate triangles (collinear or repeated corners)",
                degenerate_all.len()
            ),
        });
    }
    if !dup_all.is_empty() {
        findings.push(Finding {
            code: "IMM007",
            message: format!("{} duplicate triangles", dup_all.len()),
        });
    }
    if closed && volume < 0.0 {
        findings.push(Finding {
            code: "IMM008",
            message: format!(
                "the normals point inward (enclosed volume {volume:.6e} < 0): the body would be its complement"
            ),
        });
    }
    if diagonal > 20.0 {
        findings.push(Finding {
            code: "IMM009",
            message: format!(
                "the bounding box diagonal is {diagonal:.4} m: a file in millimetres needs scale = 1e-3"
            ),
        });
    }
    if near > 0 {
        findings.push(Finding {
            code: "IMM009",
            message: format!(
                "{near} pairs of distinct corners closer than 1e-9 of the box diagonal (not welded)"
            ),
        });
    }

    SurfaceReport {
        triangles: n,
        vertices: s.vertices.len(),
        lower,
        upper,
        diagonal,
        area,
        volume,
        components,
        closed,
        boundary_edges: boundary,
        non_manifold_edges: non_manifold,
        flipped_edges: flipped,
        n_degenerate: degenerate_all.len(),
        degenerate: degenerate_all.into_iter().take(LISTED).collect(),
        n_duplicates: dup_all.len(),
        duplicates: dup_all.into_iter().take(LISTED).collect(),
        n_intersecting: intersecting_all.len(),
        intersecting: intersecting_all.into_iter().take(LISTED).collect(),
        near_duplicate_vertices: near,
        findings,
    }
}

/// Repeated or exactly collinear corners.
fn degenerate(s: &Surface, t: usize) -> bool {
    let [i, j, k] = s.triangles[t];
    if i == j || j == k || k == i {
        return true;
    }
    let [a, b, c] = s.corners(t);
    // collinear in 3D ⇔ collinear in all three coordinate projections
    (0..3).all(|drop| orient2d(proj(a, drop), proj(b, drop), proj(c, drop)) == 0.0)
}

/// The point without coordinate `drop`.
fn proj(p: [f64; 3], drop: usize) -> [f64; 2] {
    match drop {
        0 => [p[1], p[2]],
        1 => [p[2], p[0]],
        _ => [p[0], p[1]],
    }
}

/// The coordinate to drop for a 2D view of the plane of `t` that keeps it
/// non-degenerate: the one where the projected triangle is not collinear
/// (largest normal component first).
fn drop_axis(t: &[[f64; 3]; 3]) -> usize {
    let nrm = cross(sub(t[1], t[0]), sub(t[2], t[0]));
    let mut order = [0usize, 1, 2];
    order.sort_by(|&a, &b| nrm[b].abs().total_cmp(&nrm[a].abs()));
    for d in order {
        if orient2d(proj(t[0], d), proj(t[1], d), proj(t[2], d)) != 0.0 {
            return d;
        }
    }
    order[0]
}

fn sign(x: f64) -> i8 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

/// Closed segments `p q` and `r s` in 2D intersect.
fn seg_seg_2d(p: [f64; 2], q: [f64; 2], r: [f64; 2], s: [f64; 2]) -> bool {
    let d1 = sign(orient2d(p, q, r));
    let d2 = sign(orient2d(p, q, s));
    let d3 = sign(orient2d(r, s, p));
    let d4 = sign(orient2d(r, s, q));
    if d1 * d2 < 0 && d3 * d4 < 0 {
        return true;
    }
    let on = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        // c on the closed segment a b, given collinear
        a[0].min(b[0]) <= c[0]
            && c[0] <= a[0].max(b[0])
            && a[1].min(b[1]) <= c[1]
            && c[1] <= a[1].max(b[1])
    };
    (d1 == 0 && on(p, q, r))
        || (d2 == 0 && on(p, q, s))
        || (d3 == 0 && on(r, s, p))
        || (d4 == 0 && on(r, s, q))
}

/// Point `p` in the closed 2D triangle `t`.
fn point_in_tri_2d(p: [f64; 2], t: &[[f64; 2]; 3]) -> bool {
    let a = sign(orient2d(t[0], t[1], p));
    let b = sign(orient2d(t[1], t[2], p));
    let c = sign(orient2d(t[2], t[0], p));
    !((a < 0 || b < 0 || c < 0) && (a > 0 || b > 0 || c > 0))
}

/// Closed segment `p q` and closed triangle `t` in 3D.
fn seg_tri(p: [f64; 3], q: [f64; 3], t: &[[f64; 3]; 3]) -> bool {
    let op = sign(orient3d(t[0], t[1], t[2], p));
    let oq = sign(orient3d(t[0], t[1], t[2], q));
    if op * oq > 0 {
        return false;
    }
    if op == 0 && oq == 0 {
        let d = drop_axis(t);
        let t2 = t.map(|v| proj(v, d));
        let (p2, q2) = (proj(p, d), proj(q, d));
        return point_in_tri_2d(p2, &t2)
            || point_in_tri_2d(q2, &t2)
            || (0..3).any(|k| seg_seg_2d(p2, q2, t2[k], t2[(k + 1) % 3]));
    }
    let wa = sign(orient3d(p, q, t[1], t[2]));
    let wb = sign(orient3d(p, q, t[2], t[0]));
    let wc = sign(orient3d(p, q, t[0], t[1]));
    !((wa < 0 || wb < 0 || wc < 0) && (wa > 0 || wb > 0 || wc > 0))
}

/// Exactly coplanar triangles.
fn coplanar(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> bool {
    b.iter().all(|&p| orient3d(a[0], a[1], a[2], p) == 0.0)
}

/// The direction from `v` to `p` lies in the closed wedge of the 2D
/// triangle (v, b1, b2) at v, or along one of its sides.
fn in_wedge(v: [f64; 2], b1: [f64; 2], b2: [f64; 2], p: [f64; 2]) -> bool {
    let s = sign(orient2d(v, b1, b2));
    if s == 0 {
        return false;
    }
    let o1 = sign(orient2d(v, b1, p)) * s;
    let o2 = sign(orient2d(v, p, b2)) * s;
    let same_dir = |b: [f64; 2]| {
        // p collinear with v b: on the same side of v
        (0..2).any(|k| (b[k] - v[k]) * (p[k] - v[k]) > 0.0)
    };
    (o1 > 0 && o2 > 0) || (o1 == 0 && same_dir(b1)) || (o2 == 0 && same_dir(b2))
}

/// Exact intersection of triangles `i` and `j` (both non-degenerate,
/// not duplicates).
fn intersect(s: &Surface, i: usize, j: usize) -> bool {
    let (ti, tj) = (s.triangles[i], s.triangles[j]);
    let (a, b) = (s.corners(i), s.corners(j));
    let shared: Vec<(usize, usize)> = (0..3)
        .flat_map(|x| (0..3).map(move |y| (x, y)))
        .filter(|&(x, y)| ti[x] == tj[y])
        .collect();
    match shared.len() {
        0 => {
            if coplanar(&a, &b) {
                let d = drop_axis(&a);
                let (a2, b2) = (a.map(|v| proj(v, d)), b.map(|v| proj(v, d)));
                (0..3).any(|x| {
                    (0..3).any(|y| seg_seg_2d(a2[x], a2[(x + 1) % 3], b2[y], b2[(y + 1) % 3]))
                }) || point_in_tri_2d(a2[0], &b2)
                    || point_in_tri_2d(b2[0], &a2)
            } else {
                (0..3).any(|x| seg_tri(a[x], a[(x + 1) % 3], &b))
                    || (0..3).any(|y| seg_tri(b[y], b[(y + 1) % 3], &a))
            }
        }
        1 => {
            let (x, y) = shared[0];
            let (a1, a2) = (a[(x + 1) % 3], a[(x + 2) % 3]);
            let (b1, b2) = (b[(y + 1) % 3], b[(y + 2) % 3]);
            if coplanar(&a, &b) {
                let d = drop_axis(&a);
                let v = proj(a[x], d);
                let (a1, a2, b1, b2) = (proj(a1, d), proj(a2, d), proj(b1, d), proj(b2, d));
                in_wedge(v, b1, b2, a1)
                    || in_wedge(v, b1, b2, a2)
                    || in_wedge(v, a1, a2, b1)
                    || in_wedge(v, a1, a2, b2)
            } else {
                seg_tri(a1, a2, &b) || seg_tri(b1, b2, &a)
            }
        }
        2 => {
            // shared edge: a fold when coplanar with the third corners on
            // the same side of it
            if !coplanar(&a, &b) {
                return false;
            }
            let xa = (0..3)
                .find(|&x| shared.iter().all(|&(sx, _)| sx != x))
                .unwrap_or(0);
            let yb = (0..3)
                .find(|&y| shared.iter().all(|&(_, sy)| sy != y))
                .unwrap_or(0);
            let (u, w) = (a[shared[0].0], a[shared[1].0]);
            let d = drop_axis(&a);
            let (u, w) = (proj(u, d), proj(w, d));
            let sa = sign(orient2d(u, w, proj(a[xa], d)));
            let sb = sign(orient2d(u, w, proj(b[yb], d)));
            sa * sb >= 0
        }
        _ => false,
    }
}

/// Pairs of intersecting triangles, sorted, skipping degenerate and
/// duplicate triangles (reported on their own).
fn self_intersections(s: &Surface, degenerate: &[usize], dups: &[[usize; 2]]) -> Vec<[usize; 2]> {
    let n = s.len();
    let mut skip = vec![false; n];
    for &t in degenerate {
        skip[t] = true;
    }
    for d in dups {
        skip[d[1]] = true;
    }
    let bvh = &s.bvh;
    if bvh.nodes.is_empty() {
        return Vec::new();
    }
    let mut pairs: Vec<[usize; 2]> = (0..n)
        .into_par_iter()
        .filter(|&i| !skip[i])
        .flat_map_iter(|i| {
            let t = s.corners(i);
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for p in t {
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
            let mut out = Vec::new();
            let mut stack = vec![0u32];
            while let Some(nd) = stack.pop() {
                let node = bvh.nodes[nd as usize];
                let bx = node.bx;
                if (0..3).any(|k| bx.hi[k] < lo[k] || bx.lo[k] > hi[k]) {
                    continue;
                }
                if node.count > 0 {
                    let st = node.first as usize;
                    for &j in &bvh.order[st..st + node.count as usize] {
                        let j = j as usize;
                        if j > i && !skip[j] && intersect(s, i, j) {
                            out.push([i, j]);
                        }
                    }
                } else {
                    stack.push(node.first);
                    stack.push(node.first + 1);
                }
            }
            out
        })
        .collect();
    pairs.sort_unstable();
    pairs
}

/// Pairs of distinct corners closer than `tol` (grid hashing).
fn near_duplicates(s: &Surface, tol: f64) -> usize {
    if tol <= 0.0 || s.vertices.len() < 2 {
        return 0;
    }
    let key = |p: [f64; 3]| p.map(|x| (x / tol).floor() as i64);
    let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    for (i, &v) in s.vertices.iter().enumerate() {
        grid.entry(key(v)).or_default().push(i as u32);
    }
    let mut count = 0;
    for (i, &v) in s.vertices.iter().enumerate() {
        let k = key(v);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(list) = grid.get(&[k[0] + dx, k[1] + dy, k[2] + dz]) {
                        for &j in list {
                            let j = j as usize;
                            if j > i {
                                let d = sub(s.vertices[j], v);
                                if dot(d, d) <= tol * tol {
                                    count += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    count
}
