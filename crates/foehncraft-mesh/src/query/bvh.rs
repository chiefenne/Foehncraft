//! Bounding volume hierarchy over surface faces and the closest-point
//! query (D108 part 1), for the exact wall distance (D48).
//!
//! Each face is split into the triangles of the face geometry
//! ([`crate::geometry::polygon_area_centroid`]: one triangle for a
//! triangle, otherwise the fan around the node average), so the surface
//! the query sees is the surface whose areas and centroids the solver uses.
//!
//! **Build** (host): top-down with the surface-area heuristic over binned
//! triangle centroids (MacDonald & Booth 1990; the binned form of Wald
//! 2007), 16 bins per axis, leaves of at most [`LEAF`] triangles. The
//! build orders with total orders only (bin of the centroid coordinate,
//! then triangle index), so the tree is the same on every machine. The
//! same build makes the hierarchy over cell boxes of
//! [`super::cells::CellTree`].
//!
//! **Query**: depth-first with an explicit stack, the nearer child first,
//! a subtree skipped when its box is farther than the best distance so
//! far. The closest point on a triangle follows the Voronoi-region tests
//! of Ericson (*Real-Time Collision Detection*, 2005, section 5.1.5).
//! Among faces at the same distance the smallest global face id wins, and
//! boxes at exactly the best distance are still visited, so the answer does
//! not depend on the traversal order (D93).

// coordinates indexed by axis read more clearly than zipped iterators
#![allow(clippy::needless_range_loop)]

use crate::geometry::{cross, dot, node_average, norm, scale, sub};
use crate::{Geometry, Mesh};
use rayon::prelude::*;

/// Most triangles in a leaf.
pub const LEAF: usize = 4;
/// Bins of the surface-area heuristic per axis.
const BINS: usize = 16;

/// The faces of a surface as triangles.
#[derive(Debug, Clone, Default)]
pub struct SurfaceFaces {
    /// Triangle corners.
    pub tris: Vec<[[f64; 3]; 3]>,
    /// Face (index into `global_id` and `normal`) of each triangle.
    pub tri_face: Vec<u32>,
    /// Global id of each face (ties are broken by the smallest).
    pub global_id: Vec<u64>,
    /// Local mesh face of each face.
    pub mesh_face: Vec<usize>,
    /// Unit normal of each face, into the domain.
    pub normal: Vec<[f64; 3]>,
}

impl SurfaceFaces {
    /// The boundary faces of `mesh` on the boundaries `boundaries`
    /// (indices into `mesh.boundaries`), in mesh face order.
    pub fn of_boundaries(mesh: &Mesh, geo: &Geometry, boundaries: &[usize]) -> Self {
        let mut s = SurfaceFaces::default();
        let mut list: Vec<usize> = boundaries.to_vec();
        list.sort_unstable();
        list.dedup();
        for b in list {
            for f in mesh.boundaries[b].range() {
                s.push_face(mesh, geo, f);
            }
        }
        s
    }

    /// Appends the mesh face `f` (its triangles, global id and normal).
    pub fn push_face(&mut self, mesh: &Mesh, geo: &Geometry, f: usize) {
        let id = self.global_id.len() as u32;
        let nodes = mesh.face(f);
        let p = |i: usize| mesh.points[nodes[i] as usize];
        if nodes.len() == 3 {
            self.tris.push([p(0), p(1), p(2)]);
            self.tri_face.push(id);
        } else {
            let m = node_average(&mesh.points, nodes);
            for i in 0..nodes.len() {
                self.tris.push([p(i), p((i + 1) % nodes.len()), m]);
                self.tri_face.push(id);
            }
        }
        let s = geo.face_area[f];
        let a = norm(s);
        let n = if a > 0.0 {
            scale(s, -1.0 / a)
        } else {
            [0.0; 3]
        };
        self.global_id
            .push(mesh.face_id.get(f).copied().unwrap_or(f as u64));
        self.mesh_face.push(f);
        self.normal.push(n);
    }

    /// Number of faces.
    pub fn len(&self) -> usize {
        self.global_id.len()
    }

    /// No faces.
    pub fn is_empty(&self) -> bool {
        self.global_id.is_empty()
    }
}

/// An axis-aligned box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Aabb {
    pub(crate) lo: [f64; 3],
    pub(crate) hi: [f64; 3],
}

impl Aabb {
    pub(crate) const EMPTY: Aabb = Aabb {
        lo: [f64::INFINITY; 3],
        hi: [f64::NEG_INFINITY; 3],
    };

    pub(crate) fn grow(&mut self, p: [f64; 3]) {
        for i in 0..3 {
            self.lo[i] = self.lo[i].min(p[i]);
            self.hi[i] = self.hi[i].max(p[i]);
        }
    }

    pub(crate) fn merge(&mut self, o: &Aabb) {
        for i in 0..3 {
            self.lo[i] = self.lo[i].min(o.lo[i]);
            self.hi[i] = self.hi[i].max(o.hi[i]);
        }
    }

    fn area(&self) -> f64 {
        let d = sub(self.hi, self.lo);
        if d[0] < 0.0 {
            return 0.0;
        }
        2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])
    }

    /// Squared distance from `x` to the box (0 inside).
    fn dist2(&self, x: [f64; 3]) -> f64 {
        let mut s = 0.0;
        for i in 0..3 {
            let d = if x[i] < self.lo[i] {
                self.lo[i] - x[i]
            } else if x[i] > self.hi[i] {
                x[i] - self.hi[i]
            } else {
                0.0
            };
            s += d * d;
        }
        s
    }
}

/// A node: an inner node holds its two children, a leaf a range of the
/// item order.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Node {
    pub(crate) bx: Aabb,
    /// Inner: index of the left child (the right one follows it). Leaf:
    /// start in `order`.
    pub(crate) first: u32,
    /// 0 for an inner node, else the number of items of the leaf.
    pub(crate) count: u32,
}

/// Builds the nodes over items with centroids `cent` and boxes `bx`
/// (binned surface-area heuristic, at most [`LEAF`] items per leaf);
/// returns the nodes (root first) and the items in leaf order. The result
/// depends only on the input values and their order (D93).
pub(crate) fn build_nodes(cent: &[[f64; 3]], bx: &[Aabb]) -> (Vec<Node>, Vec<u32>) {
    let n = cent.len();
    let mut order: Vec<u32> = (0..n as u32).collect();
    let mut nodes = Vec::with_capacity(2 * n.max(1));
    if n > 0 {
        nodes.push(Node {
            bx: Aabb::EMPTY,
            first: 0,
            count: 0,
        });
        // (node, start, end) still to split
        let mut todo = vec![(0usize, 0usize, n)];
        let mut scratch = Vec::new();
        while let Some((id, s, e)) = todo.pop() {
            let mut b = Aabb::EMPTY;
            for &t in &order[s..e] {
                b.merge(&bx[t as usize]);
            }
            nodes[id].bx = b;
            let split = if e - s > LEAF {
                split(&mut order[s..e], cent, bx, &mut scratch)
            } else {
                None
            };
            match split {
                Some(m) => {
                    let left = nodes.len();
                    for _ in 0..2 {
                        nodes.push(Node {
                            bx: Aabb::EMPTY,
                            first: 0,
                            count: 0,
                        });
                    }
                    nodes[id].first = left as u32;
                    nodes[id].count = 0;
                    todo.push((left + 1, s + m, e));
                    todo.push((left, s, s + m));
                }
                None => {
                    nodes[id].first = s as u32;
                    nodes[id].count = (e - s) as u32;
                }
            }
        }
    }
    (nodes, order)
}

/// The closest point of a surface to a query point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Closest {
    /// Euclidean distance (m).
    pub distance: f64,
    /// The closest point.
    pub point: [f64; 3],
    /// The face (index into the [`SurfaceFaces`]).
    pub face: usize,
    /// Global id of the face.
    pub global_id: u64,
    /// Unit normal of the face, into the domain.
    pub normal: [f64; 3],
}

/// Bounding volume hierarchy over the triangles of a [`SurfaceFaces`].
#[derive(Debug, Clone)]
pub struct Bvh {
    faces: SurfaceFaces,
    nodes: Vec<Node>,
    /// Triangles in leaf order.
    order: Vec<u32>,
}

impl Bvh {
    /// Builds the hierarchy (surface-area heuristic, binned).
    pub fn build(faces: SurfaceFaces) -> Self {
        let cent: Vec<[f64; 3]> = faces
            .tris
            .iter()
            .map(|t| {
                [
                    (t[0][0] + t[1][0] + t[2][0]) / 3.0,
                    (t[0][1] + t[1][1] + t[2][1]) / 3.0,
                    (t[0][2] + t[1][2] + t[2][2]) / 3.0,
                ]
            })
            .collect();
        let tbox: Vec<Aabb> = faces
            .tris
            .iter()
            .map(|t| {
                let mut b = Aabb::EMPTY;
                for &p in t {
                    b.grow(p);
                }
                b
            })
            .collect();
        let (nodes, order) = build_nodes(&cent, &tbox);
        Bvh {
            faces,
            nodes,
            order,
        }
    }

    /// The faces of the hierarchy.
    pub fn faces(&self) -> &SurfaceFaces {
        &self.faces
    }

    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// The closest point of the surface to `x`; `None` for an empty
    /// surface. Among faces at the same distance the one with the smallest
    /// global id.
    pub fn closest(&self, x: [f64; 3]) -> Option<Closest> {
        self.closest_within(x, f64::INFINITY)
    }

    /// The closest point of the surface to `x` at a distance of at most
    /// `radius`, or `None` when no point is that near (a capped query:
    /// boxes farther than the radius are never opened). With an infinite
    /// radius this is [`closest`](Self::closest).
    pub fn closest_within(&self, x: [f64; 3], radius: f64) -> Option<Closest> {
        if self.nodes.is_empty() {
            return None;
        }
        // best: (d², global id, triangle, point)
        let mut best_d2 = radius * radius;
        let mut best_id = u64::MAX;
        let mut best_tri = usize::MAX;
        let mut best_p = [0.0; 3];
        let mut stack: Vec<u32> = Vec::with_capacity(64);
        stack.push(0);
        while let Some(i) = stack.pop() {
            let nd = self.nodes[i as usize];
            if nd.bx.dist2(x) > best_d2 {
                continue;
            }
            if nd.count > 0 {
                let s = nd.first as usize;
                for &t in &self.order[s..s + nd.count as usize] {
                    let t = t as usize;
                    let p = closest_on_triangle(x, &self.faces.tris[t]);
                    let d = sub(x, p);
                    let d2 = dot(d, d);
                    let id = self.faces.global_id[self.faces.tri_face[t] as usize];
                    if d2 < best_d2 || (d2 == best_d2 && (id, t) < (best_id, best_tri)) {
                        best_d2 = d2;
                        best_id = id;
                        best_tri = t;
                        best_p = p;
                    }
                }
            } else {
                let (l, r) = (nd.first, nd.first + 1);
                let dl = self.nodes[l as usize].bx.dist2(x);
                let dr = self.nodes[r as usize].bx.dist2(x);
                // the nearer child is visited first (pushed last)
                if dl <= dr {
                    stack.push(r);
                    stack.push(l);
                } else {
                    stack.push(l);
                    stack.push(r);
                }
            }
        }
        if best_tri == usize::MAX {
            return None;
        }
        let face = self.faces.tri_face[best_tri] as usize;
        Some(Closest {
            distance: best_d2.sqrt(),
            point: best_p,
            face,
            global_id: best_id,
            normal: self.faces.normal[face],
        })
    }

    /// The closest point by testing every triangle (for the tests).
    pub fn brute_force(&self, x: [f64; 3]) -> Option<Closest> {
        let mut best: Option<(f64, u64, usize, [f64; 3])> = None;
        for (t, tri) in self.faces.tris.iter().enumerate() {
            let p = closest_on_triangle(x, tri);
            let d = sub(x, p);
            let d2 = dot(d, d);
            let id = self.faces.global_id[self.faces.tri_face[t] as usize];
            if best.is_none_or(|(bd, bid, bt, _)| d2 < bd || (d2 == bd && (id, t) < (bid, bt))) {
                best = Some((d2, id, t, p));
            }
        }
        best.map(|(d2, id, t, p)| {
            let face = self.faces.tri_face[t] as usize;
            Closest {
                distance: d2.sqrt(),
                point: p,
                face,
                global_id: id,
                normal: self.faces.normal[face],
            }
        })
    }
}

/// Orders `order` (a range of items, ascending on entry) by the best split
/// of the binned surface-area heuristic and returns the size of the left
/// part, or `None` when no split is cheaper than a leaf. Both parts stay
/// ascending: the split is a stable partition of an ascending range, the
/// same order as a sort by (side, index), in linear time.
fn split(
    order: &mut [u32],
    cent: &[[f64; 3]],
    tbox: &[Aabb],
    scratch: &mut Vec<u32>,
) -> Option<usize> {
    let n = order.len();
    let mut cb = Aabb::EMPTY;
    for &t in order.iter() {
        cb.grow(cent[t as usize]);
    }
    // cost of a leaf (relative units: one triangle test = 1, one box = 1)
    let mut bx = Aabb::EMPTY;
    for &t in order.iter() {
        bx.merge(&tbox[t as usize]);
    }
    let parent = bx.area().max(f64::MIN_POSITIVE);
    let mut best: Option<(f64, usize, usize)> = None; // (cost, axis, bin)
    for axis in 0..3 {
        let (lo, hi) = (cb.lo[axis], cb.hi[axis]);
        if hi <= lo {
            continue;
        }
        let w = (hi - lo) / BINS as f64;
        let bin_of = |c: f64| (((c - lo) / w) as usize).min(BINS - 1);
        let mut counts = [0usize; BINS];
        let mut boxes = [Aabb::EMPTY; BINS];
        for &t in order.iter() {
            let b = bin_of(cent[t as usize][axis]);
            counts[b] += 1;
            boxes[b].merge(&tbox[t as usize]);
        }
        // sweep: area and count of the left part of each split position
        let mut left_a = [0.0; BINS];
        let mut left_n = [0usize; BINS];
        let (mut acc, mut cnt) = (Aabb::EMPTY, 0);
        for b in 0..BINS - 1 {
            acc.merge(&boxes[b]);
            cnt += counts[b];
            left_a[b] = acc.area();
            left_n[b] = cnt;
        }
        let (mut acc, mut cnt) = (Aabb::EMPTY, 0);
        for b in (1..BINS).rev() {
            acc.merge(&boxes[b]);
            cnt += counts[b];
            let nl = left_n[b - 1];
            if nl == 0 || cnt == 0 {
                continue;
            }
            let cost = 1.0 + (left_a[b - 1] * nl as f64 + acc.area() * cnt as f64) / parent;
            if best.is_none_or(|(c, _, _)| cost < c) {
                best = Some((cost, axis, b - 1));
            }
        }
    }
    let leaf_cost = n as f64;
    match best {
        Some((cost, axis, bin)) if cost < leaf_cost || n > 4 * LEAF => {
            let (lo, hi) = (cb.lo[axis], cb.hi[axis]);
            let w = (hi - lo) / BINS as f64;
            let bin_of = |c: f64| (((c - lo) / w) as usize).min(BINS - 1);
            // a stable total order: left bins first, then by index
            debug_assert!(order.windows(2).all(|w| w[0] < w[1]));
            scratch.clear();
            let mut m = 0;
            for i in 0..n {
                let t = order[i];
                if bin_of(cent[t as usize][axis]) <= bin {
                    order[m] = t;
                    m += 1;
                } else {
                    scratch.push(t);
                }
            }
            order[m..].copy_from_slice(scratch);
            Some(m)
        }
        Some(_) => None,
        None if n > 4 * LEAF => {
            // all centroids coincide: split the index range in halves
            order.sort_unstable();
            Some(n / 2)
        }
        None => None,
    }
}

/// The closest point of triangle `t` to `x` (Ericson 2005, 5.1.5).
pub fn closest_on_triangle(x: [f64; 3], t: &[[f64; 3]; 3]) -> [f64; 3] {
    let (a, b, c) = (t[0], t[1], t[2]);
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(x, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(x, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add_scaled(a, ab, v);
    }
    let cp = sub(x, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add_scaled(a, ac, w);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add_scaled(b, sub(c, b), w);
    }
    let den = va + vb + vc;
    if den == 0.0 {
        // a degenerate triangle (zero area): the nearest of its edges
        let cands = [seg(x, a, b), seg(x, b, c), seg(x, c, a)];
        let mut best = cands[0];
        let mut bd = dot(sub(x, best), sub(x, best));
        for &p in &cands[1..] {
            let d = dot(sub(x, p), sub(x, p));
            if d < bd {
                bd = d;
                best = p;
            }
        }
        return best;
    }
    // inside the face region: the projection onto the plane of the
    // triangle, x − n (n·(x − a))/|n|² (exact for a wall in a coordinate
    // plane, where the barycentric form rounds)
    let n = cross(ab, ac);
    let s = dot(n, ap) / dot(n, n);
    [x[0] - n[0] * s, x[1] - n[1] * s, x[2] - n[2] * s]
}

fn add_scaled(a: [f64; 3], d: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] + d[0] * s, a[1] + d[1] * s, a[2] + d[2] * s]
}

/// The closest point of segment a–b to x.
fn seg(x: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let d = sub(b, a);
    let dd = dot(d, d);
    if dd == 0.0 {
        return a;
    }
    let s = (dot(sub(x, a), d) / dd).clamp(0.0, 1.0);
    add_scaled(a, d, s)
}

/// The exact wall distance of D48 at the cell centroids of `mesh`: the
/// closest point on the faces of the boundaries `walls` (indices into
/// `mesh.boundaries`) per local cell, or `None` without wall faces. The
/// cells are independent queries (any number of threads gives the same
/// result).
pub fn wall_distance(mesh: &Mesh, geo: &Geometry, walls: &[usize]) -> Option<Vec<Closest>> {
    let faces = SurfaceFaces::of_boundaries(mesh, geo, walls);
    if faces.is_empty() {
        return None;
    }
    let bvh = Bvh::build(faces);
    let nc = mesh.n_cells;
    Some(
        (0..nc)
            .into_par_iter()
            .map(|c| {
                bvh.closest(geo.cell_centroid[c])
                    .expect("a non-empty surface")
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(s: &mut u64) -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*s >> 11) as f64) / ((1u64 << 53) as f64)
    }

    #[test]
    fn closest_point_on_a_triangle_by_regions() {
        let t = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert_eq!(closest_on_triangle([0.2, 0.2, 3.0], &t), [0.2, 0.2, 0.0]);
        assert_eq!(closest_on_triangle([-1.0, -1.0, 0.0], &t), [0.0, 0.0, 0.0]);
        assert_eq!(closest_on_triangle([2.0, -0.5, 1.0], &t), [1.0, 0.0, 0.0]);
        assert_eq!(closest_on_triangle([0.5, -1.0, 0.0], &t), [0.5, 0.0, 0.0]);
        assert_eq!(closest_on_triangle([-1.0, 0.5, 0.0], &t), [0.0, 0.5, 0.0]);
        let p = closest_on_triangle([1.0, 1.0, 0.0], &t);
        assert!((p[0] - 0.5).abs() < 1e-15 && (p[1] - 0.5).abs() < 1e-15);
        // degenerate (collinear) triangle
        let d = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        assert_eq!(closest_on_triangle([1.5, 1.0, 0.0], &d), [1.5, 0.0, 0.0]);
    }

    /// A random soup of triangles: the hierarchy gives the brute-force
    /// answer bitwise, ties included.
    #[test]
    fn hierarchy_agrees_with_brute_force() {
        let mut s = 7u64;
        let mut f = SurfaceFaces::default();
        for i in 0..500 {
            let c = [lcg(&mut s), lcg(&mut s), lcg(&mut s)];
            let mut t = [[0.0; 3]; 3];
            for p in &mut t {
                for k in 0..3 {
                    p[k] = c[k] + 0.05 * (lcg(&mut s) - 0.5);
                }
            }
            f.tris.push(t);
            f.tri_face.push(i);
            f.global_id.push(1000 - i as u64);
            f.mesh_face.push(i as usize);
            f.normal.push([0.0, 0.0, 1.0]);
        }
        // a duplicate triangle with a smaller global id: ties go to it
        let t0 = f.tris[17];
        f.tris.push(t0);
        f.tri_face.push(500);
        f.global_id.push(1);
        f.mesh_face.push(500);
        f.normal.push([0.0, 0.0, 1.0]);
        let bvh = Bvh::build(f);
        assert!(bvh.n_nodes() > 1);
        for _ in 0..2000 {
            let x = [
                2.0 * lcg(&mut s) - 0.5,
                2.0 * lcg(&mut s) - 0.5,
                2.0 * lcg(&mut s) - 0.5,
            ];
            assert_eq!(bvh.closest(x), bvh.brute_force(x));
        }
        let c = bvh.closest(f_centroid(&t0)).unwrap();
        assert!(c.distance < 1e-15 && c.global_id == 1, "{c:?}");
    }

    fn f_centroid(t: &[[f64; 3]; 3]) -> [f64; 3] {
        [
            (t[0][0] + t[1][0] + t[2][0]) / 3.0,
            (t[0][1] + t[1][1] + t[2][1]) / 3.0,
            (t[0][2] + t[1][2] + t[2][2]) / 3.0,
        ]
    }
    /// Plane walls of a channel (D48 acceptance, stage 6 spec 9.2): the
    /// distance is min(y, H − y) to round-off on hexahedra, polyhedra and
    /// perturbed tetrahedra (the walls stay planes).
    #[test]
    fn channel_walls_give_the_exact_distance() {
        use crate::generate::{BoxSpec, Kind};
        for (kind, perturb) in [(Kind::Hex, 0.0), (Kind::Poly, 0.0), (Kind::Tet, 0.2)] {
            let mut spec = BoxSpec::new(kind, 3, [6, 8, 4]);
            spec.upper = [2.0, 0.5, 1.0];
            spec.perturb = perturb;
            let mesh = spec.generate().unwrap();
            let geo = Geometry::compute(&mesh);
            let walls: Vec<usize> = ["ymin", "ymax"]
                .iter()
                .map(|n| mesh.boundaries.iter().position(|b| b.name == *n).unwrap())
                .collect();
            let d = wall_distance(&mesh, &geo, &walls).unwrap();
            for (c, w) in d.iter().enumerate() {
                let y = geo.cell_centroid[c][1];
                let exact = y.min(0.5 - y);
                assert!(
                    (w.distance - exact).abs() <= 1e-14,
                    "{kind:?} cell {c}: {} vs {exact}",
                    w.distance
                );
                // the normal points into the domain
                let n = w.normal[1];
                assert!(if y < 0.25 {
                    n > 0.99
                } else if y > 0.25 {
                    n < -0.99
                } else {
                    n.abs() > 0.99
                });
            }
        }
    }

    /// The cylinder of the O-grid: the BVH gives the brute-force answer for
    /// every cell, and the distance to the faceted wall is within the
    /// faceting error of the distance to the circle.
    #[test]
    fn cylinder_wall_distance() {
        let spec = crate::cylinder::CylinderSpec {
            n: 8,
            radial: 8,
            wake: 8,
            ..Default::default()
        };
        let mesh = spec.generate().unwrap();
        let geo = Geometry::compute(&mesh);
        let b = mesh
            .boundaries
            .iter()
            .position(|b| b.name == "cylinder")
            .unwrap();
        let faces = SurfaceFaces::of_boundaries(&mesh, &geo, &[b]);
        let bvh = Bvh::build(faces);
        let r = 0.5 * spec.diameter;
        let facet = r * (1.0 - (std::f64::consts::PI / (4 * spec.n) as f64).cos());
        for c in 0..mesh.n_cells {
            let x = geo.cell_centroid[c];
            let q = bvh.closest(x).unwrap();
            assert_eq!(Some(q), bvh.brute_force(x));
            let rc = ((x[0] - spec.center[0]).powi(2) + (x[1] - spec.center[1]).powi(2)).sqrt();
            assert!(q.distance >= rc - r - 1e-14 && q.distance <= rc - r + facet + 1e-14);
        }
    }

    /// A faceted sphere (octahedron refined four times, nodes on the
    /// sphere): the hierarchy agrees with brute force, and the distance lies
    /// between that to the sphere and that plus the faceting error.
    #[test]
    fn sphere_wall_distance() {
        let mut tris: Vec<[[f64; 3]; 3]> = Vec::new();
        let o = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        for &(a, b) in &[(0, 2), (2, 1), (1, 3), (3, 0)] {
            tris.push([o[a], o[b], o[4]]);
            tris.push([o[b], o[a], o[5]]);
        }
        let unit = |p: [f64; 3]| scale(p, 1.0 / norm(p));
        for _ in 0..4 {
            let mut next = Vec::new();
            for t in &tris {
                let m = |i: usize, j: usize| {
                    unit(scale(
                        [t[i][0] + t[j][0], t[i][1] + t[j][1], t[i][2] + t[j][2]],
                        0.5,
                    ))
                };
                let (ab, bc, ca) = (m(0, 1), m(1, 2), m(2, 0));
                next.push([t[0], ab, ca]);
                next.push([ab, t[1], bc]);
                next.push([ca, bc, t[2]]);
                next.push([ab, bc, ca]);
            }
            tris = next;
        }
        let mut f = SurfaceFaces::default();
        let mut max_edge: f64 = 0.0;
        for (i, t) in tris.iter().enumerate() {
            for k in 0..3 {
                max_edge = max_edge.max(norm(sub(t[k], t[(k + 1) % 3])));
            }
            f.tris.push(*t);
            f.tri_face.push(i as u32);
            f.global_id.push(i as u64);
            f.mesh_face.push(i);
            f.normal.push([0.0; 3]);
        }
        // the circumradius of the facets bounds the gap to the sphere
        let facet = 1.0 - (1.0 - max_edge * max_edge / 3.0).sqrt();
        let bvh = Bvh::build(f);
        let mut s = 11u64;
        for _ in 0..500 {
            let x = [
                4.0 * lcg(&mut s) - 2.0,
                4.0 * lcg(&mut s) - 2.0,
                4.0 * lcg(&mut s) - 2.0,
            ];
            let q = bvh.closest(x).unwrap();
            assert_eq!(Some(q), bvh.brute_force(x));
            let exact = (norm(x) - 1.0).abs();
            assert!(
                (q.distance - exact).abs() <= facet + 1e-14,
                "{} vs {exact}",
                q.distance
            );
        }
    }
}
