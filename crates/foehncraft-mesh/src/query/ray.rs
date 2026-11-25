//! Ray casting on the surface hierarchy (D108 part 2): the first face a
//! segment or a ray meets, and whether anything lies between two points.
//!
//! The segment–triangle test is the classical one from signed volumes
//! (O'Rourke 1998, *Computational Geometry in C*, 2nd ed., §7.3; Segura &
//! Feito 2001): the segment from `p` to `q` crosses the plane of the
//! triangle when `p` and `q` lie on opposite sides of it, and its line
//! passes through the closed triangle when the three tetrahedra of the
//! line with the triangle edges have no two volumes of opposite sign. All
//! five volumes are the exact orientation predicate of
//! [`foehncraft_core::predicates::orient3d`] (Shewchuk 1997, D225), so the
//! decision is exact for the given floating-point points:
//!
//! * **Watertight.** Two triangles that share an edge (the same two corner
//!   points) evaluate the volume of that edge as exact negations, so a
//!   segment that crosses a closed surface meets at least one of its
//!   triangles; through an edge or a corner it meets every triangle there,
//!   at the same parameter, and the smallest global face id wins.
//! * **No self-hits, no offsets.** The segment is half open, (0, 1]: a
//!   start point that lies exactly in the plane of a triangle does not
//!   meet it. A ray started on an emitting wall face therefore never hits
//!   that face or a coplanar neighbour, without the usual offset along the
//!   normal.
//! * A segment that lies in the plane of a triangle does not meet it.
//!
//! Only the hit parameter and the hit point are rounded (from the same
//! volumes, as barycentric weights), the decision is not. The filtered
//! predicate is the fast path, its exact stage the fallback near zero. The
//! test was chosen over the shear-and-scale "watertight" test of graphics
//! because patents carry that title (assessment
//! `docs/assessments/graphics-algorithms-for-modules.md`, sections 3.3 and
//! 8); a floating-point Möller–Trumbore test is the benchmark reference
//! (`examples/ray_bench.rs`).
//!
//! The box test of the traversal is the slab test with the conservative
//! rounding margin of Ize 2013 ("Robust BVH ray traversal", JCGT 2(2)),
//! so it never discards a box that the exact test would enter.

use super::bvh::{Aabb, Bvh};
use crate::geometry::{dot, norm, sub};
use foehncraft_core::predicates::orient3d;

/// Unit roundoff of f64.
const U: f64 = f64::EPSILON / 2.0;
/// Relative margin of the slab test: four times γ₃ = 3u / (1 − 3u) (Ize
/// 2013 shows 2γ₃; the extra factor covers the rounded segment direction).
const MARGIN: f64 = 4.0 * (3.0 * U / (1.0 - 3.0 * U));
/// Boxes whose entry parameter exceeds the best hit by more than this are
/// not opened. The hit parameters are rounded, so a little slack keeps
/// hits at the same true parameter in the comparison by global id.
const SLACK: f64 = 1e-9;

/// The end point of a segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum End {
    /// The end point belongs to the segment, (0, 1]: a segment that ends
    /// on a triangle meets it (parcel paths, rays).
    Closed,
    /// The end point does not belong to it, (0, 1): a segment between two
    /// surface points meets only what lies between them (visibility).
    Open,
}

/// Where a segment crosses one triangle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SegmentHit {
    /// Segment parameter, in (0, 1].
    pub t: f64,
    /// Barycentric weights of the triangle corners (sum 1, none negative).
    pub weights: [f64; 3],
}

/// The first face a segment or ray meets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hit {
    /// Parameter of the hit: for a segment from `p` to `q` the fraction
    /// in (0, 1] (hit point about `p + t (q - p)`), for a ray from `o`
    /// along `d` the ray parameter (hit point about `o + t d`).
    pub t: f64,
    /// The hit point, on the triangle (barycentric form).
    pub point: [f64; 3],
    /// The face (index into the [`SurfaceFaces`](super::SurfaceFaces)).
    pub face: usize,
    /// Global id of the face.
    pub global_id: u64,
    /// Unit normal of the face, into the domain.
    pub normal: [f64; 3],
    /// The segment runs against the face normal, from the domain into the
    /// wall.
    pub entering: bool,
}

/// The exact segment–triangle test: where the segment from `p` to `q`
/// crosses the closed triangle `tri`, or `None`. The segment is half open
/// at `p`; `end` decides whether `q` belongs to it. A segment in the plane
/// of the triangle does not cross it.
pub fn segment_triangle(
    p: [f64; 3],
    q: [f64; 3],
    tri: &[[f64; 3]; 3],
    end: End,
) -> Option<SegmentHit> {
    let [a, b, c] = *tri;
    let vp = orient3d(a, b, c, p);
    if vp == 0.0 {
        return None;
    }
    let vq = orient3d(a, b, c, q);
    if vq == 0.0 {
        if end == End::Open {
            return None;
        }
    } else if (vp > 0.0) == (vq > 0.0) {
        return None;
    }
    // the line p–q against the edges: the weight of a corner is the volume
    // with the opposite edge
    let wa = orient3d(p, q, b, c);
    let wb = orient3d(p, q, c, a);
    if opposite(wa, wb) {
        return None;
    }
    let wc = orient3d(p, q, a, b);
    if opposite(wa, wc) || opposite(wb, wc) {
        return None;
    }
    // p and q on opposite sides: the line is not in the plane, so the
    // three volumes are not all zero and their sum is not zero
    let s = wa + wb + wc;
    let t = if vq == 0.0 {
        1.0
    } else {
        (vp / (vp - vq)).clamp(f64::MIN_POSITIVE, 1.0)
    };
    Some(SegmentHit {
        t,
        weights: [wa / s, wb / s, wc / s],
    })
}

/// Two values of strictly opposite sign.
fn opposite(x: f64, y: f64) -> bool {
    (x > 0.0 && y < 0.0) || (x < 0.0 && y > 0.0)
}

impl Aabb {
    /// The parameter where the segment `p + t d`, 0 ≤ t ≤ 1, enters the
    /// box, or `None` when it misses. Conservative: rounding can only make
    /// a miss look like a hit.
    pub(crate) fn enter(&self, p: [f64; 3], d: [f64; 3]) -> Option<f64> {
        let (mut t0, mut t1) = (0.0f64, 1.0f64);
        for k in 0..3 {
            if d[k] == 0.0 {
                if p[k] < self.lo[k] || p[k] > self.hi[k] {
                    return None;
                }
                continue;
            }
            let a = (self.lo[k] - p[k]) / d[k];
            let b = (self.hi[k] - p[k]) / d[k];
            let (near, far) = if a <= b { (a, b) } else { (b, a) };
            t0 = t0.max(near - near.abs() * MARGIN);
            t1 = t1.min(far + far.abs() * MARGIN);
            if t0 > t1 {
                return None;
            }
        }
        Some(t0)
    }
}

/// The best hit so far: (t, global id, triangle, weights).
type Best = (f64, u64, usize, [f64; 3]);

impl Bvh {
    /// The first face the segment from `p` to `q` meets, in (0, 1] (`p`
    /// excluded, `q` included), or `None`. Among faces met at the same
    /// parameter the one with the smallest global id.
    pub fn first_hit(&self, p: [f64; 3], q: [f64; 3]) -> Option<Hit> {
        let best = self.trace(p, q, End::Closed, false)?;
        Some(self.hit(p, q, best, best.0))
    }

    /// The first face the ray from `o` along `d` meets (`o` excluded), or
    /// `None` (also for `d = 0`). The ray is the segment from `o` to a
    /// point beyond the whole surface; `t` of the hit is the ray parameter.
    pub fn ray_first_hit(&self, o: [f64; 3], d: [f64; 3]) -> Option<Hit> {
        let len = norm(d);
        if self.nodes.is_empty() || len == 0.0 || !len.is_finite() {
            return None;
        }
        let bx = self.nodes[0].bx;
        let mut r2 = 0.0f64;
        for i in 0..8 {
            let corner = std::array::from_fn(|k| if i >> k & 1 == 0 { bx.lo[k] } else { bx.hi[k] });
            let e = sub(corner, o);
            r2 = r2.max(dot(e, e));
        }
        let scale = 2.0 * r2.sqrt() / len + 1.0;
        let q = std::array::from_fn(|k| o[k] + d[k] * scale);
        let best = self.trace(o, q, End::Closed, false)?;
        Some(self.hit(o, q, best, best.0 * scale))
    }

    /// Whether any face lies strictly between `p` and `q` (both end points
    /// excluded): the visibility test between two surface points.
    pub fn occluded(&self, p: [f64; 3], q: [f64; 3]) -> bool {
        self.trace(p, q, End::Open, true).is_some()
    }

    /// Every face the segment from `p` to `q` meets, in (0, 1], ordered
    /// by the parameter (ties by the global id). Hits at the same parameter
    /// (within the slack of the hit test) and in the same direction count once, the face with
    /// the smallest global id: a segment through an edge or a corner of the
    /// surface is one crossing (the spray datums, D263).
    pub fn all_hits(&self, p: [f64; 3], q: [f64; 3]) -> Vec<Hit> {
        let mut found: Vec<Best> = Vec::new();
        if self.nodes.is_empty() {
            return Vec::new();
        }
        let d = sub(q, p);
        let mut stack: Vec<u32> = Vec::with_capacity(64);
        if self.nodes[0].bx.enter(p, d).is_some() {
            stack.push(0);
        }
        while let Some(i) = stack.pop() {
            let nd = self.nodes[i as usize];
            if nd.count > 0 {
                let s = nd.first as usize;
                for &t in &self.order[s..s + nd.count as usize] {
                    let t = t as usize;
                    if let Some(h) = segment_triangle(p, q, &self.faces.tris[t], End::Closed) {
                        let id = self.faces.global_id[self.faces.tri_face[t] as usize];
                        found.push((h.t, id, t, h.weights));
                    }
                }
            } else {
                for c in [nd.first, nd.first + 1] {
                    if self.nodes[c as usize].bx.enter(p, d).is_some() {
                        stack.push(c);
                    }
                }
            }
        }
        found.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        let mut out: Vec<Hit> = Vec::new();
        for b in found {
            let h = self.hit(p, q, b, b.0);
            if let Some(last) = out.last()
                && h.t - last.t <= SLACK
                && h.entering == last.entering
            {
                continue;
            }
            out.push(h);
        }
        out
    }

    /// The first hit by testing every triangle (for the tests and the
    /// benchmark).
    pub fn first_hit_brute_force(&self, p: [f64; 3], q: [f64; 3]) -> Option<Hit> {
        let mut best: Option<Best> = None;
        for (t, tri) in self.faces.tris.iter().enumerate() {
            if let Some(h) = segment_triangle(p, q, tri, End::Closed) {
                let id = self.faces.global_id[self.faces.tri_face[t] as usize];
                if best.is_none_or(|b| (h.t, id, t) < (b.0, b.1, b.2)) {
                    best = Some((h.t, id, t, h.weights));
                }
            }
        }
        best.map(|b| self.hit(p, q, b, b.0))
    }

    /// The traversal: nearer child first, boxes entered after the best hit
    /// (plus [`SLACK`]) skipped; with `any` the first hit found ends it.
    fn trace(&self, p: [f64; 3], q: [f64; 3], end: End, any: bool) -> Option<Best> {
        if self.nodes.is_empty() {
            return None;
        }
        let d = sub(q, p);
        let mut best: Option<Best> = None;
        let mut stack: Vec<(u32, f64)> = Vec::with_capacity(64);
        if let Some(t) = self.nodes[0].bx.enter(p, d) {
            stack.push((0, t));
        }
        while let Some((i, t_in)) = stack.pop() {
            if best.is_some_and(|b| t_in > b.0 + SLACK) {
                continue;
            }
            let nd = self.nodes[i as usize];
            if nd.count > 0 {
                let s = nd.first as usize;
                for &t in &self.order[s..s + nd.count as usize] {
                    let t = t as usize;
                    let Some(h) = segment_triangle(p, q, &self.faces.tris[t], end) else {
                        continue;
                    };
                    let id = self.faces.global_id[self.faces.tri_face[t] as usize];
                    if best.is_none_or(|b| (h.t, id, t) < (b.0, b.1, b.2)) {
                        best = Some((h.t, id, t, h.weights));
                        if any {
                            return best;
                        }
                    }
                }
            } else {
                let (l, r) = (nd.first, nd.first + 1);
                let tl = self.nodes[l as usize].bx.enter(p, d);
                let tr = self.nodes[r as usize].bx.enter(p, d);
                // the nearer child is visited first (pushed last)
                match (tl, tr) {
                    (Some(a), Some(b)) if a <= b => {
                        stack.push((r, b));
                        stack.push((l, a));
                    }
                    (Some(a), Some(b)) => {
                        stack.push((l, a));
                        stack.push((r, b));
                    }
                    (Some(a), None) => stack.push((l, a)),
                    (None, Some(b)) => stack.push((r, b)),
                    (None, None) => {}
                }
            }
        }
        best
    }

    fn hit(&self, p: [f64; 3], q: [f64; 3], b: Best, t: f64) -> Hit {
        let (_, id, tri, w) = b;
        let [a, bb, c] = self.faces.tris[tri];
        let point = std::array::from_fn(|k| w[0] * a[k] + w[1] * bb[k] + w[2] * c[k]);
        let face = self.faces.tri_face[tri] as usize;
        let normal = self.faces.normal[face];
        Hit {
            t,
            point,
            face,
            global_id: id,
            normal,
            entering: dot(sub(q, p), normal) < 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const T: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];

    #[test]
    fn a_segment_through_the_interior() {
        let h = segment_triangle([0.2, 0.3, 1.0], [0.2, 0.3, -3.0], &T, End::Closed).unwrap();
        assert!((h.t - 0.25).abs() < 1e-15);
        let w = h.weights;
        assert!((w[0] - 0.5).abs() < 1e-15 && (w[1] - 0.2).abs() < 1e-15);
        assert!((w[2] - 0.3).abs() < 1e-15);
        // the same from below
        let h = segment_triangle([0.2, 0.3, -1.0], [0.2, 0.3, 1.0], &T, End::Closed).unwrap();
        assert!((h.t - 0.5).abs() < 1e-15);
    }

    #[test]
    fn misses_and_end_points() {
        // beside the triangle, not reaching it, on one side
        assert!(segment_triangle([0.8, 0.8, 1.0], [0.8, 0.8, -1.0], &T, End::Closed).is_none());
        assert!(segment_triangle([0.2, 0.2, 1.0], [0.2, 0.2, 0.5], &T, End::Closed).is_none());
        // starting on the plane: half open at p
        assert!(segment_triangle([0.2, 0.2, 0.0], [0.2, 0.2, -1.0], &T, End::Closed).is_none());
        // ending on the triangle: closed end meets it at t = 1, open not
        let h = segment_triangle([0.2, 0.2, 1.0], [0.2, 0.2, 0.0], &T, End::Closed).unwrap();
        assert_eq!(h.t, 1.0);
        assert!(segment_triangle([0.2, 0.2, 1.0], [0.2, 0.2, 0.0], &T, End::Open).is_none());
        // in the plane
        assert!(segment_triangle([-1.0, 0.2, 0.0], [2.0, 0.2, 0.0], &T, End::Closed).is_none());
    }

    #[test]
    fn corners_and_edges_belong_to_the_triangle() {
        let h = segment_triangle([1.0, 0.0, 1.0], [1.0, 0.0, -1.0], &T, End::Closed).unwrap();
        assert_eq!(h.weights, [0.0, 1.0, 0.0]);
        let h = segment_triangle([0.5, 0.5, 1.0], [0.5, 0.5, -1.0], &T, End::Closed).unwrap();
        assert_eq!(h.weights[0], 0.0);
        assert!((h.weights[1] - 0.5).abs() < 1e-15);
        // the neighbour across the edge b–c meets it too
        let n = [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(segment_triangle([0.5, 0.5, 1.0], [0.5, 0.5, -1.0], &n, End::Closed).is_some());
        // just beside the edge: exactly one of the two
        let p = [0.5 + 1e-17, 0.5 + 1e-16, 1.0];
        let q = [0.5, 0.5, -1.0];
        let a = segment_triangle(p, q, &T, End::Closed).is_some();
        let b = segment_triangle(p, q, &n, End::Closed).is_some();
        assert!(a || b);
    }

    #[test]
    fn the_slab_test_keeps_touching_boxes() {
        let bx = Aabb {
            lo: [0.0; 3],
            hi: [1.0; 3],
        };
        // grazing a face, along an edge, through a corner
        assert!(bx.enter([-1.0, 1.0, 0.5], [3.0, 0.0, 0.0]).is_some());
        assert!(bx.enter([-1.0, 1.0, 1.0], [3.0, 0.0, 0.0]).is_some());
        assert!(bx.enter([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]).is_some());
        assert!(
            bx.enter([-1.0, 1.0 + 1e-15, 0.5], [3.0, 0.0, 0.0])
                .is_none()
        );
        // starting inside, ending short of the box
        assert_eq!(bx.enter([0.5; 3], [0.1, 0.0, 0.0]), Some(0.0));
        assert!(bx.enter([-2.0, 0.5, 0.5], [1.0, 0.0, 0.0]).is_none());
    }
}
