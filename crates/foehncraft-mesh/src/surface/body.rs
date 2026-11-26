//! A triangle surface as the body of a cut (L13 stage 2, D227): the
//! [`ImplicitBody`] of a [`Surface`].
//!
//! - **Node values**: the signed distance (closest point, winding-number
//!   sign, negative inside), in a narrow band: a node farther from the
//!   surface than its longest edge gets ± its reach, only its sign counts.
//! - **Edge points**: the first exact crossing of the edge from its fluid
//!   end with the surface (watertight segment–triangle test, D225), the
//!   point where the edge enters the body; the zero of the linear
//!   interpolant (marching cubes) when [`SurfaceBody::exact`] is off or no
//!   crossing is found (a node within round-off of a dirty surface).
//! - **Wall point**: where a sharp edge or corner of the surface runs
//!   through the cell, the vertex of dual contouring (Ju et al. 2002) from
//!   the loop points and the normals of their triangles, kept only when it
//!   lies on the surface; elsewhere the exact closest point to the centroid
//!   of the wall loop.
//! - **Face points**: on every chord of a cut face, the crossing of the
//!   surface with the chord's perpendicular bisector in the face plane, or
//!   the point where a sharp edge pierces the face (the crossing of the two
//!   traces of its sides), again only when it lies on the surface.
//!
//! The normals are those of the closest triangle (outward, the STL
//! convention); everything is a function of the point alone, so the cut is
//! deterministic.

use super::Surface;
use crate::cut::{FEATURE_COS, ImplicitBody, chord_point, feature_point};
use crate::geometry::{add, norm, scale, sub};
use rayon::prelude::*;

/// A [`Surface`] as an [`ImplicitBody`] with its options.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceBody<'a> {
    /// The surface.
    pub surface: &'a Surface,
    /// Edge points at the exact crossings (else the linear interpolant of
    /// the node distances).
    pub exact: bool,
    /// Wall and face points on sharp edges and corners (dual contouring);
    /// else the closest point and the bisector only.
    pub features: bool,
}

impl<'a> SurfaceBody<'a> {
    /// The default body of a surface: exact crossings, sharp features.
    pub fn new(surface: &'a Surface) -> Self {
        Self {
            surface,
            exact: true,
            features: true,
        }
    }
}

impl Surface {
    /// The surface as the body of a cut with the default options.
    pub fn body(&self) -> SurfaceBody<'_> {
        SurfaceBody::new(self)
    }

    /// The outward unit normal of the triangle closest to `x`.
    pub fn normal_at(&self, x: [f64; 3]) -> Option<[f64; 3]> {
        self.bvh.closest(x).map(|c| c.normal)
    }

    /// The distance of `x` to the surface (unsigned).
    pub fn distance(&self, x: [f64; 3]) -> f64 {
        self.bvh.closest(x).map_or(f64::INFINITY, |c| c.distance)
    }

    /// The root of φ on the perpendicular bisector of the chord `p q` in
    /// the plane with normal `n` nearest to the chord's middle, from the
    /// exact crossings: the bisector segment of half-length |q − p| / 2 on
    /// each side.
    fn bisector_point(&self, p: [f64; 3], q: [f64; 3], n: [f64; 3]) -> Option<[f64; 3]> {
        let e = sub(q, p);
        let l = norm(e);
        if l == 0.0 {
            return None;
        }
        let d = crate::geometry::cross(n, scale(e, 1.0 / l));
        let dl = norm(d);
        if dl == 0.0 {
            return None;
        }
        let d = scale(d, 0.5 * l / dl);
        let c = scale(add(p, q), 0.5);
        let s0 = sub(c, d);
        let s1 = add(c, d);
        // the crossing nearest to the middle (t = ½)
        let x = self
            .crossings(s0, s1)
            .into_iter()
            .min_by(|a, b| (a.t - 0.5).abs().total_cmp(&(b.t - 0.5).abs()))?;
        // on the chord: a flat surface, nothing to bend
        ((x.t - 0.5).abs() > 1e-3).then_some(x.point)
    }
}

impl ImplicitBody for SurfaceBody<'_> {
    fn value(&self, x: [f64; 3]) -> f64 {
        self.surface.signed_distance(x)
    }

    fn values(&self, points: &[[f64; 3]], reach: &[f64]) -> Vec<f64> {
        points
            .par_iter()
            .zip(reach.par_iter())
            .map(|(&q, &r)| {
                if r > 0.0 && r.is_finite() {
                    self.surface.signed_distance_band(q, r)
                } else {
                    self.surface.signed_distance(q)
                }
            })
            .collect()
    }

    fn edge_point(&self, a: [f64; 3], b: [f64; 3], fa: f64, fb: f64) -> [f64; 3] {
        if self.exact
            && let Some(c) = self.surface.first_crossing(a, b)
        {
            return c.point;
        }
        let t = fa / (fa - fb);
        add(a, scale(sub(b, a), t))
    }

    fn project(&self, x: [f64; 3]) -> Option<[f64; 3]> {
        self.surface.bvh.closest(x).map(|c| c.point)
    }

    fn wall_point(&self, ring: &[[f64; 3]], centre: [f64; 3]) -> Option<[f64; 3]> {
        if self.features {
            let n: Vec<[f64; 3]> = ring
                .iter()
                .map(|&p| self.surface.normal_at(p).unwrap_or([0.0; 3]))
                .collect();
            if let Some(x) = feature_point(ring, &n, centre) {
                let r = ring
                    .iter()
                    .map(|&y| norm(sub(y, centre)))
                    .fold(0.0, f64::max);
                if self.surface.distance(x) <= 1e-7 * r {
                    return Some(x);
                }
            }
        }
        self.project(centre)
    }

    fn face_point(&self, p: [f64; 3], q: [f64; 3], n: [f64; 3]) -> Option<[f64; 3]> {
        if self.features {
            let s = self.surface;
            // a sharp edge through the face: the crossing of the traces of
            // its two sides; checked to lie on the surface by `chord_point`
            let (n1, n2) = (s.normal_at(p), s.normal_at(q));
            if let (Some(a), Some(b)) = (n1, n2)
                && crate::geometry::dot(a, b) < FEATURE_COS
                && let Some(x) = chord_point(|y| s.distance(y), |y| s.normal_at(y), p, q, n)
                && s.distance(x) <= 1e-7 * norm(sub(q, p))
            {
                return Some(x);
            }
        }
        self.surface.bisector_point(p, q, n)
    }
}

impl ImplicitBody for Surface {
    fn value(&self, x: [f64; 3]) -> f64 {
        self.body().value(x)
    }
    fn values(&self, points: &[[f64; 3]], reach: &[f64]) -> Vec<f64> {
        self.body().values(points, reach)
    }
    fn edge_point(&self, a: [f64; 3], b: [f64; 3], fa: f64, fb: f64) -> [f64; 3] {
        self.body().edge_point(a, b, fa, fb)
    }
    fn project(&self, x: [f64; 3]) -> Option<[f64; 3]> {
        self.body().project(x)
    }
    fn wall_point(&self, ring: &[[f64; 3]], centre: [f64; 3]) -> Option<[f64; 3]> {
        self.body().wall_point(ring, centre)
    }
    fn face_point(&self, p: [f64; 3], q: [f64; 3], n: [f64; 3]) -> Option<[f64; 3]> {
        self.body().face_point(p, q, n)
    }
}
