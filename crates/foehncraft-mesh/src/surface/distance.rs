//! Signed distance to a triangle surface and the crossings of a segment
//! with it (L13 stage 0, D225).
//!
//! The distance is the exact closest-point distance of the surface
//! hierarchy (D108 part 1, Ericson's Voronoi-region test); the sign comes
//! from the generalised winding number ([`super::winding`]), negative
//! inside. The sign is therefore defined for dirty surfaces too (holes,
//! overlaps), where a pseudo-normal sign is not.
//!
//! The narrow band version computes the distance only up to a band width
//! and returns ±band beyond it (still signed), which is what a cut or a
//! level set on a background mesh needs and costs O(log n) per point near
//! the surface, less far from it.
//!
//! The crossings of a segment are the exact segment–triangle tests of
//! [`crate::query::ray`] (watertight, orientation predicates of D225);
//! the root of φ along a mesh edge is the first crossing, the exact place
//! where the edge leaves or enters the body, not a linear interpolation of
//! two distances.

use super::Surface;
use rayon::prelude::*;

/// Where a segment crosses the surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Crossing {
    /// Segment parameter in (0, 1]: the point is about `p + t (q − p)`.
    pub t: f64,
    /// The crossing point, on the triangle.
    pub point: [f64; 3],
    /// The triangle (input order).
    pub triangle: usize,
    /// The segment goes from outside into the body (against the outward
    /// normal of the triangle).
    pub entering: bool,
}

impl Surface {
    /// The signed distance at `q` (m): negative inside the body.
    pub fn signed_distance(&self, q: [f64; 3]) -> f64 {
        match self.bvh.closest(q) {
            None => f64::INFINITY,
            Some(c) => {
                if self.inside(q) {
                    -c.distance
                } else {
                    c.distance
                }
            }
        }
    }

    /// The signed distance limited to the band `|φ| ≤ band`: points
    /// farther from the surface get ±band with the sign of their side.
    pub fn signed_distance_band(&self, q: [f64; 3], band: f64) -> f64 {
        let d = self
            .bvh
            .closest_within(q, band)
            .map_or(band, |c| c.distance.min(band));
        if self.inside(q) { -d } else { d }
    }

    /// The signed distances of many points (parallel, deterministic); with
    /// `band = Some(b)` limited to the band.
    pub fn signed_distances(&self, points: &[[f64; 3]], band: Option<f64>) -> Vec<f64> {
        points
            .par_iter()
            .map(|&q| match band {
                Some(b) => self.signed_distance_band(q, b),
                None => self.signed_distance(q),
            })
            .collect()
    }

    /// Every crossing of the segment from `p` to `q` with the surface,
    /// ordered along it (`p` excluded, `q` included). A segment through an
    /// edge or a corner of the surface counts once.
    pub fn crossings(&self, p: [f64; 3], q: [f64; 3]) -> Vec<Crossing> {
        self.bvh
            .all_hits(p, q)
            .into_iter()
            .map(|h| Crossing {
                t: h.t,
                point: h.point,
                triangle: h.face,
                entering: h.entering,
            })
            .collect()
    }

    /// The first crossing of the segment from `p` to `q`, the root of the
    /// signed distance along a mesh edge.
    pub fn first_crossing(&self, p: [f64; 3], q: [f64; 3]) -> Option<Crossing> {
        self.bvh.first_hit(p, q).map(|h| Crossing {
            t: h.t,
            point: h.point,
            triangle: h.face,
            entering: h.entering,
        })
    }
}
