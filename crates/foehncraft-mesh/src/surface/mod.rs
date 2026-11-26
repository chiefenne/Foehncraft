//! Triangle surfaces of immersed bodies (L13 stage 0, D225): the surface
//! check, generalised winding numbers, inside/outside, signed distance and
//! the crossings of a segment with the surface.
//!
//! A [`Surface`] is a triangle soup, typically an STL file: the corners are
//! welded where their coordinates are **bitwise equal** (no tolerance, so
//! nothing is silently merged), which gives the edge structure the
//! [`check()`] needs. The geometry queries run on the bounding volume
//! hierarchy of [`crate::query`] (D108): the exact closest point
//! (distance), the exact segment–triangle test (crossings, D225 predicates)
//! and the hierarchical winding number of [`winding`].
//!
//! Orientation: an outward surface (corners counter-clockwise seen from
//! outside, the STL convention) has the winding number 1 inside and 0
//! outside. The signed distance is negative inside (φ < 0 in the body), the
//! convention of the cut-cell stage (D227). Dirty input (holes,
//! self-intersections, duplicate or flipped triangles) is tolerated: the
//! winding number degrades gracefully (Jacobson et al. 2013), the check
//! reports, and nothing is repaired (owner decision Q2, 2026-11-22).

pub mod body;
pub mod check;
pub mod distance;
pub mod moments;
pub mod winding;

pub use body::SurfaceBody;
pub use check::{SurfaceReport, check};
pub use distance::Crossing;
pub use moments::{CutMoments, WallPiece};
pub use winding::{Winding, solid_angle};

use crate::format::stl::Triangle;
use crate::geometry::{cross, norm, scale, sub};
use crate::query::{Bvh, SurfaceFaces};
use std::collections::HashMap;

/// A triangle surface with welded corners and its hierarchies.
#[derive(Debug, Clone)]
pub struct Surface {
    /// The distinct corner points (bitwise-equal corners welded), in order
    /// of first appearance.
    pub vertices: Vec<[f64; 3]>,
    /// Corner indices of every triangle, in input order.
    pub triangles: Vec<[u32; 3]>,
    pub(crate) bvh: Bvh,
    pub(crate) winding: Winding,
}

impl Surface {
    /// The surface of the triangles `tris` (in metres).
    pub fn new(tris: &[Triangle]) -> Self {
        let mut index: HashMap<[u64; 3], u32> = HashMap::new();
        let mut vertices = Vec::new();
        let mut triangles = Vec::with_capacity(tris.len());
        for t in tris {
            let mut ids = [0u32; 3];
            for (k, p) in t.iter().enumerate() {
                // -0.0 and 0.0 are the same point
                let key = p.map(|x| if x == 0.0 { 0u64 } else { x.to_bits() });
                ids[k] = *index.entry(key).or_insert_with(|| {
                    vertices.push(*p);
                    (vertices.len() - 1) as u32
                });
            }
            triangles.push(ids);
        }
        let mut faces = SurfaceFaces::default();
        for (i, t) in tris.iter().enumerate() {
            faces.tris.push(*t);
            faces.tri_face.push(i as u32);
            faces.global_id.push(i as u64);
            faces.mesh_face.push(i);
            let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
            let a = norm(n);
            faces
                .normal
                .push(if a > 0.0 { scale(n, 1.0 / a) } else { [0.0; 3] });
        }
        let bvh = Bvh::build(faces);
        let winding = Winding::build(&bvh);
        Surface {
            vertices,
            triangles,
            bvh,
            winding,
        }
    }

    /// Reads an STL file (ASCII or binary) and scales its coordinates by
    /// `scale` (for example 1e-3 for a file in millimetres).
    pub fn read_stl(path: &std::path::Path, scale: f64) -> Result<Self, String> {
        let mut tris = crate::format::stl::read(path)?;
        if scale != 1.0 {
            for t in &mut tris {
                for p in t.iter_mut() {
                    *p = p.map(|x| x * scale);
                }
            }
        }
        Ok(Surface::new(&tris))
    }

    /// Number of triangles.
    pub fn len(&self) -> usize {
        self.triangles.len()
    }

    /// No triangles.
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    /// The corners of triangle `t`.
    pub fn corners(&self, t: usize) -> Triangle {
        self.triangles[t].map(|v| self.vertices[v as usize])
    }

    /// The surface hierarchy (closest point, ray casting).
    pub fn bvh(&self) -> &Bvh {
        &self.bvh
    }

    /// Lower and upper corner of the bounding box (infinite for an empty
    /// surface).
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for v in &self.vertices {
            for k in 0..3 {
                lo[k] = lo[k].min(v[k]);
                hi[k] = hi[k].max(v[k]);
            }
        }
        (lo, hi)
    }
}
