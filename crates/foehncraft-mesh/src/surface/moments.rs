//! The geometry of mesh cells cut by a triangle surface (L13 stage 0,
//! D225 item 5): for every cell the volume and first moment of its part
//! inside the body, for every face the area and first moment of its part
//! inside, and the wall pieces (the surface clipped to each cell) with
//! their area vector, centroid and surface triangle. These are the
//! apertures, volume fractions and wall faces a cut-cell finite-volume
//! method needs (stage 2, D227), on any polyhedral mesh.
//!
//! **Method: Gauss–Green by dimension.** Nothing is built as a polygon or
//! polyhedron; every quantity is a boundary integral of pieces that are
//! straight or planar:
//!
//! 1. *Edges.* The part of a segment inside the body comes from its exact
//!    crossings with the surface (orientation predicates of D225, D108 ray
//!    test): after each crossing the segment is inside when it entered.
//!    Only an edge without any crossing asks the winding number.
//! 2. *Faces.* A mesh face is split into the triangles of its geometry (one
//!    triangle, or the fan about the node average, the surfaces whose
//!    areas the solver uses). For a planar triangle `F` with centroid
//!    `c`, the divergence theorem in its plane gives
//!    `A = ½ ∮ (x − c)·ν ds` and `∫ (x − c) dA = ⅓ ∮ (x − c)(x − c)·ν ds`
//!    over the boundary of `F ∩ body`: the inside parts of the three edges
//!    (outward normal of `F`) and the segments `t ∩ F` of the surface
//!    triangles (in-plane normal: the projection of the triangle's outward
//!    normal). On each straight piece `(x − c)·ν` is constant, so each
//!    piece adds its length and midpoint only.
//! 3. *Cells.* With the cell's apex `a` (the same interior point as the
//!    tetrahedral split of the solver's geometry), `V = ⅓ ∮ (x − a)·n dA`
//!    and `∫ (x − a) dV = ¼ ∮ (x − a)(x − a)·n dA` over `cell ∩ body`: the
//!    inside parts of the face triangles from step 2, and the surface
//!    clipped to the cell. The surface is clipped to each tetrahedron
//!    (apex, face triangle) by the four planes (Sutherland–Hodgman), the
//!    corners of the triangle classified by the exact orientation
//!    predicate. The planes through the apex carry no flux (their
//!    `(x − a)·n` is zero), so the split adds nothing but the clipping.
//!
//! The pieces of a face are computed once per face, so two cells see the
//! same apertures (the cut is conservative). For a closed surface the
//! inside volumes add up to the volume the surface encloses inside the
//! mesh, and the wall pieces to its area, up to round-off (tests).
//!
//! Not handled yet: a surface lying exactly in a mesh face or through
//! exactly a mesh node (simulation of simplicity, Edelsbrunner & Mücke
//! 1990, would decide those ties consistently); the polygon pieces with
//! implicit intersection points (Attene 2020) that a mesh-building cut
//! needs. The rounded intersection points only move the moments by
//! round-off.

use super::Surface;
use crate::Mesh;
use crate::geometry::{Geometry, add, cell_apex, cross, dot, node_average, norm, scale, sub};
use foehncraft_core::predicates::orient3d;
use rayon::prelude::*;

/// A piece of the surface inside one cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallPiece {
    /// The cell.
    pub cell: usize,
    /// The surface triangle (input order).
    pub triangle: usize,
    /// Area vector, out of the body (into the fluid).
    pub area: [f64; 3],
    /// Centroid.
    pub centroid: [f64; 3],
}

/// The cut of a mesh by a surface.
#[derive(Debug, Clone, Default)]
pub struct CutMoments {
    /// Volume of each cell inside the body.
    pub cell_inside_volume: Vec<f64>,
    /// Centroid of the inside part of each cell (the cell centroid when
    /// nothing is inside).
    pub cell_inside_centroid: Vec<[f64; 3]>,
    /// Area of each face inside the body.
    pub face_inside_area: Vec<f64>,
    /// Centroid of the inside part of each face.
    pub face_inside_centroid: Vec<[f64; 3]>,
    /// The surface pieces, ordered by cell, then triangle.
    pub walls: Vec<WallPiece>,
}

/// Area and first moment `∫ x dA` of the inside part of a planar triangle.
#[derive(Debug, Clone, Copy, Default)]
struct Wet {
    area: f64,
    moment: [f64; 3],
}

impl Surface {
    /// The surface triangles whose boxes overlap the box `lo`–`hi`, in
    /// increasing order.
    pub fn triangles_in_box(&self, lo: [f64; 3], hi: [f64; 3]) -> Vec<usize> {
        let bvh = &self.bvh;
        let mut out = Vec::new();
        if bvh.nodes.is_empty() {
            return out;
        }
        let mut stack = vec![0u32];
        while let Some(i) = stack.pop() {
            let nd = bvh.nodes[i as usize];
            if (0..3).any(|k| nd.bx.hi[k] < lo[k] || nd.bx.lo[k] > hi[k]) {
                continue;
            }
            if nd.count > 0 {
                let s = nd.first as usize;
                for &t in &bvh.order[s..s + nd.count as usize] {
                    let tri = &bvh.faces.tris[t as usize];
                    if (0..3).all(|k| {
                        tri.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min) <= hi[k]
                            && tri.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max) >= lo[k]
                    }) {
                        out.push(t as usize);
                    }
                }
            } else {
                stack.push(nd.first);
                stack.push(nd.first + 1);
            }
        }
        out.sort_unstable();
        out
    }

    /// The inside intervals of the segment `p q` as parameter pairs.
    fn inside_intervals(&self, p: [f64; 3], q: [f64; 3]) -> Vec<(f64, f64)> {
        let cross = self.crossings(p, q);
        let mut state = match cross.first() {
            Some(c) => !c.entering,
            None => self.inside(scale(add(p, q), 0.5)),
        };
        let mut out = Vec::new();
        let mut start = 0.0;
        for c in &cross {
            if state && !c.entering {
                out.push((start, c.t));
            }
            if !state && c.entering {
                start = c.t;
            }
            state = c.entering;
        }
        if state {
            out.push((start, 1.0));
        }
        out
    }

    /// Area and first moment of the planar triangle `f` inside the body,
    /// with the surface triangles `cand` as the candidates.
    fn wet_triangle(&self, f: [[f64; 3]; 3], cand: &[usize]) -> Wet {
        let nf = cross(sub(f[1], f[0]), sub(f[2], f[0]));
        let len = norm(nf);
        if len == 0.0 {
            return Wet::default();
        }
        let nf = scale(nf, 1.0 / len);
        let c = scale(add(add(f[0], f[1]), f[2]), 1.0 / 3.0);
        // Σ h L and Σ h L (mid − c) over the boundary pieces
        let mut a2 = 0.0;
        let mut m3 = [0.0; 3];
        let mut piece = |x0: [f64; 3], x1: [f64; 3], nu: [f64; 3]| {
            let l = norm(sub(x1, x0));
            let mid = scale(add(x0, x1), 0.5);
            let h = dot(sub(mid, c), nu);
            a2 += h * l;
            m3 = add(m3, scale(sub(mid, c), h * l));
        };
        for k in 0..3 {
            let (p, q) = (f[k], f[(k + 1) % 3]);
            let nu = cross(sub(q, p), nf);
            let nl = norm(nu);
            if nl == 0.0 {
                continue;
            }
            let nu = scale(nu, 1.0 / nl);
            for (t0, t1) in self.inside_intervals(p, q) {
                let x0 = add(p, scale(sub(q, p), t0));
                let x1 = add(p, scale(sub(q, p), t1));
                piece(x0, x1, nu);
            }
        }
        for &t in cand {
            let tri = self.corners(t);
            let nt = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]));
            let nu = sub(nt, scale(nf, dot(nt, nf)));
            let nl = norm(nu);
            if nl <= 1e-12 * norm(nt) {
                continue; // coplanar: no segment
            }
            if !straddles(&tri, &f) {
                continue;
            }
            if let Some((x0, x1)) = tri_tri_segment(&tri, &f) {
                piece(x0, x1, scale(nu, 1.0 / nl));
            }
        }
        let area = 0.5 * a2;
        Wet {
            area,
            moment: add(scale(c, area), scale(m3, 1.0 / 3.0)),
        }
    }

    /// The cut of every cell and face of `mesh` (geometry `geo`) by the
    /// surface, in parallel; the same result on any number of threads.
    pub fn cut_moments(&self, mesh: &Mesh, geo: &Geometry) -> CutMoments {
        let nf = mesh.n_faces();
        // face triangles: (face, corners)
        let face_tris = |f: usize| -> Vec<[[f64; 3]; 3]> {
            let (nodes, _) = mesh.face_global_order(f);
            let p = |i: usize| mesh.points[nodes[i] as usize];
            if nodes.len() == 3 {
                vec![[p(0), p(1), p(2)]]
            } else {
                let m = node_average(&mesh.points, &nodes);
                (0..nodes.len())
                    .map(|i| [m, p(i), p((i + 1) % nodes.len())])
                    .collect()
            }
        };
        let bbox = |pts: &[[f64; 3]]| {
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for p in pts {
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
            (lo, hi)
        };
        // per face triangle: wet area and moment (computed once per face)
        let wet: Vec<Vec<Wet>> = (0..nf)
            .into_par_iter()
            .map(|f| {
                let tris = face_tris(f);
                tris.iter()
                    .map(|t| {
                        let (lo, hi) = bbox(t);
                        let cand = self.triangles_in_box(lo, hi);
                        if cand.is_empty() {
                            // the whole triangle on one side
                            let c = scale(add(add(t[0], t[1]), t[2]), 1.0 / 3.0);
                            let a = 0.5 * norm(cross(sub(t[1], t[0]), sub(t[2], t[0])));
                            if self.inside(c) {
                                Wet {
                                    area: a,
                                    moment: scale(c, a),
                                }
                            } else {
                                Wet::default()
                            }
                        } else {
                            self.wet_triangle(*t, &cand)
                        }
                    })
                    .collect()
            })
            .collect();
        let mut face_inside_area = vec![0.0; nf];
        let mut face_inside_centroid = geo.face_centroid.clone();
        for f in 0..nf {
            let a: f64 = wet[f].iter().map(|w| w.area).sum();
            let m = wet[f].iter().fold([0.0; 3], |s, w| add(s, w.moment));
            face_inside_area[f] = a;
            if a > 0.0 {
                face_inside_centroid[f] = scale(m, 1.0 / a);
            }
        }
        // cells
        let per_cell: Vec<(f64, [f64; 3], Vec<WallPiece>)> = (0..mesh.n_cells)
            .into_par_iter()
            .map(|c| {
                let apex = cell_apex(mesh, &geo.face_centroid, c);
                let mut pts = Vec::new();
                let mut tris: Vec<([[f64; 3]; 3], Wet)> = Vec::new();
                for (f, _) in mesh.cell_face_iter(c) {
                    for (t, w) in face_tris(f).into_iter().zip(&wet[f]) {
                        pts.extend_from_slice(&t);
                        tris.push((t, *w));
                    }
                }
                let (lo, hi) = bbox(&pts);
                let cand = self.triangles_in_box(lo, hi);
                let mut vol3 = 0.0;
                let mut mom4 = [0.0; 3];
                // faces: h A and h ∫ (x − a) dA
                for (t, w) in &tris {
                    if w.area == 0.0 {
                        continue;
                    }
                    let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
                    let nl = norm(n);
                    if nl == 0.0 {
                        continue;
                    }
                    let h = (dot(sub(t[0], apex), n) / nl).abs();
                    vol3 += h * w.area;
                    mom4 = add(mom4, scale(sub(w.moment, scale(apex, w.area)), h));
                }
                // the surface clipped to each tetrahedron
                let mut walls: Vec<WallPiece> = Vec::new();
                if !cand.is_empty() {
                    for &s in &cand {
                        let tri = self.corners(s);
                        let mut area = [0.0; 3];
                        let mut first = [0.0; 3];
                        for (t, _) in &tris {
                            let tet = [apex, t[0], t[1], t[2]];
                            if let Some((av, cen)) = clip_to_tet(&tri, &tet) {
                                area = add(area, av);
                                first = add(first, scale(cen, norm(av)));
                                let h = dot(sub(cen, apex), av);
                                vol3 += h;
                                mom4 = add(mom4, scale(sub(cen, apex), h));
                            }
                        }
                        let a = norm(area);
                        if a > 0.0 {
                            walls.push(WallPiece {
                                cell: c,
                                triangle: s,
                                area,
                                centroid: scale(first, 1.0 / a),
                            });
                        }
                    }
                }
                let v = vol3 / 3.0;
                let cen = if v > 0.0 {
                    add(apex, scale(mom4, 1.0 / (4.0 * v)))
                } else {
                    geo.cell_centroid[c]
                };
                (v, cen, walls)
            })
            .collect();
        let mut out = CutMoments {
            face_inside_area,
            face_inside_centroid,
            ..Default::default()
        };
        for (v, cen, w) in per_cell {
            out.cell_inside_volume.push(v);
            out.cell_inside_centroid.push(cen);
            out.walls.extend(w);
        }
        out
    }
}

/// Whether the surface triangle `t` has corners on both sides of the plane
/// of `f`, with the corners in the plane decided by simulation of
/// simplicity (Edelsbrunner & Mücke 1990): the surface is moved by the
/// infinitesimal `δ = (ε, ε², ε³)`, so a corner in the plane takes the
/// sign of `−N·δ`, the first nonzero component of the plane normal `N`
/// (exact, from 2D orientations) with the sign reversed. Two surface
/// triangles that share an edge lying in a mesh face then give that face
/// one segment, not two.
fn straddles(t: &[[f64; 3]; 3], f: &[[f64; 3]; 3]) -> bool {
    use foehncraft_core::predicates::orient2d;
    let yz = |p: [f64; 3]| [p[1], p[2]];
    let zx = |p: [f64; 3]| [p[2], p[0]];
    let xy = |p: [f64; 3]| [p[0], p[1]];
    let tie = || {
        for proj in [yz, zx, xy] {
            let n = orient2d(proj(f[0]), proj(f[1]), proj(f[2]));
            if n != 0.0 {
                return if n > 0.0 { -1.0 } else { 1.0 };
            }
        }
        0.0
    };
    let mut pos = false;
    let mut neg = false;
    for &p in t {
        let mut o = orient3d(f[0], f[1], f[2], p);
        if o == 0.0 {
            o = tie();
        }
        pos |= o > 0.0;
        neg |= o < 0.0;
    }
    pos && neg
}

/// The segment where two triangles meet (not coplanar), or `None`: the
/// crossings of each triangle's edges with the other, the two farthest
/// apart.
fn tri_tri_segment(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> Option<([f64; 3], [f64; 3])> {
    use crate::query::ray::{End, segment_triangle};
    let mut pts: Vec<[f64; 3]> = Vec::with_capacity(6);
    let mut add_hits = |s: &[[f64; 3]; 3], t: &[[f64; 3]; 3]| {
        for k in 0..3 {
            let (p, q) = (s[k], s[(k + 1) % 3]);
            // closed segment: both directions cover an end point in the plane
            let hit = segment_triangle(p, q, t, End::Closed)
                .map(|h| add(p, scale(sub(q, p), h.t)))
                .or_else(|| {
                    segment_triangle(q, p, t, End::Closed).map(|h| add(q, scale(sub(p, q), h.t)))
                });
            if let Some(x) = hit {
                pts.push(x);
            }
        }
    };
    add_hits(a, b);
    add_hits(b, a);
    if pts.len() < 2 {
        return None;
    }
    let mut best = (0, 0, 0.0);
    for i in 0..pts.len() {
        for j in i + 1..pts.len() {
            let d = norm(sub(pts[i], pts[j]));
            if d > best.2 {
                best = (i, j, d);
            }
        }
    }
    if best.2 == 0.0 {
        return None;
    }
    Some((pts[best.0], pts[best.1]))
}

/// The triangle `tri` clipped to the tetrahedron `tet` (Sutherland–Hodgman
/// by its four planes): area vector (the triangle's orientation) and
/// centroid of the piece, or `None` when nothing is left.
fn clip_to_tet(tri: &[[f64; 3]; 3], tet: &[[f64; 3]; 4]) -> Option<([f64; 3], [f64; 3])> {
    // reject by box
    for k in 0..3 {
        let tmin = tri.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min);
        let tmax = tri.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max);
        let kmin = tet.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min);
        let kmax = tet.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max);
        if tmax < kmin || tmin > kmax {
            return None;
        }
    }
    let mut poly: Vec<[f64; 3]> = tri.to_vec();
    for i in 0..4 {
        let (p, q, r, s) = (tet[(i + 1) % 4], tet[(i + 2) % 4], tet[(i + 3) % 4], tet[i]);
        let side = orient3d(p, q, r, s);
        if side == 0.0 {
            return None; // flat tetrahedron
        }
        let val = |x: [f64; 3]| {
            let o = orient3d(p, q, r, x);
            if side > 0.0 { o } else { -o }
        };
        let vals: Vec<f64> = poly.iter().map(|&x| val(x)).collect();
        if vals.iter().all(|&v| v >= 0.0) {
            continue;
        }
        if vals.iter().all(|&v| v <= 0.0) {
            return None;
        }
        let mut out = Vec::with_capacity(poly.len() + 1);
        for k in 0..poly.len() {
            let (a, b) = (poly[k], poly[(k + 1) % poly.len()]);
            let (va, vb) = (vals[k], vals[(k + 1) % poly.len()]);
            if va >= 0.0 {
                out.push(a);
            }
            if (va > 0.0 && vb < 0.0) || (va < 0.0 && vb > 0.0) {
                let w = va / (va - vb);
                out.push(add(a, scale(sub(b, a), w)));
            }
        }
        if out.len() < 3 {
            return None;
        }
        poly = out;
    }
    // area vector and centroid by the fan about the first corner
    let mut av = [0.0; 3];
    let mut first = [0.0; 3];
    let mut total = 0.0;
    let n_tri = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]));
    for k in 1..poly.len() - 1 {
        let s = scale(cross(sub(poly[k], poly[0]), sub(poly[k + 1], poly[0])), 0.5);
        let a = dot(s, n_tri).signum() * norm(s);
        av = add(av, s);
        let c = scale(add(add(poly[0], poly[k]), poly[k + 1]), 1.0 / 3.0);
        first = add(first, scale(c, a));
        total += a;
    }
    if total <= 0.0 {
        return None;
    }
    Some((av, scale(first, 1.0 / total)))
}
