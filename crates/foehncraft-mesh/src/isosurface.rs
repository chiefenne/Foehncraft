//! Iso-surfaces of cell fields on any polyhedral mesh: marching tetrahedra
//! (L13 stage 1, D226).
//!
//! A cell field (a solid fraction or a design density) is first given
//! values at the mesh nodes, the mean of the cells around each node
//! weighted by the inverse distance to their centroids. Every cell is split
//! into tetrahedra (cell centroid, face centroid, one face edge), the cell
//! centroid carrying the cell value and the face centroid the mean of the
//! face's node values. Marching tetrahedra (Doi & Koide 1991) then puts at
//! most two triangles into each tetrahedron, with the crossing points
//! interpolated linearly along the tetrahedron edges. A crossing point
//! belongs to the edge, which two cells share when the edge is a mesh
//! edge or the edge from a face centroid: the surface is closed across
//! cells (watertight) and, since each tetrahedron holds at most one
//! disc that stays inside it, free of self-intersections by construction.
//! Triangles are oriented with the normal from the side of the larger
//! values to the side of the smaller ones (out of a solid).
//!
//! The output order is that of the cells, so the surface is bitwise
//! reproducible (D93).
//!
//! Benchmarked against dual contouring and subgrid marching tetrahedra
//! (Baktash, Gillespie & Crane 2026) in `tools/bench/isosurface.py`: on a
//! field with values at the nodes interpolated linearly along the edges
//! each edge holds at most one crossing, so subgrid marching gives the same
//! surface as this one; dual contouring keeps sharp corners but its
//! surfaces are not always manifold.

use crate::{Geometry, Mesh};
use std::collections::HashMap;
use std::fmt::Write as _;

/// A triangle surface.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Surface {
    /// Points (m).
    pub points: Vec<[f64; 3]>,
    /// Triangles (point indices), oriented out of the larger values.
    pub triangles: Vec<[u32; 3]>,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Edge counts of a surface: edges of one triangle (open boundary of the
/// surface), two (manifold) and more (non-manifold).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EdgeCounts {
    /// Edges with one triangle.
    pub boundary: usize,
    /// Edges with two triangles.
    pub manifold: usize,
    /// Edges with three or more.
    pub non_manifold: usize,
    /// Edges whose two triangles run along them in the same direction
    /// (inconsistent orientation).
    pub misoriented: usize,
}

impl Surface {
    /// Area (m²).
    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| self.points[i as usize]);
                0.5 * norm(cross(sub(b, a), sub(c, a)))
            })
            .sum()
    }

    /// The signed volume enclosed (divergence theorem; m³): the volume of
    /// the larger values for a closed surface.
    pub fn volume(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| self.points[i as usize]);
                dot(a, cross(b, c)) / 6.0
            })
            .sum::<f64>()
    }

    /// The edge counts.
    pub fn edges(&self) -> EdgeCounts {
        let mut m: HashMap<(u32, u32), (usize, i32)> = HashMap::new();
        for t in &self.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                let e = m.entry((a.min(b), a.max(b))).or_insert((0, 0));
                e.0 += 1;
                e.1 += if a < b { 1 } else { -1 };
            }
        }
        let mut out = EdgeCounts::default();
        for (n, s) in m.values() {
            match n {
                1 => out.boundary += 1,
                2 => {
                    out.manifold += 1;
                    if *s != 0 {
                        out.misoriented += 1;
                    }
                }
                _ => out.non_manifold += 1,
            }
        }
        out
    }

    /// The smallest angle of each triangle (degrees), in triangle order.
    pub fn min_angles(&self) -> Vec<f64> {
        self.triangles
            .iter()
            .map(|t| {
                let p = t.map(|i| self.points[i as usize]);
                (0..3)
                    .map(|k| {
                        let a = sub(p[(k + 1) % 3], p[k]);
                        let b = sub(p[(k + 2) % 3], p[k]);
                        let c = (dot(a, b) / (norm(a) * norm(b))).clamp(-1.0, 1.0);
                        c.acos().to_degrees()
                    })
                    .fold(180.0f64, |x, y| if y.is_nan() { x } else { x.min(y) })
            })
            .collect()
    }

    /// ASCII STL.
    pub fn stl(&self, name: &str) -> String {
        let mut s = format!("solid {name}\n");
        for t in &self.triangles {
            let [a, b, c] = t.map(|i| self.points[i as usize]);
            let n = cross(sub(b, a), sub(c, a));
            let l = norm(n);
            let n = if l > 0.0 { n.map(|x| x / l) } else { [0.0; 3] };
            let _ = writeln!(s, "  facet normal {:e} {:e} {:e}", n[0], n[1], n[2]);
            s.push_str("    outer loop\n");
            for p in [a, b, c] {
                let _ = writeln!(s, "      vertex {:e} {:e} {:e}", p[0], p[1], p[2]);
            }
            s.push_str("    endloop\n  endfacet\n");
        }
        let _ = writeln!(s, "endsolid {name}");
        s
    }
}

/// The node values of a cell field: the mean of the cells around each node
/// weighted by the inverse distance to the cell centroids.
pub fn node_values(mesh: &Mesh, geo: &Geometry, cell: &[f64]) -> Vec<f64> {
    let np = mesh.points.len();
    let mut sum = vec![0.0; np];
    let mut wsum = vec![0.0; np];
    for (c, &v) in cell.iter().enumerate().take(mesh.n_cells) {
        let xc = geo.cell_centroid[c];
        for n in mesh.cell_nodes(c) {
            let n = n as usize;
            let d = norm(sub(mesh.points[n], xc));
            let w = 1.0 / d.max(1e-300);
            sum[n] += w * v;
            wsum[n] += w;
        }
    }
    sum.iter()
        .zip(&wsum)
        .map(|(s, w)| if *w > 0.0 { s / w } else { 0.0 })
        .collect()
}

/// The surface `cell = level` of a cell field by marching tetrahedra.
pub fn extract(mesh: &Mesh, geo: &Geometry, cell: &[f64], level: f64) -> Surface {
    let nodes = node_values(mesh, geo, cell);
    extract_with(mesh, geo, cell, &nodes, level)
}

/// As [`extract`], with given node values.
pub fn extract_with(
    mesh: &Mesh,
    geo: &Geometry,
    cell: &[f64],
    nodes: &[f64],
    level: f64,
) -> Surface {
    let np = mesh.points.len() as u64;
    let nf = mesh.n_faces() as u64;
    // the value of every face centroid: the mean of its nodes (the same
    // from both cells)
    let face_value = |f: usize| -> f64 {
        let n = mesh.face(f);
        n.iter().map(|&i| nodes[i as usize]).sum::<f64>() / n.len() as f64
    };
    let mut out = Surface::default();
    let mut at: HashMap<(u64, u64), u32> = HashMap::new();
    // a corner of a tetrahedron: its global id (nodes, then face
    // centroids, then cell centroids), position and value
    type Corner = (u64, [f64; 3], f64);
    let mut point = |a: &Corner, b: &Corner, out: &mut Surface| -> u32 {
        let key = (a.0.min(b.0), a.0.max(b.0));
        *at.entry(key).or_insert_with(|| {
            // interpolate from the lower id, so both cells get the same bits
            let (p, q) = if a.0 < b.0 { (a, b) } else { (b, a) };
            let t = ((level - p.2) / (q.2 - p.2)).clamp(0.0, 1.0);
            let x = std::array::from_fn(|i| p.1[i] + t * (q.1[i] - p.1[i]));
            out.points.push(x);
            (out.points.len() - 1) as u32
        })
    };
    for (c, &value) in cell.iter().enumerate().take(mesh.n_cells) {
        let cc: Corner = (np + nf + c as u64, geo.cell_centroid[c], value);
        for (f, _) in mesh.cell_face_iter(c) {
            let fc: Corner = (np + f as u64, geo.face_centroid[f], face_value(f));
            let fnodes = mesh.face(f);
            for e in 0..fnodes.len() {
                let i = fnodes[e] as usize;
                let j = fnodes[(e + 1) % fnodes.len()] as usize;
                let tet: [Corner; 4] = [
                    cc,
                    fc,
                    (i as u64, mesh.points[i], nodes[i]),
                    (j as u64, mesh.points[j], nodes[j]),
                ];
                let inside: Vec<usize> = (0..4).filter(|&k| tet[k].2 >= level).collect();
                let outside: Vec<usize> = (0..4).filter(|&k| tet[k].2 < level).collect();
                let mean = |ks: &[usize]| -> [f64; 3] {
                    let mut m = [0.0; 3];
                    for &k in ks {
                        for (md, td) in m.iter_mut().zip(tet[k].1) {
                            *md += td / ks.len() as f64;
                        }
                    }
                    m
                };
                let polygon: Vec<u32> = match (inside.len(), outside.len()) {
                    (1, 3) => {
                        let a = inside[0];
                        outside
                            .iter()
                            .map(|&b| point(&tet[a], &tet[b], &mut out))
                            .collect()
                    }
                    (3, 1) => {
                        let b = outside[0];
                        inside
                            .iter()
                            .map(|&a| point(&tet[a], &tet[b], &mut out))
                            .collect()
                    }
                    (2, 2) => {
                        let (a, b) = (inside[0], inside[1]);
                        let (c2, d) = (outside[0], outside[1]);
                        vec![
                            point(&tet[a], &tet[c2], &mut out),
                            point(&tet[a], &tet[d], &mut out),
                            point(&tet[b], &tet[d], &mut out),
                            point(&tet[b], &tet[c2], &mut out),
                        ]
                    }
                    _ => continue,
                };
                // orientation: the normal from the inside to the outside
                let dir = sub(mean(&outside), mean(&inside));
                let p: Vec<[f64; 3]> = polygon.iter().map(|&k| out.points[k as usize]).collect();
                let mut n = [0.0; 3];
                for k in 1..p.len() - 1 {
                    let t = cross(sub(p[k], p[0]), sub(p[k + 1], p[0]));
                    for d in 0..3 {
                        n[d] += t[d];
                    }
                }
                let flip = dot(n, dir) < 0.0;
                for k in 1..polygon.len() - 1 {
                    let (b, c3) = (polygon[k], polygon[k + 1]);
                    let t = if flip {
                        [polygon[0], c3, b]
                    } else {
                        [polygon[0], b, c3]
                    };
                    // drop triangles collapsed onto a corner at the level
                    if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
                        out.triangles.push(t);
                    }
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate::BoxSpec;

    fn mesh(kind: &str, n: usize) -> (Mesh, Geometry) {
        let spec: BoxSpec = toml::from_str(&format!(
            "type = \"{kind}\"\nlower = [-1.0, -1.0, -1.0]\nupper = [1.0, 1.0, 1.0]\ncells = [{n}, {n}, {n}]"
        ))
        .expect("spec");
        let m = spec.generate().expect("mesh");
        let g = Geometry::compute(&m);
        (m, g)
    }

    /// The surface of the sphere r = 0.6 from the signed distance at the
    /// cell centroids: closed, manifold, consistently oriented, its area
    /// and volume converge to the sphere's at second order.
    #[test]
    fn sphere_from_a_distance_field_is_closed_and_converges() {
        let r = 0.6;
        let (sa, sv) = (
            4.0 * std::f64::consts::PI * r * r,
            4.0 / 3.0 * std::f64::consts::PI * r * r * r,
        );
        for kind in ["hex", "tet", "poly"] {
            let mut errs = Vec::new();
            for n in [8, 16] {
                let (m, g) = mesh(kind, n);
                // the distance inside the sphere, positive inside
                let field: Vec<f64> = g.cell_centroid.iter().map(|x| r - norm(*x)).collect();
                // exact node values: the test is the extraction, not the
                // node interpolation
                let nodes: Vec<f64> = m.points.iter().map(|x| r - norm(*x)).collect();
                let s = extract_with(&m, &g, &field, &nodes, 0.0);
                let e = s.edges();
                assert_eq!(e.boundary, 0, "{kind} {n}: open");
                assert_eq!(e.non_manifold, 0, "{kind} {n}: non-manifold");
                assert_eq!(e.misoriented, 0, "{kind} {n}: orientation");
                errs.push(
                    ((s.area() - sa) / sa)
                        .abs()
                        .max(((s.volume() - sv) / sv).abs()),
                );
                // every point on the sphere within the chord error
                let h = 2.0 / n as f64;
                let worst = s
                    .points
                    .iter()
                    .map(|p| (norm(*p) - r).abs())
                    .fold(0.0f64, f64::max);
                assert!(worst < h * h, "{kind} {n}: {worst}");
            }
            assert!(errs[1] < 0.4 * errs[0], "{kind}: {errs:?}");
        }
    }

    /// A smeared indicator (a cell solid fraction) with node values by
    /// inverse distance: still closed and manifold, the iso-surface ½
    /// within a cell of the sphere.
    #[test]
    fn sphere_from_a_cell_fraction_lies_within_a_cell() {
        let (m, g) = mesh("hex", 16);
        let h = 2.0 / 16.0;
        let r = 0.6;
        let field: Vec<f64> = g
            .cell_centroid
            .iter()
            .map(|x| (0.5 + (r - norm(*x)) / h).clamp(0.0, 1.0))
            .collect();
        let s = extract(&m, &g, &field, 0.5);
        let e = s.edges();
        assert_eq!((e.boundary, e.non_manifold, e.misoriented), (0, 0, 0));
        let worst = s
            .points
            .iter()
            .map(|p| (norm(*p) - r).abs())
            .fold(0.0f64, f64::max);
        assert!(worst < h, "{worst}");
        // reproducible
        assert_eq!(s, extract(&m, &g, &field, 0.5));
    }
}
