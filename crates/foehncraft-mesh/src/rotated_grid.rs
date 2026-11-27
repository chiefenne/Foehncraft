//! Box meshes whose cell lines are turned in the x–y plane (`grid_angle`,
//! D337): a square-cornered grid of spacings `(x₁ − x₀)/n_x` along the
//! direction at the angle α from x and `(y₁ − y₀)/n_y` across it, cut by
//! the four sides of the box. The cells the sides cut become convex
//! polygons; a cut piece smaller than [`SMALL`] of a whole cell joins the
//! neighbour with which it shares the longest edge (the pieces of a joined
//! cell keep their faces, so two cells can share more than one face, as
//! the merged cut cells of L13, D228). The polygons are extruded over the
//! z intervals of the box (one layer in 2D), or, with `wedge`, revolved
//! about the x axis into the sector of that angle as [`crate::wedge`] does
//! for a box of hexahedra: the node (x, y, z) goes to (x, y cos φ, y sin φ)
//! with φ from −θ/2 at `zmin` to +θ/2 at `zmax`, every face planar.
//!
//! The point is a grid whose lines are not those of the box: a jet along
//! the direction α runs along the cell lines, one at α ± 45° runs across
//! them, on otherwise the same cells (the grid orientation test of Spray
//! G, gallery 337).

use crate::builder::MeshBuilder;
use crate::mesh::Mesh;
use std::collections::HashMap;

/// A cut piece of less than this fraction of a whole cell joins a
/// neighbour.
pub const SMALL: f64 = 0.3;

/// The mesh of the box `lower`–`upper` with `n` intervals (x, y and z),
/// the cell lines turned by `angle_deg` in the x–y plane; `dimension` 2
/// gives one layer in z; `wedge` revolves into a sector (y is the radius,
/// `lower[1] > 0`).
pub fn generate(
    lower: [f64; 3],
    upper: [f64; 3],
    n: [usize; 3],
    angle_deg: f64,
    wedge: Option<f64>,
) -> Result<Mesh, String> {
    let (x0, y0, x1, y1) = (lower[0], lower[1], upper[0], upper[1]);
    let hu = (x1 - x0) / n[0] as f64;
    let hv = (y1 - y0) / n[1] as f64;
    let a = angle_deg.to_radians();
    let (eu, ev) = ([a.cos(), a.sin()], [-a.sin(), a.cos()]);
    // grid nodes p(i, j) = origin + i hu eu + j hv ev, origin at the box
    // centre; the index ranges cover the box
    let (cx, cy) = (0.5 * (x0 + x1), 0.5 * (y0 + y1));
    let half = 0.5 * ((x1 - x0).hypot(y1 - y0));
    let iu = (half / hu).ceil() as i64 + 1;
    let iv = (half / hv).ceil() as i64 + 1;
    let node = |i: i64, j: i64| -> [f64; 2] {
        [
            cx + i as f64 * hu * eu[0] + j as f64 * hv * ev[0],
            cy + i as f64 * hu * eu[1] + j as f64 * hv * ev[1],
        ]
    };
    let full = hu * hv;
    let len = (x1 - x0).max(y1 - y0);
    let tol = 1e-10 * len;

    // the polygons of the cut cells, counter-clockwise, on shared vertices
    let mut verts: Vec<[f64; 2]> = Vec::new();
    let mut index: HashMap<(i64, i64), u32> = HashMap::new();
    let q = |v: f64| (v / tol).round() as i64;
    let mut vertex = |p: [f64; 2], verts: &mut Vec<[f64; 2]>| -> u32 {
        *index.entry((q(p[0]), q(p[1]))).or_insert_with(|| {
            verts.push(p);
            (verts.len() - 1) as u32
        })
    };
    let mut polys: Vec<Vec<u32>> = Vec::new();
    for i in -iu..iu {
        for j in -iv..iv {
            let sq = [node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)];
            let cut = clip(&sq, x0, y0, x1, y1);
            if cut.len() < 3 || area(&cut) < 1e-9 * full {
                continue;
            }
            let mut p: Vec<u32> = Vec::with_capacity(cut.len());
            for v in cut {
                let k = vertex(v, &mut verts);
                if p.last() != Some(&k) && p.first() != Some(&k) {
                    p.push(k);
                }
            }
            if p.len() >= 3 {
                polys.push(p);
            }
        }
    }
    let np = polys.len();
    let areas: Vec<f64> = polys
        .iter()
        .map(|p| area(&p.iter().map(|&k| verts[k as usize]).collect::<Vec<_>>()))
        .collect();

    // the edges: (sorted vertex pair) -> polygons with the orientation
    let mut edges: HashMap<(u32, u32), Vec<(u32, bool)>> = HashMap::new();
    for (c, p) in polys.iter().enumerate() {
        for k in 0..p.len() {
            let (a, b) = (p[k], p[(k + 1) % p.len()]);
            edges
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push((c as u32, a < b));
        }
    }
    for (e, cs) in &edges {
        if cs.len() > 2 {
            return Err(format!("rotated grid: edge {e:?} has {} cells", cs.len()));
        }
    }

    // the small pieces join their neighbour along the longest edge
    let mut group: Vec<usize> = (0..np).collect();
    fn root(g: &mut [usize], mut c: usize) -> usize {
        while g[c] != c {
            g[c] = g[g[c]];
            c = g[c];
        }
        c
    }
    let mut gsize = areas.clone();
    let mut order: Vec<usize> = (0..np).filter(|&c| areas[c] < SMALL * full).collect();
    order.sort_by(|&a, &b| areas[a].total_cmp(&areas[b]));
    for c in order {
        let rc = root(&mut group, c);
        if gsize[rc] >= SMALL * full {
            continue;
        }
        let mut shared: HashMap<usize, f64> = HashMap::new();
        let p = &polys[c];
        for k in 0..p.len() {
            let (a, b) = (p[k], p[(k + 1) % p.len()]);
            for &(o, _) in &edges[&(a.min(b), a.max(b))] {
                let ro = root(&mut group, o as usize);
                if ro != rc {
                    let (pa, pb) = (verts[a as usize], verts[b as usize]);
                    *shared.entry(ro).or_default() += (pb[0] - pa[0]).hypot(pb[1] - pa[1]);
                }
            }
        }
        if let Some((&to, _)) = shared
            .iter()
            .max_by(|x, y| x.1.total_cmp(y.1).then(y.0.cmp(x.0)))
        {
            group[rc] = to;
            gsize[to] += gsize[rc];
        }
    }
    let mut cell_of_group: HashMap<usize, u32> = HashMap::new();
    let mut cell2d = vec![0u32; np];
    for c in 0..np {
        let r = root(&mut group, c);
        let k = cell_of_group.len() as u32;
        cell2d[c] = *cell_of_group.entry(r).or_insert(k);
    }
    let nc2 = cell_of_group.len();

    // the nodes of the layers and the 3D cells (layer-major)
    let nz = n[2];
    let nv = verts.len();
    let mut pts: Vec<[f64; 3]> = Vec::with_capacity(nv * (nz + 1));
    for k in 0..=nz {
        let t = k as f64 / nz as f64;
        for v in &verts {
            pts.push(match wedge {
                Some(th) => {
                    let th = th.to_radians();
                    let phi = th * t - 0.5 * th;
                    [v[0], v[1] * phi.cos(), v[1] * phi.sin()]
                }
                None => [v[0], v[1], lower[2] + t * (upper[2] - lower[2])],
            });
        }
    }
    let at = |v: u32, k: usize| v + (k * nv) as u32;
    let cell = |c2: u32, k: usize| c2 + (k * nc2) as u32;
    let mut b = MeshBuilder::new(pts, nc2 * nz);
    b.allow_shared_pairs();
    for name in crate::generate::BOX_BOUNDARIES {
        b.declare_boundary(name);
    }
    // faces of constant z (the normal of a counter-clockwise polygon is +z)
    for (c, p) in polys.iter().enumerate() {
        let c2 = cell2d[c];
        for k in 0..=nz {
            let f: Vec<u32> = p.iter().map(|&v| at(v, k)).collect();
            if k == 0 {
                let r: Vec<u32> = f.iter().rev().copied().collect();
                b.add_boundary(&r, cell(c2, 0), "zmin");
            } else if k == nz {
                b.add_boundary(&f, cell(c2, nz - 1), "zmax");
            } else {
                b.add_interior(&f, cell(c2, k - 1), cell(c2, k));
            }
        }
    }
    // side faces, in a fixed order
    let mut keys: Vec<&(u32, u32)> = edges.keys().collect();
    keys.sort_unstable();
    for key in keys {
        let cs = &edges[key];
        // the polygon whose counter-clockwise loop runs a -> b
        let (own, fwd) = cs[0];
        let (a, bb) = if fwd { (key.0, key.1) } else { (key.1, key.0) };
        let side = match cs.get(1) {
            Some(&(o, _)) => {
                if cell2d[o as usize] == cell2d[own as usize] {
                    continue;
                }
                Some(cell2d[o as usize])
            }
            None => None,
        };
        let name = match side {
            Some(_) => "",
            None => {
                let (pa, pb) = (verts[a as usize], verts[bb as usize]);
                let m = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
                let d = [
                    (m[0] - x0).abs(),
                    (m[0] - x1).abs(),
                    (m[1] - y0).abs(),
                    (m[1] - y1).abs(),
                ];
                let i = (0..4).min_by(|&i, &j| d[i].total_cmp(&d[j])).unwrap();
                if d[i] > 1e-6 * len {
                    return Err(format!("rotated grid: an open edge inside the box at {m:?}"));
                }
                ["xmin", "xmax", "ymin", "ymax"][i]
            }
        };
        for k in 0..nz {
            let f = [at(a, k), at(bb, k), at(bb, k + 1), at(a, k + 1)];
            match side {
                Some(o) => b.add_interior(&f, cell(cell2d[own as usize], k), cell(o, k)),
                None => b.add_boundary(&f, cell(cell2d[own as usize], k), name),
            }
        }
    }
    b.finish()
}

/// The polygon `p` cut by the box (Sutherland and Hodgman); intersections
/// are computed from the sorted end points, so that both cells of a grid
/// edge get the same point.
fn clip(p: &[[f64; 2]], x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<[f64; 2]> {
    let mut out: Vec<[f64; 2]> = p.to_vec();
    // (axis, value, keep the side above)
    for (d, c, above) in [(0usize, x0, true), (0, x1, false), (1, y0, true), (1, y1, false)] {
        let inside = |v: &[f64; 2]| if above { v[d] >= c } else { v[d] <= c };
        let input = std::mem::take(&mut out);
        if input.is_empty() {
            break;
        }
        for k in 0..input.len() {
            let (s, e) = (input[(k + input.len() - 1) % input.len()], input[k]);
            let cross = |s: [f64; 2], e: [f64; 2]| {
                let (lo, hi) = if (s[0], s[1]) <= (e[0], e[1]) { (s, e) } else { (e, s) };
                let t = (c - lo[d]) / (hi[d] - lo[d]);
                let mut r = [lo[0] + t * (hi[0] - lo[0]), lo[1] + t * (hi[1] - lo[1])];
                r[d] = c;
                r
            };
            match (inside(&s), inside(&e)) {
                (true, true) => out.push(e),
                (true, false) => out.push(cross(s, e)),
                (false, true) => {
                    out.push(cross(s, e));
                    out.push(e);
                }
                (false, false) => {}
            }
        }
    }
    out
}

/// Signed area (counter-clockwise positive).
fn area(p: &[[f64; 2]]) -> f64 {
    let mut s = 0.0;
    for k in 0..p.len() {
        let (a, b) = (p[k], p[(k + 1) % p.len()]);
        s += a[0] * b[1] - b[0] * a[1];
    }
    0.5 * s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Geometry;

    fn total_volume(m: &Mesh) -> f64 {
        let g = Geometry::compute(m);
        g.cell_volume.iter().sum()
    }

    /// At every angle the cells fill the box, no cut piece is smaller than
    /// a part of a whole cell, and the boundaries are the box sides.
    #[test]
    fn turned_grids_fill_the_box() {
        for angle in [0.0, 17.0, 33.0, 45.0, 78.0, 90.0] {
            let m = generate([0.0, 0.0, 0.0], [2.0, 1.0, 0.1], [20, 10, 1], angle, None).unwrap();
            assert!((total_volume(&m) - 0.2).abs() < 1e-12, "angle {angle}");
            let g = Geometry::compute(&m);
            let min = g.cell_volume.iter().copied().fold(f64::INFINITY, f64::min);
            assert!(min > 0.2 * 0.01 * 0.1, "angle {angle}: smallest cell {min}");
            for name in ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"] {
                assert!(m.boundaries.iter().any(|b| b.name == name && !b.range().is_empty()));
            }
            // the whole cells are rotated squares of the spacing
            let n = m.n_cells as f64;
            assert!(n > 150.0 && n < 260.0, "angle {angle}: {n} cells");
        }
    }

    /// The revolved grid is the sector of the box: the volume of the
    /// revolved rectangle, ½ θ (y₁² − y₀²)(x₁ − x₀).
    #[test]
    fn a_turned_sector_has_the_volume_of_the_ring() {
        let th = 22.5f64;
        let m = generate([0.0, 0.1, 0.0], [2.0, 1.0, 1.0], [20, 9, 4], 33.0, Some(th)).unwrap();
        // 4 angular cells with chords: ½ sin δ (y₁² − y₀²) per unit length
        let d = th.to_radians() / 4.0;
        let chord = 4.0 * 0.5 * d.sin() * (1.0 - 0.01) * 2.0;
        let exact = chord;
        assert!((total_volume(&m) - chord).abs() < 1e-3 * exact, "{} vs {chord}", total_volume(&m));
    }
}
