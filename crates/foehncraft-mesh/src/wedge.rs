//! Axisymmetric wedge meshes (D3, D29): a 2D mesh of the (x, r) half plane,
//! one cell layer thick in z, revolved about the x axis into a wedge of a
//! small angle θ. The node (x, y, z) of the 2D mesh goes to
//! (x, y cos φ, y sin φ) with φ = θ (z − z₀)/(z₁ − z₀) − θ/2, so y is the
//! radius and the wedge is symmetric about the x–y plane.
//!
//! Nodes on the axis (y = 0) of both layers become one node: the cells
//! next to the axis lose their face on it (a face of zero area) and become
//! prisms (from hexahedra) or polyhedra with an edge on the axis, so the
//! mesh has no degenerate faces and no axis boundary. The boundary `ymin`
//! disappears when it lies on the axis (`lower[1] = 0`); with
//! `lower[1] > 0` (an annulus) it stays. `zmin` and `zmax` become the two
//! wedge planes (at −θ/2 and +θ/2); the other boundaries keep their names.
//!
//! A 3D box of hexahedra revolves the same way into a sector (D309): its
//! z intervals become the angular cells, every face stays planar (faces of constant x
//! lie in a plane x = const, faces of constant y are chords, faces of
//! constant z lie in a plane through the axis), and `zmin`, `zmax` are the
//! side planes of the sector.
//!
//! Every cell type of the box generator revolves: hexahedra, prisms on a
//! triangulation, polygonal prisms and the mixed bands, perturbed or not
//! (the perturbation moves nodes in x and y only, so the wedge planes stay
//! planar and axis nodes stay on the axis).

use crate::builder::{BoundaryFace, CellSoup};
use crate::generate::BOX_BOUNDARIES;
use crate::mesh::Mesh;
use std::collections::HashMap;

/// Revolves a 2D mesh (one layer in z, y ≥ 0) into a wedge of `angle`
/// degrees about the x axis, or a 3D box mesh into a sector whose z
/// intervals are the angular cells.
pub fn revolve(m: &Mesh, angle: f64) -> Result<Mesh, String> {
    if !(angle > 0.0 && angle <= 90.0) {
        return Err(format!(
            "the wedge angle must lie in (0°, 90°], not {angle}°"
        ));
    }
    let (mut z0, mut z1, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
    for p in &m.points {
        z0 = z0.min(p[2]);
        z1 = z1.max(p[2]);
        ymax = ymax.max(p[1]);
        if p[1] < 0.0 {
            return Err(format!(
                "a wedge mesh needs y >= 0 (y is the radius); node at {p:?}"
            ));
        }
    }
    if z1 <= z0 || z1.is_nan() || z0.is_nan() {
        return Err("the 2D mesh of a wedge needs a layer thickness in z".into());
    }
    let th = angle.to_radians();
    let on_axis = |y: f64| y <= 1e-14 * ymax;
    // the new nodes: one per axis position, two per off-axis node
    let mut map = vec![0u32; m.points.len()];
    let mut pts: Vec<[f64; 3]> = Vec::with_capacity(m.points.len());
    let mut axis: HashMap<u64, u32> = HashMap::new();
    for (i, p) in m.points.iter().enumerate() {
        if on_axis(p[1]) {
            map[i] = *axis.entry(p[0].to_bits()).or_insert_with(|| {
                pts.push([p[0], 0.0, 0.0]);
                (pts.len() - 1) as u32
            });
        } else {
            let phi = th * (p[2] - z0) / (z1 - z0) - 0.5 * th;
            map[i] = pts.len() as u32;
            pts.push([p[0], p[1] * phi.cos(), p[1] * phi.sin()]);
        }
    }
    // a face with its nodes mapped and repeated nodes removed; None when
    // fewer than three remain (the faces on the axis)
    let face = |f: usize| -> Option<Vec<u32>> {
        let mut v: Vec<u32> = Vec::with_capacity(8);
        for &n in m.face(f) {
            let k = map[n as usize];
            if v.last() != Some(&k) {
                v.push(k);
            }
        }
        while v.len() > 1 && v.first() == v.last() {
            v.pop();
        }
        (v.len() >= 3).then_some(v)
    };
    // the boundary of each boundary face that survives, by its sorted nodes
    let mut names: HashMap<Vec<u32>, String> = HashMap::new();
    let mut kept: Vec<String> = Vec::new();
    for p in &m.boundaries {
        for f in p.range() {
            if let Some(mut v) = face(f) {
                v.sort_unstable();
                names.insert(v, p.name.clone());
                if !kept.contains(&p.name) {
                    kept.push(p.name.clone());
                }
            }
        }
    }
    let mut soup = CellSoup::new(pts);
    for c in 0..m.n_cells {
        let faces: Vec<Vec<u32>> = m.cell_face_iter(c).filter_map(|(f, _)| face(f)).collect();
        soup.add_cell(faces.iter().map(|v| &v[..]));
    }
    // boundary order: that of the box, then any other
    let mut order: Vec<&str> = BOX_BOUNDARIES
        .iter()
        .copied()
        .filter(|n| kept.iter().any(|k| k == n))
        .collect();
    for k in &kept {
        if !order.contains(&k.as_str()) {
            order.push(k);
        }
    }
    let mut out = soup.finish(&order, |b: BoundaryFace| {
        let mut v = b.nodes.to_vec();
        v.sort_unstable();
        names
            .get(&v)
            .cloned()
            .unwrap_or_else(|| "unknown".to_string())
    })?;
    out.zones.clone_from(&m.zones);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Geometry;
    use crate::generate::{BoxSpec, Kind};

    fn spec(kind: Kind, perturb: f64) -> BoxSpec {
        BoxSpec {
            kind,
            dimension: 2,
            lower: [0.0, 0.0, 0.0],
            upper: [2.0, 0.5, 0.1],
            cells: [8, 6, 1],
            perturb,
            ..Default::default()
        }
    }

    /// The volume of every type is that of the planar wedge (the faces are
    /// planar polygons, so the revolved rectangle becomes a prism with the
    /// section ½ R² sin θ), every cell is closed, the wedge planes are
    /// planar and the boundaries are named.
    #[test]
    fn wedges_of_every_type() {
        let angle: f64 = 5.0;
        let th = angle.to_radians();
        for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
            for perturb in [0.0, 0.2] {
                let m2 = spec(kind, perturb).generate().unwrap();
                let m = revolve(&m2, angle).unwrap();
                let g = Geometry::compute(&m);
                assert!(g.cell_volume.iter().all(|&v| v > 0.0), "{kind:?}");
                // the revolved rectangle: the planar wedge of the radius
                // R = 0.5 has the section ½ R² sin θ, length 2
                let vol: f64 = g.cell_volume.iter().sum();
                let exact = 2.0 * 0.5 * 0.5 * 0.5 * th.sin();
                assert!(
                    (vol - exact).abs() < 1e-12 * exact,
                    "{kind:?}: {vol} {exact}"
                );
                let names: Vec<&str> = m.boundaries.iter().map(|p| p.name.as_str()).collect();
                assert_eq!(names, ["xmin", "xmax", "ymax", "zmin", "zmax"], "{kind:?}");
                // closure of every cell
                for c in 0..m.n_cells {
                    let mut s = [0.0; 3];
                    for (f, nb) in m.cell_face_iter(c) {
                        let a = g.face_area[f];
                        let sg = if nb { -1.0 } else { 1.0 };
                        for i in 0..3 {
                            s[i] += sg * a[i];
                        }
                    }
                    let n = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
                    assert!(n < 1e-15, "{kind:?} cell {c}: {n}");
                }
                // the wedge planes are planar at ∓θ/2
                for (name, sgn) in [("zmin", -1.0), ("zmax", 1.0)] {
                    let p = m.boundary(name).unwrap();
                    for f in p.range() {
                        let a = g.face_area[f];
                        let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
                        let e = [0.0, -(0.5 * th).sin(), (0.5 * th).cos() * sgn];
                        let d = (0..3).map(|i| (a[i] / n - e[i]).abs()).fold(0.0, f64::max);
                        assert!(d < 1e-12, "{kind:?} {name}: normal off by {d}");
                    }
                }
            }
        }
    }

    /// A 3D box revolves into a sector: its z intervals are the angular
    /// cells, and the volume is that of the polygonal sector, N planar
    /// wedges of the section ½ R² sin(θ/N); every cell is closed and the
    /// side planes lie at ∓θ/2 (D309).
    #[test]
    fn sectors_of_hexahedra() {
        let angle: f64 = 22.5;
        let th = angle.to_radians();
        for kind in [Kind::Hex] {
            {
                let s = BoxSpec {
                    kind,
                    dimension: 3,
                    lower: [0.0, 0.0, 0.0],
                    upper: [2.0, 0.5, 0.1],
                    cells: [6, 5, 4],
                    wedge: Some(angle),
                    ..Default::default()
                };
                let m = s.generate().unwrap();
                let g = Geometry::compute(&m);
                assert!(g.cell_volume.iter().all(|&v| v > 0.0), "{kind:?}");
                let vol: f64 = g.cell_volume.iter().sum();
                let exact = 2.0 * 4.0 * 0.5 * 0.5 * 0.5 * (th / 4.0).sin();
                assert!(
                    (vol - exact).abs() < 1e-12 * exact,
                    "{kind:?}: {vol} {exact}"
                );
                let names: Vec<&str> = m.boundaries.iter().map(|p| p.name.as_str()).collect();
                assert_eq!(names, ["xmin", "xmax", "ymax", "zmin", "zmax"], "{kind:?}");
                for c in 0..m.n_cells {
                    let mut sum = [0.0; 3];
                    for (f, nb) in m.cell_face_iter(c) {
                        let a = g.face_area[f];
                        let sg = if nb { -1.0 } else { 1.0 };
                        for i in 0..3 {
                            sum[i] += sg * a[i];
                        }
                    }
                    let n = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
                    assert!(n < 1e-15, "{kind:?} cell {c}: {n}");
                }
                for (name, sgn) in [("zmin", -1.0), ("zmax", 1.0)] {
                    let p = m.boundary(name).unwrap();
                    for f in p.range() {
                        let a = g.face_area[f];
                        let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
                        let e = [0.0, -(0.5 * th).sin(), (0.5 * th).cos() * sgn];
                        let d = (0..3).map(|i| (a[i] / n - e[i]).abs()).fold(0.0, f64::max);
                        assert!(d < 1e-12, "{kind:?} {name}: normal off by {d}");
                    }
                }
            }
        }
    }

    /// A perturbed 3D box would move the nodes on the axis apart, and
    /// cells other than hexahedra would collapse on the axis: refused.
    #[test]
    fn a_perturbed_or_tetrahedral_sector_is_refused() {
        let s = BoxSpec {
            dimension: 3,
            cells: [2, 2, 2],
            perturb: 0.1,
            wedge: Some(10.0),
            ..Default::default()
        };
        let e = s.generate().unwrap_err();
        assert!(e.contains("perturb"), "{e}");
        let s = BoxSpec {
            kind: Kind::Tet,
            dimension: 3,
            cells: [2, 2, 2],
            wedge: Some(10.0),
            ..Default::default()
        };
        let e = s.generate().unwrap_err();
        assert!(e.contains("hex"), "{e}");
    }

    /// An annulus keeps its inner boundary.
    #[test]
    fn annulus_keeps_the_inner_wall() {
        let mut s = spec(Kind::Hex, 0.0);
        s.lower[1] = 0.1;
        let m = revolve(&s.generate().unwrap(), 2.0).unwrap();
        assert!(m.boundary("ymin").is_some());
        assert_eq!(m.n_cells, 48);
    }
}
