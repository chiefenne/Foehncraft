//! The walls of a turbulent flow (stage 6 specification 4, D134; D48):
//! the wall faces with their first-cell distance, normal and roughness,
//! the near-wall cells whose ω is fixed (4.4), and the wall distance,
//! exact (D48, the BVH of stage 6a) or from the Poisson equation
//! ∇²φ = −1, φ = 0 on walls, ∂φ/∂n = 0 elsewhere, with
//! d = −|∇φ| + (|∇φ|² + 2φ)^½ (P. G. Tucker, *Differential equation-based
//! wall distance computation for DES and RANS*, J. Comput. Phys. 190,
//! 229–248 (2003); exact for a plane wall).

use crate::bc::BoundaryCondition;
use crate::device::WallStatic;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_mesh::{Geometry, Mesh};
use std::collections::BTreeMap;

/// The wall faces of a turbulent flow.
#[derive(Debug, Clone, Default)]
pub(crate) struct Walls {
    /// Device arrays of the wall faces and the near-wall cells.
    pub stat: WallStatic,
    /// Outward unit normal per wall face.
    pub normal: [Vec<f64>; 3],
    /// First-cell distance y_P = |(x_f − x_P)·n̂| per wall face, at least
    /// 10⁻³ of the cell's thickness V_P/A_f normal to the face (not of
    /// V_P^⅓, which on thin wall cells of a unit-depth 2D mesh exceeds y_P
    /// itself and lowers ω_P by an order of magnitude).
    pub y: Vec<f64>,
    /// k_s per wall face (0: smooth).
    pub ks: Vec<f64>,
    /// C_s of Cebeci & Bradshaw per wall face.
    pub cs: Vec<f64>,
    /// Weight of each wall face in the mean of its cell (its area over the
    /// smooth wall area of the cell; 0 for rough faces of SST, whose ω
    /// takes the rough-wall value as a boundary value; the k-ε family
    /// counts rough faces like smooth ones).
    pub weight: Vec<f64>,
    /// 1 in the cells with a smooth wall face (ω fixed), per cell.
    pub fix_mask: Vec<f64>,
    /// Face area per wall face.
    pub area: Vec<f64>,
    /// Boundary index per wall face.
    pub boundary: Vec<usize>,
    /// 1 on the wall faces, per boundary face.
    pub mask_b: Vec<f64>,
    /// The exact wall distance of the sample cell per wall face (with
    /// `stat.sample`, cut-cell walls).
    pub y_sample: Vec<f64>,
}

impl Walls {
    /// The wall faces of the conditions `bcs` on `mesh`; `rough_fixed`:
    /// rough faces fix the cell value as smooth ones do (the k-ε family).
    pub fn new(
        mesh: &Mesh,
        geo: &Geometry,
        bcs: &BTreeMap<String, BoundaryCondition>,
        rough_fixed: bool,
    ) -> Self {
        let nc = mesh.n_cells;
        let ni = mesh.n_interior_faces();
        let nb = mesh.n_boundary_faces();
        let mut w = Walls {
            fix_mask: vec![0.0; nc],
            mask_b: vec![0.0; nb],
            ..Default::default()
        };
        for (bi, p) in mesh.boundaries.iter().enumerate() {
            let Some(bc @ BoundaryCondition::Wall { .. }) = bcs.get(&p.name) else {
                continue;
            };
            let ks = bc.roughness();
            let cs = bc.roughness_constant();
            for f in p.range() {
                let s = geo.face_area[f];
                let a = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
                let n = [s[0] / a, s[1] / a, s[2] / a];
                let c = mesh.owner[f] as usize;
                let xf = geo.face_centroid[f];
                let xc = geo.cell_centroid[c];
                let d = (0..3).map(|i| (xf[i] - xc[i]) * n[i]).sum::<f64>().abs();
                let floor = 1e-3 * geo.cell_volume[c] / a;
                w.stat.faces.push(f as u32);
                for (i, ni_) in n.iter().enumerate() {
                    w.normal[i].push(*ni_);
                }
                w.y.push(d.max(floor));
                w.ks.push(ks);
                w.cs.push(cs);
                w.area.push(a);
                w.boundary.push(bi);
                w.mask_b[f - ni] = 1.0;
            }
        }
        // the near-wall cells (owners of smooth wall faces), ascending, and
        // their wall faces in face order (CSR)
        let mut per_cell: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
        for (slot, &f) in w.stat.faces.iter().enumerate() {
            let c = mesh.owner[f as usize];
            per_cell.entry(c).or_default().push(slot as u32);
        }
        w.weight = vec![0.0; w.stat.faces.len()];
        w.stat.ptr.push(0);
        for (c, slots) in per_cell {
            let fixes = |s: u32| rough_fixed || w.ks[s as usize] <= 0.0;
            let smooth: f64 = slots
                .iter()
                .filter(|&&s| fixes(s))
                .map(|&s| w.area[s as usize])
                .exact_sum();
            if smooth <= 0.0 {
                // only rough faces: ω from the transport equation with the
                // wall value of the roughness model
                continue;
            }
            for &s in &slots {
                if fixes(s) {
                    w.weight[s as usize] = w.area[s as usize] / smooth;
                }
            }
            w.fix_mask[c as usize] = 1.0;
            w.stat.cells.push(c);
            w.stat.slots.extend(&slots);
            w.stat.ptr.push(w.stat.slots.len() as u32);
        }
        // the device arrays address every wall face by its slot: the CSR
        // lists only the slots of the fixed cells, the face list all of them
        w
    }

    /// The walls of cut bodies (L13 stage 2c, D229).
    ///
    /// **Exact first-cell distance.** y_P of every face of a body's wall
    /// boundary is the distance of the cell centroid to the true surface
    /// (the closest point, D108), the same for all wall faces of a cut or
    /// merged cell, instead of the distance to the plane of each facet; the
    /// floor stays.
    ///
    /// **Sample cell** (`[immersed.cut] wall_sample`, default true). The
    /// cut cells along a wall are cut anywhere from a sliver to the whole
    /// cell, so their centroids lie at any distance from 0.15 to 1.5 cell
    /// widths, and their velocities follow the log law poorly (the ragged
    /// first layer: ν_t, the production and the momentum fluxes of cells at
    /// unequal heights). The wall function of each cut wall face then
    /// takes its velocity, k and distance from the cell one layer out (the
    /// cell holding the point just above the cut cell along the wall
    /// normal), whose centroid lies between one and two cell widths from
    /// the wall, and applies the wall shear stress it gives to the cut
    /// cell's face (ν_w = τ_w y_P/U_P); the cut cell keeps its own ε or ω.
    /// Wall models of immersed-boundary RANS sample the velocity at a
    /// distance from the wall in the same way (Capizzano 2011; Tamaki,
    /// Harada & Imamura 2017; Berger & Aftosmis 2012).
    ///
    /// Returns the number of faces set and the largest relative change of
    /// y_P.
    pub fn set_cut_walls(
        &mut self,
        mesh: &Mesh,
        geo: &Geometry,
        bodies: &crate::immersed::WallBodies,
    ) -> (usize, f64) {
        let (mut n, mut change) = (0, 0.0f64);
        let nw = self.stat.faces.len();
        let sample = bodies.sample();
        let loc = sample.then(|| foehncraft_mesh::locate::Locator::new(mesh, geo));
        let mut cells = vec![u32::MAX; nw];
        let mut ys = vec![0.0; nw];
        let ni = mesh.n_interior_faces();
        let neighbours = |c: usize| -> Vec<usize> {
            mesh.cell_face_iter(c)
                .filter(|&(f, _)| f < ni)
                .map(|(f, is_nb)| {
                    if is_nb {
                        mesh.owner[f] as usize
                    } else {
                        mesh.neighbour[f] as usize
                    }
                })
                .collect()
        };
        for (i, &f) in self.stat.faces.iter().enumerate() {
            let name = &mesh.boundaries[self.boundary[i]].name;
            let c = mesh.owner[f as usize] as usize;
            let xc = geo.cell_centroid[c];
            let Some(d) = bodies.distance_to(name, xc) else {
                continue;
            };
            let floor = 1e-3 * geo.cell_volume[c] / self.area[i];
            let y = d.max(floor);
            change = change.max((y - self.y[i]).abs() / y);
            self.y[i] = y;
            n += 1;
            let Some(loc) = &loc else { continue };
            // the point just above the cut cell along the inward normal
            let nn = [-self.normal[0][i], -self.normal[1][i], -self.normal[2][i]];
            let mut top = 0.0f64;
            for (g, _) in mesh.cell_face_iter(c) {
                for &p in mesh.face(g) {
                    let x = mesh.points[p as usize];
                    let h = (0..3).map(|k| (x[k] - xc[k]) * nn[k]).sum::<f64>();
                    top = top.max(h);
                }
            }
            let probe: [f64; 3] = std::array::from_fn(|k| xc[k] + nn[k] * (top + 1e-3 * (top + y)));
            // in a face neighbour, or a neighbour of one
            let mut found = None;
            let first = neighbours(c);
            for &q in &first {
                if loc.contains(q, probe) {
                    found = Some(q);
                    break;
                }
            }
            if found.is_none() {
                'outer: for &q in &first {
                    for r in neighbours(q) {
                        if r != c && loc.contains(r, probe) {
                            found = Some(r);
                            break 'outer;
                        }
                    }
                }
            }
            let Some(q) = found.filter(|&q| q < mesh.n_cells) else {
                continue;
            };
            let Some(dq) = bodies.distance_to(name, geo.cell_centroid[q]) else {
                continue;
            };
            if dq > y * (1.0 + 1e-9) {
                cells[i] = q as u32;
                ys[i] = dq;
            }
        }
        if sample && cells.iter().any(|&q| q != u32::MAX) {
            self.stat.sample = cells;
            self.y_sample = ys;
        }
        (n, change)
    }

    /// Number of wall faces.
    pub fn len(&self) -> usize {
        self.stat.faces.len()
    }

    /// No wall faces.
    pub fn is_empty(&self) -> bool {
        self.stat.faces.is_empty()
    }
}

/// d = −|∇φ| + (|∇φ|² + 2φ)^½ of the Poisson wall distance, per cell.
pub(crate) fn poisson_distance(phi: &[f64], grad: [&[f64]; 3]) -> Vec<f64> {
    phi.iter()
        .enumerate()
        .map(|(c, &p)| {
            let g2 = grad[0][c] * grad[0][c] + grad[1][c] * grad[1][c] + grad[2][c] * grad[2][c];
            let g = g2.sqrt();
            (-g + (g2 + 2.0 * p.max(0.0)).sqrt()).max(0.0)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poisson_distance_is_exact_for_a_plane_wall() {
        // channel of half height h: φ = y (2h − y)/2, ∇φ = h − y
        let h = 0.7;
        let ys = [0.01, 0.1, 0.3, 0.69];
        let phi: Vec<f64> = ys.iter().map(|y| y * (2.0 * h - y) / 2.0).collect();
        let gy: Vec<f64> = ys.iter().map(|y| h - y).collect();
        let z = vec![0.0; ys.len()];
        let d = poisson_distance(&phi, [&z, &gy, &z]);
        for (di, y) in d.iter().zip(ys) {
            assert!((di - y).abs() < 1e-14, "{di} {y}");
        }
    }

    /// The first-cell distance on a cylinder cut out of a box mesh (L13
    /// stage 2c): with the exact distance every wall face of a cut or
    /// merged cell takes the distance of the cell centroid to the circle,
    /// |r_P − R| to round-off; the distance to the plane of each facet
    /// differs from it, by the most on merged cells.
    #[test]
    fn cut_walls_take_the_exact_first_cell_distance() {
        let (r, c) = (0.1234, [0.5, 0.5]);
        let text = format!(
            r#"
case_format = "0.7"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [40, 40, 1]
[fluid]
density = 1.0
viscosity = 1e-5
pressure_values = "gauge"
[immersed]
method = "cut"
[[immersed.body]]
name = "cyl"
shape = "cylinder"
centre = [{}, {}, 0.0]
radius = {r}
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "symmetry"
[boundary.ymax]
type = "symmetry"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[boundary.cyl]
type = "wall"
"#,
            c[0], c[1]
        );
        let case = crate::case::FlowCase::from_toml(&text, std::path::Path::new(".")).unwrap();
        let (mesh, geo, _) = crate::case::mesh_for(&case, None).unwrap();
        let facet = Walls::new(&mesh, &geo, &case.boundary, true);
        let mut exact = facet.clone();
        let wb = crate::immersed::WallBodies::of(&case).unwrap().unwrap();
        let (n, change) = exact.set_cut_walls(&mesh, &geo, &wb);
        assert_eq!(n, facet.len());
        assert!(n > 100, "{n} wall faces");
        let mut per_cell: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
        let (mut worst, mut spread) = (0.0f64, 0.0f64);
        for (i, &f) in exact.stat.faces.iter().enumerate() {
            let p = mesh.owner[f as usize];
            per_cell.entry(p).or_default().push(i);
            let x = geo.cell_centroid[p as usize];
            let rp = ((x[0] - c[0]).powi(2) + (x[1] - c[1]).powi(2)).sqrt();
            let want = (rp - r).abs();
            assert!(
                (exact.y[i] - want).abs() <= 1e-12 * r,
                "{} {want}",
                exact.y[i]
            );
            worst = worst.max((facet.y[i] - want).abs() / want);
        }
        for faces in per_cell.values() {
            let y: Vec<f64> = faces.iter().map(|&i| facet.y[i]).collect();
            let lo = y.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = y.iter().copied().fold(0.0, f64::max);
            spread = spread.max((hi - lo) / hi);
            let e: Vec<f64> = faces.iter().map(|&i| exact.y[i]).collect();
            assert!(e.iter().all(|&v| v == e[0]), "one y per cell: {e:?}");
        }
        // the facets put the wall a few per cent off the true distance
        assert!(worst > 1e-3 && (change - worst).abs() < 0.5 * worst + 1e-9);
        assert!(spread > 1e-3, "{spread}");
        eprintln!(
            "cut cylinder 40²: {n} wall faces, facet distance up to {:.2} % off, \
             spread within a cell up to {:.2} %",
            100.0 * worst,
            100.0 * spread
        );
    }
}
