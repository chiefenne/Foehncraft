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
}
