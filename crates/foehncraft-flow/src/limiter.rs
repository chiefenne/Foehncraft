//! The multi-dimensional limiting process (MLP, Park & Kim 2012) for the
//! scalars and the energy equation: the vertex stencil of a mesh and the
//! limiter functions of the settings.
//!
//! A face-neighbour limiter (Barth–Jespersen, Venkatakrishnan) keeps the
//! reconstructed face values of a cell inside the range of its face
//! neighbours. On a multi-dimensional mesh the face value that flows into
//! a cell then still comes from the range of the upwind cell's neighbours,
//! which are not all neighbours of the receiving cell, so a converged
//! solution can carry new extrema (the 12 K overshoot of the hot-air nozzle
//! of the heat exchanger on tetrahedra). MLP bounds the reconstruction at
//! every vertex by the cells that share the vertex: a face value is then
//! inside the range of every cell that touches the face's vertices, which
//! is local extremum diminishing in multi-dimensions.

use crate::device::VertexStencil;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_mesh::Mesh;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// The limiter function of `limiter = "mlp"` (`limiter_function`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum LimiterFunction {
    /// ψ = min(1, Δ₊/Δ₋): strictly bounded, not differentiable (can stall a
    /// steady residual).
    BarthJespersen,
    /// Venkatakrishnan's function (second order in smooth regions).
    Venkatakrishnan,
    /// Nishikawa's R3 (2022): third order in smooth regions,
    /// differentiable.
    #[default]
    R3,
    /// Nishikawa's R4.
    R4,
    /// Nishikawa's R5.
    R5,
}

impl LimiterFunction {
    /// The code of [`foehncraft_kernels::fv::mlp_limiter`].
    pub fn code(self) -> u32 {
        match self {
            Self::BarthJespersen => foehncraft_kernels::fv::MLP_BARTH_JESPERSEN,
            Self::Venkatakrishnan => foehncraft_kernels::fv::MLP_VENKATAKRISHNAN,
            Self::R3 => 3,
            Self::R4 => 4,
            Self::R5 => 5,
        }
    }
}

/// How ε of `limiter = "mlp"` scales (`limiter_epsilon_scaling`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum EpsilonScaling {
    /// ε = c (φ_max − φ_min): one value for the domain (Wang 2000).
    #[default]
    Domain,
    /// ε_c = c (φ_max − φ_min) √(V_c / V_domain): shrinks with the cell
    /// size (as h^1.5 in 3D, h in a 2D case one cell thick; Moe et al.
    /// 2015, Setzwein et al. 2026), so the bounds tighten as the mesh is
    /// refined.
    Cell,
}

/// The factors √(V_c / V_domain) of [`EpsilonScaling::Cell`]: `volume`
/// of the own cells; on ranks (`comm`) V_domain is the exact sum over all
/// ranks (the same value as on one rank).
pub fn cell_factors(
    volume: &[f64],
    comm: Option<&foehncraft_mesh::decompose::SharedComm>,
) -> Vec<f64> {
    let total: f64 = match comm {
        Some(c) => {
            let mut s = foehncraft_core::exact::ExactSum::new();
            for &v in volume {
                s.add(v);
            }
            c.sum_exact(&s)
        }
        None => volume.iter().exact_sum(),
    };
    volume.iter().map(|v| (v / total).sqrt()).collect()
}

/// The vertex stencil of `mesh`: the cells and boundary faces of every
/// vertex and the vertices of every cell (each list sorted, no
/// duplicates). A periodic face (each side keeps its own face) counts for
/// its owner only, so the bounds at a vertex on a periodic plane come from
/// the cells on its own side: tighter, never looser. Vertices no cell uses
/// (the points of the solid cells in the flow mesh of a conjugate case)
/// keep an empty list; no cell reads them.
///
/// The part of a rank (MPI step W7): the own cells only (a processor face
/// counts for its own side, the halo cells for no vertex); the vertices of
/// the processor faces are `shared`, and the ranks that share one combine
/// their bounds into those of the whole mesh (a minimum and a maximum: the
/// same in any order).
pub fn vertex_stencil(mesh: &Mesh) -> VertexStencil {
    let nv = mesh.points.len();
    let nc = mesh.n_cells;
    let inner = mesh.n_inner_faces;
    let processor_end = mesh.periodic_start();
    let ni = mesh.n_interior_faces();
    let nf = mesh.n_faces();
    let mut of_cell: Vec<Vec<u32>> = vec![Vec::new(); nc];
    let mut bf_of_vert: Vec<Vec<u32>> = vec![Vec::new(); nv];
    let mut shared = Vec::new();
    for f in 0..nf {
        let nodes = mesh.face(f);
        let o = mesh.owner[f] as usize;
        if o < nc {
            of_cell[o].extend_from_slice(nodes);
        }
        if f < inner {
            of_cell[mesh.neighbour[f] as usize].extend_from_slice(nodes);
        } else if f < processor_end {
            // a processor face: its own side (the owner may be the halo)
            let n = mesh.neighbour[f] as usize;
            if n < nc {
                of_cell[n].extend_from_slice(nodes);
            }
            shared.extend_from_slice(nodes);
        } else if f >= ni {
            for &v in nodes {
                bf_of_vert[v as usize].push((f - ni) as u32);
            }
        }
    }
    shared.sort_unstable();
    shared.dedup();
    let shared_id = shared
        .iter()
        .map(|&v| {
            mesh.point_id
                .get(v as usize)
                .copied()
                .unwrap_or(u64::from(v))
        })
        .collect();
    let mut cells_of_vert: Vec<Vec<u32>> = vec![Vec::new(); nv];
    let mut cell_ptr = Vec::with_capacity(nc + 1);
    let mut cell_verts = Vec::new();
    cell_ptr.push(0u32);
    for (c, vs) in of_cell.iter_mut().enumerate() {
        vs.sort_unstable();
        vs.dedup();
        for &v in vs.iter() {
            cells_of_vert[v as usize].push(c as u32);
        }
        cell_verts.extend_from_slice(vs);
        cell_ptr.push(cell_verts.len() as u32);
    }
    let csr = |lists: &mut [Vec<u32>]| {
        let mut ptr = Vec::with_capacity(lists.len() + 1);
        let mut out = Vec::new();
        ptr.push(0u32);
        for l in lists.iter_mut() {
            l.sort_unstable();
            l.dedup();
            out.extend_from_slice(l);
            ptr.push(out.len() as u32);
        }
        (ptr, out)
    };
    let (vert_ptr, vert_cells) = csr(&mut cells_of_vert);
    let (bface_ptr, bfaces) = csr(&mut bf_of_vert);
    let points = [
        mesh.points.iter().map(|p| p[0]).collect(),
        mesh.points.iter().map(|p| p[1]).collect(),
        mesh.points.iter().map(|p| p[2]).collect(),
    ];
    VertexStencil {
        vert_ptr,
        vert_cells,
        bface_ptr,
        bfaces,
        cell_ptr,
        cell_verts,
        points,
        shared,
        shared_id,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_strip_of_two_hexahedra_shares_four_vertices() {
        let src = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [2.0, 1.0, 1.0]
cells = [2, 1, 1]
"#;
        let case: crate::case::FlowCase = toml::from_str(src).unwrap();
        let (mesh, _, _) = crate::case::mesh_for(&case, None).unwrap();
        let s = vertex_stencil(&mesh);
        assert_eq!(s.cell_ptr, vec![0, 8, 16]);
        let shared = (0..mesh.points.len())
            .filter(|&v| s.vert_ptr[v + 1] - s.vert_ptr[v] == 2)
            .count();
        assert_eq!(shared, 4);
        // every vertex of the strip is on the boundary
        assert!((0..mesh.points.len()).all(|v| s.bface_ptr[v + 1] > s.bface_ptr[v]));
    }
}
