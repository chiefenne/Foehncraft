//! Benchmark data: mesh, fields and matrices, built once on the host.
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use crate::matrix::{Csr, SellCs};
use crate::mesh::{FaceMesh, Ordering};
use foehncraft_core::rng;
use serde::Serialize;
use std::fmt;

/// The kernels of gate G1.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Kernel {
    /// Pass 1 of the flux loop: one flux per face.
    FaceFlux,
    /// Pass 2 of the flux loop: net flux per cell, a gather over its faces.
    CellGather,
    /// y = A x, CSR.
    SpmvCsr,
    /// y = A x, SELL-C-σ.
    SpmvSell,
}

impl Kernel {
    /// All kernels, in report order.
    pub const ALL: [Kernel; 4] = [
        Self::FaceFlux,
        Self::CellGather,
        Self::SpmvCsr,
        Self::SpmvSell,
    ];

    /// Name used in case files and reports.
    pub fn name(self) -> &'static str {
        match self {
            Self::FaceFlux => "face-flux",
            Self::CellGather => "cell-gather",
            Self::SpmvCsr => "spmv-csr",
            Self::SpmvSell => "spmv-sell",
        }
    }

    /// Parses a name from a case file.
    pub fn parse(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|k| k.name() == s)
    }

    /// Name of the GPU entry point (identical in the CUDA reference and the
    /// Rust PTX module).
    pub fn entry(self) -> &'static str {
        match self {
            Self::FaceFlux => "g1_face_flux",
            Self::CellGather => "g1_cell_gather",
            Self::SpmvCsr => "g1_spmv_csr",
            Self::SpmvSell => "g1_spmv_sell",
        }
    }
}

impl fmt::Display for Kernel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// Everything the kernels read, in float64 (D7).
pub struct Workload {
    /// Mesh connectivity.
    pub mesh: FaceMesh,
    /// Box dimensions (of the hexahedral grid under a stage-1 mesh).
    pub dims: (usize, usize, usize),
    /// Cell ordering.
    pub ordering: Ordering,
    /// Cell values of the transported scalar.
    pub phi: Vec<f64>,
    /// Diffusion coefficient per face.
    pub gamma: Vec<f64>,
    /// Mass flux per face.
    pub mass_flux: Vec<f64>,
    /// Face fluxes: the input of [`Kernel::CellGather`], computed on the
    /// host so the gather can be timed on its own.
    pub flux: Vec<f64>,
    /// Matrix in CSR form.
    pub csr: Csr<f64>,
    /// The same matrix in SELL-C-σ form.
    pub sell: SellCs<f64>,
    /// Vector x of the SpMV.
    pub x: Vec<f64>,
}

impl Workload {
    /// Builds the workload on a periodic box of `nx × ny × nz` cells.
    pub fn new(
        dims: (usize, usize, usize),
        ordering: Ordering,
        seed: u64,
        chunk: usize,
        sigma: usize,
    ) -> Self {
        let mesh = FaceMesh::periodic_box(dims.0, dims.1, dims.2, ordering, seed);
        Self::on_mesh(mesh, dims, ordering, seed, chunk, sigma)
    }

    /// Builds the workload on a stage-1 test mesh of the unit cube whose
    /// hexahedral grid has `dims` intervals.
    pub fn stage1(
        kind: foehncraft_mesh::generate::Kind,
        dims: (usize, usize, usize),
        ordering: Ordering,
        seed: u64,
        chunk: usize,
        sigma: usize,
    ) -> anyhow::Result<Self> {
        let spec = foehncraft_mesh::generate::BoxSpec::new(kind, 3, [dims.0, dims.1, dims.2]);
        let mut m = spec.generate().map_err(anyhow::Error::msg)?;
        if ordering == Ordering::Shuffled {
            let perm = foehncraft_mesh::renumber::shuffled(m.n_cells, seed);
            m = foehncraft_mesh::renumber::apply(&m, &perm).map_err(anyhow::Error::msg)?;
        }
        let mesh = FaceMesh::from_mesh(&m);
        Ok(Self::on_mesh(mesh, dims, ordering, seed, chunk, sigma))
    }

    fn on_mesh(
        mesh: FaceMesh,
        dims: (usize, usize, usize),
        ordering: Ordering,
        seed: u64,
        chunk: usize,
        sigma: usize,
    ) -> Self {
        let u = |i: usize, stream: u64| rng::uniform([seed, i as u64, stream]);
        let phi: Vec<f64> = (0..mesh.n_cells).map(|i| u(i, 1)).collect();
        let gamma: Vec<f64> = (0..mesh.n_faces()).map(|f| 0.5 + u(f, 2)).collect();
        let mass_flux: Vec<f64> = (0..mesh.n_faces()).map(|f| u(f, 3) - 0.5).collect();
        let x: Vec<f64> = (0..mesh.n_cells).map(|i| u(i, 4) - 0.5).collect();
        let flux: Vec<f64> = (0..mesh.n_faces())
            .map(|f| unsafe {
                // SAFETY: the mesh is validated by `periodic_box` and
                // `from_mesh`.
                foehncraft_kernels::face::face_flux(
                    f,
                    &mesh.owner,
                    &mesh.neighbour,
                    &phi,
                    &gamma,
                    &mass_flux,
                )
            })
            .collect();
        let csr = Csr::from_mesh(&mesh, &gamma, 1e-2);
        let sell = SellCs::from_csr(&csr, chunk, sigma);
        Self {
            mesh,
            dims,
            ordering,
            phi,
            gamma,
            mass_flux,
            flux,
            csr,
            sell,
            x,
        }
    }

    /// Number of output values of a kernel.
    pub fn output_len(&self, k: Kernel) -> usize {
        match k {
            Kernel::FaceFlux => self.mesh.n_faces(),
            _ => self.mesh.n_cells,
        }
    }

    /// Bytes a kernel must move at least (each array read or written once;
    /// gathered values counted once per distinct value). Effective bandwidth
    /// is this number over the run time; it is a lower bound on the real
    /// traffic, so it rewards good cache use.
    pub fn bytes(&self, k: Kernel) -> f64 {
        let nc = self.mesh.n_cells as f64;
        let nf = self.mesh.n_faces() as f64;
        let (i, r) = (4.0, 8.0); // u32 index, f64 value
        match k {
            Kernel::FaceFlux => nf * (2.0 * i + 2.0 * r + r) + nc * r,
            Kernel::CellGather => nc * i + 2.0 * nf * i + nf * r + nc * r,
            Kernel::SpmvCsr => nc * i + self.csr.nnz() as f64 * (i + r) + 2.0 * nc * r,
            Kernel::SpmvSell => {
                let slices = self.sell.slice_ptr.len() as f64;
                let perm = if self.sell.sigma > 1 { nc * i } else { 0.0 };
                2.0 * slices * i + self.sell.stored() as f64 * (i + r) + 2.0 * nc * r + perm
            }
        }
    }
}
