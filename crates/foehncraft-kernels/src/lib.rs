//! Compute kernels of Foehncraft, written once for every target.
//!
//! - [`face`], [`spmv`]: the kernels of gate G1 (stage 0);
//! - [`fv`]: the finite-volume operators on a scalar (stage 2): gradient,
//!   limiter, face flux coefficients, row assembly;
//! - [`flow`]: the incompressible flow solver (stage 4): momentum faces and
//!   rows, Rhie–Chow mass flux, pressure correction, boundary values;
//! - [`scalar`]: the rows of the generic scalar transport equation
//!   (stage 6, D116): bounded form, point-implicit sources, bounded BDF2;
//! - [`turb`]: the turbulence models (stage 6): invariants and the D166
//!   features, k-ω SST 2003, its wall treatments;
//! - [`particle`]: the parcel step of the Lagrangian sprays (L12): exact
//!   relaxation of the parcel velocity to the gas velocity with the drag of
//!   Schiller and Naumann, gravity and buoyancy;
//! - [`evaporation`]: droplet evaporation of the sprays (SP5): the d² law
//!   and the exact temperature relaxation with frozen coefficients, Ranz and
//!   Marshall and Abramzon and Sirignano, the vapour pressure of water;
//! - [`blas`]: vector operations and reproducible block dot products of the
//!   Krylov solvers;
//! - [`expr`]: the stack machine of the compiled case-file expressions (D53);
//! - [`adjoint`]: the transposed matrix and product for the adjoint (D59);
//! - [`reverse`]: the reverse sweep of a face kernel from local dual
//!   Jacobians, gathered over the inverse stencil (L17 stage A0, D221);
//! - [`amg`]: the smoothers, transfers and coarse solve of the algebraic
//!   multigrid V-cycle (stage 3), generic over float32 and float64;
//! - [`walk`]: the straight walk of a particle or a radiation bundle
//!   through the cells (D47, D257), resumable on another rank;
//! - [`radiation`]: Monte Carlo radiation (step 29): Owen-scrambled Sobol
//!   points, the sampling of the bundles and their transport with exact
//!   integer tallies;
//! - [`chemistry`]: gas-phase reaction rates, their analytic Jacobian and
//!   the stiff RODAS3 integrator of the reactor of a cell (stage SF2);
//! - [`realfluid`]: the Peng-Robinson equation of state of a mixture, its
//!   fugacities and departure enthalpies and the phase equilibrium at a
//!   droplet surface at high pressure (ECN Spray A, SA1);
//! - [`species`]: the properties of an ideal-gas mixture (stage SF1):
//!   NASA-7 specific heat, Chapman–Enskog viscosity, Wilke and Mathur rules.
//!
//! Each kernel is a *body*: a function that computes the result for one
//! index (one face, one cell, one matrix row). A backend supplies the loop:
//! a sequential or multi-threaded loop on the CPU (`foehncraft-bench`), or one
//! GPU thread per index (`gpu/ptx-kernels`, compiled for
//! `nvptx64-nvidia-cuda`). This is the single-source route measured by gate
//! G1 (requirements D8, D11) against a hand-written CUDA reference.
//!
//! Every body is a gather: it reads its inputs and writes one output, with no
//! atomics and a fixed summation order, so results are bitwise reproducible
//! from run to run (D31).
//!
//! # Safety
//!
//! Inside the bodies, loads are unchecked unless the `checked-indexing`
//! feature is on. The connectivity and matrix structures are validated once
//! on the host when they are built (`foehncraft-bench`: `FaceMesh::validate`,
//! `Csr::validate`, `SellCs::validate`); the bodies are `unsafe fn` and their
//! callers guarantee that every index lies inside its array.
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]
#![no_std]

pub mod adjoint;
pub mod amg;
pub mod blas;
pub mod chemistry;
pub mod entry;
pub mod evaporation;
pub mod expr;
pub mod face;
pub mod flow;
pub mod fv;
pub mod particle;
pub mod radiation;
pub mod realfluid;
pub mod reverse;
pub mod scalar;
pub mod species;
pub mod spmv;
pub mod turb;
pub mod walk;

use foehncraft_core::index::FLAG_BIT;

/// Loads `s[i]`; unchecked unless the `checked-indexing` feature is on.
///
/// # Safety
/// `i < s.len()`.
#[inline(always)]
pub(crate) unsafe fn ld<T: Copy>(s: &[T], i: usize) -> T {
    #[cfg(feature = "checked-indexing")]
    {
        s[i]
    }
    #[cfg(not(feature = "checked-indexing"))]
    {
        debug_assert!(i < s.len());
        // SAFETY: guaranteed by the caller.
        unsafe { *s.get_unchecked(i) }
    }
}

/// Splits a packed cell-to-face entry into the face index and the
/// orientation: `true` when the cell is the neighbour of the face, so the
/// face flux leaves the cell with a negative sign.
#[inline(always)]
pub fn unpack_face(entry: u32) -> (usize, bool) {
    ((entry & !FLAG_BIT) as usize, entry & FLAG_BIT != 0)
}

/// Packs a face index and its orientation seen from a cell.
#[inline(always)]
pub fn pack_face(face: u32, cell_is_neighbour: bool) -> u32 {
    debug_assert!(face & FLAG_BIT == 0);
    if cell_is_neighbour {
        face | FLAG_BIT
    } else {
        face
    }
}
