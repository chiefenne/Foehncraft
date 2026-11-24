//! The single-source GPU route: the kernel bodies of `foehncraft-kernels`
//! wrapped as PTX entry points, one GPU thread per output. The `g1_*`
//! entries are the kernels of gate G1 (stage 0); the `fv_*` and `blas_*`
//! entries in [`fv`] are the scalar operators of stage 2, the `amg_*`
//! entries in [`amg`] the multigrid V-cycle of stage 3; `radiation_transport`
//! in [`radiation`] moves the energy bundles of the radiation (step 29);
//! `steam_state` in [`steam`] evaluates water and steam states (SF7a).
//!
//! The wrappers only compute the thread index, build slices from the raw
//! device pointers and their sizes, and call the body that the CPU backends
//! also call. Argument lists are identical to the hand-written CUDA C
//! reference (`crates/foehncraft-cuda/kernels/g1_reference.cu`).
//!
//! Needs the nightly toolchain: the `ptx-kernel` ABI and the NVPTX
//! intrinsics are not stable yet. This is the open point of D11 that gate G1
//! measures.
#![no_std]
#![feature(abi_ptx, stdarch_nvptx)]
#![allow(clippy::missing_safety_doc, clippy::too_many_arguments)]

mod amg;
mod flow;
mod fv;
mod particle;
mod radiation;
mod steam;

use foehncraft_kernels::{face, spmv};
use core::arch::nvptx::{_block_dim_x, _block_idx_x, _thread_idx_x};
use core::slice::{from_raw_parts as sl, from_raw_parts_mut as sl_mut};

#[panic_handler]
fn panic(_: &core::panic::PanicInfo) -> ! {
    // SAFETY: ends the kernel with an error the host sees on synchronisation.
    unsafe { core::arch::nvptx::trap() }
}

#[inline(always)]
pub(crate) fn tid() -> usize {
    // SAFETY: special registers, always readable.
    unsafe { (_block_idx_x() * _block_dim_x() + _thread_idx_x()) as usize }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn g1_face_flux(
    n_faces: u32,
    n_cells: u32,
    owner: *const u32,
    neighbour: *const u32,
    phi: *const f64,
    gamma: *const f64,
    mass_flux: *const f64,
    flux: *mut f64,
) {
    let f = tid();
    let nf = n_faces as usize;
    if f >= nf {
        return;
    }
    // SAFETY: the host passes device arrays of the given sizes and a mesh
    // validated on the host.
    unsafe {
        let v = face::face_flux(
            f,
            sl(owner, nf),
            sl(neighbour, nf),
            sl(phi, n_cells as usize),
            sl(gamma, nf),
            sl(mass_flux, nf),
        );
        *flux.add(f) = v;
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn g1_cell_gather(
    n_cells: u32,
    n_entries: u32,
    n_faces: u32,
    cell_ptr: *const u32,
    cell_faces: *const u32,
    flux: *const f64,
    out: *mut f64,
) {
    let c = tid();
    let nc = n_cells as usize;
    if c >= nc {
        return;
    }
    // SAFETY: as above.
    unsafe {
        let v = face::cell_net_flux(
            c,
            sl(cell_ptr, nc + 1),
            sl(cell_faces, n_entries as usize),
            sl(flux, n_faces as usize),
        );
        *out.add(c) = v;
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn g1_spmv_csr(
    n: u32,
    nnz: u32,
    row_ptr: *const u32,
    col: *const u32,
    val: *const f64,
    x: *const f64,
    y: *mut f64,
) {
    let r = tid();
    let n = n as usize;
    if r >= n {
        return;
    }
    // SAFETY: as above.
    unsafe {
        let nz = nnz as usize;
        let v = spmv::csr_row(r, sl(row_ptr, n + 1), sl(col, nz), sl(val, nz), sl(x, n));
        *y.add(r) = v;
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn g1_spmv_sell(
    n: u32,
    chunk: u32,
    n_slices: u32,
    stored: u32,
    slice_ptr: *const u32,
    slice_len: *const u32,
    col: *const u32,
    val: *const f64,
    x: *const f64,
    y_sorted: *mut f64,
) {
    let r = tid();
    let n = n as usize;
    if r >= n {
        return;
    }
    // SAFETY: as above.
    unsafe {
        let (ns, st) = (n_slices as usize, stored as usize);
        let v = spmv::sell_row(
            r,
            chunk as usize,
            sl(slice_ptr, ns),
            sl(slice_len, ns),
            sl(col, st),
            sl(val, st),
            sl(x, n),
        );
        sl_mut(y_sorted, n)[r] = v;
    }
}
