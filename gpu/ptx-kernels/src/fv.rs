//! Stage 2: the finite-volume operators and the Krylov vector operations of
//! `foehncraft_kernels::{fv, blas}` as PTX entry points, one GPU thread per
//! output. The host side is `crates/foehncraft-cuda/src/fv.rs`; the
//! argument lists there must match these signatures exactly.
//!
//! Every pointer comes with its length (or a length from which it follows);
//! the host validated the mesh and built every array for it.

use crate::tid;
use foehncraft_kernels::fv::{self, Boundary, FaceInputs, SellLayout, Topology, Vec3};
use foehncraft_kernels::{blas, spmv};
use core::slice::from_raw_parts as sl;

/// Topology from raw parts.
///
/// # Safety
/// Device arrays of the given lengths.
#[inline(always)]
unsafe fn topo<'a>(
    nc: usize,
    nf: usize,
    ni: usize,
    ne: usize,
    owner: *const u32,
    neighbour: *const u32,
    cell_ptr: *const u32,
    cell_faces: *const u32,
) -> Topology<'a> {
    // SAFETY: guaranteed by the caller.
    unsafe {
        Topology {
            owner: sl(owner, nf),
            neighbour: sl(neighbour, ni),
            cell_ptr: sl(cell_ptr, nc + 1),
            cell_faces: sl(cell_faces, ne),
            n_interior: ni,
            periodic: ni,
        }
    }
}

/// # Safety
/// Device arrays of length `n` each.
#[inline(always)]
unsafe fn v3<'a>(n: usize, x: *const f64, y: *const f64, z: *const f64) -> Vec3<'a, f64> {
    // SAFETY: guaranteed by the caller.
    unsafe {
        Vec3 {
            x: sl(x, n),
            y: sl(y, n),
            z: sl(z, n),
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn fv_gradient(
    n_cells: u32,
    n_faces: u32,
    n_interior: u32,
    n_entries: u32,
    owner: *const u32,
    neighbour: *const u32,
    cell_ptr: *const u32,
    cell_faces: *const u32,
    wx: *const f64,
    wy: *const f64,
    wz: *const f64,
    phi: *const f64,
    bc_a: *const f64,
    bc_b: *const f64,
    gx: *mut f64,
    gy: *mut f64,
    gz: *mut f64,
) {
    let c = tid();
    let (nc, nf, ni, ne) = (
        n_cells as usize,
        n_faces as usize,
        n_interior as usize,
        n_entries as usize,
    );
    if c >= nc {
        return;
    }
    // SAFETY: see the module docs.
    unsafe {
        let t = topo(nc, nf, ni, ne, owner, neighbour, cell_ptr, cell_faces);
        let bc = Boundary {
            a: sl(bc_a, nf - ni),
            b: sl(bc_b, nf - ni),
            jump: &[],
        };
        let g = fv::lsq_gradient(c, &t, &v3(ne, wx, wy, wz), sl(phi, nc), &bc);
        *gx.add(c) = g[0];
        *gy.add(c) = g[1];
        *gz.add(c) = g[2];
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn fv_limiter(
    n_cells: u32,
    n_faces: u32,
    n_interior: u32,
    n_entries: u32,
    n_eps2: u32,
    owner: *const u32,
    neighbour: *const u32,
    cell_ptr: *const u32,
    cell_faces: *const u32,
    phi: *const f64,
    bc_a: *const f64,
    bc_b: *const f64,
    gx: *const f64,
    gy: *const f64,
    gz: *const f64,
    ccx: *const f64,
    ccy: *const f64,
    ccz: *const f64,
    fcx: *const f64,
    fcy: *const f64,
    fcz: *const f64,
    eps2: *const f64,
    psi: *mut f64,
) {
    let c = tid();
    let (nc, nf, ni, ne) = (
        n_cells as usize,
        n_faces as usize,
        n_interior as usize,
        n_entries as usize,
    );
    if c >= nc {
        return;
    }
    // SAFETY: see the module docs.
    unsafe {
        let t = topo(nc, nf, ni, ne, owner, neighbour, cell_ptr, cell_faces);
        let bc = Boundary {
            a: sl(bc_a, nf - ni),
            b: sl(bc_b, nf - ni),
            jump: &[],
        };
        let v = fv::venkat_limiter(
            c,
            &t,
            sl(phi, nc),
            &bc,
            &v3(nc, gx, gy, gz),
            &v3(nc, ccx, ccy, ccz),
            &v3(nf, fcx, fcy, fcz),
            sl(eps2, n_eps2 as usize),
        );
        *psi.add(c) = v;
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn fv_faces(
    n_cells: u32,
    n_faces: u32,
    n_interior: u32,
    blend: f64,
    central: f64,
    beta: f64,
    nonorth: f64,
    nonorth_limit: f64,
    owner: *const u32,
    neighbour: *const u32,
    bc_a: *const f64,
    bc_b: *const f64,
    mass_flux: *const f64,
    gamma: *const f64,
    dgeo: *const f64,
    kx: *const f64,
    ky: *const f64,
    kz: *const f64,
    weight: *const f64,
    fcx: *const f64,
    fcy: *const f64,
    fcz: *const f64,
    ccx: *const f64,
    ccy: *const f64,
    ccz: *const f64,
    phi: *const f64,
    gx: *const f64,
    gy: *const f64,
    gz: *const f64,
    psi: *const f64,
    cp: *mut f64,
    cn: *mut f64,
    bf: *mut f64,
) {
    let f = tid();
    let (nc, nf, ni) = (n_cells as usize, n_faces as usize, n_interior as usize);
    if f >= nf {
        return;
    }
    // SAFETY: see the module docs; the face kernel does not read the
    // cell-to-face lists.
    unsafe {
        let inp = FaceInputs {
            topo: Topology {
                owner: sl(owner, nf),
                neighbour: sl(neighbour, ni),
                cell_ptr: &[],
                cell_faces: &[],
                n_interior: ni,
                periodic: ni,
            },
            bc: Boundary {
                a: sl(bc_a, nf - ni),
                b: sl(bc_b, nf - ni),
                jump: &[],
            },
            mass_flux: sl(mass_flux, nf),
            gamma: sl(gamma, nf),
            dgeo: sl(dgeo, nf),
            k: v3(nf, kx, ky, kz),
            weight: sl(weight, nf),
            face_centroid: v3(nf, fcx, fcy, fcz),
            cell_centroid: v3(nc, ccx, ccy, ccz),
            shift: Vec3 {
                x: &[],
                y: &[],
                z: &[],
            },
            phi: sl(phi, nc),
            grad: v3(nc, gx, gy, gz),
            psi: sl(psi, nc),
            blend,
            central,
            beta,
            nonorth,
            nonorth_limit,
            bound: 0.0,
        };
        let (a, b, c) = fv::face_coefficients(f, &inp);
        *cp.add(f) = a;
        *cn.add(f) = b;
        *bf.add(f) = c;
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn fv_assemble(
    n_cells: u32,
    n_faces: u32,
    n_interior: u32,
    n_entries: u32,
    chunk: u32,
    n_slices: u32,
    alpha: f64,
    owner: *const u32,
    neighbour: *const u32,
    cell_ptr: *const u32,
    cell_faces: *const u32,
    cp: *const f64,
    cn: *const f64,
    bf: *const f64,
    source: *const f64,
    phi: *const f64,
    slice_ptr: *const u32,
    val: *mut f64,
    rhs: *mut f64,
    dinv: *mut f64,
) {
    let c = tid();
    let (nc, nf, ni, ne) = (
        n_cells as usize,
        n_faces as usize,
        n_interior as usize,
        n_entries as usize,
    );
    if c >= nc {
        return;
    }
    // SAFETY: see the module docs; row c is written by this thread only.
    unsafe {
        let t = topo(nc, nf, ni, ne, owner, neighbour, cell_ptr, cell_faces);
        let layout = SellLayout {
            chunk: chunk as usize,
            slice_ptr: sl(slice_ptr, n_slices as usize),
        };
        fv::assemble_row(
            c,
            &t,
            sl(cp, nf),
            sl(cn, nf),
            sl(bf, nf),
            sl(source, nc),
            sl(phi, nc),
            alpha,
            &layout,
            val,
            rhs,
            dinv,
        );
    }
}

/// out = b − A x (sliced ELL, rows in cell order).
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn fv_residual(
    n: u32,
    chunk: u32,
    n_slices: u32,
    stored: u32,
    slice_ptr: *const u32,
    slice_len: *const u32,
    col: *const u32,
    val: *const f64,
    x: *const f64,
    b: *const f64,
    out: *mut f64,
) {
    let r = tid();
    let n = n as usize;
    if r >= n {
        return;
    }
    // SAFETY: see the module docs.
    unsafe {
        let (ns, st) = (n_slices as usize, stored as usize);
        let ax = spmv::sell_row(
            r,
            chunk as usize,
            sl(slice_ptr, ns),
            sl(slice_len, ns),
            sl(col, st),
            sl(val, st),
            sl(x, n),
        );
        let v = *b.add(r) - ax;
        *out.add(r) = v;
    }
}

/// Block sums of a·b: one thread per block of `foehncraft_core::repro::BLOCK`
/// values, in the same order as on the CPU.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn blas_block_dot(
    n: u32,
    n_blocks: u32,
    a: *const f64,
    b: *const f64,
    sums: *mut f64,
) {
    let k = tid();
    if k >= n_blocks as usize {
        return;
    }
    // SAFETY: see the module docs.
    unsafe {
        let v = blas::block_dot(k, sl(a, n as usize), sl(b, n as usize));
        *sums.add(k) = v;
    }
}

/// Block minima and maxima of a.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn blas_block_min_max(
    n: u32,
    n_blocks: u32,
    a: *const f64,
    mins: *mut f64,
    maxs: *mut f64,
) {
    let k = tid();
    if k >= n_blocks as usize {
        return;
    }
    // SAFETY: see the module docs.
    unsafe {
        let (lo, hi) = blas::block_min_max(k, sl(a, n as usize));
        *mins.add(k) = lo;
        *maxs.add(k) = hi;
    }
}

/// p = r + β (p − ω v), in place.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn blas_update_p(
    n: u32,
    beta: f64,
    omega: f64,
    r: *const f64,
    v: *const f64,
    p: *mut f64,
) {
    let i = tid();
    let n = n as usize;
    if i >= n {
        return;
    }
    // SAFETY: see the module docs; the shared view of p ends before the
    // write.
    unsafe {
        let x = blas::update_p(i, sl(r, n), sl(p, n), sl(v, n), beta, omega);
        *p.add(i) = x;
    }
}

/// out = x − α z (out is not x or z).
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn blas_sub_scaled(
    n: u32,
    alpha: f64,
    x: *const f64,
    z: *const f64,
    out: *mut f64,
) {
    let i = tid();
    let n = n as usize;
    if i >= n {
        return;
    }
    // SAFETY: see the module docs.
    unsafe {
        let y = blas::sub_scaled(i, sl(x, n), sl(z, n), alpha);
        *out.add(i) = y;
    }
}

/// x += α p̂ + ω ŝ, in place.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn blas_update_x(
    n: u32,
    alpha: f64,
    omega: f64,
    ph: *const f64,
    sh: *const f64,
    x: *mut f64,
) {
    let i = tid();
    let n = n as usize;
    if i >= n {
        return;
    }
    // SAFETY: see the module docs; the shared view of x ends before the
    // write.
    unsafe {
        let y = blas::update_x(i, sl(x, n), sl(ph, n), sl(sh, n), alpha, omega);
        *x.add(i) = y;
    }
}

/// out = d x (Jacobi; out is not d or x).
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn blas_scale(n: u32, d: *const f64, x: *const f64, out: *mut f64) {
    let i = tid();
    let n = n as usize;
    if i >= n {
        return;
    }
    // SAFETY: see the module docs.
    unsafe {
        let y = blas::scale(i, sl(d, n), sl(x, n));
        *out.add(i) = y;
    }
}

/// CG update y = y + α x, in place.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn blas_axpy(n: u32, alpha: f64, x: *const f64, y: *mut f64) {
    let i = tid();
    let n = n as usize;
    if i >= n {
        return;
    }
    // SAFETY: see the module docs; y is read at index i only.
    unsafe {
        let v = blas::axpy(i, alpha, sl(x, n), core::slice::from_raw_parts(y, n));
        *y.add(i) = v;
    }
}

/// CG direction update y = x + β y, in place.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn blas_xpby(n: u32, x: *const f64, beta: f64, y: *mut f64) {
    let i = tid();
    let n = n as usize;
    if i >= n {
        return;
    }
    // SAFETY: as above.
    unsafe {
        let v = blas::xpby(i, sl(x, n), beta, core::slice::from_raw_parts(y, n));
        *y.add(i) = v;
    }
}
