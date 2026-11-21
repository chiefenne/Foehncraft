//! hypre as the reference solver (D20, feature `hypre`): BoomerAMG as the
//! preconditioner of hypre's PCG, through the C shim `hypre/shim.c` over
//! hypre's public interface. Never needed to build or run Foehncraft.
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use crate::amg::Csr;
pub use crate::poisson::HypreSettings;
use serde::Serialize;
use std::ffi::{CStr, c_char, c_double, c_int};

#[repr(C)]
struct Options {
    strong_threshold: c_double,
    coarsen_type: c_int,
    interp_type: c_int,
    relax_type: c_int,
    num_sweeps: c_int,
    agg_num_levels: c_int,
    pmax: c_int,
    max_coarse: c_int,
    print_level: c_int,
}

#[repr(C)]
#[derive(Default)]
struct Res {
    iterations: c_int,
    status: c_int,
    final_rel_res: c_double,
    setup_seconds: c_double,
    solve_seconds: c_double,
}

unsafe extern "C" {
    fn foehncraft_hypre_uses_gpu() -> c_int;
    fn foehncraft_hypre_version() -> *const c_char;
    #[allow(clippy::too_many_arguments)]
    fn foehncraft_hypre_solve(
        n: c_int,
        row_ptr: *const c_int,
        col: *const c_int,
        val: *const c_double,
        rhs: *const c_double,
        x: *mut c_double,
        tol: c_double,
        max_iter: c_int,
        o: *const Options,
        res: *mut Res,
    ) -> c_int;
}

/// Result of a hypre solve.
#[derive(Debug, Clone, Serialize)]
pub struct HypreOutcome {
    /// PCG iterations.
    pub iterations: usize,
    /// ‖r‖ / ‖b‖ reported by hypre.
    pub final_rel_res: f64,
    /// Seconds of the setup (BoomerAMG hierarchy, on hypre's device).
    pub setup_seconds: f64,
    /// Seconds of the solve.
    pub solve_seconds: f64,
    /// The solution.
    #[serde(skip)]
    pub x: Vec<f64>,
}

/// True when the linked hypre runs on a GPU.
pub fn uses_gpu() -> bool {
    // SAFETY: no arguments, returns a constant.
    unsafe { foehncraft_hypre_uses_gpu() != 0 }
}

/// Version of the linked hypre.
pub fn version() -> String {
    // SAFETY: returns a static, NUL-terminated string.
    unsafe { CStr::from_ptr(foehncraft_hypre_version()) }
        .to_string_lossy()
        .into_owned()
}

/// Device description.
pub fn device() -> String {
    format!(
        "hypre {} ({})",
        version(),
        if uses_gpu() { "GPU" } else { "CPU, 1 thread" }
    )
}

/// Solves A x = b from x = 0 to ‖r‖ ≤ `tol` ‖b‖.
pub fn solve(
    a: &Csr,
    b: &[f64],
    tol: f64,
    max_iter: usize,
    s: &HypreSettings,
) -> anyhow::Result<HypreOutcome> {
    anyhow::ensure!(
        a.n_rows == b.len() && a.n_rows < i32::MAX as usize && a.nnz() < i32::MAX as usize,
        "hypre: matrix and right-hand side do not match"
    );
    let row_ptr: Vec<c_int> = a.row_ptr.iter().map(|&p| p as c_int).collect();
    let col: Vec<c_int> = a.col.iter().map(|&c| c as c_int).collect();
    let mut x = vec![0.0; a.n_rows];
    let o = Options {
        strong_threshold: s.strong_threshold,
        coarsen_type: s.coarsen_type,
        interp_type: s.interp_type,
        relax_type: s.relax_type.unwrap_or(if uses_gpu() { 18 } else { 6 }),
        num_sweeps: s.num_sweeps,
        agg_num_levels: s.agg_num_levels,
        pmax: s.pmax,
        max_coarse: s.max_coarse,
        print_level: s.print_level,
    };
    let mut r = Res::default();
    // SAFETY: arrays of the CSR matrix with n + 1 row pointers and nnz
    // entries, vectors of length n; the shim copies them.
    unsafe {
        foehncraft_hypre_solve(
            a.n_rows as c_int,
            row_ptr.as_ptr(),
            col.as_ptr(),
            a.val.as_ptr(),
            b.as_ptr(),
            x.as_mut_ptr(),
            tol,
            max_iter as c_int,
            &o,
            &mut r,
        );
    }
    Ok(HypreOutcome {
        iterations: r.iterations as usize,
        final_rel_res: r.final_rel_res,
        setup_seconds: r.setup_seconds,
        solve_seconds: r.solve_seconds,
        x,
    })
}
