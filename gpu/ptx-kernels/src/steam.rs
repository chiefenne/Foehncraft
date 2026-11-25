//! Water and steam states of IAPWS-IF97 (stage SF7a) as the PTX entry
//! `steam_state`, one GPU thread per state. The body is
//! `foehncraft_kernels::steam::state_entry`, the one the CPU runs; the host
//! side is `crates/foehncraft-cuda/src/steam.rs`. No reductions, no
//! atomics: each thread writes the outputs of its own point.

use crate::tid;
use foehncraft_kernels::steam::{StateArgs, state_entry};

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn steam_state(a: StateArgs) {
    let i = tid();
    if (i as u64) < a.a.len && (i as u64) < a.b.len {
        // SAFETY: i lies inside `a` and `b`; the host allocates `out` with
        // OUT_WORDS entries per point.
        unsafe { state_entry(i, &a) }
    }
}
