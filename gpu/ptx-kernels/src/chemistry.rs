//! Gas-phase chemistry in the cells of a flow (steam furnace SF2b): the
//! reactor of one cell over a time step as the PTX entry `chemistry_react`,
//! one GPU thread per cell. The body is
//! `foehncraft_kernels::chemistry::react_entry`, the one the CPU threads
//! run (through `react_cell`); the host side is
//! `crates/foehncraft-cuda/src/chemistry.rs`.

use crate::tid;
use foehncraft_kernels::chemistry::{ReactArgs, react_entry};

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn chemistry_react(a: ReactArgs) {
    let i = tid();
    if (i as u64) < a.n_cells {
        // SAFETY: i < n_cells; the host sized state, work and stats for
        // n_cells rows, and each thread writes only its own rows.
        unsafe { react_entry(i, &a) }
    }
}
