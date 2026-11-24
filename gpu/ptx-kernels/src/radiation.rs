//! Radiation R1a (step 29, D259): the transport of the energy bundles as
//! the PTX entry `radiation_transport`, one GPU thread per bundle. The body
//! is `foehncraft_kernels::radiation::transport_entry`, the one the CPU
//! runs; the host side is `crates/foehncraft-cuda/src/radiation.rs`.
//!
//! The tallies are whole numbers of quanta added with 64-bit integer
//! atomics: an integer sum is the same in any order, so the result is
//! bitwise that of the CPU (the one exception of R10, D269). No
//! floating-point atomics.

use crate::tid;
use core::sync::atomic::{AtomicU64, Ordering};
use foehncraft_kernels::entry::Raw;
use foehncraft_kernels::radiation::{Tally, TransportArgs, transport_entry};

/// The tallies of the argument block, added to with integer atomics.
struct Atomics<'a> {
    a: &'a TransportArgs,
}

#[inline(always)]
fn add(r: Raw<u64>, i: usize, q: i64) {
    // SAFETY: the host allocates every tally array with the size of the
    // argument block (batches × cells or faces, three ledger entries) and
    // the transport passes indices inside it; AtomicU64 has the layout of u64.
    unsafe { AtomicU64::from_ptr(r.out().add(i)).fetch_add(q as u64, Ordering::Relaxed) };
}

impl Tally for Atomics<'_> {
    fn cell(&self, batch: u32, cell: usize, q: i64) {
        add(self.a.cell, batch as usize * self.a.n_cells as usize + cell, q);
    }

    fn wall(&self, batch: u32, face: usize, incident: i64, absorbed: i64) {
        let i = batch as usize * self.a.n_boundary as usize + face;
        add(self.a.incident, i, incident);
        add(self.a.absorbed, i, absorbed);
    }

    fn escaped(&self, q: i64) {
        add(self.a.ledger, 0, q);
    }

    fn roulette(&self, q: i64) {
        add(self.a.ledger, 1, q);
    }

    fn lost(&self, q: i64) {
        add(self.a.ledger, 2, q);
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn radiation_transport(a: TransportArgs) {
    let i = tid();
    if (i as u64) < a.n_bundles {
        // SAFETY: i < n_bundles; the host built every array of `a` for one
        // mesh (one GPU, no ranks: a bundle never reaches a processor face).
        unsafe {
            transport_entry(i, &a, &Atomics { a: &a });
        }
    }
}
