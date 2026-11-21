//! CPU backends.
//!
//! - `cpu-serial`: the single-source kernel bodies in a sequential loop;
//! - `cpu-threads`: the same bodies in a multi-threaded loop (rayon), the
//!   CPU half of the single-source route;
//! - `cpu-native`: hand-written, idiomatic multi-threaded Rust loops that do
//!   not go through the bodies; the CPU reference for G1.
//!
//! Each output value is written by exactly one iteration, so all three give
//! bitwise the same results for any number of threads (D31).
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use super::Backend;
use crate::workload::{Kernel, Workload};
use foehncraft_kernels::{face, spmv, unpack_face};
use rayon::prelude::*;
use std::time::Instant;

/// How the CPU backend runs the kernels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// Kernel bodies, one thread.
    Serial,
    /// Kernel bodies, all threads.
    Threads,
    /// Hand-written loops, all threads.
    Native,
}

/// CPU backend.
pub struct CpuBackend {
    mode: Mode,
    pool: rayon::ThreadPool,
    out: Vec<f64>,
}

impl CpuBackend {
    /// New backend; `threads = None` uses all cores.
    pub fn new(mode: Mode, threads: Option<usize>) -> anyhow::Result<Self> {
        let n = match mode {
            Mode::Serial => 1,
            _ => threads.unwrap_or(0),
        };
        let pool = rayon::ThreadPoolBuilder::new().num_threads(n).build()?;
        Ok(Self {
            mode,
            pool,
            out: Vec::new(),
        })
    }
}

/// Sequential or parallel loop over the output, by mode.
fn for_each_out(mode: Mode, out: &mut [f64], body: impl Fn(usize) -> f64 + Sync) {
    const GRAIN: usize = 1 << 16;
    match mode {
        Mode::Serial => out.iter_mut().enumerate().for_each(|(i, o)| *o = body(i)),
        _ => out
            .par_chunks_mut(GRAIN)
            .enumerate()
            .for_each(|(c, chunk)| {
                for (j, o) in chunk.iter_mut().enumerate() {
                    *o = body(c * GRAIN + j);
                }
            }),
    }
}

impl Backend for CpuBackend {
    fn name(&self) -> &'static str {
        match self.mode {
            Mode::Serial => "cpu-serial",
            Mode::Threads => "cpu-threads",
            Mode::Native => "cpu-native",
        }
    }

    fn device(&self) -> String {
        format!(
            "{} ({} threads)",
            crate::sysinfo::cpu_model(),
            self.pool.current_num_threads()
        )
    }

    fn setup(&mut self, _w: &Workload) -> anyhow::Result<()> {
        Ok(())
    }

    fn run(&mut self, w: &Workload, k: Kernel) -> anyhow::Result<f64> {
        let n = w.output_len(k);
        let mut out = std::mem::take(&mut self.out);
        out.resize(n, 0.0);
        let mode = self.mode;
        let t0 = Instant::now();
        self.pool.install(|| {
            if mode == Mode::Native {
                native(k, w, &mut out);
                return;
            }
            let m = &w.mesh;
            // SAFETY (all bodies): the mesh and matrices were validated when
            // the workload was built, and `i` ranges over the outputs.
            match k {
                Kernel::FaceFlux => for_each_out(mode, &mut out, |f| {
                    // SAFETY: as above.
                    unsafe {
                        face::face_flux(f, &m.owner, &m.neighbour, &w.phi, &w.gamma, &w.mass_flux)
                    }
                }),
                Kernel::CellGather => for_each_out(mode, &mut out, |c| {
                    // SAFETY: as above.
                    unsafe { face::cell_net_flux(c, &m.cell_ptr, &m.cell_faces, &w.flux) }
                }),
                Kernel::SpmvCsr => {
                    let a = &w.csr;
                    for_each_out(mode, &mut out, |r| {
                        // SAFETY: as above.
                        unsafe { spmv::csr_row(r, &a.row_ptr, &a.col, &a.val, &w.x) }
                    })
                }
                Kernel::SpmvSell => {
                    // Results in sorted order; scattered back by `perm` after
                    // timing when σ > 1 (identity otherwise).
                    let a = &w.sell;
                    for_each_out(mode, &mut out, |r| {
                        // SAFETY: as above.
                        unsafe {
                            spmv::sell_row(
                                r,
                                a.chunk,
                                &a.slice_ptr,
                                &a.slice_len,
                                &a.col,
                                &a.val,
                                &w.x,
                            )
                        }
                    })
                }
            }
        });
        let dt = t0.elapsed().as_secs_f64();
        self.out = out;
        Ok(dt)
    }

    fn output(&mut self, w: &Workload, k: Kernel) -> anyhow::Result<Vec<f64>> {
        if k == Kernel::SpmvSell {
            let mut y = vec![0.0; w.sell.n];
            for (r, &v) in self.out.iter().enumerate() {
                y[w.sell.perm[r] as usize] = v;
            }
            return Ok(y);
        }
        Ok(self.out.clone())
    }
}

/// Hand-written loops: plain safe Rust over slices and iterators.
fn native(k: Kernel, w: &Workload, out: &mut [f64]) {
    let m = &w.mesh;
    match k {
        Kernel::FaceFlux => out
            .par_iter_mut()
            .zip(m.owner.par_iter().zip(&m.neighbour))
            .zip(w.gamma.par_iter().zip(&w.mass_flux))
            .for_each(|((o, (&p, &n)), (&g, &mf))| {
                let (pp, pn) = (w.phi[p as usize], w.phi[n as usize]);
                *o = mf.max(0.0) * pp + mf.min(0.0) * pn - g * (pn - pp);
            }),
        Kernel::CellGather => {
            out.par_iter_mut()
                .zip(m.cell_ptr.par_windows(2))
                .for_each(|(o, r)| {
                    let mut s = 0.0;
                    for &e in &m.cell_faces[r[0] as usize..r[1] as usize] {
                        let (f, is_nb) = unpack_face(e);
                        if is_nb {
                            s -= w.flux[f];
                        } else {
                            s += w.flux[f];
                        }
                    }
                    *o = s;
                })
        }
        Kernel::SpmvCsr => {
            let a = &w.csr;
            out.par_iter_mut()
                .zip(a.row_ptr.par_windows(2))
                .for_each(|(o, r)| {
                    let (lo, hi) = (r[0] as usize, r[1] as usize);
                    let mut s = 0.0;
                    for (&j, &v) in a.col[lo..hi].iter().zip(&a.val[lo..hi]) {
                        s += v * w.x[j as usize];
                    }
                    *o = s;
                })
        }
        Kernel::SpmvSell => {
            // One slice per task, results in sorted row order (scattered
            // back by `perm` in `output`, as for the other modes).
            let a = &w.sell;
            let c = a.chunk;
            out.par_chunks_mut(c).enumerate().for_each(|(s, ys)| {
                let base = a.slice_ptr[s] as usize;
                let len = a.slice_len[s] as usize;
                for (lane, y) in ys.iter_mut().enumerate() {
                    let mut acc = 0.0;
                    for j in 0..len {
                        let k = base + j * c + lane;
                        acc += a.val[k] * w.x[a.col[k] as usize];
                    }
                    *y = acc;
                }
            });
        }
    }
}
