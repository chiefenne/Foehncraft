//! Sizing spike: cost of the two communication steps of a multi-rank solver
//! step on one machine (reduction and halo message).
//!
//! ```text
//! cargo build --release -p foehncraft-parallel --features mpi --example mpi_latency
//! mpirun -n 4 target/release/examples/mpi_latency
//! ```

use foehncraft_parallel::{Communicator, MpiComm};
use std::time::Instant;

fn main() {
    let comm = MpiComm::init().expect("MPI initialisation");
    let (rank, size) = (comm.rank(), comm.size());
    let reps = 20_000;
    comm.barrier();
    let t = Instant::now();
    let mut s = 0.0;
    for i in 0..reps {
        s += comm.sum_reproducible(i as f64 + rank as f64);
    }
    let reduce_us = t.elapsed().as_secs_f64() / reps as f64 * 1e6;
    let mut line = format!("{size} ranks: reproducible sum {reduce_us:.2} us");
    for n in [100usize, 1_000, 10_000, 100_000] {
        let send = vec![1.0; n];
        let mut recv = vec![0.0; n];
        let (dest, source) = ((rank + 1) % size, (rank + size - 1) % size);
        let reps = if n >= 100_000 { 500 } else { 5_000 };
        comm.barrier();
        let t = Instant::now();
        for _ in 0..reps {
            comm.send_recv_f64(&send, dest, &mut recv, source);
        }
        let us = t.elapsed().as_secs_f64() / reps as f64 * 1e6;
        line += &format!(", ring {n} values {us:.1} us");
    }
    if rank == 0 {
        println!("{line} (check {s:.0})");
    }
}
