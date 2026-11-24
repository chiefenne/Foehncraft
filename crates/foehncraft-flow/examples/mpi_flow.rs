//! MPI step W5 across processes: steady SIMPLEC on the ranks of a
//! decomposed mesh gives bitwise the velocity, pressure, pressure gradient,
//! limiter and residual history of one rank (D93).
//!
//! ```text
//! cargo build --release -p foehncraft-flow --features mpi --example mpi_flow
//! mpirun -n 3 target/release/examples/mpi_flow
//! ```
//!
//! Every rank solves each case of `foehncraft_flow::ranks::check_cases` on
//! the whole mesh (the reference) and on its part (recursive coordinate
//! bisection, one halo layer), then compares its own cells. The run fails
//! when any bit differs on any rank.

use foehncraft_flow::ranks::{check, check_cases, solve_one, solve_part};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::decompose::SharedComm;
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_parallel::MpiComm;
use std::sync::Arc;

fn main() {
    let comm: SharedComm = Arc::new(MpiComm::init().expect("MPI initialisation"));
    let (rank, size) = (comm.rank(), comm.size());
    let mut failures = 0.0;
    for case in check_cases().expect("check cases") {
        let mesh = case.mesh().expect("mesh");
        let geo = Geometry::compute(&mesh);
        let part = partition(&mesh, &geo, size, Method::Rcb).expect("partition");
        let reference = solve_one(&case.case, &mesh, case.iterations, Some(1)).expect("one rank");
        let mine = solve_part(
            &case.case,
            &mesh,
            &part,
            comm.clone(),
            case.iterations,
            Some(1),
        )
        .expect("ranks");
        let c = check(rank, &mine, &reference);
        println!(
            "rank {rank}/{size} {}: {} cells, fields {} / gradients {} differ, history {}",
            case.name,
            c.cells,
            c.differing_fields,
            c.differing_gradients,
            if c.same_history { "same" } else { "DIFFERENT" }
        );
        if !c.bitwise() {
            failures += 1.0;
        }
    }
    let failures = comm.sum_reproducible(failures);
    if rank == 0 {
        println!("{size} ranks: {failures} case-rank pairs differ from one rank");
    }
    if failures > 0.0 {
        std::process::exit(1);
    }
}
