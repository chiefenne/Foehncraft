//! MPI step W3 across processes: the steady scalar solver on the ranks of
//! a decomposed mesh gives bitwise the solution, gradient, limiter, matrix
//! and residual history of one rank (D93).
//!
//! ```text
//! cargo build --release -p foehncraft-fv --features mpi --example mpi_scalar
//! mpirun -n 3 target/release/examples/mpi_scalar
//! ```
//!
//! Every rank solves each case of `foehncraft_fv::ranks::check_cases` on
//! the whole mesh (the reference) and on its part (recursive coordinate
//! bisection, one halo layer), then compares its own cells. The run fails
//! when any bit differs on any rank.

use foehncraft_fv::ranks::{check, check_cases, solve_one, solve_part};
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
        let geo = Geometry::compute(&case.mesh);
        let part = partition(&case.mesh, &geo, size, Method::Rcb).expect("partition");
        let reference = solve_one(
            &case.mesh,
            &case.problem,
            &case.settings,
            case.reduction.as_ref(),
            Some(1),
        )
        .expect("one rank");
        let mine = solve_part(
            &case.mesh,
            &part,
            1,
            comm.clone(),
            &case.problem,
            &case.settings,
            case.reduction.as_ref(),
            Some(1),
        )
        .expect("ranks");
        let c = check(rank, &reference, &mine);
        println!(
            "rank {rank}/{size} {}: {} cells, φ {} / gradient {} / matrix {} differ, history {}",
            case.name,
            c.cells,
            c.differing_phi,
            c.differing_gradients,
            c.differing_matrix,
            if c.same_history { "same" } else { "DIFFERENT" }
        );
        if !c.bitwise_equal() {
            failures += 1.0;
        }
    }
    let total = comm.sum_reproducible(failures);
    if rank == 0 {
        println!("{size} ranks: {total} case-rank pairs differ from one rank");
    }
    assert_eq!(total, 0.0, "rank results differ from one rank");
}
