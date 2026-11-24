//! MPI step W5 (D93): steady SIMPLEC on 1 to 4 ranks gives bitwise the
//! velocity, pressure, pressure gradient, limiter and residual history of
//! the undecomposed run, for polyhedra, mixed cells and tetrahedra, linear
//! upwind and bounded central, CG with multigrid and with Jacobi for the
//! pressure. Ranks are threads here ([`ThreadComm`]); the example
//! `mpi_flow` runs the same check across MPI processes in CI.

use foehncraft_flow::ranks::{check, check_cases, solve_one, solve_part};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::decompose::SharedComm;
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_parallel::ThreadComm;

#[test]
fn steady_flow_on_ranks_is_bitwise_one_rank() {
    for case in check_cases().unwrap() {
        let mesh = case.mesh().unwrap();
        let reference = solve_one(&case.case, &mesh, case.iterations, Some(1)).unwrap();
        let geo = Geometry::compute(&mesh);
        for ranks in 1..=4usize {
            let part = partition(&mesh, &geo, ranks, Method::Rcb).unwrap();
            // two worker threads per rank: also independent of the threads
            let checks = ThreadComm::run(ranks, |c| {
                let rank = foehncraft_parallel::Communicator::rank(&*c);
                let comm: SharedComm = c;
                let mine =
                    solve_part(&case.case, &mesh, &part, comm, case.iterations, Some(2)).unwrap();
                check(rank, &mine, &reference)
            });
            for c in checks {
                assert!(c.bitwise(), "{} on {ranks} ranks: {c:?}", case.name);
            }
        }
    }
}
