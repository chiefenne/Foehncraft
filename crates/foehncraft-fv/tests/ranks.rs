//! MPI step W3 (D93): the steady scalar solver on 1 to 4 ranks gives
//! bitwise the solution, gradient, limiter, matrix, right-hand side and
//! residual history of the undecomposed run, for polyhedra, mixed cells and
//! tetrahedra, linear upwind with Venkatakrishnan's limiter, the bounded
//! central scheme with the Venkatakrishnan–Wang limiter, BiCGStab and CG.
//! Ranks are threads here ([`ThreadComm`]); the example `mpi_scalar` runs
//! the same check across MPI processes in CI.

use foehncraft_fv::ranks::{check, check_cases, solve_one, solve_part};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::decompose::SharedComm;
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_parallel::ThreadComm;

#[test]
fn scalar_solves_on_ranks_are_bitwise_one_rank() {
    for case in check_cases().unwrap() {
        let reference = solve_one(&case.mesh, &case.problem, &case.settings, Some(1)).unwrap();
        assert!(
            reference.outcome.outer > 2,
            "{}: too few outer iterations",
            case.name
        );
        let geo = Geometry::compute(&case.mesh);
        for ranks in 1..=4usize {
            let part = partition(&case.mesh, &geo, ranks, Method::Rcb).unwrap();
            // two worker threads per rank: also independent of the threads
            let checks = ThreadComm::run(ranks, |c| {
                let rank = foehncraft_parallel::Communicator::rank(&*c);
                let comm: SharedComm = c;
                let mine = solve_part(
                    &case.mesh,
                    &part,
                    1,
                    comm,
                    &case.problem,
                    &case.settings,
                    Some(2),
                )
                .unwrap();
                check(rank, &reference, &mine)
            });
            let cells: usize = checks.iter().map(|c| c.cells).sum();
            assert_eq!(cells, case.mesh.n_cells, "{} on {ranks} ranks", case.name);
            for c in &checks {
                assert!(
                    c.bitwise_equal(),
                    "{} on {ranks} ranks, rank {}: {c:?}",
                    case.name,
                    c.rank
                );
            }
        }
    }
}

/// A partition that cuts the mesh into scattered pieces (cells dealt out
/// round-robin): every face a processor face, many owned by halo cells.
#[test]
fn scattered_partitions_are_bitwise_one_rank() {
    let case = &check_cases().unwrap()[0];
    let reference = solve_one(&case.mesh, &case.problem, &case.settings, Some(1)).unwrap();
    let part: Vec<u32> = (0..case.mesh.n_cells as u32).map(|c| c % 3).collect();
    let checks = ThreadComm::run(3, |c| {
        let rank = foehncraft_parallel::Communicator::rank(&*c);
        let mine = solve_part(
            &case.mesh,
            &part,
            1,
            c,
            &case.problem,
            &case.settings,
            Some(1),
        )
        .unwrap();
        check(rank, &reference, &mine)
    });
    for c in &checks {
        assert!(c.bitwise_equal(), "rank {}: {c:?}", c.rank);
    }
}
