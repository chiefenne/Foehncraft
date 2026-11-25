//! MPI step W5 (D93): the flow solver on 1 to 4 ranks gives bitwise the
//! velocity, pressure, pressure gradient, limiter, scalars, temperature,
//! turbulence variables and residual history of the undecomposed run, for
//! polyhedra, mixed cells and tetrahedra: steady SIMPLEC (linear upwind,
//! bounded central, CG with multigrid and with Jacobi for the pressure),
//! transient PISO with BDF2, a passive scalar, the energy equation, the
//! k-ω SST and the k-ε model with wall functions, the coupled solver, and
//! the species of an ideal-gas mixture with mixture properties (SF1).
//! Ranks are threads here ([`ThreadComm`]); the example `mpi_flow` runs the
//! same check across MPI processes in CI.

use foehncraft_flow::ranks::{check, check_cases, solve_one, solve_part};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::decompose::SharedComm;
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_parallel::ThreadComm;

#[test]
fn flow_on_ranks_is_bitwise_one_rank() {
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

/// MPI step W6 (D93, D84): the restart written on 4 ranks is gathered into
/// one file byte for byte that of one rank, and the run continued from it
/// on 2 ranks (and 3) is bitwise the run of one rank without the
/// interruption: steady, transient, a scalar, energy, SST and the coupled
/// solver.
#[cfg(feature = "hdf5")]
#[test]
fn restart_written_on_four_ranks_continues_on_two() {
    let dir = std::env::temp_dir().join(format!("foehn-w6-restart-{}", std::process::id()));
    for case in check_cases().unwrap() {
        let steps = case.iterations;
        let (first, second) = (steps / 2, steps - steps / 2);
        for (w, r) in [(4, 2), (2, 3)] {
            let sub = dir.join(format!(
                "{}-{w}-{r}",
                case.name.replace([' ', ',', '-'], "_")
            ));
            let c =
                foehncraft_flow::ranks::restart_check(&case, first, second, w, r, &sub).unwrap();
            assert!(
                c.same_bytes,
                "{}: file written on {w} ranks differs",
                case.name
            );
            assert!(
                c.bitwise(),
                "{}: written on {w}, read on {r}: {c:?}",
                case.name
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}
