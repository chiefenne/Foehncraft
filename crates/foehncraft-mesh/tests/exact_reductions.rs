//! MPI plan step W2 (D31 (1), D93 (4)): sums over a decomposed mesh are
//! bitwise the one-rank sums, for 1 to 9 ranks, METIS and recursive
//! coordinate bisection, four cell types, 2D and 3D. The decomposition is
//! held in one process; the accumulators of the ranks are merged through
//! their integer images, as `Communicator::sum_exact` does across MPI
//! processes (example `mpi_reductions`, in CI with 1 to 4 ranks).

use foehncraft_mesh::Geometry;
use foehncraft_mesh::decompose::{decompose, reduction_check};
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::partition::{Method, partition};

fn methods() -> Vec<Method> {
    if cfg!(feature = "metis") {
        vec![Method::Rcb, Method::Metis]
    } else {
        vec![Method::Rcb]
    }
}

#[test]
fn sums_over_ranks_are_bitwise_the_one_rank_sums() {
    let mut old_differs = 0;
    for (kind, dim, cells) in [
        (Kind::Hex, 3, [9, 7, 5]),
        (Kind::Tet, 3, [6, 5, 4]),
        (Kind::Poly, 3, [6, 5, 4]),
        (Kind::Mixed, 3, [6, 6, 4]),
        (Kind::Poly, 2, [20, 12, 1]),
    ] {
        let mut spec = BoxSpec::new(kind, dim, cells);
        spec.perturb = 0.2;
        let m = spec.generate().unwrap();
        let geo = Geometry::compute(&m);
        for method in methods() {
            for parts in [1usize, 2, 3, 4, 5, 9] {
                let part = partition(&m, &geo, parts, method).unwrap();
                let locals = decompose(&m, &part, parts, 1).unwrap();
                for r in reduction_check(&m, &geo, &locals) {
                    assert!(
                        r.bitwise_equal,
                        "{kind:?} {dim}D, {method:?}, {parts} ranks: {} {:e} against {:e}",
                        r.name, r.ranks, r.one_rank
                    );
                    if r.pairwise_ulps > 0 {
                        old_differs += 1;
                    }
                }
            }
        }
    }
    // the check must bite: the old rank-order pairwise sum moves
    assert!(
        old_differs > 0,
        "the old sum never differed: the test shows nothing"
    );
}
