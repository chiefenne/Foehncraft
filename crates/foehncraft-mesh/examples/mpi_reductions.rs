//! Multi-rank MPI step W2 (D31 (1), D93 (4)): exact reductions give the
//! same bits on any number of ranks.
//!
//! ```text
//! cargo build --release -p foehncraft-mesh --features mpi --example mpi_reductions
//! mpirun -n 1 target/release/examples/mpi_reductions
//! mpirun -n 4 target/release/examples/mpi_reductions
//! ```
//!
//! Every rank builds the global mesh, partitions it (METIS, or RCB without
//! the `metis` feature), extracts its part with one halo layer and computes
//! the canonical geometry of step W1. On its own cells it then evaluates the
//! reductions a solver needs, with the threaded exact sums of
//! `foehncraft_parallel::repro` and one exact all-reduce each:
//!
//! 1. the volume of the domain Σ V,
//! 2. the volume integral Σ φ V of a field with values over 22 decades,
//! 3. a dot product (φ, ψ),
//! 4. the 2-norm of the residual r = b − A x of a diffusion operator, with x
//!    in the halo cells from the halo exchange (the owner's values),
//! 5. a boundary flux Σ φ_f |S_f| over the boundary faces,
//! 6. min and max of φ (which holds −0 and +0),
//! 7. a three-component sum in one all-reduce (the force of a wall).
//!
//! Each value is compared bit for bit with the same reduction over the
//! whole mesh on one rank, which every rank computes as the reference. The
//! old reduction (pairwise sums per rank, combined in rank order) is printed
//! next to it: it changes with the rank count, the exact one must not. The
//! example exits with status 1 when any exact value differs.

use foehncraft_core::exact::{ExactSum, max_total, min_total};
use foehncraft_core::repro::pairwise_sum;
use foehncraft_core::rng::uniform;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::Mesh;
use foehncraft_mesh::decompose::extract;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::geometry::{dot, norm, sub};
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_parallel::repro::{par_exact_dot_acc, par_exact_sum_acc};
use foehncraft_parallel::{Communicator, MpiComm, SerialComm};
use std::collections::HashSet;

/// A field keyed by global id, with magnitudes over 22 decades and signed
/// zeros, so that the order of a sum matters in the last bits.
fn phi(id: u64) -> f64 {
    if id.is_multiple_of(97) {
        return if id.is_multiple_of(2) { 0.0 } else { -0.0 };
    }
    let u = uniform([id, 11, 3]) - 0.5;
    u * 10f64.powi((id % 23) as i32 - 11)
}

fn psi(id: u64) -> f64 {
    uniform([id, 12, 3]) * 2.0 - 0.7
}

/// The reductions of one rank's part, combined over `comm`, plus the old
/// rank-order pairwise value of the field integral.
struct Values {
    exact: Vec<f64>,
    old_integral: f64,
}

const NAMES: [&str; 9] = [
    "volume",
    "integral phi V",
    "dot (phi, psi)",
    "residual norm",
    "boundary flux",
    "min phi",
    "max phi",
    "force x",
    "force y",
];

/// The reductions on the own cells `0..n_own` of `mesh` (cells beyond are
/// halo cells, whose `x` came from the exchange) and the boundary faces
/// whose global id is in `boundary`.
fn reduce(
    comm: &dyn Communicator,
    mesh: &Mesh,
    geo: &Geometry,
    n_own: usize,
    x: &[f64],
    boundary: &HashSet<u64>,
) -> Values {
    let ids = &mesh.cell_id[..n_own];
    let vol = &geo.cell_volume[..n_own];
    let f: Vec<f64> = ids.iter().map(|&g| phi(g)).collect();
    let g: Vec<f64> = ids.iter().map(|&g| psi(g)).collect();
    let fv: Vec<f64> = f.iter().zip(vol).map(|(a, v)| a * v).collect();
    // diffusion residual: each cell sums its faces in the cell's face order
    // (global face order since step W1)
    let ni = mesh.neighbour.len();
    let r: Vec<f64> = (0..n_own)
        .map(|c| {
            let mut ax = 0.0;
            for (fc, _) in mesh.cell_face_iter(c) {
                let s = geo.face_area[fc];
                let p = mesh.owner[fc] as usize;
                let (d, other) = if fc < ni {
                    let q = mesh.neighbour[fc] as usize;
                    let o = if p == c { q } else { p };
                    (sub(geo.cell_centroid[q], geo.cell_centroid[p]), x[o])
                } else {
                    (
                        sub(geo.face_centroid[fc], geo.cell_centroid[p]),
                        psi(mesh.face_id[fc]),
                    )
                };
                let a = dot(s, s) / dot(d, s).abs();
                ax += a * (x[c] - other);
            }
            vol[c] * psi(mesh.cell_id[c] ^ 0x5555) - ax
        })
        .collect();
    let own_boundary: Vec<usize> = (ni..mesh.n_faces())
        .filter(|&fc| (mesh.owner[fc] as usize) < n_own && boundary.contains(&mesh.face_id[fc]))
        .collect();
    let mut flux = ExactSum::new();
    let mut force = [ExactSum::new(), ExactSum::new()];
    for &fc in &own_boundary {
        let s = geo.face_area[fc];
        let p = phi(mesh.face_id[fc]);
        flux.add(p * norm(s));
        force[0].add(p * s[0]);
        force[1].add(p * s[1]);
    }
    let sums = comm.sum_exact_many(&[
        par_exact_sum_acc(vol),
        par_exact_sum_acc(&fv),
        par_exact_dot_acc(&f, &g),
        par_exact_dot_acc(&r, &r),
        flux,
        force[0].clone(),
        force[1].clone(),
    ]);
    let (mn, mx) = f
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
            (min_total(a, v), max_total(b, v))
        });
    Values {
        exact: vec![
            sums[0],
            sums[1],
            sums[2],
            sums[3].sqrt(),
            sums[4],
            comm.min(mn),
            comm.max(mx),
            sums[5],
            sums[6],
        ],
        old_integral: pairwise_sum(&comm.all_gather_f64(pairwise_sum(&fv))),
    }
}

fn main() {
    let comm = MpiComm::init().expect("MPI initialisation");
    let (rank, size) = (comm.rank(), comm.size());
    let kind = match std::env::var("SPIKE_KIND").as_deref() {
        Ok("hex") => Kind::Hex,
        Ok("tet") => Kind::Tet,
        Ok("mixed") => Kind::Mixed,
        _ => Kind::Poly,
    };
    let mut spec = BoxSpec::new(kind, 3, [12, 10, 8]);
    spec.perturb = 0.2;
    let global = spec.generate().expect("mesh");
    let ggeo = Geometry::compute(&global);
    let boundary: HashSet<u64> = global.face_id[global.neighbour.len()..]
        .iter()
        .copied()
        .collect();
    // the reference: the whole mesh on one rank
    let xg: Vec<f64> = global.cell_id.iter().map(|&g| psi(g ^ 0xaaaa)).collect();
    let reference = reduce(&SerialComm, &global, &ggeo, global.n_cells, &xg, &boundary);

    let method = if cfg!(feature = "metis") {
        Method::Metis
    } else {
        Method::Rcb
    };
    let part = partition(&global, &ggeo, size, method).expect("partition");
    let local = extract(&global, &part, rank, 1).expect("local mesh");
    let lgeo = local.geometry(&comm);
    let n_own = local.mesh.n_cells;
    // x on the own cells, the halo filled by the owners (D93 (3))
    let mut x: Vec<f64> = local.mesh.cell_id.iter().map(|_| f64::NAN).collect();
    for (xc, &id) in x.iter_mut().zip(&local.mesh.cell_id[..n_own]) {
        *xc = psi(id ^ 0xaaaa);
    }
    local.halo.exchange(&comm, &mut x);
    let v = reduce(&comm, &local.mesh, &lgeo, n_own, &x, &boundary);

    let mut differ = 0usize;
    if rank == 0 {
        println!("kind {kind:?}, {size} ranks, {} cells", global.n_cells);
    }
    for (i, name) in NAMES.iter().enumerate() {
        let (a, b) = (v.exact[i], reference.exact[i]);
        let same = a.to_bits() == b.to_bits();
        if !same {
            differ += 1;
        }
        if rank == 0 {
            println!(
                "  {name:<16} {a:+.17e}  {}",
                if same {
                    "bitwise = 1 rank"
                } else {
                    "DIFFERS from 1 rank"
                }
            );
        }
    }
    if rank == 0 {
        let old_same = v.old_integral.to_bits() == reference.old_integral.to_bits();
        println!(
            "  old pairwise integral {:+.17e}  {} (rank-order sum, for comparison)",
            v.old_integral,
            if old_same {
                "= 1 rank"
            } else {
                "differs (local cell order, rank order)"
            }
        );
    }
    let total = comm.sum_reproducible(differ as f64);
    if rank == 0 {
        println!(
            "exact reductions on {size} ranks: {}",
            if total == 0.0 {
                "BITWISE EQUAL TO ONE RANK"
            } else {
                "NOT bitwise equal"
            }
        );
    }
    if total != 0.0 {
        std::process::exit(1);
    }
}
