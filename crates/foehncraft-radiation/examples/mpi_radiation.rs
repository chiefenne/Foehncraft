//! Radiation R1a (D258): a radiation call on MPI ranks gives bitwise the
//! tallies of one rank.
//!
//! ```text
//! cargo build --release -p foehncraft-radiation --features mpi --example mpi_radiation
//! mpirun -n 1 target/release/examples/mpi_radiation
//! mpirun -n 4 target/release/examples/mpi_radiation
//! RAYON_NUM_THREADS=1 RADIATION_BUNDLES=2000000 mpirun -n 4 target/release/examples/mpi_radiation
//! ```
//!
//! `RADIATION_BUNDLES` sets the number of bundles (200 000 by default); the
//! example prints the time of the call on the ranks and on one rank.
//!
//! Every rank builds the global mesh (polyhedra, perturbed and turned),
//! with a temperature and an absorption coefficient that vary in space,
//! grey walls at different temperatures, a black wall, a symmetry plane
//! and an open boundary. It solves the call on the whole mesh as the
//! reference, then on its part of the partitioned mesh, where bundles that
//! reach a processor face are handed to the rank behind it. The absorbed
//! quanta per cell and per wall face are summed over the ranks with an
//! integer all-reduce and compared bit for bit with the reference, and so
//! are the energy accounts. The example exits with status 1 when anything
//! differs.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::decompose::extract;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_parallel::{Communicator, MpiComm};
use foehncraft_radiation::{Domain, Settings, Wall};
use std::collections::HashMap;
use std::time::Instant;

fn wall(name: &str) -> Wall {
    match name {
        "xmin" => Wall::grey(0.8, 1200.0),
        "xmax" => Wall::grey(0.4, 500.0),
        "ymin" => Wall::black(300.0),
        "ymax" => Wall::mirror(),
        "zmin" => Wall::grey(0.6, 900.0),
        _ => Wall::open(),
    }
}

fn main() {
    let comm = MpiComm::init().expect("MPI initialisation");
    let (rank, size) = (comm.rank(), comm.size());
    let mut spec = BoxSpec::new(Kind::Poly, 3, [16, 12, 10]);
    spec.perturb = 0.25;
    spec.upper = [2.0, 1.5, 1.0];
    let mut global = spec.generate().expect("mesh");
    let (s, c) = 0.7f64.sin_cos();
    for p in &mut global.points {
        let [x, y, z] = *p;
        *p = [c * x - s * y, s * x + c * y, z];
    }
    let geo = Geometry::compute(&global);
    let n = global.n_cells;
    let t: Vec<f64> = geo
        .cell_centroid
        .iter()
        .map(|x| 600.0 + 300.0 * (x[0] + 0.5 * x[1]).sin().abs())
        .collect();
    let kappa: Vec<f64> = geo.cell_centroid.iter().map(|x| 0.2 + 1.5 * x[2]).collect();
    let walls: Vec<Wall> = global.boundaries.iter().map(|b| wall(&b.name)).collect();
    let bundles = std::env::var("RADIATION_BUNDLES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(200_000);
    let settings = Settings {
        bundles,
        ..Settings::default()
    };

    // the reference: the whole mesh on every rank
    let start = Instant::now();
    let one = Domain::whole(&global, &geo)
        .solve(&t, &kappa, &walls, &settings)
        .expect("radiation call");
    let one_seconds = start.elapsed().as_secs_f64();
    let first = global.n_interior_faces();
    let slot: HashMap<u64, usize> = (first..global.n_faces())
        .map(|f| (global.face_id[f], f - first))
        .collect();

    // the same call on the ranks
    let method = if cfg!(feature = "metis") {
        Method::Metis
    } else {
        Method::Rcb
    };
    let part = partition(&global, &geo, size, method).expect("partition");
    let local = extract(&global, &part, rank, 1).expect("local mesh");
    let lgeo = local.geometry(&comm);
    let own = local.mesh.n_cells;
    let pos = &local.cell_position[..own];
    let lt: Vec<f64> = pos.iter().map(|&p| t[p as usize]).collect();
    let lk: Vec<f64> = pos.iter().map(|&p| kappa[p as usize]).collect();
    let lw: Vec<Wall> = local
        .mesh
        .boundaries
        .iter()
        .map(|b| wall(&b.name))
        .collect();
    // the time of the call on the ranks, from a common start to the last
    // rank's end
    comm.all_sum_i64(&mut [0]);
    let start = Instant::now();
    let mine = Domain::rank(&local, &lgeo, &comm)
        .solve(&lt, &lk, &lw, &settings)
        .expect("radiation call");
    comm.all_sum_i64(&mut [0]);
    let rank_seconds = start.elapsed().as_secs_f64();

    let mut cells = vec![0i64; 2 * n];
    for (c, &p) in pos.iter().enumerate() {
        cells[p as usize] = mine.absorbed[c];
        cells[n + p as usize] = mine.emitted[c];
    }
    comm.all_sum_i64(&mut cells);
    let nb = one.wall_absorbed.len();
    let mut faces = vec![0i64; 2 * nb];
    let lfirst = local.mesh.n_interior_faces();
    for f in 0..mine.wall_absorbed.len() {
        let g = slot[&local.mesh.face_id[lfirst + f]];
        faces[g] = mine.wall_absorbed[f];
        faces[nb + g] = mine.wall_emitted[f];
    }
    comm.all_sum_i64(&mut faces);

    let cells_differ = (0..n)
        .filter(|&c| cells[c] != one.absorbed[c] || cells[n + c] != one.emitted[c])
        .count();
    let faces_differ = (0..nb)
        .filter(|&f| faces[f] != one.wall_absorbed[f] || faces[nb + f] != one.wall_emitted[f])
        .count();
    let ledger_same = mine.ledger == one.ledger && mine.bundles == one.bundles;
    let same = cells_differ == 0 && faces_differ == 0 && ledger_same && one.ledger.closes();
    if rank == 0 {
        let l = &one.ledger;
        println!(
            "{size} ranks, {n} cells, {} bundles, {} hand-overs",
            mine.bundles, mine.hand_overs
        );
        println!(
            "  time: {rank_seconds:.2} s on {size} ranks with {} threads each, {one_seconds:.2} s on one rank",
            rayon::current_num_threads()
        );
        println!(
            "  ledger (quanta): emitted {} = gas {} + walls {} + escaped {} + roulette {} + lost {}: {}",
            l.emitted,
            l.absorbed_gas,
            l.absorbed_walls,
            l.escaped,
            l.roulette,
            l.lost,
            if l.closes() {
                "closes"
            } else {
                "DOES NOT CLOSE"
            }
        );
        println!(
            "  ledger on the ranks {}",
            if ledger_same {
                "= 1 rank"
            } else {
                "DIFFERS from 1 rank"
            }
        );
        println!("  absorbed and emitted per cell: {cells_differ} cells differ from 1 rank");
        println!("  absorbed and emitted per wall face: {faces_differ} faces differ from 1 rank");
        println!(
            "radiation on {size} ranks: {}",
            if same {
                "BITWISE EQUAL TO ONE RANK"
            } else {
                "NOT bitwise equal"
            }
        );
    }
    if !same {
        std::process::exit(1);
    }
}
