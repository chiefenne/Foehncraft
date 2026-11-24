//! Radiation R1a (D257, D258): the cell walk handed over between ranks
//! gives bitwise the walk of one rank.
//!
//! ```text
//! cargo build --release -p foehncraft-mesh --features mpi --example mpi_walk
//! mpirun -n 1 target/release/examples/mpi_walk
//! mpirun -n 4 target/release/examples/mpi_walk
//! ```
//!
//! Every rank builds the global mesh (perturbed and turned, so no wall lies
//! in a coordinate plane), partitions it (METIS, or RCB without the `metis`
//! feature) and extracts its part with one halo layer. Rays start in every
//! third cell of the global mesh, on the rank that owns the cell, in keyed
//! random directions, and walk to the wall. A walk that reaches a processor
//! face is handed to the rank behind it ([`HandOver`]); the ranks repeat
//! until the global count of walks in flight is zero.
//!
//! Each rank tallies, per cell of the global mesh, the path length in
//! 64-bit fixed point (the tally of the radiation, D107 (3)), and a hash of
//! every interval (ray, cell, the bits of λ₀ and λ₁) and of every end (ray,
//! face, the bits of λ and of the exit point). The tallies are summed over
//! the ranks with an integer all-reduce and compared bit for bit with the
//! same walk on the whole mesh, which every rank computes as the reference.
//! The example exits with status 1 when anything differs.

use foehncraft_mesh::decompose::extract;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::geometry::{norm, sub};
use foehncraft_mesh::handover::HandOver;
use foehncraft_mesh::locate::{Locator, Walk, WalkEnd, step_limit};
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_mesh::{Geometry, Mesh};
use foehncraft_parallel::{Communicator, MpiComm};

/// Fixed-point scale of the length tally: 2⁴⁰ per metre.
const SCALE: f64 = (1u64 << 40) as f64;

fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// A 40-bit hash of a record, so the sum of millions of them stays far
/// below the 64-bit limit.
fn hash(words: &[u64]) -> i64 {
    let mut h = 0u64;
    for &w in words {
        h = mix(h ^ w);
    }
    (h >> 24) as i64
}

fn direction(i: u64) -> [f64; 3] {
    let u = foehncraft_core::rng::uniform([21, i, 0]);
    let v = foehncraft_core::rng::uniform([21, i, 1]);
    let z = 2.0 * u - 1.0;
    let r = (1.0 - z * z).sqrt();
    let phi = 2.0 * std::f64::consts::PI * v;
    [r * phi.cos(), r * phi.sin(), z]
}

/// The tallies: length per global cell, then the number and hash sum of
/// the intervals and of the ends.
struct Tally {
    length: Vec<i64>,
    sums: [i64; 4],
}

impl Tally {
    fn new(n: usize) -> Self {
        Self {
            length: vec![0; n],
            sums: [0; 4],
        }
    }

    fn interval(&mut self, ray: u64, cell: u32, a: f64, b: f64, len: f64) {
        self.length[cell as usize] += ((b - a) * len * SCALE).round() as i64;
        self.sums[0] += 1;
        self.sums[1] += hash(&[ray, u64::from(cell), a.to_bits(), b.to_bits()]);
    }

    fn end(&mut self, ray: u64, face: u64, w: &Walk) {
        self.sums[2] += 1;
        self.sums[3] += hash(&[
            ray,
            face,
            w.lambda.to_bits(),
            w.x[0].to_bits(),
            w.x[1].to_bits(),
            w.x[2].to_bits(),
            u64::from(w.steps),
            u64::from(w.crossed),
        ]);
    }
}

fn rotate(m: &mut Mesh) {
    let (s, c) = 0.7f64.sin_cos();
    let (s2, c2) = 0.4f64.sin_cos();
    for p in &mut m.points {
        let [x, y, z] = *p;
        let (x, y) = (c * x - s * y, s * x + c * y);
        let (y, z) = (c2 * y - s2 * z, s2 * y + c2 * z);
        *p = [x, y, z];
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
    let mut spec = BoxSpec::new(kind, 3, [16, 12, 10]);
    spec.perturb = 0.25;
    spec.upper = [2.0, 1.5, 1.0];
    let mut global = spec.generate().expect("mesh");
    rotate(&mut global);
    let ggeo = Geometry::compute(&global);
    let n = global.n_cells;
    let limit = step_limit(n);
    let rays: Vec<(u64, usize)> = (0..n as u64)
        .filter(|i| i % 3 == 0)
        .map(|i| (i, i as usize))
        .collect();
    let target = |x: [f64; 3], ray: u64| -> [f64; 3] {
        let d = direction(ray);
        std::array::from_fn(|k| x[k] + 10.0 * d[k])
    };

    // the reference: every ray on the whole mesh
    let gloc = Locator::new(&global, &ggeo);
    let mut reference = Tally::new(n);
    for &(ray, c) in &rays {
        let from = gloc.start(c);
        let to = target(from.x, ray);
        let len = norm(sub(to, from.x));
        let mut w = gloc.start_walk(from, to);
        let end = gloc.walk(&mut w, limit, |cell, a, b| {
            reference.interval(ray, cell, a, b, len);
            true
        });
        let WalkEnd::Boundary { face } = end else {
            panic!("reference ray {ray}: {end:?}")
        };
        reference.end(ray, global.face_id[face as usize], &w);
    }

    // the same rays on the ranks
    let method = if cfg!(feature = "metis") {
        Method::Metis
    } else {
        Method::Rcb
    };
    let part = partition(&global, &ggeo, size, method).expect("partition");
    let local = extract(&global, &part, rank, 1).expect("local mesh");
    let lgeo = local.geometry(&comm);
    let loc = Locator::new(&local.mesh, &lgeo);
    let hand = HandOver::new(&local, 1);
    let own = local.mesh.n_cells;
    let pos = &local.cell_position;
    let mut local_of = vec![u32::MAX; n];
    for (c, &p) in pos[..own].iter().enumerate() {
        local_of[p as usize] = c as u32;
    }
    let mut queue: Vec<(Walk, u64)> = rays
        .iter()
        .filter(|&&(_, c)| local_of[c] != u32::MAX)
        .map(|&(ray, c)| {
            let from = loc.start(local_of[c] as usize);
            (loc.start_walk(from, target(from.x, ray)), ray)
        })
        .collect();
    let mut tally = Tally::new(n);
    let mut outbox = hand.outbox();
    let mut rounds = 0usize;
    let mut handed = 0i64;
    loop {
        let mut in_flight = [queue.len() as i64];
        comm.all_sum_i64(&mut in_flight);
        if in_flight[0] == 0 {
            break;
        }
        rounds += 1;
        for (mut w, ray) in std::mem::take(&mut queue) {
            let len = norm(sub(w.target, w.origin));
            let end = loc.walk(&mut w, limit, |cell, a, b| {
                tally.interval(ray, pos[cell as usize], a, b, len);
                true
            });
            match end {
                WalkEnd::Boundary { face } => {
                    tally.end(ray, local.mesh.face_id[face as usize], &w);
                }
                WalkEnd::Processor { .. } => {
                    hand.pack(&w, &[ray as f64], &mut outbox);
                    handed += 1;
                }
                e => panic!("rank {rank}, ray {ray}: {e:?}"),
            }
        }
        let inbox = hand.exchange(&comm, &mut outbox);
        for (w, payload) in hand.unpack(&inbox) {
            queue.push((w, payload[0] as u64));
        }
    }
    let mut sums = tally.sums.to_vec();
    sums.push(handed);
    comm.all_sum_i64(&mut sums);
    comm.all_sum_i64(&mut tally.length);

    let names = ["intervals", "interval hash", "ends", "end hash"];
    let mut differ = 0usize;
    if rank == 0 {
        println!(
            "kind {kind:?}, {size} ranks, {n} cells, {} rays, {rounds} rounds, {} hand-overs",
            rays.len(),
            sums[4]
        );
    }
    for (i, name) in names.iter().enumerate() {
        let same = sums[i] == reference.sums[i];
        differ += usize::from(!same);
        if rank == 0 {
            println!(
                "  {name:<14} {:>20}  {}",
                sums[i],
                if same {
                    "= 1 rank"
                } else {
                    "DIFFERS from 1 rank"
                }
            );
        }
    }
    let cells_differ = tally
        .length
        .iter()
        .zip(&reference.length)
        .filter(|(a, b)| a != b)
        .count();
    differ += cells_differ;
    if rank == 0 {
        println!("  length per cell: {cells_differ} cells differ from 1 rank");
        println!(
            "walk on {size} ranks: {}",
            if differ == 0 {
                "BITWISE EQUAL TO ONE RANK"
            } else {
                "NOT bitwise equal"
            }
        );
    }
    if differ != 0 {
        std::process::exit(1);
    }
}
