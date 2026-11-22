//! Feasibility spike for multi-rank MPI (sizing, not a feature): do `n` ranks
//! build the same bits as one rank?
//!
//! ```text
//! cargo build --release -p foehncraft-mesh --no-default-features \
//!     --features mpi,metis --example mpi_matrix
//! mpirun -n 2 target/release/examples/mpi_matrix
//! ```
//!
//! Every rank builds the global mesh, partitions it, extracts its part (two
//! halo layers) and computes the geometry of the local mesh with the halo
//! centroids from the exchange. It then builds the rows of a pure diffusion
//! matrix (face coefficient |S|²/(d·S), the stage-2 form without the
//! non-orthogonal correction) for its own cells and compares them bit for bit
//! with the rows of the one-rank run, keyed by global ids. Four levels:
//!
//! 1. geometry (face area vectors, face centroids, cell volumes and
//!    centroids),
//! 2. off-diagonal entries,
//! 3. the diagonal summed in the local face order (what the stage-2 assembly
//!    does today),
//! 4. the diagonal summed in global face order (D93 (1)).

use foehncraft_mesh::Geometry;
use foehncraft_mesh::Mesh;
use foehncraft_mesh::decompose::LocalMesh;
use foehncraft_mesh::decompose::extract;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::geometry::{add, cross, dot, node_average, polygon_area_centroid, scale, sub};
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_parallel::{Communicator, MpiComm};
use std::collections::HashMap;

/// One row: diagonal (local order), diagonal (global order), off-diagonals
/// keyed by the global id of the neighbour cell.
struct Row {
    diag_local: f64,
    diag_global: f64,
    off: HashMap<u64, f64>,
}

/// Rows of the diffusion matrix for the own cells of `mesh`.
fn rows(mesh: &Mesh, geo: &Geometry) -> HashMap<u64, Row> {
    let ni = mesh.neighbour.len();
    let mut out = HashMap::new();
    for c in 0..mesh.n_cells {
        let mut terms: Vec<(u64, f64)> = Vec::new(); // (global face id, coefficient)
        let mut off = HashMap::new();
        for (f, _) in mesh.cell_face_iter(c) {
            let p = mesh.owner[f] as usize;
            let s = geo.face_area[f];
            let (d, other) = if f < ni {
                let q = mesh.neighbour[f] as usize;
                (
                    sub(geo.cell_centroid[q], geo.cell_centroid[p]),
                    Some(if p == c { q } else { p }),
                )
            } else {
                (sub(geo.face_centroid[f], geo.cell_centroid[p]), None)
            };
            let a = dot(s, s) / dot(d, s);
            terms.push((mesh.face_id[f], a));
            if let Some(o) = other {
                off.insert(mesh.cell_id[o], -a);
            }
        }
        let diag_local = terms.iter().fold(0.0, |s, t| s + t.1);
        terms.sort_by_key(|t| t.0);
        let diag_global = terms.iter().fold(0.0, |s, t| s + t.1);
        out.insert(
            mesh.cell_id[c],
            Row {
                diag_local,
                diag_global,
                off,
            },
        );
    }
    out
}

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

/// Geometry of the local mesh in the canonical way (the proposed fix): every
/// face is computed with the node order and orientation of the global mesh,
/// every cell sums its faces in global face order, halo cells come from the
/// owner by exchange (D93 (1), (3)).
fn canonical_geometry(local: &LocalMesh, comm: &dyn Communicator) -> Geometry {
    let m = &local.mesh;
    let ni = m.neighbour.len();
    // global orientation: the local owner is the global owner unless flipped;
    // for a flipped face extract() reversed the nodes after the first one
    let flipped = |f: usize| -> bool {
        // global owner = lower global cell id in the generated meshes, whose
        // global order is the input order (cell_id)
        f < ni && m.cell_id[m.owner[f] as usize] > m.cell_id[m.neighbour[f] as usize]
    };
    let gnodes = |f: usize| -> Vec<u32> {
        let n = m.face(f);
        if flipped(f) {
            let mut v = vec![n[0]];
            v.extend(n[1..].iter().rev());
            v
        } else {
            n.to_vec()
        }
    };
    let mut face_area = Vec::new();
    let mut face_centroid = Vec::new();
    for f in 0..m.n_faces() {
        let (s, c) = polygon_area_centroid(&m.points, &gnodes(f));
        face_area.push(if flipped(f) { [-s[0], -s[1], -s[2]] } else { s });
        face_centroid.push(c);
    }
    let n = m.n_cells + m.n_halo;
    let mut cell_volume = vec![0.0; n];
    let mut cell_centroid = vec![[0.0; 3]; n];
    for c in 0..m.n_cells {
        // faces of the cell in global face order, with the global sense
        let mut fs: Vec<(u64, usize, bool)> = m
            .cell_face_iter(c)
            .map(|(f, _)| {
                // is_nb in the global sense: the cell is the global neighbour
                let is_nb = f < ni && {
                    let (a, b) = (
                        m.cell_id[m.owner[f] as usize],
                        m.cell_id[m.neighbour[f] as usize],
                    );
                    m.cell_id[c] == a.max(b)
                };
                (m.face_id[f], f, is_nb)
            })
            .collect();
        fs.sort_unstable_by_key(|x| x.0);
        let mut apex = [0.0; 3];
        let mut k = 0.0;
        for &(_, f, _) in &fs {
            apex = add(apex, face_centroid[f]);
            k += 1.0;
        }
        let apex = scale(apex, 1.0 / k);
        let (mut vol, mut cen) = (0.0, [0.0; 3]);
        for &(_, f, is_nb) in &fs {
            let nodes = gnodes(f);
            let kk = nodes.len();
            let fc = if kk == 3 {
                face_centroid[f]
            } else {
                node_average(&m.points, &nodes)
            };
            for i in 0..kk {
                let p = m.points[nodes[i] as usize];
                let q = m.points[nodes[(i + 1) % kk] as usize];
                let mut v = dot(sub(fc, apex), cross(sub(p, apex), sub(q, apex))) / 6.0;
                if is_nb {
                    v = -v;
                }
                vol += v;
                let tc = scale(add(add(apex, fc), add(p, q)), 0.25);
                cen = add(cen, scale(tc, v));
            }
        }
        cell_volume[c] = vol;
        cell_centroid[c] = scale(cen, 1.0 / vol);
    }
    let mut x: Vec<Vec<f64>> = (0..3)
        .map(|d| cell_centroid.iter().map(|c| c[d]).collect())
        .collect();
    let [x0, x1, x2] = &mut x[..] else {
        unreachable!()
    };
    local
        .halo
        .exchange_fields(
            comm,
            foehncraft_mesh::decompose::Transport::HostStaged,
            &mut [&mut cell_volume, x0, x1, x2],
        )
        .expect("exchange");
    for (c, cen) in cell_centroid.iter_mut().enumerate() {
        *cen = [x[0][c], x[1][c], x[2][c]];
    }
    Geometry {
        face_area,
        face_centroid,
        cell_volume,
        cell_centroid,
    }
}

/// Counts of own entities that differ bit for bit from the one-rank run:
/// area, face centroid, volume, cell centroid, off-diagonals, diagonal in
/// local face order, diagonal in global face order.
fn check(
    local: &LocalMesh,
    global: &Mesh,
    ggeo: &Geometry,
    lgeo: &Geometry,
    reference: &HashMap<u64, Row>,
) -> [usize; 7] {
    let gface: HashMap<u64, usize> = global
        .face_id
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let gcell: HashMap<u64, usize> = global
        .cell_id
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let mut bad = [0usize; 7];
    for f in 0..local.mesh.n_faces() {
        let g = gface[&local.mesh.face_id[f]];
        let flipped = f < local.mesh.neighbour.len()
            && global.cell_id[global.owner[g] as usize]
                != local.mesh.cell_id[local.mesh.owner[f] as usize];
        let sg = ggeo.face_area[g];
        let sg = if flipped {
            [-sg[0], -sg[1], -sg[2]]
        } else {
            sg
        };
        if (0..3).any(|d| !same(sg[d], lgeo.face_area[f][d])) {
            bad[0] += 1;
        }
        if (0..3).any(|d| !same(ggeo.face_centroid[g][d], lgeo.face_centroid[f][d])) {
            bad[1] += 1;
        }
    }
    for c in 0..local.mesh.n_cells + local.mesh.n_halo {
        let g = gcell[&local.mesh.cell_id[c]];
        if !same(ggeo.cell_volume[g], lgeo.cell_volume[c]) {
            bad[2] += 1;
        }
        if (0..3).any(|d| !same(ggeo.cell_centroid[g][d], lgeo.cell_centroid[c][d])) {
            bad[3] += 1;
        }
    }
    let mine = rows(&local.mesh, lgeo);
    for (id, row) in &mine {
        let r = &reference[id];
        if row.off.len() != r.off.len()
            || row
                .off
                .iter()
                .any(|(k, v)| r.off.get(k).is_none_or(|w| !same(*v, *w)))
        {
            bad[4] += 1;
        }
        if !same(row.diag_local, r.diag_local) {
            bad[5] += 1;
        }
        if !same(row.diag_global, r.diag_global) {
            bad[6] += 1;
        }
    }
    bad
}

fn main() {
    let comm = MpiComm::init().expect("MPI initialisation");
    let (rank, size) = (comm.rank(), comm.size());
    let kind = match std::env::var("SPIKE_KIND").as_deref() {
        Ok("hex") => Kind::Hex,
        _ => Kind::Poly,
    };
    let mut spec = BoxSpec::new(kind, 3, [12, 10, 8]);
    spec.perturb = 0.2;
    let global = spec.generate().expect("mesh");
    let ggeo = Geometry::compute(&global);
    let reference = rows(&global, &ggeo);
    let method = if cfg!(feature = "metis") {
        Method::Metis
    } else {
        Method::Rcb
    };
    let part = partition(&global, &ggeo, size, method).expect("partition");
    let local = extract(&global, &part, rank, 2).expect("local mesh");
    let today = check(&local, &global, &ggeo, &local.geometry(&comm), &reference);
    let canon = check(
        &local,
        &global,
        &ggeo,
        &canonical_geometry(&local, &comm),
        &reference,
    );
    let names = "area, face centroid, volume, cell centroid, off-diagonals, diagonal (local order), diagonal (global order)";
    println!(
        "rank {rank}/{size}: {} own cells, {} halo | differing [{names}] today {today:?} | canonical {canon:?}",
        local.mesh.n_cells, local.mesh.n_halo,
    );
    let total = (canon.iter().sum::<usize>() - canon[5]) as f64; // the local-order diagonal is shown for comparison only
    let sum = comm.sum_reproducible(total);
    if rank == 0 {
        println!(
            "kind {kind:?}, {size} ranks: canonical geometry {}",
            if sum == 0.0 {
                "BITWISE EQUAL TO ONE RANK"
            } else {
                "NOT bitwise equal"
            }
        );
    }
    if sum != 0.0 {
        std::process::exit(1);
    }
}
