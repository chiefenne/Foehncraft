//! Local meshes with halo cells, and the halo exchange.
//!
//! Local order of a rank's [`LocalMesh`] (owner's recommendation R1 of the
//! MFC assessment, refining D23 and D38):
//!
//! - **cells**: inner cells (all face neighbours on this rank), then the own
//!   cells that touch the halo, then the halo cells, grouped by the rank
//!   that owns them, so a halo exchange receives one contiguous block per
//!   neighbouring rank and needs no unpacking; the inner and the outer cells
//!   are each ordered by reverse Cuthill–McKee for locality, the halo cells
//!   of one rank by layer, then by their order on that rank;
//! - **faces**: inner faces (both cells on this rank), then processor faces
//!   (between an own cell and a halo cell, grouped by neighbouring rank),
//!   then the boundary faces with the same boundaries, in the same order, on
//!   every rank. Every face keeps the owner, the neighbour and the node
//!   order of the global mesh, so a processor face is owned by the halo
//!   cell where the global mesh says so.
//!
//! A solver step can therefore launch the inner faces and inner cells, wait
//! for the halo, then do the processor faces and outer cells, with plain
//! index ranges (D38).
//!
//! Every local cell, face and node keeps its global identity, the position
//! in the input mesh (R2), so output and restart do not depend on the
//! decomposition.
//!
//! **Canonical local mesh** (D93 (1), (3); MPI plan step W1). The geometry
//! and every per-cell sum of a rank are bitwise those of the one-rank run on
//! the global mesh, for any number of ranks and any partition:
//!
//! - every face keeps its global orientation: the owner and neighbour of
//!   the global mesh and its node order (MPI step W3). So the geometry of
//!   a face is bitwise the global one, and so is everything a kernel
//!   computes on the face from the values of its two cells: interpolation
//!   weights, upwind choices (also at a zero mass flux), limiters of the
//!   face. W1 stored a face flipped (owner and neighbour swapped) where
//!   the global owner had the higher local index or was a halo cell, which
//!   kept the geometry exact (a sign) but not the asymmetric face formulas
//!   (the owner weight w against 1 − w); [`Mesh::face_flip`] is all false
//!   for a local mesh since W3 and only marks it as a rank's part;
//! - every own cell lists its faces in ascending position in the global
//!   mesh, so volume, centroid and every later per-cell sum over faces add
//!   the same terms in the same order as on one rank;
//! - halo cells get their volume and centroid from the owning rank by the
//!   exchange ([`LocalMesh::geometry`]), never recomputed.
//!
//! The reference is the global mesh handed to [`extract`] (after any
//! renumbering), so one rank through [`extract`] gives the bits of the
//! undecomposed run, and runs without a decomposition do not change.
//!
//! The halo exchange (R3) gathers the send values of all fields for one
//! neighbouring rank into one buffer, sends one message per neighbour and
//! copies the received block into the contiguous halo range of each field.
//! The [`Transport`] is chosen at run time: host-staged now; GPU-aware MPI,
//! with send lists and buffers resident on the device, comes with the
//! device fields of stage 2 and is tested locally (D22). The halo depth is a
//! parameter, one layer by default: enough for the face-neighbour
//! least-squares gradients (D24) and for particles, which are handed over at
//! processor faces (D47).
//!
//! In stage 1 every rank holds the whole mesh (generated or read on every
//! rank) and extracts its part with [`extract`]; distributing a mesh read on
//! one rank follows with the parallel I/O of later stages.

use crate::geometry::Geometry;
use crate::mesh::{Boundary, Mesh};
use crate::renumber::rcm_order;
use foehncraft_kernels::unpack_face;
use foehncraft_parallel::Communicator;
use serde::{Deserialize, Serialize};

/// What a rank exchanges with one neighbouring rank.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Neighbour {
    /// The neighbouring rank.
    pub rank: usize,
    /// Own cells (local indices) the neighbour holds as halo cells, in the
    /// order of its halo: the gather list of the pack step.
    pub send: Vec<u32>,
    /// Halo cells (local indices) owned by the neighbour: a contiguous
    /// range, filled without unpacking.
    pub recv: std::ops::Range<usize>,
}

/// The halo exchange plan of one rank.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct HaloPlan {
    /// Neighbouring ranks, ascending.
    pub neighbours: Vec<Neighbour>,
}

/// How halo buffers travel between ranks; chosen at run time (R3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "kebab-case")]
pub enum Transport {
    /// Pack on the host (or copy device buffers to the host), MPI, copy
    /// back. Always available; tested in CI.
    #[default]
    HostStaged,
    /// Device buffers handed to a GPU-aware MPI library. Needs device fields
    /// (stage 2 onward); tested locally (D22).
    GpuAware,
}

/// The part of a mesh held by one rank.
#[derive(Debug, Clone)]
pub struct LocalMesh {
    /// The rank.
    pub rank: usize,
    /// The local mesh: inner cells, outer cells, halo cells; its
    /// `cell_id`, `face_id` and `point_id` are the global identities.
    pub mesh: Mesh,
    /// Halo layer of each halo cell (1 = face neighbour of an own cell).
    pub halo_layer: Vec<u8>,
    /// Number of processor faces per neighbour, in the order of
    /// `halo.neighbours`.
    pub processor_faces: Vec<usize>,
    /// The halo exchange plan.
    pub halo: HaloPlan,
}

/// Own cells of `rank` in local order (inner cells by RCM, then outer cells
/// by RCM), the number of inner cells, and the halo cells (owner rank,
/// layer, cell), all as positions in `global`.
#[allow(clippy::type_complexity)]
fn local_cells(
    global: &Mesh,
    part: &[u32],
    rank: u32,
    layers: usize,
) -> (Vec<u32>, usize, Vec<(u32, u8, u32)>) {
    let n = global.n_cells;
    let mut layer = vec![0u8; n]; // 0: not reached, MAX: own
    let own: Vec<u32> = (0..n as u32)
        .filter(|&c| part[c as usize] == rank)
        .collect();
    for &c in &own {
        layer[c as usize] = u8::MAX;
    }
    let mut halo = Vec::new();
    let mut front = own.clone();
    for l in 1..=layers {
        let mut next = Vec::new();
        for &c in &front {
            for (f, _) in global.cell_face_iter(c as usize) {
                if let Some(d) = global.across(f, c as usize)
                    && layer[d] == 0
                {
                    layer[d] = l as u8;
                    next.push(d as u32);
                    halo.push((part[d], l as u8, d as u32));
                }
            }
        }
        front = next;
    }
    // inner (no face neighbour in the halo) and outer own cells, each in RCM
    // order on the subgraph of own cells
    let mut local = vec![u32::MAX; n];
    let (mut inner, mut outer) = (Vec::new(), Vec::new());
    for &c in &own {
        let touches = global
            .cell_face_iter(c as usize)
            .any(|(f, _)| global.across(f, c as usize).is_some_and(|d| layer[d] == 1));
        if touches {
            outer.push(c);
        } else {
            inner.push(c);
        }
    }
    let mut ordered = Vec::with_capacity(own.len());
    for group in [&inner, &outer] {
        for (i, &c) in group.iter().enumerate() {
            local[c as usize] = i as u32;
        }
        let order = rcm_order(group.len(), |i, out| {
            let c = group[i] as usize;
            for (f, _) in global.cell_face_iter(c) {
                if let Some(d) = global.across(f, c)
                    && local[d] != u32::MAX
                    && group.get(local[d] as usize) == Some(&(d as u32))
                {
                    out.push(local[d] as usize);
                }
            }
        });
        ordered.extend(order.into_iter().map(|i| group[i]));
        for &c in group.iter() {
            local[c as usize] = u32::MAX;
        }
    }
    halo.sort_unstable();
    (ordered, inner.len(), halo)
}

/// The local mesh of `rank` from the global mesh and the partition.
pub fn extract(
    global: &Mesh,
    part: &[u32],
    rank: usize,
    layers: usize,
) -> Result<LocalMesh, String> {
    if global.n_halo != 0 {
        return Err("the global mesh must not have halo cells".into());
    }
    if layers == 0 || layers > 8 {
        return Err("halo layers must lie in 1..=8".into());
    }
    let r = rank as u32;
    let n = global.n_cells;
    let (own, n_inner, halo) = local_cells(global, part, r, layers);
    if own.is_empty() {
        return Err(format!("rank {rank} has no cells"));
    }
    let mut local = vec![u32::MAX; n];
    for (i, &c) in own.iter().enumerate() {
        local[c as usize] = i as u32;
    }
    let n_own = own.len();
    for (i, &(_, _, c)) in halo.iter().enumerate() {
        local[c as usize] = (n_own + i) as u32;
    }

    // faces of own cells: (class, key, local owner, local neighbour, global
    // face, reverse); class 0 inner, 1 processor (key = neighbour rank),
    // 2 + boundary for boundary faces
    let ni = global.n_interior_faces();
    let mut faces: Vec<(u32, u32, u32, u32, u32, bool)> = Vec::new();
    for f in 0..global.n_faces() {
        let p = global.owner[f] as usize;
        if f < ni {
            let q = global.neighbour[f] as usize;
            let (lp, lq) = (local[p], local[q]);
            // every face keeps the owner and neighbour of the global mesh
            // (MPI step W3), so every face formula sees its inputs in the
            // roles of the one-rank run; a processor face may be owned by
            // the halo cell
            match (part[p] == r, part[q] == r) {
                (true, true) => faces.push((0, 0, lp, lq, f as u32, false)),
                (true, false) => faces.push((1, part[q], lp, lq, f as u32, false)),
                (false, true) => faces.push((1, part[p], lp, lq, f as u32, false)),
                (false, false) => {}
            }
        } else if part[p] == r {
            let boundary = global.boundary_of(f).unwrap() as u32;
            faces.push((2 + boundary, 0, local[p], 0, f as u32, false));
        }
    }
    faces.sort_unstable();

    // nodes, in global order
    let mut node_local = vec![u32::MAX; global.n_points()];
    let mut used: Vec<u32> = Vec::new();
    for x in &faces {
        for &v in global.face(x.4 as usize) {
            if node_local[v as usize] == u32::MAX {
                node_local[v as usize] = 0;
                used.push(v);
            }
        }
    }
    used.sort_unstable();
    for (i, &v) in used.iter().enumerate() {
        node_local[v as usize] = i as u32;
    }

    let mut cell_id: Vec<u64> = own.iter().map(|&c| global.cell_id[c as usize]).collect();
    cell_id.extend(halo.iter().map(|h| global.cell_id[h.2 as usize]));
    let mut mesh = Mesh {
        points: used.iter().map(|&v| global.points[v as usize]).collect(),
        point_id: used.iter().map(|&v| global.point_id[v as usize]).collect(),
        n_cells: n_own,
        n_inner_cells: n_inner,
        n_halo: halo.len(),
        cell_zone: own.iter().map(|&c| global.cell_zone[c as usize]).collect(),
        zones: global.zones.clone(),
        cell_id,
        ..Default::default()
    };
    mesh.face_ptr.push(0);
    mesh.face_flip = faces.iter().map(|x| x.5).collect();
    let mut processor_faces: Vec<(u32, usize)> = Vec::new();
    for &(class, key, a, b, f, rev) in &faces {
        let nodes = global.face(f as usize);
        if rev {
            mesh.face_nodes.push(node_local[nodes[0] as usize]);
            mesh.face_nodes
                .extend(nodes[1..].iter().rev().map(|&v| node_local[v as usize]));
        } else {
            mesh.face_nodes
                .extend(nodes.iter().map(|&v| node_local[v as usize]));
        }
        mesh.face_ptr.push(mesh.face_nodes.len() as u32);
        mesh.face_id.push(global.face_id[f as usize]);
        mesh.owner.push(a);
        match class {
            0 => mesh.neighbour.push(b),
            1 => {
                mesh.neighbour.push(b);
                match processor_faces.last_mut() {
                    Some(x) if x.0 == key => x.1 += 1,
                    _ => processor_faces.push((key, 1)),
                }
            }
            _ => {}
        }
    }
    mesh.n_inner_faces = faces.iter().filter(|x| x.0 == 0).count();
    let mut start = mesh.neighbour.len();
    for (i, p) in global.boundaries.iter().enumerate() {
        let len = faces.iter().filter(|x| x.0 == 2 + i as u32).count();
        mesh.boundaries.push(Boundary {
            name: p.name.clone(),
            start,
            len,
        });
        start += len;
    }
    mesh.build_cell_faces();
    // every cell lists its faces in ascending position in the global mesh,
    // so per-cell sums run in global face order (D93 (1))
    for c in 0..n_own {
        let range = mesh.cell_ptr[c] as usize..mesh.cell_ptr[c + 1] as usize;
        mesh.cell_faces[range].sort_unstable_by_key(|&e| faces[unpack_face(e).0].4);
    }

    // exchange plan: for each neighbouring rank, the own cells in its halo,
    // in its halo order
    let mut neighbours: Vec<Neighbour> = Vec::new();
    let mut i = 0;
    while i < halo.len() {
        let s = halo[i].0;
        let j = i + halo[i..].iter().take_while(|h| h.0 == s).count();
        let (_, _, their_halo) = local_cells(global, part, s, layers);
        let send = their_halo
            .into_iter()
            .filter(|h| h.0 == r)
            .map(|h| local[h.2 as usize])
            .collect();
        neighbours.push(Neighbour {
            rank: s as usize,
            send,
            recv: n_own + i..n_own + j,
        });
        i = j;
    }
    let processor_faces = neighbours
        .iter()
        .map(|nb| {
            processor_faces
                .iter()
                .find(|x| x.0 as usize == nb.rank)
                .map_or(0, |x| x.1)
        })
        .collect();

    Ok(LocalMesh {
        rank,
        mesh,
        halo_layer: halo.iter().map(|h| h.1).collect(),
        processor_faces,
        halo: HaloPlan { neighbours },
    })
}

/// The local meshes of all ranks.
pub fn decompose(
    global: &Mesh,
    part: &[u32],
    parts: usize,
    layers: usize,
) -> Result<Vec<LocalMesh>, String> {
    if global.periodic.is_some() {
        return Err(
            "periodic meshes run on one partition (the decomposition follows with \
                    the multi-rank stage)"
                .into(),
        );
    }
    (0..parts)
        .map(|r| extract(global, part, r, layers))
        .collect()
}

/// The communicator shared by the solver structures of a rank: `Send +
/// Sync`, so a structure that runs its loops on worker threads can hold it
/// (only the rank's own thread calls it).
pub type SharedComm = std::sync::Arc<dyn Communicator + Send + Sync>;

/// What a rank's solver needs to work with the other ranks (MPI step W3):
/// the communicator and the halo exchange plan of its local mesh.
#[derive(Clone)]
pub struct Ranks {
    /// The communicator.
    pub comm: SharedComm,
    /// The exchange plan of the rank's local mesh.
    pub halo: HaloPlan,
}

impl std::fmt::Debug for Ranks {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Ranks")
            .field("rank", &self.comm.rank())
            .field("size", &self.comm.size())
            .field("neighbours", &self.halo.neighbours.len())
            .finish()
    }
}

impl Ranks {
    /// The link of `local`'s rank through `comm`.
    pub fn new(comm: SharedComm, local: &LocalMesh) -> Self {
        Self {
            comm,
            halo: local.halo.clone(),
        }
    }

    /// Fills the halo entries of `fields` (one value per own and halo
    /// cell) from their owners, one message per neighbouring rank.
    pub fn exchange(&self, fields: &mut [&mut [f64]]) {
        self.halo
            .exchange_fields(&*self.comm, Transport::HostStaged, fields)
            .expect("host-staged exchange");
    }
}

impl HaloPlan {
    /// Fills the halo entries of all `fields` (each one value per own and
    /// halo cell) from the neighbouring ranks: one message per neighbour
    /// carries all fields.
    pub fn exchange_fields(
        &self,
        comm: &dyn Communicator,
        transport: Transport,
        fields: &mut [&mut [f64]],
    ) -> Result<(), String> {
        if transport == Transport::GpuAware {
            return Err(
                "the GPU-aware halo transport needs device fields (stage 2 onward); \
                 use host-staged"
                    .into(),
            );
        }
        let mut send = Vec::new();
        let mut recv = Vec::new();
        for nb in &self.neighbours {
            // pack: one contiguous buffer, field by field
            send.clear();
            for f in fields.iter() {
                send.extend(nb.send.iter().map(|&c| f[c as usize]));
            }
            let n = nb.recv.len();
            recv.clear();
            recv.resize(n * fields.len(), 0.0);
            comm.send_recv_f64(&send, nb.rank, &mut recv, nb.rank);
            // contiguous halo range per field: a copy, no unpack
            for (k, f) in fields.iter_mut().enumerate() {
                f[nb.recv.clone()].copy_from_slice(&recv[k * n..(k + 1) * n]);
            }
        }
        Ok(())
    }

    /// One field, host-staged.
    pub fn exchange(&self, comm: &dyn Communicator, field: &mut [f64]) {
        self.exchange_fields(comm, Transport::HostStaged, &mut [field])
            .expect("host-staged exchange");
    }

    /// The exchange for all ranks of a decomposition held in one process
    /// (tests and tools): `fields[r]` is the field of rank `r`.
    pub fn exchange_all(plans: &[&HaloPlan], fields: &mut [Vec<f64>]) {
        for r in 0..plans.len() {
            for nb in &plans[r].neighbours {
                let other = plans[nb.rank]
                    .neighbours
                    .iter()
                    .find(|x| x.rank == r)
                    .expect("halo plans are symmetric");
                let vals: Vec<f64> = other
                    .send
                    .iter()
                    .map(|&c| fields[nb.rank][c as usize])
                    .collect();
                fields[r][nb.recv.clone()].copy_from_slice(&vals);
            }
        }
    }
}

impl LocalMesh {
    /// Geometry of the local mesh, with the halo volumes and centroids
    /// filled by one exchange of four fields: bitwise the geometry of the
    /// global mesh for every own face and own and halo cell (see the module
    /// docs).
    pub fn geometry(&self, comm: &dyn Communicator) -> Geometry {
        let mut g = Geometry::compute(&self.mesh);
        let mut x: Vec<Vec<f64>> = (0..3)
            .map(|d| g.cell_centroid.iter().map(|c| c[d]).collect())
            .collect();
        let [x0, x1, x2] = &mut x[..] else {
            unreachable!()
        };
        self.halo
            .exchange_fields(
                comm,
                Transport::HostStaged,
                &mut [&mut g.cell_volume, x0, x1, x2],
            )
            .expect("host-staged exchange");
        for (c, cen) in g.cell_centroid.iter_mut().enumerate() {
            *cen = [x[0][c], x[1][c], x[2][c]];
        }
        g
    }
}

/// What [`canonical_check`] found on one rank.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct CanonicalRank {
    /// The rank.
    pub rank: usize,
    /// Faces stored flipped against the global mesh.
    pub flipped_faces: usize,
    /// Faces whose area vector or centroid differ in any bit from the
    /// global mesh (0 when the local mesh is canonical).
    pub differing_faces: usize,
    /// Own and halo cells whose volume or centroid differ in any bit
    /// (0 when the local mesh is canonical).
    pub differing_cells: usize,
}

/// Checks that the local meshes of a decomposition held in one process are
/// canonical (MPI plan step W1, D93 (1), (3)): every own face and every own
/// and halo cell has bitwise the geometry of `global` (halo filled by
/// [`HaloPlan::exchange_all`]). `global_geo` is
/// `Geometry::compute(global)`. Entities are matched by their global
/// identities.
pub fn canonical_check(
    global: &Mesh,
    global_geo: &Geometry,
    locals: &[LocalMesh],
) -> Vec<CanonicalRank> {
    let same = |a: [f64; 3], b: [f64; 3]| (0..3).all(|d| a[d].to_bits() == b[d].to_bits());
    let gface: std::collections::HashMap<u64, usize> = global
        .face_id
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let gcell: std::collections::HashMap<u64, usize> = global
        .cell_id
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let mut geos: Vec<Geometry> = locals.iter().map(|l| Geometry::compute(&l.mesh)).collect();
    let plans: Vec<&HaloPlan> = locals.iter().map(|l| &l.halo).collect();
    let mut vol: Vec<Vec<f64>> = geos.iter().map(|g| g.cell_volume.clone()).collect();
    HaloPlan::exchange_all(&plans, &mut vol);
    for d in 0..3 {
        let mut x: Vec<Vec<f64>> = geos
            .iter()
            .map(|g| g.cell_centroid.iter().map(|c| c[d]).collect())
            .collect();
        HaloPlan::exchange_all(&plans, &mut x);
        for (g, x) in geos.iter_mut().zip(&x) {
            for (c, v) in g.cell_centroid.iter_mut().zip(x) {
                c[d] = *v;
            }
        }
    }
    locals
        .iter()
        .zip(geos.iter_mut().zip(vol))
        .map(|(l, (g, v))| {
            g.cell_volume = v;
            let m = &l.mesh;
            let mut out = CanonicalRank {
                rank: l.rank,
                ..Default::default()
            };
            for f in 0..m.n_faces() {
                let gf = gface[&m.face_id[f]];
                let s = global_geo.face_area[gf];
                let s = if m.is_flipped(f) {
                    out.flipped_faces += 1;
                    [-s[0], -s[1], -s[2]]
                } else {
                    s
                };
                if !same(s, g.face_area[f])
                    || !same(global_geo.face_centroid[gf], g.face_centroid[f])
                {
                    out.differing_faces += 1;
                }
            }
            for c in 0..m.n_cells + m.n_halo {
                let gc = gcell[&m.cell_id[c]];
                if global_geo.cell_volume[gc].to_bits() != g.cell_volume[c].to_bits()
                    || !same(global_geo.cell_centroid[gc], g.cell_centroid[c])
                {
                    out.differing_cells += 1;
                }
            }
            out
        })
        .collect()
}
/// One reduction of [`reduction_check`]: the exact value on the ranks and
/// on one rank, and the old rank-order pairwise value for comparison.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct ReductionValue {
    /// What is summed.
    pub name: String,
    /// The exact sum over the ranks (local accumulators merged as the
    /// integer all-reduce of a multi-rank run does it).
    pub ranks: f64,
    /// The exact sum over the whole mesh on one rank.
    pub one_rank: f64,
    /// Whether the two are bitwise equal (always; the run fails otherwise).
    pub bitwise_equal: bool,
    /// The sum before MPI step W2: a pairwise sum per rank, the rank sums
    /// combined pairwise in rank order.
    pub pairwise_ranks: f64,
    /// Distance of `pairwise_ranks` from the exact one-rank value in units
    /// in the last place (0: the same bits).
    pub pairwise_ulps: u64,
}

/// Distance of two finite values of the same sign in units in the last
/// place.
fn ulps(a: f64, b: f64) -> u64 {
    let key = |x: f64| {
        let b = x.to_bits() as i64;
        if b < 0 { i64::MIN - b } else { b }
    };
    key(a).abs_diff(key(b))
}

/// Exact reductions over a decomposition held in one process (MPI plan step
/// W2, D31 (1), D93 (4)): the domain volume Σ V, its first moment Σ V x and
/// the boundary area Σ |S| over the boundary faces, each summed per rank
/// over the own cells and faces into an exact accumulator, the accumulators
/// merged through their integer images (what `Communicator::sum_exact`
/// sends), and compared bit for bit with the exact sum over `global`.
pub fn reduction_check(
    global: &Mesh,
    global_geo: &Geometry,
    locals: &[LocalMesh],
) -> Vec<ReductionValue> {
    use foehncraft_core::exact::{ExactSum, WORDS};
    use foehncraft_core::repro::pairwise_sum;
    let boundary: std::collections::HashSet<u64> = global.face_id[global.neighbour.len()..]
        .iter()
        .copied()
        .collect();
    // the three per-entity terms of a mesh, own cells `0..n` only
    let terms = |m: &Mesh, g: &Geometry, n: usize| -> [Vec<f64>; 3] {
        let v = g.cell_volume[..n].to_vec();
        let vx = (0..n)
            .map(|c| g.cell_volume[c] * g.cell_centroid[c][0])
            .collect();
        let area = (m.neighbour.len()..m.n_faces())
            .filter(|&f| (m.owner[f] as usize) < n && boundary.contains(&m.face_id[f]))
            .map(|f| crate::geometry::norm(g.face_area[f]))
            .collect();
        [v, vx, area]
    };
    let reference = terms(global, global_geo, global.n_cells);
    let mut words = [[0i64; WORDS]; 3];
    let mut rank_sums = [vec![], vec![], vec![]];
    for l in locals {
        let g = Geometry::compute(&l.mesh);
        let t = terms(&l.mesh, &g, l.mesh.n_cells);
        for q in 0..3 {
            let mut acc = ExactSum::new();
            acc.add_slice(&t[q]);
            for (w, x) in words[q].iter_mut().zip(acc.to_words()) {
                *w += x;
            }
            rank_sums[q].push(pairwise_sum(&t[q]));
        }
    }
    ["volume", "first moment V x", "boundary area"]
        .iter()
        .enumerate()
        .map(|(q, name)| {
            let mut one = ExactSum::new();
            one.add_slice(&reference[q]);
            let one_rank = one.value();
            let ranks = ExactSum::from_words(&words[q]).value();
            let pairwise_ranks = pairwise_sum(&rank_sums[q]);
            ReductionValue {
                name: (*name).to_string(),
                ranks,
                one_rank,
                bitwise_equal: ranks.to_bits() == one_rank.to_bits(),
                pairwise_ranks,
                pairwise_ulps: ulps(pairwise_ranks, one_rank),
            }
        })
        .collect()
}
