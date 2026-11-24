//! Hand-over of walks between ranks (D47, D258).
//!
//! A walk ([`crate::locate::Walk`]) that reaches a processor face stops with
//! the state of the halo cell behind the face. [`HandOver::pack`] writes
//! that state with global identities (the cell's position in the global
//! mesh, the face's global identity; the edge index is the same on every
//! rank because every face keeps the node order of the global mesh) and the
//! caller's payload into the buffer for the rank that owns the cell.
//! [`HandOver::exchange`] sends one message per neighbouring rank, the
//! lengths first, and [`HandOver::unpack`] gives the walk back in the local
//! indices of the receiving rank, which goes on with the same numbers. On a
//! canonical local mesh (MPI step W1) the path and every interval of a walk
//! are therefore bitwise those of the walk on one rank.
//!
//! A driver walks all its walks, exchanges, and repeats until the global
//! count of walks in flight (an integer all-reduce) is zero.

use crate::decompose::LocalMesh;
use crate::locate::Walk;
use foehncraft_parallel::Communicator;
use std::collections::HashMap;

/// Numbers of a walk in a hand-over record; the payload follows.
pub const WALK_LEN: usize = 16;

/// Largest integer a record carries exactly (2⁵³).
const EXACT: u64 = 1 << 53;

/// The hand-over of walks for one rank's local mesh.
#[derive(Debug)]
pub struct HandOver<'a> {
    local: &'a LocalMesh,
    /// Own cell of each global cell position held by this rank.
    own: HashMap<u32, u32>,
    /// Slot in `local.halo.neighbours` of each halo cell.
    slot: Vec<u32>,
    /// Numbers per record: [`WALK_LEN`] plus the payload.
    record: usize,
}

impl<'a> HandOver<'a> {
    /// The hand-over of `local`'s rank, for records with `payload` numbers
    /// after the walk.
    pub fn new(local: &'a LocalMesh, payload: usize) -> Self {
        let n = local.mesh.n_cells;
        let own = local.cell_position[..n]
            .iter()
            .enumerate()
            .map(|(c, &p)| (p, c as u32))
            .collect();
        let mut slot = vec![u32::MAX; local.mesh.n_halo];
        for (s, nb) in local.halo.neighbours.iter().enumerate() {
            for h in nb.recv.clone() {
                slot[h - n] = s as u32;
            }
        }
        Self {
            local,
            own,
            slot,
            record: WALK_LEN + payload,
        }
    }

    /// Numbers per record.
    pub fn record_len(&self) -> usize {
        self.record
    }

    /// Empty buffers, one per neighbouring rank.
    pub fn outbox(&self) -> Vec<Vec<f64>> {
        vec![Vec::new(); self.local.halo.neighbours.len()]
    }

    /// The rank that owns the halo cell `w` stands in.
    pub fn destination(&self, w: &Walk) -> usize {
        let s = self.slot_of(w);
        self.local.halo.neighbours[s].rank
    }

    fn slot_of(&self, w: &Walk) -> usize {
        let n = self.local.mesh.n_cells;
        let c = w.cell as usize;
        assert!(
            c >= n && c < n + self.local.mesh.n_halo,
            "hand-over: walk in cell {c}, not in a halo cell"
        );
        self.slot[c - n] as usize
    }

    /// Writes walk `w` (stopped at a processor face, in the halo cell behind
    /// it) and its `payload` into the buffer of the owning rank.
    pub fn pack(&self, w: &Walk, payload: &[f64], outbox: &mut [Vec<f64>]) {
        assert_eq!(WALK_LEN + payload.len(), self.record, "payload length");
        let s = self.slot_of(w);
        let mesh = &self.local.mesh;
        let id = |v: u64| {
            assert!(v < EXACT, "hand-over: identity {v} above 2^53");
            v as f64
        };
        let out = &mut outbox[s];
        out.extend_from_slice(&w.origin);
        out.extend_from_slice(&w.target);
        out.push(w.lambda);
        out.extend_from_slice(&w.x);
        out.push(id(u64::from(self.local.cell_position[w.cell as usize])));
        out.push(id(mesh.face_id[w.face as usize]));
        out.push(f64::from(w.edge));
        out.push(f64::from(w.entered));
        out.push(f64::from(w.steps));
        out.push(f64::from(w.crossed));
        out.extend_from_slice(payload);
    }

    /// Sends every buffer of `outbox` to its rank and returns the records
    /// received, from the neighbouring ranks in ascending order. The
    /// outbox is left empty.
    pub fn exchange(&self, comm: &dyn Communicator, outbox: &mut [Vec<f64>]) -> Vec<f64> {
        let mut inbox = Vec::new();
        for (nb, out) in self.local.halo.neighbours.iter().zip(outbox.iter_mut()) {
            let mut len = [0.0];
            comm.send_recv_f64(&[out.len() as f64], nb.rank, &mut len, nb.rank);
            let mut recv = vec![0.0; len[0] as usize];
            comm.send_recv_f64(out, nb.rank, &mut recv, nb.rank);
            inbox.extend_from_slice(&recv);
            out.clear();
        }
        inbox
    }

    /// The records of `inbox`, each the walk in local indices and its
    /// payload.
    pub fn unpack<'b>(&self, inbox: &'b [f64]) -> Vec<(Walk, &'b [f64])> {
        assert!(
            inbox.len().is_multiple_of(self.record),
            "hand-over: {} numbers are no whole records of {}",
            inbox.len(),
            self.record
        );
        inbox
            .chunks_exact(self.record)
            .map(|r| (self.walk_of(r), &r[WALK_LEN..]))
            .collect()
    }

    /// The walk of one record, in local indices.
    fn walk_of(&self, r: &[f64]) -> Walk {
        let mesh = &self.local.mesh;
        let pos = r[10] as u32;
        let cell = *self
            .own
            .get(&pos)
            .unwrap_or_else(|| panic!("hand-over: cell {pos} is not on this rank"));
        let face_id = r[11] as u64;
        let face = mesh
            .cell_face_iter(cell as usize)
            .map(|(f, _)| f)
            .find(|&f| mesh.face_id[f] == face_id)
            .unwrap_or_else(|| panic!("hand-over: face {face_id} is not a face of cell {pos}"));
        Walk {
            origin: [r[0], r[1], r[2]],
            target: [r[3], r[4], r[5]],
            lambda: r[6],
            x: [r[7], r[8], r[9]],
            cell,
            face: face as u32,
            edge: r[12] as u32,
            entered: r[13] as u32,
            steps: r[14] as u32,
            crossed: r[15] as u32,
        }
    }

    /// The exchange for all ranks of a decomposition held in one process
    /// (tests and tools): `outboxes[r]` is the outbox of rank `r`; returns
    /// the inbox of every rank, in the order [`HandOver::exchange`] gives.
    pub fn exchange_all(all: &[HandOver<'_>], outboxes: &mut [Vec<Vec<f64>>]) -> Vec<Vec<f64>> {
        let mut inboxes = vec![Vec::new(); all.len()];
        for (r, h) in all.iter().enumerate() {
            for nb in &h.local.halo.neighbours {
                let back = all[nb.rank]
                    .local
                    .halo
                    .neighbours
                    .iter()
                    .position(|x| x.rank == r)
                    .expect("halo plans are symmetric");
                inboxes[r].extend_from_slice(&outboxes[nb.rank][back]);
            }
        }
        for o in outboxes.iter_mut() {
            for b in o.iter_mut() {
                b.clear();
            }
        }
        inboxes
    }
}
