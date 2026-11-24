//! The multigrid hierarchy distributed over ranks (MPI step W4): the
//! setup is replicated, the cycle distributed.
//!
//! Every rank gathers the fine matrix in global order and builds the whole
//! hierarchy, which is therefore the hierarchy of the one-rank run
//! (D93 (5)); then [`distribute`] keeps the rows each rank owns and plans
//! the exchanges of the cycle. Nothing of the plan is communicated: every
//! rank knows the hierarchy and the partition, so it can work out what
//! every other rank needs from it.
//!
//! - **Ownership.** A fine row (level 0) belongs to the rank of its cell; a
//!   coarse row i of level l + 1 to the owner of the first fine row of
//!   R_l row i (the lowest index of its aggregate).
//! - **Local layout of a level.** The own rows first (level 0 in the cell
//!   order of the rank's local mesh, coarser levels ascending), then the
//!   ghost entries the rank reads, grouped by owner, ascending.
//! - **Ghosts of level l** are the columns of the own rows of A_l (the
//!   smoother and the residual read x_l), of the own rows of R_l (the
//!   restriction reads w_l) and of the own rows of P_{l−1} (the
//!   prolongation reads x_l). The coarsest, dense level is gathered whole.
//!
//! Every own row keeps the entries of the global row in their global
//! order, only the column numbers change, so a row product adds the same
//! terms in the same order as on one rank and the cycle is bitwise the
//! one-rank cycle.

use super::csr::Csr;
use super::setup::Hierarchy;

/// What one rank sends to and receives from one other rank in the exchange
/// of a level.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Peer {
    /// The other rank.
    pub rank: usize,
    /// Local indices of own entries it reads, in its receive order.
    pub send: Vec<u32>,
    /// Local indices of the ghost entries it owns.
    pub recv: Vec<u32>,
}

/// The vector layout and the exchange of one level on one rank.
#[derive(Debug, Clone, Default)]
pub struct LevelPlan {
    /// Global index of every local entry: the own entries first.
    pub global: Vec<u32>,
    /// Number of own entries.
    pub n_own: usize,
    /// The ranks this rank exchanges with, ascending (both directions,
    /// a list may be empty).
    pub peers: Vec<Peer>,
}

/// One level of the hierarchy on one rank: the own rows with local
/// column numbers.
#[derive(Debug, Clone)]
pub struct DistLevel {
    /// Own rows of A_l, columns in the layout of level l.
    pub a: Csr,
    /// 1 / diag(A_l) of the own rows.
    pub dinv: Vec<f64>,
    /// Own rows of P_l (own rows of level l), columns in the layout of
    /// level l + 1.
    pub p: Csr,
    /// Own rows of R_l (own rows of level l + 1), columns in the layout of
    /// level l.
    pub r: Csr,
    /// Smoother steps.
    pub steps: Vec<(f64, f64)>,
}

/// The hierarchy of one rank.
#[derive(Debug, Clone)]
pub struct DistHierarchy {
    /// Levels 0 .. L−1.
    pub levels: Vec<DistLevel>,
    /// Layouts of levels 0 ..= L (the coarsest last: every entry, own or
    /// ghost).
    pub plans: Vec<LevelPlan>,
    /// Size of the coarsest level.
    pub coarse_n: usize,
    /// Dense inverse of the coarsest operator, row-major (replicated).
    pub coarse_inv: Vec<f64>,
}

/// The owner of every row of every level (levels 0 ..= L).
fn owners(h: &Hierarchy, fine_owner: &[u32]) -> Vec<Vec<u32>> {
    let mut out = vec![fine_owner.to_vec()];
    for lv in &h.levels {
        let prev = out.last().expect("level 0");
        let next = (0..lv.r.n_rows)
            .map(|i| {
                let (j, _) = lv.r.row(i).next().expect("an aggregate has a fine row");
                prev[j]
            })
            .collect();
        out.push(next);
    }
    out
}

/// The global columns of the own rows of `m` (rows `rows`).
fn columns(m: &Csr, rows: &[u32], out: &mut Vec<u32>) {
    for &i in rows {
        out.extend(m.row(i as usize).map(|(j, _)| j as u32));
    }
}

/// The ghost entries `rank` reads at every level, ascending, without its
/// own rows.
fn ghosts(h: &Hierarchy, own: &[Vec<u32>], owner: &[Vec<u32>], rank: u32) -> Vec<Vec<u32>> {
    let nl = h.levels.len();
    (0..=nl)
        .map(|l| {
            let mut g: Vec<u32> = Vec::new();
            if l == nl {
                // the coarsest level is gathered whole
                g.extend(0..h.coarse_n as u32);
            } else {
                columns(&h.levels[l].a, &own[l], &mut g);
                columns(&h.levels[l].r, &own[l + 1], &mut g);
            }
            if l > 0 {
                columns(&h.levels[l - 1].p, &own[l - 1], &mut g);
            }
            g.sort_unstable();
            g.dedup();
            g.retain(|&j| owner[l][j as usize] != rank);
            // grouped by owner, ascending within the group
            g.sort_by_key(|&j| (owner[l][j as usize], j));
            g
        })
        .collect()
}

/// The own rows of every rank at every level: level 0 for `rank` in
/// `fine_order`, every other rank and level ascending.
fn own_rows(owner: &[Vec<u32>], size: usize) -> Vec<Vec<Vec<u32>>> {
    let mut out = vec![vec![Vec::new(); owner.len()]; size];
    for (l, ow) in owner.iter().enumerate() {
        for (i, &r) in ow.iter().enumerate() {
            out[r as usize][l].push(i as u32);
        }
    }
    out
}

/// Own rows of `m` (global rows `rows`, in that order) with the columns
/// renumbered by `local` (global → local index of the column level); the
/// entries keep their global order.
fn local_rows(m: &Csr, rows: &[u32], local: &[u32], n_cols: usize) -> Csr {
    let mut row_ptr = Vec::with_capacity(rows.len() + 1);
    row_ptr.push(0);
    let (mut col, mut val) = (Vec::new(), Vec::new());
    for &i in rows {
        for (j, v) in m.row(i as usize) {
            let lj = local[j];
            debug_assert!(lj != u32::MAX, "column {j} is not in the layout");
            col.push(lj);
            val.push(v);
        }
        row_ptr.push(col.len());
    }
    Csr {
        n_rows: rows.len(),
        n_cols,
        row_ptr,
        col,
        val,
    }
}

/// The part of the replicated hierarchy `h` that `rank` of `size` ranks
/// owns, with the exchange plans. `fine_owner` is the rank of every fine
/// row (the partition, by global position), `fine_order` the global
/// positions of the rank's own cells in its local order.
pub fn distribute(
    h: &Hierarchy,
    fine_owner: &[u32],
    fine_order: &[u32],
    rank: usize,
    size: usize,
) -> DistHierarchy {
    let nl = h.levels.len();
    let owner = owners(h, fine_owner);
    let mut own = own_rows(&owner, size);
    // level 0 of this rank in the order of its local mesh
    debug_assert_eq!(own[rank][0].len(), fine_order.len());
    own[rank][0] = fine_order.to_vec();
    let all_ghosts: Vec<Vec<Vec<u32>>> = (0..size)
        .map(|q| ghosts(h, &own[q], &owner, q as u32))
        .collect();
    let n_rows = |l: usize| {
        if l < nl {
            h.levels[l].a.n_rows
        } else {
            h.coarse_n
        }
    };
    let mut plans = Vec::with_capacity(nl + 1);
    let mut local_of: Vec<Vec<u32>> = Vec::with_capacity(nl + 1);
    for l in 0..=nl {
        let mine = &own[rank][l];
        let gh = &all_ghosts[rank][l];
        let global: Vec<u32> = mine.iter().chain(gh.iter()).copied().collect();
        let mut local = vec![u32::MAX; n_rows(l)];
        for (k, &g) in global.iter().enumerate() {
            local[g as usize] = k as u32;
        }
        let mut peers = Vec::new();
        for (q, theirs) in all_ghosts.iter().enumerate() {
            if q == rank {
                continue;
            }
            // what q reads from me, in q's order; what I read from q
            let send: Vec<u32> = theirs[l]
                .iter()
                .filter(|&&j| owner[l][j as usize] == rank as u32)
                .map(|&j| local[j as usize])
                .collect();
            let recv: Vec<u32> = gh
                .iter()
                .filter(|&&j| owner[l][j as usize] == q as u32)
                .map(|&j| local[j as usize])
                .collect();
            if !send.is_empty() || !recv.is_empty() {
                peers.push(Peer {
                    rank: q,
                    send,
                    recv,
                });
            }
        }
        plans.push(LevelPlan {
            n_own: mine.len(),
            global,
            peers,
        });
        local_of.push(local);
    }
    let levels = (0..nl)
        .map(|l| {
            let lv = &h.levels[l];
            let rows = &own[rank][l];
            let next = &own[rank][l + 1];
            DistLevel {
                a: local_rows(&lv.a, rows, &local_of[l], plans[l].global.len()),
                dinv: rows.iter().map(|&i| lv.dinv[i as usize]).collect(),
                p: local_rows(&lv.p, rows, &local_of[l + 1], plans[l + 1].global.len()),
                r: local_rows(&lv.r, next, &local_of[l], plans[l].global.len()),
                steps: lv.steps.clone(),
            }
        })
        .collect();
    DistHierarchy {
        levels,
        plans,
        coarse_n: h.coarse_n,
        coarse_inv: h.coarse_inv.clone(),
    }
}
