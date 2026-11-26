//! Cell renumbering.
//!
//! The order of the cells decides how well the face and cell loops use the
//! caches and how coalesced GPU loads are: gate G1 measured three to five
//! times less bandwidth for a random order than for a lexicographic one.
//! Meshes from a mesh generator come in whatever order it produced.
//!
//! - [`rcm`]: reverse Cuthill–McKee (Cuthill & McKee, *Reducing the
//!   bandwidth of sparse symmetric matrices*, Proc. ACM National Conference
//!   1969; reversed as proposed by George, 1971): breadth-first numbering
//!   from a peripheral cell, neighbours by increasing degree, then reversed.
//!   Small matrix bandwidth and good locality.
//! - [`shuffled`]: a random order, the worst case, for benchmarks.
//!
//! [`apply`] builds the renumbered mesh; the faces are re-sorted, so the
//! upper-triangular face order of [`Mesh`] holds again.

use crate::builder::MeshBuilder;
use crate::mesh::Mesh;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;

/// Cell order.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Ordering {
    /// As generated or read.
    #[default]
    None,
    /// Reverse Cuthill–McKee.
    Rcm,
    /// Random permutation.
    Shuffle,
}

/// New number of each cell (`perm[old] = new`) by reverse Cuthill–McKee.
pub fn rcm(mesh: &Mesh) -> Vec<u32> {
    let order = rcm_order(mesh.n_cells, |c, out| {
        out.extend(
            mesh.cell_face_iter(c)
                .filter(|&(f, _)| f < mesh.n_inner_faces)
                .map(|(f, _)| mesh.across(f, c).unwrap()),
        )
    });
    let mut perm = vec![0u32; mesh.n_cells];
    for (new, &old) in order.iter().enumerate() {
        perm[old] = new as u32;
    }
    perm
}

/// Reverse Cuthill–McKee order of a graph of `n` vertices: `order[new] =
/// old`. `adj(v, out)` appends the neighbours of `v` (vertices `< n`).
pub fn rcm_order(n: usize, adj: impl Fn(usize, &mut Vec<usize>)) -> Vec<usize> {
    let mut nb = Vec::new();
    let degree: Vec<usize> = (0..n)
        .map(|c| {
            nb.clear();
            adj(c, &mut nb);
            nb.len()
        })
        .collect();
    let mut order: Vec<usize> = Vec::with_capacity(n);
    let mut seen = vec![false; n];
    let mut level = vec![usize::MAX; n];
    let mut queue = VecDeque::new();
    // last vertex of a breadth-first search from `s` (deepest level, lowest
    // degree), over unseen vertices
    let bfs_last = |s: usize, seen: &[bool], level: &mut Vec<usize>| -> usize {
        let mut visited = vec![s];
        level[s] = 0;
        let mut q = VecDeque::from([s]);
        let mut last = (s, 0usize);
        let mut nb = Vec::new();
        while let Some(c) = q.pop_front() {
            nb.clear();
            adj(c, &mut nb);
            for &d in &nb {
                if level[d] == usize::MAX && !seen[d] {
                    level[d] = level[c] + 1;
                    visited.push(d);
                    if level[d] > last.1 || (level[d] == last.1 && degree[d] < degree[last.0]) {
                        last = (d, level[d]);
                    }
                    q.push_back(d);
                }
            }
        }
        for v in visited {
            level[v] = usize::MAX;
        }
        last.0
    };
    for start in 0..n {
        if seen[start] {
            continue;
        }
        // pseudo-peripheral start: two sweeps
        let a = bfs_last(start, &seen, &mut level);
        let s = bfs_last(a, &seen, &mut level);
        let first = order.len();
        seen[s] = true;
        queue.push_back(s);
        while let Some(c) = queue.pop_front() {
            order.push(c);
            nb.clear();
            adj(c, &mut nb);
            nb.retain(|&d| !seen[d]);
            nb.sort_unstable_by_key(|&d| (degree[d], d));
            nb.dedup();
            for &d in &nb {
                seen[d] = true;
                queue.push_back(d);
            }
        }
        order[first..].reverse();
    }
    order
}

/// A random permutation (`perm[old] = new`), reproducible from `seed`.
pub fn shuffled(n: usize, seed: u64) -> Vec<u32> {
    let mut p: Vec<u32> = (0..n as u32).collect();
    for i in (1..n).rev() {
        let j = (foehncraft_core::rng::keyed([
            foehncraft_core::rng::stream::MESH,
            seed,
            i as u64,
            0x5eed,
        ]) % (i as u64 + 1)) as usize;
        p.swap(i, j);
    }
    p
}

/// The mesh with its cells renumbered (`perm[old] = new`).
pub fn apply(mesh: &Mesh, perm: &[u32]) -> Result<Mesh, String> {
    if mesh.n_halo != 0 {
        return Err("renumber the global mesh before the decomposition".into());
    }
    if mesh.periodic.is_some() {
        return Err("renumber the mesh before pairing its periodic boundaries".into());
    }
    let n = mesh.n_cells;
    let mut b = MeshBuilder::new(mesh.points.clone(), n);
    // a mesh with cells that share several faces (merged cut cells) keeps them
    b.allow_shared_pairs();
    let mut zone = vec![0u32; n];
    for c in 0..n {
        zone[perm[c] as usize] = mesh.cell_zone[c];
    }
    b.set_zones(mesh.zones.clone(), zone);
    for p in &mesh.boundaries {
        b.declare_boundary(&p.name);
    }
    for f in 0..mesh.n_interior_faces() {
        b.add_interior(
            mesh.face(f),
            perm[mesh.owner[f] as usize],
            perm[mesh.neighbour[f] as usize],
        );
    }
    for p in &mesh.boundaries {
        for f in p.range() {
            b.add_boundary(mesh.face(f), perm[mesh.owner[f] as usize], &p.name);
        }
    }
    let (mut m, order) = b.finish_with_order()?;
    // faces were added in the old order, so the insertion index is the old
    // face index; map all identities back to the input order
    m.face_id = order.iter().map(|&f| mesh.face_id[f as usize]).collect();
    m.point_id = mesh.point_id.clone();
    for (&p, &id) in perm.iter().zip(&mesh.cell_id[..n]) {
        m.cell_id[p as usize] = id;
    }
    Ok(m)
}

/// Matrix bandwidth of the cell graph: the largest |owner − neighbour|.
pub fn bandwidth(mesh: &Mesh) -> usize {
    (0..mesh.n_inner_faces)
        .map(|f| (mesh.neighbour[f] - mesh.owner[f]) as usize)
        .max()
        .unwrap_or(0)
}

/// Renumbers the cells of `mesh` in place by `ordering`.
pub fn renumber(mesh: &mut Mesh, ordering: Ordering, seed: u64) -> Result<(), String> {
    let perm = match ordering {
        Ordering::None => return Ok(()),
        Ordering::Rcm => rcm(mesh),
        Ordering::Shuffle => shuffled(mesh.n_cells, seed),
    };
    *mesh = apply(mesh, &perm)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate::{BoxSpec, Kind};
    use crate::geometry::Geometry;

    #[test]
    fn rcm_restores_a_small_bandwidth() {
        let m = BoxSpec::new(Kind::Poly, 3, [8, 7, 6]).generate().unwrap();
        let vol: f64 = Geometry::compute(&m).cell_volume.iter().sum();
        let s = apply(&m, &shuffled(m.n_cells, 3)).unwrap();
        s.validate().unwrap();
        let r = apply(&s, &rcm(&s)).unwrap();
        r.validate().unwrap();
        assert!(
            bandwidth(&s) > 4 * bandwidth(&r),
            "{} {}",
            bandwidth(&s),
            bandwidth(&r)
        );
        let vol_r: f64 = Geometry::compute(&r).cell_volume.iter().sum();
        assert!((vol - vol_r).abs() < 1e-12);
        assert_eq!(r.boundaries, m.boundaries);
        // identities still name the input cells and faces
        let g = Geometry::compute(&r);
        let g0 = Geometry::compute(&m);
        for c in 0..r.n_cells {
            assert!((g.cell_volume[c] - g0.cell_volume[r.cell_id[c] as usize]).abs() < 1e-15);
        }
        for f in 0..r.n_faces() {
            let a = g.face_area[f];
            let b = g0.face_area[r.face_id[f] as usize];
            assert!((0..3).all(|d| (a[d].abs() - b[d].abs()).abs() < 1e-15));
        }
    }
}
