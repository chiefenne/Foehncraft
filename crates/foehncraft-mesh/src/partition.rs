//! Domain decomposition (D32).
//!
//! - `metis`: multilevel k-way or recursive bisection of the cell graph
//!   (cells connected through interior faces) with METIS 5 (Karypis &
//!   Kumar, *A fast and high quality multilevel scheme for partitioning
//!   irregular graphs*, SIAM J. Sci. Comput. 20 (1998) 359–392), linked as a
//!   C library built from source (D9; feature `metis`). Recursive bisection
//!   up to 8 parts, k-way above.
//! - `rcb`: recursive coordinate bisection of the cell centroids (Berger &
//!   Bokhari, *A partitioning strategy for nonuniform problems on
//!   multiprocessors*, IEEE Trans. Comput. C-36 (1987) 570–580): always
//!   available, no C code, perfectly balanced, more cut faces.
//!
//! Both are deterministic: the same mesh and part count give the same
//! partition on every rank and in every run (D31).

use crate::geometry::Geometry;
use crate::mesh::Mesh;
use serde::{Deserialize, Serialize};

/// Partitioning method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "kebab-case")]
pub enum Method {
    /// METIS (feature `metis`).
    #[default]
    Metis,
    /// Recursive coordinate bisection.
    Rcb,
}

impl Method {
    /// Name used in case files.
    pub fn name(self) -> &'static str {
        match self {
            Self::Metis => "metis",
            Self::Rcb => "rcb",
        }
    }
}

/// Splits the local cells of `mesh` into `parts` parts: the part of each
/// cell.
pub fn partition(
    mesh: &Mesh,
    geo: &Geometry,
    parts: usize,
    method: Method,
) -> Result<Vec<u32>, String> {
    if parts == 0 {
        return Err("at least one part".into());
    }
    if parts > mesh.n_cells {
        return Err(format!("{parts} parts for {} cells", mesh.n_cells));
    }
    if parts == 1 {
        return Ok(vec![0; mesh.n_cells]);
    }
    match method {
        Method::Rcb => Ok(rcb(&geo.cell_centroid[..mesh.n_cells], parts)),
        Method::Metis => metis_partition(mesh, parts),
    }
}

#[cfg(feature = "metis")]
fn metis_partition(mesh: &Mesh, parts: usize) -> Result<Vec<u32>, String> {
    let n = mesh.n_cells;
    let mut xadj: Vec<metis::Idx> = Vec::with_capacity(n + 1);
    let mut adjncy: Vec<metis::Idx> = Vec::with_capacity(2 * mesh.n_inner_faces);
    xadj.push(0);
    for c in 0..n {
        let start = adjncy.len();
        for (f, _) in mesh.cell_face_iter(c) {
            if f < mesh.n_inner_faces {
                let o = mesh.across(f, c).unwrap() as metis::Idx;
                // one edge per neighbour, also where two cells share several
                // faces (merged cut cells)
                if !adjncy[start..].contains(&o) {
                    adjncy.push(o);
                }
            }
        }
        xadj.push(adjncy.len() as metis::Idx);
    }
    let mut part: Vec<metis::Idx> = vec![0; n];
    let g = metis::Graph::new(1, parts as metis::Idx, &xadj, &adjncy)
        .map_err(|e| format!("METIS: {e}"))?;
    let r = if parts <= 8 {
        g.part_recursive(&mut part)
    } else {
        g.part_kway(&mut part)
    };
    r.map_err(|e| format!("METIS: {e}"))?;
    Ok(part.into_iter().map(|p| p as u32).collect())
}

#[cfg(not(feature = "metis"))]
fn metis_partition(_mesh: &Mesh, _parts: usize) -> Result<Vec<u32>, String> {
    Err("built without METIS (feature `metis`); use method = \"rcb\"".into())
}

/// Recursive coordinate bisection: splits along the longest extent at the
/// position that gives each side its share of the parts.
pub fn rcb(centroids: &[[f64; 3]], parts: usize) -> Vec<u32> {
    let mut part = vec![0u32; centroids.len()];
    let mut idx: Vec<u32> = (0..centroids.len() as u32).collect();
    split(centroids, &mut idx, 0, parts, &mut part);
    part
}

fn split(x: &[[f64; 3]], idx: &mut [u32], first: usize, parts: usize, out: &mut [u32]) {
    if parts == 1 {
        for &c in idx.iter() {
            out[c as usize] = first as u32;
        }
        return;
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &c in idx.iter() {
        for d in 0..3 {
            lo[d] = lo[d].min(x[c as usize][d]);
            hi[d] = hi[d].max(x[c as usize][d]);
        }
    }
    let d = (0..3)
        .max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b])))
        .unwrap();
    let left = parts / 2;
    let k = idx.len() * left / parts;
    idx.sort_unstable_by(|&a, &b| {
        x[a as usize][d]
            .total_cmp(&x[b as usize][d])
            .then(a.cmp(&b))
    });
    let (l, r) = idx.split_at_mut(k);
    split(x, l, first, left, out);
    split(x, r, first + left, parts - left, out);
}

/// Statistics of a partition.
#[derive(Debug, Clone, Serialize)]
pub struct PartitionStats {
    /// Number of parts.
    pub parts: usize,
    /// Cells per part.
    pub cells: Vec<usize>,
    /// Largest part over the mean part size.
    pub imbalance: f64,
    /// Interior faces between different parts.
    pub cut_faces: usize,
    /// Neighbouring parts per part.
    pub neighbours: Vec<usize>,
}

impl PartitionStats {
    /// Statistics of `part` on `mesh`.
    pub fn compute(mesh: &Mesh, part: &[u32], parts: usize) -> Self {
        let mut cells = vec![0usize; parts];
        for &p in part {
            cells[p as usize] += 1;
        }
        let mut cut = 0;
        let mut pairs: Vec<(u32, u32)> = Vec::new();
        for f in 0..mesh.n_inner_faces {
            let (a, b) = (
                part[mesh.owner[f] as usize],
                part[mesh.neighbour[f] as usize],
            );
            if a != b {
                cut += 1;
                pairs.push((a, b));
                pairs.push((b, a));
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        let mut neighbours = vec![0usize; parts];
        for (a, _) in pairs {
            neighbours[a as usize] += 1;
        }
        let mean = part.len() as f64 / parts as f64;
        Self {
            parts,
            imbalance: *cells.iter().max().unwrap() as f64 / mean,
            cells,
            cut_faces: cut,
            neighbours,
        }
    }
}
