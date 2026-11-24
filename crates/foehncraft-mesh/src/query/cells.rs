//! Bounding volume hierarchy over cell boxes (D108 part 2): the search for
//! the cell that holds a point (probes, line samples, the pressure reference
//! cell, injectors, mesh-to-mesh interpolation) and for the cells whose
//! boxes overlap a box (donor candidates).
//!
//! Each cell's box is the box of its nodes, widened by a relative margin
//! [`PAD`] of its largest side, so that a point on a face, which the
//! containment test of [`crate::locate`] accepts with a relative
//! tolerance, also lies in the box. The build is the binned
//! surface-area heuristic of [`super::bvh`] (MacDonald & Booth 1990; Wald
//! 2007) over the box centroids, so a strongly graded mesh (boundary layers,
//! refinement zones) costs O(log n) per query like a uniform one, which the
//! former uniform grid of boxes, sized by the mean cell, did not (see
//! `docs/assessments/graphics-algorithms-for-modules.md`, section 3.1).
//!
//! **Determinism (D93).** The build orders with total orders only. A query
//! visits every leaf whose box holds the point and reports the candidate
//! cells; [`crate::locate::Locator::locate_search`] keeps the smallest cell
//! index that contains the point, so the answer does not depend on the
//! traversal order or the tree shape.

// coordinates indexed by axis read more clearly than zipped iterators
#![allow(clippy::needless_range_loop)]

use super::bvh::{Aabb, Node, build_nodes};
use crate::mesh::Mesh;

/// Relative margin of a cell box (of its largest side).
pub const PAD: f64 = 1e-9;

/// Hierarchy over the boxes of the local cells of a mesh.
#[derive(Debug, Clone)]
pub struct CellTree {
    nodes: Vec<Node>,
    /// Cells in leaf order.
    order: Vec<u32>,
    /// Box of each cell (with the margin), by cell index.
    boxes: Vec<Aabb>,
}

impl CellTree {
    /// Builds the hierarchy over the local cells (`0..mesh.n_cells`).
    pub fn build(mesh: &Mesh) -> Self {
        let n = mesh.n_cells;
        let mut boxes = Vec::with_capacity(n);
        let mut cent = Vec::with_capacity(n);
        for c in 0..n {
            let mut b = Aabb::EMPTY;
            for v in mesh.cell_nodes(c) {
                b.grow(mesh.points[v as usize]);
            }
            let side = (0..3).map(|d| b.hi[d] - b.lo[d]).fold(0.0, f64::max);
            let pad = PAD * side;
            for d in 0..3 {
                b.lo[d] -= pad;
                b.hi[d] += pad;
            }
            cent.push(std::array::from_fn(|d| 0.5 * (b.lo[d] + b.hi[d])));
            boxes.push(b);
        }
        let (nodes, order) = build_nodes(&cent, &boxes);
        CellTree {
            nodes,
            order,
            boxes,
        }
    }

    /// Number of cells.
    pub fn len(&self) -> usize {
        self.boxes.len()
    }

    /// `true` for a mesh without cells.
    pub fn is_empty(&self) -> bool {
        self.boxes.is_empty()
    }

    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// The box (lower and upper corner, with the margin) of cell `c`.
    pub fn cell_box(&self, c: usize) -> ([f64; 3], [f64; 3]) {
        (self.boxes[c].lo, self.boxes[c].hi)
    }

    /// Calls `f` for every cell whose box holds `x` (boxes closed), in leaf
    /// order. The caller decides containment; use the smallest index among
    /// the cells that contain `x` for a deterministic answer.
    pub fn for_each_candidate(&self, x: [f64; 3], mut f: impl FnMut(usize)) {
        self.walk(|b| holds(b, x), |c| f(c));
    }

    /// Calls `f` for every cell whose box overlaps the box `[lo, hi]`
    /// (closed boxes, touching counts), in leaf order.
    pub fn for_each_overlap(&self, lo: [f64; 3], hi: [f64; 3], mut f: impl FnMut(usize)) {
        let q = Aabb { lo, hi };
        self.walk(|b| overlaps(b, &q), |c| f(c));
    }

    /// The cells whose boxes hold `x`, ascending.
    pub fn candidates(&self, x: [f64; 3]) -> Vec<usize> {
        let mut v = Vec::new();
        self.for_each_candidate(x, |c| v.push(c));
        v.sort_unstable();
        v
    }

    fn walk(&self, test: impl Fn(&Aabb) -> bool, mut f: impl FnMut(usize)) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack: Vec<u32> = Vec::with_capacity(64);
        stack.push(0);
        while let Some(i) = stack.pop() {
            let nd = self.nodes[i as usize];
            if !test(&nd.bx) {
                continue;
            }
            if nd.count > 0 {
                let s = nd.first as usize;
                for &c in &self.order[s..s + nd.count as usize] {
                    if test(&self.boxes[c as usize]) {
                        f(c as usize);
                    }
                }
            } else {
                stack.push(nd.first + 1);
                stack.push(nd.first);
            }
        }
    }
}

#[inline]
fn holds(b: &Aabb, x: [f64; 3]) -> bool {
    (0..3).all(|d| x[d] >= b.lo[d] && x[d] <= b.hi[d])
}

#[inline]
fn overlaps(b: &Aabb, q: &Aabb) -> bool {
    (0..3).all(|d| b.lo[d] <= q.hi[d] && q.lo[d] <= b.hi[d])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate::{BoxSpec, Grading, Kind};
    use foehncraft_core::rng;

    fn brute_candidates(t: &CellTree, x: [f64; 3]) -> Vec<usize> {
        (0..t.len()).filter(|&c| holds(&t.boxes[c], x)).collect()
    }

    #[test]
    fn candidates_and_overlaps_agree_with_brute_force() {
        for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
            let mut s = BoxSpec::new(kind, 3, [9, 7, 6]);
            s.perturb = 0.25;
            s.grading.y = Grading {
                ratio: 200.0,
                symmetric: true,
            };
            let m = s.generate().unwrap();
            let t = CellTree::build(&m);
            assert_eq!(t.len(), m.n_cells);
            for i in 0..500u64 {
                let x: [f64; 3] =
                    std::array::from_fn(|d| -0.1 + 1.2 * rng::uniform([5, i, d as u64]));
                assert_eq!(t.candidates(x), brute_candidates(&t, x), "{kind:?} {x:?}");
                let h = 0.05 * rng::uniform([6, i, 0]);
                let lo = x.map(|v| v - h);
                let hi = x.map(|v| v + h);
                let mut got = Vec::new();
                t.for_each_overlap(lo, hi, |c| got.push(c));
                got.sort_unstable();
                let q = Aabb { lo, hi };
                let want: Vec<usize> = (0..t.len())
                    .filter(|&c| overlaps(&t.boxes[c], &q))
                    .collect();
                assert_eq!(got, want, "{kind:?} overlap");
            }
        }
    }

    #[test]
    fn build_is_reproducible() {
        let mut s = BoxSpec::new(Kind::Mixed, 3, [8, 8, 8]);
        s.perturb = 0.3;
        let m = s.generate().unwrap();
        let a = CellTree::build(&m);
        let b = CellTree::build(&m);
        assert_eq!(a.order, b.order);
        assert_eq!(a.nodes.len(), b.nodes.len());
        for (x, y) in a.nodes.iter().zip(&b.nodes) {
            assert_eq!((x.first, x.count), (y.first, y.count));
            assert_eq!(x.bx, y.bx);
        }
    }
}
