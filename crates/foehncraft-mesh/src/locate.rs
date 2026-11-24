//! Point location and cell-to-cell particle tracking (D47).
//!
//! Every cell is split into tetrahedra (as in [`crate::geometry`]): an
//! interior point of the cell (the apex: the average of its face centroids),
//! the fan point of a face (its node average) and one edge of that face. A
//! tetrahedron is named by (cell, face, edge), where the edge index counts
//! in the stored node order of the face, so the two cells of a face see the
//! same triangles: the decomposition is watertight across faces, also for
//! warped faces of polyhedra.
//!
//! A particle moves along a straight segment from tetrahedron to
//! tetrahedron: in the current tetrahedron the segment leaves through the
//! face with the smallest crossing parameter, computed from signed volumes,
//! and continues in the tetrahedron on the other side, which is in the same
//! cell or, through a mesh face, in the neighbouring cell. The method follows
//! the tetrahedral decomposition tracking of Macpherson, Nordin & Weller,
//! *Particle tracking in unstructured, arbitrary polyhedral meshes for use in
//! CFD and molecular dynamics*, Commun. Numer. Meth. Engng 25 (2009) 263–273.
//! A track ends at the target, at a boundary face (the caller applies the
//! wall or outlet model) or at a processor face (the particle is handed to
//! the neighbouring rank, so a halo of one cell layer suffices).
//!
//! The steps are the walk body [`foehncraft_kernels::walk`], one source for
//! the CPU and the GPU (D259) and shared by particles and radiation bundles
//! (D257). [`Locator::walk`] also reports every cell the segment passes,
//! with the interval of the segment parameter inside it (the exact path
//! length per cell); a walk stopped at a processor face goes on on the rank
//! behind it through [`crate::handover`], bitwise as on one rank (D258).
//!
//! Point location walks from a hint cell to the point; if the walk is
//! blocked (by a boundary of a non-convex domain), or no hint is known, a
//! bounding volume hierarchy over the cell boxes
//! ([`crate::query::CellTree`], D108 part 2) finds the cell: O(log n) per
//! query also on strongly graded meshes. The former uniform grid of cell
//! boxes is kept as [`Locator::locate_search_grid`], the reference of the
//! tests and of the benchmark `examples/locate_bench.rs`.

use crate::geometry::{Geometry, add, node_average, scale, sub};
use crate::mesh::Mesh;
use crate::query::CellTree;
pub use foehncraft_kernels::walk::{Walk, WalkEnd, WalkMesh};
use foehncraft_kernels::walk::{inner, sides, walk};
use std::sync::OnceLock;

/// A step limit for walks on a mesh of `n_cells` cells: 64 + 16 tetrahedra
/// per cell, at most 16 × 2²⁰ + 64. Give it the global cell count, so the
/// limit is the same on every rank.
pub fn step_limit(n_cells: usize) -> u32 {
    (64 + 16 * n_cells.min(1 << 20)) as u32
}

/// A position inside the mesh: the point and the tetrahedron that holds it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Position {
    /// Coordinates.
    pub x: [f64; 3],
    /// Cell.
    pub cell: usize,
    /// Face of the cell's tetrahedron.
    pub face: usize,
    /// Edge of that face (in the stored node order).
    pub edge: usize,
}

/// How a track ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrackEnd {
    /// The target was reached.
    Reached,
    /// The segment left the domain through boundary face `face`.
    Boundary {
        /// The face.
        face: usize,
    },
    /// The segment reached processor face `face`: hand the particle to the
    /// rank that owns the halo cell behind it.
    Processor {
        /// The face.
        face: usize,
    },
    /// The walk did not converge (degenerate geometry).
    Lost,
}

/// Result of a track.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Track {
    /// Where the particle is at the end (on the face for `Boundary` and
    /// `Processor`).
    pub position: Position,
    /// How the track ended.
    pub end: TrackEnd,
    /// Fraction of the segment travelled.
    pub fraction: f64,
    /// Number of mesh faces crossed.
    pub faces_crossed: usize,
}

/// Point location and tracking on one mesh.
pub struct Locator<'a> {
    mesh: &'a Mesh,
    apex: Vec<[f64; 3]>,
    fan: Vec<[f64; 3]>,
    tree: OnceLock<CellTree>,
    bins: OnceLock<Bins>,
}

struct Bins {
    lo: [f64; 3],
    h: [f64; 3],
    n: [usize; 3],
    ptr: Vec<u32>,
    cells: Vec<u32>,
}

impl<'a> Locator<'a> {
    /// Prepares the decomposition (the fan points of the faces and the
    /// apexes of the local cells).
    pub fn new(mesh: &'a Mesh, geo: &Geometry) -> Self {
        let fan = (0..mesh.n_faces())
            .map(|f| node_average(&mesh.points, mesh.face(f)))
            .collect();
        let apex = (0..mesh.n_cells)
            .map(|c| crate::geometry::cell_apex(mesh, &geo.face_centroid, c))
            .collect();
        Self {
            mesh,
            apex,
            fan,
            tree: OnceLock::new(),
            bins: OnceLock::new(),
        }
    }

    /// The mesh.
    pub fn mesh(&self) -> &Mesh {
        self.mesh
    }

    /// The arrays of the walk body: the mesh and its tetrahedral
    /// decomposition (see [`foehncraft_kernels::walk`]).
    pub fn walk_mesh(&self) -> WalkMesh<'_> {
        let m = self.mesh;
        WalkMesh {
            points: &m.points,
            face_ptr: &m.face_ptr,
            face_nodes: &m.face_nodes,
            owner: &m.owner,
            neighbour: &m.neighbour,
            cell_ptr: &m.cell_ptr,
            cell_faces: &m.cell_faces,
            apex: &self.apex,
            fan: &self.fan,
            n_cells: m.n_cells,
            periodic_start: m.periodic_start(),
        }
    }

    /// Vertices of tetrahedron (cell, face, edge), with its outer triangle
    /// (fan point, p, q) oriented outward from the cell: (apex, m, p, q).
    #[inline]
    fn tet(&self, cell: usize, face: usize, edge: usize) -> [[f64; 3]; 4] {
        self.walk_mesh().tet(cell, face, edge)
    }

    /// `true` if `x` lies in cell `c` (boundaries included, with a relative
    /// tolerance).
    pub fn contains(&self, c: usize, x: [f64; 3]) -> bool {
        self.find_tet(c, x).is_some()
    }

    fn find_tet(&self, c: usize, x: [f64; 3]) -> Option<(usize, usize)> {
        for (f, _) in self.mesh.cell_face_iter(c) {
            for e in 0..self.mesh.face(f).len() {
                let t = self.tet(c, f, e);
                let s = sides(&t, x);
                let vol = inner(t[1], t[2], t[3], t[0]).abs();
                let tol = -1e-12 * vol;
                if s.iter().all(|&v| v >= tol) {
                    return Some((f, e));
                }
            }
        }
        None
    }

    /// A position inside cell `c` near its apex (in its first tetrahedron).
    pub fn start(&self, c: usize) -> Position {
        let (f, _) = self.mesh.cell_face_iter(c).next().unwrap();
        let t = self.tet(c, f, 0);
        let x = scale(add(add(t[0], t[1]), add(t[2], t[3])), 0.25);
        Position {
            x,
            cell: c,
            face: f,
            edge: 0,
        }
    }

    /// The position of point `x` if it lies in cell `c`.
    pub fn position_in(&self, c: usize, x: [f64; 3]) -> Option<Position> {
        self.find_tet(c, x).map(|(face, edge)| Position {
            x,
            cell: c,
            face,
            edge,
        })
    }

    /// Moves a particle from `from` towards `to` along a straight line.
    pub fn track(&self, from: Position, to: [f64; 3]) -> Track {
        let mut w = Walk::new(
            from.x,
            from.cell as u32,
            from.face as u32,
            from.edge as u32,
            to,
        );
        let end = self.walk(&mut w, step_limit(self.mesh.n_cells), |_, _, _| true);
        let (end, cell, crossed) = match end {
            WalkEnd::Reached => (TrackEnd::Reached, w.cell, w.crossed),
            WalkEnd::Boundary { face } => (
                TrackEnd::Boundary {
                    face: face as usize,
                },
                w.cell,
                w.crossed,
            ),
            WalkEnd::Processor { face } => {
                // the particle stays in the own cell, on the face
                let f = face as usize;
                let (p, q) = (self.mesh.owner[f], self.mesh.neighbour[f]);
                let own = if p == w.cell { q } else { p };
                (TrackEnd::Processor { face: f }, own, w.crossed - 1)
            }
            WalkEnd::Lost | WalkEnd::Stopped => (TrackEnd::Lost, w.cell, w.crossed),
        };
        Track {
            position: Position {
                x: w.x,
                cell: cell as usize,
                face: w.face as usize,
                edge: w.edge as usize,
            },
            end,
            fraction: w.lambda,
            faces_crossed: crossed as usize,
        }
    }

    /// Moves walk `w` towards its target (see [`foehncraft_kernels::walk`]):
    /// `visit(cell, λ₀, λ₁)` gets every cell the segment passes with the
    /// interval of the segment parameter inside it. `max_steps` limits the
    /// tetrahedra of the whole walk ([`step_limit`] of the global cell
    /// count gives the same limit on every rank).
    ///
    /// # Panics
    /// If `w` does not stand in a tetrahedron of a local cell (cell, face
    /// of that cell, edge of that face).
    pub fn walk<F: FnMut(u32, f64, f64) -> bool>(
        &self,
        w: &mut Walk,
        max_steps: u32,
        visit: F,
    ) -> WalkEnd {
        let (c, f) = (w.cell as usize, w.face as usize);
        assert!(
            c < self.mesh.n_cells
                && self.mesh.cell_face_iter(c).any(|(g, _)| g == f)
                && (w.edge as usize) < self.mesh.face(f).len(),
            "walk: ({c}, {f}, {}) is not a tetrahedron of a local cell",
            w.edge
        );
        walk(&self.walk_mesh(), w, max_steps, visit)
    }

    /// The walk from `from` to `to`.
    pub fn start_walk(&self, from: Position, to: [f64; 3]) -> Walk {
        Walk::new(
            from.x,
            from.cell as u32,
            from.face as u32,
            from.edge as u32,
            to,
        )
    }

    /// Cell containing `x`: a walk from cell `hint`, with the search of
    /// [`locate_search`](Self::locate_search) as fallback. `None` if `x`
    /// lies outside the mesh.
    pub fn locate(&self, x: [f64; 3], hint: usize) -> Option<usize> {
        let t = self.track(self.start(hint), x);
        if t.end == TrackEnd::Reached && self.contains(t.position.cell, x) {
            return Some(t.position.cell);
        }
        self.locate_search(x)
    }

    /// Cell containing `x` by the hierarchy over the cell boxes (built on
    /// the first call): the smallest cell index among the cells that
    /// contain `x`, so the answer does not depend on the tree (D93). `None`
    /// if `x` lies outside the mesh.
    pub fn locate_search(&self, x: [f64; 3]) -> Option<usize> {
        let tree = self.cell_tree();
        let mut best = usize::MAX;
        tree.for_each_candidate(x, |c| {
            if c < best && self.contains(c, x) {
                best = c;
            }
        });
        (best != usize::MAX).then_some(best)
    }

    /// The hierarchy over the cell boxes (built on the first call).
    pub fn cell_tree(&self) -> &CellTree {
        self.tree.get_or_init(|| CellTree::build(self.mesh))
    }

    /// Cell containing `x` by the former search over a uniform grid of cell
    /// boxes (bin size twice the mean cell size): the first cell, in index
    /// order, of the bin of `x` that contains it. Kept as the reference of
    /// the tests and of the benchmark; [`locate_search`](Self::locate_search)
    /// gives the same cell.
    pub fn locate_search_grid(&self, x: [f64; 3]) -> Option<usize> {
        let b = self.bins.get_or_init(|| Bins::build(self.mesh));
        let mut idx = [0usize; 3];
        for d in 0..3 {
            let r = (x[d] - b.lo[d]) / b.h[d];
            if !(r >= -1e-9 && r <= b.n[d] as f64 + 1e-9) {
                return None;
            }
            idx[d] = (r.max(0.0) as usize).min(b.n[d] - 1);
        }
        let bin = idx[0] + b.n[0] * (idx[1] + b.n[1] * idx[2]);
        b.cells[b.ptr[bin] as usize..b.ptr[bin + 1] as usize]
            .iter()
            .map(|&c| c as usize)
            .find(|&c| self.contains(c, x))
    }
}

impl Bins {
    fn build(mesh: &Mesh) -> Self {
        let (lo, hi) = mesh.bounding_box();
        let nc = mesh.n_cells.max(1);
        let ext = sub(hi, lo);
        let vol = ext.iter().map(|e| e.max(1e-300)).product::<f64>();
        let h0 = (vol / nc as f64).cbrt() * 2.0;
        let mut n = [1usize; 3];
        let mut h = [0.0; 3];
        for d in 0..3 {
            n[d] = ((ext[d] / h0).ceil() as usize).clamp(1, 1024);
            h[d] = ext[d].max(1e-300) / n[d] as f64;
        }
        let nb = n[0] * n[1] * n[2];
        let mut pairs: Vec<(u32, u32)> = Vec::new();
        for c in 0..mesh.n_cells {
            let mut clo = [f64::INFINITY; 3];
            let mut chi = [f64::NEG_INFINITY; 3];
            for v in mesh.cell_nodes(c) {
                let p = mesh.points[v as usize];
                for d in 0..3 {
                    clo[d] = clo[d].min(p[d]);
                    chi[d] = chi[d].max(p[d]);
                }
            }
            let r =
                |d: usize, x: f64| (((x - lo[d]) / h[d]).floor().max(0.0) as usize).min(n[d] - 1);
            for k in r(2, clo[2])..=r(2, chi[2]) {
                for j in r(1, clo[1])..=r(1, chi[1]) {
                    for i in r(0, clo[0])..=r(0, chi[0]) {
                        pairs.push(((i + n[0] * (j + n[1] * k)) as u32, c as u32));
                    }
                }
            }
        }
        pairs.sort_unstable();
        let mut ptr = vec![0u32; nb + 1];
        for &(b, _) in &pairs {
            ptr[b as usize + 1] += 1;
        }
        for b in 0..nb {
            ptr[b + 1] += ptr[b];
        }
        Bins {
            lo,
            h,
            n,
            ptr,
            cells: pairs.into_iter().map(|p| p.1).collect(),
        }
    }
}
