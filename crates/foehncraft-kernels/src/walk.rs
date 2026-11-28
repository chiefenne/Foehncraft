//! Straight walk through the cells of a mesh (D47, D257): the body of the
//! particle tracker of `foehncraft_mesh::locate` and of the radiation
//! bundles (step 29), one source for the CPU and the GPU (D259).
//!
//! Every cell is split into tetrahedra: the cell's apex (the average of its
//! face centroids), the fan point of one of its faces (the node average) and
//! one edge of that face. A tetrahedron is named by (cell, face, edge), the
//! edge counted in the stored node order of the face, so the two cells of a
//! face see the same triangles and the decomposition is watertight across
//! faces, also for warped faces of polyhedra (Macpherson, Nordin & Weller,
//! Commun. Numer. Meth. Engng 25 (2009) 263–273).
//!
//! A walk follows the straight segment from `origin` to `target`: in the
//! current tetrahedron the segment leaves through the face with the
//! smallest crossing parameter, computed from signed volumes, and goes on
//! in the tetrahedron behind it, in the same cell or, through a mesh face,
//! in the neighbouring cell. Each cell the segment passes is reported once,
//! with the interval `[λ₀, λ₁]` of the segment parameter inside it, so a
//! caller gets the exact path length `(λ₁ − λ₀) |target − origin|` per cell
//! (the track-length estimator of the radiation needs nothing else).
//!
//! [`Walk`] is the whole state. A walk stops at a processor face with the
//! state of the cell behind it; on the rank that owns that cell it goes on
//! with the same operations on the same numbers, so the path and every
//! interval are bitwise those of the walk on one rank (D258), on a
//! canonical local mesh (the geometry of every face and cell bitwise the
//! global one, MPI step W1).
//!
//! The body uses checked indexing: it is called from safe code on the CPU
//! (D202), and the inputs are the arrays of a mesh.

use crate::unpack_face;

/// Tetrahedron face index meaning "none" for [`Walk::entered`].
pub const NONE: u32 = 4;

/// Edge index meaning "not found yet" while a walk crosses convex cells by
/// their planes (set again before a tetrahedron step or the end of a walk).
const NO_EDGE: u32 = u32::MAX;

/// The arrays a walk reads: one (local) mesh and its tetrahedral
/// decomposition.
#[derive(Clone, Copy)]
pub struct WalkMesh<'a> {
    /// Node coordinates.
    pub points: &'a [[f64; 3]],
    /// Start of the node list of each face (length `n_faces + 1`).
    pub face_ptr: &'a [u32],
    /// Nodes of each face, the normal pointing out of the owner.
    pub face_nodes: &'a [u32],
    /// Owner cell of each face.
    pub owner: &'a [u32],
    /// Neighbour cell of each interior face (inner, processor, periodic).
    pub neighbour: &'a [u32],
    /// Start of the face list of each local cell.
    pub cell_ptr: &'a [u32],
    /// Faces of each local cell, packed with their orientation.
    pub cell_faces: &'a [u32],
    /// Apex of each local cell.
    pub apex: &'a [[f64; 3]],
    /// Fan point of each face.
    pub fan: &'a [[f64; 3]],
    /// Number of local cells; a cell index at or above it is a halo cell.
    pub n_cells: usize,
    /// First periodic face (`neighbour.len()` without periodic faces); a
    /// periodic face ends a walk like a boundary face.
    pub periodic_start: usize,
    /// Plane of each face (RA0): the unit normal out of the owner and the
    /// offset, `n · x = d`; empty switches the plane walk off.
    pub plane: &'a [[f64; 4]],
    /// Per local cell: 1 where the cell is convex with planar faces (a
    /// walk crosses it by the planes of its faces, with no tetrahedra), 0
    /// elsewhere; empty switches the plane walk off.
    pub convex: &'a [u8],
}

/// The state of a walk: where it is and where it goes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Walk {
    /// Start of the segment.
    pub origin: [f64; 3],
    /// End of the segment.
    pub target: [f64; 3],
    /// Fraction of the segment travelled.
    pub lambda: f64,
    /// Current point, `origin + λ (target − origin)`.
    pub x: [f64; 3],
    /// Cell of the current tetrahedron (a halo cell after
    /// [`WalkEnd::Processor`]).
    pub cell: u32,
    /// Face of the current tetrahedron.
    pub face: u32,
    /// Edge of that face, in its stored node order.
    pub edge: u32,
    /// Face of the current tetrahedron the walk came through (0 the mesh
    /// face, 1 to 3 the inner faces), or [`NONE`].
    pub entered: u32,
    /// Tetrahedra visited so far.
    pub steps: u32,
    /// Mesh faces crossed so far.
    pub crossed: u32,
}

impl Walk {
    /// A walk from `x`, which lies in tetrahedron (cell, face, edge), to
    /// `target`.
    pub fn new(x: [f64; 3], cell: u32, face: u32, edge: u32, target: [f64; 3]) -> Self {
        Self {
            origin: x,
            target,
            lambda: 0.0,
            x,
            cell,
            face,
            edge,
            entered: NONE,
            steps: 0,
            crossed: 0,
        }
    }
}

/// How a call of [`walk`] ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WalkEnd {
    /// The target was reached (`λ = 1`).
    Reached,
    /// The segment left the domain through boundary (or periodic) face
    /// `face`; the walk stands on the face, in the cell inside.
    Boundary {
        /// The face.
        face: u32,
    },
    /// The segment crossed processor face `face`; the walk stands on the
    /// face, in the halo cell behind it, ready to go on on the rank that
    /// owns that cell.
    Processor {
        /// The face.
        face: u32,
    },
    /// The caller's visit returned `false`.
    Stopped,
    /// The step limit was reached or the decomposition was broken
    /// (degenerate geometry).
    Lost,
}

#[inline(always)]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline(always)]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline(always)]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline(always)]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Six times the signed volume of the tetrahedron (a, b, c, d), positive
/// when `d` lies on the inner side of the outward-oriented triangle
/// (a, b, c).
#[inline(always)]
pub fn inner(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    -dot(cross(sub(b, a), sub(c, a)), sub(d, a))
}

/// Signed inner distances (volumes) of `x` to the four faces of a
/// tetrahedron (a, m, p, q), each face oriented outward: 0 = (m, p, q) on
/// the mesh face, 1 = (a, p, m), 2 = (a, m, q), 3 = (a, q, p). All positive
/// inside.
#[inline(always)]
pub fn sides(t: &[[f64; 3]; 4], x: [f64; 3]) -> [f64; 4] {
    let [a, m, p, q] = *t;
    [
        inner(m, p, q, x),
        inner(a, p, m, x),
        inner(a, m, q, x),
        inner(a, q, p, x),
    ]
}

/// One entry of [`sides`]: the signed inner distance of `x` to face `j`
/// of the tetrahedron `t`, the same operations as in `sides`.
#[inline(always)]
pub fn side(t: &[[f64; 3]; 4], j: u32, x: [f64; 3]) -> f64 {
    let [a, m, p, q] = *t;
    match j {
        0 => inner(m, p, q, x),
        1 => inner(a, p, m, x),
        2 => inner(a, m, q, x),
        _ => inner(a, q, p, x),
    }
}

/// Edge of the fan triangle of face `g` that contains the point `x` on the
/// face (the one with the largest smallest barycentric coordinate).
fn fan_edge(m: &WalkMesh<'_>, g: usize, x: [f64; 3]) -> u32 {
    let nodes = m.face(g);
    let k = nodes.len();
    let f = m.fan[g];
    let mut best = (f64::NEG_INFINITY, 0usize);
    for e in 0..k {
        let p = m.points[nodes[e] as usize];
        let q = m.points[nodes[if e + 1 == k { 0 } else { e + 1 }] as usize];
        let n = cross(sub(p, f), sub(q, f));
        let w0 = dot(cross(sub(p, x), sub(q, x)), n);
        let w1 = dot(cross(sub(q, x), sub(f, x)), n);
        let w2 = dot(cross(sub(f, x), sub(p, x)), n);
        let w = w0.min(w1).min(w2) / dot(n, n).max(f64::MIN_POSITIVE);
        if w > best.0 {
            best = (w, e);
        }
    }
    best.1 as u32
}

/// Puts the tetrahedron state of `w` on the tetrahedron of its cell that
/// holds `w.x` (after steps through the planes of a convex cell, which keep
/// none).
fn settle(m: &WalkMesh<'_>, w: &mut Walk) {
    let cell = w.cell as usize;
    if cell >= m.n_cells {
        return;
    }
    let cf = m.cell_ptr[cell] as usize..m.cell_ptr[cell + 1] as usize;
    let mut best = (f64::NEG_INFINITY, w.face, w.edge);
    for &entry in &m.cell_faces[cf] {
        let (g, _) = unpack_face(entry);
        for e in 0..m.face(g).len() {
            let t = m.tet(cell, g, e);
            let s = sides(&t, w.x);
            let vol = inner(t[1], t[2], t[3], t[0]).abs().max(f64::MIN_POSITIVE);
            let r = s[0].min(s[1]).min(s[2]).min(s[3]) / vol;
            if r > best.0 {
                best = (r, g as u32, e as u32);
            }
        }
    }
    w.face = best.1;
    w.edge = best.2;
    w.entered = NONE;
}

impl WalkMesh<'_> {
    /// Nodes of face `f`.
    #[inline(always)]
    fn face(&self, f: usize) -> &[u32] {
        &self.face_nodes[self.face_ptr[f] as usize..self.face_ptr[f + 1] as usize]
    }

    /// Vertices of tetrahedron (cell, face, edge), its outer triangle
    /// (fan point, p, q) oriented outward from the cell: (apex, m, p, q).
    #[inline(always)]
    pub fn tet(&self, cell: usize, face: usize, edge: usize) -> [[f64; 3]; 4] {
        let nodes = self.face(face);
        let k = nodes.len();
        let mut p = self.points[nodes[edge] as usize];
        let mut q = self.points[nodes[if edge + 1 == k { 0 } else { edge + 1 }] as usize];
        if self.owner[face] as usize != cell {
            core::mem::swap(&mut p, &mut q);
        }
        [self.apex[cell], self.fan[face], p, q]
    }
}

/// Moves `w` towards its target, at most up to `max_steps` tetrahedra in
/// all (counted in [`Walk::steps`], so the limit does not depend on where a
/// walk was handed over). `visit(cell, λ₀, λ₁)` is called once for every
/// cell the segment passes, in order, also with `λ₀ = λ₁` where the segment
/// only touches a cell at an edge or a corner; returning `false` stops the
/// walk ([`WalkEnd::Stopped`]).
///
/// `w` must stand in tetrahedron (cell, face, edge) of a local cell of `m`;
/// the indices are checked.
pub fn walk<F: FnMut(u32, f64, f64) -> bool>(
    m: &WalkMesh<'_>,
    w: &mut Walk,
    max_steps: u32,
    mut visit: F,
) -> WalkEnd {
    let n_interior = m.neighbour.len();
    let d = sub(w.target, w.origin);
    // where the segment entered the current cell
    let mut start = w.lambda;
    while w.steps < max_steps {
        if (w.cell as usize) < m.convex.len() && m.convex[w.cell as usize] != 0 {
            // a convex cell with planar faces (RA0): the exit is the first
            // plane the segment passes, no tetrahedra
            w.steps += 1;
            let cell = w.cell as usize;
            let cf = m.cell_ptr[cell] as usize..m.cell_ptr[cell + 1] as usize;
            let mut best = (f64::INFINITY, usize::MAX);
            for &entry in &m.cell_faces[cf] {
                let (g, neighbour_side) = unpack_face(entry);
                let pl = m.plane[g];
                let sg = if neighbour_side { -1.0 } else { 1.0 };
                let n = [pl[0] * sg, pl[1] * sg, pl[2] * sg];
                let off = pl[3] * sg;
                // inner distances (positive inside) of the target and of
                // the current point
                let st = off - dot(n, w.target);
                if st < 0.0 {
                    let sx = off - dot(n, w.x);
                    if sx > st {
                        let s0 = sx.max(0.0);
                        let l = s0 / (s0 - st);
                        if l < best.0 {
                            best = (l, g);
                        }
                    }
                }
            }
            if best.1 == usize::MAX {
                w.x = w.target;
                w.lambda = 1.0;
                settle(m, w);
                return if visit(w.cell, start, 1.0) {
                    WalkEnd::Reached
                } else {
                    WalkEnd::Stopped
                };
            }
            let g = best.1;
            let rest = 1.0 - w.lambda;
            w.lambda += best.0.clamp(0.0, 1.0) * rest;
            w.x = add(w.origin, scale(d, w.lambda));
            if g >= n_interior || g >= m.periodic_start {
                w.face = g as u32;
                w.edge = fan_edge(m, g, w.x);
                w.entered = 0;
                return if visit(w.cell, start, w.lambda) {
                    WalkEnd::Boundary { face: w.face }
                } else {
                    WalkEnd::Stopped
                };
            }
            if !visit(w.cell, start, w.lambda) {
                settle(m, w);
                return WalkEnd::Stopped;
            }
            let (p, q) = (m.owner[g], m.neighbour[g]);
            let other = if p as usize == cell { q } else { p };
            start = w.lambda;
            w.cell = other;
            w.crossed += 1;
            w.face = g as u32;
            // the edge of the fan triangle is found only when a tetrahedron
            // step needs it (or the walk ends)
            w.edge = NO_EDGE;
            w.entered = 0;
            if other as usize >= m.n_cells {
                w.edge = fan_edge(m, g, w.x);
                return WalkEnd::Processor { face: w.face };
            }
            if w.lambda >= 1.0 {
                w.x = w.target;
                w.lambda = 1.0;
                settle(m, w);
                return if visit(w.cell, start, 1.0) {
                    WalkEnd::Reached
                } else {
                    WalkEnd::Stopped
                };
            }
            continue;
        }
        w.steps += 1;
        let cell = w.cell as usize;
        let face = w.face as usize;
        if w.edge == NO_EDGE {
            w.edge = fan_edge(m, face, w.x);
        }
        let t = m.tet(cell, face, w.edge as usize);
        let sb = sides(&t, w.target);
        // exit face: smallest crossing parameter among the faces the
        // target lies behind; the distance of the current point is
        // computed only for those faces (RA0: most tetrahedra have one
        // or none), with the same expression, so the result is bitwise
        // the one that evaluated all four
        let mut best = (f64::INFINITY, NONE);
        for j in 0..4u32 {
            let b = sb[j as usize];
            if b < 0.0 && j != w.entered {
                let a = side(&t, j, w.x);
                let s = a.max(0.0);
                let l = s / (s - b);
                if l < best.0 {
                    best = (l, j);
                }
            }
        }
        if best.1 == NONE {
            w.x = w.target;
            w.lambda = 1.0;
            return if visit(w.cell, start, 1.0) {
                WalkEnd::Reached
            } else {
                WalkEnd::Stopped
            };
        }
        let l = best.0.clamp(0.0, 1.0);
        // parameter along the whole segment
        let rest = 1.0 - w.lambda;
        w.lambda += l * rest;
        w.x = add(w.origin, scale(d, w.lambda));
        let k = m.face(face).len();
        let owner = m.owner[face] as usize == cell;
        match best.1 {
            0 => {
                // through the mesh face
                if face >= n_interior || face >= m.periodic_start {
                    return if visit(w.cell, start, w.lambda) {
                        WalkEnd::Boundary { face: w.face }
                    } else {
                        WalkEnd::Stopped
                    };
                }
                if !visit(w.cell, start, w.lambda) {
                    return WalkEnd::Stopped;
                }
                let (p, q) = (m.owner[face], m.neighbour[face]);
                let other = if p as usize == cell { q } else { p };
                start = w.lambda;
                w.cell = other;
                w.crossed += 1;
                w.entered = 0;
                if other as usize >= m.n_cells {
                    return WalkEnd::Processor { face: w.face };
                }
            }
            j @ (1 | 2) => {
                // (a, p, m) or (a, m, q): neighbouring fan triangle; p is
                // node `edge` when the cell owns the face
                let towards_p = (j == 1) == owner;
                let e = w.edge as usize;
                w.edge = if towards_p {
                    (if e == 0 { k - 1 } else { e - 1 }) as u32
                } else {
                    (if e + 1 == k { 0 } else { e + 1 }) as u32
                };
                w.entered = 3 - j; // the shared face has the other index
            }
            _ => {
                // (a, q, p): the other face of the cell with edge p-q
                let nodes = m.face(face);
                let e = w.edge as usize;
                let (u, v) = (nodes[e], nodes[(e + 1) % k]);
                let mut found = None;
                let cf = m.cell_ptr[cell] as usize..m.cell_ptr[cell + 1] as usize;
                for &entry in &m.cell_faces[cf] {
                    let (g, _) = unpack_face(entry);
                    if g == face {
                        continue;
                    }
                    let gn = m.face(g);
                    let kg = gn.len();
                    for e in 0..kg {
                        let (x0, x1) = (gn[e], gn[(e + 1) % kg]);
                        if (x0 == u && x1 == v) || (x0 == v && x1 == u) {
                            found = Some((g, e));
                        }
                    }
                }
                match found {
                    Some((g, e)) => {
                        w.face = g as u32;
                        w.edge = e as u32;
                    }
                    None => return WalkEnd::Lost,
                }
                w.entered = 3;
            }
        }
        if w.lambda >= 1.0 {
            // the target lies on the face just crossed
            w.x = w.target;
            w.lambda = 1.0;
            return if visit(w.cell, start, 1.0) {
                WalkEnd::Reached
            } else {
                WalkEnd::Stopped
            };
        }
    }
    WalkEnd::Lost
}
