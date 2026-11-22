//! The polyhedral face-based mesh (D23).
//!
//! Layout, chosen so that the kernels of gate G1 run on it unchanged:
//!
//! - **Cells** `0..n_cells` are the cells of this partition: first the
//!   inner cells `0..n_inner_cells`, whose face neighbours are all local,
//!   then the cells that touch the halo. Halo cells, the copies of
//!   neighbouring cells owned by other ranks, follow as
//!   `n_cells..n_cells + n_halo`, grouped by the rank that owns them (see
//!   [`crate::decompose`]). A solver can finish the inner cells while the
//!   halo exchange is in flight (D38).
//! - **Faces** come in three ranges:
//!   1. inner faces `0..n_inner_faces`: both cells local, `owner < neighbour`;
//!   2. processor faces `n_inner_faces..n_interior_faces()`: the neighbour is
//!      a halo cell; kept apart so that their work can wait for the halo
//!      exchange while the inner faces are computed (D38);
//!   3. boundary faces `n_interior_faces()..n_faces()`, grouped by boundary.
//!
//!   `owner` has one entry per face, `neighbour` one per interior face. The
//!   nodes of a face are ordered so that its normal (right-hand rule) points
//!   out of the owner cell.
//! - **Cell to face** lists (`cell_ptr`, `cell_faces`) for the gather loops,
//!   packed with the orientation bit ([`foehncraft_kernels::pack_face`]), faces
//!   in ascending order of the global mesh (for a mesh that is not a rank's
//!   part, ascending local order; for a local mesh of
//!   [`crate::decompose`], ascending face position in the global mesh, so
//!   every per-cell sum over faces runs in global order, D93 (1)). Only for
//!   local cells.
//!
//! Inner faces are sorted by owner, then neighbour (upper-triangular order),
//! which gives the matrix rows the order of a well-numbered mesh.
//!
//! **Periodic faces** (D29, D109 (1); [`crate::periodic`]). On a single
//! partition the faces of translational periodic boundary pairs take the
//! place of the processor faces: `n_inner_faces..n_interior_faces()`. Each
//! side keeps its own face (its nodes, area and centroid stay those of the
//! boundary face), owned by its cell and listed only in that cell's face
//! list, with the cell of the partner face as neighbour; [`Periodic`]
//! holds the shift that moves the neighbour's centroid next to the face.
//! So every flux of a pair is computed from both sides with the same
//! formula, the matrix rows couple the two cells, and the geometry of the
//! cells stays that of the mesh as read.
//!
//! **Global identities.** Every cell, face and node carries a 64-bit
//! identity, its position in the input mesh (as read or generated, before
//! any renumbering or decomposition). Output and restart files are written
//! in this order, so they do not depend on the cell order or on the number
//! of ranks (D41, D42).

use foehncraft_core::index::{GlobalIndex, MAX_LOCAL_INDEX};
use foehncraft_kernels::{pack_face, unpack_face};

/// A named set of boundary faces.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Boundary {
    /// Name, from the file (Gmsh physical group, SU2 marker) or the generator.
    pub name: String,
    /// First face of the boundary.
    pub start: usize,
    /// Number of faces.
    pub len: usize,
}

impl Boundary {
    /// Face range of the boundary.
    pub fn range(&self) -> std::ops::Range<usize> {
        self.start..self.start + self.len
    }
}

/// One translational periodic pair: the faces of boundary `a` meet those
/// of boundary `b` when moved by `translation` (x_b = x_a + T).
#[derive(Debug, Clone, PartialEq)]
pub struct PeriodicPair {
    /// Name of the first boundary.
    pub a: String,
    /// Name of the partner boundary.
    pub b: String,
    /// Translation T from `a` to `b` (m).
    pub translation: [f64; 3],
    /// The faces of side `a` (consecutive periodic faces).
    pub a_faces: std::ops::Range<usize>,
    /// The faces of side `b`, in the order of their partners on side `a`.
    pub b_faces: std::ops::Range<usize>,
}

/// The periodic faces of a mesh (see the module docs).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Periodic {
    /// The pairs.
    pub pairs: Vec<PeriodicPair>,
    /// First periodic face (`n_inner_faces`); the periodic faces end at
    /// `n_interior_faces()`.
    pub start: usize,
    /// Per periodic face (index `f − start`): the shift s that places the
    /// neighbour cell next to the face, x_N + s (−T on side `a`, +T on
    /// side `b`).
    pub shift: Vec<[f64; 3]>,
    /// Per periodic face: the partner face (absolute index).
    pub partner: Vec<u32>,
    /// Per periodic face: its pair.
    pub pair: Vec<u32>,
}

impl Periodic {
    /// Number of periodic faces.
    pub fn len(&self) -> usize {
        self.shift.len()
    }

    /// True without periodic faces.
    pub fn is_empty(&self) -> bool {
        self.shift.is_empty()
    }

    /// −1 for a face on side `a` of its pair, +1 on side `b`: the sign of
    /// its shift along the translation, so that a variable with the trend
    /// φ(x + T) = φ(x) + Δ seen across face `f` is φ_N + side(f) Δ.
    pub fn side(&self, f: usize) -> f64 {
        let i = f - self.start;
        let pr = &self.pairs[self.pair[i] as usize];
        if pr.a_faces.contains(&f) { -1.0 } else { 1.0 }
    }

    /// The pair of boundary `name`, if it is one of the sides.
    pub fn pair_of(&self, name: &str) -> Option<&PeriodicPair> {
        self.pairs.iter().find(|p| p.a == name || p.b == name)
    }
}

/// A polyhedral mesh in face-based form.
#[derive(Debug, Clone, Default)]
pub struct Mesh {
    /// Node coordinates.
    pub points: Vec<[f64; 3]>,
    /// Start of the node list of each face in `face_nodes` (length
    /// `n_faces + 1`).
    pub face_ptr: Vec<u32>,
    /// Nodes of each face, ordered so that the normal points out of the
    /// owner.
    pub face_nodes: Vec<u32>,
    /// Owner cell of each face.
    pub owner: Vec<u32>,
    /// Neighbour cell of each interior face (inner and processor faces).
    pub neighbour: Vec<u32>,
    /// Number of inner faces (both cells local).
    pub n_inner_faces: usize,
    /// Number of local cells.
    pub n_cells: usize,
    /// Number of inner cells (local cells with no halo neighbour); equal to
    /// `n_cells` without a decomposition.
    pub n_inner_cells: usize,
    /// Number of halo cells.
    pub n_halo: usize,
    /// Boundary boundaries, in face order.
    pub boundaries: Vec<Boundary>,
    /// Start of the face list of each local cell (length `n_cells + 1`).
    pub cell_ptr: Vec<u32>,
    /// Faces of each local cell, packed with their orientation.
    pub cell_faces: Vec<u32>,
    /// Cell zone of each local cell (index into `zones`).
    pub cell_zone: Vec<u32>,
    /// Names of the cell zones (at least one).
    pub zones: Vec<String>,
    /// Global identity of each local and halo cell (input order).
    pub cell_id: Vec<GlobalIndex>,
    /// Global identity of each face (input order).
    pub face_id: Vec<GlobalIndex>,
    /// Global identity of each node (input order).
    pub point_id: Vec<GlobalIndex>,
    /// Translational periodic faces (`None`: none; see the module docs).
    pub periodic: Option<Periodic>,
    /// Faces stored with the reverse orientation of the global mesh they
    /// were extracted from (empty: none, the mesh is its own global mesh).
    /// A flipped face has its owner and neighbour swapped and its nodes
    /// reversed after the first one; geometry evaluates it in the global
    /// node order and flips the sign, so every rank gets the bits of the
    /// one-rank run (D93 (3); see [`crate::decompose`]).
    pub face_flip: Vec<bool>,
}

impl Mesh {
    /// Number of faces.
    pub fn n_faces(&self) -> usize {
        self.owner.len()
    }

    /// Number of interior faces (inner and processor).
    pub fn n_interior_faces(&self) -> usize {
        self.neighbour.len()
    }

    /// Number of boundary faces.
    pub fn n_boundary_faces(&self) -> usize {
        self.n_faces() - self.n_interior_faces()
    }

    /// Number of nodes.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }

    /// Nodes of face `f`.
    #[inline]
    pub fn face(&self, f: usize) -> &[u32] {
        &self.face_nodes[self.face_ptr[f] as usize..self.face_ptr[f + 1] as usize]
    }

    /// Whether face `f` is stored reversed against the global mesh (see
    /// [`Mesh::face_flip`]).
    #[inline]
    pub fn is_flipped(&self, f: usize) -> bool {
        self.face_flip.get(f).copied().unwrap_or(false)
    }

    /// Nodes of face `f` in the order of the global mesh, and whether the
    /// face is flipped. Reversing the loop after its first node is its own
    /// inverse, so the global order is recovered exactly.
    pub fn face_global_order(&self, f: usize) -> (std::borrow::Cow<'_, [u32]>, bool) {
        let n = self.face(f);
        if self.is_flipped(f) {
            let mut v = Vec::with_capacity(n.len());
            v.push(n[0]);
            v.extend(n[1..].iter().rev());
            (std::borrow::Cow::Owned(v), true)
        } else {
            (std::borrow::Cow::Borrowed(n), false)
        }
    }

    /// Packed faces of local cell `c` (see [`foehncraft_kernels::unpack_face`]).
    #[inline]
    pub fn cell(&self, c: usize) -> &[u32] {
        &self.cell_faces[self.cell_ptr[c] as usize..self.cell_ptr[c + 1] as usize]
    }

    /// Faces of cell `c` with the flag `true` where the cell is the
    /// neighbour (the stored normal points into the cell).
    pub fn cell_face_iter(&self, c: usize) -> impl Iterator<Item = (usize, bool)> + '_ {
        self.cell(c).iter().map(|&e| unpack_face(e))
    }

    /// First periodic face: `n_interior_faces()` without periodic faces.
    #[inline]
    pub fn periodic_start(&self) -> usize {
        self.periodic
            .as_ref()
            .map_or(self.n_interior_faces(), |p| p.start)
    }

    /// True for a periodic face.
    #[inline]
    pub fn is_periodic(&self, f: usize) -> bool {
        f >= self.periodic_start() && f < self.n_interior_faces()
    }

    /// The cell on the other side of face `f` seen from cell `c`, or `None`
    /// for a boundary face and a periodic face (geometrically a boundary:
    /// its neighbour lies at the other end of the domain).
    #[inline]
    pub fn across(&self, f: usize, c: usize) -> Option<usize> {
        if f >= self.n_interior_faces() || self.is_periodic(f) {
            return None;
        }
        let (p, n) = (self.owner[f] as usize, self.neighbour[f] as usize);
        Some(if p == c { n } else { p })
    }

    /// Boundary index of boundary face `f`.
    pub fn boundary_of(&self, f: usize) -> Option<usize> {
        self.boundaries.iter().position(|p| p.range().contains(&f))
    }

    /// Boundary by name.
    pub fn boundary(&self, name: &str) -> Option<&Boundary> {
        self.boundaries.iter().find(|p| p.name == name)
    }

    /// Distinct nodes of local cell `c`, in first-seen order.
    pub fn cell_nodes(&self, c: usize) -> Vec<u32> {
        let mut v: Vec<u32> = Vec::with_capacity(16);
        for (f, _) in self.cell_face_iter(c) {
            for &n in self.face(f) {
                if !v.contains(&n) {
                    v.push(n);
                }
            }
        }
        v
    }

    /// Shape of local cell `c`, from the number of faces and their sizes.
    pub fn cell_shape(&self, c: usize) -> CellShape {
        let mut tri = 0;
        let mut quad = 0;
        let mut other = 0;
        for (f, _) in self.cell_face_iter(c) {
            match self.face(f).len() {
                3 => tri += 1,
                4 => quad += 1,
                _ => other += 1,
            }
        }
        match (tri, quad, other) {
            (4, 0, 0) => CellShape::Tetrahedron,
            (0, 6, 0) => CellShape::Hexahedron,
            (2, 3, 0) => CellShape::Prism,
            (4, 1, 0) => CellShape::Pyramid,
            _ => CellShape::Polyhedron,
        }
    }

    /// Checks the connectivity: index ranges, face ranges, orientation bits,
    /// that every face appears in the lists of its cells, that every local
    /// cell is closed in the topological sense (every edge shared by exactly
    /// two of its faces, in opposite directions) and the precondition of the
    /// unchecked loads of `foehncraft-kernels`.
    pub fn validate(&self) -> Result<(), String> {
        let nf = self.n_faces();
        let ni = self.n_interior_faces();
        let nc = self.n_cells;
        let nt = nc + self.n_halo;
        if nt > MAX_LOCAL_INDEX as usize || self.face_nodes.len() > MAX_LOCAL_INDEX as usize {
            return Err("mesh exceeds the 32-bit local index range".into());
        }
        if self.face_ptr.len() != nf + 1 || self.face_ptr[0] != 0 {
            return Err("face_ptr has the wrong length or start".into());
        }
        if *self.face_ptr.last().unwrap() as usize != self.face_nodes.len() {
            return Err("face_ptr does not end at face_nodes.len()".into());
        }
        if self.n_inner_faces > ni {
            return Err("more inner faces than interior faces".into());
        }
        for f in 0..nf {
            let (a, b) = (self.face_ptr[f], self.face_ptr[f + 1]);
            if b < a + 3 {
                return Err(format!("face {f} has fewer than 3 nodes"));
            }
            if self
                .face(f)
                .iter()
                .any(|&n| n as usize >= self.points.len())
            {
                return Err(format!("face {f}: node out of range"));
            }
            let p = self.owner[f] as usize;
            if p >= nc {
                return Err(format!("face {f}: owner {p} is not a local cell"));
            }
            if f < ni {
                let n = self.neighbour[f] as usize;
                if f < self.n_inner_faces && (n >= nc || n <= p) {
                    return Err(format!("inner face {f}: cells {p}, {n} not ordered"));
                }
                if f >= self.n_inner_faces && f < self.periodic_start() && (n < nc || n >= nt) {
                    return Err(format!(
                        "processor face {f}: neighbour {n} is not a halo cell"
                    ));
                }
            }
        }
        let mut next = ni;
        for p in &self.boundaries {
            if p.start != next {
                return Err(format!(
                    "boundary {} does not follow the previous one",
                    p.name
                ));
            }
            next += p.len;
        }
        if next != nf {
            return Err("boundaries do not cover all boundary faces".into());
        }
        if self.cell_ptr.len() != nc + 1 || self.cell_ptr[0] != 0 {
            return Err("cell_ptr has the wrong length or start".into());
        }
        if self.cell_ptr.windows(2).any(|w| w[0] > w[1])
            || *self.cell_ptr.last().unwrap() as usize != self.cell_faces.len()
        {
            return Err("cell_ptr is not a valid offset array".into());
        }
        let expected_entries = nf + self.n_inner_faces;
        if self.cell_faces.len() != expected_entries {
            return Err(format!(
                "cell_faces has {} entries, expected {expected_entries}",
                self.cell_faces.len()
            ));
        }
        if self.cell_zone.len() != nc || self.zones.is_empty() {
            return Err("cell zones missing".into());
        }
        if self
            .cell_zone
            .iter()
            .any(|&z| z as usize >= self.zones.len())
        {
            return Err("cell zone index out of range".into());
        }
        if self.cell_id.len() != nt
            || self.face_id.len() != nf
            || self.point_id.len() != self.points.len()
        {
            return Err("global identities missing".into());
        }
        if self.n_inner_cells > nc {
            return Err("more inner cells than local cells".into());
        }
        if let Some(pd) = &self.periodic {
            if pd.start != self.n_inner_faces || self.n_halo != 0 {
                return Err("periodic faces need a single partition after the inner faces".into());
            }
            let np = ni - pd.start;
            if pd.shift.len() != np || pd.partner.len() != np || pd.pair.len() != np {
                return Err("periodic arrays have the wrong length".into());
            }
            for f in pd.start..ni {
                let i = f - pd.start;
                let q = pd.partner[i] as usize;
                if q < pd.start || q >= ni || pd.partner[q - pd.start] as usize != f {
                    return Err(format!("periodic face {f}: partner {q} is not paired back"));
                }
                if self.neighbour[f] != self.owner[q] || self.owner[f] == self.neighbour[f] {
                    return Err(format!("periodic face {f}: wrong neighbour"));
                }
            }
        }
        for f in self.n_inner_faces..self.periodic_start() {
            if (self.owner[f] as usize) < self.n_inner_cells {
                return Err(format!(
                    "processor face {f} belongs to inner cell {}",
                    self.owner[f]
                ));
            }
        }
        if !self.face_flip.is_empty() && self.face_flip.len() != nf {
            return Err(format!(
                "face_flip has {} entries for {nf} faces",
                self.face_flip.len()
            ));
        }
        // a rank's part lists the faces of a cell in global face order (the
        // order of the mesh it was extracted from), not in local order
        let local_order = self.face_flip.is_empty();
        let mut edges: Vec<(u32, u32)> = Vec::with_capacity(64);
        let mut seen: Vec<usize> = Vec::with_capacity(64);
        for c in 0..nc {
            edges.clear();
            seen.clear();
            let mut last = None;
            for (f, is_nb) in self.cell_face_iter(c) {
                if f >= nf {
                    return Err(format!("cell {c}: face {f} out of range"));
                }
                if local_order && last.is_some_and(|l| l >= f) {
                    return Err(format!("cell {c}: faces not in ascending order"));
                }
                last = Some(f);
                seen.push(f);
                let expected = if is_nb {
                    if f >= ni {
                        return Err(format!("cell {c}: boundary face {f} flagged as neighbour"));
                    }
                    self.neighbour[f]
                } else {
                    self.owner[f]
                };
                if expected as usize != c {
                    return Err(format!("cell {c}: face {f} has the wrong orientation"));
                }
                let nodes = self.face(f);
                let k = nodes.len();
                for i in 0..k {
                    let (a, b) = (nodes[i], nodes[(i + 1) % k]);
                    // outward direction of the edge seen from this cell
                    edges.push(if is_nb { (b, a) } else { (a, b) });
                }
            }
            seen.sort_unstable();
            if seen.windows(2).any(|w| w[0] == w[1]) {
                return Err(format!("cell {c}: a face is listed twice"));
            }
            edges.sort_unstable();
            for w in edges.windows(2) {
                if w[0] == w[1] {
                    return Err(format!(
                        "cell {c}: edge {:?} used twice in one direction",
                        w[0]
                    ));
                }
            }
            for &(a, b) in &edges {
                if edges.binary_search(&(b, a)).is_err() {
                    return Err(format!(
                        "cell {c} is not closed: edge ({a}, {b}) has no partner"
                    ));
                }
            }
        }
        Ok(())
    }

    /// Builds the cell-to-face lists from `owner` and `neighbour`.
    pub(crate) fn build_cell_faces(&mut self) {
        let nc = self.n_cells;
        let mut count = vec![0u32; nc + 1];
        for &p in &self.owner {
            count[p as usize + 1] += 1;
        }
        for &n in &self.neighbour[..self.n_inner_faces] {
            count[n as usize + 1] += 1;
        }
        for c in 0..nc {
            count[c + 1] += count[c];
        }
        let mut next = count.clone();
        let mut cell_faces = vec![0u32; count[nc] as usize];
        // Faces are visited in ascending order, so each list is sorted.
        for f in 0..self.n_faces() {
            let p = self.owner[f] as usize;
            cell_faces[next[p] as usize] = pack_face(f as u32, false);
            next[p] += 1;
            if f < self.n_inner_faces {
                let n = self.neighbour[f] as usize;
                cell_faces[next[n] as usize] = pack_face(f as u32, true);
                next[n] += 1;
            }
        }
        self.cell_ptr = count;
        self.cell_faces = cell_faces;
    }

    /// Bounding box of the nodes.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &self.points {
            for d in 0..3 {
                lo[d] = lo[d].min(p[d]);
                hi[d] = hi[d].max(p[d]);
            }
        }
        (lo, hi)
    }

    /// Renames a boundary.
    pub fn rename_boundary(&mut self, old: &str, new: &str) -> Result<(), String> {
        if self.boundary(new).is_some() {
            return Err(format!("boundary {new} exists already"));
        }
        let p = self
            .boundaries
            .iter_mut()
            .find(|p| p.name == old)
            .ok_or_else(|| format!("no boundary named {old}"))?;
        p.name = new.to_string();
        Ok(())
    }
}

/// Shape of a cell, recognised from its faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum CellShape {
    /// 4 triangles.
    Tetrahedron,
    /// 6 quadrilaterals.
    Hexahedron,
    /// 2 triangles and 3 quadrilaterals.
    Prism,
    /// 4 triangles and 1 quadrilateral.
    Pyramid,
    /// Anything else.
    Polyhedron,
}

impl CellShape {
    /// All shapes, in report order.
    pub const ALL: [CellShape; 5] = [
        Self::Tetrahedron,
        Self::Hexahedron,
        Self::Prism,
        Self::Pyramid,
        Self::Polyhedron,
    ];

    /// Name for reports.
    pub fn name(self) -> &'static str {
        match self {
            Self::Tetrahedron => "tetrahedra",
            Self::Hexahedron => "hexahedra",
            Self::Prism => "prisms",
            Self::Pyramid => "pyramids",
            Self::Polyhedron => "polyhedra",
        }
    }
}
