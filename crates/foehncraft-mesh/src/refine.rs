//! Refinement of hexahedra into hexahedra with hanging-node faces (L13
//! stage 2f, D57): the first piece of the refinement kernel of L14.
//!
//! Every hexahedron of the input mesh is the root of a tree of hexahedra
//! halved in each of its three reference directions (in a planar mesh, one
//! cell thick, the direction across the layer stays whole). A cell is
//! refined where the caller's marking function asks for it, up to
//! [`RefineOptions::levels`] times; then the tree is balanced, 2:1 across
//! faces, edges and corners, so that a face of a cell meets at most four
//! (planar mesh: two) faces of its neighbour and an edge holds at most one
//! hanging node. Cells that are not hexahedra are kept as they are.
//!
//! **Nodes.** A node is named by the lowest entity of the input mesh it lies
//! on: an input node, a point on an edge of a hexahedron (its two input
//! nodes and the position from the lower one), a point on a face (the four
//! nodes in an order fixed by their numbers and the position in that frame)
//! or a point inside a hexahedron. Its position comes from that entity
//! alone (linear on the edge, bilinear on the face, trilinear inside), so
//! the two cells on either side of a face compute the same point bit for
//! bit, and the names give the order of the new nodes: the result does not
//! depend on threads or ranks (D93).
//!
//! **Faces.** The faces of the input mesh and the planes inside a refined
//! hexahedron are tiled by the cell faces of each side; two tilings of
//! squares halved in steps are nested, so every face of the result is the
//! smaller of two overlapping squares, between the two cells. A face is a
//! polygon with every node of the mesh on its edges (the hanging nodes of
//! D23), so every cell stays closed edge by edge. The faces of a cell face
//! of the input mesh add up to its area vector exactly (the hanging nodes
//! lie on the straight edges), and Σ S = 0 holds on every cell.

use crate::builder::MeshBuilder;
use crate::geometry::{cross, dot, norm, sub};
use crate::mesh::Mesh;
use std::collections::{BTreeSet, HashMap, HashSet};

const NONE: u32 = u32::MAX;

/// The most levels of [`refine`] (a cell is halved at most this often).
pub const MAX_LEVELS: u32 = 6;

/// Settings of [`refine`].
#[derive(Debug, Clone, Copy, Default)]
pub struct RefineOptions {
    /// Halvings of a cell at most (1 to [`MAX_LEVELS`]).
    pub levels: u32,
    /// A planar mesh (one cell thick): the reference direction of each
    /// hexahedron along this vector is not split.
    pub keep_direction: Option<[f64; 3]>,
}

/// A cell of the tree as the marking function sees it.
#[derive(Debug, Clone, Copy)]
pub struct Candidate {
    /// The input cell.
    pub parent: usize,
    /// Halvings so far (0: the input cell).
    pub level: u32,
    /// Centre (the trilinear map at the middle of the cell).
    pub centre: [f64; 3],
    /// Largest width across the split directions.
    pub size: f64,
    /// Half the longest diagonal: every point of the cell lies within it of
    /// the centre.
    pub radius: f64,
}

/// What [`refine`] did.
#[derive(Debug, Clone, Default)]
pub struct RefineReport {
    /// Cells of the input mesh.
    pub cells_in: usize,
    /// Cells of the result.
    pub cells_out: usize,
    /// Input cells that were refined.
    pub refined: usize,
    /// Input cells that are not hexahedra (kept whole).
    pub not_hexahedra: usize,
    /// Cells of the result per level (0: input cells kept whole).
    pub per_level: Vec<usize>,
    /// Cells refined by the 2:1 balance, beyond the marking.
    pub balance: usize,
    /// Interior faces between cells of different levels.
    pub hanging_faces: usize,
    /// Faces with more nodes than their corners (hanging nodes on an edge).
    pub faces_with_hanging_nodes: usize,
}

/// The refined mesh.
#[derive(Debug, Clone)]
pub struct Refined {
    /// The mesh: the input nodes first and in their order, then the new nodes.
    pub mesh: Mesh,
    /// Input cell of each cell.
    pub parent: Vec<u32>,
    /// Level of each cell.
    pub level: Vec<u8>,
    /// Statistics.
    pub report: RefineReport,
}

/// The name of a node (see the module docs); positions in units of
/// 1/2^levels of the input entity.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
enum Key {
    Corner(u32),
    Edge(u32, u32, u32),
    Face([u32; 4], u32, u32),
    Cell(u32, [u32; 3]),
}

#[derive(Clone, Copy, Debug)]
struct Leaf {
    lo: [u32; 3],
    n: [u32; 3],
    level: u8,
}

/// The eight nodes of a hexahedron in lexicographic order (index
/// i + 2 j + 4 k for the corner (i, j, k) of the unit cube), or `None`.
fn hex_corners(mesh: &Mesh, c: usize) -> Option<[u32; 8]> {
    let faces: Vec<&[u32]> = mesh.cell_face_iter(c).map(|(f, _)| mesh.face(f)).collect();
    if faces.len() != 6 || faces.iter().any(|f| f.len() != 4) {
        return None;
    }
    let b = [faces[0][0], faces[0][1], faces[0][2], faces[0][3]];
    let pos = |n: u32| b.iter().position(|&x| x == n);
    let mut t = [NONE; 4];
    for f in &faces[1..] {
        for i in 0..4 {
            let (p, q) = (f[i], f[(i + 1) % 4]);
            for (x, y) in [(p, q), (q, p)] {
                if let (Some(k), None) = (pos(x), pos(y)) {
                    if t[k] != NONE && t[k] != y {
                        return None;
                    }
                    t[k] = y;
                }
            }
        }
    }
    let mut all: Vec<u32> = b.iter().chain(t.iter()).copied().collect();
    if t.contains(&NONE) {
        return None;
    }
    all.sort_unstable();
    all.dedup();
    if all.len() != 8 {
        return None;
    }
    // the top face and the four sides must be faces of the cell
    let has = |q: [u32; 4]| {
        let mut q = q;
        q.sort_unstable();
        faces.iter().any(|f| {
            let mut g = [f[0], f[1], f[2], f[3]];
            g.sort_unstable();
            g == q
        })
    };
    if !has(t) {
        return None;
    }
    for i in 0..4 {
        let j = (i + 1) % 4;
        if !has([b[i], b[j], t[j], t[i]]) {
            return None;
        }
    }
    Some([b[0], b[1], b[3], b[2], t[0], t[1], t[3], t[2]])
}

fn corner_index(p: [u32; 3], r: u32) -> usize {
    (p[0] == r) as usize + 2 * (p[1] == r) as usize + 4 * (p[2] == r) as usize
}

fn lerp(a: [f64; 3], b: [f64; 3], s: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * s,
        a[1] + (b[1] - a[1]) * s,
        a[2] + (b[2] - a[2]) * s,
    ]
}

fn trilinear(pts: &[[f64; 3]], h: &[u32; 8], s: [f64; 3]) -> [f64; 3] {
    let p = |i: usize| pts[h[i] as usize];
    let x00 = lerp(p(0), p(1), s[0]);
    let x10 = lerp(p(2), p(3), s[0]);
    let x01 = lerp(p(4), p(5), s[0]);
    let x11 = lerp(p(6), p(7), s[0]);
    lerp(lerp(x00, x10, s[1]), lerp(x01, x11, s[1]), s[2])
}

/// A face frame: the four nodes [origin, along u, along v, opposite].
type Frame = [u32; 4];

/// The frame of a quadrilateral given as `q[a][b]` (the corner at (a, b) of
/// two free directions), and the map of a point (a, b) into it.
fn frame_of(q: [[u32; 2]; 2]) -> (Frame, impl Fn(u32, u32, u32) -> (u32, u32)) {
    let mut best = (0usize, 0usize);
    for a in 0..2 {
        for b in 0..2 {
            if q[a][b] < q[best.0][best.1] {
                best = (a, b);
            }
        }
    }
    let (oa, ob) = best;
    let along_a = q[1 - oa][ob];
    let along_b = q[oa][1 - ob];
    let swap = along_b < along_a;
    let frame = if swap {
        [q[oa][ob], along_b, along_a, q[1 - oa][1 - ob]]
    } else {
        [q[oa][ob], along_a, along_b, q[1 - oa][1 - ob]]
    };
    let map = move |a: u32, b: u32, r: u32| {
        let da = if oa == 0 { a } else { r - a };
        let db = if ob == 0 { b } else { r - b };
        if swap { (db, da) } else { (da, db) }
    };
    (frame, map)
}

fn edge_key(n0: u32, n1: u32, t: u32, r: u32) -> Key {
    if t == 0 {
        Key::Corner(n0)
    } else if t == r {
        Key::Corner(n1)
    } else if n0 < n1 {
        Key::Edge(n0, n1, t)
    } else {
        Key::Edge(n1, n0, r - t)
    }
}

/// The name of the point (u, v) of a frame.
fn frame_key(f: &Frame, u: u32, v: u32, r: u32) -> Key {
    let [o, nu, nv, w] = *f;
    match (u == 0 || u == r, v == 0 || v == r) {
        (true, _) => {
            let (a, b) = if u == 0 { (o, nv) } else { (nu, w) };
            edge_key(a, b, v, r)
        }
        (false, true) => {
            let (a, b) = if v == 0 { (o, nu) } else { (nv, w) };
            edge_key(a, b, u, r)
        }
        _ => Key::Face(*f, u, v),
    }
}

struct Tree<'a> {
    mesh: &'a Mesh,
    r: u32,
    hex: Vec<Option<[u32; 8]>>,
    split: Vec<[bool; 3]>,
    leaves: Vec<Vec<Leaf>>,
}

impl Tree<'_> {
    /// The name of the point p (units of 1/r) of hexahedron c.
    fn cell_key(&self, c: usize, p: [u32; 3]) -> Key {
        let h = self.hex[c].as_ref().expect("a hexahedron");
        let r = self.r;
        let at: Vec<bool> = p.iter().map(|&x| x == 0 || x == r).collect();
        let nb = at.iter().filter(|&&a| a).count();
        match nb {
            3 => Key::Corner(h[corner_index(p, r)]),
            2 => {
                let f = (0..3).find(|&d| !at[d]).unwrap();
                let mut p0 = p;
                p0[f] = 0;
                let mut p1 = p;
                p1[f] = r;
                edge_key(h[corner_index(p0, r)], h[corner_index(p1, r)], p[f], r)
            }
            1 => {
                let f = (0..3).find(|&d| at[d]).unwrap();
                let (g, k) = match f {
                    0 => (1, 2),
                    1 => (0, 2),
                    _ => (0, 1),
                };
                let (fr, map) = self.local_frame(c, f, p[f]);
                let _ = h;
                let (u, v) = map(p[g], p[k], r);
                frame_key(&fr, u, v, r)
            }
            _ => Key::Cell(c as u32, p),
        }
    }

    /// The frame of the face of hexahedron c at p[f] = side, and the map
    /// of its free coordinates (increasing directions) into the frame.
    fn local_frame(
        &self,
        c: usize,
        f: usize,
        side: u32,
    ) -> (Frame, impl Fn(u32, u32, u32) -> (u32, u32)) {
        let h = self.hex[c].as_ref().expect("a hexahedron");
        let r = self.r;
        let (g, k) = match f {
            0 => (1, 2),
            1 => (0, 2),
            _ => (0, 1),
        };
        let q: [[u32; 2]; 2] = std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                let mut p = [0u32; 3];
                p[f] = side;
                p[g] = a as u32 * r;
                p[k] = b as u32 * r;
                h[corner_index(p, r)]
            })
        });
        frame_of(q)
    }

    fn corners(&self, c: usize, l: &Leaf) -> [Key; 8] {
        let mut out = [Key::Corner(0); 8];
        for (i, o) in out.iter_mut().enumerate() {
            let m = [(i & 1) as u32, ((i >> 1) & 1) as u32, ((i >> 2) & 1) as u32];
            *o = self.cell_key(
                c,
                [
                    l.lo[0] + m[0] * l.n[0],
                    l.lo[1] + m[1] * l.n[1],
                    l.lo[2] + m[2] * l.n[2],
                ],
            );
        }
        out
    }

    fn children(&self, c: usize, l: &Leaf) -> Vec<Leaf> {
        let s = self.split[c];
        let half: [u32; 3] = std::array::from_fn(|d| if s[d] { l.n[d] / 2 } else { l.n[d] });
        let mut out = Vec::with_capacity(8);
        for k in 0..=(s[2] as u32) {
            for j in 0..=(s[1] as u32) {
                for i in 0..=(s[0] as u32) {
                    let m = [i, j, k];
                    out.push(Leaf {
                        lo: std::array::from_fn(|d| l.lo[d] + m[d] * half[d]),
                        n: half,
                        level: l.level + 1,
                    });
                }
            }
        }
        out
    }

    fn candidate(&self, c: usize, l: &Leaf) -> Candidate {
        let h = self.hex[c].as_ref().unwrap();
        let pts = &self.mesh.points;
        let r = self.r as f64;
        let at = |m: [f64; 3]| {
            trilinear(
                pts,
                h,
                std::array::from_fn(|d| (l.lo[d] as f64 + m[d] * l.n[d] as f64) / r),
            )
        };
        let centre = at([0.5; 3]);
        let mut size = 0.0f64;
        for d in 0..3 {
            if self.split[c][d] {
                let mut a = [0.5; 3];
                let mut b = [0.5; 3];
                a[d] = 0.0;
                b[d] = 1.0;
                size = size.max(norm(sub(at(b), at(a))));
            }
        }
        let mut radius = 0.0f64;
        for i in 0..4 {
            let m = [(i & 1) as f64, ((i >> 1) & 1) as f64, 0.0];
            let a = at(m);
            let b = at([1.0 - m[0], 1.0 - m[1], 1.0]);
            radius = radius.max(0.5 * norm(sub(b, a)));
        }
        Candidate {
            parent: c,
            level: l.level as u32,
            centre,
            size,
            radius,
        }
    }
}

fn register(maxlev: &mut HashMap<Key, u8>, keys: &[Key], level: u8) {
    for k in keys {
        let e = maxlev.entry(*k).or_insert(level);
        *e = (*e).max(level);
    }
}

/// Refines the hexahedra of `mesh` where `mark` asks for it (see the
/// module docs). `mark` sees each cell of the tree once per level, from
/// the input cells down; a cell it marks is halved.
pub fn refine(
    mesh: &Mesh,
    opt: &RefineOptions,
    mark: impl Fn(&Candidate) -> bool,
) -> Result<Refined, String> {
    if opt.levels == 0 || opt.levels > MAX_LEVELS {
        return Err(format!(
            "levels must be 1 to {MAX_LEVELS}, not {}",
            opt.levels
        ));
    }
    if mesh.periodic.is_some() || mesh.n_halo > 0 || !mesh.face_flip.is_empty() {
        return Err(
            "refinement needs the whole mesh without periodic pairs (refine before pairing)".into(),
        );
    }
    let r = 1u32 << opt.levels;
    let nc = mesh.n_cells;
    let pts = &mesh.points;
    let hex: Vec<Option<[u32; 8]>> = (0..nc).map(|c| hex_corners(mesh, c)).collect();
    let mut split = vec![[true; 3]; nc];
    if let Some(k) = opt.keep_direction {
        let kn = norm(k);
        for c in 0..nc {
            let Some(h) = &hex[c] else { continue };
            let mut found = 0;
            for d in 0..3 {
                // the mean edge along direction d
                let mut e = [0.0; 3];
                for i in 0..8 {
                    if i >> d & 1 == 0 {
                        let a = pts[h[i] as usize];
                        let b = pts[h[i | 1 << d] as usize];
                        e = [e[0] + b[0] - a[0], e[1] + b[1] - a[1], e[2] + b[2] - a[2]];
                    }
                }
                if dot(e, k).abs() > 0.99 * norm(e) * kn {
                    split[c][d] = false;
                    found += 1;
                }
            }
            if found != 1 {
                return Err(format!(
                    "cell {c}: a hexahedron of a planar mesh needs one direction across the layer"
                ));
            }
        }
    }
    let root = Leaf {
        lo: [0; 3],
        n: [r; 3],
        level: 0,
    };
    let mut tree = Tree {
        mesh,
        r,
        hex,
        split,
        leaves: (0..nc).map(|_| vec![root]).collect(),
    };
    let mut maxlev: HashMap<Key, u8> = HashMap::new();
    for c in 0..nc {
        if tree.hex[c].is_some() {
            let k = tree.corners(c, &root);
            register(&mut maxlev, &k, 0);
        }
    }
    // marking, level by level
    for level in 0..opt.levels as u8 {
        for c in 0..nc {
            if tree.hex[c].is_none() {
                continue;
            }
            let old = std::mem::take(&mut tree.leaves[c]);
            let mut new = Vec::with_capacity(old.len());
            for l in old {
                if l.level == level && mark(&tree.candidate(c, &l)) {
                    for ch in tree.children(c, &l) {
                        let k = tree.corners(c, &ch);
                        register(&mut maxlev, &k, ch.level);
                        new.push(ch);
                    }
                } else {
                    new.push(l);
                }
            }
            tree.leaves[c] = new;
        }
    }
    // 2:1 balance across faces, edges and corners: a cell is halved while a
    // point of its boundary on the lattice of the next level is a node of a
    // cell two or more levels finer
    let mut balance = 0usize;
    loop {
        let mut changed = false;
        for c in 0..nc {
            if tree.hex[c].is_none() {
                continue;
            }
            let s = tree.split[c];
            let mut i = 0;
            while i < tree.leaves[c].len() {
                let l = tree.leaves[c][i];
                let needs = (l.level as u32) + 1 < opt.levels && {
                    let mut hit = false;
                    'scan: for mk in 0..3u32 {
                        for mj in 0..3u32 {
                            for mi in 0..3u32 {
                                let m = [mi, mj, mk];
                                if (0..3).any(|d| !s[d] && m[d] == 1) {
                                    continue;
                                }
                                if m == [1, 1, 1] {
                                    continue;
                                }
                                let p: [u32; 3] = std::array::from_fn(|d| {
                                    if s[d] {
                                        l.lo[d] + m[d] * (l.n[d] / 2)
                                    } else {
                                        l.lo[d] + m[d] / 2 * l.n[d]
                                    }
                                });
                                let k = tree.cell_key(c, p);
                                if maxlev.get(&k).is_some_and(|&v| v >= l.level + 2) {
                                    hit = true;
                                    break 'scan;
                                }
                            }
                        }
                    }
                    hit
                };
                if needs {
                    let ch = tree.children(c, &l);
                    for x in &ch {
                        let k = tree.corners(c, x);
                        register(&mut maxlev, &k, x.level);
                    }
                    tree.leaves[c].swap_remove(i);
                    tree.leaves[c].extend(ch);
                    balance += 1;
                    changed = true;
                } else {
                    i += 1;
                }
            }
        }
        if !changed {
            break;
        }
    }
    build(&mut tree, &maxlev, balance)
}

fn build(tree: &mut Tree, maxlev: &HashMap<Key, u8>, balance: usize) -> Result<Refined, String> {
    let mesh = tree.mesh;
    let nc = mesh.n_cells;
    let r = tree.r;
    let pts = &mesh.points;
    // cells: per input cell, its leaves in (k, j, i) order
    for leaves in &mut tree.leaves {
        leaves.sort_by_key(|l| (l.lo[2], l.lo[1], l.lo[0]));
    }
    let mut first = Vec::with_capacity(nc + 1);
    first.push(0u32);
    for c in 0..nc {
        first.push(first[c] + tree.leaves[c].len() as u32);
    }
    let n_out = first[nc] as usize;
    // nodes: the input nodes, then the new ones in the order of their names
    let mut new_keys: BTreeSet<Key> = BTreeSet::new();
    for k in maxlev.keys() {
        if !matches!(k, Key::Corner(_)) {
            new_keys.insert(*k);
        }
    }
    let mut points = pts.clone();
    let mut index: HashMap<Key, u32> = HashMap::with_capacity(maxlev.len());
    for k in maxlev.keys() {
        if let Key::Corner(n) = k {
            index.insert(*k, *n);
        }
    }
    let rf = r as f64;
    for k in &new_keys {
        let x = match *k {
            Key::Corner(n) => pts[n as usize],
            Key::Edge(a, b, t) => lerp(pts[a as usize], pts[b as usize], t as f64 / rf),
            Key::Face([o, u, v, w], cu, cv) => {
                let s = cu as f64 / rf;
                let t = cv as f64 / rf;
                lerp(
                    lerp(pts[o as usize], pts[u as usize], s),
                    lerp(pts[v as usize], pts[w as usize], s),
                    t,
                )
            }
            Key::Cell(c, p) => trilinear(
                pts,
                tree.hex[c as usize].as_ref().unwrap(),
                std::array::from_fn(|d| p[d] as f64 / rf),
            ),
        };
        index.insert(*k, points.len() as u32);
        points.push(x);
    }
    // the input edges that carry new nodes
    let mut split_edges: HashSet<(u32, u32)> = HashSet::new();
    for k in maxlev.keys() {
        if let Key::Edge(a, b, _) = k {
            split_edges.insert((*a, *b));
        }
    }
    let mut level = vec![0u8; n_out];
    let mut parent = vec![0u32; n_out];
    for c in 0..nc {
        for (i, l) in tree.leaves[c].iter().enumerate() {
            level[first[c] as usize + i] = l.level;
            parent[first[c] as usize + i] = c as u32;
        }
    }
    let mut b = MeshBuilder::new(points, n_out);
    for p in &mesh.boundaries {
        b.declare_boundary(&p.name);
    }
    b.set_zones(
        mesh.zones.clone(),
        (0..n_out)
            .map(|i| mesh.cell_zone[parent[i] as usize])
            .collect(),
    );
    let refined = |c: usize| tree.leaves[c].len() > 1;
    let mut with_hanging = 0usize;
    // an input polygon with the new nodes on its edges
    let insert_edges = |poly: &[u32]| -> Vec<u32> {
        let mut out = Vec::with_capacity(poly.len());
        for i in 0..poly.len() {
            let (a, bb) = (poly[i], poly[(i + 1) % poly.len()]);
            out.push(a);
            if split_edges.contains(&(a.min(bb), a.max(bb))) {
                for t in 1..r {
                    if let Some(&n) = index.get(&edge_key(a, bb, t, r)) {
                        out.push(n);
                    }
                }
            }
        }
        out
    };
    // a rectangle of a frame, (u0, v0) → (u1, v0) → (u1, v1) → (u0, v1), with
    // every node on its edges
    let rect = |key: &dyn Fn(u32, u32) -> Key, u0: u32, u1: u32, v0: u32, v1: u32| -> Vec<u32> {
        let corners = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)];
        let mut out = Vec::with_capacity(8);
        for i in 0..4 {
            let (pa, pb) = (corners[i], corners[(i + 1) % 4]);
            out.push(index[&key(pa.0, pa.1)]);
            let len = pa.0.abs_diff(pb.0).max(pa.1.abs_diff(pb.1));
            for t in 1..len {
                let q = (
                    if pb.0 >= pa.0 {
                        pa.0 + t * (pb.0 != pa.0) as u32
                    } else {
                        pa.0 - t
                    },
                    if pb.1 >= pa.1 {
                        pa.1 + t * (pb.1 != pa.1) as u32
                    } else {
                        pa.1 - t
                    },
                );
                if let Some(&n) = index.get(&key(q.0, q.1)) {
                    out.push(n);
                }
            }
        }
        out
    };
    let overlap = |a: (u32, u32, u32, u32), b: (u32, u32, u32, u32)| {
        a.0 < b.1 && b.0 < a.1 && a.2 < b.3 && b.2 < a.3
    };
    let smaller = |a: (u32, u32, u32, u32), b: (u32, u32, u32, u32)| {
        if (a.1 - a.0) * (a.3 - a.2) <= (b.1 - b.0) * (b.3 - b.2) {
            a
        } else {
            b
        }
    };
    // 1. the faces of the input mesh
    let ni = mesh.n_interior_faces();
    for f in 0..mesh.n_faces() {
        let poly = mesh.face(f);
        let a = mesh.owner[f] as usize;
        let nb = (f < ni).then(|| mesh.neighbour[f] as usize);
        let bname = (f >= ni).then(|| &mesh.boundaries[mesh.boundary_of(f).unwrap()].name);
        if !refined(a) && !nb.is_some_and(refined) {
            let p = insert_edges(poly);
            if p.len() > poly.len() {
                with_hanging += 1;
            }
            match (nb, bname) {
                (Some(n), _) => b.add_interior(&p, first[a], first[n]),
                (None, Some(name)) => b.add_boundary(&p, first[a], name),
                _ => unreachable!(),
            }
            continue;
        }
        // a quadrilateral of a refined hexahedron: the tilings of its sides
        let io = (0..4).min_by_key(|&i| poly[i]).unwrap();
        let (prev, next) = (poly[(io + 3) % 4], poly[(io + 1) % 4]);
        let fr: Frame = if next < prev {
            [poly[io], next, prev, poly[(io + 2) % 4]]
        } else {
            [poly[io], prev, next, poly[(io + 2) % 4]]
        };
        // +1: the frame order is the order of the polygon (normal out of a)
        let same = next < prev;
        let tiles = |c: usize| -> Vec<((u32, u32, u32, u32), u32)> {
            if !refined(c) {
                return vec![((0, r, 0, r), first[c])];
            }
            let h = tree.hex[c].as_ref().unwrap();
            let mut set = [0u32; 4];
            set.copy_from_slice(poly);
            set.sort_unstable();
            for fx in 0..3 {
                for side in [0, r] {
                    let (lf, map) = tree.local_frame(c, fx, side);
                    let mut s2 = lf;
                    s2.sort_unstable();
                    if s2 != set {
                        continue;
                    }
                    debug_assert_eq!(lf, fr);
                    let _ = h;
                    let (g, k) = match fx {
                        0 => (1, 2),
                        1 => (0, 2),
                        _ => (0, 1),
                    };
                    let mut out = Vec::new();
                    for (i, l) in tree.leaves[c].iter().enumerate() {
                        let on = if side == 0 {
                            l.lo[fx] == 0
                        } else {
                            l.lo[fx] + l.n[fx] == r
                        };
                        if !on {
                            continue;
                        }
                        let p0 = map(l.lo[g], l.lo[k], r);
                        let p1 = map(l.lo[g] + l.n[g], l.lo[k] + l.n[k], r);
                        out.push((
                            (
                                p0.0.min(p1.0),
                                p0.0.max(p1.0),
                                p0.1.min(p1.1),
                                p0.1.max(p1.1),
                            ),
                            first[c] + i as u32,
                        ));
                    }
                    return out;
                }
            }
            unreachable!("a face of a refined hexahedron")
        };
        let key = |u: u32, v: u32| frame_key(&fr, u, v, r);
        let ta = tiles(a);
        let oriented = |mut p: Vec<u32>| {
            if !same {
                p[1..].reverse();
            }
            p
        };
        match nb {
            Some(n) => {
                let tb = tiles(n);
                for (ra, ca) in &ta {
                    for (rb, cb) in &tb {
                        if overlap(*ra, *rb) {
                            let s = smaller(*ra, *rb);
                            let p = oriented(rect(&key, s.0, s.1, s.2, s.3));
                            if p.len() > 4 {
                                with_hanging += 1;
                            }
                            b.add_interior(&p, *ca, *cb);
                        }
                    }
                }
            }
            None => {
                for (ra, ca) in &ta {
                    let p = oriented(rect(&key, ra.0, ra.1, ra.2, ra.3));
                    if p.len() > 4 {
                        with_hanging += 1;
                    }
                    b.add_boundary(&p, *ca, bname.unwrap());
                }
            }
        }
    }
    // 2. the planes inside each refined hexahedron
    for c in 0..nc {
        if !refined(c) {
            continue;
        }
        let h = tree.hex[c].as_ref().unwrap();
        let p = |i: usize| pts[h[i] as usize];
        let mean = |d: usize, side: usize| {
            let mut s = [0.0; 3];
            for i in 0..8 {
                if (i >> d) & 1 == side {
                    let x = p(i);
                    s = [s[0] + x[0], s[1] + x[1], s[2] + x[2]];
                }
            }
            s
        };
        let e: Vec<[f64; 3]> = (0..3).map(|d| sub(mean(d, 1), mean(d, 0))).collect();
        let right = dot(e[0], cross(e[1], e[2])) > 0.0;
        // (axis, position) → (cells below, cells above)
        type Side = Vec<((u32, u32, u32, u32), u32)>;
        let mut planes: std::collections::BTreeMap<(usize, u32), (Side, Side)> = Default::default();
        for (i, l) in tree.leaves[c].iter().enumerate() {
            for f in 0..3 {
                if !tree.split[c][f] {
                    continue;
                }
                let (g, k) = ((f + 1) % 3, (f + 2) % 3);
                let rr = (l.lo[g], l.lo[g] + l.n[g], l.lo[k], l.lo[k] + l.n[k]);
                let id = first[c] + i as u32;
                let hi = l.lo[f] + l.n[f];
                if hi < r {
                    planes.entry((f, hi)).or_default().0.push((rr, id));
                }
                if l.lo[f] > 0 {
                    planes.entry((f, l.lo[f])).or_default().1.push((rr, id));
                }
            }
        }
        for ((f, pos), (below, above)) in &planes {
            let (f, pos) = (*f, *pos);
            let (g, k) = ((f + 1) % 3, (f + 2) % 3);
            let key = |u: u32, v: u32| {
                let mut q = [0u32; 3];
                q[f] = pos;
                q[g] = u;
                q[k] = v;
                tree.cell_key(c, q)
            };
            for (ra, ca) in below {
                for (rb, cb) in above {
                    if overlap(*ra, *rb) {
                        let s = smaller(*ra, *rb);
                        let mut poly = rect(&key, s.0, s.1, s.2, s.3);
                        if poly.len() > 4 {
                            with_hanging += 1;
                        }
                        if !right {
                            poly[1..].reverse();
                        }
                        b.add_interior(&poly, *ca, *cb);
                    }
                }
            }
        }
    }
    let out = b.finish()?;
    let mut per_level = vec![0usize; 1 + level.iter().copied().max().unwrap_or(0) as usize];
    for &l in &level {
        per_level[l as usize] += 1;
    }
    let hanging_faces = (0..out.n_interior_faces())
        .filter(|&f| level[out.owner[f] as usize] != level[out.neighbour[f] as usize])
        .count();
    // cells of the result are numbered as given to the builder
    let report = RefineReport {
        cells_in: nc,
        cells_out: n_out,
        refined: (0..nc).filter(|&c| refined(c)).count(),
        not_hexahedra: tree.hex.iter().filter(|h| h.is_none()).count(),
        per_level,
        balance,
        hanging_faces,
        faces_with_hanging_nodes: with_hanging,
    };
    Ok(Refined {
        mesh: out,
        parent,
        level,
        report,
    })
}
