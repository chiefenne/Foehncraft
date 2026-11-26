//! Conservative cut cells (L13 stage 2, D227–D229): a body is cut out of
//! any polyhedral mesh, the fluid parts of the cut cells stay ordinary
//! polyhedral cells (D23), small ones are merged into a neighbour and the
//! body surface becomes an ordinary wall boundary.
//!
//! **Input.** The body is an implicit function φ(x) ([`ImplicitBody`]),
//! negative inside the solid, zero on its surface, positive in the fluid
//! (the surface counts as solid, as in the inside test of the penalised
//! bodies). For a signed distance the values are lengths, but any level
//! set works. The value is taken at the mesh nodes; where an edge runs from
//! a fluid node to a solid one, the crossing point is
//! [`ImplicitBody::edge_point`]: the zero of the linear interpolant
//! (marching-cubes style), or the root of the body's own function on the
//! edge, which puts the point on the true surface.
//!
//! **Faces.** Every face is clipped to its fluid part: the ring of its
//! nodes keeps the fluid nodes and puts the crossing point of each mixed
//! edge in between. A crossing point belongs to the edge, so the faces of
//! an edge share it, and the clipped faces of two cells match exactly. The
//! ring edge from the point where the ring leaves the fluid to the point
//! where it comes back is a *chord*: a straight piece of the cut surface
//! across the face.
//!
//! **Cells.** The chords of the faces of a cell, reversed, form closed
//! loops: the wall faces that close the fluid part of the cell. The cell
//! keeps its clipped faces plus these wall faces; as every polygon's area
//! vector is ½ Σ xᵢ × xᵢ₊₁ over its ring (the fan of the geometry module
//! gives exactly that), the area vectors of a cut cell sum to zero to
//! round-off, whatever the shape of the wall faces (closure, the
//! conservation of the discretisation). A cell whose fluid falls into
//! pieces that share no point becomes one cell per piece.
//!
//! **Wall faces.** A wall loop is one (generally non-planar) face, or, with
//! [`CutOptions::wall_point`], a fan of triangles around a point put on
//! the surface by [`ImplicitBody::project`] (in the spirit of the vertex of
//! dual contouring: the face follows a curved surface inside the cell). Each
//! wall face goes to the boundary of the body that is nearest to its
//! centre ([`ImplicitBody::patch`]).
//!
//! **Small cells (D228, owner decision Q3: merge).** A cut cell whose
//! volume is below [`CutOptions::min_fraction`] of its uncut cell, or that
//! has a face whose normal does not point from its centroid to the
//! neighbour's (the discretisation needs d · S > 0), is merged into the
//! neighbour with which it shares the largest face area. A merged cell is a
//! polyhedron; the faces it shares with one neighbour are joined into one
//! face along their common edges (two cells share at most one face). A
//! merge that would give a face that is not one ring is not made, the next
//! neighbour is tried, and a cell that cannot be merged at all is kept and
//! counted in the report.
//!
//! **Determinism (D93).** The cut is a serial mesh operation on the global
//! mesh, before renumbering and partitioning, in face and cell order: the
//! cut mesh is bitwise the same for any number of threads and ranks.

use crate::builder::MeshBuilder;
use crate::geometry::{add, dot, norm, polygon_area_centroid, scale, sub};
use crate::{Geometry, Mesh};
use std::collections::HashMap;

/// A body given by an implicit function φ (negative inside, positive in
/// the fluid).
pub trait ImplicitBody: Sync {
    /// φ(x).
    fn value(&self, x: [f64; 3]) -> f64;

    /// φ at the mesh nodes. Node `i` is farther than `reach[i]` from the
    /// surface only if no mesh edge from it reaches the surface, so beyond
    /// that distance only the sign of its value matters (a body can give
    /// a narrow-band distance there). The default is [`Self::value`] at
    /// every point.
    fn values(&self, points: &[[f64; 3]], reach: &[f64]) -> Vec<f64> {
        use rayon::prelude::*;
        let _ = reach;
        points.par_iter().map(|&x| self.value(x)).collect()
    }

    /// The crossing point on the edge from fluid node `a` (φ_a > 0) to
    /// solid node `b` (φ_b ≤ 0). The default is the zero of the linear
    /// interpolant.
    fn edge_point(&self, a: [f64; 3], b: [f64; 3], fa: f64, fb: f64) -> [f64; 3] {
        let t = fa / (fa - fb);
        add(a, scale(sub(b, a), t))
    }

    /// A point on the surface near `x`, for the wall point of
    /// [`CutOptions::wall_point`]; `None` keeps the wall face whole.
    fn project(&self, _x: [f64; 3]) -> Option<[f64; 3]> {
        None
    }

    /// The point a wall face is fanned around ([`CutOptions::wall_point`]):
    /// `ring` is the closed loop of the wall face (points on the surface),
    /// `centre` its centroid. The default projects the centroid onto the
    /// surface ([`Self::project`]); a body with sharp edges can put the
    /// point on the edge or corner that runs through the cell (the vertex
    /// of dual contouring, [`feature_point`]).
    fn wall_point(&self, ring: &[[f64; 3]], centre: [f64; 3]) -> Option<[f64; 3]> {
        let _ = ring;
        self.project(centre)
    }

    /// A point on the surface between the ends `p` and `q` of a chord, the
    /// straight cut across a face with unit normal `n`
    /// ([`CutOptions::face_point`]): in the plane of the face, where the
    /// surface crosses the perpendicular bisector of the chord, or where a
    /// sharp edge of the body pierces the face. `None` keeps the chord
    /// straight. The cut checks that the face stays a simple polygon inside
    /// the uncut face and drops the point otherwise.
    fn face_point(&self, p: [f64; 3], q: [f64; 3], n: [f64; 3]) -> Option<[f64; 3]> {
        let _ = (p, q, n);
        None
    }

    /// The body (index into the names of [`cut`]) whose surface is
    /// nearest to `x`.
    fn patch(&self, _x: [f64; 3]) -> usize {
        0
    }
}

/// The root of φ on the edge from `a` (φ_a > 0) to `b` (φ_b ≤ 0) by the
/// Illinois variant of regula falsi, a fixed number of steps at most, so
/// the point is the same on every machine: the point on the true surface
/// for [`ImplicitBody::edge_point`] of a body with an exact function.
pub fn edge_root<F: Fn([f64; 3]) -> f64>(
    f: F,
    a: [f64; 3],
    b: [f64; 3],
    fa: f64,
    fb: f64,
) -> [f64; 3] {
    if fb == 0.0 {
        return b;
    }
    let d = sub(b, a);
    let at = |t: f64| add(a, scale(d, t));
    let (mut t0, mut t1, mut f0, mut f1) = (0.0f64, 1.0f64, fa, fb);
    let mut side = 0i8;
    for _ in 0..80 {
        let t = (t0 * f1 - t1 * f0) / (f1 - f0);
        if !(t > t0 && t < t1) {
            break;
        }
        let ft = f(at(t));
        if ft == 0.0 {
            return at(t);
        }
        if ft > 0.0 {
            t0 = t;
            f0 = ft;
            if side == 1 {
                f1 *= 0.5;
            }
            side = 1;
        } else {
            t1 = t;
            f1 = ft;
            if side == -1 {
                f0 *= 0.5;
            }
            side = -1;
        }
        if t1 - t0 <= 1e-15 {
            break;
        }
    }
    // the end on the solid side: the surface counts as solid
    let t = (t0 * f1 - t1 * f0) / (f1 - f0);
    at(t.clamp(t0, t1))
}

/// A point on the zero level of `f` near `x` by Newton steps along the
/// gradient (central differences with step `h`); `None` if it does not
/// settle.
pub fn project_newton<F: Fn([f64; 3]) -> f64>(f: F, x: [f64; 3], h: f64) -> Option<[f64; 3]> {
    let mut p = x;
    for _ in 0..8 {
        let v = f(p);
        if v == 0.0 {
            return Some(p);
        }
        let mut g = [0.0; 3];
        for (i, gi) in g.iter_mut().enumerate() {
            let mut a = p;
            let mut b = p;
            a[i] += h;
            b[i] -= h;
            *gi = (f(a) - f(b)) / (2.0 * h);
        }
        let gg = dot(g, g);
        if gg.is_nan() || gg <= 0.0 || !gg.is_finite() {
            return None;
        }
        let step = scale(g, v / gg);
        p = sub(p, step);
        if norm(step) <= 1e-13 * (1.0 + norm(p)) {
            break;
        }
    }
    let v = f(p);
    (v.abs() <= 1e-9 * (1.0 + norm(p)) || v.abs() <= 1e-3 * h).then_some(p)
}

/// Opening angle of the surface normals at the points of a wall loop or
/// the ends of a chord above which the cut looks for a sharp edge or
/// corner (the feature test of extended marching cubes, Kobbelt et al.
/// 2001, there with cos θ < 0.9): 15°.
pub const FEATURE_COS: f64 = 0.965_925_826_289_068_3;

/// The vertex of dual contouring (Ju, Losasso, Schaefer & Warren 2002) for
/// a wall loop: the point `x` that minimises Σ (nᵢ · (x − pᵢ))² over the
/// loop points `p` with their surface normals `n`, the quadratic error
/// function (QEF), solved with the pseudo-inverse about the mass point of
/// the loop (Lindstrom 2000; Schaefer & Warren 2002, eigenvalues below
/// 1 % of the largest dropped). It lies on the sharp edge or corner where
/// the planes of the normals meet. `None` when the normals open less than
/// [`FEATURE_COS`] (a smooth piece, where the projection of the centroid is
/// better) or the point leaves the loop's ball about `centre`. The caller
/// checks that the point lies on its surface.
pub fn feature_point(p: &[[f64; 3]], n: &[[f64; 3]], centre: [f64; 3]) -> Option<[f64; 3]> {
    let k = p.len();
    if k < 2 || n.len() != k {
        return None;
    }
    let mut cmin = 1.0f64;
    for i in 0..k {
        for j in i + 1..k {
            cmin = cmin.min(dot(n[i], n[j]));
        }
    }
    if cmin.is_nan() || cmin > FEATURE_COS {
        return None;
    }
    let mut m = [0.0; 3];
    for &x in p {
        m = add(m, x);
    }
    let m = scale(m, 1.0 / k as f64);
    // AᵀA and Aᵀ(b − A m), about the mass point
    let mut ata = [[0.0f64; 3]; 3];
    let mut atb = [0.0f64; 3];
    for (x, nn) in p.iter().zip(n) {
        let r = dot(*nn, sub(*x, m));
        for a in 0..3 {
            atb[a] += nn[a] * r;
            for b in 0..3 {
                ata[a][b] += nn[a] * nn[b];
            }
        }
    }
    let (lam, vec) = sym_eigen3(ata);
    let lmax = lam.iter().fold(0.0f64, |a, &b| a.max(b));
    if lmax <= 0.0 || !lmax.is_finite() {
        return None;
    }
    let mut x = m;
    for j in 0..3 {
        if lam[j] > 1e-2 * lmax {
            let v = [vec[0][j], vec[1][j], vec[2][j]];
            x = add(x, scale(v, dot(v, atb) / lam[j]));
        }
    }
    let r = p.iter().map(|&y| norm(sub(y, centre))).fold(0.0, f64::max);
    (norm(sub(x, centre)) <= r && x.iter().all(|v| v.is_finite())).then_some(x)
}

/// Eigenvalues and eigenvectors (columns) of a symmetric 3 × 3 matrix by
/// cyclic Jacobi rotations (a fixed number of sweeps at most: the same
/// result on every machine).
pub fn sym_eigen3(a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = a;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..32 {
        let off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        let diag = a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2];
        if off <= 1e-30 * diag || off == 0.0 {
            break;
        }
        for (p, q) in [(0usize, 1usize), (0, 2), (1, 2)] {
            if a[p][q] == 0.0 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            // A ← Jᵀ A J
            for k in 0..3 {
                let (akp, akq) = (a[k][p], a[k][q]);
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for k in 0..3 {
                let (apk, aqk) = (a[p][k], a[q][k]);
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for row in &mut v {
                let (vp, vq) = (row[p], row[q]);
                row[p] = c * vp - s * vq;
                row[q] = s * vp + c * vq;
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}

/// The face point of [`ImplicitBody::face_point`] for a body given by its
/// function `f` and its unit normal `normal` (out of the body): where two
/// sharp sides meet in the face (the normals at `p` and `q`, projected into
/// the face, open more than [`FEATURE_COS`]: the crossing of their two
/// lines, if it lies on the surface), else the root of `f` on the
/// perpendicular bisector of the chord nearest to its middle. `None` if
/// the point would be on the chord (a flat surface) or no root is found.
pub fn chord_point<F, N>(f: F, normal: N, p: [f64; 3], q: [f64; 3], n: [f64; 3]) -> Option<[f64; 3]>
where
    F: Fn([f64; 3]) -> f64,
    N: Fn([f64; 3]) -> Option<[f64; 3]>,
{
    let e = sub(q, p);
    let l = norm(e);
    if l == 0.0 || !l.is_finite() {
        return None;
    }
    let eu = scale(e, 1.0 / l);
    let d = crate::geometry::cross(n, eu);
    let dl = norm(d);
    if dl == 0.0 {
        return None;
    }
    let d = scale(d, 1.0 / dl);
    let c = scale(add(p, q), 0.5);
    let in_plane = |v: [f64; 3]| {
        let w = sub(v, scale(n, dot(v, n)));
        let wl = norm(w);
        (wl > 0.0).then(|| scale(w, 1.0 / wl))
    };
    // a sharp edge through the face: the crossing of the two lines
    if let (Some(n1), Some(n2)) = (normal(p), normal(q))
        && let (Some(a1), Some(a2)) = (in_plane(n1), in_plane(n2))
        && dot(a1, a2) < FEATURE_COS
    {
        // x = c + α eu + β d with a_i · (x − end_i) = 0
        let (m11, m12, r1) = (dot(a1, eu), dot(a1, d), dot(a1, sub(p, c)));
        let (m21, m22, r2) = (dot(a2, eu), dot(a2, d), dot(a2, sub(q, c)));
        let det = m11 * m22 - m12 * m21;
        if det.abs() > 1e-12 {
            let al = (r1 * m22 - m12 * r2) / det;
            let be = (m11 * r2 - r1 * m21) / det;
            let x = add(c, add(scale(eu, al), scale(d, be)));
            if norm(sub(x, c)) <= l && f(x).abs() <= 1e-7 * l {
                return (be.abs() > 1e-3 * l).then_some(x);
            }
        }
    }
    // the root on the bisector nearest to the middle of the chord
    let fc = f(c);
    if fc == 0.0 {
        return None;
    }
    let mut best: Option<[f64; 3]> = None;
    for s in [1.0, -1.0] {
        let end = add(c, scale(d, s * 0.5 * l));
        let fe = f(end);
        if fe.is_nan() || (fe > 0.0) == (fc > 0.0) {
            continue;
        }
        let x = if fc > 0.0 {
            edge_root(&f, c, end, fc, fe)
        } else {
            edge_root(&f, end, c, fe, fc)
        };
        if best.is_none_or(|b| norm(sub(x, c)) < norm(sub(b, c))) {
            best = Some(x);
        }
    }
    best.filter(|x| dot(sub(*x, c), d).abs() > 1e-3 * l)
}

/// Nodes closer to the surface than this fraction of their shortest edge
/// are put on it.
const SNAP: f64 = 1e-8;

/// Options of [`cut`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CutOptions {
    /// Cut cells with a smaller volume fraction (of their uncut cell) are
    /// merged into a neighbour; 0 merges only cells the discretisation
    /// cannot take.
    pub min_fraction: f64,
    /// Fan each wall face around a point on the surface (needs
    /// [`ImplicitBody::wall_point`] or [`ImplicitBody::project`]).
    pub wall_point: bool,
    /// Bend each chord of a cut face through a point on the surface
    /// ([`ImplicitBody::face_point`]): the clipped face follows a curved
    /// body or a sharp edge through the face instead of cutting straight
    /// across.
    pub face_point: bool,
}

impl Default for CutOptions {
    fn default() -> Self {
        Self {
            min_fraction: 0.2,
            wall_point: false,
            face_point: false,
        }
    }
}

/// What the cut did, for the log, the limitations and `result.json`.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct CutReport {
    /// Cells of the input mesh.
    pub cells_in: usize,
    /// Cells of the cut mesh.
    pub cells_out: usize,
    /// Input cells cut by the surface (some nodes inside, some outside).
    pub cells_cut: usize,
    /// Input cells wholly inside the body, removed.
    pub cells_removed: usize,
    /// Extra cells from cut cells whose fluid falls into several pieces.
    pub cells_split: usize,
    /// Small or ill-shaped cut cells merged into a neighbour.
    pub cells_merged: usize,
    /// Small or ill-shaped cut cells that could not be merged (kept).
    pub merge_failed: usize,
    /// Smallest volume fraction (of its uncut cell) among the cut cells
    /// that were not merged.
    pub min_fraction: f64,
    /// Faces of the cut mesh where d · S ≤ 0 (the discretisation refuses
    /// them); zero on a usable mesh.
    pub bad_faces: usize,
    /// Fluid volume of the cut mesh (m³).
    pub fluid_volume: f64,
    /// Wall area per body (m²), in the order of the names.
    pub wall_area: Vec<f64>,
    /// Largest |Σ S| of a cell over the largest face area of that cell.
    pub max_closure: f64,
    /// Joined faces of merged cells that fold: their fan moves more than
    /// 1 % of the smaller cell's volume into the other cell (the merge
    /// picks a neighbour without a fold where it can).
    pub folded_faces: usize,
    /// The largest volume a joined face moves, over the smaller cell.
    pub max_fold: f64,
    /// Input cells cut with fewer points than asked for, so that the
    /// discretisation takes their faces: the face points dropped first,
    /// then the sharp-feature wall point, then the projected one (the
    /// fallback runs where a face has d · S ≤ 0).
    pub simplified: usize,
    /// Rounds of the cut (one, plus one per fallback step).
    pub rounds: usize,
}

/// The mesh with the bodies cut out and their report; the cells of the cut
/// mesh know their cell in the input mesh (`parent`).
#[derive(Debug, Clone)]
pub struct CutMesh {
    /// The cut mesh: the input boundaries first, then one wall boundary per
    /// body (named by `names`, empty ones included).
    pub mesh: Mesh,
    /// Input cell of each cut-mesh cell (for a merged cell, the cell it
    /// was merged into).
    pub parent: Vec<u32>,
    /// The input cells whose fluid makes up each cut-mesh cell, in
    /// increasing order (one for a plain or split cell, several for a
    /// merged one).
    pub covers: Vec<Vec<u32>>,
    /// The fluid area of every input face as the cut clipped it (zero for
    /// a face inside the body): the aperture times the face area.
    pub face_fluid_area: Vec<f64>,
    /// The report.
    pub report: CutReport,
}

/// A face during the cut: its ring (outward from `a`) and its two cells
/// (`b`: neighbour cell or `!boundary`).
#[derive(Debug, Clone)]
struct Face {
    ring: Vec<u32>,
    a: u32,
    b: i64,
}

impl Face {
    fn boundary(&self) -> Option<usize> {
        (self.b < 0).then(|| (-self.b - 1) as usize)
    }
}

/// Cuts the bodies of `body` out of `mesh`. `names` names the body
/// boundaries; the mesh must be a global mesh (no halo) without periodic
/// faces (they are paired after the cut).
pub fn cut(
    mesh: &Mesh,
    body: &dyn ImplicitBody,
    names: &[String],
    opt: &CutOptions,
) -> Result<CutMesh, String> {
    if mesh.n_halo != 0 {
        return Err("the cut runs on the global mesh, not on a rank's part".into());
    }
    if mesh.periodic.is_some() {
        return Err("the cut runs before the periodic pairing".into());
    }
    if names.is_empty() {
        return Err("no body to cut".into());
    }
    for n in names {
        if mesh.boundary(n).is_some() {
            return Err(format!(
                "the mesh already has a boundary named `{n}`, the name of a body"
            ));
        }
    }
    let prep = Prep::new(mesh, body)?;
    // the per-cell fallback: cells next to a face the discretisation
    // refuses lose the finest point first (face points, then the
    // sharp-feature wall point, then the projected wall point)
    let first = if opt.face_point { 0 } else { 1 };
    let last = if opt.wall_point { 3 } else { 1 };
    let mut level = vec![first; mesh.n_cells];
    let mut edges: HashMap<(u32, u32), [f64; 3]> = HashMap::new();
    let mut rounds = 0;
    loop {
        let mut out = cut_with(mesh, body, names, opt, &prep, &level, &mut edges)?;
        rounds += 1;
        let mut raised = false;
        if out.report.bad_faces > 0 && rounds < MAX_ROUNDS {
            let geo = Geometry::compute(&out.mesh);
            for (c, bad) in bad_cells(&out.mesh, &geo).into_iter().enumerate() {
                if !bad {
                    continue;
                }
                for &p in &out.covers[c] {
                    let l = &mut level[p as usize];
                    if *l < last && prep.cut[p as usize] {
                        *l += 1;
                        raised = true;
                    }
                }
            }
        }
        if !raised {
            out.report.simplified = level.iter().filter(|&&l| l > first).count();
            out.report.rounds = rounds;
            return Ok(out);
        }
    }
}

/// Rounds of the per-cell fallback of [`cut`].
const MAX_ROUNDS: usize = 8;

/// What every round of the cut shares: the uncut geometry, the node
/// values and the input cells the surface cuts.
struct Prep {
    geo0: Geometry,
    phi: Vec<f64>,
    cut: Vec<bool>,
}

impl Prep {
    fn new(mesh: &Mesh, body: &dyn ImplicitBody) -> Result<Self, String> {
        let geo0 = Geometry::compute(mesh);
        // the shortest and longest edge at every node: the snap below and
        // the reach of the node values
        let mut hmin = vec![f64::INFINITY; mesh.points.len()];
        let mut hmax = vec![0.0f64; mesh.points.len()];
        for f in 0..mesh.n_faces() {
            let nodes = mesh.face(f);
            for i in 0..nodes.len() {
                let (a, b) = (nodes[i] as usize, nodes[(i + 1) % nodes.len()] as usize);
                let l = norm(sub(mesh.points[a], mesh.points[b]));
                hmin[a] = hmin[a].min(l);
                hmin[b] = hmin[b].min(l);
                hmax[a] = hmax[a].max(l);
                hmax[b] = hmax[b].max(l);
            }
        }
        let reach: Vec<f64> = hmax
            .iter()
            .map(|&h| if h > 0.0 { 1.01 * h } else { f64::INFINITY })
            .collect();
        let phi = body.values(&mesh.points, &reach);
        if phi.len() != mesh.points.len() {
            return Err("the body gave a value count other than the node count".into());
        }
        if let Some(i) = phi.iter().position(|v| v.is_nan()) {
            return Err(format!("the body function is not a number at node {i}"));
        }
        // a node within round-off of the surface lies on it (solid, its own
        // crossing point): no crossing points a round-off away from a node
        let phi: Vec<f64> = phi
            .iter()
            .zip(&hmin)
            .map(|(&v, &h)| if v.abs() <= SNAP * h { 0.0 } else { v })
            .collect();
        let cut = (0..mesh.n_cells)
            .map(|c| {
                let (mut fl, mut so) = (false, false);
                for (f, _) in mesh.cell_face_iter(c) {
                    for &n in mesh.face(f) {
                        if phi[n as usize] > 0.0 {
                            fl = true;
                        } else {
                            so = true;
                        }
                    }
                }
                fl && so
            })
            .collect();
        Ok(Self { geo0, phi, cut })
    }
}

/// One round of [`cut`]: `level` per input cell, 0 all points, 1 no face
/// points, 2 the projected wall point only, 3 whole wall faces; a face
/// takes face points only when both its cells are at level 0. `edges`
/// keeps the crossing point of every edge between the rounds.
fn cut_with(
    mesh: &Mesh,
    body: &dyn ImplicitBody,
    names: &[String],
    opt: &CutOptions,
    prep: &Prep,
    level: &[u8],
    edges: &mut HashMap<(u32, u32), [f64; 3]>,
) -> Result<CutMesh, String> {
    let geo0 = &prep.geo0;
    let phi = &prep.phi;
    let fluid = |n: u32| phi[n as usize] > 0.0;
    let mut points = mesh.points.clone();
    let mut crossing: HashMap<(u32, u32), u32> = HashMap::new();
    // the crossing point of edge (a fluid, b solid), shared by the faces
    let mut cross_pt = |a: u32, b: u32, points: &mut Vec<[f64; 3]>| -> u32 {
        if phi[b as usize] == 0.0 {
            return b;
        }
        *crossing.entry((a, b)).or_insert_with(|| {
            let p = *edges.entry((a, b)).or_insert_with(|| {
                body.edge_point(
                    points[a as usize],
                    points[b as usize],
                    phi[a as usize],
                    phi[b as usize],
                )
            });
            points.push(p);
            (points.len() - 1) as u32
        })
    };

    // 1. clip the faces: ring outward from the owner, chords marked
    let ni = mesh.n_interior_faces();
    let nf = mesh.n_faces();
    let mut faces: Vec<Option<(Vec<u32>, Vec<bool>)>> = Vec::with_capacity(nf);
    for f in 0..nf {
        let nodes = mesh.face(f);
        let nfl = nodes.iter().filter(|&&n| fluid(n)).count();
        if nfl == nodes.len() {
            faces.push(Some((nodes.to_vec(), vec![false; nodes.len()])));
            continue;
        }
        if nfl == 0 {
            faces.push(None);
            continue;
        }
        let k = nodes.len();
        // (point, the ring edge leaving it is a chord)
        let mut ring: Vec<(u32, bool)> = Vec::with_capacity(k + 2);
        for i in 0..k {
            let (a, b) = (nodes[i], nodes[(i + 1) % k]);
            match (fluid(a), fluid(b)) {
                (true, true) => ring.push((a, false)),
                (true, false) => {
                    ring.push((a, false));
                    let x = cross_pt(a, b, &mut points);
                    ring.push((x, true));
                }
                (false, true) => {
                    let x = cross_pt(b, a, &mut points);
                    ring.push((x, false));
                }
                (false, false) => {}
            }
        }
        dedupe(&mut ring);
        if ring.len() < 3 {
            faces.push(None);
            continue;
        }
        // face points where both cells take them
        if opt.face_point
            && level[mesh.owner[f] as usize] == 0
            && (f >= ni || level[mesh.neighbour[f] as usize] == 0)
        {
            bend_chords(&mut ring, nodes, &mut points, body);
        }
        let (r, c) = ring.into_iter().unzip();
        faces.push(Some((r, c)));
    }
    let face_fluid_area: Vec<f64> = faces
        .iter()
        .map(|f| {
            f.as_ref()
                .map_or(0.0, |(r, _)| norm(polygon_area_centroid(&points, r).0))
        })
        .collect();

    // 2. the cells: wall loops from the chords, pieces
    let mut out_faces: Vec<Face> = Vec::new();
    let mut walls: Vec<(u32, Vec<u32>)> = Vec::new(); // (new cell, loop outward)
    let mut parent: Vec<u32> = Vec::new();
    let mut cut_cell: Vec<bool> = Vec::new();
    // new cell of each (input cell, first point of a ring) is found by the
    // piece of that point
    let mut face_cell: Vec<[u32; 2]> = vec![[u32::MAX; 2]; nf];
    let mut report = CutReport {
        cells_in: mesh.n_cells,
        wall_area: vec![0.0; names.len()],
        ..Default::default()
    };
    for c in 0..mesh.n_cells {
        let mut has = false;
        let mut cut_here = false;
        let mut edges: Vec<(u32, u32)> = Vec::new();
        for (f, is_nb) in mesh.cell_face_iter(c) {
            let Some((ring, chord)) = &faces[f] else {
                cut_here = true;
                continue;
            };
            has = true;
            let k = ring.len();
            for i in 0..k {
                if chord[i] {
                    cut_here = true;
                    let (p, q) = (ring[i], ring[(i + 1) % k]);
                    // outward ring of the owner runs p → q, of the
                    // neighbour q → p; the wall runs the other way
                    edges.push(if is_nb { (p, q) } else { (q, p) });
                }
            }
        }
        if !has {
            report.cells_removed += 1;
            continue;
        }
        if cut_here {
            report.cells_cut += 1;
        }
        // pieces: union-find over the points of the rings and wall edges
        let mut uf = LocalUf::default();
        for (f, _) in mesh.cell_face_iter(c) {
            if let Some((ring, _)) = &faces[f] {
                for w in ring.windows(2) {
                    uf.union(w[0], w[1]);
                }
                uf.find(ring[0]);
            }
        }
        for &(p, q) in &edges {
            uf.union(p, q);
        }
        let mut roots: Vec<u32> = Vec::new();
        let mut piece_of = |p: u32, uf: &mut LocalUf| -> usize {
            let r = uf.find(p);
            match roots.iter().position(|&x| x == r) {
                Some(i) => i,
                None => {
                    roots.push(r);
                    roots.len() - 1
                }
            }
        };
        let base = parent.len() as u32;
        for (f, is_nb) in mesh.cell_face_iter(c) {
            if let Some((ring, _)) = &faces[f] {
                let k = piece_of(ring[0], &mut uf);
                face_cell[f][usize::from(is_nb)] = base + k as u32;
            }
        }
        let loops = chain(&edges);
        let mut loop_piece = Vec::with_capacity(loops.len());
        for l in &loops {
            loop_piece.push(piece_of(l[0], &mut uf));
        }
        let n_pieces = roots.len();
        report.cells_split += n_pieces - 1;
        for _ in 0..n_pieces {
            parent.push(c as u32);
            cut_cell.push(cut_here);
        }
        for (l, k) in loops.into_iter().zip(loop_piece) {
            if l.len() >= 3 {
                walls.push((base + k as u32, l));
            }
        }
    }
    for f in 0..nf {
        let Some((ring, _)) = faces[f].take() else {
            continue;
        };
        let a = face_cell[f][0];
        if f < ni {
            let b = face_cell[f][1];
            out_faces.push(Face {
                ring,
                a,
                b: b as i64,
            });
        } else {
            let p = mesh.boundary_of(f).expect("boundary face");
            out_faces.push(Face {
                ring,
                a,
                b: -(p as i64) - 1,
            });
        }
    }
    let nb0 = mesh.boundaries.len();
    for (c, l) in walls {
        let (_, x) = polygon_area_centroid(&points, &l);
        let pc = parent[c as usize] as usize;
        let centre = if opt.wall_point && level[pc] <= 2 {
            // the fan stays inside the uncut cell: a point outside it (a
            // sharp edge that passes the cell, the projection of a bent
            // loop) falls back to the projection, then to no fan
            let ok = |m: &[f64; 3]| inside_cell(mesh, geo0, pc, *m);
            let feature = if level[pc] <= 1 {
                let ring: Vec<[f64; 3]> = l.iter().map(|&i| points[i as usize]).collect();
                body.wall_point(&ring, x).filter(ok)
            } else {
                None
            };
            feature.or_else(|| body.project(x).filter(ok))
        } else {
            None
        };
        let pid = body.patch(x).min(names.len() - 1);
        if let Some(m) = centre {
            points.push(m);
            let mi = (points.len() - 1) as u32;
            let k = l.len();
            for i in 0..k {
                out_faces.push(Face {
                    ring: vec![l[i], l[(i + 1) % k], mi],
                    a: c,
                    b: -((nb0 + pid) as i64) - 1,
                });
            }
        } else {
            out_faces.push(Face {
                ring: l,
                a: c,
                b: -((nb0 + pid) as i64) - 1,
            });
        }
    }
    let n_cells = parent.len();
    if n_cells == 0 {
        return Err("the bodies fill the whole mesh".into());
    }
    let mut bnames: Vec<String> = mesh.boundaries.iter().map(|b| b.name.clone()).collect();
    bnames.extend(names.iter().cloned());

    // 3. small cells: merge, again while a merged cell is ill-shaped
    let mut covers: Vec<Vec<u32>> = parent.iter().map(|&p| vec![p]).collect();
    let mut min_fraction = opt.min_fraction;
    for _ in 0..10 {
        let (f, p, c, had_bad, gid) = merge_pass(
            mesh,
            &points,
            &out_faces,
            &parent,
            &cut_cell,
            geo0,
            &bnames,
            min_fraction,
            &mut report.cells_merged,
            &mut report.merge_failed,
            &mut report.folded_faces,
            &mut report.max_fold,
        )?;
        let changed = p.len() != parent.len();
        let mut next: Vec<Vec<u32>> = vec![Vec::new(); p.len()];
        for (old, g) in gid.iter().enumerate() {
            next[*g as usize].extend_from_slice(&covers[old]);
        }
        for v in &mut next {
            v.sort_unstable();
            v.dedup();
        }
        covers = next;
        out_faces = f;
        parent = p;
        cut_cell = c;
        // later passes merge only the cells the discretisation cannot take
        min_fraction = 0.0;
        if !had_bad && !changed {
            break;
        }
    }
    let final_parent = parent;
    let nfin = final_parent.len();
    let mut fin_faces = out_faces;
    let fzone: Vec<u32> = final_parent
        .iter()
        .map(|&p| mesh.cell_zone[p as usize])
        .collect();
    // drop unused points
    let mut used = vec![u32::MAX; points.len()];
    let mut fpoints = Vec::new();
    for f in &mut fin_faces {
        for n in &mut f.ring {
            if used[*n as usize] == u32::MAX {
                used[*n as usize] = fpoints.len() as u32;
                fpoints.push(points[*n as usize]);
            }
            *n = used[*n as usize];
        }
    }
    let out = build(&fpoints, nfin, &fin_faces, &bnames, &mesh.zones, &fzone)?;
    let geo = Geometry::compute(&out);
    report.cells_out = out.n_cells;
    report.fluid_volume = geo.cell_volume.iter().sum();
    report.bad_faces = bad_faces(&out, &geo);
    for (k, _) in names.iter().enumerate() {
        let b = &out.boundaries[nb0 + k];
        report.wall_area[k] = b.range().map(|f| norm(geo.face_area[f])).sum();
    }
    let mut minf = f64::INFINITY;
    for c in 0..nfin {
        if cut_cell[c] {
            minf = minf.min(geo.cell_volume[c] / geo0.cell_volume[final_parent[c] as usize]);
        }
        let mut amax = 0.0f64;
        for (f, _) in out.cell_face_iter(c) {
            amax = amax.max(norm(geo.face_area[f]));
        }
        report.max_closure = report.max_closure.max(norm(geo.closure(&out, c)) / amax);
    }
    report.min_fraction = if minf.is_finite() { minf } else { 1.0 };
    Ok(CutMesh {
        mesh: out,
        parent: final_parent,
        covers,
        face_fluid_area,
        report,
    })
}

fn build(
    points: &[[f64; 3]],
    n_cells: usize,
    faces: &[Face],
    bnames: &[String],
    zones: &[String],
    zone: &[u32],
) -> Result<Mesh, String> {
    let mut b = MeshBuilder::new(points.to_vec(), n_cells);
    for n in bnames {
        b.declare_boundary(n);
    }
    b.set_zones(zones.to_vec(), zone.to_vec());
    for f in faces {
        match f.boundary() {
            Some(p) => b.add_boundary(&f.ring, f.a, &bnames[p]),
            None => b.add_interior(&f.ring, f.a, f.b as u32),
        }
    }
    b.finish()
}

/// Cells with a face whose normal does not point from the cell centroid
/// to the neighbour centroid (or the face centroid at a boundary).
fn bad_cells(mesh: &Mesh, geo: &Geometry) -> Vec<bool> {
    let mut bad = vec![false; mesh.n_cells];
    for f in 0..mesh.n_faces() {
        if !face_ok(mesh, geo, f) {
            bad[mesh.owner[f] as usize] = true;
            if f < mesh.n_interior_faces() {
                bad[mesh.neighbour[f] as usize] = true;
            }
        }
    }
    bad
}

fn bad_faces(mesh: &Mesh, geo: &Geometry) -> usize {
    (0..mesh.n_faces())
        .filter(|&f| !face_ok(mesh, geo, f))
        .count()
}

fn face_ok(mesh: &Mesh, geo: &Geometry, f: usize) -> bool {
    let p = mesh.owner[f] as usize;
    let to = if f < mesh.n_interior_faces() {
        geo.cell_centroid[mesh.neighbour[f] as usize]
    } else {
        geo.face_centroid[f]
    };
    dot(sub(to, geo.cell_centroid[p]), geo.face_area[f]) > 0.0
}

/// Puts a point of the surface into the chords of a clipped face ring
/// ([`CutOptions::face_point`]) where the face stays a simple polygon,
/// oriented as before, with the point inside the uncut face `nodes`; a
/// chord whose point fails that stays straight. The ring is computed once
/// per face, so both cells of the face see the same points.
fn bend_chords(
    ring: &mut Vec<(u32, bool)>,
    nodes: &[u32],
    points: &mut Vec<[f64; 3]>,
    body: &dyn ImplicitBody,
) {
    let (area, _) = polygon_area_centroid(points, nodes);
    let an = norm(area);
    if an == 0.0 || !an.is_finite() {
        return;
    }
    let n = scale(area, 1.0 / an);
    let u = {
        let t = if n[0].abs() < 0.6 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let w = crate::geometry::cross(t, n);
        scale(w, 1.0 / norm(w))
    };
    let v = crate::geometry::cross(n, u);
    let to2 = |x: [f64; 3]| [dot(x, u), dot(x, v)];
    let poly: Vec<[f64; 2]> = nodes.iter().map(|&i| to2(points[i as usize])).collect();
    let mut i = 0;
    while i < ring.len() {
        if !ring[i].1 {
            i += 1;
            continue;
        }
        let k = ring.len();
        let (a, b) = (ring[i].0, ring[(i + 1) % k].0);
        let Some(m) = body.face_point(points[a as usize], points[b as usize], n) else {
            i += 1;
            continue;
        };
        let m2 = to2(m);
        let mut cand: Vec<[f64; 2]> = ring.iter().map(|e| to2(points[e.0 as usize])).collect();
        cand.insert(i + 1, m2);
        if m.iter().all(|x| x.is_finite())
            && inside_2d(&poly, m2)
            && simple_at(&cand, i + 1)
            && area_2d(&cand) > 0.0
        {
            points.push(m);
            ring.insert(i + 1, ((points.len() - 1) as u32, true));
            i += 2;
        } else {
            i += 1;
        }
    }
}

/// Whether `x` lies inside cell `c` of the uncut mesh: behind the plane of
/// every face (exact for convex cells, the usual case of the cut cells).
fn inside_cell(mesh: &Mesh, geo: &Geometry, c: usize, x: [f64; 3]) -> bool {
    mesh.cell_face_iter(c).all(|(f, is_nb)| {
        let s = dot(sub(x, geo.face_centroid[f]), geo.face_area[f]);
        if is_nb { s >= 0.0 } else { s <= 0.0 }
    })
}

/// Twice the signed area of a 2D polygon.
fn area_2d(p: &[[f64; 2]]) -> f64 {
    let k = p.len();
    (0..k)
        .map(|i| {
            let (a, b) = (p[i], p[(i + 1) % k]);
            a[0] * b[1] - a[1] * b[0]
        })
        .sum()
}

/// Whether `x` lies strictly inside the polygon `p` (crossing number).
fn inside_2d(p: &[[f64; 2]], x: [f64; 2]) -> bool {
    let k = p.len();
    let mut inside = false;
    for i in 0..k {
        let (a, b) = (p[i], p[(i + 1) % k]);
        // on an edge: not strictly inside
        let cr = (b[0] - a[0]) * (x[1] - a[1]) - (b[1] - a[1]) * (x[0] - a[0]);
        let within = (x[0] - a[0]) * (x[0] - b[0]) + (x[1] - a[1]) * (x[1] - b[1]) <= 0.0;
        if cr == 0.0 && within {
            return false;
        }
        if (a[1] > x[1]) != (b[1] > x[1]) {
            let t = (x[1] - a[1]) / (b[1] - a[1]);
            if x[0] < a[0] + t * (b[0] - a[0]) {
                inside = !inside;
            }
        }
    }
    inside
}

/// Whether the two edges at vertex `j` of the polygon `p` meet no other
/// edge (touching counts as meeting).
fn simple_at(p: &[[f64; 2]], j: usize) -> bool {
    let k = p.len();
    let orient = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    };
    let on = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        orient(a, b, c) == 0.0
            && c[0] >= a[0].min(b[0])
            && c[0] <= a[0].max(b[0])
            && c[1] >= a[1].min(b[1])
            && c[1] <= a[1].max(b[1])
    };
    let meet = |a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]| {
        let (o1, o2) = (orient(a, b, c), orient(a, b, d));
        let (o3, o4) = (orient(c, d, a), orient(c, d, b));
        ((o1 > 0.0) != (o2 > 0.0) && o1 != 0.0 && o2 != 0.0)
            && ((o3 > 0.0) != (o4 > 0.0) && o3 != 0.0 && o4 != 0.0)
            || on(a, b, c)
            || on(a, b, d)
            || on(c, d, a)
            || on(c, d, b)
    };
    let jm = (j + k - 1) % k;
    for (s0, s1) in [(jm, j), (j, (j + 1) % k)] {
        for e in 0..k {
            let e1 = (e + 1) % k;
            // adjacent edges share a vertex
            if e == s0 || e == s1 || e1 == s0 || e1 == s1 {
                continue;
            }
            if meet(p[s0], p[s1], p[e], p[e1]) {
                return false;
            }
        }
    }
    true
}

/// Removes ring entries whose point repeats the next one (zero-length ring
/// edges), keeping the chord flag of the edge that leaves the point.
fn dedupe(ring: &mut Vec<(u32, bool)>) {
    loop {
        let k = ring.len();
        if k < 2 {
            return;
        }
        let Some(i) = (0..k).find(|&i| ring[i].0 == ring[(i + 1) % k].0) else {
            return;
        };
        ring.remove(i);
    }
}

/// Closed loops from directed edges whose points have as many edges in as
/// out; the first unused edge starts a loop, at a point with several edges
/// out the first given is taken.
fn chain(edges: &[(u32, u32)]) -> Vec<Vec<u32>> {
    let mut used = vec![false; edges.len()];
    let mut out = Vec::new();
    for s in 0..edges.len() {
        if used[s] {
            continue;
        }
        used[s] = true;
        let start = edges[s].0;
        let mut l = vec![start];
        let mut at = edges[s].1;
        while at != start {
            l.push(at);
            match (0..edges.len()).find(|&e| !used[e] && edges[e].0 == at) {
                Some(e) => {
                    used[e] = true;
                    at = edges[e].1;
                }
                None => break, // open chain: degenerate input, kept as is
            }
        }
        out.push(l);
    }
    out
}

fn reversed(r: &[u32]) -> Vec<u32> {
    let mut v = Vec::with_capacity(r.len());
    v.push(r[0]);
    v.extend(r[1..].iter().rev());
    v
}

/// Joins rings of the same orientation along their common edges into one
/// ring, `None` if what is left is not a single loop.
fn join(rings: &[Vec<u32>]) -> Option<Vec<u32>> {
    let mut edges: Vec<(u32, u32)> = Vec::new();
    for r in rings {
        let k = r.len();
        for i in 0..k {
            edges.push((r[i], r[(i + 1) % k]));
        }
    }
    let mut keep = vec![true; edges.len()];
    for i in 0..edges.len() {
        if !keep[i] {
            continue;
        }
        let (u, v) = edges[i];
        if let Some(j) = (0..edges.len()).find(|&j| keep[j] && j != i && edges[j] == (v, u)) {
            keep[i] = false;
            keep[j] = false;
        }
    }
    let rest: Vec<(u32, u32)> = edges
        .iter()
        .zip(&keep)
        .filter(|x| *x.1)
        .map(|x| *x.0)
        .collect();
    if rest.len() < 3 {
        return None;
    }
    let mut next: HashMap<u32, u32> = HashMap::new();
    for &(u, v) in &rest {
        if next.insert(u, v).is_some() {
            return None;
        }
    }
    let start = rest[0].0;
    let mut l = vec![start];
    let mut at = next[&start];
    while at != start {
        if l.len() > rest.len() {
            return None;
        }
        l.push(at);
        at = *next.get(&at)?;
    }
    (l.len() == rest.len()).then_some(l)
}

/// Faces that a merged cell shares with one neighbour are joined into one
/// face (two cells share at most one face); the fan of a joined face that
/// folds (an L-shaped group of merged cells around an uncut neighbour)
/// cuts through the neighbour and moves volume between the two cells. A
/// join may move at most this fraction of the smaller cell's volume.
const FOLD: f64 = 1e-2;

/// The groups (other than those in `skip`) with which the cells `un` share
/// faces whose join would move more than [`FOLD`] of the smaller volume
/// (the volume between the joined face's fan about its node average and the
/// pieces, ⅓ |Σᵢ (cᵢ − m) · Sᵢ|), and whether all shared faces join into
/// one ring at all.
#[allow(clippy::too_many_arguments)]
fn folds(
    points: &[[f64; 3]],
    un: &[usize],
    cell_faces: &[Vec<usize>],
    faces: &[Face],
    uf: &mut Uf,
    skip: &[usize],
    group_volume: &dyn Fn(usize) -> f64,
    un_volume: f64,
) -> (Vec<usize>, bool) {
    let mut by: Vec<(usize, Vec<Vec<u32>>)> = Vec::new();
    for &m in un {
        for &fi in &cell_faces[m] {
            let f = &faces[fi];
            if f.b < 0 {
                continue;
            }
            let (o, ring) = if f.a as usize == m {
                (f.b as usize, f.ring.clone())
            } else {
                (f.a as usize, reversed(&f.ring))
            };
            let og = uf.find(o);
            if skip.contains(&og) {
                continue;
            }
            match by.iter_mut().find(|x| x.0 == og) {
                Some(x) => x.1.push(ring),
                None => by.push((og, vec![ring])),
            }
        }
    }
    let mut joins = true;
    let mut out = Vec::new();
    for (og, r) in by {
        if r.len() == 1 {
            continue;
        }
        let Some(joined) = join(&r) else {
            joins = false;
            out.push(og);
            continue;
        };
        if fold_fraction(points, &r, &joined, group_volume(og).min(un_volume)) > FOLD {
            out.push(og);
        }
    }
    (out, joins)
}

/// The smallest cosine between the centroid distance d and the outward
/// area vector S over the faces of the union of the cells `un` (from the
/// union's centroid to the neighbour group's centroid, or to the face
/// centroid at a boundary): positive when the union is a cell the
/// discretisation takes.
#[allow(clippy::too_many_arguments)]
fn union_score(
    un: &[usize],
    cell_faces: &[Vec<usize>],
    faces: &[Face],
    face_geo: &[([f64; 3], [f64; 3])],
    geo: &Geometry,
    members: &[Vec<usize>],
    uf: &mut Uf,
) -> f64 {
    let centroid = |cells: &[usize]| -> [f64; 3] {
        let mut v = 0.0;
        let mut m = [0.0; 3];
        for &c in cells {
            v += geo.cell_volume[c];
            m = add(m, scale(geo.cell_centroid[c], geo.cell_volume[c]));
        }
        scale(m, 1.0 / v)
    };
    let x = centroid(un);
    let mut worst = f64::INFINITY;
    for &m in un {
        for &fi in &cell_faces[m] {
            let f = &faces[fi];
            let (s0, c0) = face_geo[fi];
            let (s, other) = if f.a as usize == m {
                (s0, f.b)
            } else {
                (scale(s0, -1.0), f.a as i64)
            };
            let to = if other >= 0 {
                let o = other as usize;
                if un.contains(&o) {
                    continue;
                }
                let og = uf.find(o);
                centroid(&members[og])
            } else {
                c0
            };
            let d = sub(to, x);
            let den = norm(d) * norm(s);
            let cos = if den > 0.0 { dot(d, s) / den } else { -1.0 };
            worst = worst.min(cos);
        }
    }
    worst
}

/// The volume a joined face moves between its two cells, over `volume`:
/// ⅓ |Σᵢ (cᵢ − m) · Sᵢ| over the pieces `rings`, m the node average of
/// the joined ring (its fan carries no flux of (x − m)).
fn fold_fraction(points: &[[f64; 3]], rings: &[Vec<u32>], joined: &[u32], volume: f64) -> f64 {
    let m = crate::geometry::node_average(points, joined);
    let mut dv = 0.0;
    for ring in rings {
        let (a, c) = polygon_area_centroid(points, ring);
        dv += dot(sub(c, m), a);
    }
    (dv / 3.0).abs() / volume.max(f64::MIN_POSITIVE)
}

/// Union-find over cells (the smaller root wins: deterministic).
struct Uf {
    p: Vec<usize>,
}

impl Uf {
    fn new(n: usize) -> Self {
        Self {
            p: (0..n).collect(),
        }
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.p[x] != x {
            self.p[x] = self.p[self.p[x]];
            x = self.p[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) -> usize {
        let (a, b) = (self.find(a), self.find(b));
        let r = a.min(b);
        self.p[a] = r;
        self.p[b] = r;
        r
    }
}

/// Union-find over the few points of one cell.
#[derive(Default)]
struct LocalUf {
    p: HashMap<u32, u32>,
}

impl LocalUf {
    fn find(&mut self, x: u32) -> u32 {
        let mut r = x;
        loop {
            let q = *self.p.entry(r).or_insert(r);
            if q == r {
                break;
            }
            r = q;
        }
        // path compression
        let mut y = x;
        while y != r {
            let q = self.p[&y];
            self.p.insert(y, r);
            y = q;
        }
        r
    }
    fn union(&mut self, a: u32, b: u32) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            let (lo, hi) = (a.min(b), a.max(b));
            self.p.insert(hi, lo);
        }
    }
}

/// The analytic bodies of the tests and the benchmark: a sphere.
#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    /// Centre.
    pub centre: [f64; 3],
    /// Radius.
    pub radius: f64,
    /// Put crossing points on the true surface (else linear).
    pub exact: bool,
}

impl ImplicitBody for Sphere {
    fn value(&self, x: [f64; 3]) -> f64 {
        norm(sub(x, self.centre)) - self.radius
    }
    fn edge_point(&self, a: [f64; 3], b: [f64; 3], fa: f64, fb: f64) -> [f64; 3] {
        if self.exact {
            edge_root(|x| self.value(x), a, b, fa, fb)
        } else {
            let t = fa / (fa - fb);
            add(a, scale(sub(b, a), t))
        }
    }
    fn project(&self, x: [f64; 3]) -> Option<[f64; 3]> {
        let d = sub(x, self.centre);
        let n = norm(d);
        (n > 0.0).then(|| add(self.centre, scale(d, self.radius / n)))
    }
    fn face_point(&self, p: [f64; 3], q: [f64; 3], n: [f64; 3]) -> Option<[f64; 3]> {
        let normal = |x: [f64; 3]| {
            let d = sub(x, self.centre);
            let l = norm(d);
            (l > 0.0).then(|| scale(d, 1.0 / l))
        };
        chord_point(|x| self.value(x), normal, p, q, n)
    }
}

/// The faces, parents and cut flags of the merged cells, whether the mesh
/// before the pass had ill-shaped cells, and the merged cell of every cell.
type Merged = (Vec<Face>, Vec<u32>, Vec<bool>, bool, Vec<u32>);

/// One pass of the small-cell merge (see the module docs): the faces,
/// parents and cut flags of the merged cells, and whether the mesh before
/// the pass had ill-shaped cells.
#[allow(clippy::too_many_arguments)]
fn merge_pass(
    mesh: &Mesh,
    points: &[[f64; 3]],
    out_faces: &[Face],
    parent: &[u32],
    cut_cell: &[bool],
    geo0: &Geometry,
    bnames: &[String],
    min_fraction: f64,
    merged: &mut usize,
    failed: &mut usize,
    folded: &mut usize,
    max_fold: &mut f64,
) -> Result<Merged, String> {
    let n_cells = parent.len();
    let zone: Vec<u32> = parent.iter().map(|&p| mesh.cell_zone[p as usize]).collect();
    let tmp = build(points, n_cells, out_faces, bnames, &mesh.zones, &zone)?;
    let geo1 = Geometry::compute(&tmp);
    let frac: Vec<f64> = (0..n_cells)
        .map(|c| geo1.cell_volume[c] / geo0.cell_volume[parent[c] as usize])
        .collect();
    let bad1 = bad_cells(&tmp, &geo1);
    let mut small: Vec<usize> = (0..n_cells)
        .filter(|&c| cut_cell[c] && (frac[c] < min_fraction || bad1[c]))
        .collect();
    small.sort_by(|&x, &y| frac[x].total_cmp(&frac[y]).then(x.cmp(&y)));
    let is_small: Vec<bool> = {
        let mut v = vec![false; n_cells];
        for &c in &small {
            v[c] = true;
        }
        v
    };
    // faces of each cell in the cut mesh (index into out_faces)
    let mut cell_faces: Vec<Vec<usize>> = vec![Vec::new(); n_cells];
    for (i, f) in out_faces.iter().enumerate() {
        cell_faces[f.a as usize].push(i);
        if f.b >= 0 {
            cell_faces[f.b as usize].push(i);
        }
    }
    let face_geo: Vec<([f64; 3], [f64; 3])> = out_faces
        .iter()
        .map(|f| polygon_area_centroid(points, &f.ring))
        .collect();
    let face_area: Vec<f64> = face_geo.iter().map(|g| norm(g.0)).collect();
    let mut uf = Uf::new(n_cells);
    let mut members: Vec<Vec<usize>> = (0..n_cells).map(|c| vec![c]).collect();
    for &s in &small {
        let g = uf.find(s);
        if g != s && members[g].len() > 1 && members[g].iter().any(|&m| !is_small[m]) {
            continue; // already part of a large enough cell
        }
        // candidate groups by shared area
        let mut cand: Vec<(usize, f64)> = Vec::new();
        for &m in &members[g] {
            for &fi in &cell_faces[m] {
                let f = &out_faces[fi];
                if f.b < 0 {
                    continue;
                }
                let o = if f.a as usize == m {
                    f.b as usize
                } else {
                    f.a as usize
                };
                let og = uf.find(o);
                if og == g {
                    continue;
                }
                match cand.iter_mut().find(|x| x.0 == og) {
                    Some(x) => x.1 += face_area[fi],
                    None => cand.push((og, face_area[fi])),
                }
            }
        }
        let small_group = |h: usize| members[h].iter().all(|&m| is_small[m]);
        cand.sort_by(|x, y| {
            small_group(x.0)
                .cmp(&small_group(y.0))
                .then(y.1.total_cmp(&x.1))
                .then(x.0.cmp(&y.0))
        });
        // the first neighbour (largest shared face first) whose shared
        // faces join into one ring and whose union with the cell is valid
        // (every face of the union with d · S > 0 from its centroid); if
        // none is, the one whose union comes closest (a join may fold:
        // `folded_faces` counts those)
        let gv = |x: usize, members: &[Vec<usize>]| -> f64 {
            members[x].iter().map(|&m| geo1.cell_volume[m]).sum()
        };
        let mut pick: Option<Vec<usize>> = None;
        let mut valid: Option<Vec<usize>> = None;
        let mut fallback: Option<(f64, Vec<usize>)> = None;
        for &(h, _) in &cand {
            let mut un: Vec<usize> = members[g].clone();
            un.extend_from_slice(&members[h]);
            let gvm = |x: usize| gv(x, &members);
            let v_un = gvm(g) + gvm(h);
            let (fold, joins) = folds(
                points,
                &un,
                &cell_faces,
                out_faces,
                &mut uf,
                &[g, h],
                &gvm,
                v_un,
            );
            if !joins {
                continue;
            }
            let score = union_score(
                &un,
                &cell_faces,
                out_faces,
                &face_geo,
                &geo1,
                &members,
                &mut uf,
            );
            if score > 0.0 {
                if fold.is_empty() {
                    pick = Some(vec![g, h]);
                    break;
                }
                if valid.is_none() {
                    valid = Some(vec![g, h]);
                }
            } else if fallback.as_ref().is_none_or(|f| score > f.0) {
                fallback = Some((score, vec![g, h]));
            }
        }
        let pick = pick.or(valid);
        let take = pick.or(fallback.map(|f| f.1));
        let mut done = false;
        if let Some(take) = take {
            let mut root = take[0];
            for &x in &take[1..] {
                root = uf.union(root, x);
            }
            for &x in &take {
                if x != root {
                    let moved = std::mem::take(&mut members[x]);
                    members[root].extend(moved);
                }
            }
            members[root].sort_unstable();
            done = true;
        }
        if done {
            *merged += 1;
        } else {
            *failed += 1;
        }
    }

    // 4. the final mesh: one cell per group, faces between groups joined
    let mut gid = vec![u32::MAX; n_cells];
    let mut final_parent = Vec::new();
    for c in 0..n_cells {
        let g = uf.find(c);
        if gid[g] == u32::MAX {
            gid[g] = final_parent.len() as u32;
            // the parent of a group is that of its largest member
            let big = *members[g]
                .iter()
                .max_by(|&&x, &&y| {
                    geo1.cell_volume[x]
                        .total_cmp(&geo1.cell_volume[y])
                        .then(y.cmp(&x))
                })
                .expect("member");
            final_parent.push(parent[big]);
        }
        gid[c] = gid[g];
    }
    let nfin = final_parent.len();
    let mut pair_faces: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    let mut order: Vec<(u32, u32)> = Vec::new();
    let mut fin_faces: Vec<Face> = Vec::new();
    for (i, f) in out_faces.iter().enumerate() {
        let a = gid[f.a as usize];
        if f.b < 0 {
            fin_faces.push(Face {
                ring: f.ring.clone(),
                a,
                b: f.b,
            });
            continue;
        }
        let b = gid[f.b as usize];
        if a == b {
            continue;
        }
        let key = (a.min(b), a.max(b));
        let e = pair_faces.entry(key).or_insert_with(|| {
            order.push(key);
            Vec::new()
        });
        e.push(i);
    }
    for key in order {
        let list = &pair_faces[&key];
        // rings outward from key.0
        let rings: Vec<Vec<u32>> = list
            .iter()
            .map(|&i| {
                let f = &out_faces[i];
                if gid[f.a as usize] == key.0 {
                    f.ring.clone()
                } else {
                    reversed(&f.ring)
                }
            })
            .collect();
        let ring = if rings.len() == 1 {
            rings.into_iter().next().expect("one")
        } else {
            let joined = join(&rings).ok_or("internal: merged faces do not form one ring")?;
            let vol = |g: u32| -> f64 {
                (0..n_cells)
                    .filter(|&c| gid[c] == g)
                    .map(|c| geo1.cell_volume[c])
                    .sum()
            };
            let fr = fold_fraction(points, &rings, &joined, vol(key.0).min(vol(key.1)));
            if fr > FOLD {
                *folded += 1;
            }
            *max_fold = max_fold.max(fr);
            joined
        };
        fin_faces.push(Face {
            ring,
            a: key.0,
            b: key.1 as i64,
        });
    }
    let mut fin_cut = vec![false; nfin];
    for c in 0..n_cells {
        if cut_cell[c] {
            fin_cut[gid[c] as usize] = true;
        }
    }
    Ok((
        fin_faces,
        final_parent,
        fin_cut,
        bad1.iter().any(|&b| b),
        gid,
    ))
}
