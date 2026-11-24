//! Datums and their reports (SP1b, D262 to D267): named surfaces and
//! volumes of a case, the gas flowing through them and the liquid of the
//! sprays crossing them.
//!
//! ```toml
//! [[datum]]
//! name = "section"            # column prefix in datums.csv
//! type = "plane"              # plane, boundary, interface, surface, box, cylinder, zone
//! point = [1.0, 0.2, 0.05]
//! normal = [1.0, 0.0, 0.0]    # the positive direction
//! # radius = 0.1              # bound: a disc ...
//! # up = [0, 1, 0]            # ... or a rectangle: height along `up`,
//! # width = 0.2               #     width along up × normal
//! # height = 0.1
//! # snap = "faces"            # exact: the mesh faces nearest to the plane
//! # histogram = { bins = 40, min = 1e-6, max = 1e-3, scale = "log" }
//! ```
//!
//! Kinds (D262):
//!
//! - `plane` (`point`, `normal`; optional bound by `radius`, or by `up`,
//!   `width`, `height`): with the default `snap = "cut"` the plane is cut
//!   through the cells, each cut polygon carrying the cell value and its
//!   least-squares gradient (D24) at its centroid; `cut_error` reports the
//!   relative difference of its mass flow against the exact flow through
//!   the faces of the same plane snapped to the mesh. With `snap = "faces"`
//!   the plane is the layer of mesh faces between the cells whose centroids
//!   lie on its two sides: exact and conservative.
//! - `boundary` (`boundary`): the faces of a boundary; positive out of the
//!   domain. Exact.
//! - `interface` (`zones = [a, b]`): the faces between two cell zones;
//!   positive from `a` to `b`. Exact.
//! - `surface` (`file`, STL; `scale`, default 1): triangles; gas values from
//!   the cell that holds each triangle's centroid (value and gradient);
//!   liquid by the exact segment-triangle test. Positive along the
//!   triangle normals (right-hand rule of the corners).
//! - `box` (`lower`, `upper`), `cylinder` (`start`, `end`, `radius`),
//!   `zone` (`zone`): volumes; the gas averages over the cells whose
//!   centroids lie inside, the parcels inside by their positions.
//!
//! Exact surfaces sum the face mass fluxes of the solver; every sum runs in
//! a fixed order, so reports are bitwise equal for any thread count.

use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use crate::simplec::Solution;
use anyhow::Result;
use foehncraft_mesh::Mesh;
use foehncraft_mesh::geometry::Geometry;
use foehncraft_particles::datum::{Bound, Datum, HistSpec, Shape, Surface, Volume};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

/// The kinds of `[[datum]] type`.
pub const KINDS: &[&str] = &[
    "plane",
    "boundary",
    "interface",
    "surface",
    "box",
    "cylinder",
    "zone",
];

/// `histogram = { … }` of a datum: the droplet size bins (D264).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct HistogramKeys {
    /// Number of bins (default 40).
    pub bins: usize,
    /// Lower edge (m, default 1e-6).
    pub min: f64,
    /// Upper edge (m, default 1e-3).
    pub max: f64,
    /// `log` (default) or `linear`.
    pub scale: String,
}

impl Default for HistogramKeys {
    fn default() -> Self {
        Self {
            bins: 40,
            min: 1e-6,
            max: 1e-3,
            scale: "log".into(),
        }
    }
}

/// A `[[datum]]` table (D262).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct DatumSection {
    /// Name: the column prefix in `datums.csv` (letters, digits, `_`, `-`).
    pub name: String,
    /// `plane`, `boundary`, `interface`, `surface`, `box`, `cylinder` or
    /// `zone`.
    #[serde(rename = "type")]
    pub kind: String,
    /// Plane: a point of it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub point: Option<[f64; 3]>,
    /// Plane: the normal (any length), the positive direction.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub normal: Option<[f64; 3]>,
    /// Plane: a disc of this radius about `point`; cylinder: its radius.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub radius: Option<f64>,
    /// Plane: a rectangle about `point`, its height along `up` (projected
    /// into the plane) …
    #[serde(skip_serializing_if = "Option::is_none")]
    pub up: Option<[f64; 3]>,
    /// … its width along up × normal (m).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub width: Option<f64>,
    /// … and its height (m).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub height: Option<f64>,
    /// Plane: `cut` (default) through the cells, or `faces`, moved to the
    /// mesh faces nearest to it (exact).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub snap: Option<String>,
    /// Boundary: its name.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub boundary: Option<String>,
    /// Interface: the two cell zones; positive from the first to the second.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub zones: Option<[String; 2]>,
    /// Surface: the STL file, relative to the case file.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub file: Option<PathBuf>,
    /// Surface: factor on the coordinates of the file (default 1).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scale: Option<f64>,
    /// Box: lower corner.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lower: Option<[f64; 3]>,
    /// Box: upper corner.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub upper: Option<[f64; 3]>,
    /// Cylinder: centre of one end.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub start: Option<[f64; 3]>,
    /// Cylinder: centre of the other end.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end: Option<[f64; 3]>,
    /// Zone: the cell zone.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub zone: Option<String>,
    /// Droplet size histogram (sprays).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub histogram: Option<HistogramKeys>,
}

impl DatumSection {
    /// The histogram bins.
    pub fn hist(&self) -> HistSpec {
        let h = self.histogram.clone().unwrap_or_default();
        HistSpec {
            bins: h.bins,
            min: h.min,
            max: h.max,
            log: h.scale != "linear",
        }
    }
}

fn err(key: &str, message: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::DATUM_CONFIG, message).key(key))
}

fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    (n.is_finite() && n > 0.0).then(|| [a[0] / n, a[1] / n, a[2] / n])
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Two unit vectors e1, e2 in the plane of `n` with e1 × e2 = n.
fn frame(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let t = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let e1 = unit(cross(t, n)).unwrap_or([1.0, 0.0, 0.0]);
    let e2 = cross(n, e1);
    (e1, e2)
}

/// Resolves the paths of the datums against the case directory.
pub fn resolve_paths(list: &mut [DatumSection], dir: &Path) {
    for d in list {
        if let Some(f) = &mut d.file {
            *f = dir.join(&*f);
        }
    }
}

/// Checks the tables that need no mesh (DAT001).
pub fn validate(case: &FlowCase) -> Result<()> {
    if !case.datums.is_empty() && !case.sliding.is_empty() {
        return Err(err(
            "datum",
            "[[datum]] with a sliding mesh is not available in this version (the faces move)"
                .into(),
        ));
    }
    let mut names = std::collections::BTreeSet::new();
    for (i, d) in case.datums.iter().enumerate() {
        let key = format!("datum[{i}]");
        if d.name.is_empty()
            || !d
                .name
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
        {
            return Err(err(
                &key,
                format!(
                    "[[datum]] name `{}`: use letters, digits, `_` and `-`",
                    d.name
                ),
            ));
        }
        if !names.insert(d.name.clone()) {
            return Err(err(
                &key,
                format!("[[datum]] {}: the name is used twice", d.name),
            ));
        }
        if !KINDS.contains(&d.kind.as_str()) {
            return Err(err(
                &key,
                format!(
                    "[[datum]] {}: unknown type `{}`; known: {}",
                    d.name,
                    d.kind,
                    KINDS.join(", ")
                ),
            ));
        }
        let h = d.hist();
        if h.bins == 0
            || !(h.min.is_finite() && h.max.is_finite() && h.max > h.min)
            || (h.log && h.min <= 0.0)
        {
            return Err(err(
                &key,
                format!(
                    "[[datum]] {}: histogram needs bins > 0 and 0 < min < max (log bins)",
                    d.name
                ),
            ));
        }
        if let Some(s) = &d.histogram
            && s.scale != "log"
            && s.scale != "linear"
        {
            return Err(err(
                &key,
                format!("[[datum]] {}: histogram scale is `log` or `linear`", d.name),
            ));
        }
    }
    Ok(())
}

/// What the gas report of a datum sums over.
#[derive(Debug, Clone)]
pub enum Facets {
    /// Mesh faces with a sign (+1: the face's mass flux, out of its owner,
    /// is along the datum's direction).
    Faces(Vec<(usize, f64)>),
    /// Pieces of a surface in cells: cell, area vector (along the positive
    /// direction), centroid.
    Cut(Vec<(usize, [f64; 3], [f64; 3])>),
    /// Cells of a volume.
    Cells(Vec<usize>),
}

/// A datum ready for reports.
#[derive(Debug, Clone)]
pub struct Built {
    /// The table.
    pub section: DatumSection,
    /// The gas side.
    pub facets: Facets,
    /// The faces of the same plane snapped to the mesh (cut planes only),
    /// for `cut_error`.
    pub snapped: Option<Vec<(usize, f64)>>,
    /// The parcel side.
    pub parcels: Datum,
}

impl Built {
    /// True for volumes.
    pub fn is_volume(&self) -> bool {
        matches!(self.facets, Facets::Cells(_))
    }
}

/// Point, unit normal, bound, and the in-plane axes e1, e2 of the bound.
type PlaneParts = ([f64; 3], [f64; 3], Bound, [f64; 3], [f64; 3]);

/// The plane of a `plane` datum: point, unit normal, bound, e1, e2 of the
/// bound.
fn plane_of(d: &DatumSection, key: &str) -> Result<PlaneParts> {
    let point = d
        .point
        .ok_or_else(|| err(key, format!("[[datum]] {}: a plane needs `point`", d.name)))?;
    let normal = d.normal.and_then(unit).ok_or_else(|| {
        err(
            key,
            format!("[[datum]] {}: a plane needs a nonzero `normal`", d.name),
        )
    })?;
    let pos = |x: Option<f64>| x.is_some_and(|v| v.is_finite() && v > 0.0);
    match (d.radius, d.up, d.width, d.height) {
        (None, None, None, None) => {
            let (e1, e2) = frame(normal);
            Ok((point, normal, Bound::None, e1, e2))
        }
        (Some(r), None, None, None) if pos(Some(r)) => {
            let (e1, e2) = frame(normal);
            Ok((point, normal, Bound::Disc { radius: r }, e1, e2))
        }
        (None, Some(up), Some(w), Some(h)) if pos(Some(w)) && pos(Some(h)) => {
            let upn = dot(up, normal);
            let e2 = unit([
                up[0] - upn * normal[0],
                up[1] - upn * normal[1],
                up[2] - upn * normal[2],
            ])
            .ok_or_else(|| {
                err(
                    key,
                    format!("[[datum]] {}: `up` is parallel to the normal", d.name),
                )
            })?;
            let e1 = cross(e2, normal);
            Ok((
                point,
                normal,
                Bound::Rectangle {
                    e1,
                    e2,
                    half_width: 0.5 * w,
                    half_height: 0.5 * h,
                },
                e1,
                e2,
            ))
        }
        _ => Err(err(
            key,
            format!(
                "[[datum]] {}: bound a plane by `radius` (> 0) or by `up`, `width` and `height` (> 0), not both",
                d.name
            ),
        )),
    }
}

/// The convex polygon of a bound in (e1, e2) coordinates, counter-clockwise;
/// a disc is the regular 128-gon of the same area.
fn bound_polygon(b: &Bound) -> Option<Vec<[f64; 2]>> {
    match *b {
        Bound::None => None,
        Bound::Disc { radius } => {
            let n = 128;
            let a = 2.0 * std::f64::consts::PI / n as f64;
            // n/2 R² sin(2π/n) = π r²
            let r = radius * (std::f64::consts::PI / (0.5 * n as f64 * a.sin())).sqrt();
            Some(
                (0..n)
                    .map(|i| {
                        let t = a * i as f64;
                        [r * t.cos(), r * t.sin()]
                    })
                    .collect(),
            )
        }
        Bound::Rectangle {
            half_width: w,
            half_height: h,
            ..
        } => Some(vec![[-w, -h], [w, -h], [w, h], [-w, h]]),
    }
}

/// Sutherland-Hodgman: `poly` clipped by the convex counter-clockwise
/// polygon `clip` (the orientation of `poly` is kept).
fn clip2(poly: &[[f64; 2]], clip: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut out = poly.to_vec();
    for i in 0..clip.len() {
        if out.is_empty() {
            break;
        }
        let a = clip[i];
        let b = clip[(i + 1) % clip.len()];
        let side = |p: [f64; 2]| (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
        let inp = std::mem::take(&mut out);
        for j in 0..inp.len() {
            let p = inp[j];
            let q = inp[(j + 1) % inp.len()];
            let (sp, sq) = (side(p), side(q));
            if sp >= 0.0 {
                out.push(p);
            }
            if (sp >= 0.0) != (sq >= 0.0) {
                let t = sp / (sp - sq);
                out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
            }
        }
    }
    out
}

/// Signed area and first moments of a 2D polygon.
fn area2(p: &[[f64; 2]]) -> (f64, [f64; 2]) {
    let (mut a, mut mx, mut my) = (0.0, 0.0, 0.0);
    for i in 0..p.len() {
        let (u, v) = (p[i], p[(i + 1) % p.len()]);
        let c = u[0] * v[1] - v[0] * u[1];
        a += c;
        mx += (u[0] + v[0]) * c;
        my += (u[1] + v[1]) * c;
    }
    (0.5 * a, [mx / 6.0, my / 6.0])
}

/// The cut of cell `c` by the plane (point `p`, unit normal `n`, in-plane
/// frame e1, e2 with e1 × e2 = n), clipped to `clip`: its area (along `n`)
/// and centroid, or `None` when the plane misses the cell.
#[allow(clippy::too_many_arguments)]
pub fn cut_cell(
    mesh: &Mesh,
    geo: &Geometry,
    c: usize,
    p: [f64; 3],
    n: [f64; 3],
    e1: [f64; 3],
    e2: [f64; 3],
    clip: Option<&[[f64; 2]]>,
) -> Option<(f64, [f64; 3])> {
    let to2 = |x: [f64; 3]| {
        let r = sub(x, p);
        [dot(r, e1), dot(r, e2)]
    };
    let xc = geo.cell_centroid[c];
    let q0 = to2(xc);
    let (mut area, mut m) = (0.0, [0.0; 2]);
    let mut any = false;
    for (f, neighbour) in mesh.cell_face_iter(c) {
        let nodes = mesh.face(f);
        let k = nodes.len();
        let pt = |i: usize| {
            let j = if neighbour { (k - i) % k } else { i };
            mesh.points[nodes[j] as usize]
        };
        let s: Vec<f64> = (0..k).map(|i| dot(sub(pt(i), p), n)).collect();
        let Some(i0) = (0..k).find(|&i| s[i] <= 0.0) else {
            continue;
        };
        if (0..k).all(|i| s[i] <= 0.0) {
            continue;
        }
        let at = |i: usize, j: usize| {
            let (a, b) = (pt(i), pt(j));
            let t = s[i] / (s[i] - s[j]);
            to2([
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ])
        };
        let mut exit: Option<[f64; 2]> = None;
        for step in 0..k {
            let i = (i0 + step) % k;
            let j = (i + 1) % k;
            if s[i] <= 0.0 && s[j] > 0.0 {
                exit = Some(at(i, j));
            } else if s[i] > 0.0 && s[j] <= 0.0 {
                let entry = at(i, j);
                if let Some(e) = exit.take() {
                    // the cap edge runs from the entry to the exit point
                    let tri = [q0, entry, e];
                    let piece = match clip {
                        Some(cl) => clip2(&tri, cl),
                        None => tri.to_vec(),
                    };
                    if piece.len() >= 3 {
                        let (a, mm) = area2(&piece);
                        area += a;
                        m[0] += mm[0];
                        m[1] += mm[1];
                        any = true;
                    }
                }
            }
        }
    }
    if !any || area.abs() <= 1e-300 {
        return None;
    }
    let (u, v) = (m[0] / area, m[1] / area);
    Some((
        area,
        [
            p[0] + u * e1[0] + v * e2[0],
            p[1] + u * e1[1] + v * e2[1],
            p[2] + u * e1[2] + v * e2[2],
        ],
    ))
}

/// The face layer of a plane: inner faces between cells whose centroids lie
/// on different sides (s < 0 against s ≥ 0), the face centroid inside the
/// bound; and the side of every cell.
fn snapped_faces(
    mesh: &Mesh,
    geo: &Geometry,
    p: [f64; 3],
    n: [f64; 3],
    bound: &Bound,
) -> (Vec<(usize, f64)>, Vec<i8>) {
    let side: Vec<i8> = (0..mesh.n_cells)
        .map(|c| {
            if dot(sub(geo.cell_centroid[c], p), n) < 0.0 {
                -1
            } else {
                1
            }
        })
        .collect();
    let mut faces = Vec::new();
    for f in 0..mesh.n_inner_faces {
        let (o, nb) = (mesh.owner[f] as usize, mesh.neighbour[f] as usize);
        if side[o] == side[nb] {
            continue;
        }
        let x = geo.face_centroid[f];
        let h = dot(sub(x, p), n);
        let xp = [x[0] - h * n[0], x[1] - h * n[1], x[2] - h * n[2]];
        if !bound.contains(p, xp) {
            continue;
        }
        faces.push((f, if side[o] < 0 { 1.0 } else { -1.0 }));
    }
    (faces, side)
}

/// The gas side of a surface triangle: pieces smaller than a quarter of the
/// size of the cell holding their centroid (midpoint subdivision, at most
/// ten levels; six for pieces whose centroid is outside the mesh, which are
/// dropped then), each with the cell of its centroid. A triangle much larger
/// than the cells is thus sampled at many points instead of one.
fn sample_triangle(
    loc: &foehncraft_mesh::locate::Locator,
    geo: &Geometry,
    t: [[f64; 3]; 3],
    out: &mut Vec<(usize, [f64; 3], [f64; 3])>,
) {
    let mut stack = vec![(t, 0u32)];
    while let Some((t, depth)) = stack.pop() {
        let x = [
            (t[0][0] + t[1][0] + t[2][0]) / 3.0,
            (t[0][1] + t[1][1] + t[2][1]) / 3.0,
            (t[0][2] + t[1][2] + t[2][2]) / 3.0,
        ];
        let edge = (0..3)
            .map(|i| {
                let d = sub(t[(i + 1) % 3], t[i]);
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .fold(0.0, f64::max);
        let cell = loc.locate_search(x);
        let h = cell.map_or(0.0, |c| geo.cell_volume[c].cbrt());
        if (cell.is_none() && depth < 6) || (depth < 10 && edge > 0.25 * h) {
            let m = |a: [f64; 3], b: [f64; 3]| {
                [
                    0.5 * (a[0] + b[0]),
                    0.5 * (a[1] + b[1]),
                    0.5 * (a[2] + b[2]),
                ]
            };
            let (ab, bc, ca) = (m(t[0], t[1]), m(t[1], t[2]), m(t[2], t[0]));
            for s in [[t[0], ab, ca], [ab, t[1], bc], [ca, bc, t[2]], [ab, bc, ca]] {
                stack.push((s, depth + 1));
            }
            continue;
        }
        if let Some(c) = cell {
            let a = cross(sub(t[1], t[0]), sub(t[2], t[0]));
            out.push((c, [0.5 * a[0], 0.5 * a[1], 0.5 * a[2]], x));
        }
    }
}

/// Builds the datums of the case on its mesh (DAT001, DAT002).
pub fn build(case: &FlowCase, mesh: &Mesh, geo: &Geometry) -> Result<Vec<Built>> {
    let mut out = Vec::new();
    for (i, d) in case.datums.iter().enumerate() {
        let key = format!("datum[{i}]");
        let name = &d.name;
        let need = |v: Option<[f64; 3]>, k: &str| {
            v.ok_or_else(|| err(&key, format!("[[datum]] {name}: a {} needs `{k}`", d.kind)))
        };
        let (facets, snapped, shape) = match d.kind.as_str() {
            "plane" => {
                let (p, n, bound, e1, e2) = plane_of(d, &key)?;
                let (faces, side) = snapped_faces(mesh, geo, p, n, &bound);
                match d.snap.as_deref().unwrap_or("cut") {
                    "cut" => {
                        let poly = bound_polygon(&bound);
                        let mut cut = Vec::new();
                        for c in 0..mesh.n_cells {
                            if let Some((a, x)) =
                                cut_cell(mesh, geo, c, p, n, e1, e2, poly.as_deref())
                            {
                                cut.push((c, [a * n[0], a * n[1], a * n[2]], x));
                            }
                        }
                        (
                            Facets::Cut(cut),
                            Some(faces),
                            Shape::Surface(Surface::Plane {
                                point: p,
                                normal: n,
                                bound,
                            }),
                        )
                    }
                    "faces" => (
                        Facets::Faces(faces),
                        None,
                        Shape::Surface(Surface::Cells {
                            side,
                            plane: Some((p, n, bound)),
                            faces: None,
                        }),
                    ),
                    s => {
                        return Err(err(
                            &key,
                            format!("[[datum]] {name}: snap `{s}`; known: cut, faces"),
                        ));
                    }
                }
            }
            "boundary" => {
                let bn = d.boundary.as_deref().ok_or_else(|| {
                    err(
                        &key,
                        format!("[[datum]] {name}: a boundary needs `boundary`"),
                    )
                })?;
                let Some(b) = mesh.boundaries.iter().position(|x| x.name == bn) else {
                    let all: Vec<&str> = mesh.boundaries.iter().map(|x| x.name.as_str()).collect();
                    return Err(err(
                        &key,
                        format!(
                            "[[datum]] {name}: the mesh has no boundary `{bn}`; it has {}",
                            all.join(", ")
                        ),
                    ));
                };
                (
                    Facets::Faces(mesh.boundaries[b].range().map(|f| (f, 1.0)).collect()),
                    None,
                    Shape::Surface(Surface::Boundary { boundary: b }),
                )
            }
            "interface" => {
                let [za, zb] = d.zones.clone().ok_or_else(|| {
                    err(
                        &key,
                        format!("[[datum]] {name}: an interface needs `zones = [a, b]`"),
                    )
                })?;
                let zi = |z: &str| {
                    mesh.zones.iter().position(|x| x == z).ok_or_else(|| {
                        err(
                            &key,
                            format!(
                                "[[datum]] {name}: the mesh has no cell zone `{z}`; it has {}",
                                mesh.zones.join(", ")
                            ),
                        )
                    })
                };
                let (a, b) = (zi(&za)? as u32, zi(&zb)? as u32);
                let side: Vec<i8> = (0..mesh.n_cells)
                    .map(|c| match mesh.cell_zone[c] {
                        z if z == a => -1,
                        z if z == b => 1,
                        _ => 0,
                    })
                    .collect();
                let mut faces = Vec::new();
                for f in 0..mesh.n_inner_faces {
                    let (o, nb) = (mesh.owner[f] as usize, mesh.neighbour[f] as usize);
                    match (side[o], side[nb]) {
                        (-1, 1) => faces.push((f, 1.0)),
                        (1, -1) => faces.push((f, -1.0)),
                        _ => {}
                    }
                }
                let mut sf = foehncraft_mesh::query::SurfaceFaces::default();
                for &(f, _) in &faces {
                    sf.push_face(mesh, geo, f);
                }
                let bvh =
                    (!faces.is_empty()).then(|| Box::new(foehncraft_mesh::query::Bvh::build(sf)));
                (
                    Facets::Faces(faces),
                    None,
                    Shape::Surface(Surface::Cells {
                        side,
                        plane: None,
                        faces: bvh,
                    }),
                )
            }
            "surface" => {
                let file = d.file.as_ref().ok_or_else(|| {
                    err(
                        &key,
                        format!("[[datum]] {name}: a surface needs `file` (STL)"),
                    )
                })?;
                let scale = d.scale.unwrap_or(1.0);
                let tris = foehncraft_mesh::format::stl::read(file)
                    .map_err(|e| err(&key, format!("[[datum]] {name}: {e}")))?;
                let mut sf = foehncraft_mesh::query::SurfaceFaces::default();
                let loc = foehncraft_mesh::locate::Locator::new(mesh, geo);
                let mut cut = Vec::new();
                for (k, t) in tris.iter().enumerate() {
                    let t = t.map(|v| v.map(|x| x * scale));
                    let a = cross(sub(t[1], t[0]), sub(t[2], t[0]));
                    let Some(nrm) = unit(a) else {
                        continue;
                    };
                    sf.tris.push(t);
                    sf.tri_face.push(sf.global_id.len() as u32);
                    sf.global_id.push(k as u64);
                    sf.mesh_face.push(usize::MAX);
                    // the direction of the hits: `entering` is against it
                    sf.normal.push(nrm);
                    sample_triangle(&loc, geo, t, &mut cut);
                }
                if sf.tris.is_empty() {
                    return Err(err(
                        &key,
                        format!("[[datum]] {name}: the surface has no triangle of nonzero area"),
                    ));
                }
                (
                    Facets::Cut(cut),
                    None,
                    Shape::Surface(Surface::Triangles {
                        bvh: Box::new(foehncraft_mesh::query::Bvh::build(sf)),
                    }),
                )
            }
            "box" | "cylinder" | "zone" => {
                let v = match d.kind.as_str() {
                    "box" => {
                        let (lo, hi) = (need(d.lower, "lower")?, need(d.upper, "upper")?);
                        if (0..3).any(|k| hi[k] <= lo[k]) {
                            return Err(err(
                                &key,
                                format!("[[datum]] {name}: upper must exceed lower"),
                            ));
                        }
                        Volume::Box {
                            lower: lo,
                            upper: hi,
                        }
                    }
                    "cylinder" => {
                        let (a, b) = (need(d.start, "start")?, need(d.end, "end")?);
                        let r =
                            d.radius
                                .filter(|r| r.is_finite() && *r > 0.0)
                                .ok_or_else(|| {
                                    err(
                                        &key,
                                        format!("[[datum]] {name}: a cylinder needs `radius` > 0"),
                                    )
                                })?;
                        if unit(sub(b, a)).is_none() {
                            return Err(err(
                                &key,
                                format!("[[datum]] {name}: start and end coincide"),
                            ));
                        }
                        Volume::Cylinder {
                            start: a,
                            end: b,
                            radius: r,
                        }
                    }
                    _ => {
                        let z = d.zone.as_deref().ok_or_else(|| {
                            err(&key, format!("[[datum]] {name}: a zone datum needs `zone`"))
                        })?;
                        let zi = mesh.zones.iter().position(|x| x == z).ok_or_else(|| {
                            err(
                                &key,
                                format!(
                                    "[[datum]] {name}: the mesh has no cell zone `{z}`; it has {}",
                                    mesh.zones.join(", ")
                                ),
                            )
                        })?;
                        Volume::Zone { zone: zi as u32 }
                    }
                };
                let cells: Vec<usize> = (0..mesh.n_cells)
                    .filter(|&c| v.contains(geo.cell_centroid[c], c, &mesh.cell_zone))
                    .collect();
                (Facets::Cells(cells), None, Shape::Volume(v))
            }
            k => return Err(err(&key, format!("[[datum]] {name}: unknown type `{k}`"))),
        };
        let empty = match &facets {
            Facets::Faces(f) => f.is_empty(),
            Facets::Cut(c) => c.is_empty(),
            Facets::Cells(c) => c.is_empty(),
        };
        if empty {
            return Err(diag::err(
                Diagnostic::new(
                    diag::DATUM_EMPTY,
                    format!("[[datum]] {name}: it meets no face or cell of the mesh"),
                )
                .key(&key),
            ));
        }
        out.push(Built {
            section: d.clone(),
            facets,
            snapped,
            parcels: Datum {
                name: name.clone(),
                shape,
                hist: d.hist(),
            },
        });
    }
    Ok(out)
}

/// One number of a report: quantity name, value.
pub type Value = (String, f64);

/// The gas side of a solution at a datum.
struct Gas<'a> {
    sol: &'a Solution<'a>,
    mesh: &'a Mesh,
    geo: &'a Geometry,
    offsets: Option<Vec<f64>>,
}

impl Gas<'_> {
    fn u_at(&self, c: usize, x: [f64; 3]) -> [f64; 3] {
        let d = sub(x, self.geo.cell_centroid[c]);
        let mut u = self.sol.u.get(c);
        for (i, ui) in u.iter_mut().enumerate() {
            let g = self.sol.grad_u[i].get(c);
            *ui += dot(g, d);
        }
        u
    }

    fn p_at(&self, c: usize, x: [f64; 3]) -> f64 {
        let d = sub(x, self.geo.cell_centroid[c]);
        let off = self.offsets.as_ref().map_or(0.0, |o| o[c]);
        self.sol.p[c] + dot(self.sol.grad_p.get(c), d) + off
    }

    /// Velocity, pressure and scalars on face `f`.
    fn face(&self, f: usize) -> ([f64; 3], f64, Vec<f64>) {
        let ni = self.mesh.n_interior_faces();
        let o = self.mesh.owner[f] as usize;
        let x = self.geo.face_centroid[f];
        let sc = |c: usize| -> Vec<f64> { self.sol.scalars.iter().map(|s| s.phi[c]).collect() };
        if f >= ni {
            let off = self.offsets.as_ref().map_or(0.0, |v| v[o]);
            (self.sol.ub.get(f - ni), self.sol.pb[f - ni] + off, sc(o))
        } else {
            let nb = self.mesh.neighbour[f] as usize;
            let (ua, ub) = (self.u_at(o, x), self.u_at(nb, x));
            let (sa, sb) = (sc(o), sc(nb));
            (
                [
                    0.5 * (ua[0] + ub[0]),
                    0.5 * (ua[1] + ub[1]),
                    0.5 * (ua[2] + ub[2]),
                ],
                0.5 * (self.p_at(o, x) + self.p_at(nb, x)),
                sa.iter().zip(&sb).map(|(a, b)| 0.5 * (a + b)).collect(),
            )
        }
    }
}

/// Sums of a surface report.
#[derive(Default)]
struct Sums {
    area: f64,
    mdot: f64,
    vdot: f64,
    mabs: f64,
    mom: [f64; 3],
    pa: f64,
    ua: [f64; 3],
    sa: Vec<f64>,
    sm: Vec<f64>,
}

impl Sums {
    fn add(&mut self, a: [f64; 3], mdot: f64, rho: f64, u: [f64; 3], p: f64, s: &[f64]) {
        let am = dot(a, a).sqrt();
        self.area += am;
        self.mdot += mdot;
        self.vdot += mdot / rho;
        self.mabs += mdot.abs();
        for (k, &uk) in u.iter().enumerate() {
            self.mom[k] += mdot * uk;
            self.ua[k] += uk * am;
        }
        self.pa += p * am;
        if self.sa.len() < s.len() {
            self.sa.resize(s.len(), 0.0);
            self.sm.resize(s.len(), 0.0);
        }
        for (i, v) in s.iter().enumerate() {
            self.sa[i] += v * am;
            self.sm[i] += v * mdot.abs();
        }
    }
}

fn face_mass_flow(sol: &Solution, faces: &[(usize, f64)]) -> f64 {
    faces.iter().map(|&(f, s)| s * sol.mass_flux[f]).sum()
}

/// The gas report of datum `b` (D262): surfaces give `area` (m²),
/// `mass_flow` (kg/s, along the datum's direction), `volume_flow` (m³/s),
/// `momentum_flow_x/y/z` (N), area averages `velocity_x/y/z_area` (m/s) and
/// `pressure_area` (Pa), for each scalar `<name>_area` and the
/// mass-flow-weighted `<name>_mass` (weights |ṁ|), and for a cut plane
/// `cut_error` (relative difference of its mass flow against the exact
/// flow through the same plane snapped to the faces). Volumes give
/// `volume` (m³) and the volume averages `velocity_x/y/z_volume`,
/// `pressure_volume`, `<name>_volume`.
pub fn gas_report(
    b: &Built,
    case: &FlowCase,
    sol: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
) -> Vec<Value> {
    let g = Gas {
        sol,
        mesh,
        geo,
        offsets: case
            .pressure_frame
            .as_ref()
            .and_then(|f| f.cell_offsets(case, mesh)),
    };
    let names: Vec<String> = sol.scalars.iter().map(|s| s.name.clone()).collect();
    let mut out: Vec<Value> = Vec::new();
    let surface = |s: Sums, out: &mut Vec<Value>| {
        let q = |a: f64, w: f64| if w > 0.0 { a / w } else { 0.0 };
        out.push(("area".into(), s.area));
        out.push(("mass_flow".into(), s.mdot));
        out.push(("volume_flow".into(), s.vdot));
        for (k, c) in ["x", "y", "z"].iter().enumerate() {
            out.push((format!("momentum_flow_{c}"), s.mom[k]));
        }
        for (k, c) in ["x", "y", "z"].iter().enumerate() {
            out.push((format!("velocity_{c}_area"), q(s.ua[k], s.area)));
        }
        out.push(("pressure_area".into(), q(s.pa, s.area)));
        for (i, n) in names.iter().enumerate() {
            out.push((
                format!("{n}_area"),
                q(s.sa.get(i).copied().unwrap_or(0.0), s.area),
            ));
            out.push((
                format!("{n}_mass"),
                q(s.sm.get(i).copied().unwrap_or(0.0), s.mabs),
            ));
        }
    };
    match &b.facets {
        Facets::Faces(faces) => {
            let mut s = Sums::default();
            for &(f, sign) in faces {
                let a = geo.face_area[f];
                let (u, p, sc) = g.face(f);
                s.add(
                    [sign * a[0], sign * a[1], sign * a[2]],
                    sign * sol.mass_flux[f],
                    sol.phys.rho_face_at(f),
                    u,
                    p,
                    &sc,
                );
            }
            surface(s, &mut out);
        }
        Facets::Cut(pieces) => {
            let mut s = Sums::default();
            for &(c, a, x) in pieces {
                let u = g.u_at(c, x);
                let rho = sol.phys.rho_at(c);
                let sc: Vec<f64> = sol.scalars.iter().map(|v| v.phi[c]).collect();
                s.add(a, rho * dot(u, a), rho, u, g.p_at(c, x), &sc);
            }
            let mdot = s.mdot;
            surface(s, &mut out);
            if let Some(faces) = &b.snapped {
                let exact = face_mass_flow(sol, faces);
                let scale = exact.abs().max(1e-300);
                out.push(("cut_error".into(), (mdot - exact).abs() / scale));
            }
        }
        Facets::Cells(cells) => {
            let (mut v, mut ua, mut pa) = (0.0, [0.0; 3], 0.0);
            let mut sa = vec![0.0; names.len()];
            for &c in cells {
                let w = geo.cell_volume[c];
                v += w;
                let u = sol.u.get(c);
                for k in 0..3 {
                    ua[k] += u[k] * w;
                }
                pa += (sol.p[c] + g.offsets.as_ref().map_or(0.0, |o| o[c])) * w;
                for (i, sc) in sol.scalars.iter().enumerate() {
                    sa[i] += sc.phi[c] * w;
                }
            }
            let q = |a: f64| if v > 0.0 { a / v } else { 0.0 };
            out.push(("volume".into(), v));
            for (k, c) in ["x", "y", "z"].iter().enumerate() {
                out.push((format!("velocity_{c}_volume"), q(ua[k])));
            }
            out.push(("pressure_volume".into(), q(pa)));
            for (i, n) in names.iter().enumerate() {
                out.push((format!("{n}_volume"), q(sa[i])));
            }
        }
    }
    out
}

/// The liquid columns of datum `i` (D263, D264), closing its report
/// interval at the spray's time: surfaces give the cumulative
/// `liquid_kg_pos`, `liquid_kg_neg` (kg along and against the direction),
/// `liquid_rate` (net kg/s over the interval), `crossings` (in the
/// interval), and the flux-weighted `d10`, `d32`, `dv50` (m) of all
/// crossings so far; volumes give `liquid_kg`, `parcels`, `d10`, `d32`,
/// `dv50` of the parcels inside now. Also returns the interval's sums
/// (surfaces) for the histogram file.
pub fn liquid_report(
    sp: &mut foehncraft_particles::Spray,
    i: usize,
    zone: &[u32],
) -> (Vec<Value>, Option<foehncraft_particles::datum::Flux>) {
    let spec = sp.datums()[i].hist;
    if matches!(sp.datums()[i].shape, Shape::Volume(_)) {
        let st = sp.moments(Some(i), zone).stats(&spec);
        return (
            vec![
                ("liquid_kg".into(), st.mass),
                ("parcels".into(), st.parcels as f64),
                ("d10".into(), st.d10),
                ("d32".into(), st.d32),
                ("dv50".into(), st.dv50),
            ],
            None,
        );
    }
    let t = sp.time;
    let tally = &mut sp.datum_tally[i];
    let (f, dt) = tally.close(&spec, t);
    let st = tally.total.all.stats(&spec);
    (
        vec![
            ("liquid_kg_pos".into(), tally.total.mass_pos),
            ("liquid_kg_neg".into(), tally.total.mass_neg),
            (
                "liquid_rate".into(),
                if dt > 0.0 { f.net() / dt } else { 0.0 },
            ),
            ("crossings".into(), f.crossings() as f64),
            ("d10".into(), st.d10),
            ("d32".into(), st.d32),
            ("dv50".into(), st.dv50),
        ],
        Some(f),
    )
}

/// Units of the columns of `datums.csv`.
pub const CSV_UNITS: &str = "step -; time s; <datum>.area m²; .mass_flow kg/s; .volume_flow m³/s; \
     .momentum_flow_* N; .velocity_*_area m/s; .pressure_area Pa; .<scalar>_area, .<scalar>_mass \
     (scalar unit; mass-flow weighted with |ṁ|); .cut_error - (cut plane against the same plane \
     snapped to the faces); volumes: .volume m³, .*_volume (volume averages); sprays: \
     .liquid_kg_pos, .liquid_kg_neg kg (cumulative along/against the datum's direction), \
     .liquid_rate kg/s (net, over the report interval), .crossings - (in the interval), .d10 .d32 \
     .dv50 m (all crossings so far, flux-weighted); volumes .liquid_kg kg, .parcels -, .d10 .d32 \
     .dv50 m (parcels inside now)";

/// The reports of a run: `datums.csv` and `spray_hist_<datum>.csv`.
pub struct Reports {
    /// The datums.
    pub built: Vec<Built>,
    csv: Option<crate::monitor::History>,
    hist: Vec<Option<crate::monitor::History>>,
    header: Option<String>,
    out: PathBuf,
    name: String,
    cont: Option<usize>,
    /// The last row's values per datum (for `result.json`).
    pub last: Vec<Vec<Value>>,
}

impl Reports {
    /// The reports of the datums `built` of a run writing into `out`.
    pub fn new(built: Vec<Built>, out: &Path, name: &str, cont: Option<usize>) -> Self {
        let n = built.len();
        Self {
            built,
            csv: None,
            hist: (0..n).map(|_| None).collect(),
            header: None,
            out: out.to_path_buf(),
            name: name.into(),
            cont,
            last: vec![Vec::new(); n],
        }
    }

    /// Writes one row at step `n`, time `t`.
    #[allow(clippy::too_many_arguments)]
    pub fn row(
        &mut self,
        n: usize,
        t: f64,
        case: &FlowCase,
        sol: &Solution,
        mesh: &Mesh,
        geo: &Geometry,
        mut spray: Option<&mut foehncraft_particles::Spray>,
    ) -> Result<()> {
        let mut head = String::from("step,time");
        let mut row = format!("{n},{t:e}");
        for (i, b) in self.built.iter().enumerate() {
            let mut vals = gas_report(b, case, sol, mesh, geo);
            if let Some(sp) = spray.as_deref_mut() {
                let (liq, flux) = liquid_report(sp, i, &mesh.cell_zone);
                vals.extend(liq);
                if let Some(f) = flux {
                    let spec = b.parcels.hist;
                    if self.hist[i].is_none() {
                        self.hist[i] = Some(crate::monitor::History::new(
                            &self.out,
                            &format!("spray_hist_{}.csv", b.section.name),
                            &crate::case::csv_preamble(
                                "spray-histogram",
                                &self.name,
                                None,
                                "step -; time s (end of the interval); bin -; d_lo m; d_hi m; \
                                 droplets - (crossings of the interval, both directions); mass kg; \
                                 bin -1 is below d_lo of bin 0, bin = bins above the last",
                            ),
                            "step,time,bin,d_lo,d_hi,droplets,mass",
                            self.cont,
                        )?);
                    }
                    let h = self.hist[i].as_mut().expect("just made");
                    let e = spec.edges();
                    let a = &f.all;
                    h.row(&format!(
                        "{n},{t:e},-1,0,{:e},{:e},{:e}",
                        e[0], a.under[0], a.under[1]
                    ));
                    for k in 0..spec.bins {
                        h.row(&format!(
                            "{n},{t:e},{k},{:e},{:e},{:e},{:e}",
                            e[k],
                            e[k + 1],
                            a.hist_n[k],
                            a.hist_m[k]
                        ));
                    }
                    h.row(&format!(
                        "{n},{t:e},{},{:e},inf,{:e},{:e}",
                        spec.bins, e[spec.bins], a.over[0], a.over[1]
                    ));
                }
            }
            for (q, v) in &vals {
                let _ = write!(head, ",{}.{q}", b.section.name);
                let _ = write!(row, ",{v:e}");
            }
            self.last[i] = vals;
        }
        if self.csv.is_none() {
            self.csv = Some(crate::monitor::History::new(
                &self.out,
                "datums.csv",
                &crate::case::csv_preamble("datums", &self.name, None, CSV_UNITS),
                &head,
                self.cont,
            )?);
            self.header = Some(head);
        }
        if let Some(c) = self.csv.as_mut() {
            c.row(&row);
        }
        Ok(())
    }

    /// Appends the buffered rows to the files.
    pub fn flush(&mut self) -> Result<()> {
        if let Some(c) = self.csv.as_mut() {
            c.flush()?;
        }
        for h in self.hist.iter_mut().flatten() {
            h.flush()?;
        }
        Ok(())
    }

    /// The summary for `result.json`.
    pub fn summary(
        &self,
        spray: Option<&foehncraft_particles::Spray>,
        zone: &[u32],
    ) -> Vec<DatumSummary> {
        self.built
            .iter()
            .enumerate()
            .map(|(i, b)| {
                let spec = b.parcels.hist;
                let liquid = spray.map(|sp| {
                    if b.is_volume() {
                        LiquidSummary {
                            mass_pos: None,
                            mass_neg: None,
                            crossings: None,
                            mean_rate: None,
                            rate_relative_error: None,
                            statistics: sp.moments(Some(i), zone).stats(&spec),
                        }
                    } else {
                        let t = &sp.datum_tally[i];
                        LiquidSummary {
                            mass_pos: Some(t.total.mass_pos),
                            mass_neg: Some(t.total.mass_neg),
                            crossings: Some(t.total.crossings()),
                            mean_rate: Some(t.mean_rate(sp.time)),
                            rate_relative_error: t.rate_relative_error(),
                            statistics: t.total.all.stats(&spec),
                        }
                    }
                });
                DatumSummary {
                    name: b.section.name.clone(),
                    kind: b.section.kind.clone(),
                    exact: matches!(b.facets, Facets::Faces(_)),
                    gas: self.last[i].iter().cloned().collect(),
                    liquid,
                }
            })
            .collect()
    }
}

/// A datum in `result.json` (D267).
#[derive(Debug, Clone, Serialize)]
pub struct DatumSummary {
    /// Name.
    pub name: String,
    /// Type.
    pub kind: String,
    /// Sums over mesh faces (exact, conservative), not a cut.
    pub exact: bool,
    /// The values of the last report, by quantity.
    pub gas: std::collections::BTreeMap<String, f64>,
    /// The liquid of the sprays.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub liquid: Option<LiquidSummary>,
}

/// The liquid at a datum over the run (D263, D264).
#[derive(Debug, Clone, Serialize)]
pub struct LiquidSummary {
    /// Cumulative mass along the direction (kg).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mass_pos: Option<f64>,
    /// Cumulative mass against it (kg).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mass_neg: Option<f64>,
    /// Crossings.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub crossings: Option<u64>,
    /// Mean net mass rate since the start of the sums (kg/s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mean_rate: Option<f64>,
    /// Relative standard error of the mean rate (ten batch means of the
    /// report intervals, D71); absent with fewer than ten intervals.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rate_relative_error: Option<f64>,
    /// Droplet statistics (flux-weighted at surfaces, of the parcels inside
    /// for volumes).
    pub statistics: foehncraft_particles::Stats,
}

#[cfg(test)]
mod tests {
    use super::*;
    use foehncraft_mesh::generate::{BoxSpec, Kind};

    fn mesh(kind: Kind) -> (Mesh, Geometry) {
        let mut s = BoxSpec::new(kind, 3, [5, 5, 5]);
        s.lower = [0.0; 3];
        s.upper = [1.0; 3];
        let m = s.generate().unwrap();
        let g = Geometry::compute(&m);
        (m, g)
    }

    /// The cut polygons of a plane through any mesh add up to the area of
    /// the plane in the box, with the centroid at its centre.
    #[test]
    fn cuts_add_up_to_the_section() {
        for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
            let (m, g) = mesh(kind);
            for (p, n) in [
                ([0.37, 0.5, 0.5], [1.0, 0.0, 0.0]),
                ([0.5, 0.5, 0.5], unit([1.0, 1.0, 0.0]).unwrap()),
            ] {
                let (e1, e2) = frame(n);
                let (mut a, mut mx) = (0.0, [0.0; 3]);
                for c in 0..m.n_cells {
                    if let Some((ac, x)) = cut_cell(&m, &g, c, p, n, e1, e2, None) {
                        a += ac;
                        for k in 0..3 {
                            mx[k] += ac * x[k];
                        }
                    }
                }
                let exact = if n[1] == 0.0 { 1.0 } else { 2f64.sqrt() };
                assert!((a - exact).abs() < 1e-12, "{kind:?}: area {a}");
                let cy = mx[2] / a;
                assert!((cy - 0.5).abs() < 1e-12, "{kind:?}: centroid z {cy}");
            }
            // a disc of radius 0.3 about the centre
            let p = [0.5, 0.5, 0.5];
            let n = [0.0, 0.0, 1.0];
            let (e1, e2) = frame(n);
            let poly = bound_polygon(&Bound::Disc { radius: 0.3 }).unwrap();
            let a: f64 = (0..m.n_cells)
                .filter_map(|c| cut_cell(&m, &g, c, p, n, e1, e2, Some(&poly)))
                .map(|x| x.0)
                .sum();
            assert!(
                (a - std::f64::consts::PI * 0.09).abs() < 1e-12,
                "{kind:?}: disc {a}"
            );
        }
    }
}
