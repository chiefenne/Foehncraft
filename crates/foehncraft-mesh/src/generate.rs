//! Built-in mesh generator for canonical test cases (D34).
//!
//! All meshes fill an axis-aligned box `lower..upper` with `cells[d]`
//! intervals per direction, graded or uniform, optionally refined by whole
//! levels and with randomly perturbed nodes:
//!
//! - `hex`: hexahedra;
//! - `tet`: each hexahedron split into six tetrahedra around its main
//!   diagonal (Kuhn's subdivision, which gives matching diagonals on shared
//!   faces);
//! - `poly`: polyhedra, the dual of the `tet` mesh: one cell per node of the
//!   tetrahedral mesh, one face per tetrahedral edge (a polygon through the
//!   centroids of the triangles and tetrahedra around the edge), boundary
//!   faces merged into one planar polygon per node and boundary; interior cells
//!   have 14 faces on a regular grid;
//! - `mixed`: three bands along x: hexahedra, prisms (each hexahedron split
//!   into two along the diagonal of its z faces) and tetrahedra. Where a
//!   quadrilateral meets two triangles it is split, so the cells along the
//!   band interfaces are polyhedra.
//!
//! A 2D mesh (`dimension = 2`) is one cell layer thick in z (D3): `hex`
//! gives hexahedra, `tet` prisms on a triangulation, `poly` prisms on the
//! dual polygons of that triangulation (mostly hexagons), `mixed` a band of
//! hexahedra and a band of prisms. The two z boundaries take the symmetry
//! condition in a 2D case.
//!
//! Boundary boundaries are named `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax`.
//!
//! **Grading.** Per direction a ratio r of the last to the first interval
//! (one-sided, r < 1 clusters at the upper end), or with `symmetric` of the
//! middle to the end intervals (clustering at both walls, for boundary
//! layers). The intervals follow a continuous exponential mapping fixed by
//! the base mesh, so a refined mesh contains the nodes of the coarser one
//! (nested refinement for order-of-accuracy studies).
//!
//! **Segments.** Instead of one grading, a direction may be a chain of
//! segments, each with its upper coordinate `to`, its number of base
//! intervals and its ratio of the last to the first interval (a flat plate:
//! intervals that shrink towards the leading edge from both sides). The
//! intervals of a direction sum to `cells`; refinement halves every one.
//!
//! **Split boundaries.** `split` cuts a box boundary at a coordinate: the
//! faces whose centroid lies above `at` along `axis` form the boundary
//! `name` (the plate after its symmetry stretch).
//!
//! **Removed boxes.** `remove` takes the cells whose centroid lies inside
//! a box out of the mesh (hexahedra, tetrahedra and mixed meshes): a
//! backward-facing step is the box lower corner to the step corner. The
//! faces that the removal exposes form the boundary `name`; nodes that no
//! cell uses are dropped.
//!
//! **Zones.** `zones` puts the cells whose centroid lies inside a box in a
//! cell zone of that name; the other cells form the zone `fluid`. The
//! solids of a conjugate heat transfer case are such zones, named in
//! `[regions]`.
//!
//! **Perturbation.** Each node moves by up to `perturb` times the smaller
//! adjacent interval in each direction, from a counter-based random number
//! keyed by the seed and the node position (D47), so it is reproducible.
//! Boundary nodes move only along the boundary, and so do the nodes on a
//! face plane of a `zones` box (the interfaces of a conjugate case stay
//! planes); in 2D the nodes move in x and y only.

use crate::builder::{BoundaryFace, CellSoup, MeshBuilder};
use crate::geometry::{add, dot, polygon_area_centroid, scale, sub};
use crate::mesh::Mesh;
use foehncraft_core::rng;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Cell type of a generated mesh.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Kind {
    /// Hexahedra.
    #[default]
    Hex,
    /// Tetrahedra (prisms in 2D).
    Tet,
    /// Polyhedra (polygonal prisms in 2D).
    Poly,
    /// Hexahedra, prisms, tetrahedra and transitional polyhedra.
    Mixed,
}

impl Kind {
    /// Name used in case files.
    pub fn name(self) -> &'static str {
        match self {
            Self::Hex => "hex",
            Self::Tet => "tet",
            Self::Poly => "poly",
            Self::Mixed => "mixed",
        }
    }
}

/// Grading of one direction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Grading {
    /// Ratio of the last to the first interval of the base mesh (with
    /// `symmetric`: of the middle to the end intervals).
    #[serde(default = "one")]
    pub ratio: f64,
    /// Cluster at both ends.
    #[serde(default)]
    pub symmetric: bool,
}

fn one() -> f64 {
    1.0
}

impl Default for Grading {
    fn default() -> Self {
        Self {
            ratio: 1.0,
            symmetric: false,
        }
    }
}

impl Grading {
    /// Node positions in [0, 1] for `n` intervals at refinement level
    /// `level` (`n` already refined; `n0 = n >> level` is the base count).
    pub fn nodes(&self, n: usize, level: u32) -> Vec<f64> {
        let n0 = (n >> level).max(1);
        let one_sided = |xi: f64, m0: usize| -> f64 {
            if (self.ratio - 1.0).abs() < 1e-14 || m0 < 2 {
                return xi;
            }
            let beta = self.ratio.ln() * m0 as f64 / (m0 as f64 - 1.0);
            ((beta * xi).exp() - 1.0) / (beta.exp() - 1.0)
        };
        (0..=n)
            .map(|i| {
                let xi = i as f64 / n as f64;
                if self.symmetric {
                    let m0 = n0 / 2;
                    if xi <= 0.5 {
                        0.5 * one_sided(2.0 * xi, m0)
                    } else {
                        1.0 - 0.5 * one_sided(2.0 * (1.0 - xi), m0)
                    }
                } else {
                    one_sided(xi, n0)
                }
            })
            .collect()
    }
}

/// Specification of a generated box mesh.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BoxSpec {
    /// Cell type.
    #[serde(default, rename = "type")]
    pub kind: Kind,
    /// 2 for a mesh one cell layer thick in z, 3 otherwise.
    #[serde(default = "three")]
    pub dimension: u8,
    /// Lower corner.
    #[serde(default = "zeros")]
    pub lower: [f64; 3],
    /// Upper corner.
    #[serde(default = "ones")]
    pub upper: [f64; 3],
    /// Intervals per direction of the base mesh (z ignored in 2D).
    #[serde(default = "default_cells")]
    pub cells: [usize; 3],
    /// Grading per direction.
    #[serde(default)]
    pub grading: GradingSpec,
    /// Refinement level: every interval is halved `refine` times.
    #[serde(default)]
    pub refine: u32,
    /// Node perturbation, as a fraction of the local interval (0 to 0.3).
    #[serde(default)]
    pub perturb: f64,
    /// Seed of the perturbation.
    #[serde(default = "default_seed")]
    pub seed: u64,
    /// Revolve the 2D mesh about the x axis into an axisymmetric wedge of
    /// this angle (degrees; y is the radius, `zmin` and `zmax` become the
    /// wedge planes, see [`crate::wedge`]). With `dimension = 3` the z
    /// intervals become the angular cells of a sector of the angle (D309),
    /// whose side planes `zmin` and `zmax` are ordinary boundaries, for
    /// example symmetry planes; a sector is of hexahedra, not perturbed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wedge: Option<f64>,
    /// Directions whose two boundaries form a periodic pair (set by the
    /// flow case from its periodic boundaries, D29): the perturbation of
    /// the nodes on the two boundaries is then the same, so their faces
    /// are translated copies.
    #[serde(skip)]
    #[schemars(skip)]
    pub periodic: [bool; 3],
    /// Graded segments per direction (replace `grading` in that direction).
    #[serde(default, skip_serializing_if = "SegmentSpec::is_empty")]
    pub segments: SegmentSpec,
    /// Boundaries cut at a coordinate.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub split: Vec<Split>,
    /// Boxes of cells taken out of the mesh.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub remove: Vec<RemoveBox>,
    /// Boxes of cells put in a cell zone of their own (the rest stays in
    /// the zone `fluid`): the solid parts of a conjugate heat transfer
    /// case, named in `[regions]`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub zones: Vec<ZoneBox>,
}

/// A box of cells of a generated mesh in a cell zone of its own: the cells
/// whose centroid lies inside the box (later boxes win where they
/// overlap).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ZoneBox {
    /// Lower corner (m; z is ignored in 2D).
    pub lower: [f64; 3],
    /// Upper corner (m; z is ignored in 2D).
    pub upper: [f64; 3],
    /// Name of the cell zone.
    pub name: String,
}

impl ZoneBox {
    fn contains(&self, dim: u8, p: [f64; 3]) -> bool {
        RemoveBox {
            lower: self.lower,
            upper: self.upper,
            name: String::new(),
        }
        .contains(dim, p)
    }
}

/// A box of cells taken out of a generated mesh.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RemoveBox {
    /// Lower corner (m; z is ignored in 2D).
    pub lower: [f64; 3],
    /// Upper corner (m; z is ignored in 2D).
    pub upper: [f64; 3],
    /// Name of the boundary that the removal exposes.
    pub name: String,
}

impl RemoveBox {
    /// The box in the dimension `dim` (2D: the whole z range).
    fn range(&self, dim: u8) -> ([f64; 3], [f64; 3]) {
        let (mut lo, mut hi) = (self.lower, self.upper);
        if dim == 2 {
            lo[2] = f64::NEG_INFINITY;
            hi[2] = f64::INFINITY;
        }
        (lo, hi)
    }

    fn contains(&self, dim: u8, p: [f64; 3]) -> bool {
        let (lo, hi) = self.range(dim);
        (0..3).all(|d| p[d] > lo[d] && p[d] < hi[d])
    }
}

/// One graded segment of a direction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Segment {
    /// Upper coordinate of the segment (m); the last one is `upper`.
    pub to: f64,
    /// Base intervals of the segment.
    pub cells: usize,
    /// Ratio of the last to the first interval (< 1: smaller at the upper
    /// end).
    #[serde(default = "one")]
    pub ratio: f64,
}

/// The segments of the three directions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SegmentSpec {
    /// x direction.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub x: Vec<Segment>,
    /// y direction.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub y: Vec<Segment>,
    /// z direction.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub z: Vec<Segment>,
}

impl SegmentSpec {
    /// No segments.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty() && self.y.is_empty() && self.z.is_empty()
    }

    /// The segments of direction `d`.
    pub fn of(&self, d: usize) -> &[Segment] {
        match d {
            0 => &self.x,
            1 => &self.y,
            _ => &self.z,
        }
    }
}

/// A box boundary cut at a coordinate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Split {
    /// The box boundary (`xmin`, …, `zmax`).
    pub boundary: String,
    /// The direction of the cut: `x`, `y` or `z`.
    pub axis: String,
    /// Coordinate of the cut (m).
    pub at: f64,
    /// Name of the part above `at`.
    pub name: String,
}

impl Split {
    fn axis_index(&self) -> Option<usize> {
        match self.axis.as_str() {
            "x" => Some(0),
            "y" => Some(1),
            "z" => Some(2),
            _ => None,
        }
    }
}

/// Grading of the three directions.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GradingSpec {
    /// x direction.
    #[serde(default)]
    pub x: Grading,
    /// y direction.
    #[serde(default)]
    pub y: Grading,
    /// z direction.
    #[serde(default)]
    pub z: Grading,
}

fn three() -> u8 {
    3
}
fn zeros() -> [f64; 3] {
    [0.0; 3]
}
fn ones() -> [f64; 3] {
    [1.0; 3]
}
fn default_cells() -> [usize; 3] {
    [8, 8, 8]
}
fn default_seed() -> u64 {
    2026
}

impl Default for BoxSpec {
    fn default() -> Self {
        Self {
            kind: Kind::Hex,
            dimension: 3,
            lower: zeros(),
            upper: ones(),
            cells: default_cells(),
            grading: GradingSpec::default(),
            refine: 0,
            perturb: 0.0,
            seed: default_seed(),
            wedge: None,
            periodic: [false; 3],
            segments: SegmentSpec::default(),
            split: Vec::new(),
            remove: Vec::new(),
            zones: Vec::new(),
        }
    }
}

/// Boundary names of a box, in face order.
pub const BOX_BOUNDARIES: [&str; 6] = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"];

impl BoxSpec {
    /// A spec with the given type, dimension and intervals.
    pub fn new(kind: Kind, dimension: u8, cells: [usize; 3]) -> Self {
        Self {
            kind,
            dimension,
            cells,
            ..Default::default()
        }
    }

    /// Intervals per direction after refinement (z = 1 in 2D).
    pub fn intervals(&self) -> [usize; 3] {
        let f = 1usize << self.refine;
        [
            self.cells[0] * f,
            self.cells[1] * f,
            if self.dimension == 2 {
                1
            } else {
                self.cells[2] * f
            },
        ]
    }

    /// Checks the spec.
    pub fn validate(&self) -> Result<(), String> {
        if self.dimension != 2 && self.dimension != 3 {
            return Err(format!("dimension must be 2 or 3, not {}", self.dimension));
        }
        for d in 0..3 {
            if self.upper[d] <= self.lower[d] {
                return Err(format!("upper[{d}] must be larger than lower[{d}]"));
            }
        }
        let n = self.intervals();
        if n.contains(&0) {
            return Err("at least one interval per direction".into());
        }
        if self.kind == Kind::Mixed && n[0] < if self.dimension == 2 { 2 } else { 3 } {
            return Err("a mixed mesh needs at least 3 intervals in x (2 in 2D)".into());
        }
        if !(0.0..=0.3).contains(&self.perturb) {
            return Err(format!(
                "perturb must lie in [0, 0.3], not {}",
                self.perturb
            ));
        }
        for (d, g) in [self.grading.x, self.grading.y, self.grading.z]
            .iter()
            .enumerate()
        {
            if !(g.ratio > 0.0 && g.ratio.is_finite()) {
                return Err(format!("grading ratio of direction {d} must be positive"));
            }
        }
        for d in 0..3 {
            let seg = self.segments.of(d);
            if seg.is_empty() {
                continue;
            }
            let name = ["x", "y", "z"][d];
            if d == 2 && self.dimension == 2 {
                return Err("segments in z: a 2D mesh is one layer thick".into());
            }
            let total: usize = seg.iter().map(|s| s.cells).sum();
            if total != self.cells[d] {
                return Err(format!(
                    "the segments in {name} have {total} intervals, cells[{d}] is {}",
                    self.cells[d]
                ));
            }
            let mut lo = self.lower[d];
            for (i, s) in seg.iter().enumerate() {
                if s.cells == 0 || !(s.ratio > 0.0 && s.ratio.is_finite()) || s.to <= lo {
                    return Err(format!(
                        "segment {i} in {name}: needs cells >= 1, a positive ratio and `to` above \
                         the segment before"
                    ));
                }
                lo = s.to;
            }
            if (lo - self.upper[d]).abs() > 1e-12 * (self.upper[d] - self.lower[d]) {
                return Err(format!(
                    "the last segment in {name} must end at upper[{d}] = {}",
                    self.upper[d]
                ));
            }
        }
        for s in &self.split {
            let Some(a) = s.axis_index() else {
                return Err(format!("split of `{}`: axis must be x, y or z", s.boundary));
            };
            if !BOX_BOUNDARIES.contains(&s.boundary.as_str()) {
                return Err(format!(
                    "split: `{}` is not a box boundary ({})",
                    s.boundary,
                    BOX_BOUNDARIES.join(", ")
                ));
            }
            if s.boundary.starts_with(["x", "y", "z"][a]) {
                return Err(format!(
                    "split of `{}`: the axis must lie in the boundary",
                    s.boundary
                ));
            }
            if BOX_BOUNDARIES.contains(&s.name.as_str())
                || self.split.iter().filter(|t| t.name == s.name).count() > 1
            {
                return Err(format!("split: the name `{}` is taken", s.name));
            }
            if !(s.at > self.lower[a] && s.at < self.upper[a]) {
                return Err(format!(
                    "split of `{}`: `at` must lie inside the box",
                    s.boundary
                ));
            }
        }
        for r in &self.remove {
            if self.kind == Kind::Poly {
                return Err("remove: not with polyhedral meshes (hex, tet or mixed)".into());
            }
            let dims = if self.dimension == 2 { 2 } else { 3 };
            if (0..dims)
                .any(|d| r.upper[d].partial_cmp(&r.lower[d]) != Some(std::cmp::Ordering::Greater))
            {
                return Err(format!(
                    "remove `{}`: upper must be larger than lower",
                    r.name
                ));
            }
            if BOX_BOUNDARIES.contains(&r.name.as_str())
                || self.split.iter().any(|t| t.name == r.name)
            {
                return Err(format!("remove: the name `{}` is taken", r.name));
            }
        }
        for (i, z) in self.zones.iter().enumerate() {
            let dims = if self.dimension == 2 { 2 } else { 3 };
            if (0..dims)
                .any(|d| z.upper[d].partial_cmp(&z.lower[d]) != Some(std::cmp::Ordering::Greater))
            {
                return Err(format!(
                    "zone `{}`: upper must be larger than lower",
                    z.name
                ));
            }
            if z.name.is_empty()
                || z.name == "fluid"
                || self.zones[..i].iter().any(|y| y.name == z.name)
            {
                return Err(format!(
                    "zones: the name `{}` is empty, `fluid` (the zone of the other cells) or \
                     given twice",
                    z.name
                ));
            }
        }
        if let Some(a) = self.wedge {
            if self.lower[1] < 0.0 {
                return Err(format!(
                    "a wedge needs lower[1] >= 0 (y is the radius), not {}",
                    self.lower[1]
                ));
            }
            if !(a > 0.0 && a <= 90.0) {
                return Err(format!("the wedge angle must lie in (0°, 90°], not {a}°"));
            }
            if self.dimension == 3 && self.kind != Kind::Hex {
                return Err(
                    "a sector (wedge with dimension = 3) needs type = \"hex\": other cells collapse on the axis"
                        .into(),
                );
            }
            if self.dimension == 3 && self.perturb != 0.0 {
                return Err(
                    "a sector (wedge with dimension = 3) cannot be perturbed: perturb would move the nodes on the axis apart"
                        .into(),
                );
            }
        }
        let est = self.estimated_cells();
        if est > foehncraft_core::index::MAX_LOCAL_INDEX as usize / 8 {
            return Err(format!(
                "about {est} cells: too many for one generated mesh"
            ));
        }
        Ok(())
    }

    /// Approximate number of cells.
    pub fn estimated_cells(&self) -> usize {
        let [a, b, c] = self.intervals();
        let hexes = a * b * c;
        match (self.kind, self.dimension) {
            (Kind::Hex, _) => hexes,
            (Kind::Tet, 3) => 6 * hexes,
            (Kind::Tet, _) => 2 * hexes,
            (Kind::Poly, 3) => (a + 1) * (b + 1) * (c + 1),
            (Kind::Poly, _) => (a + 1) * (b + 1),
            (Kind::Mixed, 3) => hexes * (1 + 2 + 6) / 3,
            (Kind::Mixed, _) => hexes * 3 / 2,
        }
    }

    /// Generates the mesh.
    pub fn generate(&self) -> Result<Mesh, String> {
        self.validate()?;
        let m = match (self.kind, self.dimension) {
            (Kind::Poly, 3) => self.poly3d(),
            (Kind::Poly, _) => self.poly2d(),
            _ => self.structured(),
        }?;
        let m = self.with_zones(m)?;
        match self.wedge {
            Some(a) => crate::wedge::revolve(&m, a),
            None => Ok(m),
        }
    }

    /// The cell zones of `zones`: a cell in the last box that holds its
    /// centroid, the rest in `fluid`.
    fn with_zones(&self, mut m: Mesh) -> Result<Mesh, String> {
        if self.zones.is_empty() {
            return Ok(m);
        }
        let geo = crate::Geometry::compute(&m);
        let mut used = vec![false; self.zones.len()];
        m.zones = std::iter::once("fluid".to_string())
            .chain(self.zones.iter().map(|z| z.name.clone()))
            .collect();
        for c in 0..m.n_cells {
            let x = geo.cell_centroid[c];
            m.cell_zone[c] = match self
                .zones
                .iter()
                .rposition(|z| z.contains(self.dimension, x))
            {
                Some(i) => {
                    used[i] = true;
                    1 + i as u32
                }
                None => 0,
            };
        }
        if let Some(i) = used.iter().position(|u| !u) {
            return Err(format!("zone `{}` holds no cell", self.zones[i].name));
        }
        Ok(m)
    }

    /// Node coordinates of the (perturbed) grid, index `i + nx1 (j + ny1 k)`.
    fn grid(&self) -> (Vec<[f64; 3]>, [usize; 3]) {
        let n = self.intervals();
        let g = [self.grading.x, self.grading.y, self.grading.z];
        let coords: Vec<Vec<f64>> = (0..3)
            .map(|d| {
                let level = if d == 2 && self.dimension == 2 {
                    0
                } else {
                    self.refine
                };
                let seg = self.segments.of(d);
                if seg.is_empty() {
                    g[d].nodes(n[d], level)
                        .into_iter()
                        .map(|s| self.lower[d] + s * (self.upper[d] - self.lower[d]))
                        .collect()
                } else {
                    let f = 1usize << level;
                    let mut x = vec![self.lower[d]];
                    let mut lo = self.lower[d];
                    for s in seg {
                        let gr = Grading {
                            ratio: s.ratio,
                            symmetric: false,
                        };
                        let nodes = gr.nodes(s.cells * f, level);
                        x.extend(nodes[1..].iter().map(|t| lo + t * (s.to - lo)));
                        lo = s.to;
                    }
                    let last = x.len() - 1;
                    x[last] = self.upper[d];
                    x
                }
            })
            .collect();
        let n1 = [n[0] + 1, n[1] + 1, n[2] + 1];
        let mut pts = Vec::with_capacity(n1[0] * n1[1] * n1[2]);
        // nodes on a face plane of a zone box move only along it, so that
        // the interfaces between the zones stay planes (stage 7d2)
        let ext: Vec<f64> = (0..3)
            .map(|d| (self.upper[d] - self.lower[d]).abs())
            .collect();
        let on_zone_plane = |d: usize, x: f64| {
            self.zones.iter().any(|z| {
                (z.lower[d] - x).abs() <= 1e-12 * ext[d] || (z.upper[d] - x).abs() <= 1e-12 * ext[d]
            })
        };
        let dims = if self.dimension == 2 { 2 } else { 3 };
        for k in 0..n1[2] {
            for j in 0..n1[1] {
                for i in 0..n1[0] {
                    let idx = [i, j, k];
                    let mut p = [coords[0][i], coords[1][j], coords[2][k]];
                    if self.perturb > 0.0 {
                        // key: node position in the 2D or 3D grid, so both
                        // layers of a 2D mesh move alike; in a periodic
                        // direction the last node layer takes the numbers
                        // of the first
                        let wrap = |e: usize, m: usize| {
                            if self.periodic[e] && m == n[e] {
                                0
                            } else {
                                m as u64
                            }
                        };
                        let kk = if dims == 2 { 0 } else { wrap(2, k) };
                        let (ki, kj) = (wrap(0, i), wrap(1, j));
                        for d in 0..dims {
                            let m = idx[d];
                            if m == 0 || m == n[d] || on_zone_plane(d, coords[d][m]) {
                                continue;
                            }
                            let h = (coords[d][m] - coords[d][m - 1])
                                .min(coords[d][m + 1] - coords[d][m]);
                            let u = rng::uniform_keyed([
                                foehncraft_core::rng::stream::MESH,
                                self.seed,
                                ki | (kj << 21) | (kk << 42),
                                d as u64,
                            ]);
                            p[d] += self.perturb * h * (2.0 * u - 1.0);
                        }
                    }
                    pts.push(p);
                }
            }
        }
        (pts, n1)
    }

    /// The boundary names in face order: the box boundaries, then the
    /// split parts.
    fn boundary_names(&self) -> Vec<String> {
        BOX_BOUNDARIES
            .iter()
            .map(|s| s.to_string())
            .chain(self.split.iter().map(|s| s.name.clone()))
            .chain(self.remove.iter().map(|r| r.name.clone()))
            .collect()
    }

    fn classifier(&self) -> impl Fn(BoundaryFace) -> String + use<> {
        let (lo, hi) = (self.lower, self.upper);
        let split: Vec<(String, usize, f64, String)> = self
            .split
            .iter()
            .filter_map(|s| Some((s.boundary.clone(), s.axis_index()?, s.at, s.name.clone())))
            .collect();
        let dim = self.dimension;
        let remove: Vec<([f64; 3], [f64; 3], String)> = self
            .remove
            .iter()
            .map(|r| {
                let (a, b) = r.range(dim);
                (a, b, r.name.clone())
            })
            .collect();
        move |b: BoundaryFace| {
            let d = (0..3)
                .max_by(|&a, &c| b.normal[a].abs().total_cmp(&b.normal[c].abs()))
                .unwrap();
            let tol = 1e-9 * (hi[d] - lo[d]);
            // a face that a removed box exposes: on one of its sides,
            // inside the box (not on the outer boundary)
            let outer = (b.centroid[d] - lo[d]).abs() < tol || (b.centroid[d] - hi[d]).abs() < tol;
            if !outer {
                for (rl, rh, name) in &remove {
                    let on =
                        (b.centroid[d] - rl[d]).abs() < tol || (b.centroid[d] - rh[d]).abs() < tol;
                    let inside = (0..3)
                        .filter(|&e| e != d)
                        .all(|e| b.centroid[e] > rl[e] - tol && b.centroid[e] < rh[e] + tol);
                    if on && inside {
                        return name.clone();
                    }
                }
            }
            let s = if b.normal[d] < 0.0 { 0 } else { 1 };
            debug_assert!(
                (b.centroid[d] - if s == 0 { lo[d] } else { hi[d] }).abs() < tol,
                "face {:?} at {:?} normal {:?} is not on the box boundary",
                b.nodes,
                b.centroid,
                b.normal
            );
            let name = BOX_BOUNDARIES[2 * d + s];
            for (bn, a, at, new) in &split {
                if bn == name && b.centroid[*a] > *at {
                    return new.clone();
                }
            }
            name.to_string()
        }
    }

    /// Hexahedra, tetrahedra, prisms and mixed meshes through the cell soup.
    fn structured(&self) -> Result<Mesh, String> {
        let (pts, n1) = self.grid();
        let n = self.intervals();
        let raw = |i: usize, j: usize, k: usize| i + n1[0] * (j + n1[1] * k);
        // the cells of the removed boxes (by the centroid of the
        // hexahedron) and the nodes that the kept cells use
        let cell = |i: usize, j: usize, k: usize| i + n[0] * (j + n[1] * k);
        let mask: Vec<bool> = if self.remove.is_empty() {
            Vec::new()
        } else {
            let mut m = vec![true; n[0] * n[1] * n[2]];
            for k in 0..n[2] {
                for j in 0..n[1] {
                    for i in 0..n[0] {
                        let (p, q) = (pts[raw(i, j, k)], pts[raw(i + 1, j + 1, k + 1)]);
                        let c = [
                            0.5 * (p[0] + q[0]),
                            0.5 * (p[1] + q[1]),
                            0.5 * (p[2] + q[2]),
                        ];
                        m[cell(i, j, k)] =
                            !self.remove.iter().any(|r| r.contains(self.dimension, c));
                    }
                }
            }
            m
        };
        let keep = |i: usize, j: usize, k: usize| mask.is_empty() || mask[cell(i, j, k)];
        let mut map: Vec<u32> = (0..pts.len() as u32).collect();
        let pts = if self.remove.is_empty() {
            pts
        } else {
            let mut used = vec![false; pts.len()];
            for k in 0..n[2] {
                for j in 0..n[1] {
                    for i in 0..n[0] {
                        if keep(i, j, k) {
                            for (a, b, e) in (0..8).map(|m| (m & 1, (m >> 1) & 1, m >> 2)) {
                                used[raw(i + a, j + b, k + e)] = true;
                            }
                        }
                    }
                }
            }
            let mut kept = Vec::new();
            for (q, p) in pts.iter().enumerate() {
                if used[q] {
                    map[q] = kept.len() as u32;
                    kept.push(*p);
                } else {
                    map[q] = u32::MAX;
                }
            }
            if kept.is_empty() {
                return Err("remove: the boxes take out every cell".into());
            }
            kept
        };
        let node = |i: usize, j: usize, k: usize| map[raw(i, j, k)];
        let mut soup = CellSoup::new(pts);
        let x_band = |i: usize| -> Kind {
            match self.kind {
                Kind::Mixed if self.dimension == 3 => {
                    if 3 * i < n[0] {
                        Kind::Hex
                    } else if 3 * i < 2 * n[0] {
                        Kind::Mixed // prisms
                    } else {
                        Kind::Tet
                    }
                }
                Kind::Mixed => {
                    if 2 * i < n[0] {
                        Kind::Hex
                    } else {
                        Kind::Mixed
                    }
                }
                Kind::Tet if self.dimension == 2 => Kind::Mixed, // prisms
                k => k,
            }
        };
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    if !keep(i, j, k) {
                        continue;
                    }
                    let v = |a: usize, b: usize, c: usize| node(i + a, j + b, k + c);
                    match x_band(i) {
                        Kind::Hex => {
                            let h = [
                                v(0, 0, 0),
                                v(1, 0, 0),
                                v(1, 1, 0),
                                v(0, 1, 0),
                                v(0, 0, 1),
                                v(1, 0, 1),
                                v(1, 1, 1),
                                v(0, 1, 1),
                            ];
                            let faces = hex_faces(&h);
                            soup.add_cell(faces.iter().map(|f| &f[..]));
                        }
                        Kind::Mixed => {
                            // two prisms, split along the (0,0)-(1,1) diagonal
                            for (a, b) in [((1, 0), (1, 1)), ((1, 1), (0, 1))] {
                                let t0 = [v(0, 0, 0), v(a.0, a.1, 0), v(b.0, b.1, 0)];
                                let t1 = [v(0, 0, 1), v(a.0, a.1, 1), v(b.0, b.1, 1)];
                                let faces = prism_faces(&t0, &t1);
                                soup.add_cell(faces.iter().map(|f| &f[..]));
                            }
                        }
                        Kind::Tet => {
                            for perm in [
                                [0, 1, 2],
                                [0, 2, 1],
                                [1, 0, 2],
                                [1, 2, 0],
                                [2, 0, 1],
                                [2, 1, 0],
                            ] {
                                let mut c = [0usize; 3];
                                let mut t = [v(0, 0, 0); 4];
                                for (s, &ax) in perm.iter().enumerate() {
                                    c[ax] = 1;
                                    t[s + 1] = v(c[0], c[1], c[2]);
                                }
                                let faces = tet_faces(&t);
                                soup.add_cell(faces.iter().map(|f| &f[..]));
                            }
                        }
                        Kind::Poly => unreachable!(),
                    }
                }
            }
        }
        soup.finish(
            &self
                .boundary_names()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            self.classifier(),
        )
    }

    /// 2D polygonal prisms: the dual of the triangulation, extruded.
    fn poly2d(&self) -> Result<Mesh, String> {
        let (grid, n1) = self.grid();
        let n = self.intervals();
        let nv = n1[0] * n1[1];
        let v2 = |i: usize, j: usize| (i + n1[0] * j) as u32;
        let pt = |v: u32| grid[v as usize];
        // triangles: two per quadrilateral, along the (0,0)-(1,1) diagonal
        let mut tris: Vec<[u32; 3]> = Vec::with_capacity(2 * n[0] * n[1]);
        for j in 0..n[1] {
            for i in 0..n[0] {
                tris.push([v2(i, j), v2(i + 1, j), v2(i + 1, j + 1)]);
                tris.push([v2(i, j), v2(i + 1, j + 1), v2(i, j + 1)]);
            }
        }
        // 2D points: triangle centroids, boundary edge midpoints, boundary
        // vertices
        let mut p2: Vec<[f64; 2]> = Vec::new();
        let xy = |p: [f64; 3]| [p[0], p[1]];
        for t in &tris {
            let c = scale(add(add(pt(t[0]), pt(t[1])), pt(t[2])), 1.0 / 3.0);
            p2.push(xy(c));
        }
        let mut links: Vec<Vec<(u32, u32)>> = vec![Vec::new(); nv]; // (other vertex, triangle)
        for (ti, t) in tris.iter().enumerate() {
            for a in 0..3 {
                let v = t[a];
                links[v as usize].push((t[(a + 1) % 3], ti as u32));
                links[v as usize].push((t[(a + 2) % 3], ti as u32));
            }
        }
        let mut mid: HashMap<(u32, u32), u32> = HashMap::new();
        let mut vert: HashMap<u32, u32> = HashMap::new();
        let mut polys: Vec<Vec<u32>> = Vec::with_capacity(nv);
        for v in 0..nv as u32 {
            let ring = walk_ring(&links[v as usize]);
            let mut poly = Vec::new();
            if ring.open {
                let id = *vert.entry(v).or_insert_with(|| {
                    p2.push(xy(pt(v)));
                    (p2.len() - 1) as u32
                });
                poly.push(id);
                let mut push_mid = |x: u32, poly: &mut Vec<u32>| {
                    let key = (v.min(x), v.max(x));
                    let id = *mid.entry(key).or_insert_with(|| {
                        p2.push(xy(scale(add(pt(v), pt(x)), 0.5)));
                        (p2.len() - 1) as u32
                    });
                    poly.push(id);
                };
                push_mid(ring.nodes[0], &mut poly);
                for &t in &ring.items {
                    poly.push(t);
                }
                push_mid(*ring.nodes.last().unwrap(), &mut poly);
            } else {
                poly.extend_from_slice(&ring.items);
            }
            // counter-clockwise
            let mut a2 = 0.0;
            for i in 0..poly.len() {
                let (p, q) = (
                    p2[poly[i] as usize],
                    p2[poly[(i + 1) % poly.len()] as usize],
                );
                a2 += p[0] * q[1] - q[0] * p[1];
            }
            if a2 < 0.0 {
                poly.reverse();
            }
            polys.push(poly);
        }
        // extrude
        let np = p2.len() as u32;
        let (z0, z1) = (self.lower[2], self.upper[2]);
        let mut pts: Vec<[f64; 3]> = p2.iter().map(|p| [p[0], p[1], z0]).collect();
        pts.extend(p2.iter().map(|p| [p[0], p[1], z1]));
        let mut soup = CellSoup::new(pts);
        let mut faces: Vec<Vec<u32>> = Vec::new();
        for poly in &polys {
            faces.clear();
            faces.push(poly.iter().rev().copied().collect());
            faces.push(poly.iter().map(|&p| p + np).collect());
            for i in 0..poly.len() {
                let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
                faces.push(vec![a, b, b + np, a + np]);
            }
            soup.add_cell(faces.iter().map(|f| &f[..]));
        }
        soup.finish(
            &self
                .boundary_names()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            self.classifier(),
        )
    }

    /// 3D polyhedra: the dual of the Kuhn tetrahedral mesh.
    fn poly3d(&self) -> Result<Mesh, String> {
        let (grid, n1) = self.grid();
        let n = self.intervals();
        let node = |i: usize, j: usize, k: usize| (i + n1[0] * (j + n1[1] * k)) as u32;
        let mut tets: Vec<[u32; 4]> = Vec::with_capacity(6 * n[0] * n[1] * n[2]);
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    for perm in [
                        [0, 1, 2],
                        [0, 2, 1],
                        [1, 0, 2],
                        [1, 2, 0],
                        [2, 0, 1],
                        [2, 1, 0],
                    ] {
                        let mut c = [0usize; 3];
                        let mut t = [node(i, j, k); 4];
                        for (s, &ax) in perm.iter().enumerate() {
                            c[ax] = 1;
                            t[s + 1] = node(i + c[0], j + c[1], k + c[2]);
                        }
                        tets.push(t);
                    }
                }
            }
        }
        dual_of_tets(
            &grid,
            &tets,
            &self
                .boundary_names()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            self.classifier(),
        )
    }
}

/// Faces of a hexahedron with nodes in the usual order (bottom loop, then
/// the top loop above it).
pub fn hex_faces(h: &[u32; 8]) -> [[u32; 4]; 6] {
    [
        [h[0], h[3], h[2], h[1]],
        [h[4], h[5], h[6], h[7]],
        [h[0], h[1], h[5], h[4]],
        [h[1], h[2], h[6], h[5]],
        [h[2], h[3], h[7], h[6]],
        [h[3], h[0], h[4], h[7]],
    ]
}

/// Faces of a tetrahedron, consistently oriented (outward for a positive
/// orientation of the nodes).
pub fn tet_faces(t: &[u32; 4]) -> [[u32; 3]; 4] {
    [
        [t[0], t[2], t[1]],
        [t[0], t[1], t[3]],
        [t[0], t[3], t[2]],
        [t[1], t[2], t[3]],
    ]
}

fn prism_faces(t0: &[u32; 3], t1: &[u32; 3]) -> [Vec<u32>; 5] {
    [
        vec![t0[0], t0[2], t0[1]],
        vec![t1[0], t1[1], t1[2]],
        vec![t0[0], t0[1], t1[1], t1[0]],
        vec![t0[1], t0[2], t1[2], t1[1]],
        vec![t0[2], t0[0], t1[0], t1[2]],
    ]
}

/// A ring of items around a centre: `nodes[i]` and `nodes[i + 1]` are
/// joined by `items[i]`; for a closed ring `nodes.len() == items.len()`
/// (the last item joins the last and the first node), for an open one
/// `nodes.len() == items.len() + 1`.
#[derive(Debug)]
struct Ring {
    nodes: Vec<u32>,
    items: Vec<u32>,
    open: bool,
}

/// Orders links `(node, item)` (each item appears twice, with its two
/// nodes) into a ring.
fn walk_ring(links: &[(u32, u32)]) -> Ring {
    // node -> items
    let mut by_node: Vec<(u32, u32)> = links.to_vec();
    by_node.sort_unstable();
    let items_of = |v: u32| -> Vec<u32> {
        let lo = by_node.partition_point(|x| x.0 < v);
        by_node[lo..]
            .iter()
            .take_while(|x| x.0 == v)
            .map(|x| x.1)
            .collect()
    };
    let other = |item: u32, v: u32| -> u32 {
        links
            .iter()
            .find(|x| x.1 == item && x.0 != v)
            .map(|x| x.0)
            .expect("item with two nodes")
    };
    // start at a node with one item (open ring), else the smallest node
    let mut start = by_node[0].0;
    let mut open = false;
    let mut i = 0;
    while i < by_node.len() {
        let v = by_node[i].0;
        let cnt = by_node[i..].iter().take_while(|x| x.0 == v).count();
        if cnt == 1 {
            start = v;
            open = true;
            break;
        }
        i += cnt;
    }
    let n_items = links.len() / 2;
    let mut nodes = vec![start];
    let mut items = Vec::with_capacity(n_items);
    let mut prev_item = u32::MAX;
    let mut v = start;
    while items.len() < n_items {
        let its = items_of(v);
        let item = *its.iter().find(|&&t| t != prev_item).unwrap_or(&its[0]);
        if items.contains(&item) {
            break;
        }
        items.push(item);
        let w = other(item, v);
        prev_item = item;
        v = w;
        if items.len() < n_items || open {
            nodes.push(w);
        }
    }
    Ring { nodes, items, open }
}

/// The polyhedral mesh of a tetrahedral mesh (for example from Gmsh): the
/// dual of [`dual_of_tets`], with the boundaries of the boundary triangles.
/// Every cell must be a tetrahedron. Nodes of the tetrahedral mesh that no
/// cell uses are not allowed (the Gmsh reader drops them).
pub fn polyhedral_from_tets(mesh: &Mesh) -> Result<Mesh, String> {
    let mut tets = Vec::with_capacity(mesh.n_cells);
    for c in 0..mesh.n_cells {
        if mesh.cell_shape(c) != crate::mesh::CellShape::Tetrahedron {
            return Err(format!(
                "cell {c} is a {}: only tetrahedral meshes are converted to polyhedra",
                mesh.cell_shape(c).name()
            ));
        }
        let v = mesh.cell_nodes(c);
        tets.push([v[0], v[1], v[2], v[3]]);
    }
    let mut boundary_of: HashMap<[u32; 3], usize> = HashMap::new();
    for (i, p) in mesh.boundaries.iter().enumerate() {
        for f in p.range() {
            let n = mesh.face(f);
            let mut k = [n[0], n[1], n[2]];
            k.sort_unstable();
            boundary_of.insert(k, i);
        }
    }
    let names: Vec<&str> = mesh.boundaries.iter().map(|p| p.name.as_str()).collect();
    let mut m = dual_of_tets(&mesh.points, &tets, &names, |b| {
        let mut k = [b.nodes[0], b.nodes[1], b.nodes[2]];
        k.sort_unstable();
        names[boundary_of[&k]].to_string()
    })?;
    // zones: the zone of the tetrahedra around a node, if they agree
    let mut zone = vec![u32::MAX; m.n_cells];
    for (c, t) in tets.iter().enumerate() {
        for &v in t {
            let z = &mut zone[v as usize];
            *z = if *z == u32::MAX || *z == mesh.cell_zone[c] {
                mesh.cell_zone[c]
            } else {
                0
            };
        }
    }
    m.zones = mesh.zones.clone();
    m.cell_zone = zone
        .into_iter()
        .map(|z| if z == u32::MAX { 0 } else { z })
        .collect();
    Ok(m)
}

/// The polyhedral dual of a tetrahedral mesh: one cell per node.
///
/// Dual nodes: tetrahedron centroids, triangle centroids, midpoints of
/// boundary edges and the boundary nodes themselves. The face between the
/// cells of nodes a and b runs through the centroids of the triangles and
/// tetrahedra around edge (a, b), closed by the edge midpoint if the edge
/// lies on the boundary. The boundary faces of a cell are the pieces of the
/// boundary triangles around its node, merged into one polygon per boundary
/// where they are coplanar and kept as one planar piece per triangle on
/// curved boundaries, so the dual has exactly the volume of the tet mesh.
pub fn dual_of_tets(
    points: &[[f64; 3]],
    tets: &[[u32; 4]],
    boundaries: &[&str],
    classify: impl Fn(BoundaryFace) -> String,
) -> Result<Mesh, String> {
    let nv = points.len();
    let pt = |v: u32| points[v as usize];
    // triangles with the tetrahedra on either side
    let mut tri_list: Vec<([u32; 3], u32)> = Vec::with_capacity(4 * tets.len());
    for (ti, t) in tets.iter().enumerate() {
        for skip in 0..4 {
            let mut k = [0u32; 3];
            let mut m = 0;
            for (a, &v) in t.iter().enumerate() {
                if a != skip {
                    k[m] = v;
                    m += 1;
                }
            }
            k.sort_unstable();
            tri_list.push((k, ti as u32));
        }
    }
    tri_list.sort_unstable();
    let mut tri_keys: Vec<[u32; 3]> = Vec::new();
    let mut tri_tets: Vec<[u32; 2]> = Vec::new();
    let mut i = 0;
    while i < tri_list.len() {
        let k = tri_list[i].0;
        if i + 1 < tri_list.len() && tri_list[i + 1].0 == k {
            tri_keys.push(k);
            tri_tets.push([tri_list[i].1, tri_list[i + 1].1]);
            i += 2;
        } else {
            tri_keys.push(k);
            tri_tets.push([tri_list[i].1, u32::MAX]);
            i += 1;
        }
    }
    drop(tri_list);
    let tri_index = |k: [u32; 3]| -> u32 {
        tri_keys
            .binary_search(&k)
            .expect("triangle of the tetrahedral mesh") as u32
    };
    let tet_centroid =
        |t: &[u32; 4]| scale(add(add(pt(t[0]), pt(t[1])), add(pt(t[2]), pt(t[3]))), 0.25);
    // dual points
    let n_tet = tets.len() as u32;
    let mut dp: Vec<[f64; 3]> = tets.iter().map(tet_centroid).collect();
    dp.extend(
        tri_keys
            .iter()
            .map(|k| scale(add(add(pt(k[0]), pt(k[1])), pt(k[2])), 1.0 / 3.0)),
    );
    let tri_pt = |t: u32| n_tet + t;
    // boundary triangles, outward orientation and boundary
    let mut btri: Vec<(u32, [u32; 3], String)> = Vec::new(); // (tri, outward nodes, boundary)
    for (t, (k, tt)) in tri_keys.iter().zip(&tri_tets).enumerate() {
        if tt[1] != u32::MAX {
            continue;
        }
        let tc = tet_centroid(&tets[tt[0] as usize]);
        let (s, c) = polygon_area_centroid(points, k);
        let mut nodes = *k;
        let mut s = s;
        if dot(s, sub(c, tc)) < 0.0 {
            nodes.swap(1, 2);
            s = scale(s, -1.0);
        }
        let ns = crate::geometry::norm(s);
        let name = classify(BoundaryFace {
            nodes: &nodes,
            centroid: c,
            normal: scale(s, 1.0 / ns),
        });
        btri.push((t as u32, nodes, name));
    }
    // boundary edge midpoints and boundary vertices
    let mut mid: HashMap<(u32, u32), u32> = HashMap::new();
    let mut bvert: HashMap<u32, u32> = HashMap::new();
    for (_, k, _) in &btri {
        for a in 0..3 {
            let (u, w) = (k[a].min(k[(a + 1) % 3]), k[a].max(k[(a + 1) % 3]));
            mid.entry((u, w)).or_insert_with(|| {
                dp.push(scale(add(pt(u), pt(w)), 0.5));
                (dp.len() - 1) as u32
            });
            bvert.entry(k[a]).or_insert_with(|| {
                dp.push(pt(k[a]));
                (dp.len() - 1) as u32
            });
        }
    }
    // edges with their tetrahedra
    let mut edge_list: Vec<(u32, u32, u32)> = Vec::with_capacity(6 * tets.len());
    for (ti, t) in tets.iter().enumerate() {
        for a in 0..4 {
            for b in a + 1..4 {
                edge_list.push((t[a].min(t[b]), t[a].max(t[b]), ti as u32));
            }
        }
    }
    edge_list.sort_unstable();
    let mut builder = MeshBuilder::new(dp.clone(), nv);
    for p in boundaries {
        builder.declare_boundary(p);
    }
    let mut links: Vec<(u32, u32)> = Vec::new();
    let mut poly: Vec<u32> = Vec::new();
    let mut i = 0;
    while i < edge_list.len() {
        let (a, b) = (edge_list[i].0, edge_list[i].1);
        let mut j = i;
        links.clear();
        while j < edge_list.len() && edge_list[j].0 == a && edge_list[j].1 == b {
            let t = &tets[edge_list[j].2 as usize];
            for &v in t {
                if v != a && v != b {
                    links.push((v, edge_list[j].2));
                }
            }
            j += 1;
        }
        let ring = walk_ring(&links);
        poly.clear();
        let tri_of = |v: u32| {
            let mut k = [a, b, v];
            k.sort_unstable();
            tri_pt(tri_index(k))
        };
        if ring.open {
            poly.push(mid[&(a, b)]);
        }
        for (m, &t) in ring.items.iter().enumerate() {
            poly.push(tri_of(ring.nodes[m]));
            poly.push(t);
        }
        if ring.open {
            poly.push(tri_of(*ring.nodes.last().unwrap()));
        }
        let (s, _) = polygon_area_centroid(&dp, &poly);
        if dot(s, sub(pt(b), pt(a))) < 0.0 {
            poly.reverse();
        }
        builder.add_interior(&poly, a, b);
        i = j;
    }
    // boundary faces: group boundary triangles by (vertex, boundary)
    let mut by_vertex: Vec<(u32, String, usize)> = Vec::with_capacity(3 * btri.len());
    for (bi, (_, k, name)) in btri.iter().enumerate() {
        for &v in k {
            by_vertex.push((v, name.clone(), bi));
        }
    }
    by_vertex.sort();
    let mut i = 0;
    while i < by_vertex.len() {
        let (v, name) = (by_vertex[i].0, by_vertex[i].1.clone());
        let mut j = i;
        links.clear();
        let mut normal = [0.0; 3];
        let mut unit: Option<[f64; 3]> = None;
        let mut planar = true;
        while j < by_vertex.len() && by_vertex[j].0 == v && by_vertex[j].1 == name {
            let (t, k, _) = &btri[by_vertex[j].2];
            for &w in k {
                if w != v {
                    links.push((w, *t));
                }
            }
            let s = polygon_area_centroid(points, k).0;
            let u = scale(s, 1.0 / crate::geometry::norm(s));
            match unit {
                None => unit = Some(u),
                Some(u0) => planar &= dot(u0, u) > 1.0 - 1e-12,
            }
            normal = add(normal, s);
            j += 1;
        }
        let m = |w: u32| mid[&(v.min(w), v.max(w))];
        if !planar {
            // curved boundary: one planar piece per triangle, so the dual
            // cells fill exactly the volume of the tetrahedral mesh
            for q in i..j {
                let (t, k, _) = &btri[by_vertex[q].2];
                let a = k.iter().position(|&w| w == v).unwrap();
                let (x, y) = (k[(a + 1) % 3], k[(a + 2) % 3]);
                let quad = [bvert[&v], m(x), tri_pt(*t), m(y)];
                builder.add_boundary(&quad, v, &name);
            }
            i = j;
            continue;
        }
        let ring = walk_ring(&links);
        poly.clear();
        if ring.open {
            poly.push(bvert[&v]);
        }
        for (q, &t) in ring.items.iter().enumerate() {
            poly.push(m(ring.nodes[q]));
            poly.push(tri_pt(t));
        }
        if ring.open {
            poly.push(m(*ring.nodes.last().unwrap()));
        }
        let (s, _) = polygon_area_centroid(&dp, &poly);
        if dot(s, normal) < 0.0 {
            poly.reverse();
        }
        builder.add_boundary(&poly, v, &name);
        i = j;
    }
    let mesh = builder.finish()?;
    debug_assert_eq!(mesh.points.len(), dp.len());
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grading_is_nested_and_has_the_ratio() {
        let g = Grading {
            ratio: 8.0,
            symmetric: false,
        };
        let x0 = g.nodes(10, 0);
        let h: Vec<f64> = x0.windows(2).map(|w| w[1] - w[0]).collect();
        assert!((h[9] / h[0] - 8.0).abs() < 1e-12);
        let x1 = g.nodes(20, 1);
        for i in 0..=10 {
            assert!((x1[2 * i] - x0[i]).abs() < 1e-14);
        }
        let s = Grading {
            ratio: 5.0,
            symmetric: true,
        };
        let x = s.nodes(12, 0);
        let h: Vec<f64> = x.windows(2).map(|w| w[1] - w[0]).collect();
        assert!((h[5] / h[0] - 5.0).abs() < 1e-12);
        assert!((h[0] - h[11]).abs() < 1e-14);
    }

    /// The TMR plate mesh: x segments meet at the leading edge with the
    /// requested end ratios, and the bottom boundary splits there.
    #[test]
    fn segments_and_split_build_the_plate_mesh() {
        let spec = BoxSpec {
            dimension: 2,
            lower: [-1.0 / 3.0, 0.0, 0.0],
            upper: [2.0, 1.0, 1.0],
            cells: [32, 6, 1],
            segments: SegmentSpec {
                x: vec![
                    Segment {
                        to: 0.0,
                        cells: 8,
                        ratio: 0.1,
                    },
                    Segment {
                        to: 2.0,
                        cells: 24,
                        ratio: 20.0,
                    },
                ],
                ..Default::default()
            },
            split: vec![Split {
                boundary: "ymin".into(),
                axis: "x".into(),
                at: 0.0,
                name: "plate".into(),
            }],
            ..Default::default()
        };
        for refine in [0, 1] {
            let mesh = BoxSpec {
                refine,
                ..spec.clone()
            }
            .generate()
            .unwrap();
            let geo = crate::Geometry::compute(&mesh);
            let n = 32usize << refine;
            assert_eq!(mesh.n_cells, n * (6 << refine));
            let names: Vec<&str> = mesh.boundaries.iter().map(|p| p.name.as_str()).collect();
            assert!(
                names.contains(&"ymin") && names.contains(&"plate"),
                "{names:?}"
            );
            let side = |name: &str| {
                let p = mesh.boundaries.iter().find(|p| p.name == name).unwrap();
                p.range()
                    .map(|f| geo.face_centroid[f][0])
                    .collect::<Vec<f64>>()
            };
            let (ahead, plate) = (side("ymin"), side("plate"));
            assert_eq!(ahead.len(), 8 << refine);
            assert_eq!(plate.len(), 24 << refine);
            assert!(ahead.iter().all(|&x| x < 0.0) && plate.iter().all(|&x| x > 0.0));
            // the plate's first and last intervals: ratio 20 (nested halving
            // keeps the end ratio of the base mesh within one level's growth)
            let mut xs = plate.clone();
            xs.sort_by(f64::total_cmp);
            let h0 = 2.0 * xs[0];
            let h1 = 2.0 * (2.0 - xs[xs.len() - 1]);
            assert!(h1 > h0 * 10.0, "first {h0}, last {h1}");
        }
    }

    /// A removed box makes a backward-facing step: the cells under the
    /// step corner go, their exposed faces form the boundary `step`, the
    /// unused nodes are dropped and the volume is that of the L shape.
    #[test]
    fn removed_box_builds_a_step() {
        for (kind, per_hex) in [(Kind::Hex, 1), (Kind::Tet, 2), (Kind::Mixed, 0)] {
            let spec = BoxSpec {
                kind,
                dimension: 2,
                lower: [-2.0, 0.0, 0.0],
                upper: [4.0, 3.0, 0.5],
                cells: [6, 3, 1],
                remove: vec![RemoveBox {
                    lower: [-2.0, 0.0, 0.0],
                    upper: [0.0, 1.0, 0.5],
                    name: "step".into(),
                }],
                ..Default::default()
            };
            let mesh = spec.generate().unwrap();
            mesh.validate().unwrap();
            if per_hex > 0 {
                assert_eq!(mesh.n_cells, per_hex * (18 - 2));
            }
            let geo = crate::Geometry::compute(&mesh);
            let vol: f64 = geo.cell_volume.iter().sum();
            assert!((vol - (18.0 - 2.0) * 0.5).abs() < 1e-12, "{kind:?}: {vol}");
            let p = mesh.boundaries.iter().find(|p| p.name == "step").unwrap();
            // the top of the step (2 faces, y = 1) and its face (1, x = 0)
            let area: f64 = p
                .range()
                .map(|f| dot(geo.face_area[f], geo.face_area[f]).sqrt())
                .sum();
            assert!((area - 3.0 * 0.5).abs() < 1e-12, "{kind:?}: {area}");
            for f in p.range() {
                let c = geo.face_centroid[f];
                assert!(
                    (c[1] - 1.0).abs() < 1e-12 && c[0] < 0.0 || c[0].abs() < 1e-12 && c[1] < 1.0,
                    "{c:?}"
                );
            }
            let ymin = mesh.boundaries.iter().find(|p| p.name == "ymin").unwrap();
            assert!(ymin.range().all(|f| geo.face_centroid[f][0] > 0.0));
            // every node is used
            let mut used = vec![false; mesh.points.len()];
            for &n in &mesh.face_nodes {
                used[n as usize] = true;
            }
            assert!(used.iter().all(|&u| u), "{kind:?}: unused nodes");
        }
        let bad = BoxSpec {
            remove: vec![RemoveBox {
                lower: [0.0; 3],
                upper: [0.5; 3],
                name: "xmin".into(),
            }],
            ..Default::default()
        };
        assert!(bad.validate().unwrap_err().contains("taken"));
    }
}
