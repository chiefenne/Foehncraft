//! Immersed bodies by Brinkman penalisation (L13 stage 1, D226).
//!
//! A body is put into the flow without meshing it: every cell gets its
//! solid fraction χ_P ∈ [0, 1] (the part of its volume inside the body), and
//! the momentum rows get the Brinkman drag (Angot, Bruneau & Fabrie 1999)
//!
//! −(μ / K) χ_P u_P V_P = −s_P u_P,
//!
//! with K a permeability (m²): as K → 0 the velocity in the solid goes to
//! zero and the body acts as a no-slip wall smeared over about one cell.
//! The sink s_P is implicit: it is added to the diagonal a_P of the momentum
//! row before the under-relaxation, inside the shared kernel
//! (`foehncraft_kernels::flow::momentum_row`), so that D_P = V_P / a_P of
//! the Rhie–Chow flux and the SIMPLEC factor see it as well. Otherwise the
//! pressure decouples in the solid. The CPU and GPU paths run the same
//! kernel; the coupled solver takes the same momentum rows.
//!
//! The force of the fluid on a body is the sum of the sink over its cells,
//! F = Σ s_P u_P (the momentum the penalisation takes out of the fluid).
//!
//! **What this is and what it is not (D226, owner decisions Q5, Q9).** A
//! topology search and a cheap first look: laminar and low-Re only. The wall
//! is smeared over about one cell, the wall shear stress and heat flux are
//! not accurate, a wall function does not fit (refused, IMM102), the
//! result depends on K and on the mesh width, the stiffness grows as K → 0,
//! and thin walls leak. Sharp walls come with the cut cells of stage 2.
//!
//! **Solid fraction.** Each cell is split into tetrahedra (cell centroid,
//! face centroid, one face edge); each tetrahedron is sampled with the
//! midpoint rule of `samples³` points in collapsed coordinates (Duffy), each
//! point weighted with its Jacobian. The inside test of the primitive
//! shapes is exact; the fraction converges with the number of samples. The
//! result is deterministic and independent of threads and ranks (each cell
//! by itself).
//!
//! Case file:
//!
//! ```toml
//! [immersed]
//! permeability = 1e-8        # K (m²); the sink is μ χ V / K
//! samples = 4                # per direction and tetrahedron (default 4)
//!
//! [[immersed.body]]
//! name = "cylinder"
//! shape = "cylinder"         # "cylinder" | "sphere" | "box"
//! centre = [0.2, 0.2, 0.0]   # cylinder: a point on the axis; sphere: centre
//! radius = 0.05
//! axis = [0.0, 0.0, 1.0]     # cylinder axis (default z)
//! # length = 0.1             # cylinder length centred on `centre`; default endless
//! # lower = [..], upper = [..]   # box corners
//!
//! [[immersed.body]]
//! name = "part"
//! shape = "stl"              # a triangle surface (L13 stage 2)
//! file = "part.stl"          # ASCII or binary, relative to the case file
//! scale = 1e-3               # factor on the coordinates (a file in mm)
//! ```
//!
//! **STL bodies.** The surface is checked when the case is read (IMM001 to
//! IMM009 of stage 0: holes, flips, self-intersections; nothing is
//! repaired) and its triangles are the body: inside is where the
//! generalised winding number exceeds ½. A penalised STL body gets the
//! exact solid fraction of every cell from the cut moments of stage 0
//! (`Surface::cut_moments`, no sampling); a cut STL body is cut with the
//! exact segment–triangle crossings, the exact closest point and the
//! dual-contouring point on its sharp edges ([`foehncraft_mesh::surface::SurfaceBody`]).

use crate::diag::{self, Diagnostic};
use anyhow::Result;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// The shape of an immersed body.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Shape {
    /// A circular cylinder about `axis` through `centre`, endless or of
    /// `length`.
    Cylinder,
    /// A sphere.
    Sphere,
    /// An axis-aligned box from `lower` to `upper`.
    Box,
    /// A solid fraction given as a field χ(x, y, z) ∈ [0, 1] (an
    /// expression, D53): a smooth density, as in the manufactured solution
    /// of the penalised equations or a design density.
    Field,
    /// A closed triangle surface from an STL file (`file`, `scale`).
    Stl,
    /// A half-space: the solid lies behind the plane through `centre`
    /// with the unit `normal` pointing into the fluid (an inclined wall
    /// across the mesh, L13 stage 2c).
    Plane,
}

/// One `[[immersed.body]]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BodySection {
    /// Name of the body (in the report, the solid-fraction field and as the
    /// `[verification] benchmark_boundary`).
    pub name: String,
    /// Shape.
    pub shape: Shape,
    /// Cylinder: a point on the axis (the middle of a finite one); sphere:
    /// the centre (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub centre: Option<[f64; 3]>,
    /// Radius of a cylinder or sphere (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub radius: Option<f64>,
    /// Axis of a cylinder (default z).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub axis: Option<[f64; 3]>,
    /// Normal of a `plane` body, pointing out of the solid into the fluid.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub normal: Option<[f64; 3]>,
    /// Length of a cylinder, centred on `centre` (m); default endless.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub length: Option<f64>,
    /// Lower corner of a box (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lower: Option<[f64; 3]>,
    /// Upper corner of a box (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper: Option<[f64; 3]>,
    /// The solid fraction χ(x, y, z) of a `field` body (an expression of
    /// x, y, z and the case constants, clipped to [0, 1]).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fraction: Option<crate::expr::Scalar>,
    /// The STL file of an `stl` body (ASCII or binary), relative to the
    /// case file.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file: Option<std::path::PathBuf>,
    /// Factor on the coordinates of the STL file (default 1; 1e-3 for a
    /// file in millimetres).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scale: Option<f64>,
    /// The fluid lies inside the shape and the solid outside it: a pipe or
    /// the outer cylinder of a gap cut out of a box mesh (cylinder, sphere
    /// and box; L13 stage 2e).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub invert: Option<bool>,
    /// The body turns rigidly with this angular velocity (rad/s, right-hand
    /// about `axis` through `centre`; a box turns about its middle when
    /// `centre` is not given): its wall faces move with ω × (x − centre),
    /// and the summary of each wall of the body gives the torque about the
    /// axis. 0 gives the torque of a body at rest. Cut bodies only (L13
    /// stage 2e).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub omega: Option<f64>,
    /// Named surface patches of an `stl` body (L13 stage 2e): `"solids"`
    /// takes the solids of an ASCII file, a table `{ <patch> = <tag> }` the
    /// 16-bit attributes of a binary file. The wall faces of the cut become
    /// one boundary per patch, `<body>-<patch>`, each with any condition of
    /// `[boundary.<name>]`; faces of triangles of no listed tag stay in the
    /// boundary `<body>`. A face takes the patch of the triangle closest to
    /// its centroid.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub patches: Option<Patches>,
}

/// The named surface patches of an `stl` body.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum Patches {
    /// `"solids"`: one patch per solid of an ASCII file, named after it.
    Solids(PatchSource),
    /// The patches of a binary file: patch name to attribute value.
    Tags(std::collections::BTreeMap<String, u16>),
}

/// Where the patches of an `stl` body come from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum PatchSource {
    /// The solids of an ASCII STL file.
    Solids,
}

/// How the bodies enter the flow.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Method {
    /// Brinkman penalisation on the uncut mesh (stage 1, D226).
    #[default]
    Penalisation,
    /// Conservative cut cells: the bodies are cut out of the mesh and their
    /// surfaces become wall boundaries named after them (stage 2, D227).
    Cut,
}

/// How a cut puts the crossing point on a mesh edge.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum EdgeRoot {
    /// The root of the body's own function on the edge: the point lies on
    /// the true surface.
    #[default]
    Exact,
    /// The zero of the linear interpolant of the node values (marching
    /// cubes).
    Linear,
}

/// `[immersed.cut]`: the settings of the cut cells.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CutSection {
    /// Cut cells below this fraction of their uncut volume are merged into
    /// a neighbour (0 to 0.9, default 0.3).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_volume_fraction: Option<f64>,
    /// The crossing point on a mesh edge (default `exact`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub edge_root: Option<EdgeRoot>,
    /// Fan every wall face around a point on the surface, so that it
    /// follows a curved body inside the cell, and a sharp edge or corner
    /// of an STL body (default true).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_point: Option<bool>,
    /// Bend every straight cut across a face through a point on the
    /// surface in the face plane (or where a sharp edge pierces the face),
    /// so that the clipped faces follow the body too (default true).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub face_point: Option<bool>,
    /// The wall function of a cut cell takes its velocity from the cell one
    /// layer out along the wall normal, whose distance does not depend on
    /// how the cell was cut (L13 stage 2c, default true); false: from the
    /// cut cell itself.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_sample: Option<bool>,
}

impl CutSection {
    /// The merge threshold.
    pub fn min_volume_fraction(&self) -> f64 {
        self.min_volume_fraction.unwrap_or(0.3)
    }

    /// Bent face chords (default true).
    pub fn face_point(&self) -> bool {
        self.face_point.unwrap_or(true)
    }

    /// Sampled wall functions on cut walls (default true).
    pub fn wall_sample(&self) -> bool {
        self.wall_sample.unwrap_or(true)
    }
}

/// `[immersed]`: bodies in the flow by Brinkman penalisation or cut cells.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImmersedSection {
    /// `penalisation` (default) or `cut`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub method: Option<Method>,
    /// The settings of the cut cells, `[immersed.cut]`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cut: Option<CutSection>,
    /// Permeability K (m²) of the Brinkman term μ χ / K; smaller is closer
    /// to a no-slip wall and stiffer. No default: the error depends on it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub permeability: Option<f64>,
    /// Sample points per direction and tetrahedron of the solid fraction
    /// (1 to 16, default 4).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub samples: Option<u32>,
    /// The bodies.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub body: Vec<BodySection>,
    /// A density-field topology search (D226): `[immersed.topology]`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub topology: Option<crate::topology::TopologySection>,
    /// Refinement of the hexahedra near the bodies before the cut (L13
    /// stage 2f), `[immersed.refine]`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub refine: Option<RefineSection>,
}

/// `[immersed.refine]`: the hexahedra of the mesh near the surfaces of the
/// bodies are halved, up to `levels` times, into hexahedra with
/// hanging-node faces (L13 stage 2f, D57); 2:1 balanced. In a planar mesh
/// (one cell thick) the cells are halved in the plane only.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RefineSection {
    /// Halvings of a cell next to a surface (1 to 6).
    pub levels: u32,
    /// Width of each level around the surfaces, in cells of that level
    /// (default 4): a cell is halved while it comes closer to a surface
    /// than `layers` times the size of its halves.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub layers: Option<f64>,
    /// Width (m) of the finest level around the surfaces, instead of
    /// `layers`: a cell of level l is halved while it comes closer to a
    /// surface than distance · 2^(levels − 1 − l), so every level is twice
    /// as wide as the next finer one, the same in metres on any mesh (a
    /// convergence study refines the same region).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub distance: Option<f64>,
}

/// Refines the hexahedra near the bodies (L13 stage 2f): a cell is halved
/// while the distance of its centre to the nearest body surface, less its
/// radius, is below the width of its level (`layers` times half its size,
/// or from `distance`).
pub fn refine_mesh(case: &crate::case::FlowCase, mesh: &Mesh) -> Result<(Mesh, Vec<Diagnostic>)> {
    let Some(rs) = &case.immersed.refine else {
        return Ok((mesh.clone(), Vec::new()));
    };
    if rs.layers.is_some() && rs.distance.is_some() {
        return Err(bad(
            "immersed.refine".into(),
            "give `layers` or `distance`, not both".into(),
        ));
    }
    if let Some(d) = rs.distance
        && !(d.is_finite() && d > 0.0)
    {
        return Err(bad(
            "immersed.refine.distance".into(),
            format!("`distance` must be positive, not {d}"),
        ));
    }
    let layers = rs.layers.unwrap_or(4.0);
    if !(1..=foehncraft_mesh::refine::MAX_LEVELS).contains(&rs.levels) {
        return Err(bad(
            "immersed.refine.levels".into(),
            format!(
                "`levels` must be 1 to {}, not {}",
                foehncraft_mesh::refine::MAX_LEVELS,
                rs.levels
            ),
        ));
    }
    if !(layers.is_finite() && layers >= 1.0) {
        return Err(bad(
            "immersed.refine.layers".into(),
            format!("`layers` must be at least 1, not {layers}"),
        ));
    }
    let bs = bodies_of(case)?;
    if bs.is_empty() {
        return Err(bad(
            "immersed.refine".into(),
            "near-body refinement needs `[[immersed.body]]` bodies".into(),
        ));
    }
    // a planar mesh: every node on one of two planes z = const
    let (lo, hi) = mesh.bounding_box();
    let tol = 1e-9 * norm(sub(hi, lo));
    let planar = mesh
        .points
        .iter()
        .all(|p| (p[2] - lo[2]).abs() <= tol || (p[2] - hi[2]).abs() <= tol);
    let opt = foehncraft_mesh::refine::RefineOptions {
        levels: rs.levels,
        keep_direction: planar.then_some([0.0, 0.0, 1.0]),
    };
    let r = foehncraft_mesh::refine::refine(mesh, &opt, |c| {
        let d = bs
            .iter()
            .map(|b| b.distance(c.centre).abs())
            .fold(f64::INFINITY, f64::min);
        let band = match rs.distance {
            Some(w) => w * f64::powi(2.0, rs.levels as i32 - 1 - c.level as i32),
            None => layers * 0.5 * c.size,
        };
        d - c.radius < band
    })
    .map_err(|e| {
        bad(
            "immersed.refine".into(),
            format!("near-body refinement: {e}"),
        )
    })?;
    let geo = foehncraft_mesh::Geometry::compute(&r.mesh);
    let (q, _) = foehncraft_mesh::quality::QualityReport::compute(&r.mesh, &geo);
    let rep = &r.report;
    let levels: Vec<String> = rep
        .per_level
        .iter()
        .enumerate()
        .map(|(l, n)| format!("level {l}: {n}"))
        .collect();
    let note = Diagnostic::new(
        diag::IMM_REFINE_REPORT,
        format!(
            "near-body refinement: {} cells into {} ({}; {} input cells refined, {} more by the \
             2:1 balance, {} not hexahedra kept whole), {} faces between levels, {} faces with \
             hanging nodes; max non-orthogonality {:.1}°, max skewness {:.3}, max closure {:.2e}",
            rep.cells_in,
            rep.cells_out,
            levels.join(", "),
            rep.refined,
            rep.balance,
            rep.not_hexahedra,
            rep.hanging_faces,
            rep.faces_with_hanging_nodes,
            q.non_orthogonality.max,
            q.skewness.max,
            q.closure.max,
        ),
    );
    Ok((r.mesh, vec![note]))
}

/// Resolves the STL files of the bodies against the case directory.
pub fn resolve_paths(s: &mut ImmersedSection, dir: &std::path::Path) {
    for b in &mut s.body {
        if let Some(f) = &mut b.file {
            *f = dir.join(&*f);
        }
    }
}

impl ImmersedSection {
    /// True without an `[immersed]` table.
    pub fn is_empty(&self) -> bool {
        self == &Self::default()
    }

    /// True with `method = "cut"`.
    pub fn is_cut(&self) -> bool {
        self.method == Some(Method::Cut)
    }

    /// Samples per direction.
    pub fn samples(&self) -> usize {
        self.samples.unwrap_or(4) as usize
    }
}

fn bad(key: String, msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::IMM_CONFIG, msg).key(key))
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

/// A checked body with its inside test.
#[derive(Debug, Clone)]
pub struct Body {
    /// Name.
    pub name: String,
    kind: Kind,
    /// Fluid inside, solid outside (`invert`).
    invert: bool,
    /// Rigid rotation of the wall (`omega`): ω, the unit axis and a point on
    /// it.
    pub rotation: Option<Rotation>,
    /// The surface patches of an STL body: the patch names (without the
    /// body prefix) and, per triangle, its patch (None: the body itself).
    patches: Option<(Vec<String>, Vec<Option<u32>>)>,
    /// The surface check of an STL body (IMM001 to IMM009).
    pub findings: Vec<foehncraft_mesh::surface::check::Finding>,
}

/// The rigid rotation of a cut body's wall (`omega`, L13 stage 2e).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rotation {
    /// Angular velocity (rad/s), right-hand about `axis`.
    pub omega: f64,
    /// Unit axis.
    pub axis: [f64; 3],
    /// A point on the axis (m).
    pub origin: [f64; 3],
}

impl Rotation {
    /// The wall velocity ω a × (x − o) at `x`.
    pub fn velocity(&self, x: [f64; 3]) -> [f64; 3] {
        let w = cross(self.axis, sub(x, self.origin));
        [self.omega * w[0], self.omega * w[1], self.omega * w[2]]
    }
}

/// An STL body: the surface and what its check found.
#[derive(Debug, Clone)]
struct StlBody {
    surface: std::sync::Arc<foehncraft_mesh::surface::Surface>,
    /// The tag of every triangle and the solid names of an ASCII file.
    tags: std::sync::Arc<(Vec<u32>, Vec<String>)>,
    /// Closed, oriented outward and free of self-intersections: its
    /// enclosed volume and area are exact references.
    clean: bool,
    volume: f64,
    area: f64,
}

#[derive(Debug, Clone)]
enum Kind {
    Cylinder {
        c: [f64; 3],
        a: [f64; 3],
        r2: f64,
        half: f64,
    },
    Sphere {
        c: [f64; 3],
        r2: f64,
    },
    Box {
        lo: [f64; 3],
        hi: [f64; 3],
    },
    Field(crate::expr::Value),
    Stl(StlBody),
    Plane {
        p: [f64; 3],
        n: [f64; 3],
    },
}

impl Body {
    /// The surface of an STL body.
    pub fn surface(&self) -> Option<&foehncraft_mesh::surface::Surface> {
        match &self.kind {
            Kind::Stl(s) => Some(&s.surface),
            _ => None,
        }
    }

    /// The solid fraction at `x`: 1 inside a primitive body (the surface
    /// counts as inside), 0 outside; the clipped field of a `field` body.
    #[inline]
    pub fn value(&self, x: [f64; 3]) -> f64 {
        match &self.kind {
            Kind::Field(v) => {
                let f = v.at(x, 0.0);
                if f.is_nan() { 0.0 } else { f.clamp(0.0, 1.0) }
            }
            _ => {
                if self.inside(x) {
                    1.0
                } else {
                    0.0
                }
            }
        }
    }

    /// The signed distance to the surface, negative inside (exact for the
    /// primitive shapes); for a `field` body the level ½ − χ.
    pub fn distance(&self, x: [f64; 3]) -> f64 {
        let d = self.shape_distance(x);
        if self.invert { -d } else { d }
    }

    /// The signed distance of the shape, before `invert`.
    fn shape_distance(&self, x: [f64; 3]) -> f64 {
        match &self.kind {
            Kind::Cylinder { c, a, r2, half } => {
                let d = sub(x, *c);
                let t = dot(d, *a);
                let p = [d[0] - t * a[0], d[1] - t * a[1], d[2] - t * a[2]];
                let dr = dot(p, p).sqrt() - r2.sqrt();
                if half.is_infinite() {
                    return dr;
                }
                let dz = t.abs() - half;
                let (er, ez) = (dr.max(0.0), dz.max(0.0));
                (er * er + ez * ez).sqrt() + dr.max(dz).min(0.0)
            }
            Kind::Sphere { c, r2 } => {
                let d = sub(x, *c);
                dot(d, d).sqrt() - r2.sqrt()
            }
            Kind::Box { lo, hi } => {
                let mut q = [0.0; 3];
                for i in 0..3 {
                    let m = 0.5 * (lo[i] + hi[i]);
                    q[i] = (x[i] - m).abs() - 0.5 * (hi[i] - lo[i]);
                }
                let o = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
                dot(o, o).sqrt() + q[0].max(q[1]).max(q[2]).min(0.0)
            }
            Kind::Field(v) => {
                let f = v.at(x, 0.0);
                if f.is_nan() {
                    0.5
                } else {
                    0.5 - f.clamp(0.0, 1.0)
                }
            }
            Kind::Stl(s) => s.surface.signed_distance(x),
            Kind::Plane { p, n } => dot(sub(x, *p), *n),
        }
    }

    /// The area of the body's surface inside the box `lo`–`hi`, where
    /// simple (as [`Self::exact_volume`]; a cylinder that spans the box:
    /// its mantle only).
    pub fn exact_area(&self, lo: [f64; 3], hi: [f64; 3]) -> Option<f64> {
        let inside = |p: [f64; 3]| (0..3).all(|i| p[i] >= lo[i] && p[i] <= hi[i]);
        match &self.kind {
            Kind::Sphere { c, r2 } => {
                let r = r2.sqrt();
                (inside(sub(*c, [r; 3])) && inside([c[0] + r, c[1] + r, c[2] + r]))
                    .then(|| 4.0 * std::f64::consts::PI * r2)
            }
            Kind::Box { lo: a, hi: b } => (inside(*a) && inside(*b)).then(|| {
                let e = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                2.0 * (e[0] * e[1] + e[1] * e[2] + e[2] * e[0])
            }),
            Kind::Cylinder { c, a, r2, half } => {
                let i = (0..3).find(|&i| a[i].abs() == 1.0)?;
                let r = r2.sqrt();
                for j in (0..3).filter(|&j| j != i) {
                    if c[j] - r < lo[j] || c[j] + r > hi[j] {
                        return None;
                    }
                }
                if c[i] - half > lo[i] || c[i] + half < hi[i] {
                    return None;
                }
                Some(2.0 * std::f64::consts::PI * r * (hi[i] - lo[i]))
            }
            Kind::Field(_) | Kind::Plane { .. } => None,
            Kind::Stl(s) => {
                let (a, b) = s.surface.bounding_box();
                (s.clean && inside(a) && inside(b)).then_some(s.area)
            }
        }
    }

    /// True if `x` lies inside the body (the surface counts as inside; a
    /// `field` body where χ ≥ ½).
    #[inline]
    pub fn inside(&self, x: [f64; 3]) -> bool {
        if self.invert {
            // the surface counts as inside both ways
            return self.shape_distance(x) >= 0.0;
        }
        match &self.kind {
            Kind::Cylinder { c, a, r2, half } => {
                let d = sub(x, *c);
                let t = dot(d, *a);
                if t.abs() > *half {
                    return false;
                }
                let p = [d[0] - t * a[0], d[1] - t * a[1], d[2] - t * a[2]];
                dot(p, p) <= *r2
            }
            Kind::Sphere { c, r2 } => {
                let d = sub(x, *c);
                dot(d, d) <= *r2
            }
            Kind::Box { lo, hi } => (0..3).all(|i| x[i] >= lo[i] && x[i] <= hi[i]),
            Kind::Field(v) => v.at(x, 0.0) >= 0.5,
            Kind::Stl(s) => s.surface.inside(x),
            Kind::Plane { p, n } => dot(sub(x, *p), *n) <= 0.0,
        }
    }

    /// The exact volume of the body inside the box `lo`–`hi`, where it is
    /// simple to give: a sphere or box inside it, a cylinder along a
    /// coordinate axis that spans it or lies inside it. None otherwise.
    pub fn exact_volume(&self, lo: [f64; 3], hi: [f64; 3]) -> Option<f64> {
        let v = self.shape_volume(lo, hi)?;
        Some(if self.invert {
            (0..3).map(|i| hi[i] - lo[i]).product::<f64>() - v
        } else {
            v
        })
    }

    /// The volume of the shape inside the box, before `invert`.
    fn shape_volume(&self, lo: [f64; 3], hi: [f64; 3]) -> Option<f64> {
        let inside = |p: [f64; 3]| (0..3).all(|i| p[i] >= lo[i] && p[i] <= hi[i]);
        match &self.kind {
            Kind::Sphere { c, r2 } => {
                let r = r2.sqrt();
                (inside(sub(*c, [r; 3])) && inside([c[0] + r, c[1] + r, c[2] + r]))
                    .then(|| 4.0 / 3.0 * std::f64::consts::PI * r * r * r)
            }
            Kind::Box { lo: a, hi: b } => {
                let v: f64 = (0..3)
                    .map(|i| (b[i].min(hi[i]) - a[i].max(lo[i])).max(0.0))
                    .product();
                Some(v)
            }
            Kind::Cylinder { c, a, r2, half } => {
                let i = (0..3).find(|&i| a[i].abs() == 1.0)?;
                let r = r2.sqrt();
                for j in (0..3).filter(|&j| j != i) {
                    if c[j] - r < lo[j] || c[j] + r > hi[j] {
                        return None;
                    }
                }
                let l = (c[i] + half).min(hi[i]) - (c[i] - half).max(lo[i]);
                Some(std::f64::consts::PI * r2 * l.max(0.0))
            }
            Kind::Field(_) | Kind::Plane { .. } => None,
            Kind::Stl(s) => {
                let (a, b) = s.surface.bounding_box();
                (s.clean && inside(a) && inside(b)).then_some(s.volume)
            }
        }
    }
}

/// The checked bodies of a case (its constants for the expressions of
/// `field` bodies).
pub fn bodies_of(case: &crate::case::FlowCase) -> Result<Vec<Body>> {
    bodies_with(&case.immersed, &case.env()?)
}

/// The checked bodies, expressions without case constants.
pub fn bodies(s: &ImmersedSection) -> Result<Vec<Body>> {
    bodies_with(s, &crate::expr::Env::default())
}

/// The checked bodies; `env` holds the constants of the expressions.
pub fn bodies_with(s: &ImmersedSection, env: &crate::expr::Env) -> Result<Vec<Body>> {
    let mut out: Vec<Body> = Vec::with_capacity(s.body.len());
    if s.body.is_empty() && s.topology.is_none() {
        return Err(bad(
            "immersed.body".into(),
            "[immersed] without a [[immersed.body]]".into(),
        ));
    }
    match s.permeability {
        Some(k) if k > 0.0 && k.is_finite() => {}
        None if s.is_cut() => {}
        _ => {
            return Err(bad(
                "immersed.permeability".into(),
                "[immersed] permeability must be given and > 0 (m²); the penalisation error \
                 depends on it, so there is no default (try 1e-4 h² for cell width h)"
                    .into(),
            ));
        }
    }
    if s.is_cut() {
        let f = s.cut.clone().unwrap_or_default().min_volume_fraction();
        if !(0.0..=0.9).contains(&f) {
            return Err(diag::err(
                Diagnostic::new(
                    diag::IMM_CUT,
                    format!("[immersed.cut] min_volume_fraction = {f}: give 0 to 0.9"),
                )
                .key("immersed.cut.min_volume_fraction"),
            ));
        }
    } else if s.cut.is_some() {
        return Err(bad(
            "immersed.cut".into(),
            "[immersed.cut] needs [immersed] method = \"cut\"".into(),
        ));
    }
    if let Some(n) = s.samples
        && !(1..=16).contains(&n)
    {
        return Err(bad(
            "immersed.samples".into(),
            format!("[immersed] samples = {n}: give 1 to 16"),
        ));
    }
    for (i, b) in s.body.iter().enumerate() {
        let key = format!("immersed.body[{i}]");
        let mut findings = Vec::new();
        if b.shape != Shape::Stl && (b.file.is_some() || b.scale.is_some()) {
            return Err(bad(
                format!("{key}.file"),
                format!(
                    "immersed body {}: `file` and `scale` belong to shape = \"stl\"",
                    b.name
                ),
            ));
        }
        if out.iter().any(|o| o.name == b.name) {
            return Err(bad(key, format!("two immersed bodies named {}", b.name)));
        }
        let need = |v: Option<[f64; 3]>, k: &str| {
            v.ok_or_else(|| {
                bad(
                    format!("{key}.{k}"),
                    format!("immersed body {}: a {:?} needs `{k}`", b.name, b.shape),
                )
            })
        };
        let radius = || match b.radius {
            Some(r) if r > 0.0 && r.is_finite() => Ok(r),
            _ => Err(bad(
                format!("{key}.radius"),
                format!("immersed body {}: `radius` must be given and > 0", b.name),
            )),
        };
        let kind = match b.shape {
            Shape::Cylinder => {
                let c = need(b.centre, "centre")?;
                let r = radius()?;
                let a = b.axis.unwrap_or([0.0, 0.0, 1.0]);
                let n = dot(a, a).sqrt();
                if n.is_nan() || n == 0.0 {
                    return Err(bad(
                        format!("{key}.axis"),
                        format!("immersed body {}: the axis has zero length", b.name),
                    ));
                }
                let half = match b.length {
                    None => f64::INFINITY,
                    Some(l) if l > 0.0 && l.is_finite() => 0.5 * l,
                    Some(_) => {
                        return Err(bad(
                            format!("{key}.length"),
                            format!("immersed body {}: `length` must be > 0", b.name),
                        ));
                    }
                };
                Kind::Cylinder {
                    c,
                    a: [a[0] / n, a[1] / n, a[2] / n],
                    r2: r * r,
                    half,
                }
            }
            Shape::Sphere => {
                let c = need(b.centre, "centre")?;
                let r = radius()?;
                Kind::Sphere { c, r2: r * r }
            }
            Shape::Box => {
                let lo = need(b.lower, "lower")?;
                let hi = need(b.upper, "upper")?;
                if (0..3).any(|i| hi[i] <= lo[i]) {
                    return Err(bad(
                        format!("{key}.upper"),
                        format!("immersed body {}: `upper` must exceed `lower`", b.name),
                    ));
                }
                Kind::Box { lo, hi }
            }
            Shape::Field => {
                let f = b.fraction.as_ref().ok_or_else(|| {
                    bad(
                        format!("{key}.fraction"),
                        format!(
                            "immersed body {}: a field body needs `fraction`, the solid \
                             fraction as an expression of x, y, z",
                            b.name
                        ),
                    )
                })?;
                Kind::Field(crate::expr::value(
                    f,
                    env,
                    crate::expr::Vars::SPACE,
                    &format!("{key}.fraction"),
                )?)
            }
            Shape::Stl => {
                let (st, f) = stl_body(b, &key)?;
                findings = f;
                Kind::Stl(st)
            }
            Shape::Plane => {
                let p = need(b.centre, "centre")?;
                let n = need(b.normal, "normal")?;
                let l = dot(n, n).sqrt();
                if !(l > 0.0 && l.is_finite()) {
                    return Err(bad(
                        format!("{key}.normal"),
                        format!("immersed body {}: the normal has zero length", b.name),
                    ));
                }
                Kind::Plane {
                    p,
                    n: [n[0] / l, n[1] / l, n[2] / l],
                }
            }
        };
        let invert = b.invert.unwrap_or(false);
        if invert && !matches!(b.shape, Shape::Cylinder | Shape::Sphere | Shape::Box) {
            return Err(bad(
                format!("{key}.invert"),
                format!(
                    "immersed body {}: `invert` is for a cylinder, sphere or box (turn the \
                     normal of a plane; an STL surface is inverted in its file)",
                    b.name
                ),
            ));
        }
        let rotation = match b.omega {
            None => None,
            Some(w) => {
                if !w.is_finite() {
                    return Err(bad(
                        format!("{key}.omega"),
                        format!("immersed body {}: `omega` must be a number", b.name),
                    ));
                }
                if !s.is_cut() {
                    return Err(bad(
                        format!("{key}.omega"),
                        format!(
                            "immersed body {}: `omega` moves the wall of a cut body \
                             ([immersed] method = \"cut\")",
                            b.name
                        ),
                    ));
                }
                let (axis, origin) = match &kind {
                    Kind::Cylinder { c, a, .. } => (*a, *c),
                    _ => {
                        let a = b.axis.unwrap_or([0.0, 0.0, 1.0]);
                        let n = dot(a, a).sqrt();
                        if !(n > 0.0 && n.is_finite()) {
                            return Err(bad(
                                format!("{key}.axis"),
                                format!("immersed body {}: the axis has zero length", b.name),
                            ));
                        }
                        let o = match (&kind, b.centre) {
                            (_, Some(c)) => c,
                            (Kind::Box { lo, hi }, None) => {
                                std::array::from_fn(|i| 0.5 * (lo[i] + hi[i]))
                            }
                            _ => need(None, "centre")?,
                        };
                        ([a[0] / n, a[1] / n, a[2] / n], o)
                    }
                };
                Some(Rotation {
                    omega: w,
                    axis,
                    origin,
                })
            }
        };
        let patches = match (&b.patches, &kind) {
            (None, _) => None,
            (Some(p), Kind::Stl(st)) => {
                if !s.is_cut() {
                    return Err(bad(
                        format!("{key}.patches"),
                        format!(
                            "immersed body {}: surface patches are walls of a cut body \
                             ([immersed] method = \"cut\")",
                            b.name
                        ),
                    ));
                }
                let (tags, solids) = &*st.tags;
                let (names, of_tag): (Vec<String>, Vec<Option<u32>>) = match p {
                    Patches::Solids(PatchSource::Solids) => {
                        if solids.is_empty() {
                            return Err(bad(
                                format!("{key}.patches"),
                                format!(
                                    "immersed body {}: patches = \"solids\" needs an ASCII STL \
                                     file (a binary file tags its triangles: give \
                                     patches = {{ <name> = <attribute> }})",
                                    b.name
                                ),
                            ));
                        }
                        // one patch per distinct solid name, in file order
                        let mut names: Vec<String> = Vec::new();
                        let mut of_solid = Vec::with_capacity(solids.len());
                        for (k, n) in solids.iter().enumerate() {
                            let n = if n.is_empty() {
                                format!("solid{}", k + 1)
                            } else {
                                n.replace(char::is_whitespace, "_")
                            };
                            let i = match names.iter().position(|m| *m == n) {
                                Some(i) => i,
                                None => {
                                    names.push(n);
                                    names.len() - 1
                                }
                            };
                            of_solid.push(Some(i as u32));
                        }
                        (names, of_solid)
                    }
                    Patches::Tags(t) => {
                        if !solids.is_empty() {
                            return Err(bad(
                                format!("{key}.patches"),
                                format!(
                                    "immersed body {}: attribute tags need a binary STL file \
                                     (an ASCII file names its solids: give patches = \"solids\")",
                                    b.name
                                ),
                            ));
                        }
                        let names: Vec<String> = t.keys().cloned().collect();
                        let top = t.values().copied().max().unwrap_or(0) as usize;
                        let mut of_tag = vec![None; top + 1];
                        for (i, v) in t.values().enumerate() {
                            if of_tag[*v as usize].is_some() {
                                return Err(bad(
                                    format!("{key}.patches"),
                                    format!(
                                        "immersed body {}: two patches with the attribute {v}",
                                        b.name
                                    ),
                                ));
                            }
                            of_tag[*v as usize] = Some(i as u32);
                        }
                        (names, of_tag)
                    }
                };
                for n in &names {
                    if n.is_empty()
                        || !n
                            .chars()
                            .all(|c| c.is_ascii_alphanumeric() || "_-".contains(c))
                    {
                        return Err(bad(
                            format!("{key}.patches"),
                            format!(
                                "immersed body {}: the patch name `{n}` must be letters, digits, \
                                 `_` or `-`",
                                b.name
                            ),
                        ));
                    }
                }
                let per_tri: Vec<Option<u32>> = tags
                    .iter()
                    .map(|&t| of_tag.get(t as usize).copied().flatten())
                    .collect();
                Some((names, per_tri))
            }
            (Some(_), _) => {
                return Err(bad(
                    format!("{key}.patches"),
                    format!(
                        "immersed body {}: `patches` belong to shape = \"stl\"",
                        b.name
                    ),
                ));
            }
        };
        out.push(Body {
            name: b.name.clone(),
            kind,
            invert,
            rotation,
            patches,
            findings,
        });
    }
    Ok(out)
}

/// Reads and checks the surface of an `stl` body (IMM001 is an error, the
/// other findings of the check are reported).
fn stl_body(
    b: &BodySection,
    key: &str,
) -> Result<(StlBody, Vec<foehncraft_mesh::surface::check::Finding>)> {
    let file = b.file.as_ref().ok_or_else(|| {
        bad(
            format!("{key}.file"),
            format!(
                "immersed body {}: an stl body needs `file`, the STL file",
                b.name
            ),
        )
    })?;
    let scale = b.scale.unwrap_or(1.0);
    if !(scale > 0.0 && scale.is_finite()) {
        return Err(bad(
            format!("{key}.scale"),
            format!("immersed body {}: `scale` must be > 0", b.name),
        ));
    }
    // the case is checked and read several times per run: the surface,
    // its hierarchy and its check are kept for the same file and scale
    let meta = std::fs::metadata(file).ok();
    let id: StlKey = (
        file.clone(),
        scale.to_bits(),
        meta.as_ref().map_or(0, |m| m.len()),
        meta.and_then(|m| m.modified().ok()),
    );
    let cache = STL_CACHE.get_or_init(Default::default);
    if let Ok(c) = cache.lock()
        && let Some((_, hit)) = c.iter().find(|(k, _)| *k == id)
    {
        return Ok(hit.clone());
    }
    let (surface, tags, solids) = foehncraft_mesh::surface::Surface::read_stl_tagged(file, scale)
        .map_err(|e| {
        diag::err(
            Diagnostic::new(
                diag::CASE_FILE,
                format!("immersed body {}: {}: {e}", b.name, file.display()),
            )
            .key(format!("{key}.file")),
        )
    })?;
    let report = foehncraft_mesh::surface::check(&surface);
    if let Some(f) = report
        .findings
        .iter()
        .find(|f| f.code == diag::IMM_SURFACE_EMPTY)
    {
        return Err(diag::err(
            Diagnostic::new(
                diag::IMM_SURFACE_EMPTY,
                format!("immersed body {}: {}", b.name, f.message),
            )
            .key(format!("{key}.file")),
        ));
    }
    let out = (
        StlBody {
            surface: std::sync::Arc::new(surface),
            tags: std::sync::Arc::new((tags, solids)),
            clean: report.clean(),
            volume: report.volume,
            area: report.area,
        },
        report.findings,
    );
    if let Ok(mut c) = cache.lock() {
        if c.len() >= 8 {
            c.remove(0);
        }
        c.push((id, out.clone()));
    }
    Ok(out)
}

/// File, scale, length and modification time of a read STL body.
type StlKey = (std::path::PathBuf, u64, u64, Option<std::time::SystemTime>);

/// The STL bodies read in this process.
#[allow(clippy::type_complexity)]
static STL_CACHE: std::sync::OnceLock<
    std::sync::Mutex<
        Vec<(
            StlKey,
            (StlBody, Vec<foehncraft_mesh::surface::check::Finding>),
        )>,
    >,
> = std::sync::OnceLock::new();

/// The settings that do not go with a penalised body (IMM102).
pub fn check_case(case: &crate::case::FlowCase) -> Result<()> {
    if case.immersed.is_empty() {
        return Ok(());
    }
    bodies_of(case)?;
    if case.immersed.is_cut() {
        let what = if case.immersed.topology.is_some() {
            Some("a topology search (it runs on the penalised density field)")
        } else if !case.mrf.is_empty() || !case.sliding.is_empty() {
            Some("rotating zones")
        } else {
            None
        };
        return match what {
            Some(w) => Err(diag::err(
                Diagnostic::new(
                    diag::IMM_CUT_UNSUPPORTED,
                    format!("[immersed] method = \"cut\" with {w} is not available yet"),
                )
                .key("immersed.method"),
            )),
            None => Ok(()),
        };
    }
    if let Some(t) = &case.immersed.topology {
        t.check()?;
    }
    let what = if !case.turbulence.is_laminar() {
        Some(
            "a turbulence model (penalised walls are laminar and low-Re only; wall functions do \
             not fit a smeared wall)",
        )
    } else if case.energy.is_some() {
        Some("the energy equation (no penalised conjugate heat transfer)")
    } else if !case.scalars.is_empty() {
        Some("passive scalars")
    } else if case.gravity.is_some() {
        Some("gravity")
    } else if !case.regions.is_empty() {
        Some("regions")
    } else if !case.mrf.is_empty() || !case.sliding.is_empty() {
        Some("rotating zones")
    } else if !case.spray.is_empty() {
        Some("sprays")
    } else if case.immersed.topology.is_some() && case.time.is_some() {
        Some("a transient run (the topology search is steady)")
    } else {
        None
    };
    match what {
        Some(w) => Err(diag::err(
            Diagnostic::new(
                diag::IMM_UNSUPPORTED,
                format!(
                    "[immersed] bodies with {w} are not available (L13 stage 1 is laminar, \
                     isothermal flow; sharp walls come with the cut cells of stage 2)"
                ),
            )
            .key("immersed"),
        )),
        None => Ok(()),
    }
}

/// The midpoint rule of `n³` points in collapsed coordinates of the unit
/// tetrahedron: barycentric-free maps (s, t, w) ∈ [0, 1]³ to
/// x = v0 + s (e1 + t (e2 + w e3)) with e1 = v1 − v0, e2 = v2 − v1,
/// e3 = v3 − v2, Jacobian ∝ s² t; returns the points as (s, st, stw) and
/// the normalised weights.
fn rule(n: usize) -> Vec<([f64; 3], f64)> {
    let mut out = Vec::with_capacity(n * n * n);
    let h = 1.0 / n as f64;
    let mut sum = 0.0;
    for i in 0..n {
        let s = (i as f64 + 0.5) * h;
        for j in 0..n {
            let t = (j as f64 + 0.5) * h;
            for k in 0..n {
                let w = (k as f64 + 0.5) * h;
                let wt = s * s * t;
                sum += wt;
                out.push(([s, s * t, s * t * w], wt));
            }
        }
    }
    for p in &mut out {
        p.1 /= sum;
    }
    out
}

/// The solid fraction χ of every local cell, one array per body: sampled
/// for the primitive and field bodies, exact for STL bodies (the inside
/// volume of the cut moments of stage 0 over the cell volume).
pub fn fractions(bodies: &[Body], samples: usize, mesh: &Mesh, geo: &Geometry) -> Vec<Vec<f64>> {
    let nc = mesh.n_cells;
    let q = rule(samples);
    let mut out = vec![vec![0.0; nc]; bodies.len()];
    let sampled: Vec<usize> = (0..bodies.len())
        .filter(|&b| bodies[b].surface().is_none())
        .collect();
    for (b, body) in bodies.iter().enumerate() {
        if let Some(s) = body.surface() {
            let m = s.cut_moments(mesh, geo);
            for (o, (&vi, &v)) in out[b]
                .iter_mut()
                .zip(m.cell_inside_volume.iter().zip(&geo.cell_volume))
            {
                *o = if v > 0.0 {
                    (vi / v).clamp(0.0, 1.0)
                } else {
                    0.0
                };
            }
        }
    }
    if sampled.is_empty() {
        return out;
    }
    let mut inside = vec![0.0; bodies.len()];
    for c in 0..nc {
        let xc = geo.cell_centroid[c];
        let mut vol = 0.0;
        inside.iter_mut().for_each(|v| *v = 0.0);
        for (f, _) in mesh.cell_face_iter(c) {
            let nodes = mesh.face(f);
            let xf = geo.face_centroid[f];
            for e in 0..nodes.len() {
                let p = mesh.points[nodes[e] as usize];
                let r = mesh.points[nodes[(e + 1) % nodes.len()] as usize];
                // the tetrahedron (centroid, face centroid, p, r)
                let (e1, e2, e3) = (sub(xf, xc), sub(p, xf), sub(r, p));
                let v = dot(sub(xf, xc), cross(sub(p, xc), sub(r, xc))).abs() / 6.0;
                if v == 0.0 {
                    continue;
                }
                vol += v;
                for &b in &sampled {
                    let body = &bodies[b];
                    let mut frac = 0.0;
                    for (z, w) in &q {
                        let x = [
                            xc[0] + z[0] * e1[0] + z[1] * e2[0] + z[2] * e3[0],
                            xc[1] + z[0] * e1[1] + z[1] * e2[1] + z[2] * e3[1],
                            xc[2] + z[0] * e1[2] + z[1] * e2[2] + z[2] * e3[2],
                        ];
                        frac += w * body.value(x);
                    }
                    inside[b] += v * frac;
                }
            }
        }
        for &b in &sampled {
            out[b][c] = if vol > 0.0 {
                (inside[b] / vol).clamp(0.0, 1.0)
            } else {
                0.0
            };
        }
    }
    out
}

/// The solid fraction of all bodies together (capped at 1).
pub fn total_fraction(per_body: &[Vec<f64>], nc: usize) -> Vec<f64> {
    (0..nc)
        .map(|c| per_body.iter().map(|b| b[c]).sum::<f64>().min(1.0))
        .collect()
}

/// The penalisation attached to a flow: the device buffer of the sink
/// s_P = μ χ_P V_P / K (kg/s) per cell.
#[derive(Debug, Clone)]
pub struct Attached {
    /// The sink on the device.
    pub sink: crate::device::Buf,
}

/// One body in the report: its solid volume and the force of the fluid on
/// it, Σ s_P u_P (N).
#[derive(Debug, Clone, Serialize)]
pub struct BodyForce {
    /// Name.
    pub name: String,
    /// Σ χ_P V_P (m³).
    pub solid_volume: f64,
    /// The exact volume inside the mesh's bounding box, where simple.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub exact_volume: Option<f64>,
    /// Cells with χ > 0.
    pub cells: usize,
    /// Force of the fluid on the body (N).
    pub force: [f64; 3],
}

/// The solid fraction χ of every local cell of the primitive and field
/// bodies (zero everywhere without bodies).
pub fn body_fraction(bs: &[Body], samples: usize, mesh: &Mesh, geo: &Geometry) -> Vec<f64> {
    if bs.is_empty() {
        return vec![0.0; mesh.n_cells];
    }
    total_fraction(&fractions(bs, samples, mesh, geo), mesh.n_cells)
}

/// Sets the penalisation of `case` on the flow (nothing without
/// `[immersed]`); returns the diagnostics (IMM103, IMM104) and the topology
/// search of `[immersed.topology]`, if any.
pub fn attach(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
) -> Result<(Vec<Diagnostic>, Option<crate::topology::Search>)> {
    if case.immersed.is_empty() {
        return Ok((Vec::new(), None));
    }
    check_case(case)?;
    if case.immersed.is_cut() {
        // the bodies are walls of the cut mesh (`cut_mesh`)
        return Ok((Vec::new(), None));
    }
    let bs = bodies_of(case)?;
    let k = case.immersed.permeability.unwrap_or(f64::INFINITY);
    let mu = case.mu();
    let nc = mesh.n_cells;
    let per = fractions(&bs, case.immersed.samples(), mesh, geo);
    let chi = if bs.is_empty() {
        vec![0.0; nc]
    } else {
        total_fraction(&per, nc)
    };
    let mut notes = surface_notes(&bs);
    let (lo, hi) = mesh.bounding_box();
    for (b, x) in bs.iter().zip(&per) {
        let vol: f64 = x.iter().zip(&geo.cell_volume).map(|(a, v)| a * v).sum();
        let cells = x.iter().filter(|&&a| a > 0.0).count();
        let half = x.iter().any(|&a| a > 0.5);
        let exact = b.exact_volume(lo, hi);
        let mut line = format!(
            "immersed body {}: solid volume {vol:.6e} m³ in {cells} cells, μ/K = {:.3e} kg/(m³ s)",
            b.name,
            mu / k
        );
        if let Some(e) = exact {
            line.push_str(&format!(
                ", exact {e:.6e} m³ ({:+.2} %)",
                100.0 * (vol - e) / e.max(f64::MIN_POSITIVE)
            ));
        }
        notes.push(Diagnostic::new(diag::IMM_BODIES, line));
        let off = exact.is_some_and(|e| e > 0.0 && ((vol - e) / e).abs() > 0.05);
        if !half || off {
            notes.push(Diagnostic::new(
                diag::IMM_RESOLUTION,
                format!(
                    "immersed body {}: the mesh does not resolve it ({}); refine the mesh \
                     around it",
                    b.name,
                    if half {
                        "solid volume more than 5 % off the exact volume"
                    } else {
                        "no cell lies more than half inside it"
                    }
                ),
            ));
        }
    }
    let search = crate::topology::Search::new(case, mesh, geo, chi.clone())?;
    match &search {
        Some(t) => {
            notes.push(Diagnostic::new(
                diag::IMM_BODIES,
                format!(
                    "topology search: {} design cells, fluid fraction {}, μ/K = {:.3e} kg/(m³ s)",
                    t.cells.len(),
                    case.immersed
                        .topology
                        .as_ref()
                        .map_or(0.0, |s| s.fluid_fraction),
                    mu / k
                ),
            ));
            fl.set_immersed(&t.sink(geo))?;
            fl.set_immersed_fields(t.solid(), t.design_chi(nc));
        }
        None => {
            let sink: Vec<f64> = chi
                .iter()
                .zip(&geo.cell_volume)
                .map(|(x, v)| mu * x * v / k)
                .collect();
            fl.set_immersed(&sink)?;
            fl.set_immersed_fields(chi, Vec::new());
        }
    }
    Ok((notes, search))
}

/// Puts the current design of a topology search on the flow.
pub fn update_design(
    fl: &mut crate::simplec::Flow,
    t: &crate::topology::Search,
    geo: &Geometry,
) -> Result<()> {
    fl.set_immersed(&t.sink(geo))?;
    fl.set_immersed_fields(t.solid(), t.design_chi(geo.cell_volume.len()));
    Ok(())
}

/// The force of the fluid on every body, Σ s_P u_P, from the host velocity;
/// with a topology search the design is the last entry (`design`), with
/// χ from `phys`.
pub fn forces(
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    u: &foehncraft_fv::disc::Soa3,
    phys: &crate::simplec::Physics,
) -> Result<Vec<BodyForce>> {
    if case.immersed.is_empty() || case.immersed.is_cut() {
        // cut bodies are wall boundaries: their forces are in the boundaries
        return Ok(Vec::new());
    }
    let bs = bodies_of(case)?;
    let k = case.immersed.permeability.unwrap_or(f64::INFINITY);
    let mu = case.mu();
    let per = fractions(&bs, case.immersed.samples(), mesh, geo);
    let (lo, hi) = mesh.bounding_box();
    let one = |name: String, x: &[f64], exact: Option<f64>| {
        let mut f = [0.0; 3];
        let mut vol = 0.0;
        for (c, &a) in x.iter().enumerate() {
            if a == 0.0 {
                continue;
            }
            let v = geo.cell_volume[c];
            vol += a * v;
            let s = mu * a * v / k;
            let uc = u.get(c);
            for i in 0..3 {
                f[i] += s * uc[i];
            }
        }
        BodyForce {
            name,
            solid_volume: vol,
            exact_volume: exact,
            cells: x.iter().filter(|&&a| a > 0.0).count(),
            force: f,
        }
    };
    let mut out: Vec<BodyForce> = bs
        .iter()
        .zip(&per)
        .map(|(b, x)| one(b.name.clone(), x, b.exact_volume(lo, hi)))
        .collect();
    if case.immersed.topology.is_some()
        && let Some(f) = &phys.immersed_fields
        && !f.1.is_empty()
    {
        out.push(one("design".into(), &f.1, None));
    }
    Ok(out)
}

/// The bodies of a cut as one implicit function: φ is the smallest signed
/// distance (narrow band for STL bodies), the crossing point on an edge is
/// the nearest entry into any body the edge ends in (the exact crossing of
/// an STL body, the root of a primitive's own function; or linear), the
/// wall and face points are those of the nearest body (dual contouring on
/// the sharp edges of an STL body, the Newton projection or the bisector
/// root of a primitive).
struct CutBodies<'a> {
    bodies: &'a [Body],
    exact: bool,
    /// Step of the gradient of the projection (m).
    h: f64,
}

impl CutBodies<'_> {
    fn nearest(&self, x: [f64; 3]) -> usize {
        let mut best = (f64::INFINITY, 0);
        for (i, b) in self.bodies.iter().enumerate() {
            let d = match b.surface() {
                Some(s) => s.distance(x),
                None => b.distance(x).abs(),
            };
            if d < best.0 {
                best = (d, i);
            }
        }
        best.1
    }

    fn stl(&self, i: usize) -> Option<foehncraft_mesh::surface::SurfaceBody<'_>> {
        self.bodies[i]
            .surface()
            .map(|s| foehncraft_mesh::surface::SurfaceBody {
                surface: s,
                exact: self.exact,
                features: true,
            })
    }

    /// The unit outward normal of primitive body `b` at `x` (central
    /// differences of its distance).
    fn normal(&self, b: &Body, x: [f64; 3]) -> Option<[f64; 3]> {
        let mut g = [0.0; 3];
        for (k, gk) in g.iter_mut().enumerate() {
            let (mut a, mut c) = (x, x);
            a[k] += self.h;
            c[k] -= self.h;
            *gk = b.distance(a) - b.distance(c);
        }
        let n = norm(g);
        (n > 0.0 && n.is_finite()).then(|| [g[0] / n, g[1] / n, g[2] / n])
    }
}

impl foehncraft_mesh::cut::ImplicitBody for CutBodies<'_> {
    fn value(&self, x: [f64; 3]) -> f64 {
        self.bodies
            .iter()
            .map(|b| b.distance(x))
            .fold(f64::INFINITY, f64::min)
    }
    fn values(&self, points: &[[f64; 3]], reach: &[f64]) -> Vec<f64> {
        use rayon::prelude::*;
        let mut phi = vec![f64::INFINITY; points.len()];
        for (i, b) in self.bodies.iter().enumerate() {
            let v: Vec<f64> = match self.stl(i) {
                Some(s) => s.values(points, reach),
                None => points.par_iter().map(|&x| b.distance(x)).collect(),
            };
            for (p, x) in phi.iter_mut().zip(v) {
                *p = p.min(x);
            }
        }
        phi
    }
    fn edge_point(&self, a: [f64; 3], b: [f64; 3], fa: f64, fb: f64) -> [f64; 3] {
        let linear = || {
            let t = fa / (fa - fb);
            [
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ]
        };
        if !self.exact {
            return linear();
        }
        // the nearest entry along the edge into a body that holds `b`
        let mut best: Option<(f64, [f64; 3])> = None;
        for (i, body) in self.bodies.iter().enumerate() {
            let x = match body.surface() {
                Some(s) => {
                    if !s.inside(b) {
                        continue;
                    }
                    match s.first_crossing(a, b) {
                        Some(c) => c.point,
                        None => continue,
                    }
                }
                None => {
                    let (da, db) = (body.distance(a), body.distance(b));
                    if db > 0.0 || da <= 0.0 {
                        continue;
                    }
                    foehncraft_mesh::cut::edge_root(|x| body.distance(x), a, b, da, db)
                }
            };
            let t = norm(sub(x, a));
            if best.is_none_or(|(bt, _)| t < bt) {
                best = Some((t, x));
            }
            let _ = i;
        }
        match best {
            Some((_, x)) => x,
            None => foehncraft_mesh::cut::edge_root(|x| self.value(x), a, b, fa, fb),
        }
    }
    fn project(&self, x: [f64; 3]) -> Option<[f64; 3]> {
        let i = self.nearest(x);
        match self.stl(i) {
            Some(s) => s.project(x),
            None => {
                let b = &self.bodies[i];
                foehncraft_mesh::cut::project_newton(|y| b.distance(y), x, self.h)
            }
        }
    }
    fn wall_point(&self, ring: &[[f64; 3]], centre: [f64; 3]) -> Option<[f64; 3]> {
        let i = self.nearest(centre);
        match self.stl(i) {
            Some(s) => s.wall_point(ring, centre),
            None => {
                let b = &self.bodies[i];
                // a sharp edge of a box: the dual-contouring point
                let n: Vec<[f64; 3]> = ring
                    .iter()
                    .map(|&p| self.normal(b, p).unwrap_or([0.0; 3]))
                    .collect();
                if let Some(x) = foehncraft_mesh::cut::feature_point(ring, &n, centre) {
                    let r = ring
                        .iter()
                        .map(|&y| norm(sub(y, centre)))
                        .fold(0.0, f64::max);
                    if b.distance(x).abs() <= 1e-7 * r {
                        return Some(x);
                    }
                }
                self.project(centre)
            }
        }
    }
    fn face_point(&self, p: [f64; 3], q: [f64; 3], n: [f64; 3]) -> Option<[f64; 3]> {
        let i = self.nearest(scale3(add3(p, q), 0.5));
        match self.stl(i) {
            Some(s) => s.face_point(p, q, n),
            None => {
                let b = &self.bodies[i];
                foehncraft_mesh::cut::chord_point(|x| b.distance(x), |x| self.normal(b, x), p, q, n)
            }
        }
    }
    fn patch(&self, x: [f64; 3]) -> usize {
        self.nearest(x)
    }
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// The solid region a cut body is kept as (L13 stage 2d, D229): the
/// `[regions.<name>] kind = "solid"` whose `zones` name the body.
pub fn solid_region_of<'a>(case: &'a crate::case::FlowCase, body: &str) -> Option<&'a str> {
    case.regions
        .iter()
        .find(|(_, r)| {
            r.kind == crate::materials::RegionKind::Solid && r.zones.iter().any(|z| z == body)
        })
        .map(|(n, _)| n.as_str())
}

/// The cut of `[immersed] method = "cut"` (L13 stage 2, D227): the mesh
/// with the bodies cut out, their surfaces as wall boundaries named after
/// the bodies, and the notes for the log (IMM203, IMM204). A body that is
/// a solid region keeps its inside as cells of a zone named by the body,
/// joined to the fluid cells by interior faces (stage 2d, D229); thin ones
/// are refused (IMM205).
pub fn cut_mesh(case: &crate::case::FlowCase, mesh: &Mesh) -> Result<(Mesh, Vec<Diagnostic>)> {
    let bs = bodies_of(case)?;
    let solid: Vec<bool> = bs
        .iter()
        .map(|b| solid_region_of(case, &b.name).is_some())
        .collect();
    for (b, &s) in bs.iter().zip(&solid) {
        if s && case.energy.is_none() {
            return Err(diag::err(
                Diagnostic::new(
                    diag::IMM_CUT_UNSUPPORTED,
                    format!(
                        "the cut body `{}` is a solid region, which needs the energy equation \
                         ([energy])",
                        b.name
                    ),
                )
                .key("regions"),
            ));
        }
    }
    let sec = case.immersed.cut.clone().unwrap_or_default();
    let (lo, hi) = mesh.bounding_box();
    let diag_len = norm(sub(hi, lo));
    let body = CutBodies {
        bodies: &bs,
        exact: sec.edge_root.unwrap_or_default() == EdgeRoot::Exact,
        h: 1e-7 * diag_len.max(f64::MIN_POSITIVE),
    };
    let names: Vec<String> = bs.iter().map(|b| b.name.clone()).collect();
    let opt = foehncraft_mesh::cut::CutOptions {
        min_fraction: sec.min_volume_fraction(),
        wall_point: sec.wall_point.unwrap_or(true),
        face_point: sec.face_point(),
    };
    let c = foehncraft_mesh::cut::cut_both(mesh, &body, &names, &solid, &opt)
        .map_err(|e| diag::err(Diagnostic::new(diag::IMM_CUT, format!("cut cells: {e}"))))?;
    let r = &c.report;
    let mut notes = surface_notes(&bs);
    notes.push(Diagnostic::new(
        diag::IMM_CUT_REPORT,
        format!(
            "cut cells: {} of {} cells cut, {} removed, {} extra from split cells, {} merged \
             (min_volume_fraction {}), smallest kept fraction {:.3}, closure {:.1e}; {} cells",
            r.cells_cut,
            r.cells_in,
            r.cells_removed,
            r.cells_split,
            r.cells_merged,
            opt.min_fraction,
            r.min_fraction,
            r.max_closure,
            r.cells_out
        ),
    ));
    for (k, (b, a)) in bs.iter().zip(&r.wall_area).enumerate() {
        let mut line = if solid[k] {
            format!(
                "cut body {} (solid region {}): interface area {a:.6e} m²",
                b.name,
                solid_region_of(case, &b.name).unwrap_or_default()
            )
        } else {
            format!("cut body {}: wall area {a:.6e} m²", b.name)
        };
        if let Some(e) = b.exact_area(lo, hi) {
            line.push_str(&format!(
                ", exact {e:.6e} m² ({:+.3} %)",
                100.0 * (a - e) / e
            ));
        }
        if solid[k] {
            line.push_str(&format!(", solid volume {:.6e} m³", r.solid_volume[k]));
            if let Some(e) = b.exact_volume(lo, hi) {
                let v = r.solid_volume[k];
                line.push_str(&format!(
                    ", exact {e:.6e} m³ ({:+.3} %)",
                    100.0 * (v - e) / e
                ));
            }
        }
        notes.push(Diagnostic::new(diag::IMM_CUT_REPORT, line));
    }
    if solid.iter().any(|&s| s) {
        notes.push(Diagnostic::new(
            diag::IMM_CUT_REPORT,
            format!(
                "cut cells both ways: {} solid cells, {} small ones merged, smallest kept \
                 fraction {:.3}",
                r.solid_cells, r.solid_merged, r.solid_min_fraction
            ),
        ));
        thin_solids(&bs, &solid, &c, mesh)?;
    }
    if r.simplified > 0 {
        notes.push(Diagnostic::new(
            diag::IMM_CUT_REPORT,
            format!(
                "cut cells: {} cells cut with fewer wall and face points so that the \
                 discretisation takes their faces ({} rounds)",
                r.simplified, r.rounds
            ),
        ));
    }
    if r.shared_pairs > 0 {
        notes.push(Diagnostic::new(
            diag::IMM_CUT_REPORT,
            format!(
                "cut cells: {} merged cells share several faces with a neighbour (faces \
                 that do not lie in one plane are kept apart, not joined)",
                r.shared_pairs
            ),
        ));
    }
    if r.merge_failed > 0 || r.bad_faces > 0 {
        notes.push(Diagnostic::new(
            diag::IMM_CUT_QUALITY,
            format!(
                "cut cells: {} small cells could not be merged, {} faces with the normal \
                 against the centroid distance; refine the mesh near the body",
                r.merge_failed, r.bad_faces
            ),
        ));
    }
    let mut out = c.mesh;
    // a boundary of the input mesh that lies inside the bodies is gone: it
    // needs no condition (the box around an inverted cylinder, stage 2e);
    // one the case gives a condition stays, empty
    let walls: Vec<&str> = bs.iter().map(|b| b.name.as_str()).collect();
    let keep = |p: &foehncraft_mesh::mesh::Boundary| {
        p.len > 0 || walls.contains(&p.name.as_str()) || case.boundary.contains_key(&p.name)
    };
    let gone: Vec<String> = out
        .boundaries
        .iter()
        .filter(|p| !keep(p))
        .map(|p| p.name.clone())
        .collect();
    if !gone.is_empty() {
        out.boundaries.retain(keep);
        notes.push(Diagnostic::new(
            diag::IMM_CUT_REPORT,
            format!(
                "cut cells: the boundaries {} lie inside the bodies and are left out",
                gone.join(", ")
            ),
        ));
    }
    split_patches(&bs, &solid, &mut out, &mut notes)?;
    Ok((out, notes))
}

/// The walls that turn with their cut body (`omega`, L13 stage 2e): the
/// boundary names of the body and its patches with the rotation.
pub fn wall_rotations(
    case: &crate::case::FlowCase,
) -> Result<std::collections::BTreeMap<String, Rotation>> {
    let mut out = std::collections::BTreeMap::new();
    if !case.immersed.is_cut() || case.immersed.body.iter().all(|b| b.omega.is_none()) {
        return Ok(out);
    }
    for b in bodies_of(case)? {
        let Some(r) = b.rotation else { continue };
        if solid_region_of(case, &b.name).is_some() {
            return Err(diag::err(
                Diagnostic::new(
                    diag::IMM_CUT_UNSUPPORTED,
                    format!(
                        "the cut body `{}` is a solid region: a turning solid is rigid motion \
                         (stage 5), `omega` moves walls only",
                        b.name
                    ),
                )
                .key("immersed.body.omega"),
            ));
        }
        for n in wall_names(&b) {
            out.insert(n, r);
        }
    }
    Ok(out)
}

/// Gives the walls of turning cut bodies their velocity ω × (x − o) in the
/// compiled boundary values (L13 stage 2e); a wall of such a body with a
/// `velocity` of its own is refused.
pub fn moving_walls(case: &crate::case::FlowCase, c: &mut crate::bc::Compiled) -> Result<()> {
    for (name, r) in wall_rotations(case)? {
        let Some(crate::bc::BoundaryCondition::Wall { velocity, .. }) = case.boundary.get(&name)
        else {
            continue;
        };
        if velocity.iter().any(|v| v.number() != Some(0.0)) {
            return Err(bad(
                format!("boundary.{name}.velocity"),
                format!(
                    "the wall `{name}` turns with its body (`omega`): give its speed by \
                     `omega` or by `velocity`, not both"
                ),
            ));
        }
        let [ax, ay, az] = r.axis;
        let [ox, oy, oz] = r.origin;
        let w = r.omega;
        let e = |a: f64, u: &str, ou: f64, b: f64, v: &str, ov: f64| {
            crate::expr::Scalar::Expr(format!(
                "({w:?}) * (({a:?}) * ({u} - ({ou:?})) - ({b:?}) * ({v} - ({ov:?})))"
            ))
        };
        let v = [
            e(ay, "z", oz, az, "y", oy),
            e(az, "x", ox, ax, "z", oz),
            e(ax, "y", oy, ay, "x", ox),
        ];
        let vals = crate::expr::vector_value(
            &v,
            &case.env()?,
            crate::expr::Vars::SPACE,
            &format!("immersed.body.omega ({name})"),
        )?;
        c.entry(name).or_default().velocity = Some(vals);
    }
    Ok(())
}

/// The boundary of every name of a cut body's walls: the body itself and
/// its surface patches `<body>-<patch>` (L13 stage 2e).
pub fn wall_names(b: &Body) -> Vec<String> {
    let mut v = vec![b.name.clone()];
    if let Some((names, _)) = &b.patches {
        v.extend(names.iter().map(|n| format!("{}-{n}", b.name)));
    }
    v
}

/// Splits the wall boundary of every STL body with `patches` into one
/// boundary per patch (L13 stage 2e): a face takes the patch of the triangle
/// closest to its centroid.
fn split_patches(
    bs: &[Body],
    solid: &[bool],
    m: &mut Mesh,
    notes: &mut Vec<Diagnostic>,
) -> Result<()> {
    if !bs.iter().any(|b| b.patches.is_some()) {
        return Ok(());
    }
    let geo = Geometry::compute(m);
    for (b, &s) in bs.iter().zip(solid) {
        let (Some((names, per_tri)), Some(surf)) = (&b.patches, b.surface()) else {
            continue;
        };
        if s {
            return Err(diag::err(
                Diagnostic::new(
                    diag::IMM_CUT_UNSUPPORTED,
                    format!(
                        "the cut body `{}` is a solid region: its surface is the region \
                         interface, which has no patches",
                        b.name
                    ),
                )
                .key("immersed.body.patches"),
            ));
        }
        let Some(range) = m.boundary(&b.name).map(|p| p.range()) else {
            continue;
        };
        // part 0 is the body itself, part k + 1 patch k
        let part: Vec<usize> = range
            .clone()
            .map(|f| {
                surf.bvh()
                    .closest(geo.face_centroid[f])
                    .and_then(|c| per_tri.get(c.global_id as usize).copied().flatten())
                    .map_or(0, |k| k as usize + 1)
            })
            .collect();
        let all = wall_names(b);
        let mut line = format!("cut body {}: surface patches", b.name);
        for (k, n) in all.iter().enumerate() {
            let a: f64 = range
                .clone()
                .zip(&part)
                .filter(|(_, p)| **p == k)
                .map(|(f, _)| norm(geo.face_area[f]))
                .sum();
            let nf = part.iter().filter(|p| **p == k).count();
            if k > 0 || nf > 0 {
                line.push_str(&format!(", {n} {nf} faces {a:.6e} m²"));
            }
        }
        m.split_boundary(&b.name, &part, &all)
            .map_err(|e| diag::err(Diagnostic::new(diag::IMM_CUT, format!("cut cells: {e}"))))?;
        let _ = names;
        notes.push(Diagnostic::new(diag::IMM_CUT_REPORT, line));
    }
    Ok(())
}

/// The refusal of thin solids (L13 stage 2d, owner decision Q4): a body kept
/// as solid cells must be about two cells thick. Its depth is the largest
/// distance of a solid cell centroid to the surface, the cell width the
/// median over its cut cells of the uncut volume over the interface area
/// in the cell (the width across the surface, in 2D and 3D alike); the body
/// is refused when the depth is below one cell width (a slab thinner than
/// two cells, a pin thinner than two cells across) or when small solid
/// cells found no solid neighbour to merge into.
fn thin_solids(
    bs: &[Body],
    solid: &[bool],
    c: &foehncraft_mesh::cut::CutMesh,
    input: &Mesh,
) -> Result<()> {
    let m = &c.mesh;
    let geo = Geometry::compute(m);
    let geo0 = Geometry::compute(input);
    let thin = |name: &str, why: String| {
        Err(diag::err(
            Diagnostic::new(
                diag::IMM_THIN_SOLID,
                format!(
                    "the cut body `{name}` is too thin for a solid region on this mesh: {why}; \
                     refine the mesh in the body, or model it as a thin wall on a region \
                     interface ([interfaces.<name>] thin_wall, D189)"
                ),
            )
            .key("regions"),
        ))
    };
    if c.report.solid_merge_failed > 0 {
        let name = bs
            .iter()
            .zip(solid)
            .find(|(_, s)| **s)
            .map_or("", |(b, _)| b.name.as_str());
        return thin(
            name,
            format!(
                "{} small solid cells found no solid neighbour to merge into",
                c.report.solid_merge_failed
            ),
        );
    }
    let nf = c.n_fluid;
    for (b, _) in bs.iter().zip(solid).filter(|(_, s)| **s) {
        // a `field` body has a level, not a distance
        if matches!(b.kind, Kind::Field(_)) {
            continue;
        }
        let Some(z) = m.zones.iter().position(|z| *z == b.name) else {
            continue;
        };
        let cells: Vec<usize> = (nf..m.n_cells)
            .filter(|&k| m.cell_zone[k] as usize == z)
            .collect();
        if cells.is_empty() {
            return thin(&b.name, "no node of the mesh lies inside it".into());
        }
        let depth = cells
            .iter()
            .map(|&k| -b.distance(geo.cell_centroid[k]))
            .fold(0.0f64, f64::max);
        let mut widths: Vec<f64> = Vec::new();
        for &k in &cells {
            let a: f64 = m
                .cell_face_iter(k)
                .filter(|&(f, _)| {
                    f < m.n_interior_faces() && {
                        let o = if m.owner[f] as usize == k {
                            m.neighbour[f]
                        } else {
                            m.owner[f]
                        };
                        (o as usize) < nf
                    }
                })
                .map(|(f, _)| norm(geo.face_area[f]))
                .sum();
            if a > 0.0 {
                let v: f64 = c.covers[k]
                    .iter()
                    .map(|&p| geo0.cell_volume[p as usize])
                    .sum();
                widths.push(v / a);
            }
        }
        if widths.is_empty() {
            continue;
        }
        widths.sort_by(f64::total_cmp);
        let h = widths[widths.len() / 2];
        if depth < h {
            return thin(
                &b.name,
                format!(
                    "its deepest solid cell lies {depth:.3e} m inside the surface, less than \
                     one cell width ({h:.3e} m), so it is less than about two cells thick"
                ),
            );
        }
    }
    Ok(())
}

/// The cut bodies whose whole surface is a wall, by the name of their
/// wall boundary: the exact first-cell distance of the wall functions and
/// the exact wall distance (L13 stage 2c, D229, D108).
#[derive(Debug, Clone)]
pub struct WallBodies {
    bodies: Vec<Body>,
    /// `[immersed.cut] wall_sample`.
    sample: bool,
    /// Wall boundary name → body.
    by_boundary: std::collections::BTreeMap<String, usize>,
}

impl WallBodies {
    /// The bodies of a cut case whose surfaces are walls (`None` without a
    /// cut or without such a body).
    pub fn of(case: &crate::case::FlowCase) -> Result<Option<Self>> {
        if !case.immersed.is_cut() {
            return Ok(None);
        }
        let bodies = bodies_of(case)?;
        let mut by_boundary = std::collections::BTreeMap::new();
        for (i, b) in bodies.iter().enumerate() {
            // a `field` body has a level, not a distance: its wall faces
            // keep the distance to the face
            if matches!(b.kind, Kind::Field(_)) {
                continue;
            }
            if let Some(r) = solid_region_of(case, &b.name) {
                // a solid region: its surface is the interface wall of the
                // fluid view (stage 2d)
                by_boundary.insert(crate::cht::interface_name(r), i);
            } else {
                for n in wall_names(b) {
                    if matches!(
                        case.boundary.get(&n),
                        Some(crate::bc::BoundaryCondition::Wall { .. })
                    ) {
                        by_boundary.insert(n, i);
                    }
                }
            }
        }
        let sample = case
            .immersed
            .cut
            .as_ref()
            .is_none_or(CutSection::wall_sample);
        Ok((!by_boundary.is_empty()).then_some(Self {
            bodies,
            sample,
            by_boundary,
        }))
    }

    /// The wall functions of the cut walls sample the cell one layer out
    /// (`[immersed.cut] wall_sample`, default true).
    pub fn sample(&self) -> bool {
        self.sample
    }

    /// True if `boundary` is the wall of a cut body.
    pub fn is_body(&self, boundary: &str) -> bool {
        self.by_boundary.contains_key(boundary)
    }

    /// The exact distance of `x` to the surface of the body whose wall
    /// boundary is `boundary` (closest point, D108).
    pub fn distance_to(&self, boundary: &str, x: [f64; 3]) -> Option<f64> {
        let &i = self.by_boundary.get(boundary)?;
        Some(self.bodies[i].distance(x).abs())
    }

    /// The exact distance of every point to the nearest wall body.
    pub fn distance(&self, points: &[[f64; 3]]) -> Vec<f64> {
        use rayon::prelude::*;
        points
            .par_iter()
            .map(|&x| {
                self.by_boundary
                    .values()
                    .map(|&i| self.bodies[i].distance(x).abs())
                    .fold(f64::INFINITY, f64::min)
            })
            .collect()
    }
}

/// The momentum faces of the cut walls under a wall function (L13 stage
/// 2c): the wall shear stress of the wall function, ν_w = τ_w y_P/U_P, is
/// carried to the cut cell only when its wall faces take the exact y_P as
/// the distance, μ_w |S| (u_w − u_P)/y_P, and no non-orthogonal
/// correction: on the facets of a cut cell the correction k · ∇u_P picks up
/// the steep wall-normal gradient of the cell (its centroid does not lie
/// over the facet centroids) and moved the wall stress of single faces by
/// up to ±50 % on the tilted plate. Laminar flows and the low-Re
/// treatment keep the face gradient with its correction.
pub fn wall_function_faces(
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    disc: &mut foehncraft_fv::disc::Discretisation,
) -> Result<()> {
    let t = &case.turbulence;
    if t.is_laminar() || t.treatment() == crate::turbulence::WallTreatment::LowRe {
        return Ok(());
    }
    let Some(wb) = WallBodies::of(case)? else {
        return Ok(());
    };
    for p in &mesh.boundaries {
        if !wb.is_body(&p.name) {
            continue;
        }
        for f in p.range() {
            let c = mesh.owner[f] as usize;
            let Some(d) = wb.distance_to(&p.name, geo.cell_centroid[c]) else {
                continue;
            };
            let a = norm(geo.face_area[f]);
            let y = d.max(1e-3 * geo.cell_volume[c] / a);
            disc.dgeo[f] = a / y;
            disc.k.x[f] = 0.0;
            disc.k.y[f] = 0.0;
            disc.k.z[f] = 0.0;
        }
    }
    Ok(())
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// The findings of the surface check of the STL bodies (IMM002 to IMM009),
/// for the log.
fn surface_notes(bs: &[Body]) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    for b in bs {
        for f in &b.findings {
            out.push(Diagnostic::new(
                f.code,
                format!("immersed body {}: {}", b.name, f.message),
            ));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cyl(r: f64) -> ImmersedSection {
        ImmersedSection {
            method: None,
            cut: None,
            permeability: Some(1e-8),
            samples: None,
            refine: None,
            body: vec![BodySection {
                name: "c".into(),
                shape: Shape::Cylinder,
                centre: Some([0.5, 0.5, 0.0]),
                radius: Some(r),
                axis: None,
                normal: None,
                length: None,
                lower: None,
                upper: None,
                fraction: None,
                file: None,
                scale: None,
                invert: None,
                omega: None,
                patches: None,
            }],
            topology: None,
        }
    }

    #[test]
    fn rule_integrates_tetrahedron_moments() {
        // the mean of the collapsed coordinates over the unit tetrahedron
        // (0,0,0) (1,0,0) (1,1,0) (1,1,1): its centroid (3/4, 1/2, 1/4)
        for n in [2, 4, 8] {
            let q = rule(n);
            let m: [f64; 3] = std::array::from_fn(|i| q.iter().map(|(z, w)| z[i] * w).sum());
            let tol = 2.0 / (n * n) as f64;
            assert!((m[0] - 0.75).abs() < tol, "{n} {m:?}");
            assert!((m[1] - 0.5).abs() < tol, "{n} {m:?}");
            assert!((m[2] - 0.25).abs() < tol, "{n} {m:?}");
        }
    }

    #[test]
    fn checks_refuse_bad_bodies() {
        let mut s = cyl(0.1);
        assert!(bodies(&s).is_ok());
        s.permeability = None;
        assert!(bodies(&s).is_err());
        s.permeability = Some(1e-8);
        s.body[0].radius = Some(0.0);
        assert!(bodies(&s).is_err());
        s.body[0].radius = Some(0.1);
        s.body[0].axis = Some([0.0; 3]);
        assert!(bodies(&s).is_err());
        s.body[0].axis = None;
        s.body.push(s.body[0].clone());
        assert!(bodies(&s).is_err());
        s.body.clear();
        assert!(bodies(&s).is_err());
    }

    #[test]
    fn solid_volume_of_a_cylinder_converges() {
        // a 2D box mesh with an endless cylinder of radius 0.2: the solid
        // volume on the mesh against π r² (depth of the mesh)
        let spec: foehncraft_mesh::generate::BoxSpec = toml::from_str(
            "dimension = 2\nlower = [0.0, 0.0, 0.0]\nupper = [1.0, 1.0, 0.1]\ncells = [20, 20, 1]",
        )
        .expect("spec");
        let mesh = spec.generate().expect("mesh");
        let geo = Geometry::compute(&mesh);
        let s = cyl(0.2);
        let bs = bodies(&s).expect("bodies");
        let exact = bs[0]
            .exact_volume(mesh.bounding_box().0, mesh.bounding_box().1)
            .expect("exact");
        let mut errs = Vec::new();
        for n in [2, 4, 8] {
            let x = fractions(&bs, n, &mesh, &geo);
            let v: f64 = x[0].iter().zip(&geo.cell_volume).map(|(a, v)| a * v).sum();
            errs.push(((v - exact) / exact).abs());
            assert!(x[0].iter().all(|a| (0.0..=1.0).contains(a)));
        }
        assert!(errs[2] < 5e-3, "{errs:?}");
        assert!(errs[2] < errs[0], "{errs:?}");
    }
}
