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
//! ```

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
    /// follows a curved body inside the cell (default true).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_point: Option<bool>,
}

impl CutSection {
    /// The merge threshold.
    pub fn min_volume_fraction(&self) -> f64 {
        self.min_volume_fraction.unwrap_or(0.3)
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
}

impl Body {
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
            Kind::Field(_) => None,
        }
    }

    /// True if `x` lies inside the body (the surface counts as inside; a
    /// `field` body where χ ≥ ½).
    #[inline]
    pub fn inside(&self, x: [f64; 3]) -> bool {
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
        }
    }

    /// The exact volume of the body inside the box `lo`–`hi`, where it is
    /// simple to give: a sphere or box inside it, a cylinder along a
    /// coordinate axis that spans it or lies inside it. None otherwise.
    pub fn exact_volume(&self, lo: [f64; 3], hi: [f64; 3]) -> Option<f64> {
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
            Kind::Field(_) => None,
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
        };
        out.push(Body {
            name: b.name.clone(),
            kind,
        });
    }
    Ok(out)
}

/// The settings that do not go with a penalised body (IMM102).
pub fn check_case(case: &crate::case::FlowCase) -> Result<()> {
    if case.immersed.is_empty() {
        return Ok(());
    }
    bodies_of(case)?;
    if case.immersed.is_cut() {
        let what = if !case.regions.is_empty() {
            Some("regions (conjugate heat transfer on cut cells is stage 2d)")
        } else if case.immersed.topology.is_some() {
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

/// The solid fraction χ of every local cell, one array per body.
pub fn fractions(bodies: &[Body], samples: usize, mesh: &Mesh, geo: &Geometry) -> Vec<Vec<f64>> {
    let nc = mesh.n_cells;
    let q = rule(samples);
    let mut out = vec![vec![0.0; nc]; bodies.len()];
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
                for (b, body) in bodies.iter().enumerate() {
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
        for (b, o) in out.iter_mut().enumerate() {
            o[c] = if vol > 0.0 {
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
    let mut notes = Vec::new();
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
/// distance, the crossing points are the exact roots on the edges (or
/// linear), the wall point is the Newton projection onto the nearest body.
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
            let d = b.distance(x).abs();
            if d < best.0 {
                best = (d, i);
            }
        }
        best.1
    }
}

impl foehncraft_mesh::cut::ImplicitBody for CutBodies<'_> {
    fn value(&self, x: [f64; 3]) -> f64 {
        self.bodies
            .iter()
            .map(|b| b.distance(x))
            .fold(f64::INFINITY, f64::min)
    }
    fn edge_point(&self, a: [f64; 3], b: [f64; 3], fa: f64, fb: f64) -> [f64; 3] {
        if self.exact {
            foehncraft_mesh::cut::edge_root(|x| self.value(x), a, b, fa, fb)
        } else {
            let t = fa / (fa - fb);
            [
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ]
        }
    }
    fn project(&self, x: [f64; 3]) -> Option<[f64; 3]> {
        let b = &self.bodies[self.nearest(x)];
        foehncraft_mesh::cut::project_newton(|y| b.distance(y), x, self.h)
    }
    fn patch(&self, x: [f64; 3]) -> usize {
        self.nearest(x)
    }
}

/// The cut of `[immersed] method = "cut"` (L13 stage 2, D227): the mesh
/// with the bodies cut out, their surfaces as wall boundaries named after
/// the bodies, and the notes for the log (IMM203, IMM204).
pub fn cut_mesh(case: &crate::case::FlowCase, mesh: &Mesh) -> Result<(Mesh, Vec<Diagnostic>)> {
    let bs = bodies_of(case)?;
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
    };
    let c = foehncraft_mesh::cut::cut(mesh, &body, &names, &opt)
        .map_err(|e| diag::err(Diagnostic::new(diag::IMM_CUT, format!("cut cells: {e}"))))?;
    let r = &c.report;
    let mut notes = vec![Diagnostic::new(
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
    )];
    for (b, a) in bs.iter().zip(&r.wall_area) {
        let mut line = format!("cut body {}: wall area {a:.6e} m²", b.name);
        if let Some(e) = b.exact_area(lo, hi) {
            line.push_str(&format!(
                ", exact {e:.6e} m² ({:+.3} %)",
                100.0 * (a - e) / e
            ));
        }
        notes.push(Diagnostic::new(diag::IMM_CUT_REPORT, line));
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
    Ok((c.mesh, notes))
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
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
            body: vec![BodySection {
                name: "c".into(),
                shape: Shape::Cylinder,
                centre: Some([0.5, 0.5, 0.0]),
                radius: Some(r),
                axis: None,
                length: None,
                lower: None,
                upper: None,
                fraction: None,
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
