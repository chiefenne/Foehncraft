//! The structural case file (TOML): materials, sections, beams, supports,
//! clamp springs (D290), springs, masses, loads and the analyses; the
//! model built from it and the results written for it.
//!
//! ```toml
//! [analysis]
//! static = true      # solve K u = f
//! modes = 5          # the lowest modes (0: none)
//!
//! [materials.spring_steel]
//! young = 206e9
//! poisson = 0.3
//! density = 7850.0
//!
//! [sections.reed]
//! shape = "rectangle"
//! width = 6e-3       # along the local y axis
//! thickness = 0.2e-3 # along the local z axis
//!
//! [[beams]]
//! name = "reed"
//! from = [0.0, 0.0, 0.0]
//! to = [0.03, 0.0, 0.0]
//! elements = 40
//! section = "reed"
//! material = "spring_steel"
//!
//! [[clamps]]           # springs per length along the clamp line
//! at = [0.0, 0.0, 0.0]
//! normal = "z"
//! line = "y"
//! k_rotation = 40.0    # N·m/rad per metre of clamp line; absent: rigid
//! ```

use crate::analysis::{self, Modal, Static};
use crate::beam::{BeamProps, Theory};
use crate::model::{DOFS, LineLoad, Model};
use crate::{Error, codes};
use serde::Deserialize;
use serde_json::{Value, json};
use std::collections::BTreeMap;

/// Format of `result.json`.
pub const RESULT_FORMAT: &str = "foehncraft.structure-result";
/// Version of the result format.
pub const RESULT_VERSION: u32 = 1;

fn default_true() -> bool {
    true
}

/// `[structure]`.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Header {
    /// Name of the case (default: the file name).
    #[serde(default)]
    pub name: Option<String>,
}

/// `[analysis]`.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnalysisKeys {
    /// Solve the static problem K u = f.
    #[serde(default, rename = "static")]
    pub statics: bool,
    /// Number of modes (0: no modal analysis).
    #[serde(default)]
    pub modes: usize,
    /// Shift σ of the eigenproblem (rad²/s²): the modes closest to it.
    /// Default 0 (the lowest modes); negative for a free structure.
    #[serde(default)]
    pub shift: f64,
}

/// `[materials.<name>]`: isotropic linear elastic.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Material {
    /// Young's modulus (Pa).
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Density (kg/m³).
    pub density: f64,
}

/// Shape of a section.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Shape {
    /// Width (local y) and thickness (local z).
    Rectangle,
    /// Diameter.
    Circle,
    /// Area, second moments, torsion constant and shear coefficients given.
    General,
}

/// `[sections.<name>]`.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Section {
    /// Shape.
    pub shape: Shape,
    /// Rectangle: width along local y (m).
    #[serde(default)]
    pub width: Option<f64>,
    /// Rectangle: thickness along local z (m).
    #[serde(default)]
    pub thickness: Option<f64>,
    /// Circle: diameter (m).
    #[serde(default)]
    pub diameter: Option<f64>,
    /// General: area (m²).
    #[serde(default)]
    pub area: Option<f64>,
    /// General: I_y (m⁴).
    #[serde(default)]
    pub iy: Option<f64>,
    /// General: I_z (m⁴).
    #[serde(default)]
    pub iz: Option<f64>,
    /// General: torsion constant (m⁴).
    #[serde(default)]
    pub torsion: Option<f64>,
    /// Shear coefficient for shear along local y (default: Cowper 1966
    /// for rectangles and circles, 5/6 otherwise).
    #[serde(default)]
    pub kappa_y: Option<f64>,
    /// Shear coefficient for shear along local z.
    #[serde(default)]
    pub kappa_z: Option<f64>,
}

/// `[[beams]]`: a straight line divided into equal elements.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BeamLine {
    /// Name (for line loads and results).
    pub name: String,
    /// Start point.
    pub from: [f64; 3],
    /// End point.
    pub to: [f64; 3],
    /// Number of elements.
    pub elements: usize,
    /// Section name.
    pub section: String,
    /// Material name.
    pub material: String,
    /// Beam theory (default Timoshenko).
    #[serde(default)]
    pub theory: Option<Theory>,
    /// Direction of the local y axis (the section's width); default: local
    /// z as close to global z as possible.
    #[serde(default)]
    pub y_axis: Option<[f64; 3]>,
}

/// Element of a plate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Default)]
pub enum ShellElement {
    /// Four-node MITC4 (Dvorkin & Bathe 1984).
    #[default]
    #[serde(rename = "mitc4")]
    Mitc4,
    /// Three-node MITC3+ (Lee, Lee & Bathe 2014), two per cell.
    #[serde(rename = "mitc3+")]
    Mitc3Plus,
}

/// `[[plates]]`: a flat parallelogram `origin + a·edge_a + b·edge_b`
/// (a, b ∈ [0, 1]) of shell elements; the normal is edge_a × edge_b.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Plate {
    /// Name (for loads and results).
    pub name: String,
    /// A corner.
    pub origin: [f64; 3],
    /// First edge from the corner.
    pub edge_a: [f64; 3],
    /// Second edge from the corner.
    pub edge_b: [f64; 3],
    /// Cells along each edge.
    pub divisions: [usize; 2],
    /// Thickness (m).
    pub thickness: f64,
    /// Material name.
    pub material: String,
    /// Element (`mitc4`, default, or `mitc3+`).
    #[serde(default)]
    pub element: ShellElement,
}

/// `[[pressures]]`: a uniform load per area on a plate.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pressure {
    /// Name of the plate.
    pub plate: String,
    /// Pressure along the plate's normal (Pa); positive pushes along it.
    #[serde(default)]
    pub value: f64,
    /// Force per area, global axes (N/m²), added to the pressure.
    #[serde(default)]
    pub force_per_area: [f64; 3],
}

/// An unknown of a node by name.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum DofName {
    /// Every unknown.
    All,
    /// Translations.
    Ux,
    /// .
    Uy,
    /// .
    Uz,
    /// Rotations.
    Rx,
    /// .
    Ry,
    /// .
    Rz,
}

/// `[[supports]]`: fixed unknowns at a node or along a line of nodes.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Support {
    /// The node (a point of the mesh).
    #[serde(default)]
    pub at: Option<[f64; 3]>,
    /// Every node on the segment between two points (a shell edge).
    #[serde(default)]
    pub along: Option<[[f64; 3]; 2]>,
    /// Fixed unknowns.
    pub fix: Vec<DofName>,
}

/// A global axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Axis {
    /// x.
    X,
    /// y.
    Y,
    /// z.
    Z,
}

impl Axis {
    fn index(self) -> usize {
        self as usize
    }
}

/// `[[clamps]]`: the mounting of a reed (D290). The clamp holds every
/// unknown of the node except the translation along `normal` and the
/// rotation about `line`, which get springs per unit length of the clamp
/// line (`width` long); a spring left out is rigid, and the result says so.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Clamp {
    /// The clamped node (a beam).
    #[serde(default)]
    pub at: Option<[f64; 3]>,
    /// The clamp line (a shell edge): every node on the segment, the
    /// springs per length shared by the lengths each node stands for.
    #[serde(default)]
    pub along: Option<[[f64; 3]; 2]>,
    /// Direction of the reed's lift.
    pub normal: Axis,
    /// Direction of the clamp line (across the reed).
    pub line: Axis,
    /// Length of the clamp line (m) of a beam clamp; default the width of
    /// the beam's rectangular section. A clamp `along` an edge has the
    /// length of its segment.
    #[serde(default)]
    pub width: Option<f64>,
    /// Translational stiffness along `normal` per metre of clamp line
    /// (N/m per m); absent: rigid.
    #[serde(default)]
    pub k_translation: Option<f64>,
    /// Rotational stiffness about `line` per metre of clamp line
    /// (N·m/rad per m); absent: rigid.
    #[serde(default)]
    pub k_rotation: Option<f64>,
}

/// `[[springs]]`: springs from a node to ground, global axes.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Spring {
    /// The node.
    pub at: [f64; 3],
    /// Translational stiffnesses (N/m).
    #[serde(default)]
    pub translation: [f64; 3],
    /// Rotational stiffnesses (N·m/rad).
    #[serde(default)]
    pub rotation: [f64; 3],
}

/// `[[masses]]`: a point mass.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointMass {
    /// The node.
    pub at: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
}

/// `[[loads]]`: a force and a moment at a node.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Load {
    /// The node.
    pub at: [f64; 3],
    /// Force (N).
    #[serde(default)]
    pub force: [f64; 3],
    /// Moment (N·m).
    #[serde(default)]
    pub moment: [f64; 3],
}

/// `[[line_loads]]`: a uniform force per length on a beam line.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LineLoadKeys {
    /// Name of the beam line.
    pub beam: String,
    /// Force per length, global axes (N/m).
    pub force_per_length: [f64; 3],
}

/// A structural case file.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    /// `[structure]`.
    #[serde(default)]
    pub structure: Header,
    /// `[analysis]`.
    #[serde(default)]
    pub analysis: AnalysisKeys,
    /// `[materials]`.
    #[serde(default)]
    pub materials: BTreeMap<String, Material>,
    /// `[sections]`.
    #[serde(default)]
    pub sections: BTreeMap<String, Section>,
    /// `[[beams]]`.
    #[serde(default)]
    pub beams: Vec<BeamLine>,
    /// `[[plates]]`.
    #[serde(default)]
    pub plates: Vec<Plate>,
    /// `[[pressures]]`.
    #[serde(default)]
    pub pressures: Vec<Pressure>,
    /// `[[supports]]`.
    #[serde(default)]
    pub supports: Vec<Support>,
    /// `[[clamps]]`.
    #[serde(default)]
    pub clamps: Vec<Clamp>,
    /// `[[springs]]`.
    #[serde(default)]
    pub springs: Vec<Spring>,
    /// `[[masses]]`.
    #[serde(default)]
    pub masses: Vec<PointMass>,
    /// `[[loads]]`.
    #[serde(default)]
    pub loads: Vec<Load>,
    /// `[[line_loads]]`.
    #[serde(default)]
    pub line_loads: Vec<LineLoadKeys>,
    /// `[points]`: named points whose results are reported.
    #[serde(default)]
    pub points: BTreeMap<String, [f64; 3]>,
    /// `[output]`.
    #[serde(default)]
    pub output: OutputKeys,
}

/// `[output]`.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputKeys {
    /// Write the mesh with displacements and mode shapes (`structure.vtk`).
    #[serde(default = "default_true")]
    pub vtk: bool,
}

impl Default for OutputKeys {
    fn default() -> Self {
        Self { vtk: true }
    }
}

fn bad(msg: impl Into<String>) -> Error {
    Error::new(codes::CASE, msg.into())
}

impl Case {
    /// Reads a case from TOML text.
    pub fn from_toml(text: &str) -> Result<Self, Error> {
        toml::from_str(text).map_err(|e| bad(format!("structural case: {e}")))
    }
}

/// How a clamp was taken (D290: the result says when it was rigid).
#[derive(Debug, Clone, serde::Serialize)]
pub struct ClampInfo {
    /// The clamped nodes.
    pub nodes: Vec<usize>,
    /// Length of the clamp line (m).
    pub width: f64,
    /// Translation along the normal: total stiffness (N/m), `None` rigid.
    pub k_translation: Option<f64>,
    /// Rotation about the line: total stiffness (N·m/rad), `None` rigid.
    pub k_rotation: Option<f64>,
}

/// The model of a case with what the results need.
#[derive(Debug, Clone)]
pub struct Built {
    /// The model.
    pub model: Model,
    /// Elements of each beam line, by name.
    pub lines: BTreeMap<String, std::ops::Range<usize>>,
    /// Shell elements of each plate, by name.
    pub plates: BTreeMap<String, std::ops::Range<usize>>,
    /// The clamps.
    pub clamps: Vec<ClampInfo>,
    /// Named points and their nodes.
    pub points: BTreeMap<String, usize>,
}

fn positive(v: Option<f64>, what: &str) -> Result<f64, Error> {
    match v {
        Some(x) if x > 0.0 && x.is_finite() => Ok(x),
        Some(x) => Err(bad(format!("{what} must be positive (is {x})"))),
        None => Err(bad(format!("{what} is missing"))),
    }
}

fn section_props(
    name: &str,
    s: &Section,
    m: &Material,
    theory: Theory,
) -> Result<BeamProps, Error> {
    let nu = m.poisson;
    if !(m.young > 0.0 && m.density > 0.0 && nu > -1.0 && nu < 0.5) {
        return Err(bad(format!(
            "material of section `{name}`: needs young > 0, density > 0 and −1 < poisson < 0.5"
        )));
    }
    let shear = m.young / (2.0 * (1.0 + nu));
    let (area, iy, iz, torsion, ky, kz) = match s.shape {
        Shape::Rectangle => {
            let b = positive(s.width, &format!("[sections.{name}] width"))?;
            let h = positive(s.thickness, &format!("[sections.{name}] thickness"))?;
            // Saint-Venant torsion of a rectangle (Roark's Formulas, table 10.7)
            let (long, short) = if b >= h { (b, h) } else { (h, b) };
            let r = short / long;
            let j = long * short.powi(3) * (1.0 / 3.0 - 0.21 * r * (1.0 - r.powi(4) / 12.0));
            let k = 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu); // Cowper 1966
            (b * h, b * h.powi(3) / 12.0, h * b.powi(3) / 12.0, j, k, k)
        }
        Shape::Circle => {
            let d = positive(s.diameter, &format!("[sections.{name}] diameter"))?;
            let pi = std::f64::consts::PI;
            let i = pi * d.powi(4) / 64.0;
            let k = 6.0 * (1.0 + nu) / (7.0 + 6.0 * nu); // Cowper 1966
            (pi * d * d / 4.0, i, i, 2.0 * i, k, k)
        }
        Shape::General => (
            positive(s.area, &format!("[sections.{name}] area"))?,
            positive(s.iy, &format!("[sections.{name}] iy"))?,
            positive(s.iz, &format!("[sections.{name}] iz"))?,
            positive(s.torsion, &format!("[sections.{name}] torsion"))?,
            5.0 / 6.0,
            5.0 / 6.0,
        ),
    };
    let kappa_y = s.kappa_y.unwrap_or(ky);
    let kappa_z = s.kappa_z.unwrap_or(kz);
    if !(kappa_y > 0.0 && kappa_z > 0.0) {
        return Err(bad(format!(
            "[sections.{name}] shear coefficients must be positive"
        )));
    }
    Ok(BeamProps {
        young: m.young,
        shear,
        density: m.density,
        area,
        iy,
        iz,
        torsion,
        kappa_y,
        kappa_z,
        theory,
    })
}

fn dofs_of(names: &[DofName]) -> Vec<usize> {
    let mut v = Vec::new();
    for n in names {
        match n {
            DofName::All => v.extend(0..DOFS),
            DofName::Ux => v.push(0),
            DofName::Uy => v.push(1),
            DofName::Uz => v.push(2),
            DofName::Rx => v.push(3),
            DofName::Ry => v.push(4),
            DofName::Rz => v.push(5),
        }
    }
    v
}

/// Builds the model of a case.
pub fn build(case: &Case) -> Result<Built, Error> {
    let mut model = Model::default();
    let mut lines = BTreeMap::new();
    let mut widths: Vec<Option<f64>> = Vec::new(); // rectangle width per node
    if case.beams.is_empty() && case.plates.is_empty() {
        return Err(bad("the case has no `[[beams]]` and no `[[plates]]`"));
    }
    // merge tolerance: end points closer than this are one node
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for b in &case.beams {
        for p in [b.from, b.to] {
            for d in 0..3 {
                lo[d] = lo[d].min(p[d]);
                hi[d] = hi[d].max(p[d]);
            }
        }
    }
    let size = crate::beam::norm(crate::beam::sub(hi, lo));
    let tol = 1e-9 * size.max(f64::MIN_POSITIVE);
    let mut ends: Vec<usize> = Vec::new();
    let node_at =
        |model: &mut Model, widths: &mut Vec<Option<f64>>, ends: &mut Vec<usize>, x: [f64; 3]| {
            for &e in ends.iter() {
                if crate::beam::norm(crate::beam::sub(model.nodes[e], x)) <= tol {
                    return e;
                }
            }
            let n = model.add_node(x);
            widths.push(None);
            ends.push(n);
            n
        };
    let mut seen = std::collections::BTreeSet::new();
    for b in &case.beams {
        if !seen.insert(b.name.clone()) {
            return Err(bad(format!("two `[[beams]]` are named `{}`", b.name)));
        }
        if b.elements == 0 {
            return Err(bad(format!(
                "beam `{}`: `elements` must be at least 1",
                b.name
            )));
        }
        let sec = case.sections.get(&b.section).ok_or_else(|| {
            bad(format!(
                "beam `{}`: no section `{}` in [sections]",
                b.name, b.section
            ))
        })?;
        let mat = case.materials.get(&b.material).ok_or_else(|| {
            bad(format!(
                "beam `{}`: no material `{}` in [materials]",
                b.name, b.material
            ))
        })?;
        let props = section_props(&b.section, sec, mat, b.theory.unwrap_or(Theory::Timoshenko))?;
        model.props.push(props);
        let pi = model.props.len() - 1;
        let width = (sec.shape == Shape::Rectangle)
            .then_some(sec.width)
            .flatten();
        let a = node_at(&mut model, &mut widths, &mut ends, b.from);
        let z = node_at(&mut model, &mut widths, &mut ends, b.to);
        if a == z {
            return Err(bad(format!("beam `{}` has no length", b.name)));
        }
        let first = model.beams.len();
        let mut prev = a;
        for i in 1..=b.elements {
            let next = if i == b.elements {
                z
            } else {
                let t = i as f64 / b.elements as f64;
                let x = [
                    b.from[0] + t * (b.to[0] - b.from[0]),
                    b.from[1] + t * (b.to[1] - b.from[1]),
                    b.from[2] + t * (b.to[2] - b.from[2]),
                ];
                widths.push(None);
                model.add_node(x)
            };
            model.add_beam(prev, next, pi, b.y_axis).ok_or_else(|| {
                bad(format!(
                    "beam `{}`: `y_axis` is parallel to the beam (or zero)",
                    b.name
                ))
            })?;
            for n in [prev, next] {
                if widths[n].is_none() {
                    widths[n] = width;
                }
            }
            prev = next;
        }
        lines.insert(b.name.clone(), first..model.beams.len());
    }
    let mut plates = BTreeMap::new();
    for p in &case.plates {
        if lines.contains_key(&p.name) || !seen.insert(p.name.clone()) {
            return Err(bad(format!("two parts are named `{}`", p.name)));
        }
        let mat = case.materials.get(&p.material).ok_or_else(|| {
            bad(format!("plate `{}`: no material `{}` in [materials]", p.name, p.material))
        })?;
        if !(mat.young > 0.0 && mat.density > 0.0 && mat.poisson > -1.0 && mat.poisson < 0.5) {
            return Err(bad(format!(
                "plate `{}`: the material needs young > 0, density > 0 and −1 < poisson < 0.5",
                p.name
            )));
        }
        if p.thickness.is_nan() || p.thickness <= 0.0 {
            return Err(bad(format!("plate `{}`: thickness must be positive", p.name)));
        }
        let [na, nb] = p.divisions;
        if na == 0 || nb == 0 {
            return Err(bad(format!("plate `{}`: `divisions` must be at least 1", p.name)));
        }
        let n = crate::beam::cross(p.edge_a, p.edge_b);
        if crate::beam::norm(n) <= 1e-12 * crate::beam::norm(p.edge_a) * crate::beam::norm(p.edge_b) {
            return Err(bad(format!("plate `{}`: edge_a and edge_b are parallel or zero", p.name)));
        }
        model.shell_props.push(crate::shell::ShellProps {
            young: mat.young,
            poisson: mat.poisson,
            density: mat.density,
            thickness: p.thickness,
            kappa: 5.0 / 6.0,
        });
        let pi = model.shell_props.len() - 1;
        let first = model.shells.len();
        let at = |a: f64, b: f64| {
            [
                p.origin[0] + a * p.edge_a[0] + b * p.edge_b[0],
                p.origin[1] + a * p.edge_a[1] + b * p.edge_b[1],
                p.origin[2] + a * p.edge_a[2] + b * p.edge_b[2],
            ]
        };
        let before = model.nodes.len();
        match p.element {
            ShellElement::Mitc4 => model.add_shell_grid(at, na, nb, pi),
            ShellElement::Mitc3Plus => model.add_tri_grid(at, na, nb, pi),
        };
        widths.extend(std::iter::repeat_n(None, model.nodes.len() - before));
        plates.insert(p.name.clone(), first..model.shells.len());
    }
    let size = model.size();
    let find = |model: &Model, x: [f64; 3], what: &str| -> Result<usize, Error> {
        let (n, d) = model.nearest(x).ok_or_else(|| bad("no nodes"))?;
        if d <= 1e-6 * size.max(f64::MIN_POSITIVE) {
            Ok(n)
        } else {
            let p = model.nodes[n];
            Err(bad(format!(
                "{what} at ({}, {}, {}) is not a node; the nearest node is at ({}, {}, {}), {d:.3e} m away",
                x[0], x[1], x[2], p[0], p[1], p[2]
            )))
        }
    };
    // the nodes on a segment, with the length each one stands for
    let along = |model: &Model, seg: [[f64; 3]; 2], what: &str| -> Result<Vec<(usize, f64)>, Error> {
        let d = crate::beam::sub(seg[1], seg[0]);
        let len = crate::beam::norm(d);
        if len.is_nan() || len <= 0.0 {
            return Err(bad(format!("{what}: `along` needs two different points")));
        }
        let tol = 1e-6 * size.max(f64::MIN_POSITIVE);
        let mut hits: Vec<(f64, usize)> = Vec::new();
        for (i, x) in model.nodes.iter().enumerate() {
            let r = crate::beam::sub(*x, seg[0]);
            let t = crate::beam::dot3(r, d) / (len * len);
            if !(-1e-9..=1.0 + 1e-9).contains(&t) {
                continue;
            }
            let off = crate::beam::sub(r, crate::beam::scale(d, t));
            if crate::beam::norm(off) <= tol {
                hits.push((t, i));
            }
        }
        if hits.len() < 2 {
            return Err(bad(format!(
                "{what}: fewer than two nodes lie on the segment from {:?} to {:?}",
                seg[0], seg[1]
            )));
        }
        hits.sort_by(|a, b| a.0.total_cmp(&b.0));
        let m = hits.len();
        Ok((0..m)
            .map(|k| {
                let lo = if k == 0 { hits[0].0 } else { 0.5 * (hits[k - 1].0 + hits[k].0) };
                let hi = if k + 1 == m { hits[m - 1].0 } else { 0.5 * (hits[k].0 + hits[k + 1].0) };
                (hits[k].1, (hi - lo) * len)
            })
            .collect())
    };
    let place = |model: &Model, at: Option<[f64; 3]>, seg: Option<[[f64; 3]; 2]>, what: &str| -> Result<Vec<(usize, f64)>, Error> {
        match (at, seg) {
            (Some(x), None) => Ok(vec![(find(model, x, what)?, 0.0)]),
            (None, Some(s)) => along(model, s, what),
            _ => Err(bad(format!("{what}: give `at` (a node) or `along` (a segment), not both"))),
        }
    };
    for s in &case.supports {
        for (n, _) in place(&model, s.at, s.along, "[[supports]]")? {
            for d in dofs_of(&s.fix) {
                model.fixed[DOFS * n + d] = true;
            }
        }
    }
    let mut clamps = Vec::new();
    for c in &case.clamps {
        let nodes = place(&model, c.at, c.along, "[[clamps]]")?;
        if c.normal == c.line {
            return Err(bad("[[clamps]]: `normal` and `line` must be different axes"));
        }
        let (width, share): (f64, Vec<f64>) = if c.along.is_some() {
            if c.width.is_some() {
                return Err(bad("[[clamps]]: a clamp `along` an edge takes its length from the segment; leave out `width`"));
            }
            (nodes.iter().map(|x| x.1).sum(), nodes.iter().map(|x| x.1).collect())
        } else {
            let n = nodes[0].0;
            match c.width.or(widths[n]) {
                Some(w) if w > 0.0 => (w, vec![w]),
                _ => {
                    return Err(bad(
                        "[[clamps]]: `width` (the length of the clamp line) is needed when the beam has no rectangular section",
                    ));
                }
            }
        };
        let tn = c.normal.index();
        let rl = 3 + c.line.index();
        for k in [c.k_translation, c.k_rotation].into_iter().flatten() {
            if k.is_nan() || k < 0.0 || !k.is_finite() {
                return Err(bad(format!("[[clamps]]: spring stiffness must be ≥ 0 (is {k})")));
            }
        }
        for (&(n, _), &l) in nodes.iter().zip(&share) {
            for d in 0..DOFS {
                let spring = if d == tn {
                    c.k_translation
                } else if d == rl {
                    c.k_rotation
                } else {
                    None
                };
                match spring {
                    Some(k) => {
                        model.fixed[DOFS * n + d] = false;
                        model.springs.push((DOFS * n + d, k * l));
                    }
                    None => model.fixed[DOFS * n + d] = true,
                }
            }
        }
        clamps.push(ClampInfo {
            nodes: nodes.iter().map(|x| x.0).collect(),
            width,
            k_translation: c.k_translation.map(|k| k * width),
            k_rotation: c.k_rotation.map(|k| k * width),
        });
    }
    for s in &case.springs {
        let n = find(&model, s.at, "[[springs]]")?;
        for d in 0..3 {
            for (off, k) in [(0, s.translation[d]), (3, s.rotation[d])] {
                if k < 0.0 || !k.is_finite() {
                    return Err(bad(format!("[[springs]]: stiffness must be ≥ 0 (is {k})")));
                }
                if k > 0.0 {
                    model.springs.push((DOFS * n + off + d, k));
                }
            }
        }
    }
    for m in &case.masses {
        let n = find(&model, m.at, "[[masses]]")?;
        if m.mass.is_nan() || m.mass < 0.0 {
            return Err(bad("[[masses]]: mass must be ≥ 0"));
        }
        model.masses.push((n, m.mass));
    }
    for l in &case.loads {
        let n = find(&model, l.at, "[[loads]]")?;
        for d in 0..3 {
            model.forces[DOFS * n + d] += l.force[d];
            model.forces[DOFS * n + 3 + d] += l.moment[d];
        }
    }
    for pr in &case.pressures {
        let r = plates
            .get(&pr.plate)
            .ok_or_else(|| bad(format!("[[pressures]]: no plate `{}`", pr.plate)))?;
        for e in r.clone() {
            model.area_loads.push(crate::model::AreaLoad {
                element: e,
                q: pr.force_per_area,
                pressure: pr.value,
            });
        }
    }
    for l in &case.line_loads {
        let r = lines
            .get(&l.beam)
            .ok_or_else(|| bad(format!("[[line_loads]]: no beam `{}`", l.beam)))?;
        for e in r.clone() {
            model.line_loads.push(LineLoad {
                element: e,
                q: l.force_per_length,
            });
        }
    }
    let mut points = BTreeMap::new();
    for (name, x) in &case.points {
        points.insert(name.clone(), find(&model, *x, &format!("[points] {name}"))?);
    }
    Ok(Built {
        model,
        lines,
        plates,
        clamps,
        points,
    })
}

/// What a run computed.
#[derive(Debug, Clone)]
pub struct Outcome {
    /// The model.
    pub built: Built,
    /// The static result.
    pub statics: Option<Static>,
    /// The modal result.
    pub modal: Option<Modal>,
}

/// Runs the analyses of a case.
pub fn run(case: &Case) -> Result<Outcome, Error> {
    let built = build(case)?;
    if !case.analysis.statics && case.analysis.modes == 0 {
        return Err(bad(
            "[analysis] asks for nothing: set `static = true` or `modes`",
        ));
    }
    let statics = if case.analysis.statics {
        Some(analysis::solve_static(&built.model)?)
    } else {
        None
    };
    let modal = if case.analysis.modes > 0 {
        Some(analysis::solve_modal(
            &built.model,
            case.analysis.modes,
            case.analysis.shift,
        )?)
    } else {
        None
    };
    Ok(Outcome {
        built,
        statics,
        modal,
    })
}

fn node_vec(u: &[f64], n: usize, off: usize) -> [f64; 3] {
    [
        u[DOFS * n + off],
        u[DOFS * n + off + 1],
        u[DOFS * n + off + 2],
    ]
}

impl Outcome {
    /// The result as JSON (`result.json`).
    pub fn to_json(&self, name: &str) -> Value {
        let m = &self.built.model;
        let free = m.fixed.iter().filter(|f| !**f).count();
        let mut warnings = Vec::new();
        let mut out = json!({
            "format": RESULT_FORMAT,
            "format_version": RESULT_VERSION,
            "name": name,
            "nodes": m.nodes.len(),
            "elements": m.beams.len() + m.shells.len(),
            "beam_elements": m.beams.len(),
            "shell_elements": m.shells.len(),
            "unknowns": m.dofs(),
            "free_unknowns": free,
            "mass": m.total_mass(),
            "clamps": self.built.clamps.iter().map(|c| json!({
                "nodes": c.nodes.len(),
                "at": m.nodes[c.nodes[0]],
                "width": c.width,
                "translation": match c.k_translation { Some(k) => json!({"rigid": false, "stiffness": k}), None => json!({"rigid": true}) },
                "rotation": match c.k_rotation { Some(k) => json!({"rigid": false, "stiffness": k}), None => json!({"rigid": true}) },
            })).collect::<Vec<_>>(),
        });
        if let Some(s) = &self.statics {
            let (mut best, mut at) = (0.0, 0);
            for n in 0..m.nodes.len() {
                let v = crate::beam::norm(node_vec(&s.u, n, 0));
                if v > best {
                    best = v;
                    at = n;
                }
            }
            let mut rf = [0.0; 3];
            for n in 0..m.nodes.len() {
                for (d, r) in rf.iter_mut().enumerate() {
                    *r += s.reactions[DOFS * n + d];
                }
            }
            let points: serde_json::Map<String, Value> = self
                .built
                .points
                .iter()
                .map(|(k, &n)| {
                    (
                        k.clone(),
                        json!({"node": n, "displacement": node_vec(&s.u, n, 0), "rotation": node_vec(&s.u, n, 3)}),
                    )
                })
                .collect();
            out["static"] = json!({
                "max_displacement": {"value": best, "node": at, "at": m.nodes[at]},
                "strain_energy": s.energy,
                "reaction_force": rf,
                "factor_entries": s.fill,
                "points": points,
            });
        }
        if let Some(md) = &self.modal {
            for w in &md.warnings {
                warnings.push(json!({"code": codes::MODES, "message": w}));
            }
            let points: serde_json::Map<String, Value> = self
                .built
                .points
                .iter()
                .map(|(k, &n)| {
                    (
                        k.clone(),
                        Value::Array(
                            md.modes
                                .iter()
                                .map(|mo| json!(node_vec(&mo.shape, n, 0)))
                                .collect(),
                        ),
                    )
                })
                .collect();
            out["modal"] = json!({
                "steps": md.steps,
                "converged": md.converged,
                "sturm": {"lo": md.sturm.lo, "hi": md.sturm.hi, "expected": md.sturm.expected, "found": md.sturm.found, "ok": md.sturm.ok()},
                "modes": md.modes.iter().enumerate().map(|(i, mo)| json!({
                    "mode": i + 1,
                    "frequency": mo.frequency(),
                    "omega": mo.lambda.max(0.0).sqrt(),
                    "effective_mass_fraction": mo.effective_mass.map(|e| e / md.total_mass),
                    "residual": mo.residual,
                })).collect::<Vec<_>>(),
                "points_mass_normalised": points,
            });
        }
        out["warnings"] = Value::Array(warnings);
        out
    }

    /// The mesh with the static displacement and the mode shapes, as a
    /// legacy VTK polydata file (ParaView, VisIt). Mode shapes are scaled
    /// to a largest translation of 1 (`mode_<n>`, and `mode_<n>_rotation`
    /// with the same factor).
    pub fn to_vtk(&self, name: &str) -> String {
        use std::fmt::Write;
        let m = &self.built.model;
        let mut s = String::new();
        let _ = writeln!(
            s,
            "# vtk DataFile Version 3.0\nFoehncraft structure {name}\nASCII\nDATASET POLYDATA"
        );
        let _ = writeln!(s, "POINTS {} double", m.nodes.len());
        for p in &m.nodes {
            let _ = writeln!(s, "{:.9e} {:.9e} {:.9e}", p[0], p[1], p[2]);
        }
        if !m.beams.is_empty() {
            let _ = writeln!(s, "LINES {} {}", m.beams.len(), 3 * m.beams.len());
            for b in &m.beams {
                let _ = writeln!(s, "2 {} {}", b.nodes[0], b.nodes[1]);
            }
        }
        if !m.shells.is_empty() {
            let size: usize = m.shells.iter().map(|e| e.used().len() + 1).sum();
            let _ = writeln!(s, "POLYGONS {} {size}", m.shells.len());
            for e in &m.shells {
                let nodes = e.used();
                let ids: Vec<String> = nodes.iter().map(|n| n.to_string()).collect();
                let _ = writeln!(s, "{} {}", nodes.len(), ids.join(" "));
            }
        }
        let _ = writeln!(s, "POINT_DATA {}", m.nodes.len());
        let vector = |s: &mut String, label: &str, u: &[f64], off: usize, scale: f64| {
            let _ = writeln!(s, "VECTORS {label} double");
            for n in 0..m.nodes.len() {
                let v = node_vec(u, n, off);
                let _ = writeln!(
                    s,
                    "{:.9e} {:.9e} {:.9e}",
                    v[0] * scale,
                    v[1] * scale,
                    v[2] * scale
                );
            }
        };
        if let Some(st) = &self.statics {
            vector(&mut s, "displacement", &st.u, 0, 1.0);
            vector(&mut s, "rotation", &st.u, 3, 1.0);
        }
        if let Some(md) = &self.modal {
            for (i, mo) in md.modes.iter().enumerate() {
                let big = (0..m.nodes.len())
                    .map(|n| crate::beam::norm(node_vec(&mo.shape, n, 0)))
                    .fold(0.0, f64::max);
                let sc = if big > 0.0 { 1.0 / big } else { 1.0 };
                vector(&mut s, &format!("mode_{}", i + 1), &mo.shape, 0, sc);
                vector(
                    &mut s,
                    &format!("mode_{}_rotation", i + 1),
                    &mo.shape,
                    3,
                    sc,
                );
            }
        }
        s
    }
}
