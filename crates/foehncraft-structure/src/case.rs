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
use crate::dynamics::{GenAlpha, ModalBasis, Rayleigh};
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
    /// Linear thermal expansion coefficient α (1/K), for `[[temperatures]]`.
    #[serde(default)]
    pub expansion: Option<f64>,
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

/// Element of a generated solid block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum SolidElement {
    /// Twenty-node hexahedra.
    #[default]
    Hex20,
    /// Ten-node tetrahedra, six per hexahedron.
    Tet10,
}

/// A generated block: `origin + [0, size]` divided into `divisions`
/// hexahedra along x, y and z.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BlockKeys {
    /// The corner with the smallest coordinates.
    pub origin: [f64; 3],
    /// Edge lengths along x, y, z (m).
    pub size: [f64; 3],
    /// Cells along x, y, z.
    pub divisions: [usize; 3],
    /// Element: `hex20` (default) or `tet10`.
    #[serde(default)]
    pub element: SolidElement,
}

/// `[[solids]]`: a part of quadratic solid elements, from a Gmsh file of
/// second order (`mesh`) or a generated block (`block`).
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SolidKeys {
    /// Name (for loads, supports and results).
    pub name: String,
    /// A Gmsh `.msh` 4.1 file (relative to the case file): ten-node
    /// tetrahedra or twenty-node hexahedra; physical surfaces of six-node
    /// triangles or eight-node quadrilaterals name the surfaces.
    #[serde(default)]
    pub mesh: Option<String>,
    /// Physical volumes of the mesh that make the part (default: all).
    #[serde(default)]
    pub volumes: Vec<String>,
    /// A generated block (surfaces `xmin` … `zmax`).
    #[serde(default)]
    pub block: Option<BlockKeys>,
    /// Material name.
    pub material: String,
}

/// `[[contacts]]`: a rigid plane (a seat, a stop) that the nodes of a part
/// may touch but not cross (MS5 part 4).
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContactKeys {
    /// Name (for the results).
    pub name: String,
    /// The part whose nodes may touch it (a beam line, a plate or a
    /// solid).
    pub part: String,
    /// A solid's named surface: only its nodes (default: every node of the
    /// part).
    #[serde(default)]
    pub surface: Option<String>,
    /// A point of the plane.
    pub point: [f64; 3],
    /// The plane's normal, pointing to the side the nodes stay on.
    pub normal: [f64; 3],
    /// Distance the nodes keep from the plane (m): half the thickness of a
    /// plate or beam, so that its face touches (default 0).
    #[serde(default)]
    pub offset: f64,
}

/// `[contact]`: settings of the contact solution.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContactSettings {
    /// Penalty stiffness over the largest stiffness of a contact node
    /// along its plane's normal (default 1000).
    #[serde(default = "default_penalty")]
    pub penalty: f64,
    /// Largest penetration of a static solution, relative to the model
    /// size (default 10⁻⁹).
    #[serde(default = "default_contact_tol")]
    pub tolerance: f64,
}

fn default_penalty() -> f64 {
    1000.0
}

fn default_contact_tol() -> f64 {
    1e-9
}

impl Default for ContactSettings {
    fn default() -> Self {
        Self {
            penalty: default_penalty(),
            tolerance: default_contact_tol(),
        }
    }
}

/// `[[pressures]]`: a uniform load per area on a plate, or on a named
/// surface of a solid.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pressure {
    /// Name of the plate.
    #[serde(default)]
    pub plate: Option<String>,
    /// Name of the solid (with `surface`).
    #[serde(default)]
    pub part: Option<String>,
    /// Named surface of the solid.
    #[serde(default)]
    pub surface: Option<String>,
    /// Pressure (Pa). Plate: along the plate's normal, positive pushes
    /// along it. Solid: positive pushes into the solid (a fluid's
    /// pressure).
    #[serde(default)]
    pub value: f64,
    /// Force per area, global axes (N/m²), added to the pressure.
    #[serde(default)]
    pub force_per_area: [f64; 3],
    /// Time history that scales the load.
    #[serde(default)]
    pub history: Option<String>,
}

/// `[[temperatures]]`: a temperature change of a beam line or a plate
/// from its stress-free state, read as a thermal strain α ΔT (the
/// material's `expansion`).
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Temperature {
    /// Name of the beam line or plate.
    pub part: String,
    /// Temperature change of the axis or mid-surface (K).
    #[serde(default)]
    pub change: f64,
    /// Difference across the thickness (K): plates, top (along the
    /// normal) minus bottom; beams, +z minus −z of the section.
    #[serde(default)]
    pub difference: f64,
    /// Beams: difference across the width, +y minus −y (K).
    #[serde(default)]
    pub difference_y: f64,
    /// Time history that scales it.
    #[serde(default)]
    pub history: Option<String>,
}

/// `[histories.<name>]`: a factor in time, a table (linear between its
/// points, held after the last) or a sine.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct History {
    /// Points [time (s), factor].
    #[serde(default)]
    pub table: Option<Vec<[f64; 2]>>,
    /// Sine: amplitude.
    #[serde(default)]
    pub amplitude: Option<f64>,
    /// Sine: frequency (Hz).
    #[serde(default)]
    pub frequency: Option<f64>,
    /// Sine: phase (degrees).
    #[serde(default)]
    pub phase: f64,
    /// Sine: offset added.
    #[serde(default)]
    pub offset: f64,
}

impl History {
    /// The factor at time `t`.
    pub fn at(&self, t: f64) -> f64 {
        if let Some(tab) = &self.table {
            if t <= tab[0][0] {
                return tab[0][1];
            }
            for w in tab.windows(2) {
                if t <= w[1][0] {
                    let f = (t - w[0][0]) / (w[1][0] - w[0][0]);
                    return w[0][1] + f * (w[1][1] - w[0][1]);
                }
            }
            return tab[tab.len() - 1][1];
        }
        let a = self.amplitude.unwrap_or(0.0);
        let f = self.frequency.unwrap_or(0.0);
        self.offset + a * (2.0 * std::f64::consts::PI * f * t + self.phase.to_radians()).sin()
    }

    fn check(&self, name: &str) -> Result<(), Error> {
        match (&self.table, self.amplitude, self.frequency) {
            (Some(tab), None, None) => {
                if tab.is_empty() || tab.windows(2).any(|w| w[1][0] <= w[0][0]) {
                    return Err(bad(format!(
                        "[histories.{name}] table: give points with increasing times"
                    )));
                }
                Ok(())
            }
            (None, Some(_), Some(f)) if f >= 0.0 => Ok(()),
            _ => Err(bad(format!(
                "[histories.{name}]: give either `table` or `amplitude` and `frequency`"
            ))),
        }
    }
}

/// Integration of a transient analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum Method {
    /// Generalised-α on the whole model.
    #[default]
    Direct,
    /// Mode superposition of the lowest `modes`.
    Modal,
}

fn default_rho() -> f64 {
    0.9
}

/// `[transient]`: a response in time from rest.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TransientKeys {
    /// End time (s).
    pub end_time: f64,
    /// Time step (s).
    pub time_step: f64,
    /// `direct` (generalised-α, default) or `modal`.
    #[serde(default)]
    pub method: Method,
    /// Spectral radius at infinite frequency of generalised-α (0 to 1,
    /// default 0.9; 1 adds no numerical damping).
    #[serde(default = "default_rho")]
    pub rho_infinity: f64,
    /// Modes kept by the modal method.
    #[serde(default)]
    pub modes: usize,
    /// Keep every n-th step in the results (default: at most 1000 kept).
    #[serde(default)]
    pub output_every: Option<usize>,
    /// Times (s) whose whole displacement field goes into `structure.vtk`
    /// (`transient_<i>`, the step nearest each).
    #[serde(default)]
    pub snapshots: Vec<f64>,
}

/// `[damping]`.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DampingKeys {
    /// Rayleigh factors [α (1/s), β (s)]: C = αM + βK.
    #[serde(default)]
    pub rayleigh: Option<[f64; 2]>,
    /// Two [frequency (Hz), ratio] pairs that set the Rayleigh factors.
    #[serde(default)]
    pub ratios: Option<[[f64; 2]; 2]>,
    /// Modal method: the same ratio for every mode (added to Rayleigh's).
    #[serde(default)]
    pub modal_ratio: f64,
}

impl DampingKeys {
    fn rayleigh(&self) -> Result<Rayleigh, Error> {
        match (self.rayleigh, self.ratios) {
            (Some(_), Some(_)) => Err(bad("[damping]: give `rayleigh` or `ratios`, not both")),
            (Some([a, b]), None) => {
                if a < 0.0 || b < 0.0 || a.is_nan() || b.is_nan() {
                    return Err(bad("[damping] rayleigh: factors must be ≥ 0"));
                }
                Ok(Rayleigh { alpha: a, beta: b })
            }
            (None, Some([[f1, z1], [f2, z2]])) => Rayleigh::from_ratios(f1, z1, f2, z2),
            (None, None) => Ok(Rayleigh::default()),
        }
    }
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

/// `[[supports]]`: fixed unknowns at a node, along a line of nodes, or on
/// a named surface of a solid (`part` and `surface`).
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Support {
    /// The solid whose `surface` is held.
    #[serde(default)]
    pub part: Option<String>,
    /// A named surface of the solid: every node on it.
    #[serde(default)]
    pub surface: Option<String>,
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
    /// Time history that scales the load (`[histories.<name>]`); none:
    /// constant.
    #[serde(default)]
    pub history: Option<String>,
}

/// `[[line_loads]]`: a uniform force per length on a beam line.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LineLoadKeys {
    /// Name of the beam line.
    pub beam: String,
    /// Force per length, global axes (N/m).
    pub force_per_length: [f64; 3],
    /// Time history that scales the load.
    #[serde(default)]
    pub history: Option<String>,
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
    /// `[[solids]]`.
    #[serde(default)]
    pub solids: Vec<SolidKeys>,
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
    /// `[[temperatures]]`.
    #[serde(default)]
    pub temperatures: Vec<Temperature>,
    /// `[histories]`.
    #[serde(default)]
    pub histories: BTreeMap<String, History>,
    /// `[transient]`.
    #[serde(default)]
    pub transient: Option<TransientKeys>,
    /// `[damping]`.
    #[serde(default)]
    pub damping: DampingKeys,
    /// `[[contacts]]`.
    #[serde(default)]
    pub contacts: Vec<ContactKeys>,
    /// `[contact]`.
    #[serde(default)]
    pub contact: ContactSettings,
    /// `[points]`: named points whose results are reported.
    #[serde(default)]
    pub points: BTreeMap<String, [f64; 3]>,
    /// `[output]`.
    #[serde(default)]
    pub output: OutputKeys,
    /// The folder of the case file: relative mesh paths start there (set
    /// by whoever reads the file; default the working directory).
    #[serde(skip)]
    pub base: Option<std::path::PathBuf>,
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

    /// Reads a case file; relative mesh paths start at its folder.
    pub fn from_file(path: &std::path::Path) -> Result<Self, Error> {
        let text =
            std::fs::read_to_string(path).map_err(|e| bad(format!("{}: {e}", path.display())))?;
        let mut c = Self::from_toml(&text)?;
        c.base = path.parent().map(|p| p.to_path_buf());
        Ok(c)
    }

    /// The mesh files the case reads, resolved.
    pub fn mesh_paths(&self) -> Vec<std::path::PathBuf> {
        self.solids
            .iter()
            .filter_map(|s| s.mesh.as_ref())
            .map(|m| match &self.base {
                Some(b) => b.join(m),
                None => std::path::PathBuf::from(m),
            })
            .collect()
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

/// A solid part of the model.
#[derive(Debug, Clone)]
pub struct SolidPart {
    /// Its elements.
    pub elements: std::ops::Range<usize>,
    /// Named surfaces: faces (model nodes, pointing out of the solid).
    pub surfaces: BTreeMap<String, Vec<Vec<usize>>>,
    /// Where the mesh came from.
    pub source: String,
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
    /// The solid parts, by name.
    pub solids: BTreeMap<String, SolidPart>,
    /// The contact planes with their nodes.
    pub contacts: Vec<crate::contact::Plane>,
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
    if case.beams.is_empty() && case.plates.is_empty() && case.solids.is_empty() {
        return Err(bad(
            "the case has no `[[beams]]`, no `[[plates]]` and no `[[solids]]`",
        ));
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
    // per beam line: α and the section's heights along local y and z
    type LineThermal = (Option<f64>, Option<(f64, f64)>);
    let mut line_thermal: BTreeMap<String, LineThermal> = BTreeMap::new();
    let mut plate_thermal: BTreeMap<String, (Option<f64>, f64)> = BTreeMap::new();
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
        let heights = match sec.shape {
            Shape::Rectangle => sec.width.zip(sec.thickness),
            Shape::Circle => sec.diameter.map(|d| (d, d)),
            Shape::General => None,
        };
        line_thermal.insert(b.name.clone(), (mat.expansion, heights));
    }
    let mut plates = BTreeMap::new();
    for p in &case.plates {
        if lines.contains_key(&p.name) || !seen.insert(p.name.clone()) {
            return Err(bad(format!("two parts are named `{}`", p.name)));
        }
        let mat = case.materials.get(&p.material).ok_or_else(|| {
            bad(format!(
                "plate `{}`: no material `{}` in [materials]",
                p.name, p.material
            ))
        })?;
        if !(mat.young > 0.0 && mat.density > 0.0 && mat.poisson > -1.0 && mat.poisson < 0.5) {
            return Err(bad(format!(
                "plate `{}`: the material needs young > 0, density > 0 and −1 < poisson < 0.5",
                p.name
            )));
        }
        if p.thickness.is_nan() || p.thickness <= 0.0 {
            return Err(bad(format!(
                "plate `{}`: thickness must be positive",
                p.name
            )));
        }
        let [na, nb] = p.divisions;
        if na == 0 || nb == 0 {
            return Err(bad(format!(
                "plate `{}`: `divisions` must be at least 1",
                p.name
            )));
        }
        let n = crate::beam::cross(p.edge_a, p.edge_b);
        if crate::beam::norm(n) <= 1e-12 * crate::beam::norm(p.edge_a) * crate::beam::norm(p.edge_b)
        {
            return Err(bad(format!(
                "plate `{}`: edge_a and edge_b are parallel or zero",
                p.name
            )));
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
        plate_thermal.insert(p.name.clone(), (mat.expansion, p.thickness));
    }
    let mut solids: BTreeMap<String, SolidPart> = BTreeMap::new();
    let mut solid_alpha: BTreeMap<String, Option<f64>> = BTreeMap::new();
    for sk in &case.solids {
        let what = format!("solid `{}`", sk.name);
        if lines.contains_key(&sk.name)
            || plates.contains_key(&sk.name)
            || !seen.insert(sk.name.clone())
        {
            return Err(bad(format!("two parts are named `{}`", sk.name)));
        }
        let mat = case.materials.get(&sk.material).ok_or_else(|| {
            bad(format!(
                "{what}: no material `{}` in [materials]",
                sk.material
            ))
        })?;
        if !(mat.young > 0.0 && mat.density > 0.0 && mat.poisson > -1.0 && mat.poisson < 0.5) {
            return Err(bad(format!(
                "{what}: the material needs young > 0, density > 0 and −1 < poisson < 0.5"
            )));
        }
        let mesh = match (&sk.mesh, &sk.block) {
            (Some(path), None) => {
                let full = match &case.base {
                    Some(b) => b.join(path),
                    None => std::path::PathBuf::from(path),
                };
                let raw = foehncraft_mesh::gmsh::read_raw(&full)
                    .map_err(|e| bad(format!("{what}: {e}")))?;
                crate::solid_mesh::from_gmsh(&raw)
                    .and_then(|m| m.select(&sk.volumes))
                    .map_err(|e| bad(format!("{what}: {}: {e}", full.display())))?
            }
            (None, Some(b)) => {
                if !sk.volumes.is_empty() {
                    return Err(bad(format!("{what}: `volumes` is for a `mesh`")));
                }
                if b.size.iter().any(|v| !(v.is_finite() && *v > 0.0))
                    || b.origin.iter().any(|v| !v.is_finite())
                {
                    return Err(bad(format!(
                        "{what}: the block's `size` must be positive and `origin` finite"
                    )));
                }
                if b.divisions.contains(&0) {
                    return Err(bad(format!(
                        "{what}: `divisions` must be at least 1 along each axis"
                    )));
                }
                if b.divisions.iter().product::<usize>() > 2_000_000 {
                    return Err(bad(format!("{what}: more than 2 million cells")));
                }
                let kind = match b.element {
                    SolidElement::Hex20 => crate::solid::Kind::Hex20,
                    SolidElement::Tet10 => crate::solid::Kind::Tet10,
                };
                crate::solid_mesh::block(b.origin, b.size, b.divisions, kind)
            }
            _ => {
                return Err(bad(format!(
                    "{what}: give `mesh` (a Gmsh file) or `block`, one of them"
                )));
            }
        };
        model.solid_props.push(crate::solid::SolidProps {
            young: mat.young,
            poisson: mat.poisson,
            density: mat.density,
        });
        let pi = model.solid_props.len() - 1;
        let off = model.nodes.len();
        for &x in &mesh.nodes {
            model.add_node(x);
            widths.push(None);
        }
        let first = model.solids.len();
        for (kind, e) in &mesh.elements {
            let v: Vec<usize> = e.iter().map(|n| n + off).collect();
            model.add_solid(*kind, &v, pi);
        }
        for e in first..model.solids.len() {
            crate::solid::volume(model.solids[e].kind, &model.solid_coords(e)).map_err(|m| {
                let c = model.nodes[model.solids[e].nodes[0]];
                bad(format!(
                    "{what}: element {} (first node at ({}, {}, {})): {m}",
                    e - first,
                    c[0],
                    c[1],
                    c[2]
                ))
            })?;
        }
        let surfaces = mesh
            .surfaces
            .into_iter()
            .map(|(k, v)| {
                (
                    k,
                    v.into_iter()
                        .map(|f| f.into_iter().map(|n| n + off).collect())
                        .collect(),
                )
            })
            .collect();
        solids.insert(
            sk.name.clone(),
            SolidPart {
                elements: first..model.solids.len(),
                surfaces,
                source: mesh.source,
            },
        );
        solid_alpha.insert(sk.name.clone(), mat.expansion);
    }
    // the faces of a named surface of a solid
    let surface = |part: &str, name: &str, what: &str| -> Result<&Vec<Vec<usize>>, Error> {
        let p = solids
            .get(part)
            .ok_or_else(|| bad(format!("{what}: no solid `{part}`")))?;
        p.surfaces.get(name).ok_or_else(|| {
            let names: Vec<String> = p.surfaces.keys().map(|k| format!("`{k}`")).collect();
            bad(format!(
                "{what}: solid `{part}` has no surface `{name}`; it has {}",
                if names.is_empty() {
                    "none".to_string()
                } else {
                    names.join(", ")
                }
            ))
        })
    };
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
    let along =
        |model: &Model, seg: [[f64; 3]; 2], what: &str| -> Result<Vec<(usize, f64)>, Error> {
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
                    let lo = if k == 0 {
                        hits[0].0
                    } else {
                        0.5 * (hits[k - 1].0 + hits[k].0)
                    };
                    let hi = if k + 1 == m {
                        hits[m - 1].0
                    } else {
                        0.5 * (hits[k].0 + hits[k + 1].0)
                    };
                    (hits[k].1, (hi - lo) * len)
                })
                .collect())
        };
    let place = |model: &Model,
                 at: Option<[f64; 3]>,
                 seg: Option<[[f64; 3]; 2]>,
                 what: &str|
     -> Result<Vec<(usize, f64)>, Error> {
        match (at, seg) {
            (Some(x), None) => Ok(vec![(find(model, x, what)?, 0.0)]),
            (None, Some(s)) => along(model, s, what),
            _ => Err(bad(format!(
                "{what}: give `at` (a node) or `along` (a segment), not both"
            ))),
        }
    };
    for s in &case.supports {
        let nodes: Vec<usize> = match (&s.part, &s.surface) {
            (None, None) => place(&model, s.at, s.along, "[[supports]]")?
                .into_iter()
                .map(|x| x.0)
                .collect(),
            (Some(p), Some(sf)) => {
                if s.at.is_some() || s.along.is_some() {
                    return Err(bad(
                        "[[supports]]: give `at`, `along` or `part` with `surface`, one of them",
                    ));
                }
                let mut v: Vec<usize> = surface(p, sf, "[[supports]]")?
                    .iter()
                    .flatten()
                    .copied()
                    .collect();
                v.sort_unstable();
                v.dedup();
                v
            }
            _ => {
                return Err(bad("[[supports]]: `part` and `surface` go together"));
            }
        };
        for n in nodes {
            for d in dofs_of(&s.fix) {
                model.fixed[DOFS * n + d] = true;
            }
        }
    }
    let mut clamps = Vec::new();
    for c in &case.clamps {
        let nodes = place(&model, c.at, c.along, "[[clamps]]")?;
        if c.normal == c.line {
            return Err(bad(
                "[[clamps]]: `normal` and `line` must be different axes",
            ));
        }
        let (width, share): (f64, Vec<f64>) = if c.along.is_some() {
            if c.width.is_some() {
                return Err(bad(
                    "[[clamps]]: a clamp `along` an edge takes its length from the segment; leave out `width`",
                ));
            }
            (
                nodes.iter().map(|x| x.1).sum(),
                nodes.iter().map(|x| x.1).collect(),
            )
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
                return Err(bad(format!(
                    "[[clamps]]: spring stiffness must be ≥ 0 (is {k})"
                )));
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
        if !pr.value.is_finite() || pr.force_per_area.iter().any(|v| !v.is_finite()) {
            return Err(bad("[[pressures]]: values must be finite"));
        }
        match (&pr.plate, &pr.part, &pr.surface) {
            (Some(plate), None, None) => {
                let r = plates
                    .get(plate)
                    .ok_or_else(|| bad(format!("[[pressures]]: no plate `{plate}`")))?;
                for e in r.clone() {
                    model.area_loads.push(crate::model::AreaLoad {
                        element: e,
                        q: pr.force_per_area,
                        pressure: pr.value,
                    });
                }
            }
            (None, Some(p), Some(sf)) => {
                for f in surface(p, sf, "[[pressures]]")? {
                    let mut nodes = [0; 8];
                    nodes[..f.len()].copy_from_slice(f);
                    model.face_loads.push(crate::model::FaceLoad {
                        nodes,
                        count: f.len(),
                        q: pr.force_per_area,
                        pressure: pr.value,
                    });
                }
            }
            _ => {
                return Err(bad(
                    "[[pressures]]: give `plate`, or `part` with `surface` (a solid)",
                ));
            }
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
    for t in &case.temperatures {
        let what = format!("[[temperatures]] `{}`", t.part);
        if [t.change, t.difference, t.difference_y]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(bad(format!("{what}: temperatures must be finite")));
        }
        let need = |a: Option<f64>| {
            a.filter(|a| a.is_finite())
                .ok_or_else(|| bad(format!("{what}: its material needs `expansion` (1/K)")))
        };
        if let Some(r) = lines.get(&t.part) {
            let (alpha, h) = line_thermal[&t.part];
            let alpha = need(alpha)?;
            let curvature = if t.difference != 0.0 || t.difference_y != 0.0 {
                let (hy, hz) = h.ok_or_else(|| {
                    bad(format!("{what}: a difference across a `general` section needs its size; use a rectangle or a circle"))
                })?;
                [alpha * t.difference_y / hy, alpha * t.difference / hz]
            } else {
                [0.0, 0.0]
            };
            for e in r.clone() {
                model.thermal.push(crate::model::ThermalLoad {
                    element: e,
                    shell: false,
                    strain: alpha * t.change,
                    curvature,
                });
            }
        } else if let Some(r) = plates.get(&t.part) {
            let (alpha, h) = plate_thermal[&t.part];
            let alpha = need(alpha)?;
            if t.difference_y != 0.0 {
                return Err(bad(format!(
                    "{what}: `difference_y` is for beams; a plate takes `difference`"
                )));
            }
            for e in r.clone() {
                model.thermal.push(crate::model::ThermalLoad {
                    element: e,
                    shell: true,
                    strain: alpha * t.change,
                    curvature: [alpha * t.difference / h, 0.0],
                });
            }
        } else if let Some(sp) = solids.get(&t.part) {
            let alpha = need(solid_alpha[&t.part])?;
            if t.difference != 0.0 || t.difference_y != 0.0 {
                return Err(bad(format!(
                    "{what}: a solid takes a uniform `change`; differences are for beams and plates"
                )));
            }
            for e in sp.elements.clone() {
                model.solid_thermal.push((e, alpha * t.change));
            }
        } else {
            return Err(bad(format!(
                "{what}: no beam line, plate or solid of that name"
            )));
        }
    }
    let mut contacts = Vec::new();
    for ck in &case.contacts {
        let what = format!("[[contacts]] `{}`", ck.name);
        if contacts.iter().any(|p: &crate::contact::Plane| p.name == ck.name) {
            return Err(bad(format!("{what}: two contacts have this name")));
        }
        let len = crate::beam::norm(ck.normal);
        if !(len > 0.0 && len.is_finite()) || ck.point.iter().any(|v| !v.is_finite()) {
            return Err(bad(format!("{what}: `normal` must be a finite non-zero vector and `point` finite")));
        }
        if !(ck.offset.is_finite() && ck.offset >= 0.0) {
            return Err(bad(format!("{what}: `offset` must be ≥ 0")));
        }
        let mut nodes: Vec<usize> = if let Some(sf) = &ck.surface {
            surface(&ck.part, sf, &what)?.iter().flatten().copied().collect()
        } else if let Some(r) = lines.get(&ck.part) {
            r.clone().flat_map(|e| model.beams[e].nodes).collect()
        } else if let Some(r) = plates.get(&ck.part) {
            r.clone()
                .flat_map(|e| model.shells[e].used().to_vec())
                .collect()
        } else if let Some(sp) = solids.get(&ck.part) {
            sp.elements
                .clone()
                .flat_map(|e| model.solids[e].used().to_vec())
                .collect()
        } else {
            return Err(bad(format!(
                "{what}: no beam line, plate or solid `{}`",
                ck.part
            )));
        };
        nodes.sort_unstable();
        nodes.dedup();
        contacts.push(crate::contact::Plane {
            name: ck.name.clone(),
            point: ck.point,
            normal: crate::beam::scale(ck.normal, 1.0 / len),
            offset: ck.offset,
            nodes,
        });
    }
    let mut points = BTreeMap::new();
    for (name, x) in &case.points {
        points.insert(name.clone(), find(&model, *x, &format!("[points] {name}"))?);
    }
    Ok(Built {
        model,
        lines,
        plates,
        solids,
        contacts,
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
    /// The transient result.
    pub transient: Option<Transient>,
    /// The contact state of the static result.
    pub contact: Option<StaticContact>,
}

/// The contact of a static solution.
#[derive(Debug, Clone)]
pub struct StaticContact {
    /// Per plane: total force (N, along its normal) and nodes in contact.
    pub planes: Vec<(f64, usize)>,
    /// Largest penetration (m).
    pub penetration: f64,
    /// Linear solves.
    pub solves: usize,
}

/// The response in time.
#[derive(Debug, Clone)]
pub struct Transient {
    /// Method.
    pub method: Method,
    /// Spectral radius of generalised-α (direct method).
    pub rho_infinity: f64,
    /// Rayleigh damping.
    pub damping: Rayleigh,
    /// Modal method: frequencies (Hz) and damping ratios of the modes kept.
    pub modes: Vec<(f64, f64)>,
    /// Time step (s).
    pub time_step: f64,
    /// Steps taken.
    pub steps: usize,
    /// Times kept (s).
    pub times: Vec<f64>,
    /// Displacements of the named points at the kept times.
    pub points: BTreeMap<String, Vec<[f64; 3]>>,
    /// Kinetic energy ½ vᵀMv at the kept times (J).
    pub kinetic: Vec<f64>,
    /// Strain energy ½ uᵀKu at the kept times (J).
    pub strain: Vec<f64>,
    /// Largest displacement of a node over every step: value, node, time.
    pub max_displacement: (f64, usize, f64),
    /// Displacements of every unknown at the end.
    pub last: Vec<f64>,
    /// The fields at the requested snapshot times: (time, every unknown).
    pub snapshots: Vec<(f64, Vec<f64>)>,
    /// Contact force of each plane (N) at the kept times.
    pub contact_force: Vec<Vec<f64>>,
    /// Largest penetration over every step (m).
    pub penetration: f64,
    /// Factorisations of the step matrix (contact).
    pub factorisations: usize,
}

/// The load vectors (free unknowns) of each history: `None` constant.
/// A load vector and the history that scales it (`None`: constant).
type LoadGroup = (Option<String>, Vec<f64>);

fn load_groups(case: &Case) -> Result<Vec<LoadGroup>, Error> {
    let mut keys: Vec<Option<String>> = vec![None];
    let used = case
        .loads
        .iter()
        .map(|l| &l.history)
        .chain(case.line_loads.iter().map(|l| &l.history))
        .chain(case.pressures.iter().map(|l| &l.history))
        .chain(case.temperatures.iter().map(|l| &l.history));
    for h in used.flatten() {
        if !case.histories.contains_key(h) {
            return Err(bad(format!(
                "no `[histories.{h}]` for a load that names it"
            )));
        }
        if !keys.contains(&Some(h.clone())) {
            keys.push(Some(h.clone()));
        }
    }
    let mut out = Vec::new();
    for key in keys {
        let mut c = case.clone();
        c.loads.retain(|l| l.history == key);
        c.line_loads.retain(|l| l.history == key);
        c.pressures.retain(|l| l.history == key);
        c.temperatures.retain(|l| l.history == key);
        let f = build(&c)?.model.assemble().f;
        out.push((key, f));
    }
    Ok(out)
}

fn transient(case: &Case, tk: &TransientKeys, built: &Built) -> Result<Transient, Error> {
    for (name, h) in &case.histories {
        h.check(name)?;
    }
    if !(tk.end_time > 0.0 && tk.time_step > 0.0 && tk.end_time.is_finite()) {
        return Err(bad(
            "[transient]: `end_time` and `time_step` must be positive",
        ));
    }
    let steps = (tk.end_time / tk.time_step - 1e-9).ceil().max(1.0) as usize;
    if steps > 1_000_000 {
        return Err(bad(format!(
            "[transient]: {steps} steps; at most 10⁶ (raise `time_step`)"
        )));
    }
    let every = tk.output_every.unwrap_or(steps.div_ceil(1000)).max(1);
    let damping = case.damping.rayleigh()?;
    let model = &built.model;
    let sys = model.assemble();
    let groups = load_groups(case)?;
    let nf = sys.free.len();
    let load = |t: f64| -> Vec<f64> {
        let mut f = vec![0.0; nf];
        for (key, g) in &groups {
            let s = key.as_ref().map_or(1.0, |k| case.histories[k].at(t));
            for (a, b) in f.iter_mut().zip(g) {
                *a += s * b;
            }
        }
        f
    };
    let dt = tk.time_step;
    let mut out = Transient {
        method: tk.method,
        rho_infinity: tk.rho_infinity,
        damping,
        modes: Vec::new(),
        time_step: dt,
        steps,
        times: Vec::new(),
        points: built
            .points
            .keys()
            .map(|k| (k.clone(), Vec::new()))
            .collect(),
        kinetic: Vec::new(),
        strain: Vec::new(),
        max_displacement: (0.0, 0, 0.0),
        last: Vec::new(),
        snapshots: Vec::new(),
        contact_force: vec![Vec::new(); built.contacts.len()],
        penetration: 0.0,
        factorisations: 0,
    };
    let snap_steps: Vec<usize> = tk
        .snapshots
        .iter()
        .map(|&t| ((t / dt).round().max(0.0) as usize).min(steps))
        .collect();
    let record = |out: &mut Transient, k: usize, t: f64, u: &[f64], v: &[f64]| {
        let full = sys.expand_zero(u);
        for n in 0..model.nodes.len() {
            let d = crate::beam::norm(node_vec(&full, n, 0));
            if d > out.max_displacement.0 {
                out.max_displacement = (d, n, t);
            }
        }
        if k.is_multiple_of(every) || k == steps {
            out.times.push(t);
            for (name, &n) in &built.points {
                if let Some(h) = out.points.get_mut(name) {
                    h.push(node_vec(&full, n, 0));
                }
            }
            out.kinetic.push(0.5 * sys.m.form(v, v));
            out.strain.push(0.5 * sys.k.form(u, u));
        }
        for &ks in &snap_steps {
            if ks == k {
                out.snapshots.push((t, full.clone()));
            }
        }
        if k == steps {
            out.last = full;
        }
    };
    if !built.contacts.is_empty() {
        if tk.method == Method::Modal {
            return Err(bad(
                "[transient] method = \"modal\" cannot take `[[contacts]]` (contact is not linear); use \"direct\"",
            ));
        }
        let c = crate::contact::Contacts::new(
            model,
            &sys,
            built.contacts.clone(),
            case.contact.penalty,
            case.contact.tolerance * model.size(),
        )?;
        let ga = GenAlpha::new(&sys, damping, dt, tk.rho_infinity)?;
        let mut cache = crate::contact::StepCache::default();
        let mut f0 = load(0.0);
        let mut st = ga.initial(0.0, vec![0.0; nf], vec![0.0; nf], &f0)?;
        record(&mut out, 0, 0.0, &st.u, &st.v);
        for f in &mut out.contact_force {
            f.push(0.0);
        }
        for k in 1..=steps {
            let f1 = load(k as f64 * dt);
            let (next, cs) = ga.step_contact(model, &c, &mut cache, &st, &f0, &f1)?;
            st = next;
            let kept = out.times.len();
            record(&mut out, k, k as f64 * dt, &st.u, &st.v);
            out.penetration = out.penetration.max(cs.penetration());
            if out.times.len() > kept {
                for (j, (f, _)) in cs.per_plane(&c).into_iter().enumerate() {
                    out.contact_force[j].push(f);
                }
            }
            f0 = f1;
        }
        out.factorisations = cache.factorisations;
        return Ok(out);
    }
    match tk.method {
        Method::Direct => {
            let ga = GenAlpha::new(&sys, damping, dt, tk.rho_infinity)?;
            let mut f0 = load(0.0);
            let mut st = ga.initial(0.0, vec![0.0; nf], vec![0.0; nf], &f0)?;
            record(&mut out, 0, 0.0, &st.u, &st.v);
            for k in 1..=steps {
                let f1 = load(k as f64 * dt);
                st = ga.step(&st, &f0, &f1);
                record(&mut out, k, k as f64 * dt, &st.u, &st.v);
                f0 = f1;
            }
        }
        Method::Modal => {
            if tk.modes == 0 {
                return Err(bad("[transient] method = \"modal\" needs `modes`"));
            }
            let basis = ModalBasis::lowest(model, tk.modes, case.damping.modal_ratio, damping)?;
            out.modes = basis
                .omega
                .iter()
                .zip(&basis.zeta)
                .map(|(w, z)| (w / (2.0 * std::f64::consts::PI), *z))
                .collect();
            let nm = basis.omega.len();
            let (mut q, mut qd) = (vec![0.0; nm], vec![0.0; nm]);
            let mut p0 = basis.project(&load(0.0));
            record(&mut out, 0, 0.0, &vec![0.0; nf], &vec![0.0; nf]);
            for k in 1..=steps {
                let p1 = basis.project(&load(k as f64 * dt));
                basis.step(&mut q, &mut qd, &p0, &p1, dt);
                record(
                    &mut out,
                    k,
                    k as f64 * dt,
                    &basis.expand(&q),
                    &basis.expand(&qd),
                );
                p0 = p1;
            }
        }
    }
    Ok(out)
}

/// Runs the analyses of a case.
pub fn run(case: &Case) -> Result<Outcome, Error> {
    let built = build(case)?;
    if !case.analysis.statics && case.analysis.modes == 0 && case.transient.is_none() {
        return Err(bad(
            "[analysis] asks for nothing: set `static = true` or `modes`, or add `[transient]`",
        ));
    }
    let (statics, contact) = if !case.analysis.statics {
        (None, None)
    } else if built.contacts.is_empty() {
        (Some(analysis::solve_static(&built.model)?), None)
    } else {
        let model = &built.model;
        let sys = model.assemble();
        let tol = case.contact.tolerance * model.size();
        let c = crate::contact::Contacts::new(
            model,
            &sys,
            built.contacts.clone(),
            case.contact.penalty,
            tol,
        )?;
        let (uf, cs, solves) = crate::contact::solve_static(model, &sys, &c, tol)?;
        let u = sys.expand(&uf);
        let ku = sys.k_all.apply(&u);
        let fc = c.force_vector(&cs, model.dofs());
        // reactions of the supports and springs; the contact forces are
        // reported per plane
        let mut reactions: Vec<f64> = (0..u.len()).map(|i| ku[i] - sys.f_all[i] - fc[i]).collect();
        for &(d, k) in &model.springs {
            reactions[d] -= k * u[d];
        }
        let energy = 0.5 * crate::sparse::dot(&u, &ku);
        (
            Some(Static {
                u,
                reactions,
                energy,
                fill: 0,
            }),
            Some(StaticContact {
                planes: cs.per_plane(&c),
                penetration: cs.penetration(),
                solves,
            }),
        )
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
    let transient = match &case.transient {
        Some(tk) => Some(transient(case, tk, &built)?),
        None => None,
    };
    Ok(Outcome {
        built,
        statics,
        modal,
        transient,
        contact,
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
    /// The von Mises stress at each node of the solids for the
    /// displacements `u` (every unknown), with the static thermal strains:
    /// the mean over the elements around the node (zero off the solids).
    pub fn solid_von_mises(&self, u: &[f64]) -> Vec<f64> {
        let m = &self.built.model;
        let mut sum = vec![0.0; m.nodes.len()];
        let mut count = vec![0usize; m.nodes.len()];
        let mut strain = vec![0.0; m.solids.len()];
        for &(e, a) in &m.solid_thermal {
            strain[e] += a;
        }
        for (e, el) in m.solids.iter().enumerate() {
            let ue: Vec<f64> = el
                .used()
                .iter()
                .flat_map(|&n| (0..3).map(move |d| u[DOFS * n + d]))
                .collect();
            let Ok(st) = crate::solid::nodal_stresses(
                el.kind,
                &m.solid_coords(e),
                &m.solid_props[el.props],
                &ue,
                strain[e],
            ) else {
                continue;
            };
            for (a, &n) in el.used().iter().enumerate() {
                sum[n] += crate::solid::von_mises(&st[a]);
                count[n] += 1;
            }
        }
        sum.iter()
            .zip(&count)
            .map(|(s, &c)| if c > 0 { s / c as f64 } else { 0.0 })
            .collect()
    }

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
            "elements": m.beams.len() + m.shells.len() + m.solids.len(),
            "beam_elements": m.beams.len(),
            "shell_elements": m.shells.len(),
            "solid_elements": m.solids.len(),
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
        if !self.built.solids.is_empty() {
            out["solids"] = Value::Array(
                self.built
                    .solids
                    .iter()
                    .map(|(k, p)| {
                        json!({
                            "name": k,
                            "elements": p.elements.len(),
                            "element": match m.solids[p.elements.start].kind {
                                crate::solid::Kind::Tet10 => "tet10",
                                crate::solid::Kind::Hex20 => "hex20",
                            },
                            "mesh": p.source,
                            "surfaces": p.surfaces.iter().map(|(n, f)| json!({"name": n, "faces": f.len()})).collect::<Vec<_>>(),
                        })
                    })
                    .collect(),
            );
        }
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
            if let Some(c) = &self.contact {
                out["static"]["contact"] = json!({
                    "planes": self.built.contacts.iter().zip(&c.planes).map(|(p, (f, n))| json!({
                        "name": p.name, "force": f, "nodes_in_contact": n, "nodes": p.nodes.len(),
                    })).collect::<Vec<_>>(),
                    "penetration": c.penetration,
                    "solves": c.solves,
                });
            }
            if !m.solids.is_empty() {
                let vm = self.solid_von_mises(&s.u);
                let (mut best, mut at) = (0.0, 0);
                for (n, &v) in vm.iter().enumerate() {
                    if v > best {
                        best = v;
                        at = n;
                    }
                }
                out["static"]["max_von_mises"] =
                    json!({"value": best, "node": at, "at": m.nodes[at]});
            }
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
        if let Some(tr) = &self.transient {
            let (v, n, t) = tr.max_displacement;
            let points: serde_json::Map<String, Value> = tr
                .points
                .iter()
                .map(|(k, h)| {
                    (
                        k.clone(),
                        json!({"node": self.built.points[k], "displacement": h}),
                    )
                })
                .collect();
            let mut j = json!({
                "method": match tr.method { Method::Direct => "direct", Method::Modal => "modal" },
                "time_step": tr.time_step,
                "steps": tr.steps,
                "damping": tr.damping,
                "times": tr.times,
                "kinetic_energy": tr.kinetic,
                "strain_energy": tr.strain,
                "max_displacement": {"value": v, "node": n, "at": m.nodes[n], "time": t},
                "points": points,
                "snapshots": tr.snapshots.iter().map(|(t, _)| *t).collect::<Vec<_>>(),
            });
            if !self.built.contacts.is_empty() {
                j["contact"] = json!({
                    "planes": self.built.contacts.iter().zip(&tr.contact_force).map(|(p, f)| json!({
                        "name": p.name, "force": f,
                    })).collect::<Vec<_>>(),
                    "penetration": tr.penetration,
                    "factorisations": tr.factorisations,
                });
            }
            match tr.method {
                Method::Direct => j["rho_infinity"] = json!(tr.rho_infinity),
                Method::Modal => {
                    j["modes"] = Value::Array(
                        tr.modes
                            .iter()
                            .map(|(f, z)| json!({"frequency": f, "damping_ratio": z}))
                            .collect(),
                    )
                }
            }
            out["transient"] = j;
        }
        out["warnings"] = Value::Array(warnings);
        out
    }

    /// The transient history as CSV: time, kinetic and strain energy, and
    /// the displacement of each named point. `None` without `[transient]`.
    pub fn to_csv(&self) -> Option<String> {
        use std::fmt::Write;
        let tr = self.transient.as_ref()?;
        let mut s = String::from("time,kinetic_energy,strain_energy");
        for k in tr.points.keys() {
            let _ = write!(s, ",{k}_ux,{k}_uy,{k}_uz");
        }
        for p in &self.built.contacts {
            let _ = write!(s, ",contact_{}", p.name);
        }
        s.push('\n');
        for (i, t) in tr.times.iter().enumerate() {
            let _ = write!(s, "{t:.9e},{:.9e},{:.9e}", tr.kinetic[i], tr.strain[i]);
            for h in tr.points.values() {
                let u = h[i];
                let _ = write!(s, ",{:.9e},{:.9e},{:.9e}", u[0], u[1], u[2]);
            }
            for f in &tr.contact_force {
                let _ = write!(s, ",{:.9e}", f.get(i).copied().unwrap_or(0.0));
            }
            s.push('\n');
        }
        Some(s)
    }

    /// The mesh with the static displacement and the mode shapes, as a
    /// legacy VTK file (ParaView, VisIt): polydata, or an unstructured grid
    /// when the model has solids. Mode shapes are scaled
    /// to a largest translation of 1 (`mode_<n>`, and `mode_<n>_rotation`
    /// with the same factor).
    pub fn to_vtk(&self, name: &str) -> String {
        use std::fmt::Write;
        let m = &self.built.model;
        let mut s = String::new();
        // with solids an unstructured grid (quadratic cells), else polydata
        let grid = !m.solids.is_empty();
        let _ = writeln!(
            s,
            "# vtk DataFile Version 3.0\nFoehncraft structure {name}\nASCII\nDATASET {}",
            if grid {
                "UNSTRUCTURED_GRID"
            } else {
                "POLYDATA"
            }
        );
        let _ = writeln!(s, "POINTS {} double", m.nodes.len());
        for p in &m.nodes {
            let _ = writeln!(s, "{:.9e} {:.9e} {:.9e}", p[0], p[1], p[2]);
        }
        if grid {
            let mut cells: Vec<(u8, Vec<usize>)> = Vec::new();
            for b in &m.beams {
                cells.push((3, b.nodes.to_vec()));
            }
            for e in &m.shells {
                let n = e.used().to_vec();
                cells.push((if n.len() == 3 { 5 } else { 9 }, n));
            }
            for e in &m.solids {
                cells.push((e.kind.vtk_type(), e.used().to_vec()));
            }
            let size: usize = cells.iter().map(|c| c.1.len() + 1).sum();
            let _ = writeln!(s, "CELLS {} {size}", cells.len());
            for (_, c) in &cells {
                let ids: Vec<String> = c.iter().map(|n| n.to_string()).collect();
                let _ = writeln!(s, "{} {}", c.len(), ids.join(" "));
            }
            let _ = writeln!(s, "CELL_TYPES {}", cells.len());
            for (t, _) in &cells {
                let _ = writeln!(s, "{t}");
            }
        }
        if !grid && !m.beams.is_empty() {
            let _ = writeln!(s, "LINES {} {}", m.beams.len(), 3 * m.beams.len());
            for b in &m.beams {
                let _ = writeln!(s, "2 {} {}", b.nodes[0], b.nodes[1]);
            }
        }
        if !grid && !m.shells.is_empty() {
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
            if grid {
                let _ = writeln!(s, "SCALARS von_mises double 1\nLOOKUP_TABLE default");
                for v in self.solid_von_mises(&st.u) {
                    let _ = writeln!(s, "{v:.9e}");
                }
            }
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
        if let Some(tr) = &self.transient {
            vector(&mut s, "transient_end", &tr.last, 0, 1.0);
            for (i, (_, u)) in tr.snapshots.iter().enumerate() {
                vector(&mut s, &format!("transient_{}", i + 1), u, 0, 1.0);
            }
        }
        s
    }
}
