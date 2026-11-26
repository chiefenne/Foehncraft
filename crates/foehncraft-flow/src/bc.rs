//! Boundary conditions of the incompressible flow solver (D29) and the check
//! of their combination before a run (D30).
//!
//! Every boundary of the mesh takes one condition, `[boundary.<boundary>]` with a
//! `type`:
//!
//! | type | velocity | pressure | mass flux |
//! |---|---|---|---|
//! | `wall` | `velocity` (default 0; a moving wall keeps only its tangential part) | extrapolated | 0 |
//! | `slip-wall` | no shear stress, zero normal component (as `symmetry`, on any surface) | extrapolated | 0 |
//! | `velocity-inlet` | `velocity`, `normal_velocity`, `profile = "analytic"`, `profile = "parabolic"` (with `max_velocity`, `profile_axes`) or `profile = "table"` (with `table`, `table_coordinate`) | extrapolated | ρ u_b·S |
//! | `mass-flow-inlet` | normal, from `mass_flow` (kg/s) spread by area | extrapolated | prescribed |
//! | `static-pressure-inlet` | normal to the face, from the flux | `pressure` | Rhie–Chow |
//! | `total-pressure-inlet` | normal to the face, from the flux | p₀ − ½ρu_n² from `total_pressure` | Rhie–Chow |
//! | `static-pressure-outlet` | extrapolated; backflow by `backflow` | `pressure`; backflow by `backflow` | Rhie–Chow |
//! | `total-pressure-outlet` | extrapolated; backflow by `backflow` | p₀ − ½ρu_n² from `total_pressure` | Rhie–Chow |
//! | `convective-outlet` | ∂u/∂t + U_c ∂u/∂n = 0 (transient runs; `convection_velocity` U_c) | extrapolated | ρ u_b·S, balanced with the inflow |
//! | `symmetry` | no shear stress, zero normal component (implicit per component) | extrapolated | 0 |
//! | `wedge` | as `symmetry`: the two side planes of an axisymmetric wedge mesh | extrapolated | 0 |
//! | `periodic` | across the pair (`partner`): the partner's cell is the neighbour, translated by the period | periodic, with a mean gradient as body force | through the pair |
//!
//! **Periodic pairs** (D29, D109 (1), stage 6 part 6p): two boundaries of
//! `type = "periodic"` that name each other as `partner` and carry
//! conformal, translated copies of the same faces (the translation from the
//! face centroids, or `translation`). The flow through the pair is driven
//! by `mass_flow` (kg/s, held by adjusting a uniform body force) or by
//! `pressure_gradient` (Pa/m, the mean pressure falling from this boundary
//! towards its partner), given on one of the two; the pressure itself is
//! periodic and its level fixed in the reference cell. Scalars cross the
//! pair unchanged or, with `periodic = "heat-flux"` in their boundary
//! sub-table, with the jump of the heat-flux mode of D109.
//!
//! **Backflow at pressure outlets** (`backflow`): with `"total-pressure"`
//! (the default) a face where the flow enters takes the prescribed pressure
//! as the total pressure of the entering fluid and lets it enter normal to
//! the face; `"zero-gradient"` keeps the static pressure and the
//! extrapolated velocity; `"suppress"` closes such a face for the outer
//! iteration (no flux; the face acts as a slip wall, with the extrapolated
//! pressure), so no fluid enters. A closed face opens again when the
//! prescribed pressure would drive the flow out of it.
//!
//! Static and total pressure conditions combine freely, for example a static
//! pressure inlet with a total pressure outlet. Turbulence values at inlets
//! and, for the backflow, at pressure outlets (`k` with `omega` or
//! `epsilon`, `turbulence_intensity`, `turbulence_length_scale`,
//! `viscosity_ratio`, `zeta`) are accepted for the turbulence models of
//! stage 6 and reported as unused in laminar flow. At a pressure outlet
//! they act where the flow enters, face by face and outer iteration by
//! outer iteration (fixed value there, zero gradient where it leaves);
//! without them the backflow takes the values of the cell.

use crate::analytic::Analytic;
use crate::diag;
use crate::expr::{Env, Scalar, Value, Vars, Vector};
use anyhow::Result;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_kernels::flow::{
    FLUX_FIXED, FLUX_PRESSURE, FLUX_PRESSURE_NO_BACKFLOW, MOM_FIXED, MOM_OUTFLOW, MOM_SYMMETRY,
    PB_EXTRAPOLATE, PB_FIXED, PB_NO_BACKFLOW, PB_OUTLET, PB_TOTAL, PB_TOTAL_NO_BACKFLOW,
    UB_CONVECTIVE, UB_EXTRAPOLATE, UB_FIXED, UB_FROM_FLUX, UB_NO_BACKFLOW, UB_OUTLET, UB_SYMMETRY,
};
use foehncraft_mesh::geometry::{add, cross, dot, norm, scale, sub};
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

/// The density the boundary conditions take: one number, or one per face
/// of the mesh (D194: several fluids in one case).
#[derive(Debug, Clone, Copy)]
pub enum Density<'a> {
    /// One density everywhere.
    One(f64),
    /// The density of every face.
    Faces(&'a [f64]),
}

impl Density<'_> {
    /// The density at face `f`.
    #[inline]
    pub fn at(&self, f: usize) -> f64 {
        match self {
            Density::One(r) => *r,
            Density::Faces(v) => v[f],
        }
    }
}

impl From<f64> for Density<'_> {
    fn from(r: f64) -> Self {
        Density::One(r)
    }
}

impl<'a> From<&'a crate::simplec::Physics> for Density<'a> {
    fn from(p: &'a crate::simplec::Physics) -> Self {
        match &p.rho_var {
            Some(r) => Density::Faces(&r.1),
            None => Density::One(p.rho),
        }
    }
}
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Inlet profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum InletProfile {
    /// The velocity of the `[analytic]` solution at the face centroid (at
    /// the time of the step in a transient case).
    Analytic,
    /// A parabola normal to the face, u_n = U_max Π_a 4 η_a (1 − η_a) with
    /// η_a = (x_a − x_a,min)/(x_a,max − x_a,min) over the boundary extent along
    /// each axis a of `profile_axes` (developed laminar inflow; two axes for
    /// a duct).
    Parabolic,
    /// A table of the velocity over one coordinate of the face centroid
    /// (`table`, `table_coordinate`), interpolated linearly (measured or
    /// computed profiles).
    Table,
}

/// The coordinate of a table profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum TableCoordinate {
    /// x of the face centroid.
    X,
    /// y of the face centroid.
    Y,
    /// z of the face centroid.
    Z,
    /// The distance of the face centroid from the axis through
    /// `table_origin` along `table_direction` (pipes, annuli, wedges).
    Radius,
}

/// What [`BoundaryCondition::Wall`]-free outlets do where the flow enters
/// (D29, backflow treatment).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Backflow {
    /// The prescribed pressure is the total pressure of the entering fluid,
    /// which enters normal to the face.
    #[default]
    TotalPressure,
    /// Static pressure kept, velocity extrapolated (zero normal gradient).
    ZeroGradient,
    /// The face is closed while the flow would enter.
    Suppress,
}

/// A velocity profile read from a table file (D29): rows of a coordinate
/// and either the normal speed into the domain (two columns) or the three
/// velocity components (four columns), the coordinate increasing. A
/// header that names the columns may add the turbulence of the profile
/// (stage 6 specification 7.1): `k` with `epsilon` or `omega`, and `zeta`
/// ([`ProfileColumns`]).
#[derive(Debug, Clone, PartialEq)]
pub struct ProfileTable {
    /// The file it was read from.
    pub path: PathBuf,
    /// SHA-256 of the file (provenance, D85).
    pub sha256: String,
    /// Coordinate values, strictly increasing.
    pub s: Vec<f64>,
    /// Values per row (the columns after the coordinate).
    pub values: Vec<Vec<f64>>,
    /// The names of the value columns from the header (trimmed, lower
    /// case); empty without a header.
    pub names: Vec<String>,
}

/// The velocity columns of a profile table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VelocityColumns {
    /// The speed into the domain, normal to the face.
    Normal(usize),
    /// The three components.
    Components([usize; 3]),
}

/// What the value columns of a profile table hold (indices into a row of
/// [`ProfileTable::values`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProfileColumns {
    /// The velocity.
    pub velocity: VelocityColumns,
    /// k (m²/s²).
    pub k: Option<usize>,
    /// ε (m²/s³).
    pub epsilon: Option<usize>,
    /// ω (1/s).
    pub omega: Option<usize>,
    /// ζ = v²/k of k-ε-ζ-f and φ of the elliptic-blending k-ε.
    pub zeta: Option<usize>,
}

impl ProfileColumns {
    /// True when the table gives k with ε or ω.
    pub fn has_turbulence(&self) -> bool {
        self.k.is_some()
    }
}

/// A loaded table profile and where its coordinate is measured.
#[derive(Debug, Clone)]
pub struct TableProfile {
    /// The table.
    pub table: Arc<ProfileTable>,
    /// Its coordinate.
    pub coordinate: TableCoordinate,
    /// A point on the axis of a radius.
    pub origin: [f64; 3],
    /// The axis of a radius.
    pub direction: Option<[f64; 3]>,
}

impl TableProfile {
    /// The values of the table row at point `x` (None without the axis of
    /// a radius).
    pub fn at(&self, x: [f64; 3]) -> Option<Vec<f64>> {
        table_position(x, self.coordinate, self.origin, self.direction).map(|s| self.table.at(s).0)
    }
}

/// The column names of a profile table.
pub const PROFILE_COLUMNS: [&str; 8] =
    ["u_n", "u_x", "u_y", "u_z", "k", "epsilon", "omega", "zeta"];

impl ProfileTable {
    /// Reads and parses a table file.
    pub fn read(path: &Path) -> std::result::Result<Self, String> {
        let bytes = std::fs::read(path).map_err(|e| format!("reading {}: {e}", path.display()))?;
        let text = String::from_utf8_lossy(&bytes);
        let mut t = Self::parse(&text).map_err(|e| format!("{}: {e}", path.display()))?;
        t.path = path.to_path_buf();
        t.sha256 = crate::engine::sha256(&bytes);
        Ok(t)
    }

    /// Parses the text of a table: `#` starts a comment; columns are
    /// separated by commas, semicolons or white space; a first line that is
    /// not numeric is a header and skipped.
    pub fn parse(text: &str) -> std::result::Result<Self, String> {
        let mut s = Vec::new();
        let mut values = Vec::new();
        let mut width = None;
        let mut names = Vec::new();
        for (i, raw) in text.lines().enumerate() {
            let line = raw.split('#').next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }
            let cols: Vec<&str> = line
                .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
                .filter(|c| !c.is_empty())
                .collect();
            let nums: std::result::Result<Vec<f64>, _> =
                cols.iter().map(|c| c.parse::<f64>()).collect();
            let Ok(nums) = nums else {
                if s.is_empty() && width.is_none() {
                    // the header
                    width = Some(cols.len());
                    names = cols[1..].iter().map(|c| c.to_lowercase()).collect();
                    continue;
                }
                return Err(format!("line {}: not a number in `{line}`", i + 1));
            };
            // a header may name the columns of other widths (the
            // turbulence of a profile, [`Self::columns`]); without one a
            // table has 2 or 4
            let named = width.is_some();
            if nums.len() < 2 || (!named && nums.len() != 2 && nums.len() != 4) {
                return Err(format!(
                    "line {}: {} columns; a table has 2 (coordinate, normal speed) or 4 \
                     (coordinate, u_x, u_y, u_z), or a header that names its columns",
                    i + 1,
                    nums.len()
                ));
            }
            if width.is_some_and(|w| w != nums.len())
                || values
                    .first()
                    .is_some_and(|v: &Vec<f64>| v.len() + 1 != nums.len())
            {
                return Err(format!(
                    "line {}: {} columns, unlike the lines before",
                    i + 1,
                    nums.len()
                ));
            }
            if nums.iter().any(|x| !x.is_finite()) {
                return Err(format!("line {}: a value that is not finite", i + 1));
            }
            if let Some(&last) = s.last()
                && nums[0] <= last
            {
                return Err(format!(
                    "line {}: the coordinate {} does not increase (after {last})",
                    i + 1,
                    nums[0]
                ));
            }
            s.push(nums[0]);
            values.push(nums[1..].to_vec());
        }
        if s.len() < 2 {
            return Err(format!("{} rows; a table needs at least 2", s.len()));
        }
        Ok(Self {
            path: PathBuf::new(),
            sha256: String::new(),
            s,
            values,
            names,
        })
    }

    /// True for a table of the normal speed (two columns).
    pub fn is_normal(&self) -> bool {
        self.values[0].len() == 1
    }

    /// What the columns of a profile table hold. A header that uses the
    /// names of [`PROFILE_COLUMNS`] names every value column: `u_n` or
    /// `u_x`, `u_y`, `u_z`, optionally `k` with `epsilon` or `omega`, and
    /// `zeta`, in any order. Without such a header two columns are the
    /// coordinate and the normal speed, four the coordinate and the
    /// components.
    pub fn columns(&self) -> std::result::Result<ProfileColumns, String> {
        let width = self.values[0].len();
        let named = self
            .names
            .iter()
            .any(|n| PROFILE_COLUMNS.contains(&n.as_str()));
        if !named {
            return match width {
                1 => Ok(ProfileColumns {
                    velocity: VelocityColumns::Normal(0),
                    k: None,
                    epsilon: None,
                    omega: None,
                    zeta: None,
                }),
                3 => Ok(ProfileColumns {
                    velocity: VelocityColumns::Components([0, 1, 2]),
                    k: None,
                    epsilon: None,
                    omega: None,
                    zeta: None,
                }),
                n => Err(format!(
                    "{} value columns without a header that names them ({})",
                    n,
                    PROFILE_COLUMNS.join(", ")
                )),
            };
        }
        let mut at: BTreeMap<&str, usize> = BTreeMap::new();
        for (i, n) in self.names.iter().enumerate() {
            let Some(&known) = PROFILE_COLUMNS.iter().find(|&&c| c == n) else {
                return Err(format!(
                    "column `{n}` of the header is not one of {}",
                    PROFILE_COLUMNS.join(", ")
                ));
            };
            if at.insert(known, i).is_some() {
                return Err(format!("column `{n}` twice in the header"));
            }
        }
        let comps = [at.get("u_x"), at.get("u_y"), at.get("u_z")];
        let velocity = match (at.get("u_n"), comps) {
            (Some(&i), [None, None, None]) => VelocityColumns::Normal(i),
            (None, [Some(&x), Some(&y), Some(&z)]) => VelocityColumns::Components([x, y, z]),
            _ => {
                return Err(
                    "the header names the velocity as `u_n` or as `u_x`, `u_y`, `u_z`, one of the two"
                        .into(),
                );
            }
        };
        let c = ProfileColumns {
            velocity,
            k: at.get("k").copied(),
            epsilon: at.get("epsilon").copied(),
            omega: at.get("omega").copied(),
            zeta: at.get("zeta").copied(),
        };
        if c.epsilon.is_some() && c.omega.is_some() {
            return Err("columns `epsilon` and `omega` together: give one (with `k`)".into());
        }
        if c.k.is_some() != (c.epsilon.is_some() || c.omega.is_some()) {
            return Err("column `k` goes with `epsilon` or `omega`".into());
        }
        if c.zeta.is_some() && c.k.is_none() {
            return Err("column `zeta` goes with `k` and `epsilon` or `omega`".into());
        }
        Ok(c)
    }

    /// The velocity columns ([`Self::columns`]; the legacy forms when the
    /// header does not name them).
    pub fn velocity_columns(&self) -> VelocityColumns {
        self.columns().map_or(
            if self.is_normal() {
                VelocityColumns::Normal(0)
            } else {
                VelocityColumns::Components([0, 1, 2])
            },
            |c| c.velocity,
        )
    }

    /// The values at coordinate `x`, linear between the rows and constant
    /// beyond the ends; the flag is true when `x` is outside the rows.
    pub fn at(&self, x: f64) -> (Vec<f64>, bool) {
        let n = self.s.len();
        if x <= self.s[0] {
            return (self.values[0].clone(), x < self.s[0]);
        }
        if x >= self.s[n - 1] {
            return (self.values[n - 1].clone(), x > self.s[n - 1]);
        }
        let i = self.s.partition_point(|&v| v <= x) - 1;
        let t = (x - self.s[i]) / (self.s[i + 1] - self.s[i]);
        let v = self.values[i]
            .iter()
            .zip(&self.values[i + 1])
            .map(|(a, b)| a + t * (b - a))
            .collect();
        (v, false)
    }
}

/// A coordinate axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum Axis {
    /// x.
    X,
    /// y.
    Y,
    /// z.
    Z,
}

/// Bounding box of the nodes of a boundary.
pub type Extent = ([f64; 3], [f64; 3]);

/// The bounding box of the nodes of boundary `p`.
/// Stores the extent of every boundary with a parabolic inlet profile on
/// the whole mesh `mesh` in its condition, for a run on ranks (MPI step
/// W7): the profile then spans the whole boundary on every rank, as on one
/// rank, and not the part of the boundary that a rank holds.
pub fn global_profile_extents(bcs: &mut BTreeMap<String, BoundaryCondition>, mesh: &Mesh) {
    for p in &mesh.boundaries {
        if let Some(BoundaryCondition::VelocityInlet {
            profile: Some(InletProfile::Parabolic),
            profile_extent,
            ..
        }) = bcs.get_mut(&p.name)
        {
            *profile_extent = Some(boundary_extent(mesh, p));
        }
    }
}

/// The bounding box of the nodes of boundary `p`.
pub fn boundary_extent(mesh: &Mesh, p: &foehncraft_mesh::Boundary) -> Extent {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for f in p.range() {
        for &n in mesh.face(f) {
            let x = mesh.points[n as usize];
            for i in 0..3 {
                lo[i] = lo[i].min(x[i]);
                hi[i] = hi[i].max(x[i]);
            }
        }
    }
    (lo, hi)
}

/// Turbulence values at an inlet (used from stage 6 on).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct InletTurbulence {
    /// Turbulence intensity (fraction, e.g. 0.05).
    pub turbulence_intensity: Option<f64>,
    /// Turbulence length scale (m).
    pub turbulence_length_scale: Option<f64>,
    /// Ratio of the eddy viscosity to the molecular one.
    pub viscosity_ratio: Option<f64>,
    /// k directly.
    pub k: Option<Scalar>,
    /// ω directly.
    pub omega: Option<Scalar>,
    /// ε directly.
    pub epsilon: Option<Scalar>,
    /// ζ = v²/k (k-ε-ζ-f).
    pub zeta: Option<f64>,
}

impl InletTurbulence {
    /// The keys that are given.
    pub fn given(&self) -> Vec<&'static str> {
        let mut v = Vec::new();
        if self.turbulence_intensity.is_some() {
            v.push("turbulence_intensity");
        }
        if self.turbulence_length_scale.is_some() {
            v.push("turbulence_length_scale");
        }
        if self.viscosity_ratio.is_some() {
            v.push("viscosity_ratio");
        }
        if self.k.is_some() {
            v.push("k");
        }
        if self.omega.is_some() {
            v.push("omega");
        }
        if self.epsilon.is_some() {
            v.push("epsilon");
        }
        if self.zeta.is_some() {
            v.push("zeta");
        }
        v
    }
}

/// One boundary condition, `[boundary.<boundary>]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(tag = "type", rename_all = "kebab-case", deny_unknown_fields)]
#[allow(clippy::large_enum_variant)] // case-file data, a few per case
pub enum BoundaryCondition {
    /// No-slip wall, at rest or moving tangentially.
    Wall {
        /// Wall velocity (m/s): numbers, or expressions of x, y, z and, in
        /// a transient case, t (D53); only the tangential part is used
        /// (normal to itself, a wall of a turning zone, `[[mrf]]` or
        /// `[[sliding]]`, moves with the zone).
        #[serde(default = "zero_vector")]
        velocity: Vector,
        /// Mass flux through the wall into the domain (kg/(m² s), D124 W3):
        /// suction (< 0) or blowing (> 0), a number or an expression; a
        /// coupled module (the wall film of D123) sets it per face.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        mass_flux: Option<Scalar>,
        /// Equivalent sand-grain roughness k_s (m) of the turbulence models
        /// (stage 6, D134; SST: the wall values of Aupoix 2015; the k-ε
        /// family: the log-law shift of Cebeci & Bradshaw 1977).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        roughness: Option<f64>,
        /// C_s of the roughness shift of Cebeci & Bradshaw (k-ε family;
        /// default 0.5).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        roughness_constant: Option<f64>,
    },
    /// Prescribed velocity.
    VelocityInlet {
        /// Velocity vector (m/s): numbers or expressions (D53).
        #[serde(default)]
        velocity: Option<Vector>,
        /// Speed normal to the face, into the domain (m/s): a number or an
        /// expression (D53).
        #[serde(default)]
        normal_velocity: Option<Scalar>,
        /// Profile instead of a constant.
        #[serde(default)]
        profile: Option<InletProfile>,
        /// Largest speed of the parabolic profile (m/s).
        #[serde(default)]
        max_velocity: Option<f64>,
        /// Axes across which the parabolic profile varies.
        #[serde(default)]
        profile_axes: Vec<Axis>,
        /// The table of a table profile (a file relative to the case file:
        /// coordinate and normal speed, or coordinate and u_x, u_y, u_z; a
        /// header that names the columns `u_n` or `u_x`, `u_y`, `u_z` may
        /// add the turbulence `k` with `epsilon` or `omega`, and `zeta`).
        #[serde(default)]
        table: Option<PathBuf>,
        /// The coordinate of the table.
        #[serde(default)]
        table_coordinate: Option<TableCoordinate>,
        /// A point on the axis of `table_coordinate = "radius"` (m).
        #[serde(default)]
        table_origin: [f64; 3],
        /// The direction of the axis of `table_coordinate = "radius"`.
        #[serde(default)]
        table_direction: Option<[f64; 3]>,
        /// The table, read with the case.
        #[serde(skip)]
        #[schemars(skip)]
        table_data: Option<Arc<ProfileTable>>,
        /// The extent of the boundary on the whole mesh, for the parabolic
        /// profile of a run on ranks (MPI step W7), where a rank holds only
        /// a part of the boundary ([`global_profile_extents`]).
        #[serde(skip)]
        #[schemars(skip)]
        profile_extent: Option<Extent>,
        /// Turbulence values (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
        /// k (m²/s²) directly: a number or an expression (D53) (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        k: Option<Scalar>,
        /// ω (1/s) directly, with `k` (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        omega: Option<Scalar>,
        /// ε (m²/s³) directly, with `k` (the k-ε family; stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        epsilon: Option<Scalar>,
        /// ζ = v²/k (k-ε-ζ-f; φ of the elliptic-blending k-ε; default 2/3)
        /// (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        zeta: Option<f64>,
    },
    /// Prescribed mass flow (kg/s) into the domain, normal to the faces,
    /// spread by face area.
    MassFlowInlet {
        /// Mass flow into the domain (kg/s): a number, or an expression of
        /// t (D53).
        mass_flow: Scalar,
        /// Turbulence values (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
        /// k (m²/s²) directly: a number or an expression (D53) (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        k: Option<Scalar>,
        /// ω (1/s) directly, with `k` (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        omega: Option<Scalar>,
        /// ε (m²/s³) directly, with `k` (the k-ε family; stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        epsilon: Option<Scalar>,
        /// ζ = v²/k (k-ε-ζ-f; φ of the elliptic-blending k-ε; default 2/3)
        /// (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        zeta: Option<f64>,
    },
    /// Prescribed static pressure, inflow normal to the face.
    StaticPressureInlet {
        /// Static pressure (Pa): a number or an expression (D53).
        pressure: Scalar,
        /// What the value means with gravity: the reduced pressure p′
        /// (default) or the static pressure p (D153, stage 7b1).
        #[serde(
            default,
            skip_serializing_if = "crate::buoyancy::PressureKind::is_default"
        )]
        pressure_kind: crate::buoyancy::PressureKind,
        /// Turbulence values (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
        /// k (m²/s²) directly: a number or an expression (D53) (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        k: Option<Scalar>,
        /// ω (1/s) directly, with `k` (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        omega: Option<Scalar>,
        /// ε (m²/s³) directly, with `k` (the k-ε family; stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        epsilon: Option<Scalar>,
        /// ζ = v²/k (k-ε-ζ-f; φ of the elliptic-blending k-ε; default 2/3)
        /// (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        zeta: Option<f64>,
    },
    /// Prescribed total pressure, inflow normal to the face.
    TotalPressureInlet {
        /// Total pressure p₀ = p + ½ρ|u|² (Pa): a number or an expression
        /// (D53).
        total_pressure: Scalar,
        /// What the value means with gravity: the reduced pressure p′
        /// (default) or the static pressure p (D153, stage 7b1).
        #[serde(
            default,
            skip_serializing_if = "crate::buoyancy::PressureKind::is_default"
        )]
        pressure_kind: crate::buoyancy::PressureKind,
        /// Turbulence values (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
        /// k (m²/s²) directly: a number or an expression (D53) (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        k: Option<Scalar>,
        /// ω (1/s) directly, with `k` (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        omega: Option<Scalar>,
        /// ε (m²/s³) directly, with `k` (the k-ε family; stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        epsilon: Option<Scalar>,
        /// ζ = v²/k (k-ε-ζ-f; φ of the elliptic-blending k-ε; default 2/3)
        /// (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        zeta: Option<f64>,
    },
    /// Prescribed static pressure, velocity extrapolated.
    StaticPressureOutlet {
        /// Static pressure (Pa): a number or an expression (D53), for
        /// example the hydrostatic pressure of a free outlet.
        #[serde(default)]
        pressure: Scalar,
        /// What the value means with gravity: the reduced pressure p′
        /// (default) or the static pressure p (D153, stage 7b1).
        #[serde(
            default,
            skip_serializing_if = "crate::buoyancy::PressureKind::is_default"
        )]
        pressure_kind: crate::buoyancy::PressureKind,
        /// Treatment of the faces where the flow enters.
        #[serde(default)]
        backflow: Backflow,
        /// Turbulence values of the backflow (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale of the backflow (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio of the backflow (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
        /// k (m²/s²) directly: a number or an expression (D53) (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        k: Option<Scalar>,
        /// ω (1/s) directly, with `k` (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        omega: Option<Scalar>,
        /// ε (m²/s³) directly, with `k` (the k-ε family; stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        epsilon: Option<Scalar>,
        /// ζ = v²/k (k-ε-ζ-f; φ of the elliptic-blending k-ε; default 2/3)
        /// (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        zeta: Option<f64>,
    },
    /// Prescribed total pressure, velocity extrapolated.
    TotalPressureOutlet {
        /// Total pressure (Pa): a number or an expression (D53).
        total_pressure: Scalar,
        /// What the value means with gravity: the reduced pressure p′
        /// (default) or the static pressure p (D153, stage 7b1).
        #[serde(
            default,
            skip_serializing_if = "crate::buoyancy::PressureKind::is_default"
        )]
        pressure_kind: crate::buoyancy::PressureKind,
        /// Treatment of the faces where the flow enters.
        #[serde(default)]
        backflow: Backflow,
        /// Turbulence values of the backflow (stage 6).
        #[serde(default)]
        turbulence_intensity: Option<f64>,
        /// Turbulence length scale of the backflow (stage 6).
        #[serde(default)]
        turbulence_length_scale: Option<f64>,
        /// Viscosity ratio of the backflow (stage 6).
        #[serde(default)]
        viscosity_ratio: Option<f64>,
        /// k (m²/s²) directly: a number or an expression (D53) (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        k: Option<Scalar>,
        /// ω (1/s) directly, with `k` (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        omega: Option<Scalar>,
        /// ε (m²/s³) directly, with `k` (the k-ε family; stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        epsilon: Option<Scalar>,
        /// ζ = v²/k (k-ε-ζ-f; φ of the elliptic-blending k-ε; default 2/3)
        /// (stage 6).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        zeta: Option<f64>,
    },
    /// Non-reflecting outlet of a transient run: the velocity is convected
    /// out of the domain, ∂u/∂t + U_c ∂u/∂n = 0.
    ConvectiveOutlet {
        /// Convection velocity U_c (m/s); by default the mean outflow
        /// speed that carries the prescribed inflow.
        #[serde(default)]
        convection_velocity: Option<f64>,
    },
    /// Wall without friction: no shear stress and no flow through it.
    SlipWall,
    /// Symmetry plane (also the two side planes of a 2D case).
    Symmetry,
    /// One of the two side planes of an axisymmetric wedge mesh (D3).
    Wedge,
    /// Translational periodic pair (D29, D109 (1)).
    Periodic {
        /// The partner boundary.
        #[serde(default)]
        partner: Option<String>,
        /// Mass flow through the pair (kg/s), from this boundary towards
        /// its partner (into the domain here); held by a uniform body
        /// force along the translation that the solver adjusts.
        #[serde(default)]
        mass_flow: Option<f64>,
        /// Mean pressure gradient (Pa/m) that drives the flow: the mean
        /// pressure falls by this much per metre from this boundary towards
        /// its partner (a uniform body force along the translation).
        #[serde(default)]
        pressure_gradient: Option<f64>,
        /// Translation from this boundary to its partner (m); by default
        /// the difference of the area-weighted face centroids.
        #[serde(default)]
        translation: Option<[f64; 3]>,
    },
}

impl BoundaryCondition {
    /// The `type` name.
    pub fn type_name(&self) -> &'static str {
        match self {
            Self::Wall { .. } => "wall",
            Self::VelocityInlet { .. } => "velocity-inlet",
            Self::MassFlowInlet { .. } => "mass-flow-inlet",
            Self::StaticPressureInlet { .. } => "static-pressure-inlet",
            Self::TotalPressureInlet { .. } => "total-pressure-inlet",
            Self::StaticPressureOutlet { .. } => "static-pressure-outlet",
            Self::TotalPressureOutlet { .. } => "total-pressure-outlet",
            Self::ConvectiveOutlet { .. } => "convective-outlet",
            Self::SlipWall => "slip-wall",
            Self::Symmetry => "symmetry",
            Self::Wedge => "wedge",
            Self::Periodic { .. } => "periodic",
        }
    }

    /// True for the conditions that fix the pressure.
    pub fn is_pressure(&self) -> bool {
        matches!(
            self,
            Self::StaticPressureInlet { .. }
                | Self::TotalPressureInlet { .. }
                | Self::StaticPressureOutlet { .. }
                | Self::TotalPressureOutlet { .. }
        )
    }

    /// True for a boundary through which fluid enters or leaves (an
    /// inlet, an outlet, a pressure boundary).
    pub fn is_open(&self) -> bool {
        matches!(
            self,
            Self::VelocityInlet { .. } | Self::MassFlowInlet { .. } | Self::ConvectiveOutlet { .. }
        ) || self.is_pressure()
    }

    /// True for a periodic boundary.
    pub fn is_periodic(&self) -> bool {
        matches!(self, Self::Periodic { .. })
    }

    /// True for a convective outlet.
    pub fn is_convective_outlet(&self) -> bool {
        matches!(self, Self::ConvectiveOutlet { .. })
    }

    /// True for a mass-flow inlet.
    pub fn is_mass_flow_inlet(&self) -> bool {
        matches!(self, Self::MassFlowInlet { .. })
    }

    /// True for a wall with a mass flux.
    pub fn has_wall_mass_flux(&self) -> bool {
        matches!(
            self,
            Self::Wall {
                mass_flux: Some(_),
                ..
            }
        )
    }

    /// True for the pressure conditions meant for outflow.
    pub fn is_pressure_outlet(&self) -> bool {
        matches!(
            self,
            Self::StaticPressureOutlet { .. } | Self::TotalPressureOutlet { .. }
        )
    }

    /// The turbulence values of an inlet or of the backflow of an outlet.
    pub fn turbulence(&self) -> InletTurbulence {
        match self {
            Self::VelocityInlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                k,
                omega,
                epsilon,
                zeta,
                ..
            }
            | Self::MassFlowInlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                k,
                omega,
                epsilon,
                zeta,
                ..
            }
            | Self::StaticPressureInlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                k,
                omega,
                epsilon,
                zeta,
                ..
            }
            | Self::TotalPressureInlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                k,
                omega,
                epsilon,
                zeta,
                ..
            }
            | Self::StaticPressureOutlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                k,
                omega,
                epsilon,
                zeta,
                ..
            }
            | Self::TotalPressureOutlet {
                turbulence_intensity,
                turbulence_length_scale,
                viscosity_ratio,
                k,
                omega,
                epsilon,
                zeta,
                ..
            } => InletTurbulence {
                turbulence_intensity: *turbulence_intensity,
                turbulence_length_scale: *turbulence_length_scale,
                viscosity_ratio: *viscosity_ratio,
                k: k.clone(),
                omega: omega.clone(),
                epsilon: epsilon.clone(),
                zeta: *zeta,
            },
            _ => InletTurbulence::default(),
        }
    }

    /// The roughness k_s of a wall (m; 0 for other conditions).
    pub fn roughness(&self) -> f64 {
        match self {
            Self::Wall { roughness, .. } => roughness.unwrap_or(0.0),
            _ => 0.0,
        }
    }

    /// C_s of the roughness shift of a wall (default 0.5).
    pub fn roughness_constant(&self) -> f64 {
        match self {
            Self::Wall {
                roughness_constant, ..
            } => roughness_constant.unwrap_or(0.5),
            _ => 0.5,
        }
    }

    /// True for the conditions with no flow through them and no shear
    /// stress (symmetry planes, slip walls, wedge planes).
    pub fn is_symmetry_like(&self) -> bool {
        matches!(self, Self::Symmetry | Self::SlipWall | Self::Wedge)
    }

    /// The table of a table profile.
    pub fn table(&self) -> Option<&ProfileTable> {
        match self {
            Self::VelocityInlet { table_data, .. } => table_data.as_deref(),
            _ => None,
        }
    }

    /// The table of a table profile with its coordinate, the origin and
    /// the direction of its axis (the turbulence columns of the inlet,
    /// stage 6 specification 7.1).
    pub fn table_profile(&self) -> Option<TableProfile> {
        match self {
            Self::VelocityInlet {
                profile: Some(InletProfile::Table),
                table_coordinate: Some(c),
                table_origin,
                table_direction,
                table_data: Some(t),
                ..
            } => Some(TableProfile {
                table: t.clone(),
                coordinate: *c,
                origin: *table_origin,
                direction: *table_direction,
            }),
            _ => None,
        }
    }

    /// The file of a table profile.
    pub fn table_path(&self) -> Option<&Path> {
        match self {
            Self::VelocityInlet {
                profile: Some(InletProfile::Table),
                table: Some(t),
                ..
            } => Some(t),
            _ => None,
        }
    }

    /// Makes the paths of the condition relative to `dir` (the directory
    /// of the case file).
    pub fn resolve_paths(&mut self, dir: &Path) {
        if let Self::VelocityInlet { table: Some(t), .. } = self {
            *t = dir.join(&*t);
        }
    }

    /// Reads the table of a table profile, `table` relative to `dir`;
    /// nothing for other conditions.
    pub fn load_table(&mut self, dir: &Path) -> Result<()> {
        if let Self::VelocityInlet {
            profile: Some(InletProfile::Table),
            table: Some(t),
            table_data,
            ..
        } = self
        {
            let path = dir.join(&*t);
            let data = ProfileTable::read(&path)
                .and_then(|d| {
                    d.columns()
                        .map_err(|e| format!("{}: {e}", path.display()))
                        .map(|_| d)
                })
                .map_err(|e| {
                    diag::err(diag::Diagnostic::new(
                        diag::BC_TABLE,
                        format!("table profile: {e}"),
                    ))
                })?;
            *table_data = Some(Arc::new(data));
        }
        Ok(())
    }
}

fn zero_vector() -> Vector {
    crate::expr::vector([0.0; 3])
}

/// The compiled values of one condition (D53): every number or expression
/// of the condition that sets a face value.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Values {
    /// Wall or inlet velocity.
    pub velocity: Option<[Value; 3]>,
    /// Inlet normal speed.
    pub normal_velocity: Option<Value>,
    /// Wall mass flux (W3).
    pub mass_flux: Option<Value>,
    /// Inlet mass flow.
    pub mass_flow: Option<Value>,
    /// Static or total pressure.
    pub pressure: Option<Value>,
}

impl Values {
    /// The variables the values read.
    pub fn uses(&self) -> Vars {
        let mut u = Vars::NONE;
        for v in self.velocity.iter().flatten() {
            u = u.with(v.uses());
        }
        for v in [
            &self.normal_velocity,
            &self.mass_flux,
            &self.mass_flow,
            &self.pressure,
        ]
        .into_iter()
        .flatten()
        {
            u = u.with(v.uses());
        }
        u
    }
}

/// The compiled values of every condition, by boundary.
pub type Compiled = BTreeMap<String, Values>;

/// A cell zone that turns about an axis: a moving reference frame (D55) or
/// the rotor of a sliding mesh (D182). Its walls move with it, also normal
/// to themselves (a blade pushes the fluid): their normal velocity is that
/// of the zone, **Ω** × (**x** − **x**₀), whatever the `velocity` of the
/// wall says.
#[derive(Debug, Clone, PartialEq)]
pub struct Turning {
    /// The cell zone (`None`: every cell).
    pub zone: Option<String>,
    /// **Ω** (rad/s, along the axis, right-handed).
    pub omega: [f64; 3],
    /// A point on the axis (m).
    pub origin: [f64; 3],
    /// The unit vector of the axis, along **Ω** (the given `axis` when
    /// Ω = 0): the circumferential direction of the turbomachinery fields
    /// (D183) is the direction of the rotation.
    pub axis: [f64; 3],
}

impl Turning {
    /// True if local cell `c` of `mesh` is in the zone.
    pub fn contains(&self, mesh: &Mesh, c: usize) -> bool {
        self.zone.as_ref().is_none_or(|z| {
            mesh.zones
                .get(mesh.cell_zone[c] as usize)
                .is_some_and(|n| n == z)
        })
    }

    /// The first zone of `turning` that holds local cell `c`.
    pub fn of<'t>(turning: &'t [Turning], mesh: &Mesh, c: usize) -> Option<&'t Turning> {
        turning.iter().find(|t| t.contains(mesh, c))
    }

    /// The velocity of the zone at `x`, **Ω** × (**x** − **x**₀).
    pub fn at(&self, x: [f64; 3]) -> [f64; 3] {
        cross(self.omega, sub(x, self.origin))
    }

    /// The velocity of the zone of local cell `c` at `x` among `turning`
    /// (zero outside all of them).
    pub fn velocity(turning: &[Turning], mesh: &Mesh, c: usize, x: [f64; 3]) -> [f64; 3] {
        Turning::of(turning, mesh, c).map_or([0.0; 3], |t| t.at(x))
    }
}

/// Where and when the conditions are evaluated: the exact solution, the
/// compiled expressions, the simulated time and the turning zones.
#[derive(Debug, Clone, Copy)]
pub struct At<'a> {
    /// The `[analytic]` solution (at the time `t`).
    pub analytic: Option<&'a Analytic>,
    /// The compiled values (an empty map: every value is a number).
    pub values: &'a Compiled,
    /// Simulated time (s).
    pub t: f64,
    /// The zones that turn (`[[mrf]]`, `[[sliding]]`; see
    /// [`FlowCase::turning`](crate::case::FlowCase::turning)).
    pub turning: &'a [Turning],
    /// The hydrostatic reference of gravity (stage 7b1): the values of
    /// boundaries with `pressure_kind = "static"` become reduced pressures.
    pub hydrostatic: Option<crate::buoyancy::Hydrostatic>,
}

static NO_VALUES: std::sync::LazyLock<Compiled> = std::sync::LazyLock::new(Compiled::new);

impl<'a> At<'a> {
    /// Numbers only, at t = 0, with the exact solution `analytic`.
    pub fn numbers(analytic: Option<&'a Analytic>) -> Self {
        Self {
            analytic,
            values: &NO_VALUES,
            t: 0.0,
            turning: &[],
            hydrostatic: None,
        }
    }

    fn of(&self, boundary: &str) -> Option<&'a Values> {
        self.values.get(boundary)
    }
}

/// A scalar of a condition at point `x`: its compiled value, or the number.
fn sval(s: &Scalar, v: Option<&Value>, x: [f64; 3], t: f64) -> f64 {
    match v {
        Some(v) => v.at(x, t),
        None => s.number().unwrap_or(f64::NAN),
    }
}

/// A vector of a condition at point `x`.
fn vval(s: &Vector, v: Option<&[Value; 3]>, x: [f64; 3], t: f64) -> [f64; 3] {
    std::array::from_fn(|i| sval(&s[i], v.map(|v| &v[i]), x, t))
}

impl BoundaryCondition {
    /// Compiles the numbers and expressions of the condition of `boundary`
    /// (D53): a transient case allows `t`; a mass flow is one number for
    /// the boundary (only `t`).
    pub fn compile(&self, boundary: &str, env: &Env, transient: bool) -> Result<Values> {
        let field = if transient { Vars::ALL } else { Vars::SPACE };
        let boundary_total = if transient { Vars::T } else { Vars::NONE };
        let key = |k: &str| format!("boundary.{boundary}.{k}");
        let one = |s: &Scalar, k: &str, a: Vars| crate::expr::value(s, env, a, &key(k));
        let mut v = Values::default();
        match self {
            BoundaryCondition::Wall {
                velocity,
                mass_flux,
                ..
            } => {
                v.velocity = Some(crate::expr::vector_value(
                    velocity,
                    env,
                    field,
                    &key("velocity"),
                )?);
                if let Some(m) = mass_flux {
                    v.mass_flux = Some(one(m, "mass_flux", field)?);
                }
            }
            BoundaryCondition::VelocityInlet {
                velocity,
                normal_velocity,
                ..
            } => {
                if let Some(u) = velocity {
                    v.velocity = Some(crate::expr::vector_value(u, env, field, &key("velocity"))?);
                }
                if let Some(u) = normal_velocity {
                    v.normal_velocity = Some(one(u, "normal_velocity", field)?);
                }
            }
            BoundaryCondition::MassFlowInlet { mass_flow, .. } => {
                v.mass_flow = Some(one(mass_flow, "mass_flow", boundary_total)?);
            }
            BoundaryCondition::StaticPressureInlet { pressure, .. }
            | BoundaryCondition::StaticPressureOutlet { pressure, .. } => {
                v.pressure = Some(one(pressure, "pressure", field)?);
            }
            BoundaryCondition::TotalPressureInlet { total_pressure, .. }
            | BoundaryCondition::TotalPressureOutlet { total_pressure, .. } => {
                v.pressure = Some(one(total_pressure, "total_pressure", field)?);
            }
            _ => {}
        }
        Ok(v)
    }

    /// The mass flux into the domain through a wall face at `x`.
    pub fn wall_mass_flux(&self, vals: Option<&Values>, x: [f64; 3], t: f64) -> f64 {
        match self {
            BoundaryCondition::Wall {
                mass_flux: Some(m), ..
            } => sval(m, vals.and_then(|v| v.mass_flux.as_ref()), x, t),
            _ => 0.0,
        }
    }

    /// `pressure_kind` of a pressure boundary (`None` for the others).
    pub fn pressure_kind(&self) -> Option<crate::buoyancy::PressureKind> {
        match self {
            BoundaryCondition::StaticPressureInlet { pressure_kind, .. }
            | BoundaryCondition::StaticPressureOutlet { pressure_kind, .. }
            | BoundaryCondition::TotalPressureInlet { pressure_kind, .. }
            | BoundaryCondition::TotalPressureOutlet { pressure_kind, .. } => Some(*pressure_kind),
            _ => None,
        }
    }

    /// The prescribed pressure or total pressure at `x`.
    pub fn pressure_at(&self, vals: Option<&Values>, x: [f64; 3], t: f64) -> f64 {
        let p = vals.and_then(|v| v.pressure.as_ref());
        match self {
            BoundaryCondition::StaticPressureInlet { pressure, .. }
            | BoundaryCondition::StaticPressureOutlet { pressure, .. } => sval(pressure, p, x, t),
            BoundaryCondition::TotalPressureInlet { total_pressure, .. }
            | BoundaryCondition::TotalPressureOutlet { total_pressure, .. } => {
                sval(total_pressure, p, x, t)
            }
            _ => 0.0,
        }
    }

    /// The mass flow of a mass-flow inlet at time `t`.
    pub fn mass_flow_at(&self, vals: Option<&Values>, t: f64) -> f64 {
        match self {
            BoundaryCondition::MassFlowInlet { mass_flow, .. } => sval(
                mass_flow,
                vals.and_then(|v| v.mass_flow.as_ref()),
                [0.0; 3],
                t,
            ),
            _ => 0.0,
        }
    }

    /// The wall velocity at `x`.
    pub fn wall_velocity(&self, vals: Option<&Values>, x: [f64; 3], t: f64) -> [f64; 3] {
        match self {
            BoundaryCondition::Wall { velocity, .. } => {
                vval(velocity, vals.and_then(|v| v.velocity.as_ref()), x, t)
            }
            _ => [0.0; 3],
        }
    }
}

/// Compiles every condition (D53).
pub fn compile_all(
    bcs: &BTreeMap<String, BoundaryCondition>,
    env: &Env,
    transient: bool,
) -> Result<Compiled> {
    bcs.iter()
        .map(|(n, b)| Ok((n.clone(), b.compile(n, env, transient)?)))
        .collect()
}

/// Severity of a finding of [`check`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Level {
    /// Information (for example how the pressure level is fixed).
    Note,
    /// Allowed but probably not intended.
    Warning,
    /// The run cannot start.
    Error,
}

/// One finding of [`check`]: what is wrong, where, and what to change.
#[derive(Debug, Clone, Serialize)]
pub struct Finding {
    /// Severity.
    pub level: Level,
    /// Diagnostic code (D78).
    pub code: &'static str,
    /// The boundary it concerns, if one.
    pub boundary: Option<String>,
    /// What is wrong and what to change.
    pub message: String,
}

impl Finding {
    /// As a diagnostic (D78).
    pub fn diagnostic(&self) -> diag::Diagnostic {
        let mut d =
            diag::Diagnostic::new(self.code, self.message.clone()).severity(match self.level {
                Level::Note => diag::Severity::Note,
                Level::Warning => diag::Severity::Warning,
                Level::Error => diag::Severity::Error,
            });
        if let Some(p) = &self.boundary {
            d = d.boundary(p.clone()).key(format!("boundary.{p}"));
        }
        d
    }
}

impl fmt::Display for Finding {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let l = match self.level {
            Level::Note => "note",
            Level::Warning => "warning",
            Level::Error => "error",
        };
        match &self.boundary {
            Some(p) => write!(f, "{l} {}: boundary `{p}`: {}", self.code, self.message),
            None => write!(f, "{l} {}: {}", self.code, self.message),
        }
    }
}

/// Levenshtein distance, for the "did you mean" hints.
pub fn edit_distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut cur = vec![i + 1; b.len() + 1];
        for (j, &cb) in b.iter().enumerate() {
            cur[j + 1] = (prev[j] + usize::from(ca != cb))
                .min(prev[j + 1] + 1)
                .min(cur[j] + 1);
        }
        prev = cur;
    }
    prev[b.len()]
}

/// The candidate closest to `word`, if it is close enough to be a typo.
pub fn did_you_mean<'a>(
    word: &str,
    candidates: impl IntoIterator<Item = &'a str>,
) -> Option<&'a str> {
    candidates
        .into_iter()
        .map(|c| (edit_distance(word, c), c))
        .filter(|&(d, c)| d <= 2.max(c.len() / 3))
        .min_by_key(|&(d, _)| d)
        .map(|(_, c)| c)
}

fn hint(word: &str, candidates: &[&str]) -> String {
    match did_you_mean(word, candidates.iter().copied()) {
        Some(c) => format!(" (did you mean `{c}`?)"),
        None => String::new(),
    }
}

/// Mass flow (kg/s) out of the domain through the faces of every boundary whose
/// condition prescribes the flux, for the initial velocity of the
/// conditions (walls and symmetry planes give 0).
fn prescribed_outflow(
    mesh: &Mesh,
    geo: &Geometry,
    bcs: &BTreeMap<String, BoundaryCondition>,
    rho: Density,
    at: At,
) -> (BTreeMap<String, f64>, f64) {
    let mut out = BTreeMap::new();
    let mut gross = 0.0;
    for p in &mesh.boundaries {
        let Some(bc) = bcs.get(&p.name) else { continue };
        let area: f64 = p.range().map(|f| norm(geo.face_area[f])).exact_sum();
        let ext = boundary_extent(mesh, p);
        let q = match bc {
            BoundaryCondition::VelocityInlet { .. } => p
                .range()
                .map(|f| {
                    let q = rho.at(f)
                        * dot(
                            inlet_velocity(bc, at.of(&p.name), geo, f, at, &ext)
                                .unwrap_or([0.0; 3]),
                            geo.face_area[f],
                        );
                    gross += q.abs();
                    q
                })
                .exact_sum(),
            BoundaryCondition::MassFlowInlet { .. } => {
                let mass_flow = bc.mass_flow_at(at.of(&p.name), at.t);
                gross += mass_flow.abs();
                if area > 0.0 { -mass_flow } else { 0.0 }
            }
            BoundaryCondition::Wall {
                mass_flux: Some(_), ..
            } => p
                .range()
                .map(|f| {
                    let q = -bc.wall_mass_flux(at.of(&p.name), geo.face_centroid[f], at.t)
                        * norm(geo.face_area[f]);
                    gross += q.abs();
                    q
                })
                .exact_sum(),
            _ => continue,
        };
        out.insert(p.name.clone(), q);
    }
    (out, gross)
}

/// The prescribed velocity of face `f` of a velocity inlet whose boundary has
/// the extent `ext` (None for other conditions).
pub fn inlet_velocity(
    bc: &BoundaryCondition,
    vals: Option<&Values>,
    geo: &Geometry,
    f: usize,
    at: At,
    ext: &Extent,
) -> Option<[f64; 3]> {
    let analytic = at.analytic;
    let BoundaryCondition::VelocityInlet {
        velocity,
        normal_velocity,
        profile,
        max_velocity,
        profile_axes,
        table_coordinate,
        table_origin,
        table_direction,
        table_data,
        profile_extent,
        ..
    } = bc
    else {
        return None;
    };
    // on ranks the extent of the whole boundary, not of this rank's part
    let ext = profile_extent.as_ref().unwrap_or(ext);
    if let Some(InletProfile::Analytic) = profile {
        return analytic.map(|a| a.velocity(geo.face_centroid[f]));
    }
    if let Some(InletProfile::Table) = profile {
        let (t, c) = (table_data.as_ref()?, (*table_coordinate)?);
        let x = table_position(geo.face_centroid[f], c, *table_origin, *table_direction)?;
        let (v, _) = t.at(x);
        let s = geo.face_area[f];
        return Some(match t.velocity_columns() {
            VelocityColumns::Normal(i) => scale(s, -v[i] / norm(s)),
            VelocityColumns::Components([a, b, c]) => [v[a], v[b], v[c]],
        });
    }
    if let Some(InletProfile::Parabolic) = profile {
        let x = geo.face_centroid[f];
        let mut un = (*max_velocity)?;
        for a in profile_axes {
            let i = *a as usize;
            let (lo, hi) = (ext.0[i], ext.1[i]);
            let eta = if hi > lo {
                (x[i] - lo) / (hi - lo)
            } else {
                0.5
            };
            un *= 4.0 * eta * (1.0 - eta);
        }
        let s = geo.face_area[f];
        return Some(scale(s, -un / norm(s)));
    }
    let x = geo.face_centroid[f];
    if let Some(v) = velocity {
        return Some(vval(v, vals.and_then(|v| v.velocity.as_ref()), x, at.t));
    }
    let s = geo.face_area[f];
    let n = scale(s, 1.0 / norm(s));
    normal_velocity.as_ref().map(|un| {
        scale(
            n,
            -sval(un, vals.and_then(|v| v.normal_velocity.as_ref()), x, at.t),
        )
    })
}

/// The table coordinate of point `x`.
pub fn table_position(
    x: [f64; 3],
    c: TableCoordinate,
    origin: [f64; 3],
    direction: Option<[f64; 3]>,
) -> Option<f64> {
    Some(match c {
        TableCoordinate::X => x[0],
        TableCoordinate::Y => x[1],
        TableCoordinate::Z => x[2],
        TableCoordinate::Radius => {
            let d = direction?;
            let l = norm(d);
            if l == 0.0 {
                return None;
            }
            let d = scale(d, 1.0 / l);
            let r = sub(x, origin);
            norm(sub(r, scale(d, dot(r, d))))
        }
    })
}

/// The plane of a planar boundary: unit normal and a point, if every face
/// normal agrees with the first to `tol` (the distance of their tips).
fn boundary_plane(
    geo: &Geometry,
    p: &foehncraft_mesh::Boundary,
    tol: f64,
) -> Option<([f64; 3], [f64; 3])> {
    if p.len == 0 {
        return None;
    }
    let s0 = geo.face_area[p.start];
    let n0 = scale(s0, 1.0 / norm(s0));
    let planar = p.range().all(|f| {
        let s = geo.face_area[f];
        norm(sub(scale(s, 1.0 / norm(s)), n0)) <= tol
    });
    planar.then(|| (n0, geo.face_centroid[p.start]))
}

/// What the check of the conditions needs to know about the case.
#[derive(Debug, Clone, Default)]
pub struct CheckContext {
    /// Laminar flow (turbulence values unused).
    pub laminar: bool,
    /// A transient case (`[time]`).
    pub transient: bool,
    /// The potential-flow start (D52).
    pub potential_flow: bool,
    /// A body force acts (`flow.body_force` not zero).
    pub body_force: bool,
    /// The turbulence model solves ζ = v²/k (k-ε-ζ-f): `zeta` at inlets
    /// is used.
    pub zeta: bool,
    /// Boundaries that are no boundaries of the solved mesh: the two sides
    /// of a sliding interface that is imprinted in every step (D182, variant
    /// B). They need no condition.
    pub interface: Vec<String>,
}

/// How the flow through a periodic pair is driven (D29).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Drive {
    /// Not driven (a body force of the case, or a decaying flow).
    None,
    /// The mass flow (kg/s) from side `a` towards side `b`.
    MassFlow(f64),
    /// The mean pressure gradient (Pa/m), the pressure falling from side
    /// `a` towards side `b`.
    PressureGradient(f64),
}

/// A periodic pair of a case: the boundaries (side `a` carries the driving,
/// so that the flow goes along the translation from `a` to `b`) and the
/// driving.
#[derive(Debug, Clone, PartialEq)]
pub struct PeriodicCase {
    /// The pair for the mesh.
    pub spec: foehncraft_mesh::periodic::PairSpec,
    /// The driving.
    pub drive: Drive,
}

/// The periodic pairs of the conditions `bcs`: every periodic boundary
/// names a periodic partner that names it back (one of the two may leave
/// `partner` out), the driving and the translation are given on at most
/// one side. Errors carry the code BC030.
pub fn periodic_pairs(bcs: &BTreeMap<String, BoundaryCondition>) -> Result<Vec<PeriodicCase>> {
    let err = |name: &str, msg: String| {
        diag::err(
            diag::Diagnostic::new(diag::BC_PERIODIC, format!("[boundary.{name}]: {msg}"))
                .key(format!("boundary.{name}")),
        )
    };
    struct Side<'a> {
        partner: Option<&'a str>,
        drive: Drive,
        translation: Option<[f64; 3]>,
    }
    let mut sides: BTreeMap<&str, Side> = BTreeMap::new();
    for (name, bc) in bcs {
        if let BoundaryCondition::Periodic {
            partner,
            mass_flow,
            pressure_gradient,
            translation,
        } = bc
        {
            let drive = match (mass_flow, pressure_gradient) {
                (Some(_), Some(_)) => {
                    return Err(err(
                        name,
                        "give either `mass_flow` or `pressure_gradient`, not both".into(),
                    ));
                }
                (Some(m), None) => Drive::MassFlow(*m),
                (None, Some(g)) => Drive::PressureGradient(*g),
                (None, None) => Drive::None,
            };
            for v in [mass_flow, pressure_gradient].into_iter().flatten() {
                if !v.is_finite() {
                    return Err(err(name, "the driving must be a finite number".into()));
                }
            }
            sides.insert(
                name,
                Side {
                    partner: partner.as_deref(),
                    drive,
                    translation: *translation,
                },
            );
        }
    }
    let mut out = Vec::new();
    let mut done: Vec<&str> = Vec::new();
    for (&name, s) in &sides {
        if done.contains(&name) {
            continue;
        }
        // the partner: named here, or the periodic boundary that names this one
        let named_by: Vec<&str> = sides
            .iter()
            .filter(|(n, o)| **n != name && o.partner == Some(name))
            .map(|(n, _)| *n)
            .collect();
        let partner = match (s.partner, named_by.as_slice()) {
            (Some(p), _) => p,
            (None, [p]) => *p,
            (None, []) => {
                return Err(err(
                    name,
                    "a periodic boundary needs a `partner` (the boundary it repeats on)".into(),
                ));
            }
            (None, _) => {
                return Err(err(
                    name,
                    format!(
                        "named as partner by several boundaries: {}",
                        named_by.join(", ")
                    ),
                ));
            }
        };
        if partner == name {
            return Err(err(
                name,
                "a boundary cannot be its own periodic partner".into(),
            ));
        }
        let Some(o) = sides.get(partner) else {
            return Err(err(
                name,
                format!(
                    "the partner `{partner}` {}",
                    if bcs.contains_key(partner) {
                        "is not of type periodic"
                    } else {
                        "has no boundary condition"
                    }
                ),
            ));
        };
        if o.partner.is_some_and(|p| p != name) {
            return Err(err(
                partner,
                format!(
                    "its partner is `{}`, but `{name}` names it",
                    o.partner.unwrap_or("")
                ),
            ));
        }
        if done.contains(&partner) {
            return Err(err(name, format!("`{partner}` is in another pair already")));
        }
        if s.drive != Drive::None && o.drive != Drive::None {
            return Err(err(
                name,
                format!("the driving is given on both `{name}` and `{partner}`: give it on one"),
            ));
        }
        if s.translation.is_some() && o.translation.is_some() {
            return Err(err(
                name,
                format!("`translation` is given on both `{name}` and `{partner}`: give it on one"),
            ));
        }
        // side a: the one with the driving (else with the translation,
        // else this one, the first in name order)
        let swap = o.drive != Drive::None || (s.drive == Drive::None && o.translation.is_some());
        let (a, b, sa, sb) = if swap {
            (partner, name, o, s)
        } else {
            (name, partner, s, o)
        };
        let translation = sa.translation.or(sb.translation.map(|t| t.map(|v| -v)));
        out.push(PeriodicCase {
            spec: foehncraft_mesh::periodic::PairSpec {
                a: a.to_string(),
                b: b.to_string(),
                translation,
            },
            drive: sa.drive,
        });
        done.push(name);
        done.push(partner);
    }
    Ok(out)
}

/// Checks the combination of boundary conditions on the whole domain (D30)
/// before the run. Errors stop the run; warnings and notes go to the log.
/// Every finding names the boundary and says what to change.
pub fn check<'a>(
    mesh: &Mesh,
    geo: &Geometry,
    bcs: &BTreeMap<String, BoundaryCondition>,
    rho: impl Into<Density<'a>>,
    at: At,
    ctx: CheckContext,
) -> Vec<Finding> {
    let analytic = at.analytic;
    let laminar = ctx.laminar;
    let mut out = Vec::new();
    let mut push = |level, code, boundary: Option<&str>, message: String| {
        out.push(Finding {
            level,
            code,
            boundary: boundary.map(str::to_string),
            message,
        })
    };
    let names: Vec<&str> = mesh.boundaries.iter().map(|p| p.name.as_str()).collect();
    // the interface boundaries that `[[sliding]]` names are on the mesh
    for n in &ctx.interface {
        if !names.contains(&n.as_str()) {
            let hint = did_you_mean(n, names.iter().copied())
                .map(|c| format!(" (did you mean `{c}`?)"))
                .unwrap_or_default();
            push(
                Level::Error,
                diag::SLIDING_INTERFACE,
                Some(n),
                format!(
                    "[[sliding]]: the mesh has no boundary `{n}` for the interface of the \
                     sliding mesh{hint}; its boundaries are {names:?}"
                ),
            );
        }
    }
    // every boundary has a condition and every condition a boundary
    for p in &mesh.boundaries {
        if !bcs.contains_key(&p.name) && !ctx.interface.contains(&p.name) {
            push(
                Level::Error,
                diag::BC_NO_CONDITION,
                Some(&p.name),
                format!(
                    "no boundary condition: add a table [boundary.{}] with a `type` \
                     (wall, slip-wall, velocity-inlet, mass-flow-inlet, static-pressure-inlet, \
                     total-pressure-inlet, static-pressure-outlet, total-pressure-outlet, \
                     convective-outlet, symmetry, wedge)",
                    p.name
                ),
            );
        }
    }
    let paired = |n: &str| {
        mesh.periodic
            .as_ref()
            .is_some_and(|p| p.pair_of(n).is_some())
    };
    if let Ok(pairs) = periodic_pairs(bcs)
        && !ctx.body_force
    {
        for pc in pairs.iter().filter(|p| p.drive == Drive::None) {
            push(
                Level::Warning,
                diag::BC_PERIODIC_UNDRIVEN,
                Some(&pc.spec.a),
                format!(
                    "nothing drives the flow through the periodic pair `{}`–`{}`: give \
                     `mass_flow` (kg/s) or `pressure_gradient` (Pa/m) on one of them, or a \
                     body force, or the flow decays to rest",
                    pc.spec.a, pc.spec.b
                ),
            );
        }
    }
    for name in bcs.keys() {
        if mesh.boundary(name).is_none() && !paired(name) {
            push(
                Level::Error,
                diag::BC_NO_BOUNDARY,
                Some(name),
                format!(
                    "the mesh has no such boundary{}; its boundaries are {}",
                    hint(name, &names),
                    names.join(", ")
                ),
            );
        }
    }
    let present: Vec<(&str, &BoundaryCondition)> = mesh
        .boundaries
        .iter()
        .filter_map(|p| bcs.get(&p.name).map(|b| (p.name.as_str(), b)))
        .collect();
    // conditions that are not available yet or need more input
    for &(name, bc) in &present {
        match bc {
            BoundaryCondition::Periodic { .. } => push(
                Level::Error,
                diag::BC_PERIODIC,
                Some(name),
                "this periodic boundary was not paired: give it a `partner` of type periodic \
                 that names it back"
                    .into(),
            ),
            BoundaryCondition::VelocityInlet {
                velocity,
                normal_velocity,
                profile,
                max_velocity,
                profile_axes,
                table,
                table_coordinate,
                table_direction,
                table_data,
                ..
            } => {
                if *profile == Some(InletProfile::Table) {
                    let missing = match (table, table_coordinate) {
                        (None, _) => Some("`table` (the file of the profile)"),
                        (_, None) => Some("`table_coordinate` (x, y, z or radius)"),
                        _ if *table_coordinate == Some(TableCoordinate::Radius)
                            && table_direction.is_none_or(|d| norm(d) == 0.0) =>
                        {
                            Some("`table_direction` (the axis of the radius)")
                        }
                        _ if table_data.is_none() => Some("the table data (the file was not read)"),
                        _ => None,
                    };
                    if let Some(m) = missing {
                        push(
                            Level::Error,
                            diag::BC_TABLE,
                            Some(name),
                            format!("`profile = \"table\"` needs {m}"),
                        );
                    }
                } else if table.is_some() || table_coordinate.is_some() || table_direction.is_some()
                {
                    push(
                        Level::Warning,
                        diag::BC_PROFILE_UNUSED,
                        Some(name),
                        "`table`, `table_coordinate` and `table_direction` are only used with \
                         `profile = \"table\"`"
                            .into(),
                    );
                }
                let given = usize::from(velocity.is_some())
                    + usize::from(normal_velocity.is_some())
                    + usize::from(profile.is_some());
                if given != 1 {
                    push(
                        Level::Error,
                        diag::BC_INLET_INPUT,
                        Some(name),
                        format!(
                            "a velocity-inlet takes exactly one of `velocity`, `normal_velocity` \
                             and `profile` ({given} given)"
                        ),
                    );
                }
                if *profile == Some(InletProfile::Analytic) && analytic.is_none() {
                    push(
                        Level::Error,
                        diag::BC_PROFILE_ANALYTIC,
                        Some(name),
                        "`profile = \"analytic\"` needs an [analytic] table".into(),
                    );
                }
                if *profile == Some(InletProfile::Parabolic)
                    && (max_velocity.is_none() || profile_axes.is_empty())
                {
                    push(
                        Level::Error,
                        diag::BC_PROFILE_PARABOLIC,
                        Some(name),
                        "`profile = \"parabolic\"` needs `max_velocity` (m/s) and \
                         `profile_axes` (for example [\"y\"])"
                            .into(),
                    );
                }
                if *profile != Some(InletProfile::Parabolic)
                    && (max_velocity.is_some() || !profile_axes.is_empty())
                {
                    push(
                        Level::Warning,
                        diag::BC_PROFILE_UNUSED,
                        Some(name),
                        "`max_velocity` and `profile_axes` are only used with \
                         `profile = \"parabolic\"`"
                            .into(),
                    );
                }
            }
            BoundaryCondition::ConvectiveOutlet {
                convection_velocity,
            } => {
                if !ctx.transient {
                    push(
                        Level::Error,
                        diag::BC_CONVECTIVE_STEADY,
                        Some(name),
                        "a convective outlet convects the flow out in time and needs a [time] \
                         table; for a steady case use static-pressure-outlet"
                            .into(),
                    );
                }
                if let Some(u) = convection_velocity
                    && !(*u > 0.0 && u.is_finite())
                {
                    push(
                        Level::Error,
                        diag::CASE_RANGE,
                        Some(name),
                        format!("convection_velocity = {u} m/s must be > 0"),
                    );
                }
            }
            BoundaryCondition::MassFlowInlet {
                mass_flow: Scalar::Number(mass_flow),
                ..
            } if *mass_flow < 0.0 => push(
                Level::Warning,
                diag::BC_NEGATIVE_MASS_FLOW,
                Some(name),
                format!(
                    "mass_flow = {mass_flow} kg/s is negative, so the flow leaves the domain \
                     here; for an outlet use static-pressure-outlet"
                ),
            ),
            _ => {}
        }
        let t = bc.turbulence().given();
        if laminar && !t.is_empty() {
            push(
                Level::Warning,
                diag::BC_TURBULENCE_UNUSED,
                Some(name),
                format!(
                    "{} unused: the flow is laminar (set `model` in [turbulence] for a \
                     turbulence model)",
                    t.join(", ")
                ),
            );
        } else if !laminar {
            if !t.is_empty() && !crate::turbulence::inlet::is_inlet(bc) && !bc.is_pressure_outlet()
            {
                push(
                    Level::Warning,
                    diag::BC_TURBULENCE_UNUSED_VALUE,
                    Some(name),
                    format!(
                        "{} unused: only inlets and the backflow of pressure outlets take \
                         turbulence values",
                        t.join(", ")
                    ),
                );
            }
            if !ctx.zeta
                && bc.turbulence().zeta.is_some()
                && (crate::turbulence::inlet::is_inlet(bc) || bc.is_pressure_outlet())
            {
                push(
                    Level::Warning,
                    diag::BC_TURBULENCE_UNUSED_VALUE,
                    Some(name),
                    "zeta unused: the turbulence model does not solve ζ = v²/k (only \
                     k-epsilon-zeta-f does)"
                        .into(),
                );
            }
            if let Some(z) = bc.turbulence().zeta
                && !(z > 0.0 && z <= 2.0)
            {
                push(
                    Level::Error,
                    diag::CASE_RANGE,
                    Some(name),
                    format!("zeta = {z} must lie in (0, 2] (ζ = v²/k ≤ 2)"),
                );
            }
            for (error, code, message) in crate::turbulence::inlet::check(bc) {
                let level = if error {
                    Level::Error
                } else if code == diag::TURB_INLET_DEFAULT {
                    Level::Note
                } else {
                    Level::Warning
                };
                push(level, code, Some(name), message);
            }
        }
    }
    // rough walls: k_s beyond the first cell centre (BC027)
    if !laminar {
        for p in &mesh.boundaries {
            let Some(bc) = bcs.get(&p.name) else { continue };
            let ks = bc.roughness();
            if !matches!(bc, BoundaryCondition::Wall { .. }) || ks <= 0.0 {
                continue;
            }
            let (mut total, mut beyond) = (0.0, 0.0);
            for f in p.range() {
                let s = geo.face_area[f];
                let a = norm(s);
                let c = mesh.owner[f] as usize;
                let y = (0..3)
                    .map(|i| (geo.face_centroid[f][i] - geo.cell_centroid[c][i]) * s[i] / a)
                    .sum::<f64>()
                    .abs();
                total += a;
                if ks > y {
                    beyond += a;
                }
            }
            if total > 0.0 && beyond > 0.1 * total {
                push(
                    Level::Warning,
                    diag::BC_ROUGHNESS,
                    Some(&p.name),
                    format!(
                        "roughness = {ks} m exceeds the first-cell distance on {:.0} % of \
                         the wall's area: the roughness elements reach beyond the first \
                         cell centre (coarsen the wall cells or check k_s)",
                        100.0 * beyond / total
                    ),
                );
            }
        }
    }
    // wall velocity normal to the wall, velocity inlets with outflow,
    // non-planar symmetry planes
    for p in &mesh.boundaries {
        let Some(bc) = bcs.get(&p.name) else { continue };
        match bc {
            BoundaryCondition::Wall { .. }
                if p.range().any(|f| {
                    norm(bc.wall_velocity(at.of(&p.name), geo.face_centroid[f], at.t)) > 0.0
                }) =>
            {
                let mut largest = 0.0f64;
                let worst = p
                    .range()
                    .map(|f| {
                        let s = geo.face_area[f];
                        let x = geo.face_centroid[f];
                        let v = bc.wall_velocity(at.of(&p.name), x, at.t);
                        largest = largest.max(norm(v));
                        // the normal velocity of a wall of a turning zone is
                        // that of the zone
                        let zone = Turning::velocity(at.turning, mesh, mesh.owner[f] as usize, x);
                        dot(sub(v, zone), s).abs() / norm(s)
                    })
                    .fold(0.0f64, f64::max);
                if worst > 1e-9 * largest {
                    push(
                        Level::Warning,
                        diag::BC_WALL_NORMAL_VELOCITY,
                        Some(&p.name),
                        format!(
                            "the wall velocity has a component of up to {worst:.3e} m/s normal to \
                             the wall; only the tangential part is used (a wall lets no mass through; \
                             a wall of a turning zone moves normal to itself with the zone)"
                        ),
                    );
                }
            }
            BoundaryCondition::VelocityInlet {
                profile: Some(InletProfile::Table),
                table_coordinate: Some(c),
                table_origin,
                table_direction,
                table_data: Some(t),
                ..
            } => {
                let outside = p
                    .range()
                    .filter(|&f| {
                        table_position(geo.face_centroid[f], *c, *table_origin, *table_direction)
                            .is_some_and(|x| t.at(x).1)
                    })
                    .count();
                if outside > 0 {
                    push(
                        Level::Warning,
                        diag::BC_TABLE_RANGE,
                        Some(&p.name),
                        format!(
                            "{outside} of {} faces lie outside the coordinate range [{}, {}] of \
                             the table {} and take its end values",
                            p.len,
                            t.s[0],
                            t.s[t.s.len() - 1],
                            t.path.display()
                        ),
                    );
                }
            }
            _ => {}
        }
        match bc {
            BoundaryCondition::VelocityInlet { profile, .. } => {
                let ext = boundary_extent(mesh, p);
                let out_faces = p
                    .range()
                    .filter(|&f| {
                        inlet_velocity(bc, at.of(&p.name), geo, f, at, &ext)
                            .is_some_and(|u| dot(u, geo.face_area[f]) > 0.0)
                    })
                    .count();
                if out_faces > 0 && *profile == Some(InletProfile::Analytic) {
                    push(
                        Level::Note,
                        diag::BC_EXACT_OUTFLOW,
                        Some(&p.name),
                        format!(
                            "the exact solution leaves the domain on {out_faces} of {} faces \
                             (the velocity is prescribed there as well)",
                            p.len
                        ),
                    );
                } else if out_faces > 0 {
                    push(
                        Level::Warning,
                        diag::BC_INLET_OUTFLOW,
                        Some(&p.name),
                        format!(
                            "the velocity points out of the domain on {out_faces} of {} faces; \
                             for an outlet use static-pressure-outlet, or reverse the velocity",
                            p.len
                        ),
                    );
                }
            }
            BoundaryCondition::Symmetry if p.len > 0 => {
                let n0 = {
                    let s = geo.face_area[p.start];
                    scale(s, 1.0 / norm(s))
                };
                let worst = p
                    .range()
                    .map(|f| {
                        let s = geo.face_area[f];
                        let n = scale(s, 1.0 / norm(s));
                        norm(sub(n, n0))
                    })
                    .fold(0.0f64, f64::max);
                if worst > 1e-6 {
                    push(
                        Level::Warning,
                        diag::BC_SYMMETRY_NOT_PLANAR,
                        Some(&p.name),
                        format!(
                            "the symmetry boundary is not planar (normals differ by up to {:.2}°); \
                             the condition is applied face by face",
                            2.0 * (worst / 2.0).asin().to_degrees()
                        ),
                    );
                }
            }
            _ => {}
        }
    }
    // the wedge planes of an axisymmetric case: two planes that meet in
    // the axis
    let wedges: Vec<&foehncraft_mesh::Boundary> = mesh
        .boundaries
        .iter()
        .filter(|p| matches!(bcs.get(&p.name), Some(BoundaryCondition::Wedge)))
        .collect();
    if !wedges.is_empty() {
        let planes: Vec<Option<([f64; 3], [f64; 3])>> = wedges
            .iter()
            .map(|p| boundary_plane(geo, p, 1e-6))
            .collect();
        if wedges.len() != 2 {
            push(
                Level::Error,
                diag::BC_WEDGE,
                Some(&wedges[0].name),
                format!(
                    "{} wedge boundary(es); an axisymmetric wedge has exactly two side planes \
                     (the other boundaries are walls, inlets, outlets or symmetry planes)",
                    wedges.len()
                ),
            );
        } else if let [Some((n1, _)), Some((n2, _))] = planes[..] {
            let axis = cross(n1, n2);
            if norm(axis) < 1e-9 {
                push(
                    Level::Error,
                    diag::BC_WEDGE,
                    Some(&wedges[0].name),
                    format!(
                        "the wedge planes `{}` and `{}` are parallel: a wedge mesh is the 2D mesh \
                         revolved by a small angle about the axis (for a plane 2D case use \
                         symmetry)",
                        wedges[0].name, wedges[1].name
                    ),
                );
            } else {
                let angle = (-dot(n1, n2)).clamp(-1.0, 1.0).acos().to_degrees();
                if angle > 10.0 {
                    push(
                        Level::Warning,
                        diag::BC_WEDGE_ANGLE,
                        Some(&wedges[0].name),
                        format!(
                            "the wedge is {angle:.2}° wide; 1° to 5° keeps the error of the \
                             planar faces small"
                        ),
                    );
                }
            }
        } else {
            for (p, pl) in wedges.iter().zip(&planes) {
                if pl.is_none() {
                    push(
                        Level::Error,
                        diag::BC_WEDGE,
                        Some(&p.name),
                        "the wedge boundary is not planar; each side plane of a wedge mesh is one \
                         plane through the axis"
                            .into(),
                    );
                }
            }
        }
    }
    // the domain as a whole: pressure level, mass balance, flow paths
    let convective: Vec<&str> = present
        .iter()
        .filter(|(_, b)| matches!(b, BoundaryCondition::ConvectiveOutlet { .. }))
        .map(|(n, _)| *n)
        .collect();
    let pressure: Vec<&str> = present
        .iter()
        .filter(|(_, b)| b.is_pressure())
        .map(|(n, _)| *n)
        .collect();
    let outlets: Vec<&str> = present
        .iter()
        .filter(|(_, b)| b.is_pressure_outlet())
        .map(|(n, _)| *n)
        .collect();
    // the scale of the balance is the flow through the faces, so that a
    // boundary that the flow enters and leaves (an exact solution) does not
    // make round-off look like an imbalance
    let (q, scale_q) = prescribed_outflow(mesh, geo, bcs, rho.into(), at);
    let net: f64 = q.values().exact_sum();
    let inflow: Vec<&str> = q
        .iter()
        .filter(|&(_, &v)| v < 0.0)
        .map(|(n, _)| n.as_str())
        .collect();
    // the convection velocity of the outlets without one: the mean outflow
    // speed that carries the prescribed inflow
    for &(name, bc) in &present {
        if let BoundaryCondition::ConvectiveOutlet {
            convection_velocity: None,
        } = bc
            && net >= -1e-12 * scale_q.max(f64::MIN_POSITIVE)
        {
            push(
                Level::Error,
                diag::BC_CONVECTIVE_VELOCITY,
                Some(name),
                "no prescribed inflow gives the mean outflow speed: set `convection_velocity` \
                 (m/s), for example the bulk velocity"
                    .into(),
            );
        }
    }
    if ctx.potential_flow && !convective.is_empty() && pressure.is_empty() {
        push(
            Level::Error,
            diag::BC_POTENTIAL_CONVECTIVE,
            Some(convective[0]),
            "the potential-flow start needs a boundary that takes up the inflow: with \
             convective outlets and no pressure boundary start without `potential_flow`"
                .into(),
        );
    }
    if pressure.is_empty() && !convective.is_empty() {
        push(
            Level::Note,
            diag::BC_REFERENCE_CELL,
            None,
            format!(
                "the convective outlet(s) {} carry the prescribed inflow (their velocity is \
                 balanced with it in every iteration); the pressure level is fixed in the \
                 reference cell ([pressure] reference_point, reference_value)",
                convective
                    .iter()
                    .map(|n| format!("`{n}`"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        );
    } else if pressure.is_empty() {
        let wall_net: f64 = q
            .iter()
            .filter(|(n, _)| {
                matches!(
                    bcs.get(n.as_str()),
                    Some(BoundaryCondition::Wall {
                        mass_flux: Some(_),
                        ..
                    })
                )
            })
            .map(|(_, v)| *v)
            .exact_sum();
        let walls: Vec<&str> = q
            .keys()
            .filter(|n| {
                matches!(
                    bcs.get(n.as_str()),
                    Some(BoundaryCondition::Wall {
                        mass_flux: Some(_),
                        ..
                    })
                )
            })
            .map(String::as_str)
            .collect();
        if scale_q > 0.0 && net.abs() > 1e-9 * scale_q && wall_net.abs() > 1e-9 * scale_q {
            push(
                Level::Error,
                diag::BC_WALL_MASS_CLOSED,
                Some(walls[0]),
                format!(
                    "the wall mass flux of {} {} a net {:.6e} kg/s {} a closed volume, \
                     and the density is constant: give the domain a pressure boundary \
                     (static-pressure-outlet), or balance the wall mass flux with the \
                     other prescribed boundaries",
                    walls
                        .iter()
                        .map(|n| format!("`{n}`"))
                        .collect::<Vec<_>>()
                        .join(", "),
                    if net < 0.0 { "adds" } else { "removes" },
                    net.abs(),
                    if net < 0.0 { "to" } else { "from" }
                ),
            );
        } else if scale_q > 0.0 && net.abs() > 1e-9 * scale_q {
            push(
                Level::Error,
                diag::BC_UNBALANCED,
                None,
                format!(
                    "the prescribed boundaries let a net {:.6e} kg/s {} the domain and there is no \
                     pressure boundary to balance it: make one boundary a static-pressure-outlet \
                     (or a total-pressure-outlet), or balance the inflow and outflow",
                    net.abs(),
                    if net < 0.0 { "into" } else { "out of" }
                ),
            );
        } else {
            push(
                Level::Note,
                diag::BC_REFERENCE_CELL,
                None,
                "no pressure boundary: the pressure level is fixed in the reference cell \
                 ([pressure] reference_point, reference_value)"
                    .into(),
            );
        }
    } else if outlets.is_empty() && convective.is_empty() && !inflow.is_empty() {
        push(
            Level::Error,
            diag::BC_OVER_SPECIFIED,
            None,
            format!(
                "the flow entering through {} can only leave through the pressure inlet(s) {}: \
                 the flow rate is over-specified on this flow path; change one of {} to \
                 static-pressure-outlet",
                inflow
                    .iter()
                    .map(|n| format!("`{n}`"))
                    .collect::<Vec<_>>()
                    .join(", "),
                pressure
                    .iter()
                    .map(|n| format!("`{n}`"))
                    .collect::<Vec<_>>()
                    .join(", "),
                pressure
                    .iter()
                    .map(|n| format!("`{n}`"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        );
    } else if inflow.is_empty() {
        let inlets: Vec<&str> = present
            .iter()
            .filter(|(_, b)| b.is_pressure() && !b.is_pressure_outlet())
            .map(|(n, _)| *n)
            .collect();
        let moving = present.iter().any(|(_, b)| {
            matches!(b, BoundaryCondition::Wall { velocity, .. }
                if velocity.iter().any(|v| v.number() != Some(0.0)))
        });
        if inlets.is_empty() && outlets.len() == 1 && !moving {
            push(
                Level::Warning,
                diag::BC_FLUID_AT_REST,
                Some(outlets[0]),
                "the only open boundary is this outlet and no boundary drives a flow; unless a \
                 body force acts, the result is the fluid at rest (add an inlet or a moving wall)"
                    .into(),
            );
        } else {
            push(
                Level::Note,
                diag::BC_PRESSURE_DRIVEN,
                None,
                "no prescribed inflow: the flow rate follows from the pressure boundaries".into(),
            );
        }
    }
    out.sort_by_key(|f| std::cmp::Reverse(f.level));
    out
}

/// Per boundary face arrays of the conditions, in boundary-face order
/// (face index − number of interior faces).
#[derive(Debug, Clone, Default)]
pub struct BoundaryArrays {
    /// Momentum kind (`MOM_*`).
    pub mom_kind: Vec<u32>,
    /// Flux kind (`FLUX_*`).
    pub flux_kind: Vec<u32>,
    /// Update rule of the face velocity (`UB_*`).
    pub ub_rule: Vec<u32>,
    /// Update rule of the face pressure (`PB_*`).
    pub pb_rule: Vec<u32>,
    /// Prescribed face velocity (and the start value of the others).
    pub ub: [Vec<f64>; 3],
    /// Prescribed pressure or total pressure.
    pub pset: Vec<f64>,
    /// Boundary of each boundary face.
    pub boundary: Vec<u32>,
    /// 1 where the velocity has a zero normal gradient in the least-squares
    /// gradient (outlets), else 0 (the face value enters).
    pub u_zero_gradient: Vec<f64>,
    /// 1 where the pressure has a zero normal gradient in the least-squares
    /// gradient (walls, symmetry planes), else 0.
    pub p_zero_gradient: Vec<f64>,
    /// Stencil point of the velocity gradient per boundary: `ZeroGradient` for
    /// the normal projection of the face centroid (outlets, symmetry
    /// planes), `FixedValue` for the face centroid.
    pub u_points: BTreeMap<String, foehncraft_fv::problem::BoundaryKind>,
    /// Stencil point of the pressure gradient per boundary (`ZeroGradient` at
    /// walls and symmetry planes).
    pub p_points: BTreeMap<String, foehncraft_fv::problem::BoundaryKind>,
    /// Convection velocity U_c (m/s) of the faces of convective outlets
    /// that give one, 0 for the default (and on all other faces).
    pub convection: Vec<f64>,
    /// The driving of the periodic pairs (D29), one per pair of the mesh.
    pub drives: Vec<PeriodicDrive>,
}

/// The driving of one periodic pair on the mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct PeriodicDrive {
    /// The driving.
    pub drive: Drive,
    /// Unit vector along the translation, from side `a` to side `b`.
    pub direction: [f64; 3],
    /// The faces of side `b` (global face indices): the flow from `a`
    /// towards `b` leaves the domain through them.
    pub b_faces: std::ops::Range<usize>,
}

impl BoundaryArrays {
    /// True when some face fixes the pressure.
    pub fn has_pressure(&self) -> bool {
        self.flux_kind
            .iter()
            .any(|&k| k == FLUX_PRESSURE || k == FLUX_PRESSURE_NO_BACKFLOW)
    }
}

/// The per-face arrays of the conditions (after [`check`] found no error).
pub fn arrays<'a>(
    mesh: &Mesh,
    geo: &Geometry,
    bcs: &BTreeMap<String, BoundaryCondition>,
    rho: impl Into<Density<'a>>,
    at: At,
) -> Result<BoundaryArrays> {
    let rho: Density = rho.into();
    let ni = mesh.n_interior_faces();
    let nb = mesh.n_boundary_faces();
    let mut a = BoundaryArrays {
        mom_kind: vec![MOM_FIXED; nb],
        flux_kind: vec![FLUX_FIXED; nb],
        ub_rule: vec![UB_FIXED; nb],
        pb_rule: vec![PB_EXTRAPOLATE; nb],
        ub: [vec![0.0; nb], vec![0.0; nb], vec![0.0; nb]],
        pset: vec![0.0; nb],
        boundary: vec![0; nb],
        u_zero_gradient: vec![0.0; nb],
        p_zero_gradient: vec![0.0; nb],
        u_points: BTreeMap::new(),
        p_points: BTreeMap::new(),
        convection: vec![0.0; nb],
        drives: Vec::new(),
    };
    if let Some(per) = &mesh.periodic {
        let cases = periodic_pairs(bcs)?;
        for pair in &per.pairs {
            let drive = cases
                .iter()
                .find(|c| c.spec.a == pair.a && c.spec.b == pair.b)
                .map_or(Drive::None, |c| c.drive);
            let t = pair.translation;
            let l = norm(t);
            a.drives.push(PeriodicDrive {
                drive,
                direction: [t[0] / l, t[1] / l, t[2] / l],
                b_faces: pair.b_faces.clone(),
            });
        }
    }
    use foehncraft_fv::problem::BoundaryKind as K;
    for (pi, p) in mesh.boundaries.iter().enumerate() {
        let bc = bcs
            .get(&p.name)
            .ok_or_else(|| anyhow::anyhow!("no boundary condition for `{}`", p.name))?;
        let area: f64 = p.range().map(|f| norm(geo.face_area[f])).exact_sum();
        let ext = boundary_extent(mesh, p);
        let outlet = bc.is_pressure_outlet();
        let wall_like = bc.is_symmetry_like()
            || matches!(
                bc,
                BoundaryCondition::Wall { .. } | BoundaryCondition::ConvectiveOutlet { .. }
            );
        let symmetry = bc.is_symmetry_like();
        a.u_points.insert(
            p.name.clone(),
            if outlet || symmetry {
                K::ZeroGradient
            } else {
                K::FixedValue
            },
        );
        a.p_points.insert(
            p.name.clone(),
            if wall_like {
                K::ZeroGradient
            } else {
                K::FixedValue
            },
        );
        for f in p.range() {
            let fb = f - ni;
            a.boundary[fb] = pi as u32;
            a.u_zero_gradient[fb] = if outlet { 1.0 } else { 0.0 };
            a.p_zero_gradient[fb] = if wall_like { 1.0 } else { 0.0 };
            let s = geo.face_area[f];
            let n = scale(s, 1.0 / norm(s));
            let mut set_ub = |v: [f64; 3]| {
                for (ub, vi) in a.ub.iter_mut().zip(v) {
                    ub[fb] = vi;
                }
            };
            let vals = at.of(&p.name);
            let xf = geo.face_centroid[f];
            match bc {
                BoundaryCondition::Wall { mass_flux, .. } => {
                    let velocity = bc.wall_velocity(vals, xf, at.t);
                    // the tangential part of the given velocity; normal to
                    // the wall it moves with its zone if that turns (a
                    // blade, D55/D182: absolute flux ρ S·V, relative flux 0)
                    let zone = Turning::velocity(at.turning, mesh, mesh.owner[f] as usize, xf);
                    let wall = add(
                        sub(velocity, scale(n, dot(velocity, n))),
                        scale(n, dot(zone, n)),
                    );
                    if mass_flux.is_some() {
                        // the normal velocity that carries the mass flux into
                        // the domain (D124 W3)
                        let un = bc.wall_mass_flux(vals, xf, at.t) / rho.at(f);
                        set_ub(sub(wall, scale(n, un)));
                    } else {
                        set_ub(wall);
                    }
                }
                BoundaryCondition::VelocityInlet { .. } => {
                    set_ub(inlet_velocity(bc, vals, geo, f, at, &ext).unwrap_or([0.0; 3]));
                }
                BoundaryCondition::MassFlowInlet { .. } => {
                    let mass_flow = bc.mass_flow_at(vals, at.t);
                    set_ub(scale(n, -mass_flow / (rho.at(f) * area)));
                }
                BoundaryCondition::StaticPressureInlet { .. } => {
                    a.ub_rule[fb] = UB_FROM_FLUX;
                    a.flux_kind[fb] = FLUX_PRESSURE;
                    a.pb_rule[fb] = PB_FIXED;
                    a.pset[fb] = bc.pressure_at(vals, xf, at.t);
                }
                BoundaryCondition::TotalPressureInlet { .. } => {
                    a.ub_rule[fb] = UB_FROM_FLUX;
                    a.flux_kind[fb] = FLUX_PRESSURE;
                    a.pb_rule[fb] = PB_TOTAL;
                    a.pset[fb] = bc.pressure_at(vals, xf, at.t);
                }
                BoundaryCondition::StaticPressureOutlet { backflow, .. } => {
                    a.mom_kind[fb] = MOM_OUTFLOW;
                    let (ub, flux, pb) = match backflow {
                        Backflow::TotalPressure => (UB_OUTLET, FLUX_PRESSURE, PB_OUTLET),
                        Backflow::ZeroGradient => (UB_EXTRAPOLATE, FLUX_PRESSURE, PB_FIXED),
                        Backflow::Suppress => {
                            (UB_NO_BACKFLOW, FLUX_PRESSURE_NO_BACKFLOW, PB_NO_BACKFLOW)
                        }
                    };
                    a.ub_rule[fb] = ub;
                    a.flux_kind[fb] = flux;
                    a.pb_rule[fb] = pb;
                    a.pset[fb] = bc.pressure_at(vals, xf, at.t);
                }
                BoundaryCondition::TotalPressureOutlet { backflow, .. } => {
                    a.mom_kind[fb] = MOM_OUTFLOW;
                    let (ub, flux, pb) = match backflow {
                        Backflow::TotalPressure => (UB_OUTLET, FLUX_PRESSURE, PB_TOTAL),
                        Backflow::ZeroGradient => (UB_EXTRAPOLATE, FLUX_PRESSURE, PB_TOTAL),
                        Backflow::Suppress => (
                            UB_NO_BACKFLOW,
                            FLUX_PRESSURE_NO_BACKFLOW,
                            PB_TOTAL_NO_BACKFLOW,
                        ),
                    };
                    a.ub_rule[fb] = ub;
                    a.flux_kind[fb] = flux;
                    a.pb_rule[fb] = pb;
                    a.pset[fb] = bc.pressure_at(vals, xf, at.t);
                }
                BoundaryCondition::ConvectiveOutlet {
                    convection_velocity,
                } => {
                    a.ub_rule[fb] = UB_CONVECTIVE;
                    a.convection[fb] = convection_velocity.unwrap_or(0.0);
                }
                BoundaryCondition::Symmetry
                | BoundaryCondition::SlipWall
                | BoundaryCondition::Wedge => {
                    a.mom_kind[fb] = MOM_SYMMETRY;
                    a.ub_rule[fb] = UB_SYMMETRY;
                }
                BoundaryCondition::Periodic { .. } => {
                    anyhow::bail!("periodic boundary `{}` was not paired", p.name)
                }
            }
            // a static pressure with gravity: the reduced pressure of the
            // face, p′_b = p_b − ρ_ref g·(x_f − x_ref) (D153)
            if let (Some(h), Some(crate::buoyancy::PressureKind::Static)) =
                (at.hydrostatic, bc.pressure_kind())
            {
                a.pset[fb] -= h.head(xf);
            }
        }
    }
    Ok(a)
}

/// Guidance: the sensible combinations and the cases each one suits (also
/// in `docs/cases.md`).
pub const GUIDANCE: &str = "\
| Combination | Suits | Pressure level |
|---|---|---|
| velocity-inlet + static-pressure-outlet (+ walls, symmetry) | internal flows with a known flow rate: channels, pipes, ducts, the Poiseuille case | outlet pressure |
| mass-flow-inlet + static-pressure-outlet | as above when the mass flow is known rather than a velocity profile | outlet pressure |
| total-pressure-inlet + static-pressure-outlet | pressure-driven flows: the flow rate is a result (fans, valves, plenum to plenum) | both |
| static-pressure-inlet + static-pressure-outlet | pressure-driven flows with a known static pressure drop | both |
| static-pressure-inlet + total-pressure-outlet | flows that exhaust into a large volume at a known total pressure | both |
| walls only, one wall moving (+ symmetry) | enclosed flows: the lid-driven cavity | reference cell |
| velocity-inlet + velocity-inlet + static-pressure-outlet | mixing of two streams: the T-junction | outlet pressure |
| velocity-inlet + convective-outlet (+ walls, slip walls), transient | wakes and vortices leaving a short domain without reflection: the cylinder wake | reference cell |
| velocity-inlet + slip-wall + static-pressure-outlet | a body in a free stream with the far field as slip walls | outlet pressure |
| wedge + wedge + velocity-inlet + static-pressure-outlet (+ walls) | axisymmetric flows on a wedge mesh: pipes, nozzles, sudden expansions | outlet pressure |
| periodic pair (`mass_flow` or `pressure_gradient`) + walls (+ symmetry, wedges, further periodic pairs) | fully developed flow in one module: channels, pipes, ribbed passages, heat transfer with the heat-flux jump | reference cell |
";
