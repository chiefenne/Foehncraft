//! The `foehn flow` case: steady incompressible laminar flow on one mesh,
//! or an order-of-accuracy sweep over refinement levels.
//!
//! Case file (TOML, D39). Every key is optional unless noted, SI units
//! (D44), unknown keys and values are errors with a "did you mean" hint:
//!
//! ```toml
//! [mesh]                         # as in `foehn mesh` (file or generate)
//! [mesh.generate]
//! type = "hex"
//! dimension = 2
//! cells = [64, 64, 1]
//! [fluid]
//! density = 1.0                  # ρ (kg/m³)
//! viscosity = 0.01               # μ (Pa s), or an expression of x, y, z (D53)
//! # material = "water"           # a material instead (D94-D102, crate::materials)
//! # temperature = 293.15         # K, where the material is evaluated
//! [flow]
//! convection = true              # false: Stokes flow
//! body_force = [0.0, 0.0, 0.0]   # N/m³, numbers or expressions of x, y, z, t
//! [constants]                    # named numbers of the expressions (D53)
//! u_max = 1.5
//! [tables.ramp]                  # a table, called as ramp(t) (crate::expr)
//! x = [0.0, 1.0]
//! y = [0.0, 1.0]
//! [initial]
//! velocity = [0.0, 0.0, 0.0]
//! pressure = 0.0
//! [pressure]                     # only used when no boundary fixes the pressure
//! reference_point = [0.5, 0.5, 0.0]
//! reference_value = 0.0
//! [boundary.ymax]                # one table per boundary (required), crate::bc
//! type = "wall"
//! velocity = [1.0, 0.0, 0.0]
//! [solver]                       # crate::simplec::Settings
//! relaxation_velocity = 0.9
//! [analytic]                     # optional exact solution, crate::analytic
//! type = "poiseuille"
//! [verification]
//! levels = [0, 1, 2]             # refine levels of [mesh.generate]: a sweep
//! expected_order = [1.8, 2.2]
//! max_error = 1e-9               # largest relative L2 velocity error
//! reference = "ghia-re100"       # centre-line data of the cavity
//! reference_tolerance = 0.02
//! [output]
//! vtkhdf = true
//! history_every = 1              # residual history: every N iterations
//! wall_csv = true                # one CSV of wall quantities per wall
//! derived = ["vorticity", "q-criterion"]  # derived fields (D72), crate::derived
//! [output.reference]             # reference values of pressure-coefficient
//! velocity = 1.0
//! [[output.lines]]               # velocity and pressure along a line
//! name = "vertical"
//! from = [0.5, 0.0, 0.0]
//! to = [0.5, 1.0, 0.0]
//! ```
//!
//! Output (D43): `flow.log` (header with what is written and how large),
//! `residuals.csv`, `boundaries.csv` (mass flow, mean pressure and force per
//! boundary), `wall_<boundary>.csv`, `line_<name>.csv`, `fields.vtkhdf`,
//! `case.resolved.toml` (the case with every default filled in), `flow.json`
//! and `flow.md` (the report), SVG charts.

use crate::analytic::{Analytic, Benchmark, Reference};
use crate::balance::{self, BalanceSummary};
use crate::bc::{self, BoundaryCondition, Level};
use crate::control::{self, Control};
use crate::diag;
use crate::divergence;
use crate::engine::Steer;
use crate::expr::{Env, Scalar, Value, Vars};
use crate::plot::{self, Series};
use crate::restart::{self, Restart};
use crate::simplec::{self, Flow, Physics, Residuals, Settings, Solution, Timings};
use crate::transient::TimeScheme;
use anyhow::{Context, Result, bail};
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_mesh::case::{MeshSection, build_mesh};
use foehncraft_mesh::geometry::{cross, dot, norm, sub};
use foehncraft_mesh::locate::Locator;
use foehncraft_mesh::renumber::renumber;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// `[fluid]`: the fluid's properties, as numbers (the shorthand of one
/// constant material, D100) or from a material (D94–D102).
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct FluidSection {
    /// Density ρ (kg/m³); default 1, or the material's. One number: a
    /// density that varies comes with the energy equation (L1).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub density: Option<Scalar>,
    /// Dynamic viscosity μ (Pa s): a number or an expression of x, y, z
    /// (D53); default 0.001, or the material's.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub viscosity: Option<Scalar>,
    /// A material of `[materials]` or of the database (`foehn materials`);
    /// `density` and `viscosity` given here replace its values.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub material: Option<String>,
    /// The temperature (K) at which the material's properties are
    /// evaluated (isothermal flow); default 293.15 K.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature: Option<f64>,
    /// The operating pressure p_op (Pa): the absolute level of the gauge
    /// pressure the solver works with, and the pressure of the ideal-gas
    /// density. With absolute pressure values (the default) each fluid part
    /// takes it from its own pressure boundaries and it is rarely given;
    /// with gauge values the default is 101 325 Pa.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub operating_pressure: Option<f64>,
    /// What the pressure values of the case file are (D239): `"absolute"`
    /// static pressures (the default from case format 0.7) or `"gauge"`
    /// values relative to the operating pressure, p = p_op + p_gauge (the
    /// meaning of files of format 0.6 and older). It applies to every
    /// pressure value of the file: boundary `pressure` and
    /// `total_pressure` (numbers and expressions), `[initial] pressure`
    /// and `[pressure] reference_value`; the outputs (fields, probes,
    /// boundary values) are written in the same convention.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pressure_values: Option<crate::pressure_values::PressureValues>,
    /// Specific heat c_p (J/(kg K)) for the energy equation: a number or
    /// an expression of T; replaces the material's.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub specific_heat: Option<Scalar>,
    /// Thermal conductivity k (W/(m K)) for the energy equation: a number
    /// or an expression of T; replaces the material's.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub conductivity: Option<Scalar>,
}

impl FluidSection {
    /// The shorthand of one constant fluid.
    pub fn constant(density: f64, viscosity: f64) -> Self {
        Self {
            density: Some(Scalar::Number(density)),
            viscosity: Some(Scalar::Number(viscosity)),
            ..Self::default()
        }
    }

    /// The conditions at which the properties are evaluated.
    pub fn conditions(&self) -> crate::materials::Conditions {
        crate::materials::Conditions {
            temperature: self
                .temperature
                .unwrap_or(crate::materials::DEFAULT_TEMPERATURE),
            pressure: self
                .operating_pressure
                .unwrap_or(crate::materials::DEFAULT_PRESSURE),
        }
    }
}

/// `[flow]`
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct FlowSection {
    /// Convection on (false: Stokes flow).
    pub convection: bool,
    /// Body force per unit volume (N/m³): three numbers, or expressions of
    /// x, y, z and, in a transient case, t (D53).
    pub body_force: crate::expr::Vector,
    /// Angular velocity Ω (rad/s) of a frame in which the whole domain is
    /// solved (stage 6, 9.8): the velocity is relative to that frame and
    /// the momentum equations take the Coriolis force −2ρ Ω × u; the
    /// centrifugal force −ρ Ω × (Ω × r) is a gradient and stays in the
    /// pressure (the solved pressure is p − ½ρ|Ω × r|²). Not with `[[mrf]]`
    /// or `[[sliding]]`.
    pub frame_rotation: [f64; 3],
}

impl Default for FlowSection {
    fn default() -> Self {
        Self {
            convection: true,
            body_force: crate::expr::vector([0.0; 3]),
            frame_rotation: [0.0; 3],
        }
    }
}

/// `[initial]`
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct InitialSection {
    /// Start velocity (m/s).
    pub velocity: [f64; 3],
    /// Start pressure (Pa).
    pub pressure: f64,
    /// Start from potential flow (D52) instead of `velocity`: the
    /// irrotational velocity that carries the inflow through the domain to
    /// the pressure boundaries, from one Laplace solve with the pressure
    /// multigrid.
    pub potential_flow: bool,
    /// Transient: start from the velocity and pressure of `[analytic]`
    /// (at `start_time`) instead of `velocity` and `pressure`, to see
    /// whether a numerical method keeps an exact solution.
    pub from_analytic: bool,
    /// Start turbulence intensity (stage 6, with `viscosity_ratio`; default:
    /// the values of the first inlet by name, else 0.01).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub turbulence_intensity: Option<f64>,
    /// Start ratio ν_t/ν (default 1 with `turbulence_intensity`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub viscosity_ratio: Option<f64>,
    /// Start k (m²/s²): a number or an expression of x, y, z (with `omega`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub k: Option<crate::expr::Scalar>,
    /// Start ω (1/s): a number or an expression of x, y, z (with `k`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub omega: Option<crate::expr::Scalar>,
    /// Start ε (m²/s³): a number or an expression of x, y, z (with `k`;
    /// the k-ε family).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub epsilon: Option<crate::expr::Scalar>,
    /// Start ζ = v²/k (k-ε-ζ-f; φ of the elliptic-blending k-ε; default
    /// 2/3, isotropy; in (0, 2]).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub zeta: Option<f64>,
}

/// `[pressure]`
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct PressureSection {
    /// Point whose cell holds the reference pressure (default: the cell
    /// nearest to the lower corner of the bounding box).
    pub reference_point: Option<[f64; 3]>,
    /// Pressure in the reference cell (Pa).
    pub reference_value: f64,
}

/// `[run]`
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct RunSection {
    /// CPU threads (default: all cores; meshes below 20 000 cells run on
    /// one).
    pub threads: Option<usize>,
    /// Where the flow loop runs: `cpu` or `cuda` (a build with the `cuda`
    /// feature; the results are bitwise those of `cpu`, D31).
    pub backend: String,
    /// PTX module of the `cuda` backend.
    pub ptx: PathBuf,
    /// Restart file to continue from (D42). A transient case may also
    /// start from the restart file of a steady run: its fields are the
    /// start fields at `time.start_time`, with no older time levels.
    pub restart: Option<PathBuf>,
    /// Write `restart.h5` every N iterations or time steps (0: never).
    pub restart_every: usize,
    /// Also write a numbered `restart_<N>.h5` at these iterations or steps.
    pub restart_at: Vec<usize>,
    /// Transient: write `restart.h5` every this much simulated time (s), at
    /// the multiples of it; the time step is shortened to land on them
    /// (D66).
    pub restart_interval: Option<f64>,
    /// Transient: also write a numbered `restart_<N>.h5` at these simulated
    /// times (s), landed on exactly (D66).
    pub restart_times: Vec<f64>,
    /// Write `restart.h5` at the end of the run.
    pub restart_final: bool,
    /// Re-read the case file every N iterations or time steps (0: only on a
    /// `.read_toml_N` control file), D40.
    pub reread_every: usize,
    /// Limits of the run (D83).
    pub limits: Limits,
    /// Also solve the flow on these numbers of ranks (MPI step W5): the
    /// ranks run as threads of this process on the parts of a recursive
    /// coordinate bisection, and the run fails unless every rank's cells
    /// are bitwise those of one rank. Steady and transient, SIMPLEC, PISO
    /// and the coupled solver, laminar, k-ω SST and k-ε, passive scalars
    /// and the energy equation of one fluid; not yet with gravity, turning
    /// zones, regions, a schedule or a restart.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub ranks: Vec<usize>,
}

/// `[run.limits]` (D83): when one is reached the run writes a restart file
/// and the final output and ends with exit code 3 (stopped, D64).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Limits {
    /// Largest wall time of the run (s), setup included.
    pub max_wall_time: Option<f64>,
    /// Largest number of time steps of the run (transient; a continued run
    /// counts from its restart).
    pub max_steps: Option<usize>,
    /// Largest number of outer iterations of the run (steady).
    pub max_iterations: Option<usize>,
    /// Largest volume of the output directory (bytes).
    pub max_output: Option<f64>,
    /// Steady: the outer iteration at which a run that has not converged is
    /// diagnosed from its residual history (D198 (d); default 2000, 0
    /// switches it off). The diagnosis goes into the run record and the
    /// log.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub diagnose_at: Option<usize>,
    /// What a diagnosis that finds the run off track does: `record`
    /// (default) or `stop` (restart and final output written, exit code 3).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub diagnosis: Option<DiagnosisAction>,
    /// Cell-hours (cells times wall hours) that all runs of the study may
    /// spend together (D198 (c)); the study is the parent of the output
    /// directory, its ledger `study_budget.jsonl`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub study_cell_hours: Option<f64>,
    /// GPU-hours that all runs of the study may spend together.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub study_gpu_hours: Option<f64>,
    /// Bytes that all runs of the study may write together.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub study_output: Option<f64>,
}

/// `[run.limits] diagnosis` (D198 (d)).
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum DiagnosisAction {
    /// Write the diagnosis into the run record and the log; go on.
    #[default]
    Record,
    /// Stop a run that is not on its way to the tolerance within its
    /// budget.
    Stop,
}

impl Limits {
    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        if let Some(t) = self.max_wall_time {
            crate::need!(
                t > 0.0,
                CASE_RANGE,
                "[run.limits] max_wall_time must be > 0"
            );
        }
        if let Some(n) = self.max_steps {
            crate::need!(n >= 1, CASE_RANGE, "[run.limits] max_steps must be >= 1");
        }
        if let Some(n) = self.max_iterations {
            crate::need!(
                n >= 1,
                CASE_RANGE,
                "[run.limits] max_iterations must be >= 1"
            );
        }
        if let Some(b) = self.max_output {
            crate::need!(b > 0.0, CASE_RANGE, "[run.limits] max_output must be > 0");
        }
        for (k, v) in [
            ("study_cell_hours", self.study_cell_hours),
            ("study_gpu_hours", self.study_gpu_hours),
            ("study_output", self.study_output),
        ] {
            if let Some(b) = v {
                crate::need!(b > 0.0, CASE_RANGE, "[run.limits] {k} must be > 0");
            }
        }
        Ok(())
    }

    /// A diagnosis off track stops the run.
    pub fn stops_at_diagnosis(&self) -> bool {
        self.diagnosis == Some(DiagnosisAction::Stop)
    }

    /// A study budget is set.
    pub fn has_study_budget(&self) -> bool {
        self.study_cell_hours.is_some()
            || self.study_gpu_hours.is_some()
            || self.study_output.is_some()
    }

    /// True when no limit is set.
    pub fn is_empty(&self) -> bool {
        *self == Self::default()
    }
}

impl Default for RunSection {
    fn default() -> Self {
        Self {
            threads: None,
            backend: "cpu".into(),
            ptx: PathBuf::from("target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx"),
            restart: None,
            restart_every: 0,
            restart_at: Vec::new(),
            restart_interval: None,
            restart_times: Vec::new(),
            restart_final: true,
            reread_every: 0,
            limits: Limits::default(),
            ranks: Vec::new(),
        }
    }
}

/// `[time]`: a transient case (D28).
#[derive(Debug, Clone, Copy, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct TimeSection {
    /// Time integration.
    pub scheme: TimeScheme,
    /// Time step Δt (s); with `adaptive` the first one.
    pub time_step: f64,
    /// End time (s).
    pub end_time: Option<f64>,
    /// Number of time steps (instead of, or as well as, `end_time`).
    pub steps: Option<usize>,
    /// Time of the initial fields (s).
    pub start_time: f64,
    /// Adapt Δt to the largest Courant number `max_cfl`.
    pub adaptive: bool,
    /// Target of the largest Courant number with `adaptive`.
    pub max_cfl: f64,
    /// Largest ratio of two successive time steps with `adaptive`.
    pub max_growth: f64,
    /// Smallest time step with `adaptive` (s).
    pub min_time_step: f64,
    /// Largest time step with `adaptive` (s).
    pub max_time_step: Option<f64>,
    /// Relative tolerance of the local time-integration error with
    /// `adaptive` (D65): the step follows from the error estimate of BDF2
    /// or Crank–Nicolson, bounded by `max_cfl` and `max_growth`, and a
    /// step whose error exceeds the tolerance is repeated with a shorter
    /// one. Off when absent.
    pub error_tolerance: Option<f64>,
    /// Absolute tolerance of the local time-integration error (m/s), added
    /// to `error_tolerance` × the RMS velocity.
    pub error_absolute: f64,
}

impl Default for TimeSection {
    fn default() -> Self {
        Self {
            scheme: TimeScheme::Bdf2,
            time_step: 0.0,
            end_time: None,
            steps: None,
            start_time: 0.0,
            adaptive: false,
            max_cfl: 1.0,
            max_growth: 1.2,
            min_time_step: 0.0,
            max_time_step: None,
            error_tolerance: None,
            error_absolute: 1e-6,
        }
    }
}

impl TimeSection {
    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        crate::need!(
            self.time_step > 0.0,
            CASE_RANGE,
            "[time] time_step must be > 0 (the first step with `adaptive`)"
        );
        crate::need!(
            self.end_time.is_some() || self.steps.is_some(),
            CASE_MISSING_KEY,
            "[time] needs `end_time` or `steps`"
        );
        if let Some(e) = self.end_time {
            crate::need!(
                e > self.start_time,
                CASE_RANGE,
                "[time] end_time must be later than start_time"
            );
        }
        crate::need!(self.max_cfl > 0.0, CASE_RANGE, "[time] max_cfl must be > 0");
        crate::diag::check_fix(
            self.max_growth > 1.0 && self.max_growth < 1.0 + std::f64::consts::SQRT_2,
            crate::diag::CASE_RANGE,
            &["time.max_growth=1.2"],
            || {
                String::from(
                    "[time] max_growth must lie in (1, 2.414): BDF2 with variable steps is zero-stable \
             only for step ratios below 1 + √2",
                )
            },
        )?;
        crate::need!(
            self.min_time_step >= 0.0,
            CASE_RANGE,
            "[time] min_time_step must be >= 0"
        );
        if let Some(m) = self.max_time_step {
            crate::need!(
                m >= self.min_time_step && m > 0.0,
                CASE_RANGE,
                "[time] max_time_step must be > 0 and >= min_time_step"
            );
        }
        if let Some(e) = self.error_tolerance {
            crate::need!(
                e > 0.0 && e < 1.0,
                CASE_RANGE,
                "[time] error_tolerance must lie in (0, 1)"
            );
            crate::need!(
                self.adaptive,
                CASE_RANGE,
                "[time] error_tolerance needs adaptive = true"
            );
            crate::need!(
                self.scheme != TimeScheme::Bdf1,
                CASE_RANGE,
                "[time] error_tolerance needs scheme bdf2 or crank-nicolson (the estimate is that \
                 of a second-order method)"
            );
        }
        crate::need!(
            self.error_absolute >= 0.0,
            CASE_RANGE,
            "[time] error_absolute must be >= 0"
        );
        Ok(())
    }
}

/// A point probe, `[[output.probes]]`.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Probe {
    /// Name (columns `<name>_u` … in `probes.csv`).
    pub name: String,
    /// The point.
    pub point: [f64; 3],
}

/// A wall probe, `[[output.wall_probes]]`: the wall shear stress and y⁺ on
/// a wall at the point nearest to `point`, linear between the two nearest
/// face centroids of `boundary` (stage 6: the skin friction of the flat
/// plate at x = 0.97).
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WallProbe {
    /// Name (in the report and the uncertainty table).
    pub name: String,
    /// The wall boundary.
    pub boundary: String,
    /// The point (m).
    pub point: [f64; 3],
}

/// A sign change of the wall shear, `[[output.wall_crossings]]`: the first
/// point after `from` along `direction` where the component of the wall
/// shear stress along `direction` changes sign on a wall, linear between
/// the two face centroids (stage 6: the reattachment length of the
/// backward-facing step, with its uncertainty in a level sweep). The
/// faces are taken in the order of their centroids along `direction`, so
/// the wall is one row of faces: a 2D case (one cell thick), or `width`.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WallCrossing {
    /// Name (in the report and the uncertainty table).
    pub name: String,
    /// The wall boundary.
    pub boundary: String,
    /// `reattachment` (the shear along `direction` turns from negative to
    /// positive) or `separation` (from positive to negative).
    #[serde(default)]
    pub kind: CrossingKind,
    /// The point the search starts from (m); faces whose centroids lie
    /// before it along `direction` are left out. The position is reported
    /// as the distance from it along `direction`.
    pub from: [f64; 3],
    /// The direction of the search (normalised).
    #[serde(default = "x_axis")]
    pub direction: [f64; 3],
    /// Leaves out the faces closer than this distance (m) to `from` along
    /// `direction`; the position is still measured from `from` (default
    /// 0). Behind a step the corner vortices change the sign of the shear
    /// close to the step face, so the search for the main reattachment
    /// starts after them.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skip: Option<f64>,
    /// Only the faces whose centroids lie within this distance (m) of the
    /// line through `from` along `direction` (default: every face of the
    /// boundary).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub width: Option<f64>,
}

/// The kind of a wall crossing.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum CrossingKind {
    /// The shear along the direction turns from negative to positive.
    #[default]
    Reattachment,
    /// The shear along the direction turns from positive to negative.
    Separation,
}

impl CrossingKind {
    /// The name in the case file.
    pub fn name(self) -> &'static str {
        match self {
            Self::Reattachment => "reattachment",
            Self::Separation => "separation",
        }
    }
}

/// Force coefficients, `[output.coefficients]`: C = 2 F / (ρ U² A) of the
/// sum of the forces on `boundaries`.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Coefficients {
    /// Boundaries whose forces are summed.
    pub boundaries: Vec<String>,
    /// Reference velocity U (m/s).
    pub velocity: f64,
    /// Reference area A (m²).
    pub area: f64,
    /// Drag direction (unit vector).
    #[serde(default = "x_axis")]
    pub drag_direction: [f64; 3],
    /// Lift direction (unit vector).
    #[serde(default = "y_axis")]
    pub lift_direction: [f64; 3],
}

fn x_axis() -> [f64; 3] {
    [1.0, 0.0, 0.0]
}

fn y_axis() -> [f64; 3] {
    [0.0, 1.0, 0.0]
}

/// `[verification]`
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct VerificationSection {
    /// Refinement levels of `[mesh.generate]`; two or more give a sweep:
    /// with `[analytic]` the observed order of accuracy, with three or
    /// more the discretisation uncertainty of the monitored quantities
    /// (D67), also without an exact solution.
    pub levels: Vec<u32>,
    /// Accepted range of the observed order on the finest pair.
    pub expected_order: Option<[f64; 2]>,
    /// Largest accepted relative L2 velocity error (needs `[analytic]`).
    pub max_error: Option<f64>,
    /// Reference data of the lid-driven cavity.
    pub reference: Option<Reference>,
    /// Largest accepted deviation from the reference, in units of the lid
    /// speed.
    pub reference_tolerance: Option<f64>,
    /// Lid speed and cavity size for the reference comparison.
    pub lid_velocity: Option<f64>,
    /// Fails the case unless it converges.
    pub require_convergence: Option<bool>,
    /// Time steps of a time-order sweep on one mesh (transient, at least
    /// three, largest first): the order follows from the differences of
    /// successive solutions.
    pub time_steps: Vec<f64>,
    /// Accepted range of the observed time order on the smallest steps.
    pub expected_time_order: Option<[f64; 2]>,
    /// Benchmark data (the cylinder of Schäfer & Turek).
    pub benchmark: Option<Benchmark>,
    /// Boundary of the cylinder of the benchmark.
    pub benchmark_boundary: Option<String>,
    /// Largest accepted distance from the reference ranges, relative to
    /// the midpoint of each range.
    pub benchmark_tolerance: Option<f64>,
    /// Start of the evaluation window of a periodic benchmark (s; default:
    /// the second half of the run).
    pub evaluation_start: Option<f64>,
    /// Accepted range of the observed order of the scalar errors on the
    /// finest pair of a mesh sweep (scalars with `exact`, D51).
    pub scalar_order: Option<[f64; 2]>,
    /// Largest accepted L2 error of every scalar with `exact` (absolute).
    pub scalar_max_error: Option<f64>,
    /// Accepted range of every scalar of the report (the passive scalars
    /// and, with the energy equation, the temperature) on every mesh: the
    /// check of a bounded convection scheme (D237).
    pub scalar_range: Option<[f64; 2]>,
}

/// A line along which velocity and pressure are sampled.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Line {
    /// Name (file `line_<name>.csv`).
    pub name: String,
    /// Start point.
    pub from: [f64; 3],
    /// End point.
    pub to: [f64; 3],
    /// Number of points.
    #[serde(default = "default_points")]
    pub points: usize,
}

fn default_points() -> usize {
    101
}

/// `[output]`
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct OutputSection {
    /// Output directory (default `output/<case name>`).
    pub directory: Option<PathBuf>,
    /// Write `fields.vtkhdf` (velocity, pressure, error).
    pub vtkhdf: bool,
    /// Write a row of `residuals.csv` every N outer iterations.
    pub history_every: usize,
    /// Write `wall_<boundary>.csv` for every wall.
    pub wall_csv: bool,
    /// Sample lines.
    pub lines: Vec<Line>,
    /// Transient: write `fields_<step>.vtkhdf` every N steps (0: none; the
    /// final fields are `fields.vtkhdf`).
    pub fields_every: usize,
    /// Transient: also write the fields at these steps.
    pub fields_at: Vec<usize>,
    /// Transient: write `fields_<step>.vtkhdf` every this much simulated
    /// time (s), at the multiples of it; the time step is shortened to
    /// land on them (D66).
    pub fields_interval: Option<f64>,
    /// Transient: also write the fields at these simulated times (s),
    /// landed on exactly (D66).
    pub fields_times: Vec<f64>,
    /// A row of `forces.csv` and `probes.csv` every N steps (transient).
    pub monitor_every: usize,
    /// Transient: a row of `forces.csv` and `probes.csv` every this much
    /// simulated time (s) instead of every `monitor_every` steps, landed
    /// on exactly (D66): equally spaced samples for spectra.
    pub monitor_interval: Option<f64>,
    /// Point probes: velocity and pressure per step (`probes.csv`).
    pub probes: Vec<Probe>,
    /// Wall probes: wall shear stress, skin friction (with
    /// `[output.reference] velocity`) and y⁺ at points of a wall, in the
    /// report (and its uncertainty table of a level sweep).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub wall_probes: Vec<WallProbe>,
    /// Wall crossings: the first point along a direction where the wall
    /// shear changes sign (reattachment or separation), in the report (and
    /// its uncertainty table of a level sweep).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub wall_crossings: Vec<WallCrossing>,
    /// Force coefficients of some boundaries (`forces.csv`, the report).
    pub coefficients: Option<Coefficients>,
    /// Derived fields written with the fields (D72), from the cell
    /// gradients: vorticity, q-criterion, lambda2, enstrophy, strain-rate,
    /// pressure-coefficient, total-pressure, relative-velocity,
    /// frame-velocity, turbomachinery, relative-total-pressure, rothalpy
    /// (with `[[mrf]]` or `[[sliding]]`: the fields of rotating machinery,
    /// D183), wall-distance, wall-shear-stress, y-plus, turbulence-features
    /// (with a turbulence model: the invariant features of D166).
    pub derived: Vec<crate::derived::Derived>,
    /// Reference values of the pressure coefficient.
    pub reference: crate::derived::ReferenceValues,
    /// The turning zone (`zone` of an `[[mrf]]` or `[[sliding]]` table) whose
    /// frame the relative velocity of every cell is taken in (D183); default:
    /// each turning zone in its own frame, the rest in the fixed one.
    pub frame: Option<String>,
    /// Pumps, fans and turbines: flow, pressure rise, head, shaft power,
    /// efficiency and Euler work between an inlet and an outlet (D184).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub machines: Vec<crate::rotor::MachineSection>,
    /// Transient: running means, RMS and covariances of the fields and the
    /// statistics of the probes and force coefficients with batch-means
    /// confidence intervals and phase averages (D71).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub statistics: Option<crate::statistics::StatisticsSection>,
}

impl Default for OutputSection {
    fn default() -> Self {
        Self {
            directory: None,
            vtkhdf: true,
            history_every: 1,
            wall_csv: true,
            lines: Vec::new(),
            fields_every: 0,
            fields_at: Vec::new(),
            fields_interval: None,
            fields_times: Vec::new(),
            monitor_every: 1,
            monitor_interval: None,
            probes: Vec::new(),
            wall_probes: Vec::new(),
            wall_crossings: Vec::new(),
            coefficients: None,
            derived: Vec::new(),
            reference: Default::default(),
            frame: None,
            machines: Vec::new(),
            statistics: None,
        }
    }
}

/// A flow case.
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct FlowCase {
    /// Format of the case file: the Foehncraft version that introduced its
    /// key layout (default `"0.5"`, the layout before this key). An older
    /// file keeps running (its old keys are mapped, the log says which);
    /// `foehn case upgrade` rewrites it in the current format.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub case_format: Option<String>,
    /// What reading the case file did: one note per key of an older format
    /// that was mapped (printed in the log; not a key).
    #[serde(skip)]
    pub notes: Vec<String>,
    /// The mesh.
    pub mesh: MeshSection,
    /// Fluid properties.
    pub fluid: FluidSection,
    /// Materials of the case (D94): new ones, or changes of database
    /// materials with `base`.
    #[serde(skip_serializing_if = "crate::materials::CaseMaterials::is_empty")]
    #[schemars(schema_with = "crate::materials::case_materials_schema")]
    pub materials: crate::materials::CaseMaterials,
    /// Regions of cell zones with a kind and a material (D112); zones in no
    /// region are fluid with the `[fluid]` material.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub regions: BTreeMap<String, crate::materials::RegionSection>,
    /// Contact resistances and thin walls between two regions of a
    /// conjugate case (stage 7d2, D157): `[interfaces.<name>]`.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub interfaces: BTreeMap<String, crate::conduction::InterfaceSection>,
    /// Named constants of the expressions (D53): numbers, or expressions
    /// of other constants.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub constants: BTreeMap<String, Scalar>,
    /// Tables of the expressions (D53), piecewise linear in one argument.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub tables: BTreeMap<String, crate::expr::TableSection>,
    /// Flow model.
    pub flow: FlowSection,
    /// Start fields.
    pub initial: InitialSection,
    /// Pressure level.
    pub pressure: PressureSection,
    /// Turbulence model (stage 6): `[turbulence]`.
    #[serde(skip_serializing_if = "crate::turbulence::TurbulenceSection::is_default")]
    pub turbulence: crate::turbulence::TurbulenceSection,
    /// Boundary conditions per boundary.
    pub boundary: BTreeMap<String, BoundaryCondition>,
    /// Passive scalars (D51), transported by the flow: `[scalars.<name>]`.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub scalars: BTreeMap<String, crate::scalar::ScalarSection>,
    /// The scalar sub-tables of the boundaries (D116),
    /// `[boundary.<boundary>.scalar.<scalar>]`, per boundary: moved here
    /// from the boundary tables by the reader (part of the case hash).
    #[serde(skip)]
    pub boundary_scalar: BTreeMap<String, BTreeMap<String, crate::scalar::ScalarBoundary>>,
    /// The energy equation of the fluid (stage 7a1): `[energy]`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub energy: Option<crate::energy::EnergySection>,
    /// The energy sub-tables of the boundaries,
    /// `[boundary.<boundary>.energy]`: moved here from the boundary tables
    /// by the reader.
    #[serde(skip)]
    pub boundary_energy: BTreeMap<String, crate::energy::EnergyBoundary>,
    /// Gravity (stage 7b1): `[gravity]`. With it the solved pressure is the
    /// reduced pressure p′ = p − ρ_ref g·(x − x_ref) (D153).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gravity: Option<crate::buoyancy::GravitySection>,
    /// The buoyancy model (stage 7b1): `[buoyancy]`.
    #[serde(
        default,
        skip_serializing_if = "crate::buoyancy::BuoyancySection::is_default"
    )]
    pub buoyancy: crate::buoyancy::BuoyancySection,
    /// Lagrangian parcels of the sprays (L12): `[spray]`.
    #[serde(skip_serializing_if = "foehncraft_particles::SpraySection::is_empty")]
    pub spray: foehncraft_particles::SpraySection,
    /// Datums of the reports (SP1b, D262): `[[datum]]` tables, named planes,
    /// boundaries, interfaces, surfaces and volumes whose gas and liquid
    /// flows go to `datums.csv`.
    #[serde(rename = "datum", skip_serializing_if = "Vec::is_empty")]
    pub datums: Vec<crate::datum::DatumSection>,
    /// The fluid view of a conjugate heat transfer case (stage 7d1): set on
    /// the case of the flow, which runs on the fluid cells only, with the
    /// full case and mesh behind it ([`crate::cht`]).
    #[serde(skip)]
    #[schemars(skip)]
    pub fluid_view: Option<std::sync::Arc<crate::cht::FluidView>>,
    /// Solver settings.
    pub solver: Settings,
    /// Exact solution.
    pub analytic: Option<Analytic>,
    /// Checks.
    pub verification: VerificationSection,
    /// Threads, device, restart and re-read schedule.
    pub run: RunSection,
    /// Output.
    pub output: OutputSection,
    /// Time stepping; present in a transient case.
    pub time: Option<TimeSection>,
    /// Conservation ledger: report schedule and tolerances (D88, D89).
    pub balance: crate::balance::BalanceSection,
    /// How the console and `flow.log` look.
    pub log: crate::console::LogSection,
    /// Settings that change at prepared iterations, steps or times
    /// (D128): `[[schedule]]` entries, each with one trigger and a `set`
    /// table of dotted keys.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub schedule: Vec<crate::schedule::ScheduleEntry>,
    /// Zones of the mesh that rotate steadily (the moving reference frame,
    /// D55): `[[mrf]]` tables.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub mrf: Vec<crate::frame::MrfSection>,
    /// A zone of the mesh that turns through the time steps of a transient
    /// run (the sliding mesh, D56, D182): the `[[sliding]]` table.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub sliding: Vec<crate::sliding::SlidingSection>,
    /// The case file (for re-reading it during a run).
    #[serde(skip)]
    pub source: Option<PathBuf>,
    /// The `--set key=value` overrides applied to the case file.
    #[serde(skip)]
    pub overrides: Vec<String>,
    /// The pressure convention of the case file and the operating pressure
    /// of every fluid part (D239): set by the reader, which turns absolute
    /// values into the gauge values the solver works with; `None` for a
    /// case built in code (gauge values).
    #[serde(skip)]
    #[schemars(skip)]
    pub pressure_frame: Option<std::sync::Arc<crate::pressure_values::Frame>>,
}

/// Adds "did you mean" to the unknown-field and unknown-variant messages of
/// the TOML parser.
pub fn explain_toml_error(msg: &str) -> String {
    let mut out = msg.to_string();
    for key in ["unknown field `", "unknown variant `"] {
        let Some(i) = msg.find(key) else { continue };
        let rest = &msg[i + key.len()..];
        let Some(j) = rest.find('`') else { continue };
        let word = &rest[..j];
        let after = &rest[j + 1..];
        let end = after.find('\n').unwrap_or(after.len());
        let cands: Vec<&str> = after[..end].split('`').skip(1).step_by(2).collect();
        if let Some(c) = bc::did_you_mean(word, cands.iter().copied()) {
            out.push_str(&format!("\n  did you mean `{c}`?"));
        }
    }
    out
}

/// `run.restart = "latest"` (`--restart latest`): the newest restart file
/// of the output directory, by the step or iteration it holds (D64).
pub const LATEST: &str = "latest";

/// The diagnostic of a TOML error: the code from its kind, the location
/// and the key path from the span in `text` (D78).
pub fn toml_diagnostic(text: &str, e: &toml::de::Error) -> diag::Diagnostic {
    let msg = e.message();
    let code = if msg.contains("unknown field") {
        diag::CASE_UNKNOWN_KEY
    } else if msg.contains("missing field") {
        diag::CASE_MISSING_KEY
    } else if msg.contains("unknown variant")
        || msg.contains("invalid type")
        || msg.contains("invalid value")
        || msg.contains("invalid length")
    {
        diag::CASE_BAD_VALUE
    } else {
        diag::CASE_SYNTAX
    };
    let mut d = diag::Diagnostic::new(code, explain_toml_error(&e.to_string()));
    if let Some(span) = e.span()
        && let Some(before) = text.get(..span.start)
    {
        let line = before.matches('\n').count() + 1;
        let col = before.rsplit('\n').next().map_or(0, |l| l.chars().count()) + 1;
        d.location = Some([line, col]);
        // the table header above the span and the key at it
        let table = before
            .lines()
            .rev()
            .map(str::trim)
            .find(|l| l.starts_with('[') && !l.starts_with("[["))
            .map(|l| l.trim_matches(|c| c == '[' || c == ']').trim().to_string());
        let at = text.get(span.clone()).unwrap_or("").trim();
        let key: String = at
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || *c == '_' || *c == '-' || *c == '.')
            .collect();
        let line_text = before.rsplit('\n').next().unwrap_or("").trim();
        let is_key = !key.is_empty() && (line_text.is_empty() || !line_text.contains('='));
        d.key = match (table, is_key) {
            (Some(t), true) => Some(format!("{t}.{key}")),
            (None, true) => Some(key),
            (Some(t), false) => {
                let k = line_text.split('=').next().unwrap_or("").trim();
                Some(if k.is_empty() { t } else { format!("{t}.{k}") })
            }
            (None, false) => {
                let k = line_text.split('=').next().unwrap_or("").trim();
                (!k.is_empty()).then(|| k.to_string())
            }
        };
    }
    d
}

/// Applies one `key=value` override (`--set`, D62) to a case table: the key
/// is a dotted path (`fluid.viscosity`, `boundary.inlet.velocity`), the
/// value a TOML value (`1e-3`, `[1.0, 0.0, 0.0]`, `true`, `"rcm"`); a value
/// that is not TOML is taken as a string (`run.backend=cuda`).
pub fn apply_set(v: &mut toml::Table, set: &str) -> Result<()> {
    let bad = |m: String| diag::err(diag::Diagnostic::new(diag::CASE_SET, m));
    let Some((key, value)) = set.split_once('=') else {
        return Err(bad(format!("--set {set}: expected key=value")));
    };
    let key = key.trim();
    let value = value.trim();
    if key.is_empty() || key.split('.').any(str::is_empty) {
        return Err(bad(format!(
            "--set {set}: the key `{key}` is not a dotted path"
        )));
    }
    let parsed = toml::from_str::<toml::Table>(&format!("v = {value}"))
        .ok()
        .and_then(|mut t| t.remove("v"))
        .unwrap_or_else(|| toml::Value::String(value.to_string()));
    let mut parts: Vec<&str> = key.split('.').collect();
    let last = parts.pop().unwrap_or_default();
    let mut cur = v;
    let mut path = String::new();
    for p in parts {
        if !path.is_empty() {
            path.push('.');
        }
        path.push_str(p);
        let next = cur
            .entry(p.to_string())
            .or_insert_with(|| toml::Value::Table(Default::default()));
        let toml::Value::Table(t) = next else {
            return Err(bad(format!(
                "--set {set}: `{path}` is a value, not a table"
            )));
        };
        cur = t;
    }
    cur.insert(last.to_string(), parsed);
    Ok(())
}

impl FlowCase {
    /// Reads and validates a case file; relative paths are taken from its
    /// directory.
    pub fn from_file(path: &Path) -> Result<Self> {
        Self::from_file_with(path, &[])
    }

    /// As [`from_file`](Self::from_file), with `key=value` overrides
    /// applied before the case is validated (`--set`, D62). A re-read of
    /// the case file during the run applies them again.
    pub fn from_file_with(path: &Path, sets: &[String]) -> Result<Self> {
        let mut c = Self::from_file_unloaded(path, sets)?;
        c.load_tables()?;
        Ok(c)
    }

    /// As [`from_file_with`](Self::from_file_with), without reading the
    /// files the case refers to (the profile tables): the caller checks
    /// their paths first (read roots, D82) and then calls
    /// [`load_tables`](Self::load_tables).
    pub fn from_file_unloaded(path: &Path, sets: &[String]) -> Result<Self> {
        Self::from_file_unloaded_inner(path, sets).map_err(diag::case_error)
    }

    fn from_file_unloaded_inner(path: &Path, sets: &[String]) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(|e| {
            diag::err(diag::Diagnostic::new(
                diag::CASE_FILE,
                format!("reading {}: {e}", path.display()),
            ))
        })?;
        let mut c = Self::from_toml_parts(&text, path.parent().unwrap_or(Path::new(".")), sets)
            .with_context(|| format!("in case file {}", path.display()))?;
        c.source = Some(path.to_path_buf());
        Ok(c)
    }

    /// Reads the profile tables of the boundary conditions (their paths are
    /// resolved against the case file already).
    pub fn load_tables(&mut self) -> Result<()> {
        self.load_tables_inner().map_err(diag::case_error)
    }

    fn load_tables_inner(&mut self) -> Result<()> {
        for b in self.boundary.values_mut() {
            b.load_table(Path::new(""))?;
        }
        for (name, t) in &mut self.tables {
            t.load(name)?;
        }
        // the user material files and the data files of the materials
        self.materials.load()?;
        crate::spray::load(&mut self.spray)?;
        // the expressions can be checked once the tables are read
        self.validate_expressions()?;
        Ok(())
    }

    /// Parses and validates a case from TOML text.
    pub fn from_toml(text: &str, dir: &Path) -> Result<Self> {
        Self::from_toml_with(text, dir, &[])
    }

    /// Parses a case from TOML text, applies the `key=value` overrides and
    /// validates it.
    pub fn from_toml_with(text: &str, dir: &Path, sets: &[String]) -> Result<Self> {
        let mut c = Self::from_toml_parts(text, dir, sets).map_err(diag::case_error)?;
        c.load_tables()?;
        Ok(c)
    }

    /// As [`from_toml_with`](Self::from_toml_with), without reading the
    /// profile tables.
    fn from_toml_parts(text: &str, dir: &Path, sets: &[String]) -> Result<Self> {
        let mut v: toml::Table =
            toml::from_str(text).map_err(|e| diag::err(toml_diagnostic(text, &e)))?;
        for s in sets {
            apply_set(&mut v, s)?;
        }
        // keys of an older case format are mapped (after the overrides, so
        // that `--set` may use the old names of an old file)
        let has_convention = |v: &toml::Table| {
            v.get("fluid")
                .and_then(toml::Value::as_table)
                .is_some_and(|f| f.contains_key("pressure_values"))
        };
        let stated = has_convention(&v);
        let notes = foehncraft_mesh::case_format::migrate(&mut v).map_err(|m| {
            diag::err(diag::Diagnostic::new(diag::CASE_BAD_VALUE, m).key("case_format"))
        })?;
        // which pressure values the file gives (D239: a value left out sits
        // at the operating pressure)
        let given = crate::pressure_values::Given::of(&v);
        // the per-physics sub-tables of the boundaries (D116)
        let boundary_scalar = crate::scalar::extract_boundary_tables(&mut v)?;
        let boundary_energy = crate::energy::extract_boundary_tables(&mut v)?;
        // the changed default of D239 alone does not need the table: the
        // text keeps the line numbers of its errors
        let defaulted = !stated && has_convention(&v);
        let mut c: Self = if sets.is_empty()
            && notes.len() == usize::from(defaulted)
            && boundary_scalar.is_empty()
            && boundary_energy.is_empty()
        {
            toml::from_str(text).map_err(|e| diag::err(toml_diagnostic(text, &e)))?
        } else {
            Self::deserialize(toml::Value::Table(v)).map_err(|e| {
                let mut d = toml_diagnostic("", &e);
                if !sets.is_empty() {
                    d.message = format!("after --set: {}", d.message);
                }
                diag::err(d)
            })?
        };
        if defaulted {
            c.fluid.pressure_values = Some(crate::pressure_values::PressureValues::Gauge);
        }
        c.case_format = Some(foehncraft_mesh::case_format::CURRENT.into());
        c.boundary_scalar = boundary_scalar;
        c.boundary_energy = boundary_energy;
        c.notes = notes;
        c.overrides = sets.to_vec();
        if let Some(f) = &c.mesh.file {
            c.mesh.file = Some(dir.join(f));
        }
        if let Some(r) = &c.run.restart
            && r.as_os_str() != LATEST
        {
            c.run.restart = Some(dir.join(r));
        }
        for b in c.boundary.values_mut() {
            b.resolve_paths(dir);
        }
        crate::expr::resolve_paths(&mut c.tables, dir);
        c.materials.resolve_paths(dir);
        c.spray.resolve_paths(dir);
        crate::datum::resolve_paths(&mut c.datums, dir);
        c.solver.reduction_defaults = crate::turbulence::reduction_defaults(&c);
        // `[energy.cycle] enabled = false` runs the case fully coupled (the
        // reference of the loop, `--set energy.cycle.enabled=false`)
        if let Some(e) = &mut c.energy
            && e.cycle.as_ref().is_some_and(|y| !y.enabled)
        {
            e.cycle = None;
        }
        // as `[energy.frozen_flow] enabled = false` (D118)
        if let Some(e) = &mut c.energy
            && e.frozen_flow.as_ref().is_some_and(|y| !y.enabled)
        {
            e.frozen_flow = None;
        }
        c.validate()?;
        // absolute pressure values become the gauge values of the solver,
        // relative to the operating pressure of each fluid part (D239)
        crate::pressure_values::resolve(&mut c, &given)?;
        Ok(c)
    }

    /// The number of the pressure or total pressure of a pressure boundary
    /// (as the solver sees it: the gauge value, D239); `None` for an
    /// expression or another kind of boundary.
    pub fn boundary_pressure(&self, name: &str) -> Option<f64> {
        match self.boundary.get(name)? {
            BoundaryCondition::StaticPressureInlet { pressure, .. }
            | BoundaryCondition::StaticPressureOutlet { pressure, .. } => pressure.number(),
            BoundaryCondition::TotalPressureInlet { total_pressure, .. }
            | BoundaryCondition::TotalPressureOutlet { total_pressure, .. } => {
                total_pressure.number()
            }
            _ => None,
        }
    }

    /// The constants and tables of the expressions.
    pub fn env(&self) -> Result<Env> {
        Env::new(&self.constants, &self.tables)
    }

    /// The resolved property program (D94 (3)): the `[fluid]` material and
    /// the regions, evaluated at the case conditions.
    pub fn properties(&self) -> Result<crate::materials::PropertySet> {
        use crate::materials::{RegionKind, resolve_material};
        let env = self.env()?;
        let cond = self.fluid.conditions();
        let mut warnings = Vec::new();
        let fluid = resolve_material(
            self.fluid.material.as_deref(),
            &self.materials,
            self.fluid.density.as_ref(),
            self.fluid.viscosity.as_ref(),
            &env,
            cond,
            &mut warnings,
        )?;
        let mut regions = Vec::new();
        let mut seen: BTreeMap<&str, &str> = BTreeMap::new();
        for (name, r) in &self.regions {
            // solids belong to the energy equation (stage 7d1, crate::cht):
            // the flow runs on the fluid cells
            if (r.kind == RegionKind::Solid && self.energy.is_some())
                || r.kind == RegionKind::Porous
            {
                continue;
            }
            crate::need!(
                r.kind == RegionKind::Fluid,
                REG_NOT_FLUID,
                "[regions.{name}] kind = \"{}\": solid regions need the energy equation \
                 ([energy], stage 7d1)",
                match r.kind {
                    RegionKind::Solid => "solid",
                    RegionKind::Porous => "porous",
                    RegionKind::Fluid => "fluid",
                }
            );
            crate::need!(
                !r.zones.is_empty(),
                REG_ZONE,
                "[regions.{name}] zones must name at least one cell zone"
            );
            for z in &r.zones {
                if let Some(other) = seen.insert(z.as_str(), name.as_str()) {
                    crate::need!(
                        false,
                        REG_ZONE,
                        "[regions.{name}] zone `{z}` is also in region `{other}`"
                    );
                }
            }
            // a region of its own temperature evaluates its material there
            // (D194); the [fluid] shorthand numbers stay those of [fluid]
            let rcond = crate::materials::Conditions {
                temperature: r.temperature.unwrap_or(cond.temperature),
                pressure: r.operating_pressure.unwrap_or(cond.pressure),
            };
            if let Some(t) = r.temperature {
                crate::need!(
                    t > 0.0 && t.is_finite(),
                    CASE_RANGE,
                    "[regions.{name}] temperature must be > 0 K, is {t}"
                );
            }
            let m = match &r.material {
                Some(mat) => resolve_material(
                    Some(mat),
                    &self.materials,
                    None,
                    None,
                    &env,
                    rcond,
                    &mut warnings,
                )?,
                None if r.temperature.is_some() || r.operating_pressure.is_some() => {
                    resolve_material(
                        self.fluid.material.as_deref(),
                        &self.materials,
                        self.fluid.density.as_ref(),
                        self.fluid.viscosity.as_ref(),
                        &env,
                        rcond,
                        &mut warnings,
                    )?
                }
                None => fluid.clone(),
            };
            if let Some(rho) = m.density.value.constant() {
                crate::need!(
                    rho > 0.0 && rho.is_finite(),
                    CASE_RANGE,
                    "[regions.{name}] the density of material `{}` must be > 0, is {rho}",
                    r.material.as_deref().unwrap_or("")
                );
            }
            regions.push((name.clone(), r.zones.clone(), m));
        }
        Ok(crate::materials::PropertySet {
            conditions: cond,
            fluid,
            regions,
            warnings,
        })
    }

    /// The density ρ (kg/m³).
    pub fn rho(&self) -> f64 {
        self.properties().map_or(f64::NAN, |p| p.density())
    }

    /// The reference viscosity μ (Pa s): the `[fluid]` value (Reynolds
    /// numbers, exact solutions).
    pub fn mu(&self) -> f64 {
        self.properties()
            .map_or(f64::NAN, |p| p.reference_viscosity())
    }

    /// The body force as compiled values.
    pub fn body_force_values(&self) -> Result<[Value; 3]> {
        let allowed = if self.time.is_some() {
            Vars::ALL
        } else {
            Vars::SPACE
        };
        crate::expr::vector_value(
            &self.flow.body_force,
            &self.env()?,
            allowed,
            "flow.body_force",
        )
    }

    /// The compiled boundary values (D53).
    pub fn compiled_boundaries(&self) -> Result<bc::Compiled> {
        bc::compile_all(&self.boundary, &self.env()?, self.time.is_some())
    }

    /// The uniform body force, when all three components are numbers.
    pub fn uniform_body_force(&self) -> Option<[f64; 3]> {
        crate::expr::numbers(&self.flow.body_force)
    }

    /// Checks every expression of the case: syntax, names, the variables
    /// each key allows (D53), and the property program (D94).
    pub fn validate_expressions(&self) -> Result<()> {
        let env = self.env()?;
        self.body_force_values()?;
        let transient = self.time.is_some();
        for (name, b) in &self.boundary {
            b.compile(name, &env, transient)?;
        }
        crate::scalar::validate(
            &self.scalars,
            &self.boundary_scalar,
            &self.boundary,
            &self.solver.scalar,
            &env,
            transient,
        )?;
        crate::energy::validate(self, &env)?;
        crate::buoyancy::validate(self)?;
        crate::density::validate(self)?;
        crate::fluids::validate(self)?;
        crate::porous::validate(self)?;
        crate::turbulence::validate(self)?;
        let p = self.properties()?;
        let rho = p.density();
        crate::need!(
            rho > 0.0 && rho.is_finite(),
            CASE_RANGE,
            "[fluid] density must be > 0"
        );
        if let Some(mu) = p.fluid.viscosity.value.constant() {
            crate::need!(
                mu > 0.0 && mu.is_finite(),
                CASE_RANGE,
                "[fluid] viscosity must be > 0"
            );
        }
        Ok(())
    }

    /// The resolved case (`case.resolved.toml`): every default filled in,
    /// the fluid's density and viscosity as resolved, and the resolved
    /// property program as comments (D94).
    pub fn resolved_toml(&self) -> Result<String> {
        let mut c = self.clone();
        let p = self.properties()?;
        if c.fluid.density.is_none() {
            c.fluid.density = Some(Scalar::Number(p.density()));
        }
        if c.fluid.viscosity.is_none()
            && let Some(mu) = p.fluid.viscosity.value.constant()
        {
            c.fluid.viscosity = Some(Scalar::Number(mu));
        }
        let mut text = if c.boundary_scalar.is_empty() && c.boundary_energy.is_empty() {
            toml::to_string_pretty(&c)?
        } else {
            // the scalar and energy sub-tables live outside the flow
            // conditions: put them back into their `[boundary.<name>]`
            // tables
            let mut v = toml::Table::try_from(&c)?;
            if let Some(toml::Value::Table(bs)) = v.get_mut("boundary") {
                for (b, t) in &c.boundary_scalar {
                    if let Some(toml::Value::Table(bt)) = bs.get_mut(b) {
                        bt.insert("scalar".into(), toml::Value::try_from(t)?);
                    }
                }
                for (b, t) in &c.boundary_energy {
                    if let Some(toml::Value::Table(bt)) = bs.get_mut(b) {
                        bt.insert("energy".into(), toml::Value::try_from(t)?);
                    }
                }
            }
            toml::to_string_pretty(&v)?
        };
        text.push('\n');
        text.push_str(&p.text());
        Ok(text)
    }

    /// Ranges and consistency of the settings (the expressions and the
    /// file tables are checked by [`load_tables`](Self::load_tables)).
    pub fn validate(&self) -> Result<()> {
        crate::frame::validate(self)?;
        crate::start::validate(self)?;
        crate::sliding::validate(self)?;
        crate::rotor::validate(self)?;
        crate::spray::validate(self)?;
        crate::datum::validate(self)?;
        if let Some(d) = &self.fluid.density
            && let Some(v) = d.number()
        {
            crate::need!(v > 0.0, CASE_RANGE, "[fluid] density must be > 0");
        }
        if let Some(d) = &self.fluid.viscosity
            && let Some(v) = d.number()
        {
            crate::need!(v > 0.0, CASE_RANGE, "[fluid] viscosity must be > 0");
        }
        if let Some(t) = self.fluid.temperature {
            crate::need!(
                t > 0.0 && t.is_finite(),
                CASE_RANGE,
                "[fluid] temperature must be > 0 (K)"
            );
        }
        if let Some(p) = self.fluid.operating_pressure {
            crate::need!(
                p > 0.0 && p.is_finite(),
                CASE_RANGE,
                "[fluid] operating_pressure must be > 0 (Pa)"
            );
        }
        if self
            .tables
            .values()
            .all(|t| t.file.is_none() || t.data.is_some())
        {
            self.validate_expressions()?;
        }
        let transient = self.time.is_some();
        self.solver.validate(transient)?;
        self.run.limits.validate()?;
        crate::schedule::validate(self)?;
        self.balance.validate()?;
        crate::logcols::validate(self)?;
        if let Some(t) = &self.time {
            t.validate()?;
        }
        crate::need!(
            self.output.history_every >= 1,
            CASE_RANGE,
            "[output] history_every must be >= 1"
        );
        crate::need!(
            self.output.monitor_every >= 1,
            CASE_RANGE,
            "[output] monitor_every must be >= 1"
        );
        // D66: output in simulated time
        let (o, r) = (&self.output, &self.run);
        for (key, v) in [
            ("output.fields_interval", o.fields_interval),
            ("output.monitor_interval", o.monitor_interval),
            ("run.restart_interval", r.restart_interval),
        ] {
            if let Some(v) = v {
                crate::need!(
                    v > 0.0 && v.is_finite(),
                    CASE_RANGE,
                    "{key} must be > 0 (s)"
                );
                crate::need!(
                    transient,
                    CASE_INCONSISTENT,
                    "{key} schedules output in simulated time and needs a [time] table"
                );
            }
        }
        for (key, v) in [
            ("output.fields_times", &o.fields_times),
            ("run.restart_times", &r.restart_times),
        ] {
            if !v.is_empty() {
                crate::need!(
                    v.iter().all(|x| x.is_finite()),
                    CASE_RANGE,
                    "{key} must hold finite times (s)"
                );
                crate::need!(
                    transient,
                    CASE_INCONSISTENT,
                    "{key} schedules output in simulated time and needs a [time] table"
                );
            }
        }
        if let Some(a) = &self.analytic
            && a.is_transient()
            && !transient
        {
            anyhow::bail!(
                "[analytic] type = \"taylor-green\" is transient and needs a [time] table"
            );
        }
        if let Some(c) = &self.output.coefficients {
            crate::need!(
                c.velocity > 0.0 && c.area > 0.0,
                CASE_RANGE,
                "[output.coefficients] velocity and area must be > 0"
            );
            crate::need!(
                !c.boundaries.is_empty(),
                CASE_INCONSISTENT,
                "[output.coefficients] boundaries must name at least one boundary"
            );
        }
        crate::derived::validate(self)?;
        if let Some(st) = &self.output.statistics {
            st.validate(transient)?;
        }
        let v = &self.verification;
        if !v.time_steps.is_empty() {
            crate::need!(
                transient,
                CASE_INCONSISTENT,
                "[verification] time_steps need a [time] table"
            );
            crate::need!(
                v.time_steps.len() >= 3,
                CASE_INCONSISTENT,
                "[verification] time_steps need at least three steps for an observed order"
            );
            crate::need!(
                v.levels.len() < 2,
                CASE_INCONSISTENT,
                "[verification] a sweep is either over levels or over time_steps"
            );
            crate::diag::check_fix(
                self.time.is_some_and(|t| !t.adaptive),
                crate::diag::CASE_INCONSISTENT,
                &["time.adaptive=false"],
                || {
                    "[verification] time_steps need a fixed time step ([time] adaptive = false)"
                        .into()
                },
            )?;
        }
        if let Some(b) = v.benchmark {
            let re = self.rho() * b.mean_velocity() * Benchmark::DIAMETER / self.mu();
            crate::need!(
                (re - b.reynolds()).abs() < 1e-6 * b.reynolds(),
                CASE_INCONSISTENT,
                "[verification] benchmark {} needs Re = ρ Ū D / μ = {} (Ū = {} m/s, D = 0.1 m), \
                 the case has {re}",
                b.label(),
                b.reynolds(),
                b.mean_velocity()
            );
            if b == Benchmark::SchaferTurek2d2 {
                crate::need!(
                    transient,
                    CASE_INCONSISTENT,
                    "[verification] benchmark schafer-turek-2d2 is periodic and needs a [time] table"
                );
            }
        }
        let v = &self.verification;
        if v.levels.len() >= 2 {
            crate::need!(
                self.mesh.generate.is_some()
                    || self.mesh.cylinder.is_some()
                    || self.mesh.annulus.is_some()
                    || self.mesh.revolved.is_some(),
                CASE_INCONSISTENT,
                "[verification] levels need a [mesh.generate], [mesh.cylinder], [mesh.annulus] or \
                 [mesh.revolved] table"
            );
        }
        if v.expected_order.is_some() {
            crate::need!(
                self.analytic.is_some(),
                CASE_INCONSISTENT,
                "[verification] expected_order needs an [analytic] solution (without one the sweep \
                 reports the observed orders of the monitored quantities, D67)"
            );
        }
        if v.max_error.is_some() {
            crate::need!(
                self.analytic.is_some(),
                CASE_INCONSISTENT,
                "[verification] max_error needs an [analytic] solution"
            );
        }
        if let Some(reference) = v.reference {
            let re = self.rho() * v.lid_velocity.unwrap_or(1.0) / self.mu();
            let want = reference.reynolds();
            crate::need!(
                (re - want).abs() < 1e-6 * want,
                CASE_INCONSISTENT,
                "[verification] reference {:?} is for Re = {want}, the case has Re = ρ U / μ = {re} \
                 (unit cavity, lid_velocity = {})",
                reference,
                v.lid_velocity.unwrap_or(1.0)
            );
        }
        for l in &self.output.lines {
            crate::need!(
                l.points >= 2,
                CASE_RANGE,
                "[[output.lines]] {}: points must be >= 2",
                l.name
            );
        }
        for w in &self.output.wall_crossings {
            crate::need!(
                norm(w.direction) > 0.0 && w.direction.iter().all(|x| x.is_finite()),
                CASE_RANGE,
                "[[output.wall_crossings]] {}: direction must be a non-zero vector",
                w.name
            );
            crate::need!(
                w.width.is_none_or(|b| b > 0.0),
                CASE_RANGE,
                "[[output.wall_crossings]] {}: width must be > 0",
                w.name
            );
            crate::need!(
                w.skip.is_none_or(|b| b >= 0.0 && b.is_finite()),
                CASE_RANGE,
                "[[output.wall_crossings]] {}: skip must be >= 0",
                w.name
            );
        }
        Ok(())
    }

    pub(crate) fn analytic(&self) -> Option<Analytic> {
        self.analytic.map(|a| a.with_fluid(self.rho(), self.mu()))
    }

    /// What the boundary-condition check needs to know about the case.
    pub fn check_context(&self) -> bc::CheckContext {
        bc::CheckContext {
            laminar: self.turbulence.is_laminar(),
            transient: self.time.is_some(),
            potential_flow: self.initial.potential_flow,
            body_force: self.uniform_body_force() != Some([0.0; 3]),
            zeta: self.turbulence.model.has_zeta(),
            interface: self
                .sliding
                .iter()
                .filter(|e| e.imprinted())
                .flat_map(|e| [e.rotor_interface.clone(), e.stator_interface.clone()])
                .flatten()
                .collect(),
        }
    }

    /// The zones that turn: the `[[mrf]]` and `[[sliding]]` entries (their
    /// walls move normal to themselves with them, [`bc::Turning`]); an entry
    /// whose speed is not valid is left out (its own check reports it).
    pub fn turning(&self) -> Vec<bc::Turning> {
        let make = |zone: &Option<String>, omega: [f64; 3], origin: [f64; 3], axis: [f64; 3]| {
            let w = norm(omega);
            let a = if w > 0.0 { omega } else { axis };
            let n = norm(a);
            bc::Turning {
                zone: zone.clone(),
                omega,
                origin,
                axis: a.map(|v| v / n),
            }
        };
        let mrf = self.mrf.iter().filter_map(|m| {
            m.omega_vector()
                .ok()
                .map(|omega| make(&m.zone, omega, m.origin, m.axis))
        });
        let sliding = self.sliding.iter().filter_map(|m| {
            m.omega_vector()
                .ok()
                .map(|omega| make(&m.zone, omega, m.origin, m.axis))
        });
        mrf.chain(sliding).collect()
    }

    /// The exact solution at time `t` (transient solutions decay in time).
    pub fn analytic_at(&self, t: f64) -> Option<Analytic> {
        self.analytic().map(|a| a.at_time(t))
    }

    /// Reynolds number ρ U L / μ: the lid speed of the cavity reference,
    /// the scales of the exact solution or the mean velocity and the
    /// cylinder diameter of the benchmark.
    pub fn reynolds(&self) -> Option<f64> {
        let rho = self.rho();
        let mu = self.mu();
        if let Some(b) = self.verification.benchmark {
            return Some(rho * b.mean_velocity() * crate::analytic::Benchmark::DIAMETER / mu);
        }
        self.verification
            .reference
            .map(|_| rho * self.verification.lid_velocity.unwrap_or(1.0) / mu)
            .or_else(|| {
                self.analytic()
                    .map(|a| rho * a.velocity_scale() * a.length_scale() / mu)
            })
    }
}

/// Error norms against the exact solution.
#[derive(Debug, Clone, Serialize)]
pub struct ErrorNorms {
    /// Cells.
    pub cells: usize,
    /// Mesh size h (cube root of the mean cell volume; square root of the
    /// mean cell area in 2D).
    pub h: f64,
    /// L2 velocity error over the velocity scale.
    pub l2: f64,
    /// Largest velocity error over the velocity scale.
    pub linf: f64,
    /// Largest pressure error over the pressure range of the exact solution
    /// (absolute when the exact pressure is constant).
    pub pressure_linf: f64,
}

/// Comparison with the reference data.
#[derive(Debug, Clone, Serialize)]
pub struct ReferenceComparison {
    /// The data.
    pub reference: String,
    /// Largest |u − u_ref| / U on the vertical centre line.
    pub u_max_dev: f64,
    /// Largest |v − v_ref| / U on the horizontal centre line.
    pub v_max_dev: f64,
    /// (y, u/U) computed at the reference points.
    pub u: Vec<(f64, f64)>,
    /// (x, v/U) computed at the reference points.
    pub v: Vec<(f64, f64)>,
    /// (y, u/U) on 201 points of the vertical centre line, for the chart.
    #[serde(skip)]
    pub u_line: Vec<(f64, f64)>,
    /// (x, v/U) on 201 points of the horizontal centre line.
    #[serde(skip)]
    pub v_line: Vec<(f64, f64)>,
}

/// Forces and fluxes of one boundary.
#[derive(Debug, Clone, Serialize)]
pub struct BoundarySummary {
    /// Name.
    pub name: String,
    /// Condition.
    pub kind: String,
    /// Area (m²).
    pub area: f64,
    /// Mass flow out of the domain (kg/s).
    pub mass_flow: f64,
    /// Area-weighted mean pressure (Pa).
    pub mean_pressure: f64,
    /// Force of the fluid on the boundary, pressure and viscous (N).
    pub force: [f64; 3],
    /// Torque of the fluid on a wall about the axis of its rotor, or of the
    /// first turning zone (N m, D184); none without a turning zone.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub torque: Option<f64>,
}

/// One solve (one level).
#[derive(Debug, Clone, Serialize)]
pub struct Solve {
    /// Mesh description.
    pub mesh: String,
    /// SHA-256 of the mesh (D85).
    pub mesh_hash: String,
    /// Where it ran.
    pub device: String,
    /// Cells.
    pub cells: usize,
    /// Outer iterations.
    pub iterations: usize,
    /// Converged to the tolerance.
    pub converged: bool,
    /// Last residuals.
    pub residuals: Residuals,
    /// Phase timings.
    pub timings: Timings,
    /// Multigrid setups.
    pub amg_setups: usize,
    /// Multigrid setups of the momentum sweeps of the coupled algorithm
    /// (D179).
    pub sweep_setups: usize,
    /// Multigrid setups of the transport equations of the turbulence
    /// model (D180).
    pub turbulence_setups: usize,
    /// Nanoseconds per cell and outer iteration (R6).
    pub ns_per_cell_iteration: f64,
    /// Error norms.
    pub errors: Option<ErrorNorms>,
    /// Reference comparison.
    pub reference: Option<ReferenceComparison>,
    /// Boundaries.
    pub boundaries: Vec<BoundarySummary>,
    /// Net mass flow out of the domain over the largest boundary mass flow (in
    /// a closed domain over ρ |u|max times the largest boundary area).
    pub mass_imbalance: f64,
    /// Force coefficients (drag, lift) of `[output.coefficients]` at the
    /// end.
    pub coefficients: Option<[f64; 2]>,
    /// Torque, shaft power and machine performance of the turning zones at
    /// the end (D184).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rotors: Option<crate::rotor::Values>,
    /// Transient run.
    pub transient: Option<crate::unsteady::TransientSummary>,
    /// The datums at the end of the run (SP1b, D267): the gas report of the
    /// last row and the liquid through each over the run.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub datums: Option<Vec<crate::datum::DatumSummary>>,
    /// The sprays at the end of the run (SP1b, D265, D267).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub spray: Option<crate::spray::SpraySummary>,
    /// Benchmark comparison.
    pub benchmark: Option<crate::monitor::BenchmarkResult>,
    /// Conservation ledger (D88, D89).
    pub balance: Option<BalanceSummary>,
    /// Cell-quality check and remediation (D50), on a mesh with poor faces.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remediation: Option<crate::simplec::RemediationReport>,
    /// The operator exactness check (D236).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub operators: Option<crate::operators::OperatorReport>,
    /// The passive scalars (D51).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub scalars: Vec<crate::scalar::ScalarSummary>,
    /// The periodic pairs (D29): mass flow and driving force.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub periodic: Vec<PeriodicSummary>,
    /// The wall probes of `[[output.wall_probes]]`.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub wall_probes: Vec<WallProbeValue>,
    /// The wall crossings of `[[output.wall_crossings]]`.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub wall_crossings: Vec<WallCrossingValue>,
    /// The turbulence model (stage 6): first-cell y⁺ per wall, ν_t/ν max.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub turbulence: Option<crate::turbulence::TurbulenceSummary>,
    /// The statistics of a transient run (D71).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub statistics: Option<crate::statistics::StatisticsSummary>,
    /// The variable density (stage 7c): p₀, the mass, the density range
    /// and the largest Mach number.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub density: Option<crate::density::Summary>,
}

/// The value of a wall probe at the end of a solve.
#[derive(Debug, Clone, Serialize)]
pub struct WallProbeValue {
    /// Name.
    pub name: String,
    /// The point on the wall the values are taken at (m).
    pub at: [f64; 3],
    /// Wall shear stress (Pa).
    pub tau: [f64; 3],
    /// |τ| (Pa).
    pub magnitude: f64,
    /// Skin friction |τ| / (½ ρ_ref U_ref²), with a reference velocity.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub skin_friction: Option<f64>,
    /// First-cell y⁺ there (with a turbulence model: of its wall treatment).
    pub y_plus: f64,
    /// The passive scalars there: wall value and diffusive flux into the
    /// domain (the Nusselt number of a heated wall, stage 6 part 6g).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub scalars: Vec<WallProbeScalar>,
}

/// A passive scalar at a wall probe.
#[derive(Debug, Clone, Serialize)]
pub struct WallProbeScalar {
    /// The scalar.
    pub name: String,
    /// Its value on the wall.
    pub wall: f64,
    /// The diffusive flux into the domain per unit area (the unit of the
    /// scalar times m/s).
    pub flux: f64,
}

/// The wall value and the diffusive flux into the domain of the scalar
/// `st` at the wall face `f`: φ_b = a φ_P + b, flux = Γ_f/ρ (|S|/d_n)
/// (φ_b − φ_P). None where the scalar has no wall data.
fn scalar_wall_face(
    fl: &Solution,
    geo: &Geometry,
    st: &crate::scalar::ScalarState,
    f: usize,
    ni: usize,
) -> Option<(f64, f64)> {
    let w = st.wall.as_ref()?;
    let fb = f - ni;
    let c = fl.disc.owner[f] as usize;
    let pp = st.phi[c];
    let pb = w.a[fb] * pp + w.b[fb];
    let ex = w.gamma[fb] / fl.phys.rho_face_at(f) * fl.disc.dgeo[f] / norm(geo.face_area[f]);
    Some((pb, ex * (pb - pp)))
}

/// The log line of a wall probe.
pub(crate) fn wall_probe_line(w: &WallProbeValue) -> String {
    format!(
        "wall probe {} at ({:.6e}, {:.6e}, {:.6e}) m: |τ_w| {:.6e} Pa{}, y⁺ {:.3}",
        w.name,
        w.at[0],
        w.at[1],
        w.at[2],
        w.magnitude,
        w.skin_friction
            .map(|c| format!(", C_f {c:.6e}"))
            .unwrap_or_default(),
        w.y_plus
    ) + &w
        .scalars
        .iter()
        .map(|s| format!(", {} wall {:.6e}, flux {:.6e}", s.name, s.wall, s.flux))
        .collect::<String>()
}

/// A wall crossing at the end of a solve.
#[derive(Debug, Clone, Serialize)]
pub struct WallCrossingValue {
    /// Name.
    pub name: String,
    /// `reattachment` or `separation`.
    pub kind: &'static str,
    /// Distance from `from` along `direction` (m); none when the shear
    /// does not change sign that way.
    pub position: Option<f64>,
    /// The point on the wall (m), between the two face centroids.
    pub at: Option<[f64; 3]>,
}

/// The log line of a wall crossing.
pub(crate) fn wall_crossing_line(w: &WallCrossingValue, c: &WallCrossing) -> String {
    match (w.position, w.at) {
        (Some(s), Some(x)) => format!(
            "wall crossing {} ({}) on {}: {:.6e} m from ({:.6e}, {:.6e}, {:.6e}) m, at \
             ({:.6e}, {:.6e}, {:.6e}) m",
            w.name, w.kind, c.boundary, s, c.from[0], c.from[1], c.from[2], x[0], x[1], x[2]
        ),
        _ => format!(
            "wall crossing {} ({}) on {}: the shear does not change sign that way",
            w.name, w.kind, c.boundary
        ),
    }
}

/// The wall crossings of `case` in the solution `fl`: per crossing the
/// faces of the boundary after `from` and `skip` (and within `width` of
/// the line),
/// ordered by the projection of their centroids on the direction (in face
/// order on ties), and the first pair whose shear along the direction
/// changes sign the given way, linear between the two.
pub(crate) fn wall_crossings(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
) -> Vec<WallCrossingValue> {
    let ni = mesh.n_interior_faces();
    let mut out = Vec::new();
    for wc in &case.output.wall_crossings {
        let n = norm(wc.direction);
        let d: [f64; 3] = std::array::from_fn(|i| wc.direction[i] / n);
        let mut faces: Vec<(f64, usize, f64)> = Vec::new();
        if let Some(p) = mesh.boundaries.iter().find(|p| p.name == wc.boundary) {
            for f in p.range() {
                let r = sub(geo.face_centroid[f], wc.from);
                let s = dot(r, d);
                if s < wc.skip.unwrap_or(0.0) {
                    continue;
                }
                if let Some(b) = wc.width {
                    let off: [f64; 3] = std::array::from_fn(|i| r[i] - s * d[i]);
                    if norm(off) > b {
                        continue;
                    }
                }
                let (_, _, tau) = wall_face_values(fl, geo, f, ni);
                faces.push((s, f, dot(tau, d)));
            }
        }
        faces.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let turns = |t0: f64, t1: f64| match wc.kind {
            CrossingKind::Reattachment => t0 < 0.0 && t1 >= 0.0,
            CrossingKind::Separation => t0 > 0.0 && t1 <= 0.0,
        };
        let hit = faces.windows(2).find(|w| turns(w[0].2, w[1].2)).map(|w| {
            let ((s0, f0, t0), (s1, f1, t1)) = (w[0], w[1]);
            let a = -t0 / (t1 - t0);
            let (x0, x1) = (geo.face_centroid[f0], geo.face_centroid[f1]);
            (
                s0 + a * (s1 - s0),
                std::array::from_fn(|i| x0[i] + a * (x1[i] - x0[i])),
            )
        });
        out.push(WallCrossingValue {
            name: wc.name.clone(),
            kind: wc.kind.name(),
            position: hit.map(|h| h.0),
            at: hit.map(|h| h.1),
        });
    }
    out
}

/// The turbulence summary of a solve, its log lines and findings.
pub(crate) fn turbulence_summary<F: FnMut(&str)>(
    fl: &mut Flow,
    mesh: &Mesh,
    lg: &mut Log<F>,
    warnings: &mut Vec<diag::Diagnostic>,
) -> Result<Option<crate::turbulence::TurbulenceSummary>> {
    let Some((sm, finds)) = crate::turbulence::summary(fl, mesh)? else {
        return Ok(None);
    };
    lg.line(&format!(
        "turbulence: {}; largest ν_t/ν {:.4e} in cell {}",
        sm.model, sm.viscosity_ratio_max, sm.viscosity_ratio_cell
    ));
    for w in &sm.walls {
        lg.line(&format!(
            "  wall {}: first-cell y⁺ min {:.3}, mean {:.3}, max {:.3}",
            w.boundary, w.min, w.mean, w.max
        ));
    }
    let h = &sm.health;
    let zeta = if h.floor_cells_max.len() > 3 {
        format!(
            ", {} at most {} cells ({}, {})",
            if sm
                .model
                .starts_with(crate::turbulence::TurbulenceModel::EllipticBlendingKEpsilon.name(),)
            {
                "φ"
            } else {
                "ζ"
            },
            h.floor_cells_max[3],
            h.floor_iterations[3],
            h.floor_iterations_late[3]
        )
    } else {
        String::new()
    };
    lg.line(&format!(
        "  floors: k at most {} cells ({} solves above 0.1 %, {} in the final tenth), {} at \
         most {} cells ({}, {}){zeta}; ν_t capped in at most {} cells",
        h.floor_cells_max[0],
        h.floor_iterations[0],
        h.floor_iterations_late[0],
        if sm.second == "epsilon" { "ε" } else { "ω" },
        h.floor_cells_max[1],
        h.floor_iterations[1],
        h.floor_iterations_late[1],
        h.capped_max
    ));
    for d in finds {
        lg.note(&d.message);
        warnings.push(d);
    }
    Ok(Some(sm))
}

/// The wall probes of `case` in the solution `fl`.
pub(crate) fn wall_probes(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
) -> Vec<WallProbeValue> {
    let ni = mesh.n_interior_faces();
    let q = case
        .output
        .reference
        .velocity
        .or(case.output.coefficients.as_ref().map(|c| c.velocity))
        .map(|u| 0.5 * case.output.reference.density.unwrap_or(case.rho()) * u * u);
    let mut out = Vec::new();
    for wp in &case.output.wall_probes {
        let Some(p) = mesh.boundaries.iter().find(|p| p.name == wp.boundary) else {
            continue;
        };
        // the two nearest face centroids, in face order on ties
        let mut d: Vec<(f64, usize)> = p
            .range()
            .map(|f| (norm(sub(geo.face_centroid[f], wp.point)), f))
            .collect();
        d.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let Some(&(_, f1)) = d.first() else { continue };
        let f2 = d.get(1).map_or(f1, |x| x.1);
        let (x1, x2) = (geo.face_centroid[f1], geo.face_centroid[f2]);
        let e = sub(x2, x1);
        let ee = dot(e, e);
        let t = if ee > 0.0 {
            (dot(sub(wp.point, x1), e) / ee).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let val = |f: usize| -> ([f64; 3], f64) {
            let (_, _, tau) = wall_face_values(fl, geo, f, ni);
            let yp = match &fl.state.turbulence {
                Some(ts) => ts.y_plus[f - ni],
                None => {
                    let c = fl.disc.owner[f] as usize;
                    let s = geo.face_area[f];
                    let y = dot(sub(geo.face_centroid[f], geo.cell_centroid[c]), s).abs() / norm(s);
                    let rho = fl.phys.rho_face_at(f);
                    let ut = (norm(tau) / rho).sqrt();
                    y * ut * rho / fl.phys.mu_at(f)
                }
            };
            (tau, yp)
        };
        let (t1, y1) = val(f1);
        let (t2, y2) = val(f2);
        let tau: [f64; 3] = std::array::from_fn(|i| (1.0 - t) * t1[i] + t * t2[i]);
        let magnitude = norm(tau);
        let scalars = fl
            .state
            .scalars
            .iter()
            .filter_map(|st| {
                let (w1, q1) = scalar_wall_face(fl, geo, st, f1, ni)?;
                let (w2, q2) = scalar_wall_face(fl, geo, st, f2, ni)?;
                Some(WallProbeScalar {
                    name: st.name.clone(),
                    wall: (1.0 - t) * w1 + t * w2,
                    flux: (1.0 - t) * q1 + t * q2,
                })
            })
            .collect();
        out.push(WallProbeValue {
            name: wp.name.clone(),
            at: std::array::from_fn(|i| x1[i] + t * e[i]),
            tau,
            magnitude,
            skin_friction: q.map(|q| magnitude / q),
            y_plus: (1.0 - t) * y1 + t * y2,
            scalars,
        });
    }
    out
}

/// The flow through one periodic pair at the end of a solve (D29).
#[derive(Debug, Clone, Serialize)]
pub struct PeriodicSummary {
    /// The boundary the flow enters through (side `a`).
    pub from: String,
    /// The boundary the flow leaves through (side `b`).
    pub to: String,
    /// Mass flow from `from` towards `to` (kg/s).
    pub mass_flow: f64,
    /// The driving force per unit volume along the translation (N/m³),
    /// the mean pressure gradient it stands for with the opposite sign
    /// (Pa/m).
    pub driving_force: f64,
}

/// The periodic pairs of `fl` (D29) for the report, with their log lines.
pub(crate) fn periodic_summaries<F: FnMut(&str)>(
    fl: &mut Flow,
    mesh: &Mesh,
    lg: &mut Log<F>,
) -> Result<Vec<PeriodicSummary>> {
    let Some(per) = &mesh.periodic else {
        return Ok(Vec::new());
    };
    let q = fl.periodic_mass_flow()?;
    let b = fl.periodic_force();
    let out: Vec<PeriodicSummary> = per
        .pairs
        .iter()
        .zip(q.iter().zip(&b))
        .map(|(p, (&q, &b))| PeriodicSummary {
            from: p.a.clone(),
            to: p.b.clone(),
            mass_flow: q,
            driving_force: b,
        })
        .collect();
    for s in &out {
        lg.line(&format!(
            "periodic {}–{}: mass flow {:.6e} kg/s, driving force {:.6e} N/m³",
            s.from, s.to, s.mass_flow, s.driving_force
        ));
    }
    Ok(out)
}

/// The cell-quality check and remediation of `fl` (D50) for the log and
/// the warnings: nothing on a mesh without poor faces; else a MESH001
/// warning with the counts and how to avoid them, and with changed schemes a
/// MESH002 note. Returns the report for [`Solve::remediation`].
pub(crate) fn report_remediation<F: FnMut(&str)>(
    fl: &Flow,
    lg: &mut Log<F>,
    warnings: &mut Vec<diag::Diagnostic>,
) -> Option<crate::simplec::RemediationReport> {
    let r = &fl.remediation;
    if r.skewed_faces == 0 {
        return None;
    }
    let m = format!(
        "mesh quality: {} faces with skewness above {} (largest {:.3}) in {} cells; \
         an optimised mesh avoids them (Gmsh: Mesh.OptimizeNetgen = 1)",
        r.skewed_faces, r.skewness_limit, r.max_skewness, r.poor_cells
    );
    lg.line(&m);
    if !warnings.iter().any(|w| w.code == diag::MESH_SKEWED) {
        warnings.push(diag::Diagnostic::new(diag::MESH_SKEWED, m));
    }
    if !r.changes.is_empty() {
        let m = format!(
            "cell-quality remediation (D50, solver.remediation = \"auto\"): {}",
            r.changes.join(", ")
        );
        lg.line(&m);
        if !warnings.iter().any(|w| w.code == diag::MESH_REMEDIATION) {
            let mut d = diag::Diagnostic::new(diag::MESH_REMEDIATION, m);
            d.key = Some("solver.remediation".into());
            warnings.push(d);
        }
    }
    Some(r.clone())
}

/// Report of a flow case.
#[derive(Debug, Clone, Serialize)]
pub struct FlowReport {
    /// Case name.
    pub name: String,
    /// Reynolds number ρ U L / μ when the case defines one (lid speed or
    /// analytic velocity scale and height), else `None`.
    pub reynolds: Option<f64>,
    /// Findings of the boundary-condition check.
    pub findings: Vec<bc::Finding>,
    /// Solves per level.
    pub solves: Vec<Solve>,
    /// Observed orders of the sweep.
    pub orders: Vec<f64>,
    /// Observed time orders of a time-step sweep.
    pub time_orders: Vec<f64>,
    /// Observed orders of the scalar errors of a mesh sweep, per scalar
    /// with `exact` (D51).
    pub scalar_orders: Vec<(String, Vec<f64>)>,
    /// Discretisation uncertainty of the monitored quantities of a sweep
    /// over three or more meshes or time steps (D67).
    pub uncertainty: Vec<crate::uncertainty::Estimate>,
    /// What the sweep of `uncertainty` varies: `meshes` or `time steps`.
    pub uncertainty_of: &'static str,
    /// The run was stopped early: `.stop`, a signal or a limit.
    pub stopped: bool,
    /// Why it stopped early.
    pub stop: Option<crate::engine::Stop>,
    /// Warnings of the run with their codes (re-reads, restarts).
    pub warnings: Vec<diag::Diagnostic>,
    /// Problems found by the checks.
    pub failures: Vec<String>,
    /// Seconds.
    pub seconds: f64,
    /// The pressure convention of the case file and the operating pressure
    /// of every fluid part (D239), for people and agents.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pressure_values: Option<serde_json::Value>,
    /// The start values of the regions and fluid parts and where they come
    /// from (D238, REG010); empty when the case starts from its uniform
    /// `[initial]` values or restarts.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub start_values: Vec<crate::start::StartValue>,
}

pub(crate) fn mesh_for(case: &FlowCase, level: Option<u32>) -> Result<(Mesh, Geometry, String)> {
    let mut sec = case.mesh.clone();
    if let (Some(l), Some(g)) = (level, sec.generate.as_mut()) {
        g.refine = l;
    }
    // periodic pairs (D29) between the opposite sides of a generated box:
    // the perturbation keeps their nodes translated copies
    let pairs = bc::periodic_pairs(&case.boundary)?;
    if let Some(g) = sec.generate.as_mut() {
        for p in &pairs {
            for (d, ends) in [("xmin", "xmax"), ("ymin", "ymax"), ("zmin", "zmax")]
                .into_iter()
                .enumerate()
            {
                let (a, b) = (p.spec.a.as_str(), p.spec.b.as_str());
                if (a, b) == ends || (b, a) == ends {
                    g.periodic[d] = true;
                }
            }
        }
    }
    if let (Some(l), Some(g)) = (level, sec.cylinder.as_mut()) {
        g.refine = l;
    }
    if let (Some(l), Some(g)) = (level, sec.annulus.as_mut()) {
        g.refine = l;
    }
    if let (Some(l), Some(g)) = (level, sec.revolved.as_mut()) {
        g.refine = l;
    }
    let (mut mesh, source, _) = build_mesh(&sec)?;
    renumber(&mut mesh, sec.renumber, sec.seed).map_err(anyhow::Error::msg)?;
    // periodic pairs (D29): the last step, after the renumbering
    let mut desc = format!("{source}: {} cells", mesh.n_cells);
    if !pairs.is_empty() {
        let specs: Vec<_> = pairs.iter().map(|p| p.spec.clone()).collect();
        mesh = foehncraft_mesh::periodic::pair(&mesh, &specs)
            .map_err(|e| diag::err(diag::Diagnostic::new(diag::BC_PERIODIC, e)))?;
        let pd = mesh.periodic.as_ref().expect("paired");
        for p in &pd.pairs {
            let t = p.translation;
            let _ = write!(
                desc,
                "; periodic {}–{} ({} face pairs, translation ({:.6e}, {:.6e}, {:.6e}) m)",
                p.a,
                p.b,
                p.a_faces.len(),
                t[0],
                t[1],
                t[2]
            );
        }
    }
    let geo = Geometry::compute(&mesh);
    Ok((mesh, geo, desc))
}

pub(crate) fn dimension(case: &FlowCase) -> u8 {
    if case.mesh.cylinder.is_some() {
        return 2;
    }
    if let Some(a) = &case.mesh.annulus {
        return a.dimension();
    }
    case.mesh.generate.as_ref().map_or(3, |g| g.dimension)
}

/// Writes log lines to the terminal callback and to `flow.log`.
pub(crate) struct Log<'a, F: FnMut(&str)> {
    pub(crate) out: &'a mut F,
    pub(crate) text: String,
    pub(crate) style: crate::console::Style,
    /// The frames printed since the last [`sync`](Self::sync), as rows of
    /// the run record.
    pending: Vec<serde_json::Value>,
}

impl<'a, F: FnMut(&str)> Log<'a, F> {
    pub(crate) fn new(out: &'a mut F, style: crate::console::Style) -> Self {
        Self {
            out,
            text: String::new(),
            style,
            pending: Vec::new(),
        }
    }

    pub(crate) fn line(&mut self, l: &str) {
        (self.out)(l);
        self.text.push_str(l);
        self.text.push('\n');
    }

    /// An information line, ` %FOEHN-I-<FILE>, ...`.
    #[track_caller]
    pub(crate) fn note(&mut self, msg: &str) {
        self.line(&crate::console::tagged("note", msg));
    }

    /// A warning line, ` %FOEHN-W-<FILE>, ...`.
    #[track_caller]
    pub(crate) fn warn(&mut self, msg: &str) {
        self.line(&crate::console::tagged("warning", msg));
    }

    /// A blank line and a frame.
    pub(crate) fn frame(&mut self, f: &crate::console::Frame) {
        self.pending.push(f.record());
        self.line("");
        for l in f.render(self.style) {
            self.line(&l);
        }
    }
}

impl<F: FnMut(&str)> Log<'_, F> {
    /// Hands the frames printed so far to the run record as `section`
    /// rows, so that a GUI shows the same blocks as the console.
    pub(crate) fn sync(&mut self, rec: &mut crate::engine::Record) {
        for v in self.pending.drain(..) {
            rec.row(&v);
        }
    }
}

/// Velocity and pressure at `x` in cell `c`, reconstructed with the cell
/// gradients (second order).
pub(crate) fn value_at(
    fl: &crate::simplec::State,
    geo: &Geometry,
    c: usize,
    x: [f64; 3],
) -> ([f64; 3], f64) {
    let r = sub(x, geo.cell_centroid[c]);
    let u = fl.u.get(c);
    let g = [
        fl.grad_u[0].get(c),
        fl.grad_u[1].get(c),
        fl.grad_u[2].get(c),
    ];
    (
        [
            u[0] + dot(g[0], r),
            u[1] + dot(g[1], r),
            u[2] + dot(g[2], r),
        ],
        fl.p[c] + dot(fl.grad_p.get(c), r),
    )
}

/// The warnings and failed checks of a run in one frame, at its end.
pub(crate) fn end_frames<F: FnMut(&str)>(lg: &mut Log<F>, rep: &FlowReport) {
    let mut rows: Vec<String> = Vec::new();
    for f in &rep.findings {
        if f.level == Level::Warning {
            rows.push(finding_line(f).trim_start().to_string());
        }
    }
    for d in &rep.warnings {
        rows.push(format!("[W] {}: {}", d.code, d.message));
    }
    for f in &rep.failures {
        rows.push(format!("[E] FAILED: {f}"));
    }
    let mut fr = crate::console::Frame::new(format!("Warnings and failed checks ({})", rows.len()))
        .summary();
    if rows.is_empty() {
        fr.push_text("none");
    }
    for r in rows {
        fr.push_text(r);
    }
    lg.frame(&fr);
}

/// A boundary-condition finding as a tagged message line.
pub(crate) fn finding_line(f: &bc::Finding) -> String {
    let t = f.to_string();
    let rest = t.split_once(' ').map_or(t.as_str(), |x| x.1);
    crate::console::tagged(
        match f.level {
            Level::Note => "note",
            Level::Warning => "warning",
            Level::Error => "error",
        },
        rest,
    )
}

pub(crate) fn locate(loc: &Locator, geo: &Geometry, x: [f64; 3]) -> usize {
    loc.locate_search(x).unwrap_or_else(|| {
        (0..geo.cell_volume.len())
            .min_by(|&a, &b| {
                let da = norm(sub(geo.cell_centroid[a], x));
                let db = norm(sub(geo.cell_centroid[b], x));
                da.total_cmp(&db)
            })
            .unwrap_or(0)
    })
}

/// Runs a case; `log` receives progress lines; files go to `out`.
pub fn run(case: &FlowCase, name: &str, out: &Path, log: impl FnMut(&str)) -> Result<FlowReport> {
    let mut steer = Steer::quiet(out);
    run_steered(case, name, out, &mut steer, log)
}

/// The steady restart file of the current state.
fn steady_restart(fl: &mut Flow, mesh: &Mesh, name: &str) -> Result<Restart> {
    Ok(Restart {
        case: name.into(),
        cell_id: mesh.cell_id[..mesh.n_cells].to_vec(),
        face_id: mesh.face_id.clone(),
        state: fl.state()?,
        iteration: fl.iteration(),
        amg_setups: fl.amg_setups,
        amg_matrix: fl.amg_matrix.clone(),
        sweep: fl.sweep.clone(),
        transient: None,
        geometry: None,
        cht: crate::cht::restart_part(fl, None, None)?,
        spray: None,
    })
}

/// The physics of a case without a mesh: the density, the reference
/// viscosity and the uniform body force (zero when it is given by
/// expressions, which become a momentum source).
pub fn physics(case: &FlowCase) -> Result<Physics> {
    let p = case.properties()?;
    Ok(Physics::constant(
        p.density(),
        p.reference_viscosity(),
        case.uniform_body_force().unwrap_or([0.0; 3]),
        case.flow.convection,
    ))
}

/// `phys` on a mesh: the viscosity per face of the property program
/// (D94), when it is not one number.
pub fn physics_on(case: &FlowCase, phys: &Physics, mesh: &Mesh, geo: &Geometry) -> Result<Physics> {
    for (name, r) in &case.regions {
        for z in &r.zones {
            if !mesh.zones.contains(z) {
                let hint = match bc::did_you_mean(z, mesh.zones.iter().map(String::as_str)) {
                    Some(c) => format!(" (did you mean `{c}`?)"),
                    None => String::new(),
                };
                return Err(diag::err(
                    diag::Diagnostic::new(
                        diag::REG_ZONE,
                        format!(
                            "[regions.{name}] zone `{z}` is not a cell zone of the mesh{hint}; \
                             the mesh has {}",
                            mesh.zones
                                .iter()
                                .map(|z| format!("`{z}`"))
                                .collect::<Vec<_>>()
                                .join(", ")
                        ),
                    )
                    .key(format!("regions.{name}.zones")),
                ));
            }
        }
    }
    let mut p = phys.clone();
    let props = case.properties()?;
    p.mu_face = props.face_viscosity(mesh, geo);
    // several fluids (D194): the density of each cell and face
    if let Some((cell, face)) = props.cell_face_density(mesh) {
        p.rho_var = Some(std::sync::Arc::new((cell, face)));
        p.zones = true;
    }
    Ok(p)
}

/// The momentum source per unit volume (N/m³) of a body force given by
/// expressions (D53), at the cell centroids at time `t`; `None` for a
/// uniform body force.
pub fn momentum_source(
    case: &FlowCase,
    geo: &Geometry,
    n_cells: usize,
    t: f64,
) -> Result<Option<[Vec<f64>; 3]>> {
    if case.uniform_body_force().is_some() {
        return Ok(None);
    }
    let v = case.body_force_values()?;
    let out: [Vec<f64>; 3] = std::array::from_fn(|i| {
        (0..n_cells)
            .map(|c| v[i].at(geo.cell_centroid[c], t))
            .collect()
    });
    if let Some((i, c)) = out
        .iter()
        .enumerate()
        .find_map(|(i, s)| s.iter().position(|x| !x.is_finite()).map(|c| (i, c)))
    {
        return Err(diag::err(
            diag::Diagnostic::new(
                diag::EXPR_VALUE,
                format!(
                    "flow.body_force[{i}] is not finite at the centroid {:?} of cell {c} (t = {t} s)",
                    geo.cell_centroid[c]
                ),
            )
            .key(format!("flow.body_force[{i}]")),
        ));
    }
    Ok(Some(out))
}

/// True when the body force or a boundary value depends on t.
pub fn time_dependent(case: &FlowCase, values: &bc::Compiled) -> Result<(bool, bool)> {
    let source = case.body_force_values()?.iter().any(|v| v.uses().time());
    let boundary = values.values().any(|v| v.uses().time());
    Ok((source, boundary))
}

/// The names of the residuals of a steady run, in the order of
/// [`residual_values`] (D198 (d)).
fn residual_names(case: &FlowCase) -> Vec<String> {
    let mut v: Vec<String> = ["momentum_x", "momentum_y", "momentum_z", "continuity"]
        .iter()
        .map(|s| (*s).to_string())
        .collect();
    if !case.turbulence.is_laminar() {
        v.extend(
            crate::turbulence::equations(case.turbulence.model)
                .iter()
                .map(|s| (*s).to_string()),
        );
    }
    v.extend(case.scalars.keys().map(|s| format!("scalar_{s}")));
    if case.energy.is_some() {
        v.push("energy".into());
    }
    v
}

/// The residuals of one iteration in the order of [`residual_names`].
fn residual_values(r: &Residuals) -> Vec<f64> {
    let mut v = r.momentum.to_vec();
    v.push(r.continuity);
    v.extend(&r.turbulence.residual);
    v.extend(&r.scalars.residual);
    v
}

/// The log frame of a diagnosis (D198 (d)).
fn diagnosis_frame(d: &crate::budget::Diagnosis) -> crate::console::Frame {
    let verdict = serde_json::to_value(d.verdict)
        .ok()
        .and_then(|v| v.as_str().map(str::to_uppercase))
        .unwrap_or_default();
    let mut f = crate::console::Frame::new(format!("Diagnosis at iteration {}", d.iteration))
        .kv("Verdict", verdict)
        .kv("Budget", format!("{} iterations", d.budget))
        .kv("Controlling", &d.controlling)
        .kv("Window", format!("last {} iterations", d.window));
    for e in &d.equations {
        let trend = serde_json::to_value(e.trend)
            .ok()
            .and_then(|v| v.as_str().map(str::to_string))
            .unwrap_or_default();
        f.push_kv(
            &e.equation,
            format!(
                "{:.2e}, {trend}, {:+.2} decades per 1000 iterations{}",
                e.residual,
                e.decades_per_1000,
                e.tolerance_at
                    .map_or(String::new(), |n| format!(", tolerance at about {n}"))
            ),
        );
    }
    f
}

/// As [`run`], with the requests, limits and run record of `steer`
/// (crate::engine).
pub fn run_steered(
    case: &FlowCase,
    name: &str,
    out: &Path,
    steer: &mut Steer,
    log: impl FnMut(&str),
) -> Result<FlowReport> {
    if case.time.is_some() {
        return crate::unsteady::run(case, name, out, steer, log);
    }
    let t0 = Instant::now();
    std::fs::create_dir_all(out)?;
    let control = Control::new(out)?;
    let mut case = case.clone();
    let mut log = log;
    let mut lg = Log::new(&mut log, case.log.style);
    let levels: Vec<Option<u32>> = if case.verification.levels.is_empty() {
        vec![None]
    } else {
        case.verification.levels.iter().map(|&l| Some(l)).collect()
    };
    let analytic = case.analytic();
    let phys = physics(&case)?;
    let values = case.compiled_boundaries()?;
    let at = bc::At {
        analytic: analytic.as_ref(),
        values: &values,
        t: 0.0,
        turning: &case.turning(),
        hydrostatic: crate::buoyancy::hydrostatic(&case),
    };
    let reynolds = case.reynolds();
    let mut rep = FlowReport {
        name: name.into(),
        reynolds,
        findings: Vec::new(),
        solves: Vec::new(),
        orders: Vec::new(),
        time_orders: Vec::new(),
        scalar_orders: Vec::new(),
        uncertainty: Vec::new(),
        uncertainty_of: "",
        stopped: false,
        stop: None,
        warnings: Vec::new(),
        failures: Vec::new(),
        seconds: 0.0,
        pressure_values: case.pressure_frame.as_ref().map(|f| f.json()),
        start_values: Vec::new(),
    };
    std::fs::write(out.join("case.resolved.toml"), case.resolved_toml()?)?;
    if !case.run.ranks.is_empty() {
        crate::ranks::supported(&case)?;
    }
    let last = levels.len() - 1;
    let single = levels.len() == 1;
    if !single && case.run.restart.is_some() {
        bail!("[run] restart continues a single run, not a sweep over levels");
    }
    let mut histories: Vec<(String, Vec<(f64, f64)>)> = Vec::new();
    let mut stopped = false;
    let case_at_start = (!case.schedule.is_empty()).then(|| case.clone());
    for (li, level) in levels.iter().enumerate() {
        // every level of a sweep runs the whole schedule from the start
        if li > 0
            && let Some(c) = &case_at_start
        {
            case = c.clone();
        }
        let (mesh, geo, desc) = mesh_for(&case, *level)?;
        let mesh_hash = crate::engine::mesh_hash(&mesh);
        // conjugate heat transfer (stage 7d1): the flow runs on the fluid
        // view, the energy equation on the full mesh
        let (fluid_case, mesh, geo) = crate::cht::view(&case, mesh, geo)?;
        case = fluid_case;
        // the fluid parts (D194): REG004, REG003
        let parts = crate::fluids::check(&case, &mesh)?;
        let nc = mesh.n_cells;
        steer.set_cells(nc, usize::from(case.run.backend != "cpu"));
        if li == 0 {
            lg.frame(&crate::logblocks::banner(steer.build.as_ref(), &case));
            for n in std::mem::take(&mut steer.notes) {
                lg.note(&n);
            }
        }
        // the density and viscosity per cell and face (D94, D194)
        let phys_on = physics_on(&case, &phys, &mesh, &geo)?;
        // boundary-condition check (D30)
        let mut findings = bc::check(
            &mesh,
            &geo,
            &case.boundary,
            &phys_on,
            at,
            case.check_context(),
        );
        findings.extend(crate::turbulence::setup_findings(&case));
        findings.extend(crate::energy::setup_findings(&case));
        findings.extend(crate::pressure_values::setup_findings(&case));
        findings.extend(crate::buoyancy::setup_findings(&case, &mesh, &geo));
        if li == 0 {
            rep.findings = findings.clone();
        }
        if findings.iter().any(|f| f.level == Level::Error) {
            let msg: Vec<String> = findings
                .iter()
                .filter(|f| f.level == Level::Error)
                .map(|f| f.to_string())
                .collect();
            lg.line("");
            for f in &findings {
                lg.line(&finding_line(f));
            }
            std::fs::write(out.join("flow.log"), &lg.text)?;
            bail!(
                "the boundary conditions cannot be run:\n  {}",
                msg.join("\n  ")
            );
        }
        let (quality, _) = foehncraft_mesh::quality::QualityReport::compute(&mesh, &geo);
        if li == 0 {
            if single {
                lg.frame(&files_frame(&case, &mesh, &desc, out));
            }
            lg.frame(&crate::logblocks::control_steady(&case));
        }
        lg.frame(&crate::logblocks::mesh(&quality, &desc));
        if li == 0 {
            lg.frame(&crate::logblocks::fluid(&case));
            lg.frame(&crate::logblocks::solution(&case, false));
            lg.frame(&crate::logblocks::boundaries(&case, &quality));
            if let Some(f) = crate::logblocks::limitations(&case) {
                lg.frame(&f);
            }
            let mut any = false;
            for n in crate::logblocks::notes(&case)
                .into_iter()
                .chain(crate::fluids::notes(&parts))
            {
                if !any {
                    lg.line("");
                    any = true;
                }
                lg.line(&n);
            }
            for f in &findings {
                if !any {
                    lg.line("");
                    any = true;
                }
                lg.line(&finding_line(f));
            }
        }
        if !single && li == last {
            lg.frame(&files_frame(&case, &mesh, &desc, out));
        }
        let arrays = bc::arrays(&mesh, &geo, &case.boundary, &phys_on, at)?;
        let disc = simplec::discretisation(&mesh, &geo, &arrays)?;
        let operators = crate::operators::report(
            &case,
            &mesh,
            &geo,
            &arrays.u_points,
            out,
            &mut lg,
            &mut rep.warnings,
        );
        let phys = phys_on;
        let loc = Locator::new(&mesh, &geo);
        let reference_cell = if arrays.has_pressure() {
            None
        } else {
            let x = case.pressure.reference_point.unwrap_or_else(|| {
                let (lo, hi) = mesh.bounding_box();
                std::array::from_fn(|i| lo[i] + 1e-6 * (hi[i] - lo[i]))
            });
            Some(locate(&loc, &geo, x))
        };
        let p0 = if reference_cell.is_some() {
            case.pressure.reference_value
        } else {
            case.initial.pressure
        };
        let mut fl = Flow::new(
            disc,
            arrays,
            phys,
            case.solver.clone(),
            case.initial.velocity,
            p0,
            reference_cell,
            &crate::simplec::Target {
                device: case.run.backend.clone(),
                threads: case.run.threads,
                ptx: case.run.ptx.clone(),
                ranks: None,
            },
        )?;
        let remediation = report_remediation(&fl, &mut lg, &mut rep.warnings);
        if let Some(src) = momentum_source(&case, &geo, nc, 0.0)? {
            fl.set_momentum_source(src)?;
        }
        crate::scalar::attach(&mut fl, &case, &mesh, &geo)?;
        crate::density::attach(&mut fl, &case)?;
        crate::buoyancy::attach(&mut fl, &case)?;
        crate::porous::attach(&mut fl, &case, &mesh)?;
        crate::turbulence::attach(&mut fl, &case, &mesh, &geo)?;
        if let Some(notes) = crate::frame::attach(&mut fl, &case, &mesh, &geo)? {
            for l in &notes.lines {
                lg.line(l);
            }
            rep.warnings.extend(notes.diagnostics);
        }
        if let Some(path) = case.run.restart.clone() {
            let r = restart::read(&path)?;
            let (r, same) =
                r.to_local(&mesh.cell_id[..nc], &mesh.face_id, mesh.n_interior_faces())?;
            if r.transient.is_some() {
                bail!(
                    "{} is the restart file of a transient run; this case is steady (no [time] table)",
                    path.display()
                );
            }
            fl.load(&r.state)?;
            crate::cht::load_restart(&mut fl, r.cht.as_ref(), &path)?;
            fl.iteration = r.iteration;
            if let Some(m) = &r.amg_matrix {
                fl.restore_amg(m, r.amg_setups)?;
            }
            if case.solver.algorithm(false) == simplec::Algorithm::Coupled
                && case.solver.coupled.momentum_sweep == crate::coupled::MomentumSweep::Multigrid
            {
                fl.restore_sweep(&r.sweep)?;
            }
            if !same {
                let m = "note: the restart file has another cell or face order; the fields were \
                         mapped by identity and the multigrid rebuilt (the continuation is not \
                         bit-identical)";
                lg.note(m.trim_start_matches("note: "));
                rep.warnings
                    .push(diag::Diagnostic::new(diag::RUN_RESTART_MAPPED, m));
            }
            lg.note(&format!(
                "restart from {} at iteration {}",
                path.display(),
                r.iteration
            ));
        } else {
            // the start values of the regions and parts (D238)
            let mut starts = Vec::new();
            if let Some((u, p, lines)) = crate::start::initial_flow(
                &case,
                &mesh,
                &geo,
                case.initial.velocity,
                p0,
                reference_cell.is_some(),
            )? {
                fl.initialise(&u, &p)?;
                starts = lines;
            }
            starts.extend(crate::start::temperature_report(&case, &mesh, &geo)?);
            starts.extend(crate::start::turbulence_report(&case));
            if li == 0 && !starts.is_empty() {
                lg.line("");
                for n in crate::start::notes(&starts) {
                    lg.line(&n);
                }
                rep.start_values = starts;
            }
        }
        if case.run.restart.is_none() && case.initial.potential_flow {
            let it = fl.potential_flow(&mesh, &geo)?;
            lg.note(&format!(
                "start from potential flow (D52): {it} CG iterations"
            ));
        }
        let first = fl.iteration();
        let mut units = String::from(
            "iteration -, momentum_x -, momentum_y -, momentum_z -, continuity - (scaled residuals), \
             bicgstab_x -, bicgstab_y -, bicgstab_z -, cg - (iterations)",
        );
        let mut head = String::from(
            "iteration,momentum_x,momentum_y,momentum_z,continuity,bicgstab_x,bicgstab_y,bicgstab_z,cg",
        );
        if !case.turbulence.is_laminar() {
            for s in crate::turbulence::equations(case.turbulence.model) {
                let _ = write!(
                    units,
                    ", {s} - (scaled residual), bicgstab_{s} - (iterations)"
                );
                let _ = write!(head, ",{s},bicgstab_{s}");
            }
        }
        for s in case.scalars.keys() {
            let _ = write!(
                units,
                ", scalar_{s} - (scaled residual), bicgstab_{s} - (iterations)"
            );
            let _ = write!(head, ",scalar_{s},bicgstab_{s}");
        }
        // the energy equation, after the passive scalars (stage 7a1)
        if case.energy.is_some() {
            units.push_str(", energy - (scaled residual), bicgstab_energy - (iterations)");
            head.push_str(",energy,bicgstab_energy");
        }
        // the Krylov iterations of the coupled algorithm: a last column
        let coupled = case.solver.algorithm(false) == crate::simplec::Algorithm::Coupled;
        if coupled {
            units.push_str(", krylov - (iterations of the coupled system)");
            head.push_str(",krylov");
        }
        let mut csv = csv_preamble("residuals", name, None, &units);
        csv.push_str(&head);
        csv.push('\n');
        let (bunits, bhead) = balance::csv_head(false, !case.turning().is_empty());
        let mut bcsv = csv_preamble("balance", name, None, &bunits);
        bcsv.push_str(&bhead);
        bcsv.push('\n');
        let nonorth = case.solver.nonorth_correction;
        let mut hist = Vec::new();
        let mut converged = false;
        let mut res = Residuals::default();
        let ts = Instant::now();
        let mut cols = crate::logcols::Columns::build(&case, &mesh, &geo, &loc, false)?;
        let mut table = cols.table("Iter", case.log.style, case.log.header_every);
        if steer.record.on() {
            steer.record.row(&cols.record_header("iteration", "Iter"));
        }
        lg.frame(&crate::console::Frame::new(format!(
            "Iterations ({}, steady)",
            crate::logblocks::spelled(&case.solver.algorithm(false)).to_uppercase()
        )));
        lg.line("");
        if !case.log.detail.rows() {
            lg.note(
                "the table is not shown ([log] detail = \"summary\"); run.jsonl holds every row",
            );
        }
        lg.sync(&mut steer.record);
        let mut sched = crate::schedule::Runner::default();
        // the iteration budget and its diagnosis (D198 (d)): reads the
        // residuals only
        let budget =
            crate::budget::IterationBudget::of(case.solver.max_iterations, &case.run.limits);
        let mut diag_hist = crate::budget::History::new(
            residual_names(&case),
            budget.diagnose_at.filter(|&n| n > first),
        );
        loop {
            if !case.schedule.is_empty() {
                let restarted = case.run.restart.is_some();
                let p = crate::schedule::Position {
                    n: fl.iteration() + 1,
                    time: 0.0,
                    tol: 0.0,
                };
                schedule_step(&mut case, &mut fl, &mut lg, steer, &mut sched, p, restarted)?;
            }
            if fl.iteration() >= case.solver.max_iterations {
                break;
            }
            let previous = res;
            fl.keep_valid()?;
            res = match fl.iterate() {
                Ok(r) => r,
                Err(e) if divergence::is_divergence(&e, &mut fl) => {
                    // D63: the last valid state, the snapshot and the report
                    lg.line(&crate::console::tagged("error", &format!("{e:#}")));
                    let at = fl.iteration();
                    let ctx = divergence::Context {
                        name,
                        case: &case,
                        at,
                        done: at.saturating_sub(first + 1),
                        time: None,
                        dt: None,
                        residuals: fl.diverged.clone().unwrap_or_default(),
                        last_residuals: (at > first + 1).then_some(previous.clone()),
                        first_courant: None,
                    };
                    let prepared = (|| -> Result<(simplec::State, Vec<&'static str>)> {
                        let failed = fl.state()?;
                        anyhow::ensure!(fl.rewind()?, "no valid state was kept");
                        let mut files = Vec::new();
                        if cfg!(feature = "hdf5") && case.output.vtkhdf {
                            let sol = Solution {
                                state: fl.state()?,
                                disc: &fl.disc,
                                phys: fl.physics(),
                            };
                            write_fields(
                                &case,
                                &out.join("fields_last_valid.vtkhdf"),
                                &mesh,
                                &geo,
                                &sol,
                                analytic.as_ref(),
                                Some((fl.iteration(), None)),
                            )?;
                            files.push("fields_last_valid.vtkhdf");
                        }
                        if single && restart::AVAILABLE {
                            let rs = steady_restart(&mut fl, &mesh, name)?;
                            restart::write(&out.join("restart.h5"), &rs)?;
                            files.push("restart.h5");
                        }
                        Ok((failed, files))
                    })();
                    let e = divergence::finish(
                        &mut fl,
                        &ctx,
                        &mesh,
                        &geo,
                        out,
                        prepared,
                        case.log.style,
                        &mut |l| lg.line(l),
                    );
                    lg.sync(&mut steer.record);
                    steer.record.flush();
                    std::fs::write(out.join("flow.log"), &lg.text)?;
                    return Err(e);
                }
                Err(e) => return Err(e),
            };
            if let Some(m) = fl.sweep_note.take() {
                lg.warn(&m);
            }
            if let Some(m) = fl.turbulence.as_ref().and_then(|t| t.take_note()) {
                lg.warn(&m);
            }
            let it = fl.iteration();
            let done = res.max() < case.solver.tolerance;
            diag_hist.push(it, &residual_values(&res));
            if diag_hist.due(it)
                && !done
                && let Some(d) =
                    crate::budget::diagnose(&diag_hist, case.solver.tolerance, budget.iterations)
            {
                lg.frame(&diagnosis_frame(&d));
                if d.verdict.off_track() {
                    lg.warn(&format!("diagnosis at iteration {it}: {}", d.advice));
                } else {
                    lg.note(&format!("diagnosis at iteration {it}: {}", d.advice));
                }
                steer.diagnosis(&d, case.run.limits.stops_at_diagnosis());
            }
            if it % case.output.history_every == 0 || done {
                let _ = write!(
                    csv,
                    "{},{:e},{:e},{:e},{:e},{},{},{},{}",
                    it,
                    res.momentum[0],
                    res.momentum[1],
                    res.momentum[2],
                    res.continuity,
                    res.momentum_iterations[0],
                    res.momentum_iterations[1],
                    res.momentum_iterations[2],
                    res.pressure_iterations
                );
                for (r, n) in res
                    .turbulence
                    .residual
                    .iter()
                    .zip(&res.turbulence.iterations)
                    .chain(res.scalars.residual.iter().zip(&res.scalars.iterations))
                {
                    let _ = write!(csv, ",{r:e},{n}");
                }
                if coupled {
                    let _ = write!(csv, ",{}", res.coupled_iterations);
                }
                csv.push('\n');
                hist.push((it as f64, res.max()));
                if steer.record.on() {
                    let mut row = serde_json::json!({
                        "type": "iteration",
                        "iteration": it,
                        "momentum": res.momentum,
                        "continuity": res.continuity,
                        "bicgstab": res.momentum_iterations,
                        "cg": res.pressure_iterations,
                        "scalars": res.scalars.residual,
                        "turbulence": res.turbulence.residual,
                        "turbulence_floors": res.turbulence.floors,
                    });
                    if coupled {
                        row["krylov"] = serde_json::json!(res.coupled_iterations);
                        if let Some(dt) = res.pseudo_step {
                            row["pseudo_step"] = serde_json::json!(dt);
                        }
                    }
                    steer.record.row(&row);
                }
            }
            if it == first + 1 || it % case.solver.log_every == 0 || done {
                lg.sync(&mut steer.record);
                steer.record.flush();
                // `log.columns` changed during the run (D214): new columns,
                // the header row again
                crate::logcols::refresh(
                    &mut cols,
                    &mut table,
                    &mut crate::logcols::Refresh {
                        case: &case,
                        mesh: &mesh,
                        geo: &geo,
                        loc: &loc,
                        transient: false,
                        table_name: "iteration",
                        first: "Iter",
                        at: it,
                    },
                    &mut |m| lg.warn(m),
                    steer,
                );
                let sol = if cols.needs_state() {
                    Some(Solution {
                        state: fl.state()?,
                        disc: &fl.disc,
                        phys: fl.physics(),
                    })
                } else {
                    None
                };
                let rotors = (cols.needs_rotors() && sol.is_some())
                    .then(|| crate::rotor::evaluate(sol.as_ref().unwrap(), &mesh, &geo, &case));
                let entries = cols.eval(
                    &crate::logcols::RowCtx {
                        index: it,
                        res: &res,
                        time: None,
                        sol: sol.as_ref(),
                        rotors: rotors.as_ref(),
                    },
                    &geo,
                );
                let cells: Vec<String> = entries.iter().map(|e| e.text.clone()).collect();
                if case.log.detail.rows() {
                    for l in table.row(it, &cells) {
                        lg.line(&l);
                    }
                }
                if steer.record.on() {
                    steer
                        .record
                        .row(&cols.record_row("iteration", it, &entries));
                }
            }
            if done {
                converged = true;
                break;
            }
            // the ledger of the current state (D88, D89), with its boundary
            // values and gradients; the run goes on from those it had (D31)
            if case.balance.every > 0 && it.is_multiple_of(case.balance.every) {
                let lagged = fl.lagged()?;
                fl.update_boundaries()?;
                fl.gradients()?;
                let sol = Solution {
                    state: fl.state()?,
                    disc: &fl.disc,
                    phys: fl.physics(),
                };
                let r = balance::rates(
                    &sol,
                    &mesh,
                    &geo,
                    &case,
                    &fl.bc.mom_kind,
                    &geo.cell_volume,
                    nonorth,
                );
                fl.set_lagged(&lagged)?;
                let b = balance::steady(&r, it);
                bcsv.push_str(&balance::csv_row(&b));
                bcsv.push('\n');
                if steer.record.on() {
                    steer.record.row(&balance::record_row(&b));
                }
                if case.log.detail.rows() {
                    lg.frame(&crate::logblocks::balance_block(
                        &format!("Balance of iteration {it}"),
                        &b,
                        false,
                    ));
                }
            }
            // control files (D40), signals and limits (D64, D83)
            let req = control.check(it, case.run.reread_every);
            if req.result && li == last {
                // as for the ledger: the run goes on from the boundary
                // values and gradients it had (D31)
                let lagged = fl.lagged()?;
                fl.update_boundaries()?;
                fl.gradients()?;
                let sol = Solution {
                    state: fl.state()?,
                    disc: &fl.disc,
                    phys: fl.physics(),
                };
                write_fields(
                    &case,
                    &out.join(format!("fields_{it:06}.vtkhdf")),
                    &mesh,
                    &geo,
                    &sol,
                    analytic.as_ref(),
                    Some((it, None)),
                )?;
                fl.set_lagged(&lagged)?;
                lg.note(&format!(
                    "fields of iteration {it} written (fields_{it:06}.vtkhdf)"
                ));
            }
            let r = &case.run;
            let restart_now = single
                && ((r.restart_every > 0 && it.is_multiple_of(r.restart_every))
                    || r.restart_at.contains(&it)
                    || steer.restart_requested());
            if restart_now && restart::AVAILABLE {
                let rs = steady_restart(&mut fl, &mesh, name)?;
                restart::write(&out.join("restart.h5"), &rs)?;
                if case.run.restart_at.contains(&it) {
                    restart::write(&out.join(format!("restart_{it:06}.h5")), &rs)?;
                }
                lg.note(&format!("iteration {it}: restart.h5 written"));
            }
            if req.reread {
                reread_case(
                    &mut case,
                    &mut fl,
                    &mut lg,
                    steer,
                    &mut rep.warnings,
                    it,
                    out,
                    sched.fired,
                )?;
            }
            apply_live_settings(&mut case, &mut fl, &mut lg, steer, &mut rep.warnings, it)?;
            if let Some(s) = steer.poll(&case.run.limits, it - first, false, req.stop) {
                lg.note(&format!("iteration {it}: {}", s.message));
                stopped = true;
                break;
            }
        }
        steer.record.flush();
        let secs = ts.elapsed().as_secs_f64();
        let iters = fl.iteration();
        {
            let wall = ts.elapsed().as_secs_f64();
            let _ = secs;
            let mut f = crate::console::Frame::new("Convergence").summary().kv(
                "Status",
                if converged {
                    format!(
                        "CONVERGED in {iters} iterations (criterion {})",
                        crate::console::sci(case.solver.tolerance, 1)
                    )
                } else if stopped {
                    format!("STOPPED after {iters} iterations")
                } else {
                    format!("NOT CONVERGED in {iters} iterations")
                },
            );
            f.push_kv("Largest scaled residual", crate::console::sci3(res.max()));
            if let Some(dt) = res.pseudo_step {
                f.push_kv(
                    "Pseudo-time step",
                    format!(
                        "{} s  ({}, last iteration)",
                        crate::console::sci3(dt),
                        if case.solver.pseudo_time.step.is_some() {
                            "given"
                        } else {
                            "automatic"
                        }
                    ),
                );
            }
            f.push_kv(
                "Solve time",
                format!(
                    "{wall:.2} s  ({:.0} ns per cell and iteration)",
                    1e9 * fl.timings.total() / (nc as f64 * (iters - first).max(1) as f64)
                ),
            );
            lg.frame(&f);
        }
        let suffix = level.map_or(String::new(), |l| format!("_level{l}"));
        std::fs::write(out.join(format!("residuals{suffix}.csv")), &csv)?;
        histories.push((
            level.map_or("residual".into(), |l| format!("level {l}")),
            hist,
        ));
        if !converged && !stopped && case.verification.require_convergence.unwrap_or(true) {
            rep.failures.push(format!(
                "{desc}: not converged in {iters} outer iterations (largest residual {:.2e})",
                res.max()
            ));
        }
        // final gradients for the output
        fl.update_boundaries()?;
        fl.gradients()?;
        // conjugate heat transfer: the fields of the output on the full mesh
        let full_out = if case.fluid_view.is_some() {
            let st = fl.state()?;
            Some(crate::cht::full_state(&mut fl, &st)?)
        } else {
            None
        };
        let fl_host = Solution {
            state: fl.state()?,
            disc: &fl.disc,
            phys: fl.physics(),
        };
        let sol = &fl_host;
        let errors = analytic.map(|a| error_norms(sol, &mesh, &geo, &a, dimension(&case)));
        let mut verification = crate::console::Frame::new("Verification").summary();
        if let Some(e) = &errors {
            verification.push_kv(
                "Velocity error L2 / Linf",
                format!(
                    "{} / {}  (relative)",
                    crate::console::sci3(e.l2),
                    crate::console::sci3(e.linf)
                ),
            );
            verification.push_kv("Pressure error Linf", crate::console::sci3(e.pressure_linf));
        }
        let reference = case
            .verification
            .reference
            .map(|r| compare_reference(sol, &mesh, &geo, &loc, r, &case));
        if let Some(r) = &reference {
            verification.push_kv(
                &format!("Reference {}", r.reference),
                format!(
                    "largest deviation u {:.4}, v {:.4}  (lid speeds)",
                    r.u_max_dev, r.v_max_dev
                ),
            );
        }
        let boundaries = boundary_summaries(sol, &mesh, &geo, &case);
        let wall_probe_values = wall_probes(&case, &mesh, &geo, sol);
        for w in &wall_probe_values {
            lg.line(&wall_probe_line(w));
        }
        let wall_crossing_values = wall_crossings(&case, &mesh, &geo, sol);
        for (w, c) in wall_crossing_values.iter().zip(&case.output.wall_crossings) {
            lg.line(&wall_crossing_line(w, c));
        }
        let total: f64 = boundaries.iter().map(|p| p.mass_flow).exact_sum();
        // scale: the largest boundary mass flow, or in a closed domain the
        // flow ρ |u|max A through the largest boundary
        let umax = (0..nc).map(|c| norm(sol.u.get(c))).fold(0.0f64, f64::max);
        let big = boundaries
            .iter()
            .map(|p| p.mass_flow.abs().max(case.rho() * umax * p.area))
            .fold(0.0f64, f64::max);
        let mass_imbalance = if big > 0.0 { total / big } else { 0.0 };
        let coefficients = crate::monitor::coefficients(&case, &boundaries);
        if let Some([cd, cl]) = coefficients {
            verification.push_kv("Force coefficients", format!("C_D {cd:.5}, C_L {cl:.5}"));
        }
        let benchmark =
            crate::unsteady::steady_benchmark_of(&case, sol, &mesh, &geo, &loc, &boundaries);
        if let Some(b) = &benchmark {
            crate::unsteady::log_benchmark(b, |l| {
                verification.push_text(l.trim_start().to_string())
            });
        }
        if !verification.is_empty() {
            lg.frame(&verification);
        }
        lg.frame(&crate::logblocks::results(&boundaries));
        let rotor_values = rotor_report(sol, &mesh, &geo, &case, &mut lg, &mut rep.warnings);
        // the ledger of the final state (D88, D89, D90 part 4)
        let bal = balance::steady(
            &balance::rates(
                sol,
                &mesh,
                &geo,
                &case,
                &fl.bc.mom_kind,
                &geo.cell_volume,
                nonorth,
            ),
            iters,
        );
        lg.frame(&crate::logblocks::balance_block(
            "Balance of the final state",
            &bal,
            false,
        ));
        lg.frame(&crate::logblocks::run_time(
            &fl.timings,
            nc,
            iters - first,
            "iterations",
        ));
        bcsv.push_str(&balance::csv_row(&bal));
        bcsv.push('\n');
        if steer.record.on() {
            steer.record.row(&balance::record_row(&bal));
        }
        std::fs::write(out.join(format!("balance{suffix}.csv")), &bcsv)?;
        let mut summary = BalanceSummary {
            largest_mass_relative: bal.mass.relative.abs(),
            last: bal,
            ..Default::default()
        };
        for d in balance::check(&summary, &case.balance) {
            lg.warn(&format!("{}: {}", d.code, d.message));
            steer.event(d.code, iters, &d.message);
            rep.warnings.push(d);
        }
        if li == last || stopped {
            std::fs::write(
                out.join("balance.svg"),
                balance::chart(name, &summary.last, false),
            )?;
            match (&case.fluid_view, &full_out, fl.cht.as_deref()) {
                (Some(v), Some(st), Some(c)) => {
                    let full = Solution {
                        state: st.clone(),
                        disc: &c.fe.disc,
                        phys: c.fe.physics(),
                    };
                    write_outputs(
                        &v.full_case,
                        out,
                        &v.mesh,
                        &v.geo,
                        &full,
                        None,
                        &boundaries,
                        &Locator::new(&v.mesh, &v.geo),
                        Some((iters, None)),
                    )?;
                }
                _ => write_outputs(
                    &case,
                    out,
                    &mesh,
                    &geo,
                    sol,
                    analytic.as_ref(),
                    &boundaries,
                    &loc,
                    Some((iters, None)),
                )?,
            }
            if let Some(r) = &reference {
                write_reference_svg(out, name, r)?;
            }
        }
        drop(fl_host);
        let scalars = crate::scalar::summaries(&mut fl, &case, &mesh, &geo, dimension(&case), 0.0)?;
        for sm in &scalars {
            lg.line(&crate::scalar::summary_line(sm));
            if let Some(d) = converged
                .then(|| crate::scalar::heat_balance_finding(sm))
                .flatten()
            {
                lg.note(&d.message);
                rep.warnings.push(d);
            }
        }
        // the thermal account of the final state (D119, D161)
        if let Some(sm) = scalars.iter().find(|s| s.energy.is_some())
            && let Some(e) = &sm.energy
        {
            summary
                .last
                .book(balance::thermal_account(e, &sm.boundaries));
        }
        let periodic = periodic_summaries(&mut fl, &mesh, &mut lg)?;
        let turbulence = turbulence_summary(&mut fl, &mesh, &mut lg, &mut rep.warnings)?;
        let density = crate::density::summary(&mut fl, &mut rep.warnings)?;
        if single && (case.run.restart_final || stopped) && restart::AVAILABLE {
            let rs = steady_restart(&mut fl, &mesh, name)?;
            restart::write(&out.join("restart.h5"), &rs)?;
        }
        rep.solves.push(Solve {
            mesh: desc,
            mesh_hash,
            device: fl.device(),
            cells: nc,
            iterations: iters,
            converged,
            residuals: res,
            ns_per_cell_iteration: 1e9 * fl.timings.total()
                / (nc as f64 * (iters - first).max(1) as f64),
            timings: fl.timings.clone(),
            amg_setups: fl.amg_setups,
            sweep_setups: fl.sweep_setups,
            turbulence_setups: fl.turbulence.as_ref().map_or(0, |t| t.reduction_setups()),
            errors,
            reference,
            boundaries,
            mass_imbalance,
            coefficients,
            rotors: rotor_values,
            transient: None,
            datums: None,
            spray: None,
            benchmark,
            balance: Some(summary),
            remediation,
            operators,
            scalars,
            periodic,
            wall_probes: wall_probe_values,
            wall_crossings: wall_crossing_values,
            turbulence,
            statistics: None,
            density,
        });
        if stopped {
            rep.stopped = true;
            rep.stop = steer.stop.clone();
            break;
        }
    }
    // checks
    let v = &case.verification;
    for s in &rep.solves {
        if let (Some(tol), Some(b)) = (v.benchmark_tolerance, &s.benchmark)
            && b.max_deviation() > tol
        {
            rep.failures.push(format!(
                "{}: largest distance from the {} ranges {:.4} above {tol}",
                s.mesh,
                b.benchmark,
                b.max_deviation()
            ));
        }
    }
    crate::scalar::check_sweep(&mut rep, v, levels.len(), &mut |l| lg.line(l));
    if levels.len() >= 2 && rep.solves.iter().all(|s| s.errors.is_some()) {
        let e: Vec<&ErrorNorms> = rep
            .solves
            .iter()
            .filter_map(|s| s.errors.as_ref())
            .collect();
        rep.orders = e
            .windows(2)
            .map(|w| (w[0].l2 / w[1].l2).ln() / (w[0].h / w[1].h).ln())
            .collect();
        lg.line(&format!("observed orders {:.3?}", rep.orders));
        if let (Some([lo, hi]), Some(&o)) = (v.expected_order, rep.orders.last())
            && !(lo..=hi).contains(&o)
        {
            rep.failures.push(format!(
                "observed order {o:.3} on the finest pair, outside [{lo}, {hi}]"
            ));
        }
        let pts: Vec<(f64, f64)> = e.iter().map(|e| (e.h, e.l2)).collect();
        let slope: Vec<(f64, f64)> = match (pts.first(), pts.last()) {
            (Some(&a), Some(&b)) => vec![(a.0, b.1 * (a.0 / b.0).powi(2)), b],
            _ => Vec::new(),
        };
        std::fs::write(
            out.join("convergence.svg"),
            plot::chart(
                &format!("{name}: velocity error over mesh size"),
                "h (m)",
                "L2 error / U",
                &[
                    Series {
                        label: "Foehncraft",
                        points: &pts,
                        markers: true,
                    },
                    Series {
                        label: "slope 2",
                        points: &slope,
                        markers: false,
                    },
                ],
                true,
                true,
            ),
        )?;
    }
    if levels.len() >= 3 && rep.solves.len() == levels.len() {
        let d = f64::from(dimension(&case));
        let sizes: Vec<f64> = rep
            .solves
            .iter()
            .map(|s| (s.cells as f64).powf(-1.0 / d))
            .collect();
        let js: Vec<serde_json::Value> = rep
            .solves
            .iter()
            .map(|s| serde_json::to_value(s).unwrap_or_default())
            .collect();
        rep.uncertainty = crate::uncertainty::sweep(&js, &sizes);
        rep.uncertainty_of = "meshes";
    }
    for s in &rep.solves {
        if let (Some(max), Some(e)) = (v.max_error, &s.errors)
            && e.l2 > max
        {
            rep.failures.push(format!(
                "{}: relative L2 velocity error {:.3e} above {max:e}",
                s.mesh, e.l2
            ));
        }
        if let (Some(tol), Some(r)) = (v.reference_tolerance, &s.reference)
            && r.u_max_dev.max(r.v_max_dev) > tol
        {
            rep.failures.push(format!(
                "{}: deviation from {} of {:.4} above {tol}",
                s.mesh,
                r.reference,
                r.u_max_dev.max(r.v_max_dev)
            ));
        }
        if s.mass_imbalance.abs() > 1e-6 {
            rep.failures.push(format!(
                "{}: net mass flow out of the domain {:.2e} of the largest boundary flow",
                s.mesh, s.mass_imbalance
            ));
        }
    }
    let series: Vec<Series> = histories
        .iter()
        .map(|(l, h)| Series {
            label: l,
            points: h,
            markers: false,
        })
        .collect();
    std::fs::write(
        out.join("residuals.svg"),
        plot::chart(
            &format!("{name}: largest residual"),
            "outer iteration",
            "residual",
            &series,
            false,
            true,
        ),
    )?;
    if !case.run.ranks.is_empty() {
        lg.line("");
        let failures = crate::ranks::report(&case, out, &mut |l| lg.line(l))?;
        rep.failures.extend(failures);
    }
    rep.seconds = t0.elapsed().as_secs_f64();
    end_frames(&mut lg, &rep);
    lg.sync(&mut steer.record);
    steer.record.flush();
    std::fs::write(out.join("flow.log"), &lg.text)?;
    std::fs::write(out.join("flow.json"), serde_json::to_string_pretty(&rep)?)?;
    std::fs::write(out.join("flow.md"), rep.to_markdown())?;
    Ok(rep)
}

/// Applies the `[[schedule]]` entries (D128) that are due before iteration
/// or step `p.n`: logs every change (` [I]` line and an `event` row of the
/// run record) and hands the new settings to the flow. Returns the case as
/// it was when something was applied.
pub(crate) fn schedule_step<F: FnMut(&str)>(
    case: &mut FlowCase,
    fl: &mut Flow,
    lg: &mut Log<F>,
    steer: &mut Steer,
    runner: &mut crate::schedule::Runner,
    p: crate::schedule::Position,
    restart: bool,
) -> Result<Option<FlowCase>> {
    let what = if case.time.is_some() {
        "step"
    } else {
        "iteration"
    };
    let Some((new, changes)) = runner.fire(case, p)? else {
        return Ok(None);
    };
    for c in &changes {
        let m = format!(
            "{what} {}: schedule ({}): {}: {} -> {}{}",
            p.n,
            c.when,
            c.key,
            c.old,
            c.new,
            if restart { " (at the restart)" } else { "" }
        );
        lg.note(&m);
        if steer.record.on() {
            steer.record.row(&serde_json::json!({
                "type": "event",
                "code": diag::RUN_SCHEDULED,
                "severity": "note",
                "message": m,
                "at": p.n,
                "key": c.key,
                "old": c.old,
                "new": c.new,
            }));
        }
        steer.remember(diag::RUN_SCHEDULED, p.n, &m);
    }
    fl.set_settings(new.solver.clone())?;
    Ok(Some(std::mem::replace(case, new)))
}

/// Applies the settings asked for through the live control port (L19) after
/// iteration or step `n`, one by one: a setting that is invalid for the case
/// is reported (`RUN024`) and left out, the others are applied (`RUN023`).
/// Every outcome is an `event` row with the equivalent `[[schedule]]` entry
/// (D128) and an entry of the status of the port. Returns the case as it was
/// when something was applied.
pub(crate) fn apply_live_settings<F: FnMut(&str)>(
    case: &mut FlowCase,
    fl: &mut Flow,
    lg: &mut Log<F>,
    steer: &mut Steer,
    warnings: &mut Vec<diag::Diagnostic>,
    n: usize,
) -> Result<Option<FlowCase>> {
    let pending = steer.take_settings();
    if pending.is_empty() {
        return Ok(None);
    }
    let (what, at) = if case.time.is_some() {
        ("step", "at_step")
    } else {
        ("iteration", "at_iteration")
    };
    let mut first = None;
    for p in pending {
        // the equivalent entry of a case file: effective from the next
        // iteration or step
        let entry = format!(
            "[[schedule]] {at} = {}, set = {{ {} = {} }}",
            n + 1,
            p.key,
            p.value
        );
        match control::with_setting(case, &p.key, &p.value) {
            Ok((new, old, now)) => {
                fl.set_settings(new.solver.clone())?;
                let prev = std::mem::replace(case, new);
                first.get_or_insert(prev);
                let m = format!("{what} {n}: control port: {}: {old} -> {now}", p.key);
                lg.note(&m);
                steer.event(diag::RUN_CONTROL_SET, n, &m);
                steer.live_note(serde_json::json!({
                    "at": n, "key": p.key, "old": old, "new": now,
                    "status": "applied", "schedule": entry,
                }));
            }
            Err(e) => {
                let m = format!(
                    "{what} {n}: control port: {} = {} not applied: {e:#}",
                    p.key, p.value
                );
                lg.warn(&m);
                steer.event(diag::RUN_CONTROL_REJECTED, n, &m);
                warnings.push(diag::Diagnostic::new(diag::RUN_CONTROL_REJECTED, m));
                steer.live_note(serde_json::json!({
                    "at": n, "key": p.key, "value": p.value.to_string(),
                    "status": "rejected", "error": format!("{e:#}"),
                }));
            }
        }
    }
    Ok(first)
}

/// Re-reads the case file at iteration or step `n` (D40) and applies what
/// may change during a run; returns the case as it was when something was
/// applied.
#[allow(clippy::too_many_arguments)]
pub(crate) fn reread_case<F: FnMut(&str)>(
    case: &mut FlowCase,
    fl: &mut Flow,
    lg: &mut Log<F>,
    steer: &mut Steer,
    warnings: &mut Vec<diag::Diagnostic>,
    n: usize,
    out: &Path,
    fired: usize,
) -> Result<Option<FlowCase>> {
    let what = if case.time.is_some() {
        "step"
    } else {
        "iteration"
    };
    let Some(src) = case.source.clone() else {
        lg.note(&format!("{what} {n}: no case file to re-read"));
        return Ok(None);
    };
    match control::reread(&src, case, fired) {
        Ok(rr) => {
            if rr.applied.is_empty() && rr.not_applied.is_empty() {
                lg.note(&format!("{what} {n}: case file re-read, nothing changed"));
            }
            for l in &rr.applied {
                let m = format!("{what} {n}: changed {l}");
                lg.note(&m);
                steer.event(diag::RUN_REREAD, n, &m);
            }
            for l in &rr.not_applied {
                let m = format!("{what} {n}: not applied (cannot change during a run): {l}");
                lg.warn(&m);
                steer.event(diag::RUN_NOT_APPLIED, n, &m);
                warnings.push(diag::Diagnostic::new(diag::RUN_NOT_APPLIED, m));
            }
            let old = if rr.applied.is_empty() {
                None
            } else {
                fl.set_settings(rr.case.solver.clone())?;
                Some(std::mem::replace(case, rr.case))
            };
            std::fs::write(
                out.join(format!("case.read_toml_{n}.toml")),
                toml::to_string_pretty(&*case)?,
            )?;
            Ok(old)
        }
        Err(e) => {
            let m =
                format!("{what} {n}: the case file was not applied, it does not validate: {e:#}");
            lg.warn(&m);
            steer.event(diag::RUN_REREAD_INVALID, n, &m);
            warnings.push(diag::Diagnostic::new(diag::RUN_REREAD_INVALID, m));
            Ok(None)
        }
    }
}

/// The files of a steady run: the inputs, the output directory and what
/// will be written there with its estimated size.
pub(crate) fn files_frame(
    case: &FlowCase,
    mesh: &Mesh,
    desc: &str,
    out: &Path,
) -> crate::console::Frame {
    use crate::console::size;
    let nc = mesh.n_cells;
    let walls = case
        .boundary
        .iter()
        .filter(|(_, b)| matches!(b, BoundaryCondition::Wall { .. }))
        .count();
    let rows = case.solver.max_iterations / case.output.history_every;
    let mut f = crate::logblocks::files(case, desc, out);
    f.push_rule();
    f.push_kv(
        "residuals.csv",
        format!(
            "up to {rows} rows (every {} iterations), about {}",
            case.output.history_every,
            size(rows * 90)
        ),
    );
    if case.output.vtkhdf {
        f.push_kv(
            "fields.vtkhdf",
            format!(
                "velocity, pressure{} on {nc} cells, about {}",
                if case.analytic.is_some() {
                    ", error"
                } else {
                    ""
                },
                size(nc * (8 * 6) + mesh.face_nodes.len() * 8 + mesh.points.len() * 24)
            ),
        );
    }
    if case.output.wall_csv && walls > 0 {
        f.push_kv(
            "wall_<boundary>.csv",
            format!("{walls} wall(s): pressure and wall shear stress per face"),
        );
    }
    for l in &case.output.lines {
        f.push_kv(
            &format!("line_{}.csv", l.name),
            format!("{} points", l.points),
        );
    }
    f.push_kv("flow.md, flow.json", "report and its data");
    f.push_kv("result.json", "summary for programs (exit code, checks)");
    f.push_kv("balance.csv, balance.svg", "mass and energy ledger");
    f.push_kv(
        "boundaries.csv",
        "area, mass flow, pressure and force per boundary",
    );
    f.push_kv("residuals.svg", "residual history");
    f
}

pub(crate) fn error_norms(
    fl: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    a: &Analytic,
    dim: u8,
) -> ErrorNorms {
    let (mut e2, mut einf, mut vol, mut pinf) = (0.0, 0.0f64, 0.0, 0.0f64);
    let (mut pmin, mut pmax) = (f64::INFINITY, f64::NEG_INFINITY);
    // periodic pairs (D29): the pressure of the solution is the periodic
    // part; the driving force β t̂ adds the mean gradient −β t̂
    let drive: Vec<([f64; 3], f64)> = mesh.periodic.as_ref().map_or(Vec::new(), |per| {
        per.pairs
            .iter()
            .zip(&fl.periodic_force)
            .map(|(pair, &b)| {
                let l = norm(pair.translation);
                (pair.translation.map(|t| t / l), b)
            })
            .collect()
    });
    let p_full = |c: usize| -> f64 {
        let x = geo.cell_centroid[c];
        drive.iter().fold(fl.p[c], |p, (t, b)| {
            p - b * (t[0] * x[0] + t[1] * x[1] + t[2] * x[2])
        })
    };
    // transient solutions and periodic pairs: the pressure level is free
    // (reference cell), so the volume means are compared
    let (mut pm, mut pem) = (0.0, 0.0);
    if a.is_transient() || a.free_pressure_level() || mesh.periodic.is_some() {
        let mut v = 0.0;
        for c in 0..mesh.n_cells {
            let w = geo.cell_volume[c];
            pm += p_full(c) * w;
            pem += a.pressure(geo.cell_centroid[c]) * w;
            v += w;
        }
        pm /= v;
        pem /= v;
    }
    for c in 0..mesh.n_cells {
        let x = geo.cell_centroid[c];
        let ue = a.velocity(x);
        let e = norm(sub(fl.u.get(c), ue));
        e2 += e * e * geo.cell_volume[c];
        einf = einf.max(e);
        vol += geo.cell_volume[c];
        let pe = a.pressure(x);
        pmin = pmin.min(pe);
        pmax = pmax.max(pe);
        pinf = pinf.max(((p_full(c) - pm) - (pe - pem)).abs());
    }
    let n = mesh.n_cells as f64;
    let h = if dim == 2 {
        let (lo, hi) = mesh.bounding_box();
        (vol / (hi[2] - lo[2]) / n).sqrt()
    } else {
        (vol / n).cbrt()
    };
    let us = a.velocity_scale().max(f64::MIN_POSITIVE);
    let prange = pmax - pmin;
    ErrorNorms {
        cells: mesh.n_cells,
        h,
        l2: (e2 / vol).sqrt() / us,
        linf: einf / us,
        pressure_linf: if prange > 0.0 { pinf / prange } else { pinf },
    }
}

pub(crate) fn compare_reference(
    fl: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    loc: &Locator,
    r: Reference,
    case: &FlowCase,
) -> ReferenceComparison {
    let (ys, us, xs, vs) = r.centre_lines();
    let lid = case.verification.lid_velocity.unwrap_or(1.0);
    let (lo, hi) = mesh.bounding_box();
    let zc = 0.5 * (lo[2] + hi[2]);
    let span = [hi[0] - lo[0], hi[1] - lo[1]];
    let mut u = Vec::new();
    let mut v = Vec::new();
    let (mut du, mut dv) = (0.0f64, 0.0f64);
    for i in 0..ys.len() {
        // keep the probe points inside the domain
        let y = ys[i].clamp(1e-9, 1.0 - 1e-9);
        let x = [lo[0] + 0.5 * span[0], lo[1] + y * span[1], zc];
        let c = locate(loc, geo, x);
        let val = value_at(fl, geo, c, x).0[0] / lid;
        let val = if ys[i] == 1.0 {
            1.0
        } else if ys[i] == 0.0 {
            0.0
        } else {
            val
        };
        du = du.max((val - us[i]).abs());
        u.push((ys[i], val));
        let xx = xs[i].clamp(1e-9, 1.0 - 1e-9);
        let x = [lo[0] + xx * span[0], lo[1] + 0.5 * span[1], zc];
        let c = locate(loc, geo, x);
        let val = if xs[i] == 1.0 || xs[i] == 0.0 {
            0.0
        } else {
            value_at(fl, geo, c, x).0[1] / lid
        };
        dv = dv.max((val - vs[i]).abs());
        v.push((xs[i], val));
    }
    let n = 200;
    let (mut u_line, mut v_line) = (Vec::new(), Vec::new());
    for i in 0..=n {
        let t = i as f64 / n as f64;
        let tt = t.clamp(1e-9, 1.0 - 1e-9);
        let x = [lo[0] + 0.5 * span[0], lo[1] + tt * span[1], zc];
        let uu = if i == n {
            1.0
        } else if i == 0 {
            0.0
        } else {
            value_at(fl, geo, locate(loc, geo, x), x).0[0] / lid
        };
        u_line.push((t, uu));
        let x = [lo[0] + tt * span[0], lo[1] + 0.5 * span[1], zc];
        let vv = if i == 0 || i == n {
            0.0
        } else {
            value_at(fl, geo, locate(loc, geo, x), x).0[1] / lid
        };
        v_line.push((t, vv));
    }
    ReferenceComparison {
        reference: r.label().into(),
        u_max_dev: du,
        v_max_dev: dv,
        u,
        v,
        u_line,
        v_line,
    }
}

pub(crate) fn wall_face_values(
    fl: &Solution,
    geo: &Geometry,
    f: usize,
    ni: usize,
) -> (f64, [f64; 3], [f64; 3]) {
    // pressure at the face (zero normal gradient: tangential extrapolation),
    // viscous force of the fluid on the face, wall shear stress
    let fb = f - ni;
    let c = fl.disc.owner[f] as usize;
    let s = geo.face_area[f];
    let r = sub(geo.face_centroid[f], geo.cell_centroid[c]);
    let ss = dot(s, s);
    let rn = dot(r, s) / ss;
    let rt = [r[0] - rn * s[0], r[1] - rn * s[1], r[2] - rn * s[2]];
    let p = fl.p[c] + dot(fl.grad_p.get(c), rt);
    let ub = fl.ub.get(fb);
    let up = fl.u.get(c);
    let k = fl.disc.k.get(f);
    let dg = fl.disc.dgeo[f];
    let mu = fl.mu_eff_at(f);
    let mut fv = [0.0; 3];
    for i in 0..3 {
        fv[i] = -mu * (dg * (ub[i] - up[i]) + dot(k, fl.grad_u[i].get(c)));
    }
    let a = ss.sqrt();
    let n = [s[0] / a, s[1] / a, s[2] / a];
    // shear: tangential part of the viscous force per area, on the wall
    let fn_ = dot(fv, n);
    let tau = [
        (fv[0] - fn_ * n[0]) / a,
        (fv[1] - fn_ * n[1]) / a,
        (fv[2] - fn_ * n[2]) / a,
    ];
    (p, fv, tau)
}

/// The torques and machine values of the end state: the log block, the
/// warnings (`ROT001`) and the values for the report (D184).
pub(crate) fn rotor_report<F: FnMut(&str)>(
    sol: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    case: &FlowCase,
    lg: &mut Log<'_, F>,
    warnings: &mut Vec<crate::diag::Diagnostic>,
) -> Option<crate::rotor::Values> {
    let v = crate::rotor::evaluate(sol, mesh, geo, case);
    if v.rotors.is_empty() {
        return None;
    }
    lg.frame(&crate::logblocks::rotors(&v));
    for d in crate::rotor::check(&v) {
        lg.warn(&format!("{}: {}", d.code, d.message));
        warnings.push(d);
    }
    Some(v)
}

/// The transposed part of the viscous force of the fluid on a moving wall
/// face (D184): −μ (∇u)ᵀ·S, with the gradient of the cell next to the face.
/// The momentum rows hold the Laplacian form ∇·(μ∇u), whose force
/// [`wall_face_values`] gives; the full stress μ(∇u + ∇uᵀ) has this term
/// besides. It is zero on a wall at rest (no slip), and on a wall that moves
/// in its own plane without curvature; on a turning cylinder it is the
/// difference between the torque of the Laplacian form (which misses the
/// −u_θ/r of the shear stress) and the true one.
pub(crate) fn wall_face_transposed(fl: &Solution, geo: &Geometry, f: usize, ni: usize) -> [f64; 3] {
    let ub = fl.ub.get(f - ni);
    if ub == [0.0; 3] {
        return [0.0; 3];
    }
    let c = fl.disc.owner[f] as usize;
    let s = geo.face_area[f];
    let mu = fl.mu_eff_at(f);
    let mut out = [0.0; 3];
    for (j, sj) in s.iter().enumerate() {
        let g = fl.grad_u[j].get(c);
        for (o, gi) in out.iter_mut().zip(g) {
            *o -= mu * sj * gi;
        }
    }
    out
}

pub(crate) fn boundary_summaries(
    fl: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    case: &FlowCase,
) -> Vec<BoundarySummary> {
    let ni = mesh.n_interior_faces();
    // the torque of a wall about the axis of its rotor, or of the first
    // turning zone (D184)
    let rotors = crate::rotor::rotors(case, mesh);
    let offsets = case
        .pressure_frame
        .as_ref()
        .and_then(|f| f.cell_offsets(case, mesh));
    mesh.boundaries
        .iter()
        .enumerate()
        .map(|(b, p)| {
            let bc = &case.boundary[&p.name];
            let wall = matches!(bc, BoundaryCondition::Wall { .. });
            let axis = (wall && !rotors.is_empty()).then(|| {
                &rotors
                    .iter()
                    .find(|r| r.walls.contains(&b))
                    .unwrap_or(&rotors[0])
                    .turning
            });
            let (mut area, mut q, mut pa, mut torque) = (0.0, 0.0, 0.0, 0.0);
            let mut force = [0.0; 3];
            for f in p.range() {
                let s = geo.face_area[f];
                let a = norm(s);
                area += a;
                q += fl.mass_flux[f];
                let (pf, fv) = if wall {
                    let (pf, fv, _) = wall_face_values(fl, geo, f, ni);
                    (pf, fv)
                } else {
                    (fl.pb[f - ni], [0.0; 3])
                };
                pa += pf * a;
                let ff = [pf * s[0] + fv[0], pf * s[1] + fv[1], pf * s[2] + fv[2]];
                for i in 0..3 {
                    force[i] += ff[i];
                }
                if let Some(t) = axis {
                    // the torque has the full stress (D184)
                    let tr = wall_face_transposed(fl, geo, f, ni);
                    let ff = [ff[0] + tr[0], ff[1] + tr[1], ff[2] + tr[2]];
                    torque += dot(t.axis, cross(sub(geo.face_centroid[f], t.origin), ff));
                }
            }
            BoundarySummary {
                name: p.name.clone(),
                kind: bc.type_name().into(),
                area,
                mass_flow: q,
                // in the case file's convention (D239): the operating
                // pressure of each face's part added; the forces stay those
                // of p − p_op
                mean_pressure: {
                    let m = if area > 0.0 { pa / area } else { 0.0 };
                    match &offsets {
                        Some(o) if area > 0.0 => {
                            let w: f64 = p
                                .range()
                                .map(|f| o[mesh.owner[f] as usize] * norm(geo.face_area[f]))
                                .sum();
                            m + w / area
                        }
                        _ => m,
                    }
                },
                force,
                torque: axis.map(|_| torque),
            }
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn write_outputs(
    case: &FlowCase,
    out: &Path,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
    analytic: Option<&Analytic>,
    boundaries: &[BoundarySummary],
    loc: &Locator,
    at: Option<(usize, Option<f64>)>,
) -> Result<()> {
    let name = case_name(case);
    let mut csv = csv_preamble(
        "boundaries",
        &name,
        at,
        "boundary -, type -, area m², mass_flow kg/s (out of the domain), mean_pressure Pa, \
         force_x force_y force_z N, torque N m (of a wall about the axis of its rotor, D184)",
    );
    csv.push_str("boundary,type,area,mass_flow,mean_pressure,force_x,force_y,force_z,torque\n");
    for p in boundaries {
        let _ = writeln!(
            csv,
            "{},{},{:e},{:e},{:e},{:e},{:e},{:e},{}",
            p.name,
            p.kind,
            p.area,
            p.mass_flow,
            p.mean_pressure,
            p.force[0],
            p.force[1],
            p.force[2],
            p.torque.map_or(String::new(), |t| format!("{t:e}"))
        );
    }
    std::fs::write(out.join("boundaries.csv"), csv)?;
    let ni = mesh.n_interior_faces();
    if case.output.wall_csv {
        for p in &mesh.boundaries {
            if !matches!(case.boundary[&p.name], BoundaryCondition::Wall { .. }) {
                continue;
            }
            let turb = fl.state.turbulence.as_ref();
            let mut csv = csv_preamble(
                "wall",
                &name,
                at,
                if turb.is_some() {
                    "x y z m (face centroid), pressure Pa, tau_x tau_y tau_z Pa (wall shear \
                     stress), y_plus - (first cell), u_tau m/s (wall treatment), per scalar \
                     <name>_wall (wall value), <name>_flux (diffusive flux into the domain per \
                     unit area, scalar unit m/s), <name>_h (exchange velocity flux/(wall - cell \
                     value), m/s: h_w/(rho c_p) of a temperature) and <name>_t_plus (u_tau/h)"
                } else {
                    "x y z m (face centroid), pressure Pa, tau_x tau_y tau_z Pa (wall shear stress)"
                },
            );
            csv.push_str(if turb.is_some() {
                "x,y,z,pressure,tau_x,tau_y,tau_z,y_plus,u_tau"
            } else {
                "x,y,z,pressure,tau_x,tau_y,tau_z"
            });
            for st in fl.state.scalars.iter().filter(|s| s.wall.is_some()) {
                let n = &st.name;
                let _ = write!(csv, ",{n}_wall,{n}_flux,{n}_h");
                if turb.is_some() {
                    let _ = write!(csv, ",{n}_t_plus");
                }
            }
            csv.push('\n');
            let mut rows: Vec<(u64, String)> = p
                .range()
                .map(|f| {
                    let x = geo.face_centroid[f];
                    let (pf, _, tau) = wall_face_values(fl, geo, f, ni);
                    let mut r = format!(
                        "{:e},{:e},{:e},{:e},{:e},{:e},{:e}",
                        x[0], x[1], x[2], pf, tau[0], tau[1], tau[2]
                    );
                    if let Some(t) = turb {
                        let _ = write!(r, ",{:e},{:e}", t.y_plus[f - ni], t.u_tau[f - ni]);
                    }
                    // scalars: wall value, flux and exchange velocity
                    let area = norm(geo.face_area[f]);
                    for st in &fl.state.scalars {
                        let Some(w) = &st.wall else { continue };
                        let Some((pb, flux)) = scalar_wall_face(fl, geo, st, f, ni) else {
                            continue;
                        };
                        let fb = f - ni;
                        let ex = w.gamma[fb] / fl.phys.rho_face_at(f) * fl.disc.dgeo[f] / area;
                        let _ = write!(r, ",{pb:e},{flux:e},{ex:e}");
                        if let Some(t) = turb {
                            let _ = write!(r, ",{:e}", t.u_tau[fb] / ex);
                        }
                    }
                    (mesh.face_id[f], r)
                })
                .collect();
            // global face order: independent of the renumbering (R2)
            rows.sort_by_key(|r| r.0);
            for (_, r) in rows {
                csv.push_str(&r);
                csv.push('\n');
            }
            std::fs::write(out.join(format!("wall_{}.csv", p.name)), csv)?;
        }
    }
    for l in &case.output.lines {
        let mut units = String::from("s m (distance from `from`), x y z m, u v w m/s, p Pa");
        let mut head = String::from("s,x,y,z,u,v,w,p");
        if let Some(t) = &fl.state.turbulence {
            if t.epsilon {
                units
                    .push_str(", k m2/s2, epsilon m2/s3, nu_t m2/s, wall_distance m (cell values)");
                head.push_str(",k,epsilon,nu_t,wall_distance");
            } else {
                units.push_str(", k m2/s2, omega 1/s, nu_t m2/s, wall_distance m (cell values)");
                head.push_str(",k,omega,nu_t,wall_distance");
            }
            if t.zeta.is_some() {
                if case.turbulence.model.is_eb() {
                    units.push_str(", phi -, alpha - (cell values)");
                    head.push_str(",phi,alpha");
                } else {
                    units.push_str(", zeta -, f 1/s (cell values)");
                    head.push_str(",zeta,f");
                }
            }
        }
        for st in &fl.state.scalars {
            let _ = write!(units, ", {} (cell value)", st.name);
            let _ = write!(head, ",{}", st.name);
        }
        // rotating machinery (D183): the relative velocity and the
        // cylindrical components about the axis of the zone
        let frames = crate::derived::Frames::new(case);
        if frames.is_some() {
            units.push_str(
                ", u_rel v_rel w_rel m/s (relative to the turning zone), c_r c_u c_a m/s \
                 (radial, circumferential, axial), w_u m/s (relative circumferential)",
            );
            head.push_str(",u_rel,v_rel,w_rel,c_r,c_u,c_a,w_u");
        }
        let mut csv = csv_preamble("line", &name, at, &units);
        csv.push_str(&head);
        csv.push('\n');
        let d = sub(l.to, l.from);
        let len = norm(d);
        for i in 0..l.points {
            let t = i as f64 / (l.points - 1) as f64;
            let x = std::array::from_fn(|k| l.from[k] + t * d[k]);
            let c = locate(loc, geo, x);
            let (u, p) = value_at(fl, geo, c, x);
            let _ = write!(
                csv,
                "{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e}",
                t * len,
                x[0],
                x[1],
                x[2],
                u[0],
                u[1],
                u[2],
                p
            );
            if let Some(t) = &fl.state.turbulence {
                let _ = write!(
                    csv,
                    ",{:e},{:e},{:e},{:e}",
                    t.k[c], t.omega[c], t.nu_t[c], t.wall_distance[c]
                );
                if let Some(z) = &t.zeta {
                    let _ = write!(csv, ",{:e},{:e}", z.zeta[c], z.f[c]);
                }
            }
            for st in &fl.state.scalars {
                let _ = write!(csv, ",{:e}", st.phi[c]);
            }
            if let Some(fr) = &frames {
                let t = fr.at(mesh, c, x, u);
                let _ = write!(
                    csv,
                    ",{:e},{:e},{:e},{:e},{:e},{:e},{:e}",
                    t.w[0], t.w[1], t.w[2], t.c_r, t.c_u, t.c_a, t.w_u
                );
            }
            csv.push('\n');
        }
        std::fs::write(out.join(format!("line_{}.csv", l.name)), csv)?;
    }
    write_fields(
        case,
        &out.join("fields.vtkhdf"),
        mesh,
        geo,
        fl,
        analytic,
        at,
    )
}

/// The case name (the stem of the case file).
pub(crate) fn case_name(case: &FlowCase) -> String {
    case.source
        .as_ref()
        .and_then(|p| p.file_stem())
        .map_or("flow".into(), |s| s.to_string_lossy().into_owned())
}

/// The comment lines at the top of a CSV file (D75): format and version,
/// the Foehncraft version, the case, the iteration or step and time of
/// the state it holds, and the units of the columns (`-`: dimensionless).
pub(crate) fn csv_preamble(
    kind: &str,
    case: &str,
    at: Option<(usize, Option<f64>)>,
    units: &str,
) -> String {
    let when = match at {
        Some((n, Some(t))) => format!("; step {n}; time {t:e} s"),
        Some((n, None)) => format!("; iteration {n}"),
        None => String::new(),
    };
    format!(
        "# format foehncraft-{kind} 1; Foehncraft {}; case {case}{when}\n# units: {units}\n",
        foehncraft_core::VERSION
    )
}

/// Units of the fields of `fields*.vtkhdf` (attribute `units`, D75).
pub const FIELD_UNITS: &str = "coordinates m; velocity_x m/s; velocity_y m/s; velocity_z m/s; \
     speed m/s; pressure Pa; velocity_error m/s";

/// Writes the fields; `at` is the iteration or step and the time (D75).
#[cfg(feature = "hdf5")]
pub(crate) fn write_fields(
    case: &FlowCase,
    file: &Path,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
    analytic: Option<&Analytic>,
    at: Option<(usize, Option<f64>)>,
) -> Result<()> {
    use foehncraft_mesh::vtkhdf::{Attr, CellField, FaceField, write_all};
    if !case.output.vtkhdf {
        return Ok(());
    }
    let derived = crate::derived::compute(case, &case.output.derived, mesh, geo, fl)?;
    let speed: Vec<f64> = (0..mesh.n_cells).map(|c| norm(fl.u.get(c))).collect();
    // absolute pressure values (D239): the solver's gauge pressure plus the
    // operating pressure of each cell's fluid part
    let offsets = case
        .pressure_frame
        .as_ref()
        .and_then(|f| f.cell_offsets(case, mesh));
    let p_out: Option<Vec<f64>> = offsets
        .as_ref()
        .map(|o| (0..mesh.n_cells).map(|c| fl.p[c] + o[c]).collect());
    let mut fields = vec![
        CellField {
            name: "velocity_x",
            values: &fl.u.x,
        },
        CellField {
            name: "velocity_y",
            values: &fl.u.y,
        },
        CellField {
            name: "velocity_z",
            values: &fl.u.z,
        },
        CellField {
            name: "speed",
            values: &speed,
        },
        CellField {
            name: "pressure",
            values: p_out.as_deref().unwrap_or(&fl.p),
        },
    ];
    // gravity (stage 7b1): `pressure` is the reduced pressure p′, the
    // static pressure p = p′ + ρ_ref g·(x − x_ref) beside it (D153)
    let hydro = crate::buoyancy::hydrostatic(case);
    let static_p: Vec<f64>;
    if let Some(h) = hydro {
        static_p = (0..mesh.n_cells)
            .map(|c| {
                let p = fl.p[c] + h.head(geo.cell_centroid[c]);
                match &offsets {
                    Some(o) => p + o[c],
                    None => p,
                }
            })
            .collect();
        fields.push(CellField {
            name: "static_pressure",
            values: &static_p,
        });
    }
    // the variable density (stage 7c): ρ(T, p₀) of the cells
    // and the density of each fluid (D194), with that of the solid regions
    let rho_out: Option<Vec<f64>> = fl.phys.rho_var.as_ref().map(|r| {
        let mut v = r.0.clone();
        if fl.phys.zones && v.len() == mesh.n_cells {
            for (name, reg) in &case.regions {
                if reg.kind != crate::materials::RegionKind::Solid {
                    continue;
                }
                let Ok(m) = crate::energy::solid(case, &Env::default(), name) else {
                    continue;
                };
                for (c, rho) in v.iter_mut().enumerate() {
                    let z = &mesh.zones[mesh.cell_zone[c] as usize];
                    if reg.zones.iter().any(|n| n == z) {
                        *rho = m.rho;
                    }
                }
            }
        }
        v
    });
    if let Some(r) = &rho_out {
        fields.push(CellField {
            name: "density",
            values: r,
        });
    }
    let err: Vec<f64>;
    if let Some(a) = analytic {
        err = (0..mesh.n_cells)
            .map(|c| norm(sub(fl.u.get(c), a.velocity(geo.cell_centroid[c]))))
            .collect();
        fields.push(CellField {
            name: "velocity_error",
            values: &err,
        });
    }
    // passive scalars (D51) and, with an exact solution, their error
    let scalar_err: Vec<(String, Vec<f64>)> = crate::scalar::field_errors(case, mesh, geo, fl, at)?;
    // the turbulence model (stage 6): k, ω, ν_t, the wall distance and
    // the modified pressure p̂ of the solver; k-ε-ζ-f also ζ, f, the time
    // scale T and v² = ζ k
    let v2: Vec<f64>;
    if let Some(t) = &fl.state.turbulence {
        for (name, values) in [
            ("k", &t.k),
            (if t.epsilon { "epsilon" } else { "omega" }, &t.omega),
            ("nu_t", &t.nu_t),
            ("turbulence_wall_distance", &t.wall_distance),
            ("modified_pressure", &t.p_modified),
        ] {
            fields.push(CellField { name, values });
        }
        if let Some(z) = &t.zeta {
            v2 = z.zeta.iter().zip(&t.k).map(|(z, k)| z * k).collect();
            // the elliptic-blending k-ε: φ and α in the slots of ζ and f
            let eb = case.turbulence.model.is_eb();
            for (name, values) in [
                (if eb { "phi" } else { "zeta" }, &z.zeta),
                (if eb { "alpha" } else { "f" }, &z.f),
                ("turbulent_time_scale", &z.t),
                ("v2", &v2),
            ] {
                fields.push(CellField { name, values });
            }
        }
    }
    for st in &fl.state.scalars {
        fields.push(CellField {
            name: &st.name,
            values: &st.phi,
        });
    }
    for (name, values) in &scalar_err {
        fields.push(CellField { name, values });
    }
    for (name, values) in &derived.cells {
        fields.push(CellField { name, values });
    }
    let faces: Vec<FaceField> = derived
        .faces
        .iter()
        .map(|(name, values)| FaceField { name, values })
        .collect();
    let name = case_name(case);
    let units = format!(
        "{FIELD_UNITS}{}{}{}",
        if hydro.is_some() {
            "; static_pressure Pa"
        } else {
            ""
        },
        match &fl.state.turbulence {
            Some(t) if t.zeta.is_some() => {
                "; k m2/s2; epsilon m2/s3; nu_t m2/s; turbulence_wall_distance m; \
                 modified_pressure Pa; zeta 1; f 1/s; turbulent_time_scale s; v2 m2/s2"
            }
            Some(t) if t.epsilon => {
                "; k m2/s2; epsilon m2/s3; nu_t m2/s; turbulence_wall_distance m; \
                 modified_pressure Pa"
            }
            Some(_) => {
                "; k m2/s2; omega 1/s; nu_t m2/s; turbulence_wall_distance m; modified_pressure Pa"
            }
            None => "",
        },
        crate::derived::Fields::units(&case.output.derived)
    );
    let mut attrs = vec![
        ("foehncraft_format", Attr::Str("foehncraft-fields")),
        ("foehncraft_format_version", Attr::U64(1)),
        ("foehncraft_version", Attr::Str(foehncraft_core::VERSION)),
        ("case", Attr::Str(&name)),
        ("units", Attr::Str(&units)),
    ];
    // the convention of `pressure` and the operating pressure (D239 (4))
    if let Some(f) = case.pressure_frame.as_deref() {
        attrs.push(("pressure_values", Attr::Str(f.pressure_values.as_str())));
        attrs.push((
            "operating_pressure",
            Attr::F64(f.parts.first().map_or(0.0, |p| p.operating_pressure)),
        ));
    }
    match at {
        Some((n, Some(t))) => {
            attrs.push(("step", Attr::U64(n as u64)));
            attrs.push(("time", Attr::F64(t)));
        }
        Some((n, None)) => attrs.push(("iteration", Attr::U64(n as u64))),
        None => {}
    }
    write_all(file, mesh, &fields, &faces, &attrs)
}

#[cfg(not(feature = "hdf5"))]
pub(crate) fn write_fields(
    _: &FlowCase,
    _: &Path,
    _: &Mesh,
    _: &Geometry,
    _: &Solution,
    _: Option<&Analytic>,
    _: Option<(usize, Option<f64>)>,
) -> Result<()> {
    Ok(())
}

pub(crate) fn write_reference_svg(out: &Path, name: &str, r: &ReferenceComparison) -> Result<()> {
    let (ys, us, xs, vs) = match r.reference.contains("1000") {
        true => Reference::GhiaRe1000.centre_lines(),
        false => Reference::GhiaRe100.centre_lines(),
    };
    let ref_u: Vec<(f64, f64)> = ys.iter().zip(us.iter()).map(|(&y, &u)| (u, y)).collect();
    let ref_v: Vec<(f64, f64)> = xs.iter().zip(vs.iter()).map(|(&x, &v)| (x, v)).collect();
    let cu: Vec<(f64, f64)> = r.u_line.iter().map(|&(y, u)| (u, y)).collect();
    let cv = r.v_line.clone();
    std::fs::write(
        out.join("centreline_u.svg"),
        plot::chart(
            &format!("{name}: u on the vertical centre line"),
            "u / U",
            "y / L",
            &[
                Series {
                    label: "Foehncraft",
                    points: &cu,
                    markers: false,
                },
                Series {
                    label: "Ghia et al. 1982",
                    points: &ref_u,
                    markers: true,
                },
            ],
            false,
            false,
        ),
    )?;
    std::fs::write(
        out.join("centreline_v.svg"),
        plot::chart(
            &format!("{name}: v on the horizontal centre line"),
            "x / L",
            "v / U",
            &[
                Series {
                    label: "Foehncraft",
                    points: &cv,
                    markers: false,
                },
                Series {
                    label: "Ghia et al. 1982",
                    points: &ref_v,
                    markers: true,
                },
            ],
            false,
            false,
        ),
    )?;
    Ok(())
}

/// Output directory of a case: `[output] directory`, `--out`, or
/// `output/<case name>`.
/// Checks the case without solving it (`foehn flow --check`): builds the
/// mesh of the first level and prints the findings of the boundary-condition
/// check with their guidance; fails on an error.
pub fn check_case(case: &FlowCase, path: &Path, mut log: impl FnMut(&str)) -> Result<()> {
    let level = case.verification.levels.first().copied();
    let (mesh, geo, desc) = mesh_for(case, level)?;
    let (fluid_case, mesh, geo) = crate::cht::view(case, mesh, geo)?;
    let case = &fluid_case;
    log(&format!("{}: {desc}", path.display()));
    for n in crate::fluids::notes(&crate::fluids::check(case, &mesh)?) {
        log(&n);
    }
    // the property program on this mesh: regions name existing zones (D112)
    let phys_on = physics_on(case, &physics(case)?, &mesh, &geo)?;
    let mut findings = bc::check(
        &mesh,
        &geo,
        &case.boundary,
        &phys_on,
        bc::At {
            analytic: case.analytic().as_ref(),
            values: &case.compiled_boundaries()?,
            t: case.time.map_or(0.0, |t| t.start_time),
            turning: &case.turning(),
            hydrostatic: crate::buoyancy::hydrostatic(case),
        },
        case.check_context(),
    );
    findings.extend(crate::turbulence::setup_findings(case));
    findings.extend(crate::energy::setup_findings(case));
    findings.extend(crate::buoyancy::setup_findings(case, &mesh, &geo));
    for f in &findings {
        log(&f.to_string());
    }
    let errors = findings.iter().filter(|f| f.level == Level::Error).count();
    if errors > 0 {
        bail!("{errors} error(s) in the boundary conditions");
    }
    log("the case can be run");
    Ok(())
}

/// The case name and output directory.
pub fn output_dir(case: &FlowCase, path: &Path, out: Option<PathBuf>) -> (String, PathBuf) {
    let name = path
        .file_stem()
        .map_or("flow".into(), |s| s.to_string_lossy().into_owned());
    let dir = out
        .or_else(|| case.output.directory.clone())
        .unwrap_or_else(|| PathBuf::from("output").join(&name));
    (name, dir)
}

/// Fails when the report lists failures.
pub fn check(rep: &FlowReport) -> Result<()> {
    if rep.failures.is_empty() {
        Ok(())
    } else {
        bail!(
            "{} check(s) failed:\n  {}",
            rep.failures.len(),
            rep.failures.join("\n  ")
        )
    }
}

impl FlowReport {
    /// Markdown report.
    pub fn to_markdown(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "# Flow case `{}`\n", self.name);
        if let Some(re) = self.reynolds {
            let _ = writeln!(s, "Reynolds number {re}\n");
        }
        if !self.findings.is_empty() {
            let _ = writeln!(s, "Boundary-condition check:\n");
            for f in &self.findings {
                let _ = writeln!(s, "- {f}");
            }
            s.push('\n');
        }
        let _ = writeln!(
            s,
            "| Mesh | Iterations | Converged | Largest residual | ns/cell/iteration | L2 error | Linf error | p error | Reference deviation |\n|---|---|---|---|---|---|---|---|---|"
        );
        for v in &self.solves {
            let _ = writeln!(
                s,
                "| {} | {} | {} | {:.2e} | {:.0} | {} | {} | {} | {} |",
                v.mesh,
                v.iterations,
                if v.converged { "yes" } else { "no" },
                v.residuals.max(),
                v.ns_per_cell_iteration,
                v.errors
                    .as_ref()
                    .map_or("–".into(), |e| format!("{:.3e}", e.l2)),
                v.errors
                    .as_ref()
                    .map_or("–".into(), |e| format!("{:.3e}", e.linf)),
                v.errors
                    .as_ref()
                    .map_or("–".into(), |e| format!("{:.3e}", e.pressure_linf)),
                v.reference.as_ref().map_or("–".into(), |r| format!(
                    "u {:.4}, v {:.4}",
                    r.u_max_dev, r.v_max_dev
                )),
            );
        }
        for v in &self.solves {
            if let Some(r) = &v.remediation {
                let _ = writeln!(
                    s,
                    "\n**Cell quality (D50)** on {}: {} faces with skewness above {} (largest {:.3}) in {} cells; an optimised mesh avoids them (Gmsh: `Mesh.OptimizeNetgen = 1`). {}",
                    v.mesh,
                    r.skewed_faces,
                    r.skewness_limit,
                    r.max_skewness,
                    r.poor_cells,
                    if r.changes.is_empty() {
                        "The schemes were kept as set.".to_string()
                    } else {
                        format!(
                            "Remediation (`solver.remediation = \"auto\"`, `\"off\"` keeps the schemes): {}.",
                            r.changes
                                .iter()
                                .map(|c| format!("`{c}`"))
                                .collect::<Vec<_>>()
                                .join(", ")
                        )
                    }
                );
            }
        }
        for v in &self.solves {
            if let Some(r) = &v.operators {
                let _ = writeln!(
                    s,
                    "\n**Operator check (D236)** on {}: {}",
                    v.mesh,
                    r.findings.join(" ")
                );
            }
        }
        if !self.orders.is_empty() {
            let _ = writeln!(s, "\nObserved orders: {:.3?}", self.orders);
        }
        if !self.time_orders.is_empty() {
            let _ = writeln!(s, "\nObserved time orders: {:.3?}", self.time_orders);
        }
        s.push_str(&crate::uncertainty::markdown(
            &self.uncertainty,
            self.uncertainty_of,
        ));
        if self.solves.iter().any(|v| v.transient.is_some()) {
            let _ = writeln!(
                s,
                "\n| Mesh | Scheme | Steps | End time (s) | Δt (s) | Largest Courant | Outer per step | CG per step | ns/cell/step |\n|---|---|---|---|---|---|---|---|---|"
            );
            for v in &self.solves {
                if let Some(t) = &v.transient {
                    let _ = writeln!(
                        s,
                        "| {} | {} | {} | {:.6e} | {:.3e}–{:.3e} | {:.3} | {:.2} | {:.1} | {:.0} |",
                        v.cells,
                        t.scheme,
                        t.steps,
                        t.time,
                        t.dt_min,
                        t.dt_max,
                        t.cfl_max,
                        t.outer_mean,
                        t.cg_mean,
                        t.ns_per_cell_step
                    );
                }
            }
        }
        if self
            .solves
            .iter()
            .any(|v| v.benchmark.is_some() || v.coefficients.is_some())
        {
            let _ = writeln!(s);
            for v in &self.solves {
                if let Some([cd, cl]) = v.coefficients {
                    let _ = writeln!(
                        s,
                        "{} cells: force coefficients C_D {cd:.5}, C_L {cl:.5}\n",
                        v.cells
                    );
                }
                if let Some(b) = &v.benchmark {
                    let _ = writeln!(
                        s,
                        "{} cells, {}{}:\n\n| Quantity | Computed | Reference | Distance |\n|---|---|---|---|",
                        v.cells,
                        b.benchmark,
                        if b.periods > 0 {
                            format!(" ({} periods)", b.periods)
                        } else {
                            String::new()
                        }
                    );
                    for x in &b.values {
                        let _ = writeln!(
                            s,
                            "| {} | {:.5} | {} – {} | {:.4} |",
                            x.range.name, x.value, x.range.lower, x.range.upper, x.deviation
                        );
                    }
                    s.push('\n');
                }
            }
        }
        for v in &self.solves {
            if let Some(b) = &v.balance {
                s.push_str(&balance::markdown(b, v.cells));
            }
        }
        if self.stopped {
            let _ = writeln!(s, "\nThe run was stopped by `.stop`.");
        }
        if let Some(v) = self.solves.last() {
            let _ = writeln!(
                s,
                "\n| Boundary | Type | Mass flow (kg/s) | Mean pressure (Pa) | Force (N) | Torque (N m) |\n|---|---|---|---|---|---|"
            );
            for p in &v.boundaries {
                let _ = writeln!(
                    s,
                    "| {} | {} | {:.6e} | {:.6e} | [{:.4e}, {:.4e}, {:.4e}] | {} |",
                    p.name,
                    p.kind,
                    p.mass_flow,
                    p.mean_pressure,
                    p.force[0],
                    p.force[1],
                    p.force[2],
                    p.torque.map_or("-".to_string(), |t| format!("{t:.4e}"))
                );
            }
            if let Some(r) = &v.rotors {
                s.push_str(&crate::rotor::markdown(r));
            }
            let t = &v.timings;
            let mut sweeps = String::new();
            if v.sweep_setups > 0 {
                let _ = write!(sweeps, ", {} of the momentum sweeps", v.sweep_setups);
            }
            if v.turbulence_setups > 0 {
                let _ = write!(
                    sweeps,
                    ", {} of the turbulence equations",
                    v.turbulence_setups
                );
            }
            let mut models = String::new();
            if t.turbulence > 0.0 {
                let _ = write!(models, ", turbulence model {:.2}", t.turbulence);
            }
            if t.scalars > 0.0 {
                let _ = write!(models, ", scalars {:.2}", t.scalars);
            }
            let _ = writeln!(
                s,
                "\nSeconds per phase (R6): gradients {:.2}, momentum assembly {:.2}, momentum solve {:.2}, pressure assembly {:.2}, pressure solve {:.2}, multigrid setup {:.2} ({} setups{sweeps}){models}",
                t.gradients,
                t.momentum_assembly,
                t.momentum_solve,
                t.pressure_assembly,
                t.pressure_solve,
                t.amg_setup,
                v.amg_setups
            );
        }
        if self.failures.is_empty() {
            let _ = writeln!(s, "\nAll checks passed ({:.1} s).", self.seconds);
        } else {
            let _ = writeln!(s, "\nFailed checks:\n");
            for f in &self.failures {
                let _ = writeln!(s, "- {f}");
            }
        }
        s
    }
}
