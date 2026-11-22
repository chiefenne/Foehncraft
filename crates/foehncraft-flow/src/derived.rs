//! Derived fields (D72): quantities computed from the solution when the
//! fields are written (`[output] derived = [...]`), or after the run from a
//! restart file (`foehn post derive`), with the same code, so both give the
//! same numbers.
//!
//! Cell fields, from the cell velocity gradient G with Gᵢⱼ = ∂uᵢ/∂xⱼ, its
//! symmetric part S = ½ (G + Gᵀ) (rate of strain) and its antisymmetric
//! part Ω = ½ (G − Gᵀ) (rate of rotation):
//!
//! | name | fields | definition | unit |
//! |---|---|---|---|
//! | `vorticity` | `vorticity_x`, `_y`, `_z`, `vorticity_magnitude` | ω = ∇ × u | 1/s |
//! | `q-criterion` | `q_criterion` | Q = ½ (‖Ω‖² − ‖S‖²); vortices where Q > 0 (Hunt, Wray and Moin 1988) | 1/s² |
//! | `lambda2` | `lambda2` | the middle eigenvalue of S² + Ω²; vortices where λ₂ < 0 (Jeong and Hussain 1995) | 1/s² |
//! | `enstrophy` | `enstrophy` | ½ ω·ω | 1/s² |
//! | `strain-rate` | `strain_rate` | √(2 S:S), the shear rate | 1/s |
//! | `pressure-coefficient` | `pressure_coefficient` | (p − p_ref) / (½ ρ_ref U_ref²) | - |
//! | `total-pressure` | `total_pressure` | p + ½ ρ u·u | Pa |
//! | `relative-velocity` | `relative_velocity_x`, `_y`, `_z`, `relative_velocity_magnitude` | **w** = **c** − **u**, the velocity relative to the turning zone of the cell (D183; **c** elsewhere) | m/s |
//! | `frame-velocity` | `frame_velocity_x`, `_y`, `_z`, `frame_velocity_magnitude` | **u** = **Ω** × (**x** − **x**₀) of the turning zone of the cell (0 elsewhere) | m/s |
//! | `turbomachinery` | `velocity_radial`, `velocity_circumferential`, `velocity_axial`, `velocity_meridional`, `relative_velocity_circumferential`, `flow_angle`, `relative_flow_angle` | c_r, c_u, c_a, c_m = (c_r² + c_a²)^½, w_u about the axis of the zone (along **Ω**), α = atan2(c_m, c_u), β = atan2(c_m, −w_u) | m/s, ° |
//! | `relative-total-pressure` | `relative_total_pressure` | p + ½ ρ w·w | Pa |
//! | `rothalpy` | `rothalpy` | I = p/ρ + ½ w·w − ½ u·u | m²/s² |
//! | `wall-distance` | `wall_distance` | the exact distance of the cell centroid to the nearest no-slip wall face (D48, a BVH closest-point search, D108) | m |
//! | `turbulence-features` | `feature_lambda1` … `feature_lambda5`, `feature_wall_reynolds`, `feature_time_ratio`, `feature_viscosity_ratio`, `feature_shielding`, `feature_q`, `feature_pressure_gradient`, `feature_curvature` | the invariant features of D166 (below); needs a turbulence model | - |
//! | `turbulence-scales` | `dissipation_rate`, `specific_dissipation_rate`, `time_scale`, `length_scale`, `viscosity_ratio` | ε (β* k ω for SST), ω (ε/(C_μ k) for the ε models), the time scale T of the model (1/(β* ω) for SST, k/ε for the k-ε family, the bounded T of ζ-f and EB), k^{3/2}/ε, ν_t/ν; needs a turbulence model | m²/s³, 1/s, s, m, - |
//! | `reynolds-stress` | `reynolds_stress_xx`, `_yy`, `_zz`, `_xy`, `_yz`, `_xz` | the Boussinesq estimate ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ of u'ᵢu'ⱼ (spec 2.3); needs a turbulence model | m²/s² |
//! | `wall-normal` | `wall_normal_x`, `_y`, `_z` | n = ∇α/\|∇α\| of the elliptic-blending k-ε (0 where ∇α vanishes), for the GGDH of stage 7; α = 0 on walls, zero gradient elsewhere | - |
//! | `curvature-correction` | `rotation_function` | f_r1 of the curvature correction of SST (D170; 1 in pure shear, above 1 where curvature or rotation destabilise, below where they stabilise); needs `turbulence.curvature_correction` | - |
//!
//! Boundary face fields, on the boundary blocks of the file (zero on boundaries
//! that are not walls), from the same wall values as `wall_<boundary>.csv`:
//!
//! | name | fields | definition | unit |
//! |---|---|---|---|
//! | `wall-shear-stress` | `wall_shear_stress_x`, `_y`, `_z`, `wall_shear_stress_magnitude` | the tangential viscous stress of the fluid on the wall | Pa |
//! | `y-plus` | `y_plus` | y⁺ = d u_τ / ν with u_τ = √(τ_w / ρ) and d the wall distance of the cell centroid | - |
//! | `turbulence-wall` | `friction_velocity`, `y_star`, `blending` | u_τ of the wall treatment, y* = d C_μ^¼ k_P^½ / ν of the wall-adjacent cell (β* for SST), and the viscous weight of the wall treatment (spec 2.3): 1 where the wall is resolved (`low-re`, the sublayer of `standard`), 1 − tanh((y⁺/5)²) of `automatic`, e^{−Γ} of `compound`, Γ = dU⁺/dy⁺ of `compound-2023`, 0 in the log layer and for `scalable`; needs a turbulence model | m/s, -, - |
//!
//! The fields of rotating machinery (D183) need an `[[mrf]]` or a
//! `[[sliding]]` table and use the notation of turbomachinery: **c** is
//! the absolute (solved) velocity, **u** = **Ω** × (**x** − **x**₀) the
//! velocity of the frame (the blade speed), **w** = **c** − **u** the
//! relative velocity. A cell in a turning zone is seen in the frame of its
//! zone, any other cell in the fixed frame (**w** = **c**, **u** = 0), so a
//! pump shows the relative flow in the impeller and the absolute flow in
//! the volute; `[output] frame = "<zone>"` takes the frame of that zone in
//! every cell. On a sliding mesh **x** is the current position of the cell.
//! The cylindrical components are taken about the axis of the cell's
//! frame (outside the zones about that of `frame`, else of the first
//! table), which points along **Ω**: the blade speed and a swirl in the
//! direction of rotation are positive. The flow angles are in degrees and
//! measured from the circumferential direction: α between **c** and **u**,
//! β between **w** and −**u** (90° for a flow without swirl relative to
//! the blades; the convention of pump design). Rothalpy is constant along
//! a relative streamline of a steady flow without losses, so its drop in a
//! rotor shows where energy is lost.
//!
//! The features of `turbulence-features` (D166; the definitions of the
//! export of D106), with the turbulence time scale τ = 1/(β* ω) = k/ε,
//! the turbulence length L = √k/(β* ω), the wall distance d of the model
//! and |S| = √(2 S:S):
//!
//! - λ₁ … λ₅: tr Ŝ², tr Ŵ², tr Ŝ³, tr Ŵ²Ŝ, tr Ŵ²Ŝ² of Ŝ = τS, Ŵ = τΩ
//!   (Pope 1975);
//! - wall Reynolds number d√k/ν (based on k, not on the cell);
//! - time-scale ratio τ |S| (turbulent over mean time scale);
//! - viscosity ratio ν_t/ν;
//! - shielding f_d = 1 − tanh((8 r_d)³), r_d = (ν_t + ν)/(‖G‖ κ² d²)
//!   (delayed DES, Spalart et al. 2006): 0 in attached boundary layers;
//! - Q = (‖Ω‖² − ‖S‖²)/(‖Ω‖² + ‖S‖²), in [−1, 1];
//! - pressure gradient along the flow u·∇p/(|u| |∇p|), in [−1, 1]
//!   (velocity in the frame of the walls);
//! - streamline curvature |u × (G u)| L/|u|³.
//!
//! Every quotient has its denominator floored at 10⁻³⁰⁰ (zero for a fluid at
//! rest). ‖·‖ is the Frobenius norm. The eigenvalues of λ₂ come from Jacobi
//! rotations, exact to round-off also where two of them are equal (as in
//! every 2D vortex). Turbulence quantities follow with the turbulence
//! models (`turbulence-scales`, `reynolds-stress`, `wall-normal`,
//! `turbulence-wall`, stage 6 D72); the wall values of the passive scalars
//! (h_w and T⁺) are the columns `<name>_h` and `<name>_t_plus` of
//! `wall_<boundary>.csv`.

use crate::bc::BoundaryCondition;
use crate::case::FlowCase;
use crate::simplec::Solution;
use anyhow::Result;
use foehncraft_mesh::geometry::{dot, norm, sub};
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// A derived field of `[output] derived`.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Derived {
    /// ω = ∇ × u, components and magnitude (1/s).
    Vorticity,
    /// Q = ½ (‖Ω‖² − ‖S‖²) (1/s²).
    QCriterion,
    /// Middle eigenvalue of S² + Ω² (1/s²).
    Lambda2,
    /// ½ ω·ω (1/s²).
    Enstrophy,
    /// √(2 S:S) (1/s).
    StrainRate,
    /// (p − p_ref) / (½ ρ_ref U_ref²) (-).
    PressureCoefficient,
    /// p + ½ ρ u·u (Pa).
    TotalPressure,
    /// **w** = **c** − **Ω** × (**x** − **x**₀), the velocity relative to
    /// the turning zone of the cell, components and magnitude (m/s, D183).
    RelativeVelocity,
    /// **u** = **Ω** × (**x** − **x**₀) of the turning zone of the cell,
    /// components and magnitude (m/s, D183).
    FrameVelocity,
    /// Cylindrical components about the axis of the zone and the flow
    /// angles α and β (m/s and degrees, D183).
    Turbomachinery,
    /// p + ½ ρ w·w (Pa, D183).
    RelativeTotalPressure,
    /// p/ρ + ½ w·w − ½ u·u (m²/s², D183).
    Rothalpy,
    /// Exact distance of the cell centroid to the nearest wall (m).
    WallDistance,
    /// Wall shear stress on the walls, components and magnitude (Pa).
    WallShearStress,
    /// y⁺ of the wall-adjacent cells, on the walls (-).
    YPlus,
    /// The invariant features of D166 (-).
    TurbulenceFeatures,
    /// ε, ω, the time and length scales and ν_t/ν of the turbulence
    /// model (stage 6, D72).
    TurbulenceScales,
    /// The Boussinesq Reynolds stress ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ, six components
    /// (m²/s², D72).
    ReynoldsStress,
    /// n = ∇α/|∇α| of the elliptic-blending k-ε (-, D72).
    WallNormal,
    /// u_τ and y* of the wall treatment, on the walls (D72).
    TurbulenceWall,
    /// f_r1 of the curvature correction of SST (-, D170).
    CurvatureCorrection,
}

impl Derived {
    /// Every derived field, in the order they are written.
    pub const ALL: [Derived; 21] = [
        Derived::Vorticity,
        Derived::QCriterion,
        Derived::Lambda2,
        Derived::Enstrophy,
        Derived::StrainRate,
        Derived::PressureCoefficient,
        Derived::TotalPressure,
        Derived::RelativeVelocity,
        Derived::FrameVelocity,
        Derived::Turbomachinery,
        Derived::RelativeTotalPressure,
        Derived::Rothalpy,
        Derived::WallDistance,
        Derived::WallShearStress,
        Derived::YPlus,
        Derived::TurbulenceFeatures,
        Derived::TurbulenceScales,
        Derived::ReynoldsStress,
        Derived::WallNormal,
        Derived::TurbulenceWall,
        Derived::CurvatureCorrection,
    ];

    /// The name in the case file.
    pub fn name(self) -> &'static str {
        match self {
            Derived::Vorticity => "vorticity",
            Derived::QCriterion => "q-criterion",
            Derived::Lambda2 => "lambda2",
            Derived::Enstrophy => "enstrophy",
            Derived::StrainRate => "strain-rate",
            Derived::PressureCoefficient => "pressure-coefficient",
            Derived::TotalPressure => "total-pressure",
            Derived::RelativeVelocity => "relative-velocity",
            Derived::FrameVelocity => "frame-velocity",
            Derived::Turbomachinery => "turbomachinery",
            Derived::RelativeTotalPressure => "relative-total-pressure",
            Derived::Rothalpy => "rothalpy",
            Derived::WallDistance => "wall-distance",
            Derived::WallShearStress => "wall-shear-stress",
            Derived::YPlus => "y-plus",
            Derived::TurbulenceFeatures => "turbulence-features",
            Derived::TurbulenceScales => "turbulence-scales",
            Derived::ReynoldsStress => "reynolds-stress",
            Derived::WallNormal => "wall-normal",
            Derived::TurbulenceWall => "turbulence-wall",
            Derived::CurvatureCorrection => "curvature-correction",
        }
    }

    /// The field from its name in the case file.
    pub fn parse(s: &str) -> Option<Derived> {
        Self::ALL.into_iter().find(|d| d.name() == s)
    }

    /// The datasets it writes, with their units.
    pub fn datasets(self) -> &'static [(&'static str, &'static str)] {
        match self {
            Derived::Vorticity => &[
                ("vorticity_x", "1/s"),
                ("vorticity_y", "1/s"),
                ("vorticity_z", "1/s"),
                ("vorticity_magnitude", "1/s"),
            ],
            Derived::QCriterion => &[("q_criterion", "1/s²")],
            Derived::Lambda2 => &[("lambda2", "1/s²")],
            Derived::Enstrophy => &[("enstrophy", "1/s²")],
            Derived::StrainRate => &[("strain_rate", "1/s")],
            Derived::PressureCoefficient => &[("pressure_coefficient", "-")],
            Derived::TotalPressure => &[("total_pressure", "Pa")],
            Derived::RelativeVelocity => &[
                ("relative_velocity_x", "m/s"),
                ("relative_velocity_y", "m/s"),
                ("relative_velocity_z", "m/s"),
                ("relative_velocity_magnitude", "m/s"),
            ],
            Derived::FrameVelocity => &[
                ("frame_velocity_x", "m/s"),
                ("frame_velocity_y", "m/s"),
                ("frame_velocity_z", "m/s"),
                ("frame_velocity_magnitude", "m/s"),
            ],
            Derived::Turbomachinery => &[
                ("velocity_radial", "m/s"),
                ("velocity_circumferential", "m/s"),
                ("velocity_axial", "m/s"),
                ("velocity_meridional", "m/s"),
                ("relative_velocity_circumferential", "m/s"),
                ("flow_angle", "°"),
                ("relative_flow_angle", "°"),
            ],
            Derived::RelativeTotalPressure => &[("relative_total_pressure", "Pa")],
            Derived::Rothalpy => &[("rothalpy", "m²/s²")],
            Derived::WallDistance => &[("wall_distance", "m")],
            Derived::WallShearStress => &[
                ("wall_shear_stress_x", "Pa"),
                ("wall_shear_stress_y", "Pa"),
                ("wall_shear_stress_z", "Pa"),
                ("wall_shear_stress_magnitude", "Pa"),
            ],
            Derived::YPlus => &[("y_plus", "-")],
            Derived::TurbulenceFeatures => &[
                ("feature_lambda1", "-"),
                ("feature_lambda2", "-"),
                ("feature_lambda3", "-"),
                ("feature_lambda4", "-"),
                ("feature_lambda5", "-"),
                ("feature_wall_reynolds", "-"),
                ("feature_time_ratio", "-"),
                ("feature_viscosity_ratio", "-"),
                ("feature_shielding", "-"),
                ("feature_q", "-"),
                ("feature_pressure_gradient", "-"),
                ("feature_curvature", "-"),
            ],
            Derived::TurbulenceScales => &[
                ("dissipation_rate", "m²/s³"),
                ("specific_dissipation_rate", "1/s"),
                ("time_scale", "s"),
                ("length_scale", "m"),
                ("viscosity_ratio", "-"),
            ],
            Derived::ReynoldsStress => &[
                ("reynolds_stress_xx", "m²/s²"),
                ("reynolds_stress_yy", "m²/s²"),
                ("reynolds_stress_zz", "m²/s²"),
                ("reynolds_stress_xy", "m²/s²"),
                ("reynolds_stress_yz", "m²/s²"),
                ("reynolds_stress_xz", "m²/s²"),
            ],
            Derived::WallNormal => &[
                ("wall_normal_x", "-"),
                ("wall_normal_y", "-"),
                ("wall_normal_z", "-"),
            ],
            Derived::TurbulenceWall => &[
                ("friction_velocity", "m/s"),
                ("y_star", "-"),
                ("blending", "-"),
            ],
            Derived::CurvatureCorrection => &[("rotation_function", "-")],
        }
    }

    /// True for the fields on the boundary faces (the boundary blocks).
    pub fn on_faces(self) -> bool {
        matches!(
            self,
            Derived::WallShearStress | Derived::YPlus | Derived::TurbulenceWall
        )
    }

    /// True for the fields that need a turbulence model.
    pub fn of_turbulence(self) -> bool {
        matches!(
            self,
            Derived::TurbulenceFeatures
                | Derived::TurbulenceScales
                | Derived::ReynoldsStress
                | Derived::WallNormal
                | Derived::TurbulenceWall
                | Derived::CurvatureCorrection
        )
    }

    /// True for the fields of rotating machinery (D183), which need a
    /// turning zone.
    pub fn of_rotors(self) -> bool {
        matches!(
            self,
            Derived::RelativeVelocity
                | Derived::FrameVelocity
                | Derived::Turbomachinery
                | Derived::RelativeTotalPressure
                | Derived::Rothalpy
        )
    }
}

/// The velocities of rotating machinery at a point (D183): absolute **c**
/// (the solved velocity), frame **u**, relative **w** = **c** − **u**, the
/// cylindrical components and the flow angles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Turbo {
    /// **w** = **c** − **u** (m/s).
    pub w: [f64; 3],
    /// **u** = **Ω** × (**x** − **x**₀) of the frame, 0 in the fixed one (m/s).
    pub u: [f64; 3],
    /// c_r, radial, away from the axis (m/s).
    pub c_r: f64,
    /// c_u, circumferential, positive in the direction of rotation (m/s).
    pub c_u: f64,
    /// c_a, axial, along **Ω** (m/s).
    pub c_a: f64,
    /// c_m = (c_r² + c_a²)^½, meridional (m/s); w_m = c_m.
    pub c_m: f64,
    /// w_u = c_u − u_u, circumferential relative (m/s).
    pub w_u: f64,
    /// α = atan2(c_m, c_u), between **c** and the circumferential direction (°).
    pub alpha: f64,
    /// β = atan2(c_m, −w_u), between **w** and the negative circumferential
    /// direction (°).
    pub beta: f64,
}

/// The velocities of rotating machinery of the velocity `c` at `x`, seen
/// in `frame` (the fixed frame for `None`), with the cylindrical
/// components about the axis of `axis`. On the axis itself the radial and
/// circumferential directions are taken as zero.
pub fn turbo(
    c: [f64; 3],
    x: [f64; 3],
    frame: Option<&crate::bc::Turning>,
    axis: &crate::bc::Turning,
) -> Turbo {
    let u = frame.map_or([0.0; 3], |t| t.at(x));
    let w = sub(c, u);
    let a = axis.axis;
    let d = sub(x, axis.origin);
    let da = dot(d, a);
    let rv = [d[0] - da * a[0], d[1] - da * a[1], d[2] - da * a[2]];
    let r = norm(rv);
    let er = if r > 0.0 { rv.map(|v| v / r) } else { [0.0; 3] };
    let et = [
        a[1] * er[2] - a[2] * er[1],
        a[2] * er[0] - a[0] * er[2],
        a[0] * er[1] - a[1] * er[0],
    ];
    let (c_r, c_u, c_a) = (dot(c, er), dot(c, et), dot(c, a));
    let c_m = c_r.hypot(c_a);
    let w_u = dot(w, et);
    Turbo {
        w,
        u,
        c_r,
        c_u,
        c_a,
        c_m,
        w_u,
        alpha: c_m.atan2(c_u).to_degrees(),
        beta: c_m.atan2(-w_u).to_degrees(),
    }
}

/// The frames of the fields of rotating machinery (D183): the turning
/// zones of `[[mrf]]` and `[[sliding]]`, and the zone of `[output] frame`.
#[derive(Debug, Clone)]
pub struct Frames {
    zones: Vec<crate::bc::Turning>,
    fixed: Option<usize>,
}

impl Frames {
    /// The frames of `case`; `None` without a turning zone.
    pub fn new(case: &FlowCase) -> Option<Frames> {
        let zones = case.turning();
        if zones.is_empty() {
            return None;
        }
        let fixed = case.output.frame.as_ref().and_then(|f| {
            zones
                .iter()
                .position(|t| t.zone.as_deref() == Some(f.as_str()))
        });
        Some(Frames { zones, fixed })
    }

    /// The frame that local cell `c` is seen in: that of `[output] frame`,
    /// else that of its turning zone, else the fixed one (`None`).
    pub fn frame(&self, mesh: &Mesh, c: usize) -> Option<&crate::bc::Turning> {
        match self.fixed {
            Some(i) => Some(&self.zones[i]),
            None => crate::bc::Turning::of(&self.zones, mesh, c),
        }
    }

    /// The values of the velocity `v` at `x` in local cell `c`.
    pub fn at(&self, mesh: &Mesh, c: usize, x: [f64; 3], v: [f64; 3]) -> Turbo {
        let frame = self.frame(mesh, c);
        turbo(v, x, frame, frame.unwrap_or(&self.zones[0]))
    }
}

/// `[output.reference]`: the reference state of the pressure coefficient.
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct ReferenceValues {
    /// Reference velocity U_ref (m/s); default `[output.coefficients]
    /// velocity`.
    pub velocity: Option<f64>,
    /// Reference pressure p_ref (Pa).
    pub pressure: f64,
    /// Reference density ρ_ref (kg/m³); default `[fluid] density`.
    pub density: Option<f64>,
}

/// The reference dynamic pressure ½ ρ_ref U_ref² of `case` and p_ref.
fn reference(case: &FlowCase) -> Option<(f64, f64)> {
    let r = &case.output.reference;
    let u = r
        .velocity
        .or(case.output.coefficients.as_ref().map(|c| c.velocity))?;
    let rho = r.density.unwrap_or(case.rho());
    Some((0.5 * rho * u * u, r.pressure))
}

/// Checks `[output] derived` and `[output.reference]`.
pub fn validate(case: &FlowCase) -> Result<()> {
    let r = &case.output.reference;
    if let Some(u) = r.velocity {
        crate::need!(
            u > 0.0 && u.is_finite(),
            CASE_RANGE,
            "[output.reference] velocity must be > 0 (m/s)"
        );
    }
    if let Some(d) = r.density {
        crate::need!(
            d > 0.0 && d.is_finite(),
            CASE_RANGE,
            "[output.reference] density must be > 0 (kg/m³)"
        );
    }
    crate::need!(
        r.pressure.is_finite(),
        CASE_RANGE,
        "[output.reference] pressure must be finite (Pa)"
    );
    if let Some(f) = &case.output.frame {
        let zones: Vec<&str> = case
            .mrf
            .iter()
            .filter_map(|m| m.zone.as_deref())
            .chain(case.sliding.iter().filter_map(|m| m.zone.as_deref()))
            .collect();
        let hint = crate::bc::did_you_mean(f, zones.iter().copied())
            .map(|z| format!(" (did you mean \"{z}\"?)"))
            .unwrap_or_default();
        crate::need!(
            zones.contains(&f.as_str()),
            CASE_INCONSISTENT,
            "[output] frame = \"{f}\" is not the `zone` of an [[mrf]] or [[sliding]] table{hint}"
        );
    }
    validate_list(case, &case.output.derived)
}

/// Checks a list of derived fields against `case` (also for
/// `foehn post derive --fields`).
pub fn validate_list(case: &FlowCase, list: &[Derived]) -> Result<()> {
    if list.contains(&Derived::WallDistance) {
        crate::need!(
            case.boundary
                .values()
                .any(|b| matches!(b, BoundaryCondition::Wall { .. })),
            CASE_INCONSISTENT,
            "derived field `wall-distance` needs at least one wall boundary"
        );
    }
    if let Some(d) = list.iter().find(|d| d.of_turbulence()) {
        crate::need!(
            !case.turbulence.is_laminar(),
            CASE_INCONSISTENT,
            "derived field `{}` needs a turbulence model ([turbulence] model)",
            d.name()
        );
    }
    if list.contains(&Derived::WallNormal) {
        crate::need!(
            case.turbulence.model.is_eb(),
            CASE_INCONSISTENT,
            "derived field `wall-normal` is the n = ∇α/|∇α| of the elliptic-blending k-ε: \
             set [turbulence] model = \"elliptic-blending-k-epsilon\""
        );
    }
    if let Some(d) = list.iter().find(|d| d.of_rotors()) {
        crate::need!(
            !case.mrf.is_empty() || !case.sliding.is_empty(),
            CASE_INCONSISTENT,
            "derived field `{}` needs a turning zone ([[mrf]] or [[sliding]], D183)",
            d.name()
        );
    }
    if list.contains(&Derived::PressureCoefficient) {
        crate::need!(
            reference(case).is_some(),
            CASE_INCONSISTENT,
            "derived field `pressure-coefficient` needs a reference velocity: set \
             [output.reference] velocity (or [output.coefficients] velocity)"
        );
    }
    Ok(())
}

/// Computed derived fields: cell fields (one value per cell) and boundary
/// face fields (one value per boundary face).
#[derive(Debug, Clone, Default)]
pub struct Fields {
    /// Cell fields, name and values.
    pub cells: Vec<(&'static str, Vec<f64>)>,
    /// Boundary face fields, name and values.
    pub faces: Vec<(&'static str, Vec<f64>)>,
}

impl Fields {
    /// The values of the field `name`.
    pub fn get(&self, name: &str) -> Option<&[f64]> {
        self.cells
            .iter()
            .chain(&self.faces)
            .find(|(n, _)| *n == name)
            .map(|(_, v)| v.as_slice())
    }

    /// The units of the datasets, `name unit; …`.
    pub fn units(list: &[Derived]) -> String {
        let mut s = String::new();
        for d in sorted(list) {
            for (n, u) in d.datasets() {
                s.push_str(&format!("; {n} {u}"));
            }
        }
        s
    }
}

/// `list` without repetitions, in the order of [`Derived::ALL`] (the file
/// does not depend on the order of the list).
pub fn sorted(list: &[Derived]) -> Vec<Derived> {
    Derived::ALL
        .into_iter()
        .filter(|d| list.contains(d))
        .collect()
}

/// The velocity gradient of cell `c`, G[i][j] = ∂uᵢ/∂xⱼ.
fn gradient(fl: &Solution, c: usize) -> [[f64; 3]; 3] {
    [
        fl.grad_u[0].get(c),
        fl.grad_u[1].get(c),
        fl.grad_u[2].get(c),
    ]
}

/// The cell quantities of one velocity gradient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kinematics {
    /// Vorticity ω (1/s).
    pub vorticity: [f64; 3],
    /// Q criterion (1/s²).
    pub q: f64,
    /// λ₂ (1/s²).
    pub lambda2: f64,
    /// Enstrophy ½ ω·ω (1/s²).
    pub enstrophy: f64,
    /// Shear rate √(2 S:S) (1/s).
    pub strain_rate: f64,
}

/// The derived quantities of the velocity gradient `g`, gᵢⱼ = ∂uᵢ/∂xⱼ.
pub fn kinematics(g: [[f64; 3]; 3]) -> Kinematics {
    let mut s = [[0.0; 3]; 3];
    let mut w = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            s[i][j] = 0.5 * (g[i][j] + g[j][i]);
            w[i][j] = 0.5 * (g[i][j] - g[j][i]);
        }
    }
    let (mut ss, mut ww) = (0.0, 0.0);
    for i in 0..3 {
        for j in 0..3 {
            ss += s[i][j] * s[i][j];
            ww += w[i][j] * w[i][j];
        }
    }
    let vorticity = [g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1]];
    // S² + Ω², symmetric
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut v = 0.0;
            for k in 0..3 {
                v += s[i][k] * s[k][j] + w[i][k] * w[k][j];
            }
            m[i][j] = v;
        }
    }
    Kinematics {
        vorticity,
        q: 0.5 * (ww - ss),
        lambda2: symmetric_eigenvalues(m)[1],
        enstrophy: 0.5 * dot(vorticity, vorticity),
        strain_rate: (2.0 * ss).sqrt(),
    }
}

/// The features of D166 of one cell (the order of the datasets of
/// [`Derived::TurbulenceFeatures`]): velocity gradient `g`, velocity `u`,
/// pressure gradient `gp`, k, ω, ν_t, wall distance `d`, ν and β*.
#[allow(clippy::too_many_arguments)]
pub fn features(
    g: [[f64; 3]; 3],
    u: [f64; 3],
    gp: [f64; 3],
    k: f64,
    omega: f64,
    nu_t: f64,
    d: f64,
    nu: f64,
    beta_star: f64,
) -> [f64; 12] {
    use foehncraft_kernels::turb::invariants as inv;
    let tiny = 1e-300;
    let tau = 1.0 / (beta_star * omega).max(tiny);
    let l = inv::pope(&g, tau);
    let (s2, w2) = inv::s2_omega2(&g);
    let k = k.max(0.0);
    let sk = k.sqrt();
    let speed = norm(u);
    let a = [
        g[0][0] * u[0] + g[0][1] * u[1] + g[0][2] * u[2],
        g[1][0] * u[0] + g[1][1] * u[1] + g[1][2] * u[2],
        g[2][0] * u[0] + g[2][1] * u[1] + g[2][2] * u[2],
    ];
    let cross = [
        u[1] * a[2] - u[2] * a[1],
        u[2] * a[0] - u[0] * a[2],
        u[0] * a[1] - u[1] * a[0],
    ];
    let length = sk * tau;
    [
        l[0],
        l[1],
        l[2],
        l[3],
        l[4],
        if d.is_finite() { d * sk / nu } else { 1e300 },
        tau * s2.sqrt(),
        nu_t / nu,
        if d.is_finite() {
            inv::shielding(&g, nu_t, nu, d)
        } else {
            1.0
        },
        (w2 - s2) / (w2 + s2).max(tiny),
        dot(u, gp) / (speed * norm(gp)).max(tiny),
        norm(cross) * length / (speed * speed * speed).max(tiny),
    ]
}

/// The eigenvalues of the symmetric matrix `a`, ascending: cyclic Jacobi
/// rotations (Golub and Van Loan, Matrix Computations, section 8.5),
/// accurate to round-off also for repeated eigenvalues, where the closed
/// form of the cubic loses half the digits.
pub fn symmetric_eigenvalues(a: [[f64; 3]; 3]) -> [f64; 3] {
    let mut m = a;
    let scale = (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .fold(0.0f64, |s, (i, j)| s + a[i][j] * a[i][j])
        .sqrt();
    for _ in 0..32 {
        let off = m[0][1] * m[0][1] + m[0][2] * m[0][2] + m[1][2] * m[1][2];
        if off.sqrt() <= f64::EPSILON * 1e-3 * scale || off == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if m[p][q] == 0.0 {
                continue;
            }
            // the rotation that zeroes m[p][q]
            let theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            let mut r = m;
            for k in 0..3 {
                r[k][p] = c * m[k][p] - s * m[k][q];
                r[k][q] = s * m[k][p] + c * m[k][q];
            }
            let mut n = r;
            for k in 0..3 {
                n[p][k] = c * r[p][k] - s * r[q][k];
                n[q][k] = s * r[p][k] + c * r[q][k];
            }
            n[p][q] = 0.0;
            n[q][p] = 0.0;
            m = n;
        }
    }
    let mut e = [m[0][0], m[1][1], m[2][2]];
    e.sort_by(f64::total_cmp);
    e
}

/// Computes the derived fields `list` of the solution `fl` (cell order and
/// boundary face order of `mesh`).
pub fn compute(
    case: &FlowCase,
    list: &[Derived],
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
) -> Result<Fields> {
    validate_list(case, list)?;
    let list = sorted(list);
    let nc = mesh.n_cells;
    let mut out = Fields::default();
    let cellwise = list.iter().any(|d| {
        matches!(
            d,
            Derived::Vorticity
                | Derived::QCriterion
                | Derived::Lambda2
                | Derived::Enstrophy
                | Derived::StrainRate
        )
    });
    let kin: Vec<Kinematics> = if cellwise {
        (0..nc).map(|c| kinematics(gradient(fl, c))).collect()
    } else {
        Vec::new()
    };
    let rho = fl.phys.rho;
    let turbo: Vec<Turbo> = match list.iter().any(|d| d.of_rotors()) {
        true => {
            let frames = Frames::new(case).ok_or_else(|| {
                crate::diag::err(crate::diag::Diagnostic::new(
                    crate::diag::CASE_INCONSISTENT,
                    "the fields of rotating machinery need a turning zone ([[mrf]] or [[sliding]])",
                ))
            })?;
            (0..nc)
                .map(|c| frames.at(mesh, c, geo.cell_centroid[c], fl.u.get(c)))
                .collect()
        }
        false => Vec::new(),
    };
    let vector = |out: &mut Fields, names: [&'static str; 4], v: &dyn Fn(&Turbo) -> [f64; 3]| {
        for (i, n) in names.into_iter().take(3).enumerate() {
            out.cells.push((n, turbo.iter().map(|t| v(t)[i]).collect()));
        }
        out.cells
            .push((names[3], turbo.iter().map(|t| norm(v(t))).collect()));
    };
    let nb = mesh.n_boundary_faces();
    let walls = (list.contains(&Derived::WallShearStress) || list.contains(&Derived::YPlus))
        .then(|| wall_values(case, mesh, geo, fl));
    for d in list {
        match d {
            Derived::Vorticity => {
                for (i, n) in ["vorticity_x", "vorticity_y", "vorticity_z"]
                    .into_iter()
                    .enumerate()
                {
                    out.cells
                        .push((n, kin.iter().map(|k| k.vorticity[i]).collect()));
                }
                out.cells.push((
                    "vorticity_magnitude",
                    kin.iter().map(|k| norm(k.vorticity)).collect(),
                ));
            }
            Derived::QCriterion => out
                .cells
                .push(("q_criterion", kin.iter().map(|k| k.q).collect())),
            Derived::Lambda2 => out
                .cells
                .push(("lambda2", kin.iter().map(|k| k.lambda2).collect())),
            Derived::Enstrophy => out
                .cells
                .push(("enstrophy", kin.iter().map(|k| k.enstrophy).collect())),
            Derived::StrainRate => out
                .cells
                .push(("strain_rate", kin.iter().map(|k| k.strain_rate).collect())),
            Derived::PressureCoefficient => {
                let (q, p0) = reference(case).expect("checked by validate_list");
                out.cells.push((
                    "pressure_coefficient",
                    fl.p.iter().map(|&p| (p - p0) / q).collect(),
                ));
            }
            Derived::TotalPressure => out.cells.push((
                "total_pressure",
                (0..nc)
                    .map(|c| {
                        let u = fl.u.get(c);
                        fl.p[c] + 0.5 * rho * dot(u, u)
                    })
                    .collect(),
            )),
            Derived::RelativeVelocity => vector(
                &mut out,
                [
                    "relative_velocity_x",
                    "relative_velocity_y",
                    "relative_velocity_z",
                    "relative_velocity_magnitude",
                ],
                &|t| t.w,
            ),
            Derived::FrameVelocity => vector(
                &mut out,
                [
                    "frame_velocity_x",
                    "frame_velocity_y",
                    "frame_velocity_z",
                    "frame_velocity_magnitude",
                ],
                &|t| t.u,
            ),
            Derived::Turbomachinery => {
                let get: [fn(&Turbo) -> f64; 7] = [
                    |t| t.c_r,
                    |t| t.c_u,
                    |t| t.c_a,
                    |t| t.c_m,
                    |t| t.w_u,
                    |t| t.alpha,
                    |t| t.beta,
                ];
                for ((n, _), g) in d.datasets().iter().zip(get) {
                    out.cells.push((n, turbo.iter().map(g).collect()));
                }
            }
            Derived::RelativeTotalPressure => out.cells.push((
                "relative_total_pressure",
                turbo
                    .iter()
                    .zip(&fl.p)
                    .map(|(t, p)| p + 0.5 * rho * dot(t.w, t.w))
                    .collect(),
            )),
            Derived::Rothalpy => out.cells.push((
                "rothalpy",
                turbo
                    .iter()
                    .zip(&fl.p)
                    .map(|(t, p)| p / rho + 0.5 * dot(t.w, t.w) - 0.5 * dot(t.u, t.u))
                    .collect(),
            )),
            Derived::WallDistance => out
                .cells
                .push(("wall_distance", wall_distance(case, mesh, geo)?)),
            Derived::WallShearStress => {
                let w = walls.as_ref().expect("computed above");
                for (i, n) in [
                    "wall_shear_stress_x",
                    "wall_shear_stress_y",
                    "wall_shear_stress_z",
                ]
                .into_iter()
                .enumerate()
                {
                    out.faces.push((n, w.iter().map(|v| v.0[i]).collect()));
                }
                out.faces.push((
                    "wall_shear_stress_magnitude",
                    w.iter().map(|v| norm(v.0)).collect(),
                ));
            }
            Derived::YPlus => {
                let w = walls.as_ref().expect("computed above");
                out.faces.push(("y_plus", w.iter().map(|v| v.1).collect()));
            }
            Derived::TurbulenceFeatures => {
                let t = fl.turbulence.as_ref().ok_or_else(|| {
                    crate::diag::err(crate::diag::Diagnostic::new(
                        crate::diag::CASE_INCONSISTENT,
                        "derived field `turbulence-features` needs the turbulence fields",
                    ))
                })?;
                // the k-ε family: ω = ε/(C_μ k), so that τ = 1/(β* ω) = k/ε
                let beta_star = case.turbulence.c_mu()?;
                let nu = fl.phys.mu / fl.phys.rho;
                let rows: Vec<[f64; 12]> = (0..nc)
                    .map(|c| {
                        let gp = fl.grad_p.get(c);
                        features(
                            gradient(fl, c),
                            fl.u.get(c),
                            gp,
                            t.k[c],
                            if t.epsilon {
                                t.omega[c] / (beta_star * t.k[c]).max(1e-300)
                            } else {
                                t.omega[c]
                            },
                            t.nu_t[c],
                            t.wall_distance.get(c).copied().unwrap_or(f64::INFINITY),
                            nu,
                            beta_star,
                        )
                    })
                    .collect();
                for (i, (n, _)) in d.datasets().iter().enumerate() {
                    out.cells.push((n, rows.iter().map(|r| r[i]).collect()));
                }
            }
            Derived::TurbulenceScales => {
                let t = turbulence_state(fl, d)?;
                let c_mu = case.turbulence.c_mu()?;
                let nu = fl.phys.mu / fl.phys.rho;
                let rows: Vec<[f64; 5]> = (0..nc)
                    .map(|c| {
                        let zt = t.zeta.as_ref().map(|z| z.t[c]);
                        scales(t.k[c], t.omega[c], t.epsilon, t.nu_t[c], zt, c_mu, nu)
                    })
                    .collect();
                for (i, (n, _)) in d.datasets().iter().enumerate() {
                    out.cells.push((n, rows.iter().map(|r| r[i]).collect()));
                }
            }
            Derived::ReynoldsStress => {
                let t = turbulence_state(fl, d)?;
                let rows: Vec<[f64; 6]> = (0..nc)
                    .map(|c| reynolds_stress(gradient(fl, c), t.k[c], t.nu_t[c]))
                    .collect();
                for (i, (n, _)) in d.datasets().iter().enumerate() {
                    out.cells.push((n, rows.iter().map(|r| r[i]).collect()));
                }
            }
            Derived::WallNormal => {
                let t = turbulence_state(fl, d)?;
                let z = t.zeta.as_ref().ok_or_else(|| {
                    crate::diag::err(crate::diag::Diagnostic::new(
                        crate::diag::CASE_INCONSISTENT,
                        "derived field `wall-normal` needs the α field of the elliptic-blending k-ε",
                    ))
                })?;
                let n = wall_normal(case, mesh, fl, &z.f);
                for (i, (name, _)) in d.datasets().iter().enumerate() {
                    out.cells.push((name, n.iter().map(|v| v[i]).collect()));
                }
            }
            Derived::CurvatureCorrection => {
                let t = turbulence_state(fl, d)?;
                let f = t.fr1.as_ref().ok_or_else(|| {
                    crate::diag::err(crate::diag::Diagnostic::new(
                        crate::diag::CASE_INCONSISTENT,
                        "derived field `curvature-correction` needs `turbulence.curvature_correction = true`",
                    ))
                })?;
                out.cells.push(("rotation_function", f.clone()));
            }
            Derived::TurbulenceWall => {
                let t = turbulence_state(fl, d)?;
                let c_mu = case.turbulence.c_mu()?;
                let u_star = c_mu.sqrt().sqrt();
                let ni = mesh.n_interior_faces();
                let mut ut = vec![0.0; nb];
                let mut ys = vec![0.0; nb];
                for p in &mesh.boundaries {
                    if !matches!(
                        case.boundary.get(&p.name),
                        Some(BoundaryCondition::Wall { .. })
                    ) {
                        continue;
                    }
                    for f in p.range() {
                        let s = geo.face_area[f];
                        let c = fl.disc.owner[f] as usize;
                        let y =
                            dot(sub(geo.face_centroid[f], geo.cell_centroid[c]), s).abs() / norm(s);
                        let nu = fl.phys.mu_at(f) / fl.phys.rho;
                        ut[f - ni] = t.u_tau[f - ni];
                        ys[f - ni] = y * u_star * t.k[c].max(0.0).sqrt() / nu;
                    }
                }
                let bl = t.wall_blending.clone().ok_or_else(|| {
                    crate::diag::err(crate::diag::Diagnostic::new(
                        crate::diag::CASE_INCONSISTENT,
                        "derived field `turbulence-wall`: the restart file has no wall blending \
                         weight (`turb_wall_blending`, written since 0.6); continue the run for \
                         one iteration to write it",
                    ))
                })?;
                out.faces.push(("friction_velocity", ut));
                out.faces.push(("y_star", ys));
                out.faces.push(("blending", bl));
            }
        }
    }
    debug_assert!(out.faces.iter().all(|(_, v)| v.len() == nb));
    Ok(out)
}

/// The turbulence fields of `fl`, for the derived field `d`.
fn turbulence_state<'a>(fl: &'a Solution, d: Derived) -> Result<&'a crate::turbulence::TurbState> {
    fl.turbulence.as_ref().ok_or_else(|| {
        crate::diag::err(crate::diag::Diagnostic::new(
            crate::diag::CASE_INCONSISTENT,
            format!("derived field `{}` needs the turbulence fields", d.name()),
        ))
    })
}

/// The scales of `turbulence-scales` of one cell (the order of its
/// datasets): ε, ω, T, L = k^{3/2}/ε and ν_t/ν, from k, the second
/// variable `w` (ε when `epsilon`, else ω), ν_t, the bounded time scale
/// `zeta_t` of ζ-f and EB, C_μ (β* for SST) and ν.
pub fn scales(
    k: f64,
    w: f64,
    epsilon: bool,
    nu_t: f64,
    zeta_t: Option<f64>,
    c_mu: f64,
    nu: f64,
) -> [f64; 5] {
    let tiny = 1e-300;
    let k = k.max(0.0);
    let (eps, omega) = if epsilon {
        (w, w / (c_mu * k).max(tiny))
    } else {
        (c_mu * k * w, w)
    };
    let t = zeta_t.unwrap_or(if epsilon {
        k / eps.max(tiny)
    } else {
        1.0 / (c_mu * omega).max(tiny)
    });
    [eps, omega, t, k * k.sqrt() / eps.max(tiny), nu_t / nu]
}

/// The Boussinesq Reynolds stress ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ of the velocity
/// gradient `g`, components xx, yy, zz, xy, yz, xz.
pub fn reynolds_stress(g: [[f64; 3]; 3], k: f64, nu_t: f64) -> [f64; 6] {
    let s = |i: usize, j: usize| 0.5 * (g[i][j] + g[j][i]);
    let iso = 2.0 / 3.0 * k.max(0.0);
    [
        iso - 2.0 * nu_t * s(0, 0),
        iso - 2.0 * nu_t * s(1, 1),
        iso - 2.0 * nu_t * s(2, 2),
        -2.0 * nu_t * s(0, 1),
        -2.0 * nu_t * s(1, 2),
        -2.0 * nu_t * s(0, 2),
    ]
}

/// n = ∇α/|∇α| per cell (0 where the gradient vanishes): the least-squares
/// gradient of the velocity stencil (D24) with α = 0 on the walls and a
/// zero gradient on every other boundary.
fn wall_normal(case: &FlowCase, mesh: &Mesh, fl: &Solution, alpha: &[f64]) -> Vec<[f64; 3]> {
    let d = fl.disc;
    let ni = d.n_interior;
    let mut wall = vec![false; mesh.n_boundary_faces()];
    for p in &mesh.boundaries {
        if matches!(
            case.boundary.get(&p.name),
            Some(BoundaryCondition::Wall { .. })
        ) {
            for f in p.range() {
                wall[f - ni] = true;
            }
        }
    }
    (0..d.n_cells)
        .map(|c| {
            let pc = alpha[c];
            let mut g = [0.0; 3];
            for k in d.cell_ptr[c] as usize..d.cell_ptr[c + 1] as usize {
                let (f, is_nb) = foehncraft_kernels::unpack_face(d.cell_faces[k]);
                let diff = if f < ni {
                    let o = if is_nb { d.owner[f] } else { d.neighbour[f] };
                    alpha[o as usize] - pc
                } else if wall[f - ni] {
                    -pc
                } else {
                    0.0
                };
                let w = d.lsq.get(k);
                for i in 0..3 {
                    g[i] += w[i] * diff;
                }
            }
            let m = norm(g);
            if m > 1e-300 {
                g.map(|v| v / m)
            } else {
                [0.0; 3]
            }
        })
        .collect()
}

/// The exact wall distance of every cell centroid (D48): the closest point
/// on the faces of the wall boundaries of `case`, by the BVH search of
/// D108.
pub fn wall_distance(case: &FlowCase, mesh: &Mesh, geo: &Geometry) -> Result<Vec<f64>> {
    let walls: Vec<usize> = mesh
        .boundaries
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            matches!(
                case.boundary.get(&p.name),
                Some(BoundaryCondition::Wall { .. })
            )
        })
        .map(|(i, _)| i)
        .collect();
    let d = foehncraft_mesh::query::wall_distance(mesh, geo, &walls).ok_or_else(|| {
        crate::diag::err(crate::diag::Diagnostic::new(
            crate::diag::CASE_INCONSISTENT,
            "derived field `wall-distance` needs at least one wall face in the mesh",
        ))
    })?;
    Ok(d.into_iter().map(|c| c.distance).collect())
}

/// What [`from_restart`] wrote.
#[derive(Debug, Clone, Serialize)]
pub struct PostReport {
    /// The file written.
    pub file: std::path::PathBuf,
    /// The restart file read.
    pub restart: std::path::PathBuf,
    /// Cells.
    pub cells: usize,
    /// Outer iterations of the solution.
    pub iteration: usize,
    /// Time step of a transient solution.
    pub step: Option<usize>,
    /// Simulated time of a transient solution (s).
    pub time: Option<f64>,
    /// The derived fields written.
    pub derived: Vec<Derived>,
    /// Their datasets.
    pub datasets: Vec<&'static str>,
}

/// `foehn post derive` (D72): reads the solution of `restart`, a restart
/// file of `case`, and writes the fields with the derived fields `list` to
/// `file`, as the run would have written them (the same file as
/// `fields.vtkhdf` of a run with `[output] derived = list` at the same
/// iteration or step).
pub fn from_restart(
    case: &FlowCase,
    restart: &std::path::Path,
    list: &[Derived],
    file: &std::path::Path,
) -> Result<PostReport> {
    validate_list(case, list)?;
    with_restart(case, restart, |l| {
        let mut c = case.clone();
        c.output.vtkhdf = true;
        c.output.derived = sorted(list);
        let at = match l.clock {
            Some(k) => Some((k.step, Some(k.time))),
            None => Some((l.iteration, None)),
        };
        crate::case::write_fields(&c, file, l.mesh, l.geo, l.sol, l.analytic.as_ref(), at)?;
        Ok(PostReport {
            file: file.to_path_buf(),
            restart: restart.to_path_buf(),
            cells: l.mesh.n_cells,
            iteration: l.iteration,
            step: l.clock.map(|k| k.step),
            time: l.clock.map(|k| k.time),
            datasets: c
                .output
                .derived
                .iter()
                .flat_map(|d| d.datasets().iter().map(|(n, _)| *n))
                .collect(),
            derived: c.output.derived,
        })
    })
}

/// The derived fields `list` of the solution in `restart`, a restart file
/// of `case`, with the mesh and geometry they belong to (for checks and
/// scripts; [`from_restart`] writes them).
pub fn evaluate(
    case: &FlowCase,
    restart: &std::path::Path,
    list: &[Derived],
) -> Result<(Mesh, Geometry, Fields)> {
    validate_list(case, list)?;
    with_restart(case, restart, |l| {
        let f = compute(case, list, l.mesh, l.geo, l.sol)?;
        Ok((l.mesh.clone(), l.geo.clone(), f))
    })
}

/// The [`WallExchange`] of the solution in `restart`, a restart file of
/// `case` (D124 W4; for checks and for the wall models of later stages).
pub fn wall_exchange_from_restart(
    case: &FlowCase,
    restart: &std::path::Path,
) -> Result<(Mesh, Geometry, WallExchange)> {
    with_restart(case, restart, |l| {
        Ok((
            l.mesh.clone(),
            l.geo.clone(),
            wall_exchange(case, l.mesh, l.geo, l.sol),
        ))
    })
}

/// A solution read back from a restart file.
struct Loaded<'a> {
    mesh: &'a Mesh,
    geo: &'a Geometry,
    sol: &'a Solution<'a>,
    analytic: Option<crate::analytic::Analytic>,
    iteration: usize,
    clock: Option<crate::transient::Clock>,
}

/// Reads `restart`, builds the mesh and the discretisation of `case` and
/// hands the solution to `f`.
fn with_restart<R>(
    case: &FlowCase,
    restart: &std::path::Path,
    f: impl FnOnce(&Loaded) -> Result<R>,
) -> Result<R> {
    use anyhow::{bail, ensure};
    ensure!(
        crate::restart::AVAILABLE,
        "post derive reads restart files and writes VTKHDF: this build has no HDF5 (feature `hdf5`)"
    );
    let level = match case.verification.levels.as_slice() {
        [] => None,
        [l] => Some(*l),
        _ => bail!(
            "a restart file holds a single run; this case is a sweep over levels \
             (choose one with --set verification.levels=[N])"
        ),
    };
    let (mesh, geo, _) = crate::case::mesh_for(case, level)?;
    let nc = mesh.n_cells;
    let r = crate::restart::read(restart)?;
    ensure!(
        r.geometry.is_none(),
        "{}: restart files of moving meshes are not read by post derive yet",
        restart.display()
    );
    let (r, _) = r.to_local(&mesh.cell_id[..nc], &mesh.face_id, mesh.n_interior_faces())?;
    let clock = r.transient.as_ref().map(|t| t.1);
    match (clock.is_some(), case.time.is_some()) {
        (true, false) => bail!(
            "{} is the restart file of a transient run; the case is steady (no [time] table)",
            restart.display()
        ),
        (false, true) => bail!(
            "{} is the restart file of a steady run; the case is transient",
            restart.display()
        ),
        _ => {}
    }
    let analytic = match clock {
        Some(c) => case.analytic_at(c.time),
        None => case.analytic(),
    };
    let phys = crate::case::physics(case)?;
    let phys = crate::case::physics_on(case, &phys, &mesh, &geo)?;
    let analytic0 = match &case.time {
        Some(t) => case.analytic_at(t.start_time),
        None => analytic,
    };
    let values = case.compiled_boundaries()?;
    let t_now = clock.map_or(case.time.map_or(0.0, |t| t.start_time), |c| c.time);
    let arrays = crate::bc::arrays(
        &mesh,
        &geo,
        &case.boundary,
        phys.rho,
        crate::bc::At {
            analytic: analytic0.as_ref(),
            values: &values,
            t: t_now,
            turning: &case.turning(),
            hydrostatic: crate::buoyancy::hydrostatic(case),
        },
    )?;
    let (disc, _) = crate::simplec::discretisation(&mesh, &geo, &arrays)?;
    // the driving of the periodic pairs (D29) is part of the body force
    let mut phys = phys;
    for (d, b) in arrays.drives.iter().zip(&r.state.periodic_force) {
        for (f, t) in phys.body_force.iter_mut().zip(d.direction) {
            *f += b * t;
        }
    }
    let sol = Solution {
        state: r.state,
        disc: &disc,
        phys,
    };
    f(&Loaded {
        mesh: &mesh,
        geo: &geo,
        sol: &sol,
        analytic,
        iteration: r.iteration,
        clock,
    })
}

/// What the wall treatment gives a module coupled to the walls (D124 W4,
/// the wall film of D123 and later CHT): per wall face the gas-side shear
/// stress and the transfer coefficients of the wall treatment.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct WallExchange {
    /// The wall faces (mesh face indices), boundary by boundary.
    pub faces: Vec<usize>,
    /// Shear stress of the fluid on the wall, tangential (Pa).
    pub shear: Vec<[f64; 3]>,
    /// The momentum transfer coefficient μ_f/Δn of the wall treatment
    /// (Pa s/m): τ = coefficient · (u_P − u_wall) tangential for a
    /// laminar wall; the turbulence models of stage 6 replace it by their
    /// wall functions.
    pub momentum: Vec<f64>,
    /// The heat transfer coefficient (W/(m² K)); with the energy equation
    /// (L1).
    pub heat: Option<Vec<f64>>,
    /// The mass transfer coefficient (kg/(m² s)); with species (L2).
    pub mass: Option<Vec<f64>>,
}

/// The [`WallExchange`] of a solution (D124 W4).
pub fn wall_exchange(case: &FlowCase, mesh: &Mesh, geo: &Geometry, fl: &Solution) -> WallExchange {
    let ni = mesh.n_interior_faces();
    let mut w = WallExchange::default();
    for p in &mesh.boundaries {
        if !matches!(
            case.boundary.get(&p.name),
            Some(BoundaryCondition::Wall { .. })
        ) {
            continue;
        }
        for f in p.range() {
            let (_, _, tau) = crate::case::wall_face_values(fl, geo, f, ni);
            let s = geo.face_area[f];
            let c = fl.disc.owner[f] as usize;
            let dn = dot(sub(geo.face_centroid[f], geo.cell_centroid[c]), s).abs() / norm(s);
            w.faces.push(f);
            w.shear.push(tau);
            w.momentum.push(fl.mu_eff_at(f) / dn);
        }
    }
    w
}

/// Wall shear stress and y⁺ per boundary face (zero off the walls).
fn wall_values(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    fl: &Solution,
) -> Vec<([f64; 3], f64)> {
    let ni = mesh.n_interior_faces();
    let mut v = vec![([0.0; 3], 0.0); mesh.n_boundary_faces()];
    for p in &mesh.boundaries {
        if !matches!(
            case.boundary.get(&p.name),
            Some(BoundaryCondition::Wall { .. })
        ) {
            continue;
        }
        for f in p.range() {
            let (_, _, tau) = crate::case::wall_face_values(fl, geo, f, ni);
            let s = geo.face_area[f];
            let c = fl.disc.owner[f] as usize;
            let y = dot(sub(geo.face_centroid[f], geo.cell_centroid[c]), s).abs() / norm(s);
            let u_tau = (norm(tau) / fl.phys.rho).sqrt();
            let nu = fl.phys.mu_at(f) / fl.phys.rho;
            // with a turbulence model: the y⁺ of its wall treatment
            let yp = match &fl.state.turbulence {
                Some(t) => t.y_plus[f - ni],
                None => y * u_tau / nu,
            };
            v[f - ni] = (tau, yp);
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * (1.0 + b.abs())
    }

    #[test]
    fn solid_body_rotation() {
        // u = −Ω y, v = Ω x
        let om = 3.0;
        let k = kinematics([[0.0, -om, 0.0], [om, 0.0, 0.0], [0.0; 3]]);
        assert_eq!(k.vorticity, [0.0, 0.0, 2.0 * om]);
        assert!(close(k.q, om * om));
        assert!(close(k.lambda2, -om * om));
        assert!(close(k.enstrophy, 2.0 * om * om));
        assert!(close(k.strain_rate, 0.0));
    }

    #[test]
    fn simple_shear_is_no_vortex() {
        // u = γ y: as much strain as rotation
        let g = 2.5;
        let k = kinematics([[0.0, g, 0.0], [0.0; 3], [0.0; 3]]);
        assert_eq!(k.vorticity, [0.0, 0.0, -g]);
        assert!(close(k.q, 0.0));
        assert!(close(k.lambda2, 0.0));
        assert!(close(k.strain_rate, g));
        assert!(close(k.enstrophy, 0.5 * g * g));
    }

    #[test]
    fn pure_strain() {
        // u = a x, v = −a y: Q = −a², λ₂ = 0 (eigenvalues a², a², 0)
        let a = 1.5;
        let k = kinematics([[a, 0.0, 0.0], [0.0, -a, 0.0], [0.0; 3]]);
        assert_eq!(k.vorticity, [0.0; 3]);
        assert!(close(k.q, -a * a));
        assert!(close(k.lambda2, a * a));
        assert!(close(k.strain_rate, 2.0 * a));
    }

    #[test]
    fn eigenvalues_of_symmetric_matrices() {
        // A = R diag(e) Rᵀ with a rotation R
        let (c, s) = (0.3f64.cos(), 0.3f64.sin());
        let (c2, s2) = (1.1f64.cos(), 1.1f64.sin());
        let rz = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let rx = [[1.0, 0.0, 0.0], [0.0, c2, -s2], [0.0, s2, c2]];
        let mul = |a: [[f64; 3]; 3], b: [[f64; 3]; 3]| {
            let mut m = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        m[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
            m
        };
        let r = mul(rz, rx);
        let mut rt = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                rt[i][j] = r[j][i];
            }
        }
        for e in [
            [-2.0, 0.5, 4.0],
            [1.0, 1.0, 3.0],
            [-1e-3, -1e-3, -1e-3],
            [0.0, 0.0, 7.0],
        ] {
            let d = [[e[0], 0.0, 0.0], [0.0, e[1], 0.0], [0.0, 0.0, e[2]]];
            let a = mul(mul(r, d), rt);
            let got = symmetric_eigenvalues(a);
            for i in 0..3 {
                assert!((got[i] - e[i]).abs() < 1e-14 * 8.0, "{got:?} {e:?}");
            }
        }
    }

    /// The wall distance in a channel between walls at y = 0 and y = 1
    /// (the ends are an inlet and an outlet, which do not count): min(y,
    /// 1 − y) at every centroid, exactly on hex cells; the polyhedral mesh
    /// gives the distance to its (plane) wall faces, the same function.
    #[test]
    fn wall_distance_of_a_channel() {
        for mesh in ["hex", "poly", "tet"] {
            let case = FlowCase::from_toml(
                &format!(
                    r#"
[mesh.generate]
type = "{mesh}"
dimension = 2
upper = [2.0, 1.0, 0.1]
cells = [8, 6, 1]
[fluid]
viscosity = 0.1
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[output]
derived = ["wall-distance"]
"#
                ),
                std::path::Path::new("."),
            )
            .unwrap();
            validate(&case).unwrap();
            let (m, g, _) = crate::case::mesh_for(&case, None).unwrap();
            let d = wall_distance(&case, &m, &g).unwrap();
            assert_eq!(d.len(), m.n_cells);
            for (c, &dc) in d.iter().enumerate() {
                let y = g.cell_centroid[c][1];
                let e = y.min(1.0 - y);
                assert!((dc - e).abs() < 1e-12, "{mesh}: cell {c}, {dc} vs {e}");
            }
        }
    }

    #[test]
    fn names_round_trip() {
        for d in Derived::ALL {
            assert_eq!(Derived::parse(d.name()), Some(d));
            let t: toml::Value = toml::Value::try_from(d).unwrap();
            assert_eq!(t.as_str(), Some(d.name()));
        }
    }

    /// The identities of the turbomachinery values (D183) on random
    /// velocities, points and frames: c = w + u, c_r² + c_u² + c_a² = |c|²,
    /// c_u − w_u = |Ω| r (the blade speed, positive for either sense of
    /// rotation), the angles from their definitions, and the fixed frame.
    #[test]
    fn turbomachinery_identities() {
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
        };
        let tol = |a: f64, b: f64, s: f64| (a - b).abs() <= 1e-12 * (1.0 + s);
        for _ in 0..200 {
            let axis = [rnd(), rnd(), rnd()];
            let n = norm(axis);
            let w0 = 50.0 * rnd();
            let omega = axis.map(|a| w0 * a / n);
            let t = crate::bc::Turning {
                zone: None,
                omega,
                origin: [rnd(), rnd(), rnd()],
                axis: omega.map(|o| o / w0.abs()),
            };
            let c = [3.0 * rnd(), 3.0 * rnd(), 3.0 * rnd()];
            let x = [rnd(), rnd(), rnd()];
            let v = turbo(c, x, Some(&t), &t);
            let s = norm(c) + norm(v.u);
            for ((w, u), c) in v.w.iter().zip(v.u).zip(c) {
                assert!(tol(w + u, c, s), "c = w + u");
            }
            assert!(tol(v.c_r.hypot(v.c_u).hypot(v.c_a), norm(c), s), "|c|");
            assert!(tol(v.c_m, v.c_r.hypot(v.c_a), s));
            // the blade speed |Ω| r, r the distance from the axis
            let d = sub(x, t.origin);
            let r = norm(foehncraft_mesh::geometry::cross(t.axis, d));
            assert!(tol(v.c_u - v.w_u, w0.abs() * r, s), "u_u");
            assert!(tol(norm(v.u), w0.abs() * r, s), "|u|");
            assert!(tol(v.alpha, v.c_m.atan2(v.c_u).to_degrees(), 180.0));
            assert!(tol(v.beta, v.c_m.atan2(-v.w_u).to_degrees(), 180.0));
            assert!((0.0..=180.0).contains(&v.alpha) && (0.0..=180.0).contains(&v.beta));
            // the fixed frame: w = c, u = 0, the same components
            let f = turbo(c, x, None, &t);
            assert_eq!(f.w, c);
            assert_eq!(f.u, [0.0; 3]);
            assert_eq!((f.c_r, f.c_u, f.c_a), (v.c_r, v.c_u, v.c_a));
            assert!(tol(f.w_u, f.c_u, s));
        }
        // a fluid turning with its frame (the negative sense): w = 0, a
        // positive swirl c_u = |Ω| r, α = 0
        let t = crate::bc::Turning {
            zone: None,
            omega: [0.0, 0.0, -4.0],
            origin: [0.0; 3],
            axis: [0.0, 0.0, -1.0],
        };
        let x = [0.3, -0.2, 0.1];
        let v = turbo(t.at(x), x, Some(&t), &t);
        assert!(norm(v.w) < 1e-15);
        assert!(tol(v.c_u, 4.0 * 0.3f64.hypot(0.2), 1.0), "{}", v.c_u);
        assert!(tol(v.alpha, 0.0, 1.0), "{}", v.alpha);
        // without relative swirl β = 90°
        let c = foehncraft_mesh::geometry::add(t.at(x), [0.0, 0.0, 2.0]);
        let v = turbo(c, x, Some(&t), &t);
        assert!(tol(v.beta, 90.0, 90.0), "{}", v.beta);
        // on the axis: no radial or circumferential direction
        let v = turbo([1.0, 2.0, 3.0], [0.0, 0.0, 5.0], Some(&t), &t);
        assert_eq!((v.c_r, v.c_u), (0.0, 0.0));
        assert!(tol(v.c_a, -3.0, 1.0));
    }
}
