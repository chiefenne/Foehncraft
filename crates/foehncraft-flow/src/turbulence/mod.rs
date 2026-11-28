//! RANS turbulence models on the flow solver (stage 6, D129–D139): the
//! k-ω SST model, 2003 version (D130, stage 6 PR 6b), with the production
//! options of D138, the low-Re and automatic wall treatments and rough
//! walls (D134); the standard and the realizable k-ε model (D131, PR 6c)
//! with scalable, standard and compound wall treatments and the rough-wall
//! shift of Cebeci & Bradshaw; the inlet values (D137), the modified
//! pressure p̂ = p + ⅔ρk and the exact or Poisson wall distance (D48).
//!
//! The two-equation models share one machinery: k and a second variable
//! (ω of SST, ε of the k-ε family), each on the scalar path of stage 6a;
//! the model enters through its cell kernel ([`ksst::sst_cell`] or
//! [`kkeps::keps_cell`]) and its wall-face kernel.
//!
//! Case file:
//!
//! ```toml
//! [turbulence]
//! model = "k-omega-sst"
//! wall_treatment = "automatic"        # or "low-re"
//! production_limiter = "menter"       # "kato-launder", "durbin", "none"
//! # ambient = { k = 1e-6, omega = 1.0 }  # sustaining terms
//! # wall_distance = "pde"             # Poisson distance instead of the exact one
//! # constants = { a1 = 0.31 }         # TURB007
//!
//! [initial]
//! turbulence_intensity = 0.01          # or k and omega (numbers or expressions)
//! viscosity_ratio = 1.0
//!
//! [boundary.inlet]
//! type = "velocity-inlet"
//! velocity = [10.0, 0.0, 0.0]
//! turbulence_intensity = 0.05          # with turbulence_length_scale or viscosity_ratio,
//! viscosity_ratio = 10.0               # or k and omega
//!
//! [solver.turbulence]
//! relaxation = 0.8
//! ```
//!
//! Coupling (stage 6 specification 6.4): after the pressure correction of
//! every outer iteration (and of every step of a transient run), with the
//! velocity gradient of the corrected velocity:
//!
//! 1. inlet values; gradients and limiters of k and ω;
//! 2. the SST sources per cell (F₁, F₂, CD_kω, P̃_k and the split of
//!    section 6.2), the wall treatment per wall face and the ω of the
//!    near-wall cells;
//! 3. the face diffusivities, then the k equation and the ω equation (the
//!    near-wall rows replaced by the fixed ω), the positivity floors;
//! 4. ν_t from the new fields, μ_eff = μ + ρν_t on the faces (the wall
//!    faces take μ_w of the wall treatment), the pressure boundaries shifted
//!    to p̂.
//!
//! The momentum equation then takes the face viscosity μ_eff and the
//! transposed part ∇·[μ_eff (∇u)ᵀ] of the stress as a deferred term; the
//! pressure of the solver is p̂, the output and the forces take the static
//! pressure p = p̂ − ⅔ρk ([`crate::simplec::Flow::state`]).

pub mod inlet;
#[cfg(test)]
mod mms;
pub mod wall;

use crate::bc::{BoundaryArrays, BoundaryCondition};
use crate::device::{
    Buf, Buf3, Eps2, FlowDevice, Lsq, Mask, ScalarFaceArgs, ScalarRowArgs, ScalarTime, SstCellOps,
    TurbFaceOps, WallCellOps, WallFaceOps,
};
use crate::device::{CurvatureOps, KepsCellOps, ZetaFCellOps};
use crate::diag::{self, Diagnostic};
use crate::expr::{self, Env, Scalar, Vars};
use anyhow::Result;
use foehncraft_fv::problem::BoundaryKind;
use foehncraft_fv::solver::{LimiterKind, Preconditioner, Scheme};
use foehncraft_kernels::turb::keps::{self as kkeps, KepsConstants, KepsOptions};
use foehncraft_kernels::turb::sst::{self as ksst, SstConstants, SstOptions};
use foehncraft_kernels::turb::wall as kwall;
use foehncraft_kernels::turb::zeta_f::{self as kzf, ZetaFConstants, ZetaFOptions};
use foehncraft_mesh::{Geometry, Mesh};
use inlet::InletSpec;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use wall::Walls;

/// `turbulence.model`.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum TurbulenceModel {
    /// No turbulence model.
    #[default]
    Laminar,
    /// k-ω SST, 2003 version (Menter, Kuntz & Langtry 2003).
    KOmegaSst,
    /// Standard k-ε (Launder & Spalding 1974).
    KEpsilon,
    /// Realizable k-ε (Shih, Liou, Shabbir, Yang & Zhu 1995).
    RealizableKEpsilon,
    /// k-ε-ζ-f (Hanjalić, Popovac & Hadžiabdić 2004).
    KEpsilonZetaF,
    /// k-ε-ζ-f₀ (Popovac 2023): ζ-f with f₀ = 0 at the wall, and the
    /// compound wall treatment of the same paper by default.
    KEpsilonZetaF0,
    /// Elliptic-blending k-ε (Billard & Laurence 2012, the "BL-v²/k"
    /// model, D133): φ = v²/k and the blending parameter α on the fields
    /// and solvers of k-ε-ζ-f, the compound wall treatment by default.
    EllipticBlendingKEpsilon,
}

impl TurbulenceModel {
    /// The second variable is ε (the k-ε family and k-ε-ζ-f).
    pub fn is_epsilon(self) -> bool {
        matches!(
            self,
            Self::KEpsilon
                | Self::RealizableKEpsilon
                | Self::KEpsilonZetaF
                | Self::KEpsilonZetaF0
                | Self::EllipticBlendingKEpsilon
        )
    }

    /// The standard or the realizable k-ε model (high-Re, with wall
    /// functions).
    pub fn is_keps(self) -> bool {
        matches!(self, Self::KEpsilon | Self::RealizableKEpsilon)
    }

    /// The model solves ζ = v²/k and the elliptic f (k-ε-ζ-f, ζ-f₀), or
    /// φ = v²/k and the elliptic α in their places (elliptic blending).
    pub fn has_zeta(self) -> bool {
        matches!(
            self,
            Self::KEpsilonZetaF | Self::KEpsilonZetaF0 | Self::EllipticBlendingKEpsilon
        )
    }

    /// The elliptic-blending k-ε.
    pub fn is_eb(self) -> bool {
        self == Self::EllipticBlendingKEpsilon
    }

    /// The label of the ζ-f family in messages.
    pub fn zeta_label(self) -> &'static str {
        if self.is_eb() {
            "the elliptic-blending k-ε"
        } else {
            "k-ε-ζ-f"
        }
    }

    /// The name in the case file.
    pub fn name(self) -> &'static str {
        match self {
            Self::Laminar => "laminar",
            Self::KOmegaSst => "k-omega-sst",
            Self::KEpsilon => "k-epsilon",
            Self::RealizableKEpsilon => "realizable-k-epsilon",
            Self::KEpsilonZetaF => "k-epsilon-zeta-f",
            Self::KEpsilonZetaF0 => "k-epsilon-zeta-f0",
            Self::EllipticBlendingKEpsilon => "elliptic-blending-k-epsilon",
        }
    }
}

/// `turbulence.wall_treatment` (D134).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum WallTreatment {
    /// The y⁺-insensitive blend of the viscous and the logarithmic layer
    /// (the default of SST).
    Automatic,
    /// The integrated form to the wall (first-cell y⁺ ≲ 1): SST k_w = 0
    /// and ω_w; k-ε-ζ-f k_w = 0, ζ_w = 0, ε_w = 2νk_P/y_P²,
    /// f_w = −2νζ_P/y_P².
    LowRe,
    /// k-ε family: the log law at y* = max(y*, 11.06) (Grotjans & Menter
    /// 1998; the default of the k-ε family).
    Scalable,
    /// k-ε family: the wall functions of Launder & Spalding (1974), first
    /// cell at 30 ≤ y⁺ ≤ 300.
    Standard,
    /// The compound treatment of Popovac & Hanjalić (2007): the default of
    /// k-ε-ζ-f (any first-cell y⁺), an option of the k-ε family (claimed
    /// for 30 ≤ y⁺ ≤ 300 there).
    Compound,
    /// The compound treatment of Popovac (2023): the harmonic blend of the
    /// viscous and the log layer (n = 10/3) with the non-equilibrium
    /// factor ψ; the default of k-ε-ζ-f₀ (any first-cell y⁺), an option
    /// of k-ε-ζ-f.
    #[serde(rename = "compound-2023")]
    Compound2023,
}

impl WallTreatment {
    /// The name in the case file.
    pub fn name(self) -> &'static str {
        match self {
            Self::Automatic => "automatic",
            Self::LowRe => "low-re",
            Self::Scalable => "scalable",
            Self::Standard => "standard",
            Self::Compound => "compound",
            Self::Compound2023 => "compound-2023",
        }
    }

    /// The treatment's id in the wall-face kernel.
    fn id(self) -> u32 {
        match self {
            Self::Automatic => kwall::AUTOMATIC,
            Self::LowRe => kwall::LOW_RE,
            Self::Scalable => kwall::SCALABLE,
            Self::Standard => kwall::STANDARD,
            Self::Compound => kwall::COMPOUND,
            Self::Compound2023 => kwall::COMPOUND_2023,
        }
    }
}

/// `turbulence.production_limiter` (D138).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum ProductionLimiter {
    /// P̃_k = min(ν_t S², 10 β* k ω) (2003); for the standard k-ε the
    /// same clip min(ν_t S², 10 ε).
    Menter,
    /// Kato & Launder (1993): P_k = ν_t S Ω, then the 2003 limiter.
    KatoLaunder,
    /// Durbin (1996): ν_t ≤ C_lim k/(√6 |S|), then the 2003 limiter.
    Durbin,
    /// P_k = ν_t S², unlimited.
    None,
}

/// `turbulence.qcr` (D170): the quadratic constitutive relation of SST.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum QcrMode {
    /// The Boussinesq stress only.
    #[default]
    None,
    /// QCR2000 (Spalart 2000): −C_cr1 (O_ik τ_jk + O_jk τ_ik).
    Qcr2000,
    /// QCR2013 (Mani, Babcock, Winkler & Spalart 2013): QCR2000 and
    /// −C_cr2 μ_t √(2 S*_mn S*_mn) δ_ij.
    Qcr2013,
}

impl QcrMode {
    /// The mode of the kernel (0: off).
    pub fn kernel(self) -> u32 {
        match self {
            Self::None => 0,
            Self::Qcr2000 => foehncraft_kernels::turb::sst::QCR_2000,
            Self::Qcr2013 => foehncraft_kernels::turb::sst::QCR_2013,
        }
    }
}

/// `turbulence.wall_distance` (D48).
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum WallDistanceMethod {
    /// The exact distance to the nearest wall face (bounding-volume
    /// hierarchy).
    #[default]
    Exact,
    /// From the Poisson equation ∇²φ = −1 (Tucker 2003).
    Pde,
}

/// `turbulence.k_sink`: the sink of the k equation of k-ε-ζ-f and
/// k-ε-ζ-f₀ (stage 6 specification 3.3.2).
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum KSink {
    /// max(k/T_K, 2νk/d²), T_K = max(k/ε, C_T (ν/ε)^½), d the wall
    /// distance: ε at the wall, so that k vanishes like y² and the results
    /// converge with the mesh; k/T_K away from walls, which keeps k
    /// positive in laminar regions.
    #[default]
    Wall,
    /// ε, as printed in the papers (Hanjalić et al. 2004; Popovac 2023):
    /// holds k at its floor at a laminar leading edge.
    Epsilon,
    /// k/T_K (the ζ-f₀ default of 6d2): k vanishes like y only, results
    /// drift with refinement below y⁺ 1.
    Kolmogorov,
}

impl KSink {
    /// The kernel code.
    pub fn code(self) -> u32 {
        match self {
            Self::Wall => foehncraft_kernels::turb::zeta_f::K_SINK_WALL,
            Self::Epsilon => foehncraft_kernels::turb::zeta_f::K_SINK_EPSILON,
            Self::Kolmogorov => foehncraft_kernels::turb::zeta_f::K_SINK_KOLMOGOROV,
        }
    }
}

/// `solver.turbulence.destruction`.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Destruction {
    /// −βω² ≈ βω*² − 2βω* ω.
    #[default]
    Newton,
    /// −βω² ≈ −βω* ω.
    Picard,
}

/// `turbulence.ambient`: the sustaining terms of Spalart & Rumsey (2007).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Ambient {
    /// Ambient k (m²/s²).
    pub k: f64,
    /// Ambient ω (1/s; SST).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub omega: Option<f64>,
    /// Ambient ε (m²/s³; the k-ε family).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub epsilon: Option<f64>,
}

/// `[turbulence]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct TurbulenceSection {
    /// Turbulence model (D12, D130–D133).
    pub model: TurbulenceModel,
    /// Momentum wall treatment (D134; SST default `automatic`, k-ε family
    /// `scalable`, k-ε-ζ-f and the elliptic-blending k-ε `compound`,
    /// k-ε-ζ-f₀ `compound-2023`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wall_treatment: Option<WallTreatment>,
    /// Production form of SST (default `menter`) and of the standard k-ε
    /// (default `durbin`: the realizability bound is inactive in the log
    /// layer and keeps the impulsive start of near-wall cells bounded; the
    /// realizable k-ε takes `none` only).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub production_limiter: Option<ProductionLimiter>,
    /// C_lim of the Durbin bound.
    pub durbin_coefficient: f64,
    /// Sustaining terms (Spalart & Rumsey 2007): ambient k and ω (SST) or
    /// ε (k-ε).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ambient: Option<Ambient>,
    /// Cap of ν_t/ν (counted, TURB008).
    pub max_viscosity_ratio: f64,
    /// Smooth min and max in the primal solution (specification 2.6).
    pub smooth_limits: bool,
    /// Wall distance (D48).
    pub wall_distance: WallDistanceMethod,
    /// SST: the curvature and rotation correction of Smirnov & Menter
    /// (2009) (stage 6 specification 3.6, D170; constants `c_r1`, `c_r2`,
    /// `c_r3`).
    pub curvature_correction: bool,
    /// SST: the quadratic constitutive relation QCR2000 or QCR2013
    /// (specification 3.6, D170; constants `c_cr1`, `c_cr2`), an explicit
    /// non-linear stress in the momentum equations.
    pub qcr: QcrMode,
    /// k-ε-ζ-f and k-ε-ζ-f₀: the sink of k (default `wall`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub k_sink: Option<KSink>,
    /// Constant overrides by their literature name (SST: `beta_star`, `a1`,
    /// `kappa`, `gamma1`, `gamma2`, `beta1`, `beta2`, `sigma_k1`,
    /// `sigma_k2`, `sigma_omega1`, `sigma_omega2`, and `c_r1`, `c_r2`,
    /// `c_r3` of the curvature correction, `c_cr1`, `c_cr2` of the QCR; k-ε: `c_mu`,
    /// `c_epsilon1`, `c_epsilon2`, `sigma_k`, `sigma_epsilon`, `a0`,
    /// `kappa`; k-ε-ζ-f and ζ-f₀: `c_mu`, `c_epsilon1`,
    /// `c_epsilon1_zeta`, `c_epsilon1_power` (1 or ½), `c_epsilon2`, `c_1`,
    /// `c_2`, `sigma_k`, `sigma_epsilon`, `sigma_zeta`, `a`, `c_t`, `c_l`,
    /// `c_eta`, `kappa`, `n_zeta` (ζ-f₀ only); the elliptic-blending k-ε: `c_mu`, `c_epsilon1` …
    /// `c_epsilon4`, `c_1`, `c_2`, `sigma_k`, `sigma_epsilon`, `sigma_phi`,
    /// `a`, `c_t`, `c_l`, `c_eta`, `kappa`, where `c_epsilon4 = 1.83`
    /// makes C*_ε2 constant); TURB007.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub constants: BTreeMap<String, f64>,
}

impl Default for TurbulenceSection {
    fn default() -> Self {
        Self {
            model: TurbulenceModel::Laminar,
            wall_treatment: None,
            production_limiter: None,
            durbin_coefficient: 1.0,
            ambient: None,
            max_viscosity_ratio: 1e5,
            smooth_limits: false,
            wall_distance: WallDistanceMethod::Exact,
            curvature_correction: false,
            qcr: QcrMode::None,
            k_sink: None,
            constants: BTreeMap::new(),
        }
    }
}

/// The names of the SST constants that `[turbulence.constants]` takes.
pub const SST_CONSTANTS: &[&str] = &[
    "beta_star",
    "a1",
    "kappa",
    "gamma1",
    "gamma2",
    "beta1",
    "beta2",
    "sigma_k1",
    "sigma_k2",
    "sigma_omega1",
    "sigma_omega2",
    "c_r1",
    "c_r2",
    "c_r3",
    "c_cr1",
    "c_cr2",
];

/// The names of the k-ε constants that `[turbulence.constants]` takes.
pub const KEPS_CONSTANTS: &[&str] = &[
    "c_mu",
    "c_epsilon1",
    "c_epsilon2",
    "sigma_k",
    "sigma_epsilon",
    "a0",
    "kappa",
];

/// The names of the k-ε-ζ-f constants that `[turbulence.constants]` takes.
pub const ZETA_F_CONSTANTS: &[&str] = &[
    "c_mu",
    "c_epsilon1",
    "c_epsilon1_zeta",
    "c_epsilon2",
    "c_1",
    "c_2",
    "sigma_k",
    "sigma_epsilon",
    "sigma_zeta",
    "a",
    "c_t",
    "c_l",
    "c_eta",
    "kappa",
    "n_zeta",
    "c_epsilon1_power",
];

/// The names of the elliptic-blending k-ε constants that
/// `[turbulence.constants]` takes (`a`: the coefficient of the
/// realizability bound of T in ν_t).
pub const EB_CONSTANTS: &[&str] = &[
    "c_mu",
    "c_epsilon1",
    "c_epsilon2",
    "c_epsilon3",
    "c_epsilon4",
    "c_1",
    "c_2",
    "sigma_k",
    "sigma_epsilon",
    "sigma_phi",
    "a",
    "c_t",
    "c_l",
    "c_eta",
    "kappa",
];

/// Relative width of the smooth switches with `smooth_limits = true`.
pub const SMOOTH_WIDTH: f64 = 1e-2;

impl TurbulenceSection {
    /// The default table (laminar flow; left out of the written case).
    pub fn is_default(&self) -> bool {
        *self == Self::default()
    }

    /// No turbulence model.
    pub fn is_laminar(&self) -> bool {
        self.model == TurbulenceModel::Laminar
    }

    /// The model constants with the overrides.
    pub fn sst_constants(&self) -> Result<SstConstants> {
        let mut c = SstConstants::default();
        for (name, &v) in &self.constants {
            let slot = match name.as_str() {
                "beta_star" => &mut c.beta_star,
                "a1" => &mut c.a1,
                "kappa" => &mut c.kappa,
                "gamma1" => &mut c.gamma1,
                "gamma2" => &mut c.gamma2,
                "beta1" => &mut c.beta1,
                "beta2" => &mut c.beta2,
                "sigma_k1" => &mut c.sigma_k1,
                "sigma_k2" => &mut c.sigma_k2,
                "sigma_omega1" => &mut c.sigma_w1,
                "sigma_omega2" => &mut c.sigma_w2,
                "c_r1" => &mut c.c_r1,
                "c_r2" => &mut c.c_r2,
                "c_r3" => &mut c.c_r3,
                "c_cr1" => &mut c.c_cr1,
                "c_cr2" => &mut c.c_cr2,
                _ => {
                    let hint = crate::bc::did_you_mean(name, SST_CONSTANTS.iter().copied())
                        .map(|h| format!(" (did you mean `{h}`?)"))
                        .unwrap_or_default();
                    return Err(diag::err(
                        Diagnostic::new(
                            diag::CASE_UNKNOWN_KEY,
                            format!(
                                "[turbulence.constants] has no constant `{name}`{hint}; the SST \
                                 constants are {}",
                                SST_CONSTANTS.join(", ")
                            ),
                        )
                        .key(format!("turbulence.constants.{name}")),
                    ));
                }
            };
            crate::need!(
                v > 0.0 && v.is_finite(),
                CASE_RANGE,
                "[turbulence.constants] {name} must be > 0, is {v}"
            );
            *slot = v;
        }
        Ok(c)
    }

    /// The constants of the k-ε family with the overrides.
    pub fn keps_constants(&self) -> Result<KepsConstants> {
        let mut c = if self.model == TurbulenceModel::RealizableKEpsilon {
            KepsConstants::realizable()
        } else {
            KepsConstants::standard()
        };
        for (name, &v) in &self.constants {
            let slot = match name.as_str() {
                "c_mu" => &mut c.c_mu,
                "c_epsilon1" => &mut c.c_e1,
                "c_epsilon2" => &mut c.c_e2,
                "sigma_k" => &mut c.sigma_k,
                "sigma_epsilon" => &mut c.sigma_e,
                "a0" => &mut c.a0,
                "kappa" => &mut c.kappa,
                _ => {
                    let hint = crate::bc::did_you_mean(name, KEPS_CONSTANTS.iter().copied())
                        .map(|h| format!(" (did you mean `{h}`?)"))
                        .unwrap_or_default();
                    return Err(diag::err(
                        Diagnostic::new(
                            diag::CASE_UNKNOWN_KEY,
                            format!(
                                "[turbulence.constants] has no constant `{name}`{hint}; the k-ε \
                                 constants are {}",
                                KEPS_CONSTANTS.join(", ")
                            ),
                        )
                        .key(format!("turbulence.constants.{name}")),
                    ));
                }
            };
            crate::need!(
                v > 0.0 && v.is_finite(),
                CASE_RANGE,
                "[turbulence.constants] {name} must be > 0, is {v}"
            );
            *slot = v;
        }
        Ok(c)
    }

    /// The constants of k-ε-ζ-f (or of the elliptic-blending k-ε) with the
    /// overrides.
    pub fn zeta_f_constants(&self) -> Result<ZetaFConstants> {
        let f0 = self.model == TurbulenceModel::KEpsilonZetaF0;
        let eb = self.model.is_eb();
        let mut c = if eb {
            ZetaFConstants::billard_laurence_2012()
        } else if f0 {
            ZetaFConstants::popovac_2023()
        } else {
            ZetaFConstants::paper()
        };
        let names = if eb { EB_CONSTANTS } else { ZETA_F_CONSTANTS };
        for (name, &v) in &self.constants {
            let slot = match name.as_str() {
                "c_mu" => &mut c.c_mu,
                "c_epsilon1" => &mut c.c_e1,
                "c_epsilon1_zeta" if !eb => &mut c.c_e1_zeta,
                "c_epsilon2" => &mut c.c_e2,
                "c_epsilon3" if eb => &mut c.c_e3,
                "c_epsilon4" if eb => &mut c.c_e4,
                "c_1" => &mut c.c1,
                "c_2" => &mut c.c2,
                "sigma_k" => &mut c.sigma_k,
                "sigma_epsilon" => &mut c.sigma_e,
                "sigma_zeta" if !eb => &mut c.sigma_z,
                "sigma_phi" if eb => &mut c.sigma_z,
                "a" => &mut c.a,
                "c_t" => &mut c.c_t,
                "c_l" => &mut c.c_l,
                "c_eta" => &mut c.c_eta,
                "kappa" => &mut c.kappa,
                "n_zeta" if f0 => &mut c.n_zeta,
                "c_epsilon1_power" if !eb => &mut c.c_e1_power,
                _ => {
                    let hint = crate::bc::did_you_mean(name, names.iter().copied())
                        .map(|h| format!(" (did you mean `{h}`?)"))
                        .unwrap_or_default();
                    return Err(diag::err(
                        Diagnostic::new(
                            diag::CASE_UNKNOWN_KEY,
                            format!(
                                "[turbulence.constants] has no constant `{name}`{hint}; the \
                                 constants of {} are {}",
                                self.model.zeta_label(),
                                names.join(", ")
                            ),
                        )
                        .key(format!("turbulence.constants.{name}")),
                    ));
                }
            };
            crate::need!(
                v > 0.0 && v.is_finite(),
                CASE_RANGE,
                "[turbulence.constants] {name} must be > 0, is {v}"
            );
            *slot = v;
        }
        if !eb {
            c.k_sink = self.k_sink.unwrap_or_default().code();
        }
        crate::need!(
            c.c_e1_power == 1.0 || c.c_e1_power == 0.5,
            CASE_RANGE,
            "[turbulence.constants] c_epsilon1_power must be 1 (C_ε1 = C_ε1′(1 + C_ε1″/ζ)) \
             or 0.5 (C_ε1″/√ζ), is {}",
            c.c_e1_power
        );
        Ok(c)
    }

    /// The run-time options of the k-ε-ζ-f kernel.
    pub fn zeta_f_options(&self, s: &TurbulenceSettings) -> ZetaFOptions {
        let (k_amb, e_amb) = self
            .ambient
            .map_or((0.0, 0.0), |a| (a.k, a.epsilon.unwrap_or(0.0)));
        ZetaFOptions {
            newton: u32::from(s.destruction == Some(Destruction::Newton)),
            smooth: if self.smooth_limits {
                SMOOTH_WIDTH
            } else {
                0.0
            },
            k_amb,
            e_amb,
            max_ratio: self.max_viscosity_ratio,
        }
    }

    /// The run-time options of the k-ε kernel.
    pub fn keps_options(&self, s: &TurbulenceSettings) -> KepsOptions {
        let (k_amb, e_amb) = self
            .ambient
            .map_or((0.0, 0.0), |a| (a.k, a.epsilon.unwrap_or(0.0)));
        let realizable = self.model == TurbulenceModel::RealizableKEpsilon;
        KepsOptions {
            variant: if realizable {
                kkeps::REALIZABLE
            } else {
                kkeps::STANDARD
            },
            limiter: match self.production_limiter.unwrap_or(if realizable {
                ProductionLimiter::None
            } else {
                ProductionLimiter::Durbin
            }) {
                ProductionLimiter::None => kkeps::LIMITER_NONE,
                ProductionLimiter::KatoLaunder => kkeps::LIMITER_KATO_LAUNDER,
                ProductionLimiter::Durbin => kkeps::LIMITER_DURBIN,
                ProductionLimiter::Menter => kkeps::LIMITER_CLIP,
            },
            newton: u32::from(s.destruction == Some(Destruction::Newton)),
            c_lim: self.durbin_coefficient,
            smooth: if self.smooth_limits {
                SMOOTH_WIDTH
            } else {
                0.0
            },
            k_amb,
            e_amb,
            max_ratio: self.max_viscosity_ratio,
        }
    }

    /// C_μ of the model (β* of SST; the equilibrium 0.09 of k-ε-ζ-f; the
    /// inlet values and the floors).
    pub fn c_mu(&self) -> Result<f64> {
        Ok(if self.model.has_zeta() {
            self.zeta_f_constants()?.c_mu_eq
        } else if self.model.is_epsilon() {
            self.keps_constants()?.c_mu
        } else {
            self.sst_constants()?.beta_star
        })
    }

    /// The run-time options of the kernel.
    pub fn sst_options(&self, s: &TurbulenceSettings) -> SstOptions {
        let (k_amb, w_amb) = self
            .ambient
            .map_or((0.0, 0.0), |a| (a.k, a.omega.unwrap_or(0.0)));
        SstOptions {
            limiter: match self.production_limiter.unwrap_or(ProductionLimiter::Menter) {
                ProductionLimiter::Menter => ksst::LIMITER_MENTER,
                ProductionLimiter::KatoLaunder => ksst::LIMITER_KATO_LAUNDER,
                ProductionLimiter::Durbin => ksst::LIMITER_DURBIN,
                ProductionLimiter::None => ksst::LIMITER_NONE,
            },
            newton: u32::from(s.destruction.unwrap_or_default() == Destruction::Newton),
            c_lim: self.durbin_coefficient,
            smooth: if self.smooth_limits {
                SMOOTH_WIDTH
            } else {
                0.0
            },
            k_amb,
            w_amb,
            max_ratio: self.max_viscosity_ratio,
        }
    }

    /// The wall treatment (the model's default when not set).
    pub fn treatment(&self) -> WallTreatment {
        self.wall_treatment
            .unwrap_or(if self.model == TurbulenceModel::KEpsilonZetaF0 {
                WallTreatment::Compound2023
            } else if self.model.has_zeta() {
                WallTreatment::Compound
            } else if self.model.is_epsilon() {
                WallTreatment::Scalable
            } else {
                WallTreatment::Automatic
            })
    }

    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        if self.is_laminar() {
            return Ok(());
        }
        let t = self.treatment();
        crate::need!(
            !self.curvature_correction || !(self.model.has_zeta() || self.model.is_epsilon()),
            CASE_INCONSISTENT,
            "turbulence.curvature_correction is an option of k-ω SST (Smirnov & Menter 2009)"
        );
        crate::need!(
            self.qcr == QcrMode::None || !(self.model.has_zeta() || self.model.is_epsilon()),
            CASE_INCONSISTENT,
            "turbulence.qcr is an option of k-ω SST (Spalart 2000; Mani et al. 2013)"
        );
        crate::need!(
            self.k_sink.is_none() || (self.model.has_zeta() && !self.model.is_eb()),
            CASE_INCONSISTENT,
            "turbulence.k_sink is an option of k-ε-ζ-f and k-ε-ζ-f₀"
        );
        if self.model.is_eb() {
            self.zeta_f_constants()?;
            crate::need!(
                matches!(t, WallTreatment::Compound | WallTreatment::LowRe),
                CASE_INCONSISTENT,
                "turbulence.wall_treatment = \"{}\" is not a treatment of the \
                 elliptic-blending k-ε; use \"compound\" (any y⁺) or \"low-re\" (y⁺ ≲ 1)",
                t.name()
            );
        }
        if self.model.has_zeta() {
            self.zeta_f_constants()?;
            crate::need!(
                matches!(
                    t,
                    WallTreatment::Compound | WallTreatment::Compound2023 | WallTreatment::LowRe
                ),
                CASE_INCONSISTENT,
                "turbulence.wall_treatment = \"{}\" is not a treatment of k-ε-ζ-f; use \
                 \"compound-2023\" or \"compound\" (any y⁺) or \"low-re\" (y⁺ ≲ 1)",
                t.name()
            );
            crate::need!(
                matches!(
                    self.production_limiter,
                    None | Some(ProductionLimiter::None)
                ),
                CASE_INCONSISTENT,
                "turbulence.production_limiter: {} bounds its time scale by \
                 realizability itself and takes no production limiter (use \"none\")",
                self.model.zeta_label()
            );
            if let Some(a) = self.ambient {
                crate::need!(
                    a.epsilon.is_some() && a.omega.is_none(),
                    CASE_INCONSISTENT,
                    "turbulence.ambient of {} takes k and epsilon",
                    self.model.zeta_label()
                );
            }
        } else if self.model.is_epsilon() {
            self.keps_constants()?;
            crate::need!(
                matches!(
                    t,
                    WallTreatment::Scalable | WallTreatment::Standard | WallTreatment::Compound
                ),
                CASE_INCONSISTENT,
                "turbulence.wall_treatment = \"{}\" is not a treatment of the k-ε family; use \
                 \"scalable\", \"standard\" or \"compound\"",
                t.name()
            );
            crate::need!(
                self.model == TurbulenceModel::KEpsilon
                    || matches!(
                        self.production_limiter,
                        None | Some(ProductionLimiter::None)
                    ),
                CASE_INCONSISTENT,
                "turbulence.production_limiter: the realizable k-ε model has its own \
                 realizable C_μ and takes no production limiter (use \"none\")"
            );
            if let Some(a) = self.ambient {
                crate::need!(
                    a.epsilon.is_some() && a.omega.is_none(),
                    CASE_INCONSISTENT,
                    "turbulence.ambient of a k-ε model takes k and epsilon"
                );
            }
        } else {
            self.sst_constants()?;
            crate::need!(
                matches!(t, WallTreatment::Automatic | WallTreatment::LowRe),
                CASE_INCONSISTENT,
                "turbulence.wall_treatment = \"{}\" is not a treatment of SST; use \
                 \"automatic\" or \"low-re\"",
                t.name()
            );
            if let Some(a) = self.ambient {
                crate::need!(
                    a.omega.is_some() && a.epsilon.is_none(),
                    CASE_INCONSISTENT,
                    "turbulence.ambient of SST takes k and omega"
                );
            }
        }
        crate::need!(
            self.durbin_coefficient > 0.0 && self.durbin_coefficient.is_finite(),
            CASE_RANGE,
            "turbulence.durbin_coefficient must be > 0, is {}",
            self.durbin_coefficient
        );
        crate::need!(
            self.max_viscosity_ratio > 0.0,
            CASE_RANGE,
            "turbulence.max_viscosity_ratio must be > 0, is {}",
            self.max_viscosity_ratio
        );
        if let Some(a) = self.ambient {
            let second = a.omega.or(a.epsilon).unwrap_or(0.0);
            crate::need!(
                a.k > 0.0 && second > 0.0 && a.k.is_finite() && second.is_finite(),
                CASE_RANGE,
                "turbulence.ambient: k and omega (or epsilon) must be > 0"
            );
        }
        Ok(())
    }
}

/// `[solver.turbulence]`: the transport equations of the turbulence model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct TurbulenceSettings {
    /// Convection scheme (`central` only in transient runs, TURB006).
    pub scheme: Scheme,
    /// Limiter of the gradients (`venkatakrishnan` and
    /// `venkatakrishnan-wang` are the same for k and ω: the threshold is
    /// scaled with the variable's range).
    pub limiter: LimiterKind,
    /// Relative threshold K of the limiter, ε = K (φ_max − φ_min) over the
    /// domain (Wang 2000). The default 0 (taken as 10⁻¹²) is the bounded
    /// Venkatakrishnan function: the dimensional threshold (K h)³ of the
    /// momentum equations is far above the differences of k (10⁻⁸ to
    /// 10⁻³ m²/s² on a flat plate) and would switch the limiter off.
    pub limiter_k: f64,
    /// Under-relaxation of the transport equations (default 0.8 steady, 1
    /// transient; with the pseudo-time step of the coupled algorithm a cap
    /// on the local step: 0.98, or 0.995 for the k-ε family with the
    /// reduction multigrid, D180).
    pub relaxation: Option<f64>,
    /// Relative tolerance per solve (default 1e-2 steady, 1e-8 transient).
    pub tolerance: Option<f64>,
    /// Largest number of BiCGStab iterations per solve.
    pub max_iter: usize,
    /// Linearisation of the destruction of ω or ε (default `newton` for
    /// ω, `picard` for ε: the Newton form takes at most half of ε away per
    /// iteration and let k collapse in decaying free streams).
    pub destruction: Option<Destruction>,
    /// Relaxation α_νt of ν_t.
    pub viscosity_relaxation: f64,
    /// Preconditioner of BiCGStab: `jacobi`; `reduction`, the reduction
    /// multigrid of the momentum sweeps (D179, D180) on the equation's own
    /// matrix, kept while the matrix changes by at most `reuse` (settings
    /// `[solver.coupled.multigrid]`); or `amg`, the aggregation multigrid
    /// of the pressure (D36) built from the symmetric part of the matrix
    /// (settings `[solver.pressure.amg]`; not with a pseudo-time step,
    /// where the symmetric part is not positive definite). Default:
    /// `reduction` for the k-ε family with the coupled algorithm and its
    /// pseudo-time step, `jacobi` otherwise.
    pub preconditioner: Option<TransportPreconditioner>,
    /// `preconditioner = "reduction"`: a hierarchy is kept while the
    /// matrix of its equation differs from the one it was built from by at
    /// most this fraction (Σ|Δa| / Σ|a|); 0 builds one for every solve.
    pub reuse: f64,
}

/// Preconditioner of the transport equations of a turbulence model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum TransportPreconditioner {
    /// Diagonal (Jacobi).
    Jacobi,
    /// One V-cycle of the smoothed-aggregation multigrid on the symmetric
    /// part of the matrix.
    Amg,
    /// One V-cycle of the reduction multigrid (D179) on the matrix.
    Reduction,
}

/// `[solver.elliptic]`: the elliptic equation of k-ε-ζ-f (f), stage 6
/// specification 6.3. The matrix of −∇·(L²∇f) + f is symmetric positive
/// definite: CG with the aggregation multigrid (its own hierarchy, rebuilt
/// for every solve, since L² changes with k and ε in every outer
/// iteration; a restart then continues bitwise without a stored matrix).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct EllipticSettings {
    /// Relative tolerance per solve (default 1e-2 steady, 1e-8 transient).
    pub tolerance: Option<f64>,
    /// Largest number of CG iterations per solve.
    pub max_iter: usize,
    /// Preconditioner of CG: `amg` (default) or `jacobi`.
    pub preconditioner: Preconditioner,
}

impl Default for EllipticSettings {
    fn default() -> Self {
        Self {
            tolerance: None,
            max_iter: 100,
            preconditioner: Preconditioner::Amg,
        }
    }
}

impl EllipticSettings {
    /// Relative tolerance of a solve.
    pub fn tol(&self, transient: bool) -> f64 {
        self.tolerance
            .unwrap_or(if transient { 1e-8 } else { 1e-2 })
    }

    /// Checks the ranges.
    pub fn validate(&self, transient: bool) -> Result<()> {
        crate::need!(
            self.tol(transient) > 0.0,
            CASE_RANGE,
            "[solver.elliptic] tolerance must be > 0"
        );
        crate::need!(
            self.max_iter >= 1,
            CASE_RANGE,
            "[solver.elliptic] max_iter must be >= 1"
        );
        Ok(())
    }
}

impl Default for TurbulenceSettings {
    fn default() -> Self {
        Self {
            scheme: Scheme::LinearUpwind,
            limiter: LimiterKind::Venkatakrishnan,
            limiter_k: 0.0,
            relaxation: None,
            tolerance: None,
            max_iter: 50,
            destruction: None,
            viscosity_relaxation: 1.0,
            preconditioner: None,
            reuse: 0.02,
        }
    }
}

impl TurbulenceSettings {
    /// α of k and ω.
    pub fn alpha(&self, transient: bool) -> f64 {
        self.relaxation.unwrap_or(if transient { 1.0 } else { 0.8 })
    }

    /// Relative tolerance of a solve.
    pub fn tol(&self, transient: bool) -> f64 {
        self.tolerance
            .unwrap_or(if transient { 1e-8 } else { 1e-2 })
    }

    /// Checks the ranges.
    pub fn validate(&self, transient: bool) -> Result<()> {
        let a = self.alpha(transient);
        crate::need!(
            a > 0.0 && a <= 1.0,
            CASE_RANGE,
            "[solver.turbulence] relaxation must be in (0, 1], is {a}"
        );
        crate::need!(
            self.limiter != LimiterKind::Mlp,
            CASE_RANGE,
            "[solver.turbulence] limiter = \"mlp\" is not built for the turbulence variables yet \
             (scalars and energy only)"
        );
        crate::need!(
            self.viscosity_relaxation > 0.0 && self.viscosity_relaxation <= 1.0,
            CASE_RANGE,
            "[solver.turbulence] viscosity_relaxation must be in (0, 1], is {}",
            self.viscosity_relaxation
        );
        crate::need!(
            self.tol(transient) > 0.0,
            CASE_RANGE,
            "[solver.turbulence] tolerance must be > 0"
        );
        crate::need!(
            self.max_iter >= 1,
            CASE_RANGE,
            "[solver.turbulence] max_iter must be >= 1"
        );
        crate::need!(
            self.limiter_k >= 0.0 && self.limiter_k.is_finite(),
            CASE_RANGE,
            "[solver.turbulence] limiter_k must be >= 0"
        );
        crate::need!(
            (0.0..=0.2).contains(&self.reuse),
            CASE_RANGE,
            "[solver.turbulence] reuse must lie in [0, 0.2], is {}",
            self.reuse
        );
        crate::need!(
            transient || self.scheme != Scheme::Central,
            TURB_UNBOUNDED_SCHEME,
            "[solver.turbulence] scheme = \"central\" is not bounded for k and ω in a steady run; \
             use \"linear-upwind\""
        );
        Ok(())
    }
}

/// The start values of the turbulence (`[initial]`, specification 7.3).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InitialTurbulence {
    /// Intensity.
    pub turbulence_intensity: Option<f64>,
    /// ν_t/ν.
    pub viscosity_ratio: Option<f64>,
    /// k directly (number or expression of x, y, z).
    pub k: Option<Scalar>,
    /// ω directly.
    pub omega: Option<Scalar>,
}

/// Checks the turbulence keys of a case.
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let t = &case.turbulence;
    t.validate()?;
    if t.is_laminar() {
        let i = &case.initial;
        crate::need!(
            i.turbulence_intensity.is_none()
                && i.viscosity_ratio.is_none()
                && i.k.is_none()
                && i.omega.is_none()
                && i.epsilon.is_none(),
            CASE_INCONSISTENT,
            "[initial] turbulence values in a laminar case (turbulence.model = \"laminar\")"
        );
        return Ok(());
    }
    let transient = case.time.is_some();
    case.solver.turbulence.validate(transient)?;
    case.solver.elliptic.validate(transient)?;
    if t.model.has_zeta() {
        for (name, bc) in &case.boundary {
            crate::need!(
                !(t.treatment() == WallTreatment::LowRe && bc.roughness() > 0.0),
                CASE_INCONSISTENT,
                "[boundary.{name}] roughness with {} needs the compound wall treatment \
                 (its log-law shift); the integration to the wall has no rough-wall form",
                t.model.zeta_label()
            );
        }
    }
    let i = &case.initial;
    crate::need!(
        !(i.omega.is_some() && i.epsilon.is_some()),
        CASE_INCONSISTENT,
        "[initial] omega and epsilon together: give one (with k)"
    );
    crate::need!(
        i.k.is_some() == (i.omega.is_some() || i.epsilon.is_some()),
        CASE_INCONSISTENT,
        "[initial] k goes with omega or epsilon"
    );
    crate::need!(
        i.k.is_none() || (i.turbulence_intensity.is_none() && i.viscosity_ratio.is_none()),
        CASE_INCONSISTENT,
        "[initial]: k with omega or epsilon, or turbulence_intensity and viscosity_ratio, \
         not both"
    );
    if let Some(z) = i.zeta {
        crate::need!(
            z > 0.0 && z <= 2.0,
            CASE_RANGE,
            "[initial] zeta = {z} must lie in (0, 2] (ζ = v²/k ≤ 2)"
        );
    }
    let env = case.env()?;
    for (key, v) in [
        ("initial.k", &i.k),
        ("initial.omega", &i.omega),
        ("initial.epsilon", &i.epsilon),
    ] {
        if let Some(v) = v {
            expr::value(v, &env, Vars::SPACE, key)?;
        }
    }
    for (name, bc) in &case.boundary {
        for (error, code, msg) in inlet::check(bc) {
            if error {
                return Err(diag::err(
                    Diagnostic::new(code, format!("[boundary.{name}]: {msg}"))
                        .key(format!("boundary.{name}")),
                ));
            }
        }
        if inlet::is_inlet(bc) || matches!(bc, BoundaryCondition::Wall { .. }) {
            InletSpec::of(
                &bc.turbulence(),
                &env,
                transient,
                &format!("boundary.{name}"),
            )?;
        }
        if let BoundaryCondition::Wall { .. } = bc {
            let ks = bc.roughness();
            crate::need!(
                ks >= 0.0 && ks.is_finite(),
                CASE_RANGE,
                "[boundary.{name}] roughness must be >= 0, is {ks}"
            );
        }
    }
    Ok(())
}

/// The run-time findings of the setup (notes and warnings for the log).
pub fn setup_findings(case: &crate::case::FlowCase) -> Vec<crate::bc::Finding> {
    use crate::bc::{Finding, Level};
    let mut out = Vec::new();
    let t = &case.turbulence;
    if t.is_laminar() {
        return out;
    }
    let mut push = |level, code, message: String| {
        out.push(Finding {
            level,
            code,
            boundary: None,
            message,
        })
    };
    if t.curvature_correction && t.production_limiter == Some(ProductionLimiter::KatoLaunder) {
        push(
            Level::Note,
            diag::TURB_CURVATURE_KATO,
            "turbulence.curvature_correction with production_limiter = \"kato-launder\": both act on \
             stagnation and swirl regions; Smirnov & Menter (2009) use the standard production"
                .to_string(),
        );
    }
    if !t.constants.is_empty() {
        push(
            Level::Warning,
            diag::TURB_CONSTANTS,
            format!(
                "{} constants changed in [turbulence.constants]: {}",
                if t.model.has_zeta() {
                    t.model.zeta_label()
                } else if t.model.is_epsilon() {
                    "k-ε"
                } else {
                    "SST"
                },
                t.constants
                    .iter()
                    .map(|(k, v)| format!("{k} = {v}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        );
    }
    if !case
        .boundary
        .values()
        .any(|b| matches!(b, BoundaryCondition::Wall { .. }))
    {
        push(
            Level::Note,
            diag::TURB_NO_WALL,
            if t.model.is_eb() {
                "no wall boundary: no wall treatment acts (α = 1, φ relaxes to 2/3)".into()
            } else if t.model.has_zeta() {
                "no wall boundary: no wall treatment acts (ζ-f relaxes ζ to 2/3)".into()
            } else if t.model.is_epsilon() {
                "no wall boundary: no wall treatment acts".into()
            } else {
                "no wall boundary: the wall distance is infinite (F₁ = 0, F₂ = 0)".into()
            },
        );
    }
    if case
        .time
        .as_ref()
        .is_some_and(|t| t.scheme == crate::transient::TimeScheme::CrankNicolson)
    {
        push(
            Level::Note,
            diag::TURB_UNBOUNDED_SCHEME,
            "Crank–Nicolson is not bounded for the turbulence variables: they take BDF2".into(),
        );
    }
    out
}

/// Results of the turbulence solves of one outer iteration.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct TurbResiduals {
    /// Scaled residuals per equation of [`equations`] (before the solves).
    pub residual: Vec<f64>,
    /// Linear iterations per equation (BiCGStab; CG for f).
    pub iterations: Vec<usize>,
    /// Cells where the positivity floor acted, per equation.
    pub floors: Vec<usize>,
    /// Cells where the ν_t cap acted.
    pub capped: usize,
}

impl TurbResiduals {
    /// No turbulence model.
    pub fn is_empty(&self) -> bool {
        self.residual.is_empty()
    }
}

/// The names of the equations, in the order of the solves and of the
/// residuals.
pub fn equations(model: TurbulenceModel) -> &'static [&'static str] {
    if model.is_eb() {
        &["k", "epsilon", "alpha", "phi"]
    } else if model.has_zeta() {
        &["k", "epsilon", "f", "zeta"]
    } else if model.is_epsilon() {
        &["k", "epsilon"]
    } else {
        &["k", "omega"]
    }
}

/// Host copy of the turbulence state (restart, output).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TurbState {
    /// The second variable is ε (the k-ε family): `omega`, `omega_old`
    /// and `omega_b` hold ε.
    pub epsilon: bool,
    /// k per cell.
    pub k: Vec<f64>,
    /// ω (SST) or ε (k-ε family, see `epsilon`) per cell.
    pub omega: Vec<f64>,
    /// ν_t per cell.
    pub nu_t: Vec<f64>,
    /// Transient: k and ω of the level before.
    pub k_old: Option<Vec<f64>>,
    /// See `k_old`.
    pub omega_old: Option<Vec<f64>>,
    /// The face viscosity μ_eff of the momentum equation.
    pub mu_eff: Vec<f64>,
    /// The boundary values of k and ω per boundary face (inlets and rough
    /// walls).
    pub k_b: Vec<f64>,
    /// See `k_b`.
    pub omega_b: Vec<f64>,
    /// u_τ per boundary face (0 off the walls).
    pub u_tau: Vec<f64>,
    /// First-cell y⁺ per boundary face (0 off the walls).
    pub y_plus: Vec<f64>,
    /// The viscous weight of the wall treatment per boundary face (0 off
    /// the walls): 1 for the integration to the wall, the weight of the
    /// automatic or compound blending, 0 in the log layer (spec 2.3, the
    /// wall field `blending`). `None` from a restart file written before
    /// it was stored.
    pub wall_blending: Option<Vec<f64>>,
    /// The modified pressure p̂ = p + ⅔ρk of the solver, per cell (the
    /// state's `p` is the static pressure).
    pub p_modified: Vec<f64>,
    /// The modified pressure per boundary face.
    pub pb_modified: Vec<f64>,
    /// The wall distance per cell.
    pub wall_distance: Vec<f64>,
    /// k-ε-ζ-f: ζ and f per cell, ζ of the level before (transient), and
    /// their boundary values per boundary face.
    pub zeta: Option<ZetaState>,
    /// The matrices (values in the order of the matrix pattern) of the
    /// last multigrid setups of the turbulence equations, per hierarchy
    /// slot (1 f, 2 k, 3 ε or ω, 4 ζ): a continued run rebuilds the same
    /// hierarchies from them (bit-identical restarts, D42). Empty without
    /// a multigrid of the turbulence equations; dropped when a restart
    /// file is mapped to another cell order.
    pub amg: Vec<(usize, Vec<f64>)>,
    /// f_r1 of the curvature correction per cell (output only; not in the
    /// restart file, recomputed in every iteration).
    pub fr1: Option<Vec<f64>>,
}

/// The ζ-f part of a [`TurbState`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ZetaState {
    /// ζ per cell.
    pub zeta: Vec<f64>,
    /// f per cell.
    pub f: Vec<f64>,
    /// ζ of the level before (transient).
    pub zeta_old: Option<Vec<f64>>,
    /// ζ per boundary face (inlet values).
    pub zeta_b: Vec<f64>,
    /// f per boundary face (wall values).
    pub f_b: Vec<f64>,
    /// The time scale T per cell (output).
    pub t: Vec<f64>,
}

/// Time terms of a step.
pub type StepTime = crate::scalar::StepTime;

/// One transported turbulence variable on the device.
#[derive(Debug, Clone, Copy)]
struct Var {
    phi: Buf,
    grad: Buf3,
    psi: Buf,
    bc_a: Buf,
    bc_b: Buf,
    /// Inlet values per boundary face (0 elsewhere).
    b_in: Buf,
    /// Wall values per boundary face (rough walls; 0 elsewhere).
    b_wall: Buf,
    lsq: Buf3,
    old: Option<(Buf, Buf, Buf)>,
    /// The sliding mesh (D182): the old levels n and n−1 of the time term
    /// with the data of the rotor cells moved along the turn.
    content: Option<(Buf, Buf)>,
    valid: Option<Buf>,
    floor: f64,
}

/// The wall treatment on the device.
#[derive(Debug, Clone, Copy)]
struct WallDev {
    normal: Buf3,
    y: Buf,
    ks: Buf,
    cs: Buf,
    prod_w: Buf,
    u_tau: Buf,
    y_plus: Buf,
    /// The y⁺ of the thermal law (y* = C_μ^¼ k_P^½ y_P/ν of the k-ε
    /// family, y⁺ of SST).
    y_star: Buf,
    omega: Buf,
    blending: Buf,
    /// The viscous weight of the wall treatment (output only).
    viscous_weight: Buf,
    pk: Buf,
    weight: Buf,
}

/// The k-ε-ζ-f part on the device: ζ, f and their coefficients.
#[derive(Debug, Clone, Copy)]
struct ZetaFDev {
    c: ZetaFConstants,
    o: ZetaFOptions,
    z: Var,
    f: Var,
    /// T and L² per cell.
    t: Buf,
    l2: Buf,
    /// ν_t/σ_ζ per cell.
    dz: Buf,
    /// ρ(ν + ν_t/σ_ζ) and ρL² per face.
    gamma_z: Buf,
    gamma_f: Buf,
    /// S_u, S_p of ζ and S_u (the right-hand side), S_p of f.
    src: [Buf; 4],
    /// The values of the transported variable before its solve (the
    /// relative lower bound, [`ZETA_F_DROP`]).
    prev: Buf,
}

/// The extra inputs of the elliptic-blending k-ε on the device: the
/// second derivatives of the velocity (the E term) and the turbulent
/// transport of k (C*_ε2), from gradients of the gradients with zero
/// normal gradient on every boundary (∂²U/∂y² = ∂p/∂x/μ ≈ 0 at a wall).
#[derive(Debug, Clone, Copy)]
struct EbDev {
    /// Σ_ijl (∂_l ∂_j U_i)² per cell.
    g2: Buf,
    /// ∇·(ν_t/σ_k ∇k) per cell.
    dkt: Buf,
    /// A gradient and the flux ν_t/σ_k ∇k per cell.
    tmp: Buf3,
    q: Buf3,
    /// The least-squares weights with every boundary point at the owner
    /// value, and the boundary coefficients a = 1, b = 0.
    lsq: Buf3,
    ones: Buf,
    zeros: Buf,
}

/// The device data of the QCR of SST (D170) with the consistent stress
/// divergence of D236: the kernel mode, the cell stress, its cell
/// gradients (zero normal gradient on the boundary) for the skewness
/// correction of the face stress, and on the cells where the
/// face-neighbour velocity gradient is not second order the quadratic
/// least-squares gradient (weights, coordinates) into `gq`.
#[derive(Debug, Clone, Copy)]
struct QcrDev {
    mode: u32,
    tau: [Buf; 6],
    gtau: [Buf3; 6],
    lsq: Buf3,
    ones: Buf,
    zeros: Buf,
    /// The quadratic fit and the velocity gradient it overwrites.
    quad: Option<(QuadDev, [Buf3; 3])>,
}

/// A quadratic least-squares fit of D236 on the device: the stencil slot,
/// the gradient weights and the coordinates of the fit.
#[derive(Debug, Clone, Copy)]
struct QuadDev {
    slot: usize,
    w: Buf3,
    dims: u32,
}

/// The device data of the curvature correction of SST (D170): the six
/// components of S_ij, their values at the last time level, DS_ij/Dt and
/// f_r1, with the zero-gradient least-squares weights of [`EbDev`].
#[derive(Debug, Clone, Copy)]
struct CurvDev {
    s: [Buf; 6],
    s_old: [Buf; 6],
    dsdt: [Buf; 6],
    fr1: Buf,
    tmp: Buf3,
    lsq: Buf3,
    ones: Buf,
    zeros: Buf,
    /// Ω^rot and the frame velocity Ω × (x − x₀) per cell (moving
    /// reference frames, D55), if any.
    frame: Option<(Buf3, Buf3)>,
    /// Ω^rot of every cell and the spin of a relative velocity (the
    /// rotating frame of the whole domain, `flow.frame_rotation`).
    rot: [f64; 3],
    spin: [f64; 3],
    /// `s_old` holds a time level.
    has_old: bool,
    /// D236: the quadratic velocity fit on the cells where the
    /// face-neighbour gradient is first order, with the gradient it
    /// overwrites (S_ij from it is second order everywhere) …
    quad: Option<(QuadDev, [Buf3; 3])>,
    /// … and the quadratic fit of S_ij over interior cells on the cells
    /// where the gradient of the gradient is not exact, with the zero
    /// boundary values it never reads and the gradients of three
    /// components it overwrites.
    sfit: Option<(QuadDev, Buf3, [Buf3; 3])>,
}

/// The curvature correction computed S_ij since the start or the last
/// restart (so that a time level of S_ij can be kept).
#[derive(Debug, Default)]
struct CurvFresh(std::cell::Cell<bool>);

impl Clone for CurvFresh {
    fn clone(&self) -> Self {
        Self(std::cell::Cell::new(self.0.get()))
    }
}

/// The pseudo-time step of the turbulence transport equations in a steady
/// iteration (`[solver.pseudo_time]`): Δτ = min(Δτ_max, c_T T) per cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PseudoStep {
    /// 1/Δτ_max.
    pub inv_dt: f64,
    /// 1/c_T.
    pub rate_factor: f64,
}

/// Default under-relaxation of the transport equations with a pseudo-time
/// step: the step takes the place of the relaxation, which only caps the
/// local Courant number at α/(1 − α) = 49, as for the momentum equations.
/// (Measured on the backward-facing step: without the cap k-ε diverges in
/// the first iterations from a uniform start; with 0.9 the turbulence lags
/// the momentum equations and the run needs several times the iterations.)
pub const PSEUDO_RELAXATION: f64 = 0.98;

/// The same default for the k-ε family with the reduction multigrid as
/// preconditioner (D180): a local Courant number of at most 199. The
/// relaxation is a local step that shrinks with the cells (α/(1 − α) h/U,
/// or about half of that times h²/ν_eff in flat cells at a wall), so it
/// makes the number of outer iterations grow with the mesh: 150 on the
/// backward-facing step with 32 800 cells against 65 with 8 200 at 0.98.
/// Above 0.98 Jacobi-BiCGStab no longer solves the equations within its
/// 50 iterations; with the multigrid one or two iterations do. (The
/// steady state of the SST model on that step is not held with accurate
/// solves of its equations, at any relaxation: SST keeps 0.98 and the
/// Jacobi preconditioner, step 22 section 7.8.)
pub const PSEUDO_RELAXATION_REDUCTION: f64 = 0.995;

/// A kept reduction hierarchy of a transport equation is rebuilt for the
/// next solve when BiCGStab needed more than this many iterations with it
/// (a fresh one needs one or two).
const REDUCTION_STALE: usize = 8;

/// The case runs the coupled algorithm with pseudo-time steps in the
/// turbulence equations (steady).
fn pseudo_steps(case: &crate::case::FlowCase) -> bool {
    let alg = case.solver.algorithm(false);
    alg == crate::simplec::Algorithm::Coupled
        && case.solver.pseudo_time.turbulence_factor(alg) > 0.0
}

/// The k-ε family with the coupled algorithm and the reduction multigrid
/// in its transport equations (D180): the case where the defaults of the
/// relaxations are 0.995 instead of 0.98.
pub fn reduction_defaults(case: &crate::case::FlowCase) -> bool {
    case.time.is_none()
        && case.turbulence.model.is_epsilon()
        && pseudo_steps(case)
        && case
            .solver
            .turbulence
            .preconditioner
            .is_none_or(|p| p == TransportPreconditioner::Reduction)
}

/// k-ε-ζ-f: one outer iteration lowers k, ε or ζ of a cell by at most
/// this factor. An unconverged Krylov iterate of the loose steady
/// tolerance can be negative where the solution is small; the absolute
/// floor alone would then put ε near zero, and T = k/ε, L ∝ (ν³/ε)^¼ and
/// P/ε in the source of f would blow up (stage 6 specification 6.4,
/// 6d amendment). Inactive at convergence and in transient runs, whose
/// steps change the fields by small factors.
pub const ZETA_F_DROP: f64 = 0.1;

/// A backflow face of an outlet: boundary face index and, where the flow
/// enters, the entering k, ω or ε and ζ.
type BackflowFace = (usize, Option<(f64, f64, f64)>);

/// The inlet faces.
#[derive(Debug, Clone)]
struct Inlets {
    specs: Vec<InletSpec>,
    /// ζ per spec (k-ε-ζ-f; default 2/3; a table column wins).
    zeta: Vec<f64>,
    /// Per inlet face: boundary face index, spec, centroid.
    faces: Vec<(u32, usize, [f64; 3])>,
    speed: bool,
    time: bool,
}

/// The turbulence model of a flow on its device: k-ω SST or the k-ε
/// family (k and the second variable, ω or ε), or k-ε-ζ-f (k, ε and ζ
/// transported, f elliptic).
pub struct Turbulence {
    /// The model settings.
    pub section: TurbulenceSection,
    /// The solver settings.
    pub settings: TurbulenceSettings,
    model: TurbulenceModel,
    /// The second variable is ε (the k-ε family).
    eps: bool,
    c: SstConstants,
    o: SstOptions,
    ke: KepsConstants,
    ko: KepsOptions,
    /// C_μ (β* of SST): inlet values and floors.
    c_mu: f64,
    treatment: u32,
    nu: f64,
    rho: f64,
    /// ν and ρ per cell and per boundary face (its owner's) of a case of
    /// several fluids (D194); `None`: `nu` and `rho` everywhere.
    fluids: Option<CellFluids>,
    nc: usize,
    ni: usize,
    nb: usize,
    k: Var,
    w: Var,
    nu_t: Buf,
    nu_t_prev: Option<Buf>,
    /// F₁ of SST, C_μ of the k-ε family.
    f1: Buf,
    pk: Buf,
    dk: Buf,
    dw: Buf,
    capped: Buf,
    src: [Buf; 4],
    d: Buf,
    gu: [Buf3; 3],
    mu_lam: Buf,
    gamma_k: Buf,
    gamma_w: Buf,
    mask_b: Buf,
    walls: Walls,
    wd: Option<WallDev>,
    fix_mask: Buf,
    fix_val: Buf,
    wall_pk: Buf,
    prod_mask: Buf,
    zero_face: Buf,
    cp: Buf,
    cn: Buf,
    bf: Buf,
    val: Buf,
    rhs: Buf,
    diag: Buf,
    order: Buf,
    work: Buf,
    inlets: Inlets,
    /// k-ε-ζ-f: ζ and f.
    zf: Option<ZetaFDev>,
    /// The elliptic-blending k-ε (with `zf`).
    eb: Option<EbDev>,
    /// The curvature correction of SST.
    curv: Option<CurvDev>,
    curv_fresh: CurvFresh,
    /// The QCR of SST.
    qcr: Option<QcrDev>,
    /// Multigrid settings of the turbulence and elliptic solves (those of
    /// the pressure).
    amg: foehncraft_fv::amg::AmgSettings,
    /// The elliptic solves (f).
    elliptic: EllipticSettings,
    /// Rebuild the multigrid hierarchies of the turbulence equations in
    /// this outer iteration (set by the caller on the schedule of the
    /// pressure, `amg_rebuild`).
    pub(crate) amg_rebuild: bool,
    /// The pseudo-time step of the transport equations in the steady
    /// iteration that follows (set by the flow loop; None: relaxation
    /// alone).
    pub(crate) pseudo: Option<PseudoStep>,
    /// ρ V/Δτ per cell of that step.
    pseudo_coef: Buf,
    /// The matrix of the last setup per hierarchy slot (index = slot).
    amg_built: std::cell::RefCell<Vec<Option<Vec<f64>>>>,
    /// The preconditioner of the transport equations in force.
    transport_pc: TransportPreconditioner,
    /// Settings of the reduction multigrid (`[solver.coupled.multigrid]`).
    reduction: foehncraft_fv::amg::reduction::ReductionSettings,
    /// Setups of the reduction multigrid of the transport equations.
    reduction_setups: std::cell::Cell<usize>,
    /// Solves that took the Jacobi preconditioner after a failed setup.
    reduction_fallbacks: std::cell::Cell<usize>,
    /// The message of the first such failure, until the run loop reports
    /// it.
    note: std::cell::RefCell<Option<String>>,
    /// Host copies of a and of the inlet b of k (for p̂ at the pressure
    /// boundaries).
    ka_host: Vec<f64>,
    kin_host: Vec<f64>,
    /// Host copies of a and of the inlet b of ω or ε and of ζ (the
    /// backflow faces of the outlets change them).
    wa_host: Vec<f64>,
    win_host: Vec<f64>,
    za_host: Vec<f64>,
    zin_host: Vec<f64>,
    /// The faces of the outlets with backflow values (BC028): where the
    /// flow enters, k, ω or ε and ζ take the given values (fixed value),
    /// elsewhere zero gradient.
    backflow: Inlets,
    /// The state of the backflow faces of the last [`Self::set_backflow`]
    /// (boundary face, entering values or `None`), applied again after the
    /// inlet values are rebuilt.
    bf_state: Vec<BackflowFace>,
    /// The pressure faces (boundary face index) and their owners.
    pressure_faces: Vec<u32>,
    pressure_owners: Vec<u32>,
    added: Option<(Buf, Buf)>,
    /// The added sources of ζ and f (k-ε-ζ-f; φ and α of EB), if any.
    added_z: Option<(Buf, Buf)>,
    /// Zeros per cell.
    zero_cell: Buf,
    t_set: f64,
    /// The wall distance per cell (host).
    pub wall_distance: Vec<f64>,
    /// Health (D69): per equation (in the order of [`equations`]) the
    /// iterations in which the floor acted on more than 0.1 % of the cells,
    /// and the most cells floored in one iteration; the most cells with a
    /// capped ν_t.
    pub health: TurbHealth,
}

/// A wall's y⁺ statistics: boundary, min, mean, max, area fraction above 2.
pub type YPlusRow = (String, f64, f64, f64, f64);

/// Floor and cap counts of a run (D69, TURB002).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct TurbHealth {
    /// Iterations with floors on more than 0.1 % of the cells, per
    /// equation (k and ω or ε; k-ε-ζ-f: k, ε, f (never floored), ζ).
    pub floor_iterations: Vec<usize>,
    /// Of those, the ones in the final tenth of the solves so far (the
    /// inexact linear solves of the first iterations may undershoot; a
    /// converging run does not need the floor at its end).
    pub floor_iterations_late: Vec<usize>,
    /// Per solve, a bit per equation: the floor acted on more than 0.1 %.
    #[serde(skip)]
    flags: Vec<u8>,
    /// Most cells floored in one iteration, per equation.
    pub floor_cells_max: Vec<usize>,
    /// Most cells with ν_t capped by `max_viscosity_ratio` in one iteration.
    pub capped_max: usize,
}

/// A new buffer with the values `v`.
fn up(dev: &mut dyn FlowDevice, v: &[f64]) -> Result<Buf> {
    let b = dev.alloc(v.len())?;
    dev.upload(b, v)?;
    Ok(b)
}

/// A transported or elliptic turbulence variable (for its boundary kinds).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VarKind {
    /// k.
    K,
    /// ω or ε.
    Second,
    /// ζ (k-ε-ζ-f).
    Zeta,
    /// f (k-ε-ζ-f).
    F,
}

/// The boundary kinds of a turbulence variable per boundary (for the
/// gradient weights) and the coefficient a per boundary face. Inlets fix
/// k, ω or ε and ζ (f has zero gradient there). Walls: SST fixes ω, and k
/// with the low-Re treatment or on rough walls; the k-ε family has zero
/// gradient of k and ε (ε is fixed in the cell); k-ε-ζ-f fixes f always,
/// and k, ε and ζ with the integration to the wall (`low-re`; with the
/// compound treatment k and ζ have zero gradient, since dk/dy = dζ/dy = 0
/// at the wall, and ε is fixed in the cell).
fn kinds(
    mesh: &Mesh,
    bcs: &BTreeMap<String, BoundaryCondition>,
    treatment: WallTreatment,
    var: VarKind,
    model: TurbulenceModel,
) -> (BTreeMap<String, BoundaryKind>, Vec<f64>) {
    let ni = mesh.n_interior_faces();
    let mut a = vec![1.0; mesh.n_boundary_faces()];
    let mut points = BTreeMap::new();
    let low_re = treatment == WallTreatment::LowRe;
    for p in &mesh.boundaries {
        let bc = bcs.get(&p.name);
        let fixed = match bc {
            Some(b) if inlet::is_inlet(b) => var != VarKind::F,
            Some(b @ BoundaryCondition::Wall { .. }) => {
                if model.has_zeta() {
                    var == VarKind::F || low_re
                } else if model.is_epsilon() {
                    false
                } else {
                    var == VarKind::Second || b.roughness() > 0.0 || low_re
                }
            }
            _ => false,
        };
        let kind = if fixed {
            BoundaryKind::FixedValue
        } else {
            BoundaryKind::ZeroGradient
        };
        points.insert(p.name.clone(), kind);
        if fixed {
            for f in p.range() {
                a[f - ni] = 0.0;
            }
        }
    }
    (points, a)
}

/// ν and ρ of a case of several fluids (D194): per cell (host and device)
/// and per boundary face (the owner's).
#[derive(Debug, Clone)]
struct CellFluids {
    buf: Buf,
    nu: Vec<f64>,
    rho: Vec<f64>,
    nu_b: Vec<f64>,
    rho_b: Vec<f64>,
}

impl Turbulence {
    /// ν of cell `c`.
    fn nu_c(&self, c: usize) -> f64 {
        self.fluids.as_ref().map_or(self.nu, |f| f.nu[c])
    }

    /// ν at boundary face `fb` (its owner's).
    fn nu_b(&self, fb: usize) -> f64 {
        self.fluids.as_ref().map_or(self.nu, |f| f.nu_b[fb])
    }

    /// ρ at boundary face `fb` (its owner's).
    fn rho_b(&self, fb: usize) -> f64 {
        self.fluids.as_ref().map_or(self.rho, |f| f.rho_b[fb])
    }

    /// The turbulence model of a case on the device `dev` of a flow.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        dev: &mut dyn FlowDevice,
        case: &crate::case::FlowCase,
        mesh: &Mesh,
        geo: &Geometry,
        disc: &foehncraft_fv::disc::Discretisation,
        bc_arrays: &BoundaryArrays,
        rho: f64,
        mu: f64,
        face_viscosity: Buf,
        wall_distance: Option<Vec<f64>>,
        cells: Option<(Vec<f64>, Vec<f64>)>,
    ) -> Result<Self> {
        let section = case.turbulence.clone();
        let settings = case.solver.turbulence.clone();
        let transient = case.time.is_some();
        let env = case.env()?;
        let (nc, nf, ni) = (disc.n_cells, disc.n_faces, disc.n_interior);
        let nb = nf - ni;
        let stored = disc.sell.stored();
        let model = section.model;
        let eps = model.is_epsilon();
        let transport_pc =
            settings
                .preconditioner
                .unwrap_or(if eps && !transient && pseudo_steps(case) {
                    TransportPreconditioner::Reduction
                } else {
                    TransportPreconditioner::Jacobi
                });
        let c = if eps {
            SstConstants::default()
        } else {
            section.sst_constants()?
        };
        let o = section.sst_options(&settings);
        // the ζ-f constants are read only for the ζ-f family: the
        // overrides of SST (`gamma2`) or k-ε (`a0`) are not ζ-f names
        let zc = if model.has_zeta() {
            section.zeta_f_constants()?
        } else {
            ZetaFConstants::default()
        };
        // the wall kernel of the k-ε family serves k-ε-ζ-f with its
        // equilibrium C_μ (compound treatment)
        let ke = if model.has_zeta() {
            zc.wall()
        } else if eps {
            section.keps_constants()?
        } else {
            KepsConstants::default()
        };
        let ko = section.keps_options(&settings);
        let c_mu = section.c_mu()?;
        let treatment = section.treatment();
        let nu = mu / rho;
        let fluids = match cells {
            Some((nu_c, rho_c)) => {
                let owner_b = &disc.owner[ni..nf];
                Some(CellFluids {
                    buf: {
                        let b = dev.alloc(nc)?;
                        dev.upload(b, &nu_c)?;
                        b
                    },
                    nu_b: owner_b.iter().map(|&o| nu_c[o as usize]).collect(),
                    rho_b: owner_b.iter().map(|&o| rho_c[o as usize]).collect(),
                    nu: nu_c,
                    rho: rho_c,
                })
            }
            None => None,
        };
        // the inlets
        let mut specs = Vec::new();
        let mut zetas = Vec::new();
        let mut faces = Vec::new();
        let mut spec_of = BTreeMap::new();
        for p in &mesh.boundaries {
            let Some(bc) = case.boundary.get(&p.name) else {
                continue;
            };
            if !inlet::is_inlet(bc) {
                continue;
            }
            let s = InletSpec::of_inlet(bc, &env, transient, &format!("boundary.{}", p.name))?;
            spec_of.insert(p.name.clone(), specs.len());
            specs.push(s);
            zetas.push(bc.turbulence().zeta.unwrap_or(2.0 / 3.0));
            for f in p.range() {
                faces.push(((f - ni) as u32, specs.len() - 1, geo.face_centroid[f]));
            }
        }
        let inlets = Inlets {
            speed: specs.iter().any(InletSpec::speed_dependent),
            time: specs.iter().any(InletSpec::time_dependent),
            specs,
            zeta: zetas,
            faces,
        };
        // the outlets with backflow values
        let mut specs = Vec::new();
        let mut zetas = Vec::new();
        let mut faces = Vec::new();
        for p in &mesh.boundaries {
            let Some(bc) = case.boundary.get(&p.name) else {
                continue;
            };
            if !inlet::has_backflow_values(bc) {
                continue;
            }
            let s = InletSpec::of(
                &bc.turbulence(),
                &env,
                transient,
                &format!("boundary.{}", p.name),
            )?;
            specs.push(s);
            zetas.push(bc.turbulence().zeta.unwrap_or(2.0 / 3.0));
            for f in p.range() {
                faces.push(((f - ni) as u32, specs.len() - 1, geo.face_centroid[f]));
            }
        }
        let backflow = Inlets {
            speed: specs.iter().any(InletSpec::speed_dependent),
            time: specs.iter().any(InletSpec::time_dependent),
            specs,
            zeta: zetas,
            faces,
        };
        // the walls
        let walls = Walls::new(mesh, geo, &case.boundary, eps);
        let var = |dev: &mut dyn FlowDevice, vk: VarKind| -> Result<Var> {
            let (points, a) = kinds(mesh, &case.boundary, treatment, vk, model);
            let lsq = crate::simplec::disc_with(mesh, geo, &points)?.lsq;
            let lsq = Buf3 {
                x: up(dev, &lsq.x)?,
                y: up(dev, &lsq.y)?,
                z: up(dev, &lsq.z)?,
            };
            let bc_a = up(dev, &a)?;
            let psi = dev.alloc(nc)?;
            dev.fill(psi, 1.0)?;
            Ok(Var {
                phi: dev.alloc(nc)?,
                grad: Buf3 {
                    x: dev.alloc(nc)?,
                    y: dev.alloc(nc)?,
                    z: dev.alloc(nc)?,
                },
                psi,
                bc_a,
                bc_b: dev.alloc(nb)?,
                b_in: dev.alloc(nb)?,
                b_wall: dev.alloc(nb)?,
                lsq,
                old: if transient && vk != VarKind::F {
                    Some((dev.alloc(nc)?, dev.alloc(nc)?, dev.alloc(nc)?))
                } else {
                    None
                },
                content: None,
                valid: None,
                floor: 0.0,
            })
        };
        let k = var(dev, VarKind::K)?;
        let w = var(dev, VarKind::Second)?;
        let zf = if model.has_zeta() {
            Some(ZetaFDev {
                c: zc,
                o: section.zeta_f_options(&settings),
                z: var(dev, VarKind::Zeta)?,
                f: var(dev, VarKind::F)?,
                t: dev.alloc(nc)?,
                l2: dev.alloc(nc)?,
                dz: dev.alloc(nc)?,
                gamma_z: dev.alloc(nf)?,
                gamma_f: dev.alloc(nf)?,
                src: [
                    dev.alloc(nc)?,
                    dev.alloc(nc)?,
                    dev.alloc(nc)?,
                    dev.alloc(nc)?,
                ],
                prev: dev.alloc(nc)?,
            })
        } else {
            None
        };
        let eb = if model.is_eb() {
            let points = mesh
                .boundaries
                .iter()
                .map(|p| (p.name.clone(), BoundaryKind::ZeroGradient))
                .collect();
            let lsq = crate::simplec::disc_with(mesh, geo, &points)?.lsq;
            let ones = dev.alloc(nb)?;
            dev.fill(ones, 1.0)?;
            Some(EbDev {
                g2: dev.alloc(nc)?,
                dkt: dev.alloc(nc)?,
                tmp: Buf3 {
                    x: dev.alloc(nc)?,
                    y: dev.alloc(nc)?,
                    z: dev.alloc(nc)?,
                },
                q: Buf3 {
                    x: dev.alloc(nc)?,
                    y: dev.alloc(nc)?,
                    z: dev.alloc(nc)?,
                },
                lsq: Buf3 {
                    x: up(dev, &lsq.x)?,
                    y: up(dev, &lsq.y)?,
                    z: up(dev, &lsq.z)?,
                },
                ones,
                zeros: dev.alloc(nb)?,
            })
        } else {
            None
        };
        // the quadratic least-squares fits on the cells where the QCR or
        // the curvature correction needs them (D236)
        let (quad, sfit) = if section.qcr != QcrMode::None || section.curvature_correction {
            let (q, qs) = crate::operators::quadratic_fits(case, mesh, geo, &bc_arrays.u_points)?;
            let mut setup =
                |q: &foehncraft_fv::exactness::QuadraticGradient| -> Result<Option<QuadDev>> {
                    if q.cells.is_empty() {
                        return Ok(None);
                    }
                    let slot = dev.quadratic_setup(&q.cells, &q.ptr, &q.idx)?;
                    Ok(Some(QuadDev {
                        slot,
                        w: Buf3 {
                            x: up(dev, &q.w.x)?,
                            y: up(dev, &q.w.y)?,
                            z: up(dev, &q.w.z)?,
                        },
                        dims: q.dims as u32,
                    }))
                };
            let quad = setup(&q)?;
            let sfit = match qs {
                Some(qs) => setup(&qs)?,
                None => None,
            };
            (quad, sfit)
        } else {
            (None, None)
        };
        let a3n = |dev: &mut dyn FlowDevice, n: usize| -> Result<Buf3> {
            Ok(Buf3 {
                x: dev.alloc(n)?,
                y: dev.alloc(n)?,
                z: dev.alloc(n)?,
            })
        };
        let curv = if section.curvature_correction {
            let points = mesh
                .boundaries
                .iter()
                .map(|p| (p.name.clone(), BoundaryKind::ZeroGradient))
                .collect();
            let lsq = crate::simplec::disc_with(mesh, geo, &points)?.lsq;
            let ones = dev.alloc(nb)?;
            dev.fill(ones, 1.0)?;
            let six = |dev: &mut dyn FlowDevice| -> Result<[Buf; 6]> {
                let mut v = [ones; 6];
                for b in &mut v {
                    *b = dev.alloc(nc)?;
                }
                Ok(v)
            };
            let spin = case.flow.frame_rotation;
            Some(CurvDev {
                s: six(dev)?,
                s_old: six(dev)?,
                dsdt: six(dev)?,
                fr1: {
                    let b = dev.alloc(nc)?;
                    dev.fill(b, 1.0)?;
                    b
                },
                tmp: Buf3 {
                    x: dev.alloc(nc)?,
                    y: dev.alloc(nc)?,
                    z: dev.alloc(nc)?,
                },
                lsq: Buf3 {
                    x: up(dev, &lsq.x)?,
                    y: up(dev, &lsq.y)?,
                    z: up(dev, &lsq.z)?,
                },
                ones,
                zeros: dev.alloc(nb)?,
                frame: None,
                rot: spin,
                spin,
                has_old: false,
                quad: match quad {
                    Some(q) => Some((q, [a3n(dev, nc)?, a3n(dev, nc)?, a3n(dev, nc)?])),
                    None => None,
                },
                sfit: match sfit {
                    Some(q) => {
                        let zb = dev.alloc(nb)?;
                        let zeros = Buf3 {
                            x: zb,
                            y: zb,
                            z: zb,
                        };
                        Some((q, zeros, [a3n(dev, nc)?, a3n(dev, nc)?, a3n(dev, nc)?]))
                    }
                    None => None,
                },
            })
        } else {
            None
        };
        let qcr = if section.qcr == QcrMode::None {
            None
        } else {
            let mut t = [face_viscosity; 6];
            for b in &mut t {
                *b = dev.alloc(nc)?;
            }
            let a3 = |dev: &mut dyn FlowDevice| -> Result<Buf3> {
                Ok(Buf3 {
                    x: dev.alloc(nc)?,
                    y: dev.alloc(nc)?,
                    z: dev.alloc(nc)?,
                })
            };
            let mut gtau = [Buf3 {
                x: face_viscosity,
                y: face_viscosity,
                z: face_viscosity,
            }; 6];
            for g in &mut gtau {
                *g = a3(dev)?;
            }
            let points = mesh
                .boundaries
                .iter()
                .map(|p| (p.name.clone(), BoundaryKind::ZeroGradient))
                .collect();
            let lsq = crate::simplec::lsq_with(mesh, geo, &points)?;
            let lsq = Buf3 {
                x: up(dev, &lsq.x)?,
                y: up(dev, &lsq.y)?,
                z: up(dev, &lsq.z)?,
            };
            let ones = dev.alloc(nb)?;
            dev.fill(ones, 1.0)?;
            // the quadratic gradient on the cells where the face-neighbour
            // gradient is first order (D236)
            let quad = match quad {
                Some(q) => Some((q, [a3(dev)?, a3(dev)?, a3(dev)?])),
                None => None,
            };
            Some(QcrDev {
                mode: section.qcr.kernel(),
                tau: t,
                gtau,
                lsq,
                ones,
                zeros: dev.alloc(nb)?,
                quad,
            })
        };
        let ka_host = kinds(mesh, &case.boundary, treatment, VarKind::K, model).1;
        let wa_host = kinds(mesh, &case.boundary, treatment, VarKind::Second, model).1;
        let za_host = kinds(mesh, &case.boundary, treatment, VarKind::Zeta, model).1;
        let a3 = |dev: &mut dyn FlowDevice| -> Result<Buf3> {
            Ok(Buf3 {
                x: dev.alloc(nc)?,
                y: dev.alloc(nc)?,
                z: dev.alloc(nc)?,
            })
        };
        let gu = [a3(dev)?, a3(dev)?, a3(dev)?];
        let mu_lam = dev.alloc(nf)?;
        dev.copy(face_viscosity, mu_lam)?;
        let mask_b = {
            let b = dev.alloc(nb)?;
            dev.upload(b, &walls.mask_b)?;
            b
        };
        let wd = if walls.is_empty() {
            None
        } else {
            let nw = walls.len();
            let wd = WallDev {
                normal: Buf3 {
                    x: up(dev, &walls.normal[0])?,
                    y: up(dev, &walls.normal[1])?,
                    z: up(dev, &walls.normal[2])?,
                },
                y: up(dev, &walls.y)?,
                ks: up(dev, &walls.ks)?,
                cs: up(dev, &walls.cs)?,
                prod_w: dev.alloc(nw)?,
                weight: up(dev, &walls.weight)?,
                u_tau: dev.alloc(nw)?,
                y_plus: dev.alloc(nw)?,
                y_star: dev.alloc(nw)?,
                omega: dev.alloc(nw)?,
                blending: dev.alloc(nw)?,
                viscous_weight: dev.alloc(nw)?,
                pk: dev.alloc(nw)?,
            };
            dev.wall_setup(&walls.stat)?;
            Some(wd)
        };
        let fix_mask = up(dev, &vec![0.0; nc])?;
        let pressure: Vec<u32> = bc_arrays
            .pb_rule
            .iter()
            .enumerate()
            .filter(|(_, r)| **r != foehncraft_kernels::flow::PB_EXTRAPOLATE)
            .map(|(fb, _)| fb as u32)
            .collect();
        let pressure_owners = pressure
            .iter()
            .map(|&fb| mesh.owner[ni + fb as usize])
            .collect();
        let mut t = Self {
            section,
            settings,
            model,
            eps,
            c,
            o,
            ke,
            ko,
            c_mu,
            treatment: treatment.id(),
            nu,
            rho,
            fluids,
            nc,
            ni,
            nb,
            k,
            w,
            nu_t: dev.alloc(nc)?,
            nu_t_prev: None,
            f1: dev.alloc(nc)?,
            pk: dev.alloc(nc)?,
            dk: dev.alloc(nc)?,
            dw: dev.alloc(nc)?,
            capped: dev.alloc(nc)?,
            src: [
                dev.alloc(nc)?,
                dev.alloc(nc)?,
                dev.alloc(nc)?,
                dev.alloc(nc)?,
            ],
            d: dev.alloc(nc)?,
            gu,
            mu_lam,
            gamma_k: dev.alloc(nf)?,
            gamma_w: dev.alloc(nf)?,
            mask_b,
            walls,
            wd,
            fix_mask,
            fix_val: dev.alloc(nc)?,
            wall_pk: up(dev, &vec![0.0; nc])?,
            prod_mask: up(dev, &vec![0.0; nc])?,
            zero_face: dev.alloc(nf)?,
            cp: dev.alloc(nf)?,
            cn: dev.alloc(nf)?,
            bf: dev.alloc(nf)?,
            val: dev.alloc(stored)?,
            rhs: dev.alloc(nc)?,
            diag: dev.alloc(nc)?,
            order: dev.alloc(nc)?,
            work: dev.alloc(nc)?,
            inlets,
            zf,
            eb,
            curv,
            curv_fresh: CurvFresh::default(),
            qcr,
            amg: foehncraft_fv::amg::AmgSettings {
                scaled_pivots: true,
                ..case.solver.pressure.amg.clone()
            },
            elliptic: case.solver.elliptic.clone(),
            amg_rebuild: true,
            pseudo: None,
            pseudo_coef: dev.alloc(nc)?,
            amg_built: std::cell::RefCell::new(Vec::new()),
            transport_pc,
            reduction: case.solver.coupled.multigrid.clone(),
            reduction_setups: std::cell::Cell::new(0),
            reduction_fallbacks: std::cell::Cell::new(0),
            note: std::cell::RefCell::new(None),
            ka_host,
            kin_host: vec![0.0; nb],
            wa_host,
            win_host: vec![0.0; nb],
            za_host,
            zin_host: vec![0.0; nb],
            backflow,
            bf_state: Vec::new(),
            pressure_faces: pressure,
            pressure_owners,
            added: None,
            added_z: None,
            zero_cell: dev.alloc(nc)?,
            t_set: f64::NAN,
            wall_distance: Vec::new(),
            health: TurbHealth::default(),
        };
        if t.settings.viscosity_relaxation < 1.0 {
            t.nu_t_prev = Some(dev.alloc(nc)?);
        }
        // the wall distance (D48), unless the caller has it (the sliding
        // mesh by imprinting keeps it per part, D182)
        // (given first: a rank's part may hold no wall face, MPI step W5)
        let d = if let Some(d) = wall_distance {
            anyhow::ensure!(
                d.len() == nc,
                "a wall distance of {} cells for {nc}",
                d.len()
            );
            d
        } else if t.walls.is_empty() {
            vec![1e10; nc]
        } else {
            match t.section.wall_distance {
                WallDistanceMethod::Exact => crate::derived::wall_distance(case, mesh, geo)?,
                WallDistanceMethod::Pde => t.poisson_distance(dev, mesh, geo, disc, case)?,
            }
        };
        dev.upload(t.d, &d)?;
        t.wall_distance = d;
        // the start values (specification 7.3)
        let (mut k0, mut w0) = t.initial_values(dev, case, mesh, geo, bc_arrays, &env)?;
        // the start values of the regions (D238)
        crate::start::turbulence(case, mesh, geo, &env, t.eps, t.c_mu, &mut k0, &mut w0)?;
        dev.upload(t.k.phi, &k0)?;
        dev.upload(t.w.phi, &w0)?;
        if let Some(z) = &mut t.zf {
            // ζ = 2/3 (isotropy) unless given, f = 0 at the start
            // (specification 7.3); α = 1, its value away from walls
            dev.fill(z.z.phi, case.initial.zeta.unwrap_or(2.0 / 3.0))?;
            dev.fill(z.f.phi, if model.is_eb() { 1.0 } else { 0.0 })?;
            z.z.floor = kzf::ZETA_MIN;
        }
        let ne = equations(model).len();
        t.health.floor_iterations = vec![0; ne];
        t.health.floor_iterations_late = vec![0; ne];
        t.health.floor_cells_max = vec![0; ne];
        let kref = k0
            .iter()
            .copied()
            .fold(0.0, f64::max)
            .max(t.inlet_k_max(bc_arrays));
        let kref = if kref > 0.0 { kref } else { 1e-6 };
        t.k.floor = 1e-10 * kref;
        // ω_min and ε_min from k_min and the viscosity ratio 10⁻²
        t.w.floor = if eps {
            t.c_mu * t.k.floor * t.k.floor / (1e-2 * nu)
        } else {
            t.k.floor / (1e-2 * nu)
        };
        Ok(t)
    }

    /// The largest inlet k (at the prescribed boundary velocity).
    fn inlet_k_max(&self, bc: &BoundaryArrays) -> f64 {
        self.inlets
            .faces
            .iter()
            .map(|&(fb, s, x)| {
                let fb = fb as usize;
                let u = [bc.ub[0][fb], bc.ub[1][fb], bc.ub[2][fb]];
                let sp = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
                self.inlets.specs[s]
                    .at(sp, x, 0.0, self.nu_b(fb), self.c_mu, self.eps)
                    .0
            })
            .fold(0.0, f64::max)
    }

    /// The start values: `[initial]`, else the values of the first inlet
    /// (by name) at its mean speed, else I = 0.01 and r = 1 with the
    /// initial speed (1 m/s when that is 0).
    fn initial_values(
        &self,
        _dev: &mut dyn FlowDevice,
        case: &crate::case::FlowCase,
        mesh: &Mesh,
        geo: &Geometry,
        bc: &BoundaryArrays,
        env: &Env,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let nc = self.nc;
        let i = &case.initial;
        let u0 = i.velocity;
        let speed0 = (u0[0] * u0[0] + u0[1] * u0[1] + u0[2] * u0[2]).sqrt();
        let speed0 = if speed0 > 0.0 { speed0 } else { 1.0 };
        let cells = &geo.cell_centroid[..nc];
        let second = match (&i.omega, &i.epsilon) {
            (Some(w), _) => Some((w, "initial.omega", false)),
            (None, Some(e)) => Some((e, "initial.epsilon", true)),
            _ => None,
        };
        if let (Some(k), Some((w, key, given_eps))) = (&i.k, second) {
            let k = expr::value(k, env, Vars::SPACE, "initial.k")?;
            let w = expr::value(w, env, Vars::SPACE, key)?;
            let kv: Vec<f64> = cells.iter().map(|&x| k.at(x, 0.0)).collect();
            let mut wv: Vec<f64> = cells.iter().map(|&x| w.at(x, 0.0)).collect();
            // ε = C_μ k ω when the model solves the other one
            if given_eps != self.eps {
                for (w, &k) in wv.iter_mut().zip(&kv) {
                    *w = if self.eps {
                        self.c_mu * k * *w
                    } else {
                        *w / (self.c_mu * k).max(1e-300)
                    };
                }
            }
            return Ok((kv, wv));
        }
        let from = |ii: f64, r: f64, sp: f64, nu: f64| {
            let k = 1.5 * (ii * sp) * (ii * sp);
            if self.eps {
                (k, self.c_mu * k * k / (r * nu))
            } else {
                (k, k / (r * nu))
            }
        };
        // the start values of one fluid part (`in_part`: its boundary
        // faces; ν its viscosity): the initial intensity and viscosity
        // ratio, else the area-free mean of its first inlet's face values
        let start = |in_part: &dyn Fn(usize) -> bool, nu: f64| -> (f64, f64) {
            if i.turbulence_intensity.is_some() || i.viscosity_ratio.is_some() {
                return from(
                    i.turbulence_intensity.unwrap_or(0.01),
                    i.viscosity_ratio.unwrap_or(1.0),
                    speed0,
                    nu,
                );
            }
            let Some(&(_, s, _)) = self.inlets.faces.iter().find(|f| in_part(f.0 as usize)) else {
                return from(0.01, 1.0, speed0, nu);
            };
            let (mut ks, mut ws, mut n) = (0.0, 0.0, 0.0);
            for &(fb, si, x) in &self.inlets.faces {
                if si != s || !in_part(fb as usize) {
                    continue;
                }
                let fb = fb as usize;
                let u = [bc.ub[0][fb], bc.ub[1][fb], bc.ub[2][fb]];
                let sp = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
                let sp = if sp > 0.0 { sp } else { speed0 };
                let (kk, ww) =
                    self.inlets.specs[s].at(sp, x, 0.0, self.nu_b(fb), self.c_mu, self.eps);
                ks += kk;
                ws += ww;
                n += 1.0;
            }
            (ks / n, ws / n)
        };
        let Some(cf) = &self.fluids else {
            let (k, w) = start(&|_| true, self.nu);
            return Ok((vec![k; nc], vec![w; nc]));
        };
        // several fluids (D194): each fluid part from its own inlet
        let parts = crate::fluids::parts(mesh);
        let ni = mesh.n_interior_faces();
        let (mut kv, mut wv) = (vec![0.0; nc], vec![0.0; nc]);
        let mut done = vec![false; parts.n];
        for c in 0..nc {
            let pc = parts.part[c];
            if done[pc as usize] {
                continue;
            }
            done[pc as usize] = true;
            let in_part = |fb: usize| parts.part[mesh.owner[ni + fb] as usize] == pc;
            let (k, w) = start(&in_part, cf.nu[c]);
            for d in c..nc {
                if parts.part[d] == pc {
                    kv[d] = k;
                    wv[d] = w;
                }
            }
        }
        Ok((kv, wv))
    }

    /// The Poisson wall distance (Tucker 2003): ∇²φ = −1, φ = 0 on the
    /// walls, zero normal gradient elsewhere, solved with the scalar
    /// machinery (non-orthogonal corrections in sweeps).
    fn poisson_distance(
        &mut self,
        dev: &mut dyn FlowDevice,
        mesh: &Mesh,
        geo: &Geometry,
        disc: &foehncraft_fv::disc::Discretisation,
        case: &crate::case::FlowCase,
    ) -> Result<Vec<f64>> {
        let nc = self.nc;
        let ni = disc.n_interior;
        let mut points = BTreeMap::new();
        let mut a = vec![1.0; self.nb];
        for p in &mesh.boundaries {
            let wall = matches!(
                case.boundary.get(&p.name),
                Some(BoundaryCondition::Wall { .. })
            );
            points.insert(
                p.name.clone(),
                if wall {
                    BoundaryKind::FixedValue
                } else {
                    BoundaryKind::ZeroGradient
                },
            );
            if wall {
                for f in p.range() {
                    a[f - ni] = 0.0;
                }
            }
        }
        let lsq = crate::simplec::disc_with(mesh, geo, &points)?.lsq;
        let lsq = Buf3 {
            x: up(dev, &lsq.x)?,
            y: up(dev, &lsq.y)?,
            z: up(dev, &lsq.z)?,
        };
        let bc_a = up(dev, &a)?;
        let bc_b = up(dev, &vec![0.0; self.nb])?;
        let gamma = up(dev, &vec![1.0; disc.n_faces])?;
        let su = match &self.fluids {
            Some(cf) => up(dev, &cf.rho.iter().map(|r| 1.0 / r).collect::<Vec<_>>())?,
            None => up(dev, &vec![1.0 / self.rho; nc])?,
        };
        let phi = dev.alloc(nc)?;
        let psi = dev.alloc(nc)?;
        dev.fill(psi, 1.0)?;
        let g = Buf3 {
            x: dev.alloc(nc)?,
            y: dev.alloc(nc)?,
            z: dev.alloc(nc)?,
        };
        for _ in 0..8 {
            dev.gradient_with(lsq, phi, bc_a, bc_b, None, g)?;
            dev.scalar_faces(&ScalarFaceArgs {
                bc_a,
                bc_b,
                jump: None,
                mass_flux: self.zero_face,
                gamma,
                phi,
                grad: g,
                psi,
                blend: 0.0,
                central: 0.0,
                beta: 0.1,
                nonorth: 1.0,
                bound: 0.0,
                cp: self.cp,
                cn: self.cn,
                bf: self.bf,
            })?;
            dev.scalar_rows(&ScalarRowArgs {
                cp: self.cp,
                cn: self.cn,
                bf: self.bf,
                mass_flux: self.zero_face,
                phi,
                su: Some(su),
                sp: None,
                alpha: 1.0,
                time: None,
                fix: None,
                positive: false,
                pseudo: None,
                capacity: None,
                val: self.val,
                rhs: self.rhs,
                diag: self.diag,
                order: self.order,
            })?;
            dev.solve(
                self.val,
                self.rhs,
                phi,
                Preconditioner::Jacobi,
                false,
                1e-12,
                2000,
            )?;
        }
        dev.gradient_with(lsq, phi, bc_a, bc_b, None, g)?;
        let p = dev.download(phi)?;
        let (gx, gy, gz) = (dev.download(g.x)?, dev.download(g.y)?, dev.download(g.z)?);
        Ok(wall::poisson_distance(&p, [&gx, &gy, &gz]))
    }

    /// Sets the wall distance (verification: a manufactured solution with
    /// its own distance field).
    pub fn set_wall_distance(&mut self, dev: &mut dyn FlowDevice, d: &[f64]) -> Result<()> {
        anyhow::ensure!(d.len() == self.nc, "wall distance: one value per cell");
        dev.upload(self.d, d)?;
        self.wall_distance = d.to_vec();
        Ok(())
    }

    /// Adds sources of k and ω per unit mass (m²/s³ and 1/s²; a
    /// manufactured solution), split by sign into the rows.
    pub fn set_added_sources(
        &mut self,
        dev: &mut dyn FlowDevice,
        sk: &[f64],
        sw: &[f64],
    ) -> Result<()> {
        let (a, b) = match self.added {
            Some(x) => x,
            None => (dev.alloc(self.nc)?, dev.alloc(self.nc)?),
        };
        dev.upload(a, sk)?;
        dev.upload(b, sw)?;
        self.added = Some((a, b));
        Ok(())
    }

    /// Sets k and ω (verification and tests).
    pub fn set_fields(&mut self, dev: &mut dyn FlowDevice, k: &[f64], w: &[f64]) -> Result<()> {
        dev.upload(self.k.phi, k)?;
        dev.upload(self.w.phi, w)
    }

    /// Adds a source of ζ (φ of EB, 1/s) and a right-hand side of the
    /// elliptic equation of f (α) (a manufactured solution); with
    /// [`Self::set_added_sources`], which they need.
    pub fn set_added_zeta_sources(
        &mut self,
        dev: &mut dyn FlowDevice,
        sz: &[f64],
        sf: &[f64],
    ) -> Result<()> {
        anyhow::ensure!(self.zf.is_some(), "added ζ sources: the model has no ζ");
        let (a, b) = match self.added_z {
            Some(x) => x,
            None => (dev.alloc(self.nc)?, dev.alloc(self.nc)?),
        };
        dev.upload(a, sz)?;
        dev.upload(b, sf)?;
        self.added_z = Some((a, b));
        Ok(())
    }

    /// Sets ζ and f (φ and α of EB; verification and tests).
    pub fn set_zeta_fields(
        &mut self,
        dev: &mut dyn FlowDevice,
        z: &[f64],
        f: &[f64],
    ) -> Result<()> {
        let Some(zf) = self.zf.as_ref() else {
            anyhow::bail!("ζ fields: the model has no ζ");
        };
        let (zp, fp) = (zf.z.phi, zf.f.phi);
        dev.upload(zp, z)?;
        dev.upload(fp, f)
    }

    /// The inlet values of k and ω for the boundary velocity `ub` at time
    /// `t` (only when they depend on the speed or on t, unless `all`).
    fn set_inlets(&mut self, dev: &mut dyn FlowDevice, ub: Buf3, t: f64, all: bool) -> Result<()> {
        if self.inlets.faces.is_empty() {
            return Ok(());
        }
        if !all && !self.inlets.speed && !(self.inlets.time && t != self.t_set) {
            return Ok(());
        }
        let idx: Vec<u32> = self.inlets.faces.iter().map(|f| f.0).collect();
        let (ux, uy, uz) = if self.inlets.speed {
            (
                dev.gather(ub.x, &idx)?,
                dev.gather(ub.y, &idx)?,
                dev.gather(ub.z, &idx)?,
            )
        } else {
            let z = vec![0.0; idx.len()];
            (z.clone(), z.clone(), z)
        };
        let mut kb = vec![0.0; self.nb];
        let mut wb = vec![0.0; self.nb];
        for (i, &(fb, s, x)) in self.inlets.faces.iter().enumerate() {
            let sp = (ux[i] * ux[i] + uy[i] * uy[i] + uz[i] * uz[i]).sqrt();
            let (k, w) =
                self.inlets.specs[s].at(sp, x, t, self.nu_b(fb as usize), self.c_mu, self.eps);
            kb[fb as usize] = k.max(self.k.floor);
            wb[fb as usize] = w.max(self.w.floor);
        }
        dev.upload(self.k.b_in, &kb)?;
        dev.upload(self.w.b_in, &wb)?;
        if let Some(z) = self.zf {
            let mut zb = vec![0.0; self.nb];
            for &(fb, s, x) in &self.inlets.faces {
                zb[fb as usize] = self.inlets.specs[s]
                    .zeta_at(x)
                    .unwrap_or(self.inlets.zeta[s]);
            }
            dev.upload(z.z.b_in, &zb)?;
            self.zin_host = zb;
        }
        self.kin_host = kb;
        self.win_host = wb;
        self.t_set = t;
        if !self.bf_state.is_empty() {
            let st = std::mem::take(&mut self.bf_state);
            self.apply_backflow(dev, &st)?;
        }
        Ok(())
    }

    /// The backflow faces of the outlets (BC028) for the fluxes
    /// `mass_flux`: a face where the flow enters (ṁ < 0) takes the given
    /// k, ω or ε and ζ as fixed values (a = 0, b the value), a face where
    /// it leaves zero gradient (a = 1, b = 0). Called at the start of every
    /// turbulence solve, so a continued run sets the same state from the
    /// same fluxes.
    fn set_backflow(
        &mut self,
        dev: &mut dyn FlowDevice,
        ub: Buf3,
        mass_flux: Buf,
        t: f64,
    ) -> Result<()> {
        if self.backflow.faces.is_empty() {
            return Ok(());
        }
        let ni = self.ni as u32;
        let fidx: Vec<u32> = self.backflow.faces.iter().map(|f| f.0 + ni).collect();
        let m = dev.gather(mass_flux, &fidx)?;
        let bidx: Vec<u32> = self.backflow.faces.iter().map(|f| f.0).collect();
        let (ux, uy, uz) = if self.backflow.speed {
            (
                dev.gather(ub.x, &bidx)?,
                dev.gather(ub.y, &bidx)?,
                dev.gather(ub.z, &bidx)?,
            )
        } else {
            let z = vec![0.0; bidx.len()];
            (z.clone(), z.clone(), z)
        };
        let mut st = Vec::with_capacity(bidx.len());
        for (i, &(fb, s, x)) in self.backflow.faces.iter().enumerate() {
            let v = if m[i] < 0.0 {
                let sp = (ux[i] * ux[i] + uy[i] * uy[i] + uz[i] * uz[i]).sqrt();
                let (k, w) = self.backflow.specs[s].at(
                    sp,
                    x,
                    t,
                    self.nu_b(fb as usize),
                    self.c_mu,
                    self.eps,
                );
                let z = self.backflow.specs[s]
                    .zeta_at(x)
                    .unwrap_or(self.backflow.zeta[s]);
                Some((k.max(self.k.floor), w.max(self.w.floor), z))
            } else {
                None
            };
            st.push((fb as usize, v));
        }
        if st == self.bf_state {
            return Ok(());
        }
        self.apply_backflow(dev, &st)
    }

    /// Writes the backflow state `st` (boundary face, entering values or
    /// `None`) into a and b of k, ω or ε and ζ and uploads them.
    fn apply_backflow(&mut self, dev: &mut dyn FlowDevice, st: &[BackflowFace]) -> Result<()> {
        self.bf_state = st.to_vec();
        for &(fb, v) in st {
            let (a, (k, w, z)) = match v {
                Some(v) => (0.0, v),
                None => (1.0, (0.0, 0.0, 0.0)),
            };
            self.ka_host[fb] = a;
            self.wa_host[fb] = a;
            self.za_host[fb] = a;
            self.kin_host[fb] = k;
            self.win_host[fb] = w;
            self.zin_host[fb] = z;
        }
        dev.upload(self.k.bc_a, &self.ka_host)?;
        dev.upload(self.w.bc_a, &self.wa_host)?;
        dev.upload(self.k.b_in, &self.kin_host)?;
        dev.upload(self.w.b_in, &self.win_host)?;
        if let Some(z) = self.zf {
            dev.upload(z.z.bc_a, &self.za_host)?;
            dev.upload(z.z.b_in, &self.zin_host)?;
        }
        Ok(())
    }

    /// The transported and elliptic variables: k, ω or ε, and with
    /// k-ε-ζ-f ζ and f.
    fn vars(&self) -> Vec<Var> {
        let mut v = vec![self.k, self.w];
        if let Some(z) = self.zf {
            v.extend([z.z, z.f]);
        }
        v
    }

    /// bc_b = inlet values + wall values.
    fn boundary_values(&self, dev: &mut dyn FlowDevice) -> Result<()> {
        for v in self.vars() {
            dev.lincomb(v.bc_b, 1.0, v.b_in, 1.0, v.b_wall)?;
        }
        Ok(())
    }

    /// Gradients and limiters of the transported variables (f: the gradient
    /// only, for the non-orthogonal correction).
    fn gradients(&self, dev: &mut dyn FlowDevice) -> Result<()> {
        let s = &self.settings;
        if let Some(z) = self.zf {
            dev.gradient_with(z.f.lsq, z.f.phi, z.f.bc_a, z.f.bc_b, None, z.f.grad)?;
        }
        let mut tr = vec![self.k, self.w];
        if let Some(z) = self.zf {
            tr.push(z.z);
        }
        for v in &tr {
            dev.gradient_with(v.lsq, v.phi, v.bc_a, v.bc_b, None, v.grad)?;
            if s.scheme != Scheme::Upwind && s.limiter != LimiterKind::None {
                let eps2 = if s.limiter == LimiterKind::BarthJespersen {
                    Eps2::Single(-1.0)
                } else {
                    let (mn, mx) = dev.min_max(v.phi)?;
                    let e = s.limiter_k.max(1e-12) * (mx - mn);
                    Eps2::Single((e * e).max(1e-300))
                };
                dev.limiter_with(v.phi, v.bc_a, v.bc_b, None, v.grad, eps2, v.psi)?;
            }
        }
        Ok(())
    }

    /// The model per cell: with `sources`, also the source split.
    fn cells(&self, dev: &mut dyn FlowDevice, sources: bool) -> Result<()> {
        if let Some(z) = self.zf {
            return dev.zeta_f_cells(&ZetaFCellOps {
                grad_u: self.gu,
                k: self.k.phi,
                e: self.w.phi,
                zeta: z.z.phi,
                f: z.f.phi,
                nu: self.nu,
                nu_cell: self.fluids.as_ref().map(|f| f.buf),
                c: z.c,
                o: z.o,
                add: if sources {
                    self.added.map(|(a, b)| {
                        let (z, f) = self.added_z.unwrap_or((self.zero_cell, self.zero_cell));
                        [a, b, z, f]
                    })
                } else {
                    None
                },
                wall: (sources && self.wd.is_some() && !self.walls.stat.cells.is_empty())
                    .then_some((self.prod_mask, self.wall_pk)),
                dist: self.d,
                nu_t: self.nu_t,
                t: z.t,
                l2: z.l2,
                pk: self.pk,
                sources: sources.then_some([
                    self.src[0],
                    self.src[1],
                    self.src[2],
                    self.src[3],
                    z.src[0],
                    z.src[1],
                    z.src[2],
                    z.src[3],
                ]),
                dk: self.dk,
                de: self.dw,
                dz: z.dz,
                capped: self.capped,
                eb: self.eb.map(|e| (self.k.grad, z.z.grad, e.g2, e.dkt)),
            });
        }
        if self.eps {
            return dev.keps_cells(&KepsCellOps {
                grad_u: self.gu,
                k: self.k.phi,
                e: self.w.phi,
                nu: self.nu,
                nu_cell: self.fluids.as_ref().map(|f| f.buf),
                c: self.ke,
                o: self.ko,
                add: if sources { self.added } else { None },
                wall: (sources && self.wd.is_some() && !self.walls.stat.cells.is_empty())
                    .then_some((self.prod_mask, self.wall_pk)),
                nu_t: self.nu_t,
                c_mu: self.f1,
                pk: self.pk,
                sources: sources.then_some(self.src),
                dk: self.dk,
                de: self.dw,
                capped: self.capped,
            });
        }
        dev.sst_cells(&SstCellOps {
            grad_u: self.gu,
            k: self.k.phi,
            w: self.w.phi,
            grad_k: self.k.grad,
            grad_w: self.w.grad,
            d: self.d,
            nu: self.nu,
            nu_cell: self.fluids.as_ref().map(|f| f.buf),
            c: self.c,
            o: self.o,
            add: if sources { self.added } else { None },
            wall: (sources && self.treatment == kwall::AUTOMATIC && self.wd.is_some())
                .then_some((self.prod_mask, self.wall_pk)),
            fr1: self.curv.map(|c| c.fr1),
            nu_t: self.nu_t,
            f1: self.f1,
            pk: self.pk,
            sources: sources.then_some(self.src),
            dk: self.dk,
            dw: self.dw,
            capped: self.capped,
        })
    }

    /// The wall treatment: u_τ, y⁺, μ_w, the wall values of rough walls and
    /// the ω of the near-wall cells.
    fn walls(
        &self,
        dev: &mut dyn FlowDevice,
        u: Buf3,
        ub: Buf3,
        grad_p: Option<Buf3>,
    ) -> Result<()> {
        let Some(wd) = self.wd else { return Ok(()) };
        dev.wall_faces(&WallFaceOps {
            normal: wd.normal,
            y: wd.y,
            ks: wd.ks,
            u,
            ub,
            nu: self.nu,
            nu_cell: self.fluids.as_ref().map(|f| f.buf),
            treatment: self.treatment,
            c: self.c,
            ke: self.ke,
            k: self.k.phi,
            cs: wd.cs,
            u_tau: wd.u_tau,
            y_plus: wd.y_plus,
            y_star: wd.y_star,
            omega: wd.omega,
            blending: wd.blending,
            viscous_weight: wd.viscous_weight,
            pk: wd.pk,
            prod_w: wd.prod_w,
            mu_eff: dev.face_viscosity(),
            k_b: self.k.b_wall,
            w_b: self.w.b_wall,
            zeta_f: self.zf.map(|z| (z.z.phi, z.f.b_wall)),
            n_zeta: self.zf.map_or(1.0, |z| z.c.n_zeta),
            // C_U of the 2023 treatment (none at the start: ψ = 1 there)
            non_eq: (self.treatment == kwall::COMPOUND_2023)
                .then_some(grad_p)
                .flatten()
                .map(|gp| (self.gu, gp)),
        })?;
        dev.wall_cells(&WallCellOps {
            omega: wd.omega,
            weight: wd.weight,
            blending: wd.blending,
            pk: wd.pk,
            prod_w: wd.prod_w,
            wall_pk: self.wall_pk,
            prod_mask: self.prod_mask,
            fix_val: self.fix_val,
            fix_mask: self.fix_mask,
            phi: self.eps.then_some(self.w.phi),
        })
    }

    /// μ_eff on the faces off the walls and the diffusivities of k and ω.
    fn faces(&self, dev: &mut dyn FlowDevice) -> Result<()> {
        dev.turb_faces(&TurbFaceOps {
            nu_t: self.nu_t,
            dk: self.dk,
            dw: self.dw,
            mu: self.mu_lam,
            wall: self.mask_b,
            mu_eff: dev.face_viscosity(),
            gamma_k: self.gamma_k,
            gamma_w: self.gamma_w,
            zeta_f: self.zf.map(|z| [z.dz, z.l2, z.gamma_z, z.gamma_f]),
            mol: self.zf.map_or(1.0, |z| z.c.mol),
        })
    }

    /// The inputs of the elliptic-blending k-ε that need gradients of
    /// gradients: Σ (∂_l ∂_j U_i)² from the gradients of the nine
    /// velocity-gradient components, and ∇·(ν_t/σ_k ∇k) from the
    /// gradients of the flux q = ν_t/σ_k ∇k (ν_t of the last evaluation;
    /// the zero-gradient boundary values of [`EbDev`]).
    fn eb_inputs(&self, dev: &mut dyn FlowDevice) -> Result<()> {
        let Some(e) = self.eb else { return Ok(()) };
        dev.fill(e.g2, 0.0)?;
        for g in self.gu {
            for c in g.all() {
                dev.gradient_with(e.lsq, c, e.ones, e.zeros, None, e.tmp)?;
                dev.squares_add(e.tmp, e.g2)?;
            }
        }
        dev.scale_vector(self.dk, self.k.grad, e.q)?;
        dev.gradient_with(e.lsq, e.q.x, e.ones, e.zeros, None, e.tmp)?;
        dev.copy(e.tmp.x, e.dkt)?;
        dev.gradient_with(e.lsq, e.q.y, e.ones, e.zeros, None, e.tmp)?;
        dev.lincomb(e.dkt, 1.0, e.dkt, 1.0, e.tmp.y)?;
        dev.gradient_with(e.lsq, e.q.z, e.ones, e.zeros, None, e.tmp)?;
        dev.lincomb(e.dkt, 1.0, e.dkt, 1.0, e.tmp.z)?;
        Ok(())
    }

    /// f_r1 of the curvature correction (D170) from the velocity gradient
    /// of the last evaluation (D236: on the cells where it is first order
    /// the quadratic fit of the velocity `u`, `ub`, and the gradients of
    /// S_ij from a quadratic fit over interior cells where the
    /// least-squares gradient of S_ij is not consistent): S_ij per cell, the gradients of its six
    /// components (zero normal gradient on every boundary), DS_ij/Dt with
    /// the velocity `u` less the frame velocity, and in a transient run
    /// the BDF1 time part from the last time level (`dt`).
    fn curvature_inputs(
        &self,
        dev: &mut dyn FlowDevice,
        u: Buf3,
        ub: Buf3,
        dt: Option<f64>,
    ) -> Result<()> {
        let Some(c) = self.curv else { return Ok(()) };
        // S_ij from a velocity gradient that is second order everywhere
        // (D236): the quadratic fit where the face-neighbour one is not
        let g = match c.quad {
            Some((q, gq)) => {
                for (a, b) in self.gu.iter().zip(&gq) {
                    for (x, y) in a.all().into_iter().zip(b.all()) {
                        dev.copy(x, y)?;
                    }
                }
                dev.quadratic_gradient(q.slot, q.w, q.dims, u, ub, gq)?;
                gq
            }
            None => self.gu,
        };
        for (i, &b) in c.s.iter().enumerate() {
            dev.strain_component(g, i as u32, b)?;
        }
        self.curv_fresh.0.set(true);
        let old = |i: usize| match dt {
            Some(dt) if c.has_old => Some((c.s[i], c.s_old[i], 1.0 / dt)),
            _ => None,
        };
        let v_frame = c.frame.map(|f| f.1);
        match c.sfit {
            // the gradients of S_ij: least squares (zero normal gradient on
            // every boundary), and on the cells where that is not
            // consistent the quadratic fit over interior cells (D236)
            Some((q, zeros, gs)) => {
                for half in [0usize, 3] {
                    for (j, gj) in gs.iter().enumerate() {
                        dev.gradient_with(c.lsq, c.s[half + j], c.ones, c.zeros, None, *gj)?;
                    }
                    let three = Buf3 {
                        x: c.s[half],
                        y: c.s[half + 1],
                        z: c.s[half + 2],
                    };
                    dev.quadratic_gradient(q.slot, q.w, q.dims, three, zeros, gs)?;
                    for (j, gj) in gs.iter().enumerate() {
                        let i = half + j;
                        dev.convective_derivative(u, v_frame, *gj, old(i), c.dsdt[i])?;
                    }
                }
            }
            None => {
                for i in 0..6 {
                    dev.gradient_with(c.lsq, c.s[i], c.ones, c.zeros, None, c.tmp)?;
                    dev.convective_derivative(u, v_frame, c.tmp, old(i), c.dsdt[i])?;
                }
            }
        }
        dev.curvature_cells(&CurvatureOps {
            grad_u: self.gu,
            dsdt: c.dsdt,
            w: self.w.phi,
            omega_rot: c.frame.map(|f| f.0),
            rot: c.rot,
            spin: c.spin,
            c: self.c,
            smooth: self.o.smooth,
            out: c.fr1,
        })
    }

    /// The non-linear stress of the QCR (D170) in the momentum face parts
    /// `bf`, from the current velocity gradient and ν_t of the last
    /// turbulence update (an explicit stress; nothing without `turbulence.qcr`).
    ///
    /// The divergence is consistent on irregular cells (D236): the velocity
    /// gradient of the stress is the quadratic least-squares fit where the
    /// face-neighbour gradient is first order (cell velocity `u`, boundary
    /// velocity `ub`), and the face stress is corrected for skewness with
    /// the cell gradients of the stress.
    pub(crate) fn add_explicit_stress(
        &self,
        dev: &mut dyn FlowDevice,
        grad_u: [Buf3; 3],
        u: Buf3,
        ub: Buf3,
        bf: Buf3,
    ) -> Result<()> {
        let Some(q) = self.qcr else {
            return Ok(());
        };
        let g = match q.quad {
            Some((QuadDev { slot, w, dims, .. }, gq)) => {
                for (a, b) in grad_u.iter().zip(&gq) {
                    for (x, y) in a.all().into_iter().zip(b.all()) {
                        dev.copy(x, y)?;
                    }
                }
                dev.quadratic_gradient(slot, w, dims, u, ub, gq)?;
                gq
            }
            None => grad_u,
        };
        dev.qcr_cells(g, self.nu_t, &self.c, q.mode, q.tau)?;
        for (t, gt) in q.tau.iter().zip(&q.gtau) {
            dev.gradient_with(q.lsq, *t, q.ones, q.zeros, None, *gt)?;
        }
        dev.explicit_stress_faces(q.tau, Some(q.gtau), bf)
    }

    /// The moving reference frames of the curvature correction (D170):
    /// Ω^rot and Ω × (x − x₀) per cell (zero outside the turning zones).
    pub(crate) fn set_curvature_frames(
        &mut self,
        dev: &mut dyn FlowDevice,
        omega: &[[f64; 3]],
        velocity: &[[f64; 3]],
    ) -> Result<()> {
        let Some(c) = self.curv.as_mut() else {
            return Ok(());
        };
        let nc = self.nc;
        let up3 = |dev: &mut dyn FlowDevice, v: &[[f64; 3]]| -> Result<Buf3> {
            let b = Buf3 {
                x: dev.alloc(nc)?,
                y: dev.alloc(nc)?,
                z: dev.alloc(nc)?,
            };
            for (k, comp) in b.all().into_iter().enumerate() {
                let x: Vec<f64> = v.iter().map(|a| a[k]).collect();
                dev.upload(comp, &x)?;
            }
            Ok(b)
        };
        c.frame = Some((up3(dev, omega)?, up3(dev, velocity)?));
        Ok(())
    }

    /// The buffers a passive scalar couples to (D135): the laminar and the
    /// effective face viscosity, the wall faces and the y⁺ of the thermal
    /// law (y* = u* y_P/ν, u* = C_μ^¼ k_P^½, for the k-ε family and
    /// ζ-f/EB; the blended u_τ for SST; section 5.2 of the stage-6
    /// specification).
    pub(crate) fn scalar_coupling(&self, dev: &dyn FlowDevice) -> crate::scalar::TurbCoupling {
        crate::scalar::TurbCoupling {
            mu: self.mu_lam,
            mu_eff: dev.face_viscosity(),
            wall: self.mask_b,
            y_plus: self.wd.map(|w| w.y_star),
            y_floor: if self.treatment == kwall::SCALABLE {
                kwall::Y_PLUS_LOG
            } else {
                0.0
            },
        }
    }

    /// The wall faces of the flow (the face index of each wall face, in the
    /// order of the per-wall-face buffers such as the y⁺ of
    /// [`scalar_coupling`](Self::scalar_coupling)) and their first-cell
    /// distances y_P (conjugate heat transfer, stage 7d5).
    pub(crate) fn wall_faces(&self) -> (&[u32], &[f64]) {
        (&self.walls.stat.faces, &self.walls.y)
    }

    /// The velocity gradient of `u` with the boundary values `ub`.
    fn velocity_gradient(&self, dev: &mut dyn FlowDevice, u: Buf3, ub: Buf3) -> Result<()> {
        for ((ui, bi), g) in u.all().into_iter().zip(ub.all()).zip(self.gu) {
            dev.gradient(Lsq::Velocity, ui, Mask::Velocity, bi, g)?;
        }
        Ok(())
    }

    /// The face viscosity of the start (before the first iteration and after
    /// a restart without turbulence state).
    pub(crate) fn start(
        &mut self,
        dev: &mut dyn FlowDevice,
        u: Buf3,
        ub: Buf3,
        bc: &BoundaryArrays,
    ) -> Result<()> {
        self.set_inlets(dev, ub, 0.0, true)?;
        self.velocity_gradient(dev, u, ub)?;
        self.walls(dev, u, ub, None)?;
        self.boundary_values(dev)?;
        self.gradients(dev)?;
        self.cells(dev, false)?;
        self.faces(dev)?;
        self.modified_pressure(dev, bc)
    }

    /// The pressure boundaries take p̂ = p + ⅔ρk with k of the face (also
    /// after a restart: the values of [`Self::start`] came from the start
    /// fields).
    pub(crate) fn modified_pressure(
        &self,
        dev: &mut dyn FlowDevice,
        bc: &BoundaryArrays,
    ) -> Result<()> {
        if self.pressure_faces.is_empty() {
            return Ok(());
        }
        let kp = dev.gather(self.k.phi, &self.pressure_owners)?;
        let mut v = bc.pset.clone();
        for (i, &fb) in self.pressure_faces.iter().enumerate() {
            let fb = fb as usize;
            let kf = self.ka_host[fb] * kp[i] + self.kin_host[fb];
            v[fb] += 2.0 / 3.0 * self.rho_b(fb) * kf;
        }
        dev.set_pressure_values(&v)
    }

    /// The turbulence solves of one outer iteration with the corrected
    /// velocity `u` and the fluxes `mass_flux`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve(
        &mut self,
        dev: &mut dyn FlowDevice,
        u: Buf3,
        ub: Buf3,
        mass_flux: Buf,
        grad_p: Buf3,
        bc: &BoundaryArrays,
        step: Option<StepTime>,
    ) -> Result<TurbResiduals> {
        let tr = step.is_some();
        let t = step.map_or(0.0, |s| s.time);
        self.set_inlets(dev, ub, t, false)?;
        self.set_backflow(dev, ub, mass_flux, t)?;
        // ν_t of the last iteration (relaxation of ν_t)
        if let Some(prev) = self.nu_t_prev {
            dev.copy(self.nu_t, prev)?;
        }
        // 1. velocity gradient of the corrected velocity, the wall
        // treatment, gradients of k and ω
        self.velocity_gradient(dev, u, ub)?;
        self.walls(dev, u, ub, Some(grad_p))?;
        self.boundary_values(dev)?;
        self.gradients(dev)?;
        self.eb_inputs(dev)?;
        self.curvature_inputs(dev, u, ub, step.map(|s| s.dt))?;
        // 2. the sources and the diffusivities from the fields at the start
        self.cells(dev, true)?;
        self.faces(dev)?;
        let s = self.settings.clone();
        // the pseudo-time step of a steady iteration: per cell
        // min(Δτ_max, c_T T) with the time scale T of the model (k/ε,
        // 1/(β* ω), T of k-ε-ζ-f), from the fields at the start
        let pseudo = match self.pseudo.filter(|_| !tr) {
            Some(p) => {
                let (factor, num, den) = match self.zf {
                    Some(z) => (p.rate_factor, None, Some(z.t)),
                    None if self.eps => (p.rate_factor, Some(self.w.phi), Some(self.k.phi)),
                    None => (p.rate_factor * self.c.beta_star, Some(self.w.phi), None),
                };
                dev.pseudo_coefficient(self.pseudo_coef, p.inv_dt, factor, num, den)?;
                Some(self.pseudo_coef)
            }
            None => None,
        };
        let alpha = if pseudo.is_some() {
            s.relaxation
                .unwrap_or(if self.transport_pc == TransportPreconditioner::Reduction {
                    PSEUDO_RELAXATION_REDUCTION
                } else {
                    PSEUDO_RELAXATION
                })
        } else {
            s.alpha(tr)
        };
        let mut out = TurbResiduals::default();
        let fix_second = (self.treatment == kwall::AUTOMATIC || self.eps)
            && self.wd.is_some()
            && !self.walls.stat.cells.is_empty();
        // k, then ω or ε (the near-wall rows of the second take the fixed
        // value)
        for (e, (v, gamma)) in [(self.k, self.gamma_k), (self.w, self.gamma_w)]
            .into_iter()
            .enumerate()
        {
            let (su, sp) = (self.src[2 * e], self.src[2 * e + 1]);
            let fix = (e == 1 && fix_second).then_some((self.fix_mask, self.fix_val));
            let r = self.transport(
                dev,
                v,
                gamma,
                (su, sp),
                fix,
                mass_flux,
                step,
                alpha,
                2 + e,
                pseudo,
            )?;
            out.residual.push(r.0);
            out.iterations.push(r.1);
            out.floors.push(r.2);
        }
        if let Some(z) = self.zf {
            // k-ε-ζ-f (specification 6.4): f from the new k and ε, then ζ
            // with the new f
            self.cells(dev, true)?;
            self.faces(dev)?;
            let r = self.elliptic(dev, z, tr)?;
            out.residual.push(r.0);
            out.iterations.push(r.1);
            // α in [0, 1] (the discrete maximum principle holds only on
            // M-matrices, and the solve stops at its tolerance: not counted
            // as a floor, since α = 1 ± the tolerance in the outer flow)
            if self.eb.is_some() {
                dev.floor(z.f.phi, 0.0, 1.0, None, self.work)?;
            }
            out.floors.push(0);
            self.cells(dev, true)?;
            let r = self.transport(
                dev,
                z.z,
                z.gamma_z,
                (z.src[0], z.src[1]),
                None,
                mass_flux,
                step,
                alpha,
                4,
                pseudo,
            )?;
            out.residual.push(r.0);
            out.iterations.push(r.1);
            out.floors.push(r.2);
        }
        // 3. ν_t from the new fields
        self.cells(dev, false)?;
        if let Some(prev) = self.nu_t_prev {
            let a = s.viscosity_relaxation;
            dev.lincomb(self.nu_t, a, self.nu_t, 1.0 - a, prev)?;
        }
        out.capped = dev.sum_sq(self.capped)?.round() as usize;
        // the floor counts are sums over all ranks, so is the share they
        // are measured against (MPI step W7)
        let nc = match dev.comm() {
            Some(c) => {
                let mut n = [self.n_cells() as i64];
                c.all_sum_i64(&mut n);
                n[0] as usize
            }
            None => self.n_cells(),
        };
        let mut flag = 0u8;
        for (e, &n) in out.floors.iter().enumerate() {
            if n * 1000 > nc {
                self.health.floor_iterations[e] += 1;
                flag |= 1 << e;
            }
            self.health.floor_cells_max[e] = self.health.floor_cells_max[e].max(n);
        }
        let h = &mut self.health;
        h.flags.push(flag);
        let tail = h.flags.len().div_ceil(10);
        for e in 0..out.floors.len() {
            h.floor_iterations_late[e] = h.flags[h.flags.len() - tail..]
                .iter()
                .filter(|&&f| f & (1 << e) != 0)
                .count();
        }
        self.health.capped_max = self.health.capped_max.max(out.capped);
        self.faces(dev)?;
        self.modified_pressure(dev, bc)?;
        Ok(out)
    }

    /// Assembles and solves one transport equation (k, ω, ε or ζ) with the
    /// sources `src` (S_u, S_p) and the face diffusivity `gamma`, then
    /// applies its bounds. Returns (scaled residual, iterations, cells
    /// bounded).
    #[allow(clippy::too_many_arguments)]
    fn transport(
        &self,
        dev: &mut dyn FlowDevice,
        v: Var,
        gamma: Buf,
        src: (Buf, Buf),
        fix: Option<(Buf, Buf)>,
        mass_flux: Buf,
        step: Option<StepTime>,
        alpha: f64,
        slot: usize,
        pseudo: Option<Buf>,
    ) -> Result<(f64, usize, usize)> {
        let s = &self.settings;
        let tr = step.is_some();
        dev.scalar_faces(&ScalarFaceArgs {
            bc_a: v.bc_a,
            bc_b: v.bc_b,
            jump: None,
            mass_flux,
            gamma,
            phi: v.phi,
            grad: v.grad,
            psi: v.psi,
            blend: if s.scheme == Scheme::Upwind { 0.0 } else { 1.0 },
            central: if s.scheme == Scheme::Central {
                1.0
            } else {
                0.0
            },
            beta: 0.1,
            nonorth: 1.0,
            bound: 0.0,
            cp: self.cp,
            cn: self.cn,
            bf: self.bf,
        })?;
        let time = match (step, v.old) {
            (Some(st), Some((o1, o2, _))) => Some(ScalarTime {
                inv_dt: 1.0 / st.dt,
                a: st.a,
                old: v.content.map_or(o1, |c| c.0),
                old2: v.content.map_or(o2, |c| c.1),
                bounds: (0.0, f64::INFINITY),
            }),
            _ => None,
        };
        dev.scalar_rows(&ScalarRowArgs {
            cp: self.cp,
            cn: self.cn,
            bf: self.bf,
            mass_flux,
            phi: v.phi,
            su: Some(src.0),
            sp: Some(src.1),
            alpha,
            time,
            fix,
            positive: true,
            pseudo,
            capacity: None,
            val: self.val,
            rhs: self.rhs,
            diag: self.diag,
            order: self.order,
        })?;
        let scale = dev.scalar_scale(v.phi, self.rhs, self.diag, alpha, pseudo)?;
        let rel_old = match self.zf {
            Some(z) => {
                dev.copy(v.phi, z.prev)?;
                Some((z.prev, ZETA_F_DROP))
            }
            None => None,
        };
        // its own hierarchy (slot 2 k, 3 ε or ω, 4 ζ): the aggregation
        // multigrid from the symmetric part, or the reduction multigrid
        let (pc, kept) = match self.transport_pc {
            TransportPreconditioner::Jacobi => (Preconditioner::Jacobi, false),
            TransportPreconditioner::Amg => {
                if self.amg_prepare(dev, slot, true)? {
                    (Preconditioner::Amg, false)
                } else {
                    (Preconditioner::Jacobi, false)
                }
            }
            TransportPreconditioner::Reduction => match self.reduction_prepare(dev, slot)? {
                Some(kept) => (Preconditioner::Amg, kept),
                None => (Preconditioner::Jacobi, false),
            },
        };
        let solved = dev.solve(self.val, self.rhs, v.phi, pc, false, s.tol(tr), s.max_iter);
        if pc == Preconditioner::Amg {
            dev.amg_select(0)?;
        }
        let (it, r0) = solved?;
        if kept && it > REDUCTION_STALE {
            // a kept hierarchy that no longer serves: rebuilt next time
            self.amg_built.borrow_mut()[slot] = None;
        }
        let hi = if self.zf.is_some_and(|z| z.z.phi == v.phi) {
            kzf::ZETA_MAX
        } else {
            f64::INFINITY
        };
        let n = dev.floor(v.phi, v.floor, hi, rel_old, self.work)?;
        Ok((if scale > 0.0 { r0 / scale } else { r0 }, it, n))
    }

    /// Selects the hierarchy slot `slot` of an equation, rebuilt from the
    /// assembled matrix on the schedule of the pressure (`amg_rebuild`) or
    /// when the slot has none yet. Between rebuilds a hierarchy of an
    /// earlier matrix of the same equation preconditions the solve (a
    /// fixed preconditioner within each solve: the Krylov method stays
    /// valid, only its rate changes).
    ///
    /// Returns false when the setup failed (the matrix is not positive
    /// definite to the precision of the setup): the solve then takes the
    /// Jacobi preconditioner, the first failure is noted for the log, and
    /// the next solve tries a new setup.
    fn amg_prepare(
        &self,
        dev: &mut dyn FlowDevice,
        slot: usize,
        symmetric_part: bool,
    ) -> Result<bool> {
        let mut built = self.amg_built.borrow_mut();
        if built.len() <= slot {
            built.resize(slot + 1, None);
        }
        if self.amg_rebuild || built[slot].is_none() {
            match dev.amg_setup_slot(slot, self.val, &self.amg, symmetric_part) {
                Ok(()) => built[slot] = Some(dev.download(self.val)?),
                Err(e) => {
                    built[slot] = None;
                    dev.amg_select(0)?;
                    let n = self.reduction_fallbacks.get() + 1;
                    self.reduction_fallbacks.set(n);
                    if n == 1 {
                        *self.note.borrow_mut() = Some(format!(
                            "the multigrid setup of a turbulence equation (slot {slot}) failed \
                             ({e:#}); its solves take the Jacobi preconditioner until a setup \
                             succeeds"
                        ));
                    }
                    return Ok(false);
                }
            }
        } else {
            dev.amg_select(slot)?;
        }
        Ok(true)
    }

    /// The reduction multigrid (D180) of an equation in slot `slot`: the
    /// kept hierarchy while the assembled matrix differs from the one it
    /// was built from by at most `reuse`, else a new one. Returns whether
    /// the hierarchy was kept, or None when the setup failed (the solve
    /// then takes the Jacobi preconditioner; the first failure is noted
    /// for the log).
    fn reduction_prepare(&self, dev: &mut dyn FlowDevice, slot: usize) -> Result<Option<bool>> {
        let mut built = self.amg_built.borrow_mut();
        if built.len() <= slot {
            built.resize(slot + 1, None);
        }
        let cur = dev.download(self.val)?;
        let comm = dev.comm();
        let keep = built[slot].as_ref().is_some_and(|m| {
            crate::coupled::matrix_change(&cur, m, None, &[], comm.as_ref()) <= self.settings.reuse
        });
        if keep {
            dev.amg_select(slot)?;
            return Ok(Some(true));
        }
        match dev.reduction_setup_slot(slot, self.val, &self.reduction) {
            Ok(()) => {
                built[slot] = Some(cur);
                self.reduction_setups.set(self.reduction_setups.get() + 1);
                Ok(Some(false))
            }
            Err(e) => {
                built[slot] = None;
                dev.amg_select(0)?;
                let n = self.reduction_fallbacks.get() + 1;
                self.reduction_fallbacks.set(n);
                if n == 1 {
                    *self.note.borrow_mut() = Some(format!(
                        "the multigrid setup of a turbulence equation failed ({e:#}); its \
                         solves take the Jacobi preconditioner until a setup succeeds"
                    ));
                }
                Ok(None)
            }
        }
    }

    /// Setups of the reduction multigrid of the transport equations.
    pub fn reduction_setups(&self) -> usize {
        self.reduction_setups.get()
    }

    /// The message of the first failed multigrid setup, once.
    pub fn take_note(&self) -> Option<String> {
        self.note.borrow_mut().take()
    }

    /// Rebuilds the hierarchies of a loaded state from their matrices
    /// (`TurbState::amg`); the slot of the f equation is not symmetrized,
    /// the others are (as in [`Self::amg_prepare`]), or carry the
    /// reduction multigrid.
    fn amg_restore(&self, dev: &mut dyn FlowDevice, amg: &[(usize, Vec<f64>)]) -> Result<()> {
        let mut built = self.amg_built.borrow_mut();
        built.clear();
        for (slot, m) in amg {
            anyhow::ensure!(
                (1..=4).contains(slot),
                "turbulence multigrid slot {slot} outside 1..=4"
            );
            dev.upload(self.val, m)?;
            if *slot != 1 && self.transport_pc == TransportPreconditioner::Reduction {
                dev.reduction_setup_slot(*slot, self.val, &self.reduction)?;
            } else {
                dev.amg_setup_slot(*slot, self.val, &self.amg, *slot != 1)?;
            }
            if built.len() <= *slot {
                built.resize(slot + 1, None);
            }
            built[*slot] = Some(m.clone());
        }
        dev.amg_select(0)
    }

    /// Assembles and solves the elliptic f equation −∇·(L²∇f) + f = R of
    /// k-ε-ζ-f (symmetric positive definite: CG, by default with the
    /// aggregation multigrid on its own hierarchy). Returns (scaled
    /// residual, iterations).
    fn elliptic(&self, dev: &mut dyn FlowDevice, z: ZetaFDev, tr: bool) -> Result<(f64, usize)> {
        let v = z.f;
        dev.gradient_with(v.lsq, v.phi, v.bc_a, v.bc_b, None, v.grad)?;
        dev.scalar_faces(&ScalarFaceArgs {
            bc_a: v.bc_a,
            bc_b: v.bc_b,
            jump: None,
            mass_flux: self.zero_face,
            gamma: z.gamma_f,
            phi: v.phi,
            grad: v.grad,
            psi: v.psi,
            blend: 0.0,
            central: 0.0,
            beta: 0.1,
            nonorth: 1.0,
            bound: 0.0,
            cp: self.cp,
            cn: self.cn,
            bf: self.bf,
        })?;
        dev.scalar_rows(&ScalarRowArgs {
            cp: self.cp,
            cn: self.cn,
            bf: self.bf,
            mass_flux: self.zero_face,
            phi: v.phi,
            su: Some(z.src[2]),
            sp: Some(z.src[3]),
            alpha: 1.0,
            time: None,
            fix: None,
            pseudo: None,
            capacity: None,
            positive: false,
            val: self.val,
            rhs: self.rhs,
            diag: self.diag,
            order: self.order,
        })?;
        let scale = dev.scalar_scale(v.phi, self.rhs, self.diag, 1.0, None)?;
        let e = &self.elliptic;
        let mut pc = e.preconditioner;
        if pc == Preconditioner::Amg && !self.amg_prepare(dev, 1, false)? {
            pc = Preconditioner::Jacobi;
        }
        let solved = dev.solve(self.val, self.rhs, v.phi, pc, true, e.tol(tr), e.max_iter);
        if pc == Preconditioner::Amg {
            dev.amg_select(0)?;
        }
        let (it, r0) = solved?;
        Ok((if scale > 0.0 { r0 / scale } else { r0 }, it))
    }

    /// Start of a time step (as the scalars: kⁿ⁻¹ ← kⁿ ← k, with the
    /// extrapolated start value of BDF2).
    pub fn step_start(
        &mut self,
        dev: &mut dyn FlowDevice,
        has_old: bool,
        extrapolate: Option<f64>,
    ) -> Result<()> {
        if let Some(c) = self.curv.as_mut() {
            // a time level of S_ij only from a computed S_ij (not in the
            // first step of a run or after a restart: no time part there)
            c.has_old = self.curv_fresh.0.get();
            if c.has_old {
                for i in 0..6 {
                    dev.copy(c.s[i], c.s_old[i])?;
                }
            }
        }
        for v in &self.vars() {
            let Some((o1, o2, keep)) = v.old else {
                continue;
            };
            dev.copy(o2, keep)?;
            if has_old {
                dev.copy(o1, o2)?;
            }
            dev.copy(v.phi, o1)?;
            if has_old && let Some(w) = extrapolate {
                dev.lincomb(v.phi, 1.0 + w, o1, -w, o2)?;
                // The extrapolated start may not fall below a fraction of
                // the last value (every model, k, ε or ω and ζ). Where a
                // field more than halves in a step (a run or a jet started
                // into gas at rest, an injection that starts or ends, a
                // valve, a schedule) the linear extrapolation turns
                // negative and the absolute floor would put ε or ω near
                // zero: ν_t jumps to its cap, and with k-ε-ζ-f T = k/ε,
                // L ∝ k^{3/2}/ε and P/ε in the source of f blow up before
                // the first solve of the step (f of ±1e7, ζ between its
                // bounds, divergence within ten steps). The bound of a
                // solve of k-ε-ζ-f; inactive where the fields change by
                // less than this factor per step.
                let rel = Some((o1, ZETA_F_DROP));
                dev.floor(v.phi, v.floor, f64::INFINITY, rel, self.work)?;
            }
        }
        Ok(())
    }

    /// Takes the last step back exactly.
    pub fn step_back(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for v in &self.vars() {
            let Some((o1, o2, keep)) = v.old else {
                continue;
            };
            dev.copy(o1, v.phi)?;
            dev.copy(o2, o1)?;
            dev.copy(keep, o2)?;
        }
        Ok(())
    }

    /// The old levels of the time terms of k, ε or ω and ζ with the data of
    /// the cells moved as the sliding mesh says (`shift`: by one step for
    /// the old level n, by two for n−1; D182, as the passive scalars): the time term of a rotor cell sees
    /// what that cell held. The start values of the step, which
    /// [`step_start`](Self::step_start) made from the levels at the places
    /// of the mesh, are made again from the moved levels in the same way
    /// (`has_old`, `extrapolate` as there), and f (no time term) starts from
    /// its moved value: the coefficients of the first outer iteration (T,
    /// L, ν_t of the sources) belong to the turbulence that the cell
    /// carries, not to the one that was at its place.
    pub fn shift_old(
        &mut self,
        dev: &mut dyn FlowDevice,
        shift: &crate::sliding::Shift,
        has_old: bool,
        extrapolate: Option<f64>,
    ) -> Result<()> {
        let nc = self.nc;
        let work = self.work;
        let moved = |dev: &mut dyn FlowDevice, src: Buf, dst: Buf, steps: usize| -> Result<()> {
            shift.moved(dev, src, dst, steps)
        };
        let shift = |dev: &mut dyn FlowDevice, v: &mut Var| -> Result<()> {
            let Some((o1, o2, _)) = v.old else {
                // (f: the iterate itself)
                return moved(dev, v.phi, v.phi, 1);
            };
            let (c1, c2) = match v.content {
                Some(c) => c,
                None => {
                    let c = (dev.alloc(nc)?, dev.alloc(nc)?);
                    v.content = Some(c);
                    c
                }
            };
            moved(dev, o1, c1, 1)?;
            moved(dev, o2, c2, 2)?;
            dev.copy(c1, v.phi)?;
            if has_old && let Some(w) = extrapolate {
                dev.lincomb(v.phi, 1.0 + w, c1, -w, c2)?;
                dev.floor(v.phi, v.floor, f64::INFINITY, None, work)?;
            }
            Ok(())
        };
        shift(dev, &mut self.k)?;
        shift(dev, &mut self.w)?;
        if let Some(z) = &mut self.zf {
            shift(dev, &mut z.z)?;
            shift(dev, &mut z.f)?;
        }
        Ok(())
    }

    /// Keeps a device copy of k and ω as the last valid state (D63).
    pub fn keep_valid(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        let nc = self.nc;
        let keep = |dev: &mut dyn FlowDevice, v: &mut Var| -> Result<()> {
            let b = match v.valid {
                Some(b) => b,
                None => dev.alloc(nc)?,
            };
            dev.copy(v.phi, b)?;
            v.valid = Some(b);
            Ok(())
        };
        keep(dev, &mut self.k)?;
        keep(dev, &mut self.w)?;
        if let Some(z) = &mut self.zf {
            keep(dev, &mut z.z)?;
            keep(dev, &mut z.f)?;
        }
        Ok(())
    }

    /// Puts the last valid k and ω back.
    pub fn rewind(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for v in &self.vars() {
            if let Some(b) = v.valid {
                dev.copy(b, v.phi)?;
            }
        }
        Ok(())
    }

    /// After a rewind of a time step: the old level of n−1 back into kⁿ.
    pub fn rewind_levels(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        for v in &self.vars() {
            if let Some((o1, o2, _)) = v.old {
                dev.copy(o2, o1)?;
            }
        }
        Ok(())
    }

    /// Host copy of the state; `p`, `pb` are the solver's p̂ (cells and
    /// boundary faces).
    /// k and the dissipation ε per cell (the eddy dissipation concept,
    /// SF3): ε of the k-ε family, β* k ω of SST.
    pub(crate) fn k_and_epsilon(&self, dev: &mut dyn FlowDevice) -> Result<(Vec<f64>, Vec<f64>)> {
        let k = dev.download(self.k.phi)?;
        let mut e = dev.download(self.w.phi)?;
        if !self.eps {
            let bs = self.section.sst_constants()?.beta_star;
            for (e, k) in e.iter_mut().zip(&k) {
                *e *= bs * k.max(0.0);
            }
        }
        Ok((k, e))
    }

    pub(crate) fn state(
        &self,
        dev: &mut dyn FlowDevice,
        p: &[f64],
        pb: &[f64],
    ) -> Result<TurbState> {
        let ni = self.ni;
        let k = dev.download(self.k.phi)?;
        let omega = dev.download(self.w.phi)?;
        let k_b = dev.download(self.k.bc_b)?;
        let omega_b = dev.download(self.w.bc_b)?;
        let mut u_tau = vec![0.0; self.nb];
        let mut y_plus = vec![0.0; self.nb];
        let mut wall_blending = vec![0.0; self.nb];
        if let Some(wd) = self.wd {
            let ut = dev.download(wd.u_tau)?;
            let yp = dev.download(wd.y_plus)?;
            let vw = dev.download(wd.viscous_weight)?;
            for (i, &f) in self.walls.stat.faces.iter().enumerate() {
                u_tau[f as usize - ni] = ut[i];
                y_plus[f as usize - ni] = yp[i];
                wall_blending[f as usize - ni] = vw[i];
            }
        }
        Ok(TurbState {
            epsilon: self.eps,
            nu_t: dev.download(self.nu_t)?,
            k_old: match self.k.old {
                Some((o1, _, _)) => Some(dev.download(o1)?),
                None => None,
            },
            omega_old: match self.w.old {
                Some((o1, _, _)) => Some(dev.download(o1)?),
                None => None,
            },
            mu_eff: dev.download(dev.face_viscosity())?,
            k_b,
            omega_b,
            u_tau,
            y_plus,
            wall_blending: Some(wall_blending),
            p_modified: p.to_vec(),
            pb_modified: pb.to_vec(),
            wall_distance: self.wall_distance.clone(),
            zeta: match self.zf {
                Some(z) => Some(ZetaState {
                    zeta: dev.download(z.z.phi)?,
                    f: dev.download(z.f.phi)?,
                    zeta_old: match z.z.old {
                        Some((o1, _, _)) => Some(dev.download(o1)?),
                        None => None,
                    },
                    zeta_b: dev.download(z.z.bc_b)?,
                    f_b: dev.download(z.f.bc_b)?,
                    t: dev.download(z.t)?,
                }),
                None => None,
            },
            amg: self
                .amg_built
                .borrow()
                .iter()
                .enumerate()
                .filter_map(|(s, m)| m.as_ref().map(|m| (s, m.clone())))
                .collect(),
            k,
            omega,
            fr1: match self.curv {
                Some(c) => Some(dev.download(c.fr1)?),
                None => None,
            },
        })
    }

    /// k on the boundary faces, a k_P + b.
    pub(crate) fn k_face(&self, k: &[f64], k_b: &[f64], owner: &[u32]) -> Vec<f64> {
        (0..self.nb)
            .map(|fb| self.ka_host[fb] * k[owner[self.ni + fb] as usize] + k_b[fb])
            .collect()
    }

    /// Loads a host state (a restart). A state of the other family (ω for a
    /// k-ε case, or ε for SST) is converted with ε = C_μ k ω.
    pub(crate) fn load(&mut self, dev: &mut dyn FlowDevice, s: &TurbState) -> Result<()> {
        self.curv_fresh.0.set(false);
        let converted;
        let s = if s.epsilon == self.eps {
            s
        } else {
            let cm = self.c_mu;
            let to = |k: &[f64], v: &[f64]| -> Vec<f64> {
                k.iter()
                    .zip(v)
                    .map(|(&k, &v)| {
                        if self.eps {
                            cm * k * v
                        } else {
                            v / (cm * k).max(1e-300)
                        }
                    })
                    .collect()
            };
            let mut c = s.clone();
            c.epsilon = self.eps;
            c.omega = to(&s.k, &s.omega);
            if let (Some(ko), Some(wo)) = (&s.k_old, &s.omega_old) {
                c.omega_old = Some(to(ko, wo));
            }
            c.omega_b = to(&s.k_b, &s.omega_b);
            converted = c;
            &converted
        };
        dev.upload(self.k.phi, &s.k)?;
        dev.upload(self.w.phi, &s.omega)?;
        dev.upload(self.nu_t, &s.nu_t)?;
        if let (Some((o1, _, _)), Some(o)) = (self.k.old, &s.k_old) {
            dev.upload(o1, o)?;
        }
        if let (Some((o1, _, _)), Some(o)) = (self.w.old, &s.omega_old) {
            dev.upload(o1, o)?;
        }
        dev.upload(dev.face_viscosity(), &s.mu_eff)?;
        // inlet and wall parts of the boundary values
        let mut kin = vec![0.0; self.nb];
        let mut win = vec![0.0; self.nb];
        let mut kw = vec![0.0; self.nb];
        let mut ww = vec![0.0; self.nb];
        for &(fb, _, _) in &self.inlets.faces {
            let fb = fb as usize;
            kin[fb] = s.k_b[fb];
            win[fb] = s.omega_b[fb];
        }
        for fb in 0..self.nb {
            if self.walls.mask_b[fb] > 0.0 {
                kw[fb] = s.k_b[fb];
                ww[fb] = s.omega_b[fb];
            }
        }
        dev.upload(self.k.b_in, &kin)?;
        dev.upload(self.w.b_in, &win)?;
        dev.upload(self.k.b_wall, &kw)?;
        dev.upload(self.w.b_wall, &ww)?;
        self.kin_host = kin;
        self.win_host = win;
        if let Some(z) = self.zf {
            match &s.zeta {
                Some(zs) => {
                    dev.upload(z.z.phi, &zs.zeta)?;
                    dev.upload(z.f.phi, &zs.f)?;
                    if let (Some((o1, _, _)), Some(o)) = (z.z.old, &zs.zeta_old) {
                        dev.upload(o1, o)?;
                    }
                    let mut zin = vec![0.0; self.nb];
                    for &(fb, _, _) in &self.inlets.faces {
                        zin[fb as usize] = zs.zeta_b[fb as usize];
                    }
                    let fw: Vec<f64> = (self.walls.mask_b.iter().zip(&zs.f_b))
                        .map(|(&m, &f)| if m > 0.0 { f } else { 0.0 })
                        .collect();
                    dev.upload(z.z.b_in, &zin)?;
                    dev.upload(z.f.b_wall, &fw)?;
                    self.zin_host = zin;
                }
                None => {
                    // a restart of a two-equation model: ζ = 2/3, f = 0
                    // (α = 1)
                    dev.fill(z.z.phi, 2.0 / 3.0)?;
                    dev.fill(z.f.phi, if self.model.is_eb() { 1.0 } else { 0.0 })?;
                    if let Some((o1, _, _)) = z.z.old {
                        dev.fill(o1, 2.0 / 3.0)?;
                    }
                }
            }
        }
        // the backflow faces of the outlets: the stored boundary values
        // tell which faces took entering fluid (k_b > 0 there, 0 where the
        // flow left)
        if !self.backflow.faces.is_empty() {
            let zb = s.zeta.as_ref().map(|z| &z.zeta_b);
            let st: Vec<_> = self
                .backflow
                .faces
                .iter()
                .map(|&(fb, _, _)| {
                    let fb = fb as usize;
                    let v = (s.k_b[fb] > 0.0)
                        .then(|| (s.k_b[fb], s.omega_b[fb], zb.map_or(2.0 / 3.0, |z| z[fb])));
                    (fb, v)
                })
                .collect();
            self.apply_backflow(dev, &st)?;
        }
        self.amg_restore(dev, &s.amg)?;
        self.boundary_values(dev)?;
        if let Some(wd) = self.wd {
            let ut: Vec<f64> = self
                .walls
                .stat
                .faces
                .iter()
                .map(|&f| s.u_tau[f as usize - self.ni])
                .collect();
            dev.upload(wd.u_tau, &ut)?;
        }
        Ok(())
    }

    /// The modified pressures to load into the solver.
    pub(crate) fn pressure_of(s: &TurbState) -> (&[f64], &[f64]) {
        (&s.p_modified, &s.pb_modified)
    }

    /// Per wall boundary: the area-weighted min, mean and max of y⁺ and
    /// the area fraction outside the claimed band (low-Re: above 2;
    /// automatic: above 300), for the report and TURB001.
    pub fn y_plus_summary(&self, dev: &mut dyn FlowDevice, mesh: &Mesh) -> Result<Vec<YPlusRow>> {
        let Some(wd) = self.wd else {
            return Ok(Vec::new());
        };
        let yp = dev.download(wd.y_plus)?;
        Ok(self.y_plus_rows(&self.walls, &yp, mesh))
    }

    /// [`Self::y_plus_summary`] of a run on ranks (MPI step W7): the y⁺ of
    /// the wall faces of this rank gathered into the boundary faces of the
    /// whole mesh `mesh` and summed over its walls in the order of one rank.
    /// Collective.
    fn y_plus_summary_on(
        &self,
        dev: &mut dyn FlowDevice,
        mesh: &Mesh,
        geo: &Geometry,
        bcs: &BTreeMap<String, BoundaryCondition>,
        g: &mut crate::ranks::Fields<'_>,
    ) -> Result<Vec<YPlusRow>> {
        use crate::restart::FieldMap as _;
        let Some(wd) = self.wd else {
            return Ok(Vec::new());
        };
        let yp = dev.download(wd.y_plus)?;
        let ni = self.ni;
        let mut b = vec![0.0; self.nb];
        for (i, &f) in self.walls.stat.faces.iter().enumerate() {
            b[f as usize - ni] = yp[i];
        }
        let b = g.boundary(&b)?;
        let walls = Walls::new(mesh, geo, bcs, self.eps);
        let gi = mesh.n_interior_faces();
        let yp: Vec<f64> = walls
            .stat
            .faces
            .iter()
            .map(|&f| b[f as usize - gi])
            .collect();
        Ok(self.y_plus_rows(&walls, &yp, mesh))
    }

    /// Per wall boundary of `walls`: min, area mean and max of `yp` (one
    /// value per wall face) and the share of the area outside the band of
    /// the wall treatment.
    fn y_plus_rows(&self, walls: &Walls, yp: &[f64], mesh: &Mesh) -> Vec<YPlusRow> {
        let mut per: BTreeMap<usize, (f64, f64, f64, f64, f64)> = BTreeMap::new();
        for (i, &b) in walls.boundary.iter().enumerate() {
            let a = walls.area[i];
            let e = per
                .entry(b)
                .or_insert((f64::INFINITY, 0.0, f64::NEG_INFINITY, 0.0, 0.0));
            e.0 = e.0.min(yp[i]);
            e.1 += a * yp[i];
            e.2 = e.2.max(yp[i]);
            e.3 += a;
            // area outside the claimed band (automatic: above the log
            // layer, y⁺ > 300; validated on the flat plate from 0.24 to 85)
            let out = match self.treatment {
                kwall::LOW_RE => yp[i] > 2.0,
                // ζ-f with the compound treatment: any y⁺ (validated from
                // 0.5 to 100 on the flat plate)
                kwall::COMPOUND | kwall::COMPOUND_2023 if self.zf.is_some() => false,
                kwall::AUTOMATIC => yp[i] > 300.0,
                kwall::SCALABLE => yp[i] < 30.0,
                kwall::STANDARD => !(30.0..=300.0).contains(&yp[i]),
                // compound with the high-Re k-ε: below y⁺ 30 the cells of
                // the buffer layer have no near-wall damping (6c finding)
                _ => !(30.0..=300.0).contains(&yp[i]),
            };
            if out {
                e.4 += a;
            }
        }
        per.into_iter()
            .map(|(b, e)| {
                (
                    mesh.boundaries[b].name.clone(),
                    e.0,
                    e.1 / e.3,
                    e.2,
                    e.4 / e.3,
                )
            })
            .collect()
    }

    /// The model's name for the log.
    pub fn label(&self) -> String {
        if self.zf.is_some() {
            return format!(
                "{}, {} wall treatment",
                self.model.name(),
                match self.treatment {
                    kwall::LOW_RE => "low-Re",
                    kwall::COMPOUND_2023 => "compound (2023)",
                    _ => "compound",
                }
            );
        }
        if self.eps {
            return format!(
                "{} k-epsilon, {} wall treatment{}",
                if self.model == TurbulenceModel::RealizableKEpsilon {
                    "realizable"
                } else {
                    "standard"
                },
                match self.treatment {
                    kwall::STANDARD => "standard",
                    kwall::COMPOUND => "compound",
                    _ => "scalable",
                },
                match self.ko.limiter {
                    kkeps::LIMITER_KATO_LAUNDER => ", Kato-Launder production",
                    kkeps::LIMITER_DURBIN => ", Durbin-bounded",
                    kkeps::LIMITER_CLIP => ", production clipped at 10 epsilon",
                    _ => "",
                }
            );
        }
        format!(
            "k-omega SST (2003), {} wall treatment, {} production",
            match self.treatment {
                kwall::LOW_RE => "low-Re",
                _ => "automatic",
            },
            match self.o.limiter {
                ksst::LIMITER_KATO_LAUNDER => "Kato-Launder",
                ksst::LIMITER_DURBIN => "Durbin-bounded",
                ksst::LIMITER_NONE => "unlimited",
                _ => "Menter-limited",
            }
        )
    }

    /// Number of cells.
    pub fn n_cells(&self) -> usize {
        self.nc
    }

    /// The ν_t cap that acted in the last iteration, ν_t/ν max and the
    /// cell (for the extremes block).
    pub fn ratio_max(&self, dev: &mut dyn FlowDevice) -> Result<(f64, usize)> {
        self.ratio_max_on(dev, None)
    }

    /// [`Self::ratio_max`]; on ranks (`g`, MPI step W7) over ν_t gathered
    /// into the cells of the whole mesh, with the cell index there.
    /// Collective on ranks.
    fn ratio_max_on(
        &self,
        dev: &mut dyn FlowDevice,
        g: Option<&mut crate::ranks::Fields<'_>>,
    ) -> Result<(f64, usize)> {
        let mut v = dev.download(self.nu_t)?;
        if let Some(g) = g {
            use crate::restart::FieldMap as _;
            anyhow::ensure!(
                self.fluids.is_none(),
                "the turbulence summary of several fluids does not run on ranks yet (MPI step W7)"
            );
            v = g.cells(&v)?;
        }
        let (mut m, mut at) = (0.0, 0);
        for (c, &x) in v.iter().enumerate() {
            if x > m {
                m = x;
                at = c;
            }
        }
        Ok((m / self.nu_c(at), at))
    }
}

/// The first-cell y⁺ of one wall boundary.
#[derive(Debug, Clone, Serialize)]
pub struct WallYPlus {
    /// Boundary.
    pub boundary: String,
    /// Smallest.
    pub min: f64,
    /// Area-weighted mean.
    pub mean: f64,
    /// Largest.
    pub max: f64,
}

/// The turbulence model in the report of a solve.
#[derive(Debug, Clone, Serialize)]
pub struct TurbulenceSummary {
    /// Model and options.
    pub model: String,
    /// The second variable, `omega` or `epsilon`.
    pub second: &'static str,
    /// First-cell y⁺ per wall.
    pub walls: Vec<WallYPlus>,
    /// Largest ν_t/ν and its cell.
    pub viscosity_ratio_max: f64,
    /// The cell of `viscosity_ratio_max`.
    pub viscosity_ratio_cell: usize,
    /// Floor and cap counts of the run.
    pub health: TurbHealth,
}

/// The summary of the turbulence model of `fl`, the findings it gives
/// (TURB001) and the log lines.
pub fn summary(
    fl: &mut crate::simplec::Flow,
    mesh: &Mesh,
) -> Result<Option<(TurbulenceSummary, Vec<Diagnostic>)>> {
    summary_on(fl, mesh, None)
}

/// [`summary`]; on ranks (MPI step W7) `mesh` is the whole mesh and the
/// wall and cell values are gathered into its order, with its geometry and
/// the boundary conditions of the case. Collective on ranks.
pub(crate) fn summary_on(
    fl: &mut crate::simplec::Flow,
    mesh: &Mesh,
    on: Option<(
        &crate::ranks::OnRanks,
        &Geometry,
        &BTreeMap<String, BoundaryCondition>,
    )>,
) -> Result<Option<(TurbulenceSummary, Vec<Diagnostic>)>> {
    let Some(t) = fl.turbulence.as_ref() else {
        return Ok(None);
    };
    let dev = fl.dev.as_mut();
    let (walls, (rmax, cell)) = match on {
        Some((o, geo, bcs)) => {
            let mut g = o.fields();
            (
                t.y_plus_summary_on(dev, mesh, geo, bcs, &mut g)?,
                t.ratio_max_on(dev, Some(&mut g))?,
            )
        }
        None => (t.y_plus_summary(dev, mesh)?, t.ratio_max(dev)?),
    };
    let mut finds = Vec::new();
    for (e, name) in equations(t.model).iter().enumerate() {
        let n = t.health.floor_iterations_late[e];
        if n > 0 {
            finds.push(Diagnostic::new(
                diag::TURB_FLOOR,
                format!(
                    "the positivity floor of {name} acted on more than 0.1 % of the cells in \
                     {n} of the final tenth of the solves, at most {} cells in one (a \
                     converging run does not need it at its end: check the mesh near the \
                     boundary-layer edge and the inlet values)",
                    t.health.floor_cells_max[e]
                ),
            ));
        }
    }
    for (b, _, _, _, out) in &walls {
        if *out > 0.1 {
            let text = if t.eps && t.treatment != kwall::LOW_RE {
                let (band, advice) = match t.treatment {
                    kwall::SCALABLE => (
                        "below 30",
                        "the scalable wall functions take it without failure, but their \
                         accuracy is claimed for y⁺ ≥ 30 only (a flat plate at y⁺ 1 to 10 \
                         stays within 3 %); resolved walls need ζ-f or EB",
                    ),
                    kwall::STANDARD => (
                        "outside 30 … 300",
                        "the standard wall functions need the first cell in the log layer; \
                         the scalable treatment takes any y⁺ without failure",
                    ),
                    _ => (
                        "outside 30 … 300",
                        "with the high-Re k-ε models the compound treatment is claimed for \
                         30 ≤ y⁺ ≤ 300: below, the cells of the buffer layer have no \
                         near-wall damping and the wall shear stress is overpredicted (a \
                         flat plate at y⁺ 10 by 21 %, at y⁺ 1 by 58 %); resolved walls need \
                         ζ-f or EB",
                    ),
                };
                format!(
                    "wall `{b}`: first-cell y⁺ {band} on {:.0} % of the area: {advice}",
                    100.0 * out
                )
            } else if t.treatment == kwall::LOW_RE {
                format!(
                    "wall `{b}`: first-cell y⁺ above 2 on {:.0} % of the area with the low-Re \
                     treatment (it needs y⁺ ≲ 1; the {} treatment takes any y⁺)",
                    100.0 * out,
                    if t.model.is_eb() {
                        "compound"
                    } else if t.zf.is_some() {
                        "compound-2023 or compound"
                    } else {
                        "automatic"
                    }
                )
            } else {
                format!(
                    "wall `{b}`: first-cell y⁺ above 300 on {:.0} % of the area, above the log \
                     layer that the automatic treatment assumes (validated for y⁺ 0.24 to 85 on \
                     the flat plate, gallery 68): refine the wall cells",
                    100.0 * out
                )
            };
            finds.push(Diagnostic::new(diag::TURB_Y_PLUS, text));
        }
    }
    Ok(Some((
        TurbulenceSummary {
            model: t.label(),
            second: equations(t.model)[1],
            walls: walls
                .into_iter()
                .map(|(boundary, min, mean, max, _)| WallYPlus {
                    boundary,
                    min,
                    mean,
                    max,
                })
                .collect(),
            viscosity_ratio_max: rmax,
            viscosity_ratio_cell: cell,
            health: t.health.clone(),
        },
        finds,
    )))
}

/// Adds the turbulence model of `case` to the flow solver `fl` (no-op for a
/// laminar case); call after [`crate::scalar::attach`] and before loading a
/// restart.
pub fn attach(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
) -> Result<()> {
    attach_with(fl, case, mesh, geo, None)
}

/// [`attach`] with the wall distance per cell given (`None`: the model
/// computes it as `turbulence.wall_distance` says).
pub fn attach_with(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    wall_distance: Option<Vec<f64>>,
) -> Result<()> {
    if case.turbulence.is_laminar() {
        return Ok(());
    }
    let fv = fl.dev.face_viscosity();
    // several fluids (D194): ν and ρ per cell
    let cells = if crate::fluids::zone_fluids(case, mesh).1.len() > 1 {
        let mu = case.properties()?.cell_viscosity(mesh, geo);
        let rho: Vec<f64> = (0..mesh.n_cells).map(|c| fl.phys.rho_at(c)).collect();
        let nu: Vec<f64> = mu.iter().zip(&rho).map(|(m, r)| m / r).collect();
        Some((nu, rho))
    } else {
        None
    };
    let t = Turbulence::new(
        fl.dev.as_mut(),
        case,
        mesh,
        geo,
        &fl.disc,
        &fl.bc,
        fl.phys.rho,
        fl.phys.mu,
        fv,
        wall_distance,
        cells,
    )?;
    fl.add_turbulence(t)
}
