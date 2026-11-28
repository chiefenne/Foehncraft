//! Two fluids by the Euler–Euler two-fluid model (L3, stage MP3, D352 to
//! D355): dispersed gas bubbles in a liquid, one momentum equation per
//! phase and one pressure.
//!
//! **Equations** (D352). Per phase k (l liquid, g gas), with p the pressure
//! less the hydrostatic part of the liquid, p = p_true − ρ_l g·x:
//!
//! α_k ρ_k (∂u_k/∂t + u_k·∇u_k) = −α_k ∇p + ∇·(α_k μ_k,eff ∇u_k)
//!     + α_k (ρ_k − ρ_l) g + M_k,
//!
//! ∂α_g/∂t + ∇·(α_g u_g) = 0, α_l = 1 − α_g, and the volume continuity
//! ∇·(α_l u_l + α_g u_g) = 0 for the pressure. The interphase forces M_g =
//! −M_l (drag, lift, wall lubrication, turbulent dispersion) are those of
//! [`foehncraft_kernels::twofluid`] (D353).
//!
//! **Discretisation.** One device holds both phases: the liquid in the
//! fields of the [`Flow`], the gas in a second set of the same buffers
//! that is swapped in for its operations, and the material arrays of a
//! phase (α_k ρ_k per cell, ρ_k per face, α_k,f μ_k,eff per face, the cell
//! force b_k and the face offsets) are copied into the device's material
//! buffers before them. So the momentum rows, the Rhie–Chow flux and the
//! pressure coefficients of the single fluid serve both phases, and a case
//! without `[two_fluid]` runs bit for bit as before.
//!
//! - The momentum equations are in the non-conservative form: the time
//!   term α_kⁿ⁺¹ ρ_k (u − uⁿ)/Δt and the convection with the phase mass
//!   flux α_k,f ṁ_k less u_P Σ_f α_k,f ṁ_k (the continuity error of the
//!   phase, an implicit sink), so a cell's velocity does not take the
//!   change of its phase fraction.
//! - The pressure term α_k ∇p and the cell force α_k b_k enter the rows;
//!   D of the phase is α_k V/a_P, so the Rhie–Chow flux of the phase, the
//!   pressure coefficient and the velocity correction carry α_k.
//! - Drag is implicit in the rows of both phases (K V on the diagonal, K V
//!   u of the other phase on the right, the latest one; the gas first).
//! - The pressure equation is the volume continuity Σ_k α_k,f ṁ_k/ρ_k = 0
//!   with the coefficients Σ_k α_k,f c_k,f/ρ_k; both fluxes and both
//!   velocities take the correction.
//! - α_g advances explicitly at the start of a step with the gas volume
//!   flux of the step before, first-order upwind with sub-steps of a
//!   bounded Courant number (D354); the face fractions α_k,f of the step
//!   are upwind of the phase fluxes.
//! - Turbulence (D355): the liquid's model on the liquid velocity and mass
//!   flux (ρ_l), the gas viscosity μ_g + μ_t ρ_g/ρ_l, the liquid's
//!   viscosity plus the bubble-induced μ_BIT of Sato.
//!
//! The time step is backward Euler; the outer iterations of
//! `[solver] outer_iterations` couple the phases.

use crate::device::{Buf, Buf3, EliminateBufs, Lsq, Mask, MaterialBufs, TwoFluidBufs, TwoFluidOp, VofBufs, VofOp};
use crate::diag::{self, Diagnostic};
use crate::simplec::{Fields, Flow, Residuals, Work};
use anyhow::{Context, Result};
use foehncraft_core::exact::ExactSumIter;
use foehncraft_kernels::twofluid as k;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::time::Instant;

/// The lift model.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Lift {
    /// Tomiyama et al. (2002).
    #[default]
    Tomiyama,
    /// The constant `lift_coefficient`.
    Constant,
    /// No lift.
    None,
}

/// The wall lubrication model.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum WallLubrication {
    /// Hosokawa et al. (2002), in the form of Rzehak et al. (2012).
    #[default]
    Hosokawa,
    /// None.
    None,
}

/// `[two_fluid.gas]`: the dispersed phase.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GasSection {
    /// Name (output).
    #[serde(default = "gas_name")]
    pub name: String,
    /// Density ρ_g (kg/m³).
    pub density: f64,
    /// Dynamic viscosity μ_g (Pa s).
    pub viscosity: f64,
    /// Bubble diameter d (m), one size.
    pub diameter: f64,
}

fn gas_name() -> String {
    "gas".into()
}

/// One `[two_fluid.inflow.<boundary>]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InflowSection {
    /// α_g of the inflow.
    pub alpha: f64,
    /// Gas velocity of the inflow (m/s); default: the liquid's plus the
    /// terminal rise velocity against gravity.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub velocity: Option<[f64; 3]>,
}

/// `[two_fluid]`: gas bubbles dispersed in the liquid of `[fluid]`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TwoFluidSection {
    /// The gas.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gas: Option<GasSection>,
    /// Surface tension σ (N/m) of the Eötvös numbers.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub surface_tension: Option<f64>,
    /// α_g at the start (default 0).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial: Option<f64>,
    /// The inflow of gas per boundary (default: none).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub inflow: BTreeMap<String, InflowSection>,
    /// Lift (default Tomiyama).
    #[serde(default)]
    pub lift: Lift,
    /// C_L of `lift = "constant"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lift_coefficient: Option<f64>,
    /// Damp the lift near walls after Shaver & Podowski (2015): none
    /// closer than a bubble radius, full beyond a diameter (default true).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lift_wall_damping: Option<bool>,
    /// Wall lubrication (default Hosokawa).
    #[serde(default)]
    pub wall_lubrication: WallLubrication,
    /// C_TD of the turbulent dispersion of Burns et al. (default 1; 0
    /// none; needs a turbulence model).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub turbulent_dispersion: Option<f64>,
    /// C_μb of the bubble-induced viscosity of Sato et al. (default 0.6
    /// with a turbulence model, 0 without).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bubble_viscosity: Option<f64>,
    /// The largest outflow Courant number of a cell in one α sub-step
    /// (default 0.5).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_courant: Option<f64>,
}

impl TwoFluidSection {
    /// True without `[two_fluid]`.
    pub fn is_empty(&self) -> bool {
        self.gas.is_none()
    }

    fn courant(&self) -> f64 {
        self.max_courant.unwrap_or(0.5)
    }
}

fn bad(key: &str, msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::TWO_FLUID_CONFIG, msg).key(key.to_string()))
}

/// The ranges of `[two_fluid]` and the settings it does not work with yet
/// (TWF001, TWF002).
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let t = &case.two_fluid;
    let Some(g) = &t.gas else {
        return Ok(());
    };
    let pos = |x: f64| x > 0.0 && x.is_finite();
    if !(pos(g.density) && pos(g.viscosity) && pos(g.diameter)) {
        return Err(bad(
            "two_fluid.gas",
            "[two_fluid.gas] density, viscosity and diameter must be > 0".into(),
        ));
    }
    if g.density >= case.rho() {
        return Err(bad(
            "two_fluid.gas.density",
            format!(
                "[two_fluid.gas] density {} must be below that of the liquid in [fluid] ({})",
                g.density,
                case.rho()
            ),
        ));
    }
    let s = t.surface_tension.unwrap_or(0.0);
    if !(s >= 0.0 && s.is_finite()) {
        return Err(bad(
            "two_fluid.surface_tension",
            format!("[two_fluid] surface_tension {s} must be >= 0"),
        ));
    }
    let a0 = t.initial.unwrap_or(0.0);
    if !(0.0..1.0).contains(&a0) {
        return Err(bad(
            "two_fluid.initial",
            format!("[two_fluid] initial {a0} must be in [0, 1)"),
        ));
    }
    for (name, i) in &t.inflow {
        if !(0.0..1.0).contains(&i.alpha) {
            return Err(bad(
                "two_fluid.inflow",
                format!("[two_fluid.inflow.{name}] alpha {} must be in [0, 1)", i.alpha),
            ));
        }
        if !case.boundary.contains_key(name) {
            return Err(bad(
                "two_fluid.inflow",
                format!("[two_fluid.inflow] names the boundary `{name}`, which the case does not have"),
            ));
        }
    }
    if t.lift == Lift::Constant && t.lift_coefficient.is_none() {
        return Err(bad(
            "two_fluid.lift_coefficient",
            "[two_fluid] lift = \"constant\" needs lift_coefficient".into(),
        ));
    }
    let c = t.courant();
    if !(c > 0.0 && c <= 1.0) {
        return Err(bad(
            "two_fluid.max_courant",
            format!("[two_fluid] max_courant {c} must be in (0, 1]"),
        ));
    }
    for (key, v) in [
        ("turbulent_dispersion", t.turbulent_dispersion),
        ("bubble_viscosity", t.bubble_viscosity),
    ] {
        if let Some(v) = v
            && !(v >= 0.0 && v.is_finite())
        {
            return Err(bad(
                &format!("two_fluid.{key}"),
                format!("[two_fluid] {key} {v} must be >= 0"),
            ));
        }
    }
    let mut no = Vec::new();
    match &case.time {
        None => no.push("a steady case (the two-fluid model needs [time])"),
        Some(tm) => {
            if tm.scheme != crate::transient::TimeScheme::Bdf1 {
                no.push("a time scheme other than bdf1 ([time] scheme = \"bdf1\")");
            }
            if tm.adaptive || tm.error_tolerance.is_some() {
                no.push("an adaptive time step");
            }
        }
    }
    if !case.vof.is_empty() {
        no.push("[vof]");
    }
    if case.energy.is_some() {
        no.push("[energy] (wall boiling, MP4, adds it)");
    }
    if !case.scalars.is_empty() {
        no.push("passive scalars or species");
    }
    if case.buoyancy.model != crate::buoyancy::Model::None {
        no.push("a [buoyancy] model (gravity acts through [gravity] alone)");
    }
    if !case.regions.is_empty() {
        no.push("regions or several fluids");
    }
    if !case.spray.is_empty() {
        no.push("sprays");
    }
    if !case.immersed.is_empty() {
        no.push("immersed bodies");
    }
    if case.solver.algorithm == Some(crate::simplec::Algorithm::Coupled) {
        no.push("the coupled algorithm");
    }
    if !case.mrf.is_empty() || !case.sliding.is_empty() {
        no.push("rotating zones or sliding interfaces");
    }
    if case.mesh_motion.is_some() {
        no.push("a moving or deforming mesh");
    }
    if case.run.restart.is_some() {
        no.push("a restart");
    }
    if case.boundary.values().any(|b| b.is_periodic()) {
        no.push("periodic boundaries");
    }
    if case.boundary.values().any(|b| b.is_convective_outlet()) {
        no.push("convective outlets");
    }
    if case.uniform_body_force() != Some([0.0; 3]) {
        no.push("[flow] body_force (gravity acts through [gravity])");
    }
    if !no.is_empty() {
        return Err(diag::err(Diagnostic::new(
            diag::TWO_FLUID_UNSUPPORTED,
            format!(
                "[two_fluid] does not work with {} in this version",
                no.join(", ")
            ),
        )));
    }
    Ok(())
}

/// The material arrays of one phase (copied into the device's material
/// buffers before its operations).
#[derive(Debug, Clone, Copy)]
struct Material {
    rho: Buf,
    rho_face: Buf,
    mu: Buf,
    body: Buf3,
    jump: Buf,
}

/// What one phase keeps besides its fields.
#[derive(Debug, Clone, Copy)]
struct Phase {
    mat: Material,
    /// α_k per cell (clipped) and per face.
    alpha: Buf,
    alpha_f: Buf,
    /// The convecting flux α_k,f ṁ_k.
    conv: Buf,
    /// The implicit sink K V − Σ_f α_k,f ṁ_k per cell.
    sink: Buf,
    /// The explicit interphase source per cell (N).
    src: Buf3,
    /// α_k (∇p − b_k) of the rows.
    gp: Buf3,
    /// uⁿ.
    u_old: Buf3,
    /// The old-level Rhie–Chow difference ṁ − ρ ū·S.
    x_old: Buf,
    /// c_f of the pressure correction.
    cf: Buf,
    /// ρ_k.
    rho: f64,
}

/// The totals of a run (for `result.json` and `two_fluid.csv`).
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema, PartialEq)]
pub struct Totals {
    /// Gas volume at the start (m³).
    pub volume_start: f64,
    /// Gas volume now (m³).
    pub volume: f64,
    /// Gas volume that left through the boundaries (m³, net).
    pub volume_out: f64,
    /// Gas volume clipped to keep 0 ≤ α ≤ 1 (m³, signed sum).
    pub clipped_volume: f64,
    /// Sub-steps of α in the last step.
    pub substeps: usize,
    /// Smallest and largest α_g after the last step.
    pub alpha_range: [f64; 2],
    /// The largest momentum residual of the gas in the last outer
    /// iteration.
    pub gas_residual: f64,
}

/// The two-fluid state of a flow.
#[derive(Debug)]
pub struct TwoFluid {
    /// The settings.
    pub section: TwoFluidSection,
    /// The constants of the kernels.
    pub k: k::Constants<f64>,
    gas: Phase,
    liq: Phase,
    /// The gas fields (swapped with the flow's for the gas operations).
    gas_f: Fields,
    gas_w: Work,
    /// α_g of the inflow per boundary face, and 1 − that.
    alpha_b: Buf,
    alpha_lb: Buf,
    alpha_new: Buf,
    flux: Buf,
    clipped: Buf,
    work: Buf,
    grad_alpha: Buf3,
    gnorm: Buf,
    wall_y: Buf,
    /// The cell size h (the dispersion limiter, D354).
    cell_h: Buf,
    wall_n: Buf3,
    kv: Buf,
    fv: Buf3,
    mu_bit: Buf,
    disp: Buf,
    /// μ_l + μ_t per face (the turbulence model's face viscosity).
    mu_eff: Buf,
    /// The material of the turbulence model: ρ_l everywhere.
    turb: Material,
    tmp: Buf,
    tmp3: Buf3,
    /// The liquid velocity of the gas rows (the elimination).
    u_lag: Buf3,
    cv: Buf,
    mv: Buf,
    volume: Vec<f64>,
    /// The totals.
    pub totals: Totals,
    /// α_g per cell of the end of the last step (host).
    pub alpha_host: Vec<f64>,
    /// The gas velocity of the end of the last step (host).
    pub gas_velocity: [Vec<f64>; 3],
}

/// The terminal rise speed (m/s) of a bubble of the constants `c` in still
/// liquid: K(U) U = (ρ_l − ρ_g) |g| α, by fixed-point iteration.
pub fn terminal_speed(c: &k::Constants<f64>) -> f64 {
    let g = (c.g[0] * c.g[0] + c.g[1] * c.g[1] + c.g[2] * c.g[2]).sqrt();
    if g == 0.0 {
        return 0.0;
    }
    let a = 1e-3;
    let mut u: f64 = 0.1;
    for _ in 0..200 {
        let kk = k::drag_k(c, a, u);
        let n = a * (c.rho_l - c.rho_g) * g / kk;
        if (n - u).abs() < 1e-12 * u.max(1e-12) {
            u = n;
            break;
        }
        u = 0.5 * (u + n);
    }
    u
}

impl TwoFluid {
    /// The constants of `case`.
    pub fn constants(case: &crate::case::FlowCase) -> Result<k::Constants<f64>> {
        let t = &case.two_fluid;
        let g = t.gas.as_ref().context("[two_fluid.gas]")?;
        let turbulent = case.turbulence.model != crate::turbulence::TurbulenceModel::Laminar;
        Ok(k::Constants {
            rho_l: case.rho(),
            rho_g: g.density,
            mu_l: case.mu(),
            mu_g: g.viscosity,
            sigma: t.surface_tension.unwrap_or(0.0),
            d: g.diameter,
            g: case.gravity.as_ref().map(|g| g.vector).unwrap_or([0.0; 3]),
            lift: match t.lift {
                Lift::Tomiyama => k::LIFT_TOMIYAMA,
                Lift::Constant => k::LIFT_CONSTANT,
                Lift::None => k::LIFT_NONE,
            },
            c_lift: t.lift_coefficient.unwrap_or(0.0),
            lift_damping: u32::from(t.lift_wall_damping.unwrap_or(true)),
            wall: match t.wall_lubrication {
                WallLubrication::Hosokawa => k::WALL_HOSOKAWA,
                WallLubrication::None => k::WALL_NONE,
            },
            c_td: if turbulent {
                t.turbulent_dispersion.unwrap_or(1.0)
            } else {
                0.0
            },
            sato: t
                .bubble_viscosity
                .unwrap_or(if turbulent { 0.6 } else { 0.0 }),
            alpha_min: 1e-6,
        })
    }

    /// The gas fraction of the cells and the material arrays of both phases
    /// from the current α_g (device).
    fn refresh_fractions(&self, fl: &mut Flow) -> Result<()> {
        let dev = fl.dev.as_mut();
        // α_l = 1 − α_g
        dev.lincomb(self.liq.alpha, -1.0, self.gas.alpha, 0.0, self.gas.alpha)?;
        dev.subtract(self.liq.alpha, -1.0)?;
        // α_k ρ_k with the residual fraction
        for ph in [&self.gas, &self.liq] {
            dev.copy(ph.alpha, ph.mat.rho)?;
            dev.floor(ph.mat.rho, self.k.alpha_min, 1.0, None, self.work)?;
            dev.lincomb(ph.mat.rho, ph.rho, ph.mat.rho, 0.0, ph.mat.rho)?;
        }
        Ok(())
    }

    /// Copies the material of `m` into the device's buffers.
    fn select(fl: &mut Flow, m: &Material, d: &MaterialBufs) -> Result<()> {
        let dev = fl.dev.as_mut();
        dev.copy(m.rho, d.rho)?;
        dev.copy(m.rho_face, d.rho_face)?;
        dev.copy(m.mu, d.mu)?;
        if let Some(b) = d.body {
            for (s, t) in m.body.all().into_iter().zip(b.all()) {
                dev.copy(s, t)?;
            }
        }
        if let Some(j) = d.jump {
            dev.copy(m.jump, j)?;
        }
        Ok(())
    }

    fn swap(&mut self, fl: &mut Flow) {
        std::mem::swap(&mut fl.f, &mut self.gas_f);
        std::mem::swap(&mut fl.w, &mut self.gas_w);
    }

    /// The face fractions α_k,f upwind of the current phase fluxes.
    fn face_fractions(&self, fl: &mut Flow) -> Result<()> {
        let (ml, mg) = (fl.f.mass_flux, self.gas_f.mass_flux);
        let dev = fl.dev.as_mut();
        for (ph, m, ab) in [(&self.gas, mg, self.alpha_b), (&self.liq, ml, self.alpha_lb)] {
            dev.lincomb(self.flux, 1.0 / ph.rho, m, 0.0, m)?;
            dev.vof(
                VofOp::Faces,
                &VofBufs {
                    flux: Some(self.flux),
                    alpha: Some(ph.alpha),
                    grad: Some(self.grad_alpha),
                    alpha_b: Some(ab),
                    scheme: foehncraft_kernels::vof::SCHEME_UPWIND,
                    dt: 0.0,
                    beta: 1.0,
                    out: Some(ph.alpha_f),
                    ..Default::default()
                },
            )?;
        }
        Ok(())
    }

    /// Advances α_g over `dt` with the gas volume flux of the last step.
    fn advance(&mut self, fl: &mut Flow, dt: f64) -> Result<()> {
        let ni = fl.disc.n_interior;
        let mg = self.gas_f.mass_flux;
        let dev = fl.dev.as_mut();
        dev.lincomb(self.flux, 1.0 / self.gas.rho, mg, 0.0, mg)?;
        dev.vof(
            VofOp::OutflowRate,
            &VofBufs {
                flux: Some(self.flux),
                out: Some(self.work),
                ..Default::default()
            },
        )?;
        let rate = dev.min_max(self.work)?.1;
        let co = rate * dt;
        let n = if co.is_finite() && co > 0.0 {
            ((co / self.section.courant()).ceil() as usize).max(1)
        } else {
            1
        };
        anyhow::ensure!(
            n <= 100_000,
            "[two_fluid]: {n} sub-steps of α in one time step (outflow Courant number {co:.3e}); \
             the time step is far too long for the gas velocity"
        );
        let h = dt / n as f64;
        let fl_host = dev.download(self.flux)?;
        let mut clip = 0.0;
        for _ in 0..n {
            dev.vof(
                VofOp::Faces,
                &VofBufs {
                    flux: Some(self.flux),
                    alpha: Some(self.gas.alpha),
                    grad: Some(self.grad_alpha),
                    alpha_b: Some(self.alpha_b),
                    scheme: foehncraft_kernels::vof::SCHEME_UPWIND,
                    dt: h,
                    beta: 1.0,
                    out: Some(self.gas.alpha_f),
                    ..Default::default()
                },
            )?;
            dev.vof(
                VofOp::Update,
                &VofBufs {
                    flux: Some(self.flux),
                    alpha_f: Some(self.gas.alpha_f),
                    alpha: Some(self.gas.alpha),
                    dt: h,
                    out: Some(self.alpha_new),
                    out2: Some(self.clipped),
                    ..Default::default()
                },
            )?;
            dev.copy(self.alpha_new, self.gas.alpha)?;
            let c = dev.download(self.clipped)?;
            clip += c.iter().copied().exact_sum();
            let af = dev.download(self.gas.alpha_f)?;
            self.totals.volume_out += (ni..fl_host.len())
                .map(|f| fl_host[f] * af[f] * h)
                .exact_sum();
        }
        self.totals.clipped_volume += clip;
        self.totals.substeps = n;
        Ok(())
    }

    /// The interphase arrays of the current fields: K V, the lift and wall
    /// forces, μ_BIT, D, the cell forces b_k and the face arrays.
    fn interphase(&self, fl: &mut Flow) -> Result<()> {
        let nu_t = fl.turbulence.as_ref().map(|t| t.nu_t_buf());
        let (lf, gf) = (fl.f, self.gas_f);
        let dev = fl.dev.as_mut();
        dev.vof(
            VofOp::Gradient,
            &VofBufs {
                alpha: Some(self.gas.alpha),
                out: Some(self.gnorm),
                out3: Some(self.grad_alpha),
                ..Default::default()
            },
        )?;
        dev.two_fluid(
            TwoFluidOp::Cells,
            &TwoFluidBufs {
                alpha: Some(self.gas.alpha),
                ul: Some(lf.u),
                ug: Some(gf.u),
                grad_ul: Some(lf.grad_u),
                wall_y: Some(self.wall_y),
                wall_n: Some(self.wall_n),
                nu_t,
                grad_alpha: Some(self.grad_alpha),
                cell_h: Some(self.cell_h),
                kv: Some(self.kv),
                fv: Some(self.fv),
                mu_bit: Some(self.mu_bit),
                disp: Some(self.disp),
                body_g: Some(self.gas.mat.body),
                body_l: Some(self.liq.mat.body),
                ..Default::default()
            },
            &self.k,
        )?;
        dev.two_fluid(
            TwoFluidOp::Faces,
            &TwoFluidBufs {
                alpha: Some(self.gas.alpha),
                alpha_b: Some(self.alpha_b),
                alpha_gf: Some(self.gas.alpha_f),
                alpha_lf: Some(self.liq.alpha_f),
                mu_eff: Some(self.mu_eff),
                mu_bit: Some(self.mu_bit),
                disp: Some(self.disp),
                mu_g: Some(self.gas.mat.mu),
                mu_l: Some(self.liq.mat.mu),
                off_g: Some(self.gas.mat.jump),
                off_l: Some(self.liq.mat.jump),
                ..Default::default()
            },
            &self.k,
        )
    }

    /// The momentum equation of the phase in the flow's fields (`gas`:
    /// the gas is swapped in), solved with the drag implicit on its own
    /// diagonal and the other phase lagged; `w.d`, `w.ds`, `w.dc` stay
    /// V/a of the row (the elimination adds α_k and the drag).
    fn momentum(
        &self,
        fl: &mut Flow,
        gas: bool,
        dt: f64,
        alpha: f64,
        d: &MaterialBufs,
        res: &mut Residuals,
    ) -> Result<()> {
        // the other phase's velocity: the gas slot holds the liquid while the
        // gas is swapped in, else the gas
        let other = self.gas_f.u;
        let ph = if gas { self.gas } else { self.liq };
        Self::select(fl, &ph.mat, d)?;
        let t = Instant::now();
        fl.update_boundaries()?;
        fl.gradients()?;
        if alpha < 1.0 {
            fl.relaxation_term(alpha)?;
        }
        let (f, w) = (fl.f, fl.w);
        {
            let dev = fl.dev.as_mut();
            // the convecting flux α_k,f ṁ_k and the sink K V − Σ α_k,f ṁ_k
            dev.product(ph.conv, 1.0, ph.alpha_f, f.mass_flux)?;
            dev.net_outflow(ph.conv, self.tmp)?;
            dev.lincomb(ph.sink, 1.0, self.kv, -1.0, self.tmp)?;
            // α_k (∇p − b_k)
            for ((g, b), o) in f
                .grad_p
                .all()
                .into_iter()
                .zip(ph.mat.body.all())
                .zip(self.tmp3.all())
            {
                dev.lincomb(o, 1.0, g, -1.0, b)?;
            }
            // (with the residual fraction of the cell density, as the drag:
            // an empty cell carries the bubbles' terminal slip)
            dev.scale_vector(ph.mat.rho, self.tmp3, ph.gp)?;
            for g in ph.gp.all() {
                dev.lincomb(g, 1.0 / ph.rho, g, 0.0, g)?;
            }
            // the drag of the other phase and the lift and wall forces
            let sign = if gas { 1.0 } else { -1.0 };
            for ((s, u), q) in ph.src.all().into_iter().zip(other.all()).zip(self.fv.all()) {
                dev.product(s, 1.0, self.kv, u)?;
                dev.lincomb(s, 1.0, s, sign, q)?;
            }
        }
        fl.momentum_assembly_faces(ph.conv)?;
        let time = crate::device::TimeTerms {
            a: [1.0, 1.0, 0.0],
            dt,
            u_old: ph.u_old,
            u_old2: ph.u_old,
            explicit: (0.0, ph.u_old),
        };
        fl.dev.momentum_rows(&crate::device::MomentumRowArgs {
            cp: w.cp,
            cn: w.cn,
            bf: w.bf,
            grad_p: ph.gp,
            body_force: [0.0; 3],
            u: f.u,
            alpha,
            theta: 1.0,
            time: Some(time),
            sink: Some(ph.sink),
            val: w.val,
            rhs: w.rhs,
            d: w.d,
            ds: w.ds,
            dc: w.dc,
            diag: w.diag,
        })?;
        for (r, q) in w.rhs.all().into_iter().zip(ph.src.all()) {
            fl.dev.lincomb(r, 1.0, r, 1.0, q)?;
        }
        if fl.settings.algorithm(true) == crate::simplec::Algorithm::Simple {
            fl.dev.copy(w.d, w.dc)?;
        }
        let scale = fl.dev.residual_scale(f.u, w.rhs, w.diag, alpha)?;
        fl.timings.momentum_assembly += t.elapsed().as_secs_f64();
        let t = Instant::now();
        fl.momentum_solve(scale, res)?;
        fl.timings.momentum_solve += t.elapsed().as_secs_f64();
        Ok(())
    }

    /// The predicted flux of the phase in the flow's fields into `w.m`
    /// and its pressure coefficients, after the elimination (the D of
    /// `w` carry α_k).
    fn phase_flux(fl: &mut Flow, ph: &Phase, d: &MaterialBufs, alpha: f64) -> Result<()> {
        Self::select(fl, &ph.mat, d)?;
        fl.update_boundaries()?;
        let w = fl.w;
        let relax = if alpha < 1.0 { w.relax } else { w.zero_face };
        // (the time-step term of the flux with the old-level difference)
        let a = fl.flux_args(Some(w.d), relax, Some((w.ds, ph.x_old, alpha)));
        fl.dev.mass_flux(&a, w.m)?;
        fl.dev.pressure_faces(w.dc, ph.cf)
    }

    /// The old-level Rhie–Chow difference ṁ − ρ ū·S of the phase in the
    /// flow's fields into `x_old`.
    fn finish_phase(fl: &mut Flow, ph: &Phase, d: &MaterialBufs) -> Result<()> {
        Self::select(fl, &ph.mat, d)?;
        fl.update_boundaries()?;
        fl.gradients()?;
        let a = fl.flux_args(None, fl.w.zero_face, None);
        fl.dev.mass_flux(&a, ph.x_old)?;
        fl.dev
            .scaled_difference(ph.x_old, 1.0, fl.f.mass_flux, ph.x_old)
    }

    /// The host copies of α_g and the gas velocity.
    fn download(&mut self, fl: &mut Flow) -> Result<()> {
        let dev = fl.dev.as_mut();
        self.alpha_host = dev.download(self.gas.alpha)?;
        let u = self.gas_f.u;
        self.gas_velocity = [dev.download(u.x)?, dev.download(u.y)?, dev.download(u.z)?];
        self.totals.volume = self
            .alpha_host
            .iter()
            .zip(&self.volume)
            .map(|(a, v)| a * v)
            .exact_sum();
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &a in &self.alpha_host {
            lo = lo.min(a);
            hi = hi.max(a);
        }
        self.totals.alpha_range = [lo, hi];
        fl.phys.vof_alpha = Some(std::sync::Arc::new(self.alpha_host.clone()));
        fl.phys.gas_velocity = Some(std::sync::Arc::new(self.gas_velocity.clone()));
        Ok(())
    }
}

/// Attaches `[two_fluid]` of `case` to the flow (after the start fields
/// and the turbulence model; nothing without `[two_fluid]`).
pub fn attach(fl: &mut Flow, case: &crate::case::FlowCase, mesh: &Mesh, geo: &Geometry) -> Result<()> {
    if case.two_fluid.is_empty() {
        return Ok(());
    }
    let sec = case.two_fluid.clone();
    let kc = TwoFluid::constants(case)?;
    let (nc, nf, ni) = (fl.disc.n_cells, fl.disc.n_faces, fl.disc.n_interior);
    let nb = nf - ni;
    // the inflow α and the flags of the face arrays
    let mut ab = vec![0.0; nb];
    let mut inflow = vec![0u32; nb];
    let mut gas_ub: Vec<Option<[f64; 3]>> = vec![None; nb];
    let mut is_in = vec![false; nb];
    for b in &mesh.boundaries {
        if let Some(i) = sec.inflow.get(&b.name) {
            for f in b.range() {
                if f >= ni && f - ni < nb {
                    ab[f - ni] = i.alpha;
                    inflow[f - ni] = 1;
                    gas_ub[f - ni] = i.velocity;
                    is_in[f - ni] = true;
                }
            }
        }
    }
    // the wall distance and the wall normal
    let walls: Vec<usize> = mesh
        .boundaries
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            matches!(
                case.boundary.get(&p.name),
                Some(crate::bc::BoundaryCondition::Wall { .. })
            )
        })
        .map(|(i, _)| i)
        .collect();
    let (wy, wn) = match foehncraft_mesh::query::wall_distance(mesh, geo, &walls) {
        Some(c) => {
            let mut n = [vec![0.0; nc], vec![0.0; nc], vec![0.0; nc]];
            for (i, q) in c.iter().enumerate().take(nc) {
                let x = geo.cell_centroid[i];
                let d = [x[0] - q.point[0], x[1] - q.point[1], x[2] - q.point[2]];
                let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                let v = if l > 0.0 {
                    [d[0] / l, d[1] / l, d[2] / l]
                } else {
                    q.normal
                };
                for j in 0..3 {
                    n[j][i] = v[j];
                }
            }
            (c.iter().take(nc).map(|q| q.distance).collect(), n)
        }
        None => (
            vec![f64::INFINITY; nc],
            [vec![0.0; nc], vec![0.0; nc], vec![0.0; nc]],
        ),
    };
    let ut = terminal_speed(&kc);
    let gm = (kc.g[0] * kc.g[0] + kc.g[1] * kc.g[1] + kc.g[2] * kc.g[2]).sqrt();
    let up: [f64; 3] = if gm > 0.0 {
        [-kc.g[0] / gm * ut, -kc.g[1] / gm * ut, -kc.g[2] / gm * ut]
    } else {
        [0.0; 3]
    };
    let a0 = sec.initial.unwrap_or(0.0);
    let volume = fl.disc.cell_volume[..nc].to_vec();
    let (gas_f, gas_w) = fl.phase_slot()?;
    let lf = fl.f;
    let dev = fl.dev.as_mut();
    // the device material buffers with a cell force and offsets
    dev.set_body_force(&[vec![0.0; nc], vec![0.0; nc], vec![0.0; nc]])?;
    dev.set_balanced_force(&vec![0.0; nf], &vec![0.0; nb])?;
    dev.two_fluid_inflow(&inflow)?;
    let cell = |dev: &mut dyn crate::device::FlowDevice| dev.alloc(nc);
    let face = |dev: &mut dyn crate::device::FlowDevice| dev.alloc(nf);
    let cell3 = |dev: &mut dyn crate::device::FlowDevice| -> Result<Buf3> {
        Ok(Buf3 {
            x: dev.alloc(nc)?,
            y: dev.alloc(nc)?,
            z: dev.alloc(nc)?,
        })
    };
    let phase = |dev: &mut dyn crate::device::FlowDevice, rho: f64| -> Result<Phase> {
        let rf = face(dev)?;
        dev.fill(rf, rho)?;
        Ok(Phase {
            mat: Material {
                rho: cell(dev)?,
                rho_face: rf,
                mu: face(dev)?,
                body: cell3(dev)?,
                jump: face(dev)?,
            },
            alpha: cell(dev)?,
            alpha_f: face(dev)?,
            conv: face(dev)?,
            sink: cell(dev)?,
            src: cell3(dev)?,
            gp: cell3(dev)?,
            u_old: cell3(dev)?,
            x_old: face(dev)?,
            cf: face(dev)?,
            rho,
        })
    };
    let gas = phase(dev, kc.rho_g)?;
    let liq = phase(dev, kc.rho_l)?;
    let turb = {
        let r = cell(dev)?;
        dev.fill(r, kc.rho_l)?;
        let rf = face(dev)?;
        dev.fill(rf, kc.rho_l)?;
        Material {
            rho: r,
            rho_face: rf,
            mu: face(dev)?,
            body: cell3(dev)?,
            jump: face(dev)?,
        }
    };
    for b in turb.body.all() {
        dev.fill(b, 0.0)?;
    }
    dev.fill(turb.jump, 0.0)?;
    for ph in [&gas, &liq] {
        dev.fill(ph.x_old, 0.0)?;
    }
    let alpha_b = dev.alloc(nb)?;
    dev.upload(alpha_b, &ab)?;
    let alpha_lb = dev.alloc(nb)?;
    dev.upload(alpha_lb, &ab.iter().map(|a| 1.0 - a).collect::<Vec<_>>())?;
    dev.fill(gas.alpha, a0)?;
    let wall_y = cell(dev)?;
    dev.upload(wall_y, &wy)?;
    // h: the shortest distance from a cell centre to a neighbour's
    let mut h = vec![f64::INFINITY; mesh.n_cells];
    for (f, &n) in mesh.neighbour.iter().enumerate() {
        let (p, n) = (mesh.owner[f] as usize, n as usize);
        if p < mesh.n_cells && n < mesh.n_cells {
            let (a, b) = (geo.cell_centroid[p], geo.cell_centroid[n]);
            let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
            h[p] = h[p].min(d);
            h[n] = h[n].min(d);
        }
    }
    for x in &mut h {
        if !x.is_finite() {
            *x = 0.0;
        }
    }
    let cell_h = cell(dev)?;
    dev.upload(cell_h, &h)?;
    let wall_n = cell3(dev)?;
    for (b, v) in wall_n.all().into_iter().zip(&wn) {
        dev.upload(b, v)?;
    }
    // μ_l on every face until the turbulence model sets μ_l + μ_t
    let mu_eff = face(dev)?;
    dev.fill(mu_eff, kc.mu_l)?;
    // the gas: the liquid's start velocity plus the terminal rise; at the
    // inflows the given gas velocity
    for (i, (g, l)) in gas_f.u.all().into_iter().zip(lf.u.all()).enumerate() {
        dev.copy(l, g)?;
        dev.subtract(g, -up[i])?;
    }
    let ubl = [
        dev.download(lf.ub_set.x)?,
        dev.download(lf.ub_set.y)?,
        dev.download(lf.ub_set.z)?,
    ];
    for i in 0..3 {
        let v: Vec<f64> = (0..nb)
            .map(|fb| match gas_ub[fb] {
                Some(u) => u[i],
                None if is_in[fb] => ubl[i][fb] + up[i],
                None => ubl[i][fb],
            })
            .collect();
        dev.upload(gas_f.ub_set.all()[i], &v)?;
    }
    let tf = TwoFluid {
        section: sec,
        k: kc,
        gas,
        liq,
        gas_f,
        gas_w,
        alpha_b,
        alpha_lb,
        alpha_new: cell(dev)?,
        flux: face(dev)?,
        clipped: cell(dev)?,
        work: cell(dev)?,
        grad_alpha: cell3(dev)?,
        gnorm: cell(dev)?,
        wall_y,
        cell_h,
        wall_n,
        kv: cell(dev)?,
        fv: cell3(dev)?,
        mu_bit: cell(dev)?,
        disp: cell(dev)?,
        mu_eff,
        turb,
        tmp: cell(dev)?,
        tmp3: cell3(dev)?,
        u_lag: cell3(dev)?,
        cv: face(dev)?,
        mv: face(dev)?,
        volume,
        totals: Totals::default(),
        alpha_host: vec![a0; nc],
        gas_velocity: [Vec::new(), Vec::new(), Vec::new()],
    };
    let mut tf = Box::new(tf);
    tf.refresh_fractions(fl)?;
    // the start fluxes: ρ_k ū_k·S of both phases
    let d = fl.dev.material_bufs();
    for gas in [true, false] {
        let ph = if gas { tf.gas } else { tf.liq };
        if gas {
            tf.swap(fl);
        }
        TwoFluid::select(fl, &ph.mat, &d)?;
        fl.update_boundaries()?;
        fl.gradients()?;
        let a = fl.flux_args(None, fl.w.zero_face, None);
        let m = fl.f.mass_flux;
        fl.dev.mass_flux(&a, m)?;
        if gas {
            tf.swap(fl);
        }
    }
    tf.face_fractions(fl)?;
    tf.download(fl)?;
    tf.totals.volume_start = tf.totals.volume;
    fl.two_fluid = Some(tf);
    Ok(())
}

/// Debug: cell 80 of several buffers (FOEHN_TF_DEBUG).
fn probe(fl: &mut Flow, tag: &str, bufs: &[(&str, Buf)]) {
    if std::env::var_os("FOEHN_TF_DEBUG").is_none() {
        return;
    }
    let mut line = format!("  @{tag}");
    for (n, b) in bufs {
        if let Ok(v) = fl.dev.download(*b) {
            line += &format!(" {n}={:.4e}", v[80.min(v.len() - 1)]);
        }
    }
    eprintln!("{line}");
}

/// Debug: the finiteness and range of a buffer (FOEHN_TF_DEBUG).
fn chk(fl: &mut Flow, name: &str, b: Buf) {
    if std::env::var_os("FOEHN_TF_DEBUG").is_none() {
        return;
    }
    if let Ok(v) = fl.dev.download(b) {
        let bad = v.iter().filter(|x| !x.is_finite()).count();
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &x in v.iter().filter(|x| x.is_finite()) {
            lo = lo.min(x);
            hi = hi.max(x);
        }
        eprintln!("  {name:<14} n {:6} bad {bad:6} [{lo:.4e}, {hi:.4e}]", v.len());
    }
}

/// One backward-Euler time step of size `dt` of a two-fluid flow.
pub fn step(
    tr: &mut crate::transient::Transient,
    fl: &mut Flow,
    dt: f64,
) -> Result<crate::transient::StepReport> {
    let mut tf = fl.two_fluid.take().context("a two-fluid flow")?;
    let r = step_inner(tr, fl, &mut tf, dt);
    if let Err(e) = &r {
        if std::env::var_os("FOEHN_TF_DEBUG").is_some() {
            eprintln!("two-fluid step: {e:#}");
        }
    }
    fl.two_fluid = Some(tf);
    r
}

fn step_inner(
    tr: &mut crate::transient::Transient,
    fl: &mut Flow,
    tf: &mut TwoFluid,
    dt: f64,
) -> Result<crate::transient::StepReport> {
    anyhow::ensure!(dt > 0.0 && dt.is_finite(), "time step {dt} must be > 0");
    let s = fl.settings.clone();
    let alpha = s.alpha_u(true);
    let alpha_p = s.alpha_p(true);
    let d = fl.dev.material_bufs();
    // 1. α_g of the new level, the face fractions, the old levels
    tf.advance(fl, dt)?;
    tf.refresh_fractions(fl)?;
    tf.face_fractions(fl)?;
    for (ph, u) in [(tf.gas, tf.gas_f.u), (tf.liq, fl.f.u)] {
        for (a, b) in u.all().into_iter().zip(ph.u_old.all()) {
            fl.dev.copy(a, b)?;
        }
    }
    let levels = tr.clock.levels;
    if let Some(t) = fl.turbulence.as_mut() {
        t.step_start(fl.dev.as_mut(), levels >= 1, None)?;
    }
    let scalar_step = fl.turbulence.as_ref().map(|_| crate::scalar::StepTime {
        time: tr.clock.time + dt,
        dt,
        a: crate::scalar::bdf(false, 1.0),
    });
    let mut first = Residuals::default();
    let mut last = Residuals::default();
    let mut outer = 0;
    tr.outer_log.clear();
    for kk in 1..=s.outer_iterations {
        outer = kk;
        fl.iteration += 1;
        let mut res = Residuals {
            iteration: fl.iteration,
            ..Default::default()
        };
        // the liquid's gradients for the lift (the gas's follow in its
        // momentum)
        TwoFluid::select(fl, &tf.liq.mat, &d)?;
        fl.update_boundaries()?;
        fl.gradients()?;
        tf.interphase(fl)?;
        for (n, b) in [
            ("alpha_g", tf.gas.alpha), ("alpha_gf", tf.gas.alpha_f), ("alpha_lf", tf.liq.alpha_f),
            ("rho_g", tf.gas.mat.rho), ("kv", tf.kv), ("fv.y", tf.fv.y), ("mu_bit", tf.mu_bit),
            ("disp", tf.disp), ("body_g.y", tf.gas.mat.body.y), ("mu_g", tf.gas.mat.mu),
            ("mu_l", tf.liq.mat.mu), ("off_g", tf.gas.mat.jump), ("off_l", tf.liq.mat.jump),
            ("ug.y", tf.gas_f.u.y), ("mg", tf.gas_f.mass_flux), ("ml", fl.f.mass_flux),
        ] {
            chk(fl, n, b);
        }
        // the gas
        let mut rg = Residuals::default();
        tf.swap(fl);
        let r = tf.momentum(fl, true, dt, alpha, &d, &mut rg);
        tf.swap(fl);
        r?;
        for (a, b) in fl.f.u.all().into_iter().zip(tf.u_lag.all()) {
            fl.dev.copy(a, b)?;
        }
        for (n, b) in [("g d", tf.gas_w.d), ("g dc", tf.gas_w.dc), ("g m", tf.gas_w.m), ("g cf", tf.gas.cf), ("g u.y", tf.gas_f.u.y)] {
            chk(fl, n, b);
        }
        // the liquid
        tf.momentum(fl, false, dt, alpha, &d, &mut res)?;
        // the partial elimination of the drag, then the predicted fluxes
        let (gw, lw) = (tf.gas_w, fl.w);
        fl.dev.two_fluid_eliminate(&EliminateBufs {
            alpha: tf.gas.alpha,
            kv: tf.kv,
            diag_g: gw.diag,
            diag_l: lw.diag,
            ug: tf.gas_f.u,
            ul: fl.f.u,
            ul_lag: tf.u_lag,
            d_g: [gw.d, gw.ds, gw.dc],
            d_l: [lw.d, lw.ds, lw.dc],
        })?;
        let (g, l) = (tf.gas, tf.liq);
        tf.swap(fl);
        let r = TwoFluid::phase_flux(fl, &g, &d, alpha);
        tf.swap(fl);
        r?;
        TwoFluid::phase_flux(fl, &l, &d, alpha)?;
        for (n, b) in [("l d", fl.w.d), ("l dc", fl.w.dc), ("l m", fl.w.m), ("l cf", tf.liq.cf), ("l u.y", fl.f.u.y)] {
            chk(fl, n, b);
        }
        tf.totals.gas_residual = rg.momentum.iter().copied().fold(0.0, f64::max);
        // the pressure correction of the volume continuity
        let t = Instant::now();
        let (gw, w) = (tf.gas_w, fl.w);
        {
            let dev = fl.dev.as_mut();
            let (g, l) = (&tf.gas, &tf.liq);
            dev.product(tf.cv, 1.0 / l.rho, l.alpha_f, l.cf)?;
            dev.product(tf.flux, 1.0 / g.rho, g.alpha_f, g.cf)?;
            dev.lincomb(tf.cv, 1.0, tf.cv, 1.0, tf.flux)?;
            dev.product(tf.mv, 1.0 / l.rho, l.alpha_f, w.m)?;
            dev.product(tf.flux, 1.0 / g.rho, g.alpha_f, gw.m)?;
            dev.lincomb(tf.mv, 1.0, tf.mv, 1.0, tf.flux)?;
            dev.net_outflow(tf.mv, w.defect)?;
            dev.throughflow(tf.mv, w.through)?;
            let dd = dev.sum_sq(w.defect)?.sqrt();
            let tt = dev.sum_sq(w.through)?.sqrt();
            res.continuity = dd / tt.max(f64::MIN_POSITIVE);
            dev.pressure_rows(tf.cv, tf.mv, w.pval, w.prhs)?;
        }
        for (n, b) in [("cv", tf.cv), ("mv", tf.mv), ("pval", fl.w.pval), ("prhs", fl.w.prhs)] {
            chk(fl, n, b);
        }
        if let Some(r) = fl.reference_cell {
            let pos = fl.disc.sell.pos(r, 0);
            fl.dev.scale_entry(w.pval, pos, 2.0)?;
        }
        fl.timings.pressure_assembly += t.elapsed().as_secs_f64();
        fl.pressure_amg(kk == 1)?;
        let t = Instant::now();
        fl.dev.fill(w.pc, 0.0)?;
        let (it, _) = fl.dev.solve(
            w.pval,
            w.prhs,
            w.pc,
            s.pressure.preconditioner,
            true,
            s.pressure_tol(true),
            s.pressure.max_iter,
        )?;
        res.pressure_iterations = it;
        if let Some(shift) = fl.reference_value(w.pc)? {
            fl.dev.subtract(w.pc, shift)?;
        }
        fl.timings.pressure_solve += t.elapsed().as_secs_f64();
        // the corrections of both phases
        let t = Instant::now();
        {
            let dev = fl.dev.as_mut();
            dev.gradient(
                Lsq::Pressure,
                w.pc,
                Mask::PressureCorrection,
                w.zero_boundary,
                w.gpc,
            )?;
            dev.correct_fluxes(tf.liq.cf, w.m, w.pc, w.mc)?;
            dev.correct_velocity(w.dc, w.gpc, fl.f.u)?;
            dev.correct_fluxes(tf.gas.cf, gw.m, w.pc, gw.mc)?;
            dev.correct_velocity(gw.dc, w.gpc, tf.gas_f.u)?;
            dev.lincomb(fl.f.p, 1.0, fl.f.p, alpha_p, w.pc)?;
        }
        std::mem::swap(&mut fl.f.mass_flux, &mut fl.w.mc);
        std::mem::swap(&mut tf.gas_f.mass_flux, &mut tf.gas_w.mc);
        fl.timings.pressure_assembly += t.elapsed().as_secs_f64();
        // the turbulence of the liquid (ρ_l), then its face viscosity μ_l +
        // μ_t for the phases
        if res.finite() && fl.turbulence.is_some() {
            TwoFluid::select(fl, &tf.turb, &d)?;
            fl.dev.fill(d.mu, tf.k.mu_l)?;
            if let Some(t) = fl.turbulence.as_mut() {
                t.amg_rebuild = kk == 1;
            }
            fl.solve_scalars(scalar_step, &mut res)?;
            fl.dev.copy(d.mu, tf.mu_eff)?;
        }
        if !res.finite() || !tf.totals.gas_residual.is_finite() {
            fl.diverged = Some(res.clone());
            return Err(diag::err(Diagnostic::new(
                diag::RUN_DIVERGED,
                format!(
                    "the solution diverged in time step {} (t = {:.6e} s, Δt = {dt:.3e} s; \
                     liquid residuals {:?}, gas {:.3e}, continuity {:.3e})",
                    tr.clock.step + 1,
                    tr.clock.time + dt,
                    res.momentum,
                    tf.totals.gas_residual,
                    res.continuity
                ),
            )));
        }
        if std::env::var_os("FOEHN_TF_DEBUG").is_some() {
            eprintln!("  #it {kk} liq {:?} gas {:?} cont {:.3e}", res.momentum, rg.momentum, res.continuity);
        }
        for (a, b) in res.momentum.iter_mut().zip(rg.momentum) {
            *a = a.max(b);
        }
        if kk == 1 {
            first = res.clone();
        } else {
            first.pressure_iterations += res.pressure_iterations;
        }
        last = res.clone();
        tr.outer_log.push(res.clone());
        if kk > 1 && res.max() < s.outer_tolerance {
            break;
        }
    }
    // the old-level differences of the next step
    let t = Instant::now();
    let (g, l) = (tf.gas, tf.liq);
    tf.swap(fl);
    let r = TwoFluid::finish_phase(fl, &g, &d);
    tf.swap(fl);
    r?;
    TwoFluid::finish_phase(fl, &l, &d)?;
    tf.download(fl)?;
    let cfl = tr.courant(fl, dt)?;
    fl.timings.gradients += t.elapsed().as_secs_f64();
    tr.clock.time += dt;
    tr.clock.step += 1;
    tr.clock.dt_prev2 = tr.clock.dt_prev;
    tr.clock.dt_prev = dt;
    tr.clock.levels = (tr.clock.levels + 1).min(2);
    Ok(crate::transient::StepReport {
        step: tr.clock.step,
        time: tr.clock.time,
        dt,
        residuals: first,
        last,
        outer,
        cfl,
        error: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn air_water() -> k::Constants<f64> {
        k::Constants {
            rho_l: 998.2,
            rho_g: 1.2,
            mu_l: 1.0e-3,
            mu_g: 1.8e-5,
            sigma: 0.072,
            d: 4.0e-3,
            g: [0.0, -9.81, 0.0],
            lift: k::LIFT_TOMIYAMA,
            c_lift: 0.0,
            lift_damping: 1,
            wall: k::WALL_HOSOKAWA,
            c_td: 1.0,
            sato: 0.6,
            alpha_min: 1e-6,
        }
    }

    /// The lift damping is 0 inside a bubble radius, 1 beyond a diameter,
    /// smooth and monotone between.
    #[test]
    fn lift_damping_near_walls() {
        assert_eq!(k::lift_wall_damping(0.3), 0.0);
        assert_eq!(k::lift_wall_damping(0.5), 0.0);
        assert_eq!(k::lift_wall_damping(1.0), 1.0);
        assert_eq!(k::lift_wall_damping(f64::INFINITY), 1.0);
        assert!((k::lift_wall_damping(0.75_f64) - 0.5).abs() < 1e-15);
        let mut last = 0.0;
        for i in 0..=100 {
            let v: f64 = k::lift_wall_damping(0.5 + 0.005 * f64::from(i));
            assert!(v >= last);
            last = v;
        }
    }

    /// The partial elimination solves the drag pair exactly and gives the
    /// separate rows' D without drag.
    #[test]
    fn elimination_solves_the_drag_pair() {
        let (kv, ag, al) = (3.0_f64, 5.0_f64, 40.0_f64);
        let (ug0, ul0, lag) = (0.7_f64, 0.4_f64, 0.45_f64);
        let (ug, ul) = k::eliminate(kv, ag, al, ug0, ul0, lag);
        let hg = ag * ug0 - kv * lag;
        let hl = al * ul0 - kv * ug0;
        assert!((ag * ug - kv * ul - hg).abs() < 1e-12);
        assert!((al * ul - kv * ug - hl).abs() < 1e-12);
        let (dg, dl): (f64, f64) = k::eliminate_d(2.0, 0.0, 0.25, 0.1, 0.3);
        assert!((dg - 0.025).abs() < 1e-15 && (dl - 0.225).abs() < 1e-15);
    }

    /// The distorted regime of Ishii & Zuber: C_D = ⅔ √Eo gives the rise
    /// speed U = √2 (g Δρ σ)^¼ / √ρ_l, independent of d.
    #[test]
    fn terminal_speed_of_a_distorted_bubble() {
        let c = air_water();
        let u = terminal_speed(&c);
        let exact = 2f64.sqrt() * (9.81 * (998.2 - 1.2) * 0.072f64).powf(0.25) / 998.2f64.sqrt();
        assert!((u / exact - 1.0).abs() < 1e-9, "{u} vs {exact}");
        assert!((u - 0.2306).abs() < 1e-3, "{u}");
    }

    /// Stokes drag at a vanishing slip: K = 18 μ_m α / d².
    #[test]
    fn drag_tends_to_stokes() {
        let c = air_water();
        let a = 0.01;
        let kk = k::drag_k(&c, a, 0.0);
        let mu_m = c.mu_l / (1.0 - a);
        assert!((kk / (18.0 * mu_m * a / (c.d * c.d)) - 1.0).abs() < 1e-12);
    }

    /// C_L of Tomiyama: positive for 4 mm air bubbles in water (towards
    /// the wall in up-flow), negative for 7 mm ones.
    #[test]
    fn lift_changes_sign_with_size() {
        let mut c = air_water();
        assert!(k::lift_tomiyama(&c, 900.0) > 0.2);
        c.d = 7.0e-3;
        assert!(k::lift_tomiyama(&c, 1500.0) < 0.0);
    }
}
