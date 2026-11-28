//! Wall boiling of the two-fluid model (L3, stage MP4, D357, D358): the
//! energy equations of the liquid and the vapour, the heat flux
//! partitioning of Kurul & Podowski at heated walls and the heat and mass
//! transfer at the bubble interfaces; the kernels are those of
//! [`foehncraft_kernels::boiling`].
//!
//! **Equations.** Each phase k has a temperature equation of constant
//! properties (the saturation state of the system pressure by IAPWS-IF97
//! for the vapour, the given c_p and conductivity for the liquid),
//!
//! α_k ρ_k c_p,k (∂T_k/∂t + u_k·∇T_k) = ∇·(α_k (k_k + c_p,k μ_t,k/Pr_t)
//!     ∇T_k) + S_k,
//!
//! with S_k of [`foehncraft_kernels::boiling::energy_sources`]. The vapour
//! made per volume, Γ (kg/(m³ s)), enters α_g as Γ/ρ_g after the advection
//! of a step and the volume continuity as Γ (1/ρ_g − 1/ρ_l).
//!
//! **Discretisation.** Per step: (1) after the advection of α_g, the wall
//! split of every heated face and the phase change of every cell with the
//! temperatures and slip of the step start (the rate lags one step, as
//! MP2); (2) the outer iterations of the two momentum equations with the
//! continuity source; (3) the two temperature equations with the final
//! fluxes, backward Euler, upwind convection (the scalar rows of the
//! single fluid with the heat capacity α_k ρ_k c_p,k), heated walls as
//! cell sources, the other walls adiabatic.

use super::{Phase, TwoFluid};
use crate::device::{BoilingBufs, BoilingOp, Buf, Buf3, ScalarFaceArgs, ScalarRowArgs, ScalarTime};
use crate::diag::{self, Diagnostic};
use crate::simplec::Flow;
use anyhow::{Context, Result};
use foehncraft_core::exact::ExactSumIter;
use foehncraft_kernels::boiling as kb;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// `[two_fluid.boiling]`: boiling at heated walls and the energy equations
/// of both phases (MP4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BoilingSection {
    /// System pressure (Pa): the saturation state (T_sat, h_fg) and the
    /// saturated vapour by IAPWS-IF97.
    pub pressure: f64,
    /// Specific heat of the liquid c_p,l (J/(kg K)).
    pub liquid_specific_heat: f64,
    /// Thermal conductivity of the liquid k_l (W/(m K)).
    pub liquid_conductivity: f64,
    /// The liquid temperature at the start (K; default: the first
    /// inflow temperature, else T_sat).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial_temperature: Option<f64>,
    /// The liquid temperature (K) per inflow boundary.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub inflow_temperature: BTreeMap<String, f64>,
    /// The heat flux into the fluid (W/m²) per heated wall boundary.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub heat_flux: BTreeMap<String, f64>,
    /// Pr_t of the turbulent conductivity (default 0.85).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub turbulent_prandtl: Option<f64>,
    /// m and p of the site density (m ΔT_sup)^p of Lemmert & Chawla
    /// (default 185 and 1.805).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nucleation_m: Option<f64>,
    /// See `nucleation_m`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nucleation_p: Option<f64>,
    /// d_ref (m), ΔT_ref (K) and d_max (m) of the departure diameter of
    /// Tolubinsky & Kostanchuk (default 0.6 mm, 45 K, 1.4 mm).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub departure_reference: Option<f64>,
    /// See `departure_reference`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub departure_subcooling: Option<f64>,
    /// See `departure_reference`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub departure_max: Option<f64>,
    /// The influence factor K of the bubble area (default 4).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub influence_factor: Option<f64>,
    /// The largest bubble area fraction (default 1).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub area_max: Option<f64>,
    /// y⁺ of the reference liquid temperature of the quenching and the
    /// subcooling (default 250; 0: the near-wall cell).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_y_plus: Option<f64>,
    /// The bulk bubble diameter (m) at the subcooling `subcooling_cold`
    /// (K) and at `subcooling_hot`, linear between (default 0.1 mm at
    /// 13.5 K and 2 mm at −5 K, Kurul & Podowski).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub diameter_cold: Option<f64>,
    /// See `diameter_cold`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subcooling_cold: Option<f64>,
    /// See `diameter_cold`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub diameter_hot: Option<f64>,
    /// See `diameter_cold`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subcooling_hot: Option<f64>,
    /// Nu of the vapour side of the interface (default 2).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vapour_nusselt: Option<f64>,
    /// The dryout break point α_dry of the near-wall vapour fraction: above
    /// it a growing part of the wall heats the vapour (default 0.9, Yoon
    /// et al. 2017; 1 turns the dryout off).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dryout_void: Option<f64>,
}

fn bad(key: &str, msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::TWO_FLUID_CONFIG, msg).key(key.to_string()))
}

/// The saturation state of `[two_fluid.boiling]`.
fn saturation(b: &BoilingSection) -> Result<crate::steam::Saturation> {
    crate::steam::saturation_p(b.pressure).ok_or_else(|| {
        bad(
            "two_fluid.boiling.pressure",
            format!(
                "[two_fluid.boiling] pressure {} Pa is outside the saturation line of water \
                 (611.2 Pa to 22.06 MPa)",
                b.pressure
            ),
        )
    })
}

/// The ranges of `[two_fluid.boiling]`.
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let t = &case.two_fluid;
    let Some(b) = &t.boiling else {
        return Ok(());
    };
    let sat = saturation(b)?;
    let pos = |x: f64| x > 0.0 && x.is_finite();
    if !(pos(b.liquid_specific_heat) && pos(b.liquid_conductivity)) {
        return Err(bad(
            "two_fluid.boiling",
            "[two_fluid.boiling] liquid_specific_heat and liquid_conductivity must be > 0".into(),
        ));
    }
    if let Some(g) = &t.gas {
        let rv = sat.vapour.rho;
        if ((g.density - rv) / rv).abs() > 0.02 {
            return Err(bad(
                "two_fluid.gas.density",
                format!(
                    "[two_fluid.gas] density {} differs from that of saturated steam at \
                     {} Pa ({rv:.4} kg/m³, IAPWS-IF97); set it to that value",
                    g.density, b.pressure
                ),
            ));
        }
    }
    if case.turbulence.model == crate::turbulence::TurbulenceModel::Laminar {
        return Err(bad(
            "two_fluid.boiling",
            "[two_fluid.boiling] needs a turbulence model with k (the wall functions of the \
             heat flux partitioning)"
                .into(),
        ));
    }
    for (name, q) in &b.heat_flux {
        if !q.is_finite() {
            return Err(bad(
                "two_fluid.boiling.heat_flux",
                format!("[two_fluid.boiling.heat_flux] {name} = {q} must be finite"),
            ));
        }
        if !matches!(
            case.boundary.get(name),
            Some(crate::bc::BoundaryCondition::Wall { .. })
        ) {
            return Err(bad(
                "two_fluid.boiling.heat_flux",
                format!("[two_fluid.boiling.heat_flux] names `{name}`, which is not a wall"),
            ));
        }
    }
    for (name, v) in &b.inflow_temperature {
        if !pos(*v) {
            return Err(bad(
                "two_fluid.boiling.inflow_temperature",
                format!("[two_fluid.boiling.inflow_temperature] {name} = {v} K must be > 0"),
            ));
        }
        if !case.boundary.contains_key(name) {
            return Err(bad(
                "two_fluid.boiling.inflow_temperature",
                format!(
                    "[two_fluid.boiling.inflow_temperature] names the boundary `{name}`, which \
                     the case does not have"
                ),
            ));
        }
    }
    let k = constants(case, b, &sat);
    for (key, v) in [
        ("turbulent_prandtl", b.turbulent_prandtl.unwrap_or(0.85)),
        ("nucleation_m", k.site_m),
        ("nucleation_p", k.site_p),
        ("departure_reference", k.dw_ref),
        ("departure_subcooling", k.dw_dt),
        ("departure_max", k.dw_max),
        ("influence_factor", k.influence),
        ("area_max", k.area_max),
        ("diameter_cold", k.d_cold),
        ("diameter_hot", k.d_hot),
        ("vapour_nusselt", k.nu_g),
    ] {
        if !pos(v) {
            return Err(bad(
                &format!("two_fluid.boiling.{key}"),
                format!("[two_fluid.boiling] {key} {v} must be > 0"),
            ));
        }
    }
    if !(k.alpha_dry > 0.0 && k.alpha_dry <= 1.0) {
        return Err(bad(
            "two_fluid.boiling.dryout_void",
            format!(
                "[two_fluid.boiling] dryout_void {} must be in (0, 1]",
                k.alpha_dry
            ),
        ));
    }
    if !(k.y_ref >= 0.0 && k.y_ref.is_finite()) {
        return Err(bad(
            "two_fluid.boiling.reference_y_plus",
            format!(
                "[two_fluid.boiling] reference_y_plus {} must be >= 0",
                k.y_ref
            ),
        ));
    }
    if k.dt_cold <= k.dt_hot {
        return Err(bad(
            "two_fluid.boiling.subcooling_cold",
            format!(
                "[two_fluid.boiling] subcooling_cold {} K must be above subcooling_hot {} K",
                k.dt_cold, k.dt_hot
            ),
        ));
    }
    Ok(())
}

/// The kernel constants of `b` at the saturation state `sat`.
fn constants(
    case: &crate::case::FlowCase,
    b: &BoilingSection,
    sat: &crate::steam::Saturation,
) -> kb::Constants<f64> {
    let g = case
        .gravity
        .as_ref()
        .map(|g| g.vector)
        .unwrap_or([0.0, 0.0, 0.0]);
    let v = &sat.vapour;
    kb::Constants {
        t_sat: sat.t,
        h_fg: sat.h_fg(),
        rho_l: case.rho(),
        cp_l: b.liquid_specific_heat,
        k_l: b.liquid_conductivity,
        mu_l: case.mu(),
        rho_g: case.two_fluid.gas.as_ref().map_or(v.rho, |g| g.density),
        cp_g: v.cp,
        k_g: crate::steam::conductivity(v),
        mu_g: case.two_fluid.gas.as_ref().map_or(1e-5, |g| g.viscosity),
        g: (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt(),
        site_m: b.nucleation_m.unwrap_or(185.0),
        site_p: b.nucleation_p.unwrap_or(1.805),
        dw_ref: b.departure_reference.unwrap_or(0.6e-3),
        dw_dt: b.departure_subcooling.unwrap_or(45.0),
        dw_max: b.departure_max.unwrap_or(1.4e-3),
        cole_cd: 1.0,
        influence: b.influence_factor.unwrap_or(4.0),
        area_max: b.area_max.unwrap_or(1.0),
        wait: 0.8,
        y_ref: b.reference_y_plus.unwrap_or(250.0),
        d_cold: b.diameter_cold.unwrap_or(0.1e-3),
        dt_cold: b.subcooling_cold.unwrap_or(13.5),
        d_hot: b.diameter_hot.unwrap_or(2.0e-3),
        dt_hot: b.subcooling_hot.unwrap_or(-5.0),
        nu_g: b.vapour_nusselt.unwrap_or(2.0),
        alpha_dry: b.dryout_void.unwrap_or(0.9),
        alpha_min: 1e-6,
    }
}

/// The surface tension of `[two_fluid.boiling]` at T_sat (N/m, IAPWS).
pub fn surface_tension(b: &BoilingSection) -> Option<f64> {
    crate::steam::saturation_p(b.pressure).map(|s| crate::steam::surface_tension(s.t))
}

/// The totals of the boiling (for `result.json` and `two_fluid.csv`).
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema, PartialEq)]
pub struct BoilingTotals {
    /// Saturation temperature (K) and enthalpy of evaporation (J/kg).
    pub t_sat: f64,
    /// See `t_sat`.
    pub h_fg: f64,
    /// Heat into the fluid through the heated walls (W).
    pub wall_heat: f64,
    /// Of that, the evaporation q_e (W).
    pub wall_evaporation_heat: f64,
    /// Of that, the heat into the vapour at dry wall (W).
    pub wall_vapour_heat: f64,
    /// Vapour made at the walls (kg/s) and at the interfaces (kg/s, net;
    /// negative: condensation).
    pub wall_vapour: f64,
    /// See `wall_vapour`.
    pub interface_vapour: f64,
    /// The largest wall temperature (K).
    pub max_wall_temperature: f64,
    /// Enthalpy flowing out less in through the boundaries (W), liquid
    /// enthalpy c_p,l T, vapour c_p,l T_sat + h_fg + c_p,g (T − T_sat).
    pub enthalpy_out: f64,
    /// Iterations of the two temperature solves of the last step.
    pub energy_iterations: [usize; 2],
}

/// The boiling state of a two-fluid flow.
#[derive(Debug)]
pub struct Boil {
    /// The kernel constants.
    pub k: kb::Constants<f64>,
    prt: f64,
    b: BoilingBufs,
    t_l_old: Buf,
    t_g_old: Buf,
    bc_l: (Buf, Buf),
    bc_g: (Buf, Buf),
    conv: Buf,
    gam: Buf,
    grad: Buf3,
    psi: Buf,
    /// Per heated slot: the face, area and heat flux (host).
    slot_face: Vec<usize>,
    slot_area: Vec<f64>,
    slot_q: Vec<f64>,
    /// Liquid temperature of the inflow faces (host, per boundary face;
    /// NaN elsewhere).
    t_in: Vec<f64>,
    /// The totals.
    pub totals: BoilingTotals,
}

/// Attaches `[two_fluid.boiling]` (after the two-fluid state).
pub(super) fn attach(
    fl: &mut Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
) -> Result<Option<Box<Boil>>> {
    let Some(b) = &case.two_fluid.boiling else {
        return Ok(None);
    };
    let sat = saturation(b)?;
    let k = constants(case, b, &sat);
    let (nc, nf, ni) = (fl.disc.n_cells, fl.disc.n_faces, fl.disc.n_interior);
    let nb = nf - ni;
    // the heated slots, ordered by their near-wall cell
    let mut slots: Vec<(u32, usize, f64, f64, f64)> = Vec::new();
    for p in &mesh.boundaries {
        if let Some(&q) = b.heat_flux.get(&p.name) {
            for f in p.range() {
                let c = mesh.owner[f];
                let s = geo.face_area[f];
                let a = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
                let n = [s[0] / a, s[1] / a, s[2] / a];
                let (xc, xf) = (geo.cell_centroid[c as usize], geo.face_centroid[f]);
                let y = ((xf[0] - xc[0]) * n[0] + (xf[1] - xc[1]) * n[1] + (xf[2] - xc[2]) * n[2])
                    .abs();
                slots.push((c, f, a, y, q));
            }
        }
    }
    slots.sort_by_key(|s| (s.0, s.1));
    let mut wcell: Vec<u32> = Vec::new();
    let mut wptr: Vec<u32> = vec![0];
    for (i, s) in slots.iter().enumerate() {
        if wcell.last() != Some(&s.0) {
            if !wcell.is_empty() {
                wptr.push(i as u32);
            }
            wcell.push(s.0);
        }
    }
    if !wcell.is_empty() {
        wptr.push(slots.len() as u32);
    }
    let cells: Vec<u32> = slots.iter().map(|s| s.0).collect();
    let areas: Vec<f64> = slots.iter().map(|s| s.2).collect();
    let ys: Vec<f64> = slots.iter().map(|s| s.3).collect();
    let qs: Vec<f64> = slots.iter().map(|s| s.4).collect();
    // the inflow temperatures and the boundary rules of both phases
    let mut t_in = vec![f64::NAN; nb];
    for p in &mesh.boundaries {
        if let Some(&t) = b.inflow_temperature.get(&p.name) {
            for f in p.range() {
                if f >= ni {
                    t_in[f - ni] = t;
                }
            }
        }
    }
    let t0 = b
        .initial_temperature
        .or_else(|| b.inflow_temperature.values().next().copied())
        .unwrap_or(sat.t);
    let dev = fl.dev.as_mut();
    dev.boiling_slots(&cells, &areas, &ys, &qs, &wcell, &wptr)?;
    let ns = slots.len().max(1);
    let cell = |dev: &mut dyn crate::device::FlowDevice| -> Result<Buf> {
        let x = dev.alloc(nc)?;
        dev.fill(x, 0.0)?;
        Ok(x)
    };
    let slot = |dev: &mut dyn crate::device::FlowDevice| -> Result<Buf> {
        let x = dev.alloc(ns)?;
        dev.fill(x, 0.0)?;
        Ok(x)
    };
    let cell3 = |dev: &mut dyn crate::device::FlowDevice| -> Result<Buf3> {
        Ok(Buf3 {
            x: cell(dev)?,
            y: cell(dev)?,
            z: cell(dev)?,
        })
    };
    let bc = |dev: &mut dyn crate::device::FlowDevice, fixed: &dyn Fn(usize) -> Option<f64>| {
        let a: Vec<f64> = (0..nb)
            .map(|i| if fixed(i).is_some() { 0.0 } else { 1.0 })
            .collect();
        let v: Vec<f64> = (0..nb).map(|i| fixed(i).unwrap_or(0.0)).collect();
        let (ba, bb) = (dev.alloc(nb)?, dev.alloc(nb)?);
        dev.upload(ba, &a)?;
        dev.upload(bb, &v)?;
        anyhow::Ok((ba, bb))
    };
    let bc_l = bc(dev, &|i| t_in[i].is_finite().then_some(t_in[i]))?;
    let bc_g = bc(dev, &|i| t_in[i].is_finite().then_some(sat.t))?;
    let t_l = cell(dev)?;
    dev.fill(t_l, t0)?;
    let t_g = cell(dev)?;
    dev.fill(t_g, sat.t)?;
    let bufs = BoilingBufs {
        t_l,
        t_g,
        tk: cell(dev)?,
        t_w: slot(dev)?,
        q_l: slot(dev)?,
        m_e: slot(dev)?,
        q_e: slot(dev)?,
        area: slot(dev)?,
        q_g: slot(dev)?,
        q_wall: cell(dev)?,
        m_wall: cell(dev)?,
        q_wall_g: cell(dev)?,
        alpha: cell(dev)?,
        ul: cell3(dev)?,
        ug: cell3(dev)?,
        alpha_out: cell(dev)?,
        diam: cell(dev)?,
        gamma_i: cell(dev)?,
        gamma: cell(dev)?,
        h_l: cell(dev)?,
        h_g: cell(dev)?,
        su_l: cell(dev)?,
        sp_l: cell(dev)?,
        su_g: cell(dev)?,
        sp_g: cell(dev)?,
        cap_l: cell(dev)?,
        cap_g: cell(dev)?,
        src: cell(dev)?,
    };
    // the bubble size of the start temperature
    let d0 = kb::bubble_diameter(&k, t0);
    dev.fill(bufs.diam, d0)?;
    let boil = Boil {
        k,
        prt: b.turbulent_prandtl.unwrap_or(0.85),
        t_l_old: cell(dev)?,
        t_g_old: cell(dev)?,
        bc_l,
        bc_g,
        conv: dev.alloc(nf)?,
        gam: dev.alloc(nf)?,
        grad: cell3(dev)?,
        psi: cell(dev)?,
        b: bufs,
        slot_face: slots.iter().map(|s| s.1).collect(),
        slot_area: areas,
        slot_q: qs,
        t_in,
        totals: BoilingTotals {
            t_sat: sat.t,
            h_fg: sat.h_fg(),
            ..Default::default()
        },
    };
    dev.fill(boil.psi, 1.0)?;
    dev.set_mass_source(Some(boil.b.src))?;
    Ok(Some(Box::new(boil)))
}

impl Boil {
    /// The bubble diameter per cell (the forces of the two fluids).
    pub fn diameter(&self) -> Buf {
        self.b.diam
    }

    /// The continuity source per cell (m³/s).
    pub fn source(&self) -> Buf {
        self.b.src
    }
}

/// Step (1): the wall split and the phase change of the step `dt` from α*
/// (in the gas's α), the temperatures and velocities of the step start;
/// α_g of the new level into the gas's α.
pub(super) fn phase_change(fl: &mut Flow, tf: &mut TwoFluid, dt: f64) -> Result<()> {
    let Some(bo) = tf.boil.as_mut() else {
        return Ok(());
    };
    let tk = fl
        .turbulence
        .as_ref()
        .map(|t| t.k_buf())
        .context("turbulence")?;
    let dev = fl.dev.as_mut();
    let mut b = bo.b;
    dev.copy(tk, b.tk)?;
    dev.copy(tf.gas.alpha, b.alpha)?;
    for (s, t) in fl.f.u.all().into_iter().zip(b.ul.all()) {
        dev.copy(s, t)?;
    }
    for (s, t) in tf.gas_f.u.all().into_iter().zip(b.ug.all()) {
        dev.copy(s, t)?;
    }
    dev.copy(b.t_l, bo.t_l_old)?;
    dev.copy(b.t_g, bo.t_g_old)?;
    b.t_l = bo.t_l_old;
    b.t_g = bo.t_g_old;
    dev.boiling(BoilingOp::Wall, &b, &bo.k, dt)?;
    dev.boiling(BoilingOp::WallCells, &b, &bo.k, dt)?;
    dev.boiling(BoilingOp::Cells, &b, &bo.k, dt)?;
    dev.copy(b.alpha_out, tf.gas.alpha)?;
    Ok(())
}

/// Step (3): the temperature equations of both phases over `dt` with the
/// final fluxes of the step.
pub(super) fn energy(fl: &mut Flow, tf: &mut TwoFluid, dt: f64) -> Result<()> {
    let Some(mut bo) = tf.boil.take() else {
        return Ok(());
    };
    let r = energy_inner(fl, tf, &mut bo, dt);
    tf.boil = Some(bo);
    r
}

fn energy_inner(fl: &mut Flow, tf: &TwoFluid, bo: &mut Boil, dt: f64) -> Result<()> {
    let w = fl.w;
    let mu_l = tf.k.mu_l;
    let k = bo.k;
    let lm = fl.f.mass_flux;
    let gm = tf.gas_f.mass_flux;
    let dev = fl.dev.as_mut();
    let phases: [(
        &Phase,
        Buf,
        Buf,
        Buf,
        Buf,
        Buf,
        Buf,
        (Buf, Buf),
        f64,
        f64,
        f64,
    ); 2] = [
        (
            &tf.liq, lm, bo.b.t_l, bo.t_l_old, bo.b.su_l, bo.b.sp_l, bo.b.cap_l, bo.bc_l, k.cp_l,
            k.k_l, 1.0,
        ),
        (
            &tf.gas,
            gm,
            bo.b.t_g,
            bo.t_g_old,
            bo.b.su_g,
            bo.b.sp_g,
            bo.b.cap_g,
            bo.bc_g,
            k.cp_g,
            k.k_g,
            k.rho_g / k.rho_l,
        ),
    ];
    for (i, (ph, m, t, t_old, su, sp, cap, (ba, bb), cp, kk, s)) in phases.into_iter().enumerate() {
        // the convecting heat-capacity flux c_p α_f ṁ and the conductance
        // α_f (k + c_p (μ_t s)/Pr_t)
        dev.product(bo.conv, cp, ph.alpha_f, m)?;
        dev.copy(tf.mu_eff, bo.gam)?;
        dev.subtract(bo.gam, mu_l)?;
        dev.lincomb(bo.gam, cp * s / bo.prt, bo.gam, 0.0, bo.gam)?;
        dev.subtract(bo.gam, -kk)?;
        dev.product(bo.gam, 1.0, ph.alpha_f, bo.gam)?;
        dev.scalar_faces(&ScalarFaceArgs {
            bc_a: ba,
            bc_b: bb,
            jump: None,
            mass_flux: bo.conv,
            gamma: bo.gam,
            phi: t,
            grad: bo.grad,
            psi: bo.psi,
            blend: 0.0,
            central: 0.0,
            beta: 1.0,
            nonorth: 0.0,
            bound: 0.0,
            cp: w.cp,
            cn: w.cn,
            bf: w.bf.x,
        })?;
        dev.scalar_rows(&ScalarRowArgs {
            cp: w.cp,
            cn: w.cn,
            bf: w.bf.x,
            mass_flux: bo.conv,
            phi: t,
            su: Some(su),
            sp: Some(sp),
            alpha: 1.0,
            time: Some(ScalarTime {
                inv_dt: 1.0 / dt,
                a: [1.0, -1.0, 0.0],
                old: t_old,
                old2: t_old,
                bounds: (f64::NEG_INFINITY, f64::INFINITY),
            }),
            fix: None,
            positive: false,
            pseudo: None,
            capacity: Some(cap),
            val: w.val,
            rhs: w.rhs.x,
            diag: w.diag.x,
            order: w.defect,
        })?;
        let (it, _) = dev.solve(
            w.val,
            w.rhs.x,
            t,
            foehncraft_fv::solver::Preconditioner::Jacobi,
            false,
            1e-10,
            1000,
        )?;
        bo.totals.energy_iterations[i] = it;
    }
    Ok(())
}

/// The host copies, the totals and the output fields of the boiling.
pub(super) fn download(fl: &mut Flow, tf: &mut TwoFluid) -> Result<()> {
    let Some(bo) = tf.boil.as_mut() else {
        return Ok(());
    };
    let k = bo.k;
    let ni = fl.disc.n_interior;
    let dev = fl.dev.as_mut();
    let tl = dev.download(bo.b.t_l)?;
    let tg = dev.download(bo.b.t_g)?;
    let gam = dev.download(bo.b.gamma)?;
    let gi = dev.download(bo.b.gamma_i)?;
    let diam = dev.download(bo.b.diam)?;
    let vol = &fl.disc.cell_volume;
    let n = bo.slot_face.len();
    let (tw, qe) = if n > 0 {
        let mut a = dev.download(bo.b.t_w)?;
        let mut b = dev.download(bo.b.q_e)?;
        a.truncate(n);
        b.truncate(n);
        (a, b)
    } else {
        (Vec::new(), Vec::new())
    };
    let (me, qg) = if n > 0 {
        let mut a = dev.download(bo.b.m_e)?;
        let mut b = dev.download(bo.b.q_g)?;
        a.truncate(n);
        b.truncate(n);
        (a, b)
    } else {
        (Vec::new(), Vec::new())
    };
    let t = &mut bo.totals;
    t.wall_heat = bo
        .slot_q
        .iter()
        .zip(&bo.slot_area)
        .map(|(q, a)| q * a)
        .exact_sum();
    t.wall_evaporation_heat = qe.iter().zip(&bo.slot_area).map(|(q, a)| q * a).exact_sum();
    t.wall_vapour = me.iter().copied().exact_sum();
    t.wall_vapour_heat = qg.iter().copied().exact_sum();
    t.interface_vapour = gi.iter().zip(vol).map(|(g, v)| g * v).exact_sum();
    t.max_wall_temperature = tw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    // the enthalpy through the boundaries
    let ml = dev.download(fl.f.mass_flux)?;
    let mg = dev.download(tf.gas_f.mass_flux)?;
    let alf = dev.download(tf.liq.alpha_f)?;
    let agf = dev.download(tf.gas.alpha_f)?;
    let owner = &fl.disc.owner;
    let mut parts = Vec::new();
    for f in ni..ml.len() {
        let c = owner[f] as usize;
        let tin = bo.t_in[f - ni];
        let (tlb, tgb) = if tin.is_finite() && ml[f] < 0.0 {
            (tin, k.t_sat)
        } else {
            (tl[c], tg[c])
        };
        parts.push(alf[f] * ml[f] * k.cp_l * tlb);
        let hg = k.cp_l * k.t_sat + k.h_fg + k.cp_g * (tgb - k.t_sat);
        parts.push(agf[f] * mg[f] * hg);
    }
    t.enthalpy_out = parts.into_iter().exact_sum();
    let rate: Vec<f64> = gam;
    fl.phys.boiling_fields = Some(std::sync::Arc::new(vec![
        ("liquid_temperature".to_string(), tl),
        ("vapour_temperature".to_string(), tg),
        ("bubble_diameter".to_string(), diam),
        ("vapour_generation".to_string(), rate),
    ]));
    Ok(())
}
