//! Conservation ledger (D88), balance reports (D89) and the kinetic-energy
//! form of the second-law check (D90 part 4).
//!
//! The ledger books, for each conserved quantity, the storage (the time
//! derivative of its integral over the domain), the boundary fluxes per
//! boundary and the volume sources, and closes it:
//!
//! imbalance = storage − Σ boundary fluxes into the domain − Σ sources.
//!
//! **Mass** (kg/s). The density is constant, so the storage is zero and
//! the imbalance is the net mass flow out of the domain. The interior face
//! fluxes cancel exactly, so it is the sum of the continuity residuals of
//! all cells: a check of the pressure-correction tolerance.
//!
//! **Kinetic energy** (W). The momentum equation dotted with u gives, for
//! constant ρ and μ and the Laplacian form of the viscous term that the
//! solver discretises,
//!
//! d/dt ∫ ½ρ|u|² dV = Σ_boundaries [ −∮ ½ρ|u|² u·dS − ∮ p u·dS + ∮ (μ∇u·dS)·u ]
//!                    + ∫ f·u dV − ∫ μ ∇u:∇u dV,
//!
//! boundary power (convected kinetic energy, pressure work, viscous work),
//! the power of the body force and the viscous dissipation Φ = μ ∇u:∇u,
//! which is never negative. The discrete solution satisfies this only up
//! to the dissipation of the schemes (upwind-biased convection, the
//! interpolations, the pressure–velocity coupling) and, in transient runs,
//! of the time integration. The ledger computes every physical term from
//! the solution and books the remainder as the **numerical dissipation**;
//! the second law asks that it is not negative (D90 part 4), and it falls
//! as the mesh is refined.
//!
//! The ledger splits it in two. The momentum rows hold the pressure as
//! the volume source −∇p_P V_P with the least-squares cell gradient, so
//! their energy budget has −Σ u_P·∇p_P V_P where the exact one has the
//! boundary pressure work; the difference, the discrete −∫ p ∇·u dV of the
//! cell velocities, is the **pressure part**. The momentum interpolation
//! of the face fluxes (Rhie–Chow) makes it a dissipation on regular cells;
//! on distorted cells it holds a consistency error of either sign that
//! falls with the mesh. The rest is the **scheme part**, of the convection
//! and time schemes: positive for upwind convection, of either sign for
//! the linear-upwind correction near singular points.
//!
//! **Mass in transient runs.** The discrete continuity equation holds at
//! the new time level, so the mass flows of a step are those at its end.
//!
//! In a transient run the storage of kinetic energy over a step is
//! (Kⁿ⁺¹ − Kⁿ)/Δt and its rates are the mean of their values at tⁿ and tⁿ⁺¹ (trapezoidal rule,
//! second order), so that the cumulative ledger of the run is
//! ΔK = ∫ (boundary power + body-force power − Φ − numerical dissipation) dt.
//!
//! All sums are pairwise sums over the cells or faces in a fixed order
//! (D31): the ledger is bitwise reproducible and does not depend on the
//! thread count or the device. It is computed from the fields on the host
//! and changes nothing of the solution.
//!
//! The volumes are passed per evaluation (D91): with a deforming mesh
//! (L22) each time level brings its own; today they are the volumes of the
//! fixed mesh.

use crate::case::{FlowCase, wall_face_values};
use crate::diag::{self, Diagnostic};
use crate::simplec::Solution;
use foehncraft_core::repro::pairwise_sum;
use foehncraft_kernels::flow::{MOM_FIXED, MOM_SYMMETRY};
use foehncraft_mesh::geometry::{dot, norm};
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// What an exceeded tolerance does.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum BalanceAction {
    /// A warning with a code (D78) in the log and the result.
    #[default]
    Warn,
    /// Stop the run as at a limit: restart file, final output, exit code 3.
    Stop,
}

/// `[balance]` (D89).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct BalanceSection {
    /// The full balance table in the log every N outer iterations (steady)
    /// or time steps (transient); 0: at the end only. A transient run also
    /// gives one compact balance line with every progress line.
    pub every: usize,
    /// Largest relative mass imbalance of one time step, or of the final
    /// state of a steady run (-).
    pub tolerance_step: f64,
    /// Largest relative mass imbalance over the run (-).
    pub tolerance_run: f64,
    /// Largest negative numerical dissipation of kinetic energy, relative
    /// to the largest term of the ledger, of the final steady state or over
    /// the run (-). Below it the second-law check warns; a transient step
    /// below it counts as negative. The default leaves room for the
    /// consistency errors of distorted cells, which fall with the mesh.
    pub tolerance_energy: f64,
    /// What an exceeded tolerance does: `warn`, or `stop` (a transient run
    /// stops at the step as at a limit, exit code 3; a steady run checks
    /// its final state only, so `stop` acts as `warn` there).
    pub action: BalanceAction,
}

impl Default for BalanceSection {
    fn default() -> Self {
        Self {
            every: 0,
            tolerance_step: 1e-8,
            tolerance_run: 1e-4,
            tolerance_energy: 1e-2,
            action: BalanceAction::Warn,
        }
    }
}

impl BalanceSection {
    /// Checks the ranges.
    pub fn validate(&self) -> anyhow::Result<()> {
        crate::need!(
            self.tolerance_step > 0.0,
            CASE_RANGE,
            "[balance] tolerance_step must be > 0"
        );
        crate::need!(
            self.tolerance_run > 0.0,
            CASE_RANGE,
            "[balance] tolerance_run must be > 0"
        );
        crate::need!(
            self.tolerance_energy > 0.0,
            CASE_RANGE,
            "[balance] tolerance_energy must be > 0"
        );
        Ok(())
    }
}

/// The boundary terms of one boundary.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct BoundaryTerms {
    /// Boundary name.
    pub boundary: String,
    /// Boundary condition.
    pub kind: String,
    /// Mass flow out of the domain (kg/s).
    pub mass_flow: f64,
    /// Kinetic energy carried out of the domain, Σ ṁ ½|u_b|² (W).
    pub convective_power: f64,
    /// Work of the pressure on the outflow, Σ p_b u_b·S (W).
    pub pressure_power: f64,
    /// Work of the viscous stress on the fluid, Σ (μ ∇u·S)·u_b (W).
    pub viscous_power: f64,
}

impl BoundaryTerms {
    /// Power into the fluid: viscous work − convected kinetic energy −
    /// pressure work (W).
    pub fn power(&self) -> f64 {
        self.viscous_power - self.convective_power - self.pressure_power
    }
}

/// The instantaneous terms of one state.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Rates {
    /// Kinetic energy in the domain, ∫ ½ρ|u|² dV (J).
    pub kinetic_energy: f64,
    /// Boundary terms per boundary.
    pub boundaries: Vec<BoundaryTerms>,
    /// Power of the body force, ∫ f·u dV (W).
    pub body_force_power: f64,
    /// Power of the pressure source of the momentum rows on the cell
    /// velocities, −Σ u_P·∇p_P V_P (W). With the boundary pressure work it
    /// gives the pressure part of the numerical dissipation.
    pub pressure_source_power: f64,
    /// Viscous dissipation: of the discrete viscous operator, which
    /// converges to ∫ μ ∇u:∇u dV (W).
    pub viscous_dissipation: f64,
    /// Scale of the mass flows: the largest boundary mass flow, or in a closed
    /// domain ρ |u|max times the largest boundary area (kg/s).
    pub mass_scale: f64,
    /// Torque and shaft power of the rotors, machine values (D184).
    pub rotors: crate::rotor::Values,
    /// The mass in the domain, Σ ρ V (kg), with a variable density (stage
    /// 7c; its storage enters the mass account of a time step).
    pub mass: Option<f64>,
    /// The pressure dilatation ∫ p ∇·u dV of a variable density (W, stage
    /// 7c): the physical exchange between kinetic and internal energy,
    /// with ∇·u from the face fluxes, Σ_f ṁ_f/ρ_f per cell.
    pub pressure_dilatation: Option<f64>,
}

/// The shaft work of one rotor (D184): −Ω·T, the work of its walls on the
/// fluid when they turn with the zone, booked in the wall terms of its walls.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct ShaftTerms {
    /// The rotor (its zone).
    pub rotor: String,
    /// Its walls.
    pub walls: Vec<String>,
    /// Shaft power (W), or shaft work over a run (J).
    pub power: f64,
    /// Of it the work of the transposed viscous stress, which is not in the
    /// wall terms (the ledger has the Laplacian form of the viscous term).
    pub transposed: f64,
}

impl ShaftTerms {
    fn of(r: &Rates) -> Vec<Self> {
        r.rotors
            .rotors
            .iter()
            .map(|v| Self {
                rotor: v.name.clone(),
                walls: v.walls.clone(),
                power: v.power,
                transposed: v.transposed_work,
            })
            .collect()
    }
}

/// The terms of a state on the host. `volume` holds the cell volumes of
/// the time level of the state (D91); `nonorth` is the non-orthogonal
/// correction of the solver, which the viscous boundary flux follows.
pub fn rates(
    sol: &Solution,
    mesh: &Mesh,
    geo: &Geometry,
    case: &FlowCase,
    mom_kind: &[u32],
    volume: &[f64],
    nonorth: bool,
) -> Rates {
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let rho = sol.phys.rho;
    let f = sol.phys.body_force;
    let mut ke = vec![0.0; nc];
    let mut body = vec![0.0; nc];
    let mut press = vec![0.0; nc];
    let mut umax = 0.0f64;
    for c in 0..nc {
        let u = sol.u.get(c);
        let v = volume[c];
        let uu = dot(u, u);
        umax = umax.max(uu.sqrt());
        ke[c] = 0.5 * sol.phys.rho_at(c) * uu * v;
        body[c] = match &sol.phys.source {
            Some(src) => dot([src[0][c], src[1][c], src[2][c]], u) * v,
            None => dot(f, u) * v,
        };
        // the work of the buoyancy force of the momentum rows (stage 7b1)
        if let Some(bf) = &sol.phys.buoyancy {
            body[c] += dot([bf[0][c], bf[1][c], bf[2][c]], u) * v;
        }
        // the pressure source of the momentum rows, −∇p_P V_P
        // (foehncraft_kernels::flow::momentum_row)
        press[c] = -dot(sol.grad_p.get(c), u) * v;
    }
    let d = sol.disc;
    let on = if nonorth { 1.0 } else { 0.0 };
    // the dissipation of the discrete viscous operator, face by face:
    // F_f·(u_N − u_P) with the diffusive flux F_f of the momentum rows
    // (foehncraft_kernels::flow::momentum_face, non-orthogonal and
    // skewness corrections included); it converges to ∫ μ ∇u:∇u dV and
    // is what the discrete momentum equation dotted with u removes
    let mut phi = vec![0.0; d.n_faces];
    for (fc, ph) in phi.iter_mut().enumerate().take(ni) {
        let (p, n) = (d.owner[fc] as usize, d.neighbour[fc] as usize);
        let (up, un) = (sol.u.get(p), sol.u.get(n));
        let w = d.weight[fc];
        let xf = d.face_centroid.get(fc);
        let (xp, xn) = (d.cell_centroid.get(p), d.neighbour_centroid(fc));
        let rp = [xf[0] - xp[0], xf[1] - xp[1], xf[2] - xp[2]];
        let rm = [
            xf[0] - 0.5 * (xn[0] + xf[0] - rp[0]),
            xf[1] - 0.5 * (xn[1] + xf[1] - rp[1]),
            xf[2] - 0.5 * (xn[2] + xf[2] - rp[2]),
        ];
        let k = d.k.get(fc);
        let mu = sol.mu_eff_at(fc);
        let dd = mu * d.dgeo[fc];
        let mut e = 0.0;
        for i in 0..3 {
            let (gp, gn) = (sol.grad_u[i].get(p), sol.grad_u[i].get(n));
            let gf = [
                w * gp[0] + (1.0 - w) * gn[0],
                w * gp[1] + (1.0 - w) * gn[1],
                w * gp[2] + (1.0 - w) * gn[2],
            ];
            let skew = dot([gn[0] - gp[0], gn[1] - gp[1], gn[2] - gp[2]], rm);
            let flux = dd * (un[i] - up[i]) + on * (mu * dot(k, gf) + dd * skew);
            e += flux * (un[i] - up[i]);
        }
        // a periodic pair carries the face twice, once from each side
        *ph = if fc >= d.periodic { 0.5 * e } else { e };
    }
    let mut boundaries = Vec::with_capacity(mesh.boundaries.len());
    let mut scale = 0.0f64;
    // in a turning zone the kinetic energy is carried by the convecting
    // flux, relative to the zone (zero on a wall that turns with it)
    let turning = case.turning();
    for p in &mesh.boundaries {
        let bc = &case.boundary[&p.name];
        let wall = matches!(bc, crate::bc::BoundaryCondition::Wall { .. });
        let n = p.len;
        let (mut m, mut kc, mut pw, mut vw, mut area) = (
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
        );
        for fc in p.range() {
            let fb = fc - ni;
            let s = geo.face_area[fc];
            let ub = sol.ub.get(fb);
            let mf = sol.mass_flux[fc];
            m.push(mf);
            area.push(norm(s));
            let conv = if turning.is_empty() {
                mf
            } else {
                mf - crate::rotor::frame_flux(&turning, mesh, geo, fc, rho)
            };
            kc.push(conv * 0.5 * dot(ub, ub));
            let pf = if wall {
                wall_face_values(sol, geo, fc, ni).0
            } else {
                sol.pb[fb]
            };
            pw.push(pf * dot(ub, s));
            let kind = mom_kind[fb];
            let c = d.owner[fc] as usize;
            let up = sol.u.get(c);
            let mu = sol.mu_eff_at(fc);
            let dd = mu * d.dgeo[fc];
            if kind == MOM_FIXED {
                // the diffusive flux of the momentum rows
                // (foehncraft_kernels::flow::momentum_face): work on the
                // fluid F·u_b, dissipation F·(u_b − u_P)
                let k = d.k.get(fc);
                let (mut w, mut e) = (0.0, 0.0);
                for i in 0..3 {
                    let flux = dd * (ub[i] - up[i]) + on * mu * dot(k, sol.grad_u[i].get(c));
                    w += flux * ub[i];
                    e += flux * (ub[i] - up[i]);
                }
                vw.push(w);
                phi[fc] = e;
            } else if kind == MOM_SYMMETRY {
                // the normal velocity goes to zero at the face: no work,
                // dissipation μ d (u_P·n̂)²
                let un = dot(up, s) / norm(s);
                vw.push(0.0);
                phi[fc] = dd * un * un;
            } else {
                // outflow: zero normal gradient, no viscous flux
                vw.push(0.0);
            }
        }
        let t = BoundaryTerms {
            boundary: p.name.clone(),
            kind: bc.type_name().into(),
            mass_flow: pairwise_sum(&m),
            convective_power: pairwise_sum(&kc),
            pressure_power: pairwise_sum(&pw),
            viscous_power: pairwise_sum(&vw),
        };
        scale = scale
            .max(t.mass_flow.abs())
            .max(rho * umax * pairwise_sum(&area));
        boundaries.push(t);
    }
    // the variable density (stage 7c): the mass in the domain and the
    // pressure dilatation Σ_P p_P Σ_f ṁ_f/ρ_f
    let (mass, pressure_dilatation) = match &sol.phys.rho_var {
        Some(r) => {
            let m: Vec<f64> = (0..nc).map(|c| r.0[c] * volume[c]).collect();
            let mut div = vec![0.0; nc];
            for fc in 0..d.n_faces {
                let q = sol.mass_flux[fc] / r.1[fc];
                div[d.owner[fc] as usize] += q;
                if fc < ni {
                    div[d.neighbour[fc] as usize] -= q;
                }
            }
            let pd: Vec<f64> = (0..nc).map(|c| sol.p[c] * div[c]).collect();
            (Some(pairwise_sum(&m)), Some(pairwise_sum(&pd)))
        }
        None => (None, None),
    };
    Rates {
        mass,
        pressure_dilatation,
        kinetic_energy: pairwise_sum(&ke),
        boundaries,
        body_force_power: pairwise_sum(&body),
        pressure_source_power: pairwise_sum(&press),
        viscous_dissipation: pairwise_sum(&phi),
        mass_scale: scale,
        rotors: crate::rotor::evaluate(sol, mesh, geo, case),
    }
}

/// What a term of an account is (D119). Every kind but [`TermKind::Part`]
/// is a contribution to the quantity and enters the closure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum TermKind {
    /// A flow through a boundary (convective, diffusive, work).
    Boundary,
    /// A named volume source (positive) or sink (negative).
    Source,
    /// A physical dissipation or other loss to another quantity
    /// (negative), booked as a pair when the partner is in the ledger.
    Dissipation,
    /// One side of an exchange pair between two accounts: the same face
    /// or volume sum with opposite signs (D88, D117 (1)).
    Exchange,
    /// A breakdown of other terms, reported and not closed (the pressure
    /// source of the momentum rows, whose difference to the boundary
    /// pressure work is the pressure part of the numerical dissipation).
    Part,
}

/// How a term is booked over a time step.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum StepRule {
    /// The mean of its values at the two ends of the step (trapezoidal
    /// rule, second order).
    #[default]
    Trapezoid,
    /// Its value at the end of the step: a term of an equation that holds
    /// at the new time level (the mass flows of the discrete continuity
    /// equation).
    End,
}

/// One term of an account (D119): `{quantity, unit, region, term kind,
/// name, exchange-pair id}`, the quantity, unit and region being those of
/// its [`Account`]. The value is the contribution to the quantity: a
/// rate (per second) for a state, an integral over a run.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Term {
    /// Name of the term ("outflow", "viscous work", "body force", …).
    pub name: String,
    /// Its kind.
    pub kind: TermKind,
    /// The boundary of a boundary term; the terms of one boundary are
    /// summed together, in their order, before the boundaries are summed.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub boundary: Option<String>,
    /// The id of an exchange pair: the other account has the term with the
    /// same id and the opposite sign.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pair: Option<u32>,
    /// How it is booked over a time step.
    #[serde(skip)]
    pub rule: StepRule,
    /// Its contribution to the quantity (unit of the account per second,
    /// or the unit over a run).
    pub value: f64,
}

impl Term {
    /// A term of a boundary.
    pub fn boundary(name: &str, boundary: &str, value: f64) -> Self {
        Self {
            name: name.into(),
            kind: TermKind::Boundary,
            boundary: Some(boundary.into()),
            pair: None,
            rule: StepRule::Trapezoid,
            value,
        }
    }

    /// A volume term of the given kind.
    pub fn volume(name: &str, kind: TermKind, value: f64) -> Self {
        Self {
            name: name.into(),
            kind,
            boundary: None,
            pair: None,
            rule: StepRule::Trapezoid,
            value,
        }
    }

    /// The same term with another step rule.
    pub fn rule(mut self, rule: StepRule) -> Self {
        self.rule = rule;
        self
    }

    /// The same term as one side of the exchange pair `id`.
    pub fn paired(mut self, id: u32) -> Self {
        self.pair = Some(id);
        self
    }
}

/// The ledger of one conserved quantity in one region (D88, D119): its
/// content and its terms, with one closure routine ([`Account::close`])
/// for every quantity a module books.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct Account {
    /// The quantity ("mass", "kinetic energy", later "thermal energy",
    /// species, charge, …).
    pub quantity: String,
    /// Unit of the content (kg, J, …); the terms are this unit per second,
    /// or this unit over a run.
    pub unit: String,
    /// The region ("domain" for the whole mesh).
    pub region: String,
    /// The content at the end (the integral of the quantity), if the
    /// account knows it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub content: Option<f64>,
    /// The scale the imbalance is measured against (for mass: the largest
    /// mass flow), if the account has one; otherwise the largest term.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scale: Option<f64>,
    /// The terms in their booking order.
    pub terms: Vec<Term>,
}

/// The closure of an account: imbalance = storage − Σ contributions.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize)]
pub struct Closure {
    /// Storage: the change of the content per time (or over the run).
    pub storage: f64,
    /// Σ of the boundary terms: per boundary in the order of its terms,
    /// then a pairwise sum over the boundaries (D31).
    pub boundary: f64,
    /// Σ of all contributions: the boundary sum, then the volume terms
    /// added one after the other in their order.
    pub contributions: f64,
    /// storage − contributions.
    pub imbalance: f64,
    /// The largest magnitude of the storage and of every term (parts
    /// included).
    pub largest: f64,
    /// imbalance over the scale, or over the largest term (-).
    pub relative: f64,
}

impl Account {
    /// A new account of the whole domain.
    pub fn new(quantity: &str, unit: &str) -> Self {
        Self {
            quantity: quantity.into(),
            unit: unit.into(),
            region: "domain".into(),
            ..Default::default()
        }
    }

    /// The first term of this name (of the boundary, if one is given).
    pub fn term(&self, name: &str, boundary: Option<&str>) -> Option<&Term> {
        self.terms
            .iter()
            .find(|t| t.name == name && t.boundary.as_deref() == boundary)
    }

    /// The value of a term, zero if the account has none of this name.
    pub fn value(&self, name: &str, boundary: Option<&str>) -> f64 {
        self.term(name, boundary).map_or(0.0, |t| t.value)
    }

    /// The boundaries in the order of their first term.
    pub fn boundaries(&self) -> Vec<&str> {
        let mut out: Vec<&str> = Vec::new();
        for t in &self.terms {
            if let Some(b) = t.boundary.as_deref()
                && !out.contains(&b)
            {
                out.push(b);
            }
        }
        out
    }

    /// The sum of the terms of each boundary (in the order of the terms),
    /// in the order of [`Account::boundaries`].
    pub fn boundary_sums(&self) -> Vec<f64> {
        self.boundaries()
            .iter()
            .map(|b| {
                let mut s = 0.0;
                let mut first = true;
                for t in self
                    .terms
                    .iter()
                    .filter(|t| t.kind != TermKind::Part && t.boundary.as_deref() == Some(*b))
                {
                    s = if first { t.value } else { s + t.value };
                    first = false;
                }
                s
            })
            .collect()
    }

    /// The one closure routine of the ledger (D119): the boundary terms
    /// summed per boundary and pairwise over the boundaries, the volume
    /// terms added in their order, imbalance = storage − contributions.
    /// Every sum has a fixed order, so the closure is bitwise
    /// reproducible (D31).
    pub fn close(&self, storage: f64) -> Closure {
        let boundary = pairwise_sum(&self.boundary_sums());
        let mut contributions = boundary;
        for t in self
            .terms
            .iter()
            .filter(|t| t.kind != TermKind::Part && t.boundary.is_none())
        {
            contributions += t.value;
        }
        let imbalance = storage - contributions;
        let mut largest = storage.abs();
        for t in &self.terms {
            largest = largest.max(t.value.abs());
        }
        let denom = self.scale.unwrap_or(largest);
        Closure {
            storage,
            boundary,
            contributions,
            imbalance,
            largest,
            relative: if denom > 0.0 { imbalance / denom } else { 0.0 },
        }
    }

    /// The account with every value (terms, scale) multiplied by `a`.
    fn scaled(&self, a: f64) -> Self {
        let mut o = self.clone();
        for t in &mut o.terms {
            t.value *= a;
        }
        o.scale = o.scale.map(|s| a * s);
        o
    }

    /// The account of a time step from the accounts at its two ends: each
    /// term by its [`StepRule`], the content and the scale of the end.
    fn over_step(a: &Self, b: &Self) -> Self {
        let mut o = b.clone();
        for (t, ta) in o.terms.iter_mut().zip(&a.terms) {
            if t.rule == StepRule::Trapezoid {
                t.value = 0.5 * ta.value + 0.5 * t.value;
            }
        }
        o
    }

    /// Adds `dt` times the step account `s` to this (run) account.
    fn accumulate(&mut self, s: &Self, dt: f64) {
        for (t, ts) in self.terms.iter_mut().zip(&s.terms) {
            t.value += dt * ts.value;
        }
        if let (Some(acc), Some(v)) = (self.scale.as_mut(), s.scale) {
            *acc += dt * v;
        }
        self.content = s.content;
    }
}

/// The sum of the two sides of every exchange pair over a set of accounts:
/// each is zero to the bit when the two modules book the same face or volume
/// sum with opposite signs (D88, D117 (1)). Pairs whose partner is missing
/// are returned with the one value they have.
pub fn exchange_residuals(accounts: &[Account]) -> Vec<(u32, f64)> {
    let mut out: Vec<(u32, f64)> = Vec::new();
    for a in accounts {
        for t in &a.terms {
            if let Some(id) = t.pair {
                match out.iter_mut().find(|(i, _)| *i == id) {
                    Some((_, v)) => *v += t.value,
                    None => out.push((id, t.value)),
                }
            }
        }
    }
    out
}

/// Names of the mass and kinetic-energy terms.
mod names {
    pub const MASS: &str = "mass";
    pub const KINETIC: &str = "kinetic energy";
    pub const OUTFLOW: &str = "outflow";
    pub const VISCOUS_WORK: &str = "viscous work";
    pub const CONVECTED: &str = "convected";
    pub const PRESSURE_WORK: &str = "pressure work";
    pub const BODY_FORCE: &str = "body force";
    pub const VISCOUS_DISSIPATION: &str = "viscous dissipation";
    pub const PRESSURE_SOURCE: &str = "pressure source";
    pub const PRESSURE_DILATATION: &str = "pressure dilatation";
}

impl Rates {
    /// The terms of the state as accounts (D119): mass and kinetic energy
    /// of the domain. Every term is a contribution to the quantity, so the
    /// outflows and the dissipation carry a minus sign; negating is exact,
    /// so the closures give the numbers of the fixed terms bitwise.
    pub fn accounts(&self) -> Vec<Account> {
        use names::*;
        let mut mass = Account::new(MASS, "kg");
        mass.scale = Some(self.mass_scale);
        let mut ke = Account::new(KINETIC, "J");
        ke.content = Some(self.kinetic_energy);
        for p in &self.boundaries {
            let b = p.boundary.as_str();
            mass.terms
                .push(Term::boundary(OUTFLOW, b, -p.mass_flow).rule(StepRule::End));
            ke.terms
                .push(Term::boundary(VISCOUS_WORK, b, p.viscous_power));
            ke.terms
                .push(Term::boundary(CONVECTED, b, -p.convective_power));
            ke.terms
                .push(Term::boundary(PRESSURE_WORK, b, -p.pressure_power));
        }
        ke.terms.push(Term::volume(
            BODY_FORCE,
            TermKind::Source,
            self.body_force_power,
        ));
        ke.terms.push(Term::volume(
            VISCOUS_DISSIPATION,
            TermKind::Dissipation,
            -self.viscous_dissipation,
        ));
        ke.terms.push(Term::volume(
            PRESSURE_SOURCE,
            TermKind::Part,
            self.pressure_source_power,
        ));
        if let Some(pd) = self.pressure_dilatation {
            ke.terms
                .push(Term::volume(PRESSURE_DILATATION, TermKind::Exchange, pd));
        }
        mass.content = self.mass;
        vec![mass, ke]
    }
}

/// Mass ledger (kg/s, or kg over a run).
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct MassBalance {
    /// Storage (zero for constant density).
    pub storage: f64,
    /// Net mass flow out of the domain.
    pub outflow: f64,
    /// storage − inflow = storage + outflow.
    pub imbalance: f64,
    /// Imbalance over the scale of the mass flows (-).
    pub relative: f64,
}

/// Kinetic-energy ledger (W, or J over a run).
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct EnergyBalance {
    /// Kinetic energy in the domain at the end (J).
    pub energy: f64,
    /// Storage: its change per time (W) or over the run (J).
    pub storage: f64,
    /// Power of all boundaries into the fluid.
    pub boundary: f64,
    /// Of it the shaft power of the rotors, −Ω·T (D184).
    pub shaft: f64,
    /// Power of the body force.
    pub body_force: f64,
    /// Viscous dissipation (physical, ≥ 0).
    pub viscous_dissipation: f64,
    /// Numerical dissipation: boundary + body force − viscous dissipation −
    /// storage; ≥ 0 by the second law (D90 part 4). It is the sum of the
    /// two parts below.
    pub numerical_dissipation: f64,
    /// Its pressure part: the boundary pressure work minus the power of the
    /// pressure source −∇p_P V_P of the momentum rows on the cell
    /// velocities, the discrete −∫ p ∇·u dV of the cell velocities. The
    /// momentum interpolation of the face fluxes makes it a dissipation;
    /// it is not sign-definite on distorted cells.
    pub pressure_dissipation: f64,
    /// Its scheme part: of the convection and time schemes (the rest).
    pub scheme_dissipation: f64,
    /// Numerical dissipation over the largest term (-).
    pub relative: f64,
}

/// One closed ledger: of a state (steady), a time step or a whole run.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct Balance {
    /// Outer iteration (steady).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub iteration: Option<usize>,
    /// Time step (transient; the last one for a run).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub step: Option<usize>,
    /// Time at the end (s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub time: Option<f64>,
    /// Length of the step or of the run (s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub duration: Option<f64>,
    /// Mass.
    pub mass: MassBalance,
    /// Kinetic energy.
    pub kinetic_energy: EnergyBalance,
    /// Boundary terms per boundary (rates; integrals over a run).
    pub boundaries: Vec<BoundaryTerms>,
    /// The shaft work per rotor, in the terms of its walls (D184).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub shafts: Vec<ShaftTerms>,
    /// The accounts the balance is closed from (D119): every quantity
    /// with its terms. Not yet written to `result.json` (the fixed fields
    /// above are); later quantities report through them.
    #[serde(skip)]
    pub accounts: Vec<Account>,
}

/// The balance of a set of accounts: the closure of each, and the fixed
/// views of mass and kinetic energy (with their storages) that the
/// reports print.
fn close(
    accounts: Vec<Account>,
    storage_ke: f64,
    storage_mass: f64,
    shafts: &[ShaftTerms],
) -> (MassBalance, EnergyBalance, Vec<BoundaryTerms>, Vec<Account>) {
    use names::*;
    let find = |q: &str| {
        accounts
            .iter()
            .find(|a| a.quantity == q)
            .expect("the mass and kinetic-energy accounts are booked")
    };
    let (mass, ke) = (find(MASS), find(KINETIC));
    let ck = ke.close(storage_ke);
    // the boundary terms per boundary, in the sign of the reports
    let boundaries: Vec<BoundaryTerms> = ke
        .boundaries()
        .into_iter()
        .map(|b| BoundaryTerms {
            boundary: b.into(),
            kind: String::new(),
            mass_flow: -mass.value(OUTFLOW, Some(b)),
            convective_power: -ke.value(CONVECTED, Some(b)),
            pressure_power: -ke.value(PRESSURE_WORK, Some(b)),
            viscous_power: ke.value(VISCOUS_WORK, Some(b)),
        })
        .collect();
    let wall_pressure = pairwise_sum(
        &boundaries
            .iter()
            .map(|p| p.pressure_power)
            .collect::<Vec<_>>(),
    );
    let source = ke.value(PRESSURE_SOURCE, None);
    let pressure = match ke.term(PRESSURE_DILATATION, None) {
        // the physical part of the discrete −∫ p ∇·u dV is not numerical
        // (a variable density, stage 7c)
        Some(t) => -(source + wall_pressure) + t.value,
        None => -(source + wall_pressure),
    };
    // the views keep the signs of the reports, signed zeros included: the
    // outflow and the numerical dissipation (= −imbalance of the closure)
    // are formed from the terms in the sign the tables print
    let out = pairwise_sum(&boundaries.iter().map(|p| p.mass_flow).collect::<Vec<_>>());
    let num = ck.contributions - storage_ke;
    let scale = mass.scale.unwrap_or(0.0);
    (
        {
            // storage + outflow (no storage: the outflow, signed zero kept)
            let imb = if storage_mass == 0.0 {
                out
            } else {
                storage_mass + out
            };
            MassBalance {
                storage: storage_mass,
                outflow: out,
                imbalance: imb,
                relative: if scale > 0.0 { imb / scale } else { 0.0 },
            }
        },
        EnergyBalance {
            energy: ke.content.unwrap_or(0.0),
            storage: storage_ke,
            boundary: ck.boundary,
            shaft: pairwise_sum(&shafts.iter().map(|t| t.power).collect::<Vec<_>>()),
            body_force: ke.value(BODY_FORCE, None),
            viscous_dissipation: -ke.value(VISCOUS_DISSIPATION, None),
            numerical_dissipation: num,
            pressure_dissipation: pressure,
            scheme_dissipation: num - pressure,
            relative: if ck.largest > 0.0 {
                num / ck.largest
            } else {
                0.0
            },
        },
        boundaries,
        accounts,
    )
}

/// The ledger of a set of accounts with the boundary kinds of the rates.
fn balance_of(
    accounts: Vec<Account>,
    storage_ke: f64,
    storage_mass: f64,
    shafts: Vec<ShaftTerms>,
    kinds: &[BoundaryTerms],
) -> Balance {
    let (mass, kinetic_energy, mut boundaries, accounts) =
        close(accounts, storage_ke, storage_mass, &shafts);
    for (b, k) in boundaries.iter_mut().zip(kinds) {
        b.kind = k.kind.clone();
    }
    Balance {
        mass,
        kinetic_energy,
        boundaries,
        shafts,
        accounts,
        ..Default::default()
    }
}

/// The ledger of a steady state (no storage).
pub fn steady(r: &Rates, iteration: usize) -> Balance {
    Balance {
        iteration: Some(iteration),
        ..balance_of(r.accounts(), 0.0, 0.0, ShaftTerms::of(r), &r.boundaries)
    }
}

/// The ledger of a transient run: one balance per step and the totals.
#[derive(Debug, Clone, Default)]
pub struct Ledger {
    prev: Option<Rates>,
    start_energy: f64,
    start_time: f64,
    time: f64,
    step: usize,
    /// The accounts integrated over the run.
    run: Vec<Account>,
    shafts: Vec<ShaftTerms>,
    kinds: Vec<BoundaryTerms>,
    energy: f64,
    /// The mass in the domain at the start and now (a variable density).
    start_mass: Option<f64>,
    mass: Option<f64>,
    /// A step counts as negative below −tolerance (relative).
    pub tolerance: f64,
    /// Steps booked.
    pub steps: usize,
    /// Steps with a negative numerical dissipation (below −tolerance).
    pub negative_steps: usize,
    /// Largest relative mass imbalance of a step.
    pub largest_mass_relative: f64,
}

impl Ledger {
    /// Starts at the state of time `t` (start fields or a restart). A step
    /// counts as negative when its relative numerical dissipation is below
    /// −`tolerance` (round-off is not a finding).
    pub fn start(r: Rates, t: f64, step: usize, tolerance: f64) -> Self {
        Self {
            start_mass: r.mass,
            mass: r.mass,
            tolerance,
            start_energy: r.kinetic_energy,
            energy: r.kinetic_energy,
            start_time: t,
            time: t,
            step,
            run: r.accounts().iter().map(|a| a.scaled(0.0)).collect(),
            shafts: ShaftTerms::of(&r)
                .into_iter()
                .map(|t| ShaftTerms {
                    power: 0.0,
                    transposed: 0.0,
                    ..t
                })
                .collect(),
            kinds: r.boundaries.clone(),
            prev: Some(r),
            ..Default::default()
        }
    }

    /// Books the step to the state `r` at time `t`.
    pub fn step(&mut self, r: Rates, step: usize, t: f64) -> Balance {
        let prev = self.prev.take().unwrap_or_else(|| r.clone());
        let dt = t - self.time;
        let half = 0.5;
        let accounts: Vec<Account> = prev
            .accounts()
            .iter()
            .zip(&r.accounts())
            .map(|(a, b)| Account::over_step(a, b))
            .collect();
        let shafts: Vec<ShaftTerms> = ShaftTerms::of(&prev)
            .into_iter()
            .zip(ShaftTerms::of(&r))
            .map(|(a, b)| ShaftTerms {
                power: half * (a.power + b.power),
                transposed: half * (a.transposed + b.transposed),
                ..b
            })
            .collect();
        let storage = if dt > 0.0 {
            (r.kinetic_energy - prev.kinetic_energy) / dt
        } else {
            0.0
        };
        let storage_mass = match (r.mass, prev.mass) {
            (Some(b), Some(a)) if dt > 0.0 => (b - a) / dt,
            _ => 0.0,
        };
        // a variable density (stage 7c): the mass is known to the round-off
        // of its sum, a floor of the scale where nothing flows (a closed box
        // at rest)
        let mut accounts = accounts;
        if let (Some(m), true) = (r.mass, dt > 0.0)
            && let Some(a) = accounts.iter_mut().find(|a| a.quantity == names::MASS)
        {
            a.scale = a.scale.map(|x| x.max(1e3 * f64::EPSILON * m.abs() / dt));
        }
        for (acc, s) in self.run.iter_mut().zip(&accounts) {
            acc.accumulate(s, dt);
        }
        for (acc, t) in self.shafts.iter_mut().zip(&shafts) {
            acc.power += dt * t.power;
            acc.transposed += dt * t.transposed;
        }
        let b = balance_of(accounts, storage, storage_mass, shafts, &r.boundaries);
        self.mass = r.mass;
        self.energy = r.kinetic_energy;
        self.time = t;
        self.step = step;
        self.steps += 1;
        if b.kinetic_energy.relative < -self.tolerance {
            self.negative_steps += 1;
        }
        self.largest_mass_relative = self.largest_mass_relative.max(b.mass.relative.abs());
        self.kinds = r.boundaries.clone();
        self.prev = Some(r);
        Balance {
            step: Some(step),
            time: Some(t),
            duration: Some(dt),
            ..b
        }
    }

    /// The ledger of the run so far: integrals over time (kg, J).
    pub fn run(&self) -> Balance {
        let mut accounts = self.run.clone();
        for a in &mut accounts {
            if a.quantity == names::KINETIC {
                a.content = Some(self.energy);
            }
        }
        Balance {
            step: Some(self.step),
            time: Some(self.time),
            duration: Some(self.time - self.start_time),
            ..balance_of(
                accounts,
                self.energy - self.start_energy,
                match (self.mass, self.start_mass) {
                    (Some(b), Some(a)) => b - a,
                    _ => 0.0,
                },
                self.shafts.clone(),
                &self.kinds,
            )
        }
    }
}

/// The ledger of a solve in the report and in `result.json`.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct BalanceSummary {
    /// The last balance: of the final state (steady) or of the last step
    /// (transient), in kg/s and W.
    pub last: Balance,
    /// Transient: the whole run, in kg and J.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run: Option<Balance>,
    /// Transient: steps booked.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub steps: Option<usize>,
    /// Transient: steps whose numerical dissipation is negative (a health
    /// number: anti-diffusive steps, D90).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub negative_steps: Option<usize>,
    /// Largest relative mass imbalance of a step or of the final state.
    pub largest_mass_relative: f64,
}

/// The tolerance checks of D89 on a finished ledger: the warnings with
/// their codes.
pub fn check(s: &BalanceSummary, set: &BalanceSection) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    if s.largest_mass_relative > set.tolerance_step {
        out.push(
            Diagnostic::new(
                diag::BAL_MASS_STEP,
                format!(
                    "the relative mass imbalance reached {:.2e}, above [balance] tolerance_step {:.1e} \
                     (tighten [solver.pressure] tolerance)",
                    s.largest_mass_relative, set.tolerance_step
                ),
            )
            .key("balance.tolerance_step"),
        );
    }
    let whole = s.run.as_ref().unwrap_or(&s.last);
    if s.run.is_some() && whole.mass.relative.abs() > set.tolerance_run {
        out.push(
            Diagnostic::new(
                diag::BAL_MASS_RUN,
                format!(
                    "the mass imbalance of the run is {:.2e} of the mass flows, above [balance] \
                     tolerance_run {:.1e}",
                    whole.mass.relative, set.tolerance_run
                ),
            )
            .key("balance.tolerance_run"),
        );
    }
    if whole.kinetic_energy.relative < -set.tolerance_energy {
        out.push(
            Diagnostic::new(
                diag::BAL_NEGATIVE_DISSIPATION,
                format!(
                    "the numerical dissipation of kinetic energy is negative: {:.3e} {u} \
                     ({:.2e} of the largest term, below [balance] tolerance_energy {:.1e}; \
                     pressure part {:+.3e} {u}, schemes {:+.3e} {u}); the discretisation adds \
                     energy (D90)",
                    whole.kinetic_energy.numerical_dissipation,
                    whole.kinetic_energy.relative,
                    set.tolerance_energy,
                    whole.kinetic_energy.pressure_dissipation,
                    whole.kinetic_energy.scheme_dissipation,
                    u = if s.run.is_some() { "J" } else { "W" },
                ),
            )
            .key("balance.tolerance_energy"),
        );
    }
    out
}

/// One compact line for the progress log (transient).
pub fn compact(b: &Balance) -> String {
    let shaft = if b.shafts.is_empty() {
        String::new()
    } else {
        format!(" (shaft {:+.3e} W)", b.kinetic_energy.shaft)
    };
    format!(
        "      balance  mass {:+.2e}  kinetic energy: dK/dt {:+.3e} W, boundary {:+.3e} W{shaft}, \
         dissipation {:.3e} W, numerical {:+.3e} W ({:+.1e})",
        b.mass.relative,
        b.kinetic_energy.storage,
        b.kinetic_energy.boundary + b.kinetic_energy.body_force,
        b.kinetic_energy.viscous_dissipation,
        b.kinetic_energy.numerical_dissipation,
        b.kinetic_energy.relative
    )
}

/// The full table (log, report): per quantity, boundary and term.
pub fn table(title: &str, b: &Balance, integrated: bool) -> String {
    let (m, e) = if integrated {
        ("kg", "J")
    } else {
        ("kg/s", "W")
    };
    // the rows are indented as the title is
    let ind = &title[..title.len() - title.trim_start().len()];
    let mut s = format!("{title}\n");
    let row = |s: &mut String, name: &str, v: f64, unit: &str| {
        s.push_str(&format!("{ind}  {name:<44} {v:>+14.6e} {unit}\n"));
    };
    s.push_str(&format!("{ind}  mass\n"));
    for p in &b.boundaries {
        row(
            &mut s,
            &format!("  {} ({}) out", p.boundary, p.kind),
            p.mass_flow,
            m,
        );
    }
    row(&mut s, "  storage", b.mass.storage, m);
    row(
        &mut s,
        &format!("  imbalance ({:+.2e})", b.mass.relative),
        b.mass.imbalance,
        m,
    );
    s.push_str(&format!("{ind}  kinetic energy\n"));
    let wall_rows = |s: &mut String, p: &BoundaryTerms, pre: &str| {
        if p.viscous_power != 0.0 {
            row(
                s,
                &format!("{pre}{} viscous work in", p.boundary),
                p.viscous_power,
                e,
            );
        }
        if p.convective_power != 0.0 {
            row(
                s,
                &format!("{pre}{} convected out", p.boundary),
                p.convective_power,
                e,
            );
        }
        if p.pressure_power != 0.0 {
            row(
                s,
                &format!("{pre}{} pressure work out", p.boundary),
                p.pressure_power,
                e,
            );
        }
    };
    // the shaft work of a rotor (D184) with the terms of its walls
    for t in &b.shafts {
        row(
            &mut s,
            &format!("  shaft work of rotor {} in (-Ω·T)", t.rotor),
            t.power,
            e,
        );
        for p in b
            .boundaries
            .iter()
            .filter(|p| t.walls.contains(&p.boundary))
        {
            wall_rows(&mut s, p, "    of it ");
        }
        if t.transposed != 0.0 {
            row(
                &mut s,
                "    of it transposed viscous stress",
                t.transposed,
                e,
            );
        }
    }
    for p in &b.boundaries {
        if !b.shafts.iter().any(|t| t.walls.contains(&p.boundary)) {
            wall_rows(&mut s, p, "  ");
        }
    }
    let k = &b.kinetic_energy;
    if k.body_force != 0.0 {
        row(&mut s, "  body force", k.body_force, e);
    }
    row(&mut s, "  viscous dissipation", k.viscous_dissipation, e);
    row(&mut s, "  storage", k.storage, e);
    row(
        &mut s,
        &format!("  numerical dissipation ({:+.2e})", k.relative),
        k.numerical_dissipation,
        e,
    );
    row(&mut s, "    of it pressure", k.pressure_dissipation, e);
    row(&mut s, "    of it schemes", k.scheme_dissipation, e);
    s
}

/// The ledger in the report (`flow.md`): one table per solve.
pub fn markdown(s: &BalanceSummary, cells: usize) -> String {
    let (b, unit_m, unit_e, what) = match &s.run {
        Some(r) => (r, "kg", "J", "over the run"),
        None => (&s.last, "kg/s", "W", "of the final state"),
    };
    let mut o = format!(
        "\n{cells} cells, conservation ledger {what}{}:\n\n| Term | Value |\n|---|---|\n",
        match (s.steps, s.negative_steps) {
            (Some(n), Some(k)) => format!(" ({n} steps, {k} with negative numerical dissipation)"),
            _ => String::new(),
        }
    );
    let mut row =
        |name: String, v: f64, u: &str| o.push_str(&format!("| {name} | {v:+.6e} {u} |\n"));
    row("mass imbalance".into(), b.mass.imbalance, unit_m);
    row(
        "largest relative mass imbalance (-)".into(),
        s.largest_mass_relative,
        "",
    );
    for t in &b.shafts {
        row(
            format!(
                "shaft work of rotor {} (−Ω·T, in the power of {})",
                t.rotor,
                t.walls.join(", ")
            ),
            t.power,
            unit_e,
        );
    }
    for p in &b.boundaries {
        if p.power() != 0.0 {
            row(
                format!("power into the fluid through {}", p.boundary),
                p.power(),
                unit_e,
            );
        }
    }
    let k = &b.kinetic_energy;
    if k.body_force != 0.0 {
        row("body force".into(), k.body_force, unit_e);
    }
    row("viscous dissipation".into(), k.viscous_dissipation, unit_e);
    row(
        "storage (change of kinetic energy)".into(),
        k.storage,
        unit_e,
    );
    row(
        "numerical dissipation".into(),
        k.numerical_dissipation,
        unit_e,
    );
    row("numerical dissipation, relative (-)".into(), k.relative, "");
    row(
        "of it pressure (momentum interpolation)".into(),
        k.pressure_dissipation,
        unit_e,
    );
    row(
        "of it convection and time schemes".into(),
        k.scheme_dissipation,
        unit_e,
    );
    o
}

const CSV_TERMS: &str = "mass_outflow,mass_relative,kinetic_energy,storage,boundary_power,\
body_force_power,viscous_dissipation,numerical_dissipation,numerical_relative,pressure_dissipation,\
scheme_dissipation";

/// The units line (of the D75 preamble) and the header of
/// `balance.csv`: per time step with the time (transient), or per outer
/// iteration (steady); with rotors (D184) a last column `shaft_power`.
pub fn csv_head(transient: bool, shaft: bool) -> (String, String) {
    let mut units = String::from(
        "mass_outflow kg/s, mass_relative -, kinetic_energy J, storage boundary_power \
         body_force_power viscous_dissipation numerical_dissipation W, numerical_relative -, \
         pressure_dissipation scheme_dissipation W",
    );
    let mut terms = String::from(CSV_TERMS);
    if shaft {
        units.push_str(", shaft_power W (of the rotors, in boundary_power)");
        terms.push_str(",shaft_power");
    }
    if transient {
        (
            format!("step -, time s, {units}"),
            format!("step,time,{terms}"),
        )
    } else {
        (
            format!("iteration -, {units}"),
            format!("iteration,{terms}"),
        )
    }
}

/// A row of `balance.csv` (with the time when the balance has one).
pub fn csv_row(b: &Balance) -> String {
    let k = &b.kinetic_energy;
    let mut s = format!("{}", b.step.or(b.iteration).unwrap_or(0));
    if let Some(t) = b.time {
        s.push_str(&format!(",{t:e}"));
    }
    for v in [
        b.mass.outflow,
        b.mass.relative,
        k.energy,
        k.storage,
        k.boundary,
        k.body_force,
        k.viscous_dissipation,
        k.numerical_dissipation,
        k.relative,
        k.pressure_dissipation,
        k.scheme_dissipation,
    ] {
        s.push_str(&format!(",{v:e}"));
    }
    if !b.shafts.is_empty() {
        s.push_str(&format!(",{:e}", k.shaft));
    }
    s
}

/// The `balance` row of the run record (D69).
pub fn record_row(b: &Balance) -> serde_json::Value {
    let k = &b.kinetic_energy;
    let mut v = serde_json::json!({
        "type": "balance",
        "mass_outflow": b.mass.outflow,
        "mass_relative": b.mass.relative,
        "kinetic_energy": k.energy,
        "storage": k.storage,
        "boundary_power": k.boundary,
        "body_force_power": k.body_force,
        "viscous_dissipation": k.viscous_dissipation,
        "numerical_dissipation": k.numerical_dissipation,
        "numerical_relative": k.relative,
        "pressure_dissipation": k.pressure_dissipation,
        "scheme_dissipation": k.scheme_dissipation,
    });
    if let Some(m) = v.as_object_mut() {
        if !b.shafts.is_empty() {
            m.insert("shaft_power".into(), k.shaft.into());
        }
        if let Some(s) = b.step {
            m.insert("step".into(), s.into());
        }
        if let Some(i) = b.iteration {
            m.insert("iteration".into(), i.into());
        }
        if let Some(t) = b.time {
            m.insert("time".into(), t.into());
        }
    }
    v
}

/// The end-of-run chart: bars of the kinetic-energy terms (W, or J over a
/// transient run) with the numerical dissipation.
pub fn chart(name: &str, b: &Balance, integrated: bool) -> String {
    let e = if integrated { "J" } else { "W" };
    let mut bars: Vec<(String, f64)> = Vec::new();
    for p in &b.boundaries {
        let w = p.power();
        if w != 0.0 {
            bars.push((format!("{} (power in)", p.boundary), w));
        }
    }
    let k = &b.kinetic_energy;
    if k.body_force != 0.0 {
        bars.push(("body force".into(), k.body_force));
    }
    bars.push(("viscous dissipation".into(), -k.viscous_dissipation));
    bars.push(("numerical: pressure".into(), -k.pressure_dissipation));
    bars.push(("numerical: schemes".into(), -k.scheme_dissipation));
    bars.push(("storage".into(), -k.storage));
    crate::plot::bars(
        &format!(
            "{name}: kinetic energy{}",
            if integrated { " over the run" } else { "" }
        ),
        &format!("{e} (sources +, sinks −)"),
        &bars,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rates(ke: f64, lid: f64, phi: f64) -> Rates {
        Rates {
            kinetic_energy: ke,
            boundaries: vec![BoundaryTerms {
                boundary: "lid".into(),
                kind: "wall".into(),
                viscous_power: lid,
                ..Default::default()
            }],
            body_force_power: 0.0,
            pressure_source_power: 0.0,
            viscous_dissipation: phi,
            mass_scale: 1.0,
            rotors: Default::default(),
        }
    }

    #[test]
    fn steady_ledger_closes() {
        let b = steady(&rates(2.0, 1.0, 0.9), 7);
        assert_eq!(b.kinetic_energy.storage, 0.0);
        assert!((b.kinetic_energy.numerical_dissipation - 0.1).abs() < 1e-15);
        assert_eq!(b.mass.imbalance, 0.0);
    }

    #[test]
    fn transient_ledger_integrates_with_the_trapezoidal_rule() {
        // K(t) = t², power in 2t + 1, dissipation 1: exact, no numerical
        // dissipation, for the linear rates
        let mut l = Ledger::start(rates(0.0, 1.0, 1.0), 0.0, 0, 1e-4);
        let mut t: f64 = 0.0;
        for n in 1..=10 {
            t = 0.1 * n as f64;
            let b = l.step(rates(t * t, 2.0 * t + 1.0, 1.0), n, t);
            assert!(b.kinetic_energy.numerical_dissipation.abs() < 1e-12);
        }
        let r = l.run();
        assert!((r.kinetic_energy.storage - t * t).abs() < 1e-12);
        assert!((r.boundaries[0].viscous_power - (t * t + t)).abs() < 1e-12);
        assert!(r.kinetic_energy.numerical_dissipation.abs() < 1e-12);
        assert_eq!(l.negative_steps, 0);
    }

    #[test]
    fn negative_dissipation_and_mass_are_flagged() {
        let s = BalanceSummary {
            last: steady(&rates(1.0, 1.0, 1.1), 3),
            largest_mass_relative: 1e-6,
            ..Default::default()
        };
        let d = check(&s, &BalanceSection::default());
        let codes: Vec<&str> = d.iter().map(|d| d.code).collect();
        assert_eq!(codes, [diag::BAL_MASS_STEP, diag::BAL_NEGATIVE_DISSIPATION]);
    }

    /// A thermal-like account of a later module: two boundaries with two
    /// terms each, a source, a part that is not closed and one side of an
    /// exchange pair, closed by the one routine.
    fn heat_account() -> Account {
        let mut a = Account::new("thermal energy", "J");
        a.terms = vec![
            Term::boundary("conduction", "hot", 3.0),
            Term::boundary("convected", "hot", -0.5),
            Term::boundary("conduction", "cold", -1.25),
            Term::volume("heater", TermKind::Source, 0.75),
            Term::volume("to the solid", TermKind::Exchange, -0.25).paired(7),
            Term::volume("of it radiation", TermKind::Part, 100.0),
        ];
        a
    }

    #[test]
    fn one_closure_routine_closes_any_account() {
        let a = heat_account();
        assert_eq!(a.boundaries(), ["hot", "cold"]);
        assert_eq!(a.boundary_sums(), [2.5, -1.25]);
        let c = a.close(1.5);
        assert_eq!(c.boundary, 1.25);
        // 1.25 + 0.75 − 0.25; the part is reported, not closed
        assert_eq!(c.contributions, 1.75);
        assert_eq!(c.imbalance, -0.25);
        assert_eq!(c.largest, 100.0);
        assert_eq!(c.relative, -0.25 / 100.0);
        let mut b = a.clone();
        b.scale = Some(2.0);
        assert_eq!(b.close(1.5).relative, -0.125);
        assert_eq!(a.value("conduction", Some("cold")), -1.25);
        assert_eq!(a.value("missing", None), 0.0);
    }

    #[test]
    fn exchange_pairs_cancel_to_the_bit() {
        let gas = heat_account();
        let mut solid = Account::new("thermal energy", "J");
        solid.region = "solid".into();
        solid.terms = vec![Term::volume("from the gas", TermKind::Exchange, 0.25).paired(7)];
        assert_eq!(exchange_residuals(&[gas.clone(), solid]), [(7, 0.0)]);
        assert_eq!(exchange_residuals(&[gas]), [(7, -0.25)]);
    }

    #[test]
    fn step_rules_and_run_integrals() {
        // a trapezoidal term and an end-of-step term over two steps
        let at = |v: f64| {
            let mut a = Account::new("q", "kg");
            a.scale = Some(v);
            a.terms = vec![
                Term::volume("source", TermKind::Source, v),
                Term::boundary("outflow", "out", -v).rule(StepRule::End),
            ];
            a
        };
        let s1 = Account::over_step(&at(1.0), &at(3.0));
        assert_eq!(s1.terms[0].value, 2.0);
        assert_eq!(s1.terms[1].value, -3.0);
        let mut run = at(1.0).scaled(0.0);
        run.accumulate(&s1, 0.5);
        run.accumulate(&Account::over_step(&at(3.0), &at(5.0)), 0.25);
        assert_eq!(run.terms[0].value, 0.5 * 2.0 + 0.25 * 4.0);
        assert_eq!(run.terms[1].value, -(0.5 * 3.0 + 0.25 * 5.0));
        assert_eq!(run.scale, Some(0.5 * 3.0 + 0.25 * 5.0));
    }

    #[test]
    fn mass_and_kinetic_energy_are_accounts() {
        let r = rates(2.0, 1.0, 0.9);
        let acc = r.accounts();
        assert_eq!(acc.len(), 2);
        assert_eq!(acc[0].quantity, "mass");
        assert_eq!(acc[1].quantity, "kinetic energy");
        assert_eq!(acc[1].content, Some(2.0));
        let b = steady(&r, 1);
        // the fixed views and the closure of the account agree
        assert_eq!(
            b.kinetic_energy.numerical_dissipation,
            -b.accounts[1].close(0.0).imbalance
        );
        assert_eq!(b.boundaries[0].kind, "wall");
        // a closed domain reports +0 outflow, as before the accounts
        assert!(b.mass.outflow.is_sign_positive());
    }
}
