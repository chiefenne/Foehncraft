//! The discrete steady residual R(q; x, p) as one callable (L17 stage A0,
//! decision Q4 of the adjoint plan, docs/requirements/step27-adjoint-solver.md).
//!
//! The discrete adjoint (D220) differentiates the equations a converged
//! steady run satisfies, not the iteration that got there. SIMPLEC and the
//! coupled solver converge to the same discrete equations (see
//! [`crate::coupled`]), so one function serves both:
//!
//! ```text
//! R_u,i = b_i − A_i u_i                          (momentum, per component)
//! R_c   = Σ_f ± ṁ_f                              (continuity, per cell)
//! R_ṁ   = (ṁ_f − ṁ*_f(u, p, ṁ)) / α_u            (face flux, per face)
//! ```
//!
//! with A, b the unrelaxed momentum rows assembled from the current fields
//! (convection with the current face fluxes ṁ, the pressure gradient on
//! the right-hand side) and ṁ* the Rhie–Chow flux of the current fields
//! with the relaxation term (1 − α_u)(ṁ − ρ ū·S) of the current fields
//! (Majumdar 1988).
//!
//! **State vector (A0 decision).** The face fluxes are part of the state:
//! q = (u, p, ṁ). The Rhie–Chow flux reads the momentum diagonal D = α V/a_P,
//! and a_P holds the convection coefficients of ṁ, so ṁ is defined only
//! implicitly by R_ṁ = 0 and is not a function of (u, p) alone. Written
//! out, ṁ − ṁ* = α (ṁ − ρ ū·S) + ρ D̄ |S|²/(d·S) [jump], and D̄/α = V/a_P
//! does not hold α: dividing by α makes R_ṁ independent of the relaxation,
//! so the converged state and the residual are the same for every α_u and
//! for both algorithms.
//!
//! Scope of A0: steady, laminar, constant density, one rank. Turbulence,
//! energy and passive scalars (their equations are not in R yet), the
//! variable density and time steps are refused with a diagnostic.
//! Boundary values, gradients and limiters are brought up to date for the
//! evaluation and put back afterwards (`Flow::lagged`), so a call between
//! two outer iterations does not change the run.

use crate::diag::{self, Diagnostic};
use crate::simplec::Flow;
use anyhow::Result;
use foehncraft_kernels::flow as k;

/// The residual of the steady discrete equations at the current state.
#[derive(Debug, Clone)]
pub struct SteadyResidual {
    /// R_u per component and cell (N, the units of the momentum rows).
    pub momentum: [Vec<f64>; 3],
    /// R_c per cell (kg/s): the net outflow of the current face fluxes.
    pub continuity: Vec<f64>,
    /// R_ṁ per face (kg/s): (ṁ − ṁ*)/α_u.
    pub mass_flux: Vec<f64>,
    /// The residual scale of the momentum rows, ‖a_P u‖ + ‖b‖ over the
    /// three components (as the run's own residuals).
    pub momentum_scale: f64,
    /// The scale of the fluxes, ‖ṁ‖ over all faces.
    pub flux_scale: f64,
}

impl SteadyResidual {
    /// The three norms ‖R_u‖/scale, ‖R_c‖/‖ṁ‖, ‖R_ṁ‖/‖ṁ‖ (2-norms).
    pub fn norms(&self) -> [f64; 3] {
        let n = |v: &[f64]| v.iter().map(|a| a * a).sum::<f64>().sqrt();
        let mu = self.momentum.iter().map(|v| n(v).powi(2)).sum::<f64>().sqrt();
        let fs = self.flux_scale.max(f64::MIN_POSITIVE);
        [
            mu / self.momentum_scale.max(f64::MIN_POSITIVE),
            n(&self.continuity) / fs,
            n(&self.mass_flux) / fs,
        ]
    }
}

fn unsupported(what: &str) -> anyhow::Error {
    diag::err(Diagnostic::new(
        diag::CASE_INCONSISTENT,
        format!(
            "the steady residual of the adjoint (L17 stage A0) covers steady laminar flow of \
             constant density on one rank; this run has {what}"
        ),
    ))
}

impl Flow {
    /// The residual R(q) of the steady discrete equations at the current
    /// fields (see the module documentation). The run continues unchanged
    /// after the call.
    pub fn steady_residual(&mut self) -> Result<SteadyResidual> {
        if self.transient {
            return Err(unsupported("time steps"));
        }
        if self.turbulence.is_some() {
            return Err(unsupported("a turbulence model"));
        }
        if self.scalars.is_some() {
            return Err(unsupported("transported scalars (energy or passive scalars)"));
        }
        if self.density.is_some() {
            return Err(unsupported("a variable density"));
        }
        if self.comm.is_some() {
            return Err(unsupported("several ranks"));
        }
        let lagged = self.lagged()?;
        let r = self.evaluate_steady_residual();
        self.set_lagged(&lagged)?;
        r
    }

    fn evaluate_steady_residual(&mut self) -> Result<SteadyResidual> {
        let alpha = self.settings.alpha_u(false);
        self.update_boundaries()?;
        self.gradients()?;
        // the relaxation term of the current fields into w.relax, then the
        // unrelaxed momentum rows with the current fluxes
        self.relaxation_term(alpha)?;
        let conv = self.convecting_flux()?;
        let scale = self.momentum_assembly(conv, alpha, 1.0, None)?;
        let (f, w) = (self.f, self.w);
        let sym = self.bc.mom_kind.contains(&k::MOM_SYMMETRY);
        let tmp = w.gpc.x;
        let mut momentum: [Vec<f64>; 3] = Default::default();
        for (i, out) in momentum.iter_mut().enumerate() {
            if sym {
                self.dev.set_diagonal(w.val, w.diag.all()[i])?;
            }
            self.dev.matrix_product(w.val, f.u.all()[i], tmp)?;
            self.dev.lincomb(tmp, 1.0, w.rhs.all()[i], -1.0, tmp)?;
            *out = self.dev.download(tmp)?;
        }
        // the Rhie–Chow flux of the current fields and R_ṁ
        let a = self.flux_args(Some(w.d), w.relax, None);
        self.dev.mass_flux(&a, w.m)?;
        self.close_backflow(w.m)?;
        let m = self.dev.download(f.mass_flux)?;
        let ms = self.dev.download(w.m)?;
        let mass_flux = m.iter().zip(&ms).map(|(a, b)| (a - b) / alpha).collect();
        self.dev.net_outflow(f.mass_flux, w.defect)?;
        let continuity = self.dev.download(w.defect)?;
        let flux_scale = m.iter().map(|a| a * a).sum::<f64>().sqrt();
        Ok(SteadyResidual {
            momentum,
            continuity,
            mass_flux,
            momentum_scale: scale,
            flux_scale,
        })
    }
}

/// A steady laminar flow solver set up from `case` at its start fields,
/// for the residual studies of stage A0 (the example `adjoint_a0`): the
/// mesh, the boundary arrays and the physics of the case, the reference
/// cell 0 when no boundary fixes the pressure. The full run setup (wall
/// distance, sources, scalars, remediation report) is
/// [`crate::case::run`]'s; this one refuses what the residual does not
/// cover yet.
pub fn steady_flow(case: &crate::case::FlowCase, threads: Option<usize>) -> Result<Flow> {
    if case.time.is_some() {
        return Err(unsupported("a [time] table"));
    }
    let (mesh, geo, _) = crate::case::mesh_for(case, None)?;
    let phys = crate::case::physics(case)?;
    let arrays = crate::bc::arrays(
        &mesh,
        &geo,
        &case.boundary,
        &phys,
        crate::bc::At::numbers(None),
    )?;
    let disc = crate::simplec::discretisation(&mesh, &geo, &arrays)?;
    let reference = (!arrays.has_pressure()).then_some(0);
    Flow::new(
        disc,
        arrays,
        phys,
        case.solver.clone(),
        case.initial.velocity,
        case.initial.pressure,
        reference,
        &crate::simplec::Target {
            device: "cpu".into(),
            threads,
            ptx: Default::default(),
            ranks: None,
        },
    )
}
