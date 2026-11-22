//! The coupled pressure–velocity algorithm (L11, D19, D176) and the
//! pseudo-time steps of a steady run (D177, D178): the momentum equations
//! and the continuity equation of one outer iteration are solved as one
//! linear system for the increments of the velocity and the pressure, on
//! the discretisation of [`crate::simplec`].
//!
//! One outer iteration, from the fields u, p and the face fluxes ṁ:
//!
//! 1. boundary values, gradients, limiters and the momentum rows A as in
//!    SIMPLEC, with the pseudo-time term ρ V/Δτ (u − uᵏ) on the diagonal
//!    (uᵏ: the current iterate, so the term vanishes at convergence) and the
//!    relaxation α_u;
//! 2. the residuals of the current fields: r_u = b − A u − V ∇p per
//!    component and r_c = −Σ_f ±ṁ*_f, with ṁ* the Rhie–Chow fluxes of u
//!    and p ([`foehncraft_kernels::flow::face_mass_flux`]);
//! 3. the increments δu, δp from
//!
//!    ```text
//!    | A  G | |δu|   |r_u|        G δp = V ∇(δp),
//!    | D  L | |δp| = |r_c|        D δu = Σ_f ± ρ (δu)_f·S,
//!    ```
//!
//!    L δp = Σ_f c_f (δp_P − δp_N) with c_f = ρ D̄_f |S|²/(d·S),
//!    D = V/(a_P/α_u + ρ V/Δτ): the linearisation of the same flux formula
//!    (the interpolated pressure gradient of the Rhie–Chow term and the
//!    skewness corrections stay at the old fields; they are in r_c);
//! 4. u += δu, p += δp, ṁ = ṁ* + ρ (δu)_f·S − c_f (δp_N − δp_P), then one
//!    pressure correction on the remaining continuity defect: the face
//!    fluxes satisfy continuity to the tolerance of `[solver.pressure]`
//!    after every outer iteration, as with SIMPLEC;
//! 5. turbulence and scalars as in SIMPLEC.
//!
//! The residuals are those of the segregated algorithm and the face flux
//! depends neither on α_u (Majumdar 1988) nor on Δτ (the time-step term of
//! D27 with the current iterate as the old level), so a converged result
//! is the one of SIMPLEC. The discrete equations are the same too: what the
//! convection scheme conserves (L39) it conserves with either algorithm.
//!
//! **Pseudo-time steps.** Under-relaxation is a time step of its own,
//! Δτ_P = α/(1 − α) ρV/a_P per cell: it shrinks with the cell, and the
//! number of outer iterations grows with the mesh. The momentum equations
//! therefore take one step Δτ for all cells, by default `factor` × L / U
//! with L the volume of the domain over the area of its boundaries
//! (without symmetry planes: the half height of a plane channel, a quarter
//! of the diameter of a pipe, and no function of the thickness of a mesh
//! of one cell layer) and U the larger of the largest prescribed boundary
//! speed (inlets, moving walls) and the root mean square of the cell
//! speeds; α_u (close to 1) only caps the local Courant number at
//! α/(1 − α). The turbulence equations take
//! min(Δτ, c_T T) per cell with the time scale T of the model
//! ([`crate::turbulence`]): their sources are stiff on that scale.
//!
//! **Linear solve.** The operator is applied matrix-free with the
//! operations of the device (the momentum matrix, the least-squares or
//! Gauss gradient with the mask of a pressure correction, the flux of a
//! velocity increment, the pressure matrix), so the system needs no block
//! matrix. Flexible GMRES (Saad 1993) with right preconditioning; the
//! preconditioner is one SIMPLE-type sweep on the increments (Patankar &
//! Spalding 1972; as a preconditioner: Klaij & Vuik 2013):
//!
//! ```text
//! z_u = Ã⁻¹ r_u,   Ŝ z_p = r_p − D z_u,   z_u −= D̂ ∇z_p,
//! ```
//!
//! with Ã⁻¹ BiCGStab per component (the outer iteration needs these sweeps
//! to about 10⁻²: looser sweeps make the Krylov count and the outer
//! iteration erratic), D̂ = min(D^C, V/((1 − α) a_P/α + ρV/(α Δτ))) (the
//! SIMPLEC diagonal of Van Doormaal & Raithby 1984, bounded by its value
//! for fluxes that satisfy continuity), and Ŝ the pressure matrix with the
//! coefficient ρ (D̄ + D̂̄)_f |S|²/(d·S) (the Rhie–Chow part L plus the
//! compact form of −D D̂ G), solved by CG with the multigrid of the
//! pressure.

use crate::device::{Buf, Buf3, Lsq, Mask, TimeTerms};
use crate::simplec::{Algorithm, Flow, Residuals};
use crate::turbulence::PseudoStep;
use anyhow::Result;
use foehncraft_fv::amg::ReductionSettings;
use foehncraft_fv::solver::Preconditioner;
use foehncraft_kernels::flow as k;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// `[solver.coupled]`: the linear solve of the coupled algorithm.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct CoupledSettings {
    /// Relative reduction of the residual of the coupled system per outer
    /// iteration.
    pub tolerance: f64,
    /// Largest number of Krylov iterations per outer iteration.
    pub max_iter: usize,
    /// The momentum sweeps of the preconditioner: cycles of a reduction
    /// multigrid (`multigrid`, D179) or BiCGStab with the Jacobi
    /// preconditioner (`bicgstab`).
    pub momentum_sweep: MomentumSweep,
    /// Multigrid cycles per momentum sweep.
    pub momentum_cycles: usize,
    /// A multigrid hierarchy of the momentum sweeps is kept while the
    /// momentum matrix differs from the one it was built from by at most
    /// this fraction (the sum of the absolute changes of the coefficients
    /// over the sum of their absolute values); 0: a setup in every outer
    /// iteration.
    pub momentum_reuse: f64,
    /// Relative tolerance of a BiCGStab momentum sweep.
    pub momentum_tolerance: f64,
    /// Largest number of BiCGStab iterations of a momentum sweep.
    pub momentum_max_iter: usize,
    /// Relative tolerance of the pressure solve of the preconditioner. The
    /// Krylov iteration corrects what this solve leaves: 0.3 needs the
    /// same Krylov iterations as 0.1 with a fifth to a quarter fewer CG
    /// iterations (requirements step 22, section 7.9).
    pub pressure_tolerance: f64,
    /// Largest number of CG iterations of that pressure solve.
    pub pressure_max_iter: usize,
    /// The multigrid hierarchy of the pressure matrix Ŝ is built anew when
    /// Ŝ differs from the matrix of its last setup by more than this
    /// fraction (the measure of `momentum_reuse`), besides every
    /// `amg_rebuild` iterations; 0: whenever Ŝ has changed; 1: only every
    /// `amg_rebuild` iterations (requirements step 22, section 7.9).
    pub pressure_reuse: f64,
    /// The multigrid of the momentum sweeps.
    pub multigrid: ReductionSettings,
}

impl Default for CoupledSettings {
    fn default() -> Self {
        Self {
            tolerance: 0.1,
            max_iter: 30,
            momentum_sweep: MomentumSweep::Multigrid,
            momentum_cycles: 1,
            momentum_reuse: 0.02,
            momentum_tolerance: 0.01,
            momentum_max_iter: 10,
            pressure_tolerance: 0.3,
            pressure_max_iter: 20,
            pressure_reuse: 0.05,
            multigrid: ReductionSettings::default(),
        }
    }
}

/// The momentum sweep of the preconditioner of the coupled algorithm.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum MomentumSweep {
    /// Cycles of the reduction multigrid of the momentum matrix (D179);
    /// a hierarchy is kept while the matrix changes little
    /// (`momentum_reuse`).
    #[default]
    Multigrid,
    /// BiCGStab with the Jacobi preconditioner, to `momentum_tolerance` or
    /// `momentum_max_iter`.
    Bicgstab,
}

/// First multigrid slot of the momentum sweeps (one per velocity
/// component; 0: the pressure, 1 to 4: the turbulence equations).
const MOMENTUM_SLOT: usize = 5;

/// A velocity component whose part of a vector holds less than this
/// fraction of the squared norm of the three components takes no sweep
/// (the third component of a planar case, zero up to rounding).
const NEGLIGIBLE: f64 = 1e-20;

/// A kept hierarchy whose cycle leaves more than this fraction of the
/// momentum residual of the first Krylov vector of an outer iteration is
/// rebuilt at once (a fresh one leaves 0.1 to 0.4 of it).
const STALE: f64 = 0.5;

/// The multigrid hierarchies of the momentum sweeps (D179) that outlive an
/// outer iteration: the hierarchy (0 to 2, in the multigrid slots 5 to 7)
/// each velocity component uses, and the momentum matrix each hierarchy
/// was built from. A restart file holds them, so that a continued run
/// builds the same hierarchies.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SweepState {
    /// Hierarchy of each velocity component (`None`: none yet, or the
    /// matrix has changed too much since its setup).
    pub slot: [Option<usize>; 3],
    /// The matrix (in the storage of the discretisation) each hierarchy
    /// was built from.
    pub matrix: [Option<Vec<f64>>; 3],
}

/// Σ|a − b| / Σ|a| of two matrices in the same storage; `diag` replaces
/// the diagonal of `a` (at the positions `dpos`).
pub(crate) fn matrix_change(a: &[f64], b: &[f64], diag: Option<&[f64]>, dpos: &[usize]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for (x, y) in a.iter().zip(b) {
        num += (x - y).abs();
        den += x.abs();
    }
    if let Some(d) = diag {
        for (&k, &di) in dpos.iter().zip(d) {
            num += (di - b[k]).abs() - (a[k] - b[k]).abs();
            den += di.abs() - a[k].abs();
        }
    }
    if den > 0.0 { num / den } else { f64::INFINITY }
}

impl CoupledSettings {
    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        for (name, t) in [
            ("tolerance", self.tolerance),
            ("momentum_tolerance", self.momentum_tolerance),
            ("pressure_tolerance", self.pressure_tolerance),
        ] {
            crate::need!(
                t > 0.0 && t < 1.0,
                CASE_RANGE,
                "[solver.coupled] {name} must lie in (0, 1)"
            );
        }
        for (name, n) in [
            ("max_iter", self.max_iter),
            ("momentum_max_iter", self.momentum_max_iter),
            ("pressure_max_iter", self.pressure_max_iter),
        ] {
            crate::need!(n >= 1, CASE_RANGE, "[solver.coupled] {name} must be >= 1");
        }
        crate::need!(
            self.max_iter <= 200,
            CASE_RANGE,
            "[solver.coupled] max_iter must be <= 200 (the Krylov vectors are stored)"
        );
        crate::need!(
            (1..=10).contains(&self.momentum_cycles),
            CASE_RANGE,
            "[solver.coupled] momentum_cycles must lie in 1..=10"
        );
        crate::need!(
            (0.0..=0.2).contains(&self.momentum_reuse),
            CASE_RANGE,
            "[solver.coupled] momentum_reuse must lie in [0, 0.2]"
        );
        crate::need!(
            (0.0..=1.0).contains(&self.pressure_reuse),
            CASE_RANGE,
            "[solver.coupled] pressure_reuse must lie in [0, 1]"
        );
        if let Err(e) = self.multigrid.validate() {
            crate::need!(false, CASE_RANGE, "[solver.coupled.multigrid]: {e}");
        }
        Ok(())
    }
}

/// `[solver.pseudo_time]`: the pseudo-time steps of a steady run (D177,
/// D178). A converged result does not depend on them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct PseudoTimeSettings {
    /// The step Δτ in s, the same in every cell. Default: automatic,
    /// `factor` × L / U in every outer iteration, with L the volume of the
    /// domain over the area of its boundaries (without symmetry planes and
    /// periodic boundaries) and U the larger of the largest prescribed
    /// boundary speed (inlets, moving walls) and the root mean square of
    /// the cell speeds.
    pub step: Option<f64>,
    /// Factor on the automatic step.
    pub factor: f64,
    /// The momentum equations of the coupled algorithm take the step;
    /// false: under-relaxation only. (SIMPLEC and SIMPLE relax only.)
    pub momentum: bool,
    /// c_T: the turbulence equations take the step min(Δτ, c_T T) per cell,
    /// T the time scale of the model (k/ε; 1/(β* ω); the bounded scale of
    /// k-ε-ζ-f). 0: under-relaxation only. Default 0.5 with the coupled
    /// algorithm, 0 with SIMPLEC and SIMPLE.
    pub turbulence: Option<f64>,
}

impl Default for PseudoTimeSettings {
    fn default() -> Self {
        Self {
            step: None,
            factor: 1.0,
            momentum: true,
            turbulence: None,
        }
    }
}

impl PseudoTimeSettings {
    /// c_T in force with the algorithm `alg` (0: none).
    pub fn turbulence_factor(&self, alg: Algorithm) -> f64 {
        self.turbulence
            .unwrap_or(if alg == Algorithm::Coupled { 0.5 } else { 0.0 })
    }

    /// Checks the ranges.
    pub fn validate(&self) -> Result<()> {
        if let Some(dt) = self.step {
            crate::need!(
                dt > 0.0 && dt.is_finite(),
                CASE_RANGE,
                "[solver.pseudo_time] step must be > 0"
            );
        }
        crate::need!(
            self.factor > 0.0 && self.factor.is_finite(),
            CASE_RANGE,
            "[solver.pseudo_time] factor must be > 0"
        );
        if let Some(c) = self.turbulence {
            crate::need!(
                c >= 0.0 && c.is_finite(),
                CASE_RANGE,
                "[solver.pseudo_time] turbulence must be >= 0 (0: no pseudo-time step)"
            );
        }
        Ok(())
    }
}

/// A vector of the coupled system: three velocity components and the
/// pressure, one value per cell each.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Vec4 {
    pub u: Buf3,
    pub p: Buf,
}

impl Vec4 {
    fn all(self) -> [Buf; 4] {
        [self.u.x, self.u.y, self.u.z, self.p]
    }
}

/// Buffers of the coupled solve.
#[derive(Debug)]
pub(crate) struct CoupledWork {
    /// −V per cell.
    neg_vol: Buf,
    /// 1 per cell.
    one: Buf,
    /// ρ Σ|S| over the interior faces plus half of that of the boundary
    /// faces (the scale of the continuity row; ρ per face with several
    /// fluids).
    area: f64,
    /// Several fluids (D194): ρ and 1/ρ per cell and a cell work buffer,
    /// for the pseudo-time term of D̂ per cell.
    zones: Option<(Buf, Buf, Buf)>,
    /// The Rhie–Chow pressure matrix L and its face coefficients.
    pval: Buf,
    cf: Buf,
    /// D + D̂.
    dsum: Buf,
    /// D̂.
    dhat: Buf,
    /// Right-hand side, solution, work vector.
    r: Vec4,
    x: Vec4,
    w: Vec4,
    /// Gradient, cell and face work.
    g: Buf3,
    t: Buf,
    mf: Buf,
    /// Y = ṁ − ρ ū·S of the current fields per face.
    y: Buf,
    /// Krylov basis and preconditioned basis, grown as the solves need
    /// them (a solve that ends after k iterations uses k + 1 and k).
    v: Vec<Vec4>,
    z: Vec<Vec4>,
    /// Storage position of the diagonal of every row of the momentum
    /// matrix.
    dpos: Vec<usize>,
    /// The hierarchies of the momentum sweeps built in this outer
    /// iteration (the others were kept from earlier ones).
    mg_fresh: [bool; 3],
    /// The first sweeps of this outer iteration are still to come (they
    /// check the kept hierarchies).
    mg_first: bool,
    /// A multigrid setup failed in this outer iteration: BiCGStab does its
    /// sweeps.
    mg_failed: bool,
}

/// Counters of one coupled solve.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Counts {
    pub krylov: usize,
    pub momentum: [usize; 3],
    pub pressure: usize,
    /// Reached relative residual.
    pub reduction: f64,
}

impl Flow {
    fn vec4(&mut self) -> Result<Vec4> {
        let nc = self.disc.n_cells;
        Ok(Vec4 {
            u: Buf3 {
                x: self.dev.alloc(nc)?,
                y: self.dev.alloc(nc)?,
                z: self.dev.alloc(nc)?,
            },
            p: self.dev.alloc(nc)?,
        })
    }

    /// The buffers of the coupled solve, allocated at the first coupled
    /// iteration.
    fn coupled_work(&mut self) -> Result<()> {
        if self.coupled.is_some() {
            return Ok(());
        }
        let (nc, nf) = (self.disc.n_cells, self.disc.n_faces);
        let stored = self.disc.sell.stored();
        let neg_vol = self.dev.alloc(nc)?;
        let v: Vec<f64> = self.disc.cell_volume.iter().map(|v| -v).collect();
        self.dev.upload(neg_vol, &v)?;
        let one = self.dev.alloc(nc)?;
        self.dev.fill(one, 1.0)?;
        let d = &self.disc;
        let mut area = 0.0;
        for f in 0..nf {
            let a =
                (d.face_area.x[f].powi(2) + d.face_area.y[f].powi(2) + d.face_area.z[f].powi(2))
                    .sqrt();
            let a = if f < d.n_interior { a } else { 0.5 * a };
            area += if self.phys.zones {
                self.phys.rho_face_at(f) * a
            } else {
                a
            };
        }
        let area = if self.phys.zones {
            area
        } else {
            self.phys.rho * area
        };
        let zones = match (&self.phys.rho_var, self.phys.zones) {
            (Some(r), true) => {
                let (rho, inv, work) = (
                    self.dev.alloc(nc)?,
                    self.dev.alloc(nc)?,
                    self.dev.alloc(nc)?,
                );
                let i: Vec<f64> = r.0.iter().map(|v| 1.0 / v).collect();
                self.dev.upload(rho, &r.0)?;
                self.dev.upload(inv, &i)?;
                Some((rho, inv, work))
            }
            _ => None,
        };
        let cw = CoupledWork {
            neg_vol,
            one,
            area,
            zones,
            pval: self.dev.alloc(stored)?,
            cf: self.dev.alloc(nf)?,
            dsum: self.dev.alloc(nc)?,
            dhat: self.dev.alloc(nc)?,
            r: self.vec4()?,
            x: self.vec4()?,
            w: self.vec4()?,
            g: Buf3 {
                x: self.dev.alloc(nc)?,
                y: self.dev.alloc(nc)?,
                z: self.dev.alloc(nc)?,
            },
            t: self.dev.alloc(nc)?,
            mf: self.dev.alloc(nf)?,
            y: self.dev.alloc(nf)?,
            v: Vec::new(),
            z: Vec::new(),
            dpos: (0..nc).map(|r| self.disc.sell.pos(r, 0)).collect(),
            mg_fresh: [false; 3],
            mg_first: true,
            mg_failed: false,
        };
        self.coupled = Some(cw);
        Ok(())
    }

    /// The Krylov vectors v_j and z_j (allocated at their first use).
    fn krylov_vectors(&mut self, j: usize) -> Result<(Vec4, Vec4)> {
        let (nv, nz) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (cw.v.len(), cw.z.len())
        };
        for _ in nv..=j {
            let v = self.vec4()?;
            self.coupled.as_mut().expect("coupled work").v.push(v);
        }
        for _ in nz..=j {
            let z = self.vec4()?;
            self.coupled.as_mut().expect("coupled work").z.push(z);
        }
        let cw = self.coupled.as_ref().expect("coupled work");
        Ok((cw.v[j], cw.z[j]))
    }

    fn dot4(&mut self, a: Vec4, b: Vec4) -> Result<f64> {
        let mut s = 0.0;
        for (x, y) in a.all().into_iter().zip(b.all()) {
            s += self.dev.dot(x, y)?;
        }
        Ok(s)
    }

    fn norm4(&mut self, a: Vec4) -> Result<f64> {
        let mut s = 0.0;
        for x in a.all() {
            s += self.dev.sum_sq(x)?;
        }
        Ok(s.sqrt())
    }

    /// out = a x + b y (out may be x or y).
    fn lincomb4(&mut self, out: Vec4, a: f64, x: Vec4, b: f64, y: Vec4) -> Result<()> {
        for ((o, x), y) in out.all().into_iter().zip(x.all()).zip(y.all()) {
            self.dev.lincomb(o, a, x, b, y)?;
        }
        Ok(())
    }

    /// The momentum matrix of component `i` in `w.val` (the diagonal of the
    /// symmetry planes differs per component).
    fn momentum_matrix(&mut self, i: usize, sym: bool) -> Result<()> {
        if sym {
            self.dev.set_diagonal(self.w.val, self.w.diag.all()[i])?;
        }
        Ok(())
    }

    /// Before the sweeps of an outer iteration: the hierarchies of the
    /// momentum sweeps whose matrix has changed by more than
    /// `momentum_reuse` since their setup are given up (they are rebuilt
    /// when a sweep needs them).
    fn momentum_hierarchies_check(&mut self, sym: bool) -> Result<()> {
        let reuse = self.settings.coupled.momentum_reuse;
        {
            let cw = self.coupled.as_mut().expect("coupled work");
            cw.mg_fresh = [false; 3];
            cw.mg_first = true;
            cw.mg_failed = false;
        }
        if self.settings.coupled.momentum_sweep != MomentumSweep::Multigrid {
            return Ok(());
        }
        if self.sweep.slot.iter().any(Option::is_some) {
            let val = self.dev.download(self.w.val)?;
            for i in 0..3 {
                let Some(slot) = self.sweep.slot[i] else {
                    continue;
                };
                let diag = if sym {
                    Some(self.dev.download(self.w.diag.all()[i])?)
                } else {
                    None
                };
                let cw = self.coupled.as_ref().expect("coupled work");
                let kept = match &self.sweep.matrix[slot] {
                    Some(m) if m.len() == val.len() => {
                        matrix_change(&val, m, diag.as_deref(), &cw.dpos) <= reuse
                    }
                    _ => false,
                };
                if !kept {
                    self.sweep.slot[i] = None;
                }
            }
        }
        self.sweep_release();
        Ok(())
    }

    /// Forgets the matrices of the hierarchies no component uses.
    fn sweep_release(&mut self) {
        for slot in 0..3 {
            if !self.sweep.slot.contains(&Some(slot)) {
                self.sweep.matrix[slot] = None;
            }
        }
    }

    /// Makes the multigrid hierarchy of the momentum matrix of component
    /// `i` (in `w.val`) the active one. A hierarchy kept from an earlier
    /// outer iteration serves while the matrix has changed little
    /// ([`Self::momentum_hierarchies_check`]); otherwise it is built at the
    /// first sweep that needs it, so that a component that never takes a
    /// sweep costs no setup. Components with the same matrix share a
    /// hierarchy (all three without symmetry planes; with them the
    /// diagonals differ). Returns `false` when a setup failed: BiCGStab
    /// then does the sweeps of this outer iteration.
    fn momentum_hierarchy(&mut self, i: usize) -> Result<bool> {
        if self.coupled.as_ref().expect("coupled work").mg_failed {
            return Ok(false);
        }
        if let Some(slot) = self.sweep.slot[i] {
            self.dev.amg_select(MOMENTUM_SLOT + slot)?;
            return Ok(true);
        }
        let t = Instant::now();
        let cur = self.dev.download(self.w.val)?;
        // a hierarchy built in this outer iteration from the same matrix
        let fresh = self.coupled.as_ref().expect("coupled work").mg_fresh;
        let shared = (0..3).filter(|&j| j != i && fresh[j]).find_map(|j| {
            let slot = self.sweep.slot[j]?;
            (self.sweep.matrix[slot].as_ref() == Some(&cur)).then_some(slot)
        });
        let slot = match shared {
            Some(slot) => {
                self.dev.amg_select(MOMENTUM_SLOT + slot)?;
                slot
            }
            None => {
                let slot = (0..3)
                    .find(|k| !self.sweep.slot.contains(&Some(*k)))
                    .expect("three components, three hierarchies");
                let ms = self.settings.coupled.multigrid.clone();
                let built = self
                    .dev
                    .reduction_setup_slot(MOMENTUM_SLOT + slot, self.w.val, &ms);
                if let Err(e) = built {
                    self.timings.amg_setup += t.elapsed().as_secs_f64();
                    self.coupled.as_mut().expect("coupled work").mg_failed = true;
                    self.sweep_fallbacks += 1;
                    if self.sweep_fallbacks == 1 {
                        self.sweep_note = Some(format!(
                            "outer iteration {}: the multigrid setup of the momentum sweep \
                             failed ({e:#}); BiCGStab does the sweeps of such iterations",
                            self.iteration
                        ));
                    }
                    self.dev.amg_select(0)?;
                    return Ok(false);
                }
                self.sweep_setups += 1;
                self.sweep.matrix[slot] = Some(cur);
                slot
            }
        };
        self.sweep.slot[i] = Some(slot);
        self.coupled.as_mut().expect("coupled work").mg_fresh[i] = true;
        self.timings.amg_setup += t.elapsed().as_secs_f64();
        Ok(true)
    }

    /// Rebuilds the hierarchies of the momentum sweeps of a restart from
    /// their matrices, so that the run continues bit-identically.
    pub fn restore_sweep(&mut self, sweep: &SweepState) -> Result<()> {
        let stored = self.disc.sell.stored();
        let ms = self.settings.coupled.multigrid.clone();
        for (slot, m) in sweep.matrix.iter().enumerate() {
            let Some(m) = m else {
                continue;
            };
            anyhow::ensure!(
                m.len() == stored,
                "restart: the matrix of the momentum multigrid has {} values, the mesh needs \
                 {stored}",
                m.len()
            );
            self.dev.upload(self.w.val, m)?;
            self.dev
                .reduction_setup_slot(MOMENTUM_SLOT + slot, self.w.val, &ms)?;
        }
        self.dev.amg_select(0)?;
        self.sweep = sweep.clone();
        self.sweep_release();
        Ok(())
    }

    /// out = K x: the coupled operator, the continuity row scaled with
    /// `sc`.
    fn coupled_operator(&mut self, x: Vec4, out: Vec4, sc: f64, sym: bool) -> Result<()> {
        let (g, t, mf, neg_vol, pval) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (cw.g, cw.t, cw.mf, cw.neg_vol, cw.pval)
        };
        let w = self.w;
        self.dev.gradient(
            Lsq::Pressure,
            x.p,
            Mask::PressureCorrection,
            w.zero_boundary,
            g,
        )?;
        for i in 0..3 {
            self.momentum_matrix(i, sym)?;
            self.dev
                .matrix_product(w.val, x.u.all()[i], out.u.all()[i])?;
        }
        // + V ∇x_p
        self.dev.correct_velocity(neg_vol, g, out.u)?;
        self.dev.delta_flux(x.u, mf)?;
        self.dev.net_outflow(mf, t)?;
        self.dev.matrix_product(pval, x.p, out.p)?;
        self.dev.lincomb(out.p, sc, out.p, sc, t)
    }

    /// z = M⁻¹ r: one SIMPLE-type sweep on the increments.
    fn coupled_preconditioner(
        &mut self,
        r: Vec4,
        z: Vec4,
        sc: f64,
        sym: bool,
        n: &mut Counts,
    ) -> Result<()> {
        let (g, t, mf, dhat) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (cw.g, cw.t, cw.mf, cw.dhat)
        };
        let w = self.w;
        let s = self.settings.coupled.clone();
        let tm = Instant::now();
        let setup = self.timings.amg_setup;
        let mut part = [0.0; 3];
        for (p, b) in part.iter_mut().zip(r.u.all()) {
            *p = self.dev.sum_sq(b)?;
        }
        let total: f64 = part.iter().sum();
        let multigrid = s.momentum_sweep == MomentumSweep::Multigrid;
        let first = {
            let cw = self.coupled.as_mut().expect("coupled work");
            std::mem::replace(&mut cw.mg_first, false)
        };
        for (i, count) in n.momentum.iter_mut().enumerate() {
            let (ri, zi) = (r.u.all()[i], z.u.all()[i]);
            self.dev.fill(zi, 0.0)?;
            if part[i] <= NEGLIGIBLE * total {
                continue;
            }
            self.momentum_matrix(i, sym)?;
            if multigrid && self.momentum_hierarchy(i)? {
                self.dev
                    .multigrid_cycles(w.val, ri, zi, s.momentum_cycles)?;
                *count += s.momentum_cycles;
                let kept = !self.coupled.as_ref().expect("coupled work").mg_fresh[i];
                if first && kept {
                    // a kept hierarchy must still reduce the residual
                    self.dev.matrix_product(w.val, zi, t)?;
                    self.dev.lincomb(t, 1.0, ri, -1.0, t)?;
                    let left = self.dev.sum_sq(t)?;
                    // (a NaN counts as not reduced)
                    let reduced = left <= STALE * STALE * part[i];
                    if !reduced {
                        self.sweep.slot[i] = None;
                        self.sweep_release();
                        if self.momentum_hierarchy(i)? {
                            self.dev
                                .multigrid_cycles(w.val, ri, zi, s.momentum_cycles)?;
                            *count += s.momentum_cycles;
                        } else {
                            self.dev.fill(zi, 0.0)?;
                            let (it, _) = self.dev.solve(
                                w.val,
                                ri,
                                zi,
                                Preconditioner::Jacobi,
                                false,
                                s.momentum_tolerance,
                                s.momentum_max_iter,
                            )?;
                            *count += it;
                        }
                    }
                }
            } else {
                let (it, _) = self.dev.solve(
                    w.val,
                    ri,
                    zi,
                    Preconditioner::Jacobi,
                    false,
                    s.momentum_tolerance,
                    s.momentum_max_iter,
                )?;
                *count += it;
            }
        }
        if multigrid {
            self.dev.amg_select(0)?;
        }
        self.dev.sync()?;
        // the setups are counted under the multigrid setup
        self.timings.momentum_solve +=
            tm.elapsed().as_secs_f64() - (self.timings.amg_setup - setup);
        let tm = Instant::now();
        self.dev.delta_flux(z.u, mf)?;
        self.dev.net_outflow(mf, t)?;
        self.dev.lincomb(t, 1.0 / sc, r.p, -1.0, t)?;
        self.dev.fill(z.p, 0.0)?;
        let (it, _) = self.dev.solve(
            w.pval,
            t,
            z.p,
            self.settings.pressure.preconditioner,
            true,
            s.pressure_tolerance,
            s.pressure_max_iter,
        )?;
        n.pressure += it;
        if let Some(rc) = self.reference_cell {
            let shift = self.dev.value(z.p, rc)?;
            self.dev.subtract(z.p, shift)?;
        }
        self.dev.gradient(
            Lsq::Pressure,
            z.p,
            Mask::PressureCorrection,
            w.zero_boundary,
            g,
        )?;
        self.dev.correct_velocity(dhat, g, z.u)?;
        self.dev.sync()?;
        self.timings.pressure_solve += tm.elapsed().as_secs_f64();
        Ok(())
    }

    /// Flexible GMRES on K x = r from x = 0 (right preconditioning,
    /// modified Gram–Schmidt, Givens rotations); the solution in `cw.x`.
    fn coupled_solve(&mut self, sc: f64, sym: bool) -> Result<Counts> {
        let (r, x, wv) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (cw.r, cw.x, cw.w)
        };
        let s = self.settings.coupled.clone();
        let m = s.max_iter;
        let mut n = Counts::default();
        for b in x.all() {
            self.dev.fill(b, 0.0)?;
        }
        let beta = self.norm4(r)?;
        if beta == 0.0 || !beta.is_finite() {
            n.reduction = if beta == 0.0 { 0.0 } else { f64::NAN };
            return Ok(n);
        }
        let (v0, _) = self.krylov_vectors(0)?;
        self.lincomb4(v0, 1.0 / beta, r, 0.0, r)?;
        // Hessenberg matrix by columns, rotations, right-hand side
        let mut h = vec![vec![0.0; m + 1]; m];
        let (mut cs, mut sn) = (vec![0.0; m], vec![0.0; m]);
        let mut gk = vec![0.0; m + 1];
        gk[0] = beta;
        let mut kk = 0;
        for j in 0..m {
            let (vj, zj) = self.krylov_vectors(j)?;
            self.coupled_preconditioner(vj, zj, sc, sym, &mut n)?;
            self.coupled_operator(zj, wv, sc, sym)?;
            for (i, hi) in h[j].iter_mut().enumerate().take(j + 1) {
                let vi = self.coupled.as_ref().expect("coupled work").v[i];
                *hi = self.dot4(wv, vi)?;
                self.lincomb4(wv, 1.0, wv, -*hi, vi)?;
            }
            let hn = self.norm4(wv)?;
            h[j][j + 1] = hn;
            for i in 0..j {
                let (a, b) = (h[j][i], h[j][i + 1]);
                h[j][i] = cs[i] * a + sn[i] * b;
                h[j][i + 1] = -sn[i] * a + cs[i] * b;
            }
            let (a, b) = (h[j][j], h[j][j + 1]);
            let d = a.hypot(b);
            if d == 0.0 || !d.is_finite() {
                break;
            }
            cs[j] = a / d;
            sn[j] = b / d;
            h[j][j] = d;
            h[j][j + 1] = 0.0;
            gk[j + 1] = -sn[j] * gk[j];
            gk[j] *= cs[j];
            kk = j + 1;
            n.krylov = kk;
            n.reduction = gk[j + 1].abs() / beta;
            if n.reduction <= s.tolerance || hn == 0.0 || kk == m {
                break;
            }
            let (vn, _) = self.krylov_vectors(j + 1)?;
            self.lincomb4(vn, 1.0 / hn, wv, 0.0, wv)?;
        }
        // back substitution, x = Σ y_j z_j
        let mut y = vec![0.0; kk];
        for i in (0..kk).rev() {
            let mut a = gk[i];
            for j in i + 1..kk {
                a -= h[j][i] * y[j];
            }
            y[i] = a / h[i][i];
        }
        for (j, yj) in y.iter().enumerate() {
            let zj = self.coupled.as_ref().expect("coupled work").z[j];
            self.lincomb4(x, 1.0, x, *yj, zj)?;
        }
        Ok(n)
    }

    /// The pressure matrix Ŝ (in `w.pval`) differs from the one of the last
    /// multigrid setup by more than `reuse`: with D̂ it changes in every
    /// outer iteration, and a hierarchy from an older matrix needs many
    /// more CG iterations (requirements step 22, section 7.9).
    fn pressure_changed(&mut self, reuse: f64) -> Result<bool> {
        if self.settings.pressure.preconditioner != Preconditioner::Amg
            || reuse >= 1.0
            || self.amg_matrix.is_none()
        {
            return Ok(false);
        }
        let cur = self.dev.download(self.w.pval)?;
        let old = self.amg_matrix.as_deref().unwrap_or_default();
        Ok(old.len() != cur.len() || matrix_change(&cur, old, None, &[]) > reuse)
    }

    /// The length L of the automatic pseudo-time step: the volume of the
    /// domain over the area of its boundary faces without the symmetry
    /// planes (periodic faces are interior faces); (ΣV)^(1/3) for a domain
    /// without such boundaries.
    pub(crate) fn pseudo_length(&mut self) -> f64 {
        if let Some(l) = self.pseudo_length {
            return l;
        }
        let d = &self.disc;
        let volume: f64 = d.cell_volume.iter().sum();
        let mut area = 0.0;
        for (fb, kind) in self.bc.mom_kind.iter().enumerate() {
            if *kind != k::MOM_SYMMETRY {
                let f = d.n_interior + fb;
                area += (d.face_area.x[f].powi(2)
                    + d.face_area.y[f].powi(2)
                    + d.face_area.z[f].powi(2))
                .sqrt();
            }
        }
        let l = if area > 0.0 {
            volume / area
        } else {
            volume.cbrt()
        };
        self.pseudo_length = Some(l);
        l
    }

    /// The largest prescribed boundary speed: the faces that keep their
    /// velocity (inlets, moving walls). Symmetry planes and outlets take
    /// their velocity from the cells and do not count; on a mesh of one
    /// layer the symmetry planes would carry the whole field.
    fn prescribed_speed(&mut self) -> f64 {
        if let Some(u) = self.pseudo_speed {
            return u;
        }
        let b = &self.bc;
        let mut sq: f64 = 0.0;
        for (fb, rule) in b.ub_rule.iter().enumerate() {
            if *rule == k::UB_FIXED {
                sq = sq.max(b.ub[0][fb].powi(2) + b.ub[1][fb].powi(2) + b.ub[2][fb].powi(2));
            }
        }
        let u = sq.sqrt();
        self.pseudo_speed = Some(u);
        u
    }

    /// The pseudo-time step Δτ of this outer iteration (`[solver.pseudo_time]`):
    /// the given one, else `factor` × L / U with L of
    /// [`pseudo_length`](Self::pseudo_length) and U the larger of the
    /// largest prescribed boundary speed
    /// ([`prescribed_speed`](Self::prescribed_speed)) and the root mean
    /// square of the cell speeds; none in a fluid at rest.
    pub(crate) fn pseudo_step(&mut self) -> Result<Option<f64>> {
        let s = self.settings.pseudo_time.clone();
        if let Some(dt) = s.step {
            return Ok(Some(dt));
        }
        let length = self.pseudo_length();
        let bnd = self.prescribed_speed();
        let mut cells = 0.0;
        for u in self.f.u.all() {
            cells += self.dev.sum_sq(u)?;
        }
        // with buoyancy also the buoyancy velocity (stage 7b1): a fluid at
        // rest in a stratification has no other speed
        let buoyant = self.buoyancy.as_ref().map_or(0.0, |b| b.speed(length));
        let speed = (cells / self.disc.n_cells as f64)
            .sqrt()
            .max(bnd)
            .max(buoyant);
        Ok((speed > 0.0 && speed.is_finite()).then(|| s.factor * length / speed))
    }

    /// Hands the turbulence equations their pseudo-time step of this outer
    /// iteration, min(Δτ, c_T T) per cell (`step`: Δτ where the iteration
    /// has computed it already).
    pub(crate) fn turbulence_step(&mut self, step: Option<Option<f64>>) -> Result<()> {
        if self.turbulence.is_none() {
            return Ok(());
        }
        let c = self
            .settings
            .pseudo_time
            .turbulence_factor(self.settings.algorithm(false));
        let pseudo = if c > 0.0 {
            let step = match step {
                Some(s) => s,
                None => self.pseudo_step()?,
            };
            Some(PseudoStep {
                inv_dt: step.map_or(0.0, |dt| 1.0 / dt),
                rate_factor: 1.0 / c,
            })
        } else {
            None
        };
        if let Some(t) = self.turbulence.as_mut() {
            t.pseudo = pseudo;
        }
        Ok(())
    }

    /// One steady outer iteration of the coupled algorithm; returns the
    /// residuals of the fields before the solve.
    pub(crate) fn iterate_coupled(&mut self) -> Result<Residuals> {
        self.iteration += 1;
        let mut res = Residuals {
            iteration: self.iteration,
            ..Default::default()
        };
        let s = self.settings.clone();
        let alpha = s.alpha_u(false);
        self.coupled_work()?;
        let (r, x, t, mf, cf, pval, dsum, one, area, dhat, y, zones) = {
            let cw = self.coupled.as_ref().expect("coupled work");
            (
                cw.r, cw.x, cw.t, cw.mf, cw.cf, cw.pval, cw.dsum, cw.one, cw.area, cw.dhat, cw.y,
                cw.zones,
            )
        };
        let tm = Instant::now();
        self.update_boundaries()?;
        self.gradients()?;
        let step = self.pseudo_step()?;
        let dt = step.filter(|_| s.pseudo_time.momentum);
        res.pseudo_step = dt;
        let tau = dt.map_or(0.0, |dt| self.phys.rho / dt);
        self.dev.sync()?;
        self.timings.gradients += tm.elapsed().as_secs_f64();

        // Y = ṁ − ρ ū·S of the current fields: the relaxation term of the
        // Rhie–Chow flux is (1 − α) Y, its pseudo-time term γ_f Y (D27 with
        // the current iterate as the old level)
        let tm = Instant::now();
        let (f, w) = (self.f, self.w);
        let a = self.flux_args(None, w.zero_face, None);
        self.dev.mass_flux(&a, y)?;
        self.dev.scaled_difference(y, 1.0, f.mass_flux, y)?;
        self.dev.lincomb(w.relax, 1.0 - alpha, y, 0.0, y)?;
        // momentum rows and the residual of the current fields
        let time = dt.map(|_| TimeTerms {
            a: [1.0, 1.0, 0.0],
            dt: dt.unwrap_or(1.0),
            u_old: f.u,
            u_old2: f.u,
            explicit: (0.0, f.u),
        });
        let conv = self.convecting_flux()?;
        let scale = self.momentum_assembly(conv, alpha, 1.0, time)?;
        let sym = self.bc.mom_kind.contains(&k::MOM_SYMMETRY);
        for (i, residual) in res.momentum.iter_mut().enumerate() {
            self.momentum_matrix(i, sym)?;
            self.dev.matrix_product(w.val, f.u.all()[i], t)?;
            self.dev
                .lincomb(r.u.all()[i], 1.0, w.rhs.all()[i], -1.0, t)?;
            *residual = self.dev.sum_sq(r.u.all()[i])?.sqrt() / scale;
        }
        self.dev.sync()?;
        self.timings.momentum_assembly += tm.elapsed().as_secs_f64();

        // Rhie–Chow fluxes of the current fields, the continuity residual,
        // the matrix L and the matrix Ŝ of the preconditioner
        let tm = Instant::now();
        let a = self.flux_args(Some(w.d), w.relax, dt.map(|_| (w.ds, y, alpha)));
        self.dev.mass_flux(&a, w.m)?;
        self.close_backflow(w.m)?;
        res.continuity = self.continuity(w.m)?;
        self.dev.pressure_faces(w.d, cf)?;
        self.dev.pressure_rows(cf, w.m, pval, r.p)?;
        match zones {
            // several fluids: 1/(a/D^s + ρ b) = (1/ρ) / (a/(ρ D^s) + b),
            // and the bound D^C as ρ D^C before the division
            Some((rho, inv, work)) if tau > 0.0 => {
                let b = tau / (self.phys.rho * alpha);
                self.dev.product(work, 1.0, rho, w.ds)?;
                self.dev.product(dhat, 1.0, rho, w.dc)?;
                self.dev.schur_diagonal(dhat, 1.0 - alpha, b, work, dhat)?;
                self.dev.product(dhat, 1.0, inv, dhat)?;
            }
            _ => self
                .dev
                .schur_diagonal(dhat, 1.0 - alpha, tau / alpha, w.ds, w.dc)?,
        }
        self.dev.lincomb(dsum, 1.0, w.d, 1.0, dhat)?;
        self.dev.pressure_faces(dsum, w.cf)?;
        self.dev.pressure_rows(w.cf, w.m, w.pval, w.prhs)?;
        if let Some(rc) = self.reference_cell {
            let pos = self.disc.sell.pos(rc, 0);
            self.dev.scale_entry(pval, pos, 2.0)?;
            self.dev.scale_entry(w.pval, pos, 2.0)?;
        }
        self.dev.sync()?;
        self.timings.pressure_assembly += tm.elapsed().as_secs_f64();
        let rebuild = (self.iteration - 1).is_multiple_of(s.amg_rebuild);
        let changed = self.pressure_changed(s.coupled.pressure_reuse)?;
        self.pressure_amg(rebuild || changed)?;

        // the continuity row in the units of the momentum rows: times the
        // velocity Σ a_P / (ρ Σ|S|), the mean diagonal coefficient over the
        // face areas (the mean normal speed where convection dominates,
        // ν/h in a fluid at rest)
        let sc = self.dev.dot(w.diag.x, one)? / area;
        self.dev.lincomb(r.p, sc, r.p, 0.0, r.p)?;
        let tm = Instant::now();
        self.momentum_hierarchies_check(sym)?;
        self.timings.amg_setup += tm.elapsed().as_secs_f64();
        let n = self.coupled_solve(sc, sym)?;
        res.momentum_iterations = n.momentum;
        res.pressure_iterations = n.pressure;
        res.coupled_iterations = n.krylov;

        // the increments
        let tm = Instant::now();
        if let Some(rc) = self.reference_cell {
            let shift = self.dev.value(x.p, rc)?;
            self.dev.subtract(x.p, shift)?;
        }
        for (u, d) in f.u.all().into_iter().zip(x.u.all()) {
            self.dev.lincomb(u, 1.0, u, 1.0, d)?;
        }
        self.dev.correct_fluxes(cf, w.m, x.p, w.mc)?;
        self.dev.delta_flux(x.u, mf)?;
        self.dev.lincomb(w.mc, 1.0, w.mc, 1.0, mf)?;
        self.dev.lincomb(f.p, 1.0, f.p, 1.0, x.p)?;
        // one pressure correction on the remaining continuity defect
        self.dev.pressure_rows(w.cf, w.mc, w.pval, w.prhs)?;
        if let Some(rc) = self.reference_cell {
            let pos = self.disc.sell.pos(rc, 0);
            self.dev.scale_entry(w.pval, pos, 2.0)?;
        }
        self.dev.fill(w.pc, 0.0)?;
        let (it, _) = self.dev.solve(
            w.pval,
            w.prhs,
            w.pc,
            s.pressure.preconditioner,
            true,
            s.pressure_tol(false),
            s.pressure.max_iter,
        )?;
        res.pressure_iterations += it;
        if let Some(rc) = self.reference_cell {
            let shift = self.dev.value(w.pc, rc)?;
            self.dev.subtract(w.pc, shift)?;
        }
        self.dev.correct_fluxes(w.cf, w.mc, w.pc, w.m)?;
        self.dev.copy(w.m, w.mc)?;
        self.dev.gradient(
            Lsq::Pressure,
            w.pc,
            Mask::PressureCorrection,
            w.zero_boundary,
            w.gpc,
        )?;
        self.dev.correct_velocity(dhat, w.gpc, f.u)?;
        self.dev.lincomb(f.p, 1.0, f.p, 1.0, w.pc)?;
        self.drive_mass_flow()?;
        self.dev.sync()?;
        self.timings.pressure_assembly += tm.elapsed().as_secs_f64();
        if !res.finite() || !n.reduction.is_finite() {
            self.diverged = Some(res.clone());
            return Err(crate::diag::err(crate::diag::Diagnostic::new(
                crate::diag::RUN_DIVERGED,
                format!(
                    "the solution diverged in outer iteration {} (residuals {:?}, continuity {:.3e})",
                    self.iteration, res.momentum, res.continuity
                ),
            )));
        }
        std::mem::swap(&mut self.f.mass_flux, &mut self.w.mc);
        if let Some(t) = self.turbulence.as_mut() {
            t.amg_rebuild = rebuild;
        }
        self.turbulence_step(Some(step))?;
        self.solve_scalars(None, &mut res)?;
        Ok(res)
    }
}
