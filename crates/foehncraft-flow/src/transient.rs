//! Transient incompressible flow: the PISO algorithm (Issa 1986) with
//! backward differences of first and second order (BDF1, BDF2 with variable
//! time steps) or Crank–Nicolson (D28), on the state of [`Flow`].
//!
//! One time step from tⁿ to tⁿ⁺¹ = tⁿ + Δt:
//!
//! 1. the old levels move on (uⁿ⁻¹ ← uⁿ, ṁⁿ⁻¹ ← ṁⁿ); for Crank–Nicolson the
//!    spatial operator R(uⁿ) of the old level is evaluated with the old flux;
//! 2. per outer iteration (one for plain PISO):
//!    - boundary values (the prescribed ones of tⁿ⁺¹), gradients, limiters;
//!    - the momentum predictor: the rows of [`crate::simplec`] with the time
//!      term ρV a₀/Δt on the diagonal and the old levels on the right-hand
//!      side, convected by the flux of the latest iterate;
//!    - in the first outer iteration of a second-order step the explicit
//!      quantities are extrapolated to tⁿ⁺¹ from the levels n and n−1
//!      (φ* = (1 + ω)φⁿ − ωφⁿ⁻¹): the convecting flux and the velocity
//!      gradients and limiters of the deferred corrections, so that plain
//!      PISO keeps second order (`extrapolation`; a steady state is not
//!      changed by it);
//!    - the predicted Rhie–Chow flux with the time-step term γ_f Y_f of
//!      `foehncraft_kernels::flow::face_mass_flux` (D27: a steady state of the
//!      time stepping is exactly the steady solution, for any Δt);
//!    - `correctors` pressure corrections: Σ c_f (p'_P − p'_N) = −Σ ±ṁ*
//!      with c_f from a coefficient D, the fluxes corrected to continuity,
//!      u' = δ − D∇p', p += p'; between two correctors the neighbour
//!      correction δ (the part the corrector neglects) is added to u and its
//!      flux to ṁ. `solver.corrector`: Issa's D = V/(a_P/α) with
//!      δ = H(u')/a_P (plain PISO by default), or the SIMPLEC form
//!      D^C = V/(a_P/α − Σ|a_N|) (Van Doormaal & Raithby 1984) with
//!      δ = Σ_N |a_N| (u'_N − u'_P)/a_P, since D^C already counts
//!      Σ|a_N| u'_P (by default with outer iterations). Where a_P is almost
//!      all neighbour coupling (a fine pitch about an axis at high
//!      viscosity: the thin wedge cells, whose a_P is 50 times the time
//!      term in gallery 103) Issa's D underestimates the velocity response
//!      to the pressure by a_P/(a_P − Σ|a_N|); with two correctors the
//!      pressure there grew by a factor of 2.4 in every outer iteration. In
//!      plain PISO Issa's form damps more: Crank–Nicolson at Courant number
//!      4 settles with it and not with the SIMPLEC form;
//! 3. boundary values and gradients of the new level, the old-level
//!    difference Xⁿ⁺¹ = ṁ − ρ ū·S of the Rhie–Chow time-step term and the
//!    largest Courant number.
//!
//! Backward differences: du/dt ≈ (a₀ uⁿ⁺¹ + a₁ uⁿ + a₂ uⁿ⁻¹)/Δt with, for
//! ω = Δtₙ/Δtₙ₋₁, a₀ = (1 + 2ω)/(1 + ω), a₁ = −(1 + ω), a₂ = ω²/(1 + ω)
//! (BDF2; the first step of a run is BDF1: a₀ = 1, a₁ = −1, a₂ = 0).
//! Crank–Nicolson: a₀ = 1, a₁ = −1, θ = ½ on the spatial terms; the pressure
//! is then the one of tⁿ⁺¹ᐟ². Its first step is a BDF1 step as well:
//! Crank–Nicolson does not damp the part of the initial fields that is not
//! consistent with the discrete equations (the face fluxes, the pressure),
//! which otherwise leaves an error of first order in Δt (the start of
//! Rannacher 1984; measured with the Taylor–Green vortex).

use crate::device::{Buf, Buf3, TimeTerms};
use crate::simplec::{Flow, Residuals};
use anyhow::{Result, ensure};
use foehncraft_fv::solver::Scheme;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// A module that takes part in the outer iterations of a time step
/// ([`Transient::step_with`]).
pub trait OuterHook {
    /// The most outer iterations it may ask for.
    fn max_outer(&self) -> usize;
    /// Before the momentum assembly of outer iteration `k` (the gradients
    /// of the current iterate are set).
    fn before_momentum(&mut self, fl: &mut Flow, k: usize) -> Result<()>;
    /// After the solves of outer iteration `k`: whether it wants another
    /// outer iteration.
    fn after_solve(&mut self, fl: &mut Flow, k: usize) -> Result<bool>;
}

/// Time integration (D28).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum TimeScheme {
    /// First-order backward differences (backward Euler).
    Bdf1,
    /// Second-order backward differences with variable steps (default);
    /// the first step is BDF1.
    #[default]
    Bdf2,
    /// Crank–Nicolson, second order; the first step is BDF1.
    CrankNicolson,
}

impl TimeScheme {
    /// Name for the log.
    pub fn label(self) -> &'static str {
        match self {
            Self::Bdf1 => "BDF1 (backward Euler)",
            Self::Bdf2 => "BDF2",
            Self::CrankNicolson => "Crank–Nicolson",
        }
    }
}

/// The clock of a transient run (part of the restart).
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq)]
pub struct Clock {
    /// Time of the current level (s).
    pub time: f64,
    /// Time steps done.
    pub step: usize,
    /// Size of the last step (s), 0 before the first.
    pub dt_prev: f64,
    /// Steps behind the current level, up to 2: 0 at the start, 1 after
    /// the first step (an old level exists: BDF2, Crank–Nicolson), 2 later
    /// (the old level is itself a solution of the discrete equations: the
    /// extrapolation, which would otherwise carry an initial field that is
    /// not divergence-free into the convecting flux).
    pub levels: usize,
    /// Size of the step before the last one (s), 0 while uⁿ⁻² is not a
    /// level of this run (the first two steps, after a rewind): the
    /// predictor of the error estimate (D65) needs it.
    #[serde(default)]
    pub dt_prev2: f64,
}

/// Report of one time step.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct StepReport {
    /// Step number (from 1).
    pub step: usize,
    /// Time at the end of the step (s).
    pub time: f64,
    /// Time step (s).
    pub dt: f64,
    /// Residuals of the first outer iteration (momentum before the solve,
    /// continuity of the predicted fluxes); pressure iterations of all
    /// correctors and outer iterations.
    pub residuals: Residuals,
    /// Residuals at the start of the last outer iteration.
    pub last: Residuals,
    /// Outer iterations done.
    pub outer: usize,
    /// Largest Courant number of the new level.
    pub cfl: f64,
    /// Estimated local error of the time integration relative to the
    /// tolerance (D65; > 1: the step is too long), when the estimate is on
    /// and three levels exist.
    pub error: Option<f64>,
}

/// Tolerances of the error-based time step (D65).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ErrorControl {
    /// Relative tolerance on the velocity.
    pub relative: f64,
    /// Absolute tolerance on the velocity (m/s).
    pub absolute: f64,
}

/// Work of the error estimate (D65): the predictor and what a rejected
/// step needs to be repeated exactly.
#[derive(Debug, Clone, Copy)]
struct ErrorWork {
    /// Quadratic extrapolation of uⁿ⁻², uⁿ⁻¹, uⁿ to tⁿ⁺¹.
    pred: Buf3,
    /// uⁿ⁻² at the start of the step.
    u2: Buf3,
    /// Xⁿ⁻¹ at the start of the step.
    x_old: Buf,
    /// Gⁿ⁻¹ at the start of the step.
    g_old: [Buf; N_EXPLICIT],
    /// The current explicit buffers and those of Gⁿ⁻¹ at the start.
    ids: [Buf; N_EXPLICIT],
    old_ids: [Buf; N_EXPLICIT],
    /// The clock at the start.
    clock: Clock,
}

/// Old levels and work buffers of the time integration.
#[derive(Debug, Clone, Copy)]
struct Levels {
    u_old: Buf3,
    u_old2: Buf3,
    m_old: Buf,
    m_old2: Buf,
    x_cur: Buf,
    x_old: Buf,
    y: Buf,
    conv: Buf,
    r_old: Buf3,
    ucorr: Buf3,
    delta: Buf3,
    dflux: Buf,
    cfl: Buf,
    /// Velocity gradients and limiters of level n−1.
    g_old: [Buf; N_EXPLICIT],
    /// Boundary velocities of levels n and n−1 of the convective outlets
    /// (per boundary face; only with convective outlets).
    ub_old: Option<(Buf3, Buf3)>,
}

/// Number of cell fields in [`explicit_buffers`].
const N_EXPLICIT: usize = 12;

/// The velocity gradients and limiters of `f` (the explicit parts of the
/// deferred corrections).
fn explicit_buffers(f: &crate::simplec::Fields) -> [Buf; N_EXPLICIT] {
    let g = &f.grad_u;
    [
        g[0].x, g[0].y, g[0].z, g[1].x, g[1].y, g[1].z, g[2].x, g[2].y, g[2].z, f.psi.x, f.psi.y,
        f.psi.z,
    ]
}

fn set_explicit_buffers(f: &mut crate::simplec::Fields, b: [Buf; N_EXPLICIT]) {
    let b3 = |i: usize| Buf3 {
        x: b[i],
        y: b[i + 1],
        z: b[i + 2],
    };
    f.grad_u = [b3(0), b3(3), b3(6)];
    f.psi = b3(9);
}

/// The transient driver on a [`Flow`].
pub struct Transient {
    /// Time integration.
    pub scheme: TimeScheme,
    /// The clock.
    pub clock: Clock,
    /// Error estimate of the time integration (D65), off with `None`.
    pub error_control: Option<ErrorControl>,
    /// Residuals of every outer iteration of the last step (for the log;
    /// no result depends on them).
    pub outer_log: Vec<Residuals>,
    l: Levels,
    e: Option<ErrorWork>,
    /// The sliding mesh (D182): how the old levels of the time terms move.
    content: Option<Content>,
}

/// The sliding mesh moves the cells of the rotor, so the time term of a
/// moving cell must see what the cell held at the old times, and not what
/// the cell of the mesh that is now at its place held: the old levels of the
/// time terms are the Eulerian levels (the fields at the places of the mesh)
/// with the data of the rotor cells shifted along the turn, one step for the
/// level n and two for the level n−1 (D182).
#[derive(Debug, Clone)]
struct Content {
    /// Entry i of the shifted array is entry `perm1[i]` of the old one
    /// (one step), `perm2` for two steps.
    perm1: Vec<usize>,
    perm2: Vec<usize>,
    /// The shifted velocities of the levels n and n−1.
    c1: Buf3,
    c2: Buf3,
}

/// Host copy of the old levels (for a restart).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OldLevels {
    /// uⁿ⁻¹ per component.
    pub u_old: [Vec<f64>; 3],
    /// ṁⁿ⁻¹.
    pub m_old: Vec<f64>,
    /// Xⁿ = ṁⁿ − ρ ūⁿ·S.
    pub x_cur: Vec<f64>,
    /// Xⁿ⁻¹.
    pub x_old: Vec<f64>,
    /// Velocity gradients (9) and limiters (3) of level n−1.
    pub g_old: Vec<Vec<f64>>,
    /// uⁿ⁻¹ of the step before (uⁿ⁻² of the next step), for the predictor
    /// of the error estimate (D65); only with the estimate on.
    #[serde(default)]
    pub u_old2: Option<[Vec<f64>; 3]>,
    /// Boundary velocity u_bⁿ⁻¹ per boundary face, for the backward
    /// difference of the convective outlets (D29); only with them.
    #[serde(default)]
    pub ub_old: Option<[Vec<f64>; 3]>,
}

impl OldLevels {
    /// Bytes held.
    pub fn bytes(&self) -> usize {
        8 * (self.u_old.iter().map(Vec::len).sum::<usize>()
            + self.m_old.len()
            + self.x_cur.len()
            + self.x_old.len()
            + self.g_old.iter().map(Vec::len).sum::<usize>()
            + self
                .u_old2
                .as_ref()
                .map_or(0, |u| u.iter().map(Vec::len).sum::<usize>())
            + self
                .ub_old
                .as_ref()
                .map_or(0, |u| u.iter().map(Vec::len).sum::<usize>()))
    }
}

impl Transient {
    /// Buffers on the device of `fl`, the clock at `start_time`. Call
    /// [`start`](Self::start) for new start fields, or
    /// [`restore`](Self::restore) after a restart.
    pub fn new(fl: &mut Flow, scheme: TimeScheme, start_time: f64) -> Result<Self> {
        ensure!(
            fl.transient,
            "the flow solver was created for a steady case (no [time] table)"
        );
        let nc = fl.disc.n_cells;
        let nf = fl.disc.n_faces;
        let nb = nf - fl.disc.n_interior;
        let has_conv = fl.conv.is_some();
        let dev = &mut fl.dev;
        let ub_old = if has_conv {
            let mut b3 = || -> Result<Buf3> {
                Ok(Buf3 {
                    x: dev.alloc(nb)?,
                    y: dev.alloc(nb)?,
                    z: dev.alloc(nb)?,
                })
            };
            Some((b3()?, b3()?))
        } else {
            None
        };
        let mut c3 = || -> Result<Buf3> {
            Ok(Buf3 {
                x: dev.alloc(nc)?,
                y: dev.alloc(nc)?,
                z: dev.alloc(nc)?,
            })
        };
        let (u_old, u_old2, r_old, ucorr, delta) = (c3()?, c3()?, c3()?, c3()?, c3()?);
        let l = Levels {
            u_old,
            u_old2,
            r_old,
            ucorr,
            delta,
            m_old: dev.alloc(nf)?,
            m_old2: dev.alloc(nf)?,
            x_cur: dev.alloc(nf)?,
            x_old: dev.alloc(nf)?,
            y: dev.alloc(nf)?,
            conv: dev.alloc(nf)?,
            dflux: dev.alloc(nf)?,
            cfl: dev.alloc(nc)?,
            g_old: {
                let mut g = [Buf(0); N_EXPLICIT];
                for b in &mut g {
                    *b = dev.alloc(nc)?;
                }
                g
            },
            ub_old,
        };
        Ok(Self {
            scheme,
            clock: Clock {
                time: start_time,
                ..Default::default()
            },
            error_control: None,
            outer_log: Vec::new(),
            l,
            e: None,
            content: None,
        })
    }

    /// The time terms see the old levels with the data of the cells moved
    /// as `perm` says (entry i of the moved array is entry `perm[i]` of
    /// the Eulerian one), once for the level n and twice for the level
    /// n−1: the sliding mesh of D182 ([`crate::sliding`]).
    pub fn set_sliding(&mut self, fl: &mut Flow, perm: Vec<usize>) -> Result<()> {
        ensure!(
            perm.len() == fl.disc.n_cells,
            "the shift of the cells has {} entries for {} cells",
            perm.len(),
            fl.disc.n_cells
        );
        let nc = fl.disc.n_cells;
        let perm2: Vec<usize> = perm.iter().map(|&i| perm[i]).collect();
        let mut c3 = || -> Result<Buf3> {
            Ok(Buf3 {
                x: fl.dev.alloc(nc)?,
                y: fl.dev.alloc(nc)?,
                z: fl.dev.alloc(nc)?,
            })
        };
        let (c1, c2) = (c3()?, c3()?);
        self.content = Some(Content {
            perm1: perm,
            perm2,
            c1,
            c2,
        });
        Ok(())
    }

    /// Boundary values, gradients and Xⁿ of new start fields.
    pub fn start(&mut self, fl: &mut Flow) -> Result<()> {
        self.finish_level(fl)
    }

    /// Boundary values, gradients and Xⁿ of the current fields.
    fn finish_level(&mut self, fl: &mut Flow) -> Result<()> {
        fl.update_boundaries()?;
        fl.gradients()?;
        let a = fl.flux_args(None, fl.w.zero_face, None);
        fl.dev.mass_flux(&a, self.l.x_cur)?;
        fl.dev
            .scaled_difference(self.l.x_cur, 1.0, fl.f.mass_flux, self.l.x_cur)
    }

    /// The old levels on the host.
    pub fn old_levels(&self, fl: &mut Flow) -> Result<OldLevels> {
        let l = self.l;
        Ok(OldLevels {
            u_old: [
                fl.dev.download(l.u_old.x)?,
                fl.dev.download(l.u_old.y)?,
                fl.dev.download(l.u_old.z)?,
            ],
            m_old: fl.dev.download(l.m_old)?,
            x_cur: fl.dev.download(l.x_cur)?,
            x_old: fl.dev.download(l.x_old)?,
            g_old: l
                .g_old
                .iter()
                .map(|&b| fl.dev.download(b))
                .collect::<Result<_>>()?,
            u_old2: if self.error_control.is_some() && self.clock.dt_prev2 > 0.0 {
                Some([
                    fl.dev.download(l.u_old2.x)?,
                    fl.dev.download(l.u_old2.y)?,
                    fl.dev.download(l.u_old2.z)?,
                ])
            } else {
                None
            },
            ub_old: match l.ub_old {
                Some((b, _)) => Some([
                    fl.dev.download(b.x)?,
                    fl.dev.download(b.y)?,
                    fl.dev.download(b.z)?,
                ]),
                None => None,
            },
        })
    }

    /// Restores the old levels and the clock (after [`Flow::load`]).
    pub fn restore(&mut self, fl: &mut Flow, o: &OldLevels, clock: Clock) -> Result<()> {
        let l = self.l;
        fl.dev.upload(l.u_old.x, &o.u_old[0])?;
        fl.dev.upload(l.u_old.y, &o.u_old[1])?;
        fl.dev.upload(l.u_old.z, &o.u_old[2])?;
        fl.dev.upload(l.m_old, &o.m_old)?;
        fl.dev.upload(l.x_cur, &o.x_cur)?;
        fl.dev.upload(l.x_old, &o.x_old)?;
        ensure!(
            o.g_old.len() == l.g_old.len(),
            "the restart file holds {} old gradient fields, expected {}",
            o.g_old.len(),
            l.g_old.len()
        );
        for (&b, g) in l.g_old.iter().zip(&o.g_old) {
            fl.dev.upload(b, g)?;
        }
        self.clock = clock;
        match &o.u_old2 {
            Some(u) => {
                for (&b, v) in l.u_old2.all().iter().zip(u) {
                    fl.dev.upload(b, v)?;
                }
            }
            // no uⁿ⁻²: the predictor waits a step
            None => self.clock.dt_prev2 = 0.0,
        }
        if let Some((b, _)) = l.ub_old {
            let Some(u) = &o.ub_old else {
                anyhow::bail!(
                    "the restart file holds no old boundary velocity of the convective outlets \
                     (written by a case without them)"
                );
            };
            for (&bi, v) in b.all().iter().zip(u) {
                fl.dev.upload(bi, v)?;
            }
        }
        Ok(())
    }

    /// After a step that failed: the state of level n back (u, p and ṁ
    /// from [`Flow::rewind`]), the old velocity and flux of level n−1 back
    /// into the old levels, and Xⁿ, the boundary values and gradients
    /// recomputed. The clock still stands at level n. A restart file of
    /// this state continues the run, though not bit-identically: the old
    /// gradients of the deferred corrections are those of level n (D63).
    pub fn rewind(&mut self, fl: &mut Flow) -> Result<bool> {
        if !fl.rewind()? {
            return Ok(false);
        }
        let l = self.l;
        self.hold_convective(fl)?;
        if self.clock.levels >= 1 {
            for (a, b) in l.u_old2.all().into_iter().zip(l.u_old.all()) {
                fl.dev.copy(a, b)?;
            }
            fl.dev.copy(l.m_old2, l.m_old)?;
            if let Some(sc) = fl.scalars.as_mut() {
                sc.rewind_levels(fl.dev.as_mut())?;
            }
            if let Some((sc, dev)) = crate::cht::energy_fields(fl) {
                sc.rewind_levels(dev)?;
            }
            if let Some(t) = fl.turbulence.as_mut() {
                t.rewind_levels(fl.dev.as_mut())?;
            }
        }
        // uⁿ⁻² is gone
        self.clock.dt_prev2 = 0.0;
        fl.density_step_back()?;
        self.finish_level(fl)?;
        Ok(true)
    }

    /// After a step was taken back: the convective outlets hold u_bⁿ (the
    /// level the step started from), and their old levels move back.
    fn hold_convective(&mut self, fl: &mut Flow) -> Result<()> {
        if let Some((ub1, ub2)) = self.l.ub_old {
            fl.hold_convective(ub1)?;
            for (a, b) in ub2.all().into_iter().zip(ub1.all()) {
                fl.dev.copy(a, b)?;
            }
        }
        Ok(())
    }

    /// Takes back the step just done, after its error estimate exceeded the
    /// tolerance (D65): every level, the old gradients, the clock and the
    /// state of level n return exactly as before the step, so that the
    /// step repeated with a shorter Δt is the step a run with that Δt
    /// would have taken. The prescribed boundary values of tⁿ are set
    /// beforehand.
    pub fn reject(&mut self, fl: &mut Flow) -> Result<()> {
        let Some(e) = self.e else {
            anyhow::bail!("no step to take back");
        };
        let l = self.l;
        self.clock = e.clock;
        if self.clock.levels >= 1 {
            for (a, b) in l.u_old2.all().into_iter().zip(l.u_old.all()) {
                fl.dev.copy(a, b)?;
            }
            fl.dev.copy(l.m_old2, l.m_old)?;
        }
        for (a, b) in e.u2.all().into_iter().zip(l.u_old2.all()) {
            fl.dev.copy(a, b)?;
        }
        fl.dev.copy(l.x_old, l.x_cur)?;
        fl.dev.copy(e.x_old, l.x_old)?;
        if let Some(sc) = fl.scalars.as_mut() {
            sc.step_back(fl.dev.as_mut())?;
        }
        if let Some((sc, dev)) = crate::cht::energy_fields(fl) {
            sc.step_back(dev)?;
        }
        if let Some(t) = fl.turbulence.as_mut() {
            t.step_back(fl.dev.as_mut())?;
        }
        set_explicit_buffers(&mut fl.f, e.ids);
        self.l.g_old = e.old_ids;
        for (&a, &b) in e.g_old.iter().zip(&e.old_ids) {
            fl.dev.copy(a, b)?;
        }
        // the convective outlets hold u_bⁿ; u, p, ṁ of level n with their
        // boundary values and gradients
        self.hold_convective(fl)?;
        anyhow::ensure!(fl.rewind()?, "no valid state was kept");
        fl.density_step_back()?;
        fl.dev.sync()
    }

    /// The buffers of the error estimate, allocated at first use.
    fn error_work(&mut self, fl: &mut Flow) -> Result<ErrorWork> {
        if let Some(e) = self.e {
            return Ok(e);
        }
        let (nc, nf) = (fl.disc.n_cells, fl.disc.n_faces);
        let dev = &mut fl.dev;
        let mut c3 = || -> Result<Buf3> {
            Ok(Buf3 {
                x: dev.alloc(nc)?,
                y: dev.alloc(nc)?,
                z: dev.alloc(nc)?,
            })
        };
        let (pred, u2) = (c3()?, c3()?);
        let mut g_old = [Buf(0); N_EXPLICIT];
        for b in &mut g_old {
            *b = fl.dev.alloc(nc)?;
        }
        let e = ErrorWork {
            pred,
            u2,
            x_old: fl.dev.alloc(nf)?,
            g_old,
            ids: [Buf(0); N_EXPLICIT],
            old_ids: [Buf(0); N_EXPLICIT],
            clock: self.clock,
        };
        self.e = Some(e);
        Ok(e)
    }

    /// Largest Courant number of the current fluxes for a step `dt`.
    pub fn courant(&mut self, fl: &mut Flow, dt: f64) -> Result<f64> {
        // the flux that carries the quantities: the relative one in a
        // rotating zone (frames, D55, D182)
        let m = fl.convecting_flux()?;
        fl.dev.courant(m, dt, self.l.cfl)?;
        Ok(fl.dev.min_max(self.l.cfl)?.1)
    }

    /// One time step of size `dt`. The prescribed boundary values of the new
    /// time are set beforehand ([`Flow::set_boundary_velocity`]).
    pub fn step(&mut self, fl: &mut Flow, dt: f64) -> Result<StepReport> {
        self.step_with(fl, dt, None)
    }

    /// [`Transient::step`] with a module that takes part in every outer
    /// iteration (the parcel passes of a two-way coupled spray, D252): it
    /// sets its sources before each momentum assembly and may ask for more
    /// outer iterations than `solver.outer_iterations`.
    pub fn step_with(
        &mut self,
        fl: &mut Flow,
        dt: f64,
        mut hook: Option<&mut dyn OuterHook>,
    ) -> Result<StepReport> {
        ensure!(dt > 0.0 && dt.is_finite(), "time step {dt} must be > 0");
        let s = fl.settings.clone();
        let (l, w, f) = (self.l, fl.w, fl.f);
        let second = self.scheme != TimeScheme::Bdf1 && self.clock.levels >= 1;
        let omega = if second { dt / self.clock.dt_prev } else { 1.0 };
        let (a0, a1, a2, theta) = match self.scheme {
            TimeScheme::CrankNicolson if self.clock.levels >= 1 => (1.0, -1.0, 0.0, 0.5),
            TimeScheme::Bdf2 if second => (
                (1.0 + 2.0 * omega) / (1.0 + omega),
                -(1.0 + omega),
                omega * omega / (1.0 + omega),
                1.0,
            ),
            _ => (1.0, -1.0, 0.0, 1.0),
        };
        let alpha = s.alpha_u(true);
        let alpha_p = s.alpha_p(true);
        let extrapolate = second && s.extrapolation && self.clock.levels >= 2;
        let t = Instant::now();

        // 0. error estimate (D65): the state of the start for a rejection,
        // and the predictor P = L₀uⁿ⁻² + L₁uⁿ⁻¹ + L₂uⁿ, the quadratic
        // through the last three levels at tⁿ⁺¹
        let control = match self.error_control {
            Some(c) if self.scheme != TimeScheme::Bdf1 => Some(c),
            _ => None,
        };
        // (from the fourth step on: uⁿ⁻² is then a level after the first
        // step, whose BDF1 error and the start fields that are not a
        // solution of the discrete equations would count as error)
        let predict = control.is_some()
            && self.clock.levels >= 2
            && self.clock.dt_prev2 > 0.0
            && self.clock.step >= 3;
        if control.is_some() {
            let mut e = self.error_work(fl)?;
            e.ids = explicit_buffers(&fl.f);
            e.old_ids = self.l.g_old;
            e.clock = self.clock;
            for (a, b) in l.u_old2.all().into_iter().zip(e.u2.all()) {
                fl.dev.copy(a, b)?;
            }
            fl.dev.copy(l.x_old, e.x_old)?;
            for (&a, &b) in self.l.g_old.iter().zip(&e.g_old) {
                fl.dev.copy(a, b)?;
            }
            if predict {
                let (h, h1, h2) = (dt, self.clock.dt_prev, self.clock.dt_prev2);
                let c2 = (h + h1) * (h + h1 + h2) / (h1 * (h1 + h2));
                let c1 = -h * (h + h1 + h2) / (h1 * h2);
                let c0 = h * (h + h1) / (h2 * (h1 + h2));
                for ((p, (a, b)), c) in e
                    .pred
                    .all()
                    .into_iter()
                    .zip(f.u.all().into_iter().zip(l.u_old.all()))
                    .zip(l.u_old2.all())
                {
                    fl.dev.lincomb(p, c2, a, c1, b)?;
                    fl.dev.lincomb(p, 1.0, p, c0, c)?;
                }
            }
            self.e = Some(e);
        }

        // 1. old levels
        if self.clock.levels >= 1 {
            for (a, b) in l.u_old.all().into_iter().zip(l.u_old2.all()) {
                fl.dev.copy(a, b)?;
            }
            fl.dev.copy(l.m_old, l.m_old2)?;
        }
        for (a, b) in f.u.all().into_iter().zip(l.u_old.all()) {
            fl.dev.copy(a, b)?;
        }
        fl.dev.copy(f.mass_flux, l.m_old)?;
        // the passive scalars: φⁿ⁻¹ ← φⁿ ← φ; BDF2 (also under
        // Crank–Nicolson, which does not keep them bounded) from the second
        // step, BDF1 in the first and with `bdf1`
        let bdf2 = self.scheme != TimeScheme::Bdf1 && self.clock.levels >= 1;
        if let Some(t) = fl.turbulence.as_mut() {
            t.step_start(
                fl.dev.as_mut(),
                self.clock.levels >= 1,
                bdf2.then(|| dt / self.clock.dt_prev),
            )?;
        }
        let scalar_step = if fl.scalars.is_some() || fl.turbulence.is_some() || fl.cht.is_some() {
            if let Some(sc) = fl.scalars.as_mut() {
                sc.step_start(
                    fl.dev.as_mut(),
                    self.clock.levels >= 1,
                    bdf2.then(|| dt / self.clock.dt_prev),
                )?;
            }
            // the temperature of a conjugate case (stage 7d6): its levels
            // in the energy view, the time terms of the fluid and the
            // solids in one matrix
            if let Some((sc, dev)) = crate::cht::energy_fields(fl) {
                sc.step_start(
                    dev,
                    self.clock.levels >= 1,
                    bdf2.then(|| dt / self.clock.dt_prev),
                )?;
            }
            Some(crate::scalar::StepTime {
                time: self.clock.time + dt,
                dt,
                a: crate::scalar::bdf(second, omega),
            })
        } else {
            None
        };
        // the convective outlets: u_bⁿ⁻¹ ← u_bⁿ ← u_b, then c and β of the
        // backward difference of this step
        if let Some((ub1, ub2)) = l.ub_old {
            for ((a, b), c) in ub1.all().into_iter().zip(ub2.all()).zip(f.ub.all()) {
                fl.dev.copy(a, b)?;
                fl.dev.copy(c, a)?;
            }
            fl.step_convective(ub1, ub2, [a0, a1, a2], dt)?;
        }
        let ew = 1.0 - theta;
        if ew != 0.0 {
            // R(uⁿ) with the old flux, gradients and boundary values
            let m = fl.relative_flux(f.mass_flux)?;
            fl.momentum_assembly_faces(m)?;
            fl.dev.momentum_operator(w.cp, w.cn, w.bf, f.u, l.r_old)?;
        }
        // G* = (1 + ω)Gⁿ − ωGⁿ⁻¹ for the velocity gradients and limiters,
        // into the buffers of Gⁿ⁻¹, which become the current ones; Gⁿ⁻¹ ←
        // Gⁿ. The gradients of level n are those of the end of the last
        // step (with the boundary values of tⁿ), so the extrapolation is
        // consistent also for time-dependent boundary values. The pressure
        // is not extrapolated: p* = 2pⁿ − pⁿ⁻¹ amplifies the pressure modes
        // that the corrections do not remove on skewed meshes (unstable on
        // the cylinder O-grid), and the corrections recover pⁿ⁺¹ anyway.
        let cur = explicit_buffers(&fl.f);
        let old = self.l.g_old;
        if extrapolate {
            for (&a, &b) in cur.iter().zip(&old) {
                fl.dev.lincomb(b, -omega, b, 1.0 + omega, a)?;
            }
            set_explicit_buffers(&mut fl.f, old);
            self.l.g_old = cur;
        } else {
            for (&a, &b) in cur.iter().zip(&old) {
                fl.dev.copy(a, b)?;
            }
        }
        // Y = Σ (−a_k/a₀) X^k
        fl.dev.lincomb(l.y, -a1 / a0, l.x_cur, -a2 / a0, l.x_old)?;
        // the convecting flux of the new level (also for Crank–Nicolson: the
        // implicit half is R(uⁿ⁺¹) with ṁⁿ⁺¹, the explicit one R(uⁿ) with ṁⁿ,
        // as in the later outer iterations)
        // (the relative flux ṁ − ρ S·V of a rotating zone, D182: the frame
        // flux is the same at every level)
        let conv = if extrapolate {
            fl.dev
                .lincomb(l.conv, 1.0 + omega, l.m_old, -omega, l.m_old2)?;
            if let Some(ff) = fl.frame_flux() {
                fl.dev.lincomb(l.conv, 1.0, l.conv, -1.0, ff)?;
            }
            l.conv
        } else {
            fl.relative_flux(f.mass_flux)?
        };
        // the old levels of the time term: those of the cells that moved
        // (the sliding mesh)
        let (u_old, u_old2) = match &self.content {
            Some(c) => {
                for (src, dst, perm) in [(l.u_old, c.c1, &c.perm1), (l.u_old2, c.c2, &c.perm2)] {
                    for (&s, &d) in src.all().iter().zip(&dst.all()) {
                        let v = fl.dev.download(s)?;
                        let m: Vec<f64> = perm.iter().map(|&i| v[i]).collect();
                        fl.dev.upload(d, &m)?;
                    }
                }
                (c.c1, c.c2)
            }
            None => (l.u_old, l.u_old2),
        };
        // (level n−1 of the scalars is shifted in the same way)
        if let (Some(c), Some(sc)) = (&self.content, fl.scalars.as_mut()) {
            sc.shift_old(fl.dev.as_mut(), &c.perm1, &c.perm2)?;
        }
        // (and those of the turbulence variables)
        if let (Some(c), Some(tu)) = (&self.content, fl.turbulence.as_mut()) {
            tu.shift_old(
                fl.dev.as_mut(),
                &c.perm1,
                &c.perm2,
                self.clock.levels >= 1,
                bdf2.then(|| dt / self.clock.dt_prev),
            )?;
        }
        let time = TimeTerms {
            a: [a0, -a1, -a2],
            dt,
            u_old,
            u_old2,
            explicit: (ew, l.r_old),
        };
        // the variable density (stage 7c): ρⁿ⁻¹ ← ρⁿ ← ρ, the mass source
        // of the continuity with the coefficients of this step
        fl.density_step_start(time.a, dt)?;
        fl.dev.sync()?;
        fl.timings.momentum_assembly += t.elapsed().as_secs_f64();

        let mut first = Residuals::default();
        let mut last = Residuals::default();
        let mut outer = 0;
        self.outer_log.clear();
        let rebuild_step = self.clock.step.is_multiple_of(s.amg_rebuild);
        let max_outer = s
            .outer_iterations
            .max(hook.as_ref().map_or(0, |h| h.max_outer()));
        for k in 1..=max_outer {
            outer = k;
            fl.iteration += 1;
            let mut res = Residuals {
                iteration: fl.iteration,
                ..Default::default()
            };
            let t = Instant::now();
            fl.update_boundaries()?;
            // first outer iteration of an extrapolated step: the velocity
            // gradients are G*, the pressure gradient is that of pⁿ
            if !(k == 1 && extrapolate) {
                fl.gradients()?;
            }
            fl.dev.sync()?;
            fl.timings.gradients += t.elapsed().as_secs_f64();
            if let Some(h) = hook.as_deref_mut() {
                h.before_momentum(fl, k)?;
            }
            let t = Instant::now();
            let relax = if alpha < 1.0 {
                fl.relaxation_term(alpha)?;
                w.relax
            } else {
                w.zero_face
            };
            let c = if k == 1 { conv } else { fl.convecting_flux()? };
            // the central correction reads cell velocities, not gradients:
            // extrapolate them as well, u* = (1 + ω)uⁿ − ωuⁿ⁻¹, into the
            // corrector scratch (free until the correctors)
            let central = k == 1 && extrapolate && s.momentum.scheme == Scheme::Central;
            if central {
                for ((d, a), b) in l
                    .delta
                    .all()
                    .into_iter()
                    .zip(l.u_old.all())
                    .zip(l.u_old2.all())
                {
                    fl.dev.lincomb(d, 1.0 + omega, a, -omega, b)?;
                }
                fl.central_u = Some(l.delta);
            }
            let scale = fl.momentum_assembly(c, alpha, theta, Some(time));
            fl.central_u = None;
            let scale = scale?;
            fl.dev.sync()?;
            fl.timings.momentum_assembly += t.elapsed().as_secs_f64();
            let t = Instant::now();
            fl.momentum_solve(scale, &mut res)?;
            fl.dev.sync()?;
            fl.timings.momentum_solve += t.elapsed().as_secs_f64();

            // predicted fluxes with the time-step term
            let t = Instant::now();
            let a = fl.flux_args(Some(w.d), relax, Some((w.ds, l.y, alpha)));
            fl.dev.mass_flux(&a, w.m)?;
            fl.close_backflow(w.m)?;
            res.continuity = fl.continuity(w.m)?;
            fl.dev.sync()?;
            fl.timings.pressure_assembly += t.elapsed().as_secs_f64();

            for j in 1..=s.correctors {
                let rebuild = rebuild_step && k == 1 && j == 1;
                res.pressure_iterations += fl.pressure_correction(w.m, rebuild)?;
                let t = Instant::now();
                fl.corrections(w.m, alpha_p)?;
                if j < s.correctors {
                    // u' = δ − D∇p' of this corrector (δ = 0 in the first)
                    if j == 1 {
                        for b in l.ucorr.all() {
                            fl.dev.fill(b, 0.0)?;
                        }
                    }
                    fl.dev.correct_velocity(w.dc, w.gpc, l.ucorr)?;
                    // δ = H(u')/a_P (Issa) or Σ_N |a_N| (u'_N − u'_P)/a_P
                    // (SIMPLEC type), added to u and, as a flux, to ṁ
                    let own = if s.simplec_corrector() { 1.0 } else { 0.0 };
                    for i in 0..3 {
                        let (x, dg, dl) = (l.ucorr.all()[i], w.diag.all()[i], l.delta.all()[i]);
                        fl.dev.neighbour_correction(w.val, x, dg, dl, own)?;
                    }
                    for (u, d) in f.u.all().into_iter().zip(l.delta.all()) {
                        fl.dev.lincomb(u, 1.0, u, 1.0, d)?;
                    }
                    fl.dev.delta_flux(l.delta, l.dflux)?;
                    fl.dev.lincomb(w.m, 1.0, w.mc, 1.0, l.dflux)?;
                    for (d, u) in l.delta.all().into_iter().zip(l.ucorr.all()) {
                        fl.dev.copy(d, u)?;
                    }
                } else {
                    std::mem::swap(&mut fl.f.mass_flux, &mut fl.w.mc);
                }
                fl.dev.sync()?;
                fl.timings.pressure_assembly += t.elapsed().as_secs_f64();
            }
            if res.finite() {
                if let Some(t) = fl.turbulence.as_mut() {
                    t.amg_rebuild = rebuild_step && k == 1;
                }
                fl.solve_scalars(scalar_step, &mut res)?;
            }
            if !res.finite() {
                fl.diverged = Some(res.clone());
                return Err(crate::diag::err(crate::diag::Diagnostic::new(
                    crate::diag::RUN_DIVERGED,
                    format!(
                        "the solution diverged in time step {} (t = {:.6e} s, Δt = {dt:.3e} s; \
                         residuals {:?}, continuity {:.3e})",
                        self.clock.step + 1,
                        self.clock.time + dt,
                        res.momentum,
                        res.continuity
                    ),
                )));
            }
            if k == 1 {
                first = res.clone();
            } else {
                first.pressure_iterations += res.pressure_iterations;
            }
            last = res.clone();
            self.outer_log.push(res.clone());
            let flow_done = k >= s.outer_iterations || (k > 1 && res.max() < s.outer_tolerance);
            let more = match hook.as_deref_mut() {
                Some(h) => h.after_solve(fl, k)?,
                None => false,
            };
            if flow_done && !more {
                break;
            }
        }
        // the mass flux buffer may have moved: f is stale from here
        let t = Instant::now();
        fl.dev.copy(self.l.x_cur, self.l.x_old)?;
        self.finish_level(fl)?;
        let cfl = self.courant(fl, dt)?;
        // the local error of the corrector from its difference to the
        // predictor (Milne's device): with the step constants of uⁿ⁺¹ − u(tⁿ⁺¹)
        // = K h³u'''/6 and P − u(tⁿ⁺¹) = −Π h³u'''/6, the error is
        // K/(K + Π)·(uⁿ⁺¹ − P)
        let error = match (control, self.e) {
            (Some(c), Some(e)) if predict => {
                let (h, h1, h2) = (dt, self.clock.dt_prev, self.clock.dt_prev2);
                let k = match self.scheme {
                    TimeScheme::CrankNicolson => -0.5,
                    _ => (1.0 + omega) * (1.0 + omega) / (omega * (1.0 + 2.0 * omega)),
                };
                let pi = (h + h1) * (h + h1 + h2) / (h * h);
                let u = fl.f.u;
                let (mut dd, mut uu) = (0.0, 0.0);
                for ((d, a), p) in l.delta.all().into_iter().zip(u.all()).zip(e.pred.all()) {
                    fl.dev.lincomb(d, 1.0, a, -1.0, p)?;
                    dd += fl.dev.sum_sq(d)?;
                    uu += fl.dev.sum_sq(a)?;
                }
                let n = 3.0 * fl.disc.n_cells as f64;
                let rms = |x: f64| (x / n).sqrt();
                Some((k / (k + pi)).abs() * rms(dd) / (c.absolute + c.relative * rms(uu)))
            }
            _ => None,
        };
        fl.dev.sync()?;
        fl.timings.gradients += t.elapsed().as_secs_f64();
        self.clock.time += dt;
        self.clock.step += 1;
        self.clock.dt_prev2 = self.clock.dt_prev;
        self.clock.dt_prev = dt;
        self.clock.levels = (self.clock.levels + 1).min(2);
        Ok(StepReport {
            step: self.clock.step,
            time: self.clock.time,
            dt,
            residuals: first,
            last,
            outer,
            cfl,
            error,
        })
    }

    /// One step of size `dt` of the energy equation alone, with the flow
    /// held (the frozen-flow march, D118): velocity, pressure, mass flux
    /// and turbulence keep their values and their time levels, the
    /// temperatures of fluid and solids take the step with the held flux
    /// (one matrix, D158). `restart`: BDF1, for the first step after the
    /// step size has changed by more than BDF2 takes (from the flow's step
    /// to the march's); later steps are BDF2 as the case asks. The energy
    /// equation is solved up to `iterations` times, until the largest
    /// change of a temperature between two solves is below 10⁻⁸ of the
    /// temperature span (its properties follow T). The flow's levels are
    /// those of its last coupled step, and the next coupled step continues
    /// from them: its step ratio to the long step of the march is small,
    /// which makes it BDF1-like, and the flow had settled.
    pub fn energy_step(
        &mut self,
        fl: &mut Flow,
        dt: f64,
        restart: bool,
        iterations: u32,
    ) -> Result<StepReport> {
        ensure!(dt > 0.0 && dt.is_finite(), "time step {dt} must be > 0");
        ensure!(
            fl.cht.is_some(),
            "the energy-only step of the frozen flow needs a conjugate case"
        );
        let second = !restart && self.scheme != TimeScheme::Bdf1 && self.clock.levels >= 1;
        let omega = if second { dt / self.clock.dt_prev } else { 1.0 };
        let t = Instant::now();
        let mut prev = match crate::cht::energy_fields(fl) {
            Some((sc, dev)) => {
                sc.step_start(dev, self.clock.levels >= 1, second.then_some(omega))?;
                sc.energy_temperature(dev)?
            }
            None => Vec::new(),
        };
        let step = crate::scalar::StepTime {
            time: self.clock.time + dt,
            dt,
            a: crate::scalar::bdf(second, omega),
        };
        let mut first = Residuals::default();
        let mut last = Residuals::default();
        let mut outer = 0;
        self.outer_log.clear();
        for k in 1..=iterations.max(1) as usize {
            outer = k;
            fl.iteration += 1;
            let mut res = Residuals {
                iteration: fl.iteration,
                ..Default::default()
            };
            crate::cht::solve(fl, Some(step), &mut res)?;
            if !res.finite() {
                fl.diverged = Some(res.clone());
                return Err(crate::diag::err(crate::diag::Diagnostic::new(
                    crate::diag::RUN_DIVERGED,
                    format!(
                        "the energy equation diverged in the frozen-flow step {} (t = {:.6e} s, \
                         Δt = {dt:.3e} s; residuals {:?})",
                        self.clock.step + 1,
                        self.clock.time + dt,
                        res.scalars.residual
                    ),
                )));
            }
            let now = match crate::cht::energy_fields(fl) {
                Some((sc, dev)) => sc.energy_temperature(dev)?,
                None => Vec::new(),
            };
            let (lo, hi) = now
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| {
                    (a.min(x), b.max(x))
                });
            let change = now
                .iter()
                .zip(&prev)
                .fold(0.0_f64, |a, (x, y)| a.max((x - y).abs()));
            prev = now;
            if k == 1 {
                first = res.clone();
            }
            last = res.clone();
            self.outer_log.push(res);
            if k > 1 && change <= 1e-8 * (hi - lo).max(1.0) {
                break;
            }
        }
        fl.dev.sync()?;
        fl.timings.scalars += t.elapsed().as_secs_f64();
        let cfl = self.courant(fl, dt)?;
        self.clock.time += dt;
        self.clock.step += 1;
        self.clock.dt_prev2 = self.clock.dt_prev;
        self.clock.dt_prev = dt;
        self.clock.levels = (self.clock.levels + 1).min(2);
        Ok(StepReport {
            step: self.clock.step,
            time: self.clock.time,
            dt,
            residuals: first,
            last,
            outer,
            cfl,
            error: None,
        })
    }
}
