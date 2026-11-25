//! Structural dynamics: Rayleigh damping, the generalised-α method of
//! Chung & Hulbert (1993), mode superposition with modal damping ratios,
//! and the Craig–Bampton reduction (1968) of a substructure.
//!
//! Everything works on the free unknowns of an assembled [`System`]; the
//! fixed unknowns stay at zero (prescribed motion in time comes with the
//! coupling of MS6).
//!
//! **Generalised-α.** With αm = (2ρ∞ − 1)/(ρ∞ + 1), αf = ρ∞/(ρ∞ + 1),
//! γ = ½ − αm + αf and β = ¼(1 − αm + αf)², the balance
//! M a(n+1−αm) + C v(n+1−αf) + K u(n+1−αf) = f(n+1−αf), with
//! x(n+1−α) = (1 − α) x(n+1) + α x(n), and the Newmark updates give one
//! linear system for a(n+1) whose matrix
//! S = (1 − αm) M + (1 − αf) γ Δt C + (1 − αf) β Δt² K
//! is factorised once. The method is second order, unconditionally
//! stable, and damps the frequencies far above 1/Δt to the spectral radius
//! ρ∞ (1: no numerical damping, the trapezoidal rule, which conserves the
//! energy of an undamped linear system exactly).
//!
//! **Mode superposition.** Each modal equation
//! q'' + 2ζω q' + ω² q = φᵀ f(t) is integrated exactly for a force linear
//! within each step (Nigam & Jennings 1969), so the only error is the
//! truncation to the modes kept.

use crate::analysis::{self, PIVOT_TOL, dof_name};
use crate::eigen;
use crate::ldlt::Ldlt;
use crate::model::{Model, System};
use crate::sparse::{Assembly, SymCsc, dot};
use crate::{Error, codes};

/// Rayleigh damping C = α M + β K.
#[derive(Debug, Clone, Copy, Default, PartialEq, serde::Serialize)]
pub struct Rayleigh {
    /// Mass factor α (1/s).
    pub alpha: f64,
    /// Stiffness factor β (s).
    pub beta: f64,
}

impl Rayleigh {
    /// The damping ratio at the circular frequency ω: α/(2ω) + βω/2.
    pub fn ratio(&self, omega: f64) -> f64 {
        self.alpha / (2.0 * omega) + self.beta * omega / 2.0
    }

    /// α and β that give the ratio ζ1 at f1 and ζ2 at f2 (Hz, f1 ≠ f2).
    pub fn from_ratios(f1: f64, z1: f64, f2: f64, z2: f64) -> Result<Self, Error> {
        let bad = |m: String| Error::new(codes::CASE, m);
        if !(f1 > 0.0 && f2 > 0.0) || (f1 - f2).abs() <= 1e-12 * f1.max(f2) {
            return Err(bad(format!(
                "Rayleigh damping needs two different positive frequencies (got {f1} and {f2} Hz)"
            )));
        }
        if !(z1 >= 0.0 && z2 >= 0.0) {
            return Err(bad("damping ratios must not be negative".into()));
        }
        let (w1, w2) = (
            2.0 * std::f64::consts::PI * f1,
            2.0 * std::f64::consts::PI * f2,
        );
        // 2ζω = α + βω² at both frequencies
        let alpha = 2.0 * w1 * w2 * (z1 * w2 - z2 * w1) / (w2 * w2 - w1 * w1);
        let beta = 2.0 * (z2 * w2 - z1 * w1) / (w2 * w2 - w1 * w1);
        let r = Rayleigh { alpha, beta };
        if alpha < -1e-12 * alpha.abs().max(1.0) || beta < 0.0 {
            return Err(bad(format!(
                "the ratios {z1} at {f1} Hz and {z2} at {f2} Hz give a negative Rayleigh factor (α = {alpha:.4e}, β = {beta:.4e}); choose ratios that a C = αM + βK can give"
            )));
        }
        Ok(r)
    }

    /// C x.
    fn apply(&self, k: &SymCsc, m: &SymCsc, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; x.len()];
        if self.alpha != 0.0 {
            for (o, v) in out.iter_mut().zip(m.apply(x)) {
                *o += self.alpha * v;
            }
        }
        if self.beta != 0.0 {
            for (o, v) in out.iter_mut().zip(k.apply(x)) {
                *o += self.beta * v;
            }
        }
        out
    }
}

/// Displacement, velocity and acceleration on the free unknowns at a time.
#[derive(Debug, Clone)]
pub struct State {
    /// Time (s).
    pub t: f64,
    /// Displacements.
    pub u: Vec<f64>,
    /// Velocities.
    pub v: Vec<f64>,
    /// Accelerations.
    pub a: Vec<f64>,
}

impl State {
    /// Kinetic energy ½ vᵀ M v and strain energy ½ uᵀ K u.
    pub fn energies(&self, sys: &System) -> (f64, f64) {
        (
            0.5 * sys.m.form(&self.v, &self.v),
            0.5 * sys.k.form(&self.u, &self.u),
        )
    }
}

/// The generalised-α integrator of one system and one time step.
pub struct GenAlpha<'a> {
    sys: &'a System,
    damping: Rayleigh,
    dt: f64,
    am: f64,
    af: f64,
    gamma: f64,
    beta: f64,
    /// S / ck = K + (cm/ck) M, and its factor.
    base: SymCsc,
    s: Ldlt,
    s_scale: f64,
}

impl<'a> GenAlpha<'a> {
    /// The integrator for the step `dt` and the spectral radius `rho_inf`
    /// at infinite frequency (0 to 1).
    pub fn new(sys: &'a System, damping: Rayleigh, dt: f64, rho_inf: f64) -> Result<Self, Error> {
        if !(dt > 0.0 && dt.is_finite()) {
            return Err(Error::new(
                codes::CASE,
                format!("the time step must be positive (got {dt})"),
            ));
        }
        if !(0.0..=1.0).contains(&rho_inf) {
            return Err(Error::new(
                codes::CASE,
                format!("`rho_infinity` must lie between 0 and 1 (got {rho_inf})"),
            ));
        }
        let am = (2.0 * rho_inf - 1.0) / (rho_inf + 1.0);
        let af = rho_inf / (rho_inf + 1.0);
        let gamma = 0.5 - am + af;
        let beta = 0.25 * (1.0 - am + af).powi(2);
        // S = cm M + ck K
        let cm = (1.0 - am) + (1.0 - af) * gamma * dt * damping.alpha;
        let ck = (1.0 - af) * gamma * dt * damping.beta + (1.0 - af) * beta * dt * dt;
        let base = sys.k.plus_scaled(cm / ck, &sys.m);
        let f = Ldlt::factor(&base, PIVOT_TOL).map_err(|z| {
            Error::new(
                codes::SINGULAR,
                format!(
                    "the matrix of the time step is singular at {} (pivot {:.3e})",
                    sys.free
                        .get(z.index)
                        .map_or("?".into(), |&d| format!("unknown {d}")),
                    z.pivot
                ),
            )
        })?;
        Ok(GenAlpha {
            sys,
            damping,
            dt,
            am,
            af,
            gamma,
            beta,
            base,
            s: f,
            s_scale: ck,
        })
    }

    /// The time step.
    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// The state at `t0` from displacements, velocities and the load: the
    /// acceleration from M a = f − C v − K u. A unit of the model without
    /// mass (the drilling rotation of a flat shell) gets the acceleration
    /// of a mass 10⁻¹² of the largest.
    pub fn initial(&self, t0: f64, u: Vec<f64>, v: Vec<f64>, f: &[f64]) -> Result<State, Error> {
        let sys = self.sys;
        let n = sys.free.len();
        let ku = sys.k.apply(&u);
        let cv = self.damping.apply(&sys.k, &sys.m, &v);
        let r: Vec<f64> = (0..n).map(|i| f[i] - ku[i] - cv[i]).collect();
        let a = if r.iter().all(|x| *x == 0.0) {
            vec![0.0; n]
        } else {
            let mut asm = Assembly::new(n);
            let dm = sys.m.diagonal();
            let eps = 1e-12 * dm.iter().fold(0.0f64, |a, b| a.max(b.abs()));
            for (i, &d) in dm.iter().enumerate() {
                if d.abs() <= eps {
                    asm.add(i, i, eps.max(f64::MIN_POSITIVE));
                }
            }
            let mm = sys.m.plus_scaled(1.0, &asm.finish());
            let fm = Ldlt::factor(&mm, PIVOT_TOL).map_err(|z| {
                Error::new(
                    codes::SINGULAR,
                    format!(
                        "the mass matrix is singular at unknown {} (pivot {:.3e})",
                        z.index, z.pivot
                    ),
                )
            })?;
            fm.solve(&r)
        };
        Ok(State { t: t0, u, v, a })
    }

    /// The predictors of a step: ũ, ṽ, the α-point displacement at a = 0
    /// and the right-hand side of (cm M + ck K) a = rhs.
    fn predict(
        &self,
        s: &State,
        f_n: &[f64],
        f_next: &[f64],
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let sys = self.sys;
        let (dt, am, af, g, b) = (self.dt, self.am, self.af, self.gamma, self.beta);
        let n = s.u.len();
        let ut: Vec<f64> = (0..n)
            .map(|i| s.u[i] + dt * s.v[i] + dt * dt * (0.5 - b) * s.a[i])
            .collect();
        let vt: Vec<f64> = (0..n).map(|i| s.v[i] + dt * (1.0 - g) * s.a[i]).collect();
        let ua: Vec<f64> = (0..n).map(|i| (1.0 - af) * ut[i] + af * s.u[i]).collect();
        let va: Vec<f64> = (0..n).map(|i| (1.0 - af) * vt[i] + af * s.v[i]).collect();
        let ma = sys.m.apply(&s.a);
        let ku = sys.k.apply(&ua);
        let cv = self.damping.apply(&sys.k, &sys.m, &va);
        let rhs: Vec<f64> = (0..n)
            .map(|i| (1.0 - af) * f_next[i] + af * f_n[i] - am * ma[i] - cv[i] - ku[i])
            .collect();
        (ut, vt, ua, rhs)
    }

    /// One step with unilateral contact (penalty at the α-point, an active
    /// set loop in the step; see [`crate::contact`]).
    pub fn step_contact(
        &self,
        model: &crate::model::Model,
        c: &crate::contact::Contacts,
        cache: &mut crate::contact::StepCache,
        s: &State,
        f_n: &[f64],
        f_next: &[f64],
    ) -> Result<(State, crate::contact::State), Error> {
        let (dt, af, g, b) = (self.dt, self.af, self.gamma, self.beta);
        let (ut, vt, ua, rhs) = self.predict(s, f_n, f_next);
        let cu = (1.0 - af) * b * dt * dt;
        let (a, mut cs) = crate::contact::step(
            model,
            self.sys,
            c,
            &self.base,
            self.s_scale,
            cu,
            &ua,
            &rhs,
            cache,
        )?;
        let n = a.len();
        let u: Vec<f64> = (0..n).map(|i| ut[i] + b * dt * dt * a[i]).collect();
        let v = (0..n).map(|i| vt[i] + g * dt * a[i]).collect();
        // the force acts at the α-point; the gap reported is the one at the
        // end of the step, where the displacement is written out
        cs.gap = (0..c.pairs.len()).map(|p| c.gap(p, &u)).collect();
        Ok((
            State {
                t: s.t + dt,
                u,
                v,
                a,
            },
            cs,
        ))
    }

    /// One step from `s` (load `f_n`) to `s.t + dt` (load `f_next`).
    pub fn step(&self, s: &State, f_n: &[f64], f_next: &[f64]) -> State {
        let sys = self.sys;
        let (dt, am, af, g, b) = (self.dt, self.am, self.af, self.gamma, self.beta);
        let n = s.u.len();
        let ut: Vec<f64> = (0..n)
            .map(|i| s.u[i] + dt * s.v[i] + dt * dt * (0.5 - b) * s.a[i])
            .collect();
        let vt: Vec<f64> = (0..n).map(|i| s.v[i] + dt * (1.0 - g) * s.a[i]).collect();
        // right-hand side
        let ua: Vec<f64> = (0..n).map(|i| (1.0 - af) * ut[i] + af * s.u[i]).collect();
        let va: Vec<f64> = (0..n).map(|i| (1.0 - af) * vt[i] + af * s.v[i]).collect();
        let ma = sys.m.apply(&s.a);
        let ku = sys.k.apply(&ua);
        let cv = self.damping.apply(&sys.k, &sys.m, &va);
        let rhs: Vec<f64> = (0..n)
            .map(|i| (1.0 - af) * f_next[i] + af * f_n[i] - am * ma[i] - cv[i] - ku[i])
            .collect();
        let a: Vec<f64> = self
            .s
            .solve(&rhs)
            .into_iter()
            .map(|x| x / self.s_scale)
            .collect();
        let u = (0..n).map(|i| ut[i] + b * dt * dt * a[i]).collect();
        let v = (0..n).map(|i| vt[i] + g * dt * a[i]).collect();
        State {
            t: s.t + dt,
            u,
            v,
            a,
        }
    }
}

/// A reduced model of the lowest modes, for mode superposition.
#[derive(Debug, Clone)]
pub struct ModalBasis {
    /// Circular frequencies ω (rad/s).
    pub omega: Vec<f64>,
    /// Damping ratio of each mode.
    pub zeta: Vec<f64>,
    /// M-normalised shapes on the free unknowns.
    pub shapes: Vec<Vec<f64>>,
}

impl ModalBasis {
    /// The lowest `count` modes of the system with the ratio `zeta` each,
    /// plus the Rayleigh ratio of `rayleigh` at each mode's frequency.
    pub fn lowest(
        model: &Model,
        count: usize,
        zeta: f64,
        rayleigh: Rayleigh,
    ) -> Result<Self, Error> {
        let md = analysis::solve_modal(model, count, 0.0)?;
        let sys = model.assemble();
        let mut omega = Vec::new();
        let mut zs = Vec::new();
        let mut shapes = Vec::new();
        for m in &md.modes {
            let w = m.lambda.max(0.0).sqrt();
            if w <= 0.0 {
                return Err(Error::new(
                    codes::CASE,
                    "mode superposition needs a supported structure (a rigid-body mode at 0 Hz was found)".to_string(),
                ));
            }
            let z = zeta
                + if rayleigh == Rayleigh::default() {
                    0.0
                } else {
                    rayleigh.ratio(w)
                };
            if z >= 1.0 {
                return Err(Error::new(
                    codes::CASE,
                    format!(
                        "mode at {:.4} Hz has a damping ratio {z:.3} ≥ 1: mode superposition takes underdamped modes",
                        m.frequency()
                    ),
                ));
            }
            omega.push(w);
            zs.push(z);
            shapes.push(sys.restrict(&m.shape));
        }
        Ok(ModalBasis {
            omega,
            zeta: zs,
            shapes,
        })
    }

    /// The modal forces φᵀ f.
    pub fn project(&self, f: &[f64]) -> Vec<f64> {
        self.shapes.iter().map(|p| dot(p, f)).collect()
    }

    /// The free-unknown vector Σ q_i φ_i.
    pub fn expand(&self, q: &[f64]) -> Vec<f64> {
        let n = self.shapes.first().map_or(0, |s| s.len());
        let mut out = vec![0.0; n];
        for (qi, p) in q.iter().zip(&self.shapes) {
            for (o, x) in out.iter_mut().zip(p) {
                *o += qi * x;
            }
        }
        out
    }

    /// One exact step of every modal equation for a modal force linear
    /// from `p0` to `p1` over `dt`: (q, q') to their values at the end.
    pub fn step(&self, q: &mut [f64], qd: &mut [f64], p0: &[f64], p1: &[f64], dt: f64) {
        for i in 0..q.len() {
            let (w, z) = (self.omega[i], self.zeta[i]);
            let wd = w * (1.0 - z * z).sqrt();
            let a = p0[i];
            let b = (p1[i] - p0[i]) / dt;
            // particular solution for a + bτ
            let qp0 = a / (w * w) - 2.0 * z * b / (w * w * w);
            let qpd = b / (w * w);
            let e = q[i] - qp0;
            let ed = qd[i] - qpd;
            let decay = (-z * w * dt).exp();
            let (s, c) = (wd * dt).sin_cos();
            let c2 = (ed + z * w * e) / wd;
            let x = decay * (e * c + c2 * s);
            let xd = decay * ((-z * w) * (e * c + c2 * s) + (-e * wd * s + c2 * wd * c));
            q[i] = x + qp0 + qpd * dt;
            qd[i] = xd + qpd;
        }
    }
}

/// A Craig–Bampton reduced model: the boundary unknowns kept, the interior
/// represented by its constraint modes and its lowest fixed-interface
/// modes.
#[derive(Debug, Clone)]
pub struct Reduced {
    /// Model unknowns kept (the boundary), in order.
    pub boundary: Vec<usize>,
    /// Number of fixed-interface modes kept.
    pub modes: usize,
    /// Circular frequencies of the fixed-interface modes (rad/s).
    pub interface_omega: Vec<f64>,
    /// Reduced stiffness, dense row-major, size (boundary + modes)².
    pub k: Vec<f64>,
    /// Reduced mass, dense row-major.
    pub m: Vec<f64>,
    /// Columns of the transformation on the free unknowns of the model:
    /// the constraint modes, then the fixed-interface modes.
    pub basis: Vec<Vec<f64>>,
}

impl Reduced {
    /// Size of the reduced model.
    pub fn n(&self) -> usize {
        self.boundary.len() + self.modes
    }
}

/// The Craig–Bampton reduction (1968) of the free unknowns of `model` to
/// the model unknowns `boundary` and `modes` fixed-interface modes.
pub fn craig_bampton(model: &Model, boundary: &[usize], modes: usize) -> Result<Reduced, Error> {
    let sys = model.assemble();
    let nf = sys.free.len();
    let bad = |m: String| Error::new(codes::CASE, m);
    let mut is_b = vec![false; nf];
    let mut bidx = Vec::new();
    for &d in boundary {
        let i = sys
            .map
            .get(d)
            .copied()
            .flatten()
            .ok_or_else(|| bad(format!("boundary unknown {d} is fixed or not in the model")))?;
        if is_b[i] {
            return Err(bad(format!("boundary unknown {d} is listed twice")));
        }
        is_b[i] = true;
        bidx.push(i);
    }
    let interior: Vec<usize> = (0..nf).filter(|&i| !is_b[i]).collect();
    let ni = interior.len();
    if modes > ni {
        return Err(bad(format!(
            "only {ni} interior unknowns: fewer modes than {modes}"
        )));
    }
    let mut pos = vec![usize::MAX; nf];
    for (k, &i) in interior.iter().enumerate() {
        pos[i] = k;
    }
    let sub = |a: &SymCsc| {
        let mut asm = Assembly::new(ni);
        for j in 0..a.n {
            for p in a.colptr[j]..a.colptr[j + 1] {
                let i = a.rowidx[p];
                if pos[i] != usize::MAX && pos[j] != usize::MAX {
                    asm.add(pos[i].min(pos[j]), pos[i].max(pos[j]), a.vals[p]);
                }
            }
        }
        asm.finish()
    };
    let (kii, mii) = (sub(&sys.k), sub(&sys.m));
    let fk = Ldlt::factor(&kii, PIVOT_TOL).map_err(|z| {
        Error::new(
            codes::SINGULAR,
            format!(
                "the interior is free to move at {} with the boundary held: add a support or a boundary unknown",
                dof_name(model, sys.free[interior[z.index]])
            ),
        )
    })?;
    let mut basis = Vec::new();
    // constraint modes: unit boundary motion, interior in static balance
    for &bi in &bidx {
        let mut e = vec![0.0; nf];
        e[bi] = 1.0;
        let ke = sys.k.apply(&e);
        let r: Vec<f64> = interior.iter().map(|&i| -ke[i]).collect();
        let x = fk.solve(&r);
        for (k, &i) in interior.iter().enumerate() {
            e[i] = x[k];
        }
        basis.push(e);
    }
    let mut interface_omega = Vec::new();
    if modes > 0 {
        let pairs = eigen::lanczos(&kii, &mii, modes, 0.0, 1e-10).map_err(|_| {
            Error::new(
                codes::MODES,
                "the fixed-interface modes did not converge".into(),
            )
        })?;
        for (val, vec) in pairs.values.iter().zip(&pairs.vectors).take(modes) {
            interface_omega.push(val.max(0.0).sqrt());
            let mut e = vec![0.0; nf];
            for (k, &i) in interior.iter().enumerate() {
                e[i] = vec[k];
            }
            basis.push(e);
        }
    }
    let n = basis.len();
    let kb: Vec<Vec<f64>> = basis.iter().map(|b| sys.k.apply(b)).collect();
    let mb: Vec<Vec<f64>> = basis.iter().map(|b| sys.m.apply(b)).collect();
    let mut k = vec![0.0; n * n];
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            k[i * n + j] = dot(&basis[i], &kb[j]);
            m[i * n + j] = dot(&basis[i], &mb[j]);
        }
    }
    // symmetrise the rounding
    for i in 0..n {
        for j in 0..i {
            let (a, b) = (
                0.5 * (k[i * n + j] + k[j * n + i]),
                0.5 * (m[i * n + j] + m[j * n + i]),
            );
            k[i * n + j] = a;
            k[j * n + i] = a;
            m[i * n + j] = b;
            m[j * n + i] = b;
        }
    }
    Ok(Reduced {
        boundary: boundary.to_vec(),
        modes: interface_omega.len(),
        interface_omega,
        k,
        m,
        basis,
    })
}

/// Eigenvalues of the dense symmetric pencil (K, M), M positive definite,
/// ascending: Cholesky of M and cyclic Jacobi on L⁻¹ K L⁻ᵀ.
pub fn dense_eigenvalues(k: &[f64], m: &[f64], n: usize) -> Option<Vec<f64>> {
    // Cholesky M = L Lᵀ
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = m[i * n + j];
            for p in 0..j {
                s -= l[i * n + p] * l[j * n + p];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i * n + i] = s.sqrt();
            } else {
                l[i * n + j] = s / l[j * n + j];
            }
        }
    }
    // A = L⁻¹ K L⁻ᵀ: solve L X = K, then L Aᵀ = Xᵀ
    let mut x = k.to_vec();
    for c in 0..n {
        for i in 0..n {
            let mut s = x[i * n + c];
            for p in 0..i {
                s -= l[i * n + p] * x[p * n + c];
            }
            x[i * n + c] = s / l[i * n + i];
        }
    }
    let mut a = vec![0.0; n * n];
    for r in 0..n {
        for i in 0..n {
            let mut s = x[r * n + i];
            for p in 0..i {
                s -= l[i * n + p] * a[r * n + p];
            }
            a[r * n + i] = s / l[i * n + i];
        }
    }
    for _sweep in 0..100 {
        let mut off = 0.0;
        let mut diag = 0.0;
        for i in 0..n {
            diag += a[i * n + i] * a[i * n + i];
            for j in 0..i {
                off += a[i * n + j] * a[i * n + j];
            }
        }
        if off <= 1e-30 * diag.max(f64::MIN_POSITIVE) {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq.abs() <= 1e-300 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for r in 0..n {
                    let (arp, arq) = (a[r * n + p], a[r * n + q]);
                    a[r * n + p] = c * arp - s * arq;
                    a[r * n + q] = s * arp + c * arq;
                }
                for r in 0..n {
                    let (apr, aqr) = (a[p * n + r], a[q * n + r]);
                    a[p * n + r] = c * apr - s * aqr;
                    a[q * n + r] = s * apr + c * aqr;
                }
            }
        }
    }
    let mut ev: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    ev.sort_by(f64::total_cmp);
    Some(ev)
}
