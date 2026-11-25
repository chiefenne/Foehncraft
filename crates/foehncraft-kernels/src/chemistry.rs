//! Gas-phase chemistry of an ideal-gas mixture (steam furnace stage SF2,
//! the gas-phase slice of the chemistry program D120): the rates of a
//! reaction mechanism, their analytic Jacobian, and a stiff integrator of
//! the constant-pressure adiabatic reactor of one cell.
//!
//! Every function is cell-local, allocation-free (the caller hands in the
//! work array of [`work_len`]) and runs its loops in the fixed order of the
//! species and reactions, so a CPU core and a GPU thread give the same bits
//! (D31) and a rank integrates its cells as one rank does (D93). The
//! transcendental functions are the in-house ones of
//! [`foehncraft_core::math`] (D93 (6)).
//!
//! # Rate laws
//!
//! Reaction j has the rate of progress (mol/(m³ s))
//!
//! q_j = k_f,j Π_k φ(C_k, a'_kj) − k_r,j Π_k φ(C_k, a''_kj),
//!
//! with the concentrations C_k = ρ Y_k/W_k and the orders a' (forward) and
//! a'' (reverse). An elementary reaction has the stoichiometric
//! coefficients as orders; a global one has its own forward orders, and
//! when it is reversible the reverse orders a'' = a' − ν' + ν'', so that
//! the ratio of the two products is the equilibrium quotient Π C^(ν''−ν')
//! and the reaction stops at the chemical equilibrium.
//!
//! - k_f = A T^b exp(−T_a/T) (modified Arrhenius);
//! - three-body reactions multiply k_f by [M] = Σ α_k C_k (efficiencies α,
//!   the default for species not listed);
//! - fall-off reactions blend the low- and high-pressure limits:
//!   k_f = k_∞ P_r/(1 + P_r) F, P_r = k₀[M]/k_∞, with F = 1 (Lindemann) or
//!   the broadening of Troe (Gilbert, Luther & Troe 1983): log F =
//!   log F_c/(1 + ((log P_r + c)/(n − 0.14 (log P_r + c)))²), c = −0.4 −
//!   0.67 log F_c, n = 0.75 − 1.27 log F_c, F_c = (1 − a) exp(−T/T₃) +
//!   a exp(−T/T₁) + exp(−T₂/T);
//! - reversible reactions take k_r = k_f/K_c with K_c = K_p (p°/(R T))^Δν
//!   and ln K_p = −Σ ν_k g_k°/(R T) from the NASA-7 polynomials (Kee,
//!   Coltrin & Glarborg 2003, chapter 9).
//!
//! An integer order n ≥ 1 is the power Cⁿ (a slightly negative C of a
//! round-off then gives a rate of the right sign). A fractional order
//! above 1 takes max(C, 0)^a. An order below 1 (fractional or negative,
//! as in global mechanisms) would have an infinite derivative or rate at
//! C = 0; it is regularised as φ = C (max(C, 0) + ε)^(a − 1), which is
//! C^a for C ≫ ε, goes to zero linearly at C = 0 and is smooth, with
//! ε = 10⁻¹² of the total concentration (a mole fraction of 10⁻¹²).
//!
//! # Reactor
//!
//! The constant-pressure adiabatic reactor of a cell, with the state
//! z = (Y₁ … Y_N, T):
//!
//! dY_k/dt = W_k ω̇_k/ρ,  dT/dt = −Σ_k H_k ω̇_k/(ρ c_p),
//!
//! ω̇_k = Σ_j ν_kj q_j, H_k the molar enthalpy (formation included), c_p
//! the specific heat of the mixture. The Jacobian ∂f/∂z is analytic (D120,
//! D59), with the dependence of ρ and of the concentrations on Y and T.
//!
//! The integrator is the L-stable four-stage Rosenbrock method RODAS3 of
//! order 3 with an embedded method of order 2 (Sandu, Verwer, Blom,
//! Spee, Carmichael & Potra 1997, "Benchmarking stiff ODE solvers for
//! atmospheric chemistry problems II: Rosenbrock solvers", Atmos. Environ.
//! 31, 3459–3472), with step-size control on the embedded error. Each
//! step solves its linear systems with one LU factorisation of
//! I/(γh) − J by Gaussian elimination with partial pivoting in the fixed
//! column order 1 … N + 1: the pivot choice depends only on the bits of
//! the matrix, so the elimination is the same on every target (the fixed
//! elimination order of the SF2 plan). Rosenbrock stages are linear in the
//! right-hand side, so element masses are conserved to round-off.

use foehncraft_core::math;

use crate::species::{R, Species};

/// The standard pressure p° of the thermodynamic data (Pa).
pub const P_STANDARD: f64 = 101_325.0;

/// An elementary reaction: k_f = A T^b exp(−T_a/T).
pub const ELEMENTARY: u32 = 0;
/// A three-body reaction: k_f [M].
pub const THREE_BODY: u32 = 1;
/// A fall-off reaction with the Lindemann form (F = 1).
pub const LINDEMANN: u32 = 2;
/// A fall-off reaction with the Troe broadening.
pub const TROE: u32 = 3;

/// One reaction of a mechanism, flat (`repr(C)`, for the GPU as well).
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Reaction {
    /// [`ELEMENTARY`], [`THREE_BODY`], [`LINDEMANN`] or [`TROE`].
    pub kind: u32,
    /// 1: reversible with k_r = k_f/K_c; 0: forward only.
    pub reversible: u32,
    /// The range of the reaction's species terms in the term array.
    pub terms: [u32; 2],
    /// The range of its third-body efficiencies in the efficiency array.
    pub efficiencies: [u32; 2],
    /// The efficiency of the species not listed.
    pub default_efficiency: f64,
    /// A (SI units: m, mol, s), b, T_a = E_a/R (K); the high-pressure limit
    /// of a fall-off reaction.
    pub rate: [f64; 3],
    /// The low-pressure limit k₀ of a fall-off reaction.
    pub low: [f64; 3],
    /// Troe: a, T₃, T₁, T₂ (K); T₂ is used only when `has_t2` is 1.
    pub troe: [f64; 4],
    /// 1 when the Troe form has its fourth parameter T₂.
    pub has_t2: u32,
}

/// One species of a reaction.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Term {
    /// The species index.
    pub species: u32,
    /// Its stoichiometric coefficient among the reactants (ν').
    pub reactant: f64,
    /// Among the products (ν'').
    pub product: f64,
    /// The order of the forward rate (a').
    pub order: f64,
}

/// A third-body efficiency.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Efficiency {
    /// The species index.
    pub species: u32,
    /// Its efficiency α.
    pub value: f64,
}

/// A mechanism as the kernels read it.
#[derive(Debug, Clone, Copy)]
pub struct Program<'a> {
    /// The species (NASA-7 data and molar masses).
    pub species: &'a [Species],
    /// The reactions.
    pub reactions: &'a [Reaction],
    /// The species terms of all reactions.
    pub terms: &'a [Term],
    /// The third-body efficiencies of all reactions.
    pub efficiencies: &'a [Efficiency],
}

impl Program<'_> {
    /// The number of species.
    #[inline(always)]
    pub fn n(&self) -> usize {
        self.species.len()
    }
}

/// The factor φ(C, a) of a concentration in a rate and its derivative
/// dφ/dC (see the module documentation).
#[inline(always)]
pub fn phi(c: f64, a: f64, eps: f64) -> (f64, f64) {
    if a == 0.0 {
        (1.0, 0.0)
    } else if a == 1.0 {
        (c, 1.0)
    } else if a == 2.0 {
        (c * c, 2.0 * c)
    } else if a == 3.0 {
        (c * c * c, 3.0 * c * c)
    } else if a > 1.0 {
        let cp = c.max(0.0);
        let p = math::pow(cp, a - 1.0);
        (cp * p, a * p)
    } else {
        let s = c.max(0.0) + eps;
        let p = math::pow(s, a - 1.0);
        let d = if c > 0.0 {
            p + (a - 1.0) * c * p / s
        } else {
            p
        };
        (c * p, d)
    }
}

impl Species {
    /// s°/R at `t` (the entropy at p°), held at the limits of the data.
    #[inline(always)]
    pub fn s_over_r(&self, t: f64) -> f64 {
        let t = t.max(self.t_range[0]).min(self.t_range[2]);
        let a = if t <= self.t_range[1] {
            &self.low
        } else {
            &self.high
        };
        a[0] * math::ln(t)
            + t * (a[1] + t * (a[2] / 2.0 + t * (a[3] / 3.0 + t * a[4] / 4.0)))
            + a[6]
    }

    /// d(c_p/R)/dT at `t` (zero outside the data range, where c_p is held).
    #[inline(always)]
    pub fn dcp_over_r_dt(&self, t: f64) -> f64 {
        if t < self.t_range[0] || t > self.t_range[2] {
            return 0.0;
        }
        let a = if t <= self.t_range[1] {
            &self.low
        } else {
            &self.high
        };
        a[1] + t * (2.0 * a[2] + t * (3.0 * a[3] + t * 4.0 * a[4]))
    }
}

/// k = A T^b exp(−T_a/T) and d ln k/dT.
#[inline(always)]
fn arrhenius(r: &[f64; 3], t: f64, ln_t: f64) -> (f64, f64) {
    let k = r[0] * math::exp(r[1] * ln_t - r[2] / t);
    (k, (r[1] + r[2] / t) / t)
}

/// The forward rate coefficient of reaction `rx` with its derivatives
/// (k, ∂k/∂T, ∂k/∂[M]) at the third-body concentration `m` (the factor
/// [M] of a three-body reaction included).
#[inline(always)]
fn forward(rx: &Reaction, t: f64, ln_t: f64, m: f64) -> (f64, f64, f64) {
    let (k, dl) = arrhenius(&rx.rate, t, ln_t);
    match rx.kind {
        THREE_BODY => (k * m, k * m * dl, k),
        LINDEMANN | TROE => {
            let (k0, dl0) = arrhenius(&rx.low, t, ln_t);
            if !(k > 0.0) || !(k0 > 0.0) || !(m > 0.0) {
                return (0.0, 0.0, if k > 0.0 { k0 } else { 0.0 });
            }
            let pr = k0 * m / k;
            // ln k_f = ln k + ln P_r − ln(1 + P_r) + ln F
            let dlpr = 1.0 / (1.0 + pr); // d ln k_f/d ln P_r without F
            let (lnf, dlnf_dlpr, dlnf_dt) = if rx.kind == TROE {
                troe(&rx.troe, rx.has_t2 == 1, t, pr)
            } else {
                (0.0, 0.0, 0.0)
            };
            let kf = k * pr / (1.0 + pr) * math::exp(lnf);
            let g = dlpr + dlnf_dlpr; // d ln k_f/d ln P_r
            // d ln P_r/dT = dl0 − dl; d ln P_r/d[M] = 1/[M]
            let dt = kf * (dl + g * (dl0 - dl) + dlnf_dt);
            let dm = kf * g / m;
            (kf, dt, dm)
        }
        _ => (k, k * dl, 0.0),
    }
}

/// ln F of the Troe form with d ln F/d ln P_r and ∂ ln F/∂T at fixed P_r.
#[inline(always)]
fn troe(p: &[f64; 4], has_t2: bool, t: f64, pr: f64) -> (f64, f64, f64) {
    let (a, t3, t1, t2) = (p[0], p[1], p[2], p[3]);
    let mut fc = 0.0;
    let mut dfc = 0.0;
    if t3 != 0.0 {
        let e = math::exp(-t / t3);
        fc += (1.0 - a) * e;
        dfc -= (1.0 - a) * e / t3;
    }
    if t1 != 0.0 {
        let e = math::exp(-t / t1);
        fc += a * e;
        dfc -= a * e / t1;
    }
    if has_t2 {
        let e = math::exp(-t2 / t);
        fc += e;
        dfc += e * t2 / (t * t);
    }
    if !(fc > 0.0) {
        return (0.0, 0.0, 0.0);
    }
    let ln10 = core::f64::consts::LN_10;
    let l = math::log10(fc);
    let lp = math::log10(pr.max(1e-300));
    let c = -0.4 - 0.67 * l;
    let n = 0.75 - 1.27 * l;
    let u = lp + c;
    let den = n - 0.14 * u;
    let x = u / den;
    let q = 1.0 + x * x;
    let log_f = l / q;
    // dX/d lp = n/den²; dX/dL through c and n
    let dx_dlp = n / (den * den);
    let dx_dl = (-0.67 * den - u * (-1.27 + 0.14 * 0.67)) / (den * den);
    let dlogf_dx = -l * 2.0 * x / (q * q);
    let dlogf_dlp = dlogf_dx * dx_dlp;
    let dlogf_dl = 1.0 / q + dlogf_dx * dx_dl;
    // d lp/d ln P_r = 1/ln 10; dL/dT = F_c'/(F_c ln 10); ln F = ln 10 log F
    let dlnf_dlnpr = dlogf_dlp;
    let dlnf_dt = dlogf_dl * dfc / fc;
    (ln10 * log_f, dlnf_dlnpr, dlnf_dt)
}

/// ln K_c of a reversible reaction and d ln K_c/dT, from g/(R T) and
/// h/(R T) of its species (`grt`, `hrt`).
#[inline(always)]
fn equilibrium(prog: &Program, rx: &Reaction, t: f64, grt: &[f64], hrt: &[f64]) -> (f64, f64) {
    let mut lk = 0.0;
    let mut dlk = 0.0;
    let mut dnu = 0.0;
    for j in rx.terms[0] as usize..rx.terms[1] as usize {
        let tm = prog.terms[j];
        let nu = tm.product - tm.reactant;
        let k = tm.species as usize;
        lk -= nu * grt[k];
        dlk += nu * hrt[k];
        dnu += nu;
    }
    // K_c = K_p (p°/(R T))^Δν
    lk += dnu * math::ln(P_STANDARD / (R * t));
    (lk, (dlk - dnu) / t)
}

/// The length of the work array of [`rates`] for `n` species.
pub const fn rates_work_len(n: usize) -> usize {
    4 * n
}

/// The molar production rates ω̇_k (mol/(m³ s)) of the mixture with the
/// concentrations `c` (mol/m³) at the temperature `t`: written to `wdot`.
/// With `dwdc` (N × N, row k = species k), also ∂ω̇_k/∂C_m, and with
/// `dwdt` ∂ω̇_k/∂T at fixed concentrations. `work` holds
/// [`rates_work_len`] values.
pub fn rates(
    prog: &Program,
    t: f64,
    c: &[f64],
    wdot: &mut [f64],
    mut dwdc: Option<&mut [f64]>,
    mut dwdt: Option<&mut [f64]>,
    work: &mut [f64],
) {
    let n = prog.n();
    let (grt, rest) = work.split_at_mut(n);
    let (hrt, _) = rest.split_at_mut(n);
    let ln_t = math::ln(t);
    let mut ctot = 0.0;
    for k in 0..n {
        let h = prog.species[k].h_over_rt(t);
        hrt[k] = h;
        grt[k] = h - prog.species[k].s_over_r(t);
        ctot += c[k];
        wdot[k] = 0.0;
    }
    if let Some(d) = dwdc.as_deref_mut() {
        for v in d[..n * n].iter_mut() {
            *v = 0.0;
        }
    }
    if let Some(d) = dwdt.as_deref_mut() {
        for v in d[..n].iter_mut() {
            *v = 0.0;
        }
    }
    let eps = 1e-12 * ctot.abs().max(1e-300);
    for rx in prog.reactions {
        let (t0, t1) = (rx.terms[0] as usize, rx.terms[1] as usize);
        // third-body concentration
        let third = rx.kind != ELEMENTARY;
        let mut m = 0.0;
        if third {
            m = rx.default_efficiency * ctot;
            for e in &prog.efficiencies[rx.efficiencies[0] as usize..rx.efficiencies[1] as usize] {
                m += (e.value - rx.default_efficiency) * c[e.species as usize];
            }
        }
        let (kf, dkf_dt, dkf_dm) = forward(rx, t, ln_t, m);
        let (inv_kc, dlnkc) = if rx.reversible == 1 {
            let (lk, dlk) = equilibrium(prog, rx, t, grt, hrt);
            (math::exp((-lk).max(-700.0).min(700.0)), dlk)
        } else {
            (0.0, 0.0)
        };
        let kr = kf * inv_kc;
        let dkr_dt = dkf_dt * inv_kc - kr * dlnkc;
        // the two products
        let mut pf = 1.0;
        let mut pr = 1.0;
        for tm in &prog.terms[t0..t1] {
            let ck = c[tm.species as usize];
            pf *= phi(ck, tm.order, eps).0;
            if rx.reversible == 1 {
                pr *= phi(ck, tm.order - tm.reactant + tm.product, eps).0;
            }
        }
        if rx.reversible == 0 {
            pr = 0.0;
        }
        let q = kf * pf - kr * pr;
        for tm in &prog.terms[t0..t1] {
            wdot[tm.species as usize] += (tm.product - tm.reactant) * q;
        }
        if let Some(d) = dwdt.as_deref_mut() {
            let dq = dkf_dt * pf - dkr_dt * pr;
            for tm in &prog.terms[t0..t1] {
                d[tm.species as usize] += (tm.product - tm.reactant) * dq;
            }
        }
        if let Some(d) = dwdc.as_deref_mut() {
            // ∂q/∂C_m through the products
            for (i, tmi) in prog.terms[t0..t1].iter().enumerate() {
                let mut df = 1.0;
                let mut dr = 1.0;
                for (j, tmj) in prog.terms[t0..t1].iter().enumerate() {
                    let cj = c[tmj.species as usize];
                    let of = tmj.order;
                    let or = tmj.order - tmj.reactant + tmj.product;
                    if i == j {
                        df *= phi(cj, of, eps).1;
                        dr *= phi(cj, or, eps).1;
                    } else {
                        df *= phi(cj, of, eps).0;
                        dr *= phi(cj, or, eps).0;
                    }
                }
                let dq = kf * df - if rx.reversible == 1 { kr * dr } else { 0.0 };
                if dq == 0.0 {
                    continue;
                }
                let col = tmi.species as usize;
                for tm in &prog.terms[t0..t1] {
                    d[tm.species as usize * n + col] += (tm.product - tm.reactant) * dq;
                }
            }
            // ∂q/∂C_m through [M]
            if third && dkf_dm != 0.0 {
                let g = dkf_dm * (pf - inv_kc * pr);
                if g != 0.0 {
                    let effs = &prog.efficiencies
                        [rx.efficiencies[0] as usize..rx.efficiencies[1] as usize];
                    for tm in &prog.terms[t0..t1] {
                        let nu = (tm.product - tm.reactant) * g;
                        if nu == 0.0 {
                            continue;
                        }
                        let row = tm.species as usize * n;
                        for col in 0..n {
                            d[row + col] += nu * rx.default_efficiency;
                        }
                        for e in effs {
                            d[row + e.species as usize] += nu * (e.value - rx.default_efficiency);
                        }
                    }
                }
            }
        }
    }
}

/// The density (kg/m³) and the concentrations (mol/m³) of the mass
/// fractions `y` at `p` and `t`; returns (ρ, W̄).
#[inline(always)]
pub fn concentrations(sp: &[Species], p: f64, t: f64, y: &[f64], c: &mut [f64]) -> (f64, f64) {
    let mut s = 0.0;
    for k in 0..sp.len() {
        s += y[k] / sp[k].molar_mass;
    }
    let wbar = 1.0 / s;
    let rho = p * wbar / (R * t);
    for k in 0..sp.len() {
        c[k] = rho * y[k] / sp[k].molar_mass;
    }
    (rho, wbar)
}

/// The length of the work array of [`reactor_rhs`] for `n` species.
pub const fn rhs_work_len(n: usize) -> usize {
    rates_work_len(n) + 3 * n
}

/// The right-hand side f of the reactor at the state `z` = (Y, T) and the
/// pressure `p`; with `jac` ((N+1) × (N+1), row-major, row i = ∂f_i/∂z)
/// also its Jacobian. Returns false if the state has no meaning (T ≤ 0,
/// a mixture without molar mass, a value that is not finite).
pub fn reactor_rhs(
    prog: &Program,
    p: f64,
    z: &[f64],
    f: &mut [f64],
    jac: Option<&mut [f64]>,
    work: &mut [f64],
) -> bool {
    let n = prog.n();
    let nz = n + 1;
    let t = z[n];
    if !(t > 0.0 && t.is_finite()) {
        return false;
    }
    let (c, rest) = work.split_at_mut(n);
    let (wdot, rest) = rest.split_at_mut(n);
    let (dwdt, rw) = rest.split_at_mut(n);
    let (rho, wbar) = concentrations(prog.species, p, t, &z[..n], c);
    if !(rho > 0.0 && rho.is_finite()) {
        return false;
    }
    let sp = prog.species;
    let mut cp = 0.0;
    let mut dcp = 0.0;
    for k in 0..n {
        cp += z[k] * sp[k].cp(t);
        dcp += z[k] * R * sp[k].dcp_over_r_dt(t) / sp[k].molar_mass;
    }
    if !(cp > 0.0) {
        return false;
    }
    match jac {
        None => {
            rates(prog, t, c, wdot, None, None, rw);
            let mut a = 0.0;
            for k in 0..n {
                f[k] = sp[k].molar_mass * wdot[k] / rho;
                a += R * t * sp[k].h_over_rt(t) * wdot[k];
            }
            f[n] = -a / (rho * cp);
        }
        Some(jac) => {
            // ∂ω̇/∂C into the top-left N × N block with the row stride N,
            // spread to the stride N + 1 afterwards
            rates(prog, t, c, wdot, Some(&mut jac[..n * n]), Some(dwdt), rw);
            for k in (0..n).rev() {
                for m in (0..n).rev() {
                    jac[k * nz + m] = jac[k * n + m];
                }
            }
            let mut a = 0.0;
            for k in 0..n {
                f[k] = sp[k].molar_mass * wdot[k] / rho;
                a += R * t * sp[k].h_over_rt(t) * wdot[k];
            }
            let ft = -a / (rho * cp);
            f[n] = ft;
            // dω̇_k/dY_m = ∂ω̇_k/∂C_m ρ/W_m − W̄/W_m g_k, g_k = Σ_l ∂ω̇_k/∂C_l C_l;
            // dω̇_k/dT = ∂ω̇_k/∂T − g_k/T
            for k in 0..n {
                let row = k * nz;
                let mut g = 0.0;
                for l in 0..n {
                    g += jac[row + l] * c[l];
                }
                for m in 0..n {
                    let wm = sp[m].molar_mass;
                    jac[row + m] = jac[row + m] * rho / wm - wbar / wm * g;
                }
                jac[row + n] = dwdt[k] - g / t;
            }
            // the temperature row from the species rows (still dω̇/dz)
            let rho_cp = rho * cp;
            for m in 0..nz {
                let mut s = 0.0;
                for k in 0..n {
                    s += R * t * sp[k].h_over_rt(t) * jac[k * nz + m];
                }
                if m < n {
                    let wm = sp[m].molar_mass;
                    jac[n * nz + m] = -s / rho_cp - ft * (-wbar / wm + sp[m].cp(t) / cp);
                } else {
                    let mut cw = 0.0;
                    for k in 0..n {
                        cw += R * sp[k].cp_over_r(t) * wdot[k];
                    }
                    jac[n * nz + n] = -(cw + s) / rho_cp - ft * (-1.0 / t + dcp / cp);
                }
            }
            // the species rows: ∂f_k/∂Y_m = W_k/ρ dω̇_k/dY_m + f_k W̄/W_m,
            // ∂f_k/∂T = W_k/ρ dω̇_k/dT + f_k/T
            for k in 0..n {
                let row = k * nz;
                let wk = sp[k].molar_mass / rho;
                for m in 0..n {
                    jac[row + m] = wk * jac[row + m] + f[k] * wbar / sp[m].molar_mass;
                }
                jac[row + n] = wk * jac[row + n] + f[k] / t;
            }
        }
    }
    f.iter().take(nz).all(|v| v.is_finite())
}

/// Settings of the integrator.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Settings {
    /// Relative tolerance.
    pub rtol: f64,
    /// Absolute tolerance of the mass fractions.
    pub atol: f64,
    /// Absolute tolerance of the temperature (K).
    pub atol_t: f64,
    /// The first step (s); 0: the whole interval, cut by the control.
    pub h0: f64,
    /// The largest number of steps (accepted and rejected).
    pub max_steps: u32,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            rtol: 1e-6,
            atol: 1e-12,
            atol_t: 1e-6,
            h0: 0.0,
            max_steps: 100_000,
        }
    }
}

/// The result of [`integrate`].
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Stats {
    /// Accepted steps.
    pub steps: u32,
    /// Rejected steps.
    pub rejected: u32,
    /// The last accepted step (s), a good first step of the next call.
    pub h_last: f64,
    /// 0: reached the end; 1: too many steps; 2: the step became too
    /// small or the state lost its meaning.
    pub status: u32,
}

/// The length of the work array of [`integrate`] for `n` species.
pub const fn work_len(n: usize) -> usize {
    let nz = n + 1;
    2 * nz * nz + 8 * nz + rhs_work_len(n)
}

/// LU factorisation of the `n` × `n` row-major matrix `a` in place with
/// partial pivoting (fixed column order; ties keep the upper row); the
/// row exchanges are written to `piv`. Returns false if a pivot is zero.
pub fn lu(a: &mut [f64], n: usize, piv: &mut [f64]) -> bool {
    for k in 0..n {
        let mut p = k;
        let mut best = a[k * n + k].abs();
        for i in k + 1..n {
            let v = a[i * n + k].abs();
            if v > best {
                best = v;
                p = i;
            }
        }
        piv[k] = p as f64;
        if !(best > 0.0) || !best.is_finite() {
            return false;
        }
        if p != k {
            for j in 0..n {
                a.swap(k * n + j, p * n + j);
            }
        }
        let d = 1.0 / a[k * n + k];
        for i in k + 1..n {
            let l = a[i * n + k] * d;
            a[i * n + k] = l;
            if l != 0.0 {
                for j in k + 1..n {
                    a[i * n + j] -= l * a[k * n + j];
                }
            }
        }
    }
    true
}

/// Solves `lu` x = b in place in `b` with the factors and exchanges of
/// [`lu`].
pub fn lu_solve(a: &[f64], n: usize, piv: &[f64], b: &mut [f64]) {
    for k in 0..n {
        let p = piv[k] as usize;
        if p != k {
            b.swap(k, p);
        }
    }
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= a[i * n + j] * b[j];
        }
        b[i] = s;
    }
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in i + 1..n {
            s -= a[i * n + j] * b[j];
        }
        b[i] = s / a[i * n + i];
    }
}

/// The RODAS3 coefficients (Sandu et al. 1997, table 4, in the form
/// (I/(γh) − J) K_i = f(z + Σ_j a_ij K_j) + Σ_j (c_ij/h) K_j).
const GAMMA: f64 = 0.5;
const A21: f64 = 0.0;
const A31: f64 = 2.0;
const A32: f64 = 0.0;
const A41: f64 = 2.0;
const A42: f64 = 0.0;
const A43: f64 = 1.0;
const C21: f64 = 4.0;
const C31: f64 = 1.0;
const C32: f64 = -1.0;
const C41: f64 = 1.0;
const C42: f64 = -1.0;
const C43: f64 = -8.0 / 3.0;
const M1: f64 = 2.0;
const M2: f64 = 0.0;
const M3: f64 = 1.0;
const M4: f64 = 1.0;
/// The embedded error estimate is K₄ (E = (0, 0, 0, 1)).
const ORDER_LOW: f64 = 3.0;

/// Integrates the reactor of the state `z` = (Y, T) at the pressure `p`
/// over the time `dt` (s). `work` holds [`work_len`] values.
pub fn integrate(
    prog: &Program,
    p: f64,
    z: &mut [f64],
    dt: f64,
    set: &Settings,
    work: &mut [f64],
) -> Stats {
    let n = prog.n();
    let nz = n + 1;
    let (jac, rest) = work.split_at_mut(nz * nz);
    let (mat, rest) = rest.split_at_mut(nz * nz);
    let (k1, rest) = rest.split_at_mut(nz);
    let (k2, rest) = rest.split_at_mut(nz);
    let (k3, rest) = rest.split_at_mut(nz);
    let (k4, rest) = rest.split_at_mut(nz);
    let (f0, rest) = rest.split_at_mut(nz);
    let (zt, rest) = rest.split_at_mut(nz);
    let (ft, rest) = rest.split_at_mut(nz);
    let (piv, rw) = rest.split_at_mut(nz);
    let mut st = Stats::default();
    if !(dt > 0.0) {
        return st;
    }
    let mut t = 0.0;
    let mut h = if set.h0 > 0.0 { set.h0.min(dt) } else { dt };
    let hmin = 1e-14 * dt;
    loop {
        if st.steps + st.rejected >= set.max_steps {
            st.status = 1;
            return st;
        }
        if t + h > dt || dt - (t + h) < 1e-12 * dt {
            h = dt - t;
        }
        // f and J at the start of the step
        if !reactor_rhs(prog, p, z, f0, Some(jac), rw) {
            st.status = 2;
            return st;
        }
        loop {
            // the matrix I/(γh) − J
            let d = 1.0 / (GAMMA * h);
            for i in 0..nz * nz {
                mat[i] = -jac[i];
            }
            for i in 0..nz {
                mat[i * nz + i] += d;
            }
            let ok = lu(mat, nz, piv);
            let mut good = ok;
            if ok {
                // stage 1
                k1[..nz].copy_from_slice(&f0[..nz]);
                lu_solve(mat, nz, piv, k1);
                // stage 2 (A21 = 0: the same point)
                for i in 0..nz {
                    zt[i] = z[i] + A21 * k1[i];
                    k2[i] = f0[i] + C21 / h * k1[i];
                }
                lu_solve(mat, nz, piv, k2);
                // stage 3
                for i in 0..nz {
                    zt[i] = z[i] + A31 * k1[i] + A32 * k2[i];
                }
                good = reactor_rhs(prog, p, zt, ft, None, rw);
                if good {
                    for i in 0..nz {
                        k3[i] = ft[i] + (C31 * k1[i] + C32 * k2[i]) / h;
                    }
                    lu_solve(mat, nz, piv, k3);
                    // stage 4
                    for i in 0..nz {
                        zt[i] = z[i] + A41 * k1[i] + A42 * k2[i] + A43 * k3[i];
                    }
                    good = reactor_rhs(prog, p, zt, ft, None, rw);
                }
                if good {
                    for i in 0..nz {
                        k4[i] = ft[i] + (C41 * k1[i] + C42 * k2[i] + C43 * k3[i]) / h;
                    }
                    lu_solve(mat, nz, piv, k4);
                }
            }
            // the error of the step
            let mut err = f64::INFINITY;
            if good {
                let mut s = 0.0;
                for i in 0..nz {
                    let znew = z[i] + M1 * k1[i] + M2 * k2[i] + M3 * k3[i] + M4 * k4[i];
                    zt[i] = znew;
                    let atol = if i == n { set.atol_t } else { set.atol };
                    let sc = atol + set.rtol * z[i].abs().max(znew.abs());
                    let e = k4[i] / sc;
                    s += e * e;
                }
                err = math::sqrt(s / nz as f64);
                if !(zt[n] > 0.0) || !err.is_finite() {
                    err = f64::INFINITY;
                }
            }
            let fac = if err.is_finite() && err > 0.0 {
                (0.9 * math::pow(err, -1.0 / ORDER_LOW)).max(0.2).min(6.0)
            } else if err == 0.0 {
                6.0
            } else {
                0.2
            };
            if err <= 1.0 {
                z[..nz].copy_from_slice(&zt[..nz]);
                t += h;
                st.steps += 1;
                st.h_last = h;
                if t >= dt * (1.0 - 1e-12) {
                    return st;
                }
                h *= fac;
                break;
            }
            st.rejected += 1;
            h *= fac.min(1.0);
            if h < hmin || st.steps + st.rejected >= set.max_steps {
                st.status = if h < hmin { 2 } else { 1 };
                return st;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phi_is_smooth_below_order_one() {
        let eps = 1e-12;
        for a in [0.5, -0.3, 0.25, -0.75] {
            for c in [-1e-12, 0.0, 1e-13, 1e-6, 0.3] {
                let (v, d) = phi(c, a, eps);
                let hh = 1e-6 * (c.abs() + eps);
                let fd = (phi(c + hh, a, eps).0 - phi(c - hh, a, eps).0) / (2.0 * hh);
                assert!(
                    (d - fd).abs() <= 1e-5 * (d.abs() + fd.abs()) + 1e-9,
                    "{a} {c} {d} {fd}"
                );
                if c > 1e-6 {
                    assert!((v / math::pow(c, a) - 1.0).abs() < 1e-5);
                }
            }
        }
    }

    #[test]
    fn lu_solves_with_pivoting() {
        let mut a = [0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 4.0, -1.0, 3.0];
        let x = [1.0, -2.0, 0.5];
        let mut b = [0.0; 3];
        for i in 0..3 {
            for j in 0..3 {
                b[i] += a[i * 3 + j] * x[j];
            }
        }
        let mut piv = [0.0; 3];
        assert!(lu(&mut a, 3, &mut piv));
        lu_solve(&a, 3, &piv, &mut b);
        for i in 0..3 {
            assert!((b[i] - x[i]).abs() < 1e-14);
        }
    }

    /// Two species with the same data, A → B with k = 10/s: Y_A = e^(−10 t)
    /// at a constant temperature. Fixed steps (huge tolerances) show the
    /// third order of RODAS3; the adaptive control meets its tolerance.
    #[test]
    fn rodas3_has_order_three_and_controls_its_error() {
        let sp = Species {
            molar_mass: 0.028,
            t_range: [200.0, 1000.0, 3500.0],
            low: [3.5, 0.0, 0.0, 0.0, 0.0, -1000.0, 3.0],
            high: [3.5, 0.0, 0.0, 0.0, 0.0, -1000.0, 3.0],
            well_depth: 100.0,
            diameter: 3.6,
            dipole: 0.0,
            geometry: 1,
        };
        let species = [sp, sp];
        let reactions = [Reaction {
            kind: ELEMENTARY,
            reversible: 0,
            terms: [0, 2],
            efficiencies: [0, 0],
            default_efficiency: 0.0,
            rate: [10.0, 0.0, 0.0],
            low: [0.0; 3],
            troe: [0.0; 4],
            has_t2: 0,
        }];
        let terms = [
            Term {
                species: 0,
                reactant: 1.0,
                product: 0.0,
                order: 1.0,
            },
            Term {
                species: 1,
                reactant: 0.0,
                product: 1.0,
                order: 0.0,
            },
        ];
        let prog = Program {
            species: &species,
            reactions: &reactions,
            terms: &terms,
            efficiencies: &[],
        };
        let mut work = [0.0; work_len(2)];
        let fixed = Settings {
            rtol: 1e30,
            atol: 1e30,
            atol_t: 1e30,
            h0: 0.0,
            max_steps: 10,
        };
        let err = |steps: usize, work: &mut [f64]| {
            let mut z = [1.0, 0.0, 1000.0];
            let h = 0.5 / steps as f64;
            for _ in 0..steps {
                let st = integrate(&prog, 1e5, &mut z, h, &fixed, work);
                assert_eq!(st.steps, 1);
            }
            assert!((z[2] - 1000.0).abs() < 1e-9);
            assert!((z[0] + z[1] - 1.0).abs() < 1e-14);
            (z[0] - math::exp(-5.0)).abs()
        };
        let (e1, e2, e3) = (err(10, &mut work), err(20, &mut work), err(40, &mut work));
        let (o1, o2) = ((e1 / e2).log2(), (e2 / e3).log2());
        assert!(
            o1 > 2.8 && o2 > 2.8,
            "orders {o1} {o2} ({e1:e} {e2:e} {e3:e})"
        );
        // adaptive
        let mut z = [1.0, 0.0, 1000.0];
        let st = integrate(
            &prog,
            1e5,
            &mut z,
            0.5,
            &Settings {
                rtol: 1e-8,
                atol: 1e-14,
                ..Default::default()
            },
            &mut work,
        );
        assert_eq!(st.status, 0);
        assert!((z[0] / math::exp(-5.0) - 1.0).abs() < 1e-6, "{}", z[0]);
    }
}
