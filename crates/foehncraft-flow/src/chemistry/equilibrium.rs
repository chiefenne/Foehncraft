//! Chemical equilibrium of an ideal-gas mixture by minimisation of the
//! Gibbs energy (SF2), with the same NASA-7 data as the kinetics, so that
//! the end state of a reactor and the equilibrium are comparable to the
//! integrator tolerance.
//!
//! At fixed T and p the mole numbers n_j (mol/kg) minimise
//! G/(R T) = Σ_j n_j (g_j°/(R T) + ln(n_j/n) + ln(p/p°)) subject to the
//! element masses Σ_j a_ij n_j = b_i. The iteration is the reduced Newton
//! method of Gordon & McBride (1994), "Computer program for calculation of
//! complex chemical equilibrium compositions and applications", NASA
//! RP-1311, eqs. (2.24) and (2.26) for the Lagrange multipliers π_i and
//! Δ ln n, the update Δ ln n_j = −μ_j/(R T) + Σ_i a_ij π_i + Δ ln n
//! (eq. 2.18), and the damping of section 3.3 (eqs. 3.1, 3.2). The mole
//! numbers are carried as logarithms, so trace species stay positive. At
//! fixed h and p the temperature is found by a bracketed secant
//! iteration on h(T) of the equilibrium mixture.

use super::mechanism::Mechanism;
use anyhow::Result;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_core::math;
use foehncraft_kernels::chemistry::P_STANDARD;
use foehncraft_kernels::species::R;

/// An equilibrium state.
#[derive(Debug, Clone, PartialEq, serde::Serialize, schemars::JsonSchema)]
pub struct State {
    /// Temperature (K).
    pub temperature: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Mole fractions in the order of the mechanism's species.
    pub mole_fractions: Vec<f64>,
    /// Mass fractions in the same order.
    pub mass_fractions: Vec<f64>,
    /// Newton iterations of the last composition solve.
    pub iterations: usize,
}

/// The element amounts b_i (mol of atoms per kg) of the mass fractions
/// `y`.
pub fn elements_of(m: &Mechanism, y: &[f64]) -> Vec<f64> {
    m.atoms
        .iter()
        .map(|a| {
            a.iter()
                .zip(y)
                .zip(&m.data)
                .map(|((a, y), s)| a * y / s.molar_mass)
                .sum()
        })
        .collect()
}

/// The specific enthalpy (J/kg) of the mass fractions `y` at `t`.
pub fn enthalpy(m: &Mechanism, y: &[f64], t: f64) -> f64 {
    y.iter()
        .zip(&m.data)
        .map(|(y, s)| y * R * t * s.h_over_rt(t) / s.molar_mass)
        .sum()
}

/// The equilibrium at the temperature `t` and pressure `p` of the mixture
/// with the mass fractions `y`; `start` (ln n_j) warm-starts the solve.
pub fn tp(m: &Mechanism, y: &[f64], t: f64, p: f64, start: Option<&mut Vec<f64>>) -> Result<State> {
    let b0 = elements_of(m, y);
    let mut own = Vec::new();
    let ln_n = match start {
        Some(s) => s,
        None => &mut own,
    };
    let it = solve_tp(m, &b0, t, p, ln_n)?;
    Ok(state(m, ln_n, t, p, it))
}

fn state(m: &Mechanism, ln_n: &[f64], t: f64, p: f64, iterations: usize) -> State {
    let n: Vec<f64> = ln_n.iter().map(|l| math::exp(*l)).collect();
    let tot: f64 = n.iter().sum();
    let mass: f64 = n.iter().zip(&m.data).map(|(n, s)| n * s.molar_mass).sum();
    State {
        temperature: t,
        pressure: p,
        mole_fractions: n.iter().map(|v| v / tot).collect(),
        mass_fractions: n
            .iter()
            .zip(&m.data)
            .map(|(n, s)| n * s.molar_mass / mass)
            .collect(),
        iterations,
    }
}

/// The ln n_j of a species that cannot exist (an element it needs is
/// absent).
const ABSENT: f64 = -1e300;

fn solve_tp(m: &Mechanism, b0: &[f64], t: f64, p: f64, ln_n: &mut Vec<f64>) -> Result<usize> {
    let ns = m.len();
    // the elements present and the species that can exist
    let present: Vec<usize> = (0..b0.len()).filter(|&i| b0[i] > 1e-300).collect();
    let can = |j: usize| (0..b0.len()).all(|i| present.contains(&i) || m.atoms[i][j] == 0.0);
    let g0: Vec<f64> = m
        .data
        .iter()
        .map(|s| s.h_over_rt(t) - s.s_over_r(t) + math::ln(p / P_STANDARD))
        .collect();
    if ln_n.len() != ns || ln_n.iter().all(|v| *v <= ABSENT) {
        ln_n.clear();
        ln_n.extend((0..ns).map(|j| {
            if can(j) {
                math::ln(0.1 / ns as f64)
            } else {
                ABSENT
            }
        }));
    } else {
        for (j, v) in ln_n.iter_mut().enumerate() {
            if !can(j) {
                *v = ABSENT;
            } else if *v < -700.0 {
                *v = -700.0;
            }
        }
    }
    let ne = present.len();
    let ln_tot: f64 = ln_n
        .iter()
        .filter(|v| **v > ABSENT)
        .map(|v| math::exp(*v))
        .exact_sum();
    let mut ln_tot = math::ln(ln_tot);
    let dim = ne + 1;
    let mut a = vec![0.0; dim * dim];
    let mut r = vec![0.0; dim];
    let mut piv = vec![0.0; dim];
    for iter in 1..=500 {
        let n: Vec<f64> = ln_n
            .iter()
            .map(|v| if *v > ABSENT { math::exp(*v) } else { 0.0 })
            .collect();
        let tot = math::exp(ln_tot);
        let mu: Vec<f64> = (0..ns)
            .map(|j| {
                if n[j] > 0.0 {
                    g0[j] + ln_n[j] - ln_tot
                } else {
                    0.0
                }
            })
            .collect();
        a.fill(0.0);
        r.fill(0.0);
        for (ii, &i) in present.iter().enumerate() {
            let ai = &m.atoms[i];
            let mut bi = 0.0;
            let mut s = 0.0;
            for j in 0..ns {
                bi += ai[j] * n[j];
                s += ai[j] * n[j] * mu[j];
            }
            for (kk, &k) in present.iter().enumerate() {
                let ak = &m.atoms[k];
                let mut g = 0.0;
                for j in 0..ns {
                    g += ai[j] * ak[j] * n[j];
                }
                a[ii * dim + kk] = g;
            }
            a[ii * dim + ne] = bi;
            a[ne * dim + ii] = bi;
            r[ii] = b0[i] - bi + s;
        }
        let sn: f64 = n.iter().sum();
        a[ne * dim + ne] = sn - tot;
        r[ne] = tot - sn + n.iter().zip(&mu).map(|(n, m)| n * m).sum::<f64>();
        if !foehncraft_kernels::chemistry::lu(&mut a, dim, &mut piv) {
            return Err(super::mechanism::bad(
                "equilibrium",
                "the equilibrium matrix is singular (an element without species?)",
            ));
        }
        foehncraft_kernels::chemistry::lu_solve(&a, dim, &piv, &mut r);
        let dln_tot = r[ne];
        let mut dln = vec![0.0; ns];
        for j in 0..ns {
            if n[j] == 0.0 {
                continue;
            }
            let mut s = -mu[j] + dln_tot;
            for (ii, &i) in present.iter().enumerate() {
                s += m.atoms[i][j] * r[ii];
            }
            dln[j] = s;
        }
        // damping (Gordon & McBride 1994, eqs. 3.1, 3.2)
        let mut big = 5.0 * dln_tot.abs();
        let mut lam2 = f64::INFINITY;
        for j in 0..ns {
            if n[j] == 0.0 {
                continue;
            }
            let rel = ln_n[j] - ln_tot;
            if rel > -18.420_681 {
                big = big.max(dln[j].abs());
            } else if dln[j] >= 0.0 && dln[j] - dln_tot != 0.0 {
                lam2 = lam2.min(((-rel - 9.210_340_4) / (dln[j] - dln_tot)).abs());
            }
        }
        let lam1 = if big > 2.0 { 2.0 / big } else { 1.0 };
        let lam = lam1.min(lam2).min(1.0);
        let mut conv = (dln_tot.abs()) < 5e-9;
        for j in 0..ns {
            if n[j] == 0.0 {
                continue;
            }
            ln_n[j] += lam * dln[j];
            if n[j] * dln[j].abs() / sn > 5e-9 {
                conv = false;
            }
        }
        ln_tot += lam * dln_tot;
        if conv && lam == 1.0 {
            // the element balance as well
            let n2: Vec<f64> = ln_n
                .iter()
                .map(|v| if *v > ABSENT { math::exp(*v) } else { 0.0 })
                .collect();
            let ok = present.iter().all(|&i| {
                let bi: f64 = (0..ns).map(|j| m.atoms[i][j] * n2[j]).sum();
                (bi - b0[i]).abs() <= 1e-10 * b0[i].abs().max(1e-300)
            });
            if ok {
                return Ok(iter);
            }
        }
    }
    Err(super::mechanism::bad(
        "equilibrium",
        format!("the equilibrium at T = {t} K did not converge in 500 iterations"),
    ))
}

/// The equilibrium at the specific enthalpy and pressure of the mixture
/// `y` at `t0` (the adiabatic flame state at constant pressure).
pub fn hp(m: &Mechanism, y: &[f64], t0: f64, p: f64) -> Result<State> {
    let h0 = enthalpy(m, y, t0);
    let mut ln_n = Vec::new();
    let f = |t: f64, ln_n: &mut Vec<f64>| -> Result<(f64, State)> {
        let s = tp(m, y, t, p, Some(ln_n))?;
        Ok((enthalpy(m, &s.mass_fractions, t) - h0, s))
    };
    // bracket: h(T) of the equilibrium mixture increases with T
    let (mut lo, mut hi) = (t0.min(300.0).max(200.0), 6000.0f64.max(t0));
    let (mut flo, _) = f(lo, &mut ln_n)?;
    let (mut fhi, _) = f(hi, &mut ln_n)?;
    if flo > 0.0 || fhi < 0.0 {
        return Err(super::mechanism::bad(
            "equilibrium",
            format!("no adiabatic equilibrium temperature between {lo} and {hi} K"),
        ));
    }
    let mut t = t0.max(lo).min(hi);
    let mut last = None;
    for _ in 0..200 {
        let (ft, s) = f(t, &mut ln_n)?;
        if ft.abs() <= 1e-9 * h0.abs().max(1e3) || hi - lo < 1e-9 {
            return Ok(s);
        }
        if ft < 0.0 {
            lo = t;
            flo = ft;
        } else {
            hi = t;
            fhi = ft;
        }
        // secant (regula falsi with the Illinois correction), else bisection
        let mut next = lo - flo * (hi - lo) / (fhi - flo);
        if !(next > lo && next < hi) {
            next = 0.5 * (lo + hi);
        }
        if last == Some(ft < 0.0) {
            next = 0.5 * (next + 0.5 * (lo + hi));
        }
        last = Some(ft < 0.0);
        t = next;
    }
    Err(super::mechanism::bad(
        "equilibrium",
        "the adiabatic equilibrium temperature did not converge",
    ))
}
