//! The clamp fit (D290): the spring per length of a clamp that makes a
//! natural frequency of the model equal to a measured one.
//!
//! A stiffer spring can only raise every natural frequency (Rayleigh's
//! principle: K grows, M stays), so the frequency is monotone in the
//! spring and the root is unique. The fit brackets it by factors of ten
//! from 1 N/rad (or N/m²), then refines on ln k by the Illinois variant of
//! regula falsi (Dowell & Jarratt 1971), to a frequency within 10⁻¹⁰ of
//! the measured one. A measured frequency at or above the rigidly clamped
//! model cannot be met by any spring: the model is too soft (length,
//! thickness, modulus) and the fit says so.

use crate::analysis;
use crate::case::{self, Case};
use crate::{Error, codes};
use foehncraft_core::math;

/// Which spring of the clamp is fitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Which {
    /// The rotational spring about the clamp line.
    Rotation,
    /// The translational spring along the normal.
    Translation,
}

/// The result of a fit.
#[derive(Debug, Clone, serde::Serialize, schemars::JsonSchema)]
pub struct Fit {
    /// The fitted stiffness per metre of clamp line.
    pub k_per_length: f64,
    /// The total stiffness (times the clamp width).
    pub k_total: f64,
    /// The frequency of the model with the fitted spring (Hz).
    pub frequency: f64,
    /// The frequency with this spring rigid (Hz): the upper bound.
    pub frequency_rigid: f64,
    /// Model evaluations (modal solves).
    pub evaluations: usize,
}

fn frequency(case: &Case, mode: usize) -> Result<f64, Error> {
    let built = case::build(case)?;
    let m = analysis::solve_modal(&built.model, mode, case.analysis.shift)?;
    m.modes
        .get(mode - 1)
        .map(analysis::Mode::frequency)
        .ok_or_else(|| Error::new(codes::MODES, format!("mode {mode} not found")))
}

/// Fits the spring `which` of clamp `clamp` (0-based) so that mode `mode`
/// (1-based) has the frequency `measured` (Hz).
pub fn fit(
    case: &Case,
    clamp: usize,
    which: Which,
    mode: usize,
    measured: f64,
) -> Result<Fit, Error> {
    if clamp >= case.clamps.len() {
        return Err(Error::new(
            codes::CASE,
            format!(
                "the case has {} clamp(s); clamp {clamp} does not exist",
                case.clamps.len()
            ),
        ));
    }
    if mode == 0 || !(measured > 0.0 && measured.is_finite()) {
        return Err(Error::new(
            codes::CASE,
            "the fit needs a mode number ≥ 1 and a positive measured frequency".into(),
        ));
    }
    let with = |k: Option<f64>| {
        let mut c = case.clone();
        match which {
            Which::Rotation => c.clamps[clamp].k_rotation = k,
            Which::Translation => c.clamps[clamp].k_translation = k,
        }
        c
    };
    let mut evals = 0;
    let mut f_of = |x: f64| -> Result<f64, Error> {
        evals += 1;
        frequency(&with(Some(math::exp(x))), mode)
    };
    let f_rigid = frequency(&with(None), mode)?;
    if measured >= f_rigid * (1.0 - 1e-9) {
        return Err(Error::new(
            codes::FIT,
            format!(
                "the measured {measured} Hz is at or above the rigidly clamped model ({f_rigid:.6} Hz): no spring can reach it; the model is too soft (free length, thickness, modulus)"
            ),
        ));
    }
    // bracket on x = ln k by factors of ten
    let step = math::ln(10.0);
    let mut x0 = 0.0;
    let mut g0 = f_of(x0)? - measured;
    let mut x1 = x0;
    let mut g1 = g0;
    for _ in 0..40 {
        x1 = if g0 < 0.0 { x0 + step } else { x0 - step };
        g1 = f_of(x1)? - measured;
        if (g1 < 0.0) != (g0 < 0.0) {
            break;
        }
        x0 = x1;
        g0 = g1;
    }
    if (g1 < 0.0) == (g0 < 0.0) {
        return Err(Error::new(
            codes::FIT,
            format!(
                "no spring between 1e-40 and 1e40 gives {measured} Hz (mode {mode}): the mode does not depend on this spring"
            ),
        ));
    }
    // Illinois regula falsi: b is the newest point, a the other end of
    // the bracket; a retained end has its value halved
    let (mut a, mut fa, mut b, mut fb) = (x0, g0, x1, g1);
    let mut x = b;
    let mut fx = fb;
    for _ in 0..200 {
        x = (a * fb - b * fa) / (fb - fa);
        fx = f_of(x)? - measured;
        if fx.abs() <= 1e-10 * measured || (b - a).abs() < 1e-14 {
            break;
        }
        if (fx < 0.0) != (fb < 0.0) {
            a = b;
            fa = fb;
        } else {
            fa *= 0.5;
        }
        b = x;
        fb = fx;
    }
    let built = case::build(&with(Some(math::exp(x))))?;
    let width = built.clamps[clamp].width;
    let k = math::exp(x);
    Ok(Fit {
        k_per_length: k,
        k_total: k * width,
        frequency: measured + fx,
        frequency_rigid: f_rigid,
        evaluations: evals + 1,
    })
}
