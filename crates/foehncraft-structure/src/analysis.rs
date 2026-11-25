//! Static and modal analyses of a [`Model`].

use crate::eigen::{self, EigenError};
use crate::ldlt::Ldlt;
use crate::model::{DOF_NAMES, DOFS, Model, System};
use crate::sparse::dot;
use crate::{Error, codes};

/// Relative pivot tolerance of the factorisations.
const PIVOT_TOL: f64 = 1e-12;

fn dof_name(model: &Model, d: usize) -> String {
    let n = d / DOFS;
    let x = model.nodes.get(n).copied().unwrap_or_default();
    format!(
        "{} of node {n} at ({:.6}, {:.6}, {:.6})",
        DOF_NAMES[d % DOFS],
        x[0],
        x[1],
        x[2]
    )
}

/// The result of a static analysis.
#[derive(Debug, Clone)]
pub struct Static {
    /// Displacements and rotations of every unknown.
    pub u: Vec<f64>,
    /// Reactions at every unknown: K u − f of the supports plus the forces
    /// of the springs to ground; zero elsewhere to rounding.
    pub reactions: Vec<f64>,
    /// Strain energy ½ uᵀ K u.
    pub energy: f64,
    /// Entries of the factor L.
    pub fill: usize,
}

/// Solves K u = f.
pub fn solve_static(model: &Model) -> Result<Static, Error> {
    let sys = model.assemble();
    let f = Ldlt::factor(&sys.k, PIVOT_TOL).map_err(|z| {
        Error::new(
            codes::SINGULAR,
            format!(
                "the stiffness is singular at {} (pivot {:.3e}): the structure is free to move there; add a support or a clamp",
                dof_name(model, sys.free[z.index]),
                z.pivot
            ),
        )
    })?;
    let uf = f.solve(&sys.f);
    let u = sys.expand(&uf);
    let ku = sys.k_all.apply(&u);
    let mut reactions: Vec<f64> = ku.iter().zip(&sys.f_all).map(|(a, b)| a - b).collect();
    // springs to ground (clamp springs among them) hold the structure too:
    // their forces on it count as reactions
    for &(d, k) in &model.springs {
        reactions[d] -= k * u[d];
    }
    let energy = 0.5 * dot(&uf, &sys.k.apply(&uf));
    Ok(Static {
        u,
        reactions,
        energy,
        fill: f.nnz(),
    })
}

/// The Sturm check of a modal analysis: the number of eigenvalues in
/// [lo, hi] from the inertia of K − σM at both ends, against the number
/// of modes found there.
#[derive(Debug, Clone, Copy, serde::Serialize)]
pub struct SturmCheck {
    /// Lower end (rad²/s²).
    pub lo: f64,
    /// Upper end (rad²/s²).
    pub hi: f64,
    /// Eigenvalues in [lo, hi] by the inertia.
    pub expected: usize,
    /// Modes found in [lo, hi].
    pub found: usize,
}

impl SturmCheck {
    /// No mode was missed.
    pub fn ok(&self) -> bool {
        self.expected == self.found
    }
}

/// One mode.
#[derive(Debug, Clone)]
pub struct Mode {
    /// Eigenvalue ω² (rad²/s²).
    pub lambda: f64,
    /// M-normalised shape on every unknown.
    pub shape: Vec<f64>,
    /// Effective mass in x, y, z (kg).
    pub effective_mass: [f64; 3],
    /// Residual estimate of Lanczos.
    pub residual: f64,
}

impl Mode {
    /// Natural frequency (Hz); a slightly negative λ of a free body is 0.
    pub fn frequency(&self) -> f64 {
        self.lambda.max(0.0).sqrt() / (2.0 * std::f64::consts::PI)
    }
}

/// The result of a modal analysis.
#[derive(Debug, Clone)]
pub struct Modal {
    /// Modes by increasing frequency.
    pub modes: Vec<Mode>,
    /// Lanczos steps.
    pub steps: usize,
    /// The Sturm check.
    pub sturm: SturmCheck,
    /// Every mode met the residual tolerance.
    pub converged: bool,
    /// Total translational mass (kg).
    pub total_mass: f64,
    /// Warnings (code STR003).
    pub warnings: Vec<String>,
}

fn eigen_error(model: &Model, sys: &System, e: EigenError, shift: f64) -> Error {
    match e {
        EigenError::Singular(z) => Error::new(
            codes::SINGULAR,
            format!(
                "K − σM is singular at {} for the shift σ = {shift:.6e} (pivot {:.3e}): the shift is an eigenvalue, or the structure is free there; for a free structure give a negative `shift`",
                dof_name(model, sys.free[z.index]),
                z.pivot
            ),
        ),
        EigenError::Tridiagonal => Error::new(
            codes::MODES,
            "the tridiagonal eigenvalue iteration did not converge".to_string(),
        ),
    }
}

/// Number of eigenvalues below σ: negative pivots of K − σM.
fn below(sys: &System, sigma: f64) -> Option<usize> {
    let a = sys.k.plus_scaled(-sigma, &sys.m);
    Ldlt::factor(&a, PIVOT_TOL)
        .ok()
        .map(|f| f.negative_pivots())
}

/// The `count` modes closest to `shift` (rad²/s²; 0 for the lowest modes
/// of a supported structure, a small negative value for a free one).
pub fn solve_modal(model: &Model, count: usize, shift: f64) -> Result<Modal, Error> {
    let sys = model.assemble();
    let nf = sys.free.len();
    if count == 0 || count > nf {
        return Err(Error::new(
            codes::CASE,
            format!("`modes` must be between 1 and the number of free unknowns ({nf})"),
        ));
    }
    let tol = 1e-10;
    let mut want = count;
    let mut warnings = Vec::new();
    let mut pairs;
    let mut sturm;
    loop {
        pairs = eigen::lanczos(&sys.k, &sys.m, want, shift, tol)
            .map_err(|e| eigen_error(model, &sys, e, shift))?;
        // keep the `count` closest to the shift
        let mut idx: Vec<usize> = (0..pairs.values.len()).collect();
        idx.sort_by(|&a, &b| {
            (pairs.values[a] - shift)
                .abs()
                .total_cmp(&(pairs.values[b] - shift).abs())
        });
        idx.truncate(count);
        idx.sort_unstable();
        let vals: Vec<f64> = idx.iter().map(|&i| pairs.values[i]).collect();
        let (vmin, vmax) = (vals[0], vals[vals.len() - 1]);
        let gap = 1e-6 * (vmax.abs() + vmin.abs()) + 1e-12;
        let lo = vmin - gap;
        let hi = vmax + gap;
        let expected = match (below(&sys, lo), below(&sys, hi)) {
            (Some(a), Some(b)) => b.saturating_sub(a),
            _ => vals.len(),
        };
        let found = vals.len();
        sturm = SturmCheck {
            lo,
            hi,
            expected,
            found,
        };
        pairs = eigen::EigenPairs {
            values: vals,
            vectors: idx.iter().map(|&i| pairs.vectors[i].clone()).collect(),
            residuals: idx.iter().map(|&i| pairs.residuals[i]).collect(),
            steps: pairs.steps,
            converged: pairs.converged,
        };
        if sturm.ok() || want >= nf || want >= 4 * count + 20 {
            break;
        }
        want = (want * 2).min(nf);
    }
    if !sturm.ok() {
        warnings.push(format!(
            "the Sturm count finds {} eigenvalues between {:.6e} and {:.6e} rad²/s², Lanczos {}: a mode may be missing",
            sturm.expected, sturm.lo, sturm.hi, sturm.found
        ));
    }
    if !pairs.converged {
        warnings.push(format!(
            "Lanczos stopped after {} steps with a residual estimate above {tol:e}",
            pairs.steps
        ));
    }
    // effective masses: Γ_d = φᵀ M r_d with r_d the unit translation in d
    let rig: Vec<Vec<f64>> = (0..3)
        .map(|d| {
            let mut r = vec![0.0; model.dofs()];
            for n in 0..model.nodes.len() {
                r[DOFS * n + d] = 1.0;
            }
            sys.restrict(&r)
        })
        .collect();
    let modes = pairs
        .values
        .iter()
        .zip(&pairs.vectors)
        .zip(&pairs.residuals)
        .map(|((&lambda, v), &residual)| {
            let mv = sys.m.apply(v);
            let mut em = [0.0; 3];
            for d in 0..3 {
                let g = dot(&mv, &rig[d]);
                em[d] = g * g;
            }
            Mode {
                lambda,
                shape: sys.expand(v),
                effective_mass: em,
                residual,
            }
        })
        .collect();
    Ok(Modal {
        modes,
        steps: pairs.steps,
        sturm,
        converged: pairs.converged,
        total_mass: model.total_mass(),
        warnings,
    })
}
