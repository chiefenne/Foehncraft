//! Unilateral contact of nodes with rigid planes: a reed on its seat and
//! against its stop (MS5 part 4).
//!
//! A plane through `point` with the unit `normal` pointing to the free
//! side lets a node be anywhere with the gap g = (x + u − point)·n − offset
//! ≥ 0 (`offset`: half the thickness of a plate or beam, so that the face
//! touches, not the mid-surface). The contact force λ ≥ 0 pushes the node
//! along n, and λ g = 0.
//!
//! - **Static**: the augmented Lagrangian method (Uzawa iterations on a
//!   penalty; Simo & Laursen 1992, *Computers & Structures* 42; Wriggers,
//!   *Computational Contact Mechanics*, 2nd ed. 2006, section 6.3). For a
//!   given λ̄ the contact force of a pair is max(0, λ̄ − k g); an active
//!   set loop finds the pairs where it is positive (K plus the penalty on
//!   them, factored once per set), then λ̄ takes the force and the loop
//!   repeats until the largest penetration is below the tolerance: the
//!   answer is the exact one, not the penalty's.
//! - **In time** (generalised-α): the penalty force k max(0, −g) at the
//!   end of the step, its α-average (1 − α_f) λ_{n+1} + α_f λ_n entering
//!   the balance like the other forces (the form of generalised-α with
//!   interpolated forces; Erlicher, Bonaventura & Bursi 2002,
//!   *Computational Mechanics* 28). Holding the α-point on the plane
//!   instead would place the end of the step α_f/(1 − α_f) of a step's
//!   travel inside it and feed energy into the next impact. An active set
//!   loop runs in each step (the factor is kept while the set does not
//!   change). The penetration is about the contact force over k; the
//!   result reports the largest one at the end of a step. The impact is
//!   elastic: at the stiffness of a penalty on a node of a shell it lasts
//!   far less than a step, and the numerical damping of generalised-α
//!   (ρ∞ < 1) takes the energy of the high frequencies it excites.
//!
//! The penalty k is `penalty` times the largest stiffness of a contact
//! node along its plane's normal (nᵀ diag(K) n), so that a stiff
//! direction of the node that the plane does not touch (the axial one of
//! a beam) does not spoil the conditioning.

use crate::analysis::{PIVOT_TOL, dof_name};
use crate::ldlt::Ldlt;
use crate::model::{DOFS, Model, System};
use crate::sparse::{Assembly, SymCsc};
use crate::{Error, codes};

/// A rigid plane.
#[derive(Debug, Clone)]
pub struct Plane {
    /// Name (for the results).
    pub name: String,
    /// A point of the plane.
    pub point: [f64; 3],
    /// Unit normal, pointing to the free side.
    pub normal: [f64; 3],
    /// Distance kept between a node and the plane (m).
    pub offset: f64,
    /// The nodes that may touch it.
    pub nodes: Vec<usize>,
}

/// One node against one plane.
#[derive(Debug, Clone, Copy)]
pub struct Pair {
    /// Node.
    pub node: usize,
    /// Plane (index into the planes).
    pub plane: usize,
    /// Gap at u = 0 (and at the prescribed values of fixed unknowns).
    pub g0: f64,
    /// Free indices of the node's translations with the normal's
    /// components: g = g0 + Σ n_d u_d.
    pub terms: [(Option<usize>, f64); 3],
}

/// The contact pairs of a model on its free unknowns.
#[derive(Debug, Clone)]
pub struct Contacts {
    /// The planes.
    pub planes: Vec<Plane>,
    /// The pairs.
    pub pairs: Vec<Pair>,
    /// The penalty stiffness (N/m).
    pub k: f64,
}

/// The contact state of a solution.
#[derive(Debug, Clone, Default)]
pub struct State {
    /// Contact force of each pair (N, along the plane's normal).
    pub force: Vec<f64>,
    /// Gap of each pair (m; negative: penetration).
    pub gap: Vec<f64>,
}

impl State {
    /// Total contact force on each plane's nodes (N, along its normal), and
    /// the number of nodes in contact.
    pub fn per_plane(&self, c: &Contacts) -> Vec<(f64, usize)> {
        let mut out = vec![(0.0, 0); c.planes.len()];
        for (p, pair) in c.pairs.iter().enumerate() {
            if self.force[p] > 0.0 {
                out[pair.plane].0 += self.force[p];
                out[pair.plane].1 += 1;
            }
        }
        out
    }

    /// The largest penetration (m, ≥ 0).
    pub fn penetration(&self) -> f64 {
        self.gap.iter().fold(0.0f64, |a, &g| a.max(-g))
    }
}

impl Contacts {
    /// The pairs of `planes` on the system of `model`; `penalty` scales the
    /// largest translational diagonal of K at the contact nodes. A node
    /// that starts inside a plane (by more than `tol`) is an error.
    pub fn new(
        model: &Model,
        sys: &System,
        planes: Vec<Plane>,
        penalty: f64,
        tol: f64,
    ) -> Result<Self, Error> {
        if !(penalty > 0.0 && penalty.is_finite()) {
            return Err(Error::new(
                codes::CASE,
                format!("[contact] penalty must be positive (is {penalty})"),
            ));
        }
        let diag = sys.k.diagonal();
        let mut kmax = 0.0f64;
        let mut pairs = Vec::new();
        for (j, pl) in planes.iter().enumerate() {
            for &n in &pl.nodes {
                let x = model.nodes[n];
                let mut g0 = -pl.offset;
                let mut terms = [(None, 0.0); 3];
                // the node's stiffness along the normal, nᵀ diag(K) n
                let mut kn = 0.0;
                for d in 0..3 {
                    let dof = DOFS * n + d;
                    g0 += (x[d] - pl.point[d] + sys.prescribed[dof]) * pl.normal[d];
                    terms[d] = (sys.map[dof], pl.normal[d]);
                    if let Some(i) = sys.map[dof] {
                        kn += pl.normal[d] * pl.normal[d] * diag[i];
                    }
                }
                kmax = kmax.max(kn);
                if g0 < -tol {
                    return Err(Error::new(
                        codes::CASE,
                        format!(
                            "contact `{}`: the node at ({}, {}, {}) starts {g0:.3e} m inside the plane",
                            pl.name, x[0], x[1], x[2]
                        ),
                    ));
                }
                pairs.push(Pair {
                    node: n,
                    plane: j,
                    g0,
                    terms,
                });
            }
        }
        Ok(Contacts {
            planes,
            pairs,
            k: penalty * kmax.max(f64::MIN_POSITIVE),
        })
    }

    /// The gap of pair `p` for the free unknowns `u`.
    pub fn gap(&self, p: usize, u: &[f64]) -> f64 {
        let pr = &self.pairs[p];
        pr.g0
            + pr
                .terms
                .iter()
                .filter_map(|&(i, n)| i.map(|i| n * u[i]))
                .sum::<f64>()
    }

    /// `Σ_active w_p Bᵀ B` on the free unknowns.
    fn penalty_matrix(&self, nf: usize, active: &[bool], w: impl Fn(usize) -> f64) -> SymCsc {
        let mut a = Assembly::new(nf);
        for (p, pr) in self.pairs.iter().enumerate() {
            if !active[p] {
                continue;
            }
            let dofs: Vec<Option<usize>> = pr.terms.iter().map(|t| t.0).collect();
            let mut b = [0.0; 9];
            let wp = w(p);
            for r in 0..3 {
                for c in 0..3 {
                    b[3 * r + c] = wp * pr.terms[r].1 * pr.terms[c].1;
                }
            }
            a.add_block(&dofs, &b);
        }
        a.finish()
    }

    /// Adds `Bᵀ λ` of pair `p` to `f`.
    fn add_force(&self, p: usize, lambda: f64, f: &mut [f64]) {
        for &(i, n) in &self.pairs[p].terms {
            if let Some(i) = i {
                f[i] += n * lambda;
            }
        }
    }

    /// The contact forces as a vector on every unknown of the model.
    pub fn force_vector(&self, st: &State, ndofs: usize) -> Vec<f64> {
        let mut f = vec![0.0; ndofs];
        for (p, pr) in self.pairs.iter().enumerate() {
            let n = self.planes[pr.plane].normal;
            for d in 0..3 {
                f[DOFS * pr.node + d] += n[d] * st.force[p];
            }
        }
        f
    }
}

fn factor(model: &Model, sys: &System, a: &SymCsc) -> Result<Ldlt, Error> {
    Ldlt::factor(a, PIVOT_TOL).map_err(|z| {
        Error::new(
            codes::SINGULAR,
            format!(
                "the stiffness is singular at {} (pivot {:.3e}): the structure is free to move there and no contact holds it; add a support or a clamp",
                dof_name(model, sys.free[z.index]),
                z.pivot
            ),
        )
    })
}

/// The static solution with contact: free displacements and the contact
/// state. `tol` is the largest penetration accepted (m).
pub fn solve_static(
    model: &Model,
    sys: &System,
    c: &Contacts,
    tol: f64,
) -> Result<(Vec<f64>, State, usize), Error> {
    let nf = sys.free.len();
    let np = c.pairs.len();
    let mut lam = vec![0.0; np];
    // start with the pairs that touch
    let mut active: Vec<bool> = c.pairs.iter().map(|p| p.g0 <= 0.0).collect();
    let mut solves = 0;
    let mut cached: Option<(Vec<bool>, Ldlt)> = None;
    for _outer in 0..100 {
        let mut u = vec![0.0; nf];
        let mut settled = false;
        for _inner in 0..100 {
            if cached.as_ref().is_none_or(|(a, _)| *a != active) {
                let k = sys.k.plus_scaled(1.0, &c.penalty_matrix(nf, &active, |_| c.k));
                cached = Some((active.clone(), factor(model, sys, &k)?));
            }
            let mut f = sys.f.clone();
            for p in 0..np {
                if active[p] {
                    c.add_force(p, lam[p] - c.k * c.pairs[p].g0, &mut f);
                }
            }
            u = cached.as_ref().map(|(_, l)| l.solve(&f)).unwrap_or_default();
            solves += 1;
            let next: Vec<bool> = (0..np)
                .map(|p| lam[p] - c.k * c.gap(p, &u) > 0.0)
                .collect();
            if next == active {
                settled = true;
                break;
            }
            active = next;
        }
        if !settled {
            return Err(Error::new(
                codes::CASE,
                "contact: the set of nodes in contact does not settle (100 solves)".into(),
            ));
        }
        let gaps: Vec<f64> = (0..np).map(|p| c.gap(p, &u)).collect();
        let force: Vec<f64> = (0..np)
            .map(|p| (lam[p] - c.k * gaps[p]).max(0.0))
            .collect();
        let pen = gaps.iter().fold(0.0f64, |a, &g| a.max(-g));
        if pen <= tol {
            return Ok((u, State { force, gap: gaps }, solves));
        }
        lam = force;
    }
    Err(Error::new(
        codes::CASE,
        "contact: the penetration does not fall below the tolerance in 100 updates; raise `[contact] penalty`".into(),
    ))
}

/// The active set and its factor for one step of generalised-α, and the
/// contact forces of the last step.
#[derive(Default)]
pub struct StepCache {
    active: Vec<bool>,
    factor: Option<Ldlt>,
    last: Vec<f64>,
    /// Factorisations so far.
    pub factorisations: usize,
}

/// The end of a generalised-α step: the displacement at a = 0 (the
/// predictor), its coefficient of the acceleration, and α_f.
pub struct StepEnd<'a> {
    /// Displacement at a = 0.
    pub u0: &'a [f64],
    /// ∂u_{n+1}/∂a (β Δt²).
    pub cu: f64,
    /// α_f.
    pub af: f64,
}

/// One generalised-α step with contact: `base` is S/ck (K + cm/ck M),
/// `ck` its scale, `at` the end of the step and `rhs` the right-hand side
/// without contact (the same as the step without contact). Returns the
/// acceleration and the contact state at the end of the step.
#[allow(clippy::too_many_arguments)]
pub fn step(
    model: &Model,
    sys: &System,
    c: &Contacts,
    base: &SymCsc,
    ck: f64,
    at: &StepEnd<'_>,
    rhs: &[f64],
    cache: &mut StepCache,
) -> Result<(Vec<f64>, State), Error> {
    let (u0, cu, af) = (at.u0, at.cu, at.af);
    let nf = sys.free.len();
    let np = c.pairs.len();
    if cache.active.len() != np {
        cache.active = vec![false; np];
        cache.factor = None;
        cache.last = vec![0.0; np];
    }
    // start from the pairs that are closed at the predicted end
    let mut active: Vec<bool> = (0..np).map(|p| c.gap(p, u0) < 0.0).collect();
    for _ in 0..50 {
        if cache.factor.is_none() || cache.active != active {
            let a = base.plus_scaled(
                1.0,
                &c.penalty_matrix(nf, &active, |_| (1.0 - af) * cu * c.k / ck),
            );
            cache.factor = Some(factor(model, sys, &a)?);
            cache.active = active.clone();
            cache.factorisations += 1;
        }
        let mut r = rhs.to_vec();
        for p in 0..np {
            // α_f λ_n, and (1 − α_f)(−k g) at a = 0 on the active pairs;
            // the part in a is in the matrix
            let mut f0 = af * cache.last[p];
            if active[p] {
                f0 -= (1.0 - af) * c.k * c.gap(p, u0);
            }
            if f0 != 0.0 {
                c.add_force(p, f0, &mut r);
            }
        }
        let a: Vec<f64> = cache
            .factor
            .as_ref()
            .map(|f| f.solve(&r))
            .unwrap_or_default()
            .into_iter()
            .map(|x| x / ck)
            .collect();
        let u: Vec<f64> = (0..nf).map(|i| u0[i] + cu * a[i]).collect();
        let gaps: Vec<f64> = (0..np).map(|p| c.gap(p, &u)).collect();
        let next: Vec<bool> = gaps.iter().map(|&g| g < 0.0).collect();
        if next == active {
            let force: Vec<f64> = gaps.iter().map(|&g| (-c.k * g).max(0.0)).collect();
            cache.last.clone_from(&force);
            return Ok((a, State { force, gap: gaps }));
        }
        active = next;
    }
    Err(Error::new(
        codes::CASE,
        "contact: the set of nodes in contact does not settle within a time step (50 solves); take a smaller `time_step`".into(),
    ))
}
