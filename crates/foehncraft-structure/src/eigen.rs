//! Eigenvalues: the symmetric tridiagonal QL method and shift-invert
//! Lanczos for K φ = λ M φ.
//!
//! **Tridiagonal.** Implicit QL with Wilkinson's shift and the eigenvectors
//! accumulated (Bowdler, Martin, Reinsch & Wilkinson 1968, *Handbook for
//! Automatic Computation* II/3; Parlett, *The Symmetric Eigenvalue
//! Problem*, 1998, ch. 8).
//!
//! **Lanczos.** The operator is (K − σM)⁻¹M, symmetric in the M inner
//! product; its largest eigenvalues θ = 1/(λ − σ) are the λ closest to the
//! shift σ (Ericsson & Ruhe 1980; Grimes, Lewis & Simon 1994). Every new
//! vector is M-orthogonalised against all earlier ones twice (full
//! reorthogonalisation, Parlett ch. 13), which is cheap for the small
//! models of this crate and keeps the basis orthogonal to rounding. The
//! start vector is the operator applied to a fixed pseudo-random vector,
//! so it carries no component of a massless unknown. A Ritz pair is
//! accepted when its residual estimate |β_j s_j| is below `tol·|θ|`.

use crate::ldlt::{Ldlt, ZeroPivot};
use crate::sparse::{SymCsc, dot};

/// √(a² + b²) without overflow.
fn hypot(a: f64, b: f64) -> f64 {
    let (a, b) = (a.abs(), b.abs());
    let (big, small) = if a > b { (a, b) } else { (b, a) };
    if big == 0.0 {
        return 0.0;
    }
    let r = small / big;
    big * (1.0 + r * r).sqrt()
}

/// Eigenvalues and eigenvectors of the symmetric tridiagonal matrix with
/// diagonal `d` and off-diagonal `e` (`e[i]` couples `i` and `i + 1`).
/// Returns the eigenvalues (unsorted) and the eigenvectors as columns of
/// `z` (row major, `n × n`), or `None` if an eigenvalue takes more than
/// 60 iterations.
pub fn tridiagonal(d: &[f64], e: &[f64]) -> Option<(Vec<f64>, Vec<f64>)> {
    let n = d.len();
    let mut d = d.to_vec();
    let mut e: Vec<f64> = (0..n).map(|i| if i + 1 < n { e[i] } else { 0.0 }).collect();
    let mut z = vec![0.0; n * n];
    for i in 0..n {
        z[i * n + i] = 1.0;
    }
    for l in 0..n {
        let mut iter = 0;
        loop {
            let mut m = l;
            while m + 1 < n {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= f64::EPSILON * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iter += 1;
            if iter > 60 {
                return None;
            }
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = hypot(g, 1.0);
            g = d[m] - d[l] + e[l] / (g + if g >= 0.0 { r } else { -r });
            let (mut s, mut c, mut p) = (1.0, 1.0, 0.0);
            let mut deflated = false;
            let mut i = m;
            while i > l {
                i -= 1;
                let f = s * e[i];
                let b = c * e[i];
                r = hypot(f, g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    deflated = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
                for k in 0..n {
                    let f = z[k * n + i + 1];
                    z[k * n + i + 1] = s * z[k * n + i] + c * f;
                    z[k * n + i] = c * z[k * n + i] - s * f;
                }
            }
            if deflated {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    Some((d, z))
}

/// Why the eigen solve failed.
#[derive(Debug, Clone, PartialEq)]
pub enum EigenError {
    /// K − σM is singular: the shift is an eigenvalue (or K and M share a
    /// null space).
    Singular(ZeroPivot),
    /// The tridiagonal QL did not converge.
    Tridiagonal,
}

/// Eigenpairs from [`lanczos`].
#[derive(Debug, Clone)]
pub struct EigenPairs {
    /// Eigenvalues λ = ω², ascending.
    pub values: Vec<f64>,
    /// M-normalised eigenvectors, one per value.
    pub vectors: Vec<Vec<f64>>,
    /// Residual estimate of each pair, |β s| / |θ|.
    pub residuals: Vec<f64>,
    /// Lanczos steps taken.
    pub steps: usize,
    /// Every requested pair met the tolerance.
    pub converged: bool,
}

fn start_vector(n: usize, seed: u64) -> Vec<f64> {
    (0..n)
        .map(|i| foehncraft_core::rng::uniform([0x5354_5255_4354, seed, i as u64]) - 0.5)
        .collect()
}

/// The `count` eigenpairs of K φ = λ M φ closest to `shift`.
pub fn lanczos(
    k: &SymCsc,
    m: &SymCsc,
    count: usize,
    shift: f64,
    tol: f64,
) -> Result<EigenPairs, EigenError> {
    let n = k.n;
    let op = k.plus_scaled(-shift, m);
    let f = Ldlt::factor(&op, 1e-13).map_err(EigenError::Singular)?;
    let apply = |q: &[f64]| f.solve(&m.apply(q));
    let max_steps = n;
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut mq: Vec<Vec<f64>> = Vec::new();
    let mut alpha: Vec<f64> = Vec::new();
    let mut beta: Vec<f64> = Vec::new();
    let orth = |r: &mut Vec<f64>, q: &[Vec<f64>], mq: &[Vec<f64>]| {
        for _ in 0..2 {
            for (qi, mqi) in q.iter().zip(mq) {
                let c = dot(mqi, r);
                r.iter_mut().zip(qi).for_each(|(a, b)| *a -= c * b);
            }
        }
    };
    let mut seed = 0;
    let mut r = apply(&start_vector(n, seed));
    let mut best: Option<EigenPairs> = None;
    let mut next_check = count.max(1);
    let mut beta_max: f64 = 0.0;
    while q.len() < max_steps {
        let mr = m.apply(&r);
        let mut b = dot(&r, &mr).max(0.0).sqrt();
        beta_max = beta_max.max(b);
        let mut restart = b.is_nan() || b <= 1e-10 * beta_max || b == 0.0;
        let mut tries = 0;
        let mut rr = std::mem::take(&mut r);
        let mut mrr = mr;
        while restart {
            // the Krylov space is invariant: a new direction, orthogonal
            // to the basis, starts the next block (β = 0 between them)
            tries += 1;
            seed += 1;
            if tries > 3 {
                break;
            }
            let mut v = apply(&start_vector(n, seed));
            orth(&mut v, &q, &mq);
            let mv = m.apply(&v);
            let bv = dot(&v, &mv).max(0.0).sqrt();
            if bv > 1e-8 * beta_max.max(f64::MIN_POSITIVE) {
                rr = v;
                mrr = mv;
                b = bv;
                restart = false;
                if !q.is_empty() {
                    beta.push(0.0);
                }
            }
        }
        if restart {
            break; // the space is exhausted (massless unknowns)
        }
        if !q.is_empty() && beta.len() < q.len() {
            beta.push(b);
        }
        let qj: Vec<f64> = rr.iter().map(|v| v / b).collect();
        let mqj: Vec<f64> = mrr.iter().map(|v| v / b).collect();
        let mut u = apply(&qj);
        let a = dot(&mqj, &u);
        if let (Some(qp), Some(&bp)) = (q.last(), beta.last()) {
            u.iter_mut().zip(qp).for_each(|(x, y)| *x -= bp * y);
        }
        u.iter_mut().zip(&qj).for_each(|(x, y)| *x -= a * y);
        q.push(qj);
        mq.push(mqj);
        alpha.push(a);
        orth(&mut u, &q, &mq);
        r = u;
        let j = q.len();
        if j >= next_check || j == max_steps {
            let mr = m.apply(&r);
            let bj = dot(&r, &mr).max(0.0).sqrt();
            let pairs = ritz(&q, &alpha, &beta, bj, count, shift, tol)?;
            let done = pairs.converged;
            best = Some(pairs);
            if done {
                break;
            }
            next_check = j + 4.max(count / 4);
        }
    }
    let mut pairs = match best {
        Some(p) if p.steps == q.len() => p,
        _ => {
            let bj = if r.len() == n {
                dot(&r, &m.apply(&r)).max(0.0).sqrt()
            } else {
                0.0 // the Krylov space was exhausted: the Ritz pairs are exact
            };
            ritz(&q, &alpha, &beta, bj, count, shift, tol)?
        }
    };
    pairs.steps = q.len();
    Ok(pairs)
}

fn ritz(
    q: &[Vec<f64>],
    alpha: &[f64],
    beta: &[f64],
    beta_last: f64,
    count: usize,
    shift: f64,
    tol: f64,
) -> Result<EigenPairs, EigenError> {
    let j = alpha.len();
    let (theta, s) =
        tridiagonal(alpha, &beta[..j.saturating_sub(1)]).ok_or(EigenError::Tridiagonal)?;
    let mut idx: Vec<usize> = (0..j).collect();
    idx.sort_by(|&a, &b| theta[b].abs().total_cmp(&theta[a].abs()));
    idx.truncate(count);
    let n = q.first().map_or(0, Vec::len);
    let mut out: Vec<(f64, Vec<f64>, f64)> = Vec::new();
    let mut converged = idx.len() == count;
    for &i in &idx {
        let th = theta[i];
        let res = (beta_last * s[(j - 1) * j + i]).abs() / th.abs().max(f64::MIN_POSITIVE);
        if res.is_nan() || res > tol {
            converged = false;
        }
        let mut v = vec![0.0; n];
        for (l, ql) in q.iter().enumerate() {
            let c = s[l * j + i];
            v.iter_mut().zip(ql).for_each(|(a, b)| *a += c * b);
        }
        out.push((shift + 1.0 / th, v, res));
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(EigenPairs {
        values: out.iter().map(|o| o.0).collect(),
        residuals: out.iter().map(|o| o.2).collect(),
        vectors: out.into_iter().map(|o| o.1).collect(),
        steps: j,
        converged,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Assembly;

    #[test]
    fn tridiagonal_of_the_path() {
        // eigenvalues of tridiag(-1, 2, -1) of order n: 2 − 2 cos(kπ/(n+1))
        let n = 30;
        let (mut v, z) = tridiagonal(&vec![2.0; n], &vec![-1.0; n - 1]).unwrap();
        // eigenvector check before sorting
        for c in 0..n {
            let x: Vec<f64> = (0..n).map(|r| z[r * n + c]).collect();
            for r in 0..n {
                let mut t = 2.0 * x[r];
                if r > 0 {
                    t -= x[r - 1];
                }
                if r + 1 < n {
                    t -= x[r + 1];
                }
                assert!((t - v[c] * x[r]).abs() < 1e-12);
            }
        }
        v.sort_by(f64::total_cmp);
        for (k, l) in v.iter().enumerate() {
            let e = 2.0 - 2.0 * (((k + 1) as f64) * std::f64::consts::PI / (n as f64 + 1.0)).cos();
            assert!((l - e).abs() < 1e-12, "{k}: {l} {e}");
        }
    }

    #[test]
    fn lanczos_finds_the_lowest_of_a_chain() {
        // spring-mass chain fixed at one end: K = tridiag, M = I·m;
        // λ_k = (k/m)·4 sin²((2k−1)π/(2(2n+1)))
        let n = 200;
        let (ks, ms) = (3.0, 2.0);
        let mut k = Assembly::new(n);
        let mut m = Assembly::new(n);
        for i in 0..n {
            k.add(i, i, if i + 1 < n { 2.0 * ks } else { ks });
            if i + 1 < n {
                k.add(i, i + 1, -ks);
            }
            m.add(i, i, ms);
        }
        let p = lanczos(&k.finish(), &m.finish(), 8, 0.0, 1e-10).unwrap();
        assert!(p.converged);
        for (i, l) in p.values.iter().enumerate() {
            let a = ((2 * i + 1) as f64) * std::f64::consts::PI / (2.0 * (2 * n + 1) as f64);
            let e = ks / ms * 4.0 * a.sin().powi(2);
            assert!((l / e - 1.0).abs() < 1e-9, "{i}: {l} {e}");
        }
    }
}
