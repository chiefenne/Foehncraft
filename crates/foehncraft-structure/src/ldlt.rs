//! Sparse direct solver: A = P L D Lᵀ Pᵀ for symmetric matrices that need
//! not be positive definite (K − σM of a shifted eigenproblem).
//!
//! The ordering P is reverse Cuthill–McKee (Cuthill & McKee 1969, George
//! 1971), which keeps the fill inside a small envelope for the slender
//! meshes of reeds and beams. The factorisation is the up-looking
//! algorithm of Davis (*Direct Methods for Sparse Linear Systems*, SIAM
//! 2006, ch. 4; Davis 2005, ACM TOMS 31): row k of L comes from a sparse
//! triangular solve whose pattern is the reach of row k of A in the
//! elimination tree, so only the true fill is stored. No pivoting: the
//! pivots D are the inertia of A (Sylvester's law), which gives the
//! Sturm count of the eigenproblem (how many eigenvalues lie below the
//! shift).

use crate::sparse::SymCsc;

/// A zero pivot: the matrix is singular (to the relative tolerance) at
/// this unknown of the original numbering.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZeroPivot {
    /// Unknown (original numbering).
    pub index: usize,
    /// The pivot found.
    pub pivot: f64,
}

/// The factors.
#[derive(Debug, Clone)]
pub struct Ldlt {
    n: usize,
    /// `perm[k]`: original index of the k-th unknown of the ordering.
    perm: Vec<usize>,
    /// Columns of L (strictly lower), in the ordering.
    lp: Vec<usize>,
    li: Vec<usize>,
    lx: Vec<f64>,
    d: Vec<f64>,
}

/// Reverse Cuthill–McKee ordering of a symmetric pattern: each connected
/// part from a pseudo-peripheral start (the last level of a breadth-first
/// search from the lowest-degree node, repeated while the depth grows).
pub fn rcm(adj: &[Vec<usize>]) -> Vec<usize> {
    let n = adj.len();
    let mut order = Vec::with_capacity(n);
    let mut placed = vec![false; n];
    let deg = |i: usize| adj[i].len();
    let levels = |s: usize, mark: &mut Vec<usize>, stamp: usize| -> (usize, usize) {
        // returns (depth, a lowest-degree node of the last level)
        let mut front = vec![s];
        mark[s] = stamp;
        let mut depth = 0;
        let mut last = s;
        while !front.is_empty() {
            last = *front.iter().min_by_key(|&&i| (deg(i), i)).unwrap_or(&s);
            let mut next = Vec::new();
            for &i in &front {
                for &j in &adj[i] {
                    if mark[j] != stamp {
                        mark[j] = stamp;
                        next.push(j);
                    }
                }
            }
            if !next.is_empty() {
                depth += 1;
            }
            front = next;
        }
        (depth, last)
    };
    let mut mark = vec![usize::MAX; n];
    let mut stamp = 0;
    for root in 0..n {
        if placed[root] {
            continue;
        }
        // lowest degree node of this part
        let mut part = vec![root];
        let mut seen = vec![root];
        placed[root] = true;
        while let Some(i) = part.pop() {
            for &j in &adj[i] {
                if !placed[j] {
                    placed[j] = true;
                    part.push(j);
                    seen.push(j);
                }
            }
        }
        let mut start = *seen.iter().min_by_key(|&&i| (deg(i), i)).unwrap_or(&root);
        stamp += 1;
        let (mut depth, mut cand) = levels(start, &mut mark, stamp);
        for _ in 0..8 {
            stamp += 1;
            let (d2, c2) = levels(cand, &mut mark, stamp);
            if d2 > depth {
                start = cand;
                depth = d2;
                cand = c2;
            } else {
                break;
            }
        }
        // Cuthill–McKee from `start`, neighbours by increasing degree
        stamp += 1;
        let first = order.len();
        order.push(start);
        mark[start] = stamp;
        let mut head = first;
        while head < order.len() {
            let i = order[head];
            head += 1;
            let mut nb: Vec<usize> = adj[i]
                .iter()
                .copied()
                .filter(|&j| mark[j] != stamp)
                .collect();
            nb.sort_by_key(|&j| (deg(j), j));
            for j in nb {
                mark[j] = stamp;
                order.push(j);
            }
        }
        order[first..].reverse();
    }
    order
}

impl Ldlt {
    /// Factorises `a` with the RCM ordering of its pattern. A pivot below
    /// `rel_tol` times the diagonal of A at that unknown is a zero pivot.
    pub fn factor(a: &SymCsc, rel_tol: f64) -> Result<Self, ZeroPivot> {
        let perm = rcm(&a.adjacency());
        Self::factor_with(a, perm, rel_tol)
    }

    /// Factorises with a given ordering (`perm[k]` is the original index
    /// of the k-th unknown).
    pub fn factor_with(a: &SymCsc, perm: Vec<usize>, rel_tol: f64) -> Result<Self, ZeroPivot> {
        let n = a.n;
        let mut pinv = vec![0usize; n];
        for (k, &i) in perm.iter().enumerate() {
            pinv[i] = k;
        }
        // permuted upper triangle C = P A Pᵀ by columns
        let mut cnt = vec![0usize; n + 1];
        for j in 0..n {
            for p in a.colptr[j]..a.colptr[j + 1] {
                let (pi, pj) = (pinv[a.rowidx[p]], pinv[j]);
                cnt[pi.max(pj) + 1] += 1;
            }
        }
        for k in 0..n {
            cnt[k + 1] += cnt[k];
        }
        let cp = cnt.clone();
        let mut next = cnt;
        let nnz = cp[n];
        let mut ci = vec![0usize; nnz];
        let mut cx = vec![0.0f64; nnz];
        for j in 0..n {
            for p in a.colptr[j]..a.colptr[j + 1] {
                let (pi, pj) = (pinv[a.rowidx[p]], pinv[j]);
                let c = pi.max(pj);
                let q = next[c];
                next[c] += 1;
                ci[q] = pi.min(pj);
                cx[q] = a.vals[p];
            }
        }
        let diag_of = |k: usize| -> f64 {
            let mut s = 0.0;
            for p in cp[k]..cp[k + 1] {
                if ci[p] == k {
                    s += cx[p];
                }
            }
            s
        };
        // symbolic: elimination tree and column counts of L
        let mut parent = vec![usize::MAX; n];
        let mut flag = vec![usize::MAX; n];
        let mut lnz = vec![0usize; n];
        for k in 0..n {
            flag[k] = k;
            for &i0 in &ci[cp[k]..cp[k + 1]] {
                let mut i = i0;
                if i < k {
                    while flag[i] != k {
                        if parent[i] == usize::MAX {
                            parent[i] = k;
                        }
                        lnz[i] += 1;
                        flag[i] = k;
                        i = parent[i];
                    }
                }
            }
        }
        let mut lp = vec![0usize; n + 1];
        for k in 0..n {
            lp[k + 1] = lp[k] + lnz[k];
        }
        let mut li = vec![0usize; lp[n]];
        let mut lx = vec![0.0f64; lp[n]];
        let mut d = vec![0.0f64; n];
        // numeric
        let mut y = vec![0.0f64; n];
        let mut pattern = vec![0usize; n];
        let mut fill = vec![0usize; n];
        for k in 0..n {
            let mut top = n;
            flag[k] = k;
            for p in cp[k]..cp[k + 1] {
                let mut i = ci[p];
                y[i] += cx[p];
                let mut len = 0;
                while flag[i] != k {
                    pattern[len] = i;
                    len += 1;
                    flag[i] = k;
                    i = parent[i];
                }
                while len > 0 {
                    top -= 1;
                    len -= 1;
                    pattern[top] = pattern[len];
                }
            }
            d[k] = y[k];
            y[k] = 0.0;
            for &i in &pattern[top..n] {
                let yi = y[i];
                y[i] = 0.0;
                let end = lp[i] + fill[i];
                for p in lp[i]..end {
                    y[li[p]] -= lx[p] * yi;
                }
                let l = yi / d[i];
                d[k] -= l * yi;
                li[end] = k;
                lx[end] = l;
                fill[i] += 1;
            }
            let scale = diag_of(k).abs();
            if !d[k].is_finite() || d[k].abs() <= rel_tol * scale {
                return Err(ZeroPivot {
                    index: perm[k],
                    pivot: d[k],
                });
            }
        }
        Ok(Self {
            n,
            perm,
            lp,
            li,
            lx,
            d,
        })
    }

    /// Order.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Entries of L (the fill).
    pub fn nnz(&self) -> usize {
        self.lp[self.n]
    }

    /// Number of negative pivots: the eigenvalues of A below zero.
    pub fn negative_pivots(&self) -> usize {
        self.d.iter().filter(|&&v| v < 0.0).count()
    }

    /// Solves `A x = b`.
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x: Vec<f64> = (0..n).map(|k| b[self.perm[k]]).collect();
        for j in 0..n {
            let xj = x[j];
            for p in self.lp[j]..self.lp[j + 1] {
                x[self.li[p]] -= self.lx[p] * xj;
            }
        }
        x.iter_mut().zip(&self.d).for_each(|(v, d)| *v /= d);
        for j in (0..n).rev() {
            let mut s = x[j];
            for p in self.lp[j]..self.lp[j + 1] {
                s -= self.lx[p] * x[self.li[p]];
            }
            x[j] = s;
        }
        let mut out = vec![0.0; n];
        for k in 0..n {
            out[self.perm[k]] = x[k];
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Assembly;

    fn laplace(n: usize, shift: f64) -> SymCsc {
        let mut a = Assembly::new(n);
        for i in 0..n {
            a.add(i, i, 2.0 - shift);
            if i + 1 < n {
                a.add(i, i + 1, -1.0);
            }
        }
        // a long-range coupling, so the ordering matters
        a.add(0, n - 1, -0.5);
        a.add(0, 0, 0.5);
        a.add(n - 1, n - 1, 0.5);
        a.finish()
    }

    #[test]
    fn solves_and_counts_inertia() {
        let a = laplace(50, 0.0);
        let f = Ldlt::factor(&a, 1e-13).unwrap();
        assert_eq!(f.negative_pivots(), 0);
        let xs: Vec<f64> = (0..50).map(|i| (i as f64 * 0.37).sin()).collect();
        let b = a.apply(&xs);
        let x = f.solve(&b);
        for (u, v) in x.iter().zip(&xs) {
            assert!((u - v).abs() < 1e-10, "{u} {v}");
        }
        // inertia: the path Laplacian tridiag(−1, 2, −1) has the eigenvalues
        // 2 − 2 cos(kπ/(n+1)); below the shift σ = 1.3 lie exactly those
        // with cos(kπ/(n+1)) > 0.35
        let n = 50;
        let mut p = Assembly::new(n);
        for i in 0..n {
            p.add(i, i, 2.0 - 1.3);
            if i + 1 < n {
                p.add(i, i + 1, -1.0);
            }
        }
        let s = p.finish();
        let g = Ldlt::factor(&s, 1e-13).unwrap();
        let exact = (1..=n)
            .filter(|&k| {
                2.0 - 2.0 * (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos() < 1.3
            })
            .count();
        assert_eq!(g.negative_pivots(), exact);
        let b = s.apply(&xs);
        let x = g.solve(&b);
        for (u, v) in x.iter().zip(&xs) {
            assert!((u - v).abs() < 1e-8, "{u} {v}");
        }
    }

    #[test]
    fn singular_is_found() {
        let mut a = Assembly::new(3);
        a.add(0, 0, 1.0);
        a.add(1, 1, 1.0);
        a.add(0, 1, -1.0);
        a.add(2, 2, 1.0);
        let e = Ldlt::factor(&a.finish(), 1e-12).unwrap_err();
        assert!(e.index < 2);
    }
}
