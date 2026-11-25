//! Symmetric sparse matrices: assembly from element blocks and the
//! product with a vector.
//!
//! The matrices of the structure are symmetric; [`Assembly`] collects the
//! entries of the upper triangle as triplets (duplicates are summed in a
//! fixed order, so the result does not depend on how the triplets were
//! grouped), and [`SymCsc`] stores the upper triangle by columns: column
//! `j` holds the rows `i ≤ j` in increasing order.

/// Triplets of a symmetric matrix (upper triangle).
#[derive(Debug, Clone, Default)]
pub struct Assembly {
    n: usize,
    rows: Vec<usize>,
    cols: Vec<usize>,
    vals: Vec<f64>,
}

impl Assembly {
    /// An empty `n × n` matrix.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            ..Default::default()
        }
    }

    /// Order.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Adds `v` to entry `(i, j)` and, by symmetry, to `(j, i)`: one call
    /// per pair (do not add both).
    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        if v == 0.0 {
            return;
        }
        let (r, c) = if i <= j { (i, j) } else { (j, i) };
        self.rows.push(r);
        self.cols.push(c);
        self.vals.push(v);
    }

    /// Adds a dense symmetric block `b` (`dofs.len()²`, row major) at the
    /// global indices `dofs`; `None` drops a row and column (a fixed DOF).
    pub fn add_block(&mut self, dofs: &[Option<usize>], b: &[f64]) {
        let m = dofs.len();
        for a in 0..m {
            let Some(i) = dofs[a] else { continue };
            for c in 0..m {
                let Some(j) = dofs[c] else { continue };
                if a < c {
                    // the pair (a, c) stands for (c, a) too; two local
                    // indices on one global index meet on its diagonal
                    let v = b[a * m + c];
                    self.add(i, j, if i == j { 2.0 * v } else { v });
                } else if a == c {
                    self.add(i, j, b[a * m + c]);
                }
            }
        }
    }

    /// The compressed matrix.
    pub fn finish(&self) -> SymCsc {
        let n = self.n;
        let mut order: Vec<usize> = (0..self.vals.len()).collect();
        order.sort_by_key(|&k| (self.cols[k], self.rows[k]));
        let mut colptr = vec![0usize; n + 1];
        let mut rowidx = Vec::new();
        let mut vals: Vec<f64> = Vec::new();
        let mut last: Option<(usize, usize)> = None;
        for &k in &order {
            let key = (self.cols[k], self.rows[k]);
            if last == Some(key) {
                if let Some(v) = vals.last_mut() {
                    *v += self.vals[k];
                }
            } else {
                rowidx.push(key.1);
                vals.push(self.vals[k]);
                colptr[key.0 + 1] += 1;
                last = Some(key);
            }
        }
        for j in 0..n {
            colptr[j + 1] += colptr[j];
        }
        SymCsc {
            n,
            colptr,
            rowidx,
            vals,
        }
    }
}

/// Upper triangle of a symmetric matrix by columns.
#[derive(Debug, Clone)]
pub struct SymCsc {
    /// Order.
    pub n: usize,
    /// Column starts (`n + 1`).
    pub colptr: Vec<usize>,
    /// Row of each entry (`i ≤ j`, increasing within a column).
    pub rowidx: Vec<usize>,
    /// Values.
    pub vals: Vec<f64>,
}

impl SymCsc {
    /// `y = A x`.
    pub fn mul(&self, x: &[f64], y: &mut [f64]) {
        y.iter_mut().for_each(|v| *v = 0.0);
        for j in 0..self.n {
            for p in self.colptr[j]..self.colptr[j + 1] {
                let i = self.rowidx[p];
                let a = self.vals[p];
                y[i] += a * x[j];
                if i != j {
                    y[j] += a * x[i];
                }
            }
        }
    }

    /// `A x` as a new vector.
    pub fn apply(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.n];
        self.mul(x, &mut y);
        y
    }

    /// Diagonal entries.
    pub fn diagonal(&self) -> Vec<f64> {
        let mut d = vec![0.0; self.n];
        for (j, dj) in d.iter_mut().enumerate() {
            let e = self.colptr[j + 1];
            if e > self.colptr[j] && self.rowidx[e - 1] == j {
                *dj = self.vals[e - 1];
            }
        }
        d
    }

    /// `self + s · other`, on the union of the two patterns.
    pub fn plus_scaled(&self, s: f64, other: &SymCsc) -> SymCsc {
        let mut a = Assembly::new(self.n);
        for m in [(1.0, self), (s, other)] {
            for j in 0..m.1.n {
                for p in m.1.colptr[j]..m.1.colptr[j + 1] {
                    a.rows.push(m.1.rowidx[p]);
                    a.cols.push(j);
                    a.vals.push(m.0 * m.1.vals[p]);
                }
            }
        }
        a.finish()
    }

    /// `xᵀ A y`.
    pub fn form(&self, x: &[f64], y: &[f64]) -> f64 {
        let ay = self.apply(y);
        dot(x, &ay)
    }

    /// The adjacency of the pattern (without the diagonal), symmetric.
    pub fn adjacency(&self) -> Vec<Vec<usize>> {
        let mut adj = vec![Vec::new(); self.n];
        for j in 0..self.n {
            for p in self.colptr[j]..self.colptr[j + 1] {
                let i = self.rowidx[p];
                if i != j {
                    adj[i].push(j);
                    adj[j].push(i);
                }
            }
        }
        for a in &mut adj {
            a.sort_unstable();
            a.dedup();
        }
        adj
    }
}

/// `xᵀ y`, summed in index order.
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}
