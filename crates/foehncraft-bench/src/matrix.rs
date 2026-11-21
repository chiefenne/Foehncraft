//! Sparse matrices for the SpMV benchmark (D37): CSR and SELL-C-σ.
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use crate::mesh::FaceMesh;
use foehncraft_core::Real;
use foehncraft_kernels::unpack_face;

/// Compressed sparse row matrix.
#[derive(Debug, Clone)]
pub struct Csr<T> {
    /// Number of rows (and columns).
    pub n: usize,
    /// Row offsets, length `n + 1`.
    pub row_ptr: Vec<u32>,
    /// Column indices, ascending within a row.
    pub col: Vec<u32>,
    /// Values.
    pub val: Vec<T>,
}

impl<T: Real> Csr<T> {
    /// Matrix of a diffusion operator on the face graph of `mesh`:
    /// a_PN = −γ_f for each face, a_PP = Σ γ_f + `shift`. Symmetric and, for
    /// `shift > 0`, positive definite — the shape of a pressure matrix.
    pub fn from_mesh(mesh: &FaceMesh, gamma: &[T], shift: T) -> Self {
        let n = mesh.n_cells;
        let mut row_ptr = Vec::with_capacity(n + 1);
        let mut col = Vec::with_capacity(n + 2 * mesh.n_faces());
        let mut val = Vec::with_capacity(n + 2 * mesh.n_faces());
        row_ptr.push(0u32);
        let mut row: Vec<(u32, T)> = Vec::new();
        for c in 0..n {
            row.clear();
            let mut diag = shift;
            for k in mesh.cell_ptr[c] as usize..mesh.cell_ptr[c + 1] as usize {
                let (f, is_nb) = unpack_face(mesh.cell_faces[k]);
                let other = if is_nb {
                    mesh.owner[f]
                } else {
                    mesh.neighbour[f]
                };
                diag += gamma[f];
                row.push((other, -gamma[f]));
            }
            row.push((c as u32, diag));
            row.sort_unstable_by_key(|e| e.0);
            for &(j, v) in &row {
                col.push(j);
                val.push(v);
            }
            row_ptr.push(col.len() as u32);
        }
        let m = Self {
            n,
            row_ptr,
            col,
            val,
        };
        m.validate().expect("assembled CSR matrix is valid");
        m
    }

    /// Number of stored entries.
    pub fn nnz(&self) -> usize {
        self.col.len()
    }

    /// Checks the preconditions of `foehncraft_kernels::spmv::csr_row`.
    pub fn validate(&self) -> Result<(), String> {
        if self.row_ptr.len() != self.n + 1 || self.row_ptr[0] != 0 {
            return Err("row_ptr has the wrong length or start".into());
        }
        if self.row_ptr.windows(2).any(|w| w[0] > w[1])
            || *self.row_ptr.last().unwrap() as usize != self.col.len()
            || self.col.len() != self.val.len()
        {
            return Err("row_ptr is not a valid offset array".into());
        }
        if self.col.iter().any(|&j| j as usize >= self.n) {
            return Err("column index out of range".into());
        }
        Ok(())
    }
}

/// SELL-C-σ matrix (Kreutzer et al. 2014).
#[derive(Debug, Clone)]
pub struct SellCs<T> {
    /// Number of rows.
    pub n: usize,
    /// Slice height C.
    pub chunk: usize,
    /// Sorting window σ (1: no sorting).
    pub sigma: usize,
    /// Start of each slice in `col` and `val`.
    pub slice_ptr: Vec<u32>,
    /// Length (longest row) of each slice.
    pub slice_len: Vec<u32>,
    /// Column indices, column-major within a slice; padding repeats a valid
    /// column with value zero.
    pub col: Vec<u32>,
    /// Values.
    pub val: Vec<T>,
    /// Original row of each sorted row.
    pub perm: Vec<u32>,
}

impl<T: Real> SellCs<T> {
    /// Converts a CSR matrix.
    pub fn from_csr(a: &Csr<T>, chunk: usize, sigma: usize) -> Self {
        assert!(chunk >= 1 && sigma >= 1);
        let n = a.n;
        let len = |r: usize| (a.row_ptr[r + 1] - a.row_ptr[r]) as usize;
        let mut perm: Vec<u32> = (0..n as u32).collect();
        if sigma > 1 {
            for w in perm.chunks_mut(sigma) {
                // stable: equal lengths keep their order
                w.sort_by_key(|&r| std::cmp::Reverse(len(r as usize)));
            }
        }
        let n_slices = n.div_ceil(chunk);
        let mut slice_ptr = Vec::with_capacity(n_slices);
        let mut slice_len = Vec::with_capacity(n_slices);
        let mut col = Vec::new();
        let mut val = Vec::new();
        for s in 0..n_slices {
            let rows = s * chunk..((s + 1) * chunk).min(n);
            let width = rows
                .clone()
                .map(|r| len(perm[r] as usize))
                .max()
                .unwrap_or(0);
            slice_ptr.push(col.len() as u32);
            slice_len.push(width as u32);
            let base = col.len();
            col.resize(base + width * chunk, 0);
            val.resize(base + width * chunk, T::ZERO);
            for lane in 0..chunk {
                let r = s * chunk + lane;
                // rows past the end of the matrix (last slice) are padding
                let orig = if r < n {
                    perm[r] as usize
                } else {
                    perm[n - 1] as usize
                };
                let lo = a.row_ptr[orig] as usize;
                let l = if r < n { len(orig) } else { 0 };
                for j in 0..width {
                    let k = base + j * chunk + lane;
                    if j < l {
                        col[k] = a.col[lo + j];
                        val[k] = a.val[lo + j];
                    } else {
                        col[k] = orig as u32;
                    }
                }
            }
        }
        let m = Self {
            n,
            chunk,
            sigma,
            slice_ptr,
            slice_len,
            col,
            val,
            perm,
        };
        m.validate().expect("converted SELL matrix is valid");
        m
    }

    /// Stored entries including padding.
    pub fn stored(&self) -> usize {
        self.col.len()
    }

    /// Checks the preconditions of `foehncraft_kernels::spmv::sell_row`.
    pub fn validate(&self) -> Result<(), String> {
        let n_slices = self.n.div_ceil(self.chunk);
        if self.slice_ptr.len() != n_slices || self.slice_len.len() != n_slices {
            return Err("slice arrays have the wrong length".into());
        }
        for s in 0..n_slices {
            let end = self.slice_ptr[s] as usize + self.slice_len[s] as usize * self.chunk;
            if end > self.col.len() {
                return Err(format!("slice {s} exceeds the storage"));
            }
        }
        if self.col.len() != self.val.len() || self.col.iter().any(|&j| j as usize >= self.n) {
            return Err("column index out of range".into());
        }
        if self.perm.len() != self.n {
            return Err("permutation has the wrong length".into());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::{FaceMesh, Ordering};

    #[test]
    fn sell_reproduces_csr_product() {
        let mesh = FaceMesh::periodic_box(5, 4, 3, Ordering::Shuffled, 3);
        let gamma: Vec<f64> = (0..mesh.n_faces()).map(|f| 1.0 + f as f64 * 1e-3).collect();
        let a = Csr::from_mesh(&mesh, &gamma, 0.5);
        assert_eq!(a.nnz(), 7 * mesh.n_cells);
        let x: Vec<f64> = (0..a.n).map(|i| (i as f64).sin()).collect();
        let mut y = vec![0.0; a.n];
        for (r, yr) in y.iter_mut().enumerate() {
            // SAFETY: r < n; `from_mesh` builds a valid CSR matrix.
            *yr = unsafe { foehncraft_kernels::spmv::csr_row(r, &a.row_ptr, &a.col, &a.val, &x) };
        }
        for (chunk, sigma) in [(1, 1), (4, 1), (32, 1), (8, 16)] {
            let s = SellCs::from_csr(&a, chunk, sigma);
            let mut ys = vec![0.0; a.n];
            for r in 0..a.n {
                // SAFETY: r < n; `from_csr` builds a valid SELL-C matrix.
                ys[s.perm[r] as usize] = unsafe {
                    foehncraft_kernels::spmv::sell_row(
                        r,
                        s.chunk,
                        &s.slice_ptr,
                        &s.slice_len,
                        &s.col,
                        &s.val,
                        &x,
                    )
                };
            }
            // same entries in the same order: identical bits
            assert_eq!(y, ys, "chunk {chunk}, sigma {sigma}");
        }
    }
}
