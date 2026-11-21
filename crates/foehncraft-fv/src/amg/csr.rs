//! Host-side sparse matrices of the multigrid setup: compressed sparse rows
//! in float64 for the setup algebra, and the sliced-ELL form in which the
//! backends store and apply the hierarchy.
//!
//! Every product and transpose works row by row with a fixed order of the
//! entries (columns ascending), so the setup is deterministic and does not
//! depend on the number of threads (D31).
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use crate::disc::CHUNK;
use foehncraft_core::Real;
use rayon::prelude::*;

/// A sparse matrix in compressed-sparse-row form, columns ascending in each
/// row.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Csr {
    /// Rows.
    pub n_rows: usize,
    /// Columns.
    pub n_cols: usize,
    /// Start of each row (length `n_rows + 1`).
    pub row_ptr: Vec<usize>,
    /// Column of each entry.
    pub col: Vec<u32>,
    /// Value of each entry.
    pub val: Vec<f64>,
}

impl Csr {
    /// Builds a matrix from rows of (column, value) pairs; each row is sorted
    /// by column and duplicate columns are summed in their given order.
    pub fn from_rows(n_cols: usize, mut rows: Vec<Vec<(u32, f64)>>) -> Self {
        let n_rows = rows.len();
        // sort (stable: duplicates keep their order) and merge each row
        rows.par_iter_mut().with_min_len(4096).for_each(|r| {
            r.sort_by_key(|e| e.0);
            let mut w = 0;
            for k in 0..r.len() {
                if w > 0 && r[w - 1].0 == r[k].0 {
                    r[w - 1].1 += r[k].1;
                } else {
                    r[w] = r[k];
                    w += 1;
                }
            }
            r.truncate(w);
        });
        let nnz = rows.iter().map(Vec::len).sum();
        let mut row_ptr = Vec::with_capacity(n_rows + 1);
        row_ptr.push(0);
        let mut col = Vec::with_capacity(nnz);
        let mut val = Vec::with_capacity(nnz);
        for r in rows {
            for (c, v) in r {
                col.push(c);
                val.push(v);
            }
            row_ptr.push(col.len());
        }
        Self {
            n_rows,
            n_cols,
            row_ptr,
            col,
            val,
        }
    }

    /// Number of stored entries.
    pub fn nnz(&self) -> usize {
        self.col.len()
    }

    /// Entries of row `i`.
    pub fn row(&self, i: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (lo, hi) = (self.row_ptr[i], self.row_ptr[i + 1]);
        self.col[lo..hi]
            .iter()
            .zip(&self.val[lo..hi])
            .map(|(&c, &v)| (c as usize, v))
    }

    /// Diagonal entries (zero where a row has none).
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.n_rows)
            .map(|i| self.row(i).find(|&(c, _)| c == i).map_or(0.0, |e| e.1))
            .collect()
    }

    /// y = A x (rows in parallel, each in column order).
    pub fn spmv(&self, x: &[f64]) -> Vec<f64> {
        (0..self.n_rows)
            .into_par_iter()
            .with_min_len(4096)
            .map(|i| self.row(i).map(|(c, v)| v * x[c]).sum())
            .collect()
    }

    /// The symmetric part ½(A + Aᵀ) (the pattern is the union of the two;
    /// every entry is the exact mean, so a symmetric matrix is returned
    /// bitwise unchanged). The multigrid of a non-symmetric transport
    /// matrix (convection) is built from it (stage 6 PR 6d, D36): the
    /// V-cycle then preconditions BiCGStab on the full matrix.
    pub fn symmetric_part(&self) -> Csr {
        let t = self.transpose();
        let rows = (0..self.n_rows)
            .map(|i| {
                let mut r = Vec::new();
                let (mut a, mut b) = (self.row(i).peekable(), t.row(i).peekable());
                loop {
                    match (a.peek().copied(), b.peek().copied()) {
                        (None, None) => break,
                        (Some((ca, va)), Some((cb, vb))) if ca == cb => {
                            r.push((ca as u32, 0.5 * va + 0.5 * vb));
                            a.next();
                            b.next();
                        }
                        (Some((ca, va)), Some((cb, _))) if ca < cb => {
                            r.push((ca as u32, 0.5 * va));
                            a.next();
                        }
                        (Some((ca, va)), None) => {
                            r.push((ca as u32, 0.5 * va));
                            a.next();
                        }
                        (_, Some((cb, vb))) => {
                            r.push((cb as u32, 0.5 * vb));
                            b.next();
                        }
                    }
                }
                r
            })
            .collect();
        Csr::from_rows(self.n_cols, rows)
    }

    /// The transpose; the entries of a column keep the row order.
    pub fn transpose(&self) -> Csr {
        let mut count = vec![0usize; self.n_cols + 1];
        for &c in &self.col {
            count[c as usize + 1] += 1;
        }
        for i in 0..self.n_cols {
            count[i + 1] += count[i];
        }
        let row_ptr = count.clone();
        let mut next = count;
        let mut col = vec![0u32; self.nnz()];
        let mut val = vec![0.0; self.nnz()];
        for i in 0..self.n_rows {
            for (c, v) in self.row(i) {
                let k = next[c];
                col[k] = i as u32;
                val[k] = v;
                next[c] += 1;
            }
        }
        Csr {
            n_rows: self.n_cols,
            n_cols: self.n_rows,
            row_ptr,
            col,
            val,
        }
    }

    /// The product self · b (Gustavson's row-by-row algorithm with a dense
    /// accumulator per row; the entries of a row are summed in the order of
    /// the row of `self`, then of the rows of `b`).
    pub fn matmul(&self, b: &Csr) -> Csr {
        assert_eq!(self.n_cols, b.n_rows, "matmul: inner dimensions");
        Csr::build_rows(
            self.n_rows,
            b.n_cols,
            || RowAccumulator::new(b.n_cols),
            |acc, i, col, val| {
                for (k, a) in self.row(i) {
                    for (j, v) in b.row(k) {
                        acc.add(j as u32, a * v);
                    }
                }
                acc.drain_sorted(col, val);
            },
        )
    }

    /// Builds a matrix row by row: `fill(scratch, i, col, val)` appends the
    /// entries of row `i`, columns ascending and each once, to `col` and
    /// `val`. Blocks of rows are built in parallel into flat arrays (one
    /// scratch per worker, no allocation per row) and joined in row order,
    /// so the result does not depend on the threads (D31).
    pub fn build_rows<S, I, F>(n_rows: usize, n_cols: usize, init: I, fill: F) -> Csr
    where
        I: Fn() -> S + Sync + Send,
        F: Fn(&mut S, usize, &mut Vec<u32>, &mut Vec<f64>) + Sync + Send,
    {
        const BLOCK: usize = 512;
        let blocks: Vec<(Vec<u32>, Vec<u32>, Vec<f64>)> = (0..n_rows.div_ceil(BLOCK))
            .into_par_iter()
            .map_init(init, |scratch, b| {
                let (lo, hi) = (b * BLOCK, ((b + 1) * BLOCK).min(n_rows));
                let mut len = Vec::with_capacity(hi - lo);
                let mut col = Vec::new();
                let mut val = Vec::new();
                for i in lo..hi {
                    let start = col.len();
                    fill(scratch, i, &mut col, &mut val);
                    len.push((col.len() - start) as u32);
                }
                (len, col, val)
            })
            .collect();
        let nnz = blocks.iter().map(|b| b.1.len()).sum();
        let mut out = Csr {
            n_rows,
            n_cols,
            row_ptr: Vec::with_capacity(n_rows + 1),
            col: Vec::with_capacity(nnz),
            val: Vec::with_capacity(nnz),
        };
        out.row_ptr.push(0);
        let mut end = 0usize;
        for (len, col, val) in blocks {
            for l in len {
                end += l as usize;
                out.row_ptr.push(end);
            }
            out.col.extend_from_slice(&col);
            out.val.extend_from_slice(&val);
        }
        out
    }

    /// Largest |a_ij − a_ji| relative to the largest |a_ij|.
    pub fn asymmetry(&self) -> f64 {
        if self.n_rows != self.n_cols {
            return f64::INFINITY;
        }
        let t = self.transpose();
        let scale = self.val.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        if scale == 0.0 {
            return 0.0;
        }
        let mut d = 0.0f64;
        for i in 0..self.n_rows {
            // both rows have ascending columns: merge
            let (mut a, mut b) = (self.row(i).peekable(), t.row(i).peekable());
            loop {
                match (a.peek().copied(), b.peek().copied()) {
                    (None, None) => break,
                    (Some((ca, va)), Some((cb, vb))) if ca == cb => {
                        d = d.max((va - vb).abs());
                        a.next();
                        b.next();
                    }
                    (Some((ca, va)), Some((cb, _))) if ca < cb => {
                        d = d.max(va.abs());
                        a.next();
                    }
                    (Some((_, va)), None) => {
                        d = d.max(va.abs());
                        a.next();
                    }
                    (_, Some((_, vb))) => {
                        d = d.max(vb.abs());
                        b.next();
                    }
                }
            }
        }
        d / scale
    }
}

/// The dense accumulator of one row of a sparse product or sum: entries
/// of the same column are added in the order they arrive (the sums of
/// [`Csr::from_rows`]), and the row leaves sorted by column.
pub struct RowAccumulator {
    pos: Vec<u32>,
    acc: Vec<(u32, f64)>,
}

impl RowAccumulator {
    /// An accumulator for rows of `n_cols` columns.
    pub fn new(n_cols: usize) -> Self {
        Self {
            pos: vec![u32::MAX; n_cols],
            acc: Vec::new(),
        }
    }

    /// Adds `v` to column `j` of the row.
    #[inline]
    pub fn add(&mut self, j: u32, v: f64) {
        let p = &mut self.pos[j as usize];
        if *p == u32::MAX {
            *p = self.acc.len() as u32;
            self.acc.push((j, v));
        } else {
            self.acc[*p as usize].1 += v;
        }
    }

    /// Appends the row, columns ascending, to `col` and `val` and empties
    /// the accumulator.
    pub fn drain_sorted(&mut self, col: &mut Vec<u32>, val: &mut Vec<f64>) {
        for e in &self.acc {
            self.pos[e.0 as usize] = u32::MAX;
        }
        // the columns are distinct: an unstable sort gives the same order
        self.acc.sort_unstable_by_key(|e| e.0);
        for &(j, v) in &self.acc {
            col.push(j);
            val.push(v);
        }
        self.acc.clear();
    }
}

/// A sliced-ELL matrix (SELL-C-σ, C = [`CHUNK`], σ = 1) owned by the host,
/// in the precision `T` of the hierarchy.
#[derive(Debug, Clone, Default)]
pub struct SellMatrix<T> {
    /// Rows.
    pub n_rows: usize,
    /// Columns.
    pub n_cols: usize,
    /// Slice height C.
    pub chunk: usize,
    /// Start of each slice.
    pub slice_ptr: Vec<u32>,
    /// Slots per row of each slice.
    pub slice_len: Vec<u32>,
    /// Column per slot (padding: column 0, value 0).
    pub col: Vec<u32>,
    /// Value per slot.
    pub val: Vec<T>,
    /// Real entries.
    pub nnz: usize,
}

impl<T: Real> SellMatrix<T> {
    /// Converts a CSR matrix (values rounded to `T`).
    pub fn from_csr(a: &Csr) -> Self {
        let n = a.n_rows;
        let n_slices = n.div_ceil(CHUNK).max(1);
        let mut slice_ptr = Vec::with_capacity(n_slices);
        let mut slice_len = Vec::with_capacity(n_slices);
        let mut total = 0usize;
        for s in 0..n_slices {
            let len = (s * CHUNK..((s + 1) * CHUNK).min(n))
                .map(|r| a.row_ptr[r + 1] - a.row_ptr[r])
                .max()
                .unwrap_or(0)
                .max(1);
            slice_ptr.push(total as u32);
            slice_len.push(len as u32);
            total += len * CHUNK;
        }
        let mut col = vec![0u32; total];
        let mut val = vec![T::ZERO; total];
        for r in 0..n {
            let s = r / CHUNK;
            for (j, (c, v)) in a.row(r).enumerate() {
                let k = slice_ptr[s] as usize + j * CHUNK + r % CHUNK;
                col[k] = c as u32;
                val[k] = T::from_f64(v);
            }
        }
        Self {
            n_rows: n,
            n_cols: a.n_cols,
            chunk: CHUNK,
            slice_ptr,
            slice_len,
            col,
            val,
            nnz: a.nnz(),
        }
    }

    /// Borrowed form for the kernel bodies.
    pub fn view(&self) -> foehncraft_kernels::amg::Sell<'_, T> {
        foehncraft_kernels::amg::Sell {
            chunk: self.chunk,
            slice_ptr: &self.slice_ptr,
            slice_len: &self.slice_len,
            col: &self.col,
            val: &self.val,
        }
    }

    /// Stored slots including padding.
    pub fn stored(&self) -> usize {
        self.col.len()
    }

    /// Bytes of the stored matrix (values and columns).
    pub fn bytes(&self) -> usize {
        self.stored() * (T::BYTES + 4) + 8 * self.slice_ptr.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lap1d(n: usize) -> Csr {
        let rows = (0..n)
            .map(|i| {
                let mut r = vec![(i as u32, 2.0)];
                if i > 0 {
                    r.push((i as u32 - 1, -1.0));
                }
                if i + 1 < n {
                    r.push((i as u32 + 1, -1.0));
                }
                r
            })
            .collect();
        Csr::from_rows(n, rows)
    }

    /// The product of two matrices with scattered, repeated columns equals
    /// the one of the reference (rows of pairs through `from_rows`), bit
    /// for bit, on any number of threads: `build_rows` joins blocks of rows
    /// in row order and `RowAccumulator` sums in arrival order.
    #[test]
    fn flat_products_equal_the_reference_bitwise() {
        let n = 1500; // more than two blocks of rows
        let row = |i: usize, shift: usize| -> Vec<(u32, f64)> {
            (0..1 + (i * 7 + shift) % 5)
                .map(|k| {
                    let c = (i * 31 + k * 97 + shift * 13) % n;
                    (c as u32, 0.1 + ((i + k * 3 + shift) % 11) as f64 * 0.37)
                })
                .collect()
        };
        let a = Csr::from_rows(n, (0..n).map(|i| row(i, 1)).collect());
        let b = Csr::from_rows(n, (0..n).map(|i| row(i, 2)).collect());
        let mut reference: Vec<Vec<(u32, f64)>> = Vec::new();
        for i in 0..n {
            let mut r = Vec::new();
            for (k, x) in a.row(i) {
                for (j, y) in b.row(k) {
                    r.push((j as u32, x * y));
                }
            }
            reference.push(r);
        }
        let reference = Csr::from_rows(n, reference);
        assert_eq!(a.matmul(&b), reference);
        for threads in [1, 3] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            assert_eq!(pool.install(|| a.matmul(&b)), reference);
        }
    }

    #[test]
    fn transpose_matmul_symmetry() {
        let a = lap1d(7);
        assert_eq!(a.transpose(), a);
        assert_eq!(a.asymmetry(), 0.0);
        let a2 = a.matmul(&a);
        assert_eq!(a2.row(3).collect::<Vec<_>>().len(), 5);
        let x: Vec<f64> = (0..7).map(|i| i as f64).collect();
        let y1 = a2.spmv(&x);
        let y2 = a.spmv(&a.spmv(&x));
        assert_eq!(y1, y2);
    }

    /// ½(A + Aᵀ) of a convection–diffusion matrix: upwind rows with an
    /// entry whose transpose is missing (one-sided coupling) keep half of
    /// it on both sides; the result is symmetric with the diagonal of A.
    #[test]
    fn symmetric_part_of_a_convection_matrix() {
        let rows = vec![
            vec![(0, 3.0), (1, -1.0)],
            vec![(0, -2.0), (1, 4.0), (2, -1.5)],
            vec![(2, 2.0)],
        ];
        let a = Csr::from_rows(3, rows);
        assert!(a.asymmetry() > 0.0);
        let s = a.symmetric_part();
        assert_eq!(s.asymmetry(), 0.0);
        assert_eq!(s.diagonal(), a.diagonal());
        let r1: Vec<_> = s.row(1).collect();
        assert_eq!(r1, vec![(0, -1.5), (1, 4.0), (2, -0.75)]);
        let r2: Vec<_> = s.row(2).collect();
        assert_eq!(r2, vec![(1, -0.75), (2, 2.0)]);
        // a symmetric matrix is its own symmetric part
        assert_eq!(lap1d(5).symmetric_part(), lap1d(5));
    }

    #[test]
    fn sell_round_trip() {
        let a = lap1d(70);
        let s = SellMatrix::<f64>::from_csr(&a);
        let x: Vec<f64> = (0..70).map(|i| (i as f64).sin()).collect();
        let y = a.spmv(&x);
        for (i, yi) in y.iter().enumerate() {
            // SAFETY: built by from_csr.
            let v = unsafe { s.view().row(i, &x) };
            assert_eq!(v, *yi);
        }
    }
}
