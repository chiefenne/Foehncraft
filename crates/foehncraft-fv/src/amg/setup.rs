//! Setup of the smoothed-aggregation hierarchy (see the module docs of
//! [`crate::amg`]). Everything here runs on the host in float64; the
//! backends round the result to the precision of the hierarchy when they
//! load it.

use super::csr::Csr;
use super::{Aggregation, AmgSettings, Precision, Smoother};
use anyhow::{Result, bail, ensure};
use rayon::prelude::*;
use serde::Serialize;
use std::time::Instant;

/// Marks a node without an aggregate.
const NONE: u32 = u32::MAX;

/// One level of the hierarchy (all but the coarsest).
#[derive(Debug, Clone)]
pub struct Level {
    /// The level operator A_l.
    pub a: Csr,
    /// 1 / diag(A_l).
    pub dinv: Vec<f64>,
    /// Prolongation P_l from level l + 1 to level l.
    pub p: Csr,
    /// Restriction R_l = P_lᵀ.
    pub r: Csr,
    /// Estimated λ_max(D⁻¹A_l).
    pub lambda_max: f64,
    /// Smoother steps (c₁, c₂): d ← c₁ d + c₂ D⁻¹(b − A x), x ← x + d.
    pub steps: Vec<(f64, f64)>,
}

/// The multigrid hierarchy on the host, in float64.
#[derive(Debug, Clone)]
pub struct Hierarchy {
    /// Levels 0 .. L−1 (level 0 is the fine matrix).
    pub levels: Vec<Level>,
    /// Size of the coarsest level L.
    pub coarse_n: usize,
    /// Dense inverse of the coarsest operator, row-major.
    pub coarse_inv: Vec<f64>,
    /// Precision in which the backends store and apply it.
    pub precision: Precision,
    /// Sizes and complexities.
    pub info: HierarchyInfo,
}

/// Sizes and complexities of a hierarchy, for the reports.
#[derive(Debug, Clone, Default, Serialize)]
pub struct HierarchyInfo {
    /// Rows per level, the coarsest last.
    pub rows: Vec<usize>,
    /// Non-zeros per level, the coarsest (dense) last.
    pub nnz: Vec<usize>,
    /// Σ rows / rows of the fine level.
    pub grid_complexity: f64,
    /// Σ nnz / nnz of the fine level (coarsest counted as its sparse
    /// Galerkin operator).
    pub operator_complexity: f64,
    /// Seconds of the host setup.
    pub setup_seconds: f64,
    /// Bytes of the hierarchy in its precision (matrices, diagonals, dense
    /// inverse).
    pub bytes: usize,
    /// Precision name.
    pub precision: String,
}

/// A deterministic start vector for the power iteration (a hash of the
/// index, so it is not orthogonal to the dominant eigenvector).
fn start_vector(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let mut z = (i as u64).wrapping_add(0x9E37_79B9_7F4A_7C15);
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            0.5 + (z >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect()
}

/// Estimate of λ_max(D⁻¹A) by `iters` power iterations with the Rayleigh
/// quotient (x, A x) / (x, D x) of the generalised problem A x = λ D x,
/// capped by the Gershgorin bound max_i Σ_j |a_ij| / a_ii (and the bound
/// itself when the quotient is not positive).
pub fn lambda_max(a: &Csr, diag: &[f64], iters: usize) -> f64 {
    let n = a.n_rows;
    let gersh = (0..n)
        .map(|i| a.row(i).map(|(_, v)| v.abs()).sum::<f64>() / diag[i].abs())
        .fold(0.0f64, f64::max);
    let mut x = start_vector(n);
    let mut lam = 0.0;
    for _ in 0..iters {
        let ax = a.spmv(&x);
        let xax = dot(&x, &ax);
        let xdx = {
            let mut s = foehncraft_core::exact::ExactSum::new();
            for (a, d) in x.iter().zip(diag) {
                s.add(a * a * d);
            }
            s.value()
        };
        lam = xax / xdx;
        let mut y: Vec<f64> = ax.iter().zip(diag).map(|(v, d)| v / d).collect();
        let nrm = y.iter().map(|v| v * v).sum::<f64>().sqrt();
        if nrm == 0.0 {
            break;
        }
        y.iter_mut().for_each(|v| *v /= nrm);
        x = y;
    }
    let lam = if lam > 0.0 && lam.is_finite() {
        lam.min(gersh)
    } else {
        gersh
    };
    lam.max(f64::MIN_POSITIVE)
}

/// Strong neighbours of every row: j ≠ i with |a_ij| ≥ θ √(|a_ii a_jj|).
fn strength(a: &Csr, diag: &[f64], theta: f64) -> Vec<Vec<(u32, f64)>> {
    (0..a.n_rows)
        .into_par_iter()
        .with_min_len(4096)
        .map(|i| {
            a.row(i)
                .filter(|&(j, v)| {
                    j != i && v != 0.0 && v.abs() >= theta * (diag[i] * diag[j]).abs().sqrt()
                })
                .map(|(j, v)| (j as u32, v.abs()))
                .collect()
        })
        .collect()
}

/// Three-phase greedy aggregation (Vaněk, Mandel, Brezina 1996) in row
/// order. Returns the aggregate of every node (`NONE` for isolated nodes)
/// and the number of aggregates.
pub fn aggregate(strong: &[Vec<(u32, f64)>]) -> (Vec<u32>, usize) {
    let n = strong.len();
    let mut agg = vec![NONE; n];
    let mut k = 0u32;
    // phase 1: whole untouched neighbourhoods
    for i in 0..n {
        if strong[i].is_empty() || agg[i] != NONE {
            continue;
        }
        if strong[i].iter().all(|&(j, _)| agg[j as usize] == NONE) {
            agg[i] = k;
            for &(j, _) in &strong[i] {
                agg[j as usize] = k;
            }
            k += 1;
        }
    }
    // phase 2: join the strongest neighbour aggregated in phase 1
    let phase1 = agg.clone();
    for i in 0..n {
        if agg[i] != NONE || strong[i].is_empty() {
            continue;
        }
        let mut best: Option<(f64, u32)> = None;
        for &(j, s) in &strong[i] {
            let a = phase1[j as usize];
            if a != NONE && best.is_none_or(|(bs, _)| s > bs) {
                best = Some((s, a));
            }
        }
        if let Some((_, a)) = best {
            agg[i] = a;
        }
    }
    // phase 3: the rest with their unaggregated strong neighbours
    for i in 0..n {
        if agg[i] != NONE || strong[i].is_empty() {
            continue;
        }
        agg[i] = k;
        for &(j, _) in &strong[i] {
            if agg[j as usize] == NONE {
                agg[j as usize] = k;
            }
        }
        k += 1;
    }
    (agg, k as usize)
}

/// Tentative prolongator: column k is the normalised indicator of
/// aggregate k (the constant near-null space of the pressure equation).
fn tentative(agg: &[u32], n_agg: usize) -> Csr {
    let mut size = vec![0usize; n_agg];
    for &a in agg {
        if a != NONE {
            size[a as usize] += 1;
        }
    }
    let rows = agg
        .iter()
        .map(|&a| {
            if a == NONE {
                Vec::new()
            } else {
                vec![(a, 1.0 / (size[a as usize] as f64).sqrt())]
            }
        })
        .collect();
    Csr::from_rows(n_agg, rows)
}

/// Filtered matrix A_F: strong off-diagonals kept, weak ones lumped onto the
/// diagonal (so that A_F keeps the row sums of A). A row whose lumped
/// diagonal would not be positive keeps its own diagonal instead: on a
/// coarse level the Galerkin product has positive off-diagonals, and a row
/// with a zero or negative sum (next to the large diagonal of a reference
/// cell, say) gave D_F⁻¹ = ∞ in the prolongator smoothing and a hierarchy of
/// NaN.
fn filtered(a: &Csr, strong: &[Vec<(u32, f64)>]) -> Csr {
    let rows = (0..a.n_rows)
        .into_par_iter()
        .with_min_len(4096)
        .map(|i| {
            let (mut own, mut lumped) = (0.0, 0.0);
            let mut r = Vec::new();
            let mut s = strong[i].iter().map(|e| e.0 as usize).peekable();
            for (j, v) in a.row(i) {
                while s.peek().is_some_and(|&sj| sj < j) {
                    s.next();
                }
                if j == i {
                    own += v;
                } else if s.peek() == Some(&j) {
                    r.push((j as u32, v));
                } else {
                    lumped += v;
                }
            }
            let diag = if own + lumped > 0.0 {
                own + lumped
            } else {
                own
            };
            r.push((i as u32, diag));
            r
        })
        .collect();
    Csr::from_rows(a.n_cols, rows)
}

/// Smoothed prolongator (I − ω D_F⁻¹ A_F) P_tent.
fn smooth_prolongator(af: &Csr, pt: &Csr, weight: f64, iters: usize) -> Csr {
    let d = af.diagonal();
    let lam = lambda_max(af, &d, iters);
    let omega = weight / lam;
    let ap = af.matmul(pt);
    let rows = (0..pt.n_rows)
        .into_par_iter()
        .with_min_len(4096)
        .map(|i| {
            let mut r: Vec<(u32, f64)> = pt.row(i).map(|(c, v)| (c as u32, v)).collect();
            let s = omega / d[i];
            r.extend(ap.row(i).map(|(c, v)| (c as u32, -s * v)));
            r
        })
        .collect();
    Csr::from_rows(pt.n_cols, rows)
}

/// Smoother coefficients of a level.
fn smoother_steps(s: &AmgSettings, lambda: f64) -> Vec<(f64, f64)> {
    match s.smoother {
        Smoother::Jacobi => vec![(0.0, s.jacobi_weight / lambda); s.degree],
        Smoother::Chebyshev => {
            let (lo, hi) = (
                s.chebyshev_interval[0] * lambda,
                s.chebyshev_interval[1] * lambda,
            );
            let theta = 0.5 * (hi + lo);
            let delta = 0.5 * (hi - lo);
            let sigma = theta / delta;
            let mut rho_prev = 1.0 / sigma;
            let mut steps = vec![(0.0, 1.0 / theta)];
            for _ in 1..s.degree {
                let rho = 1.0 / (2.0 * sigma - rho_prev);
                steps.push((rho * rho_prev, 2.0 * rho / delta));
                rho_prev = rho;
            }
            steps
        }
    }
}

/// Dense inverse of a symmetric positive definite matrix by its Cholesky
/// factorisation; fails when a pivot is not positive.
///
/// A = L Lᵀ, M = L⁻¹ and A⁻¹ = Mᵀ M, each of the three about n³/6
/// multiply-adds over contiguous rows of the lower triangle; the upper
/// triangle of the result is the mirror of the lower one (exactly
/// symmetric). The coarsest level is rebuilt with every setup, so this is
/// most of a setup on a small mesh, and of every per-step setup of a
/// sliding mesh by imprinting (D182).
pub fn dense_spd_inverse(a: &Csr) -> Result<Vec<f64>> {
    dense_spd_inverse_with(a, false)
}

/// [`dense_spd_inverse`] with the choice of the pivot test: `scaled` tests
/// each pivot against its own diagonal (independent of the scale of the
/// other rows), else against the largest diagonal.
pub fn dense_spd_inverse_with(a: &Csr, scaled: bool) -> Result<Vec<f64>> {
    let n = a.n_rows;
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for (j, v) in a.row(i) {
            if j <= i {
                l[i * n + j] = v;
            }
        }
    }
    let scale = (0..n).map(|i| l[i * n + i].abs()).fold(0.0f64, f64::max);
    // lower Cholesky factor in place, row by row: L[i][j] for j < i from
    // the dot product of the rows i and j, then the pivot
    for i in 0..n {
        let (done, rest) = l.split_at_mut(i * n);
        let row = &mut rest[..n];
        for j in 0..i {
            let rj = &done[j * n..j * n + j];
            let s = row[j] - dot(&row[..j], rj);
            row[j] = s / done[j * n + j];
        }
        // `scaled`: the pivot against its own diagonal, so that the test
        // does not depend on the scale of the other rows. The coarsest
        // level of an equation with a wide range of coefficients (the
        // elliptic equation of ζ-f on a fine mesh: 10¹⁴ between diagonals)
        // is positive definite although a pivot is far below the largest
        // diagonal.
        let d = row[i] - dot(&row[..i], &row[..i]);
        let floor = if scaled { row[i].abs() } else { scale };
        if d <= 1e-13 * floor || !d.is_finite() {
            bail!(
                "the coarsest multigrid level ({n} unknowns) is singular or not positive \
                 definite: the pressure level is not fixed (fix it with a fixed-value \
                 boundary or a reference cell, D30), or the matrix is not symmetric \
                 positive definite"
            );
        }
        row[i] = d.sqrt();
    }
    // M = L⁻¹, lower triangular, row by row: M[i] = (e_i − Σ_k<i L[i][k] M[k]) / L[i][i]
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        let (done, rest) = m.split_at_mut(i * n);
        let row = &mut rest[..n];
        let li = &l[i * n..i * n + i];
        for (k, &lik) in li.iter().enumerate() {
            if lik != 0.0 {
                axpy(-lik, &done[k * n..k * n + k + 1], &mut row[..=k]);
            }
        }
        let d = 1.0 / l[i * n + i];
        row[i] = 1.0;
        for v in &mut row[..=i] {
            *v *= d;
        }
    }
    // A⁻¹ = Mᵀ M: (A⁻¹)[i][j] = Σ_k M[k][i] M[k][j], lower triangle j ≤ i
    let mut inv = vec![0.0; n * n];
    for k in 0..n {
        let mk = &m[k * n..k * n + k + 1];
        for i in 0..=k {
            let c = mk[i];
            if c != 0.0 {
                axpy(c, &mk[..=i], &mut inv[i * n..i * n + i + 1]);
            }
        }
    }
    for i in 0..n {
        for j in i + 1..n {
            inv[i * n + j] = inv[j * n + i];
        }
    }
    Ok(inv)
}

/// Exact dot product (D93 (5): the eigenvalue estimate uses exact
/// reductions).
fn dot(a: &[f64], b: &[f64]) -> f64 {
    foehncraft_core::exact::exact_dot(a, b)
}

fn axpy(c: f64, x: &[f64], y: &mut [f64]) {
    for (y, x) in y.iter_mut().zip(x) {
        *y += c * x;
    }
}

/// Below this many rows the setup runs on one thread of its own: its
/// parallel loops are short, and handing each of them to the pool costs
/// more than the loop (the hierarchy is the same for any number of
/// threads, D31).
const SERIAL_ROWS: usize = 20_000;

/// Builds the hierarchy of the symmetric matrix `a`.
pub fn build(a: Csr, s: &AmgSettings) -> Result<Hierarchy> {
    thread_local! {
        // one helper thread per calling thread (parallel callers do not
        // queue behind each other)
        static ONE: Option<rayon::ThreadPool> = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .ok();
    }
    if a.n_rows < SERIAL_ROWS && rayon::current_num_threads() > 1 {
        return ONE.with(|one| match one {
            Some(pool) => pool.install(|| build_on(a, s)),
            None => build_on(a, s),
        });
    }
    build_on(a, s)
}

fn build_on(a: Csr, s: &AmgSettings) -> Result<Hierarchy> {
    s.validate()?;
    let t0 = Instant::now();
    ensure!(
        a.n_rows == a.n_cols && a.n_rows > 0,
        "AMG: square matrix needed"
    );
    let asym = a.asymmetry();
    ensure!(
        asym <= 1e-10,
        "AMG with CG needs a symmetric matrix (relative asymmetry {asym:.2e}); \
         use `linear_solver = \"bicgstab\"` for convection"
    );
    let mut levels = Vec::new();
    let mut cur = a;
    let mut theta = s.strength;
    while cur.n_rows > s.max_coarse && levels.len() + 1 < s.max_levels {
        let diag = cur.diagonal();
        ensure!(
            diag.iter().all(|&d| d > 0.0),
            "AMG: a diagonal entry is not positive on level {}",
            levels.len()
        );
        let strong = strength(&cur, &diag, theta);
        let (agg, n_agg) = aggregate(&strong);
        if n_agg == 0 || n_agg as f64 > 0.9 * cur.n_rows as f64 {
            break;
        }
        let pt = tentative(&agg, n_agg);
        let p = match s.aggregation {
            Aggregation::Plain => pt,
            Aggregation::Smoothed => smooth_prolongator(
                &filtered(&cur, &strong),
                &pt,
                s.prolongation_weight,
                s.power_iterations,
            ),
        };
        let r = p.transpose();
        let coarse = r.matmul(&cur.matmul(&p));
        let lambda = lambda_max(&cur, &diag, s.power_iterations);
        levels.push(Level {
            dinv: diag.iter().map(|d| 1.0 / d).collect(),
            steps: smoother_steps(s, lambda),
            lambda_max: lambda,
            a: cur,
            p,
            r,
        });
        cur = coarse;
        theta *= 0.5;
    }
    ensure!(
        cur.n_rows <= 4000,
        "AMG: coarsening stalled at {} unknowns (more than the dense limit of 4000); \
         lower `strength`",
        cur.n_rows
    );
    let coarse_inv = dense_spd_inverse_with(&cur, s.scaled_pivots)?;
    let tb = match s.precision {
        Precision::Mixed => 4,
        Precision::Double => 8,
    };
    let mut info = HierarchyInfo {
        precision: s.precision.name().into(),
        ..HierarchyInfo::default()
    };
    for l in &levels {
        info.rows.push(l.a.n_rows);
        info.nnz.push(l.a.nnz());
        info.bytes += (l.a.nnz() + l.p.nnz() + l.r.nnz()) * (tb + 4) + l.a.n_rows * tb;
    }
    info.rows.push(cur.n_rows);
    info.nnz.push(cur.nnz());
    info.bytes += cur.n_rows * cur.n_rows * tb;
    info.grid_complexity = info.rows.iter().sum::<usize>() as f64 / info.rows[0] as f64;
    info.operator_complexity = info.nnz.iter().sum::<usize>() as f64 / info.nnz[0] as f64;
    info.setup_seconds = t0.elapsed().as_secs_f64();
    Ok(Hierarchy {
        levels,
        coarse_n: cur.n_rows,
        coarse_inv,
        precision: s.precision,
        info,
    })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// 2D five-point Laplacian on an m × m grid with Dirichlet boundaries.
    pub(crate) fn lap2d(m: usize) -> Csr {
        let id = |i: usize, j: usize| (i * m + j) as u32;
        let rows = (0..m * m)
            .map(|k| {
                let (i, j) = (k / m, k % m);
                let mut r = vec![(k as u32, 4.0)];
                if i > 0 {
                    r.push((id(i - 1, j), -1.0));
                }
                if i + 1 < m {
                    r.push((id(i + 1, j), -1.0));
                }
                if j > 0 {
                    r.push((id(i, j - 1), -1.0));
                }
                if j + 1 < m {
                    r.push((id(i, j + 1), -1.0));
                }
                r
            })
            .collect();
        Csr::from_rows(m * m, rows)
    }

    #[test]
    fn aggregation_covers_all_connected_nodes() {
        let a = lap2d(20);
        let d = a.diagonal();
        let (agg, n) = aggregate(&strength(&a, &d, 0.08));
        assert!(agg.iter().all(|&g| g != NONE));
        assert!(n > 20 && n < 100, "{n} aggregates");
    }

    /// A row whose weak off-diagonals sum to more than its diagonal (here
    /// next to two rows with large diagonals) keeps a positive diagonal in
    /// the filtered matrix, and the smoothed prolongator stays finite.
    #[test]
    fn filtered_diagonal_stays_positive() {
        let a = Csr::from_rows(
            3,
            vec![
                vec![(0, 1.0), (1, -0.9), (2, -0.9)],
                vec![(0, -0.9), (1, 1e4)],
                vec![(0, -0.9), (2, 1e4)],
            ],
        );
        let strong = strength(&a, &a.diagonal(), 0.25);
        assert!(strong.iter().all(Vec::is_empty), "all weak");
        let af = filtered(&a, &strong);
        assert!(
            af.diagonal().iter().all(|&d| d > 0.0),
            "{:?}",
            af.diagonal()
        );
        let pt = tentative(&[0, 0, 0], 1);
        let p = smooth_prolongator(&af, &pt, 4.0 / 3.0, 10);
        assert!(p.val.iter().all(|v| v.is_finite()), "{:?}", p.val);
    }

    #[test]
    fn lambda_max_of_laplacian() {
        let a = lap2d(30);
        let lam = lambda_max(&a, &a.diagonal(), 50);
        // λ_max(D⁻¹A) = 1 + cos(π/(m+1)) < 2
        assert!(lam > 1.9 && lam <= 2.0, "{lam}");
    }

    #[test]
    fn dense_inverse() {
        for m in [1, 2, 6, 15] {
            let a = lap2d(m);
            let inv = dense_spd_inverse(&a).unwrap();
            let n = a.n_rows;
            for c in 0..n {
                let col: Vec<f64> = (0..n).map(|i| inv[i * n + c]).collect();
                let e = a.spmv(&col);
                for (i, v) in e.iter().enumerate() {
                    let want = if i == c { 1.0 } else { 0.0 };
                    assert!((v - want).abs() < 1e-12, "m {m}: A A⁻¹[{i}][{c}] = {v}");
                    assert_eq!(inv[i * n + c], inv[c * n + i]);
                }
            }
        }
    }

    /// A positive definite matrix whose diagonals span 10¹⁴ (the elliptic
    /// equation of ζ-f on a fine mesh) has its inverse; a singular one
    /// (a pure Neumann Laplacian) is refused at any scale.
    #[test]
    fn dense_inverse_is_scale_independent() {
        let d = [1e-7, 1.0, 1e7];
        let b = [[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]];
        let rows = (0..3)
            .map(|i| (0..3).map(|j| (j as u32, d[i] * b[i][j] * d[j])).collect())
            .collect();
        let a = Csr::from_rows(3, rows);
        assert!(dense_spd_inverse(&a).is_err(), "the strict test refuses it");
        let inv = dense_spd_inverse_with(&a, true).unwrap();
        // A = D B D, so A⁻¹ = D⁻¹ B⁻¹ D⁻¹ with B⁻¹ = ¼ [[3 2 1] [2 4 2] [1 2 3]]
        let b_inv = [[3.0, 2.0, 1.0], [2.0, 4.0, 2.0], [1.0, 2.0, 3.0]];
        for i in 0..3 {
            for j in 0..3 {
                let got = inv[i * 3 + j] * d[i] * d[j];
                let want = b_inv[i][j] / 4.0;
                assert!((got - want).abs() < 1e-12, "[{i}][{j}]: {got} {want}");
            }
        }
        for scale in [1e-8, 1.0, 1e8] {
            let n = [[1.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 1.0]];
            let rows = (0..3)
                .map(|i| (0..3).map(|j| (j as u32, scale * n[i][j])).collect())
                .collect();
            assert!(
                dense_spd_inverse(&Csr::from_rows(3, rows)).is_err(),
                "{scale}"
            );
        }
    }

    #[test]
    fn hierarchy_complexity() {
        let h = build(lap2d(64), &AmgSettings::default()).unwrap();
        assert!(h.levels.len() >= 2);
        assert!(h.coarse_n <= 256);
        assert!(h.info.operator_complexity < 2.0, "{:?}", h.info);
        // Galerkin operators stay symmetric
        for l in &h.levels {
            assert!(l.a.asymmetry() < 1e-12);
        }
    }

    #[test]
    fn singular_coarse_level_is_reported() {
        // pure Neumann Laplacian: constant null space
        let mut a = lap2d(8);
        let n = a.n_rows;
        for i in 0..n {
            let off: f64 = a.row(i).filter(|&(j, _)| j != i).map(|e| e.1).sum();
            let k = a.row_ptr[i]..a.row_ptr[i + 1];
            for kk in k {
                if a.col[kk] as usize == i {
                    a.val[kk] = -off;
                }
            }
        }
        let err = build(a, &AmgSettings::default()).unwrap_err();
        assert!(err.to_string().contains("singular"), "{err}");
    }
}
