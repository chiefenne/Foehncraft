//! A dense LU factorisation with partial pivoting for the Newton steps of
//! a circuit (a few hundred unknowns at most).

/// Solve A x = b in place (A row-major n × n, overwritten); `None` when A
/// is singular to working precision.
pub fn solve(a: &mut [f64], b: &mut [f64], n: usize) -> Option<()> {
    let scale = a
        .iter()
        .fold(0.0_f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    for k in 0..n {
        let (mut piv, mut big) = (k, a[k * n + k].abs());
        for i in k + 1..n {
            let v = a[i * n + k].abs();
            if v > big {
                big = v;
                piv = i;
            }
        }
        if big <= 1e-300 * scale || !big.is_finite() {
            return None;
        }
        if piv != k {
            for j in 0..n {
                a.swap(k * n + j, piv * n + j);
            }
            b.swap(k, piv);
        }
        let d = a[k * n + k];
        for i in k + 1..n {
            let f = a[i * n + k] / d;
            if f != 0.0 {
                for j in k..n {
                    a[i * n + j] -= f * a[k * n + j];
                }
                b[i] -= f * b[k];
            }
        }
    }
    for k in (0..n).rev() {
        let mut s = b[k];
        for j in k + 1..n {
            s -= a[k * n + j] * b[j];
        }
        b[k] = s / a[k * n + k];
    }
    Some(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn solves_a_pivoting_system() {
        let mut a = vec![0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 2.0, 0.0, 3.0];
        let x = [1.0, -2.0, 3.0];
        let mut b: Vec<f64> = (0..3)
            .map(|i| (0..3).map(|j| a[i * 3 + j] * x[j]).sum())
            .collect();
        super::solve(&mut a, &mut b, 3).unwrap();
        for i in 0..3 {
            assert!((b[i] - x[i]).abs() < 1e-14);
        }
    }
}
