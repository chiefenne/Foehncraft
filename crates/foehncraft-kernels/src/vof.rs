//! Kernel bodies of the volume-of-fluid method (L3, stage MP1): the face
//! value of the volume fraction α by THINC, CICSAM or upwind, the explicit
//! conservative update of α, the smoothing, the Green–Gauss gradient, the
//! interface curvature and the surface-tension face jump of the
//! balanced-force continuum surface force.
//!
//! Like [`crate::flow`], every body computes the outputs of one face or one
//! cell from its inputs (a gather, no atomics), so a CPU loop and one GPU
//! thread per index give bitwise the same results (D31).
//!
//! **Normalised variable.** The face schemes work in the normalised
//! variable of the donor–acceptor pair (Leonard 1991; on unstructured
//! meshes with the virtual far-upwind value of Darwish & Moukalled 2003,
//! α_UU = α_D − 2 ∇α_U·d, d from the donor U to the acceptor D, clipped to
//! [0, 1]): α̃ = (α − α_UU)/(α_D − α_UU). A donor outside the monotone
//! range 0 < α̃_U < 1 takes the upwind value.
//!
//! **THINC** (Xiao, Honma & Kono 2005): the donor cell holds the profile
//! H(x) = ½ (1 + tanh(β (x − x̃))) on its unit length along d, with x̃ from
//! the cell average α̃_U; the face value is the average of H over the part
//! of the cell that crosses the face in the step, (1/Co) ∫_{1−Co}^{1} H dx,
//! so that α̃_f Co ≤ α̃_U and (1 − α̃_f) Co ≤ 1 − α̃_U: bounded for a face
//! Courant number Co ≤ 1 (the explicit update needs the sum of the outflow
//! Courant numbers of a cell below 1, which the sub-steps of the caller
//! keep). The steepness follows the angle between the interface normal and
//! d (Xiao, Ii & Chen 2011): β_f = β |n̂·d̂| + 0.01.
//!
//! **CICSAM** (Ubbink & Issa 1999): the blend of Hyper-C and
//! ULTIMATE-QUICKEST by the angle factor γ_f = min(k_γ (cos θ_f + 1)/2, 1).

use crate::fv::{Topology, Vec3};
use crate::{ld, unpack_face};
use foehncraft_core::Real;

/// Face scheme of α: first-order upwind.
pub const SCHEME_UPWIND: u32 = 0;
/// Face scheme of α: THINC with the time-integrated face value.
pub const SCHEME_THINC: u32 = 1;
/// Face scheme of α: CICSAM.
pub const SCHEME_CICSAM: u32 = 2;

#[inline(always)]
fn dot<T: Real>(a: [T; 3], b: [T; 3]) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
fn sub<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline(always)]
fn clamp01<T: Real>(x: T) -> T {
    x.max_of(T::ZERO).min_of(T::ONE)
}

/// ln cosh x without overflow: |x| + ln(1 + e^{−2|x|}) − ln 2.
#[inline(always)]
pub fn ln_cosh<T: Real>(x: T) -> T {
    let a = x.abs();
    a + (T::ONE + (-(a + a)).exp_of()).ln_of() - T::from_f64(core::f64::consts::LN_2)
}

/// The normalised THINC face value for the donor average `c` (0 < c < 1),
/// the face Courant number `co` (0 ≤ co ≤ 1) and the steepness `beta`.
pub fn thinc_normalised<T: Real>(c: T, co: T, beta: T) -> T {
    let half = T::from_f64(0.5);
    // x̃ from the cell average: cosh(β − a)/cosh(a) = Q = e^{β(2c−1)},
    // a = β x̃, so tanh a = (cosh β − Q)/sinh β
    let eb = beta.exp_of();
    let (ch, sh) = (half * (eb + T::ONE / eb), half * (eb - T::ONE / eb));
    let q = (beta * (c + c - T::ONE)).exp_of();
    let lim = T::ONE - T::from_f64(1e-15);
    let y = ((ch - q) / sh).max_of(-lim).min_of(lim);
    let a = half * ((T::ONE + y) / (T::ONE - y)).ln_of();
    if co < T::from_f64(1e-9) {
        // the point value H(1)
        return half * (T::ONE + (beta - a).tanh_of());
    }
    // ½ + (1/(2βCo)) ln(cosh(β − a)/cosh(β(1 − Co) − a))
    let v = half + (ln_cosh(beta - a) - ln_cosh(beta * (T::ONE - co) - a)) / (beta * co + beta * co);
    clamp01(v)
}

/// The normalised CICSAM face value for the donor value `c` (0 < c < 1),
/// the face Courant number `co` and cos θ_f of the interface normal and d.
pub fn cicsam_normalised<T: Real>(c: T, co: T, cos_theta: T) -> T {
    let co = co.max_of(T::from_f64(1e-12));
    let hc = (c / co).min_of(T::ONE);
    let e = T::from_f64(8.0);
    let uq = ((e * co * c + (T::ONE - co) * (T::from_f64(6.0) * c + T::from_f64(3.0))) / e)
        .min_of(hc);
    let g = (T::from_f64(0.5) * (cos_theta.abs() + T::ONE)).min_of(T::ONE);
    g * hc + (T::ONE - g) * uq
}

/// The face value of α between the donor U (value `au`, gradient `gu`) and
/// the acceptor D (value `ad`), `d` from U to D, face Courant number `co`
/// (|F| Δt / V_U), steepness `beta` of THINC.
pub fn face_value<T: Real>(scheme: u32, au: T, ad: T, gu: [T; 3], d: [T; 3], co: T, beta: T) -> T {
    if scheme == SCHEME_UPWIND {
        return au;
    }
    let auu = clamp01(ad - (dot(gu, d) + dot(gu, d)));
    let den = ad - auu;
    if den.abs() < T::from_f64(1e-12) {
        return au;
    }
    let c = (au - auu) / den;
    let eps = T::from_f64(1e-6);
    if !(c > eps && c < T::ONE - eps) {
        return au;
    }
    let co = co.min_of(T::ONE);
    let gn = dot(gu, gu).sqrt_of();
    let dn = dot(d, d).sqrt_of();
    let cos = if gn > T::ZERO && dn > T::ZERO {
        dot(gu, d) / (gn * dn)
    } else {
        T::ZERO
    };
    let ct = if scheme == SCHEME_THINC {
        thinc_normalised(c, co, beta * cos.abs() + T::from_f64(0.01))
    } else {
        cicsam_normalised(c, co, cos)
    };
    clamp01(auu + ct * den)
}

/// Inputs of [`face_alpha`].
pub struct FaceIn<'a, T> {
    /// Connectivity.
    pub topo: &'a Topology<'a>,
    /// Cell centroids.
    pub cell_centroid: &'a Vec3<'a, T>,
    /// Cell volumes.
    pub volume: &'a [T],
    /// Volume flux F_f (m³/s) out of the owner.
    pub flux: &'a [T],
    /// α per cell.
    pub alpha: &'a [T],
    /// ∇α per cell.
    pub grad: &'a Vec3<'a, T>,
    /// α of the inflow per boundary face.
    pub alpha_b: &'a [T],
    /// Face scheme.
    pub scheme: u32,
    /// The sub-step Δt.
    pub dt: T,
    /// β of THINC.
    pub beta: T,
}

/// The face value α_f of face `f` (interior faces between the donor and
/// the acceptor of the flux; boundary faces: the owner's value for an
/// outflow, `alpha_b` for an inflow).
///
/// # Safety
/// The arrays have the lengths of the mesh; `f` is a face.
#[inline(always)]
pub unsafe fn face_alpha<T: Real>(f: usize, a: &FaceIn<'_, T>) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let p = ld(a.topo.owner, f) as usize;
        let fl = ld(a.flux, f);
        if f >= a.topo.n_interior {
            return if fl >= T::ZERO {
                ld(a.alpha, p)
            } else {
                ld(a.alpha_b, f - a.topo.n_interior)
            };
        }
        let n = ld(a.topo.neighbour, f) as usize;
        let (u, dn) = if fl >= T::ZERO { (p, n) } else { (n, p) };
        let d = sub(a.cell_centroid.at(dn), a.cell_centroid.at(u));
        let co = fl.abs() * a.dt / ld(a.volume, u);
        face_value(
            a.scheme,
            ld(a.alpha, u),
            ld(a.alpha, dn),
            a.grad.at(u),
            d,
            co,
            a.beta,
        )
    }
}

/// The explicit update of cell `c`: α* = α − Δt/V Σ ±F_f α_f; returns
/// (α clipped to [0, 1], the clipped part (α* − α_clipped) V in m³).
///
/// # Safety
/// As [`face_alpha`]; `c` is a local cell.
#[inline(always)]
pub unsafe fn update<T: Real>(
    c: usize,
    topo: &Topology<'_>,
    volume: &[T],
    flux: &[T],
    alpha_f: &[T],
    alpha: &[T],
    dt: T,
) -> (T, T) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut s = T::ZERO;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let q = ld(flux, f) * ld(alpha_f, f);
            if is_nb {
                s -= q;
            } else {
                s += q;
            }
            k += 1;
        }
        let v = ld(volume, c);
        let a = ld(alpha, c) - dt * s / v;
        let ac = clamp01(a);
        (ac, (a - ac) * v)
    }
}

/// Sum of the outflow Courant numbers Σ_out |F_f| Δt / V of cell `c` for
/// Δt = 1 (the caller scales it).
///
/// # Safety
/// As [`update`].
#[inline(always)]
pub unsafe fn outflow_rate<T: Real>(c: usize, topo: &Topology<'_>, volume: &[T], flux: &[T]) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut s = T::ZERO;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let q = ld(flux, f);
            if (is_nb && q < T::ZERO) || (!is_nb && q > T::ZERO) {
                s += q.abs();
            }
            k += 1;
        }
        s / ld(volume, c)
    }
}

/// The smoothed α̃ of cell `c`: ½ α_P + ½ Σ_f |S_f| α_f / Σ_f |S_f| with the
/// linear face value (the owner value on boundary faces).
///
/// # Safety
/// As [`update`]; `area` has an entry per face.
#[inline(always)]
pub unsafe fn smooth<T: Real>(
    c: usize,
    topo: &Topology<'_>,
    area: &Vec3<'_, T>,
    weight: &[T],
    alpha: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let (mut s, mut w) = (T::ZERO, T::ZERO);
        let mut k = lo;
        while k < hi {
            let (f, _) = unpack_face(ld(topo.cell_faces, k));
            let sf = area.at(f);
            let m = dot(sf, sf).sqrt_of();
            let af = if f < topo.n_interior {
                let p = ld(topo.owner, f) as usize;
                let n = ld(topo.neighbour, f) as usize;
                let wf = ld(weight, f);
                wf * ld(alpha, p) + (T::ONE - wf) * ld(alpha, n)
            } else {
                ld(alpha, c)
            };
            s += m * af;
            w += m;
            k += 1;
        }
        let h = T::from_f64(0.5);
        h * ld(alpha, c) + h * s / w
    }
}

/// The Green–Gauss gradient of cell `c` with the linear face value (the
/// owner value on boundary faces): (1/V) Σ ±φ_f S_f.
///
/// # Safety
/// As [`smooth`].
#[inline(always)]
pub unsafe fn gauss_gradient<T: Real>(
    c: usize,
    topo: &Topology<'_>,
    area: &Vec3<'_, T>,
    weight: &[T],
    volume: &[T],
    phi: &[T],
) -> [T; 3] {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut g = [T::ZERO; 3];
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let sf = area.at(f);
            let pf = if f < topo.n_interior {
                let p = ld(topo.owner, f) as usize;
                let n = ld(topo.neighbour, f) as usize;
                let wf = ld(weight, f);
                wf * ld(phi, p) + (T::ONE - wf) * ld(phi, n)
            } else {
                ld(phi, c)
            };
            let s = if is_nb { -pf } else { pf };
            g[0] += s * sf[0];
            g[1] += s * sf[1];
            g[2] += s * sf[2];
            k += 1;
        }
        let v = ld(volume, c);
        [g[0] / v, g[1] / v, g[2] / v]
    }
}

/// The unit normal n̂ = g/(|g| + δ).
#[inline(always)]
fn unit<T: Real>(g: [T; 3], delta: T) -> [T; 3] {
    let m = dot(g, g).sqrt_of() + delta;
    [g[0] / m, g[1] / m, g[2] / m]
}

/// The curvature κ = −∇·n̂ of cell `c` from the gradient `grad` of the
/// smoothed α: (−1/V) Σ ±n̂_f·S_f with n̂_f the unit vector of the linearly
/// interpolated gradient; boundary faces of kind 0 (walls, symmetry: a
/// contact angle of 90°) add nothing, others take the owner's n̂.
///
/// # Safety
/// As [`smooth`]; `wall` has an entry per boundary face.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub unsafe fn curvature<T: Real>(
    c: usize,
    topo: &Topology<'_>,
    area: &Vec3<'_, T>,
    weight: &[T],
    volume: &[T],
    grad: &Vec3<'_, T>,
    wall: &[u32],
    delta: T,
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut s = T::ZERO;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let sf = area.at(f);
            let nf = if f < topo.n_interior {
                let p = ld(topo.owner, f) as usize;
                let n = ld(topo.neighbour, f) as usize;
                let wf = ld(weight, f);
                let (gp, gn) = (grad.at(p), grad.at(n));
                let gf = [
                    wf * gp[0] + (T::ONE - wf) * gn[0],
                    wf * gp[1] + (T::ONE - wf) * gn[1],
                    wf * gp[2] + (T::ONE - wf) * gn[2],
                ];
                dot(unit(gf, delta), sf)
            } else if ld(wall, f - topo.n_interior) == 0 {
                T::ZERO
            } else {
                dot(unit(grad.at(c), delta), sf)
            };
            if is_nb {
                s -= nf;
            } else {
                s += nf;
            }
            k += 1;
        }
        -s / ld(volume, c)
    }
}

/// The surface-tension face offset σ κ_f (α_N − α_P) of interior face `f`
/// (0 on boundary faces), κ_f the mean of the two cell curvatures weighted
/// by |∇α̃| (0 when both are 0).
///
/// # Safety
/// As [`smooth`]; `kappa` and `gnorm` have an entry per cell.
#[inline(always)]
pub unsafe fn tension<T: Real>(
    f: usize,
    topo: &Topology<'_>,
    alpha: &[T],
    kappa: &[T],
    gnorm: &[T],
    sigma: T,
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        if f >= topo.n_interior {
            return T::ZERO;
        }
        let p = ld(topo.owner, f) as usize;
        let n = ld(topo.neighbour, f) as usize;
        let (wp, wn) = (ld(gnorm, p), ld(gnorm, n));
        let w = wp + wn;
        if !(w > T::ZERO) {
            return T::ZERO;
        }
        let kf = (wp * ld(kappa, p) + wn * ld(kappa, n)) / w;
        sigma * kf * (ld(alpha, n) - ld(alpha, p))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ∫_{1−Co}^{1} H dx / Co by the midpoint rule.
    fn thinc_quadrature(c: f64, co: f64, beta: f64) -> f64 {
        // x̃ by bisection on the cell average
        let avg = |xt: f64| {
            let n = 20000;
            (0..n)
                .map(|i| {
                    let x = (i as f64 + 0.5) / n as f64;
                    0.5 * (1.0 + (beta * (x - xt)).tanh())
                })
                .sum::<f64>()
                / n as f64
        };
        let (mut lo, mut hi) = (-10.0, 10.0);
        for _ in 0..200 {
            let m = 0.5 * (lo + hi);
            if avg(m) > c { lo = m } else { hi = m }
        }
        let xt = 0.5 * (lo + hi);
        let n = 20000;
        (0..n)
            .map(|i| {
                let x = 1.0 - co + co * (i as f64 + 0.5) / n as f64;
                0.5 * (1.0 + (beta * (x - xt)).tanh())
            })
            .sum::<f64>()
            / n as f64
    }

    #[test]
    fn thinc_matches_quadrature() {
        for &c in &[0.05, 0.3, 0.5, 0.7, 0.97] {
            for &co in &[0.0, 0.1, 0.5, 1.0] {
                for &b in &[1.0, 3.5, 6.0] {
                    let a: f64 = thinc_normalised(c, co, b);
                    let q = if co == 0.0 {
                        thinc_quadrature(c, 1e-6, b)
                    } else {
                        thinc_quadrature(c, co, b)
                    };
                    assert!((a - q).abs() < 2e-6, "c {c} co {co} beta {b}: {a} vs {q}");
                    // boundedness of the explicit update
                    assert!(a * co <= c + 1e-12 && (1.0 - a) * co <= 1.0 - c + 1e-12);
                }
            }
        }
        // a whole cell crossing gives back the cell average
        let a: f64 = thinc_normalised(0.3, 1.0, 3.5);
        assert!((a - 0.3).abs() < 1e-12);
    }

    #[test]
    fn cicsam_is_bounded() {
        for &c in &[0.01, 0.2, 0.5, 0.8, 0.99] {
            for &co in &[0.01, 0.2, 0.5, 1.0] {
                for &ct in &[0.0, 0.5, 1.0] {
                    let a: f64 = cicsam_normalised(c, co, ct);
                    assert!(a >= c - 1e-12 && a <= 1.0 + 1e-12);
                    assert!(a * co <= c + 1e-12, "c {c} co {co} {a}");
                }
            }
        }
    }

    #[test]
    fn face_value_falls_back_to_upwind() {
        // a donor outside the monotone range
        let a: f64 = face_value(SCHEME_THINC, 1.0, 1.0, [0.0; 3], [1.0, 0.0, 0.0], 0.2, 3.5);
        assert_eq!(a, 1.0);
        let b: f64 = face_value(SCHEME_UPWIND, 0.4, 1.0, [0.6, 0.0, 0.0], [1.0, 0.0, 0.0], 0.2, 3.5);
        assert_eq!(b, 0.4);
        // a sharpening face value for a jump at the face
        let c: f64 = face_value(SCHEME_THINC, 0.5, 1.0, [0.5, 0.0, 0.0], [1.0, 0.0, 0.0], 0.2, 3.5);
        assert!(c > 0.5 && c <= 1.0, "{c}");
    }
}
