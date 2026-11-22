//! Finite-volume operators on a scalar φ (stage 2): least-squares gradient
//! (D24), Venkatakrishnan limiters (D25), face flux coefficients with linear
//! upwind convection and over-relaxed non-orthogonal diffusion (D25, D26),
//! and the gather-based row assembly (D31).
//!
//! Data are structure-of-arrays slices (one slice per vector component), so
//! the same bodies run in a CPU loop and as one GPU thread per index.
//!
//! **Boundary faces.** Every boundary face carries two coefficients `a`, `b`
//! of its face value, φ_f = a φ_P + b:
//!
//! - fixed value φ_b: a = 0, b = φ_b; the gradient stencil point is the face
//!   centroid;
//! - fixed normal gradient g (zero for symmetry and zero-gradient faces):
//!   a = 1, b = g d_n with d_n = (x_f − x_P)·n̂; the stencil point is the
//!   projection x_P + d_n n̂, so the gradient of a linear field stays exact.
//!
//! With the diffusion factor of the face, D = Γ |S|²/(d·S), the same formula
//! −D (φ_f − φ_P) gives the diffusive flux of both kinds, because
//! D g d_n = Γ g |S| for a fixed gradient.
//!
//! **Face flux.** The flux of φ out of the owner P through face f is linear
//! in the cell values, F_f = c_P φ_P + c_N φ_N + b_f. The implicit part is
//! first-order upwind convection and the orthogonal part of the diffusion;
//! `b_f` carries the deferred corrections: the linear-upwind
//! reconstruction ṁ (φ_f − φ_U), blended by `blend` (0 upwind, 1 linear
//! upwind), or with `central` > 0 the bounded central scheme (see
//! [`central_weight`]) for that fraction of it, and the non-orthogonal
//! diffusion −Γ k·∇φ_f (Jasak 1996, the
//! over-relaxed decomposition S = Δ + k with Δ = d |S|²/(d·S)), switched by
//! `nonorth` and limited by `nonorth_limit`.

use crate::{ld, unpack_face};
use foehncraft_core::Real;

/// Three components of a vector field, one slice each.
#[derive(Clone, Copy)]
pub struct Vec3<'a, T> {
    /// x components.
    pub x: &'a [T],
    /// y components.
    pub y: &'a [T],
    /// z components.
    pub z: &'a [T],
}

impl<T: Real> Vec3<'_, T> {
    /// Vector `i`.
    ///
    /// # Safety
    /// `i` is inside all three slices.
    #[inline(always)]
    pub unsafe fn at(&self, i: usize) -> [T; 3] {
        // SAFETY: guaranteed by the caller.
        unsafe { [ld(self.x, i), ld(self.y, i), ld(self.z, i)] }
    }
}

#[inline(always)]
fn dot<T: Real>(a: [T; 3], b: [T; 3]) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
fn sub<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Face-based connectivity of a partition (see `foehncraft_mesh::Mesh`).
#[derive(Clone, Copy)]
pub struct Topology<'a> {
    /// Owner cell of each face.
    pub owner: &'a [u32],
    /// Neighbour cell of each interior face.
    pub neighbour: &'a [u32],
    /// Start of the face list of each local cell.
    pub cell_ptr: &'a [u32],
    /// Faces of each local cell, packed with the orientation bit.
    pub cell_faces: &'a [u32],
    /// Number of interior (inner and processor) faces; boundary faces
    /// follow.
    pub n_interior: usize,
    /// First periodic face (D29, D109 (1)): the interior faces
    /// `periodic..n_interior` are periodic, seen from their owner only,
    /// with the partner face's cell as neighbour; `n_interior` without
    /// periodic faces.
    pub periodic: usize,
}

#[inline(always)]
fn add<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// The centroid of the neighbour `n` of interior face `f` as seen from the
/// face: moved by the shift of a periodic face (`shift` has one entry per
/// periodic face), unchanged otherwise.
///
/// # Safety
/// `n` is inside `cell_centroid`, and for a periodic face `f − periodic`
/// inside `shift`.
#[inline(always)]
pub unsafe fn neighbour_centroid<T: Real>(
    f: usize,
    n: usize,
    topo: &Topology,
    cell_centroid: &Vec3<T>,
    shift: &Vec3<T>,
) -> [T; 3] {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let x = cell_centroid.at(n);
        if f >= topo.periodic {
            add(x, shift.at(f - topo.periodic))
        } else {
            x
        }
    }
}

/// Boundary coefficients of the face values, φ_f = a φ_P + b, one entry per
/// boundary face (index `f − n_interior`), and the jump of the periodic
/// faces.
#[derive(Clone, Copy)]
pub struct Boundary<'a, T> {
    /// Coefficient of the owner value.
    pub a: &'a [T],
    /// Constant part.
    pub b: &'a [T],
    /// Jump b of the affine map of the periodic faces (D109 (1), a = 1):
    /// the neighbour value seen across periodic face `f` is
    /// φ_N + jump\[f − periodic\]; empty: no jump (plain periodic).
    pub jump: &'a [T],
}

impl<T: Real> Boundary<'_, T> {
    /// The neighbour value `phi_n` of interior face `f` as seen from its
    /// owner: with the jump on a periodic face.
    ///
    /// # Safety
    /// For a periodic face with a jump, `f − topo.periodic` is inside
    /// `jump`.
    #[inline(always)]
    pub unsafe fn across(&self, f: usize, topo: &Topology, phi_n: T) -> T {
        if f >= topo.periodic && !self.jump.is_empty() {
            // SAFETY: guaranteed by the caller.
            phi_n + unsafe { ld(self.jump, f - topo.periodic) }
        } else {
            phi_n
        }
    }

    /// Face value of boundary face `fb` for owner value `phi_p`.
    ///
    /// # Safety
    /// `fb` is inside both slices.
    #[inline(always)]
    pub unsafe fn value(&self, fb: usize, phi_p: T) -> T {
        // SAFETY: guaranteed by the caller.
        unsafe { ld(self.a, fb) * phi_p + ld(self.b, fb) }
    }
}

/// Weighted least-squares gradient of cell `c` (D24):
/// ∇φ_c = Σ_k w_k (φ_k − φ_c) over the cell's faces in stored order, with
/// the weight vectors `w` precomputed on the host per cell-face entry,
/// w_k = G⁻¹ ω_k d_k, G = Σ ω_k d_k d_kᵀ.
///
/// # Safety
/// `c` is a local cell of `topo`, every index stored in `topo` is inside the
/// arrays it addresses, `w` has one entry per cell-face entry, `bc` one per
/// boundary face.
#[inline(always)]
pub unsafe fn lsq_gradient<T: Real>(
    c: usize,
    topo: &Topology,
    w: &Vec3<T>,
    phi: &[T],
    bc: &Boundary<T>,
) -> [T; 3] {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let pc = ld(phi, c);
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut g = [T::ZERO; 3];
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let d = if f < topo.n_interior {
                let other = if is_nb {
                    ld(topo.owner, f)
                } else {
                    ld(topo.neighbour, f)
                };
                bc.across(f, topo, ld(phi, other as usize)) - pc
            } else {
                bc.value(f - topo.n_interior, pc) - pc
            };
            let wk = w.at(k);
            g[0] += wk[0] * d;
            g[1] += wk[1] * d;
            g[2] += wk[2] * d;
            k += 1;
        }
        g
    }
}

/// Venkatakrishnan's limiter function of one face (Venkatakrishnan 1995):
/// for the unlimited increment Δ₋ = ∇φ·(x_f − x_P) and the allowed
/// increment Δ₊ (to the neighbourhood maximum or minimum),
///
/// ψ = (Δ₊² + ε² + 2 Δ₋ Δ₊) / (Δ₊² + 2 Δ₋² + Δ₋ Δ₊ + ε²),
///
/// clipped to 1. With ε = 0 the limited increment ψ Δ₋ never exceeds Δ₊, so
/// no new extremum is created; ε > 0 lets smooth regions keep the unlimited
/// gradient.
#[inline(always)]
pub fn venkat_psi<T: Real>(dplus: T, dminus: T, eps2: T) -> T {
    let num = dplus * dplus + eps2 + T::from_f64(2.0) * dminus * dplus;
    let den = dplus * dplus + T::from_f64(2.0) * dminus * dminus + dminus * dplus + eps2;
    (num / den).min_of(T::ONE)
}

/// Barth–Jespersen's limiter function of one face: ψ = min(1, Δ₊/Δ₋) for
/// the unlimited increment Δ₋ ≠ 0 and the allowed increment Δ₊ of the
/// same sign (or zero): the largest ψ with ψ Δ₋ inside the neighbourhood
/// bounds. Not differentiable, which can stall a steady residual (Berger,
/// Aftosmis & Murman 2005); the smooth Venkatakrishnan form is the default.
#[inline(always)]
pub fn barth_jespersen_psi<T: Real>(dplus: T, dminus: T) -> T {
    (dplus / dminus).max_of(T::ZERO).min_of(T::ONE)
}

/// Nishikawa's limiter functions R3, R4, R5 (Nishikawa 2022, AIAA
/// 2022-1374, eqs. (5.5)–(5.7)) with the ε modification of Venkatakrishnan:
/// for a = |Δ₊|, b = |Δ₋|,
///
/// ψ = ((aᵖ + εₚ) + a Sₚ) / ((aᵖ + εₚ) + b (aᵖ⁻¹ + Sₚ)) for a ≤ 2b, else 1,
///
/// S₃ = 4b², S₄ = 2b (a² − 2b (a − 2b)), S₅ = 8b² (a² − 2b (a − b)).
/// Differentiable, equal to 1 with zero slope where the unlimited value is
/// bounded (y = a/b ≥ 2), and of order p in smooth regions (the
/// Venkatakrishnan function is of order 2). `eps_p` = εₚ has the units of
/// φᵖ; `p` outside 3..=5 is taken as 3.
#[inline(always)]
pub fn nishikawa_psi<T: Real>(dplus: T, dminus: T, eps_p: T, p: u32) -> T {
    let a = dplus.abs();
    let b = dminus.abs();
    let two = T::from_f64(2.0);
    if a > two * b {
        return T::ONE;
    }
    let (ap, ap1, s) = match p {
        4 => {
            let a2 = a * a;
            (a2 * a2, a2 * a, two * b * (a2 - two * b * (a - two * b)))
        }
        5 => {
            let a2 = a * a;
            (
                a2 * a2 * a,
                a2 * a2,
                T::from_f64(8.0) * b * b * (a2 - two * b * (a - b)),
            )
        }
        _ => (a * a * a, a * a, T::from_f64(4.0) * b * b),
    };
    let num = ap + eps_p + a * s;
    let den = ap + eps_p + b * (ap1 + s);
    if den > T::ZERO {
        (num / den).max_of(T::ZERO).min_of(T::ONE)
    } else {
        T::ONE
    }
}

/// Minimum and maximum over the multi-dimensional limiting stencil of
/// vertex `v` (Park & Kim 2012, the MLP condition): the values of the cells
/// that share the vertex (`vert_ptr`, `vert_cells`) and the boundary values
/// of the boundary faces that contain it (`bface_ptr`, `bfaces`, indices of
/// boundary faces).
///
/// # Safety
/// The lists are consistent CSR arrays of the mesh; `phi` has an entry per
/// cell, the boundary arrays one per boundary face and `topo.owner` one per
/// face.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub unsafe fn vertex_bounds<T: Real>(
    v: usize,
    vert_ptr: &[u32],
    vert_cells: &[u32],
    bface_ptr: &[u32],
    bfaces: &[u32],
    topo: &Topology,
    phi: &[T],
    bc: &Boundary<T>,
) -> (T, T) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(vert_ptr, v) as usize;
        let hi = ld(vert_ptr, v + 1) as usize;
        if lo >= hi {
            // a vertex of no cell (no cell reads its bounds)
            return (T::ZERO, T::ZERO);
        }
        let first = ld(phi, ld(vert_cells, lo) as usize);
        let (mut mn, mut mx) = (first, first);
        let mut k = lo + 1;
        while k < hi {
            let x = ld(phi, ld(vert_cells, k) as usize);
            mn = mn.min_of(x);
            mx = mx.max_of(x);
            k += 1;
        }
        let lo = ld(bface_ptr, v) as usize;
        let hi = ld(bface_ptr, v + 1) as usize;
        k = lo;
        while k < hi {
            let fb = ld(bfaces, k) as usize;
            let p = ld(topo.owner, topo.n_interior + fb) as usize;
            let x = bc.value(fb, ld(phi, p));
            mn = mn.min_of(x);
            mx = mx.max_of(x);
            k += 1;
        }
        (mn, mx)
    }
}

/// The limiter function of [`mlp_limiter`].
pub const MLP_BARTH_JESPERSEN: u32 = 0;
/// [`venkat_psi`] with ε² = `eps`².
pub const MLP_VENKATAKRISHNAN: u32 = 2;

/// Limiter value ψ_c of cell `c` with the multi-dimensional limiting
/// process (MLP, Park & Kim 2012): the linear reconstruction at each vertex
/// v of the cell, Δ₋ = ∇φ_c·(x_v − x_c), is held inside the bounds of the
/// vertex ([`vertex_bounds`], `vmin`, `vmax`), Δ₊ = φ_v^max − φ_c or
/// φ_v^min − φ_c, with the limiter function `func`
/// ([`MLP_BARTH_JESPERSEN`], [`MLP_VENKATAKRISHNAN`] or 3, 4, 5 for
/// [`nishikawa_psi`]) and the threshold `eps` in the units of φ (ε of the
/// solution, Wang 2000: dimensionally consistent on any mesh size), εₚ =
/// epsᵖ. The stencil of the vertices bounds the face values of every
/// neighbour that sees the cell, so a local extremum does not grow (LED in
/// multi-dimensions, Park & Kim 2012).
///
/// # Safety
/// `cell_ptr`/`cell_verts` are a consistent CSR array of the vertices of
/// every cell, `points` and `vmin`/`vmax` have an entry per vertex, `phi`,
/// `grad`, `cell_centroid` one per cell.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub unsafe fn mlp_limiter<T: Real>(
    c: usize,
    cell_ptr: &[u32],
    cell_verts: &[u32],
    points: &Vec3<T>,
    vmin: &[T],
    vmax: &[T],
    phi: &[T],
    grad: &Vec3<T>,
    cell_centroid: &Vec3<T>,
    func: u32,
    eps: T,
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let pc = ld(phi, c);
        let g = grad.at(c);
        let xc = cell_centroid.at(c);
        let e2 = eps * eps;
        let ep = match func {
            4 => e2 * e2,
            5 => e2 * e2 * eps,
            _ => e2 * eps,
        };
        let lo = ld(cell_ptr, c) as usize;
        let hi = ld(cell_ptr, c + 1) as usize;
        let mut psi = T::ONE;
        let mut k = lo;
        while k < hi {
            let v = ld(cell_verts, k) as usize;
            let dm = dot(g, sub(points.at(v), xc));
            let dp = if dm > T::ZERO {
                ld(vmax, v) - pc
            } else {
                ld(vmin, v) - pc
            };
            let p = if dm == T::ZERO {
                T::ONE
            } else {
                match func {
                    MLP_BARTH_JESPERSEN => barth_jespersen_psi(dp, dm),
                    MLP_VENKATAKRISHNAN => venkat_psi(dp, dm, e2),
                    _ => nishikawa_psi(dp, dm, ep, func),
                }
            };
            psi = psi.min_of(p);
            k += 1;
        }
        psi
    }
}

/// Limiter value ψ_c of cell `c`: the minimum of [`venkat_psi`] over its
/// faces, with the bounds φ_min, φ_max from the cell and its face
/// neighbours (boundary faces count with their face value).
///
/// `eps2` is ε² per cell (Venkatakrishnan: (K h_c)³ with h_c the cube root
/// of the cell volume), or a single value when `eps2.len() == 1`
/// (Venkatakrishnan–Wang: ε² = (K' (φ_max − φ_min)_domain)², Wang 2000).
/// A single negative value selects Barth–Jespersen (Barth & Jespersen
/// 1989): ψ = min(1, Δ₊/Δ₋) per face, [`barth_jespersen_psi`].
///
/// # Safety
/// As [`lsq_gradient`]; `grad` and `cell_centroid` have an entry for `c`,
/// `face_centroid` one per face, `eps2` one per local cell or exactly one.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub unsafe fn venkat_limiter<T: Real>(
    c: usize,
    topo: &Topology,
    phi: &[T],
    bc: &Boundary<T>,
    grad: &Vec3<T>,
    cell_centroid: &Vec3<T>,
    face_centroid: &Vec3<T>,
    eps2: &[T],
) -> T {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let pc = ld(phi, c);
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let (mut mn, mut mx) = (pc, pc);
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            let v = if f < topo.n_interior {
                let other = if is_nb {
                    ld(topo.owner, f)
                } else {
                    ld(topo.neighbour, f)
                };
                bc.across(f, topo, ld(phi, other as usize))
            } else {
                bc.value(f - topo.n_interior, pc)
            };
            mn = mn.min_of(v);
            mx = mx.max_of(v);
            k += 1;
        }
        let e2 = if eps2.len() == 1 {
            ld(eps2, 0)
        } else {
            ld(eps2, c)
        };
        let g = grad.at(c);
        let xc = cell_centroid.at(c);
        let mut psi = T::ONE;
        k = lo;
        while k < hi {
            let (f, _) = unpack_face(ld(topo.cell_faces, k));
            let dm = dot(g, sub(face_centroid.at(f), xc));
            let p = if dm > T::ZERO {
                if e2 < T::ZERO {
                    barth_jespersen_psi(mx - pc, dm)
                } else {
                    venkat_psi(mx - pc, dm, e2)
                }
            } else if dm < T::ZERO {
                if e2 < T::ZERO {
                    barth_jespersen_psi(mn - pc, dm)
                } else {
                    venkat_psi(mn - pc, dm, e2)
                }
            } else {
                T::ONE
            };
            psi = psi.min_of(p);
            k += 1;
        }
        psi
    }
}

/// Per-face inputs of [`face_coefficients`].
#[derive(Clone, Copy)]
pub struct FaceInputs<'a, T> {
    /// Connectivity.
    pub topo: Topology<'a>,
    /// Boundary face values.
    pub bc: Boundary<'a, T>,
    /// Mass flux ṁ_f (kg/s), positive from owner to neighbour (out of the
    /// domain on a boundary face).
    pub mass_flux: &'a [T],
    /// Diffusion coefficient Γ_f on the face (kg/(m s)).
    pub gamma: &'a [T],
    /// Geometric diffusion factor |S|²/(d·S) per face (m), d from owner to
    /// neighbour or to the boundary face centroid.
    pub dgeo: &'a [T],
    /// Non-orthogonal part k = S − d |S|²/(d·S) per face (m²).
    pub k: Vec3<'a, T>,
    /// Weight of the owner in the face interpolation of the gradient.
    pub weight: &'a [T],
    /// Face centroids.
    pub face_centroid: Vec3<'a, T>,
    /// Cell centroids (local and halo cells).
    pub cell_centroid: Vec3<'a, T>,
    /// Shift of the neighbour centroid per periodic face
    /// ([`neighbour_centroid`]).
    pub shift: Vec3<'a, T>,
    /// Cell values.
    pub phi: &'a [T],
    /// Cell gradients.
    pub grad: Vec3<'a, T>,
    /// Limiter value per cell (1 for no limiter).
    pub psi: &'a [T],
    /// Blending of the linear-upwind correction: 0 first-order upwind, 1
    /// linear upwind.
    pub blend: T,
    /// Fraction of the bounded central scheme on interior faces (0: none;
    /// the rest takes the linear-upwind correction of `blend`).
    pub central: T,
    /// β of the normalised-variable bounding of the central value (0: no
    /// bounding), see [`central_weight`].
    pub beta: T,
    /// 1 to apply the non-orthogonal diffusion correction, 0 to omit it.
    pub nonorth: T,
    /// Limit of the non-orthogonal correction as a multiple of the
    /// orthogonal part (Jasak 1996); a negative value means no limit.
    pub nonorth_limit: T,
    /// 1 to clip the linear-upwind face value to the interval of the two
    /// cell values of the face, 0 to leave it ([`clip_increment`]).
    pub bound: T,
}

/// The linear-upwind increment `ho` = φ_f − φ_U clipped so that the face
/// value lies between the upwind and downwind values (`diff` = φ_D − φ_U):
/// the convection boundedness criterion of the normalised-variable
/// schemes (Gaskell & Lau 1988; Jasak, Weller & Gosman 1999) applied to
/// the face value of the reconstruction. A limiter alone bounds the face
/// value by the neighbourhood of the upwind cell, which on a
/// multi-dimensional mesh reaches cells that are not neighbours of the
/// downwind one; with the clip an inflow face never carries a value beyond
/// that of the receiving cell, so the converged cell values make no new
/// extrema (with a limiter that vanishes at a local extremum).
/// `bound` ≤ 0 leaves `ho` unchanged.
#[inline(always)]
pub fn clip_increment<T: Real>(ho: T, diff: T, bound: T) -> T {
    if bound <= T::ZERO {
        ho
    } else {
        ho.max_of(diff.min_of(T::ZERO)).min_of(diff.max_of(T::ZERO))
    }
}

/// With bounding (β > 0) the central value `c` clipped to the interval of
/// the upwind and downwind values `u`, `d` (the skewness correction can
/// leave it); without bounding unchanged.
#[inline(always)]
pub fn bounded_between<T: Real>(c: T, u: T, d: T, beta: T) -> T {
    if beta <= T::ZERO {
        c
    } else {
        c.max_of(u.min_of(d)).min_of(u.max_of(d))
    }
}

/// Weight γ ∈ [0, 1] of the central face value in the bounded central
/// scheme: φ_f = φ_U + γ (φ_f^C − φ_U), with φ_f^C the linear
/// interpolation. The normalised variable of the upwind cell,
/// φ̃_U = 1 − (φ_D − φ_U)/(2 ∇φ_U·d) with d = x_D − x_U (the far-upwind
/// value is extrapolated with the upwind gradient; Jasak, Weller & Gosman
/// 1999, after Leonard's normalised-variable diagram), decides: in the
/// monotone range β ≤ φ̃_U < 1 the central value (γ = 1), at a local
/// extremum (φ̃_U ≤ 0 or φ̃_U ≥ 1) upwind (γ = 0), in between the smooth
/// switch γ = φ̃_U/β. With the central value clipped to [φ_U, φ_D]
/// ([`bounded_between`]) the face value lies between φ_U and φ_D; this
/// bounds a 1D profile. On distorted multi-dimensional meshes the kernels
/// also scale the correction with the limiter of the upwind cell (as the
/// linear-upwind one), which with a Venkatakrishnan limiter keeps the
/// oblique step of the catalogue inside [0, 1]. `diff` = φ_D − φ_U,
/// `gd` = ∇φ_U·d; β ≤ 0 gives the plain central scheme (γ = 1).
#[inline(always)]
pub fn central_weight<T: Real>(diff: T, gd: T, beta: T) -> T {
    if beta <= T::ZERO {
        return T::ONE;
    }
    let den = T::from_f64(2.0) * gd;
    if den == T::ZERO {
        // a uniform field (central = upwind) or an extremum at U
        return if diff == T::ZERO { T::ONE } else { T::ZERO };
    }
    let r = T::ONE - diff / den;
    if r <= T::ZERO || r >= T::ONE {
        T::ZERO
    } else if r >= beta {
        T::ONE
    } else {
        r / beta
    }
}

/// Coefficients `(c_P, c_N, b_f)` of the flux out of the owner through face
/// `f`, F_f = c_P φ_P + c_N φ_N + b_f (`c_N = 0` on boundary faces).
///
/// # Safety
/// `f` is a face of the topology and every array of `inp` is long enough
/// for the indices it is read at.
#[inline(always)]
pub unsafe fn face_coefficients<T: Real>(f: usize, inp: &FaceInputs<T>) -> (T, T, T) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let t = &inp.topo;
        let p = ld(t.owner, f) as usize;
        let m = ld(inp.mass_flux, f);
        let gam = ld(inp.gamma, f);
        let dd = gam * ld(inp.dgeo, f);
        let mp = m.max_of(T::ZERO);
        let mm = m.min_of(T::ZERO);
        let xf = inp.face_centroid.at(f);
        let phi_p = ld(inp.phi, p);
        let gp = inp.grad.at(p);
        let kf = inp.k.at(f);
        if f < t.n_interior {
            let n = ld(t.neighbour, f) as usize;
            let gn = inp.grad.at(n);
            let xn = neighbour_centroid(f, n, t, &inp.cell_centroid, &inp.shift);
            let phi_n = inp.bc.across(f, t, ld(inp.phi, n));
            // linear-upwind deferred correction ṁ (φ_f − φ_U)
            let ho = if m >= T::ZERO {
                let ho = ld(inp.psi, p) * dot(gp, sub(xf, inp.cell_centroid.at(p)));
                clip_increment(ho, phi_n - phi_p, inp.bound)
            } else {
                let ho = ld(inp.psi, n) * dot(gn, sub(xf, xn));
                clip_increment(ho, phi_p - phi_n, inp.bound)
            };
            let w = ld(inp.weight, f);
            let wn = T::ONE - w;
            let gf = [
                w * gp[0] + wn * gn[0],
                w * gp[1] + wn * gn[1],
                w * gp[2] + wn * gn[2],
            ];
            let mut b = if inp.central > T::ZERO {
                // bounded central: the linear interpolation to x_w on the
                // line P–N, corrected along the face with the
                // interpolated gradient
                let xp = inp.cell_centroid.at(p);
                let xw = [
                    w * xp[0] + wn * xn[0],
                    w * xp[1] + wn * xn[1],
                    w * xp[2] + wn * xn[2],
                ];
                let phi_c = w * phi_p + wn * phi_n + dot(gf, sub(xf, xw));
                let (phi_u, phi_d, gu, d, psi_u) = if m >= T::ZERO {
                    (phi_p, phi_n, gp, sub(xn, xp), ld(inp.psi, p))
                } else {
                    (phi_n, phi_p, gn, sub(xp, xn), ld(inp.psi, n))
                };
                // the limiter of the upwind cell (1 without one) scales
                // the central correction as it scales the linear-upwind one
                let g = psi_u * central_weight(phi_d - phi_u, dot(gu, d), inp.beta);
                let phi_c = bounded_between(phi_c, phi_u, phi_d, inp.beta);
                m * ((T::ONE - inp.central) * inp.blend * ho + inp.central * g * (phi_c - phi_u))
            } else {
                inp.blend * m * ho
            };
            // non-orthogonal correction −Γ k·∇φ_f, ∇φ_f interpolated
            let corr = gam * dot(kf, gf);
            let orth = dd * (phi_n - phi_p);
            b -= inp.nonorth * limit_corr(corr, orth, inp.nonorth_limit);
            if f >= t.periodic && !inp.bc.jump.is_empty() {
                // the jump of the neighbour value in the implicit part
                b += (mm - dd) * ld(inp.bc.jump, f - t.periodic);
            }
            (mp + dd, mm - dd, b)
        } else {
            let fb = f - t.n_interior;
            let a = ld(inp.bc.a, fb);
            let bb = ld(inp.bc.b, fb);
            let one_a = T::ONE - a;
            let ho = clip_increment(
                ld(inp.psi, p) * dot(gp, sub(xf, inp.cell_centroid.at(p))),
                inp.bc.value(fb, phi_p) - phi_p,
                inp.bound,
            );
            let corr = gam * dot(kf, gp);
            let orth = dd * (inp.bc.value(fb, phi_p) - phi_p);
            // the central fraction takes the boundary value as the face
            // value of the outflow, ṁ⁺ (φ_b − φ_P)
            let lud = if inp.central > T::ZERO {
                (T::ONE - inp.central) * inp.blend
            } else {
                inp.blend
            };
            let mut b = mm * bb - dd * bb + lud * mp * ho
                - inp.nonorth * one_a * limit_corr(corr, orth, inp.nonorth_limit);
            if inp.central > T::ZERO {
                b += inp.central * mp * (inp.bc.value(fb, phi_p) - phi_p);
            }
            (mp + a * mm + one_a * dd, T::ZERO, b)
        }
    }
}

/// Limits the non-orthogonal correction `corr` to `lim` times the size of
/// the orthogonal part `orth`; `lim < 0` leaves it unchanged.
#[inline(always)]
fn limit_corr<T: Real>(corr: T, orth: T, lim: T) -> T {
    if lim < T::ZERO {
        corr
    } else {
        let bound = lim * orth.abs();
        corr.max_of(-bound).min_of(bound)
    }
}

/// Matrix storage of the assembly: sliced ELL (SELL-C-σ with σ = 1, rows in
/// cell order), slot 0 of each row is the diagonal, slots 1.. the interior
/// faces of the cell in stored order.
#[derive(Clone, Copy)]
pub struct SellLayout<'a> {
    /// Slice height C.
    pub chunk: usize,
    /// Start of each slice in `col`/`val`.
    pub slice_ptr: &'a [u32],
}

impl SellLayout<'_> {
    /// Storage position of slot `j` of row `r`.
    ///
    /// # Safety
    /// `r / chunk < slice_ptr.len()`.
    #[inline(always)]
    pub unsafe fn pos(&self, r: usize, j: usize) -> usize {
        let s = r / self.chunk;
        // SAFETY: guaranteed by the caller.
        unsafe { ld(self.slice_ptr, s) as usize + j * self.chunk + (r - s * self.chunk) }
    }
}

/// Assembly of row `c` by a gather over the cell's faces (D31): the row of
/// the conservation equation Σ_f ±F_f = S_c (out of the cell positive),
/// with implicit under-relaxation α (Patankar 1980):
/// a_P/α φ_P + Σ a_N φ_N = rhs + (1 − α)/α a_P φ_P.
///
/// Writes the diagonal and the off-diagonals of the row into `val` (SELL
/// slots), the right-hand side into `rhs[c]` and 1/a_P(relaxed) into
/// `dinv[c]`. Faces are summed in the stored order, so the result is the
/// same for any loop order over the rows.
///
/// # Safety
/// As [`lsq_gradient`]; `cp`, `cn`, `bf` have one entry per face; `val`
/// holds the SELL storage of `layout` with room for every interior face of
/// the cell; `rhs`, `dinv`, `source`, `phi` have an entry for `c`. No other
/// thread writes row `c`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub unsafe fn assemble_row<T: Real>(
    c: usize,
    topo: &Topology,
    cp: &[T],
    cn: &[T],
    bf: &[T],
    source: &[T],
    phi: &[T],
    alpha: T,
    layout: &SellLayout,
    val: *mut T,
    rhs: *mut T,
    dinv: *mut T,
) {
    // SAFETY: guaranteed by the caller.
    unsafe {
        let lo = ld(topo.cell_ptr, c) as usize;
        let hi = ld(topo.cell_ptr, c + 1) as usize;
        let mut diag = T::ZERO;
        let mut b = ld(source, c);
        let mut j = 1;
        let mut k = lo;
        while k < hi {
            let (f, is_nb) = unpack_face(ld(topo.cell_faces, k));
            if is_nb {
                diag -= ld(cn, f);
                b += ld(bf, f);
                *val.add(layout.pos(c, j)) = -ld(cp, f);
                j += 1;
            } else {
                diag += ld(cp, f);
                b -= ld(bf, f);
                if f < topo.n_interior {
                    *val.add(layout.pos(c, j)) = ld(cn, f);
                    j += 1;
                }
            }
            k += 1;
        }
        let dr = diag / alpha;
        b += (dr - diag) * ld(phi, c);
        *val.add(layout.pos(c, 0)) = dr;
        *rhs.add(c) = b;
        *dinv.add(c) = T::ONE / dr;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_clip_keeps_the_face_value_between_the_cells() {
        // off: the high-order increment as it is
        assert_eq!(clip_increment(0.7, 0.5, 0.0), 0.7);
        // inside [0, diff]
        assert_eq!(clip_increment(0.3, 0.5, 1.0), 0.3);
        // past the downwind cell, either sign
        assert_eq!(clip_increment(0.7, 0.5, 1.0), 0.5);
        assert_eq!(clip_increment(-0.7, -0.5, 1.0), -0.5);
        // against the upwind direction: back to the upwind value
        assert_eq!(clip_increment(-0.2, 0.5, 1.0), 0.0);
    }

    #[test]
    fn barth_jespersen_is_the_ratio_clipped_to_one() {
        assert_eq!(barth_jespersen_psi(0.25, 1.0), 0.25);
        assert_eq!(barth_jespersen_psi(3.0, 1.0), 1.0);
        assert_eq!(barth_jespersen_psi(-0.5, -1.0), 0.5);
        assert_eq!(barth_jespersen_psi(0.0, 1.0), 0.0);
    }

    /// Nishikawa's R3, R4 and R5 without ε: bounded (ψ y⁻¹ ≤ 1, so the
    /// reconstruction stays inside the bound), one from y = 2 on and
    /// differentiable there, below Barth–Jespersen at y = 1.
    #[test]
    fn nishikawa_functions_are_bounded_and_smooth() {
        for (p, at_one) in [(3, 5.0 / 6.0), (4, 0.875), (5, 0.9)] {
            let psi = |y: f64| nishikawa_psi(y, 1.0, 0.0, p);
            for i in 1..=3000 {
                let y = i as f64 * 1e-3;
                let v = psi(y);
                assert!((0.0..=1.0).contains(&v), "R{p} ψ({y}) = {v}");
                assert!(v <= y + 1e-15, "R{p} ψ({y}) = {v} > y");
                if y > 2.0 {
                    assert_eq!(v, 1.0);
                }
            }
            assert!(
                (psi(1.0) - at_one).abs() < 1e-12,
                "R{p} ψ(1) = {}",
                psi(1.0)
            );
            let slope = (psi(2.0) - psi(2.0 - 1e-4)) / 1e-4;
            assert!(slope.abs() < 1e-4, "R{p} slope at 2: {slope}");
            // the signs of Δ₊, Δ₋ do not matter, ε > 0 keeps ψ in [0, 1]
            assert_eq!(nishikawa_psi(-0.5, -1.0, 0.0, p), psi(0.5));
            let v = nishikawa_psi(0.0, 1.0, 1e-3, p);
            assert!(v > 0.0 && v < 1e-2, "R{p} {v}");
        }
    }
}
