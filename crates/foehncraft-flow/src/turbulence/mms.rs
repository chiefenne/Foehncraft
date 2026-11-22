//! Manufactured solutions of the turbulence models (stage 6 specification
//! 9.2, row 1; construction after Eça, Hoekstra, Hay & Pelletier 2007).
//!
//! The velocity is prescribed and frozen, a divergence-free cellular flow
//! in the unit square, u = U₀ (sin πx cos πy, −cos πx sin πy, 0), which
//! has no flux through the sides. k, the second variable (ω or ε), ζ (φ of
//! the elliptic-blending model), f (α) and the wall distance d are smooth
//! manufactured fields; every source term of the model is active and none
//! of its switches (the limiters, the realizability bounds, the
//! Kolmogorov bounds of T and L, the maxima of ν_t and arg₁) changes
//! branch in the domain, so the source is smooth
//! (`manufactured_fields_stay_on_one_branch`). The sides are velocity
//! inlets with the exact k and ω/ε as expressions (Dirichlet); ζ is
//! uniform on them (the inlet value of ζ is a number) and f has zero
//! gradient there, as the exact f does. z is a symmetry pair (2D).
//!
//! The source added to each equation is the residual of the continuous
//! equation for the exact fields,
//!
//!   S_φ = c u·∇φ − ∇·(Γ_φ ∇φ) + λ φ − Q_φ,
//!
//! with c = 1, λ = 0 for the transport equations and c = 0, λ = 1 (f,
//! Γ = L²) or 1/L² (α, Γ = 1) for the elliptic one, and Γ_φ and the source
//! Q_φ of the model written out here from the published equations
//! (independently of the kernels). The derivatives of the analytic fields
//! come from fourth-order central differences with a step of 10⁻³
//! (truncation 10⁻¹², round-off 10⁻¹⁰, far below the discretisation
//! errors), ∇Γ from the same differences of Γ.
//!
//! Only the turbulence equations are solved (`Flow::turbulence_iteration`)
//! until their residuals fall below 10⁻¹¹; the L2 error of every variable
//! over the cells against the exact cell-centre values, relative to the
//! largest exact value, gives the observed order between levels,
//! p = ln(e₁/e₂)/ln(h₁/h₂) with h = (V/N)^½.
//!
//! The ζ-f₀ form (`k-epsilon-zeta-f0`, PR 184) takes its own coefficients,
//! the sink of k, C_ε1 with 0.009/ζ and the f₀ equation. The curvature
//! correction of SST runs in a second prescribed flow ([`Flow::Mixed`])
//! whose strain axes turn along the streamlines, so that DS_ij/Dt acts.
//! The QCR changes the momentum equations, not the turbulence equations:
//! [`qcr_errors`] compares the explicit stress force of the solver with
//! the divergence of the written-out stress.

#![allow(clippy::needless_range_loop)]

use crate::case::FlowCase;
use std::path::Path;

/// A point of the plane.
type P = [f64; 2];

const U0: f64 = 0.25;
const NU: f64 = 1e-5;
const PI: f64 = std::f64::consts::PI;

/// The prescribed flow of a model's manufactured solution.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Flow {
    /// u = U₀ (sin πx cos πy, −cos πx sin πy): S_xy = 0, no DS_ij/Dt of
    /// the rotation term.
    Cellular,
    /// The stream function of [`Flow::Cellular`] with the modes
    /// (2, 1) and (1, 2) added (β = 1): a strain tensor that turns along
    /// the streamlines, so that every term of r̃ of the curvature
    /// correction is active.
    Mixed,
}

thread_local! {
    /// The flow of the test running in this thread.
    static FLOW: std::cell::Cell<Flow> = const { std::cell::Cell::new(Flow::Cellular) };
}

const BETA: f64 = 1.0;

fn vel(x: P) -> [f64; 2] {
    let (sx, cx) = (PI * x[0]).sin_cos();
    let (sy, cy) = (PI * x[1]).sin_cos();
    match FLOW.get() {
        Flow::Cellular => [U0 * sx * cy, -U0 * cx * sy],
        Flow::Mixed => {
            let (s2x, c2x) = (2.0 * PI * x[0]).sin_cos();
            let (s2y, c2y) = (2.0 * PI * x[1]).sin_cos();
            [
                U0 * (sx * cy + BETA * s2x * cy + BETA * sx * c2y),
                -U0 * (cx * sy + 2.0 * BETA * c2x * sy + 0.5 * BETA * cx * s2y),
            ]
        }
    }
}

/// The velocity components as expressions of the case file.
fn vel_text(flow: Flow) -> [String; 2] {
    match flow {
        Flow::Cellular => [
            format!("{U0} * sin(pi * x) * cos(pi * y)"),
            format!("-{U0} * cos(pi * x) * sin(pi * y)"),
        ],
        Flow::Mixed => [
            format!(
                "{U0} * (sin(pi * x) * cos(pi * y) + {BETA} * sin(2 * pi * x) * cos(pi * y) \
                 + {BETA} * sin(pi * x) * cos(2 * pi * y))"
            ),
            format!(
                "-{U0} * (cos(pi * x) * sin(pi * y) + {} * cos(2 * pi * x) * sin(pi * y) \
                 + {} * cos(pi * x) * sin(2 * pi * y))",
                2.0 * BETA,
                0.5 * BETA
            ),
        ],
    }
}

fn k_exact(x: P) -> f64 {
    0.05 * (1.0 + 0.3 * (2.0 * x[0] + 0.5).sin() * (1.5 * x[1] + 0.2).cos())
}
const K_TEXT: &str = "0.05 * (1 + 0.3 * sin(2 * x + 0.5) * cos(1.5 * y + 0.2))";

fn w_exact(x: P) -> f64 {
    10.0 * (1.0 + 0.25 * (1.7 * x[0] + 0.3).cos() * (2.0 * x[1] + 0.4).sin())
}
const W_TEXT: &str = "10 * (1 + 0.25 * cos(1.7 * x + 0.3) * sin(2 * y + 0.4))";

fn e_exact(x: P) -> f64 {
    0.09 * (1.0 + 0.25 * (1.7 * x[0] + 0.3).cos() * (2.0 * x[1] + 0.4).sin())
}
const E_TEXT: &str = "0.09 * (1 + 0.25 * cos(1.7 * x + 0.3) * sin(2 * y + 0.4))";

/// ζ (φ): 0.4 on the sides.
fn z_exact(x: P) -> f64 {
    0.4 + 0.15 * (PI * x[0]).sin() * (PI * x[1]).sin()
}
const Z_SIDE: f64 = 0.4;

/// f of ζ-f: zero normal derivative on the sides.
fn f_exact(x: P) -> f64 {
    0.2 * (1.0 + 0.3 * (PI * x[0]).cos() * (PI * x[1]).cos())
}

/// α of the elliptic-blending model: in (0, 1), zero normal derivative on
/// the sides.
fn a_exact(x: P) -> f64 {
    0.6 + 0.2 * (PI * x[0]).cos() * (PI * x[1]).cos()
}

/// The wall distance: F₁ of SST goes from about 1 (k-ω) near y = 0 to
/// about 0 (k-ε) near y = 1.
fn d_exact(x: P) -> f64 {
    0.06 + 0.6 * x[1] + 0.1 * x[0]
}

const H: f64 = 1e-3;

/// ∂f/∂x_i, fourth order.
fn d1(f: &dyn Fn(P) -> f64, x: P, i: usize) -> f64 {
    let at = |s: f64| {
        let mut y = x;
        y[i] += s * H;
        f(y)
    };
    (-at(2.0) + 8.0 * at(1.0) - 8.0 * at(-1.0) + at(-2.0)) / (12.0 * H)
}

fn grad(f: &dyn Fn(P) -> f64, x: P) -> [f64; 2] {
    [d1(f, x, 0), d1(f, x, 1)]
}

fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// ∇²f, fourth order.
fn lap(f: &dyn Fn(P) -> f64, x: P) -> f64 {
    let mut s = 0.0;
    for i in 0..2 {
        let at = |t: f64| {
            let mut y = x;
            y[i] += t * H;
            f(y)
        };
        s += (-at(2.0) + 16.0 * at(1.0) - 30.0 * at(0.0) + 16.0 * at(-1.0) - at(-2.0))
            / (12.0 * H * H);
    }
    s
}

/// The velocity gradient g_ij = ∂u_i/∂x_j (2D).
fn grad_u(x: P) -> [[f64; 2]; 2] {
    let mut g = [[0.0; 2]; 2];
    for (i, row) in g.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = d1(&|y| vel(y)[i], x, j);
        }
    }
    g
}

/// S = (2 S_ij S_ij)^½.
fn strain(x: P) -> f64 {
    let g = grad_u(x);
    let mut s2 = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            let sij = 0.5 * (g[i][j] + g[j][i]);
            s2 += 2.0 * sij * sij;
        }
    }
    s2.sqrt()
}

/// One equation of a model at a point: c u·∇φ − ∇·(Γ∇φ) + λφ = Q.
#[derive(Clone, Copy)]
struct Eq {
    conv: f64,
    gamma: f64,
    lin: f64,
    q: f64,
}

impl Eq {
    fn transport(gamma: f64, q: f64) -> Self {
        Self {
            conv: 1.0,
            gamma,
            lin: 0.0,
            q,
        }
    }
}

/// A model: its variables with their exact fields, and its equations at
/// a point.
trait Model {
    /// The case key of the second variable.
    fn second(&self) -> &'static str;
    fn fields(&self) -> Vec<fn(P) -> f64>;
    fn eqs(&self, x: P) -> Vec<Eq>;
    /// The prescribed flow.
    fn flow(&self) -> Flow {
        Flow::Cellular
    }
    /// The `turbulence.model` of the case, when it is not the name of the
    /// test.
    fn key(&self) -> Option<&'static str> {
        None
    }
    /// Further keys of the `[turbulence]` table.
    fn extra(&self) -> &'static str {
        ""
    }
}

/// S_ij of the prescribed flow.
fn sij(x: P, i: usize, j: usize) -> f64 {
    let g = grad_u(x);
    0.5 * (g[i][j] + g[j][i])
}

/// The rotation function f_r1 of the curvature correction of Smirnov &
/// Menter (2009), c_r1 = 1, c_r2 = 2, c_r3 = 1, from the definitions
/// f_r1 = max[0, min(f_rotation, 1.25)],
/// f_rotation = (1 + c_r1) 2r*/(1 + r*) [1 − c_r3 atan(c_r2 r̃)] − c_r1,
/// r* = S/Ω, r̃ = 2 Ω_ik S_jk DS_ij/Dt/(Ω D³), D² = max(S², 0.09 ω²), in
/// the inertial frame, with DS_ij/Dt = u·∇S_ij of the steady flow.
fn fr1_exact(x: P) -> f64 {
    let g = grad_u(x);
    let (mut s2, mut o2) = (0.0, 0.0);
    let mut s = [[0.0; 2]; 2];
    let mut w = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            s[i][j] = 0.5 * (g[i][j] + g[j][i]);
            w[i][j] = 0.5 * (g[i][j] - g[j][i]);
            s2 += 2.0 * s[i][j] * s[i][j];
            o2 += 2.0 * w[i][j] * w[i][j];
        }
    }
    let (sm, om) = (s2.sqrt(), o2.sqrt());
    let rs = sm / om;
    let d2 = s2.max(0.09 * w_exact(x).powi(2));
    let u = vel(x);
    let mut num = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            let ds = dot(u, grad(&|y| sij(y, i, j), x));
            let a: f64 = (0..2).map(|k| w[i][k] * s[j][k]).sum();
            num += a * ds;
        }
    }
    let rt = 2.0 * num / (om * d2 * d2.sqrt());
    let fr = 2.0 * 2.0 * rs / (1.0 + rs) * (1.0 - (2.0 * rt).atan()) - 1.0;
    fr.clamp(0.0, 1.25)
}

/// k-ω SST 2003 (Menter, Kuntz & Langtry 2003; the production limiter
/// 10 β* k ω, P̃_k in both equations); with `curv` the production of both
/// equations times f_r1 (Smirnov & Menter 2009), after the limiter.
struct Sst {
    curv: bool,
}

struct SstPoint {
    nu_t: f64,
    f1: f64,
    pk: f64,
    cross: f64,
}

fn sst_point(x: P, curv: bool) -> SstPoint {
    let (bs, a1, sw2) = (0.09, 0.31, 0.856);
    let (k, w, d) = (k_exact(x), w_exact(x), d_exact(x));
    let s = strain(x);
    let cross = 2.0 * sw2 * dot(grad(&k_exact, x), grad(&w_exact, x)) / w;
    let cd = cross.max(1e-10);
    let t1 = k.sqrt() / (bs * w * d);
    let t2 = 500.0 * NU / (d * d * w);
    let arg1 = t1.max(t2).min(4.0 * sw2 * k / (cd * d * d));
    let f1 = arg1.powi(4).tanh();
    let arg2 = (2.0 * t1).max(t2);
    let f2 = (arg2 * arg2).tanh();
    let nu_t = a1 * k / (a1 * w).max(s * f2);
    let pk = (nu_t * s * s).min(10.0 * bs * k * w);
    let pk = if curv { fr1_exact(x) * pk } else { pk };
    SstPoint {
        nu_t,
        f1,
        pk,
        cross,
    }
}

impl Model for Sst {
    fn second(&self) -> &'static str {
        "omega"
    }
    fn key(&self) -> Option<&'static str> {
        Some("k-omega-sst")
    }
    fn flow(&self) -> Flow {
        if self.curv {
            Flow::Mixed
        } else {
            Flow::Cellular
        }
    }
    fn extra(&self) -> &'static str {
        if self.curv {
            "curvature_correction = true\n"
        } else {
            ""
        }
    }
    fn fields(&self) -> Vec<fn(P) -> f64> {
        vec![k_exact, w_exact]
    }
    fn eqs(&self, x: P) -> Vec<Eq> {
        let p = sst_point(x, self.curv);
        let (k, w) = (k_exact(x), w_exact(x));
        let b = |a: f64, c: f64| p.f1 * a + (1.0 - p.f1) * c;
        let (sk, sw) = (b(0.85, 1.0), b(0.5, 0.856));
        let (gamma, beta) = (b(5.0 / 9.0, 0.44), b(0.075, 0.0828));
        vec![
            Eq::transport(NU + sk * p.nu_t, p.pk - 0.09 * k * w),
            Eq::transport(
                NU + sw * p.nu_t,
                gamma / p.nu_t * p.pk - beta * w * w + (1.0 - p.f1) * p.cross,
            ),
        ]
    }
}

/// k-ω SST with the quadratic constitutive relation (`turbulence.qcr`):
/// the k and ω equations are those of [`Sst`] (the production keeps the
/// linear stress); the stress enters the momentum equations, checked by
/// [`qcr_orders`].
struct SstQcr {
    c2013: bool,
}

impl Model for SstQcr {
    fn second(&self) -> &'static str {
        "omega"
    }
    fn fields(&self) -> Vec<fn(P) -> f64> {
        vec![k_exact, w_exact]
    }
    fn eqs(&self, x: P) -> Vec<Eq> {
        Sst { curv: false }.eqs(x)
    }
    fn key(&self) -> Option<&'static str> {
        Some("k-omega-sst")
    }
    fn extra(&self) -> &'static str {
        if self.c2013 {
            "qcr = \"qcr2013\"\n"
        } else {
            "qcr = \"qcr2000\"\n"
        }
    }
}

/// The standard k-ε model (Launder & Spalding 1974), with the default
/// Durbin bound ν_t ≤ k/(√3 S) (inactive here).
struct StandardKe;

impl Model for StandardKe {
    fn second(&self) -> &'static str {
        "epsilon"
    }
    fn fields(&self) -> Vec<fn(P) -> f64> {
        vec![k_exact, e_exact]
    }
    fn eqs(&self, x: P) -> Vec<Eq> {
        let (k, e) = (k_exact(x), e_exact(x));
        let s = strain(x);
        let nu_t = (0.09 * k * k / e).min(k / (3f64.sqrt() * s));
        let pk = nu_t * s * s;
        vec![
            Eq::transport(NU + nu_t, pk - e),
            Eq::transport(NU + nu_t / 1.3, e / k * (1.44 * pk - 1.92 * e)),
        ]
    }
}

/// The realizable k-ε model (Shih, Liou, Shabbir, Yang & Zhu 1995): C_μ
/// from the strain and the rotation, C₁ = max(0.43, η/(η + 5)) (0.43
/// here, η < 3.8), the destruction C₂ ε²/(k + √(νε)).
struct RealizableKe;

fn realizable_nu_t(x: P) -> f64 {
    let (k, e) = (k_exact(x), e_exact(x));
    let g = grad_u(x);
    let mut sij = [[0.0; 3]; 3];
    let mut wij = [[0.0; 3]; 3];
    for i in 0..2 {
        for j in 0..2 {
            sij[i][j] = 0.5 * (g[i][j] + g[j][i]);
            wij[i][j] = 0.5 * (g[i][j] - g[j][i]);
        }
    }
    let (mut ss, mut ww, mut sss) = (0.0, 0.0, 0.0);
    for i in 0..3 {
        for j in 0..3 {
            ss += sij[i][j] * sij[i][j];
            ww += wij[i][j] * wij[i][j];
            for l in 0..3 {
                sss += sij[i][j] * sij[j][l] * sij[l][i];
            }
        }
    }
    let st = ss.sqrt();
    let wt = (sss / (st * st * st)).clamp(-1.0 / 6f64.sqrt(), 1.0 / 6f64.sqrt());
    let phi = (6f64.sqrt() * wt).acos() / 3.0;
    let a_s = 6f64.sqrt() * phi.cos();
    let u_star = (ss + ww).sqrt();
    let c_mu = 1.0 / (4.0 + a_s * u_star * k / e);
    c_mu * k * k / e
}

impl Model for RealizableKe {
    fn second(&self) -> &'static str {
        "epsilon"
    }
    fn fields(&self) -> Vec<fn(P) -> f64> {
        vec![k_exact, e_exact]
    }
    fn eqs(&self, x: P) -> Vec<Eq> {
        let (k, e) = (k_exact(x), e_exact(x));
        let s = strain(x);
        let nu_t = realizable_nu_t(x);
        let pk = nu_t * s * s;
        let eta = s * k / e;
        let c1 = 0.43f64.max(eta / (eta + 5.0));
        vec![
            Eq::transport(NU + nu_t, pk - e),
            Eq::transport(
                NU + nu_t / 1.2,
                c1 * s * e - 1.9 * e * e / (k + (NU * e).sqrt()),
            ),
        ]
    }
}

/// The k-ε-ζ-f model of Hanjalić, Popovac & Hadžiabdić (2004):
/// T = max[min(k/ε, a/(√6 C_μ S ζ)), C_T (ν/ε)^½],
/// L = C_L max[min(k^{3/2}/ε, k^½/(√6 C_μ S ζ)), C_η ν^¾/ε^¼],
/// ν_t = C_μ ζ k T, C_ε1 = 1.4 (1 + 0.012/ζ); the f equation in the
/// divergence form −∇·(L²∇f) + f of the solver (the paper writes L²∇²f;
/// the two differ by ∇L²·∇f).
struct ZetaF {
    /// The f₀ form of Popovac (2023) with the choices of 6d2 and #141
    /// (`k-epsilon-zeta-f0`).
    f0: bool,
}

/// The coefficient of ν_t, ν_t = C_μ ζ k T, and the coefficient a of the
/// realizability bound of T and L (a = 0.6 √2 of our ζ-f₀ default: the
/// bound written with S = (2 S_ij S_ij)^½).
fn zeta_f_coefficients(f0: bool) -> (f64, f64) {
    if f0 {
        (0.09 / 0.41, 0.6 * 2f64.sqrt())
    } else {
        (0.22, 0.6)
    }
}

fn zeta_f_scales(x: P, f0: bool) -> (f64, f64) {
    let (k, e, z) = (k_exact(x), e_exact(x), z_exact(x));
    let s = strain(x);
    let (cm, a) = zeta_f_coefficients(f0);
    let den = 6f64.sqrt() * cm * s * z;
    let t = (k / e).min(a / den).max(6.0 * (NU / e).sqrt());
    let l = 0.36
        * (k.powf(1.5) / e)
            .min(k.sqrt() / den)
            .max(85.0 * (NU.powi(3) / e).powf(0.25));
    (t, l)
}

impl Model for ZetaF {
    fn second(&self) -> &'static str {
        "epsilon"
    }
    fn fields(&self) -> Vec<fn(P) -> f64> {
        vec![k_exact, e_exact, z_exact, f_exact]
    }
    fn eqs(&self, x: P) -> Vec<Eq> {
        let (k, e, z, f) = (k_exact(x), e_exact(x), z_exact(x), f_exact(x));
        let s = strain(x);
        let (t, l) = zeta_f_scales(x, self.f0);
        let nu_t = zeta_f_coefficients(self.f0).0 * z * k * t;
        let pk = nu_t * s * s;
        if self.f0 {
            // Popovac (2023): ν_t = (c_μ/ζ_eq) ζ k T; C_ε1 = 1.4 (1 +
            // 0.009/ζ); the sink of k max(k/T_K, 2νk/d²) with T_K =
            // max(k/ε, C_T (ν/ε)^½); the ζ equation
            // f₀ − ζ [P_k/ε + (N_ζ − 1)]/T and the f₀ equation
            // −L²∇²f₀ + f₀ = [(C_f1 + C_f2 P_k/ε)(2/3 − ζ) + (N_ζ − 1) ζ]/T
            // with N_ζ = 2
            let ce1 = 1.4 * (1.0 + 0.009 / z);
            let t_k = (k / e).max(6.0 * (NU / e).sqrt());
            let sink = (1.0 / t_k).max(2.0 * NU / d_exact(x).powi(2));
            return vec![
                Eq::transport(NU + nu_t, pk - sink * k),
                Eq::transport(NU + nu_t / 1.3, (ce1 * pk - 1.9 * e) / t),
                Eq::transport(NU + nu_t / 1.2, f - z * (pk / e + 1.0) / t),
                Eq {
                    conv: 0.0,
                    gamma: l * l,
                    lin: 1.0,
                    q: ((0.4 + 0.65 * pk / e) * (2.0 / 3.0 - z) + z) / t,
                },
            ];
        }
        let ce1 = 1.4 * (1.0 + 0.012 / z);
        vec![
            Eq::transport(NU + nu_t, pk - e),
            Eq::transport(NU + nu_t / 1.3, (ce1 * pk - 1.9 * e) / t),
            Eq::transport(NU + nu_t / 1.2, f - z / k * pk),
            Eq {
                conv: 0.0,
                gamma: l * l,
                lin: 1.0,
                q: (0.4 + 0.65 * pk / e) * (2.0 / 3.0 - z) / t,
            },
        ]
    }
}

/// The elliptic-blending k-ε of Billard & Laurence (2012), their Appendix
/// A: T = (k²/ε² + C_T² ν/ε)^½, L = C_L (k³/ε² + C_η² ν^{3/2}/ε^{1/2})^½,
/// ν_t = C_μ φ k min(T, T_lim) (T here), the halved viscous diffusion,
/// E = C_ε3 (1 − α)³ (k/ε) 2ν ν_t (∂_l ∂_j U_i)², C*_ε2 with
/// |∇·(ν_t/σ_k ∇k)/ε|^{3/2}, and α − L²∇²α = 1 with L² outside the
/// divergence.
struct Eb;

fn eb_point(x: P) -> (f64, f64, f64) {
    let (k, e, ph) = (k_exact(x), e_exact(x), z_exact(x));
    let t = ((k / e).powi(2) + 16.0 * NU / e).sqrt();
    let l = 0.164 * (k.powi(3) / (e * e) + 75.0 * 75.0 * NU.powf(1.5) / e.sqrt()).sqrt();
    let s = strain(x);
    let t_lim = 0.6 / (6f64.sqrt() * 0.22 * ph * s / 2f64.sqrt());
    let nu_t = 0.22 * ph * k * t.min(t_lim);
    (t, l, nu_t)
}

/// Σ (∂_l ∂_j U_i)², from second differences of the velocity.
fn velocity_hessian2(x: P) -> f64 {
    let mut g2 = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            for l in 0..2 {
                let v = d1(&|y| d1(&|z| vel(z)[i], y, j), x, l);
                g2 += v * v;
            }
        }
    }
    g2
}

impl Model for Eb {
    fn second(&self) -> &'static str {
        "epsilon"
    }
    fn fields(&self) -> Vec<fn(P) -> f64> {
        vec![k_exact, e_exact, z_exact, a_exact]
    }
    fn eqs(&self, x: P) -> Vec<Eq> {
        let (k, e, ph, al) = (k_exact(x), e_exact(x), z_exact(x), a_exact(x));
        let s = strain(x);
        let (t, l, nu_t) = eb_point(x);
        let pk = nu_t * s * s;
        let b3 = (1.0 - al).powi(3);
        // E = C_ε3 (1 − α)³ (k/ε) 2ν ν_t G
        let big_e = 2.3 * b3 * (k / e) * 2.0 * NU * nu_t * velocity_hessian2(x);
        // ∇·(ν_t/σ_k ∇k), σ_k = 1
        let gk = grad(&k_exact, x);
        let gnt = grad(&|y| eb_point(y).2, x);
        let dkt = dot(gnt, gk) + nu_t * lap(&k_exact, x);
        let xx = (dkt / e).abs();
        let a3 = al.powi(3);
        let ce2 = 1.83 + a3 * (0.4 - 1.83) * xx.powf(1.5).tanh();
        let fw = -0.5 * e * ph / k;
        let fh = -(0.7 + 0.9 * pk / e) * (ph - 2.0 / 3.0) / t;
        let cross = 2.0 / k * nu_t * dot(grad(&z_exact, x), gk);
        vec![
            Eq::transport(0.5 * NU + nu_t, pk - e - big_e),
            Eq::transport(0.5 * NU + nu_t / 1.5, (1.44 * pk - ce2 * e) / t),
            Eq::transport(
                0.5 * NU + nu_t,
                (1.0 - a3) * fw + a3 * fh - ph / k * pk + cross,
            ),
            Eq {
                conv: 0.0,
                gamma: 1.0,
                lin: 1.0 / (l * l),
                q: 1.0 / (l * l),
            },
        ]
    }
}

/// S_φ = c u·∇φ − ∇·(Γ∇φ) + λφ − Q of every equation at a point.
fn added_source(m: &dyn Model, x: P) -> Vec<f64> {
    let u = vel(x);
    let eqs = m.eqs(x);
    m.fields()
        .into_iter()
        .enumerate()
        .map(|(e, phi)| {
            let g = grad(&phi, x);
            let q = eqs[e];
            let gg = grad(&|y| m.eqs(y)[e].gamma, x);
            let div = dot(gg, g) + q.gamma * lap(&phi, x);
            q.conv * dot(u, g) - div + q.lin * phi(x) - q.q
        })
        .collect()
}

fn model_of(name: &str) -> Box<dyn Model> {
    match name {
        "k-omega-sst" => Box::new(Sst { curv: false }),
        "k-omega-sst-curvature" => Box::new(Sst { curv: true }),
        "k-omega-sst-qcr2000" => Box::new(SstQcr { c2013: false }),
        "k-omega-sst-qcr2013" => Box::new(SstQcr { c2013: true }),
        "k-epsilon" => Box::new(StandardKe),
        "realizable-k-epsilon" => Box::new(RealizableKe),
        "k-epsilon-zeta-f" => Box::new(ZetaF { f0: false }),
        "k-epsilon-zeta-f0" => Box::new(ZetaF { f0: true }),
        "elliptic-blending-k-epsilon" => Box::new(Eb),
        _ => unreachable!(),
    }
}

fn case_text(model: &str, mesh: &str, m: &dyn Model) -> String {
    let key = m.key().unwrap_or(model);
    let second = m.second();
    let st = if second == "omega" { W_TEXT } else { E_TEXT };
    let zeta = if m.fields().len() > 2 {
        format!("zeta = {Z_SIDE}\n")
    } else {
        String::new()
    };
    let [ux, uy] = vel_text(m.flow());
    let side = |name: &str| {
        format!(
            "[boundary.{name}]\ntype = \"velocity-inlet\"\n\
             velocity = [\"{ux}\", \"{uy}\", 0.0]\n\
             k = \"{K_TEXT}\"\n{second} = \"{st}\"\n{zeta}"
        )
    };
    format!(
        r#"
[mesh.generate]
type = "{mesh}"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [16, 16, 1]
perturb = 0.15
[fluid]
viscosity = {NU}
[flow]
convection = true
[turbulence]
model = "{key}"
{extra}
{xmin}{xmax}{ymin}{ymax}
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver.turbulence]
limiter = "none"
relaxation = 0.9
tolerance = 1e-3
max_iter = 200
[output]
vtkhdf = false
wall_csv = false
"#,
        extra = m.extra(),
        xmin = side("xmin"),
        xmax = side("xmax"),
        ymin = side("ymin"),
        ymax = side("ymax"),
    )
}

/// Size h, L2 and L∞ errors per variable (relative to the largest exact
/// value), iterations.
struct Level {
    h: f64,
    l2: Vec<f64>,
    linf: Vec<f64>,
    iterations: usize,
    /// Cell centroid and the error of k relative to its largest value, per
    /// cell (for the gallery figure of the curvature correction).
    cells: Vec<([f64; 2], f64)>,
}

/// The flow of a level with the manufactured fields set (10 % off) and the
/// added sources.
struct Built {
    fl: crate::simplec::Flow,
    xc: Vec<P>,
    vol: Vec<f64>,
    exact: Vec<Vec<f64>>,
    n: usize,
}

fn build(case: &FlowCase, level: u32, m: &dyn Model) -> Built {
    let (mesh, geo, _) = crate::case::mesh_for(case, Some(level)).unwrap();
    let phys = crate::case::physics(case).unwrap();
    let values = case.compiled_boundaries().unwrap();
    let at = crate::bc::At {
        analytic: None,
        values: &values,
        t: 0.0,
        turning: &case.turning(),
    };
    let arrays = crate::bc::arrays(&mesh, &geo, &case.boundary, phys.rho, at).unwrap();
    let disc = crate::simplec::discretisation(&mesh, &geo, &arrays).unwrap();
    let phys = crate::case::physics_on(case, &phys, &mesh, &geo).unwrap();
    let nc = disc.0.n_cells;
    let xc: Vec<P> = (0..nc)
        .map(|c| [disc.0.cell_centroid.x[c], disc.0.cell_centroid.y[c]])
        .collect();
    let vol = disc.0.cell_volume.clone();
    let mut fl = crate::simplec::Flow::new(
        disc,
        arrays,
        phys,
        case.solver.clone(),
        [0.0; 3],
        0.0,
        Some(0),
        &crate::simplec::Target {
            device: "cpu".into(),
            threads: Some(1),
            ptx: Default::default(),
        },
    )
    .unwrap();
    crate::turbulence::attach(&mut fl, case, &mesh, &geo).unwrap();
    let mut u = foehncraft_fv::disc::Soa3::zeros(nc);
    for (c, &x) in xc.iter().enumerate() {
        let v = vel(x);
        u.x[c] = v[0];
        u.y[c] = v[1];
    }
    fl.initialise(&u, &vec![0.0; nc]).unwrap();
    let fields = m.fields();
    let n = fields.len();
    let mut src = vec![vec![0.0; nc]; n];
    for (c, &x) in xc.iter().enumerate() {
        for (e, v) in added_source(m, x).into_iter().enumerate() {
            src[e][c] = v;
        }
    }
    let exact: Vec<Vec<f64>> = fields
        .iter()
        .map(|f| xc.iter().map(|&x| f(x)).collect())
        .collect();
    // start 10 % off the exact fields
    let start: Vec<Vec<f64>> = exact
        .iter()
        .enumerate()
        .map(|(e, v)| {
            let s = if e % 2 == 0 { 0.9 } else { 1.1 };
            v.iter().map(|x| s * x).collect()
        })
        .collect();
    {
        let (t, dev) = fl.turbulence_mut().unwrap();
        // without walls the model takes d = 10¹⁰; the manufactured d
        // replaces it (SST and the sink of k of ζ-f₀ read it)
        let d: Vec<f64> = xc.iter().map(|&x| d_exact(x)).collect();
        t.set_wall_distance(dev, &d).unwrap();
        t.set_added_sources(dev, &src[0], &src[1]).unwrap();
        t.set_fields(dev, &start[0], &start[1]).unwrap();
        if n == 4 {
            t.set_added_zeta_sources(dev, &src[2], &src[3]).unwrap();
            t.set_zeta_fields(dev, &start[2], &start[3]).unwrap();
        }
    }
    Built {
        fl,
        xc,
        vol,
        exact,
        n,
    }
}

fn run_level(case: &FlowCase, level: u32, m: &dyn Model) -> Level {
    let Built {
        mut fl,
        vol,
        exact,
        n,
        xc,
        ..
    } = build(case, level, m);
    let nc = vol.len();
    let mut it = 0;
    loop {
        let r = fl.turbulence_iteration().unwrap();
        it += 1;
        let worst = r.residual.iter().cloned().fold(0.0, f64::max);
        if worst < 1e-11 {
            break;
        }
        assert!(
            it < 2500,
            "level {level}: residuals {:?} after {it} iterations",
            r.residual
        );
    }
    let st = fl.state().unwrap().turbulence.unwrap();
    let mut num = vec![st.k, st.omega];
    if let Some(z) = st.zeta {
        num.push(z.zeta);
        num.push(z.f);
    }
    let vt: f64 = vol.iter().sum();
    let mut l2 = vec![0.0; n];
    let mut linf = vec![0.0f64; n];
    for e in 0..n {
        let scale = exact[e].iter().cloned().fold(0.0, f64::max);
        for c in 0..nc {
            let err = num[e][c] - exact[e][c];
            l2[e] += vol[c] * err * err;
            linf[e] = linf[e].max(err.abs());
        }
        l2[e] = (l2[e] / vt).sqrt() / scale;
        linf[e] /= scale;
    }
    let kmax = exact[0].iter().cloned().fold(0.0, f64::max);
    Level {
        h: (vt / nc as f64).sqrt(),
        l2,
        linf,
        iterations: it,
        cells: (0..nc)
            .map(|c| (xc[c], (num[0][c] - exact[0][c]) / kmax))
            .collect(),
    }
}

fn sci(v: &[f64]) -> String {
    let s: Vec<String> = v.iter().map(|x| format!("{x:.3e}")).collect();
    format!("[{}]", s.join(", "))
}

/// The observed orders of the L2 errors of every variable over `levels`.
fn orders(model: &str, mesh: &str, levels: &[u32]) -> Vec<Vec<f64>> {
    let m = model_of(model);
    FLOW.set(m.flow());
    let case = FlowCase::from_toml(&case_text(model, mesh, m.as_ref()), Path::new(".")).unwrap();
    let runs: Vec<Level> = levels
        .iter()
        .map(|&l| run_level(&case, l, m.as_ref()))
        .collect();
    for r in &runs {
        eprintln!(
            "{model} {mesh}: h {:.4e} L2 {} Linf {} ({} iterations)",
            r.h,
            sci(&r.l2),
            sci(&r.linf),
            r.iterations
        );
    }
    runs.windows(2)
        .map(|w| {
            let lh = (w[0].h / w[1].h).ln();
            (0..w[0].l2.len())
                .map(|e| (w[0].l2[e] / w[1].l2[e]).ln() / lh)
                .collect()
        })
        .collect()
}

fn check(model: &str, mesh: &str, band: (f64, f64)) {
    let o = orders(model, mesh, &[0, 1, 2]);
    eprintln!("{model} {mesh}: orders {o:.3?}");
    for (e, p) in o.last().unwrap().iter().enumerate() {
        assert!(
            *p > band.0 && *p < band.1,
            "{model} {mesh}: order of variable {e} {p:.3} outside {band:?} ({o:.3?})"
        );
    }
}

/// τ_ij/ρ of the quadratic constitutive relation less the Boussinesq stress
/// (Spalart 2000; Mani et al. 2013) in the plane, from the definitions
/// τ_ij^QCR − τ_ij = −C_cr1 (O_ik τ_jk + O_jk τ_ik),
/// O_ik = 2 W_ik/(∂_m u_n ∂_m u_n)^½, W_ik = ½(∂u_i/∂x_k − ∂u_k/∂x_i),
/// τ_ij = 2 ν_t S_ij, C_cr1 = 0.3, and for QCR2013 the further
/// −C_cr2 ν_t (2 S_mn S_mn)^½ δ_ij, C_cr2 = 2.5.
fn qcr_tau(x: P, c2013: bool) -> [[f64; 2]; 2] {
    let g = grad_u(x);
    let nu_t = sst_point(x, false).nu_t;
    let mut gg = 0.0;
    let mut s = [[0.0; 2]; 2];
    let mut w = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            gg += g[i][j] * g[i][j];
            s[i][j] = 0.5 * (g[i][j] + g[j][i]);
            w[i][j] = 0.5 * (g[i][j] - g[j][i]);
        }
    }
    let norm = gg.sqrt();
    let o = |i: usize, k: usize| 2.0 * w[i][k] / norm;
    let tau = |i: usize, k: usize| 2.0 * nu_t * s[i][k];
    let mut out = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                out[i][j] -= 0.3 * (o(i, k) * tau(j, k) + o(j, k) * tau(i, k));
            }
        }
    }
    if c2013 {
        let s2: f64 = (0..2)
            .map(|i| (0..2).map(|j| s[i][j] * s[i][j]).sum::<f64>())
            .sum();
        for i in 0..2 {
            out[i][i] -= 2.5 * nu_t * (2.0 * s2).sqrt();
        }
    }
    out
}

/// ∂_j τ_ij of the QCR stress.
fn qcr_force(x: P, c2013: bool) -> [f64; 2] {
    let mut f = [0.0; 2];
    for (i, fi) in f.iter_mut().enumerate() {
        for j in 0..2 {
            *fi += d1(&|y| qcr_tau(y, c2013)[i][j], x, j);
        }
    }
    f
}

/// L2 error of the explicit stress force of the QCR on `mesh`, over the
/// cells away from the sides (|x − ½| ≤ 0.35 in both directions) and from the centre lines of the flow, where S = 0, relative
/// to the largest exact force, per level: the stress and its divergence
/// as the momentum equations take them, from the velocity gradient and the
/// ν_t of the solver on the manufactured fields (the cell-value stress on a
/// boundary face is first order, so the cells next to the sides are left
/// out).
fn qcr_errors(model: &str, mesh: &str, perturb: f64, levels: &[u32]) -> Vec<(f64, f64)> {
    let m = model_of(model);
    FLOW.set(m.flow());
    let c2013 = model.ends_with("2013");
    let text = case_text(model, mesh, m.as_ref())
        .replace("perturb = 0.15", &format!("perturb = {perturb}"));
    let case = FlowCase::from_toml(&text, Path::new(".")).unwrap();
    levels
        .iter()
        .map(|&level| {
            let Built {
                mut fl,
                xc,
                vol,
                exact,
                ..
            } = build(&case, level, m.as_ref());
            {
                // the exact fields as the start; converged for ν_t
                let (t, dev) = fl.turbulence_mut().unwrap();
                t.set_fields(dev, &exact[0], &exact[1]).unwrap();
            }
            for it in 0.. {
                let r = fl.turbulence_iteration().unwrap();
                if r.residual.iter().cloned().fold(0.0, f64::max) < 1e-11 {
                    break;
                }
                assert!(it < 2500, "level {level}: no convergence");
            }
            let force = fl.explicit_stress_force().unwrap();
            // FOEHN_MMS_QCR_DUMP=<folder>: the force error of every cell
            // (x, y, error of f_x and f_y, exact f_x and f_y) for the
            // gallery figure (tools/plot/qcr_consistency.py)
            let mut dump = std::env::var_os("FOEHN_MMS_QCR_DUMP").map(|d| {
                let mut s = String::from("x,y,ex,ey,fx,fy\n");
                s.reserve(64 * xc.len());
                (std::path::PathBuf::from(d), s)
            });
            let (mut num, mut vt, mut big) = (0.0, 0.0, 0.0f64);
            for (c, &x) in xc.iter().enumerate() {
                // S = 0 on the centre lines: |S| of the isotropic stress of
                // QCR2013 has a kink there, left out within ±0.1
                let kink = (PI * x[0]).cos().abs() * (PI * x[1]).cos().abs() < 0.3;
                if (x[0] - 0.5).abs() > 0.35 || (x[1] - 0.5).abs() > 0.35 || kink {
                    continue;
                }
                let f = qcr_force(x, c2013);
                let e = (force[c][0] - f[0]).powi(2) + (force[c][1] - f[1]).powi(2);
                if let Some((_, s)) = dump.as_mut() {
                    s.push_str(&format!(
                        "{},{},{},{},{},{}\n",
                        x[0],
                        x[1],
                        force[c][0] - f[0],
                        force[c][1] - f[1],
                        f[0],
                        f[1]
                    ));
                }
                num += vol[c] * e;
                vt += vol[c];
                big = big.max(f[0].abs()).max(f[1].abs());
            }
            if let Some((d, s)) = dump {
                let name = format!("{model}_{mesh}_{perturb}_{level}.csv");
                std::fs::write(d.join(name), s).unwrap();
            }
            let nc = vol.len() as f64;
            let h = (vol.iter().sum::<f64>() / nc).sqrt();
            (h, (num / vt).sqrt() / big)
        })
        .collect()
}

/// The order of the force error on meshes without perturbation.
fn check_qcr(model: &str, mesh: &str, band: (f64, f64)) {
    let e = qcr_errors(model, mesh, 0.0, &[0, 1, 2]);
    eprintln!("{model} {mesh}: (h, L2 error of the force) {e:?}");
    let p = (e[1].1 / e[2].1).ln() / (e[1].0 / e[2].0).ln();
    assert!(
        p > band.0 && p < band.1,
        "{model} {mesh}: order {p:.3} outside {band:?} ({e:?})"
    );
}

/// The manufactured fields keep every switch of the models on one branch
/// and let F₁ of SST cover its range.
#[test]
fn manufactured_fields_stay_on_one_branch() {
    let (mut f1_lo, mut f1_hi) = (1.0f64, 0.0f64);
    for i in 0..=20 {
        for j in 0..=20 {
            let x = [i as f64 / 20.0, j as f64 / 20.0];
            let (k, w, d, e, z) = (k_exact(x), w_exact(x), d_exact(x), e_exact(x), z_exact(x));
            let s = strain(x);
            // SST: ν_t = k/ω, the production limiter inactive, the √k
            // term of arg₁ wins over 500ν/(d²ω) and the CD_kω term
            let p = sst_point(x, false);
            f1_lo = f1_lo.min(p.f1);
            f1_hi = f1_hi.max(p.f1);
            assert!(0.31 * w > s, "a₁ω > S F₂ at {x:?}");
            assert!(p.nu_t * s * s < 10.0 * 0.09 * k * w, "limiter at {x:?}");
            let t1 = k.sqrt() / (0.09 * w * d);
            assert!(t1 > 500.0 * NU / (d * d * w), "arg₁ at {x:?}");
            assert!(
                t1 < 4.0 * 0.856 * k / (p.cross.max(1e-10) * d * d),
                "CD at {x:?}"
            );
            // k-ε: the Durbin bound and C₁ = 0.43 of the realizable model
            assert!(
                0.09 * k * k / e < k / (3f64.sqrt() * s.max(1e-30)),
                "Durbin at {x:?}"
            );
            assert!(s * k / e < 3.77, "η at {x:?}");
            // ζ-f: T = k/ε and L = C_L k^{3/2}/ε (no bound acts)
            for f0 in [false, true] {
                let (t, l) = zeta_f_scales(x, f0);
                assert!((t - k / e).abs() < 1e-15 * t, "T at {x:?}");
                assert!((l - 0.36 * k.powf(1.5) / e).abs() < 1e-15 * l, "L at {x:?}");
            }
            // ζ-f₀: the sink of k is k/T_K = ε (2ν/d² is smaller)
            assert!(2.0 * NU / (d * d) < e / k, "sink of k at {x:?}");
            // EB: T < T_lim, α in (0, 1), φ < 2
            let (t, _, nu_t) = eb_point(x);
            assert!(
                (nu_t - 0.22 * z * k * t).abs() < 1e-15 * nu_t,
                "T_lim at {x:?}"
            );
            assert!(a_exact(x) > 0.0 && a_exact(x) < 1.0 && z < 2.0);
        }
    }
    assert!(f1_lo < 0.05 && f1_hi > 0.95, "F₁ in [{f1_lo}, {f1_hi}]");
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn sst_manufactured_solution_orders() {
    check("k-omega-sst", "hex", (1.8, 2.2));
    check("k-omega-sst", "poly", (1.6, 2.4));
    check("k-omega-sst", "tet", (1.6, 2.4));
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn standard_k_epsilon_manufactured_solution_orders() {
    check("k-epsilon", "hex", (1.8, 2.2));
    check("k-epsilon", "poly", (1.6, 2.4));
    check("k-epsilon", "tet", (1.6, 2.4));
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn realizable_k_epsilon_manufactured_solution_orders() {
    check("realizable-k-epsilon", "hex", (1.8, 2.2));
    check("realizable-k-epsilon", "poly", (1.6, 2.4));
    check("realizable-k-epsilon", "tet", (1.6, 2.4));
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn zeta_f_manufactured_solution_orders() {
    check("k-epsilon-zeta-f", "hex", (1.8, 2.2));
    check("k-epsilon-zeta-f", "poly", (1.6, 2.4));
    check("k-epsilon-zeta-f", "tet", (1.6, 2.4));
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn elliptic_blending_manufactured_solution_orders() {
    check("elliptic-blending-k-epsilon", "hex", (1.8, 2.2));
    check("elliptic-blending-k-epsilon", "poly", (1.6, 2.4));
    check("elliptic-blending-k-epsilon", "tet", (1.6, 2.4));
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn zeta_f0_manufactured_solution_orders() {
    check("k-epsilon-zeta-f0", "hex", (1.8, 2.2));
    check("k-epsilon-zeta-f0", "poly", (1.6, 2.4));
    check("k-epsilon-zeta-f0", "tet", (1.6, 2.4));
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn curvature_correction_manufactured_solution_orders() {
    check("k-omega-sst-curvature", "hex", (1.8, 2.2));
    check("k-omega-sst-curvature", "poly", (1.6, 2.4));
}

/// On tetrahedra the curvature correction converges since D236: S_ij comes
/// from a second-order velocity gradient (the quadratic fit where the
/// face-neighbour one is first order) and its gradient from a quadratic fit
/// of S_ij on the cells where the least-squares one is not consistent. (Before, the error of k and ω stayed at about 1.3 % and 0.7 % at
/// all three levels.) Regular tetrahedra and 15 % perturbed ones.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn curvature_correction_on_tetrahedra_converges() {
    let m = model_of("k-omega-sst-curvature");
    FLOW.set(m.flow());
    for perturb in [0.0, 0.15] {
        let text = case_text("k-omega-sst-curvature", "tet", m.as_ref())
            .replace("perturb = 0.15", &format!("perturb = {perturb}"));
        let case = FlowCase::from_toml(&text, Path::new(".")).unwrap();
        let runs: Vec<Level> = (0..3).map(|l| run_level(&case, l, m.as_ref())).collect();
        for (level, r) in runs.iter().enumerate() {
            eprintln!(
                "curvature tet perturb {perturb}: h {:.4e} L2 {}",
                r.h,
                sci(&r.l2)
            );
            // FOEHN_MMS_CURV_DUMP=<folder>: the error of k per cell for the
            // gallery figure (tools/plot/curvature_consistency.py)
            if let Some(d) = std::env::var_os("FOEHN_MMS_CURV_DUMP") {
                let mut s = String::from("x,y,ek\n");
                for (x, e) in &r.cells {
                    s.push_str(&format!("{},{},{}\n", x[0], x[1], e));
                }
                let name = format!("curvature_tet_{perturb}_{level}.csv");
                std::fs::write(std::path::Path::new(&d).join(name), s).unwrap();
            }
        }
        for e in 0..2 {
            let p = (runs[1].l2[e] / runs[2].l2[e]).ln() / (runs[1].h / runs[2].h).ln();
            assert!(
                p > 1.6
                    && p < 2.4
                    && runs[2].l2[e] < runs[1].l2[e]
                    && runs[1].l2[e] < runs[0].l2[e],
                "perturb {perturb}, variable {e}: order {p:.3}"
            );
        }
    }
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn qcr_stress_manufactured_solution_orders() {
    for model in ["k-omega-sst-qcr2000", "k-omega-sst-qcr2013"] {
        check_qcr(model, "hex", (1.8, 2.2));
        check_qcr(model, "poly", (1.8, 2.2));
        // first order on tetrahedra: the midpoint rule of the face sum on
        // irregular cells (the face stress itself is second order, D236)
        check_qcr(model, "tet", (0.8, 1.3));
    }
}

/// On meshes with randomly perturbed vertices (hexahedra and tetrahedra,
/// 15 % of the cell size) the divergence of the explicit stress converges
/// (D236): the quadratic least-squares velocity gradient on the cells where
/// the face-neighbour gradient is first order and the skewness-corrected
/// face stress. (Before, Green-Gauss with the linear interpolation of cell
/// values left an error of about 5 % of the largest force at every level.)
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn qcr_stress_on_perturbed_meshes_converges() {
    for mesh in ["hex", "tet"] {
        let e = qcr_errors("k-omega-sst-qcr2000", mesh, 0.15, &[0, 1, 2]);
        eprintln!("perturbed {mesh}: {e:?}");
        let p = (e[1].1 / e[2].1).ln() / (e[1].0 / e[2].0).ln();
        assert!(p > 0.9 && p < 2.2, "perturbed {mesh}: order {p:.3} ({e:?})");
        assert!(
            e[2].1 < e[1].1 && e[1].1 < e[0].1,
            "perturbed {mesh}: {e:?}"
        );
    }
}
