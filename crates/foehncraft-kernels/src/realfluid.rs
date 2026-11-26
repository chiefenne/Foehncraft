//! Real-fluid thermodynamics of the sprays at high pressure (ECN Spray A,
//! stage SA1, D311, D312): the Peng-Robinson equation of state of a mixture,
//! its fugacity coefficients and departure enthalpies, and the phase
//! equilibrium at the surface of a droplet.
//!
//! **Equation of state** (Peng and Robinson 1976):
//! p = R T/(v − b) − a(T)/(v² + 2 b v − b²), for a component
//! a = Ω_a R² T_c²/p_c · [1 + κ (1 − √(T/T_c))]², b = Ω_b R T_c/p_c, with
//! κ = 0.37464 + 1.54226 ω − 0.26992 ω² (Peng and Robinson 1976) for
//! ω ≤ 0.491 and κ = 0.379642 + 1.48503 ω − 0.164423 ω² + 0.016666 ω³
//! (Robinson and Peng 1978) for heavier components. Mixtures by the van der
//! Waals one-fluid rules a = Σ Σ xᵢ xⱼ (1 − kᵢⱼ) √(aᵢ aⱼ), b = Σ xᵢ bᵢ. In
//! Z = p v/(R T), with A = a p/(R T)², B = b p/(R T):
//! Z³ − (1 − B) Z² + (A − 3B² − 2B) Z − (A B − B² − B³) = 0.
//!
//! **Fugacity and departure** (from the residual Helmholtz energy of the
//! equation, own derivation, the standard result):
//! ln φᵢ = (bᵢ/b)(Z − 1) − ln(Z − B) − A/(2√2 B) (2 Σⱼ xⱼ aᵢⱼ/a − bᵢ/b) Λ,
//! h − h^ig = R T (Z − 1) + (T da/dT − a)/(2√2 b) Λ,
//! Λ = ln[(Z + (1 + √2) B)/(Z + (1 − √2) B)].
//! The partial molar departure enthalpy of a component is
//! h̄ᵢ − hᵢ^ig = −R T² (∂ ln φᵢ/∂T)_{p,x}, taken by a centred difference
//! with a fixed step (the same bits on every target).
//!
//! **Surface equilibrium** (D311). At the droplet temperature T and the
//! thermodynamic pressure p the fugacity of every component is the same in
//! the surface liquid x and the vapour y next to it: yᵢ = Kᵢ xᵢ,
//! Kᵢ = φᵢ^l(x)/φᵢ^v(y). Component 0 is the fuel; the others are the gases of
//! the ambient, which keep their far-field ratios rⱼ in the vapour
//! (yⱼ = s rⱼ, s = 1 − y₀). For given K the compositions follow in closed
//! form from Σ x = 1: s = (K₀ − 1)/(K₀ Q − 1), Q = Σⱼ rⱼ/Kⱼ, x₀ = (1 − s)/K₀,
//! xⱼ = s rⱼ/Kⱼ (own derivation). The K values are found by successive
//! substitution from a start value (Wilson's K, or the solution at a nearby
//! temperature) and finished by Newton's method on ln K with a difference
//! Jacobian; both run a bounded number of steps. When the K values collapse
//! to 1 (the trivial solution, x = y) the droplet surface is at or above the
//! critical mixing state and no two-phase solution exists (D314).
//!
//! Everything is a function of its inputs with the in-house `exp`, `ln`,
//! `sqrt`, `acos` and `cos` of `foehncraft_core::math`, so CPU and GPU agree
//! bit for bit (D93 (6)).

use foehncraft_core::math;

/// Molar gas constant (J/(mol K)), exact in the SI since 2019.
pub const R: f64 = 8.314_462_618;

/// The largest number of components of a mixture.
pub const MAX: usize = 6;

const SQRT2: f64 = core::f64::consts::SQRT_2;

/// Ω_b of the Peng-Robinson equation. At T_c and p_c the cubic in Z has a
/// triple root Z_c: 1 − B = 3 Z_c, A − 3B² − 2B = 3 Z_c², A B − B² − B³ =
/// Z_c³ with A = Ω_a, B = Ω_b; solved to 50 digits (own calculation, the
/// test checks the three conditions).
pub const OMEGA_B: f64 = 0.077_796_073_903_888_46;
/// Ω_a of the Peng-Robinson equation (as [`OMEGA_B`]).
pub const OMEGA_A: f64 = 0.457_235_528_921_382_2;
/// The critical compressibility of the Peng-Robinson equation (as
/// [`OMEGA_B`]).
pub const Z_CRITICAL: f64 = 0.307_401_308_698_703_85;

/// One component: its critical point, acentric factor and molar mass.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Component {
    /// Critical temperature (K).
    pub tc: f64,
    /// Critical pressure (Pa).
    pub pc: f64,
    /// Acentric factor.
    pub omega: f64,
    /// Molar mass (kg/mol).
    pub molar_mass: f64,
    /// Mathias-Copeman coefficients c₁, c₂, c₃ of α(T) (Mathias and
    /// Copeman 1983), fitted to a measured vapour pressure; all zero: the
    /// κ(ω) of Peng and Robinson.
    pub mc: [f64; 3],
    /// Critical molar volume (m³/mol) for the dense-fluid transport of
    /// Chung et al. (1988) (D313); 0: estimated by [`Component::vc`].
    pub vc: f64,
    /// Dipole moment (debye) for Chung et al. (1988).
    pub dipole: f64,
    /// Association factor κ of Chung et al. (1988) (0 for all but strongly
    /// associating fluids).
    pub assoc: f64,
}

impl Component {
    /// A component with the κ(ω) of Peng and Robinson.
    pub const fn new(tc: f64, pc: f64, omega: f64, molar_mass: f64) -> Self {
        Self {
            tc,
            pc,
            omega,
            molar_mass,
            mc: [0.0; 3],
            vc: 0.0,
            dipole: 0.0,
            assoc: 0.0,
        }
    }

    /// The critical molar volume (m³/mol): the given one, or Z_c R T_c/p_c
    /// with Z_c = 0.291 − 0.080 ω (Pitzer et al. 1955; 0.5 % below the
    /// measured V_c of nitrogen, 2.4 % below that of n-dodecane).
    #[inline]
    pub fn vc(&self) -> f64 {
        if self.vc > 0.0 {
            self.vc
        } else {
            (0.291 - 0.080 * self.omega) * R * self.tc / self.pc
        }
    }

    /// κ(ω): Peng and Robinson (1976) up to ω = 0.491, Robinson and Peng
    /// (1978) above.
    #[inline]
    pub fn kappa(&self) -> f64 {
        let w = self.omega;
        if w <= 0.491 {
            0.374_64 + 1.542_26 * w - 0.269_92 * w * w
        } else {
            0.379_642 + 1.485_03 * w - 0.164_423 * w * w + 0.016_666 * w * w * w
        }
    }

    /// The co-volume b (m³/mol).
    #[inline]
    pub fn b(&self) -> f64 {
        OMEGA_B * R * self.tc / self.pc
    }

    /// The attraction a(T) (Pa m⁶/mol²) and its derivative da/dT:
    /// a = a_c α, √α = 1 + κ s (Peng and Robinson) or
    /// √α = 1 + c₁ s + c₂ s² + c₃ s³ below T_c and 1 + c₁ s above (Mathias
    /// and Copeman 1983), s = 1 − √(T/T_c).
    #[inline]
    pub fn a(&self, t: f64) -> (f64, f64) {
        let ac = OMEGA_A * R * R * self.tc * self.tc / self.pc;
        let s = 1.0 - math::sqrt(t / self.tc);
        // √α and d√α/ds
        let (m, dm) = if self.mc == [0.0; 3] {
            let k = self.kappa();
            (1.0 + k * s, k)
        } else if t < self.tc {
            let [c1, c2, c3] = self.mc;
            (
                1.0 + s * (c1 + s * (c2 + s * c3)),
                c1 + s * (2.0 * c2 + 3.0 * s * c3),
            )
        } else {
            (1.0 + self.mc[0] * s, self.mc[0])
        };
        // ds/dT = −1/(2 √(T T_c))
        (ac * m * m, -ac * m * dm / math::sqrt(t * self.tc))
    }
}

/// Which root of the cubic a phase takes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Root {
    /// The smallest Z above B.
    Liquid,
    /// The largest Z.
    Vapour,
    /// The root of the lower Gibbs energy.
    Stable,
}

/// A mixture of up to [`MAX`] components with its binary interaction
/// parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mixture {
    /// Number of components.
    pub n: usize,
    /// The components; component 0 is the fuel of a droplet.
    pub c: [Component; MAX],
    /// kᵢⱼ, symmetric, zero on the diagonal.
    pub k: [[f64; MAX]; MAX],
}

/// The state of one phase of fixed composition.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Phase {
    /// Compressibility factor Z = p v/(R T).
    pub z: f64,
    /// ln φᵢ.
    pub ln_phi: [f64; MAX],
    /// Molar departure enthalpy h − h^ig (J/mol).
    pub h_dep: f64,
    /// Molar volume (m³/mol).
    pub v: f64,
    /// Number of real roots of the cubic above B (1 or 3).
    pub roots: u32,
}

/// ∛x for any sign.
#[inline]
fn cbrt(x: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else if x > 0.0 {
        math::exp(math::ln(x) / 3.0)
    } else {
        -math::exp(math::ln(-x) / 3.0)
    }
}

/// The real roots of Z³ + c2 Z² + c1 Z + c0 = 0 in ascending order (the
/// count in the second value), each polished by two Newton steps.
pub fn cubic_roots(c2: f64, c1: f64, c0: f64) -> ([f64; 3], usize) {
    let p = c1 - c2 * c2 / 3.0;
    let q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0;
    let shift = -c2 / 3.0;
    let disc = 0.25 * q * q + p * p * p / 27.0;
    let mut r = [0.0; 3];
    let n;
    if disc > 0.0 {
        let s = math::sqrt(disc);
        r[0] = cbrt(-0.5 * q + s) + cbrt(-0.5 * q - s) + shift;
        n = 1;
    } else if p == 0.0 {
        r[0] = shift;
        n = 1;
    } else {
        let m = 2.0 * math::sqrt(-p / 3.0);
        let arg = (3.0 * q / (p * m)).clamp(-1.0, 1.0);
        let phi = math::acos(arg) / 3.0;
        let third = 2.094_395_102_393_195_5; // 2π/3
        let mut t = [
            m * math::cos(phi) + shift,
            m * math::cos(phi - third) + shift,
            m * math::cos(phi - 2.0 * third) + shift,
        ];
        // ascending, by a fixed network of compare-exchanges
        if t[0] > t[1] {
            t.swap(0, 1);
        }
        if t[1] > t[2] {
            t.swap(1, 2);
        }
        if t[0] > t[1] {
            t.swap(0, 1);
        }
        r = t;
        n = 3;
    }
    // two Newton steps, each kept only when it lowers the residual (near
    // a double or triple root the derivative vanishes and a step would
    // amplify the rounding)
    let f = |z: f64| ((z + c2) * z + c1) * z + c0;
    for z in r.iter_mut().take(n) {
        for _ in 0..2 {
            let fz = f(*z);
            let d = (3.0 * *z + 2.0 * c2) * *z + c1;
            if d != 0.0 {
                let zn = *z - fz / d;
                if f(zn).abs() < fz.abs() {
                    *z = zn;
                }
            }
        }
    }
    (r, n)
}

impl Mixture {
    /// A mixture of the first `n` of `c` with all kᵢⱼ zero.
    pub fn new(c: &[Component]) -> Self {
        let mut m = Self {
            n: c.len().min(MAX),
            c: [Component::default(); MAX],
            k: [[0.0; MAX]; MAX],
        };
        for (i, ci) in c.iter().take(MAX).enumerate() {
            m.c[i] = *ci;
        }
        m
    }

    /// Sets kᵢⱼ = kⱼᵢ.
    pub fn set_k(&mut self, i: usize, j: usize, k: f64) {
        if i < MAX && j < MAX && i != j {
            self.k[i][j] = k;
            self.k[j][i] = k;
        }
    }

    /// Molar mass of the composition x (kg/mol).
    pub fn molar_mass(&self, x: &[f64; MAX]) -> f64 {
        let mut m = 0.0;
        for (i, xi) in x.iter().enumerate().take(self.n) {
            m += xi * self.c[i].molar_mass;
        }
        m
    }

    /// The state of a phase of composition `x` at `t` (K) and `p` (Pa).
    pub fn phase(&self, x: &[f64; MAX], t: f64, p: f64, root: Root) -> Phase {
        let n = self.n;
        let mut ai = [0.0; MAX];
        let mut dai = [0.0; MAX];
        let mut bi = [0.0; MAX];
        for i in 0..n {
            let (a, da) = self.c[i].a(t);
            ai[i] = a;
            dai[i] = da;
            bi[i] = self.c[i].b();
        }
        // Σⱼ xⱼ aᵢⱼ, a, da/dT, b
        let mut sa = [0.0; MAX];
        let mut am = 0.0;
        let mut dam = 0.0;
        let mut bm = 0.0;
        for i in 0..n {
            let mut s = 0.0;
            for j in 0..n {
                let g = math::sqrt(ai[i] * ai[j]);
                let aij = (1.0 - self.k[i][j]) * g;
                s += x[j] * aij;
                // d√(aᵢ aⱼ)/dT = ½ √(aᵢ aⱼ) (aᵢ'/aᵢ + aⱼ'/aⱼ)
                dam += x[i]
                    * x[j]
                    * (1.0 - self.k[i][j])
                    * 0.5
                    * g
                    * (dai[i] / ai[i] + dai[j] / ai[j]);
            }
            sa[i] = s;
            am += x[i] * s;
            bm += x[i] * bi[i];
        }
        let rt = R * t;
        let a_ = am * p / (rt * rt);
        let b_ = bm * p / rt;
        let (r, nr) = cubic_roots(
            -(1.0 - b_),
            a_ - 3.0 * b_ * b_ - 2.0 * b_,
            -(a_ * b_ - b_ * b_ - b_ * b_ * b_),
        );
        // the roots above B
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        let mut count = 0;
        for &z in r.iter().take(nr) {
            if z > b_ {
                lo = lo.min(z);
                hi = hi.max(z);
                count += 1;
            }
        }
        if count == 0 {
            // cannot happen for a physical state; take the largest root
            lo = r[nr - 1].max(b_ * (1.0 + 1e-12));
            hi = lo;
            count = 1;
        }
        let gdep = |z: f64| {
            let lam = math::ln((z + (1.0 + SQRT2) * b_) / (z + (1.0 - SQRT2) * b_));
            z - 1.0 - math::ln(z - b_) - a_ / (2.0 * SQRT2 * b_) * lam
        };
        let z = match root {
            Root::Liquid => lo,
            Root::Vapour => hi,
            Root::Stable => {
                if count > 1 && gdep(lo) < gdep(hi) {
                    lo
                } else {
                    hi
                }
            }
        };
        let lam = math::ln((z + (1.0 + SQRT2) * b_) / (z + (1.0 - SQRT2) * b_));
        let lnzb = math::ln(z - b_);
        let mut ln_phi = [0.0; MAX];
        for i in 0..n {
            ln_phi[i] = bi[i] / bm * (z - 1.0)
                - lnzb
                - a_ / (2.0 * SQRT2 * b_) * (2.0 * sa[i] / am - bi[i] / bm) * lam;
        }
        Phase {
            z,
            ln_phi,
            h_dep: rt * (z - 1.0) + (t * dam - am) / (2.0 * SQRT2 * bm) * lam,
            v: z * rt / p,
            roots: count,
        }
    }

    /// Partial molar departure enthalpies h̄ᵢ − hᵢ^ig (J/mol) of a phase of
    /// composition `x`, −R T² ∂ln φᵢ/∂T by a centred difference of 0.01 K.
    pub fn partial_h_dep(&self, x: &[f64; MAX], t: f64, p: f64, root: Root) -> [f64; MAX] {
        let h = 0.01;
        let up = self.phase(x, t + h, p, root);
        let dn = self.phase(x, t - h, p, root);
        let mut out = [0.0; MAX];
        for (i, o) in out.iter_mut().enumerate().take(self.n) {
            *o = -R * t * t * (up.ln_phi[i] - dn.ln_phi[i]) / (2.0 * h);
        }
        out
    }

    /// Wilson's estimate of the K values, ln Kᵢ = ln(p_c,i/p) + 5.373
    /// (1 + ωᵢ)(1 − T_c,i/T) (Wilson 1968).
    pub fn wilson(&self, t: f64, p: f64) -> [f64; MAX] {
        let mut lk = [0.0; MAX];
        for (i, l) in lk.iter_mut().enumerate().take(self.n) {
            let c = &self.c[i];
            *l = math::ln(c.pc / p) + 5.373 * (1.0 + c.omega) * (1.0 - c.tc / t);
        }
        lk
    }

    /// Pure-component saturation pressure (Pa) of component `i` at `t` by
    /// the equal-fugacity iteration p ← p φ^l/φ^v from Wilson's estimate;
    /// `None` at or above its critical temperature or without convergence.
    pub fn saturation_pressure(&self, i: usize, t: f64) -> Option<f64> {
        let c = self.c[i];
        if t >= c.tc {
            return None;
        }
        let pure = Mixture::new(&[c]);
        let x = {
            let mut x = [0.0; MAX];
            x[0] = 1.0;
            x
        };
        let mut p = c.pc * math::exp(5.373 * (1.0 + c.omega) * (1.0 - c.tc / t));
        for _ in 0..200 {
            let l = pure.phase(&x, t, p, Root::Liquid);
            let v = pure.phase(&x, t, p, Root::Vapour);
            if l.roots < 2 {
                // one root: move towards the two-phase pressure range
                p *= if l.z < Z_CRITICAL { 0.9 } else { 1.1 };
                continue;
            }
            let pn = p * math::exp(l.ln_phi[0] - v.ln_phi[0]);
            if (pn / p - 1.0).abs() < 1e-13 {
                return Some(pn);
            }
            p = pn;
        }
        None
    }
}

/// The compositions of liquid and vapour for given K values (module
/// documentation); `None` when they leave the physical range (the fuel is
/// not less volatile than the gas: boiling or no droplet).
#[inline]
fn split(n: usize, lnk: &[f64; MAX], r: &[f64; MAX]) -> Option<([f64; MAX], [f64; MAX])> {
    let k0 = math::exp(lnk[0]);
    let mut q = 0.0;
    for j in 1..n {
        q += r[j] / math::exp(lnk[j]);
    }
    let den = k0 * q - 1.0;
    if den == 0.0 {
        return None;
    }
    let s = (k0 - 1.0) / den;
    if !(s > 0.0 && s <= 1.0) {
        return None;
    }
    let mut x = [0.0; MAX];
    let mut y = [0.0; MAX];
    y[0] = 1.0 - s;
    x[0] = y[0] / k0;
    for j in 1..n {
        y[j] = s * r[j];
        x[j] = y[j] / math::exp(lnk[j]);
    }
    Some((x, y))
}

/// The phase equilibrium at a droplet surface.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Surface {
    /// Liquid composition at the surface (mole fractions; component 0 is
    /// the fuel, the rest the dissolved gases).
    pub x: [f64; MAX],
    /// Vapour composition next to the surface (mole fractions).
    pub y: [f64; MAX],
    /// ln Kᵢ = ln(yᵢ/xᵢ).
    pub ln_k: [f64; MAX],
    /// Compressibility of the surface liquid and of the vapour.
    pub z_l: f64,
    /// Compressibility of the vapour.
    pub z_v: f64,
    /// Steps taken (substitution and Newton).
    pub steps: u32,
}

/// Outcome of [`surface`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SurfaceState {
    /// Two phases in equilibrium.
    TwoPhase(Surface),
    /// The K values collapsed to 1: at or above the critical mixing state.
    Critical,
    /// The fuel boils at this temperature and pressure (its vapour fraction
    /// at the surface would reach 1).
    Boiling,
    /// No convergence in the bounded number of steps.
    NotConverged,
}

/// Steps of successive substitution before Newton takes over.
const SUBSTITUTIONS: u32 = 30;
/// Newton steps at most.
const NEWTON: u32 = 30;
/// Convergence of ln K, relative to 1 + max |ln Kᵢ|. The residual of a
/// cold fuel carries round-off from ln(Z − B) of its nearly incompressible
/// liquid, and its Newton steps stall at that floor (about 2e-11 at 267 K
/// and 6 MPa, more below): a step under [`STALL`] that is no longer half
/// the one before also ends the iteration. An absolute 1e-12 without that
/// rule failed at single temperatures and ended the surface table there
/// as if critical.
const TOL: f64 = 1e-12;
/// The largest relative Newton step accepted as the round-off floor.
const STALL: f64 = 1e-8;
/// ‖ln K‖² below which the solution is taken as trivial (x = y).
const TRIVIAL: f64 = 1e-10;

/// F, the liquid and vapour phases and their compositions x and y.
type Residual = ([f64; MAX], Phase, Phase, [f64; MAX], [f64; MAX]);

/// The residual F(ln K) = ln K − (ln φ^l(x) − ln φ^v(y)) and the phases.
#[inline]
fn residual(mix: &Mixture, t: f64, p: f64, lnk: &[f64; MAX], r: &[f64; MAX]) -> Option<Residual> {
    let (x, y) = split(mix.n, lnk, r)?;
    let l = mix.phase(&x, t, p, Root::Liquid);
    let v = mix.phase(&y, t, p, Root::Vapour);
    let mut f = [0.0; MAX];
    for i in 0..mix.n {
        f[i] = lnk[i] - (l.ln_phi[i] - v.ln_phi[i]);
    }
    Some((f, l, v, x, y))
}

/// Solves the small dense system J d = −f (Gaussian elimination with
/// partial pivoting, fixed order); `None` when singular.
#[allow(clippy::needless_range_loop)]
fn solve(n: usize, j: &mut [[f64; MAX]; MAX], f: &[f64; MAX]) -> Option<[f64; MAX]> {
    let mut b = [0.0; MAX];
    for i in 0..n {
        b[i] = -f[i];
    }
    for c in 0..n {
        let mut piv = c;
        for r in c + 1..n {
            if j[r][c].abs() > j[piv][c].abs() {
                piv = r;
            }
        }
        if j[piv][c] == 0.0 {
            return None;
        }
        j.swap(c, piv);
        b.swap(c, piv);
        for r in c + 1..n {
            let m = j[r][c] / j[c][c];
            for k in c..n {
                j[r][k] -= m * j[c][k];
            }
            b[r] -= m * b[c];
        }
    }
    let mut d = [0.0; MAX];
    for c in (0..n).rev() {
        let mut s = b[c];
        for k in c + 1..n {
            s -= j[c][k] * d[k];
        }
        d[c] = s / j[c][c];
    }
    Some(d)
}

/// The phase equilibrium at the surface of a fuel droplet (component 0) at
/// `t` (K) and `p` (Pa), the other components keeping the ratios `r` (r₀
/// unused, Σ rⱼ = 1) in the vapour. `start` holds ln K to start from
/// (Wilson's values when `None`).
pub fn surface(
    mix: &Mixture,
    t: f64,
    p: f64,
    r: &[f64; MAX],
    start: Option<&[f64; MAX]>,
) -> SurfaceState {
    let n = mix.n;
    let mut lnk = match start {
        Some(s) => *s,
        None => mix.wilson(t, p),
    };
    let mut steps = 0;
    let scale = |lk: &[f64; MAX]| {
        let mut m: f64 = 0.0;
        for v in lk.iter().take(n) {
            m = m.max(v.abs());
        }
        1.0 + m
    };
    let trivial = |lk: &[f64; MAX]| {
        let mut s = 0.0;
        for v in lk.iter().take(n) {
            s += v * v;
        }
        s < TRIVIAL
    };
    // successive substitution
    let mut done = false;
    while steps < SUBSTITUTIONS {
        steps += 1;
        let Some((f, ..)) = residual(mix, t, p, &lnk, r) else {
            return if math::exp(lnk[0]) >= 1.0 {
                SurfaceState::Boiling
            } else {
                SurfaceState::NotConverged
            };
        };
        let mut dmax: f64 = 0.0;
        for i in 0..n {
            lnk[i] -= f[i];
            dmax = dmax.max(f[i].abs());
        }
        if trivial(&lnk) {
            return SurfaceState::Critical;
        }
        if dmax < TOL * scale(&lnk) {
            done = true;
            break;
        }
    }
    // Newton with a difference Jacobian
    let mut it = 0;
    let mut last = f64::INFINITY;
    while !done && it < NEWTON {
        it += 1;
        steps += 1;
        let Some((f, ..)) = residual(mix, t, p, &lnk, r) else {
            return SurfaceState::NotConverged;
        };
        let mut jac = [[0.0; MAX]; MAX];
        for c in 0..n {
            let h = 1e-7 * (1.0 + lnk[c].abs());
            let mut lp = lnk;
            lp[c] += h;
            let Some((fp, ..)) = residual(mix, t, p, &lp, r) else {
                return SurfaceState::NotConverged;
            };
            for rr in 0..n {
                jac[rr][c] = (fp[rr] - f[rr]) / h;
            }
        }
        let Some(d) = solve(n, &mut jac, &f) else {
            return SurfaceState::NotConverged;
        };
        let mut dmax: f64 = 0.0;
        for i in 0..n {
            lnk[i] += d[i];
            dmax = dmax.max(d[i].abs());
        }
        if trivial(&lnk) {
            return SurfaceState::Critical;
        }
        let sc = scale(&lnk);
        if dmax < TOL * sc || (dmax < STALL * sc && dmax > 0.5 * last) {
            done = true;
        }
        last = dmax;
    }
    if !done {
        return SurfaceState::NotConverged;
    }
    let Some((_, l, v, x, y)) = residual(mix, t, p, &lnk, r) else {
        return SurfaceState::NotConverged;
    };
    if y[0] >= 0.999_999 {
        return SurfaceState::Boiling;
    }
    SurfaceState::TwoPhase(Surface {
        x,
        y,
        ln_k: lnk,
        z_l: l.z,
        z_v: v.z,
        steps,
    })
}

/// The bubble point of the liquid `x` at `t`: the pressure (Pa) and the
/// vapour composition, by the iteration p ← p Σ Kᵢ xᵢ, y ← K x/Σ K x from
/// `p_start` (Pa); `None` without convergence in 500 steps or at the
/// trivial solution. Used by the checks against measured bubble points
/// (VA.2).
pub fn bubble_pressure(
    mix: &Mixture,
    t: f64,
    x: &[f64; MAX],
    p_start: f64,
) -> Option<(f64, [f64; MAX])> {
    let n = mix.n;
    let mut p = p_start;
    let lk = mix.wilson(t, p);
    let mut y = [0.0; MAX];
    let mut sy = 0.0;
    for i in 0..n {
        y[i] = math::exp(lk[i]) * x[i];
        sy += y[i];
    }
    for v in y.iter_mut().take(n) {
        *v /= sy;
    }
    for _ in 0..500 {
        let l = mix.phase(x, t, p, Root::Liquid);
        let v = mix.phase(&y, t, p, Root::Vapour);
        let mut s = 0.0;
        let mut yn = [0.0; MAX];
        for i in 0..n {
            yn[i] = math::exp(l.ln_phi[i] - v.ln_phi[i]) * x[i];
            s += yn[i];
        }
        let mut dy: f64 = 0.0;
        let mut trivial = 0.0;
        for i in 0..n {
            yn[i] /= s;
            dy = dy.max((yn[i] - y[i]).abs());
            trivial += (yn[i] - x[i]) * (yn[i] - x[i]);
        }
        if trivial < 1e-14 {
            return None;
        }
        let pn = p * s;
        y = yn;
        if (pn / p - 1.0).abs() < 1e-11 && dy < 1e-11 {
            return Some((pn, y));
        }
        p = pn;
    }
    None
}

/// The outcome of [`Mixture::stability`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Stability {
    /// The smallest modified tangent-plane distance tm of the trial
    /// phases; negative when the phase splits.
    pub tm: f64,
    /// The composition of the trial phase of that minimum (mole
    /// fractions; the feed when the phase is stable).
    pub w: [f64; MAX],
    /// The trial phase is less dense than the feed (a vapour-like trial).
    pub lighter: bool,
    /// Substitution steps of both trials.
    pub steps: u32,
}

/// tm below which a phase is taken as unstable (the round-off of the sums
/// is near 1e-15).
pub const STABILITY_TOL: f64 = 1e-10;
/// Substitution steps of one trial phase at most.
const STABILITY_STEPS: u32 = 300;
/// Substitution steps of the flash before Newton takes over.
const FLASH_SUBSTITUTIONS: u32 = 60;

impl Stability {
    /// Whether the phase splits into two.
    pub fn unstable(&self) -> bool {
        self.tm < -STABILITY_TOL
    }
}

/// Two phases in equilibrium, from [`Mixture::flash`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Split {
    /// Mole fraction of the feed in the lighter phase.
    pub beta: f64,
    /// Composition of the denser phase (mole fractions).
    pub x: [f64; MAX],
    /// Composition of the lighter phase (mole fractions).
    pub y: [f64; MAX],
    /// ln Kᵢ = ln(yᵢ/xᵢ).
    pub ln_k: [f64; MAX],
    /// The denser phase.
    pub l: Phase,
    /// The lighter phase.
    pub v: Phase,
    /// Steps taken (substitution and Newton).
    pub steps: u32,
}

/// Outcome of [`Mixture::flash`] (kept by value: the kernels allocate
/// nothing).
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Flash {
    /// The feed is stable as one phase (the root of the lower Gibbs
    /// energy).
    One(Phase),
    /// The feed splits into two phases.
    Two(Split),
    /// The stability test found a split but the flash did not converge in
    /// its bounded number of steps.
    NotConverged,
}

/// The lighter-phase fraction β of the feed `z` for the K values `lnk`: the
/// root of the Rachford-Rice function Σ zᵢ (Kᵢ − 1)/(1 + β (Kᵢ − 1)) between
/// its poles 1/(1 − K_max) and 1/(1 − K_min) (β outside 0 to 1 is kept, a
/// negative flash), by Newton steps kept inside a shrinking bracket;
/// `None` when all K are on one side of 1.
fn rachford_rice(n: usize, z: &[f64; MAX], lnk: &[f64; MAX]) -> Option<f64> {
    let mut k = [1.0; MAX];
    let mut kmin = f64::INFINITY;
    let mut kmax: f64 = 0.0;
    for i in 0..n {
        if z[i] > 0.0 {
            k[i] = math::exp(lnk[i]);
            kmin = kmin.min(k[i]);
            kmax = kmax.max(k[i]);
        }
    }
    if !(kmax > 1.0 && kmin < 1.0) {
        return None;
    }
    let mut lo = 1.0 / (1.0 - kmax);
    let mut hi = 1.0 / (1.0 - kmin);
    let mut b: f64 = 0.5;
    for _ in 0..200 {
        let mut g = 0.0;
        let mut dg = 0.0;
        for i in 0..n {
            if z[i] > 0.0 {
                let c = k[i] - 1.0;
                let den = 1.0 + b * c;
                g += z[i] * c / den;
                dg -= z[i] * c * c / (den * den);
            }
        }
        if g == 0.0 {
            return Some(b);
        }
        // g falls with β: a positive g puts the root to the right
        if g > 0.0 {
            lo = b;
        } else {
            hi = b;
        }
        let mut bn = b - g / dg;
        if !(bn > lo && bn < hi) {
            bn = 0.5 * (lo + hi);
        }
        if (bn - b).abs() <= 1e-15 * (1.0 + b.abs()) {
            return Some(bn);
        }
        b = bn;
    }
    Some(b)
}

/// F, β, x, y and the denser and lighter phases of a flash.
type FlashResidual = ([f64; MAX], f64, [f64; MAX], [f64; MAX], Phase, Phase);

impl Mixture {
    /// The pressure (Pa) of the composition `x` at `t` (K) and the molar
    /// volume `v` (m³/mol), explicit in the equation of state.
    pub fn pressure(&self, x: &[f64; MAX], t: f64, v: f64) -> f64 {
        let n = self.n;
        let mut am = 0.0;
        let mut bm = 0.0;
        let mut ai = [0.0; MAX];
        for i in 0..n {
            ai[i] = self.c[i].a(t).0;
            bm += x[i] * self.c[i].b();
        }
        for i in 0..n {
            for j in 0..n {
                am += x[i] * x[j] * (1.0 - self.k[i][j]) * math::sqrt(ai[i] * ai[j]);
            }
        }
        R * t / (v - bm) - am / (v * v + 2.0 * bm * v - bm * bm)
    }

    /// The stability test of the phase `z` at `t` (K) and `p` (Pa) by the
    /// tangent plane of Michelsen (1982): with dᵢ = ln zᵢ + ln φᵢ(z) (the
    /// root of the lower Gibbs energy), a trial phase of mole numbers W is
    /// moved by ln Wᵢ ← dᵢ − ln φᵢ(W/ΣW) from a vapour-like (Wᵢ = Kᵢ zᵢ) and
    /// a liquid-like (Wᵢ = zᵢ/Kᵢ) start with Wilson's K. The modified
    /// distance tm = 1 + Σ Wᵢ (ln Wᵢ + ln φᵢ(w) − dᵢ − 1) is negative only
    /// if the tangent-plane distance of w = W/ΣW is (tm = 1 + θ (ln θ − 1) +
    /// θ TPD(w) with θ = ΣW, own rearrangement), so any negative value on
    /// the way proves a split. A trial stops when it returns to the feed or
    /// moves by less than 1e-10 in ln W, after at most 300 steps each.
    /// Components with zᵢ = 0 stay absent; a single component is stable.
    pub fn stability(&self, z: &[f64; MAX], t: f64, p: f64) -> Stability {
        let n = self.n;
        let feed = self.phase(z, t, p, Root::Stable);
        let mut d = [0.0; MAX];
        let mut active = 0;
        for i in 0..n {
            if z[i] > 0.0 {
                d[i] = math::ln(z[i]) + feed.ln_phi[i];
                active += 1;
            }
        }
        let mut best = Stability {
            tm: 0.0,
            w: *z,
            lighter: false,
            steps: 0,
        };
        if active < 2 {
            return best;
        }
        let lk = self.wilson(t, p);
        for sign in [1.0, -1.0] {
            let mut lw = [0.0; MAX];
            for i in 0..n {
                if z[i] > 0.0 {
                    lw[i] = math::ln(z[i]) + sign * lk[i];
                }
            }
            for _ in 0..STABILITY_STEPS {
                best.steps += 1;
                let mut wv = [0.0; MAX];
                let mut s = 0.0;
                for i in 0..n {
                    if z[i] > 0.0 {
                        wv[i] = math::exp(lw[i]);
                        s += wv[i];
                    }
                }
                let mut w = [0.0; MAX];
                for i in 0..n {
                    w[i] = wv[i] / s;
                }
                let ph = self.phase(&w, t, p, Root::Stable);
                let mut tm = 1.0;
                let mut dmax: f64 = 0.0;
                let mut back = 0.0;
                let mut next = lw;
                for i in 0..n {
                    if z[i] > 0.0 {
                        tm += wv[i] * (lw[i] + ph.ln_phi[i] - d[i] - 1.0);
                        next[i] = d[i] - ph.ln_phi[i];
                        dmax = dmax.max((next[i] - lw[i]).abs());
                        let r = math::ln(w[i] / z[i]);
                        back += r * r;
                    }
                }
                if tm < best.tm {
                    best.tm = tm;
                    best.w = w;
                    best.lighter = ph.v > feed.v;
                }
                if back < 1e-8 || dmax < 1e-10 {
                    break;
                }
                lw = next;
            }
        }
        best
    }

    /// The residual F(ln K) = ln K − (ln φ^l(x) − ln φ^v(y)) of a flash of
    /// the feed `z`, with β from Rachford-Rice and both phases on the root
    /// of the lower Gibbs energy.
    fn flash_residual(
        &self,
        z: &[f64; MAX],
        t: f64,
        p: f64,
        lnk: &[f64; MAX],
    ) -> Option<FlashResidual> {
        let n = self.n;
        let beta = rachford_rice(n, z, lnk)?;
        let mut x = [0.0; MAX];
        let mut y = [0.0; MAX];
        let mut sx = 0.0;
        let mut sy = 0.0;
        for i in 0..n {
            if z[i] > 0.0 {
                let k = math::exp(lnk[i]);
                x[i] = z[i] / (1.0 + beta * (k - 1.0));
                y[i] = k * x[i];
                sx += x[i];
                sy += y[i];
            }
        }
        for i in 0..n {
            x[i] /= sx;
            y[i] /= sy;
        }
        let l = self.phase(&x, t, p, Root::Stable);
        let v = self.phase(&y, t, p, Root::Stable);
        let mut f = [0.0; MAX];
        for i in 0..n {
            if z[i] > 0.0 {
                f[i] = lnk[i] - (l.ln_phi[i] - v.ln_phi[i]);
            }
        }
        Some((f, beta, x, y, l, v))
    }

    /// The isothermal flash of the feed `z` at `t` (K) and `p` (Pa): the
    /// stability test first; on a split, K from the trial phase (or from
    /// `start`, ln K of a nearby state), successive substitution ln K ←
    /// ln φ^l(x) − ln φ^v(y) with β from Rachford-Rice, finished by Newton's
    /// method on ln K with a difference Jacobian, each a bounded number of
    /// steps. When the K values collapse to 1 or β leaves 0 to 1 the feed
    /// is one phase.
    pub fn flash(&self, z: &[f64; MAX], t: f64, p: f64, start: Option<&[f64; MAX]>) -> Flash {
        let n = self.n;
        let one = || Flash::One(self.phase(z, t, p, Root::Stable));
        let st = self.stability(z, t, p);
        if !st.unstable() {
            return one();
        }
        let mut lnk = [0.0; MAX];
        if let Some(s) = start {
            lnk = *s;
        } else {
            for i in 0..n {
                if z[i] > 0.0 {
                    let r = math::ln(st.w[i] / z[i]);
                    lnk[i] = if st.lighter { r } else { -r };
                }
            }
        }
        let trivial = |lk: &[f64; MAX]| {
            let mut s = 0.0;
            for i in 0..n {
                if z[i] > 0.0 {
                    s += lk[i] * lk[i];
                }
            }
            s < TRIVIAL
        };
        let mut steps = st.steps;
        let mut done = false;
        for _ in 0..FLASH_SUBSTITUTIONS {
            steps += 1;
            let Some((f, ..)) = self.flash_residual(z, t, p, &lnk) else {
                return one();
            };
            let mut dmax: f64 = 0.0;
            for i in 0..n {
                lnk[i] -= f[i];
                dmax = dmax.max(f[i].abs());
            }
            if trivial(&lnk) {
                return one();
            }
            if dmax < TOL {
                done = true;
                break;
            }
        }
        let mut it = 0;
        while !done && it < NEWTON {
            it += 1;
            steps += 1;
            let Some((f, ..)) = self.flash_residual(z, t, p, &lnk) else {
                return Flash::NotConverged;
            };
            let mut jac = [[0.0; MAX]; MAX];
            for c in 0..n {
                if z[c] <= 0.0 {
                    jac[c][c] = 1.0;
                    continue;
                }
                let h = 1e-7 * (1.0 + lnk[c].abs());
                let mut lp = lnk;
                lp[c] += h;
                let Some((fp, ..)) = self.flash_residual(z, t, p, &lp) else {
                    return Flash::NotConverged;
                };
                for r in 0..n {
                    jac[r][c] = (fp[r] - f[r]) / h;
                }
            }
            let Some(dl) = solve(n, &mut jac, &f) else {
                return Flash::NotConverged;
            };
            let mut dmax: f64 = 0.0;
            for i in 0..n {
                lnk[i] += dl[i];
                dmax = dmax.max(dl[i].abs());
            }
            if trivial(&lnk) {
                return one();
            }
            if dmax < TOL {
                done = true;
            }
        }
        if !done {
            return Flash::NotConverged;
        }
        let Some((_, beta, x, y, l, v)) = self.flash_residual(z, t, p, &lnk) else {
            return one();
        };
        if !(beta > 0.0 && beta < 1.0) {
            return one();
        }
        // the lighter phase is y: swap when the roots say otherwise
        let s = if v.v >= l.v {
            Split {
                beta,
                x,
                y,
                ln_k: lnk,
                l,
                v,
                steps,
            }
        } else {
            let mut nk = [0.0; MAX];
            for i in 0..n {
                nk[i] = -lnk[i];
            }
            Split {
                beta: 1.0 - beta,
                x: y,
                y: x,
                ln_k: nk,
                l: v,
                v: l,
                steps,
            }
        };
        Flash::Two(s)
    }
}

/// How a mixing line treats the phases (D319, VA.4).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mixing {
    /// One phase of the equation of state at every mixture fraction, the
    /// root of the lower Gibbs energy (no phase split: "frozen").
    Frozen,
    /// Phase equilibrium of the equation of state (isobaric and
    /// isenthalpic flash).
    Equilibrium,
    /// Ideal gas and ideal solution: the liquid is pure fuel with the
    /// enthalpy h^ig − L_sat of the equation of state (L_sat = h_v − h_l at
    /// the pure fuel's saturation pressure), the vapour fraction of the fuel
    /// over it is p_sat/p (Raoult's law), no gas dissolves.
    Ideal,
}

/// A state of a mixing line.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MixingPoint {
    /// Mixture fraction Z: the mass fraction that comes from the fuel
    /// stream.
    pub z: f64,
    /// Temperature (K).
    pub t: f64,
    /// Specific enthalpy of the mixture, Z h_fuel + (1 − Z) h_gas (J/kg).
    pub h: f64,
    /// The enthalpy of the state found minus `h` (J/kg): the error of the
    /// energy balance.
    pub residual: f64,
    /// Mass fraction of the denser phase (0 for one phase).
    pub liquid: f64,
    /// Density of the mixture (kg/m³), both phases together.
    pub rho: f64,
    /// Two phases (equilibrium and ideal), or the one phase of the frozen
    /// line fails the stability test.
    pub split: bool,
}

/// The two streams of an adiabatic mixing line at one pressure: the fuel
/// (component 0, pure) and the ambient gas.
#[derive(Clone, Copy)]
pub struct Streams<'a> {
    /// The mixture of the fuel and the gases.
    pub mix: &'a Mixture,
    /// The ideal-gas molar enthalpy hᵢ^ig(T) of component i (J/mol), with
    /// the same reference for all components (the enthalpies of formation
    /// cancel in the balance; only differences of temperature enter).
    pub ideal: &'a dyn Fn(usize, f64) -> f64,
    /// Pressure (Pa).
    pub p: f64,
    /// Temperature of the fuel (K).
    pub t_fuel: f64,
    /// Temperature of the gas (K).
    pub t_gas: f64,
    /// Mole fractions of the ambient gas (component 0 zero).
    pub gas: [f64; MAX],
}

impl Streams<'_> {
    /// Molar enthalpy of a phase: Σ xᵢ hᵢ^ig(T) + h − h^ig.
    fn molar_h(&self, x: &[f64; MAX], t: f64, ph: &Phase) -> f64 {
        let mut h = ph.h_dep;
        for (i, &xi) in x.iter().enumerate().take(self.mix.n) {
            if xi > 0.0 {
                h += xi * (self.ideal)(i, t);
            }
        }
        h
    }

    /// Latent heat (J/mol) of the pure fuel at `t` at its saturation
    /// pressure and that pressure; `None` at or above its critical
    /// temperature.
    fn latent(&self, t: f64) -> Option<(f64, f64)> {
        let pure = Mixture::new(&[self.mix.c[0]]);
        let ps = pure.saturation_pressure(0, t)?;
        let mut x = [0.0; MAX];
        x[0] = 1.0;
        let l = pure.phase(&x, t, ps, Root::Liquid);
        let v = pure.phase(&x, t, ps, Root::Vapour);
        Some((v.h_dep - l.h_dep, ps))
    }

    /// Specific enthalpies (J/kg) of the fuel and the gas stream.
    pub fn stream_enthalpies(&self, model: Mixing) -> (f64, f64) {
        let m = self.mix;
        let mut f = [0.0; MAX];
        f[0] = 1.0;
        let mf = m.c[0].molar_mass;
        let mg = m.molar_mass(&self.gas);
        let hf_ig = (self.ideal)(0, self.t_fuel);
        let mut hg_ig = 0.0;
        for i in 1..m.n {
            if self.gas[i] > 0.0 {
                hg_ig += self.gas[i] * (self.ideal)(i, self.t_gas);
            }
        }
        match model {
            Mixing::Ideal => {
                let l = self.latent(self.t_fuel).map_or(0.0, |(l, _)| l);
                ((hf_ig - l) / mf, hg_ig / mg)
            }
            Mixing::Frozen | Mixing::Equilibrium => {
                let pf = m.phase(&f, self.t_fuel, self.p, Root::Stable);
                let pg = m.phase(&self.gas, self.t_gas, self.p, Root::Stable);
                ((hf_ig + pf.h_dep) / mf, (hg_ig + pg.h_dep) / mg)
            }
        }
    }

    /// Mole fractions of the mixture at the mixture fraction `z`.
    pub fn composition(&self, z: f64) -> [f64; MAX] {
        let m = self.mix;
        let mg = m.molar_mass(&self.gas);
        let nf = z / m.c[0].molar_mass;
        let ng = (1.0 - z) / mg;
        let mut x = [0.0; MAX];
        x[0] = nf / (nf + ng);
        for (xi, gi) in x.iter_mut().zip(self.gas.iter()).take(m.n).skip(1) {
            *xi = ng * gi / (nf + ng);
        }
        x
    }

    /// Specific enthalpy (J/kg), liquid mass fraction, density and the
    /// split flag of the mixture `x` at `t`.
    fn state(&self, x: &[f64; MAX], t: f64, model: Mixing) -> (f64, f64, f64, bool) {
        let m = self.mix;
        let mm = m.molar_mass(x);
        match model {
            Mixing::Frozen => {
                let ph = m.phase(x, t, self.p, Root::Stable);
                (self.molar_h(x, t, &ph) / mm, 0.0, mm / ph.v, false)
            }
            Mixing::Equilibrium => match m.flash(x, t, self.p, None) {
                Flash::Two(s) => {
                    let h = s.beta * self.molar_h(&s.y, t, &s.v)
                        + (1.0 - s.beta) * self.molar_h(&s.x, t, &s.l);
                    let ml = (1.0 - s.beta) * m.molar_mass(&s.x);
                    let v = s.beta * s.v.v + (1.0 - s.beta) * s.l.v;
                    (h / mm, ml / mm, mm / v, true)
                }
                Flash::One(ph) => (self.molar_h(x, t, &ph) / mm, 0.0, mm / ph.v, false),
                Flash::NotConverged => {
                    let ph = m.phase(x, t, self.p, Root::Stable);
                    (self.molar_h(x, t, &ph) / mm, 0.0, mm / ph.v, false)
                }
            },
            Mixing::Ideal => {
                let mut h = 0.0;
                for (i, &xi) in x.iter().enumerate().take(m.n) {
                    if xi > 0.0 {
                        h += xi * (self.ideal)(i, t);
                    }
                }
                let v_gas = R * t / self.p;
                let Some((l, ps)) = self.latent(t) else {
                    return (h / mm, 0.0, mm / v_gas, false);
                };
                let y0 = ps / self.p;
                if x[0] <= y0 || y0 >= 1.0 {
                    return (h / mm, 0.0, mm / v_gas, false);
                }
                // moles per mole of mixture: gas 1 − x₀, fuel vapour
                // y₀/(1 − y₀) of that, the rest of the fuel liquid
                let nv = y0 / (1.0 - y0) * (1.0 - x[0]);
                let nl = x[0] - nv;
                let mut one = [0.0; MAX];
                one[0] = 1.0;
                let pure = Mixture::new(&[m.c[0]]);
                let lv = pure.phase(&one, t, self.p, Root::Liquid).v;
                let v = (1.0 - nl) * v_gas + nl * lv;
                let ml = nl * m.c[0].molar_mass;
                ((h - nl * l) / mm, ml / mm, mm / v, true)
            }
        }
    }

    /// The state of the mixing line at the mixture fraction `z`: the
    /// temperature at which the enthalpy of the mixture equals Z h_fuel +
    /// (1 − Z) h_gas, by bisection between the stream temperatures (widened
    /// when needed) down to adjacent floating-point numbers, so the energy
    /// balance holds to round-off wherever h(T) is continuous.
    pub fn point(&self, z: f64, model: Mixing) -> MixingPoint {
        let (hf, hg) = self.stream_enthalpies(model);
        let target = z * hf + (1.0 - z) * hg;
        let x = self.composition(z);
        let f = |t: f64| self.state(&x, t, model).0 - target;
        let mut lo = self.t_fuel.min(self.t_gas);
        let mut hi = self.t_fuel.max(self.t_gas);
        for _ in 0..40 {
            if f(lo) <= 0.0 || lo <= 20.0 {
                break;
            }
            lo -= 10.0;
        }
        for _ in 0..40 {
            if f(hi) >= 0.0 {
                break;
            }
            hi += 10.0;
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi {
                break;
            }
            if f(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let t = if f(lo).abs() <= f(hi).abs() { lo } else { hi };
        let (h, liquid, rho, two) = self.state(&x, t, model);
        let split = match model {
            Mixing::Frozen => self.mix.stability(&x, t, self.p).unstable(),
            Mixing::Equilibrium | Mixing::Ideal => two,
        };
        MixingPoint {
            z,
            t,
            h: target,
            residual: h - target,
            liquid,
            rho,
            split,
        }
    }

    /// The smallest mixture fraction at which the line splits (two phases;
    /// for the frozen line, a phase that fails the stability test): the
    /// saturation mixture fraction Z_sat of a mixing-limited liquid length
    /// (Siebers 1999, B = Z_sat/(1 − Z_sat)). A scan in steps of 0.02 and
    /// 40 bisections; `None` if the line does not split.
    pub fn saturation(&self, model: Mixing) -> Option<f64> {
        let mut prev = 0.0;
        for k in 1..50 {
            let z = f64::from(k) / 50.0;
            if self.point(z, model).split {
                let (mut lo, mut hi) = (prev, z);
                for _ in 0..40 {
                    let mid = 0.5 * (lo + hi);
                    if self.point(mid, model).split {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                return Some(0.5 * (lo + hi));
            }
            prev = z;
        }
        None
    }
}

/// The coefficients a, b, c, d of E₁ … E₁₀ of the dense-fluid viscosity of
/// Chung et al. (1988), Eᵢ = aᵢ + bᵢ ω + cᵢ μ_r⁴ + dᵢ κ (D313; see
/// `docs/verification/DATA_SOURCES.md`, section 3h).
const CHUNG_E: [[f64; 4]; 10] = [
    [6.324_02, 50.411_9, -51.680_1, 1_189.02],
    [0.121_02e-2, -0.115_36e-2, -0.625_71e-2, 0.037_283],
    [5.283_46, 254.209, -168.481, 3_898.27],
    [6.622_63, 38.095_7, -8.464_14, 31.417_8],
    [19.745_4, 7.630_34, -14.354_4, 31.526_7],
    [-1.899_92, -12.536_7, 4.985_29, -18.150_7],
    [24.274_5, 3.449_45, -11.291_3, 69.346_6],
    [0.797_16, 1.117_64, 0.012_348, -4.116_61],
    [-0.238_16, 0.067_695, -0.816_3, 4.025_28],
    [0.068_629, 0.347_93, 0.592_56, -0.726_63],
];

/// The coefficients of B₁ … B₇ of the dense-fluid thermal conductivity of
/// Chung et al. (1988), as [`CHUNG_E`].
const CHUNG_B: [[f64; 4]; 7] = [
    [2.416_6, 0.748_24, -0.918_58, 121.72],
    [-0.509_24, -1.509_4, -49.991, 69.983],
    [6.610_7, 5.620_7, 64.760, 27.039],
    [14.543, -8.913_9, -5.637_9, 74.344],
    [0.792_74, 0.820_19, -0.693_69, 6.317_3],
    [-5.863_4, 12.801, 9.589_3, 65.529],
    [91.089, 128.11, -54.217, 523.81],
];

/// The reduced collision integral Ω_v(T*) of Neufeld, Janzen and Aziz
/// (1972) with its correction term, as Chung et al. (1988) use it.
#[inline]
pub fn neufeld_omega_v(ts: f64) -> f64 {
    let pw = math::exp(-0.148_74 * math::ln(ts));
    1.161_45 * pw + 0.524_87 * math::exp(-0.773_20 * ts) + 2.161_78 * math::exp(-2.437_87 * ts)
        - 6.435e-4 / pw * math::sin(18.032_3 * math::exp(-0.768_3 * math::ln(ts)) - 7.273_71)
}

/// The pseudo-pure fluid of a mixture by the mixing rules of Chung et al.
/// (1988): σ_m³ = Σᵢ Σⱼ xᵢ xⱼ σᵢⱼ³, (ε/k)_m σ_m³ = Σ Σ xᵢ xⱼ (ε/k)ᵢⱼ σᵢⱼ³,
/// ω_m σ_m³ = Σ Σ xᵢ xⱼ ωᵢⱼ σᵢⱼ³, M_m = [Σ Σ xᵢ xⱼ (ε/k)ᵢⱼ σᵢⱼ² √Mᵢⱼ /
/// ((ε/k)_m σ_m²)]², μ_m⁴ = σ_m³ Σ Σ xᵢ xⱼ μᵢ² μⱼ²/σᵢⱼ³, κ_m = Σ Σ xᵢ xⱼ κᵢⱼ
/// with σᵢ = 0.809 V_c,ᵢ^⅓, (ε/k)ᵢ = T_c,ᵢ/1.2593, σᵢⱼ = √(σᵢ σⱼ),
/// (ε/k)ᵢⱼ = √((ε/k)ᵢ (ε/k)ⱼ), ωᵢⱼ = (ωᵢ + ωⱼ)/2, Mᵢⱼ = 2 Mᵢ Mⱼ/(Mᵢ + Mⱼ),
/// κᵢⱼ = √(κᵢ κⱼ); then T_c = 1.2593 (ε/k)_m, V_c = (σ_m/0.809)³.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChungFluid {
    /// T_c (K).
    pub tc: f64,
    /// V_c (cm³/mol, the unit of the correlation).
    pub vc: f64,
    /// ω.
    pub omega: f64,
    /// M (g/mol).
    pub m: f64,
    /// The reduced dipole moment μ_r = 131.3 μ/√(V_c T_c).
    pub mur: f64,
    /// κ.
    pub kappa: f64,
}

impl ChungFluid {
    /// T* = 1.2593 T/T_c, F_c = 1 − 0.2756 ω + 0.059035 μ_r⁴ + κ and Ω_v.
    #[inline]
    fn reduced(&self, t: f64) -> (f64, f64, f64) {
        let ts = 1.2593 * t / self.tc;
        let fc = 1.0 - 0.2756 * self.omega
            + 0.059_035 * self.mur * self.mur * self.mur * self.mur
            + self.kappa;
        (ts, fc, neufeld_omega_v(ts))
    }

    /// The coefficients Eᵢ or Bᵢ of this fluid.
    #[inline]
    fn coefficients<const N: usize>(&self, tab: &[[f64; 4]; N]) -> [f64; N] {
        let m4 = self.mur * self.mur * self.mur * self.mur;
        let mut e = [0.0; N];
        for (ei, r) in e.iter_mut().zip(tab.iter()) {
            *ei = r[0] + r[1] * self.omega + r[2] * m4 + r[3] * self.kappa;
        }
        e
    }

    /// y = ρ V_c/6, G₁ = (1 − y/2)/(1 − y)³ and G₂ (or H₂) of the table `c`.
    #[inline]
    fn g2(y: f64, c: &[f64]) -> (f64, f64) {
        let g1 = (1.0 - 0.5 * y) / ((1.0 - y) * (1.0 - y) * (1.0 - y));
        // (1 − e^(−c₄ y))/y, by its series below y = 1e-8
        let q = if y > 1e-8 {
            (1.0 - math::exp(-c[3] * y)) / y
        } else {
            c[3] * (1.0 - 0.5 * c[3] * y)
        };
        let g2 =
            (c[0] * q + c[1] * g1 * math::exp(c[4] * y) + c[2] * g1) / (c[0] * c[3] + c[1] + c[2]);
        (g1, g2)
    }

    /// The low-pressure viscosity η° (Pa s): 40.785 F_c √(M T)/(V_c^⅔ Ω_v)
    /// µP (Chung et al. 1984, 1988).
    pub fn viscosity0(&self, t: f64) -> f64 {
        let (_, fc, ov) = self.reduced(t);
        40.785e-7 * fc * math::sqrt(self.m * t) / (cbrt(self.vc) * cbrt(self.vc) * ov)
    }

    /// The viscosity (Pa s) at `t` and the molar density `rho` (mol/m³):
    /// η = η* 36.344 √(M T_c)/V_c^⅔ µP, η* = (√T*/Ω_v) F_c (1/G₂ + E₆ y) +
    /// η**, η** = E₇ y² G₂ exp(E₈ + E₉/T* + E₁₀/T*²).
    pub fn viscosity(&self, t: f64, rho: f64) -> f64 {
        self.viscosity0(t) * self.viscosity_factor(t, rho)
    }

    /// η/η° at `t` and `rho` (mol/m³): (1/G₂ + E₆ y) + η** Ω_v/(√T* F_c).
    pub fn viscosity_factor(&self, t: f64, rho: f64) -> f64 {
        let (ts, fc, ov) = self.reduced(t);
        let e = self.coefficients(&CHUNG_E);
        let y = rho * self.vc * 1e-6 / 6.0;
        let (_, g2) = Self::g2(y, &e);
        let ess = e[6] * y * y * g2 * math::exp(e[7] + e[8] / ts + e[9] / (ts * ts));
        1.0 / g2 + e[5] * y + ess * ov / (math::sqrt(ts) * fc)
    }

    /// The low-pressure thermal conductivity λ° (W/(m K)) with the
    /// ideal-gas molar heat capacity `cv` (J/(mol K)): λ° M'/(η° c_v) =
    /// 3.75 Ψ/(c_v/R), that is λ° = 31.2 η° Ψ/M' (η° in Pa s, M' in kg/mol;
    /// Chung et al. 1984), Ψ = 1 + α (0.215 + 0.28288 α − 1.061 β +
    /// 0.26665 Z)/(0.6366 + β Z + 1.061 α β), α = c_v/R − 3/2, β = 0.7862 −
    /// 0.7109 ω + 1.3168 ω², Z = 2 + 10.5 T_r².
    pub fn conductivity0(&self, t: f64, cv: f64) -> f64 {
        let tr = t / self.tc;
        let a = cv / R - 1.5;
        let b = 0.7862 - 0.7109 * self.omega + 1.3168 * self.omega * self.omega;
        let z = 2.0 + 10.5 * tr * tr;
        // 0.6366 is Chung's printed constant, not 2/π
        #[allow(clippy::approx_constant)]
        let psi = 1.0
            + a * (0.215 + 0.282_88 * a - 1.061 * b + 0.266_65 * z)
                / (0.6366 + b * z + 1.061 * a * b);
        3.75 * R * self.viscosity0(t) * psi / (self.m * 1e-3)
    }

    /// The thermal conductivity (W/(m K)) at `t`, `rho` (mol/m³) with the
    /// ideal-gas `cv` (J/(mol K)): λ = λ° (1/H₂ + B₆ y) + q B₇ y² √T_r H₂,
    /// q = 3.586e-3 √(T_c/M')/V_c^⅔.
    pub fn conductivity(&self, t: f64, rho: f64, cv: f64) -> f64 {
        self.conductivity0(t, cv) * self.conductivity_factor(t, rho, cv)
    }

    /// λ/λ° at `t`, `rho` (mol/m³) and the ideal-gas `cv` (J/(mol K)).
    pub fn conductivity_factor(&self, t: f64, rho: f64, cv: f64) -> f64 {
        let b = self.coefficients(&CHUNG_B);
        let y = rho * self.vc * 1e-6 / 6.0;
        let (_, h2) = Self::g2(y, &b);
        let q = 3.586e-3 * math::sqrt(self.tc / (self.m * 1e-3)) / (cbrt(self.vc) * cbrt(self.vc));
        let lp = q * b[6] * y * y * math::sqrt(t / self.tc) * h2;
        1.0 / h2 + b[5] * y + lp / self.conductivity0(t, cv)
    }
}

impl Mixture {
    /// The pseudo-pure fluid of Chung et al. (1988) of the composition `x`
    /// (see [`ChungFluid`]).
    pub fn chung(&self, x: &[f64; MAX]) -> ChungFluid {
        let n = self.n;
        let mut sig = [0.0; MAX];
        let mut eps = [0.0; MAX];
        for i in 0..n {
            sig[i] = 0.809 * cbrt(self.c[i].vc() * 1e6);
            eps[i] = self.c[i].tc / 1.2593;
        }
        let (mut s3, mut e3, mut w3, mut mm, mut d4, mut k) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for i in 0..n {
            for j in 0..n {
                let xx = x[i] * x[j];
                let sij = math::sqrt(sig[i] * sig[j]);
                let sij3 = sij * sij * sij;
                let eij = math::sqrt(eps[i] * eps[j]);
                let (ci, cj) = (&self.c[i], &self.c[j]);
                let mij =
                    2.0 * ci.molar_mass * cj.molar_mass / (ci.molar_mass + cj.molar_mass) * 1e3;
                s3 += xx * sij3;
                e3 += xx * eij * sij3;
                w3 += xx * 0.5 * (ci.omega + cj.omega) * sij3;
                mm += xx * eij * sij * sij * math::sqrt(mij);
                d4 += xx * ci.dipole * ci.dipole * cj.dipole * cj.dipole / sij3;
                k += xx * math::sqrt(ci.assoc * cj.assoc);
            }
        }
        let sm = cbrt(s3);
        let em = e3 / s3;
        let tc = 1.2593 * em;
        let vc = (sm / 0.809) * (sm / 0.809) * (sm / 0.809);
        let m = (mm / (em * sm * sm)) * (mm / (em * sm * sm));
        let dm = math::sqrt(math::sqrt(s3 * d4));
        ChungFluid {
            tc,
            vc,
            omega: w3 / s3,
            m,
            mur: 131.3 * dm / math::sqrt(vc * tc),
            kappa: k,
        }
    }

    /// The factor on the low-pressure ρD of a dilute solute in this mixture
    /// at `p` (Pa) after Riazi and Whitson (1993), ρD/(ρD)° = 1.07 (η/η°)^n,
    /// n = −0.27 − 0.38 ω + (−0.05 + 0.1 ω) p_r with ω = Σ xᵢ ωᵢ and
    /// p_r = p/Σ xᵢ p_c,ᵢ, given the dense-fluid factor `eta_ratio` = η/η°;
    /// here without their constant 1.07, so that the factor is 1 at low
    /// pressure and the film keeps the diffusivity of the low-pressure
    /// model there (D313).
    pub fn riazi_whitson(&self, x: &[f64; MAX], p: f64, eta_ratio: f64) -> f64 {
        let (mut w, mut pc) = (0.0, 0.0);
        for (xi, c) in x.iter().zip(&self.c[..self.n]) {
            w += xi * c.omega;
            pc += xi * c.pc;
        }
        let n = -0.27 - 0.38 * w + (-0.05 + 0.1 * w) * p / pc;
        math::exp(n * math::ln(eta_ratio))
    }

    /// The molar heat capacity departure c_p − c_p^ig (J/(mol K)) of a phase
    /// of composition `x` at `t` and `p`: ∂(h − h^ig)/∂T at constant p and
    /// x by a centred difference of 0.01 K.
    pub fn cp_departure(&self, x: &[f64; MAX], t: f64, p: f64, root: Root) -> f64 {
        let h = 0.01;
        (self.phase(x, t + h, p, root).h_dep - self.phase(x, t - h, p, root).h_dep) / (2.0 * h)
    }
}

/// Critical constants and acentric factors of the gases of the ambient
/// (species names of the flow's species data). Critical points from the
/// reference equations of state; acentric factors from the tables of
/// Poling, Prausnitz and O'Connell (2001), appendix A (see
/// `docs/verification/DATA_SOURCES.md`, section 3e). Critical volumes
/// (the inverse critical densities), the dipole moment of water from the
/// NIST Chemistry WebBook fluid tables, the association factor of water
/// from Chung et al. (1988) (D313, section 3h).
pub mod gases {
    use super::Component;

    /// N₂: Span, Lemmon, Jacobsen, Wagner, Yokozeki (2000).
    pub const N2: Component = Component {
        vc: 1.0 / 11_183.9,
        ..Component::new(126.192, 3.3958e6, 0.0372, 0.028_013_4)
    };
    /// O₂: Schmidt and Wagner (1985).
    pub const O2: Component = Component {
        vc: 1.0 / 13_630.0,
        ..Component::new(154.581, 5.043e6, 0.0222, 0.031_998_8)
    };
    /// CO₂: Span and Wagner (1996).
    pub const CO2: Component = Component {
        vc: 1.0 / 10_624.9,
        ..Component::new(304.1282, 7.3773e6, 0.2239, 0.044_009_5)
    };
    /// H₂O: Wagner and Pruss (2002).
    pub const H2O: Component = Component {
        vc: 1.0 / 17_873.728,
        dipole: 1.855,
        assoc: 0.076,
        ..Component::new(647.096, 22.064e6, 0.3443, 0.018_015_268)
    };
    /// Ar: Tegeler, Span, Wagner (1999).
    pub const AR: Component = Component {
        vc: 1.0 / 13_407.4,
        ..Component::new(150.687, 4.863e6, -0.0022, 0.039_948)
    };

    /// The data of a species by its name in the species data, if known.
    pub fn by_name(name: &str) -> Option<Component> {
        match name {
            "N2" => Some(N2),
            "O2" => Some(O2),
            "CO2" => Some(CO2),
            "H2O" => Some(H2O),
            "AR" | "Ar" => Some(AR),
            _ => None,
        }
    }
}

/// High-pressure data of the fuels: critical point and acentric factor
/// with the vapour pressure equation they were published with, the
/// Mathias-Copeman coefficients that make the Peng-Robinson saturation
/// pressure follow that equation (own fit, reproduced by the test
/// `fit_mathias_copeman_reproduces_the_dodecane_coefficients` of
/// `foehncraft-particles`), and the binary interaction parameters with the
/// ambient gases fitted to measured phase equilibria (see
/// `docs/verification/DATA_SOURCES.md`, section 3e).
pub mod fuels {
    use super::Component;
    use crate::species::Species;

    /// A fuel of the high-pressure equilibrium.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct Fuel {
        /// Name as the liquid section gives it.
        pub name: &'static str,
        /// Critical point, acentric factor, molar mass and α coefficients.
        pub component: Component,
        /// k_ij with the gases of `super::gases`, by species name; a gas
        /// that is not listed takes k_ij = 0 unless the case sets one.
        pub k: &'static [(&'static str, f64)],
        /// The vapour as an ideal gas: NASA-7 polynomials and
        /// Lennard-Jones data, for the enthalpies of a mixing line.
        pub ideal_gas: Species,
    }

    /// n-dodecane: T_c = 658.65 K, p_c = 1.83 MPa and ω = 0.573 of Ambrose
    /// and Walton (1989), whose Wagner equation (their eq. 3, table 1) the
    /// Mathias-Copeman α follows within 1 % from 300 K to 655 K;
    /// k_ij(N₂) = 0.205 fitted to the bubble points of García-Córdova,
    /// Justo-García, García-Flores and García-Sánchez (2011) up to 20 MPa
    /// (344 K to 594 K; mean deviation 2.7 % in pressure). The ideal gas:
    /// c_p of Scott (1974) by his CH₂ increment from n-nonane and n-decane
    /// (NIST WebBook), within 0.11 % of the reference equation of Lemmon
    /// and Huber (2004) at 300 to 700 K; the NASA-7 fit (within 0.28 %),
    /// a₆ from Δ_f H° = −290.9 kJ/mol, a₇ from S° = 622.50 J/(mol K)
    /// (NIST WebBook), Lennard-Jones data by Tee, Gotoh and Stewart (1966):
    /// written by `tools/data/dodecane_ideal_gas.py`. Critical volume
    /// 1/(1.33 mol/l), the critical density of Lemmon and Huber (2004) (NIST
    /// WebBook), for the dense-fluid transport (D313).
    pub const N_DODECANE: Fuel = Fuel {
        name: "n-dodecane",
        component: Component {
            mc: [1.208_668, -0.307_470, 0.885_193],
            vc: 1.0 / 1_330.0,
            ..Component::new(658.65, 1.83e6, 0.573, 0.170_334)
        },
        k: &[("N2", 0.205)],
        ideal_gas: Species {
            molar_mass: 0.170_334,
            t_range: [273.15, 1000.0, 1500.0],
            low: [
                8.188_509_377_2e+00,
                6.192_592_395_7e-02,
                1.344_050_375_8e-04,
                -2.169_816_522_9e-07,
                8.846_824_494_5e-11,
                -4.098_147_709_6e+04,
                5.519_800_597_5e+00,
            ],
            high: [
                -1.754_171_898_3e+02,
                7.161_773_054_2e-01,
                -7.908_492_430_7e-04,
                4.051_955_614_1e-07,
                -7.910_037_037_7e-11,
                1.886_044_560_8e+03,
                9.166_985_466_4e+02,
            ],
            well_depth: 585.22,
            diameter: 7.6447,
            dipole: 0.0,
            geometry: 2,
        },
    };

    /// The fuel of a liquid name, if its high-pressure data are known.
    pub fn by_name(name: &str) -> Option<Fuel> {
        match name {
            "n-dodecane" | "dodecane" | "NC12H26" | "C12H26" => Some(N_DODECANE),
            _ => None,
        }
    }

    impl Fuel {
        /// k_ij of the fuel with the gas `gas` (species name), 0 if unknown.
        pub fn k_with(&self, gas: &str) -> f64 {
            self.k
                .iter()
                .find(|(g, _)| *g == gas)
                .map_or(0.0, |(_, k)| *k)
        }
    }
}

#[cfg(test)]
mod tests {
    extern crate std;
    use super::*;

    /// n-dodecane, critical point and acentric factor as quoted by Lemmon
    /// and Huber (2004) (test data only; the solver takes the liquid
    /// library's values).
    const C12: Component = Component::new(658.1, 1.817e6, 0.574, 0.170_334);

    fn binary(k: f64) -> Mixture {
        let mut m = Mixture::new(&[C12, gases::N2]);
        m.set_k(0, 1, k);
        m
    }

    #[test]
    fn the_omegas_are_the_roots_of_the_critical_conditions() {
        // η_c = 1/(1 + ∛(4 − √8) + ∛(4 + √8)) = Ω_b/Z_c ... checked through
        // the cubic: at T_c and p_c the cubic has the triple root Z_c
        let c = gases::N2;
        let m = Mixture::new(&[c]);
        let mut x = [0.0; MAX];
        x[0] = 1.0;
        let a = OMEGA_A;
        let b = OMEGA_B;
        // Z³ − (1 − B)Z² + (A − 3B² − 2B)Z − (AB − B² − B³) = (Z − Z_c)³
        let zc = Z_CRITICAL;
        assert!((1.0 - b - 3.0 * zc).abs() < 1e-15);
        assert!((a - 3.0 * b * b - 2.0 * b - 3.0 * zc * zc).abs() < 1e-15);
        assert!((a * b - b * b - b * b * b - zc * zc * zc).abs() < 1e-15);
        let ph = m.phase(&x, c.tc, c.pc, Root::Stable);
        assert!((ph.z - zc).abs() < 2e-5, "{}", ph.z);
    }

    #[test]
    fn cubic_roots_are_exact_for_a_known_polynomial() {
        // (Z − 0.1)(Z − 0.5)(Z − 0.9)
        let (r, n) = cubic_roots(-1.5, 0.59, -0.045);
        assert_eq!(n, 3);
        for (a, b) in r.iter().zip([0.1, 0.5, 0.9]) {
            assert!((a - b).abs() < 1e-14, "{a} {b}");
        }
        // Z³ − 1 = 0 has one real root
        let (r, n) = cubic_roots(0.0, 0.0, -1.0);
        assert_eq!(n, 1);
        assert!((r[0] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn fugacity_coefficients_agree_with_the_derivative_of_the_gibbs_energy() {
        // ln φᵢ = ∂(n g_dep/RT)/∂nᵢ: check by a difference in the mole
        // numbers at fixed T, p (the gas root, one phase)
        let m = binary(0.1561);
        let (t, p) = (700.0, 6e6);
        let n0 = [0.3, 0.7];
        let g = |nn: [f64; 2]| {
            let tot = nn[0] + nn[1];
            let mut x = [0.0; MAX];
            x[0] = nn[0] / tot;
            x[1] = nn[1] / tot;
            let ph = m.phase(&x, t, p, Root::Vapour);
            // n g_dep/RT = n Σ xᵢ ln φᵢ
            tot * (x[0] * ph.ln_phi[0] + x[1] * ph.ln_phi[1])
        };
        let mut x = [0.0; MAX];
        x[0] = 0.3;
        x[1] = 0.7;
        let ph = m.phase(&x, t, p, Root::Vapour);
        for i in 0..2 {
            let h = 1e-6;
            let mut up = n0;
            up[i] += h;
            let mut dn = n0;
            dn[i] -= h;
            let d = (g(up) - g(dn)) / (2.0 * h);
            assert!((d - ph.ln_phi[i]).abs() < 1e-8, "{i}: {d} {}", ph.ln_phi[i]);
        }
        // and the mixture: Σ xᵢ h̄ᵢ = h_dep
        let hb = m.partial_h_dep(&x, t, p, Root::Vapour);
        let sum = x[0] * hb[0] + x[1] * hb[1];
        assert!(
            (sum - ph.h_dep).abs() < 1e-6 * ph.h_dep.abs().max(1.0),
            "{sum} {}",
            ph.h_dep
        );
    }

    #[test]
    fn the_pure_saturation_pressure_has_equal_fugacities() {
        let m = Mixture::new(&[C12]);
        let ps = m.saturation_pressure(0, 500.0).unwrap();
        let mut x = [0.0; MAX];
        x[0] = 1.0;
        let l = m.phase(&x, 500.0, ps, Root::Liquid);
        let v = m.phase(&x, 500.0, ps, Root::Vapour);
        assert!((l.ln_phi[0] - v.ln_phi[0]).abs() < 1e-11);
        assert!(l.z < 0.1 && v.z > 0.8, "{} {}", l.z, v.z);
        assert!(m.saturation_pressure(0, 700.0).is_none());
    }

    #[test]
    fn the_surface_equilibrium_tends_to_raoult_at_low_pressure() {
        // VA.3: at 0.1 MPa the gas is nearly ideal and dissolves little;
        // y₀ p_sat⁻¹ p − 1 (the Poynting factor and the gas non-ideality)
        // is a few per cent and falls with p
        let m = binary(0.1561);
        let mut r = [0.0; MAX];
        r[1] = 1.0;
        for t in [320.0, 380.0, 440.0] {
            let ps = Mixture::new(&[C12]).saturation_pressure(0, t).unwrap();
            let enh = |p: f64| {
                let SurfaceState::TwoPhase(s) = surface(&m, t, p, &r, None) else {
                    panic!("no two-phase solution at {t} K, {p} Pa");
                };
                assert!(s.x[1] < 0.01, "{}", s.x[1]);
                s.y[0] / (ps / p) - 1.0
            };
            let e1 = enh(1e5);
            let e2 = enh(0.5e5);
            assert!(e1.abs() < 0.03, "{t} K: {e1}");
            assert!(e2.abs() < e1.abs(), "{t} K: {e1} {e2}");
        }
    }

    #[test]
    fn the_surface_equilibrium_at_spray_a_pressure_is_above_raoult() {
        // own check of D311: at 6 MPa the fuel vapour fraction is 1.5 to 2
        // times Raoult's value and N₂ dissolves (section 2 of the page)
        let m = binary(0.1561);
        let mut r = [0.0; MAX];
        r[1] = 1.0;
        let mut start = None;
        let mut last = [0.0; MAX];
        for t in [400.0, 450.0, 500.0, 550.0, 600.0, 620.0, 640.0] {
            let st = surface(&m, t, 6e6, &r, start.as_ref());
            let SurfaceState::TwoPhase(s) = st else {
                panic!("{t} K: {st:?}");
            };
            // the fugacities agree
            let l = m.phase(&s.x, t, 6e6, Root::Liquid);
            let v = m.phase(&s.y, t, 6e6, Root::Vapour);
            for i in 0..2 {
                let fl = l.ln_phi[i] + s.x[i].ln();
                let fv = v.ln_phi[i] + s.y[i].ln();
                assert!((fl - fv).abs() < 1e-10, "{t} K, {i}: {fl} {fv}");
            }
            let ps = Mixture::new(&[C12]).saturation_pressure(0, t).unwrap();
            let e = s.y[0] / (ps / 6e6);
            assert!((1.4..2.2).contains(&e), "{t} K: enhancement {e}");
            assert!(s.x[1] > 0.05 && s.x[1] < 0.3, "{t} K: x_N2 {}", s.x[1]);
            last = s.ln_k;
            start = Some(last);
        }
        // the same solution from Wilson's start
        let SurfaceState::TwoPhase(w) = surface(&m, 550.0, 6e6, &r, None) else {
            panic!()
        };
        let SurfaceState::TwoPhase(c) = surface(&m, 550.0, 6e6, &r, Some(&last)) else {
            panic!()
        };
        assert!((w.y[0] - c.y[0]).abs() < 1e-10);
    }

    #[test]
    fn far_above_the_mixing_critical_point_there_is_no_surface() {
        let m = binary(0.1561);
        let mut r = [0.0; MAX];
        r[1] = 1.0;
        let st = surface(&m, 700.0, 6e6, &r, None);
        assert!(
            matches!(
                st,
                SurfaceState::Critical | SurfaceState::Boiling | SurfaceState::NotConverged
            ),
            "{st:?}"
        );
        // below the fuel's critical pressure the droplet boils instead
        let st = surface(&m, 700.0, 1e6, &r, None);
        assert!(!matches!(st, SurfaceState::TwoPhase(_)), "{st:?}");
    }

    #[test]
    fn rachford_rice_finds_the_vapour_fraction_of_a_known_split() {
        // K = (2, 1/2), z = (1/2, 1/2): 1/(1 + β) = 1/(2 − β), β = 1/2
        let mut z = [0.0; MAX];
        z[0] = 0.5;
        z[1] = 0.5;
        let mut lk = [0.0; MAX];
        lk[0] = 2.0f64.ln();
        lk[1] = 0.5f64.ln();
        let b = rachford_rice(2, &z, &lk).unwrap();
        assert!((b - 0.5).abs() < 1e-15, "{b}");
        // all K above 1: no split
        lk[1] = 0.1;
        assert!(rachford_rice(2, &z, &lk).is_none());
    }

    #[test]
    fn the_pressure_is_the_inverse_of_the_phase_volume() {
        let m = binary(0.205);
        let mut x = [0.0; MAX];
        x[0] = 0.3;
        x[1] = 0.7;
        for (t, p, root) in [(700.0, 6e6, Root::Vapour), (400.0, 6e6, Root::Liquid)] {
            let ph = m.phase(&x, t, p, root);
            let back = m.pressure(&x, t, ph.v);
            assert!((back / p - 1.0).abs() < 1e-9, "{t} K: {back} {p}");
        }
    }

    #[test]
    fn the_stability_test_tells_one_phase_from_two() {
        let m = binary(0.205);
        let mut z = [0.0; MAX];
        // pure fuel and pure nitrogen: one component, stable
        z[0] = 1.0;
        assert!(!m.stability(&z, 450.0, 6e6).unstable());
        // hot nitrogen with a little fuel vapour: stable
        z[0] = 0.02;
        z[1] = 0.98;
        assert!(!m.stability(&z, 800.0, 6e6).unstable());
        // half fuel, half nitrogen at 450 K: liquid and gas
        z[0] = 0.5;
        z[1] = 0.5;
        let st = m.stability(&z, 450.0, 6e6);
        assert!(st.unstable(), "{st:?}");
    }

    #[test]
    fn a_binary_flash_gives_the_phases_of_the_surface_equilibrium() {
        // two components at fixed T and p: the coexisting compositions do
        // not depend on the feed, and they are those of the droplet
        // surface with nitrogen as the only gas
        let m = binary(0.205);
        let (t, p) = (500.0, 6e6);
        let mut r = [0.0; MAX];
        r[1] = 1.0;
        let SurfaceState::TwoPhase(s) = surface(&m, t, p, &r, None) else {
            panic!()
        };
        for z0 in [0.2, 0.5, 0.8] {
            let mut z = [0.0; MAX];
            z[0] = z0;
            z[1] = 1.0 - z0;
            let Flash::Two(f) = m.flash(&z, t, p, None) else {
                panic!("{z0}: no split");
            };
            assert!(
                (f.x[0] - s.x[0]).abs() < 1e-9,
                "{z0}: x {} {}",
                f.x[0],
                s.x[0]
            );
            assert!(
                (f.y[0] - s.y[0]).abs() < 1e-9,
                "{z0}: y {} {}",
                f.y[0],
                s.y[0]
            );
            // mass balance and the denser phase second
            for (i, zi) in z.iter().enumerate().take(2) {
                let back = f.beta * f.y[i] + (1.0 - f.beta) * f.x[i];
                assert!((back - zi).abs() < 1e-12, "{z0}, {i}: {back}");
            }
            assert!(f.l.v < f.v.v);
        }
        // a feed outside the two-phase range is one phase
        let mut z = [0.0; MAX];
        z[0] = 0.5 * s.y[0];
        z[1] = 1.0 - z[0];
        assert!(matches!(m.flash(&z, t, p, None), Flash::One(_)));
    }
}
