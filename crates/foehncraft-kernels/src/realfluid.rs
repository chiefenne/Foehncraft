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
/// Convergence of ln K.
const TOL: f64 = 1e-12;
/// ‖ln K‖² below which the solution is taken as trivial (x = y).
const TRIVIAL: f64 = 1e-10;

/// F, the liquid and vapour phases and their compositions x and y.
type Residual = ([f64; MAX], Phase, Phase, [f64; MAX], [f64; MAX]);

/// The residual F(ln K) = ln K − (ln φ^l(x) − ln φ^v(y)) and the phases.
#[inline]
fn residual(
    mix: &Mixture,
    t: f64,
    p: f64,
    lnk: &[f64; MAX],
    r: &[f64; MAX],
) -> Option<Residual> {
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
        if dmax < TOL {
            done = true;
            break;
        }
    }
    // Newton with a difference Jacobian
    let mut it = 0;
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
        if dmax < TOL {
            done = true;
        }
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

/// Critical constants and acentric factors of the gases of the ambient
/// (species names of the flow's species data). Critical points from the
/// reference equations of state; acentric factors from the tables of
/// Poling, Prausnitz and O'Connell (2001), appendix A (see
/// `docs/verification/DATA_SOURCES.md`, section 3e).
pub mod gases {
    use super::Component;

    /// N₂: Span, Lemmon, Jacobsen, Wagner, Yokozeki (2000).
    pub const N2: Component = Component::new(126.192, 3.3958e6, 0.0372, 0.028_013_4);
    /// O₂: Schmidt and Wagner (1985).
    pub const O2: Component = Component::new(154.581, 5.043e6, 0.0222, 0.031_998_8);
    /// CO₂: Span and Wagner (1996).
    pub const CO2: Component = Component::new(304.1282, 7.3773e6, 0.2239, 0.044_009_5);
    /// H₂O: Wagner and Pruss (2002).
    pub const H2O: Component = Component::new(647.096, 22.064e6, 0.3443, 0.018_015_268);
    /// Ar: Tegeler, Span, Wagner (1999).
    pub const AR: Component = Component::new(150.687, 4.863e6, -0.0022, 0.039_948);

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
    }

    /// n-dodecane: T_c = 658.65 K, p_c = 1.83 MPa and ω = 0.573 of Ambrose
    /// and Walton (1989), whose Wagner equation (their eq. 3, table 1) the
    /// Mathias-Copeman α follows within 1 % from 300 K to 655 K;
    /// k_ij(N₂) = 0.205 fitted to the bubble points of García-Córdova,
    /// Justo-García, García-Flores and García-Sánchez (2011) up to 20 MPa
    /// (344 K to 594 K; mean deviation 2.7 % in pressure).
    pub const N_DODECANE: Fuel = Fuel {
        name: "n-dodecane",
        component: Component {
            mc: [1.208_668, -0.307_470, 0.885_193],
            ..Component::new(658.65, 1.83e6, 0.573, 0.170_334)
        },
        k: &[("N2", 0.205)],
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
}
