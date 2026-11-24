//! Water and steam states of IAPWS-IF97 (stage SF7a of the steam-furnace
//! path, D94, D102): one body for the CPU and the GPU.
//!
//! - **Thermodynamics:** the Industrial Formulation 1997 for the
//!   Thermodynamic Properties of Water and Steam (IAPWS R7-97(2012);
//!   Wagner et al., J. Eng. Gas Turbines Power 122 (2000) 150): region 1
//!   (liquid, Gibbs energy), region 2 (vapour, Gibbs energy, ideal and
//!   residual parts), region 3 (near-critical, Helmholtz energy in density
//!   and temperature), region 4 (saturation line, both directions), region 5
//!   (high-temperature steam up to 2273.15 K, the 2007 revision of its
//!   coefficients) and the boundary B23 between regions 2 and 3. Range:
//!   273.15 K ≤ T ≤ 1073.15 K at p ≤ 100 MPa, and 1073.15 K < T ≤ 2273.15 K
//!   at p ≤ 50 MPa.
//! - **Region 3 at (p, T)** is solved for the density from the Helmholtz
//!   form itself (a scan of the isotherm on the side of the stable phase,
//!   then bisection with Newton steps), so the states are consistent with
//!   the basic equation to round-off; the backward equations v(p, T) of
//!   IAPWS SR5-05 are not needed.
//! - **(p, h), (p, s), (p, x) and (T, x) states** are solved from the
//!   forward equations (Newton steps on c_p and c_p/T kept inside a
//!   bracket); two-phase states carry the equilibrium quality, the
//!   homogeneous void fraction and the homogeneous mixture density. The
//!   backward equations T(p, h), T(p, s) of IF97 are not used: the forward
//!   solve is exact to the basic equations, which a 1D circuit marching in
//!   h (SF7b) needs for its energy closure.
//! - **Transport:** the viscosity of IAPWS R12-08 (Huber et al., J. Phys.
//!   Chem. Ref. Data 38 (2009) 101) with μ̄₂ = 1 (the release's
//!   recommendation for industrial use); the thermal conductivity of IAPWS
//!   R15-11 (Huber et al., J. Phys. Chem. Ref. Data 41 (2012) 033102) with
//!   its critical enhancement evaluated from the IF97 derivatives (the
//!   industrial route of the release; ζ at the reference temperature 1.5 T_c
//!   from a table of IAPWS-95, [`ZETA_R_INV`]); the surface tension of IAPWS
//!   R1-76(2014).
//! - **Two-phase mixtures:** homogeneous density and void fraction, and the
//!   mixture viscosity of McAdams, Cicchitti or Dukler ([`MixtureViscosity`]).
//!
//! Units are SI throughout: Pa, K, kg/m³, J/kg, J/(kg K), Pa s, W/(m K),
//! N/m. The bodies use the in-house `exp`, `ln` and `sqrt` of
//! `foehncraft_core::math` and fixed iteration rules, so a state evaluates
//! bit for bit the same on the CPU and the GPU (D31, D93 (6)). Every
//! function is total: a state outside the range of the formulation is
//! `None` (or NaN for the plain property functions), never a panic.

use foehncraft_core::math::{atan, exp, ln, sqrt};

/// Specific gas constant of IF97 (J/(kg K)).
pub const R: f64 = 461.526;
/// Critical temperature (K).
pub const T_CRIT: f64 = 647.096;
/// Critical pressure (Pa).
pub const P_CRIT: f64 = 22.064e6;
/// Critical density (kg/m³).
pub const RHO_CRIT: f64 = 322.0;
/// Lowest temperature of IF97 (K).
pub const T_MIN: f64 = 273.15;
/// Highest temperature of regions 1 to 3 (K).
pub const T_MAX_2: f64 = 1073.15;
/// Highest temperature of region 5 (K).
pub const T_MAX_5: f64 = 2273.15;
/// Highest pressure of regions 1 to 3 (Pa).
pub const P_MAX: f64 = 100e6;
/// Highest pressure of region 5 (Pa).
pub const P_MAX_5: f64 = 50e6;
/// The temperature of the boundary of regions 1 and 3 (K).
pub const T_13: f64 = 623.15;

/// x^n for an integer n (repeated squaring; a negative n gives 1/x^|n|).
#[inline(always)]
fn powi(x: f64, n: i32) -> f64 {
    let mut m = if n < 0 { -n } else { n } as u32;
    let mut b = x;
    let mut r = 1.0;
    while m > 0 {
        if m & 1 == 1 {
            r *= b;
        }
        b *= b;
        m >>= 1;
    }
    if n < 0 { 1.0 / r } else { r }
}

/// A Gibbs energy γ(π, τ) and its derivatives.
#[derive(Clone, Copy, Debug, Default)]
struct Gibbs {
    g: f64,
    gp: f64,
    gpp: f64,
    gt: f64,
    gtt: f64,
    gpt: f64,
}

/// A Helmholtz energy φ(δ, τ) and its derivatives.
#[derive(Clone, Copy, Debug, Default)]
struct Helmholtz {
    f: f64,
    fd: f64,
    fdd: f64,
    ft: f64,
    ftt: f64,
    fdt: f64,
}

// ---------------------------------------------------------------- region 1

const I1: [i32; 34] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29,
    30, 31, 32,
];
const J1: [i32; 34] = [
    -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11,
    -6, -29, -31, -38, -39, -40, -41,
];
const N1: [f64; 34] = [
    0.14632971213167,
    -0.84548187169114,
    -0.37563603672040e1,
    0.33855169168385e1,
    -0.95791963387872,
    0.15772038513228,
    -0.16616417199501e-1,
    0.81214629983568e-3,
    0.28319080123804e-3,
    -0.60706301565874e-3,
    -0.18990068218419e-1,
    -0.32529748770505e-1,
    -0.21841717175414e-1,
    -0.52838357969930e-4,
    -0.47184321073267e-3,
    -0.30001780793026e-3,
    0.47661393906987e-4,
    -0.44141845330846e-5,
    -0.72694996297594e-15,
    -0.31679644845054e-4,
    -0.28270797985312e-5,
    -0.85205128120103e-9,
    -0.22425281908000e-5,
    -0.65171222895601e-6,
    -0.14341729937924e-12,
    -0.40516996860117e-6,
    -0.12734301741641e-8,
    -0.17424871230634e-9,
    -0.68762131295531e-18,
    0.14478307828521e-19,
    0.26335781662795e-22,
    -0.11947622640071e-22,
    0.18228094581404e-23,
    -0.93537087292458e-25,
];

fn gibbs1(pi: f64, tau: f64) -> Gibbs {
    let a = 7.1 - pi;
    let b = tau - 1.222;
    let mut r = Gibbs::default();
    for k in 0..34 {
        let (i, j, n) = (I1[k], J1[k], N1[k]);
        let (fi, fj) = (i as f64, j as f64);
        let ai = powi(a, i);
        let ai1 = powi(a, i - 1);
        let bj = powi(b, j);
        let bj1 = powi(b, j - 1);
        r.g += n * ai * bj;
        r.gp -= n * fi * ai1 * bj;
        r.gpp += n * fi * (fi - 1.0) * powi(a, i - 2) * bj;
        r.gt += n * ai * fj * bj1;
        r.gtt += n * ai * fj * (fj - 1.0) * powi(b, j - 2);
        r.gpt -= n * fi * ai1 * fj * bj1;
    }
    r
}

// ---------------------------------------------------------------- region 2

const J0_2: [i32; 9] = [0, 1, -5, -4, -3, -2, -1, 2, 3];
const N0_2: [f64; 9] = [
    -0.96927686500217e1,
    0.10086655968018e2,
    -0.56087911283020e-2,
    0.71452738081455e-1,
    -0.40710498223928,
    0.14240819171444e1,
    -0.43839511319450e1,
    -0.28408632460772,
    0.21268463753307e-1,
];
const I2: [i32; 43] = [
    1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10,
    16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24,
];
const J2: [i32; 43] = [
    0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4,
    10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58,
];
const N2: [f64; 43] = [
    -0.17731742473213e-2,
    -0.17834862292358e-1,
    -0.45996013696365e-1,
    -0.57581259083432e-1,
    -0.50325278727930e-1,
    -0.33032641670203e-4,
    -0.18948987516315e-3,
    -0.39392777243355e-2,
    -0.43797295650573e-1,
    -0.26674547914087e-4,
    0.20481737692309e-7,
    0.43870667284435e-6,
    -0.32277677238570e-4,
    -0.15033924542148e-2,
    -0.40668253562649e-1,
    -0.78847309559367e-9,
    0.12790717852285e-7,
    0.48225372718507e-6,
    0.22922076337661e-5,
    -0.16714766451061e-10,
    -0.21171472321355e-2,
    -0.23895741934104e2,
    -0.59059564324270e-17,
    -0.12621808899101e-5,
    -0.38946842435739e-1,
    0.11256211360459e-10,
    -0.82311340897998e1,
    0.19809712802088e-7,
    0.10406965210174e-18,
    -0.10234747095929e-12,
    -0.10018179379511e-8,
    -0.80882908646985e-10,
    0.10693031879409,
    -0.33662250574171,
    0.89185845355421e-24,
    0.30629316876232e-12,
    -0.42002467698208e-5,
    -0.59056029685639e-25,
    0.37826947613457e-5,
    -0.12768608934681e-14,
    0.73087610595061e-28,
    0.55414715350778e-16,
    -0.94369707241210e-6,
];

fn gibbs2(pi: f64, tau: f64) -> Gibbs {
    let mut r = Gibbs {
        g: ln(pi),
        gp: 1.0 / pi,
        gpp: -1.0 / (pi * pi),
        ..Gibbs::default()
    };
    for k in 0..9 {
        let (j, n) = (J0_2[k], N0_2[k]);
        let fj = j as f64;
        r.g += n * powi(tau, j);
        r.gt += n * fj * powi(tau, j - 1);
        r.gtt += n * fj * (fj - 1.0) * powi(tau, j - 2);
    }
    let b = tau - 0.5;
    for k in 0..43 {
        let (i, j, n) = (I2[k], J2[k], N2[k]);
        let (fi, fj) = (i as f64, j as f64);
        let pi_i = powi(pi, i);
        let pi_i1 = powi(pi, i - 1);
        let bj = powi(b, j);
        let bj1 = powi(b, j - 1);
        r.g += n * pi_i * bj;
        r.gp += n * fi * pi_i1 * bj;
        r.gpp += n * fi * (fi - 1.0) * powi(pi, i - 2) * bj;
        r.gt += n * pi_i * fj * bj1;
        r.gtt += n * pi_i * fj * (fj - 1.0) * powi(b, j - 2);
        r.gpt += n * fi * pi_i1 * fj * bj1;
    }
    r
}

// ---------------------------------------------------------------- region 3

/// I of region 3; entry 0 is the logarithmic term n₁ ln δ.
const I3: [i32; 40] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6,
    6, 7, 8, 9, 9, 10, 10, 11,
];
const J3: [i32; 40] = [
    0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1,
    3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26,
];
const N3: [f64; 40] = [
    0.10658070028513e1,
    -0.15732845290239e2,
    0.20944396974307e2,
    -0.76867707878716e1,
    0.26185947787954e1,
    -0.28080781148620e1,
    0.12053369696517e1,
    -0.84566812812502e-2,
    -0.12654315477714e1,
    -0.11524407806681e1,
    0.88521043984318,
    -0.64207765181607,
    0.38493460186671,
    -0.85214708824206,
    0.48972281541877e1,
    -0.30502617256965e1,
    0.39420536879154e-1,
    0.12558408424308,
    -0.27999329698710,
    0.13899799569460e1,
    -0.20189915023570e1,
    -0.82147637173963e-2,
    -0.47596035734923,
    0.43984074473500e-1,
    -0.44476435428739,
    0.90572070719733,
    0.70522450087967,
    0.10770512626332,
    -0.32913623258954,
    -0.50871062041158,
    -0.22175400873096e-1,
    0.94260751665092e-1,
    0.16436278447961,
    -0.13503372241348e-1,
    -0.14834345352472e-1,
    0.57922953628084e-3,
    0.32308904703711e-2,
    0.80964802996215e-4,
    -0.16557679795037e-3,
    -0.44923899061815e-4,
];

fn helmholtz3(delta: f64, tau: f64) -> Helmholtz {
    let n1 = N3[0];
    let mut r = Helmholtz {
        f: n1 * ln(delta),
        fd: n1 / delta,
        fdd: -n1 / (delta * delta),
        ..Helmholtz::default()
    };
    for k in 1..40 {
        let (i, j, n) = (I3[k], J3[k], N3[k]);
        let (fi, fj) = (i as f64, j as f64);
        let di = powi(delta, i);
        let di1 = powi(delta, i - 1);
        let tj = powi(tau, j);
        let tj1 = powi(tau, j - 1);
        r.f += n * di * tj;
        r.fd += n * fi * di1 * tj;
        r.fdd += n * fi * (fi - 1.0) * powi(delta, i - 2) * tj;
        r.ft += n * di * fj * tj1;
        r.ftt += n * di * fj * (fj - 1.0) * powi(tau, j - 2);
        r.fdt += n * fi * di1 * fj * tj1;
    }
    r
}

/// The pressure (Pa) and (∂p/∂ρ)_T of region 3 at (ρ, T).
fn p3(rho: f64, t: f64) -> (f64, f64) {
    let d = rho / RHO_CRIT;
    let f = helmholtz3(d, T_CRIT / t);
    (
        rho * R * t * d * f.fd,
        R * t * (2.0 * d * f.fd + d * d * f.fdd),
    )
}

// ---------------------------------------------------------------- region 5

const J0_5: [i32; 6] = [0, 1, -3, -2, -1, 2];
const N0_5: [f64; 6] = [
    -0.13179983674201e2,
    0.68540841634434e1,
    -0.24805148933466e-1,
    0.36901534980333,
    -0.31161318213925e1,
    -0.32961626538917,
];
const I5: [i32; 6] = [1, 1, 1, 2, 2, 3];
const J5: [i32; 6] = [1, 2, 3, 3, 9, 7];
const N5: [f64; 6] = [
    0.15736404855259e-2,
    0.90153761673944e-3,
    -0.50270077677648e-2,
    0.22440037409485e-5,
    -0.41163275453471e-5,
    0.37919454822955e-7,
];

fn gibbs5(pi: f64, tau: f64) -> Gibbs {
    let mut r = Gibbs {
        g: ln(pi),
        gp: 1.0 / pi,
        gpp: -1.0 / (pi * pi),
        ..Gibbs::default()
    };
    for k in 0..6 {
        let (j, n) = (J0_5[k], N0_5[k]);
        let fj = j as f64;
        r.g += n * powi(tau, j);
        r.gt += n * fj * powi(tau, j - 1);
        r.gtt += n * fj * (fj - 1.0) * powi(tau, j - 2);
    }
    for k in 0..6 {
        let (i, j, n) = (I5[k], J5[k], N5[k]);
        let (fi, fj) = (i as f64, j as f64);
        let pi_i = powi(pi, i);
        let pi_i1 = powi(pi, i - 1);
        let tj = powi(tau, j);
        let tj1 = powi(tau, j - 1);
        r.g += n * pi_i * tj;
        r.gp += n * fi * pi_i1 * tj;
        r.gpp += n * fi * (fi - 1.0) * powi(pi, i - 2) * tj;
        r.gt += n * pi_i * fj * tj1;
        r.gtt += n * pi_i * fj * (fj - 1.0) * powi(tau, j - 2);
        r.gpt += n * fi * pi_i1 * fj * tj1;
    }
    r
}

// ------------------------------------------------- region 4 and boundary B23

const N4: [f64; 10] = [
    0.11670521452767e4,
    -0.72421316703206e6,
    -0.17073846940092e2,
    0.12020824702470e5,
    -0.32325550322333e7,
    0.14915108613530e2,
    -0.48232657361591e4,
    0.40511340542057e6,
    -0.23855557567849,
    0.65017534844798e3,
];

/// The saturation pressure p_s(T) (Pa) of region 4, for 273.15 K ≤ T ≤
/// T_c; NaN outside.
pub fn p_sat(t: f64) -> f64 {
    if !(T_MIN..=T_CRIT).contains(&t) {
        return f64::NAN;
    }
    let n = &N4;
    let th = t + n[8] / (t - n[9]);
    let a = th * th + n[0] * th + n[1];
    let b = n[2] * th * th + n[3] * th + n[4];
    let c = n[5] * th * th + n[6] * th + n[7];
    let x = 2.0 * c / (-b + sqrt(b * b - 4.0 * a * c));
    let x2 = x * x;
    x2 * x2 * 1e6
}

/// The saturation temperature T_s(p) (K) of region 4, for p_s(273.15 K) ≤
/// p ≤ p_c; NaN outside.
pub fn t_sat(p: f64) -> f64 {
    if !(P_SAT_MIN..=P_CRIT).contains(&p) {
        return f64::NAN;
    }
    let n = &N4;
    let beta = sqrt(sqrt(p * 1e-6));
    let e = beta * beta + n[2] * beta + n[5];
    let f = n[0] * beta * beta + n[3] * beta + n[6];
    let g = n[1] * beta * beta + n[4] * beta + n[7];
    let d = 2.0 * g / (-f - sqrt(f * f - 4.0 * e * g));
    let s = n[9] + d;
    0.5 * (s - sqrt(s * s - 4.0 * (n[8] + n[9] * d)))
}

/// The saturation pressure at 273.15 K (Pa), the lowest of region 4
/// (rounded down, so that `t_sat(p_sat(273.15))` is defined).
pub const P_SAT_MIN: f64 = 611.212_677;

const NB23: [f64; 5] = [
    0.34805185628969e3,
    -0.11671859879975e1,
    0.10192970039326e-2,
    0.57254459862746e3,
    0.13918839778870e2,
];

/// The pressure of the boundary B23 between regions 2 and 3 (Pa) at T (K).
pub fn p_b23(t: f64) -> f64 {
    (NB23[0] + NB23[1] * t + NB23[2] * t * t) * 1e6
}

/// The temperature of the boundary B23 (K) at p (Pa).
pub fn t_b23(p: f64) -> f64 {
    NB23[3] + sqrt((p * 1e-6 - NB23[4]) / NB23[2])
}

// ------------------------------------------------------------------ states

/// The phase of a state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum Phase {
    /// Liquid (below the saturation temperature, or below T_c above p_c).
    Liquid = 1,
    /// Vapour (above the saturation temperature, below p_c).
    Vapour = 2,
    /// A two-phase mixture on the saturation line.
    TwoPhase = 3,
    /// Supercritical fluid (p ≥ p_c and T ≥ T_c).
    Supercritical = 4,
}

/// A state of water or steam (SI units).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct State {
    /// The IF97 region (1, 2, 3, 5; 4 for a two-phase state).
    pub region: u8,
    /// The phase.
    pub phase: Phase,
    /// Pressure (Pa).
    pub p: f64,
    /// Temperature (K).
    pub t: f64,
    /// Density (kg/m³); of the homogeneous mixture in a two-phase state.
    pub rho: f64,
    /// Specific enthalpy (J/kg).
    pub h: f64,
    /// Specific internal energy (J/kg).
    pub u: f64,
    /// Specific entropy (J/(kg K)).
    pub s: f64,
    /// Isobaric specific heat (J/(kg K)); NaN in a two-phase state.
    pub cp: f64,
    /// Isochoric specific heat (J/(kg K)); NaN in a two-phase state.
    pub cv: f64,
    /// Speed of sound (m/s); NaN in a two-phase state.
    pub w: f64,
    /// (∂ρ/∂p)_T (s²/m²); NaN in a two-phase state.
    pub drho_dp: f64,
    /// The equilibrium quality x = (h − h′)/(h″ − h′) below p_c (below 0
    /// for a subcooled liquid, above 1 for a superheated vapour); NaN at
    /// and above p_c.
    pub x: f64,
    /// The homogeneous void fraction (0 for a liquid, 1 for a vapour, NaN
    /// above p_c).
    pub alpha: f64,
}

impl State {
    /// The specific volume (m³/kg).
    pub fn v(&self) -> f64 {
        1.0 / self.rho
    }
}

const NAN_STATE: State = State {
    region: 0,
    phase: Phase::Liquid,
    p: f64::NAN,
    t: f64::NAN,
    rho: f64::NAN,
    h: f64::NAN,
    u: f64::NAN,
    s: f64::NAN,
    cp: f64::NAN,
    cv: f64::NAN,
    w: f64::NAN,
    drho_dp: f64::NAN,
    x: f64::NAN,
    alpha: f64::NAN,
};

fn from_gibbs(region: u8, p: f64, t: f64, ps: f64, ts: f64, g: Gibbs) -> State {
    let pi = p / ps;
    let tau = ts / t;
    let v = pi * g.gp * R * t / p;
    let a = g.gp - tau * g.gpt;
    let cp = -R * tau * tau * g.gtt;
    let w2 = R * t * g.gp * g.gp / (a * a / (tau * tau * g.gtt) - g.gpp);
    // v = (R T / p*) γ_π, so (∂v/∂p)_T = R T γ_ππ / p*²
    let dv_dp = R * t * g.gpp / (ps * ps);
    State {
        region,
        phase: Phase::Liquid,
        p,
        t,
        rho: 1.0 / v,
        h: R * t * tau * g.gt,
        u: R * t * (tau * g.gt - pi * g.gp),
        s: R * (tau * g.gt - g.g),
        cp,
        cv: cp + R * a * a / g.gpp,
        w: sqrt(w2),
        drho_dp: -dv_dp / (v * v),
        x: f64::NAN,
        alpha: f64::NAN,
    }
}

/// Region 1 at (p, T), without range checks.
pub fn region1(p: f64, t: f64) -> State {
    from_gibbs(1, p, t, 16.53e6, 1386.0, gibbs1(p / 16.53e6, 1386.0 / t))
}

/// Region 2 at (p, T), without range checks.
pub fn region2(p: f64, t: f64) -> State {
    from_gibbs(2, p, t, 1e6, 540.0, gibbs2(p / 1e6, 540.0 / t))
}

/// Region 5 at (p, T), without range checks.
pub fn region5(p: f64, t: f64) -> State {
    from_gibbs(5, p, t, 1e6, 1000.0, gibbs5(p / 1e6, 1000.0 / t))
}

/// Region 3 at (ρ, T), without range checks.
pub fn region3(rho: f64, t: f64) -> State {
    let d = rho / RHO_CRIT;
    let tau = T_CRIT / t;
    let f = helmholtz3(d, tau);
    let dpd = 2.0 * d * f.fd + d * d * f.fdd;
    let a = d * f.fd - d * tau * f.fdt;
    let cv = -R * tau * tau * f.ftt;
    State {
        region: 3,
        phase: Phase::Liquid,
        p: rho * R * t * d * f.fd,
        t,
        rho,
        h: R * t * (tau * f.ft + d * f.fd),
        u: R * t * tau * f.ft,
        s: R * (tau * f.ft - f.f),
        cp: cv + R * a * a / dpd,
        cv,
        w: sqrt(R * t * (dpd - a * a / (tau * tau * f.ftt))),
        drho_dp: 1.0 / (R * t * dpd),
        x: f64::NAN,
        alpha: f64::NAN,
    }
}

/// The density of region 3 at (p, T) on the liquid side (`liquid`) or the
/// vapour side of the isotherm; NaN when the isotherm has no root there.
///
/// The isotherm is scanned in steps of 2 kg/m³ from a high density down
/// (liquid) or a low density up (vapour), stopping at the first sign change
/// of p₃(ρ) − p or where (∂p/∂ρ)_T turns negative (the spinodal); the root
/// is then bracketed and found by bisection with Newton steps. Above T_c the
/// isotherm is scanned from below.
pub fn rho3(p: f64, t: f64, liquid: bool) -> f64 {
    const STEP: f64 = 2.0;
    let liquid = liquid && t < T_CRIT;
    // the far end of the scan: well above any density of region 3 for the
    // liquid side, below the ideal-gas density for the vapour side
    let (mut a, dir) = if liquid {
        (820.0, -1.0)
    } else {
        (0.5 * p / (R * t), 1.0)
    };
    let (pa, da) = p3(a, t);
    let mut fa = pa - p;
    if da <= 0.0 || fa * dir > 0.0 {
        // the start is not on the stable branch of the side asked for
        return f64::NAN;
    }
    let mut b = a;
    let mut found = false;
    for _ in 0..512 {
        b = a + dir * STEP;
        if b <= 0.0 {
            return f64::NAN;
        }
        let (pb, db) = p3(b, t);
        let fb = pb - p;
        if fb * fa <= 0.0 && db > 0.0 {
            found = true;
            break;
        }
        if db <= 0.0 {
            // passed the spinodal: find it, the root lies before it
            let (mut lo, mut hi) = (a, b);
            for _ in 0..80 {
                let m = 0.5 * (lo + hi);
                if p3(m, t).1 > 0.0 {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            let fs = p3(lo, t).0 - p;
            if fs * fa > 0.0 {
                return f64::NAN;
            }
            b = lo;
            found = true;
            break;
        }
        a = b;
        fa = fb;
    }
    if !found {
        return f64::NAN;
    }
    // bisection with Newton steps inside the bracket [a, b]
    let (mut lo, mut hi) = if a < b { (a, b) } else { (b, a) };
    let mut rho = 0.5 * (lo + hi);
    let f_lo = p3(lo, t).0 - p;
    for _ in 0..200 {
        let (pr, dr) = p3(rho, t);
        let f = pr - p;
        if f == 0.0 {
            return rho;
        }
        if (f < 0.0) == (f_lo < 0.0) {
            lo = rho;
        } else {
            hi = rho;
        }
        let mut next = rho - f / dr;
        if !(next > lo && next < hi) || dr <= 0.0 {
            next = 0.5 * (lo + hi);
        }
        if (next - rho).abs() <= 1e-13 * rho || hi - lo <= 1e-13 * rho {
            return next;
        }
        rho = next;
    }
    rho
}

/// The IF97 region of (p, T) for a single-phase state (1, 2, 3 or 5); 0
/// outside the range. On the saturation line below 623.15 K the liquid
/// (region 1) is returned.
pub fn region_pt(p: f64, t: f64) -> u8 {
    if !(p > 0.0 && p <= P_MAX && (T_MIN..=T_MAX_5).contains(&t)) {
        return 0;
    }
    if t > T_MAX_2 {
        return if p <= P_MAX_5 { 5 } else { 0 };
    }
    if t <= T_13 {
        if p >= p_sat(t) { 1 } else { 2 }
    } else if t <= 863.15 && p > p_b23(t) {
        3
    } else {
        2
    }
}

/// The phase of a single-phase state at (p, T).
fn phase_pt(p: f64, t: f64) -> Phase {
    if p >= P_CRIT {
        if t >= T_CRIT {
            Phase::Supercritical
        } else {
            Phase::Liquid
        }
    } else if t >= T_CRIT || p < p_sat(t) {
        Phase::Vapour
    } else {
        Phase::Liquid
    }
}

/// The single-phase state at (p, T) without its quality (cheaper).
pub fn props_pt(p: f64, t: f64) -> Option<State> {
    let r = region_pt(p, t);
    let mut s = match r {
        1 => region1(p, t),
        2 => region2(p, t),
        5 => region5(p, t),
        3 => {
            let liquid = t < T_CRIT && p >= p_sat(t);
            let rho = rho3(p, t, liquid);
            if !rho.is_finite() {
                return None;
            }
            let mut s = region3(rho, t);
            // keep the pressure asked for (equal to round-off)
            s.p = p;
            s
        }
        _ => return None,
    };
    s.phase = phase_pt(p, t);
    s.alpha = match s.phase {
        Phase::Liquid if p < P_CRIT => 0.0,
        Phase::Vapour => 1.0,
        _ => f64::NAN,
    };
    Some(s)
}

/// The single-phase state at (p, T) with its equilibrium quality.
pub fn state_pt(p: f64, t: f64) -> Option<State> {
    let mut s = props_pt(p, t)?;
    if p < P_CRIT
        && let Some(sat) = saturation_p(p)
    {
        s.x = (s.h - sat.liquid.h) / (sat.vapour.h - sat.liquid.h);
    }
    Some(s)
}

/// Saturated liquid and vapour.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Saturation {
    /// Saturation temperature (K).
    pub t: f64,
    /// Saturation pressure (Pa).
    pub p: f64,
    /// The saturated liquid (x = 0).
    pub liquid: State,
    /// The saturated vapour (x = 1).
    pub vapour: State,
}

impl Saturation {
    /// The enthalpy of evaporation h″ − h′ (J/kg).
    pub fn h_fg(&self) -> f64 {
        self.vapour.h - self.liquid.h
    }
}

fn saturation(p: f64, t: f64) -> Option<Saturation> {
    let (mut l, mut v) = if t <= T_13 {
        (region1(p, t), region2(p, t))
    } else {
        let rl = rho3(p, t, true);
        let rv = rho3(p, t, false);
        if !(rl.is_finite() && rv.is_finite()) {
            return None;
        }
        (region3(rl, t), region3(rv, t))
    };
    l.p = p;
    v.p = p;
    l.phase = Phase::Liquid;
    v.phase = Phase::Vapour;
    l.x = 0.0;
    v.x = 1.0;
    l.alpha = 0.0;
    v.alpha = 1.0;
    Some(Saturation {
        t,
        p,
        liquid: l,
        vapour: v,
    })
}

/// The saturation state at the pressure p (Pa), p_s(273.15 K) ≤ p < p_c.
pub fn saturation_p(p: f64) -> Option<Saturation> {
    if !(P_SAT_MIN..P_CRIT).contains(&p) {
        return None;
    }
    saturation(p, t_sat(p))
}

/// The saturation state at the temperature t (K), 273.15 K ≤ t < T_c.
pub fn saturation_t(t: f64) -> Option<Saturation> {
    if !(T_MIN..T_CRIT).contains(&t) {
        return None;
    }
    saturation(p_sat(t), t)
}

/// The two-phase state of quality x (0 ≤ x ≤ 1) on the saturation line
/// `sat`: homogeneous mixture density and void fraction, linear h, u, s.
pub fn mixture(sat: &Saturation, x: f64) -> State {
    let (l, g) = (&sat.liquid, &sat.vapour);
    let v = (1.0 - x) / l.rho + x / g.rho;
    let mix = |a: f64, b: f64| a + x * (b - a);
    State {
        region: 4,
        phase: Phase::TwoPhase,
        p: sat.p,
        t: sat.t,
        rho: 1.0 / v,
        h: mix(l.h, g.h),
        u: mix(l.u, g.u),
        s: mix(l.s, g.s),
        x,
        alpha: x / g.rho / v,
        ..NAN_STATE
    }
}

/// The state at (p, x) for 0 ≤ x ≤ 1 below p_c.
pub fn state_px(p: f64, x: f64) -> Option<State> {
    if !(0.0..=1.0).contains(&x) {
        return None;
    }
    let sat = saturation_p(p)?;
    Some(endpoint(&sat, x))
}

/// The state at (T, x) for 0 ≤ x ≤ 1 below T_c.
pub fn state_tx(t: f64, x: f64) -> Option<State> {
    if !(0.0..=1.0).contains(&x) {
        return None;
    }
    let sat = saturation_t(t)?;
    Some(endpoint(&sat, x))
}

fn endpoint(sat: &Saturation, x: f64) -> State {
    if x == 0.0 {
        sat.liquid
    } else if x == 1.0 {
        sat.vapour
    } else {
        mixture(sat, x)
    }
}

/// Solves `f(T) = target` for T in [lo, hi], where `f` returns the value
/// and its derivative and increases with T; Newton steps kept inside the
/// bracket, bisection otherwise. `None` when the target is outside
/// [f(lo), f(hi)].
fn solve_t(
    target: f64,
    mut lo: f64,
    mut hi: f64,
    f: impl Fn(f64) -> Option<(f64, f64)>,
) -> Option<f64> {
    let (f_lo, _) = f(lo)?;
    let (f_hi, _) = f(hi)?;
    let scale = 1e-12 * (target.abs() + 1.0);
    if target < f_lo - scale || target > f_hi + scale {
        return None;
    }
    if (target - f_lo).abs() <= scale {
        return Some(lo);
    }
    if (target - f_hi).abs() <= scale {
        return Some(hi);
    }
    // start by linear interpolation
    let mut t = lo + (target - f_lo) / (f_hi - f_lo) * (hi - lo);
    for _ in 0..200 {
        let (v, d) = f(t)?;
        let r = v - target;
        if r == 0.0 {
            return Some(t);
        }
        if r < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
        let mut next = t - r / d;
        if !(next > lo && next < hi) || !(d > 0.0) {
            next = 0.5 * (lo + hi);
        }
        if (next - t).abs() <= 1e-12 * t || hi - lo <= 1e-12 * t {
            return Some(next);
        }
        t = next;
    }
    Some(t)
}

/// The upper temperature limit at p (K).
fn t_max(p: f64) -> f64 {
    if p <= P_MAX_5 { T_MAX_5 } else { T_MAX_2 }
}

/// What a (p, ·) solve matches: the enthalpy or the entropy.
#[derive(Clone, Copy)]
enum By {
    H,
    S,
}

fn by_value(s: &State, by: By) -> (f64, f64) {
    match by {
        By::H => (s.h, s.cp),
        By::S => (s.s, s.cp / s.t),
    }
}

fn state_p_by(p: f64, target: f64, by: By) -> Option<State> {
    if !(p > 0.0 && p <= P_MAX) {
        return None;
    }
    let f = |t: f64| props_pt(p, t).map(|s| by_value(&s, by));
    let solve = |lo: f64, hi: f64| -> Option<State> {
        let t = solve_t(target, lo, hi, f)?;
        state_pt(p, t)
    };
    if p < P_CRIT && p >= P_SAT_MIN {
        let sat = saturation_p(p)?;
        let (fl, _) = by_value(&sat.liquid, by);
        let (fg, _) = by_value(&sat.vapour, by);
        if target >= fl && target <= fg {
            let x = (target - fl) / (fg - fl);
            return Some(endpoint(&sat, x));
        }
        if target < fl {
            // the liquid side up to just below T_s
            let mut s = solve(T_MIN, sat.t * (1.0 - 1e-14))?;
            s.x = (s.h - sat.liquid.h) / sat.h_fg();
            return Some(s);
        }
        let mut s = solve(sat.t * (1.0 + 1e-14), t_max(p))?;
        s.x = (s.h - sat.liquid.h) / sat.h_fg();
        return Some(s);
    }
    solve(T_MIN, t_max(p))
}

/// The state at (p, h): single phase or two-phase.
pub fn state_ph(p: f64, h: f64) -> Option<State> {
    state_p_by(p, h, By::H)
}

/// The state at (p, s): single phase or two-phase.
pub fn state_ps(p: f64, s: f64) -> Option<State> {
    state_p_by(p, s, By::S)
}

// --------------------------------------------------------------- transport

const H0: [f64; 4] = [1.67752, 2.20462, 0.6366564, -0.241605];
/// H_ij of IAPWS R12-08, i = 0..5 (rows), j = 0..6 (columns).
const H1: [[f64; 7]; 6] = [
    [
        5.20094e-1,
        2.22531e-1,
        -2.81378e-1,
        1.61913e-1,
        -3.25372e-2,
        0.0,
        0.0,
    ],
    [
        8.50895e-2,
        9.99115e-1,
        -9.06851e-1,
        2.57399e-1,
        0.0,
        0.0,
        0.0,
    ],
    [-1.08374, 1.88797, -7.72479e-1, 0.0, 0.0, 0.0, 0.0],
    [
        -2.89555e-1,
        1.26613,
        -4.89837e-1,
        0.0,
        6.98452e-2,
        0.0,
        -4.35673e-3,
    ],
    [0.0, 0.0, -2.57040e-1, 0.0, 0.0, 8.72102e-3, 0.0],
    [0.0, 1.20573e-1, 0.0, 0.0, 0.0, 0.0, -5.93264e-4],
];

/// The dilute-gas and residual factors μ̄₀, μ̄₁ (μ in units of 1e-6 Pa s).
fn mu01(rho: f64, t: f64) -> (f64, f64) {
    let tb = t / T_CRIT;
    let rb = rho / RHO_CRIT;
    let mut d = 0.0;
    let mut ti = 1.0;
    for h in H0 {
        d += h / ti;
        ti *= tb;
    }
    let mu0 = 100.0 * sqrt(tb) / d;
    let a = 1.0 / tb - 1.0;
    let b = rb - 1.0;
    let mut sum = 0.0;
    let mut ai = 1.0;
    for row in &H1 {
        let mut bj = 1.0;
        let mut inner = 0.0;
        for &h in row {
            inner += h * bj;
            bj *= b;
        }
        sum += ai * inner;
        ai *= a;
    }
    (mu0, exp(rb * sum))
}

/// The dynamic viscosity (Pa s) at (ρ, T) by IAPWS R12-08 with μ̄₂ = 1
/// (the release's form for industrial use; the critical enhancement is
/// below 2 % outside 645.91 K < T < 650.77 K and 245.8 < ρ < 405.3 kg/m³).
pub fn viscosity(rho: f64, t: f64) -> f64 {
    let (mu0, mu1) = mu01(rho, t);
    1e-6 * mu0 * mu1
}

const L0: [f64; 5] = [
    2.443221e-3,
    1.323095e-2,
    6.770357e-3,
    -3.454586e-3,
    4.096266e-4,
];
/// L_ij of IAPWS R15-11, i = 0..4 (rows), j = 0..5 (columns).
const L1: [[f64; 6]; 5] = [
    [
        1.60397357,
        -0.646013523,
        0.111443906,
        0.102997357,
        -0.0504123634,
        0.00609859258,
    ],
    [
        2.33771842,
        -2.78843778,
        1.53616167,
        -0.463045512,
        0.0832827019,
        -0.00719201245,
    ],
    [
        2.19650529,
        -4.54580785,
        3.55777244,
        -1.40944978,
        0.275418278,
        -0.0205938816,
    ],
    [
        -1.21051378,
        1.60812989,
        -0.621178141,
        0.0716373224,
        0.0,
        0.0,
    ],
    [
        -2.7203370,
        4.57586331,
        -3.18369245,
        1.1168348,
        -0.19268305,
        0.012913842,
    ],
];

/// 1/ζ(T̄_R, ρ̄) with ζ = (∂ρ̄/∂p̄)_T̄ at the reference temperature T_R =
/// 1.5 T_c of the critical enhancement of IAPWS R15-11, at ρ̄ = 0, 0.025,
/// …, 3.2 (ρ̄ = ρ/322 kg/m³, p̄ = p/22.064 MPa). Evaluated from IAPWS-95
/// (Wagner & Pruss, J. Phys. Chem. Ref. Data 31 (2002) 387) with CoolProp
/// 8.1 (MIT licence), as the release prescribes; IF97 does not reach the
/// liquid densities at 970.644 K, so the release's industrial form also
/// takes ζ_R from a fit (its Table 6). Interpolated linearly, clamped.
pub const ZETA_R_INV: [f64; 129] = [
    6.537623745,
    6.400861755,
    6.266520538,
    6.135459630,
    6.008397616,
    5.885929161,
    5.768536455,
    5.656600515,
    5.550412319,
    5.450183704,
    5.356057942,
    5.268119954,
    5.186406054,
    5.110913179,
    5.041607527,
    4.978432552,
    4.921316275,
    4.870177870,
    4.824933521,
    4.785501519,
    4.751806639,
    4.723783786,
    4.701380958,
    4.684561566,
    4.673306143,
    4.667613523,
    4.667501518,
    4.673007174,
    4.684186656,
    4.701114824,
    4.723884556,
    4.752605866,
    4.787404880,
    4.828422688,
    4.875814146,
    4.929746622,
    4.990398741,
    5.057959136,
    5.132625204,
    5.214601885,
    5.304100445,
    5.401337245,
    5.506532493,
    5.619908943,
    5.741690557,
    5.872101126,
    6.011362905,
    6.159695325,
    6.317313896,
    6.484429446,
    6.661247830,
    6.847970278,
    7.044794493,
    7.251916539,
    7.469533507,
    7.697846820,
    7.937065947,
    8.187412255,
    8.449122663,
    8.722452820,
    9.007679603,
    9.305102846,
    9.615046373,
    9.937858511,
    1.027391237e1,
    1.062360616e1,
    1.098736386e1,
    1.136563623e1,
    1.175890237e1,
    1.216767140e1,
    1.259248428e1,
    1.303391522e1,
    1.349257251e1,
    1.396909848e1,
    1.446416850e1,
    1.497848897e1,
    1.551279440e1,
    1.606784377e1,
    1.664441632e1,
    1.724330691e1,
    1.786532130e1,
    1.851127118e1,
    1.918196930e1,
    1.987822465e1,
    2.060083781e1,
    2.135059642e1,
    2.212827103e1,
    2.293461112e1,
    2.377034167e1,
    2.463616008e1,
    2.553273363e1,
    2.646069753e1,
    2.742065352e1,
    2.841316916e1,
    2.943877767e1,
    3.049797847e1,
    3.159123829e1,
    3.271899286e1,
    3.388164913e1,
    3.507958797e1,
    3.631316727e1,
    3.758272543e1,
    3.888858502e1,
    4.023105678e1,
    4.161044355e1,
    4.302704438e1,
    4.448115854e1,
    4.597308937e1,
    4.750314801e1,
    4.907165686e1,
    5.067895270e1,
    5.232538954e1,
    5.401134098e1,
    5.573720226e1,
    5.750339180e1,
    5.931035236e1,
    6.115855170e1,
    6.304848285e1,
    6.498066393e1,
    6.695563756e1,
    6.897396996e1,
    7.103624961e1,
    7.314308569e1,
    7.529510621e1,
    7.749295593e1,
    7.973729405e1,
    8.202879184e1,
    8.436813008e1,
    8.675599652e1,
];

/// The thermal conductivity (W/(m K)) of a single-phase state by IAPWS
/// R15-11, with the critical enhancement λ̄₂ from the IF97 derivatives of
/// the state (the industrial route of the release); NaN in a two-phase
/// state.
pub fn conductivity(s: &State) -> f64 {
    if s.phase == Phase::TwoPhase || !s.cp.is_finite() {
        return f64::NAN;
    }
    let (rho, t) = (s.rho, s.t);
    let tb = t / T_CRIT;
    let rb = rho / RHO_CRIT;
    let mut d = 0.0;
    let mut ti = 1.0;
    for l in L0 {
        d += l / ti;
        ti *= tb;
    }
    let l0 = sqrt(tb) / d;
    let a = 1.0 / tb - 1.0;
    let b = rb - 1.0;
    let mut sum = 0.0;
    let mut ai = 1.0;
    for row in &L1 {
        let mut bj = 1.0;
        let mut inner = 0.0;
        for &l in row {
            inner += l * bj;
            bj *= b;
        }
        sum += ai * inner;
        ai *= a;
    }
    let l1 = exp(rb * sum);
    1e-3 * (l0 * l1 + lambda2(s, tb, rb))
}

/// λ̄₂ of IAPWS R15-11 (in units of 1e-3 W/(m K)).
fn lambda2(s: &State, tb: f64, rb: f64) -> f64 {
    const LAMBDA: f64 = 177.8514;
    const QD: f64 = 1.0 / 0.40e-9;
    const NU: f64 = 0.630;
    const GAMMA: f64 = 1.239;
    const XI0: f64 = 0.13e-9;
    const GAMMA0: f64 = 0.06;
    const TBR: f64 = 1.5;
    // the gas constant of the release
    const R_REL: f64 = 461.51805;
    let zeta = s.drho_dp * P_CRIT / RHO_CRIT;
    let zeta_r = 1.0 / zeta_r_inv(rb);
    let dchi = rb * (zeta - zeta_r * TBR / tb);
    if !(dchi > 0.0) {
        return 0.0;
    }
    let xi = XI0 * exp(NU / GAMMA * ln(dchi / GAMMA0));
    let y = QD * xi;
    let cpb = s.cp / R_REL;
    let kappa = s.cp / s.cv;
    let (mu0, mu1) = mu01(s.rho, s.t);
    let mub = mu0 * mu1;
    let z = if y < 1.2e-7 {
        0.0
    } else {
        let ik = 1.0 / kappa;
        2.0 / (core::f64::consts::PI * y)
            * ((1.0 - ik) * atan(y) + ik * y
                - (1.0 - exp(-1.0 / (1.0 / y + y * y / (3.0 * rb * rb)))))
    };
    LAMBDA * rb * cpb * tb / mub * z
}

fn zeta_r_inv(rb: f64) -> f64 {
    let n = ZETA_R_INV.len();
    let x = (rb / 0.025).clamp(0.0, (n - 1) as f64);
    let i = (x as usize).min(n - 2);
    let w = x - i as f64;
    ZETA_R_INV[i] + w * (ZETA_R_INV[i + 1] - ZETA_R_INV[i])
}

/// The surface tension (N/m) of water against its vapour at T (K) by
/// IAPWS R1-76(2014); NaN outside 273.15 K ≤ T ≤ T_c.
pub fn surface_tension(t: f64) -> f64 {
    if !(T_MIN..=T_CRIT).contains(&t) {
        return f64::NAN;
    }
    let tau = 1.0 - t / T_CRIT;
    if tau == 0.0 {
        return 0.0;
    }
    235.8e-3 * exp(1.256 * ln(tau)) * (1.0 - 0.625 * tau)
}

/// A mixture viscosity of a homogeneous two-phase flow.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum MixtureViscosity {
    /// 1/μ = x/μ″ + (1 − x)/μ′ (McAdams, Woods & Heroman 1942).
    McAdams = 0,
    /// μ = x μ″ + (1 − x) μ′ (Cicchitti et al. 1960).
    Cicchitti = 1,
    /// μ = ρ_m [x μ″/ρ″ + (1 − x) μ′/ρ′] (Dukler, Wicks & Cleveland 1964).
    Dukler = 2,
}

/// The mixture viscosity (Pa s) at the quality x on the saturation line.
pub fn mixture_viscosity(sat: &Saturation, x: f64, model: MixtureViscosity) -> f64 {
    let (l, g) = (&sat.liquid, &sat.vapour);
    let ml = viscosity(l.rho, l.t);
    let mg = viscosity(g.rho, g.t);
    match model {
        MixtureViscosity::McAdams => 1.0 / (x / mg + (1.0 - x) / ml),
        MixtureViscosity::Cicchitti => x * mg + (1.0 - x) * ml,
        MixtureViscosity::Dukler => {
            let v = (1.0 - x) / l.rho + x / g.rho;
            (x * mg / g.rho + (1.0 - x) * ml / l.rho) / v
        }
    }
}

// ------------------------------------------------------------ batch entry

/// The number of outputs per point of [`state_entry`].
pub const OUT_WORDS: usize = 8;

/// How the points of a batch are given.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u32)]
pub enum Input {
    /// a = p (Pa), b = T (K).
    PT = 0,
    /// a = p (Pa), b = h (J/kg).
    PH = 1,
}

/// The outputs of one point of a batch, in this order: ρ, h, T, c_p, μ, k,
/// x, region (as a number); NaN for a point outside the range.
pub fn state_outputs(input: Input, a: f64, b: f64) -> [f64; OUT_WORDS] {
    let s = match input {
        Input::PT => state_pt(a, b),
        Input::PH => state_ph(a, b),
    };
    match s {
        Some(s) => {
            let single = s.phase != Phase::TwoPhase;
            let mu = if single {
                viscosity(s.rho, s.t)
            } else {
                f64::NAN
            };
            [
                s.rho,
                s.h,
                s.t,
                s.cp,
                mu,
                conductivity(&s),
                s.x,
                s.region as f64,
            ]
        }
        None => [f64::NAN; OUT_WORDS],
    }
}

/// The argument block of the batch entry: the inputs a, b of each point and
/// the outputs, [`OUT_WORDS`] per point.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct StateArgs {
    /// How the inputs are given ([`Input`] as a number).
    pub input: u32,
    /// The first input of each point.
    pub a: crate::entry::Raw<f64>,
    /// The second input of each point.
    pub b: crate::entry::Raw<f64>,
    /// The outputs (written).
    pub out: crate::entry::Raw<f64>,
}

/// The batch body: the state of point `i`.
///
/// # Safety
/// `i` is below the length of `a` and `b`; `out` holds `OUT_WORDS` entries
/// per point and no other thread writes the entries of point `i`.
#[inline(always)]
pub unsafe fn state_entry(i: usize, args: &StateArgs) {
    let input = if args.input == 1 {
        Input::PH
    } else {
        Input::PT
    };
    // SAFETY: the caller guarantees that `i` lies inside `a` and `b`.
    let (a, b) = unsafe { (crate::ld(args.a.get(), i), crate::ld(args.b.get(), i)) };
    let o = state_outputs(input, a, b);
    let base = args.out.ptr as *mut f64;
    for (k, v) in o.iter().enumerate() {
        // SAFETY: the caller guarantees `OUT_WORDS` entries of point `i`
        // inside `out`, written by this index only.
        unsafe { *base.add(i * OUT_WORDS + k) = *v };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entry::Raw;

    /// The device path (argument block, entry body) run on the CPU gives
    /// bitwise the outputs of the plain call, as the GPU does.
    #[test]
    fn batch_entry_is_bitwise_the_state_outputs() {
        const N: usize = 36;
        let ps = [0.005e6, 0.1e6, 4e6, 18e6, 25e6, 60e6];
        let ts = [280.0, 450.0, 640.0, 660.0, 900.0, 1800.0];
        let mut a = [0.0f64; N];
        let mut b = [0.0f64; N];
        for i in 0..N {
            a[i] = ps[i / 6];
            b[i] = ts[i % 6];
        }
        let mut out = [0.0f64; N * OUT_WORDS];
        let args = StateArgs {
            input: Input::PT as u32,
            a: Raw::new(a.as_ptr() as u64, N),
            b: Raw::new(b.as_ptr() as u64, N),
            out: Raw::new(out.as_mut_ptr() as u64, N * OUT_WORDS),
        };
        for i in 0..N {
            // SAFETY: i < N, the length of a and b; out holds OUT_WORDS
            // entries per point.
            unsafe { state_entry(i, &args) };
        }
        let mut finite = 0;
        for i in 0..N {
            let want = state_outputs(Input::PT, a[i], b[i]);
            for k in 0..OUT_WORDS {
                assert_eq!(out[i * OUT_WORDS + k].to_bits(), want[k].to_bits());
            }
            finite += usize::from(want[0].is_finite());
        }
        // 60 MPa at 1800 K lies outside region 5
        assert_eq!(finite, N - 1);
    }
}
