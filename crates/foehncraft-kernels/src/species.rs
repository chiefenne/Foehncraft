//! Mixture properties of an ideal-gas mixture (steam furnace stage SF1):
//! the bodies that turn the composition and the temperature of one cell or
//! face into its molar mass, specific heat, viscosity and conductivity.
//!
//! Every function is cell-local, allocation-free and sums in the fixed
//! order of the species, so the CPU and a GPU thread give the same bits
//! (D31) and a rank computes its cells as one rank does (D93).
//!
//! - Specific heat and enthalpy: NASA 7-coefficient polynomials (McBride,
//!   Gordon & Reno 1993; the GRI-Mech 3.0 THERMO data), two ranges, held at
//!   the value of the nearest limit outside the data range.
//! - Pure-species viscosity: Chapman–Enskog with the Lennard-Jones
//!   collision integral Ω(2,2)* of Neufeld, Janzen & Aziz (1972), with
//!   Brokaw's (1969) correction 0.2 δ*²/T* for polar molecules, δ* =
//!   μ_D²/(2 ε σ³) (Poling, Prausnitz & O'Connell 2001, eqs. 9-3.9, 9-4.3,
//!   9-4.4).
//! - Pure-species conductivity: for molecules the mean of the Eucken
//!   relation λ M/(μ C_v) = 1 + 2.25 R/C_v and the modified Eucken relation
//!   1.32 + 1.77 R/C_v, which bracket measured values (Poling et al. 2001,
//!   section 10-3, eqs. 10-3.3, 10-3.4); λ = 15/4 R μ/M for atoms (eq.
//!   10-3.1).
//! - Mixture viscosity: Wilke (1950); mixture conductivity: Mathur, Tondon
//!   & Saxena (1967), λ = ½ (Σ x_i λ_i + 1/Σ (x_i/λ_i)).
//! - Binary diffusion coefficients (for Lewis-number checks and the
//!   mixture-averaged model of SF2h): Chapman–Enskog with Ω(1,1)* of
//!   Neufeld et al. (1972), Poling et al. 2001, eq. 11-3.2.

use foehncraft_core::Real;

/// The universal gas constant R (J/(mol K)), CODATA 2018.
pub const R: f64 = 8.314_462_618;

/// The Boltzmann constant in erg/K (cgs, for δ*).
const K_B_CGS: f64 = 1.380_649e-16;

/// One species of a mixture: its NASA-7 polynomials and its Lennard-Jones
/// transport data.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Species {
    /// Molar mass (kg/mol).
    pub molar_mass: f64,
    /// Lower, middle and upper temperature of the two ranges (K).
    pub t_range: [f64; 3],
    /// Coefficients a₁…a₇ of the lower range.
    pub low: [f64; 7],
    /// Coefficients a₁…a₇ of the upper range.
    pub high: [f64; 7],
    /// ε/k_B of the Lennard-Jones potential (K).
    pub well_depth: f64,
    /// σ of the Lennard-Jones potential (Å).
    pub diameter: f64,
    /// Dipole moment (Debye).
    pub dipole: f64,
    /// 0 for an atom, 1 for a linear and 2 for a nonlinear molecule.
    pub geometry: u32,
}

impl Species {
    /// The coefficients of the range of `t` and `t` held inside the data
    /// range.
    #[inline(always)]
    fn range(&self, t: f64) -> (&[f64; 7], f64) {
        let t = t.max(self.t_range[0]).min(self.t_range[2]);
        if t <= self.t_range[1] {
            (&self.low, t)
        } else {
            (&self.high, t)
        }
    }

    /// c_p/R at `t`.
    #[inline(always)]
    pub fn cp_over_r(&self, t: f64) -> f64 {
        let (a, t) = self.range(t);
        a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4])))
    }

    /// h/(R T) at `t` inside the data range (the enthalpy of formation
    /// included).
    #[inline(always)]
    pub fn h_over_rt(&self, t: f64) -> f64 {
        let (a, t) = self.range(t);
        a[0] + t * (a[1] / 2.0 + t * (a[2] / 3.0 + t * (a[3] / 4.0 + t * a[4] / 5.0))) + a[5] / t
    }

    /// Specific heat c_p (J/(kg K)) at `t`.
    #[inline(always)]
    pub fn cp(&self, t: f64) -> f64 {
        R * self.cp_over_r(t) / self.molar_mass
    }

    /// Dynamic viscosity of the pure gas (Pa s) at `t`.
    #[inline(always)]
    pub fn viscosity(&self, t: f64) -> f64 {
        let ts = t / self.well_depth;
        let mut om = omega22(ts);
        if self.dipole > 0.0 {
            let d = self.dipole * 1e-18;
            let s = self.diameter * 1e-8;
            let delta = d * d / (2.0 * self.well_depth * K_B_CGS * s * s * s);
            om += 0.2 * delta * delta / ts;
        }
        // 26.69 μP √(M T)/(σ² Ω) with M in g/mol, σ in Å
        2.6693e-6 * sqrt(1000.0 * self.molar_mass * t) / (self.diameter * self.diameter * om)
    }

    /// Thermal conductivity of the pure gas (W/(m K)) at `t` with its
    /// viscosity `mu`.
    #[inline(always)]
    pub fn conductivity(&self, t: f64, mu: f64) -> f64 {
        let cv = R * (self.cp_over_r(t) - 1.0); // molar C_v
        if self.geometry == 0 {
            3.75 * R * mu / self.molar_mass
        } else {
            // the mean of the Eucken (λM/(μC_v) = 1 + 2.25 R/C_v) and the
            // modified Eucken (1.32 + 1.77 R/C_v) values, between which
            // measured conductivities of polyatomic gases lie
            mu / self.molar_mass * (1.16 * cv + 2.01 * R)
        }
    }
}

/// √x on every target (no_std).
#[inline(always)]
fn sqrt(x: f64) -> f64 {
    x.sqrt_of()
}

/// exp on every target (no_std).
#[inline(always)]
fn exp(x: f64) -> f64 {
    x.exp_of()
}

/// x^y = exp(y ln x) for x > 0.
#[inline(always)]
fn powf(x: f64, y: f64) -> f64 {
    exp(y * x.ln_of())
}

/// The reduced collision integral Ω(2,2)* at T* = k_B T/ε (Neufeld, Janzen
/// & Aziz 1972), 0.3 ≤ T* ≤ 100 within 0.1 %.
#[inline(always)]
pub fn omega22(ts: f64) -> f64 {
    1.16145 * powf(ts, -0.14874) + 0.52487 * exp(-0.77320 * ts) + 2.16178 * exp(-2.43787 * ts)
}

/// The reduced collision integral Ω(1,1)* at T* (Neufeld et al. 1972).
#[inline(always)]
pub fn omega11(ts: f64) -> f64 {
    1.06036 * powf(ts, -0.15610)
        + 0.19300 * exp(-0.47635 * ts)
        + 1.03587 * exp(-1.52996 * ts)
        + 1.76474 * exp(-3.89411 * ts)
}

/// The binary diffusion coefficient D_ij (m²/s) of two species at `t` and
/// the pressure `p` (Pa), Chapman–Enskog with the combining rules σ_ij =
/// (σ_i + σ_j)/2, ε_ij = √(ε_i ε_j) (Poling et al. 2001, eq. 11-3.2).
#[inline(always)]
pub fn binary_diffusion(a: &Species, b: &Species, t: f64, p: f64) -> f64 {
    let m = 2.0 / (1.0 / (1000.0 * a.molar_mass) + 1.0 / (1000.0 * b.molar_mass));
    let s = 0.5 * (a.diameter + b.diameter);
    let e = sqrt(a.well_depth * b.well_depth);
    // 0.00266 T^1.5/(p[bar] √M σ² Ω_D) cm²/s
    let d = 0.00266 * t * sqrt(t) / (p * 1e-5 * sqrt(m) * s * s * omega11(t / e));
    d * 1e-4
}

/// The properties of a mixture at one point.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Mixture {
    /// Mean molar mass (kg/mol).
    pub molar_mass: f64,
    /// Specific heat c_p (J/(kg K)).
    pub cp: f64,
    /// Dynamic viscosity (Pa s).
    pub viscosity: f64,
    /// Thermal conductivity (W/(m K)).
    pub conductivity: f64,
}

/// The largest number of species of a mixture of SF1 (the 53 species of
/// GRI-Mech 3.0 with room).
pub const MAX_SPECIES: usize = 64;

/// The properties of the mixture of the species `sp` with the mass
/// fractions `y` at the temperature `t`; the work arrays hold the mole
/// fractions and the pure-species viscosities and conductivities. Mass
/// fractions below zero count as zero, and the fractions are taken as
/// they are (the caller closes them, Σ y = 1).
#[inline(always)]
pub fn mixture(sp: &[Species], y: &[f64], t: f64) -> Mixture {
    let n = sp.len();
    let mut x = [0.0; MAX_SPECIES];
    let mut mu = [0.0; MAX_SPECIES];
    let mut k = [0.0; MAX_SPECIES];
    // 1/M = Σ y_i/M_i
    let mut inv_m = 0.0;
    let mut cp = 0.0;
    for i in 0..n {
        let yi = y[i].max(0.0);
        inv_m += yi / sp[i].molar_mass;
        cp += yi * sp[i].cp(t);
    }
    let m = 1.0 / inv_m;
    for i in 0..n {
        x[i] = y[i].max(0.0) * m / sp[i].molar_mass;
        mu[i] = sp[i].viscosity(t);
        k[i] = sp[i].conductivity(t, mu[i]);
    }
    // Wilke (1950)
    let mut visc = 0.0;
    for i in 0..n {
        if x[i] == 0.0 {
            continue;
        }
        let mut den = 0.0;
        for j in 0..n {
            if x[j] == 0.0 {
                continue;
            }
            let r = sqrt(mu[i] / mu[j]) * sqrt(sqrt(sp[j].molar_mass / sp[i].molar_mass));
            let phi =
                (1.0 + r) * (1.0 + r) / sqrt(8.0 * (1.0 + sp[i].molar_mass / sp[j].molar_mass));
            den += x[j] * phi;
        }
        visc += x[i] * mu[i] / den;
    }
    // Mathur, Tondon & Saxena (1967)
    let (mut s1, mut s2) = (0.0, 0.0);
    for i in 0..n {
        if x[i] == 0.0 {
            continue;
        }
        s1 += x[i] * k[i];
        s2 += x[i] / k[i];
    }
    Mixture {
        molar_mass: m,
        cp,
        viscosity: visc,
        conductivity: 0.5 * (s1 + 1.0 / s2),
    }
}
