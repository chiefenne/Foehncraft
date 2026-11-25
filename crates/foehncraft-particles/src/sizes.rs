//! Droplet size distributions of the injectors (D242).
//!
//! Both distributions describe how the injected liquid **mass** is spread
//! over the droplet diameter, as sprays are usually specified. Every parcel
//! carries the same mass, so drawing the diameter of each parcel from the
//! mass distribution gives the right mass in every size range; the number
//! of droplets of a parcel follows from its mass and diameter.
//!
//! - **Rosin-Rammler** (Rosin and Rammler 1933): the mass fraction of
//!   droplets larger than d is Y(d) = exp(−(d/X)ⁿ), with the size parameter
//!   X (`mean`) and the spread n (`spread`).
//! - **Log-normal**: ln d of the mass is normal with the mass median
//!   diameter (`median`) and the geometric standard deviation σ_g
//!   (`sigma`, > 1).
//!
//! Both are truncated to [`min`, `max`]. A diameter is the inverse of the
//! truncated cumulative distribution at a uniform number from the
//! counter-based generator keyed by the parcel identity (D121), so it does
//! not depend on the thread count or the order of the parcels.
//!
//! The inverse of the normal distribution is the rational approximation of
//! Acklam (2003, relative error below 1.15e-9) refined by one Halley step on
//! erfc (Abramowitz and Stegun 26.2.2), which brings it to round-off.

use anyhow::{Result, bail};
use foehncraft_core::math;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// A droplet size distribution of the liquid mass, `size = { … }` of an
/// injector.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "distribution", rename_all = "kebab-case", deny_unknown_fields)]
pub enum SizeDistribution {
    /// Rosin-Rammler: mass fraction above d is exp(−(d/mean)^spread).
    RosinRammler {
        /// Size parameter X (m): 63.2 % of the mass is in droplets smaller
        /// than X.
        mean: f64,
        /// Spread n (−), typically 1.5 to 4.
        spread: f64,
        /// Smallest diameter (m).
        min: f64,
        /// Largest diameter (m).
        max: f64,
    },
    /// Log-normal distribution of the mass over ln d.
    LogNormal {
        /// Mass median diameter (m).
        median: f64,
        /// Geometric standard deviation σ_g (−), larger than 1.
        sigma: f64,
        /// Smallest diameter (m).
        min: f64,
        /// Largest diameter (m).
        max: f64,
    },
}

impl SizeDistribution {
    /// Checks the values.
    pub fn validate(&self) -> Result<()> {
        let pos = |x: f64| x.is_finite() && x > 0.0;
        let (min, max) = self.bounds();
        match *self {
            Self::RosinRammler { mean, spread, .. } => {
                if !pos(mean) {
                    bail!("mean must be positive and finite");
                }
                if !pos(spread) {
                    bail!("spread must be positive and finite");
                }
            }
            Self::LogNormal { median, sigma, .. } => {
                if !pos(median) {
                    bail!("median must be positive and finite");
                }
                if !(sigma.is_finite() && sigma > 1.0) {
                    bail!("sigma must be larger than 1");
                }
            }
        }
        if !(min.is_finite() && min >= 0.0) || !pos(max) || max <= min {
            bail!("min and max must satisfy 0 <= min < max");
        }
        if self.cdf(max) - self.cdf(min) <= 1e-12 {
            bail!("[min, max] holds no mass of the distribution");
        }
        Ok(())
    }

    /// The untruncated cumulative mass fraction below `d`.
    pub fn cdf(&self, d: f64) -> f64 {
        match *self {
            Self::RosinRammler { mean, spread, .. } => {
                if d <= 0.0 {
                    0.0
                } else {
                    -math::exp(-math::pow(d / mean, spread)) + 1.0
                }
            }
            Self::LogNormal { median, sigma, .. } => {
                if d <= 0.0 {
                    0.0
                } else {
                    let z = math::ln(d / median) / math::ln(sigma);
                    normal_cdf(z)
                }
            }
        }
    }

    /// The bounds [min, max].
    pub fn bounds(&self) -> (f64, f64) {
        match *self {
            Self::RosinRammler { min, max, .. } | Self::LogNormal { min, max, .. } => (min, max),
        }
    }

    /// The diameter at the truncated cumulative mass fraction `u` in [0, 1).
    pub fn quantile(&self, u: f64) -> f64 {
        let (min, max) = self.bounds();
        let (f0, f1) = (self.cdf(min), self.cdf(max));
        let f = f0 + u * (f1 - f0);
        let d = match *self {
            Self::RosinRammler { mean, spread, .. } => {
                // f = 1 − exp(−(d/X)ⁿ)  →  d = X (−ln(1 − f))^(1/n)
                let g = (1.0 - f).max(f64::MIN_POSITIVE);
                mean * math::pow(-math::ln(g), 1.0 / spread)
            }
            Self::LogNormal { median, sigma, .. } => {
                let z = normal_quantile(f.clamp(1e-300, 1.0 - 1e-16));
                median * math::exp(z * math::ln(sigma))
            }
        };
        d.clamp(min.max(f64::MIN_POSITIVE), max)
    }

    /// The mass-weighted mean of dᵏ over the truncated distribution, by
    /// Gauss-Legendre quadrature of the quantile (for the checks of V28.14).
    pub fn mass_moment(&self, k: f64) -> f64 {
        // ∫₀¹ d(u)ᵏ du on 4000 midpoints of u: smooth, exact enough for the
        // sampling checks (relative 1e-7)
        let n = 4000;
        let mut s = 0.0;
        for i in 0..n {
            let u = (i as f64 + 0.5) / n as f64;
            s += math::pow(self.quantile(u), k);
        }
        s / n as f64
    }
}

/// Φ(z), the standard normal distribution.
pub fn normal_cdf(z: f64) -> f64 {
    0.5 * math::erfc(-z / std::f64::consts::SQRT_2)
}

/// Φ⁻¹(p) for 0 < p < 1 (Acklam 2003 with one Halley step).
pub fn normal_quantile(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.38357751867269e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    const LOW: f64 = 0.02425;
    let tail = |q: f64| {
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    };
    let x = if p < LOW {
        tail(math::sqrt(-2.0 * math::ln(p)))
    } else if p <= 1.0 - LOW {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        -tail(math::sqrt(-2.0 * math::ln(1.0 - p)))
    };
    // Halley step
    let e = normal_cdf(x) - p;
    let u = e * math::sqrt(2.0 * std::f64::consts::PI) * math::exp(0.5 * x * x);
    x - u / (1.0 + 0.5 * x * u)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_quantile_inverts_the_cdf() {
        for &p in &[
            1e-12,
            1e-6,
            0.01,
            0.02425,
            0.1,
            0.3,
            0.5,
            0.7,
            0.9,
            0.99,
            1.0 - 1e-9,
        ] {
            let z = normal_quantile(p);
            let back = normal_cdf(z);
            assert!(((back - p) / p).abs() < 1e-12, "p {p}: Φ(Φ⁻¹) = {back}");
        }
        assert!(normal_quantile(0.5).abs() < 1e-15);
        assert!((normal_quantile(0.975) - 1.959963984540054).abs() < 1e-12);
    }

    #[test]
    fn quantiles_of_rosin_rammler_follow_the_formula() {
        let s = SizeDistribution::RosinRammler {
            mean: 100e-6,
            spread: 3.0,
            min: 0.0,
            max: 1.0,
        };
        // untruncated: mass fraction below X is 1 − 1/e
        let u = 1.0 - (-1.0f64).exp();
        assert!((s.quantile(u) / 100e-6 - 1.0).abs() < 1e-12);
        let ln = SizeDistribution::LogNormal {
            median: 50e-6,
            sigma: 1.8,
            min: 0.0,
            max: 1.0,
        };
        assert!((ln.quantile(0.5) / 50e-6 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn truncation_keeps_the_bounds() {
        let s = SizeDistribution::LogNormal {
            median: 50e-6,
            sigma: 2.0,
            min: 20e-6,
            max: 120e-6,
        };
        s.validate().unwrap();
        assert!((s.quantile(0.0) - 20e-6).abs() < 1e-18);
        assert!(s.quantile(0.999_999_999) <= 120e-6);
        let mut last = 0.0;
        for i in 0..100 {
            let d = s.quantile(i as f64 / 100.0);
            assert!(d >= last);
            last = d;
        }
    }
}
