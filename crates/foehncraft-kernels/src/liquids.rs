//! Properties of the spray liquids (D320): pure single-component liquids on
//! their saturation line, as plain functions of the temperature (no
//! allocation, `f64` in and out) for the parcel kernels on CPU and GPU.
//!
//! The coefficients are our own fits (written by
//! `tools/data/liquid_properties.py`) to the saturation tables of the NIST
//! Chemistry WebBook (NIST Standard Reference Database 69, "Thermophysical
//! Properties of Fluid Systems", E. W. Lemmon, M. O. McLinden and D. G.
//! Friend); the tables themselves are not in the repository. With
//! t = 1 − T/T_c and T_r = T/T_c:
//!
//! - vapour pressure (Wagner form): ln(p/p_c) = (a t + b t^1.5 + c t^2.5 + d t^5)/T_r;
//! - density: ρ = ρ_c (1 + Σ a_k t^e_k), e = 0.35, 2/3, 1, 4/3, 5/3, 16/3;
//! - specific heat: c_l = a + b/t + c t + d t²;
//! - latent heat (Watson form): L = t^0.38 (a + b t + c t² + d t³);
//! - surface tension: σ = t^1.26 (a + b t + c t² + d t³);
//! - viscosity: ln μ = Σ a_k T_r^e_k, e = 0, −1, −2, −3, −4, 1.
//!
//! The largest relative error of each fit is written with each liquid.
//! Below `t_min` the functions are evaluated at `t_min`; the specific heat
//! and the viscosity are fitted up to 0.95 T_c, density, latent heat and
//! surface tension up to 0.98 T_c, the vapour pressure to T_c, and every
//! function is evaluated at 0.999 T_c at most.
//!
//! These are properties of the pure liquid at its own vapour pressure:
//! with Raoult's law (an ideal gas over the droplet) they hold at low
//! pressure. At engine pressures the dissolved gas and the real-gas phase
//! equilibrium matter (the high-pressure layer of the ECN Spray A work,
//! D310 to D319), which replaces the surface vapour fraction, not these
//! pure-liquid functions.

use foehncraft_core::Real;

/// The data of one liquid: critical constants and fit coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LiquidData {
    /// Name of the liquid (`[spray.liquid] library`).
    pub name: &'static str,
    /// Molar mass (kg/mol).
    pub molar_mass: f64,
    /// Critical temperature (K), of the NIST table.
    pub tc: f64,
    /// Critical pressure (Pa), of the NIST table.
    pub pc: f64,
    /// Critical density (kg/m³), of the NIST table.
    pub rho_c: f64,
    /// Critical molar volume (m³/mol), M/ρ_c.
    pub vc: f64,
    /// Critical compressibility p_c V_c/(R T_c).
    pub zc: f64,
    /// Acentric factor −log₁₀(p_sat(0.7 T_c)/p_c) − 1 of the vapour-pressure fit.
    pub omega: f64,
    /// Normal boiling point (K) of the vapour-pressure fit.
    pub tb: f64,
    /// Lowest temperature of the fits (K, the triple point of the table).
    pub t_min: f64,
    /// Vapour pressure a, b, c, d.
    pub psat: [f64; 4],
    /// Density a₁ … a₆.
    pub density: [f64; 6],
    /// Specific heat a, b, c, d (J/(kg K)).
    pub cp: [f64; 4],
    /// Latent heat a, b, c, d (J/kg).
    pub latent: [f64; 4],
    /// Surface tension a, b, c, d (N/m).
    pub sigma: [f64; 4],
    /// Viscosity a₁ … a₆ (ln of Pa s).
    pub viscosity: [f64; 6],
}

const DENSITY_EXP: [f64; 6] = [0.35, 2.0 / 3.0, 1.0, 4.0 / 3.0, 5.0 / 3.0, 16.0 / 3.0];
const VISCOSITY_EXP: [f64; 6] = [0.0, -1.0, -2.0, -3.0, -4.0, 1.0];

#[inline(always)]
fn powf(x: f64, e: f64) -> f64 {
    (e * x.ln_of()).exp_of()
}

impl LiquidData {
    /// t = 1 − T/T_c and T_r at `t`, the temperature kept inside the fits.
    #[inline(always)]
    fn reduced(&self, t: f64) -> (f64, f64) {
        let t = t.max_of(self.t_min).min_of(0.999 * self.tc);
        let tr = t / self.tc;
        (1.0 - tr, tr)
    }

    /// Saturation (vapour) pressure (Pa) at `t` (K); p_c at and above T_c.
    #[inline]
    pub fn saturation_pressure(&self, t: f64) -> f64 {
        if t >= self.tc {
            return self.pc;
        }
        let (th, tr) = self.reduced(t);
        let a = &self.psat;
        let s = th.sqrt_of();
        let th2 = th * th;
        let sum = a[0] * th + a[1] * th * s + a[2] * th2 * s + a[3] * th2 * th2 * th;
        self.pc * (sum / tr).exp_of()
    }

    /// Density of the saturated liquid (kg/m³) at `t` (K).
    #[inline]
    pub fn density(&self, t: f64) -> f64 {
        let (th, _) = self.reduced(t);
        let mut s = 1.0;
        for k in 0..6 {
            s += self.density[k] * powf(th, DENSITY_EXP[k]);
        }
        self.rho_c * s
    }

    /// Specific heat of the saturated liquid (J/(kg K)) at `t` (K).
    #[inline]
    pub fn specific_heat(&self, t: f64) -> f64 {
        let (th, _) = self.reduced(t);
        let th = th.max_of(0.05);
        let a = &self.cp;
        a[0] + a[1] / th + a[2] * th + a[3] * th * th
    }

    /// Latent heat of vaporisation (J/kg) at `t` (K); 0 at and above T_c.
    #[inline]
    pub fn latent_heat(&self, t: f64) -> f64 {
        if t >= self.tc {
            return 0.0;
        }
        let (th, _) = self.reduced(t);
        let a = &self.latent;
        powf(th, 0.38) * (a[0] + th * (a[1] + th * (a[2] + th * a[3])))
    }

    /// Surface tension (N/m) at `t` (K); 0 at and above T_c.
    #[inline]
    pub fn surface_tension(&self, t: f64) -> f64 {
        if t >= self.tc {
            return 0.0;
        }
        let (th, _) = self.reduced(t);
        let a = &self.sigma;
        powf(th, 1.26) * (a[0] + th * (a[1] + th * (a[2] + th * a[3])))
    }

    /// Dynamic viscosity of the saturated liquid (Pa s) at `t` (K).
    #[inline]
    pub fn viscosity(&self, t: f64) -> f64 {
        let (th, tr) = self.reduced(t);
        let _ = th;
        let tr = tr.min_of(0.95);
        let mut s = 0.0;
        for k in 0..6 {
            s += self.viscosity[k] * powf(tr, VISCOSITY_EXP[k]);
        }
        s.exp_of()
    }
}

/// The liquids of the library, by the name of `[spray.liquid] library`.
pub fn by_name(name: &str) -> Option<&'static LiquidData> {
    ALL.iter().copied().find(|l| l.name == name)
}

/// Every liquid of the library.
pub const ALL: [&LiquidData; 3] = [&WATER, &N_HEPTANE, &N_DODECANE];

// The data below is written by tools/data/liquid_properties.py from the NIST WebBook
// saturation tables (SRD 69); do not edit by hand.
/// water: largest relative errors of the fits against the NIST table:
/// vapour pressure 6.9e-04 to 0.9 Tc, 6.9e-04 to 647 K.
/// density 9.4e-04 to 0.9 Tc, 9.4e-04 to 634 K.
/// specific heat (Tr < 0.95) 1.7e-03 to 0.9 Tc, 1.7e-03 to 615 K.
/// latent heat 9.6e-03 to 0.9 Tc, 9.6e-03 to 634 K.
/// surface tension 2.0e-03 to 0.9 Tc, 2.0e-03 to 634 K.
/// viscosity (Tr < 0.95) 9.5e-03 to 0.9 Tc, 9.5e-03 to 615 K.
#[allow(clippy::excessive_precision)]
pub const WATER: LiquidData = LiquidData {
    name: "water",
    molar_mass: 0.01801528,
    tc: 647.096,
    pc: 22064000.0,
    rho_c: 322.000000,
    vc: 5.594807e-05,
    zc: 0.229439,
    omega: 0.344305,
    tb: 373.1228,
    t_min: 273.179,
    psat: [-7.8563650847e+00, 1.8606030329e+00, -2.2386050055e+00, -2.1524913600e+00],
    density: [2.2128392741e+00, 9.5908058297e-01, -1.0301346641e+00, 1.9980404191e+00, -1.5983718763e+00, -2.0730607339e+00],
    cp: [3.8612249115e+03, 2.3210357941e+02, -1.5419314736e+03, 2.5603785651e+03],
    latent: [3.0641516310e+06, 1.9278376845e+06, -7.1821202092e+06, 6.8931853632e+06],
    sigma: [2.3958691563e-01, -1.6954867782e-01, 5.8638378155e-02, -5.4583644193e-02],
    viscosity: [1.4123940799e+02, -1.9029023379e+02, 1.1936679900e+02, -3.6712655786e+01, 4.4967876724e+00, -4.7857194958e+01],
};
/// n-heptane: largest relative errors of the fits against the NIST table:
/// vapour pressure 2.6e-03 to 0.9 Tc, 2.6e-03 to 540 K.
/// density 2.0e-04 to 0.9 Tc, 2.0e-04 to 529 K.
/// specific heat (Tr < 0.95) 1.5e-02 to 0.9 Tc, 1.5e-02 to 513 K.
/// latent heat 9.1e-03 to 0.9 Tc, 9.1e-03 to 529 K.
/// surface tension 9.5e-03 to 0.9 Tc, 9.5e-03 to 529 K.
/// viscosity (Tr < 0.95) 1.3e-02 to 0.9 Tc, 1.3e-02 to 513 K.
#[allow(clippy::excessive_precision)]
pub const N_HEPTANE: LiquidData = LiquidData {
    name: "n-heptane",
    molar_mass: 0.100202,
    tc: 540.2,
    pc: 2735760.0,
    rho_c: 233.471000,
    vc: 4.291839e-04,
    zc: 0.261416,
    omega: 0.349216,
    tb: 371.5638,
    t_min: 182.55,
    psat: [-7.7839490985e+00, 1.8822945254e+00, -2.8459050921e+00, -3.4972770672e+00],
    density: [2.1167244226e+00, -2.4056596457e-02, 1.6222523127e+00, -2.3295453632e+00, 1.5103409153e+00, 1.7522719171e-01],
    cp: [3.1550972887e+03, 4.0842088455e+01, -2.7249691818e+03, 1.2553109560e+03],
    latent: [4.6943733505e+05, 2.5522413925e+05, -8.6203069664e+05, 8.9925318973e+05],
    sigma: [5.5046170197e-02, 1.3375121376e-02, -6.2222953589e-02, 5.8826107660e-02],
    viscosity: [1.3203078790e+02, -1.6884677266e+02, 1.0097705493e+02, -2.9396695474e+01, 3.3366035836e+00, -4.8114926137e+01],
};
/// n-dodecane: largest relative errors of the fits against the NIST table:
/// vapour pressure 2.7e-03 to 0.9 Tc, 2.7e-03 to 658 K.
/// density 3.5e-04 to 0.9 Tc, 3.5e-04 to 645 K.
/// specific heat (Tr < 0.95) 5.8e-03 to 0.9 Tc, 5.8e-03 to 625 K.
/// latent heat 3.3e-03 to 0.9 Tc, 4.6e-03 to 645 K.
/// surface tension 3.6e-02 to 0.9 Tc, 3.6e-02 to 645 K.
/// viscosity (Tr < 0.95) 1.4e-02 to 0.9 Tc, 1.4e-02 to 625 K.
#[allow(clippy::excessive_precision)]
pub const N_DODECANE: LiquidData = LiquidData {
    name: "n-dodecane",
    molar_mass: 0.170335,
    tc: 657.957,
    pc: 1813970.0,
    rho_c: 230.423000,
    vc: 7.392274e-04,
    zc: 0.245119,
    omega: 0.574880,
    tb: 489.4521,
    t_min: 263.6,
    psat: [-8.9969852229e+00, 2.4854924713e+00, -4.7738684248e+00, -5.0992392466e+00],
    density: [3.1079944383e+00, -7.2231182966e+00, 2.2232448651e+01, -2.7195526377e+01, 1.2496513338e+01, -3.4089484162e-01],
    cp: [3.5069960956e+03, 2.3266841502e+01, -2.5386476771e+03, 2.2036152784e+02],
    latent: [3.8973673142e+05, 3.2541121968e+05, -8.9420165030e+05, 9.3639435675e+05],
    sigma: [6.9395132876e-02, -1.3364163449e-01, 3.8265751460e-01, -3.4829386370e-01],
    viscosity: [3.5820097678e+02, -4.6239907431e+02, 2.9063562653e+02, -9.0073493199e+01, 1.1006906645e+01, -1.1761633174e+02],
};

#[cfg(test)]
mod tests {
    use super::*;

    /// (T K, p Pa, ρ kg/m³, c_l J/(kg K), L J/kg, σ N/m, μ Pa s): points of
    /// the NIST WebBook saturation tables (SRD 69).
    const NIST: [(&LiquidData, [f64; 7]); 6] = [
        (&WATER, [347.39, 37384.8, 975.267, 4192.84, 2322477.0, 0.0637203, 0.000381182]),
        (&WATER, [449.433, 919920.0, 890.938, 4390.58, 2027225.0, 0.0430174, 0.000153739]),
        (&N_HEPTANE, [306.491, 9119.27, 672.387, 2278.86, 359792.0, 0.0189135, 0.000356591]),
        (&N_HEPTANE, [449.87, 642897.0, 527.145, 2979.67, 249115.0, 0.00585302, 0.000105924]),
        (&N_DODECANE, [381.888, 3028.96, 683.244, 2521.56, 315942.7, 0.0177905, 0.000466005]),
        (&N_DODECANE, [450.57, 36340.6, 628.401, 2801.09, 279149.0, 0.0125638, 0.000261062]),
    ];

    #[test]
    fn the_fits_meet_the_nist_tables() {
        for (l, [t, p, rho, cp, lat, sig, mu]) in NIST {
            let e = |f: f64, r: f64| (f / r - 1.0).abs();
            assert!(e(l.saturation_pressure(t), p) < 3e-3, "{} p {t}", l.name);
            assert!(e(l.density(t), rho) < 2e-3, "{} rho {t}", l.name);
            assert!(e(l.specific_heat(t), cp) < 1.5e-2, "{} cp {t}", l.name);
            assert!(e(l.latent_heat(t), lat) < 1e-2, "{} L {t}", l.name);
            assert!(e(l.surface_tension(t), sig) < 4e-2, "{} sigma {t}", l.name);
            assert!(e(l.viscosity(t), mu) < 1.5e-2, "{} mu {t}", l.name);
        }
    }

    #[test]
    fn the_normal_boiling_points_and_critical_ends() {
        // NIST WebBook: T_b 373.12 K (water), 371.5 K (n-heptane), 489.4 K (n-dodecane)
        for (l, tb) in [(&WATER, 373.12), (&N_HEPTANE, 371.5), (&N_DODECANE, 489.4)] {
            assert!((l.tb - tb).abs() < 0.5, "{} {}", l.name, l.tb);
            assert!((l.saturation_pressure(l.tb) / 101_325.0 - 1.0).abs() < 1e-4);
            assert_eq!(l.saturation_pressure(l.tc), l.pc);
            assert_eq!(l.latent_heat(l.tc + 1.0), 0.0);
            assert!(l.latent_heat(0.99 * l.tc) > 0.0 && l.surface_tension(0.99 * l.tc) > 0.0);
        }
        // the water fit agrees with the IAPWS-95 auxiliary equation of the evaporation
        for t in [280.0, 320.0, 373.124, 450.0, 600.0] {
            let r = WATER.saturation_pressure(t) / crate::evaporation::water::saturation_pressure(t);
            assert!((r - 1.0).abs() < 2e-3, "{t} {r}");
        }
        assert!(by_name("n-heptane").is_some() && by_name("benzene").is_none());
    }

    #[test]
    fn the_functions_are_smooth_and_monotone_where_they_must_be() {
        for l in ALL {
            let mut prev = (0.0, f64::INFINITY, f64::INFINITY);
            let n = 400;
            for i in 0..=n {
                let t = l.t_min + (0.97 * l.tc - l.t_min) * i as f64 / n as f64;
                let now = (l.saturation_pressure(t), l.density(t), l.surface_tension(t));
                assert!(now.0 > prev.0, "{} p at {t}", l.name);
                assert!(now.1 < prev.1, "{} rho at {t}", l.name);
                assert!(now.2 < prev.2, "{} sigma at {t}", l.name);
                prev = now;
            }
        }
    }
}
