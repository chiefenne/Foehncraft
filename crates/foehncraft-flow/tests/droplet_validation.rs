//! Single-droplet evaporation against measurements (D320 validation):
//! n-heptane droplets in still nitrogen at 0.1 MPa, Chauveau, Birouk,
//! Halter and Gökalp (2019), the microgravity rates K0(T) of their eq. (4).

use foehncraft_flow::species::database;
use foehncraft_kernels::evaporation::{self as ev, EvapInput};
use foehncraft_kernels::liquids::N_HEPTANE;
use foehncraft_kernels::species as ks;

const P0: f64 = 101_325.0;
const R: f64 = 8.314_462_618;

/// K (mm²/s) of a droplet of d0 at t0 in still N2 at tg, integrated with the
/// kernel; `rule`: the film state Tp + rule (Tg − Tp) (0 = droplet, 1 = far field).
fn rate(tg: f64, d0: f64, t0: f64, rule: f64, model: u32) -> (f64, f64) {
    let db = database();
    let n2 = db["N2"].kernel();
    let hep = db["NC7H16"].kernel();
    let (mut d, mut tp) = (d0, t0);
    let mut t = 0.0;
    let mut hist = Vec::new();
    let rho_l = N_HEPTANE.density(t0);
    let cp_l = N_HEPTANE.specific_heat(t0);
    while d > 0.05 * d0 {
        let psat = N_HEPTANE.saturation_pressure(tp);
        let ys = ev::surface_fraction(psat, P0, hep.molar_mass, n2.molar_mass);
        let tf = tp + rule * (tg - tp);
        let yf = ys * (1.0 - rule);
        let sp = [hep, n2];
        let mix = ks::mixture(&sp, &[yf, 1.0 - yf], tf);
        let rho = P0 * mix.molar_mass / (R * tf);
        let diff = ks::binary_diffusion(&hep, &n2, tf, P0);
        let dt = 1e-4 * (d / d0).powi(2) + 1e-7;
        let s = ev::evaporation_step(&EvapInput {
            tp,
            tg,
            yg: 0.0,
            p0: P0,
            m_vapour: hep.molar_mass,
            m_gas: n2.molar_mass,
            psat,
            ys_surface: -1.0,
            latent: N_HEPTANE.latent_heat(tp),
            re: 0.0,
            pr: mix.viscosity * mix.cp / mix.conductivity,
            sc: mix.viscosity / (rho * diff),
            rho_g: rho,
            diff,
            k_g: mix.conductivity,
            cp_g: mix.cp,
            cp_v: hep.cp(tf),
            d,
            rho_l,
            cp_l,
            dt,
            model,
        });
        d = s.d;
        tp = s.tp;
        t += dt;
        hist.push((t, (d / d0).powi(2), tp));
        if d == 0.0 {
            break;
        }
    }
    // the slope of d²/d0² between 0.8 and 0.3 (the quasi-steady part)
    let a = hist.iter().find(|h| h.1 < 0.8).unwrap();
    let b = hist.iter().find(|h| h.1 < 0.3).unwrap();
    let k = (a.1 - b.1) / (b.0 - a.0) * d0 * d0 * 1e6;
    (k, b.2)
}

#[test]
fn n_heptane_droplets_against_chauveau_2019() {
    for tg in [473.0, 548.0, 623.0, 698.0, 748.0] {
        let k0 = 3.6552e-4 * tg - 0.1078; // mm²/s, eq. (4), microgravity, 14 µm cross fibre
        for (name, rule, m) in [
            ("far field, RM", 1.0, ev::model::RANZ_MARSHALL),
            ("1/3 rule, RM", 1.0 / 3.0, ev::model::RANZ_MARSHALL),
        ] {
            let (k, tw) = rate(tg, 370e-6, 300.0, rule, m);
            eprintln!("T {tg} K {name}: K {k:.4} mm2/s, measured {k0:.4} ({:+.1} %), wet bulb {tw:.1} K", 100.0 * (k / k0 - 1.0));
        }
    }
}

#[test]
fn print_heptane_diffusivity() {
    let db = database();
    let (h, n) = (db["NC7H16"].kernel(), db["N2"].kernel());
    for t in [300.0, 500.0, 748.0] {
        eprintln!("D {t}: {:.4e} m2/s, mu_N2 {:.3e}, k_N2 {:.4}", ks::binary_diffusion(&h, &n, t, P0), n.viscosity(t), n.conductivity(t, n.viscosity(t)));
    }
}
