//! The adiabatic mixing line of n-dodecane and nitrogen at Spray A
//! conditions (SA1, D319, checks VA.4 and VA.5): the data of the gallery
//! image.
//!
//! ```text
//! cargo run --release -p foehncraft-particles --example spray_a_mixing_line -- output/spray_a_mixing_line
//! python3 -m tools.plot.spray_a_mixing_line output/spray_a_mixing_line
//! ```
//!
//! Writes into the output directory
//!
//! - `line.csv`: the mixing temperature, the liquid mass fraction and the
//!   density over the mixture fraction of 363 K n-dodecane in 900 K
//!   nitrogen at 6 MPa, frozen (one phase of the equation of state),
//!   equilibrium (Peng-Robinson flash) and ideal (ideal gas, Raoult's law),
//!   with the energy-balance error of every point;
//! - `liquid_length.csv`: the saturation mixture fraction of each model
//!   and the liquid length of Siebers (1999) from it at the measured
//!   conditions, with the measured liquid lengths.

use std::fmt::Write as _;

use anyhow::{Context, Result};
use foehncraft_kernels::realfluid::{MAX, Mixing, Mixture, Streams, fuels::N_DODECANE, gases};
use foehncraft_particles::{ideal_gas_species, siebers_liquid_length};

/// A measured Spray A liquid length: ambient temperature (K) and density
/// (kg/m³), liquid length (m); diffused back-illumination, injector
/// 210675, 150 MPa, 0 % O₂ (ECN, ecn.sandia.gov/data/dbi675/).
const MEASURED: [(f64, f64, f64); 4] = [
    (900.0, 22.8, 11.7e-3),
    (700.0, 22.8, 16.3e-3),
    (1200.0, 22.8, 7.4e-3),
    (900.0, 7.6, 20.6e-3),
];
/// Liquid density of n-dodecane at 363.15 K and 6 MPa (kg/m³), Lemmon and
/// Huber (2004), NIST WebBook fluid tables.
const RHO_FUEL: f64 = 703.71;
/// Nozzle of Spray A: diameter (m) and area contraction, Pickett et al.,
/// SAE 2011-01-0901 (C_a = 0.95).
const D: f64 = 0.090e-3;
const CA: f64 = 0.95;

fn main() -> Result<()> {
    let out = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "output/spray_a_mixing_line".into());
    std::fs::create_dir_all(&out).with_context(|| format!("create {out}"))?;
    let mut mix = Mixture::new(&[N_DODECANE.component, gases::N2]);
    mix.set_k(0, 1, N_DODECANE.k_with("N2"));
    let n2 = ideal_gas_species(
        include_str!("../../foehncraft-flow/data/species-gri30.toml"),
        "N2",
    )?;
    let fuel = N_DODECANE.ideal_gas;
    let ideal = |i: usize, t: f64| {
        let s = if i == 0 { &fuel } else { &n2 };
        foehncraft_kernels::species::R * t * s.h_over_rt(t)
    };
    let mut gas = [0.0; MAX];
    gas[1] = 1.0;
    let line = Streams {
        mix: &mix,
        ideal: &ideal,
        p: 6e6,
        t_fuel: 363.0,
        t_gas: 900.0,
        gas,
    };
    let models = [
        ("frozen", Mixing::Frozen),
        ("equilibrium", Mixing::Equilibrium),
        ("ideal", Mixing::Ideal),
    ];
    let mut s = String::from("model,Z,T,liquid,rho,split,residual\n");
    let mut worst: f64 = 0.0;
    for (name, model) in models {
        let (hf, hg) = line.stream_enthalpies(model);
        for k in 0..=200 {
            let z = f64::from(k) / 200.0;
            let pt = line.point(z, model);
            worst = worst.max(pt.residual.abs() / (hg - hf).abs());
            writeln!(
                s,
                "{name},{z},{},{:e},{},{},{:e}",
                pt.t,
                pt.liquid,
                pt.rho,
                u8::from(pt.split),
                pt.residual
            )?;
        }
    }
    std::fs::write(format!("{out}/line.csv"), s)?;
    println!("largest energy-balance error over (h_gas - h_fuel): {worst:.1e}");

    let mut l = String::from("model,T_gas,rho_gas,p_MPa,Z_sat,B,L_mm,theta_deg,measured_mm\n");
    for (t_gas, rho_gas, measured) in MEASURED {
        // the pressure of nitrogen at the measured density
        let mut x = [0.0; MAX];
        x[1] = 1.0;
        let p = mix.pressure(&x, t_gas, gases::N2.molar_mass / rho_gas);
        let at = Streams { p, t_gas, ..line };
        for (name, model) in [
            ("equilibrium", Mixing::Equilibrium),
            ("ideal", Mixing::Ideal),
        ] {
            let z = at.saturation(model).context("no saturation")?;
            for theta in [21.5, 18.7] {
                let len = siebers_liquid_length(z, RHO_FUEL, rho_gas, CA, D, theta);
                writeln!(
                    l,
                    "{name},{t_gas},{rho_gas},{},{z},{},{},{theta},{}",
                    p / 1e6,
                    z / (1.0 - z),
                    len * 1e3,
                    measured * 1e3
                )?;
            }
            println!(
                "{t_gas:6.0} K {rho_gas:5.1} kg/m3 ({:.2} MPa) {name:12}: Z_sat {z:.4}",
                p / 1e6
            );
        }
    }
    std::fs::write(format!("{out}/liquid_length.csv"), l)?;
    println!("wrote {out}/line.csv, liquid_length.csv");
    Ok(())
}
