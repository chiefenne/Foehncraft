//! The Peng-Robinson surface of an n-dodecane droplet in nitrogen (ECN
//! Spray A, SA1, D311, D312): the data of the gallery image.
//!
//! ```text
//! cargo run --release -p foehncraft-particles --example spray_a_surface -- output/spray_a_surface
//! python3 -m tools.plot.spray_a_surface output/spray_a_surface
//! ```
//!
//! Writes into the output directory
//!
//! - `surface.csv`: the surface tables at 1, 3, 6 and 8 MPa (fuel vapour
//!   fraction, its Raoult value y = p_sat/p with the same saturation
//!   pressure, dissolved nitrogen, latent heat ratio) with the end of each
//!   table (critical mixing or boiling);
//! - `bubble.csv`: the bubble pressure over the dissolved nitrogen on the
//!   seven isotherms of García-Córdova et al. (2011) up to 20 MPa;
//! - `measured.csv`: the 21 bubble points of that paper that the test
//!   `tests/high_pressure.rs` quotes (doi 10.1021/je1012372).

use std::fmt::Write as _;

use anyhow::{Context, Result};
use foehncraft_kernels::realfluid::{self as rf, MAX, Mixture, fuels::N_DODECANE, gases};
use foehncraft_particles::{SurfaceTable, TableEnd};

/// (T in K, p in kPa, x_N2): table 2 of García-Córdova, Justo-García,
/// García-Flores and García-Sánchez (2011), J. Chem. Eng. Data 56,
/// 1555-1564, the subset of the test.
const MEASURED: [(f64, f64, f64); 21] = [
    (344.4, 1430.0, 0.0204),
    (344.4, 5030.0, 0.0669),
    (344.4, 12380.0, 0.1435),
    (410.7, 2170.0, 0.0339),
    (410.7, 5020.0, 0.0754),
    (410.7, 12980.0, 0.1698),
    (463.9, 2050.0, 0.0372),
    (463.9, 5040.0, 0.0875),
    (463.9, 12500.0, 0.1945),
    (503.4, 2160.0, 0.0436),
    (503.4, 5310.0, 0.1041),
    (503.4, 12440.0, 0.2213),
    (532.9, 2030.0, 0.0423),
    (532.9, 6010.0, 0.1271),
    (532.9, 12450.0, 0.2462),
    (562.1, 2130.0, 0.0435),
    (562.1, 5960.0, 0.1315),
    (562.1, 11970.0, 0.2582),
    (593.5, 2270.0, 0.0504),
    (593.5, 5070.0, 0.1308),
    (593.5, 10960.0, 0.2837),
];

fn main() -> Result<()> {
    let out = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "output/spray_a_surface".into());
    std::fs::create_dir_all(&out).with_context(|| format!("create {out}"))?;
    let mut mix = Mixture::new(&[N_DODECANE.component, gases::N2]);
    mix.set_k(0, 1, N_DODECANE.k_with("N2"));
    let fuel = Mixture::new(&[N_DODECANE.component]);
    let mut r = [0.0; MAX];
    r[1] = 1.0;

    let mut s = String::from("p_MPa,T,y_fuel,y_raoult,x_N2,latent_ratio,end\n");
    for p in [1e6, 3e6, 6e6, 8e6] {
        let tab = SurfaceTable::build(&mix, &r, p, 300.0, 700.0, 1.0)?;
        let end = match tab.end {
            TableEnd::Critical => "critical",
            TableEnd::Boiling => "boiling",
            TableEnd::Range => "range",
        };
        for i in 0..tab.t.len() {
            let t = tab.t[i];
            let raoult = fuel.saturation_pressure(0, t).map_or(f64::NAN, |ps| ps / p);
            writeln!(
                s,
                "{},{t},{:e},{raoult:e},{:e},{:e},{}",
                p / 1e6,
                tab.y_fuel[i],
                tab.dissolved[i],
                tab.latent_ratio[i],
                if i + 1 == tab.t.len() { end } else { "" }
            )?;
        }
        println!(
            "{:4.1} MPa: two-phase surface ends at {:.2} K ({end})",
            p / 1e6,
            tab.t_end
        );
    }
    std::fs::write(format!("{out}/surface.csv"), s)?;

    let mut b = String::from("T,x_N2,p_kPa\n");
    for t in [344.4, 410.7, 463.9, 503.4, 532.9, 562.1, 593.5] {
        let mut p = 1e6;
        let mut xn = 0.002;
        while xn < 0.40 {
            let mut x = [0.0; MAX];
            x[0] = 1.0 - xn;
            x[1] = xn;
            match rf::bubble_pressure(&mix, t, &x, p) {
                Some((pb, _)) if pb <= 20.5e6 => {
                    writeln!(b, "{t},{xn},{}", pb / 1e3)?;
                    p = pb;
                }
                _ => break,
            }
            xn += 0.002;
        }
    }
    std::fs::write(format!("{out}/bubble.csv"), b)?;

    let mut m = String::from("T,p_kPa,x_N2\n");
    for (t, p, x) in MEASURED {
        writeln!(m, "{t},{p},{x}")?;
    }
    std::fs::write(format!("{out}/measured.csv"), m)?;
    println!("wrote {out}/surface.csv, bubble.csv, measured.csv");
    Ok(())
}
