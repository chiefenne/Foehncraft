//! The drag and the evaporation of a distorted drop with the Corrected
//! Distortion model (Nguyen, Dahms, Pickett and Tagliante 2022; ECN Spray
//! A, SA4, D316): the data of the factor chart of the gallery.
//!
//! ```text
//! cargo run --release -p foehncraft-particles --example spray_a_distortion -- output/spray_a_distortion
//! python3 -m tools.plot.spray_a_distortion docs/figures/gallery output/spray_a_distortion RUN_DYNAMIC RUN_CORRECTED
//! ```
//!
//! Writes `factors.csv` into the output directory: over the TAB distortion
//! y from −1 (prolate) to 1 (flattened), at the drop Reynolds numbers 20,
//! 100 and 500 and the viscosity ratios λ = μ_l/μ_g of 2.5, 10 and 40, the
//! drag of the Corrected Distortion model over that of the solid sphere,
//! the dynamic drag 1 + 2.632 y of the standard model (Liu, Mather and
//! Reitz 1993) and the Sherwood factor at Sc = 2 (the film of a Spray A
//! drop, n-dodecane in nitrogen) with its parts: the internal circulation
//! of the viscous sphere, the shape ratio and the surface ratio.

use std::fmt::Write as _;

use anyhow::{Context, Result};
use foehncraft_kernels::breakup::dynamic_drag;
use foehncraft_kernels::distortion as cdm;

fn main() -> Result<()> {
    let out = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "output/spray_a_distortion".into());
    std::fs::create_dir_all(&out).with_context(|| format!("cannot create {out}"))?;
    let sc = 2.0;
    let mut s = String::from(
        "re,lambda,y,drag_corrected,drag_dynamic,sherwood,circulation,shape,surface\n",
    );
    for re in [20.0, 100.0, 500.0] {
        for lambda in [2.5, 10.0, 40.0] {
            for k in 0..=80 {
                let y = -1.0 + 2.0 * f64::from(k) / 80.0;
                let g = cdm::shape(y, cdm::C_B);
                writeln!(
                    s,
                    "{re},{lambda},{y},{},{},{},{},{},{}",
                    cdm::drag_factor(re, lambda, y),
                    dynamic_drag(y),
                    cdm::transfer_factor(re, sc, lambda, y),
                    cdm::transfer_viscous_ratio(re, sc, lambda),
                    cdm::transfer_shape_ratio(re, sc, g.phi, g.phi_perp),
                    g.area,
                )?;
            }
        }
    }
    let path = format!("{out}/factors.csv");
    std::fs::write(&path, s).with_context(|| format!("cannot write {path}"))?;
    eprintln!("wrote {path}");
    Ok(())
}
