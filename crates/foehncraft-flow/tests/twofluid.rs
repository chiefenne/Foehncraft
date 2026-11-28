//! The two-fluid (Euler–Euler) model (L3, stage MP3, D352 to D356): the
//! drag balance of a uniform bubble column, the convergence the partial
//! elimination of the drag gives, and the refusals.

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_twofluid_{name}_{}", std::process::id()))
}

/// The columns of a csv file of the run (comment lines skipped).
fn table(path: &Path) -> (Vec<String>, Vec<Vec<f64>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let head = lines
        .next()
        .unwrap()
        .split(',')
        .map(str::to_string)
        .collect();
    let rows = lines
        .map(|l| {
            l.split(',')
                .map(|v| v.parse::<f64>().unwrap_or(f64::NAN))
                .collect()
        })
        .collect();
    (head, rows)
}

/// Air bubbles of 4 mm rise with water in a uniform column: the slip is
/// the drag balance U_T √(1 − α) of the Ishii–Zuber distorted regime,
/// U_T = √2 (g Δρ σ)^¼/√ρ_l; the gas fraction is the inflow's j_G over
/// the gas speed; the steps converge in a few outer iterations.
#[test]
fn bubble_column_drag_balance() {
    let text = std::fs::read_to_string(cases().join("two_fluid_column.toml")).unwrap();
    let case = FlowCase::from_toml(&text, &cases()).unwrap();
    let dir = out("column");
    run(&case, "column", &dir, |_| {}).unwrap();
    let (head, rows) = table(&dir.join("line_axis.csv"));
    let col = |n: &str| head.iter().position(|h| h == n).unwrap();
    let (cy, cv, ca, cg) = (col("y"), col("v"), col("alpha_gas"), col("v_gas"));
    let (rho_l, rho_g, sigma, g) = (998.2_f64, 1.2_f64, 0.072_f64, 9.81_f64);
    let ut = 2.0_f64.sqrt() * (g * (rho_l - rho_g) * sigma).powf(0.25) / rho_l.sqrt();
    // the column, filled after 0.28 s (the front rises at 0.72 m/s)
    let mut n = 0;
    for r in rows.iter().filter(|r| r[cy] > 0.02 && r[cy] < 0.18) {
        let a = r[ca];
        let slip = r[cg] - r[cv];
        let exact = ut * (1.0 - a).sqrt();
        assert!(
            (slip - exact).abs() < 2e-3 * exact,
            "y {:.3}: slip {slip:.5} against {exact:.5} (alpha {a:.4})",
            r[cy]
        );
        // α = j_G/u_g with j_G = 0.05 (0.5 + U_T) of the inflow
        let ag = 0.05 * (0.5 + ut) / r[cg];
        assert!(
            (a - ag).abs() < 0.01 * ag,
            "alpha {a:.5} against j_G/u_g {ag:.5}"
        );
        n += 1;
    }
    assert!(n > 20, "too few points in the filled column ({n})");
    let (head, rows) = table(&dir.join("history.csv"));
    let outer = head.iter().position(|h| h == "outer").unwrap();
    let last = rows.last().unwrap();
    assert!(
        last[outer] <= 3.0,
        "the steps take {} outer iterations",
        last[outer]
    );
}

fn refused(text: &str) -> String {
    let e = match FlowCase::from_toml(text, &cases()) {
        Err(e) => e,
        Ok(c) => match run(&c, "refused", &out("refused"), |_| {}) {
            Err(e) => e,
            Ok(_) => panic!("the case was accepted"),
        },
    };
    diag::of_error(&e).code.to_string()
}

#[test]
fn refusals() {
    // (the line ends of a Windows checkout are CRLF; the edits below
    // match LF)
    let base = std::fs::read_to_string(cases().join("two_fluid_column.toml"))
        .unwrap()
        .replace("\r\n", "\n");
    // steady runs, BDF2 and an adaptive step are not supported
    let steady = base.replace(
        "[time]\nscheme = \"bdf1\"\ntime_step = 2.0e-3\nsteps = 200\n",
        "",
    );
    assert_eq!(refused(&steady), "TWF002");
    let bdf2 = base.replace("scheme = \"bdf1\"", "scheme = \"bdf2\"");
    assert_eq!(refused(&bdf2), "TWF002");
    // a gas without a size, and a lift constant without its coefficient
    let no_d = base.replace("diameter = 4.0e-3\n", "diameter = 0.0\n");
    assert_eq!(refused(&no_d), "TWF001");
    let cl = base.replace("lift = \"none\"", "lift = \"constant\"");
    assert_eq!(refused(&cl), "TWF001");
    // with the volume of fluid
    let vof = format!(
        "{base}\n[vof]\n[[vof.phase]]\nname = \"a\"\ndensity = 1.0\nviscosity = 1.0\n[[vof.phase]]\nname = \"b\"\ndensity = 1.0\nviscosity = 1.0\n"
    );
    assert_eq!(refused(&vof), "TWF002");
}

/// The Bartolomei B-7 tube (MP4) cut to 0.3 m on a coarse mesh, with the
/// liquid entering 15 K below saturation so that it boils at once.
fn boiling_case() -> String {
    std::fs::read_to_string(cases().join("boiling_pipe_bartolomei_b7.toml"))
        .unwrap()
        .replace("\r\n", "\n")
        .replace("upper = [1.5, 0.006, 0.001]", "upper = [0.3, 0.006, 0.001]")
        .replace("cells = [300, 12, 1]", "cells = [30, 12, 1]")
        .replace("xmin = 466.0", "xmin = 542.4")
        .replace("steps = 2000", "steps = 60")
        .replace("fields_every = 500", "fields_every = 60")
}

/// Wall boiling (L3, stage MP4, D357 to D359): the heat the walls take is
/// the prescribed flux, part of it makes vapour at the wall, the wall stays
/// a few kelvin above saturation, and the run's volume account of both
/// phases with the vapour made closes (the run's own check, RUN007).
#[test]
fn wall_boiling_heats_and_makes_vapour() {
    let text = boiling_case();
    let case = FlowCase::from_toml(&text, &cases()).unwrap();
    let dir = out("boiling");
    run(&case, "boiling", &dir, |_| {}).unwrap();
    let (head, rows) = table(&dir.join("two_fluid.csv"));
    let col = |n: &str| head.iter().position(|h| h == n).unwrap();
    let last = rows.last().unwrap();
    // the wall of the 2° wedge: r θ L
    let area = 0.006 * 2.0_f64.to_radians() * 0.3;
    let q = last[col("wall_heat")];
    assert!(
        (q - 1.13e6 * area).abs() < 1e-3 * 1.13e6 * area,
        "wall heat {q} W against {} W",
        1.13e6 * area
    );
    let qe = last[col("wall_evaporation_heat")];
    assert!(qe > 0.0 && qe < q, "evaporation {qe} W of {q} W");
    assert!(last[col("wall_vapour")] > 0.0);
    let t_sat = 557.42;
    let tw = last[col("max_wall_temperature")];
    assert!(
        tw > t_sat && tw < t_sat + 30.0,
        "wall temperature {tw} K (T_sat {t_sat} K)"
    );
    assert!(last[col("volume")] > 0.0, "no vapour in the tube");
    for n in ["energy_iterations_liquid", "energy_iterations_vapour"] {
        assert!(last[col(n)] < 100.0, "{n} = {}", last[col(n)]);
    }
}

#[test]
fn boiling_refusals() {
    let base = boiling_case();
    // without k for the wall functions
    let laminar = base.replace("model = \"k-omega-sst\"", "model = \"laminar\"");
    assert_eq!(refused(&laminar), "TWF001");
    // a heat flux on a boundary that is not a wall
    let inlet = base.replace(
        "[two_fluid.boiling.heat_flux]\nymax",
        "[two_fluid.boiling.heat_flux]\nxmin",
    );
    assert_eq!(refused(&inlet), "TWF001");
    // a vapour that is not saturated steam at the pressure
    let light = base.replace("density = 35.5932", "density = 30.0");
    assert_eq!(refused(&light), "TWF001");
    // with [energy] as well
    let energy = format!("{base}\n[energy]\ninitial = 500.0\n");
    assert_eq!(refused(&energy), "TWF002");
}
