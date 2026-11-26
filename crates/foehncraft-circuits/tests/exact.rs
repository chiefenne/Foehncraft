//! Exact checks of the circuits (SF7b): laminar flow split of parallel
//! tubes, the energy balance of a heated tube, the pressure balance of a
//! single-phase riser, and the closure of a natural circulation loop.

use foehncraft_circuits::case::Case;
use foehncraft_circuits::output::{self, Solution};
use foehncraft_kernels::steam;
use std::f64::consts::PI;
use std::path::{Path, PathBuf};

pub fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/circuits")
}

pub fn solve(name: &str) -> (Case, Solution) {
    let case = Case::from_file(&cases().join(name)).unwrap();
    let sol = output::run(&case).unwrap_or_else(|e| panic!("{name}: {e}"));
    (case, sol)
}

#[test]
fn laminar_parallel_tubes_split_as_hagen_poiseuille() {
    let (_, sol) = solve("parallel_tubes_laminar.toml");
    let (a, b) = (sol.flow("short").unwrap(), sol.flow("long").unwrap());
    // L/D^4 ratio 2: the short tube takes twice the flow of the long one
    assert!((a / b - 2.0).abs() < 1e-6, "split {}", a / b);
    assert!((a + b - 1.5e-4).abs() < 1e-14, "sum {}", a + b);
    // the pressure drop of Hagen-Poiseuille with the inlet viscosity
    let s = steam::state_pt(1e6, 300.0).unwrap();
    let mu = steam::viscosity(s.rho, s.t);
    let d: f64 = 5e-3;
    let dp = 128.0 * mu * 2.0 * a / (PI * s.rho * d.powi(4));
    let dp_num = sol.pressure("inlet").unwrap() - 1e6;
    assert!((dp_num / dp - 1.0).abs() < 1e-6, "dp {dp_num} vs {dp}");
    assert!(sol.history.len() < 15);
}

#[test]
fn heated_tube_meets_the_energy_and_pressure_balances() {
    let (_, sol) = solve("heated_riser_single_phase.toml");
    let tube = sol.tube("riser").unwrap();
    let m = 0.3;
    let q = 2.0e5 * PI * 0.02 * 10.0;
    // a junction inflow given by its temperature takes the enthalpy at the
    // pressure level of the case (the outlet's 16 MPa here)
    let h_in = steam::state_pt(16e6, 520.0).unwrap().h;
    // energy: exact to round-off
    assert!(
        (tube.h_end - (h_in + q / m)).abs() < 1e-6,
        "h_end {} vs {}",
        tube.h_end,
        h_in + q / m
    );
    // pressure: inlet = outlet + friction + gravity + acceleration
    let dp = sol.pressure("inlet").unwrap() - 16e6;
    assert!((dp - (tube.dp_fric + tube.dp_grav + tube.dp_acc)).abs() < 1e-4);
    // gravity of the liquid column: between rho_out g L and rho_in g L
    let rho_in = steam::state_ph(16e6, h_in).unwrap().rho;
    let rho_out = steam::state_ph(16e6, tube.h_end).unwrap().rho;
    let g = foehncraft_circuits::G * 10.0;
    assert!(tube.dp_grav < rho_in * g && tube.dp_grav > rho_out * g);
    // acceleration of the expanding liquid: G^2 (1/rho_out - 1/rho_in)
    let gf = m / (PI * 0.01 * 0.01);
    let acc = gf * gf * (1.0 / rho_out - 1.0 / rho_in);
    assert!(
        (tube.dp_acc / acc - 1.0).abs() < 2e-3,
        "acc {} vs {acc}",
        tube.dp_acc
    );
    // single phase all along: Gnielinski wall temperature T + q/h, below T_sat
    let t_sat = steam::t_sat(16e6);
    for s in &tube.segments {
        assert_eq!(s.regime, "liquid");
        assert!(s.t_wall > s.t && s.t_wall < t_sat);
    }
}

#[test]
fn drum_loop_closes_its_mass_energy_and_pressure_balances() {
    let (_, sol) = solve("drum_natural_circulation.toml");
    let d = sol.drum("drum").unwrap();
    let down = sol.flow("downcomer").unwrap();
    let risers = sol.flow("risers").unwrap() * 30.0;
    // mass: what goes down comes up
    assert!(
        (down - risers).abs() < 1e-8 * down,
        "down {down} risers {risers}"
    );
    assert!((d.riser_flow - down).abs() < 1e-8 * down);
    // energy: all the heat makes steam out of feedwater
    let q: f64 = 30.0 * 1.5e5 * PI * 0.05 * 19.0;
    let sat = steam::saturation_p(7e6).unwrap();
    let h_feed = steam::state_pt(7e6, 490.0).unwrap().h;
    let steam_exact = q / (sat.vapour.h - h_feed);
    assert!(
        (d.steam / steam_exact - 1.0).abs() < 1e-6,
        "steam {} vs {steam_exact}",
        d.steam
    );
    // momentum around the loop: the drops of the three tubes add to zero
    let total: f64 = ["downcomer", "risers", "relief"]
        .iter()
        .map(|t| {
            let m = sol.tube(t).unwrap();
            m.dp_fric + m.dp_grav + m.dp_acc
        })
        .sum();
    assert!(total.abs() < 1e-3, "loop {total} Pa");
    // a sound natural circulation: ratio of a drum boiler at 7 MPa
    assert!(
        d.circulation_ratio > 3.0 && d.circulation_ratio < 40.0,
        "CR {}",
        d.circulation_ratio
    );
}
