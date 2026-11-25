//! Stage SF7a: IAPWS-IF97 water and steam in a flow. Superheated steam
//! heated in a duct with slip walls (`steam_heated_duct.toml`) meets the
//! exact balances of a one-dimensional low-Mach flow with the IF97 states:
//! the same mass flux at both ends, the outlet temperature of the IF97
//! enthalpy balance, and the momentum flux difference of the expanding
//! steam.

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::steam;
use std::path::{Path, PathBuf};

const P: f64 = 4e6;

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn steam_duct_meets_the_if97_balances() {
    let text = std::fs::read_to_string(cases().join("steam_heated_duct.toml"))
        .unwrap()
        .replace("\r\n", "\n");
    let case = FlowCase::from_toml(&text, &cases()).unwrap();
    let dir = std::env::temp_dir().join(format!("foehncraft_steam_duct_{}", std::process::id()));
    let rep = run(&case, "steam_duct", &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let s = &rep.solves[0];
    assert!(s.converged);
    let b = |n: &str| s.boundaries.iter().find(|b| b.name == n).unwrap();
    let (inl, outl) = (b("xmin"), b("xmax"));
    let area = 0.1 * 0.01;
    let s_in = steam::props_pt(P, 600.0).unwrap();
    let g = s_in.rho * 2.0;
    // the density curve is sampled every 0.25 K at most: 1e-7 relative
    assert!(
        ((-inl.mass_flow) - g * area).abs() < 1e-7 * g * area,
        "inlet mass flow {} against {}",
        -inl.mass_flow,
        g * area
    );
    assert!((outl.mass_flow + inl.mass_flow).abs() < 1e-12 * g * area);
    // the outlet temperature of the IF97 enthalpy balance: 750 K
    let q = 5.6844892012526296e7 * 0.2;
    let h_out = s_in.h + q / g;
    let s_out = steam::state_ph(P, h_out).unwrap();
    assert!((s_out.t - 750.0).abs() < 1e-6, "IF97 outlet {}", s_out.t);
    let t = s.scalars.iter().find(|x| x.name == "temperature").unwrap();
    let t_out = t
        .boundaries
        .iter()
        .find(|x| x.boundary == "xmax")
        .and_then(|x| x.mixed_mean)
        .unwrap();
    // the energy equation integrates the sampled IF97 c_p: a few mK
    assert!(
        (t_out - s_out.t).abs() < 0.01,
        "outlet T {t_out} against {}",
        s_out.t
    );
    let e = t.energy.as_ref().unwrap();
    assert!(e.relative.abs() < 1e-8, "heat imbalance {}", e.relative);
    // momentum: p_in − p_out = G² (1/ρ_out − 1/ρ_in)
    let dp = g * g * (1.0 / s_out.rho - 1.0 / s_in.rho);
    let dp_run = inl.mean_pressure - outl.mean_pressure;
    assert!((dp_run - dp).abs() < 1e-4 * dp, "Δp {dp_run} against {dp}");
}
