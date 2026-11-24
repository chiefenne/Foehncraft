//! Stage 7, PR 7a1: the energy equation of the fluid (D150, D151),
//! docs/requirements/stage7-energy-buoyancy-cht.md, checks V7.1 and V7.2.
//!
//! The case files of the gallery (`conduction_mms.toml`,
//! `couette_viscous_heating.toml` on three meshes) also run in CI; the
//! tests here use two levels.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_energy_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

/// The code and message a case text is refused with.
fn refused(text: &str) -> (String, String) {
    let e = match FlowCase::from_toml(text, &cases()) {
        Err(e) => e,
        Ok(c) => match run(&c, "refused", &out("refused"), |_| {}) {
            Err(e) => e,
            Ok(_) => panic!("the case was accepted"),
        },
    };
    (diag::of_error(&e).code.to_string(), format!("{e:#}"))
}

/// V7.1: conduction with k(T) refreshed from T; the manufactured
/// solution converges at second order and the heat balance closes.
#[test]
fn conduction_with_temperature_dependent_k_is_second_order() {
    let t = text("conduction_mms.toml")
        .replace("levels = [0, 1, 2]", "levels = [0, 1]")
        .replace("cells = [16, 16, 1]", "cells = [8, 8, 1]");
    let r = run_case(&parse(&t), "mms");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let o = &r.scalar_orders[0];
    assert_eq!(o.0, "temperature");
    assert!((1.8..=2.3).contains(&o.1[0]), "{o:?}");
    for s in &r.solves {
        let e = s.scalars[0].energy.as_ref().unwrap();
        // no flow: the conduction out of the walls takes the source
        assert!(e.convection.abs() < 1e-9, "{e:?}");
        assert!(e.relative < 5e-3, "{e:?}");
    }
}

/// V7.2: plane Couette flow heated by its own friction (Br = 1): the
/// parabolic temperature at second order, each wall takes half of the
/// friction power and the balance closes.
#[test]
fn couette_viscous_heating_is_second_order() {
    let t = text("couette_viscous_heating.toml").replace("levels = [0, 1, 2]", "levels = [0, 1]");
    let r = run_case(&parse(&t), "couette");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let o = &r.scalar_orders[0].1;
    assert!((o[0] - 2.0).abs() < 0.05, "{o:?}");
    let s = r.solves.last().unwrap();
    let sc = &s.scalars[0];
    let b = |n: &str| sc.boundaries.iter().find(|b| b.boundary == n).unwrap();
    // μ U²/h = 0.1 W/m² of friction, half into each wall: −0.05 W/m²
    for w in ["ymin", "ymax"] {
        assert!((b(w).flux + 0.05).abs() < 1e-6, "{w}: {}", b(w).flux);
    }
    let e = sc.energy.as_ref().unwrap();
    // Φ V = 0.1 W/m³ · 0.2 m³
    assert!((e.source - 0.02).abs() < 1e-9, "{e:?}");
    assert!(e.relative < 1e-6, "{e:?}");
}

/// A Couette flow of air between a wall at 300 K and a moving one at
/// 600 K: with μ(T) of Sutherland's law the shear is still uniform, but
/// the velocity profile bends, so it is no longer the linear one of the
/// constant viscosity; with `temperature_dependent = false` it is.
#[test]
fn viscosity_follows_the_temperature() {
    let base = text("couette_viscous_heating.toml")
        .replace("levels = [0, 1, 2]", "levels = [0]")
        .replace("scalar_order = [1.8, 2.2]\n", "")
        .replace(
            "density = 1.0\nviscosity = 0.1          # Re = ρ U h / μ = 10\nspecific_heat = 1000.0\nconductivity = 0.05",
            "material = \"air\"\ndensity = 0.001     # Re = 55 with the μ of air",
        )
        .replace(
            "mass_flow = 0.05         # ρ U h/2 times the depth 0.1 m: no pressure gradient",
            "pressure_gradient = 0.0  # pure shear: τ uniform, du/dy = τ/μ(T)",
        )
        .replace("tolerance = 1e-10\nmax_iterations = 2000", "tolerance = 1e-9\nmax_iterations = 3000")
        .replace("viscous_heating = true\n", "")
        .replace(
            "exact = \"T0 + mu_f * U^2 / (2 * k_f) * y * (1 - y)\"\n",
            "",
        )
        .replace(
            "[boundary.ymax.energy]\ntemperature = 300.0",
            "[boundary.ymax.energy]\ntemperature = 600.0",
        )
        + "\n[analytic]\ntype = \"couette\"\nwall_velocity = 1.0\n";
    let err = |dependent: bool| {
        let t = base.replace(
            "initial = 300.0\n",
            &format!("initial = 300.0\ntemperature_dependent = {dependent}\n"),
        );
        let r = run_case(&parse(&t), &format!("mu_{dependent}"));
        let s = &r.solves[0];
        assert!(s.converged, "dependent {dependent}");
        let t = &s.scalars[0];
        assert!(t.max > 550.0 && t.min < 350.0, "{} {}", t.min, t.max);
        s.errors.as_ref().unwrap().l2
    };
    let (fixed, varying) = (err(false), err(true));
    assert!(
        fixed < 1e-6,
        "constant μ: the linear profile, error {fixed}"
    );
    assert!(varying > 1e-2, "μ(T): a bent profile, error {varying}");
}

/// Convection to an ambient: a wall with h (T_amb − T_w) next to a wall
/// held at T₁ in a fluid at rest gives the series resistance of the slab
/// and the film, q = (T_amb − T₁)/(L/k + 1/h), exactly (a linear profile).
#[test]
fn convection_condition_is_the_series_resistance() {
    let t = text("conduction_mms.toml")
        .replace("levels = [0, 1, 2]", "levels = [0]")
        .replace("scalar_order = [1.8, 2.2]\n", "")
        .replace("perturb = 0.15", "perturb = 0.0")
        .replace("conductivity = \"1 + 0.01 * (T - 300)\"", "conductivity = 2.0")
        .replace(
            "source = \"2 * pi^2 * A * sin(pi * x) * sin(pi * y) * (1 + c * A * sin(pi * x) * sin(pi * y)) - c * A^2 * pi^2 * (cos(pi * x)^2 * sin(pi * y)^2 + sin(pi * x)^2 * cos(pi * y)^2)\"\n",
            "",
        )
        .replace("exact = \"T0 + A * sin(pi * x) * sin(pi * y)\"\n", "")
        .replace(
            "[boundary.ymax.energy]\ntemperature = 300.0",
            "[boundary.ymax.energy]\nheat_transfer_coefficient = 4.0\nambient_temperature = 400.0",
        )
        .replace("[boundary.xmin.energy]\ntemperature = 300.0\n", "")
        .replace("[boundary.xmax.energy]\ntemperature = 300.0\n", "");
    let r = run_case(&parse(&t), "convection");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let sc = &r.solves[0].scalars[0];
    let b = |n: &str| sc.boundaries.iter().find(|b| b.boundary == n).unwrap();
    // L = 1 m, k = 2, h = 4: q = 100/(0.5 + 0.25) = 133.33 W/m²
    let q = 100.0 / (1.0 / 2.0 + 1.0 / 4.0);
    assert!((b("ymax").flux - q).abs() < 1e-6 * q, "{}", b("ymax").flux);
    assert!((b("ymin").flux + q).abs() < 1e-6 * q, "{}", b("ymin").flux);
    // the wall temperature T_amb − q/h
    assert!(
        (b("ymax").mean - (400.0 - q / 4.0)).abs() < 1e-6,
        "{}",
        b("ymax").mean
    );
}

/// A heat flow in W spread over the boundary, and the heat flux of a wall
/// reported as given.
#[test]
fn heat_flow_is_spread_over_the_boundary() {
    let t = text("conduction_mms.toml")
        .replace("levels = [0, 1, 2]", "levels = [0]")
        .replace("scalar_order = [1.8, 2.2]\n", "")
        .replace("cells = [16, 16, 1]", "cells = [8, 8, 1]")
        .replace("perturb = 0.15", "perturb = 0.0")
        .replace("conductivity = \"1 + 0.01 * (T - 300)\"", "conductivity = 2.0")
        .replace(
            "source = \"2 * pi^2 * A * sin(pi * x) * sin(pi * y) * (1 + c * A * sin(pi * x) * sin(pi * y)) - c * A^2 * pi^2 * (cos(pi * x)^2 * sin(pi * y)^2 + sin(pi * x)^2 * cos(pi * y)^2)\"\n",
            "",
        )
        .replace("exact = \"T0 + A * sin(pi * x) * sin(pi * y)\"\n", "")
        .replace(
            "[boundary.ymin.energy]\ntemperature = 300.0",
            "[boundary.ymin.energy]\nheat_flow = 0.5",
        );
    let r = run_case(&parse(&t), "heat_flow");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let sc = &r.solves[0].scalars[0];
    let b = sc.boundaries.iter().find(|b| b.boundary == "ymin").unwrap();
    // 0.5 W over 1 m × 0.05 m
    assert!((b.flux - 10.0).abs() < 1e-9, "{}", b.flux);
    assert!((b.heat_flow.unwrap() - 0.5).abs() < 1e-9);
    let e = sc.energy.as_ref().unwrap();
    assert!(e.relative < 1e-6, "{e:?}");
}

/// The diagnostics of the energy keys and the thermal set-up (D156).
#[test]
fn energy_keys_are_checked() {
    let base = text("couette_viscous_heating.toml");
    // an inlet needs a temperature
    let inlet = base
        .replace(
            "[boundary.xmin]\ntype = \"periodic\"\npartner = \"xmax\"\nmass_flow = 0.05         # ρ U h/2 times the depth 0.1 m: no pressure gradient\n[boundary.xmax]\ntype = \"periodic\"\n",
            "[boundary.xmin]\ntype = \"velocity-inlet\"\nvelocity = [1.0, 0.0, 0.0]\n[boundary.xmax]\ntype = \"static-pressure-outlet\"\npressure = 0.0\n",
        );
    assert!(inlet.contains("velocity-inlet"));
    let (code, msg) = refused(&inlet);
    assert_eq!(code, "THM011", "{msg}");
    // two conditions on one wall
    let two = base.replace(
        "[boundary.ymin.energy]\ntemperature = 300.0",
        "[boundary.ymin.energy]\ntemperature = 300.0\nheat_flux = 1.0",
    );
    let (code, msg) = refused(&two);
    assert_eq!(code, "THM010", "{msg}");
    // an energy table without [energy]
    let (a, b) = (
        base.find("[energy]").unwrap(),
        base.find("[boundary.xmin]").unwrap(),
    );
    let lone = format!("{}{}", &base[..a], &base[b..]);
    let (code, msg) = refused(&lone);
    assert_eq!(code, "THM009", "{msg}");
    // a temperature on a symmetry plane
    let sym = base.replace(
        "[boundary.zmin]\ntype = \"symmetry\"\n",
        "[boundary.zmin]\ntype = \"symmetry\"\n[boundary.zmin.energy]\ntemperature = 300.0\n",
    );
    assert_ne!(sym, base);
    let (code, msg) = refused(&sym);
    assert_eq!(code, "THM003", "{msg}");
    // a temperature at or below 0 K, a heat transfer coefficient ≤ 0
    let cold = base.replace(
        "[boundary.ymin.energy]\ntemperature = 300.0",
        "[boundary.ymin.energy]\ntemperature = -5.0",
    );
    let (code, msg) = refused(&cold);
    assert_eq!(code, "THM004", "{msg}");
    let h0 = base.replace(
        "[boundary.ymin.energy]\ntemperature = 300.0",
        "[boundary.ymin.energy]\nheat_transfer_coefficient = 0.0\nambient_temperature = 300.0",
    );
    let (code, msg) = refused(&h0);
    assert_eq!(code, "THM004", "{msg}");
    // both walls adiabatic, steady: viscous heating with no way out
    let shut = base
        .replace("[boundary.ymin.energy]\ntemperature = 300.0\n", "")
        .replace("[boundary.ymax.energy]\ntemperature = 300.0\n", "");
    let (code, msg) = refused(&shut);
    assert_eq!(code, "THM002", "{msg}");
    // ... and without it, no temperature level
    let (code, msg) = refused(&shut.replace("viscous_heating = true\n", ""));
    assert_eq!(code, "THM001", "{msg}");
    // no specific heat
    let no_cp = base.replace("specific_heat = 1000.0\n", "");
    let (code, msg) = refused(&no_cp);
    assert_eq!(code, "MAT007", "{msg}");
    // a passive scalar of the same name
    let twice = base + "\n[scalars.temperature]\ndiffusivity = 1e-3\n";
    let (code, msg) = refused(&twice);
    assert_eq!(code, "ENE006", "{msg}");
}

/// A channel heated from below, linear upwind: the heat balance of the
/// report counts the enthalpy the scheme carries out of the outlet (its
/// reconstructed face value), so it closes as with upwind.
#[test]
fn heat_balance_closes_with_inflow_and_outflow() {
    let t = r#"
case_format = "0.6"
[mesh.generate]
type = "hex"
dimension = 2
upper = [2.0, 0.5, 0.05]
cells = [40, 10, 1]
[fluid]
density = 1.0
viscosity = 0.01
specific_heat = 1000.0
conductivity = 14.0
[energy]
initial = 300.0
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmin.energy]
temperature = 300.0
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "wall"
[boundary.ymin.energy]
heat_flux = 1000.0
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
tolerance = 1e-9
max_iterations = 1000
[solver.energy]
relaxation = 1.0
tolerance = 1e-12
max_iter = 300
[output]
wall_csv = false
vtkhdf = false
"#;
    for scheme in ["linear-upwind", "upwind"] {
        let case = t.replace(
            "max_iter = 300\n",
            &format!("max_iter = 300\nscheme = \"{scheme}\"\n"),
        );
        let r = run_case(&parse(&case), &format!("balance_{scheme}"));
        let s = &r.solves[0];
        assert!(s.converged, "{scheme}");
        let e = s.scalars[0].energy.as_ref().unwrap();
        // 1000 W/m² over 2 m × 0.05 m
        assert!((e.conduction - 100.0).abs() < 2.0, "{scheme}: {e:?}");
        assert!(e.relative < 1e-6, "{scheme}: {e:?}");
    }
}

/// D119, D161: the thermal account of the ledger is the heat balance of
/// the energy equation as data (per boundary and source rows), closes to
/// the same number, and is part of `result.json` next to mass and kinetic
/// energy.
#[test]
fn the_thermal_account_is_the_heat_balance_as_data() {
    let t = text("couette_viscous_heating.toml").replace("levels = [0, 1, 2]", "levels = [0]");
    let r = run_case(&parse(&t), "ledger");
    let s = r.solves.last().unwrap();
    let e = s.scalars[0].energy.as_ref().unwrap();
    let bal = &s.balance.as_ref().unwrap().last;
    let q: Vec<&str> = bal.ledger.iter().map(|a| a.quantity.as_str()).collect();
    assert_eq!(q, ["mass", "kinetic energy", "thermal energy"]);
    let th = &bal.ledger[2];
    // the sign of the closure is storage − contributions
    assert!((th.imbalance + e.imbalance).abs() < 1e-12, "{th:?} {e:?}");
    let src = th.terms.iter().find(|t| t.name == "heat sources").unwrap();
    assert_eq!(src.value, e.source);
    for w in ["ymin", "ymax"] {
        let t = th
            .terms
            .iter()
            .find(|t| t.name == "conduction" && t.boundary.as_deref() == Some(w));
        let b = s.scalars[0]
            .boundaries
            .iter()
            .find(|b| b.boundary == w)
            .unwrap();
        assert_eq!(t.unwrap().value, b.heat_flow.unwrap());
    }
    // the run record carries it
    let j = serde_json::to_value(&r.solves).unwrap();
    assert_eq!(
        j[0]["balance"]["last"]["ledger"][2]["quantity"],
        "thermal energy"
    );
}
