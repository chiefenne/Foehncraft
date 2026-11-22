//! Stage 7, PR 7c: the low-Mach variable density ρ(T, p₀) (D151), checks
//! V7c.1 to V7c.3 of the stage 7c design: the exact balances of a heated
//! duct, the thermodynamic pressure of a closed volume heated uniformly,
//! and the mass and energy of a closed box heated on one side.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

/// R/M of the ideal gas of the cases (J/(kg K)): the molar mass of air.
const R_GAS: f64 = 8.314462618 / 0.0289647;
const CP: f64 = 1005.0;
const P_OP: f64 = 101_325.0;

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_vd_{name}_{}", std::process::id()))
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
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

/// V7c.1: ideal gas heated by 100 K in a duct with slip walls. The mass
/// flux is the same at both ends, the heat ṁ c_p ΔT equals the source,
/// and the inlet pressure is the momentum flux difference G²(1/ρ_out −
/// 1/ρ_in) of the expanding gas.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn heated_duct_meets_the_one_dimensional_balances() {
    let rep = run_case(&parse(&text("vd_heated_duct.toml")), "duct");
    let s = &rep.solves[0];
    assert!(s.converged);
    let b = |n: &str| s.boundaries.iter().find(|b| b.name == n).unwrap();
    let (inl, outl) = (b("xmin"), b("xmax"));
    let rho_in = P_OP / (R_GAS * 300.0);
    let g = rho_in * 1.0; // ρ u at the inlet (kg/(m² s))
    let area = 0.1 * 0.01;
    assert!(
        ((-inl.mass_flow) - g * area).abs() < 1e-12 * g * area,
        "inlet mass flow {} against {}",
        -inl.mass_flow,
        g * area
    );
    assert!(
        (outl.mass_flow + inl.mass_flow).abs() < 1e-12 * g * area,
        "mass flows {} and {}",
        inl.mass_flow,
        outl.mass_flow
    );
    // the outlet temperature: Q = ṁ c_p ΔT
    let t = s.scalars.iter().find(|x| x.name == "temperature").unwrap();
    let t_out = t
        .boundaries
        .iter()
        .find(|x| x.boundary == "xmax")
        .and_then(|x| x.mixed_mean)
        .unwrap();
    let q = 5.9131e5 * 0.2 * area;
    let dt_exact = q / (g * area * CP);
    assert!(
        ((t_out - 300.0) - dt_exact).abs() < 1e-4 * dt_exact,
        "outlet ΔT {} against {dt_exact}",
        t_out - 300.0
    );
    let e = t.energy.as_ref().unwrap();
    assert!(e.relative.abs() < 1e-8, "heat imbalance {}", e.relative);
    // momentum: p_in − p_out = G² (1/ρ_out − 1/ρ_in)
    let rho_out = P_OP / (R_GAS * t_out);
    let dp = g * g * (1.0 / rho_out - 1.0 / rho_in);
    assert!(
        (inl.mean_pressure - outl.mean_pressure - dp).abs() < 2e-3 * dp,
        "pressure difference {} against {dp}",
        inl.mean_pressure - outl.mean_pressure
    );
    let d = s.density.as_ref().unwrap();
    assert_eq!(d.p0, P_OP);
    assert!(!d.closed_volume);
    assert!(d.mach_max.unwrap() < 0.01);
    assert!((d.rho_min - rho_out).abs() < 1e-3 * rho_out);
    // the kinetic-energy ledger books the pressure dilatation as physics
    let ke = &s.balance.as_ref().unwrap().last.kinetic_energy;
    assert!(ke.relative > -1e-2, "numerical dissipation {}", ke.relative);
}

/// V7c.2: a closed box of ideal gas at rest heated uniformly warms with
/// c_v: p₀(t) = p₀(0) T(t)/T₀, T(t) = T₀ + q''' t/(ρ₀ c_v), linear in t, so
/// exact at any Δt; the mass stays to round-off.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn closed_box_heated_uniformly_warms_with_c_v() {
    let rep = run_case(&parse(&text("vd_closed_box_heating.toml")), "box");
    let s = &rep.solves[0];
    let rho0 = P_OP / (R_GAS * 300.0);
    let cv = CP - R_GAS;
    let t_end = 300.0 + 1.0e4 * 5.0 / (rho0 * cv);
    let p0 = P_OP * t_end / 300.0;
    let d = s.density.as_ref().unwrap();
    assert!(d.closed_volume);
    assert!((d.p0 - p0).abs() < 1e-10 * p0, "p₀ {} against {p0}", d.p0);
    assert!(
        (d.mass - d.mass_start).abs() < 1e-13 * d.mass_start,
        "mass {} against {}",
        d.mass,
        d.mass_start
    );
    let t = s.scalars.iter().find(|x| x.name == "temperature").unwrap();
    assert!(
        (t.max - t_end).abs() < 1e-8 * t_end,
        "T {} against {t_end}",
        t.max
    );
    assert!((t.min - t_end).abs() < 1e-8 * t_end);
    // without the pressure work the gas warms with c_p
    let none = text("vd_closed_box_heating.toml").replace(
        "thermodynamic_pressure = \"closed-volume\"",
        "thermodynamic_pressure = \"closed-volume\"\npressure_work = \"none\"",
    );
    let r2 = run_case(&parse(&none), "box_none");
    let t2 = r2.solves[0]
        .scalars
        .iter()
        .find(|x| x.name == "temperature")
        .unwrap()
        .max;
    let t_cp = 300.0 + 1.0e4 * 5.0 / (rho0 * CP);
    assert!((t2 - t_cp).abs() < 1e-8 * t_cp, "T {t2} against {t_cp}");
}

/// V7c.3: a closed box heated through one wall: the expansion drives a
/// flow, the mass stays to round-off, the mass ledger of the steps closes
/// with the storage, and p₀ is that of the mass at the final temperature.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn closed_box_heated_on_one_side_keeps_its_mass() {
    let rep = run_case(&parse(&text("vd_closed_box_side_heated.toml")), "side");
    let s = &rep.solves[0];
    let d = s.density.as_ref().unwrap();
    assert!(
        (d.mass - d.mass_start).abs() < 1e-13 * d.mass_start,
        "mass {} against {}",
        d.mass,
        d.mass_start
    );
    assert!(d.p0 > P_OP, "p₀ {} did not rise", d.p0);
    let bal = s.balance.as_ref().unwrap();
    assert!(
        bal.largest_mass_relative < 1e-8,
        "mass imbalance of a step {}",
        bal.largest_mass_relative
    );
    let e = s
        .scalars
        .iter()
        .find(|x| x.name == "temperature")
        .unwrap()
        .energy
        .as_ref()
        .unwrap();
    assert!(e.relative.abs() < 1e-6, "heat imbalance {}", e.relative);
}

/// The settings that have no low-Mach solution or are not built yet are
/// refused with their codes.
#[test]
fn variable_density_settings_are_checked() {
    let duct = text("vd_heated_duct.toml");
    // a constant density form
    let (c, m) = refused(&duct.replace(
        "density = { form = \"ideal-gas\", source = \"ideal gas, rho = p M / (R T)\" }",
        "density = 1.2",
    ));
    assert_eq!(c, "ENE009", "{m}");
    // a closed volume in an open domain
    let (c, m) = refused(&duct.replace(
        "density = \"variable\"",
        "density = \"variable\"\nthermodynamic_pressure = \"closed-volume\"",
    ));
    assert_eq!(c, "ENE009", "{m}");
    // a closed domain with a constant p₀
    let closed = text("vd_closed_box_heating.toml");
    let (c, m) = refused(&closed.replace(
        "thermodynamic_pressure = \"closed-volume\"",
        "thermodynamic_pressure = \"constant\"",
    ));
    assert_eq!(c, "BUO005", "{m}");
    // closed-volume without the variable density
    let (c, m) = refused(&closed.replace("density = \"variable\"\n", ""));
    assert_eq!(c, "ENE009", "{m}");
    // with the Boussinesq model
    let (c, m) = refused(&duct.replace(
        "[energy]",
        "[gravity]\nvector = [0.0, -9.81, 0.0]\n\n[buoyancy]\nmodel = \"boussinesq\"\n\n[energy]",
    ));
    assert_eq!(c, "ENE009", "{m}");
    // not built yet: passive scalars
    let (c, m) = refused(&duct.replace(
        "[boundary.xmin]",
        "[scalars.c]\ndiffusivity = 1e-3\ninitial = 0.0\n\n[boundary.xmin]",
    ));
    assert_eq!(c, "ENE006", "{m}");
}

/// V7c.4: the manufactured solution (ρ = 1/T, a divergence-free mass flux,
/// the full stress with ⅔ of the trace, a closed volume) at second order
/// for velocity and temperature on nested hexahedra and polygons; SIMPLEC
/// gives the solution of the coupled algorithm.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn manufactured_solution_converges_at_second_order() {
    let t = text("vd_mms_steady.toml");
    for (kind, band) in [("hex", [1.8, 2.2]), ("poly", [1.8, 2.5])] {
        let rep = run_case(
            &parse(&t.replace("type = \"hex\"", &format!("type = \"{kind}\""))),
            &format!("mms_{kind}"),
        );
        let u = *rep.orders.last().unwrap();
        let tt = *rep.scalar_orders[0].1.last().unwrap();
        for (what, o) in [("velocity", u), ("temperature", tt)] {
            assert!(
                o >= band[0] && o <= band[1],
                "{kind}: {what} order {o} outside {band:?} ({:?}, {:?})",
                rep.orders,
                rep.scalar_orders
            );
        }
    }
    // SIMPLEC on the coarsest mesh: the same discrete solution
    let one = |alg: &str| {
        let c = t
            .replace("levels = [0, 1, 2]", "levels = [0]")
            .replace("algorithm = \"coupled\"", &format!("algorithm = \"{alg}\""));
        let rep = run_case(&parse(&c), &format!("mms_{alg}"));
        rep.solves[0].errors.as_ref().unwrap().l2
    };
    let (a, b) = (one("coupled"), one("simplec"));
    assert!((a - b).abs() < 1e-6 * a, "coupled {a}, SIMPLEC {b}");
}
