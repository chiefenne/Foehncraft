//! The water and steam circuits behind the walls (steam furnace stage SF7c,
//! D350): a slab cooled by a water tube against an independent march of
//! the same heat balance, a boiling tube behind the slab (the nonlinear
//! wall of flow boiling), the exchange pair of the ledger, and the checks
//! of `[circuits]`.

use foehncraft_flow::case::{FlowCase, run};
use std::f64::consts::PI;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "foehncraft_water_walls_{name}_{}",
        std::process::id()
    ))
}

fn load(file: &str, sets: &[&str]) -> FlowCase {
    let sets: Vec<String> = sets.iter().map(|s| s.to_string()).collect();
    FlowCase::from_file_with(&cases().join(file), &sets).unwrap()
}

/// The rows of `circuit/walls.csv`: x, q into the water (W/m²), T_wall of
/// the 3D face (K), U (W/(m² K)), T_water (K).
fn walls(dir: &Path) -> Vec<[f64; 5]> {
    let text = std::fs::read_to_string(dir.join("circuit/walls.csv")).unwrap();
    text.lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split(',').collect();
            let v = |i: usize| c[i].parse::<f64>().unwrap();
            [v(3), v(9), v(10), v(7), v(8)]
        })
        .collect()
}

/// The march of the heat balance of the tube behind the slab, independent
/// of the coupled run: m dh/dx = w (T_top − T(h)) / (t/k + φ/h_in(h)),
/// with the IF97 states, Gnielinski for h_in and φ = w/(π d), by RK4 on
/// 4000 steps. Returns the outlet temperature (K).
fn slab_reference() -> f64 {
    use foehncraft_circuits::correlations::nusselt_tube;
    use foehncraft_circuits::fluid::Phase1;
    use foehncraft_kernels::steam;
    let (p, m, d, w, t, k, t_top, len) = (10.02e6, 0.05, 0.02, 0.05, 0.05, 50.0, 700.0, 2.0);
    let g = m / (PI * d * d / 4.0);
    let phi = w / (PI * d);
    let rate = |h: f64| -> f64 {
        let s = steam::state_ph(p, h).unwrap();
        let ph = Phase1::of(&s);
        let hin = nusselt_tube(g * d / ph.mu, ph.pr()) * ph.k / d;
        w * (t_top - s.t) / (t / k + phi / hin) / m
    };
    let mut h = steam::state_pt(p, 450.0).unwrap().h;
    let n = 4000;
    let dx = len / n as f64;
    for _ in 0..n {
        let k1 = rate(h);
        let k2 = rate(h + 0.5 * dx * k1);
        let k3 = rate(h + 0.5 * dx * k2);
        let k4 = rate(h + dx * k3);
        h += dx / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    steam::state_ph(p, h).unwrap().t
}

/// The slab of `circuit_wall_slab.toml`: the coupling converges below its
/// tolerance 1e-6 in a few circuit solves; every face sends the heat flux of
/// the series resistance of the slab and the water side, (T_top − T_water) /
/// (t/k + 1/U), with the face's own U and T_water (the slab is a straight
/// conduction path: only the axial conduction of the slab, under 10⁻³,
/// separates them away from the adiabatic ends); the outlet temperature is that of the independent march
/// to 0.01 K (40 segments; measured 1.2·10⁻⁴ K); the
/// circuit closes its energy to round-off and the exchange pair of the
/// ledger is the coupling residual.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn slab_cooled_by_a_tube_meets_the_heat_balance_march() {
    let case = load("circuit_wall_slab.toml", &["verification.levels=[0]"]);
    let dir = out("slab");
    let rep = run(&case, "slab", &dir, |_| {}).unwrap();
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let s = &rep.solves[0];
    let c = s.circuits.as_ref().expect("the circuit summary");
    assert!(c.converged && c.residual < 1e-6, "{c:?}");
    assert!(c.solves <= 12 && c.failed == 0, "{c:?}");
    assert!(c.circuit_closure.abs() < 1e-10, "{c:?}");
    let t_out = c.tubes[0].outlet_temperature;
    let t_ref = slab_reference();
    eprintln!(
        "slab: outlet {t_out} K, march {t_ref} K, {} solves",
        c.solves
    );
    assert!(
        (t_out - t_ref).abs() < 0.01,
        "outlet {t_out} K, reference {t_ref} K"
    );
    assert!(t_out - 450.0 > 60.0, "the water warms: {t_out}");
    // away from the adiabatic ends of the slab (two thicknesses), where
    // the slab conducts along x as well
    for [x, q, _tw, u, t_water] in walls(&dir)
        .into_iter()
        .filter(|r| (0.1..1.9).contains(&r[0]))
    {
        let series = (700.0 - t_water) / (0.05 / 50.0 + 1.0 / u);
        assert!(
            (q - series).abs() < 3e-4 * series,
            "x = {x}: q {q}, series {series}"
        );
    }
    // the ledger: the circuit's account closes, the pair is the lag
    let bal = &s.balance.as_ref().unwrap().last;
    let acc = bal
        .ledger
        .iter()
        .find(|a| a.region == "circuit")
        .expect("the circuit account");
    assert!(acc.relative.abs() < 1e-10, "{acc:?}");
    let pair = bal
        .pairs
        .iter()
        .find(|p| p.pair == foehncraft_flow::water_walls::PAIR_WALLS)
        .expect("the pair");
    assert!(
        pair.residual.abs() <= 2e-6 * c.wall_heat,
        "{pair:?}, wall heat {}",
        c.wall_heat
    );
    assert!(dir.join("circuit/result.json").exists());
    assert!(dir.join("circuit/profiles.csv").exists());
    let _ = std::fs::remove_dir_all(&dir);
}

/// The same slab on a boiling tube: water at 7 MPa enters 5 K below
/// saturation and leaves as wet steam. The wall of flow boiling (Chen)
/// makes the water side nonlinear in the heat flux; the coupling still
/// converges, the water stays at the saturation temperature over most of
/// the tube, and every face still meets the series resistance.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn slab_on_a_boiling_tube_converges() {
    let case = load(
        "circuit_wall_slab.toml",
        &[
            "verification.levels=[0]",
            r#"circuits.file="../circuits/wall_slab_boiling.toml""#,
        ],
    );
    let dir = out("boiling");
    let rep = run(&case, "boiling", &dir, |_| {}).unwrap();
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let c = rep.solves[0].circuits.as_ref().unwrap();
    assert!(c.converged && c.failed == 0, "{c:?}");
    assert!(c.circuit_closure.abs() < 1e-9, "{c:?}");
    let t = &c.tubes[0];
    assert!(t.outlet_quality > 0.05, "{t:?}");
    let t_sat = foehncraft_kernels::steam::saturation_p(7.0e6).unwrap().t;
    let rows = walls(&dir);
    let boiling = rows.iter().filter(|r| (r[4] - t_sat).abs() < 1.0).count();
    assert!(boiling > rows.len() / 2, "{boiling} of {}", rows.len());
    for [x, q, _tw, u, t_water] in rows.into_iter().filter(|r| (0.1..1.9).contains(&r[0])) {
        let series = (700.0 - t_water) / (0.05 / 50.0 + 1.0 / u);
        assert!(
            (q - series).abs() < 1e-3 * series,
            "x = {x}: q {q}, series {series}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// The checks of `[circuits]` (CIR005).
#[test]
fn circuits_table_is_checked() {
    let fails = |edit: &dyn Fn(&mut FlowCase), needle: &str| {
        let mut c = load("circuit_wall_slab.toml", &[]);
        edit(&mut c);
        let e = c.validate().expect_err("must fail");
        let text = format!("{e:#}");
        assert!(text.contains(needle), "{needle}: {text}");
        let d = foehncraft_flow::diag::of_error(&e);
        assert_eq!(d.code, "CIR005", "{text}");
    };
    let wall = |c: &mut FlowCase| c.circuits.as_mut().unwrap().walls[0].clone();
    fails(
        &|c| c.circuits.as_mut().unwrap().walls[0].tube = "riser".into(),
        "no tube `riser`",
    );
    fails(
        &|c| c.circuits.as_mut().unwrap().walls[0].boundary = "zmin".into(),
        "is not a wall",
    );
    fails(
        &|c| c.circuits.as_mut().unwrap().walls[0].boundary = "none".into(),
        "is not a boundary",
    );
    fails(
        &|c| c.circuits.as_mut().unwrap().walls[0].end = [0.0, 0.0, 0.025],
        "two different",
    );
    fails(
        &|c| c.circuits.as_mut().unwrap().walls[0].fouling = -1.0,
        "must be a number",
    );
    fails(
        &|c| {
            let w = wall(c);
            c.circuits.as_mut().unwrap().walls.push(w);
        },
        "on two tubes",
    );
    fails(&|c| c.circuits.as_mut().unwrap().relaxation = 0.0, "(0, 1]");
    fails(&|c| c.circuits.as_mut().unwrap().every = 0, "≥ 1");
    fails(
        &|c| c.circuits.as_mut().unwrap().file = "nothing.toml".into(),
        "nothing.toml",
    );
    let with = |sets: &[&str], needle: &str| {
        let sets: Vec<String> = sets.iter().map(|s| s.to_string()).collect();
        let path = cases().join("circuit_wall_slab.toml");
        let e = match FlowCase::from_file_with(&path, &sets) {
            Ok(c) => c.validate().expect_err("must fail"),
            Err(e) => e,
        };
        let text = format!("{e:#}");
        assert!(text.contains(needle), "{sets:?}: {text}");
    };
    with(
        &["boundary.ymin.energy.temperature=400.0"],
        "remove the energy table",
    );
    with(
        &["time.end_time=1.0", "time.time_step=0.1"],
        "transient circuits",
    );
}

/// A tube axis that does not span the boundary is refused when the run
/// maps the faces.
#[test]
fn tube_axis_must_span_the_wall() {
    let mut case = load("circuit_wall_slab.toml", &[]);
    case.circuits.as_mut().unwrap().walls[0].end = [1.0, 0.0, 0.025];
    let dir = out("span");
    let e = run(&case, "span", &dir, |_| {}).expect_err("faces beyond the axis");
    let text = format!("{e:#}");
    assert!(text.contains("outside the tube axis"), "{text}");
    let _ = std::fs::remove_dir_all(&dir);
}
