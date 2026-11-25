//! The frozen-flow march of a transient conjugate run (D118,
//! docs/requirements/stage7-energy-buoyancy-cht.md §9f).

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_frz_{name}_{}", std::process::id()))
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

/// The gallery case on a coarser mesh to `end` s, with a flow step of 1 s;
/// `viscosity`: the viscosity that follows the temperature (else the
/// constant one of 300 K); `frozen`: with the march (else fully coupled).
fn heated_wall(end: f64, viscosity: bool, frozen: bool) -> String {
    let mut t = text("cht_heated_wall_frozen.toml")
        .replace("cells = [80, 30, 1]", "cells = [40, 15, 1]")
        .replace(
            "time_step = 0.5                      # the flow's step (Courant number about 10)",
            "time_step = 1.0",
        )
        .replace("end_time = 3000.0", &format!("end_time = {end:.1}"))
        .replace(
            "time_step = 10.0                     # s, the step of the march with the flow held",
            "time_step = 20.0",
        );
    if !viscosity {
        t = t.replace(
            "viscosity = { form = \"power-law\", value_ref = 0.01, t_ref = 300.0, exponent = -12.0 }",
            "viscosity = 0.01",
        );
    }
    if !frozen {
        t = t.replace(
            "[energy.frozen_flow]\n",
            "[energy.frozen_flow]\nenabled = false\n",
        );
    }
    t
}

/// The rows of a CSV file of the output with its header.
fn csv(dir: &Path, file: &str) -> (Vec<String>, Vec<Vec<f64>>) {
    let text = std::fs::read_to_string(dir.join(file)).unwrap();
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
                .map(|v| v.parse().unwrap_or(f64::NAN))
                .collect()
        })
        .collect();
    (head, rows)
}

fn col(head: &[String], name: &str) -> usize {
    head.iter().position(|h| h == name).unwrap()
}

/// The largest difference of probe `name` of run `a` against run `b` over
/// the times of `a` (b interpolated linearly), and the difference at the
/// end.
fn history_difference(a: &Path, b: &Path, name: &str) -> (f64, f64) {
    let (h, ra) = csv(a, "probes.csv");
    let (_, rb) = csv(b, "probes.csv");
    let k = col(&h, name);
    let at = |t: f64| -> f64 {
        let i = rb.partition_point(|r| r[1] < t).clamp(1, rb.len() - 1);
        let (r0, r1) = (&rb[i - 1], &rb[i]);
        let w = (t - r0[1]) / (r1[1] - r0[1]);
        r0[k] + w * (r1[k] - r0[k])
    };
    let mut worst = 0.0_f64;
    for r in &ra {
        worst = worst.max((r[k] - at(r[1])).abs());
    }
    let last = ra.last().unwrap();
    (worst, (last[k] - at(last[1])).abs())
}

fn frozen(r: &FlowReport) -> foehncraft_flow::frozen::FrozenSummary {
    r.solves[0]
        .transient
        .as_ref()
        .unwrap()
        .frozen_flow
        .clone()
        .unwrap()
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

const PROBES: [&str; 3] = ["heated_face_T", "interface_T", "channel_T"];

/// A flow that does not depend on the temperature: held once after it has
/// settled, never refreshed, and the march with steps twenty times the
/// flow's follows the fully coupled heat-up of the wall (the difference is
/// the time error of the long steps); the heat balance closes in every
/// step, frozen or coupled.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn frozen_flow_follows_the_fully_coupled_heat_up() {
    let (da, db) = (out("const"), out("const_ref"));
    let case = parse(&heated_wall(1200.0, false, true));
    let ra = run(&case, "frozen", &da, |_| {}).unwrap();
    run(
        &parse(&heated_wall(1200.0, false, false)),
        "coupled",
        &db,
        |_| {},
    )
    .unwrap();
    let f = frozen(&ra);
    assert!(!f.flow_follows_temperature);
    let at = f.frozen_at.expect("the flow was held");
    assert!(at < 100.0, "held at {at} s");
    assert_eq!(f.refreshes, 0);
    assert!(f.frozen_steps > 50, "{f:?}");
    assert!(f.coupled_steps < 100, "{f:?}");
    let tr = ra.solves[0].transient.as_ref().unwrap();
    assert!(
        tr.heat_relative_max.unwrap() < 1e-9,
        "{:?}",
        tr.heat_relative_max
    );
    for p in PROBES {
        let (worst, end) = history_difference(&da, &db, p);
        assert!(worst < 0.1, "{p}: {worst} K over the history");
        assert!(end < 0.01, "{p}: {end} K at the end");
    }
    let lims: Vec<String> = foehncraft_flow::limitations::of_run(&case, &ra)
        .into_iter()
        .map(|l| l.code.to_string())
        .collect();
    assert!(lims.iter().any(|c| c == diag::LIM_FROZEN), "{lims:?}");
}

/// A viscosity that follows the temperature: the flow is refreshed by
/// coupled steps whenever a fluid temperature has moved by more than
/// `refresh_change`, and the history stays within a few hundredths of a
/// kelvin of the fully coupled run.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn viscosity_that_follows_the_temperature_refreshes_the_flow() {
    let (da, db) = (out("visc"), out("visc_ref"));
    let ra = run(
        &parse(&heated_wall(1200.0, true, true)),
        "frozen",
        &da,
        |_| {},
    )
    .unwrap();
    run(
        &parse(&heated_wall(1200.0, true, false)),
        "coupled",
        &db,
        |_| {},
    )
    .unwrap();
    let f = frozen(&ra);
    assert!(f.flow_follows_temperature);
    assert!(f.refreshes >= 3, "{f:?}");
    assert!(
        f.phases
            .iter()
            .any(|p| p.kind == "coupled" && p.reason == "temperature"),
        "{:?}",
        f.phases
    );
    assert!(f.fluid_change <= 2.0 + 1.0, "{}", f.fluid_change);
    let tr = ra.solves[0].transient.as_ref().unwrap();
    assert!(
        tr.heat_relative_max.unwrap() < 1e-9,
        "{:?}",
        tr.heat_relative_max
    );
    for p in PROBES {
        let (worst, _) = history_difference(&da, &db, p);
        assert!(worst < 0.15, "{p}: {worst} K over the history");
    }
}

/// A run continued from a restart file written at `t_split` is
/// bit-identical to the run without it (probes and heat balance of the
/// last step, the phases of the march).
#[cfg(feature = "hdf5")] // restart files are HDF5
fn continues_bit_identically(label: &str, text: &str, t_split: f64) {
    let whole = out(&format!("{label}_whole"));
    let (a, b) = (out(&format!("{label}_a")), out(&format!("{label}_b")));
    let rw = run(&parse(text), label, &whole, |_| {}).unwrap();
    let mut first = parse(text);
    first.time.as_mut().unwrap().end_time = Some(t_split);
    first.run.restart_final = true;
    run(&first, label, &a, |_| {}).unwrap();
    let mut second = parse(text);
    second.run.restart = Some(a.join("restart.h5"));
    let rb = run(&second, label, &b, |_| {}).unwrap();
    for f in ["probes.csv", "heat.csv"] {
        let (hw, w) = csv(&whole, f);
        let (hb, v) = csv(&b, f);
        assert_eq!(hw, hb, "{label}: {f}");
        let (lw, lb) = (w.last().unwrap(), v.last().unwrap());
        assert!(
            lw.iter()
                .zip(lb)
                .all(|(x, y)| x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan())),
            "{label}: {f}\n{lw:?}\n{lb:?}"
        );
    }
    let (fw, fb) = (frozen(&rw), frozen(&rb));
    assert_eq!(fw.refreshes, fb.refreshes, "{label}");
    assert_eq!(fw.frozen_steps, fb.frozen_steps, "{label}");
    assert_eq!(fw.phases, fb.phases, "{label}");
    for d in [whole, a, b] {
        let _ = std::fs::remove_dir_all(&d);
    }
}

/// The state of the march goes into restart.h5 (`cht_frozen_state`): a
/// run restarted in a frozen phase, and in the coupled steps of a refresh,
/// continues bit for bit as the run without the restart.
#[cfg(feature = "hdf5")] // restart files are HDF5
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn frozen_flow_restart_continues_bit_identically() {
    let base = heated_wall(600.0, true, true);
    let probe = out("rs_probe");
    let f = frozen(&run(&parse(&base), "probe", &probe, |_| {}).unwrap());
    let _ = std::fs::remove_dir_all(&probe);
    // a time on the grid of the steps: inside the longest frozen phase,
    // and two coupled steps into a refresh
    let long = f
        .phases
        .iter()
        .filter(|p| p.kind == "frozen")
        .max_by_key(|p| p.steps)
        .expect("a frozen phase");
    assert!(long.steps >= 2, "{:?}", f.phases);
    let step = f.time_step;
    continues_bit_identically(
        "frozen",
        &base,
        long.start + step * f64::from(long.steps / 2),
    );
    let refresh = f
        .phases
        .iter()
        .find(|p| p.kind == "coupled" && p.reason != "start" && p.steps > 3)
        .expect("a refresh");
    continues_bit_identically("refresh", &base, refresh.start + 2.0);
}

/// Settings that cannot be used are refused with FRZ001.
#[test]
fn wrong_frozen_settings_are_refused() {
    let base = heated_wall(100.0, false, true);
    for (what, t) in [
        (
            "below the flow step",
            base.replace("time_step = 20.0", "time_step = 0.5"),
        ),
        (
            "with the acceleration",
            base.replace(
                "[energy.frozen_flow]\n",
                "[energy.acceleration]\ninertia_scale = 0.5\n\n[energy.frozen_flow]\n",
            ),
        ),
        (
            "a start word",
            base.replace(
                "[energy.frozen_flow]\n",
                "[energy.frozen_flow]\nstart = \"soon\"\n",
            ),
        ),
        (
            "a tolerance",
            base.replace(
                "[energy.frozen_flow]\n",
                "[energy.frozen_flow]\nsettle_tolerance = 2.0\n",
            ),
        ),
        (
            "an adaptive step",
            base.replace("[time]\n", "[time]\nadaptive = true\n"),
        ),
    ] {
        let (code, msg) = refused(&t);
        assert_eq!(code, diag::FRZ_SETTING, "{what}: {msg}");
    }
}
