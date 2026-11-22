//! Thermal acceleration of solid regions in transient conjugate heat
//! transfer (D213, docs/requirements/stage7-energy-buoyancy-cht.md §9e,
//! checks V7.41 and V7.42).

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_acc_{name}_{}", std::process::id()))
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

/// The heated wall on a coarser mesh.
fn coarse(t: &str) -> String {
    t.replace("cells = [80, 30, 1]", "cells = [40, 15, 1]")
}

fn temperature(r: &FlowReport) -> &foehncraft_flow::scalar::ScalarSummary {
    r.solves[0]
        .scalars
        .iter()
        .find(|s| s.name == "temperature")
        .unwrap()
}

fn boundary_mean(r: &FlowReport, n: &str) -> f64 {
    temperature(r)
        .boundaries
        .iter()
        .find(|b| b.boundary == n)
        .unwrap()
        .mean
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
        .map(|l| l.split(',').map(|v| v.parse().unwrap()).collect())
        .collect();
    (head, rows)
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

/// V7.41: the heated wall with `inertia_scale = "auto"`. The measured time
/// constant is that of the case (C/G ≈ 290 s, the case file's estimate
/// 300 s), the wall settles and is restored within the run, the end state
/// is the steady one (as the unaccelerated run after 3000 s), and the
/// balance of the solved equation closes in every step, with the
/// inertia-scaling term exactly (1 − f)/f of the solved wall storage.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn auto_acceleration_reaches_the_steady_wall_sooner() {
    let steady = {
        let dir = out("steady");
        let r = run(
            &parse(&coarse(&text("cht_heated_wall.toml"))),
            "steady",
            &dir,
            |_| {},
        )
        .unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        r
    };
    let case = parse(&coarse(&text("cht_heated_wall_accelerated.toml")));
    let dir = out("auto");
    let r = run(&case, "auto", &dir, |_| {}).unwrap();
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let tr = r.solves[0].transient.as_ref().unwrap();
    let a = tr.acceleration.as_ref().unwrap();
    assert_eq!(a.mode, "auto");
    let w = &a.regions[0];
    let tau = w.time_constant.unwrap();
    // C/G of the case: ρ c_p V of the wall cells over about 54 W/(m² K) × A
    // (300 s for the wall of 0.2 m; the coarse mesh's wall cells reach to
    // 0.24 m)
    assert!(
        (tau / (300.0 * w.heat_capacity / 2560.0) - 1.0).abs() < 0.1,
        "{a:?}"
    );
    // the scaled time constant is one window (10 steps of 5 s)
    assert!((w.factor * tau / 50.0 - 1.0).abs() < 1e-9, "{a:?}");
    assert!(a.settled_at.is_some() && a.restored_at.is_some(), "{a:?}");
    assert!(
        a.restored_at.unwrap() < 600.0 && a.history_physical,
        "{a:?}"
    );
    assert!(w.inertia_scaling_heat > 0.0);
    assert!(r.warnings.iter().all(|d| d.code != diag::ACC_UNRESTORED));
    // the end state is the steady one
    let s = temperature(&steady);
    let dt = s.max - s.min;
    for n in ["ymin", "interface-wall", "outlet"] {
        let (x, y) = (boundary_mean(&steady, n), boundary_mean(&r, n));
        assert!((x - y).abs() < 2e-3 * dt, "{n}: {y} against the steady {x}");
    }
    // the ledger: every step balanced; while accelerated the inertia
    // scaling is (1 − f)/f of the wall's solved storage, zero after
    assert!(tr.heat_relative_max.unwrap() < 1e-9, "{tr:?}");
    let (head, rows) = csv(&dir, "heat.csv");
    let col = |n: &str| head.iter().position(|h| h == n).unwrap();
    let (cs, ci, cr, ct) = (
        col("wall_storage"),
        col("inertia_scaling"),
        col("relative"),
        col("time"),
    );
    let ratio = (1.0 - w.factor) / w.factor;
    let (t0, t1) = (a.accelerated_from.unwrap(), a.restored_at.unwrap());
    for row in &rows {
        assert!(row[cr] < 1e-9, "{row:?}");
        let t = row[ct];
        if t > t0 + 1e-6 && t < t1 - 1e-6 {
            assert!(
                (row[ci] - ratio * row[cs]).abs() <= 1e-12 * row[cs].abs().max(1e-30),
                "{row:?}"
            );
        } else if t <= t0 + 1e-6 || t > t1 + 1e-6 {
            assert_eq!(row[ci], 0.0, "{row:?}");
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// V7.42: a fixed factor heats the wall about 1/f times faster at the start
/// (the true storage of the first step against the true wall's), and
/// `restore` at a time brings the true ρ c_p back then.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn fixed_factor_scales_the_heat_up_and_restores_at_a_time() {
    let base = coarse(&text("cht_heated_wall_transient.toml"))
        .replace("end_time = 3000.0", "end_time = 100.0");
    let acc = base.replace(
        "[energy]\ninitial = 300.0",
        "[energy.acceleration]\ninertia_scale = 0.25\nrestore = 50.0\n\n[energy]\ninitial = 300.0",
    );
    // the true storage of the wall in the first step (W): the same heat
    // goes into the wall, which warms 1/f times faster
    let first = |t: &str, name: &str| -> (FlowReport, f64) {
        let dir = out(name);
        let r = run(&parse(t), name, &dir, |_| {}).unwrap();
        let (head, rows) = csv(&dir, "heat.csv");
        let c = |n: &str| head.iter().position(|h| h == n);
        let s =
            rows[0][c("wall_storage").unwrap()] + c("inertia_scaling").map_or(0.0, |i| rows[0][i]);
        let _ = std::fs::remove_dir_all(&dir);
        (r, s)
    };
    let (_, s1) = first(&base, "plain");
    let (r, s2) = first(&acc, "fixed");
    // (a little less than 4: the faster wall passes a little more of the
    // heater's power on to the flow already in the first step)
    let x = s2 / s1 / 4.0;
    assert!(x > 0.9 && x <= 1.0, "{s2} against {s1}");
    let a = r.solves[0]
        .transient
        .as_ref()
        .unwrap()
        .acceleration
        .clone()
        .unwrap();
    assert_eq!(a.mode, "factor");
    assert_eq!(a.regions[0].factor, 0.25);
    assert_eq!(a.restore_started, Some(50.0));
    assert_eq!(a.restored_at, Some(50.0));
    assert!(a.settled_at.is_none());
}

/// ACC002: a run that ends while the solids are still accelerated says so,
/// and its history is not physical.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn unrestored_run_warns() {
    let t = coarse(&text("cht_heated_wall_accelerated.toml"))
        .replace("end_time = 600.0", "end_time = 100.0");
    let dir = out("short");
    let r = run(&parse(&t), "short", &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    assert!(
        r.warnings.iter().any(|d| d.code == diag::ACC_UNRESTORED),
        "{:?}",
        r.warnings
    );
    let a = r.solves[0]
        .transient
        .as_ref()
        .unwrap()
        .acceleration
        .clone()
        .unwrap();
    assert!(!a.history_physical && a.restored_at.is_none(), "{a:?}");
}

/// ACC001: settings that cannot be used.
#[test]
fn wrong_settings_are_refused() {
    let t = text("cht_heated_wall_accelerated.toml");
    for (from, to, what) in [
        ("inertia_scale = \"auto\"", "inertia_scale = 1.5", "(0, 1]"),
        (
            "inertia_scale = \"auto\"",
            "inertia_scale = \"fast\"",
            "\"auto\"",
        ),
        (
            "inertia_scale = \"auto\"",
            "inertia_scale = 0.1\nregions = [\"fluid\"]",
            "not a solid region",
        ),
        (
            "inertia_scale = \"auto\"",
            "inertia_scale = 0.1\nrestore = \"later\"",
            "\"settled\"",
        ),
        (
            "inertia_scale = \"auto\"",
            "inertia_scale = 0.1\nwindow_steps = 1",
            "window_steps",
        ),
    ] {
        let (code, msg) = refused(&t.replace(&format!("\n{from}"), &format!("\n{to}")));
        assert_eq!(code, "ACC001", "{msg}");
        assert!(msg.contains(what), "{msg}");
    }
    // a steady case
    let steady = text("cht_heated_wall.toml").replace(
        "[energy]",
        "[energy.acceleration]\ninertia_scale = 0.1\n\n[energy]",
    );
    let (code, msg) = refused(&steady);
    assert_eq!(code, "ACC001", "{msg}");
    assert!(msg.contains("steady"), "{msg}");
}
