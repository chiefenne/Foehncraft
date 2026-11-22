//! Cycle-averaged coupling of a periodic flow with its solids (D204,
//! docs/requirements/stage7-energy-buoyancy-cht.md §9d, checks V7.39 and
//! V7.40).

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_cyc_{name}_{}", std::process::id()))
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

/// The pulsating channel of the gallery on a coarser mesh.
fn coarse() -> String {
    text("cht_pulsating_channel_cycle.toml").replace("cells = [80, 20, 1]", "cells = [40, 10, 1]")
}

/// The same case fully coupled (`[energy.cycle] enabled = false`) until
/// `end` s.
fn coupled(end: f64) -> String {
    coarse()
        .replace("period = 50.0 ", "enabled = false\nperiod = 50.0 ")
        .replace(
            "end_time = 3000.0                     # an upper bound: the loop stops by itself",
            &format!("end_time = {end:.1}"),
        )
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

/// The mean of column `name` over the last period of a CSV (time in the
/// second column).
fn last_period_mean(dir: &Path, file: &str, name: &str, period: f64) -> f64 {
    let (h, rows) = csv(dir, file);
    let (k, end) = (col(&h, name), rows.last().unwrap()[1]);
    let sel: Vec<f64> = rows
        .iter()
        .filter(|r| r[1] > end - period + 1e-9)
        .map(|r| r[k])
        .collect();
    sel.iter().sum::<f64>() / sel.len() as f64
}

fn cycle(r: &FlowReport) -> foehncraft_flow::cycle::CycleSummary {
    r.solves[0]
        .transient
        .as_ref()
        .unwrap()
        .cycle_coupling
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

const PROBES: [&str; 4] = [
    "heated_face_T",
    "wall_upstream_T",
    "wall_middle_T",
    "wall_downstream_T",
];

/// V7.39 and V7.40: the pulsating laminar channel (inlet velocity and
/// temperature periodic) over a thick heated wall. The loop converges in a
/// few updates; its interface heat flow and wall temperatures are those of
/// the fully coupled transient run in its periodic state (within 1 %); the
/// heat balance of every step of the fluid phases and of every solid phase
/// closes to the solver tolerance.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn loop_matches_the_fully_coupled_periodic_state() {
    let dir_ref = out("coupled");
    let r_ref = run(&parse(&coupled(2500.0)), "coupled", &dir_ref, |_| {}).unwrap();
    let tr_ref = r_ref.solves[0].transient.as_ref().unwrap();
    assert!(tr_ref.cycle_coupling.is_none());
    let q_ref = last_period_mean(&dir_ref, "heat.csv", "interface-wall_heat", 50.0);
    let t_ref: Vec<f64> = PROBES
        .iter()
        .map(|p| last_period_mean(&dir_ref, "probes.csv", p, 50.0))
        .collect();
    let _ = std::fs::remove_dir_all(&dir_ref);
    // the fully coupled run is periodic: all of the 160 W go into the fluid
    assert!((q_ref - 160.0).abs() < 1e-3 * 160.0, "{q_ref}");
    let dt_scale = t_ref[0] - 300.0;
    assert!(dt_scale > 10.0, "{t_ref:?}");

    let dir = out("loop");
    let r = run(&parse(&coarse()), "loop", &dir, |_| {}).unwrap();
    let c = cycle(&r);
    assert!(c.converged, "{c:?}");
    assert!(c.updates.len() <= 5, "{} updates", c.updates.len());
    assert!(c.cycles <= 14, "{} cycles", c.cycles);
    let last = c.updates.last().unwrap();
    let q = last.heat_flow[0].1;
    assert!((q - q_ref).abs() < 0.01 * q_ref, "{q} against {q_ref}");
    // the held solids at the end are the loop's answer
    let (h, rows) = csv(&dir, "probes.csv");
    let end = rows.last().unwrap();
    for (p, tr) in PROBES.iter().zip(&t_ref) {
        let tl = end[col(&h, p)];
        assert!(
            (tl - tr).abs() < 0.01 * dt_scale,
            "{p}: loop {tl}, fully coupled {tr} (ΔT {dt_scale})"
        );
    }
    // V7.40: the ledgers
    let tr = r.solves[0].transient.as_ref().unwrap();
    assert!(tr.heat_relative_max.unwrap() < 1e-9, "{tr:?}");
    for u in &c.updates {
        assert!(u.solid_relative < 1e-9, "{u:?}");
    }
    // the outputs
    let (hc, rc) = csv(&dir, "cycle_coupling.csv");
    assert_eq!(rc.len(), c.updates.len());
    assert!(hc.contains(&"interface-wall_heat".to_string()));
    let (hf, rf) = csv(&dir, "cycle_faces.csv");
    assert_eq!(rf.len(), c.faces);
    let k = col(&hf, "heat_flux");
    let a = col(&hf, "area");
    let into_solid: f64 = rf.iter().map(|x| x[k] * x[a]).sum();
    assert!((into_solid + q).abs() < 1e-9 * q.abs(), "{into_solid} {q}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// The classic solid phase (`solid_model = "robin"`, the Robin pair of the
/// first fluid cell) contracts slowly in the same channel: after four
/// updates it is still far from the answer (CYC004), while its solid
/// balance closes.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn robin_solid_phase_is_slower_and_says_so() {
    let t = coarse()
        .replace("max_updates = 6", "max_updates = 4")
        .replace(
            "max_updates = 4\n",
            "max_updates = 4\nsolid_model = \"robin\"\n",
        );
    let dir = out("robin");
    let r = run(&parse(&t), "robin", &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let c = cycle(&r);
    assert!(!c.converged, "{c:?}");
    assert!(
        r.warnings.iter().any(|d| d.code == diag::CYC_MAX_UPDATES),
        "{:?}",
        r.warnings
    );
    let last = c.updates.last().unwrap();
    assert!(last.wall_change_max > 0.5, "{last:?}");
    for u in &c.updates {
        assert!(u.solid_relative < 1e-9, "{u:?}");
    }
}

/// CYC001: settings that cannot be used.
#[test]
fn wrong_settings_are_refused() {
    let t = text("cht_pulsating_channel_cycle.toml");
    let steady = format!(
        "{}\n[energy.cycle]\nperiod = 50.0\n",
        text("cht_heated_wall.toml")
    );
    let bad = [
        (steady, "steady case"),
        (
            t.replace("period = 50.0 ", "period = 51.0 "),
            "whole number",
        ),
        (
            t.replace("max_updates = 6", "max_updates = 6\nupdate = \"newton\""),
            "update must be",
        ),
        (
            t.replace("max_updates = 6", "max_updates = 6\nsolid_model = \"flux\""),
            "solid_model must be",
        ),
        (
            t.replace("max_updates = 6", "max_updates = 0"),
            "max_updates",
        ),
        (
            t.replace("max_updates = 6", "max_updates = 6\nrelaxation = 1.5"),
            "relaxation",
        ),
    ];
    for (case, what) in bad {
        let (code, msg) = refused(&case);
        assert_eq!(code, "CYC001", "{what}: {msg}");
        assert!(msg.contains(what), "{what}: {msg}");
    }
    let both = t.replace(
        "\n[initial]",
        "\n[energy.acceleration]\ninertia_scale = \"auto\"\n\n[initial]",
    );
    let (code, msg) = refused(&both);
    assert_eq!(code, "CYC001", "{msg}");
    assert!(msg.contains("[energy.acceleration]"), "{msg}");
}

/// Runs `text` once through and once in two parts (to `t_split` with the
/// restart file written, then continued from it) and checks that the
/// continuation is bit-identical: the last row of the probes, the heat
/// flow over the interfaces and, with a loop, every update.
fn continues_bit_identically(label: &str, text: &str, t_split: f64) {
    let whole = out(&format!("{label}_whole"));
    let (a, b) = (out(&format!("{label}_a")), out(&format!("{label}_b")));
    run(&parse(text), label, &whole, |_| {}).unwrap();
    let mut first = parse(text);
    first.time.as_mut().unwrap().end_time = Some(t_split);
    first.run.restart_final = true;
    run(&first, label, &a, |_| {}).unwrap();
    let mut second = parse(text);
    second.run.restart = Some(a.join("restart.h5"));
    run(&second, label, &b, |_| {}).unwrap();
    let mut files = vec!["probes.csv", "heat.csv"];
    if whole.join("cycle_coupling.csv").exists() {
        files.extend(["cycle_coupling.csv", "cycle_faces.csv"]);
    }
    for f in files {
        let (hw, rw) = csv(&whole, f);
        let (hb, rb) = csv(&b, f);
        assert_eq!(hw, hb, "{label}: {f}");
        let (lw, lb) = (rw.last().unwrap(), rb.last().unwrap());
        assert!(
            lw.iter()
                .zip(lb)
                .all(|(x, y)| x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan())),
            "{label}: {f}\n{lw:?}\n{lb:?}"
        );
        if f == "cycle_coupling.csv" {
            assert_eq!(rw.len(), rb.len(), "{label}: the updates");
        }
    }
    for d in [whole, a, b] {
        let _ = std::fs::remove_dir_all(&d);
    }
}

/// The CHT restart (§9d): the temperatures of the solids and their old
/// time level go into restart.h5, and with `[energy.cycle]` the state of
/// the loop (the averages of the running fluid phase, the history of the
/// quasi-Newton update, the held means, the skin swing and the fully
/// coupled cycles at the end). A run continued from it is bit-identical to
/// one run through: a plain transient conjugate run, the loop restarted in
/// the middle of a fluid phase, and the loop with the skin swing restarted
/// in the middle of a fluid phase and in the fully coupled cycles.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn cht_restart_continues_bit_identically() {
    let plain = text("cht_heated_wall_transient.toml")
        .replace("cells = [80, 30, 1]", "cells = [40, 15, 1]")
        .replace("end_time = 3000.0", "end_time = 200.0");
    continues_bit_identically("plain", &plain, 100.0);
    continues_bit_identically("loop", &coarse(), 175.0);
    let skin = text("cht_pulsating_channel_skin.toml")
        .replace("cells = [80, 80, 1]", "cells = [40, 40, 1]")
        .replace("max_updates = 10", "max_updates = 4")
        .replace("resolved_cycles = 3", "resolved_cycles = 2");
    // four updates: the loop ends at 500 s, the fully coupled cycles run
    // to 600 s
    continues_bit_identically("skin_phase", &skin, 425.0);
    continues_bit_identically("skin_resolved", &skin, 575.0);
}

/// The steady CHT restart: a steady conjugate run continued from its
/// restart file starts from the temperatures of the solids it ended with.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn steady_cht_restart_keeps_the_solid_temperatures() {
    let t = text("cht_heated_wall.toml").replace("cells = [80, 30, 1]", "cells = [40, 15, 1]");
    let (a, b) = (out("steady_a"), out("steady_b"));
    let mut first = parse(&t);
    first.run.restart_final = true;
    let r1 = run(&first, "steady", &a, |_| {}).unwrap();
    let mut second = parse(&t);
    second.run.restart = Some(a.join("restart.h5"));
    let r2 = run(&second, "steady", &b, |_| {}).unwrap();
    // converged already: the continuation stops at its first iteration
    // (the count goes on from the restart file)
    let (n1, n2) = (r1.solves[0].iterations, r2.solves[0].iterations);
    assert!(
        n2 <= n1 + 2,
        "{n2} iterations after the restart, {n1} before"
    );
    let (h, ra) = csv(&a, "line_wall_top.csv");
    let (_, rb) = csv(&b, "line_wall_top.csv");
    let k = col(&h, "temperature");
    assert_eq!(ra.len(), rb.len());
    for (x, y) in ra.iter().zip(&rb) {
        assert!((x[k] - y[k]).abs() < 1e-6, "{x:?} {y:?}");
    }
    for d in [a, b] {
        let _ = std::fs::remove_dir_all(&d);
    }
}

/// The skin swing (`skin = "effusivity"`) on a thin skin: the loop with
/// the swing converges; the fully coupled cycles after it give the heat
/// flow of the heater and a swing at the surface; the summary and
/// cycle_faces.csv carry the swing.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn skin_swing_and_resolved_cycles() {
    let t = text("cht_pulsating_channel_skin.toml")
        .replace("cells = [80, 80, 1]", "cells = [40, 40, 1]");
    let dir = out("skin");
    let r = run(&parse(&t), "skin", &dir, |_| {}).unwrap();
    let c = cycle(&r);
    assert!(c.converged, "{:?}", c.updates.last());
    assert!(
        !r.warnings.iter().any(|d| d.code == diag::CYC_SKIN),
        "{:?}",
        r.warnings
    );
    let sk = c.skin.as_ref().unwrap();
    assert!(sk.harmonics >= 5, "{sk:?}");
    assert!(sk.swing_max > 1.0, "{sk:?}");
    let last = c.updates.last().unwrap();
    assert!(last.swing_change.unwrap() <= 0.01, "{last:?}");
    let res = c.resolved.as_ref().unwrap();
    assert_eq!(res.cycles, 3);
    let q = res.heat_flow[0].1;
    assert!((q - 160.0).abs() < 0.01 * 160.0, "{q}");
    // the swing of the resolved cycles against the loop's estimate
    assert!(
        (res.swing_max - sk.swing_max).abs() < 0.1 * res.swing_max,
        "{res:?} {sk:?}"
    );
    let (h, rows) = csv(&dir, "cycle_faces.csv");
    let (lo, hi) = (
        col(&h, "resolved_t_wall_min"),
        col(&h, "resolved_t_wall_max"),
    );
    assert!(rows.iter().all(|x| x[hi] >= x[lo]));
    let _ = std::fs::remove_dir_all(&dir);
}

/// CYC003: the skin swing on the thick-skin wall of the gallery case, whose
/// skin depth reaches through the wall, warns.
#[test]
fn skin_deeper_than_the_wall_warns() {
    let t = coarse()
        .replace("max_updates = 6", "max_updates = 6\nskin = \"effusivity\"")
        .replace(
            "end_time = 3000.0                     # an upper bound: the loop stops by itself",
            "end_time = 5.0",
        );
    let dir = out("thick");
    let r = run(&parse(&t), "thick", &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let w = r
        .warnings
        .iter()
        .find(|d| d.code == diag::CYC_SKIN)
        .unwrap_or_else(|| panic!("{:?}", r.warnings));
    assert!(w.message.contains("wall"), "{}", w.message);
}
