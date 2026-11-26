//! Stage 6, part 6b: the k-ω SST model (D130), its wall treatments (D134),
//! inlet values and checks (D137), docs/verification/stage6.md.
//!
//! The flat plate and the channel of the gallery (`tmr_plate_sst.toml`,
//! `channel_sst.toml`) are the validation cases; the tests here use small
//! meshes.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_sst_{name}_{}", std::process::id()))
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

/// A case file with `\n` line ends (a Windows checkout has `\r\n`).
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

/// The channel of the gallery, small: 24 cells across, few iterations.
fn small_channel(extra: &str) -> String {
    text("channel_sst.toml")
        .replace("cells = [2, 96, 1]", "cells = [2, 24, 1]")
        .replace("ratio = 60.0", "ratio = 8.0")
        .replace("max_iterations = 3000", "max_iterations = 40")
        + "\n[verification]\nrequire_convergence = false\n"
        + extra
}

/// Homogeneous decay in a box at rest without walls (specification 9.2):
/// F₁ = 0, so ω = ω₀/(1 + β₂ω₀t) and k = k₀(1 + β₂ω₀t)^(−β*/β₂). BDF2 is
/// second order in time (the first step is BDF1).
#[cfg(feature = "hdf5")]
#[test]
fn homogeneous_decay_is_second_order_in_time() {
    use foehncraft_flow::restart;
    let (k0, w0, t_end): (f64, f64, f64) = (1.0, 10.0, 0.5);
    let (beta, beta_star): (f64, f64) = (0.0828, 0.09);
    let exact_w = w0 / (1.0 + beta * w0 * t_end);
    let exact_k = k0 * (1.0 + beta * w0 * t_end).powf(-beta_star / beta);
    let mut errors = Vec::new();
    for n in [10usize, 20, 40] {
        let c = parse(&format!(
            r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.5]
cells = [2, 2, 1]
[fluid]
viscosity = 1e-3
[turbulence]
model = "k-omega-sst"
[initial]
velocity = [0.0, 0.0, 0.0]
k = {k0}
omega = {w0}
[boundary.xmin]
type = "symmetry"
[boundary.xmax]
type = "symmetry"
[boundary.ymin]
type = "symmetry"
[boundary.ymax]
type = "symmetry"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver.turbulence]
tolerance = 1e-12
[time]
scheme = "bdf2"
time_step = {dt}
steps = {n}
[output]
vtkhdf = false
wall_csv = false
"#,
            dt = t_end / n as f64
        ));
        let d = out(&format!("decay{n}"));
        let r = run(&c, "decay", &d, |_| {}).unwrap();
        assert!(r.failures.is_empty(), "{:?}", r.failures);
        let rs = restart::read(&d.join("restart.h5")).unwrap();
        let t = rs.state.turbulence.unwrap();
        let _ = std::fs::remove_dir_all(&d);
        // uniform: every cell the same (to round-off)
        assert!(
            t.k.iter().all(|&k| (k - t.k[0]).abs() < 1e-14 * t.k[0]),
            "{:?}",
            t.k
        );
        errors.push((
            (t.k[0] - exact_k).abs() / exact_k,
            (t.omega[0] - exact_w).abs() / exact_w,
        ));
    }
    for i in 1..errors.len() {
        let pk = (errors[i - 1].0 / errors[i].0).log2();
        let pw = (errors[i - 1].1 / errors[i].1).log2();
        eprintln!("decay: errors {:?}, orders k {pk:.3}, ω {pw:.3}", errors[i]);
        assert!((1.8..=2.2).contains(&pk), "order of k {pk}");
        assert!((1.8..=2.2).contains(&pw), "order of ω {pw}");
    }
    assert!(errors[2].0 < 1e-3 && errors[2].1 < 1e-3, "{errors:?}");
}

/// A turbulent run continued from a restart file ends bit-identically to
/// the run without the interruption (D42): k, ω, ν_t and p̂ are part of the
/// restart file.
#[cfg(feature = "hdf5")]
#[test]
fn sst_restart_continues_bit_identically() {
    use foehncraft_flow::restart;
    let case = |extra: &str| {
        parse(&small_channel(&format!(
            "[time]\nscheme = \"bdf2\"\ntime_step = 0.002\nsteps = 10\n{extra}\n"
        )))
    };
    let da = out("restart_full");
    let ra = run(&case("[run]\nrestart_at = [5]"), "full", &da, |_| {}).unwrap();
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let from = da.join("restart_000005.h5");
    let db = out("restart_cont");
    let rb = run(
        &case(&format!("[run]\nrestart = '{}'", from.display())),
        "cont",
        &db,
        |_| {},
    )
    .unwrap();
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let a = restart::read(&da.join("restart.h5")).unwrap();
    let b = restart::read(&db.join("restart.h5")).unwrap();
    let (ta, tb) = (a.state.turbulence.unwrap(), b.state.turbulence.unwrap());
    assert!(ta.k_old.is_some());
    assert_eq!(ta, tb, "the continued turbulence state differs");
    assert_eq!(a.state.u.x, b.state.u.x);
    assert_eq!(a.state.p, b.state.p);
    assert_eq!(a.state.periodic_force, b.state.periodic_force);
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// One and four threads give the same bits (D31): residuals of every
/// equation, the wall shear stress and the forces, on a plate mesh above
/// the 20 000 cells below which one thread is used.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn sst_results_are_bitwise_equal_on_1_and_4_threads() {
    let mk = |threads: usize| {
        parse(
            &(text("tmr_plate_sst.toml")
                .replace("cells = [32, 24, 1]", "cells = [160, 128, 1]")
                .replace("cells = 8, ratio", "cells = 40, ratio")
                .replace("cells = 24, ratio", "cells = 120, ratio")
                .replace("levels = [1]", "levels = [0]")
                .replace("max_iterations = 2000", "max_iterations = 8")
                .replace("derived = [\"turbulence-features\"]", "")
                + &format!("\n[run]\nthreads = {threads}\n")),
        )
    };
    let mut c1 = mk(1);
    c1.verification.require_convergence = Some(false);
    let mut c4 = mk(4);
    c4.verification.require_convergence = Some(false);
    let a = run_case(&c1, "threads1");
    let b = run_case(&c4, "threads4");
    assert!(a.solves[0].cells >= 20_000);
    let (ra, rb) = (&a.solves[0].residuals, &b.solves[0].residuals);
    assert_eq!(ra.momentum.map(f64::to_bits), rb.momentum.map(f64::to_bits));
    assert_eq!(ra.continuity.to_bits(), rb.continuity.to_bits());
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&ra.turbulence.residual), bits(&rb.turbulence.residual));
    let (wa, wb) = (&a.solves[0].wall_probes[0], &b.solves[0].wall_probes[0]);
    assert_eq!(wa.magnitude.to_bits(), wb.magnitude.to_bits());
    for (x, y) in a.solves[0].boundaries.iter().zip(&b.solves[0].boundaries) {
        assert_eq!(
            x.force.map(f64::to_bits),
            y.force.map(f64::to_bits),
            "{}",
            x.name
        );
    }
}

/// The channel at Re_τ 395 on a coarse mesh: the driving force gives
/// u_τ, the run stays positive (no floors late in the run) and logs its
/// turbulence summary; the derived features are finite.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn coarse_channel_runs_and_reports() {
    let c = parse(
        &small_channel("")
            .replace("max_iterations = 40", "max_iterations = 400")
            .replace(
                "wall_csv = true",
                "wall_csv = true\nderived = [\"turbulence-features\", \"y-plus\"]",
            ),
    );
    let r = run_case(&c, "coarse_channel");
    let s = &r.solves[0];
    let t = s.turbulence.as_ref().expect("turbulence summary");
    assert_eq!(t.health.floor_iterations_late, [0, 0], "{:?}", t.health);
    assert_eq!(t.walls.len(), 2);
    assert!(t.walls[0].mean > 0.0);
    assert!(
        s.residuals
            .turbulence
            .residual
            .iter()
            .all(|r| r.is_finite())
    );
    assert!(
        !r.warnings.iter().any(|w| w.code == "TURB002"),
        "{:?}",
        r.warnings
    );
}

/// Turbulence values and keys that are wrong are refused with their codes;
/// a laminar case with turbulence values warns.
#[test]
fn bad_turbulence_input_is_refused() {
    let t = small_channel("");
    // k without ω
    let (code, msg) = refused(&t.replace(
        "viscosity_ratio = 100.0",
        "viscosity_ratio = 100.0\nk = 1.0",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // an unknown constant
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"low-re\"",
        "wall_treatment = \"low-re\"\nconstants = { beta_strar = 0.09 }",
    ));
    assert_eq!(code, "CASE002", "{msg}");
    assert!(msg.contains("beta_star"), "{msg}");
    // central differencing of k and ω in a steady run
    let (code, msg) = refused(&format!("{t}\n[solver.turbulence]\nscheme = \"central\"\n"));
    assert_eq!(code, "TURB006", "{msg}");
    // a negative limiter constant
    let (code, msg) = refused(&format!("{t}\n[solver.turbulence]\nlimiter_k = -1.0\n"));
    assert_eq!(code, "CASE004", "{msg}");
}

/// Inlet turbulence: conflicting keys are refused (BC029), a missing set
/// gives the default note (TURB009); backflow values at a pressure outlet
/// are used (no BC028), ζ there with SST is not (BC028).
#[test]
fn inlet_turbulence_checks() {
    let base = text("tmr_plate_sst.toml")
        .replace(
            "levels = [1]\nrequire_convergence = true",
            "levels = [0]\nrequire_convergence = false",
        )
        .replace("max_iterations = 2000", "max_iterations = 1");
    let (code, msg) = refused(&base.replace(
        "k = 2.25e-7\nomega = 125.0\n[boundary.xmax]",
        "k = 2.25e-7\nomega = 125.0\nturbulence_intensity = 0.01\n[boundary.xmax]",
    ));
    assert_eq!(code, "BC029", "{msg}");
    let r = run_case(
        &parse(&base.replace(
            "k = 2.25e-7\nomega = 125.0\n[boundary.xmax]",
            "[boundary.xmax]",
        )),
        "inlet_default",
    );
    assert!(
        r.findings.iter().any(|f| f.code == "TURB009"),
        "{:?}",
        r.findings
    );
    let r = run_case(
        &parse(&base.replace(
            "[boundary.ymax]\ntype = \"static-pressure-outlet\"",
            "[boundary.ymax]\ntype = \"static-pressure-outlet\"\nturbulence_intensity = 0.01",
        )),
        "outlet_values",
    );
    assert!(
        !r.findings.iter().any(|f| f.code == "BC028"),
        "{:?}",
        r.findings
    );
    let r = run_case(
        &parse(&base.replace(
            "[boundary.ymax]\ntype = \"static-pressure-outlet\"",
            "[boundary.ymax]\ntype = \"static-pressure-outlet\"\nturbulence_intensity = 0.01\nzeta = 0.5",
        )),
        "outlet_zeta",
    );
    assert!(
        r.findings.iter().any(|f| f.code == "BC028"),
        "{:?}",
        r.findings
    );
}

/// The smooth switches of specification 2.6 (`turbulence.smooth_limits`,
/// for the adjoint) change the converged flat-plate C_f (wall probe at
/// x = 0.97) and the drag of the plate by less than 0.1 % against the
/// published hard min and max (level 1 of `tmr_plate_sst.toml`, 64 × 48
/// cells: 0.090 % for C_f; level 0, 32 × 24 cells, gives 0.1002 %).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn smooth_limits_change_the_plate_friction_below_a_tenth_percent() {
    let mk = |smooth: bool| {
        parse(
            &text("tmr_plate_sst.toml")
                .replace("derived = [\"turbulence-features\"]", "")
                .replace(
                    "wall_treatment = \"low-re\"",
                    &format!("wall_treatment = \"low-re\"\nsmooth_limits = {smooth}"),
                ),
        )
    };
    let hard = run_case(&mk(false), "plate_hard");
    let smooth = run_case(&mk(true), "plate_smooth");
    let (h, s) = (&hard.solves[0], &smooth.solves[0]);
    assert!(hard.failures.is_empty(), "{:?}", hard.failures);
    assert!(smooth.failures.is_empty(), "{:?}", smooth.failures);
    let (th, ts) = (h.wall_probes[0].magnitude, s.wall_probes[0].magnitude);
    let d_cf = (ts - th).abs() / th;
    let drag = |r: &foehncraft_flow::case::Solve| {
        r.boundaries
            .iter()
            .find(|b| b.name == "plate")
            .map(|b| b.force[0])
            .unwrap()
    };
    let (dh, ds) = (drag(h), drag(s));
    let d_cd = (ds - dh).abs() / dh.abs();
    eprintln!("C_f relative change {d_cf:.3e}, drag {d_cd:.3e} (tau hard {th:e})");
    assert!(d_cf < 1e-3, "C_f changes by {d_cf:e}");
    assert!(d_cd < 1e-3, "drag changes by {d_cd:e}");
}

/// A duct between two pressure outlets, the left one at the higher
/// pressure, so the flow enters through it (BC028, spec 7): with backflow
/// values the entering faces take them (k_b = k, a = 0), without them the
/// cell's values; the right outlet, where the flow leaves, stays zero
/// gradient. A run continued from a restart file of step 5 ends bit for
/// bit like the run without the interruption: the backflow state is set
/// again from the stored boundary values, and p̂ = p + ⅔ρk of the pressure
/// boundaries from the loaded k (before this, from the start fields).
#[cfg(feature = "hdf5")]
#[test]
fn outlet_backflow_takes_the_given_turbulence_values() {
    use foehncraft_flow::restart;
    let duct = |values: &str, extra: &str| {
        let t = format!(
            r#"case_format = "0.6"
[mesh.generate]
type = "hex"
dimension = 2
upper = [4.0, 1.0, 0.1]
cells = [16, 12, 1]
[fluid]
density = 1.0
viscosity = 1.0e-3
[turbulence]
model = "k-omega-sst"
wall_treatment = "low-re"
[initial]
velocity = [0.0, 0.0, 0.0]
k = 1.0e-5
omega = 1.0
[boundary.xmin]
type = "static-pressure-outlet"
pressure = 0.5
{values}
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[time]
scheme = "bdf2"
time_step = 0.05
steps = 10
[verification]
require_convergence = false
[output]
vtkhdf = false
{extra}
"#
        );
        parse(&t)
    };
    let values = "k = 1.0e-2\nomega = 20.0";
    let ni_nb = |d: &Path| restart::read(&d.join("restart.h5")).unwrap();
    // with values: the left faces take k = 10⁻², the right ones 0
    let da = out("backflow_full");
    let ra = run(
        &duct(values, "[run]\nrestart_at = [5]"),
        "full",
        &da,
        |_| {},
    )
    .unwrap();
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    assert!(
        !ra.findings.iter().any(|f| f.code == "BC028"),
        "{:?}",
        ra.findings
    );
    let a = ni_nb(&da);
    let (mesh, _, _) =
        foehncraft_flow::derived::evaluate(&duct(values, ""), &da.join("restart.h5"), &[]).unwrap();
    let ni = mesh.n_interior_faces();
    let nc = mesh.n_cells;
    let local = |r: restart::Restart| {
        r.to_local(&mesh.cell_id[..nc], &mesh.face_id, ni)
            .unwrap()
            .0
            .state
            .turbulence
            .unwrap()
    };
    let ta = local(ni_nb(&da));
    let faces = |name: &str| {
        mesh.boundaries
            .iter()
            .find(|p| p.name == name)
            .unwrap()
            .range()
            .map(|f| f - ni)
            .collect::<Vec<_>>()
    };
    for fb in faces("xmin") {
        assert_eq!(ta.k_b[fb], 1.0e-2, "k of the backflow");
        assert_eq!(ta.omega_b[fb], 20.0, "ω of the backflow");
    }
    for fb in faces("xmax") {
        assert_eq!((ta.k_b[fb], ta.omega_b[fb]), (0.0, 0.0), "zero gradient");
    }
    // the entering k reaches the first cells (far above the start 10⁻⁵)
    let k_first: f64 = faces("xmin")
        .iter()
        .map(|&fb| ta.k[mesh.owner[ni + fb] as usize])
        .fold(f64::INFINITY, f64::min);
    // without values: the cell's k
    let dn = out("backflow_none");
    let rn = run(&duct("", ""), "none", &dn, |_| {}).unwrap();
    assert!(rn.failures.is_empty(), "{:?}", rn.failures);
    let tn = local(ni_nb(&dn));
    for fb in faces("xmin") {
        assert_eq!(tn.k_b[fb], 0.0, "zero gradient without backflow values");
    }
    let k_none: f64 = faces("xmin")
        .iter()
        .map(|&fb| tn.k[mesh.owner[ni + fb] as usize])
        .fold(0.0, f64::max);
    assert!(
        k_first > 10.0 * k_none,
        "k of the first cells {k_first:e} with values, {k_none:e} without"
    );
    // the continued run
    let from = da.join("restart_000005.h5");
    let db = out("backflow_cont");
    let rb = run(
        &duct(values, &format!("[run]\nrestart = '{}'", from.display())),
        "cont",
        &db,
        |_| {},
    )
    .unwrap();
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let b = ni_nb(&db);
    assert_eq!(
        a.state.turbulence, b.state.turbulence,
        "the continued turbulence state differs"
    );
    assert_eq!(a.state.u.x, b.state.u.x);
    assert_eq!(a.state.p, b.state.p);
    for d in [&da, &dn, &db] {
        let _ = std::fs::remove_dir_all(d);
    }
}

/// The SST constants of `[turbulence.constants]` reach the run (they were
/// checked against the ζ-f names at the start of the run and `gamma2` was
/// refused): homogeneous decay with β₂ = 0.1 follows ω = ω₀/(1 + β₂ω₀t).
#[cfg(feature = "hdf5")]
#[test]
fn sst_constant_overrides_reach_the_run() {
    use foehncraft_flow::restart;
    let (w0, t_end, beta): (f64, f64, f64) = (10.0, 0.5, 0.1);
    let c = parse(&format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.5]
cells = [2, 2, 1]
[fluid]
viscosity = 1e-3
[turbulence]
model = "k-omega-sst"
[turbulence.constants]
beta2 = {beta}
gamma2 = 0.44
[initial]
velocity = [0.0, 0.0, 0.0]
k = 1.0
omega = {w0}
[boundary.xmin]
type = "symmetry"
[boundary.xmax]
type = "symmetry"
[boundary.ymin]
type = "symmetry"
[boundary.ymax]
type = "symmetry"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver.turbulence]
tolerance = 1e-12
[time]
scheme = "bdf2"
time_step = {dt}
steps = 40
[output]
vtkhdf = false
wall_csv = false
"#,
        dt = t_end / 40.0
    ));
    let d = out("constants");
    let r = run(&c, "constants", &d, |_| {}).unwrap();
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let rs = restart::read(&d.join("restart.h5")).unwrap();
    let t = rs.state.turbulence.unwrap();
    let _ = std::fs::remove_dir_all(&d);
    let exact = w0 / (1.0 + beta * w0 * t_end);
    assert!(
        (t.omega[0] - exact).abs() < 1e-3 * exact,
        "{} vs {exact}",
        t.omega[0]
    );
}
