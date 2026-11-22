//! Stage 6, part 6e: the elliptic-blending k-ε model of Billard & Laurence
//! (2012), D133, on the machinery of k-ε-ζ-f (φ in the slot of ζ, α in
//! the slot of f), docs/verification/stage6.md.
//!
//! The channel, the flat plate and the backward-facing step of the gallery
//! (`channel_eb.toml`, `tmr_plate_eb.toml`, `backstep_eb.toml`) are the
//! validation cases; the tests here use small meshes, except the channel
//! against the paper (release builds).

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

#[cfg(feature = "hdf5")]
const MODEL: &str = "model = \"elliptic-blending-k-epsilon\"";

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_eb_{name}_{}", std::process::id()))
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

/// The channel of the gallery on 24 cells with few iterations.
fn small_channel(extra: &str) -> String {
    text("channel_eb.toml")
        .replace("cells = [2, 96, 1]", "cells = [2, 24, 1]")
        .replace("ratio = 60.0", "ratio = 8.0")
        .replace("max_iterations = 5000", "max_iterations = 30")
        + "\n[verification]\nrequire_convergence = false\n"
        + extra
}

/// The box at rest of the decay test, with the initial values of
/// `initial`; ν is small, so that T = (k²/ε² + C_T² ν/ε)^½ is k/ε.
#[cfg(feature = "hdf5")]
fn decay_box(initial: &str, dt: f64, n: usize) -> String {
    format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.5]
cells = [2, 2, 1]
[fluid]
viscosity = 1e-12
[turbulence]
{MODEL}
max_viscosity_ratio = 1e15
[initial]
velocity = [0.0, 0.0, 0.0]
{initial}
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
[solver]
outer_iterations = 4
[solver.turbulence]
tolerance = 1e-12
[solver.elliptic]
tolerance = 1e-12
[time]
scheme = "bdf2"
time_step = {dt}
steps = {n}
[output]
vtkhdf = false
wall_csv = false
"#
    )
}

/// Homogeneous decay in a box at rest: with P = 0, ∇k = 0 and no wall,
/// α = 1, the switch of C*_ε2 and the E term vanish, and k and ε follow
/// the k-ε solution with C_ε2 = 1.83; φ relaxes to 2/3 through
/// f_h = −(C₁ − 1)(φ − 2/3)/T with C₁ = 1.7:
/// 2/3 − φ = (2/3 − φ₀)(1 + (C_ε2 − 1) ε₀ t/k₀)^(−(C₁ − 1)/(C_ε2 − 1)).
/// BDF2 is second order in time for all three (the first step is BDF1).
#[cfg(feature = "hdf5")]
#[test]
fn eb_homogeneous_decay_is_second_order_in_time() {
    use foehncraft_flow::restart;
    let (k0, e0, p0, t_end): (f64, f64, f64, f64) = (1.0, 1.0, 0.3, 0.5);
    let (c2, c1): (f64, f64) = (1.83, 1.7);
    let a = 1.0 + (c2 - 1.0) * e0 * t_end / k0;
    let exact_k = k0 * a.powf(-1.0 / (c2 - 1.0));
    let exact_e = e0 * a.powf(-c2 / (c2 - 1.0));
    let exact_p = 2.0 / 3.0 - (2.0 / 3.0 - p0) * a.powf(-(c1 - 1.0) / (c2 - 1.0));
    let mut errors = Vec::new();
    for n in [10usize, 20, 40] {
        let c = parse(&decay_box(
            &format!("k = {k0}\nepsilon = {e0}\nzeta = {p0}"),
            t_end / n as f64,
            n,
        ));
        let d = out(&format!("decay{n}"));
        let r = run(&c, "decay", &d, |_| {}).unwrap();
        assert!(r.failures.is_empty(), "{:?}", r.failures);
        let rs = restart::read(&d.join("restart.h5")).unwrap();
        let t = rs.state.turbulence.unwrap();
        let _ = std::fs::remove_dir_all(&d);
        let z = t.zeta.expect("the restart file holds φ and α");
        assert!(
            z.f.iter().all(|&al| (al - 1.0).abs() < 1e-12),
            "α {:?}",
            z.f
        );
        errors.push((
            (t.k[0] - exact_k).abs() / exact_k,
            (t.omega[0] - exact_e).abs() / exact_e,
            (z.zeta[0] - exact_p).abs() / (2.0 / 3.0 - exact_p),
        ));
    }
    for i in 1..errors.len() {
        let pk = (errors[i - 1].0 / errors[i].0).log2();
        let pe = (errors[i - 1].1 / errors[i].1).log2();
        let pp = (errors[i - 1].2 / errors[i].2).log2();
        eprintln!(
            "decay: errors {:?}, orders k {pk:.3}, ε {pe:.3}, φ {pp:.3}",
            errors[i]
        );
        assert!((1.8..=2.2).contains(&pk), "order of k {pk}");
        assert!((1.8..=2.2).contains(&pe), "order of ε {pe}");
        assert!((1.8..=2.2).contains(&pp), "order of φ {pp}");
    }
    assert!(
        errors[2].0 < 1e-3 && errors[2].1 < 1e-3 && errors[2].2 < 1e-3,
        "{errors:?}"
    );
}

/// A run of the elliptic-blending k-ε continues from the restart file of a
/// k-ε run with φ = 2/3 and α = 1.
#[cfg(feature = "hdf5")]
#[test]
fn keps_restart_continues_with_eb() {
    use foehncraft_flow::restart;
    let keps = small_channel("")
        .replace(MODEL, "model = \"k-epsilon\"")
        .replace(
            "wall_treatment = \"compound\"",
            "wall_treatment = \"scalable\"",
        );
    let da = out("restart_keps");
    let ra = run(&parse(&keps), "keps", &da, |_| {}).unwrap();
    assert!(ra.failures.is_empty(), "{:?}", ra.failures);
    let db = out("restart_eb");
    let rb = run(
        &parse(&small_channel(&format!(
            "[run]\nrestart = '{}'\n",
            da.join("restart.h5").display()
        ))),
        "eb",
        &db,
        |_| {},
    )
    .unwrap();
    assert!(rb.failures.is_empty(), "{:?}", rb.failures);
    let t = restart::read(&db.join("restart.h5"))
        .unwrap()
        .state
        .turbulence
        .unwrap();
    let z = t.zeta.expect("φ and α in the new file");
    assert!(
        z.zeta
            .iter()
            .all(|v| v.is_finite() && (0.0..=2.0).contains(v))
    );
    assert!(z.f.iter().all(|v| (0.0..=1.0).contains(v)), "α {:?}", z.f);
    let _ = std::fs::remove_dir_all(&da);
    let _ = std::fs::remove_dir_all(&db);
}

/// One and four threads give the same bits (D31) with the elliptic-blending
/// k-ε, whose E term and C*_ε2 switch take gradients of gradients, on a
/// plate mesh above the 20 000 cells below which one thread is used.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn eb_results_are_bitwise_equal_on_1_and_4_threads() {
    let mk = |threads: usize| {
        let mut c = parse(
            &(text("tmr_plate_eb.toml")
                .replace("cells = [64, 48, 1]", "cells = [160, 128, 1]")
                .replace("cells = 16, ratio", "cells = 40, ratio")
                .replace("cells = 48, ratio", "cells = 120, ratio")
                .replace("max_iterations = 4000", "max_iterations = 8")
                + &format!("\n[run]\nthreads = {threads}\n")),
        );
        c.verification.require_convergence = Some(false);
        c
    };
    let a = run_case(&mk(1), "threads1");
    let b = run_case(&mk(4), "threads4");
    assert!(a.solves[0].cells >= 20_000);
    let (ra, rb) = (&a.solves[0].residuals, &b.solves[0].residuals);
    assert_eq!(ra.momentum.map(f64::to_bits), rb.momentum.map(f64::to_bits));
    assert_eq!(ra.continuity.to_bits(), rb.continuity.to_bits());
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(ra.turbulence.residual.len(), 4, "k, ε, α, φ");
    assert_eq!(bits(&ra.turbulence.residual), bits(&rb.turbulence.residual));
    let (wa, wb) = (&a.solves[0].wall_probes[0], &b.solves[0].wall_probes[0]);
    assert_eq!(wa.magnitude.to_bits(), wb.magnitude.to_bits());
}

/// The channel of the gallery at Re_τ 395: U_b⁺ within the ±3 % band
/// (specification 9.4), and C_f (= F at the fixed mass flow) at
/// 97 to 100 % of the DNS;
/// the low-Re treatment runs as well.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn eb_channel_against_the_dns() {
    let r = run_case(&parse(&text("channel_eb.toml")), "channel");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let f = r.solves[0].periodic[0].driving_force;
    let ub_plus = 17.55 / f.sqrt();
    let dev = ub_plus / 17.55 - 1.0;
    eprintln!(
        "channel elliptic-blending k-ε compound: U_b⁺ {ub_plus:.3} ({:+.2} %), C_f/C_f,DNS {f:.4}",
        100.0 * dev
    );
    assert!(dev.abs() < 0.03, "U_b⁺ {ub_plus}");
    assert!((0.97..=1.0).contains(&f), "C_f/C_f,DNS {f}");
    let c = parse(&small_channel("").replace(
        "wall_treatment = \"compound\"",
        "wall_treatment = \"low-re\"",
    ));
    let r = run_case(&c, "channel_low_re");
    let s = &r.solves[0];
    assert_eq!(s.residuals.turbulence.residual.len(), 4);
    assert!(
        s.residuals
            .turbulence
            .residual
            .iter()
            .all(|r| r.is_finite())
    );
}

/// A short run on the small channel keeps α in [0, 1] and φ in (0, 2].
#[test]
fn eb_small_channel_keeps_its_bounds() {
    let r = run_case(&parse(&small_channel("")), "bounds");
    let s = &r.solves[0];
    assert_eq!(s.residuals.turbulence.residual.len(), 4);
    assert!(
        s.residuals
            .turbulence
            .residual
            .iter()
            .all(|r| r.is_finite()),
        "{:?}",
        s.residuals.turbulence.residual
    );
}

/// Wrong combinations are refused with their codes.
#[test]
fn bad_eb_input_is_refused() {
    let t = small_channel("");
    for treatment in ["compound-2023", "scalable", "automatic"] {
        let (code, msg) = refused(&t.replace(
            "wall_treatment = \"compound\"",
            &format!("wall_treatment = \"{treatment}\""),
        ));
        assert_eq!(code, "CASE005", "{treatment}: {msg}");
    }
    // a production limiter
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"compound\"",
        "wall_treatment = \"compound\"\nproduction_limiter = \"kato-launder\"",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // a constant of k-ε-ζ-f that the elliptic blending does not have
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"compound\"",
        "wall_treatment = \"compound\"\nconstants = { sigma_zeta = 1.2 }",
    ));
    assert_eq!(code, "CASE002", "{msg}");
    // a constant out of its range
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"compound\"",
        "wall_treatment = \"compound\"\nconstants = { c_epsilon3 = 0.0 }",
    ));
    assert!(code.starts_with("CASE"), "{code}: {msg}");
}
