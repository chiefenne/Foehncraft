//! Stage 6, part 6d2: the k-ε-ζ-f₀ model of Popovac (2023) and the
//! compound wall treatment of the same paper (`compound-2023`), D132,
//! D134, docs/verification/stage6.md.
//!
//! The channel and the flat plate of the gallery (`channel_zeta_f0.toml`,
//! `tmr_plate_zeta_f0.toml`) and the backward-facing step
//! (`backstep_zeta_f0.toml`) are the validation cases.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_zeta_f0_{name}_{}", std::process::id()))
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
    text("channel_zeta_f0.toml")
        .replace("cells = [2, 48, 1]", "cells = [2, 24, 1]")
        .replace("ratio = 27.32", "ratio = 8.0")
        .replace("algorithm = \"coupled\"\n", "")
        .replace("max_iterations = 1000", "max_iterations = 30")
        + "\n[verification]\nrequire_convergence = false\n"
        + extra
}

/// The channel of the gallery: U_b⁺ from the driving force within the
/// ±3 % band (specification 9.4).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn zeta_f0_channel_meets_its_band() {
    let r = run_case(&parse(&text("channel_zeta_f0.toml")), "channel");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let f = r.solves[0].periodic[0].driving_force;
    let ub_plus = 17.55 / f.sqrt();
    let dev = ub_plus / 17.55 - 1.0;
    eprintln!(
        "channel k-ε-ζ-f₀: U_b⁺ {ub_plus:.3} ({:+.2} %)",
        100.0 * dev
    );
    assert!(dev.abs() < 0.03, "U_b⁺ {ub_plus}");
}

/// The flat plate from the TMR free stream: the boundary layer turns
/// turbulent on its own, every field converges (no positivity floor
/// holds k), and C_f(0.97) lies within 3 % of White's correlation
/// (specification 9.3).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn zeta_f0_plate_converges_from_the_free_stream() {
    let r = run_case(&parse(&text("tmr_plate_zeta_f0.toml")), "plate");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0];
    assert!(s.converged, "not converged");
    assert!(
        s.residuals.turbulence.floors.iter().all(|&n| n == 0),
        "{:?}",
        s.residuals.turbulence.floors
    );
    let cf = s.wall_probes[0].skin_friction.unwrap();
    let dev = cf / 2.875e-3 - 1.0;
    eprintln!(
        "plate k-ε-ζ-f₀: C_f {cf:.4e} ({:+.2} % from White)",
        100.0 * dev
    );
    assert!(dev.abs() < 0.03, "C_f {cf}");
}

/// One and four threads give the same bits (D31) with ζ-f₀ and the 2023
/// treatment (its non-equilibrium term reads the velocity and pressure
/// gradients of the near-wall cells), on a plate mesh above the 20 000
/// cells below which one thread is used.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn zeta_f0_results_are_bitwise_equal_on_1_and_4_threads() {
    let mk = |threads: usize| {
        let mut c = parse(
            &(text("tmr_plate_zeta_f0.toml")
                .replace("cells = [64, 48, 1]", "cells = [160, 128, 1]")
                .replace("cells = 16, ratio", "cells = 40, ratio")
                .replace("cells = 48, ratio", "cells = 120, ratio")
                .replace("max_iterations = 4000", "max_iterations = 8")
                + "\n[solver.turbulence]\npreconditioner = \"amg\"\n"
                + &format!("\n[run]\nthreads = {threads}\n")),
        );
        c.verification.require_convergence = Some(false);
        c
    };
    let a = run_case(&mk(1), "threads1");
    let b = run_case(&mk(4), "threads4");
    let (ra, rb) = (&a.solves[0].residuals, &b.solves[0].residuals);
    assert_eq!(ra.momentum.map(f64::to_bits), rb.momentum.map(f64::to_bits));
    assert_eq!(ra.continuity.to_bits(), rb.continuity.to_bits());
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(ra.turbulence.residual.len(), 4, "k, ε, f₀, ζ");
    assert_eq!(bits(&ra.turbulence.residual), bits(&rb.turbulence.residual));
    let (wa, wb) = (&a.solves[0].wall_probes[0], &b.solves[0].wall_probes[0]);
    assert_eq!(wa.magnitude.to_bits(), wb.magnitude.to_bits());
}

/// The 2023 treatment also serves ζ-f, and ζ-f₀ takes the 2007 treatment
/// and the integration to the wall; the runs stay finite.
#[test]
fn zeta_f0_treatment_combinations_run() {
    for t in [
        small_channel(""),
        small_channel("").replace(
            "wall_treatment = \"compound-2023\"",
            "wall_treatment = \"low-re\"",
        ),
        small_channel("").replace(
            "wall_treatment = \"compound-2023\"",
            "wall_treatment = \"compound\"",
        ),
        small_channel("").replace(
            "model = \"k-epsilon-zeta-f0\"",
            "model = \"k-epsilon-zeta-f\"",
        ),
    ] {
        let r = run_case(&parse(&t), "combination");
        let s = &r.solves[0];
        assert_eq!(s.residuals.turbulence.residual.len(), 4);
        assert!(
            s.residuals
                .turbulence
                .residual
                .iter()
                .all(|r| r.is_finite()),
            "{t}"
        );
    }
}

/// Wrong combinations are refused with their codes.
#[test]
fn bad_zeta_f0_input_is_refused() {
    let t = small_channel("");
    // a wall function of the k-ε family
    let (code, msg) = refused(&t.replace(
        "wall_treatment = \"compound-2023\"",
        "wall_treatment = \"scalable\"",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // the 2023 treatment with the high-Re k-ε model
    let (code, msg) = refused(&t.replace("model = \"k-epsilon-zeta-f0\"", "model = \"k-epsilon\""));
    assert_eq!(code, "CASE005", "{msg}");
    // N_ζ belongs to ζ-f₀
    let (code, msg) = refused(&t.replace(
        "model = \"k-epsilon-zeta-f0\"",
        "model = \"k-epsilon-zeta-f\"\nconstants = { n_zeta = 2.0 }",
    ));
    assert_eq!(code, "CASE002", "{msg}");
    // the sink of k is an option of the ζ-f models, not of EB (#141)
    let (code, msg) = refused(&t.replace(
        "model = \"k-epsilon-zeta-f0\"\nwall_treatment = \"compound-2023\"",
        "model = \"elliptic-blending-k-epsilon\"\nwall_treatment = \"compound\"\nk_sink = \"epsilon\"",
    ));
    assert_eq!(code, "CASE005", "{msg}");
    // C_ε1″/ζ^power with the power 1 or ½ only
    let (code, msg) = refused(&t.replace(
        "model = \"k-epsilon-zeta-f0\"",
        "model = \"k-epsilon-zeta-f0\"\nconstants = { c_epsilon1_power = 0.7 }",
    ));
    assert_eq!(code, "CASE004", "{msg}");
    let r = run_case(
        &parse(&t.replace(
            "model = \"k-epsilon-zeta-f0\"",
            "model = \"k-epsilon-zeta-f0\"\nk_sink = \"kolmogorov\"\n\
             constants = { n_zeta = 2.0, c_epsilon1_zeta = 0.012, c_epsilon1_power = 0.5 }",
        )),
        "constants",
    );
    assert!(r.failures.is_empty(), "{:?}", r.failures);
}

/// A jet started into gas at rest (a coarse quarter of a round jet from a
/// square nozzle with walls around it, BDF2): the extrapolated start of
/// k, ε and ζ may not fall below a tenth of the last value. Without that
/// bound the start of ζ went to its floor where the jet front arrives, f₀
/// turned large and the run diverged in step 11 (fine box) or earlier.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn zeta_f0_jet_starts_into_gas_at_rest() {
    let t = r#"case_format = "0.7"
[mesh.generate]
type = "hex"
dimension = 3
lower = [0.0, 0.0, 0.0]
upper = [30.0, 10.0, 10.0]
cells = [40, 18, 18]
segments = { x = [ { to = 5.0, cells = 15, ratio = 2.0 }, { to = 30.0, cells = 25, ratio = 3.0 } ], y = [ { to = 0.443, cells = 3 }, { to = 10.0, cells = 15, ratio = 20.0 } ], z = [ { to = 0.443, cells = 3 }, { to = 10.0, cells = 15, ratio = 20.0 } ] }
split = [ { boundary = "xmin", axis = "y", at = 0.443, name = "wall_y" }, { boundary = "xmin", axis = "z", at = 0.443, name = "wall_z" } ]
[fluid]
density = 1.0
viscosity = 1.0e-5
pressure_values = "gauge"
[turbulence]
model = "k-epsilon-zeta-f0"
wall_treatment = "compound-2023"
[initial]
velocity = [0.0, 0.0, 0.0]
turbulence_intensity = 0.05
viscosity_ratio = 10.0
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
turbulence_intensity = 0.05
viscosity_ratio = 10.0
[boundary.wall_y]
type = "wall"
[boundary.wall_z]
type = "wall"
[boundary.ymin]
type = "symmetry"
[boundary.zmin]
type = "symmetry"
[boundary.ymax]
type = "total-pressure-inlet"
total_pressure = 0.0
turbulence_intensity = 0.01
viscosity_ratio = 1.0
[boundary.zmax]
type = "total-pressure-inlet"
total_pressure = 0.0
turbulence_intensity = 0.01
viscosity_ratio = 1.0
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[time]
scheme = "bdf2"
time_step = 0.2
end_time = 6.0
[solver]
tolerance = 1e-6
max_iterations = 50
outer_iterations = 3
corrector = "simplec"
[output]
vtkhdf = false
"#;
    let dir = out("jet_at_rest");
    let r = run(&parse(t), "jet_at_rest", &dir, |_| {});
    let _ = std::fs::remove_dir_all(&dir);
    let r = r.unwrap_or_else(|e| panic!("{e:#}"));
    assert!(r.failures.is_empty(), "{:?}", r.failures);
}
