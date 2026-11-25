//! The steady spray mode (D260, stage SP2s; checks of
//! `docs/requirements/step28-sprays.md` section 5.9).
//!
//! - V28.10: the steady mode gives the same bits on 1 and 4 threads (the
//!   stream pass is a pure function, the records sorted per cell), with and
//!   without the eddy-lifetime dispersion.
//! - The books: the ledger pairs of the final pass are zero to the bit, the
//!   parcel accounts and the stream mass book close to round-off, the gas
//!   thermal account closes with the spray term; the run is converged only
//!   with the source residual below its tolerance (D260 (4)).
//! - One stream in a uniform gas: its residence time against the exact
//!   trajectory of the drag law (RK4 with small steps).
//! - The refusals of D260 (6).

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_kernels::particle::drag_correction;
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "foehncraft_steady_spray_{name}_{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&d);
    d
}

fn examples() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

/// The steady spray-cooling duct, coarser and with fewer streams.
fn duct(iterations: usize) -> String {
    std::fs::read_to_string(examples().join("spray_cooling_steady.toml"))
        .unwrap()
        .replace("cells = [100, 20, 1]", "cells = [40, 8, 1]")
        .replace("streams = 400", "streams = 60")
        .replace(
            "max_iterations = 2000",
            &format!("max_iterations = {iterations}"),
        )
}

fn rows(d: &Path, file: &str) -> String {
    std::fs::read_to_string(d.join(file))
        .unwrap()
        .lines()
        .skip(1)
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn the_steady_duct_converges_with_closed_books() {
    let case = FlowCase::from_toml(&duct(2000), &examples()).unwrap();
    let d = out("books");
    let r = run(&case, "books", &d, |_| {}).unwrap();
    let s = r.solves.last().unwrap();
    assert!(s.converged, "the steady spray duct did not converge");
    let sp = s.spray_steady.as_ref().expect("steady spray summary");
    assert!(sp.two_way && sp.converged);
    assert!(sp.source_residual.max() <= sp.keys.tolerance);
    assert!(sp.passes > 2);
    let b = &sp.books;
    assert_eq!(b.streams, 60);
    assert_eq!(b.lost_streams + b.unfinished_streams, 0);
    // the mass book of the streams
    let left: f64 = b.escaped_mass.iter().sum::<f64>() + b.deposited_mass.iter().sum::<f64>();
    assert!(
        (b.injected_mass - 0.01).abs() < 1e-15,
        "{}",
        b.injected_mass
    );
    assert!((left - b.injected_mass).abs() < 1e-15);
    // the momentum and enthalpy books of the streams
    for k in 0..3 {
        assert!(
            sp.momentum_imbalance[k].abs() < 1e-14,
            "{:?}",
            sp.momentum_imbalance
        );
    }
    assert!(
        sp.heat_imbalance.abs() < 1e-12 * b.enthalpy[0],
        "{}",
        sp.heat_imbalance
    );
    // the droplets take heat from the air
    assert!(sp.pass.gas_heat < -1000.0, "{}", sp.pass.gas_heat);
    let bal = &s.balance.as_ref().unwrap().last;
    for id in [2801, 2802, 2803, 2804] {
        let p = bal.pairs.iter().find(|p| p.pair == id).expect("pair");
        assert_eq!(p.residual, 0.0, "pair {id}");
    }
    let acc = |q: &str| bal.ledger.iter().find(|a| a.quantity == q).unwrap();
    for q in [
        "parcel momentum x",
        "parcel momentum y",
        "parcel momentum z",
        "parcel enthalpy",
    ] {
        assert!(acc(q).relative.abs() < 1e-12, "{q}: {}", acc(q).relative);
    }
    let th = acc("thermal energy");
    assert!(th.relative.abs() < 1e-4, "thermal account {}", th.relative);
    let heat = th
        .terms
        .iter()
        .find(|t| t.name == "spray heat")
        .unwrap()
        .value;
    assert_eq!(heat, sp.pass.gas_heat);
    // spray.csv: one row per pass, the final one with relaxation 0
    let csv = rows(&d, "spray.csv");
    let n = csv.lines().filter(|l| !l.starts_with('#')).count() - 1;
    assert_eq!(n as u32, sp.passes);
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_steady_mode_gives_the_same_bits_on_any_number_of_threads() {
    // V28.10, with the eddy-lifetime dispersion (keyed by the stream)
    let plain = duct(60);
    let dispersed = plain.replace(
        "wall_mode = \"stick\"",
        "wall_mode = \"stick\"\ndispersion = \"eddy-lifetime\"",
    );
    for (k, text) in [plain, dispersed].iter().enumerate() {
        let text = text.replace(
            "max_iterations = 60",
            "max_iterations = 60\n[verification]\nrequire_convergence = false",
        );
        let case = FlowCase::from_toml(&text, &examples()).unwrap();
        let runs: Vec<(String, String)> = [1, 4]
            .iter()
            .map(|&n| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(n)
                    .build()
                    .unwrap()
                    .install(|| {
                        let name = format!("threads{k}_{n}");
                        let d = out(&name);
                        run(&case, &name, &d, |_| {}).unwrap();
                        let r = (rows(&d, "spray.csv"), rows(&d, "residuals.csv"));
                        let _ = std::fs::remove_dir_all(&d);
                        r
                    })
            })
            .collect();
        assert_eq!(runs[0].0, runs[1].0, "spray.csv differs");
        assert_eq!(runs[0].1, runs[1].1, "residuals.csv differs");
        assert!(runs[0].0.lines().count() > 6);
    }
}

#[test]
fn another_seed_draws_other_streams() {
    let text = duct(30).replace(
        "max_iterations = 30",
        "max_iterations = 30\n[verification]\nrequire_convergence = false",
    );
    let go = |text: &str, name: &str| {
        let case = FlowCase::from_toml(text, &examples()).unwrap();
        let d = out(name);
        let r = run(&case, name, &d, |_| {}).unwrap();
        let _ = std::fs::remove_dir_all(&d);
        r.solves.last().unwrap().spray_steady.clone().unwrap()
    };
    let a = go(&text, "seed0");
    let b = go(
        &text.replace("max_time = 5.0", "max_time = 5.0\nseed = 7"),
        "seed7",
    );
    assert_ne!(a.pass.gas_heat, b.pass.gas_heat);
    // the same physics: the heat within a few per cent
    let rel = (a.pass.gas_heat - b.pass.gas_heat).abs() / a.pass.gas_heat.abs();
    assert!(rel < 0.1, "{rel}");
}

/// A 1 m channel with gas at U = 2 m/s (slip walls, laminar), one way.
fn channel(streams: u32) -> String {
    format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 0.1, 0.1]
cells = [10, 2, 1]
[fluid]
density = 1.2
viscosity = 1.8e-5
[initial]
velocity = [2.0, 0.0, 0.0]
[boundary.xmin]
type = "velocity-inlet"
velocity = [2.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "slip-wall"
[boundary.ymax]
type = "slip-wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[spray]
max_courant = 0.01
[spray.liquid]
density = 998.0
[[spray.injector]]
name = "jet"
position = [0.1, 0.05, 0.05]
direction = [1.0, 0.0, 0.0]
speed = 0.0
mass_flow = 1e-4
diameter = 50e-6
[spray.steady]
streams = {streams}
[output]
vtkhdf = false
wall_csv = false
[solver]
tolerance = 1e-10
max_iterations = 200
"#
    )
}

/// The time a droplet of 50 µm released at rest at x = 0.1 needs to reach
/// x = 1 in the uniform gas: dv/dt = (U − v) f(Re) / τ_St, RK4.
fn residence_exact() -> f64 {
    let (u, rho, mu, rho_p, d) = (2.0, 1.2, 1.8e-5, 998.0, 50e-6);
    let tau = rho_p * d * d / (18.0 * mu);
    let acc = |v: f64| {
        let re = rho * d * (u - v).abs() / mu;
        (u - v) * drag_correction(re) / tau
    };
    let (mut t, mut x, mut v) = (0.0, 0.1, 0.0);
    let h = 1e-6;
    loop {
        let k1 = (v, acc(v));
        let k2 = (v + 0.5 * h * k1.1, acc(v + 0.5 * h * k1.1));
        let k3 = (v + 0.5 * h * k2.1, acc(v + 0.5 * h * k2.1));
        let k4 = (v + h * k3.1, acc(v + h * k3.1));
        let xn = x + h / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
        let vn = v + h / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
        if xn >= 1.0 {
            // the last part of the step, linear in x
            return t + h * (1.0 - x) / (xn - x);
        }
        (t, x, v) = (t + h, xn, vn);
    }
}

#[test]
fn a_stream_crosses_a_uniform_gas_on_its_exact_path() {
    let case = FlowCase::from_toml(&channel(1), Path::new(".")).unwrap();
    let d = out("channel");
    let r = run(&case, "channel", &d, |_| {}).unwrap();
    let s = r.solves.last().unwrap();
    let sp = s.spray_steady.as_ref().unwrap();
    assert!(!sp.two_way && sp.passes == 1);
    let b = &sp.books;
    assert_eq!(b.streams, 1);
    let outlet = sp.boundaries.iter().find(|x| x.boundary == "xmax").unwrap();
    assert_eq!(outlet.escaped, 1e-4);
    let exact = residence_exact();
    let rel = (b.residence_max - exact).abs() / exact;
    eprintln!("residence {} against {exact} ({rel:e})", b.residence_max);
    assert!(
        rel < 1e-3,
        "residence {} against {exact} ({rel:e})",
        b.residence_max
    );
    // one way: no source in the gas; the drag on the stream is what the
    // books leave, ṁ times the exit speed (the gas speed, nearly)
    assert_eq!(sp.pass.gas_heat, 0.0);
    let f = b.momentum.from_gas[0];
    assert!((f - 1e-4 * 2.0).abs() < 1e-3 * 2e-4, "{f}");
    if cfg!(feature = "hdf5") {
        assert!(d.join("spray").join("streams.vtkhdf").exists());
    }
    let _ = std::fs::remove_dir_all(&d);
}

fn refused(text: &str) -> String {
    match FlowCase::from_toml(text, Path::new(".")) {
        Err(e) => format!("{e:#}"),
        Ok(case) => format!(
            "{:#}",
            run(&case, "refused", &out("refused"), |_| {}).unwrap_err()
        ),
    }
}

#[test]
fn the_steady_mode_refuses_what_has_no_steady_state() {
    // D260 (6): [spray.steady] with [time]
    let e = refused(&channel(1).replace(
        "[output]",
        "[time]\ntime_step = 0.01\nend_time = 0.1\n[output]",
    ));
    assert!(
        e.contains("SPR002") || e.contains("[spray.steady] is the steady spray mode"),
        "{e}"
    );
    // an injector that ends
    let e = refused(&channel(1).replace("diameter = 50e-6", "diameter = 50e-6\nend = 1.0"));
    assert!(e.contains("start and end have no meaning"), "{e}");
    // a transient case still needs parcels per second
    let e = refused(
        &channel(1)
            .replace("[spray.steady]\nstreams = 1\n", "")
            .replace(
                "[output]",
                "[time]\ntime_step = 0.01\nend_time = 0.1\n[output]",
            ),
    );
    assert!(e.contains("parcels_per_second must be positive"), "{e}");
    // two-way with the coupled algorithm
    let e = refused(
        &channel(1)
            .replace(
                "[spray.liquid]",
                "[spray.coupling]\nmode = \"two-way\"\n[spray.liquid]",
            )
            .replace(
                "max_iterations = 200",
                "max_iterations = 200\nalgorithm = \"coupled\"",
            ),
    );
    assert!(e.contains("needs SIMPLE or SIMPLEC"), "{e}");
}
