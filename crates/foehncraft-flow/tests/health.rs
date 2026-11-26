//! Stage 5 part 4b: the conservation ledger (D88, D89), the kinetic-energy
//! form of the second-law check (D90 part 4), the divergence report (D63)
//! and the restart geometry (D91); docs/verification/stage5.md.
#![cfg(feature = "hdf5")]

use foehncraft_flow::balance::BalanceSummary;
use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use foehncraft_flow::restart::{self, RestartGeometry};
use foehncraft_fv::disc::Soa3;
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("foehncraft_health_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

fn case(text: &str) -> FlowCase {
    FlowCase::from_toml(text, Path::new(".")).unwrap()
}

/// The lid-driven cavity at Re 100 on `n` x `n` hexahedra; `extra` adds
/// tables or keys at the end.
fn cavity(n: usize, extra: &str) -> FlowCase {
    case(&format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [{n}, {n}, 1]
[fluid]
viscosity = 0.01
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[output]
vtkhdf = false
wall_csv = false
{extra}
"#
    ))
}

fn taylor_green(n: usize, dt: f64) -> FlowCase {
    case(&format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [{n}, {n}, 1]
[fluid]
viscosity = 0.01
[boundary.xmin]
type = "velocity-inlet"
profile = "analytic"
[boundary.xmax]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymin]
type = "velocity-inlet"
profile = "analytic"
[boundary.ymax]
type = "velocity-inlet"
profile = "analytic"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[analytic]
type = "taylor-green"
[output]
vtkhdf = false
wall_csv = false
[time]
scheme = "bdf2"
time_step = {dt}
end_time = 0.2
"#
    ))
}

fn balance(r: &FlowReport) -> &BalanceSummary {
    r.solves.last().unwrap().balance.as_ref().unwrap()
}

/// D90 part 4: the lid power is the viscous dissipation plus the
/// numerical dissipation; the numerical dissipation is positive and its
/// share falls as the mesh is refined. The corner singularities of the
/// cavity keep its absolute value near constant (the lid power grows
/// with the mesh): the pressure part, the dissipation of the momentum
/// interpolation, stays positive and near constant there. The
/// Taylor–Green test below shows the fall of both parts on a smooth flow.
#[test]
fn cavity_energy_budget_on_three_meshes() {
    let mut rel = Vec::new();
    for n in [16, 32, 64] {
        let d = out(&format!("cavity{n}"));
        let r = run(&cavity(n, ""), "cavity", &d, |_| {}).unwrap();
        assert!(r.solves[0].converged);
        let b = &balance(&r).last;
        let k = &b.kinetic_energy;
        let lid = b.boundaries.iter().find(|p| p.boundary == "ymax").unwrap();
        assert!(lid.viscous_power > 0.0);
        // a closed domain: no mass flow, no convected energy, no pressure work
        assert_eq!(b.mass.imbalance, 0.0);
        let other: f64 = b
            .boundaries
            .iter()
            .filter(|p| p.boundary != "ymax")
            .map(|p| p.power().abs())
            .sum();
        assert!(other < 1e-12 * lid.viscous_power, "{other}");
        assert!(
            (lid.power() - k.viscous_dissipation - k.numerical_dissipation).abs()
                < 1e-12 * lid.power()
        );
        assert!(k.numerical_dissipation > 0.0, "{k:?}");
        assert!(k.relative > 0.0 && k.relative < 0.05, "{k:?}");
        assert!(
            (k.pressure_dissipation + k.scheme_dissipation - k.numerical_dissipation).abs()
                < 1e-12 * lid.power()
        );
        assert!(k.pressure_dissipation > 0.0, "{k:?}");
        rel.push(k.relative);
        assert!(d.join("balance.csv").exists() && d.join("balance.svg").exists());
        let _ = std::fs::remove_dir_all(&d);
    }
    assert!(rel[0] > rel[1] && rel[1] > rel[2], "{rel:?}");
}

/// The transient ledger on the Taylor–Green vortex: the mass ledger closes
/// to round-off, the kinetic energy that the vortex loses is the viscous
/// plus the numerical dissipation over the run, and both parts of the
/// numerical dissipation are positive and fall by more than a factor 3 per
/// halving of the mesh and the time step.
#[test]
fn taylor_green_ledger_falls_with_the_mesh() {
    let mut num = Vec::new();
    for n in [16usize, 32] {
        let d = out(&format!("tg{n}"));
        let r = run(&taylor_green(n, 0.32 / n as f64), "tg", &d, |_| {}).unwrap();
        let s = balance(&r);
        assert!(s.largest_mass_relative < 1e-12, "{s:?}");
        let run = s.run.as_ref().unwrap();
        let k = &run.kinetic_energy;
        assert!(k.storage < 0.0);
        assert!(
            (k.boundary - k.viscous_dissipation - k.numerical_dissipation - k.storage).abs()
                < 1e-12 * k.viscous_dissipation
        );
        assert!(k.numerical_dissipation > 0.0, "{k:?}");
        assert!(
            k.pressure_dissipation > 0.0 && k.scheme_dissipation > 0.0,
            "{k:?}"
        );
        assert_eq!(s.steps, Some(n * 5 / 8));
        num.push((
            k.numerical_dissipation,
            k.pressure_dissipation,
            k.scheme_dissipation,
        ));
        let csv = std::fs::read_to_string(d.join("balance.csv")).unwrap();
        assert!(csv.starts_with("# format foehncraft-balance 1"));
        assert_eq!(
            csv.lines().filter(|l| !l.starts_with('#')).count(),
            1 + n * 5 / 8
        );
        let _ = std::fs::remove_dir_all(&d);
    }
    assert!(num[0].0 > 3.0 * num[1].0, "{num:?}");
    assert!(num[0].1 > 3.0 * num[1].1, "{num:?}");
    assert!(num[0].2 > 3.0 * num[1].2, "{num:?}");
}

/// The balance tables during a steady run (`[balance] every`) evaluate
/// the state without changing it.
#[test]
fn balance_every_leaves_the_solution_unchanged() {
    let fields = |extra: &str, name: &str| {
        let d = out(name);
        let mut c = cavity(16, extra);
        c.run.restart_final = true;
        run(&c, "cavity", &d, |_| {}).unwrap();
        let r = restart::read(&d.join("restart.h5")).unwrap();
        let _ = std::fs::remove_dir_all(&d);
        (r.state.u.x, r.state.p, r.iteration)
    };
    assert_eq!(
        fields("", "every0"),
        fields("[balance]\nevery = 7", "every7")
    );
}

/// D63: a steady run with relaxation factors that make SIMPLE diverge
/// ends with RUN014, the last valid state, the snapshot and the report,
/// whose most likely cause is the relaxation with a fix that converges.
#[test]
fn steady_divergence_gives_the_report() {
    let d = out("diverge");
    let c = cavity(
        32,
        "[solver]\nalgorithm = \"simple\"\nmax_iterations = 200\n[solver.momentum]\nrelaxation = 0.99\n[solver.pressure]\nrelaxation = 1.0",
    );
    let mut c = c;
    c.fluid.viscosity = Some(0.001.into());
    let e = run(&c, "cavity", &d, |_| {}).unwrap_err();
    let dg = diag::of_error(&e);
    assert_eq!(dg.code, diag::RUN_DIVERGED, "{e:#}");
    let fix = dg.fix.unwrap().set;
    assert!(
        fix.iter()
            .any(|f| f.starts_with("solver.momentum.relaxation="))
    );
    for f in [
        "divergence.md",
        "divergence.json",
        "divergence.vtkhdf",
        "restart.h5",
        "flow.log",
    ] {
        assert!(d.join(f).exists(), "{f}");
    }
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(d.join("divergence.json")).unwrap()).unwrap();
    assert_eq!(j["format"], "foehncraft.divergence");
    assert_eq!(j["kind"], "steady");
    assert_eq!(j["causes"][0]["code"], diag::DIV_RELAXATION);
    assert!(j["failing_cells"].as_u64().unwrap() >= 1);
    let _ = std::fs::remove_dir_all(&d);
}

/// D63, transient: a time step far above the stable range; the most likely
/// cause is the time step, with a smaller one as the fix.
#[test]
fn transient_divergence_points_to_the_time_step() {
    let d = out("diverge_t");
    let c = cavity(
        32,
        "[solver]\ncorrectors = 1\n[time]\nscheme = \"crank-nicolson\"\ntime_step = 5.0\nsteps = 60",
    );
    let e = run(&c, "cavity", &d, |_| {}).unwrap_err();
    let dg = diag::of_error(&e);
    assert_eq!(dg.code, diag::RUN_DIVERGED, "{e:#}");
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(d.join("divergence.json")).unwrap()).unwrap();
    assert_eq!(j["kind"], "transient");
    assert_eq!(j["causes"][0]["code"], diag::DIV_TIME_STEP);
    assert!(
        j["causes"][0]["fix"][0]
            .as_str()
            .unwrap()
            .starts_with("time.time_step=")
    );
    // the restart file of the last valid state reads back
    let r = restart::read(&d.join("restart.h5")).unwrap();
    assert!(r.transient.is_some());
    let _ = std::fs::remove_dir_all(&d);
}

/// D91: the optional geometry datasets of the restart file round-trip.
#[test]
fn restart_geometry_round_trip() {
    let d = out("geometry");
    let mut c = cavity(8, "");
    c.run.restart_final = true;
    run(&c, "cavity", &d, |_| {}).unwrap();
    let mut r = restart::read(&d.join("restart.h5")).unwrap();
    assert!(r.geometry.is_none());
    let n = r.cell_id.len();
    let g = RestartGeometry {
        point_id: vec![3, 1, 2],
        points: Soa3 {
            x: vec![0.0, 1.0, 2.0],
            y: vec![0.5, 1.5, 2.5],
            z: vec![0.0; 3],
        },
        volumes: vec![vec![1.0; n], vec![0.5; n]],
        swept: vec![vec![0.25; r.face_id.len()]],
        volume_ratio: Some([0.5, 2.0]),
    };
    r.geometry = Some(g.clone());
    let f = d.join("geometry.h5");
    restart::write(&f, &r).unwrap();
    let back = restart::read(&f).unwrap().geometry.unwrap();
    assert_eq!(back.point_id, g.point_id);
    assert_eq!(back.points.y, g.points.y);
    assert_eq!(back.volumes, g.volumes);
    assert_eq!(back.swept, g.swept);
    assert_eq!(back.volume_ratio, g.volume_ratio);
    let _ = std::fs::remove_dir_all(&d);
}
