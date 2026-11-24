//! Two temperatures in porous regions (D271, porous plan stage PM2): the
//! gas and the matrix exchange heat through h a, against exact solutions
//! (checks V3 to V6 of the porous plan).

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::porous_energy::PorousHeat;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn run_case(text: &str, name: &str) -> FlowReport {
    let dir = std::env::temp_dir().join(format!(
        "foehncraft_porous_energy_{name}_{}",
        std::process::id()
    ));
    let case = FlowCase::from_toml(text, &cases()).unwrap();
    let rep = run(&case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn porous(rep: &FlowReport) -> PorousHeat {
    let s = &rep.solves[0];
    let e = s
        .scalars
        .iter()
        .find_map(|x| x.energy.as_ref())
        .expect("an energy balance");
    e.porous[0].clone()
}

/// A closed box filled with the matrix, gas at rest at 400 K, matrix at
/// 300 K.
fn box_case(scheme: &str, dt: f64, end: f64) -> String {
    format!(
        r#"
case_format = "0.6"

[mesh.generate]
type = "hex"
dimension = 2
upper = [0.1, 0.1, 0.01]
cells = [4, 4, 1]
zones = [{{ name = "bed", lower = [0.0, 0.0, 0.0], upper = [0.1, 0.1, 0.01] }}]

[materials.ceramic]
kind = "solid"
description = "a ceramic matrix (chosen values)"
density = 2300.0
specific_heat = 900.0
conductivity = 2.0

[fluid]
density = 1.2
viscosity = 1.8e-5
specific_heat = 1005.0
conductivity = 0.026

[regions.bed]
kind = "porous"
zones = ["bed"]
[regions.bed.porous]
viscous = 1.0e6
porosity = 0.5
solid = "ceramic"
heat_transfer_coefficient = 20.0
specific_area = 150.0
solid_initial = 300.0
source = "chosen values for the exact check"

[energy]
initial = 400.0

[boundary.xmin]
type = "wall"
[boundary.xmin.energy]
heat_flux = 0.0
[boundary.xmax]
type = "wall"
[boundary.xmax.energy]
heat_flux = 0.0
[boundary.ymin]
type = "wall"
[boundary.ymin.energy]
heat_flux = 0.0
[boundary.ymax]
type = "wall"
[boundary.ymax.energy]
heat_flux = 0.0
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"

[solver.energy]
tolerance = 1e-13

[time]
scheme = "{scheme}"
time_step = {dt}
end_time = {end}

[solver]
outer_iterations = 3

[verification]
require_convergence = false

[output]
wall_csv = false
"#
    )
}

const CG: f64 = 0.5 * 1.2 * 1005.0;
const CS: f64 = 0.5 * 2300.0 * 900.0;
const HA: f64 = 20.0 * 150.0;

/// V3: no flow, no gradients: the heat of gas and matrix is kept, and the
/// difference T_g − T_s falls by 1/(1 + λ Δt) in each implicit Euler step,
/// λ = h a (1/C_g + 1/C_s); with small steps that is exp(−λ t).
#[test]
fn gas_and_matrix_relax_to_each_other_exactly() {
    let lambda = HA * (1.0 / CG + 1.0 / CS);
    let (dt, steps) = (0.02, 10);
    let p = porous(&run_case(&box_case("bdf1", dt, dt * steps as f64), "relax"));
    let mean0 = (CG * 400.0 + CS * 300.0) / (CG + CS);
    let mean = (CG * p.gas_mean + CS * p.matrix_mean) / (CG + CS);
    assert!(
        (mean - mean0).abs() < 1e-9 * mean0,
        "heat not kept: {mean} against {mean0}"
    );
    let d = p.gas_mean - p.matrix_mean;
    let exact = 100.0 / (1.0 + lambda * dt).powi(steps);
    assert!((d - exact).abs() < 1e-8 * 100.0, "{d} against {exact}");
    assert!(p.matrix_max - p.matrix_min < 1e-9, "{p:?}");
}

/// e^{−z} I₀(z) (the series; z ≤ 40 here).
fn i0e(z: f64) -> f64 {
    let (mut term, mut sum) = (1.0, 1.0);
    let q = 0.25 * z * z;
    for k in 1..400 {
        term *= q / (k * k) as f64;
        sum += term;
        if term < 1e-17 * sum {
            break;
        }
    }
    sum * (-z).exp()
}

/// Schumann's J(x, y) = 1 − ∫₀ˣ e^{−y−s} I₀(2√(y s)) ds (Simpson).
fn schumann_j(x: f64, y: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    let n = 4000;
    let h = x / n as f64;
    let f = |s: f64| {
        let z = 2.0 * (y * s).sqrt();
        (-(y + s) + z).exp() * i0e(z)
    };
    let mut sum = f(0.0) + f(x);
    for i in 1..n {
        sum += if i % 2 == 1 { 4.0 } else { 2.0 } * f(i as f64 * h);
    }
    1.0 - sum * h / 3.0
}

/// The gas and matrix temperature of Schumann's bed (0 before the step, 1
/// at the inlet after it) at ξ = h a x/(ρ c u), η = h a (t − ε x/u)/C_s.
fn schumann(xi: f64, eta: f64) -> (f64, f64) {
    if eta < 0.0 {
        return (0.0, 0.0);
    }
    (schumann_j(xi, eta), 1.0 - schumann_j(eta, xi))
}

#[test]
fn schumann_function_has_its_limits() {
    // the gas at the front falls as e^{−ξ}; the matrix at the inlet warms
    // as 1 − e^{−η}; J(x, y) + J(y, x) = 1 + e^{−x−y} I₀(2√(x y))
    assert!((schumann(2.0, 0.0).0 - (-2.0f64).exp()).abs() < 1e-10);
    assert!((schumann(0.0, 1.5).1 - (1.0 - (-1.5f64).exp())).abs() < 1e-10);
    let (x, y) = (1.3f64, 2.1f64);
    let z = 2.0 * (x * y).sqrt();
    let id = 1.0 + (-x - y + z).exp() * i0e(z);
    assert!((schumann_j(x, y) + schumann_j(y, x) - id).abs() < 1e-10);
}

/// `porous_schumann.toml`: a bed of length 1 m without conduction, air at
/// 0.1 m/s, its inlet temperature stepped from 300 K to 400 K at t = 0.
fn schumann_case(cells: usize, dt: f64, end: f64) -> String {
    schumann_case_h(cells, dt, end, 10.0)
}

fn schumann_case_h(cells: usize, dt: f64, end: f64, h: f64) -> String {
    let dx = 1.0 / cells as f64;
    let text = std::fs::read_to_string(cases().join("porous_schumann.toml"))
        .unwrap()
        .replace("\r\n", "\n");
    let swap = |t: String, from: &str, to: String| {
        assert!(t.contains(from), "{from}");
        t.replace(from, &to)
    };
    let t = swap(
        text,
        "cells = [100, 1, 1]",
        format!("cells = [{cells}, 1, 1]"),
    );
    let t = swap(
        t,
        "heat_transfer_coefficient = 10.0",
        format!("heat_transfer_coefficient = {h:?}"),
    );
    let t = swap(t, "time_step = 15.0", format!("time_step = {dt:?}"));
    let t = swap(t, "end_time = 3000.0", format!("end_time = {end:?}"));
    let t = swap(
        t,
        "fields_times = [1000.0, 2000.0]\n",
        "vtkhdf = false\n".to_string(),
    );
    swap(
        t,
        "from = [0.005, 0.05, 0.005], to = [0.995, 0.05, 0.005], points = 100",
        format!(
            "from = [{:?}, 0.05, 0.005], to = [{:?}, 0.05, 0.005], points = {cells}",
            0.5 * dx,
            1.0 - 0.5 * dx
        ),
    )
}

/// The cell values along the bed: (x, gas, matrix).
fn schumann_run(cells: usize, dt: f64, end: f64) -> Vec<(f64, f64, f64)> {
    line_run(&schumann_case(cells, dt, end), &format!("schumann{cells}"))
}

/// V6: the heat balance of every step closes (the storage of gas and
/// matrix is what the inlet brings minus what the outlet takes) to round-off.
fn line_run(text: &str, name: &str) -> Vec<(f64, f64, f64)> {
    let dir = std::env::temp_dir().join(format!(
        "foehncraft_porous_energy_{name}_{}",
        std::process::id()
    ));
    let case = FlowCase::from_toml(text, &cases()).unwrap();
    let rep = run(&case, name, &dir, |_| {}).unwrap();
    let h = rep.solves[0].transient.as_ref().unwrap();
    assert!(h.heat_relative_max.unwrap() < 1e-9, "heat balance {h:?}");
    let e = rep.solves[0]
        .scalars
        .iter()
        .find_map(|x| x.energy.as_ref())
        .unwrap();
    assert!(
        e.storage_regions
            .iter()
            .any(|r| r.region == "bed matrix" && r.storage > 0.0),
        "{e:?}"
    );
    let text = std::fs::read_to_string(dir.join("line_axis.csv")).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let (cx, ct, cs) = (col("x"), col("temperature"), col("matrix_temperature"));
    lines
        .map(|l| {
            let v: Vec<f64> = l.split(',').map(|x| x.parse().unwrap()).collect();
            (v[cx], v[ct], v[cs])
        })
        .collect()
}

/// V4 (Schumann 1929): the gas and matrix temperatures of a bed without
/// conduction after a step of the inlet temperature, exact in x and t;
/// the error of the cell values falls with the square of the cell size.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn a_bed_without_conduction_follows_schumann() {
    let (rho_c_u, cs, ha, eps, u) = (1.2 * 1005.0 * 0.1, 0.6 * 2300.0 * 900.0, 603.0, 0.4, 0.1);
    let end = 3000.0;
    let mut errs = Vec::new();
    for cells in [25, 50, 100] {
        let prof = schumann_run(cells, 15.0, end);
        assert_eq!(prof.len(), cells);
        let (mut eg, mut es) = (0.0, 0.0);
        for &(x, tg, ts) in &prof {
            let (g, s) = schumann(ha * x / rho_c_u, ha * (end - eps * x / u) / cs);
            eg += ((tg - 300.0) / 100.0 - g).abs() / cells as f64;
            es += ((ts - 300.0) / 100.0 - s).abs() / cells as f64;
        }
        eprintln!("{cells} cells: L1 error gas {eg:.3e}, matrix {es:.3e}");
        errs.push((eg, es));
    }
    // 5.2e-5 and 4.7e-5 of the step on 100 cells, orders 1.96 to 2.01
    assert!(errs[2].0 < 1e-4 && errs[2].1 < 1e-4, "{errs:?}");
    for w in errs.windows(2) {
        let (og, os) = ((w[0].0 / w[1].0).log2(), (w[0].1 / w[1].1).log2());
        assert!(
            og > 1.8 && os > 1.8,
            "observed orders {og}, {os} ({errs:?})"
        );
    }
}

/// V5: with h a → large the gas and the matrix stay together and the
/// front moves with the mixture's speed ρ_g c_g u / (ε ρ_g c_g + (1 − ε)
/// ρ_s c_s): its middle (half the step) lies where the heat brought in puts
/// it.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn a_bed_with_a_large_exchange_is_in_equilibrium() {
    let end = 3000.0;
    let prof = line_run(&schumann_case_h(100, 15.0, end, 1.0e5), "equilibrium");
    let dmax = prof.iter().map(|p| (p.1 - p.2).abs()).fold(0.0, f64::max);
    assert!(dmax < 0.1, "gas and matrix apart by {dmax} K");
    let cg = 0.4 * 1.2 * 1005.0;
    let speed = 1.2 * 1005.0 * 0.1 / (cg + 0.6 * 2300.0 * 900.0);
    // the middle of the front: where the temperature crosses 350 K
    let i = prof.iter().position(|p| p.2 < 350.0).unwrap();
    let (a, b) = (prof[i - 1], prof[i]);
    let x = a.0 + (350.0 - a.2) / (b.2 - a.2) * (b.0 - a.0);
    assert!(
        (x - speed * end).abs() < 0.01,
        "front at {x} m, {} m expected",
        speed * end
    );
}
