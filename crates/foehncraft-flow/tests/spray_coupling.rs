//! Two-way coupling of the sprays (SP2, D247 to D252, checks V28.1 to
//! V28.4 of `docs/requirements/step28-sprays.md`).
//!
//! The test bed is one cell of gas with slip walls (and a pressure outlet
//! across the motion, which fixes the pressure level) and a cloud of
//! parcels shot into it at its centre along x: the gas of the cell is a
//! lumped mass (no shear, no flux, the least-squares gradient is zero by
//! symmetry), so gas and parcels relax towards one common velocity with
//! the mixture momentum conserved. The reference is the same two-mass
//! relaxation integrated by RK4 with small steps.
//!
//! The one thing the lumped picture leaves out: a slip wall keeps the
//! normal velocity zero, so the viscous stress of the normal component on
//! the two end walls slows the gas at the rate 4 μ / (ρ L²). With the cell
//! 10 m long (and 1 m³) that is 6·10⁻⁷ per second; the checks of the
//! mixture momentum allow for it, the checks of the exchange books do not
//! (they are exact).

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_kernels::particle::{drag_correction, nusselt};
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("foehncraft_coupling_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

const RHO_G: f64 = 1.2;
const MU: f64 = 1.8e-5;
const RHO_L: f64 = 998.0;
const D: f64 = 50e-6;
const V0: f64 = 1.0;

/// The cell case: gas at rest in 1 m³, liquid of mass `z` times the gas
/// mass shot along x at V0 in the first 1e-6 s.
fn case(z: f64, dt: f64, end: f64, passes: &str) -> String {
    let m_liquid = z * RHO_G;
    let duration = 1e-6;
    format!(
        r#"
[mesh.generate]
type = "hex"
dimension = 3
upper = [10.0, 0.316227766016838, 0.316227766016838]
cells = [1, 1, 1]
[fluid]
density = {RHO_G}
viscosity = {MU}
[initial]
velocity = [0.0, 0.0, 0.0]
[boundary.xmin]
type = "slip-wall"
[boundary.xmax]
type = "slip-wall"
[boundary.ymin]
type = "slip-wall"
[boundary.ymax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.zmin]
type = "slip-wall"
[boundary.zmax]
type = "slip-wall"
[time]
scheme = "bdf1"
time_step = {dt}
end_time = {end}
[spray]
[spray.coupling]
mode = "two-way"
passes = {passes}
tolerance = 1e-10
max_passes = 30
[spray.liquid]
density = {RHO_L}
[[spray.injector]]
name = "cloud"
position = [5.0, 0.158113883008419, 0.158113883008419]
direction = [1.0, 0.0, 0.0]
speed = {V0}
mass_flow = {mf}
diameter = {D}
parcels_per_second = 1e8
start = 0.0
end = {duration}
[output]
vtkhdf = false
wall_csv = false
[solver]
tolerance = 1e-12
max_iterations = 200
[verification]
require_convergence = false
"#,
        mf = m_liquid / duration
    )
}

struct Run {
    /// Gas velocity x at the end.
    u: f64,
    /// Momentum the gas received from the parcels over the run.
    gas_impulse: f64,
    /// Parcel momentum x and mass at the end.
    p_parcels: f64,
    m_parcels: f64,
    /// The rows of spray.csv: (gas force x).
    force: Vec<f64>,
    csv: String,
    pair_residual: f64,
    parcel_imbalance: f64,
    injected: f64,
}

fn go(name: &str, text: &str) -> Run {
    let d = out(name);
    let case = FlowCase::from_toml(text, Path::new(".")).unwrap();
    let r = run(&case, name, &d, |_| {}).unwrap();
    let s = r.solves.last().unwrap();
    let bal = s.balance.as_ref().expect("balance");
    let run = bal.run.as_ref().expect("run ledger");
    let acc = |q: &str| run.ledger.iter().find(|a| a.quantity == q).unwrap();
    let gas = acc("momentum x");
    let par = acc("parcel momentum x");
    let csv = std::fs::read_to_string(d.join("spray.csv")).unwrap();
    let mut lines = csv.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = head.iter().position(|h| *h == "gas_force_x").unwrap();
    let force = lines
        .map(|l| l.split(',').nth(col).unwrap().parse::<f64>().unwrap())
        .collect();
    let sp = s.spray.as_ref().unwrap();
    let injected = par
        .terms
        .iter()
        .find(|t| t.name == "injected")
        .map_or(0.0, |t| t.value);
    let pair = run
        .pairs
        .iter()
        .find(|p| p.pair == 2801)
        .map_or(f64::NAN, |p| p.residual);
    let out = Run {
        gas_impulse: sp.coupling.as_ref().unwrap().gas_impulse[0],
        u: gas.content.unwrap() / RHO_G,
        p_parcels: par.content.unwrap(),
        m_parcels: sp.in_flow_kg,
        force,
        csv: csv
            .lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n"),
        pair_residual: pair,
        parcel_imbalance: par.imbalance,
        injected,
    };
    let _ = std::fs::remove_dir_all(&d);
    out
}

/// The lumped two-mass reference: gas mass M at rest, liquid mass m at V0,
/// drag of Schiller and Naumann; RK4 to time t.
fn reference(z: f64, t: f64) -> (f64, f64) {
    let (big, m) = (RHO_G, z * RHO_G);
    let (mut u, mut v) = (0.0f64, V0);
    let n = 200_000;
    let h = t / n as f64;
    let rhs = |u: f64, v: f64| {
        let re = RHO_G * D * (u - v).abs() / MU;
        let tau = RHO_L * D * D / (18.0 * MU * drag_correction(re));
        let f = m * (v - u) / tau;
        (f / big, -f / m)
    };
    for _ in 0..n {
        let (a1, b1) = rhs(u, v);
        let (a2, b2) = rhs(u + 0.5 * h * a1, v + 0.5 * h * b1);
        let (a3, b3) = rhs(u + 0.5 * h * a2, v + 0.5 * h * b2);
        let (a4, b4) = rhs(u + h * a3, v + h * b3);
        u += h / 6.0 * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
        v += h / 6.0 * (b1 + 2.0 * b2 + 2.0 * b3 + b4);
    }
    (u, v)
}

#[test]
fn a_cloud_and_its_gas_reach_one_velocity_at_any_step() {
    // V28.2, V28.3: loadings 0.1 to 10, steps of 0.1 to 100 response times
    // (τ_p ≈ 7.7 ms for 50 µm water in air)
    for (z, dt) in [(0.1f64, 1e-3f64), (1.0, 0.05), (10.0, 0.5), (1.0, 0.77)] {
        let end = 40.0 * 7.7e-3 * (1.0 + 1.0 / z).max(1.0) + 3.0 * dt;
        let end = (end / dt).ceil() * dt;
        let r = go("cloud", &case(z, dt, end, "\"converged\""));
        let p0 = z * RHO_G * V0;
        // the mixture keeps the momentum it was given
        assert!(
            (r.injected - p0).abs() < 1e-12 * p0,
            "injected {} against {p0}",
            r.injected
        );
        // the exchange books: what the gas received plus what the parcels
        // hold is what was injected, to round-off
        assert!(
            (r.gas_impulse + r.p_parcels - p0).abs() < 1e-12 * p0,
            "Z {z} dt {dt}: books {} + {}",
            r.gas_impulse,
            r.p_parcels
        );
        let total = RHO_G * r.u + r.p_parcels;
        assert!(
            (total - p0).abs() < 1e-5 * p0,
            "Z {z} dt {dt}: gas {} + parcels {} = {total} against {p0}",
            RHO_G * r.u,
            r.p_parcels
        );
        // one common velocity: the equilibrium of the mixture
        let ue = p0 / (RHO_G * (1.0 + z));
        let v = r.p_parcels / r.m_parcels;
        assert!(
            (r.u - ue).abs() < 1e-5 * ue,
            "Z {z} dt {dt}: gas {} against {ue}",
            r.u
        );
        assert!(
            (v - r.u).abs() < 1e-6 * ue,
            "Z {z} dt {dt}: parcels {v} against {}",
            r.u
        );
        // no overshoot: the gas is pushed one way only
        assert!(
            r.force.iter().all(|&f| f >= -1e-12 * p0),
            "Z {z} dt {dt}: the gas force changes sign {:?}",
            r.force
        );
        // the pair of the exchange cancels to the bit; the parcel account
        // closes to round-off
        assert_eq!(r.pair_residual, 0.0, "Z {z} dt {dt}");
        assert!(
            r.parcel_imbalance.abs() < 1e-12 * p0,
            "{}",
            r.parcel_imbalance
        );
    }
}

#[test]
fn the_coupling_converges_at_first_order_in_the_step() {
    // V28.4: Z = 1 at t = 10 ms, Δt = 2, 1, 0.5 ms
    let z = 1.0;
    let t = 0.01;
    let (u_ref, _) = reference(z, t);
    let err: Vec<f64> = [2e-3, 1e-3, 5e-4]
        .iter()
        .map(|&dt| {
            let r = go("order", &case(z, dt, t, "\"converged\""));
            (r.u - u_ref).abs()
        })
        .collect();
    let p1 = (err[0] / err[1]).log2();
    let p2 = (err[1] / err[2]).log2();
    eprintln!("errors {err:?}, observed orders {p1:.3} {p2:.3} (u_ref {u_ref})");
    assert!(p1 > 0.9 && p2 > 0.9, "orders {p1} {p2}, errors {err:?}");
    assert!(err[1] < 0.05 * u_ref, "error at 1 ms {}", err[1]);
}

#[test]
fn one_pass_is_bounded_and_conservative_too() {
    // D252: a fixed single pass is the linearised implicit coupling: still
    // conservative to round-off and without overshoot at Δt = 100 τ_p
    let z = 3.0;
    // (five steps: the cloud drifts with the gas at 0.75 m/s and stays
    // clear of the end wall)
    let r = go("one", &case(z, 0.77, 3.85, "1"));
    let p0 = z * RHO_G * V0;
    assert!(
        (r.gas_impulse + r.p_parcels - p0).abs() < 1e-12 * p0,
        "{} + {} against {p0}",
        r.gas_impulse,
        r.p_parcels
    );
    let total = RHO_G * r.u + r.p_parcels;
    assert!((total - p0).abs() < 1e-5 * p0, "{total} {p0}");
    assert!(r.force.iter().all(|&f| f >= -1e-12 * p0), "{:?}", r.force);
    assert_eq!(r.pair_residual, 0.0);
}

#[test]
fn passes_give_the_same_bits_on_any_number_of_threads() {
    // V28.1: a cone of parcels (keyed directions), 1 and 4 threads
    // (and with the heat exchange of SP2b)
    for text in [
        case(1.0, 2e-3, 0.02, "\"converged\""),
        heat_case(1.0, 2e-3, 0.02, "\"converged\""),
    ] {
        let text = text.replace("speed = 1", "half_angle = 40.0\nspeed = 1");
        let runs: Vec<String> = [1, 4]
            .iter()
            .map(|&n| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(n)
                    .build()
                    .unwrap()
                    .install(|| go(&format!("threads{n}"), &text).csv)
            })
            .collect();
        assert_eq!(runs[0], runs[1]);
    }
}

// Sensible heat (SP2b): the same cell of gas, now at 400 K, with the cloud
// injected at 300 K. Gas and parcels relax to one temperature with the
// sensible enthalpy of the mixture kept (the walls are adiabatic).

const T_GAS: f64 = 400.0;
const T_LIQ: f64 = 300.0;
const CP_G: f64 = 1006.0;
const CP_L: f64 = 4181.0;

fn heat_case(z: f64, dt: f64, end: f64, passes: &str) -> String {
    case(z, dt, end, passes)
        .replace(
            "[fluid]\n",
            &format!("[fluid]\nspecific_heat = {CP_G}\nconductivity = 0.026\n"),
        )
        .replace(
            "[boundary.xmin]\n",
            &format!("[energy]\ninitial = {T_GAS}\n[boundary.xmin]\n"),
        )
        .replace(
            "name = \"cloud\"\n",
            &format!("name = \"cloud\"\ntemperature = {T_LIQ}\n"),
        )
}

struct HeatRun {
    t_gas: f64,
    t_parcels: f64,
    gas_heat: f64,
    h_parcels: f64,
    heat: Vec<f64>,
    pair: f64,
    parcel_imbalance: f64,
    m_parcels: f64,
}

fn go_heat(name: &str, text: &str) -> HeatRun {
    let d = out(name);
    let case = FlowCase::from_toml(text, Path::new(".")).unwrap();
    let r = run(&case, name, &d, |_| {}).unwrap();
    let s = r.solves.last().unwrap();
    let bal = s.balance.as_ref().expect("balance");
    let run = bal.run.as_ref().expect("run ledger");
    let acc = run
        .ledger
        .iter()
        .find(|a| a.quantity == "parcel enthalpy")
        .expect("parcel enthalpy account");
    let t = s
        .scalars
        .iter()
        .find(|x| x.name == "temperature")
        .expect("temperature");
    let csv = std::fs::read_to_string(d.join("spray.csv")).unwrap();
    let mut lines = csv.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let (ch, ct) = (col("gas_heat_w"), col("t_mean_k"));
    let rows: Vec<Vec<f64>> = lines
        .map(|l| {
            l.split(',')
                .map(|x| x.parse::<f64>().unwrap_or(f64::NAN))
                .collect()
        })
        .collect();
    let sp = s.spray.as_ref().unwrap();
    let out = HeatRun {
        t_gas: 0.5 * (t.min + t.max),
        t_parcels: rows.last().unwrap()[ct],
        gas_heat: sp.coupling.as_ref().unwrap().gas_heat,
        h_parcels: acc.content.unwrap(),
        heat: rows.iter().map(|r| r[ch]).collect(),
        pair: run
            .pairs
            .iter()
            .find(|p| p.pair == 2804)
            .map_or(f64::NAN, |p| p.residual),
        parcel_imbalance: acc.imbalance,
        m_parcels: sp.in_flow_kg,
    };
    let _ = std::fs::remove_dir_all(&d);
    out
}

#[test]
fn a_cloud_and_its_gas_reach_one_temperature_at_any_step() {
    // τ_h ≈ 33 ms for 50 µm water at Nu ≈ 2; steps of 0.03 to 30 τ_h
    for (z, dt) in [(0.1f64, 1e-3f64), (1.0, 0.05), (1.0, 0.77)] {
        let end = 40.0 * 0.033 * (1.0 + CP_L * z / CP_G).max(1.0) + 3.0 * dt;
        let end = (end / dt).ceil() * dt;
        let r = go_heat("heat", &heat_case(z, dt, end, "\"converged\""));
        let (mg, ml) = (RHO_G * CP_G, z * RHO_G * CP_L);
        let h0 = mg * T_GAS + ml * T_LIQ;
        let te = h0 / (mg + ml);
        // the exchange books: the gas lost what the parcels gained
        assert!(
            (r.h_parcels - (ml * T_LIQ - r.gas_heat)).abs() < 1e-12 * h0,
            "Z {z} dt {dt}: parcels {} against {}",
            r.h_parcels,
            ml * T_LIQ - r.gas_heat
        );
        // the mixture keeps its enthalpy (the gas side through the energy
        // equation, to its solver tolerance)
        let total = mg * r.t_gas + CP_L * r.m_parcels * r.t_parcels;
        assert!(
            (total - h0).abs() < 1e-9 * h0,
            "Z {z} dt {dt}: {total} against {h0}"
        );
        // one common temperature
        assert!(
            (r.t_gas - te).abs() < 1e-9 * te,
            "Z {z} dt {dt}: gas {} against {te}",
            r.t_gas
        );
        assert!(
            (r.t_parcels - r.t_gas).abs() < 1e-9 * te,
            "Z {z} dt {dt}: parcels {} against gas {}",
            r.t_parcels,
            r.t_gas
        );
        // the gas only loses heat: no overshoot
        assert!(
            r.heat.iter().all(|&q| q <= 1e-12 * h0),
            "Z {z} dt {dt}: {:?}",
            r.heat
        );
        assert_eq!(r.pair, 0.0, "Z {z} dt {dt}");
        assert!(
            r.parcel_imbalance.abs() < 1e-12 * h0,
            "{}",
            r.parcel_imbalance
        );
    }
}

/// The lumped reference with heat: velocities and temperatures of gas
/// (mass M) and liquid (mass m), drag of Schiller and Naumann, Nusselt
/// number of Ranz and Marshall; RK4 to time t. Returns the gas temperature.
fn reference_heat(z: f64, t: f64) -> f64 {
    let (big, m) = (RHO_G, z * RHO_G);
    let pr = MU * CP_G / 0.026;
    let mut y = [0.0f64, V0, T_GAS, T_LIQ];
    let n = 200_000;
    let h = t / n as f64;
    let rhs = |y: [f64; 4]| {
        let re = RHO_G * D * (y[0] - y[1]).abs() / MU;
        let tau = RHO_L * D * D / (18.0 * MU * drag_correction(re));
        let f = m * (y[1] - y[0]) / tau;
        let tau_h = RHO_L * CP_L * D * D / (6.0 * nusselt(re, pr) * 0.026);
        let q = m * CP_L * (y[2] - y[3]) / tau_h;
        [f / big, -f / m, -q / (big * CP_G), q / (m * CP_L)]
    };
    let add = |a: [f64; 4], b: [f64; 4], s: f64| {
        [
            a[0] + s * b[0],
            a[1] + s * b[1],
            a[2] + s * b[2],
            a[3] + s * b[3],
        ]
    };
    for _ in 0..n {
        let k1 = rhs(y);
        let k2 = rhs(add(y, k1, 0.5 * h));
        let k3 = rhs(add(y, k2, 0.5 * h));
        let k4 = rhs(add(y, k3, h));
        for i in 0..4 {
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    y[2]
}

#[test]
fn the_heat_coupling_converges_at_first_order_in_the_step() {
    // Z = 1 at t = 20 ms (τ_h ≈ 33 ms), Δt = 4, 2, 1 ms
    let z = 1.0;
    let t = 0.02;
    let t_ref = reference_heat(z, t);
    let err: Vec<f64> = [4e-3, 2e-3, 1e-3]
        .iter()
        .map(|&dt| {
            let r = go_heat("heat_order", &heat_case(z, dt, t, "\"converged\""));
            (r.t_gas - t_ref).abs()
        })
        .collect();
    let p1 = (err[0] / err[1]).log2();
    let p2 = (err[1] / err[2]).log2();
    eprintln!("errors {err:?}, observed orders {p1:.3} {p2:.3} (T_ref {t_ref})");
    assert!(p1 > 0.9 && p2 > 0.9, "orders {p1} {p2}, errors {err:?}");
    assert!(err[2] < 0.05 * (T_GAS - t_ref), "error at 1 ms {}", err[2]);
}
