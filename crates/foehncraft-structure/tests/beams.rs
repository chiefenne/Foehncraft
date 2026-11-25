//! Verification of the beams, the static and modal analyses and the clamp
//! springs against analytic beam theory (MS5; plan section 6.4, rows
//! "cantilever beam" and "clamp spring calibration").
//!
//! References: Timoshenko & Gere, *Mechanics of Materials*; Blevins,
//! *Formulas for Natural Frequency and Mode Shape* (1979), table 8-1
//! (β_n L of the cantilever and of the free–free beam); Han, Benaroya &
//! Wei (1999), *J. Sound Vib.* 225, the frequency equation of the pinned
//! Timoshenko beam; the beam with a rotational root spring from its four
//! boundary conditions (solved here by bisection).

#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use foehncraft_structure::case::{self, Case};
use foehncraft_structure::clamp::{self, Which};
use foehncraft_structure::codes;
use std::f64::consts::PI;

const E: f64 = 206e9;
const NU: f64 = 0.3;
const RHO: f64 = 7850.0;
/// A reed of spring steel: 30 mm free length, 6 mm wide, 0.2 mm thick.
const L: f64 = 0.03;
const B: f64 = 6e-3;
const H: f64 = 0.2e-3;

fn reed(theory: &str, elements: usize, extra: &str) -> Case {
    let text = format!(
        r#"
[analysis]
static = true
modes = 5
[materials.steel]
young = {E}
poisson = {NU}
density = {RHO}
[sections.reed]
shape = "rectangle"
width = {B}
thickness = {H}
[[beams]]
name = "reed"
from = [0.0, 0.0, 0.0]
to = [{L}, 0.0, 0.0]
elements = {elements}
section = "reed"
material = "steel"
theory = "{theory}"
[points]
tip = [{L}, 0.0, 0.0]
{extra}
"#
    );
    Case::from_toml(&text).unwrap()
}

const RIGID: &str = r#"
[[supports]]
at = [0.0, 0.0, 0.0]
fix = ["all"]
"#;

fn ei() -> f64 {
    E * B * H.powi(3) / 12.0
}
fn kga() -> f64 {
    let g = E / (2.0 * (1.0 + NU));
    let kappa = 10.0 * (1.0 + NU) / (12.0 + 11.0 * NU);
    kappa * g * B * H
}
fn tip(out: &case::Outcome) -> (f64, f64) {
    let s = out.statics.as_ref().unwrap();
    let n = out.built.points["tip"];
    (s.u[6 * n + 2], s.u[6 * n + 4])
}
fn rel(a: f64, b: f64) -> f64 {
    ((a - b) / b).abs()
}

#[test]
fn cantilever_tip_load_is_exact() {
    let p = -0.05; // N, downwards
    for (theory, shear) in [("euler-bernoulli", false), ("timoshenko", true)] {
        for n in [1, 7, 40] {
            let c = reed(
                theory,
                n,
                &format!("{RIGID}\n[[loads]]\nat = [{L}, 0.0, 0.0]\nforce = [0.0, 0.0, {p}]\n"),
            );
            let out = case::run(&c).unwrap();
            let (w, ry) = tip(&out);
            let mut exact = p * L.powi(3) / (3.0 * ei());
            if shear {
                exact += p * L / kga();
            }
            // rounding: the reed's axial stiffness is 10⁸ times its bending one
            assert!(rel(w, exact) < 1e-8, "{theory} {n}: {w} {exact}");
            // θy = −dw/dx; dw/dx at the tip = P L²/(2EI)
            assert!(
                rel(-ry, p * L * L / (2.0 * ei())) < 1e-8,
                "{theory} {n}: {ry}"
            );
            // the reaction balances the load
            let rf = out.to_json("t")["static"]["reaction_force"][2]
                .as_f64()
                .unwrap();
            assert!((rf + p).abs() < 1e-8 * p.abs(), "{rf}");
        }
    }
}

#[test]
fn cantilever_distributed_load_is_exact_at_the_nodes() {
    let q = -20.0; // N/m
    for (theory, shear) in [("euler-bernoulli", false), ("timoshenko", true)] {
        let c = reed(
            theory,
            12,
            &format!(
                "{RIGID}\n[[line_loads]]\nbeam = \"reed\"\nforce_per_length = [0.0, 0.0, {q}]\n"
            ),
        );
        let out = case::run(&c).unwrap();
        let s = out.statics.as_ref().unwrap();
        let m = &out.built.model;
        for (nd, x) in m.nodes.iter().enumerate() {
            let x = x[0];
            // w(x) = q x²(6L² − 4Lx + x²)/(24EI) (+ q(Lx − x²/2)/(κGA))
            let mut exact = q * x * x * (6.0 * L * L - 4.0 * L * x + x * x) / (24.0 * ei());
            if shear {
                exact += q * (L * x - 0.5 * x * x) / kga();
            }
            let w = s.u[6 * nd + 2];
            assert!(
                (w - exact).abs() <= 1e-10 * (q * L.powi(4) / (8.0 * ei())).abs(),
                "{theory} x={x}: {w} {exact}"
            );
        }
    }
}

/// β_n L of the clamped–free beam (Blevins 1979, table 8-1).
const CANTILEVER: [f64; 5] = [
    1.875_104_068_711_961,
    4.694_091_132_974_175,
    7.854_757_438_237_613,
    10.995_540_734_875_47,
    14.137_168_391_046_47,
];

fn eb_frequency(beta_l: f64) -> f64 {
    beta_l * beta_l / (L * L) * (ei() / (RHO * B * H)).sqrt() / (2.0 * PI)
}

/// The flapping (x–z) modes of the reed, the 0.2 mm direction: modes
/// whose shape is mostly w. The torsion modes (from 1.8 kHz, no
/// translation of the axis) and the in-plane modes (30 times the flapping
/// ones) lie between and above them.
fn flap_modes(out: &case::Outcome) -> Vec<f64> {
    let md = out.modal.as_ref().unwrap();
    md.modes
        .iter()
        .filter(|m| {
            let s = &m.shape;
            let n = s.len() / 6;
            let w: f64 = (0..n).map(|i| s[6 * i + 2].powi(2)).sum();
            let v: f64 = (0..n).map(|i| s[6 * i + 1].powi(2)).sum();
            let t: f64 = (0..n).map(|i| (s[6 * i + 3] * B).powi(2)).sum();
            w > v && w > t
        })
        .map(|m| m.frequency())
        .collect()
}

#[test]
fn cantilever_frequencies_converge_at_fourth_order() {
    let mut errs = Vec::new();
    for n in [10, 20, 40] {
        let mut c = reed("euler-bernoulli", n, RIGID);
        c.analysis.modes = 10;
        let out = case::run(&c).unwrap();
        let md = out.modal.as_ref().unwrap();
        assert!(md.sturm.ok(), "{:?}", md.sturm);
        assert!(md.converged);
        let f = flap_modes(&out);
        assert_eq!(f.len(), 5, "{f:?}");
        let e: Vec<f64> = (0..5)
            .map(|i| rel(f[i], eb_frequency(CANTILEVER[i])))
            .collect();
        errs.push(e);
    }
    // 40 elements: every one of the first five within 2e-5 (mode 5: 1.1e-5)
    for (i, e) in errs[2].iter().enumerate() {
        assert!(*e < 2e-5, "mode {}: {e:e}", i + 1);
    }
    // observed order from 20 to 40 elements, modes 2 to 5 (mode 1 is at
    // rounding level already)
    for i in 1..5 {
        let order = (errs[1][i] / errs[2][i]).log2();
        assert!((3.8..4.2).contains(&order), "mode {}: order {order}", i + 1);
    }
    // the modes are pure flapping: all effective mass in z
    let out = case::run(&reed("euler-bernoulli", 40, RIGID)).unwrap();
    let j = out.to_json("t");
    let em = j["modal"]["modes"][0]["effective_mass_fraction"]
        .as_array()
        .unwrap()
        .clone();
    // the classical 61.3 % of the cantilever's first mode
    assert!((em[2].as_f64().unwrap() - 0.6131).abs() < 1e-3, "{em:?}");
    assert!(em[0].as_f64().unwrap() < 1e-12 && em[1].as_f64().unwrap() < 1e-12);
}

/// The frequency equation of the pinned Timoshenko beam in mode n: the
/// lower root ω² of
/// (ρAω² − κGAk²)(ρIω² − EIk² − κGA) − (κGAk)² = 0, k = nπ/L.
fn pinned_timoshenko(
    n: usize,
    len: f64,
    e: f64,
    g: f64,
    rho: f64,
    a: f64,
    i: f64,
    kappa: f64,
) -> f64 {
    let k = n as f64 * PI / len;
    let s = kappa * g * a;
    // a2 x² + a1 x + a0 = 0 with x = ω²
    let a2 = rho * a * rho * i;
    let a1 = -(rho * a * (e * i * k * k + s) + rho * i * s * k * k);
    let a0 = s * k * k * (e * i * k * k + s) - s * s * k * k;
    let disc = (a1 * a1 - 4.0 * a2 * a0).sqrt();
    let x = (-a1 - disc) / (2.0 * a2);
    x.sqrt() / (2.0 * PI)
}

#[test]
fn thick_pinned_timoshenko_beam() {
    // the two-node element converges at second order in the frequency here
    // (the Euler–Bernoulli limit at fourth): 120 elements give 7e-6
    let mut errs: Vec<Vec<f64>> = Vec::new();
    for ne in [60usize, 120] {
        // a round bar, L/d = 10: shear and rotary inertia lower the second
        // mode by about 7 % against Euler–Bernoulli
        let (len, d) = (1.0, 0.1);
        let text = format!(
            r#"
[analysis]
modes = 4
[materials.steel]
young = {E}
poisson = {NU}
density = {RHO}
[sections.bar]
shape = "circle"
diameter = {d}
[[beams]]
name = "bar"
from = [0.0, 0.0, 0.0]
to = [{len}, 0.0, 0.0]
elements = {ne}
section = "bar"
material = "steel"
[[supports]]
at = [0.0, 0.0, 0.0]
fix = ["ux", "uy", "uz", "rx"]
[[supports]]
at = [{len}, 0.0, 0.0]
fix = ["ux", "uy", "uz", "rx"]
"#
        );
        let out = case::run(&Case::from_toml(&text).unwrap()).unwrap();
        let md = out.modal.as_ref().unwrap();
        assert!(md.sturm.ok(), "{:?}", md.sturm);
        let g = E / (2.0 * (1.0 + NU));
        let a = PI * d * d / 4.0;
        let i = PI * d.powi(4) / 64.0;
        let kappa = 6.0 * (1.0 + NU) / (7.0 + 6.0 * NU);
        // two planes: each frequency twice
        for (m, mode) in md.modes.iter().enumerate() {
            let n = m / 2 + 1;
            let exact = pinned_timoshenko(n, len, E, g, RHO, a, i, kappa);
            let f = mode.frequency();
            if m == 0 {
                errs.push(Vec::new());
            }
            errs.last_mut().unwrap().push(rel(f, exact));
            let eb = (n as f64 * PI / len).powi(2) * (E * i / (RHO * a)).sqrt() / (2.0 * PI);
            assert!(f < eb, "Timoshenko lies below Euler–Bernoulli");
        }
    }
    for m in 0..4 {
        assert!(errs[1][m] < 1e-5, "mode {}: {:e}", m + 1, errs[1][m]);
        let order = (errs[0][m] / errs[1][m]).log2();
        assert!((1.9..2.1).contains(&order), "mode {}: order {order}", m + 1);
    }
}

/// The first root β L of the cantilever with a rotational root spring k
/// (N·m/rad): w(0) = 0, EI w''(0) = k w'(0), w''(L) = w'''(L) = 0, with
/// w = A(cosh − cos) + B sinh + D sin; the 3 × 3 determinant in (A, B, D).
fn spring_root_beta(k: f64) -> f64 {
    let det = |bl: f64| {
        let beta = bl / L;
        let (ch, sh, c, s) = (bl.cosh(), bl.sinh(), bl.cos(), bl.sin());
        let m = [
            [2.0 * ei() * beta, -k, -k],
            [ch + c, sh, -s],
            [sh - s, ch, -c],
        ];
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    // scan (0, 1.875] for the sign change, then bisect
    let mut a = 1e-3;
    let mut fa = det(a);
    let mut b = a;
    while b < CANTILEVER[0] + 1e-3 {
        b = a + 1e-3;
        let fb = det(b);
        if (fb < 0.0) != (fa < 0.0) {
            break;
        }
        a = b;
        fa = fb;
    }
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        let fm = det(m);
        if (fm < 0.0) == (fa < 0.0) {
            a = m;
            fa = fm;
        } else {
            b = m;
        }
    }
    0.5 * (a + b)
}

fn clamp_keys(kr: Option<f64>, kt: Option<f64>) -> String {
    let mut s = "[[clamps]]\nat = [0.0, 0.0, 0.0]\nnormal = \"z\"\nline = \"y\"\n".to_string();
    if let Some(k) = kr {
        s += &format!("k_rotation = {k}\n");
    }
    if let Some(k) = kt {
        s += &format!("k_translation = {k}\n");
    }
    s
}

#[test]
fn clamp_springs_static_is_exact() {
    // per metre of clamp line; the clamp line is the 6 mm width
    let (kr, kt) = (40.0, 2.0e6);
    let p = -0.05;
    let c = reed(
        "euler-bernoulli",
        10,
        &format!(
            "{}\n[[loads]]\nat = [{L}, 0.0, 0.0]\nforce = [0.0, 0.0, {p}]\n",
            clamp_keys(Some(kr), Some(kt))
        ),
    );
    let out = case::run(&c).unwrap();
    let (w, _) = tip(&out);
    let exact = p * L.powi(3) / (3.0 * ei()) + p * L * L / (kr * B) + p / (kt * B);
    assert!(rel(w, exact) < 1e-10, "{w} {exact}");
    let j = out.to_json("t");
    assert_eq!(j["clamps"][0]["rotation"]["rigid"], false);
    assert!(
        rel(
            j["clamps"][0]["rotation"]["stiffness"].as_f64().unwrap(),
            kr * B
        ) < 1e-15
    );
    // a clamp without springs is rigid, and the result says so
    let out = case::run(&reed("euler-bernoulli", 4, &clamp_keys(None, None))).unwrap();
    let j = out.to_json("t");
    assert_eq!(j["clamps"][0]["rotation"]["rigid"], true);
    assert_eq!(j["clamps"][0]["translation"]["rigid"], true);
}

#[test]
fn clamp_spring_frequency_against_the_frequency_equation() {
    // from a soft to a nearly rigid mounting: k L/EI from 0.3 to 300
    for klei in [0.3, 3.0, 30.0, 300.0] {
        let k = klei * ei() / L; // total N·m/rad
        let c = reed("euler-bernoulli", 40, &clamp_keys(Some(k / B), None));
        let out = case::run(&c).unwrap();
        let f = flap_modes(&out)[0];
        let exact = eb_frequency(spring_root_beta(k));
        assert!(rel(f, exact) < 1e-8, "kL/EI = {klei}: {f} {exact}");
        assert!(f < eb_frequency(CANTILEVER[0]));
    }
}

#[test]
fn clamp_fit_recovers_the_spring() {
    let k_true = 3.0 * ei() / L / B; // per metre of clamp line
    let measured = eb_frequency(spring_root_beta(k_true * B));
    let c = reed("euler-bernoulli", 40, &clamp_keys(Some(1.0), None));
    let fit = clamp::fit(&c, 0, Which::Rotation, 1, measured).unwrap();
    assert!(rel(fit.frequency, measured) < 1e-9, "{fit:?}");
    assert!(
        rel(fit.k_per_length, k_true) < 1e-5,
        "{} {k_true}",
        fit.k_per_length
    );
    assert!(fit.frequency_rigid > measured);
    assert!(fit.evaluations < 40, "{}", fit.evaluations);
    // above the rigid clamp: no spring can do it
    let e = clamp::fit(&c, 0, Which::Rotation, 1, 1.01 * fit.frequency_rigid).unwrap_err();
    assert_eq!(e.code, codes::FIT);
}

#[test]
fn free_free_beam_with_a_negative_shift() {
    // six rigid-body modes at zero, then the free–free bending modes
    // β₁L = 4.730 (Blevins table 8-1)
    let mut c = reed("euler-bernoulli", 40, "");
    c.analysis.statics = false;
    c.analysis.modes = 8;
    // well below the first elastic λ = (2π·1.2 kHz)², far from 0 in the
    // scale of K
    c.analysis.shift = -1e4;
    let out = case::run(&c).unwrap();
    let md = out.modal.as_ref().unwrap();
    assert!(md.sturm.ok(), "{:?}", md.sturm);
    let f: Vec<f64> = md.modes.iter().map(|m| m.frequency()).collect();
    let first = eb_frequency(CANTILEVER[0]);
    for (i, fi) in f.iter().take(6).enumerate() {
        assert!(*fi < 1e-3 * first, "rigid mode {}: {fi}", i + 1);
    }
    let ff = eb_frequency(4.730_040_744_862_704);
    assert!(rel(f[6], ff) < 1e-6, "{} {ff}", f[6]);
    // the rigid modes carry the whole mass in each direction
    for d in 0..3 {
        let s: f64 = md.modes.iter().take(6).map(|m| m.effective_mass[d]).sum();
        assert!(rel(s, md.total_mass) < 1e-9, "direction {d}: {s}");
    }
}

#[test]
fn unsupported_static_case_is_refused_with_its_code() {
    let c = reed("timoshenko", 4, "");
    let e = case::run(&c).unwrap_err();
    assert_eq!(e.code, codes::SINGULAR);
    assert!(e.message.contains("free to move"), "{}", e.message);
}

#[test]
fn case_errors_carry_the_case_code() {
    let e = Case::from_toml("[analysis]\nstatics = true\n").unwrap_err();
    assert_eq!(e.code, codes::CASE);
    let c = reed(
        "timoshenko",
        4,
        "[[supports]]\nat = [0.001, 0.0, 0.0]\nfix = [\"all\"]\n",
    );
    let e = case::run(&c).unwrap_err();
    assert_eq!(e.code, codes::CASE);
    assert!(e.message.contains("not a node"), "{}", e.message);
}
