//! Verification of the dynamics (MS5 part 3): generalised-α against the
//! exact response of a damped oscillator, its order, its energy and its
//! high-frequency damping; mode superposition; Rayleigh damping from two
//! ratios; the Craig–Bampton reduction.
//!
//! References: Chung & Hulbert (1993), *J. Appl. Mech.* 60 (generalised-α:
//! second order, spectral radius ρ∞ at infinite frequency); Nigam &
//! Jennings (1969), *Bull. Seismol. Soc. Am.* 59 (exact integration of a
//! piecewise linear load); Craig & Bampton (1968), *AIAA J.* 6; the step
//! response of a damped oscillator, u = F/k [1 − e^(−ζωt)(cos ω_d t +
//! ζω/ω_d sin ω_d t)].

#![allow(clippy::needless_range_loop)]

use foehncraft_structure::case::{self, Case};
use foehncraft_structure::dynamics::{
    self, GenAlpha, ModalBasis, Rayleigh, craig_bampton, dense_eigenvalues,
};
use foehncraft_structure::model::{DOFS, Model};
use foehncraft_structure::{analysis, codes};
use std::f64::consts::PI;

fn rel(a: f64, b: f64) -> f64 {
    ((a - b) / b).abs()
}

/// One mass on one spring along x: m = 1 kg, f = 10 Hz.
fn oscillator(force: f64) -> (Model, f64) {
    let mut m = Model::default();
    m.add_node([0.0; 3]);
    for d in 1..DOFS {
        m.fixed[d] = true;
    }
    let w = 2.0 * PI * 10.0;
    m.springs.push((0, w * w));
    m.masses.push((0, 1.0));
    m.forces[0] = force;
    (m, w)
}

fn step_exact(t: f64, f: f64, k: f64, w: f64, z: f64) -> f64 {
    let wd = w * (1.0 - z * z).sqrt();
    f / k * (1.0 - (-z * w * t).exp() * ((wd * t).cos() + z * w / wd * (wd * t).sin()))
}

fn gen_alpha_end(rho: f64, z: f64, steps: usize, t_end: f64) -> f64 {
    let (model, w) = oscillator(3.0);
    let sys = model.assemble();
    let damp = Rayleigh {
        alpha: 2.0 * z * w,
        beta: 0.0,
    };
    let ga = GenAlpha::new(&sys, damp, t_end / steps as f64, rho).unwrap();
    let f = sys.f.clone();
    let mut s = ga.initial(0.0, vec![0.0], vec![0.0], &f).unwrap();
    for _ in 0..steps {
        s = ga.step(&s, &f, &f);
    }
    s.u[0]
}

#[test]
fn generalised_alpha_converges_at_second_order_to_the_damped_step_response() {
    let (z, t_end) = (0.05, 0.23);
    let w = 2.0 * PI * 10.0;
    let exact = step_exact(t_end, 3.0, w * w, w, z);
    for rho in [1.0, 0.8, 0.5, 0.0] {
        let e: Vec<f64> = [40, 80, 160, 320]
            .iter()
            .map(|&n| (gen_alpha_end(rho, z, n, t_end) - exact).abs())
            .collect();
        let order = (e[2] / e[3]).log2();
        eprintln!("ρ∞ = {rho}: errors {e:?}, order {order:.3}");
        assert!((1.9..=2.1).contains(&order), "ρ∞ = {rho}: {e:?}");
        assert!(e[3] < 1e-2 * exact.abs(), "ρ∞ = {rho}: {e:?}");
    }
}

#[test]
fn generalised_alpha_damps_high_frequencies_to_rho_infinity() {
    // a step 10⁴ times the period: the amplification matrix of (u, v, a)
    // has the triple eigenvalue −ρ∞ (Chung & Hulbert 1993), so its
    // characteristic polynomial is (λ + ρ∞)³
    for rho in [0.0, 0.5, 0.8, 1.0] {
        let (model, _) = oscillator(0.0);
        let sys = model.assemble();
        let dt = 1000.0;
        let ga = GenAlpha::new(&sys, Rayleigh::default(), dt, rho).unwrap();
        let mut amp = [[0.0; 3]; 3];
        for c in 0..3 {
            let mut x = [0.0; 3];
            x[c] = 1.0;
            let s0 = dynamics::State {
                t: 0.0,
                u: vec![x[0]],
                v: vec![x[1] / dt],
                a: vec![x[2] / (dt * dt)],
            };
            let s1 = ga.step(&s0, &[0.0], &[0.0]);
            amp[0][c] = s1.u[0];
            amp[1][c] = s1.v[0] * dt;
            amp[2][c] = s1.a[0] * dt * dt;
        }
        let tr = amp[0][0] + amp[1][1] + amp[2][2];
        let minors = amp[0][0] * amp[1][1] - amp[0][1] * amp[1][0] + amp[0][0] * amp[2][2]
            - amp[0][2] * amp[2][0]
            + amp[1][1] * amp[2][2]
            - amp[1][2] * amp[2][1];
        let det = amp[0][0] * (amp[1][1] * amp[2][2] - amp[1][2] * amp[2][1])
            - amp[0][1] * (amp[1][0] * amp[2][2] - amp[1][2] * amp[2][0])
            + amp[0][2] * (amp[1][0] * amp[2][1] - amp[1][1] * amp[2][0]);
        eprintln!("ρ∞ = {rho}: trace {tr:.6}, minors {minors:.6}, det {det:.6}");
        assert!((tr + 3.0 * rho).abs() < 1e-3, "{rho}: {tr}");
        assert!((minors - 3.0 * rho * rho).abs() < 1e-3, "{rho}: {minors}");
        assert!((det + rho.powi(3)).abs() < 1e-3, "{rho}: {det}");
    }
}

fn cantilever(elements: usize, planar: bool) -> (Model, usize) {
    let text = format!(
        r#"
[analysis]
static = true
[materials.steel]
young = 206e9
poisson = 0.3
density = 7850.0
[sections.reed]
shape = "rectangle"
width = 6e-3
thickness = 0.2e-3
[[beams]]
name = "reed"
from = [0.0, 0.0, 0.0]
to = [0.03, 0.0, 0.0]
elements = {elements}
section = "reed"
material = "steel"
theory = "euler-bernoulli"
[[supports]]
at = [0.0, 0.0, 0.0]
fix = ["all"]
[points]
tip = [0.03, 0.0, 0.0]
"#
    );
    let b = case::build(&Case::from_toml(&text).unwrap()).unwrap();
    let mut m = b.model;
    if planar {
        // the flap alone: uz and ry
        for n in 0..m.nodes.len() {
            for d in [0, 1, 3, 5] {
                m.fixed[DOFS * n + d] = true;
            }
        }
    }
    let tip = b.points["tip"];
    (m, tip)
}

#[test]
fn trapezoidal_rule_conserves_the_energy() {
    // ρ∞ = 1, no damping, a constant tip force from rest: ½vᵀMv + ½uᵀKu − fᵀu
    // stays 0 to rounding
    let (mut m, tip) = cantilever(20, false);
    m.forces[DOFS * tip + 2] = 0.05;
    let sys = m.assemble();
    let ga = GenAlpha::new(&sys, Rayleigh::default(), 2e-5, 1.0).unwrap();
    let f = sys.f.clone();
    let n = f.len();
    let mut s = ga.initial(0.0, vec![0.0; n], vec![0.0; n], &f).unwrap();
    let mut worst: f64 = 0.0;
    let mut largest: f64 = 0.0;
    for _ in 0..2000 {
        s = ga.step(&s, &f, &f);
        let (kin, strain) = s.energies(&sys);
        let work: f64 = f.iter().zip(&s.u).map(|(a, b)| a * b).sum();
        worst = worst.max((kin + strain - work).abs());
        largest = largest.max(strain);
    }
    eprintln!("energy error {worst:.3e} of the largest strain energy {largest:.3e}");
    assert!(worst < 1e-10 * largest, "{worst} {largest}");
}

#[test]
fn mode_superposition_is_exact_per_mode_and_meets_the_direct_method() {
    // the oscillator: the modal step is exact for any step length
    let (model, w) = oscillator(3.0);
    let z = 0.05;
    let basis = ModalBasis::lowest(&model, 1, z, Rayleigh::default()).unwrap();
    assert!(rel(basis.omega[0], w) < 1e-12);
    let sys = model.assemble();
    let p = basis.project(&sys.f);
    let (mut q, mut qd) = (vec![0.0], vec![0.0]);
    let dt = 0.037;
    for k in 1..=10 {
        basis.step(&mut q, &mut qd, &p, &p, dt);
        let u = basis.expand(&q)[0];
        let exact = step_exact(k as f64 * dt, 3.0, w * w, w, z);
        assert!(
            (u - exact).abs() < 1e-12 * exact.abs().max(3.0 / (w * w)),
            "{u} {exact}"
        );
    }
    // a ramp: the same end state with 1 and with 1000 steps
    let ramp = |steps: usize| {
        let (mut q, mut qd) = (vec![0.0], vec![0.0]);
        let h = 0.05 / steps as f64;
        for k in 0..steps {
            let p0 = vec![p[0] * (k as f64 * h) / 0.05];
            let p1 = vec![p[0] * ((k + 1) as f64 * h) / 0.05];
            basis.step(&mut q, &mut qd, &p0, &p1, h);
        }
        (q[0], qd[0])
    };
    let (a, b) = (ramp(1), ramp(1000));
    assert!(
        rel(a.0, b.0) < 1e-10 && rel(a.1, b.1) < 1e-10,
        "{a:?} {b:?}"
    );

    // the planar cantilever with every mode: the modal and the direct
    // responses to a step tip force meet as Δt → 0
    let (mut m, tip) = cantilever(8, true);
    m.forces[DOFS * tip + 2] = 0.05;
    let sys = m.assemble();
    let nf = sys.free.len();
    let basis = ModalBasis::lowest(
        &m,
        nf,
        0.0,
        Rayleigh {
            alpha: 20.0,
            beta: 1e-6,
        },
    )
    .unwrap();
    let ray = Rayleigh {
        alpha: 20.0,
        beta: 1e-6,
    };
    let t_end = 0.01;
    let pm = basis.project(&sys.f);
    let (mut q, mut qd) = (vec![0.0; nf], vec![0.0; nf]);
    basis.step(&mut q, &mut qd, &pm, &pm, t_end);
    let iz = sys.map[DOFS * tip + 2].unwrap();
    let modal = basis.expand(&q)[iz];
    let direct = |steps: usize| {
        let ga = GenAlpha::new(&sys, ray, t_end / steps as f64, 1.0).unwrap();
        let mut s = ga
            .initial(0.0, vec![0.0; nf], vec![0.0; nf], &sys.f)
            .unwrap();
        for _ in 0..steps {
            s = ga.step(&s, &sys.f, &sys.f);
        }
        s.u[iz]
    };
    let (d1, d2) = (direct(4000), direct(8000));
    let order = ((d1 - modal) / (d2 - modal)).abs().log2();
    eprintln!("tip: modal {modal:.9e}, direct {d1:.9e}, {d2:.9e}, order {order:.2}");
    assert!(rel(d2, modal) < 1e-5, "{d2} {modal}");
    assert!((1.8..=2.2).contains(&order), "{order}");
}

#[test]
fn rayleigh_factors_from_two_ratios() {
    let r = Rayleigh::from_ratios(100.0, 0.02, 2000.0, 0.03).unwrap();
    assert!(rel(r.ratio(2.0 * PI * 100.0), 0.02) < 1e-12);
    assert!(rel(r.ratio(2.0 * PI * 2000.0), 0.03) < 1e-12);
    // a ratio that falls too fast for αM + βK
    let e = Rayleigh::from_ratios(100.0, 0.05, 200.0, 0.001).unwrap_err();
    assert_eq!(e.code, codes::CASE);
    assert_eq!(
        Rayleigh::from_ratios(100.0, 0.02, 100.0, 0.03)
            .unwrap_err()
            .code,
        codes::CASE
    );
    // the oscillator damped by the fitted factors decays at its ratio
    let (model, w) = oscillator(0.0);
    let sys = model.assemble();
    let f0 = w / (2.0 * PI);
    let r = Rayleigh::from_ratios(f0, 0.02, 4.0 * f0, 0.05).unwrap();
    let ga = GenAlpha::new(&sys, r, 1.0 / (f0 * 400.0), 1.0).unwrap();
    let mut s = ga.initial(0.0, vec![1.0], vec![0.0], &[0.0]).unwrap();
    // peaks one period apart: ln(u1/u2) = 2πζ/√(1−ζ²)
    let mut peaks = Vec::new();
    let mut prev = (s.u[0], s.u[0]);
    for _ in 0..4000 {
        s = ga.step(&s, &[0.0], &[0.0]);
        if prev.1 > prev.0 && prev.1 > s.u[0] {
            peaks.push(prev.1);
        }
        prev = (prev.1, s.u[0]);
    }
    let delta = (peaks[1] / peaks[5]).ln() / 4.0;
    let zeta = delta / (4.0 * PI * PI + delta * delta).sqrt();
    eprintln!("decay ratio {zeta:.6}");
    assert!(rel(zeta, 0.02) < 2e-3, "{zeta}");
}

#[test]
fn craig_bampton_is_exact_with_every_mode_and_a_ritz_bound_with_few() {
    let (m, tip) = cantilever(20, true);
    let boundary = [DOFS * tip + 2, DOFS * tip + 4];
    let sys = m.assemble();
    let nf = sys.free.len();
    let full = analysis::solve_modal(&m, 8, 0.0).unwrap();
    let full: Vec<f64> = full.modes.iter().map(|x| x.lambda).collect();
    // every interior mode: the eigenvalues of the full model
    let all = craig_bampton(&m, &boundary, nf - 2).unwrap();
    let ev = dense_eigenvalues(&all.k, &all.m, all.n()).unwrap();
    for i in 0..8 {
        assert!(
            rel(ev[i], full[i]) < 1e-8,
            "mode {i}: {} {}",
            ev[i],
            full[i]
        );
    }
    // a few modes: above the full ones (Rayleigh–Ritz), converging
    let mut prev_err = f64::INFINITY;
    for keep in [2, 4, 8] {
        let r = craig_bampton(&m, &boundary, keep).unwrap();
        let ev = dense_eigenvalues(&r.k, &r.m, r.n()).unwrap();
        let err = rel(ev[2], full[2]);
        eprintln!("{keep} interior modes: third eigenvalue off by {err:.3e}");
        for i in 0..3 {
            assert!(
                ev[i] >= full[i] * (1.0 - 1e-10),
                "mode {i}: {} {}",
                ev[i],
                full[i]
            );
        }
        assert!(err < prev_err);
        prev_err = err;
    }
    // the constraint modes make the static response at the boundary exact
    let r = craig_bampton(&m, &boundary, 3).unwrap();
    let mut mt = m.clone();
    mt.forces[DOFS * tip + 2] = 0.05;
    let st = analysis::solve_static(&mt).unwrap();
    // boundary block of K: [[k_uu, k_ur], [k_ru, k_rr]] · x = [0.05, 0]
    let n = r.n();
    let (a, b, c) = (r.k[0], r.k[1], r.k[n + 1]);
    let w = 0.05 * c / (a * c - b * b);
    assert!(
        rel(w, st.u[DOFS * tip + 2]) < 1e-10,
        "{w} {}",
        st.u[DOFS * tip + 2]
    );
    // the constraint modes are K-orthogonal to the fixed-interface modes
    for i in 0..2 {
        for j in 2..n {
            assert!(r.k[i * n + j].abs() < 1e-8 * r.k[0].abs(), "{i} {j}");
        }
    }
    // errors
    assert_eq!(craig_bampton(&m, &[0], 2).unwrap_err().code, codes::CASE);
    assert_eq!(
        dynamics::craig_bampton(&m, &[boundary[0], boundary[0]], 2)
            .unwrap_err()
            .code,
        codes::CASE
    );
}
