//! Transposed matrices of the flow solver (requirement D59, item 4): Aᵀ in
//! the storage of A and the product Aᵀx without storing Aᵀ, both gathers
//! without atomics, on the momentum and pressure-correction matrices of a
//! running SIMPLEC solve; the linear solvers and the multigrid on Aᵀ.

use crate::case::{FlowCase, mesh_for};
use crate::checkpoint::{FlowRun, reverse};
use crate::device::Buf;
use crate::simplec::{Flow, Physics, Target};
use crate::transient::Transient;
use crate::{bc, simplec};
use foehncraft_fv::solver::Preconditioner;
use std::path::Path;

/// A lid-driven cavity at Re 100 on perturbed polygons after a few outer
/// iterations, so the momentum matrix carries convection (non-symmetric).
fn cavity(threads: Option<usize>) -> Flow {
    let (case, mesh, geo) = cavity_case("");
    let phys = physics(&case);
    let arrays = bc::arrays(&mesh, &geo, &case.boundary, phys.rho, bc::At::numbers(None)).unwrap();
    let disc = simplec::discretisation(&mesh, &geo, &arrays).unwrap();
    let mut fl = Flow::new(
        disc,
        arrays,
        phys,
        case.solver.clone(),
        [0.0; 3],
        0.0,
        Some(0),
        &Target {
            device: "cpu".into(),
            threads,
            ptx: Default::default(),
            ranks: None,
        },
    )
    .unwrap();
    for _ in 0..5 {
        fl.iterate().unwrap();
    }
    fl
}

fn physics(case: &FlowCase) -> Physics {
    Physics::constant(case.rho(), case.mu(), [0.0; 3], true)
}

fn cavity_case(extra: &str) -> (FlowCase, foehncraft_mesh::Mesh, foehncraft_mesh::Geometry) {
    let text = format!(
        r#"
[mesh.generate]
type = "poly"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [14, 14, 1]
perturb = 0.15
seed = 59
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
{extra}
"#
    );
    let case = FlowCase::from_toml(&text, Path::new(".")).unwrap();
    let (mesh, geo, _) = mesh_for(&case, None).unwrap();
    (case, mesh, geo)
}

/// Deterministic values in [−1, 1).
fn values(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..n)
        .map(|_| {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (s >> 11) as f64 / (1u64 << 52) as f64 - 1.0
        })
        .collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Aᵀ of `val` as a new buffer, with the identities of the transpose.
fn transposed(fl: &mut Flow, val: Buf) -> Buf {
    let nc = fl.disc.n_cells;
    let stored = fl.disc.sell.stored();
    let dev = &mut fl.dev;
    let t = dev.alloc(stored).unwrap();
    dev.transpose(val, t).unwrap();
    // transposing twice gives A back, bit for bit
    let tt = dev.alloc(stored).unwrap();
    dev.transpose(t, tt).unwrap();
    assert_eq!(dev.download(tt).unwrap(), dev.download(val).unwrap());
    // ⟨A v, w⟩ = ⟨v, Aᵀ w⟩ to round-off
    let (v, w) = (dev.alloc(nc).unwrap(), dev.alloc(nc).unwrap());
    let (hv, hw) = (values(nc, 1), values(nc, 2));
    dev.upload(v, &hv).unwrap();
    dev.upload(w, &hw).unwrap();
    let (av, atw, atw2) = (
        dev.alloc(nc).unwrap(),
        dev.alloc(nc).unwrap(),
        dev.alloc(nc).unwrap(),
    );
    dev.matrix_product(val, v, av).unwrap();
    dev.transposed_product(val, w, atw).unwrap();
    dev.matrix_product(t, w, atw2).unwrap();
    let (av, atw, atw2) = (
        dev.download(av).unwrap(),
        dev.download(atw).unwrap(),
        dev.download(atw2).unwrap(),
    );
    // the product without storing Aᵀ is the product with the stored Aᵀ, bit
    // for bit (same entries, same order)
    assert_eq!(atw, atw2);
    let (l, r) = (dot(&av, &hw), dot(&hv, &atw));
    let size: f64 = av.iter().zip(&hw).map(|(a, b)| (a * b).abs()).sum();
    assert!(
        (l - r).abs() <= 1e-14 * size,
        "<Av, w> = {l:e}, <v, A^T w> = {r:e}"
    );
    t
}

#[test]
fn transposed_matrices_of_a_simplec_solve() {
    let mut fl = cavity(None);
    let (val, pval) = (fl.w.val, fl.w.pval);
    // momentum: convection makes A non-symmetric
    let t = transposed(&mut fl, val);
    assert_ne!(fl.dev.download(t).unwrap(), fl.dev.download(val).unwrap());
    // pressure correction: symmetric, Aᵀ = A bit for bit
    let tp = transposed(&mut fl, pval);
    assert_eq!(fl.dev.download(tp).unwrap(), fl.dev.download(pval).unwrap());
}

#[test]
fn transposed_products_do_not_depend_on_the_threads() {
    let run = |threads| {
        let mut fl = cavity(Some(threads));
        let nc = fl.disc.n_cells;
        let val = fl.w.val;
        let x = fl.dev.alloc(nc).unwrap();
        fl.dev.upload(x, &values(nc, 3)).unwrap();
        let out = fl.dev.alloc(nc).unwrap();
        fl.dev.transposed_product(val, x, out).unwrap();
        fl.dev.download(out).unwrap()
    };
    let a = run(1);
    let b = run(4);
    assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
}

/// The solvers of the primal run on Aᵀ unchanged: CG with the multigrid on
/// the transposed pressure-correction matrix (symmetric, the multigrid of
/// stage 3 needs symmetry), BiCGStab with Jacobi on the transposed momentum
/// matrix. The residual is checked with the product that never stores Aᵀ.
#[test]
fn the_linear_solvers_solve_the_transposed_systems() {
    let mut fl = cavity(None);
    let nc = fl.disc.n_cells;
    let amg = fl.settings.pressure.amg.clone();
    let cases = [
        (fl.w.pval, Preconditioner::Amg, true),
        (fl.w.val, Preconditioner::Jacobi, false),
    ];
    for (val, pc, cg) in cases {
        let t = transposed(&mut fl, val);
        let dev = &mut fl.dev;
        let (b, x, r) = (
            dev.alloc(nc).unwrap(),
            dev.alloc(nc).unwrap(),
            dev.alloc(nc).unwrap(),
        );
        let hb = values(nc, 4);
        dev.upload(b, &hb).unwrap();
        if pc == Preconditioner::Amg {
            dev.amg_setup(t, &amg).unwrap();
        }
        let (it, r0) = dev.solve(t, b, x, pc, cg, 1e-10, 500).unwrap();
        assert!(it > 0 && it < 500, "{it} iterations on A^T");
        dev.transposed_product(val, x, r).unwrap();
        let hr = dev.download(r).unwrap();
        let res: f64 = hr
            .iter()
            .zip(&hb)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        assert!(res <= 1e-9 * r0, "|b - A^T x| = {res:e}, |b| = {r0:e}");
    }
}

/// Binomial checkpointing of a transient run (D59, item 6): the reverse
/// sweep of 13 BDF2 steps with 2 snapshots meets every level with the
/// bits of the plain forward run, so a step recomputed from a snapshot is
/// the step of the run.
#[test]
fn checkpointed_reverse_sweep_meets_the_forward_levels_bitwise() {
    let (case, mesh, geo) =
        cavity_case("[time]\nscheme = \"bdf2\"\ntime_step = 0.05\nsteps = 13\n");
    let time = case.time.unwrap();
    let make = || {
        let phys = physics(&case);
        let arrays =
            bc::arrays(&mesh, &geo, &case.boundary, phys.rho, bc::At::numbers(None)).unwrap();
        let disc = simplec::discretisation(&mesh, &geo, &arrays).unwrap();
        let mut fl = Flow::on(
            disc,
            arrays,
            phys,
            case.solver.clone(),
            true,
            Some(0),
            &Target {
                device: "cpu".into(),
                threads: None,
                ptx: Default::default(),
                ranks: None,
            },
        )
        .unwrap();
        let mut tr = Transient::new(&mut fl, time.scheme, 0.0).unwrap();
        let nc = mesh.n_cells;
        fl.initialise(&foehncraft_fv::disc::Soa3::zeros(nc), &vec![0.0; nc])
            .unwrap();
        tr.start(&mut fl).unwrap();
        (fl, tr)
    };
    let n = 13;
    let fields = |fl: &mut Flow| {
        let s = fl.state().unwrap();
        let mut v = s.u.x.clone();
        v.extend(&s.u.y);
        v.extend(&s.p);
        v.extend(&s.mass_flux);
        v.iter().map(|x| x.to_bits()).collect::<Vec<u64>>()
    };
    let (mut fl, mut tr) = make();
    let mut levels = vec![fields(&mut fl)];
    for _ in 0..n {
        tr.step(&mut fl, time.time_step).unwrap();
        levels.push(fields(&mut fl));
    }
    let (mut fl, mut tr) = make();
    let mut run = FlowRun {
        fl: &mut fl,
        tr: &mut tr,
        dt: time.time_step,
    };
    let mut met = 0;
    let st = reverse(&mut run, n, 2, |r, i| {
        assert!(
            fields(r.fl) == levels[i],
            "level {i} differs from the forward run"
        );
        assert_eq!(r.tr.clock.step, i);
        met += 1;
        Ok(())
    })
    .unwrap();
    assert_eq!(met, n);
    assert!(st.max_snapshots <= 2);
    // 13 steps with 2 snapshots: t = 3, 3·13 − β(4, 2) = 24 primal steps
    assert_eq!(st.steps, 24);
}

/// D65: a step taken back after its error estimate (`Transient::reject`)
/// leaves every level as before, so that the repeated step and the steps
/// after it are bitwise those of a run that took the shorter step at once;
/// for BDF2 and Crank–Nicolson.
#[test]
fn rejected_step_repeats_bitwise() {
    use crate::transient::{ErrorControl, TimeScheme};
    for scheme in [TimeScheme::Bdf2, TimeScheme::CrankNicolson] {
        let (case, mesh, geo) =
            cavity_case("[time]\nscheme = \"bdf2\"\ntime_step = 0.05\nsteps = 9\n");
        let make = || {
            let phys = physics(&case);
            let arrays =
                bc::arrays(&mesh, &geo, &case.boundary, phys.rho, bc::At::numbers(None)).unwrap();
            let disc = simplec::discretisation(&mesh, &geo, &arrays).unwrap();
            let mut fl = Flow::on(
                disc,
                arrays,
                phys,
                case.solver.clone(),
                true,
                Some(0),
                &Target {
                    device: "cpu".into(),
                    threads: None,
                    ptx: Default::default(),
                    ranks: None,
                },
            )
            .unwrap();
            let mut tr = Transient::new(&mut fl, scheme, 0.0).unwrap();
            tr.error_control = Some(ErrorControl {
                relative: 1e-3,
                absolute: 1e-6,
            });
            let nc = mesh.n_cells;
            fl.initialise(&foehncraft_fv::disc::Soa3::zeros(nc), &vec![0.0; nc])
                .unwrap();
            tr.start(&mut fl).unwrap();
            (fl, tr)
        };
        let fields = |fl: &mut Flow| {
            let s = fl.state().unwrap();
            let mut v = s.u.x.clone();
            v.extend(&s.u.y);
            v.extend(&s.p);
            v.extend(&s.mass_flux);
            v.extend(&s.grad_u[0].x);
            v.extend(&s.psi.y);
            v.iter().map(|x| x.to_bits()).collect::<Vec<u64>>()
        };
        let steps = [0.05, 0.06, 0.07, 0.08, 0.05];
        let (short, after) = (0.03, 0.04);
        // the run that takes the long step, rejects it, repeats it short
        let (mut fl, mut tr) = make();
        for &h in &steps {
            tr.step(&mut fl, h).unwrap();
        }
        fl.keep_valid().unwrap();
        let r = tr.step(&mut fl, 0.2).unwrap();
        assert!(r.error.is_some(), "no error estimate in step {}", r.step);
        tr.reject(&mut fl).unwrap();
        assert_eq!(tr.clock.step, steps.len());
        fl.keep_valid().unwrap();
        let a = tr.step(&mut fl, short).unwrap();
        let fa = fields(&mut fl);
        let a2 = tr.step(&mut fl, after).unwrap();
        let fa2 = fields(&mut fl);
        // the run that takes the short step at once
        let (mut fl, mut tr) = make();
        for &h in &steps {
            tr.step(&mut fl, h).unwrap();
        }
        fl.keep_valid().unwrap();
        let b = tr.step(&mut fl, short).unwrap();
        assert!(
            fields(&mut fl) == fa,
            "{scheme:?}: the repeated step differs"
        );
        let b2 = tr.step(&mut fl, after).unwrap();
        assert!(fields(&mut fl) == fa2, "{scheme:?}: the step after differs");
        assert_eq!(a.error.map(f64::to_bits), b.error.map(f64::to_bits));
        assert_eq!(a2.error.map(f64::to_bits), b2.error.map(f64::to_bits));
        assert_eq!(a2.time.to_bits(), b2.time.to_bits());
    }
}

/// D52: the potential-flow start in a channel with a uniform inflow is the
/// uniform flow (exact for the linear potential on hexahedra, to the
/// discretisation on perturbed polygons), its face fluxes are continuous,
/// and the run after it builds its own multigrid.
#[test]
fn potential_flow_start_in_a_channel() {
    for (kind, tol) in [("hex", 1e-9), ("poly", 0.025)] {
        let text = format!(
            r#"
[mesh.generate]
type = "{kind}"
dimension = 2
upper = [2.0, 1.0, 0.1]
cells = [24, 12, 1]
perturb = {}
[fluid]
viscosity = 0.1
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
"#,
            if kind == "poly" { 0.15 } else { 0.0 }
        );
        let case = FlowCase::from_toml(&text, Path::new(".")).unwrap();
        let (mesh, geo, _) = mesh_for(&case, None).unwrap();
        let phys = physics(&case);
        let arrays =
            bc::arrays(&mesh, &geo, &case.boundary, phys.rho, bc::At::numbers(None)).unwrap();
        let disc = simplec::discretisation(&mesh, &geo, &arrays).unwrap();
        let mut fl = Flow::new(
            disc,
            arrays,
            phys,
            case.solver.clone(),
            [0.0; 3],
            0.0,
            None,
            &Target::default(),
        )
        .unwrap();
        let it = fl.potential_flow(&mesh, &geo).unwrap();
        assert!(it > 0 && fl.amg_setups == 0);
        let s = fl.state().unwrap();
        let nc = mesh.n_cells;
        let (mut worst, mut sq) = (0.0f64, 0.0);
        for c in 0..nc {
            let e = ((s.u.x[c] - 1.0).powi(2) + s.u.y[c].powi(2) + s.u.z[c].powi(2)).sqrt();
            worst = worst.max(e);
            sq += e * e;
        }
        let rms = (sq / nc as f64).sqrt();
        eprintln!("{kind}: deviation from the uniform flow largest {worst:.3e}, RMS {rms:.3e}");
        assert!(rms < tol, "{kind}: RMS deviation {rms:e}");
        // continuity of the face fluxes: Σ ±ṁ per cell against the inflow
        let mut div = vec![0.0; nc];
        for (f, &m) in s.mass_flux.iter().enumerate() {
            div[mesh.owner[f] as usize] += m;
            if f < mesh.n_interior_faces() {
                div[mesh.neighbour[f] as usize] -= m;
            }
        }
        let inflow = 1.0 * 1.0 * 0.1;
        let d = div.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        assert!(d < 1e-8 * inflow, "{kind}: cell imbalance {d:e}");
        fl.iterate().unwrap();
        assert_eq!(fl.amg_setups, 1);
    }
}

// ---------------------------------------------------------------------
// L17 stage A0: the steady residual R(q) as one callable (residual.rs).
// ---------------------------------------------------------------------

/// The perturbed-polygon cavity of `cavity_case` with `extra` settings,
/// iterated `n` outer iterations.
fn cavity_run(extra: &str, n: usize) -> Flow {
    let (case, mesh, geo) = cavity_case(extra);
    let phys = physics(&case);
    let arrays = bc::arrays(&mesh, &geo, &case.boundary, phys.rho, bc::At::numbers(None)).unwrap();
    let disc = simplec::discretisation(&mesh, &geo, &arrays).unwrap();
    let mut fl = Flow::new(
        disc,
        arrays,
        phys,
        case.solver.clone(),
        [0.0; 3],
        0.0,
        Some(0),
        &Target {
            device: "cpu".into(),
            threads: Some(1),
            ptx: Default::default(),
            ranks: None,
        },
    )
    .unwrap();
    for _ in 0..n {
        fl.iterate().unwrap();
    }
    fl
}

const SIMPLEC: &str = "[solver.momentum]\nrelaxation = 0.8\n[solver.pressure]\ntolerance = 1e-10\n";
const COUPLED: &str = "[solver]\nalgorithm = \"coupled\"\n[solver.pressure]\ntolerance = 1e-10\n";

/// SIMPLEC and the coupled solver converge to states where the one
/// residual R(q) vanishes (to the convergence of the run), and R of an
/// early iterate is far from zero: the callable is the residual both
/// algorithms solve.
#[test]
fn steady_residual_vanishes_for_both_algorithms() {
    let early = cavity_run(SIMPLEC, 5).steady_residual().unwrap().norms();
    println!(
        "SIMPLEC after 5 iterations: |R_u| {:.2e}, |R_c| {:.2e}, |R_m| {:.2e}",
        early[0], early[1], early[2]
    );
    assert!(early[0] > 1e-3 && early[2] > 1e-3);
    for (name, extra, n) in [("SIMPLEC", SIMPLEC, 600), ("coupled", COUPLED, 150)] {
        let mut fl = cavity_run(extra, n);
        let r = fl.steady_residual().unwrap().norms();
        println!(
            "{name} after {n} iterations: |R_u| {:.2e}, |R_c| {:.2e}, |R_m| {:.2e}",
            r[0], r[1], r[2]
        );
        for (a, b) in r.iter().zip(&early) {
            assert!(
                *a < 1e-6 && *a < 1e-4 * b,
                "{name}: {r:?} against {early:?}"
            );
        }
    }
}

/// R_ṁ = (ṁ − ṁ*)/α_u does not depend on the relaxation: the same state
/// gives the same residual for α_u = 0.8 and 0.5 (to round-off).
#[test]
fn steady_residual_is_independent_of_the_relaxation() {
    let mut fl = cavity_run(SIMPLEC, 20);
    let a = fl.steady_residual().unwrap();
    fl.settings.momentum.relaxation = Some(0.5);
    let b = fl.steady_residual().unwrap();
    let rel = |x: &[f64], y: &[f64]| {
        let d: f64 = x
            .iter()
            .zip(y)
            .map(|(p, q)| (p - q).powi(2))
            .sum::<f64>()
            .sqrt();
        d / x.iter().map(|p| p * p).sum::<f64>().sqrt()
    };
    let e = (
        rel(&a.mass_flux, &b.mass_flux),
        rel(&a.momentum[0], &b.momentum[0]),
    );
    println!("relative change of R_m {:.1e}, of R_u,x {:.1e}", e.0, e.1);
    assert!(e.0 < 1e-10 && e.1 < 1e-10, "{e:?}");
}

/// A residual call between two outer iterations leaves the run unchanged,
/// bit for bit (the output-schedule rule of AGENTS.md).
#[test]
fn steady_residual_does_not_change_the_run() {
    for extra in [SIMPLEC, COUPLED] {
        let mut a = cavity_run(extra, 8);
        let mut b = cavity_run(extra, 8);
        b.steady_residual().unwrap();
        for _ in 0..4 {
            a.iterate().unwrap();
            b.iterate().unwrap();
        }
        let (ua, ub) = (
            a.dev.download(a.f.u.x).unwrap(),
            b.dev.download(b.f.u.x).unwrap(),
        );
        let (pa, pb) = (
            a.dev.download(a.f.p).unwrap(),
            b.dev.download(b.f.p).unwrap(),
        );
        assert!(ua == ub && pa == pb, "a residual call changed the run");
    }
}
