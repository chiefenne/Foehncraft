//! Breakup of the drops (SP6, D300 to D307) in a gas at rest that the test
//! holds:
//!
//! - blobs shed child parcels whose drops have the KH child size, and the
//!   liquid is kept to round-off at every split;
//! - the results are the same bits on one and on four threads with the
//!   parcels stored in another order (D247), and a restart continues with
//!   the same bits;
//! - the KH radius of a blob follows the breakup equation da/dt =
//!   −(a − B₀Λ(a))/τ(a), integrated independently;
//! - Rayleigh-Taylor breakup acts beyond the breakup length only;
//! - the dynamic drag of distorted drops slows the spray down;
//! - a rate table gives each parcel its exact share of the liquid and the
//!   speed ṁ/(ρ_l C_a A_hole) (D306);
//! - KH-ACT (D307): the nozzle turbulence strips a young blob at
//!   C_T,CAV √k exactly, the liquid and the bits are kept, and the
//!   nozzle keys are checked.

use foehncraft_kernels::breakup as brk;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::spray::child_id;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, BreakupModel, DragLaw, Env, Gas, GasState, Injector, Parcels,
    Spray, SpraySection,
};

const RHO_G: f64 = 20.0;
const MU_G: f64 = 1.8e-5;
const RHO_L: f64 = 700.0;
const SIGMA: f64 = 0.02;
const MU_L: f64 = 5e-4;

/// Gas at rest.
struct Still(f64);

impl Gas for Still {
    fn at(&self, _: usize, _: [f64; 3]) -> GasState {
        GasState::isothermal([0.0; 3], self.0, MU_G)
    }
}

fn channel() -> (foehncraft_mesh::Mesh, Geometry) {
    let mut s = BoxSpec::new(Kind::Hex, 3, [40, 4, 4]);
    s.lower = [0.0; 3];
    s.upper = [0.08, 0.02, 0.02];
    let m = s.generate().unwrap();
    let g = Geometry::compute(&m);
    (m, g)
}

fn cfg(rho_l: f64) -> SpraySection {
    let mut c = SpraySection {
        injectors: vec![Injector {
            name: "hole".into(),
            position: [0.002, 0.01, 0.01],
            direction: [1.0, 0.0, 0.0],
            half_angle: 5.0,
            mass_flow: 1e-3,
            parcels_per_second: 4e5,
            hole_diameter: Some(1e-4),
            area_contraction: 0.8,
            ..Default::default()
        }],
        ..SpraySection::default()
    };
    c.liquid.density = rho_l;
    c.liquid.surface_tension = Some(SIGMA);
    c.liquid.viscosity = Some(MU_L);
    c.breakup.model = BreakupModel::KhRt;
    c
}

fn spray(env: &Env, c: &SpraySection) -> Spray {
    let rules = vec![
        BoundaryRule {
            mode: BoundaryMode::Escape,
            e_normal: 1.0,
            e_tangential: 1.0,
        };
        env.mesh.boundaries.len()
    ];
    Spray::new(c, env, rules, [0.0; 3]).unwrap()
}

fn reversed(p: &Parcels) -> Parcels {
    let mut q = p.clone();
    q.id.reverse();
    q.x.reverse();
    q.v.reverse();
    q.d.reverse();
    q.count.reverse();
    q.born.reverse();
    q.temp.reverse();
    q.eddy.reverse();
    q.brk.reverse();
    q.cell.reverse();
    q
}

fn mass(p: &Parcels, rho_l: f64) -> f64 {
    let mut ix: Vec<usize> = (0..p.len()).collect();
    ix.sort_by_key(|&i| p.id[i]);
    ix.iter()
        .map(|&i| p.count[i] * rho_l * std::f64::consts::PI / 6.0 * p.d[i].powi(3))
        .sum()
}

/// The parcels sorted by identity: (id, x, v, d, count).
type Row = (u64, [f64; 3], [f64; 3], f64, f64);

fn rows(p: &Parcels) -> Vec<Row> {
    let mut r: Vec<Row> = (0..p.len())
        .map(|i| (p.id[i], p.x[i], p.v[i], p.d[i], p.count[i]))
        .collect();
    r.sort_by_key(|x| x.0);
    r
}

fn run(threads: usize, shuffle: bool, steps: usize, c: &SpraySection) -> (Spray, Vec<Row>) {
    let (m, geo) = channel();
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let mut sp = spray(&env, c);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let dt = 5e-6;
    for n in 0..steps {
        if shuffle {
            sp.parcels = reversed(&sp.parcels);
        }
        pool.install(|| sp.advance(&env, &Still(RHO_G), n as f64 * dt, dt).unwrap());
    }
    let r = rows(&sp.parcels);
    (sp, r)
}

#[test]
fn blobs_shed_children_and_keep_the_liquid() {
    let c = cfg(RHO_L);
    let (sp, _) = run(1, false, 40, &c);
    let p = &sp.parcels;
    let kids: Vec<usize> = (0..p.len()).filter(|&i| p.id[i] & (1 << 47) != 0).collect();
    let blobs: Vec<usize> = (0..p.len()).filter(|&i| p.brk[i].blob).collect();
    assert!(kids.len() > 10, "{} children", kids.len());
    assert!(!blobs.is_empty());
    // every child shed is still in the flow (the gas is at rest and the run
    // short)
    assert_eq!(sp.tally.child_parcels as usize, kids.len());
    // the children are much smaller than the blob (√0.8 × 100 µm)
    let blob_d = 0.8f64.sqrt() * 1e-4;
    for &i in &kids {
        assert!(p.d[i] < 0.3 * blob_d, "child of {} m", p.d[i]);
        assert!(p.count[i] > 1.0);
        // the identity names the parent's injection index
        let k = p.id[i] & 0xFFFF_FFFF;
        let j = ((p.id[i] >> 32) & 0x7FFF) as u32;
        assert_eq!(child_id(k, j), p.id[i]);
    }
    // the liquid: in the flow plus out of the domain is what was injected
    let out: f64 = sp.tally.escaped_mass.iter().sum::<f64>() + sp.tally.lost_mass;
    let inflow = mass(p, RHO_L) + out;
    let rel = (inflow - sp.tally.injected_mass).abs() / sp.tally.injected_mass;
    assert!(rel < 1e-13, "liquid off by {rel}");
    // every child carries about shed_fraction of its blob's liquid (between
    // one and two shares when the sub-steps end on the share)
    for &i in &kids {
        let m0 = p.brk[i].m0;
        let parent_m0 = 1e-3 / 4e5;
        assert!(
            m0 >= 0.03 * parent_m0 * (1.0 - 1e-9),
            "child share {}",
            m0 / parent_m0
        );
        assert!(m0 <= 0.09 * parent_m0, "child share {}", m0 / parent_m0);
    }
}

#[test]
fn same_bits_on_any_threads_and_order() {
    let mut c = cfg(RHO_L);
    c.breakup.drag = DragLaw::Dynamic;
    c.breakup.breakup_length = Some(2.0);
    let (a, ra) = run(1, false, 30, &c);
    let (_, rb) = run(4, true, 30, &c);
    assert!(a.tally.child_parcels > 0 && a.tally.rt_breakups > 0);
    assert_eq!(ra.len(), rb.len());
    for (x, y) in ra.iter().zip(&rb) {
        assert_eq!(x.0, y.0);
        for k in 0..3 {
            assert_eq!(x.1[k].to_bits(), y.1[k].to_bits());
            assert_eq!(x.2[k].to_bits(), y.2[k].to_bits());
        }
        assert_eq!(x.3.to_bits(), y.3.to_bits());
        assert_eq!(x.4.to_bits(), y.4.to_bits());
    }
}

#[test]
fn a_restart_continues_with_the_same_bits() {
    let mut c = cfg(RHO_L);
    c.breakup.drag = DragLaw::Dynamic;
    let (m, geo) = channel();
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let dt = 5e-6;
    let mut a = spray(&env, &c);
    for n in 0..15 {
        a.advance(&env, &Still(RHO_G), n as f64 * dt, dt).unwrap();
    }
    let saved = a.state().unwrap();
    for n in 15..30 {
        a.advance(&env, &Still(RHO_G), n as f64 * dt, dt).unwrap();
    }
    let mut b = spray(&env, &c);
    b.restore(&env, &saved).unwrap();
    for n in 15..30 {
        b.advance(&env, &Still(RHO_G), n as f64 * dt, dt).unwrap();
    }
    assert_eq!(rows(&a.parcels), rows(&b.parcels));
    assert_eq!(a.parcels.brk.len(), a.parcels.len());
    for (x, y) in a.parcels.brk.iter().zip(&b.parcels.brk) {
        assert_eq!(x, y);
    }
}

/// A blob of a very heavy liquid slows down little in 0.2 ms: its KH
/// radius (the radius it has stripped to, children included) follows
/// da/dt = −(a − B₀Λ(a))/τ(a) with the slip u(t) of the blob, which the
/// test integrates independently by fourth-order Runge-Kutta (the slip
/// interpolated linearly between the step ends).
#[test]
fn a_blob_strips_as_the_breakup_equation() {
    let rho_l = 1e6;
    let mut c = cfg(rho_l);
    c.injectors[0].speed = 150.0;
    c.injectors[0].half_angle = 0.0;
    c.injectors[0].parcels_per_second = 1e6;
    c.injectors[0].end = Some(1e-6); // one parcel, born at 0.5 µs
    c.breakup.b1 = 2.0;
    c.breakup.breakup_length = Some(1e6); // no RT
    let (m, geo) = channel();
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let mut sp = spray(&env, &c);
    let dt = 2e-6;
    let blob_a0 = 0.5 * 0.8f64.sqrt() * 1e-4;
    let mut hist = vec![(0.5e-6, 150.0)];
    let mut a = blob_a0;
    let mut checked = 0;
    let mut stripped = 0.0;
    for n in 0..100 {
        sp.advance(&env, &Still(RHO_G), n as f64 * dt, dt).unwrap();
        let t = (n + 1) as f64 * dt;
        let p = &sp.parcels;
        let Some(i) = (0..p.len()).find(|&i| p.id[i] & (1 << 47) == 0) else {
            break;
        };
        let speed = p.v[i][0];
        let (t0, u0) = *hist.last().unwrap();
        hist.push((t, speed));
        let a_kh = 0.5 * p.d[i] * p.brk[i].kh;
        // reference: RK4 of the breakup equation from t0 to t
        let f = |a: f64, tt: f64| {
            let u = u0 + (speed - u0) * (tt - t0) / (t - t0);
            let w = brk::kh_wave(&brk::DropInput {
                a,
                u,
                rho_l,
                rho_g: RHO_G,
                sigma: SIGMA,
                mu_l: MU_L,
            });
            let rc = 0.61 * w.lambda;
            if rc >= a {
                0.0
            } else {
                -(a - rc) / brk::kh_time(a, &w, 2.0)
            }
        };
        let nsub = 2000;
        let h = (t - t0) / nsub as f64;
        for j in 0..nsub {
            let tt = t0 + j as f64 * h;
            let k1 = f(a, tt);
            let k2 = f(a + 0.5 * h * k1, tt + 0.5 * h);
            let k3 = f(a + 0.5 * h * k2, tt + 0.5 * h);
            let k4 = f(a + h * k3, tt + h);
            a += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        let err = (a_kh - a).abs() / (blob_a0 - a).max(1e-12);
        assert!(
            err < 0.02,
            "t {t}: KH radius {a_kh} vs {a} (stripped part off by {err})"
        );
        stripped = 1.0 - (a / blob_a0).powi(3);
        checked += 1;
    }
    assert!(checked > 20);
    // the test covers a good part of the blob's life
    assert!(stripped > 0.5, "only {stripped} of the blob stripped");
}

#[test]
fn rayleigh_taylor_acts_beyond_the_breakup_length_only() {
    let mut c = cfg(RHO_L);
    // far: no RT inside the domain
    c.breakup.breakup_length = Some(1e3);
    let (far, _) = run(1, false, 30, &c);
    assert_eq!(far.tally.rt_breakups, 0);
    // at the nozzle
    c.breakup.breakup_length = Some(0.0);
    let (near, _) = run(1, false, 30, &c);
    assert!(near.tally.rt_breakups > 0);
    // RT drops are smaller than the blob
    let p = &near.parcels;
    let blob_d = 0.8f64.sqrt() * 1e-4;
    assert!(p.d.iter().all(|&d| d < blob_d));
}

#[test]
fn distorted_drops_feel_more_drag() {
    let mut c = cfg(RHO_L);
    c.breakup.model = BreakupModel::None;
    let (sphere, _) = run(1, false, 30, &c);
    c.breakup.drag = DragLaw::Dynamic;
    let (dynamic, _) = run(1, false, 30, &c);
    let far = |s: &Spray| s.parcels.x.iter().map(|x| x[0]).fold(0.0, f64::max);
    assert!(
        far(&dynamic) < far(&sphere),
        "{} vs {}",
        far(&dynamic),
        far(&sphere)
    );
    // the blobs are distorted to the limit y = 1 (We ≫ 12)
    let p = &dynamic.parcels;
    assert!(p.brk.iter().all(|b| b.y <= 1.0));
    assert!(p.brk.iter().any(|b| b.y == 1.0));
}

#[test]
fn a_rate_table_gives_exact_shares_and_speeds() {
    let mut c = cfg(RHO_L);
    c.breakup.model = BreakupModel::None;
    let inj = &mut c.injectors[0];
    inj.mass_flow = 0.0;
    inj.half_angle = 0.0;
    inj.rate = vec![[0.0, 0.0], [2e-5, 1e-3], [8e-5, 1e-3], [1e-4, 0.0]];
    inj.parcels_per_second = 3e5;
    let total = 0.5 * 2e-5 * 1e-3 + 6e-5 * 1e-3 + 0.5 * 2e-5 * 1e-3;
    let (sp, _) = run(1, false, 40, &c);
    let rel = (sp.tally.injected_mass - total).abs() / total;
    assert!(
        rel < 1e-12,
        "injected {} vs {total}",
        sp.tally.injected_mass
    );
    // the speed of a parcel of the plateau: ṁ/(ρ C_a A)
    let u = 1e-3 / (RHO_L * 0.8 * std::f64::consts::PI * 0.25 * 1e-8);
    let inj = &c.injectors[0];
    assert!((inj.exit_speed(RHO_L, 1e-3) - u).abs() < 1e-9 * u);
    // the injection ends with the table
    assert!((inj.end_time() - 1e-4).abs() < 1e-18);
    assert_eq!(sp.tally.injected_parcels, 30);
}

#[test]
fn breakup_needs_the_liquid_properties() {
    let mut c = cfg(RHO_L);
    c.liquid.surface_tension = None;
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("surface_tension"), "{e}");
    let mut c = cfg(RHO_L);
    c.liquid.viscosity = None;
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("viscosity"), "{e}");
    let mut c = cfg(RHO_L);
    c.steady = Some(Default::default());
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("steady"), "{e}");
    let mut c = cfg(RHO_L);
    c.injectors[0].rate = vec![[0.0, 1e-3], [1e-4, 1e-3]];
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("either mass_flow or rate"), "{e}");
}

/// The configuration of `cfg` with KH-ACT (D307) and the nozzle
/// turbulence estimated from the hole (Huh and Gosman), with cavities.
fn act_cfg() -> SpraySection {
    let mut c = cfg(RHO_L);
    c.breakup.model = BreakupModel::KhAct;
    c.liquid.vapour_pressure = Some(5e4);
    let inj = &mut c.injectors[0];
    inj.hole_length = Some(4e-4);
    inj.discharge_coefficient = Some(0.7);
    c
}

/// KH-ACT (D307): a young blob is stripped by the nozzle turbulence at
/// C_T,CAV √k = C_T,CAV L_T/τ_T, and the sum of the exact sub-step
/// integrals is the exact integral, whatever the sub-steps: the stripped
/// radius is a₀ − C_T,CAV ∫₀^age √k dt to round-off (Kelvin-Helmholtz made
/// slow by a large B₁, no cavity at C_a = 1, no Rayleigh-Taylor).
#[test]
fn kh_act_turbulence_strips_at_the_root_of_k() {
    let mut c = cfg(RHO_L);
    c.breakup.model = BreakupModel::KhAct;
    c.breakup.b1 = 1e6;
    c.breakup.c_tcav = 0.7;
    c.breakup.breakup_length = Some(1e6);
    let inj = &mut c.injectors[0];
    inj.speed = 150.0;
    inj.half_angle = 0.0;
    inj.parcels_per_second = 1e6;
    inj.end = Some(1e-6); // one parcel
    inj.area_contraction = 1.0;
    inj.nozzle_k = Some(25.0);
    inj.nozzle_epsilon = Some(1e7);
    c.validate().unwrap();
    let turb = brk::NozzleTurbulence {
        k0: 25.0,
        eps0: 1e7,
    };
    let (m, geo) = channel();
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let mut sp = spray(&env, &c);
    let dt = 1e-6;
    let a0 = 0.5e-4;
    let mut checked = 0;
    for n in 0..40 {
        sp.advance(&env, &Still(RHO_G), n as f64 * dt, dt).unwrap();
        let t = (n + 1) as f64 * dt;
        let p = &sp.parcels;
        let Some(i) = (0..p.len()).find(|&i| p.id[i] & (1 << 47) == 0 && p.brk[i].blob) else {
            break;
        };
        let a = 0.5 * p.d[i] * p.brk[i].kh;
        let age = t - p.born[i];
        if a <= 1.01 * turb.length(age) {
            break;
        }
        let exact = a0 - 0.7 * turb.sqrt_k_integral(0.0, age);
        assert!(
            (a - exact).abs() < 1e-12 * a0,
            "age {age}: radius {a} vs {exact}"
        );
        checked += 1;
    }
    assert!(checked > 5, "only {checked} steps checked");
    assert!(sp.tally.act_children > 0);
    assert_eq!(sp.tally.act_children, sp.tally.child_parcels);
    // every child has the integral length scale of the turbulence at the
    // age of its blob when it was shed: at least L_T0, at most the blob
    let p = &sp.parcels;
    for i in (0..p.len()).filter(|&i| p.id[i] & (1 << 47) != 0) {
        assert!(p.d[i] >= 2.0 * turb.length(0.0) * (1.0 - 1e-12));
        assert!(p.d[i] < 2.0 * a0);
    }
}

/// KH-ACT keeps the liquid to round-off and gives the same bits on any
/// number of threads and parcel order (D247), with the turbulence of the
/// Huh-Gosman estimate and the cavities.
#[test]
fn kh_act_keeps_the_liquid_and_the_bits() {
    let c = act_cfg();
    c.validate().unwrap();
    let (a, ra) = run(1, false, 30, &c);
    let (b, rb) = run(4, true, 30, &c);
    assert!(a.tally.act_children > 0);
    assert!(a.tally.child_parcels >= a.tally.act_children);
    assert_eq!(a.tally.act_children, b.tally.act_children);
    assert_eq!(ra.len(), rb.len());
    for (x, y) in ra.iter().zip(&rb) {
        assert_eq!(x.0, y.0);
        for k in 0..3 {
            assert_eq!(x.1[k].to_bits(), y.1[k].to_bits());
            assert_eq!(x.2[k].to_bits(), y.2[k].to_bits());
        }
        assert_eq!(x.3.to_bits(), y.3.to_bits());
        assert_eq!(x.4.to_bits(), y.4.to_bits());
    }
    let p = &a.parcels;
    let out: f64 = a.tally.escaped_mass.iter().sum::<f64>() + a.tally.lost_mass;
    let rel = (mass(p, RHO_L) + out - a.tally.injected_mass).abs() / a.tally.injected_mass;
    assert!(rel < 1e-13, "liquid off by {rel}");
    // the turbulence atomises the blobs faster than the KH waves alone
    let (r, _) = run(1, false, 30, &cfg(RHO_L));
    let blob_liquid = |s: &Spray| -> f64 {
        let p = &s.parcels;
        (0..p.len())
            .filter(|&i| p.brk[i].blob)
            .map(|i| p.count[i] * RHO_L * std::f64::consts::PI / 6.0 * p.d[i].powi(3))
            .sum()
    };
    assert!(
        blob_liquid(&a) < 0.8 * blob_liquid(&r),
        "blobs keep {} with KH-ACT, {} with KH-RT",
        blob_liquid(&a),
        blob_liquid(&r)
    );
}

#[test]
fn kh_act_needs_the_nozzle_turbulence() {
    let mut c = act_cfg();
    c.injectors[0].hole_length = None;
    c.injectors[0].discharge_coefficient = None;
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("turbulence of the liquid at the hole"), "{e}");
    let mut c = act_cfg();
    c.injectors[0].discharge_coefficient = None;
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("needs both hole_length"), "{e}");
    let mut c = act_cfg();
    c.injectors[0].discharge_coefficient = Some(0.85);
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("area_contraction"), "{e}");
    let mut c = act_cfg();
    c.injectors[0].area_contraction = 0.95;
    c.injectors[0].discharge_coefficient = Some(0.9);
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("no wall friction"), "{e}");
    let mut c = act_cfg();
    c.injectors[0].nozzle_k = Some(100.0);
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("together"), "{e}");
    // kh-rt leaves the nozzle turbulence unused
    let mut c = act_cfg();
    c.breakup.model = BreakupModel::KhRt;
    c.validate().unwrap();
}
