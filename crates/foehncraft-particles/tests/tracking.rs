//! Acceptance of the parcel layer (SP1, D240 to D246): the parcel against an
//! independent integration of its equation of motion, the contract of the
//! tracker (every parcel in the cell it is assigned to), the mass books of
//! the injectors and the walls, and bitwise equal results for any number of
//! threads.

use foehncraft_kernels::particle::drag_correction;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, Env, Gas, GasState, Injector, Spray, SpraySection,
};

const RHO_G: f64 = 1.2;
const MU_G: f64 = 1.8e-5;
const RHO_P: f64 = 1000.0;

struct Field<F: Fn(usize, [f64; 3]) -> [f64; 3] + Sync>(F);

impl<F: Fn(usize, [f64; 3]) -> [f64; 3] + Sync> Gas for Field<F> {
    fn at(&self, cell: usize, x: [f64; 3]) -> GasState {
        GasState {
            u: (self.0)(cell, x),
            rho: RHO_G,
            mu: MU_G,
        }
    }
}

fn rules(n: usize, mode: BoundaryMode) -> Vec<BoundaryRule> {
    vec![
        BoundaryRule {
            mode,
            e_normal: 1.0,
            e_tangential: 1.0,
        };
        n
    ]
}

fn injector(position: [f64; 3], direction: [f64; 3], speed: f64, d: f64) -> Injector {
    Injector {
        name: "nozzle".into(),
        position,
        direction,
        speed,
        half_angle: 0.0,
        mass_flow: 1e-3,
        diameter: d,
        size: None,
        parcels_per_second: 50.0,
        start: 0.0,
        end: None,
    }
}

/// dv/dt of the equation of motion, the same model as the kernel (drag of
/// Schiller and Naumann, gravity, buoyancy), written out again here.
fn accel(v: [f64; 3], u: [f64; 3], d: f64, g: [f64; 3]) -> [f64; 3] {
    let rel = ((u[0] - v[0]).powi(2) + (u[1] - v[1]).powi(2) + (u[2] - v[2]).powi(2)).sqrt();
    let re = RHO_G * d * rel / MU_G;
    let tau = RHO_P * d * d / (18.0 * MU_G * drag_correction(re));
    let mut a = [0.0; 3];
    for k in 0..3 {
        a[k] = (u[k] - v[k]) / tau + g[k] * (1.0 - RHO_G / RHO_P);
    }
    a
}

/// Reference: RK4 with small steps, uniform gas.
fn rk4(
    x0: [f64; 3],
    v0: [f64; 3],
    u: [f64; 3],
    d: f64,
    g: [f64; 3],
    t: f64,
) -> ([f64; 3], [f64; 3]) {
    let n = 20000;
    let h = t / n as f64;
    let (mut x, mut v) = (x0, v0);
    for _ in 0..n {
        let f = |v: [f64; 3]| accel(v, u, d, g);
        let k1v = f(v);
        let k1x = v;
        let v2 = [
            v[0] + 0.5 * h * k1v[0],
            v[1] + 0.5 * h * k1v[1],
            v[2] + 0.5 * h * k1v[2],
        ];
        let k2v = f(v2);
        let v3 = [
            v[0] + 0.5 * h * k2v[0],
            v[1] + 0.5 * h * k2v[1],
            v[2] + 0.5 * h * k2v[2],
        ];
        let k3v = f(v3);
        let v4 = [v[0] + h * k3v[0], v[1] + h * k3v[1], v[2] + h * k3v[2]];
        let k4v = f(v4);
        for k in 0..3 {
            x[k] += h / 6.0 * (k1x[k] + 2.0 * v2[k] + 2.0 * v3[k] + v4[k]);
            v[k] += h / 6.0 * (k1v[k] + 2.0 * k2v[k] + 2.0 * k3v[k] + k4v[k]);
        }
    }
    (x, v)
}

fn mesh_of(kind: Kind, dim: u8, cells: [usize; 3], upper: [f64; 3]) -> foehncraft_mesh::Mesh {
    let mut s = BoxSpec::new(kind, dim, cells);
    s.lower = [0.0; 3];
    s.upper = upper;
    s.perturb = 0.2;
    s.generate().unwrap()
}

/// A droplet released into a uniform stream with gravity follows the
/// reference integration of the same equation, on every mesh type.
#[test]
fn parcel_follows_the_reference_integration() {
    let u = [0.3, -0.1, 0.05];
    let g = [0.0, 0.0, -9.81];
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        let m = mesh_of(kind, 3, [6, 6, 6], [1.0, 1.0, 1.0]);
        let geo = Geometry::compute(&m);
        let loc = Locator::new(&m, &geo);
        let env = Env {
            mesh: &m,
            geo: &geo,
            loc: &loc,
        };
        // 40 µm: tau ≈ 0.12 s; Re of order 1
        let d = 40e-6;
        let mut inj = injector([0.2, 0.5, 0.9], [1.0, 0.0, 0.0], 0.5, d);
        inj.parcels_per_second = 1.0; // one parcel at t = 0.5
        let cfg = SpraySection {
            injectors: vec![inj],
            ..SpraySection::default()
        };
        let mut sp = Spray::new(
            &cfg,
            &env,
            rules(m.boundaries.len(), BoundaryMode::Escape),
            g,
        )
        .unwrap();
        let gas = Field(|_, _| u);
        let dt = 0.01;
        let mut t = 0.0;
        // the parcel is born at 0.5 s: advance to 0.5 + 0.2 s
        for _ in 0..70 {
            sp.advance(&env, &gas, t, dt).unwrap();
            t += dt;
        }
        assert_eq!(
            sp.parcels.len(),
            1,
            "{kind:?}: the parcel must still be in the box"
        );
        let (xr, vr) = rk4([0.2, 0.5, 0.9], [0.5, 0.0, 0.0], u, d, g, t - 0.5);
        let x = sp.parcels.x[0];
        let v = sp.parcels.v[0];
        let len = (xr[0] - 0.2).hypot(xr[2] - 0.9);
        for k in 0..3 {
            assert!(
                (x[k] - xr[k]).abs() < 2e-3 * len.max(1e-3),
                "{kind:?}: x[{k}] = {} against {}",
                x[k],
                xr[k]
            );
            assert!(
                (v[k] - vr[k]).abs() < 2e-3 * (vr[0].abs() + vr[2].abs()),
                "{kind:?}: v[{k}] = {} against {}",
                v[k],
                vr[k]
            );
        }
        // the cell holds the parcel
        let c = sp.parcels.cell[0] as usize;
        assert!(loc.contains(c, x), "{kind:?}: parcel outside its cell");
    }
}

/// Mass books: released = in the flow + stuck + escaped + lost, with every
/// wall mode, and the tracker contract holds for every parcel.
#[test]
fn mass_books_close_and_parcels_stay_in_their_cells() {
    for mode in [
        BoundaryMode::Stick,
        BoundaryMode::Rebound,
        BoundaryMode::Escape,
    ] {
        for kind in [Kind::Hex, Kind::Tet, Kind::Poly] {
            let m = mesh_of(kind, 3, [5, 5, 5], [1.0, 1.0, 1.0]);
            let geo = Geometry::compute(&m);
            let loc = Locator::new(&m, &geo);
            let env = Env {
                mesh: &m,
                geo: &geo,
                loc: &loc,
            };
            let mut cfg = SpraySection::default();
            let mut inj = injector([0.5, 0.5, 0.5], [1.0, 0.3, 0.2], 3.0, 150e-6);
            inj.half_angle = 25.0;
            inj.parcels_per_second = 120.0;
            cfg.injectors = vec![inj];
            let mut sp = Spray::new(
                &cfg,
                &env,
                rules(m.boundaries.len(), mode),
                [0.0, 0.0, -9.81],
            )
            .unwrap();
            // a swirling gas keeps the parcels moving through many cells
            let gas = Field(|_, x| [-(x[1] - 0.5), x[0] - 0.5, 0.0]);
            let dt = 0.01;
            for s in 0..250 {
                sp.advance(&env, &gas, s as f64 * dt, dt).unwrap();
                for i in 0..sp.parcels.len() {
                    let c = sp.parcels.cell[i] as usize;
                    assert!(
                        loc.contains(c, sp.parcels.x[i]),
                        "{mode:?} {kind:?} step {s}: parcel {i} is outside cell {c}"
                    );
                }
            }
            let t = &sp.tally;
            let stuck: f64 = t.deposited_mass.iter().sum();
            let gone: f64 = t.escaped_mass.iter().sum();
            let total = sp.mass_in_flow() + stuck + gone + t.lost_mass;
            assert!(
                (t.injected_mass - total).abs() < 1e-12 * t.injected_mass,
                "{mode:?} {kind:?}: released {} against {}",
                t.injected_mass,
                total
            );
            assert_eq!(t.lost_parcels, 0, "{mode:?} {kind:?}: parcels were lost");
            match mode {
                BoundaryMode::Stick => assert!(stuck > 0.0),
                BoundaryMode::Escape => assert!(gone > 0.0),
                _ => {
                    assert!(t.rebounds > 0);
                    assert_eq!(stuck + gone, 0.0);
                }
            }
        }
    }
}

/// A droplet that bounces elastically in a closed box keeps its speed
/// relative to still gas only up to the drag: with a very large droplet and
/// no gravity the speed is nearly constant, and the parcel never leaves.
#[test]
fn elastic_bounce_keeps_the_parcel_inside() {
    let m = mesh_of(Kind::Hex, 3, [4, 4, 4], [1.0, 1.0, 1.0]);
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let mut cfg = SpraySection::default();
    let mut inj = injector([0.3, 0.4, 0.5], [1.0, 0.7, 0.4], 2.0, 3e-3);
    inj.parcels_per_second = 1.0;
    inj.end = Some(1.0);
    cfg.injectors = vec![inj];
    let mut sp = Spray::new(
        &cfg,
        &env,
        rules(m.boundaries.len(), BoundaryMode::Rebound),
        [0.0; 3],
    )
    .unwrap();
    let gas = Field(|_, _| [0.0; 3]);
    for s in 0..400 {
        sp.advance(&env, &gas, s as f64 * 0.005, 0.005).unwrap();
    }
    assert_eq!(sp.parcels.len(), 1);
    assert!(sp.tally.rebounds >= 3, "rebounds: {}", sp.tally.rebounds);
    let v = sp.parcels.v[0];
    let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    assert!(speed > 0.0 && speed < 2.0);
    for xk in &sp.parcels.x[0] {
        assert!((0.0..=1.0).contains(xk));
    }
}

/// The same bits for one, two and four threads.
#[test]
fn threads_do_not_change_a_bit() {
    let m = mesh_of(Kind::Tet, 3, [5, 5, 5], [1.0, 1.0, 1.0]);
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let run = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let mut cfg = SpraySection::default();
            let mut inj = injector([0.5, 0.5, 0.2], [0.0, 0.2, 1.0], 4.0, 60e-6);
            inj.half_angle = 30.0;
            inj.parcels_per_second = 800.0;
            cfg.injectors = vec![inj];
            let mut sp = Spray::new(
                &cfg,
                &env,
                rules(m.boundaries.len(), BoundaryMode::Rebound),
                [0.0, 0.0, -9.81],
            )
            .unwrap();
            let gas = Field(|_, x| [0.3 * (x[1] - 0.5), -0.3 * (x[0] - 0.5), 0.1]);
            for s in 0..60 {
                sp.advance(&env, &gas, s as f64 * 0.01, 0.01).unwrap();
            }
            (sp.state().unwrap(), sp.parcels.len())
        })
    };
    let one = run(1);
    assert!(one.1 > 50);
    assert!(one == run(2), "2 threads differ");
    assert!(one == run(4), "4 threads differ");
}

/// A restart gives the bits of the uninterrupted run.
#[test]
fn restart_state_continues_bitwise() {
    let m = mesh_of(Kind::Poly, 3, [5, 5, 5], [1.0, 1.0, 1.0]);
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let mut cfg = SpraySection::default();
    let mut inj = injector([0.5, 0.5, 0.2], [0.0, 0.2, 1.0], 4.0, 60e-6);
    inj.half_angle = 30.0;
    inj.parcels_per_second = 500.0;
    cfg.injectors = vec![inj];
    let make = || {
        Spray::new(
            &cfg,
            &env,
            rules(m.boundaries.len(), BoundaryMode::Rebound),
            [0.0, 0.0, -9.81],
        )
        .unwrap()
    };
    let gas = Field(|_, x| [0.3 * (x[1] - 0.5), -0.3 * (x[0] - 0.5), 0.1]);
    let mut a = make();
    for s in 0..40 {
        a.advance(&env, &gas, s as f64 * 0.01, 0.01).unwrap();
    }
    let mut b = make();
    for s in 0..20 {
        b.advance(&env, &gas, s as f64 * 0.01, 0.01).unwrap();
    }
    let saved = b.state().unwrap();
    let mut c = make();
    c.restore(&env, &saved).unwrap();
    for s in 20..40 {
        c.advance(&env, &gas, s as f64 * 0.01, 0.01).unwrap();
    }
    assert!(a.state().unwrap() == c.state().unwrap());
}

/// The injector outside the mesh is refused.
#[test]
fn injector_outside_the_mesh_is_refused() {
    let m = mesh_of(Kind::Hex, 3, [3, 3, 3], [1.0, 1.0, 1.0]);
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let cfg = SpraySection {
        injectors: vec![injector([3.0, 0.5, 0.5], [1.0, 0.0, 0.0], 1.0, 1e-4)],
        ..SpraySection::default()
    };
    let r = Spray::new(
        &cfg,
        &env,
        rules(m.boundaries.len(), BoundaryMode::Stick),
        [0.0; 3],
    );
    assert!(r.err().unwrap().to_string().contains("outside the mesh"));
}
