//! The parcel pass of the two-way coupling (SP2, D247 to D249, V28.1,
//! V28.6): the per-cell sums of the exchange records do not depend on the
//! storage order of the parcels or the number of threads, every sub-step
//! lies in one cell, and the slopes of a parcel add up to at most its mass;
//! with the eddies of the dispersion (SP2c) too.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, CouplingMode, Dispersion, Env, Gas, GasState, Injector, Parcels,
    Spray, SpraySection,
};

struct Swirl;

impl Gas for Swirl {
    fn at(&self, _cell: usize, x: [f64; 3]) -> GasState {
        // a hot gas whose temperature varies in space: the heat exchange
        // takes part in the order independence too (SP2b)
        GasState {
            u: [0.5 - x[1], x[0] - 0.5, 0.1],
            rho: 1.2,
            mu: 1.8e-5,
            t: 500.0 + 100.0 * x[0] - 50.0 * x[2],
            cp: 1030.0,
            k: 0.04,
            // turbulence that varies in space: the eddies of the
            // dispersion take part too (SP2c)
            tke: 0.02 + 0.03 * x[1],
            eps: 0.1 + 0.2 * x[0],
        }
    }
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
    q.cell.reverse();
    q
}

#[test]
fn sources_do_not_depend_on_parcel_order_or_threads() {
    for (kind, dispersion) in [
        (Kind::Hex, Dispersion::None),
        (Kind::Tet, Dispersion::None),
        (Kind::Poly, Dispersion::None),
        (Kind::Hex, Dispersion::EddyLifetime),
        (Kind::Tet, Dispersion::EddyLifetime),
        (Kind::Poly, Dispersion::EddyLifetime),
    ] {
        let mut s = BoxSpec::new(kind, 3, [5, 5, 5]);
        s.lower = [0.0; 3];
        s.upper = [1.0; 3];
        s.perturb = 0.2;
        let m = s.generate().unwrap();
        let geo = Geometry::compute(&m);
        let loc = Locator::new(&m, &geo);
        let env = Env {
            mesh: &m,
            geo: &geo,
            loc: &loc,
        };
        let mut cfg = SpraySection {
            injectors: vec![Injector {
                name: "a".into(),
                position: [0.3, 0.5, 0.2],
                direction: [1.0, 0.2, 0.5],
                speed: 3.0,
                half_angle: 30.0,
                mass_flow: 1e-3,
                diameter: 60e-6,
                size: None,
                parcels_per_second: 2000.0,
                start: 0.0,
                end: None,
                temperature: None,
            }],
            ..SpraySection::default()
        };
        cfg.coupling.mode = CouplingMode::TwoWay;
        cfg.dispersion = dispersion;
        let rules = vec![
            BoundaryRule {
                mode: BoundaryMode::Rebound,
                e_normal: 0.8,
                e_tangential: 0.9,
            };
            m.boundaries.len()
        ];
        let mut sp = Spray::new(&cfg, &env, rules, [0.0, 0.0, -9.81]).unwrap();
        let dt = 0.01;
        for n in 0..10 {
            sp.advance(&env, &Swirl, n as f64 * dt, dt).unwrap();
        }
        assert!(sp.parcels.len() > 100, "{kind:?}: {}", sp.parcels.len());
        let t0 = 10.0 * dt;
        let a = sp.pass(&env, &Swirl, t0, dt).unwrap();
        // every sub-step lies in its cell: the records of a parcel name the
        // cells it passes; the slopes of a parcel are at most its mass
        let mut by_parcel: std::collections::BTreeMap<u64, f64> = Default::default();
        let mut heat_by_parcel: std::collections::BTreeMap<u64, f64> = Default::default();
        for r in a.records() {
            assert!(r.slope >= 0.0 && r.dt >= 0.0 && r.heat_slope >= 0.0);
            *by_parcel.entry(r.id).or_default() += r.slope;
            *heat_by_parcel.entry(r.id).or_default() += r.heat_slope;
        }
        for (i, &id) in sp.parcels.id.iter().enumerate() {
            let mass = sp.parcel_mass(i);
            let s = by_parcel.get(&id).copied().unwrap_or(0.0);
            assert!(
                s <= mass * (1.0 + 1e-12),
                "{kind:?}: slopes {s} above the mass {mass}"
            );
            let h = heat_by_parcel.get(&id).copied().unwrap_or(0.0);
            assert!(
                h <= mass * 4181.0 * (1.0 + 1e-12),
                "{kind:?}: heat slopes {h}"
            );
        }
        // the storage order reversed, and on one thread: the same bits
        let mut other = Spray::new(&cfg, &env, sp_rules(&m), [0.0, 0.0, -9.81]).unwrap();
        other.restore(&env, &sp.state().unwrap()).unwrap();
        other.parcels = reversed(&other.parcels);
        let b = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(|| other.pass(&env, &Swirl, t0, dt).unwrap());
        assert_eq!(a.sources, b.sources, "{kind:?} {dispersion:?}");
        // and the committed parcels carry the same momentum after the
        // closing step
        let du = |c: u32| ([1e-3 * c as f64, -2e-3, 0.5e-3], 0.2 - 0.01 * c as f64);
        let ba = sp.commit(&env, a, Some(&du));
        let bb = other.commit(&env, b, Some(&du));
        assert_eq!(ba, bb, "{kind:?} {dispersion:?}");
        let mut pa: Vec<_> = (0..sp.parcels.len())
            .map(|i| (sp.parcels.id[i], sp.parcels.x[i], sp.parcels.eddy[i]))
            .collect();
        let mut pb: Vec<_> = (0..other.parcels.len())
            .map(|i| {
                (
                    other.parcels.id[i],
                    other.parcels.x[i],
                    other.parcels.eddy[i],
                )
            })
            .collect();
        pa.sort_by_key(|p| p.0);
        pb.sort_by_key(|p| p.0);
        assert_eq!(pa, pb, "{kind:?} {dispersion:?}");
        if dispersion == Dispersion::EddyLifetime {
            assert!(
                pa.iter().all(|p| p.2.n > 0),
                "{kind:?}: a parcel without eddy"
            );
        }
        assert_eq!(sp.tally.lost_parcels, 0, "{kind:?} {dispersion:?}");
        for k in 0..3 {
            assert!(
                ba.imbalance()[k].abs() < 1e-12 * (ba.end[k].abs() + ba.injected[k].abs() + 1e-30),
                "{kind:?} {dispersion:?} {k}: {:?}",
                ba
            );
        }
        // the parcels took heat from the gas, and the enthalpy books close
        assert!(ba.enthalpy[2] > 0.0, "{kind:?}: {:?}", ba.enthalpy);
        assert!(
            ba.heat_imbalance().abs() < 1e-12 * ba.enthalpy[1],
            "{kind:?}: {:?}",
            ba.enthalpy
        );
    }
}

fn sp_rules(m: &foehncraft_mesh::Mesh) -> Vec<BoundaryRule> {
    vec![
        BoundaryRule {
            mode: BoundaryMode::Rebound,
            e_normal: 0.8,
            e_tangential: 0.9,
        };
        m.boundaries.len()
    ]
}

/// Gas that pushes towards the plane x = 0.5 from both sides, cell by cell
/// (the inflow towards the axis of a spray that entrains the gas).
struct Converging(Vec<[f64; 3]>);

impl Gas for Converging {
    fn at(&self, cell: usize, _x: [f64; 3]) -> GasState {
        let ux = if self.0[cell][0] < 0.5 { 0.5 } else { -0.5 };
        GasState::isothermal([ux, 0.0, 0.0], 1.2, 1.8e-5)
    }
}

#[test]
fn a_parcel_on_a_face_pushed_from_both_sides_is_not_lost() {
    // the nozzle sits on the face x = 0.5 between two cells; the gas of each
    // cell pushes the parcels back across it: two-way tracking must not end
    // every sub-step at that face after no time
    let mut s = BoxSpec::new(Kind::Hex, 3, [4, 4, 4]);
    s.lower = [0.0; 3];
    s.upper = [1.0; 3];
    let m = s.generate().unwrap();
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let mut cfg = SpraySection {
        injectors: vec![Injector {
            name: "a".into(),
            position: [0.5, 0.9, 0.4],
            direction: [0.0, -1.0, 0.0],
            speed: 2.0,
            half_angle: 0.0,
            mass_flow: 1e-3,
            diameter: 100e-6,
            size: None,
            parcels_per_second: 2000.0,
            start: 0.0,
            end: None,
            temperature: None,
        }],
        ..SpraySection::default()
    };
    cfg.coupling.mode = CouplingMode::TwoWay;
    let gas = Converging(geo.cell_centroid.clone());
    let mut sp = Spray::new(&cfg, &env, sp_rules(&m), [0.0, 0.0, 0.0]).unwrap();
    for n in 0..20 {
        sp.advance(&env, &gas, n as f64 * 0.01, 0.01).unwrap();
    }
    assert!(sp.parcels.len() > 100, "{}", sp.parcels.len());
    assert_eq!(sp.tally.lost_parcels, 0);
}
