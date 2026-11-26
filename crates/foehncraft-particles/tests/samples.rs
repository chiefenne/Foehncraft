//! Injection from measured droplets (D322): the sample table, the draw by
//! mass, the start of the parcels on their radius within the sector and the
//! keys that the table replaces.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, Env, Gas, GasState, Injector, Samples, Spray, SpraySection,
};
use std::sync::Arc;

const TABLE: &str = "# two measured droplets and one of no weight
r, d, u_axial, u_radial, u_tangential, weight
0.002, 20e-6, 30.0, 2.0, 0.5, 3.0
0.004, 40e-6, 10.0, -1.0, 0.0, 1.0
0.006, 10e-6, 50.0, 0.0, 0.0, 0.0
";

struct Still;

impl Gas for Still {
    fn at(&self, _cell: usize, _x: [f64; 3]) -> GasState {
        GasState::isothermal([0.0; 3], 1.2, 1.8e-5)
    }
}

#[test]
fn a_sample_table_reads_and_draws_by_mass() {
    let s = Samples::parse(TABLE).unwrap();
    assert_eq!(s.rows.len(), 3);
    assert_eq!(s.cumulative, vec![0.75, 1.0, 1.0]);
    assert_eq!(s.rows[0].u, [30.0, 2.0, 0.5]);
    assert_eq!(s.pick(0.0).d, 20e-6);
    assert_eq!(s.pick(0.7499).d, 20e-6);
    assert_eq!(s.pick(0.75).d, 40e-6);
    assert_eq!(s.pick(0.999_999).d, 40e-6);
    // without weights every row has the same share; white space separates
    let s = Samples::parse("r d u_axial u_radial\n0 1e-5 1 0\n0.001 2e-5 2 0\n").unwrap();
    assert_eq!(s.cumulative, vec![0.5, 1.0]);
    for (bad, why) in [
        ("r, d, u_axial\n0, 1e-5, 1\n", "u_radial"),
        (
            "r, d, u_axial, u_radial, v\n0, 1e-5, 1, 0, 0\n",
            "unknown column",
        ),
        ("r, d, u_axial, u_radial\n-1, 1e-5, 1, 0\n", "r must be"),
        ("r, d, u_axial, u_radial\n0, 0, 1, 0\n", "d must be"),
        ("r, d, u_axial, u_radial\n0, 1e-5, 1\n", "3 values"),
        (
            "r, d, u_axial, u_radial, weight\n0, 1e-5, 1, 0, 0\n",
            "no droplet",
        ),
        ("", "empty"),
    ] {
        let e = format!("{:#}", Samples::parse(bad).unwrap_err());
        assert!(e.contains(why), "{bad:?}: {e}");
    }
}

fn sampled_injector() -> Injector {
    Injector {
        name: "plane".into(),
        position: [0.0, 0.0, 0.01],
        direction: [0.0, 0.0, 1.0],
        mass_flow: 1e-3,
        parcels_per_second: 1e9,
        samples: Some("table.csv".into()),
        sector: 90.0,
        flow_share: 0.25,
        sample_data: Some(Arc::new(Samples::parse(TABLE).unwrap())),
        ..Default::default()
    }
}

#[test]
fn sampled_parcels_start_on_their_radius_in_the_sector() {
    let mut b = BoxSpec::new(Kind::Hex, 3, [4, 4, 4]);
    b.lower = [0.0; 3];
    b.upper = [0.01, 0.01, 0.05];
    let m = b.generate().unwrap();
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let cfg = SpraySection {
        injectors: vec![sampled_injector()],
        ..SpraySection::default()
    };
    cfg.validate().unwrap();
    let rules = vec![
        BoundaryRule {
            mode: BoundaryMode::Escape,
            e_normal: 1.0,
            e_tangential: 1.0,
        };
        m.boundaries.len()
    ];
    let mut sp = Spray::new(&cfg, &env, rules, [0.0; 3]).unwrap();
    // 4000 parcels in 4 µs: they move a few tenths of a millimetre at most
    sp.advance(&env, &Still, 0.0, 4e-6).unwrap();
    let p = &sp.parcels;
    assert_eq!(p.len(), 4000);
    let (mut small, mut mass) = (0usize, 0.0);
    for i in 0..p.len() {
        let [x, y, z] = p.x[i];
        let v = p.v[i];
        let (r, phi) = (x.hypot(y), y.atan2(x));
        assert!(x >= 0.0 && y >= 0.0, "outside the sector: {x} {y}");
        let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        // the radial and tangential parts in the frame of the azimuth
        let (ur, ut) = (
            v[0] * phi.cos() + v[1] * phi.sin(),
            -v[0] * phi.sin() + v[1] * phi.cos(),
        );
        if p.d[i] == 20e-6 {
            small += 1;
            assert!((r - 0.002).abs() < 2e-5, "r = {r}");
            assert!((speed / (30.0f64.powi(2) + 4.0 + 0.25).sqrt() - 1.0).abs() < 0.02);
            assert!((ur / 2.0 - 1.0).abs() < 0.05 && (ut / 0.5 - 1.0).abs() < 0.05);
        } else {
            assert_eq!(p.d[i], 40e-6, "the row of no weight is drawn");
            assert!((r - 0.004).abs() < 2e-5, "r = {r}");
            assert!(v[2] > 9.0 && ur < -0.9);
        }
        assert!(z > 0.01);
        mass += cfg.liquid.density * std::f64::consts::PI / 6.0 * p.d[i].powi(3) * p.count[i];
    }
    // three quarters of the mass in the first row, by number of parcels
    let share = small as f64 / p.len() as f64;
    assert!((share - 0.75).abs() < 0.03, "share {share}");
    // the sector carries its share of the mass flow
    assert!(
        (mass / (0.25 * 1e-3 * 4e-6) - 1.0).abs() < 1e-9,
        "mass {mass}"
    );
}

#[test]
fn the_sample_table_replaces_the_droplet_keys() {
    let check = |f: &dyn Fn(&mut Injector), why: &str| {
        let mut i = sampled_injector();
        f(&mut i);
        let cfg = SpraySection {
            injectors: vec![i],
            ..SpraySection::default()
        };
        let e = format!("{:#}", cfg.validate().unwrap_err());
        assert!(e.contains(why), "{e}");
    };
    check(&|i| i.speed = 10.0, "leave out speed");
    check(&|i| i.half_angle = 5.0, "leave out half_angle");
    check(&|i| i.diameter = 1e-5, "leave out diameter");
    check(&|i| i.flow_share = 1.0, "flow_share = 0.25");
    check(&|i| i.sector = 0.0, "sector must be");
    check(
        &|i| {
            i.samples = None;
            i.sample_data = None;
            i.diameter = 1e-5;
            i.speed = 1.0;
        },
        "sector needs samples",
    );
}
