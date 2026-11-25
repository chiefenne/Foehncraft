//! Acceptance of the spray datums and statistics (SP1b, D262 to D265):
//! the liquid through a plane closes against the boundaries and the parcels
//! in the flow (V28.13), sizes drawn from a distribution have its moments
//! (V28.14), the statistics of the sums equal a direct calculation
//! (V28.15), and the sums are bitwise equal on 1 to 4 threads and after a
//! restart.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_mesh::query::{Bvh, SurfaceFaces};
use foehncraft_particles::{
    Bound, BoundaryMode, BoundaryRule, Datum, Env, Gas, GasState, HistSpec, Injector, Shape,
    SizeDistribution, Spray, SpraySection, Surface, Volume,
};

struct Uniform([f64; 3]);

impl Gas for Uniform {
    fn at(&self, _: usize, _: [f64; 3]) -> GasState {
        GasState {
            u: self.0,
            rho: 1.2,
            mu: 1.8e-5,
        }
    }
}

fn rr() -> SizeDistribution {
    SizeDistribution::RosinRammler {
        mean: 80e-6,
        spread: 2.5,
        min: 10e-6,
        max: 300e-6,
    }
}

fn injector(rate: f64) -> Injector {
    Injector {
        name: "nozzle".into(),
        position: [0.05, 0.5, 0.5],
        direction: [1.0, 0.0, 0.0],
        speed: 2.0,
        half_angle: 20.0,
        mass_flow: 1e-3,
        diameter: 0.0,
        size: Some(rr()),
        parcels_per_second: rate,
        start: 0.0,
        end: None,
    }
}

fn square_at(x: f64) -> Bvh {
    let mut s = SurfaceFaces::default();
    let t = [
        [[x, -0.1, -0.1], [x, 1.1, -0.1], [x, 1.1, 1.1]],
        [[x, -0.1, -0.1], [x, 1.1, 1.1], [x, -0.1, 1.1]],
    ];
    for (k, tri) in t.into_iter().enumerate() {
        s.tris.push(tri);
        s.tri_face.push(k as u32);
        s.global_id.push(k as u64);
        s.mesh_face.push(usize::MAX);
        s.normal.push([1.0, 0.0, 0.0]);
    }
    Bvh::build(s)
}

struct Setup {
    mesh: foehncraft_mesh::Mesh,
    geo: Geometry,
}

fn setup(kind: Kind) -> Setup {
    let mut s = BoxSpec::new(kind, 3, [10, 4, 4]);
    s.lower = [0.0; 3];
    s.upper = [1.0; 3];
    let mesh = s.generate().unwrap();
    let geo = Geometry::compute(&mesh);
    Setup { mesh, geo }
}

/// The spray of the checks: lateral boundaries mirror, the ends let the
/// parcels go; datums: a plane at x = 0.5, the same plane as cell sides,
/// a triangulated square at x = 0.7, the outlet, and a box.
fn spray(st: &Setup, env: &Env) -> Spray {
    let m = &st.mesh;
    let rules: Vec<BoundaryRule> = m
        .boundaries
        .iter()
        .map(|b| BoundaryRule {
            mode: if b.name == "xmin" || b.name == "xmax" {
                BoundaryMode::Escape
            } else {
                BoundaryMode::Mirror
            },
            e_normal: 1.0,
            e_tangential: 1.0,
        })
        .collect();
    let cfg = SpraySection {
        injectors: vec![injector(2000.0)],
        ..SpraySection::default()
    };
    let mut sp = Spray::new(&cfg, env, rules, [0.0; 3]).unwrap();
    let hist = HistSpec::default();
    let side: Vec<i8> = (0..m.n_cells)
        .map(|c| {
            if st.geo.cell_centroid[c][0] < 0.5 {
                -1
            } else {
                1
            }
        })
        .collect();
    let xmax = m.boundaries.iter().position(|b| b.name == "xmax").unwrap();
    let n = [1.0, 0.0, 0.0];
    sp.set_datums(vec![
        Datum {
            name: "plane".into(),
            shape: Shape::Surface(Surface::Plane {
                point: [0.5, 0.5, 0.5],
                normal: n,
                bound: Bound::None,
            }),
            hist,
        },
        Datum {
            name: "faces".into(),
            shape: Shape::Surface(Surface::Cells {
                side,
                plane: Some(([0.5, 0.5, 0.5], n, Bound::None)),
                faces: None,
            }),
            hist,
        },
        Datum {
            name: "stl".into(),
            shape: Shape::Surface(Surface::Triangles {
                bvh: Box::new(square_at(0.7)),
            }),
            hist,
        },
        Datum {
            name: "outlet".into(),
            shape: Shape::Surface(Surface::Boundary { boundary: xmax }),
            hist,
        },
        Datum {
            name: "box".into(),
            shape: Shape::Volume(Volume::Box {
                lower: [0.0, 0.0, 0.0],
                upper: [0.3, 1.0, 1.0],
            }),
            hist,
        },
    ]);
    sp
}

const DT: f64 = 0.01;

fn run(sp: &mut Spray, env: &Env, from: usize, to: usize) {
    let gas = Uniform([1.5, 0.2, -0.1]);
    for s in from..to {
        sp.advance(env, &gas, s as f64 * DT, DT).unwrap();
    }
}

/// V28.13: the liquid through a plane equals what left through the outlet
/// plus what is in the flow beyond the plane (no evaporation); the plane
/// and the same plane on the cell sides agree within one parcel mass; the
/// outlet datum books what the mass book calls escaped.
#[test]
fn liquid_through_datums_closes() {
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly] {
        let st = setup(kind);
        let loc = Locator::new(&st.mesh, &st.geo);
        let env = Env {
            mesh: &st.mesh,
            geo: &st.geo,
            loc: &loc,
        };
        let mut sp = spray(&st, &env);
        run(&mut sp, &env, 0, 100);
        let xmax = st
            .mesh
            .boundaries
            .iter()
            .position(|b| b.name == "xmax")
            .unwrap();
        let escaped = sp.tally.escaped_mass[xmax];
        let beyond = |x0: f64| -> f64 {
            (0..sp.parcels.len())
                .filter(|&i| sp.parcels.x[i][0] >= x0)
                .map(|i| sp.parcel_mass(i))
                .sum()
        };
        let m_parcel = 1e-3 / 2000.0;
        let t = &sp.datum_tally;
        let plane = t[0].total.net();
        let expect = escaped + beyond(0.5);
        assert!(
            escaped > 0.0 && beyond(0.5) > 0.0,
            "{kind:?}: the check needs both"
        );
        assert!(
            (plane - expect).abs() <= 1e-12 * expect,
            "{kind:?}: plane {plane:e} against {expect:e}"
        );
        // the cell sides: what crossed is what left plus what is in the
        // cells of the far side; on the hexahedra the face layer is the
        // plane itself
        let far: f64 = (0..sp.parcels.len())
            .filter(|&i| st.geo.cell_centroid[sp.parcels.cell[i] as usize][0] >= 0.5)
            .map(|i| sp.parcel_mass(i))
            .sum();
        let faces = t[1].total.net();
        assert!(
            (faces - escaped - far).abs() <= 1e-12 * faces,
            "{kind:?}: faces {faces:e} against {:e}",
            escaped + far
        );
        if kind == Kind::Hex {
            assert!(
                (faces - plane).abs() <= m_parcel * 1.0001,
                "{kind:?}: faces {faces:e} against the plane {plane:e}"
            );
        }
        let stl = t[2].total.net();
        let expect = escaped + beyond(0.7);
        assert!(
            (stl - expect).abs() <= m_parcel * 1.0001,
            "{kind:?}: surface {stl:e} against {expect:e}"
        );
        assert_eq!(t[3].total.mass_pos, escaped, "{kind:?}: outlet datum");
        assert_eq!(t[3].total.mass_neg, 0.0);
        // the box holds the parcels near the nozzle
        let inside: f64 = (0..sp.parcels.len())
            .filter(|&i| sp.parcels.x[i][0] <= 0.3)
            .map(|i| sp.parcel_mass(i))
            .sum();
        let mb = sp.moments(Some(4), &st.mesh.cell_zone).mass;
        assert!(
            (mb - inside).abs() <= 1e-12 * inside,
            "{kind:?}: box {mb} against {inside}"
        );
    }
}

/// V28.14 and V28.15: the droplets in the flow, drawn from the
/// Rosin-Rammler distribution, have its Sauter diameter within the sampling
/// error; the statistics of the sums equal a direct calculation from the
/// parcels.
#[test]
fn statistics_follow_the_distribution_and_the_parcels() {
    let st = setup(Kind::Hex);
    let loc = Locator::new(&st.mesh, &st.geo);
    let env = Env {
        mesh: &st.mesh,
        geo: &st.geo,
        loc: &loc,
    };
    let cfg = SpraySection {
        injectors: vec![{
            let mut i = injector(200_000.0);
            i.speed = 0.0;
            i.position = [0.5, 0.5, 0.5];
            i
        }],
        ..SpraySection::default()
    };
    let rules = vec![
        BoundaryRule {
            mode: BoundaryMode::Stick,
            e_normal: 1.0,
            e_tangential: 1.0
        };
        st.mesh.boundaries.len()
    ];
    let mut sp = Spray::new(&cfg, &env, rules, [0.0; 3]).unwrap();
    // still gas, no gravity: 20 000 parcels stay where they were born
    sp.advance(&env, &Uniform([0.0; 3]), 0.0, 0.1).unwrap();
    let n = sp.parcels.len();
    assert_eq!(n, 20_000);
    let spec = HistSpec::default();
    let s = sp.moments(None, &st.mesh.cell_zone).stats(&spec);
    // direct sums over the parcels, in the order of their identities
    let (mut s2, mut s3, mut s0, mut s1, mut s4) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for i in 0..n {
        let (d, c) = (sp.parcels.d[i], sp.parcels.count[i]);
        s0 += c;
        s1 += c * d;
        s2 += c * d * d;
        s3 += c * d * d * d;
        s4 += c * d * d * d * d;
    }
    let close = |a: f64, b: f64| (a / b - 1.0).abs() < 1e-11;
    assert!(close(s.d10, s1 / s0) && close(s.d32, s3 / s2) && close(s.d43, s4 / s3));
    assert!(close(s.d20, (s2 / s0).sqrt()) && close(s.d30, (s3 / s0).cbrt()));
    // equal-mass parcels: D32 = 1 / E_mass[1/d]
    let dist = rr();
    let inv = dist.mass_moment(-1.0);
    let expect = 1.0 / inv;
    let sigma = ((dist.mass_moment(-2.0) - inv * inv).sqrt() / inv) / (n as f64).sqrt();
    assert!(
        (s.d32 / expect - 1.0).abs() < 4.0 * sigma,
        "D32 {} against {expect} (relative sampling error {sigma:e})",
        s.d32
    );
    // Dv50 is the mass median of the distribution within one bin and the
    // sampling error
    let med = dist.quantile(0.5);
    let bin = (1e-3f64 / 1e-6).powf(1.0 / 40.0);
    assert!(
        s.dv50 / med < bin * 1.01 && med / s.dv50 < bin * 1.01,
        "Dv50 {} against {med}",
        s.dv50
    );
}

/// The datum sums are the same bits on 1, 2 and 4 threads and after a
/// restart in the middle.
#[test]
fn datum_sums_are_bitwise_on_threads_and_restart() {
    let st = setup(Kind::Poly);
    let loc = Locator::new(&st.mesh, &st.geo);
    let env = Env {
        mesh: &st.mesh,
        geo: &st.geo,
        loc: &loc,
    };
    let sums = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let mut sp = spray(&st, &env);
            run(&mut sp, &env, 0, 40);
            serde_json::to_string(&sp.datum_tally).unwrap()
        })
    };
    let one = sums(1);
    assert!(one == sums(2), "2 threads differ");
    assert!(one == sums(4), "4 threads differ");
    let mut sp = spray(&st, &env);
    run(&mut sp, &env, 0, 20);
    let saved = sp.state().unwrap();
    let mut again = spray(&st, &env);
    again.restore(&env, &saved).unwrap();
    run(&mut again, &env, 20, 40);
    assert!(
        serde_json::to_string(&again.datum_tally).unwrap() == one,
        "the restart differs"
    );
}
