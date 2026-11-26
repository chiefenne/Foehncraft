//! The diagnostics of the ECN modelling standard (`[spray.ecn]`, D315):
//! liquid length from the projected liquid volume (PLV = ∫ LVF dy across
//! the spray) and vapour penetration, on parcels and cells built by hand so
//! that the exact answer is known.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::spray::{
    format_mirrors, line_of_sight, mirror_image_count, parse_mirrors, plv_lengths,
    plv_lengths_mirrored, plv_profile, plv_profile_mirrored,
};
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, Ecn, Env, Injector, MirrorPlane, Spray, SpraySection,
    injector_stats,
};

const RHO: f64 = 700.0;

fn injector(direction: [f64; 3]) -> Injector {
    Injector {
        name: "a".into(),
        position: [0.01, 0.02, 0.03],
        direction,
        speed: 500.0,
        half_angle: 0.0,
        mass_flow: 1e-3,
        diameter: 1e-5,
        size: None,
        parcels_per_second: 1e6,
        start: 0.0,
        end: None,
        temperature: None,
        ..Injector::default()
    }
}

/// Parcels on a line along the axis of `inj` from the nozzle to `len`,
/// `per_bin` of them in every bin of size `bin`, at the bin's quarter
/// points across the axis (all inside the bin across), each with the
/// liquid volume v(s) = v0 (1 − s/len): a PLV that falls linearly.
fn line(inj: &Injector, bin: f64, len: f64, v0: f64) -> (Vec<u64>, Vec<[f64; 3]>, Vec<f64>) {
    let n = (len / bin).round() as usize;
    let a = inj.direction;
    let (mut id, mut x, mut m) = (Vec::new(), Vec::new(), Vec::new());
    let mut k = 0u64;
    for i in 0..n {
        for j in 0..4 {
            let s = (i as f64 + 0.125 + 0.25 * j as f64) * bin;
            id.push(k);
            k += 1;
            // a small offset across, inside the bin [0, bin) across
            x.push([
                inj.position[0] + s * a[0],
                inj.position[1] + s * a[1],
                inj.position[2] + s * a[2],
            ]);
            m.push(RHO * v0 * (1.0 - s / len));
        }
    }
    (id, x, m)
}

#[test]
fn the_liquid_length_of_a_linear_plv_profile_is_exact() {
    let inj = injector([1.0, 0.0, 0.0]);
    let e = Ecn {
        plv: vec![2e-7, 2e-6],
        bin: 2.5e-4,
        ..Ecn::default()
    };
    let len = 0.02;
    // the PLV of a bin: four parcels of v0 (1 − s/len) over bin²; at the
    // bin centre s_c the mean of the four is v0 (1 − s_c/len) exactly
    let v0 = 4e-6 * e.bin * e.bin / 4.0;
    let (id, x, m) = line(&inj, e.bin, len, v0);
    let own: Vec<usize> = (0..id.len()).collect();
    let prof = plv_profile(&e, &inj, RHO, &own, &id, &x, &m);
    assert_eq!(prof.len(), 80);
    for (k, p) in prof.iter().enumerate() {
        let sc = (k as f64 + 0.5) * e.bin;
        let exact = 4e-6 * (1.0 - sc / len);
        assert!(
            (p / exact - 1.0).abs() < 1e-12,
            "bin {k}: {p} against {exact}"
        );
    }
    // PLV(s) = 4e-6 (1 − s/len) falls to 2e-6 at s = len/2, to 2e-7 at
    // 0.95 len: linear between the bin centres, so exact
    let l = plv_lengths(&e, &inj, RHO, &own, &id, &x, &m);
    assert!((l[0] - 0.95 * len).abs() < 1e-12, "{l:?}");
    assert!((l[1] - 0.5 * len).abs() < 1e-12, "{l:?}");
}

#[test]
fn the_volume_along_the_line_of_sight_is_summed_and_across_it_is_not() {
    let inj = injector([0.0, 0.0, 1.0]);
    let e = Ecn {
        plv: vec![1e-6],
        bin: 1e-3,
        line_of_sight: Some([1.0, 0.0, 0.0]),
        ..Ecn::default()
    };
    // at s = 2.5 mm two parcels 3 mm apart along the line of sight (x):
    // one bin; at s = 4.5 mm two parcels 3 mm apart across it (y): two bins
    let p = inj.position;
    let x = vec![
        [p[0] - 1.5e-3, p[1] + 1e-4, p[2] + 2.5e-3],
        [p[0] + 1.5e-3, p[1] + 1e-4, p[2] + 2.5e-3],
        [p[0], p[1] + 1e-4, p[2] + 4.5e-3],
        [p[0], p[1] + 3.1e-3, p[2] + 4.5e-3],
    ];
    let v = 0.6e-6 * e.bin * e.bin;
    let m = vec![RHO * v; 4];
    let id = vec![0, 1, 2, 3];
    let own = vec![0, 1, 2, 3];
    let prof = plv_profile(&e, &inj, RHO, &own, &id, &x, &m);
    assert!((prof[2] - 1.2e-6).abs() < 1e-18, "{prof:?}");
    assert!((prof[4] - 0.6e-6).abs() < 1e-18, "{prof:?}");
    // so the 1e-6 threshold is reached at 2.5 mm, not at 4.5 mm
    let l = plv_lengths(&e, &inj, RHO, &own, &id, &x, &m);
    assert!(l[0] > 2.5e-3 && l[0] < 3.5e-3, "{l:?}");
}

#[test]
fn the_default_line_of_sight_is_across_the_axis() {
    let e = Ecn::default();
    for a in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
        let l = line_of_sight(&e, a);
        let d = l[0] * a[0] + l[1] * a[1] + l[2] * a[2];
        assert!(d.abs() < 1e-15 && (l[0] * l[0] + l[1] * l[1] + l[2] * l[2] - 1.0).abs() < 1e-15);
    }
    // a slanted axis: the given direction loses its part along the axis
    let s = 0.5f64.sqrt();
    let a = [s, s, 0.0];
    let e = Ecn {
        line_of_sight: Some([1.0, 0.0, 0.0]),
        ..Ecn::default()
    };
    let l = line_of_sight(&e, a);
    assert!((l[0] - s).abs() < 1e-15 && (l[1] + s).abs() < 1e-15 && l[2].abs() < 1e-15);
}

#[test]
fn the_liquid_length_does_not_depend_on_the_parcel_order() {
    let inj = injector([0.0, 1.0, 0.0]);
    let mut cfg = SpraySection {
        injectors: vec![inj.clone()],
        ecn: Some(Ecn::default()),
        ..SpraySection::default()
    };
    cfg.liquid.density = RHO;
    let (id, x, m) = line(&inj, 2.5e-4, 0.01, 3e-6 * 2.5e-4 * 2.5e-4 / 4.0);
    let a = injector_stats(&cfg, &id, &x, &m, None);
    let mut ix: Vec<usize> = (0..id.len()).collect();
    ix.reverse();
    let id2: Vec<u64> = ix.iter().map(|&i| id[i]).collect();
    let x2: Vec<[f64; 3]> = ix.iter().map(|&i| x[i]).collect();
    let m2: Vec<f64> = ix.iter().map(|&i| m[i]).collect();
    let b = injector_stats(&cfg, &id2, &x2, &m2, None);
    assert_eq!(a[0].liquid_length.len(), 2);
    for (p, q) in a[0].liquid_length.iter().zip(&b[0].liquid_length) {
        assert_eq!(p.to_bits(), q.to_bits());
    }
    // 3e-6 falls to 2e-6 at a third of the line, to 2e-7 at 14/15 of it
    assert!((a[0].liquid_length[1] - 0.01 / 3.0).abs() < 1e-12);
    assert!((a[0].liquid_length[0] - 0.01 * 14.0 / 15.0).abs() < 1e-12);
    // without [spray.ecn] there is none
    cfg.ecn = None;
    assert!(
        injector_stats(&cfg, &id, &x, &m, None)[0]
            .liquid_length
            .is_empty()
    );
}

#[test]
fn the_vapour_penetration_is_the_farthest_cell_at_the_threshold() {
    let mut b = BoxSpec::new(Kind::Hex, 3, [1, 1, 1]);
    // around the nozzle at (0.01, 0.02, 0.03)
    b.lower = [0.0; 3];
    b.upper = [0.05; 3];
    let m = b.generate().unwrap();
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    // a slanted axis: the distance is the part along it
    let s = 0.5f64.sqrt();
    let inj = injector([s, s, 0.0]);
    let cfg = SpraySection {
        injectors: vec![inj.clone()],
        ecn: Some(Ecn::default()),
        ..SpraySection::default()
    };
    let rules = vec![
        BoundaryRule {
            mode: BoundaryMode::Rebound,
            e_normal: 0.5,
            e_tangential: 0.5,
        };
        m.boundaries.len()
    ];
    let mut sp = Spray::new(&cfg, &env, rules, [0.0; 3]).unwrap();
    let p = inj.position;
    let along = |d: f64, off: f64| [p[0] + d * s - off * s, p[1] + d * s + off * s, p[2]];
    let centre = vec![
        along(1e-3, 0.0),
        along(2e-3, 0.0),
        along(3e-3, 0.0),    // below the threshold
        along(2.5e-3, 4e-3), // off the axis, at the threshold
        along(-6e-3, 0.0),   // behind the nozzle
    ];
    let y = vec![0.01, 0.002, 0.000_9, 1e-3, 0.05];
    sp.set_vapour_penetration(&y, &centre);
    let v = sp.injector_stats()[0].vapour_penetration.unwrap();
    assert!((v - 2.5e-3).abs() < 1e-15, "{v}");
    // no cell at the threshold: zero
    sp.set_vapour_penetration(&[0.0; 5], &centre);
    assert_eq!(sp.injector_stats()[0].vapour_penetration, Some(0.0));
}

#[test]
fn the_ecn_table_is_checked() {
    let inj = injector([1.0, 0.0, 0.0]);
    let base = SpraySection {
        injectors: vec![inj],
        ecn: Some(Ecn::default()),
        ..SpraySection::default()
    };
    assert!(base.validate().is_ok());
    let bad = |e: Ecn| {
        let mut c = base.clone();
        c.ecn = Some(e);
        c.validate().unwrap_err().to_string()
    };
    let e = bad(Ecn {
        plv: vec![2e-7, -1.0],
        ..Ecn::default()
    });
    assert!(e.contains("plv"), "{e}");
    let e = bad(Ecn {
        line_of_sight: Some([2.0, 0.0, 0.0]),
        ..Ecn::default()
    });
    assert!(e.contains("line_of_sight") && e.contains("`a`"), "{e}");
    let e = bad(Ecn {
        vapour_fraction: 1.0,
        ..Ecn::default()
    });
    assert!(e.contains("vapour_fraction"), "{e}");
    let mut c = base.clone();
    c.steady = Some(Default::default());
    let e = c.validate().unwrap_err().to_string();
    assert!(e.contains("transient"), "{e}");
}

fn plane(normal: [f64; 3], offset: f64) -> MirrorPlane {
    MirrorPlane { normal, offset }
}

#[test]
fn the_mirror_planes_close_into_the_images_of_the_spray() {
    let y0 = plane([0.0, 1.0, 0.0], 0.0);
    let z0 = plane([0.0, 0.0, 1.0], 0.0);
    assert_eq!(mirror_image_count(&[]), Some(1));
    assert_eq!(mirror_image_count(&[y0]), Some(2));
    // the same plane twice (a symmetry plane split into two boundaries)
    assert_eq!(mirror_image_count(&[y0, y0]), Some(2));
    assert_eq!(mirror_image_count(&[y0, z0]), Some(4));
    // a sector of 60° between two planes through the x axis: 6 images
    let (s, c) = (60f64.to_radians().sin(), 60f64.to_radians().cos());
    assert_eq!(mirror_image_count(&[y0, plane([0.0, -c, s], 0.0)]), Some(6));
    // 50° = 5 · 180°/18: 36 images that cover the circle five times (the
    // spray refuses them, see below)
    let (s, c) = (50f64.to_radians().sin(), 50f64.to_radians().cos());
    assert_eq!(
        mirror_image_count(&[y0, plane([0.0, -c, s], 0.0)]),
        Some(36)
    );
    // one radian is no rational part of 180°: no finite set of images
    let (s, c) = (1f64.sin(), 1f64.cos());
    assert_eq!(mirror_image_count(&[y0, plane([0.0, -c, s], 0.0)]), None);
    // the text of the parcel file reads back to the bit
    let m = vec![vec![y0, plane([0.0, -0.6, 0.8], 1.25e-3)], vec![], vec![z0]];
    assert_eq!(parse_mirrors(&format_mirrors(&m)).unwrap(), m);
    assert!(parse_mirrors("0 1 0").is_err());
}

/// Parcels of a round spray along +x from `p`, inside a cone of 10°, with
/// masses that fall along the axis, in the quarter y > 0, z > 0 of the
/// planes y = 0 and z = 0 (the nozzle lies 1 µm off both).
fn quarter(p: [f64; 3]) -> (Vec<u64>, Vec<[f64; 3]>, Vec<f64>) {
    let (mut id, mut x, mut m) = (Vec::new(), Vec::new(), Vec::new());
    let mut k = 0u64;
    for i in 0..400 {
        let s = 0.01 * (f64::from(i) + 0.5) / 400.0;
        for j in 0..6 {
            let r = 0.17 * s * (f64::from(j) + 0.5) / 6.0;
            let th = std::f64::consts::FRAC_PI_2 * (f64::from((i * 7 + j * 3) % 11) + 0.5) / 11.0;
            let pos = [p[0] + s, r * th.cos(), r * th.sin()];
            if pos[1] > 0.0 && pos[2] > 0.0 {
                id.push(k);
                x.push(pos);
                m.push(RHO * 2e-13 * (1.0 - s / 0.012));
                k += 1;
            }
        }
    }
    (id, x, m)
}

#[test]
fn a_quarter_spray_with_its_mirror_images_is_the_whole_spray() {
    // the Spray A layout: the axis on two symmetry planes, seen along y
    let mut inj = injector([1.0, 0.0, 0.0]);
    inj.position = [0.0002, 1e-6, 1e-6];
    let e = Ecn {
        plv: vec![2e-8, 2e-7],
        ..Ecn::default()
    };
    let (id, x, m) = quarter(inj.position);
    let own: Vec<usize> = (0..id.len()).collect();
    let planes = [plane([0.0, 1.0, 0.0], 0.0), plane([0.0, 0.0, 1.0], 0.0)];
    // the whole spray built by hand: the quarter and its three images
    let (mut id4, mut x4, mut m4) = (Vec::new(), Vec::new(), Vec::new());
    for (q, (sy, sz)) in [(1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]
        .into_iter()
        .enumerate()
    {
        for i in 0..id.len() {
            id4.push(id[i] + 1_000_000 * q as u64);
            x4.push([x[i][0], sy * x[i][1], sz * x[i][2]]);
            m4.push(m[i]);
        }
    }
    let own4: Vec<usize> = (0..id4.len()).collect();
    let whole = plv_profile(&e, &inj, RHO, &own4, &id4, &x4, &m4);
    let mirrored = plv_profile_mirrored(&e, &inj, RHO, &own, &id, &x, &m, &planes);
    assert_eq!(whole.len(), mirrored.len());
    for (a, b) in whole.iter().zip(&mirrored) {
        assert!((a - b).abs() <= 1e-12 * a.abs(), "{a} {b}");
    }
    let lw = plv_lengths(&e, &inj, RHO, &own4, &id4, &x4, &m4);
    let lm = plv_lengths_mirrored(&e, &inj, RHO, &own, &id, &x, &m, &planes);
    for (a, b) in lw.iter().zip(&lm) {
        assert!((a - b).abs() <= 1e-12 * a, "{a} {b}");
    }
    // seen along y, the quarter alone holds half the liquid of each line
    // of sight: its projection is half (from 2 mm on, where every parcel
    // lies more than the 1 µm of the nozzle off the planes and so in the
    // bin of its image), and its liquid lengths are shorter
    let plain = plv_profile(&e, &inj, RHO, &own, &id, &x, &m);
    for (a, b) in whole.iter().zip(&plain).skip(8) {
        assert!((b / a - 0.5).abs() < 1e-12, "{a} {b}");
    }
    let lq = plv_lengths(&e, &inj, RHO, &own, &id, &x, &m);
    assert!(lq[0] < lw[0] && lq[1] < lw[1], "{lq:?} {lw:?}");
}

#[test]
fn the_spray_finds_the_mirror_planes_through_the_axis() {
    let mut b = BoxSpec::new(Kind::Hex, 3, [4, 2, 2]);
    b.lower = [0.0; 3];
    b.upper = [0.01, 0.002, 0.002];
    let m = b.generate().unwrap();
    let geo = Geometry::compute(&m);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let mut inj = injector([1.0, 0.0, 0.0]);
    inj.position = [0.0002, 1e-6, 1e-6];
    inj.hole_diameter = Some(9e-5);
    let rules = |mirror: &[&str]| -> Vec<BoundaryRule> {
        m.boundaries
            .iter()
            .map(|b| BoundaryRule {
                mode: if mirror.contains(&b.name.as_str()) {
                    BoundaryMode::Mirror
                } else {
                    BoundaryMode::Stick
                },
                e_normal: 0.5,
                e_tangential: 0.5,
            })
            .collect()
    };
    let names = |inj: &Injector, mirror: &[&str]| -> Vec<String> {
        let cfg = SpraySection {
            injectors: vec![inj.clone()],
            ecn: Some(Ecn::default()),
            ..SpraySection::default()
        };
        let sp = Spray::new(&cfg, &env, rules(mirror), [0.0; 3]).unwrap();
        sp.ecn_mirrors()[0].iter().map(|x| x.0.clone()).collect()
    };
    // ymax mirrors too, but the axis is 2 mm away from it
    assert_eq!(names(&inj, &["ymin", "zmin", "ymax"]), ["ymin", "zmin"]);
    let sp = Spray::new(
        &SpraySection {
            injectors: vec![inj.clone()],
            ..SpraySection::default()
        },
        &env,
        rules(&["ymin", "zmin"]),
        [0.0; 3],
    )
    .unwrap();
    let planes = &sp.ecn_mirrors()[0];
    assert_eq!(planes[0].1.normal, [0.0, -1.0, 0.0]);
    assert_eq!(planes[0].1.offset, 0.0);
    // a wall is no mirror
    assert_eq!(names(&inj, &["zmin"]), ["zmin"]);
    // half a millimetre off the plane is beyond half the hole
    inj.position[1] = 5e-4;
    assert_eq!(names(&inj, &["ymin", "zmin"]), ["zmin"]);
    // an axis across a plane (+y) does not lie in it
    inj.position = [0.005, 1e-6, 1e-6];
    inj.direction = [0.0, 1.0, 0.0];
    assert_eq!(names(&inj, &["ymin", "zmin"]), ["zmin"]);
}

#[test]
fn a_sector_must_tile_the_whole_spray() {
    // a box sheared into a wedge between the x-y plane (zmin) and a plane
    // through the x axis at `deg` to it (ymin), the nozzle on the x axis
    let wedge = |deg: f64| {
        let mut b = BoxSpec::new(Kind::Hex, 3, [4, 2, 2]);
        b.lower = [0.0; 3];
        b.upper = [0.01, 0.002, 0.002];
        let mut m = b.generate().unwrap();
        let (s, c) = (deg.to_radians().sin(), deg.to_radians().cos());
        for p in &mut m.points {
            *p = [p[0], p[1] + p[2] * c, p[2] * s];
        }
        m
    };
    let mut inj = injector([1.0, 0.0, 0.0]);
    // inside the wedge, 1 µm off both planes for every angle tried here
    inj.position = [0.0002, 4e-6, 1e-6];
    inj.hole_diameter = Some(9e-5);
    let cfg = SpraySection {
        injectors: vec![inj.clone()],
        ecn: Some(Ecn::default()),
        ..SpraySection::default()
    };
    let run = |deg: f64| {
        let m = wedge(deg);
        let geo = Geometry::compute(&m);
        let loc = Locator::new(&m, &geo);
        let env = Env {
            mesh: &m,
            geo: &geo,
            loc: &loc,
        };
        let rules: Vec<BoundaryRule> = m
            .boundaries
            .iter()
            .map(|b| BoundaryRule {
                mode: if ["ymin", "zmin"].contains(&b.name.as_str()) {
                    BoundaryMode::Mirror
                } else {
                    BoundaryMode::Stick
                },
                e_normal: 0.5,
                e_tangential: 0.5,
            })
            .collect();
        Spray::new(&cfg, &env, rules, [0.0; 3]).map(|sp| {
            let planes: Vec<MirrorPlane> = sp.ecn_mirrors()[0].iter().map(|x| x.1).collect();
            mirror_image_count(&planes)
        })
    };
    // 90° and 60°: the images tile the spray (4 and 6 of them)
    assert_eq!(run(90.0).unwrap(), Some(4));
    assert_eq!(run(60.0).unwrap(), Some(6));
    // 50° makes 36 images that cover the spray five times; one radian
    // makes no finite set: both refused
    for deg in [50.0, 1f64.to_degrees()] {
        let e = run(deg).unwrap_err().to_string();
        assert!(e.contains("do not tile the whole spray"), "{deg}: {e}");
    }
}
