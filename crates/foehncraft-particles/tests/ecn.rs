//! The diagnostics of the ECN modelling standard (`[spray.ecn]`, D315):
//! liquid length from the projected liquid volume (PLV = ∫ LVF dy across
//! the spray) and vapour penetration, on parcels and cells built by hand so
//! that the exact answer is known.

use foehncraft_particles::spray::{line_of_sight, plv_lengths, plv_profile};
use foehncraft_particles::{Ecn, Injector, SpraySection, injector_stats};

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
