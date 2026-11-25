//! Two fluids by the volume-of-fluid method (L3, stage MP1): the advection
//! of α by THINC, CICSAM and upwind on the device kernels of the run
//! (Zalesak's slotted disk, one revolution; conservation to round-off),
//! the refusals of the first version.

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_vof_{name}_{}", std::process::id()))
}

/// Zalesak's slotted disk (Zalesak 1979): radius 15, slot 5 wide and 25
/// long, in [0, 100]² on `n` × `n` cells, turned about (50, 50) once in
/// T = 628.
fn zalesak(scheme: &str, n: usize, kind: &str) -> FlowCase {
    let text = format!(
        r#"
case_format = "0.7"
[mesh.generate]
type = "{kind}"
dimension = 2
upper = [100.0, 100.0, 1.0]
cells = [{n}, {n}, 1]
[fluid]
density = 1.0
viscosity = 1.0
[vof]
scheme = "{scheme}"
[[vof.phase]]
name = "a"
density = 1.0
viscosity = 1.0
[[vof.phase]]
name = "b"
density = 1.0
viscosity = 1.0
[[vof.region]]
shape = "cylinder"
centre = [50.0, 75.0, 0.5]
radius = 15.0
[[vof.region]]
shape = "box"
lower = [47.5, 60.0, -1.0]
upper = [52.5, 85.0, 2.0]
alpha = 0.0
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[time]
time_step = 1.0
steps = 1
"#
    );
    FlowCase::from_toml(&text, &cases()).unwrap()
}

/// The relative L1 error Σ|α − α₀| V / Σ α₀ V after one revolution, the
/// relative change of the phase volume and α's range.
fn revolve(scheme: &str, n: usize, kind: &str) -> (f64, f64, [f64; 2]) {
    let case = zalesak(scheme, n, kind);
    let w = std::f64::consts::PI / 314.0;
    let u = move |x: [f64; 3], _t: f64| [w * (50.0 - x[1]), w * (x[0] - 50.0), 0.0];
    // one revolution in 1256 steps of 0.5 (sub-steps keep Co ≤ 0.25)
    let r = foehncraft_flow::vof::advect(&case, &u, 0.5, 1256, None).unwrap();
    let v = &r.volume;
    let e: f64 = (0..v.len())
        .map(|c| (r.end[c] - r.start[c]).abs() * v[c])
        .sum();
    let v0: f64 = (0..v.len()).map(|c| r.start[c] * v[c]).sum();
    let t = &r.totals;
    // conservation: the change of the volume is the clipped part alone
    // (walls: nothing leaves)
    let dv = (t.volume - t.volume_start + t.clipped_volume) / t.volume_start;
    assert!(dv.abs() < 1e-12, "{scheme}: volume not conserved: {dv:e}");
    (e / v0, (t.volume - t.volume_start) / t.volume_start, t.alpha_range)
}

#[test]
fn zalesak_disk_thinc_cicsam_upwind() {
    let mut errs = Vec::new();
    for s in ["thinc", "cicsam", "upwind"] {
        let (e, dv, range) = revolve(s, 100, "hex");
        eprintln!("zalesak {s}: L1 {e:.4}, volume change {dv:.2e}, alpha {range:?}");
        assert!(range[0] >= 0.0 && range[1] <= 1.0);
        errs.push(e);
    }
    // the sharpening schemes beat upwind by far; THINC is the sharpest
    assert!(errs[0] < errs[1] && errs[1] < errs[2], "{errs:?}");
    assert!(errs[0] < 0.12, "THINC L1 error {}", errs[0]);
}

#[test]
fn zalesak_disk_thinc_on_polygons() {
    let (e, _, range) = revolve("thinc", 100, "poly");
    eprintln!("zalesak thinc poly: L1 {e:.4}, alpha {range:?}");
    assert!(e < 0.25, "THINC L1 error on polygons {e}");
}

/// The code a case text is refused with.
fn refused(text: &str) -> String {
    let e = match FlowCase::from_toml(text, &cases()) {
        Err(e) => e,
        Ok(c) => match run(&c, "refused", &out("refused"), |_| {}) {
            Err(e) => e,
            Ok(_) => panic!("the case was accepted"),
        },
    };
    diag::of_error(&e).code.to_string()
}

#[test]
fn refusals() {
    let base = r#"
case_format = "0.7"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.1]
cells = [4, 4, 1]
[fluid]
density = 1.0
viscosity = 1.0
[vof]
[[vof.phase]]
name = "a"
density = 1000.0
viscosity = 1e-3
[[vof.phase]]
name = "b"
density = 1.0
viscosity = 1e-5
[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
"#;
    // steady
    assert_eq!(refused(base), "VOF002");
    let tr = format!("{base}[time]\ntime_step = 0.01\nsteps = 1\n");
    // turbulence
    assert_eq!(
        refused(&format!("{tr}[turbulence]\nmodel = \"k-omega-sst\"\n")),
        "VOF002"
    );
    // one phase only
    let one = tr.replacen(
        "[[vof.phase]]\nname = \"b\"\ndensity = 1.0\nviscosity = 1e-5\n",
        "",
        1,
    );
    assert_eq!(refused(&one), "VOF001");
    // a bad Courant number
    assert_eq!(
        refused(&tr.replace("[vof]\n", "[vof]\nmax_courant = 2.0\n")),
        "VOF001"
    );
    // an inflow on a boundary the case does not have
    assert_eq!(
        refused(&tr.replace("[boundary.xmin]", "[vof.inflow]\ninlet = 1.0\n[boundary.xmin]")),
        "VOF001"
    );
}

/// The computer-graphics benchmark of the disk (owner rule: interface
/// topics benchmark graphics algorithms too): a semi-Lagrangian level set
/// (Staniforth & Côté 1991; Enright, Fedkiw et al. 2002 use it in graphics)
/// with Catmull-Rom interpolation and the exact back-trace of the rotation,
/// no reinitialisation. φ is the signed distance of the slotted disk made
/// by CSG; α of a cell is clamp(½ − φ/h); the error is measured as for VOF.
fn level_set_disk(n: usize) -> f64 {
    let h = 100.0 / n as f64;
    let x = |i: usize| (i as f64 + 0.5) * h;
    let sd = |px: f64, py: f64| {
        let disk = ((px - 50.0).powi(2) + (py - 75.0).powi(2)).sqrt() - 15.0;
        // the slot as a box: negative inside
        let dx = (px - 50.0).abs() - 2.5;
        let dy = (py - 72.5).abs() - 12.5;
        let slot = dx.max(dy).min(0.0) + (dx.max(0.0).powi(2) + dy.max(0.0).powi(2)).sqrt();
        disk.max(-slot)
    };
    let mut phi: Vec<f64> = (0..n * n).map(|c| sd(x(c % n), x(c / n))).collect();
    let alpha = |p: &[f64]| -> Vec<f64> { p.iter().map(|v| (0.5 - v / h).clamp(0.0, 1.0)).collect() };
    let a0 = alpha(&phi);
    let cr = |p: [f64; 4], t: f64| {
        let v = 0.5
            * (2.0 * p[1]
                + (p[2] - p[0]) * t
                + (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3]) * t * t
                + (3.0 * (p[1] - p[2]) + p[3] - p[0]) * t * t * t);
        // clamp to the neighbours (no new extrema)
        v.clamp(p[1].min(p[2]), p[1].max(p[2]))
    };
    let w = std::f64::consts::PI / 314.0;
    let (dt, steps) = (0.5, 1256);
    let (s, c) = ((w * dt).sin(), (w * dt).cos());
    for _ in 0..steps {
        let old = phi.clone();
        let at = |i: isize, j: isize| {
            let i = i.clamp(0, n as isize - 1) as usize;
            let j = j.clamp(0, n as isize - 1) as usize;
            old[i + n * j]
        };
        for (k, p) in phi.iter_mut().enumerate() {
            // back-trace: turn the point back by ω dt about (50, 50)
            let (rx, ry) = (x(k % n) - 50.0, x(k / n) - 50.0);
            let (bx, by) = (50.0 + c * rx + s * ry, 50.0 - s * rx + c * ry);
            let (fx, fy) = (bx / h - 0.5, by / h - 0.5);
            let (i0, j0) = (fx.floor() as isize, fy.floor() as isize);
            let (tx, ty) = (fx - i0 as f64, fy - j0 as f64);
            let row = |j: isize| cr([at(i0 - 1, j), at(i0, j), at(i0 + 1, j), at(i0 + 2, j)], tx);
            *p = cr([row(j0 - 1), row(j0), row(j0 + 1), row(j0 + 2)], ty);
        }
    }
    let a = alpha(&phi);
    let e: f64 = a.iter().zip(&a0).map(|(p, q)| (p - q).abs()).sum();
    e / a0.iter().sum::<f64>()
}

#[test]
fn zalesak_disk_level_set_benchmark() {
    let e = level_set_disk(100);
    eprintln!("zalesak semi-Lagrangian level set: L1 {e:.4}");
    // the graphics reference: kept as a number to compare with, not a bar
    // THINC has to pass
    assert!(e.is_finite() && e < 0.5, "level set L1 error {e}");
}
