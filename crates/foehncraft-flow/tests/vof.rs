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
    // conservation: the change of the volume is the clipped part and what
    // crossed the boundary (the rotation has a normal velocity on the box
    // sides; α there is all but zero)
    let dv = (t.volume - t.volume_start + t.clipped_volume + t.volume_out) / t.volume_start;
    assert!(dv.abs() < 1e-12, "{scheme}: volume not conserved: {dv:e}");
    (
        e / v0,
        (t.volume - t.volume_start) / t.volume_start,
        t.alpha_range,
    )
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
    assert!(errs[0] < 0.06, "THINC L1 error {}", errs[0]);
}

#[test]
fn zalesak_disk_thinc_on_polygons() {
    let (e, _, range) = revolve("thinc", 100, "poly");
    eprintln!("zalesak thinc poly: L1 {e:.4}, alpha {range:?}");
    assert!(e < 0.11, "THINC L1 error on polygons {e}");
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
        refused(&tr.replace(
            "[boundary.xmin]",
            "[vof.inflow]\ninlet = 1.0\n[boundary.xmin]"
        )),
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
    let alpha =
        |p: &[f64]| -> Vec<f64> { p.iter().map(|v| (0.5 - v / h).clamp(0.0, 1.0)).collect() };
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

/// Martin & Moyce (1952), Phil. Trans. R. Soc. A 244, 312-324, part IV,
/// table 2 (n² = 2, a = 2¼ in): the surge front Z = z/a at T = n t √(g/a)
/// (docs/verification/DATA_SOURCES.md section 3h).
const MARTIN_MOYCE: [(f64, f64); 12] = [
    (1.19, 1.44),
    (1.56, 1.89),
    (1.87, 2.33),
    (2.26, 2.78),
    (2.57, 3.22),
    (2.88, 3.67),
    (3.25, 4.11),
    (3.59, 4.56),
    (3.89, 5.00),
    (4.20, 5.44),
    (4.54, 5.89),
    (4.85, 6.33),
];

/// The collapse of a water column in air (density ratio 828) on h = a/8:
/// the surge front against the measurements, the water volume conserved.
#[test]
fn dam_break_martin_moyce() {
    let text = std::fs::read_to_string(cases().join("vof_dam_break.toml")).unwrap();
    let text = text
        .replace("cells = [168, 45, 1]", "cells = [84, 22, 1]")
        .replace("time_step = 2.5e-4", "time_step = 5e-4")
        .replace("steps = 1200", "steps = 560")
        .replace("fields_every = 100", "fields_every = 0");
    let case = FlowCase::from_toml(&text, &cases()).unwrap();
    let dir = out("dam");
    run(&case, "dam", &dir, |_| {}).unwrap();
    let csv = std::fs::read_to_string(dir.join("vof.csv")).unwrap();
    let mut lines = csv.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let (ct, cx, cv, cc) = (col("time"), col("x_max"), col("volume"), col("clipped"));
    let rows: Vec<Vec<f64>> = lines
        .map(|l| {
            l.split(',')
                .map(|v| v.parse::<f64>().unwrap_or(f64::NAN))
                .collect()
        })
        .collect();
    let a: f64 = 0.05715;
    let h = 0.6 / 84.0;
    let tz: Vec<(f64, f64)> = rows
        .iter()
        .map(|r| (r[ct] * (2.0 * 9.81 / a).sqrt(), (r[cx] + 0.5 * h) / a))
        .collect();
    let mut worst: f64 = 0.0;
    for (t, z) in MARTIN_MOYCE {
        let k = tz.iter().position(|p| p.0 >= t).unwrap();
        let (p, q) = (tz[k - 1], tz[k]);
        let zs = p.1 + (q.1 - p.1) * (t - p.0) / (q.0 - p.0);
        let d = (zs - z) / z;
        eprintln!(
            "dam break: T {t:.2}: Z measured {z:.2}, computed {zs:.2} ({:+.1} %)",
            100.0 * d
        );
        worst = worst.max(d.abs());
    }
    // the front runs ahead of the measurements by 5 to 17 % on h = a/16
    // (the gate of the experiment took time to lift); a/8 alike
    assert!(
        worst < 0.2,
        "surge front off the measurements by {worst:.3}"
    );
    let (first, last) = (&rows[0], rows.last().unwrap());
    let dv = (last[cv] - first[cv] + last[cc]) / first[cv];
    assert!(dv.abs() < 1e-9, "water volume not conserved: {dv:e}");
}

/// The static drop: the Laplace jump σ/R within 15 %, the spurious
/// currents below 0.3 m/s (balanced CSF, density ratio 828).
#[cfg(feature = "hdf5")]
#[test]
fn static_drop_laplace_jump() {
    let case = FlowCase::from_toml(
        &std::fs::read_to_string(cases().join("vof_static_drop.toml")).unwrap(),
        &cases(),
    )
    .unwrap();
    let dir = out("drop");
    run(&case, "drop", &dir, |_| {}).unwrap();
    let f = hdf5_metno::File::open(dir.join("fields.vtkhdf")).unwrap();
    let g = f.group("VTKHDF/internal").unwrap();
    let read = |n: &str| {
        g.dataset(&format!("CellData/{n}"))
            .unwrap()
            .read_raw::<f64>()
            .unwrap()
    };
    let (a, p, s) = (read("volume_fraction"), read("pressure"), read("speed"));
    let mean = |pick: &dyn Fn(f64) -> bool| {
        let v: Vec<f64> = (0..a.len()).filter(|&c| pick(a[c])).map(|c| p[c]).collect();
        v.iter().sum::<f64>() / v.len() as f64
    };
    let dp = mean(&|x| x > 0.99) - mean(&|x| x < 0.01);
    let umax = s.iter().copied().fold(0.0, f64::max);
    eprintln!("static drop: dp {dp:.2} Pa (exact 35), spurious {umax:.3} m/s");
    assert!((dp - 35.0).abs() < 0.15 * 35.0, "Laplace jump {dp}");
    assert!(umax < 0.3, "spurious currents {umax}");
}

/// The columns of `vof.csv` by name, one row per step.
fn vof_rows(dir: &Path) -> (Vec<String>, Vec<Vec<f64>>) {
    let csv = std::fs::read_to_string(dir.join("vof.csv")).unwrap();
    let mut lines = csv.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<String> = lines
        .next()
        .unwrap()
        .split(',')
        .map(str::to_string)
        .collect();
    let rows = lines
        .map(|l| {
            l.split(',')
                .map(|v| v.parse::<f64>().unwrap_or(f64::NAN))
                .collect()
        })
        .collect();
    (head, rows)
}

/// Runs a phase-change example with `edits` and returns the time and the
/// liquid volume per step, and the last row by column name.
fn phase_run(
    file: &str,
    name: &str,
    edits: &[(&str, &str)],
) -> (Vec<(f64, f64)>, std::collections::BTreeMap<String, f64>) {
    let mut text = std::fs::read_to_string(cases().join(file)).unwrap();
    for (a, b) in edits {
        assert!(text.contains(a), "{file}: no `{a}`");
        text = text.replace(a, b);
    }
    let case = FlowCase::from_toml(&text, &cases()).unwrap();
    let dir = out(name);
    let rep = run(&case, name, &dir, |_| {}).unwrap();
    // the volume account of the run (net boundary volume flow against the
    // volume made) is one of its checks
    assert!(rep.failures.is_empty(), "{name}: {:?}", rep.failures);
    let (head, rows) = vof_rows(&dir);
    let col = |n: &str| head.iter().position(|h| h == n).unwrap();
    let (ct, cv) = (col("time"), col("volume"));
    let series = rows.iter().map(|r| (r[ct], r[cv])).collect();
    let last = rows.last().unwrap();
    let named = head.iter().cloned().zip(last.iter().copied()).collect();
    (series, named)
}

/// The worst relative deviation of a computed length from the exact one
/// over the rows after the first `skip`.
fn worst(
    series: &[(f64, f64)],
    skip: usize,
    len: impl Fn(f64) -> f64,
    exact: impl Fn(f64) -> f64,
) -> f64 {
    series
        .iter()
        .skip(skip)
        .map(|&(t, v)| (len(v) / exact(t) - 1.0).abs())
        .fold(0.0, f64::max)
}

/// The Stefan problem (MP2): the vapour layer on a wall 10 K above T_sat
/// grows as δ = 2λ√(α_v t) (water at 1 atm); the liquid volume closes
/// with the outflow and the evaporated mass.
#[test]
fn phase_change_stefan() {
    let (s, last) = phase_run(
        "vof_stefan.toml",
        "stefan",
        &[
            ("steps = 661", "steps = 300"),
            ("fields_every = 100", "fields_every = 0"),
        ],
    );
    let (l, a) = (1.0e-3, 4.0e-10);
    let (lam, av) = (0.067_741_985_626_831_74, 1.978_878_679_494_709e-5);
    let w = worst(&s, 20, |v| l - v / a, |t| 2.0 * lam * (av * t).sqrt());
    eprintln!("stefan: worst deviation of δ {:.2} %", 100.0 * w);
    assert!(w < 0.02, "vapour layer off the exact δ by {w:.4}");
    // V₀ = V + V_out + m/ρ_l, V_out = m (1/ρ_v − 1/ρ_l) (the continuity
    // takes the volume the phase change makes)
    let (rl, rv) = (958.372729, 0.597623116);
    let (v0, m) = (9.0e-4 * a, last["evaporated"]);
    let liquid = (last["volume"] + last["volume_out"] + m / rl - v0) / v0;
    let made = last["volume_out"] / (m * (1.0 / rv - 1.0 / rl)) - 1.0;
    eprintln!("stefan: liquid account {liquid:e}, volume made {made:e}");
    assert!(liquid.abs() < 1e-6, "liquid volume account {liquid:e}");
    assert!(
        made.abs() < 1e-3,
        "outflow against the volume made {made:e}"
    );
}

/// The sucking problem (MP2): the vapour layer grows into a liquid 2 K
/// superheated, δ = 2β√(α_l t); the liquid carries its thermal layer.
#[test]
fn phase_change_sucking() {
    let (s, _) = phase_run(
        "vof_sucking.toml",
        "sucking",
        &[
            ("steps = 800", "steps = 400"),
            ("fields_every = 200", "fields_every = 0"),
        ],
    );
    let (l, a) = (1.2e-3, 9.0e-12);
    let (beta, al) = (3.389_359_871_358_738_7, 1.675_804_574_892_131_6e-7);
    let w = worst(&s, 0, |v| l - v / a, |t| 2.0 * beta * (al * t).sqrt());
    eprintln!("sucking: worst deviation of δ {:.2} %", 100.0 * w);
    assert!(w < 0.025, "vapour layer off the exact δ by {w:.4}");
}

/// Scriven's bubble (MP2) on a thin sector of the sphere: R = 2β√(α_l t)
/// in water 1.25 K superheated, from 0.1 mm to 0.2 mm.
#[test]
fn phase_change_scriven() {
    let (s, _) = phase_run(
        "vof_scriven.toml",
        "scriven",
        &[
            ("steps = 1447", "steps = 540"),
            ("fields_every = 250", "fields_every = 0"),
        ],
    );
    // the sector (examples/meshes/sphere_sector.py): a square pyramid of
    // half angle 2°, solid angle 4 tan²θ/n³, n = √(1 + 2 tan²θ)
    let t2 = 2.0f64.to_radians().tan();
    let om = 4.0 * t2 * t2 / (1.0 + 2.0 * t2 * t2).powf(1.5);
    let (r0, r1) = (2.0e-6f64, 5.0e-4f64);
    let total = om * (r1.powi(3) - r0.powi(3)) / 3.0;
    let radius = |v: f64| (3.0 * (total - v + om * r0.powi(3) / 3.0) / om).cbrt();
    let (beta, al) = (4.066_930_541_557_649_5, 1.675_804_574_892_131_6e-7);
    let w = worst(&s, 0, radius, |t| 2.0 * beta * (al * t).sqrt());
    eprintln!("scriven: worst deviation of R {:.2} %", 100.0 * w);
    assert!(w < 0.04, "bubble radius off the exact R by {w:.4}");
}
