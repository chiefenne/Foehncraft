//! Verification of the unilateral contact with rigid planes (MS5 part 4):
//! a reed on its seat and against its stop.
//!
//! - Static, augmented Lagrangian: a cantilever pressed onto a seat below
//!   its tip is the propped cantilever with a settled support, whose tip
//!   force is P − 3EI g/L³ (Timoshenko & Gere, *Mechanics of Materials*);
//!   a block on a seat equals the block on rollers.
//! - In time: a mass on a spring driven by a constant force against a
//!   stop. With a rigid stop and an elastic impact the motion is
//!   x = F/k (1 − cos ωt) up to the stop at t₁, then mirrored: back at
//!   x = 0, at rest, at 2t₁ (the energy is conserved). The penalty contact
//!   converges to it as the step and the penetration shrink.
//! - Case keys: a reed lifted by a pressure pulse hits its stop and falls
//!   back onto its seat.

#![allow(clippy::needless_range_loop)]

use foehncraft_structure::case::{self, Case};
use foehncraft_structure::codes;
use foehncraft_structure::contact::{self, Contacts, Plane, StepCache};
use foehncraft_structure::dynamics::{GenAlpha, Rayleigh};
use foehncraft_structure::model::{DOFS, Model};

fn rel(a: f64, b: f64) -> f64 {
    ((a - b) / b).abs()
}

const E: f64 = 206e9;
const L: f64 = 0.03;
const W: f64 = 6e-3;
const H: f64 = 0.2e-3;

fn reed(extra: &str) -> Case {
    Case::from_toml(&format!(
        r#"
[analysis]
static = true
[materials.steel]
young = {E}
poisson = 0.3
density = 7850.0
[sections.reed]
shape = "rectangle"
width = {W}
thickness = {H}
[[beams]]
name = "reed"
from = [0.0, 0.0, 0.0]
to = [{L}, 0.0, 0.0]
elements = 20
section = "reed"
material = "steel"
theory = "euler-bernoulli"
[[supports]]
at = [0.0, 0.0, 0.0]
fix = ["all"]
[points]
tip = [{L}, 0.0, 0.0]
{extra}
"#
    ))
    .unwrap()
}

#[test]
fn cantilever_pressed_onto_a_seat_is_the_propped_cantilever() {
    let ei = E * W * H.powi(3) / 12.0;
    let g = 1.0e-3; // the seat 1 mm below the tip
    for p in [0.05, 0.2, 1.0] {
        let c = reed(&format!(
            r#"
[[loads]]
at = [{L}, 0.0, 0.0]
force = [0.0, 0.0, {}]
[[contacts]]
name = "seat"
part = "reed"
point = [0.0, 0.0, {}]
normal = [0.0, 0.0, 1.0]
"#,
            -p, -g
        ));
        let out = case::run(&c).unwrap();
        let s = out.statics.as_ref().unwrap();
        let n = out.built.points["tip"];
        let w = s.u[DOFS * n + 2];
        let free = -p * L.powi(3) / (3.0 * ei);
        let j = out.to_json("c");
        let force = j["static"]["contact"]["planes"][0]["force"]
            .as_f64()
            .unwrap();
        if free > -g {
            // the tip does not reach the seat
            assert!(rel(w, free) < 1e-9, "{p}: {w} {free}");
            assert_eq!(force, 0.0);
        } else {
            let r = p - 3.0 * ei * g / L.powi(3);
            eprintln!(
                "P = {p} N: tip {w:.9e} (seat {:.9e}), contact force {force:.9e} (exact {r:.9e})",
                -g
            );
            assert!((w + g).abs() < 1e-9 * L, "{p}: {w}");
            assert!(rel(force, r) < 1e-6, "{p}: {force} {r}");
            assert_eq!(j["static"]["contact"]["planes"][0]["nodes_in_contact"], 1);
            // the clamp takes the rest of the load
            let rf = j["static"]["reaction_force"][2].as_f64().unwrap();
            assert!(rel(rf, p - r) < 1e-6, "{rf}");
        }
    }
    // pulled away from the seat: no contact, the free answer
    let c = reed(
        r#"
[[loads]]
at = [0.03, 0.0, 0.0]
force = [0.0, 0.0, 0.5]
[[contacts]]
name = "seat"
part = "reed"
point = [0.0, 0.0, 0.0]
normal = [0.0, 0.0, 1.0]
"#,
    );
    let out = case::run(&c).unwrap();
    let w = out.statics.as_ref().unwrap().u[DOFS * out.built.points["tip"] + 2];
    assert!(rel(w, 0.5 * L.powi(3) / (3.0 * ei)) < 1e-9);
}

#[test]
fn block_on_a_seat_is_the_block_held_where_it_touches() {
    // a hex20 block pressed onto a seat by a uniform pressure. On rollers
    // the corner nodes of the quad8 faces carry a pulling reaction (the
    // consistent load of a uniform pressure on a serendipity face is −1/12
    // at a corner, 1/3 at a mid-side node), so on a seat they lift off: the
    // unilateral answer is the block held by rollers at the nodes that
    // touch, with every one of their reactions pushing
    let block = |extra: &str| {
        Case::from_toml(&format!(
            r#"
[analysis]
static = true
[materials.steel]
young = 200e9
poisson = 0.3
density = 7800.0
[[solids]]
name = "b"
material = "steel"
block = {{ origin = [0.0, 0.0, 0.0], size = [0.02, 0.01, 0.005], divisions = [4, 2, 1] }}
[[supports]]
part = "b"
surface = "xmin"
fix = ["ux"]
[[supports]]
part = "b"
surface = "ymin"
fix = ["uy"]
[[pressures]]
part = "b"
surface = "zmax"
value = 1e6
[contact]
tolerance = 1e-13
{extra}
"#
        ))
        .unwrap()
    };
    let seat = case::run(&block(
        "[[contacts]]\nname = \"seat\"\npart = \"b\"\nsurface = \"zmin\"\npoint = [0.0, 0.0, 0.0]\nnormal = [0.0, 0.0, 1.0]\n",
    ))
    .unwrap();
    let m = &seat.built.model;
    let b = &seat.statics.as_ref().unwrap().u;
    let big = b.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let bottom = &seat.built.contacts[0].nodes;
    assert_eq!(bottom.len(), 15 + 22);
    let mut touching = Vec::new();
    for &n in bottom {
        let w = b[DOFS * n + 2];
        // within the tolerance, 10⁻¹³ of the model size
        assert!(w > -2.0e-15, "node {n} inside the seat: {w}");
        if w < 1e-6 * big {
            touching.push(n);
        }
    }
    // the mid-side nodes touch, the corners of the quad8 faces lift off
    assert_eq!(touching.len(), 22, "{touching:?}");
    let j = seat.to_json("s");
    let c = &j["static"]["contact"]["planes"][0];
    assert_eq!(c["nodes_in_contact"], 22);
    let f = c["force"].as_f64().unwrap();
    assert!(rel(f, 1e6 * 0.02 * 0.01) < 1e-9, "{f}");
    // the same block on rollers at the touching nodes
    let rollers: String = touching
        .iter()
        .map(|&n| {
            let x = m.nodes[n];
            format!(
                "[[supports]]\nat = [{:e}, {:e}, {:e}]\nfix = [\"uz\"]\n",
                x[0], x[1], x[2]
            )
        })
        .collect();
    let held = case::run(&block(&rollers)).unwrap();
    let hs = held.statics.as_ref().unwrap();
    // equal up to the penetration left by the tolerance (2·10⁻¹⁵ m)
    for i in 0..b.len() {
        assert!(
            (hs.u[i] - b[i]).abs() < 1e-6 * big,
            "{i}: {} {}",
            hs.u[i],
            b[i]
        );
    }
    for &n in &touching {
        // the rollers push the block up
        assert!(hs.reactions[DOFS * n + 2] > 0.0, "{n}");
    }
    eprintln!(
        "block on a seat: corners lift by up to {:.3e} m ({:.1} % of the largest displacement), {} solves",
        bottom
            .iter()
            .map(|&n| b[DOFS * n + 2])
            .fold(0.0f64, f64::max),
        100.0
            * bottom
                .iter()
                .map(|&n| b[DOFS * n + 2])
                .fold(0.0f64, f64::max)
            / big,
        j["static"]["contact"]["solves"]
    );
    // a node starting inside the plane is refused
    let e = case::run(&block(
        "[[contacts]]\nname = \"seat\"\npart = \"b\"\npoint = [0.0, 0.0, 0.001]\nnormal = [0.0, 0.0, 1.0]\n",
    ))
    .err()
    .unwrap();
    assert_eq!(e.code, codes::CASE);
    assert!(e.message.contains("inside"), "{e}");
}

/// A mass on a spring along z, pushed by a constant force towards a stop.
fn oscillator(k: f64, m: f64) -> Model {
    let mut model = Model::default();
    model.add_node([0.0, 0.0, 0.0]);
    model.springs.push((2, k));
    model.masses.push((0, m));
    model.forces[2] = 1.0;
    model
}

#[test]
fn impact_on_a_stop_converges_to_the_elastic_bounce() {
    let (k, m, f): (f64, f64, f64) = (400.0, 0.01, 1.0); // ω = 200 rad/s
    let w = (k / m).sqrt();
    let g = 0.75 * f / k; // the stop below the free swing 2F/k
    let t1 = (1.0 - g * k / f).acos() / w; // the impact
    let v1 = w * ((f / k).powi(2) - (g - f / k).powi(2)).sqrt();
    let model = oscillator(k, m);
    let sys = model.assemble();
    // the penalty model in continuous time: inside the stop the mass swings
    // about x_e on the springs k + k_c with the impact speed v₁, for τ; then
    // the free motion mirrors the first one: back at x = 0, at rest, at
    // 2t₁ + τ (the energy is kept). As k_c grows, τ and the penetration
    // vanish: the rigid elastic bounce.
    let exact = |kc: f64| {
        let xe = (f + kc * g) / (k + kc);
        let w2 = ((k + kc) / m).sqrt();
        let y0 = g - xe;
        let a = (y0 * y0 + (v1 / w2).powi(2)).sqrt();
        let tau = 2.0 * (v1 / w2).atan2(y0) / w2;
        (xe + a - g, tau)
    };
    let run = |penalty: f64, steps: usize| {
        let stop = Plane {
            name: "stop".into(),
            point: [0.0, 0.0, g],
            normal: [0.0, 0.0, -1.0],
            offset: 0.0,
            nodes: vec![0],
        };
        let c = Contacts::new(&model, &sys, vec![stop], penalty, 1e-12).unwrap();
        let (deep, tau) = exact(c.k);
        let dt = (2.0 * t1 + tau) / steps as f64;
        let ga = GenAlpha::new(&sys, Rayleigh::default(), dt, 1.0).unwrap();
        let mut st = ga.initial(0.0, vec![0.0], vec![0.0], &sys.f).unwrap();
        let mut cache = StepCache::default();
        let mut top = 0.0f64;
        let mut hit = false;
        for _ in 0..steps {
            let (next, cs) = ga
                .step_contact(&model, &c, &mut cache, &st, &sys.f, &sys.f)
                .unwrap();
            hit |= cs.force[0] > 0.0;
            top = top.max(next.u[0]);
            st = next;
        }
        assert!(hit);
        let err = (st.u[0].abs() + st.v[0].abs() / w) / g;
        eprintln!(
            "penalty {penalty:.0e}, {steps} steps ({:.1} in contact): |x| + |v|/ω at 2t1 + τ = {err:.3e} of the gap; penetration {:.4e} m (penalty model {deep:.4e} m); {} factorisations",
            tau / dt,
            top - g,
            cache.factorisations
        );
        (err, top - g, deep)
    };
    // the steps converge to the penalty model
    let mut errs = Vec::new();
    for steps in [200, 400, 800, 1600] {
        let (err, pen, deep) = run(1e3, steps);
        assert!(rel(pen, deep) < 0.05, "{pen} {deep}");
        errs.push(err);
    }
    for i in 1..errs.len() {
        assert!(errs[i] < errs[i - 1], "{errs:?}");
    }
    assert!(errs[3] < 0.01, "{errs:?}");
    // a stiffer penalty approaches the rigid stop: the penetration falls
    // as 1/√penalty
    let mut pens = Vec::new();
    for (penalty, steps) in [(1e2, 1600), (1e4, 1600), (1e6, 16000)] {
        let (_, pen, deep) = run(penalty, steps);
        assert!(rel(pen, deep) < 0.05, "{pen} {deep}");
        pens.push(pen / g);
    }
    assert!(pens[2] < 2e-3 && pens[0] / pens[2] > 50.0, "{pens:?}");
}

#[test]
fn reed_hits_its_stop_and_falls_back_on_its_seat() {
    // a pressure pulse lifts the reed against its stop 1 mm above the seat;
    // then it falls back onto the seat; the plate is 0.2 mm thick, so its
    // faces keep 0.1 mm from the planes
    let text = r#"
[analysis]
static = false
[materials.steel]
young = 206e9
poisson = 0.3
density = 7850.0
[[plates]]
name = "reed"
origin = [0.0, -3e-3, 0.0]
edge_a = [0.03, 0.0, 0.0]
edge_b = [0.0, 6e-3, 0.0]
divisions = [20, 4]
thickness = 0.2e-3
material = "steel"
[[supports]]
along = [[0.0, -3e-3, 0.0], [0.0, 3e-3, 0.0]]
fix = ["all"]
[[pressures]]
plate = "reed"
value = 2000.0
history = "pulse"
[histories.pulse]
table = [[0.0, 0.0], [0.0005, 1.0], [0.002, 1.0], [0.0025, 0.0]]
[[contacts]]
name = "seat"
part = "reed"
point = [0.0, 0.0, -0.1e-3]
normal = [0.0, 0.0, 1.0]
offset = 0.1e-3
[[contacts]]
name = "stop"
part = "reed"
point = [0.0, 0.0, 1.1e-3]
normal = [0.0, 0.0, -1.0]
offset = 0.1e-3
[transient]
end_time = 0.012
time_step = 1e-5
rho_infinity = 0.8
[damping]
ratios = [[150.0, 0.03], [1000.0, 0.03]]
[points]
tip = [0.03, 0.0, 0.0]
"#;
    let c = Case::from_toml(text).unwrap();
    let out = case::run(&c).unwrap();
    let j = out.to_json("r");
    let tr = &j["transient"];
    let w: Vec<f64> = tr["points"]["tip"]["displacement"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v[2].as_f64().unwrap())
        .collect();
    let pen = tr["contact"]["penetration"].as_f64().unwrap();
    let max = w.iter().cloned().fold(f64::MIN, f64::max);
    let min = w.iter().cloned().fold(f64::MAX, f64::min);
    let fs: Vec<f64> = tr["contact"]["planes"][1]["force"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    let fseat: Vec<f64> = tr["contact"]["planes"][0]["force"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    eprintln!(
        "tip from {min:.4e} to {max:.4e} m, penetration {pen:.3e} m, stop force up to {:.3} N, seat force up to {:.3} N, {} factorisations",
        fs.iter().cloned().fold(0.0, f64::max),
        fseat.iter().cloned().fold(0.0, f64::max),
        tr["contact"]["factorisations"]
    );
    assert!(max <= 1e-3 + pen + 1e-12 && max > 0.99e-3, "{max}");
    assert!(min >= -pen - 1e-12, "{min}");
    assert!(fs.iter().any(|&f| f > 0.0), "the reed reaches its stop");
    assert!(fseat.iter().any(|&f| f > 0.0), "the reed lands on its seat");
    assert!(pen < 2e-5, "{pen}");
    let csv = out.to_csv().unwrap();
    assert!(
        csv.lines()
            .next()
            .unwrap()
            .ends_with(",contact_seat,contact_stop")
    );
    // the modal method does not take contact
    let mut m = c.clone();
    m.transient.as_mut().unwrap().method = case::Method::Modal;
    m.transient.as_mut().unwrap().modes = 4;
    let e = case::run(&m).err().unwrap();
    assert!(e.message.contains("direct"), "{e}");
    let _ = contact::State::default();
}
