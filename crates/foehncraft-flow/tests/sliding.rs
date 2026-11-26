//! The sliding mesh (D56, D182, stage 6 parts 6 and 7 of 6g).
//!
//! Variant A, conformal stepping: the rotor of a ring mesh turns by whole
//! sectors per time step. The checks: the shift of the rotor cells (a
//! bijection, N of them are the identity), the refusals of the case file,
//! laminar Taylor–Couette flow from rest against the exact solution on a
//! series of meshes (second order, Δt ∝ h), a fluid in solid-body rotation
//! that the scheme keeps, a rotor with a shift of two sectors, and a passive
//! scalar.
//!
//! Variant B, imprinting: the interface of a rotor and a stator with
//! different numbers of sectors is imprinted at the rotor angle of every
//! (free) time step. The same flows, with the refusals of its case file, and
//! the comparison with variant A on a mesh whose two pitches are equal.
//!
//! Turbulence on both variants: the spin-up of turbulent Taylor–Couette
//! flow (k-ω SST) on the sliding meshes against the same flow on a fixed
//! mesh, and a patch of turbulence (k-ε-ζ-f) that the fluid carries with
//! the rotor through half a turn against the same patch at rest.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use foehncraft_flow::frame::{self, MrfSection};
use foehncraft_flow::sliding::{SlidingSection, plan};
use foehncraft_fv::Discretisation;
use foehncraft_fv::problem::{Boundary, Problem};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::annulus::AnnulusSpec;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_sliding_{name}_{}", std::process::id()))
}

fn try_run(case: &FlowCase, name: &str) -> anyhow::Result<FlowReport> {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {});
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    try_run(case, name).unwrap()
}

fn text() -> String {
    std::fs::read_to_string(cases().join("taylor_couette_sliding.toml"))
        .unwrap()
        .replace("\r\n", "\n")
}

fn parse(t: &str) -> FlowCase {
    FlowCase::from_toml(t, &cases()).unwrap()
}

fn code(e: &anyhow::Error) -> String {
    diag::of_error(e).code.to_string()
}

fn l2(r: &FlowReport) -> f64 {
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0];
    assert!(s.converged);
    s.errors.as_ref().expect("an analytic solution").l2
}

/// The case at refinement `level` (16 · 2^level sectors) with the time step
/// of one sector per step, to the time `end`.
fn at_level(base: &FlowCase, level: u32, shift: usize, end: f64) -> FlowCase {
    let mut c = base.clone();
    c.verification.max_error = None;
    c.mesh.annulus.as_mut().unwrap().refine = level;
    c.sliding[0].shift = shift;
    let sectors = 16usize << level;
    let dt = c.sliding[0].time_step(sectors).unwrap();
    let mut t = c.time.unwrap();
    t.time_step = dt;
    t.steps = Some((end / dt).round() as usize);
    t.end_time = None;
    c.time = Some(t);
    c
}

fn discretisation(mesh: &foehncraft_mesh::Mesh, geo: &Geometry) -> Discretisation {
    let mut prob = Problem::default();
    for p in &mesh.boundaries {
        prob.boundary
            .insert(p.name.clone(), Boundary::zero_gradient());
    }
    Discretisation::new(mesh, geo, &prob).unwrap()
}

/// The rotor of a ring of 24 sectors: found from the geometry, shifted by
/// whole sectors (a bijection that returns after N steps).
#[test]
fn the_rotor_cells_are_shifted_by_whole_sectors() {
    let mesh = AnnulusSpec {
        n_theta: 24,
        n_r: 6,
        rotor_layers: 3,
        radial_ratio: 2.0,
        ..Default::default()
    }
    .generate()
    .unwrap();
    let geo = Geometry::compute(&mesh);
    let d = discretisation(&mesh, &geo);
    let entry = SlidingSection {
        zone: Some("rotor".into()),
        omega: Some(10.0),
        ..Default::default()
    };
    let cell_zone = frame::assign_cells(&[entry.frame()], &mesh).unwrap();
    let rotor: Vec<bool> = cell_zone.iter().map(Option::is_some).collect();
    let p = plan(&entry, &rotor, &d).unwrap();
    assert_eq!((p.sectors, p.rotor_cells), (24, 24 * 3));
    assert!((p.angle - std::f64::consts::TAU / 24.0).abs() < 1e-14);
    // a permutation of the cells that fixes the stator and every ring
    let mut seen = vec![false; p.cell.len()];
    for (i, &c) in p.cell.iter().enumerate() {
        assert!(!seen[c], "cell {c} twice");
        seen[c] = true;
        assert_eq!(rotor[i], rotor[c], "the shift keeps the rotor a rotor");
        if !rotor[i] {
            assert_eq!(i, c);
        } else {
            let (a, b) = (d.cell_centroid.get(i), d.cell_centroid.get(c));
            let (ra, rb) = (a[0].hypot(a[1]), b[0].hypot(b[1]));
            assert!((ra - rb).abs() < 1e-12 && (a[2] - b[2]).abs() < 1e-12);
        }
    }
    // N steps are a whole turn: the identity
    let mut at: Vec<usize> = (0..p.cell.len()).collect();
    for _ in 0..24 {
        at = at.iter().map(|&i| p.cell[i]).collect();
    }
    assert!(at.iter().enumerate().all(|(i, &c)| i == c));
    // two sectors per step: the shift twice
    let two = SlidingSection {
        shift: 2,
        ..entry.clone()
    };
    let p2 = plan(&two, &rotor, &d).unwrap();
    assert!(
        p2.cell
            .iter()
            .enumerate()
            .all(|(i, &c)| c == p.cell[p.cell[i]])
    );
    assert!((p2.angle - 2.0 * p.angle).abs() < 1e-14);
    // a whole turn per step is no step
    let whole = SlidingSection {
        shift: 24,
        ..entry.clone()
    };
    assert_eq!(code(&plan(&whole, &rotor, &d).unwrap_err()), "SLD002");
    // an axis that does not go through the centre: no cell has an image
    let off = SlidingSection {
        origin: [0.013, 0.0, 0.0],
        ..entry.clone()
    };
    let e = plan(&off, &rotor, &d).unwrap_err();
    assert_eq!(code(&e), "SLD004", "{e:#}");
    // the angular velocity sets the sense only: the shift is the same set of
    // pairs turned the other way
    let back = SlidingSection {
        omega: Some(-10.0),
        ..entry
    };
    let pb = plan(&back, &rotor, &d).unwrap();
    assert!(
        pb.cell.iter().enumerate().all(|(i, &c)| p.cell[c] == i),
        "the opposite turn is the inverse shift"
    );
}

/// The checks of the tables and of the time stepping.
#[test]
fn the_case_is_checked() {
    let valid = |t: &str| -> Result<(), String> {
        FlowCase::from_toml(t, &cases())
            .and_then(|c| c.validate())
            .map_err(|e| code(&e))
    };
    let base = text();
    valid(&base).unwrap();
    // a steady case
    let no_time = {
        let i = base.find("\n[time]\n").unwrap();
        let j = base.find("\n[[sliding]]\n").unwrap();
        format!("{}{}", &base[..i], &base[j..])
    };
    assert_eq!(valid(&no_time).unwrap_err(), "SLD001");
    // the entry
    let two = base.clone() + "\n[[sliding]]\nomega = 1.0\n";
    assert_eq!(valid(&two).unwrap_err(), "SLD002");
    let both = base.replace(
        "omega = 10.0\n\n[boundary",
        "omega = 10.0\nrpm = 5.0\n\n[boundary",
    );
    assert_eq!(valid(&both).unwrap_err(), "SLD002");
    let none = base.replace("omega = 10.0\n\n[boundary", "\n[boundary");
    assert_eq!(valid(&none).unwrap_err(), "SLD002");
    let stopped = base.replace("omega = 10.0\n\n[boundary", "omega = 0.0\n\n[boundary");
    assert_eq!(valid(&stopped).unwrap_err(), "SLD002");
    let zero = base.replace(
        "omega = 10.0\n\n[boundary",
        "omega = 10.0\nshift = 0\n\n[boundary",
    );
    assert_eq!(valid(&zero).unwrap_err(), "SLD002");
    // the stepping
    let adaptive = base.replace("[time]\n", "[time]\nadaptive = true\n");
    assert_eq!(valid(&adaptive).unwrap_err(), "SLD003");
    let sweep = base.replace("[verification]\n", "[verification]\nlevels = [0, 1]\n");
    assert_eq!(valid(&sweep).unwrap_err(), "SLD003");
}

/// What the mesh decides: the time step, the end time, the zone and the
/// boundaries.
#[test]
fn the_mesh_and_the_step_are_checked() {
    let base = text();
    let refused = |t: String| -> (String, String) {
        let c = parse(&t);
        let e = try_run(&c, "refused").unwrap_err();
        (code(&e), format!("{e:#}"))
    };
    // the time step that the sectors and the speed give is named
    let (c, m) = refused(base.replace("time_step = 0.019634954084936207", "time_step = 0.02"));
    assert_eq!(c, "SLD003", "{m}");
    assert!(m.contains("1.9634954085e-2"), "{m}");
    // an end time that is not a whole number of steps
    let (c, _) = refused(base.replace("steps = 510", "end_time = 10.0"));
    assert_eq!(c, "SLD003");
    // the zone
    let (c, m) = refused(base.replace("zone = \"rotor\"", "zone = \"rotr\""));
    assert_eq!(c, "SLD002", "{m}");
    assert!(m.contains("did you mean `rotor`"), "{m}");
    // a pressure boundary on the rotor
    let t = base.replace(
        "[boundary.inner]\ntype = \"wall\"\nvelocity = [\"-spin * y\", \"spin * x\", 0.0]\n",
        "[boundary.inner]\ntype = \"static-pressure-outlet\"\npressure = 0.0\n",
    );
    let (c, m) = refused(t);
    assert_eq!(c, "SLD004", "{m}");
    assert!(m.contains("boundary `inner`"), "{m}");
}

/// Taylor–Couette flow from rest: the error against the exact solution at
/// the end of the spin-up falls with the second power of the sector angle
/// and the time step (128, 512 and 2 048 cells).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn taylor_couette_converges_with_second_order() {
    let base = parse(&text());
    let e: Vec<f64> = (0..3)
        .map(|l| l2(&run_case(&at_level(&base, l, 1, 10.0), &format!("tc{l}"))))
        .collect();
    let order: Vec<f64> = e.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
    eprintln!("sliding Taylor–Couette: L2 {e:?}, orders {order:?}");
    assert!(e[1] < 0.02 && e[2] < 0.004, "{e:?}");
    assert!(
        order.iter().all(|o| (1.9..2.7).contains(o)),
        "orders {order:?}"
    );
}

/// Two sectors per step: the same flow with a time step twice as long and
/// four times the error of the angle (second order in the step).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn two_sectors_per_step() {
    let base = parse(&text());
    let one = l2(&run_case(&at_level(&base, 1, 1, 10.0), "shift1"));
    let two = l2(&run_case(&at_level(&base, 1, 2, 10.0), "shift2"));
    eprintln!("L2: one sector per step {one:e}, two {two:e}");
    assert!(two < 0.08 && two > one, "{one} {two}");
}

/// A fluid in solid-body rotation with the whole mesh turning, started from
/// the exact solution: the scheme keeps it (the pressure gradient turns the
/// vectors of the cells that move), with an error that does not grow.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn solid_body_rotation_is_kept() {
    let t = text()
        .replace("zone = \"rotor\"\n", "")
        .replace("rotor_layers = 4\n", "")
        .replace(
            "[boundary.outer]\ntype = \"wall\"\n",
            "[boundary.outer]\ntype = \"wall\"\nvelocity = [\"-spin * y\", \"spin * x\", 0.0]\n",
        )
        .replace(
            "inner_omega = 10.0",
            "inner_omega = 10.0\nouter_omega = 10.0",
        )
        .replace("[flow]\n", "[initial]\nfrom_analytic = true\n\n[flow]\n")
        .replace("max_error = 0.02", "max_error = 0.1");
    let base = parse(&t);
    // 0.6 s and 3 s: the error stays where it settled (the first two steps
    // of a run are backward differences of first order)
    let short = l2(&run_case(&at_level(&base, 1, 1, 0.6), "solid_short"));
    let long = l2(&run_case(&at_level(&base, 1, 1, 3.0), "solid_long"));
    eprintln!("solid body, L2: {short:e} at 0.6 s, {long:e} at 3 s");
    assert!(short < 0.03 && long < 0.03, "{short} {long}");
    assert!(long < 1.5 * short, "{short} {long}");
}

/// A passive scalar between the cylinders (θ = 1 inside, 0 outside, pure
/// conduction) with the rotor turning through it: the old level of a moved
/// cell is what that cell held, and the scalar ends in the exact solution.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn a_passive_scalar_ends_in_the_exact_solution() {
    let t = text()
        + "\n[scalars.theta]\ndiffusivity = 1e-3\ninitial = 0.0\n\
           exact = \"log(0.1 / sqrt(x^2 + y^2)) / log(2.0)\"\n\
           [boundary.inner.scalar.theta]\nvalue = 1.0\n\
           [boundary.outer.scalar.theta]\nvalue = 0.0\n";
    let base = parse(&t);
    let r = run_case(&at_level(&base, 1, 1, 12.0), "theta");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let e = r.solves[0].scalars[0].l2.expect("exact solution of theta");
    eprintln!("theta L2 {e:e}");
    assert!(e < 1.5e-2, "{e}");
}

/// The rotating zone as a moving reference frame of the same case is
/// steady only: the table of a sliding mesh needs a time, the table of the
/// frame must not have one (D55, D182).
#[test]
fn a_frame_and_a_sliding_mesh_do_not_mix() {
    let t = text() + "\n[[mrf]]\nzone = \"rotor\"\nomega = 10.0\n";
    let e = FlowCase::from_toml(&t, &cases())
        .and_then(|c| c.validate())
        .unwrap_err();
    assert_eq!(code(&e), "MRF001", "{e:#}");
    let _ = MrfSection::default();
}

// ---------------------------------------------------------------------
// variant B: the interface is imprinted in every step (D182)
// ---------------------------------------------------------------------

fn imprint_text() -> String {
    std::fs::read_to_string(cases().join("taylor_couette_imprint.toml"))
        .unwrap()
        .replace("\r\n", "\n")
}

/// The imprint case at refinement `level` (the sectors of the rotor and the
/// stator double with it) with a time step proportional to the cell size
/// (0.04 s on the coarsest mesh), to the time `end`.
fn imprint_at(base: &FlowCase, level: u32, end: f64) -> FlowCase {
    let mut c = base.clone();
    c.verification.max_error = None;
    c.mesh.annulus.as_mut().unwrap().refine = level;
    let dt = 0.04 / f64::from(1u32 << level);
    let mut t = c.time.unwrap();
    t.time_step = dt;
    t.steps = Some((end / dt).round() as usize);
    t.end_time = None;
    c.time = Some(t);
    c
}

/// The refusals of the case file that need no mesh.
#[test]
fn the_imprint_case_is_checked() {
    let valid = |t: &str| -> Result<(), String> {
        FlowCase::from_toml(t, &cases())
            .and_then(|c| c.validate())
            .map_err(|e| code(&e))
    };
    let base = imprint_text();
    valid(&base).unwrap();
    // the two boundaries go together
    let one = base.replace("stator_interface = \"stator_interface\"\n", "");
    assert_eq!(valid(&one).unwrap_err(), "SLD002");
    // the rotor is a zone, the time step is free (shift belongs to variant A)
    let no_zone = base.replace("zone = \"rotor\"\n", "");
    assert_eq!(valid(&no_zone).unwrap_err(), "SLD002");
    let shift = base.replace(
        "omega = 10.0\nrotor_interface",
        "omega = 10.0\nshift = 2\nrotor_interface",
    );
    assert_eq!(valid(&shift).unwrap_err(), "SLD002");
    // fully implicit in the step
    let extrapolated = base.replace("extrapolation = false", "extrapolation = true");
    assert_eq!(valid(&extrapolated).unwrap_err(), "SLD003");
    let one_outer = base.replace("outer_iterations = 4", "outer_iterations = 1");
    assert_eq!(valid(&one_outer).unwrap_err(), "SLD003");
    let adaptive = base.replace("[time]\n", "[time]\nadaptive = true\n");
    assert_eq!(valid(&adaptive).unwrap_err(), "SLD003");
    // a restart is not available yet
    let restart = base.replace("[output]\n", "[run]\nrestart = \"x.h5\"\n\n[output]\n");
    assert_eq!(valid(&restart).unwrap_err(), "SLD005");
}

/// What the mesh decides: the names of the interface, the zone and the
/// shape of the interface.
#[test]
fn the_imprint_mesh_is_checked() {
    let base = imprint_text();
    let refused = |t: String| -> (String, String) {
        let c = parse(&t);
        let e = try_run(&c, "imprint_refused").unwrap_err();
        (code(&e), format!("{e:#}"))
    };
    // an interface that is not on the mesh
    let (c, m) = refused(base.replace(
        "rotor_interface = \"rotor_interface\"",
        "rotor_interface = \"rotor_interfac\"",
    ));
    assert!(
        m.contains("SLD005") && m.contains("did you mean `rotor_interface`"),
        "{c}: {m}"
    );
    // a conformal mesh has no interface boundaries
    let (_, m) = refused(base.replace("n_theta_stator = 24\n", ""));
    assert!(m.contains("SLD005") && m.contains("has no boundary"), "{m}");
    // the zone, with the nearest name
    let (c, m) = refused(base.replace("zone = \"rotor\"", "zone = \"rotr\""));
    assert_eq!(c, "SLD002", "{m}");
    assert!(m.contains("did you mean `rotor`"), "{m}");
    // an interface boundary gets no condition
    let (c, m) = refused(format!(
        "{base}\n[boundary.rotor_interface]\ntype = \"wall\"\n"
    ));
    assert_eq!(c, "SLD005", "{m}");
    assert!(m.contains("interface boundary"), "{m}");
}

/// Taylor–Couette flow from rest with 2n sectors in the rotor and 3n in the
/// stator, a time step that no sector decides: the error against the exact
/// solution falls with the second power of the mesh width (160, 640 and
/// 5 120 cells, Δt ∝ h).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn taylor_couette_with_unequal_pitch_converges_with_second_order() {
    let base = parse(&imprint_text());
    let e: Vec<f64> = (0..3)
        .map(|l| l2(&run_case(&imprint_at(&base, l, 10.0), &format!("itc{l}"))))
        .collect();
    let order: Vec<f64> = e.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
    eprintln!("imprinted Taylor–Couette: L2 {e:?}, orders {order:?}");
    assert!(e[1] < 0.02 && e[2] < 0.004, "{e:?}");
    assert!(
        order.iter().all(|o| (1.9..2.8).contains(o)),
        "orders {order:?}"
    );
    // the mass that crosses the interface is the mass that leaves one side
    // and enters the other: the balance closes to round-off
    let r = run_case(&imprint_at(&base, 1, 2.0), "itc_mass");
    assert!(
        r.solves[0].mass_imbalance.abs() < 1e-12,
        "{}",
        r.solves[0].mass_imbalance
    );
}

/// With the same number of sectors on both sides and the time step of
/// variant A, imprinting and conformal stepping solve the same problem on
/// the same mesh: the errors agree (the imprint pieces are the faces of the
/// mesh at every step that turns the rotor by whole sectors).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn equal_pitch_agrees_with_conformal_stepping() {
    let conformal = l2(&run_case(
        &at_level(&parse(&text()), 1, 1, 10.0),
        "pitch_conformal",
    ));
    let mut c = parse(&imprint_text());
    c.verification.max_error = None;
    let n = {
        let a = c.mesh.annulus.as_mut().unwrap();
        a.refine = 1;
        a.n_theta_stator = a.n_theta;
        a.n_theta << a.refine
    };
    let dt = 2.0 * std::f64::consts::PI / (n as f64 * 10.0);
    let mut t = c.time.unwrap();
    t.time_step = dt;
    t.steps = Some((10.0 / dt).round() as usize);
    c.time = Some(t);
    let imprinted = l2(&run_case(&c, "pitch_imprinted"));
    eprintln!("equal pitch: conformal {conformal:e}, imprinted {imprinted:e}");
    assert!(
        (imprinted / conformal - 1.0).abs() < 0.1,
        "{conformal} {imprinted}"
    );
}

/// A fluid in solid-body rotation with both cylinders turning, started
/// from the exact solution: the rotor mesh turns through the fixed stator
/// mesh and the fluid, which does not move relative to the rotor, flows
/// through the stator cells; the error does not grow.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn solid_body_rotation_is_kept_by_imprinting() {
    let t = imprint_text()
        .replace(
            "[boundary.outer]\ntype = \"wall\"\n",
            "[boundary.outer]\ntype = \"wall\"\nvelocity = [\"-spin * y\", \"spin * x\", 0.0]\n",
        )
        .replace(
            "inner_omega = 10.0",
            "inner_omega = 10.0\nouter_omega = 10.0",
        )
        .replace("[flow]\n", "[initial]\nfrom_analytic = true\n\n[flow]\n");
    let base = parse(&t);
    let short = l2(&run_case(&imprint_at(&base, 1, 0.6), "isolid_short"));
    let long = l2(&run_case(&imprint_at(&base, 1, 3.0), "isolid_long"));
    eprintln!("imprinted solid body, L2: {short:e} at 0.6 s, {long:e} at 3 s");
    assert!(short < 0.03 && long < 0.03, "{short} {long}");
    assert!(long < 1.5 * short, "{short} {long}");
}

/// A passive scalar between the cylinders with the rotor turning through it:
/// the scalar of a cell is carried through the rebuilt solver, and it ends
/// in the exact solution.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn a_passive_scalar_ends_in_the_exact_solution_with_imprinting() {
    let t = imprint_text()
        + "\n[scalars.theta]\ndiffusivity = 1e-3\ninitial = 0.0\n\
           exact = \"log(0.1 / sqrt(x^2 + y^2)) / log(2.0)\"\n\
           [boundary.inner.scalar.theta]\nvalue = 1.0\n\
           [boundary.outer.scalar.theta]\nvalue = 0.0\n";
    let base = parse(&t);
    let r = run_case(&imprint_at(&base, 1, 12.0), "itheta");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let e = r.solves[0].scalars[0].l2.expect("exact solution of theta");
    eprintln!("imprinted theta L2 {e:e}");
    assert!(e < 1.5e-2, "{e}");
}

/// The outputs of a run on imprinted meshes: a probe in the rotor finds the
/// cell that is there now (not the cell that was there at the start), a
/// probe in the stator the stator's, and the fields of a mesh whose
/// interface cells are polyhedra are written.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn probes_and_fields_follow_the_rotor() {
    let t = imprint_text()
        + "\n[[output.probes]]\nname = \"rot\"\npoint = [0.06, 0.0, 0.005]\n\
           \n[[output.probes]]\nname = \"sta\"\npoint = [0.0, 0.09, 0.005]\n\
           \n[[output.lines]]\nname = \"radial\"\nfrom = [0.052, 0.0, 0.005]\n\
           to = [0.098, 0.0, 0.005]\npoints = 24\n";
    let t = t.replace("vtkhdf = false", "vtkhdf = true");
    let mut c = parse(&t);
    c.verification.max_error = None;
    let dir = out("outputs");
    let rep = run(&c, "outputs", &dir, |_| {}).unwrap();
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let csv = std::fs::read_to_string(dir.join("probes.csv")).unwrap();
    let last: Vec<f64> = csv
        .lines()
        .rfind(|l| !l.starts_with('#') && !l.is_empty())
        .unwrap()
        .split(',')
        .map(|x| x.parse().unwrap())
        .collect();
    let fields = dir.join("fields.vtkhdf").exists();
    let line = std::fs::read_to_string(dir.join("line_radial.csv")).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    // u_θ = A r + B / r of Taylor–Couette flow (r₁ = 0.05, r₂ = 0.10, Ω = 10)
    let exact = |r: f64| -10.0 * 0.0025 / 0.0075 * r + 10.0 * 0.0025 * 0.01 / 0.0075 / r;
    // columns: step, time, rot u v w p u_rel v_rel w_rel, sta u v w p
    // u_rel v_rel w_rel (D183: the velocity relative to the turning zone)
    let (rot_v, sta_u) = (last[3], last[9]);
    assert_eq!(last.len(), 16, "{csv}");
    assert!(
        (last[7] - (rot_v - 10.0 * 0.06)).abs() < 1e-12,
        "rotor v_rel"
    );
    assert_eq!((last[6], last[13], last[14]), (last[2], sta_u, last[10]));
    // the line: c_u − w_u is the blade speed Ω r in the rotor, 0 outside
    let mut head = None;
    for l in line.lines().filter(|l| !l.starts_with('#')) {
        if head.is_none() {
            head = Some(l.to_string());
            assert_eq!(l, "s,x,y,z,u,v,w,p,u_rel,v_rel,w_rel,c_r,c_u,c_a,w_u");
            continue;
        }
        let v: Vec<f64> = l.split(',').map(|x| x.parse().unwrap()).collect();
        let r = v[1].hypot(v[2]);
        let blade = if r < 0.075 { 10.0 * r } else { 0.0 };
        assert!((v[12] - v[14] - blade).abs() < 1e-12, "r = {r}: {l}");
        assert!(
            (v[12] - v[5]).abs() < 1e-12 && (v[11] - v[4]).abs() < 1e-12,
            "{l}"
        );
        assert!((v[12] / exact(r) - 1.0).abs() < 0.05, "c_u at r = {r}: {l}");
    }
    eprintln!(
        "probes: rotor v {rot_v} (exact {}), stator u {sta_u} (exact {})",
        exact(0.06),
        -exact(0.09)
    );
    assert!((rot_v / exact(0.06) - 1.0).abs() < 0.05, "{rot_v}");
    assert!((sta_u / -exact(0.09) - 1.0).abs() < 0.05, "{sta_u}");
    if cfg!(feature = "hdf5") {
        assert!(fields, "fields.vtkhdf");
    }
}

// ---------------------------------------------------------------------
// turbulence on the sliding mesh
// ---------------------------------------------------------------------

/// The rows of a sample line of a finished run (column name → value).
fn line_rows(
    case: &FlowCase,
    name: &str,
    line: &str,
) -> Vec<std::collections::HashMap<String, f64>> {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    assert!(rep.failures.is_empty(), "{name}: {:?}", rep.failures);
    let text = std::fs::read_to_string(dir.join(format!("line_{line}.csv"))).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    let mut rows = text.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<String> = rows.next().unwrap().split(',').map(String::from).collect();
    rows.map(|l| {
        head.iter()
            .cloned()
            .zip(l.split(',').map(|v| v.parse::<f64>().unwrap()))
            .collect()
    })
    .collect()
}

/// The largest difference of column `q` between two lines, relative to the
/// largest |q| of the first.
fn deviation(
    a: &[std::collections::HashMap<String, f64>],
    b: &[std::collections::HashMap<String, f64>],
    q: &str,
) -> f64 {
    let scale = a.iter().map(|r| r[q].abs()).fold(0.0, f64::max);
    a.iter()
        .zip(b)
        .map(|(x, y)| (x[q] - y[q]).abs())
        .fold(0.0, f64::max)
        / scale
}

/// Turbulent Taylor–Couette flow (k-ω SST, Re = Ω r₁ (r₂ − r₁)/ν = 12 500)
/// spun up from rest for 2 s: on the conformal sliding mesh and on the
/// imprinted one (32 rotor and 48 stator sectors) the profiles of u_θ, k
/// and ν_t across the gap are those of the fixed mesh with the turning wall,
/// to the discretisation error (1.6 %, 1.5 % and 4 % of the largest value
/// on the line here; about 0.5 %, 1 % and 1.5 % on the mesh of half the
/// width and half the step).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn turbulent_spin_up_agrees_with_a_fixed_mesh() {
    let t = text()
        .replace(
            "n_r = 8\nrotor_layers = 4\n",
            "n_r = 12\nrotor_layers = 6\nradial_ratio = 4.0\nsymmetric = true\n",
        )
        .replace("viscosity = 5e-4", "viscosity = 2e-6")
        .replace(
            "[flow]\n",
            "[turbulence]\nmodel = \"k-omega-sst\"\nwall_treatment = \"automatic\"\n\n\
             [initial]\nturbulence_intensity = 0.05\nviscosity_ratio = 10.0\n\n[flow]\n",
        )
        .replace("steps = 510", "steps = 102")
        .replace(
            "[solver]\ntolerance = 1e-10\n",
            "[solver]\ntolerance = 1e-9\nextrapolation = false\nouter_iterations = 4\n",
        )
        .replace("[verification]\nmax_error = 0.02\n", "")
        .replace(
            "[output]\n",
            "[output]\nlines = [{ name = \"gap\", from = [0.0505, 0.0, 0.005], \
             to = [0.0995, 0.0, 0.005], points = 40 }]\n",
        );
    let mut base = parse(&t);
    base.analytic = None;
    assert!(!base.turbulence.is_laminar() && base.output.lines.len() == 1);
    let mut fixed = base.clone();
    fixed.sliding.clear();
    let mut imprinted = base.clone();
    imprinted.mesh.annulus.as_mut().unwrap().n_theta_stator = 24;
    imprinted.sliding[0].rotor_interface = Some("rotor_interface".into());
    imprinted.sliding[0].stator_interface = Some("stator_interface".into());
    let f = line_rows(&fixed, "turb_fixed", "gap");
    for (name, case) in [("turb_conformal", &base), ("turb_imprinted", &imprinted)] {
        let s = line_rows(case, name, "gap");
        let d = [
            deviation(&f, &s, "v"),
            deviation(&f, &s, "k"),
            deviation(&f, &s, "nu_t"),
        ];
        eprintln!("{name}: largest deviation from the fixed mesh: u_θ, k, ν_t {d:?}");
        assert!(d[0] < 0.02 && d[1] < 0.03 && d[2] < 0.08, "{name}: {d:?}");
    }
}

/// A patch of turbulence (k-ε-ζ-f) in a fluid in solid-body rotation with
/// the walls, carried through half a turn: on the conformal sliding mesh
/// (all cells turn, the old levels of k, ε and ζ move with the cells) and on
/// the imprinted one (the patch in the rotor) it arrives at the opposite
/// side as the same patch decays at rest, and leaves its start empty. The
/// patch lags the turn by the velocity error of the coarse mesh (3 % here);
/// its place on a line is the centroid of k above the line's least value.
/// The imprinted run takes half the time step (see below).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn turbulence_is_carried_by_the_rotor() {
    let blob = "1e-6 + 1e-3 * exp(-((x - 0.065)^2 + y^2) / 2.5e-5)";
    let lines = "[output]\nlines = [\
        { name = \"start\", from = [0.065, -0.02, 0.005], to = [0.065, 0.02, 0.005], points = 41 }, \
        { name = \"half\", from = [-0.065, 0.02, 0.005], to = [-0.065, -0.02, 0.005], points = 41 }]\n";
    let turb = format!(
        "[turbulence]\nmodel = \"k-epsilon-zeta-f0\"\nwall_treatment = \"compound-2023\"\n\n\
         [initial]\nfrom_analytic = true\nk = \"{blob}\"\nepsilon = \"0.9 * ({blob})\"\n\n[flow]\n"
    );
    let solid =
        |t: String| -> String {
            t.replace(
            "[boundary.outer]\ntype = \"wall\"\n",
            "[boundary.outer]\ntype = \"wall\"\nvelocity = [\"-spin * y\", \"spin * x\", 0.0]\n",
        )
        .replace("inner_omega = 10.0", "inner_omega = 10.0\nouter_omega = 10.0")
        .replace("[flow]\n", &turb)
        .replace("[output]\n", lines)
        .replace("[verification]\nmax_error = 0.02\n", "")
        .replace("[verification]\nmax_error = 0.03\n", "")
        };
    let conformal = parse(
        &solid(text())
            .replace("zone = \"rotor\"\n", "")
            .replace("rotor_layers = 4\n", "")
            .replace("steps = 510", "steps = 16"),
    );
    // half the step of the conformal run (Courant number 0.85): at 1.7 the
    // rotor of the imprinted mesh drives ζ of k-ε-ζ-f₀ to its cap and the
    // run diverges in step 8 with any sink of k that is ε at the wall
    // (`"wall"` since #141, the papers' `"epsilon"`); k/T_K
    // (`"kolmogorov"`) held at 1.7 (finding of #141)
    let imprinted = parse(&solid(imprint_text()).replace(
        "time_step = 0.02\nsteps = 500",
        "time_step = 0.009817477042468103\nsteps = 32",
    ));
    let mut rest = conformal.clone();
    rest.sliding.clear();
    rest.analytic = None;
    rest.initial.from_analytic = false;
    for w in ["inner", "outer"] {
        if let Some(foehncraft_flow::bc::BoundaryCondition::Wall { velocity, .. }) =
            rest.boundary.get_mut(w)
        {
            *velocity = foehncraft_flow::expr::vector([0.0; 3]);
        }
    }
    // the largest k and the centroid of k above its least value on the line:
    // the decayed patch is a ring (a dip at its centre), so two cells share
    // the largest k and the place of the largest value is a tie that the last
    // digit decides; the centroid is the place of the patch
    let peak = |rows: &[std::collections::HashMap<String, f64>]| -> (f64, f64) {
        let low = rows.iter().map(|r| r["k"]).fold(f64::INFINITY, f64::min);
        let high = rows.iter().map(|r| r["k"]).fold(0.0, f64::max);
        let (m, w) = rows.iter().fold((0.0, 0.0), |a, r| {
            (a.0 + r["s"] * (r["k"] - low), a.1 + r["k"] - low)
        });
        (if w > 0.0 { m / w } else { 0.0 }, high)
    };
    let (s0, k0) = peak(&line_rows(&rest, "blob_rest", "start"));
    for (name, case) in [
        ("blob_conformal", &conformal),
        ("blob_imprinted", &imprinted),
    ] {
        let (s, k) = peak(&line_rows(case, name, "half"));
        let (_, left) = peak(&line_rows(case, name, "start"));
        eprintln!(
            "{name}: k at rest {k0:.3e} at s = {s0:.4}, turned {k:.3e} at s = {s:.4}, left at the \
             start {left:.2e}"
        );
        assert!((k / k0 - 1.0).abs() < 0.25, "{name}: {k} vs {k0}");
        assert!((s - s0).abs() < 0.01, "{name}: at {s} instead of {s0}");
        assert!(left < 0.1 * k0, "{name}: {left} stayed at the start");
    }
}

// ---------------------------------------------------------------------
// several rotors by imprinting (D185)
// ---------------------------------------------------------------------

fn two_rotors_text() -> String {
    std::fs::read_to_string(cases().join("couette_two_rotors.toml"))
        .unwrap()
        .replace("\r\n", "\n")
}

/// The checks of several rotors: each needs the interface by imprinting
/// and its own zone, and a zone is turned by a frame or a sliding mesh.
#[test]
fn several_rotors_are_checked() {
    let valid = |t: &str| -> Result<(), String> {
        FlowCase::from_toml(t, &cases())
            .and_then(|c| c.validate())
            .map_err(|e| code(&e))
    };
    let base = two_rotors_text();
    valid(&base).unwrap();
    // two rotors in conformal stepping
    let conformal = base
        .replace("rotor_interface = \"rotor_interface_outer\"\n", "")
        .replace("stator_interface = \"stator_interface_outer\"\n", "");
    assert_eq!(valid(&conformal).unwrap_err(), "SLD002");
    // the same zone twice
    let twice = base.replace("zone = \"rotor_outer\"", "zone = \"rotor\"");
    assert_eq!(valid(&twice).unwrap_err(), "SLD002");
    // a frame on the cells of a sliding rotor, or on all cells
    let frame = format!("{base}\n[[mrf]]\nzone = \"rotor_outer\"\nomega = 1.0\n");
    assert_eq!(valid(&frame).unwrap_err(), "MRF001");
    let all = format!("{base}\n[[mrf]]\nomega = 1.0\n");
    assert_eq!(valid(&all).unwrap_err(), "MRF001");
}

/// Couette flow between two cylinders that turn in opposite senses, each
/// a sliding rotor by imprinting with its own speed (D185), from rest: the
/// end is the exact solution u_θ = A r + B/r, the torque of the fluid on
/// the inner cylinder is 4π μ L r₁² r₂² (ω₂ − ω₁)/(r₂² − r₁²) and that on the
/// outer one is opposite (each is reported about the axis along its own
/// Ω, so both are negative: both rotors drive the fluid), the shaft powers
/// add up to the dissipation, and the angle of each rotor is its own.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn two_rotors_end_in_the_couette_flow_of_their_speeds() {
    let c = parse(&two_rotors_text());
    let r = run_case(&c, "two_rotors");
    let e = l2(&r);
    let s = &r.solves[0];
    let v = s.rotors.as_ref().expect("rotors");
    let (inner, outer) = (v.rotor("rotor").unwrap(), v.rotor("rotor_outer").unwrap());
    let (r1, r2, depth, mu, w1, w2) = (0.05_f64, 0.10_f64, 0.01_f64, 5e-4_f64, 10.0, -5.0);
    let t1 = 4.0 * std::f64::consts::PI * mu * depth * r1 * r1 * r2 * r2 * (w2 - w1)
        / (r2 * r2 - r1 * r1);
    eprintln!(
        "two rotors: L2 {e:e}; torque inner {:e}, outer {:e} (exact {t1:e} each about its own Ω); power {:e} + {:e}",
        inner.torque, outer.torque, inner.power, outer.power
    );
    assert!(e < 0.03, "{e}");
    assert_eq!(inner.walls, ["inner"]);
    assert_eq!(outer.walls, ["outer"]);
    assert!((inner.torque / t1 - 1.0).abs() < 0.02, "{}", inner.torque);
    assert!((outer.torque / t1 - 1.0).abs() < 0.02, "{}", outer.torque);
    // the shaft powers are |Ω| times the torques, and their sum is the work
    // done on the fluid, which the steady state dissipates
    assert!((inner.power + w1 * inner.torque).abs() < 1e-12 * inner.power);
    assert!((outer.power + w2.abs() * outer.torque).abs() < 1e-12 * outer.power);
    let exact_power = -t1 * (w1 - w2);
    let p = inner.power + outer.power;
    assert!((p / exact_power - 1.0).abs() < 0.02, "{p} {exact_power}");
    let b = s.balance.as_ref().unwrap();
    let k = &b.last.kinetic_energy;
    // (the ledger of the last step against the values at its end)
    assert!((k.shaft / p - 1.0).abs() < 1e-4, "{} {p}", k.shaft);
    assert!(
        (k.viscous_dissipation / p - 1.0).abs() < 0.05,
        "{} {p}",
        k.viscous_dissipation
    );
}

/// Interfaces of any rotationally symmetric shape (D186): uniform axial flow
/// between slip walls through a turning rotor.
///
/// A pot whose two discs pass through the axis, away from every wall: the
/// uniform flow is exact on every imprinted mesh (the cells close, the fluxes
/// through the pieces are conservative), to the level of the skewed wedge
/// cells at the axis. The stage (disc, cylinder, disc on the hub) and its
/// cones: where the interface meets the hub, the hub of the rotor is a
/// polygon of 24 sides and the hub of the stator one of 34, and the small
/// step between them stops a little of the flow along the hub; the error is
/// that of the polygons and falls with the square of the pitch.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn uniform_flow_crosses_discs_cones_and_the_axis() {
    let t = std::fs::read_to_string(cases().join("axial_stage_uniform.toml"))
        .unwrap()
        .replace("\r\n", "\n");
    let base = parse(&t);
    let pot = parse(
        &t.replace("r = [0.05, 0.09, 0.10]", "r = [0.0, 0.04, 0.10]")
            .replace("[boundary.inner]\ntype = \"slip-wall\"\n", "")
            .replace("[boundary.inner_rotor]\ntype = \"slip-wall\"\n", ""),
    );
    let e = l2(&run_case(&pot, "uniform_pot"));
    eprintln!("uniform flow through the pot: L2 {e:e}");
    assert!(e < 1e-5, "pot: {e}");
    for cone in [0.0, 20.0] {
        let mut errors = Vec::new();
        for f in [1, 2] {
            let mut c = base.clone();
            let m = c.mesh.revolved.as_mut().unwrap();
            m.cone = cone;
            m.n_theta *= f;
            m.n_theta_stator *= f;
            c.verification.max_error = None;
            errors.push(l2(&run_case(&c, &format!("uniform_stage_{f}"))));
        }
        eprintln!("uniform flow through the stage, cone {cone}: L2 {errors:?}");
        assert!(errors[0] < 3e-3, "cone {cone}: {errors:?}");
        assert!(errors[0] / errors[1] > 2.8, "cone {cone}: {errors:?}");
    }
}

/// Solid-body rotation in a closed vessel with a rotor pot whose interface
/// passes through the axis (D186): every wall turns with the rotor, and the
/// imprinted meshes keep the exact flow through half a revolution.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn solid_body_rotation_is_kept_by_a_pot_through_the_axis() {
    let t = std::fs::read_to_string(cases().join("pot_solid_body.toml"))
        .unwrap()
        .replace("\r\n", "\n");
    let r = run_case(&parse(&t), "pot_solid");
    let e = l2(&r);
    eprintln!("pot solid body, L2 {e:e}");
    assert!(e < 0.03, "{e}");
}

/// The same at ten times the viscosity: the wedge cells at the axis couple
/// to their neighbours around the axis about 90 times more strongly than
/// to their old level. With Issa's corrector (D = V/a_P) the pressure grew
/// in every outer iteration there and the run left the exact solution
/// within 30 steps (L2 error 2.0, checked below); the SIMPLEC-type
/// corrector keeps it (gallery 103).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn solid_body_rotation_is_kept_at_high_viscosity_on_the_axis_wedges() {
    let t = std::fs::read_to_string(cases().join("pot_solid_body_viscous.toml"))
        .unwrap()
        .replace("\r\n", "\n");
    let r = run_case(&parse(&t), "pot_solid_viscous");
    let e = l2(&r);
    eprintln!("pot solid body at nu = 1e-2, L2 {e:e}");
    assert!(e < 0.01, "{e}");
    // Issa's corrector forced: the axisymmetric pressure mode of the wedges
    // grows from round-off
    let mut c = parse(&t);
    c.solver.corrector = foehncraft_flow::simplec::PisoCorrector::Issa;
    c.verification.max_error = None;
    let issa = try_run(&c, "pot_solid_viscous_issa").map(|r| l2(&r));
    eprintln!("the same with Issa's corrector: {issa:?}");
    assert!(!matches!(issa, Ok(e) if e < 0.1), "{issa:?}");
}
