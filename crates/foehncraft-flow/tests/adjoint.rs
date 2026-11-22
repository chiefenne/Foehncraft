//! Adjoint-ready kernels (requirement D59, stage 5): the kernel bodies of
//! the scalar operators and of the flow solver run on dual numbers, and
//! their exact derivatives agree with central finite differences of the
//! same bodies on `f64`.
//!
//! Every input array of a kernel (geometry, states, gradients, limiters,
//! boundary values, properties) gets a random direction; the test compares
//! the directional derivative of every output of every face (or cell, or
//! boundary face) once for each input array alone and once for all of them
//! together. The data are random on a perturbed mixed mesh of hexahedra,
//! prisms, tetrahedra and polyhedra, with every boundary kind and update
//! rule present.
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use foehncraft_core::{Dual, Real};
use foehncraft_fv::Discretisation;
use foehncraft_fv::problem::{Boundary, Problem};
use foehncraft_kernels::flow::{
    self as k, BoundaryInputs, FaceGeometry, FluxInputs, MomentumInputs, VelocityGradient,
};
use foehncraft_kernels::fv::{self, FaceInputs, Topology, Vec3};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::case::{MeshSection, build_mesh};
use std::collections::BTreeMap;

type D = Dual<1>;

/// Step of the central differences (relative to the scale of each array).
const H: f64 = 1e-6;
/// Allowed difference between the dual and the finite-difference
/// derivative, relative to the larger of the two (or to 10⁻³ of the
/// output, below which the central difference is round-off).
const TOL: f64 = 2e-6;

/// A small random-number generator (xorshift64*), so the data do not
/// depend on another crate.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let x = self.0.wrapping_mul(0x2545_F491_4F6C_DD1D);
        (x >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Uniform in [a, b).
    fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.next()
    }
    /// ± a value in [a, b).
    fn signed(&mut self, a: f64, b: f64) -> f64 {
        let v = self.range(a, b);
        if self.next() < 0.5 { -v } else { v }
    }
    fn pick(&mut self, n: u32) -> u32 {
        ((self.next() * n as f64) as u32).min(n - 1)
    }
}

/// The mesh, its discretisation and the boundary kinds.
struct Fixture {
    disc: Discretisation,
    mom_kind: Vec<u32>,
    flux_kind: Vec<u32>,
    ub_rule: Vec<u32>,
    pb_rule: Vec<u32>,
}

impl Fixture {
    fn topo(&self) -> Topology<'_> {
        self.disc.topology()
    }
}

/// All input arrays by name, in `f64`.
type Data = BTreeMap<String, Vec<f64>>;

/// The input arrays of type `T`: value `x + s h d` for finite differences
/// (`s` = ±1) or the dual number `x + d ε`.
struct Lifted<T>(BTreeMap<String, Vec<T>>);

impl<T: Real> Lifted<T> {
    fn s(&self, n: &str) -> &[T] {
        self.0.get(n).unwrap_or_else(|| panic!("no array {n}"))
    }
    fn v3(&self, n: &str) -> Vec3<'_, T> {
        Vec3 {
            x: self.s(&format!("{n}.x")),
            y: self.s(&format!("{n}.y")),
            z: self.s(&format!("{n}.z")),
        }
    }
    fn grad(&self, n: &str) -> VelocityGradient<'_, T> {
        VelocityGradient {
            x: self.v3(&format!("{n}x")),
            y: self.v3(&format!("{n}y")),
            z: self.v3(&format!("{n}z")),
        }
    }
    fn geo(&self) -> FaceGeometry<'_, T> {
        FaceGeometry {
            area: self.v3("area"),
            face_centroid: self.v3("fc"),
            cell_centroid: self.v3("cc"),
            dgeo: self.s("dgeo"),
            k: self.v3("k"),
            weight: self.s("weight"),
            shift: Vec3 {
                x: &[],
                y: &[],
                z: &[],
            },
        }
    }
}

fn fixture() -> (Fixture, Data, Data) {
    let sec: MeshSection = toml::from_str(
        r#"
[generate]
type = "mixed"
dimension = 3
lower = [0.0, 0.0, 0.0]
upper = [1.0, 1.0, 1.0]
cells = [4, 4, 3]
perturb = 0.2
seed = 59
"#,
    )
    .unwrap();
    let (mesh, _, _) = build_mesh(&sec).unwrap();
    let geo = Geometry::compute(&mesh);
    let mut prob = Problem::default();
    for (i, p) in mesh.boundaries.iter().enumerate() {
        let b = if i % 2 == 0 {
            Boundary::value(1.0)
        } else {
            Boundary::zero_gradient()
        };
        prob.boundary.insert(p.name.clone(), b);
    }
    let disc = Discretisation::new(&mesh, &geo, &prob).unwrap();
    let (nc, nf, ni) = (disc.n_cells, disc.n_faces, disc.n_interior);
    let nb = nf - ni;
    let mut r = Rng(0x05ee_dd59);
    let fx = Fixture {
        mom_kind: (0..nb).map(|_| r.pick(3)).collect(),
        flux_kind: (0..nb).map(|_| r.pick(3)).collect(),
        ub_rule: (0..nb).map(|_| r.pick(7)).collect(),
        pb_rule: (0..nb).map(|_| r.pick(6)).collect(),
        disc,
    };
    let d = &fx.disc;
    let mut x = Data::new();
    let put3 = |x: &mut Data, n: &str, v: [&Vec<f64>; 3]| {
        for (c, a) in ["x", "y", "z"].iter().zip(v) {
            x.insert(format!("{n}.{c}"), a.clone());
        }
    };
    put3(
        &mut x,
        "area",
        [&d.face_area.x, &d.face_area.y, &d.face_area.z],
    );
    put3(
        &mut x,
        "fc",
        [&d.face_centroid.x, &d.face_centroid.y, &d.face_centroid.z],
    );
    put3(
        &mut x,
        "cc",
        [&d.cell_centroid.x, &d.cell_centroid.y, &d.cell_centroid.z],
    );
    put3(&mut x, "k", [&d.k.x, &d.k.y, &d.k.z]);
    put3(&mut x, "lsq", [&d.lsq.x, &d.lsq.y, &d.lsq.z]);
    x.insert("dgeo".into(), d.dgeo.clone());
    x.insert("weight".into(), d.weight.clone());
    x.insert("bca".into(), d.bc_a.clone());
    let mut rand = |n: usize, a: f64, b: f64, signed: bool| -> Vec<f64> {
        (0..n)
            .map(|_| {
                if signed {
                    r.signed(a, b)
                } else {
                    r.range(a, b)
                }
            })
            .collect()
    };
    // mass fluxes bounded away from zero, so no upwind switch lies within
    // the finite-difference step
    x.insert("m".into(), rand(nf, 0.2, 1.0, true));
    x.insert("mu".into(), rand(nf, 0.5, 1.5, false));
    x.insert("gamma".into(), rand(nf, 0.5, 1.5, false));
    x.insert("relax".into(), rand(nf, -0.1, 0.1, false));
    x.insert("time".into(), rand(nf, -0.1, 0.1, false));
    x.insert("cf".into(), rand(nf, 0.5, 1.5, false));
    for n in ["phi", "p", "pc", "d", "ds", "dc"] {
        let signed = matches!(n, "phi" | "p" | "pc");
        x.insert(n.into(), rand(nc, 0.5, 1.5, signed));
    }
    x.insert("psis".into(), rand(nc, 0.3, 1.0, false));
    x.insert("eps2".into(), rand(nc, 1e-3, 1e-2, false));
    for n in ["u", "gp", "psi", "delta"] {
        for c in ["x", "y", "z"] {
            let v = if n == "psi" {
                rand(nc, 0.3, 1.0, false)
            } else {
                rand(nc, 0.2, 1.0, true)
            };
            x.insert(format!("{n}.{c}"), v);
        }
    }
    for i in ["x", "y", "z"] {
        for c in ["x", "y", "z"] {
            x.insert(format!("gu{i}.{c}"), rand(nc, 0.2, 1.0, true));
        }
    }
    for c in ["x", "y", "z"] {
        x.insert(format!("ub.{c}"), rand(nb, 0.2, 1.0, true));
    }
    x.insert("pb".into(), rand(nb, 0.5, 1.5, true));
    x.insert("pset".into(), rand(nb, 0.5, 1.5, true));
    x.insert("bcb".into(), rand(nb, 0.5, 1.5, true));
    // the convective outlet: the old-level part c and the factor beta
    for c in ["x", "y", "z"] {
        x.insert(format!("convc.{c}"), rand(nb, 0.2, 1.0, true));
    }
    x.insert("beta".into(), rand(nb, 0.1, 2.0, false));
    // the face density (D150): a field, so its derivative is checked too
    x.insert("rho_face".into(), rand(nf, 1.0, 1.4, false));
    // closed (suppressed backflow) is 0 or 1, a switch
    x.insert("closed".into(), (0..nb).map(|i| (i % 2) as f64).collect());
    // random directions, scaled by the size of each array (geometry
    // arrays are small on a fine mesh)
    let mut dir = Data::new();
    for (n, v) in &x {
        let scale = (v.iter().map(|a| a * a).sum::<f64>() / v.len().max(1) as f64).sqrt();
        let scale = if scale > 0.0 { scale } else { 1.0 };
        dir.insert(
            n.clone(),
            (0..v.len()).map(|_| scale * r.signed(0.1, 1.0)).collect(),
        );
    }
    // a is 0 or 1 (a switch, not a value): no direction
    dir.get_mut("bca").unwrap().fill(0.0);
    dir.get_mut("closed").unwrap().fill(0.0);
    (fx, x, dir)
}

/// The inputs moved by `s h d` (only the arrays in `only`, or all).
fn shifted(x: &Data, dir: &Data, s: f64, only: Option<&str>) -> Lifted<f64> {
    Lifted(
        x.iter()
            .map(|(n, v)| {
                let on = only.is_none_or(|o| o == n);
                let d = &dir[n];
                let w = v
                    .iter()
                    .zip(d)
                    .map(|(a, b)| if on { a + s * H * b } else { *a })
                    .collect();
                (n.clone(), w)
            })
            .collect(),
    )
}

/// The inputs as dual numbers with derivative `d` (only the arrays in
/// `only`, or all).
fn dual(x: &Data, dir: &Data, only: Option<&str>) -> Lifted<D> {
    Lifted(
        x.iter()
            .map(|(n, v)| {
                let on = only.is_none_or(|o| o == n);
                let w = v
                    .iter()
                    .zip(&dir[n])
                    .map(|(&a, &b)| D::new(a, [if on { b } else { 0.0 }]))
                    .collect();
                (n.clone(), w)
            })
            .collect(),
    )
}

/// The outputs of a kernel for index `i`.
type Eval<T> = fn(&Fixture, &Lifted<T>, usize) -> Vec<T>;

/// Compares dual and finite-difference derivatives of `kernel` over the
/// indices `0..n`, for each input array alone and for all together.
/// Returns the number of derivatives compared and the largest relative
/// difference.
fn check(
    name: &str,
    n: usize,
    fx: &Fixture,
    x: &Data,
    dir: &Data,
    f: Eval<f64>,
    fd: Eval<D>,
) -> (usize, f64) {
    let mut sets: Vec<Option<&str>> = x.keys().map(|k| Some(k.as_str())).collect();
    sets.push(None);
    let mut count = 0;
    let mut worst = 0.0f64;
    let l0 = shifted(x, dir, 0.0, None);
    let lall = dual(x, dir, None);
    for only in sets {
        let (lp, lm, ld) = (
            shifted(x, dir, 1.0, only),
            shifted(x, dir, -1.0, only),
            dual(x, dir, only),
        );
        for i in 0..n {
            let (yp, ym, yd) = (f(fx, &lp, i), f(fx, &lm, i), fd(fx, &ld, i));
            let y0 = f(fx, &l0, i);
            // the dual value is the f64 value (equal up to the sign of zero,
            // which f64::max leaves open)
            for (a, b) in yd.iter().zip(&y0) {
                assert!(a.v == *b, "{name}: value of index {i}: {} against {b}", a.v);
            }
            // the round-off of the central difference is about
            // ε Σ|terms| / h, some 1e-10 of the output or of its derivative
            // in all directions (when the output itself cancels)
            let big = |v: &mut dyn Iterator<Item = f64>| v.fold(1e-300f64, |m, a| m.max(a.abs()));
            let floor = (1e-3 * big(&mut y0.iter().copied()))
                .max(1e-4 * big(&mut fd(fx, &lall, i).iter().map(|v| v.d[0])));
            for (j, (a, (p, m))) in yd.iter().zip(yp.iter().zip(&ym)).enumerate() {
                let b = (p - m) / (2.0 * H);
                let scale = a.d[0].abs().max(b.abs()).max(floor);
                let e = (a.d[0] - b).abs() / scale;
                worst = worst.max(e);
                assert!(
                    e < TOL,
                    "{name}, index {i}, output {j}, direction {}: dual {:e}, central difference {b:e}",
                    only.unwrap_or("all"),
                    a.d[0]
                );
                count += 1;
            }
        }
    }
    println!("{name}: {count} derivatives, largest relative difference {worst:.2e}");
    (count, worst)
}

fn gradient<T: Real>(fx: &Fixture, l: &Lifted<T>, c: usize) -> Vec<T> {
    let bc = fv::Boundary {
        a: l.s("bca"),
        b: l.s("bcb"),
        jump: &[],
    };
    // SAFETY: c is a cell of the mesh; all arrays have their full length.
    unsafe { fv::lsq_gradient(c, &fx.topo(), &l.v3("lsq"), l.s("phi"), &bc).to_vec() }
}

fn limiter<T: Real>(fx: &Fixture, l: &Lifted<T>, c: usize) -> Vec<T> {
    let bc = fv::Boundary {
        a: l.s("bca"),
        b: l.s("bcb"),
        jump: &[],
    };
    // SAFETY: as above.
    unsafe {
        vec![fv::venkat_limiter(
            c,
            &fx.topo(),
            l.s("phi"),
            &bc,
            &l.v3("gp"),
            &l.v3("cc"),
            &l.v3("fc"),
            l.s("eps2"),
        )]
    }
}

fn scalar_face<T: Real>(fx: &Fixture, l: &Lifted<T>, f: usize) -> Vec<T> {
    let inp = FaceInputs {
        topo: fx.topo(),
        bc: fv::Boundary {
            a: l.s("bca"),
            b: l.s("bcb"),
            jump: &[],
        },
        mass_flux: l.s("m"),
        gamma: l.s("gamma"),
        dgeo: l.s("dgeo"),
        k: l.v3("k"),
        weight: l.s("weight"),
        face_centroid: l.v3("fc"),
        cell_centroid: l.v3("cc"),
        shift: Vec3 {
            x: &[],
            y: &[],
            z: &[],
        },
        phi: l.s("phi"),
        grad: l.v3("gp"),
        psi: l.s("psis"),
        blend: T::ONE,
        central: T::ZERO,
        beta: T::ZERO,
        nonorth: T::ONE,
        nonorth_limit: T::from_f64(-1.0),
    };
    // SAFETY: f is a face of the mesh.
    let (a, b, c) = unsafe { fv::face_coefficients(f, &inp) };
    vec![a, b, c]
}

fn momentum_face<T: Real>(fx: &Fixture, l: &Lifted<T>, f: usize) -> Vec<T> {
    momentum_face_with(fx, l, f, T::ZERO)
}

/// The blended central scheme (D25), unbounded (β = 0, a smooth function
/// of the inputs).
fn momentum_face_central<T: Real>(fx: &Fixture, l: &Lifted<T>, f: usize) -> Vec<T> {
    momentum_face_with(fx, l, f, T::from_f64(0.7))
}

fn momentum_face_with<T: Real>(fx: &Fixture, l: &Lifted<T>, f: usize, central: T) -> Vec<T> {
    let inp = MomentumInputs {
        topo: fx.topo(),
        geo: l.geo(),
        mass_flux: l.s("m"),
        mu: l.s("mu"),
        u: l.v3("u"),
        uc: l.v3("u"),
        grad: l.grad("gu"),
        psi: l.v3("psi"),
        ub: l.v3("ub"),
        mom_kind: &fx.mom_kind,
        blend: T::ONE,
        central,
        beta: T::ZERO,
        nonorth: T::ONE,
        convection: T::ONE,
        transposed: T::ZERO,
        trace_free: T::ZERO,
    };
    // SAFETY: f is a face of the mesh.
    let (cp, cn, b) = unsafe { k::momentum_face(f, &inp) };
    vec![cp, cn, b[0], b[1], b[2]]
}

fn flux_inputs<'a, T: Real>(
    fx: &'a Fixture,
    l: &'a Lifted<T>,
    transient: bool,
) -> FluxInputs<'a, T> {
    FluxInputs {
        topo: fx.topo(),
        geo: l.geo(),
        rho_face: l.s("rho_face"),
        u: l.v3("u"),
        grad: l.grad("gu"),
        p: l.s("p"),
        grad_p: l.v3("gp"),
        pb: l.s("pb"),
        pset: l.s("pset"),
        ub: l.v3("ub"),
        flux_kind: &fx.flux_kind,
        mom_kind: &fx.mom_kind,
        d: l.s("d"),
        relax: l.s("relax"),
        transient,
        ds: l.s("ds"),
        time: l.s("time"),
        alpha: T::from_f64(0.7),
        psi: l.v3("u"),
        bounded: false,
        jump_offset: &[],
    }
}

fn mass_flux<T: Real>(fx: &Fixture, l: &Lifted<T>, f: usize) -> Vec<T> {
    // SAFETY: f is a face of the mesh.
    unsafe { vec![k::face_mass_flux(f, &flux_inputs(fx, l, false))] }
}

fn mass_flux_transient<T: Real>(fx: &Fixture, l: &Lifted<T>, f: usize) -> Vec<T> {
    // SAFETY: as above.
    unsafe { vec![k::face_mass_flux(f, &flux_inputs(fx, l, true))] }
}

fn pressure_face<T: Real>(fx: &Fixture, l: &Lifted<T>, f: usize) -> Vec<T> {
    let geo = l.geo();
    // SAFETY: as above.
    unsafe {
        vec![k::pressure_face(
            f,
            &fx.topo(),
            &geo,
            &fx.flux_kind,
            l.s("rho_face"),
            l.s("dc"),
            l.s("closed"),
        )]
    }
}

fn correct_flux<T: Real>(fx: &Fixture, l: &Lifted<T>, f: usize) -> Vec<T> {
    // SAFETY: as above.
    unsafe {
        vec![k::correct_flux(
            f,
            &fx.topo(),
            l.s("cf"),
            l.s("m"),
            l.s("pc"),
        )]
    }
}

fn delta_flux<T: Real>(fx: &Fixture, l: &Lifted<T>, f: usize) -> Vec<T> {
    let geo = l.geo();
    // SAFETY: as above.
    unsafe {
        vec![k::delta_flux(
            f,
            &fx.topo(),
            &geo,
            &fx.flux_kind,
            l.s("rho_face"),
            &l.v3("delta"),
            l.s("closed"),
        )]
    }
}

fn boundary_inputs<'a, T: Real>(fx: &'a Fixture, l: &'a Lifted<T>) -> BoundaryInputs<'a, T> {
    BoundaryInputs {
        topo: fx.topo(),
        geo: l.geo(),
        rho_face: l.s("rho_face"),
        mass_flux: l.s("m"),
    }
}

fn boundary_velocity<T: Real>(fx: &Fixture, l: &Lifted<T>, fb: usize) -> Vec<T> {
    let ub = l.v3("ub");
    // SAFETY: fb is a boundary face.
    unsafe {
        k::boundary_velocity(
            fb,
            &boundary_inputs(fx, l),
            fx.ub_rule[fb],
            ub.at(fb),
            (l.v3("convc").at(fb), l.s("beta")[fb]),
            &l.v3("u"),
            &l.grad("gu"),
        )
        .to_vec()
    }
}

fn boundary_pressure<T: Real>(fx: &Fixture, l: &Lifted<T>, fb: usize) -> Vec<T> {
    // SAFETY: as above.
    unsafe {
        vec![k::boundary_pressure(
            fb,
            &boundary_inputs(fx, l),
            fx.pb_rule[fb],
            l.s("pset")[fb],
            l.s("p"),
            &l.v3("gp"),
            &[],
        )]
    }
}

macro_rules! jacobian_test {
    ($test:ident, $kernel:ident, $n:expr) => {
        #[test]
        fn $test() {
            let (fx, x, dir) = fixture();
            let n = $n(&fx.disc);
            let (count, _) = check(
                stringify!($kernel),
                n,
                &fx,
                &x,
                &dir,
                $kernel::<f64>,
                $kernel::<D>,
            );
            assert!(count > 0);
        }
    };
}

fn cells(d: &Discretisation) -> usize {
    d.n_cells
}
fn faces(d: &Discretisation) -> usize {
    d.n_faces
}
fn boundary_faces(d: &Discretisation) -> usize {
    d.n_faces - d.n_interior
}

jacobian_test!(jacobian_of_the_least_squares_gradient, gradient, cells);
jacobian_test!(jacobian_of_the_venkatakrishnan_limiter, limiter, cells);
jacobian_test!(jacobian_of_the_scalar_face_flux, scalar_face, faces);
jacobian_test!(jacobian_of_the_momentum_face_flux, momentum_face, faces);
jacobian_test!(
    jacobian_of_the_central_momentum_face_flux,
    momentum_face_central,
    faces
);
jacobian_test!(jacobian_of_the_rhie_chow_mass_flux, mass_flux, faces);
jacobian_test!(
    jacobian_of_the_rhie_chow_mass_flux_in_a_time_step,
    mass_flux_transient,
    faces
);
jacobian_test!(
    jacobian_of_the_pressure_correction_coefficient,
    pressure_face,
    faces
);
jacobian_test!(jacobian_of_the_flux_correction, correct_flux, faces);
jacobian_test!(jacobian_of_the_velocity_increment_flux, delta_flux, faces);
jacobian_test!(
    jacobian_of_the_boundary_velocity,
    boundary_velocity,
    boundary_faces
);
jacobian_test!(
    jacobian_of_the_boundary_pressure,
    boundary_pressure,
    boundary_faces
);

/// The fixture covers every boundary kind and update rule.
#[test]
fn fixture_has_every_boundary_kind() {
    let (fx, x, _) = fixture();
    let d = &fx.disc;
    println!(
        "{} cells, {} faces ({} boundary), {} input arrays",
        d.n_cells,
        d.n_faces,
        d.n_faces - d.n_interior,
        x.len()
    );
    for (v, n) in [
        (&fx.mom_kind, 3),
        (&fx.flux_kind, 3),
        (&fx.ub_rule, 7),
        (&fx.pb_rule, 6),
    ] {
        for kind in 0..n {
            assert!(v.contains(&kind), "kind {kind} missing");
        }
    }
    assert!(fx.disc.bc_a.contains(&0.0) && fx.disc.bc_a.contains(&1.0));
}
