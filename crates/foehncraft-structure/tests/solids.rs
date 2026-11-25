//! Verification of the quadratic solids (MS5 part 3b): the ten-node
//! tetrahedron and the twenty-node hexahedron.
//!
//! - Quadrature, shape functions and face orientation.
//! - Rigid motions free of force, exact mass, the patch test on distorted
//!   meshes (Irons & Razzaque 1972): a uniform stress is reproduced to
//!   rounding.
//! - A cantilever against Timoshenko beam theory (deflection and first
//!   frequency).
//! - The thick-walled cylinder under internal pressure against Lamé's
//!   solution (Timoshenko & Goodier, *Theory of Elasticity*, 1970, §28),
//!   with ν = 0.3 and the nearly incompressible ν = 0.49 of the MacNeal &
//!   Harder test set (*Finite Elements in Analysis and Design* 1, 1985).
//! - Free and held thermal expansion.
//! - The Gmsh reader: a block written as a second-order `.msh` 4.1 file
//!   gives the same answer as the generated block.

#![allow(clippy::needless_range_loop)]

use foehncraft_structure::analysis;
use foehncraft_structure::case::{self, Case};
use foehncraft_structure::codes;
use foehncraft_structure::model::DOFS;
use foehncraft_structure::solid::{self, Kind, SolidProps};

fn rel(a: f64, b: f64) -> f64 {
    ((a - b) / b).abs()
}

fn fact(n: u32) -> f64 {
    (1..=n).map(f64::from).product()
}

#[test]
fn quadrature_and_shape_functions() {
    // the tetrahedron rule integrates every monomial of degree ≤ 5
    for a in 0..=5u32 {
        for b in 0..=(5 - a) {
            for c in 0..=(5 - a - b) {
                let sum: f64 = solid::rule(Kind::Tet10)
                    .iter()
                    .map(|(p, w)| {
                        w * p[0].powi(a as i32) * p[1].powi(b as i32) * p[2].powi(c as i32)
                    })
                    .sum();
                let exact = fact(a) * fact(b) * fact(c) / fact(a + b + c + 3);
                assert!(rel(sum, exact) < 1e-13, "x^{a} y^{b} z^{c}: {sum} {exact}");
            }
        }
    }
    for kind in [Kind::Tet10, Kind::Hex20] {
        let nodes = kind.natural_nodes();
        for (j, pj) in nodes.iter().enumerate() {
            let (n, _) = solid::shape(kind, *pj);
            for i in 0..kind.nodes() {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((n[i] - want).abs() < 1e-14, "{kind:?} N{i} at node {j}");
            }
        }
        let (n, d) = solid::shape(kind, [0.21, 0.13, 0.37]);
        let s: f64 = n.iter().sum();
        assert!((s - 1.0).abs() < 1e-14);
        for k in 0..3 {
            let s: f64 = d.iter().map(|g| g[k]).sum();
            assert!(s.abs() < 1e-13);
        }
        // derivatives against differences
        let p = [0.21, 0.13, 0.37];
        let h = 1e-6;
        for k in 0..3 {
            let mut a = p;
            let mut b = p;
            a[k] += h;
            b[k] -= h;
            let (na, _) = solid::shape(kind, a);
            let (nb, _) = solid::shape(kind, b);
            for i in 0..kind.nodes() {
                let fd = (na[i] - nb[i]) / (2.0 * h);
                assert!((fd - d[i][k]).abs() < 1e-8, "{kind:?} dN{i}/dξ{k}");
            }
        }
        // the faces point out of the element
        let c: [f64; 3] = {
            let mut c = [0.0; 3];
            for p in &nodes[..kind.corners()] {
                for k in 0..3 {
                    c[k] += p[k] / kind.corners() as f64;
                }
            }
            c
        };
        for f in kind.faces() {
            let x: Vec<[f64; 3]> = f.iter().map(|&a| nodes[a]).collect();
            let nrm = solid::face_normal(&x);
            let fc = x[0];
            let out: f64 = (0..3).map(|k| nrm[k] * (fc[k] - c[k])).sum();
            assert!(out > 0.0, "{kind:?} face {f:?}");
        }
    }
}

#[test]
fn face_loads_are_the_textbook_values() {
    // a uniform load on a flat six-node triangle: nothing at the corners,
    // a third at each edge node; on an eight-node quadrilateral −1/12 at
    // the corners and 1/3 at the edge nodes
    let tri = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.5, 0.0],
        [0.0, 0.5, 0.0],
    ];
    let (f, a) = solid::face_load(&tri, [0.0, 0.0, 3.0], 0.0);
    assert!((a - 1.0).abs() < 1e-14);
    for i in 0..6 {
        let want = if i < 3 { 0.0 } else { 3.0 / 3.0 };
        assert!((f[3 * i + 2] - want).abs() < 1e-14, "{i}: {}", f[3 * i + 2]);
    }
    let quad = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 3.0, 0.0],
        [0.0, 3.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 1.5, 0.0],
        [1.0, 3.0, 0.0],
        [0.0, 1.5, 0.0],
    ];
    // a pressure pushes against the normal (here +z)
    let (f, a) = solid::face_load(&quad, [0.0; 3], 2.0);
    assert!((a - 6.0).abs() < 1e-13);
    for i in 0..8 {
        let want = if i < 4 { -6.0 / 12.0 } else { 6.0 / 3.0 };
        assert!(
            (f[3 * i + 2] + 2.0 * want).abs() < 1e-13,
            "{i}: {}",
            f[3 * i + 2]
        );
    }
}

const STEEL: SolidProps = SolidProps {
    young: 200e9,
    poisson: 0.3,
    density: 7800.0,
};

/// A distorted element: the natural nodes moved by a smooth map (so the
/// element stays valid and curved).
fn distorted(kind: Kind) -> Vec<[f64; 3]> {
    kind.natural_nodes()
        .iter()
        .map(|p| {
            [
                1.3 * p[0] + 0.1 * p[1] * p[1] + 0.05 * p[2],
                0.9 * p[1] + 0.1 * p[0] * p[2],
                1.1 * p[2] - 0.08 * p[0] * p[0] + 0.03 * p[1],
            ]
        })
        .collect()
}

#[test]
fn rigid_motions_are_free_and_the_mass_is_exact() {
    for kind in [Kind::Tet10, Kind::Hex20] {
        let x = distorted(kind);
        let (k, m, _) = solid::matrices(kind, &x, &STEEL, 0.0).unwrap();
        let nd = 3 * kind.nodes();
        let kmax = k.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        // three translations and three small rotations
        for mode in 0..6 {
            let u: Vec<f64> = (0..nd)
                .map(|i| {
                    let (a, c) = (i / 3, i % 3);
                    let p = x[a];
                    match mode {
                        0..=2 => f64::from(c == mode),
                        3 => [0.0, -p[2], p[1]][c],
                        4 => [p[2], 0.0, -p[0]][c],
                        _ => [-p[1], p[0], 0.0][c],
                    }
                })
                .collect();
            for i in 0..nd {
                let r: f64 = (0..nd).map(|j| k[i * nd + j] * u[j]).sum();
                assert!(r.abs() < 1e-11 * kmax, "{kind:?} mode {mode}: {r}");
            }
        }
        // the mass of each direction is ρ V
        let v = solid::volume(kind, &x).unwrap();
        let mx: f64 = (0..kind.nodes())
            .flat_map(|a| (0..kind.nodes()).map(move |b| (a, b)))
            .map(|(a, b)| m[(3 * a) * nd + 3 * b])
            .sum();
        assert!(rel(mx, STEEL.density * v) < 1e-13, "{kind:?}");
        // an inverted element is refused
        let mut bad = x.clone();
        for p in &mut bad {
            p[2] = -p[2];
        }
        assert!(solid::matrices(kind, &bad, &STEEL, 0.0).is_err());
    }
}

fn block_case(element: &str, divisions: [usize; 3], size: [f64; 3], extra: &str) -> Case {
    let text = format!(
        r#"
[analysis]
static = true
[materials.steel]
young = 200e9
poisson = 0.3
density = 7800.0
expansion = 1.2e-5
[[solids]]
name = "b"
material = "steel"
block = {{ origin = [0.0, 0.0, 0.0], size = {size:?}, divisions = {divisions:?}, element = "{element}" }}
{extra}
"#
    );
    Case::from_toml(&text).unwrap()
}

const SYMMETRY: &str = r#"
[[supports]]
part = "b"
surface = "xmin"
fix = ["ux"]
[[supports]]
part = "b"
surface = "ymin"
fix = ["uy"]
[[supports]]
part = "b"
surface = "zmin"
fix = ["uz"]
"#;

#[test]
fn patch_test_on_distorted_meshes() {
    // tension σ along x on a block whose inner nodes are moved: the
    // displacement is linear and exact
    let sigma = 50e6;
    for element in ["hex20", "tet10"] {
        let extra = format!(
            "{SYMMETRY}[[pressures]]\npart = \"b\"\nsurface = \"xmax\"\nvalue = {}\n",
            -sigma
        );
        let c = block_case(element, [3, 2, 2], [0.03, 0.02, 0.02], &extra);
        let mut built = case::build(&c).unwrap();
        let m = &mut built.model;
        let on_side = |p: [f64; 3]| {
            let s = [0.03, 0.02, 0.02];
            (0..3).any(|d| p[d].abs() < 1e-12 || (p[d] - s[d]).abs() < 1e-12)
        };
        for (i, p) in m.nodes.iter_mut().enumerate() {
            if !on_side(*p) {
                let t = i as f64;
                p[0] += 1.5e-3 * (0.7 * t).sin();
                p[1] += 1.0e-3 * (1.3 * t).cos();
                p[2] += 1.0e-3 * (0.4 * t).sin();
            }
        }
        m.check_solids().unwrap();
        let s = analysis::solve_static(m).unwrap();
        let e = sigma / 200e9;
        for (n, p) in m.nodes.iter().enumerate() {
            let want = [e * p[0], -0.3 * e * p[1], -0.3 * e * p[2]];
            for d in 0..3 {
                assert!(
                    (s.u[DOFS * n + d] - want[d]).abs() < 1e-9 * e * 0.03,
                    "{element} node {n} dof {d}: {} {}",
                    s.u[DOFS * n + d],
                    want[d]
                );
            }
        }
    }
}

#[test]
fn cantilever_against_beam_theory() {
    // 100 × 10 × 10 mm, clamped at x = 0, a shear force on the tip face;
    // ν = 0 so the clamp adds no Poisson restraint
    let (l, h, force): (f64, f64, f64) = (0.1, 0.01, 100.0);
    let e = 200e9;
    let i = h.powi(4) / 12.0;
    let g = e / 2.0;
    let kappa = 5.0 / 6.0;
    let timoshenko = force * l.powi(3) / (3.0 * e * i) + force * l / (kappa * g * h * h);
    // the Timoshenko beam (shear and rotary inertia) of the same bar, 40
    // elements, for the frequencies
    let beam = Case::from_toml(&format!(
        r#"
[analysis]
modes = 1
[materials.steel]
young = {e}
poisson = 0.0
density = 7800.0
[sections.sq]
shape = "rectangle"
width = {h}
thickness = {h}
kappa_y = {kappa}
kappa_z = {kappa}
[[beams]]
name = "bar"
from = [0.0, 0.0, 0.0]
to = [{l}, 0.0, 0.0]
elements = 40
section = "sq"
material = "steel"
[[supports]]
at = [0.0, 0.0, 0.0]
fix = ["all"]
"#
    ))
    .unwrap();
    let f_beam = case::run(&beam).unwrap().modal.unwrap().modes[0].frequency();
    for (element, div, tol) in [("hex20", [10, 1, 1], 3e-3), ("tet10", [20, 2, 2], 3e-3)] {
        let extra = format!(
            r#"
[[supports]]
part = "b"
surface = "xmin"
fix = ["all"]
[[pressures]]
part = "b"
surface = "xmax"
force_per_area = [0.0, 0.0, {}]
[points]
tip = [{l}, {}, 0.0]
"#,
            force / (h * h),
            h / 2.0,
        );
        let mut c = block_case(element, div, [l, h, h], &extra);
        c.materials.get_mut("steel").unwrap().poisson = 0.0;
        c.analysis.modes = 2;
        let out = case::run(&c).unwrap();
        let s = out.statics.as_ref().unwrap();
        let n = out.built.points["tip"];
        let w = s.u[DOFS * n + 2];
        eprintln!(
            "{element}: tip {w:.6e}, Timoshenko {timoshenko:.6e} ({:+.3} %)",
            100.0 * (w / timoshenko - 1.0)
        );
        assert!(rel(w, timoshenko) < tol, "{element}: {w} {timoshenko}");
        // the reaction carries the force
        let rf = out.to_json("c")["static"]["reaction_force"][2]
            .as_f64()
            .unwrap();
        assert!(rel(-rf, force) < 1e-9, "{rf}");
        // the first two modes (bending in y and z, equal by symmetry)
        // below Euler–Bernoulli's 1.875104² / (2π L²) √(EI/ρA) by the
        // shear and rotary inertia of the Timoshenko beam
        let md = out.modal.as_ref().unwrap();
        let eb = 1.875_104_068_711_961f64.powi(2) / (2.0 * std::f64::consts::PI * l * l)
            * (e * i / (7800.0 * h * h)).sqrt();
        let f1 = md.modes[0].frequency();
        let f2 = md.modes[1].frequency();
        eprintln!(
            "{element}: f1 {f1:.3} Hz, f2 {f2:.3} Hz, Timoshenko beam {f_beam:.3} Hz, Euler–Bernoulli {eb:.3} Hz"
        );
        assert!(rel(f1, f2) < 1e-3, "{f1} {f2}");
        assert!(rel(f1, f_beam) < 3e-3, "{f1} {f_beam}");
        assert!(f1 < eb, "{f1} {eb}");
    }
}

/// The quarter of a thick-walled cylinder (radii a and b, plane strain)
/// from a block in (r, θ, z): the nodes are mapped after the build.
fn cylinder(element: &str, nu: f64, div: [usize; 3]) -> (case::Built, f64, f64, f64) {
    let (a, b, p) = (0.01, 0.02, 100e6);
    let extra = format!(
        r#"
[[supports]]
part = "b"
surface = "ymin"
fix = ["uy"]
[[supports]]
part = "b"
surface = "ymax"
fix = ["ux"]
[[supports]]
part = "b"
surface = "zmin"
fix = ["uz"]
[[supports]]
part = "b"
surface = "zmax"
fix = ["uz"]
[[pressures]]
part = "b"
surface = "xmin"
value = {p}
"#
    );
    let mut c = block_case(
        element,
        div,
        [b - a, std::f64::consts::FRAC_PI_2, 0.002],
        &extra,
    );
    c.solids[0].block.as_mut().unwrap().origin = [a, 0.0, 0.0];
    c.materials.get_mut("steel").unwrap().poisson = nu;
    let mut built = case::build(&c).unwrap();
    for x in &mut built.model.nodes {
        let (r, t) = (x[0], x[1]);
        x[0] = r * t.cos();
        x[1] = r * t.sin();
    }
    (built, a, b, p)
}

#[test]
fn thick_cylinder_against_lame() {
    for (element, nu, div, tol) in [
        ("hex20", 0.3, [4, 8, 1], 2e-4),
        ("tet10", 0.3, [4, 8, 1], 3e-3),
        ("hex20", 0.49, [4, 8, 1], 2e-3),
        ("tet10", 0.49, [4, 8, 1], 1e-2),
    ] {
        let (built, a, b, p) = cylinder(element, nu, div);
        let m = &built.model;
        m.check_solids().unwrap();
        let s = analysis::solve_static(m).unwrap();
        let e = 200e9;
        let lame = |r: f64| {
            (1.0 + nu) * p * a * a / (e * (b * b - a * a)) * ((1.0 - 2.0 * nu) * r + b * b / r)
        };
        let mut worst = 0.0f64;
        for (n, x) in m.nodes.iter().enumerate() {
            let r = (x[0] * x[0] + x[1] * x[1]).sqrt();
            let ur = (s.u[DOFS * n] * x[0] + s.u[DOFS * n + 1] * x[1]) / r;
            worst = worst.max(rel(ur, lame(r)));
        }
        eprintln!("{element} ν = {nu}: largest error of u_r {:.3e}", worst);
        assert!(worst < tol, "{element} ν {nu}: {worst}");
    }
}

#[test]
fn thermal_expansion_free_and_held() {
    for element in ["hex20", "tet10"] {
        // free: u = α ΔT x and no stress
        let extra = format!("{SYMMETRY}[[temperatures]]\npart = \"b\"\nchange = 40.0\n");
        let c = block_case(element, [2, 2, 1], [0.02, 0.01, 0.005], &extra);
        let out = case::run(&c).unwrap();
        let s = out.statics.as_ref().unwrap();
        let m = &out.built.model;
        for (n, p) in m.nodes.iter().enumerate() {
            for d in 0..3 {
                let want = 1.2e-5 * 40.0 * p[d];
                assert!((s.u[DOFS * n + d] - want).abs() < 1e-15, "{element}");
            }
        }
        let vm = out.to_json("t")["static"]["max_von_mises"]["value"]
            .as_f64()
            .unwrap();
        assert!(vm < 1e-3, "{element}: {vm}");
        // held on every side: σ = −E α ΔT / (1 − 2ν) on each axis
        let held = format!(
            "{SYMMETRY}[[supports]]\npart = \"b\"\nsurface = \"xmax\"\nfix = [\"ux\"]\n\
             [[supports]]\npart = \"b\"\nsurface = \"ymax\"\nfix = [\"uy\"]\n\
             [[supports]]\npart = \"b\"\nsurface = \"zmax\"\nfix = [\"uz\"]\n\
             [[temperatures]]\npart = \"b\"\nchange = 40.0\n"
        );
        let c = block_case(element, [2, 2, 1], [0.02, 0.01, 0.005], &held);
        let out = case::run(&c).unwrap();
        let s = out.statics.as_ref().unwrap();
        let m = &out.built.model;
        let sigma = 200e9 * 1.2e-5 * 40.0 / (1.0 - 0.6);
        let fx: f64 = (0..m.nodes.len())
            .filter(|&n| (m.nodes[n][0] - 0.02).abs() < 1e-12)
            .map(|n| s.reactions[DOFS * n])
            .sum();
        // the support at x = 0.02 pushes back along −x
        assert!(rel(fx, -sigma * 0.01 * 0.005) < 1e-10, "{element}: {fx}");
        assert!(s.u.iter().all(|u| u.abs() < 1e-15));
    }
}

/// A mesh written as a Gmsh `.msh` 4.1 ASCII file of second order: the
/// elements in Gmsh's node order, the block's sides as physical surfaces.
fn to_msh(m: &foehncraft_structure::solid_mesh::SolidMesh) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    let names: Vec<&String> = m.surfaces.keys().collect();
    let _ = writeln!(s, "$MeshFormat\n4.1 0 8\n$EndMeshFormat");
    let _ = writeln!(s, "$PhysicalNames\n{}", names.len() + 1);
    for (i, n) in names.iter().enumerate() {
        let _ = writeln!(s, "2 {} \"{n}\"", i + 1);
    }
    let _ = writeln!(s, "3 100 \"body\"\n$EndPhysicalNames");
    let _ = writeln!(s, "$Entities\n0 0 {} 1", names.len());
    for i in 0..names.len() {
        let _ = writeln!(s, "{} 0 0 0 1 1 1 1 {} 0", i + 1, i + 1);
    }
    let _ = writeln!(s, "1 0 0 0 1 1 1 1 100 0\n$EndEntities");
    let nn = m.nodes.len();
    let _ = writeln!(s, "$Nodes\n1 {nn} 1 {nn}\n3 1 0 {nn}");
    for i in 0..nn {
        let _ = writeln!(s, "{}", i + 1);
    }
    for p in &m.nodes {
        let _ = writeln!(s, "{} {} {}", p[0], p[1], p[2]);
    }
    let _ = writeln!(s, "$EndNodes");
    let faces: usize = m.surfaces.values().map(|v| v.len()).sum();
    let total = faces + m.elements.len();
    let _ = writeln!(s, "$Elements\n{} {total} 1 {total}", names.len() + 1);
    let mut tag = 1;
    for (i, n) in names.iter().enumerate() {
        let f = &m.surfaces[*n];
        let t = if f[0].len() == 8 { 16 } else { 9 };
        let _ = writeln!(s, "2 {} {t} {}", i + 1, f.len());
        for face in f {
            // reversed: the reader must orient it from the element
            let r = solid::flip_face(face);
            let ids: Vec<String> = r.iter().map(|x| (x + 1).to_string()).collect();
            let _ = writeln!(s, "{tag} {}", ids.join(" "));
            tag += 1;
        }
    }
    let (kind, _) = &m.elements[0];
    let (t, perm): (u8, &[usize]) = match kind {
        Kind::Tet10 => (11, &solid::GMSH_TET10),
        Kind::Hex20 => (17, &solid::GMSH_HEX20),
    };
    let _ = writeln!(s, "3 1 {t} {}", m.elements.len());
    for (_, e) in &m.elements {
        // Gmsh node g holds our node i with perm[i] = g
        let mut g = vec![0; e.len()];
        for (i, &p) in perm.iter().enumerate() {
            g[p] = e[i];
        }
        let ids: Vec<String> = g.iter().map(|x| (x + 1).to_string()).collect();
        let _ = writeln!(s, "{tag} {}", ids.join(" "));
        tag += 1;
    }
    let _ = writeln!(s, "$EndElements");
    s
}

#[test]
fn gmsh_file_gives_the_same_answer_as_the_block() {
    let dir = std::env::temp_dir().join(format!("fc-solids-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for (element, kind) in [("hex20", Kind::Hex20), ("tet10", Kind::Tet10)] {
        let size = [0.03, 0.01, 0.004];
        let div = [6, 2, 1];
        let mesh = foehncraft_structure::solid_mesh::block([0.0; 3], size, div, kind);
        let file = dir.join(format!("{element}.msh"));
        std::fs::write(&file, to_msh(&mesh)).unwrap();
        let loads = r#"
[[supports]]
part = "b"
surface = "xmin"
fix = ["all"]
[[pressures]]
part = "b"
surface = "zmax"
value = 2e5
[points]
tip = [0.03, 0.01, 0.004]
"#;
        let block = block_case(element, div, size, loads);
        let mut from_file = block.clone();
        from_file.solids[0].block = None;
        from_file.solids[0].mesh = Some(format!("{element}.msh"));
        from_file.solids[0].volumes = vec!["body".into()];
        from_file.base = Some(dir.clone());
        let a = case::run(&block).unwrap();
        let b = case::run(&from_file).unwrap();
        let ua = a.statics.as_ref().unwrap().u[DOFS * a.built.points["tip"] + 2];
        let ub = b.statics.as_ref().unwrap().u[DOFS * b.built.points["tip"] + 2];
        eprintln!("{element}: block {ua:.9e}, Gmsh {ub:.9e}");
        assert!(ua < 0.0, "the pressure pushes down: {ua}");
        assert!(rel(ub, ua) < 1e-12, "{element}: {ua} {ub}");
        let j = b.to_json("g");
        assert_eq!(j["solids"][0]["mesh"], "MSH 4.1 ASCII");
        assert_eq!(j["solid_elements"], a.built.model.solids.len());
        // a volume that does not exist, a surface that does not exist
        let mut c = from_file.clone();
        c.solids[0].volumes = vec!["nope".into()];
        let e = case::run(&c).err().unwrap();
        assert_eq!(e.code, codes::CASE);
        assert!(e.message.contains("body"), "{e}");
        let mut c = from_file.clone();
        c.pressures[0].surface = Some("top".into());
        let e = case::run(&c).err().unwrap();
        assert!(e.message.contains("zmax"), "{e}");
    }
    // a linear mesh is refused with the fix
    let lin = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 4 1 4\n3 1 0 4\n1\n2\n3\n4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n$Elements\n1 1 1 1\n3 1 4 1\n1 1 2 3 4\n$EndElements\n";
    std::fs::write(dir.join("lin.msh"), lin).unwrap();
    let mut c = block_case("hex20", [1, 1, 1], [1.0, 1.0, 1.0], "");
    c.solids[0].block = None;
    c.solids[0].mesh = Some("lin.msh".into());
    c.base = Some(dir.clone());
    let e = case::run(&c).err().unwrap();
    assert!(e.message.contains("second order"), "{e}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn solid_case_errors() {
    for (extra, want) in [
        (
            "[[temperatures]]\npart = \"b\"\ndifference = 3.0\n",
            "uniform",
        ),
        ("[[supports]]\npart = \"b\"\nfix = [\"all\"]\n", "together"),
        (
            "[[pressures]]\npart = \"b\"\nsurface = \"xmin\"\nplate = \"p\"\n",
            "plate",
        ),
    ] {
        let c = block_case("hex20", [1, 1, 1], [1.0, 1.0, 1.0], extra);
        let e = case::run(&c).err().unwrap();
        assert_eq!(e.code, codes::CASE);
        assert!(e.message.contains(want), "{extra}: {e}");
    }
    let mut c = block_case("hex20", [1, 0, 1], [1.0, 1.0, 1.0], "");
    assert!(case::run(&c).err().unwrap().message.contains("divisions"));
    c.solids[0].block.as_mut().unwrap().divisions = [1, 1, 1];
    c.solids[0].mesh = Some("x.msh".into());
    assert!(case::run(&c).err().unwrap().message.contains("one of them"));
}
