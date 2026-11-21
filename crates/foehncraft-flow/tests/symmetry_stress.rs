//! The full viscous stress on symmetry planes and wedge planes (stage 6
//! part 6g): with the transposed part of the stress (turbulent flow) the
//! normal stress on a symmetry face is 2μ ∂u_n/∂n. On the planes of an
//! axisymmetric wedge that part is half of the hoop stress −2μ v/r²;
//! without it the cells next to the axis of the impinging jet had v/r too
//! high and Nu too low, more so on finer meshes.
//!
//! A linear velocity field with u_n = 0 and no shear on the planes has a
//! constant stress, so the discrete viscous operator of the momentum rows
//! must vanish in every cell, with the transposed part and without it.
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use foehncraft_fv::Discretisation;
use foehncraft_fv::problem::{Boundary, Problem};
use foehncraft_kernels::flow::{
    self as k, FaceGeometry, MOM_FIXED, MOM_SYMMETRY, MomentumInputs, OperatorInputs,
    VelocityGradient,
};
use foehncraft_kernels::fv::Vec3;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::case::{MeshSection, build_mesh};

fn v3<'a>(x: &'a [f64], y: &'a [f64], z: &'a [f64]) -> Vec3<'a, f64> {
    Vec3 { x, y, z }
}

/// The largest |R_i| of the viscous operator over all cells, relative to
/// the largest face flux, for the linear field u = u0 + G x (G constant),
/// with `symmetry` the boundaries taken as symmetry planes and the others
/// fixed at the exact field.
fn residual(mesh_toml: &str, symmetry: &[&str], g: [[f64; 3]; 3], transposed: f64) -> f64 {
    let sec: MeshSection = toml::from_str(mesh_toml).unwrap();
    let (mesh, _, _) = build_mesh(&sec).unwrap();
    let geo = Geometry::compute(&mesh);
    let mut prob = Problem::default();
    for p in &mesh.boundaries {
        prob.boundary
            .insert(p.name.clone(), Boundary::zero_gradient());
    }
    let d = Discretisation::new(&mesh, &geo, &prob).unwrap();
    let (nc, nf, ni) = (d.n_cells, d.n_faces, d.n_interior);
    let u0 = [0.7, 0.0, 0.0];
    let field = |x: [f64; 3]| -> [f64; 3] {
        let mut u = u0;
        for (i, ui) in u.iter_mut().enumerate() {
            for (j, xj) in x.iter().enumerate() {
                *ui += g[i][j] * xj;
            }
        }
        u
    };
    let mut mom_kind = vec![MOM_FIXED; nf - ni];
    for p in &mesh.boundaries {
        if symmetry.contains(&p.name.as_str()) {
            for f in p.range() {
                mom_kind[f - ni] = MOM_SYMMETRY;
            }
        }
    }
    let cell = |c: usize| field(d.cell_centroid.get(c));
    let (mut ux, mut uy, mut uz) = (vec![0.0; nc], vec![0.0; nc], vec![0.0; nc]);
    for c in 0..nc {
        [ux[c], uy[c], uz[c]] = cell(c);
    }
    let (mut bx, mut by, mut bz) = (vec![0.0; nf - ni], vec![0.0; nf - ni], vec![0.0; nf - ni]);
    for f in ni..nf {
        [bx[f - ni], by[f - ni], bz[f - ni]] = field(d.face_centroid.get(f));
    }
    // the exact gradient in every cell
    let gc: Vec<Vec<f64>> = (0..9).map(|n| vec![g[n / 3][n % 3]; nc]).collect();
    let ones = vec![1.0; nc];
    let (zero, mu) = (vec![0.0; nf], vec![1.3; nf]);
    let inp = MomentumInputs {
        topo: d.topology(),
        geo: FaceGeometry {
            area: v3(&d.face_area.x, &d.face_area.y, &d.face_area.z),
            face_centroid: v3(&d.face_centroid.x, &d.face_centroid.y, &d.face_centroid.z),
            cell_centroid: v3(&d.cell_centroid.x, &d.cell_centroid.y, &d.cell_centroid.z),
            dgeo: &d.dgeo,
            k: v3(&d.k.x, &d.k.y, &d.k.z),
            weight: &d.weight,
            shift: v3(&[], &[], &[]),
        },
        mass_flux: &zero,
        mu: &mu,
        u: v3(&ux, &uy, &uz),
        uc: v3(&ux, &uy, &uz),
        grad: VelocityGradient {
            x: v3(&gc[0], &gc[1], &gc[2]),
            y: v3(&gc[3], &gc[4], &gc[5]),
            z: v3(&gc[6], &gc[7], &gc[8]),
        },
        psi: v3(&ones, &ones, &ones),
        ub: v3(&bx, &by, &bz),
        mom_kind: &mom_kind,
        blend: 1.0,
        central: 0.0,
        beta: 0.0,
        nonorth: 1.0,
        convection: 0.0,
        transposed,
        trace_free: 0.0,
    };
    let (mut cp, mut cn) = (vec![0.0; nf], vec![0.0; nf]);
    let (mut fx, mut fy, mut fz) = (vec![0.0; nf], vec![0.0; nf], vec![0.0; nf]);
    for f in 0..nf {
        // SAFETY: f is a face of the mesh, every array has its entries.
        let (a, b, c) = unsafe { k::momentum_face(f, &inp) };
        (cp[f], cn[f]) = (a, b);
        [fx[f], fy[f], fz[f]] = c;
    }
    let op = OperatorInputs {
        topo: d.topology(),
        cp: &cp,
        cn: &cn,
        bf: v3(&fx, &fy, &fz),
        mom_kind: &mom_kind,
        face_area: v3(&d.face_area.x, &d.face_area.y, &d.face_area.z),
        u: v3(&ux, &uy, &uz),
    };
    // the size of one face flux: μ |G| |S|
    let gn = g.iter().flatten().map(|a| a * a).sum::<f64>().sqrt();
    let smax = (0..nf)
        .map(|f| {
            let s = d.face_area.get(f);
            (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt()
        })
        .fold(0.0, f64::max);
    let mut worst = 0.0f64;
    for c in 0..nc {
        // SAFETY: c is a cell of the mesh.
        let r = unsafe { k::momentum_operator(c, &op) };
        worst = worst.max(r.iter().fold(0.0, |m, a| m.max(a.abs())));
    }
    worst / (1.3 * gn * smax)
}

/// A wedge of 5° about the x axis (perturbed hexahedra and the prisms on
/// the axis): the radial field v = a r, u = −2a x, constant stress
/// 2μ diag(−2a, a, a); the hoop stress needs the transposed part on the
/// wedge planes.
#[test]
fn the_stress_of_a_radial_flow_balances_on_a_wedge() {
    let mesh = r#"
[generate]
type = "hex"
dimension = 2
lower = [0.0, 0.0, 0.0]
upper = [1.0, 0.5, 0.1]
cells = [6, 5, 1]
perturb = 0.15
seed = 3
wedge = 5.0
"#;
    let a = 0.8;
    let g = [[-2.0 * a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]];
    let lap = residual(mesh, &["zmin", "zmax"], g, 0.0);
    let full = residual(mesh, &["zmin", "zmax"], g, 1.0);
    eprintln!("wedge: Laplacian {lap:.2e}, full stress {full:.2e}");
    assert!(lap < 1e-10, "Laplacian form {lap}");
    assert!(full < 1e-10, "full stress {full}");
}

/// A plane stagnation flow onto the symmetry plane y = 0 of a box of
/// perturbed hexahedra: u = a x, v = −a y; the normal stress on the plane
/// is 2μ ∂v/∂y.
#[test]
fn the_stress_of_a_stagnation_flow_balances_on_a_symmetry_plane() {
    let mesh = r#"
[generate]
type = "hex"
dimension = 3
lower = [0.0, 0.0, 0.0]
upper = [1.0, 0.5, 0.4]
cells = [5, 4, 3]
perturb = 0.15
seed = 5
"#;
    let a = 0.6;
    let g = [[a, 0.0, 0.0], [0.0, -a, 0.0], [0.0, 0.0, 0.0]];
    let lap = residual(mesh, &["ymin"], g, 0.0);
    let full = residual(mesh, &["ymin"], g, 1.0);
    eprintln!("symmetry plane: Laplacian {lap:.2e}, full stress {full:.2e}");
    assert!(lap < 1e-10, "Laplacian form {lap}");
    assert!(full < 1e-10, "full stress {full}");
}
