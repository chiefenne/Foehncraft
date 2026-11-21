//! The cell-quality remediation (D50): a pipe in raw Gmsh tetrahedra
//! (`examples/meshes/pipe_gmsh_raw.geo`, faces with a skewness up to 1.23)
//! diverges with the least-squares pressure gradient; with
//! `solver.remediation = "auto"` the solver switches to the Gauss pressure
//! gradient, the Venkatakrishnan limiter and the bounded face velocity of
//! the mass flux, reports it and converges.
// D202: a named unsafe module (tools/check/unsafe_modules.py);
// every unsafe block has a `// SAFETY:` comment.
#![allow(unsafe_code)]

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("foehncraft_quality_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

/// Developing laminar pipe flow at Re = 100 on the raw Gmsh mesh.
fn pipe(remediation: &str) -> FlowCase {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/meshes");
    FlowCase::from_toml(
        &format!(
            r#"
[mesh]
file = "pipe_gmsh_raw.msh"
renumber = "rcm"
[fluid]
viscosity = 1e-3
[boundary.inlet_main]
type = "velocity-inlet"
normal_velocity = 1.0
[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.wall]
type = "wall"
[solver]
tolerance = 1e-6
max_iterations = 400
remediation = "{remediation}"
[output]
vtkhdf = false
wall_csv = false
"#
        ),
        &dir,
    )
    .unwrap()
}

#[test]
fn raw_gmsh_tetrahedra_converge_with_the_remediation() {
    let d = out("auto");
    let r = run(&pipe("auto"), "pipe", &d, |_| {}).unwrap();
    let s = &r.solves[0];
    assert!(s.converged, "{} iterations", s.iterations);
    let q = s.remediation.as_ref().expect("a remediation report");
    assert!(q.skewed_faces > 0 && q.poor_cells > 0, "{q:?}");
    assert!(q.max_skewness > 0.9, "{q:?}");
    assert_eq!(
        q.changes,
        [
            "solver.pressure_gradient: \"least-squares\" -> \"gauss\"",
            "solver.momentum.limiter: \"none\" -> \"venkatakrishnan\"",
            "solver.face_velocity: \"interpolated\" -> \"bounded\""
        ]
    );
    let codes: Vec<&str> = r.warnings.iter().map(|w| w.code).collect();
    assert!(codes.contains(&diag::MESH_SKEWED), "{codes:?}");
    assert!(codes.contains(&diag::MESH_REMEDIATION), "{codes:?}");
    // the report says what was changed and how to avoid it
    let md = std::fs::read_to_string(d.join("flow.md")).unwrap();
    assert!(md.contains("Cell quality (D50)") && md.contains("Mesh.OptimizeNetgen"));
    // the mass balance of the converged state
    assert!(s.mass_imbalance < 1e-6, "{}", s.mass_imbalance);
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn raw_gmsh_tetrahedra_fail_without_the_remediation() {
    let d = out("off");
    let e = match run(&pipe("off"), "pipe", &d, |_| {}) {
        Ok(r) => panic!(
            "the least-squares pressure gradient ran through (converged: {})",
            r.solves.iter().all(|s| s.converged)
        ),
        Err(e) => e,
    };
    // a divergence (the pressure-correction matrix stops being positive
    // definite before a value is non-finite), not a missing pressure level
    assert_eq!(diag::of_error(&e).code, diag::RUN_DIVERGED, "{e:#}");
    // the report names the skewed faces first and the fix
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(d.join("divergence.json")).unwrap()).unwrap();
    let first = &j["causes"][0];
    assert_eq!(first["code"], "DIV003", "{j:#}");
    assert!(
        first["fix"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f == "solver.remediation=\"auto\""),
        "{first:#}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_remediation_reports_poor_faces_but_changes_nothing_when_off() {
    let d = out("off_report");
    let mut c = pipe("off");
    c.solver.max_iterations = 1;
    let r = run(&c, "pipe", &d, |_| {}).unwrap();
    let q = r.solves[0].remediation.as_ref().unwrap();
    assert!(q.skewed_faces > 0 && q.changes.is_empty(), "{q:?}");
    let codes: Vec<&str> = r.warnings.iter().map(|w| w.code).collect();
    assert!(codes.contains(&diag::MESH_SKEWED) && !codes.contains(&diag::MESH_REMEDIATION));
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn the_gauss_weights_give_the_exact_gradient_of_a_linear_field_without_skewness() {
    use foehncraft_fv::disc::Discretisation;
    use foehncraft_fv::problem::{Boundary, Problem};
    use foehncraft_mesh::Geometry;
    use foehncraft_mesh::generate::{BoxSpec, Grading, Kind};
    let mut spec = BoxSpec {
        kind: Kind::Hex,
        upper: [2.0, 1.0, 0.5],
        cells: [7, 5, 3],
        ..Default::default()
    };
    spec.grading.x = Grading {
        ratio: 3.0,
        symmetric: false,
    };
    spec.grading.y = Grading {
        ratio: 4.0,
        symmetric: true,
    };
    let m = spec.generate().unwrap();
    let g = Geometry::compute(&m);
    let mut prob = Problem::default();
    for p in &m.boundaries {
        prob.boundary.insert(p.name.clone(), Boundary::value(0.0));
    }
    let d = Discretisation::new(&m, &g, &prob).unwrap();
    let w = foehncraft_flow::simplec::gauss_weights(&d);
    let a = [0.7, -1.3, 2.1];
    let phi = |x: [f64; 3]| 0.4 + a[0] * x[0] + a[1] * x[1] + a[2] * x[2];
    let ni = d.n_interior;
    for c in 0..d.n_cells {
        let pc = phi(g.cell_centroid[c]);
        let mut grad = [0.0; 3];
        for k in d.cell_ptr[c] as usize..d.cell_ptr[c + 1] as usize {
            let (f, is_nb) = foehncraft_kernels::unpack_face(d.cell_faces[k]);
            let other = if f >= ni {
                phi(g.face_centroid[f])
            } else if is_nb {
                phi(g.cell_centroid[d.owner[f] as usize])
            } else {
                phi(g.cell_centroid[d.neighbour[f] as usize])
            };
            for (i, gi) in grad.iter_mut().enumerate() {
                *gi += [w.x[k], w.y[k], w.z[k]][i] * (other - pc);
            }
        }
        for i in 0..3 {
            assert!((grad[i] - a[i]).abs() < 1e-12, "cell {c}: {grad:?}");
        }
    }
    // and no face of this mesh is skewed
    let skew = foehncraft_flow::simplec::face_skewness(&d);
    assert!(skew.iter().all(|&s| s < 1e-12));
}

#[test]
fn the_remediation_changes_only_what_the_case_does_not_set() {
    use foehncraft_flow::simplec::{
        FaceVelocity, PressureGradient, Remediation, Settings, TransposedStress,
    };
    use foehncraft_fv::solver::{LimiterKind, Scheme};
    let s = Settings::default();
    // a mesh without poor faces, or the remediation off: nothing changes
    assert!(s.remedied(0.9).1.is_empty());
    let off = Settings {
        remediation: Remediation::Off,
        ..Settings::default()
    };
    assert!(off.remedied(1.3).1.is_empty());
    // poor faces: the Gauss pressure gradient, the limiter and the bounded
    // face velocity of the mass flux
    let (r, ch) = s.remedied(1.3);
    assert_eq!(ch.len(), 3);
    assert_eq!(r.pressure_gradient, PressureGradient::Gauss);
    assert_eq!(r.momentum.limiter, LimiterKind::Venkatakrishnan);
    assert_eq!(r.face_velocity, FaceVelocity::Bounded);
    // the transposed stress stays as set
    assert_eq!(r.transposed_stress, TransposedStress::Full);
    // first-order upwind needs no limiter; a limiter of the case stays
    let mut up = Settings::default();
    up.momentum.scheme = Scheme::Upwind;
    assert_eq!(up.remedied(1.3).1.len(), 2);
    let mut wang = Settings::default();
    wang.momentum.limiter = LimiterKind::VenkatakrishnanWang;
    let (r, ch) = wang.remedied(1.3);
    assert_eq!(ch.len(), 2);
    assert_eq!(r.momentum.limiter, LimiterKind::VenkatakrishnanWang);
    // a stricter limit treats milder faces
    let strict = Settings {
        skewness_limit: 0.5,
        ..Settings::default()
    };
    assert_eq!(strict.remedied(0.6).1.len(), 3);
}

/// The discretisation of the raw Gmsh pipe (faces with a skewness up to
/// 1.23).
fn raw_pipe() -> foehncraft_fv::Discretisation {
    use foehncraft_fv::problem::{Boundary, Problem};
    use foehncraft_mesh::Geometry;
    use foehncraft_mesh::case::{MeshSection, build_mesh};
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/meshes/pipe_gmsh_raw.msh");
    let sec: MeshSection = toml::from_str(&format!("file = {:?}", path.to_str().unwrap())).unwrap();
    let (mesh, _, _) = build_mesh(&sec).unwrap();
    let geo = Geometry::compute(&mesh);
    let mut prob = Problem::default();
    for p in &mesh.boundaries {
        prob.boundary
            .insert(p.name.clone(), Boundary::zero_gradient());
    }
    foehncraft_fv::Discretisation::new(&mesh, &geo, &prob).unwrap()
}

fn v3<'a>(x: &'a [f64], y: &'a [f64], z: &'a [f64]) -> foehncraft_kernels::fv::Vec3<'a, f64> {
    foehncraft_kernels::fv::Vec3 { x, y, z }
}

/// The face flux of `face_mass_flux` without the pressure terms (D = 0)
/// for cell velocities `u`, cell gradients `g` (∇u_x, ∇u_y, ∇u_z per
/// component, each three arrays) and the limiter `psi`.
fn plain_fluxes(
    d: &foehncraft_fv::Discretisation,
    u: &[Vec<f64>; 3],
    g: &[Vec<f64>; 9],
    psi: &[Vec<f64>; 3],
    bounded: bool,
) -> Vec<f64> {
    use foehncraft_kernels::flow::{
        self as k, FLUX_FIXED, FaceGeometry, FluxInputs, MOM_FIXED, VelocityGradient,
    };
    let (nc, nf, ni) = (d.n_cells, d.n_faces, d.n_interior);
    let zc = vec![0.0; nc];
    let zf = vec![0.0; nf];
    let zb = vec![0.0; nf - ni];
    let rho_face = vec![1.0; nf];
    let flux_kind = vec![FLUX_FIXED; nf - ni];
    let mom_kind = vec![MOM_FIXED; nf - ni];
    let inp = FluxInputs {
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
        rho_face: &rho_face,
        u: v3(&u[0], &u[1], &u[2]),
        grad: VelocityGradient {
            x: v3(&g[0], &g[1], &g[2]),
            y: v3(&g[3], &g[4], &g[5]),
            z: v3(&g[6], &g[7], &g[8]),
        },
        p: &zc,
        grad_p: v3(&zc, &zc, &zc),
        pb: &zb,
        pset: &zb,
        ub: v3(&zb, &zb, &zb),
        flux_kind: &flux_kind,
        mom_kind: &mom_kind,
        d: &zc,
        relax: &zf,
        transient: false,
        ds: &[],
        time: &[],
        alpha: 1.0,
        psi: v3(&psi[0], &psi[1], &psi[2]),
        bounded,
    };
    // SAFETY: every face of the mesh, every array has its entries.
    (0..ni)
        .map(|f| unsafe { k::face_mass_flux(f, &inp) })
        .collect()
}

#[test]
fn the_bounded_face_velocity_is_exact_for_a_linear_field_and_bounded_on_skewed_faces() {
    let d = raw_pipe();
    let (nc, ni) = (d.n_cells, d.n_interior);
    let skew = foehncraft_flow::simplec::face_skewness(&d);
    assert!(skew.iter().any(|&s| s > 0.9));
    // u = u0 + G x with the exact gradient: both forms give u(x_f)·S
    let u0 = [0.3, -0.2, 1.1];
    let gm = [[0.5, -1.0, 0.2], [0.7, 0.1, -0.4], [-0.3, 0.9, -0.6]];
    let at = |x: [f64; 3], i: usize| u0[i] + gm[i][0] * x[0] + gm[i][1] * x[1] + gm[i][2] * x[2];
    let u: [Vec<f64>; 3] =
        std::array::from_fn(|i| (0..nc).map(|c| at(d.cell_centroid.get(c), i)).collect());
    let g: [Vec<f64>; 9] = std::array::from_fn(|n| vec![gm[n / 3][n % 3]; nc]);
    let ones: [Vec<f64>; 3] = std::array::from_fn(|_| vec![1.0; nc]);
    let a = plain_fluxes(&d, &u, &g, &ones, false);
    let b = plain_fluxes(&d, &u, &g, &ones, true);
    for f in 0..ni {
        let (xf, s) = (d.face_centroid.get(f), d.face_area.get(f));
        let exact: f64 = (0..3).map(|i| at(xf, i) * s[i]).sum();
        let scale = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
        assert!(
            (a[f] - exact).abs() < 1e-12 * scale,
            "face {f}: {} {exact}",
            a[f]
        );
        assert!(
            (b[f] - exact).abs() < 1e-12 * scale,
            "face {f}: {} {exact}",
            b[f]
        );
    }
    // a cell-to-cell oscillation with a gradient of the size the cells see:
    // with the limiter at zero the bounded face velocity is the weighted
    // mean of the two cell values; the interpolated gradient moves the face
    // value of a skewed face out of the range of the two cells
    let osc: [Vec<f64>; 3] = std::array::from_fn(|i| {
        (0..nc)
            .map(|c| if (c + i) % 2 == 0 { 1.0 } else { -1.0 })
            .collect()
    });
    let zeros: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; nc]);
    let big: [Vec<f64>; 9] = std::array::from_fn(|n| {
        (0..nc)
            .map(|c| if (c + n) % 3 == 0 { 40.0 } else { -40.0 })
            .collect()
    });
    let lim = plain_fluxes(&d, &osc, &big, &zeros, true);
    let unl = plain_fluxes(&d, &osc, &big, &ones, false);
    let mut out_of_range = 0;
    for f in 0..ni {
        let (p, n) = (d.owner[f] as usize, d.neighbour[f] as usize);
        let s = d.face_area.get(f);
        let w = d.weight[f];
        let mean: f64 = (0..3)
            .map(|i| (w * osc[i][p] + (1.0 - w) * osc[i][n]) * s[i])
            .sum();
        assert!((lim[f] - mean).abs() < 1e-12, "face {f}");
        let bound: f64 = (0..3).map(|i| s[i].abs()).sum();
        if skew[f] > 0.9 && unl[f].abs() > bound {
            out_of_range += 1;
        }
    }
    assert!(
        out_of_range > 0,
        "the interpolated gradient kept every skewed face in range"
    );
}

#[test]
fn the_trace_free_transposed_stress_differs_only_by_the_divergence() {
    use foehncraft_kernels::flow::{
        self as k, FaceGeometry, MOM_FIXED, MomentumInputs, VelocityGradient,
    };
    let d = raw_pipe();
    let (nc, nf, ni) = (d.n_cells, d.n_faces, d.n_interior);
    // a linear field with divergence 0.9: the extra flux −μ (∇·u) S sums to
    // zero over the faces of a cell, so the momentum balance of every cell
    // is the same with both forms; for ∇·u = 0 every face flux is
    let run = |gm: [[f64; 3]; 3], trace_free: f64| -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let at = |x: [f64; 3], i: usize| gm[i][0] * x[0] + gm[i][1] * x[1] + gm[i][2] * x[2];
        let u: [Vec<f64>; 3] =
            std::array::from_fn(|i| (0..nc).map(|c| at(d.cell_centroid.get(c), i)).collect());
        let ub: [Vec<f64>; 3] =
            std::array::from_fn(|i| (ni..nf).map(|f| at(d.face_centroid.get(f), i)).collect());
        let g: [Vec<f64>; 9] = std::array::from_fn(|n| vec![gm[n / 3][n % 3]; nc]);
        let ones = vec![1.0; nc];
        let (zero, mu) = (vec![0.0; nf], vec![1.3; nf]);
        let mom_kind = vec![MOM_FIXED; nf - ni];
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
            u: v3(&u[0], &u[1], &u[2]),
            uc: v3(&u[0], &u[1], &u[2]),
            grad: VelocityGradient {
                x: v3(&g[0], &g[1], &g[2]),
                y: v3(&g[3], &g[4], &g[5]),
                z: v3(&g[6], &g[7], &g[8]),
            },
            psi: v3(&ones, &ones, &ones),
            ub: v3(&ub[0], &ub[1], &ub[2]),
            mom_kind: &mom_kind,
            blend: 1.0,
            central: 0.0,
            beta: 0.0,
            nonorth: 1.0,
            convection: 0.0,
            transposed: 1.0,
            trace_free,
        };
        let faces: Vec<[f64; 3]> = (0..nf)
            .map(|f| {
                // SAFETY: f is a face of the mesh, every array has its entries.
                let (cp, cn, b) = unsafe { k::momentum_face(f, &inp) };
                let p = d.owner[f] as usize;
                let (up, un) = if f < ni {
                    let n = d.neighbour[f] as usize;
                    ([u[0][p], u[1][p], u[2][p]], [u[0][n], u[1][n], u[2][n]])
                } else {
                    ([u[0][p], u[1][p], u[2][p]], [0.0; 3])
                };
                std::array::from_fn(|i| cp * up[i] + cn * un[i] + b[i])
            })
            .collect();
        let mut cells = vec![[0.0; 3]; nc];
        for (f, fl) in faces.iter().enumerate() {
            let p = d.owner[f] as usize;
            for i in 0..3 {
                cells[p][i] += fl[i];
            }
            if f < ni {
                let n = d.neighbour[f] as usize;
                for i in 0..3 {
                    cells[n][i] -= fl[i];
                }
            }
        }
        (faces, cells)
    };
    let free = [[0.5, -1.0, 0.2], [0.7, 0.1, -0.4], [-0.3, 0.9, -0.6]];
    let (a, _) = run(free, 0.0);
    let (b, _) = run(free, 1.0);
    for f in 0..nf {
        for i in 0..3 {
            assert!((a[f][i] - b[f][i]).abs() < 1e-12, "face {f}");
        }
    }
    let div = [[0.5, -1.0, 0.2], [0.7, 0.1, -0.4], [-0.3, 0.9, 0.3]];
    let (_, ca) = run(div, 0.0);
    let (_, cb) = run(div, 1.0);
    for c in 0..nc {
        for i in 0..3 {
            assert!(
                (ca[c][i] - cb[c][i]).abs() < 1e-12,
                "cell {c}: {:?} {:?}",
                ca[c],
                cb[c]
            );
        }
    }
}
