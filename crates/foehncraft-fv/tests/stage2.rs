//! Acceptance checks of stage 2 (docs/requirements.md, stage plan row 2).
//!
//! The order-of-accuracy sweep runs in release builds (CI runs the tests in
//! debug and in release); the longer sweep is `#[ignore]`d:
//! `cargo test --release -p foehncraft-fv -- --ignored` (R7).

use foehncraft_fv::backend::{Backend, V, cpu::CpuBackend};
use foehncraft_fv::solver::{LimiterKind, Settings};
use foehncraft_fv::verify::*;
use foehncraft_fv::{Discretisation, Problem};
use foehncraft_kernels::unpack_face;
use foehncraft_mesh::generate::Kind;
use std::path::Path;

const KINDS: [Kind; 4] = [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed];
const NO_PTX: &str = "";

#[test]
fn gradient_of_a_linear_field_is_exact_on_every_mesh_type() {
    for kind in KINDS {
        for dim in [2u8, 3] {
            for perturb in [0.0, 0.2] {
                let (m, g) = test_mesh(kind, dim, 4, 0, perturb).unwrap();
                for grad_boundaries in [&[][..], &["xmin", "ymax"][..]] {
                    let e = gradient_linear(&m, &g, dim, grad_boundaries).unwrap();
                    assert!(
                        e < 1e-13,
                        "{} {dim}D perturb {perturb} {grad_boundaries:?}: error {e:e}",
                        kind.name()
                    );
                }
            }
        }
    }
}

/// With K = 0 the limited reconstruction φ_P + ψ ∇φ·(x_f − x_P) of every
/// face lies inside the range of the cell and its neighbours, for a rough
/// random field on every mesh type.
#[test]
fn limited_face_values_create_no_new_extrema() {
    for kind in KINDS {
        let (m, g) = test_mesh(kind, 3, 4, 0, 0.2).unwrap();
        let p = mms_problem(3);
        let d = Discretisation::new(&m, &g, &p).unwrap();
        let phi: Vec<f64> = (0..m.n_cells)
            .map(|c| foehncraft_core::rng::uniform([c as u64, 7, 9]))
            .collect();
        let mut be = CpuBackend::new(d.clone(), None).unwrap();
        be.upload(V::Phi, &phi).unwrap();
        be.gradient().unwrap();
        be.limiter(foehncraft_fv::backend::Limiter::Venkatakrishnan { k: 0.0 })
            .unwrap();
        let (gx, gy, gz, psi) = (
            be.download(V::Gx).unwrap(),
            be.download(V::Gy).unwrap(),
            be.download(V::Gz).unwrap(),
            be.download(V::Psi).unwrap(),
        );
        let ni = m.n_interior_faces();
        for c in 0..m.n_cells {
            let mut lo = phi[c];
            let mut hi = phi[c];
            for (f, is_nb) in m.cell_face_iter(c) {
                let v = if f < ni {
                    phi[if is_nb { m.owner[f] } else { m.neighbour[f] } as usize]
                } else {
                    d.bc_a[f - ni] * phi[c] + d.bc_b[f - ni]
                };
                lo = lo.min(v);
                hi = hi.max(v);
            }
            for &e in m.cell(c) {
                let (f, _) = unpack_face(e);
                let r = d.face_centroid.get(f);
                let x = g.cell_centroid[c];
                let fv = phi[c]
                    + psi[c]
                        * (gx[c] * (r[0] - x[0]) + gy[c] * (r[1] - x[1]) + gz[c] * (r[2] - x[2]));
                let tol = 1e-12 * (hi - lo).max(1.0);
                assert!(
                    fv >= lo - tol && fv <= hi + tol,
                    "{} cell {c}: {fv} outside [{lo}, {hi}]",
                    kind.name()
                );
            }
        }
    }
}

/// The converged oblique-step solution stays inside [0, 1] with both
/// limiters at K = 0, and overshoots without a limiter.
#[test]
fn oblique_step_stays_bounded_with_the_limiters() {
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly] {
        let (m, g) = test_mesh(kind, 2, 12, 0, 0.1).unwrap();
        let s0 = Settings {
            max_outer: 400,
            ..Default::default()
        };
        let (mn, mx, _, _) = oblique_step(&m, &g, &s0, "cpu", Path::new(NO_PTX)).unwrap();
        assert!(
            mn < -1e-3 || mx > 1.0 + 1e-3,
            "{}: no overshoot unlimited",
            kind.name()
        );
        for lim in [
            LimiterKind::Venkatakrishnan,
            LimiterKind::VenkatakrishnanWang,
        ] {
            let s = Settings {
                limiter: lim,
                limiter_k: 0.0,
                ..s0.clone()
            };
            let (mn, mx, o, _) = oblique_step(&m, &g, &s, "cpu", Path::new(NO_PTX)).unwrap();
            assert!(o.converged, "{} {lim:?}: not converged", kind.name());
            assert!(
                mn > -1e-7 && mx < 1.0 + 1e-6,
                "{} {lim:?}: range [{mn}, {mx}]",
                kind.name()
            );
        }
    }
}

fn sweep(levels3: &[u32], levels2: &[u32]) {
    let s = Settings::default();
    let mut failures = Vec::new();
    for kind in KINDS {
        for (dim, levels) in [(3u8, levels3), (2, levels2)] {
            let r = order_sweep(kind, dim, 4, levels, 0.15, &s, "cpu", Path::new(NO_PTX)).unwrap();
            let orders = observed_orders(&r, dim);
            let last = *orders.last().unwrap();
            eprintln!(
                "{} {dim}D cells {:?} orders {orders:.2?}",
                kind.name(),
                r.iter().map(|x| x.cells).collect::<Vec<_>>()
            );
            assert!(r.iter().all(|x| x.converged));
            if !(1.8..=2.2).contains(&last) {
                failures.push(format!("{} {dim}D: {last:.3}", kind.name()));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "observed order outside 1.8-2.2: {failures:?}"
    );
}

/// Convection–diffusion with a manufactured solution: observed order 1.8–2.2
/// on hexahedral, tetrahedral, polyhedral and mixed meshes, 2D and 3D,
/// perturbed nodes (acceptance test of stage 2).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn order_of_accuracy_on_all_mesh_types() {
    sweep(&[1, 2, 3], &[2, 3, 4]);
}

/// One level more (about 1.6 million tetrahedra).
#[test]
#[ignore = "slow: cargo test --release -p foehncraft-fv -- --ignored"]
fn order_of_accuracy_finer() {
    sweep(&[2, 3, 4], &[3, 4, 5]);
}

/// Bitwise the same results for 1 to 4 threads and from run to run (D31).
#[test]
fn results_are_bitwise_reproducible() {
    for kind in [Kind::Poly, Kind::Mixed] {
        let (m, g) = test_mesh(kind, 3, 4, 0, 0.15).unwrap();
        let p: Problem = mms_problem(3);
        let s = Settings {
            limiter: LimiterKind::Venkatakrishnan,
            limiter_k: 0.3,
            max_outer: 30,
            ..Default::default()
        };
        let mut reference: Option<Vec<u64>> = None;
        for threads in [1, 2, 3, 4, 4] {
            let (_, phi) = manufactured(
                "x",
                &m,
                &g,
                3,
                &p,
                &s,
                "cpu",
                Some(threads),
                Path::new(NO_PTX),
            )
            .unwrap();
            let bits: Vec<u64> = phi.iter().map(|x| x.to_bits()).collect();
            match &reference {
                None => reference = Some(bits),
                Some(r) => assert!(r == &bits, "{}: {threads} threads differ", kind.name()),
            }
        }
    }
}

/// BiCGStab with Jacobi reduces the residual of a pure diffusion system by
/// the requested factor.
#[test]
fn bicgstab_reaches_its_tolerance() {
    let (m, g) = test_mesh(Kind::Tet, 3, 4, 0, 0.1).unwrap();
    let mut p = mms_problem(3);
    p.velocity = [0.0; 3];
    let d = Discretisation::new(&m, &g, &p).unwrap();
    let mut be = CpuBackend::new(d, None).unwrap();
    be.upload(V::Phi, &vec![0.0; m.n_cells]).unwrap();
    be.gradient().unwrap();
    be.faces(foehncraft_fv::backend::FaceSettings {
        blend: 0.0,
        central: 0.0,
        beta: 0.0,
        nonorth: false,
        nonorth_limit: -1.0,
    })
    .unwrap();
    be.assemble(1.0).unwrap();
    be.residual(V::Phi, V::Rhs, V::R).unwrap();
    let r0 = be.dot(V::R, V::R).unwrap().sqrt();
    let it = foehncraft_fv::solver::bicgstab(
        &mut be,
        foehncraft_fv::solver::Preconditioner::Jacobi,
        500,
        1e-10,
        r0,
    )
    .unwrap();
    be.residual(V::Phi, V::Rhs, V::R).unwrap();
    let r1 = be.dot(V::R, V::R).unwrap().sqrt();
    assert!(
        it < 500 && r1 <= 2e-10 * r0,
        "{it} iterations, {r1:e} / {r0:e}"
    );
}

/// Two cells that share several faces (the pieces of a merged cut cell
/// around a bend, L13 D228) are one face to the two-point flux: over the
/// pieces of each pair Σ |S|²/(d·S) and Σ k are those of the summed area
/// vector, all pieces take the pair's interpolation weight, and the
/// least-squares gradient of a linear field stays exact.
#[test]
fn faces_of_one_cell_pair_act_as_one_face() {
    use foehncraft_mesh::cut::{CutOptions, Sphere, cut};
    use foehncraft_mesh::generate::BoxSpec;
    use foehncraft_mesh::geometry::{add, dot, norm, scale, sub};
    let mut s = BoxSpec::new(Kind::Poly, 3, [16, 16, 16]);
    s.lower = [0.0; 3];
    s.upper = [1.0; 3];
    let m = s.generate().unwrap();
    let body = Sphere {
        centre: [0.5, 0.5, 0.5],
        radius: 0.2,
        exact: true,
    };
    let opts = CutOptions {
        min_fraction: 0.3,
        wall_point: true,
        face_point: true,
    };
    let c = cut(&m, &body, &["body".to_string()], &opts).unwrap();
    assert!(c.report.shared_pairs > 0, "{:?}", c.report);
    let mesh = c.mesh;
    let geo = foehncraft_mesh::Geometry::compute(&mesh);
    let mut p = Problem {
        velocity: [0.0; 3],
        ..Default::default()
    };
    for b in &mesh.boundaries {
        p.boundary.insert(
            b.name.clone(),
            foehncraft_fv::problem::Boundary::zero_gradient(),
        );
    }
    let d = Discretisation::new(&mesh, &geo, &p).unwrap();
    let ni = mesh.n_interior_faces();
    let (mut f, mut pairs) = (0, 0);
    while f < ni {
        let mut e = f + 1;
        while e < ni && mesh.owner[e] == mesh.owner[f] && mesh.neighbour[e] == mesh.neighbour[f] {
            e += 1;
        }
        if e - f > 1 {
            pairs += 1;
            let dd = sub(
                geo.cell_centroid[mesh.neighbour[f] as usize],
                geo.cell_centroid[mesh.owner[f] as usize],
            );
            let st = (f..e).fold([0.0; 3], |a, g| add(a, geo.face_area[g]));
            let want = dot(st, st) / dot(dd, st);
            let got: f64 = (f..e).map(|g| d.dgeo[g]).sum();
            assert!(
                (got / want - 1.0).abs() < 1e-12,
                "faces {f}..{e}: {got} vs {want}"
            );
            let k = (f..e).fold([0.0; 3], |a, g| add(a, d.k.get(g)));
            let kw = sub(st, scale(dd, want));
            assert!(norm(sub(k, kw)) < 1e-12 * norm(st), "faces {f}..{e}");
            for g in f + 1..e {
                assert_eq!(d.weight[g], d.weight[f]);
            }
            assert!((0.0..=1.0).contains(&d.weight[f]));
        }
        f = e;
    }
    assert!(pairs > 0);
    // the least-squares gradient of φ = a · x (zero normal gradient at
    // the boundaries does not hold for it: interior cells only)
    let a = [0.3, -1.1, 0.7];
    let phi: Vec<f64> = (0..mesh.n_cells)
        .map(|c| dot(a, geo.cell_centroid[c]))
        .collect();
    for cell in 0..mesh.n_cells {
        let lo = mesh.cell_ptr[cell] as usize;
        let hi = mesh.cell_ptr[cell + 1] as usize;
        let mut inner = true;
        let mut g = [0.0; 3];
        for (i, &en) in mesh.cell_faces[lo..hi].iter().enumerate() {
            let (fc, is_nb) = unpack_face(en);
            if fc >= ni {
                inner = false;
                break;
            }
            let o = if is_nb {
                mesh.owner[fc]
            } else {
                mesh.neighbour[fc]
            } as usize;
            g = add(g, scale(d.lsq.get(lo + i), phi[o] - phi[cell]));
        }
        if inner {
            assert!(norm(sub(g, a)) < 1e-9, "cell {cell}: {g:?}");
        }
    }
}
