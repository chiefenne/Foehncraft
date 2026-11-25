//! Verification of the radiation call R1a (step 29, section 3): cases 1 to
//! 4 against exact solutions, the exact energy account, and bitwise the
//! same result on any number of threads and ranks.

use foehncraft_mesh::decompose::extract;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_mesh::{Geometry, Mesh};
use foehncraft_parallel::{Communicator, ThreadComm};
use foehncraft_radiation::exact::{
    cube_wall_flux_mean, parallel_rectangles, slab_source, slab_wall_flux,
};
use foehncraft_radiation::{CpuDevicePath, Domain, Points, SIGMA, Settings, Solution, Wall};

fn mesh(kind: Kind, n: [usize; 3], upper: [f64; 3], perturb: f64) -> Mesh {
    let mut s = BoxSpec::new(kind, 3, n);
    s.upper = upper;
    s.perturb = perturb;
    s.generate().unwrap()
}

/// One wall per boundary of the mesh, by name.
fn walls(m: &Mesh, wall: impl Fn(&str) -> Wall) -> Vec<Wall> {
    m.boundaries.iter().map(|b| wall(&b.name)).collect()
}

/// Boundary faces (index from the first boundary face) of the boundary
/// `name`.
fn faces_of(m: &Mesh, name: &str) -> std::ops::Range<usize> {
    let b = m.boundaries.iter().find(|b| b.name == name).unwrap();
    let first = m.n_interior_faces();
    b.start - first..b.start - first + b.len
}

fn area(g: &Geometry, m: &Mesh, f: usize) -> f64 {
    let a = g.face_area[m.n_interior_faces() + f];
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Area-weighted mean of a wall quantity over a boundary, and the error of
/// the mean from the errors of the faces (taken as independent).
fn wall_mean(m: &Mesh, g: &Geometry, name: &str, v: &[f64], e: &[f64]) -> (f64, f64) {
    let (mut s, mut a, mut err) = (0.0, 0.0, 0.0);
    for f in faces_of(m, name) {
        let af = area(g, m, f);
        s += v[f] * af;
        err += (e[f] * af).powi(2);
        a += af;
    }
    (s / a, err.sqrt() / a)
}

fn check_ledger(s: &Solution) {
    assert!(s.ledger.closes(), "{:?}", s.ledger);
    assert_eq!(s.ledger.lost, 0, "{:?}", s.ledger);
    assert!(s.ledger.emitted > 0);
}

/// Case 1: isothermal grey slab between cold black plates (the side walls
/// symmetry planes): wall flux σT⁴ (1 − 2 E₃(τ)) and the source of every
/// layer, for τ from 0.01 to 10, within three standard errors.
#[test]
fn case1_isothermal_slab() {
    let (t, l) = (1000.0f64, 1.0);
    let et4 = SIGMA * t.powi(4);
    for (kappa, n) in [(0.01, 10), (0.1, 10), (1.0, 20), (5.0, 40), (10.0, 50)] {
        let m = mesh(Kind::Hex, [n, 2, 2], [l, 0.5, 0.5], 0.0);
        let g = Geometry::compute(&m);
        let w = walls(&m, |name| match name {
            "xmin" | "xmax" => Wall::black(0.0),
            _ => Wall::mirror(),
        });
        let s = Settings {
            bundles: 200_000,
            ..Settings::default()
        };
        let sol = Domain::whole(&m, &g)
            .solve(&vec![t; m.n_cells], &vec![kappa; m.n_cells], &w, &s)
            .unwrap();
        check_ledger(&sol);
        assert_eq!(sol.ledger.escaped, 0);
        // flux into each plate (net = absorbed, the plates emit nothing)
        let exact = et4 * slab_wall_flux(kappa * l);
        for name in ["xmin", "xmax"] {
            let (q, e) = wall_mean(&m, &g, name, &sol.wall_net, &sol.wall_net_error);
            let dev = (q - exact) / exact;
            eprintln!(
                "slab κ={kappa} {name}: {q:.6e} exact {exact:.6e} ({dev:+.2e}, σ {:.1e})",
                e / exact
            );
            assert!((q - exact).abs() < 3.0 * e, "{name}: {q} against {exact}");
        }
        // source of each layer of cells
        let dx = l / n as f64;
        let mut z2 = 0.0;
        for i in 0..n {
            let (x0, x1) = (i as f64 * dx, (i + 1) as f64 * dx);
            let (mut sv, mut vol, mut err) = (0.0, 0.0, 0.0);
            for c in 0..m.n_cells {
                let x = g.cell_centroid[c][0];
                if x > x0 && x < x1 {
                    let v = g.cell_volume[c];
                    sv += sol.source[c] * v;
                    err += (sol.source_error[c] * v).powi(2);
                    vol += v;
                }
            }
            let (q, e) = (sv / vol, err.sqrt() / vol);
            let exact = et4 * slab_source(kappa, l, x0, x1);
            z2 += ((q - exact) / e).powi(2);
            assert!(
                (q - exact).abs() < 3.0 * e,
                "κ {kappa}, layer {i}: {q} against {exact} (σ {e})"
            );
        }
        eprintln!("slab κ={kappa}: mean z² of the layers {:.2}", z2 / n as f64);
    }
}

/// Case 2: isothermal enclosure (black, then grey walls at the gas
/// temperature): no net source in any cell, no net flux into any wall.
#[test]
fn case2_isothermal_enclosure() {
    let t = 800.0f64;
    let m = mesh(Kind::Tet, [6, 5, 4], [1.2, 1.0, 0.8], 0.2);
    let g = Geometry::compute(&m);
    let kappa: Vec<f64> = g.cell_centroid.iter().map(|c| 0.5 + c[0]).collect();
    for e in [1.0, 0.5] {
        let w = walls(&m, |_| Wall::grey(e, t));
        let sol = Domain::whole(&m, &g)
            .solve(&vec![t; m.n_cells], &kappa, &w, &Settings::default())
            .unwrap();
        check_ledger(&sol);
        let scale = 4.0 * SIGMA * t.powi(4);
        let (mut z2, mut rel) = (0.0, 0.0);
        for c in 0..m.n_cells {
            z2 += (sol.source[c] / sol.source_error[c]).powi(2);
            rel += sol.source[c].abs() / (scale * kappa[c]);
        }
        let n = m.n_cells as f64;
        eprintln!(
            "enclosure ε={e}: mean |source| / emission {:.2e}, mean z² {:.2}",
            rel / n,
            z2 / n
        );
        assert!(rel / n < 0.05, "ε {e}: mean relative source {}", rel / n);
        assert!(
            z2 / n < 3.0,
            "ε {e}: the batch errors are too small: z² {}",
            z2 / n
        );
        // gas and walls in balance as a whole
        let gas: i64 = sol.absorbed.iter().sum::<i64>() - sol.emitted.iter().sum::<i64>();
        let emitted = sol.ledger.emitted as f64;
        assert!(
            (gas as f64).abs() < 0.01 * emitted,
            "ε {e}: gas balance {gas}"
        );
    }
}

/// Case 3: transparent medium between directly opposed black plates (the
/// other boundaries open): the power the cold plate absorbs is the
/// configuration factor times the power of the hot plate.
#[test]
fn case3_parallel_plates() {
    let t = 1000.0;
    for (c, a, b) in [(1.0, 1.0, 1.0), (0.5, 2.0, 1.0)] {
        let m = mesh(Kind::Poly, [3, 4, 4], [c, a, b], 0.2);
        let g = Geometry::compute(&m);
        let w = walls(&m, |name| match name {
            "xmin" => Wall::black(t),
            "xmax" => Wall::black(0.0),
            _ => Wall::open(),
        });
        let s = Settings {
            bundles: 200_000,
            ..Settings::default()
        };
        let sol = Domain::whole(&m, &g)
            .solve(&vec![t; m.n_cells], &vec![0.0; m.n_cells], &w, &s)
            .unwrap();
        check_ledger(&sol);
        assert_eq!(sol.ledger.absorbed_gas, 0);
        let hot: i64 = faces_of(&m, "xmin").map(|f| sol.wall_emitted[f]).sum();
        let cold: i64 = faces_of(&m, "xmax").map(|f| sol.wall_absorbed[f]).sum();
        let back: i64 = faces_of(&m, "xmin").map(|f| sol.wall_absorbed[f]).sum();
        assert_eq!(back, 0, "nothing comes back to the hot plate");
        assert_eq!(hot, sol.ledger.emitted);
        let f = cold as f64 / hot as f64;
        let exact = parallel_rectangles(a, b, c);
        let sigma = (exact * (1.0 - exact) / sol.bundles as f64).sqrt();
        eprintln!("plates {a} x {b} at {c}: F {f:.5} exact {exact:.5} (σ {sigma:.1e})");
        assert!((f - exact).abs() < 4.0 * sigma, "F {f} against {exact}");
    }
}

/// The two coordinates in the wall of a point on a wall of the unit cube,
/// and the wall's axis.
fn in_wall(name: &str, x: [f64; 3]) -> (usize, [f64; 2]) {
    match &name[..1] {
        "x" => (0, [x[1], x[2]]),
        "y" => (1, [x[0], x[2]]),
        _ => (2, [x[0], x[1]]),
    }
}

/// Case 5: unit cube of isothermal grey gas (κL = 1) between cold black
/// walls. The flux into every wall face against the exact flux averaged
/// over the face (`exact::cube_wall_flux`, the gas emission integrated over
/// the five other walls), and the mean flux of every wall, on a hex mesh;
/// the wall means on a perturbed polyhedral mesh.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn case5_cube_enclosure() {
    let (t, kl) = (1000.0f64, 1.0);
    let et4 = SIGMA * t.powi(4);
    let mean = cube_wall_flux_mean(kl, [0.0, 1.0], [0.0, 1.0], 8);
    let names = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"];
    // per face on 10³ hexahedra
    let n = 10;
    let m = mesh(Kind::Hex, [n; 3], [1.0; 3], 0.0);
    let g = Geometry::compute(&m);
    let w = walls(&m, |_| Wall::black(0.0));
    // 32 batches: the error of each face is then close to normal (8 batches
    // give Student's t with 7 degrees of freedom, with wide tails)
    let s = Settings {
        bundles: 1_000_000,
        batches: 32,
        ..Settings::default()
    };
    let sol = Domain::whole(&m, &g)
        .solve(&vec![t; m.n_cells], &vec![kl; m.n_cells], &w, &s)
        .unwrap();
    check_ledger(&sol);
    assert_eq!(sol.ledger.escaped, 0);
    // the exact face means of one wall (every wall is the same)
    let h = 1.0 / n as f64;
    let exact: Vec<f64> = (0..n * n)
        .map(|k| {
            let (i, j) = ((k % n) as f64, (k / n) as f64);
            et4 * cube_wall_flux_mean(kl, [i * h, (i + 1.0) * h], [j * h, (j + 1.0) * h], 1)
        })
        .collect();
    let first = m.n_interior_faces();
    let (mut z2, mut zmax, mut count, mut above3) = (0.0f64, 0.0f64, 0, 0);
    for name in names {
        for f in faces_of(&m, name) {
            let (_, [a, b]) = in_wall(name, g.face_centroid[first + f]);
            let k = (a / h) as usize + n * (b / h) as usize;
            let z = (sol.wall_net[f] - exact[k]) / sol.wall_net_error[f];
            z2 += z * z;
            zmax = zmax.max(z.abs());
            above3 += usize::from(z.abs() > 3.0);
            count += 1;
        }
        let (q, e) = wall_mean(&m, &g, name, &sol.wall_net, &sol.wall_net_error);
        let ex = et4 * mean;
        eprintln!(
            "cube hex {name}: {q:.6e} exact {ex:.6e} ({:+.2e}, σ {:.1e})",
            (q - ex) / ex,
            e / ex
        );
        assert!((q - ex).abs() < 3.0 * e, "{name}: {q} against {ex}");
    }
    let z2 = z2 / count as f64;
    eprintln!(
        "cube hex: {count} wall faces, mean z² {z2:.2}, largest |z| {zmax:.2}, {above3} beyond 3σ"
    );
    assert!(
        z2 < 1.3 && zmax < 4.5 && above3 <= count / 100,
        "mean z² {z2}, largest |z| {zmax}, {above3} beyond 3σ"
    );
    // the mean of every wall on a perturbed polyhedral mesh
    let m = mesh(Kind::Poly, [6, 6, 6], [1.0; 3], 0.25);
    let g = Geometry::compute(&m);
    let w = walls(&m, |_| Wall::black(0.0));
    let s = Settings {
        bundles: 300_000,
        ..Settings::default()
    };
    let sol = Domain::whole(&m, &g)
        .solve(&vec![t; m.n_cells], &vec![kl; m.n_cells], &w, &s)
        .unwrap();
    check_ledger(&sol);
    for name in names {
        let (q, e) = wall_mean(&m, &g, name, &sol.wall_net, &sol.wall_net_error);
        let ex = et4 * mean;
        eprintln!(
            "cube poly {name}: {q:.6e} exact {ex:.6e} ({:+.2e}, σ {:.1e})",
            (q - ex) / ex,
            e / ex
        );
        assert!((q - ex).abs() < 3.0 * e, "{name}: {q} against {ex}");
    }
}

/// Case 6: Sobol against pseudo-random emission. The error of the slab
/// wall flux (case 1, τ = 1) and of the configuration factor of the plates
/// (case 3) against the number of bundles, as the root mean square over
/// independent calls: the scrambled Sobol points converge faster and are
/// several times more accurate at the largest number of bundles.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn case6_sobol_against_random() {
    let t = 1000.0f64;
    let et4 = SIGMA * t.powi(4);
    let slab = mesh(Kind::Hex, [10, 2, 2], [1.0, 0.5, 0.5], 0.0);
    let gs = Geometry::compute(&slab);
    let ws = walls(&slab, |name| match name {
        "xmin" | "xmax" => Wall::black(0.0),
        _ => Wall::mirror(),
    });
    let plates = mesh(Kind::Hex, [2, 4, 4], [1.0; 3], 0.0);
    let gp = Geometry::compute(&plates);
    let wp = walls(&plates, |name| match name {
        "xmin" => Wall::black(t),
        "xmax" => Wall::black(0.0),
        _ => Wall::open(),
    });
    let exact_slab = et4 * slab_wall_flux(1.0);
    let exact_plates = parallel_rectangles(1.0, 1.0, 1.0);
    let calls = 8;
    // root-mean-square relative error over the calls
    let rms = |points: Points, bundles: u64| -> [f64; 2] {
        let (mut e_slab, mut e_plates) = (0.0, 0.0);
        for call in 0..calls {
            let s = Settings {
                bundles,
                call,
                points,
                ..Settings::default()
            };
            let sol = Domain::whole(&slab, &gs)
                .solve(&vec![t; slab.n_cells], &vec![1.0; slab.n_cells], &ws, &s)
                .unwrap();
            let (q, _) = wall_mean(&slab, &gs, "xmin", &sol.wall_net, &sol.wall_net_error);
            e_slab += ((q - exact_slab) / exact_slab).powi(2);
            let sol = Domain::whole(&plates, &gp)
                .solve(
                    &vec![t; plates.n_cells],
                    &vec![0.0; plates.n_cells],
                    &wp,
                    &s,
                )
                .unwrap();
            let hot: i64 = faces_of(&plates, "xmin").map(|f| sol.wall_emitted[f]).sum();
            let cold: i64 = faces_of(&plates, "xmax")
                .map(|f| sol.wall_absorbed[f])
                .sum();
            let f = cold as f64 / hot as f64;
            e_plates += ((f - exact_plates) / exact_plates).powi(2);
        }
        [e_slab, e_plates].map(|e| (e / calls as f64).sqrt())
    };
    let sizes = [4_000u64, 16_000, 64_000, 256_000];
    let mut table = Vec::new();
    for &n in &sizes {
        let (s, r) = (rms(Points::Sobol, n), rms(Points::Random, n));
        eprintln!(
            "N {n:>7}: slab Sobol {:.2e} random {:.2e}; plates Sobol {:.2e} random {:.2e}",
            s[0], r[0], s[1], r[1]
        );
        table.push((s, r));
    }
    // observed order of the error against N, from the first to the last
    let order = |a: f64, b: f64| -(b / a).ln() / (sizes[3] as f64 / sizes[0] as f64).ln();
    for (case, k) in [("slab", 0), ("plates", 1)] {
        let (s0, r0) = table[0];
        let (s3, r3) = table[3];
        let (os, or) = (order(s0[k], s3[k]), order(r0[k], r3[k]));
        eprintln!("{case}: order Sobol {os:.2}, random {or:.2}");
        assert!(
            s3[k] < 0.3 * r3[k],
            "{case}: Sobol {} against random {}",
            s3[k],
            r3[k]
        );
        assert!(
            os > or + 0.05,
            "{case}: order Sobol {os} against random {or}"
        );
    }
}

/// A case with everything: a turned polyhedral mesh, varying temperature
/// and absorption, grey walls at different temperatures, a symmetry plane
/// and an open boundary.
fn mixed_case() -> (Mesh, Vec<f64>, Vec<f64>, Vec<Wall>) {
    let mut m = mesh(Kind::Poly, [8, 6, 5], [2.0, 1.5, 1.0], 0.25);
    let (s, c) = 0.7f64.sin_cos();
    for p in &mut m.points {
        let [x, y, z] = *p;
        *p = [c * x - s * y, s * x + c * y, z];
    }
    let g = Geometry::compute(&m);
    let t: Vec<f64> = g
        .cell_centroid
        .iter()
        .map(|x| 600.0 + 300.0 * (x[0] + 0.5 * x[1]).sin().abs())
        .collect();
    let kappa: Vec<f64> = g.cell_centroid.iter().map(|x| 0.2 + 1.5 * x[2]).collect();
    let w = walls(&m, |name| match name {
        "xmin" => Wall::grey(0.8, 1200.0),
        "xmax" => Wall::grey(0.4, 500.0),
        "ymin" => Wall::black(300.0),
        "ymax" => Wall::mirror(),
        "zmin" => Wall::grey(0.6, 900.0),
        _ => Wall::open(),
    });
    (m, t, kappa, w)
}

/// Case 4a: the energy account closes to the quantum, and the call gives
/// bitwise the same result on 1, 2 and 4 threads.
#[test]
fn case4_ledger_and_threads() {
    let (m, t, kappa, w) = mixed_case();
    let g = Geometry::compute(&m);
    let s = Settings {
        bundles: 40_000,
        ..Settings::default()
    };
    let run = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| Domain::whole(&m, &g).solve(&t, &kappa, &w, &s).unwrap())
    };
    let one = run(1);
    check_ledger(&one);
    assert!(one.ledger.escaped > 0 && one.ledger.absorbed_walls > 0 && one.ledger.absorbed_gas > 0);
    assert_ne!(one.ledger.roulette, 0, "roulette played");
    for threads in [2, 4] {
        assert_eq!(run(threads), one, "{threads} threads");
    }
    // the data path of the GPU (flat records, the body of the PTX entry,
    // tallies per launch) on the CPU gives bitwise the same, also when the
    // bundles go in several launches
    let device = Domain::whole(&m, &g)
        .on_device(Box::new(CpuDevicePath))
        .unwrap();
    assert_eq!(device.device_name(), "cpu-device-path");
    let launches = Settings { launch: 5000, ..s };
    assert_eq!(device.solve(&t, &kappa, &w, &launches).unwrap(), one);
    // another call number gives other bundles
    let other = Domain::whole(&m, &g)
        .solve(&t, &kappa, &w, &Settings { call: 1, ..s })
        .unwrap();
    check_ledger(&other);
    assert_ne!(other.absorbed, one.absorbed);
}

/// Case 4b: on 2, 3 and 4 ranks (threads as ranks here; the example
/// `mpi_radiation` runs the same check across MPI processes in CI) every
/// tally is bitwise the tally of one rank.
#[test]
fn case4_ranks_are_bitwise_one_rank() {
    let (m, t, kappa, w) = mixed_case();
    let g = Geometry::compute(&m);
    let s = Settings {
        bundles: 40_000,
        ..Settings::default()
    };
    let one = Domain::whole(&m, &g).solve(&t, &kappa, &w, &s).unwrap();
    check_ledger(&one);
    let first = m.n_interior_faces();
    let mut cases = vec![(2, Method::Rcb), (3, Method::Rcb), (4, Method::Rcb)];
    if cfg!(feature = "metis") {
        cases.push((3, Method::Metis));
    }
    for (ranks, method) in cases {
        let part = partition(&m, &g, ranks, method).unwrap();
        let parts = ThreadComm::run(ranks, |c| {
            let comm: &dyn Communicator = &*c;
            let local = extract(&m, &part, comm.rank(), 1).unwrap();
            let lg = local.geometry(comm);
            let own = local.mesh.n_cells;
            let pos: Vec<usize> = local.cell_position[..own]
                .iter()
                .map(|&p| p as usize)
                .collect();
            let lt: Vec<f64> = pos.iter().map(|&p| t[p]).collect();
            let lk: Vec<f64> = pos.iter().map(|&p| kappa[p]).collect();
            // walls by boundary name: the local mesh may hold other boundaries
            let lw: Vec<Wall> = local
                .mesh
                .boundaries
                .iter()
                .map(|b| w[m.boundaries.iter().position(|x| x.name == b.name).unwrap()])
                .collect();
            // one device per run: not on ranks
            assert!(
                Domain::rank(&local, &lg, comm)
                    .on_device(Box::new(CpuDevicePath))
                    .is_err()
            );
            let sol = Domain::rank(&local, &lg, comm)
                .solve(&lt, &lk, &lw, &s)
                .unwrap();
            let lfirst = local.mesh.n_interior_faces();
            let faces: Vec<u64> = (lfirst..local.mesh.n_faces())
                .map(|f| local.mesh.face_id[f])
                .collect();
            (pos, faces, sol)
        });
        let mut seen = 0;
        let mut handed = 0;
        for (pos, faces, sol) in &parts {
            assert_eq!(sol.ledger, one.ledger, "{ranks} ranks, {method:?}");
            assert_eq!(sol.bundles, one.bundles);
            handed = sol.hand_overs;
            for (c, &p) in pos.iter().enumerate() {
                assert_eq!(sol.absorbed[c], one.absorbed[p], "cell {p}, {ranks} ranks");
                assert_eq!(sol.emitted[c], one.emitted[p]);
                assert_eq!(sol.source[c].to_bits(), one.source[p].to_bits());
                assert_eq!(sol.source_error[c].to_bits(), one.source_error[p].to_bits());
            }
            for (f, &id) in faces.iter().enumerate() {
                let gf = (first..m.n_faces()).find(|&x| m.face_id[x] == id).unwrap() - first;
                assert_eq!(sol.wall_absorbed[f], one.wall_absorbed[gf], "face {id}");
                assert_eq!(sol.wall_emitted[f], one.wall_emitted[gf]);
                assert_eq!(sol.wall_net[f].to_bits(), one.wall_net[gf].to_bits());
            }
            seen += pos.len();
        }
        assert_eq!(seen, m.n_cells);
        assert!(handed > 1000, "{handed} hand-overs on {ranks} ranks");
    }
}
