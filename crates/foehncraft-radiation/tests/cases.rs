//! Verification of the radiation call R1a (step 29, section 3): cases 1 to
//! 4 against exact solutions, the exact energy account, and bitwise the
//! same result on any number of threads and ranks.

use foehncraft_mesh::decompose::extract;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_mesh::{Geometry, Mesh};
use foehncraft_parallel::{Communicator, ThreadComm};
use foehncraft_radiation::exact::{parallel_rectangles, slab_source, slab_wall_flux};
use foehncraft_radiation::{Domain, SIGMA, Settings, Solution, Wall};

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
        let sol = Domain::whole(&m, &g).solve(&vec![t; m.n_cells], &vec![kappa; m.n_cells], &w, &s);
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
        let sol =
            Domain::whole(&m, &g).solve(&vec![t; m.n_cells], &kappa, &w, &Settings::default());
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
        let sol = Domain::whole(&m, &g).solve(&vec![t; m.n_cells], &vec![0.0; m.n_cells], &w, &s);
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
            .install(|| Domain::whole(&m, &g).solve(&t, &kappa, &w, &s))
    };
    let one = run(1);
    check_ledger(&one);
    assert!(one.ledger.escaped > 0 && one.ledger.absorbed_walls > 0 && one.ledger.absorbed_gas > 0);
    assert_ne!(one.ledger.roulette, 0, "roulette played");
    for threads in [2, 4] {
        assert_eq!(run(threads), one, "{threads} threads");
    }
    // another call number gives other bundles
    let other = Domain::whole(&m, &g).solve(&t, &kappa, &w, &Settings { call: 1, ..s });
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
    let one = Domain::whole(&m, &g).solve(&t, &kappa, &w, &s);
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
            let sol = Domain::rank(&local, &lg, comm).solve(&lt, &lk, &lw, &s);
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
