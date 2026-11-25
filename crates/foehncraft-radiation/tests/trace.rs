//! The traced paths of a call are the paths of its bundles: over all
//! bundles, the quanta the paths leave in the cells, on the walls, through
//! the openings and to Russian roulette are the call's tallies, to the
//! quantum.

use foehncraft_mesh::decompose::extract;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_mesh::{Geometry, Mesh};
use foehncraft_parallel::{Communicator, ThreadComm};
use foehncraft_radiation::{Domain, Settings, Source, Trace, Wall};

/// Turned, perturbed polyhedra with grey, black, mirror and open walls,
/// and a gas whose temperature and absorption vary (case 4 of `cases.rs`).
fn mixed_case() -> (Mesh, Vec<f64>, Vec<f64>, Vec<Wall>) {
    let mut s = BoxSpec::new(Kind::Poly, 3, [8, 6, 5]);
    s.upper = [2.0, 1.5, 1.0];
    s.perturb = 0.25;
    let mut m = s.generate().unwrap();
    let (s, c) = 0.7f64.sin_cos();
    for p in &mut m.points {
        let [x, y, z] = *p;
        *p = [c * x - s * y, s * x + c * y, z];
    }
    let g = Geometry::compute(&m);
    let t = g
        .cell_centroid
        .iter()
        .map(|x| 600.0 + 300.0 * (x[0] + 0.5 * x[1]).sin().abs())
        .collect();
    let kappa = g.cell_centroid.iter().map(|x| 0.2 + 1.5 * x[2]).collect();
    let w = m
        .boundaries
        .iter()
        .map(|b| match b.name.as_str() {
            "xmin" => Wall::grey(0.8, 1200.0),
            "xmax" => Wall::grey(0.4, 500.0),
            "ymin" => Wall::black(300.0),
            "ymax" => Wall::mirror(),
            "zmin" => Wall::grey(0.6, 900.0),
            _ => Wall::open(),
        })
        .collect();
    (m, t, kappa, w)
}

#[test]
fn traced_paths_carry_the_tallies_of_the_call() {
    let (m, t, kappa, w) = mixed_case();
    let g = Geometry::compute(&m);
    let s = Settings {
        bundles: 20_000,
        ..Settings::default()
    };
    let d = Domain::whole(&m, &g);
    let sol = d.solve(&t, &kappa, &w, &s).unwrap();
    let traced = d.trace(&t, &kappa, &w, &s, |_| u32::MAX).unwrap();
    assert_eq!(traced.quantum, sol.quantum);
    assert_eq!(traced.paths.len() as i64, sol.bundles);

    let nb = sol.wall_absorbed.len();
    let (mut cells, mut walls) = (vec![0i64; m.n_cells], vec![0i64; nb]);
    let (mut escaped, mut roulette, mut lost) = (0i64, 0i64, 0i64);
    let mut kinds = [0usize; 6];
    for p in &traced.paths {
        for &(c, q) in &p.cells {
            cells[c as usize] += q;
        }
        for &(f, _, q) in &p.walls {
            walls[f as usize] += q;
        }
        escaped += p.escaped;
        roulette += p.roulette;
        lost += p.lost;
        // every quantum of the bundle is accounted for
        let left: i64 = p.cells.iter().map(|c| c.1).sum::<i64>()
            + p.walls.iter().map(|w| w.2).sum::<i64>()
            + p.escaped
            + p.roulette
            + p.lost;
        assert_eq!(left, p.w0, "bundle {:?} {}", p.source, p.index);
        // a path ends where its bundle ends, and only there
        let (last, rest) = p.points.split_last().expect("a path has points");
        assert!(
            matches!(last.at, Trace::Absorbed | Trace::Escaped | Trace::Lost),
            "{last:?}"
        );
        match last.at {
            Trace::Absorbed => assert_eq!(last.w, 0),
            Trace::Escaped => assert_eq!(last.w, p.escaped),
            _ => {}
        }
        assert!(
            rest.iter()
                .all(|x| matches!(x.at, Trace::Cell | Trace::Wall | Trace::Mirror))
        );
        assert!(rest.iter().all(|x| x.w > 0));
        for x in &p.points {
            kinds[x.at as usize] += 1;
        }
    }
    assert_eq!(cells, sol.absorbed);
    assert_eq!(walls, sol.wall_absorbed);
    assert_eq!(escaped, sol.ledger.escaped);
    assert_eq!(roulette, sol.ledger.roulette);
    assert_eq!(lost, sol.ledger.lost);
    // every kind of event occurs in this case (no bundle is lost)
    for (k, n) in kinds.iter().enumerate().take(5) {
        assert!(*n > 0, "no point of kind {k}");
    }
    assert_eq!(kinds[Trace::Lost as usize], 0);
}

#[test]
fn chosen_bundles_only_and_not_on_ranks() {
    let (m, t, kappa, w) = mixed_case();
    let g = Geometry::compute(&m);
    let s = Settings {
        bundles: 5_000,
        ..Settings::default()
    };
    let d = Domain::whole(&m, &g);
    // the first two bundles of the cells below z = 0.3 m and of no wall
    let low = |c: usize| g.cell_centroid[c][2] < 0.3;
    let traced = d
        .trace(&t, &kappa, &w, &s, |src| match src {
            Source::Cell(c) if low(c) => 2,
            _ => 0,
        })
        .unwrap();
    assert!(!traced.paths.is_empty());
    for p in &traced.paths {
        let Source::Cell(c) = p.source else {
            panic!("a wall bundle was traced")
        };
        assert!(low(c) && p.index < 2);
    }
    // the same bundles in a full trace have the same paths
    let all = d.trace(&t, &kappa, &w, &s, |_| u32::MAX).unwrap();
    for p in &traced.paths {
        assert!(all.paths.contains(p));
    }

    // a rank's part of the mesh is refused (before any message)
    let part = partition(&m, &g, 2, Method::Rcb).unwrap();
    let errors = ThreadComm::run(2, |c| {
        let comm: &dyn Communicator = &*c;
        let local = extract(&m, &part, comm.rank(), 1).unwrap();
        let lg = local.geometry(comm);
        let own = local.mesh.n_cells;
        let lw: Vec<Wall> = local
            .mesh
            .boundaries
            .iter()
            .map(|b| w[m.boundaries.iter().position(|x| x.name == b.name).unwrap()])
            .collect();
        Domain::rank(&local, &lg, comm)
            .trace(&vec![900.0; own], &vec![0.5; own], &lw, &s, |_| 1)
            .unwrap_err()
            .to_string()
    });
    for e in errors {
        assert!(e.contains("whole mesh"), "{e}");
    }
}
