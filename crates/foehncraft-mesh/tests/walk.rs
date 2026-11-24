//! The cell walk (D257) and its hand-over between ranks (D258): exact
//! lengths per cell, the exit point of the ray cast, and bitwise the walk
//! of one rank on every decomposition.

use foehncraft_core::rng;
use foehncraft_mesh::decompose::{LocalMesh, extract};
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::geometry::{norm, sub};
use foehncraft_mesh::handover::HandOver;
use foehncraft_mesh::locate::{Locator, Walk, WalkEnd, step_limit};
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_mesh::query::{Bvh, SurfaceFaces};
use foehncraft_mesh::{Geometry, Mesh};

/// Turns every point about an oblique axis, so no wall lies in a
/// coordinate plane.
fn rotate(m: &mut Mesh) {
    let (s, c) = 0.7f64.sin_cos();
    let (s2, c2) = 0.4f64.sin_cos();
    for p in &mut m.points {
        let [x, y, z] = *p;
        let (x, y) = (c * x - s * y, s * x + c * y);
        let (y, z) = (c2 * y - s2 * z, s2 * y + c2 * z);
        *p = [x, y, z];
    }
}

fn mesh(kind: Kind, n: [usize; 3], perturb: f64) -> Mesh {
    let mut s = BoxSpec::new(kind, 3, n);
    s.perturb = perturb;
    s.upper = [2.0, 1.0, 0.5];
    s.generate().unwrap()
}

/// A unit direction from three keyed uniform numbers.
fn direction(i: u64) -> [f64; 3] {
    let u = rng::uniform([21, i, 0]);
    let v = rng::uniform([21, i, 1]);
    let z = 2.0 * u - 1.0;
    let r = (1.0 - z * z).sqrt();
    let phi = 2.0 * std::f64::consts::PI * v;
    [r * phi.cos(), r * phi.sin(), z]
}

#[test]
fn intervals_tile_the_path_and_end_where_the_ray_cast_ends() {
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        let mut m = mesh(kind, [7, 5, 4], 0.25);
        rotate(&mut m);
        let g = Geometry::compute(&m);
        let loc = Locator::new(&m, &g);
        let all: Vec<usize> = (0..m.boundaries.len()).collect();
        let bvh = Bvh::build(SurfaceFaces::of_boundaries(&m, &g, &all));
        let limit = step_limit(m.n_cells);
        for i in 0..400u64 {
            let c = (i as usize * 7919) % m.n_cells;
            let from = loc.start(c);
            let d = direction(i);
            // far outside the box
            let to: [f64; 3] = std::array::from_fn(|k| from.x[k] + 10.0 * d[k]);
            let mut w = loc.start_walk(from, to);
            let mut seen: Vec<(u32, f64, f64)> = Vec::new();
            let end = loc.walk(&mut w, limit, |cell, a, b| {
                seen.push((cell, a, b));
                true
            });
            let WalkEnd::Boundary { face } = end else {
                panic!("{kind:?} ray {i}: {end:?}");
            };
            assert!(face as usize >= m.n_interior_faces());
            // the intervals tile [0, λ] without gap or overlap
            assert_eq!(seen[0].0 as usize, c);
            assert_eq!(seen[0].1, 0.0);
            for p in seen.windows(2) {
                assert_eq!(p[0].2.to_bits(), p[1].1.to_bits(), "{kind:?} ray {i}");
                assert!(p[0].2 >= p[0].1);
            }
            assert_eq!(seen.last().unwrap().2.to_bits(), w.lambda.to_bits());
            assert_eq!(seen.len(), w.crossed as usize + 1);
            // the exit point is the first hit of the ray on the walls
            let h = bvh.ray_first_hit(from.x, d).expect("closed surface");
            let e = norm(sub(w.x, h.point));
            assert!(
                e < 1e-12,
                "{kind:?} ray {i}: exit {:?} against {:?}",
                w.x,
                h.point
            );
            let len = norm(sub(to, from.x));
            assert!((w.lambda * len - h.t).abs() < 1e-12, "{kind:?} ray {i}");
        }
    }
}

#[test]
fn lengths_per_cell_are_the_chords_of_a_hex_mesh() {
    // cells 0.25 x 0.2 x 0.125; a ray along x crosses whole cells
    let m = mesh(Kind::Hex, [8, 5, 4], 0.0);
    let g = Geometry::compute(&m);
    let loc = Locator::new(&m, &g);
    let x0 = [0.01, 0.37, 0.31];
    let c = loc.locate_search(x0).unwrap();
    let from = loc.position_in(c, x0).unwrap();
    let mut w = loc.start_walk(from, [5.0, 0.37, 0.31]);
    let mut lengths = Vec::new();
    let end = loc.walk(&mut w, step_limit(m.n_cells), |_, a, b| {
        lengths.push((b - a) * 4.99);
        true
    });
    assert!(matches!(end, WalkEnd::Boundary { .. }));
    let lengths: Vec<f64> = lengths.into_iter().filter(|&l| l > 0.0).collect();
    assert_eq!(lengths.len(), 8);
    assert!((lengths[0] - 0.24).abs() < 1e-14, "{}", lengths[0]);
    for l in &lengths[1..] {
        assert!((l - 0.25).abs() < 1e-14, "{l}");
    }
}

/// One walk's record: (ray, cell position, λ₀ bits, λ₁ bits) per cell, and
/// its end (ray, face identity, the bits of λ and x, steps, crossings).
#[derive(Default, Debug, PartialEq)]
struct Log {
    cells: Vec<(u64, u32, u64, u64)>,
    ends: Vec<(u64, u64, [u64; 4], u32, u32)>,
}

impl Log {
    fn sort(&mut self) {
        self.cells.sort_unstable();
        self.ends.sort_unstable();
    }
}

fn rays(m: &Mesh) -> impl Iterator<Item = (u64, usize, [f64; 3])> + '_ {
    (0..m.n_cells as u64)
        .filter(|i| i % 3 == 0)
        .map(|i| (i, i as usize, direction(i)))
}

/// Walks every ray on the whole mesh.
fn one_rank(m: &Mesh, g: &Geometry) -> Log {
    let loc = Locator::new(m, g);
    let limit = step_limit(m.n_cells);
    let mut log = Log::default();
    for (ray, c, d) in rays(m) {
        let from = loc.start(c);
        let to = std::array::from_fn(|k| from.x[k] + 10.0 * d[k]);
        let mut w = loc.start_walk(from, to);
        let end = loc.walk(&mut w, limit, |cell, a, b| {
            log.cells.push((ray, cell, a.to_bits(), b.to_bits()));
            true
        });
        let WalkEnd::Boundary { face } = end else {
            panic!("{end:?}")
        };
        log.ends.push((
            ray,
            m.face_id[face as usize],
            [w.lambda, w.x[0], w.x[1], w.x[2]].map(f64::to_bits),
            w.steps,
            w.crossed,
        ));
    }
    log.sort();
    log
}

/// Walks every ray on the local meshes of a decomposition held in one
/// process, handing walks over at processor faces until none is left.
fn ranks(global: &Mesh, locals: &[LocalMesh]) -> Log {
    let geos: Vec<Geometry> = locals.iter().map(|l| Geometry::compute(&l.mesh)).collect();
    let locs: Vec<Locator> = locals
        .iter()
        .zip(&geos)
        .map(|(l, g)| Locator::new(&l.mesh, g))
        .collect();
    let hand: Vec<HandOver> = locals.iter().map(|l| HandOver::new(l, 1)).collect();
    let limit = step_limit(global.n_cells);
    let gg = Geometry::compute(global);
    let gloc = Locator::new(global, &gg);
    // each ray starts on the rank that owns its cell
    let mut queue: Vec<Vec<(Walk, u64)>> = vec![Vec::new(); locals.len()];
    for (ray, c, d) in rays(global) {
        let (r, lc) = locals
            .iter()
            .enumerate()
            .find_map(|(r, l)| {
                l.cell_position[..l.mesh.n_cells]
                    .iter()
                    .position(|&p| p as usize == c)
                    .map(|lc| (r, lc))
            })
            .unwrap();
        let from = locs[r].start(lc);
        assert_eq!(from.x, gloc.start(c).x, "start of ray {ray}");
        let to = std::array::from_fn(|k| from.x[k] + 10.0 * d[k]);
        queue[r].push((locs[r].start_walk(from, to), ray));
    }
    let mut log = Log::default();
    let mut outboxes: Vec<Vec<Vec<f64>>> = hand.iter().map(|h| h.outbox()).collect();
    let mut rounds = 0;
    while queue.iter().any(|q| !q.is_empty()) {
        rounds += 1;
        assert!(rounds < 100, "hand-over does not end");
        for r in 0..locals.len() {
            let pos = &locals[r].cell_position;
            for (mut w, ray) in std::mem::take(&mut queue[r]) {
                let end = locs[r].walk(&mut w, limit, |cell, a, b| {
                    log.cells
                        .push((ray, pos[cell as usize], a.to_bits(), b.to_bits()));
                    true
                });
                match end {
                    WalkEnd::Boundary { face } => log.ends.push((
                        ray,
                        locals[r].mesh.face_id[face as usize],
                        [w.lambda, w.x[0], w.x[1], w.x[2]].map(f64::to_bits),
                        w.steps,
                        w.crossed,
                    )),
                    WalkEnd::Processor { .. } => {
                        hand[r].pack(&w, &[ray as f64], &mut outboxes[r]);
                    }
                    e => panic!("{e:?}"),
                }
            }
        }
        let inboxes = HandOver::exchange_all(&hand, &mut outboxes);
        for (r, inbox) in inboxes.iter().enumerate() {
            for (w, payload) in hand[r].unpack(inbox) {
                queue[r].push((w, payload[0] as u64));
            }
        }
    }
    assert!(rounds > 2, "{rounds} rounds: walks were handed over");
    log.sort();
    log
}

#[test]
fn walks_handed_over_between_ranks_are_the_walks_of_one_rank() {
    for kind in [Kind::Tet, Kind::Poly, Kind::Mixed] {
        let mut m = mesh(kind, [9, 6, 5], 0.25);
        rotate(&mut m);
        let g = Geometry::compute(&m);
        let reference = one_rank(&m, &g);
        assert!(
            reference.cells.len() > 4 * reference.ends.len(),
            "{} {}",
            reference.cells.len(),
            reference.ends.len()
        );
        let mut cases = vec![(2, Method::Rcb), (3, Method::Rcb), (4, Method::Rcb)];
        if cfg!(feature = "metis") {
            cases.push((3, Method::Metis));
        }
        for (parts, method) in cases {
            let part = partition(&m, &g, parts, method).unwrap();
            let locals: Vec<LocalMesh> = (0..parts)
                .map(|r| extract(&m, &part, r, 1).unwrap())
                .collect();
            let log = ranks(&m, &locals);
            assert_eq!(
                log.ends, reference.ends,
                "{kind:?}, {parts} ranks, {method:?}"
            );
            assert_eq!(
                log.cells, reference.cells,
                "{kind:?}, {parts} ranks, {method:?}"
            );
        }
    }
}
