//! Point search without a hint cell: the hierarchy over cell boxes
//! (`Locator::locate_search`, D108 part 2) against the former uniform grid
//! of cell boxes (`Locator::locate_search_grid`).
//!
//! ```text
//! cargo run --release -p foehncraft-mesh --example locate_bench -- [mesh files]
//! ```
//!
//! Generated meshes (uniform and strongly graded hexahedra, graded
//! tetrahedra and polyhedra) are always run; mesh files given on the
//! command line are added. Queries: every cell centroid (at most 200 000)
//! and as many random points in the bounding box (some outside a non-convex
//! domain). Prints one markdown row per mesh: build times, query times,
//! the mean number of candidate cells per query of the hierarchy, and the
//! number of queries where the two searches disagree (must be 0).

use foehncraft_core::rng;
use foehncraft_mesh::generate::{BoxSpec, Grading, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_mesh::{Geometry, Mesh};
use std::time::Instant;

fn main() {
    let mut meshes: Vec<(String, Mesh)> = Vec::new();
    let make = |kind: Kind, n: [usize; 3], ratio: f64| {
        let mut s = BoxSpec::new(kind, 3, n);
        s.perturb = 0.1;
        s.grading.y = Grading {
            ratio,
            symmetric: true,
        };
        s.generate().expect("mesh")
    };
    meshes.push(("hex 64³ uniform".into(), make(Kind::Hex, [64, 64, 64], 1.0)));
    meshes.push((
        "hex 64³ graded 1e4 in y".into(),
        make(Kind::Hex, [64, 64, 64], 1e4),
    ));
    meshes.push((
        "tet 24³ graded 1e3 in y".into(),
        make(Kind::Tet, [24, 24, 24], 1e3),
    ));
    meshes.push((
        "poly 20³ graded 1e3 in y".into(),
        make(Kind::Poly, [20, 20, 20], 1e3),
    ));
    for path in std::env::args().skip(1) {
        let f = foehncraft_mesh::format::read(&path, None).expect("mesh file");
        let name = std::path::Path::new(&path)
            .file_name()
            .map_or(path.clone(), |n| n.to_string_lossy().into_owned());
        meshes.push((name, f.mesh));
    }
    println!(
        "| Mesh | Cells | Grid build (s) | Tree build (s) | Queries | Grid search (µs/query) | Tree search (µs/query) | Speed-up | Candidates per query | Disagreements |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|");
    for (name, m) in &meshes {
        let g = Geometry::compute(m);
        let (lo, hi) = m.bounding_box();
        let nq = m.n_cells.min(200_000);
        let step = (m.n_cells / nq).max(1);
        let mut q: Vec<[f64; 3]> = (0..nq).map(|i| g.cell_centroid[i * step]).collect();
        for i in 0..nq as u64 {
            q.push(std::array::from_fn(|d| {
                lo[d] + (hi[d] - lo[d]) * rng::uniform([17, i, d as u64])
            }));
        }
        let loc = Locator::new(m, &g);
        // builds (first call)
        let t = Instant::now();
        let _ = loc.locate_search_grid(q[0]);
        let grid_build = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let _ = loc.cell_tree();
        let tree_build = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let a: Vec<Option<usize>> = q.iter().map(|&x| loc.locate_search_grid(x)).collect();
        let grid_q = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let b: Vec<Option<usize>> = q.iter().map(|&x| loc.locate_search(x)).collect();
        let tree_q = t.elapsed().as_secs_f64();
        let differ = a.iter().zip(&b).filter(|(x, y)| x != y).count();
        let mut cand = 0usize;
        for &x in &q {
            loc.cell_tree().for_each_candidate(x, |_| cand += 1);
        }
        let n = q.len() as f64;
        println!(
            "| {name} | {} | {grid_build:.3} | {tree_build:.3} | {} | {:.2} | {:.2} | {:.1} | {:.1} | {differ} |",
            m.n_cells,
            q.len(),
            1e6 * grid_q / n,
            1e6 * tree_q / n,
            grid_q / tree_q,
            cand as f64 / n
        );
    }
}
