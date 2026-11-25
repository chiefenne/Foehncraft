//! Verification case 6 of the radiation (R1a): the error of scrambled
//! Sobol against pseudo-random emission points, against the number of
//! bundles, on the slab (case 1, τ = 1, the wall flux) and on the parallel
//! plates (case 3, the configuration factor). Prints a CSV table of the
//! root-mean-square relative error over 16 independent calls per point.
//!
//! ```text
//! cargo run --release -p foehncraft-radiation --example sobol_study > output/sobol.csv
//! ```

use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::{Geometry, Mesh};
use foehncraft_radiation::exact::{parallel_rectangles, slab_wall_flux};
use foehncraft_radiation::{Domain, Points, SIGMA, Settings, Wall};

fn mesh(n: [usize; 3], upper: [f64; 3]) -> Mesh {
    let mut s = BoxSpec::new(Kind::Hex, 3, n);
    s.upper = upper;
    s.generate().expect("mesh")
}

fn walls(m: &Mesh, wall: impl Fn(&str) -> Wall) -> Vec<Wall> {
    m.boundaries.iter().map(|b| wall(&b.name)).collect()
}

/// Area-weighted sum of `v` over the boundary `name`, and its area.
fn on(m: &Mesh, g: &Geometry, name: &str, v: impl Fn(usize) -> f64) -> (f64, f64) {
    let b = m
        .boundaries
        .iter()
        .find(|b| b.name == name)
        .expect("boundary");
    let first = m.n_interior_faces();
    b.range().fold((0.0, 0.0), |(s, a), f| {
        let af = foehncraft_mesh::geometry::norm(g.face_area[f]);
        (s + v(f - first) * af, a + af)
    })
}

fn main() {
    let t = 1000.0f64;
    let slab = mesh([10, 2, 2], [1.0, 0.5, 0.5]);
    let gs = Geometry::compute(&slab);
    let ws = walls(&slab, |name| match name {
        "xmin" | "xmax" => Wall::black(0.0),
        _ => Wall::mirror(),
    });
    let plates = mesh([2, 4, 4], [1.0; 3]);
    let gp = Geometry::compute(&plates);
    let wp = walls(&plates, |name| match name {
        "xmin" => Wall::black(t),
        "xmax" => Wall::black(0.0),
        _ => Wall::open(),
    });
    let exact_slab = SIGMA * t.powi(4) * slab_wall_flux(1.0);
    let exact_plates = parallel_rectangles(1.0, 1.0, 1.0);
    let calls = 16u64;
    println!("bundles,slab_sobol,slab_random,plates_sobol,plates_random");
    for k in 0..9 {
        let n = 2_000u64 << k;
        let mut row = vec![n as f64];
        let mut errors = [[0.0f64; 2]; 2];
        for (j, points) in [Points::Sobol, Points::Random].into_iter().enumerate() {
            for call in 0..calls {
                let s = Settings {
                    bundles: n,
                    call,
                    points,
                    ..Settings::default()
                };
                let sol = Domain::whole(&slab, &gs)
                    .solve(&vec![t; slab.n_cells], &vec![1.0; slab.n_cells], &ws, &s)
                    .expect("slab");
                let (q, a) = on(&slab, &gs, "xmin", |f| sol.wall_net[f]);
                errors[0][j] += (q / a / exact_slab - 1.0).powi(2);
                let sol = Domain::whole(&plates, &gp)
                    .solve(
                        &vec![t; plates.n_cells],
                        &vec![0.0; plates.n_cells],
                        &wp,
                        &s,
                    )
                    .expect("plates");
                let (cold, _) = on(&plates, &gp, "xmax", |f| sol.wall_absorbed[f] as f64);
                let (hot, _) = on(&plates, &gp, "xmin", |f| sol.wall_emitted[f] as f64);
                // the faces have equal areas, so the area weights cancel
                errors[1][j] += (cold / hot / exact_plates - 1.0).powi(2);
            }
        }
        for case in errors {
            for e in case {
                row.push((e / calls as f64).sqrt());
            }
        }
        let text: Vec<String> = row.iter().map(|v| format!("{v:e}")).collect();
        println!("{}", text.join(","));
    }
}
