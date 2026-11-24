//! Radiation R1a on the GPU (D259): every case once with the bundles on the
//! CPU threads and once on the GPU, the two results compared bit for bit
//! (sources, errors, wall fluxes, tallies, ledger), with the time per
//! bundle on both.
//!
//! ```text
//! python3 gpu/ptx-kernels/build.py
//! cargo run --release -p foehncraft-cuda --features cuda --example radiation_gpu [PTX]
//! ```
//!
//! `python3 -m tools.bench.radiation_cloud` runs it on a cloud GPU (Colab:
//! `tools/notebooks/radiation_colab.ipynb`). Prints a Markdown table and
//! exits with status 1 when a GPU result differs from its CPU result.

use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::{Geometry, Mesh};
use foehncraft_radiation::{Domain, Settings, Solution, Wall};
use std::path::PathBuf;
use std::time::Instant;

struct Case {
    name: &'static str,
    mesh: Mesh,
    t: Vec<f64>,
    kappa: Vec<f64>,
    walls: Vec<Wall>,
    bundles: u64,
}

fn mesh(kind: Kind, n: [usize; 3], upper: [f64; 3], perturb: f64) -> Mesh {
    let mut s = BoxSpec::new(kind, 3, n);
    s.upper = upper;
    s.perturb = perturb;
    s.generate().expect("mesh")
}

fn walls(m: &Mesh, f: impl Fn(&str) -> Wall) -> Vec<Wall> {
    m.boundaries.iter().map(|b| f(&b.name)).collect()
}

fn cases() -> Vec<Case> {
    let mut v = Vec::new();
    // case 1: grey slab between cold black plates, mirror sides
    let m = mesh(Kind::Hex, [40, 2, 2], [1.0, 0.5, 0.5], 0.0);
    v.push(Case {
        name: "slab, τ = 1 (case 1)",
        t: vec![1000.0; m.n_cells],
        kappa: vec![1.0; m.n_cells],
        walls: walls(&m, |n| match n {
            "xmin" | "xmax" => Wall::black(0.0),
            _ => Wall::mirror(),
        }),
        mesh: m,
        bundles: 1_000_000,
    });
    // case 2: isothermal enclosure, grey walls, tetrahedra
    let m = mesh(Kind::Tet, [12, 10, 8], [1.2, 1.0, 0.8], 0.2);
    let g = Geometry::compute(&m);
    v.push(Case {
        name: "isothermal enclosure (case 2)",
        t: vec![800.0; m.n_cells],
        kappa: g.cell_centroid.iter().map(|c| 0.5 + c[0]).collect(),
        walls: walls(&m, |_| Wall::grey(0.5, 800.0)),
        mesh: m,
        bundles: 1_000_000,
    });
    // case 3: transparent medium between opposed plates, open sides
    let m = mesh(Kind::Poly, [6, 8, 8], [1.0, 1.0, 1.0], 0.2);
    v.push(Case {
        name: "parallel plates (case 3)",
        t: vec![1000.0; m.n_cells],
        kappa: vec![0.0; m.n_cells],
        walls: walls(&m, |n| match n {
            "xmin" => Wall::black(1000.0),
            "xmax" => Wall::black(0.0),
            _ => Wall::open(),
        }),
        mesh: m,
        bundles: 1_000_000,
    });
    // everything at once on a larger turned polyhedral mesh
    let mut m = mesh(Kind::Poly, [40, 30, 25], [2.0, 1.5, 1.0], 0.25);
    let (s, c) = 0.7f64.sin_cos();
    for p in &mut m.points {
        let [x, y, z] = *p;
        *p = [c * x - s * y, s * x + c * y, z];
    }
    let g = Geometry::compute(&m);
    v.push(Case {
        name: "mixed walls, polyhedra",
        t: g.cell_centroid
            .iter()
            .map(|x| 600.0 + 300.0 * (x[0] + 0.5 * x[1]).sin().abs())
            .collect(),
        kappa: g.cell_centroid.iter().map(|x| 0.2 + 1.5 * x[2]).collect(),
        walls: walls(&m, |n| match n {
            "xmin" => Wall::grey(0.8, 1200.0),
            "xmax" => Wall::grey(0.4, 500.0),
            "ymin" => Wall::black(300.0),
            "ymax" => Wall::mirror(),
            "zmin" => Wall::grey(0.6, 900.0),
            _ => Wall::open(),
        }),
        mesh: m,
        bundles: 4_000_000,
    });
    v
}

fn main() {
    foehncraft_cuda::register();
    let ptx = std::env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from("target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx")
        });
    let threads = rayon::current_num_threads();
    println!(
        "| case | cells | bundles | CPU ({threads} threads) ns/bundle | GPU ns/bundle | speed-up | ledger closes | bitwise equal |"
    );
    println!("|---|---:|---:|---:|---:|---:|---|---|");
    let mut failed = false;
    let mut device = String::new();
    for c in cases() {
        let g = Geometry::compute(&c.mesh);
        let s = Settings {
            bundles: c.bundles,
            ..Settings::default()
        };
        let run = |d: Domain<'_>| -> (Solution, f64) {
            let t0 = Instant::now();
            let sol = d
                .solve(&c.t, &c.kappa, &c.walls, &s)
                .expect("radiation call");
            (sol, t0.elapsed().as_secs_f64())
        };
        let (cpu, t_cpu) = run(Domain::whole(&c.mesh, &g));
        let gpu = match Domain::whole(&c.mesh, &g).on_gpu(0, &ptx) {
            Ok(d) => d,
            Err(e) => {
                eprintln!("GPU: {e:#}");
                std::process::exit(1);
            }
        };
        if device.is_empty() {
            device = gpu.device_name();
        }
        let (gpu, t_gpu) = run(gpu);
        let equal = gpu == cpu;
        let closes = cpu.ledger.closes() && gpu.ledger.closes();
        failed |= !equal || !closes;
        let per = |t: f64| t * 1e9 / cpu.bundles as f64;
        println!(
            "| {} | {} | {} | {:.0} | {:.0} | {:.1} | {} | {} |",
            c.name,
            c.mesh.n_cells,
            cpu.bundles,
            per(t_cpu),
            per(t_gpu),
            t_cpu / t_gpu,
            if closes { "yes" } else { "**no**" },
            if equal { "yes" } else { "**no**" }
        );
    }
    println!(
        "\nDevice: {device}. Times include the emission on the host and the transfers.\n{}",
        if failed {
            "FAILED: a GPU result differs from the CPU result or a ledger does not close"
        } else {
            "all GPU results bitwise equal to the CPU results"
        }
    );
    if failed {
        std::process::exit(1);
    }
}
