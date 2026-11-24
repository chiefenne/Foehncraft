//! Strong scaling of the flow solver on ranks (MPI step W5): one steady
//! laminar channel on hexahedra, a fixed number of SIMPLEC outer
//! iterations, on 1 rank and on several ranks (threads of this process,
//! [`ThreadComm`]) with the cores shared among them; prints the wall time
//! of the iterations (setup excluded) and checks every run bitwise against
//! the first.
//!
//! ```text
//! cargo run --release -p foehncraft-flow --example ranks_scaling -- 200 100 50 10
//! ```
//!
//! Arguments: cells in x, y, z and the outer iterations (default 100 50 25
//! and 10).

use foehncraft_flow::case::FlowCase;
use foehncraft_flow::ranks::{check, solve_one, solve_part};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::decompose::SharedComm;
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_parallel::ThreadComm;

fn main() -> anyhow::Result<()> {
    let a: Vec<usize> = std::env::args().skip(1).filter_map(|x| x.parse().ok()).collect();
    let (nx, ny, nz, its) = match a.as_slice() {
        [x, y, z, i] => (*x, *y, *z, *i),
        _ => (100, 50, 25, 10),
    };
    let text = format!(
        r#"
        case_format = "0.7"
        [mesh.generate]
        type = "hex"
        dimension = 3
        upper = [4.0, 2.0, 1.0]
        cells = [{nx}, {ny}, {nz}]
        [fluid]
        density = 1.0
        viscosity = 0.01
        [boundary.xmin]
        type = "velocity-inlet"
        velocity = [1.0, 0.1, 0.0]
        [boundary.xmax]
        type = "static-pressure-outlet"
        pressure = 101325.0
        [boundary.ymin]
        type = "wall"
        [boundary.ymax]
        type = "wall"
        velocity = [0.5, 0.0, 0.0]
        [boundary.zmin]
        type = "wall"
        [boundary.zmax]
        type = "wall"
        [solver.momentum]
        limiter = "venkatakrishnan"
        "#
    );
    let case = FlowCase::from_toml(&text, std::path::Path::new("."))?;
    let mesh = foehncraft_flow::ranks::CheckCase {
        name: "scaling",
        case: case.clone(),
        iterations: its,
    }
    .mesh()?;
    let geo = Geometry::compute(&mesh);
    let cores = std::thread::available_parallelism().map_or(1, |c| c.get());
    println!(
        "{} cells, {its} outer iterations, {cores} cores",
        mesh.n_cells
    );
    let one = solve_one(&case, &mesh, its, Some(cores))?;
    println!(
        "{:>6} {:>8} {:>12} {:>14} {:>10}",
        "ranks", "threads", "seconds", "s/iteration", "bitwise"
    );
    println!(
        "{:>6} {:>8} {:>12.3} {:>14.4} {:>10}",
        "none", cores, one.seconds, one.seconds / its as f64, "reference"
    );
    let mut counts = vec![1usize];
    let mut r = 2;
    while r <= cores {
        counts.push(r);
        r *= 2;
    }
    let mut ok = true;
    for ranks in counts {
        let threads = (cores / ranks).max(1);
        let part = partition(&mesh, &geo, ranks, Method::Rcb).map_err(anyhow::Error::msg)?;
        let solves = ThreadComm::run(ranks, |c| {
            let comm: SharedComm = c;
            solve_part(&case, &mesh, &part, comm, its, Some(threads))
        });
        let mut seconds: f64 = 0.0;
        let mut bitwise = true;
        for (rank, s) in solves.into_iter().enumerate() {
            let s = s?;
            seconds = seconds.max(s.seconds);
            bitwise &= check(rank, &s, &one).bitwise();
        }
        ok &= bitwise;
        println!(
            "{ranks:>6} {threads:>8} {seconds:>12.3} {:>14.4} {:>10}",
            seconds / its as f64,
            if bitwise { "yes" } else { "NO" }
        );
    }
    anyhow::ensure!(ok, "a run on ranks differs from one rank");
    Ok(())
}
