//! The box furnace of the radiation gallery (R1a, D255): a flame of
//! prescribed temperature and absorption in a chamber with refractory
//! walls and a row of cooled tubes under the roof.
//!
//! ```text
//! pip install gmsh
//! python3 examples/meshes/furnace.py --out output/furnace.msh
//! cargo run --release -p foehncraft-radiation --example furnace -- \
//!     output/furnace.msh output/furnace [BUNDLES]
//! ```
//!
//! The flame leaves the burner disc in the front wall along x and rises
//! slowly; its temperature peaks at 1 950 K about 1.2 m from the burner
//! and falls to the 1 250 K of the flue gas, its absorption coefficient
//! (soot) rises from the 0.3 1/m of the gas to 1.1 1/m in the core. The
//! walls are grey (ε 0.75) at 1 150 K, the tubes grey (ε 0.85) at 620 K,
//! the burner disc grey (ε 0.8) at 1 100 K; the flue gas leaves through the
//! opening in the back wall (open boundary).
//!
//! Writes into the output directory `furnace.vtkhdf` (the source, its
//! error, the temperature and the absorption coefficient in the cells; the
//! incident and net flux on the walls), `cells.csv` and `walls.csv` (for
//! `tools/plot/radiation_furnace.py`) and prints the energy account and the
//! power the tubes take up.

use foehncraft_mesh::Geometry;
use foehncraft_radiation::output::{cell_csv, wall_csv};
use foehncraft_radiation::{Domain, Settings, Wall};
use std::time::Instant;

/// Gas temperature (K) and absorption coefficient (1/m) at a point.
fn flame(x: [f64; 3]) -> (f64, f64) {
    let s = x[0].max(0.0);
    // the axis rises slowly from the burner centre
    let axis = [1.0, 0.9 + 0.06 * s * s];
    let r2 = (x[1] - axis[0]).powi(2) + (x[2] - axis[1]).powi(2);
    let width = 0.22 + 0.16 * s;
    let along = (s / 1.2) * (1.0 - s / 1.2).exp(); // 1 at s = 1.2 m
    let core = along * (-r2 / (width * width)).exp();
    (1250.0 + 700.0 * core, 0.3 + 0.8 * core)
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    anyhow::ensure!(
        args.len() >= 3,
        "usage: furnace MESH.msh OUT_DIR [BUNDLES] (mesh from examples/meshes/furnace.py)"
    );
    let bundles: u64 = args.get(3).map_or(Ok(4_000_000), |b| b.parse())?;
    let out = std::path::Path::new(&args[2]);
    std::fs::create_dir_all(out)?;

    let mesh = foehncraft_mesh::format::read(&args[1], None)
        .map_err(anyhow::Error::msg)?
        .mesh;
    let geo = Geometry::compute(&mesh);
    let (t, kappa): (Vec<f64>, Vec<f64>) = geo.cell_centroid.iter().map(|&x| flame(x)).unzip();
    let walls: Vec<Wall> = mesh
        .boundaries
        .iter()
        .map(|b| match b.name.as_str() {
            "burner" => Ok(Wall::grey(0.8, 1100.0)),
            "front" | "walls" => Ok(Wall::grey(0.75, 1150.0)),
            "tubes" => Ok(Wall::grey(0.85, 620.0)),
            "outlet" => Ok(Wall::open()),
            other => Err(anyhow::anyhow!("unknown boundary {other}")),
        })
        .collect::<anyhow::Result<_>>()?;
    let s = Settings {
        bundles,
        ..Settings::default()
    };
    let start = Instant::now();
    let sol = Domain::whole(&mesh, &geo).solve(&t, &kappa, &walls, &s)?;
    let seconds = start.elapsed().as_secs_f64();

    let l = &sol.ledger;
    let q = sol.quantum;
    let kw = |n: i64| n as f64 * q * 1e-3;
    println!(
        "furnace: {} cells, {} bundles on {} threads in {seconds:.2} s ({:.0} ns per bundle)",
        mesh.n_cells,
        sol.bundles,
        rayon::current_num_threads(),
        seconds * 1e9 / sol.bundles as f64
    );
    println!(
        "  emitted {:.1} kW = gas {:.1} + walls {:.1} + escaped {:.1} + roulette {:.3} kW; closes: {}",
        kw(l.emitted),
        kw(l.absorbed_gas),
        kw(l.absorbed_walls),
        kw(l.escaped),
        kw(l.roulette),
        l.closes()
    );
    let first = mesh.n_interior_faces();
    println!("  boundary   area (m²)  incident (kW)  net in (kW)  net (kW/m²)");
    for b in &mesh.boundaries {
        let (mut inc, mut net, mut err, mut area) = (0.0, 0.0, 0.0, 0.0);
        for f in b.range() {
            let a = foehncraft_mesh::geometry::norm(geo.face_area[f]);
            inc += sol.wall_incident[f - first] * a;
            net += sol.wall_net[f - first] * a;
            err += (sol.wall_net_error[f - first] * a).powi(2);
            area += a;
        }
        println!(
            "  {:<8} {area:10.2} {:14.1} {:8.1} ± {:4.1} {:10.1}",
            b.name,
            inc * 1e-3,
            net * 1e-3,
            err.sqrt() * 1e-3,
            net * 1e-3 / area
        );
    }
    std::fs::write(
        out.join("cells.csv"),
        cell_csv(
            &mesh,
            &geo,
            &sol,
            &[("temperature", &t), ("absorption", &kappa)],
        ),
    )?;
    std::fs::write(out.join("walls.csv"), wall_csv(&mesh, &geo, &sol))?;
    #[cfg(feature = "hdf5")]
    {
        use foehncraft_mesh::vtkhdf::CellField;
        foehncraft_radiation::output::write_vtkhdf(
            out.join("furnace.vtkhdf"),
            &mesh,
            &sol,
            &[
                CellField {
                    name: "temperature",
                    values: &t,
                },
                CellField {
                    name: "absorption coefficient",
                    values: &kappa,
                },
            ],
        )?;
    }
    println!("  wrote {}", out.display());
    Ok(())
}
