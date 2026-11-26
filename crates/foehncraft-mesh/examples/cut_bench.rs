//! Benchmark of the cut cells (L13 stage 2, maintainer rule: benchmark the
//! computer-graphics algorithms too): a sphere cut out of hexahedral,
//! tetrahedral and polyhedral box meshes with three ways of placing the
//! cut surface:
//!
//! - `linear`: crossing points at the zero of the linear interpolant of
//!   the node values (marching cubes / marching tetrahedra);
//! - `exact`: crossing points at the root of the body's own function on
//!   the edge (on the true surface);
//! - `exact+point`: as `exact`, and every wall face fanned around a point
//!   projected onto the surface (the vertex of dual contouring).
//!
//! Prints the relative error of the fluid volume and the wall area, their
//! observed orders and the time of the cut, as a Markdown table.
//!
//! ```sh
//! cargo run --release -p foehncraft-mesh --example cut_bench
//! ```

use foehncraft_mesh::cut::{CutOptions, Sphere, cut};
use foehncraft_mesh::generate::{BoxSpec, Kind};
use std::f64::consts::PI;
use std::time::Instant;

fn main() {
    let r = 0.55;
    let centre = [0.031, -0.017, 0.012];
    let v_body = 4.0 / 3.0 * PI * r * r * r;
    let a_body = 4.0 * PI * r * r;
    println!(
        "| mesh | n | cells | method | volume error | order | area error | order | merged | cut (ms) |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|");
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly] {
        let meshes: Vec<_> = [8usize, 16, 32]
            .iter()
            .map(|&n| {
                let mut s = BoxSpec::new(kind, 3, [n, n, n]);
                s.lower = [-1.0; 3];
                s.upper = [1.0; 3];
                (n, s.generate().expect("mesh"))
            })
            .collect();
        for (name, exact, wall_point) in [
            ("linear", false, false),
            ("exact", true, false),
            ("exact+point", true, true),
        ] {
            let mut prev: Option<(f64, f64)> = None;
            for (n, m) in &meshes {
                let body = Sphere {
                    centre,
                    radius: r,
                    exact,
                };
                let t = Instant::now();
                let c = cut(
                    m,
                    &body,
                    &["sphere".to_string()],
                    &CutOptions {
                        min_fraction: 0.3,
                        wall_point,
                    },
                )
                .expect("cut");
                let ms = t.elapsed().as_secs_f64() * 1e3;
                let ev = (c.report.fluid_volume - (8.0 - v_body)).abs() / v_body;
                let ea = (c.report.wall_area[0] - a_body).abs() / a_body;
                let (ov, oa) = match prev {
                    Some((pv, pa)) => (
                        format!("{:.2}", (pv / ev).log2()),
                        format!("{:.2}", (pa / ea).log2()),
                    ),
                    None => ("".into(), "".into()),
                };
                prev = Some((ev, ea));
                println!(
                    "| {} | {n} | {} | {name} | {ev:.2e} | {ov} | {ea:.2e} | {oa} | {} | {ms:.0} |",
                    kind.name(),
                    c.mesh.n_cells,
                    c.report.cells_merged
                );
            }
        }
    }
}
