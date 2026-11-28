//! Gas-phase chemistry of the cells on the GPU (SF2b): the reactors of a
//! row of cells (stoichiometric methane and air from 900 to 2400 K,
//! fresh, half and fully burnt) integrated over one flow time step on the
//! CPU threads and on the GPU, the end states and step counts compared bit
//! for bit, with the time per cell on both, for GRI-Mech 3.0 and the
//! global mechanism of Jones & Lindstedt. A second row of each mechanism
//! gives every cell its own interval (the fine-structure residence times
//! τ* of the eddy dissipation concept, SF3: 10 µs to 1 ms over the row).
//!
//! ```text
//! python3 gpu/ptx-kernels/build.py
//! cargo run --release -p foehncraft-cuda --features cuda --example chemistry_gpu [PTX]
//! ```
//!
//! Prints a Markdown table and exits with status 1 when a GPU result
//! differs from its CPU result.

use foehncraft_flow::chemistry::Mechanism;
use foehncraft_flow::chemistry::device::react_cpu;
use foehncraft_kernels::chemistry::{Settings, Stats};
use std::path::PathBuf;
use std::time::Instant;

fn states(m: &Mechanism, cells: usize) -> Vec<f64> {
    let n = m.len();
    let mut z = Vec::with_capacity(cells * (n + 1));
    for c in 0..cells {
        let s = c as f64 / (cells - 1).max(1) as f64;
        // burnt fraction 0, 0.5, 1 in turn; temperature over the row
        let b = (c % 3) as f64 / 2.0;
        let mut row = vec![0.0; n + 1];
        let set = |row: &mut [f64], name: &str, y: f64| {
            if let Some(k) = m.index(name) {
                row[k] = y;
            }
        };
        set(&mut row, "CH4", 0.05519 * (1.0 - b));
        set(&mut row, "O2", 0.2200 * (1.0 - b));
        set(&mut row, "CO2", 0.1514 * b);
        set(&mut row, "H2O", 0.1240 * b);
        let rest: f64 = row[..n].iter().sum();
        set(&mut row, "N2", 1.0 - rest);
        row[n] = 900.0 + 1500.0 * s;
        z.extend(row);
    }
    z
}

fn main() {
    let ptx = std::env::args().nth(1).map_or_else(
        || PathBuf::from("target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx"),
        PathBuf::from,
    );
    foehncraft_cuda::register();
    let set = Settings::default();
    let (p, dt, tmin, cells) = (101_325.0, 1e-5, 400.0, 4096);
    println!("| mechanism | interval | cells | CPU µs/cell | GPU µs/cell | bitwise |");
    println!("|---|---|---|---|---|---|");
    let mut ok = true;
    // the residence times of the fine structures (SF3), per cell
    let taus: Vec<f64> = (0..cells)
        .map(|c| 1e-5 * 100f64.powf(c as f64 / (cells - 1) as f64))
        .collect();
    for (name, per_cell) in [
        ("gri30", false),
        ("jones-lindstedt-4", false),
        ("gri30", true),
        ("jones-lindstedt-4", true),
    ] {
        let dt_cell = per_cell.then_some(&taus[..]);
        let m = Mechanism::load(name, None).expect("mechanism");
        let prog = m.program();
        let start = states(&m, cells);
        let mut cpu = start.clone();
        let mut st_cpu = vec![Stats::default(); cells];
        let t0 = Instant::now();
        react_cpu(&prog, p, dt, dt_cell, tmin, &set, &mut cpu, &mut st_cpu);
        let t_cpu = t0.elapsed().as_secs_f64();
        let mut gpu = start.clone();
        let mut st_gpu = vec![Stats::default(); cells];
        let mut dev = foehncraft_flow::chemistry::device::create_gpu("cuda", &prog, 0, &ptx)
            .expect("GPU device");
        let t0 = Instant::now();
        dev.react(p, dt, dt_cell, tmin, &set, &mut gpu, &mut st_gpu)
            .expect("GPU launch");
        let t_gpu = t0.elapsed().as_secs_f64();
        let same = cpu
            .iter()
            .zip(&gpu)
            .all(|(a, b)| a.to_bits() == b.to_bits())
            && st_cpu == st_gpu;
        if !same {
            let worst = cpu
                .iter()
                .zip(&gpu)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            eprintln!("{name}: the GPU differs from the CPU (largest difference {worst:e})");
        }
        ok &= same;
        println!(
            "| {name} | {} | {cells} | {:.1} | {:.1} | {} |",
            if per_cell { "τ* per cell" } else { "Δt" },
            1e6 * t_cpu / cells as f64,
            1e6 * t_gpu / cells as f64,
            if same { "yes" } else { "**no**" }
        );
    }
    if !ok {
        std::process::exit(1);
    }
}
