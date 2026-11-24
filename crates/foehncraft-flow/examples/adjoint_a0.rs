//! L17 stage A0: the one steady residual R(q) = (R_u, R_c, R_ṁ) of
//! SIMPLEC and the coupled solver along both runs of one case, and the
//! fields of the converged state (gallery 267).
//!
//! ```text
//! cargo run --release -p foehncraft-flow --example adjoint_a0 -- \
//!     examples/cases/flow/adjoint_residual_cavity.toml output/adjoint_a0
//! python3 -m tools.plot.adjoint_a0 docs/figures/gallery/adjoint_a0_residual.png output/adjoint_a0
//! ```
//!
//! Writes `history_simplec.csv` and `history_coupled.csv` (outer iteration
//! and the three relative norms of R after it) and `cells.csv` (cell
//! centroid, velocity and the momentum residual of the start fields and of
//! the converged SIMPLEC state per cell). Optional third argument: the
//! largest number of outer iterations per run (default 2 000).

use foehncraft_flow::case::FlowCase;
use foehncraft_flow::residual::steady_flow;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

fn history(case: &Path, sets: &[String], max: usize, out: &Path) -> anyhow::Result<usize> {
    let case = FlowCase::from_file_with(case, sets)?;
    let mut fl = steady_flow(&case, Some(1))?;
    let mut csv = String::from("iteration,momentum,continuity,mass_flux\n");
    let mut n = 0;
    for it in 1..=max {
        fl.iterate()?;
        let r = fl.steady_residual()?.norms();
        writeln!(csv, "{it},{:.6e},{:.6e},{:.6e}", r[0], r[1], r[2])?;
        n = it;
        if r.iter().all(|&v| v < 1e-11) {
            break;
        }
    }
    std::fs::write(out, csv)?;
    Ok(n)
}

fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let (case_path, out) = match a.as_slice() {
        [c, o, ..] => (PathBuf::from(c), PathBuf::from(o)),
        _ => anyhow::bail!("usage: adjoint_a0 <case.toml> <output dir> [max iterations]"),
    };
    let max: usize = a.get(2).map(|s| s.parse()).transpose()?.unwrap_or(2000);
    std::fs::create_dir_all(&out)?;
    let ns = history(&case_path, &[], max, &out.join("history_simplec.csv"))?;
    let nc = history(
        &case_path,
        &["solver.algorithm=\"coupled\"".to_string()],
        max,
        &out.join("history_coupled.csv"),
    )?;
    println!("SIMPLEC: {ns} outer iterations, coupled: {nc}");

    // the cells: start residual, converged state and its residual
    let case = FlowCase::from_file(&case_path)?;
    let mut fl = steady_flow(&case, Some(1))?;
    fl.iterate()?;
    let r0 = fl.steady_residual()?;
    for _ in 1..ns {
        fl.iterate()?;
    }
    let r1 = fl.steady_residual()?;
    let st = fl.state()?;
    let c = &fl.disc.cell_centroid;
    let mag = |r: &foehncraft_flow::residual::SteadyResidual, i: usize| {
        (r.momentum[0][i].powi(2) + r.momentum[1][i].powi(2) + r.momentum[2][i].powi(2)).sqrt()
    };
    let mut csv = String::from("x,y,ux,uy,r_start,r_converged\n");
    for i in 0..fl.disc.n_cells {
        writeln!(
            csv,
            "{:.9e},{:.9e},{:.9e},{:.9e},{:.6e},{:.6e}",
            c.x[i],
            c.y[i],
            st.u.x[i],
            st.u.y[i],
            mag(&r0, i),
            mag(&r1, i)
        )?;
    }
    std::fs::write(out.join("cells.csv"), csv)?;
    println!("norms after 1 iteration {:?}", r0.norms());
    println!("norms at the end        {:?}", r1.norms());
    Ok(())
}
