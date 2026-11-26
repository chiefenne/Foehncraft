//! `circuits.run`: the one-dimensional water and steam tube circuits of a
//! boiler (steam furnace stage SF7b; crate `foehncraft-circuits`).

use crate::{Action, Category, Ctx, Effects, Tool};
use anyhow::Result;
use foehncraft_circuits::case::Case;
use foehncraft_circuits::output;
use foehncraft_flow::diag::{self, Diagnostic};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// The diagnostic of an error of the circuits (its code is one of
/// `CIR001` to `CIR004`, registered in `foehncraft_flow::diag`).
fn diagnostic(e: foehncraft_circuits::Error) -> anyhow::Error {
    let code = diag::info(e.code).map_or(diag::GENERAL, |c| c.code);
    diag::err(Diagnostic::new(code, e.message))
}

fn read_case(path: &Path, ctx: &mut Ctx<'_>) -> Result<Case> {
    ctx.roots.check_read(path, "the circuit case")?;
    let case = Case::from_file(path).map_err(diagnostic)?;
    if let Some(t) = case.chf_table_path() {
        ctx.roots.check_read(&t, "the CHF table")?;
    }
    Ok(case)
}

/// Input of `circuits.run`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunInput {
    /// The circuit case file (TOML).
    pub case: PathBuf,
    /// Output directory (default `output/<case name>`).
    #[serde(default)]
    pub out: Option<PathBuf>,
}

impl RunInput {
    fn name(&self, case: &Case) -> String {
        case.circuit.name.clone().unwrap_or_else(|| {
            self.case
                .file_stem()
                .map_or("circuit".into(), |s| s.to_string_lossy().into_owned())
        })
    }

    fn dir(&self, case: &Case) -> PathBuf {
        self.out
            .clone()
            .unwrap_or_else(|| PathBuf::from("output").join(self.name(case)))
    }
}

/// `circuits.run`.
pub struct Run;

impl Tool for Run {
    const ID: &'static str = "circuits.run";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Solve a water and steam tube circuit of a boiler (tubes, headers, pressure nodes, steam drums with natural circulation; two-phase pressure drop, flow boiling wall temperatures, critical heat flux margin) and write `result.json` and `profiles.csv`; returns the result.";
    const CATEGORY: Category = Category::Circuits;
    const ACTION: Action = Action::StartRun;
    type Input = RunInput;
    type Output = serde_json::Value;

    fn call(a: RunInput, ctx: &mut Ctx<'_>) -> Result<serde_json::Value> {
        let c = read_case(&a.case, ctx)?;
        let name = a.name(&c);
        let dir = a.dir(&c);
        ctx.roots.check_write(&dir)?;
        let sol = output::run(&c).map_err(diagnostic)?;
        let json = sol.to_json(&name, &c);
        std::fs::create_dir_all(&dir)?;
        std::fs::write(
            dir.join("result.json"),
            serde_json::to_string_pretty(&json)? + "\n",
        )?;
        if c.output.profiles.unwrap_or(true) {
            std::fs::write(dir.join("profiles.csv"), sol.to_csv())?;
        }
        ctx.line(&format!(
            "{name}: converged in {} Newton iterations",
            sol.history.len() - 1
        ));
        for n in &sol.model.nodes {
            if let Some(d) = sol.drum(&n.name) {
                ctx.line(&format!(
                    "  drum {}: steam {:.6} kg/s, riser flow {:.6} kg/s, circulation ratio {:.3}",
                    n.name, d.steam, d.riser_flow, d.circulation_ratio
                ));
            }
        }
        for (t, m) in sol.model.tubes.iter().zip(&sol.marches) {
            let chf =
                output::min_chf(m).map_or(String::new(), |(_, r)| format!(", CHF ratio {r:.3}"));
            ctx.line(&format!(
                "  tube {}: {:.6e} kg/s each, dp {:.6e} Pa{chf}",
                t.name,
                sol.x[t.im],
                m.dp_fric + m.dp_grav + m.dp_acc
            ));
        }
        for w in &sol.warnings {
            ctx.line(&format!("{}: {}", w.code, w.message));
        }
        ctx.line(&format!("results in {}", dir.display()));
        Ok(json)
    }

    fn effects(a: &RunInput, ctx: &mut Ctx<'_>) -> Result<Effects> {
        let c = read_case(&a.case, ctx)?;
        Ok(Effects {
            writes: vec![a.dir(&c)],
            ..Default::default()
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::{Ctx, call};
    use foehncraft_flow::diag;
    use foehncraft_flow::engine::BuildInfo;
    use serde_json::json;
    use std::path::PathBuf;

    #[test]
    fn circuit_codes_are_registered() {
        use foehncraft_circuits::codes;
        for c in [codes::CASE, codes::CONVERGENCE, codes::RANGE, codes::CHF] {
            assert!(
                diag::info(c).is_some(),
                "{c} is not in foehncraft_flow::diag"
            );
        }
    }

    #[test]
    fn the_drum_loop_runs_through_the_tool() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
        let tmp = std::env::temp_dir().join(format!("foehn-circuits-{}", std::process::id()));
        let mut ctx = Ctx::new(BuildInfo::default());
        let case = root.join("examples/cases/circuits/drum_natural_circulation.toml");
        let r = call("circuits.run", json!({"case": case, "out": tmp}), &mut ctx);
        assert_eq!(r.status, "finished", "{:?}", r.diagnostics);
        assert_eq!(r.result["format"], "foehncraft.circuit-result");
        assert!(tmp.join("profiles.csv").is_file());
        let cr = r.result["nodes"][0]["drum"]["circulation_ratio"]
            .as_f64()
            .unwrap();
        assert!(cr > 3.0, "{cr}");
    }
}
