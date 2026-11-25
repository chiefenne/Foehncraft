//! `structure.run` and `structure.fitclamp`: the linear structural solver
//! (L53, MS5; crate `foehncraft-structure`).

use crate::{Action, Category, Ctx, Effects, Tool};
use anyhow::Result;
use foehncraft_flow::diag::{self, Diagnostic};
use foehncraft_structure::case::{self, Case};
use foehncraft_structure::clamp::{self, Which};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// The diagnostic of an error of the structural solver (its code is one of
/// `STR001` to `STR004`, registered in `foehncraft_flow::diag`).
fn diagnostic(e: foehncraft_structure::Error) -> anyhow::Error {
    let code = diag::info(e.code).map_or(diag::GENERAL, |c| c.code);
    diag::err(Diagnostic::new(code, e.message))
}

fn read_case(path: &Path, ctx: &mut Ctx<'_>) -> Result<Case> {
    ctx.roots.check_read(path, "the structural case")?;
    let text = std::fs::read_to_string(path).map_err(|e| {
        diag::err(Diagnostic::new(
            diag::CASE_FILE,
            format!("{}: {e}", path.display()),
        ))
    })?;
    let mut case = Case::from_toml(&text).map_err(diagnostic)?;
    // mesh files of solids are read relative to the case file
    case.base = Some(path.parent().map_or_else(PathBuf::new, Path::to_path_buf));
    for mesh in case.mesh_paths() {
        ctx.roots.check_read(&mesh, "the mesh of a solid")?;
    }
    Ok(case)
}

/// Input of `structure.run`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunInput {
    /// The structural case file (TOML).
    pub case: PathBuf,
    /// Output directory (default `output/<case name>`).
    #[serde(default)]
    pub out: Option<PathBuf>,
}

impl RunInput {
    fn name(&self, case: &Case) -> String {
        case.structure.name.clone().unwrap_or_else(|| {
            self.case
                .file_stem()
                .map_or("structure".into(), |s| s.to_string_lossy().into_owned())
        })
    }

    fn dir(&self, case: &Case) -> PathBuf {
        self.out
            .clone()
            .unwrap_or_else(|| PathBuf::from("output").join(self.name(case)))
    }
}

/// `structure.run`.
pub struct Run;

impl Tool for Run {
    const ID: &'static str = "structure.run";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Run a structural case (beams, shells, quadratic solids from Gmsh or a block, clamp springs, temperatures; static, modal and transient analyses) and write `result.json`, `structure.vtk` and, for a transient, `transient.csv`; returns the result.";
    const CATEGORY: Category = Category::Structure;
    const ACTION: Action = Action::StartRun;
    type Input = RunInput;
    type Output = serde_json::Value;

    fn call(a: RunInput, ctx: &mut Ctx<'_>) -> Result<serde_json::Value> {
        let c = read_case(&a.case, ctx)?;
        let name = a.name(&c);
        let dir = a.dir(&c);
        ctx.roots.check_write(&dir)?;
        let out = case::run(&c).map_err(diagnostic)?;
        let json = out.to_json(&name);
        std::fs::create_dir_all(&dir)?;
        std::fs::write(
            dir.join("result.json"),
            serde_json::to_string_pretty(&json)? + "\n",
        )?;
        if let Some(csv) = out.to_csv() {
            std::fs::write(dir.join("transient.csv"), csv)?;
        }
        if c.output.vtk {
            std::fs::write(dir.join("structure.vtk"), out.to_vtk(&name))?;
        }
        if let Some(m) = &out.modal {
            ctx.line(&format!(
                "{name}: {} modes, Sturm check {}",
                m.modes.len(),
                if m.sturm.ok() { "passed" } else { "FAILED" }
            ));
            for (i, mo) in m.modes.iter().enumerate() {
                ctx.line(&format!("  mode {:>3}  {:>14.6} Hz", i + 1, mo.frequency()));
            }
        }
        if out.statics.is_some() {
            ctx.line(&format!(
                "{name}: static, largest displacement {:.6e} m",
                json["static"]["max_displacement"]["value"]
                    .as_f64()
                    .unwrap_or(0.0)
            ));
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

/// Input of `structure.fitclamp`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FitInput {
    /// The structural case file (TOML) with the clamp.
    pub case: PathBuf,
    /// The clamp, by its place in `[[clamps]]` (from 0).
    #[serde(default)]
    pub clamp: usize,
    /// The spring to fit: `rotation` (about the clamp line) or
    /// `translation` (along the normal); the other keeps its case value.
    #[serde(default = "rotation")]
    pub spring: String,
    /// The mode whose frequency was measured (from 1).
    #[serde(default = "one")]
    pub mode: usize,
    /// The measured natural frequency (Hz).
    pub frequency: f64,
}

fn rotation() -> String {
    "rotation".into()
}
fn one() -> usize {
    1
}

/// `structure.fitclamp`.
pub struct FitClamp;

impl Tool for FitClamp {
    const ID: &'static str = "structure.fitclamp";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Fit the clamp spring per length (D290) so that a mode of the structural case has a measured natural frequency.";
    const CATEGORY: Category = Category::Structure;
    const ACTION: Action = Action::Read;
    type Input = FitInput;
    type Output = clamp::Fit;

    fn call(a: FitInput, ctx: &mut Ctx<'_>) -> Result<clamp::Fit> {
        let c = read_case(&a.case, ctx)?;
        let which = match a.spring.as_str() {
            "rotation" => Which::Rotation,
            "translation" => Which::Translation,
            s => {
                return Err(diag::err(Diagnostic::new(
                    diag::TOOL_INPUT,
                    format!("`spring` is `rotation` or `translation`, not `{s}`"),
                )));
            }
        };
        let f = clamp::fit(&c, a.clamp, which, a.mode, a.frequency).map_err(diagnostic)?;
        ctx.line(&format!(
            "clamp {} {}: k = {:.6e} per metre of clamp line ({:.6e} in total), mode {} at {:.6} Hz (rigid clamp {:.6} Hz), {} modal solves",
            a.clamp, a.spring, f.k_per_length, f.k_total, a.mode, f.frequency, f.frequency_rigid, f.evaluations
        ));
        Ok(f)
    }
}

#[cfg(test)]
mod tests {
    use crate::{Ctx, call};
    use foehncraft_flow::diag;
    use foehncraft_flow::engine::BuildInfo;
    use serde_json::json;
    use std::path::PathBuf;

    fn root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
    }

    #[test]
    fn structure_codes_are_registered() {
        use foehncraft_structure::codes;
        for c in [codes::CASE, codes::SINGULAR, codes::MODES, codes::FIT] {
            assert!(
                diag::info(c).is_some(),
                "{c} is not in foehncraft_flow::diag"
            );
        }
    }

    #[test]
    fn catalogue_cases_run_and_the_fit_matches() {
        let tmp = std::env::temp_dir().join(format!("foehn-structure-{}", std::process::id()));
        let mut ctx = Ctx::new(BuildInfo::default());
        for name in [
            "reed_beam",
            "reed_beam_clamp_springs",
            "reed_plate",
            "reed_plate_pulse",
            "reed_plate_hot",
        ] {
            let case = root().join(format!("examples/cases/structure/{name}.toml"));
            let out = tmp.join(name);
            let r = call("structure.run", json!({"case": case, "out": out}), &mut ctx);
            assert_eq!(r.status, "finished", "{name}: {:?}", r.diagnostics);
            let res = &r.result;
            if !res["modal"].is_null() {
                assert_eq!(res["modal"]["sturm"]["ok"], true, "{name}");
            }
            assert!(out.join("structure.vtk").is_file());
            assert_eq!(
                out.join("transient.csv").is_file(),
                name == "reed_plate_pulse"
            );
            let rigid = res["clamps"][0]["rotation"]["rigid"].as_bool().unwrap();
            assert_eq!(rigid, name == "reed_beam" || name == "reed_plate_hot");
        }
        // the fit of examples/cases/tools/structure_fit_clamp.json
        let case = root().join("examples/cases/structure/reed_beam_clamp_springs.toml");
        let r = call(
            "structure.fitclamp",
            json!({"case": case, "frequency": 150.0}),
            &mut ctx,
        );
        assert_eq!(r.status, "finished", "{:?}", r.diagnostics);
        let f = r.result["frequency"].as_f64().unwrap();
        assert!((f - 150.0).abs() < 1e-6, "{f}");
        // a clamp along a shell edge fits the same way
        let plate = root().join("examples/cases/structure/reed_plate.toml");
        let r = call(
            "structure.fitclamp",
            json!({"case": plate, "frequency": 150.0}),
            &mut ctx,
        );
        assert_eq!(r.status, "finished", "{:?}", r.diagnostics);
        let f = r.result["frequency"].as_f64().unwrap();
        assert!((f - 150.0).abs() < 1e-6, "{f}");
        // above the rigid clamp: STR004
        let r = call(
            "structure.fitclamp",
            json!({"case": case, "frequency": 400.0}),
            &mut ctx,
        );
        assert_eq!(r.diagnostics[0].code, "STR004");
        let _ = std::fs::remove_dir_all(&tmp);
    }
}
