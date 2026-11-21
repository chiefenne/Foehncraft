//! `post.derive`: derived fields from a restart file (D72).

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::diag;
use foehncraft_flow::engine;
use schemars::JsonSchema;
use serde::Deserialize;
use std::path::{Path, PathBuf};

fn latest() -> PathBuf {
    PathBuf::from("latest")
}

/// Input of `post.derive`, the options of `foehn post derive`.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DeriveInput {
    /// Flow case file (TOML) of the run.
    pub case: PathBuf,
    /// Restart file, or `latest` for the newest one of the output directory.
    #[serde(default = "latest")]
    pub restart: PathBuf,
    /// Derived fields (default: `output.derived` of the case, or every
    /// field when that is empty).
    #[serde(default)]
    pub fields: Option<Vec<String>>,
    /// File to write (default: `derived.vtkhdf` for `restart.h5`,
    /// `derived_<N>.vtkhdf` for `restart_<N>.h5`, next to the restart file).
    #[serde(default)]
    pub file: Option<PathBuf>,
    /// Case keys to override, `key=value` (D62).
    #[serde(default)]
    pub set: Vec<String>,
    /// Output directory of the run (overrides `output.directory`), where
    /// `latest` looks.
    #[serde(default)]
    pub out: Option<PathBuf>,
}

/// `post.derive`.
pub struct Derive;

impl Tool for Derive {
    const ID: &'static str = "post.derive";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Write the fields with derived fields (vorticity, Q, lambda2, wall shear stress, y+ ...) from a restart file of a case (D72).";
    const CATEGORY: Category = Category::Post;
    const ACTION: Action = Action::WriteFiles;
    type Input = DeriveInput;
    type Output = serde_json::Value;

    fn call(a: DeriveInput, ctx: &mut Ctx<'_>) -> Result<serde_json::Value> {
        use foehncraft_flow::case::output_dir;
        use foehncraft_flow::derived::{self, Derived};
        let roots = ctx.roots.clone();
        let case = engine::load(&a.case, &a.set, &roots)?;
        case.validate()?;
        let (_, dir) = output_dir(&case, &a.case, a.out.clone());
        let restart = if a.restart == Path::new("latest") {
            engine::latest_restart(&dir)?.ok_or_else(|| {
                diag::err(diag::Diagnostic::new(
                    diag::RUN_NO_RESTART,
                    format!(
                        "--restart latest: no restart file (restart.h5, restart_<N>.h5) in {}",
                        dir.display()
                    ),
                ))
            })?
        } else {
            a.restart.clone()
        };
        roots.check_read(&restart, "the restart file")?;
        let list: Vec<Derived> = match &a.fields {
            Some(names) => names
                .iter()
                .map(|n| {
                    Derived::parse(n.trim()).ok_or_else(|| {
                        let all: Vec<&str> = Derived::ALL.iter().map(|d| d.name()).collect();
                        diag::err(diag::Diagnostic::new(
                            diag::CASE_BAD_VALUE,
                            format!("unknown derived field `{n}`; known: {}", all.join(", ")),
                        ))
                    })
                })
                .collect::<Result<_>>()?,
            None if !case.output.derived.is_empty() => case.output.derived.clone(),
            // every field; the pressure coefficient only with a reference
            None => Derived::ALL
                .into_iter()
                .filter(|d| {
                    *d != Derived::PressureCoefficient
                        || derived::validate_list(&case, &[*d]).is_ok()
                })
                .collect(),
        };
        let file = a.file.clone().unwrap_or_else(|| {
            let stem = restart
                .file_stem()
                .map_or("restart".into(), |s| s.to_string_lossy().into_owned());
            let name = match stem.strip_prefix("restart") {
                Some(rest) => format!("derived{rest}.vtkhdf"),
                None => format!("{stem}_derived.vtkhdf"),
            };
            restart.with_file_name(name)
        });
        if let Some(d) = file.parent().filter(|d| !d.as_os_str().is_empty()) {
            roots.check_write(d)?;
        }
        let rep = derived::from_restart(&case, &restart, &list, &file)?;
        let at = match (rep.step, rep.time) {
            (Some(n), Some(t)) => format!("step {n}, t = {t} s"),
            _ => format!("iteration {}", rep.iteration),
        };
        ctx.line(&format!(
            "{}: {} cells at {at} from {}",
            rep.file.display(),
            rep.cells,
            rep.restart.display()
        ));
        ctx.line(&format!(
            "  derived: {}",
            rep.derived
                .iter()
                .map(|d| d.name())
                .collect::<Vec<_>>()
                .join(", ")
        ));
        Ok(serde_json::to_value(&rep)?)
    }
}
