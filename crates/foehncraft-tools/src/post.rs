//! `post.derive`: derived fields from a restart file (D72).

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::diag;
use foehncraft_flow::engine;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

fn latest() -> PathBuf {
    PathBuf::from("latest")
}

/// Input of `post.derive`, the options of `foehn post derive`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
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

    /// The output directory of the run, where the fields are written.
    fn effects(a: &DeriveInput, ctx: &mut Ctx<'_>) -> Result<crate::Effects> {
        let case = engine::load(&a.case, &a.set, &ctx.roots)?;
        let (_, dir) = foehncraft_flow::case::output_dir(&case, &a.case, a.out.clone());
        let writes = match &a.file {
            Some(f) => vec![f.clone()],
            None => vec![dir],
        };
        Ok(crate::Effects {
            writes,
            ..Default::default()
        })
    }

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

/// Input of `datum.report`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DatumReportInput {
    /// Flow case file (TOML) with the `[[datum]]` tables.
    pub case: PathBuf,
    /// Restart file, or `latest` for the newest one of the output directory.
    #[serde(default = "latest")]
    pub restart: PathBuf,
    /// Only these datums (default: all of the case).
    #[serde(default)]
    pub datums: Vec<String>,
    /// Case keys to override, `key=value` (D62); datums can be added this
    /// way only in the case file.
    #[serde(default)]
    pub set: Vec<String>,
    /// Output directory of the run (overrides `output.directory`), where
    /// `latest` looks.
    #[serde(default)]
    pub out: Option<PathBuf>,
}

/// `datum.report`.
pub struct DatumReport;

impl Tool for DatumReport {
    const ID: &'static str = "datum.report";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Gas flow through the datums of a case (mass, volume and momentum flow, area and mass-flow averages, cut error) from a restart file, without a rerun (D262, D267).";
    const CATEGORY: Category = Category::Post;
    const ACTION: Action = Action::Read;
    type Input = DatumReportInput;
    type Output = serde_json::Value;

    fn call(a: DatumReportInput, ctx: &mut Ctx<'_>) -> Result<serde_json::Value> {
        let roots = ctx.roots.clone();
        let case = engine::load(&a.case, &a.set, &roots)?;
        case.validate()?;
        let restart = restart_of(&case, &a.case, &a.restart, a.out.clone())?;
        roots.check_read(&restart, "the restart file")?;
        let rep = foehncraft_flow::derived::datum_report(&case, &restart, &a.datums)?;
        for d in &rep.datums {
            ctx.line(&format!(
                "{} ({}{}):",
                d.name,
                d.kind,
                if d.exact { ", exact" } else { "" }
            ));
            for (k, v) in &d.values {
                ctx.line(&format!("  {k:<24} {v:.6e}"));
            }
        }
        Ok(serde_json::to_value(&rep)?)
    }
}

fn restart_of(
    case: &foehncraft_flow::case::FlowCase,
    file: &Path,
    restart: &Path,
    out: Option<PathBuf>,
) -> Result<PathBuf> {
    let (_, dir) = foehncraft_flow::case::output_dir(case, file, out);
    if restart == Path::new("latest") {
        engine::latest_restart(&dir)?.ok_or_else(|| {
            diag::err(diag::Diagnostic::new(
                diag::RUN_NO_RESTART,
                format!(
                    "--restart latest: no restart file (restart.h5, restart_<N>.h5) in {}",
                    dir.display()
                ),
            ))
        })
    } else {
        Ok(restart.to_path_buf())
    }
}

/// Input of `spray.statistics` and `spray.penetration`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SprayFileInput {
    /// Flow case file (TOML) of the run: the liquid, the injectors and the
    /// penetration settings.
    pub case: PathBuf,
    /// Parcel file (`spray/parcels_<step>.vtkhdf`).
    pub file: PathBuf,
    /// Histogram bins (default 40).
    #[serde(default)]
    pub bins: Option<usize>,
    /// Smallest diameter of the histogram (m, default 1e-6).
    #[serde(default)]
    pub min: Option<f64>,
    /// Largest diameter of the histogram (m, default 1e-3).
    #[serde(default)]
    pub max: Option<f64>,
    /// Linear bins instead of logarithmic ones.
    #[serde(default)]
    pub linear: bool,
    /// Only the parcels inside this box or cylinder datum of the case.
    #[serde(default)]
    pub datum: Option<String>,
    /// Case keys to override, `key=value` (D62).
    #[serde(default)]
    pub set: Vec<String>,
}

fn spray_file(
    a: &SprayFileInput,
    ctx: &mut Ctx<'_>,
) -> Result<foehncraft_flow::spray::FileStatistics> {
    let roots = ctx.roots.clone();
    let case = engine::load(&a.case, &a.set, &roots)?;
    roots.check_read(&a.file, "the parcel file")?;
    let mut spec = foehncraft_flow::spray::HistSpec::default();
    if let Some(b) = a.bins {
        spec.bins = b;
    }
    if let Some(x) = a.min {
        spec.min = x;
    }
    if let Some(x) = a.max {
        spec.max = x;
    }
    spec.log = !a.linear;
    anyhow::ensure!(
        spec.bins > 0 && spec.max > spec.min && (!spec.log || spec.min > 0.0),
        "histogram: bins > 0 and 0 < min < max"
    );
    foehncraft_flow::spray::file_statistics(&case, &a.file, spec, a.datum.as_deref())
}

/// `spray.statistics`.
pub struct SprayStatistics;

impl Tool for SprayStatistics {
    const ID: &'static str = "spray.statistics";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Droplet statistics of a parcel file: D10, D20, D30, D32, D43, Dv10, Dv50, Dv90, span, velocities, histogram and Rosin-Rammler fit, from the same moment sums as the run (D264, D267).";
    const CATEGORY: Category = Category::Post;
    const ACTION: Action = Action::Read;
    type Input = SprayFileInput;
    type Output = serde_json::Value;

    fn call(a: SprayFileInput, ctx: &mut Ctx<'_>) -> Result<serde_json::Value> {
        let r = spray_file(&a, ctx)?;
        let s = &r.statistics;
        ctx.line(&format!(
            "{}: {} parcels, {:.4e} droplets, {:.4e} kg",
            r.file.display(),
            s.parcels,
            s.droplets,
            s.mass
        ));
        ctx.line(&format!(
            "  D10 {:.4e}  D32 {:.4e}  D43 {:.4e}  Dv10 {:.4e}  Dv50 {:.4e}  Dv90 {:.4e} m  span {:.3}",
            s.d10, s.d32, s.d43, s.dv10, s.dv50, s.dv90, s.span
        ));
        Ok(serde_json::to_value(&r)?)
    }
}

/// `spray.penetration`.
pub struct SprayPenetration;

impl Tool for SprayPenetration {
    const ID: &'static str = "spray.penetration";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Liquid penetration per injector of a parcel file: the distance within which the share `penetration_fraction` of the liquid lies, and the farthest parcel (D265, D267).";
    const CATEGORY: Category = Category::Post;
    const ACTION: Action = Action::Read;
    type Input = SprayFileInput;
    type Output = serde_json::Value;

    fn call(a: SprayFileInput, ctx: &mut Ctx<'_>) -> Result<serde_json::Value> {
        let r = spray_file(&a, ctx)?;
        for i in &r.injectors {
            ctx.line(&format!(
                "{}: {} parcels, {:.4e} kg, penetration {:.4e} m, farthest {:.4e} m",
                i.name, i.parcels, i.mass, i.penetration, i.farthest
            ));
        }
        Ok(serde_json::json!({
            "file": r.file,
            "time": r.time,
            "step": r.step,
            "injectors": r.injectors,
            "definitions": r.definitions,
        }))
    }
}
