//! `materials.*`: the material database (D94).

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::engine;
use foehncraft_flow::materials as m;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::PathBuf;

fn default_temperature() -> f64 {
    m::DEFAULT_TEMPERATURE
}

fn default_pressure() -> f64 {
    m::DEFAULT_PRESSURE
}

/// Input of `materials.list`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ListInput {
    /// Temperature (K) at which the properties are evaluated.
    #[serde(default = "default_temperature")]
    pub temperature: f64,
    /// Operating pressure (Pa) of the ideal-gas densities.
    #[serde(default = "default_pressure")]
    pub pressure: f64,
    /// User material database files, searched before the built-in database
    /// like `[materials] include` of a case; the files of `FOEHN_MATERIALS`
    /// are read as well.
    #[serde(default)]
    pub include: Vec<PathBuf>,
}

impl Default for ListInput {
    fn default() -> Self {
        Self {
            temperature: m::DEFAULT_TEMPERATURE,
            pressure: m::DEFAULT_PRESSURE,
            include: Vec::new(),
        }
    }
}

/// Input of `materials.show`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ShowInput {
    /// The material.
    pub name: String,
    /// Temperature (K) at which the properties are evaluated.
    #[serde(default = "default_temperature")]
    pub temperature: f64,
    /// Operating pressure (Pa) of the ideal-gas densities.
    #[serde(default = "default_pressure")]
    pub pressure: f64,
    /// User material database files, as for `materials.list`.
    #[serde(default)]
    pub include: Vec<PathBuf>,
}

/// Materials evaluated at a temperature and a pressure.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct MaterialsOutput {
    /// `foehncraft.materials`.
    pub format: &'static str,
    /// Format version of the engine.
    pub format_version: &'static str,
    /// SHA-256 of the built-in database.
    pub database_sha256: String,
    /// Each material: kind, description, the properties with their form
    /// and origin, the data range, sources and warnings.
    pub materials: Vec<Value>,
}

fn evaluate(
    name: Option<&str>,
    temperature: f64,
    pressure: f64,
    include: &[PathBuf],
    ctx: &mut Ctx<'_>,
) -> Result<MaterialsOutput> {
    anyhow::ensure!(
        temperature > 0.0 && pressure > 0.0,
        "--temperature and --pressure must be > 0"
    );
    for f in include {
        ctx.roots.check_read(f, "the material database file")?;
    }
    let cond = m::Conditions {
        temperature,
        pressure,
    };
    let mut lib = m::CaseMaterials {
        include: include.to_vec(),
        ..Default::default()
    };
    lib.load()?;
    let names: Vec<String> = match name {
        Some(n) => vec![n.to_string()],
        None => lib.names(),
    };
    let mut list = Vec::new();
    for n in &names {
        list.push(m::describe_in(n, &lib, cond)?);
    }
    text(name.is_none(), temperature, pressure, &list, ctx);
    Ok(MaterialsOutput {
        format: "foehncraft.materials",
        format_version: engine::FORMAT_VERSION,
        database_sha256: m::database_hash(),
        materials: list,
    })
}

fn text(all: bool, temperature: f64, pressure: f64, list: &[Value], ctx: &mut Ctx<'_>) {
    let num = |d: &Value, k: &str| -> Option<f64> {
        d[k]["value"]
            .as_f64()
            .or_else(|| d[k]["tensor"][0][0].as_f64())
    };
    let cell = |v: Option<f64>, w: usize, p: usize| match v {
        Some(x) => format!("{x:>w$.p$e}"),
        None => format!("{:>w$}", "–"),
    };
    if all {
        ctx.line(&format!(
            "material database ({} materials) at T = {} K, p = {} Pa",
            list.len(),
            temperature,
            pressure
        ));
        ctx.line(&format!(
            "  {:<18} {:<7} {:>11} {:>11} {:>11} {:>11}  {:<16} description",
            "name", "kind", "ρ (kg/m³)", "μ (Pa s)", "c_p (J/kgK)", "k (W/mK)", "range (K)"
        ));
        for d in list {
            let range = d["temperature_range"]
                .as_array()
                .map(|r| format!("{}–{}", r[0], r[1]))
                .unwrap_or_default();
            ctx.line(&format!(
                "  {:<18} {:<7} {} {} {} {}  {:<16} {}",
                d["name"].as_str().unwrap_or(""),
                d["kind"].as_str().unwrap_or(""),
                cell(num(d, "density"), 11, 4),
                cell(num(d, "viscosity"), 11, 3),
                cell(num(d, "specific_heat"), 11, 4),
                cell(num(d, "conductivity"), 11, 3),
                range,
                d["description"].as_str().unwrap_or("")
            ));
        }
    } else if let Some(d) = list.first() {
        ctx.line(&format!(
            "{}: {}",
            d["name"].as_str().unwrap_or(""),
            d["description"].as_str().unwrap_or("")
        ));
        for what in ["density", "viscosity", "specific_heat", "conductivity"] {
            if d[what].is_null() {
                continue;
            }
            let v = match num(d, what) {
                Some(x) => format!("{x}"),
                None => "–".into(),
            };
            ctx.line(&format!(
                "  {what:<14} {v} ({})",
                d[what]["origin"].as_str().unwrap_or("")
            ));
        }
        for w in d["warnings"].as_array().into_iter().flatten() {
            ctx.line(&format!("  warning: {}", w.as_str().unwrap_or("")));
        }
    }
}

/// `materials.list`.
pub struct List;

impl Tool for List {
    const ID: &'static str = "materials.list";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str =
        "Every material of the database with its properties at a temperature and a pressure.";
    const CATEGORY: Category = Category::Materials;
    const ACTION: Action = Action::Read;
    type Input = ListInput;
    type Output = MaterialsOutput;

    fn call(a: ListInput, ctx: &mut Ctx<'_>) -> Result<MaterialsOutput> {
        evaluate(None, a.temperature, a.pressure, &a.include, ctx)
    }
}

/// `materials.show`.
pub struct Show;

impl Tool for Show {
    const ID: &'static str = "materials.show";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str =
        "One material with its property forms, values, data range and sources.";
    const CATEGORY: Category = Category::Materials;
    const ACTION: Action = Action::Read;
    type Input = ShowInput;
    type Output = MaterialsOutput;

    fn call(a: ShowInput, ctx: &mut Ctx<'_>) -> Result<MaterialsOutput> {
        evaluate(Some(&a.name), a.temperature, a.pressure, &a.include, ctx)
    }
}
