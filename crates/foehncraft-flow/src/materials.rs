//! Material properties (D94–D102): density and viscosity (stage 5), the
//! thermal data c_p and k with solids and user database files (stage 7,
//! PR 7.0c).
//!
//! Three layers (D94):
//!
//! 1. the **database**: named materials with their data, validity ranges
//!    and sources, built into the binary from `data/materials.toml`
//!    (air, water, an engine oil, ethylene glycol, glycerol), read-only;
//! 2. the **case materials** `[materials.<name>]`: new materials, or
//!    changes of a database material with `base = "<name>"`, and the
//!    assignment of materials to the fluid (`[fluid] material`) and to
//!    regions of cell zones (`[regions.<name>]`, D112);
//! 3. the **resolved property program** ([`PropertySet`]): every property
//!    of every region as one value or one compiled expression (D53), with
//!    its source, written into the resolved case and hashed into the run
//!    record (D85). The solver sees only this; the database is not read
//!    again during the run, and properties do not change during a run
//!    (D100).
//!
//! Forms (D95) evaluated in this version: a number, an expression (of x,
//! y, z for the viscosity; `T` is the case temperature), `constant`,
//! `polynomial` in T, `sutherland`, `power-law`, `table` in T (linear or
//! the monotone cubic of Steffen 1990, clamped with a warning outside the
//! table, D96) and `ideal-gas` for the density. The flow is isothermal
//! and incompressible until the energy equation (L1): every form is
//! evaluated once at `[fluid] temperature` (default 293.15 K) and the
//! operating pressure `[fluid] operating_pressure` (default 101 325 Pa),
//! the density must be one number for the whole fluid (a varying density
//! comes with L1), and only the viscosity may vary in space.
//!
//! Thermal data (stage 7, D163): every material may carry a specific heat
//! c_p (any form above, or NASA-7/9 polynomials with the molar mass) and a
//! thermal conductivity k (any form, or orthotropic: three principal values
//! along two given axes). Solids need ρ, c_p and k (MAT006). They are
//! evaluated by [`thermal_properties`] at any temperature; the energy
//! equation reads them from PR 7a1 on. User database files
//! (`[materials] include = [..]`, and the files named by the environment
//! variable `FOEHN_MATERIALS`) are searched after the case's own materials
//! and before the built-in database (D97).
//!
//! `[fluid] density` and `viscosity` stay valid as the shorthand of one
//! constant material on the whole mesh (D100), so every older case keeps
//! running with the same numbers.

use crate::diag::{self, Diagnostic};
use crate::expr::{self, Env, Scalar, Value, Vars};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::sync::{Arc, OnceLock};

/// Universal gas constant (J/(mol K), CODATA 2018, exact).
pub const R_GAS: f64 = 8.314_462_618;
/// The default case temperature (K).
pub const DEFAULT_TEMPERATURE: f64 = 293.15;
/// The default operating pressure (Pa).
pub const DEFAULT_PRESSURE: f64 = 101_325.0;

/// Interpolation of a property table (D96).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Interpolation {
    /// Piecewise linear.
    #[default]
    Linear,
    /// Monotone piecewise cubic (Steffen 1990): no overshoot between the
    /// points, continuous first derivative.
    MonotoneCubic,
}

/// A property given by a form (D95).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(tag = "form", rename_all = "kebab-case", deny_unknown_fields)]
pub enum Form {
    /// A constant.
    Constant {
        /// The value (SI).
        value: f64,
        /// Where the value comes from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// Σ aᵢ Tⁱ (Horner form), T in K.
    Polynomial {
        /// a₀, a₁, … (SI).
        coefficients: Vec<f64>,
        /// Validity range in T (K).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        range: Option<[f64; 2]>,
        /// Where the coefficients come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// μ = μ_ref (T/T_ref)^1.5 (T_ref + S)/(T + S) (Sutherland 1893), for
    /// the viscosity of a gas.
    Sutherland {
        /// μ_ref (Pa s).
        mu_ref: f64,
        /// T_ref (K).
        t_ref: f64,
        /// S (K).
        s: f64,
        /// Validity range in T (K).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        range: Option<[f64; 2]>,
        /// Where the coefficients come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// value_ref (T/T_ref)^exponent.
    PowerLaw {
        /// The value at T_ref (SI).
        value_ref: f64,
        /// T_ref (K).
        t_ref: f64,
        /// The exponent.
        exponent: f64,
        /// Validity range in T (K).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        range: Option<[f64; 2]>,
        /// Where the coefficients come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// A table in T (K), clamped with a warning outside (D96).
    Table {
        /// Temperatures (K), strictly increasing.
        #[serde(default)]
        temperature: Vec<f64>,
        /// Values (SI).
        #[serde(default)]
        values: Vec<f64>,
        /// A file of two columns (T, value) instead of the arrays.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        file: Option<PathBuf>,
        /// Interpolation.
        #[serde(default)]
        interpolation: Interpolation,
        /// Where the data come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// An expression (D53); `T` is the case temperature.
    Expression {
        /// The expression.
        expression: String,
        /// Where it comes from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// ρ = p M / (R T) at the operating pressure (for the density).
    IdealGas {
        /// Molar mass (kg/mol); default the material's `molar_mass`.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        molar_mass: Option<f64>,
        /// Where it comes from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// NASA 7-coefficient polynomials of the specific heat (the CHEMKIN
    /// form; McBride, Zehe & Gordon 2002): c_p/R = a₁ + a₂T + a₃T² + a₄T³ +
    /// a₅T⁴ per temperature range, a₆ and a₇ (enthalpy and entropy
    /// constants) kept; c_p = (R/M)(c_p/R) with the material's molar mass.
    /// Only for `specific_heat`.
    Nasa7 {
        /// The range limits (K), increasing: n + 1 values for n ranges
        /// (the CHEMKIN record: T_low, T_mid, T_high).
        #[serde(default)]
        temperature: Vec<f64>,
        /// a₁…a₇ per range, from the lowest range up.
        #[serde(default)]
        coefficients: Vec<Vec<f64>>,
        /// A CHEMKIN THERMO file to read the record of `species` from
        /// (relative to the file that names it), instead of the arrays.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        file: Option<PathBuf>,
        /// The species name in `file`.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        species: Option<String>,
        /// Where the coefficients come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
    /// NASA 9-coefficient polynomials of the specific heat (McBride, Zehe
    /// & Gordon 2002): c_p/R = a₁T⁻² + a₂T⁻¹ + a₃ + a₄T + a₅T² + a₆T³ + a₇T⁴
    /// per range (b₁, b₂ may follow and are kept). Only for
    /// `specific_heat`.
    Nasa9 {
        /// The range limits (K), increasing: n + 1 values for n ranges.
        #[serde(default)]
        temperature: Vec<f64>,
        /// a₁…a₇ (and b₁, b₂) per range, from the lowest range up.
        #[serde(default)]
        coefficients: Vec<Vec<f64>>,
        /// Where the coefficients come from.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        source: Option<String>,
    },
}

impl Form {
    /// The form's name.
    pub fn name(&self) -> &'static str {
        match self {
            Form::Constant { .. } => "constant",
            Form::Polynomial { .. } => "polynomial",
            Form::Sutherland { .. } => "sutherland",
            Form::PowerLaw { .. } => "power-law",
            Form::Table { .. } => "table",
            Form::Expression { .. } => "expression",
            Form::IdealGas { .. } => "ideal-gas",
            Form::Nasa7 { .. } => "nasa-7",
            Form::Nasa9 { .. } => "nasa-9",
        }
    }

    fn source(&self) -> Option<&str> {
        match self {
            Form::Constant { source, .. }
            | Form::Polynomial { source, .. }
            | Form::Sutherland { source, .. }
            | Form::PowerLaw { source, .. }
            | Form::Table { source, .. }
            | Form::Expression { source, .. }
            | Form::IdealGas { source, .. }
            | Form::Nasa7 { source, .. }
            | Form::Nasa9 { source, .. } => source.as_deref(),
        }
    }

    fn range(&self) -> Option<[f64; 2]> {
        match self {
            Form::Polynomial { range, .. }
            | Form::Sutherland { range, .. }
            | Form::PowerLaw { range, .. } => *range,
            Form::Table { temperature, .. }
            | Form::Nasa7 { temperature, .. }
            | Form::Nasa9 { temperature, .. }
                if temperature.len() >= 2 =>
            {
                Some([temperature[0], temperature[temperature.len() - 1]])
            }
            _ => None,
        }
    }
}

/// A property: a number, an expression string, or a form.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum Property {
    /// A number (SI).
    Number(f64),
    /// An expression (D53).
    Expr(String),
    /// A form of D95.
    Form(Form),
}

impl From<Scalar> for Property {
    fn from(s: Scalar) -> Self {
        match s {
            Scalar::Number(v) => Property::Number(v),
            Scalar::Expr(e) => Property::Expr(e),
        }
    }
}

/// The kind of a material.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum MaterialKind {
    /// A gas.
    Gas,
    /// A liquid.
    Liquid,
    /// A solid (with the energy equation, L1).
    Solid,
}

/// A thermal conductivity: one value (any property form) or orthotropic.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum Conductivity {
    /// The same k in every direction (W/(m K)).
    Isotropic(Property),
    /// Three principal values along two given axes (a circuit board, a
    /// composite, a laminated core).
    Orthotropic(Box<Orthotropic>),
}

/// The `form` of an orthotropic conductivity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum OrthotropicForm {
    /// Three principal values along two given axes.
    Orthotropic,
}

/// An orthotropic conductivity K = Σᵢ kᵢ eᵢ eᵢᵀ with e₁ = `axes[0]`, e₂ the
/// part of `axes[1]` normal to e₁ and e₃ = e₁ × e₂ (all normalised):
/// `{ form = "orthotropic", values = [k1, k2, k3], axes = [[..], [..]] }`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Orthotropic {
    /// Always `orthotropic`.
    pub form: OrthotropicForm,
    /// k₁, k₂, k₃ along the axes (W/(m K)), each a number, an expression
    /// or any property form.
    pub values: [Property; 3],
    /// The first two principal axes in the mesh frame (the second is made
    /// normal to the first; not parallel to it). Left out: MAT008.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub axes: Option<[[f64; 3]; 2]>,
    /// Where the values come from.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
}

/// A material, `[materials.<name>]` of the case or of the database.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct MaterialSection {
    /// A material (of the case or the database) this one starts from;
    /// the keys given here replace its keys.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub base: Option<String>,
    /// One line on what it is.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// Gas, liquid or solid.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub kind: Option<MaterialKind>,
    /// Temperature range of the data (K).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature_range: Option<[f64; 2]>,
    /// Molar mass (kg/mol).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub molar_mass: Option<f64>,
    /// Where the molar mass comes from.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub molar_mass_source: Option<String>,
    /// Density ρ (kg/m³).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub density: Option<Property>,
    /// Dynamic viscosity μ (Pa s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub viscosity: Option<Property>,
    /// Specific heat at constant pressure c_p (J/(kg K)); with the energy
    /// equation, and always for solids.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub specific_heat: Option<Property>,
    /// Thermal conductivity k (W/(m K)), isotropic or orthotropic; with
    /// the energy equation, and always for solids.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub conductivity: Option<Conductivity>,
}

/// `[materials]` of a case: the case's own materials by name and the user
/// database files to search after them (D97).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
pub struct CaseMaterials {
    /// User database files in the format of the built-in database
    /// (`format = "foehncraft.materials"`, `[materials.<name>]` tables),
    /// relative to the case file. They are searched in this order after
    /// the case's own materials, then the files of `FOEHN_MATERIALS`, then
    /// the built-in database; a material cannot be named `include`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub include: Vec<PathBuf>,
    /// The case's materials: new ones, or changes of database materials
    /// with `base`.
    #[serde(flatten)]
    pub named: BTreeMap<String, MaterialSection>,
    /// The materials read from the user database files, file by file (not
    /// a key; [`CaseMaterials::load`]).
    #[serde(skip)]
    #[schemars(skip)]
    pub files: Vec<(PathBuf, BTreeMap<String, MaterialSection>)>,
}

impl From<BTreeMap<String, MaterialSection>> for CaseMaterials {
    fn from(named: BTreeMap<String, MaterialSection>) -> Self {
        Self {
            named,
            ..Default::default()
        }
    }
}

impl CaseMaterials {
    /// True without own materials and include files.
    pub fn is_empty(&self) -> bool {
        self.include.is_empty() && self.named.is_empty()
    }

    /// Resolves the include paths against the case directory.
    pub fn resolve_paths(&mut self, dir: &std::path::Path) {
        for f in &mut self.include {
            *f = dir.join(&*f);
        }
    }

    /// Reads the include files and the files of `FOEHN_MATERIALS`, and the
    /// NASA and table data files that the materials name (MAT009 when a
    /// file does not read).
    pub fn load(&mut self) -> Result<()> {
        let mut paths = self.include.clone();
        if let Some(v) = std::env::var_os("FOEHN_MATERIALS") {
            paths.extend(std::env::split_paths(&v).filter(|p| !p.as_os_str().is_empty()));
        }
        self.files.clear();
        for p in paths {
            let text = std::fs::read_to_string(&p).map_err(|e| {
                mat_err(
                    diag::MAT_FILE,
                    "materials.include",
                    format!("material file {}: {e}", p.display()),
                )
            })?;
            let mut f: DatabaseFile = toml::from_str(&text).map_err(|e| {
                mat_err(
                    diag::MAT_FILE,
                    "materials.include",
                    format!("material file {}: {e}", p.display()),
                )
            })?;
            if f.format != "foehncraft.materials" {
                return Err(mat_err(
                    diag::MAT_FILE,
                    "materials.include",
                    format!(
                        "material file {}: format `{}`, expected `foehncraft.materials`",
                        p.display(),
                        f.format
                    ),
                ));
            }
            let dir = p.parent().unwrap_or(std::path::Path::new("")).to_path_buf();
            for (n, m) in &mut f.materials {
                load_data_files(m, &dir, &format!("material {n} ({})", p.display()))?;
            }
            self.files.push((p, f.materials));
        }
        for (n, m) in &mut self.named {
            load_data_files(m, std::path::Path::new(""), &format!("materials.{n}"))?;
        }
        Ok(())
    }

    /// The search layers: own materials, the files, the built-in database.
    fn layers(&self) -> Vec<(&BTreeMap<String, MaterialSection>, String)> {
        let mut v = vec![(&self.named, "case".to_string())];
        for (p, m) in &self.files {
            v.push((m, p.display().to_string()));
        }
        v.push((database(), "database".to_string()));
        v
    }

    /// Every material name of every layer (the case's first).
    pub fn names(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for (m, _) in self.layers() {
            for n in m.keys() {
                if !out.contains(n) {
                    out.push(n.clone());
                }
            }
        }
        out
    }
}

/// Reads the data files a material names (NASA-7 CHEMKIN records, table
/// files) into its forms; `dir` resolves relative paths.
fn load_data_files(m: &mut MaterialSection, dir: &std::path::Path, owner: &str) -> Result<()> {
    let mut forms: Vec<(&str, &mut Form)> = Vec::new();
    for (what, p) in [
        ("density", m.density.as_mut()),
        ("viscosity", m.viscosity.as_mut()),
        ("specific_heat", m.specific_heat.as_mut()),
    ] {
        if let Some(Property::Form(f)) = p {
            forms.push((what, f));
        }
    }
    match m.conductivity.as_mut() {
        Some(Conductivity::Isotropic(Property::Form(f))) => forms.push(("conductivity", f)),
        Some(Conductivity::Orthotropic(o)) => {
            for p in o.values.iter_mut() {
                if let Property::Form(f) = p {
                    forms.push(("conductivity", f));
                }
            }
        }
        _ => {}
    }
    for (what, f) in forms {
        let key = format!("{owner}.{what}");
        let read = |file: &PathBuf| -> Result<String> {
            let path = dir.join(file);
            std::fs::read_to_string(&path).map_err(|e| {
                mat_err(
                    diag::MAT_FILE,
                    &key,
                    format!("{key}: data file {}: {e}", path.display()),
                )
            })
        };
        match f {
            Form::Nasa7 {
                temperature,
                coefficients,
                file,
                species,
                source,
            } => {
                let Some(path) = file.take() else { continue };
                let text = read(&path)?;
                let sp = species.clone().ok_or_else(|| {
                    mat_err(
                        diag::MAT_FORM,
                        &key,
                        format!("{key}: a NASA-7 `file` needs `species`"),
                    )
                })?;
                let (t, c) = nasa7_from_chemkin(&text, &sp).map_err(|e| {
                    mat_err(
                        diag::MAT_FILE,
                        &key,
                        format!("{key}: {}: {e:#}", dir.join(&path).display()),
                    )
                })?;
                *temperature = t;
                *coefficients = c;
                if source.is_none() {
                    *source = Some(format!("NASA-7 record {sp} of {}", path.display()));
                }
            }
            Form::Table {
                temperature,
                values,
                file,
                ..
            } => {
                let Some(path) = file.take() else { continue };
                let text = read(&path)?;
                let (t, v) = two_columns(&text).map_err(|e| {
                    mat_err(
                        diag::MAT_FILE,
                        &key,
                        format!("{key}: {}: {e}", dir.join(&path).display()),
                    )
                })?;
                *temperature = t;
                *values = v;
            }
            _ => {}
        }
    }
    Ok(())
}

/// The rows `T value` of a two-column text file (`#` comments, blank
/// lines, commas or white space between the columns).
fn two_columns(text: &str) -> std::result::Result<(Vec<f64>, Vec<f64>), String> {
    let (mut t, mut v) = (Vec::new(), Vec::new());
    for (i, line) in text.lines().enumerate() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|c| !c.is_empty())
            .collect();
        let num = |c: &str| {
            c.parse::<f64>()
                .map_err(|_| format!("line {}: `{c}` is not a number", i + 1))
        };
        // a header line of names before the data
        if t.is_empty() && !cols.is_empty() && cols.iter().all(|c| c.parse::<f64>().is_err()) {
            continue;
        }
        if cols.len() != 2 {
            return Err(format!("line {}: two columns (T, value) expected", i + 1));
        }
        t.push(num(cols[0])?);
        v.push(num(cols[1])?);
    }
    Ok((t, v))
}

/// The NASA-7 record of `species` in a CHEMKIN THERMO text (fixed columns:
/// line 1 with the name in columns 1–18, T_low, T_high and T_mid in
/// columns 46–55, 56–65 and 66–73 and `1` in column 80; lines 2–4 with
/// five numbers of 15 columns each: a₁…a₇ of the upper range, then a₁…a₇
/// of the lower range). A `THERMO` header line with default temperatures
/// supplies T_mid when the record leaves it blank. Returns the range
/// limits and the coefficients from the lowest range up.
#[allow(clippy::type_complexity)]
pub fn nasa7_from_chemkin(text: &str, species: &str) -> Result<(Vec<f64>, Vec<Vec<f64>>)> {
    let lines: Vec<&str> = text.lines().collect();
    let field = |l: &str, a: usize, b: usize| -> String {
        l.chars()
            .skip(a)
            .take(b - a)
            .collect::<String>()
            .trim()
            .to_string()
    };
    let parse = |s: &str| -> Result<f64> {
        s.replace(['D', 'd'], "E")
            .parse::<f64>()
            .map_err(|_| anyhow::anyhow!("`{s}` is not a number"))
    };
    // default T_mid of the THERMO header (line after `THERMO`)
    let mut t_mid_default = None;
    for (i, l) in lines.iter().enumerate() {
        if l.trim_start().to_ascii_uppercase().starts_with("THERMO") {
            if let Some(h) = lines.get(i + 1) {
                let v: Vec<f64> = h
                    .split_whitespace()
                    .filter_map(|x| x.parse::<f64>().ok())
                    .collect();
                if v.len() >= 3 {
                    t_mid_default = Some(v[1]);
                }
            }
            break;
        }
    }
    let i = lines
        .iter()
        .position(|l| {
            l.len() >= 80
                && l.chars().nth(79) == Some('1')
                && l.split_whitespace().next() == Some(species)
        })
        .ok_or_else(|| anyhow::anyhow!("no NASA-7 record of `{species}`"))?;
    anyhow::ensure!(
        i + 3 < lines.len(),
        "the record of `{species}` has fewer than four lines"
    );
    let l1 = lines[i];
    let t_low = parse(&field(l1, 45, 55))?;
    let t_high = parse(&field(l1, 55, 65))?;
    let mid = field(l1, 65, 73);
    let t_mid = if mid.is_empty() {
        t_mid_default.ok_or_else(|| anyhow::anyhow!("no T_mid for `{species}`"))?
    } else {
        parse(&mid)?
    };
    let mut a = Vec::with_capacity(15);
    for (k, line) in lines[i + 1..i + 4].iter().enumerate() {
        let n = if k == 2 { 4 } else { 5 };
        for j in 0..n {
            a.push(parse(&field(line, 15 * j, 15 * (j + 1)))?);
        }
    }
    anyhow::ensure!(
        t_low < t_mid && t_mid < t_high,
        "the ranges of `{species}` are not increasing: {t_low}, {t_mid}, {t_high}"
    );
    Ok((
        vec![t_low, t_mid, t_high],
        vec![a[7..14].to_vec(), a[0..7].to_vec()],
    ))
}

/// The kind of a region (D112).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum RegionKind {
    /// A fluid region: the flow is solved on it.
    #[default]
    Fluid,
    /// A solid (conjugate heat transfer, with L1).
    Solid,
    /// A porous zone: the gas flows through a fixed matrix (D270).
    Porous,
}

/// `[regions.<name>]`: cell zones grouped into a region with a kind and a
/// material (D112).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct RegionSection {
    /// Fluid, solid or porous.
    pub kind: RegionKind,
    /// The cell zones of the region (mesh zone names).
    pub zones: Vec<String>,
    /// The material; default the material of `[fluid]`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub material: Option<String>,
    /// The first two principal axes of an orthotropic conductivity or of a
    /// porous resistance in this region (stage 7d2, L48; D270): one material
    /// oriented differently in every part; default the axes of the material.
    /// A porous region with direction-dependent values needs them, along
    /// the mesh axes in this version (POR002).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub axes: Option<[[f64; 3]; 2]>,
    /// The matrix of a porous region (D270): its resistances, or the
    /// channels of a honeycomb; the first of `axes` runs along the
    /// channels.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub porous: Option<crate::porous::PorousSection>,
    /// The temperature (K) at which the region's material is evaluated
    /// where its properties stay constant (D194): the density without the
    /// variable density, and every property of an isothermal run; default
    /// `[fluid] temperature`. Water at 300 K beside air at 450 K in one
    /// heat exchanger.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature: Option<f64>,
    /// The operating pressure (Pa) of a fluid region (D238 (3), D239): the
    /// absolute level of its gauge pressure and the pressure at which its
    /// material is evaluated; default `[fluid] operating_pressure`. With
    /// absolute pressure values the reader sets it from the region's own
    /// fluid part (the water at 1 bar beside the air at 1.5 bar).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub operating_pressure: Option<f64>,
    /// The start values of the region (D238): velocity, pressure,
    /// temperature and turbulence of a fluid region, the temperature of a
    /// solid one; before the case-wide `[initial]` and `[energy] initial`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial: Option<crate::start::RegionInitial>,
}

/// The JSON Schema of `[materials]`: the include list and, for every other
/// key, a material (D77: the keys of the case file are documented).
pub fn case_materials_schema(generator: &mut schemars::SchemaGenerator) -> schemars::Schema {
    let m = generator.subschema_for::<MaterialSection>();
    schemars::json_schema!({
        "type": "object",
        "description": "Materials of the case (D94): new ones, or changes of database materials with `base`. The key `include` lists user database files (D97).",
        "properties": {
            "include": {
                "type": "array",
                "items": {"type": "string"},
                "description": "User database files in the format of the built-in database, relative to the case file; searched after the case's own materials and before the built-in database."
            }
        },
        "additionalProperties": m
    })
}

/// The database file.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct DatabaseFile {
    format: String,
    #[allow(dead_code)]
    format_version: String,
    materials: BTreeMap<String, MaterialSection>,
}

const DATABASE_TOML: &str = include_str!("../data/materials.toml");

/// The built-in material database.
pub fn database() -> &'static BTreeMap<String, MaterialSection> {
    static DB: OnceLock<BTreeMap<String, MaterialSection>> = OnceLock::new();
    DB.get_or_init(|| {
        let f: DatabaseFile =
            toml::from_str(DATABASE_TOML).expect("the built-in material database parses");
        f.materials
    })
}

/// SHA-256 of the built-in database (the run record names it, D85).
pub fn database_hash() -> String {
    crate::engine::sha256(DATABASE_TOML.as_bytes())
}

/// The inputs of the property resolution that the case gives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conditions {
    /// The isothermal temperature (K).
    pub temperature: f64,
    /// The operating pressure (Pa).
    pub pressure: f64,
}

impl Default for Conditions {
    fn default() -> Self {
        Self {
            temperature: DEFAULT_TEMPERATURE,
            pressure: DEFAULT_PRESSURE,
        }
    }
}

fn mat_err(code: &'static str, key: &str, m: String) -> anyhow::Error {
    diag::err(Diagnostic::new(code, m).key(key))
}

/// The material `name` with its `base` chain merged (the case's materials
/// first, then the user database files in order, then the built-in
/// database).
pub fn lookup(name: &str, lib: &CaseMaterials) -> Result<MaterialSection> {
    Ok(lookup_layer(name, lib)?.0)
}

/// As [`lookup`], with the layer of `name` itself (0: the case, then the
/// files, the database last) and its label.
fn lookup_layer(name: &str, lib: &CaseMaterials) -> Result<(MaterialSection, usize, String)> {
    let layers = lib.layers();
    let mut seen: Vec<(String, usize)> = Vec::new();
    let mut chain: Vec<MaterialSection> = Vec::new();
    let mut found: Option<(usize, String)> = None;
    let mut cur = name.to_string();
    // a material's `base` is searched from its own layer on; a `base` of
    // its own name means the next layers (a case material that changes a
    // database one of the same name)
    let mut start = 0;
    loop {
        if seen.contains(&(cur.clone(), start)) || seen.len() > 32 {
            let names: Vec<&str> = seen.iter().map(|s| s.0.as_str()).collect();
            return Err(mat_err(
                diag::MAT_BASE,
                &format!("materials.{name}.base"),
                format!(
                    "the materials {} → {cur} refer to each other through `base` in a cycle",
                    names.join(" → ")
                ),
            ));
        }
        seen.push((cur.clone(), start));
        let Some(i) = (start..layers.len()).find(|&i| layers[i].0.contains_key(&cur)) else {
            let names = lib.names();
            let hint = match crate::bc::did_you_mean(&cur, names.iter().map(String::as_str)) {
                Some(c) => format!(" (did you mean `{c}`?)"),
                None => String::new(),
            };
            return Err(mat_err(
                diag::MAT_UNKNOWN,
                "fluid.material",
                format!("unknown material `{cur}`{hint}; `foehn materials` lists the database"),
            ));
        };
        let m = &layers[i].0[&cur];
        if found.is_none() {
            found = Some((i, layers[i].1.clone()));
        }
        chain.push(m.clone());
        match &m.base {
            Some(b) => {
                start = if *b == cur { i + 1 } else { i };
                cur = b.clone();
            }
            None => break,
        }
    }
    // merge from the root up
    let mut out = MaterialSection::default();
    for m in chain.into_iter().rev() {
        out.description = m.description.or(out.description);
        out.kind = m.kind.or(out.kind);
        out.temperature_range = m.temperature_range.or(out.temperature_range);
        out.molar_mass = m.molar_mass.or(out.molar_mass);
        out.molar_mass_source = m.molar_mass_source.or(out.molar_mass_source);
        out.density = m.density.or(out.density);
        out.viscosity = m.viscosity.or(out.viscosity);
        out.specific_heat = m.specific_heat.or(out.specific_heat);
        out.conductivity = m.conductivity.or(out.conductivity);
    }
    out.base = None;
    let (layer, label) = found.unwrap_or((0, String::new()));
    Ok((out, layer, label))
}

/// One resolved property.
#[derive(Debug, Clone, PartialEq)]
pub struct Resolved {
    /// The value (a number, or a program in x, y, z).
    pub value: Value,
    /// How it was obtained (form, T, source), for the resolved case.
    pub origin: String,
}

/// The properties of one material at the case conditions.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialProperties {
    /// The material (`None`: the `[fluid]` shorthand without a material).
    pub material: Option<String>,
    /// Density.
    pub density: Resolved,
    /// Viscosity.
    pub viscosity: Resolved,
}

/// Monotone cubic interpolation of Steffen (1990), *A simple method for
/// monotonic interpolation in one dimension*, A&A 239, 443–450.
pub fn steffen(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    // a diverging run: the NaN goes on to the residual check, which reports
    // the divergence (D202), instead of an index out of bounds here
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    let i = xs.partition_point(|&v| v <= x) - 1;
    if n == 2 {
        return ys[0] + (x - xs[0]) / (xs[1] - xs[0]) * (ys[1] - ys[0]);
    }
    let h = |j: usize| xs[j + 1] - xs[j];
    let s = |j: usize| (ys[j + 1] - ys[j]) / h(j);
    let sign = |v: f64| {
        if v > 0.0 {
            1.0
        } else if v < 0.0 {
            -1.0
        } else {
            0.0
        }
    };
    let slope = |j: usize| -> f64 {
        if j == 0 {
            let p = s(0) * (1.0 + h(0) / (h(0) + h(1))) - s(1) * h(0) / (h(0) + h(1));
            if p * s(0) <= 0.0 {
                0.0
            } else if p.abs() > 2.0 * s(0).abs() {
                2.0 * s(0)
            } else {
                p
            }
        } else if j == n - 1 {
            let (a, b) = (n - 2, n - 3);
            let p = s(a) * (1.0 + h(a) / (h(a) + h(b))) - s(b) * h(a) / (h(a) + h(b));
            if p * s(a) <= 0.0 {
                0.0
            } else if p.abs() > 2.0 * s(a).abs() {
                2.0 * s(a)
            } else {
                p
            }
        } else {
            let (s0, s1) = (s(j - 1), s(j));
            let p = (s0 * h(j) + s1 * h(j - 1)) / (h(j - 1) + h(j));
            (sign(s0) + sign(s1)) * s0.abs().min(s1.abs()).min(0.5 * p.abs())
        }
    };
    let (y0, y1) = (slope(i), slope(i + 1));
    let hi = h(i);
    let si = s(i);
    let a = (y0 + y1 - 2.0 * si) / (hi * hi);
    let b = (3.0 * si - 2.0 * y0 - y1) / hi;
    let d = x - xs[i];
    ((a * d + b) * d + y0) * d + ys[i]
}

/// Evaluates a property of `material` (name for messages) at the
/// conditions. `space` allows x, y, z (the viscosity).
#[allow(clippy::too_many_arguments)]
fn evaluate(
    p: &Property,
    what: &str,
    owner: &str,
    m: &MaterialSection,
    env: &Env,
    cond: Conditions,
    space: bool,
    warnings: &mut Vec<Diagnostic>,
) -> Result<Resolved> {
    let key = format!("{owner}.{what}");
    let t = cond.temperature;
    let mut tenv = env.clone();
    tenv.constants.insert("T".into(), t);
    let allowed = if space { Vars::SPACE } else { Vars::NONE };
    let warn_range = |range: Option<[f64; 2]>, warnings: &mut Vec<Diagnostic>| -> f64 {
        if let Some([lo, hi]) = range
            && (t < lo || t > hi)
        {
            warnings.push(
                Diagnostic::new(
                    diag::MAT_RANGE,
                    format!(
                        "{key}: T = {t} K is outside the data range [{lo}, {hi}] K; \
                         the property is taken at {} K",
                        t.clamp(lo, hi)
                    ),
                )
                .key(key.clone()),
            );
            return t.clamp(lo, hi);
        }
        t
    };
    let num = |v: f64, origin: String| Resolved {
        value: Value::Number(v),
        origin,
    };
    let r = match p {
        Property::Number(v) => num(*v, format!("{v:?}")),
        Property::Expr(s) => {
            let v = expr::value(&Scalar::Expr(s.clone()), &tenv, allowed, &key).map_err(|e| {
                if !space && format!("{e:#}").contains("but this key takes no variable") {
                    mat_err(
                        diag::MAT_FORM,
                        &key,
                        format!(
                            "{key}: `{s}` varies in space; the density is one number for the \
                             whole fluid until the energy equation (L1) brings a variable density"
                        ),
                    )
                } else {
                    e
                }
            })?;
            let origin = format!("expression \"{s}\" at T = {t} K");
            match v.constant() {
                Some(c) => num(c, origin),
                None => Resolved { value: v, origin },
            }
        }
        Property::Form(f) => {
            let src = f.source().map(|s| format!("; {s}")).unwrap_or_default();
            let tt = warn_range(f.range(), warnings);
            let origin = |extra: &str| format!("{}{extra} at T = {tt} K{src}", f.name());
            match f {
                Form::Constant { value, .. } => num(*value, origin("")),
                Form::Polynomial { coefficients, .. } => {
                    let mut v = 0.0;
                    for c in coefficients.iter().rev() {
                        v = v * tt + c;
                    }
                    num(v, origin(""))
                }
                Form::Sutherland {
                    mu_ref, t_ref, s, ..
                } => {
                    let r = tt / t_ref;
                    num(
                        mu_ref * r * foehncraft_core::math::sqrt(r) * (t_ref + s) / (tt + s),
                        origin(""),
                    )
                }
                Form::PowerLaw {
                    value_ref,
                    t_ref,
                    exponent,
                    ..
                } => num(
                    value_ref * foehncraft_core::math::pow(tt / t_ref, *exponent),
                    origin(""),
                ),
                Form::Table {
                    temperature,
                    values,
                    file,
                    interpolation,
                    ..
                } => {
                    if file.is_some() {
                        return Err(mat_err(
                            diag::MAT_FORM,
                            &key,
                            format!(
                                "{key}: property tables from files come with the importers of D103; give `temperature` and `values`"
                            ),
                        ));
                    }
                    check_table(temperature, values, &key)?;
                    let v = match interpolation {
                        Interpolation::Linear => linear(temperature, values, tt),
                        Interpolation::MonotoneCubic => steffen(temperature, values, tt),
                    };
                    let name = match interpolation {
                        Interpolation::Linear => " (linear)",
                        Interpolation::MonotoneCubic => " (monotone-cubic)",
                    };
                    num(v, origin(name))
                }
                Form::Expression { expression, .. } => {
                    let inner = Property::Expr(expression.clone());
                    let mut r = evaluate(&inner, what, owner, m, env, cond, space, warnings)?;
                    r.origin.push_str(&src);
                    r
                }
                Form::IdealGas { molar_mass, .. } => {
                    let mm = molar_mass.or(m.molar_mass).ok_or_else(|| {
                        mat_err(
                            diag::MAT_FORM,
                            &key,
                            format!("{key}: the ideal-gas density needs `molar_mass`"),
                        )
                    })?;
                    num(
                        cond.pressure * mm / (R_GAS * t),
                        format!(
                            "ideal-gas, p = {} Pa, M = {mm} kg/mol, T = {t} K{src}",
                            cond.pressure
                        ),
                    )
                }
                Form::Nasa7 {
                    temperature,
                    coefficients,
                    ..
                }
                | Form::Nasa9 {
                    temperature,
                    coefficients,
                    ..
                } => {
                    if what != "specific_heat" {
                        return Err(mat_err(
                            diag::MAT_FORM,
                            &key,
                            format!("{key}: the {} form gives the specific heat only", f.name()),
                        ));
                    }
                    let nine = matches!(f, Form::Nasa9 { .. });
                    check_nasa(temperature, coefficients, nine, &key)?;
                    let mm = m.molar_mass.ok_or_else(|| {
                        mat_err(
                            diag::MAT_FORM,
                            &key,
                            format!(
                                "{key}: the {} form needs the material's `molar_mass`",
                                f.name()
                            ),
                        )
                    })?;
                    num(
                        R_GAS / mm * nasa_cp_over_r(temperature, coefficients, nine, tt),
                        origin(&format!(", M = {mm} kg/mol")),
                    )
                }
            }
        }
    };
    if let Some(v) = r.value.constant()
        && !(v > 0.0 && v.is_finite())
    {
        return Err(mat_err(
            diag::MAT_VALUE,
            &key,
            format!("{key} = {v} at T = {t} K must be positive and finite"),
        ));
    }
    Ok(r)
}

/// c_p/R of NASA-7 (`nine` false) or NASA-9 polynomials at `t` (inside the
/// ranges; the range of the upper limit at a shared limit).
pub fn nasa_cp_over_r(temperature: &[f64], coefficients: &[Vec<f64>], nine: bool, t: f64) -> f64 {
    let n = coefficients.len();
    let i = (0..n).find(|&i| t <= temperature[i + 1]).unwrap_or(n - 1);
    let a = &coefficients[i];
    if nine {
        a[0] / (t * t) + a[1] / t + a[2] + t * (a[3] + t * (a[4] + t * (a[5] + t * a[6])))
    } else {
        a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4])))
    }
}

fn check_nasa(t: &[f64], c: &[Vec<f64>], nine: bool, key: &str) -> Result<()> {
    let need = if nine { [7, 9] } else { [7, 7] };
    let ok = !c.is_empty()
        && t.len() == c.len() + 1
        && t.iter().all(|x| x.is_finite() && *x > 0.0)
        && t.windows(2).all(|w| w[1] > w[0])
        && c.iter()
            .all(|r| need.contains(&r.len()) && r.iter().all(|x| x.is_finite()));
    if ok {
        return Ok(());
    }
    Err(mat_err(
        diag::MAT_FORM,
        key,
        format!(
            "{key}: NASA-{} needs n ranges with {} coefficients each and n + 1 increasing \
             positive `temperature` limits (has {} ranges, {} limits)",
            if nine { 9 } else { 7 },
            if nine { "7 or 9" } else { "7" },
            c.len(),
            t.len()
        ),
    ))
}

fn linear(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    let i = xs.partition_point(|&v| v <= x) - 1;
    ys[i] + (x - xs[i]) / (xs[i + 1] - xs[i]) * (ys[i + 1] - ys[i])
}

fn check_table(t: &[f64], v: &[f64], key: &str) -> Result<()> {
    if t.len() < 2 || t.len() != v.len() {
        return Err(mat_err(
            diag::MAT_FORM,
            key,
            format!(
                "{key}: a table needs `temperature` and `values` of the same length, at least 2 \
                 (has {} and {})",
                t.len(),
                v.len()
            ),
        ));
    }
    if t.iter().chain(v).any(|x| !x.is_finite()) || t.windows(2).any(|w| w[1] <= w[0]) {
        return Err(mat_err(
            diag::MAT_FORM,
            key,
            format!("{key}: the temperatures must be finite and increase strictly"),
        ));
    }
    Ok(())
}

/// How messages name a material of layer `layer` (0: the case).
fn owner_label(name: &str, layer: usize, label: &str) -> String {
    match (layer, label) {
        (0, _) => format!("materials.{name}"),
        (_, "database") => format!("material {name}"),
        (_, file) => format!("material {name} ({file})"),
    }
}

/// The properties of a material (`None`: the shorthand of `[fluid]`), with
/// the `[fluid]` density and viscosity replacing the material's.
pub fn resolve_material(
    name: Option<&str>,
    case_materials: &CaseMaterials,
    fluid_density: Option<&Scalar>,
    fluid_viscosity: Option<&Scalar>,
    env: &Env,
    cond: Conditions,
    warnings: &mut Vec<Diagnostic>,
) -> Result<MaterialProperties> {
    let (m, owner) = match name {
        Some(n) => {
            let (m, layer, label) = lookup_layer(n, case_materials)?;
            (m, owner_label(n, layer, &label))
        }
        None => (MaterialSection::default(), "fluid".to_string()),
    };
    if m.kind == Some(MaterialKind::Solid) {
        return Err(mat_err(
            diag::REG_NOT_FLUID,
            "fluid.material",
            format!(
                "`{}` is a solid; solids come with the energy equation (L1)",
                name.unwrap_or("")
            ),
        ));
    }
    let pick = |given: Option<&Scalar>, own: &Option<Property>, default: f64| -> (Property, bool) {
        match (given, own) {
            (Some(s), _) => (Property::from(s.clone()), true),
            (None, Some(p)) => (p.clone(), false),
            (None, None) => (Property::Number(default), false),
        }
    };
    let (rho_p, rho_given) = pick(fluid_density, &m.density, 1.0);
    let (mu_p, mu_given) = pick(fluid_viscosity, &m.viscosity, 1e-3);
    if name.is_some() {
        for (what, own, given) in [
            ("density", &m.density, rho_given),
            ("viscosity", &m.viscosity, mu_given),
        ] {
            if own.is_none() && !given {
                return Err(mat_err(
                    diag::MAT_FORM,
                    &format!("{owner}.{what}"),
                    format!("{owner} has no {what}; give it there or in [fluid] {what}"),
                ));
            }
        }
    }
    let rho_owner = if rho_given {
        "fluid".to_string()
    } else {
        owner.clone()
    };
    let mu_owner = if mu_given {
        "fluid".to_string()
    } else {
        owner.clone()
    };
    let density = evaluate(
        &rho_p, "density", &rho_owner, &m, env, cond, false, warnings,
    )?;
    let viscosity = evaluate(&mu_p, "viscosity", &mu_owner, &m, env, cond, true, warnings)?;
    Ok(MaterialProperties {
        material: name.map(str::to_string),
        density,
        viscosity,
    })
}

/// A resolved thermal conductivity.
#[derive(Debug, Clone, PartialEq)]
pub enum ResolvedConductivity {
    /// One value.
    Isotropic(Resolved),
    /// Principal values, the orthonormal axes (rows e₁, e₂, e₃) and the
    /// tensor K = Σᵢ kᵢ eᵢ eᵢᵀ in the mesh frame.
    Orthotropic {
        /// k₁, k₂, k₃ (W/(m K)).
        principal: [f64; 3],
        /// e₁, e₂, e₃.
        axes: [[f64; 3]; 3],
        /// K.
        tensor: [[f64; 3]; 3],
        /// How it was obtained.
        origin: String,
    },
}

impl ResolvedConductivity {
    /// The conductivity tensor (k I for an isotropic value).
    pub fn tensor(&self) -> [[f64; 3]; 3] {
        match self {
            ResolvedConductivity::Isotropic(r) => {
                let k = r.value.constant().unwrap_or(f64::NAN);
                [[k, 0.0, 0.0], [0.0, k, 0.0], [0.0, 0.0, k]]
            }
            ResolvedConductivity::Orthotropic { tensor, .. } => *tensor,
        }
    }
}

/// The thermal data of one material at one temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct Thermal {
    /// The material.
    pub material: String,
    /// Gas, liquid or solid.
    pub kind: Option<MaterialKind>,
    /// Density ρ.
    pub density: Option<Resolved>,
    /// Specific heat c_p.
    pub specific_heat: Option<Resolved>,
    /// Thermal conductivity k.
    pub conductivity: Option<ResolvedConductivity>,
}

impl Thermal {
    /// MAT007 unless c_p and k are both given (a fluid in a run with the
    /// energy equation; solids are checked by [`thermal_properties`]).
    pub fn require(&self) -> Result<()> {
        let missing: Vec<&str> = [
            ("specific_heat", self.specific_heat.is_none()),
            ("conductivity", self.conductivity.is_none()),
        ]
        .iter()
        .filter(|x| x.1)
        .map(|x| x.0)
        .collect();
        if missing.is_empty() {
            return Ok(());
        }
        Err(mat_err(
            diag::MAT_THERMAL,
            &format!("materials.{}", self.material),
            format!(
                "material `{}` has no {}; the energy equation needs c_p and k of every fluid \
                 (give them in [materials.{}], with `base = \"{}\"` to keep the rest)",
                self.material,
                missing.join(" and no "),
                self.material,
                self.material
            ),
        ))
    }
}

fn normalised(v: [f64; 3]) -> Option<[f64; 3]> {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    (n > 0.0 && n.is_finite()).then(|| [v[0] / n, v[1] / n, v[2] / n])
}

/// The orthonormal axes of an orthotropic conductivity (MAT008 without
/// both axes, or with parallel ones).
pub(crate) fn axes_of(o: &Orthotropic, key: &str) -> Result<[[f64; 3]; 3]> {
    let bad = |m: String| mat_err(diag::MAT_AXES, key, m);
    let Some([a1, a2]) = o.axes else {
        return Err(bad(format!(
            "{key}: an orthotropic conductivity needs `axes` (the directions of k₁ and k₂ in \
             the mesh frame)"
        )));
    };
    let e1 = normalised(a1).ok_or_else(|| bad(format!("{key}: `axes[0]` has no length")))?;
    let d = a2[0] * e1[0] + a2[1] * e1[1] + a2[2] * e1[2];
    let e2 = normalised([a2[0] - d * e1[0], a2[1] - d * e1[1], a2[2] - d * e1[2]])
        .filter(|_| {
            let n2 = (a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2]).sqrt();
            (a2[0] - d * e1[0])
                .hypot(a2[1] - d * e1[1])
                .hypot(a2[2] - d * e1[2])
                > 1e-6 * n2
        })
        .ok_or_else(|| {
            bad(format!(
                "{key}: `axes[1]` is parallel to `axes[0]` or has no length"
            ))
        })?;
    let e3 = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    Ok([e1, e2, e3])
}

/// The thermal data of material `name` at the temperature of `cond`
/// (forms evaluated there, clamped with MAT003 outside their ranges).
/// A solid must have ρ, c_p and k (MAT006); an orthotropic k needs its
/// axes (MAT008).
pub fn thermal_properties(
    name: &str,
    lib: &CaseMaterials,
    env: &Env,
    cond: Conditions,
    warnings: &mut Vec<Diagnostic>,
) -> Result<Thermal> {
    let (m, layer, label) = lookup_layer(name, lib)?;
    let owner = owner_label(name, layer, &label);
    if m.kind == Some(MaterialKind::Solid) {
        let missing: Vec<&str> = [
            ("density", m.density.is_none()),
            ("specific_heat", m.specific_heat.is_none()),
            ("conductivity", m.conductivity.is_none()),
        ]
        .iter()
        .filter(|x| x.1)
        .map(|x| x.0)
        .collect();
        if !missing.is_empty() {
            return Err(mat_err(
                diag::MAT_SOLID,
                &owner,
                format!(
                    "{owner} is a solid without {}: a solid needs density, specific_heat and \
                     conductivity",
                    missing.join(", ")
                ),
            ));
        }
    }
    let one = |p: &Property, what: &str, w: &mut Vec<Diagnostic>| {
        evaluate(p, what, &owner, &m, env, cond, false, w)
    };
    let density = match &m.density {
        Some(p) => Some(one(p, "density", warnings)?),
        None => None,
    };
    let specific_heat = match &m.specific_heat {
        Some(p) => Some(one(p, "specific_heat", warnings)?),
        None => None,
    };
    let conductivity = match &m.conductivity {
        None => None,
        Some(Conductivity::Isotropic(p)) => Some(ResolvedConductivity::Isotropic(one(
            p,
            "conductivity",
            warnings,
        )?)),
        Some(Conductivity::Orthotropic(o)) => {
            let key = format!("{owner}.conductivity");
            let axes = axes_of(o, &key)?;
            let mut k = [0.0; 3];
            let mut origins = Vec::new();
            for (i, p) in o.values.iter().enumerate() {
                let r = one(p, &format!("conductivity.values[{i}]"), warnings)?;
                k[i] = r.value.constant().unwrap_or(f64::NAN);
                origins.push(r.origin);
            }
            let mut t = [[0.0; 3]; 3];
            for (ki, e) in k.iter().zip(&axes) {
                for a in 0..3 {
                    for b in 0..3 {
                        t[a][b] += ki * e[a] * e[b];
                    }
                }
            }
            let src = o
                .source
                .as_deref()
                .map(|s| format!("; {s}"))
                .unwrap_or_default();
            Some(ResolvedConductivity::Orthotropic {
                principal: k,
                axes,
                tensor: t,
                origin: format!("orthotropic [{}]{src}", origins.join(" | ")),
            })
        }
    };
    Ok(Thermal {
        material: name.to_string(),
        kind: m.kind,
        density,
        specific_heat,
        conductivity,
    })
}

/// A property as a function of the temperature, for the energy equation
/// (stage 7a1): c_p(T), k(T) and μ(T) per cell and face, refreshed every
/// outer iteration. The forms of D95 are evaluated directly (clamped to
/// their data range, as `evaluate` does with MAT003 at the case
/// temperature); an expression in `T` is sampled every
/// [`Curve::SAMPLE_STEP`] K from [`Curve::SAMPLE_RANGE`] and interpolated
/// linearly between the samples.
#[derive(Debug, Clone, PartialEq)]
pub struct Curve {
    kind: CurveKind,
    /// The data range (K), if the form has one.
    pub range: Option<[f64; 2]>,
    /// How it was obtained (form, source), for the log.
    pub origin: String,
}

#[derive(Debug, Clone, PartialEq)]
enum CurveKind {
    Constant(f64),
    Polynomial(Vec<f64>),
    Sutherland {
        mu_ref: f64,
        t_ref: f64,
        s: f64,
    },
    PowerLaw {
        value_ref: f64,
        t_ref: f64,
        exponent: f64,
    },
    Table {
        t: Vec<f64>,
        v: Vec<f64>,
        cubic: bool,
    },
    Nasa {
        t: Vec<f64>,
        c: Vec<Vec<f64>>,
        nine: bool,
        r: f64,
    },
    IdealGas {
        pm_over_r: f64,
    },
    Sampled {
        t0: f64,
        dt: f64,
        v: Vec<f64>,
    },
}

impl Curve {
    /// The step of the samples of an expression in `T` (K).
    pub const SAMPLE_STEP: f64 = 0.5;
    /// The range of the samples of an expression in `T` (K).
    pub const SAMPLE_RANGE: [f64; 2] = [1.0, 4000.0];

    /// A constant.
    pub fn constant(v: f64) -> Self {
        Self {
            kind: CurveKind::Constant(v),
            range: None,
            origin: format!("{v:?}"),
        }
    }

    /// True for a constant (the energy module then keeps the value).
    pub fn is_constant(&self) -> bool {
        matches!(self.kind, CurveKind::Constant(_))
    }

    /// p M / R (kg K/m³) of an ideal-gas density curve, at the pressure it
    /// was built for (`None` for any other form).
    pub fn ideal_gas(&self) -> Option<f64> {
        match self.kind {
            CurveKind::IdealGas { pm_over_r } => Some(pm_over_r),
            _ => None,
        }
    }

    /// The value at the temperature `t` (K), with `t` clamped to the data
    /// range.
    pub fn at(&self, t: f64) -> f64 {
        let t = match self.range {
            Some([lo, hi]) => t.clamp(lo, hi),
            None => t,
        };
        match &self.kind {
            CurveKind::Constant(v) => *v,
            CurveKind::Polynomial(c) => {
                let mut v = 0.0;
                for a in c.iter().rev() {
                    v = v * t + a;
                }
                v
            }
            CurveKind::Sutherland { mu_ref, t_ref, s } => {
                let r = t / t_ref;
                mu_ref * r * foehncraft_core::math::sqrt(r) * (t_ref + s) / (t + s)
            }
            CurveKind::PowerLaw {
                value_ref,
                t_ref,
                exponent,
            } => value_ref * foehncraft_core::math::pow(t / t_ref, *exponent),
            CurveKind::Table { t: xs, v, cubic } => {
                if *cubic {
                    steffen(xs, v, t)
                } else {
                    linear(xs, v, t)
                }
            }
            CurveKind::Nasa { t: lim, c, nine, r } => r * nasa_cp_over_r(lim, c, *nine, t),
            CurveKind::IdealGas { pm_over_r } => pm_over_r / t,
            CurveKind::Sampled { t0, dt, v } => {
                let x = ((t - t0) / dt).clamp(0.0, (v.len() - 1) as f64);
                let i = (x.floor() as usize).min(v.len() - 2);
                let w = x - i as f64;
                v[i] + w * (v[i + 1] - v[i])
            }
        }
    }
}

/// The curve of property `p` (`what` of material `m`, named `owner` in
/// messages) over T: checked as `evaluate` checks it at the case
/// temperature `cond` (MAT002, MAT005), so that a curve that is built
/// evaluates everywhere.
pub fn curve(
    p: &Property,
    what: &str,
    owner: &str,
    m: &MaterialSection,
    env: &Env,
    cond: Conditions,
) -> Result<Curve> {
    // the checks and the origin text of the case temperature
    let mut warnings = Vec::new();
    let r = evaluate(p, what, owner, m, env, cond, false, &mut warnings)?;
    let origin = r.origin;
    let sampled = |text: &str| -> Result<Curve> {
        let [lo, hi] = Curve::SAMPLE_RANGE;
        let n = ((hi - lo) / Curve::SAMPLE_STEP).round() as usize + 1;
        let mut v = Vec::with_capacity(n);
        // samples where the expression is not positive are kept: the
        // energy module refuses such a value where the temperature reaches
        // it (MAT005)
        let key = format!("{owner}.{what}");
        let mut tenv = env.clone();
        for i in 0..n {
            let t = lo + i as f64 * Curve::SAMPLE_STEP;
            tenv.constants.insert("T".into(), t);
            let x = expr::value(&Scalar::Expr(text.to_string()), &tenv, Vars::NONE, &key)?;
            v.push(x.constant().unwrap_or(f64::NAN));
        }
        Ok(Curve {
            kind: CurveKind::Sampled {
                t0: lo,
                dt: Curve::SAMPLE_STEP,
                v,
            },
            range: Some(Curve::SAMPLE_RANGE),
            origin: origin.clone(),
        })
    };
    let kind = match p {
        Property::Number(v) => CurveKind::Constant(*v),
        Property::Expr(s) => return sampled(s),
        Property::Form(f) => match f {
            Form::Constant { value, .. } => CurveKind::Constant(*value),
            Form::Polynomial { coefficients, .. } => CurveKind::Polynomial(coefficients.clone()),
            Form::Sutherland {
                mu_ref, t_ref, s, ..
            } => CurveKind::Sutherland {
                mu_ref: *mu_ref,
                t_ref: *t_ref,
                s: *s,
            },
            Form::PowerLaw {
                value_ref,
                t_ref,
                exponent,
                ..
            } => CurveKind::PowerLaw {
                value_ref: *value_ref,
                t_ref: *t_ref,
                exponent: *exponent,
            },
            Form::Table {
                temperature,
                values,
                interpolation,
                ..
            } => CurveKind::Table {
                t: temperature.clone(),
                v: values.clone(),
                cubic: *interpolation == Interpolation::MonotoneCubic,
            },
            Form::Expression { expression, .. } => return sampled(expression),
            Form::IdealGas { molar_mass, .. } => CurveKind::IdealGas {
                pm_over_r: cond.pressure * molar_mass.or(m.molar_mass).unwrap_or(f64::NAN) / R_GAS,
            },
            Form::Nasa7 {
                temperature,
                coefficients,
                ..
            }
            | Form::Nasa9 {
                temperature,
                coefficients,
                ..
            } => CurveKind::Nasa {
                t: temperature.clone(),
                c: coefficients.clone(),
                nine: matches!(f, Form::Nasa9 { .. }),
                r: R_GAS / m.molar_mass.unwrap_or(f64::NAN),
            },
        },
    };
    let range = match p {
        Property::Form(f) => f.range(),
        _ => None,
    };
    Ok(Curve {
        kind,
        range,
        origin,
    })
}

/// The resolved property program of a case (D94 (3)): the default fluid
/// and one entry per region.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertySet {
    /// Conditions of the evaluation.
    pub conditions: Conditions,
    /// The `[fluid]` material (on every zone outside the regions).
    pub fluid: MaterialProperties,
    /// `(region name, zones, properties)` in name order.
    pub regions: Vec<(String, Vec<String>, MaterialProperties)>,
    /// Range warnings (MAT003).
    pub warnings: Vec<Diagnostic>,
}

impl PropertySet {
    /// The density of the `[fluid]` material: the reference value (and the
    /// value in every cell unless [`density_varies`](Self::density_varies)).
    pub fn density(&self) -> f64 {
        self.fluid.density.value.constant().unwrap_or(f64::NAN)
    }

    /// True when a fluid region has a density of its own (D194: several
    /// fluids in one case).
    pub fn density_varies(&self) -> bool {
        let r0 = self.density().to_bits();
        self.regions
            .iter()
            .any(|r| r.2.density.value.constant().map(f64::to_bits) != Some(r0))
    }

    /// The density of each cell zone of `mesh`.
    fn zone_density(&self, mesh: &foehncraft_mesh::Mesh) -> Vec<f64> {
        mesh.zones
            .iter()
            .map(|z| {
                self.regions
                    .iter()
                    .find(|r| r.1.iter().any(|n| n == z))
                    .map_or(&self.fluid.density, |r| &r.2.density)
                    .value
                    .constant()
                    .unwrap_or(f64::NAN)
            })
            .collect()
    }

    /// The density per cell and per face of `mesh` (D194; `None` when it is
    /// one number everywhere, which keeps the solver's constant path bit for
    /// bit). Fluids of different densities never share a face (REG004);
    /// the face takes the owner's value.
    pub fn cell_face_density(&self, mesh: &foehncraft_mesh::Mesh) -> Option<(Vec<f64>, Vec<f64>)> {
        if !self.density_varies() {
            return None;
        }
        let zr = self.zone_density(mesh);
        let cell: Vec<f64> = (0..mesh.n_cells)
            .map(|c| zr[mesh.cell_zone[c] as usize])
            .collect();
        let face: Vec<f64> = mesh.owner.iter().map(|&o| cell[o as usize]).collect();
        Some((cell, face))
    }

    /// The viscosity of the `[fluid]` material, as the reference value of
    /// Reynolds numbers and exact solutions; for a viscosity that varies in
    /// space, its value at the origin.
    pub fn reference_viscosity(&self) -> f64 {
        self.fluid.viscosity.value.at([0.0; 3], 0.0)
    }

    /// True when the viscosity is not one number on the whole mesh.
    pub fn viscosity_varies(&self) -> bool {
        let mu0 = self.fluid.viscosity.value.constant();
        mu0.is_none()
            || self
                .regions
                .iter()
                .any(|r| r.2.viscosity.value.constant().map(f64::to_bits) != mu0.map(f64::to_bits))
    }

    /// The program as text: what the resolved case shows and the case hash
    /// includes.
    pub fn text(&self) -> String {
        let mut s = String::new();
        let c = self.conditions;
        let _ = writeln!(
            s,
            "# resolved property program (D94): T = {} K, p = {} Pa, database {}",
            c.temperature,
            c.pressure,
            &database_hash()[..16]
        );
        let mut one = |label: &str, m: &MaterialProperties| {
            let mat = m.material.as_deref().unwrap_or("(fluid shorthand)");
            for (what, r) in [("density", &m.density), ("viscosity", &m.viscosity)] {
                let v = match &r.value {
                    Value::Number(v) => format!("{v:?}"),
                    Value::Program(p) => format!("\"{}\"", p.source),
                };
                let _ = writeln!(s, "#   {label} [{mat}] {what} = {v}  ({})", r.origin);
            }
        };
        one("fluid", &self.fluid);
        for (name, zones, m) in &self.regions {
            one(&format!("region {name} {zones:?}"), m);
        }
        s
    }

    /// The viscosity per cell of a mesh (zone by zone; expressions at the
    /// cell centroids).
    pub fn cell_viscosity(
        &self,
        mesh: &foehncraft_mesh::Mesh,
        geo: &foehncraft_mesh::Geometry,
    ) -> Vec<f64> {
        let zone_mu = self.zone_values(mesh);
        (0..mesh.n_cells)
            .map(|c| zone_mu[mesh.cell_zone[c] as usize].at(geo.cell_centroid[c], 0.0))
            .collect()
    }

    /// The viscosity of each cell zone of `mesh`.
    fn zone_values<'a>(&'a self, mesh: &foehncraft_mesh::Mesh) -> Vec<&'a Value> {
        mesh.zones
            .iter()
            .map(|z| {
                self.regions
                    .iter()
                    .find(|r| r.1.iter().any(|n| n == z))
                    .map_or(&self.fluid.viscosity.value, |r| &r.2.viscosity.value)
            })
            .collect()
    }

    /// The viscosity per face (`None` when it is one number everywhere,
    /// which keeps the solver's constant path bit for bit). Inside a zone
    /// an expression is evaluated at the face centroid; between zones the
    /// harmonic mean of the two cell values (D99).
    pub fn face_viscosity(
        &self,
        mesh: &foehncraft_mesh::Mesh,
        geo: &foehncraft_mesh::Geometry,
    ) -> Option<Arc<Vec<f64>>> {
        if !self.viscosity_varies() {
            return None;
        }
        let zone_mu = self.zone_values(mesh);
        let ni = mesh.n_interior_faces();
        let mut mu = Vec::with_capacity(mesh.n_faces());
        let cell = |c: usize| zone_mu[mesh.cell_zone[c] as usize].at(geo.cell_centroid[c], 0.0);
        for f in 0..mesh.n_faces() {
            let o = mesh.owner[f] as usize;
            let zo = mesh.cell_zone[o];
            if f < ni {
                let n = mesh.neighbour[f] as usize;
                if mesh.cell_zone[n] != zo {
                    let (a, b) = (cell(o), cell(n));
                    mu.push(2.0 * a * b / (a + b));
                    continue;
                }
            }
            mu.push(zone_mu[zo as usize].at(geo.face_centroid[f], 0.0));
        }
        Some(Arc::new(mu))
    }
}

/// The JSON description of a database material at the conditions (for
/// `foehn materials`).
pub fn describe(name: &str, cond: Conditions) -> Result<serde_json::Value> {
    describe_in(name, &CaseMaterials::default(), cond)
}

/// As [`describe`], for a material of the case materials or user database
/// files `lib` as well.
pub fn describe_in(name: &str, lib: &CaseMaterials, cond: Conditions) -> Result<serde_json::Value> {
    let (m, _, label) = lookup_layer(name, lib)?;
    let mut warnings = Vec::new();
    let env = Env::default();
    let th = thermal_properties(name, lib, &env, cond, &mut warnings)?;
    let prop = |r: &Resolved| {
        serde_json::json!({
            "value": r.value.constant(),
            "origin": r.origin,
        })
    };
    // fluids: density and viscosity as the flow takes them
    let (density, viscosity) = if m.kind == Some(MaterialKind::Solid) {
        (th.density.as_ref().map(prop), None)
    } else {
        let p = resolve_material(Some(name), lib, None, None, &env, cond, &mut warnings)?;
        (Some(prop(&p.density)), Some(prop(&p.viscosity)))
    };
    let conductivity = th.conductivity.as_ref().map(|k| match k {
        ResolvedConductivity::Isotropic(r) => prop(r),
        ResolvedConductivity::Orthotropic {
            principal,
            axes,
            tensor,
            origin,
        } => serde_json::json!({
            "principal": principal,
            "axes": axes,
            "tensor": tensor,
            "origin": origin,
        }),
    });
    let form = |p: &Option<Property>| {
        p.as_ref().map(|d| match d {
            Property::Form(f) => f.name(),
            Property::Number(_) => "number",
            Property::Expr(_) => "expression",
        })
    };
    let kform = m.conductivity.as_ref().map(|k| match k {
        Conductivity::Isotropic(p) => form(&Some(p.clone())).unwrap_or("number"),
        Conductivity::Orthotropic(_) => "orthotropic",
    });
    warnings.dedup_by(|a, b| a.message == b.message);
    Ok(serde_json::json!({
        "name": name,
        "from": label,
        "description": m.description,
        "kind": m.kind,
        "temperature_range": m.temperature_range,
        "molar_mass": m.molar_mass,
        "temperature": cond.temperature,
        "operating_pressure": cond.pressure,
        "density": density,
        "viscosity": viscosity,
        "specific_heat": th.specific_heat.as_ref().map(prop),
        "conductivity": conductivity,
        "forms": {
            "density": form(&m.density),
            "viscosity": form(&m.viscosity),
            "specific_heat": form(&m.specific_heat),
            "conductivity": kform,
        },
        "warnings": warnings.iter().map(|w| w.message.clone()).collect::<Vec<_>>(),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at(name: &str, t: f64) -> MaterialProperties {
        let mut w = Vec::new();
        resolve_material(
            Some(name),
            &CaseMaterials::default(),
            None,
            None,
            &Env::default(),
            Conditions {
                temperature: t,
                pressure: DEFAULT_PRESSURE,
            },
            &mut w,
        )
        .unwrap()
    }

    #[test]
    fn database_parses_and_resolves() {
        for name in database().keys() {
            let m = at(name, 300.0);
            let (r, mu) = (
                m.density.value.constant().unwrap(),
                m.viscosity.value.constant().unwrap(),
            );
            assert!(r > 0.0 && mu > 0.0, "{name}: {r} {mu}");
        }
        // air at 300 K: ρ = 1.1766, μ = 1.846e-5 (White, Table A.2: 1.846e-5)
        let a = at("air", 300.0);
        let rho = a.density.value.constant().unwrap();
        let mu = a.viscosity.value.constant().unwrap();
        assert!((rho - 1.1766).abs() < 1e-3, "{rho}");
        assert!((mu - 1.846e-5).abs() < 0.005e-5, "{mu}");
        // water at the table points is exact
        let w = at("water", 293.15);
        assert_eq!(w.density.value.constant(), Some(998.21));
        assert_eq!(w.viscosity.value.constant(), Some(1.0016e-3));
    }

    #[test]
    fn steffen_is_monotone_and_exact_at_points() {
        let x = [0.0, 1.0, 2.0, 3.0, 4.0];
        let y = [0.0, 0.1, 5.0, 5.1, 5.2];
        for (xi, yi) in x.iter().zip(&y) {
            assert_eq!(steffen(&x, &y, *xi), *yi);
        }
        let mut prev = -1.0;
        for i in 0..=400 {
            let v = steffen(&x, &y, 4.0 * i as f64 / 400.0);
            assert!(v >= prev - 1e-15, "not monotone at {i}: {v} < {prev}");
            assert!((0.0..=5.2).contains(&v));
            prev = v;
        }
        assert!(steffen(&x, &y, f64::NAN).is_nan());
    }

    #[test]
    fn range_warning_and_base() {
        let mut w = Vec::new();
        let mut case = BTreeMap::new();
        case.insert(
            "thin-oil".to_string(),
            MaterialSection {
                base: Some("engine-oil".into()),
                viscosity: Some(Property::Number(0.05)),
                ..Default::default()
            },
        );
        let p = resolve_material(
            Some("thin-oil"),
            &CaseMaterials::from(case),
            None,
            None,
            &Env::default(),
            Conditions {
                temperature: 260.0,
                pressure: DEFAULT_PRESSURE,
            },
            &mut w,
        )
        .unwrap();
        assert_eq!(p.viscosity.value.constant(), Some(0.05));
        assert_eq!(p.density.value.constant(), Some(899.1));
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].code, diag::MAT_RANGE);
        // a cycle
        let mut case = BTreeMap::new();
        case.insert(
            "a".to_string(),
            MaterialSection {
                base: Some("b".into()),
                ..Default::default()
            },
        );
        case.insert(
            "b".to_string(),
            MaterialSection {
                base: Some("a".into()),
                ..Default::default()
            },
        );
        let e = resolve_material(
            Some("a"),
            &CaseMaterials::from(case),
            None,
            None,
            &Env::default(),
            Conditions::default(),
            &mut w,
        )
        .unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MAT_BASE);
        // an unknown name with a hint
        let e = resolve_material(
            Some("watr"),
            &CaseMaterials::default(),
            None,
            None,
            &Env::default(),
            Conditions::default(),
            &mut w,
        )
        .unwrap_err();
        assert!(format!("{e:#}").contains("did you mean `water`"));
    }

    #[test]
    fn a_varying_density_is_refused() {
        let mut w = Vec::new();
        let e = resolve_material(
            None,
            &CaseMaterials::default(),
            Some(&Scalar::Expr("1000 + x".into())),
            None,
            &Env::default(),
            Conditions::default(),
            &mut w,
        )
        .unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MAT_FORM);
        // T in an expression is the case temperature
        let p = resolve_material(
            None,
            &CaseMaterials::default(),
            None,
            Some(&Scalar::Expr("1e-3 * T / 293.15".into())),
            &Env::default(),
            Conditions::default(),
            &mut w,
        )
        .unwrap();
        assert_eq!(p.viscosity.value.constant(), Some(1e-3));
    }

    fn lib_with(name: &str, m: MaterialSection) -> CaseMaterials {
        let mut named = BTreeMap::new();
        named.insert(name.to_string(), m);
        CaseMaterials::from(named)
    }

    fn th(lib: &CaseMaterials, name: &str, t: f64) -> Result<Thermal> {
        let mut w = Vec::new();
        thermal_properties(
            name,
            lib,
            &Env::default(),
            Conditions {
                temperature: t,
                pressure: DEFAULT_PRESSURE,
            },
            &mut w,
        )
    }

    #[test]
    fn nasa_forms_give_cp_and_are_continuous_at_the_mid_limit() {
        // N2 (NASA-7, McBride, Zehe & Gordon 2002 / GRI-Mech 3.0 THERMO):
        // c_p/R at 300 K is 3.5 (a1 + a2 T + ...), the table value
        // c_p = 1040 J/(kg K) at 300 K (Incropera Table A.4: 1.039 kJ/(kg K)).
        let lo = vec![
            3.298_677,
            1.408_240_4e-3,
            -3.963_222e-6,
            5.641_515e-9,
            -2.444_854e-12,
            -1020.8999,
            3.950_372,
        ];
        let hi = vec![
            2.926_640,
            1.487_976_8e-3,
            -5.684_761_e-7,
            1.009_704_e-10,
            -6.753_351e-15,
            -922.7977,
            5.980_528,
        ];
        let m = MaterialSection {
            molar_mass: Some(0.028_013_4),
            specific_heat: Some(Property::Form(Form::Nasa7 {
                temperature: vec![200.0, 1000.0, 5000.0],
                coefficients: vec![lo.clone(), hi.clone()],
                file: None,
                species: None,
                source: Some("GRI-Mech 3.0 THERMO".into()),
            })),
            ..Default::default()
        };
        let lib = lib_with("n2", m);
        let cp = |t: f64| {
            th(&lib, "n2", t)
                .unwrap()
                .specific_heat
                .unwrap()
                .value
                .constant()
                .unwrap()
        };
        assert!((cp(300.0) - 1039.0).abs() < 3.0, "{}", cp(300.0));
        // the two ranges agree at T_mid to a part in a hundred
        let (a, b) = (cp(1000.0 - 1e-9), cp(1000.0 + 1e-9));
        assert!((a - b).abs() / a < 0.01, "{a} {b}");
        // NASA-9 with the same c_p/R shape: a3 only, constant c_p/R = 3.5
        let m9 = MaterialSection {
            molar_mass: Some(0.028_013_4),
            specific_heat: Some(Property::Form(Form::Nasa9 {
                temperature: vec![200.0, 6000.0],
                coefficients: vec![vec![0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
                source: None,
            })),
            ..Default::default()
        };
        let lib9 = lib_with("c", m9);
        let c = th(&lib9, "c", 400.0).unwrap().specific_heat.unwrap();
        assert!((c.value.constant().unwrap() - R_GAS / 0.028_013_4 * 3.5).abs() < 1e-9);
        // a NASA form for another property, and one without molar mass
        let bad = MaterialSection {
            density: Some(Property::Form(Form::Nasa9 {
                temperature: vec![200.0, 6000.0],
                coefficients: vec![vec![0.0; 9]],
                source: None,
            })),
            ..Default::default()
        };
        let e = th(&lib_with("b", bad), "b", 300.0).unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MAT_FORM);
    }

    #[test]
    fn chemkin_records_are_read() {
        // an N2 record in the fixed CHEMKIN columns (the NASA/Burcat
        // coefficients for 200-1000-6000 K, public data)
        let p = |x: &str| format!("{x:>15}");
        let l1 = format!(
            "{:<18}{:<6}{:<20}{}{:>10}{:>10}{:>8}{:<6}1",
            "N2", "G 8/02", "N   2", "G", "200.000", "6000.000", "1000.000", ""
        );
        let l2 = [
            "2.95257637E+00",
            "1.39690040E-03",
            "-4.92631603E-07",
            "7.86010195E-11",
            "-4.60755204E-15",
        ]
        .map(p)
        .concat()
            + &format!("{:>5}", "2");
        let l3 = [
            "-9.23948688E+02",
            "5.87188762E+00",
            "3.53100528E+00",
            "-1.23660988E-04",
            "-5.02999433E-07",
        ]
        .map(p)
        .concat()
            + &format!("{:>5}", "3");
        let l4 = [
            "2.43530612E-09",
            "-1.40881235E-12",
            "-1.04697628E+03",
            "2.96747038E+00",
        ]
        .map(p)
        .concat()
            + &format!("{:>20}", "4");
        let rec = format!("THERMO\n   300.000  1000.000  5000.000\n{l1}\n{l2}\n{l3}\n{l4}\nEND\n");
        let rec = rec.as_str();
        let (t, c) = nasa7_from_chemkin(rec, "N2").unwrap();
        assert_eq!(t, vec![200.0, 1000.0, 6000.0]);
        assert_eq!(c.len(), 2);
        assert_eq!(c[0][0], 3.53100528);
        assert_eq!(c[1][0], 2.95257637);
        let cpr = nasa_cp_over_r(&t, &c, false, 300.0);
        assert!((cpr - 3.5).abs() < 0.01, "{cpr}");
        assert!(nasa7_from_chemkin(rec, "O2").is_err());
    }

    #[test]
    fn orthotropic_k_builds_the_rotated_tensor() {
        let m = MaterialSection {
            kind: Some(MaterialKind::Solid),
            density: Some(Property::Number(7650.0)),
            specific_heat: Some(Property::Number(460.0)),
            conductivity: Some(Conductivity::Orthotropic(Box::new(Orthotropic {
                form: OrthotropicForm::Orthotropic,
                values: [
                    Property::Number(25.0),
                    Property::Number(25.0),
                    Property::Number(1.0),
                ],
                // the sheets lie in a plane tilted 30 degrees about x
                axes: Some([[1.0, 0.0, 0.0], [0.0, 0.866_025_403_784_438_6, 0.5]]),
                source: None,
            }))),
            ..Default::default()
        };
        let lib = lib_with("core", m.clone());
        let k = th(&lib, "core", 300.0).unwrap().conductivity.unwrap();
        let t = k.tensor();
        // e3 = e1 x e2 = (0, -0.5, 0.866): K = 25 I - 24 e3 e3^T
        assert!((t[0][0] - 25.0).abs() < 1e-12);
        assert!(
            (t[1][1] - (25.0 - 24.0 * 0.25)).abs() < 1e-12,
            "{}",
            t[1][1]
        );
        assert!((t[2][2] - (25.0 - 24.0 * 0.75)).abs() < 1e-12);
        assert!((t[1][2] - (24.0 * 0.5 * 0.866_025_403_784_438_6)).abs() < 1e-12);
        assert_eq!(t[1][2], t[2][1]);
        // the trace is the sum of the principal values
        assert!((t[0][0] + t[1][1] + t[2][2] - 51.0).abs() < 1e-12);
        // MAT008: no axes, parallel axes
        let mut no = m.clone();
        if let Some(Conductivity::Orthotropic(o)) = no.conductivity.as_mut() {
            o.axes = None;
        }
        let e = th(&lib_with("core", no), "core", 300.0).unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MAT_AXES);
        let mut par = m;
        if let Some(Conductivity::Orthotropic(o)) = par.conductivity.as_mut() {
            o.axes = Some([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        }
        let e = th(&lib_with("core", par), "core", 300.0).unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MAT_AXES);
    }

    #[test]
    fn a_solid_needs_its_thermal_data_and_a_fluid_when_energy_is_on() {
        let solid = MaterialSection {
            kind: Some(MaterialKind::Solid),
            density: Some(Property::Number(2700.0)),
            ..Default::default()
        };
        let e = th(&lib_with("s", solid), "s", 300.0).unwrap_err();
        let d = diag::of_error(&e);
        assert_eq!(d.code, diag::MAT_SOLID);
        assert!(d.message.contains("specific_heat"), "{}", d.message);
        // a fluid without c_p and k is fine for the flow, MAT007 for the energy equation
        let t = th(&CaseMaterials::default(), "engine-oil", 300.0).unwrap();
        let e = t.require().unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MAT_THERMAL);
    }

    #[test]
    fn user_database_files_are_searched_in_order_before_the_built_in_one() {
        let dir = std::env::temp_dir().join(format!("foehn-mat-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let a = dir.join("a.toml");
        std::fs::write(
            &a,
            r#"format = "foehncraft.materials"
format_version = "1"
[materials.my-oil]
kind = "liquid"
density = 850.0
viscosity = 0.04
specific_heat = 1900.0
conductivity = 0.14
[materials.water]
base = "water"
description = "water with a different c_p"
specific_heat = 4000.0
"#,
        )
        .unwrap();
        let mut lib = CaseMaterials {
            include: vec![a.clone()],
            ..Default::default()
        };
        lib.load().unwrap();
        let t = th(&lib, "my-oil", 300.0).unwrap();
        assert_eq!(t.specific_heat.unwrap().value.constant(), Some(1900.0));
        t_require_ok(&th(&lib, "my-oil", 300.0).unwrap());
        // a file entry of a database name changes it (base of the same name)
        let w = th(&lib, "water", 300.0).unwrap();
        assert_eq!(w.specific_heat.unwrap().value.constant(), Some(4000.0));
        let p = resolve_material(
            Some("water"),
            &lib,
            None,
            None,
            &Env::default(),
            Conditions::default(),
            &mut Vec::new(),
        )
        .unwrap();
        assert_eq!(p.density.value.constant(), Some(998.21));
        // the case's own material wins over the file's
        let mut lib2 = lib.clone();
        lib2.named.insert(
            "my-oil".into(),
            MaterialSection {
                base: Some("my-oil".into()),
                specific_heat: Some(Property::Number(2000.0)),
                ..Default::default()
            },
        );
        let t2 = th(&lib2, "my-oil", 300.0).unwrap();
        assert_eq!(t2.specific_heat.unwrap().value.constant(), Some(2000.0));
        assert_eq!(t2.density.unwrap().value.constant(), Some(850.0));
        // MAT009: a missing file, a wrong format
        let mut bad = CaseMaterials {
            include: vec![dir.join("nope.toml")],
            ..Default::default()
        };
        assert_eq!(
            diag::of_error(&bad.load().unwrap_err()).code,
            diag::MAT_FILE
        );
        let wf = dir.join("wf.toml");
        std::fs::write(
            &wf,
            "format = \"other\"\nformat_version = \"1\"\n[materials]\n",
        )
        .unwrap();
        bad.include = vec![wf];
        assert_eq!(
            diag::of_error(&bad.load().unwrap_err()).code,
            diag::MAT_FILE
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    fn t_require_ok(t: &Thermal) {
        t.require().unwrap();
    }

    #[test]
    fn nasa7_and_tables_read_their_data_files() {
        let dir = std::env::temp_dir().join(format!("foehn-mat2-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("k.csv"), "# T, k\nT,k\n300, 15.0\n400, 17.0\n").unwrap();
        let f = dir.join("m.toml");
        std::fs::write(
            &f,
            r#"format = "foehncraft.materials"
format_version = "1"
[materials.steel-x]
kind = "solid"
density = 7900.0
specific_heat = 500.0
[materials.steel-x.conductivity]
form = "table"
file = "k.csv"
"#,
        )
        .unwrap();
        let mut lib = CaseMaterials {
            include: vec![f],
            ..Default::default()
        };
        lib.load().unwrap();
        let k = th(&lib, "steel-x", 350.0).unwrap().conductivity.unwrap();
        assert!((k.tensor()[0][0] - 16.0).abs() < 1e-12);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn water_and_air_have_specific_heat_and_conductivity() {
        // water at 293.15 K (IAPWS-95): c_p = 4184.05 J/(kg K); k = 0.598 W/(m K)
        let w = th(&CaseMaterials::default(), "water", 293.15).unwrap();
        let cp = w.specific_heat.unwrap().value.constant().unwrap();
        assert!((cp - 4184.05).abs() < 0.5, "{cp}");
        let k = w.conductivity.unwrap().tensor()[0][0];
        assert!((k - 0.598).abs() < 0.002, "{k}");
        // air at 300 K: c_p = 1006.4 J/(kg K), k = 0.0264 W/(m K)
        let a = th(&CaseMaterials::default(), "air", 300.0).unwrap();
        let cp = a.specific_heat.unwrap().value.constant().unwrap();
        assert!((cp - 1006.4).abs() < 0.5, "{cp}");
        let k = a.conductivity.unwrap().tensor()[0][0];
        assert!((k - 0.02638).abs() < 2e-5, "{k}");
        // the material is a gas: no solid checks; the tables are not
        // clamped inside their range
        let mut warn = Vec::new();
        th_warn("water", 293.15, &mut warn);
        assert!(warn.is_empty(), "{warn:?}");
    }

    fn th_warn(name: &str, t: f64, w: &mut Vec<Diagnostic>) {
        thermal_properties(
            name,
            &CaseMaterials::default(),
            &Env::default(),
            Conditions {
                temperature: t,
                pressure: DEFAULT_PRESSURE,
            },
            w,
        )
        .unwrap();
    }
}
