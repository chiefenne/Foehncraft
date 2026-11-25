//! Species transport of an ideal-gas mixture (steam furnace stage SF1).
//!
//! N species with mass fractions Y_i: N − 1 of them are transported on
//! the generic scalar equation of D116 with the density of the flow,
//!
//! ρ ∂Y_i/∂t + ∇·(ṁ Y_i) − Y_i ∇·ṁ = ∇·((ρ D_i + μ_t/Sc_t) ∇Y_i),
//!
//! and the last one, the *closure species* (in air N₂), is Y_N = 1 − Σ Y_i,
//! so Σ Y = 1 holds in every cell to round-off. The face diffusivity is
//! ρ D_i = μ/Sc_i (`diffusion = "schmidt"`) or λ/(c_p Le_i)
//! (`diffusion = "lewis"`); the turbulent part is μ_t/Sc_t. The diffusion
//! flux of the closure species is minus the sum of the others, so the mass
//! of the mixture is conserved. Convection is bounded by the vertex
//! limiter (MLP, D237) by default; the mass fractions are never clipped
//! after a solve.
//!
//! With `[energy] density = "variable"` the density follows the mixture,
//! ρ = p₀ M/(R T) with 1/M = Σ Y_i/M_i (stage 7c), and with
//! `properties = "mixture"` (the default) the energy equation and the flow
//! take the mixture's c_p = Σ Y_i c_p,i (NASA-7), λ (Mathur et al. 1967)
//! and μ (Wilke 1950) from kinetic theory
//! ([`foehncraft_kernels::species`]). The species data are those of
//! GRI-Mech 3.0 (53 species).
//!
//! The ledger (D119) gets one account per species and one per chemical
//! element: per boundary the convected and the diffusive mass flow, and
//! the continuity defect −Σ Y_P (∇·ṁ)_P V_P that the non-conservative form
//! carries, so that each account closes to the solver tolerance.
//!
//! Case file:
//!
//! ```toml
//! [species]
//! names = ["CH4", "O2", "N2"]        # the last one is the closure species
//! initial = { mole_fractions = { O2 = 0.21, N2 = 0.79 } }
//! diffusion = "schmidt"             # or "lewis"
//! schmidt = 0.7                     # or one per species: { CH4 = 0.71 }
//! turbulent_schmidt = 0.7
//!
//! [boundary.fuel.species]
//! mole_fractions = { CH4 = 1.0 }
//! ```

use crate::diag::{self, Diagnostic};
use crate::expr::Scalar;
use crate::scalar::{ScalarBoundary, ScalarSection, ThermalWallFunction};
use anyhow::Result;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_kernels::species as ks;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::OnceLock;

/// The prefix of the scalar field of a transported species (`Y_CH4`).
pub const PREFIX: &str = "Y_";

/// The scalar field name of species `s`.
pub fn field_name(s: &str) -> String {
    format!("{PREFIX}{s}")
}

/// The species data file (GRI-Mech 3.0, written by
/// `tools/data/gri30_species.py`).
pub const DATA: &str = include_str!("../data/species-gri30.toml");

/// One species of the data file.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Record {
    /// Atoms per molecule by element.
    pub composition: BTreeMap<String, f64>,
    /// Molar mass (kg/mol).
    pub molar_mass: f64,
    /// The limits of the two NASA-7 ranges (K).
    pub temperature: Vec<f64>,
    /// The NASA-7 coefficients of each range.
    pub coefficients: Vec<Vec<f64>>,
    /// `atom`, `linear` or `nonlinear`.
    pub geometry: String,
    /// ε/k_B (K).
    pub well_depth: f64,
    /// σ (Å).
    pub diameter: f64,
    /// Dipole moment (Debye).
    #[serde(default)]
    pub dipole: f64,
    /// Polarizability (Å³).
    #[serde(default)]
    pub polarizability: f64,
    /// Rotational relaxation number at 298 K.
    #[serde(default)]
    pub rotational_relaxation: f64,
}

#[derive(Deserialize)]
struct DataFile {
    species: BTreeMap<String, Record>,
}

/// The species the solver knows, by name.
pub fn database() -> &'static BTreeMap<String, Record> {
    static DB: OnceLock<BTreeMap<String, Record>> = OnceLock::new();
    DB.get_or_init(|| {
        // the file is part of the build; a test reads it (an error here is
        // a broken data file, reported by `Resolved::of`)
        toml::from_str::<DataFile>(DATA)
            .map(|d| d.species)
            .unwrap_or_default()
    })
}

impl Record {
    /// The kernel form of the species.
    pub fn kernel(&self) -> ks::Species {
        let c = |i: usize| -> [f64; 7] {
            let mut a = [0.0; 7];
            for (k, v) in self.coefficients[i].iter().take(7).enumerate() {
                a[k] = *v;
            }
            a
        };
        ks::Species {
            molar_mass: self.molar_mass,
            t_range: [
                self.temperature[0],
                self.temperature[1],
                self.temperature[2],
            ],
            low: c(0),
            high: c(1),
            well_depth: self.well_depth,
            diameter: self.diameter,
            dipole: self.dipole,
            geometry: match self.geometry.as_str() {
                "atom" => 0,
                "linear" => 1,
                _ => 2,
            },
        }
    }
}

/// A composition: mole fractions (numbers) or mass fractions (numbers or
/// expressions of x, y, z and, in a transient case, t). A species left out
/// has the fraction 0; the closure species may be left out and takes the
/// rest.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Composition {
    /// Mole fractions x_i by species.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mole_fractions: Option<BTreeMap<String, f64>>,
    /// Mass fractions Y_i by species: numbers or expressions (D53).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mass_fractions: Option<BTreeMap<String, Scalar>>,
}

/// `[species] diffusion`: the molecular diffusivity of the species.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum DiffusionModel {
    /// ρ D_i = μ/Sc_i.
    #[default]
    Schmidt,
    /// ρ D_i = λ/(c_p Le_i).
    Lewis,
}

/// `[species] properties`: where c_p, λ and μ come from.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum PropertyModel {
    /// The mixture from kinetic theory: c_p = Σ Y_i c_p,i (NASA-7), μ by
    /// Wilke (1950), λ by Mathur, Tondon & Saxena (1967).
    #[default]
    Mixture,
    /// The fluid's own c_p, λ and μ (only the molar mass, and with it the
    /// density, follows the mixture).
    Fluid,
}

/// A number for every species, or one per species.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum PerSpecies {
    /// The same for all.
    All(f64),
    /// By species; a species left out takes 0.7 (Schmidt) or 1 (Lewis).
    Each(BTreeMap<String, f64>),
}

impl PerSpecies {
    /// The value of species `s` (`default` where a table leaves it out).
    pub fn of(&self, s: &str, default: f64) -> f64 {
        match self {
            Self::All(v) => *v,
            Self::Each(m) => m.get(s).copied().unwrap_or(default),
        }
    }
}

fn default_sc() -> PerSpecies {
    PerSpecies::All(0.7)
}

fn default_le() -> PerSpecies {
    PerSpecies::All(1.0)
}

fn default_sct() -> f64 {
    0.7
}

/// `[species]`: the species of an ideal-gas mixture (SF1).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SpeciesSection {
    /// The species, names of the GRI-Mech 3.0 set (`CH4`, `O2`, `N2`,
    /// `H2O`, `CO2`, `CO`, `H2`, `AR`, …).
    pub names: Vec<String>,
    /// The closure species, Y = 1 − Σ of the others (default: the last
    /// name).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub closure: Option<String>,
    /// The start composition.
    pub initial: Composition,
    /// The molecular diffusivity model.
    #[serde(default)]
    pub diffusion: DiffusionModel,
    /// Schmidt numbers Sc_i = μ/(ρ D_i) of `diffusion = "schmidt"`.
    #[serde(default = "default_sc")]
    pub schmidt: PerSpecies,
    /// Lewis numbers Le_i = λ/(ρ c_p D_i) of `diffusion = "lewis"`.
    #[serde(default = "default_le")]
    pub lewis: PerSpecies,
    /// Turbulent Schmidt number Sc_t (with a turbulence model).
    #[serde(default = "default_sct")]
    pub turbulent_schmidt: f64,
    /// Where c_p, λ and μ come from.
    #[serde(default)]
    pub properties: PropertyModel,
    /// The exact mass fractions of transported species, for the
    /// verification (expressions of x, y, z and, in a transient case, t):
    /// the report gives the L2 and largest error as for a scalar.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub exact: BTreeMap<String, Scalar>,
}

/// Moves the `[boundary.<name>.species]` sub-tables out of the boundary
/// tables (as the scalar sub-tables of D116).
pub fn extract_boundary_tables(v: &mut toml::Table) -> Result<BTreeMap<String, Composition>> {
    let mut out = BTreeMap::new();
    let Some(toml::Value::Table(bs)) = v.get_mut("boundary") else {
        return Ok(out);
    };
    for (name, b) in bs.iter_mut() {
        let toml::Value::Table(t) = b else { continue };
        let Some(s) = t.remove("species") else {
            continue;
        };
        let c: Composition = s.try_into().map_err(|e: toml::de::Error| {
            diag::err(
                Diagnostic::new(
                    diag::SPECIES,
                    format!("[boundary.{name}.species]: {}", e.message()),
                )
                .key(format!("boundary.{name}.species")),
            )
        })?;
        out.insert(name.clone(), c);
    }
    Ok(out)
}

fn bad(key: impl Into<String>, msg: impl Into<String>) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::SPECIES, msg).key(key))
}

/// The species of a case, resolved: the closure species last.
#[derive(Debug, Clone, PartialEq)]
pub struct Resolved {
    /// Names, the closure species last.
    pub names: Vec<String>,
    /// The kernel data, in the order of `names`.
    pub data: Vec<ks::Species>,
    /// The elements of the mixture and the mass of each per kg of each
    /// species (kg/kg), in the order of `names`.
    pub elements: Vec<(String, Vec<f64>)>,
    /// Sc_i or Le_i per species (by the diffusion model).
    pub number: Vec<f64>,
    /// The section.
    pub section: SpeciesSection,
}

/// Standard atomic weights (kg/mol), as in GRI-Mech 3.0.
pub fn atomic_weight(e: &str) -> f64 {
    match e {
        "H" => 1.00794e-3,
        "C" => 12.0107e-3,
        "O" => 15.9994e-3,
        "N" => 14.0067e-3,
        "Ar" => 39.948e-3,
        _ => f64::NAN,
    }
}

impl Resolved {
    /// The species of `s`, checked (SPE001).
    pub fn of(s: &SpeciesSection) -> Result<Self> {
        let db = database();
        if db.is_empty() {
            return Err(bad("species", "the species data file could not be read"));
        }
        if s.names.len() < 2 {
            return Err(bad(
                "species.names",
                "[species] names needs at least two species (a single gas is a [fluid])",
            ));
        }
        if s.names.len() > ks::MAX_SPECIES {
            return Err(bad(
                "species.names",
                format!("[species] names has more than {} species", ks::MAX_SPECIES),
            ));
        }
        let mut seen = std::collections::BTreeSet::new();
        for n in &s.names {
            if !seen.insert(n) {
                return Err(bad(
                    "species.names",
                    format!("species `{n}` is named twice"),
                ));
            }
            if !db.contains_key(n) {
                let known: Vec<&str> = db.keys().map(String::as_str).collect();
                return Err(bad(
                    "species.names",
                    format!(
                        "species `{n}` is not in the species data (GRI-Mech 3.0); known: {}",
                        known.join(", ")
                    ),
                ));
            }
        }
        let closure = s
            .closure
            .clone()
            .unwrap_or_else(|| s.names[s.names.len() - 1].clone());
        if !s.names.contains(&closure) {
            return Err(bad(
                "species.closure",
                format!("the closure species `{closure}` is not in [species] names"),
            ));
        }
        let mut names: Vec<String> = s.names.iter().filter(|n| **n != closure).cloned().collect();
        names.push(closure);
        let recs: Vec<&Record> = names.iter().map(|n| &db[n]).collect();
        let mut els: Vec<String> = Vec::new();
        for r in &recs {
            for e in r.composition.keys() {
                if !els.contains(e) {
                    els.push(e.clone());
                }
            }
        }
        els.sort();
        let elements = els
            .iter()
            .map(|e| {
                let w = atomic_weight(e);
                (
                    e.clone(),
                    recs.iter()
                        .map(|r| r.composition.get(e).copied().unwrap_or(0.0) * w / r.molar_mass)
                        .collect(),
                )
            })
            .collect();
        let (model, default) = match s.diffusion {
            DiffusionModel::Schmidt => (&s.schmidt, 0.7),
            DiffusionModel::Lewis => (&s.lewis, 1.0),
        };
        let key = match s.diffusion {
            DiffusionModel::Schmidt => "species.schmidt",
            DiffusionModel::Lewis => "species.lewis",
        };
        if let PerSpecies::Each(m) = model {
            for k in m.keys() {
                if !s.names.contains(k) {
                    return Err(bad(key, format!("`{k}` is not in [species] names")));
                }
            }
        }
        for k in s.exact.keys() {
            if !names[..names.len() - 1].contains(k) {
                return Err(bad(
                    "species.exact",
                    format!("[species] exact: `{k}` is not a transported species"),
                ));
            }
        }
        let number: Vec<f64> = names.iter().map(|n| model.of(n, default)).collect();
        if let Some(v) = number.iter().find(|v| !(**v > 0.0 && v.is_finite())) {
            return Err(bad(key, format!("{key} must be > 0, not {v}")));
        }
        if !(s.turbulent_schmidt > 0.0 && s.turbulent_schmidt.is_finite()) {
            return Err(bad(
                "species.turbulent_schmidt",
                "[species] turbulent_schmidt must be > 0",
            ));
        }
        Ok(Self {
            data: recs.iter().map(|r| r.kernel()).collect(),
            names,
            elements,
            number,
            section: s.clone(),
        })
    }

    /// The number of species.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Never empty (at least two species).
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// The transported species (all but the closure).
    pub fn transported(&self) -> &[String] {
        &self.names[..self.names.len() - 1]
    }

    /// The mass fractions of the transported species of a composition
    /// (numbers or expressions); `key` names it in errors.
    pub fn mass_fractions(&self, c: &Composition, key: &str) -> Result<Vec<Scalar>> {
        let n = self.len();
        let closure = &self.names[n - 1];
        let check_names = |m: &mut dyn Iterator<Item = &String>| -> Result<()> {
            for k in m {
                if !self.names.contains(k) {
                    return Err(bad(key, format!("{key}: `{k}` is not in [species] names")));
                }
            }
            Ok(())
        };
        match (&c.mole_fractions, &c.mass_fractions) {
            (Some(_), Some(_)) | (None, None) => Err(bad(
                key,
                format!("{key}: give one of `mole_fractions` and `mass_fractions`"),
            )),
            (Some(x), None) => {
                check_names(&mut x.keys())?;
                let mut xs: Vec<f64> = self
                    .names
                    .iter()
                    .map(|s| x.get(s).copied().unwrap_or(0.0))
                    .collect();
                if let Some((s, v)) = x.iter().find(|(_, v)| !(**v >= 0.0 && **v <= 1.0)) {
                    return Err(bad(key, format!("{key}: the mole fraction of {s} is {v}")));
                }
                let rest: f64 = xs[..n - 1].iter().sum();
                close(&mut xs[n - 1], rest, x.contains_key(closure), key, closure)?;
                let mm: f64 = xs
                    .iter()
                    .zip(&self.data)
                    .map(|(x, d)| x * d.molar_mass)
                    .sum();
                Ok(xs[..n - 1]
                    .iter()
                    .zip(&self.data)
                    .map(|(x, d)| Scalar::Number(x * d.molar_mass / mm))
                    .collect())
            }
            (None, Some(y)) => {
                check_names(&mut y.keys())?;
                let mut out = Vec::with_capacity(n - 1);
                let mut sum = 0.0;
                let mut numbers = true;
                for s in self.transported() {
                    let v = y.get(s).cloned().unwrap_or(Scalar::Number(0.0));
                    match v.number() {
                        Some(x) if !(0.0..=1.0).contains(&x) => {
                            return Err(bad(
                                key,
                                format!("{key}: the mass fraction of {s} is {x}"),
                            ));
                        }
                        Some(x) => sum += x,
                        None => numbers = false,
                    }
                    out.push(v);
                }
                if numbers {
                    let mut last = y.get(closure).and_then(Scalar::number).unwrap_or(0.0);
                    close(&mut last, sum, y.contains_key(closure), key, closure)?;
                } else if y.contains_key(closure) {
                    return Err(bad(
                        key,
                        format!(
                            "{key}: with expressions leave out the closure species {closure}; \
                             it takes 1 − Σ of the others"
                        ),
                    ));
                }
                Ok(out)
            }
        }
    }
}

/// The closure fraction: `v` (given or not) against 1 − `rest`.
fn close(v: &mut f64, rest: f64, given: bool, key: &str, closure: &str) -> Result<()> {
    let left = 1.0 - rest;
    if given {
        if (*v + rest - 1.0).abs() > 1e-6 {
            return Err(bad(
                key,
                format!("{key}: the fractions sum to {}, not 1", *v + rest),
            ));
        }
    } else if left < -1e-12 {
        return Err(bad(
            key,
            format!(
                "{key}: the fractions of the others sum to {rest} > 1, so the closure species \
                 {closure} would be negative"
            ),
        ));
    }
    *v = left.max(0.0);
    Ok(())
}

/// The scalar sections and boundary sub-tables of a case with its species
/// added as scalars `Y_<name>` (the case's own scalars unchanged).
#[allow(clippy::type_complexity)]
pub fn scalar_sections(
    case: &crate::case::FlowCase,
) -> Result<(
    BTreeMap<String, ScalarSection>,
    BTreeMap<String, BTreeMap<String, ScalarBoundary>>,
)> {
    let mut sc = case.scalars.clone();
    let mut tb = case.boundary_scalar.clone();
    let Some(s) = &case.species else {
        return Ok((sc, tb));
    };
    let r = Resolved::of(s)?;
    let ini = r.mass_fractions(&s.initial, "species.initial")?;
    let mut bvals = BTreeMap::new();
    for (b, c) in &case.boundary_species {
        if !case.boundary.contains_key(b) {
            return Err(bad(
                format!("boundary.{b}.species"),
                format!("[boundary.{b}.species]: there is no boundary `{b}`"),
            ));
        }
        bvals.insert(
            b.clone(),
            r.mass_fractions(c, &format!("boundary.{b}.species"))?,
        );
    }
    for (i, n) in r.transported().iter().enumerate() {
        let f = field_name(n);
        if sc.contains_key(&f) {
            return Err(bad(
                format!("scalars.{f}"),
                format!("[scalars.{f}] has the name of the mass fraction of species {n}"),
            ));
        }
        sc.insert(
            f.clone(),
            ScalarSection {
                // the face diffusivity follows the mixture every outer
                // iteration; this is only its start
                diffusivity: Scalar::Number(1e-5),
                initial: ini[i].clone(),
                source: Scalar::Number(0.0),
                bounds: Some([0.0, 1.0]),
                exact: s.exact.get(n).cloned(),
                turbulent_prandtl: s.turbulent_schmidt,
                turbulent_prandtl_model: crate::scalar::PrandtlModel::Constant,
                wall_function: ThermalWallFunction::None,
            },
        );
        for (b, v) in &bvals {
            tb.entry(b.clone()).or_default().insert(
                f.clone(),
                ScalarBoundary {
                    value: Some(v[i].clone()),
                    ..Default::default()
                },
            );
        }
    }
    Ok((sc, tb))
}

/// The checks of `[species]` (SPE001, SPE002).
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let Some(s) = &case.species else {
        if let Some(b) = case.boundary_species.keys().next() {
            return Err(bad(
                format!("boundary.{b}.species"),
                format!("[boundary.{b}.species] without the table [species]"),
            ));
        }
        return Ok(());
    };
    scalar_sections(case)?;
    let unsupported = |key: &str, what: &str| {
        Err(diag::err(
            Diagnostic::new(
                diag::SPECIES_UNSUPPORTED,
                format!("[species] {what} (comes with a later furnace stage)"),
            )
            .key(key.to_string()),
        ))
    };
    if case.fluid_view.is_some() || !case.regions.is_empty() {
        return unsupported(
            "regions",
            "with regions (solids, several fluids, porous regions): species run in one \
             fluid region; conjugate walls come with SF5",
        );
    }
    // every inlet needs its composition
    for (b, bc) in &case.boundary {
        use crate::bc::BoundaryCondition as B;
        let inlet = matches!(
            bc,
            B::VelocityInlet { .. }
                | B::MassFlowInlet { .. }
                | B::StaticPressureInlet { .. }
                | B::TotalPressureInlet { .. }
        );
        if inlet && !case.boundary_species.contains_key(b) {
            return Err(bad(
                format!("boundary.{b}.species"),
                format!(
                    "the inlet `{b}` needs its composition: add [boundary.{b}.species] with \
                     `mole_fractions` or `mass_fractions`"
                ),
            ));
        }
    }
    let variable = crate::density::wanted(case);
    if s.properties == PropertyModel::Mixture && case.energy.is_none() {
        return Err(bad(
            "species.properties",
            "[species] properties = \"mixture\" needs the energy equation ([energy]) for the \
             temperature of the properties; without it set properties = \"fluid\" (an \
             isothermal mixture with the fluid's density and viscosity)",
        ));
    }
    if s.diffusion == DiffusionModel::Lewis && case.energy.is_none() {
        return Err(bad(
            "species.diffusion",
            "[species] diffusion = \"lewis\" needs the energy equation ([energy]) for λ and c_p",
        ));
    }
    if variable {
        let m = fluid_molar_mass(case)?;
        if m.is_none() {
            return Err(bad(
                "fluid.material",
                "the variable density of a mixture follows its molar mass: give the fluid a \
                 material with `molar_mass` and `density = { form = \"ideal-gas\" }`",
            ));
        }
    }
    Ok(())
}

/// The molar mass of the fluid's ideal-gas density form, if it has one.
pub fn fluid_molar_mass(case: &crate::case::FlowCase) -> Result<Option<f64>> {
    let Some(n) = &case.fluid.material else {
        return Ok(None);
    };
    let m = crate::materials::lookup(n, &case.materials)?;
    let ideal = m
        .density
        .as_ref()
        .is_some_and(|p| format!("{p:?}").to_ascii_lowercase().contains("ideal"));
    Ok(m.molar_mass.filter(|_| ideal))
}

/// The mixture state of the last refresh, on the host.
#[derive(Debug, Clone)]
pub struct Link {
    /// The species.
    pub res: Resolved,
    /// The index of each transported species' field in the scalars.
    pub fields: Vec<usize>,
    /// The temperature without the energy equation (K).
    pub t_fluid: f64,
    /// The fluid's μ, λ, c_p (properties = "fluid"; μ also without energy).
    pub fluid_props: (f64, f64, f64),
    /// The molar mass of the fluid's density form (kg/mol), with a variable
    /// density.
    pub m_fluid: Option<f64>,
    owner: Vec<u32>,
    neighbour: Vec<u32>,
    weight: Vec<f64>,
    /// Y per species (closure last) per cell (own and halo cells).
    pub y_cell: Vec<Vec<f64>>,
    /// Y per species per boundary face.
    pub y_bface: Vec<Vec<f64>>,
    /// The mixture per cell.
    pub cell: Vec<ks::Mixture>,
    /// The mixture per face.
    pub face: Vec<ks::Mixture>,
    /// The smallest closure fraction before it is held at zero (the sum of
    /// the transported fractions above 1 by this much).
    pub closure_min: f64,
}

impl Link {
    /// A new link for `res` with the faces of the mesh.
    pub fn new(
        res: Resolved,
        fields: Vec<usize>,
        disc: &foehncraft_fv::disc::Discretisation,
        t_fluid: f64,
        fluid_props: (f64, f64, f64),
        m_fluid: Option<f64>,
    ) -> Self {
        let ni = disc.n_interior;
        Self {
            res,
            fields,
            t_fluid,
            fluid_props,
            m_fluid,
            owner: disc.owner.clone(),
            neighbour: disc.neighbour[..ni].to_vec(),
            weight: disc.weight.clone(),
            y_cell: Vec::new(),
            y_bface: Vec::new(),
            cell: Vec::new(),
            face: Vec::new(),
            closure_min: 0.0,
        }
    }

    /// The owner cell of boundary face `fb`.
    pub fn boundary_owner(&self, fb: usize) -> usize {
        self.owner[self.neighbour.len() + fb] as usize
    }

    /// The properties model.
    pub fn properties(&self) -> PropertyModel {
        self.res.section.properties
    }

    /// Recomputes the composition and the mixture per cell and face from
    /// the transported fractions `y` per cell (own and halo cells) and
    /// per boundary face `yb`, the temperature per cell `t` and per
    /// boundary face `tb` (`None`: the fluid temperature).
    pub fn compute(&mut self, y: &[Vec<f64>], yb: &[Vec<f64>], t: Option<(&[f64], &[f64])>) {
        let n = self.res.len();
        let nc = y.first().map_or(0, Vec::len);
        let nf = self.owner.len();
        let ni = self.neighbour.len();
        let nb = nf - ni;
        let close = |v: &mut Vec<Vec<f64>>, len: usize, min: &mut f64| {
            let mut last = vec![0.0; len];
            for (c, l) in last.iter_mut().enumerate() {
                let s: f64 = (0..n - 1).map(|i| v[i][c]).sum();
                let r = 1.0 - s;
                *min = min.min(r);
                *l = r.max(0.0);
            }
            v.push(last);
        };
        let mut min = 0.0f64;
        self.y_cell = y.to_vec();
        close(&mut self.y_cell, nc, &mut min);
        self.y_bface = yb.to_vec();
        close(&mut self.y_bface, nb, &mut min);
        self.closure_min = min;
        let tf = self.t_fluid;
        let tc = |c: usize| t.map_or(tf, |(t, _)| t[c]);
        let mut yy = [0.0; ks::MAX_SPECIES];
        let sp = &self.res.data;
        self.cell = (0..nc)
            .map(|c| {
                for (v, y) in yy.iter_mut().zip(&self.y_cell) {
                    *v = y[c];
                }
                ks::mixture(sp, &yy[..n], tc(c))
            })
            .collect();
        self.face = (0..nf)
            .map(|f| {
                let tfv = if f < ni {
                    let (o, nn) = (self.owner[f] as usize, self.neighbour[f] as usize);
                    let w = self.weight[f];
                    for (v, y) in yy.iter_mut().zip(&self.y_cell) {
                        *v = w * y[o] + (1.0 - w) * y[nn];
                    }
                    w * tc(o) + (1.0 - w) * tc(nn)
                } else {
                    for (v, y) in yy.iter_mut().zip(&self.y_bface) {
                        *v = y[f - ni];
                    }
                    t.map_or(tf, |(_, tb)| tb[f - ni])
                };
                ks::mixture(sp, &yy[..n], tfv)
            })
            .collect();
    }

    /// ρ D_i per face of transported species `i` (μ/Sc_i or λ/(c_p Le_i):
    /// no density needed).
    pub fn gamma(&self, i: usize) -> Vec<f64> {
        let num = self.res.number[i];
        let fluid = self.properties() == PropertyModel::Fluid;
        let (mu_f, k_f, cp_f) = self.fluid_props;
        self.face
            .iter()
            .map(|m| match self.res.section.diffusion {
                DiffusionModel::Schmidt => {
                    let mu = if fluid { mu_f } else { m.viscosity };
                    mu / num
                }
                DiffusionModel::Lewis => {
                    let (k, cp) = if fluid {
                        (k_f, cp_f)
                    } else {
                        (m.conductivity, m.cp)
                    };
                    k / (cp * num)
                }
            })
            .collect()
    }

    /// The factor M/M_fluid of the density per cell and per boundary face,
    /// with a variable density.
    pub fn molar_factor(&self) -> Option<(Vec<f64>, Vec<f64>)> {
        let mf = self.m_fluid?;
        let ni = self.neighbour.len();
        Some((
            self.cell.iter().map(|m| m.molar_mass / mf).collect(),
            self.face[ni..].iter().map(|m| m.molar_mass / mf).collect(),
        ))
    }
}

/// The mass of species `i` per kg of mixture at every face of a boundary
/// summary: see [`accounts`].
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct SpeciesFlow {
    /// The species.
    pub species: String,
    /// Per boundary: (name, convected, diffusive) mass flow into the
    /// domain (kg/s).
    pub boundaries: Vec<(String, f64, f64)>,
    /// −Σ Y_P (∇·ṁ)_P V_P (kg/s): the continuity defect of the
    /// non-conservative form (zero when the continuity is converged).
    pub defect: f64,
    /// The storage Σ ρ V ∂Y/∂t of a time step (kg/s), 0 steady.
    pub storage: f64,
    /// The mass an evaporating spray gives the species (kg/s, SP5): the
    /// vapour S_v, the mixture the continuity source S_c; the continuity
    /// defect is then taken against S_c.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub spray: f64,
}

fn is_zero(x: &f64) -> bool {
    *x == 0.0
}

/// The ledger accounts of the species and of the elements (D119): one
/// account `species <name>` per species and one `element <E>` per
/// element, in kg.
pub fn accounts(flows: &[SpeciesFlow], res: &Resolved) -> Vec<crate::balance::Account> {
    use crate::balance::{Account, Term, TermKind};
    let book = |q: String, weights: &dyn Fn(usize) -> f64| {
        let mut a = Account::new(&q, "kg");
        if let Some(first) = flows.first() {
            for (bi, (b, _, _)) in first.boundaries.iter().enumerate() {
                let conv: f64 = flows
                    .iter()
                    .enumerate()
                    .map(|(i, s)| weights(i) * s.boundaries[bi].1)
                    .exact_sum();
                let dif: f64 = flows
                    .iter()
                    .enumerate()
                    .map(|(i, s)| weights(i) * s.boundaries[bi].2)
                    .exact_sum();
                a.terms.push(Term::boundary("convection", b, conv));
                a.terms.push(Term::boundary("diffusion", b, dif));
            }
        }
        let defect: f64 = flows
            .iter()
            .enumerate()
            .map(|(i, s)| weights(i) * s.defect)
            .exact_sum();
        a.terms
            .push(Term::volume("continuity defect", TermKind::Source, defect));
        if flows.iter().any(|s| s.spray != 0.0) {
            let spray: f64 = flows
                .iter()
                .enumerate()
                .map(|(i, s)| weights(i) * s.spray)
                .exact_sum();
            a.terms
                .push(Term::volume("spray vapour", TermKind::Source, spray));
        }
        a.storage = Some(
            flows
                .iter()
                .enumerate()
                .map(|(i, s)| weights(i) * s.storage)
                .exact_sum(),
        );
        a
    };
    let mut out = Vec::new();
    for (i, s) in flows.iter().enumerate() {
        out.push(book(format!("species {}", s.species), &|j| {
            if j == i { 1.0 } else { 0.0 }
        }));
    }
    for (e, w) in &res.elements {
        out.push(book(format!("element {e}"), &|j| w[j]));
    }
    out
}

/// The species part of the report: per species its extremes and its
/// mass flows, and the closure of each element account.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct SpeciesSummary {
    /// Per species (closure last): (name, min Y, max Y).
    pub extremes: Vec<(String, f64, f64)>,
    /// The mass flows per species.
    pub flows: Vec<SpeciesFlow>,
    /// Per element: (name, inflow kg/s, imbalance kg/s, relative).
    pub elements: Vec<(String, f64, f64, f64)>,
    /// The smallest closure fraction 1 − Σ Y_i before it is held at zero.
    pub closure_min: f64,
}

/// One log line per element of the summary.
pub fn summary_lines(s: &SpeciesSummary) -> Vec<String> {
    let mut out = Vec::new();
    for (n, lo, hi) in &s.extremes {
        out.push(format!("species {n:<6} Y in [{lo:.6e}, {hi:.6e}]"));
    }
    for (e, inflow, imb, rel) in &s.elements {
        out.push(if *inflow > 0.0 {
            format!(
                "element {e:<2} inflow {inflow:.6e} kg/s, imbalance {imb:.3e} kg/s ({rel:.2e} of the inflow)"
            )
        } else {
            format!("element {e:<2} no inflow, imbalance {imb:.3e} kg/s")
        });
    }
    out
}

/// The species of a solved state for the report and the ledger: the mass
/// flows of every species per boundary (convected as the scheme convects
/// them, diffusive as the assembly), the continuity defect and the storage
/// of the last time step; `None` without species.
pub fn summary(
    fl: &mut crate::simplec::Flow,
    mesh: &foehncraft_mesh::Mesh,
    geo: &foehncraft_mesh::Geometry,
    transient: bool,
) -> Result<Option<SpeciesSummary>> {
    let Some(sc) = fl.scalars.as_ref() else {
        return Ok(None);
    };
    let Some(link) = sc.species() else {
        return Ok(None);
    };
    let res = link.res.clone();
    let closure_min = link.closure_min;
    let data = sc.species_host(fl.dev.as_mut())?;
    let old = if transient {
        sc.species_old(fl.dev.as_mut())?
    } else {
        None
    };
    let m = fl.dev.download(fl.f.mass_flux)?;
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let nf = mesh.n_faces();
    // the net mass outflow of every own cell, Σ_f ṁ_f out
    let mut div = vec![0.0; nc];
    for (f, &mf) in m.iter().enumerate().take(nf) {
        let o = mesh.owner[f] as usize;
        if o < nc {
            div[o] += mf;
        }
        if f < ni {
            let n = mesh.neighbour[f] as usize;
            if n < nc {
                div[n] -= mf;
            }
        }
    }
    // an evaporating spray (SP5): per spray cell (cell, a, b, S_c); the
    // vapour receives a − b Y, the mixture S_c
    let spray: Vec<(usize, f64, f64, f64)> = sc
        .spray_species()
        .map(|s| {
            s.rows
                .iter()
                .map(|&(c, a, b, m, _)| (c as usize, a, b, m))
                .collect()
        })
        .unwrap_or_default();
    let vapour_field = sc.spray_species().map(|s| s.vapour.clone());
    let rho = |c: usize| fl.phys.rho_at(c);
    let n = res.len();
    let mut flows: Vec<SpeciesFlow> = Vec::with_capacity(n);
    let mut extremes = Vec::with_capacity(n);
    // the closure: the mixture's total minus the transported
    let mut total = SpeciesFlow {
        species: res.names[n - 1].clone(),
        ..Default::default()
    };
    for p in &mesh.boundaries {
        let conv: f64 = p.range().map(|f| -m[f]).exact_sum();
        total.boundaries.push((p.name.clone(), conv, 0.0));
    }
    total.defect = div.iter().exact_sum();
    if !spray.is_empty() {
        total.defect -= spray.iter().map(|r| r.3).exact_sum();
        total.spray = spray.iter().map(|r| r.3).exact_sum();
    }
    for (k, (phi, a, b, gamma, out)) in data.iter().enumerate() {
        let mut s = SpeciesFlow {
            species: res.names[k].clone(),
            ..Default::default()
        };
        for p in &mesh.boundaries {
            let (mut conv, mut dif) = (Vec::new(), Vec::new());
            for f in p.range() {
                let fb = f - ni;
                let o = mesh.owner[f] as usize;
                let pf = a[fb] * phi[o] + b[fb];
                let v = if m[f] > 0.0 {
                    out.at(geo, f, o, phi[o], pf)
                } else {
                    pf
                };
                conv.push(-m[f] * v);
                dif.push(gamma[f] * fl.disc.dgeo[f] * (pf - phi[o]));
            }
            s.boundaries.push((
                p.name.clone(),
                conv.iter().exact_sum(),
                dif.iter().exact_sum(),
            ));
        }
        s.defect = (0..nc).map(|c| phi[c] * div[c]).exact_sum();
        if !spray.is_empty() {
            s.defect -= spray.iter().map(|r| phi[r.0] * r.3).exact_sum();
            if vapour_field.as_deref() == Some(field_name(&res.names[k]).as_str()) {
                s.spray = spray.iter().map(|r| r.1 - r.2 * phi[r.0]).exact_sum();
            }
        }
        if let Some((st, levels)) = &old {
            let (o1, o2) = &levels[k];
            s.storage = (0..nc)
                .map(|c| {
                    rho(c) * geo.cell_volume[c] / st.dt
                        * (st.a[0] * phi[c] + st.a[1] * o1[c] + st.a[2] * o2[c])
                })
                .exact_sum();
        }
        let (lo, hi) = phi[..nc]
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &v| {
                (l.min(v), h.max(v))
            });
        extremes.push((res.names[k].clone(), lo, hi));
        flows.push(s);
    }
    // the closure species: total minus the sum of the others
    let mut cl = total.clone();
    for (bi, e) in cl.boundaries.iter_mut().enumerate() {
        e.1 -= flows.iter().map(|s| s.boundaries[bi].1).exact_sum();
        e.2 = -flows.iter().map(|s| s.boundaries[bi].2).exact_sum();
    }
    cl.defect -= flows.iter().map(|s| s.defect).exact_sum();
    cl.storage = -flows.iter().map(|s| s.storage).exact_sum();
    cl.spray -= flows.iter().map(|s| s.spray).exact_sum();
    let (lo, hi) = link.y_cell[n - 1][..nc]
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &v| {
            (l.min(v), h.max(v))
        });
    extremes.push((res.names[n - 1].clone(), lo, hi));
    flows.push(cl);
    let accounts = accounts(&flows, &res);
    let elements = accounts
        .iter()
        .filter_map(|a| {
            let e = a.quantity.strip_prefix("element ")?;
            let c = a.close(a.storage.unwrap_or(0.0));
            let inflow: f64 = a
                .terms
                .iter()
                .filter(|t| t.boundary.is_some() && t.value > 0.0)
                .map(|t| t.value)
                .sum();
            let rel = if inflow > 0.0 {
                c.imbalance.abs() / inflow
            } else {
                c.relative
            };
            Some((e.to_string(), inflow, c.imbalance, rel))
        })
        .collect();
    Ok(Some(SpeciesSummary {
        extremes,
        flows,
        elements,
        closure_min,
    }))
}

/// The SPE003 warning of a summary, if the closure species went negative.
pub fn closure_finding(s: &SpeciesSummary) -> Option<Diagnostic> {
    (s.closure_min < -1e-4).then(|| {
        Diagnostic::new(
            diag::SPECIES_CLOSURE,
            format!(
                "the transported mass fractions summed to 1 + {:.3e} in some cells; the closure \
                 species was held at zero there",
                -s.closure_min
            ),
        )
        .key("species")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn db(n: &str) -> ks::Species {
        database()[n].kernel()
    }

    #[test]
    fn the_data_file_has_the_53_species_of_gri_mech() {
        let d = database();
        assert_eq!(d.len(), 53);
        assert!((d["CH4"].molar_mass - 0.01604246).abs() < 1e-8);
        assert!((d["N2"].molar_mass - 0.0280134).abs() < 1e-8);
    }

    /// Pure-gas viscosities and conductivities at 300 K against measured
    /// values (Poling, Prausnitz & O'Connell 2001, appendix; Kadoya et al.
    /// 1985 for N₂): μ within 3 %, λ within 5 % (the mean of the Eucken
    /// and the modified Eucken relations, Poling et al. 2001, section 10-3).
    #[test]
    fn pure_gas_properties_match_measured_values() {
        // (species, μ measured Pa s, λ measured W/(m K)) at 300 K, 1 atm
        let cases = [
            ("N2", 1.79e-5, 0.0260),
            ("O2", 2.07e-5, 0.0267),
            ("CH4", 1.12e-5, 0.0343),
            ("AR", 2.27e-5, 0.0177),
            ("H2", 8.96e-6, 0.183),
            ("CO2", 1.50e-5, 0.0166),
        ];
        for (n, mu_m, k_m) in cases {
            let s = db(n);
            let mu = s.viscosity(300.0);
            let k = s.conductivity(300.0, mu);
            assert!(
                (mu / mu_m - 1.0).abs() < 0.03,
                "{n}: μ = {mu:e}, measured {mu_m:e}"
            );
            assert!((k / k_m - 1.0).abs() < 0.05, "{n}: λ = {k}, measured {k_m}");
        }
    }

    /// Air (x_O2 = 0.21, x_N2 = 0.78, x_Ar = 0.01) at 300 K: M = 28.97
    /// g/mol, c_p = 1005 J/(kg K), μ = 1.85e-5 Pa s, λ = 0.0263 W/(m K)
    /// (Lemmon et al. 2004, measured and correlated); c_p within 0.3 %, μ
    /// within 3 %, λ within 4 %.
    #[test]
    fn air_as_a_mixture() {
        let sp = [db("O2"), db("N2"), db("AR")];
        let x = [0.21, 0.78, 0.01];
        let mm: f64 = x.iter().zip(&sp).map(|(x, s)| x * s.molar_mass).sum();
        let y: Vec<f64> = x
            .iter()
            .zip(&sp)
            .map(|(x, s)| x * s.molar_mass / mm)
            .collect();
        let m = ks::mixture(&sp, &y, 300.0);
        assert!((m.molar_mass - 0.028_97).abs() < 2e-5, "{}", m.molar_mass);
        assert!((m.cp / 1005.0 - 1.0).abs() < 0.003, "{}", m.cp);
        assert!(
            (m.viscosity / 1.85e-5 - 1.0).abs() < 0.03,
            "{}",
            m.viscosity
        );
        assert!(
            (m.conductivity / 0.0263 - 1.0).abs() < 0.04,
            "{}",
            m.conductivity
        );
    }

    /// NASA-7: the enthalpy of formation of CH₄ at 298.15 K, −74.53 kJ/mol
    /// (Active Thermochemical Tables, Ruscic et al.), within 0.2 %, and
    /// c_p continuous at the middle temperature.
    #[test]
    fn nasa_polynomials_give_the_formation_enthalpy() {
        let s = db("CH4");
        let h = s.h_over_rt(298.15) * ks::R * 298.15;
        assert!((h / -74_530.0 - 1.0).abs() < 2e-3, "{h}");
        let t = s.t_range[1];
        let lo = s.cp_over_r(t);
        let hi = s.cp_over_r(t * (1.0 + 1e-12));
        assert!((lo - hi).abs() < 1e-3 * lo);
    }

    /// The binary diffusion coefficient of CH₄ in air at 298 K and 1 atm:
    /// measured 2.10e-5 to 2.20e-5 m²/s (Poling et al. 2001, table 11-2,
    /// CH₄–air 0.210 cm²/s; Marrero & Mason 1972), within 6 %.
    #[test]
    fn binary_diffusion_of_methane_in_nitrogen() {
        let d = ks::binary_diffusion(&db("CH4"), &db("N2"), 298.0, 101_325.0);
        assert!((d / 2.15e-5 - 1.0).abs() < 0.06, "{d:e}");
    }

    #[test]
    fn mole_fractions_become_mass_fractions_and_close() {
        let s = SpeciesSection {
            names: vec!["CH4".into(), "O2".into(), "N2".into()],
            closure: None,
            initial: Composition {
                mole_fractions: Some([("O2".to_string(), 0.21)].into_iter().collect()),
                mass_fractions: None,
            },
            diffusion: DiffusionModel::Schmidt,
            schmidt: default_sc(),
            lewis: default_le(),
            turbulent_schmidt: 0.7,
            properties: PropertyModel::Mixture,
            exact: BTreeMap::new(),
        };
        let r = Resolved::of(&s).unwrap();
        assert_eq!(r.names, ["CH4", "O2", "N2"]);
        let y = r.mass_fractions(&s.initial, "species.initial").unwrap();
        let y_o2 = y[1].number().unwrap();
        assert!((y_o2 - 0.2330).abs() < 2e-3, "{y_o2}");
        assert_eq!(y[0].number(), Some(0.0));
        // a sum above one is an error
        let c = Composition {
            mole_fractions: Some(
                [("O2".to_string(), 0.5), ("CH4".to_string(), 0.6)]
                    .into_iter()
                    .collect(),
            ),
            mass_fractions: None,
        };
        assert!(r.mass_fractions(&c, "x").is_err());
        // the element masses of each species sum to one
        for j in 0..r.len() {
            let s: f64 = r.elements.iter().map(|(_, w)| w[j]).sum();
            assert!((s - 1.0).abs() < 1e-6, "{} {s}", r.names[j]);
        }
    }
}
