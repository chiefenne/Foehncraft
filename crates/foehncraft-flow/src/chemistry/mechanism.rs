//! Reaction mechanisms (SF2): the built-in ones and the importer of the
//! YAML mechanism format of D103 (the format of the Cantera library).
//!
//! A mechanism is read once into the flat program of the kernels
//! ([`foehncraft_kernels::chemistry`]): species with their NASA-7 and
//! transport data, reactions with SI rate constants, species terms and
//! third-body efficiencies. The importer reads the documented YAML keys
//! (`units`, `phases`, `species` with `NASA7` thermo and `gas` transport,
//! `reactions` with `equation`, `rate-constant`, `type: three-body`,
//! `type: falloff` with `low-P-rate-constant`, `high-P-rate-constant`,
//! `Troe`, `efficiencies`, `default-efficiency`, `orders`,
//! `nonreactant-orders`, `duplicate`); no library is linked. A phase
//! whose species have no record in the file takes them from the built-in
//! species database of SF1 (GRI-Mech 3.0 thermo and transport), so a global
//! mechanism needs only its reactions. Forms not implemented yet
//! (pressure-dependent Arrhenius, Chebyshev, SRI fall-off, chemically
//! activated reactions, a specific collider `(+N2)`) are refused with
//! CHM001, as is a reaction that does not conserve the elements.

use crate::diag::{self, Diagnostic};
use crate::species::{self, Record};
use anyhow::Result;
use foehncraft_kernels::chemistry as kc;
use foehncraft_kernels::species as ks;
use std::collections::BTreeMap;
use yaml_rust2::{Yaml, YamlLoader};

/// The built-in mechanisms: (name, description, YAML text).
pub const BUILT_IN: &[(&str, &str, &str)] = &[
    (
        "gri30",
        "GRI-Mech 3.0 (Smith et al. 1999): 53 species, 325 reactions, natural gas with NOx",
        include_str!("../../data/mechanisms/gri30.yaml"),
    ),
    (
        "westbrook-dryer-2",
        "Two-step global CH4 mechanism of Westbrook & Dryer (1981)",
        include_str!("../../data/mechanisms/westbrook-dryer-2.yaml"),
    ),
    (
        "jones-lindstedt-4",
        "Four-step global CH4 mechanism of Jones & Lindstedt (1988)",
        include_str!("../../data/mechanisms/jones-lindstedt-4.yaml"),
    ),
];

pub(crate) fn bad(key: impl Into<String>, msg: impl Into<String>) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::CHEMISTRY, msg).key(key))
}

/// A mechanism resolved into the flat program of the kernels.
#[derive(Debug, Clone, PartialEq)]
pub struct Mechanism {
    /// Its name (the built-in name or the file).
    pub name: String,
    /// The description of the file.
    pub description: String,
    /// The species names.
    pub species: Vec<String>,
    /// The species records (composition, thermo, transport).
    pub records: Vec<Record>,
    /// The kernel form of the species.
    pub data: Vec<ks::Species>,
    /// The elements, sorted.
    pub elements: Vec<String>,
    /// Atoms of each element per molecule of each species
    /// (`atoms[e][k]`).
    pub atoms: Vec<Vec<f64>>,
    /// The reactions.
    pub reactions: Vec<kc::Reaction>,
    /// Their species terms.
    pub terms: Vec<kc::Term>,
    /// Their third-body efficiencies.
    pub efficiencies: Vec<kc::Efficiency>,
    /// The equation of each reaction, as written in the file.
    pub equations: Vec<String>,
}

impl Mechanism {
    /// The kernel view.
    pub fn program(&self) -> kc::Program<'_> {
        kc::Program {
            species: &self.data,
            reactions: &self.reactions,
            terms: &self.terms,
            efficiencies: &self.efficiencies,
        }
    }

    /// The index of species `s`.
    pub fn index(&self, s: &str) -> Option<usize> {
        self.species.iter().position(|n| n == s)
    }

    /// The number of species.
    pub fn len(&self) -> usize {
        self.species.len()
    }

    /// Never empty.
    pub fn is_empty(&self) -> bool {
        self.species.is_empty()
    }

    /// A built-in mechanism by name, or a mechanism file (YAML) by path.
    pub fn load(spec: &str, base: Option<&std::path::Path>) -> Result<Self> {
        if let Some((n, _, text)) = BUILT_IN.iter().find(|(n, _, _)| *n == spec) {
            return Self::from_yaml(text, n);
        }
        let looks_like_file =
            spec.ends_with(".yaml") || spec.ends_with(".yml") || spec.contains('/');
        if !looks_like_file {
            let known: Vec<&str> = BUILT_IN.iter().map(|b| b.0).collect();
            return Err(bad(
                "mechanism",
                format!(
                    "unknown mechanism `{spec}`: give a built-in one ({}) or a YAML mechanism file",
                    known.join(", ")
                ),
            ));
        }
        let path = match base {
            Some(b) if std::path::Path::new(spec).is_relative() => b.join(spec),
            _ => std::path::PathBuf::from(spec),
        };
        let text = std::fs::read_to_string(&path).map_err(|e| {
            bad(
                "mechanism",
                format!("the mechanism file {} cannot be read: {e}", path.display()),
            )
        })?;
        Self::from_yaml(&text, spec)
    }

    /// Reads a mechanism from the YAML text of D103.
    pub fn from_yaml(text: &str, name: &str) -> Result<Self> {
        let docs = YamlLoader::load_from_str(text)
            .map_err(|e| bad("mechanism", format!("{name}: not a YAML file: {e}")))?;
        let Some(doc) = docs.first() else {
            return Err(bad("mechanism", format!("{name}: the file is empty")));
        };
        let units = Units::of(&doc["units"], name)?;
        // the species of the (first) phase, else all records
        let records_yaml: BTreeMap<String, &Yaml> = doc["species"]
            .as_vec()
            .map(|v| {
                v.iter()
                    .filter_map(|s| Some((s["name"].as_str()?.to_string(), s)))
                    .collect()
            })
            .unwrap_or_default();
        let phase = doc["phases"].as_vec().and_then(|p| p.first());
        let names: Vec<String> = match phase.and_then(|p| p["species"].as_vec()) {
            Some(v) => v
                .iter()
                .map(|s| {
                    s.as_str().map(str::to_string).ok_or_else(|| {
                        bad(
                            "mechanism",
                            format!(
                                "{name}: the species of the phase must be a list of names \
                                 (species from other files are not read)"
                            ),
                        )
                    })
                })
                .collect::<Result<_>>()?,
            None => doc["species"]
                .as_vec()
                .map(|v| {
                    v.iter()
                        .filter_map(|s| s["name"].as_str().map(str::to_string))
                        .collect()
                })
                .unwrap_or_default(),
        };
        if names.is_empty() {
            return Err(bad("mechanism", format!("{name}: no species")));
        }
        if names.len() > ks::MAX_SPECIES {
            return Err(bad(
                "mechanism",
                format!("{name}: more than {} species", ks::MAX_SPECIES),
            ));
        }
        let db = species::database();
        let mut records = Vec::with_capacity(names.len());
        for s in &names {
            let r = match records_yaml.get(s) {
                Some(y) => record(y, s, name)?,
                None => db.get(s).cloned().ok_or_else(|| {
                    bad(
                        "mechanism",
                        format!(
                            "{name}: species `{s}` has no record in the file and is not in the \
                             built-in species database"
                        ),
                    )
                })?,
            };
            records.push(r);
        }
        let mut elements: Vec<String> = Vec::new();
        for r in &records {
            for e in r.composition.keys() {
                if !elements.contains(e) {
                    elements.push(e.clone());
                }
            }
        }
        elements.sort();
        let atoms = elements
            .iter()
            .map(|e| {
                records
                    .iter()
                    .map(|r| r.composition.get(e).copied().unwrap_or(0.0))
                    .collect()
            })
            .collect();
        let mut m = Self {
            name: name.to_string(),
            description: doc["description"].as_str().unwrap_or("").trim().to_string(),
            data: records.iter().map(Record::kernel).collect(),
            species: names,
            records,
            elements,
            atoms,
            reactions: Vec::new(),
            terms: Vec::new(),
            efficiencies: Vec::new(),
            equations: Vec::new(),
        };
        if let Some(rs) = doc["reactions"].as_vec() {
            for (i, r) in rs.iter().enumerate() {
                m.reaction(r, &units)
                    .map_err(|e| bad("mechanism", format!("{name}: reaction {}: {e}", i + 1)))?;
            }
        }
        Ok(m)
    }

    fn reaction(&mut self, r: &Yaml, units: &Units) -> std::result::Result<(), String> {
        let eq = r["equation"]
            .as_str()
            .ok_or_else(|| "no `equation`".to_string())?
            .to_string();
        let kind = r["type"].as_str().unwrap_or("elementary");
        let parsed = parse_equation(&eq)?;
        let kind = match (kind, parsed.third_body) {
            ("elementary", Third::None) => kc::ELEMENTARY,
            ("three-body", Third::M) | ("elementary", Third::M) => kc::THREE_BODY,
            ("falloff", Third::Falloff) => {
                if !r["SRI"].is_badvalue() {
                    return Err("SRI fall-off is not implemented yet".into());
                }
                if r["Troe"].is_badvalue() {
                    kc::LINDEMANN
                } else {
                    kc::TROE
                }
            }
            ("falloff", Third::Collider(c)) | ("three-body", Third::Collider(c)) => {
                return Err(format!(
                    "a specific collider `{c}` is not implemented yet; write `M` with efficiencies"
                ));
            }
            (k, t) => {
                return Err(format!(
                    "the rate form `{k}` (third body {t:?}) is not implemented yet (implemented: \
                     elementary, three-body, falloff with Lindemann or Troe)"
                ));
            }
        };
        // species terms, merged per species
        let mut terms: BTreeMap<usize, kc::Term> = BTreeMap::new();
        let mut order_of_species: Vec<usize> = Vec::new();
        let mut add = |s: &str, nu: f64, reactant: bool| -> std::result::Result<(), String> {
            let k = self
                .index(s)
                .ok_or_else(|| format!("species `{s}` is not in the phase"))?;
            let t = terms.entry(k).or_insert_with(|| {
                order_of_species.push(k);
                kc::Term {
                    species: k as u32,
                    reactant: 0.0,
                    product: 0.0,
                    order: 0.0,
                }
            });
            if reactant {
                t.reactant += nu;
                t.order += nu;
            } else {
                t.product += nu;
            }
            Ok(())
        };
        for (s, nu) in &parsed.reactants {
            add(s, *nu, true)?;
        }
        for (s, nu) in &parsed.products {
            add(s, *nu, false)?;
        }
        // explicit orders
        if let Some(h) = r["orders"].as_hash() {
            let nonreactant = r["nonreactant-orders"].as_bool().unwrap_or(false);
            for (k, v) in h {
                let s = k.as_str().ok_or("an order without a species name")?;
                let a = number(v).ok_or_else(|| format!("the order of `{s}` is not a number"))?;
                let idx = self
                    .index(s)
                    .ok_or_else(|| format!("the order of `{s}`: not a species of the phase"))?;
                match terms.get_mut(&idx) {
                    Some(t) if t.reactant > 0.0 || nonreactant => t.order = a,
                    Some(_) | None if !nonreactant => {
                        return Err(format!(
                            "an order for `{s}`, which is not a reactant: add \
                             `nonreactant-orders: true`"
                        ));
                    }
                    _ => {
                        order_of_species.push(idx);
                        terms.insert(
                            idx,
                            kc::Term {
                                species: idx as u32,
                                reactant: 0.0,
                                product: 0.0,
                                order: a,
                            },
                        );
                    }
                }
            }
        }
        // the elements must balance
        for (e, at) in self.elements.iter().zip(&self.atoms) {
            let d: f64 = terms
                .values()
                .map(|t| (t.product - t.reactant) * at[t.species as usize])
                .sum();
            if d.abs() > 1e-9 {
                return Err(format!(
                    "`{eq}` does not conserve the element {e} (Δ = {d})"
                ));
            }
        }
        let order_sum: f64 = terms.values().map(|t| t.order).sum();
        // a reversible reaction with orders reverses with a' − ν' + ν''
        // (kernel module documentation)
        let reversible = parsed.reversible;
        let conv = |key: &str, extra: f64| -> std::result::Result<[f64; 3], String> {
            let y = &r[key];
            if y.is_badvalue() {
                return Err(format!("no `{key}`"));
            }
            units.arrhenius(y, order_sum + extra)
        };
        let (rate, low) = match kind {
            kc::ELEMENTARY => (conv("rate-constant", 0.0)?, [0.0; 3]),
            kc::THREE_BODY => (conv("rate-constant", 1.0)?, [0.0; 3]),
            _ => (
                conv("high-P-rate-constant", 0.0)?,
                conv("low-P-rate-constant", 1.0)?,
            ),
        };
        if rate[0] < 0.0 && !r["negative-A"].as_bool().unwrap_or(false) {
            return Err("a negative pre-exponential factor without `negative-A: true`".into());
        }
        let (mut troe, mut has_t2) = ([0.0; 4], 0u32);
        if kind == kc::TROE {
            let t = &r["Troe"];
            let g = |k: &str| number(&t[k]);
            troe = [
                g("A").ok_or("Troe without A")?,
                g("T3").ok_or("Troe without T3")?,
                g("T1").ok_or("Troe without T1")?,
                g("T2").unwrap_or(0.0),
            ];
            has_t2 = u32::from(g("T2").is_some());
        }
        let e0 = self.efficiencies.len() as u32;
        let default_efficiency = if kind == kc::ELEMENTARY {
            0.0
        } else {
            number(&r["default-efficiency"]).unwrap_or(1.0)
        };
        if kind != kc::ELEMENTARY
            && let Some(h) = r["efficiencies"].as_hash()
        {
            for (k, v) in h {
                let s = k.as_str().ok_or("an efficiency without a species name")?;
                let a =
                    number(v).ok_or_else(|| format!("the efficiency of `{s}` is not a number"))?;
                // species outside the phase do not take part
                if let Some(idx) = self.index(s) {
                    self.efficiencies.push(kc::Efficiency {
                        species: idx as u32,
                        value: a,
                    });
                }
            }
        }
        let t0 = self.terms.len() as u32;
        for k in order_of_species {
            self.terms.push(terms[&k]);
        }
        self.reactions.push(kc::Reaction {
            kind,
            reversible: u32::from(reversible),
            terms: [t0, self.terms.len() as u32],
            efficiencies: [e0, self.efficiencies.len() as u32],
            default_efficiency,
            rate,
            low,
            troe,
            has_t2,
        });
        self.equations.push(eq);
        Ok(())
    }
}

/// A number of a YAML value (an integer or a real).
pub(crate) fn number(y: &Yaml) -> Option<f64> {
    match y {
        Yaml::Real(s) => s.parse().ok(),
        Yaml::Integer(i) => Some(*i as f64),
        _ => None,
    }
}

/// The units of a file.
#[derive(Debug, Clone, Copy)]
struct Units {
    /// m³/mol of one length³/quantity.
    volume: f64,
    /// s of one time unit.
    time: f64,
    /// J/mol of one activation-energy unit; NaN for kelvin.
    energy: f64,
}

impl Units {
    fn of(y: &Yaml, name: &str) -> Result<Self> {
        // the defaults of the format: m, kmol, s, J/kmol
        let get = |k: &str, d: &str| y[k].as_str().unwrap_or(d).to_string();
        let length = match get("length", "m").as_str() {
            "m" => 1.0,
            "cm" => 1e-2,
            "mm" => 1e-3,
            u => return Err(bad("mechanism", format!("{name}: length unit `{u}`"))),
        };
        let quantity = match get("quantity", "kmol").as_str() {
            "mol" => 1.0,
            "kmol" => 1e3,
            u => return Err(bad("mechanism", format!("{name}: quantity unit `{u}`"))),
        };
        let time = match get("time", "s").as_str() {
            "s" => 1.0,
            "ms" => 1e-3,
            u => return Err(bad("mechanism", format!("{name}: time unit `{u}`"))),
        };
        let energy = energy_unit(&get("activation-energy", "J/kmol"))
            .ok_or_else(|| bad("mechanism", format!("{name}: activation-energy unit")))?;
        Ok(Self {
            volume: length * length * length / quantity,
            time,
            energy,
        })
    }

    /// [A, b, T_a] in SI of a rate constant whose concentration orders sum
    /// to `orders`.
    fn arrhenius(&self, y: &Yaml, orders: f64) -> std::result::Result<[f64; 3], String> {
        let a = number(&y["A"]).ok_or("A is not a number")?;
        let b = number(&y["b"]).unwrap_or(0.0);
        let ta = match &y["Ea"] {
            Yaml::String(s) => {
                // "<value> <unit>"
                let mut it = s.split_whitespace();
                let v: f64 = it
                    .next()
                    .and_then(|v| v.parse().ok())
                    .ok_or("Ea is not a number")?;
                let u = energy_unit(it.next().unwrap_or("")).ok_or("Ea unit")?;
                to_kelvin(v, u)
            }
            e => to_kelvin(number(e).unwrap_or(0.0), self.energy),
        };
        let a_si = a * libm_pow(self.volume, orders - 1.0) / self.time;
        Ok([a_si, b, ta])
    }
}

fn libm_pow(x: f64, y: f64) -> f64 {
    foehncraft_core::math::pow(x, y)
}

fn energy_unit(u: &str) -> Option<f64> {
    Some(match u {
        "J/mol" => 1.0,
        "kJ/mol" => 1e3,
        "J/kmol" => 1e-3,
        "cal/mol" => 4.184,
        "kcal/mol" => 4184.0,
        "K" => f64::NAN,
        _ => return None,
    })
}

fn to_kelvin(v: f64, unit: f64) -> f64 {
    if unit.is_nan() { v } else { v * unit / ks::R }
}

/// The third body of an equation.
#[derive(Debug, Clone, PartialEq)]
enum Third {
    None,
    /// `+ M`.
    M,
    /// `(+M)`.
    Falloff,
    /// `(+N2)` or `+ N2` as a named collider.
    Collider(String),
}

#[derive(Debug, Clone, PartialEq)]
struct Equation {
    reactants: Vec<(String, f64)>,
    products: Vec<(String, f64)>,
    reversible: bool,
    third_body: Third,
}

/// Parses `2 O + M <=> O2 + M`, `O + CO (+M) <=> CO2 (+M)`,
/// `CH4 + 1.5 O2 => CO + 2 H2O`.
fn parse_equation(eq: &str) -> std::result::Result<Equation, String> {
    let (lhs, rhs, reversible) = if let Some((l, r)) = eq.split_once("<=>") {
        (l, r, true)
    } else if let Some((l, r)) = eq.split_once("=>") {
        (l, r, false)
    } else if let Some((l, r)) = eq.split_once('=') {
        (l, r, true)
    } else {
        return Err(format!("`{eq}` has no `<=>`, `=>` or `=`"));
    };
    let mut third = Third::None;
    let mut side = |s: &str| -> std::result::Result<Vec<(String, f64)>, String> {
        let mut s = s.to_string();
        // (+M) or (+X)
        if let Some(a) = s.find("(+") {
            let b = s[a..].find(')').ok_or("an unclosed `(+`")? + a;
            let c = s[a + 2..b].trim().to_string();
            third = if c == "M" {
                Third::Falloff
            } else {
                Third::Collider(c)
            };
            s.replace_range(a..=b, "");
        }
        let mut out = Vec::new();
        for part in s.split(" + ") {
            let p = part.trim();
            if p.is_empty() {
                return Err(format!("an empty term in `{eq}`"));
            }
            if p == "M" {
                if third == Third::None {
                    third = Third::M;
                }
                continue;
            }
            let (nu, name) = match p.split_once(' ') {
                Some((a, b)) => match a.parse::<f64>() {
                    Ok(v) => (v, b.trim()),
                    Err(_) => return Err(format!("the term `{p}` in `{eq}`")),
                },
                None => {
                    // a coefficient written without a space (`2O`)
                    let digits: String = p
                        .chars()
                        .take_while(|c| c.is_ascii_digit() || *c == '.')
                        .collect();
                    if digits.is_empty() || digits.len() == p.len() {
                        (1.0, p)
                    } else {
                        (
                            digits.parse().map_err(|_| format!("the term `{p}`"))?,
                            &p[digits.len()..],
                        )
                    }
                }
            };
            if !(nu > 0.0) {
                return Err(format!("the coefficient of `{name}` in `{eq}`"));
            }
            out.push((name.to_string(), nu));
        }
        Ok(out)
    };
    let reactants = side(lhs)?;
    let products = side(rhs)?;
    Ok(Equation {
        reactants,
        products,
        reversible,
        third_body: third,
    })
}

/// A species record of a YAML file.
fn record(y: &Yaml, s: &str, file: &str) -> Result<Record> {
    let fail = |m: &str| bad("mechanism", format!("{file}: species `{s}`: {m}"));
    let mut composition = BTreeMap::new();
    for (k, v) in y["composition"]
        .as_hash()
        .ok_or_else(|| fail("no composition"))?
    {
        let e = k.as_str().ok_or_else(|| fail("an element name"))?;
        composition.insert(
            e.to_string(),
            number(v).ok_or_else(|| fail("an atom count"))?,
        );
    }
    let mut molar_mass = 0.0;
    for (e, n) in &composition {
        let w = species::atomic_weight(e);
        if w.is_nan() {
            return Err(fail(&format!(
                "the element {e} is not known (H, C, O, N, Ar)"
            )));
        }
        molar_mass += n * w;
    }
    let th = &y["thermo"];
    if th["model"].as_str() != Some("NASA7") {
        return Err(fail("only NASA7 thermo is read"));
    }
    let temperature: Vec<f64> = th["temperature-ranges"]
        .as_vec()
        .ok_or_else(|| fail("no temperature-ranges"))?
        .iter()
        .filter_map(number)
        .collect();
    let coefficients: Vec<Vec<f64>> = th["data"]
        .as_vec()
        .ok_or_else(|| fail("no thermo data"))?
        .iter()
        .map(|r| {
            r.as_vec()
                .map(|v| v.iter().filter_map(number).collect())
                .unwrap_or_default()
        })
        .collect();
    let mut temperature = temperature;
    let mut coefficients = coefficients;
    if temperature.len() == 2 && coefficients.len() == 1 {
        // one range: the same coefficients above
        temperature = vec![temperature[0], temperature[1], temperature[1]];
        coefficients.push(coefficients[0].clone());
    }
    if temperature.len() != 3
        || coefficients.len() != 2
        || coefficients.iter().any(|c| c.len() != 7)
    {
        return Err(fail("NASA7 data need two ranges of 7 coefficients"));
    }
    let tr = &y["transport"];
    let g = |k: &str| number(&tr[k]);
    Ok(Record {
        composition,
        molar_mass,
        temperature,
        coefficients,
        geometry: tr["geometry"].as_str().unwrap_or("nonlinear").to_string(),
        well_depth: g("well-depth").unwrap_or(100.0),
        diameter: g("diameter").unwrap_or(3.5),
        dipole: g("dipole").unwrap_or(0.0),
        polarizability: g("polarizability").unwrap_or(0.0),
        rotational_relaxation: g("rotational-relaxation").unwrap_or(0.0),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equations_parse() {
        let e = parse_equation("2 O + M <=> O2 + M").unwrap();
        assert_eq!(e.reactants, [("O".to_string(), 2.0)]);
        assert_eq!(e.third_body, Third::M);
        assert!(e.reversible);
        let e = parse_equation("O + CO (+M) <=> CO2 (+M)").unwrap();
        assert_eq!(e.third_body, Third::Falloff);
        assert_eq!(e.products, [("CO2".to_string(), 1.0)]);
        let e = parse_equation("CH4 + 1.5 O2 => CO + 2 H2O").unwrap();
        assert!(!e.reversible);
        assert_eq!(e.reactants[1], ("O2".to_string(), 1.5));
        let e = parse_equation("H + O2 + O2 <=> HO2 + O2").unwrap();
        assert_eq!(e.reactants.len(), 3);
        assert_eq!(e.third_body, Third::None);
        assert!(parse_equation("A + B").is_err());
        assert!(matches!(
            parse_equation("H + O2 (+N2) <=> HO2 (+N2)")
                .unwrap()
                .third_body,
            Third::Collider(_)
        ));
    }

    #[test]
    fn the_built_in_mechanisms_load() {
        let g = Mechanism::load("gri30", None).unwrap();
        assert_eq!(g.len(), 53);
        assert_eq!(g.reactions.len(), 325);
        assert_eq!(
            g.reactions.iter().filter(|r| r.kind == kc::TROE).count(),
            26
        );
        assert_eq!(
            g.reactions
                .iter()
                .filter(|r| r.kind == kc::LINDEMANN)
                .count(),
            3
        );
        // reaction 3, O + H2 <=> H + OH: A = 3.87e4 cm³/(mol s) = 3.87e-2 m³/(mol s)
        assert!((g.reactions[2].rate[0] / 3.87e-2 - 1.0).abs() < 1e-12);
        assert!((g.reactions[2].rate[2] - 6260.0 * 4.184 / ks::R).abs() < 1e-9);
        // reaction 1, 2 O + M: A = 1.2e17 cm⁶/(mol² s) = 1.2e5 m⁶/(mol² s)
        assert!((g.reactions[0].rate[0] / 1.2e5 - 1.0).abs() < 1e-12);
        // Jones & Lindstedt in m, kmol, s give the cm, mol, s values of Si
        // et al. (2020, table 1): 7.82e13 and 1.21e18 for the orders 1.75
        let jl = Mechanism::load("jones-lindstedt-4", None).unwrap();
        let cgs = |a: f64| a * foehncraft_core::math::pow(1e-6, 0.75);
        assert!((jl.reactions[0].rate[0] / cgs(7.82e13) - 1.0).abs() < 2e-3);
        assert!((jl.reactions[3].rate[0] / cgs(1.21e18) - 1.0).abs() < 2e-3);
        assert!((jl.reactions[1].rate[0] / (3.0e11 * 1e-6) - 1.0).abs() < 2e-3);
        for n in ["westbrook-dryer-2", "jones-lindstedt-4"] {
            let m = Mechanism::load(n, None).unwrap();
            assert!(m.reactions.len() >= 3, "{n}");
        }
        assert!(Mechanism::load("nonsense", None).is_err());
    }

    #[test]
    fn unbalanced_and_unsupported_reactions_are_refused() {
        let base = "units: {length: cm, quantity: mol, activation-energy: cal/mol}\n\
                    phases:\n- name: g\n  species: [H2, O2, H2O, OH, N2]\n\
                    reactions:\n";
        let bad_eq =
            format!("{base}- equation: H2 + O2 => H2O\n  rate-constant: {{A: 1.0, b: 0, Ea: 0}}\n");
        let e = Mechanism::from_yaml(&bad_eq, "t").unwrap_err().to_string();
        assert!(e.contains("element O"), "{e}");
        let plog = format!(
            "{base}- equation: H2 + O2 <=> 2 OH\n  type: pressure-dependent-Arrhenius\n  rate-constants: []\n"
        );
        let e = Mechanism::from_yaml(&plog, "t").unwrap_err().to_string();
        assert!(e.contains("not implemented"), "{e}");
        assert_eq!(
            diag::of_error(&Mechanism::from_yaml(&plog, "t").unwrap_err()).code,
            diag::CHEMISTRY
        );
    }
}
