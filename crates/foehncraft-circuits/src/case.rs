//! The circuit case file (TOML): models, nodes and tubes.

use crate::{Error, codes};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// The two-phase friction model.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Friction {
    /// The multiplier of Friedel (1979) on the whole flow as liquid.
    #[default]
    Friedel,
    /// The homogeneous model: Churchill's friction factor with the
    /// McAdams mixture viscosity and the homogeneous density.
    Homogeneous,
}

/// The void fraction model.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Void {
    /// Zuber–Findlay drift flux with Rouhani–Axelsson (1970).
    #[default]
    DriftFlux,
    /// No slip between the phases.
    Homogeneous,
}

/// The flow boiling heat transfer correlation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Boiling {
    /// Chen (1966) with the subcooled form of Butterworth (1977).
    #[default]
    Chen,
    /// Gungor & Winterton (1986).
    GungorWinterton,
}

/// The vapour of subcooled boiling.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Subcooled {
    /// Onset of significant void by Saha & Zuber (1974), true quality by
    /// the profile of Levy (1967).
    #[default]
    SahaZuber,
    /// Thermal equilibrium: no vapour below x_eq = 0.
    None,
}

/// The critical heat flux source.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Chf {
    /// Bowring (1972), built in.
    #[default]
    Bowring,
    /// A look-up table in the layout of the 2006 table (`chf_table`).
    Table,
    /// No critical heat flux check.
    None,
}

/// The `[circuit]` table: name and models.
#[derive(Clone, Debug, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Header {
    /// Name of the case (default: the file stem).
    #[serde(default)]
    pub name: Option<String>,
    /// Two-phase friction model.
    #[serde(default)]
    pub friction: Friction,
    /// Void fraction model.
    #[serde(default)]
    pub void: Void,
    /// Flow boiling heat transfer.
    #[serde(default)]
    pub boiling: Boiling,
    /// Vapour in subcooled boiling.
    #[serde(default)]
    pub subcooled: Subcooled,
    /// Critical heat flux source.
    #[serde(default)]
    pub chf: Chf,
    /// The CSV file of the CHF look-up table (`chf = "table"`), relative
    /// to the case file: header `p_kpa,g,x,chf_kw`.
    #[serde(default)]
    pub chf_table: Option<String>,
    /// The smallest CHF ratio q_CHF/q accepted without CIR004 (default 1.3).
    #[serde(default = "default_chf_margin")]
    pub chf_margin: f64,
    /// Segments of every tube unless the tube sets its own (default 40).
    #[serde(default = "default_segments")]
    pub segments: usize,
    /// Newton tolerance on the scaled residuals (default 1e-10).
    #[serde(default = "default_tolerance")]
    pub tolerance: f64,
    /// Newton iterations at most (default 100).
    #[serde(default = "default_iterations")]
    pub max_iterations: usize,
}

fn default_chf_margin() -> f64 {
    1.3
}
fn default_segments() -> usize {
    40
}
fn default_tolerance() -> f64 {
    1e-10
}
fn default_iterations() -> usize {
    100
}

impl Default for Header {
    fn default() -> Self {
        Self {
            name: None,
            friction: Friction::default(),
            void: Void::default(),
            boiling: Boiling::default(),
            subcooled: Subcooled::default(),
            chf: Chf::default(),
            chf_table: None,
            chf_margin: default_chf_margin(),
            segments: default_segments(),
            tolerance: default_tolerance(),
            max_iterations: default_iterations(),
        }
    }
}

/// The kind of a node.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum NodeKind {
    /// A header or junction: pressure and enthalpy follow from the
    /// balances of mass and energy (an optional `mass_flow` enters here).
    Junction,
    /// A fixed pressure; water leaving it into a tube has the given
    /// temperature or enthalpy.
    Pressure,
    /// A steam drum at a fixed pressure: the steam of the risers leaves,
    /// the same mass of feedwater (`feed_temperature`) mixes with the
    /// separated water that goes to the downcomers.
    Drum,
}

/// A node of the circuit.
#[derive(Clone, Debug, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Node {
    /// Name, used by the tubes.
    pub name: String,
    /// Kind of node.
    pub kind: NodeKind,
    /// Pressure (Pa) of a `pressure` or `drum` node; the start value of a
    /// junction.
    #[serde(default)]
    pub pressure: Option<f64>,
    /// Temperature (K) of the water that enters the circuit here (a
    /// `pressure` node or the `mass_flow` of a junction). It is turned
    /// into an enthalpy at the node's `pressure`, or for a junction
    /// without one at the mean of the fixed pressures of the case; give
    /// `enthalpy` to fix the energy exactly.
    #[serde(default)]
    pub temperature: Option<f64>,
    /// Specific enthalpy (J/kg) instead of the temperature.
    #[serde(default)]
    pub enthalpy: Option<f64>,
    /// Mass flow (kg/s) that enters the circuit at a junction (negative:
    /// leaves it).
    #[serde(default)]
    pub mass_flow: Option<f64>,
    /// Feedwater temperature (K) of a drum.
    #[serde(default)]
    pub feed_temperature: Option<f64>,
    /// Height of the node (m), for the report only; the tubes carry the
    /// rise.
    #[serde(default)]
    pub elevation: Option<f64>,
}

/// A tube (or a bank of identical parallel tubes) between two nodes.
#[derive(Clone, Debug, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Tube {
    /// Name.
    pub name: String,
    /// The node where positive flow enters.
    pub from: String,
    /// The node where positive flow leaves.
    pub to: String,
    /// Length (m).
    pub length: f64,
    /// Inner diameter (m).
    pub diameter: f64,
    /// Rise from `from` to `to` (m), between −length and length (default 0).
    #[serde(default)]
    pub rise: f64,
    /// Absolute roughness (m) (default 0, smooth).
    #[serde(default)]
    pub roughness: f64,
    /// Identical tubes in parallel (default 1).
    #[serde(default = "one")]
    pub count: usize,
    /// Segments of the march (default: `[circuit] segments`).
    #[serde(default)]
    pub segments: Option<usize>,
    /// Uniform heat flux on the inner wall (W/m²).
    #[serde(default)]
    pub heat_flux: Option<f64>,
    /// Heat of one tube (W) instead of the heat flux.
    #[serde(default)]
    pub heat: Option<f64>,
    /// Axial shape of the heat flux: pairs [s, factor] with s from 0 at
    /// `from` to 1 at `to`, linear between; scaled so that the mean is one.
    #[serde(default)]
    pub heat_shape: Option<Vec<[f64; 2]>>,
    /// Unheated length at the `from` end (m) (default 0).
    #[serde(default)]
    pub unheated_inlet: f64,
    /// Extra loss coefficient K (per tube), Δp = K G²/(2ρ) at the inlet.
    #[serde(default)]
    pub loss: f64,
    /// Start value of the mass flow of one tube (kg/s).
    #[serde(default)]
    pub start_flow: Option<f64>,
}

fn one() -> usize {
    1
}

/// The `[output]` table.
#[derive(Clone, Debug, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Output {
    /// Write `profiles.csv` with every segment (default true).
    #[serde(default)]
    pub profiles: Option<bool>,
}

/// A circuit case.
#[derive(Clone, Debug, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Case {
    /// Name and models.
    #[serde(default)]
    pub circuit: Header,
    /// The nodes.
    pub nodes: Vec<Node>,
    /// The tubes.
    pub tubes: Vec<Tube>,
    /// Outputs.
    #[serde(default)]
    pub output: Output,
    /// The folder of the case file, for relative paths (not a key).
    #[serde(skip)]
    #[schemars(skip)]
    pub base: Option<std::path::PathBuf>,
}

impl Case {
    /// Read a case from TOML text.
    pub fn from_toml(text: &str) -> Result<Self, Error> {
        let c: Case = toml::from_str(text).map_err(|e| Error::new(codes::CASE, e.to_string()))?;
        c.check()?;
        Ok(c)
    }

    /// Read a case file; relative paths are taken from its folder.
    pub fn from_file(path: &std::path::Path) -> Result<Self, Error> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::new(codes::CASE, format!("{}: {e}", path.display())))?;
        let mut c = Self::from_toml(&text)?;
        c.base = Some(
            path.parent()
                .map_or_else(Default::default, |p| p.to_path_buf()),
        );
        Ok(c)
    }

    /// The path of the CHF table, if any.
    pub fn chf_table_path(&self) -> Option<std::path::PathBuf> {
        self.circuit.chf_table.as_ref().map(|t| match &self.base {
            Some(b) => b.join(t),
            None => t.into(),
        })
    }

    fn check(&self) -> Result<(), Error> {
        let bad = |m: String| Err(Error::new(codes::CASE, m));
        if self.nodes.is_empty() || self.tubes.is_empty() {
            return bad("a circuit needs nodes and tubes".into());
        }
        let mut names = std::collections::HashSet::new();
        let mut fixed = 0;
        for n in &self.nodes {
            if !names.insert(n.name.as_str()) {
                return bad(format!("node `{}` twice", n.name));
            }
            let has_state = n.temperature.is_some() || n.enthalpy.is_some();
            if n.temperature.is_some() && n.enthalpy.is_some() {
                return bad(format!(
                    "node `{}`: give temperature or enthalpy, not both",
                    n.name
                ));
            }
            match n.kind {
                NodeKind::Pressure | NodeKind::Drum => {
                    fixed += 1;
                    match n.pressure {
                        Some(p) if p > 0.0 => {}
                        _ => {
                            return bad(format!(
                                "node `{}`: a positive `pressure` is needed",
                                n.name
                            ));
                        }
                    }
                    if n.kind == NodeKind::Pressure && !has_state {
                        return bad(format!(
                            "node `{}`: a pressure node needs the `temperature` or `enthalpy` of the water it delivers",
                            n.name
                        ));
                    }
                    if n.kind == NodeKind::Drum {
                        if n.feed_temperature.is_none() {
                            return bad(format!("drum `{}`: `feed_temperature` is needed", n.name));
                        }
                        if n.pressure.is_some_and(|p| p >= crate::fluid::P_CRIT) {
                            return bad(format!(
                                "drum `{}`: a drum needs a pressure below the critical pressure",
                                n.name
                            ));
                        }
                    }
                    if n.mass_flow.is_some() {
                        return bad(format!(
                            "node `{}`: `mass_flow` belongs to a junction",
                            n.name
                        ));
                    }
                }
                NodeKind::Junction => {
                    if n.mass_flow.is_some_and(|m| m > 0.0) && !has_state {
                        return bad(format!(
                            "junction `{}`: the inflow `mass_flow` needs its `temperature` or `enthalpy`",
                            n.name
                        ));
                    }
                    if n.feed_temperature.is_some() {
                        return bad(format!(
                            "junction `{}`: `feed_temperature` belongs to a drum",
                            n.name
                        ));
                    }
                }
            }
        }
        if fixed == 0 {
            return bad(
                "the circuit needs a `pressure` or `drum` node (the pressure level)".into(),
            );
        }
        let mut tnames = std::collections::HashSet::new();
        for t in &self.tubes {
            if !tnames.insert(t.name.as_str()) {
                return bad(format!("tube `{}` twice", t.name));
            }
            for end in [&t.from, &t.to] {
                if !names.contains(end.as_str()) {
                    return bad(format!("tube `{}`: no node `{end}`", t.name));
                }
            }
            if t.from == t.to {
                return bad(format!(
                    "tube `{}`: `from` and `to` are the same node",
                    t.name
                ));
            }
            if !(t.length > 0.0 && t.diameter > 0.0) {
                return bad(format!(
                    "tube `{}`: length and diameter must be positive",
                    t.name
                ));
            }
            if t.rise.abs() > t.length * (1.0 + 1e-12) {
                return bad(format!(
                    "tube `{}`: |rise| is larger than the length",
                    t.name
                ));
            }
            if t.roughness < 0.0 || t.count == 0 || t.loss < 0.0 {
                return bad(format!(
                    "tube `{}`: roughness, count and loss must not be negative (count ≥ 1)",
                    t.name
                ));
            }
            if t.heat.is_some() && t.heat_flux.is_some() {
                return bad(format!(
                    "tube `{}`: give heat or heat_flux, not both",
                    t.name
                ));
            }
            if !(0.0..t.length).contains(&t.unheated_inlet) {
                return bad(format!(
                    "tube `{}`: unheated_inlet must be in [0, length)",
                    t.name
                ));
            }
            if let Some(s) = &t.heat_shape
                && (s.len() < 2
                    || s[0][0] != 0.0
                    || s[s.len() - 1][0] != 1.0
                    || s.windows(2).any(|w| w[1][0] <= w[0][0])
                    || s.iter().any(|v| v[1] < 0.0))
            {
                return bad(format!(
                    "tube `{}`: heat_shape needs increasing s from 0 to 1 and factors ≥ 0",
                    t.name
                ));
            }
            if t.segments == Some(0) {
                return bad(format!("tube `{}`: segments ≥ 1", t.name));
            }
        }
        if self.circuit.segments == 0 || self.circuit.chf_margin <= 0.0 {
            return bad("[circuit]: segments ≥ 1 and chf_margin > 0".into());
        }
        if self.circuit.chf == Chf::Table && self.circuit.chf_table.is_none() {
            return bad("[circuit]: chf = \"table\" needs `chf_table`".into());
        }
        Ok(())
    }
}
