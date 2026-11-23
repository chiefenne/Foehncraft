//! Pressure values absolute by default, gauge optional (D239).
//!
//! The case file speaks the language of the application: from case format
//! 0.7 its pressure values are absolute static pressures
//! (`[fluid] pressure_values = "absolute"`, the default). The solver keeps
//! working with the gauge pressure p − p_op, relative to an operating
//! pressure p_op, so that round-off and convergence do not change.
//!
//! The reader turns the case into its *gauge twin* once, in
//! [`resolve`]:
//!
//! - every fluid part (fluid cells connected without a solid between them,
//!   D194) takes its operating pressure from its own pressure boundaries:
//!   the smallest number of its static-pressure outlets, else of its
//!   total-pressure outlets, else of its pressure inlets; a part without a
//!   pressure boundary from its `[pressure] reference_value` or
//!   `[initial] pressure`, else 101 325 Pa. An explicit `[fluid]
//!   operating_pressure` or `[regions.<name>] operating_pressure` wins;
//! - every pressure value of the file (boundary `pressure` and
//!   `total_pressure`, numbers and expressions, `[initial] pressure`,
//!   `[pressure] reference_value`) becomes p − p_op of its part; for a
//!   number within a factor 2 of p_op the difference is exact (Sterbenz's
//!   lemma), so an absolute case and its gauge twin run the same bits;
//! - `[fluid] operating_pressure` and the operating pressure of the fluid
//!   regions are set, so that each gas's density comes from the absolute
//!   pressure of its own part, and `pressure_values` becomes `"gauge"`:
//!   `case.resolved.toml` is the gauge twin.
//!
//! The [`Frame`] keeps what the file said: the outputs (the `pressure`
//! field, probes, boundary values) go back to absolute values by adding
//! the operating pressure of each cell's part, and the log and `flow.json`
//! state the convention and the operating pressure of every part (PRS005).
//! Files of format 0.6 and older mean gauge: the case-format migration puts
//! `pressure_values = "gauge"` into them (a changed default), so they keep
//! their meaning and their results bit for bit.

use crate::bc::BoundaryCondition;
use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use crate::expr::Scalar;
use crate::materials::{DEFAULT_PRESSURE, RegionKind};
use anyhow::Result;
use foehncraft_mesh::Mesh;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

/// What the pressure values of a case file are.
#[derive(
    Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "lowercase")]
pub enum PressureValues {
    /// Absolute static pressures (the default from case format 0.7).
    #[default]
    Absolute,
    /// Values relative to the operating pressure, p = p_op + p_gauge (the
    /// meaning of files of format 0.6 and older).
    Gauge,
}

impl PressureValues {
    /// `"absolute"` or `"gauge"`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Absolute => "absolute",
            Self::Gauge => "gauge",
        }
    }
}

/// Which pressure values the case file gives (a value left out sits at the
/// operating pressure of its part, gauge 0).
#[derive(Debug, Clone, Copy, Default)]
pub struct Given {
    /// `[initial] pressure`.
    pub initial: bool,
    /// `[pressure] reference_value`.
    pub reference: bool,
}

impl Given {
    /// From the table of the case file.
    pub fn of(t: &toml::Table) -> Self {
        let has = |a: &str, b: &str| {
            t.get(a)
                .and_then(toml::Value::as_table)
                .is_some_and(|s| s.contains_key(b))
        };
        Self {
            initial: has("initial", "pressure"),
            reference: has("pressure", "reference_value"),
        }
    }
}

/// One fluid part: its operating pressure and where it comes from.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct PartFrame {
    /// The operating pressure p_op (Pa).
    pub operating_pressure: f64,
    /// Where p_op comes from (`boundary <name>`, `[fluid]
    /// operating_pressure`, …).
    pub source: String,
    /// The boundaries of the part.
    pub boundaries: Vec<String>,
    /// The fluid regions of the part (`fluid`: cells outside the regions).
    pub regions: Vec<String>,
    /// The smallest absolute pressure value given on the part (Pa), for
    /// the check of PRS002.
    #[serde(skip)]
    pub smallest: Option<f64>,
}

/// The pressure convention of a case file and the operating pressure of
/// every fluid part.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct Frame {
    /// What the values of the file are.
    pub pressure_values: PressureValues,
    /// The fluid parts (one for a case of one connected fluid).
    pub parts: Vec<PartFrame>,
}

impl Frame {
    /// True when the file's values are absolute.
    pub fn absolute(&self) -> bool {
        self.pressure_values == PressureValues::Absolute
    }

    /// The operating pressure of the part of `boundary` (the first part
    /// when no part has it).
    pub fn of_boundary(&self, boundary: &str) -> f64 {
        self.parts
            .iter()
            .find(|p| p.boundaries.iter().any(|b| b == boundary))
            .or(self.parts.first())
            .map_or(DEFAULT_PRESSURE, |p| p.operating_pressure)
    }

    /// The offset that turns the solver's gauge pressure of each cell of
    /// `mesh` into the value of the case file's convention: `None` with
    /// gauge values (the outputs stay as they are, bit for bit). The parts
    /// of `mesh` (connected cells, without the cells of solid regions) are
    /// matched to the parts of the frame by their boundaries, then by their
    /// regions.
    pub fn cell_offsets(&self, case: &FlowCase, mesh: &Mesh) -> Option<Vec<f64>> {
        if !self.absolute() {
            return None;
        }
        if self.parts.len() == 1 {
            return Some(vec![self.parts[0].operating_pressure; mesh.n_cells]);
        }
        let comp = components(case, mesh);
        let mut op = vec![self.parts[0].operating_pressure; comp.n];
        let mut found = vec![false; comp.n];
        for b in &mesh.boundaries {
            for f in b.range() {
                let k = comp.of[mesh.owner[f] as usize];
                if k == u32::MAX || found[k as usize] {
                    continue;
                }
                if let Some(p) = self
                    .parts
                    .iter()
                    .find(|p| p.boundaries.iter().any(|n| *n == b.name))
                {
                    op[k as usize] = p.operating_pressure;
                    found[k as usize] = true;
                }
            }
        }
        let region_of = zone_regions(case, mesh);
        for c in 0..mesh.n_cells {
            let k = comp.of[c];
            if k == u32::MAX || found[k as usize] {
                continue;
            }
            let r = &region_of[mesh.cell_zone[c] as usize];
            if let Some(p) = self.parts.iter().find(|p| p.regions.iter().any(|n| n == r)) {
                op[k as usize] = p.operating_pressure;
                found[k as usize] = true;
            }
        }
        Some(
            (0..mesh.n_cells)
                .map(|c| match comp.of[c] {
                    u32::MAX => 0.0,
                    k => op[k as usize],
                })
                .collect(),
        )
    }

    /// The offset of the boundary values of `boundary`: its part's
    /// operating pressure with absolute values, else 0.
    pub fn boundary_offset(&self, boundary: &str) -> f64 {
        if self.absolute() {
            self.of_boundary(boundary)
        } else {
            0.0
        }
    }

    /// The frame as JSON (`flow.json`).
    pub fn json(&self) -> serde_json::Value {
        serde_json::to_value(self).unwrap_or_default()
    }

    /// The PRS005 lines of the log.
    pub fn notes(&self) -> Vec<String> {
        let mut v = vec![format!(
            "[PRS005] pressure values are {} (D239){}",
            self.pressure_values.as_str(),
            if self.absolute() {
                "; the solver works with p - p_op, the outputs are absolute"
            } else {
                "; p = p_op + p_gauge"
            }
        )];
        for (k, p) in self.parts.iter().enumerate() {
            v.push(format!(
                "[PRS005]   part {}{}: operating pressure {} Pa ({})",
                k + 1,
                if p.regions.is_empty() || self.parts.len() == 1 {
                    String::new()
                } else {
                    format!(" (regions {})", p.regions.join(", "))
                },
                p.operating_pressure,
                p.source
            ));
        }
        v
    }

    /// The gauge value of a pressure key set at run time (`[[schedule]]`,
    /// D128) in the file's convention: `None` when `key` is not a pressure
    /// value or the values are gauge.
    pub fn gauge_value(&self, key: &str, v: &toml::Value) -> Option<toml::Value> {
        if !self.absolute() {
            return None;
        }
        let op = if let Some(rest) = key.strip_prefix("boundary.") {
            let (b, k) = rest.rsplit_once('.')?;
            if k != "pressure" && k != "total_pressure" {
                return None;
            }
            self.of_boundary(b)
        } else if key == "initial.pressure" || key == "pressure.reference_value" {
            self.parts.first()?.operating_pressure
        } else {
            return None;
        };
        match v {
            toml::Value::Float(x) => Some(toml::Value::Float(x - op)),
            toml::Value::Integer(i) => Some(toml::Value::Float(*i as f64 - op)),
            toml::Value::String(s) => Some(toml::Value::String(minus(s, op))),
            _ => None,
        }
    }
}

/// An expression minus the operating pressure.
fn minus(s: &str, op: f64) -> String {
    format!("({s}) - {op:?}")
}

fn scalar_minus(s: &Scalar, op: f64) -> Scalar {
    match s {
        Scalar::Number(v) => Scalar::Number(v - op),
        Scalar::Expr(e) => Scalar::Expr(minus(e, op)),
    }
}

/// The pressure value of a pressure boundary and how it ranks as the
/// source of the operating pressure (0: static-pressure outlet, …).
fn pressure_of(bc: &BoundaryCondition) -> Option<(&Scalar, u8)> {
    match bc {
        BoundaryCondition::StaticPressureOutlet { pressure, .. } => Some((pressure, 0)),
        BoundaryCondition::TotalPressureOutlet { total_pressure, .. } => Some((total_pressure, 1)),
        BoundaryCondition::StaticPressureInlet { pressure, .. } => Some((pressure, 2)),
        BoundaryCondition::TotalPressureInlet { total_pressure, .. } => Some((total_pressure, 3)),
        _ => None,
    }
}

fn pressure_of_mut(bc: &mut BoundaryCondition) -> Option<&mut Scalar> {
    match bc {
        BoundaryCondition::StaticPressureOutlet { pressure, .. }
        | BoundaryCondition::StaticPressureInlet { pressure, .. } => Some(pressure),
        BoundaryCondition::TotalPressureOutlet { total_pressure, .. }
        | BoundaryCondition::TotalPressureInlet { total_pressure, .. } => Some(total_pressure),
        _ => None,
    }
}

/// The connected parts of the fluid cells of `mesh` (cells of solid regions
/// are in none: `u32::MAX`).
struct Components {
    of: Vec<u32>,
    n: usize,
}

fn zone_regions(case: &FlowCase, mesh: &Mesh) -> Vec<String> {
    mesh.zones
        .iter()
        .map(|z| {
            case.regions
                .iter()
                .find(|(_, r)| r.zones.contains(z))
                .map_or_else(|| "fluid".to_string(), |(n, _)| n.clone())
        })
        .collect()
}

fn components(case: &FlowCase, mesh: &Mesh) -> Components {
    let solid: Vec<bool> = mesh
        .zones
        .iter()
        .map(|z| {
            case.regions
                .values()
                .any(|r| r.kind != RegionKind::Fluid && r.zones.contains(z))
        })
        .collect();
    let nc = mesh.n_cells;
    let is_fluid = |c: usize| {
        !solid
            .get(mesh.cell_zone[c] as usize)
            .copied()
            .unwrap_or(false)
    };
    let mut up: Vec<u32> = (0..nc as u32).collect();
    fn root(up: &mut [u32], mut c: u32) -> u32 {
        while up[c as usize] != c {
            let g = up[up[c as usize] as usize];
            up[c as usize] = g;
            c = g;
        }
        c
    }
    for f in 0..mesh.n_interior_faces() {
        let (o, n) = (mesh.owner[f] as usize, mesh.neighbour[f] as usize);
        if !(is_fluid(o) && is_fluid(n)) {
            continue;
        }
        let (a, b) = (root(&mut up, o as u32), root(&mut up, n as u32));
        if a != b {
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            up[hi as usize] = lo;
        }
    }
    let mut id = vec![u32::MAX; nc];
    let mut of = vec![u32::MAX; nc];
    let mut n = 0u32;
    for c in 0..nc {
        if !is_fluid(c) {
            continue;
        }
        let r = root(&mut up, c as u32) as usize;
        if id[r] == u32::MAX {
            id[r] = n;
            n += 1;
        }
        of[c] = id[r];
    }
    Components { of, n: n as usize }
}

/// A group of boundaries and regions that belong to one fluid part.
#[derive(Debug, Clone, Default)]
struct Group {
    boundaries: Vec<String>,
    regions: Vec<String>,
}

/// The fluid parts of the case: one group for a case without regions or
/// with one candidate operating pressure; else from the mesh.
fn groups(case: &FlowCase) -> Result<Vec<Group>> {
    let all = || Group {
        boundaries: case.boundary.keys().cloned().collect(),
        regions: Vec::new(),
    };
    let numbers: Vec<u64> = {
        let mut v: Vec<u64> = case
            .boundary
            .values()
            .filter_map(pressure_of)
            .filter_map(|(s, _)| s.number().map(f64::to_bits))
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let expressions = case
        .boundary
        .values()
        .filter_map(pressure_of)
        .any(|(s, _)| s.number().is_none());
    if case.regions.is_empty()
        || case.fluid.operating_pressure.is_some()
            && case
                .regions
                .values()
                .all(|r| r.operating_pressure.is_none())
        || (numbers.len() <= 1
            && !expressions
            && case
                .regions
                .values()
                .all(|r| r.operating_pressure.is_none()))
    {
        return Ok(vec![all()]);
    }
    // the parts need the mesh: the first level of a sweep (the parts are
    // the same on every level)
    let (mesh, _, _) = crate::case::mesh_for(case, case.verification.levels.first().copied())?;
    let comp = components(case, &mesh);
    let region_of = zone_regions(case, &mesh);
    let mut g = vec![Group::default(); comp.n];
    for c in 0..mesh.n_cells {
        let k = comp.of[c];
        if k == u32::MAX {
            continue;
        }
        let r = &region_of[mesh.cell_zone[c] as usize];
        let gr = &mut g[k as usize];
        if !gr.regions.contains(r) {
            gr.regions.push(r.clone());
        }
    }
    for b in &mesh.boundaries {
        let mut seen = vec![false; comp.n];
        for f in b.range() {
            let k = comp.of[mesh.owner[f] as usize];
            if k != u32::MAX {
                seen[k as usize] = true;
            }
        }
        for (k, s) in seen.iter().enumerate() {
            if *s {
                g[k].boundaries.push(b.name.clone());
            }
        }
    }
    for x in &mut g {
        x.regions.sort();
    }
    if g.is_empty() {
        return Ok(vec![all()]);
    }
    Ok(g)
}

fn not_positive(what: &str, v: f64, key: String) -> anyhow::Error {
    diag::err(
        Diagnostic::new(
            diag::PRS_NOT_POSITIVE,
            format!(
                "{what} = {v} Pa is an absolute pressure (D239: `[fluid] pressure_values = \
                 \"absolute\"`, the default from case format 0.7) and must be > 0; for a value \
                 relative to the operating pressure set `pressure_values = \"gauge\"`"
            ),
        )
        .key(key),
    )
}

/// Turns the absolute pressure values of `case` into the gauge values of
/// the solver and records the [`Frame`] (see the module documentation).
/// A gauge case keeps its values; its frame states the operating pressure.
pub fn resolve(case: &mut FlowCase, given: &Given) -> Result<()> {
    let written = case.fluid.pressure_values.unwrap_or_default();
    for (name, r) in &case.regions {
        if let Some(p) = r.operating_pressure
            && !(p > 0.0 && p.is_finite())
        {
            return Err(not_positive(
                &format!("[regions.{name}] operating_pressure"),
                p,
                format!("regions.{name}.operating_pressure"),
            ));
        }
    }
    if written == PressureValues::Gauge {
        let op = case.fluid.conditions().pressure;
        case.pressure_frame = Some(Arc::new(Frame {
            pressure_values: written,
            parts: vec![PartFrame {
                operating_pressure: op,
                source: if case.fluid.operating_pressure.is_some() {
                    "[fluid] operating_pressure".into()
                } else {
                    "default".into()
                },
                boundaries: case.boundary.keys().cloned().collect(),
                regions: Vec::new(),
                smallest: None,
            }],
        }));
        return Ok(());
    }
    // the absolute values must be > 0 (PRS001)
    for (name, bc) in &case.boundary {
        if let Some((s, _)) = pressure_of(bc)
            && let Some(v) = s.number()
            && v <= 0.0
        {
            let k = if matches!(
                bc,
                BoundaryCondition::TotalPressureInlet { .. }
                    | BoundaryCondition::TotalPressureOutlet { .. }
            ) {
                "total_pressure"
            } else {
                "pressure"
            };
            return Err(not_positive(
                &format!("[boundary.{name}] {k}"),
                v,
                format!("boundary.{name}.{k}"),
            ));
        }
    }
    if given.initial && case.initial.pressure <= 0.0 {
        return Err(not_positive(
            "[initial] pressure",
            case.initial.pressure,
            "initial.pressure".into(),
        ));
    }
    if given.reference && case.pressure.reference_value <= 0.0 {
        return Err(not_positive(
            "[pressure] reference_value",
            case.pressure.reference_value,
            "pressure.reference_value".into(),
        ));
    }
    let groups = groups(case)?;
    let mut parts = Vec::with_capacity(groups.len());
    for g in &groups {
        let region_op = g
            .regions
            .iter()
            .filter_map(|r| case.regions.get(r).and_then(|r| r.operating_pressure))
            .next();
        let mut best: Option<(u8, f64, &str)> = None;
        let mut smallest: Option<f64> = None;
        let mut pressure_boundary = false;
        for b in &g.boundaries {
            let Some((s, rank)) = case.boundary.get(b).and_then(pressure_of) else {
                continue;
            };
            pressure_boundary = true;
            let Some(v) = s.number() else { continue };
            smallest = Some(smallest.map_or(v, |m: f64| m.min(v)));
            let better = match best {
                None => true,
                Some((r, x, _)) => rank < r || (rank == r && v < x),
            };
            if better {
                best = Some((rank, v, b.as_str()));
            }
        }
        let (op, source) = if let Some(p) = region_op {
            (p, "[regions] operating_pressure".to_string())
        } else if let Some(p) = case.fluid.operating_pressure {
            (p, "[fluid] operating_pressure".to_string())
        } else if let Some((_, v, b)) = best {
            (v, format!("boundary {b}"))
        } else if pressure_boundary {
            return Err(diag::err(
                Diagnostic::new(
                    diag::PRS_NO_OPERATING,
                    format!(
                        "the pressure boundaries of fluid part {} ({}) give their values as \
                         expressions only, so its operating pressure cannot be taken from them \
                         (D239): give `[fluid] operating_pressure` (or `[regions.<name>] \
                         operating_pressure`) in Pa",
                        parts.len() + 1,
                        g.boundaries.join(", ")
                    ),
                )
                .key("fluid.operating_pressure"),
            ));
        } else if given.reference {
            (
                case.pressure.reference_value,
                "[pressure] reference_value".to_string(),
            )
        } else if given.initial {
            (case.initial.pressure, "[initial] pressure".to_string())
        } else {
            (
                DEFAULT_PRESSURE,
                "default (no pressure boundary)".to_string(),
            )
        };
        if given.initial {
            smallest =
                Some(smallest.map_or(case.initial.pressure, |m| m.min(case.initial.pressure)));
        }
        parts.push(PartFrame {
            operating_pressure: op,
            source,
            boundaries: g.boundaries.clone(),
            regions: g.regions.clone(),
            smallest,
        });
    }
    let frame = Frame {
        pressure_values: written,
        parts,
    };
    // the gauge twin: every value relative to its part's p_op
    let names: Vec<String> = case.boundary.keys().cloned().collect();
    for name in names {
        let op = frame.of_boundary(&name);
        if let Some(s) = case.boundary.get_mut(&name).and_then(pressure_of_mut) {
            *s = scalar_minus(s, op);
        }
    }
    let first = frame.parts[0].operating_pressure;
    let same = frame
        .parts
        .iter()
        .all(|p| p.operating_pressure.to_bits() == first.to_bits());
    if given.initial {
        if !same {
            return Err(diag::err(
                Diagnostic::new(
                    diag::PRS_INITIAL_PARTS,
                    format!(
                        "[initial] pressure = {} Pa is one absolute start value for fluid parts \
                         of different operating pressures ({}): leave it out, so that every \
                         part starts at its own operating pressure",
                        case.initial.pressure,
                        frame
                            .parts
                            .iter()
                            .map(|p| format!("{} Pa", p.operating_pressure))
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                )
                .key("initial.pressure"),
            ));
        }
        case.initial.pressure -= first;
    }
    if given.reference {
        // the reference cell serves a part without a pressure boundary
        let op = frame
            .parts
            .iter()
            .find(|p| {
                !p.boundaries
                    .iter()
                    .any(|b| case.boundary.get(b).and_then(pressure_of).is_some())
            })
            .map_or(first, |p| p.operating_pressure);
        case.pressure.reference_value -= op;
    }
    // the operating pressure of the [fluid] (cells outside the regions, or
    // the first part) and of the regions of the other parts
    let fluid_op = frame
        .parts
        .iter()
        .find(|p| p.regions.is_empty() || p.regions.iter().any(|r| r == "fluid"))
        .map_or(first, |p| p.operating_pressure);
    case.fluid.operating_pressure = Some(fluid_op);
    for p in &frame.parts {
        if p.operating_pressure.to_bits() == fluid_op.to_bits() {
            continue;
        }
        for r in &p.regions {
            // a region in two parts of different pressures keeps the
            // [fluid] value
            let elsewhere = frame.parts.iter().any(|q| {
                q.operating_pressure.to_bits() != p.operating_pressure.to_bits()
                    && q.regions.contains(r)
            });
            if let Some(sec) = case.regions.get_mut(r)
                && sec.kind == RegionKind::Fluid
                && sec.operating_pressure.is_none()
                && !elsewhere
            {
                sec.operating_pressure = Some(p.operating_pressure);
            }
        }
    }
    case.fluid.pressure_values = Some(PressureValues::Gauge);
    case.pressure_frame = Some(Arc::new(frame));
    Ok(())
}

/// The setup findings of the pressure values (PRS002): an absolute value
/// below 1 kPa in a gas looks like a gauge value.
pub fn setup_findings(case: &FlowCase) -> Vec<crate::bc::Finding> {
    let Some(f) = case.pressure_frame.as_deref() else {
        return Vec::new();
    };
    if !f.absolute() {
        return Vec::new();
    }
    let Ok(props) = case.properties() else {
        return Vec::new();
    };
    let gas = props.density() < 50.0
        || props
            .regions
            .iter()
            .any(|r| r.2.density.value.constant().is_some_and(|d| d < 50.0));
    let mut out = Vec::new();
    for (k, p) in f.parts.iter().enumerate() {
        if let Some(s) = p.smallest
            && s < 1000.0
            && gas
        {
            out.push(crate::bc::Finding {
                level: crate::bc::Level::Warning,
                code: diag::PRS_LOOKS_GAUGE,
                boundary: None,
                message: format!(
                    "fluid part {}: an absolute pressure of {s} Pa in a gas looks like a gauge \
                     value; set `[fluid] pressure_values = \"gauge\"` if the values are gauge \
                     values",
                    k + 1
                ),
            });
        }
    }
    out
}

/// The operating pressure per region, for the report.
pub fn by_region(f: &Frame) -> BTreeMap<String, f64> {
    let mut m = BTreeMap::new();
    for p in &f.parts {
        for r in &p.regions {
            m.insert(r.clone(), p.operating_pressure);
        }
    }
    m
}
