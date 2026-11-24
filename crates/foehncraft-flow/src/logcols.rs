//! The columns of the iteration and time-step tables (`[log] columns`, D125).
//!
//! A column is named `<kind>.<target>.<quantity>`:
//!
//! | Name | Columns |
//! |---|---|
//! | `residual.momentum`, `residual.continuity` | `ResU ResV ResW`, `ResM` |
//! | `residual.scalars` | one column per passive scalar (`Res<name>`) |
//! | `residual.turbulence` | the turbulence equations (`ResK ResOm`) |
//! | `linear.momentum`, `linear.pressure` | solver iterations (`U V W`, `P`) |
//! | `linear.coupled` | Krylov iterations of the coupled algorithm (`K`) |
//! | `probe.<name>.<q>` | value at an `[[output.probes]]` point |
//! | `cell.<index>.<q>` | value in one cell (debugging) |
//! | `min.<q>`, `max.<q>` | extreme over the domain and the cell it occurs in |
//! | `boundary.<name>.mass_flow` | mass flow out of a boundary (kg/s) |
//! | `rotor.<zone>.torque`, `.power` | torque of the fluid on a rotor (N m), its shaft power (W, D184) |
//! | `machine.<name>.flow`, `.pressure_rise`, `.head`, `.efficiency` | an `[[output.machines]]` entry (D184) |
//! | `time.simulated`, `.dt`, `.courant`, `.outer`, `.error` | transient runs |
//! | `time.elapsed`, `time.per_iteration` | wall-clock seconds |
//!
//! `<q>` is `velocity_x`, `velocity_y`, `velocity_z`, `velocity_magnitude`
//! or `pressure`. One formatter turns these into the table of the console
//! and `flow.log` and into the `columns` and `table_row` rows of `run.jsonl`,
//! so a GUI shows what the terminal shows. Evaluating a column only reads
//! the fields; no result changes with the choice of columns.

use crate::console::{Col, Group, Style, Table, sci};
use crate::diag::CASE_BAD_VALUE;
use crate::simplec::{Residuals, Solution};
use anyhow::Result;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_mesh::locate::Locator;
use foehncraft_mesh::{Geometry, Mesh};
use serde_json::{Value, json};
use std::ops::Range;
use std::time::Instant;

/// A field quantity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quantity {
    /// `velocity_x`.
    Ux,
    /// `velocity_y`.
    Uy,
    /// `velocity_z`.
    Uz,
    /// `velocity_magnitude`.
    Umag,
    /// `pressure`.
    P,
}

impl Quantity {
    const ALL: [(&'static str, Quantity); 5] = [
        ("velocity_x", Quantity::Ux),
        ("velocity_y", Quantity::Uy),
        ("velocity_z", Quantity::Uz),
        ("velocity_magnitude", Quantity::Umag),
        ("pressure", Quantity::P),
    ];

    fn parse(s: &str) -> Option<Self> {
        Self::ALL.iter().find(|(n, _)| *n == s).map(|(_, q)| *q)
    }

    /// Short header name.
    fn short(self) -> &'static str {
        match self {
            Self::Ux => "Ux",
            Self::Uy => "Uy",
            Self::Uz => "Uz",
            Self::Umag => "U",
            Self::P => "P",
        }
    }

    fn unit(self) -> &'static str {
        if self == Self::P { "Pa" } else { "m/s" }
    }

    fn of(self, u: [f64; 3], p: f64) -> f64 {
        match self {
            Self::Ux => u[0],
            Self::Uy => u[1],
            Self::Uz => u[2],
            Self::Umag => (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt(),
            Self::P => p,
        }
    }
}

/// One selectable column set.
#[derive(Debug, Clone, PartialEq)]
pub enum Spec {
    /// `residual.momentum`.
    ResMomentum,
    /// `residual.continuity`.
    ResContinuity,
    /// `residual.scalars`: the names of the case's scalars (filled in by
    /// [`Columns::build`]).
    ResScalars(Vec<String>),
    /// `residual.turbulence`: the residuals of the turbulence equations,
    /// by name (k and ω, k and ε, or k, ε, f and ζ; filled in by
    /// [`Columns::build`]).
    ResTurbulence(Vec<&'static str>),
    /// `linear.momentum`.
    LinMomentum,
    /// `linear.pressure`.
    LinPressure,
    /// `linear.coupled`.
    LinCoupled,
    /// `probe.<name>.<q>`.
    Probe(String, Quantity),
    /// `cell.<index>.<q>`.
    Cell(usize, Quantity),
    /// `min.<q>`.
    Min(Quantity),
    /// `max.<q>`.
    Max(Quantity),
    /// `boundary.<name>.mass_flow`.
    BoundaryMassFlow(String),
    /// `rotor.<zone>.<quantity>` (D184).
    Rotor(String, RotorQ),
    /// `machine.<name>.<quantity>` (D184).
    Machine(String, MachineQ),
    /// `time.simulated`.
    Simulated,
    /// `time.dt`.
    Dt,
    /// `time.courant`.
    Courant,
    /// `time.outer`.
    Outer,
    /// `time.error`.
    Error,
    /// `time.elapsed`.
    Elapsed,
    /// `time.per_iteration`.
    PerIteration,
}

/// A quantity of a rotor (D184).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RotorQ {
    /// `torque`: of the fluid on its walls (N m).
    Torque,
    /// `power`: the shaft power (W).
    Power,
    /// `angle`: the angle the rotor has turned since the start, |Ω| t
    /// (degrees; transient runs, the sliding mesh).
    Angle,
}

/// A quantity of a machine of `[[output.machines]]` (D184).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MachineQ {
    /// `flow`: volume flow (m³/s).
    Flow,
    /// `pressure_rise`: total pressure rise (Pa).
    PressureRise,
    /// `head`: H (m).
    Head,
    /// `efficiency`: hydraulic efficiency (-).
    Efficiency,
}

impl MachineQ {
    fn parse(s: &str) -> Option<Self> {
        Some(match s {
            "flow" => Self::Flow,
            "pressure_rise" => Self::PressureRise,
            "head" => Self::Head,
            "efficiency" => Self::Efficiency,
            _ => return None,
        })
    }

    fn id(self) -> &'static str {
        match self {
            Self::Flow => "flow",
            Self::PressureRise => "pressure_rise",
            Self::Head => "head",
            Self::Efficiency => "efficiency",
        }
    }
}

/// The names that `[log] columns` accepts, for the schema and the hints.
pub const KINDS: &[&str] = &[
    "residual.momentum",
    "residual.continuity",
    "residual.scalars",
    "residual.turbulence",
    "linear.momentum",
    "linear.coupled",
    "linear.pressure",
    "probe.<name>.<quantity>",
    "cell.<index>.<quantity>",
    "min.<quantity>",
    "max.<quantity>",
    "boundary.<name>.mass_flow",
    "rotor.<zone>.<torque|power|angle>",
    "machine.<name>.<flow|pressure_rise|head|efficiency>",
    "time.simulated",
    "time.dt",
    "time.courant",
    "time.outer",
    "time.error",
    "time.elapsed",
    "time.per_iteration",
];

fn bad(name: &str, why: &str) -> anyhow::Error {
    let valid = if why.contains("quantity") {
        format!("quantities: {}", Quantity::ALL.map(|(n, _)| n).join(", "))
    } else {
        format!("valid names: {}", KINDS.join(", "))
    };
    crate::diag::err(
        crate::diag::Diagnostic::new(
            CASE_BAD_VALUE,
            format!("[log] columns: `{name}`: {why}; {valid}"),
        )
        .key("log.columns"),
    )
}

impl Spec {
    /// Parses a column name.
    pub fn parse(name: &str) -> Result<Self> {
        let parts: Vec<&str> = name.split('.').collect();
        let quantity = |s: &str| Quantity::parse(s).ok_or_else(|| bad(name, "unknown quantity"));
        // `<kind>.<target>.<last>` with a target that may contain dots
        let target = || parts[1..parts.len() - 1].join(".");
        Ok(match (parts[0], parts.len()) {
            ("residual", 2) if parts[1] == "momentum" => Self::ResMomentum,
            ("residual", 2) if parts[1] == "continuity" => Self::ResContinuity,
            ("residual", 2) if parts[1] == "scalars" => Self::ResScalars(Vec::new()),
            ("residual", 2) if parts[1] == "turbulence" => Self::ResTurbulence(Vec::new()),
            ("linear", 2) if parts[1] == "momentum" => Self::LinMomentum,
            ("linear", 2) if parts[1] == "pressure" => Self::LinPressure,
            ("linear", 2) if parts[1] == "coupled" => Self::LinCoupled,
            ("min", 2) => Self::Min(quantity(parts[1])?),
            ("max", 2) => Self::Max(quantity(parts[1])?),
            ("probe", n) if n >= 3 => Self::Probe(target(), quantity(parts[n - 1])?),
            ("cell", n) if n >= 3 => {
                let i = target()
                    .parse::<usize>()
                    .map_err(|_| bad(name, "the cell index is not a number"))?;
                Self::Cell(i, quantity(parts[n - 1])?)
            }
            ("boundary", n) if n >= 3 && parts[n - 1] == "mass_flow" => {
                Self::BoundaryMassFlow(target())
            }
            ("rotor", n) if n >= 3 => match parts[n - 1] {
                "torque" => Self::Rotor(target(), RotorQ::Torque),
                "power" => Self::Rotor(target(), RotorQ::Power),
                "angle" => Self::Rotor(target(), RotorQ::Angle),
                _ => {
                    return Err(bad(name, "a rotor column is `torque`, `power` or `angle`"));
                }
            },
            ("machine", n) if n >= 3 => match MachineQ::parse(parts[n - 1]) {
                Some(q) => Self::Machine(target(), q),
                None => {
                    return Err(bad(
                        name,
                        "a machine column is `flow`, `pressure_rise`, `head` or `efficiency`",
                    ));
                }
            },
            ("time", 2) => match parts[1] {
                "simulated" => Self::Simulated,
                "dt" => Self::Dt,
                "courant" => Self::Courant,
                "outer" => Self::Outer,
                "error" => Self::Error,
                "elapsed" => Self::Elapsed,
                "per_iteration" => Self::PerIteration,
                _ => return Err(bad(name, "unknown time column")),
            },
            _ => return Err(bad(name, "unknown column")),
        })
    }

    fn transient_only(&self) -> bool {
        matches!(
            self,
            Self::Simulated | Self::Dt | Self::Courant | Self::Outer | Self::Error
        )
    }

    fn needs_state(&self) -> bool {
        matches!(
            self,
            Self::Probe(..)
                | Self::Cell(..)
                | Self::Min(_)
                | Self::Max(_)
                | Self::BoundaryMassFlow(_)
                | Self::Rotor(..)
                | Self::Machine(..)
        )
    }

    fn needs_rotors(&self) -> bool {
        matches!(self, Self::Rotor(..) | Self::Machine(..))
    }
}

/// Checks the names of `[log] columns` against the case (syntax, probe and
/// boundary names, and the kind of run); the cell indices are checked
/// against the mesh when the run starts.
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let Some(names) = &case.log.columns else {
        return Ok(());
    };
    let transient = case.time.is_some();
    for n in names {
        let spec = Spec::parse(n)?;
        match &spec {
            Spec::Probe(p, _) => crate::need!(
                case.output.probes.iter().any(|q| &q.name == p),
                CASE_BAD_VALUE,
                "[log] columns: `{n}`: no probe `{p}` in [[output.probes]]"
            ),
            Spec::BoundaryMassFlow(b) => crate::need!(
                case.boundary.contains_key(b),
                CASE_BAD_VALUE,
                "[log] columns: `{n}`: no boundary `{b}` in the case"
            ),
            Spec::Rotor(r, _) => crate::need!(
                rotor_names(case).contains(r),
                CASE_BAD_VALUE,
                "[log] columns: `{n}`: no turning zone `{r}` ([[mrf]] or [[sliding]] zones: {})",
                rotor_names(case).join(", ")
            ),
            Spec::Machine(m, _) => crate::need!(
                case.output.machines.iter().any(|x| &x.name == m),
                CASE_BAD_VALUE,
                "[log] columns: `{n}`: no `{m}` in [[output.machines]]"
            ),
            s if s.transient_only() => crate::need!(
                transient,
                CASE_BAD_VALUE,
                "[log] columns: `{n}` needs a transient run ([time] table)"
            ),
            _ => {}
        }
    }
    crate::need!(
        (1..=8).contains(&case.log.digits),
        CASE_RANGE,
        "[log] digits must be between 1 and 8"
    );
    Ok(())
}

/// The names of the rotors of `case` (D184).
fn rotor_names(case: &crate::case::FlowCase) -> Vec<String> {
    case.mrf
        .iter()
        .map(|m| crate::rotor::rotor_name(&m.zone))
        .chain(
            case.sliding
                .iter()
                .map(|s| crate::rotor::rotor_name(&s.zone)),
        )
        .collect()
}

/// One column of the table.
#[derive(Debug, Clone)]
pub struct Column {
    /// Name in the case and the run record.
    pub id: String,
    /// Short name in the header row.
    pub name: String,
    /// Group title above the column.
    pub group: String,
    /// Unit (empty: none).
    pub unit: &'static str,
    /// Width in the table.
    pub width: usize,
}

/// A table cell: the text of the console and the value of the record.
#[derive(Debug, Clone)]
pub struct Entry {
    /// As printed.
    pub text: String,
    /// As recorded.
    pub value: Value,
}

/// Time quantities of a step.
#[derive(Debug, Clone, Copy, Default)]
pub struct TimeInfo {
    /// Simulated time (s).
    pub time: f64,
    /// Time step (s).
    pub dt: f64,
    /// Courant number.
    pub courant: f64,
    /// Outer iterations of the step.
    pub outer: usize,
    /// Error of the step control, if one is set.
    pub error: Option<f64>,
}

/// What a row is evaluated from.
pub struct RowCtx<'a> {
    /// Iteration or step.
    pub index: usize,
    /// Residuals and linear iterations.
    pub res: &'a Residuals,
    /// The step, in a transient run.
    pub time: Option<TimeInfo>,
    /// The fields (needed when [`Columns::needs_state`]).
    pub sol: Option<&'a Solution<'a>>,
    /// Torques, powers and machine values (needed when
    /// [`Columns::needs_rotors`]).
    pub rotors: Option<&'a crate::rotor::Values>,
}

/// What [`refresh`] needs to build the columns again.
pub struct Refresh<'a> {
    /// The case in force.
    pub case: &'a crate::case::FlowCase,
    /// The mesh.
    pub mesh: &'a Mesh,
    /// Its geometry.
    pub geo: &'a Geometry,
    /// The cell finder of the probes.
    pub loc: &'a Locator<'a>,
    /// A transient run.
    pub transient: bool,
    /// `iteration` or `step` (the `columns` row of the run record).
    pub table_name: &'a str,
    /// The title of the counter column.
    pub first: &'a str,
    /// The iteration or step.
    pub at: usize,
}

/// Builds the columns and the table again when `log.columns` changed during
/// the run (D214: `[[schedule]]`, `run.set`, a re-read): the next row
/// prints the header again and a `columns` row goes to the run record.
/// Columns that cannot be built (a cell the mesh does not have) are
/// reported through `warn` and the table stays as it was.
pub fn refresh(
    cols: &mut Columns,
    table: &mut Table,
    r: &mut Refresh,
    warn: &mut dyn FnMut(&str),
    steer: &mut crate::engine::Steer,
) {
    if cols.names == r.case.log.columns {
        return;
    }
    match Columns::build(r.case, r.mesh, r.geo, r.loc, r.transient) {
        Ok(c) => {
            *cols = c;
            *table = cols.table(r.first, r.case.log.style, r.case.log.header_every);
            if steer.record.on() {
                steer.record.row(&cols.record_header(r.table_name, r.first));
            }
        }
        Err(e) => {
            warn(&format!(
                "{} {}: the new log.columns are not used ({e:#}); the table keeps its columns",
                r.table_name, r.at
            ));
            cols.names = r.case.log.columns.clone();
        }
    }
}

/// The table columns of a run.
pub struct Columns {
    /// `[log] columns` as the columns were built from (D214: a change
    /// during the run rebuilds them).
    names: Option<Vec<String>>,
    specs: Vec<Spec>,
    columns: Vec<Column>,
    digits: usize,
    transient: bool,
    probes: Option<crate::monitor::Probes>,
    boundaries: Vec<(String, Range<usize>)>,
    start: Instant,
    last: Option<(Instant, usize)>,
}

impl Columns {
    /// The columns of the case (`[log] columns`, or the default set) for
    /// this mesh. The Locator finds the cells of the probes.
    pub fn build(
        case: &crate::case::FlowCase,
        mesh: &Mesh,
        geo: &Geometry,
        loc: &Locator,
        transient: bool,
    ) -> Result<Self> {
        let digits = case.log.digits;
        let with_error = case.time.is_some_and(|t| t.error_tolerance.is_some());
        let mut specs: Vec<Spec> = match &case.log.columns {
            Some(names) => names
                .iter()
                .map(|n| Spec::parse(n))
                .collect::<Result<_>>()?,
            None if transient => {
                let mut v = vec![
                    Spec::Simulated,
                    Spec::Dt,
                    Spec::Courant,
                    Spec::ResMomentum,
                    Spec::ResContinuity,
                    Spec::Outer,
                    Spec::LinPressure,
                ];
                if with_error {
                    v.push(Spec::Error);
                }
                v
            }
            None if case.solver.algorithm(false) == crate::simplec::Algorithm::Coupled => vec![
                Spec::ResMomentum,
                Spec::ResContinuity,
                Spec::LinCoupled,
                Spec::LinMomentum,
                Spec::LinPressure,
            ],
            None => vec![
                Spec::ResMomentum,
                Spec::ResContinuity,
                Spec::LinMomentum,
                Spec::LinPressure,
            ],
        };
        // a turning zone: the torque of its rotor, the head and efficiency
        // of the machines (D184)
        if case.log.columns.is_none() {
            specs.extend(
                rotor_names(case)
                    .into_iter()
                    .map(|r| Spec::Rotor(r, RotorQ::Torque)),
            );
            // the angle of every sliding rotor (D184): they may turn at
            // different speeds
            if transient {
                specs.extend(
                    case.sliding
                        .iter()
                        .map(|m| Spec::Rotor(crate::rotor::rotor_name(&m.zone), RotorQ::Angle)),
                );
            }
            for m in &case.output.machines {
                specs.push(Spec::Machine(m.name.clone(), MachineQ::Head));
                specs.push(Spec::Machine(m.name.clone(), MachineQ::Efficiency));
            }
        }
        if case.log.columns.is_none() && !case.turbulence.is_laminar() {
            let at = specs
                .iter()
                .position(|s| *s == Spec::ResContinuity)
                .map_or(specs.len(), |i| i + 1);
            specs.insert(at, Spec::ResTurbulence(Vec::new()));
        }
        let names: Vec<String> = case.scalars.keys().cloned().collect();
        if case.log.columns.is_none() && !names.is_empty() {
            let at = specs
                .iter()
                .position(|s| *s == Spec::ResContinuity)
                .map_or(specs.len(), |i| i + 1);
            specs.insert(at, Spec::ResScalars(Vec::new()));
        }
        for s in &mut specs {
            if let Spec::ResScalars(v) = s {
                v.clone_from(&names);
            }
        }
        for s in &specs {
            if let Spec::Cell(i, _) = s {
                crate::need!(
                    *i < mesh.n_cells,
                    CASE_RANGE,
                    "[log] columns: cell {i} does not exist (the mesh has {} cells)",
                    mesh.n_cells
                );
            }
            if s.transient_only() {
                crate::need!(
                    transient,
                    CASE_BAD_VALUE,
                    "[log] columns: a time column needs a transient run"
                );
            }
        }
        let eqs = crate::turbulence::equations(case.turbulence.model);
        for s in &mut specs {
            if let Spec::ResTurbulence(e) = s {
                *e = eqs.to_vec();
            }
        }
        let probes = specs
            .iter()
            .any(|s| matches!(s, Spec::Probe(..)))
            .then(|| crate::monitor::Probes::new(case, &[], loc, geo));
        let boundaries = mesh
            .boundaries
            .iter()
            .map(|p| (p.name.clone(), p.range()))
            .collect();
        let columns = specs
            .iter()
            .flat_map(|s| layout(s, digits, transient))
            .collect();
        Ok(Self {
            names: case.log.columns.clone(),
            specs,
            columns,
            digits,
            transient,
            probes,
            boundaries,
            start: Instant::now(),
            last: None,
        })
    }

    /// Whether a row needs the fields (and so a download from the device).
    pub fn needs_state(&self) -> bool {
        self.specs.iter().any(Spec::needs_state)
    }

    /// Whether a row needs the torques and machine values (D184).
    pub fn needs_rotors(&self) -> bool {
        self.specs.iter().any(Spec::needs_rotors)
    }

    /// The columns, in order.
    pub fn columns(&self) -> &[Column] {
        &self.columns
    }

    /// The table: a counter column and groups of adjacent columns with the
    /// same group title.
    pub fn table(&self, first: &str, style: Style, header_every: usize) -> Table {
        let mut groups: Vec<Group> = Vec::new();
        for c in &self.columns {
            let col = Col::int(&c.name, c.width);
            match groups.last_mut() {
                Some(g) if g.title == c.group => g.cols.push(col),
                _ => groups.push(Group::new(&c.group, vec![col])),
            }
        }
        // a group title that does not fit above its columns widens the last
        // one (the default groups all fit)
        for g in &mut groups {
            let need = g.title.chars().count() + 4;
            let have: usize = g.cols.iter().map(|c| c.width + 1).sum::<usize>() + 1;
            if let Some(last) = g.cols.last_mut()
                && have < need
            {
                last.width += need - have;
            }
        }
        Table::new(Col::int(first, 6), groups, style, header_every)
    }

    /// The `columns` row of the run record.
    pub fn record_header(&self, table: &str, first: &str) -> Value {
        json!({
            "type": "columns",
            "table": table,
            "first": first,
            "columns": self.columns.iter().map(|c| json!({
                "id": c.id, "name": c.name, "group": c.group, "unit": c.unit,
            })).collect::<Vec<_>>(),
        })
    }

    /// The `table_row` row of the run record.
    pub fn record_row(&self, table: &str, index: usize, entries: &[Entry]) -> Value {
        json!({
            "type": "table_row",
            "table": table,
            "index": index,
            "values": self.columns.iter().zip(entries)
                .map(|(c, e)| (c.id.clone(), e.value.clone()))
                .collect::<serde_json::Map<_, _>>(),
        })
    }

    /// The cells of a row.
    pub fn eval(&mut self, ctx: &RowCtx, geo: &Geometry) -> Vec<Entry> {
        let d = self.digits;
        let num = |v: f64| Entry {
            text: sci(v, d - 1),
            value: json_num(v),
        };
        let int = |n: usize| Entry {
            text: n.to_string(),
            value: json!(n),
        };
        let now = Instant::now();
        let elapsed = now.duration_since(self.start).as_secs_f64();
        let per_it = self.last.map_or(f64::NAN, |(t, i)| {
            let n = ctx.index.saturating_sub(i).max(1);
            now.duration_since(t).as_secs_f64() / n as f64
        });
        self.last = Some((now, ctx.index));
        let mut out = Vec::with_capacity(self.columns.len());
        let ti = ctx.time.unwrap_or_default();
        for s in &self.specs {
            match s {
                Spec::ResMomentum => out.extend(ctx.res.momentum.iter().map(|&r| num(r))),
                Spec::ResContinuity => out.push(num(ctx.res.continuity)),
                Spec::ResScalars(names) => out
                    .extend((0..names.len()).map(|i| {
                        num(ctx.res.scalars.residual.get(i).copied().unwrap_or(f64::NAN))
                    })),
                Spec::ResTurbulence(names) => out.extend((0..names.len()).map(|i| {
                    num(ctx
                        .res
                        .turbulence
                        .residual
                        .get(i)
                        .copied()
                        .unwrap_or(f64::NAN))
                })),
                Spec::LinMomentum => {
                    out.extend(ctx.res.momentum_iterations.iter().map(|&n| int(n)))
                }
                Spec::LinPressure => out.push(int(ctx.res.pressure_iterations)),
                Spec::LinCoupled => out.push(int(ctx.res.coupled_iterations)),
                Spec::Probe(name, q) => {
                    let v = ctx.sol.zip(self.probes.as_ref()).and_then(|(sol, p)| {
                        let i = p.index(name)?;
                        let (u, pr) = p.values(sol, geo)[i];
                        Some(q.of(u, pr))
                    });
                    out.push(v.map_or_else(dash, num));
                }
                Spec::Cell(c, q) => {
                    let v = ctx.sol.map(|s| q.of(s.u.get(*c), s.p[*c]));
                    out.push(v.map_or_else(dash, num));
                }
                Spec::Min(q) | Spec::Max(q) => {
                    let max = matches!(s, Spec::Max(_));
                    match ctx.sol {
                        Some(sol) => {
                            let (v, c) = extreme(sol, *q, max, geo.cell_volume.len());
                            out.push(num(v));
                            out.push(Entry {
                                text: c.to_string(),
                                value: json!({"cell": c, "x": geo.cell_centroid[c]}),
                            });
                        }
                        None => out.extend([dash(), dash()]),
                    }
                }
                Spec::BoundaryMassFlow(name) => {
                    let v = ctx.sol.and_then(|sol| {
                        let (_, r) = self.boundaries.iter().find(|(n, _)| n == name)?;
                        Some(r.clone().map(|f| sol.mass_flux[f]).exact_sum())
                    });
                    out.push(v.map_or_else(dash, num));
                }
                Spec::Rotor(name, q) => {
                    let v = ctx
                        .rotors
                        .and_then(|v| v.rotor(name))
                        .and_then(|r| match q {
                            RotorQ::Torque => Some(r.torque),
                            RotorQ::Power => Some(r.power),
                            RotorQ::Angle => ctx.time.map(|t| (r.omega * t.time).to_degrees()),
                        });
                    out.push(v.map_or_else(dash, |v| {
                        if *q == RotorQ::Angle {
                            Entry {
                                text: format!("{v:.2}"),
                                value: json_num(v),
                            }
                        } else {
                            num(v)
                        }
                    }));
                }
                Spec::Machine(name, q) => {
                    let v = ctx.rotors.and_then(|v| v.machine(name)).map(|m| match q {
                        MachineQ::Flow => m.volume_flow,
                        MachineQ::PressureRise => m.total_pressure_rise,
                        MachineQ::Head => m.head,
                        MachineQ::Efficiency => m.efficiency,
                    });
                    out.push(v.map_or_else(dash, num));
                }
                Spec::Simulated => out.push(Entry {
                    text: sci(ti.time, d + 1),
                    value: json_num(ti.time),
                }),
                Spec::Dt => out.push(num(ti.dt)),
                Spec::Courant => out.push(Entry {
                    text: format!("{:.3}", ti.courant),
                    value: json_num(ti.courant),
                }),
                Spec::Outer => out.push(int(ti.outer)),
                Spec::Error => out.push(match ti.error {
                    Some(e) => Entry {
                        text: format!("{e:.2}"),
                        value: json_num(e),
                    },
                    None => dash(),
                }),
                Spec::Elapsed => out.push(Entry {
                    text: format!("{elapsed:.1}"),
                    value: json_num(elapsed),
                }),
                Spec::PerIteration => out.push(if per_it.is_nan() { dash() } else { num(per_it) }),
            }
        }
        let _ = self.transient;
        out
    }
}

fn dash() -> Entry {
    Entry {
        text: "-".into(),
        value: Value::Null,
    }
}

fn json_num(v: f64) -> Value {
    if v.is_finite() {
        json!(v)
    } else {
        // JSON has no NaN or Inf
        Value::String(sci(v, 2))
    }
}

/// Extreme of a quantity over the cells and the first cell it occurs in
/// (NaN cells are skipped; all NaN gives NaN in cell 0).
fn extreme(sol: &Solution, q: Quantity, max: bool, n: usize) -> (f64, usize) {
    let mut best = (f64::NAN, 0usize);
    for c in 0..n {
        let v = q.of(sol.u.get(c), sol.p[c]);
        if v.is_nan() {
            continue;
        }
        if best.0.is_nan() || (max && v > best.0) || (!max && v < best.0) {
            best = (v, c);
        }
    }
    best
}

/// The table columns of one spec.
fn layout(s: &Spec, d: usize, transient: bool) -> Vec<Column> {
    let w = d + 5;
    let ws = d + 6; // with room for a sign
    let col = |id: String, name: &str, group: &str, unit: &'static str, width: usize| Column {
        id,
        name: name.into(),
        group: group.into(),
        unit,
        width,
    };
    match s {
        Spec::ResMomentum => ["U", "V", "W"]
            .iter()
            .map(|c| {
                col(
                    format!("residual.momentum.{}", c.to_lowercase()),
                    &format!("Res{c}"),
                    "scaled residuals",
                    "",
                    w,
                )
            })
            .collect(),
        Spec::ResContinuity => vec![col(
            "residual.continuity".into(),
            "ResM",
            "scaled residuals",
            "",
            w,
        )],
        Spec::ResScalars(names) => names
            .iter()
            .map(|n| {
                let short: String = format!("Res{n}").chars().take(w).collect();
                col(
                    format!("residual.scalars.{n}"),
                    &short,
                    "scaled residuals",
                    "",
                    w,
                )
            })
            .collect(),
        Spec::ResTurbulence(names) => names
            .iter()
            .map(|&id| {
                (
                    id,
                    match id {
                        "k" => "ResK",
                        "epsilon" => "ResEps",
                        "omega" => "ResOm",
                        "zeta" => "ResZeta",
                        "phi" => "ResPhi",
                        "alpha" => "ResAlpha",
                        _ => "ResF",
                    },
                )
            })
            .map(|(id, n)| {
                col(
                    format!("residual.turbulence.{id}"),
                    n,
                    "scaled residuals",
                    "",
                    w,
                )
            })
            .collect(),
        Spec::LinMomentum => {
            let g = if transient {
                "solvers"
            } else {
                "linear iterations"
            };
            ["U", "V", "W"]
                .iter()
                .map(|c| col(format!("linear.momentum.{}", c.to_lowercase()), c, g, "", 4))
                .collect()
        }
        Spec::LinPressure => {
            if transient {
                vec![col("linear.pressure".into(), "CG", "solvers", "", 5)]
            } else {
                vec![col(
                    "linear.pressure".into(),
                    "P",
                    "linear iterations",
                    "",
                    5,
                )]
            }
        }
        Spec::LinCoupled => vec![col(
            "linear.coupled".into(),
            "K",
            if transient {
                "solvers"
            } else {
                "linear iterations"
            },
            "",
            4,
        )],
        Spec::Probe(n, q) => vec![col(
            format!("probe.{n}.{}", name_of(*q)),
            q.short(),
            &format!("probe {n}"),
            q.unit(),
            ws,
        )],
        Spec::Cell(i, q) => vec![col(
            format!("cell.{i}.{}", name_of(*q)),
            q.short(),
            &format!("cell {i}"),
            q.unit(),
            ws,
        )],
        Spec::Min(q) | Spec::Max(q) => {
            let (k, suffix) = if matches!(s, Spec::Max(_)) {
                ("max", "max")
            } else {
                ("min", "min")
            };
            vec![
                col(
                    format!("{k}.{}", name_of(*q)),
                    &format!("{}{suffix}", q.short()),
                    "extremes",
                    q.unit(),
                    ws,
                ),
                col(
                    format!("{k}.{}.cell", name_of(*q)),
                    "cell",
                    "extremes",
                    "",
                    8,
                ),
            ]
        }
        Spec::BoundaryMassFlow(n) => vec![col(
            format!("boundary.{n}.mass_flow"),
            &n.chars().take(14).collect::<String>(),
            "flow kg/s",
            "kg/s",
            ws.max(n.chars().count().min(14)),
        )],
        Spec::Rotor(n, q) => {
            let (id, name, unit) = match q {
                RotorQ::Torque => ("torque", "Torque", "N m"),
                RotorQ::Power => ("power", "Power", "W"),
                RotorQ::Angle => ("angle", "Angle", "deg"),
            };
            // (the zone in the heading: several rotors turn at their own
            // speeds, D185)
            let head = format!("{name} {}", n.chars().take(14).collect::<String>());
            let w = ws.max(head.chars().count());
            vec![col(
                format!("rotor.{n}.{id}"),
                &head,
                &format!("rotor {n}"),
                unit,
                w,
            )]
        }
        Spec::Machine(n, q) => {
            let (name, unit) = match q {
                MachineQ::Flow => ("Q", "m³/s"),
                MachineQ::PressureRise => ("dpt", "Pa"),
                MachineQ::Head => ("H", "m"),
                MachineQ::Efficiency => ("eta", "-"),
            };
            vec![col(
                format!("machine.{n}.{}", q.id()),
                name,
                &format!("machine {n}"),
                unit,
                ws,
            )]
        }
        Spec::Simulated => vec![col("time.simulated".into(), "Time [s]", "time", "s", d + 7)],
        Spec::Dt => vec![col("time.dt".into(), "Dt [s]", "time", "s", w)],
        Spec::Courant => vec![col("time.courant".into(), "Cour.", "time", "", 5)],
        Spec::Outer => vec![col("time.outer".into(), "Outer", "solvers", "", 5)],
        Spec::Error => vec![col("time.error".into(), "Error", "solvers", "", 6)],
        Spec::Elapsed => vec![col("time.elapsed".into(), "Elapsed", "run [s]", "s", 9)],
        Spec::PerIteration => vec![col(
            "time.per_iteration".into(),
            "s/iter",
            "run [s]",
            "s",
            w,
        )],
    }
}

fn name_of(q: Quantity) -> &'static str {
    Quantity::ALL
        .iter()
        .find(|(_, x)| *x == q)
        .map_or("", |(n, _)| n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_parse() {
        assert_eq!(
            Spec::parse("probe.p.one.velocity_x").unwrap(),
            Spec::Probe("p.one".into(), Quantity::Ux)
        );
        assert_eq!(
            Spec::parse("cell.17.pressure").unwrap(),
            Spec::Cell(17, Quantity::P)
        );
        assert_eq!(
            Spec::parse("max.velocity_magnitude").unwrap(),
            Spec::Max(Quantity::Umag)
        );
        assert_eq!(
            Spec::parse("boundary.outlet.mass_flow").unwrap(),
            Spec::BoundaryMassFlow("outlet".into())
        );
        assert!(Spec::parse("max.speed").is_err());
        assert!(Spec::parse("cell.x.pressure").is_err());
        assert!(Spec::parse("residual").is_err());
        assert!(Spec::parse("nothing.here").is_err());
    }

    #[test]
    fn the_default_steady_layout_is_the_classic_table() {
        let c: Vec<Column> = [
            Spec::ResMomentum,
            Spec::ResContinuity,
            Spec::LinMomentum,
            Spec::LinPressure,
        ]
        .iter()
        .flat_map(|s| layout(s, 3, false))
        .collect();
        let names: Vec<&str> = c.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, ["ResU", "ResV", "ResW", "ResM", "U", "V", "W", "P"]);
        assert!(c[..4].iter().all(|c| c.width == 8));
    }
}
