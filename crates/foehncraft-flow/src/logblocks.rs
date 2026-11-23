//! The framed blocks of the log (D43): what was run, with which files, on
//! which mesh and fluid, with which settings and boundary conditions, and
//! at the end what came out. The text comes from the case and the mesh;
//! nothing here computes a result.

use crate::bc::BoundaryCondition;
use crate::case::FlowCase;
use crate::console::{Frame, clock, count, hostname, pretty_date, sci, sci3, sci5, size};
use crate::engine::BuildInfo;
use foehncraft_mesh::quality::QualityReport;
use serde::Serialize;
use std::path::Path;

/// A commit as seven characters, with `-dirty` kept.
pub fn short_commit(c: &str) -> String {
    match c.strip_suffix("-dirty") {
        Some(h) => format!("{}-dirty", h.chars().take(7).collect::<String>()),
        None => c.chars().take(7).collect(),
    }
}

/// The name of a value as the case file spells it (`linear-upwind`).
pub fn spelled<T: Serialize>(v: &T) -> String {
    serde_json::to_value(v)
        .ok()
        .and_then(|j| j.as_str().map(str::to_string))
        .unwrap_or_default()
}

pub use crate::banner_art::BANNER_ART;

/// The version, the date and the host.
pub fn banner(build: Option<&BuildInfo>, case: &FlowCase) -> Frame {
    let mut f = Frame::new("Foehncraft").art(&BANNER_ART).kv(
        "Version",
        match build {
            Some(b) => format!(
                "{}  (build {})",
                foehncraft_core::VERSION,
                short_commit(&b.commit)
            ),
            None => foehncraft_core::VERSION.to_string(),
        },
    );
    f.push_kv("Date and time", pretty_date(&crate::engine::run_date()));
    f.push_kv("Host name", hostname());
    if let Some(b) = build {
        f.push_kv("Platform", format!("{}, {}", b.platform, b.cpu));
        if !b.features.is_empty() {
            f.push_kv("Features", b.features.join(", "));
        }
    }
    f.push_kv(
        "Threads / backend",
        format!(
            "{} / {}",
            case.run.threads.unwrap_or_else(rayon::current_num_threads),
            case.run.backend
        ),
    );
    f
}

/// The input files and the output directory; the driver appends the files
/// it will write.
pub fn files(case: &FlowCase, desc: &str, out: &Path) -> Frame {
    let mut f = Frame::new("Input and output files").kv(
        "Case file",
        case.source
            .as_ref()
            .map_or("none (built from a string)".to_string(), |p| {
                p.display().to_string()
            }),
    );
    f.push_kv(
        "Overrides (--set)",
        if case.overrides.is_empty() {
            "none".to_string()
        } else {
            case.overrides.join("\n")
        },
    );
    f.push_kv("Mesh", desc);
    f.push_kv(
        "Material",
        case.fluid
            .material
            .clone()
            .unwrap_or_else(|| "constant properties of the case file".into()),
    );
    f.push_kv(
        "Restart from",
        case.run
            .restart
            .as_ref()
            .map_or("none".to_string(), |p| p.display().to_string()),
    );
    for (name, t) in &case.tables {
        let _ = t;
        f.push_kv("Table", name);
    }
    f.push_rule();
    f.push_kv("Output directory", out.display());
    f.push_kv("Resolved case", "case.resolved.toml");
    f.push_kv("Log / run record", "flow.log / run.jsonl");
    f
}

/// The flow category of the run-control frame: laminar (Stokes without
/// convection) or turbulent with the model and its wall treatment.
fn flow_category(case: &FlowCase) -> String {
    let t = &case.turbulence;
    if !t.is_laminar() {
        format!(
            "incompressible, turbulent ({}, {})",
            t.model.name(),
            t.treatment().name()
        )
    } else if case.flow.convection {
        "incompressible, laminar".to_string()
    } else {
        "incompressible, laminar, Stokes (convection off)".to_string()
    }
}

/// What a run of a steady case does.
pub fn control_steady(case: &FlowCase) -> Frame {
    let s = &case.solver;
    let alg = spelled(&s.algorithm(false)).to_uppercase();
    let mut f = Frame::new("Run control")
        .kv("Run mode", "steady")
        .kv("Flow category", flow_category(case))
        .kv("Algorithm", alg)
        .kv("Max. iterations", s.max_iterations)
        .kv(
            "Convergence criterion",
            format!("{}  (every scaled residual)", sci(s.tolerance, 1)),
        );
    f.push_kv(
        "Initial field",
        if case.run.restart.is_some() {
            "restart file".to_string()
        } else if case.initial.potential_flow {
            "potential flow (D52)".to_string()
        } else {
            format!(
                "u = ({}, {}, {}) m/s, p = {} Pa",
                case.initial.velocity[0],
                case.initial.velocity[1],
                case.initial.velocity[2],
                case.initial.pressure
            )
        },
    );
    cadence(&mut f, case, "iterations");
    f
}

fn cadence(f: &mut Frame, case: &FlowCase, unit: &str) {
    if !case.output.probes.is_empty() {
        f.push_kv(
            "Monitoring points",
            case.output
                .probes
                .iter()
                .map(|p| {
                    format!(
                        "{} at ({}, {}, {}) m",
                        p.name, p.point[0], p.point[1], p.point[2]
                    )
                })
                .collect::<Vec<_>>()
                .join("\n"),
        );
    }
    f.push_kv("Log every", format!("{} {unit}", case.solver.log_every));
    f.push_kv(
        "History every",
        format!("{} {unit}", case.output.history_every),
    );
    f.push_kv(
        "Balance table every",
        if case.balance.every == 0 {
            "run end only".to_string()
        } else {
            format!("{} {unit}", case.balance.every)
        },
    );
    if case.run.reread_every > 0 {
        f.push_kv(
            "Case file re-read",
            format!("every {} {unit}", case.run.reread_every),
        );
    }
    f.push_kv("Control files", ".stop, .result_N, .read_toml_N");
    // the changes prepared in `[[schedule]]` (D128)
    for (when, list) in crate::schedule::rows(case) {
        f.push_kv(&format!("Scheduled: {when}"), list);
    }
}

/// What a run of a transient case does.
pub fn control_transient(case: &FlowCase, dt: f64, steps: usize) -> Frame {
    let t = case.time.unwrap_or_default();
    let s = &case.solver;
    let mut f = Frame::new("Run control")
        .kv("Run mode", "transient")
        .kv("Flow category", flow_category(case))
        .kv("Algorithm", "PISO")
        .kv("Time integration", t.scheme.label())
        .kv(
            "Time step",
            if t.adaptive {
                format!("{} s, adaptive (Courant number {})", sci(dt, 4), t.max_cfl)
            } else {
                format!("{} s", sci(dt, 4))
            },
        )
        .kv("Start time", format!("{} s", t.start_time))
        .kv(
            "End",
            match (t.end_time, t.steps) {
                (Some(e), _) => format!("t = {e} s  ({steps} steps)"),
                (None, Some(n)) => format!("{n} steps"),
                _ => "open".to_string(),
            },
        )
        .kv(
            "Correctors / outer iterations",
            format!("{} / {}", s.correctors, s.outer_iterations),
        );
    f.push_kv(
        "Initial field",
        if case.run.restart.is_some() {
            "restart file".to_string()
        } else if case.initial.potential_flow {
            "potential flow (D52)".to_string()
        } else {
            format!(
                "u = ({}, {}, {}) m/s, p = {} Pa",
                case.initial.velocity[0],
                case.initial.velocity[1],
                case.initial.velocity[2],
                case.initial.pressure
            )
        },
    );
    cadence(&mut f, case, "steps");
    f
}

/// The mesh: counts, size and quality.
pub fn mesh(q: &QualityReport, desc: &str) -> Frame {
    let mut f = Frame::new("Mesh information");
    f.push_kv("Source", desc);
    let shapes = q
        .shapes
        .iter()
        .filter(|(_, n)| *n > 0)
        .map(|(s, n)| format!("{} {s}", count(*n)))
        .collect::<Vec<_>>()
        .join(", ");
    f.push_kv("No. of cells", format!("{}  ({shapes})", count(q.cells)));
    f.push_kv(
        "No. of faces",
        format!(
            "{}  (interior {}, boundary {})",
            count(q.faces),
            count(q.interior_faces),
            count(q.boundary_faces)
        ),
    );
    f.push_kv("No. of nodes", count(q.points));
    let (lo, hi) = q.bounding_box;
    f.push_kv(
        "Dimensions [m]",
        format!(
            "{} x {} x {}",
            sci(hi[0] - lo[0], 3),
            sci(hi[1] - lo[1], 3),
            sci(hi[2] - lo[2], 3)
        ),
    );
    f.push_kv("Volume [m3]", sci(q.total_volume, 4));
    f.push_rule();
    f.push_text(format!(
        "{:<26}{:>11}{:>12}{:>12}",
        "Quality", "Min", "Mean", "Max"
    ));
    let row = |name: &str, s: &foehncraft_mesh::quality::Stat, d: usize| {
        format!(
            "{name:<26}{:>11.d$}{:>12.d$}{:>12.d$}",
            s.min,
            s.mean,
            s.max,
            d = d
        )
    };
    f.push_text(row("Non-orthogonality [deg]", &q.non_orthogonality, 1));
    f.push_text(row("Skewness [-]", &q.skewness, 3));
    f.push_text(row("Aspect ratio [-]", &q.aspect_ratio, 2));
    f.push_text(row("Volume ratio [-]", &q.volume_ratio, 2));
    if q.non_orthogonal_above_70 > 0 || q.negative_volumes > 0 || q.inverted_faces > 0 {
        f.push_rule();
        if q.non_orthogonal_above_70 > 0 {
            f.push_kv(
                "Faces above 70 deg",
                format!("{}  (non-orthogonality)", q.non_orthogonal_above_70),
            );
        }
        if q.negative_volumes > 0 {
            f.push_kv("Negative volumes", q.negative_volumes);
        }
        if q.inverted_faces > 0 {
            f.push_kv("Inverted faces", q.inverted_faces);
        }
    }
    f
}

/// The fluid.
pub fn fluid(case: &FlowCase) -> Frame {
    let mut f = Frame::new("Fluid properties").kv(
        "Material",
        case.fluid
            .material
            .clone()
            .unwrap_or_else(|| "constant (case file)".into()),
    );
    if let Ok(p) = case.properties() {
        f.push_kv("Density", format!("{}  [kg/m3]", sci(p.density(), 4)));
        f.push_kv(
            "Dynamic viscosity",
            format!("{}  [Pa s]", sci(p.reference_viscosity(), 4)),
        );
    }
    if let Some(re) = case.reynolds() {
        f.push_kv("Reynolds number", sci(re, 3));
    }
    if let Some(t) = case.fluid.temperature {
        f.push_kv("Temperature", format!("{t}  [K]"));
    }
    f.push_kv(
        "Convection",
        if case.flow.convection {
            "on"
        } else {
            "off (Stokes)"
        },
    );
    let g = &case.flow.body_force;
    match crate::expr::numbers(g) {
        Some(n) if n == [0.0; 3] => {}
        Some(n) => f.push_kv(
            "Body force [N/m3]",
            format!("({}, {}, {})", n[0], n[1], n[2]),
        ),
        None => f.push_kv("Body force [N/m3]", "expressions (see the case file)"),
    }
    f
}

/// The "solution flags": one column per equation.
pub fn solution(case: &FlowCase, transient: bool) -> Frame {
    let s = &case.solver;
    let m = &s.momentum;
    let p = &s.pressure;
    let w = 12;
    let row = |label: &str, u: &str, pr: &str| format!("{label:<24}{u:>w$}{u:>w$}{u:>w$}{pr:>w$}");
    let mut f = Frame::new("Solution settings");
    f.push_text(format!(
        "{:<24}{:>w$}{:>w$}{:>w$}{:>w$}",
        "Variable", "U", "V", "W", "P"
    ));
    f.push_rule();
    f.push_text(row(
        "Linear solver",
        "BiCGStab",
        &format!(
            "CG+{}",
            match spelled(&p.preconditioner).as_str() {
                "amg" => "AMG".to_string(),
                o => o.to_string(),
            }
        ),
    ));
    f.push_text(row(
        "Max. solver iterations",
        &m.max_iter.to_string(),
        &p.max_iter.to_string(),
    ));
    f.push_text(row(
        "Solver tolerance",
        &sci(s.momentum_tol(transient), 0),
        &sci(s.pressure_tol(transient), 0),
    ));
    f.push_text(row(
        "Relaxation factor",
        &format!("{:.2}", s.alpha_u(transient)),
        &format!("{:.2}", s.alpha_p(transient)),
    ));
    f.push_text(row("Differencing scheme", &short(&spelled(&m.scheme)), "-"));
    f.push_text(row("Gradient limiter", &short(&spelled(&m.limiter)), "-"));
    f.push_rule();
    f.push_kv("Gradients", "weighted least squares");
    f.push_kv(
        "Non-orthogonal correction",
        if s.nonorth_correction { "on" } else { "off" },
    );
    f.push_kv(
        "Multigrid rebuild",
        format!(
            "every {} {}",
            s.amg_rebuild,
            if transient { "steps" } else { "iterations" }
        ),
    );
    f
}

/// Short names for the columns of the solution settings.
fn short(name: &str) -> String {
    match name {
        "linear-upwind" => "LUD".into(),
        "upwind" => "UD".into(),
        "central" => "CD".into(),
        "venkatakrishnan" => "Venkat.".into(),
        "venkatakrishnan-wang" => "Venkat-W".into(),
        "barth-jespersen" => "BJ".into(),
        o => o.to_string(),
    }
}

/// What a boundary condition sets, from its serialised form.
pub fn describe(bc: &BoundaryCondition) -> (String, String) {
    let v = toml::Value::try_from(bc).ok();
    let Some(toml::Value::Table(t)) = v else {
        return (String::new(), String::new());
    };
    let kind = t
        .get("type")
        .and_then(toml::Value::as_str)
        .unwrap_or("")
        .to_string();
    let vals = t
        .iter()
        .filter(|(k, _)| *k != "type")
        .map(|(k, v)| format!("{k} = {}", inline(v)))
        .collect::<Vec<_>>()
        .join(", ");
    (kind, vals)
}

fn inline(v: &toml::Value) -> String {
    match v {
        toml::Value::String(s) => format!("\"{s}\""),
        toml::Value::Float(x) => format!("{x}"),
        toml::Value::Array(a) => {
            format!("[{}]", a.iter().map(inline).collect::<Vec<_>>().join(", "))
        }
        toml::Value::Table(t) => format!(
            "{{ {} }}",
            t.iter()
                .map(|(k, v)| format!("{k} = {}", inline(v)))
                .collect::<Vec<_>>()
                .join(", ")
        ),
        o => o.to_string(),
    }
}

/// The boundaries of the mesh with their conditions.
pub fn boundaries(case: &FlowCase, q: &QualityReport) -> Frame {
    let mut f = Frame::new("Boundary conditions");
    f.push_text(format!(
        "{:<16}{:<18}{:>9}{:>13}",
        "Boundary", "Type", "Faces", "Area [m2]"
    ));
    f.push_rule();
    for (name, faces, area) in &q.boundaries {
        let (kind, vals) = match case.boundary.get(name) {
            Some(b) => describe(b),
            None => ("(none)".into(), String::new()),
        };
        f.push_text(format!(
            "{name:<16}{kind:<18}{faces:>9}{:>13}",
            sci(*area, 3)
        ));
        if !vals.is_empty() {
            f.push_text(format!("    {vals}"));
        }
    }
    f
}

/// The `[I]` lines of the notes of reading the case.
pub fn notes(case: &FlowCase) -> Vec<String> {
    case.notes
        .iter()
        .map(|n| crate::console::tagged("note", n))
        .chain(
            case.pressure_frame
                .iter()
                .flat_map(|f| f.notes())
                .map(|n| crate::console::tagged("note", &n)),
        )
        .collect()
}

/// A line of the run time analysis with its share.
pub fn time_rows(f: &mut Frame, rows: &[(&str, f64)], total: f64) {
    for (name, s) in rows {
        f.push_kv(
            name,
            format!(
                "{:>10.2} s  {:>5.1} %",
                s,
                if total > 0.0 { 100.0 * s / total } else { 0.0 }
            ),
        );
    }
}

/// The end of the output list: what a file is about.
pub fn file_size(bytes: usize) -> String {
    size(bytes)
}

/// Values of the end-of-run summary that are used by both drivers.
pub fn quantity(v: f64, unit: &str) -> String {
    format!("{}  [{unit}]", sci5(v))
}

/// `clock` and `sci3` re-exported for the drivers.
pub fn elapsed(secs: f64) -> String {
    clock(secs)
}

/// A residual cell of a table.
pub fn cell(v: f64) -> String {
    sci3(v)
}

/// A balance table (D88-D90) in a frame.
pub fn balance_block(title: &str, b: &crate::balance::Balance, integrated: bool) -> Frame {
    let (head, detail) = match title.split_once(", ") {
        Some((h, d)) => (h, Some(d)),
        None => (title, None),
    };
    let text = crate::balance::table(head, b, integrated);
    let mut f = Frame::new(head).summary();
    if let Some(d) = detail {
        f.push_text(d.to_string());
        f.push_rule();
    }
    for l in text.lines().skip(1) {
        f.push_text(l.strip_prefix("  ").unwrap_or(l).to_string());
    }
    f
}

/// Mass flow, mean pressure and force per boundary, with the torque of
/// the walls about the axis of a turning zone (D184).
pub fn results(bs: &[crate::case::BoundarySummary]) -> Frame {
    let mut f = Frame::new("Results per boundary").summary();
    let torque = bs.iter().any(|b| b.torque.is_some());
    f.push_text(format!(
        "{:<16}{:>11}{:>11}{:>11}{:>11}{:>11}{}",
        "Boundary",
        "Mass flow",
        "Mean p",
        "Force x",
        "Force y",
        "Force z",
        if torque {
            format!("{:>11}", "Torque")
        } else {
            String::new()
        }
    ));
    f.push_text(format!(
        "{:<16}{:>11}{:>11}{:>11}{:>11}{:>11}{}",
        "",
        "[kg/s]",
        "[Pa]",
        "[N]",
        "[N]",
        "[N]",
        if torque {
            format!("{:>11}", "[N m]")
        } else {
            String::new()
        }
    ));
    f.push_rule();
    for b in bs {
        f.push_text(format!(
            "{:<16}{:>11}{:>11}{:>11}{:>11}{:>11}{}",
            b.name,
            sci3(b.mass_flow),
            sci3(b.mean_pressure),
            sci3(b.force[0]),
            sci3(b.force[1]),
            sci3(b.force[2]),
            if torque {
                format!("{:>11}", b.torque.map_or("-".to_string(), sci3))
            } else {
                String::new()
            }
        ));
    }
    f
}

/// Torque, shaft power and the performance of the machines (D184).
pub fn rotors(v: &crate::rotor::Values) -> Frame {
    let mut f = Frame::new("Rotors and machines").summary();
    for r in &v.rotors {
        f.push_text(format!("Rotor {}  ({})", r.name, r.walls.join(", ")));
        f.push_kv("  Speed", format!("{} rad/s", sci3(r.omega)));
        f.push_kv(
            "  Torque of the fluid on it",
            format!("{} N m", sci3(r.torque)),
        );
        f.push_kv("  Shaft power  (-Ω·T)", format!("{} W", sci3(r.power)));
        let a = &r.angular_momentum;
        f.push_kv(
            "  Angular momentum",
            format!(
                "rotor {}, walls {}, openings {}, out {} N m; imbalance {:.1e} of the largest",
                sci3(a.rotor),
                sci3(a.walls),
                sci3(a.openings),
                sci3(a.outflow),
                a.relative
            ),
        );
    }
    for m in &v.machines {
        f.push_rule();
        f.push_text(format!("Machine {}  (rotor {})", m.name, m.rotor));
        f.push_kv("  Volume flow", format!("{} m³/s", sci3(m.volume_flow)));
        f.push_kv(
            "  Total pressure rise",
            format!(
                "{} Pa  (total to static {} Pa)",
                sci3(m.total_pressure_rise),
                sci3(m.total_to_static_pressure_rise)
            ),
        );
        f.push_kv("  Head", format!("{} m", sci3(m.head)));
        f.push_kv("  Shaft power", format!("{} W", sci3(m.shaft_power)));
        f.push_kv(
            "  Hydraulic power",
            format!("{} W", sci3(m.hydraulic_power)),
        );
        f.push_kv(
            "  Efficiency",
            format!(
                "{:.4}  (total to static {:.4})",
                m.efficiency, m.total_to_static_efficiency
            ),
        );
        f.push_kv(
            "  Euler work",
            format!("{} W  (head {} m)", sci3(m.euler_work), sci3(m.euler_head)),
        );
    }
    f
}

/// The run time analysis of one solve.
pub fn run_time(t: &crate::simplec::Timings, cells: usize, iterations: usize, unit: &str) -> Frame {
    let total = t.total();
    let mut f = Frame::new("Run time analysis").summary();
    let mut rows = vec![
        ("Gradients and limiters", t.gradients),
        ("Momentum assembly", t.momentum_assembly),
        ("Momentum solver", t.momentum_solve),
        ("Pressure assembly", t.pressure_assembly),
        ("Pressure solver", t.pressure_solve),
        ("Multigrid setup", t.amg_setup),
    ];
    // the phases a case without a turbulence model or scalars does not have
    if t.turbulence > 0.0 {
        rows.push(("Turbulence model", t.turbulence));
    }
    if t.scalars > 0.0 {
        rows.push(("Passive scalars", t.scalars));
    }
    time_rows(&mut f, &rows, total);
    f.push_rule();
    f.push_kv("Total", format!("{total:>10.2} s"));
    f.push_kv(
        "Per cell and iteration",
        format!(
            "{:.0} ns  ({iterations} {unit}, {} cells)",
            1e9 * total / (cells as f64 * iterations.max(1) as f64),
            count(cells)
        ),
    );
    f
}
