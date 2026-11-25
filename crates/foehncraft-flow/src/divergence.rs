//! Divergence report (D63).
//!
//! A run checks the residuals of every outer iteration (steady) or
//! corrector (transient) for values that are not finite; that check was
//! there before and costs nothing. Before each iteration or step the
//! solver keeps a device copy of u, p and ṁ ([`Flow::keep_valid`], four
//! fields against the dozens of passes of an iteration). When the check
//! fails, the run
//!
//! 1. takes the failed fields to the host, puts the last valid state back
//!    ([`Flow::rewind`], [`crate::transient::Transient::rewind`]) and writes
//!    it as `fields_last_valid.vtkhdf` and `restart.h5` (a restart from it
//!    continues the run with changed settings; transient, not
//!    bit-identically);
//! 2. writes `divergence.vtkhdf`, the failed fields on the failing cells
//!    and three layers of neighbours, with the field `divergence_layer`
//!    (0 failing, 1–3 neighbour layers) and the last valid fields beside
//!    the failed ones;
//! 3. writes the report `divergence.md` and `divergence.json` (format
//!    `foehncraft.divergence`): the failing cells with their global ID and
//!    rank, the variable that failed first, the cell nearest to failing in
//!    the last valid state, and the likely causes ranked by their evidence
//!    with what to change, as `--set` arguments where one fits (codes
//!    `DIV001` … `DIV006`);
//! 4. ends with exit code 1 and the code `RUN014`.
//!
//! The bounded quantities of the later models (turbulence, species,
//! volume fractions) get their range warnings at the progress lines with
//! those models; the laminar flow of stage 5 has none.

use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use crate::simplec::{Flow, Residuals, Solution, State};
use foehncraft_mesh::geometry::{norm, sub};
use foehncraft_mesh::{Geometry, Mesh};
use serde::Serialize;
use std::fmt::Write as _;
use std::path::Path;

/// Format name of `divergence.json` (D80).
pub const FORMAT: &str = "foehncraft.divergence";
/// Format version of `divergence.json`.
pub const FORMAT_VERSION: &str = "1.0";
/// Neighbour layers in the snapshot.
pub const LAYERS: i32 = 3;
/// Failing cells listed in the report (the snapshot holds all).
pub const LISTED: usize = 20;

/// Where the run stood when it failed.
#[derive(Debug, Clone)]
pub struct Context<'a> {
    /// Case name.
    pub name: &'a str,
    /// The case in force.
    pub case: &'a FlowCase,
    /// Outer iteration (steady) or time step (transient) that failed.
    pub at: usize,
    /// Iterations or steps of this run before the failure.
    pub done: usize,
    /// Time of the last valid state (transient).
    pub time: Option<f64>,
    /// The time step that failed (transient).
    pub dt: Option<f64>,
    /// Residuals of the iteration or corrector that failed.
    pub residuals: Residuals,
    /// Residuals of the last valid iteration or step.
    pub last_residuals: Option<Residuals>,
    /// Largest Courant number of the first step of this run (transient).
    pub first_courant: Option<f64>,
}

/// A failing cell.
#[derive(Debug, Clone, Serialize)]
pub struct FailedCell {
    /// Global cell ID (`cell_id` of the mesh files).
    pub cell: u64,
    /// Rank that holds it.
    pub rank: usize,
    /// Centroid (m).
    pub centroid: [f64; 3],
    /// Velocity in the failed state (null: not finite) (m/s).
    pub velocity: [f64; 3],
    /// Pressure in the failed state (Pa).
    pub pressure: f64,
    /// Speed in the last valid state (m/s).
    pub last_valid_speed: f64,
}

/// The cell nearest to failing in the last valid state.
#[derive(Debug, Clone, Serialize)]
pub struct NearestCell {
    /// Global cell ID.
    pub cell: u64,
    /// Rank.
    pub rank: usize,
    /// Centroid (m).
    pub centroid: [f64; 3],
    /// What made it the nearest: `courant` (transient) or `speed`.
    pub measure: &'static str,
    /// Its value: Courant number (-) or speed (m/s).
    pub value: f64,
    /// Distance to the nearest failing cell (m).
    pub distance_to_failure: f64,
}

/// A likely cause.
#[derive(Debug, Clone, Serialize)]
pub struct Cause {
    /// Code (`DIV001` …).
    pub code: &'static str,
    /// Short name.
    pub cause: &'static str,
    /// Weight of the evidence, 0 … 1; the causes are sorted by it.
    pub score: f64,
    /// What points to it.
    pub evidence: String,
    /// What to change.
    pub advice: String,
    /// The change as `--set` arguments, where one fits.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub fix: Vec<String>,
}

/// The report (`divergence.json`).
#[derive(Debug, Clone, Serialize)]
pub struct Report {
    /// [`FORMAT`].
    pub format: &'static str,
    /// [`FORMAT_VERSION`].
    pub format_version: &'static str,
    /// Foehncraft version.
    pub foehncraft_version: &'static str,
    /// Case name.
    pub case: String,
    /// `steady` or `transient`.
    pub kind: &'static str,
    /// Outer iteration or time step that failed.
    pub at: usize,
    /// Iterations or steps of this run before it.
    pub done: usize,
    /// Time of the last valid state (s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub time: Option<f64>,
    /// The failed time step (s).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub time_step: Option<f64>,
    /// `velocity` or `pressure`: the variable that failed first (the
    /// momentum predictor runs before the pressure correction).
    pub first_variable: &'static str,
    /// Residuals of the failed iteration or corrector.
    pub residuals: Residuals,
    /// Number of failing cells.
    pub failing_cells: usize,
    /// What makes a cell failing: a velocity or pressure that is not
    /// finite, a mass flux that is not finite at one of its faces, or (when
    /// nothing is non-finite yet) a speed above a tenth of the largest.
    pub failing_criterion: &'static str,
    /// The first failing cells by global ID (all of them are in the
    /// snapshot).
    pub cells: Vec<FailedCell>,
    /// The cell nearest to failing in the last valid state.
    pub nearest: Option<NearestCell>,
    /// Likely causes, most likely first.
    pub causes: Vec<Cause>,
    /// Files written.
    pub files: Vec<String>,
}

impl Report {
    /// The diagnostic of the run: `RUN014`, with the fix of the most
    /// likely cause.
    pub fn diagnostic(&self) -> Diagnostic {
        let mut d = Diagnostic::new(
            diag::RUN_DIVERGED,
            format!(
                "the solution diverged in {} {} ({} failing cells, {} first); report in \
                 divergence.md{}",
                if self.kind == "steady" {
                    "outer iteration"
                } else {
                    "time step"
                },
                self.at,
                self.failing_cells,
                self.first_variable,
                self.causes.first().map_or(String::new(), |c| format!(
                    "; most likely cause: {} ({})",
                    c.cause, c.code
                ))
            ),
        );
        if let Some(c) = self.causes.iter().find(|c| !c.fix.is_empty()) {
            d = d.fix(c.fix.iter().cloned());
        }
        d
    }

    /// The text report.
    pub fn markdown(&self) -> String {
        let mut s = format!(
            "# Divergence of `{}`\n\nThe {} run diverged in {} {}",
            self.case,
            self.kind,
            if self.kind == "steady" {
                "outer iteration"
            } else {
                "time step"
            },
            self.at
        );
        if let (Some(t), Some(dt)) = (self.time, self.time_step) {
            let _ = write!(s, " (from t = {t:.6e} s with Δt = {dt:.3e} s)");
        }
        let _ = writeln!(
            s,
            ", after {} {} of this run. The {} failed first ({}).\n",
            self.done,
            if self.kind == "steady" {
                "iterations"
            } else {
                "steps"
            },
            self.first_variable,
            if self.residuals == Residuals::default() {
                "a linear solver broke down before the residuals of the iteration were complete"
                    .to_string()
            } else {
                format!(
                    "residuals momentum {:.3e} {:.3e} {:.3e}, continuity {:.3e}",
                    self.residuals.momentum[0],
                    self.residuals.momentum[1],
                    self.residuals.momentum[2],
                    self.residuals.continuity
                )
            }
        );
        let _ = writeln!(s, "## Likely causes\n");
        if self.causes.is_empty() {
            let _ = writeln!(
                s,
                "None of the checks found evidence; look at the failing cells in \
                 `divergence.vtkhdf`.\n"
            );
        }
        for (i, c) in self.causes.iter().enumerate() {
            let _ = writeln!(
                s,
                "{}. **{}** (`{}`, evidence {:.2}): {}. {}{}",
                i + 1,
                c.cause,
                c.code,
                c.score,
                c.evidence,
                c.advice,
                if c.fix.is_empty() {
                    String::new()
                } else {
                    format!(
                        " Try `{}`.",
                        c.fix
                            .iter()
                            .map(|f| format!("--set {f}"))
                            .collect::<Vec<_>>()
                            .join(" ")
                    )
                }
            );
        }
        let _ = writeln!(
            s,
            "\n## Failing cells\n\n{} cells{}{}.\n\n\
             | Cell (global ID) | Rank | Centroid (m) | Velocity (m/s) | Pressure (Pa) | Last valid speed (m/s) |\n|---|---|---|---|---|---|",
            self.failing_cells,
            self.failing_criterion,
            if self.failing_cells > self.cells.len() {
                format!("; the first {} by ID", self.cells.len())
            } else {
                String::new()
            }
        );
        for c in &self.cells {
            let _ = writeln!(
                s,
                "| {} | {} | {:.4e}, {:.4e}, {:.4e} | {:.3e}, {:.3e}, {:.3e} | {:.3e} | {:.4e} |",
                c.cell,
                c.rank,
                c.centroid[0],
                c.centroid[1],
                c.centroid[2],
                c.velocity[0],
                c.velocity[1],
                c.velocity[2],
                c.pressure,
                c.last_valid_speed
            );
        }
        if let Some(n) = &self.nearest {
            let _ = writeln!(
                s,
                "\nNearest to failing in the last valid state: cell {} (rank {}) at {:.4e}, {:.4e}, \
                 {:.4e} m, {} {:.4e}, {:.3e} m from the nearest failing cell.",
                n.cell,
                n.rank,
                n.centroid[0],
                n.centroid[1],
                n.centroid[2],
                if n.measure == "courant" {
                    "Courant number"
                } else {
                    "speed (m/s)"
                },
                n.value,
                n.distance_to_failure
            );
        }
        let _ = writeln!(s, "\n## Files\n");
        for f in &self.files {
            let _ = writeln!(s, "- `{f}`");
        }
        s
    }
}

/// Courant number of each cell for a step `dt`: Δt Σ|ṁ_f| / (2 ρ V).
fn courant(mesh: &Mesh, geo: &Geometry, mass_flux: &[f64], rho: f64, dt: f64) -> Vec<f64> {
    let nc = mesh.n_cells;
    let ni = mesh.periodic_start();
    let mut s = vec![0.0; nc];
    for (f, &m) in mass_flux.iter().enumerate().take(mesh.n_faces()) {
        let a = m.abs();
        s[mesh.owner[f] as usize] += a;
        if f < ni {
            let n = mesh.neighbour[f] as usize;
            if n < nc {
                s[n] += a;
            }
        }
    }
    (0..nc)
        .map(|c| dt * s[c] / (2.0 * rho * geo.cell_volume[c]))
        .collect()
}

/// Layers from the failing cells: 0 failing, 1 … `LAYERS` neighbours, −1
/// beyond.
fn layers(mesh: &Mesh, failing: &[usize]) -> Vec<i32> {
    let nc = mesh.n_cells;
    let mut layer = vec![-1i32; nc];
    let mut front: Vec<usize> = failing.to_vec();
    for &c in failing {
        layer[c] = 0;
    }
    for l in 1..=LAYERS {
        let mut next = Vec::new();
        for &c in &front {
            for (f, _) in mesh.cell_face_iter(c) {
                if let Some(d) = mesh.across(f, c)
                    && d < nc
                    && layer[d] < 0
                {
                    layer[d] = l;
                    next.push(d);
                }
            }
        }
        next.sort_unstable();
        front = next;
    }
    layer
}

/// Whether the error `e` of an iteration or step is a divergence: the
/// check of the residuals (`RUN014`), a breakdown of a linear solver, or
/// another error without a code of its own while the fields hold values
/// that are not finite. Downloads the fields only for the
/// second kind.
pub(crate) fn is_divergence(e: &anyhow::Error, fl: &mut Flow) -> bool {
    match diag::of_error(e).code {
        diag::RUN_DIVERGED => true,
        diag::GENERAL if format!("{e:#}").contains("breakdown") => true,
        diag::GENERAL => fl.state().is_ok_and(|s| {
            (0..s.p.len()).any(|c| !finite3(s.u.get(c)) || !s.p[c].is_finite())
                || s.mass_flux.iter().any(|m| !m.is_finite())
        }),
        _ => false,
    }
}

fn finite3(v: [f64; 3]) -> bool {
    v.iter().all(|x| x.is_finite())
}

/// Analyses the failure: `failed` holds the fields after the failed
/// iteration, `valid` the last valid state (with its gradients).
pub fn analyse(
    ctx: &Context,
    mesh: &Mesh,
    geo: &Geometry,
    mom_kind: &[u32],
    failed: &State,
    valid: &Solution,
) -> (Report, Vec<i32>) {
    let case = ctx.case;
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let transient = case.time.is_some();
    let mut failing: Vec<usize> = (0..nc)
        .filter(|&c| !finite3(failed.u.get(c)) || !failed.p[c].is_finite())
        .collect();
    let mut criterion = " with a velocity or pressure that is not finite";
    if failing.is_empty() {
        // the fields are finite but a flux is not: the cells of those faces
        for f in 0..mesh.n_faces() {
            if !failed.mass_flux[f].is_finite() {
                failing.push(mesh.owner[f] as usize);
                if f < ni {
                    failing.push(mesh.neighbour[f] as usize);
                }
            }
        }
        failing.retain(|&c| c < nc);
        failing.sort_unstable();
        failing.dedup();
        criterion = " at a face whose mass flux is not finite";
    }
    let mut blew_up = false;
    if failing.is_empty() {
        // nothing is non-finite (a linear solver broke down first): the
        // cells where the speed exceeds a tenth of its largest value
        let speed: Vec<f64> = (0..nc).map(|c| norm(failed.u.get(c))).collect();
        let top = speed.iter().copied().fold(0.0f64, f64::max);
        failing = (0..nc).filter(|&c| speed[c] >= 0.1 * top).collect();
        blew_up = top > 0.0;
        criterion = " where the speed exceeds a tenth of its largest value (no value is \
                     non-finite yet)";
    }
    let layer = layers(mesh, &failing);
    let near: Vec<usize> = (0..nc).filter(|&c| (0..=1).contains(&layer[c])).collect();
    // the momentum predictor runs first: non-finite momentum residuals, or
    // (a linear solver broke down before the residuals were complete)
    // non-finite velocities, point to the velocity
    let first_variable = if ctx.residuals.momentum.iter().any(|r| !r.is_finite())
        || (ctx.residuals.finite() && (blew_up || (0..nc).any(|c| !finite3(failed.u.get(c)))))
    {
        "velocity"
    } else {
        "pressure"
    };
    let mut by_id: Vec<usize> = failing.clone();
    by_id.sort_by_key(|&c| mesh.cell_id[c]);
    let cells: Vec<FailedCell> = by_id
        .iter()
        .take(LISTED)
        .map(|&c| FailedCell {
            cell: mesh.cell_id[c],
            rank: 0,
            centroid: geo.cell_centroid[c],
            velocity: failed.u.get(c),
            pressure: failed.p[c],
            last_valid_speed: norm(valid.u.get(c)),
        })
        .collect();
    let rho = case.rho();
    let dt = ctx.dt.unwrap_or(0.0);
    let co = if transient {
        Some(courant(mesh, geo, &valid.mass_flux, rho, dt))
    } else {
        None
    };
    let measure: Vec<f64> = match &co {
        Some(v) => v.clone(),
        None => (0..nc).map(|c| norm(valid.u.get(c))).collect(),
    };
    let nearest = (0..nc)
        .max_by(|&a, &b| measure[a].total_cmp(&measure[b]))
        .map(|c| {
            let x = geo.cell_centroid[c];
            let d = failing
                .iter()
                .map(|&f| norm(sub(geo.cell_centroid[f], x)))
                .fold(f64::INFINITY, f64::min);
            NearestCell {
                cell: mesh.cell_id[c],
                rank: 0,
                centroid: x,
                measure: if transient { "courant" } else { "speed" },
                value: measure[c],
                distance_to_failure: d,
            }
        });

    let s = &case.solver;
    let mut causes = Vec::new();
    // DIV001: the time step. The Courant number of the first step of
    // this run measures the setting; the one of the last valid state shows
    // how far the flow had run away before a value became non-finite.
    if let (Some(co), Some(t)) = (&co, &case.time) {
        let local = near.iter().map(|&c| co[c]).fold(0.0f64, f64::max);
        let global = co.iter().copied().fold(0.0f64, f64::max);
        let first = ctx.first_courant.unwrap_or(global);
        let setting = first.max(if ctx.done < 3 { global } else { 0.0 });
        let score = if setting > t.max_cfl {
            (0.5 + 0.25 * (setting / t.max_cfl).log2()).min(1.0)
        } else if global > 10.0 * t.max_cfl {
            0.5
        } else {
            (0.3 * global / t.max_cfl).min(0.3)
        };
        let (advice, fix) = if t.adaptive {
            let m = (0.5 * t.max_cfl).max(0.05);
            (
                format!(
                    "Halve the largest Courant number of the adaptive step (now {}).",
                    t.max_cfl
                ),
                vec![format!("time.max_cfl={}", round2(m))],
            )
        } else {
            let f = if setting > t.max_cfl {
                t.max_cfl / setting
            } else {
                0.5
            };
            (
                format!(
                    "Take a smaller time step (now {dt:.3e} s), or let the step adapt to \
                     [time] max_cfl = {} (adaptive = true).",
                    t.max_cfl
                ),
                vec![format!("time.time_step={}", sig3(dt * f.min(0.5)))],
            )
        };
        causes.push(Cause {
            code: diag::DIV_TIME_STEP,
            cause: "time step",
            score,
            evidence: format!(
                "Courant number {} in the first step of the run, {} in the last valid state \
                 with the failed step ({} at the failing cells)",
                num(first),
                num(global),
                num(local)
            ),
            advice,
            fix,
        });
    }
    // DIV002: relaxation; SIMPLE pairs the factors as α_p ≈ 1 − α_u
    let (au, ap) = (s.alpha_u(transient), s.alpha_p(transient));
    let simple = s.algorithm(transient) == crate::simplec::Algorithm::Simple;
    if !transient || s.outer_iterations > 1 {
        let score = if transient {
            0.2
        } else {
            let mut x = ((au - 0.5) / 0.4).clamp(0.0, 1.0) * 0.6;
            if simple {
                x = x.max(((au + ap - 1.0) / 0.5).clamp(0.0, 1.0) * 0.8);
            }
            x.max(if ctx.done < 50 { 0.3 } else { 0.1 })
        };
        let nu = round2((0.7 * au).max(0.1));
        let mut fix = vec![format!("solver.momentum.relaxation={nu}")];
        let mut advice = "Relax the velocity more strongly.".to_string();
        if simple {
            let np = round2((0.7 * ap).min(1.0 - nu).max(0.05));
            if np < ap {
                fix.push(format!("solver.pressure.relaxation={np}"));
                advice = "Relax more strongly; SIMPLE wants α_p ≈ 1 − α_u.".into();
            }
        }
        causes.push(Cause {
            code: diag::DIV_RELAXATION,
            cause: "under-relaxation",
            score,
            evidence: format!(
                "relaxation factors velocity {au}, pressure {ap} ({})",
                match s.algorithm(transient) {
                    crate::simplec::Algorithm::Simplec => "SIMPLEC",
                    crate::simplec::Algorithm::Simple => "SIMPLE",
                    crate::simplec::Algorithm::Piso => "PISO",
                    crate::simplec::Algorithm::Coupled => "coupled",
                }
            ),
            advice,
            fix,
        });
    }
    // DIV003: cell quality
    {
        let (_, q) = foehncraft_mesh::quality::QualityReport::compute(mesh, geo);
        let theta = near
            .iter()
            .map(|&c| q.non_orthogonality[c])
            .fold(0.0f64, f64::max);
        let ratio = near
            .iter()
            .map(|&c| q.volume_ratio[c])
            .fold(0.0f64, f64::max);
        let skew = near.iter().map(|&c| q.skewness[c]).fold(0.0f64, f64::max);
        let skewed = skew > s.skewness_limit;
        // skewed faces anywhere in the mesh with the remediation off: the
        // least-squares pressure gradient on them can drive the whole
        // iteration away, also where the values run away first
        let off = s.remediation == crate::simplec::Remediation::Off;
        let poor: Vec<f64> = q
            .skewness
            .iter()
            .copied()
            .filter(|&k| k > s.skewness_limit)
            .collect();
        let poor_max = poor.iter().copied().fold(0.0f64, f64::max);
        let untreated = off && !poor.is_empty();
        let score = ((theta - 50.0) / 30.0)
            .clamp(0.0, 1.0)
            .max(((ratio - 5.0) / 15.0).clamp(0.0, 1.0))
            .max(if skewed { 0.6 } else { 0.0 })
            .max(if untreated { 0.5 } else { 0.0 });
        if score > 0.0 {
            let mut fix = Vec::new();
            if !s.nonorth_correction && theta > 40.0 {
                fix.push("solver.nonorth_correction=true".to_string());
            }
            if (skewed || untreated) && off {
                fix.push("solver.remediation=\"auto\"".to_string());
            }
            let mut evidence = format!(
                "at the failing cells and their neighbours: non-orthogonality up to \
                 {theta:.1}°, skewness up to {skew:.2}, volume ratio up to {ratio:.2}"
            );
            if untreated {
                evidence += &format!(
                    "; in the mesh {} cells with a face skewness above {} (largest {poor_max:.2}) \
                     and solver.remediation = \"off\"",
                    poor.len(),
                    s.skewness_limit
                );
            }
            causes.push(Cause {
                code: diag::DIV_CELL_QUALITY,
                cause: "cell quality",
                score,
                evidence,
                advice: if skewed || untreated {
                    "Faces with a skewness above `solver.skewness_limit` (D50): optimise the \
                     mesh (Gmsh: `Mesh.OptimizeNetgen = 1`), or let the cell-quality \
                     remediation treat them (`solver.remediation = \"auto\"`)."
                        .into()
                } else {
                    "Improve the mesh there (see `foehn mesh --quality`), or relax more \
                     strongly."
                        .into()
                },
                fix,
            });
        }
    }
    // DIV004: a boundary next to the failing cells
    {
        let mut hits: Vec<(String, String, bool)> = Vec::new();
        for p in &mesh.boundaries {
            let mut touch = false;
            let mut backflow = false;
            for f in p.range() {
                let c = mesh.owner[f] as usize;
                if c < nc && (0..=1).contains(&layer[c]) {
                    touch = true;
                    if mom_kind[f - ni] == foehncraft_kernels::flow::MOM_OUTFLOW
                        && valid.mass_flux[f] < 0.0
                    {
                        backflow = true;
                    }
                }
            }
            if touch {
                let kind = case
                    .boundary
                    .get(&p.name)
                    .map_or("?".to_string(), |b| b.type_name().to_string());
                hits.push((p.name.clone(), kind, backflow));
            }
        }
        if !hits.is_empty() {
            let back = hits.iter().any(|h| h.2);
            let names: Vec<String> = hits
                .iter()
                .map(|(n, k, b)| format!("{n} ({k}{})", if *b { ", backflow" } else { "" }))
                .collect();
            causes.push(Cause {
                code: diag::DIV_BOUNDARY,
                cause: "boundary condition",
                score: if back { 0.8 } else { 0.4 },
                evidence: format!("failing cells at the boundaries {}", names.join(", ")),
                advice: if back {
                    "Flow enters through an outflow boundary: move it further downstream or \
                     prescribe the pressure there; see the boundary-condition check (`foehn flow \
                     --check`)."
                        .into()
                } else {
                    "Check the condition of those boundaries against the guidance of `foehn flow \
                     --check`."
                        .into()
                },
                fix: Vec::new(),
            });
        }
    }
    // DIV005: the initial state
    if ctx.done < 10 {
        let score = if ctx.done < 3 { 0.7 } else { 0.4 };
        let fix = match &case.time {
            Some(t) if !t.adaptive => vec![format!("time.time_step={}", sig3(0.25 * dt))],
            Some(_) => vec!["time.max_growth=1.1".to_string()],
            None => vec![format!(
                "solver.momentum.relaxation={}",
                round2((0.5 * au).max(0.1))
            )],
        };
        causes.push(Cause {
            code: diag::DIV_INITIAL_STATE,
            cause: "initial state",
            score,
            evidence: format!(
                "the run failed after {} {}",
                ctx.done,
                if transient { "steps" } else { "iterations" }
            ),
            advice: "Start more gently: a smaller first step or stronger relaxation, or an \
                     initial velocity closer to the flow ([initial] velocity)."
                .into(),
            fix,
        });
    }
    // DIV006: a linear solve at its limit; in the last valid iteration it
    // weighs more than in the failed one, where a solver at its limit is
    // as much a consequence as a cause
    {
        let mut hit = Vec::new();
        let mut fix: Vec<String> = Vec::new();
        let mut score = 0.0f64;
        for (r, w, when) in [
            (Some(&ctx.residuals), 0.3, "failed"),
            (ctx.last_residuals.as_ref(), 0.6, "last valid"),
        ] {
            let Some(r) = r else { continue };
            if r.momentum_iterations
                .iter()
                .any(|&n| n >= s.momentum.max_iter)
            {
                hit.push(format!(
                    "momentum BiCGStab at its limit of {} iterations ({when} iteration)",
                    s.momentum.max_iter
                ));
                let f = format!("solver.momentum.max_iter={}", 2 * s.momentum.max_iter);
                if !fix.contains(&f) {
                    fix.push(f);
                }
                score = score.max(w);
            }
            if r.pressure_iterations >= s.pressure.max_iter {
                hit.push(format!(
                    "pressure CG at its limit of {} iterations ({when} iteration)",
                    s.pressure.max_iter
                ));
                let f = format!("solver.pressure.max_iter={}", 2 * s.pressure.max_iter);
                if !fix.contains(&f) {
                    fix.push(f);
                }
                score = score.max(w);
            }
        }
        if !hit.is_empty() {
            causes.push(Cause {
                code: diag::DIV_LINEAR_SOLVER,
                cause: "linear solve",
                score,
                evidence: hit.join("; "),
                advice: "Let the linear solver iterate further.".into(),
                fix,
            });
        }
    }
    causes.sort_by(|a, b| b.score.total_cmp(&a.score));
    (
        Report {
            format: FORMAT,
            format_version: FORMAT_VERSION,
            foehncraft_version: foehncraft_core::VERSION,
            case: ctx.name.into(),
            kind: if transient { "transient" } else { "steady" },
            at: ctx.at,
            done: ctx.done,
            time: ctx.time,
            time_step: ctx.dt,
            first_variable,
            residuals: ctx.residuals.clone(),
            failing_cells: failing.len(),
            failing_criterion: criterion,
            cells,
            nearest,
            causes,
            files: Vec::new(),
        },
        layer,
    )
}

fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

/// A number for the text: fixed for moderate values, else exponential.
fn num(x: f64) -> String {
    if x.abs() >= 1e-3 && x.abs() < 1e4 {
        format!("{x:.3}")
    } else {
        format!("{x:.3e}")
    }
}

fn sig3(x: f64) -> String {
    format!("{x:.2e}")
}

/// The snapshot `divergence.vtkhdf`: the failing cells and their
/// neighbour layers as a mesh of their own (or the whole mesh when they
/// fill it), with the failed and the last valid fields.
#[cfg(feature = "hdf5")]
fn write_snapshot(
    file: &Path,
    name: &str,
    mesh: &Mesh,
    failed: &State,
    valid: &Solution,
    layer: &[i32],
) -> anyhow::Result<()> {
    use foehncraft_mesh::vtkhdf::{Attr, CellField, write_with};
    let nc = mesh.n_cells;
    let inside: Vec<u32> = layer.iter().map(|&l| u32::from(l < 0)).collect();
    let sub = if inside.contains(&1) {
        foehncraft_mesh::decompose::extract(mesh, &inside, 0, 1).ok()
    } else {
        None
    };
    let (m, pick): (&Mesh, Vec<usize>) = match &sub {
        Some(l) => {
            let mut pos = std::collections::HashMap::with_capacity(nc);
            for (c, id) in mesh.cell_id[..nc].iter().enumerate() {
                pos.insert(*id, c);
            }
            let pick = l.mesh.cell_id[..l.mesh.n_cells]
                .iter()
                .map(|id| pos[id])
                .collect();
            (&l.mesh, pick)
        }
        None => (mesh, (0..nc).collect()),
    };
    let take = |v: &dyn Fn(usize) -> f64| -> Vec<f64> { pick.iter().map(|&c| v(c)).collect() };
    let cols: Vec<(&str, Vec<f64>)> = vec![
        ("divergence_layer", take(&|c| f64::from(layer[c]))),
        ("velocity_x", take(&|c| failed.u.x[c])),
        ("velocity_y", take(&|c| failed.u.y[c])),
        ("velocity_z", take(&|c| failed.u.z[c])),
        ("pressure", take(&|c| failed.p[c])),
        ("last_valid_velocity_x", take(&|c| valid.u.x[c])),
        ("last_valid_velocity_y", take(&|c| valid.u.y[c])),
        ("last_valid_velocity_z", take(&|c| valid.u.z[c])),
        ("last_valid_pressure", take(&|c| valid.p[c])),
        ("cell_id", take(&|c| mesh.cell_id[c] as f64)),
    ];
    let fields: Vec<CellField> = cols
        .iter()
        .map(|(n, v)| CellField { name: n, values: v })
        .collect();
    write_with(
        file,
        m,
        &fields,
        &[
            (
                "foehncraft_format",
                Attr::Str("foehncraft-divergence-fields"),
            ),
            ("foehncraft_format_version", Attr::U64(1)),
            ("foehncraft_version", Attr::Str(foehncraft_core::VERSION)),
            ("case", Attr::Str(name)),
            (
                "units",
                Attr::Str(
                    "coordinates m; velocity m/s; pressure Pa; divergence_layer -; cell_id -",
                ),
            ),
        ],
    )
}

#[cfg(not(feature = "hdf5"))]
fn write_snapshot(
    _: &Path,
    _: &str,
    _: &Mesh,
    _: &State,
    _: &Solution,
    _: &[i32],
) -> anyhow::Result<()> {
    Ok(())
}

/// Writes the report and the snapshot into `out` and returns the report;
/// `files` are the files the caller already wrote (last valid state,
/// restart).
pub fn write(
    out: &Path,
    mut report: Report,
    mesh: &Mesh,
    failed: &State,
    valid: &Solution,
    layer: &[i32],
    files: &[&str],
) -> anyhow::Result<Report> {
    report.files = files.iter().map(|s| s.to_string()).collect();
    if cfg!(feature = "hdf5") {
        write_snapshot(
            &out.join("divergence.vtkhdf"),
            &report.case,
            mesh,
            failed,
            valid,
            layer,
        )?;
        report.files.push("divergence.vtkhdf".into());
    }
    report.files.push("divergence.md".into());
    report.files.push("divergence.json".into());
    std::fs::write(out.join("divergence.md"), report.markdown())?;
    std::fs::write(
        out.join("divergence.json"),
        serde_json::to_string_pretty(&report)? + "\n",
    )?;
    Ok(report)
}

/// The last steps of D63, for both drivers, after they took the failed
/// fields to the host, put the last valid state back into `fl`
/// ([`Flow::rewind`], [`crate::transient::Transient::rewind`]) and wrote
/// its fields and restart file (`prepared`: the failed fields and the
/// names of those files): the analysis, the snapshot and the report.
/// Returns the error that ends the run (code `RUN014`, with the fix of the
/// most likely cause).
#[allow(clippy::too_many_arguments)]
pub(crate) fn finish(
    fl: &mut Flow,
    ctx: &Context,
    mesh: &Mesh,
    geo: &Geometry,
    out: &Path,
    prepared: anyhow::Result<(State, Vec<&'static str>)>,
    style: crate::console::Style,
    line: &mut dyn FnMut(&str),
) -> anyhow::Error {
    let valid = match &prepared {
        Ok(_) => fl.state().ok().map(|state| Solution {
            state,
            disc: &fl.disc,
            phys: fl.phys.clone(),
        }),
        Err(_) => None,
    };
    finish_from(
        valid.as_ref(),
        &fl.bc.mom_kind,
        ctx,
        mesh,
        geo,
        out,
        prepared,
        style,
        line,
    )
}

/// [`finish`] with the last valid state `valid` (on ranks gathered into the
/// global order, MPI step W7; `None` when it could not be read) and the
/// momentum boundary kinds of `mesh`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn finish_from(
    valid: Option<&Solution<'_>>,
    mom_kind: &[u32],
    ctx: &Context,
    mesh: &Mesh,
    geo: &Geometry,
    out: &Path,
    prepared: anyhow::Result<(State, Vec<&'static str>)>,
    style: crate::console::Style,
    line: &mut dyn FnMut(&str),
) -> anyhow::Error {
    let r = (|| -> anyhow::Result<Report> {
        let (failed, files) = prepared?;
        let valid =
            valid.ok_or_else(|| anyhow::anyhow!("the last valid state could not be read"))?;
        let (rep, layer) = analyse(ctx, mesh, geo, mom_kind, &failed, valid);
        write(out, rep, mesh, &failed, valid, &layer, &files)
    })();
    match r {
        Ok(rep) => {
            let unit = if rep.kind == "steady" {
                "outer iteration"
            } else {
                "time step"
            };
            let mut f = crate::console::Frame::new("Divergence")
                .summary()
                .kv(
                    "Diverged in",
                    format!("{unit} {}, {} first", rep.at, rep.first_variable),
                )
                .kv("Failing cells", rep.failing_cells)
                .kv("Report", "divergence.md, divergence.json");
            if !rep.causes.is_empty() {
                f.push_rule();
            }
            for c in rep.causes.iter().take(3) {
                f.push_kv(
                    &format!("{} {} ({:.2})", c.code, c.cause, c.score),
                    &c.evidence,
                );
                if !c.fix.is_empty() {
                    f.push_kv("    try", format!("--set {}", c.fix.join(" --set ")));
                }
            }
            line("");
            for l in f.render(style) {
                line(&l);
            }
            diag::err(rep.diagnostic())
        }
        Err(e) => diag::err(Diagnostic::new(
            diag::RUN_DIVERGED,
            format!("the solution diverged; the divergence report failed: {e:#}"),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounding_of_the_fixes() {
        assert_eq!(round2(0.7 * 0.7), 0.49);
        assert_eq!(sig3(0.000123456), "1.23e-4");
    }

    #[test]
    fn courant_of_a_uniform_flux() {
        use foehncraft_mesh::generate::{BoxSpec, Kind};
        let m = BoxSpec::new(Kind::Hex, 3, [2, 1, 1]).generate().unwrap();
        let g = Geometry::compute(&m);
        let mf = vec![1.0; m.n_faces()];
        let co = courant(&m, &g, &mf, 1.0, 0.5);
        // each cell has six faces with |ṁ| = 1: 0.5 · 6 / (2 V)
        for (c, co) in co.iter().enumerate() {
            assert!((co - 1.5 / g.cell_volume[c]).abs() < 1e-12);
        }
    }
}
