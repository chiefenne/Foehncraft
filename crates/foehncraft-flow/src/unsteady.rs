//! The transient `foehn flow` case (a case file with a `[time]` table):
//! the time loop with PISO ([`crate::transient`]), the adaptive time step,
//! restart files (D42), control files and re-reading the case file (D40),
//! probes, force histories and field series (D43), and the transient
//! verification: order sweeps in space (levels) and time (time steps)
//! against the Taylor–Green vortex, and the cylinder benchmark of Schäfer &
//! Turek (1996).
//!
//! Output in addition to the steady one: `history.csv` (per step: time, Δt,
//! Courant number, residuals, outer and CG iterations), `forces.csv` (force
//! per wall and the coefficients), `probes.csv`, `fields_<step>.vtkhdf`,
//! `restart.h5`.

use crate::analytic::Benchmark;
use crate::balance::{self, BalanceAction, BalanceSummary, Ledger};
use crate::bc::{self, BoundaryCondition, InletProfile, Level, boundary_extent, inlet_velocity};
use crate::case::{
    ErrorNorms, FlowCase, FlowReport, Log, Solve, boundary_summaries, dimension, error_norms,
    locate, mesh_for, write_fields, write_outputs,
};
use crate::control::Control;
use crate::engine::Steer;
use crate::monitor::{self, History, Probes};
use crate::plot::{self, Series};
use crate::restart::{self, Restart};
use crate::simplec::{self, Flow, Solution, Target};
use crate::transient::{ErrorControl, StepReport, TimeScheme, Transient};
use foehncraft_core::exact::ExactSumIter as _;

/// Times of an output schedule in simulated time (D66): the multiples of
/// `interval` and the listed `times`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Schedule<'a> {
    pub interval: Option<f64>,
    pub times: &'a [f64],
}

impl Schedule<'_> {
    /// The first scheduled time later than t + tol.
    pub fn next(&self, t: f64, tol: f64) -> Option<f64> {
        let a = self.interval.map(|i| (((t + tol) / i).floor() + 1.0) * i);
        self.times
            .iter()
            .copied()
            .filter(|&x| x > t + tol)
            .chain(a)
            .reduce(f64::min)
    }

    /// Whether a scheduled time lies in (t0, t1], to within tol.
    pub fn fired(&self, t0: f64, t1: f64, tol: f64) -> bool {
        self.next(t0, tol).is_some_and(|x| x <= t1 + tol)
    }

    /// Scheduled times in (t0, t1] (for the estimate of the output).
    pub fn count(&self, t0: f64, t1: f64) -> usize {
        let a = self.interval.map_or(0, |i| {
            ((t1 / i).floor() - (t0 / i).floor()).max(0.0) as usize
        });
        a + self.times.iter().filter(|&&x| x > t0 && x <= t1).count()
    }
}

/// The schedules of fields, monitors and restart files in simulated time.
fn schedules(case: &FlowCase) -> [Schedule<'_>; 3] {
    let (o, r) = (&case.output, &case.run);
    [
        Schedule {
            interval: o.fields_interval,
            times: &o.fields_times,
        },
        Schedule {
            interval: o.monitor_interval,
            times: &[],
        },
        Schedule {
            interval: r.restart_interval,
            times: &r.restart_times,
        },
    ]
}

/// Numbered field files and restart files that the schedules in simulated
/// time ask for over a run of `steps` steps of `dt` (the estimate, D81).
pub(crate) fn scheduled_outputs(case: &FlowCase, steps: usize, dt: f64) -> (usize, usize) {
    let Some(t) = case.time else {
        return (0, 0);
    };
    let (t0, t1) = (
        t.start_time,
        t.end_time.unwrap_or(t.start_time + steps as f64 * dt),
    );
    let numbered = Schedule {
        interval: None,
        times: &case.run.restart_times,
    };
    (schedules(case)[0].count(t0, t1), numbered.count(t0, t1))
}

/// Repetitions of one step at most when its error estimate exceeds the
/// tolerance (D65); the last try is kept.
const MAX_REJECTIONS: usize = 8;
use anyhow::{Context, Result, bail};
use foehncraft_fv::disc::Soa3;
use foehncraft_mesh::geometry::norm;
use foehncraft_mesh::locate::Locator;
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use std::path::Path;
use std::time::Instant;

/// Summary of a transient run.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TransientSummary {
    /// Time integration.
    pub scheme: String,
    /// Time steps done in this run.
    pub steps: usize,
    /// Step number at the end (with a restart: counted from the start of
    /// the first run).
    pub last_step: usize,
    /// Time at the end (s).
    pub time: f64,
    /// Smallest and largest time step (s).
    pub dt_min: f64,
    /// See `dt_min`.
    pub dt_max: f64,
    /// Largest Courant number.
    pub cfl_max: f64,
    /// Mean outer iterations per step.
    pub outer_mean: f64,
    /// Mean CG iterations per step.
    pub cg_mean: f64,
    /// Steps repeated with a shorter Δt because their time-integration
    /// error exceeded the tolerance (D65).
    #[serde(default)]
    pub rejected: usize,
    /// Restart file the run continued from.
    pub restart_from: Option<String>,
    /// Ended by `.stop`.
    pub stopped: bool,
    /// Nanoseconds per cell and time step.
    pub ns_per_cell_step: f64,
    /// The energy equation (stage 7d6): the largest relative heat
    /// imbalance of a step (|imbalance| over the largest heat flow, the
    /// `relative` column of `heat.csv`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_relative_max: Option<f64>,
    /// The thermal acceleration of the solids (D213).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub acceleration: Option<crate::acceleration::AccelerationSummary>,
    /// The cycle-averaged coupling loop of a periodic conjugate case
    /// (D204).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cycle_coupling: Option<crate::cycle::CycleSummary>,
    /// The frozen-flow march of a transient conjugate case (D118).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub frozen_flow: Option<crate::frozen::FrozenSummary>,
    /// The deforming mesh (D91, stage MS0).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mesh_motion: Option<crate::motion::MotionSummary>,
}

/// The heat balance of one time step (stage 7d6): the balance of the
/// energy equation and the heat flow through each interface of a
/// conjugate case (W, into the fluid, the order of `heat.csv`).
#[derive(Debug, Clone, Serialize)]
struct StepHeat {
    #[serde(flatten)]
    balance: crate::scalar::EnergyBalance,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    interfaces: Vec<f64>,
}

/// Time-dependent boundary faces: boundary index and the face.
struct AnalyticFaces {
    faces: Vec<(usize, usize)>,
    ub: [Vec<f64>; 3],
}

impl AnalyticFaces {
    fn new(case: &FlowCase, mesh: &foehncraft_mesh::Mesh, ub: &[Vec<f64>; 3]) -> Self {
        let ni = mesh.n_interior_faces();
        let mut faces = Vec::new();
        for p in &mesh.boundaries {
            if let Some(BoundaryCondition::VelocityInlet {
                profile: Some(InletProfile::Analytic),
                ..
            }) = case.boundary.get(&p.name)
            {
                faces.extend(p.range().map(|f| (f - ni, f)));
            }
        }
        Self {
            faces,
            ub: ub.clone(),
        }
    }

    /// Sets the analytic velocity of time `t` on the faces.
    fn update(
        &mut self,
        case: &FlowCase,
        mesh: &foehncraft_mesh::Mesh,
        geo: &foehncraft_mesh::Geometry,
        fl: &mut Flow,
        t: f64,
    ) -> Result<()> {
        let Some(a) = case.analytic_at(t) else {
            return Ok(());
        };
        if self.faces.is_empty() || !a.is_transient() {
            return Ok(());
        }
        for p in &mesh.boundaries {
            let Some(b) = case.boundary.get(&p.name) else {
                continue;
            };
            let ext = boundary_extent(mesh, p);
            for f in p.range() {
                if let Some((fb, _)) = self.faces.iter().find(|x| x.1 == f) {
                    let at = bc::At::numbers(Some(&a));
                    let v = inlet_velocity(b, None, geo, f, at, &ext).unwrap_or([0.0; 3]);
                    for (ub, vi) in self.ub.iter_mut().zip(v) {
                        ub[*fb] = vi;
                    }
                }
            }
        }
        fl.set_boundary_velocity(&self.ub)
    }
}

/// Boundary values and the body force that depend on t (D53): evaluated
/// on the host at the new time level before each step and uploaded.
struct TimeValues {
    values: bc::Compiled,
    boundary: bool,
    source: bool,
}

impl TimeValues {
    fn new(case: &FlowCase) -> Result<Self> {
        let values = case.compiled_boundaries()?;
        let (source, boundary) = crate::case::time_dependent(case, &values)?;
        Ok(Self {
            values,
            boundary,
            source,
        })
    }

    /// The values of time `t` on the faces and in the cells.
    fn update(
        &self,
        case: &FlowCase,
        mesh: &foehncraft_mesh::Mesh,
        geo: &foehncraft_mesh::Geometry,
        fl: &mut Flow,
        t: f64,
    ) -> Result<()> {
        if self.boundary {
            let a = case.analytic_at(t);
            let arr = bc::arrays(
                mesh,
                geo,
                &case.boundary,
                &fl.phys,
                bc::At {
                    analytic: a.as_ref(),
                    values: &self.values,
                    t,
                    turning: &case.turning(),
                    hydrostatic: crate::buoyancy::hydrostatic(case),
                },
            )?;
            fl.set_boundary_velocity(&arr.ub)?;
            fl.set_boundary_pressure(&arr.pset)?;
        }
        if self.source
            && let Some(src) = crate::case::momentum_source(case, geo, mesh.n_cells, t)?
        {
            fl.set_momentum_source(src)?;
        }
        Ok(())
    }
}

/// A run of the transient sweep: a refinement level and a time step.
#[derive(Debug, Clone, Copy)]
struct RunSpec {
    level: Option<u32>,
    dt: f64,
}

/// The result of one run that the sweeps need.
struct RunOut {
    solve: Solve,
    /// The start values of the regions and parts (D238).
    start_values: Vec<crate::start::StartValue>,
    u: Soa3,
    /// Final scalar fields, in the order of the case's scalar sections
    /// (on the full mesh in a conjugate case).
    scalars: Vec<Vec<f64>>,
    /// The cell volumes of `scalars`.
    scalar_volume: Vec<f64>,
    volume: Vec<f64>,
    stopped: bool,
}

/// Runs a transient case; `log` receives progress lines, `steer` carries
/// the requests, limits and run record (crate::engine).
pub fn run(
    case: &FlowCase,
    name: &str,
    out: &Path,
    steer: &mut Steer,
    log: impl FnMut(&str),
) -> Result<FlowReport> {
    let t0 = Instant::now();
    std::fs::create_dir_all(out)?;
    let control = Control::new(out)?;
    let mut log = log;
    let mut lg = Log::new(&mut log, case.log.style);
    let time = case.time.context("a transient case needs a [time] table")?;
    let v = &case.verification;
    let specs: Vec<RunSpec> = if !v.time_steps.is_empty() {
        let level = v.levels.first().copied();
        v.time_steps
            .iter()
            .map(|&dt| RunSpec { level, dt })
            .collect()
    } else if v.levels.is_empty() {
        vec![RunSpec {
            level: None,
            dt: time.time_step,
        }]
    } else {
        v.levels
            .iter()
            .map(|&l| RunSpec {
                level: Some(l),
                dt: time.time_step,
            })
            .collect()
    };
    if specs.len() > 1 && case.run.restart.is_some() {
        bail!("[run] restart continues a single run, not a sweep over levels or time steps");
    }
    let analytic = case.analytic_at(time.start_time);
    let reynolds = case.reynolds();
    let mut rep = FlowReport {
        name: name.into(),
        reynolds,
        findings: Vec::new(),
        solves: Vec::new(),
        orders: Vec::new(),
        time_orders: Vec::new(),
        scalar_orders: Vec::new(),
        uncertainty: Vec::new(),
        uncertainty_of: "",
        stopped: false,
        stop: None,
        warnings: Vec::new(),
        failures: Vec::new(),
        seconds: 0.0,
        pressure_values: case.pressure_frame.as_ref().map(|f| f.json()),
        start_values: Vec::new(),
    };
    std::fs::write(out.join("case.resolved.toml"), case.resolved_toml()?)?;
    if !case.run.ranks.is_empty() {
        crate::ranks::supported(case)?;
    }
    if steer.ranks.is_some() {
        crate::ranks::supported_run(case)?;
        anyhow::ensure!(
            specs.len() == 1,
            "a sweep over time steps does not run on ranks yet (MPI step W7)"
        );
    }
    let mut outs: Vec<RunOut> = Vec::new();
    let last = specs.len() - 1;
    for (ri, spec) in specs.iter().enumerate() {
        let (mesh, geo, desc) = mesh_for(case, spec.level)?;
        // conjugate heat transfer (stage 7d6): the flow runs on the fluid
        // view, the energy equation on the full mesh
        let (view_case, mesh, geo) = crate::cht::view(case, mesh, geo)?;
        let case = &view_case;
        let parts = crate::fluids::check(case, &mesh)?;
        let values = case.compiled_boundaries()?;
        let phys_on = crate::case::physics_on(case, &crate::case::physics(case)?, &mesh, &geo)?;
        let mut findings = bc::check(
            &mesh,
            &geo,
            &case.boundary,
            &phys_on,
            bc::At {
                analytic: analytic.as_ref(),
                values: &values,
                t: case.time.map_or(0.0, |t| t.start_time),
                turning: &case.turning(),
                hydrostatic: crate::buoyancy::hydrostatic(case),
            },
            case.check_context(),
        );
        findings.extend(crate::turbulence::setup_findings(case));
        findings.extend(crate::energy::setup_findings(case));
        findings.extend(crate::pressure_values::setup_findings(case));
        findings.extend(crate::buoyancy::setup_findings(case, &mesh, &geo));
        steer.set_cells(mesh.n_cells, usize::from(case.run.backend != "cpu"));
        if ri == 0 {
            rep.findings = findings.clone();
            lg.frame(&crate::logblocks::banner(steer.build.as_ref(), case));
            for n in std::mem::take(&mut steer.notes) {
                lg.note(&n);
            }
        }
        if findings.iter().any(|f| f.level == Level::Error) {
            let msg: Vec<String> = findings
                .iter()
                .filter(|f| f.level == Level::Error)
                .map(|f| f.to_string())
                .collect();
            lg.line("");
            for f in &findings {
                lg.line(&crate::case::finding_line(f));
            }
            std::fs::write(out.join("flow.log"), &lg.text)?;
            bail!(
                "the boundary conditions cannot be run:\n  {}",
                msg.join("\n  ")
            );
        }
        let (quality, _) = foehncraft_mesh::quality::QualityReport::compute(&mesh, &geo);
        if ri == 0 {
            let t = case.time.unwrap_or_default();
            let steps = t.steps.unwrap_or_else(|| {
                t.end_time
                    .map_or(0, |e| ((e - t.start_time) / spec.dt).ceil() as usize)
            });
            if specs.len() == 1 {
                lg.frame(&files_frame(case, &mesh, &desc, out, spec.dt));
            }
            lg.frame(&crate::logblocks::control_transient(case, spec.dt, steps));
        }
        lg.frame(&crate::logblocks::mesh(&quality, &desc));
        if ri == 0 {
            lg.frame(&crate::logblocks::fluid(case));
            lg.frame(&crate::logblocks::solution(case, true));
            lg.frame(&crate::logblocks::boundaries(case, &quality));
            if let Some(f) = crate::logblocks::limitations(case) {
                lg.frame(&f);
            }
            let mut any = false;
            for n in crate::logblocks::notes(case)
                .into_iter()
                .chain(crate::fluids::notes(&parts))
            {
                if !any {
                    lg.line("");
                    any = true;
                }
                lg.line(&n);
            }
            for f in &findings {
                if !any {
                    lg.line("");
                    any = true;
                }
                lg.line(&crate::case::finding_line(f));
            }
        }
        if specs.len() > 1 && ri == last {
            lg.frame(&files_frame(case, &mesh, &desc, out, spec.dt));
        }
        let label = if specs.len() > 1 {
            format!("{desc}, Δt = {:e} s", spec.dt)
        } else {
            desc.clone()
        };
        let r = run_one(
            case,
            name,
            out,
            &mut lg,
            &control,
            steer,
            &mut rep.warnings,
            &mesh,
            &geo,
            &label,
            spec.dt,
            ri == last,
            specs.len() == 1,
        )?;
        rep.stopped |= r.stopped;
        if rep.start_values.is_empty() {
            rep.start_values = r.start_values.clone();
        }
        outs.push(r);
        if rep.stopped {
            rep.stop = steer.stop.clone();
            break;
        }
    }
    // sweeps
    let v = &case.verification;
    if !v.time_steps.is_empty() && outs.len() == specs.len() {
        // differences of successive solutions on one mesh
        let diffs: Vec<f64> = outs
            .windows(2)
            .map(|w| {
                let (a, b) = (&w[0], &w[1]);
                let (mut e2, mut vol) = (0.0, 0.0);
                for c in 0..a.volume.len() {
                    let d = [
                        a.u.x[c] - b.u.x[c],
                        a.u.y[c] - b.u.y[c],
                        a.u.z[c] - b.u.z[c],
                    ];
                    e2 += (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) * a.volume[c];
                    vol += a.volume[c];
                }
                (e2 / vol).sqrt()
            })
            .collect();
        let dts = &v.time_steps;
        // the scalars: the same differences of successive solutions
        for (k, sm) in outs[0].solve.scalars.iter().enumerate() {
            let d: Vec<f64> = outs
                .windows(2)
                .map(|w| {
                    let (a, b) = (&w[0], &w[1]);
                    let (mut e2, mut vol) = (0.0, 0.0);
                    for c in 0..a.scalar_volume.len() {
                        let x = a.scalars[k][c] - b.scalars[k][c];
                        e2 += x * x * a.scalar_volume[c];
                        vol += a.scalar_volume[c];
                    }
                    (e2 / vol).sqrt()
                })
                .collect();
            let o: Vec<f64> = d
                .windows(2)
                .enumerate()
                .map(|(i, d)| (d[0] / d[1]).ln() / (dts[i] / dts[i + 1]).ln())
                .collect();
            lg.note(&format!(
                "time sweep: scalar {}: observed time orders {o:.3?}",
                sm.name
            ));
            if let (Some([lo, hi]), Some(&last)) = (v.scalar_order, o.last())
                && !(lo..=hi).contains(&last)
            {
                rep.failures.push(format!(
                    "scalar {}: observed time order {last:.3} on the smallest steps, outside [{lo}, {hi}]",
                    sm.name
                ));
            }
            rep.scalar_orders.push((format!("{} (time)", sm.name), o));
        }
        // a fluid at rest (the conduction of a conjugate case) gives no
        // velocity differences and no order of them
        if diffs.iter().any(|&d| d > 0.0) {
            rep.time_orders = diffs
                .windows(2)
                .enumerate()
                .map(|(i, d)| (d[0] / d[1]).ln() / (dts[i] / dts[i + 1]).ln())
                .collect();
            lg.note(&format!(
                "time sweep: differences of successive solutions [{}], observed time orders {:.3?}",
                diffs
                    .iter()
                    .map(|d| format!("{d:.3e}"))
                    .collect::<Vec<_>>()
                    .join(", "),
                rep.time_orders
            ));
            if let (Some([lo, hi]), Some(&o)) = (v.expected_time_order, rep.time_orders.last())
                && !(lo..=hi).contains(&o)
            {
                rep.failures.push(format!(
                    "observed time order {o:.3} on the smallest steps, outside [{lo}, {hi}]"
                ));
            }
            let pts: Vec<(f64, f64)> = dts[1..]
                .iter()
                .copied()
                .zip(diffs.iter().copied())
                .collect();
            let slope: Vec<(f64, f64)> = match (pts.first(), pts.last()) {
                (Some(&a), Some(&b)) => vec![(a.0, b.1 * (a.0 / b.0).powi(2)), b],
                _ => Vec::new(),
            };
            std::fs::write(
                out.join("time_convergence.svg"),
                plot::chart(
                    &format!("{name}: change of the solution over the time step"),
                    "Δt (s)",
                    "‖u(Δt) − u(2Δt)‖ (m/s)",
                    &[
                        Series {
                            label: "Foehncraft",
                            points: &pts,
                            markers: true,
                        },
                        Series {
                            label: "slope 2",
                            points: &slope,
                            markers: false,
                        },
                    ],
                    true,
                    true,
                ),
            )?;
        }
    }
    if specs.len() >= 3 && outs.len() == specs.len() {
        let js: Vec<serde_json::Value> = outs
            .iter()
            .map(|o| serde_json::to_value(&o.solve).unwrap_or_default())
            .collect();
        let (sizes, of): (Vec<f64>, &'static str) = if v.time_steps.is_empty() {
            let d = f64::from(dimension(case));
            (
                outs.iter()
                    .map(|o| (o.solve.cells as f64).powf(-1.0 / d))
                    .collect(),
                "meshes",
            )
        } else {
            (specs.iter().map(|s| s.dt).collect(), "time steps")
        };
        rep.uncertainty = crate::uncertainty::sweep(&js, &sizes);
        rep.uncertainty_of = of;
    }
    // the velocity against [analytic] (the scalars are checked below)
    if v.levels.len() >= 2
        && v.time_steps.is_empty()
        && outs.len() == specs.len()
        && outs.iter().all(|o| o.solve.errors.is_some())
    {
        let e: Vec<&ErrorNorms> = outs
            .iter()
            .filter_map(|s| s.solve.errors.as_ref())
            .collect();
        rep.orders = e
            .windows(2)
            .map(|w| (w[0].l2 / w[1].l2).ln() / (w[0].h / w[1].h).ln())
            .collect();
        lg.note(&format!("observed orders in space {:.3?}", rep.orders));
        if let (Some([lo, hi]), Some(&o)) = (v.expected_order, rep.orders.last())
            && !(lo..=hi).contains(&o)
        {
            rep.failures.push(format!(
                "observed order {o:.3} on the finest pair, outside [{lo}, {hi}]"
            ));
        }
        let pts: Vec<(f64, f64)> = e.iter().map(|e| (e.h, e.l2)).collect();
        let slope: Vec<(f64, f64)> = match (pts.first(), pts.last()) {
            (Some(&a), Some(&b)) => vec![(a.0, b.1 * (a.0 / b.0).powi(2)), b],
            _ => Vec::new(),
        };
        std::fs::write(
            out.join("convergence.svg"),
            plot::chart(
                &format!("{name}: velocity error over mesh size"),
                "h (m)",
                "L2 error / U",
                &[
                    Series {
                        label: "Foehncraft",
                        points: &pts,
                        markers: true,
                    },
                    Series {
                        label: "slope 2",
                        points: &slope,
                        markers: false,
                    },
                ],
                true,
                true,
            ),
        )?;
    }
    for o in &outs {
        let s = &o.solve;
        if let (Some(max), Some(e)) = (v.max_error, &s.errors)
            && e.l2 > max
        {
            rep.failures.push(format!(
                "{}: relative L2 velocity error {:.3e} above {max:e}",
                s.mesh, e.l2
            ));
        }
        if s.mass_imbalance.abs() > 1e-6 {
            rep.failures.push(format!(
                "{}: net mass flow out of the domain {:.2e} of the largest boundary flow",
                s.mesh, s.mass_imbalance
            ));
        }
        if let Some(f) = s.chemistry.as_ref().and_then(|c| c.failure()) {
            rep.failures.push(format!("{}: {f}", s.mesh));
        }
        if let (Some(tol), Some(r)) = (v.reference_tolerance, &s.reference)
            && r.u_max_dev.max(r.v_max_dev) > tol
        {
            rep.failures.push(format!(
                "{}: deviation from {} of {:.4} above {tol}",
                s.mesh,
                r.reference,
                r.u_max_dev.max(r.v_max_dev)
            ));
        }
        if let (Some(tol), Some(b)) = (v.benchmark_tolerance, &s.benchmark)
            && b.max_deviation() > tol
        {
            rep.failures.push(format!(
                "{}: largest distance from the {} ranges {:.4} above {tol}",
                s.mesh,
                b.benchmark,
                b.max_deviation()
            ));
        }
        if v.benchmark == Some(Benchmark::SchaferTurek2d2) && s.benchmark.is_none() && !o.stopped {
            rep.failures.push(format!(
                "{}: fewer than two full periods of the lift in the evaluation window",
                s.mesh
            ));
        }
    }
    rep.solves = outs.into_iter().map(|o| o.solve).collect();
    // the scalars with an exact solution: the largest error, the orders
    // of a sweep over levels (those of a time sweep are above)
    let levels = if v.time_steps.is_empty() {
        v.levels.len()
    } else {
        0
    };
    crate::scalar::check_sweep(&mut rep, v, levels, &mut |l| lg.line(l));
    if !case.run.ranks.is_empty() {
        lg.line("");
        let failures = crate::ranks::report(case, out, &mut |l| lg.line(l))?;
        rep.failures.extend(failures);
    }
    rep.seconds = t0.elapsed().as_secs_f64();
    crate::case::end_frames(&mut lg, &rep);
    lg.sync(&mut steer.record);
    steer.record.flush();
    std::fs::write(out.join("flow.log"), &lg.text)?;
    std::fs::write(out.join("flow.json"), serde_json::to_string_pretty(&rep)?)?;
    std::fs::write(out.join("flow.md"), rep.to_markdown())?;
    Ok(rep)
}

/// The files of a transient run: the inputs, the output directory and what
/// will be written there with its estimated size.
fn files_frame(
    case: &FlowCase,
    mesh: &foehncraft_mesh::Mesh,
    desc: &str,
    out: &std::path::Path,
    dt: f64,
) -> crate::console::Frame {
    use crate::console::size;
    let t = case.time.unwrap_or_default();
    let steps = t.steps.unwrap_or_else(|| {
        t.end_time
            .map_or(0, |e| ((e - t.start_time) / dt).ceil() as usize)
    });
    let nc = mesh.n_cells;
    let o = &case.output;
    let mut f = crate::logblocks::files(case, desc, out);
    f.push_rule();
    f.push_kv(
        "history.csv",
        format!(
            "{} rows, about {}",
            steps / o.history_every + 1,
            size(steps / o.history_every * 120)
        ),
    );
    let walls = case
        .boundary
        .values()
        .filter(|b| matches!(b, BoundaryCondition::Wall { .. }))
        .count();
    if walls > 0 || o.coefficients.is_some() {
        f.push_kv(
            "forces.csv",
            format!(
                "{} rows (every {} steps), forces on {walls} wall(s){}",
                steps / o.monitor_every + 1,
                o.monitor_every,
                if o.coefficients.is_some() {
                    " and the force coefficients"
                } else {
                    ""
                }
            ),
        );
    }
    if !o.probes.is_empty() {
        f.push_kv(
            "probes.csv",
            format!(
                "{} probe(s), {} rows",
                o.probes.len(),
                steps / o.monitor_every + 1
            ),
        );
    }
    let bytes = nc * (8 * 6) + mesh.face_nodes.len() * 8 + mesh.points.len() * 24;
    if o.vtkhdf {
        let series = steps.checked_div(o.fields_every).unwrap_or(0)
            + o.fields_at.len()
            + scheduled_outputs(case, steps, dt).0;
        f.push_kv(
            "fields.vtkhdf (final)",
            format!(
                "about {} each{}",
                size(bytes),
                if series > 0 {
                    format!(", and {series} fields_<step>.vtkhdf")
                } else {
                    String::new()
                }
            ),
        );
    }
    let r = &case.run;
    if r.restart_final
        || r.restart_every > 0
        || !r.restart_at.is_empty()
        || r.restart_interval.is_some()
        || !r.restart_times.is_empty()
    {
        f.push_kv(
            "restart.h5",
            format!(
                "about {}{}",
                size(nc * 8 * 22 + mesh.n_faces() * 8 * 4),
                if r.restart_at.is_empty() {
                    String::new()
                } else {
                    format!(", and restart_<step>.h5 at {:?}", r.restart_at)
                }
            ),
        );
    }
    f.push_kv("flow.md, flow.json", "report and its data");
    f.push_kv("result.json", "summary for programs (exit code, checks)");
    f.push_kv("balance.csv, balance.svg", "mass and energy ledger");
    f.push_kv(
        "boundaries.csv",
        "area, mass flow, pressure and force per boundary",
    );
    f
}

/// The restart content of the current state.
#[allow(clippy::too_many_arguments)]
fn restart_of(
    fl: &mut Flow,
    tr: &Transient,
    mesh: &foehncraft_mesh::Mesh,
    name: &str,
    dt_next: f64,
    cycle: Option<&crate::cycle::Controller>,
    frozen: Option<&crate::frozen::Controller>,
    on: Option<&crate::ranks::OnRanks>,
) -> Result<Restart> {
    let r = restart_part(fl, tr, mesh, name, dt_next, cycle, frozen)?;
    match on {
        // on ranks gathered into the global order (the bytes of the
        // one-rank file, W6; collective)
        Some(o) => o.restart(r),
        None => Ok(r),
    }
}

/// The restart of this rank's part (`mesh`: the mesh of its solver).
fn restart_part(
    fl: &mut Flow,
    tr: &Transient,
    mesh: &foehncraft_mesh::Mesh,
    name: &str,
    dt_next: f64,
    cycle: Option<&crate::cycle::Controller>,
    frozen: Option<&crate::frozen::Controller>,
) -> Result<Restart> {
    let cycle = match cycle {
        Some(c) => Some(c.save()?),
        None => None,
    };
    let frozen = match frozen {
        Some(c) => Some(c.save()?),
        None => None,
    };
    Ok(Restart {
        cht: crate::cht::restart_part(fl, cycle, frozen)?,
        spray: None,
        case: name.into(),
        cell_id: mesh.cell_id[..mesh.n_cells].to_vec(),
        face_id: mesh.face_id.clone(),
        state: fl.state()?,
        iteration: fl.iteration(),
        amg_setups: fl.amg_setups,
        amg_matrix: fl.amg_matrix.clone(),
        sweep: Default::default(),
        transient: Some((tr.old_levels(fl)?, tr.clock, dt_next, tr.scheme)),
        geometry: None,
    })
}

/// The mesh and the geometry of the run at this moment: those of the
/// rotor angle with a sliding mesh by imprinting, the base mesh otherwise.
fn current<'a>(
    work: &'a Option<crate::sliding::Working>,
    deform: &'a Option<crate::motion::Deforming>,
    mesh: &'a foehncraft_mesh::Mesh,
    geo: &'a foehncraft_mesh::Geometry,
) -> (&'a foehncraft_mesh::Mesh, &'a foehncraft_mesh::Geometry) {
    match (work, deform) {
        (Some(w), _) => (&w.mesh, &w.geo),
        (None, Some(d)) => (&d.mesh, &d.geo),
        (None, None) => (mesh, geo),
    }
}

/// A new flow solver on `mesh` built as at the start of the run, with the
/// settings of `fl` (its fields not yet loaded) and the density of the
/// fluid.
#[allow(clippy::too_many_arguments)]
fn new_flow(
    case: &FlowCase,
    tv: &TimeValues,
    target: &Target,
    reference_cell: Option<usize>,
    (mesh, geo): (&foehncraft_mesh::Mesh, &foehncraft_mesh::Geometry),
    t: f64,
    fl: &Flow,
    tr: &Transient,
    wall_distance: Option<Vec<f64>>,
) -> Result<(Flow, Transient, f64)> {
    let phys = crate::case::physics(case)?;
    let phys = crate::case::physics_on(case, &phys, mesh, geo)?;
    let rho = phys.rho;
    let analytic = case.analytic_at(t);
    let arrays = bc::arrays(
        mesh,
        geo,
        &case.boundary,
        &phys,
        bc::At {
            analytic: analytic.as_ref(),
            values: &tv.values,
            t,
            turning: &case.turning(),
            hydrostatic: crate::buoyancy::hydrostatic(case),
        },
    )?;
    let disc = simplec::discretisation(mesh, geo, &arrays)?;
    let mut nfl = Flow::on(
        disc,
        arrays,
        phys,
        fl.settings.clone(),
        true,
        reference_cell,
        target,
    )?;
    if let Some(src) = crate::case::momentum_source(case, geo, mesh.n_cells, t)? {
        nfl.set_momentum_source(src)?;
    }
    crate::scalar::attach(&mut nfl, case, mesh, geo)?;
    crate::density::attach(&mut nfl, case)?;
    crate::buoyancy::attach(&mut nfl, case)?;
    crate::turbulence::attach_with(&mut nfl, case, mesh, geo, wall_distance)?;
    crate::frame::attach(&mut nfl, case, mesh, geo)?;
    let scheme = case.time.context("[time]")?.scheme;
    let ntr = Transient::new(&mut nfl, scheme, tr.clock.time)?;
    Ok((nfl, ntr, rho))
}

/// The flow solver on the mesh of the next rotor angle (sliding mesh by
/// imprinting, D182 variant B), or on the deformed mesh (no `map`): built
/// again on `next` as at the start of the run, with the cell fields, the
/// old time levels and the fluxes of the faces that stay taken over from
/// `fl` and `tr` (the faces that the imprint made new start from the
/// interpolated flux). `t` is the time of the end of the step.
///
/// A deforming mesh (stage MS0) passes no `map`: the faces are the same,
/// and the variable density goes along with its mass and old levels.
#[allow(clippy::too_many_arguments)]
fn rebuild_flow(
    case: &FlowCase,
    tv: &TimeValues,
    target: &Target,
    reference_cell: Option<usize>,
    (mesh, geo): (&foehncraft_mesh::Mesh, &foehncraft_mesh::Geometry),
    map: Option<&crate::sliding::FaceMap>,
    t: f64,
    fl: &mut Flow,
    tr: &Transient,
    wall_distance: Option<Vec<f64>>,
) -> Result<(Flow, Transient)> {
    let (mut nfl, mut ntr, rho) = new_flow(
        case,
        tv,
        target,
        reference_cell,
        (mesh, geo),
        t,
        fl,
        tr,
        wall_distance,
    )?;
    let state = fl.state()?;
    let old = tr.old_levels(fl)?;
    let (state, old) = match map {
        Some(map) => crate::sliding::carry(state, old, map, &nfl.disc, rho, nfl.phys.mu),
        None => (state, old),
    };
    take_over(fl, tr, &mut nfl, &mut ntr, &state, &old, map.is_none())?;
    Ok((nfl, ntr))
}

/// The sliding mesh by imprinting on ranks (MPI step W7): the solver of
/// this rank on its part of `next`, the mesh of the next rotor angle, and
/// the part itself. The partition stays that of the start (the cells keep
/// their numbers); the state and the old levels are gathered into the
/// global order, carried onto the faces of `next` there by the code of
/// one rank (bitwise its result) and handed back to the ranks. Collective.
#[allow(clippy::too_many_arguments)]
fn rebuild_on_ranks(
    case: &FlowCase,
    tv: &TimeValues,
    target: &Target,
    reference_cell: Option<usize>,
    next: &crate::sliding::Working,
    map: &crate::sliding::FaceMap,
    t: f64,
    fl: &mut Flow,
    tr: &Transient,
    wall_distance: Option<Vec<f64>>,
    o: &crate::ranks::OnRanks,
) -> Result<(crate::ranks::OnRanks, Flow, Transient)> {
    let (mesh, geo) = (&next.mesh, &next.geo);
    // the boundary arrays of the whole mesh, for the output and the carry
    let phys = crate::case::physics(case)?;
    let gphys = crate::case::physics_on(case, &phys, mesh, geo)?;
    let analytic = case.analytic_at(t);
    let at = bc::At {
        analytic: analytic.as_ref(),
        values: &tv.values,
        t,
        turning: &case.turning(),
        hydrostatic: crate::buoyancy::hydrostatic(case),
    };
    let garrays = bc::arrays(mesh, geo, &case.boundary, &gphys, at)?;
    let no =
        crate::ranks::OnRanks::imprinted(o.comm.clone(), mesh, geo, garrays, Some(o.part.clone()))?;
    let (mut nfl, mut ntr, rho) = new_flow(
        case,
        tv,
        &no.target(target),
        reference_cell,
        (&no.local.mesh, &no.geo),
        t,
        fl,
        tr,
        wall_distance.map(|d| no.own(&d)),
    )?;
    // the state and the old levels of the whole mesh (collective)
    let mut r = restart_part(fl, tr, &o.local.mesh, "", 0.0, None, None)?;
    r.amg_matrix = None;
    let mut r = o.restart(r)?;
    let (old, clock, dt_next, scheme) = r.transient.take().context("the old levels")?;
    let (state, old) = crate::sliding::carry(r.state.clone(), old, map, &no.disc, rho, nfl.phys.mu);
    r.state = state;
    r.transient = Some((old, clock, dt_next, scheme));
    r.cell_id = mesh.cell_id[..mesh.n_cells].to_vec();
    r.face_id = mesh.face_id.clone();
    // this rank's part of them
    let mut r = no.scatter(r)?;
    let (old, ..) = r.transient.take().context("the old levels")?;
    take_over(fl, tr, &mut nfl, &mut ntr, &r.state, &old, false)?;
    Ok((no, nfl, ntr))
}

/// Loads the carried `state` and old levels `old` into the new solver and
/// takes over what goes on from the old one (`density`: the variable
/// density goes along, a deforming mesh).
fn take_over(
    fl: &mut Flow,
    tr: &Transient,
    nfl: &mut Flow,
    ntr: &mut Transient,
    state: &crate::simplec::State,
    old: &crate::transient::OldLevels,
    density: bool,
) -> Result<()> {
    nfl.load(state)?;
    if density && let Some(mut d) = fl.density.take() {
        let tc = nfl.fluid_temperature()?.context("the temperature")?;
        let tb = nfl.fluid_boundary_temperature(&tc)?;
        d.set_geometry(&nfl.disc, &tc, &tb);
        nfl.density = Some(d);
        nfl.refresh_density()?;
    }
    // (the floor and cap counts of the turbulence model go on, D69)
    if let (Some(a), Some(b)) = (fl.turbulence.as_ref(), nfl.turbulence.as_mut()) {
        b.health = a.health.clone();
    }
    ntr.restore(nfl, old, tr.clock)?;
    nfl.iteration = fl.iteration;
    // (the multigrid setups restart with the new matrix: run_one adds up the
    // setups of the solvers it replaces)
    nfl.timings = fl.timings.clone();
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_one<F: FnMut(&str)>(
    case0: &FlowCase,
    name: &str,
    out: &Path,
    lg: &mut Log<F>,
    control: &Control,
    steer: &mut Steer,
    warnings: &mut Vec<crate::diag::Diagnostic>,
    mesh: &foehncraft_mesh::Mesh,
    geo: &foehncraft_mesh::Geometry,
    desc: &str,
    dt0: f64,
    last: bool,
    single: bool,
) -> Result<RunOut> {
    let mut case = case0.clone();
    let t_setup = Instant::now();
    let (mesh0, geo0) = (mesh, geo);
    let time0 = case.time.context("[time]")?;
    // the sliding mesh by imprinting (D182 variant B): the mesh of the rotor
    // angle of this step replaces the base mesh in every step
    let imprint = crate::sliding::attach_imprint(&case, mesh0)?;
    // multigrid setups of the solvers that the imprint replaced (amg, sweep)
    let mut setups_past = (0usize, 0usize);
    let mut work: Option<crate::sliding::Working> = match &imprint {
        Some((imp, notes)) => {
            let w = imp.working(mesh0, geo0, time0.start_time)?;
            for l in &notes.lines {
                lg.line(l);
            }
            let defect = imp.closure(&w);
            if defect > 1e-8 {
                let m = format!(
                    "the cells of the imprinted interface close to {defect:.1e} only (|sum of S| / \
                     sum of |S|): the faces of the interface are not those of a mesh extruded \
                     along the axis, the mesh flux of the rigid turn leaves a mass source of \
                     this size"
                );
                lg.warn(&format!("{}: {m}", crate::diag::SLIDING_CLOSURE));
                warnings.push(crate::diag::Diagnostic::new(
                    crate::diag::SLIDING_CLOSURE,
                    m,
                ));
            }
            Some(w)
        }
        None => None,
    };
    // the deforming mesh (D91, stage MS0): the nodes at the start time; the
    // mesh of the end of each step replaces it before the step
    let mut deform = crate::motion::Deforming::new(&case, mesh0, time0.start_time)?;
    let moving = imprint.is_some() || deform.is_some();
    if let Some(d) = &deform {
        lg.line(&format!(
            "deforming mesh (D91, stage MS0): the nodes move by [mesh_motion] displacement; volume \
             {:.6e} m³ at t = {:.6e} s; the solver is built again on the moved mesh in every step",
            d.volume(),
            time0.start_time
        ));
    }
    let (mesh, geo) = current(&work, &deform, mesh0, geo0);
    let phys = crate::case::physics(&case)?;
    let phys = crate::case::physics_on(&case, &phys, mesh, geo)?;
    let analytic0 = case.analytic_at(time0.start_time);
    let tv = TimeValues::new(&case)?;
    let arrays = bc::arrays(
        mesh,
        geo,
        &case.boundary,
        &phys,
        bc::At {
            analytic: analytic0.as_ref(),
            values: &tv.values,
            t: time0.start_time,
            turning: &case.turning(),
            hydrostatic: crate::buoyancy::hydrostatic(&case),
        },
    )?;
    // a run on ranks (MPI step W7): this rank's part, its own arrays and
    // physics; the global ones stay for the output
    let mut on = match steer.ranks.clone() {
        Some(c) => {
            crate::ranks::supported_mesh(&case, mesh)?;
            Some(if work.is_some() {
                // the geometry of the rotor angle (rigid turn of the base
                // mesh), the partition kept for the whole run
                crate::ranks::OnRanks::imprinted(c, mesh, geo, arrays.clone(), None)?
            } else {
                crate::ranks::OnRanks::new(c, mesh, geo, arrays.clone())?
            })
        }
        None => None,
    };
    let lead = on.as_ref().is_none_or(crate::ranks::OnRanks::lead);
    if let Some(o) = &on {
        o.write_partition(out, mesh)?;
        let own = o.local.mesh.n_cells;
        let fewest = o.comm.min(own as f64) as usize;
        let most = o.comm.max(own as f64) as usize;
        lg.note(&format!(
            "solved on {} ranks ({} communicator), {fewest} to {most} cells each; rank 0 \
             writes the output, computed on the fields of all ranks in the order of the \
             whole mesh (MPI step W7)",
            o.comm.size(),
            o.comm.name()
        ));
    }
    let (arrays, phys) = match &on {
        Some(o) => {
            let lp = crate::case::physics(&case)?;
            let lp = crate::case::physics_on(&case, &lp, &o.local.mesh, &o.geo)?;
            let la = bc::arrays(
                &o.local.mesh,
                &o.geo,
                &case.boundary,
                &lp,
                bc::At {
                    analytic: analytic0.as_ref(),
                    values: &tv.values,
                    t: time0.start_time,
                    turning: &case.turning(),
                    hydrostatic: crate::buoyancy::hydrostatic(&case),
                },
            )?;
            (la, lp)
        }
        None => (arrays, phys),
    };
    // the mesh and geometry of this rank's solver
    let (smesh, sgeo) = match &on {
        Some(o) => (&o.local.mesh, &o.geo),
        None => (mesh, geo),
    };
    // the boundary kinds of the whole mesh, for the ledger
    let mom_kind = on
        .as_ref()
        .map_or_else(|| arrays.mom_kind.clone(), |o| o.arrays.mom_kind.clone());
    let ub0 = arrays.ub.clone();
    let disc = simplec::discretisation(smesh, sgeo, &arrays)?;
    let operators = crate::operators::report(
        &case,
        mesh,
        geo,
        on.as_ref().map_or(&arrays.u_points, |o| &o.arrays.u_points),
        out,
        lg,
        warnings,
    );
    // (on the base mesh: a probe or the reference cell is found there, a
    // point inside a sliding rotor is a limit of variant B)
    let loc = Locator::new(mesh0, geo0);
    let mut spray = crate::spray::create(&case, mesh0, geo0, &loc)?;
    // chemistry in the cells (SF2b): the reaction step of each time step
    let mut chem = crate::chemistry::flow::Chemistry::of(&case, &geo0.cell_volume)?;
    // the datums of the reports (SP1b, D262): the parcels see the same ones
    let datum_built = crate::datum::build(&case, mesh0, geo0)?;
    if let Some(sp) = spray.as_mut() {
        sp.set_datums(datum_built.iter().map(|b| b.parcels.clone()).collect());
    }
    let reference_cell = if on
        .as_ref()
        .map_or(arrays.has_pressure(), |o| o.arrays.has_pressure())
    {
        None
    } else {
        let x = case.pressure.reference_point.unwrap_or_else(|| {
            let (lo, hi) = mesh.bounding_box();
            std::array::from_fn(|i| lo[i] + 1e-6 * (hi[i] - lo[i]))
        });
        Some(locate(&loc, geo0, x))
    };
    let target = Target {
        device: case.run.backend.clone(),
        threads: case.run.threads,
        ptx: case.run.ptx.clone(),
        ranks: on.as_ref().and_then(|o| o.target(&Target::default()).ranks),
    };
    let mut fl = Flow::on(
        disc,
        arrays,
        phys,
        case.solver.clone(),
        true,
        reference_cell,
        &target,
    )?;
    if let Some(o) = &on {
        // the numbers of the whole mesh, not of this rank's part
        let skew = simplec::face_skewness(&o.disc);
        let lim = fl.remediation.skewness_limit;
        (fl.remediation.skewed_faces, fl.remediation.poor_cells) =
            simplec::skew_counts(&o.disc, &skew, lim);
    }
    let remediation = crate::case::report_remediation(&fl, lg, warnings);
    let nc = mesh.n_cells;
    if let Some(src) = crate::case::momentum_source(&case, sgeo, smesh.n_cells, time0.start_time)? {
        fl.set_momentum_source(src)?;
    }
    crate::scalar::attach(&mut fl, &case, smesh, sgeo)?;
    crate::density::attach(&mut fl, &case)?;
    crate::buoyancy::attach(&mut fl, &case)?;
    crate::porous::attach(&mut fl, &case, smesh)?;
    crate::vof::attach(&mut fl, &case, smesh, sgeo)?;
    for d in crate::immersed::attach(&mut fl, &case, smesh, sgeo)?.0 {
        if d.severity == crate::diag::Severity::Warning {
            lg.warn(&format!("{}: {}", d.code, d.message));
            warnings.push(d);
        } else {
            lg.note(&d.message);
        }
    }
    // the wall distance of the sliding mesh by imprinting: kept per part and
    // updated for the rotor angle of every step (D182, D48)
    let mut turning_distance = match (&imprint, &work) {
        (Some((imp, _)), Some(w)) => match imp.wall_distance(&case, mesh0, geo0) {
            Some(mut twd) => {
                lg.line(&format!(
                    "wall distance (D48) of the turning rotor: kept per part, {} of {} cells \
                     nearer to the edge of their part than to its walls search the other part \
                     in every step",
                    twd.candidates(),
                    mesh0.n_cells
                ));
                let d = twd.at_angles(&w.angles, &w.geo.cell_centroid[..w.mesh.n_cells]);
                Some((twd, d))
            }
            None => None,
        },
        _ => None,
    };
    match &on {
        // the wall distance of the whole mesh, each rank its cells (W5)
        Some(o) if !case.turbulence.is_laminar() => {
            let d = match turning_distance.as_mut() {
                Some((_, d)) => std::mem::take(d),
                None => crate::derived::wall_distance(&case, mesh, geo)?,
            };
            crate::turbulence::attach_with(&mut fl, &case, smesh, sgeo, Some(o.own(&d)))?;
        }
        _ => crate::turbulence::attach_with(
            &mut fl,
            &case,
            smesh,
            sgeo,
            turning_distance.as_mut().map(|(_, d)| std::mem::take(d)),
        )?,
    }
    // the sliding mesh (D182): the rotating zone and the shift of its data
    // (on ranks the plan of the whole mesh)
    let slide =
        crate::sliding::attach(&case, mesh, on.as_ref().map_or(&fl.disc, |o| &o.disc), dt0)?;
    // the rotating zones: the frames of `[[mrf]]` (a transient frozen rotor,
    // D185), the grid flux of the sliding rotors (D182) and the rotating
    // frame of the whole domain (`flow.frame_rotation`)
    if (slide.is_some() || !case.mrf.is_empty() || case.flow.frame_rotation != [0.0; 3])
        && let Some(notes) = crate::frame::attach(&mut fl, &case, smesh, sgeo)?
    {
        for l in &notes.lines {
            lg.line(l);
        }
        warnings.extend(notes.diagnostics);
    }
    let slide = slide.map(|(slide, notes)| {
        for l in &notes.lines {
            lg.line(l);
        }
        slide
    });
    let mut tr = Transient::new(&mut fl, time0.scheme, time0.start_time)?;
    if let Some(sl) = &slide {
        sl.install(&mut fl, &mut tr, on.as_ref())?;
    }
    let mut af = AnalyticFaces::new(&case, smesh, &ub0);
    let mut dt = dt0;
    let mut restart_from = None;
    let mut start_values = Vec::new();
    let mut cycle_state: Option<Vec<u8>> = None;
    let mut frozen_state: Option<Vec<u8>> = None;
    let restart_file = if single {
        case.run.restart.clone()
    } else {
        None
    };
    if let Some(path) = &restart_file {
        let r = restart::read(path)?;
        let (r, same) = match &on {
            // every rank its part, by global identity (W6)
            Some(o) => (o.scatter(r)?, true),
            None => r.to_local(&mesh.cell_id[..nc], &mesh.face_id, mesh.n_interior_faces())?,
        };
        if let Some(sp) = spray.as_mut() {
            match &r.spray {
                Some(b) => {
                    let env = foehncraft_particles::Env {
                        mesh: mesh0,
                        geo: geo0,
                        loc: &loc,
                    };
                    sp.restore(&env, b).map_err(|e| {
                        crate::diag::err(crate::diag::Diagnostic::new(
                            crate::diag::SPRAY_CONFIG,
                            format!(
                                "{}: the spray state of the restart file: {e:#}",
                                path.display()
                            ),
                        ))
                    })?;
                }
                None => {
                    lg.note("the restart file holds no spray state: the run starts without parcels")
                }
            }
        }
        if let Some((old, clock, dt_next, scheme)) = r.transient.clone() {
            if scheme != time0.scheme {
                lg.note(&format!(
                    "the restart file was written with {}, the case runs {}",
                    scheme.label(),
                    time0.scheme.label()
                ));
            }
            fl.load(&r.state)?;
            tr.restore(&mut fl, &old, clock)?;
            crate::cht::load_restart(&mut fl, r.cht.as_ref(), path)?;
            cycle_state = r.cht.as_ref().and_then(|c| c.cycle.clone());
            frozen_state = r.cht.as_ref().and_then(|c| c.frozen.clone());
            if (cycle_state.is_some() || frozen_state.is_some()) && !same {
                bail!(
                    "{}: the restart file holds the state of the {}, which follows the cell \
                     and face order of the run that wrote it; the case has another cell or face \
                     order",
                    path.display(),
                    if cycle_state.is_some() {
                        "cycle-averaged loop"
                    } else {
                        "frozen-flow march"
                    }
                );
            }
            fl.iteration = r.iteration;
            if let Some(m) = &r.amg_matrix {
                fl.restore_amg(m, r.amg_setups)?;
            }
            if !same {
                let m = "note: the restart file has another cell or face order; the fields were mapped \
                     by identity and the multigrid rebuilt (the continuation is not bit-identical)";
                lg.note(m.trim_start_matches("note: "));
                warnings.push(crate::diag::Diagnostic::new(
                    crate::diag::RUN_RESTART_MAPPED,
                    m,
                ));
            }
            if !time0.adaptive {
                dt = dt0;
            } else {
                dt = dt_next;
            }
            lg.note(&format!(
                "restart from {} at step {}, t = {:.6e} s",
                path.display(),
                clock.step,
                clock.time
            ));
            restart_from = Some(path.display().to_string());
        } else {
            // a steady restart file: its fields (velocity, pressure,
            // turbulence) are the start fields of the run at `start_time`,
            // with no older time levels, as after an initial condition
            fl.load(&r.state)?;
            crate::cht::load_restart(&mut fl, r.cht.as_ref(), path)?;
            af.update(&case, smesh, sgeo, &mut fl, time0.start_time)?;
            tv.update(&case, smesh, sgeo, &mut fl, time0.start_time)?;
            tr.start(&mut fl)?;
            lg.note(&format!(
                "start from the steady solution of {} (iteration {}) at t = {:.6e} s",
                path.display(),
                r.iteration,
                time0.start_time
            ));
        }
    } else {
        let (u, p) = match analytic0.filter(|a| a.is_transient() || case.initial.from_analytic) {
            Some(a) => {
                let mut u = Soa3::zeros(nc);
                let mut p = vec![0.0; nc];
                for (c, pc) in p.iter_mut().enumerate() {
                    let x = geo.cell_centroid[c];
                    u.set(c, a.velocity(x));
                    *pc = a.pressure(x);
                }
                (u, p)
            }
            None => {
                let u0 = case.initial.velocity;
                (
                    Soa3 {
                        x: vec![u0[0]; nc],
                        y: vec![u0[1]; nc],
                        z: vec![u0[2]; nc],
                    },
                    vec![
                        if reference_cell.is_some() {
                            case.pressure.reference_value
                        } else {
                            case.initial.pressure
                        };
                        nc
                    ],
                )
            }
        };
        // the start values of the regions and parts (D238)
        let analytic_start = analytic0
            .filter(|a| a.is_transient() || case.initial.from_analytic)
            .is_some();
        let (u, p) = match (!analytic_start)
            .then(|| {
                crate::start::initial_flow(
                    &case,
                    mesh,
                    geo,
                    case.initial.velocity,
                    p[0],
                    reference_cell.is_some(),
                )
            })
            .transpose()?
            .flatten()
        {
            Some((u, p, lines)) => {
                start_values = lines;
                (u, p)
            }
            None => (u, p),
        };
        start_values.extend(crate::start::temperature_report(&case, mesh, geo)?);
        start_values.extend(crate::start::turbulence_report(&case));
        if !start_values.is_empty() {
            for n in crate::start::notes(&start_values) {
                lg.note(n.as_str());
            }
        }
        af.update(&case, smesh, sgeo, &mut fl, time0.start_time)?;
        tv.update(&case, smesh, sgeo, &mut fl, time0.start_time)?;
        match &on {
            // computed on the whole mesh, each rank its cells (W5)
            Some(o) => {
                let u = Soa3 {
                    x: o.own(&u.x),
                    y: o.own(&u.y),
                    z: o.own(&u.z),
                };
                fl.initialise(&u, &o.own(&p))?;
            }
            None => fl.initialise(&u, &p)?,
        }
        if case.initial.potential_flow {
            let it = fl.potential_flow(smesh, sgeo)?;
            lg.note(&format!(
                "start from potential flow (D52): {it} CG iterations"
            ));
        }
        tr.start(&mut fl)?;
    }
    let setup = t_setup.elapsed().as_secs_f64();
    // the fields of the output: on ranks gathered into the global order
    // (collective), with the global discretisation
    macro_rules! out_sol {
        () => {
            match &on {
                Some(o) => o.solution(&mut fl)?,
                None => Solution {
                    state: fl.state()?,
                    disc: &fl.disc,
                    phys: fl.physics(),
                },
            }
        };
    }
    let first_step = tr.clock.step;

    // monitors
    let bench = case.verification.benchmark;
    let extra = benchmark_points(&case, mesh);
    let probes = Probes::new(&case, &extra, &loc, geo0);
    // conjugate heat transfer (stage 7d6): the probes of the output on the
    // full mesh, with the temperature in the solids (velocity and pressure
    // zero there, as in the fields)
    let full_probes = fl.cht.as_deref().map(|c| {
        let v = c.view.clone();
        let loc = Locator::new(&v.mesh, &v.geo);
        (Probes::new(&case, &extra, &loc, &v.geo), v)
    });
    let probe_temperature = case.energy.is_some();
    let cont = restart_from.as_ref().map(|_| first_step);
    let hist_dir = out;
    let mut hunits = String::from(
        "step -, time s, dt s, cfl - (largest Courant number), momentum_x momentum_y \
         momentum_z continuity - (scaled residuals), outer cg - (iterations)",
    );
    let mut hhead =
        String::from("step,time,dt,cfl,momentum_x,momentum_y,momentum_z,continuity,outer,cg");
    if !case.turbulence.is_laminar() {
        for s in crate::turbulence::equations(case.turbulence.model) {
            let _ = write!(hunits, ", {s} - (scaled residual)");
            let _ = write!(hhead, ",{s}");
        }
    }
    for s in case.scalars.keys() {
        let _ = write!(
            hunits,
            ", scalar_{s} - (scaled residual), first_order_{s} - (cells at first order)"
        );
        let _ = write!(hhead, ",scalar_{s},first_order_{s}");
    }
    if case.energy.is_some() {
        hunits
            .push_str(", energy - (scaled residual), first_order_energy - (cells at first order)");
        hhead.push_str(",energy,first_order_energy");
    }
    let vd_csv = crate::density::wanted(&case);
    if vd_csv {
        hunits.push_str(
            ", p0 Pa (thermodynamic pressure), mass kg (in the domain), rho_min rho_max kg/m³",
        );
        hhead.push_str(",p0,mass,rho_min,rho_max");
    }
    if deform.is_some() {
        hunits.push_str(", volume m³ (of the deforming mesh, stage MS0)");
        hhead.push_str(",volume");
    }
    let mut history = History::new(
        hist_dir,
        &suffix_name("history", desc, last),
        &crate::case::csv_preamble("history", name, None, &hunits),
        &hhead,
        cont,
    )?;
    // the heat balance of every step (stage 7d6): what the boundaries and
    // the sources bring in, what each region stores, what goes through
    // the interfaces of a conjugate case
    let heat_regions: Vec<String> = match fl.cht.as_deref() {
        Some(c) => c.fe.scalars.as_ref().map(|s| s.energy_storage_names()),
        None => fl.scalars.as_ref().map(|s| s.energy_storage_names()),
    }
    .unwrap_or_default();
    let heat_interfaces: Vec<String> = fl
        .cht
        .as_deref()
        .map(|c| c.view.interfaces.clone())
        .unwrap_or_default();
    let mut heat_relative_max: Option<f64> = None;
    let accel_csv = case
        .energy
        .as_ref()
        .is_some_and(|e| e.acceleration.is_some());
    let mut heat_csv = if last && !heat_regions.is_empty() {
        let mut head = String::from("step,time,conduction,convection,source,storage");
        for r in &heat_regions {
            let _ = write!(head, ",{r}_storage");
        }
        for i in &heat_interfaces {
            let _ = write!(head, ",{i}_heat");
        }
        if accel_csv {
            head.push_str(",inertia_scaling");
        }
        head.push_str(",imbalance,relative");
        Some(History::new(
            out,
            "heat.csv",
            &crate::case::csv_preamble(
                "heat",
                name,
                None,
                "step -, time s, conduction W (into the domain through the boundaries), \
                 convection W (sensible enthalpy carried in by the flow), source W (q''' and \
                 viscous heating), storage W (stored in the step, all regions), \
                 <region>_storage W, <interface>_heat W (into the fluid), inertia_scaling W \
                 (with [energy.acceleration]: the heat the true ρ c_p of the accelerated solids \
                 would have stored in addition), imbalance W \
                 (conduction + convection + source - storage), relative - (|imbalance| over \
                 the largest heat flow)",
            ),
            &head,
            cont,
        )?)
    } else {
        None
    };
    let walls: Vec<String> = mesh
        .boundaries
        .iter()
        .filter(|p| {
            matches!(
                case.boundary.get(&p.name),
                Some(BoundaryCondition::Wall { .. })
            )
        })
        .map(|p| p.name.clone())
        .collect();
    let mut fh = String::from("step,time");
    // rotating machinery (D184): torque per wall and per rotor, shaft power,
    // the machines
    let rotor_names: Vec<String> = crate::rotor::rotors(&case, mesh)
        .into_iter()
        .map(|r| r.name)
        .collect();
    for w in &walls {
        let _ = write!(fh, ",{w}_fx,{w}_fy,{w}_fz");
        if !rotor_names.is_empty() {
            let _ = write!(fh, ",{w}_torque");
        }
    }
    for r in &rotor_names {
        let _ = write!(fh, ",{r}_torque,{r}_power");
        if case
            .sliding
            .iter()
            .any(|m| &crate::rotor::rotor_name(&m.zone) == r)
        {
            let _ = write!(fh, ",{r}_angle");
        }
    }
    for m in &case.output.machines {
        let n = &m.name;
        let _ = write!(fh, ",{n}_flow,{n}_dpt,{n}_head,{n}_power,{n}_efficiency");
    }
    if case.output.coefficients.is_some() {
        fh.push_str(",cd,cl");
    }
    let monitors = last && (!walls.is_empty() || !probes.names.is_empty());
    let mut forces = if monitors {
        Some(History::new(
            out,
            "forces.csv",
            &crate::case::csv_preamble(
                "forces",
                name,
                None,
                "step -, time s, <wall>_fx <wall>_fy <wall>_fz N (force of the fluid on the wall), <wall>_torque N m (about the axis of its rotor), <rotor>_torque N m, <rotor>_power W (shaft power), <rotor>_angle deg, <machine>_flow m3/s, _dpt Pa, _head m, _power W, _efficiency -, cd cl -",
            ),
            &fh,
            cont,
        )?)
    } else {
        None
    };
    // rotating machinery (D183): the velocity relative to the turning zone
    let frames = crate::derived::Frames::new(&case);
    let mut ph = String::from("step,time");
    for n in &probes.names {
        let _ = write!(ph, ",{n}_u,{n}_v,{n}_w,{n}_p");
        if frames.is_some() {
            let _ = write!(ph, ",{n}_u_rel,{n}_v_rel,{n}_w_rel");
        }
        if probe_temperature {
            let _ = write!(ph, ",{n}_T");
        }
    }
    let mut probe_units =
        String::from("step -, time s, <probe>_u <probe>_v <probe>_w m/s, <probe>_p Pa");
    if frames.is_some() {
        probe_units.push_str(
            ", <probe>_u_rel <probe>_v_rel <probe>_w_rel m/s (relative to the turning zone)",
        );
    }
    if probe_temperature {
        probe_units.push_str(", <probe>_T K (the value of the cell)");
    }
    let mut probe_csv = if monitors && !probes.names.is_empty() {
        Some(History::new(
            out,
            "probes.csv",
            &crate::case::csv_preamble("probes", name, None, &probe_units),
            &ph,
            cont,
        )?)
    } else {
        None
    };
    // benchmark histories
    let (mut bt, mut bcd, mut bcl, mut bdp) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    // the conservation ledger (D88, D89, D90 part 4), booked every step
    // from the state at the start or the restart
    let nonorth = case.solver.nonorth_correction;
    let (bunits, bhead) = balance::csv_head(true, !case.turning().is_empty());
    let mut balance_csv = History::new(
        out,
        &suffix_name("balance", desc, last),
        &crate::case::csv_preamble("balance", name, None, &bunits),
        &bhead,
        cont,
    )?;
    let mut spray_csv = if spray.is_some() {
        Some(History::new(
            out,
            &suffix_name("spray", desc, last),
            &crate::case::csv_preamble("spray", name, None, crate::spray::CSV_UNITS),
            &crate::spray::csv_head(spray.as_ref().context("spray")?),
            cont,
        )?)
    } else {
        None
    };
    // two fluids by the volume of fluid (L3): `vof.csv`
    let mut vof_csv = if case.vof.is_empty() {
        None
    } else {
        Some(History::new(
            out,
            &suffix_name("vof", desc, last),
            &crate::case::csv_preamble("vof", name, None, crate::vof::CSV_UNITS),
            crate::vof::CSV_HEAD,
            cont,
        )?)
    };
    let mut chem_csv = match (&chem, &case.species) {
        (Some(_), Some(sp)) => {
            let names = crate::species::Resolved::of(sp)?.names;
            Some(History::new(
                out,
                &suffix_name("chemistry", desc, last),
                &crate::case::csv_preamble(
                    "chemistry",
                    name,
                    None,
                    crate::chemistry::flow::Chemistry::CSV_UNITS,
                ),
                &crate::chemistry::flow::Chemistry::csv_head(&names),
                cont,
            )?)
        }
        _ => None,
    };
    let mut spray_lost_warned = false;
    // the two-way coupling over the run and of the last step (SP2)
    let mut spray_coupling: Option<crate::spray::CouplingSummary> = None;
    let mut spray_step: Option<crate::spray::StepCoupling> = None;
    let mut coupling_warned = false;
    let mut limiter_warned = false;
    let mut reports = (last && !datum_built.is_empty())
        .then(|| crate::datum::Reports::new(datum_built, out, name, cont));
    let mut series = last.then(|| crate::series::Series::new(out, name, cont));
    // the statistics of D71 (the final solve of a sweep only)
    let stat_cfg = case.output.statistics.clone().filter(|_| last);
    let stat_start = stat_cfg
        .as_ref()
        .map(|c| c.start.unwrap_or(time0.start_time).max(tr.clock.time));
    if let (Some(c), Some(s0)) = (&stat_cfg, stat_start)
        && restart_from.is_some()
        && c.start.unwrap_or(time0.start_time) < s0
    {
        lg.note(&format!(
            "statistics (D71) start at the restart time {s0:.6e} s: the averages \
             of the earlier run are not continued"
        ));
    }
    let mut stat_series = stat_cfg.as_ref().map(|_| {
        let mut names = Vec::new();
        for n in &probes.names {
            for c in ["u", "v", "w", "p"] {
                names.push(format!("probes.{n}.{c}"));
            }
        }
        if case.output.coefficients.is_some() {
            names.push("coefficients.cd".into());
            names.push("coefficients.cl".into());
        }
        crate::statistics::Series::new(names)
    });
    let mut stat_fields = stat_cfg.as_ref().filter(|c| c.fields).map(|_| {
        let mut names: Vec<String> = ["velocity_x", "velocity_y", "velocity_z", "pressure"]
            .map(String::from)
            .into();
        let mut pairs = vec![
            (0, 1, "covariance_velocity_xy".to_string()),
            (1, 2, "covariance_velocity_yz".to_string()),
            (0, 2, "covariance_velocity_xz".to_string()),
        ];
        if !case.turbulence.is_laminar() {
            names.push("k".into());
        }
        for sname in case.scalars.keys() {
            let i = names.len();
            names.push(sname.clone());
            for (j, c) in ["x", "y", "z"].into_iter().enumerate() {
                pairs.push((j, i, format!("covariance_velocity_{c}_{sname}")));
            }
        }
        crate::statistics::FieldMoments::new(names, pairs, mesh.n_cells)
    });
    let mut ledger = {
        let sol = out_sol!();
        let mut r = balance::rates(&sol, mesh, geo, &case, &mom_kind, &geo.cell_volume, nonorth);
        if let Some(sp) = spray.as_ref().filter(|s| s.config().coupling.two_way()) {
            let gm = crate::spray::gas_momentum(&sol, &geo.cell_volume);
            crate::spray::book(&mut r, sp, None, dt, gm, case.energy.is_some());
        }
        Ledger::start(r, tr.clock.time, first_step, case.balance.tolerance_energy)
    };
    let mut last_balance = None;
    let mut first_cfl = None;
    let mut balance_warned = false;

    let ts = Instant::now();
    let mut cols = crate::logcols::Columns::build(&case, mesh, geo, &loc, true)?;
    let mut table = cols.table("Step", case.log.style, case.log.header_every);
    if steer.record.on() {
        steer.record.row(&cols.record_header("step", "Step"));
    }
    lg.frame(&crate::console::Frame::new(format!(
        "Time steps ({}, {})",
        crate::logblocks::spelled(&case.solver.algorithm(true)).to_uppercase(),
        case.time.unwrap_or_default().scheme.label()
    )));
    lg.line("");
    if !case.log.detail.rows() {
        lg.note("the table is not shown ([log] detail = \"summary\"); run.jsonl holds every row");
    }
    lg.sync(&mut steer.record);
    let mut stopped = false;
    let (mut dt_min, mut dt_max, mut cfl_max) = (f64::INFINITY, 0.0f64, 0.0f64);
    let (mut outer_sum, mut cg_sum, mut nsteps) = (0usize, 0usize, 0usize);
    let mut rejected = 0usize;
    // the largest closure of the space conservation law of a step (MS0)
    let mut gcl_worst = 0.0_f64;
    let mut last_rep = StepReport::default();
    let mut sched = crate::schedule::Runner::default();
    // the thermal acceleration of the solids (D213)
    let mut accel = match case.energy.as_ref().and_then(|e| e.acceleration.clone()) {
        Some(a) => crate::acceleration::Controller::new(&mut fl, &a, tr.clock.time)?,
        None => None,
    };
    // the cycle-averaged coupling of a periodic flow with its solids (D204)
    let mut cycle = match case.energy.as_ref().and_then(|e| e.cycle.clone()) {
        Some(c) => Some(crate::cycle::Controller::new(
            &mut fl,
            &c,
            tr.clock.time,
            dt,
            crate::case::dimension(&case),
            case.time.as_ref().and_then(|t| t.end_time),
        )?),
        None => None,
    };
    match (&mut cycle, &cycle_state) {
        (Some(c), Some(b)) => c.restore(&mut fl, b, tr.clock.time + dt)?,
        (None, Some(_)) => lg.note(
            "the restart file holds the state of a cycle-averaged loop; the case runs without \
             [energy.cycle], so the run continues fully coupled",
        ),
        _ => {}
    }
    // the frozen-flow march of a transient conjugate run (D118)
    let mut frozen = match case.energy.as_ref().and_then(|e| e.frozen_flow.clone()) {
        Some(f) => {
            let c = crate::frozen::Controller::new(
                &mut fl,
                &f,
                case.buoyancy.model != crate::buoyancy::Model::None,
                tv.boundary || tv.source || !af.faces.is_empty(),
                tr.clock.time,
            )?;
            for l in c.start_lines() {
                lg.line(&l);
            }
            Some(c)
        }
        None => None,
    };
    match (&mut frozen, &frozen_state) {
        (Some(c), Some(b)) => c.restore(b)?,
        (None, Some(_)) => lg.note(
            "the restart file holds the state of a frozen-flow march; the case runs without \
             [energy.frozen_flow], so the run continues fully coupled",
        ),
        _ => {}
    }
    loop {
        // [[schedule]] (D128): the entries due before this step
        if !case.schedule.is_empty() {
            let restarted = case.run.restart.is_some();
            let p = crate::schedule::Position {
                n: tr.clock.step + 1,
                time: tr.clock.time,
                tol: 1e-6 * dt,
            };
            if let Some(old) =
                crate::case::schedule_step(&mut case, &mut fl, lg, steer, &mut sched, p, restarted)?
                && let (Some(a), Some(b)) = (old.time, case.time)
                && a.time_step != b.time_step
                && !b.adaptive
            {
                dt = b.time_step;
            }
        }
        let t = case.time.context("[time]")?;
        let step = tr.clock.step;
        if t.steps.is_some_and(|n| step >= n) {
            break;
        }
        // the frozen-flow march (D118): an energy-only step of its own size
        let frozen_step = frozen.as_ref().is_some_and(|f| f.frozen());
        let mut h = match &frozen {
            Some(f) if frozen_step => f.time_step(),
            _ => dt,
        };
        // D66: land on the next scheduled output time; when it lies within
        // two steps, in two equal steps, so that no step gets shorter than
        // half the planned one (not with a sliding mesh, whose step is the
        // time that the rotor needs for its sectors: the outputs of a
        // schedule in time come with the step that reaches them)
        let t0 = tr.clock.time;
        let tol = 1e-6 * h;
        if slide.is_none()
            && let Some(te) = schedules(&case)
                .iter()
                .filter_map(|s| s.next(t0, tol))
                .reduce(f64::min)
        {
            let rest = te - t0;
            if rest <= h + tol {
                h = rest;
            } else if rest < 2.0 * h {
                h = 0.5 * rest;
            }
        }
        // D128: land on the time of the next scheduled change as well
        if slide.is_none()
            && let Some(te) = sched.next_time(&case, t0, tol)
        {
            let rest = te - t0;
            if rest <= h + tol {
                h = rest;
            } else if rest < 2.0 * h {
                h = 0.5 * rest;
            }
        }
        if let Some(end) = t.end_time {
            let rest = end - tr.clock.time;
            if slide.is_some() {
                // whole steps (checked when the sliding mesh was set up)
                if rest < 0.5 * dt {
                    break;
                }
            } else {
                if rest <= 1e-9 * dt {
                    break;
                }
                if h > rest || rest - h < 1e-6 * h {
                    h = rest;
                }
            }
        }
        // D182 variant B: the mesh of the end of the step, the solver built
        // on it
        if let Some((imp, _)) = &imprint {
            let next = imp.working(mesh0, geo0, t0 + h)?;
            let prev = work.as_ref().context("the mesh of the rotor angle")?;
            let map = crate::sliding::face_map(prev, &next);
            let d = turning_distance.as_mut().map(|(twd, _)| {
                twd.at_angles(&next.angles, &next.geo.cell_centroid[..next.mesh.n_cells])
            });
            let (nfl, ntr) = match on.as_ref() {
                None => rebuild_flow(
                    &case,
                    &tv,
                    &target,
                    reference_cell,
                    (&next.mesh, &next.geo),
                    Some(&map),
                    t0 + h,
                    &mut fl,
                    &tr,
                    d,
                )?,
                Some(o) => {
                    let (no, nfl, ntr) = rebuild_on_ranks(
                        &case,
                        &tv,
                        &target,
                        reference_cell,
                        &next,
                        &map,
                        t0 + h,
                        &mut fl,
                        &tr,
                        d,
                        o,
                    )?;
                    on = Some(no);
                    (nfl, ntr)
                }
            };
            setups_past.0 += fl.amg_setups;
            setups_past.1 += fl.sweep_setups;
            fl = nfl;
            tr = ntr;
            work = Some(next);
        }
        // D91 (stage MS0): the deforming mesh at the end of the step, the
        // solver built on it, the grid of the step attached
        if let Some(d) = deform.as_mut() {
            let grid = d.advance(t0 + h)?;
            let (nfl, ntr) = rebuild_flow(
                &case,
                &tv,
                &target,
                reference_cell,
                (&d.mesh, &d.geo),
                None,
                t0 + h,
                &mut fl,
                &tr,
                None,
            )?;
            setups_past.0 += fl.amg_setups;
            setups_past.1 += fl.sweep_setups;
            fl = nfl;
            tr = ntr;
            gcl_worst = gcl_worst.max(d.closure(&grid));
            fl.attach_grid(grid)?;
        }
        let (mesh, geo) = current(&work, &deform, mesh0, geo0);
        let (smesh, sgeo) = match &on {
            Some(o) => (&o.local.mesh, &o.geo),
            None => (mesh, geo),
        };
        // D65: the error estimate, with the tolerances of now (re-read)
        tr.error_control = t
            .error_tolerance
            .filter(|_| t.adaptive)
            .map(|r| ErrorControl {
                relative: r,
                absolute: t.error_absolute,
            });
        let min_dt = t.min_time_step;
        let mut tries = 0;
        let mut spray_passes = None;
        let stepped = loop {
            af.update(&case, smesh, sgeo, &mut fl, t0 + h)?;
            tv.update(&case, smesh, sgeo, &mut fl, t0 + h)?;
            fl.keep_valid()?;
            let two_way = spray.as_ref().filter(|s| s.config().coupling.two_way());
            let r = match &frozen {
                Some(f) if frozen_step => {
                    tr.energy_step(&mut fl, h, f.restart_levels(), f.iterations())
                }
                // the parcel passes inside the step (SP2, D252)
                _ if two_way.is_some() => {
                    let sp = two_way.context("spray")?;
                    let env = foehncraft_particles::Env {
                        mesh,
                        geo,
                        loc: &loc,
                    };
                    let mut hk = crate::spray::TwoWay::new(
                        sp,
                        env,
                        geo,
                        t0,
                        h,
                        crate::spray::VapourLink::of(&case),
                    );
                    let r = tr.step_with(&mut fl, h, Some(&mut hk));
                    fl.spray_momentum = None;
                    spray_passes = Some(hk.finish());
                    r
                }
                // the reaction step of the cells (SF2b)
                _ if chem.is_some() => {
                    let ch = chem.as_mut().context("chemistry")?;
                    let mut hk = crate::chemistry::flow::Step::new(ch, h, t0 + h);
                    tr.step_with(&mut fl, h, Some(&mut hk))
                }
                _ => tr.step(&mut fl, h),
            };
            // a step whose error exceeds the tolerance is taken back and
            // repeated with a shorter one (at most MAX_REJECTIONS times, and
            // not below min_time_step); with the exponent ½ rather than ⅓,
            // since the estimate falls only as Δt² when Δt gets small
            // against the step before (the error constant of variable-step
            // BDF2 grows as 1/ω)
            if let Ok(rp) = &r
                && let Some(e) = rp.error
                && e > 1.0
                && tries < MAX_REJECTIONS
                && h > min_dt
            {
                tries += 1;
                rejected += 1;
                let hn = (h * (0.9 * e.powf(-0.5)).max(0.2)).max(min_dt);
                lg.note(&format!(
                    "step {} rejected: error {e:.2} of the tolerance at Δt {h:.3e} s, \
                     repeated with Δt {hn:.3e} s",
                    rp.step
                ));
                af.update(&case, smesh, sgeo, &mut fl, t0)?;
                tv.update(&case, smesh, sgeo, &mut fl, t0)?;
                tr.reject(&mut fl)?;
                h = hn;
                continue;
            }
            break r;
        };
        let rp = match stepped {
            Ok(r) => r,
            Err(e) if crate::divergence::is_divergence(&e, &mut fl) => {
                // D63: the last valid state, the snapshot and the report
                lg.line(&crate::console::tagged("error", &format!("{e:#}")));
                let ctx = crate::divergence::Context {
                    name,
                    case: &case,
                    at: tr.clock.step + 1,
                    done: nsteps,
                    time: Some(tr.clock.time),
                    dt: Some(h),
                    residuals: fl.diverged.clone().unwrap_or_default(),
                    last_residuals: (nsteps > 0).then_some(last_rep.residuals.clone()),
                    first_courant: first_cfl,
                };
                let prepared = (|| -> Result<(crate::simplec::State, Vec<&'static str>)> {
                    let failed = out_sol!().state;
                    anyhow::ensure!(tr.rewind(&mut fl)?, "no valid state was kept");
                    let mut files = Vec::new();
                    if last && cfg!(feature = "hdf5") && case.output.vtkhdf {
                        let path = out.join("fields_last_valid.vtkhdf");
                        let at = Some((tr.clock.step, Some(tr.clock.time)));
                        if case.fluid_view.is_some() {
                            let st = fl.state()?;
                            let full = crate::cht::full_state(&mut fl, &st)?;
                            write_full_fields(&fl, &full, &path, at)?;
                        } else {
                            let sol = out_sol!();
                            if lead {
                                write_fields(
                                    &case,
                                    &path,
                                    mesh,
                                    geo,
                                    &sol,
                                    case.analytic_at(tr.clock.time).as_ref(),
                                    at,
                                )?;
                            }
                        }
                        files.push("fields_last_valid.vtkhdf");
                    }
                    if last && restart::AVAILABLE && !moving {
                        let mut rs = restart_of(
                            &mut fl,
                            &tr,
                            smesh,
                            name,
                            h,
                            cycle.as_ref(),
                            frozen.as_ref(),
                            on.as_ref(),
                        )?;
                        rs.spray = crate::spray::state(spray.as_ref())?;
                        if lead {
                            restart::write(&out.join("restart.h5"), &rs)?;
                        }
                        files.push("restart.h5");
                    }
                    Ok((failed, files))
                })();
                let e = match &on {
                    None => crate::divergence::finish(
                        &mut fl,
                        &ctx,
                        mesh,
                        geo,
                        out,
                        prepared,
                        case.log.style,
                        &mut |l| lg.line(l),
                    ),
                    // on ranks the last valid state gathered (W7)
                    Some(o) => {
                        let valid = match &prepared {
                            Ok(_) => o.solution(&mut fl).ok(),
                            Err(_) => None,
                        };
                        crate::divergence::finish_from(
                            valid.as_ref(),
                            &mom_kind,
                            &ctx,
                            mesh,
                            geo,
                            out,
                            prepared,
                            case.log.style,
                            &mut |l| lg.line(l),
                        )
                    }
                };
                history.flush()?;
                balance_csv.flush()?;
                if let Some(f) = heat_csv.as_mut() {
                    f.flush()?;
                }
                lg.sync(&mut steer.record);
                steer.record.flush();
                std::fs::write(out.join("flow.log"), &lg.text)?;
                return Err(e);
            }
            Err(e) => return Err(e),
        };
        last_rep = rp.clone();
        if nsteps == 0 {
            first_cfl = Some(rp.cfl);
        }
        nsteps += 1;
        outer_sum += rp.outer;
        cg_sum += rp.residuals.pressure_iterations;
        dt_min = dt_min.min(h);
        dt_max = dt_max.max(h);
        cfl_max = cfl_max.max(rp.cfl);
        let n = rp.step;
        let r = &rp.residuals;
        // the heat balance of the step (stage 7d6)
        let heat = match heat_csv.as_mut() {
            Some(hc) if n.is_multiple_of(case.output.history_every) => {
                let s = crate::scalar::summaries_on(
                    &mut fl,
                    &case,
                    mesh,
                    geo,
                    crate::case::dimension(&case),
                    rp.time,
                    on.as_ref(),
                )?;
                let t = s.into_iter().find(|s| s.name == crate::energy::FIELD);
                let hb = t.as_ref().and_then(|t| t.energy.clone()).map(|e| StepHeat {
                    interfaces: heat_interfaces
                        .iter()
                        .map(|i| {
                            t.iter()
                                .flat_map(|t| &t.boundaries)
                                .find(|b| &b.boundary == i)
                                .and_then(|b| b.heat_flow)
                                .unwrap_or(f64::NAN)
                        })
                        .collect(),
                    balance: e,
                });
                if let Some(h) = &hb {
                    let e = &h.balance;
                    let mut row = format!(
                        "{n},{:e},{:e},{:e},{:e},{:e}",
                        rp.time,
                        e.conduction,
                        e.convection,
                        e.source,
                        e.storage.unwrap_or(0.0)
                    );
                    for r in &heat_regions {
                        let v = e
                            .storage_regions
                            .iter()
                            .find(|x| &x.region == r)
                            .map_or(0.0, |x| x.storage);
                        let _ = write!(row, ",{v:e}");
                    }
                    for q in &h.interfaces {
                        let _ = write!(row, ",{q:e}");
                    }
                    if accel_csv {
                        let _ = write!(row, ",{:e}", e.inertia_scaling.unwrap_or(0.0));
                    }
                    let _ = write!(row, ",{:e},{:e}", e.imbalance, e.relative);
                    hc.row(&row);
                    heat_relative_max = Some(heat_relative_max.unwrap_or(0.0).max(e.relative));
                }
                hb
            }
            _ => None,
        };
        if let Some(acc) = accel.as_mut() {
            for l in acc.after_step(&mut fl, rp.time, rp.time - t0)? {
                lg.line(&l);
                steer.event(crate::diag::ACC_PHASE, n, &l);
            }
        }
        if let Some(fz) = frozen.as_mut() {
            for l in fz.after_step(&mut fl, &rp, frozen_step)? {
                lg.line(&l);
                steer.event(crate::diag::FRZ_PHASE, n, &l);
            }
        }
        let mut cycle_done = false;
        if let Some(cy) = cycle.as_mut() {
            for l in cy.after_step(&mut fl, rp.time, rp.time - t0)? {
                lg.line(&l);
                steer.event(crate::diag::CYC_PHASE, n, &l);
            }
            // a transient solid phase: the solids have marched past the
            // end of the fluid phase, the clock follows (whole periods)
            if let Some(j) = cy.take_jump() {
                tr.clock.time += j;
            }
            cycle_done = cy.done();
        }
        if n.is_multiple_of(case.output.history_every) {
            let mut sc = String::new();
            for x in &r.turbulence.residual {
                let _ = write!(sc, ",{x:e}");
            }
            for (x, k) in r.scalars.residual.iter().zip(&r.scalars.first_order) {
                let _ = write!(sc, ",{x:e},{k}");
            }
            if let (true, Some(d)) = (vd_csv, fl.density.as_deref()) {
                let (lo, hi) = d
                    .rho
                    .iter()
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &r| {
                        (a.min(r), b.max(r))
                    });
                let _ = write!(sc, ",{:e},{:e},{lo:e},{hi:e}", d.p0, d.total_mass());
            }
            if let Some(d) = &deform {
                let _ = write!(sc, ",{:e}", d.volume());
            }
            history.row(&format!(
                "{n},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{},{}{sc}",
                rp.time,
                h,
                rp.cfl,
                r.momentum[0],
                r.momentum[1],
                r.momentum[2],
                r.continuity,
                rp.outer,
                r.pressure_iterations
            ));
            if steer.record.on() {
                let mut row = serde_json::json!({
                    "type": "step",
                    "step": n,
                    "time": rp.time,
                    "dt": h,
                    "cfl": rp.cfl,
                    "momentum": r.momentum,
                    "continuity": r.continuity,
                    "outer": rp.outer,
                    "cg": r.pressure_iterations,
                    "scalars": r.scalars.residual,
                    "turbulence": r.turbulence.residual,
                    "turbulence_floors": r.turbulence.floors,
                });
                if let Some(e) = rp.error {
                    row["error"] = e.into();
                }
                if let Some(hb) = &heat {
                    row["heat"] = serde_json::to_value(hb)?;
                }
                steer.record.row(&row);
            }
        }
        if let (Some(v), Some(csv)) = (fl.vof.as_deref(), vof_csv.as_mut()) {
            csv.row(&v.csv_row(n, tr.clock.time, &geo.cell_centroid));
        }
        if let (Some(ch), Some(csv)) = (chem.as_ref(), chem_csv.as_mut())
            && let Some(row) = ch.csv_row(n)
        {
            csv.row(&row);
        }
        if let (Some(sp), Some(csv)) = (spray.as_mut(), spray_csv.as_mut()) {
            let sol = out_sol!();
            let (lost, file) = match spray_passes.take() {
                Some(done) => {
                    drop(sol);
                    let (lost, file, c) = crate::spray::step_coupled(
                        &case, sp, &mut fl, done, mesh, geo, &loc, n, out, name, csv,
                    )?;
                    spray_coupling
                        .get_or_insert_with(Default::default)
                        .add(&c, h);
                    if c.unconverged && !coupling_warned {
                        coupling_warned = true;
                        let m = format!(
                            "step {n}: the spray passes ended at [spray.coupling] max_passes with \
                             the coupling lag {:.2e} above the tolerance (more passes or a smaller \
                             step; `spray.csv` counts such steps)",
                            c.lag
                        );
                        lg.warn(&format!("{}: {m}", crate::diag::SPRAY_COUPLING));
                        steer.event(crate::diag::SPRAY_COUPLING, n, &m);
                    }
                    if c.limiter_cells > 0 && !limiter_warned {
                        limiter_warned = true;
                        let m = format!(
                            "step {n}: the saturation limiter (D253) held back {:.3e} kg/s of \
                             vapour in {} spray cells (more passes or a smaller step; \
                             `spray.csv` counts every step)",
                            c.limiter_mass, c.limiter_cells
                        );
                        lg.warn(&format!("{}: {m}", crate::diag::SPRAY_LIMITER));
                        steer.event(crate::diag::SPRAY_LIMITER, n, &m);
                    }
                    spray_step = Some(c);
                    (lost, file)
                }
                None => {
                    let heat = match fl.scalars.as_ref() {
                        Some(sc) => sc
                            .energy_gas(fl.dev.as_mut())?
                            .map(|(t, cp, k)| crate::spray::GasHeat { t, cp, k }),
                        None => None,
                    };
                    // the vapour the parcels evaporate into (SP4)
                    let vapour = match (crate::spray::VapourLink::of(&case), fl.scalars.as_ref()) {
                        (Some(l), Some(sc)) => crate::spray::GasVapour::from_parts(
                            sc,
                            fl.dev.as_mut(),
                            &sol.phys,
                            fl.density.as_deref(),
                            &l,
                            heat.as_ref(),
                        )?,
                        _ => None,
                    };
                    crate::spray::step(
                        &case,
                        sp,
                        &sol,
                        mesh,
                        geo,
                        &loc,
                        n,
                        rp.time,
                        h,
                        out,
                        name,
                        csv,
                        heat.as_ref(),
                        vapour.as_ref(),
                    )?
                }
            };
            if let (Some(f), Some(se)) = (file, series.as_mut()) {
                se.add(n, rp.time, "parcels", &f)?;
            }
            if lost > 0 && !spray_lost_warned {
                spray_lost_warned = true;
                lg.warn(&format!(
                    "{}: step {n}: {lost} parcel(s) lost (the walk failed or [spray] max_substeps \
                     was exceeded); `spray.csv` counts them",
                    crate::diag::SPRAY_LOST
                ));
            }
        }
        let bal = {
            let sol = out_sol!();
            let r = balance::rates(&sol, mesh, geo, &case, &mom_kind, &geo.cell_volume, nonorth);
            let mut r = r;
            if let Some(sp) = spray.as_ref().filter(|s| s.config().coupling.two_way()) {
                let gm = crate::spray::gas_momentum(&sol, &geo.cell_volume);
                crate::spray::book(
                    &mut r,
                    sp,
                    spray_step.as_ref(),
                    h,
                    gm,
                    case.energy.is_some(),
                );
            }
            // the tolerance may change during the run (D40)
            ledger.tolerance = case.balance.tolerance_energy;
            ledger.step(r, n, rp.time)
        };
        if let Some(s0) = stat_start
            && rp.time - 0.5 * h >= s0
        {
            let sol = out_sol!();
            if let Some(m) = stat_fields.as_mut() {
                let mut x: Vec<&[f64]> = vec![&sol.u.x, &sol.u.y, &sol.u.z, &sol.p];
                if let Some(t) = &sol.state.turbulence {
                    x.push(&t.k);
                }
                for st in &sol.state.scalars {
                    x.push(&st.phi);
                }
                m.add(h, &x);
            }
            if let Some(ser) = stat_series.as_mut() {
                let here = moving.then(|| probes.relocated(mesh, geo));
                let pr = here.as_ref().unwrap_or(&probes);
                let mut row = Vec::with_capacity(ser.names.len());
                for (u, p) in pr.values(&sol, geo) {
                    row.extend([u[0], u[1], u[2], p]);
                }
                if case.output.coefficients.is_some() {
                    let b = boundary_summaries(&sol, mesh, geo, &case);
                    let [cd, cl] = monitor::coefficients(&case, &b).unwrap_or([f64::NAN; 2]);
                    row.extend([cd, cl]);
                }
                ser.add(rp.time, h, row);
            }
        }
        if n.is_multiple_of(case.output.history_every) {
            balance_csv.row(&balance::csv_row(&bal));
            if steer.record.on() {
                steer.record.row(&balance::record_row(&bal));
            }
        }
        if bal.mass.relative.abs() > case.balance.tolerance_step && !balance_warned {
            balance_warned = true;
            let m = format!(
                "step {n}: the relative mass imbalance {:.2e} is above [balance] tolerance_step {:.1e}",
                bal.mass.relative, case.balance.tolerance_step
            );
            lg.warn(&format!("{}: {m}", crate::diag::BAL_MASS_STEP));
            steer.event(crate::diag::BAL_MASS_STEP, n, &m);
            if case.balance.action == BalanceAction::Stop {
                steer.stop_for(crate::diag::BAL_MASS_STEP, format!("stopped: {m}"));
            }
        }
        let log_now = n == first_step + 1 || n.is_multiple_of(case.solver.log_every);
        if log_now {
            lg.sync(&mut steer.record);
            steer.record.flush();
            // `log.columns` changed during the run (D214): new columns, the
            // header row again
            crate::logcols::refresh(
                &mut cols,
                &mut table,
                &mut crate::logcols::Refresh {
                    case: &case,
                    mesh,
                    geo,
                    loc: &loc,
                    transient: true,
                    table_name: "step",
                    first: "Step",
                    at: n,
                },
                &mut |m| lg.warn(m),
                steer,
            );
            let sol = if cols.needs_state() {
                Some(out_sol!())
            } else {
                None
            };
            let rotors = (cols.needs_rotors() && sol.is_some())
                .then(|| crate::rotor::evaluate(sol.as_ref().unwrap(), mesh, geo, &case));
            let entries = cols.eval(
                &crate::logcols::RowCtx {
                    index: n,
                    res: r,
                    time: Some(crate::logcols::TimeInfo {
                        time: rp.time,
                        dt: h,
                        courant: rp.cfl,
                        outer: rp.outer,
                        error: rp.error,
                    }),
                    sol: sol.as_ref(),
                    rotors: rotors.as_ref(),
                },
                geo,
            );
            let cells: Vec<String> = entries.iter().map(|e| e.text.clone()).collect();
            let detail = case.log.detail;
            // the box of the time step, then the table row
            if detail.rows() {
                lg.line("");
                let tb = crate::console::Frame::new(format!(
                    "{n:>7}. TIME STEP     TIME = {}     DT = {}",
                    crate::console::sci(rp.time, 5),
                    crate::console::sci(h, 5)
                ))
                .plain()
                .text(format!(
                    "{}COURANT NUMBER = {:.3}     OUTER ITERATIONS = {}",
                    " ".repeat(9),
                    rp.cfl,
                    rp.outer
                ));
                for l in tb.render(case.log.style) {
                    lg.line(&l);
                }
                lg.line("");
            }
            // the residuals of every outer iteration of this step
            if detail.outer() {
                let mut outer_table = crate::console::Table::new(
                    crate::console::Col::int("Outer", 6),
                    vec![
                        crate::console::Group::new(
                            "scaled residuals",
                            ["ResU", "ResV", "ResW", "ResM"]
                                .map(crate::console::Col::num)
                                .into(),
                        ),
                        crate::console::Group::new(
                            "linear iterations",
                            vec![
                                crate::console::Col::int("U", 4),
                                crate::console::Col::int("V", 4),
                                crate::console::Col::int("W", 4),
                                crate::console::Col::int("P", 4),
                            ],
                        ),
                    ],
                    case.log.style,
                    0,
                );
                for (k, o) in tr.outer_log.iter().enumerate() {
                    let cells: Vec<String> = o
                        .momentum
                        .iter()
                        .chain([&o.continuity])
                        .map(|v| crate::console::sci3(*v))
                        .chain(o.momentum_iterations.iter().map(ToString::to_string))
                        .chain([o.pressure_iterations.to_string()])
                        .collect();
                    for l in outer_table.row(k + 1, &cells) {
                        lg.line(&l);
                    }
                }
                lg.line("");
            }
            if detail.rows() {
                for l in table.row(n, &cells) {
                    lg.line(&l);
                }
            }
            if steer.record.on() {
                steer.record.row(&cols.record_row("step", n, &entries));
            }
            if detail.rows() {
                lg.line(&balance::compact(&bal));
            }
            if detail.rows() && case.balance.every > 0 && n.is_multiple_of(case.balance.every) {
                lg.frame(&crate::logblocks::balance_block(
                    &format!("Balance of the run up to step {n}"),
                    &ledger.run(),
                    true,
                ));
            }
            // the histories reach the disk with every log line, so a run
            // that is killed keeps them up to there
            history.flush()?;
            balance_csv.flush()?;
            if let Some(f) = heat_csv.as_mut() {
                f.flush()?;
            }
            if let Some(f) = forces.as_mut() {
                f.flush()?;
            }
            if let Some(f) = probe_csv.as_mut() {
                f.flush()?;
            }
            if let Some(r) = reports.as_mut() {
                r.flush()?;
            }
        }
        last_balance = Some(bal);
        let mut req = control.check(n, case.run.reread_every);
        // on ranks rank 0 reads the requests and decides for all
        let asked = steer.restart_requested();
        let asked = if let Some(o) = &on {
            let flags = o.agree(
                u32::from(req.result)
                    | u32::from(req.reread) << 1
                    | u32::from(req.stop) << 2
                    | u32::from(asked) << 3,
            );
            req.result = flags & 1 != 0;
            req.reread = flags & 2 != 0;
            req.stop = flags & 4 != 0;
            flags & 8 != 0
        } else {
            asked
        };
        let [sf, sm, sr] = schedules(&case);
        let tol = 1e-6 * h;
        let at_monitor = monitors
            && if sm.interval.is_some() {
                sm.fired(t0, rp.time, tol)
            } else {
                n.is_multiple_of(case.output.monitor_every)
            };
        // the datum reports on the schedule of the monitors (D262)
        let datum_now = if sm.interval.is_some() {
            sm.fired(t0, rp.time, tol)
        } else {
            n.is_multiple_of(case.output.monitor_every)
        };
        if datum_now && let Some(r) = reports.as_mut() {
            let sol = out_sol!();
            r.row(n, rp.time, &case, &sol, mesh, geo, spray.as_mut())?;
        }
        let o = &case.output;
        let fields_now = last
            && ((o.fields_every > 0 && n.is_multiple_of(o.fields_every))
                || o.fields_at.contains(&n)
                || sf.fired(t0, rp.time, tol)
                || req.result);
        let r = &case.run;
        let restart_numbered = r.restart_at.contains(&n)
            || Schedule {
                interval: None,
                times: &r.restart_times,
            }
            .fired(t0, rp.time, tol);
        let restart_now = last
            && ((r.restart_every > 0 && n.is_multiple_of(r.restart_every))
                || restart_numbered
                || sr.fired(t0, rp.time, tol)
                || asked);
        // the fields of a conjugate case are written on the whole mesh, the
        // probes are found there
        let full_fields = match (fields_now || at_monitor, &case.fluid_view) {
            (true, Some(_)) => {
                let st = fl.state()?;
                Some(crate::cht::full_state(&mut fl, &st)?)
            }
            _ => None,
        };
        if at_monitor || fields_now {
            let sol = out_sol!();
            if at_monitor {
                let boundaries = boundary_summaries(&sol, mesh, geo, &case);
                let mut row = format!("{n},{:e}", rp.time);
                for w in &walls {
                    if let Some(p) = boundaries.iter().find(|p| &p.name == w) {
                        let _ = write!(row, ",{:e},{:e},{:e}", p.force[0], p.force[1], p.force[2]);
                        if !rotor_names.is_empty() {
                            let _ = write!(row, ",{:e}", p.torque.unwrap_or(0.0));
                        }
                    }
                }
                if !rotor_names.is_empty() {
                    let rv = crate::rotor::evaluate(&sol, mesh, geo, &case);
                    for r in &rotor_names {
                        let v = rv.rotor(r);
                        let _ = write!(
                            row,
                            ",{:e},{:e}",
                            v.map_or(f64::NAN, |v| v.torque),
                            v.map_or(f64::NAN, |v| v.power)
                        );
                        if case
                            .sliding
                            .iter()
                            .any(|m| &crate::rotor::rotor_name(&m.zone) == r)
                        {
                            let _ = write!(
                                row,
                                ",{:e}",
                                v.map_or(f64::NAN, |v| (v.omega * rp.time).to_degrees())
                            );
                        }
                    }
                    for m in &case.output.machines {
                        match rv.machine(&m.name) {
                            Some(x) => {
                                let _ = write!(
                                    row,
                                    ",{:e},{:e},{:e},{:e},{:e}",
                                    x.volume_flow,
                                    x.total_pressure_rise,
                                    x.head,
                                    x.shaft_power,
                                    x.efficiency
                                );
                            }
                            None => row.push_str(",nan,nan,nan,nan,nan"),
                        }
                    }
                }
                let cc = monitor::coefficients(&case, &boundaries);
                if let Some([cd, cl]) = cc {
                    let _ = write!(row, ",{cd:e},{cl:e}");
                }
                if let Some(fh) = forces.as_mut() {
                    fh.row(&row);
                }
                let here = moving.then(|| probes.relocated(mesh, geo));
                let pr = here.as_ref().unwrap_or(&probes);
                let (vals, temps) = match (&full_probes, &full_fields) {
                    (Some((fp, v)), Some(st)) => (fp.values(st, &v.geo), fp.temperatures(st)),
                    _ => (pr.values(&sol, geo), pr.temperatures(&sol)),
                };
                if let Some(pc) = probe_csv.as_mut() {
                    let rel = frames.as_ref().map(|f| pr.relative(f, mesh, &vals));
                    let mut row = format!("{n},{:e}", rp.time);
                    for (k, (u, p)) in vals.iter().enumerate() {
                        let _ = write!(row, ",{:e},{:e},{:e},{p:e}", u[0], u[1], u[2]);
                        if let Some(r) = &rel {
                            let _ = write!(row, ",{:e},{:e},{:e}", r[k][0], r[k][1], r[k][2]);
                        }
                        if probe_temperature {
                            let _ = write!(row, ",{:e}", temps[k]);
                        }
                    }
                    pc.row(&row);
                }
                if let (Some(_), Some([cd, cl])) =
                    (bench, benchmark_coefficients(&case, &boundaries, &[], mesh))
                {
                    let (i, j) = (
                        probes.index("benchmark_front").unwrap_or(0),
                        probes.index("benchmark_back").unwrap_or(0),
                    );
                    bt.push(rp.time);
                    bcd.push(cd);
                    bcl.push(cl);
                    bdp.push(vals[i].1 - vals[j].1);
                }
            }
            if fields_now {
                let path = out.join(format!("fields_{n:06}.vtkhdf"));
                match &full_fields {
                    Some(st) => write_full_fields(&fl, st, &path, Some((n, Some(rp.time))))?,
                    None if !lead => {}
                    None => write_fields(
                        &case,
                        &path,
                        mesh,
                        geo,
                        &sol,
                        case.analytic_at(rp.time).as_ref(),
                        Some((n, Some(rp.time))),
                    )?,
                }
                lg.note(&format!(
                    "fields of step {n} written (fields_{n:06}.vtkhdf)"
                ));
                if let Some(se) = series.as_mut() {
                    se.add(n, rp.time, "fields", &path)?;
                }
            }
        }
        // next time step, from the step just taken (h, which is shorter
        // than dt when it was cut to end at end_time): the Courant number
        // scales with the step, and the growth limit holds for the ratio
        // of two successive steps
        let tnow = case.time.context("[time]")?;
        if tnow.adaptive {
            dt = h;
            let mut grow = if rp.cfl > 0.0 {
                (tnow.max_cfl / rp.cfl).min(tnow.max_growth)
            } else {
                tnow.max_growth
            };
            // D65: the step of the error estimate, Δt (0.9/E)^⅓ for an
            // error of third order in Δt, inside the Courant and growth
            // bounds
            if let Some(e) = rp.error.filter(|_| tnow.error_tolerance.is_some()) {
                let f = if e > 0.0 {
                    0.9 * e.powf(-1.0 / 3.0)
                } else {
                    tnow.max_growth
                };
                grow = grow.min(f.max(0.2));
            }
            dt = (dt * grow).max(tnow.min_time_step);
            if let Some(m) = tnow.max_time_step {
                dt = dt.min(m);
            }
        }
        if restart_now && restart::AVAILABLE && !moving {
            let mut rs = restart_of(
                &mut fl,
                &tr,
                smesh,
                name,
                dt,
                cycle.as_ref(),
                frozen.as_ref(),
                on.as_ref(),
            )?;
            rs.spray = crate::spray::state(spray.as_ref())?;
            if lead {
                restart::write(&out.join("restart.h5"), &rs)?;
            }
            if lead && restart_numbered {
                restart::write(&out.join(format!("restart_{n:06}.h5")), &rs)?;
            }
            if asked {
                lg.note(&format!("step {n}: restart.h5 written (SIGUSR1)"));
            }
        }
        if req.reread
            && let Some(old) = crate::case::reread_case(
                &mut case,
                &mut fl,
                lg,
                steer,
                warnings,
                n,
                out,
                sched.fired,
            )?
            && let (Some(a), Some(b)) = (old.time, case.time)
            && a.time_step != b.time_step
            && !b.adaptive
        {
            dt = b.time_step;
        }
        if let Some(old) =
            crate::case::apply_live_settings(&mut case, &mut fl, lg, steer, warnings, n)?
            && let (Some(a), Some(b)) = (old.time, case.time)
            && a.time_step != b.time_step
            && !b.adaptive
        {
            dt = b.time_step;
        }
        let stop_now = steer.poll(&case.run.limits, nsteps, true, req.stop);
        if let Some(s) = crate::ranks::agree_stop(on.as_ref(), stop_now) {
            lg.note(&format!("step {n}: {}", s.message));
            stopped = true;
            break;
        }
        if cycle_done {
            lg.note(&format!(
                "step {n}: the cycle-averaged loop ended ([energy.cycle]); the run stops"
            ));
            break;
        }
    }
    // the mesh of the last rotor angle (D182 variant B)
    let (mesh, geo) = current(&work, &deform, mesh0, geo0);
    history.flush()?;
    balance_csv.flush()?;
    if let Some(f) = spray_csv.as_mut() {
        f.flush()?;
    }
    if let Some(f) = vof_csv.as_mut() {
        f.flush()?;
    }
    if let Some(t) = fl.vof.as_deref().map(|v| &v.totals)
        && t.clipped_volume.abs() > 1e-6 * t.volume_start.abs()
    {
        lg.warn(&format!(
            "{}: {:.3e} m³ of phase 1 clipped over the run ({:.2e} of its volume); \
             `vof.csv` counts it per step",
            crate::diag::VOF_CLIPPED,
            t.clipped_volume,
            t.clipped_volume / t.volume_start
        ));
    }
    if let Some(f) = chem_csv.as_mut() {
        f.flush()?;
    }
    if let Some(sp) = spray.as_ref() {
        let f = crate::spray::write_parcels(sp, out, name, tr.clock.step)?;
        if let (Some(f), Some(se)) = (f, series.as_mut()) {
            se.add(tr.clock.step, tr.clock.time, "parcels", &f)?;
        }
        crate::spray::finish(sp)?;
    }
    if let Some(r) = reports.as_mut() {
        r.flush()?;
    }
    let datum_summary = reports
        .as_ref()
        .map(|r| r.summary(spray.as_ref(), &mesh.cell_zone));
    let spray_summary = spray.as_ref().map(|sp| crate::spray::SpraySummary {
        coupling: spray_coupling.clone(),
        ..crate::spray::summary(sp, &mesh.cell_zone)
    });
    if let Some(f) = heat_csv.as_mut() {
        f.flush()?;
    }
    if let Some(f) = forces.as_mut() {
        f.flush()?;
    }
    if let Some(f) = probe_csv.as_mut() {
        f.flush()?;
    }
    let secs = ts.elapsed().as_secs_f64();
    let steps_total = tr.clock.step;
    {
        let mut f = crate::console::Frame::new("Time stepping").summary().kv(
            "Status",
            format!(
                "{} {nsteps} steps to t = {} s",
                if stopped { "STOPPED after" } else { "DONE:" },
                crate::console::sci(tr.clock.time, 5)
            ),
        );
        f.push_kv(
            "Time step",
            format!(
                "{} to {} s",
                crate::console::sci(if nsteps > 0 { dt_min } else { 0.0 }, 3),
                crate::console::sci(dt_max, 3)
            ),
        );
        f.push_kv("Largest Courant number", format!("{cfl_max:.3}"));
        f.push_kv(
            "Solve time",
            format!(
                "{secs:.2} s  (setup {setup:.2} s, {:.0} ns per cell and step)",
                1e9 * secs / (nc as f64 * nsteps.max(1) as f64)
            ),
        );
        if rejected > 0 {
            f.push_kv(
                "Repeated steps",
                format!("{rejected}  (time-integration error above the tolerance)"),
            );
        }
        lg.frame(&f);
    }
    if last && (case.run.restart_final || stopped) && restart::AVAILABLE && !moving {
        let (smesh, _) = match &on {
            Some(o) => (&o.local.mesh, &o.geo),
            None => (mesh, geo),
        };
        let mut rs = restart_of(
            &mut fl,
            &tr,
            smesh,
            name,
            dt,
            cycle.as_ref(),
            frozen.as_ref(),
            on.as_ref(),
        )?;
        rs.spray = crate::spray::state(spray.as_ref())?;
        if lead {
            restart::write(&out.join("restart.h5"), &rs)?;
        }
    }
    let scalars = crate::scalar::summaries_on(
        &mut fl,
        &case,
        mesh,
        geo,
        crate::case::dimension(&case),
        tr.clock.time,
        on.as_ref(),
    )?;
    let periodic = crate::case::periodic_summaries(&mut fl, mesh, lg)?;
    let turbulence = crate::case::turbulence_summary(
        &mut fl,
        mesh,
        on.as_ref().map(|o| (o, geo, &case)),
        lg,
        warnings,
    )?;
    // the thermal acceleration (D213): not restored (ACC002), a skin
    // through the wall while accelerated (ACC003)
    let accel_summary = match accel.as_mut() {
        Some(acc) => {
            let thick = acc.skin_check(&mut fl)?;
            if !thick.is_empty() {
                let d = crate::diag::Diagnostic::new(
                    crate::diag::ACC_SKIN,
                    format!(
                        "[energy.acceleration]: while accelerated the skin depth of the period was \
                         larger than the thickness V/A of {}; the swing reached through the wall \
                         until the restore",
                        thick.join(", ")
                    ),
                );
                lg.line(&format!("{}: {}", d.code, d.message));
                steer.event(d.code, steps_total, &d.message);
                warnings.push(d);
            }
            if !acc.restored() {
                let d = crate::diag::Diagnostic::new(
                    crate::diag::ACC_UNRESTORED,
                    "[energy.acceleration]: the run ended before the true ρ c_p of the solids was \
                     restored; the final state is that of the accelerated solids and the history \
                     is not physical (run longer, or set `restore` to a time)"
                        .to_string(),
                );
                lg.warn(&format!("{}: {}", d.code, d.message));
                steer.event(d.code, steps_total, &d.message);
                warnings.push(d);
            }
            Some(acc.summary(tr.clock.time))
        }
        None => None,
    };
    // the frozen-flow march (D118): its phases, a flow that did not settle
    // (FRZ002), time-dependent flow values (FRZ004)
    let frozen_summary = match frozen.as_mut() {
        Some(fz) => {
            let sm = fz.summary();
            for d in fz.diagnostics.drain(..) {
                lg.warn(&format!("{}: {}", d.code, d.message));
                steer.event(d.code, steps_total, &d.message);
                warnings.push(d);
            }
            lg.line(&format!(
                "frozen flow: {} frozen steps over {:.6e} s, {} coupled steps over {:.6e} s, {} \
                 refreshes, {:.2} energy solves per frozen step",
                sm.frozen_steps,
                sm.frozen_time,
                sm.coupled_steps,
                sm.coupled_time,
                sm.refreshes,
                sm.solves_mean
            ));
            Some(sm)
        }
        None => None,
    };
    // the cycle-averaged loop (D204): its updates and the faces of the
    // last update, the warnings of the loop (CYC002, CYC004)
    let cycle_summary = match cycle.as_mut() {
        Some(cy) => {
            let sm = cy.summary();
            for d in cy.diagnostics.drain(..) {
                lg.warn(&format!("{}: {}", d.code, d.message));
                steer.event(d.code, steps_total, &d.message);
                warnings.push(d);
            }
            if last {
                let (head, rows) = cy.csv();
                let mut f = History::new(
                    out,
                    "cycle_coupling.csv",
                    &crate::case::csv_preamble(
                        "cycle_coupling",
                        name,
                        None,
                        "update -, time s (end of the fluid phase), cycles - (of the fluid \
                         phase), periodic - (1: within cycle_tolerance), cycle_change - \
                         (relative change of the mean interface heat flow, last two cycles), \
                         <interface>_heat W (cycle mean, into the fluid), heat_change - \
                         (against the fluid phase before), wall_change_max K, wall_change_mean K \
                         (solid side of the interface faces, solid phase against the update \
                         before), wall_mean K, solid_iterations -, relaxation - (of the \
                         update), solid_imbalance W, solid_relative - (balance of the solids with \
                         the Robin heat or the defect of the mean-flow model), swing_max K \
                         (skin: largest surface swing amplitude), swing_change K (skin: largest \
                         change of a harmonic of the swing)",
                    ),
                    &head,
                    None,
                )?;
                for r in rows {
                    f.row(&r);
                }
                f.flush()?;
                let centroids = fl
                    .cht
                    .as_deref()
                    .map(|c| c.view.geo.face_centroid.clone())
                    .unwrap_or_default();
                let (head, rows) = cy.faces_csv(&centroids);
                let mut f = History::new(
                    out,
                    "cycle_faces.csv",
                    &crate::case::csv_preamble(
                        "cycle_faces",
                        name,
                        None,
                        "face - (of the full mesh), interface -, x y z m (face centroid), area \
                         m², h_fluid W/(m² K) (cycle mean of the fluid side), conductance W/(m² \
                         K) (fluid side, contact and solid half cell in series), t_ref K, \
                         heat_flux W/m² (cycle mean, into the solid), t_wall K (solid side \
                         after the last solid phase); with skin = \"effusivity\": \
                         t_wall_amplitude K (half the range of the estimated surface swing), \
                         t_wall_min K, t_wall_max K; with resolved_cycles: \
                         resolved_t_wall_min K, resolved_t_wall_max K (over the last fully \
                         coupled cycle)",
                    ),
                    &head,
                    None,
                )?;
                for r in rows {
                    f.row(&r);
                }
                f.flush()?;
            }
            Some(sm)
        }
        None => None,
    };
    // conjugate heat transfer: the fields of the output on the full mesh
    let full_out = if case.fluid_view.is_some() {
        let st = fl.state()?;
        Some(crate::cht::full_state(&mut fl, &st)?)
    } else {
        None
    };
    let density = crate::density::summary(&mut fl, warnings)?;
    let species_sum = crate::species::summary(&mut fl, mesh, geo, true)?;
    let species_res = fl
        .scalars()
        .and_then(|s| s.species())
        .map(|l| l.res.clone());
    let sol = out_sol!();
    let analytic_end = case.analytic_at(tr.clock.time);
    let errors = analytic_end.map(|a| error_norms(&sol, mesh, geo, &a, dimension(&case)));
    let mut verification = crate::console::Frame::new("Verification").summary();
    if let Some(e) = &errors {
        verification.push_kv(
            &format!("Error at t = {} s", crate::console::sci(tr.clock.time, 4)),
            format!(
                "velocity L2 {}, Linf {} (relative), pressure Linf {}",
                crate::console::sci3(e.l2),
                crate::console::sci3(e.linf),
                crate::console::sci3(e.pressure_linf)
            ),
        );
    }
    // the imprinted sliding mesh (D182): the cells of the rotor have moved,
    // points are found on the mesh of the end of the run
    let loc_end = moving.then(|| Locator::new(mesh, geo));
    let loc_end = loc_end.as_ref().unwrap_or(&loc);
    let reference = case
        .verification
        .reference
        .map(|r| crate::case::compare_reference(&sol, mesh, geo, loc_end, r, &case));
    if let Some(r) = &reference {
        verification.push_kv(
            &format!("Reference {}", r.reference),
            format!(
                "largest deviation u {:.4}, v {:.4}  (lid speeds)",
                r.u_max_dev, r.v_max_dev
            ),
        );
        if last {
            crate::case::write_reference_svg(out, name, r)?;
        }
    }
    let boundaries = boundary_summaries(&sol, mesh, geo, &case);
    let rotor_values = crate::case::rotor_report(&sol, mesh, geo, &case, lg, warnings);
    let wall_probe_values = crate::case::wall_probes(&case, mesh, geo, &sol);
    let wall_crossing_values = crate::case::wall_crossings(&case, mesh, geo, &sol);
    let total: f64 = boundaries.iter().map(|p| p.mass_flow).exact_sum()
        // less the mass an evaporating spray gives the gas (SP4)
        - fl
            .spray_mass
            .as_ref()
            .map_or(0.0, |r| r.iter().map(|x| x.1).exact_sum());
    let umax = (0..nc).map(|c| norm(sol.u.get(c))).fold(0.0f64, f64::max);
    let big = boundaries
        .iter()
        .map(|p| p.mass_flow.abs().max(case.rho() * umax * p.area))
        .fold(0.0f64, f64::max);
    // in a reacting flow the mass in the domain keeps changing (a flame
    // turns dense fresh gas into light burnt gas): the check is then the
    // mass account of the last step, outflow plus storage (D119); other
    // variable-density cases keep the net outflow (a closed box has no
    // flow to scale the account by)
    let mass_imbalance = match (&last_balance, case.chemistry.is_some()) {
        (Some(b), true) => b.mass.relative,
        _ if big > 0.0 => total / big,
        _ => 0.0,
    };
    let coefficients = monitor::coefficients(&case, &boundaries);
    if let Some([cd, cl]) = coefficients {
        verification.push_kv("Force coefficients", format!("C_D {cd:.5}, C_L {cl:.5}"));
    }
    // the chart of the coefficients shows the evaluation window of a periodic
    // benchmark (the impulsive start is orders of magnitude larger)
    let mut chart_start = f64::NEG_INFINITY;
    let benchmark = match bench {
        Some(b @ Benchmark::SchaferTurek2d1) => {
            benchmark_coefficients(&case, &boundaries, &[], mesh).map(|[cd, cl]| {
                let vals = if moving {
                    probes.relocated(mesh, geo).values(&sol, geo)
                } else {
                    probes.values(&sol, geo)
                };
                let (i, j) = (
                    probes.index("benchmark_front").unwrap_or(0),
                    probes.index("benchmark_back").unwrap_or(0),
                );
                monitor::steady_benchmark(b, cd, cl, vals[i].1 - vals[j].1)
            })
        }
        Some(b @ Benchmark::SchaferTurek2d2) => {
            let start = case
                .verification
                .evaluation_start
                .unwrap_or(0.5 * (time0.start_time + tr.clock.time));
            chart_start = start;
            monitor::periodic_benchmark(b, &bt, &bcd, &bcl, &bdp, start)
        }
        None => None,
    };
    if let Some(b) = &benchmark {
        log_benchmark(b, |l| verification.push_text(l.trim_start().to_string()));
    }
    if !verification.is_empty() {
        lg.frame(&verification);
    }
    lg.frame(&crate::logblocks::results(&boundaries));
    let run_balance = ledger.run();
    lg.frame(&crate::logblocks::balance_block(
        &format!(
            "Balance of the run, t = {} to {} s ({} steps, {} with negative numerical dissipation)",
            crate::console::sci(
                run_balance.time.unwrap_or(0.0) - run_balance.duration.unwrap_or(0.0),
                4
            ),
            crate::console::sci(run_balance.time.unwrap_or(0.0), 4),
            ledger.steps,
            ledger.negative_steps
        ),
        &run_balance,
        true,
    ));
    lg.frame(&crate::logblocks::run_time(
        &fl.timings,
        nc,
        nsteps,
        "steps",
    ));
    if steer.record.on() {
        let mut row = balance::record_row(&run_balance);
        if let Some(m) = row.as_object_mut() {
            m.insert("type".into(), "balance_run".into());
        }
        steer.record.row(&row);
    }
    let mut last_balance = last_balance.unwrap_or_else(|| run_balance.clone());
    // the thermal account of the last step (D119, D161), from the heat
    // balance of the energy equation
    if let Some(sm) = scalars.iter().find(|s| s.energy.is_some())
        && let Some(e) = &sm.energy
    {
        last_balance.book(balance::thermal_account(e, &sm.boundaries));
    }
    // the species and element accounts of the last step (SF1, D119)
    if let Some(sp) = species_sum {
        for l in crate::species::summary_lines(&sp) {
            lg.line(&l);
        }
        if let Some(d) = crate::species::closure_finding(&sp) {
            lg.warn(&format!("{}: {}", d.code, d.message));
            warnings.push(d);
        }
        if let Some(r) = &species_res {
            for a in crate::species::accounts(&sp.flows, r) {
                last_balance.book(a);
            }
        }
    }
    let balance_summary = BalanceSummary {
        last: last_balance,
        largest_mass_relative: ledger.largest_mass_relative,
        steps: Some(ledger.steps),
        negative_steps: Some(ledger.negative_steps),
        run: Some(run_balance),
    };
    for d in balance::check(&balance_summary, &case.balance) {
        if d.code == crate::diag::BAL_MASS_STEP && balance_warned {
            warnings.push(d);
            continue;
        }
        lg.warn(&format!("{}: {}", d.code, d.message));
        steer.event(d.code, steps_total, &d.message);
        warnings.push(d);
    }
    if last || stopped {
        if let Some(rb) = &balance_summary.run {
            std::fs::write(out.join("balance.svg"), balance::chart(name, rb, true))?;
        }
        match (&full_out, fl.cht.as_deref()) {
            (Some(st), Some(c)) => {
                let v = &c.view;
                let full = Solution {
                    state: st.clone(),
                    disc: &c.fe.disc,
                    phys: c.fe.physics(),
                };
                write_outputs(
                    &v.full_case,
                    out,
                    &v.mesh,
                    &v.geo,
                    &full,
                    None,
                    &boundaries,
                    &Locator::new(&v.mesh, &v.geo),
                    Some((tr.clock.step, Some(tr.clock.time))),
                )?;
            }
            _ if !lead => {}
            _ => write_outputs(
                &case,
                out,
                mesh,
                geo,
                &sol,
                analytic_end.as_ref(),
                &boundaries,
                loc_end,
                Some((tr.clock.step, Some(tr.clock.time))),
            )?,
        }
        let ff = out.join("fields.vtkhdf");
        if let Some(se) = series.as_mut()
            && ff.exists()
        {
            se.add(tr.clock.step, tr.clock.time, "fields", &ff)?;
        }
        if !bt.is_empty() {
            let window = |c: &[f64]| -> Vec<(f64, f64)> {
                bt.iter()
                    .copied()
                    .zip(c.iter().copied())
                    .filter(|&(t, _)| t >= chart_start)
                    .collect()
            };
            let (cd, cl) = (window(&bcd), window(&bcl));
            let title = if chart_start.is_finite() {
                format!("{name}: force coefficients from t = {chart_start} s")
            } else {
                format!("{name}: force coefficients")
            };
            std::fs::write(
                out.join("coefficients.svg"),
                plot::chart(
                    &title,
                    "t (s)",
                    "C",
                    &[
                        Series {
                            label: "C_D",
                            points: &cd,
                            markers: false,
                        },
                        Series {
                            label: "C_L",
                            points: &cl,
                            markers: false,
                        },
                    ],
                    false,
                    false,
                ),
            )?;
        }
    }
    for sm in &scalars {
        lg.line(&crate::scalar::summary_line(sm));
    }
    for (w, c) in wall_crossing_values.iter().zip(&case.output.wall_crossings) {
        lg.line(&crate::case::wall_crossing_line(w, c));
    }
    let statistics = match (&stat_cfg, stat_start, &stat_series) {
        (Some(c), Some(s0), Some(ser)) => Some(crate::statistics::finish(
            &case,
            c,
            s0,
            ser,
            stat_fields.as_ref(),
            out,
            name,
            mesh,
            (tr.clock.step, tr.clock.time),
            lg,
            warnings,
        )?),
        _ => None,
    };
    // the scalars of the sweeps: those of the full mesh in a conjugate case
    let (scalar_fields, scalar_volume): (Vec<Vec<f64>>, Vec<f64>) =
        match (&full_out, fl.cht.as_deref()) {
            (Some(st), Some(c)) => (
                st.scalars.iter().map(|x| x.phi.clone()).collect(),
                c.view.geo.cell_volume[..c.view.mesh.n_cells].to_vec(),
            ),
            _ => (
                sol.state.scalars.iter().map(|x| x.phi.clone()).collect(),
                geo.cell_volume.clone(),
            ),
        };
    let immersed = crate::immersed::forces(&case, mesh, geo, &sol.u, &sol.phys)?;
    let chemistry = match chem.as_ref() {
        Some(ch) => {
            let s = ch.summary(&case, &boundaries)?;
            match s.flame_speed {
                Some(v) => {
                    let w = s.flame_speed_window.unwrap_or_default();
                    lg.line(&format!(
                        "chemistry {}: flame speed {v:.5} m/s (fuel consumption, mean over \
                         {:.4e} … {:.4e} s), heat release {:.4e} W",
                        s.mechanism, w[0], w[1], s.heat_release
                    ));
                    verification.push_kv("Flame speed", format!("{v:.5} m/s"));
                }
                None => lg.line(&format!(
                    "chemistry {}: heat release {:.4e} W",
                    s.mechanism, s.heat_release
                )),
            }
            Some(s)
        }
        None => None,
    };
    let solve = Solve {
        chemistry,
        immersed,
        topology: None,
        density,
        scalars,
        periodic,
        wall_probes: wall_probe_values,
        wall_crossings: wall_crossing_values,
        turbulence,
        statistics,
        remediation,
        operators,
        balance: Some(balance_summary),
        mesh: desc.into(),
        mesh_hash: crate::engine::mesh_hash(mesh),
        device: fl.device(),
        cells: nc,
        iterations: nsteps,
        converged: !stopped,
        residuals: last_rep.residuals,
        ns_per_cell_iteration: 1e9 * secs / (nc as f64 * nsteps.max(1) as f64),
        timings: fl.timings.clone(),
        amg_setups: fl.amg_setups + setups_past.0,
        sweep_setups: fl.sweep_setups + setups_past.1,
        turbulence_setups: fl.turbulence.as_ref().map_or(0, |t| t.reduction_setups()),
        errors,
        reference,
        boundaries,
        mass_imbalance,
        coefficients,
        rotors: rotor_values,
        transient: Some(TransientSummary {
            scheme: time0.scheme.label().into(),
            steps: nsteps,
            last_step: steps_total,
            time: tr.clock.time,
            dt_min: if nsteps > 0 { dt_min } else { 0.0 },
            dt_max,
            cfl_max,
            outer_mean: outer_sum as f64 / nsteps.max(1) as f64,
            cg_mean: cg_sum as f64 / nsteps.max(1) as f64,
            rejected,
            restart_from,
            stopped,
            ns_per_cell_step: 1e9 * secs / (nc as f64 * nsteps.max(1) as f64),
            heat_relative_max,
            acceleration: accel_summary,
            cycle_coupling: cycle_summary,
            frozen_flow: frozen_summary,
            mesh_motion: deform.as_ref().map(|d| crate::motion::MotionSummary {
                volume: d.volume(),
                volume_ratio_min: d.volume_ratio.0,
                volume_ratio_max: d.volume_ratio.1,
                closure_max: gcl_worst,
            }),
        }),
        datums: datum_summary,
        spray: spray_summary,
        spray_steady: None,
        benchmark,
    };
    let volume = geo.cell_volume.clone();
    Ok(RunOut {
        u: sol.state.u.clone(),
        scalars: scalar_fields,
        scalar_volume,
        solve,
        volume,
        stopped,
        start_values,
    })
}

/// C_D and C_L of the benchmark cylinder: 2 F / (ρ Ū² D H) with H the depth
/// of the slab.
/// The fields of a conjugate case on the whole mesh (stage 7d6): the state
/// `st` of [`crate::cht::full_state`] with the full case, mesh and the
/// energy view's discretisation.
fn write_full_fields(
    fl: &Flow,
    st: &crate::simplec::State,
    path: &Path,
    at: Option<(usize, Option<f64>)>,
) -> Result<()> {
    let c = fl
        .cht
        .as_deref()
        .context("the fields of a conjugate case without its energy view")?;
    let v = &c.view;
    let sol = Solution {
        state: st.clone(),
        disc: &c.fe.disc,
        phys: c.fe.physics(),
    };
    write_fields(&v.full_case, path, &v.mesh, &v.geo, &sol, None, at)
}

/// The front and back points of the benchmark (Δp), as extra probes.
pub(crate) fn benchmark_points(
    case: &FlowCase,
    mesh: &foehncraft_mesh::Mesh,
) -> Vec<(String, [f64; 3])> {
    let (lo, hi) = mesh.bounding_box();
    let zc = 0.5 * (lo[2] + hi[2]);
    match case.verification.benchmark {
        Some(_) => vec![
            (
                "benchmark_front".into(),
                [Benchmark::FRONT[0], Benchmark::FRONT[1], zc],
            ),
            (
                "benchmark_back".into(),
                [Benchmark::BACK[0], Benchmark::BACK[1], zc],
            ),
        ],
        None => Vec::new(),
    }
}

/// The steady benchmark (2D-1) of a converged solution.
pub(crate) fn steady_benchmark_of(
    case: &FlowCase,
    sol: &Solution,
    mesh: &foehncraft_mesh::Mesh,
    geo: &foehncraft_mesh::Geometry,
    loc: &Locator,
    boundaries: &[crate::case::BoundarySummary],
    immersed: &[crate::immersed::BodyForce],
) -> Option<monitor::BenchmarkResult> {
    let b = case.verification.benchmark?;
    let [cd, cl] = benchmark_coefficients(case, boundaries, immersed, mesh)?;
    let probes = Probes::new(case, &benchmark_points(case, mesh), loc, geo);
    let vals = probes.values(sol, geo);
    let (i, j) = (
        probes.index("benchmark_front")?,
        probes.index("benchmark_back")?,
    );
    Some(monitor::steady_benchmark(b, cd, cl, vals[i].1 - vals[j].1))
}

/// Log lines of a benchmark result.
pub(crate) fn log_benchmark(b: &monitor::BenchmarkResult, mut line: impl FnMut(&str)) {
    line(&format!("      {}:", b.benchmark));
    for v in &b.values {
        line(&format!(
            "        {:10} {:.5}  reference [{}, {}]  distance {:.4}",
            v.range.name, v.value, v.range.lower, v.range.upper, v.deviation
        ));
    }
}

/// The drag and lift coefficients of the benchmark body: the boundary
/// `benchmark_boundary` (default `cylinder`), or else the immersed body of
/// that name (L13 stage 1).
fn benchmark_coefficients(
    case: &FlowCase,
    boundaries: &[crate::case::BoundarySummary],
    immersed: &[crate::immersed::BodyForce],
    mesh: &foehncraft_mesh::Mesh,
) -> Option<[f64; 2]> {
    let b = case.verification.benchmark?;
    let name = case
        .verification
        .benchmark_boundary
        .clone()
        .unwrap_or_else(|| "cylinder".into());
    let force = match boundaries.iter().find(|p| p.name == name) {
        Some(p) => p.force,
        None => immersed.iter().find(|b| b.name == name)?.force,
    };
    let (lo, hi) = mesh.bounding_box();
    let depth = hi[2] - lo[2];
    let u = b.mean_velocity();
    let q = 0.5 * case.rho() * u * u * Benchmark::DIAMETER * depth;
    Some([force[0] / q, force[1] / q])
}

fn suffix_name(base: &str, _desc: &str, last: bool) -> String {
    if last {
        format!("{base}.csv")
    } else {
        format!("{base}_run.csv")
    }
}

/// The time scheme label (for reports).
pub fn scheme_label(s: TimeScheme) -> &'static str {
    s.label()
}
