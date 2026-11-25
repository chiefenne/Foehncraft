//! The steady spray mode of a flow case (D260, stage SP2s): parcel streams
//! through the steady gas field, the source-in-cell iteration of Crowe,
//! Sharma and Stock (1977) with the rules of `docs/requirements/
//! step28-sprays.md` section 5.9.
//!
//! With two-way coupling a stream pass runs every `[spray.steady] every`
//! outer iterations of the steady solver. The pass is the same pure
//! function as a transient pass (D247): streams from their injectors
//! through the frozen gas, records sorted and summed per cell (D249). The
//! gas holds its sources in value-plus-slope form (D250) as
//! S(u) = A − β u (momentum) and Q(T) = B − G T (heat), and relaxes them
//! between passes, A ← A + α (A_new − A) and the same for β, B and G. There
//! is no closing step: no parcel state is carried from pass to pass. The
//! *source residual* (the sources the gas holds against those of the last
//! pass, at the gas state of that pass) is part of the convergence test.

use crate::balance::{Account, Rates, Term, TermKind, names};
use crate::case::FlowCase;
use crate::spray::{FlowGas, GasHeat, MomentumSource, PAIR_HEAT, PAIR_MOMENTUM};
use anyhow::Result;
use foehncraft_mesh::Mesh;
use foehncraft_mesh::geometry::Geometry;
use foehncraft_particles::{Env, Spray, Steady, StreamBooks, StreamPass};

/// Units of `spray.csv` in a steady run.
pub const CSV_UNITS: &str = "iteration - (the outer iteration the pass ran before); pass -; \
     relaxation - (α of the pass; 1 for the first pass after a restart); source_residual - \
     (the larger of the momentum and heat parts, D260 (4)); residual_momentum -, residual_heat \
     - (Σ |S_held − S_pass| / Σ |S_pass| over the spray cells, at the gas state of the pass); \
     streams -; injected_kg_s, escaped_kg_s, deposited_kg_s, lost_kg_s, unfinished_kg_s kg/s \
     (liquid mass flows of the pass); residence_mean_s s (mass-flow mean time a stream is \
     followed), residence_max_s s; gas_force_x/y/z N (momentum the gas receives per second from \
     the pass), gas_heat_w W (heat the gas receives per second from the pass); coupling_max - \
     (largest β / ṁ of a spray cell: the slope of the momentum source over the mass flow \
     through the cell, Δt/τ with the residence time of the gas in the cell, D260 (5)); \
     coupling_heat_max - (largest G / (ṁ c_p)); liquid_fraction_max - (largest liquid volume \
     fraction of a cell); streams_per_cell_max -; substeps -";

/// The header of `spray.csv` in a steady run.
pub const CSV_HEAD: &str = "iteration,pass,relaxation,source_residual,residual_momentum,\
     residual_heat,streams,injected_kg_s,escaped_kg_s,deposited_kg_s,lost_kg_s,unfinished_kg_s,\
     residence_mean_s,residence_max_s,gas_force_x,gas_force_y,gas_force_z,gas_heat_w,\
     coupling_max,coupling_heat_max,liquid_fraction_max,streams_per_cell_max,substeps";

/// The sources the gas holds, per spray cell, ascending in the cell:
/// S(u) = A − β u (N), Q(T) = B − G T (W).
#[derive(Debug, Clone, Default, PartialEq)]
struct Held {
    cells: Vec<u32>,
    a: Vec<[f64; 3]>,
    beta: Vec<f64>,
    b: Vec<f64>,
    g: Vec<f64>,
}

impl Held {
    /// The momentum source of cell row `j` at the velocity `u`.
    fn momentum(&self, j: usize, u: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|k| self.a[j][k] - self.beta[j] * u[k])
    }

    /// The heat source of cell row `j` at the temperature `t`.
    fn heat(&self, j: usize, t: f64) -> f64 {
        self.b[j] - self.g[j] * t
    }
}

/// The source residual of a pass (D260 (4)).
#[derive(Debug, Clone, Copy, Default, PartialEq, serde::Serialize)]
pub struct SourceResidual {
    /// Σ |S_held − S_pass| / Σ |S_pass| of the momentum sources.
    pub momentum: f64,
    /// The same of the heat sources (0 without an energy equation).
    pub heat: f64,
}

impl SourceResidual {
    /// The larger part.
    pub fn max(&self) -> f64 {
        self.momentum.max(self.heat)
    }
}

/// The numbers of a pass for `spray.csv` and the summary.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct PassNumbers {
    /// Momentum the gas receives per second from the pass, Σ E (N).
    pub gas_force: [f64; 3],
    /// Heat the gas receives per second from the pass, Σ H (W).
    pub gas_heat: f64,
    /// Largest β / ṁ of a spray cell (D260 (5)).
    pub coupling_max: f64,
    /// Largest G / (ṁ c_p) of a spray cell.
    pub coupling_heat_max: f64,
    /// Largest liquid volume fraction of a cell.
    pub liquid_fraction_max: f64,
    /// Most streams through one cell.
    pub streams_per_cell_max: u32,
}

/// The steady spray in `result.json` (D260).
#[derive(Debug, Clone, serde::Serialize)]
pub struct SteadySummary {
    /// "steady".
    pub mode: &'static str,
    /// Two-way coupling.
    pub two_way: bool,
    /// Stream passes of the run (the final one included).
    pub passes: u32,
    /// The keys of `[spray.steady]`.
    pub keys: Steady,
    /// The source residual of the final pass (two-way).
    pub source_residual: SourceResidual,
    /// It is at or below `[spray.steady] tolerance` (always true one way).
    pub converged: bool,
    /// The books of the final pass (mass flows kg/s, momentum flows N,
    /// enthalpy flows W).
    pub books: StreamBooks,
    /// Per boundary: the liquid that leaves through it and sticks on it
    /// (kg/s).
    pub boundaries: Vec<BoundaryLiquid>,
    /// The numbers of the final pass.
    pub pass: PassNumbers,
    /// The momentum books of the streams closed: injected + weight +
    /// walls + removed + from the gas (N; zero to round-off).
    pub momentum_imbalance: [f64; 3],
    /// The same of the enthalpy (W).
    pub heat_imbalance: f64,
    /// How the numbers are defined.
    pub definitions: &'static str,
}

/// The liquid at one boundary (kg/s).
#[derive(Debug, Clone, serde::Serialize)]
pub struct BoundaryLiquid {
    /// The boundary.
    pub boundary: String,
    /// Leaves through it.
    pub escaped: f64,
    /// Sticks on it.
    pub deposited: f64,
}

/// The definitions written with the summary.
pub const DEFINITIONS: &str = "steady mode (D260): a stream stands for the droplets per second \
     of mass_flow / streams of its injector and is followed through the frozen steady gas until \
     it leaves, sticks or reaches max_time; mass flows in kg/s, momentum flows in N, enthalpy \
     flows (sensible, from 0 K) in W; the gas holds the sources relaxed between passes; the \
     source residual is Σ |S_held − S_pass| / Σ |S_pass| over the spray cells at the gas state \
     of the pass, the larger of momentum and heat";

/// The steady spray of a run.
pub struct SteadySpray {
    /// The spray (injectors, liquid, rules).
    pub sp: Spray,
    keys: Steady,
    two_way: bool,
    held: Held,
    /// Passes run.
    pub passes: u32,
    pass_at: Option<usize>,
    force: bool,
    full_first: bool,
    /// The source residual of the last pass.
    pub residual: Option<SourceResidual>,
    csv: String,
}

impl SteadySpray {
    /// The steady spray of a case without `[time]`, or `None` without
    /// `[spray]`; `restarted`: the first pass takes its sources in full
    /// (the relaxed sources are not in the restart file).
    pub fn create(
        case: &FlowCase,
        mesh: &Mesh,
        geo: &Geometry,
        loc: &foehncraft_mesh::locate::Locator,
        restarted: bool,
    ) -> Result<Option<Self>> {
        let Some(sp) = crate::spray::create(case, mesh, geo, loc)? else {
            return Ok(None);
        };
        let keys = case.spray.steady.unwrap_or_default();
        Ok(Some(Self {
            two_way: case.spray.coupling.two_way(),
            sp,
            keys,
            held: Held::default(),
            passes: 0,
            pass_at: None,
            force: false,
            full_first: restarted,
            residual: None,
            csv: String::new(),
        }))
    }

    /// The keys.
    pub fn keys(&self) -> &Steady {
        &self.keys
    }

    /// True with two-way coupling.
    pub fn two_way(&self) -> bool {
        self.two_way
    }

    /// Runs a pass before outer iteration `it` when one is due: the first
    /// iteration of the run, then every `every` iterations, and at once
    /// when the gas has converged on sources that are not (D260 (2), (4)).
    pub fn before_iteration(
        &mut self,
        fl: &mut crate::simplec::Flow,
        env: &Env,
        it: usize,
        first: usize,
    ) -> Result<()> {
        if !self.two_way {
            return Ok(());
        }
        let due =
            self.passes == 0 || self.force || (it - first - 1).is_multiple_of(self.keys.every);
        if !due {
            return Ok(());
        }
        self.force = false;
        let heat = GasHeat::of(fl)?;
        let pass = {
            let sol = crate::simplec::Solution {
                state: fl.state()?,
                disc: &fl.disc,
                phys: fl.physics(),
            };
            let gas = FlowGas {
                sol: &sol,
                geo: env.geo,
                heat: heat.as_ref(),
            };
            let pass = self.sp.streams(env, &gas, &self.keys, false)?;
            let alpha = if self.passes == 0 && self.full_first {
                1.0
            } else {
                self.keys.relaxation
            };
            let r = self.relax(&pass, &sol, heat.as_ref(), alpha);
            let n = numbers(&pass, &sol, env.geo, heat.as_ref());
            self.passes += 1;
            self.residual = Some(r);
            self.csv_row(it, alpha, &r, &pass.books, &n);
            pass
        };
        drop(pass);
        self.apply(fl)?;
        self.pass_at = Some(it);
        Ok(())
    }

    /// Whether the spray lets the run stop after outer iteration `it`
    /// whose gas residuals are below the tolerance (`gas_done`): a pass
    /// ran before this iteration and its source residual is at or below
    /// the tolerance. Otherwise, with the gas converged, the next
    /// iteration runs a pass at once.
    pub fn settled(&mut self, it: usize, gas_done: bool) -> bool {
        if !self.two_way {
            return true;
        }
        let ok = self.pass_at == Some(it)
            && self
                .residual
                .is_some_and(|r| r.max() <= self.keys.tolerance);
        if gas_done && !ok {
            self.force = true;
        }
        ok
    }

    /// The source residual against the held sources, then the relaxation
    /// of the held sources towards those of `pass` (α).
    fn relax(
        &mut self,
        pass: &StreamPass,
        sol: &crate::simplec::Solution,
        heat: Option<&GasHeat>,
        alpha: f64,
    ) -> SourceResidual {
        let s = &pass.sources;
        let old = std::mem::take(&mut self.held);
        let mut new = Held::default();
        let (mut dm, mut sm, mut dh, mut sh) = (0.0, 0.0, 0.0, 0.0);
        // the union of the cells, ascending
        let (mut i, mut j) = (0usize, 0usize);
        while i < old.cells.len() || j < s.cells.len() {
            let c = match (old.cells.get(i), s.cells.get(j)) {
                (Some(&a), Some(&b)) => a.min(b),
                (Some(&a), None) => a,
                (None, Some(&b)) => b,
                (None, None) => break,
            };
            let u = sol.state.u.get(c as usize);
            let t = heat.map_or(0.0, |h| h.t[c as usize]);
            let (a0, beta0, b0, g0) = if old.cells.get(i) == Some(&c) {
                i += 1;
                (old.a[i - 1], old.beta[i - 1], old.b[i - 1], old.g[i - 1])
            } else {
                ([0.0; 3], 0.0, 0.0, 0.0)
            };
            let (e, beta, h, g) = if s.cells.get(j) == Some(&c) {
                j += 1;
                let k = j - 1;
                (s.impulse[k], s.slope[k], s.heat[k], s.heat_slope[k])
            } else {
                ([0.0; 3], 0.0, 0.0, 0.0)
            };
            // the pass in the held form: A = E + β u*, B = H + G T*
            let a1: [f64; 3] = std::array::from_fn(|k| e[k] + beta * u[k]);
            let b1 = h + g * t;
            for k in 0..3 {
                dm += (a0[k] - beta0 * u[k] - e[k]).abs();
                sm += e[k].abs();
            }
            if heat.is_some() {
                dh += (b0 - g0 * t - h).abs();
                sh += h.abs();
            }
            new.cells.push(c);
            new.a
                .push(std::array::from_fn(|k| a0[k] + alpha * (a1[k] - a0[k])));
            new.beta.push(beta0 + alpha * (beta - beta0));
            new.b.push(b0 + alpha * (b1 - b0));
            new.g.push(g0 + alpha * (g - g0));
        }
        self.held = new;
        let ratio = |d: f64, s: f64| {
            if s > 0.0 {
                d / s
            } else if d > 0.0 {
                1.0
            } else {
                0.0
            }
        };
        SourceResidual {
            momentum: ratio(dm, sm),
            heat: ratio(dh, sh),
        }
    }

    /// Gives the gas the held sources.
    fn apply(&self, fl: &mut crate::simplec::Flow) -> Result<()> {
        let h = &self.held;
        fl.spray_momentum = Some(MomentumSource {
            cells: h.cells.clone(),
            e: h.a.clone(),
            beta: h.beta.clone(),
            ustar: vec![[0.0; 3]; h.cells.len()],
        });
        if let Some(sc) = fl.scalars.as_mut() {
            let rows: Vec<(u32, f64, f64)> = h
                .cells
                .iter()
                .enumerate()
                .map(|(j, &c)| (c, h.b[j], h.g[j]))
                .collect();
            sc.set_energy_spray(fl.dev.as_mut(), Some(rows))?;
        }
        Ok(())
    }

    fn csv_row(
        &mut self,
        it: usize,
        alpha: f64,
        r: &SourceResidual,
        b: &StreamBooks,
        n: &PassNumbers,
    ) {
        self.csv.push_str(&format!(
            "{it},{},{alpha:e},{:e},{:e},{:e},{},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},\
             {:e},{:e},{:e},{:e},{},{}\n",
            self.passes,
            r.max(),
            r.momentum,
            r.heat,
            b.streams,
            b.injected_mass,
            b.escaped_mass.iter().sum::<f64>(),
            b.deposited_mass.iter().sum::<f64>(),
            b.lost_mass,
            b.unfinished_mass,
            b.residence_mean,
            b.residence_max,
            n.gas_force[0],
            n.gas_force[1],
            n.gas_force[2],
            n.gas_heat,
            n.coupling_max,
            n.coupling_heat_max,
            n.liquid_fraction_max,
            n.streams_per_cell_max,
            b.substeps
        ));
    }

    /// The final pass on the converged gas (with the paths): its books,
    /// its source residual against the held sources (no relaxation: the
    /// gas is final), its `spray.csv` row and the summary.
    pub fn finish(
        &mut self,
        sol: &crate::simplec::Solution,
        env: &Env,
        heat: Option<&GasHeat>,
        it: usize,
    ) -> Result<(SteadySummary, StreamPass)> {
        let gas = FlowGas {
            sol,
            geo: env.geo,
            heat,
        };
        let pass = self.sp.streams(env, &gas, &self.keys, true)?;
        let n = numbers(&pass, sol, env.geo, heat);
        let r = if self.two_way {
            // the residual without changing what the gas holds
            let held = self.held.clone();
            let r = self.relax(&pass, sol, heat, 0.0);
            self.held = held;
            r
        } else {
            SourceResidual::default()
        };
        self.passes += 1;
        self.residual = Some(r);
        self.csv_row(it, 0.0, &r, &pass.books, &n);
        let b = &pass.books;
        let summary = SteadySummary {
            mode: "steady",
            two_way: self.two_way,
            passes: self.passes,
            keys: self.keys,
            source_residual: r,
            converged: !self.two_way || r.max() <= self.keys.tolerance,
            books: b.clone(),
            boundaries: self
                .sp
                .boundary_names()
                .iter()
                .enumerate()
                .map(|(k, name)| BoundaryLiquid {
                    boundary: name.clone(),
                    escaped: b.escaped_mass[k],
                    deposited: b.deposited_mass[k],
                })
                .collect(),
            pass: n,
            momentum_imbalance: b.momentum.imbalance(),
            heat_imbalance: b.heat_imbalance(),
            definitions: DEFINITIONS,
        };
        Ok((summary, pass))
    }

    /// The rows of `spray.csv`.
    pub fn csv(&self, name: &str) -> String {
        let mut s = crate::case::csv_preamble("spray", name, None, CSV_UNITS);
        s.push_str(CSV_HEAD);
        s.push('\n');
        s.push_str(&self.csv);
        s
    }

    /// The ledger rows of the steady spray on the balance of the final
    /// state (D119, D244, D260 (4)): the pairs hold the sources of the
    /// final pass on both sides (zero to the bit); the gas accounts add
    /// the difference between what the gas holds and that pass as the
    /// `spray source residual`, so they close with what the gas received.
    pub fn book(
        &self,
        r: &mut Rates,
        pass: &StreamPass,
        sol: &crate::simplec::Solution,
        heat: Option<&GasHeat>,
    ) {
        if !self.two_way {
            return;
        }
        let s = &pass.sources;
        let b = &pass.books;
        let mut gas = [0.0; 3];
        for e in &s.impulse {
            for k in 0..3 {
                gas[k] += e[k];
            }
        }
        // what the gas holds at the final state, less the final pass
        let (mut held, mut power) = ([0.0; 3], 0.0);
        let (mut held_q, mut gas_q) = (0.0, 0.0);
        for (j, &c) in self.held.cells.iter().enumerate() {
            let u = sol.state.u.get(c as usize);
            let f = self.held.momentum(j, u);
            for k in 0..3 {
                held[k] += f[k];
                power += u[k] * f[k];
            }
            if let Some(h) = heat {
                held_q += self.held.heat(j, h.t[c as usize]);
            }
        }
        for q in &s.heat {
            gas_q += q;
        }
        let account = |quantity: &str, region: &str, terms: Vec<Term>| {
            let mut a = Account::new(
                quantity,
                if quantity.contains("momentum") {
                    "kg m/s"
                } else {
                    "J"
                },
            );
            a.region = region.into();
            a.terms = terms;
            a
        };
        for k in 0..3 {
            r.book(account(
                names::MOMENTUM[k],
                "domain",
                vec![
                    Term::volume("spray drag", TermKind::Exchange, gas[k]).paired(PAIR_MOMENTUM[k]),
                    Term::volume("spray source residual", TermKind::Source, held[k] - gas[k]),
                ],
            ));
            let m = &b.momentum;
            r.book(account(
                &format!("parcel {}", names::MOMENTUM[k]),
                "parcels",
                vec![
                    Term::volume("injected", TermKind::Source, m.injected[k]),
                    Term::volume("weight", TermKind::Source, m.weight[k]),
                    Term::volume("walls", TermKind::Source, m.walls[k]),
                    Term::volume("left the flow", TermKind::Source, m.removed[k]),
                    Term::volume("from the gas", TermKind::Exchange, m.from_gas[k])
                        .paired(PAIR_MOMENTUM[k]),
                ],
            ));
        }
        r.book(account(
            names::KINETIC,
            "domain",
            vec![Term::volume("spray drag power", TermKind::Source, power)],
        ));
        if heat.is_none() {
            return;
        }
        r.book(account(
            names::THERMAL,
            "domain",
            vec![
                Term::volume("spray heat", TermKind::Exchange, gas_q).paired(PAIR_HEAT),
                Term::volume("spray source residual", TermKind::Source, held_q - gas_q),
            ],
        ));
        r.book(account(
            "parcel enthalpy",
            "parcels",
            vec![
                Term::volume("injected", TermKind::Source, b.enthalpy[0]),
                Term::volume("left the flow", TermKind::Source, b.enthalpy[1]),
                Term::volume("from the gas", TermKind::Exchange, b.enthalpy[2]).paired(PAIR_HEAT),
            ],
        ));
    }

    /// Writes `spray/streams.vtkhdf`; returns its path (`None` in a build
    /// without HDF5 or without paths).
    pub fn write_streams(
        &self,
        pass: &StreamPass,
        out: &std::path::Path,
        name: &str,
        it: usize,
    ) -> Result<Option<std::path::PathBuf>> {
        #[cfg(feature = "hdf5")]
        {
            let Some(p) = &pass.paths else {
                return Ok(None);
            };
            let dir = out.join("spray");
            std::fs::create_dir_all(&dir)?;
            let path = dir.join("streams.vtkhdf");
            self.sp.write_streams(&path, name, it as u64, p)?;
            Ok(Some(path))
        }
        #[cfg(not(feature = "hdf5"))]
        {
            let _ = (pass, out, name, it);
            Ok(None)
        }
    }
}

/// The health numbers of a pass (D254 (3), (5) with the residence time of
/// the gas in the cell in place of the time step, D260 (5)).
fn numbers(
    pass: &StreamPass,
    sol: &crate::simplec::Solution,
    geo: &Geometry,
    heat: Option<&GasHeat>,
) -> PassNumbers {
    let s = &pass.sources;
    let mut n = PassNumbers::default();
    for (j, e) in s.impulse.iter().enumerate() {
        for k in 0..3 {
            n.gas_force[k] += e[k];
        }
        n.gas_heat += s.heat[j];
    }
    if !s.cells.is_empty() {
        // the mass flow through each spray cell, ½ Σ |ṁ_f| over its faces
        let d = sol.disc;
        let mut through = vec![0.0; s.cells.len()];
        for (f, &m) in sol.state.mass_flux.iter().enumerate().take(d.n_faces) {
            let m = 0.5 * m.abs();
            if let Some(j) = s.find(d.owner[f]) {
                through[j] += m;
            }
            if f < d.n_interior
                && let Some(j) = s.find(d.neighbour[f])
            {
                through[j] += m;
            }
        }
        for (j, &c) in s.cells.iter().enumerate() {
            let q = through[j].max(1e-300);
            n.coupling_max = n.coupling_max.max(s.slope[j] / q);
            if let Some(h) = heat {
                n.coupling_heat_max = n
                    .coupling_heat_max
                    .max(s.heat_slope[j] / (q * h.cp[c as usize]));
            }
        }
    }
    for &(c, k, v) in &pass.loading {
        n.streams_per_cell_max = n.streams_per_cell_max.max(k);
        n.liquid_fraction_max = n.liquid_fraction_max.max(v / geo.cell_volume[c as usize]);
    }
    n
}

/// The log frame of the steady spray.
pub fn frame(s: &SteadySummary) -> crate::console::Frame {
    use crate::console::sci3;
    let b = &s.books;
    let mut f = crate::console::Frame::new("Spray (steady mode, D260)").summary();
    f.push_kv(
        "Streams",
        format!(
            "{} ({} per injector), {} passes{}",
            b.streams,
            s.keys.streams,
            s.passes,
            if s.two_way {
                format!(
                    ", every {} iterations, relaxation {}",
                    s.keys.every, s.keys.relaxation
                )
            } else {
                " (one way: the final pass)".into()
            }
        ),
    );
    f.push_kv(
        "Liquid",
        format!(
            "injected {} kg/s, left {} kg/s, stuck {} kg/s, lost {} kg/s, unfinished {} kg/s",
            sci3(b.injected_mass),
            sci3(b.escaped_mass.iter().sum::<f64>()),
            sci3(b.deposited_mass.iter().sum::<f64>()),
            sci3(b.lost_mass),
            sci3(b.unfinished_mass)
        ),
    );
    f.push_kv(
        "Residence time",
        format!(
            "mean {} s, longest {} s",
            sci3(b.residence_mean),
            sci3(b.residence_max)
        ),
    );
    if s.two_way {
        let p = &s.pass;
        f.push_kv(
            "To the gas",
            format!(
                "force [{}, {}, {}] N, heat {} W",
                sci3(p.gas_force[0]),
                sci3(p.gas_force[1]),
                sci3(p.gas_force[2]),
                sci3(p.gas_heat)
            ),
        );
        f.push_kv(
            "Source residual",
            format!(
                "{} (momentum {}, heat {}), tolerance {}",
                sci3(s.source_residual.max()),
                sci3(s.source_residual.momentum),
                sci3(s.source_residual.heat),
                sci3(s.keys.tolerance)
            ),
        );
        f.push_kv(
            "Coupling",
            format!(
                "β/ṁ up to {}, G/(ṁ c_p) up to {}, liquid fraction up to {}",
                sci3(p.coupling_max),
                sci3(p.coupling_heat_max),
                sci3(p.liquid_fraction_max)
            ),
        );
    }
    f
}
