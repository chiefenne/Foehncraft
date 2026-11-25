//! The steady spray mode (D260, stage SP2s): parcel *streams* through a
//! frozen steady gas field.
//!
//! A stream stands for a number of droplets per second. It starts at its
//! injector with the direction and size a parcel of a transient run with
//! the same key would get (D121) and is followed until it leaves the flow,
//! sticks, is lost or reaches `[spray.steady] max_time`. Its mass is a mass
//! flow (kg/s), so the exchange records of its sub-steps hold momentum
//! flows (N) and heat flows (W), and their sorted per-cell sums (D249) are
//! the sources of the gas per second. A stream pass is a pure function of
//! the gas and the keys: it carries nothing from pass to pass, so every
//! pass on the same gas gives the same bits, whatever the thread count
//! (D247, D260 (3) and (4)).

use super::{Env, Fate, Gas, One, Spray, TrackOpts, direction, frame, mass_of, unit};
use crate::config::Steady;
use crate::coupling::{CellSources, Record};
use anyhow::{Result, bail};
use foehncraft_core::rng::{self, stream};
use rayon::prelude::*;
use serde::Serialize;

/// The momentum flows of the streams (N), per component: the parcel side
/// of the ledger of a steady spray. injected + weight + walls + removed +
/// from the gas = 0 to round-off.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize)]
pub struct MomentumFlows {
    /// Momentum flow of the streams at their injectors.
    pub injected: [f64; 3],
    /// Weight (gravity less buoyancy) along the paths.
    pub weight: [f64; 3],
    /// Walls (rebounds, mirrors).
    pub walls: [f64; 3],
    /// What the streams carry at the end of their paths, with the sign of
    /// a loss (left, stuck, lost or unfinished).
    pub removed: [f64; 3],
    /// Received from the gas (the opposite of what the gas receives,
    /// summed per cell in cell order; one way, where the gas books nothing,
    /// what the other terms leave).
    pub from_gas: [f64; 3],
}

impl MomentumFlows {
    /// The sum of the terms per component: zero to round-off.
    pub fn imbalance(&self) -> [f64; 3] {
        std::array::from_fn(|k| {
            self.injected[k] + self.weight[k] + self.walls[k] + self.removed[k] + self.from_gas[k]
        })
    }
}

/// The books of one stream pass: mass flows (kg/s) per fate, the momentum
/// flows (N) and the sensible enthalpy flows (W).
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct StreamBooks {
    /// Streams followed.
    pub streams: u64,
    /// Liquid mass flow of the injectors (kg/s).
    pub injected_mass: f64,
    /// Mass flow that left through each boundary of the mesh (kg/s).
    pub escaped_mass: Vec<f64>,
    /// Mass flow that sticks on each boundary of the mesh (kg/s).
    pub deposited_mass: Vec<f64>,
    /// Mass flow of streams that were lost (failed walk or more than
    /// `max_substeps` sub-steps).
    pub lost_mass: f64,
    /// Streams lost.
    pub lost_streams: u64,
    /// Mass flow of streams still in the flow at `max_time`.
    pub unfinished_mass: f64,
    /// Streams still in the flow at `max_time`.
    pub unfinished_streams: u64,
    /// Bounces off boundaries.
    pub rebounds: u64,
    /// Sub-steps of all streams.
    pub substeps: u64,
    /// Mass-flow mean of the time a stream is followed (s).
    pub residence_mean: f64,
    /// Longest time a stream was followed (s).
    pub residence_max: f64,
    /// The momentum flows.
    pub momentum: MomentumFlows,
    /// Sensible enthalpy flows Σ ṁ c_p T (W): [injected, removed (with the
    /// sign of a loss), from the gas (one way: what the other two leave)].
    pub enthalpy: [f64; 3],
    /// Liquid mass flow that evaporated (kg/s, the streams' own loss, net
    /// of condensation; SP5).
    #[serde(default)]
    pub evaporated_mass: f64,
    /// Streams that evaporated completely.
    #[serde(default)]
    pub evaporated_streams: u64,
    /// Liquid mass flow the saturation limiter of D253 held back in streams
    /// that the frozen gas evaporated completely (kg/s): it ends the pass
    /// as liquid where the stream ended (SP5).
    #[serde(default)]
    pub limited_mass: f64,
}

impl StreamBooks {
    /// injected + removed + from the gas of the enthalpy: zero to
    /// round-off.
    pub fn heat_imbalance(&self) -> f64 {
        self.enthalpy.iter().sum()
    }
}

/// A point of the stream paths of a pass, for the paths file.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Streams {
    /// Stream identity (injector in the top 16 bits, key in the low 48).
    pub id: Vec<u64>,
    /// Time since the injector (s).
    pub t: Vec<f64>,
    /// Position (m).
    pub x: Vec<[f64; 3]>,
    /// Velocity (m/s).
    pub v: Vec<[f64; 3]>,
    /// Droplet diameter (m).
    pub d: Vec<f64>,
    /// Droplet temperature (K).
    pub temp: Vec<f64>,
    /// Droplets per second of the stream.
    pub rate: Vec<f64>,
}

/// One stream pass (D260).
pub struct StreamPass {
    /// The per-cell sums of the records (two-way coupling; empty one way):
    /// momentum flows (N), slopes (kg/s), heat flows (W), heat slopes (W/K).
    pub sources: CellSources,
    /// The books.
    pub books: StreamBooks,
    /// Per cell with records: the cell, the streams that pass it and the
    /// liquid volume it holds (m³, Σ of the liquid volume flow of a stream
    /// times the time it spends in the cell), ascending in the cell.
    pub loading: Vec<(u32, u32, f64)>,
    /// The paths, when asked for.
    pub paths: Option<Streams>,
    /// The streams as tracked, in the order of their identities (for the
    /// limiter, [`Spray::limit_streams`]).
    outcomes: Vec<(One, super::Outcome)>,
}

impl Spray {
    /// The streams of the injectors (D260): `streams` per injector, stream
    /// k with the key `seed · streams + k`, its direction and size drawn as
    /// those of a transient parcel with that key, its droplets per second
    /// from the injector's mass flow shared equally.
    fn stream_release(&self, env: &Env, st: &Steady) -> Result<Vec<One>> {
        let rho_p = self.cfg.liquid.density;
        let n = u64::from(st.streams);
        let mut out = Vec::new();
        for (ii, inj) in self.cfg.injectors.iter().enumerate() {
            let axis = unit(inj.direction);
            let (e1, e2) = frame(axis);
            let cos_a = inj.half_angle.to_radians().cos();
            let m_stream = inj.mass_flow / n as f64;
            let Some(cell) = env.loc.locate_search(inj.position) else {
                bail!(
                    "[[spray.injector]] {}: the position lies outside the mesh",
                    inj.name
                );
            };
            for k in 0..n {
                let key = st.seed.wrapping_mul(n).wrapping_add(k) & ((1u64 << 48) - 1);
                let id = ((ii as u64) << 48) | key;
                let dir = direction(id, axis, e1, e2, cos_a);
                let d = match &inj.size {
                    Some(sz) => {
                        sz.quantile(rng::uniform_keyed([stream::SPRAY_INJECTION, id, 1, 0]))
                    }
                    None => inj.diameter,
                };
                out.push(One {
                    id,
                    x: inj.position,
                    v: [dir[0] * inj.speed, dir[1] * inj.speed, dir[2] * inj.speed],
                    d,
                    count: m_stream / (rho_p * std::f64::consts::PI / 6.0 * d * d * d),
                    born: 0.0,
                    temp: inj.temperature.unwrap_or(self.cfg.liquid.temperature),
                    eddy: super::Eddy::default(),
                    cell,
                    face: None,
                });
            }
        }
        Ok(out)
    }

    /// One stream pass through the frozen gas `gas` (D260) with the keys
    /// `st`; `paths` keeps the path of every stream. The spray is not
    /// changed.
    pub fn streams(
        &self,
        env: &Env,
        gas: &dyn Gas,
        st: &Steady,
        paths: bool,
    ) -> Result<StreamPass> {
        let fresh = self.stream_release(env, st)?;
        let opts = TrackOpts {
            paths,
            steady: true,
        };
        let mut outcomes: Vec<(One, super::Outcome)> = fresh
            .par_iter()
            .map(|&p| (p, self.track(env, gas, p, st.max_time, st.max_time, opts)))
            .collect();
        // the books in the order of the stream identities (D93)
        outcomes.sort_by_key(|(p, _)| p.id);
        Ok(self.stream_pass(env, outcomes, paths))
    }

    /// The books, the sources and the paths of the tracked streams
    /// `outcomes` (sorted by identity).
    fn stream_pass(
        &self,
        env: &Env,
        outcomes: Vec<(One, super::Outcome)>,
        paths: bool,
    ) -> StreamPass {
        let rho_p = self.cfg.liquid.density;
        let cp_p = self.cfg.liquid.specific_heat;
        let nb = env.mesh.boundaries.len();
        let mut b = StreamBooks {
            streams: outcomes.len() as u64,
            escaped_mass: vec![0.0; nb],
            deposited_mass: vec![0.0; nb],
            ..Default::default()
        };
        let mut residence = 0.0;
        for (p0, o) in &outcomes {
            let m = mass_of(p0.d, p0.count, rho_p);
            b.injected_mass += m;
            b.rebounds += o.rebounds;
            b.substeps += o.substeps;
            residence += m * o.used;
            b.residence_max = b.residence_max.max(o.used);
            b.evaporated_mass += o.vapour;
            let me = mass_of(o.p.d, o.p.count, rho_p);
            let mf = &mut b.momentum;
            for k in 0..3 {
                mf.injected[k] += m * p0.v[k];
                mf.weight[k] += o.weight[k];
                mf.walls[k] += o.walls[k];
                mf.removed[k] -= me * o.p.v[k];
            }
            b.enthalpy[0] += m * cp_p * p0.temp;
            b.enthalpy[1] -= me * cp_p * o.p.temp;
            match o.fate {
                Fate::Active => {
                    b.unfinished_mass += me;
                    b.unfinished_streams += 1;
                }
                Fate::Stuck(f) => {
                    if let Some(k) = env.mesh.boundary_of(f) {
                        b.deposited_mass[k] += me;
                    }
                }
                Fate::Escaped(f) => {
                    if let Some(k) = env.mesh.boundary_of(f) {
                        b.escaped_mass[k] += me;
                    }
                }
                Fate::Lost => {
                    b.lost_mass += me;
                    b.lost_streams += 1;
                }
                Fate::Evaporated if me > 0.0 => b.limited_mass += me,
                Fate::Evaporated => b.evaporated_streams += 1,
            }
        }
        if b.injected_mass > 0.0 {
            b.residence_mean = residence / b.injected_mass;
        }
        let two_way = self.cfg.coupling.two_way();
        let (sources, loading) = if two_way {
            let all: Vec<Record> = outcomes
                .iter()
                .flat_map(|(_, o)| o.records.iter().copied())
                .collect();
            let s = CellSources::from_records(&all);
            // the liquid in each cell: volume flow times the time spent
            let mut dwell: Vec<(u32, u64, f64)> = Vec::with_capacity(all.len());
            for (p0, o) in &outcomes {
                let q = p0.count * std::f64::consts::PI / 6.0 * p0.d.powi(3);
                dwell.extend(o.records.iter().map(|r| (r.cell, p0.id, q * r.dt)));
            }
            dwell.sort_by_key(|x| (x.0, x.1));
            let mut loading: Vec<(u32, u32, f64)> = Vec::new();
            let mut last: Option<(u32, u64)> = None;
            for (c, id, v) in dwell {
                match loading.last_mut() {
                    Some(l) if l.0 == c => {
                        if last != Some((c, id)) {
                            l.1 += 1;
                        }
                        l.2 += v;
                    }
                    _ => loading.push((c, 1, v)),
                }
                last = Some((c, id));
            }
            (s, loading)
        } else {
            (CellSources::default(), Vec::new())
        };
        if two_way {
            for (j, imp) in sources.impulse.iter().enumerate() {
                for (f, i) in b.momentum.from_gas.iter_mut().zip(imp) {
                    *f -= i;
                }
                b.enthalpy[2] -= sources.heat[j];
            }
        } else {
            // one way the gas books nothing: what the streams received from
            // it is what the other terms leave
            let m = &mut b.momentum;
            for k in 0..3 {
                m.from_gas[k] = -(m.injected[k] + m.weight[k] + m.walls[k] + m.removed[k]);
            }
            b.enthalpy[2] = -(b.enthalpy[0] + b.enthalpy[1]);
        }
        let paths = paths.then(|| {
            let mut s = Streams::default();
            for (p0, o) in &outcomes {
                for q in &o.path {
                    s.id.push(p0.id);
                    s.t.push(q.t);
                    s.x.push(q.x);
                    s.v.push(q.v);
                    s.d.push(o.p.d);
                    s.temp.push(q.temp);
                    s.rate.push(p0.count);
                }
            }
            s
        });
        StreamPass {
            sources,
            books: b,
            loading,
            paths,
            outcomes,
        }
    }

    /// The saturation limiter of D253 on a stream pass (D260 (5)): the
    /// evaporating records of each cell with a factor `factor(cell)` below
    /// 1 keep that share of their vapour; the rest stays liquid in its
    /// stream, at the temperature of the record, with the momentum it
    /// carried, and the books and the sources are made again. Returns the
    /// cells limited and the mass flow held back (kg/s).
    pub fn limit_streams(
        &self,
        env: &Env,
        pass: &mut StreamPass,
        factor: &dyn Fn(u32) -> f64,
    ) -> (u32, f64) {
        let cells = pass
            .sources
            .cells
            .iter()
            .filter(|&&c| factor(c) < 1.0)
            .count() as u32;
        if cells == 0 {
            return (0, 0.0);
        }
        let mut outcomes = std::mem::take(&mut pass.outcomes);
        let mut back = 0.0;
        for (_, o) in &mut outcomes {
            back += self.hold_back(o, factor, false);
        }
        let paths = pass.paths.is_some();
        *pass = self.stream_pass(env, outcomes, paths);
        (cells, back)
    }

    /// Writes the stream paths to a VTKHDF `PolyData` file of points
    /// (format `foehncraft.streams`): one point at the start of every
    /// sub-step and at the end of each path.
    #[cfg(feature = "hdf5")]
    pub fn write_streams(
        &self,
        path: &std::path::Path,
        case: &str,
        iteration: u64,
        s: &Streams,
    ) -> Result<()> {
        use foehncraft_mesh::vtkhdf::{Attr, PointField, write_points};
        let vel: Vec<f64> = s.v.iter().flatten().copied().collect();
        let id: Vec<f64> = s.id.iter().map(|&i| i as f64).collect();
        let injector: Vec<f64> = s.id.iter().map(|&i| (i >> 48) as f64).collect();
        let fields = [
            PointField {
                name: "velocity",
                values: &vel,
                components: 3,
            },
            PointField {
                name: "diameter",
                values: &s.d,
                components: 1,
            },
            PointField {
                name: "temperature",
                values: &s.temp,
                components: 1,
            },
            PointField {
                name: "droplet_rate",
                values: &s.rate,
                components: 1,
            },
            PointField {
                name: "residence_time",
                values: &s.t,
                components: 1,
            },
            PointField {
                name: "stream_id",
                values: &id,
                components: 1,
            },
            PointField {
                name: "injector",
                values: &injector,
                components: 1,
            },
        ];
        let attrs = [
            ("format", Attr::Str("foehncraft.streams")),
            ("format_version", Attr::U64(1)),
            ("case", Attr::Str(case)),
            ("iteration", Attr::U64(iteration)),
            ("liquid", Attr::Str(&self.cfg.liquid.name)),
            ("liquid_density", Attr::F64(self.cfg.liquid.density)),
            ("units", Attr::Str("m, m/s, K, s, droplets/s")),
            ("points", Attr::U64(s.id.len() as u64)),
        ];
        write_points(path, &s.x, &fields, &attrs)
    }
}
