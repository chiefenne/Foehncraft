//! The parcels and their tracking through the mesh.

use crate::config::{BoundaryMode, Injector, SpraySection};
use anyhow::{Result, bail};
use foehncraft_core::rng::{self, stream};
use foehncraft_kernels::particle::{ParcelInput, parcel_step};
use foehncraft_mesh::geometry::Geometry;
use foehncraft_mesh::locate::{Locator, TrackEnd};
use foehncraft_mesh::mesh::Mesh;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// The gas at a point: what the parcel step needs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasState {
    /// Velocity (m/s).
    pub u: [f64; 3],
    /// Density (kg/m³).
    pub rho: f64,
    /// Dynamic viscosity (Pa s).
    pub mu: f64,
}

/// The gas field seen by the parcels (one-way coupling in SP1).
pub trait Gas: Sync {
    /// The gas state at `x`, a point of cell `cell`.
    fn at(&self, cell: usize, x: [f64; 3]) -> GasState;
}

/// What a boundary of the mesh does to a parcel (the case's choice, resolved
/// per boundary of the mesh).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryRule {
    /// The mode.
    pub mode: BoundaryMode,
    /// Restitution of the normal velocity.
    pub e_normal: f64,
    /// Restitution of the tangential velocity.
    pub e_tangential: f64,
}

/// The mesh with its geometry and the point locator.
#[derive(Clone, Copy)]
pub struct Env<'a> {
    /// The mesh.
    pub mesh: &'a Mesh,
    /// Its geometry.
    pub geo: &'a Geometry,
    /// Point location and tracking.
    pub loc: &'a Locator<'a>,
}

/// The parcels in structure-of-arrays form (D47). A parcel stands for
/// `count` identical droplets.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Parcels {
    /// Identity: injector index in the top 16 bits, injection index in the
    /// low 48 (independent of any partition, D117).
    pub id: Vec<u64>,
    /// Position (m).
    pub x: Vec<[f64; 3]>,
    /// Velocity (m/s).
    pub v: Vec<[f64; 3]>,
    /// Droplet diameter (m).
    pub d: Vec<f64>,
    /// Number of droplets the parcel stands for.
    pub count: Vec<f64>,
    /// Time of injection (s).
    pub born: Vec<f64>,
    /// Cell that holds the parcel (rebuilt from the position on a restart).
    #[serde(skip)]
    pub cell: Vec<u32>,
}

impl Parcels {
    /// Number of parcels.
    pub fn len(&self) -> usize {
        self.id.len()
    }

    /// True without parcels.
    pub fn is_empty(&self) -> bool {
        self.id.is_empty()
    }

    fn push(&mut self, p: &One) {
        self.id.push(p.id);
        self.x.push(p.x);
        self.v.push(p.v);
        self.d.push(p.d);
        self.count.push(p.count);
        self.born.push(p.born);
        self.cell.push(p.cell as u32);
    }

    fn get(&self, i: usize) -> One {
        One {
            id: self.id[i],
            x: self.x[i],
            v: self.v[i],
            d: self.d[i],
            count: self.count[i],
            born: self.born[i],
            cell: self.cell[i] as usize,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct One {
    id: u64,
    x: [f64; 3],
    v: [f64; 3],
    d: f64,
    count: f64,
    born: f64,
    cell: usize,
}

/// How the tracking of one parcel ended in one time step.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Fate {
    /// Still in the flow.
    Active,
    /// Stuck on a boundary face.
    Stuck(usize),
    /// Left through a boundary face.
    Escaped(usize),
    /// Not tracked any further: the walk failed or needed too many sub-steps.
    Lost,
}

struct Outcome {
    p: One,
    fate: Fate,
    rebounds: u64,
    substeps: u64,
}

/// Mass and parcel counts, per boundary and in total (kg, parcels).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Tally {
    /// Parcels released.
    pub injected_parcels: u64,
    /// Liquid mass released (kg).
    pub injected_mass: f64,
    /// Liquid mass stuck on each boundary of the mesh (kg).
    pub deposited_mass: Vec<f64>,
    /// Liquid mass that left through each boundary of the mesh (kg).
    pub escaped_mass: Vec<f64>,
    /// Parcels stuck or escaped per boundary of the mesh.
    pub boundary_parcels: Vec<u64>,
    /// Liquid mass of parcels that were lost (kg).
    pub lost_mass: f64,
    /// Parcels lost.
    pub lost_parcels: u64,
    /// Bounces off boundaries.
    pub rebounds: u64,
    /// Sub-steps of all parcels.
    pub substeps: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct InjectorState {
    name: String,
    released: u64,
}

/// The sprays of a run.
pub struct Spray {
    cfg: SpraySection,
    rules: Vec<BoundaryRule>,
    boundary_names: Vec<String>,
    gravity: [f64; 3],
    injectors: Vec<InjectorState>,
    /// The parcels now in the flow.
    pub parcels: Parcels,
    /// The running totals.
    pub tally: Tally,
    /// Time of the last advance (s).
    pub time: f64,
}

/// What `state()` writes and `restore()` reads (restart, D42).
#[derive(Serialize, Deserialize)]
struct Saved {
    format: String,
    version: u32,
    time: f64,
    injectors: Vec<InjectorState>,
    parcels: Parcels,
    tally: Tally,
}

const FORMAT: &str = "foehncraft.spray-state";

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn mass_of(d: f64, count: f64, rho_p: f64) -> f64 {
    count * rho_p * std::f64::consts::PI / 6.0 * d * d * d
}

impl Spray {
    /// Prepares the sprays: `rules` has one entry per boundary of the mesh
    /// (`mesh.boundaries`), `gravity` is the vector of `[gravity]` (zero
    /// without it). Each injector must lie inside the mesh.
    pub fn new(
        cfg: &SpraySection,
        env: &Env,
        rules: Vec<BoundaryRule>,
        gravity: [f64; 3],
    ) -> Result<Self> {
        cfg.validate()?;
        let nb = env.mesh.boundaries.len();
        if rules.len() != nb {
            bail!("{} boundary rules for {nb} boundaries", rules.len());
        }
        for i in &cfg.injectors {
            if env.loc.locate_search(i.position).is_none() {
                bail!(
                    "[[spray.injector]] {}: the position {:?} lies outside the mesh",
                    i.name,
                    i.position
                );
            }
        }
        Ok(Self {
            cfg: cfg.clone(),
            rules,
            boundary_names: env.mesh.boundaries.iter().map(|b| b.name.clone()).collect(),
            gravity,
            injectors: cfg
                .injectors
                .iter()
                .map(|i| InjectorState {
                    name: i.name.clone(),
                    released: 0,
                })
                .collect(),
            parcels: Parcels::default(),
            tally: Tally {
                deposited_mass: vec![0.0; nb],
                escaped_mass: vec![0.0; nb],
                boundary_parcels: vec![0; nb],
                ..Default::default()
            },
            time: 0.0,
        })
    }

    /// The case's section.
    pub fn config(&self) -> &SpraySection {
        &self.cfg
    }

    /// The names of the boundaries, in the order of the tally.
    pub fn boundary_names(&self) -> &[String] {
        &self.boundary_names
    }

    /// Liquid mass in the flow now (kg).
    pub fn mass_in_flow(&self) -> f64 {
        let r = self.cfg.liquid.density;
        (0..self.parcels.len())
            .map(|i| mass_of(self.parcels.d[i], self.parcels.count[i], r))
            .sum()
    }

    /// Advances the sprays from `t0` to `t0 + dt`: releases the parcels that
    /// are due, tracks every parcel through the mesh. The gas field is the
    /// one given (taken as constant over the step). A parcel's path depends
    /// on nothing but itself and the gas, so any number of threads gives
    /// the same bits (D93); the tallies are summed in parcel order.
    pub fn advance(&mut self, env: &Env, gas: &dyn Gas, t0: f64, dt: f64) -> Result<()> {
        let t1 = t0 + dt;
        let mut fresh = Vec::new();
        self.release(env, t1, &mut fresh)?;
        // existing parcels first, in storage order, then the new ones
        let n_old = self.parcels.len();
        for p in &fresh {
            self.parcels.push(p);
        }
        let rho_p = self.cfg.liquid.density;
        let this: &Spray = self;
        let results: Vec<Outcome> = (0..this.parcels.len())
            .into_par_iter()
            .map(|i| {
                let p = this.parcels.get(i);
                let time = if i < n_old { dt } else { t1 - p.born };
                this.track(env, gas, p, time.max(0.0))
            })
            .collect();
        let mut kept = Parcels::default();
        for o in results {
            self.tally.rebounds += o.rebounds;
            self.tally.substeps += o.substeps;
            let m = mass_of(o.p.d, o.p.count, rho_p);
            match o.fate {
                Fate::Active => kept.push(&o.p),
                Fate::Stuck(f) => {
                    if let Some(b) = env.mesh.boundary_of(f) {
                        self.tally.deposited_mass[b] += m;
                        self.tally.boundary_parcels[b] += 1;
                    }
                }
                Fate::Escaped(f) => {
                    if let Some(b) = env.mesh.boundary_of(f) {
                        self.tally.escaped_mass[b] += m;
                        self.tally.boundary_parcels[b] += 1;
                    }
                }
                Fate::Lost => {
                    self.tally.lost_mass += m;
                    self.tally.lost_parcels += 1;
                }
            }
        }
        self.parcels = kept;
        self.time = t1;
        Ok(())
    }

    /// The parcels due before `t1`: parcel k of an injector is born at
    /// start + (k + 1/2)/rate, so the mean mass flow is the one asked for.
    fn release(&mut self, env: &Env, t1: f64, out: &mut Vec<One>) -> Result<()> {
        let rho_p = self.cfg.liquid.density;
        for (ii, inj) in self.cfg.injectors.iter().enumerate() {
            let rate = inj.parcels_per_second;
            let end = inj.end.unwrap_or(f64::INFINITY);
            let axis = unit(inj.direction);
            let (e1, e2) = frame(axis);
            let cos_a = inj.half_angle.to_radians().cos();
            let m_parcel = inj.mass_flow / rate;
            let d = inj.diameter;
            let count = m_parcel / (rho_p * std::f64::consts::PI / 6.0 * d * d * d);
            loop {
                let k = self.injectors[ii].released;
                let born = inj.start + (k as f64 + 0.5) / rate;
                if born > t1 || born > end {
                    break;
                }
                self.injectors[ii].released += 1;
                let id = ((ii as u64) << 48) | k;
                let dir = direction(id, axis, e1, e2, cos_a);
                let Some(cell) = env.loc.locate_search(inj.position) else {
                    bail!(
                        "[[spray.injector]] {}: the position lies outside the mesh",
                        inj.name
                    );
                };
                out.push(One {
                    id,
                    x: inj.position,
                    v: [dir[0] * inj.speed, dir[1] * inj.speed, dir[2] * inj.speed],
                    d,
                    count,
                    born,
                    cell,
                });
                self.tally.injected_parcels += 1;
                self.tally.injected_mass += m_parcel;
            }
        }
        Ok(())
    }

    fn rule_of(&self, env: &Env, face: usize) -> BoundaryRule {
        match env.mesh.boundary_of(face) {
            Some(b) => self.rules[b],
            None => BoundaryRule {
                mode: BoundaryMode::Escape,
                e_normal: 1.0,
                e_tangential: 1.0,
            },
        }
    }

    /// Tracks one parcel for `time` seconds.
    fn track(&self, env: &Env, gas: &dyn Gas, mut p: One, time: f64) -> Outcome {
        let rho_p = self.cfg.liquid.density;
        let mut left = time;
        let mut rebounds = 0u64;
        let mut substeps = 0u64;
        let mut fate = Fate::Active;
        let tiny = 1e-14 * time.max(1e-300);
        while left > tiny {
            substeps += 1;
            if substeps as usize > self.cfg.max_substeps {
                fate = Fate::Lost;
                break;
            }
            let gs = gas.at(p.cell, p.x);
            let h = env.geo.cell_volume[p.cell].cbrt();
            let input = |dt: f64| ParcelInput {
                v: p.v,
                u: gs.u,
                rho_g: gs.rho,
                mu_g: gs.mu,
                d: p.d,
                rho_p,
                g: self.gravity,
                dt,
            };
            let mut sub = left;
            let mut s = parcel_step(&input(sub));
            let len = norm(s.dx);
            if len > self.cfg.max_courant * h {
                sub = left * self.cfg.max_courant * h / len;
                s = parcel_step(&input(sub));
            }
            let target = [p.x[0] + s.dx[0], p.x[1] + s.dx[1], p.x[2] + s.dx[2]];
            let start = env.loc.position_in(p.cell, p.x).or_else(|| {
                env.loc
                    .locate_search(p.x)
                    .and_then(|c| env.loc.position_in(c, p.x))
            });
            let Some(start) = start else {
                fate = Fate::Lost;
                break;
            };
            let tr = env.loc.track(start, target);
            match tr.end {
                TrackEnd::Reached => {
                    p.x = target;
                    p.v = s.v;
                    p.cell = tr.position.cell;
                    left -= sub;
                }
                TrackEnd::Boundary { face } => {
                    // the time to the face: the chord fraction of the
                    // sub-step (the sub-step is a fraction of a cell)
                    let thit = tr.fraction * sub;
                    let sh = parcel_step(&input(thit));
                    p.x = tr.position.x;
                    p.v = sh.v;
                    p.cell = tr.position.cell;
                    left -= thit;
                    let rule = self.rule_of(env, face);
                    match rule.mode {
                        BoundaryMode::Stick => {
                            fate = Fate::Stuck(face);
                            break;
                        }
                        BoundaryMode::Escape => {
                            fate = Fate::Escaped(face);
                            break;
                        }
                        BoundaryMode::Rebound | BoundaryMode::Mirror => {
                            let a = env.geo.face_area[face];
                            let n = unit(a);
                            let vn = p.v[0] * n[0] + p.v[1] * n[1] + p.v[2] * n[2];
                            let (en, et) = if rule.mode == BoundaryMode::Mirror {
                                (1.0, 1.0)
                            } else {
                                (rule.e_normal, rule.e_tangential)
                            };
                            for k in 0..3 {
                                let vt = p.v[k] - vn * n[k];
                                p.v[k] = et * vt - en * vn * n[k];
                            }
                            // a hair inside the cell, so the next track
                            // starts in it
                            let eps = 1e-6 * h;
                            for k in 0..3 {
                                p.x[k] -= eps * n[k];
                            }
                            rebounds += 1;
                        }
                    }
                }
                TrackEnd::Processor { .. } | TrackEnd::Lost => {
                    // the walk failed: one search for the target
                    match env.loc.locate_search(target) {
                        Some(c) => {
                            p.x = target;
                            p.v = s.v;
                            p.cell = c;
                            left -= sub;
                        }
                        None => {
                            fate = Fate::Lost;
                            break;
                        }
                    }
                }
            }
        }
        Outcome {
            p,
            fate,
            rebounds,
            substeps,
        }
    }

    /// The state for a restart file: the parcels, the counters of the
    /// injectors and the tallies (JSON, `foehncraft.spray-state`).
    pub fn state(&self) -> Result<Vec<u8>> {
        Ok(serde_json::to_vec(&Saved {
            format: FORMAT.into(),
            version: 1,
            time: self.time,
            injectors: self.injectors.clone(),
            parcels: self.parcels.clone(),
            tally: self.tally.clone(),
        })?)
    }

    /// Restores the state of [`state`](Self::state); the cells of the
    /// parcels are found again from their positions (the mesh may be
    /// numbered differently).
    pub fn restore(&mut self, env: &Env, bytes: &[u8]) -> Result<()> {
        let s: Saved = serde_json::from_slice(bytes)?;
        if s.format != FORMAT || s.version != 1 {
            bail!("the restart holds a spray state of format {} version {}", s.format, s.version);
        }
        let names: Vec<&str> = s.injectors.iter().map(|i| i.name.as_str()).collect();
        let mine: Vec<&str> = self.injectors.iter().map(|i| i.name.as_str()).collect();
        if names != mine {
            bail!("the restart's injectors {names:?} differ from the case's {mine:?}");
        }
        if s.tally.deposited_mass.len() != self.boundary_names.len() {
            bail!("the restart's spray tally has another number of boundaries");
        }
        let mut parcels = s.parcels;
        parcels.cell = Vec::with_capacity(parcels.len());
        for x in &parcels.x {
            match env.loc.locate_search(*x) {
                Some(c) => parcels.cell.push(c as u32),
                None => bail!("a parcel of the restart lies outside the mesh at {x:?}"),
            }
        }
        self.injectors = s.injectors;
        self.parcels = parcels;
        self.tally = s.tally;
        self.time = s.time;
        Ok(())
    }

    /// Writes the parcels to a VTKHDF `PolyData` file (D75 attributes).
    #[cfg(feature = "hdf5")]
    pub fn write(&self, path: &std::path::Path, case: &str, step: u64) -> Result<()> {
        use foehncraft_mesh::vtkhdf::{Attr, PointField, write_points};
        let n = self.parcels.len();
        let vel: Vec<f64> = self.parcels.v.iter().flatten().copied().collect();
        let id: Vec<f64> = self.parcels.id.iter().map(|&i| i as f64).collect();
        let age: Vec<f64> = self.parcels.born.iter().map(|b| self.time - b).collect();
        let fields = [
            PointField {
                name: "velocity",
                values: &vel,
                components: 3,
            },
            PointField {
                name: "diameter",
                values: &self.parcels.d,
                components: 1,
            },
            PointField {
                name: "droplets",
                values: &self.parcels.count,
                components: 1,
            },
            PointField {
                name: "age",
                values: &age,
                components: 1,
            },
            PointField {
                name: "parcel_id",
                values: &id,
                components: 1,
            },
        ];
        let attrs = [
            ("format", Attr::Str("foehncraft.parcels")),
            ("format_version", Attr::U64(1)),
            ("case", Attr::Str(case)),
            ("step", Attr::U64(step)),
            ("time", Attr::F64(self.time)),
            ("liquid", Attr::Str(&self.cfg.liquid.name)),
            ("units", Attr::Str("m, m/s, s, droplets")),
            ("parcels", Attr::U64(n as u64)),
        ];
        write_points(path, &self.parcels.x, &fields, &attrs)
    }
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Two unit vectors orthogonal to `a` and to each other.
fn frame(a: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    // the axis least aligned with a
    let t = if a[0].abs() <= a[1].abs() && a[0].abs() <= a[2].abs() {
        [1.0, 0.0, 0.0]
    } else if a[1].abs() <= a[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let c = [
        a[1] * t[2] - a[2] * t[1],
        a[2] * t[0] - a[0] * t[2],
        a[0] * t[1] - a[1] * t[0],
    ];
    let e1 = unit(c);
    let e2 = [
        a[1] * e1[2] - a[2] * e1[1],
        a[2] * e1[0] - a[0] * e1[2],
        a[0] * e1[1] - a[1] * e1[0],
    ];
    (e1, e2)
}

/// A direction in the cone about `axis`, uniform in solid angle, from the
/// counter-based generator keyed by the parcel identity (D121).
fn direction(id: u64, axis: [f64; 3], e1: [f64; 3], e2: [f64; 3], cos_a: f64) -> [f64; 3] {
    if cos_a >= 1.0 {
        return axis;
    }
    let u1 = rng::uniform_keyed([stream::SPRAY_INJECTION, id, 0, 0]);
    let u2 = rng::uniform_keyed([stream::SPRAY_INJECTION, id, 0, 1]);
    let ct = 1.0 - u1 * (1.0 - cos_a);
    let st = (1.0 - ct * ct).max(0.0).sqrt();
    let phi = 2.0 * std::f64::consts::PI * u2;
    let (sp, cp) = phi.sin_cos();
    let mut d = [0.0; 3];
    for k in 0..3 {
        d[k] = ct * axis[k] + st * (cp * e1[k] + sp * e2[k]);
    }
    d
}

impl Injector {
    /// Droplets per parcel at the injector.
    pub fn droplets_per_parcel(&self, rho_p: f64) -> f64 {
        let m = self.mass_flow / self.parcels_per_second;
        m / (rho_p * std::f64::consts::PI / 6.0 * self.diameter.powi(3))
    }
}
