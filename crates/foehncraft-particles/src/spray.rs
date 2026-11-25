//! The parcels and their tracking through the mesh.

use crate::config::{BoundaryMode, Dispersion, Injector, Penetration, SpraySection};
use crate::coupling::{CellSources, Record, StepBooks};
use crate::datum::{Crossing, Datum, DatumTally, Moments, Shape, Surface, plane_crossing};
use crate::equilibrium::{SurfacePoint, SurfaceTable};
use anyhow::{Result, bail};
use foehncraft_core::Real;
use foehncraft_core::rng::{self, stream};
use foehncraft_kernels::evaporation::{self as evap, EvapInput, EvapStep};
use foehncraft_kernels::particle::{HeatInput, ParcelInput, eddy_times, heat_step, parcel_step};
use foehncraft_mesh::geometry::Geometry;
use foehncraft_mesh::locate::{Locator, Position, TrackEnd, WalkEnd, step_limit};
use foehncraft_mesh::mesh::Mesh;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

mod steady;
pub use steady::{MomentumFlows, StreamBooks, StreamPass, Streams};

/// The gas at a point: what the parcel step needs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasState {
    /// Velocity (m/s).
    pub u: [f64; 3],
    /// Density (kg/m³).
    pub rho: f64,
    /// Dynamic viscosity (Pa s).
    pub mu: f64,
    /// Temperature (K); NaN without an energy equation (no heat exchange).
    pub t: f64,
    /// Specific heat at constant pressure (J/(kg K)).
    pub cp: f64,
    /// Thermal conductivity (W/(m K)); zero: no heat exchange.
    pub k: f64,
    /// Turbulent kinetic energy (m²/s²); zero: no turbulent dispersion.
    pub tke: f64,
    /// Its dissipation rate ε (m²/s³).
    pub eps: f64,
    /// Mass fraction of the liquid's vapour in the gas (SP5); NaN without
    /// the vapour species (no evaporation).
    pub y_v: f64,
    /// Thermodynamic pressure p₀ (Pa).
    pub p0: f64,
    /// Molar mass of the gas without the vapour (kg/mol).
    pub m_gas: f64,
    /// Diffusivity of the vapour in the gas (m²/s).
    pub diff: f64,
}

impl GasState {
    /// A gas without temperature: the parcels keep theirs (no heat
    /// exchange).
    pub fn isothermal(u: [f64; 3], rho: f64, mu: f64) -> Self {
        Self {
            u,
            rho,
            mu,
            t: f64::NAN,
            cp: 1006.0,
            k: 0.0,
            tke: 0.0,
            eps: 0.0,
            y_v: f64::NAN,
            p0: 101_325.0,
            m_gas: 0.028_965,
            diff: 0.0,
        }
    }

    /// True when the gas carries the vapour of the liquid (evaporation).
    pub fn has_vapour(&self) -> bool {
        self.y_v.is_finite() && self.diff > 0.0 && self.p0 > 0.0 && self.has_heat()
    }

    /// True when the gas exchanges heat with the parcels.
    pub fn has_heat(&self) -> bool {
        self.t.is_finite() && self.k > 0.0 && self.cp > 0.0
    }
}

/// The vapour of the liquid (SP4, SP5): its molar mass and the NASA-7
/// polynomials of the species data, for the latent heat L(T) = h_v(T) −
/// c_l T + const (Kirchhoff's law with the latent heat of the case at one
/// temperature) and the specific heat of the vapour.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vapour {
    /// The species of the data file.
    pub species: foehncraft_kernels::species::Species,
}

impl Vapour {
    /// Water vapour with the NASA-7 data of GRI-Mech 3.0 (the `H2O` of the
    /// species data of the flow solver), for tests and tools without the
    /// data file.
    pub fn water() -> Self {
        Self {
            species: foehncraft_kernels::species::Species {
                molar_mass: 0.018_015_28,
                t_range: [200.0, 1000.0, 3500.0],
                low: [
                    4.19864056,
                    -0.0020364341,
                    6.52040211e-06,
                    -5.48797062e-09,
                    1.77197817e-12,
                    -30293.7267,
                    -0.849032208,
                ],
                high: [
                    3.03399249,
                    0.00217691804,
                    -1.64072518e-07,
                    -9.7041987e-11,
                    1.68200992e-14,
                    -30004.2971,
                    4.9667701,
                ],
                well_depth: 572.4,
                diameter: 2.605,
                dipole: 1.844,
                geometry: 2,
            },
        }
    }

    /// Enthalpy of the vapour (J/kg, the reference of the species data).
    pub fn h(&self, t: f64) -> f64 {
        foehncraft_kernels::species::R * t * self.species.h_over_rt(t) / self.species.molar_mass
    }

    /// Specific heat of the vapour (J/(kg K)).
    pub fn cp(&self, t: f64) -> f64 {
        self.species.cp(t)
    }
}

/// The eddy a parcel is in (eddy-lifetime dispersion, D243): the
/// fluctuation of the gas velocity it sees, the time left until it enters
/// the next eddy, and the number of eddies it has entered (the event word
/// of its random numbers, stream `SPRAY_DISPERSION`, D121).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct Eddy {
    /// Fluctuation of the gas velocity (m/s).
    pub u: [f64; 3],
    /// Time left in the eddy (s); zero or less: a new eddy at the next
    /// sub-step.
    pub left: f64,
    /// Eddies entered so far.
    pub n: u64,
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
    /// Droplet temperature (K); empty in restart files written before
    /// SP2b (filled with the injection temperature on restore).
    #[serde(default)]
    pub temp: Vec<f64>,
    /// The eddy each parcel is in (turbulent dispersion, D243); empty in
    /// restart files written before SP2c (every parcel starts a new eddy).
    #[serde(default)]
    pub eddy: Vec<Eddy>,
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
        self.temp.push(p.temp);
        self.eddy.push(p.eddy);
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
            temp: self.temp[i],
            eddy: self.eddy[i],
            cell: self.cell[i] as usize,
            face: None,
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
    temp: f64,
    eddy: Eddy,
    cell: usize,
    /// Where the walk stands when the parcel was stopped on an interior
    /// face (two-way tracking): the tetrahedron of the next cell.
    face: Option<Position>,
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
    /// Evaporated completely (SP5).
    Evaporated,
}

struct Outcome {
    p: One,
    fate: Fate,
    rebounds: u64,
    substeps: u64,
    cross: Vec<Crossing>,
    /// The exchange records (two-way coupling), in sub-step order.
    records: Vec<Record>,
    /// Impulse of the weight over the step (kg m/s).
    weight: [f64; 3],
    /// Impulse of the walls (rebounds, mirrors).
    walls: [f64; 3],
    /// Time the parcel was followed (s).
    used: f64,
    /// Liquid the parcel lost to vapour over the step (kg, its own change;
    /// negative: condensation).
    vapour: f64,
    /// The path (steady streams with [`TrackOpts::paths`]): the state at
    /// the start of every sub-step and at the end.
    path: Vec<PathPoint>,
}

/// A point of the path of a stream (steady mode, D260).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PathPoint {
    /// Time since the start of the stream (s).
    pub(crate) t: f64,
    pub(crate) x: [f64; 3],
    pub(crate) v: [f64; 3],
    pub(crate) temp: f64,
}

/// How [`Spray::track`] follows a parcel.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct TrackOpts {
    /// Keep the path.
    pub(crate) paths: bool,
    /// A steady stream (D260): the slope of a record is the local one of
    /// its visit of the cell, m (1 − Π e) over the sub-steps of the visit
    /// ([`crate::coupling::set_slopes_local`]).
    pub(crate) steady: bool,
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
    /// Liquid mass that evaporated (kg; net of condensation, SP5).
    #[serde(default)]
    pub evaporated_mass: f64,
    /// Parcels that evaporated completely.
    #[serde(default)]
    pub evaporated_parcels: u64,
    /// Parcels at or above the critical mixing temperature of their
    /// surface at the end of the last step (Peng-Robinson surface, D314:
    /// no latent heat, diffusive mixing); a health number, not a total.
    #[serde(default)]
    pub critical_parcels: u64,
    /// Their liquid mass (kg).
    #[serde(default)]
    pub critical_mass: f64,
}

/// The parcels released in a step, booked when the pass is committed.
#[derive(Default)]
struct Release {
    fresh: Vec<One>,
    masses: Vec<f64>,
    released: Vec<u64>,
}

/// One parcel pass (D247): the parcels at the end of the step and their
/// exchange records, not yet committed ([`Spray::pass`], [`Spray::commit`]).
pub struct Pass {
    release: Release,
    outcomes: Vec<Outcome>,
    /// The sorted per-cell sums of the exchange records (two-way coupling;
    /// empty one way).
    pub sources: CellSources,
    /// End of the step (s).
    pub t1: f64,
    /// Length of the step (s).
    pub dt: f64,
}

impl Pass {
    /// The exchange records of the pass, parcel by parcel.
    pub fn records(&self) -> impl Iterator<Item = &Record> {
        self.outcomes.iter().flat_map(|o| o.records.iter())
    }

    /// Sub-steps of all parcels in the pass.
    pub fn substeps(&self) -> u64 {
        self.outcomes.iter().map(|o| o.substeps).sum()
    }
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
    /// The datums the parcels report to (D262).
    datums: Vec<Datum>,
    /// Running sums of each datum (aligned with the datums; volumes keep
    /// empty ones).
    pub datum_tally: Vec<DatumTally>,
    /// A two-dimensional case in the x-y plane: the eddies have no z
    /// component (D243).
    planar: bool,
    /// The vapour of the liquid (evaporation, SP5).
    vapour: Option<Vapour>,
    /// The Peng-Robinson surface states over the droplet temperature at
    /// p₀ (`equilibrium = "peng-robinson"`, D311), from
    /// [`Spray::set_surface_pressure`].
    surface: Option<SurfaceTable>,
    /// The vapour penetration per injector of the last gas the flow gave
    /// (`[spray.ecn]`, D315), from [`Spray::set_vapour_penetration`].
    vapour_penetration: Vec<Option<f64>>,
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
    /// The sums of the datums (SP1b); absent in files of SP1.
    #[serde(default)]
    datums: Vec<DatumTally>,
}

const FORMAT: &str = "foehncraft.spray-state";

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn mass_of(d: f64, count: f64, rho_p: f64) -> f64 {
    count * rho_p * std::f64::consts::PI / 6.0 * d * d * d
}

/// ∫ m(t) dt over a sub-step of length `dt` of a droplet of mass `m0` and
/// diameter `d0` that follows the d² law d² = d₀² − K t to `d1`:
/// m₀ (d₀²/K) (2/5) (1 − (d₁/d₀)⁵).
fn mass_integral(m0: f64, d0: f64, d1: f64, k: f64, dt: f64) -> f64 {
    let s = k * dt / (d0 * d0);
    if s.abs() < 1e-8 || d0 <= 0.0 {
        return m0 * dt;
    }
    let r = d1 / d0;
    let r5 = r * r * r * r * r;
    m0 * d0 * d0 / k * 0.4 * (1.0 - r5)
}

/// A standard normal number from two keyed uniforms (Box and Muller 1958).
fn normal(key: [u64; 4]) -> f64 {
    let u1 = rng::uniform_keyed(key);
    let u2 = rng::uniform_keyed([key[0], key[1], key[2], key[3] + 1]);
    // 1 − u1 lies in (0, 1]
    (-2.0 * (1.0 - u1).ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

impl Spray {
    /// Marks a two-dimensional case in the x-y plane: the eddies of the
    /// dispersion get no z component.
    pub fn set_planar(&mut self, planar: bool) {
        self.planar = planar;
    }

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
            planar: false,
            vapour: None,
            surface: None,
            vapour_penetration: Vec::new(),
            datums: Vec::new(),
            datum_tally: Vec::new(),
            tally: Tally {
                deposited_mass: vec![0.0; nb],
                escaped_mass: vec![0.0; nb],
                boundary_parcels: vec![0; nb],
                ..Default::default()
            },
            time: 0.0,
        })
    }

    /// Sets the vapour of the liquid (the species data of
    /// `[spray.liquid] vapour`): needed for the evaporation (SP5).
    pub fn set_vapour(&mut self, v: Option<Vapour>) {
        self.vapour = v;
    }

    /// True when the droplets evaporate (`[spray.evaporation]` and the
    /// vapour species set).
    pub fn evaporates(&self) -> bool {
        self.cfg.evaporation.on() && self.vapour.is_some()
    }

    /// Latent heat of the liquid at `t` (J/kg): L(T_ref) + h_v(T) −
    /// h_v(T_ref) − c_l (T − T_ref) (Kirchhoff's law with the vapour
    /// enthalpy of the species data and the constant c_l of the liquid);
    /// with the Peng-Robinson surface times L_PR(T, p₀)/L_PR,sat(T) (D312),
    /// 0 at and above the critical mixing state (D314). The gas receives
    /// the vapour at c_l T + L, so the books close either way.
    pub fn latent_heat(&self, t: f64) -> f64 {
        let l = &self.cfg.liquid;
        let low = match &self.vapour {
            Some(v) => {
                l.latent_heat + v.h(t)
                    - v.h(l.latent_heat_temperature)
                    - l.specific_heat * (t - l.latent_heat_temperature)
            }
            None => l.latent_heat,
        };
        match &self.surface {
            Some(tab) => low * tab.at(t).latent_ratio,
            None => low,
        }
    }

    /// The Peng-Robinson surface state at the droplet temperature `t`, or
    /// `None` with Raoult's law.
    pub fn surface_state(&self, t: f64) -> Option<SurfacePoint> {
        self.surface.as_ref().map(|tab| tab.at(t))
    }

    /// The Peng-Robinson surface table, if built.
    pub fn surface_table(&self) -> Option<&SurfaceTable> {
        self.surface.as_ref()
    }

    /// Builds the Peng-Robinson surface table for the thermodynamic
    /// pressure `p0` (Pa) of the run, when `[spray.evaporation]
    /// equilibrium = "peng-robinson"` (D311); nothing with Raoult's law.
    /// p₀ is one number per run here (a closed volume is refused with
    /// evaporation), so the table is a pure function of the case and a
    /// restart rebuilds the same bits.
    pub fn set_surface_pressure(&mut self, p0: f64) -> Result<()> {
        use crate::config::Equilibrium;
        use foehncraft_kernels::realfluid::{self as rf, MAX};
        let ev = &self.cfg.evaporation;
        if ev.equilibrium != Equilibrium::PengRobinson {
            self.surface = None;
            return Ok(());
        }
        let pr = &ev.peng_robinson;
        let key = "[spray.evaporation.peng_robinson]";
        let consts = [
            pr.critical_temperature,
            pr.critical_pressure,
            pr.acentric_factor,
            pr.molar_mass,
        ];
        let (fuel, defaults): (rf::Component, Option<rf::fuels::Fuel>) = match consts {
            [Some(tc), Some(pc), Some(w), Some(m)] => {
                if !(tc > 0.0 && pc > 0.0 && m > 0.0 && w.is_finite()) {
                    bail!(
                        "{key}: critical_temperature {tc} K, critical_pressure {pc} Pa, acentric_factor {w}, molar_mass {m} kg/mol"
                    );
                }
                (rf::Component::new(tc, pc, w, m), None)
            }
            [None, None, None, None] => {
                let name = pr.fuel.as_deref().unwrap_or(&self.cfg.liquid.name);
                let f = rf::fuels::by_name(name).ok_or_else(|| {
                    anyhow::anyhow!(
                        "{key}: no high-pressure data for the liquid \"{name}\"; known: \
                             \"n-dodecane\"; or give critical_temperature, critical_pressure, \
                             acentric_factor and molar_mass"
                    )
                })?;
                (f.component, Some(f))
            }
            _ => bail!(
                "{key}: give all four of critical_temperature, critical_pressure, \
                     acentric_factor and molar_mass, or none (then the data of `fuel`)"
            ),
        };
        if pr.ambient.is_empty() || pr.ambient.len() + 1 > MAX {
            bail!("{key}: ambient needs 1 to {} gases", MAX - 1);
        }
        let mut comps = vec![fuel];
        let mut r = [0.0; MAX];
        let mut sum = 0.0;
        for (j, (name, &x)) in pr.ambient.iter().enumerate() {
            let g = rf::gases::by_name(name).ok_or_else(|| {
                anyhow::anyhow!(
                    "{key}: ambient gas \"{name}\" unknown; known: N2, O2, CO2, H2O, AR"
                )
            })?;
            if !(x > 0.0 && x.is_finite()) {
                bail!("{key}: ambient {name} = {x}: mole fractions are positive");
            }
            comps.push(g);
            r[j + 1] = x;
            sum += x;
        }
        for v in r.iter_mut().skip(1) {
            *v /= sum;
        }
        for name in pr.k_ij.keys() {
            if !pr.ambient.contains_key(name) {
                bail!("{key}: k_ij.{name} names no gas of ambient");
            }
        }
        let mut mix = rf::Mixture::new(&comps);
        for (j, name) in pr.ambient.keys().enumerate() {
            let k = match pr.k_ij.get(name) {
                Some(&k) => k,
                None => defaults.map_or(0.0, |f| f.k_with(name)),
            };
            mix.set_k(0, j + 1, k);
        }
        if !(pr.table_step > 0.0 && pr.table_step <= 10.0) {
            bail!("{key}: table_step {} K: between 0 and 10 K", pr.table_step);
        }
        let t_lo = (self.cfg.liquid.temperature - 60.0).clamp(150.0, 250.0);
        let t_hi = fuel.tc + 50.0;
        let tab = SurfaceTable::build(&mix, &r, p0, t_lo, t_hi, pr.table_step)?;
        self.surface = Some(tab);
        Ok(())
    }

    /// Enthalpy of the vapour at `t` in the reference of the parcels (J/kg):
    /// c_l T + L(T), so the liquid at 0 K has none.
    pub fn vapour_enthalpy(&self, t: f64) -> f64 {
        self.cfg.liquid.specific_heat * t + self.latent_heat(t)
    }

    /// Specific heat of the vapour at `t` (J/(kg K)); 0 without one.
    pub fn vapour_cp(&self, t: f64) -> f64 {
        self.vapour.map_or(0.0, |v| v.cp(t))
    }

    /// Saturation pressure of the liquid at `t` (Pa); with the
    /// Peng-Robinson surface the partial pressure y_f p₀ of the fuel in the
    /// equilibrium vapour at the surface (D311), which the Spalding numbers
    /// then take unchanged.
    pub fn saturation_pressure(&self, t: f64) -> f64 {
        if let Some(tab) = &self.surface {
            return tab.at(t).y_fuel * tab.p0;
        }
        match self.cfg.liquid.saturation {
            crate::config::Saturation::Water => evap::water::saturation_pressure(t),
        }
    }

    /// The evaporation of a droplet of diameter `d` and temperature `tp`
    /// over `dt` in the gas `gs` at the particle Reynolds number `re`, or
    /// `None` when it does not evaporate (no model, no vapour in the gas).
    pub fn evaporation(
        &self,
        gs: &GasState,
        d: f64,
        tp: f64,
        re: f64,
        dt: f64,
    ) -> Option<EvapStep> {
        let v = self.vapour.as_ref()?;
        if !self.cfg.evaporation.on() || !gs.has_vapour() || d <= 0.0 {
            return None;
        }
        let l = &self.cfg.liquid;
        Some(evap::evaporation_step(&EvapInput {
            tp,
            tg: gs.t,
            yg: gs.y_v.clamp(0.0, 1.0 - 1e-12),
            p0: gs.p0,
            m_vapour: v.species.molar_mass,
            m_gas: gs.m_gas,
            psat: self.saturation_pressure(tp),
            latent: self.latent_heat(tp),
            re,
            pr: gs.mu * gs.cp / gs.k,
            sc: gs.mu / (gs.rho * gs.diff),
            rho_g: gs.rho,
            diff: gs.diff,
            k_g: gs.k,
            cp_g: gs.cp,
            cp_v: v.cp(0.5 * (tp + gs.t)),
            d,
            rho_l: l.density,
            cp_l: l.specific_heat,
            dt,
            model: self.cfg.evaporation.code(),
        }))
    }

    /// Sets the datums the parcels report to; their sums start at the
    /// spray's time.
    pub fn set_datums(&mut self, datums: Vec<Datum>) {
        self.datum_tally = datums
            .iter()
            .map(|d| DatumTally::new(&d.name, &d.hist, self.time))
            .collect();
        self.datums = datums;
    }

    /// The datums.
    pub fn datums(&self) -> &[Datum] {
        &self.datums
    }

    /// The liquid mass of a parcel.
    pub fn parcel_mass(&self, i: usize) -> f64 {
        mass_of(
            self.parcels.d[i],
            self.parcels.count[i],
            self.cfg.liquid.density,
        )
    }

    /// The moment sums of the parcels in the flow, of all (`None`) or of
    /// those inside volume datum `datum`; summed in the order of the parcel
    /// identities, so the bits do not depend on the storage order.
    /// `zone` is the cell zone of each cell.
    pub fn moments(&self, datum: Option<usize>, zone: &[u32]) -> Moments {
        let spec = datum.map_or_else(crate::datum::HistSpec::default, |i| self.datums[i].hist);
        let mut m = Moments::new(&spec);
        for i in self.by_id() {
            if let Some(di) = datum {
                let Shape::Volume(v) = &self.datums[di].shape else {
                    continue;
                };
                if !v.contains(self.parcels.x[i], self.parcels.cell[i] as usize, zone) {
                    continue;
                }
            }
            m.add(
                &spec,
                self.parcels.d[i],
                self.parcels.count[i],
                self.parcel_mass(i),
                self.parcels.v[i],
            );
        }
        m
    }

    /// Parcel indices in the order of their identities.
    fn by_id(&self) -> Vec<usize> {
        let mut ix: Vec<usize> = (0..self.parcels.len()).collect();
        ix.sort_by_key(|&i| self.parcels.id[i]);
        ix
    }

    /// Per injector: parcels and liquid in the flow, the liquid
    /// penetration and the farthest parcel (D265), with the definitions of
    /// `[spray] penetration` and `penetration_fraction`.
    pub fn injector_stats(&self) -> Vec<InjectorStats> {
        let masses: Vec<f64> = (0..self.parcels.len())
            .map(|i| self.parcel_mass(i))
            .collect();
        let released: Vec<u64> = self.injectors.iter().map(|i| i.released).collect();
        let mut st = injector_stats(
            &self.cfg,
            &self.parcels.id,
            &self.parcels.x,
            &masses,
            Some(&released),
        );
        for (s, v) in st.iter_mut().zip(&self.vapour_penetration) {
            s.vapour_penetration = *v;
        }
        st
    }

    /// Sets the vapour penetration of each injector from the fuel mass
    /// fraction `y` of the gas cells with the centroids `centre` (D315,
    /// `[spray.ecn]`): the farthest distance along the injector's axis,
    /// from the nozzle, of a cell centre where y ≥ `vapour_fraction`; 0
    /// when no cell reaches it. With several injectors the vapour of all
    /// of them counts for each. Nothing without `[spray.ecn]`.
    pub fn set_vapour_penetration(&mut self, y: &[f64], centre: &[[f64; 3]]) {
        let Some(e) = &self.cfg.ecn else {
            self.vapour_penetration.clear();
            return;
        };
        self.vapour_penetration = self
            .cfg
            .injectors
            .iter()
            .map(|inj| {
                let a = unit(inj.direction);
                let mut far = 0.0f64;
                for (c, &yc) in y.iter().enumerate().take(centre.len()) {
                    if yc >= e.vapour_fraction {
                        let r = [
                            centre[c][0] - inj.position[0],
                            centre[c][1] - inj.position[1],
                            centre[c][2] - inj.position[2],
                        ];
                        far = far.max(r[0] * a[0] + r[1] * a[1] + r[2] * a[2]);
                    }
                }
                Some(far)
            })
            .collect();
    }

    /// Centre of mass of the liquid in the flow (m); zero without parcels.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let mut c = [0.0; 3];
        let mut m = 0.0;
        for i in self.by_id() {
            let mi = self.parcel_mass(i);
            m += mi;
            for (k, ck) in c.iter_mut().enumerate() {
                *ck += mi * self.parcels.x[i][k];
            }
        }
        if m > 0.0 {
            for ck in &mut c {
                *ck /= m;
            }
        }
        c
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
        let pass = self.pass(env, gas, t0, dt)?;
        self.commit(env, pass, None);
        Ok(())
    }

    /// One parcel pass from `t0` to `t0 + dt` (D247): the parcels of the
    /// old time level and those released in the step, tracked through the
    /// frozen gas `gas`. The spray is not changed: a pass can be repeated
    /// (D252) and gives the same parcels for the same gas, whatever the
    /// order of the parcels or the number of threads. With two-way coupling
    /// the pass carries the exchange records and their sorted per-cell sums
    /// (D249).
    pub fn pass(&self, env: &Env, gas: &dyn Gas, t0: f64, dt: f64) -> Result<Pass> {
        let t1 = t0 + dt;
        let release = self.release(env, t1)?;
        let n_old = self.parcels.len();
        let fresh = &release.fresh;
        let outcomes: Vec<Outcome> = (0..n_old + fresh.len())
            .into_par_iter()
            .map(|i| {
                let (p, time) = if i < n_old {
                    (self.parcels.get(i), dt)
                } else {
                    let p = fresh[i - n_old];
                    (p, t1 - p.born)
                };
                self.track(env, gas, p, time.max(0.0), t1, TrackOpts::default())
            })
            .collect();
        let sources = if self.cfg.coupling.two_way() {
            let all: Vec<Record> = outcomes
                .iter()
                .flat_map(|o| o.records.iter().copied())
                .collect();
            CellSources::from_records(&all)
        } else {
            CellSources::default()
        };
        Ok(Pass {
            release,
            outcomes,
            sources,
            t1,
            dt,
        })
    }

    /// Commits pass `pass` (D250, D252): the closing step, when `du` gives
    /// per cell the change u_new − u* of the gas velocity and T_new − T*
    /// of the gas temperature against the frozen state of the pass (every
    /// record r gives its parcel slope_r (u_new − u*) of momentum and
    /// heat_slope_r (T_new − T*) of enthalpy, the opposite of what the gas
    /// received), then the books, the release and the fates. Returns the
    /// momentum and enthalpy books of the step.
    pub fn commit(
        &mut self,
        env: &Env,
        pass: Pass,
        du: Option<&dyn Fn(u32) -> crate::coupling::Delta>,
    ) -> StepBooks {
        let rho_p = self.cfg.liquid.density;
        let cp_p = self.cfg.liquid.specific_heat;
        let mut books = StepBooks {
            start: self.momentum(),
            ..Default::default()
        };
        books.enthalpy[0] = self.enthalpy();
        books.mass[0] = self.mass_in_flow_by_id();
        let Pass {
            release,
            mut outcomes,
            sources,
            t1,
            ..
        } = pass;
        // the release of the step (counters, mass book)
        for (ii, &r) in release.released.iter().enumerate() {
            self.injectors[ii].released = r;
        }
        for (p, &m) in release.fresh.iter().zip(&release.masses) {
            self.tally.injected_parcels += 1;
            self.tally.injected_mass += m;
            let mp = mass_of(p.d, p.count, rho_p);
            for k in 0..3 {
                books.injected[k] += mp * p.v[k];
            }
            books.enthalpy[3] += mp * cp_p * p.temp;
            books.mass[3] += mp;
        }
        // the closing step: the parcel receives what the gas gave up
        // (momentum, heat and, with evaporation, the vapour the gas did not
        // take with the enthalpy it carries, SP5)
        if let Some(du) = du {
            for o in &mut outcomes {
                if o.records.is_empty() {
                    continue;
                }
                let m = mass_of(o.p.d, o.p.count, rho_p);
                if m <= 0.0 {
                    // evaporated completely: no slopes (D250), nothing to
                    // give back
                    continue;
                }
                let mut dp = [0.0; 3];
                let mut dh = 0.0;
                let mut dm = 0.0;
                for r in &o.records {
                    let d = du(r.cell);
                    for (p, u) in dp.iter_mut().zip(d.u) {
                        *p += r.slope * u;
                    }
                    dh += r.heat_slope * d.t;
                    if r.vapour_slope != 0.0 {
                        let x = r.vapour_slope * d.y;
                        dm += x;
                        dh += x * r.liquid_h;
                    }
                }
                if dm == 0.0 {
                    for (v, d) in o.p.v.iter_mut().zip(dp) {
                        *v += d / m;
                    }
                    o.p.temp += dh / (m * cp_p);
                } else {
                    // the mass changes: momentum and enthalpy are kept
                    // exactly, the velocity, temperature and diameter follow
                    let m1 = (m + dm).max(0.0);
                    o.vapour -= dm;
                    if m1 > 0.0 {
                        for (v, d) in o.p.v.iter_mut().zip(dp) {
                            *v = (m * *v + d) / m1;
                        }
                        o.p.temp = (m * cp_p * o.p.temp + dh) / (m1 * cp_p);
                        o.p.d = (m1 / (o.p.count * rho_p * std::f64::consts::PI / 6.0)).cbrt();
                    } else {
                        o.p.d = 0.0;
                        o.fate = Fate::Evaporated;
                    }
                }
            }
            for (j, &c) in sources.cells.iter().enumerate() {
                let d = du(c);
                for ((b, i), dk) in books.from_gas.iter_mut().zip(sources.impulse[j]).zip(d.u) {
                    *b -= i - sources.slope[j] * dk;
                }
                books.enthalpy[2] -=
                    sources.heat[j] - sources.heat_slope[j] * d.t - sources.vapour_slope_h[j] * d.y;
                books.mass[2] -= sources.vapour[j] - sources.vapour_slope[j] * d.y;
            }
        } else {
            for (j, imp) in sources.impulse.iter().enumerate() {
                for (b, i) in books.from_gas.iter_mut().zip(imp) {
                    *b -= i;
                }
                books.enthalpy[2] -= sources.heat[j];
                books.mass[2] -= sources.vapour[j];
            }
        }
        let mut kept = Parcels::default();
        let mut cross: Vec<Crossing> = Vec::new();
        for o in &outcomes {
            cross.extend_from_slice(&o.cross);
        }
        // D263: sorted by (datum, parcel, sub-step, order), summed in turn
        cross.sort_by_key(|c| c.key());
        for c in &cross {
            let di = c.datum as usize;
            let m = mass_of(c.d, c.count, rho_p);
            let spec = self.datums[di].hist;
            self.datum_tally[di].book(&spec, c, m);
        }
        // the books in the order of the parcel identities
        let mut order: Vec<usize> = (0..outcomes.len()).collect();
        order.sort_by_key(|&i| outcomes[i].p.id);
        for &i in &order {
            let o = &outcomes[i];
            let m = mass_of(o.p.d, o.p.count, rho_p);
            for k in 0..3 {
                books.gravity[k] += o.weight[k];
                books.walls[k] += o.walls[k];
                if o.fate == Fate::Active {
                    books.end[k] += m * o.p.v[k];
                } else {
                    books.removed[k] -= m * o.p.v[k];
                }
            }
            if o.fate == Fate::Active {
                books.enthalpy[1] += m * cp_p * o.p.temp;
                books.mass[1] += m;
            } else {
                books.enthalpy[4] -= m * cp_p * o.p.temp;
                books.mass[4] -= m;
            }
        }
        // the liquid that evaporated, in the order of the identities
        for &i in &order {
            self.tally.evaporated_mass += outcomes[i].vapour;
        }
        if !self.cfg.coupling.two_way() {
            // one way the gas books nothing: the vapour is what the parcels
            // lost
            books.mass[2] = books.mass[1] - books.mass[0] - books.mass[3] - books.mass[4];
        }
        self.tally.critical_parcels = 0;
        self.tally.critical_mass = 0.0;
        for o in outcomes {
            self.tally.rebounds += o.rebounds;
            self.tally.substeps += o.substeps;
            let m = mass_of(o.p.d, o.p.count, rho_p);
            match o.fate {
                Fate::Active => {
                    if self.surface_state(o.p.temp).is_some_and(|s| !s.two_phase) {
                        self.tally.critical_parcels += 1;
                        self.tally.critical_mass += m;
                    }
                    kept.push(&o.p)
                }
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
                Fate::Evaporated => {
                    self.tally.evaporated_parcels += 1;
                }
            }
        }
        self.parcels = kept;
        self.time = t1;
        books
    }

    /// The vapour that brings a gas of mass `mass` in the state `gs` to
    /// adiabatic saturation with this liquid (kg; 0 when the gas is
    /// saturated or does not carry the vapour): Δ = M (Y_s − Y*)/(1 − Y_s)
    /// at the temperature T where M c_p (T* − T) = Δ L(T) (D253).
    pub fn saturation_capacity(&self, gs: &GasState, mass: f64) -> f64 {
        let Some(v) = self.vapour else {
            return 0.0;
        };
        if !gs.has_vapour() || mass <= 0.0 {
            return 0.0;
        }
        let ys = |t: f64| {
            foehncraft_kernels::evaporation::surface_fraction(
                self.saturation_pressure(t),
                gs.p0,
                v.species.molar_mass,
                gs.m_gas,
            )
        };
        let delta = |t: f64| {
            let s = ys(t).min(0.999);
            (mass * (s - gs.y_v) / (1.0 - s)).max(0.0)
        };
        let g = |t: f64| mass * gs.cp * (gs.t - t) - delta(t) * self.latent_heat(t);
        let (mut lo, mut hi) = (200.0f64.min(gs.t), gs.t);
        if g(hi) >= 0.0 {
            return 0.0;
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if g(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        delta(0.5 * (lo + hi))
    }

    /// The saturation limiter of D253 on a pass: the evaporating records of
    /// each cell with a factor `factor(cell)` below 1 keep that share of
    /// their vapour; the rest goes back to its parcel as liquid at the
    /// temperature of the record, with the momentum it carried, before the
    /// gas sees the sources. Returns the
    /// cells limited and the mass held back (kg).
    pub fn limit_vapour(&self, pass: &mut Pass, factor: &dyn Fn(u32) -> f64) -> (u32, f64) {
        let (mut cells, mut back) = (0u32, 0.0);
        for &c in &pass.sources.cells {
            if factor(c) < 1.0 {
                cells += 1;
            }
        }
        if cells == 0 {
            return (0, 0.0);
        }
        let mut per: Vec<(u64, f64)> = Vec::new();
        for o in &mut pass.outcomes {
            let x = self.hold_back(o, factor, true);
            if x > 0.0 {
                per.push((o.p.id, x));
            }
        }
        per.sort_by_key(|p| p.0);
        for (_, x) in per {
            back += x;
        }
        let all: Vec<Record> = pass.records().copied().collect();
        pass.sources = CellSources::from_records(&all);
        (cells, back)
    }

    /// One parcel or stream of the limiter of D253: each evaporating record
    /// with a factor below 1 keeps that share of its vapour; the rest goes
    /// back to the parcel as liquid at the temperature of the record (the
    /// enthalpy c_l T_r, so the gas keeps the latent heat it would have
    /// given) with the momentum it carried; `revive` makes a parcel that had
    /// evaporated completely active again (transient). Returns the mass
    /// held back.
    fn hold_back(&self, o: &mut Outcome, factor: &dyn Fn(u32) -> f64, revive: bool) -> f64 {
        let rho_p = self.cfg.liquid.density;
        let cp_p = self.cfg.liquid.specific_heat;
        let (mut x, mut xh) = (0.0, 0.0);
        let v_end = o.p.v;
        for r in &mut o.records {
            let f = factor(r.cell);
            if f >= 1.0 || r.vapour <= 0.0 {
                continue;
            }
            let dx = (1.0 - f) * r.vapour;
            r.vapour -= dx;
            r.heat -= dx * r.liquid_h;
            r.vapour_slope *= f;
            for (i, v) in r.impulse.iter_mut().zip(v_end) {
                *i -= dx * v;
            }
            x += dx;
            xh += dx * r.liquid_h;
        }
        if x == 0.0 {
            return 0.0;
        }
        let m = mass_of(o.p.d, o.p.count, rho_p);
        let m1 = m + x;
        o.p.temp = (m * cp_p * o.p.temp + xh) / (m1 * cp_p);
        o.p.d = (m1 / (o.p.count * rho_p * std::f64::consts::PI / 6.0)).cbrt();
        o.vapour -= x;
        if revive && o.fate == Fate::Evaporated {
            o.fate = Fate::Active;
        }
        x
    }

    /// Liquid mass of the parcels in the flow (kg), summed in the order of
    /// their identities.
    pub fn mass_in_flow_by_id(&self) -> f64 {
        let mut m = 0.0;
        for i in self.by_id() {
            m += self.parcel_mass(i);
        }
        m
    }

    /// Sensible enthalpy Σ m c_p T of the parcels in the flow (J, from 0 K),
    /// summed in the order of their identities.
    pub fn enthalpy(&self) -> f64 {
        let cp = self.cfg.liquid.specific_heat;
        let mut h = 0.0;
        for i in self.by_id() {
            h += self.parcel_mass(i) * cp * self.parcels.temp[i];
        }
        h
    }

    /// Momentum of the parcels in the flow (kg m/s), summed in the order of
    /// their identities.
    pub fn momentum(&self) -> [f64; 3] {
        let mut p = [0.0; 3];
        for i in self.by_id() {
            let m = self.parcel_mass(i);
            for (k, pk) in p.iter_mut().enumerate() {
                *pk += m * self.parcels.v[i][k];
            }
        }
        p
    }

    /// Per cell with parcels: the cell, the number of parcels and the
    /// liquid volume (m³), ascending in the cell (health numbers of D254
    /// (5)); summed in the order of the parcel identities.
    pub fn cell_loading(&self) -> Vec<(u32, u32, f64)> {
        let mut v: Vec<(u32, u64, f64)> = (0..self.parcels.len())
            .map(|i| {
                let d = self.parcels.d[i];
                (
                    self.parcels.cell[i],
                    self.parcels.id[i],
                    self.parcels.count[i] * std::f64::consts::PI / 6.0 * d * d * d,
                )
            })
            .collect();
        v.sort_by_key(|x| (x.0, x.1));
        let mut out: Vec<(u32, u32, f64)> = Vec::new();
        for (c, _, vol) in v {
            match out.last_mut() {
                Some(l) if l.0 == c => {
                    l.1 += 1;
                    l.2 += vol;
                }
                _ => out.push((c, 1, vol)),
            }
        }
        out
    }

    /// The parcels due before `t1`: parcel k of an injector is born at
    /// start + (k + 1/2)/rate, so the mean mass flow is the one asked for.
    /// Nothing is booked here (the pass may be repeated); [`Spray::commit`]
    /// books the release.
    fn release(&self, env: &Env, t1: f64) -> Result<Release> {
        let rho_p = self.cfg.liquid.density;
        let mut out = Release {
            released: self.injectors.iter().map(|i| i.released).collect(),
            ..Default::default()
        };
        for (ii, inj) in self.cfg.injectors.iter().enumerate() {
            let rate = inj.parcels_per_second;
            let end = inj.end.unwrap_or(f64::INFINITY);
            let axis = unit(inj.direction);
            let (e1, e2) = frame(axis);
            let cos_a = inj.half_angle.to_radians().cos();
            let m_parcel = inj.mass_flow / rate;
            loop {
                let k = out.released[ii];
                let born = inj.start + (k as f64 + 0.5) / rate;
                if born > t1 || born > end {
                    break;
                }
                out.released[ii] += 1;
                let id = ((ii as u64) << 48) | k;
                let dir = direction(id, axis, e1, e2, cos_a);
                let d = match &inj.size {
                    Some(sz) => {
                        sz.quantile(rng::uniform_keyed([stream::SPRAY_INJECTION, id, 1, 0]))
                    }
                    None => inj.diameter,
                };
                let count = m_parcel / (rho_p * std::f64::consts::PI / 6.0 * d * d * d);
                let Some(cell) = env.loc.locate_search(inj.position) else {
                    bail!(
                        "[[spray.injector]] {}: the position lies outside the mesh",
                        inj.name
                    );
                };
                out.fresh.push(One {
                    id,
                    x: inj.position,
                    v: [dir[0] * inj.speed, dir[1] * inj.speed, dir[2] * inj.speed],
                    d,
                    count,
                    born,
                    temp: inj.temperature.unwrap_or(self.cfg.liquid.temperature),
                    eddy: Eddy::default(),
                    cell,
                    face: None,
                });
                out.masses.push(m_parcel);
            }
        }
        Ok(out)
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

    /// The crossings of the datum surfaces by the chord from `x0` (cell
    /// `c0`) to `x1` (cell `c1`) of sub-step `sub` of parcel `p`, which
    /// starts at time `t0` and lasts `dt`; `vel(τ)` is the velocity τ after
    /// the start.
    #[allow(clippy::too_many_arguments)]
    fn crossings(
        &self,
        out: &mut Vec<Crossing>,
        p: &One,
        sub: u64,
        (x0, c0): ([f64; 3], usize),
        (x1, c1): ([f64; 3], usize),
        t0: f64,
        dt: f64,
        vel: &dyn Fn(f64) -> [f64; 3],
    ) {
        let lerp = |f: f64| {
            [
                x0[0] + f * (x1[0] - x0[0]),
                x0[1] + f * (x1[1] - x0[1]),
                x0[2] + f * (x1[2] - x0[2]),
            ]
        };
        let mut push = |datum: usize, f: f64, dir: i8| {
            let k = out.len() as u32;
            out.push(Crossing {
                datum: datum as u32,
                id: p.id,
                sub,
                k,
                dir,
                t: t0 + f * dt,
                x: lerp(f),
                v: vel(f * dt),
                d: p.d,
                count: p.count,
            });
        };
        for (di, d) in self.datums.iter().enumerate() {
            let Shape::Surface(s) = &d.shape else {
                continue;
            };
            match s {
                Surface::Plane {
                    point,
                    normal,
                    bound,
                } => {
                    if let Some((f, dir)) = plane_crossing(*point, *normal, bound, x0, x1) {
                        push(di, f, dir);
                    }
                }
                Surface::Cells { side, plane, faces } => {
                    let (a, b) = (side[c0], side[c1]);
                    if a == 0 || b == 0 || a == b {
                        continue;
                    }
                    let f = match (plane, faces) {
                        (Some((pt, n, bound)), _) => {
                            let s0 = crate::datum::dot(crate::datum::sub(x0, *pt), *n);
                            let s1 = crate::datum::dot(crate::datum::sub(x1, *pt), *n);
                            let f = if s0 != s1 {
                                (s0 / (s0 - s1)).clamp(0.0, 1.0)
                            } else {
                                0.5
                            };
                            let x = lerp(f);
                            let h = crate::datum::dot(crate::datum::sub(x, *pt), *n);
                            let xp = [x[0] - h * n[0], x[1] - h * n[1], x[2] - h * n[2]];
                            if !bound.contains(*pt, xp) {
                                continue;
                            }
                            f
                        }
                        (None, Some(bvh)) => bvh.first_hit(x0, x1).map_or(0.5, |h| h.t),
                        (None, None) => 0.5,
                    };
                    push(di, f, b);
                }
                Surface::Boundary { .. } => {}
                Surface::Triangles { bvh } => {
                    for h in bvh.all_hits(x0, x1) {
                        push(di, h.t, if h.entering { -1 } else { 1 });
                    }
                }
            }
        }
    }

    /// Books the arrival of parcel `p` at boundary `b` (`dir` +1 into the
    /// boundary, −1 back from it) at time `t`.
    fn boundary_crossing(
        &self,
        out: &mut Vec<Crossing>,
        p: &One,
        sub: u64,
        b: usize,
        dir: i8,
        t: f64,
    ) {
        for (di, d) in self.datums.iter().enumerate() {
            if let Shape::Surface(Surface::Boundary { boundary }) = &d.shape
                && *boundary == b
            {
                let k = out.len() as u32;
                out.push(Crossing {
                    datum: di as u32,
                    id: p.id,
                    sub,
                    k,
                    dir,
                    t,
                    x: p.x,
                    v: p.v,
                    d: p.d,
                    count: p.count,
                });
            }
        }
    }

    /// Tracks one parcel for `time` seconds, to the end `t1` of the step.
    ///
    /// One way (SP1): sub-steps end at `max_courant` cell sizes and at
    /// boundary faces. Two way (SP2, D247 (a), D251): they end in addition
    /// at every interior face, so that each sub-step lies in one cell and
    /// sees that cell's frozen gas, and when the particle Reynolds number
    /// would change by more than `max_reynolds_change`; every sub-step
    /// writes its exchange record (D249).
    /// Puts the parcel into a new eddy (eddy-lifetime model, D243): a
    /// Gaussian fluctuation of variance 2k/3 per component, drawn from the
    /// stream `SPRAY_DISPERSION` keyed by the parcel and its eddy count, so
    /// every pass of a step draws the same numbers; the eddy lasts the
    /// interaction time of Gosman and Ioannides (1983).
    fn new_eddy(&self, p: &mut One, gs: &GasState, rho_p: f64) {
        p.eddy.n += 1;
        let sigma = (2.0 / 3.0 * gs.tke).sqrt();
        let mut u = [0.0; 3];
        let comps = if self.planar { 2 } else { 3 };
        for (k, uk) in u.iter_mut().enumerate().take(comps) {
            *uk = sigma * normal([stream::SPRAY_DISPERSION, p.id, p.eddy.n, 2 * k as u64]);
        }
        let ug = [gs.u[0] + u[0], gs.u[1] + u[1], gs.u[2] + u[2]];
        let s = parcel_step(&ParcelInput {
            v: p.v,
            u: ug,
            rho_g: gs.rho,
            mu_g: gs.mu,
            d: p.d,
            rho_p,
            g: self.gravity,
            dt: 0.0,
        });
        let slip = norm([ug[0] - p.v[0], ug[1] - p.v[1], ug[2] - p.v[2]]);
        p.eddy.u = u;
        p.eddy.left = eddy_times(gs.tke, gs.eps, s.tau, slip).interaction;
    }

    fn track(
        &self,
        env: &Env,
        gas: &dyn Gas,
        mut p: One,
        time: f64,
        t1: f64,
        opts: TrackOpts,
    ) -> Outcome {
        let rho_p = self.cfg.liquid.density;
        let two_way = self.cfg.coupling.two_way();
        let re_frac = self.cfg.coupling.max_reynolds_change;
        let mut mass;
        let mut vapour_total = 0.0;
        let cp_p = self.cfg.liquid.specific_heat;
        let g_eff = |rho_g: f64| 1.0 - rho_g / rho_p;
        let mut cross = Vec::new();
        let mut records: Vec<Record> = Vec::new();
        let mut weight = [0.0; 3];
        let mut walls = [0.0; 3];
        let any_datum = !self.datums.is_empty();
        let mut left = time;
        let mut rebounds = 0u64;
        let mut substeps = 0u64;
        let mut fate = Fate::Active;
        let tiny = 1e-14 * time.max(1e-300);
        let dispersion = self.cfg.dispersion == Dispersion::EddyLifetime;
        // the time of the eddy is spent with the sub-steps (D251: a
        // sub-step ends at the end of an eddy)
        let mut last_left = left;
        let mut path: Vec<PathPoint> = Vec::new();
        while left > tiny {
            if p.d <= 0.0 {
                // evaporated completely (SP5)
                fate = Fate::Evaporated;
                break;
            }
            mass = mass_of(p.d, p.count, rho_p);
            p.eddy.left -= last_left - left;
            last_left = left;
            if opts.paths {
                path.push(PathPoint {
                    t: time - left,
                    x: p.x,
                    v: p.v,
                    temp: p.temp,
                });
            }
            substeps += 1;
            if substeps as usize > self.cfg.max_substeps {
                fate = Fate::Lost;
                break;
            }
            let mut gs = gas.at(p.cell, p.x);
            let h = env.geo.cell_volume[p.cell].cbrt();
            if dispersion && gs.tke > 0.0 && gs.eps > 0.0 {
                // an eddy whose rest is below what the step time can
                // resolve is over (else its sub-steps would not advance)
                if p.eddy.left <= tiny {
                    self.new_eddy(&mut p, &gs, rho_p);
                }
                for k in 0..3 {
                    gs.u[k] += p.eddy.u[k];
                }
            } else {
                // no turbulence here: the parcel sees the mean gas
                p.eddy.u = [0.0; 3];
                p.eddy.left = 0.0;
            }
            let d_start = p.d;
            let input = |dt: f64| ParcelInput {
                v: p.v,
                u: gs.u,
                rho_g: gs.rho,
                mu_g: gs.mu,
                d: d_start,
                rho_p,
                g: self.gravity,
                dt,
            };
            let mut sub = left;
            if dispersion && p.eddy.left > 0.0 && p.eddy.left < sub {
                sub = p.eddy.left;
            }
            let mut s = parcel_step(&input(sub));
            let len = norm(s.dx);
            if len > self.cfg.max_courant * h {
                sub = sub * self.cfg.max_courant * h / len;
                s = parcel_step(&input(sub));
            }
            if two_way {
                // D251: the coefficients are frozen over the sub-step, so
                // it ends before the Reynolds number moves too far
                let re_end = |v: [f64; 3]| {
                    let w = [gs.u[0] - v[0], gs.u[1] - v[1], gs.u[2] - v[2]];
                    gs.rho * p.d * norm(w) / gs.mu
                };
                let lim = re_frac * s.re.max(1.0);
                let mut halvings = 0;
                while (re_end(s.v) - s.re).abs() > lim && halvings < 40 {
                    sub *= 0.5;
                    s = parcel_step(&input(sub));
                    halvings += 1;
                }
            }
            // D251 with evaporation (SP5): a sub-step ends before d² has
            // fallen by more than max_mass_change and before the droplet
            // temperature moves by more than max_temperature_change; a
            // droplet below min_diameter whose life ends inside the
            // sub-step finishes there (the d² law is exact with frozen
            // coefficients, so the end of the life is found exactly)
            let evaporating = self.evaporates() && gs.has_vapour();
            if evaporating && let Some(e) = self.evaporation(&gs, p.d, p.temp, s.re, sub) {
                let ev = &self.cfg.evaporation;
                let d02 = p.d * p.d;
                let mut finishing = false;
                let cut = if e.k > 0.0 && e.life <= sub && p.d < ev.min_diameter {
                    finishing = true;
                    e.life
                } else if e.k.abs() * sub > ev.max_mass_change * d02 {
                    ev.max_mass_change * d02 / e.k.abs()
                } else {
                    sub
                };
                if cut < sub {
                    sub = cut;
                    s = parcel_step(&input(sub));
                }
                let mut halvings = 0;
                while !finishing && halvings < 40 {
                    let Some(e) = self.evaporation(&gs, p.d, p.temp, s.re, sub) else {
                        break;
                    };
                    if (e.tp - p.temp).abs() <= ev.max_temperature_change {
                        break;
                    }
                    sub *= 0.5;
                    s = parcel_step(&input(sub));
                    halvings += 1;
                }
            }
            let target = [p.x[0] + s.dx[0], p.x[1] + s.dx[1], p.x[2] + s.dx[2]];
            // the face the previous sub-step ended on (two-way tracking)
            let on_face = p.face.map(|q| q.face);
            let start = p.face.take().or_else(|| {
                env.loc.position_in(p.cell, p.x).or_else(|| {
                    env.loc
                        .locate_search(p.x)
                        .and_then(|c| env.loc.position_in(c, p.x))
                })
            });
            let Some(start) = start else {
                fate = Fate::Lost;
                break;
            };
            let t_sub = t1 - left;
            let vel = |tau: f64| parcel_step(&input(tau)).v;
            // the exchange of a sub-step of length dt that ends with the
            // velocity v1 (D249: the parcel's own change less its weight);
            // returns the parcel temperature and diameter at its end and the
            // liquid it lost to vapour (the exact relaxations of SP2b and
            // the d² law of SP5)
            let book = |records: &mut Vec<Record>,
                        weight: &mut [f64; 3],
                        p: &One,
                        dt: f64,
                        v1: [f64; 3],
                        tau: f64|
             -> (f64, f64, f64) {
                let ge = g_eff(gs.rho);
                let slip = [gs.u[0] - p.v[0], gs.u[1] - p.v[1], gs.u[2] - p.v[2]];
                let re = gs.rho * p.d * norm(slip) / gs.mu;
                let ev = if evaporating {
                    self.evaporation(&gs, p.d, p.temp, re, dt)
                } else {
                    None
                };
                let mut imp = [0.0; 3];
                let m1 = match &ev {
                    None => {
                        for k in 0..3 {
                            let w = mass * self.gravity[k] * ge * dt;
                            weight[k] += w;
                            imp[k] = -(mass * (v1[k] - p.v[k]) - w);
                        }
                        mass
                    }
                    Some(e) => {
                        // the weight of the shrinking droplet, ∫ m(t) dt
                        let mi = mass_integral(mass, p.d, e.d, e.k, dt);
                        let m1 = mass_of(e.d, p.count, rho_p);
                        for k in 0..3 {
                            let w = self.gravity[k] * ge * mi;
                            weight[k] += w;
                            // the vapour leaves with the droplet's velocity
                            imp[k] = -(m1 * v1[k] - mass * p.v[k] - w);
                        }
                        m1
                    }
                };
                let (tp1, heat_decay) = match &ev {
                    Some(e) => (e.tp, e.decay),
                    None if gs.has_heat() => {
                        let h = heat_step(&HeatInput {
                            tp: p.temp,
                            tg: gs.t,
                            re,
                            pr: gs.mu * gs.cp / gs.k,
                            k_g: gs.k,
                            d: p.d,
                            rho_p,
                            cp_p,
                            dt,
                        });
                        (h.tp, h.decay)
                    }
                    None => (p.temp, 1.0),
                };
                let vapour = mass - m1;
                if two_way {
                    // the vapour slope: the chord to the saturation state
                    // at the droplet surface (D250), capped so that the
                    // closing step cannot take more than the parcel holds
                    let (vapour_slope, heat) = match &ev {
                        Some(e) => {
                            let gap = e.ys - gs.y_v;
                            let mut b = if e.d > 0.0 && gap.abs() > 1e-12 {
                                vapour / gap
                            } else {
                                0.0
                            };
                            if !(b > 0.0 && b.is_finite()) {
                                b = 0.0;
                            }
                            (b, -(m1 * cp_p * tp1 - mass * cp_p * p.temp))
                        }
                        None => (0.0, -mass * cp_p * (tp1 - p.temp)),
                    };
                    records.push(Record {
                        cell: p.cell as u32,
                        id: p.id,
                        sub: substeps as u32,
                        impulse: imp,
                        slope: 0.0,
                        dt,
                        decay: (-(dt / tau)).exp_of(),
                        heat,
                        heat_slope: 0.0,
                        heat_decay,
                        vapour,
                        vapour_slope,
                        vapour_y: if evaporating { gs.y_v } else { 0.0 },
                        liquid_h: if opts.steady || !evaporating {
                            cp_p * tp1
                        } else {
                            // the droplet's heat follows a change of its
                            // evaporation only as far as it relaxes in the
                            // step: the latent heat with the weight e_h
                            let hl = cp_p * tp1;
                            hl + heat_decay * (self.vapour_enthalpy(tp1) - hl)
                        },
                    });
                }
                (tp1, ev.map_or(p.d, |e| e.d), vapour)
            };
            // A parcel that moves along a face, with the gas of the two
            // cells pushing it towards the face from both sides, would end
            // every sub-step at that face after no time and chatter between
            // the cells: a sub-step that starts on a face and goes straight
            // back through it takes its full length instead (booked to the
            // cell it starts in, the path below).
            let ff = if two_way {
                Some(first_face(env, start, target))
            } else {
                None
            };
            let ff = ff.filter(|f| {
                !matches!(f, FirstFace::Interior { lambda, next, .. }
                    if on_face == Some(next.face) && *lambda < 1e-9)
            });
            if let Some(ff) = ff {
                match ff {
                    FirstFace::Reached(cell) => {
                        let bk = book(&mut records, &mut weight, &p, sub, s.v, s.tau);
                        p.temp = bk.0;
                        p.d = bk.1;
                        vapour_total += bk.2;
                        if any_datum {
                            self.crossings(
                                &mut cross,
                                &p,
                                substeps,
                                (p.x, p.cell),
                                (target, cell),
                                t_sub,
                                sub,
                                &vel,
                            );
                        }
                        p.x = target;
                        p.v = s.v;
                        p.cell = cell;
                        left -= sub;
                        continue;
                    }
                    FirstFace::Interior { lambda, x, next } => {
                        let tf =
                            time_at_fraction(&|t| parcel_step(&input(t)).dx, s.dx, lambda, sub);
                        let sf = parcel_step(&input(tf));
                        let bk = book(&mut records, &mut weight, &p, tf, sf.v, sf.tau);
                        p.temp = bk.0;
                        p.d = bk.1;
                        vapour_total += bk.2;
                        if any_datum {
                            self.crossings(
                                &mut cross,
                                &p,
                                substeps,
                                (p.x, p.cell),
                                (x, next.cell),
                                t_sub,
                                tf,
                                &vel,
                            );
                        }
                        p.x = x;
                        p.v = sf.v;
                        p.cell = next.cell;
                        p.face = Some(next);
                        left -= tf;
                        continue;
                    }
                    FirstFace::Boundary { face, lambda, x } => {
                        let tf =
                            time_at_fraction(&|t| parcel_step(&input(t)).dx, s.dx, lambda, sub);
                        let sf = parcel_step(&input(tf));
                        let bk = book(&mut records, &mut weight, &p, tf, sf.v, sf.tau);
                        p.temp = bk.0;
                        p.d = bk.1;
                        vapour_total += bk.2;
                        if any_datum {
                            self.crossings(
                                &mut cross,
                                &p,
                                substeps,
                                (p.x, p.cell),
                                (x, p.cell),
                                t_sub,
                                tf,
                                &vel,
                            );
                        }
                        p.x = x;
                        p.v = sf.v;
                        left -= tf;
                        let m_now = mass_of(p.d, p.count, rho_p);
                        match self.hit_boundary(
                            env,
                            &mut p,
                            face,
                            h,
                            &mut cross,
                            substeps,
                            t1 - left,
                            &mut walls,
                            m_now,
                        ) {
                            Some(f) => {
                                fate = f;
                                break;
                            }
                            None => {
                                rebounds += 1;
                                continue;
                            }
                        }
                    }
                    FirstFace::Lost => {
                        fate = Fate::Lost;
                        break;
                    }
                }
            }
            let tr = env.loc.track(start, target);
            match tr.end {
                TrackEnd::Reached => {
                    let bk = book(&mut records, &mut weight, &p, sub, s.v, s.tau);
                    p.temp = bk.0;
                    p.d = bk.1;
                    vapour_total += bk.2;
                    if any_datum {
                        self.crossings(
                            &mut cross,
                            &p,
                            substeps,
                            (p.x, p.cell),
                            (target, tr.position.cell),
                            t_sub,
                            sub,
                            &vel,
                        );
                    }
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
                    let bk = book(&mut records, &mut weight, &p, thit, sh.v, sh.tau);
                    p.temp = bk.0;
                    p.d = bk.1;
                    vapour_total += bk.2;
                    if any_datum {
                        self.crossings(
                            &mut cross,
                            &p,
                            substeps,
                            (p.x, p.cell),
                            (tr.position.x, tr.position.cell),
                            t_sub,
                            thit,
                            &vel,
                        );
                    }
                    p.x = tr.position.x;
                    p.v = sh.v;
                    p.cell = tr.position.cell;
                    left -= thit;
                    let m_now = mass_of(p.d, p.count, rho_p);
                    match self.hit_boundary(
                        env,
                        &mut p,
                        face,
                        h,
                        &mut cross,
                        substeps,
                        t1 - left,
                        &mut walls,
                        m_now,
                    ) {
                        Some(f) => {
                            fate = f;
                            break;
                        }
                        None => rebounds += 1,
                    }
                }
                TrackEnd::Processor { .. } | TrackEnd::Lost => {
                    // the walk failed: one search for the target
                    match env.loc.locate_search(target) {
                        Some(c) => {
                            let bk = book(&mut records, &mut weight, &p, sub, s.v, s.tau);
                            p.temp = bk.0;
                            p.d = bk.1;
                            vapour_total += bk.2;
                            if any_datum {
                                self.crossings(
                                    &mut cross,
                                    &p,
                                    substeps,
                                    (p.x, p.cell),
                                    (target, c),
                                    t_sub,
                                    sub,
                                    &vel,
                                );
                            }
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
        p.eddy.left -= last_left - left;
        p.face = None;
        if fate == Fate::Active && p.d <= 0.0 {
            fate = Fate::Evaporated;
        }
        // the slopes with the mass at the end (SP5: the parcel may have
        // lost liquid; one that evaporated completely has none, D250)
        let m_end = mass_of(p.d, p.count, rho_p);
        // a stream has no closing step: its vapour slopes stay the chords
        // to saturation, the stabilising part of the held sources (D260)
        if !opts.steady {
            crate::coupling::cap_vapour_slopes(&mut records, m_end);
        }
        if opts.steady {
            crate::coupling::set_slopes_local(&mut records, m_end, m_end * cp_p);
        } else {
            crate::coupling::set_slopes(&mut records, m_end, m_end * cp_p);
        }
        let used = time - left.max(0.0);
        if opts.paths {
            path.push(PathPoint {
                t: used,
                x: p.x,
                v: p.v,
                temp: p.temp,
            });
        }
        Outcome {
            p,
            fate,
            rebounds,
            substeps,
            cross,
            records,
            weight,
            walls,
            used,
            vapour: vapour_total,
            path,
        }
    }

    /// A parcel `p` at boundary face `face`: books the arrival at the
    /// datums, then sticks, escapes, or bounces back (adding the impulse of
    /// the wall to `walls`). Returns the fate when the parcel leaves the
    /// flow, `None` after a bounce.
    #[allow(clippy::too_many_arguments)]
    fn hit_boundary(
        &self,
        env: &Env,
        p: &mut One,
        face: usize,
        h: f64,
        cross: &mut Vec<Crossing>,
        substeps: u64,
        t: f64,
        walls: &mut [f64; 3],
        mass: f64,
    ) -> Option<Fate> {
        let any_datum = !self.datums.is_empty();
        let rule = self.rule_of(env, face);
        let bnd = env.mesh.boundary_of(face);
        if any_datum && let Some(b) = bnd {
            self.boundary_crossing(cross, p, substeps, b, 1, t);
        }
        match rule.mode {
            BoundaryMode::Stick => Some(Fate::Stuck(face)),
            BoundaryMode::Escape => Some(Fate::Escaped(face)),
            BoundaryMode::Rebound | BoundaryMode::Mirror => {
                let a = env.geo.face_area[face];
                let n = unit(a);
                let vn = p.v[0] * n[0] + p.v[1] * n[1] + p.v[2] * n[2];
                let (en, et) = if rule.mode == BoundaryMode::Mirror {
                    (1.0, 1.0)
                } else {
                    (rule.e_normal, rule.e_tangential)
                };
                for ((vk, nk), wk) in p.v.iter_mut().zip(n).zip(walls.iter_mut()) {
                    let vt = *vk - vn * nk;
                    let v1 = et * vt - en * vn * nk;
                    *wk += mass * (v1 - *vk);
                    *vk = v1;
                }
                // a hair inside the cell, so the next track starts in it
                let eps = 1e-6 * h;
                for (xk, nk) in p.x.iter_mut().zip(n) {
                    *xk -= eps * nk;
                }
                p.face = None;
                if any_datum && let Some(b) = bnd {
                    self.boundary_crossing(cross, p, substeps, b, -1, t);
                }
                None
            }
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
            datums: self.datum_tally.clone(),
        })?)
    }

    /// Restores the state of [`state`](Self::state); the cells of the
    /// parcels are found again from their positions (the mesh may be
    /// numbered differently).
    pub fn restore(&mut self, env: &Env, bytes: &[u8]) -> Result<()> {
        let s: Saved = serde_json::from_slice(bytes)?;
        if s.format != FORMAT || s.version != 1 {
            bail!(
                "the restart holds a spray state of format {} version {}",
                s.format,
                s.version
            );
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
        if parcels.temp.len() != parcels.len() {
            // a restart written before SP2b: the injection temperature
            parcels.temp = parcels
                .id
                .iter()
                .map(|&id| {
                    let inj = &self.cfg.injectors[(id >> 48) as usize];
                    inj.temperature.unwrap_or(self.cfg.liquid.temperature)
                })
                .collect();
        }
        if parcels.eddy.len() != parcels.len() {
            // a restart written before SP2c: every parcel starts a new eddy
            parcels.eddy = vec![Eddy::default(); parcels.len()];
        }
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
        // the sums of the datums of the same name continue; a new datum
        // starts at the restart time
        for (d, t) in self.datums.iter().zip(self.datum_tally.iter_mut()) {
            *t = s
                .datums
                .iter()
                .find(|x| x.name == d.name && x.total.all.hist_n.len() == d.hist.bins)
                .cloned()
                .unwrap_or_else(|| DatumTally::new(&d.name, &d.hist, s.time));
        }
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
        let injector: Vec<f64> = self.parcels.id.iter().map(|&i| (i >> 48) as f64).collect();
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
                name: "temperature",
                values: &self.parcels.temp,
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
            PointField {
                name: "injector",
                values: &injector,
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
            ("liquid_density", Attr::F64(self.cfg.liquid.density)),
            ("units", Attr::Str("m, m/s, K, s, droplets")),
            ("parcels", Attr::U64(n as u64)),
        ];
        write_points(path, &self.parcels.x, &fields, &attrs)
    }
}

/// Per injector of `cfg`: parcels and liquid in the flow, the liquid
/// penetration and the farthest parcel (D265), from the identities,
/// positions and masses of parcels (in the flow, or of a parcel file);
/// `released` the counters of the injectors when known.
pub fn injector_stats(
    cfg: &SpraySection,
    id: &[u64],
    x: &[[f64; 3]],
    mass: &[f64],
    released: Option<&[u64]>,
) -> Vec<InjectorStats> {
    let frac = cfg.penetration_fraction;
    cfg.injectors
        .iter()
        .enumerate()
        .map(|(ii, inj)| {
            let axis = unit(inj.direction);
            let mut pts: Vec<(f64, u64, f64)> = (0..id.len())
                .filter(|&i| (id[i] >> 48) as usize == ii)
                .map(|i| {
                    let r = [
                        x[i][0] - inj.position[0],
                        x[i][1] - inj.position[1],
                        x[i][2] - inj.position[2],
                    ];
                    let dist = match cfg.penetration {
                        Penetration::Axial => r[0] * axis[0] + r[1] * axis[1] + r[2] * axis[2],
                        Penetration::Radial => norm(r),
                    };
                    (dist, id[i], mass[i])
                })
                .collect();
            pts.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            let total: f64 = pts.iter().map(|p| p.2).sum();
            let mut cum = 0.0;
            let mut pen = 0.0;
            for p in &pts {
                cum += p.2;
                pen = p.0;
                if cum >= frac * total {
                    break;
                }
            }
            let liquid_length = match &cfg.ecn {
                Some(e) => {
                    let own: Vec<usize> = (0..id.len())
                        .filter(|&i| (id[i] >> 48) as usize == ii)
                        .collect();
                    plv_lengths(e, inj, cfg.liquid.density, &own, id, x, mass)
                }
                None => Vec::new(),
            };
            InjectorStats {
                name: inj.name.clone(),
                released: released.and_then(|r| r.get(ii).copied()).unwrap_or(0),
                parcels: pts.len() as u64,
                mass: total,
                penetration: if pts.is_empty() { 0.0 } else { pen },
                farthest: pts.last().map_or(0.0, |p| p.0),
                liquid_length,
                vapour_penetration: None,
            }
        })
        .collect()
}

/// The unit line of sight of `[spray.ecn]` for the unit axis `a`: the
/// given direction, or the coordinate axis most nearly perpendicular to
/// `a` (the first of equals), with its part along `a` removed.
pub fn line_of_sight(e: &crate::config::Ecn, a: [f64; 3]) -> [f64; 3] {
    let l = e.line_of_sight.unwrap_or_else(|| {
        let mut k = 0;
        for j in 1..3 {
            if a[j].abs() < a[k].abs() {
                k = j;
            }
        }
        let mut v = [0.0; 3];
        v[k] = 1.0;
        v
    });
    let d = l[0] * a[0] + l[1] * a[1] + l[2] * a[2];
    unit([l[0] - d * a[0], l[1] - d * a[1], l[2] - d * a[2]])
}

/// The projected liquid volume of the parcels `own` of injector `inj`
/// (D315): the liquid volume m/ρ_l of each parcel in square bins of size
/// `e.bin` in the plane of the axis and the direction across both the axis
/// and the line of sight (the volume along the line of sight sums into
/// one bin), over the bin area. Returns per axial bin (from the nozzle,
/// index 0 = [0, bin)) the largest PLV across (m), summed in the order of
/// the parcel identities.
pub fn plv_profile(
    e: &crate::config::Ecn,
    inj: &Injector,
    rho_l: f64,
    own: &[usize],
    id: &[u64],
    x: &[[f64; 3]],
    mass: &[f64],
) -> Vec<f64> {
    let a = unit(inj.direction);
    let l = line_of_sight(e, a);
    // across the axis and the line of sight
    let b = unit([
        a[1] * l[2] - a[2] * l[1],
        a[2] * l[0] - a[0] * l[2],
        a[0] * l[1] - a[1] * l[0],
    ]);
    let mut ix = own.to_vec();
    ix.sort_by_key(|&i| id[i]);
    let mut bins: BTreeMap<(i64, i64), f64> = BTreeMap::new();
    for i in ix {
        let r = [
            x[i][0] - inj.position[0],
            x[i][1] - inj.position[1],
            x[i][2] - inj.position[2],
        ];
        let s = r[0] * a[0] + r[1] * a[1] + r[2] * a[2];
        if s < 0.0 {
            continue;
        }
        let t = r[0] * b[0] + r[1] * b[1] + r[2] * b[2];
        let key = ((s / e.bin).floor() as i64, (t / e.bin).floor() as i64);
        *bins.entry(key).or_insert(0.0) += mass[i] / rho_l;
    }
    let n = bins.keys().map(|k| k.0 + 1).max().unwrap_or(0).max(0) as usize;
    let mut prof = vec![0.0f64; n];
    let area = e.bin * e.bin;
    for (&(i, _), &v) in &bins {
        let p = &mut prof[i as usize];
        *p = p.max(v / area);
    }
    prof
}

/// The liquid length per PLV threshold of `e` (m, D315): the farthest
/// axial position where the PLV profile of [`plv_profile`] falls to the
/// threshold, linear between the bin centres of the last bin at or above
/// it and the next one (beyond the last bin the profile is zero); 0 when
/// no bin reaches it.
pub fn plv_lengths(
    e: &crate::config::Ecn,
    inj: &Injector,
    rho_l: f64,
    own: &[usize],
    id: &[u64],
    x: &[[f64; 3]],
    mass: &[f64],
) -> Vec<f64> {
    let prof = plv_profile(e, inj, rho_l, own, id, x, mass);
    e.plv
        .iter()
        .map(|&th| match prof.iter().rposition(|&p| p >= th) {
            None => 0.0,
            Some(k) => {
                let centre = (k as f64 + 0.5) * e.bin;
                let next = prof.get(k + 1).copied().unwrap_or(0.0);
                let w = (prof[k] - th) / (prof[k] - next);
                centre + w.clamp(0.0, 1.0) * e.bin
            }
        })
        .collect()
}

/// The moment sums of parcels given as arrays (a parcel file), summed in
/// the order of the identities, as the run sums them.
pub fn moments_of(
    spec: &crate::datum::HistSpec,
    id: &[u64],
    d: &[f64],
    count: &[f64],
    v: &[[f64; 3]],
    rho_p: f64,
) -> Moments {
    let mut ix: Vec<usize> = (0..id.len()).collect();
    ix.sort_by_key(|&i| id[i]);
    let mut m = Moments::new(spec);
    for i in ix {
        m.add(spec, d[i], count[i], mass_of(d[i], count[i], rho_p), v[i]);
    }
    m
}

/// The numbers of one injector (D265).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct InjectorStats {
    /// Name.
    pub name: String,
    /// Parcels released so far.
    pub released: u64,
    /// Parcels in the flow.
    pub parcels: u64,
    /// Liquid in the flow (kg).
    pub mass: f64,
    /// Liquid penetration (m): the distance within which the share
    /// `penetration_fraction` of the injector's liquid in the flow lies.
    pub penetration: f64,
    /// Distance of the farthest parcel (m); outlier-sensitive.
    pub farthest: f64,
    /// Liquid length per PLV threshold of `[spray.ecn]` (m, D315); empty
    /// without the table.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub liquid_length: Vec<f64>,
    /// Vapour penetration of `[spray.ecn]` (m, D315): set by the flow from
    /// the fuel mass fraction of the gas.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vapour_penetration: Option<f64>,
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

/// Where a straight segment from a position first leaves its cell.
enum FirstFace {
    /// The target lies in the cell (or on its far side, in that cell).
    Reached(usize),
    /// Through an interior face at the segment fraction `lambda`, at `x`;
    /// `next` is the walk's tetrahedron in the cell behind it.
    Interior {
        lambda: f64,
        x: [f64; 3],
        next: Position,
    },
    /// Through boundary face `face` at `lambda`, at `x`.
    Boundary {
        face: usize,
        lambda: f64,
        x: [f64; 3],
    },
    /// The walk failed.
    Lost,
}

/// The walk of [`Locator::walk`] stopped at the first mesh face the segment
/// from `start` to `target` crosses (two-way tracking, D247 (a)).
fn first_face(env: &Env, start: Position, target: [f64; 3]) -> FirstFace {
    let mut w = env.loc.start_walk(start, target);
    let end = env
        .loc
        .walk(&mut w, step_limit(env.mesh.n_cells), |_, _, l1| l1 >= 1.0);
    match end {
        WalkEnd::Reached => FirstFace::Reached(w.cell as usize),
        WalkEnd::Boundary { face } => FirstFace::Boundary {
            face: face as usize,
            lambda: 1.0,
            x: w.x,
        },
        WalkEnd::Stopped => {
            let f = w.face as usize;
            if env.mesh.boundary_of(f).is_some() || f >= env.mesh.neighbour.len() {
                return FirstFace::Boundary {
                    face: f,
                    lambda: w.lambda,
                    x: w.x,
                };
            }
            let (a, b) = (env.mesh.owner[f], env.mesh.neighbour[f]);
            let other = if a == w.cell { b } else { a } as usize;
            if other >= env.mesh.n_cells {
                return FirstFace::Lost;
            }
            FirstFace::Interior {
                lambda: w.lambda,
                x: w.x,
                next: Position {
                    x: w.x,
                    cell: other,
                    face: f,
                    edge: w.edge as usize,
                },
            }
        }
        WalkEnd::Processor { .. } | WalkEnd::Lost => FirstFace::Lost,
    }
}

/// The time in [0, `sub`] at which the displacement `dx(t)` projected on the
/// chord `chord` (= dx(sub)) reaches the fraction `lambda` of its length:
/// bisection, 50 halvings (pure, the same bits on every thread).
fn time_at_fraction(dx: &dyn Fn(f64) -> [f64; 3], chord: [f64; 3], lambda: f64, sub: f64) -> f64 {
    let l2 = chord[0] * chord[0] + chord[1] * chord[1] + chord[2] * chord[2];
    if l2.is_nan() || l2 <= 0.0 || lambda <= 0.0 {
        return 0.0;
    }
    if lambda >= 1.0 {
        return sub;
    }
    let goal = lambda * l2;
    let (mut lo, mut hi) = (0.0, sub);
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        let d = dx(mid);
        let f = d[0] * chord[0] + d[1] * chord[1] + d[2] * chord[2];
        if f < goal {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}
