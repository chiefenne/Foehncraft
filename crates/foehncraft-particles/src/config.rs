//! The `[spray]` table of the case file (L12, D240 to D246).

use anyhow::{Result, bail};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use foehncraft_kernels::breakup::{Cavitation, NozzleTurbulence};

/// What a parcel does when it reaches a boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum BoundaryMode {
    /// The parcel stays on the wall: it leaves the flow and its mass is
    /// booked as deposited on that boundary (the wall film of D123 takes
    /// over from here later).
    #[default]
    Stick,
    /// The parcel bounces back with the restitution coefficients of the
    /// boundary.
    Rebound,
    /// The parcel leaves the domain (inlets, outlets, open boundaries).
    Escape,
    /// Mirror reflection (symmetry planes and slip walls): the normal
    /// velocity changes sign, nothing is lost; the normal fluctuation of
    /// the parcel's eddy changes sign with it (the parcel goes on in the
    /// mirror image of the flow, D308).
    Mirror,
}

/// A boundary's behaviour for parcels, `[spray.boundary.<name>]`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct SprayBoundary {
    /// What a parcel does there; by default `wall_mode` on walls, `mirror`
    /// on symmetry planes and slip walls, `escape` on every other boundary.
    pub mode: Option<BoundaryMode>,
    /// Restitution of the normal velocity in `rebound` (0 to 1, default 1).
    pub restitution_normal: Option<f64>,
    /// Restitution of the tangential velocity in `rebound` (0 to 1,
    /// default 1: no friction).
    pub restitution_tangential: Option<f64>,
}

/// The liquid of the parcels: a constant density and specific heat, and
/// for the evaporation (SP4, SP5) the species it evaporates into, its
/// vapour pressure law and its latent heat at one temperature.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Liquid {
    /// Name, for the log and the output.
    pub name: String,
    /// Density (kg/m³).
    pub density: f64,
    /// Specific heat (J/(kg K)), constant (SP2b: sensible heat).
    pub specific_heat: f64,
    /// Temperature of the injected liquid (K), unless an injector gives
    /// its own.
    pub temperature: f64,
    /// The species of `[species]` the liquid evaporates into (SP4), for
    /// example `"H2O"`; needed with `[spray.evaporation]`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vapour: Option<String>,
    /// The vapour pressure law: `water` (the auxiliary equation of the
    /// IAPWS-95 formulation, Wagner and Pruss 1993) or `library` (the
    /// Wagner-form fit of the `library` liquid, which sets it).
    pub saturation: Saturation,
    /// Latent heat of vaporisation (J/kg) at `latent_heat_temperature`
    /// (default 2.4417e6 J/kg, water at 298.15 K, IAPWS-95 tables). At
    /// other temperatures L(T) = h_v(T) − h_l(T) with the vapour enthalpy of
    /// the species data (NASA-7) and h_l = c_l T (Kirchhoff's law), so the
    /// energy of the gas and of the parcels close with one property set.
    pub latent_heat: f64,
    /// The temperature of `latent_heat` (K, default 298.15).
    pub latent_heat_temperature: f64,
    /// Surface tension against the gas (N/m), constant; needed with
    /// `[spray.breakup]` (D300).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub surface_tension: Option<f64>,
    /// Dynamic viscosity of the liquid (Pa s), constant; needed with
    /// `[spray.breakup]` (D300).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub viscosity: Option<f64>,
    /// Vapour pressure at the injection temperature (Pa), constant: the
    /// collapse time of the cavities of `kh-act` (D307); without it the
    /// cavities only burst.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vapour_pressure: Option<f64>,
    /// A liquid of the property library (D320): `"water"`, `"n-heptane"`,
    /// `"n-dodecane"` (fits to the NIST WebBook saturation tables) or
    /// `"iso-octane"` (fits to published correlations, see
    /// docs/verification/DATA_SOURCES.md). It
    /// sets the vapour pressure law, the density and the specific heat (at
    /// `property_temperature`) and the latent heat (at
    /// `latent_heat_temperature`); leave `density`, `specific_heat`,
    /// `latent_heat` and `saturation` out with it. It also gives
    /// `surface_tension`, `viscosity` (at `property_temperature`) and
    /// `vapour_pressure` (at `temperature`) of the breakup unless they are
    /// set.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub library: Option<String>,
    /// The temperature (K) of the constant density and specific heat of a
    /// library liquid (default: `temperature`, the injected liquid's).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub property_temperature: Option<f64>,
}

impl Liquid {
    /// The library data of the liquid, if it names one.
    pub fn data(&self) -> Option<&'static foehncraft_kernels::liquids::LiquidData> {
        self.library
            .as_deref()
            .and_then(foehncraft_kernels::liquids::by_name)
    }

    /// Takes the properties of a library liquid (D320): the vapour pressure
    /// law, ρ and c_l at the property temperature and L at
    /// `latent_heat_temperature`. Refuses a library liquid with its own
    /// density, specific heat, latent heat or vapour pressure law.
    pub fn apply_library(&mut self) -> Result<()> {
        let Some(name) = self.library.clone() else {
            if self.property_temperature.is_some() {
                bail!("[spray.liquid] property_temperature needs a library liquid (`library`)");
            }
            return Ok(());
        };
        let Some(d) = foehncraft_kernels::liquids::by_name(&name) else {
            let names: Vec<&str> = foehncraft_kernels::liquids::ALL
                .iter()
                .map(|l| l.name)
                .collect();
            bail!(
                "[spray.liquid] library = \"{name}\" is not a liquid of the library; the liquids are {}",
                names.join(", ")
            );
        };
        let def = Liquid::default();
        let tp = self.property_temperature.unwrap_or(self.temperature);
        if !(tp.is_finite() && tp > 0.0) {
            bail!("[spray.liquid] property_temperature must be positive and finite (K)");
        }
        // a value given in the file must be the default or the library's
        // own (a second application, a re-read of the case)
        for (k, v, dv, lib) in [
            ("density", self.density, def.density, d.density(tp)),
            (
                "specific_heat",
                self.specific_heat,
                def.specific_heat,
                d.specific_heat(tp),
            ),
            (
                "latent_heat",
                self.latent_heat,
                def.latent_heat,
                d.latent_heat(self.latent_heat_temperature),
            ),
        ] {
            if v != dv && v != lib {
                bail!(
                    "[spray.liquid] {k} is given with library = \"{name}\": a library liquid takes it from the library"
                );
            }
        }
        self.density = d.density(tp);
        self.specific_heat = d.specific_heat(tp);
        self.latent_heat = d.latent_heat(self.latent_heat_temperature);
        // the constants of the breakup (D300, D307) unless the file gives
        // them: σ and μ at the property temperature, the vapour pressure at
        // the injected liquid's temperature
        self.surface_tension.get_or_insert(d.surface_tension(tp));
        self.viscosity.get_or_insert(d.viscosity(tp));
        self.vapour_pressure
            .get_or_insert(d.saturation_pressure(self.temperature));
        // water keeps the IAPWS-95 auxiliary equation as its vapour pressure
        self.saturation = if name == "water" {
            Saturation::Water
        } else {
            Saturation::Library
        };
        if self.name == def.name && name != "water" {
            self.name = name;
        }
        Ok(())
    }
}

/// The vapour pressure law of a liquid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Saturation {
    /// Water: Wagner and Pruss (1993), eq. 2.5 of Wagner and Pruss (2002).
    #[default]
    Water,
    /// The Wagner-form fit of the library liquid (`library`, D320); set by
    /// `library`, not written by hand.
    Library,
}

impl Default for Liquid {
    fn default() -> Self {
        Self {
            name: "water".into(),
            density: 998.0,
            specific_heat: 4181.0,
            temperature: 293.15,
            vapour: None,
            saturation: Saturation::Water,
            latent_heat: 2.4417e6,
            latent_heat_temperature: 298.15,
            surface_tension: None,
            viscosity: None,
            vapour_pressure: None,
            library: None,
            property_temperature: None,
        }
    }
}

/// The evaporation model of the droplets (SP5).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum EvaporationModel {
    /// No evaporation (default).
    #[default]
    None,
    /// Classical film theory with the Ranz-Marshall correlations (Spalding
    /// 1953; Ranz and Marshall 1952).
    RanzMarshall,
    /// Abramzon and Sirignano (1989): film thickening by the outflow of
    /// the vapour.
    AbramzonSirignano,
}

/// The phase equilibrium at the droplet surface (D311).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Equilibrium {
    /// Raoult's law with the saturation pressure of the liquid (default):
    /// right at low pressure.
    #[default]
    Raoult,
    /// The Peng-Robinson equation of state: the fugacities of the fuel and
    /// of the ambient gases equal in the surface liquid and the vapour,
    /// the latent heat from the same equation (Spray A, D311, D312, D314).
    PengRobinson,
}

/// `[spray.evaporation.peng_robinson]`: the data of the high-pressure
/// surface equilibrium (D311).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct PengRobinson {
    /// The fuel of the built-in high-pressure data (`"n-dodecane"`);
    /// default: the name of `[spray.liquid]`. Not needed when
    /// `critical_temperature`, `critical_pressure`, `acentric_factor` and
    /// `molar_mass` are given.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fuel: Option<String>,
    /// Critical temperature of the liquid (K), with the three keys below
    /// instead of `fuel` (then with the κ(ω) of Peng and Robinson).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub critical_temperature: Option<f64>,
    /// Critical pressure of the liquid (Pa).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub critical_pressure: Option<f64>,
    /// Acentric factor of the liquid.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub acentric_factor: Option<f64>,
    /// Molar mass of the liquid (kg/mol).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub molar_mass: Option<f64>,
    /// Mole fractions of the ambient gases far from the droplets, by species
    /// name (N2, O2, CO2, H2O, AR), normalised; their ratios hold in the
    /// vapour at the surface (default `{ N2 = 1.0 }`).
    pub ambient: BTreeMap<String, f64>,
    /// Binary interaction parameters k_ij of the fuel with the gases, by
    /// species name; a gas not listed takes the built-in value of the fuel
    /// (n-dodecane + N2: 0.205, fitted to measured bubble points) or 0.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub k_ij: BTreeMap<String, f64>,
    /// Spacing of the surface table over the droplet temperature (K,
    /// default 1).
    pub table_step: f64,
}

impl Default for PengRobinson {
    fn default() -> Self {
        Self {
            fuel: None,
            critical_temperature: None,
            critical_pressure: None,
            acentric_factor: None,
            molar_mass: None,
            ambient: BTreeMap::from([("N2".to_string(), 1.0)]),
            k_ij: BTreeMap::new(),
            table_step: 1.0,
        }
    }
}

/// The reference state of the gas film around a droplet at which the
/// evaporation takes its gas properties (D321).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Film {
    /// The 1/3 rule (Hubbard, Denny and Mills 1975; Yuen and Chen 1976):
    /// T_m = T_p + (T_g − T_p)/3, Y_m = Y_s + (Y_g − Y_s)/3, with the
    /// mixture properties of the vapour and the carrier gas from the
    /// species data (kinetic theory) at that state (default).
    #[default]
    OneThird,
    /// The properties of the gas at the droplet as the flow has them (the
    /// SP5 behaviour; the vapour c_p at the mean temperature).
    FarField,
}

/// `[spray.evaporation]` (SP5, D251).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Evaporation {
    /// `none` (default), `ranz-marshall` or `abramzon-sirignano`.
    pub model: EvaporationModel,
    /// The film state of the gas properties: `one-third` (default, the
    /// 1/3 rule with kinetic-theory mixture properties) or `far-field`
    /// (D321).
    pub film: Film,
    /// A sub-step ends when d² would fall by more than this fraction of its
    /// value at the start (default 0.1, D251).
    pub max_mass_change: f64,
    /// A sub-step ends when the droplet temperature would change by more
    /// than this (K, default 2).
    pub max_temperature_change: f64,
    /// A droplet smaller than this evaporates in one sub-step when its life
    /// ends inside it (m, default 1e-6).
    pub min_diameter: f64,
    /// The surface equilibrium: `raoult` (default) or `peng-robinson`
    /// (D311: high ambient pressure, ECN Spray A).
    pub equilibrium: Equilibrium,
    /// The data of `equilibrium = "peng-robinson"`.
    pub peng_robinson: PengRobinson,
}

impl Default for Evaporation {
    fn default() -> Self {
        Self {
            model: EvaporationModel::None,
            film: Film::OneThird,
            max_mass_change: 0.1,
            max_temperature_change: 2.0,
            min_diameter: 1e-6,
            equilibrium: Equilibrium::Raoult,
            peng_robinson: PengRobinson::default(),
        }
    }
}

impl Evaporation {
    /// True when the droplets evaporate.
    pub fn on(&self) -> bool {
        self.model != EvaporationModel::None
    }

    /// The kernel code of the model.
    pub fn code(&self) -> u32 {
        match self.model {
            EvaporationModel::AbramzonSirignano => {
                foehncraft_kernels::evaporation::model::ABRAMZON_SIRIGNANO
            }
            _ => foehncraft_kernels::evaporation::model::RANZ_MARSHALL,
        }
    }
}

/// One injector, `[[injector]]` (D246).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Injector {
    /// Name (in the log and the output).
    pub name: String,
    /// Where the parcels start (m). It must lie inside the mesh.
    pub position: [f64; 3],
    /// Axis of the spray (any length, it is normalised).
    pub direction: [f64; 3],
    /// Speed of the droplets at the injector (m/s). With `hole_diameter`
    /// it may be left out: the speed then follows from the mass flow,
    /// U = ṁ/(ρ_l C_a A_hole) (D306).
    #[serde(default)]
    pub speed: f64,
    /// Half angle of the spray cone (degrees); 0 is a straight jet. The
    /// directions are uniform in the cone's solid angle (counter-based random
    /// numbers keyed by the parcel identity, D121).
    #[serde(default)]
    pub half_angle: f64,
    /// Mass flow of the liquid (kg/s), constant; leave it out with `rate`.
    #[serde(default)]
    pub mass_flow: f64,
    /// Droplet diameter (m) when every droplet has the same size; leave it
    /// out with `size`.
    #[serde(default)]
    pub diameter: f64,
    /// Droplet size distribution (D242): Rosin-Rammler or log-normal, as
    /// the distribution of the liquid mass over the diameter. Every parcel
    /// carries the same mass; its diameter is drawn from the distribution
    /// with the counter-based generator keyed by the parcel identity
    /// (D121), so the spray is the same on every run and thread count.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size: Option<crate::sizes::SizeDistribution>,
    /// Parcels released per second. A parcel stands for the droplets of
    /// `mass_flow / parcels_per_second` of liquid. Needed in a transient
    /// case; the steady mode (D260) uses `[spray.steady] streams` instead.
    #[serde(default)]
    pub parcels_per_second: f64,
    /// Start of the injection (s).
    #[serde(default)]
    pub start: f64,
    /// End of the injection (s); default: the run's end.
    #[serde(default)]
    pub end: Option<f64>,
    /// Temperature of the injected droplets (K); default: the liquid's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature: Option<f64>,
    /// Diameter of the nozzle hole (m), for blob injection (D306): without
    /// `diameter` and `size` the parcels start as blobs of the diameter
    /// √C_a d_hole (Reitz and Diwakar 1987), and without `speed` they leave
    /// at U = ṁ/(ρ_l C_a A_hole). The breakup length of `[spray.breakup]`
    /// is measured in hole diameters.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hole_diameter: Option<f64>,
    /// Area contraction C_a of the jet at the hole exit (0 to 1, default 1;
    /// the effective flow area over the geometric one).
    #[serde(default = "one", skip_serializing_if = "is_one")]
    pub area_contraction: f64,
    /// The mass flow over time (D306): pairs [t, ṁ] (s from `start`,
    /// kg/s), linear in between and zero outside; replaces `mass_flow`. A
    /// parcel carries the exact integral of the table over its share of
    /// the injection time, and the injection ends at the table's last time
    /// unless `end` comes first.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub rate: Vec<[f64; 2]>,
    /// The share of the injector's liquid that enters the domain (0 to 1,
    /// default 1): a quarter model of a round jet about its axis takes
    /// 0.25. `mass_flow` and `rate` are the whole injector's; the parcels
    /// carry this share of them, and the speed from the hole follows from
    /// the whole flow.
    #[serde(default = "one", skip_serializing_if = "is_one")]
    pub flow_share: f64,
    /// Length of the nozzle hole (m), for the nozzle turbulence of
    /// `kh-act` (D307): with `discharge_coefficient` the force balance of
    /// Huh and Gosman (1991) estimates k₀ and ε₀ at the exit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hole_length: Option<f64>,
    /// Discharge coefficient C_d = C_a C_v of the hole (0 to 1, at most
    /// `area_contraction`): the mass flow over that of the ideal Bernoulli
    /// jet, for the nozzle turbulence estimate (D307).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub discharge_coefficient: Option<f64>,
    /// Loss coefficient K_c of the inlet of the hole, for the nozzle
    /// turbulence estimate (D307; default 0.5, the textbook value of a
    /// sharp-edged pipe entrance; a rounded inlet loses less).
    #[serde(default = "half", skip_serializing_if = "is_half")]
    pub inlet_loss: f64,
    /// Turbulent kinetic energy k₀ of the liquid at the hole exit
    /// (m²/s², with `nozzle_epsilon`), for example from a simulation of
    /// the nozzle flow; replaces the estimate (D307).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nozzle_k: Option<f64>,
    /// Its dissipation rate ε₀ (m²/s³, with `nozzle_k`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nozzle_epsilon: Option<f64>,
}

fn one() -> f64 {
    1.0
}

fn half() -> f64 {
    0.5
}

fn is_half(x: &f64) -> bool {
    *x == 0.5
}

/// The linear piece `w` of a rate table at `t`.
fn piece(w: &[[f64; 2]], t: f64) -> f64 {
    let ([t0, m0], [t1, m1]) = (w[0], w[1]);
    if t1 <= t0 {
        return m1;
    }
    m0 + (m1 - m0) * (t - t0) / (t1 - t0)
}

fn is_one(x: &f64) -> bool {
    *x == 1.0
}

impl Default for Injector {
    fn default() -> Self {
        Self {
            name: String::new(),
            position: [0.0; 3],
            direction: [1.0, 0.0, 0.0],
            speed: 0.0,
            half_angle: 0.0,
            mass_flow: 0.0,
            diameter: 0.0,
            size: None,
            parcels_per_second: 0.0,
            start: 0.0,
            end: None,
            temperature: None,
            hole_diameter: None,
            area_contraction: 1.0,
            rate: Vec::new(),
            flow_share: 1.0,
            hole_length: None,
            discharge_coefficient: None,
            inlet_loss: 0.5,
            nozzle_k: None,
            nozzle_epsilon: None,
        }
    }
}

impl Injector {
    /// The mass flow at time `t` since the start of the injection (kg/s).
    pub fn mass_flow_at(&self, t: f64) -> f64 {
        if self.rate.is_empty() {
            return self.mass_flow;
        }
        let r = &self.rate;
        if t < r[0][0] || t > r[r.len() - 1][0] {
            return 0.0;
        }
        for w in r.windows(2) {
            if t <= w[1][0] {
                return piece(w, t);
            }
        }
        r[r.len() - 1][1]
    }

    /// The liquid injected between `a` and `b` (s since the start of the
    /// injection, kg): the exact integral of the rate table (trapezoids of
    /// the linear pieces), or ṁ (b − a) at a constant mass flow.
    pub fn mass_between(&self, a: f64, b: f64) -> f64 {
        if self.rate.is_empty() {
            return self.mass_flow * (b - a);
        }
        let mut m = 0.0;
        for w in self.rate.windows(2) {
            let (t0, t1) = (w[0][0], w[1][0]);
            let lo = a.max(t0);
            let hi = b.min(t1);
            if hi > lo {
                m += 0.5 * (piece(w, lo) + piece(w, hi)) * (hi - lo);
            }
        }
        m
    }

    /// The end of the injection (s): `end`, else the last time of the rate
    /// table, else never.
    pub fn end_time(&self) -> f64 {
        let table = self
            .rate
            .last()
            .map_or(f64::INFINITY, |r| self.start + r[0]);
        self.end.unwrap_or(f64::INFINITY).min(table)
    }

    /// The flow area of the jet at the hole exit, C_a π d²/4 (m²), with a
    /// hole diameter.
    pub fn flow_area(&self) -> Option<f64> {
        self.hole_diameter
            .map(|d| self.area_contraction * std::f64::consts::PI * 0.25 * d * d)
    }

    /// The blob diameter √C_a d_hole (m), with a hole diameter.
    pub fn blob_diameter(&self) -> Option<f64> {
        self.hole_diameter.map(|d| self.area_contraction.sqrt() * d)
    }

    /// The turbulence of the liquid at the hole exit for the jet speed `u`
    /// (D307): `nozzle_k` and `nozzle_epsilon`, else the estimate of Huh
    /// and Gosman from the hole length and the discharge coefficient.
    pub fn nozzle_turbulence(&self, u: f64) -> Option<NozzleTurbulence> {
        if let (Some(k0), Some(eps0)) = (self.nozzle_k, self.nozzle_epsilon) {
            return Some(NozzleTurbulence { k0, eps0 });
        }
        NozzleTurbulence::huh_gosman(
            u,
            self.hole_length?,
            self.hole_diameter?,
            self.discharge_coefficient?,
            self.inlet_loss,
        )
    }

    /// The cavity of the hole (D307) for a liquid of density `rho_l` and
    /// vapour pressure `p_v`, when the jet is contracted (C_a < 1).
    pub fn cavitation(&self, rho_l: f64, p_v: Option<f64>) -> Option<Cavitation> {
        Cavitation::of(0.5 * self.hole_diameter?, self.area_contraction, rho_l, p_v)
    }
}

/// The breakup model of the drops (SP6, D300).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum BreakupModel {
    /// The drops keep their size (default).
    #[default]
    None,
    /// Kelvin-Helmholtz stripping of the injected blobs into child parcels
    /// and of the drops in place, Rayleigh-Taylor breakup beyond the
    /// breakup length (Reitz 1987; Patterson and Reitz 1998; Beale and
    /// Reitz 1999), D301, D302.
    KhRt,
    /// KH-RT whose blobs are also stripped by the turbulence and the
    /// cavitation the nozzle leaves in the liquid (Som and Aggarwal 2010,
    /// D307): the mechanism with the largest ratio of length to time scale
    /// strips the blob; needs the nozzle turbulence of each hole (D307).
    KhAct,
}

/// The drag law of the drops (D303).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum DragLaw {
    /// A rigid sphere (Schiller and Naumann 1935; default).
    #[default]
    Sphere,
    /// The drag of a distorted drop, C_D = C_D,sphere (1 + 2.632 y) with the
    /// TAB distortion y of every parcel (O'Rourke and Amsden 1987; Liu,
    /// Mather and Reitz 1993).
    Dynamic,
}

/// `[spray.breakup]` (SP6, D300 to D305): the breakup of the drops and the
/// drag of distorted drops. The constants are those of the model papers by
/// default; a validation case calibrates `b1` against its measurements.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Breakup {
    /// `none` (default), `kh-rt` or `kh-act`.
    pub model: BreakupModel,
    /// `sphere` (default) or `dynamic`.
    pub drag: DragLaw,
    /// KH child size constant B₀: r_c = B₀ Λ (default 0.61, Reitz 1987).
    pub b0: f64,
    /// KH time constant B₁: τ = 3.726 B₁ a/(Λ Ω) (default 10, Reitz
    /// 1987); it stands for the disturbance level the nozzle leaves on the
    /// jet and is the constant a validation case calibrates.
    pub b1: f64,
    /// RT size constant C_RT: the drops break into drops of the diameter
    /// C_RT Λ_RT (default 1).
    pub c_rt: f64,
    /// RT time constant C_τ: the wave grows for C_τ/Ω_RT (default 1).
    pub c_tau: f64,
    /// Breakup length constant C_bl: RT acts beyond L_b = C_bl (ρ_l/ρ_g)^½
    /// d_hole from the injector (default B₁/2, Beale and Reitz 1999).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub breakup_length: Option<f64>,
    /// A blob sheds a child parcel whenever its stripped liquid reaches
    /// this share of its injected liquid (default 0.03, Patterson and Reitz
    /// 1998).
    pub shed_fraction: f64,
    /// `kh-act`: the rate constant C_T,CAV of the turbulence and cavitation
    /// stripping, da/dt = −C_T,CAV L/τ (default 1, D307).
    pub c_tcav: f64,
}

impl Default for Breakup {
    fn default() -> Self {
        Self {
            model: BreakupModel::None,
            drag: DragLaw::Sphere,
            b0: 0.61,
            b1: 10.0,
            c_rt: 1.0,
            c_tau: 1.0,
            breakup_length: None,
            shed_fraction: 0.03,
            c_tcav: 1.0,
        }
    }
}

impl Breakup {
    /// True when the drops break up.
    pub fn on(&self) -> bool {
        self.model != BreakupModel::None
    }

    /// True when the nozzle turbulence and cavitation strip the blobs
    /// (KH-ACT, D307).
    pub fn act(&self) -> bool {
        self.model == BreakupModel::KhAct
    }

    /// True when the parcels carry the TAB distortion (dynamic drag).
    pub fn distortion(&self) -> bool {
        self.drag == DragLaw::Dynamic
    }

    /// True when the surface tension and viscosity are needed.
    pub fn needs_properties(&self) -> bool {
        self.on() || self.distortion()
    }

    /// The breakup length constant C_bl.
    pub fn c_bl(&self) -> f64 {
        self.breakup_length.unwrap_or(0.5 * self.b1)
    }
}

/// `[spray]`: the Lagrangian parcels of the case.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct SpraySection {
    /// A file with the injectors and the liquid (`[[injector]]` and
    /// `[liquid]`, D246), relative to the case file. The tables of the file
    /// are read into the case (`case.resolved.toml` holds them inline);
    /// defining `injector` or `liquid` in the case as well is an error
    /// (D114).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub include: Option<PathBuf>,
    /// The liquid of the parcels.
    pub liquid: Liquid,
    /// The injectors.
    #[serde(rename = "injector", skip_serializing_if = "Vec::is_empty")]
    pub injectors: Vec<Injector>,
    /// What a parcel does at a wall (default `stick`).
    pub wall_mode: BoundaryMode,
    /// Restitution of the normal velocity at `rebound` boundaries.
    pub restitution_normal: f64,
    /// Restitution of the tangential velocity at `rebound` boundaries.
    pub restitution_tangential: f64,
    /// Behaviour per boundary.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub boundary: BTreeMap<String, SprayBoundary>,
    /// Largest path of a parcel in one sub-step in cell sizes (the gas
    /// state is taken at the start of the sub-step; default 0.5).
    pub max_courant: f64,
    /// Largest number of sub-steps of one parcel in one time step (default
    /// 100000); a parcel that needs more is counted as lost.
    pub max_substeps: usize,
    /// A parcel file `spray/parcels_<step>.vtkhdf` every this many time
    /// steps (default 0: only at the end of the run).
    pub write_every: usize,
    /// Liquid penetration (D265): the share of an injector's liquid mass in
    /// the flow that lies within the penetration length (default 0.95).
    pub penetration_fraction: f64,
    /// How the penetration length is measured from the injector: `axial`
    /// (along the injector axis, default) or `radial` (straight distance).
    pub penetration: Penetration,
    /// The coupling of the parcels with the gas (SP2, D247 to D252).
    pub coupling: Coupling,
    /// Turbulent dispersion of the parcels (D243): `none` (default) or
    /// `eddy-lifetime` (Gosman and Ioannides 1983); `random-walk` follows
    /// behind the same key. Needs a turbulence model.
    pub dispersion: Dispersion,
    /// The steady spray mode (D260, SP2s): the keys of a case without
    /// `[time]` (the defaults when the table is left out); refused in a
    /// transient case.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub steady: Option<Steady>,
    /// The evaporation of the droplets (SP5).
    pub evaporation: Evaporation,
    /// The breakup of the drops and the drag of distorted drops (SP6).
    pub breakup: Breakup,
    /// The diagnostics of the Engine Combustion Network modelling standard
    /// (D315): liquid length from the projected liquid volume and vapour
    /// penetration, per injector in `spray.csv` and `result.json`; off
    /// without the table.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ecn: Option<Ecn>,
}

/// `[spray.ecn]`: the liquid length and the vapour penetration as the ECN
/// modelling standard defines them (ecn.sandia.gov, Modeling Standards,
/// D315). Liquid length: the farthest axial position where the projected
/// liquid volume PLV = ∫ LVF dy across the spray (along the line of sight)
/// falls to a threshold. Vapour penetration: the farthest axial distance
/// from the nozzle where the fuel mass fraction is `vapour_fraction`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Ecn {
    /// PLV thresholds (m, liquid volume per projected area: 2e-7 m is the
    /// ECN 0.2·10⁻³ mm³/mm²); default `[2e-7, 2e-6]`, the two values of the
    /// standard. One liquid length per threshold.
    pub plv: Vec<f64>,
    /// Size of the square bins of the projection (m), along the axis and
    /// across it (default 2.5e-4, the smallest cell the standard
    /// recommends for RANS Lagrangian sprays).
    pub bin: f64,
    /// The line of sight (any direction not along the injector axis; its
    /// part across the axis is used); default: the coordinate axis most
    /// nearly perpendicular to the injector axis.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub line_of_sight: Option<[f64; 3]>,
    /// Fuel mass fraction of the vapour penetration (default 1e-3, the
    /// 0.1 % of the standard); needs evaporation.
    pub vapour_fraction: f64,
}

impl Default for Ecn {
    fn default() -> Self {
        Self {
            plv: vec![2e-7, 2e-6],
            bin: 2.5e-4,
            line_of_sight: None,
            vapour_fraction: 1e-3,
        }
    }
}

/// `[spray.steady]`: the steady spray mode (D260, SP2s). A parcel is a
/// *stream*: it stands for a number of droplets per second, starts at its
/// injector and is followed through the frozen steady gas until it leaves,
/// sticks or reaches `max_time`. With two-way coupling a stream pass runs
/// every `every` outer iterations and the gas takes the sources relaxed
/// between passes, S ← S + α (S_new − S).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Steady {
    /// Streams per injector (default 200). Their directions and sizes are
    /// drawn as those of the parcels of a transient run (keyed by the
    /// stream, D121), so every pass follows the same streams.
    pub streams: u32,
    /// A stream pass every this many outer iterations (default 10).
    pub every: usize,
    /// Relaxation α of the sources between passes (default 0.5; 1: the
    /// sources of the last pass as they are).
    pub relaxation: f64,
    /// The source residual below which the spray counts as converged
    /// (default 1e-3): the difference between the relaxed sources the gas
    /// holds and those of the last pass, over the sources of the last pass
    /// (L1 sums over the spray cells, momentum and heat, the larger).
    pub tolerance: f64,
    /// Longest time a stream is followed (s, default 10); a stream still
    /// in the flow then is counted as unfinished (SPR007).
    pub max_time: f64,
    /// The key of the streams (default 0): another seed draws other
    /// directions, sizes and eddies, for the statistical error of the one
    /// sample (D260 (3)).
    pub seed: u64,
}

impl Default for Steady {
    fn default() -> Self {
        Self {
            streams: 200,
            every: 10,
            relaxation: 0.5,
            tolerance: 1e-3,
            max_time: 10.0,
            seed: 0,
        }
    }
}

/// The turbulent dispersion model of the parcels (D243).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Dispersion {
    /// The parcels see the mean gas velocity.
    #[default]
    None,
    /// Eddy lifetime (Gosman and Ioannides 1983): a parcel sees the mean
    /// gas velocity plus an isotropic Gaussian fluctuation of variance
    /// 2k/3 per component, kept for the eddy interaction time, the shorter
    /// of the eddy lifetime and the time the parcel needs to cross the
    /// eddy.
    EddyLifetime,
    /// Continuous random walk (D243, answered 2026-11-24): follows; refused
    /// in this version.
    RandomWalk,
}

/// Whether the gas feels the parcels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum CouplingMode {
    /// The parcels feel the gas, the gas does not feel them (SP1).
    #[default]
    OneWay,
    /// Momentum exchange in both directions (SP2): the drag on the parcels
    /// is a source of the gas momentum, in value-plus-slope form with the
    /// closing step (D250) and passes inside the time step (D252).
    TwoWay,
}

/// The number of parcel passes in a time step (D252): a number, or
/// `"converged"` (passes until the coupling lag is below `tolerance`, at
/// most `max_passes`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(untagged)]
pub enum Passes {
    /// A fixed number of passes (at least 1).
    Fixed(u32),
    /// `"converged"`.
    Named(PassesName),
}

/// The named form of [`Passes`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum PassesName {
    /// Passes until the coupling lag is below the tolerance.
    Converged,
}

impl Default for Passes {
    fn default() -> Self {
        Self::Named(PassesName::Converged)
    }
}

/// `[spray.coupling]` (D252, decided by the maintainer 2026-11-24: passes
/// = "converged", max_passes = 3).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Coupling {
    /// `one-way` (default) or `two-way`.
    pub mode: CouplingMode,
    /// Parcel passes per time step: a number, or `"converged"`.
    pub passes: Passes,
    /// The coupling lag below which the passes of a step end (largest
    /// change of the gas velocity of a spray cell between the frozen state
    /// of the pass and the solve, over the larger of the largest slip
    /// velocity of the spray cells and the largest change of their gas
    /// velocity in the step; default 1e-3).
    pub tolerance: f64,
    /// The largest number of passes with `"converged"` (default 3).
    pub max_passes: u32,
    /// A sub-step ends when the particle Reynolds number would change by
    /// more than this fraction of max(Re, 1) (D251; default 0.2).
    pub max_reynolds_change: f64,
}

impl Default for Coupling {
    fn default() -> Self {
        Self {
            mode: CouplingMode::OneWay,
            passes: Passes::default(),
            tolerance: 1e-3,
            max_passes: 3,
            max_reynolds_change: 0.2,
        }
    }
}

impl Coupling {
    /// True with two-way coupling.
    pub fn two_way(&self) -> bool {
        self.mode == CouplingMode::TwoWay
    }

    /// The largest number of passes of a step, and whether they end on the
    /// coupling lag.
    pub fn pass_limit(&self) -> (u32, bool) {
        match self.passes {
            Passes::Fixed(n) => (n, false),
            Passes::Named(PassesName::Converged) => (self.max_passes, true),
        }
    }
}

/// How the liquid penetration is measured (D265).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Penetration {
    /// Distance along the injector axis.
    #[default]
    Axial,
    /// Straight distance from the injector position.
    Radial,
}

/// The file of `include`: the injectors and the liquid.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct SprayFile {
    liquid: Option<Liquid>,
    #[serde(rename = "injector")]
    injectors: Vec<Injector>,
}

impl Default for SpraySection {
    fn default() -> Self {
        Self {
            include: None,
            liquid: Liquid::default(),
            injectors: Vec::new(),
            wall_mode: BoundaryMode::Stick,
            restitution_normal: 1.0,
            restitution_tangential: 1.0,
            boundary: BTreeMap::new(),
            max_courant: 0.5,
            max_substeps: 100_000,
            write_every: 0,
            penetration_fraction: 0.95,
            penetration: Penetration::Axial,
            coupling: Coupling::default(),
            dispersion: Dispersion::None,
            steady: None,
            evaporation: Evaporation::default(),
            breakup: Breakup::default(),
            ecn: None,
        }
    }
}

impl SpraySection {
    /// True when the case has no spray.
    pub fn is_empty(&self) -> bool {
        self.include.is_none() && self.injectors.is_empty()
    }

    /// Resolves the include path against the case directory.
    pub fn resolve_paths(&mut self, dir: &Path) {
        if let Some(f) = &mut self.include {
            *f = dir.join(&*f);
        }
    }

    /// Reads the include file into the section (D246) and checks the values.
    pub fn load(&mut self) -> Result<()> {
        if let Some(path) = self.include.take() {
            let text = std::fs::read_to_string(&path)
                .map_err(|e| anyhow::anyhow!("spray file {}: {e}", path.display()))?;
            let f: SprayFile = toml::from_str(&text)
                .map_err(|e| anyhow::anyhow!("spray file {}: {e}", path.display()))?;
            if !f.injectors.is_empty() && !self.injectors.is_empty() {
                bail!(
                    "[spray] injector is defined in the case and in {} (a key defined twice is an error, D114)",
                    path.display()
                );
            }
            if f.liquid.is_some() && self.liquid != Liquid::default() {
                bail!(
                    "[spray.liquid] is defined in the case and in {} (a key defined twice is an error, D114)",
                    path.display()
                );
            }
            if let Some(l) = f.liquid {
                self.liquid = l;
            }
            self.injectors.extend(f.injectors);
        }
        self.liquid.apply_library()?;
        self.validate()
    }

    /// The nozzle keys of the KH-ACT breakup of injector `i` named `n`
    /// (D307).
    fn validate_nozzle(&self, i: &Injector, n: &str) -> Result<()> {
        let pos = |x: f64| x.is_finite() && x > 0.0;
        for (k, v) in [
            ("hole_length", i.hole_length),
            ("nozzle_k", i.nozzle_k),
            ("nozzle_epsilon", i.nozzle_epsilon),
        ] {
            if let Some(v) = v
                && !pos(v)
            {
                bail!("[[spray.injector]] {n}: {k} must be positive and finite");
            }
        }
        if let Some(cd) = i.discharge_coefficient
            && !(cd > 0.0 && cd <= i.area_contraction)
        {
            bail!(
                "[[spray.injector]] {n}: discharge_coefficient must be in (0, area_contraction] \
                 (C_d = C_a C_v with C_v at most 1; area_contraction is {})",
                i.area_contraction
            );
        }
        if !(i.inlet_loss.is_finite() && i.inlet_loss >= 0.0) {
            bail!("[[spray.injector]] {n}: inlet_loss must be zero or positive");
        }
        if i.nozzle_k.is_some() != i.nozzle_epsilon.is_some() {
            bail!("[[spray.injector]] {n}: give nozzle_k and nozzle_epsilon together");
        }
        let given = i.nozzle_k.is_some();
        let estimate = i.hole_length.is_some() || i.discharge_coefficient.is_some();
        if given && estimate {
            bail!(
                "[[spray.injector]] {n}: give either nozzle_k and nozzle_epsilon or \
                 hole_length and discharge_coefficient, not both"
            );
        }
        if estimate && (i.hole_length.is_none() || i.discharge_coefficient.is_none()) {
            bail!(
                "[[spray.injector]] {n}: the nozzle turbulence estimate needs both hole_length \
                 and discharge_coefficient"
            );
        }
        if (given || estimate) && i.hole_diameter.is_none() {
            bail!("[[spray.injector]] {n}: the nozzle turbulence needs hole_diameter");
        }
        // the nozzle describes the hole, whatever the model: kh-rt leaves
        // its turbulence unused
        if !self.breakup.act() {
            return Ok(());
        }
        if i.hole_diameter.is_some() && !given && !estimate {
            bail!(
                "[[spray.injector]] {n}: kh-act needs the turbulence of the liquid at the hole \
                 exit: hole_length and discharge_coefficient (the estimate of Huh and Gosman), \
                 or nozzle_k and nozzle_epsilon"
            );
        }
        if let (Some(cd), true) = (i.discharge_coefficient, estimate)
            && 1.0 / (cd * cd) - i.inlet_loss - 1.0 <= 0.0
        {
            bail!(
                "[[spray.injector]] {n}: with discharge_coefficient {cd} and inlet_loss {} the \
                 losses leave no wall friction for the nozzle turbulence estimate \
                 (1/C_d² − K_c − 1 must be positive); give nozzle_k and nozzle_epsilon",
                i.inlet_loss
            );
        }
        Ok(())
    }

    /// Checks the values.
    pub fn validate(&self) -> Result<()> {
        let pos = |x: f64| x.is_finite() && x > 0.0;
        if self.dispersion == Dispersion::RandomWalk {
            bail!(
                "[spray] dispersion = \"random-walk\" is not available in this version; use \"eddy-lifetime\""
            );
        }
        if matches!(self.liquid.saturation, Saturation::Library) && self.liquid.data().is_none() {
            bail!("[spray.liquid] saturation = \"library\" needs a library liquid (`library`)");
        }
        if !pos(self.liquid.density) {
            bail!("[spray.liquid] density must be positive and finite");
        }
        if !pos(self.liquid.specific_heat) {
            bail!("[spray.liquid] specific_heat must be positive and finite");
        }
        if !pos(self.liquid.temperature) {
            bail!("[spray.liquid] temperature must be positive and finite (K)");
        }
        if !pos(self.liquid.latent_heat) {
            bail!("[spray.liquid] latent_heat must be positive and finite (J/kg)");
        }
        if !pos(self.liquid.latent_heat_temperature) {
            bail!("[spray.liquid] latent_heat_temperature must be positive and finite (K)");
        }
        let ev = &self.evaporation;
        if !(ev.max_mass_change > 0.0 && ev.max_mass_change < 1.0) {
            bail!("[spray.evaporation] max_mass_change must be in (0, 1)");
        }
        if !pos(ev.max_temperature_change) {
            bail!("[spray.evaporation] max_temperature_change must be positive (K)");
        }
        if !(ev.min_diameter.is_finite() && ev.min_diameter >= 0.0) {
            bail!("[spray.evaporation] min_diameter must be zero or positive (m)");
        }
        if ev.on() && self.liquid.vapour.is_none() {
            bail!(
                "[spray.evaporation] needs the species the liquid evaporates into: \
                 [spray.liquid] vapour = \"H2O\" (a species of [species])"
            );
        }
        let br = &self.breakup;
        for (k, v) in [
            ("b0", br.b0),
            ("b1", br.b1),
            ("c_rt", br.c_rt),
            ("c_tau", br.c_tau),
            ("c_tcav", br.c_tcav),
        ] {
            if !pos(v) {
                bail!("[spray.breakup] {k} must be positive and finite");
            }
        }
        if let Some(c) = br.breakup_length
            && !(c.is_finite() && c >= 0.0)
        {
            bail!("[spray.breakup] breakup_length must be zero or positive");
        }
        if !(br.shed_fraction >= 1e-4 && br.shed_fraction < 1.0) {
            bail!("[spray.breakup] shed_fraction must be in [1e-4, 1)");
        }
        if br.needs_properties() {
            let what = if br.on() {
                "the breakup"
            } else {
                "the dynamic drag"
            };
            if !self.liquid.surface_tension.is_some_and(pos) {
                bail!(
                    "[spray.breakup] {what} needs the surface tension of the liquid: \
                     [spray.liquid] surface_tension = 0.02 (N/m, positive)"
                );
            }
            if !self.liquid.viscosity.is_some_and(pos) {
                bail!(
                    "[spray.breakup] {what} needs the viscosity of the liquid: \
                     [spray.liquid] viscosity = 1e-3 (Pa s, positive)"
                );
            }
            if self.steady.is_some() {
                bail!(
                    "[spray.breakup] is not available in the steady spray mode in this version \
                     (a stream cannot shed child streams yet); run the case transient"
                );
            }
        }
        for (k, v) in [
            ("surface_tension", self.liquid.surface_tension),
            ("viscosity", self.liquid.viscosity),
            ("vapour_pressure", self.liquid.vapour_pressure),
        ] {
            if let Some(v) = v
                && !pos(v)
            {
                bail!("[spray.liquid] {k} must be positive and finite");
            }
        }
        for inj in &self.injectors {
            if let Some(t) = inj.temperature
                && !pos(t)
            {
                bail!(
                    "[[injector]] `{}`: temperature must be positive and finite (K)",
                    inj.name
                );
            }
        }
        for (k, v) in [
            ("restitution_normal", self.restitution_normal),
            ("restitution_tangential", self.restitution_tangential),
        ] {
            if !(0.0..=1.0).contains(&v) {
                bail!("[spray] {k} must be between 0 and 1");
            }
        }
        if !pos(self.max_courant) {
            bail!("[spray] max_courant must be positive");
        }
        if !(self.penetration_fraction > 0.0 && self.penetration_fraction <= 1.0) {
            bail!("[spray] penetration_fraction must be in (0, 1]");
        }
        if let Some(e) = &self.ecn {
            if self.steady.is_some() {
                bail!(
                    "[spray.ecn] needs the transient spray: the streams of [spray.steady] carry \
                     rates, not the liquid in the flow"
                );
            }
            if e.plv.iter().any(|&v| !pos(v)) {
                bail!("[spray.ecn] plv thresholds must be positive (m, 2e-7 = 0.2e-3 mm3/mm2)");
            }
            if !pos(e.bin) {
                bail!("[spray.ecn] bin must be positive (m)");
            }
            if !(e.vapour_fraction > 0.0 && e.vapour_fraction < 1.0) {
                bail!("[spray.ecn] vapour_fraction must be in (0, 1)");
            }
            if let Some(l) = e.line_of_sight {
                for inj in &self.injectors {
                    let a = inj.direction;
                    let na = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
                    let nl = (l[0] * l[0] + l[1] * l[1] + l[2] * l[2]).sqrt();
                    let c = (a[0] * l[0] + a[1] * l[1] + a[2] * l[2]) / (na * nl);
                    if !(nl > 0.0 && c.abs() < 0.99) {
                        bail!(
                            "[spray.ecn] line_of_sight must not lie along the axis of injector `{}`",
                            inj.name
                        );
                    }
                }
            }
        }
        let c = &self.coupling;
        if let Passes::Fixed(0) = c.passes {
            bail!("[spray.coupling] passes must be at least 1");
        }
        if c.max_passes == 0 {
            bail!("[spray.coupling] max_passes must be at least 1");
        }
        if !pos(c.tolerance) {
            bail!("[spray.coupling] tolerance must be positive");
        }
        if !pos(c.max_reynolds_change) {
            bail!("[spray.coupling] max_reynolds_change must be positive");
        }
        if self.max_substeps == 0 {
            bail!("[spray] max_substeps must be at least 1");
        }
        if let Some(st) = &self.steady {
            if st.streams == 0 {
                bail!("[spray.steady] streams must be at least 1");
            }
            if st.every == 0 {
                bail!("[spray.steady] every must be at least 1");
            }
            if !(st.relaxation > 0.0 && st.relaxation <= 1.0) {
                bail!("[spray.steady] relaxation must be in (0, 1]");
            }
            if !pos(st.tolerance) {
                bail!("[spray.steady] tolerance must be positive");
            }
            if !pos(st.max_time) {
                bail!("[spray.steady] max_time must be positive and finite (s)");
            }
        }
        let mut names = std::collections::BTreeSet::new();
        for i in &self.injectors {
            let n = &i.name;
            if !names.insert(n.clone()) {
                bail!("[[spray.injector]] {n}: the name is used twice");
            }
            let dir2: f64 = i.direction.iter().map(|x| x * x).sum();
            if !(dir2.is_finite() && dir2 > 0.0) {
                bail!("[[spray.injector]] {n}: direction must not be zero");
            }
            if !i.position.iter().all(|x| x.is_finite()) {
                bail!("[[spray.injector]] {n}: position must be finite");
            }
            if !(i.speed.is_finite() && i.speed >= 0.0) {
                bail!("[[spray.injector]] {n}: speed must be zero or positive");
            }
            if let Some(h) = i.hole_diameter
                && !pos(h)
            {
                bail!("[[spray.injector]] {n}: hole_diameter must be positive and finite (m)");
            }
            if !(i.area_contraction > 0.0 && i.area_contraction <= 1.0) {
                bail!("[[spray.injector]] {n}: area_contraction must be in (0, 1]");
            }
            if !(i.flow_share > 0.0 && i.flow_share <= 1.0) {
                bail!("[[spray.injector]] {n}: flow_share must be in (0, 1]");
            }
            if i.area_contraction != 1.0 && i.hole_diameter.is_none() {
                bail!("[[spray.injector]] {n}: area_contraction needs hole_diameter");
            }
            self.validate_nozzle(i, n)?;
            if !i.rate.is_empty() {
                if i.mass_flow != 0.0 {
                    bail!("[[spray.injector]] {n}: give either mass_flow or rate, not both");
                }
                if i.rate.len() < 2 {
                    bail!("[[spray.injector]] {n}: rate needs at least two points [t, kg/s]");
                }
                for w in i.rate.windows(2) {
                    if w[1][0] < w[0][0] {
                        bail!("[[spray.injector]] {n}: the times of rate must not decrease");
                    }
                }
                if i.rate
                    .iter()
                    .any(|r| !(r[0].is_finite() && r[1].is_finite() && r[1] >= 0.0))
                {
                    bail!(
                        "[[spray.injector]] {n}: rate holds finite times and mass flows of zero or more"
                    );
                }
                if i.rate[0][0] < 0.0 {
                    bail!(
                        "[[spray.injector]] {n}: the times of rate start at zero or later (s from start)"
                    );
                }
                let total = i.mass_between(i.rate[0][0], i.rate[i.rate.len() - 1][0]);
                if total.is_nan() || total <= 0.0 {
                    bail!("[[spray.injector]] {n}: rate injects no liquid");
                }
            }

            if !(0.0..180.0).contains(&i.half_angle) {
                bail!("[[spray.injector]] {n}: half_angle must be in [0, 180) degrees");
            }
            match &i.size {
                Some(sz) => {
                    if i.diameter != 0.0 {
                        bail!("[[spray.injector]] {n}: give either diameter or size, not both");
                    }
                    sz.validate()
                        .map_err(|e| anyhow::anyhow!("[[spray.injector]] {n}: size: {e}"))?;
                }
                None => {
                    if !pos(i.diameter) && i.hole_diameter.is_none() {
                        bail!(
                            "[[spray.injector]] {n}: diameter must be positive and finite (or give a size distribution, or hole_diameter for blobs)"
                        );
                    }
                }
            }
            if i.rate.is_empty() && !pos(i.mass_flow) {
                bail!(
                    "[[spray.injector]] {n}: mass_flow must be positive and finite (or give rate)"
                );
            }
            if !i.rate.is_empty() && self.steady.is_some() {
                bail!(
                    "[[spray.injector]] {n}: rate is a transient injection; the steady mode takes mass_flow"
                );
            }
            // zero: a steady case (streams instead of parcels); a
            // transient case needs it (checked with the case)
            if !(i.parcels_per_second.is_finite() && i.parcels_per_second >= 0.0) {
                bail!("[[spray.injector]] {n}: parcels_per_second must be zero or positive");
            }
            if !i.start.is_finite() || i.start < 0.0 {
                bail!("[[spray.injector]] {n}: start must be zero or positive");
            }
            if let Some(e) = i.end
                && !(e.is_finite() && e > i.start)
            {
                bail!("[[spray.injector]] {n}: end must be after start");
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod library_tests {
    use super::*;

    #[test]
    fn a_library_liquid_takes_its_properties_and_checks_its_keys() {
        let mut l = Liquid {
            library: Some("n-heptane".into()),
            temperature: 300.0,
            vapour: Some("NC7H16".into()),
            ..Liquid::default()
        };
        l.apply_library().unwrap();
        let d = &foehncraft_kernels::liquids::N_HEPTANE;
        assert_eq!(l.density, d.density(300.0));
        assert_eq!(l.saturation, Saturation::Library);
        assert_eq!(l.name, "n-heptane");
        assert_eq!(l.surface_tension, Some(d.surface_tension(300.0)));
        assert_eq!(l.viscosity, Some(d.viscosity(300.0)));
        assert_eq!(l.vapour_pressure, Some(d.saturation_pressure(300.0)));
        // applying twice (a re-read of the case) keeps the values
        let once = l.clone();
        l.apply_library().unwrap();
        assert_eq!(l, once);
        // water keeps its IAPWS-95 vapour pressure
        let mut w = Liquid {
            library: Some("water".into()),
            ..Liquid::default()
        };
        w.apply_library().unwrap();
        assert_eq!(w.saturation, Saturation::Water);
        w.apply_library().unwrap();
        // refused: an unknown liquid (the names are listed), a value given
        // beside the library, a property temperature without a library
        let e = Liquid {
            library: Some("benzene".into()),
            ..Liquid::default()
        }
        .apply_library()
        .unwrap_err()
        .to_string();
        assert!(e.contains("iso-octane") && e.contains("n-dodecane"), "{e}");
        let e = Liquid {
            library: Some("n-heptane".into()),
            density: 700.0,
            ..Liquid::default()
        }
        .apply_library()
        .unwrap_err()
        .to_string();
        assert!(e.contains("density"), "{e}");
        assert!(
            Liquid {
                property_temperature: Some(320.0),
                ..Liquid::default()
            }
            .apply_library()
            .is_err()
        );
    }
}
