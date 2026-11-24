//! The `[spray]` table of the case file (L12, D240 to D246).

use anyhow::{Result, bail};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

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
    /// velocity changes sign, nothing is lost.
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

/// The liquid of the parcels. SP1 takes a constant density; the property
/// set of the liquids (vapour pressure, latent heat, …) comes with the
/// evaporation (SP4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Liquid {
    /// Name, for the log and the output.
    pub name: String,
    /// Density (kg/m³).
    pub density: f64,
}

impl Default for Liquid {
    fn default() -> Self {
        Self {
            name: "water".into(),
            density: 998.0,
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
    /// Speed of the droplets at the injector (m/s).
    pub speed: f64,
    /// Half angle of the spray cone (degrees); 0 is a straight jet. The
    /// directions are uniform in the cone's solid angle (counter-based random
    /// numbers keyed by the parcel identity, D121).
    #[serde(default)]
    pub half_angle: f64,
    /// Mass flow of the liquid (kg/s).
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
    /// `mass_flow / parcels_per_second` of liquid.
    pub parcels_per_second: f64,
    /// Start of the injection (s).
    #[serde(default)]
    pub start: f64,
    /// End of the injection (s); default: the run's end.
    #[serde(default)]
    pub end: Option<f64>,
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
        self.validate()
    }

    /// Checks the values.
    pub fn validate(&self) -> Result<()> {
        let pos = |x: f64| x.is_finite() && x > 0.0;
        if !pos(self.liquid.density) {
            bail!("[spray.liquid] density must be positive and finite");
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
                    if !pos(i.diameter) {
                        bail!(
                            "[[spray.injector]] {n}: diameter must be positive and finite (or give a size distribution)"
                        );
                    }
                }
            }
            for (k, v) in [
                ("mass_flow", i.mass_flow),
                ("parcels_per_second", i.parcels_per_second),
            ] {
                if !pos(v) {
                    bail!("[[spray.injector]] {n}: {k} must be positive and finite");
                }
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
