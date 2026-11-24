//! The `[spray]` table of the case file (L12, D240 to D246).

use anyhow::{Result, bail};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// What a parcel does when it reaches a boundary.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema,
)]
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
    /// Droplet diameter (m): constant in this version.
    pub diameter: f64,
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
            for (k, v) in [
                ("mass_flow", i.mass_flow),
                ("diameter", i.diameter),
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
