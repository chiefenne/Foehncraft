//! `spray.devices`: the injector and nozzle devices of a case or an
//! injector file, their frames and holes in case coordinates (D346).

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::engine;
use foehncraft_particles::SpraySection;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Input of `spray.devices`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DevicesInput {
    /// A flow case file (TOML) with `[spray]`, or an injector file (the
    /// file of `[spray] include`: `[liquid]`, `[[injector]]`, `[[device]]`).
    pub file: PathBuf,
    /// Case keys to override, `key=value` (D62; case files only).
    #[serde(default)]
    pub set: Vec<String>,
}

/// One hole in the output of `spray.devices`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct HoleOut {
    /// Name of the hole.
    pub name: String,
    /// Centre of the exit (m).
    pub exit: [f64; 3],
    /// Unit direction of the hole axis.
    pub direction: [f64; 3],
    /// Angle from the device axis (degrees).
    pub polar: f64,
    /// Azimuth about the device axis (degrees).
    pub azimuth: f64,
    /// Diameter (m), when the device gives it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub diameter: Option<f64>,
    /// The injectors that name this hole.
    pub injectors: Vec<String>,
}

/// One device in the output of `spray.devices`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct DeviceOut {
    /// Name of the device.
    pub name: String,
    /// The surface (`none`, `standard-gasoline`, `standard-diesel`, `stl`).
    pub geometry: String,
    /// Triangles of the placed surface.
    pub triangles: usize,
    /// The mounting point (m).
    pub mount: [f64; 3],
    /// Zero of the azimuth (unit).
    pub e1: [f64; 3],
    /// e₃ × e₁ (unit).
    pub e2: [f64; 3],
    /// The injector axis (unit).
    pub axis: [f64; 3],
    /// The holes, in the order of the file (holes, then patterns).
    pub holes: Vec<HoleOut>,
}

/// One injector in the output of `spray.devices`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct InjectorOut {
    /// Name (`<name>/<hole>` for an injector of `hole = "all"`).
    pub name: String,
    /// Its device, if it names one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub device: Option<String>,
    /// Its hole, if it names one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hole: Option<String>,
    /// Where the parcels start (m).
    pub position: [f64; 3],
    /// The spray axis (as given or from the hole).
    pub direction: [f64; 3],
}

/// Output of `spray.devices`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct DevicesOutput {
    /// The devices.
    pub devices: Vec<DeviceOut>,
    /// The injectors after the expansion of `hole = "all"`, in the order of
    /// the run (the `injector` index of the parcel files).
    pub injectors: Vec<InjectorOut>,
    /// The SPR010 warnings the run would give.
    pub warnings: Vec<String>,
}

/// The output of `spray.devices` for a loaded `[spray]` section.
pub fn describe(sp: &SpraySection) -> Result<DevicesOutput> {
    let mut devices = Vec::new();
    for d in &sp.devices {
        let f = d.frame()?;
        let holes = d
            .placed_holes()?
            .into_iter()
            .map(|h| HoleOut {
                injectors: sp
                    .injectors
                    .iter()
                    .filter(|i| {
                        i.device.as_deref() == Some(&d.name) && i.hole.as_deref() == Some(&h.name)
                    })
                    .map(|i| i.name.clone())
                    .collect(),
                name: h.name,
                exit: h.exit,
                direction: h.direction,
                polar: h.polar,
                azimuth: h.azimuth,
                diameter: h.diameter,
            })
            .collect();
        devices.push(DeviceOut {
            name: d.name.clone(),
            geometry: serde_json::to_value(d.geometry)?
                .as_str()
                .unwrap_or("none")
                .to_string(),
            triangles: d.surface()?.len(),
            mount: f.origin,
            e1: f.e1,
            e2: f.e2,
            axis: f.e3,
            holes,
        });
    }
    let injectors = sp
        .injectors
        .iter()
        .map(|i| InjectorOut {
            name: i.name.clone(),
            device: i.device.clone(),
            hole: i.hole.clone(),
            position: i.position,
            direction: i.direction,
        })
        .collect();
    Ok(DevicesOutput {
        devices,
        injectors,
        warnings: sp.device_warnings.clone(),
    })
}

fn section(a: &DevicesInput, ctx: &Ctx<'_>) -> Result<SpraySection> {
    let roots = ctx.roots.clone();
    roots.check_read(&a.file, "the case or injector file")?;
    let text = std::fs::read_to_string(&a.file)?;
    let table: toml::Table = toml::from_str(&text)?;
    Ok(
        if table.contains_key("spray") || table.contains_key("mesh") {
            engine::load(&a.file, &a.set, &roots)?.spray
        } else {
            anyhow::ensure!(a.set.is_empty(), "set applies to case files only");
            SpraySection::from_file(&a.file)?
        },
    )
}

/// Input of `spray.place`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PlaceInput {
    /// A flow case file (TOML) with `[spray]`, or an injector file.
    pub file: PathBuf,
    /// Folder to write `spray/devices.vtu` into (as a run would).
    pub out: PathBuf,
    /// Case keys to override, `key=value` (D62; case files only).
    #[serde(default)]
    pub set: Vec<String>,
}

/// `spray.place`.
pub struct Place;

impl Tool for Place {
    const ID: &'static str = "spray.place";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Write the injector devices of a case or an injector file as a run would, `<out>/spray/devices.vtu` (placed surfaces, hole axes, injector cones; D343, D345), to look at a mounting before meshing or running; returns the devices as `spray.devices` does.";
    const CATEGORY: Category = Category::Spray;
    const ACTION: Action = Action::WriteFiles;
    type Input = PlaceInput;
    type Output = serde_json::Value;

    fn call(a: PlaceInput, ctx: &mut Ctx<'_>) -> Result<serde_json::Value> {
        let d = DevicesInput {
            file: a.file.clone(),
            set: a.set.clone(),
        };
        let sp = section(&d, ctx)?;
        ctx.roots.check_write(&a.out)?;
        let path = foehncraft_flow::spray::write_devices(&sp, &a.out)?;
        let out = describe(&sp)?;
        match &path {
            Some(p) => ctx.line(&format!("wrote {}", p.display())),
            None => ctx.line("no [[device]]: nothing written"),
        }
        let mut v = serde_json::to_value(&out)?;
        v["file"] = serde_json::to_value(&path)?;
        Ok(v)
    }
}

/// `spray.devices`.
pub struct Devices;

impl Tool for Devices {
    const ID: &'static str = "spray.devices";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "The injector and nozzle devices of a case or an injector file: each device's frame (mount, e1, e2, axis), every hole's exit and direction in case coordinates, the injectors after `hole = \"all\"` and the SPR010 warnings, without running the flow (D340 to D342, D346).";
    const CATEGORY: Category = Category::Spray;
    const ACTION: Action = Action::Read;
    type Input = DevicesInput;
    type Output = DevicesOutput;

    fn call(a: DevicesInput, ctx: &mut Ctx<'_>) -> Result<DevicesOutput> {
        let sp = section(&a, ctx)?;
        let out = describe(&sp)?;
        for d in &out.devices {
            ctx.line(&format!(
                "device {}: mount {:?}, axis {:?}, {} holes",
                d.name,
                d.mount,
                d.axis,
                d.holes.len()
            ));
            for h in &d.holes {
                ctx.line(&format!(
                    "  hole {}: exit {:?}, direction {:?}, injectors {:?}",
                    h.name, h.exit, h.direction, h.injectors
                ));
            }
        }
        for w in &out.warnings {
            ctx.line(&format!("SPR010: {w}"));
        }
        Ok(out)
    }
}
