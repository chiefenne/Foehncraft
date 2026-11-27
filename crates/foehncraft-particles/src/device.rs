//! Injector and nozzle devices (D340 to D342): the frame of a mounted
//! injector and its holes given in that frame.
//!
//! A device is a right-handed frame: the mounting point `mount` (the tip of
//! the injector), the injector `axis` into the chamber (e₃) and a
//! `reference` direction whose part normal to the axis is the zero of the
//! azimuth (e₁); e₂ = e₃ × e₁. `clocking` turns the whole device about its
//! axis. A hole is given by its polar angle θ from the axis, its azimuth φ,
//! the radius of its exit centre from the axis and the setback of the exit
//! along the axis:
//!
//! - exit x = mount + setback e₃ + radius (cos φₑ e₁ + sin φₑ e₂),
//! - direction n = cos θ e₃ + sin θ (cos φ e₁ + sin φ e₂),
//!
//! with φₑ = `exit_azimuth` (default φ). A ring of equal holes is a
//! `[[device.pattern]]`. An `[[injector]]` names a device and a hole; its
//! `position` and `direction` are then computed (or, when given as well,
//! checked against the hole). Nothing here changes an injector that names
//! no device.

use anyhow::{Result, bail};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// The surface of a device (D343).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Geometry {
    /// No surface: the frame and the holes only, drawn as lines.
    #[default]
    None,
    /// The standard gasoline tip of `assets/injectors/` (8 holes, 37°,
    /// 0.165 mm; ECN Spray G numbers, our drawing). Without holes of its
    /// own the device takes these.
    StandardGasoline,
    /// The standard diesel tip of `assets/injectors/` (7 holes, 77.5°,
    /// 0.11 mm on a sac cone; our drawing).
    StandardDiesel,
    /// The user's own STL, `file`.
    Stl,
}

/// The two standard tips: their STL (millimetres, the frame of D343) and
/// their holes. Keep in step with `tools/data/injector_standard.py`.
struct Standard {
    stl: &'static [u8],
    count: usize,
    polar: f64,
    /// Exit circle radius (mm).
    exit: f64,
    /// Exit z on the tip (mm, negative: behind the tip point).
    exit_z: f64,
    /// Hole diameter (mm).
    hole: f64,
}

fn standard(g: Geometry) -> Option<Standard> {
    match g {
        // dome of radius 4 mm through the tip: z(1 mm) = -(4 - √15)
        Geometry::StandardGasoline => Some(Standard {
            stl: include_bytes!("../../../assets/injectors/standard-gasoline.stl"),
            count: 8,
            polar: 37.0,
            exit: 1.0,
            exit_z: -(4.0 - 15f64.sqrt()),
            hole: 0.165,
        }),
        // sac dome of radius 1.2 mm: z(0.65 mm) = -(1.2 - √(1.44 - 0.4225))
        Geometry::StandardDiesel => Some(Standard {
            stl: include_bytes!("../../../assets/injectors/standard-diesel.stl"),
            count: 7,
            polar: 77.5,
            exit: 0.65,
            exit_z: -(1.2 - (1.44f64 - 0.4225).sqrt()),
            hole: 0.11,
        }),
        _ => None,
    }
}

/// One hole of a device, `[[device.hole]]` (D341).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Hole {
    /// Name of the hole (an injector names it with `hole`).
    pub name: String,
    /// Angle of the hole axis from the device axis (degrees, 0 to 180).
    #[serde(default)]
    pub polar: f64,
    /// Azimuth of the hole axis about the device axis (degrees, from e₁
    /// towards e₂).
    #[serde(default)]
    pub azimuth: f64,
    /// Distance of the hole exit centre from the device axis (m, default 0).
    #[serde(default)]
    pub radius: f64,
    /// Distance of the hole exit along the axis from `mount` into the
    /// chamber (m, default 0; negative: inside the body).
    #[serde(default)]
    pub setback: f64,
    /// Azimuth of the exit centre (degrees) when it lies on another plane
    /// than the hole axis; default: `azimuth`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exit_azimuth: Option<f64>,
    /// Diameter of the hole (m), for the picture only; the flow keys stay
    /// on the injector (`hole_diameter`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub diameter: Option<f64>,
}

/// A ring of equal holes, `[[device.pattern]]` (D341): `count` holes at
/// `first_azimuth` + k·360°/`count`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Pattern {
    /// Number of holes (1 or more).
    pub count: usize,
    /// Angle of every hole axis from the device axis (degrees).
    #[serde(default)]
    pub polar: f64,
    /// Distance of the exit centres from the axis (m, default 0).
    #[serde(default)]
    pub radius: f64,
    /// Distance of the exits along the axis from `mount` (m, default 0).
    #[serde(default)]
    pub setback: f64,
    /// Azimuth of the first hole (degrees, default 0).
    #[serde(default)]
    pub first_azimuth: f64,
    /// Names of the holes (default `1` to `count`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub names: Vec<String>,
    /// Diameter of the holes (m), for the picture only.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub diameter: Option<f64>,
}

/// A mounted injector or nozzle, `[[device]]` (D340).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Device {
    /// Name (an injector names it with `device`).
    pub name: String,
    /// The tip of the injector, the origin of the device frame (m).
    pub mount: [f64; 3],
    /// The injector axis, from the body into the chamber (any length, it
    /// is normalised).
    pub axis: [f64; 3],
    /// A direction not along `axis`: its part normal to the axis is the
    /// zero of the azimuth (e₁). Default: the global axis least aligned
    /// with `axis` (x for an injector along z), the rule of `sector`
    /// (D322).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference: Option<[f64; 3]>,
    /// Turns the whole device about its axis (degrees, default 0), so one
    /// description serves several mounting positions.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub clocking: f64,
    /// Holes one by one.
    #[serde(default, rename = "hole", skip_serializing_if = "Vec::is_empty")]
    pub holes: Vec<Hole>,
    /// Rings of equal holes.
    #[serde(default, rename = "pattern", skip_serializing_if = "Vec::is_empty")]
    pub patterns: Vec<Pattern>,
    /// The surface drawn in the output (D343): `none` (default),
    /// `standard-gasoline`, `standard-diesel` or `stl` with `file`. A
    /// standard tip without holes of its own brings its holes.
    #[serde(default, skip_serializing_if = "is_none_geometry")]
    pub geometry: Geometry,
    /// The STL of `geometry = "stl"` (ASCII or binary), relative to the file
    /// that names it. Drawn in its own frame: the tip at the origin, the
    /// injector axis along +z into the chamber, the azimuth zero along +x
    /// (or as `stl_origin`, `stl_axis`, `stl_reference` say).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file: Option<PathBuf>,
    /// Factor on the coordinates of `file` (default 1; 1e-3 for a file in
    /// millimetres).
    #[serde(default = "one", skip_serializing_if = "is_one")]
    pub scale: f64,
    /// The tip point in the coordinates of `file` (default the origin).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stl_origin: Option<[f64; 3]>,
    /// The injector axis in the coordinates of `file` (default +z).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stl_axis: Option<[f64; 3]>,
    /// The zero of the azimuth in the coordinates of `file` (default +x).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stl_reference: Option<[f64; 3]>,
    /// Length of the lines that show the hole axes in the output (m;
    /// default five hole diameters where the hole has one, else 1 mm).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub glyph_length: Option<f64>,
}

fn is_none_geometry(g: &Geometry) -> bool {
    *g == Geometry::None
}

fn one() -> f64 {
    1.0
}

fn is_one(x: &f64) -> bool {
    *x == 1.0
}

fn is_zero(x: &f64) -> bool {
    *x == 0.0
}

/// The right-handed frame of a device: origin, e₁, e₂, e₃ (the axis).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Frame {
    /// The mounting point (m).
    pub origin: [f64; 3],
    /// Zero of the azimuth.
    pub e1: [f64; 3],
    /// e₃ × e₁.
    pub e2: [f64; 3],
    /// The injector axis.
    pub e3: [f64; 3],
}

/// A hole in case coordinates.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlacedHole {
    /// Name of the hole.
    pub name: String,
    /// Centre of the exit (m).
    pub exit: [f64; 3],
    /// Unit direction of the hole axis.
    pub direction: [f64; 3],
    /// Polar angle (degrees).
    pub polar: f64,
    /// Azimuth (degrees).
    pub azimuth: f64,
    /// Diameter (m) when given.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub diameter: Option<f64>,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn finite(a: [f64; 3]) -> bool {
    a.iter().all(|x| x.is_finite())
}

fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = dot(a, a).sqrt();
    (n.is_finite() && n > 0.0).then(|| [a[0] / n, a[1] / n, a[2] / n])
}

fn combo(o: [f64; 3], terms: &[(f64, [f64; 3])]) -> [f64; 3] {
    let mut x = o;
    for (s, v) in terms {
        for k in 0..3 {
            x[k] += s * v[k];
        }
    }
    x
}

impl Device {
    /// The frame of the device; fails on zero or non-finite vectors and a
    /// reference along the axis.
    pub fn frame(&self) -> Result<Frame> {
        let n = &self.name;
        if !finite(self.mount) {
            bail!("[[device]] {n}: mount must be finite");
        }
        let Some(e3) = finite(self.axis).then(|| unit(self.axis)).flatten() else {
            bail!("[[device]] {n}: axis must be finite and not zero");
        };
        let r = match self.reference {
            Some(r) => {
                if !finite(r) {
                    bail!("[[device]] {n}: reference must be finite");
                }
                r
            }
            None => {
                let mut best = 0;
                for k in 1..3 {
                    if e3[k].abs() < e3[best].abs() {
                        best = k;
                    }
                }
                let mut t = [0.0; 3];
                t[best] = 1.0;
                t
            }
        };
        let rn = dot(r, r).sqrt();
        let along = dot(r, e3);
        let normal = combo(r, &[(-along, e3)]);
        // a reference within about 0.06° of the axis has no usable normal part
        let Some(e1) = unit(normal).filter(|_| dot(normal, normal).sqrt() > 1e-3 * rn) else {
            bail!("[[device]] {n}: reference must not lie along the axis");
        };
        if !self.clocking.is_finite() {
            bail!("[[device]] {n}: clocking must be finite (degrees)");
        }
        let (s, c) = self.clocking.to_radians().sin_cos();
        let e2 = cross(e3, e1);
        let e1 = combo([0.0; 3], &[(c, e1), (s, e2)]);
        let e2 = cross(e3, e1);
        Ok(Frame {
            origin: self.mount,
            e1,
            e2,
            e3,
        })
    }

    /// All holes of the device in case coordinates: the `hole` tables in
    /// order, then the patterns in order.
    pub fn placed_holes(&self) -> Result<Vec<PlacedHole>> {
        let f = self.frame()?;
        let n = &self.name;
        let mut list: Vec<Hole> = self.holes.clone();
        let std_pattern;
        let patterns = if self.holes.is_empty() && self.patterns.is_empty() {
            match standard(self.geometry) {
                Some(st) => {
                    std_pattern = [Pattern {
                        count: st.count,
                        polar: st.polar,
                        radius: st.exit * 1e-3,
                        setback: st.exit_z * 1e-3,
                        first_azimuth: 0.0,
                        names: Vec::new(),
                        diameter: Some(st.hole * 1e-3),
                    }];
                    &std_pattern[..]
                }
                None => &self.patterns[..],
            }
        } else {
            &self.patterns[..]
        };
        for (p_i, p) in patterns.iter().enumerate() {
            if p.count == 0 {
                bail!("[[device]] {n}: pattern {} has count 0", p_i + 1);
            }
            if !p.names.is_empty() && p.names.len() != p.count {
                bail!(
                    "[[device]] {n}: pattern {} gives {} names for {} holes",
                    p_i + 1,
                    p.names.len(),
                    p.count
                );
            }
            for k in 0..p.count {
                let name = p
                    .names
                    .get(k)
                    .cloned()
                    .unwrap_or_else(|| (k + 1).to_string());
                list.push(Hole {
                    name,
                    polar: p.polar,
                    azimuth: p.first_azimuth + 360.0 * k as f64 / p.count as f64,
                    radius: p.radius,
                    setback: p.setback,
                    exit_azimuth: None,
                    diameter: p.diameter,
                });
            }
        }
        if list.is_empty() {
            bail!("[[device]] {n}: give at least one [[device.hole]] or [[device.pattern]]");
        }
        let mut names = std::collections::BTreeSet::new();
        let mut out = Vec::with_capacity(list.len());
        for h in list {
            let hn = &h.name;
            if hn.is_empty() || hn == "all" {
                bail!("[[device]] {n}: a hole must have a name other than \"all\"");
            }
            if !names.insert(hn.clone()) {
                bail!("[[device]] {n}: hole {hn} is defined twice");
            }
            if !(h.polar.is_finite() && (0.0..=180.0).contains(&h.polar)) {
                bail!("[[device]] {n}: hole {hn}: polar must be in [0, 180] degrees");
            }
            if !(h.azimuth.is_finite() && h.exit_azimuth.is_none_or(f64::is_finite)) {
                bail!("[[device]] {n}: hole {hn}: azimuth must be finite (degrees)");
            }
            if !(h.radius.is_finite() && h.radius >= 0.0) {
                bail!("[[device]] {n}: hole {hn}: radius must be zero or positive (m)");
            }
            if !h.setback.is_finite() {
                bail!("[[device]] {n}: hole {hn}: setback must be finite (m)");
            }
            if let Some(d) = h.diameter
                && !(d.is_finite() && d > 0.0)
            {
                bail!("[[device]] {n}: hole {hn}: diameter must be positive (m)");
            }
            let (st, ct) = h.polar.to_radians().sin_cos();
            let (sp, cp) = h.azimuth.to_radians().sin_cos();
            let (se, ce) = h.exit_azimuth.unwrap_or(h.azimuth).to_radians().sin_cos();
            let exit = combo(
                f.origin,
                &[
                    (h.setback, f.e3),
                    (h.radius * ce, f.e1),
                    (h.radius * se, f.e2),
                ],
            );
            let direction = combo([0.0; 3], &[(ct, f.e3), (st * cp, f.e1), (st * sp, f.e2)]);
            out.push(PlacedHole {
                name: h.name,
                exit,
                direction,
                polar: h.polar,
                azimuth: h.azimuth,
                diameter: h.diameter,
            });
        }
        Ok(out)
    }
}

impl Device {
    /// The surface of the device placed in case coordinates (D343): the
    /// triangles of the standard tip or of `file`, moved from their own
    /// frame into the device frame. Empty with `geometry = "none"`.
    pub fn surface(&self) -> Result<Vec<[[f64; 3]; 3]>> {
        let n = &self.name;
        let f = self.frame()?;
        let (tris, scale) = match self.geometry {
            Geometry::None => {
                if self.file.is_some() {
                    bail!("[[device]] {n}: file needs geometry = \"stl\"");
                }
                return Ok(Vec::new());
            }
            Geometry::Stl => {
                let Some(p) = &self.file else {
                    bail!("[[device]] {n}: geometry = \"stl\" needs file");
                };
                let t = foehncraft_mesh::format::stl::read(p)
                    .map_err(|e| anyhow::anyhow!("[[device]] {n}: {e}"))?;
                (t, self.scale)
            }
            g => {
                if self.file.is_some() {
                    bail!("[[device]] {n}: file needs geometry = \"stl\" (a standard tip has its own)");
                }
                let st = standard(g).expect("a standard tip");
                let t = foehncraft_mesh::format::stl::parse(st.stl)
                    .map_err(|e| anyhow::anyhow!("[[device]] {n}: standard tip: {e}"))?;
                // the standard files are in millimetres
                (t, 1e-3 * self.scale)
            }
        };
        if !(scale.is_finite() && scale > 0.0) {
            bail!("[[device]] {n}: scale must be positive and finite");
        }
        // the frame of the file
        let o = self.stl_origin.unwrap_or([0.0; 3]);
        let a = self.stl_axis.unwrap_or([0.0, 0.0, 1.0]);
        let r = self.stl_reference.unwrap_or([1.0, 0.0, 0.0]);
        let Some(g3) = finite(a).then(|| unit(a)).flatten() else {
            bail!("[[device]] {n}: stl_axis must be finite and not zero");
        };
        let rn = dot(r, r).sqrt();
        let normal = combo(r, &[(-dot(r, g3), g3)]);
        let Some(g1) = (finite(o) && finite(r))
            .then(|| unit(normal).filter(|_| dot(normal, normal).sqrt() > 1e-3 * rn))
            .flatten()
        else {
            bail!("[[device]] {n}: stl_origin and stl_reference must be finite, stl_reference not along stl_axis");
        };
        let g2 = cross(g3, g1);
        let place = |p: [f64; 3]| {
            let d = combo(p, &[(-1.0, o)]);
            let (x1, x2, x3) = (dot(d, g1) * scale, dot(d, g2) * scale, dot(d, g3) * scale);
            combo(f.origin, &[(x1, f.e1), (x2, f.e2), (x3, f.e3)])
        };
        Ok(tris.iter().map(|t| t.map(place)).collect())
    }

    /// The length of the line that shows hole `h` in the output.
    pub fn glyph(&self, h: &PlacedHole) -> f64 {
        self.glyph_length
            .unwrap_or_else(|| h.diameter.map_or(1e-3, |d| 5.0 * d))
    }
}

/// How far an injector's own `position` and `direction` lie from its hole
/// (D342): Some(message) when the position is more than 5 % of the exit
/// radius plus 10 µm from the exit (plus `offset`), or the direction more
/// than 0.5° from the hole axis.
pub fn mismatch(
    hole: &PlacedHole,
    radius: f64,
    offset: [f64; 3],
    position: [f64; 3],
    direction: [f64; 3],
) -> Option<String> {
    let want = combo(hole.exit, &[(1.0, offset)]);
    let d = combo(position, &[(-1.0, want)]);
    let dist = dot(d, d).sqrt();
    let tol = 0.05 * radius + 1e-5;
    let angle = unit(direction)
        .map(|u| dot(u, hole.direction).clamp(-1.0, 1.0).acos().to_degrees())
        .unwrap_or(180.0);
    let mut m = Vec::new();
    if dist > tol {
        m.push(format!(
            "position {dist:.3e} m from the exit of hole {} (tolerance {tol:.1e} m)",
            hole.name
        ));
    }
    if angle > 0.5 {
        m.push(format!(
            "direction {angle:.2}° from the axis of hole {} (tolerance 0.5°)",
            hole.name
        ));
    }
    (!m.is_empty()).then(|| m.join(", "))
}

/// The radius of a hole's exit from the device axis, for the tolerance of
/// [`mismatch`].
pub fn exit_radius(f: &Frame, h: &PlacedHole) -> f64 {
    let d = combo(h.exit, &[(-1.0, f.origin)]);
    let a = dot(d, f.e3);
    let r = combo(d, &[(-a, f.e3)]);
    dot(r, r).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spray_g() -> Device {
        Device {
            name: "g".into(),
            mount: [0.0; 3],
            axis: [1.0, 0.0, 0.0],
            reference: Some([0.0, 1.0, 0.0]),
            clocking: 0.0,
            holes: vec![],
            patterns: vec![Pattern {
                count: 8,
                polar: 33.0,
                radius: 1e-3,
                setback: 1e-4,
                first_azimuth: -11.25,
                names: vec![],
                diameter: Some(1.65e-4),
            }],
            geometry: Geometry::None,
            file: None,
            scale: 1.0,
            stl_origin: None,
            stl_axis: None,
            stl_reference: None,
            glyph_length: None,
        }
    }

    /// VD.4: a placed standard tip has its tip point at `mount` and its +z
    /// along the axis; a file in another frame given by `stl_*` places to
    /// the same triangles.
    #[test]
    fn placement_of_the_surface() {
        let mut d = spray_g();
        d.mount = [0.01, -0.02, 0.03];
        d.axis = [0.2, 0.3, -0.9];
        d.reference = None;
        d.clocking = 25.0;
        d.geometry = Geometry::StandardGasoline;
        let f = d.frame().unwrap();
        let t = d.surface().unwrap();
        assert!(t.len() > 1000);
        // the tip: the point of the largest height along the axis
        let mut top = f64::NEG_INFINITY;
        let mut at = [0.0; 3];
        for p in t.iter().flatten() {
            let h = dot(combo(*p, &[(-1.0, d.mount)]), f.e3);
            if h > top {
                top = h;
                at = *p;
            }
        }
        assert!(top.abs() < 1e-12, "{top}");
        let off = combo(at, &[(-1.0, d.mount)]);
        assert!(dot(off, off).sqrt() < 1e-12);
        // the same file written in another frame: origin (5, 6, 7) mm,
        // axis -y, reference +z, scale 1
        let raw = foehncraft_mesh::format::stl::parse(standard(Geometry::StandardGasoline).unwrap().stl).unwrap();
        let moved: Vec<[[f64; 3]; 3]> = raw
            .iter()
            .map(|t| t.map(|p| [5.0 - p[1], 6.0 - p[2], 7.0 + p[0]]))
            .collect();
        let dir = std::env::temp_dir().join(format!("device_stl_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("moved.stl");
        let mut text = String::from("solid moved\n");
        for t in &moved {
            text.push_str("facet normal 0 0 0\nouter loop\n");
            for p in t {
                text.push_str(&format!("vertex {:e} {:e} {:e}\n", p[0], p[1], p[2]));
            }
            text.push_str("endloop\nendfacet\n");
        }
        text.push_str("endsolid moved\n");
        std::fs::write(&path, text).unwrap();
        let mut e = d.clone();
        e.geometry = Geometry::Stl;
        e.file = Some(path);
        e.scale = 1e-3;
        e.stl_origin = Some([5.0, 6.0, 7.0]);
        e.stl_axis = Some([0.0, -1.0, 0.0]);
        e.stl_reference = Some([0.0, 0.0, 1.0]);
        let u = e.surface().unwrap();
        assert_eq!(u.len(), t.len());
        for (a, b) in t.iter().zip(&u) {
            for k in 0..3 {
                for i in 0..3 {
                    assert!((a[k][i] - b[k][i]).abs() < 1e-12);
                }
            }
        }
    }

    /// A standard tip without holes brings its own: 8 holes 37° from the
    /// axis on the exit circle of 1 mm, the exits on the tip's dome.
    #[test]
    fn standard_tip_brings_its_holes() {
        let mut d = spray_g();
        d.patterns.clear();
        d.geometry = Geometry::StandardGasoline;
        let hs = d.placed_holes().unwrap();
        assert_eq!(hs.len(), 8);
        assert!((hs[0].polar - 37.0).abs() < 1e-12);
        assert!((hs[0].exit[0] + (4.0 - 15f64.sqrt()) * 1e-3).abs() < 1e-15);
        d.geometry = Geometry::StandardDiesel;
        assert_eq!(d.placed_holes().unwrap().len(), 7);
        d.geometry = Geometry::None;
        assert!(d.placed_holes().is_err());
    }

    #[test]
    fn frame_is_right_handed_and_orthonormal() {
        let mut d = spray_g();
        d.axis = [0.3, -0.2, 0.9];
        d.reference = None;
        d.clocking = 17.0;
        let f = d.frame().unwrap();
        for (a, b) in [(f.e1, f.e2), (f.e2, f.e3), (f.e1, f.e3)] {
            assert!(dot(a, b).abs() < 1e-15);
        }
        for e in [f.e1, f.e2, f.e3] {
            assert!((dot(e, e) - 1.0).abs() < 1e-15);
        }
        let c = cross(f.e1, f.e2);
        assert!((dot(c, f.e3) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn spray_g_hole_one_matches_the_hand_written_file() {
        // examples/cases/flow/spray_g_g1_injectors.toml (D309)
        let pos = [1.000000e-04, 9.81175461e-04, -1.93128751e-04];
        let dir = [0.83867057, 0.53417395, -0.10625380];
        let h = &spray_g().placed_holes().unwrap()[0];
        assert_eq!(h.name, "1");
        let d = combo(pos, &[(-1.0, h.exit)]);
        assert!(dot(d, d).sqrt() < 2.01e-6, "{d:?}");
        for k in 0..3 {
            assert!((h.direction[k] - dir[k]).abs() < 1e-8);
        }
        // the 2 µm into the sector, normal to the plume plane
        let nrm = [
            0.0,
            (11.25f64).to_radians().sin(),
            (11.25f64).to_radians().cos(),
        ];
        let x = combo(h.exit, &[(2e-6, nrm)]);
        for k in 0..3 {
            assert!((x[k] - pos[k]).abs() < 1e-11, "{k}: {} {}", x[k], pos[k]);
        }
    }

    #[test]
    fn pattern_turns_by_equal_steps() {
        let hs = spray_g().placed_holes().unwrap();
        assert_eq!(hs.len(), 8);
        for (k, h) in hs.iter().enumerate() {
            let a = (k as f64 * 45.0).to_radians();
            let (s, c) = a.sin_cos();
            // hole 1 turned by k·45° about x
            let d0 = hs[0].direction;
            let want = [d0[0], c * d0[1] - s * d0[2], s * d0[1] + c * d0[2]];
            for i in 0..3 {
                assert!((h.direction[i] - want[i]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn refusals() {
        let mut d = spray_g();
        d.reference = Some([2.0, 0.0, 0.0]);
        assert!(d.frame().is_err());
        let mut d = spray_g();
        d.axis = [0.0; 3];
        assert!(d.frame().is_err());
        let mut d = spray_g();
        d.patterns[0].count = 0;
        assert!(d.placed_holes().is_err());
        let mut d = spray_g();
        d.holes.push(Hole {
            name: "1".into(),
            polar: 0.0,
            azimuth: 0.0,
            radius: 0.0,
            setback: 0.0,
            exit_azimuth: None,
            diameter: None,
        });
        assert!(d.placed_holes().is_err());
    }
}
