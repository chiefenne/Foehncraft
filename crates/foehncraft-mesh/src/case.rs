//! The mesh case file and the `foehn mesh` run: build or read a mesh,
//! renumber it, check it, partition it, test point location, write the
//! report and the VTKHDF file.
//!
//! Case file (TOML, D39 style: every key optional, unknown keys are errors):
//!
//! ```toml
//! [mesh]
//! file = "channel.msh"     # a mesh file (Gmsh .msh 4.1, SU2 .su2; see
//!                          # crate::format), relative to the case file, or:
//! [mesh.generate]          # the built-in generator (crate::generate::BoxSpec)
//! type = "poly"
//! cells = [16, 16, 16]
//! # or [mesh.cylinder]: a channel with a cylinder (crate::cylinder::CylinderSpec)
//! # or [mesh.annulus]: the gap between two cylinders (crate::annulus::AnnulusSpec)
//! # or [mesh.revolved]: a block grid turned about the axis, rotor and stator
//! #   with an interface of revolution (crate::revolved::RevolvedSpec)
//! ```

use crate::Geometry;
use crate::annulus::AnnulusSpec;
use crate::cylinder::CylinderSpec;
use crate::decompose::{HaloPlan, decompose};
use crate::generate::BoxSpec;
use crate::locate::{Locator, TrackEnd};
use crate::mesh::Mesh;
use crate::partition::{Method, PartitionStats, partition};
use crate::quality::QualityReport;
use crate::renumber::{Ordering, bandwidth, renumber};
use crate::revolved::RevolvedSpec;
use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// A mesh case.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct MeshCase {
    /// Format of the case file (`case_format`, see [`crate::case_format`]).
    pub case_format: Option<String>,
    /// The mesh.
    pub mesh: MeshSection,
    /// Domain decomposition.
    pub partition: PartitionSection,
    /// Checks beyond the quality report.
    pub check: CheckSection,
    /// Output.
    pub output: OutputSection,
}

/// `[mesh]`
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct MeshSection {
    /// Mesh file (relative to the case file): Gmsh `.msh` 4.1, SU2 `.su2`, the
    /// Fluent face-based `.msh`, CGNS `.cgns` (HDF5), VTK `.vtk` / `.vtu` or
    /// an OpenFOAM case / `polyMesh` directory (docs/cases/meshes.md). The
    /// boundaries are the file's physical groups, markers, patches or zones;
    /// rename them with `boundaries`.
    pub file: Option<PathBuf>,
    /// Format of `file`: `gmsh`, `su2`, `fluent_msh`, `cgns`, `vtk` or
    /// `openfoam` (default: from the extension or directory, else from the
    /// content; `.msh` is Gmsh or the Fluent format, told apart by the
    /// content).
    pub format: Option<String>,
    /// Cell layer thickness in z when a 2D file (SU2 `NDIME= 2`) is extruded
    /// into one layer, in the units of the file (before `scale`); default the
    /// median edge length of the 2D cells.
    pub thickness: Option<f64>,
    /// Join coincident nodes of a multi-zone file (SU2 with `NZONE` > 1, CGNS
    /// with several zones). Off: the zones are kept apart, each a cell zone,
    /// and the interface is two boundaries, as `[[sliding]]` wants it. On:
    /// touching zones become one connected mesh.
    pub merge_nodes: bool,
    /// Factor applied to the coordinates of a file (0.001: mm to m).
    pub scale: f64,
    /// Built-in generator.
    pub generate: Option<BoxSpec>,
    /// Built-in channel with a cylinder.
    pub cylinder: Option<CylinderSpec>,
    /// Built-in annulus (the gap between two coaxial cylinders).
    pub annulus: Option<AnnulusSpec>,
    /// Built-in mesh of revolution: a block grid in (r, z) turned about the
    /// z axis, a rotor and a stator with an interface of cylinders, discs and
    /// cones (D186).
    pub revolved: Option<RevolvedSpec>,
    /// Convert a tetrahedral mesh into polyhedra (the dual mesh, one cell
    /// per node).
    pub polyhedral: bool,
    /// Cell order: `none`, `rcm` or `shuffle`.
    pub renumber: Ordering,
    /// Seed of `shuffle`.
    pub seed: u64,
    /// Boundary renaming, `old = "new"`.
    pub boundaries: BTreeMap<String, String>,
}

impl Default for MeshSection {
    fn default() -> Self {
        Self {
            file: None,
            format: None,
            thickness: None,
            merge_nodes: false,
            scale: 1.0,
            generate: None,
            cylinder: None,
            annulus: None,
            revolved: None,
            polyhedral: false,
            renumber: Ordering::None,
            seed: 2026,
            boundaries: BTreeMap::new(),
        }
    }
}

/// `[partition]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PartitionSection {
    /// Number of parts (1: no decomposition).
    pub parts: usize,
    /// `metis` or `rcb`.
    pub method: Method,
    /// Halo layers.
    pub halo_layers: usize,
}

impl Default for PartitionSection {
    fn default() -> Self {
        Self {
            parts: 1,
            method: Method::Metis,
            halo_layers: 1,
        }
    }
}

/// `[check]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CheckSection {
    /// Random points for the point-location test (0: no test).
    pub locate_points: usize,
    /// Seed of the points.
    pub seed: u64,
}

impl Default for CheckSection {
    fn default() -> Self {
        Self {
            locate_points: 1000,
            seed: 2026,
        }
    }
}

/// `[output]`
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct OutputSection {
    /// Output directory (default `output/<case name>`).
    pub directory: Option<PathBuf>,
    /// Write `mesh.vtkhdf`.
    pub vtkhdf: bool,
    /// Also write the mesh as `mesh.vtu` (VTK XML, polyhedra, with the names
    /// of the zones and boundaries; any tool reads it, and the format layer
    /// reads it back).
    pub vtu: bool,
}

impl Default for OutputSection {
    fn default() -> Self {
        Self {
            directory: None,
            vtkhdf: true,
            vtu: false,
        }
    }
}

/// A mesh file that cannot be read: missing, of an unknown format,
/// truncated or damaged (D202). It carries the message of the reader; the
/// flow layer gives it the diagnostic code `MESH003`.
#[derive(Debug, Clone)]
pub struct ReadError(pub String);

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ReadError {}

impl MeshCase {
    /// Reads a case file.
    pub fn from_file(path: &Path) -> Result<Self> {
        let text =
            std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let (mut case, _notes): (Self, _) =
            crate::case_format::read(&text).with_context(|| format!("{}", path.display()))?;
        case.case_format = Some(crate::case_format::CURRENT.into());
        if let Some(f) = &case.mesh.file
            && f.is_relative()
        {
            case.mesh.file = Some(path.parent().unwrap_or(Path::new(".")).join(f));
        }
        Ok(case)
    }

    /// A case that reads a mesh file (any format of [`crate::format`](mod@crate::format)) with
    /// the default settings.
    pub fn for_file(path: &Path) -> Self {
        let mut c = Self::default();
        c.mesh.file = Some(path.to_path_buf());
        c
    }
}

/// Result of the point-location test.
#[derive(Debug, Clone, Serialize, Default)]
pub struct LocateReport {
    /// Random points in the bounding box.
    pub points: usize,
    /// Points inside the mesh (found by the bin search).
    pub inside: usize,
    /// Walks that reached their point.
    pub walks_reached: usize,
    /// Walks blocked by the boundary of a non-convex domain (the bin search
    /// then finds the cell).
    pub walks_blocked: usize,
    /// Walks that ended in a cell that does not contain the point (must be
    /// zero).
    pub mismatches: usize,
    /// Mean number of faces crossed per walk.
    pub faces_per_walk: f64,
    /// Walks per second.
    pub walks_per_second: f64,
}

/// What a mesh run produced.
#[derive(Debug, Clone, Serialize)]
pub struct MeshReport {
    /// Case name.
    pub name: String,
    /// Where the mesh came from.
    pub source: String,
    /// Reader warnings.
    pub warnings: Vec<String>,
    /// Cell order and matrix bandwidth before and after renumbering.
    pub renumber: (String, usize, usize),
    /// Quality report.
    pub quality: QualityReport,
    /// Partition statistics.
    pub partition: Option<PartitionStats>,
    /// Halo cells and processor faces per rank.
    pub halo: Vec<(usize, usize, usize)>,
    /// Point-location test.
    pub locate: Option<LocateReport>,
    /// Seconds per step.
    pub timings: Vec<(String, f64)>,
}

/// Builds or reads the mesh of a case.
pub fn build_mesh(case: &MeshSection) -> Result<(Mesh, String, Vec<String>)> {
    let given = [
        case.file.is_some(),
        case.generate.is_some(),
        case.cylinder.is_some(),
        case.annulus.is_some(),
        case.revolved.is_some(),
    ];
    if given.iter().filter(|&&g| g).count() > 1 {
        bail!(
            "[mesh]: give one of `file`, `[mesh.generate]`, `[mesh.cylinder]`, `[mesh.annulus]` and \
             `[mesh.revolved]`"
        );
    }
    if let Some(spec) = &case.cylinder {
        let mut mesh = spec.generate().map_err(anyhow::Error::msg)?;
        for (old, new) in &case.boundaries {
            mesh.rename_boundary(old, new).map_err(anyhow::Error::msg)?;
        }
        return Ok((mesh, spec.describe(), Vec::new()));
    }
    if let Some(spec) = &case.annulus {
        let mut mesh = spec.generate().map_err(anyhow::Error::msg)?;
        for (old, new) in &case.boundaries {
            mesh.rename_boundary(old, new).map_err(anyhow::Error::msg)?;
        }
        return Ok((mesh, spec.describe(), Vec::new()));
    }
    if let Some(spec) = &case.revolved {
        let mut mesh = spec.generate().map_err(anyhow::Error::msg)?;
        for (old, new) in &case.boundaries {
            mesh.rename_boundary(old, new).map_err(anyhow::Error::msg)?;
        }
        return Ok((mesh, spec.describe(), Vec::new()));
    }
    let (mut mesh, source, warnings) = match (&case.file, &case.generate) {
        (Some(_), Some(_)) => bail!("[mesh]: give either `file` or `[mesh.generate]`, not both"),
        (None, None) => {
            bail!(
                "[mesh]: give `file`, a `[mesh.generate]`, a `[mesh.cylinder]`, a `[mesh.annulus]` or \
                 a `[mesh.revolved]` table"
            )
        }
        (Some(f), None) => {
            let options = crate::format::ReadOptions {
                thickness: case.thickness,
                merge_nodes: case.merge_nodes,
            };
            let g = crate::format::read_with(f, case.format.as_deref(), &options)
                .map_err(|e| anyhow::Error::new(ReadError(e)))?;
            let mut mesh = g.mesh;
            if case.scale != 1.0 {
                for p in &mut mesh.points {
                    *p = p.map(|x| x * case.scale);
                }
            }
            (mesh, format!("{} ({})", f.display(), g.format), g.warnings)
        }
        (None, Some(spec)) => {
            let mesh = spec.generate().map_err(anyhow::Error::msg)?;
            let src = format!(
                "generated: {} {}D, {} x {} x {} intervals, refine {}, perturb {}",
                spec.kind.name(),
                spec.dimension,
                spec.intervals()[0],
                spec.intervals()[1],
                spec.intervals()[2],
                spec.refine,
                spec.perturb
            );
            (mesh, src, Vec::new())
        }
    };
    let mut source = source;
    if case.polyhedral {
        let n = mesh.n_cells;
        mesh = crate::generate::polyhedral_from_tets(&mesh).map_err(anyhow::Error::msg)?;
        source += &format!(
            "; converted to polyhedra ({n} tetrahedra to {} cells)",
            mesh.n_cells
        );
    }
    for (old, new) in &case.boundaries {
        mesh.rename_boundary(old, new).map_err(anyhow::Error::msg)?;
    }
    Ok((mesh, source, warnings))
}

/// Runs a mesh case; `log` receives progress lines. Returns the report
/// and the files written.
pub fn run(
    case: &MeshCase,
    name: &str,
    log: &mut dyn FnMut(&str),
) -> Result<(MeshReport, Vec<PathBuf>)> {
    let mut timings = Vec::new();
    let t = Instant::now();
    let (mut mesh, source, warnings) = build_mesh(&case.mesh)?;
    timings.push(("build or read".to_string(), t.elapsed().as_secs_f64()));
    log(&format!("mesh: {source}"));
    for w in &warnings {
        log(&format!("warning: {w}"));
    }
    let bw0 = bandwidth(&mesh);
    let t = Instant::now();
    renumber(&mut mesh, case.mesh.renumber, case.mesh.seed).map_err(anyhow::Error::msg)?;
    timings.push(("renumber".to_string(), t.elapsed().as_secs_f64()));
    let bw1 = bandwidth(&mesh);
    let t = Instant::now();
    mesh.validate()
        .map_err(|e| anyhow::anyhow!("invalid mesh: {e}"))?;
    let geo = Geometry::compute(&mesh);
    let (quality, cq) = QualityReport::compute(&mesh, &geo);
    timings.push((
        "validate, geometry, quality".to_string(),
        t.elapsed().as_secs_f64(),
    ));
    log(&format!(
        "{} cells, {} faces, {} nodes; total volume {:.12e}; {}",
        quality.cells,
        quality.faces,
        quality.points,
        quality.total_volume,
        if quality.ok() { "valid" } else { "INVALID" }
    ));

    // partition
    let mut rank = vec![0.0; mesh.n_cells];
    let mut stats = None;
    let mut halo = Vec::new();
    if case.partition.parts > 1 {
        let t = Instant::now();
        let p = &case.partition;
        let part = partition(&mesh, &geo, p.parts, p.method).map_err(anyhow::Error::msg)?;
        let st = PartitionStats::compute(&mesh, &part, p.parts);
        let locals = decompose(&mesh, &part, p.parts, p.halo_layers).map_err(anyhow::Error::msg)?;
        // check the plans with a simulated exchange of the global numbers
        let plans: Vec<&HaloPlan> = locals.iter().map(|l| &l.halo).collect();
        let mut fields: Vec<Vec<f64>> = locals
            .iter()
            .map(|l| {
                let mut v = vec![-1.0; l.mesh.cell_id.len()];
                for (x, &id) in v.iter_mut().zip(&l.mesh.cell_id[..l.mesh.n_cells]) {
                    *x = id as f64;
                }
                v
            })
            .collect();
        HaloPlan::exchange_all(&plans, &mut fields);
        for (l, f) in locals.iter().zip(&fields) {
            l.mesh
                .validate()
                .map_err(|e| anyhow::anyhow!("rank {}: {e}", l.rank))?;
            anyhow::ensure!(
                f.iter().zip(&l.mesh.cell_id).all(|(&v, &g)| v == g as f64),
                "rank {}: halo exchange delivered wrong cells",
                l.rank
            );
            halo.push((l.rank, l.mesh.n_halo, l.processor_faces.iter().sum()));
        }
        for (c, &r) in part.iter().enumerate() {
            rank[c] = r as f64;
        }
        log(&format!(
            "partition: {} parts ({}), imbalance {:.3}, {} cut faces",
            p.parts,
            p.method.name(),
            st.imbalance,
            st.cut_faces
        ));
        stats = Some(st);
        timings.push(("partition and halo".to_string(), t.elapsed().as_secs_f64()));
    }

    // point location
    let mut locate = None;
    if case.check.locate_points > 0 {
        let t = Instant::now();
        let loc = Locator::new(&mesh, &geo);
        let (lo, hi) = mesh.bounding_box();
        let mut r = LocateReport {
            points: case.check.locate_points,
            ..Default::default()
        };
        let mut crossed = 0usize;
        let mut walk_time = 0.0;
        for i in 0..case.check.locate_points as u64 {
            let mut x = [0.0; 3];
            for d in 0..3 {
                x[d] = lo[d]
                    + (hi[d] - lo[d])
                        * foehncraft_core::rng::uniform_keyed([
                            foehncraft_core::rng::stream::MESH,
                            case.check.seed,
                            i,
                            d as u64,
                        ]);
            }
            if loc.locate_search(x).is_none() {
                continue;
            }
            r.inside += 1;
            let hint = (foehncraft_core::rng::keyed([
                foehncraft_core::rng::stream::MESH,
                case.check.seed,
                i,
                7,
            ]) % mesh.n_cells as u64) as usize;
            let tw = Instant::now();
            let tr = loc.track(loc.start(hint), x);
            walk_time += tw.elapsed().as_secs_f64();
            match tr.end {
                TrackEnd::Reached => {
                    r.walks_reached += 1;
                    crossed += tr.faces_crossed;
                    if !loc.contains(tr.position.cell, x) {
                        r.mismatches += 1;
                    }
                }
                _ => r.walks_blocked += 1,
            }
        }
        r.faces_per_walk = crossed as f64 / r.walks_reached.max(1) as f64;
        r.walks_per_second = (r.walks_reached + r.walks_blocked) as f64 / walk_time.max(1e-12);
        log(&format!(
            "point location: {} of {} points inside, {} walks reached, {} blocked, {} mismatches",
            r.inside, r.points, r.walks_reached, r.walks_blocked, r.mismatches
        ));
        anyhow::ensure!(r.mismatches == 0, "point location found the wrong cell");
        timings.push(("point location".to_string(), t.elapsed().as_secs_f64()));
        locate = Some(r);
    }

    let report = MeshReport {
        name: name.to_string(),
        source,
        warnings,
        renumber: (format!("{:?}", case.mesh.renumber).to_lowercase(), bw0, bw1),
        quality,
        partition: stats,
        halo,
        locate,
        timings,
    };

    // output
    let dir = case
        .output
        .directory
        .clone()
        .unwrap_or_else(|| PathBuf::from("output").join(name));
    std::fs::create_dir_all(&dir).with_context(|| format!("creating {}", dir.display()))?;
    let mut files = Vec::new();
    let md = dir.join("mesh_report.md");
    std::fs::write(&md, report.to_markdown())?;
    files.push(md);
    let js = dir.join("mesh_report.json");
    std::fs::write(&js, serde_json::to_string_pretty(&report)?)?;
    files.push(js);
    if case.output.vtu {
        let p = dir.join("mesh.vtu");
        crate::format::vtk::write_vtu(&mesh, &p, true).map_err(anyhow::Error::msg)?;
        files.push(p);
    }
    if case.output.vtkhdf {
        #[cfg(feature = "hdf5")]
        {
            use crate::vtkhdf::{CellField, write};
            let p = dir.join("mesh.vtkhdf");
            let mut fields = vec![
                CellField {
                    name: "volume",
                    values: &geo.cell_volume[..mesh.n_cells],
                },
                CellField {
                    name: "non_orthogonality",
                    values: &cq.non_orthogonality,
                },
                CellField {
                    name: "skewness",
                    values: &cq.skewness,
                },
                CellField {
                    name: "aspect_ratio",
                    values: &cq.aspect_ratio,
                },
                CellField {
                    name: "volume_ratio",
                    values: &cq.volume_ratio,
                },
                CellField {
                    name: "closure",
                    values: &cq.closure,
                },
            ];
            if case.partition.parts > 1 {
                fields.push(CellField {
                    name: "rank",
                    values: &rank,
                });
            }
            write(&p, &mesh, &fields)?;
            files.push(p);
        }
        #[cfg(not(feature = "hdf5"))]
        {
            let _ = (&cq, &rank);
            log("VTKHDF output skipped: built without the feature `hdf5`");
        }
    }
    Ok((report, files))
}

impl MeshReport {
    /// The report as Markdown.
    pub fn to_markdown(&self) -> String {
        let mut s = self.quality.to_markdown(&self.name);
        let _ = writeln!(s, "\n## Source\n\n{}\n", self.source);
        for w in &self.warnings {
            let _ = writeln!(s, "- warning: {w}");
        }
        let _ = writeln!(
            s,
            "\nCell order: {} (matrix bandwidth {} before, {} after)",
            self.renumber.0, self.renumber.1, self.renumber.2
        );
        if let Some(p) = &self.partition {
            let _ = writeln!(
                s,
                "\n## Partition\n\n{} parts, imbalance {:.3}, {} cut faces\n\n| rank | cells | halo cells | processor faces | neighbours |\n|---:|---:|---:|---:|---:|",
                p.parts, p.imbalance, p.cut_faces
            );
            for (r, h, f) in &self.halo {
                let _ = writeln!(
                    s,
                    "| {r} | {} | {h} | {f} | {} |",
                    p.cells[*r], p.neighbours[*r]
                );
            }
        }
        if let Some(l) = &self.locate {
            let _ = writeln!(
                s,
                "\n## Point location (D47)\n\n{} random points in the bounding box, {} inside the mesh; walks from a random cell: {} reached the point, {} blocked by the boundary (found by the bin search), **{} in a wrong cell**; {:.1} faces crossed per walk, {:.0} walks per second",
                l.points,
                l.inside,
                l.walks_reached,
                l.walks_blocked,
                l.mismatches,
                l.faces_per_walk,
                l.walks_per_second
            );
        }
        let _ = writeln!(s, "\n## Timings\n\n| step | seconds |\n|---|---:|");
        for (n, t) in &self.timings {
            let _ = writeln!(s, "| {n} | {t:.3} |");
        }
        s
    }
}
