//! `mesh.process`: build or read a mesh, check, partition and write it.

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Input of `mesh.process`, the options of `foehn mesh`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MeshInput {
    /// Mesh case file (TOML) or a mesh file (Gmsh or Fluent `.msh`, SU2,
    /// CGNS, VTK, an OpenFOAM case or polyMesh directory).
    pub input: PathBuf,
    /// Number of parts (overrides `partition.parts`).
    #[serde(default)]
    pub parts: Option<usize>,
    /// Partitioning method: `metis` or `rcb` (overrides `partition.method`).
    #[serde(default)]
    pub method: Option<String>,
    /// Cell order: `none`, `rcm` or `shuffle` (overrides `mesh.renumber`).
    #[serde(default)]
    pub renumber: Option<String>,
    /// Random points of the point-location test (overrides
    /// `check.locate_points`).
    #[serde(default)]
    pub locate: Option<usize>,
    /// Do not write the VTKHDF file.
    #[serde(default)]
    pub no_vtkhdf: bool,
    /// Also write the mesh as `mesh.vtu`.
    #[serde(default)]
    pub vtu: bool,
    /// Output directory (overrides `output.directory`).
    #[serde(default)]
    pub out: Option<PathBuf>,
}

/// Output of `mesh.process`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct MeshOutput {
    /// Every check of the quality report passed.
    pub valid: bool,
    /// The files written.
    pub files: Vec<PathBuf>,
}

fn bad(msg: &str) -> anyhow::Error {
    anyhow::anyhow!("{msg}")
}

/// `mesh.process`.
pub struct Process;

impl Tool for Process {
    const ID: &'static str = "mesh.process";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Build or read a mesh, check it, partition it and write the quality report and a VTKHDF file.";
    const CATEGORY: Category = Category::Mesh;
    const ACTION: Action = Action::WriteFiles;
    type Input = MeshInput;
    type Output = MeshOutput;

    /// The output directory, when the input names one.
    fn effects(a: &MeshInput, _: &mut Ctx<'_>) -> Result<crate::Effects> {
        Ok(crate::Effects {
            writes: a.out.iter().cloned().collect(),
            ..Default::default()
        })
    }

    fn call(a: MeshInput, ctx: &mut Ctx<'_>) -> Result<MeshOutput> {
        use foehncraft_mesh::case::MeshCase;
        ctx.roots.check_read(&a.input, "the mesh input")?;
        let mut case = if foehncraft_mesh::format::is_mesh_file(&a.input) {
            MeshCase::for_file(&a.input)
        } else {
            MeshCase::from_file(&a.input).map_err(foehncraft_flow::diag::case_error)?
        };
        let name = a
            .input
            .file_stem()
            .map_or("mesh".into(), |s| s.to_string_lossy().into_owned());
        if let Some(p) = a.parts {
            case.partition.parts = p;
        }
        if let Some(m) = &a.method {
            case.partition.method = match m.as_str() {
                "metis" => foehncraft_mesh::partition::Method::Metis,
                "rcb" => foehncraft_mesh::partition::Method::Rcb,
                _ => return Err(bad("--method: metis or rcb")),
            };
        }
        if let Some(r) = &a.renumber {
            use foehncraft_mesh::renumber::Ordering;
            case.mesh.renumber = match r.as_str() {
                "none" => Ordering::None,
                "rcm" => Ordering::Rcm,
                "shuffle" => Ordering::Shuffle,
                _ => return Err(bad("--renumber: none, rcm or shuffle")),
            };
        }
        if let Some(n) = a.locate {
            case.check.locate_points = n;
        }
        if a.no_vtkhdf {
            case.output.vtkhdf = false;
        }
        if a.vtu {
            case.output.vtu = true;
        }
        if a.out.is_some() {
            case.output.directory = a.out.clone();
        }
        if let Some(d) = &case.output.directory {
            ctx.roots.check_write(d)?;
        }
        let (report, files) = foehncraft_mesh::case::run(&case, &name, &mut |l| ctx.line(l))?;
        for f in &files {
            ctx.line(&format!("wrote {}", f.display()));
        }
        Ok(MeshOutput {
            valid: report.quality.ok(),
            files,
        })
    }

    fn exit_code(o: &MeshOutput) -> i32 {
        if o.valid { 0 } else { 2 }
    }
}

/// Input of `mesh.surface`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SurfaceInput {
    /// Triangle surface (STL, ASCII or binary).
    pub input: PathBuf,
    /// Factor from the file's length unit to metres (1e-3 for a file in
    /// millimetres).
    #[serde(default = "one")]
    pub scale: f64,
}

fn one() -> f64 {
    1.0
}

/// Output of `mesh.surface`.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct SurfaceOutput {
    /// Closed, outward, no self-intersections, no degenerate or duplicate
    /// triangles: usable as a solid body as it is.
    pub clean: bool,
    /// The report (counts, first defects, findings IMM001 to IMM009).
    pub report: foehncraft_mesh::surface::SurfaceReport,
}

/// `mesh.surface`: the surface check of L13 stage 0 (D225).
pub struct SurfaceCheck;

impl Tool for SurfaceCheck {
    const ID: &'static str = "mesh.surface";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "Check a triangle surface (STL) before it becomes an immersed body: closed, orientation, holes, non-manifold edges, self-intersections, degenerate and duplicate triangles, units (exact predicates, no repair).";
    const CATEGORY: Category = Category::Mesh;
    const ACTION: Action = Action::Read;
    type Input = SurfaceInput;
    type Output = SurfaceOutput;

    fn call(a: SurfaceInput, ctx: &mut Ctx<'_>) -> Result<SurfaceOutput> {
        ctx.roots.check_read(&a.input, "the surface")?;
        if !(a.scale.is_finite() && a.scale > 0.0) {
            return Err(bad("scale: a positive number"));
        }
        let s =
            foehncraft_mesh::surface::Surface::read_stl(&a.input, a.scale).map_err(|e| bad(&e))?;
        let report = foehncraft_mesh::surface::check(&s);
        ctx.line(&format!(
            "{} triangles, {} corners, {} pieces, area {:.6e} m², volume {:.6e} m³, {}",
            report.triangles,
            report.vertices,
            report.components,
            report.area,
            report.volume,
            if report.closed { "closed" } else { "open" }
        ));
        for f in &report.findings {
            ctx.line(&format!("{} {}", f.code, f.message));
        }
        Ok(SurfaceOutput {
            clean: report.clean(),
            report,
        })
    }

    fn exit_code(o: &SurfaceOutput) -> i32 {
        if o.report.findings.iter().any(|f| f.code == "IMM001") {
            2
        } else {
            0
        }
    }
}
