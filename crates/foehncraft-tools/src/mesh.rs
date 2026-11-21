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
    /// Mesh case file (TOML) or a mesh file (Gmsh or industry `.msh`, SU2,
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
