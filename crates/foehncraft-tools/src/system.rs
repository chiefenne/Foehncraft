//! `system.*`: build, platform, backends and the catalogue itself.

use crate::{Action, Catalogue, Category, Ctx, NoInput, Tool};
use anyhow::Result;
use schemars::JsonSchema;
use serde::Serialize;

/// The build and the platform (`foehn info`).
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct InfoOutput {
    /// `foehncraft.info`.
    pub format: &'static str,
    /// Format version of the engine.
    pub format_version: &'static str,
    /// Version of Foehncraft.
    pub version: &'static str,
    /// Git commit of the build.
    pub commit: String,
    /// Operating system and architecture.
    pub platform: String,
    /// CPU model.
    pub cpu: String,
    /// CPU threads of the thread pool.
    pub threads: usize,
    /// MPI ranks of this process.
    pub ranks: usize,
    /// Cargo features of the build.
    pub features: Vec<String>,
    /// Kernel backends (benchmarks, scalar and Poisson studies).
    pub backends: Vec<String>,
    /// Backends of the flow loop.
    pub flow_backends: Vec<String>,
}

/// `system.info`.
pub struct Info;

impl Tool for Info {
    const ID: &'static str = "system.info";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str =
        "Version, commit, platform, threads, features and backends of this build.";
    const CATEGORY: Category = Category::System;
    const ACTION: Action = Action::Read;
    type Input = NoInput;
    type Output = InfoOutput;

    fn call(_: NoInput, ctx: &mut Ctx<'_>) -> Result<InfoOutput> {
        use foehncraft_parallel::Communicator as _;
        let b = ctx.build.clone();
        let threads = rayon::current_num_threads();
        let comm = foehncraft_parallel::SerialComm;
        let backends: Vec<String> = foehncraft_bench::backend::available()
            .into_iter()
            .map(Into::into)
            .collect();
        ctx.line(&format!("Foehncraft {}", foehncraft_core::VERSION));
        ctx.line(&format!("commit       {}", b.commit));
        ctx.line(&format!("platform     {}", b.platform));
        ctx.line(&format!("cpu          {}", b.cpu));
        ctx.line(&format!("threads      {threads}"));
        ctx.line("precision    f64 (reference path, D7)");
        ctx.line(&format!(
            "ranks        {} ({} communicator; MPI runs start with mpirun)",
            comm.size(),
            comm.name()
        ));
        ctx.line(&format!(
            "features     {}",
            if b.features.is_empty() {
                "none".to_string()
            } else {
                b.features.join(", ")
            }
        ));
        ctx.line(&format!("backends     {}", backends.join(", ")));
        Ok(InfoOutput {
            format: "foehncraft.info",
            format_version: foehncraft_flow::engine::FORMAT_VERSION,
            version: foehncraft_core::VERSION,
            commit: b.commit,
            platform: b.platform,
            cpu: b.cpu,
            threads,
            ranks: comm.size(),
            features: b.features,
            backends,
            flow_backends: foehncraft_flow::device::available()
                .into_iter()
                .map(Into::into)
                .collect(),
        })
    }
}

/// The backends of this build.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct BackendsOutput {
    /// Kernel backends compiled into this build.
    pub backends: Vec<String>,
}

/// `system.backends`.
pub struct Backends;

impl Tool for Backends {
    const ID: &'static str = "system.backends";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "The kernel backends compiled into this build.";
    const CATEGORY: Category = Category::System;
    const ACTION: Action = Action::Read;
    type Input = NoInput;
    type Output = BackendsOutput;

    fn call(_: NoInput, ctx: &mut Ctx<'_>) -> Result<BackendsOutput> {
        let backends: Vec<String> = foehncraft_bench::backend::available()
            .into_iter()
            .map(Into::into)
            .collect();
        for b in &backends {
            ctx.line(b);
        }
        Ok(BackendsOutput { backends })
    }
}

/// `system.commands`.
pub struct Commands;

impl Tool for Commands {
    const ID: &'static str = "system.commands";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str = "The catalogue of every tool with its input and output schemas.";
    const CATEGORY: Category = Category::System;
    const ACTION: Action = Action::Read;
    type Input = NoInput;
    type Output = Catalogue;

    fn call(_: NoInput, _: &mut Ctx<'_>) -> Result<Catalogue> {
        Ok(crate::catalogue())
    }
}
