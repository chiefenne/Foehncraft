//! Repository tasks of Foehncraft (D203), run with `cargo xtask <task>`.
//!
//! - `layering`: the crate layering. Every dependency between the crates
//!   of the repository must be in [`ALLOWED`] (lower layers never depend
//!   on higher ones), and the front ends (the command line, the Python
//!   module) reach the engine only through the tool registry
//!   `foehncraft-tools` (D196): their sources must not call the engine
//!   entry points of [`ENGINE_CALLS`].

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

/// Each crate with the repository crates it may depend on (normal and
/// build dependencies; dev dependencies are free).
const ALLOWED: &[(&str, &[&str])] = &[
    ("foehncraft-core", &[]),
    ("foehncraft-kernels", &["foehncraft-core"]),
    ("foehncraft-parallel", &["foehncraft-core"]),
    (
        "foehncraft-mesh",
        &[
            "foehncraft-core",
            "foehncraft-kernels",
            "foehncraft-parallel",
        ],
    ),
    (
        "foehncraft-radiation",
        &[
            "foehncraft-core",
            "foehncraft-kernels",
            "foehncraft-parallel",
            "foehncraft-mesh",
        ],
    ),
    (
        "foehncraft-fv",
        &[
            "foehncraft-core",
            "foehncraft-kernels",
            "foehncraft-parallel",
            "foehncraft-mesh",
        ],
    ),
    // The structural solver (L53): numerics on the core crate; the mesh
    // crate only for its Gmsh reader (second-order solids).
    (
        "foehncraft-structure",
        &["foehncraft-core", "foehncraft-mesh"],
    ),
    // The Lagrangian parcels (L12): the parcel kernel and the tracker on the
    // mesh locator; the gas enters through a trait, so no `foehncraft-fv`.
    (
        "foehncraft-particles",
        &["foehncraft-core", "foehncraft-kernels", "foehncraft-mesh"],
    ),
    (
        "foehncraft-flow",
        &[
            "foehncraft-core",
            "foehncraft-kernels",
            "foehncraft-parallel",
            "foehncraft-mesh",
            "foehncraft-fv",
            "foehncraft-particles",
        ],
    ),
    (
        "foehncraft-bench",
        &["foehncraft-core", "foehncraft-kernels", "foehncraft-mesh"],
    ),
    (
        "foehncraft-cuda",
        &[
            "foehncraft-core",
            "foehncraft-kernels",
            "foehncraft-fv",
            "foehncraft-flow",
            "foehncraft-bench",
            "foehncraft-radiation",
        ],
    ),
    (
        "foehncraft-tools",
        &[
            "foehncraft-core",
            "foehncraft-parallel",
            "foehncraft-mesh",
            "foehncraft-fv",
            "foehncraft-flow",
            "foehncraft-bench",
            "foehncraft-structure",
        ],
    ),
    // The front ends. `foehncraft-flow` for its types (diagnostics, build
    // information, interrupts), `foehncraft-cuda` to register the GPU
    // backends, `foehncraft-fv` only for the hypre library directory its
    // build script reports to the CLI's build script (no calls: the
    // crate paths of FORBIDDEN_PATHS are refused in the sources).
    (
        "foehncraft-cli",
        &[
            "foehncraft-tools",
            "foehncraft-flow",
            "foehncraft-cuda",
            "foehncraft-fv",
        ],
    ),
    ("foehncraft-py", &["foehncraft-tools", "foehncraft-flow"]),
    ("xtask", &[]),
];

/// The front ends and their sources.
const FRONT_ENDS: &[&str] = &["crates/foehncraft-cli/src", "crates/foehncraft-py/src"];

/// Crates whose items front ends never use directly (every function of
/// them is a tool).
const FORBIDDEN_PATHS: &[&str] = &[
    "foehncraft_fv::",
    "foehncraft_mesh::",
    "foehncraft_bench::",
    "foehncraft_parallel::",
    "foehncraft_kernels::",
    "foehncraft_radiation::",
    "foehncraft_structure::",
];

/// Engine entry points that front ends call only through a tool.
const ENGINE_CALLS: &[&str] = &[
    "engine::load(",
    "engine::run(",
    "engine::check_file(",
    "engine::check_case(",
    "engine::estimate(",
    "engine::read_result(",
    "engine::latest_restart(",
    "derived::from_restart(",
    "materials::describe(",
    "materials::describe_in(",
    "m::describe_in(",
    "expr::compile(",
    "schema::json_schema(",
    "case_format::upgrade_text(",
    "cf::upgrade_text(",
    "foehncraft_mesh::case::run(",
];

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .map_or_else(|| PathBuf::from("."), Path::to_path_buf)
}

/// (crate, dependency) pairs of the repository crates in `manifest`.
fn edges(manifest: &Path) -> Result<Vec<(String, String)>, String> {
    let out = Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()))
        .args([
            "metadata",
            "--format-version",
            "1",
            "--no-deps",
            "--offline",
        ])
        .arg("--manifest-path")
        .arg(manifest)
        .output()
        .map_err(|e| format!("cargo metadata: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "cargo metadata {}: {}",
            manifest.display(),
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    let m: serde_json::Value =
        serde_json::from_slice(&out.stdout).map_err(|e| format!("cargo metadata: {e}"))?;
    let mut v = Vec::new();
    for p in m["packages"].as_array().into_iter().flatten() {
        let name = p["name"].as_str().unwrap_or_default().to_string();
        for d in p["dependencies"].as_array().into_iter().flatten() {
            let dep = d["name"].as_str().unwrap_or_default();
            let local = d.get("path").is_some_and(|x| !x.is_null());
            if local && d["kind"].as_str() != Some("dev") {
                v.push((name.clone(), dep.to_string()));
            }
        }
    }
    Ok(v)
}

fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    let mut entries: Vec<PathBuf> = rd.filter_map(|e| e.ok().map(|e| e.path())).collect();
    entries.sort();
    for p in entries {
        if p.is_dir() {
            rust_files(&p, out);
        } else if p.extension().is_some_and(|e| e == "rs") {
            out.push(p);
        }
    }
}

/// The findings of the layering check (empty: it passes).
fn layering(root: &Path) -> Result<Vec<String>, String> {
    let mut found = Vec::new();
    let mut all = edges(&root.join("Cargo.toml"))?;
    all.extend(edges(&root.join("crates/foehncraft-py/Cargo.toml"))?);
    let seen: BTreeSet<(String, String)> = all.into_iter().collect();
    for (krate, dep) in &seen {
        match ALLOWED.iter().find(|(k, _)| k == krate) {
            None => found.push(format!(
                "{krate}: not in the layering table of xtask/src/main.rs"
            )),
            Some((_, ok)) if !ok.contains(&dep.as_str()) => found.push(format!(
                "{krate} depends on {dep}, which its layer does not allow"
            )),
            Some(_) => {}
        }
    }
    for fe in FRONT_ENDS {
        let mut files = Vec::new();
        rust_files(&root.join(fe), &mut files);
        for f in files {
            let text = std::fs::read_to_string(&f).map_err(|e| format!("{}: {e}", f.display()))?;
            for (n, line) in text.lines().enumerate() {
                let code = line.split("//").next().unwrap_or("");
                for p in FORBIDDEN_PATHS {
                    if code.contains(p) {
                        let rel = f.strip_prefix(root).unwrap_or(&f);
                        found.push(format!(
                            "{}:{}: the front end uses `{}` directly; go through foehncraft-tools (D196)",
                            rel.display(),
                            n + 1,
                            p.trim_end_matches(':')
                        ));
                    }
                }
                for call in ENGINE_CALLS {
                    if code.contains(call) {
                        let rel = f.strip_prefix(root).unwrap_or(&f);
                        found.push(format!(
                            "{}:{}: the front end calls `{}` directly; call its tool of foehncraft-tools (D196)",
                            rel.display(),
                            n + 1,
                            call.trim_end_matches('(')
                        ));
                    }
                }
            }
        }
    }
    Ok(found)
}

fn main() -> ExitCode {
    let task = std::env::args().nth(1).unwrap_or_default();
    match task.as_str() {
        "layering" => match layering(&root()) {
            Ok(f) if f.is_empty() => {
                println!(
                    "crate layering: {} crates, every dependency allowed; the front ends call the engine through foehncraft-tools",
                    ALLOWED.len()
                );
                ExitCode::SUCCESS
            }
            Ok(f) => {
                eprintln!("crate layering: {} finding(s)", f.len());
                for x in f {
                    eprintln!("  {x}");
                }
                ExitCode::FAILURE
            }
            Err(e) => {
                eprintln!("crate layering: {e}");
                ExitCode::FAILURE
            }
        },
        _ => {
            eprintln!("usage: cargo xtask layering");
            ExitCode::from(2)
        }
    }
}
