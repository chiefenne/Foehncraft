//! MPI step W7 (D93, D84): a run on 2 to 4 ranks (threads of this
//! process) writes the output files of the run on one rank byte for byte:
//! the residuals, the ledger, the probes, the log columns of the run
//! record, the fields and the restart file; the result summary differs
//! only in the rank count of its provenance. Part 1: steady laminar runs;
//! part 2: passive scalars, the energy equation, k-ω SST and k-ε, the
//! rotating frame and transient runs; then the rotating zones of `[[mrf]]`
//! (steady and the transient frozen rotor).

use foehncraft_flow::engine::{self, RunOptions};
use std::path::Path;

const CASE: &str = r#"
case_format = "0.7"
[mesh.generate]
type = "poly"
dimension = 2
upper = [3.0, 1.0, 0.1]
cells = [30, 10, 1]
perturb = 0.3
[mesh]
renumber = "rcm"
[fluid]
density = 1.0
viscosity = 0.01
[initial]
potential_flow = true
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 101325.0
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
velocity = [-1.0, 0.0, 0.0]
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
tolerance = 1e-30
max_iterations = 30
log_every = 5
[verification]
require_convergence = false
[solver.momentum]
limiter = "venkatakrishnan"
[balance]
every = 10
[[output.probes]]
name = "eddy"
point = [2.5, 0.8, 0.05]
[log]
columns = ["residual.momentum", "probe.eddy.velocity_x", "max.velocity_magnitude", "boundary.xmax.mass_flow"]
[run]
restart_every = 20
"#;

fn run(dir: &Path, ranks: usize) -> serde_json::Value {
    run_text(CASE, dir, ranks)
}

fn run_text(text: &str, dir: &Path, ranks: usize) -> serde_json::Value {
    let mut case = foehncraft_flow::case::FlowCase::from_toml(text, &std::env::temp_dir()).unwrap();
    case.run.rank_count = ranks;
    case.run.threads = Some(ranks);
    let out = engine::run(&case, "w7", dir, &RunOptions::default(), |_| {}).unwrap();
    assert_eq!(
        out.exit_code(),
        0,
        "{ranks} ranks: {:?}",
        out.summary.diagnostics
    );
    let mut v = serde_json::to_value(&out.summary).unwrap();
    assert_eq!(v["provenance"]["ranks"], ranks);
    v["provenance"] = serde_json::Value::Null;
    // the hashes of the files without timings or settings
    if let Some(files) = v["files"].as_array_mut() {
        files.retain(|f| {
            !matches!(
                f["name"].as_str().unwrap_or(""),
                "flow.log"
                    | "flow.md"
                    | "flow.json"
                    | "run.jsonl"
                    | "case.resolved.toml"
                    | "ranks.vtkhdf"
            )
        });
    }
    v
}

/// The rows of the run record that do not carry timings.
fn record(dir: &Path) -> Vec<String> {
    std::fs::read_to_string(dir.join("run.jsonl"))
        .unwrap()
        .lines()
        .filter_map(|l| {
            let v: serde_json::Value = serde_json::from_str(l).ok()?;
            let kind = v["type"].as_str()?.to_string();
            matches!(
                kind.as_str(),
                "iteration" | "step" | "balance" | "balance_run" | "table_row"
            )
            .then(|| l.to_string())
        })
        .collect()
}

#[test]
fn steady_run_on_ranks_writes_the_files_of_one_rank() {
    let base = std::env::temp_dir().join(format!("foehn-w7-run-{}", std::process::id()));
    let one = base.join("1");
    let reference = run(&one, 1);
    let files = [
        "residuals.csv",
        "balance.csv",
        "probes.csv",
        "boundaries.csv",
        #[cfg(feature = "hdf5")]
        "fields.vtkhdf",
        #[cfg(feature = "hdf5")]
        "restart.h5",
    ];
    for n in 2..=4 {
        let dir = base.join(n.to_string());
        let summary = run(&dir, n);
        for f in files {
            let (a, b) = (one.join(f), dir.join(f));
            if !a.exists() {
                assert!(!b.exists(), "{f} only on {n} ranks");
                continue;
            }
            assert!(
                std::fs::read(&a).unwrap() == std::fs::read(&b).unwrap(),
                "{f} on {n} ranks differs from one rank"
            );
        }
        assert_eq!(record(&one), record(&dir), "run.jsonl rows on {n} ranks");
        assert_eq!(
            summary["solves"], reference["solves"],
            "result.json solves on {n} ranks"
        );
        assert_eq!(
            summary["files"].as_array().map(Vec::len),
            reference["files"].as_array().map(Vec::len)
        );
    }
    let _ = std::fs::remove_dir_all(&base);
}

/// Runs `text` on one rank and on each of `ranks` ranks and checks that
/// every file of `files` (where one rank writes it), the rows of the run
/// record and the solves of the result summary are those of one rank.
fn same_as_one_rank(tag: &str, text: &str, files: &[&str], ranks: &[usize]) {
    let base = std::env::temp_dir().join(format!("foehn-w7-{tag}-{}", std::process::id()));
    let one = base.join("1");
    let reference = run_text(text, &one, 1);
    for &n in ranks {
        let dir = base.join(n.to_string());
        let summary = run_text(text, &dir, n);
        for f in files {
            let (a, b) = (one.join(f), dir.join(f));
            assert!(
                a.exists(),
                "{tag}: one rank wrote no {f}: {:?}",
                std::fs::read_dir(&one)
                    .map(|d| d
                        .filter_map(|e| e.ok())
                        .map(|e| e.file_name())
                        .collect::<Vec<_>>())
                    .unwrap_or_default()
            );
            assert!(
                std::fs::read(&a).unwrap() == std::fs::read(&b).unwrap(),
                "{tag}: {f} on {n} ranks differs from one rank"
            );
        }
        assert_eq!(
            record(&one),
            record(&dir),
            "{tag}: run.jsonl rows on {n} ranks"
        );
        assert_eq!(
            summary["solves"], reference["solves"],
            "{tag}: result.json solves on {n} ranks"
        );
        assert_eq!(
            summary["files"].as_array().map(Vec::len),
            reference["files"].as_array().map(Vec::len),
            "{tag}: the number of files on {n} ranks"
        );
    }
    let _ = std::fs::remove_dir_all(&base);
}

/// The channel of [`CASE`] with a passive scalar and the energy equation
/// (a heated lower wall), the turbulence model `model`, the `[output]`
/// keys `output` and the tables `extra`.
fn physics_case(model: &str, output: &str, extra: &str) -> String {
    let mut text = CASE
        .replace(
            "viscosity = 0.01\n",
            "viscosity = 0.01\nspecific_heat = 1000.0\nconductivity = 0.05\n",
        )
        .replace(
            "[[output.probes]]",
            &format!("[output]\n{output}\n[[output.probes]]"),
        );
    if model != "laminar" {
        text = text.replace(
            "potential_flow = true\n",
            "potential_flow = true\nturbulence_intensity = 0.05\nviscosity_ratio = 10.0\n",
        );
        text.push_str(&format!("[turbulence]\nmodel = \"{model}\"\n"));
    }
    format!(
        r#"{text}
[energy]
initial = 300.0
[boundary.xmin.energy]
temperature = 300.0
[boundary.ymin.energy]
temperature = 330.0
[scalars.c]
diffusivity = 0.01
initial = 0.0
[boundary.xmin.scalar.c]
value = 1.0
[boundary.ymin.scalar.c]
value = 0.0
{extra}
"#
    )
}

const STEADY_FILES: &[&str] = &[
    "residuals.csv",
    "balance.csv",
    "boundaries.csv",
    "wall_ymin.csv",
    "wall_ymax.csv",
    #[cfg(feature = "hdf5")]
    "fields.vtkhdf",
    #[cfg(feature = "hdf5")]
    "restart.h5",
];

const TRANSIENT_FILES: &[&str] = &[
    "history.csv",
    "heat.csv",
    "balance.csv",
    "probes.csv",
    "forces.csv",
    "boundaries.csv",
    #[cfg(feature = "hdf5")]
    "fields.vtkhdf",
    #[cfg(feature = "hdf5")]
    "fields_000004.vtkhdf",
    #[cfg(feature = "hdf5")]
    "restart.h5",
];

#[test]
fn steady_scalars_energy_and_sst_on_ranks_write_the_files_of_one_rank() {
    same_as_one_rank(
        "sst",
        &physics_case("k-omega-sst", "", ""),
        STEADY_FILES,
        &[2, 3],
    );
}

#[test]
fn steady_k_epsilon_in_a_rotating_frame_on_ranks_writes_the_files_of_one_rank() {
    same_as_one_rank(
        "keps-frame",
        &physics_case(
            "k-epsilon",
            "",
            "[flow]\nframe_rotation = [0.0, 0.0, 0.4]\n",
        ),
        STEADY_FILES,
        &[2, 4],
    );
}

#[test]
fn transient_run_on_ranks_writes_the_files_of_one_rank() {
    let text = physics_case(
        "laminar",
        "fields_every = 4",
        "[time]\nscheme = \"bdf2\"\ntime_step = 0.05\nsteps = 8\n",
    )
    .replace("restart_every = 20", "restart_every = 4");
    same_as_one_rank("transient", &text, TRANSIENT_FILES, &[2, 3]);
}

#[test]
fn adaptive_transient_sst_run_on_ranks_writes_the_files_of_one_rank() {
    let text = physics_case(
        "k-omega-sst",
        "fields_every = 4",
        "[time]\nscheme = \"bdf2\"\ntime_step = 0.02\nsteps = 6\nadaptive = true\n\
         max_cfl = 2.0\n",
    );
    same_as_one_rank("transient-sst", &text, TRANSIENT_FILES, &[2]);
}

#[test]
fn refuses_what_does_not_run_on_ranks_yet() {
    let mut case = foehncraft_flow::case::FlowCase::from_toml(
        &format!(
            "{CASE}\n[time]\nscheme = \"bdf2\"\ntime_step = 0.05\nsteps = 2\n\
             [output.statistics]\nfields = true\n"
        ),
        &std::env::temp_dir(),
    )
    .unwrap();
    case.run.rank_count = 2;
    let err = foehncraft_flow::ranks::supported_run(&case).unwrap_err();
    assert!(err.to_string().contains("statistics"), "{err}");
}

/// Laminar Taylor–Couette flow with the inner half of the gap a rotating
/// zone (`[[mrf]]`, D55): the frame flux on the faces between ranks needs
/// the zone of the halo cells, the torque and the shaft power of the rotor
/// are log columns. `extra` is appended (a `[time]` table, a turbulence
/// model).
fn mrf_case(extra: &str) -> String {
    format!(
        r#"
case_format = "0.7"
[mesh.annulus]
r_inner = 0.05
r_outer = 0.10
depth = 0.01
n_theta = 30
n_r = 8
rotor_layers = 4
[mesh]
renumber = "rcm"
[fluid]
density = 1.0
viscosity = 5e-4
specific_heat = 1000.0
conductivity = 0.05
[constants]
spin = 10.0
[flow]
convection = true
[initial]
pressure = 101325.0
[[mrf]]
zone = "rotor"
axis = [0.0, 0.0, 1.0]
omega = 10.0
[boundary.inner]
type = "wall"
velocity = ["-spin * y", "spin * x", 0.0]
[boundary.inner.energy]
temperature = 320.0
[boundary.outer]
type = "wall"
[boundary.outer.energy]
temperature = 300.0
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[energy]
initial = 300.0
[solver]
tolerance = 1e-30
max_iterations = 30
log_every = 5
[verification]
require_convergence = false
[balance]
every = 10
[[output.probes]]
name = "gap"
point = [0.0, 0.07, 0.005]
[log]
columns = ["residual.momentum", "probe.gap.velocity_y", "rotor.rotor.torque", "rotor.rotor.power"]
[run]
restart_every = 10
{extra}
"#
    )
}

#[test]
fn steady_mrf_zone_on_ranks_writes_the_files_of_one_rank() {
    same_as_one_rank(
        "mrf",
        &mrf_case(""),
        &[
            "residuals.csv",
            "balance.csv",
            "boundaries.csv",
            "wall_inner.csv",
            "wall_outer.csv",
            #[cfg(feature = "hdf5")]
            "fields.vtkhdf",
            #[cfg(feature = "hdf5")]
            "restart.h5",
        ],
        &[2, 3, 4],
    );
}

#[test]
fn transient_frozen_rotor_with_sst_on_ranks_writes_the_files_of_one_rank() {
    let text = mrf_case(
        "[time]\nscheme = \"bdf2\"\ntime_step = 0.02\nsteps = 6\n\
         [turbulence]\nmodel = \"k-omega-sst\"\n[output]\nfields_every = 3\n",
    )
    .replace(
        "pressure = 101325.0\n",
        "pressure = 101325.0\nturbulence_intensity = 0.05\nviscosity_ratio = 10.0\n",
    )
    .replace("restart_every = 10", "restart_every = 3");
    same_as_one_rank(
        "mrf-transient",
        &text,
        &[
            "history.csv",
            "balance.csv",
            "probes.csv",
            "boundaries.csv",
            #[cfg(feature = "hdf5")]
            "fields.vtkhdf",
            #[cfg(feature = "hdf5")]
            "fields_000003.vtkhdf",
            #[cfg(feature = "hdf5")]
            "restart.h5",
        ],
        &[2, 3],
    );
}

/// The Taylor–Couette gap with the inner half a sliding rotor in conformal
/// stepping (D182): 32 sectors, one sector per step, so a rotor cell lands
/// on a cell of another rank and the old levels of its time terms come
/// from there. A passive scalar is carried by the rotor; `extra` adds a
/// turbulence model.
fn sliding_case(extra: &str) -> String {
    format!(
        r#"
case_format = "0.7"
[mesh.annulus]
r_inner = 0.05
r_outer = 0.10
depth = 0.01
n_theta = 32
n_r = 8
rotor_layers = 4
[mesh]
renumber = "rcm"
[fluid]
density = 1.0
viscosity = 5e-4
[constants]
spin = 10.0
[flow]
convection = true
[time]
scheme = "bdf2"
time_step = 0.019634954084936207
steps = 8
start_time = 0.0
[[sliding]]
zone = "rotor"
axis = [0.0, 0.0, 1.0]
omega = 10.0
[boundary.inner]
type = "wall"
velocity = ["-spin * y", "spin * x", 0.0]
[boundary.inner.scalar.dye]
value = 1.0
[boundary.outer]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[scalars.dye]
diffusivity = 1e-4
initial = 0.0
[solver]
tolerance = 1e-10
max_iterations = 3000
[balance]
every = 4
[[output.probes]]
name = "gap"
point = [0.0, 0.07, 0.005]
[log]
columns = ["probe.gap.velocity_y", "rotor.rotor.angle", "rotor.rotor.torque"]
[output]
fields_every = 4
[run]
restart_every = 4
{extra}
"#
    )
}

const SLIDING_FILES: &[&str] = &[
    "history.csv",
    "balance.csv",
    "probes.csv",
    "boundaries.csv",
    #[cfg(feature = "hdf5")]
    "fields.vtkhdf",
    #[cfg(feature = "hdf5")]
    "fields_000004.vtkhdf",
    #[cfg(feature = "hdf5")]
    "restart.h5",
];

#[test]
fn sliding_mesh_in_conformal_stepping_on_ranks_writes_the_files_of_one_rank() {
    same_as_one_rank("sliding", &sliding_case(""), SLIDING_FILES, &[2, 3, 4]);
}

#[test]
fn turbulent_sliding_mesh_on_ranks_writes_the_files_of_one_rank() {
    let text = sliding_case("[turbulence]\nmodel = \"k-omega-sst\"\n").replace(
        "[flow]\n",
        "[initial]\nturbulence_intensity = 0.05\nviscosity_ratio = 10.0\n[flow]\n",
    );
    same_as_one_rank("sliding-sst", &text, SLIDING_FILES, &[3]);
}

/// The sliding mesh by imprinting (D182 variant B): a rotor of 32 sectors
/// and a stator of 48, the mesh built again at the rotor angle of every
/// step.
fn imprint_case(extra: &str) -> String {
    sliding_case(extra)
        .replace("n_theta = 32\n", "n_theta = 32\nn_theta_stator = 48\n")
        .replace(
            "[[sliding]]\nzone = \"rotor\"\n",
            "[[sliding]]\nzone = \"rotor\"\nrotor_interface = \"rotor_interface\"\n\
             stator_interface = \"stator_interface\"\n",
        )
        .replace(
            "[solver]\n",
            "[solver]\nextrapolation = false\nouter_iterations = 3\n",
        )
}

/// (a moving mesh writes no restart file yet)
const IMPRINT_FILES: &[&str] = &[
    "history.csv",
    "balance.csv",
    "probes.csv",
    "boundaries.csv",
    "forces.csv",
    "wall_inner.csv",
    #[cfg(feature = "hdf5")]
    "fields.vtkhdf",
    #[cfg(feature = "hdf5")]
    "fields_000004.vtkhdf",
];

#[test]
fn sliding_mesh_by_imprinting_on_ranks_writes_the_files_of_one_rank() {
    same_as_one_rank("imprint", &imprint_case(""), IMPRINT_FILES, &[2, 3]);
}

#[test]
fn turbulent_sliding_mesh_by_imprinting_on_ranks_writes_the_files_of_one_rank() {
    let text = imprint_case("[turbulence]\nmodel = \"k-omega-sst\"\n").replace(
        "[flow]\n",
        "[initial]\nturbulence_intensity = 0.05\nviscosity_ratio = 10.0\n[flow]\n",
    );
    same_as_one_rank("imprint-sst", &text, IMPRINT_FILES, &[4]);
}

#[test]
fn the_check_of_ranks_refuses_a_sliding_mesh() {
    let mut case =
        foehncraft_flow::case::FlowCase::from_toml(&imprint_case(""), &std::env::temp_dir())
            .unwrap();
    case.run.rank_count = 2;
    foehncraft_flow::ranks::supported_run(&case).unwrap();
    case.run.rank_count = 1;
    case.run.ranks = vec![2];
    let err = foehncraft_flow::ranks::supported(&case).unwrap_err();
    assert!(err.to_string().contains("sliding"), "{err}");
}
