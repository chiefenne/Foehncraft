//! MPI step W7 (D93, D84): a run on 2 to 4 ranks (threads of this
//! process) writes the output files of the run on one rank byte for byte:
//! the residuals, the ledger, the probes, the log columns of the run
//! record, the fields and the restart file; the result summary differs
//! only in the rank count of its provenance. Part 1: steady laminar runs;
//! part 2: passive scalars, the energy equation, k-ω SST and k-ε, the
//! rotating frame and transient runs.

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
    // a transient run on a variable density books no storage in the mass
    // check of its verification (RUN007, a gap of stage 7c): its files
    // still compare
    let checked = text.contains("[time]") && text.contains("density = \"variable\"");
    assert!(
        out.exit_code() == 0
            || checked
                && out
                    .summary
                    .diagnostics
                    .iter()
                    .all(|d| d.severity != foehncraft_flow::diag::Severity::Error
                        || d.code == "RUN007"),
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
    let mut case = foehncraft_flow::case::FlowCase::from_toml(CASE, &std::env::temp_dir()).unwrap();
    case.run.rank_count = 2;
    case.run.ranks = vec![2];
    let err = foehncraft_flow::ranks::supported_run(&case).unwrap_err();
    assert!(
        err.to_string().contains("the check of [run] ranks"),
        "{err}"
    );
}

/// The channel of [`CASE`] with a three-species mixture (CH₄ in air at the
/// inlet, air in the channel; N₂ the closure species) and the tables
/// `extra`. `variable`: an ideal gas whose density, c_p, λ and μ follow the
/// mixture and the temperature (stage 7c, SF1), with the energy equation
/// and a heated lower wall; otherwise the fluid's constant properties.
fn species_case(variable: bool, extra: &str) -> String {
    let mut text = CASE.replace(
        "[[output.probes]]",
        "[output]\nwall_csv = false\n[[output.probes]]",
    );
    let energy = if variable {
        text = text
            .replace(
                "[fluid]\ndensity = 1.0\nviscosity = 0.01\n",
                "[materials.gas]\nkind = \"gas\"\nmolar_mass = 0.0289647\n\
                 density = { form = \"ideal-gas\", source = \"ideal gas\" }\n\
                 specific_heat = 1005.0\nconductivity = 0.0263\nviscosity = 1.85e-5\n\
                 [fluid]\nmaterial = \"gas\"\ntemperature = 300.0\n\
                 operating_pressure = 101325.0\npressure_values = \"gauge\"\n",
            )
            .replace("outlet\"\npressure = 101325.0", "outlet\"\npressure = 0.0")
            .replace("velocity = [1.0, 0.0, 0.0]", "velocity = [0.02, 0.0, 0.0]")
            .replace(
                "velocity = [-1.0, 0.0, 0.0]",
                "velocity = [-0.02, 0.0, 0.0]",
            );
        "[energy]\ninitial = 300.0\ndensity = \"variable\"\n\
         [boundary.xmin.energy]\ntemperature = 300.0\n\
         [boundary.ymin.energy]\ntemperature = 400.0\n\
         [boundary.ymax.energy]\nheat_flux = 0.0\n"
    } else {
        ""
    };
    let properties = if variable {
        ""
    } else {
        "properties = \"fluid\"\n"
    };
    format!(
        r#"{text}
{energy}
[species]
names = ["CH4", "O2", "N2"]
initial = {{ mole_fractions = {{ O2 = 0.21, N2 = 0.79 }} }}
{properties}schmidt = 0.7
[boundary.xmin.species]
mole_fractions = {{ CH4 = 0.3, O2 = 0.147, N2 = 0.553 }}
[solver.scalar]
limiter = "mlp"
{extra}
"#
    )
}

#[test]
fn steady_species_on_ranks_write_the_files_of_one_rank() {
    same_as_one_rank(
        "species",
        &species_case(false, ""),
        &[
            "residuals.csv",
            "balance.csv",
            "boundaries.csv",
            #[cfg(feature = "hdf5")]
            "fields.vtkhdf",
            #[cfg(feature = "hdf5")]
            "restart.h5",
        ],
        &[2, 3],
    );
}

#[test]
fn steady_species_on_a_variable_density_on_ranks_write_the_files_of_one_rank() {
    same_as_one_rank(
        "species-rho",
        &species_case(true, ""),
        &[
            "residuals.csv",
            "balance.csv",
            "boundaries.csv",
            #[cfg(feature = "hdf5")]
            "fields.vtkhdf",
        ],
        &[2, 4],
    );
}

const SPECIES_TRANSIENT_FILES: &[&str] = &[
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
];

/// The transient form of [`species_case`]: BDF2, 6 steps, fields every 3.
fn species_transient(variable: bool, extra: &str) -> String {
    species_case(
        variable,
        &format!("[time]\nscheme = \"bdf2\"\ntime_step = 0.5\nsteps = 6\n{extra}"),
    )
    .replace("[output]\n", "[output]\nfields_every = 3\n")
    .replace("restart_every = 20", "restart_every = 3")
}

#[test]
fn transient_species_with_statistics_on_ranks_write_the_files_of_one_rank() {
    let text = species_transient(
        false,
        "[output.statistics]\nstart = 1.0\nbatches = 2\nfields = true\n",
    );
    let mut files = SPECIES_TRANSIENT_FILES.to_vec();
    files.push("statistics.csv");
    #[cfg(feature = "hdf5")]
    files.push("statistics.vtkhdf");
    same_as_one_rank("species-statistics", &text, &files, &[2, 3]);
}

#[test]
fn transient_species_on_a_variable_density_on_ranks_write_the_files_of_one_rank() {
    same_as_one_rank(
        "species-rho-transient",
        &species_transient(true, ""),
        SPECIES_TRANSIENT_FILES,
        &[2, 3],
    );
}
