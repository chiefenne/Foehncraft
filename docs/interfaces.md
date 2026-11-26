# Interfaces

What scripts, studies and agents can rely on when they drive `foehn`: the
command line, the exit codes and signals, the files a run writes, the
diagnostic codes and the engine API behind the command line. Everything
here is versioned: a file or a JSON object carries `format` and
`format_version`, and a change that breaks a reader raises the major
version. The case-file keys are in [keys.md](keys.md) (generated from the
code by `foehn schema --markdown`).

## Command line (D62, D80)

| Option | Effect |
|---|---|
| `--json` | One JSON object on standard output at the end (the result summary, the check report, the command's report or an error object). Log lines go to standard error. |
| `--quiet` | No log lines on the terminal (the files are written as usual). |
| `--set KEY=VALUE` | Changes one case key before the run, repeatable. The value is read as TOML (`--set fluid.viscosity=1e-3`, `--set 'time.adaptive=false'`); if it is not valid TOML it is taken as a string (`--set run.backend=cuda`). Applied in order, after the case file and before the validation. The overrides are recorded in the result and in the run record and survive a re-read of the case file (D40). |
| `--restart latest` | Continues from the newest restart file in the output directory: the largest step or iteration in the file wins, `restart.h5` wins a tie. No file is `RUN011`. |
| `--read-root DIR`, `--write-root DIR` | Confines the files of a run (D82): every file read (case, mesh, restart) must lie below a read root and the output directory below a write root. A path outside is `CASE007`. Without roots nothing is confined. |
| `--check` | Validates the case, checks the boundary conditions and estimates the run (D81) without running it. Exit code 0 when the case can run, 1 when it has an error. |
| `foehn schema` | The JSON Schema of the case file (draft 2020-12, with `x-unit` and `x-runtime`). `--markdown` prints the key reference instead. |
| `foehn materials [NAME]` | The built-in material database (D94): every material with its density and viscosity at `--temperature` (K, default 293.15) and `--pressure` (Pa, default 101 325), or one material with its forms, data range and sources. With `--json` one object (`foehncraft.materials`). |
| `foehn case resolve CASE` | The case the way a run reads it: defaults filled in, `--set` overrides applied, as TOML (the `case.resolved.toml` of a run); with `--json` one object (`foehncraft.resolved`). |
| `foehn diag codes` | The diagnostic codes below with their severity and meaning; with `--json` one object (`foehncraft.codes`). |
| `foehn expression SOURCE` | Compiles an expression (D53) with `[constants]` and `[tables]` of `--env FILE` and evaluates it at each `--point x,y,z,t`; `--allowed` names the variables it may use (default `xyzt`). |
| `foehn commands` | The catalogue of the tool registry (D196) as JSON (`foehn.tool-catalogue`); `--markdown` prints [tools.md](tools.md). |
| `foehn call ID` | Runs a tool from JSON (`--input FILE`, `-` for standard input, default no fields) and prints the result envelope; see [Tool registry](#tool-registry-d196). |
| `foehn run --script FILE` | Replays the tool calls of a journal or a run record in order (`--keep-going` past a call that ends differently from its record); exit 1 when one does; see [Journal and replay](#journal-and-replay-d197). |
| `foehn mcp --read-root DIR --write-root DIR` | Serves every tool to an assistant as an MCP server on standard input and output (L18); refuses to start without both roots; see [MCP server](#mcp-server-l18) and [agents.md](agents.md). |
| `foehn post derive CASE` | Writes the fields with derived fields (D72) from a restart file of the case (`--restart`, default `latest`), the file the run writes with `[output] derived` at that iteration or step (byte for byte). `--fields a,b` picks the fields (default `output.derived`, or all), `--file` the output (default `derived.vtkhdf`, `derived_<N>.vtkhdf` for `restart_<N>.h5`, next to the restart file). `--set`, `--out`, `--read-root` and `--write-root` as for `flow`. |

## Exit codes (D64)

| Code | Status | Meaning |
|---|---|---|
| 0 | `finished` | The run ended as planned and every check of the case passed. |
| 1 | `error` | The case could not be read or run, or the solver failed: diverged (`RUN014`, with the divergence report below), a file could not be written. |
| 2 | `failed-checks` | The run ended, but a check failed: not converged (`RUN006`), an order, error, reference or benchmark band, a mass balance (`RUN007`). |
| 3 | `stopped` | Stopped on request before its end: the `.stop` file, a signal, a limit of `[run.limits]`, or an exceeded balance tolerance with `[balance] action = "stop"` (`BAL001`). The restart file and the output of the last iteration or step are written; `--restart latest` continues bitwise as if the run had not stopped. |

## Signals (D64, unix)

| Signal | Effect |
|---|---|
| `SIGUSR1` | Writes `restart.h5` after the current iteration or step; the run goes on. |
| `SIGTERM`, `SIGINT`, `SIGUSR2` | Stops after the current iteration or step with the restart file and the final output (exit 3, `RUN002`). |
| second `SIGINT` | Ends at once (exit 3) without writing anything more. |

A batch system that sends SIGTERM before it kills a job (or SIGUSR2 a set
time before the end) thus gets a clean stop that the next job continues
with `--restart latest`.

## Limits (D83)

`[run.limits]` stops a run with exit code 3 when it reaches a budget:

| Key | Unit | Code |
|---|---|---|
| `max_wall_time` | s | `RUN003` |
| `max_steps` (transient) | | `RUN004` |
| `max_iterations` (steady, outer iterations) | | `RUN004` |
| `max_output` (bytes in the output directory) | B | `RUN005` |

The limits may change during a run (a re-read of the case file, D40).

## Files of a run

| File | Format | Contents |
|---|---|---|
| `result.json` | `foehncraft.result` | The result summary (D62): status, exit code, stop, overrides, provenance, diagnostics, one entry per solve (mesh hash, iterations, residuals, mass imbalance, boundaries, coefficients, errors, transient summary), observed orders, the discretisation uncertainty of a sweep of three or more meshes or time steps (`uncertainty`, D67: per quantity the values, ratios, order, extrapolated value, band, GCI and verdict), the output files with sizes and SHA-256 hashes, and `limitations` (D199: `code` LIM001 to LIM006, `feature`, `message`; below). No timings: two identical runs give the same bytes. |
| `run.jsonl` | `foehncraft.run_record` | The run record (D69), one JSON object per line, flushed at every log line, so it can be followed while the run goes on (`tail -f`). |
| `flow.log`, `flow.md`, `flow.json` | | The log and the report for people, with timings. |
| `fields.vtkhdf`, `fields_N.vtkhdf` | VTKHDF 2.x | Fields, with the root attributes `foehncraft_format`, `foehncraft_format_version`, `foehncraft_version`, `case`, `units` and `step`/`time` or `iteration`. The derived fields of `[output] derived` (D72) are cell datasets of the `internal` block, as the passive scalars (`[scalars.<name>]`, one dataset per scalar under its name, and `<name>_error` with an `exact` solution); the wall shear stress and y⁺ are cell datasets of the boundary blocks (one value per boundary face, zero off walls). `units` names every dataset. |
| `derived.vtkhdf`, `derived_N.vtkhdf` | VTKHDF 2.x | Written by `foehn post derive`; the same layout and attributes as `fields.vtkhdf`. |
| `restart.h5` | Foehncraft restart | With `foehncraft_version` and `units`; transient files carry the attribute `dt_prev2` (the step before the last; read as 0 from older files) and, with the error estimate of the time step (D65), the datasets `u_old2_x`, `u_old2_y`, `u_old2_z`, and with a convective outlet (D29) the old boundary velocity `ub_old_x`, `ub_old_y`, `ub_old_z`; optional datasets `point_id`, `points_x`, `points_y`, `points_z` and `volume_<k>` for meshes that move (D91), left out for fixed meshes. |
| `*.csv` | `foehncraft-<kind>` | Residuals, history, forces, probes, boundaries, wall and line data, each with a comment preamble (below). |
| `balance.csv` | `foehncraft-balance` | The conservation ledger (D88, D89, below): one row per time step (every `output.history_every`) or per balance table of a steady run and its final state. |
| `heat.csv` | `foehncraft-heat` | Transient runs with the energy equation (stage 7d6): one row per time step (every `output.history_every`), `step,time,conduction,convection,source,storage`, the heat stored per region `<region>_storage` (`fluid` first), the heat flow into the fluid through each interface `interface-<region>_heat`, `imbalance` (conduction + convection + source − storage, W) and `relative` (over the largest heat flow). |
| `balance.svg` | SVG | Bars of the kinetic-energy terms at the end: of the final state (steady, W) or over the run (transient, J). `tools/plot/balance_chart.py` draws the time histories of a transient `balance.csv`. |
| `spray.csv` | `foehncraft-spray` | Runs with `[spray]` (SP1, D265): one row per time step, the mass book (`parcels`, `injected_kg`, `in_flow_kg`, `deposited_kg`, `escaped_kg`, `lost_kg`, `lost_parcels`, `rebounds`, `substeps`), the centre of mass of the liquid in the flow `com_x`, `com_y`, `com_z`, its Sauter diameter `d32_m`, and per injector `<name>.parcels`, `<name>.in_flow_kg`, `<name>.penetration_m` (the length within which `[spray] penetration_fraction` of its liquid lies, along the injector axis or radial) and `<name>.farthest_m`. |
| `spray/parcels_N.vtkhdf` | VTKHDF 2.x PolyData, `foehncraft.parcels` | The parcels of step N as vertices: point fields `velocity`, `diameter`, `droplets` (droplets per parcel), `age`, `parcel_id` and `injector` (the index of the injector in the order of the case); root attributes `case`, `step`, `time`, `liquid`, `liquid_density`, `parcels`, `units`. Read by `foehncraft.read_parcels` and `foehn call post.droplets`. |
| `datums.csv` | `foehncraft-datums` | Runs with `[[datum]]` (SP1b, D262 to D267): one row every `output.monitor_every` steps and at the end, the column group `<datum>.<quantity>` of every datum. Surfaces: `area`, `mass_flow` (along the datum's direction), `volume_flow`, `momentum_flow_x/y/z`, area averages `velocity_*_area`, `pressure_area`, `<scalar>_area` and mass-flow averages `<scalar>_mass`, and for a cut plane `cut_error` (relative difference to the same plane moved to the faces); volumes: `volume` and `*_volume`. With sprays: `liquid_kg_pos`, `liquid_kg_neg` (cumulative mass through the surface each way, from the crossing records, the same bits on any thread count), `liquid_rate` (net over the interval), `crossings`, `d10`, `d32`, `dv50` of the interval; volumes `liquid_kg` and `parcels`. The units are in the preamble. |
| `spray_hist_<datum>.csv` | `foehncraft-spray-histogram` | Runs with sprays and datums: per interval and bin `d_lo`, `d_hi`, `droplets`, `mass` of the crossings (surfaces) or the parcels inside (volumes); bin −1 and bin `bins` hold what lies below and above the histogram. `[[datum]] histogram` sets the bins (default 40 logarithmic bins from 1 µm to 1 mm). |
| `series.json` | `foehncraft.series` | Transient runs (D268): one entry per output step, `{"step", "time", "files": {"fields": …, "parcels": …}}`, paths relative to the folder, sorted by step, rewritten through a temporary file after every output so a running case can be followed; a restart keeps the entries before its first step. The time player of MurmelCraft reads it. |
| `case.read_toml_N.toml` | TOML | The settings in force after a re-read (D40). |
| `divergence.md`, `divergence.json` | `foehncraft.divergence` | After a divergence (D63, below): the report for people and for programs. |
| `divergence.vtkhdf` | VTKHDF 2.x | The failing cells and three layers of neighbours as a mesh of their own, with the failed fields, the last valid fields, `divergence_layer` (0 failing, 1–3 neighbours) and the global `cell_id`. |
| `fields_last_valid.vtkhdf` | VTKHDF 2.x | After a divergence: the fields of the last valid iteration or step (`restart.h5` holds the same state). |

### Rows of `run.jsonl`

| `type` | Fields |
|---|---|
| `header` | `format`, `format_version`, `case`, `date`, `provenance`, `threads`, `overrides`, `case_file`, `settings` (the resolved case) |
| `iteration` | `iteration`, `momentum` (3), `continuity`, `bicgstab` (3), `cg`, `scalars` (one residual per passive scalar), `turbulence` (residuals of k and ω or ε, empty when laminar), `turbulence_floors` (cells where the positivity floor of k and ω or ε acted) |
| `step` | `step`, `time`, `dt`, `cfl`, `momentum` (3), `continuity`, `outer`, `cg`, `scalars`, `turbulence`, `turbulence_floors`; with the error estimate of the time step (D65) also `error`, the estimated local error over the tolerance; with the energy equation also `heat` (the row of `heat.csv`: `conduction`, `convection`, `source`, `storage`, `storage_regions`, `imbalance`, `relative`, and `interfaces`, the heat flow of each interface in the order of `heat.csv`) |
| `columns` | (1.1) `table` (`iteration` or `step`), `first` (name of the counter column), `columns`: per column `id`, `name`, `group`, `unit`; one row per table, before its first `table_row` |
| `table_row` | (1.1) `table`, `index` (iteration or step), `values`: column id to value (a number, `null` for a column without a value, a string for NaN or Inf; an extreme carries `{"cell", "x"}`) |
| `section` | (1.1) `title`, `rows` with `label`, `value`, `unit`: one framed block of the log (files, run control, mesh, fluid, settings, boundaries, balances, results, run time) |
| `event` | `code`, `severity`, `message`, `at` (iteration or step); a change by `[[schedule]]` (`RUN015`) also has `key`, `old`, `new` (1.1) |
| `balance` | `iteration` (steady) or `step` and `time` (transient), `mass_outflow`, `mass_relative`, `kinetic_energy`, `storage`, `boundary_power`, `body_force_power`, `viscous_dissipation`, `numerical_dissipation`, `numerical_relative`, `pressure_dissipation`, `scheme_dissipation` |
| `balance_run` | the same fields for the whole transient run (time integrals: kg, J) |
| `budget` | (1.1, D198) after the header: `iterations` (steady: the iteration budget with `iterations`, `source`, `diagnose_at`, `diagnosis`), `estimate` (the D81 estimate: `cells`, `wall_time`, `cell_hours`, `gpu_hours`, `output_bytes`) and, with a study budget, `study` (`ledger`, `runs`, and per budget `cell_hours`, `gpu_hours`, `output`: `budget`, `spent`, `remaining`, `estimate`, `fits`) |
| `diagnosis` | (1.1, D198) a steady run at `[run.limits] diagnose_at` that has not converged: `iteration`, `window`, `budget`, `tolerance`, `controlling` (equation with the largest residual), `verdict` (`on-track`, `slow`, `stalled`, `rising`, `oscillating`), `equations` (per equation `equation`, `residual`, `decades_per_1000`, `scatter`, `trend`, `tolerance_at`), `advice` |
| `spend` | (1.1, D198) with a study budget, before `result`: `cell_hours`, `gpu_hours`, `output_bytes` of this run, as written to the study ledger |
| `timings` | `seconds`, `cells`, `steps`, `ns_per_cell_step`, `solves` (per solve: `mesh`, `device`, `timings`, `ns_per_cell_iteration`, `ns_per_cell_step`) |
| `result` | `result` (the object of `result.json`) |
| `end` | `status`, `exit_code` |

A reader skips a `type` it does not know; new types and new fields do not
raise the major version.

### Self-describing files (D75)

Every CSV starts with comment lines that a reader may skip (`#`):

```
# format foehncraft-history 1; Foehncraft 0.5.0; case cylinder_re100; step 400; time 40 s
# units: time s, Fx N, Fy N, Fz N
time,Fx,Fy,Fz
```

### Byte-identical output (D84)

Two runs of the same case with the same build give byte-identical files,
on any thread count and on CPU and GPU (D31), except the measured files
`flow.log`, `flow.md`, `flow.json` and `run.jsonl`, which hold timings.
`result.json` holds nothing measured and is identical too; it lists every
other output file with its size (`bytes`) and BLAKE3 hash (`blake3`;
`sha256` in results of format 1.x, before 2026-11-24, and from a build
without the feature `blake3`, which the C-free build of D87 leaves out). BLAKE3 hashes a
file of many gigabytes in pieces on every core: SHA-256 took 17 s for a
5 GB file on a CPU without SHA instructions. HDF5 files are written
without object timestamps. Dates (in reports and in `run.jsonl`) come from
`SOURCE_DATE_EPOCH` when it is set.

## Conservation ledger (D88, D89, D90 part 4)

Every run books mass and kinetic energy. Per quantity: the storage (time
derivative of the integral over the domain), the boundary terms per boundary,
the volume sources, and the imbalance, relative to the largest term:

| Term | Mass (kg/s) | Kinetic energy (W) |
|---|---|---|
| storage | 0 (constant density) | (Kⁿ⁺¹ − Kⁿ)/Δt, K = ∫ ½ρ\|u\|² dV |
| boundary, per boundary | mass flow out | viscous work in, Σ (μ∇u·S)·u_b; convected out, Σ ṁ ½\|u_b\|²; pressure work out, Σ p_b u_b·S |
| sources | | body force, ∫ f·u dV |
| sinks | | viscous dissipation of the discrete viscous operator (converges to ∫ μ∇u:∇u dV) |
| imbalance | net outflow | numerical dissipation = boundary + sources − viscous dissipation − storage |
| of it | | pressure part: boundary pressure work − Σ u_P·∇p_P V_P (the momentum rows hold the pressure as the volume source −∇p_P V_P), the discrete −∫ p ∇·u dV of the cell velocities; scheme part: the rest, of the convection and time schemes |

The second law asks that the numerical dissipation is not negative (D90
part 4). The pressure part is the dissipation of the momentum
interpolation of the face fluxes; on distorted cells it holds a
consistency error of either sign that falls with the mesh. Transient runs
use the trapezoidal rule for the kinetic-energy rates over a step, the
mass flows at the end of the step (the discrete continuity equation holds
at the new level), and add a ledger of the whole run (time integrals, kg
and J, from the start or the restart). All sums are pairwise sums in a fixed
order (D31). `[balance]` sets the cadence and the tolerances:

| Key | Default | Meaning |
|---|---|---|
| `every` | 0 | The full table in the log every N iterations or steps (0: at the end). Transient runs also log one compact line with every progress line. |
| `tolerance_step` | 1e-8 | Relative mass imbalance of a step, or of the final steady state: `BAL001`. |
| `tolerance_run` | 1e-4 | Relative mass imbalance over the run: `BAL002`. |
| `tolerance_energy` | 1e-2 | Numerical dissipation below −tolerance relative to the largest term, of the final steady state or over the run: `BAL003`; a transient step below it counts as negative. |
| `action` | `warn` | `stop` ends a transient run at the step with exit code 3. |

`balance.csv` columns: `step,time` (transient) or `iteration` (steady),
then `mass_outflow,mass_relative,kinetic_energy,storage,boundary_power,
body_force_power,viscous_dissipation,numerical_dissipation,numerical_relative,
pressure_dissipation,scheme_dissipation`.
`result.json` carries per solve `balance`: `last` (the final state or the
last step, with the terms per boundary), `run` (transient), `steps`,
`negative_steps` (steps whose numerical dissipation is below
−`tolerance_energy`:
a health number, not an error) and `largest_mass_relative`.
Each balance also carries `ledger` (D119): one entry per account
(`mass`, `kinetic energy`, with the energy equation `thermal energy`, and
the accounts a module books) with `quantity`, `unit`, `region`, `content`,
`storage`, `boundary`, `contributions`, `imbalance`, `relative` and the
`terms` (`name`, `kind` = `boundary`/`source`/`dissipation`/`exchange`/
`part`, `boundary`, `pair`, `value`), and `pairs` (`pair`, `residual`), the
sum of the two sides of every exchange pair. The thermal account is in
`last` only. Both are additive; the fixed fields above are unchanged.

## Divergence report (D63)

When a residual of an outer iteration or corrector is not finite, or a
linear solver breaks down, the run puts the state of the last valid
iteration or step back (kept as a device copy of u, p and ṁ before every
iteration or step) and writes `fields_last_valid.vtkhdf`, `restart.h5` of
that state (a transient restart continues, not bit-identically: the old
gradients of the deferred corrections are those of the last valid level),
`divergence.vtkhdf`, `divergence.md` and `divergence.json`, then ends with
exit code 1 and `RUN014`, whose fix is the one of the most likely cause.
`divergence.json` (format `foehncraft.divergence`, version 1.0):

| Key | Contents |
|---|---|
| `kind`, `at`, `done`, `time`, `time_step` | steady or transient, the failed iteration or step, the iterations or steps before it, the time and step of the last valid state |
| `first_variable` | `velocity` (the momentum predictor failed) or `pressure` |
| `residuals` | of the failed iteration or corrector (zero when a linear solver broke down before they were complete) |
| `failing_cells`, `failing_criterion`, `cells` | the count, what made them fail, and the first 20 by global ID with rank, centroid, failed velocity and pressure, last valid speed |
| `nearest` | the cell nearest to failing in the last valid state: largest Courant number (transient) or speed (steady), and its distance to the failure |
| `causes` | most likely first: `code` (`DIV001` time step, `DIV002` relaxation, `DIV003` cell quality, `DIV004` boundary condition, `DIV005` initial state, `DIV006` linear solve), `score` 0–1, `evidence`, `advice` and `fix` (`--set` arguments) |
| `files` | the files written |

The bounded quantities of later models (turbulence, species, volume
fractions) get range warnings at the progress lines with those models;
the laminar flow has none.

## Provenance (D85)

The `provenance` object in `result.json`, in the header of `run.jsonl` and
in `foehn info --json`: `version`, `commit` (`-dirty` with uncommitted
changes), `features`, `backends`, `backend`, `platform`, `cpu`, `gpu`,
`ranks`, `gpus`, `reproducible_mode`, `case_hash` (the SHA-256 of the
resolved case without `run.threads`, `run.backend`, `run.ptx` and
`output.directory`, which do not change the results, with the contents of
the profile and expression tables and the resolved property program),
`properties_hash` (the SHA-256 of the resolved property program of D94
alone) and `materials_database` (the SHA-256 of the built-in material
database). Two results with the
same `case_hash` and `commit` are comparable bit for bit.

## JSON objects of the command line

| `format` | From |
|---|---|
| `foehncraft.result` | `foehn flow --json` (also `result.json`) |
| `foehncraft.check` | `foehn flow --check --json`: `case_file`, `ok`, `diagnostics`, `estimate` (cells, faces, steps, host and device memory, output bytes and files, wall time and its basis) |
| `foehncraft.error` | any command with `--json` that fails before a run: `command`, `exit_code`, `diagnostics` |
| `foehncraft.info` | `foehn info --json` |
| `foehncraft.command` | the other commands with `--json` (`mesh`, `scalar`, `poisson`, ...): `command` and the report; `post derive`: `file`, `restart`, `cells`, `iteration`, `step`, `time`, `derived`, `datasets` |
| `foehncraft.divergence` | `divergence.json` after a divergence (above) |
| `foehncraft.verification-grades` | `docs/verification/grades.json` (D199, generated): per feature `id`, `title`, `grade` (`verified`, `partial`, `unverified`), `why` (the reasons for a lower grade), `gaps`, `range` (validated range: `y_plus_first`, `reynolds`, `prandtl`, `courant`, `skewness`, `mach`) |

## Limitations of a run (D199)

Every result summary (`result.json`, `foehn flow --json`, the Python
`Result.limitations`) carries a list `limitations`, computed at the end of
the run from the case, the report and the grades of the
[verification status](verification/status.md) page compiled into the
binary (`docs/verification/grades.json`). Each entry has `code`,
`feature` (an id of `docs/verification/features.toml`, where one applies)
and `message`:

```json
"limitations": [
  {"code": "LIM001", "feature": "mesh.gmsh", "message": "Gmsh v4 reader is graded partial: ..."},
  {"code": "LIM004", "message": "one mesh only: the discretisation uncertainty is not estimated ..."}
]
```

| Code | Limitation | Taken from |
|---|---|---|
| `LIM001` | a feature in use is graded partial or unverified | the switches of the case against the grades |
| `LIM002` | a feature in use ran outside its validated range | `range` of the feature against the run: mean first-cell y⁺ of the walls, Reynolds number, Prandtl number (energy, passive scalars), largest Courant number, largest face skewness |
| `LIM003` | a steady solve did not settle | a steady solve that stopped at its iteration limit before its tolerance, when no check already failed for it |
| `LIM004` | one mesh only | no level sweep, so no discretisation uncertainty (D67) |
| `LIM005` | a GPU not verified | a GPU not in the verified devices (Tesla T4) |
| `LIM006` | automatic mesh remediation active | the remediation of D50 changed the schemes (MESH002) |
| `LIM007` | the history before t = … s is not physical (thermal acceleration) | `[energy.acceleration]` (D213): the run reports `physical_from` |
| `LIM008` | an operator of a feature in use is not exact on N of M cells | the operator check (D236): `MESH005` with `turbulence.curvature_correction`, `MESH006` with `turbulence.qcr` |
| `LIM009` | the solids carry the cycle-mean state (or its history with `solid = "transient"`), the swing over a cycle only as an estimate | `[energy.cycle]` (D204) |
| `LIM010` | the flow was held during the march (frozen flow) | `[energy.frozen_flow]` (D118): the run reports the frozen steps and refreshes |

Limitations change neither `status` nor the exit code (D199 Q1, Q2): a
run with limitations is `finished` with exit code 0. The human text says
so: the last frame of the log reads "Run finished with N limitations" and
lists them, and `flow.md` starts with a "Limitations" section.

## Diagnostics (D78)

A diagnostic is an object:

```json
{
  "code": "CASE002",
  "severity": "error",
  "key": "mesh.generate.cels",
  "location": [4, 1],
  "message": "unknown key `cels` in [mesh.generate]; did you mean `cells`?",
  "fix": null
}
```

`key` is the case key (`section.key`), `boundary` the boundary of a
boundary-condition finding, `location` line and column in the case file
(for syntax errors and unknown keys), and `fix` a list of `--set`
arguments that resolve it, where one exists (`{"set": ["time.adaptive=false"]}`).
A code keeps its meaning; a retired code stays in the table.

| Code | Severity | Title |
|---|---|---|
| `CASE001` | error | The case file is not valid TOML. |
| `CASE002` | error | A key that the case file does not have (with the nearest key as a hint). |
| `CASE003` | error | A value of the wrong type or an unknown choice. |
| `CASE004` | error | A value outside its range. |
| `CASE005` | error | Settings that do not fit together. |
| `CASE006` | error | A file that cannot be read: the case, the mesh or the restart file. |
| `CASE007` | error | A path outside the read roots or the write roots (D82). |
| `CASE008` | error | A `--set key=value` that cannot be applied. |
| `CASE009` | error | A key that the case needs is missing. |
| `CASE010` | error | A case file that its reader refuses for a reason without a more specific code yet (the message says which); D202. |
| `SCHED001` | error | An entry of `[[schedule]]` that is not valid: its trigger, its order, a key that cannot change during a run, or a value that does not validate (D128). |
| `BC001` | error | A mesh boundary without a boundary condition. |
| `BC002` | error | A boundary condition for a boundary the mesh does not have. |
| `BC003` | error | A condition that a later stage brings. |
| `BC004` | error | A velocity inlet with none or more than one of `velocity`, `normal_velocity` and `profile`. |
| `BC005` | error | `profile = "analytic"` without an `[analytic]` table. |
| `BC006` | error | `profile = "parabolic"` without `max_velocity` or `profile_axes`. |
| `BC007` | warning | `max_velocity` or `profile_axes` without a parabolic profile. |
| `BC008` | warning | A mass-flow inlet with a negative mass flow. |
| `BC009` | warning | Turbulence values in a laminar case. |
| `BC010` | warning | A wall velocity with a component normal to the wall. |
| `BC011` | note | The exact solution leaves the domain through a velocity inlet. |
| `BC012` | warning | A velocity inlet whose velocity points out of the domain. |
| `BC013` | warning | A symmetry boundary that is not planar. |
| `BC014` | error | A net mass flow through the prescribed boundaries and no pressure boundary to balance it. |
| `BC015` | note | No pressure boundary: the pressure level is fixed in the reference cell. |
| `BC016` | error | A flow path whose flow rate is over-specified. |
| `BC017` | warning | Nothing drives a flow: the result is the fluid at rest. |
| `BC018` | note | No prescribed inflow: the flow rate follows from the pressure boundaries. |
| `BC019` | error | A table inlet profile that is missing, unreadable or not a valid table. |
| `BC020` | warning | Faces of a table inlet outside the range of the table (they take the end values). |
| `BC021` | error | A convective outlet in a steady case. |
| `BC022` | error | A convective outlet without `convection_velocity` and no prescribed inflow to take it from. |
| `BC023` | error | Wedge boundaries that do not form one wedge about an axis. |
| `BC024` | warning | A wedge wider than 10° (the axisymmetric error grows with the angle). |
| `BC025` | error | A potential-flow start with a convective outlet and no pressure boundary. |
| `BC026` | error | Walls that add or remove mass in a closed volume (no pressure boundary, no uniform thermodynamic pressure). |
| `BC027` | warning | A wall roughness k_s larger than the first-cell distance on more than 10 % of the wall's area. |
| `BC028` | warning | A turbulence value that the run does not use (on a boundary other than an inlet or a pressure outlet, whose backflow takes the values; `zeta` with a model other than k-ε-ζ-f, ζ-f₀ or the elliptic-blending k-ε); `omega` with an ε model and `epsilon` with SST are converted (ε = C_μ k ω); BC009 stays for laminar cases. |
| `BC029` | error | Two incompatible turbulence specifications on one inlet (for example `turbulence_length_scale` and `viscosity_ratio`, `k` with intensity, or `omega` with `epsilon`). |
| `BC030` | error | A periodic boundary without a periodic partner that names it back, a pair whose faces do not match, or a driving given twice. |
| `BC031` | warning | A periodic pair with no driving (`mass_flow`, `pressure_gradient`) and no body force: the flow decays to rest. |
| `SCAL001` | error | A `[boundary.<name>.scalar.<scalar>]` table for an unknown scalar or boundary, without exactly one of `value` and `flux`, or with a condition its boundary does not take (a flux at an inlet, anything at an outlet, `periodic` away from a periodic pair, the heat-flux mode with fixed values). |
| `TURB001` | warning | First-cell y⁺ outside the claimed range of the wall treatment on more than 10 % of a wall's area: above 2 with `low-re`, above 300 (above the log layer) with `automatic`, below 30 with `scalable`, outside 30 to 300 with `standard` and `compound` of the k-ε family (per boundary with the area-weighted min, mean and max). |
| `TURB002` | warning | A positivity floor of k, ω or ε acted on more than 0.1 % of the cells in an iteration of the final tenth of the solves (equation and cell count); the run record counts all iterations (`turbulence_floors`). |
| `TURB003` | warning | Inlet turbulence outside the usual range: I > 0.5, or a viscosity ratio below 10⁻² or above 10⁵. |
| `TURB004` | note | Kato–Launder production in a rotating zone: Ω is the absolute-frame vorticity. |
| `TURB005` | note | No wall boundary: the wall distance is infinite (F₁ = 0, F₂ = 0; the ζ-f and EB near-wall terms vanish). |
| `TURB006` | error | A convection scheme or time scheme that is not bounded for a turbulence variable (central in a steady run, Crank–Nicolson). |
| `TURB007` | warning | Model constants changed in `[turbulence.constants]`. |
| `TURB008` | warning | The ν_t cap `max_viscosity_ratio` acted (cell count). |
| `TURB009` | note | Inlet turbulence defaulted to I = 0.05, r = 10 (boundary named). |
| `TURB010` | note | The curvature correction of SST together with Kato–Launder production: both act on stagnation and swirl regions (the paper uses the standard production). |
| `EXPR001` | error | An expression that cannot be parsed. |
| `EXPR002` | error | An unknown or reserved name in an expression, a constant or a table (with the nearest name as a hint). |
| `EXPR003` | error | An expression that uses a variable its key does not take (`t` in a steady case or in a property, x, y, z in a value for a whole boundary). |
| `EXPR004` | error | A `[tables]` entry that is missing, unreadable or not a valid table. |
| `EXPR005` | error | An expression or constant whose value is not finite where it is evaluated. |
| `MAT001` | error | A material that neither the case nor the database defines (with the nearest name as a hint). |
| `MAT002` | error | A property form that this version does not evaluate, or a property that the flow cannot take yet (a density that varies). |
| `MAT003` | warning | A property evaluated outside the validity range of its data (clamped to the range). |
| `MAT004` | error | Materials whose `base` refer to each other in a cycle. |
| `MAT005` | error | A property value that is not positive and finite. |
| `MAT006` | error | A solid material without density, specific heat or thermal conductivity. |
| `MAT007` | error | A fluid without specific heat or thermal conductivity in a run with the energy equation. |
| `MAT008` | error | An orthotropic thermal conductivity without both axes, or with parallel axes. |
| `MAT009` | error | A user material database file (`[materials] include`, `FOEHN_MATERIALS`) or a property data file (NASA-7 record, table file) that cannot be read or has the wrong format. |
| `MAT010` | error | A water or steam state (`materials.steam`, `foehn steam`) outside IAPWS-IF97, or asked with a set of inputs other than (p, T), (p, h), (p, s), (p, x), (T, x), p or T alone. |
| `SPR001` | error | A `[spray]` value that is out of range, an injector outside the mesh, a spray file (`[spray] include`) that cannot be read or defines `injector` or `liquid` twice, or a `[spray.boundary.<name>]` for a boundary the mesh does not have. |
| `SPR002` | error | The steady spray mode (D260): `[spray.steady]` in a case with `[time]`, or an injector with `start` or `end` in a steady case (its mass flow must not change in time). |
| `SPR003` | error | A `[spray]` in a case with MRF zones, sliding interfaces or a periodic boundary, or two-way coupling with a setting it does not work with yet (solid regions, the frozen-flow march, an algorithm other than PISO), or the steady spray mode with solid regions or, two-way, with the coupled algorithm (not in this version). |
| `SPR004` | warning | Parcels lost: the walk through the mesh failed or a parcel needed more than `[spray] max_substeps` sub-steps in a time step (counted in `spray.csv`). |
| `SPR005` | error | The mass book of the sprays does not close: released liquid differs from the liquid in the flow, stuck, left and lost by more than round-off. |
| `SPR006` | warning | Two-way spray coupling (SP2, D252): the parcel passes of a time step ended at `[spray.coupling] max_passes` with the coupling lag above `tolerance` (fix: more passes or a smaller time step; `spray.csv` counts such steps). |
| `SPR007` | warning | Steady spray mode (D260): streams still in the flow after `[spray.steady] max_time` (caught in a recirculation, or a domain the droplets need longer to cross); their mass flow is counted as unfinished in `spray.csv` and `result.json` (fix: a larger `max_time`). |
| `SPR008` | warning | Steady spray mode (D260 (4)): the source residual of the final stream pass is above `[spray.steady] tolerance`, so the run is not reported as converged (fix: more outer iterations, or a smaller `[spray.steady] relaxation` if the residual swings). |
| `SPR009` | warning | The saturation limiter of D253 held vapour back: the evaporation of a spray cell was above what brings the gas the cell holds and the gas that flows in to adiabatic saturation (transient: in a step; steady: in the final stream pass). Counted in `spray.csv` (`limiter_cells`, `limiter_kg_s`) and `result.json` (fix: more passes or a smaller time step; steady, more outer iterations or a smaller `[spray.steady] relaxation`). |
| `IMM001` | error | A triangle surface (STL) with no triangles or with corners whose coordinates are not finite. |
| `IMM002` | warning | An open triangle surface: edges with only one triangle (holes, open borders). Inside and outside still come from the generalised winding number, which closes small holes smoothly; nothing is repaired. |
| `IMM003` | warning | A triangle surface with edges shared by more than two triangles (non-manifold: touching bodies, duplicate or internal triangles). |
| `IMM004` | warning | A triangle surface whose orientation is not consistent: edges whose two triangles run them in the same direction (some triangles flipped). |
| `IMM005` | warning | A self-intersecting triangle surface: pairs of triangles that cross or touch away from their shared corners and edges (exact orientation predicates, D225). The winding number still decides inside and outside. |
| `IMM006` | warning | Degenerate triangles: corners that are repeated or exactly collinear (zero area). |
| `IMM007` | warning | Duplicate triangles: two triangles with the same three corners. |
| `IMM008` | warning | A closed surface whose normals point inward (negative enclosed volume): as given, the body is its complement (fix: flip the triangles in the CAD export). |
| `IMM009` | note | A surface hint: a bounding box larger than 20 m (a file in millimetres needs `scale = 1e-3`), or distinct corners closer than 1e-9 of the box diagonal (welding candidates, not welded). |
| `IMM101` | error | An `[immersed]` value that does not fit: a permeability that is not > 0, a body without the keys of its shape (`centre` and `radius`, `lower` and `upper`), a radius or length that is not > 0, a zero cylinder axis, two bodies with the same name, a `field` body without `fraction`, no body at all, or an `[immersed.topology]` value out of range (an empty design region, a fluid fraction not in (0, 1], a negative filter radius). |
| `IMM102` | error | An `[immersed]` body in a case it does not fit (L13 stage 1, D226): a turbulence model (penalised walls are laminar and low-Re only; a wall function on a penalised wall is refused), the energy equation (no penalised conjugate heat transfer, Q9), passive scalars, gravity, regions or several fluids, rotating zones, or sprays; a topology search in a transient run. |
| `IMM103` | warning | An immersed body that the mesh does not resolve: no cell lies more than half inside it, or its solid volume on the mesh differs from its exact volume inside the mesh by more than 5 %. |
| `IMM104` | note | The immersed bodies of the case: their solid volume on the mesh against the exact volume, the cells they touch and the sink coefficient μ/K; the design cells of a topology search. |
| `IMM201` | error | The cut of `[immersed] method = "cut"` failed: the bodies fill the whole mesh, a body has the name of a mesh boundary, the mesh is a rank's part or periodic before the cut, or an `[immersed.cut]` value is out of range (`min_volume_fraction` not in [0, 0.9]). |
| `IMM202` | error | Cut cells (L13 stage 2) with a setting they do not have yet: regions or conjugate heat transfer (stage 2d), a topology search (it needs the penalised density field), rotating zones or a moving mesh (rigid motion is stage 5). |
| `IMM203` | note | The cut of the immersed bodies: cells cut, removed, split into pieces and merged, the smallest kept volume fraction, the wall area of every body against its exact area where simple, and the closure of the cut cells. |
| `IMM204` | warning | Cut cells the merge could not fix: a small cut cell without a neighbour to merge into, or faces whose normal does not point from the cell centroid to its neighbour (the run refuses them). Fix: refine the mesh near the body, or move the body by a small fraction of a cell. |
| `DAT001` | error | A `[[datum]]` table that cannot be used: a name used twice or with other characters than letters, digits, `_` and `-`, an unknown type, a missing or invalid key of its type (plane point and normal, bound, boundary, zones, STL file, box corners, cylinder, zone), or bad histogram bins. |
| `DAT002` | error | A `[[datum]]` that meets no face or cell of the mesh (a plane outside the domain, an interface between zones that do not touch, an empty volume). |
| `ENE006` | error | An energy setting this version does not have: a passive scalar named `temperature` beside `[energy]`, solid regions with a sliding mesh, a restart with `[energy.acceleration]`, passive scalars or a periodic mesh, or with the variable density (stage 7c) passive scalars, a restart, a periodic mesh, a convective outlet, a wall mass flux or a rotating zone. |
| `ENE008` | error | An `[[energy.sources]]` entry for a region the case does not have, or without exactly one of `power` and `power_density`. |
| `ENE009` | error | Variable-density settings that do not fit together: a density that does not depend on T, `density = "variable"` with the Boussinesq model, a closed volume with an inlet or outlet, a closed volume of a liquid, or `closed-volume` without the variable density. |
| `CHT004` | error | An `[interfaces.<name>]` table that is not on a region interface (regions unknown, the same, two fluids, or not touching), or without exactly one of `contact_resistance` (≥ 0) and `thin_wall` (thickness > 0 with `conductivity` or `material`). |
| `ACC001` | error | An `[energy.acceleration]` setting that cannot be used: a steady case, no solid region, a region that is not a solid, `inertia_scale` outside (0, 1] and not `"auto"`, `restore` neither `"settled"` nor a time ≥ 0, or a `period`, `window_steps`, `settle_tolerance` or `target_windows` out of range. |
| `ACC002` | warning | The run ended before the true heat capacity of the accelerated solids was restored: the final state is that of the accelerated solids and the history is not physical. |
| `ACC003` | note | The accelerated solids' skin depth of the period, √(2k/(f ρ c_p ω)), is larger than the thickness V/A of a region (V its volume, A its interface area): while accelerated the swing of the period reaches through the wall; it goes away with the restore. |
| `ACC004` | note | A step of the thermal acceleration: the conductances and factors measured, the solids settled, the restore started or finished. |
| `CYC001` | error | An `[energy.cycle]` setting that cannot be used: a steady case, no solid region, together with `[energy.acceleration]`, an adaptive time step, a `period` that is not a whole number (≥ 2) of time steps, or a count, tolerance or relaxation out of range. |
| `CYC002` | warning | A fluid phase of the cycle-averaged loop ended at `max_cycles` while the mean interface heat flow still changed between the last two cycles by more than `cycle_tolerance`: the averages are not those of a periodic flow. |
| `CYC003` | warning | `skin = "effusivity"` with a skin depth of the period √(2k/(ρ c_p ω)) that is not small against the thickness V/A of a solid region (larger than a third of it): the swing reaches through the wall there and the semi-infinite solution does not hold; `resolved_cycles` gives the exact swing. |
| `CYC004` | warning | The cycle-averaged loop stopped at `max_updates` before the wall temperatures and the interface heat flows settled within `tolerance` and `heat_tolerance`. |
| `CYC005` | note | A step of the cycle-averaged loop: a fluid phase averaged, a solid phase solved, the wall temperatures updated, the loop converged. |
| `FRZ001` | error | An `[energy.frozen_flow]` setting that cannot be used: a steady case, no solid region, together with `[energy.acceleration]` or `[energy.cycle]`, the variable density, an adaptive time step, a `time_step` below the `[time]` step, or a `start`, window, tolerance or count out of range. |
| `FRZ002` | warning | The flow of a frozen-flow run did not settle by the window test (an unsteady flow, or a tolerance too tight): the run stayed fully coupled with the `[time]` step from the last freeze or refresh on. |
| `FRZ003` | note | A step of the frozen-flow march: the flow settled and was held, or it was refreshed by coupled steps because the temperatures it depends on had moved. |
| `FRZ004` | warning | A frozen-flow run whose flow boundary values or body force depend on time: the held flow does not follow them between refreshes (`refresh_every`). |
| `ENE004` | warning | A Mach number above 0.3 in a run with the variable density: the low-Mach form loses its validity (the density does not follow the dynamic pressure). |
| `ENE007` | warning | A steady run that converged while its heat balance does not close (imbalance above 1 % of the largest heat flow): the temperature residual is small but the energy solve is not converged; tighten `[solver.energy] tolerance` and raise `max_iter`. |
| `THM001` | error | A steady run whose boundaries leave the temperature level undetermined: no inlet and no `temperature` or convection boundary, and no net heat input. |
| `THM002` | error | A steady run with a net heat input (heat flux, heat flow, source or viscous heating) and no way for the heat to leave: no inlet and no `temperature` or convection boundary. |
| `THM003` | error | A thermal condition on a boundary that cannot carry it: a flux or convection at an inlet, any condition at an outlet, a periodic boundary, a symmetry plane or a wedge. |
| `THM004` | error | A thermal value out of range: a temperature or ambient temperature ≤ 0 K, a heat-transfer coefficient ≤ 0, a heat flow that is not finite. |
| `THM005` | warning | Walls without a `[boundary.<name>.energy]` table in a run with the energy equation: they are adiabatic (listed once). |
| `THM009` | error | A `[boundary.<name>.energy]` table without `[energy]`: the energy equation is solved only with an `[energy]` table. |
| `THM010` | error | A `[boundary.<name>.energy]` table for an unknown boundary, without exactly one condition (`temperature`, `heat_flux`, `heat_flow`, or `heat_transfer_coefficient` with `ambient_temperature`), or with `heat_transfer_coefficient` and `ambient_temperature` apart. |
| `THM011` | error | An inlet without `[boundary.<name>.energy] temperature` in a run with the energy equation. |
| `BUO001` | error | `[gravity]` with a vector of zero length, or a buoyancy model without `[gravity]`. |
| `BUO002` | warning | The Boussinesq model outside its range: β ΔT above 0.1 over the temperatures of the case (the density changes by more than 10 %). |
| `BUO003` | error | The Boussinesq model without what it needs: no energy equation, no expansion coefficient (the fluid's density does not depend on T and `[buoyancy] expansion_coefficient` is not given), or a reference value out of range. |
| `BUO004` | warning | `[gravity]` with the energy equation but no buoyancy model (the temperature drives no flow); as information: a turbulence model with buoyancy, whose production of k by buoyancy comes in stage 7e1. |
| `BUO005` | error | A variable density in a closed domain (no inlet, no outlet) without `[energy] thermodynamic_pressure = "closed-volume"`: the mass in the domain would not be conserved. |
| `BUO006` | warning | A pressure boundary with `pressure_kind = "static"` over a height span whose hydrostatic difference ρ_ref g Δz is above the dynamic pressure: a uniform static value drives a flow through the opening. |
| `REG001` | error | A region that is not a fluid where it cannot be: a solid region without `[energy]` (solids are solved by the energy equation, stage 7d1), or a solid material as the `[fluid]`. |
| `REG002` | error | A region that names a cell zone the mesh does not have, or a zone in two regions. |
| `REG003` | error | A fluid part (fluid cells connected without a solid between them) without a pressure boundary in a case of several fluid parts: give it a pressure outlet or inlet; a closed part of its own comes with stage MF-B. |
| `REG004` | error | Two fluids of different materials (or temperatures) share a face: that is a two-phase interface (multiphase, L3); put a solid region or a thin wall between them, or give both regions the same fluid. |
| `REG005` | note | The fluid parts of the case: their number, regions and fluids. |
| `REG006` | error | A setting that does not yet work with several fluids in one case (D194): the variable density, gravity, a rotating zone or frame, passive scalars, a convective outlet or a periodic mass-flow drive. |
| `REG007` | error | A key of `[regions.<name>.initial]` that does not fit: a solid region takes `temperature` only, turbulence start values come as `k` with one of `omega` or `epsilon`, and a start temperature needs `[energy]`. |
| `REG008` | error | A start temperature, k, omega or epsilon of `[regions.<name>.initial]` that is not > 0. |
| `POR001` | error | A porous region without its matrix (`[regions.<name>.porous]` with `viscous`/`inertial` or `honeycomb`), contradicting or invalid values, direction-dependent values without `axes`, a matrix in a region that is not porous, or with `[energy]` a matrix without `solid`, porosity, h or a. |
| `POR002` | error | A setting that does not yet work with porous regions (porous plan, D279): with the energy equation a cycle-averaged or accelerated energy mode or a density that varies, passive scalars, a turbulence model, a rotating zone, the GPU (stage PM5), or direction-dependent resistances whose axes are turned against the mesh axes. |
| `SPE001` | error | Invalid `[species]` settings or compositions: an unknown species (the GRI-Mech 3.0 set), a species named twice, a closure species not in `names`, fractions that do not sum to 1 or leave the closure species negative, both or neither of `mole_fractions` and `mass_fractions`, an inlet without `[boundary.<name>.species]`, `properties = "mixture"` or `diffusion = "lewis"` without `[energy]`, or a variable density without an ideal-gas fluid with its molar mass. |
| `SPE002` | error | A setting that does not yet work with species (stage SF1): regions (solids, several fluids, porous regions); conjugate walls come with SF5. |
| `SPE003` | warning | The transported mass fractions summed to more than 1 somewhere by more than 1e-4, so the closure species was held at zero there: the solve is not converged, the tolerance of `[solver.scalar]` is too loose, or the limiter lets extrema pass (`limiter_epsilon`). |
| `CHM001` | error | An invalid chemistry mechanism or reactor input: an unknown built-in mechanism, a mechanism file that does not parse, a reaction that does not conserve the elements, a species without thermodynamic data, an unsupported rate form (pressure-dependent Arrhenius, Chebyshev, SRI fall-off, a specific collider) or unit, or a reactor input with an invalid composition, temperature, pressure or time. |
| `CHM002` | warning | The stiff chemistry integrator of a reactor did not reach the end of its interval within `max_steps`, or its step became too small: raise `max_steps`, loosen `rtol`, or check the mechanism and the state. |
| `CHM003` | error | A setting that does not yet work with chemistry in the flow (stage SF2b): a steady case, a turbulence model (SF3), sprays, or ranks. |
| `REG009` | error | `temperature = "inlet"` in `[regions.<name>.initial]` for a region whose fluid part has no inlet with `[boundary.<name>.energy] temperature`. |
| `REG010` | note | The start values of the regions and fluid parts and where they come from: given, a boundary value, a Laplace solve between the pressure boundaries of a part, or the inlets. |
| `PRS001` | error | An absolute pressure value (`[fluid] pressure_values = "absolute"`, the default from case format 0.7, D239) that is zero or negative: a boundary `pressure` or `total_pressure`, `[initial] pressure`, `[pressure] reference_value` or an operating pressure; set `pressure_values = "gauge"` for values relative to the operating pressure. |
| `PRS002` | warning | An absolute pressure below 1 kPa in a gas (density below 50 kg/m³): it looks like a gauge value; set `[fluid] pressure_values = "gauge"` if it is one. |
| `PRS003` | error | A fluid part whose operating pressure cannot be taken from its pressure boundaries (absolute values, D239): their values are expressions only; give `[fluid] operating_pressure` or `[regions.<name>] operating_pressure`. |
| `PRS004` | error | A case-wide absolute `[initial] pressure` in a case whose fluid parts have different operating pressures: leave it out (each part starts at its own operating pressure) or give the start pressure per region. |
| `PRS005` | note | The pressure convention of the case file (absolute or gauge) and the operating pressure of every fluid part with where it comes from (D239). |
| `MRF001` | error | An `[[mrf]]` table and a `[[sliding]]` table that turn the same cells (one of them without a `zone`, or both with the same one): a zone is turned by a moving reference frame (steady or transient, D185) or by a sliding mesh (D182), not by both. |
| `MRF002` | error | An `[[mrf]]` table with no or two angular velocities (`omega`, `rpm`), a zero axis, a cell zone the mesh does not have (with the nearest name as a hint), no `zone` where there are several tables, or a zone that two tables rotate. |
| `MRF003` | warning | The interface or a boundary of a rotating zone is not a surface of revolution about its axis: the relative flux of the cells next to it does not close, and the solution there is not that of a frozen rotor (cell count and worst closure error). |
| `SLD001` | error | A `[[sliding]]` table in a steady case: the sliding mesh moves in time (D182); a steady case with a rotating zone uses `[[mrf]]` (D55). |
| `SLD002` | error | A `[[sliding]]` table with no or two angular velocities (`omega`, `rpm`), a zero speed or axis, a `shift` of 0 or of a whole turn, a cell zone the mesh does not have (with the nearest name as a hint), or a second table (one rotor per case so far); with the interface imprinted: `rotor_interface` without `stator_interface`, no `zone`, or a `shift` (the time step is free). |
| `SLD003` | error | The time stepping does not fit the sliding mesh: the time step is not `shift` × 2π/(sectors × Ω) (the right value is given), the step is adaptive or error controlled, `end_time` is not a whole number of steps, or `[verification]` asks for a sweep over levels or time steps; with the interface imprinted the scheme must be fully implicit (`extrapolation = false`, at least 2 `outer_iterations`). |
| `SLD004` | error | The rotor cannot be moved by whole sectors: its cells have no image when it is turned by one sector (it is not a ring of identical sectors about the axis), or a boundary of the rotor is not a wall, slip wall, symmetry plane or wedge. |
| `SLD005` | error | The interface of a sliding mesh by imprinting cannot be imprinted: a boundary or cell zone of `[[sliding]]` is missing or wrong, the two interface boundaries are not on one circular cylinder about the axis or do not cover the same surface, rotor and stator cells share a face or a node, a `[boundary]` table names an interface boundary, or a restart is asked for (not yet available). |
| `SLD006` | warning | The cells next to the imprinted interface of a sliding mesh do not close (|Σ S|/Σ |S| above 1e-8): the faces of the interface are not those of a mesh extruded along the axis, and the mesh flux of the rigid turn leaves a small mass source. |
| `MOV001` | error | A `[mesh_motion]` table that cannot be used: a steady case (the mesh moves in time), a displacement expression that does not compile, or a case with what the deforming mesh of stage MS0 does not work with yet (turning zones, regions, sprays, immersed bodies, periodic boundaries, the check of `[run] ranks`, an error-controlled time step, a sweep over levels), or a restart file whose nodes the displacement does not put there or whose node, cell or face order is another. |
| `MOV002` | error | The prescribed displacement of `[mesh_motion]` turns a cell inside out (its volume is no longer positive); the message gives the time, the cell and its place. |
| `MOV003` | error | The deforming mesh moves a boundary that is not a wall out of its own plane (a face sweeps a volume): only walls move with the mesh; the nodes of symmetry planes, wedges, inlets and outlets may only slide along them. |
| `MOV004` | error | The deforming mesh changes the volume of a closed domain (no pressure boundary) of constant density: an incompressible fluid cannot be compressed; give it a pressure boundary or the variable density of a closed volume. |
| `VOF001` | error | A `[vof]` value that does not fit: not exactly two `[[vof.phase]]` tables, a density or viscosity that is not > 0, a negative surface tension, `thinc_beta` outside (0, 20], `max_courant` outside (0, 1], more than 10 smoothing passes, an α outside [0, 1] (`initial`, a region, an inflow), a region without the keys of its shape, or an inflow for a boundary the case does not have. |
| `VOF002` | error | A `[vof]` case with a setting the first version (MP1) does not have: a steady case, a turbulence model, the energy equation, passive scalars or species, a `[buoyancy]` model, regions or several fluids, sprays, immersed bodies, the coupled algorithm, rotating zones or sliding interfaces, a moving mesh, periodic boundaries, ranks, or a restart. |
| `VOF003` | warning | The volume of phase 1 clipped to keep 0 ≤ α ≤ 1 is above 10⁻⁶ of the phase volume over the run (fix: a smaller `[vof] max_courant` or converged pressure corrections; `vof.csv` counts it per step). |
| `ROT001` | warning | A wall of a rotor does not turn with its zone (a stationary shroud in a rotating zone): its torque is in that of the rotor, but the shaft power −Ω·T is not its work on the fluid; give the walls of the rotor in `walls` of the `[[mrf]]` or `[[sliding]]` table. |
| `ROT002` | error | A `walls` entry of `[[mrf]]` or `[[sliding]]` that is not a wall of the case, or an `[[output.machines]]` table without a turning zone, with a `rotor` that is no turning zone (or none where there are several), without an inlet or an outlet, with a boundary the case does not have, or with a name used twice. |
| `MESH001` | warning | Faces whose skewness exceeds `solver.skewness_limit` (poor cells, typically tetrahedra of a mesh that was not optimised). |
| `MESH002` | note | Cell-quality remediation (D50): the run uses other schemes than set because of poor cells. |
| `MESH003` | error | A mesh file that cannot be read: missing, of an unknown format, truncated or damaged (the reader says where); D202. |
| `MESH004` | warning | The operator check (D236): the least-squares gradient of a linear field is not exact on some cells (a degenerate or badly conditioned stencil); the cells are listed worst first. |
| `MESH005` | warning | The operator check (D236): with `turbulence.curvature_correction`, the velocity Hessian (DS/Dt) is not consistent for a quadratic field on some cells (error above `solver.operator_tolerance`), so the correction does not converge with the mesh there; the cells are listed worst first. |
| `MESH006` | warning | The operator check (D236): with `turbulence.qcr`, the velocity gradient is first order on some cells (error above `solver.operator_tolerance` for a quadratic field), so the divergence of the explicit QCR stress does not converge with the mesh there; the cells are listed worst first. |
| `RUN001` | note | Stopped by the `.stop` control file. |
| `RUN002` | note | Stopped by a signal (SIGTERM, SIGINT or SIGUSR2). |
| `RUN003` | warning | Stopped at the wall-time limit (`[run.limits] max_wall_time`). |
| `RUN004` | warning | Stopped at the step or iteration limit (`[run.limits] max_steps`, `max_iterations`). |
| `RUN005` | warning | Stopped at the output-volume limit (`[run.limits] max_output`). |
| `RUN006` | error | A steady solve that did not converge to the tolerance. |
| `RUN007` | error | A verification check of the case failed (order, error, reference or benchmark band, mass balance). |
| `RUN008` | warning | A change of the re-read case file that cannot be applied during a run. |
| `RUN009` | warning | A re-read case file that does not validate; the run goes on with the settings it had. |
| `RUN010` | warning | A restart file in another cell or face order: the continuation is not bit-identical. |
| `RUN011` | error | `--restart latest` found no restart file in the output directory. |
| `RUN012` | note | No benchmark of this hardware: the estimate has no wall time. |
| `RUN013` | error | An output file could not be written. |
| `RUN014` | error | The solution diverged (a value that is not finite); see `divergence.md` (D63). |
| `RUN015` | note | A setting changed by an entry of `[[schedule]]` (D128): key, old and new value. |
| `RUN016` | note | The diagnosis of a steady run at `[run.limits] diagnose_at` (default 2000 iterations): the trend of each residual, the verdict and advice (D198). |
| `RUN017` | warning | Stopped by the diagnosis: the steady run is not on its way to the tolerance within its iteration budget (`[run.limits] diagnosis = "stop"`). |
| `RUN018` | warning | Stopped at a study budget (`[run.limits] study_cell_hours`, `study_gpu_hours`, `study_output`). |
| `RUN019` | error | A study budget is already spent by the earlier runs of the study (`study_budget.jsonl`): the run does not start. |
| `RUN020` | warning | The estimate of the run (D81) does not fit into what is left of a study budget. |
| `RUN021` | warning | A steady iteration budget above 10 000 iterations: a steady case needs hundreds, at most a few thousand; look at the run after 2000 to 3000 iterations instead. |
| `RUN022` | note | Stopped through the live control port (`run.stop`): restart file, final output, exit code 3 (L19). |
| `RUN023` | note | A setting changed through the live control port (`run.set`): key, old and new value, with the equivalent `[[schedule]]` entry (L19). |
| `RUN024` | warning | A setting asked for through the live control port that cannot be applied (not valid for the case at that point); the run goes on with the settings it had. |
| `RUN025` | note | A setting changed by re-reading the case file (`.read_toml_N`, `[run] reread_every`): key, old and new value (D40, D214). |
| `DIV001` | note | Divergence cause: the Courant number near the failing cells is large; reduce the time step. |
| `DIV002` | note | Divergence cause: under-relaxation factors too large for the case. |
| `DIV003` | note | Divergence cause: poor cells (skewness, non-orthogonality, volume ratio) at the failing cells. |
| `DIV004` | note | Divergence cause: a boundary condition next to the failing cells. |
| `DIV005` | note | Divergence cause: the run failed within its first steps or iterations, from the initial state. |
| `DIV006` | note | Divergence cause: a linear solve that stopped at its iteration limit. |
| `BAL001` | warning | A relative mass imbalance of a step or iteration above `[balance] tolerance_step` (D89). |
| `BAL002` | warning | A relative mass imbalance over the run above `[balance] tolerance_run` (D89). |
| `BAL003` | warning | A numerical dissipation of kinetic energy below `[balance] tolerance_energy`: the discretisation adds energy (D90). |
| `STAT001` | warning | The averaging window of a statistic of `[output.statistics]` is too short: the batch means are correlated (lag-1 autocorrelation above 0.3), so the confidence interval of the mean is too narrow; run longer or use fewer batches (D71). |
| `LIM001` | note | A feature in use is graded partial or unverified on the verification status page (docs/verification/status.md); the message names its first gap. |
| `LIM002` | note | A feature in use ran outside its validated range (first-cell y⁺, Reynolds, Prandtl or Courant number, skewness), as `range` of docs/verification/features.toml gives it. |
| `LIM003` | note | A steady solve that stopped before its tolerance (iteration limit): the result may still change with more iterations. |
| `LIM004` | note | One mesh only: the discretisation uncertainty is not estimated (a level sweep gives the GCI, D67). |
| `LIM005` | note | The GPU backend ran on a device that is not in the list of verified devices (Tesla T4 so far). |
| `LIM006` | note | Automatic mesh remediation (D50, MESH002) changed the schemes on poor faces: first order there. |
| `LIM007` | note | A transient conjugate run with thermal acceleration of its solids (D213): the history before `physical_from` is not physical (the solids' time ran faster, their swing was too large); the end state and later cycle means are. |
| `LIM008` | note | An operator a feature in use needs is not exact on some cells (MESH005, MESH006, D236): the feature does not converge with the mesh there. |
| `LIM009` | note | A cycle-averaged conjugate run (D204): the solids carry the steady temperatures of the cycle-mean heat transfer (with `solid = "transient"` a history under it, marched in large steps), and each fluid phase ran with the wall temperatures held: the swing of the wall temperature over a cycle is in the result only as the estimate of `skin` or from `resolved_cycles`. |
| `LIM010` | note | A frozen-flow run (D118): the flow was held during the march, so its history is that of the settled flow at the temperatures of the last freeze or refresh; exact for a flow that does not depend on the temperature. |
| `TOOL001` | error | A tool id that the registry does not have (with the nearest id as a hint); `foehn commands` lists them. |
| `TOOL002` | error | A tool input that does not match the input schema of the tool (unknown field, wrong type, missing field). |
| `TOOL003` | warning | The journal of tool calls (D197) could not be written: the calls run, but the journal misses them. |
| `TOOL004` | error | A script to replay (`foehn run --script`) that cannot be read: not a file of JSON lines, or without calls. |
| `TOOL005` | warning | A replayed call (`foehn run --script`) that ended with another status or exit code than in its record; the replay stops there unless `--keep-going`. |
| `TOOL006` | error | A series of `result.page` that the run directory does not have, or a cursor that does not belong to it; the call without `series` lists the series. |
| `TOOL007` | error | A request to the live control port of a run (`run.status`, `run.stop`, `run.set`) that cannot be made or is refused: no run listens in the directory, the token is wrong, or the setting cannot be changed during a run (L19). |
| `MCP001` | error | `foehn mcp` started without `--read-root` and `--write-root`, or with an action table (`--policy`) inside a write root, where a tool could rewrite it. |
| `MCP002` | error | A run started through MCP without `[run.limits] max_wall_time` and `max_output`; the fix sets them with `set`. |
| `MCP003` | error | A call that the action table of the project (`foehn-mcp.toml`) refuses: start a run, exceed the estimate caps, write outside the write roots, overwrite results, publish. |
| `MCP004` | error | A call that the person declined when the server asked for approval. |
| `MCP005` | error | A call that needs the person's approval (`ask`), from a host that cannot ask (no MCP elicitation); the person allows it in `foehn-mcp.toml` or uses a host that asks. |
| `MCP006` | error | The action table `foehn-mcp.toml` cannot be read or written: not TOML, an unknown action kind or setting. |
| `STR001` | error | A structural case or model that is not valid: an unknown key or a wrong value, a point that is not a node, a beam without length, a missing material or section, non-positive properties. |
| `STR002` | error | The stiffness of a structural case is singular: the structure can move freely (a missing support or clamp, or a mechanism); the message names the first free unknown. For the modes of a free structure give a negative `shift`. |
| `STR003` | warning | A modal analysis may have missed modes: the Sturm count (the inertia of K − σM) finds more eigenvalues in the range than Lanczos returned, or Lanczos stopped above its tolerance. |
| `STR004` | error | The clamp fit (D290) cannot reach the measured frequency: it lies at or above the frequency of the rigidly clamped model (the model is too soft), or the mode does not depend on the spring. |
| `CIR001` | error | A circuit case (SF7b) that is not valid: an unknown key or a wrong value, a tube to a node that does not exist, a circuit without a pressure or drum node, a drum above the critical pressure, a CHF table that is not complete. |
| `CIR002` | error | The Newton solve of a circuit did not converge or stalled: give start flows (`start_flow`) near the expected ones, check that every part of the circuit reaches a pressure or drum node, and that the flow can carry the heat. |
| `CIR003` | warning | A state of a circuit left the range of the IAPWS-IF97 properties (the pressure fell below 1 kPa, or the enthalpy of a tube reached 4.1 MJ/kg): too much flow for the pressure level, or too little flow for the heat. |
| `CIR004` | warning | The critical heat flux ratio q_CHF/q of a tube segment is below the margin of the case (`chf_margin`, default 1.3); below 1 the wall is beyond the CHF and its temperature is that of vapour convection (Dougall-Rohsenow). |
| `GEN001` | error | An error without a more specific code. |

## Engine API (D86)

The tool registry (`foehncraft-tools`, above) sits on `foehncraft_flow::engine`,
and the command line and the Python module call only the registry
(`cargo xtask layering` checks it); a Rust
program (and the Python package, D61, through its native module) drives
a run the same way:

```rust
use foehncraft_flow::engine::{self, RunOptions, Roots};

let case = engine::load(path, &["fluid.viscosity=1e-3".into()], &Roots::default())?;
let opts = RunOptions::default();               // interrupts, build info, roots
let stop = opts.interrupts.clone();              // stop.request_stop(StopSignal::Api) from another thread
let out = engine::run(&case, "cavity", dir, &opts, |line| eprintln!("{line}"))?;
std::process::exit(out.exit_code());             // out.summary is the result.json object
```

| Item | Purpose |
|---|---|
| `load(path, sets, roots)` | Reads a case with overrides, inside the roots. |
| `check_file(path, sets, roots, dir, build)` | The check report (never fails: every problem is a diagnostic). |
| `estimate(case, dir, build)` | The estimate (D81). |
| `run(case, name, dir, opts, log)` | Runs; writes `result.json` and `run.jsonl`; returns the outcome (report and summary). |
| `Interrupts` | `request_stop`, `request_restart`, and the flags a signal handler sets. |
| `Observer::new(f)`, `RunOptions::observer` | `f` sees every row of the run record (`run.jsonl`) as it is written: the hook of the Python callbacks; it only reads, so results do not change. |
| `latest_restart(dir)` | The newest restart file. |
| `read_result(dir)` | `result.json` of an earlier run. |
| `diag::of_error(&e)` | The diagnostic of an error. |

## Tool registry (D196)

Every engine function a front end can call is a tool of the crate
`foehncraft-tools`, with a stable id (`area.verb`, dotted lowercase), a
contract version, a typed input and output whose JSON Schemas come from
the same structs, a category and an action kind (`read`, `write-files`,
`start-run`; D198). The catalogue is [tools.md](tools.md) (generated by
`foehn commands --markdown`; a test fails when it is out of date). The
clap commands build the typed input of their tool and call it, so their
text, `--json` objects and exit codes stay as above; `foehn call` runs any
tool by id from JSON and prints the envelope:

```json
{"format": "foehn.tool-result", "version": 1, "tool": "case.check",
 "tool_version": 1, "status": "finished", "exit_code": 0,
 "result": { ... }, "diagnostics": [], "limitations": []}
```

`result` is the object the clap command prints with `--json` (for
commands whose object has no `format`, the same fields without the
`foehncraft.command` wrapper). `status` and `exit_code` follow the exit
codes above. A call that fails has `result: null` and its diagnostic in
`diagnostics`: an unknown id is `TOOL001` (with the nearest id), an input
that does not match the schema `TOOL002`. `limitations` carries the
limitations of the tool's output (D199, above): those of the result
summary for `flow.run` and `result.read`, empty for a tool whose output
has none; like the summary's, they change neither `status` nor the exit
code. A breaking change of a tool's input or output raises its version.
Human text of a tool goes to standard error (`--quiet` silences it).

## Journal and replay (D197)

With `--journal FILE` on any command (or the environment variable
`FOEHN_JOURNAL`, which the Python package reads too), every tool call is
appended to the journal as one JSON line (`foehn.journal`): a header with
the version and commit, then per call `seq`, `date` (when it ended), `via`
(`cli`, `python`), `tool`, `tool_version`, `input` (the typed input: the
input after parsing, defaults filled in), `status`, `exit_code`,
`seconds`, and for a failed call its `diagnostic`. A call refused before it
ran (`TOOL001`, `TOOL002`) keeps its input as it came. Calls that a tool
makes itself (the calls of a replay) are not journaled again. A journal
that cannot be written is the warning `TOOL003`; the call stands.

A run started by a call writes that call (tool, version, typed input) into
its run record (`run.jsonl`) as the row `{"type": "call", ...}` after the
header, so the record says how the run was started, also without a
journal.

`foehn run --script FILE` (the tool `journal.replay`) runs the call rows of
a journal or a run record again, in order, in the working directory of the
replay and inside its roots (`--read-root`, `--write-root`). Refused calls
are skipped. Each call is compared with its record (status and exit code;
a row without them is expected to finish); the replay stops at the first
call that differs (`TOOL005`) unless `--keep-going`, and exits 1 when one
did. The output is `foehn.replay`: every call with its recorded and
replayed status, and the counts replayed, skipped and different. A script
that is not JSON lines or has no calls is `TOOL004`.

## MCP server (L18)

`foehn mcp` speaks the Model Context Protocol on standard input and
output: JSON-RPC 2.0, one message per line, nothing else on standard
output; log lines go to standard error (none with `--quiet`). It speaks
the protocol versions 2025-11-25, 2025-06-18, 2025-03-26 and 2024-11-05
and answers with the client's version when it has it, else 2025-06-18.
The guide for people is [agents.md](agents.md).

| Message | Answer |
|---|---|
| `initialize` | server `foehn` with its version, the capability `tools` and short instructions; the client's capability `elicitation` is noted |
| `tools/list` | every tool of the registry in catalogue order: `name` = the tool id (`--names underscore`: the first dot as `_`), `title`, `description` (summary and action kind), `inputSchema` (the input schema of [tools.md](tools.md)), `annotations` (`readOnlyHint` for tools of kind read, `destructiveHint` otherwise, `openWorldHint` false) |
| `tools/call` | the envelope `foehn.tool-result` of the call as `structuredContent` (protocol 2025-06-18 or newer) and as text after a first line with the status, exit code, limitations and the first diagnostic with its fix; `isError` when the status is `error` |
| `notifications/cancelled` | the running call is asked to stop (a run writes a restart file and the final output, exit code 3); a queued call is dropped |
| `ping`, `resources/list`, `prompts/list` | empty answers |
| server → client `notifications/progress` | for a call with a `progressToken`: the outer iteration or time step of a run, at most one a second |
| server → client `elicitation/create` | the approval question of a call (below); the answer `allow-once`, `allow-always`, `refuse-once` or `refuse-always` |

Calls that write or start a run run one after another on one worker.
Calls of the kinds `read` and `control` do not wait for it: they run at
once on a thread of their own, so a running `flow.run` can be watched
(`run.status`, `result.page`) and steered (`run.set`, `run.stop`, see
[the control port](#live-control-port-l19)). Before a call of a tool that
is not of kind read or control, the server asks the tool what it would write
(`effects`: the output folder or files, for a run its estimate, D81, and
its `[run.limits]`) and the action table (`--policy`, default
`foehn-mcp.toml`, outside every write root) what to do. The kinds are
`start-run`, `exceed-estimate` (the estimate over a cap of `[estimate]`:
`max_cells`, `max_host_memory`, `max_output`, `max_wall_time`),
`write-outside-root`, `overwrite-results` and `publish`, each `default`,
`allow`, `ask` or `refuse` (defaults: ask, ask, refuse, ask, refuse). An
"always" answer is written into the table; the table is read at every
call. A run without `[run.limits] max_wall_time` and `max_output` is
refused (`MCP002`) with the fix; a refused call (`MCP003`), a declined
question (`MCP004`) and a question the host cannot ask (`MCP005`) run
nothing. Every call is journaled with `"via": "mcp"` (D197); an approved
call carries `approval` (`kinds`, `by`: `person` or `action table`,
`answer`); `foehn run --script` skips the refused calls.

### `result.page`

The histories of a run, paged by cursor (D198 (e)): `{"dir": …}` lists
the series of the run folder (`available`: the CSV files by name without
`.csv`, and `record` for the rows of `run.jsonl`, each with its row count
and columns); `{"dir": …, "series": "residuals", "limit": 200}` returns
`foehn.result-page` with `columns`, `rows` (numbers where a cell is one),
`first`, `total` and `next_cursor` (`"residuals:200"`, none on the last
page), which the next call passes as `cursor`. For `record`, `types`
filters the rows (`iteration`, `table_row`, …). An unknown series or a
cursor of another series is `TOOL006`. Fields are never paged: their
files are named in the result summary.

## Live control port (L19)

A run started with `foehn flow --control-port`, or by any `flow.run`
through `foehn mcp`, listens on a TCP port of the loopback interface
(`127.0.0.1` only, an operating-system port) and writes two files into its
output directory: `control.json` (`foehn.control`: `host`, `port`, `pid`,
`token_file`) and `control.token` (128 random bits as hexadecimal text; on
Unix readable by its user only, mode 0600). Both are removed when the run
ends, so a directory without `control.json` has no run that can be
steered. The tools `run.status`, `run.stop` and `run.set` ([tools.md](tools.md))
find the port there; they need the directory inside the read roots
(`run.status`) or the write roots (`run.stop`, `run.set`).

The protocol is one JSON object per line, one request and one answer per
connection, the token in every request:

| Request | Answer |
|---|---|
| `{"token": …, "op": "status"}` | `foehn.run-status`: `state` (`running`, `finished`), `pid`, `wall_time`, `stop_requested`, `last` (the last `iteration` or `step` row of the run record: residuals, linear iterations), `control` (what the run did with the requests so far), `pending` |
| `{"token": …, "op": "stop"}` | the run stops at the end of the current iteration or step, as at `.stop`: restart file, final output, exit code 3, `RUN022` |
| `{"token": …, "op": "restart"}` | `restart.h5` at the end of the current iteration or step, the run goes on |
| `{"token": …, "op": "set", "set": {"key": value}}` | the settings are queued (`queued`); only the keys a re-read of the case file may change (D40), never `run.limits.*`: a run's limits are the person's, an assistant cannot raise its own |

Everything is applied after an iteration or step, like the control files
(D76), so a result never depends on when a request arrived within one. A
queued setting is checked against the case when it is applied: one that
does not fit (for example a relaxation below 1 in a transient case with one
outer iteration) is `RUN024` and left out, the run goes on with what it
had. Each outcome is an `event` row of the run record (`RUN023`, `RUN024`)
and an entry of `control` in the status, with the equivalent
`[[schedule]]` entry (`at_iteration` or `at_step` = the next iteration or
step), so that a person can reproduce a steered run from its case file.
A wrong or missing token, a bad request and a run that is gone are
`TOOL007`.

The live calls are in the journal like every call; `foehn run --script`
skips them (they act on a process that is gone), and a replayed `flow.run`
runs to its end, not to the step where an assistant stopped it: a stopped
run ends with exit code 3, so the replay reports that call as different
unless the run record's `RUN022` stop is reproduced with a `max_steps`
limit or a schedule.

## Python package (D61)

`python/foehncraft` on the same engine (native module
`foehncraft._native`, `crates/foehncraft-py`), or on `foehn --json`
(`backend="cli"`). The native functions exchange the JSON formats of the
command line: `check(path, sets, out)` (`foehncraft.check`),
`run(path, sets, out, observer, log)` (`foehncraft.result`),
`material_database(name, temperature, pressure)` (`foehncraft.materials`),
`case_schema()` (D77), `codes()` (D78), `resolved(path, sets)` (the
resolved case as TOML), `expression(source, env, allowed, points, base)`
(values of an expression); errors are `FoehncraftError(code, message)`.
Each of these functions is a one-line wrapper of its tool of the
registry (D196); `call(id, input, observer, log)` runs any tool by id
from JSON text and returns the envelope, `commands()` the catalogue. In
the package: `foehncraft.call("case.check", case=path)` and
`foehncraft.commands()` (through `foehn call` without the native module).
The Python API is in [recipe 33](recipes/python.md#python-api).

## C-free core (D87)

`foehncraft-flow` without default features (no HDF5) depends on no C
library: `python3 -m tools.check.c_free` checks the dependency tree in CI. HDF5
(output and restart) and METIS (partitioning) are optional.
