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
other output file with its SHA-256 hash. HDF5 files are written
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
| `ENE006` | error | An energy setting this version does not have: a passive scalar named `temperature` beside `[energy]`, solid regions with a sliding mesh, a restart, passive scalars or a periodic mesh. |
| `ENE008` | error | An `[[energy.sources]]` entry for a region the case does not have, or without exactly one of `power` and `power_density`. |
| `CHT004` | error | An `[interfaces.<name>]` table that is not on a region interface (regions unknown, the same, two fluids, or not touching), or without exactly one of `contact_resistance` (≥ 0) and `thin_wall` (thickness > 0 with `conductivity` or `material`). |
| `ENE007` | warning | A steady run that converged while its heat balance does not close (imbalance above 1 % of the largest heat flow): the temperature residual is small but the energy solve is not converged; tighten `[solver.energy] tolerance` and raise `max_iter`. |
| `THM001` | error | A steady run whose boundaries leave the temperature level undetermined: no inlet and no `temperature` or convection boundary, and no net heat input. |
| `THM002` | error | A steady run with a net heat input (heat flux, heat flow, source or viscous heating) and no way for the heat to leave: no inlet and no `temperature` or convection boundary. |
| `THM003` | error | A thermal condition on a boundary that cannot carry it: a flux or convection at an inlet, any condition at an outlet, a periodic boundary, a symmetry plane or a wedge. |
| `THM004` | error | A thermal value out of range: a temperature or ambient temperature ≤ 0 K, a heat-transfer coefficient ≤ 0, a heat flow that is not finite. |
| `THM005` | warning | Walls without a `[boundary.<name>.energy]` table in a run with the energy equation: they are adiabatic (listed once). |
| `THM009` | error | A `[boundary.<name>.energy]` table without `[energy]`: the energy equation is solved only with an `[energy]` table. |
| `THM010` | error | A `[boundary.<name>.energy]` table for an unknown boundary, without exactly one condition (`temperature`, `heat_flux`, `heat_flow`, or `heat_transfer_coefficient` with `ambient_temperature`), or with `heat_transfer_coefficient` and `ambient_temperature` apart. |
| `THM011` | error | An inlet without `[boundary.<name>.energy] temperature` in a run with the energy equation. |
| `REG001` | error | A region that is not a fluid where it cannot be: a solid region without `[energy]` (solids are solved by the energy equation, stage 7d1), a porous zone (not in this version), or a solid material as the `[fluid]`. |
| `REG002` | error | A region that names a cell zone the mesh does not have, or a zone in two regions. |
| `MRF001` | error | An `[[mrf]]` table and a `[[sliding]]` table that turn the same cells (one of them without a `zone`, or both with the same one): a zone is turned by a moving reference frame (steady or transient, D185) or by a sliding mesh (D182), not by both. |
| `MRF002` | error | An `[[mrf]]` table with no or two angular velocities (`omega`, `rpm`), a zero axis, a cell zone the mesh does not have (with the nearest name as a hint), no `zone` where there are several tables, or a zone that two tables rotate. |
| `MRF003` | warning | The interface or a boundary of a rotating zone is not a surface of revolution about its axis: the relative flux of the cells next to it does not close, and the solution there is not that of a frozen rotor (cell count and worst closure error). |
| `SLD001` | error | A `[[sliding]]` table in a steady case: the sliding mesh moves in time (D182); a steady case with a rotating zone uses `[[mrf]]` (D55). |
| `SLD002` | error | A `[[sliding]]` table with no or two angular velocities (`omega`, `rpm`), a zero speed or axis, a `shift` of 0 or of a whole turn, a cell zone the mesh does not have (with the nearest name as a hint), or a second table (one rotor per case so far); with the interface imprinted: `rotor_interface` without `stator_interface`, no `zone`, or a `shift` (the time step is free). |
| `SLD003` | error | The time stepping does not fit the sliding mesh: the time step is not `shift` × 2π/(sectors × Ω) (the right value is given), the step is adaptive or error controlled, `end_time` is not a whole number of steps, or `[verification]` asks for a sweep over levels or time steps; with the interface imprinted the scheme must be fully implicit (`extrapolation = false`, at least 2 `outer_iterations`). |
| `SLD004` | error | The rotor cannot be moved by whole sectors: its cells have no image when it is turned by one sector (it is not a ring of identical sectors about the axis), or a boundary of the rotor is not a wall, slip wall, symmetry plane or wedge. |
| `SLD005` | error | The interface of a sliding mesh by imprinting cannot be imprinted: a boundary or cell zone of `[[sliding]]` is missing or wrong, the two interface boundaries are not on one circular cylinder about the axis or do not cover the same surface, rotor and stator cells share a face or a node, a `[boundary]` table names an interface boundary, or a restart is asked for (not yet available). |
| `SLD006` | warning | The cells next to the imprinted interface of a sliding mesh do not close (|Σ S|/Σ |S| above 1e-8): the faces of the interface are not those of a mesh extruded along the axis, and the mesh flux of the rigid turn leaves a small mass source. |
| `ROT001` | warning | A wall of a rotor does not turn with its zone (a stationary shroud in a rotating zone): its torque is in that of the rotor, but the shaft power −Ω·T is not its work on the fluid; give the walls of the rotor in `walls` of the `[[mrf]]` or `[[sliding]]` table. |
| `ROT002` | error | A `walls` entry of `[[mrf]]` or `[[sliding]]` that is not a wall of the case, or an `[[output.machines]]` table without a turning zone, with a `rotor` that is no turning zone (or none where there are several), without an inlet or an outlet, with a boundary the case does not have, or with a name used twice. |
| `MESH001` | warning | Faces whose skewness exceeds `solver.skewness_limit` (poor cells, typically tetrahedra of a mesh that was not optimised). |
| `MESH002` | note | Cell-quality remediation (D50): the run uses other schemes than set because of poor cells. |
| `MESH003` | error | A mesh file that cannot be read: missing, of an unknown format, truncated or damaged (the reader says where); D202. |
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
| `TOOL001` | error | A tool id that the registry does not have (with the nearest id as a hint); `foehn commands` lists them. |
| `TOOL002` | error | A tool input that does not match the input schema of the tool (unknown field, wrong type, missing field). |
| `TOOL003` | warning | The journal of tool calls (D197) could not be written: the calls run, but the journal misses them. |
| `TOOL004` | error | A script to replay (`foehn run --script`) that cannot be read: not a file of JSON lines, or without calls. |
| `TOOL005` | warning | A replayed call (`foehn run --script`) that ended with another status or exit code than in its record; the replay stops there unless `--keep-going`. |
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
limitations of a run (D199, below) once the run tools of D196 part 2
fill it from the result summary. A breaking change of a tool's input or output raises its version.
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
