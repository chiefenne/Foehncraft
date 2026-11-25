# Tool inputs (`examples/cases/tools/`)

Part of the [Case catalogue](../cases.md): JSON inputs of `foehn call` (the tool registry, D196).

Each file is the input of one tool of [tools.md](../tools.md), with paths
relative to the repository root. Run with
`foehn call <tool> --input examples/cases/tools/<file>`; the result
envelope is on standard output ([gallery 171](../gallery/runs.md#tool-registry),
[recipe 171](../recipes/runs.md#tools)).

| File | Tool | Shows |
|---|---|---|
| `check_cavity.json` | `case.check` | the check of a flow case with a `set` override: findings and the estimate |
| `run_cavity.json` | `flow.run` | the lid-driven cavity at Re = 100 run into `output/tools/cavity` (439 iterations, about 2 s) |
| `read_cavity.json` | `result.read` | the result summary of that run |
| `derive_cavity.json` | `post.derive` | vorticity and the Q-criterion from its restart file |
| `water_at_330.json` | `materials.show` | water from the material database at 330 K |
| `ramp_expression.json` | `expression.eval` | a ramped parabolic inlet profile with constants, at three points |
| `scalar_mms.json` | `study.scalar` | the manufactured-solution sweep of scalar transport on polygons (observed order about 2) |
| `cavity_session.jsonl` | `journal.replay` | a script of four calls (check, run, derive, read) for `foehn run --script` (D197); not an input of `foehn call` but its `script` |
| `equilibrium_methane_air.json` | `chemistry.equilibrium` | the adiabatic flame state of stoichiometric methane and air from 298.15 K with GRI-Mech 3.0 (T = 2224.9 K, SF2) |
| `reactor_methane_air.json` | `chemistry.reactor` | the constant-pressure reactor of `reactor_methane_air_gri30.toml`: ignition after 3.4 ms and the end state at the equilibrium, the trace in `output/tools/reactor/reactor.csv` |
| `page_cavity.json` | `result.page` | the first 100 rows of the residual history of that run; the answer's `next_cursor` gives the next page |
| `structure_reed_beam.json` | `structure.run` | the reed valve beam with a rigid clamp: static deflection and eight modes into `output/reed_beam` ([gallery 295](../gallery/structure.md#reed-beam)) |
| `structure_reed_beam_clamp_springs.json` | `structure.run` | the same reed on clamp springs (D290) into `output/reed_beam_clamp_springs` |
| `structure_reed_plate.json` | `structure.run` | the reed petal as MITC4 shells on clamp springs along its edge into `output/reed_plate` ([gallery 297](../gallery/structure.md#reed-plate)) |
| `structure_reed_plate_pulse.json` | `structure.run` | the reed petal struck by a pressure pulse, in time, into `output/reed_plate_pulse` ([gallery 299](../gallery/structure.md#reed-dynamics)) |
| `structure_reed_plate_hot.json` | `structure.run` | the hot reed petal: thermal bending into `output/reed_plate_hot` |
| `structure_reed_stop.json` | `structure.run` | the curved reed stop as twenty-node hexahedra read from Gmsh into `output/reed_stop` ([gallery 300](../gallery/structure.md#reed-stop)) |
| `structure_reed_valve_stop.json` | `structure.run` | the reed petal between its seat and its stop, lifted by a pressure pulse, into `output/reed_valve_stop` ([gallery 302](../gallery/structure.md#reed-valve-contact)) |
| `structure_fit_clamp.json` | `structure.fitclamp` | the rotational clamp spring that puts the reed's first mode at a measured 150 Hz |

The folder `mcp/` holds what an assistant session through `foehn mcp`
(L18) needs ([gallery 175](../gallery/runs.md#mcp),
[recipe 175](../recipes/runs.md#mcp), [agents.md](../agents.md)):

| File | Shows |
|---|---|
| `mcp/foehn-mcp.toml` | the action table of a project with the default of every kind written out and the estimate caps (D198 (b)) |
| `mcp/heat_sink_session.jsonl` | the calls of a scripted assistant on the finned heat sink: check, a run refused for want of `[run.limits]`, the run with the fix approved once, the series and the residual history; `python3 -m tools.plot.mcp_session` sends them to `foehn mcp` |

The folder `control/` holds the inputs of the live control of a running case
(L19, [gallery 177](../gallery/runs.md#control-port),
[recipe 177](../recipes/runs.md#control-port)). They name the output folder
of the run of the recipe, `output/tools/control/live`:

| File | Tool | Shows |
|---|---|---|
| `control/run_status.json` | `run.status` | the live state of the running cavity: the last step with its residuals, the wall time, what the run did with earlier changes |
| `control/run_set.json` | `run.set` | the time step of the transient cavity changed from 0.005 s to 0.01 s while it runs |
| `control/run_set_refused.json` | `run.set` | a relaxation below 1 that a transient case with one outer iteration refuses (`RUN024`): shown in the status, the run goes on |
| `control/run_stop.json` | `run.stop` | stops the run at the end of the current step: restart file, final output, exit code 3 (`RUN022`) |

The folder `mcp/heat_sink_study/` holds the recorded study of an assistant
([gallery 181](../gallery/runs.md#heat-sink-study),
[recipe 181](../recipes/runs.md#heat-sink-study),
[agents.md](../agents.md#a-recorded-study-the-fin-pitch-of-a-heat-sink)):

| File | Shows |
|---|---|
| `mcp/heat_sink_study/foehn-mcp.toml` | the action table of the study: start runs allowed, everything else refused, a cap of 40 000 cells |
| `mcp/heat_sink_study/journal.jsonl` | the session as `foehn mcp` journaled it: two checks, three refused runs, two live status calls, four sweep runs, the reads; `foehn run --script` replays it |
| `mcp/heat_sink_study/results.csv` | what the assistant read from each run: iterations, the mean and maximum temperature of the base top, the pressure drop, the cell-hours spent |
| `mcp/heat_sink_study/session.csv` | the 27 calls of the session with start, duration, outcome and refusal code; `tools/plot/heat_sink_agent.py` draws the timeline from it |
