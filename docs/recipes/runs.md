# Run control and health

Part of the [Recipes](../recipes.md): restart, stop and steer runs, drive runs from a script or an agent, balances and divergence reports.

| # | Recipe |
|---|---|
| 22 | [Restart a run, stop it and steer it while it runs](#restart) |
| 26 | [Drive runs from a script or an agent](#run-contract) |
| 27 | [Check the balances and read a divergence report](#health) |
| 37 | [Read the log and upgrade a case file](#log) |
| 38 | [Choose the columns of the log table](#log-columns) |
| 39 | [Prepare changes of settings before the run](#schedule) |
| 166 | [Set an iteration budget and a budget per study](#budgets) |
| 171 | [Call any tool from a script or an agent](#tools) |
| 172 | [Record a session and replay it](#journal) |
| 175 | [Connect an assistant through the MCP server](#mcp) |
| 177 | [Watch and steer a running case](#control-port) |
| 181 | [Let an assistant run a parameter study](#heat-sink-study) |

---

<a id="restart"></a>
## 22. Restart a run, stop it and steer it while it runs

Every run writes `restart.h5` at its end (`run.restart_final`), and on a
schedule if asked. A run continued from it gives bitwise the result of
the run without the interruption (same mesh, same build):

```toml
[run]
restart_every = 500          # restart.h5 every 500 steps
restart_at = [2000, 4000]    # also restart_002000.h5, restart_004000.h5
```

```
target/release/foehn flow my_case.toml
# continue to a later end time, into another directory:
target/release/foehn flow my_case.toml --restart output/my_case/restart.h5 --end-time 20 --out output/my_case_2
```

or, in the case file, `[run] restart = "output/my_case/restart.h5"` (a
path in the case file is relative to the case file, on the command line
to the current directory).

While a case runs, files in its output directory steer it (the solver
never deletes them; remove them before the next run):

```
touch output/my_case/.stop            # write restart.h5 and the results, then end
touch output/my_case/.result_1500     # write fields_001500.vtkhdf at step 1500
touch output/my_case/.read_toml_1200  # re-read my_case.toml at step 1200
```

A re-read applies the output and restart schedule, the stop criteria, the
relaxation, the linear-solver tolerances and the time-step settings; the
log lists each change (`key: old -> new`) and each change it cannot apply
during a run. `run.reread_every = N` re-reads every N steps without a
control file. An intermediate result only looks at the run: with or
without `.result_N` files (and with any `[balance] every`) the residuals
and the final result are the same bit by bit, so a history of a
monitored value taken from such files is the history of the plain run.

[↑ back to top](#run-control-and-health)

---

<a id="run-contract"></a>
## 26. Drive runs from a script or an agent

`foehn` keeps a contract for programs that run it (D62, D64, D78; the
details are in [interfaces.md](../interfaces.md)): one JSON object on
standard output with `--json`, an exit code that says how the run ended,
codes for every finding, and a result file without timings. From the
repository root, after `cargo build --release`:

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml > result.json; echo "exit $?"
python3 -c "import json; r = json.load(open('result.json')); print(r['status'], r['solves'][0]['iterations'])"
```

| Exit code | Status | Next step of a script |
|---|---|---|
| 0 | `finished` | use the result |
| 1 | `error` | read `diagnostics`: code, key, message and a fix where one exists |
| 2 | `failed-checks` | the run ended, a check failed (`RUN006` not converged, `RUN007` a check of the case) |
| 3 | `stopped` | a signal, `.stop` or a limit: continue with `--restart latest` |

**Change keys without editing the file**: `--set KEY=VALUE`, repeatable,
the value in TOML (a bare word is a string). The overrides are in the
result and survive a re-read of the case file.

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --set fluid.viscosity=0.005 --set solver.momentum.relaxation=0.8 --out output/re200
```

**Check before running** (D81): the findings with codes and an estimate
of memory, output and, where this hardware has a benchmark, wall time:

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --check
```

A finding with a fix gives the `--set` arguments that resolve it:

```
target/release/foehn --json flow examples/cases/flow/taylor_green_cfl10.toml --check --set time.max_growth=3
```

**Budgets** (D83): `[run.limits]` (`max_wall_time` in s,
`max_iterations`, `max_steps`, `max_output` in bytes) stops a run with its
restart file and results and exit code 3. Continue it with
`--restart latest`, bitwise as if it had not stopped:

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --set run.limits.max_iterations=300 --out output/budget; echo "exit $?"
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --restart latest --out output/budget; echo "exit $?"
```

**Signals** (unix): SIGUSR1 writes `restart.h5` and the run goes on;
SIGTERM, SIGINT and SIGUSR2 stop it cleanly (exit 3); a second Ctrl-C
ends at once. A batch job that is sent SIGTERM before its end is
continued by the next job with `--restart latest`:

```
target/release/foehn --quiet flow examples/cases/flow/cavity_re100_transient.toml --out output/job & pid=$!
sleep 5; kill -USR1 $pid      # a restart file now
sleep 5; kill -TERM $pid; wait $pid; echo "exit $?"
target/release/foehn --quiet flow examples/cases/flow/cavity_re100_transient.toml --restart latest --out output/job
```

**Follow a run**: `run.jsonl` gets one JSON object per iteration or step,
flushed at every log line (`tail -f output/job/run.jsonl`).

**Confine the files** (D82): `--read-root DIR` and `--write-root DIR`
refuse every path outside (`CASE007`), for agents that must not touch the
rest of the machine:

```
target/release/foehn --json flow examples/cases/flow/cavity_re100_budget.toml --read-root examples --write-root output
```

**Compare results**: two identical runs give byte-identical outputs
(D84); `result.json` lists their BLAKE3 hashes, and
`provenance.case_hash` with `provenance.commit` says whether two results
are comparable. The keys and their units: `foehn schema` (JSON Schema)
and [keys.md](../keys.md).

[↑ back to top](#run-control-and-health)

---

<a id="health"></a>
## 27. Check the balances and read a divergence report

Every flow run books mass and kinetic energy in a ledger (D88, D89;
details in [interfaces.md](../interfaces.md#conservation-ledger-d88-d89-d90-part-4)).
From the repository root, after `cargo build --release`:

```
target/release/foehn flow examples/cases/flow/cavity_re100_energy.toml --out output/energy
```

The log shows the full table every 100 iterations (`[balance] every`) and
at the end; for the steady cavity the lid power goes into the viscous
dissipation and, about 1.5 %, into the numerical dissipation, which the
second law wants positive (D90 part 4); its two parts, the pressure part
of the momentum interpolation and the part of the convection scheme,
follow in the table. The same
table is in `flow.md`, the history in `balance.csv`, the bars in
`balance.svg`, the values under `solves[].balance` in `result.json`.

A transient run logs one compact balance line with every progress line
and a ledger of the whole run at the end; draw the time histories with

```
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml --set "verification.time_steps=[]" --out output/tg
python3 -m tools.plot.balance_chart output/tg
```

**Tolerances**: `[balance] tolerance_step` (1e-8, the relative mass
imbalance of a step or of the final steady state), `tolerance_run` (1e-4,
the mass over the run) and `tolerance_energy` (1e-2, the negative
numerical dissipation relative to the largest term) give the warnings
`BAL001`–`BAL003`; `action = "stop"` ends a transient run at the
step with exit code 3 instead (here with a tolerance no run meets: the
Taylor–Green vortex has inflow and outflow through every side, and its
mass imbalance is at round-off, about 1e-17):

```
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml --set "verification.time_steps=[]" --set balance.tolerance_step=1e-30 --set 'balance.action="stop"' --out output/stop; echo "exit $?"
```

**The energy budget over the mesh** (the gallery chart):

```
python3 -m tools.plot.energy_budget_chart output/energy_budget.svg
```

**A divergence** (D63) ends the run with exit code 1 and `RUN014`, after
it wrote the last valid state (`fields_last_valid.vtkhdf`, `restart.h5`),
the failing cells with three layers of neighbours (`divergence.vtkhdf`,
threshold on `divergence_layer` in ParaView) and the report
`divergence.md` / `divergence.json` with the likely causes ranked and what
to change. Provoke one: SIMPLE with no relaxation

```
target/release/foehn flow examples/cases/flow/cavity_re1000.toml --set 'solver.algorithm="simple"' --set solver.momentum.relaxation=0.99 --set solver.pressure.relaxation=1.0 --out output/diverge; echo "exit $?"
cat output/diverge/divergence.md
```

The first cause is `DIV002` (relaxation) with the fix
`--set solver.momentum.relaxation=0.69 --set solver.pressure.relaxation=0.31`,
which runs. A time step far too large for a transient case points to
`DIV001` with a smaller time step:

```
target/release/foehn flow examples/cases/flow/cavity_re100_transient.toml --set time.adaptive=false --set time.time_step=5.0 --set solver.correctors=1 --set 'time.scheme="crank-nicolson"' --out output/diverge_t; echo "exit $?"
```

With `--json`, the `RUN014` diagnostic in `result.json` carries the fix
of the most likely cause as `--set` arguments, so a script can retry.

[↑ back to top](#run-control-and-health)

<a id="log"></a>
## 37. Read the log and upgrade a case file

The log (gallery 27, D125) starts with the version and the build, then the blocks that show what was set; read it top down: **Input and output files** (the case file and every `--set`), **Run control**, **Mesh information** (counts and quality), **Fluid properties**, **Solution settings** (one column per equation), **Boundary conditions**. The table follows (every `solver.log_every` iterations or steps; the balance lines follow `balance.every`), then the summaries and the closing frame with the exit code.

Frame style and the case-file version:

```
target/release/foehn flow examples/cases/flow/cavity_re100.toml --set 'log.style="box"'
target/release/foehn case upgrade examples/cases/flow/old_case.toml          # list the changes
target/release/foehn case upgrade --write examples/cases/flow/old_case.toml  # rewrite, comments kept
target/release/foehn case upgrade --check examples/cases/flow/*.toml         # exit 1 if one is old
```

A case file without `case_format` is read as 0.5: its flat solver keys (`solver.relaxation_velocity`, `solver.momentum_tolerance`, ...) map to the tables `[solver.momentum]` and `[solver.pressure]`, and `patch` maps to `boundary`; each mapping is listed in the log. Scripts that read `patches.csv` or the field `patches` of `result.json` read `boundaries.csv` and `boundaries` now.

Keys: `case_format`, `log.style`, `solver.log_every`, `balance.every`.

[↑ back to top](#run-control-and-health)

---

<a id="log-columns"></a>
## 38. Choose the columns of the log table

Which residuals, which monitoring values and which extremes the table shows is a list in the case file (gallery 28):

```toml
[[output.probes]]
name = "vortex"
point = [0.08, 0.08, 0.0078125]

[log]
digits = 3                       # significant digits of the numbers
columns = [
  "residual.momentum",           # ResU ResV ResW
  "residual.continuity",         # ResM
  "probe.vortex.velocity_x",     # a named probe
  "cell.2000.pressure",          # the value in one cell (debugging)
  "max.velocity_magnitude",      # the largest value and its cell
  "min.pressure",                # the lowest value and its cell
  "boundary.ymax.mass_flow",     # mass flow out of a boundary (kg/s)
  "time.per_iteration",          # wall-clock seconds per iteration
]
```

Names are `<kind>.<target>.<quantity>`; the quantities are `velocity_x`, `velocity_y`, `velocity_z`, `velocity_magnitude` and `pressure`. Linear iterations: `linear.momentum`, `linear.pressure`. Transient runs add `time.simulated`, `time.dt`, `time.courant`, `time.outer` and `time.error`. A wrong name stops `--check` with the list of valid names (code `CASE003`); a cell that the mesh does not have stops the run before the first iteration. Probe values are reconstructed with the cell gradients of the last iteration, as the probes of a transient run; they need the fields from the device at every table line, which costs a download each. Without `[log] columns` the default table is shown.

```
target/release/foehn flow examples/cases/flow/cavity_re100_columns.toml
target/release/foehn flow examples/cases/flow/cavity_re100.toml --set 'log.columns=["residual.momentum","max.velocity_magnitude"]'
target/release/foehn schema --markdown | grep '^| `log\.'
```

Keys: `log.columns`, `log.digits`, `log.header_every`, `solver.log_every`, `output.probes`.

[↑ back to top](#run-control-and-health)

---

<a id="schedule"></a>
## 39. Prepare changes of settings before the run

Start with cautious settings and let them change at fixed points, instead of editing the case file during the run (gallery 29):

```toml
[solver.momentum]
relaxation = 0.5

[[schedule]]
at_iteration = 100                 # steady: iteration; transient: at_step or at_time
set = { solver.momentum.relaxation = 0.9, solver.tolerance = 1e-9 }

[[schedule]]
at_time = 0.05                     # transient: the run lands a step on this time
set = { time.time_step = 0.005 }
```

(One case uses one kind of trigger: the second entry belongs to a transient case.) The value of an entry holds from its point on, until a later entry sets the same key. List the entries in ascending order. Keys are dotted paths of the case file; the table `set = { solver.momentum.relaxation = 0.9 }` and the quoted form `{ "solver.momentum.relaxation" = 0.9 }` are the same. `foehn flow --check` reads the whole schedule and stops with the entry and the reason when a trigger, a key or a value is wrong (code `SCHED001`): a key that cannot change during a run (mesh, boundaries, schemes, fluid) names the key, a misspelled key the nearest one, a value that does not validate (a relaxation above 1) the entry.

During the run each change is a ` [I]` line and an `event` row `RUN015` of `run.jsonl`. A re-read of the case file (`.read_toml_N`, `[run] reread_every`) keeps the entries that have fired in force over the values of the file; the schedule itself cannot be re-read: change it by stopping the run and restarting (`--restart latest` applies the entries that lie before the restart point at once). Without `[[schedule]]` the run, its files and its log are bitwise those of a run without the table.

```
target/release/foehn flow examples/cases/flow/cavity_re100_schedule.toml
target/release/foehn flow examples/cases/flow/cavity_re100_schedule.toml --check
grep schedule output/cavity_re100_schedule/flow.log
```

Keys: `schedule[].at_iteration`, `schedule[].at_step`, `schedule[].at_time`, `schedule[].set`.

[↑ back to top](#run-control-and-health)

---

<a id="budgets"></a>
## 166. Set an iteration budget and a budget per study

The iteration budget of a steady case is `solver.max_iterations` (or a
smaller `[run.limits] max_iterations`); keep it at about 3000, a steady case
needs hundreds, at most a few thousand (a budget above 10 000 draws
`RUN021` in the check). Look at the diagnosis instead of running longer
(gallery 166):

```toml
[solver]
max_iterations = 3000

[run.limits]
diagnose_at = 2000        # the default; 0 switches the diagnosis off
diagnosis = "stop"        # default "record": write it and go on
```

```
target/release/foehn flow case.toml --check        # "iteration budget 3000 (solver.max_iterations); diagnosis at iteration 2000 (stop)"
target/release/foehn flow case.toml
grep diagnosis output/case/flow.log
```

Read the verdict: `on-track` needs nothing; `slow` falls but would reach
the tolerance only after the budget (the residual is no stop measure: if
the monitored quantities have settled, the run is done, otherwise ramp the
relaxation up with `[[schedule]]` or use the coupled algorithm); `stalled`
and `oscillating` point at the relaxation, the mesh or an unsteady flow;
`rising` asks for a lower relaxation before the run diverges. A run stopped
by the diagnosis continues with `--restart latest` after the change.

For a study, give each variant the same budgets; the ledger is the
`study_budget.jsonl` next to the variant directories:

```python
import foehncraft as fc
from foehncraft.study import Study, grid

base = fc.Case.load("examples/cases/flow/cavity_re100.toml")
base = base.with_values({"run.limits.study_cell_hours": 50.0, "run.limits.study_output": 2e9})
rows = Study(base, grid({"fluid.viscosity": [0.01, 0.005, 0.002]}), "output/study").run()
```

Each run shows its estimate and what the study has left in the `budget`
row of `run.jsonl` and in the log; a variant that spends the rest stops
with exit code 3 (`RUN018`), the ones after it do not start (`RUN019`), and
`foehn flow --check` reports a spent study before anything runs. To raise a
budget, raise the key (it may also change during a run through
`.read_toml_N`); to start the study's account again, delete
`study_budget.jsonl`.

Keys: `run.limits.diagnose_at`, `run.limits.diagnosis`,
`run.limits.study_cell_hours`, `run.limits.study_gpu_hours`,
`run.limits.study_output`.

[↑ back to top](#run-control-and-health)

---

<a id="tools"></a>
## 171. Call any tool from a script or an agent

Every function of `foehn` is a tool with an id (D196); [tools.md](../tools.md)
lists them with their input fields. From the repository root, after
`cargo build --release`:

```
target/release/foehn commands > catalogue.json                     # ids, schemas, action kinds
target/release/foehn call system.info                              # no input: no fields
target/release/foehn call case.check --input examples/cases/tools/check_cavity.json
echo '{"name": "water", "temperature": 330}' | target/release/foehn call materials.show --input -
target/release/foehn call flow.run --input examples/cases/tools/run_cavity.json; echo "exit $?"
target/release/foehn call result.read --input examples/cases/tools/read_cavity.json
target/release/foehn call post.derive --input examples/cases/tools/derive_cavity.json
target/release/foehn call expression.eval --input examples/cases/tools/ramp_expression.json
```

The envelope on standard output has `status` and `exit_code` (the exit
codes of [recipe 26](#run-contract)), `result` (the tool's output, the
object of the command's `--json`) and `diagnostics` (the error of a call
that failed: `TOOL001` an unknown id, `TOOL002` an input that does not
match the schema, or the code of the tool's error). Human text goes to
standard error; `--quiet` silences it. `--read-root` and `--write-root`
confine the files of the call (D82).

From Python, with the native module in process, else through `foehn call`:

```python
import foehncraft as fc

cat = fc.commands()
[t["id"] for t in cat["tools"] if t["action"] == "read"]       # what an agent may call without asking
r = fc.call("case.estimate", case="examples/cases/flow/cavity_re100.toml")
r["result"]["cells"], r["result"]["host_memory"]
r = fc.call("flow.run", case="examples/cases/flow/cavity_re100.toml",
            out="output/tools/cavity", observer=lambda row: False)   # sees every row of run.jsonl
```

A new engine function becomes a tool in `crates/foehncraft-tools` (a type
implementing `Tool`, one line in `registry()`), then gets its clap command
and Python wrapper; regenerate the catalogue with
`target/release/foehn commands --markdown > docs/tools.md`, and run
`cargo xtask layering`.

[↑ back to top](#run-control-and-health)

---

<a id="journal"></a>
## 172. Record a session and replay it

Every tool call can go into a journal (D197). From the repository root,
after `cargo build --release`:

```
export FOEHN_JOURNAL=output/tools/session.jsonl                 # or --journal FILE on each command
target/release/foehn flow examples/cases/flow/cavity_re100.toml --check
target/release/foehn flow examples/cases/flow/cavity_re100.toml --out output/tools/cavity
target/release/foehn post derive examples/cases/flow/cavity_re100.toml --out output/tools/cavity --fields vorticity
unset FOEHN_JOURNAL
target/release/foehn run --script output/tools/session.jsonl; echo "exit $?"   # 0: every call as recorded
```

Each line of the journal is one call: `tool`, `input` (typed, every
default filled in), `status`, `exit_code`, `seconds` and `via` (`cli`,
`python`). A run started by a call also holds that call in its run record,
so one run can be replayed from its output directory:

```
target/release/foehn run --script output/tools/cavity/run.jsonl
```

A script can be written by hand: one JSON object per line with
`"type": "call"`, `tool` and `input`, as in
`examples/cases/tools/cavity_session.jsonl`; without a recorded status a
call is expected to finish.

```
target/release/foehn run --script examples/cases/tools/cavity_session.jsonl
target/release/foehn --json run --script examples/cases/tools/cavity_session.jsonl --keep-going   # foehn.replay
```

The replay runs in the working directory and inside `--read-root` and
`--write-root` (D82); it stops at the first call that ends otherwise than
recorded (`TOOL005`, exit 1) unless `--keep-going`. From Python:

```python
import os
import foehncraft as fc

os.environ["FOEHN_JOURNAL"] = "output/tools/py_session.jsonl"
fc.call("case.check", case="examples/cases/flow/cavity_re100.toml")
del os.environ["FOEHN_JOURNAL"]
fc.call("journal.replay", script="output/tools/py_session.jsonl")["result"]["different"]   # 0
```

[↑ back to top](#run-control-and-health)

---

<a id="mcp"></a>
## 175. Connect an assistant through the MCP server

`foehn mcp` serves every tool to an assistant over MCP on standard input
and output (L18); [agents.md](../agents.md) is the full guide with the
safety advice. From the repository root, after `cargo build --release`:

```
mkdir -p output/agent/runs
target/release/foehn mcp --read-root examples/cases --read-root output/agent \
    --write-root output/agent/runs --policy output/agent/foehn-mcp.toml
```

Both roots are required (`MCP001`); the action table must lie outside the
write roots. Register the same command with the assistant host (absolute
paths, working folder the repository root). Then let it work, or try the
protocol by hand: the server reads one JSON-RPC message per line.

```
printf '%s\n' \
 '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"me","version":"0"}}}' \
 '{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"case.check","arguments":{"case":"examples/cases/flow/cavity_re100.toml","out":"output/agent/runs/cavity"}}}' \
 | target/release/foehn --quiet mcp --read-root examples/cases --write-root output/agent/runs \
     --policy output/agent/foehn-mcp.toml
```

Runs need `[run.limits] max_wall_time` and `max_output` (`MCP002`, with the
fix as `set`), and each kind of action follows the table:

```toml
# output/agent/foehn-mcp.toml
[actions]
start-run = "ask"            # the host asks you; "allow" for a host without elicitation
overwrite-results = "ask"
write-outside-root = "refuse"

[estimate]
max_cells = 2000000          # above it, a run is the kind exceed-estimate
```

A host that cannot ask gets `MCP005` for a call of kind `ask`; set the
kind to `"allow"` (the table is read at every call). Histories come page
by page:

```
target/release/foehn call result.page --input examples/cases/tools/page_cavity.json   # after recipe 171's run
```

and the session replays from its journal (default
`foehn-mcp-journal.jsonl` beside the table):

```
target/release/foehn run --script output/agent/foehn-mcp-journal.jsonl \
    --read-root examples/cases --read-root output/agent --write-root output/agent/runs
```

The scripted session of gallery 175 on the heat sink:

```
python3 -m tools.plot.mcp_session docs/figures/gallery/mcp_session.svg target/release/foehn
```

[↑ back to top](#run-control-and-health)

---

<a id="control-port"></a>
## 177. Watch and steer a running case

A run opens a control port when you ask for it (L19): `foehn flow
--control-port`, or any `flow.run` through `foehn mcp`. From the repository
root, after `cargo build --release`, start the transient cavity with a
small time step so that it lasts long enough:

```
mkdir -p output/tools/control
target/release/foehn flow examples/cases/flow/cavity_re100_transient.toml \
    --out output/tools/control/live --control-port \
    --set time.time_step=0.005 --set time.steps=1000000 \
    --set run.limits.max_wall_time=120 --set run.limits.max_output=2e9
```

In a second terminal ask the run how it is doing, change the time step,
look again, and stop it. The inputs are in
`examples/cases/tools/control/`:

```
target/release/foehn call run.status --input examples/cases/tools/control/run_status.json
target/release/foehn call run.set    --input examples/cases/tools/control/run_set.json
target/release/foehn call run.set    --input examples/cases/tools/control/run_set_refused.json
target/release/foehn call run.status --input examples/cases/tools/control/run_status.json
target/release/foehn call run.stop   --input examples/cases/tools/control/run_stop.json
```

The status lists under `control` what the run did with each change: the
time step is `applied` after the current step, with the `[[schedule]]` entry
that repeats it from a case file; the relaxation of 0.7 is `rejected`,
because a transient case with one outer iteration would lag the solution
(`RUN024`). The stop ends the first terminal's run with exit code 3, a
restart file and the result summary. Only run-time settings can be changed
(the keys a re-read of the case file may change, D40), never `run.limits.*`.
The port exists only while the run goes: `control.json` and `control.token`
are in the output folder, the token readable by its user only, and both
disappear when the run ends. An assistant does the same through
[`foehn mcp`](../agents.md#7-watch-and-steer-a-running-case); the figure of
[gallery 177](../gallery/runs.md#control-port) is made by
`python3 -m tools.plot.control_session docs/figures/gallery/control_port.svg`.

[↑ back to top](#run-control-and-health)

---

<a id="heat-sink-study"></a>
## 181. Let an assistant run a parameter study

A study through `foehn mcp` needs a folder for its runs and an action table
that says what the assistant may start. From the repository root, after
`cargo build --release`:

```
mkdir -p output/agent/study
cat > output/agent/foehn-mcp.toml <<'TABLE'
[actions]
start-run = "allow"          # or "ask" with a host that can ask you
exceed-estimate = "refuse"
write-outside-root = "refuse"
overwrite-results = "refuse"

[estimate]
max_cells = 40000            # a run above this is refused, not asked
TABLE
target/release/foehn mcp --read-root examples/cases --read-root output/agent/study \
    --write-root output/agent/study --policy output/agent/foehn-mcp.toml \
    --journal output/agent/journal.jsonl
```

Register that command with the assistant host ([agents.md](../agents.md)),
and give it the task: "vary the fin pitch of
`examples/cases/flow/cht_heat_sink_turbulent.toml` by overriding
`mesh.generate.zones`, check each case first, run each into
`output/agent/study/pitch_<mm>` with run limits, and report the mean
temperature of the base top and the pressure drop with the limitations of
every result." The run limits and the study budget are `set` overrides the
assistant adds itself (`run.limits.max_wall_time`, `max_output`,
`study_cell_hours`). While a run goes it can read `run.status` and
`result.page`, and `run.stop` it.

To see what a session looks like, replay the recorded one
([gallery 181](../gallery/runs.md#heat-sink-study)); it takes about 8 minutes:

```
target/release/foehn run --script examples/cases/tools/mcp/heat_sink_study/journal.jsonl \
    --read-root examples/cases --read-root output/agent/study --write-root output/agent/study
```

The replay skips the refused and the live calls and reports any call that
ends otherwise than recorded. Compare `output/agent/study/pitch_40/residuals.csv`
with the one of your own run: the files are byte-identical.
Without a host of your own, `tools/agent/stdio_host.py` keeps `foehn mcp`
open behind a socket and sends one call per command (its docstring has the
usage); it is what the recorded session used.

[↑ back to top](#run-control-and-health)
