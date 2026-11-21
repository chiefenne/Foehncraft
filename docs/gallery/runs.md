# Run control, health and adjoint readiness

Part of the [Example gallery](../gallery.md): adjoint-ready kernels, runs that a program drives, balances and the divergence report.

| # | Example | Stage | Feature |
|---|---|---|---|
| 17 | [Adjoint-ready kernels: derivatives and checkpointing](#adjoint-ready) | 5 | dual numbers through the kernel bodies, transposed matrices, binomial checkpointing |
| 18 | [A run that a program can drive](#run-contract) | 5 | `--json`, exit codes, `--set`, limits, signals, `--restart latest`, diagnostic codes, estimate before the run, byte-identical files |
| 19 | [Balances and the divergence report](#health) | 5 | mass and kinetic-energy ledger, numerical dissipation of the schemes (second law), divergence report with ranked causes |
| 27 | [A log that can be circulated](#log) | 5 | frames, blocks of settings, tables with headers, summaries, closing frame with exit code (D125) |
| 28 | [A table with the columns you choose](#log-columns) | 5 | `[log] columns`: probes, cells, extremes with their cells, boundary mass flow; the same content in `run.jsonl` (D125) |
| 29 | [Settings that change during the run](#schedule) | 5 | `[[schedule]]`: under-relaxation, tolerances and cadences change at prepared iterations, steps or times; every change in the log and in `run.jsonl` (D128) |
| 166 | [Iteration budget, diagnosis and study budgets](#budgets) | 7 | the iteration budget before the run, a steady run diagnosed at 2000 iterations (trend, verdict, advice, stop), cell-hour, GPU-hour and output budgets per study beside the D81 estimate (D198) |
| 171 | [Every function as a tool: the catalogue an agent reads](#tool-registry) | 7 | one registry of typed tools with ids, versions, JSON Schemas and action kinds; `foehn commands`, `foehn call`, `fc.call`; front ends checked by `cargo xtask layering` (D196, D203) |

---

<a id="adjoint-ready"></a>
## 17. Adjoint-ready kernels: derivatives and checkpointing

<img src="../figures/gallery/checkpoint_steps.svg" alt="Primal steps per time step of a reverse sweep against the snapshots in memory" width="680" height="420">

Preparation for a discrete adjoint solver (D59, L17). The kernel bodies are
generic over the floating-point type, so they also run on dual numbers and
return the exact derivatives of their results; matrices can be applied
transposed as a gather; and a transient run can be reversed with few
states in memory.

The chart shows the cost of reversing a transient run with binomial
checkpointing (Griewank & Walther 2000): with 20 snapshots, the reverse
sweep of 10 000 time steps runs each step 3.8 times on average and never
more than 4 times; with 40 snapshots, 100 000 steps cost 3.9 primal steps
each (values above 12 drawn at 12). The schedule is optimal: its count of
recomputed steps equals that of a dynamic program over all positions of
the snapshots.

| Check | Result |
|---|---|
| Jacobians of 11 face, cell and boundary kernels from dual numbers, every input array, against central differences | 320 580 derivatives, largest relative difference 5.0e-7 |
| Shape derivatives of areas, centroids and volumes from the node coordinates | 165 024 derivatives, 2.8e-7; total volume unchanged by interior nodes |
| ⟨Av, w⟩ = ⟨v, Aᵀw⟩ on the momentum and pressure matrices of a SIMPLEC solve | to 10⁻¹⁴; Aᵀ without storing it equals the stored Aᵀ bit for bit, on 1 and 4 threads |
| CG + multigrid and BiCGStab on Aᵀ | converged to 10⁻¹⁰ |
| Reverse sweep of 13 BDF2 steps with 2 snapshots | every level bitwise equal to the forward run, 24 primal steps (the optimum) |

| | |
|---|---|
| Code | `foehncraft_core::Dual`, `foehncraft_kernels::adjoint`, `foehncraft_flow::checkpoint`, `foehncraft_mesh::geometry` (generic over `Coord`) |
| Tests | `crates/foehncraft-flow/tests/adjoint.rs`, `crates/foehncraft-flow/src/adjoint_tests.rs`, `checkpoint.rs`, `geometry.rs` |
| Report | [stage 5, part 3](../verification/stage5.md#adjoint-ready-kernels-part-3) |

```
cargo test --release -p foehncraft-flow --test adjoint -- --nocapture   # Jacobians, one line per kernel
cargo test --release -p foehncraft-flow --lib -- adjoint_tests checkpoint   # transposed matrices, checkpointing
python3 -m tools.plot.checkpoint_chart checkpoint_steps.svg
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="run-contract"></a>
## 18. A run that a program can drive

<img src="../figures/gallery/estimate.svg" alt="Memory and file sizes estimated by foehn flow --check against eleven runs" width="640" height="460">

A parametric study, a batch job or an agent drives `foehn` without
reading its log: it sets keys with `--set`, reads one JSON object, acts on
the exit code and continues a stopped run with `--restart latest`. Before
it starts a run, `--check` tells it what the run will cost. The chart
shows that estimate against eleven runs of the catalogue, from 128 to
28 227 cells, hexahedra to polyhedra: peak memory within 6 %, field and
restart files within 9 %. The memory model is a base of 14 MB plus
1960 bytes per interior face and 710 per boundary face.

```
$ foehn --json flow cavity_re100_budget.toml --set run.limits.max_iterations=300; echo "exit $?"
{
  "format": "foehncraft.result",
  "format_version": "1.0",
  "status": "stopped",
  "exit_code": 3,
  "stop": { "code": "RUN004", "kind": "iterations", "message": "stopped at the limit of 300 outer iterations ([run.limits] max_iterations)" },
  ...
exit 3
$ foehn --json flow cavity_re100_budget.toml --restart latest; echo "exit $?"
  ...
  "status": "finished",
exit 0
```

| Check | Result |
|---|---|
| Stop at a limit or on SIGTERM, then `--restart latest` | fields and restart file bitwise equal to the run without the stop (steady and transient) |
| Two identical runs | `result.json`, `fields.vtkhdf`, `restart.h5` and every CSV byte for byte the same; only the log, the reports and `run.jsonl` (timings) differ |
| Exit codes | 0 finished, 1 error (as a JSON object with code and key), 2 failed check (`RUN006`, `RUN007`), 3 stopped (`RUN001`–`RUN005`) |
| A typo in the case file | `CASE002` with the key `mesh.generate.cels`, line 4, column 1 and "did you mean `cells`" |
| The 20 flow cases of the catalogue, this build against main | every field, CSV and JSON value bitwise equal; the resolved case gains the empty `[run.limits]` table |

| | |
|---|---|
| Code | `foehncraft_flow::engine` (run, check, estimate, limits, result, run record), `diag` (codes), `schema` (JSON Schema), `foehncraft-cli` |
| Tests | `crates/foehncraft-cli/tests/run_contract.rs`, unit tests in `engine.rs`, `diag.rs`, `schema.rs` |
| Formats | [interfaces.md](../interfaces.md), [keys.md](../keys.md) |
| Case | `examples/cases/flow/cavity_re100_budget.toml` |
| Recipe | [26](../recipes/runs.md#run-contract) |

```
cargo test --release -p foehncraft-cli --test run_contract
python3 -m tools.plot.estimate_chart estimate.svg
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="health"></a>
## 19. Balances and the divergence report

<img src="../figures/gallery/energy_budget.svg" alt="Numerical dissipation of kinetic energy over the mesh size, lid-driven cavity and Taylor–Green vortex" width="640" height="440">

Every run books mass and kinetic energy: what flows and works through
each boundary, what the body force adds, what viscosity dissipates, and what
is left, the dissipation of the numerical schemes. The second law asks
that this rest is not negative. In the steady lid-driven cavity the lid
works on the fluid and everything ends up as dissipation; the schemes take
about 1.5 % of it, falling slowly with the mesh because the two corners
where the lid meets the walls are singular. The ledger splits the rest in
two: the pressure part, which the momentum interpolation of the face
fluxes dissipates, and the part of the convection and time schemes. In
the cavity the pressure part is nearly all of it and stays of the same
size at the singular corners, while the linear-upwind correction gives
back a little, less with every refinement. On the smooth Taylor–Green
vortex both parts are positive: the numerical dissipation falls with
order 2.5, its pressure part with order 3.

```
$ foehn flow cavity_re100_energy.toml
      balance of the final state
        kinetic energy
          ymax viscous work in                         +2.919438e-3 W
          viscous dissipation                          +2.877098e-3 W
          storage                                       +0.000000e0 W
          numerical dissipation (+1.45e-2)             +4.234029e-5 W
            of it pressure                             +5.115046e-5 W
            of it schemes                              -8.810171e-6 W
```

When a run diverges it keeps its last valid state and says why:

```
$ foehn flow cavity_re1000.toml --set 'solver.algorithm="simple"' --set solver.momentum.relaxation=0.99 --set solver.pressure.relaxation=1.0
      DIVERGED in outer iteration 16: 1 failing cells, velocity first; report in divergence.md
        DIV002 under-relaxation (0.80): relaxation factors velocity 0.99, pressure 1 (SIMPLE); try --set solver.momentum.relaxation=0.69 --set solver.pressure.relaxation=0.31
        DIV006 linear solve (0.60): momentum BiCGStab at its limit of 50 iterations; try --set solver.momentum.max_iter=100
Error RUN014: the solution diverged in outer iteration 16 ...
```

| Check | Result |
|---|---|
| Cavity Re 100, 16² to 128² | numerical dissipation 2.05 %, 1.59 %, 1.45 %, 1.33 % of the lid power, positive; pressure part 7.6e-5 to 5.0e-5 W, scheme part −3.6e-5 to −4.4e-6 W; mass ledger zero |
| Taylor–Green, 16² to 64², fixed Courant number | numerical dissipation 4.1e-4, 7.4e-5, 1.2e-5 J over 0.4 s (pressure part 1.3e-4, 1.6e-5, 1.4e-6 J); mass imbalance at round-off |
| Steady divergence (SIMPLE, no relaxation) | `RUN014`, the first cause is the relaxation; its fix runs |
| Transient divergence (time step 100 times too large) | the first cause is the time step, with a smaller one as the fix |
| The catalogue against part 4a | fields, CSV files and report values bitwise equal; one new warning, `BAL003` on the Courant-10 Taylor–Green sweep |

| | |
|---|---|
| Code | `foehncraft_flow::balance` (ledger, tables, CSV, chart), `divergence` (analysis, snapshot, report), `simplec::Flow::keep_valid`, `rewind` |
| Tests | `crates/foehncraft-flow/tests/health.rs`, unit tests in `balance.rs`, `divergence.rs` |
| Formats | [interfaces.md](../interfaces.md#conservation-ledger-d88-d89-d90-part-4) |
| Case | `examples/cases/flow/cavity_re100_energy.toml` |
| Recipe | [27](../recipes/runs.md#health) |

```
cargo test --release -p foehncraft-flow --test health
python3 -m tools.plot.energy_budget_chart energy_budget.svg
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="log"></a>
## 27. A log that can be circulated

Part of what a run says about itself (D125): the log shows what the case file set, so a log mailed on its own is enough to see the run. The log opens with the name of the code in block letters ("Foehn" over "craft", the Alpha figlet font by Lennert Stock); every block is a box of asterisks as wide as its content, centred on the page (no colour, no fixed width: the banner is 133 columns, the other blocks fit their text); the console and `flow.log` are identical. Excerpt of the transient cavity (Re 100, 64 x 64, BDF2, three outer iterations), blocks between the banner and the table, most time steps and the lower summaries left out:

```
*************************************************************************************************************************************
*                                                                                                                                   *
*             _____                   _______                   _____                    _____                    _____             *
*            /\    \                 /::\    \                 /\    \                  /\    \                  /\    \            *
*           /::\    \               /::::\    \               /::\    \                /::\____\                /::\____\           *
*          /::::\    \             /::::::\    \             /::::\    \              /:::/    /               /::::|   |           *
*         /::::::\    \           /::::::::\    \           /::::::\    \            /:::/    /               /:::::|   |           *
*        /:::/\:::\    \         /:::/~~\:::\    \         /:::/\:::\    \          /:::/    /               /::::::|   |           *
*       /:::/__\:::\    \       /:::/    \:::\    \       /:::/__\:::\    \        /:::/____/               /:::/|::|   |           *
*      /::::\   \:::\    \     /:::/    / \:::\    \     /::::\   \:::\    \      /::::\    \              /:::/ |::|   |           *
*     /::::::\   \:::\    \   /:::/____/   \:::\____\   /::::::\   \:::\    \    /::::::\    \   _____    /:::/  |::|   | _____     *
*    /:::/\:::\   \:::\    \ |:::|    |     |:::|    | /:::/\:::\   \:::\    \  /:::/\:::\    \ /\    \  /:::/   |::|   |/\    \    *
*   /:::/  \:::\   \:::\____\|:::|____|     |:::|    |/:::/__\:::\   \:::\____\/:::/  \:::\    /::\____\/:: /    |::|   /::\____\   *
*   \::/    \:::\   \::/    / \:::\    \   /:::/    / \:::\   \:::\   \::/    /\::/    \:::\  /:::/    /\::/    /|::|  /:::/    /   *
*    \/____/ \:::\   \/____/   \:::\    \ /:::/    /   \:::\   \:::\   \/____/  \/____/ \:::\/:::/    /  \/____/ |::| /:::/    /    *
*             \:::\    \        \:::\    /:::/    /     \:::\   \:::\    \               \::::::/    /           |::|/:::/    /     *
*              \:::\____\        \:::\__/:::/    /       \:::\   \:::\____\               \::::/    /            |::::::/    /      *
*               \::/    /         \::::::::/    /         \:::\   \::/    /               /:::/    /             |:::::/    /       *
*                \/____/           \::::::/    /           \:::\   \/____/               /:::/    /              |::::/    /        *
*                                   \::::/    /             \:::\    \                  /:::/    /               /:::/    /         *
*                                    \::/____/               \:::\____\                /:::/    /               /:::/    /          *
*                                     ~~                      \::/    /                \::/    /                \::/    /           *
*                                                              \/____/                  \/____/                  \/____/            *
*                                                                                                                                   *
*               _____                    _____                    _____                    _____                _____               *
*              /\    \                  /\    \                  /\    \                  /\    \              /\    \              *
*             /::\    \                /::\    \                /::\    \                /::\    \            /::\    \             *
*            /::::\    \              /::::\    \              /::::\    \              /::::\    \           \:::\    \            *
*           /::::::\    \            /::::::\    \            /::::::\    \            /::::::\    \           \:::\    \           *
*          /:::/\:::\    \          /:::/\:::\    \          /:::/\:::\    \          /:::/\:::\    \           \:::\    \          *
*         /:::/  \:::\    \        /:::/__\:::\    \        /:::/__\:::\    \        /:::/__\:::\    \           \:::\    \         *
*        /:::/    \:::\    \      /::::\   \:::\    \      /::::\   \:::\    \      /::::\   \:::\    \          /::::\    \        *
*       /:::/    / \:::\    \    /::::::\   \:::\    \    /::::::\   \:::\    \    /::::::\   \:::\    \        /::::::\    \       *
*      /:::/    /   \:::\    \  /:::/\:::\   \:::\____\  /:::/\:::\   \:::\    \  /:::/\:::\   \:::\    \      /:::/\:::\    \      *
*     /:::/____/     \:::\____\/:::/  \:::\   \:::|    |/:::/  \:::\   \:::\____\/:::/  \:::\   \:::\____\    /:::/  \:::\____\     *
*     \:::\    \      \::/    /\::/   |::::\  /:::|____|\::/    \:::\  /:::/    /\::/    \:::\   \::/    /   /:::/    \::/    /     *
*      \:::\    \      \/____/  \/____|:::::\/:::/    /  \/____/ \:::\/:::/    /  \/____/ \:::\   \/____/   /:::/    / \/____/      *
*       \:::\    \                    |:::::::::/    /            \::::::/    /            \:::\    \      /:::/    /               *
*        \:::\    \                   |::|\::::/    /              \::::/    /              \:::\____\    /:::/    /                *
*         \:::\    \                  |::| \::/____/               /:::/    /                \::/    /    \::/    /                 *
*          \:::\    \                 |::|  ~|                    /:::/    /                  \/____/      \/____/                  *
*           \:::\    \                |::|   |                   /:::/    /                                                         *
*            \:::\____\               \::|   |                  /:::/    /                                                          *
*             \::/    /                \:|   |                  \::/    /                                                           *
*              \/____/                  \|___|                   \/____/                                                            *
*                                                                                                                                   *
*************************************************************************************************************************************
*   Version                            0.5.0  (build b2191c6-dirty)                                                                 *
*   Date and time                      2026-11-16 22:21:29 UTC                                                                      *
*   Host name                          vm                                                                                           *
*   Platform                           linux-x86_64, Intel(R) Xeon(R) Processor @ 2.30GHz                                           *
*   Features                           hdf5, metis                                                                                  *
*   Threads / backend                  4 / cpu                                                                                      *
*************************************************************************************************************************************

  ...

                                   **********************************************************************
                                   *         1. TIME STEP     TIME = 1.00000E-01     DT = 1.00000E-01   *
                                   *            COURANT NUMBER = 4.762     OUTER ITERATIONS = 3         *
                                   **********************************************************************

         -- scaled residuals ----------------- -- linear iterations
   Outer      ResU     ResV     ResW     ResM      U    V    W    P 
       1  1.00E+00 0.00E+00 0.00E+00 6.33E-02     37    0    0   22
       2  1.37E-02 1.03E-02 0.00E+00 1.56E-02     43   42    0   23
       3  2.94E-03 3.07E-03 0.00E+00 5.57E-03     41   43    0   22

         -- time ------------------- -- scaled residuals ----------------- -- solvers --
    Step    Time [s]   Dt [s] Cour.       ResU     ResV     ResW     ResM   Outer    CG 
       1  1.0000E-01 1.00E-01 4.762   1.00E+00 0.00E+00 0.00E+00 6.33E-02       3    67
      balance  mass +0.00e0  kinetic energy: dK/dt +9.912e-4 W, boundary +1.308e-2 W, dissipation 1.206e-2 W, numerical +3.204e-5 W (+2.4e-3)

  ...

                       **********************************************************************************************
                       *                                                                                            *
                       *   Time stepping:                                                                           *
                       *                                                                                            *
                       *   Status                             DONE: 500 steps to t = 5.00000E+01 s                  *
                       *   Time step                          1.000E-01 to 1.000E-01 s                              *
                       *   Largest Courant number             6.095                                                 *
                       *   Solve time                         21.74 s  (setup 0.02 s, 10618 ns per cell and step)   *
                       *                                                                                            *
                       **********************************************************************************************

          ************************************************************************************************************************
          *                                                                                                                      *
          *   Run finished: every check passed                                                                                   *
          *                                                                                                                      *
          *   Exit code                          0                                                                               *
          *   Elapsed time                       0:00:22                                                                         *
          *   Results in                         /tmp/claude-0/-home-user/8486083e-90f2-5e7f-b60b-e659997a40c3/scratchpad/run5   *
          *   Report                             flow.md, result.json                                                            *
          *                                                                                                                      *
          ************************************************************************************************************************
```

Before the table, the log lists the input and output files, the run control, the mesh information with its quality, the fluid properties, the solution settings as a matrix per equation and the boundary conditions with their values. After it come the summaries: convergence or time stepping, verification, results per boundary, the ledger, the run time and the warnings. Each time step that is logged starts with a box (step, time, time step, Courant number, outer iterations), followed by the residuals and the linear iterations of every outer iteration of that step (the heavier the transient run, the more there is to watch: *Residuen-gaffen*, see [sayings](../sayings.md)), the summaries are boxes with a heading, blank rows and dashed rules. Messages read ` %FOEHN-I-CASE, text` (`I`, `W` or `E`, then the capitalised name of the source file that writes the message), so `%FOEHN` finds them in the log and in the code; a failed run ends with `Run FAILED: <code>`, the reason and the exit code. `[log] style = "box"` draws the boxes with box characters.

| | |
|---|---|
| Code | `crates/foehncraft-flow/src/console.rs`, `logblocks.rs`, `case.rs`, `unsteady.rs` |
| Tests | `crates/foehncraft-flow/src/console.rs` (frames, tables), `tests/health.rs`, `tests/stage5.rs` |
| Case files | every file of `examples/cases/flow/` |

```
target/release/foehn flow examples/cases/flow/cavity_re100_transient.toml
less output/cavity_re100_transient/flow.log
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="log-columns"></a>
## 28. A table with the columns you choose

The table of iterations or time steps shows what the case asks for (D125). The cavity at Re = 100 with a probe in the lower-left vortex, the pressure in one cell, the largest velocity and the lowest pressure of the domain with the cell where each occurs, and the mass flow through the lid (`examples/cases/flow/cavity_re100_columns.toml`):

```
         -- scaled residuals --------------------- -- probe vortex -- cell 2000 -- extremes ----------------------------- -- mass flow -- run [s]
    Iter       ResU      ResV      ResW      ResM           Ux            P         Umax     cell       Pmin     cell         ymax      s/iter
       1  1.000E+00 0.000E+00 0.000E+00 1.453E-01   -6.842E-02   -2.911E-03    6.253E-01      248 -2.707E-01      240    0.000E+00           -
      10  3.704E-03 4.362E-03 0.000E+00 7.042E-03   -1.873E-01   -9.042E-03    7.868E-01      250 -2.703E-01      240    0.000E+00   1.796E-02
     ...
```

(The numbers are those of the 16 × 16 test cavity of the test suite, to show the layout; the case file runs on 64 × 64 cells.) A group heads the columns of one kind and target; an extreme gives two columns, its value and its cell. The table is as wide as the columns need; `[log] digits` sets the significant digits. The columns are read-only: the fields are not changed, so the run is bitwise the one without `[log] columns`.

The same content goes to `run.jsonl` (format 1.1): a `columns` row describes the table once, a `table_row` row carries the values of every line (an extreme carries its cell and the cell centre), and a `section` row carries every framed block of the log with label, value and unit. A monitor reads these rows and draws the blocks and the table without parsing text.

| | |
|---|---|
| Code | `crates/foehncraft-flow/src/logcols.rs`, `console.rs`, `case.rs`, `unsteady.rs` |
| Tests | `crates/foehncraft-flow/src/logcols.rs` (names), `crates/foehncraft-cli/tests/run_contract.rs` (`configurable_log_columns_change_no_result`) |
| Case files | `flow/cavity_re100_columns.toml` |

```
target/release/foehn flow examples/cases/flow/cavity_re100_columns.toml
tail -f output/cavity_re100_columns/run.jsonl      # the rows of a monitor
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="schedule"></a>
## 29. Settings that change during the run

A whole run can be set up before it starts (D128). `[[schedule]]` lists the settings that change and the iteration, step or time from which the new value holds (`examples/cases/flow/cavity_re100_schedule.toml`):

```toml
[solver.momentum]
relaxation = 0.5               # cautious start

[[schedule]]
at_iteration = 100
set = { solver.momentum.relaxation = 0.9 }

[[schedule]]
at_iteration = 200
set = { solver.log_every = 10 }
```

The "Run control" frame lists the changes before the run, and every change is a line in the log when it takes effect:

```
 Scheduled: iteration 100          solver.momentum.relaxation = 0.9
 Scheduled: iteration 200          solver.log_every = 10
 ...
 [I] iteration 100: schedule (iteration 100): solver.momentum.relaxation: 0.5 -> 0.9
```

The same change is an `event` row of `run.jsonl` with the code `RUN015` and the fields `key`, `old` and `new`, so a monitor marks the iteration where the setting changed. The case file alone says what was run (the entries are part of the resolved case and of its hash), and a restart continues the schedule: entries before the restart point apply at once.

Which keys may be scheduled: those that a re-read of the case file may change during a run (D40) - relaxation factors, linear-solver tolerances and iteration limits, the stop and convergence criteria, the time step and its adaptive controls, the output, restart and balance cadences, `solver.log_every` - without the guard rails `run.limits.*` and `run.reread_every`. Anything that belongs to the state the run starts from (mesh, boundary values, schemes, fluid) is refused when the case is read, with the reason. A transient case uses `at_step` or `at_time` (the run lands a step on the time); steady `at_iteration`.

| | |
|---|---|
| Code | `crates/foehncraft-flow/src/schedule.rs`, `case.rs` (`schedule_step`), `unsteady.rs`, `control.rs` |
| Tests | `crates/foehncraft-flow/tests/schedule.rs`, `crates/foehncraft-cli/tests/run_contract.rs` (`schedule_changes_are_events_of_the_run_record`) |
| Case files | `flow/cavity_re100_schedule.toml` |

```
target/release/foehn flow examples/cases/flow/cavity_re100_schedule.toml
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="budgets"></a>
## 166. Iteration budget, diagnosis and study budgets

<img src="../figures/gallery/diagnosis_cylinder.svg" alt="Residuals of a steady run of the cylinder at Re 400 that stall; the diagnosis at iteration 2000 fits their trends over the last 500 iterations, finds them stalled and stops the run before the rest of its 3000-iteration budget" width="840" height="426">

A steady run must not go on for 10 000 iterations when it is not going
anywhere (D198 (d)). The iteration budget is a visible setting: `foehn flow
--check` and the log show it with the key it comes from
(`solver.max_iterations` or `[run.limits] max_iterations`), and at
`[run.limits] diagnose_at` (default 2000) a steady run that has not
converged is diagnosed from its residual history. For each equation the
trend over the last 500 iterations is fitted (decades per 1000 iterations
and the scatter about it), classified as falling, stalled, rising or
oscillating, and for a falling residual the iteration is projected at which
it reaches the tolerance. The verdict (on track, slow, stalled, rising,
oscillating) and advice go into the log and into the run record
(`run.jsonl`, row type `diagnosis`, event `RUN016`); with `diagnosis =
"stop"` a run off track stops there with its restart file and results (exit
code 3, `RUN017`).

The chart is `examples/cases/flow/cylinder_steady_diagnosis.toml`: the
cylinder in a channel of Schäfer & Turek at Re = 400 (24 576 cells), run
steady on purpose. The flow sheds vortices, so the residuals stall at
1e-2: the diagnosis says so at iteration 2000 (about 0.00 decades per 1000
iterations for momentum and continuity), advises to run it transient, and
stops it a third of its budget early. The same case at Re = 100 converges
in 993 iterations and is never diagnosed. The diagnosis only reads the
residuals: the fields and residuals of a run with and without it are
bitwise equal.

```
 Diagnosis at iteration 2000
   Verdict          STALLED
   Budget           3000 iterations
   Controlling      momentum_x
   Window           last 500 iterations
   momentum_x       1.27e-2, stalled, +0.00 decades per 1000 iterations
   momentum_y       1.06e-2, stalled, +0.00 decades per 1000 iterations
   continuity       1.72e-3, stalled, -0.01 decades per 1000 iterations
 %FOEHN-W-CASE, diagnosis at iteration 2000: momentum_x stalls at 1.3e-2: ... A stall this far
   above the tolerance often means the flow has no steady solution (vortex shedding): run it transient
```

**Study budgets** (D198 (c)). `[run.limits] study_cell_hours`,
`study_gpu_hours` and `study_output` bound what all runs of one study spend
together. A study is the set of output directories under one parent, as the
variants of `foehncraft.study.Study` are; each finished run adds a line with
its spend to `study_budget.jsonl` there. The run record shows, in its
`budget` row right after the header, the D81 estimate of the run (cells,
wall time, cell-hours, GPU-hours, output) beside each budget, what the
earlier runs spent and what is left. A run whose estimate does not fit warns
(`RUN020`), a run that spends the rest stops as at a limit (`RUN018`, exit
code 3), and the next run of a spent study does not start (`RUN019`, also
reported by `--check`). A cell-hour is one cell for one hour of wall time,
independent of the threads. Without these keys nothing changes.

| | |
|---|---|
| Code | `crates/foehncraft-flow/src/budget.rs`, `engine.rs` (`start_budgets`, `Steer::poll`, the estimate), `case.rs` (`[run.limits]`, the steady loop) |
| Tests | `crates/foehncraft-cli/tests/budgets.rs` (budget in `--check`, diagnosis recorded with bitwise equal fields, stop by the diagnosis, study ledger, stop and refusal), unit tests in `budget.rs` |
| Case files | `flow/cylinder_steady_diagnosis.toml` (CI set `long`, expected exit code 3) |
| Chart | `python3 -m tools.plot.diagnosis_chart docs/figures/gallery/diagnosis_cylinder.svg output/cylinder_steady_diagnosis` |

```
target/release/foehn flow examples/cases/flow/cylinder_steady_diagnosis.toml --check
target/release/foehn flow examples/cases/flow/cylinder_steady_diagnosis.toml
grep '"type":"diagnosis"' output/cylinder_steady_diagnosis/run.jsonl
```

[↑ back to top](#run-control-health-and-adjoint-readiness)

---

<a id="tool-registry"></a>
## 171. Every function as a tool: the catalogue an agent reads

<img src="../figures/gallery/tool_catalogue.svg" alt="The tools of the registry by category, coloured by action kind: read, writes files, starts a run" width="760" height="560">

Everything a front end can do is a tool of one registry (D196): a stable
id (`case.check`), a contract version, a typed input and output whose JSON
Schemas come from the same structs, and an action kind that says whether
it only reads, writes files or starts a run (the approvals of an agent
host follow it, D198). The chart is the catalogue of this version, drawn
from `foehn commands`: 19 tools in 10 categories, 12 that only read, 3
that write files and 4 that start a run. The command line, the Python
package and later the MCP server (L18) and the desktop app (L7) reach
exactly these tools; `cargo xtask layering` fails when a front end calls
the engine past them.

An agent learns the tools from the catalogue and calls any of them with
JSON. The result is one envelope with the status and exit code of D64,
the tool's output and, for a call that failed, its diagnostic:

```
$ foehn call case.check --input examples/cases/tools/check_cavity.json
{
  "format": "foehn.tool-result",
  "version": 1,
  "tool": "case.check",
  "tool_version": 1,
  "status": "finished",
  "exit_code": 0,
  "result": {
    "format": "foehncraft.check",
    ...
    "ok": true,
    "diagnostics": [],
    "estimate": { "cells": 4096, ... }
  },
  "diagnostics": [],
  "limitations": []
}
$ foehn call case.chek
Error TOOL001: unknown tool `case.chek`; did you mean `case.check`? (`foehn commands` lists the tools)
```

The same call from Python, in process:

```python
import foehncraft as fc

r = fc.call("flow.run", case="examples/cases/flow/cavity_re100.toml", out="output/tools/cavity")
r["status"], r["result"]["solves"][0]["iterations"]      # ('finished', 439)
fc.call("post.derive", case="examples/cases/flow/cavity_re100.toml",
        out="output/tools/cavity", fields=["vorticity", "q-criterion"])
```

`foehn flow`, `foehn mesh`, `foehn materials` and the other commands are
the same tools with flags: their text, `--json` objects, exit codes and
results are unchanged (a test compares each command with `foehn call` of
its tool, and the runs give byte-identical files).

| | |
|---|---|
| Code | `crates/foehncraft-tools/` (registry, tools), `crates/foehncraft-cli/src/main.rs` (`commands`, `call`), `crates/foehncraft-py/src/lib.rs`, `python/foehncraft/tools.py`, `xtask/` (layering) |
| Tests | `crates/foehncraft-tools/src/lib.rs` (ids, schemas, envelope, `docs/tools.md` current), `crates/foehncraft-cli/tests/tools_contract.rs`, `python/tests/test_tools.py` (every tool once through `call`) |
| Case files | `tools/*.json` (inputs of `foehn call`) |
| Catalogue | [tools.md](../tools.md) |

```
target/release/foehn commands | head -40
target/release/foehn call flow.run --input examples/cases/tools/run_cavity.json
python3 -m tools.plot.tools_chart docs/figures/gallery/tool_catalogue.svg target/release/foehn
```

[↑ back to top](#run-control-health-and-adjoint-readiness)
