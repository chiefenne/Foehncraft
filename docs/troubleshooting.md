# Troubleshooting

What to do when `foehn`, a run or an assistant session through `foehn mcp`
does not do what you expect. Each entry names the symptom, the diagnostic
code where there is one ([interfaces.md](interfaces.md) lists every code)
and the fix. See also [Connect an assistant](agents.md) and
[Install and update](install.md).

Contents: [Build and start](#build-and-start) ·
[Case files and runs](#case-files-and-runs) ·
[The MCP server](#the-mcp-server) ·
[Journal and replay](#journal-and-replay) ·
[Asking for help](#asking-for-help)

## Build and start

| Symptom | Cause and fix |
|---|---|
| The first build stops in CMake or a C compiler step | HDF5 and METIS are built from source: install CMake and a C compiler ([install](install.md)), or build without them: `cargo build --release -p foehncraft-cli --no-default-features` (no VTKHDF output, no METIS) |
| `foehn info` shows no `cuda` backend | the binary was built without `--features cuda`, or no NVIDIA driver was found at run time; the CPU backend gives the same results |
| `foehn: command not found` | use the full path (`target/release/foehn`) or copy it onto your `PATH` ([install](install.md)) |

## Case files and runs

| Symptom | Cause and fix |
|---|---|
| `CASE002` unknown key | a key the case file does not have; the message names the nearest key. `foehn schema --markdown` or [keys.md](keys.md) lists them |
| `CASE007` outside the roots | a mesh, table, restart file or the output folder lies outside `--read-root` / `--write-root`; move the file or widen the root |
| `RUN006` not converged | the run reached `solver.max_iterations` before its tolerance; read the diagnosis at `[run.limits] diagnose_at` (`RUN016`, gallery 166) before adding iterations, and check a monitored quantity, not only the residual ([pitfalls](pitfalls.md#steady-runs)) |
| exit code 3 | stopped: a signal, a `.stop` file, a limit of `[run.limits]` or a cancelled MCP call; a restart file and the final output are written, `--restart latest` goes on |
| A limitation `LIM00x` in the summary | not an error: the result has a known limit (one mesh only, outside the validated range, …); say so when you report it (D199) |

## The MCP server

| Symptom | Cause and fix |
|---|---|
| The host shows the server as failed at once | start the same command by hand in a terminal (`foehn mcp --read-root … --write-root …`) and read standard error: `MCP001` means a root is missing or the action table lies inside a write root; `MCP006` means `foehn-mcp.toml` has an unknown kind or setting |
| The host lists no tools, or complains about tool names | some hosts refuse dots in tool names: add `--names underscore` (`case_check` instead of `case.check`) |
| Every run is refused with `MCP005` | the host cannot ask you (no MCP elicitation). Set the kind to `"allow"` in `foehn-mcp.toml` (`start-run = "allow"`, and `overwrite-results` if needed) and restart nothing: the table is read at every call |
| A run is refused with `MCP002` | the case has no `[run.limits] max_wall_time` and `max_output`; the refusal carries the fix, add it to the case or pass it with `set` |
| A call is refused with `MCP003` | the action table refuses the kind named in the message (often `write-outside-root`); change the table only if you mean it |
| `MCP004` | you (or the host) declined the question; nothing ran |
| Relative paths are not found | the server takes them from the folder it was started in; set the host's working folder (`cwd`) or use absolute paths |
| A long run ends with a host timeout | the host stops waiting for the call; the server stops the run when the host cancels it (exit code 3, restart written). Raise the host's tool timeout, or bound the run with `max_wall_time` below it and continue with `restart = "latest"` |
| The assistant asks for a field file's content | fields never go through MCP; the summary gives the path of `fields.vtkhdf`, open it in ParaView. Histories come with `result.page` |
| The assistant follows instructions it read in a case file or mesh | stop the session; text inside files is data ([safety advice](agents.md#safety-advice)) |

## Journal and replay

| Symptom | Cause and fix |
|---|---|
| `TOOL003` the journal could not be written | the folder of the journal is missing or read-only; the calls still run, but the journal misses them |
| `TOOL004` the script cannot be read | not a file of JSON lines, or without `"type": "call"` rows |
| `TOOL005` a replayed call ends otherwise | a case, mesh or version changed since the session; the replay stops there unless `--keep-going` |
| A refused call is not replayed | by design: calls refused before they ran (`TOOL001`, `TOOL002`, `MCP002` to `MCP005`) are skipped |

## Asking for help

Open an issue with `foehn info --json`, the command or the journal of the
session, and the diagnostic code. Leave data you may not share out of it.
