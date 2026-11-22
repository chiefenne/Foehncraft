# Connect an assistant to Foehncraft

`foehn mcp` serves every tool of Foehncraft to an assistant (an AI agent
host) over the Model Context Protocol (MCP), on standard input and output
(L18). The assistant sees exactly the tools of `foehn commands`
([tools.md](tools.md)) with their input schemas: it can check a case, see
what a run will cost, start the run, read the result summary with its
limitations and page through the residual and monitor histories. Every
call goes into a journal, so you can replay the whole session yourself
with `foehn run --script`.

This page connects an assistant step by step. When something does not
work, see [troubleshooting](troubleshooting.md); to install or update
`foehn`, see [install and update](install.md).

Contents: [Before you start](#before-you-start) ·
[1. A folder of its own](#1-a-folder-of-its-own) ·
[2. Start the server by hand once](#2-start-the-server-by-hand-once) ·
[3. Register it with the assistant](#3-register-it-with-the-assistant) ·
[4. The action table](#4-the-action-table) ·
[5. A first session](#5-a-first-session) ·
[6. Replay the session](#6-replay-the-session) ·
[Safety advice](#safety-advice) ·
[What the server does and does not do](#what-the-server-does-and-does-not-do)

## Before you start

- A `foehn` binary built from this repository (`cargo build --release -p
  foehncraft-cli`, see [install](install.md)). `foehn mcp --help` must work.
- An assistant that can start a local MCP server (a "stdio" server). Hosts
  that support **MCP elicitation** can ask you before a run starts; with
  other hosts you allow runs in the action table (step 4).

## 1. A folder of its own

Give the assistant its own project folder, with copies of the cases and
meshes it should work on, and an empty folder for its runs:

```
mkdir -p ~/cfd/heatsink/runs
cp examples/cases/flow/cht_heat_sink_turbulent.toml ~/cfd/heatsink/
cd ~/cfd/heatsink
```

The assistant may **read** below the read roots and **write** only below
the write roots (D82, D198 (a)). Here: read `~/cfd/heatsink`, write
`~/cfd/heatsink/runs`. Never point a write root at your home folder, at a
folder with results you cannot run again, or at the repository.

## 2. Start the server by hand once

```
cd ~/cfd/heatsink
foehn mcp --read-root ~/cfd/heatsink --write-root ~/cfd/heatsink/runs
```

It prints one line on standard error (the roots, the action table, the
journal) and then waits for MCP messages on standard input; end it with
Ctrl-D. Without `--read-root` and `--write-root` it refuses to start
(`MCP001`). Relative paths in tool calls are taken from the folder the
server starts in, so start it in the project folder.

Options:

| Option | Meaning |
|---|---|
| `--read-root DIR` | may be read (repeatable; required) |
| `--write-root DIR` | may be written (repeatable; required) |
| `--policy FILE` | the action table (default `foehn-mcp.toml` in the working folder); it must lie outside every write root, so no tool can rewrite it |
| `--journal FILE` | the journal (default `foehn-mcp-journal.jsonl` beside the action table, or `FOEHN_JOURNAL`) |
| `--names underscore` | tool names `case_check` instead of `case.check`, for hosts that refuse dots in tool names |
| `--quiet` | no log lines on standard error |

## 3. Register it with the assistant

Every host has its own place for MCP servers; the entry is always the
same command with its arguments and the working folder. Use absolute
paths.

A command-line assistant that registers servers with a command, for
example:

```
claude mcp add foehn -- /path/to/foehn mcp \
    --read-root /home/me/cfd/heatsink --write-root /home/me/cfd/heatsink/runs \
    --policy /home/me/cfd/heatsink/foehn-mcp.toml
```

(started from `~/cfd/heatsink`, so that relative paths resolve there).

A desktop assistant with a JSON file of servers:

```json
{
  "mcpServers": {
    "foehn": {
      "command": "/path/to/foehn",
      "args": ["mcp",
               "--read-root", "/home/me/cfd/heatsink",
               "--write-root", "/home/me/cfd/heatsink/runs",
               "--policy", "/home/me/cfd/heatsink/foehn-mcp.toml"],
      "cwd": "/home/me/cfd/heatsink"
    }
  }
}
```

When the host has no working-folder setting, give the case paths to the
assistant as absolute paths. After a restart of the host, the assistant
lists about two dozen tools whose names start with `case.`, `flow.`,
`result.` and so on.

## 4. The action table

`foehn-mcp.toml` says, for each kind of action, whether a call runs,
asks you first or is refused (D198 (b)). Without the file every kind has
its default:

| Kind | When | Default |
|---|---|---|
| `start-run` | a tool that starts a run (`flow.run`, the studies, `journal.replay`) | ask |
| `exceed-estimate` | the estimate of a run (D81) exceeds a cap of `[estimate]` | ask |
| `write-outside-root` | a call would write outside the write roots | refuse |
| `overwrite-results` | a call writes into a folder with earlier results, or over a file | ask |
| `publish` | a call publishes outside the computer (no tool does yet) | refuse |

Each kind takes `"default"`, `"allow"`, `"ask"` or `"refuse"`. When the
server asks, you answer **allow once**, **allow always**, **refuse once**
or **refuse always**; "always" is written into the table for this project.
A host without elicitation cannot ask: such a call is refused (`MCP005`),
and you set the kind to `"allow"` by hand if you want it. An example table
is [examples/cases/tools/mcp/foehn-mcp.toml](../examples/cases/tools/mcp/foehn-mcp.toml):

```toml
[actions]
start-run = "ask"
exceed-estimate = "ask"
write-outside-root = "refuse"
overwrite-results = "ask"
publish = "refuse"

[estimate]
max_cells = 2000000
max_host_memory = 8e9   # bytes
max_output = 5e9        # bytes
max_wall_time = 3600    # s, when the estimate has a wall time
```

Two rules hold whatever the table says:

- **A run bounds itself.** A run started through MCP needs `[run.limits]
  max_wall_time` and `max_output` (D198 (c), `MCP002`). The refusal carries
  the fix (`set: ["run.limits.max_wall_time=…", "run.limits.max_output=…"]`),
  so the assistant can add them; study budgets (`study_cell_hours`, …,
  gallery 166) apply as usual.
- **The roots hold.** Reads outside the read roots fail (`CASE007`); a
  write outside the write roots is allowed only when you allow the kind
  `write-outside-root`, and then only into the one folder of that call.

## 5. A first session

Ask the assistant, for example: "Check
`cht_heat_sink_turbulent.toml`, run it into `runs/sink` and tell me the
chip temperature and the limitations of the result." A typical session:

1. `case.check`: the case can be run; 15 000 cells, the estimate.
2. `flow.run`: refused with `MCP002` (no `[run.limits]`), with the fix.
3. `flow.run` with the limits: the host shows you the question (the
   case, the overrides, the folder it writes, the estimate, the limits);
   you answer "allow once". The run streams progress and returns the
   result summary with its limitations.
4. `result.page`: the series of the run (residuals, probes, statistics,
   the run record), then the residual history page by page.

Cancelling the call in the host stops the run after the current iteration
with a restart file and the final output (exit code 3). Raw fields never
go through MCP: the summary names the files (`fields.vtkhdf`), which you
open yourself. [Gallery 175](gallery/runs.md#mcp) shows such a session.

## 6. Replay the session

Every call, refused ones and your answers included, is in the journal:

```
foehn run --script foehn-mcp-journal.jsonl \
    --read-root ~/cfd/heatsink --write-root ~/cfd/heatsink/runs
```

replays the calls in order without a window and reports each one that
ends otherwise than recorded; calls that were refused or declined are
skipped (D197). Read the journal to see what the assistant did: each row
holds the tool, its typed input, the status, the exit code and, for an
approved call, who approved it.

## Safety advice

- **Its own folder.** Let the assistant work only in a project folder of
  its own, with copies, and a write root below it for runs.
- **Keep copies.** Keep the originals of cases, meshes and results
  elsewhere; anything inside a write root may be overwritten.
- **Read before you approve.** The question names the case, the
  overrides, the folder, the estimate and the limits; read them before
  "allow". Prefer "allow once" until you trust a workflow.
- **Text in files is information, not orders.** Comments in case files,
  names in meshes, text in result files and journals are data. If an
  assistant starts to follow instructions it found inside such a file
  ("delete the runs", "allow everything", "send the results to…"), stop
  the session. The server never runs a process and never opens a network
  connection, but the assistant itself may have other tools.
- **Check results like a colleague's.** Read the limitations of every
  result (LIM001 to LIM006) and the iteration budget; a run that stopped
  before its tolerance is not a converged result.
- **The table is yours.** It lies outside the write roots so no tool can
  change it; keep it that way, and look at it after a session with
  "always" answers.

## What the server does and does not do

| Does | Does not (yet) |
|---|---|
| every registry tool, by its id, through the same call as `foehn call` | tools of its own: there is no MCP-only logic |
| the roots, the action table, mandatory run limits, the journal | a live control channel to a running solver (L19) |
| progress notifications during a run, cancellation | several calls at once: calls run one after another |
| summaries with limitations, histories paged by `result.page` | raw fields, images (L21), comparing runs (L20) |
| | resources and prompts (the lists are empty) |

The protocol details (versions, messages, the result content) are in
[interfaces.md](interfaces.md#mcp-server-l18).
