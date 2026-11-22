# Console log: review of the redesign and four additions (D214)

Status: reviewed 2026-11-22. Docs only; the build of the additions waits for the maintainer's go.

## 1. What is built

The redesign proposed on 2026-11-16 (proposal v3) was decided the same day and merged (#50 to #55, D125 to D128): frames of asterisks (or box characters with `[log] style = "box"`), the banner, header blocks for files, run control, mesh, fluid properties, solution settings as a matrix per equation and boundary conditions, tables with repeating headers, summaries, the closing frame with the exit code, `%FOEHN-I-<FILE>` message tags, `[log] columns` with extremes and cell ids, `run.jsonl` 1.1 (`columns`, `table_row`, `section` rows), the case format version with `foehn case upgrade`, and `[[schedule]]`. Since then the closing frame also carries the limitations (D199, "Run finished with N limitations"), and runs have iteration budgets, a diagnosis, the MCP server and the live control port (L19, events RUN022 to RUN024).

Rule that stays: the console text is a view. `run.jsonl` is the structured record, `result.json` and the other machine files change only additively, and no function exists only in the text (agents-first rule, D77 to D87, L18 to L21).

## 2. Decisions of the maintainer (2026-11-22)

1. **Detail level.** New `[log] detail`: `full` (default, the log as it is), `steps` (a logged time step shows its box and the step row, not the table of every outer iteration), `summary` (header blocks, summaries and the closing frame only). `run.jsonl` keeps the full data at every level. Meant for long transient runs and for agents.
2. **Columns change during a run.** `log.columns` becomes a runtime key: it can be set by a `[[schedule]]` entry and by `run.set` of the control port. The header row is printed again, a `columns` row goes to `run.jsonl`, and the change is an event with its equivalent schedule entry (as for every `run.set`). `run.limits.*` stay out of reach.
3. **Limitations up front.** A "Limitations" block after the settings blocks lists the limitations known before the run (graded features, one mesh, and so on); the closing frame keeps its list, which may add those only known at the end. `result.json` is unchanged.
4. **Events summary.** A block "Events" before the closing frame lists every change made during the run with step or iteration, time and code: schedule firings (RUN015), control-port changes and rejections (RUN023, RUN024), stops (RUN022), thermal acceleration events (ACC004), diagnoses. The rows stay in `run.jsonl`; the block lets a steered run be reproduced from its log.

## 3. Not decided, small, done with the build

- A UTF-8 check: the box style falls back to the asterisk style on a terminal that is not UTF-8 (the asterisk style and the banner are plain ASCII already).
- The proposal's `step_summary_every` and `header_every` stay as built (`header_every` exists; the step summary is the ledger line).

## 4. Build

One PR, after the current code threads finish and on the maintainer's go: `[log] detail`, the schedulable `log.columns`, the two blocks, tests that default logs stay byte-identical, gallery entry (extends 27), recipe, keys and schema.
