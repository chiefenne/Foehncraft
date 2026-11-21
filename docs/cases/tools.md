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
