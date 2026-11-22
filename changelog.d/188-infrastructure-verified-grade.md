- **Infrastructure features can reach verified (D199).** Restart, logs and
  run record, Python package, command line, tool registry, journal and
  replay and the MCP server are graded through contract evidence: a row
  marked `infrastructure = true` lists its contracts (bitwise, schema,
  equivalence, replay, control, diagnostics, round-trip) and each needs a
  test of kind `contract`. Three new rows (tool registry, journal and
  replay, MCP); 48 of 97 features are verified, up from 41 of 93. New test:
  every row of `run.jsonl` follows its schema.
