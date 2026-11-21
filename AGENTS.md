# Rules for agents and contributors

This file holds the standing rules of the project (decision D200 in
[docs/requirements.md](docs/requirements.md)). It is written for coding
agents first and applies to every contributor. **MUST** rules are never
broken; **SHOULD** rules are followed unless the pull request says why not.
Lessons from earlier work (what not to try again) are in
[docs/pitfalls.md](docs/pitfalls.md); the documentation conventions in
detail are in [CONTRIBUTING.md](CONTRIBUTING.md).

## Scope and approval

- **MUST** start code changes only after the maintainer's go for the task.
  Documentation of decisions already taken needs no extra go.
- **MUST** ask before any further work on solver speed (coupled solver,
  linear solvers, multigrid): the speed work of step 22 is finished
  ([pitfalls](docs/pitfalls.md#coupled-solver)).
- **MUST** keep every function reachable without a window: scriptable,
  with structured (TOML, JSON) inputs and outputs, no interactive-only
  paths. The engine is driven by the CLI, the Python package, the MCP
  server and later the desktop app through the same calls (D86, D196).
- **MUST** add every engine function a front end calls as a tool of
  `foehncraft-tools` (D196) and call it from the front ends through the
  registry; `cargo xtask layering` fails on a direct engine call in the
  CLI or the Python module, and on a crate dependency against the layers.
  Regenerate `docs/tools.md` with `foehn commands --markdown`.
- **SHOULD** flag known flaws of industry behaviour when a feature mirrors
  it, so the maintainer can decide; do not copy flaws.

## Clean room and naming (D201)

- **MUST NOT** read copyleft CFD source code at all.
- **MAY** read permissively licensed source (MIT, BSD, Apache);
  **MUST NOT** copy it.
- **MUST** build methods from papers, documented file formats and
  published reference data, and cite them.
- **MUST** get a decision of its own before using numbers produced by other
  solvers (for example CFD planes in a paper, reference values of other
  codes).
- **MUST NOT** name commercial CFD solvers, or other codes we assess,
  anywhere: code, docs, commit messages, pull requests. Say "industry
  standards" or "modern development codes".
- **MUST** keep data under restrictive licences out of the repository;
  only the published numbers go in (see
  [docs/verification/DATA_SOURCES.md](docs/verification/DATA_SOURCES.md)).
- **MUST** say "reacting flows", never "combustion".

## Pull requests

- **MUST** name branches with keywords of the task
  (`docs-agents-guide-pitfalls`, not random characters).
- **MUST** fill in the pull request template, including **what ran and
  what did not, and why** (D203).
- **MUST** write the changelog entry as a fragment
  `changelog.d/<PR number>-<short-name>.md`, never in `CHANGELOG.md`
  ([changelog.d/README.md](changelog.d/README.md)).
- **MUST NOT** edit the "by number" tables of `docs/gallery.md` and
  `docs/recipes.md`; CI regenerates them after the merge.
- **MUST NOT** merge main into a pull request whose CI is green; merge main
  only for a real conflict.
- **MUST** squash-merge when CI is green and delete the branch; then watch
  the main run of the merge commit (macOS and Windows are tested only
  there) and fix forward at once when it fails.
- **SHOULD** open the pull request as a draft early, so its number is known
  for the changelog fragment and the gallery and recipe numbers.

## Every feature is documented (D18)

- **MUST** add, in the same pull request, a gallery example, case files in
  `examples/cases/` listed in the case catalogue, and a recipe; gallery and
  recipe entries take the pull request's number.
- **MUST** give every gallery example at least one image, with no
  exception.
- **MUST** give every model or physics feature 1 to 3 colourful field
  images (velocity, pressure, temperature, vorticity, a scalar) in its
  gallery entry. Core-solver topics (linear solvers, speed, CI) need none
  beyond the one image.
- **MUST** write images with their size (`width` and `height`); run
  `python -m tools.docs.images` to add them.
- **MUST** look at every rendered image, full size and as a thumbnail.
  **MUST NOT** publish images with four-fold or otherwise symbol-like
  patterns: rotors in published images have 3 or 5 blades (CONTRIBUTING.md,
  "Images must not resemble symbols").
- **MUST** leave the README ("What works today", status line, commands,
  crate table) and the docs describing main as it is after the merge.
- **MUST** update the "Implemented" note of the decision or stage row in
  `docs/requirements.md`.

## Code

- **MUST** keep results independent of the output schedule: an observer
  between two outer iterations (results, balances, monitors) restores the
  lagged state it found (`Flow::lagged` / `Flow::set_lagged`).
- **MUST** keep the bitwise-reproducibility tests passing (thread counts,
  AMG, CPU backends).
- **MUST** keep CUDA code in the CUDA backend crate (D60).
- **MUST** use `unsafe` only in the named GPU and SIMD modules, each block
  with a `// SAFETY:` comment (D202).
- **SHOULD** avoid `unwrap` and `panic!` in new code of the case,
  mesh-reader and engine crates and return a diagnostic code instead;
  fix old code when it is touched (D202).
- **MUST NOT** assert exact iteration counts of iterative solvers in tests;
  quote counts to about 10 %.
- **MUST** keep the exact reference bands of verification cases; do not
  widen a band to make a test pass.
- **MUST NOT** skip, disable or mark as ignored a test to get CI green.

## Runs

- **MUST** check a steady run after 2 000 to 3 000 iterations: what goes
  wrong, or which settings (relaxation, tolerance, mesh, model switches)
  must change. Never let a test run go to 10 000 or more iterations.
- **MUST** judge steady convergence by a monitored quantity, not by the
  residual alone ([pitfalls](docs/pitfalls.md#steady-runs)).
- **SHOULD** run a new case file from the command line before quoting its
  iteration count or run time, and measure run times on a quiet machine.

## Checks to run

Run these before pushing. They are the PR gate of
[.github/workflows/ci.yml](.github/workflows/ci.yml); CI runs them again.

Rust, for any code change:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked --features foehncraft-cuda/cuda,foehncraft-cli/cuda
cargo test --workspace --locked --profile ci
cargo run --locked --profile ci -p foehncraft-cli -- bench g1 examples/cases/bench/g1_ci_small.toml
cargo doc --workspace --no-deps --locked
```

Repository checks:

```sh
python3 -m tools.check.ci_cases --check   # every flow case is in a CI case set
python3 -m tools.check.backend            # CUDA only in the CUDA backend crate
cargo xtask layering                      # crate layering, front ends through foehncraft-tools
python3 -m tools.check.c_free             # build without C libraries
```

Documentation, for any docs change:

```sh
python3 -m tools.docs.indexes            # tables can be generated (no number twice)
git checkout docs/gallery.md docs/recipes.md docs/applications.md   # do not commit them
python3 -m tools.docs.images --check      # images carry their size
python3 -m tools.release.changelog --print > /dev/null   # fragments parse
```

On a pull request the index tables are not current yet, so
`tools.docs.indexes --check` fails there by design; it is meant for main.

After changing dependencies:

```sh
cargo update -w --manifest-path crates/foehncraft-py/Cargo.toml
```

Flow cases of a CI set run with
`python -m tools.check.ci_cases <set> target/ci/foehn`; the long sets run
on main and in the nightly run.
