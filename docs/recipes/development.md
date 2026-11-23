# Development: tests, CI and releases

Part of the [Recipes](../recipes.md): MPI tests, checking a change the way CI does, releases, derivative checks of the kernels.

| # | Recipe |
|---|---|
| 4 | [Run the MPI test with several ranks on one machine](#mpi) |
| 5 | [Check a change the way CI does](#ci-local) |
| 6 | [Make a release (tag and GitHub release)](#release) |
| 25 | [Check the derivatives of the kernels (adjoint readiness)](#adjoint-ready) |
| 184 | [Check a turbulence model with a manufactured solution](#mms-turbulence) |

---

<a id="mpi"></a>
## 4. Run the MPI test with several ranks on one machine

Needs an MPI library and libclang (Ubuntu:
`sudo apt install libopenmpi-dev openmpi-bin libclang-dev`; macOS:
`brew install open-mpi`).

```
cargo test -p foehncraft-parallel --features mpi
cargo build --release -p foehncraft-parallel --features mpi --example mpi_smoke
mpirun --oversubscribe -n 4 target/release/examples/mpi_smoke
```

Prints `MPI smoke test passed` after a ring exchange (the step of a halo
exchange) and a reproducible sum that gives the same bits on every rank.

[↑ back to top](#development-tests-ci-and-releases)

---

<a id="ci-local"></a>
## 5. Check a change the way CI does

```
cargo fmt --all --check
python3 -m tools.check.ci_cases --check                               # every flow case is in a CI set
cargo clippy --workspace --all-targets -- -D warnings
cargo clippy --workspace --all-targets --features foehncraft-cuda/cuda,foehncraft-cli/cuda -- -D warnings
python3 -m tools.check.backend                          # D60: CUDA only in crates/foehncraft-cuda and gpu/
cargo test --workspace --profile ci                     # what a pull request runs (--release nightly)
cargo test -p foehncraft-mesh --no-default-features      # without HDF5 and METIS
cargo test -p foehncraft-fv --no-default-features
python3 gpu/ptx-kernels/build.py && python3 -m tools.check.kernels     # PTX, and R10: no atomics, no f32, no FMA
cargo build --profile ci -p foehncraft-cli
for s in linux-1 linux-2 linux-3; do                    # add "long" for the cylinder at Re 100
  python3 -m tools.check.ci_cases "$s" target/ci/foehn || break
done
pip install -r tools/requirements-ci.txt
python3 -m tools.check.vtkhdf output/*/mesh.vtkhdf output/*/fields.vtkhdf output/*/*/fields.vtkhdf
```

The workflow is `.github/workflows/ci.yml`. It has three tiers, so that
pull requests never wait for macOS runners (GitHub Free runs 20 jobs and
only 5 macOS jobs at once per account):

- **Pull requests** (gate, about 5 minutes): Linux only. The `ci` profile
  (release optimisation without LTO, results bitwise equal to `release`),
  all tests including the bitwise-reproducibility tests, the case catalogue
  in three parallel jobs, hypre, MPI, PTX and the Python package. A pull
  request that touches platform-sensitive files (`Cargo.toml`, `Cargo.lock`,
  `build.rs`, the toolchain, `.github/`, the command line, MPI, VTKHDF,
  restart, run control, console) also runs the tests, quick cases and MPI on
  macOS and Windows. Code jobs are skipped when only `docs/`, Markdown files,
  `assets/` or `LICENSE` change.
- **Main** (every push to main whose change is not docs-only): what the
  gate did not run: tests and quick cases on macOS and Windows, MPI on
  macOS, the long cases (cylinder at Re 100) and the debug-build tests on
  Linux. A newer merge cancels the older main run; with squash merges the
  newer commit contains the older one. When `Cargo.lock` or the toolchain
  changes, main also builds the dependencies in the `ci` profile, so that
  pull requests start from a warm cache. After a merge, watch this run and
  fix forward at once if macOS or Windows fail.
- **Full tier**: every night at 02:07 Vienna time (01:07 UTC), a pull
  request with the label `full-ci`, or by hand (Actions, CI, *Run
  workflow*): the `release` profile, every job on every platform, debug
  tests, the long cases, and in nightly and manual runs also the tests
  marked slow (`--ignored`).

The check "CI OK" sums up all jobs.

[↑ back to top](#development-tests-ci-and-releases)

---

<a id="release"></a>
## 6. Make a release (tag and GitHub release)

A release is a Git tag `v<version>` (D46; `v0.1` to `v0.7`, then `v1.0.0`) plus a GitHub release. Steps:

1. A pull request runs `python -m tools.release.changelog` (moves the
   fragments of `changelog.d/` into `## Unreleased`), then moves that section to
   `## <version> (YYYY-MM-DD) – ...` and adds a row to the summary table
   (`## Summary`, at the top of `CHANGELOG.md` below `## Unreleased`).
2. After the merge: GitHub, Actions, **Release**, *Run workflow*, version
   `1.0.1` (the number without the `v`), first with *dry run* checked (prints the notes), then without.

The workflow refuses a version without a changelog section, and creates
only the release if the tag already exists. The assistant starts it through
the GitHub API in the same way.

The release badge of the README updates by itself: the last step of the
workflow waits until img.shields.io reports the new tag and then purges
GitHub's cached copies of the README images (GitHub serves them through its
image proxy, which otherwise keeps the old badge for hours or days). If the
badge still shows the old version, run that purge by hand:

```
curl -s https://github.com/chiefenne/Foehncraft | grep -oE 'https://camo\.githubusercontent\.com/[0-9a-f]+/[0-9a-f]+' | sort -u | xargs -n1 curl -s -X PURGE
```

---

<a id="adjoint-ready"></a>
## 25. Check the derivatives of the kernels (adjoint readiness)

The kernels are prepared for a discrete adjoint solver (D59): they run on
dual numbers, matrices apply transposed, and transient runs reverse with
binomial checkpointing. The checks are tests; from the repository root:

```
cargo test --release -p foehncraft-flow --test adjoint -- --nocapture
cargo test --release -p foehncraft-flow --lib -- adjoint_tests checkpoint
cargo test --release -p foehncraft-mesh shape_derivatives -- --nocapture
```

The first prints one line per kernel, for example
`momentum_face: 101700 derivatives, largest relative difference 5.28e-8`.

**Derivatives of a kernel in your own code.** Instantiate the body on
`foehncraft_core::Dual<N>` instead of `f64`: give the inputs you
differentiate by the derivative 1 in one of the N directions
(`Dual::variable(x, k)`), all others as constants (`Dual::constant(x)`),
and read the N derivatives of each result from its `.d` (the value is
`.v`, the same bits as the `f64` result). Leave the limiter array ψ
constant for the frozen-limiter linearisation; switches such as the upwind
direction follow the values.

**Transposed matrices** on a flow device: `transpose(val, out)` writes Aᵀ
in the storage of A (for the solvers and the multigrid),
`transposed_product(val, x, out)` computes Aᵀx without storing it.

**Reversing a transient run**: implement `foehncraft_flow::checkpoint::Primal`
(or use `FlowRun` for the flow solver with a fixed time step) and call
`checkpoint::reverse(run, steps, snapshots, adjoint)`; the adjoint step i
is called for i = n − 1, …, 0 with level i current. `schedule(n, s)` lists
the actions, `python3 -m tools.plot.checkpoint_chart out.svg` draws the cost
against the snapshots.

[↑ back to top](#development-tests-ci-and-releases)

---

<a id="mms-turbulence"></a>
## 184. Check a turbulence model with a manufactured solution

Run the order checks of every turbulence model (SST, k-ε, realizable
k-ε, k-ε-ζ-f, k-ε-ζ-f₀, elliptic blending, the curvature correction and
the QCR stress); they take under a minute in a release build and print the
errors per level and the observed orders:

```
cargo test --release -p foehncraft-flow --lib turbulence::mms -- --nocapture
```

To add a model: write a `Model` in `crates/foehncraft-flow/src/turbulence/mms.rs`
with its fields and its equations at a point (the published equations,
not the kernel code), list it in `model_of`, and add a test that calls
`check(model, mesh, band)`. Keep every switch of the model on one branch
in the manufactured fields (`manufactured_fields_stay_on_one_branch`);
a continuous clip, as of f_r1, is allowed. A term that needs a different
flow selects it with `Model::flow`. The errors of the three levels give
p = ln(e₁/e₂)/ln(h₁/h₂); a wrong constant in the model equations drops p
to 0.2 to 0.7. The figure of gallery 184 comes from
`python3 -m tools.plot.mms_turbulence docs/figures/gallery/mms_turbulence.png`.

[↑ back to top](#development-tests-ci-and-releases)
