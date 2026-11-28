# Development: tests, CI and releases

Part of the [Recipes](../recipes.md): MPI tests, checking a change the way CI does, releases, derivative checks of the kernels, a shared install under /opt.

| # | Recipe |
|---|---|
| 4 | [Run the MPI test with several ranks on one machine](#mpi) |
| 5 | [Check a change the way CI does](#ci-local) |
| 6 | [Make a release (tag and GitHub release)](#release) |
| 25 | [Check the derivatives of the kernels (adjoint readiness)](#adjoint-ready) |
| 184 | [Check a turbulence model with a manufactured solution](#mms-turbulence) |
| 350 | [Install for all users of a Linux machine, several versions side by side](#shared-install) |

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

- **Pull requests** (gate): Linux only. The `ci` profile
  (release optimisation without LTO, results bitwise equal to `release`),
  all tests including the bitwise-reproducibility tests, the tests without
  HDF5 and METIS (a job of its own), the case catalogue
  in three parallel jobs, hypre, MPI, PTX and the Python package. A pull
  request that touches platform-sensitive files (`Cargo.toml`, `Cargo.lock`,
  `build.rs`, the toolchain, `.github/`, the command line, MPI, VTKHDF,
  restart, run control, console) also runs the tests, quick cases and MPI on
  macOS and Windows. Code jobs are skipped when only `docs/`, Markdown files,
  `assets/` or `LICENSE` change.
- **Main** (every push to main whose change is not docs-only): what the
  gate did not run: tests and quick cases on macOS and Windows, MPI on
  macOS, the long cases (cylinder at Re 100) and the tests with debug
  assertions on Linux (a job of its own, profile `ci-debug`: debug
  assertions and overflow checks at opt-level 1). A newer merge cancels the older main run; with squash merges the
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

The test job runs the tests with [cargo-nextest](https://nexte.st), in two
halves per platform side by side (`--partition count:1/2` and `2/2`): the
tests of all test programs share one pool (cargo test runs the programs one
after another), each in a process of its own. Both halves compile the same
build; the tests take about 30 CPU minutes. The job without HDF5 and METIS
runs its three crates in one nextest pool too. Locally
`cargo test --workspace --profile ci` runs the same tests;
`cargo nextest run --workspace --cargo-profile ci` runs them as CI does.
The three Linux case sets are balanced on their run times (about 8 minutes
each, `tools/check/ci_cases.py`).

The MPI job builds all its example programs in one `cargo build` with the
`ci` profile, so that the crates are compiled once with one feature set.

The integration tests of `foehncraft-mesh` are one test program
(`tests/main.rs`, one `mod` line per file), those of `foehncraft-flow` two
(`flow` in `tests/main.rs` and `flow2` in `tests/main2.rs`, each file in one
of them), not one program per file: each program compiled and linked the
solver's generic code again, and the test build was the longest step of
CI. Flow has two halves because the compiler front end works on one crate
with one thread: the single program compiled for nine minutes, the two
halves compile side by side. Run one file with a module filter, for example
`cargo test --release -p foehncraft-flow --test flow energy::` (or
`--test flow2` for a file listed in `tests/main2.rs`); a new test file needs
its `mod` line, in the shorter of the two flow lists
(`python3 -m tools.check.test_modules` checks it).
The G1 smoke run uses the `foehn` binary that the tests already built.

The verification status page and its chart
(`docs/verification/status.md`, `docs/figures/verification/grades.svg`) are
generated on main after each merge (workflow `docs-index.yml`), like the
index tables; pull requests commit only `features.toml` and `grades.json`.

[↑ back to top](#development-tests-ci-and-releases)

---

<a id="release"></a>
## 6. Make a release (tag and GitHub release)

A release is a Git tag `v<version>` (D46; `v0.1` to `v0.7`, then `v1.0.0`, `v1.0.1`, `v1.1.0` and `v1.2.0`) plus a GitHub release. Steps:

1. A pull request runs `python -m tools.release.changelog` (moves the
   fragments of `changelog.d/` into `## Unreleased`), then moves that section to
   `## <version> (YYYY-MM-DD) – ...` and adds a row to the summary table
   (`## Summary`, at the top of `CHANGELOG.md` below `## Unreleased`).
2. After the merge: GitHub, Actions, **Release**, *Run workflow*, version
   `1.2.0` (the number without the `v`), first with *dry run* checked (prints the notes), then without.

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
cargo test --release -p foehncraft-flow --test flow adjoint:: -- --nocapture
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

---

<a id="shared-install"></a>
## 350. Install for all users of a Linux machine, several versions side by side

As root, the release programs with a Python venv per version under
`/opt/Foehncraft/<version>/` (details: [Shared install on Linux](../install.md#shared-install-on-linux)):

```
curl -fsSLO https://raw.githubusercontent.com/chiefenne/Foehncraft/main/tools/install-shared.sh
bash install-shared.sh --version 1.2.0 --with-source
bash install-shared.sh --version 1.3.0 --with-source --no-current   # next to it, default stays 1.2.0
bash install-shared.sh --list
bash install-shared.sh --activate 1.3.0                              # new default
ln -s /opt/Foehncraft/use.sh /etc/profile.d/foehncraft.sh            # default on at login
```

A source build for this machine instead (Debian 13; Rust goes to
`/opt/Foehncraft/rust`, not into root's home):

```
apt install build-essential cmake git curl python3 python3-venv libglu1-mesa libxcursor1 libxft2 libxinerama1 libgomp1
bash install-shared.sh --version 1.2.0 --build --install-rust --native
```

Users:

```
source /opt/Foehncraft/use.sh
foehncraft-use              # the installed versions
foehncraft-use 1.3.0        # this shell runs 1.3.0
python3 $FOEHNCRAFT_SRC/examples/meshes/mixer.py --out output/mixer_mesh
# or: module use /opt/Foehncraft/modulefiles; module load foehncraft/1.3.0
```

`--dry-run` prints every step without changing anything; `--python system`
uses the distribution's `python3-*` packages instead of a venv.

[↑ back to top](#development-tests-ci-and-releases)
