# Verification

Part of the [Recipes](../recipes.md): the discretisation uncertainty of a sweep (D67) and the generated verification manual (D68).

| # | Recipe |
|---|---|
| 34 | [Estimate the discretisation uncertainty of a mesh or time-step sweep](#gci) |
| 35 | [Regenerate the verification manual](#manual) |
| 165 | [Grade a feature on the verification status page](#status) |
| 170 | [Read and act on the limitations of a run](#limitations) |
| 271 | [Evaluate the steady residual and check the reverse sweep of the adjoint](#adjoint-a0) |

---

<a id="gci"></a>
## 34. Estimate the discretisation uncertainty of a mesh or time-step sweep

A sweep of three or more meshes (`[verification] levels`, the refinement
levels of `[mesh.generate]` or `[mesh.cylinder]`) or of three or more
time steps (`[verification] time_steps`) reports the uncertainty of its
monitored quantities after Celik et al. (2008). No exact solution is
needed: a case without `[analytic]` may have levels too (an
`expected_order` still needs one, it checks the error norms).

```toml
[mesh.cylinder]
refine = 0                     # the coarsest level

[verification]
levels = [0, 1, 2]             # 6 144, 24 576 and 98 304 cells
benchmark = "schafer-turek-2d1"
benchmark_tolerance = 0.08
```

```
target/release/foehn flow examples/cases/flow/cylinder_re20_gci.toml
grep -A 25 "Discretisation uncertainty" output/cylinder_re20_gci/flow.md
python3 -c "import json; [print(e['quantity'], e['order'], e['extrapolated'], e['band'], e['verdict']) for e in json.load(open('output/cylinder_re20_gci/result.json'))['uncertainty']]"
```

From Python:

```python
from foehncraft import run
res = run("examples/cases/flow/couette_viscosity.toml")
e = res.uncertainty["boundaries.ymin.force_x"]
print(e["order"], e["extrapolated"], e["band"], e["gci"], e["verdict"])
```

**Reading the table.** The cell size is h ∝ N^(−1/d) (N cells, d = 2 for
a one-layer mesh) or the time step; φ1 is the finest solution. `Order`
is the observed order p of the three solutions, `Extrapolated` the
Richardson value φ1 + (φ1 − φ2)/(r^p − 1), `± band` the uncertainty
U = 1.25 |φ1 − φ2| / (r^p − 1) of φ1 in its own unit and `GCI` the band
relative to |φ1| (left out where |φ1| is smaller than its band, for
example the mean pressure of a cavity, whose level is arbitrary). The
verdict:

| Verdict | Ratio R = (φ1 − φ2)/(φ2 − φ3) | Meaning |
|---|---|---|
| `monotonic` | 0 < R < 1 | the solutions converge from one side: order, extrapolation and band are meaningful |
| `oscillatory` | −1 < R < 0 | they converge and change direction: the order is that of the magnitudes, the band an estimate |
| `divergent` | \|R\| ≥ 1 | the differences do not fall: no order, refine further or look for a singularity |
| `unchanged` | φ1 = φ2 = φ3 to round-off | nothing to estimate (a quantity the scheme reproduces exactly) |

Quantities at a singularity converge slowly or not at all: the force on
the lid of the lid-driven cavity and its viscous dissipation grow with
every refinement (the shear stress is singular in the two upper corners),
and the report says so with an order near zero. The error norms against
an exact solution are not in the table; their observed orders are
reported as before (`orders`, `time_orders`).

Keys: `verification.levels`, `verification.time_steps`,
`verification.benchmark`, `verification.reference`.

<a id="manual"></a>
## 35. Regenerate the verification manual

The manual ([docs/verification/manual.md](../verification/manual.md),
D68) is generated from the case files and the reports of a catalogue run:
every check of every case with its band and its value, the uncertainty
tables of the sweeps and a chart of the observed orders. Regenerate it
for a release and after a change of a scheme:

```
cargo build --release -p foehncraft-cli
for s in linux-1 linux-2 linux-3 long; do python3 -m tools.check.ci_cases $s target/release/foehn; done
python3 tools/release/verification_manual.py output
git diff --stat docs/verification docs/figures/verification
```

The script exits with 1 when a case failed and lists the cases without a
report (for example when the `long` set was not run). It writes no dates
and no timings, so the same reports give the same files.

---

<a id="status"></a>
## 165. Grade a feature on the verification status page

The status page ([docs/verification/status.md](../verification/status.md),
D199) gives every feature a grade, verified, partial or unverified,
computed from the feature table
[docs/verification/features.toml](../verification/features.toml) and the
last verification manual. The grade is never written by hand: add the
evidence and the gaps, and the script computes it.

A new feature, or new evidence for one, is a row or an item in the table:

```toml
[[feature]]
id = "bc.mass-flow-inlet"
area = "Boundary conditions"
title = "Mass-flow inlet"
switch = 'type = "mass-flow-inlet"'
evidence = [
  # kind: exact | manufactured | reference-data | consistency
  { kind = "exact", case = "poiseuille_mass_flow", what = "Poiseuille flow from its mass flow, second order" },
  { kind = "consistency", test = "crates/foehncraft-flow/tests/bcs.rs::checks_of_the_new_conditions", what = "input checks" },
]
gaps = []                       # every missing check or band not met
range = { reynolds = [0.0, 2300.0] }
```

`case` is a case of the catalogue, `test` a Rust `#[test]` or pytest
function as `path::function`, `run` a recorded run outside CI with the
`doc` page that records it. A model (`model = true`) needs reference data
to be verified. Then regenerate and check:

```
python3 -m tools.release.verification_status          # writes status.md and grades.json
python3 -m tools.release.verification_status --check  # what CI runs
git diff docs/verification/status.md
```

The check refuses a case, test, page or anchor that does not exist, a
catalogue case that checks no feature, and a page that is out of date.
`tools/release/verification_manual.py` regenerates the page as well, so a
new manual run updates the grades (a case that failed makes its features
unverified).

---

<a id="limitations"></a>
## 170. Read and act on the limitations of a run

Every run lists what to know before trusting its result (D199). They are
in the last frame of the log ("Run finished with N limitations"), at the
top of `flow.md`, and as data:

```
foehn --json --quiet flow examples/cases/flow/pipe_gmsh_raw.toml | python3 -c "import json,sys; [print(l['code'], l['message']) for l in json.load(sys.stdin)['limitations']]"
```

```python
from foehncraft import run
res = run("examples/cases/flow/pipe_gmsh_raw.toml")
for lim in res.limitations:
    print(lim["code"], lim.get("feature"), lim["message"])
```

What each code asks for:

| Code | Means | What to do |
|---|---|---|
| LIM001 | a feature in use is graded partial or unverified | read its row on the [verification status](../verification/status.md) page; check the result against your own data, or add a case that raises the grade (recipe 165) |
| LIM002 | outside the validated range (first-cell y⁺, Re, Pr, Courant number, skewness) | change the mesh or the model so the run is inside the range, or treat the result as an extrapolation |
| LIM003 | a steady solve stopped before its tolerance | raise `solver.max_iterations`, check the monitors (AGENTS.md: diagnose at 2 000 to 3 000 iterations) |
| LIM004 | one mesh only | run a level sweep with `[verification] levels` for the GCI (recipe 34) |
| LIM005 | a GPU that is not verified | run the GPU checks on it (`tools/notebooks/stage6_colab.ipynb`), then add it to `VERIFIED_DEVICES` |
| LIM006 | automatic remediation of poor cells | improve the mesh (optimised tetrahedra, prisms at walls); the remediated faces are first order |

Limitations never change the exit code (0 stays 0), so scripts that test
the exit code keep working; an agent reads the list.

<a id="adjoint-a0"></a>
## 271. Evaluate the steady residual and check the reverse sweep of the adjoint

Stage A0 of the adjoint (L17, gallery 271). The residual of a steady
laminar run, from Rust:

```rust
use foehncraft_flow::{case::FlowCase, residual::steady_flow};
let case = FlowCase::from_file("examples/cases/flow/adjoint_residual_cavity.toml".as_ref())?;
let mut fl = steady_flow(&case, Some(1))?;
for _ in 0..300 {
    fl.iterate()?;
}
let r = fl.steady_residual()?; // R_u, R_c, R_ṁ per cell and face
println!("{:?}", r.norms());   // relative norms; the run is not changed
```

Turbulence, transported scalars, a variable density, time steps and
several ranks are refused with CASE005 (they come with later stages).

The residual history of both algorithms and the picture of gallery 271:

```sh
cargo run --release -p foehncraft-flow --example adjoint_a0 -- \
    examples/cases/flow/adjoint_residual_cavity.toml output/adjoint_a0
python3 -m tools.plot.adjoint_a0 docs/figures/gallery/adjoint_a0_residual.png output/adjoint_a0
```

The checks of the residual and of the reverse sweep (dot product, central
differences, 1 against 4 threads, the cost table):

```sh
cargo test --release -p foehncraft-flow --lib steady_residual
cargo test --release -p foehncraft-flow --test adjoint reverse -- --nocapture
FOEHN_A0_CELLS=40,40,30 cargo test --release -p foehncraft-flow --test adjoint reverse_sweep_cost -- --nocapture
```

The compiler reverse mode (option B) needs a nightly toolchain and is not
part of the workspace: see `tools/spikes/l17-enzyme-momentum-face/README.md`.
