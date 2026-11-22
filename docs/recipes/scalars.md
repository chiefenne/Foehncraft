# Passive scalars and the wall distance

Part of the [Recipes](../recipes.md): passive scalars carried by the flow (D51, D116), their boundary conditions, bounds and exact solutions, their bounded convection (D237), and the exact wall distance (D48). The keys are on the [key page](../cases/flow-keys.md#passive-scalars-d51-d116).

| # | Recipe |
|---|---|
| 63 | [Carry a passive scalar and compute a Nusselt number](#passive-scalars) |
| 211 | [Keep a scalar or the temperature within its bounds on tetrahedra and polyhedra](#mlp-limiter) |

---

<a id="passive-scalars"></a>
## 63. Carry a passive scalar and compute a Nusselt number

**A scalar.** Each `[scalars.<name>]` table adds one transport equation,
ρ ∂φ/∂t + ∇·(ṁ φ) − φ ∇·ṁ − ∇·(ρΓ ∇φ) = ρ S, solved after the pressure
correction of every outer iteration with the current mass fluxes:

```toml
[scalars.theta]
diffusivity = 1e-3         # Γ (m²/s): ν/Pr, ν/Sc; a number or an expression of x, y, z
initial = 0.0              # start value, also an expression
source = 0.0               # S per unit mass, expression of x, y, z (and t)
bounds = [0.0, 1.0]        # optional: bounded BDF2, extremes in the log
exact = "G * x + ..."      # optional: L2 error, observed order of a sweep
```

A name is letters, digits and `_`; the names of the output fields, of the
derived fields and of the expression language are refused, among them
`T`, `k`, `omega` and `Y`, which belong to the energy equation, the
turbulence models and the species.

**Boundary conditions.** Every boundary kind has a default: zero flux on
walls, slip walls, symmetry planes and wedges, the value 0 at inlets, zero
gradient at outlets (no diffusion across them). A
`[boundary.<boundary>.scalar.<scalar>]` table changes it with exactly one
of the two keys:

```toml
[boundary.inlet.scalar.theta]
value = "h * (2 * (y / h)^3 - (y / h)^4 - y / h)"   # inlets and walls
[boundary.wall.scalar.theta]
flux = "q"          # Γ ∂φ/∂n into the fluid; for a temperature q_w/(ρ c_p) in K m/s
```

A flux at an inlet, anything at an outlet, an unknown scalar or boundary,
or both keys at once is `SCAL001`.

**Settings.** `[solver.scalar]` takes the convection scheme
(`linear-upwind` by default, `upwind`, `central`), the limiter, the
under-relaxation (0.9 steady, 1 transient), the tolerance and the
BiCGStab iterations of each solve, and `sweeps`, the solves per outer
iteration. Transient runs integrate the scalars with BDF2 (also under
Crank–Nicolson, whose scalars would not stay bounded); with `bounds` the
cells where the BDF2 old value leaves them take a first-order step
(`first_order_<name>` in `history.csv`).

**Output.** `residuals.csv` and `history.csv` get the columns
`scalar_<name>` (and the BiCGStab iterations or the first-order cells),
the console table one `Res<name>` column per scalar (`residual.scalars`
in `[log] columns`), `fields.vtkhdf` a dataset `<name>` and, with
`exact`, `<name>_error`, `line_<name>.csv` a column with the cell value,
the restart file the scalar and its old level. The report lists per
scalar its extremes, the errors and, per boundary, the area mean, the
mixed mean (weighted with the mass flux) and the mean diffusive flux.

**The Nusselt number of the heated channel.** Run
`channel_heated.toml` and take the means from `result.json`:

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/channel_heated.toml --out output/heated
python3 - <<'PY'
import json
r = json.load(open("output/heated/result.json"))
q, gamma, G, L = 1e-3, 1e-3, 2e-3, 2.0
for s in r["solves"]:
    b = {x["boundary"]: x for x in s["scalars"][0]["boundaries"]}
    tw = 0.5 * (b["ymin"]["mean"] + b["ymax"]["mean"])
    tb = b["xmin"]["mixed_mean"] + G * L / 2      # bulk value at x = L/2
    print(s["cells"], 2 * q / (gamma * (tw - tb)))  # Nu with D_h = 2 h; exact 140/17
PY
```

**The wall distance.** `[output] derived = ["wall-distance"]` writes the
distance of every cell centroid to the nearest face of the `wall`
boundaries (dataset `wall_distance`, m); a case without a wall is
refused. `foehn post derive CASE --fields wall-distance` writes it from a
restart file.

[↑ back to top](#passive-scalars-and-the-wall-distance)

---

<a id="mlp-limiter"></a>
## 211. Keep a scalar or the temperature within its bounds on tetrahedra and polyhedra

**When.** A scalar or the temperature has a sharp front (an inlet value
meeting another, little diffusion) on tetrahedra or converted polyhedra,
or it is a bounded quantity (a mixture fraction, a tracer) whose
overshoot matters. The face limiters (`venkatakrishnan`, the default,
and `barth-jespersen`) can leave the bounds there by several percent of
the range; on hexahedra aligned with the flow they rarely do, and a
smooth solution is more accurate without the vertex bounds.

```toml
[solver.scalar]                # or [solver.energy] for the temperature
limiter = "mlp"                # bounds at the vertices (Park & Kim 2012)
limiter_function = "r3"        # Nishikawa's R3 (default); r4, r5 are less dissipative
limiter_epsilon = 0.01         # ε = 0.01 (φ_max − φ_min): overshoot within about ε
limiter_relaxation = 0.5       # damps a limit cycle of a steady run
```

`limiter = "mlp"` is refused for the momentum and the turbulence
equations (`CASE004`). Check the extremes in the report (`flow.json`,
per scalar `min` and `max`; the energy equation as `temperature`) or let
the case fail outside a range:

```toml
[verification]
scalar_range = [-0.005, 1.005]
```

**If a steady run stalls.** The residual of the scalar stops at 10⁻⁴ to
10⁻³ while the flow has converged: the limiter value of a few cells
switches back and forth. Keep `limiter_epsilon` above 0 (0.01 is
enough), lower `limiter_relaxation` to 0.25, and keep
`limiter_epsilon_scaling = "domain"`: the per-cell ε (`"cell"`) shrinks
with the mesh and stalls like ε = 0 on fine meshes. Strict bounds
(`limiter_epsilon = 0`) suit transient runs.

**The oblique step** (gallery 211) with the four limiters:

```
cargo build --release -p foehncraft-cli
C=examples/cases/flow/scalar_oblique_step_mlp.toml
target/release/foehn flow $C --out output/limiters/mlp
target/release/foehn flow $C --out output/limiters/mlp0 --set solver.scalar.limiter_epsilon=0.0
target/release/foehn flow $C --out output/limiters/venk --set 'solver.scalar.limiter="venkatakrishnan"'
target/release/foehn flow $C --out output/limiters/bj --set 'solver.scalar.limiter="barth-jespersen"'
python3 - <<'PY'
import json
for v in ["venk", "bj", "mlp", "mlp0"]:
    s = json.load(open(f"output/limiters/{v}/flow.json"))["solves"][0]
    c = s["scalars"][0]
    print(f"{v:5s} c in [{c['min']:.4f}, {c['max']:.4f}], {s['iterations']} iterations")
PY
```

[↑ back to top](#passive-scalars-and-the-wall-distance)
