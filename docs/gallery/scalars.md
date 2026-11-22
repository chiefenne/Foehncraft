# Passive scalars and the wall distance

Part of the [Example gallery](../gallery.md): passive scalars carried by the flow (D51) with the generic scalar transport equation of stage 6 (D116), the bounded BDF2 (D28 (1)), the bounded convection of the scalars and the energy equation (D237) and the exact wall distance (D48, D108).

| # | Example | Stage | Feature |
|---|---|---|---|
| 63 | [A heated channel and a bounded front](#passive-scalars) | 6 | `[scalars.<name>]`, `[boundary.<name>.scalar.<scalar>]`, `bounds`, `exact`, `wall-distance` |
| 211 | [A step that stays a step: the MLP vertex limiter](#mlp-limiter) | 7 | `limiter = "mlp"`, `limiter_function`, `limiter_epsilon`, `limiter_epsilon_scaling`, `limiter_relaxation`, `verification.scalar_range` |

---

<a id="passive-scalars"></a>
## 63. A heated channel and a bounded front

<img src="../figures/gallery/passive_scalars.svg" alt="Left: the temperature profile across the heated channel at x = 1 m, 64 x 32 cells against the exact quartic. Middle: the Nusselt number on three meshes converging to 140/17. Right: the head of a front of a bounded scalar, plain BDF2 overshooting the upper bound by 3.3e-4, the bounded BDF2 staying below it" width="800" height="298">

**The laminar heated channel.** `channel_heated.toml` carries a passive
temperature θ in the Poiseuille flow of `poiseuille_hex.toml` (channel
height h = 1 m, mean velocity 1 m/s, Péclet number 1000). Both walls take
the same heat flux, `flux = "q"` in their `[boundary.<wall>.scalar.theta]`
tables; the inlet takes the thermally developed profile. The exact
solution is θ = G x + h (2η³ − η⁴ − η) with η = y/h and the axial gradient
G = 2q/(U h); the run computes the L2 error against the `exact` expression
of the scalar and the observed order over three meshes:

| Cells | L2 error of θ (K) | Nu | Nu − 140/17 |
|---|---|---|---|
| 512 | 1.10e-3 | 8.3132 | +0.95 % |
| 2 048 | 2.74e-4 | 8.2536 | +0.22 % |
| 8 192 | 6.23e-5 | 8.2398 | +0.055 % |

The observed orders of the error are 2.01 and 2.14 (the case requires 1.8
to 2.2). The Nusselt number with the hydraulic diameter 2h,
Nu = 2 h q / (Γ (T̄_w − T_b)), uses the area mean of θ on the walls and
the mixed mean at the inlet from the report (`result.json`, per boundary
`mean`, `mixed_mean` and `flux`), recipe 63.

**A bounded front.** `scalar_front_bdf2.toml` lets a scalar with the
bounds [0, 1] enter a channel in a uniform flow, first-order upwind in
space, BDF2 in time, Courant number 0.8. BDF2 steps from the old value
φ* = (4 cⁿ − cⁿ⁻¹)/3, which exceeds 1 at the head of the front where c
still grows: without bounds (`--set 'scalars.c.bounds=[-1e9, 1e9]'`) the
solution overshoots to 1.00033. With `bounds = [0, 1]` the cells where φ*
leaves the bounds take a backward-Euler step for that step (one or two
cells per step here, counted in `history.csv` as `first_order_c`), and
the largest value is 1 − 7e-11. Smooth scalars keep second order in time:
the Taylor–Green vortex carrying a scalar gives the observed time orders
1.98 with BDF2 and 2.02 with Crank–Nicolson (`tests/scalars.rs`).

**The wall distance.** `[output] derived = ["wall-distance"]` writes the
exact distance of every cell centroid to the nearest wall face, found
with a bounding volume hierarchy over the wall triangles (binned surface
area heuristic, closest-point queries, ties broken by the global face
id, so the field is the same on any number of threads and partitions).
The turbulence models of stage 6 take it from the same search.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/channel_heated.toml --out output/heated
target/release/foehn flow examples/cases/flow/channel_heated.toml --out output/heated_mid \
    --set 'verification.levels=[1]' \
    --set 'output.lines=[{name="mid",from=[1.0,0.015625,0.05],to=[1.0,0.984375,0.05],points=32}]'
target/release/foehn flow examples/cases/flow/scalar_front_bdf2.toml --out output/front
target/release/foehn flow examples/cases/flow/scalar_front_bdf2.toml --out output/front_free \
    --set 'scalars.c.bounds=[-1e9, 1e9]'
python3 -m tools.plot.scalar_chart docs/figures/gallery/passive_scalars.svg \
    output/heated output/heated_mid output/front output/front_free
```

[↑ back to top](#passive-scalars-and-the-wall-distance)

---

<a id="mlp-limiter"></a>
## 211. A step that stays a step: the MLP vertex limiter

<img src="../figures/gallery/scalar_limiters.png" alt="Top: the scalar c of the oblique step on perturbed polyhedra with three limiters, a yellow region above a dark blue one along a line at 31 degrees; with Venkatakrishnan's limiter many cells along the step are outlined in red and cyan where c leaves the interval from 0 to 1, with Barth–Jespersen and with MLP almost none. Bottom left: the profile across the step at x = 0.8, the face limiters overshooting near the corners of the step, MLP staying within the bounds. Bottom right: the largest excursion beyond the bounds in percent of the range for the oblique step, the T-junction and the heat exchanger on a log scale, MLP below 0.25 % in all three, the face limiters up to 15.75 %" width="800" height="573">

**Why a face limiter is not enough.** Linear upwind with the limiter of
Venkatakrishnan or Barth–Jespersen keeps the reconstructed face values of
a cell inside the range of its face neighbours. On a multi-dimensional
mesh the face value that flows into a cell then comes from the range of
the upwind cell's neighbours, which are not all neighbours of the
receiving cell, so the converged solution can carry new extrema. On the
tetrahedra of the shell and tube heat exchanger (application 203) the hot
air reached 585.75 K behind a nozzle with an inlet temperature of
573.15 K.

**The vertex limiter.** `limiter = "mlp"` in `[solver.scalar]` or
`[solver.energy]` bounds the reconstruction at every vertex by the
smallest and largest value of all cells that share the vertex (the
multi-dimensional limiting process of Park & Kim 2012). The limiter value
is Nishikawa's R3 function (2022), smooth and third order in smooth
regions, with a threshold ε in the units of the field,
ε = c (φ_max − φ_min) (`limiter_epsilon` = c, 0.01): a new extremum
stays within about ε, c = 0 is strictly bounded. `limiter_relaxation`
(0.5) damps a two-iteration cycle of the limiter value that stalls a
steady run in a mixing layer.

**The oblique step.** `scalar_oblique_step_mlp.toml`: a uniform flow
(1, 0.6) m/s across the unit square, c = 1 entering through x = 0 and
c = 0 through y = 0, nearly no diffusion (Γ = 10⁻⁹), 576 polyhedra
perturbed by 30 %. The exact solution is a step along the line at 31°
from the corner.

| Limiter | c min | c max | Iterations |
|---|---|---|---|
| Venkatakrishnan (faces, the default) | −0.0175 | 1.0178 | 24 |
| Barth–Jespersen (faces) | 0 | 1 | 30 |
| MLP, R3, c = 0.01 | 0 | 1.0021 | 34 |
| MLP, R3, c = 0 | 0 | 1 | 39 |

Barth–Jespersen is bounded on this mesh, but on the tetrahedra of the
heat exchanger it is the worst of all (617.29 K). The excursions beyond
the bounds in percent of the range (bottom right of the figure):

| Case | Venkatakrishnan | Barth–Jespersen | MLP, c = 0.01 | MLP, c = 0 |
|---|---|---|---|---|
| oblique step, polyhedra | 1.78 | 0 | 0.21 | 0 |
| T-junction θ, coarse polyhedra (gallery 196) | 0.66 | 0.03 | 0.19 | stalls |
| heat exchanger, air inlet, tetrahedra (restart, 150 iterations) | 4.50 | 15.75 | 0.21 | 0.09 |

**What it costs.** On a smooth manufactured solution MLP keeps second
order; on the finest polyhedral mesh the L2 error with R3 is 2.5 times
the unlimited one (9.13·10⁻⁴ against 3.67·10⁻⁴; R5 1.6 times), on
tetrahedra 4 % more. ε = 0 and the per-cell ε
(`limiter_epsilon_scaling = "cell"`, shrinking with the mesh) can stall a
steady run: the T-junction does not converge with them from the medium
mesh on.

**Where it is used.** Venkatakrishnan stays the default. Of thirteen
heated and scalar cases on hexahedra, eleven give the same extremes with
both limiters to within 0.3 % of the range, and the two with an exact
smooth solution are better without the vertex bounds. MLP is set in the
cases on tetrahedra and converted polyhedra with sharp fronts or a
bounded quantity: the T-junction (`tjunction_mixing.toml`, where the
scalar now converges with the mesh, gallery 196), the two heat-exchanger
cases (`[solver.energy]`) and the tracer of the stirred tank.

| | |
|---|---|
| Code | `foehncraft_flow::limiter` (vertex stencil), kernels `VertexBounds` and `MlpLimiter` (`foehncraft_kernels::fv::mlp_limiter`, `nishikawa_psi`), D237 |
| Case file | `flow/scalar_oblique_step_mlp.toml`; fails unless c stays in [−0.005, 1.005] (`verification.scalar_range`) |
| Tests | `crates/foehncraft-flow/tests/scalars.rs` (`mlp_limiter_bounds_the_oblique_step`, `mlp_limiter_is_refused_for_the_momentum`) |
| Recipe | [211](../recipes/scalars.md#mlp-limiter) |

```
C=examples/cases/flow/scalar_oblique_step_mlp.toml
target/release/foehn flow $C --out output/limiters/mlp
target/release/foehn flow $C --out output/limiters/mlp0 --set solver.scalar.limiter_epsilon=0.0
target/release/foehn flow $C --out output/limiters/venk --set 'solver.scalar.limiter="venkatakrishnan"'
target/release/foehn flow $C --out output/limiters/bj --set 'solver.scalar.limiter="barth-jespersen"'
python3 -m tools.plot.scalar_limiters docs/figures/gallery output/limiters/venk output/limiters/bj \
    output/limiters/mlp output/limiters/mlp0
```

[↑ back to top](#passive-scalars-and-the-wall-distance)
