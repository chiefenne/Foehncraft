# Boundary conditions of the flow solver

Part of the [Example gallery](../gallery.md): slip walls, axisymmetric wedges, inlet profiles from tables, backflow treatment at pressure outlets and the convective outlet (D29), periodic pairs with the heat-flux jump of the scalars (D29, D109 (1)), and absolute or gauge pressure values with an operating pressure per fluid part (D239).

| # | Example | Stage | Feature |
|---|---|---|---|
| 21 | [Boundary conditions: slip walls, wedges, tables, backflow, convective outlets](#boundary-conditions) | 5 | slip wall, axisymmetric wedge mesh, inlet profile from a table, backflow treatment at pressure outlets, convective outlet |
| 67 | [Periodic modules: a channel and a short heated pipe](#periodic) | 6 | `type = "periodic"`, `mass_flow`, `pressure_gradient`, `[boundary.<name>.scalar.<scalar>] periodic = "heat-flux"`, `bulk` |
| 220 | [Absolute pressures: two fluids at 1 bar and 1.5 bar](#absolute-pressure) | 7 | `[fluid] pressure_values`, `operating_pressure`, `[regions.<name>] operating_pressure`, case format 0.7 |

---

<a id="boundary-conditions"></a>
## 21. Boundary conditions: slip walls, wedges, tables, backflow, convective outlets

<img src="../figures/gallery/pipe_wedge_speed.png" alt="Hagen–Poiseuille flow on a wedge of 20°, coloured by speed" width="800" height="343">

**Axisymmetric flow on a wedge.** A 2D mesh of the half plane y ≥ 0 is
revolved about the x axis into a wedge one cell thick
(`mesh.generate.wedge`); the nodes on the axis merge, so the cells next to
it are prisms with an edge on the axis and the mesh has no axis boundary. The
picture is the pipe at 20° to show the shape; the verification runs at
0.5°, where Hagen–Poiseuille flow converges at second order (2.07 and 2.04
on 256 to 4 096 cells). Every cell type of the box generator revolves,
perturbed or not. The planar faces make a wedge of θ a slice of a
polygonal pipe of 360°/θ sides; its result differs from the round pipe by
about θ²/8 (9e-4 at 5°), so applications take 2–5° and a convergence
study 0.5–1°.

<img src="../figures/gallery/cylinder_outlets.svg" alt="Lift coefficient over the last second at Re = 100: long channel, short channel with a pressure outlet, short channel with a convective outlet" width="760" height="420">

**The outlet close behind a body.** The cylinder of [gallery 15](flow.md#cylinder) at Re = 100
with the outlet 5.5 D behind the cylinder instead of 20 D: the vortices
cross the outlet, so fluid enters through it in parts of every period.
A static pressure outlet with the default backflow treatment (the
prescribed pressure taken as the total pressure of the entering fluid)
gives the lift of the long channel on the same near-wake cells; a
convective outlet (∂u/∂t + U_c ∂u/∂n = 0, U_c the mean outflow speed)
gives C_L,max 1.6 % lower. Drag, Strouhal number and Δp stay inside the
reference ranges with both.

| Case | Cells | C_D,max | C_L,max | St | Δp (Pa) |
|---|---:|---:|---:|---:|---:|
| reference | | 3.22–3.24 | 0.99–1.01 | 0.295–0.305 | 2.46–2.50 |
| long channel ([gallery 15](flow.md#cylinder)) | 6 144 | 3.2319 | 0.9712 | 0.2955 | 2.4819 |
| short, pressure outlet, `backflow = "total-pressure"` | 4 928 | 3.2336 | 0.9713 | 0.2951 | 2.4825 |
| short, convective outlet | 4 928 | 3.2261 | 0.9555 | 0.2952 | 2.4779 |

| Check | Result |
|---|---|
| Slip walls, uniform flow on hex, tet, poly and mixed cells perturbed by 20 % | velocity error 1.5e-11, pressure 1.0e-11 (exact to the solver tolerance) |
| Inlet profile from a table (rows at the face centres) | the run with the exact parabola, velocity error 1.592e-3 for both |
| Reversed flow through a pressure outlet, `backflow = "total-pressure"` | uniform flow and the static pressure p − ½ρU² exact to 3e-11 |
| A side outlet above the channel pressure, `backflow = "suppress"` | the faces close, act as a slip wall, no mass crosses (exactly 0), uniform flow exact to 1e-11 |
| Convective outlet, no pressure boundary | outlet flow = inflow to round-off in every step; the start-up disturbance leaves (error 4e-5 after two flow-through times) |
| Convective outlet, restart and rejected steps | bit-identical continuation, the old boundary levels in the restart file |

| | |
|---|---|
| Code | `foehncraft_mesh::wedge`, `foehncraft_flow::bc` (`Backflow`, `ProfileTable`, the new condition types and checks), `simplec::Convective`, kernels `flow::close_backflow`, `UB_OUTLET`, `UB_CONVECTIVE`, `UB_NO_BACKFLOW`, `PB_OUTLET`, `PB_NO_BACKFLOW` |
| Tests | `tests/bcs.rs`, `tests/stage5.rs` (`convective_outlet_restarts_bit_identically`), `wedge.rs` (`wedges_of_every_type`) |
| Case files | `flow/channel_slip.toml`, `poiseuille_table.toml`, `pipe_wedge.toml`, `cylinder_re100_short.toml`, `cylinder_re100_convective.toml` |
| Recipe | [29](../recipes/boundary-conditions.md#boundary-conditions) |

```
target/release/foehn flow examples/cases/flow/pipe_wedge.toml
target/release/foehn flow examples/cases/flow/pipe_wedge.toml --set mesh.generate.wedge=20.0 --set 'verification.levels=[1]' --set mesh.generate.perturb=0.0 --set output.vtkhdf=true --out output/wedge20
python3 -m tools.mesh.render_mesh output/wedge20/fields.vtkhdf pipe_wedge.png --field speed --size 1400x600 --zoom 1.6
target/release/foehn flow examples/cases/flow/cylinder_re100.toml
target/release/foehn flow examples/cases/flow/cylinder_re100_short.toml
target/release/foehn flow examples/cases/flow/cylinder_re100_convective.toml
python3 -m tools.plot.outlet_chart cylinder_outlets.svg
```

[↑ back to top](#boundary-conditions-of-the-flow-solver)

---

<a id="periodic"></a>
## 67. Periodic modules: a channel and a short heated pipe

<img src="../figures/gallery/periodic.svg" alt="Left: the temperature across the short periodic heated pipe at x = 0.5 m, 16 x 16 cells against the exact profile. Middle: the Nusselt numbers of the periodic heated channel and pipe on three meshes approaching 140/17 and 48/11. Right: the driving force that holds the mass flow of the periodic channel and pipe, converging to the exact pressure gradients" width="800" height="298">

**One module instead of a long duct.** A fully developed flow repeats
itself along the duct, so one module with its inlet and outlet as a
periodic pair replaces the entry length. `poiseuille_periodic.toml` is
the channel of `poiseuille_hex.toml` with `xmin` and `xmax` as the pair,
driven by its mass flow:

```toml
[boundary.xmin]
type = "periodic"
partner = "xmax"
mass_flow = 0.1          # kg/s from xmin towards xmax
[boundary.xmax]
type = "periodic"
```

The faces of the two boundaries are matched by their centroids; each side
keeps its own faces, and the neighbour of a face is the cell behind its
partner, seen shifted by the period, in every face loop, gradient and
limiter. A uniform body force β along the translation drives the flow; the
solver adjusts it after every pressure correction so that the mass flow
through the pair is the prescribed one (or holds it fixed with
`pressure_gradient`). The velocity error falls at order 1.98 and 1.94 on
three perturbed hex meshes (the open channel: 1.98 and 1.89), and β
converges to the exact gradient 12 μU/h² = 1.2 Pa/m, +0.09 % on the
finest mesh. Driven by that gradient instead, the mass flow converges to
0.1 kg/s at second order (0.78 %, 0.20 %, 0.05 % on uniform meshes).

**The heat-flux jump** (D109 (1)). A temperature heated through the walls
is not periodic: it rises by Δ = Q/ṁ from one end of the module to the
other, Q the heat the walls put in. With `periodic = "heat-flux"` on one
side of the pair the scalar crosses it with that jump, computed from the
wall fluxes and the current mass flow in every iteration, and `bulk` fixes
its level (the mixing-cup value on that boundary):

```toml
[boundary.xmin.scalar.theta]
periodic = "heat-flux"
bulk = -0.24285714285714285   # −17/70 h at x = 0
```

`channel_heated_periodic.toml` is the heated channel of gallery 63 on one
module; the periodic solution is the exact one of the open channel:

| Case | Cells | L2 error of θ (K) | Nu | exact | Nu − exact |
|---|---|---|---|---|---|
| channel | 512 | 3.61e-3 | 8.4647 | 140/17 = 8.2353 | +2.79 % |
| | 2 048 | 6.87e-4 | 8.2793 | | +0.53 % |
| | 8 192 | 1.69e-4 | 8.2460 | | +0.13 % |
| pipe | 64 | 4.61e-3 | 4.5094 | 48/11 = 4.3636 | +3.34 % |
| | 256 | 9.38e-4 | 4.3919 | | +0.65 % |
| | 1 024 | 1.93e-4 | 4.3691 | | +0.13 % |

**The short heated pipe** (D111). `pipe_heated_periodic.toml` is a module
one diameter long on the axisymmetric wedge of `pipe_wedge.toml`: the
laminar pipe with a uniform wall heat flux, whose Nusselt number 48/11
otherwise needs a pipe of 50 diameters and more for the thermal entry
length. The observed orders of the temperature error are 2.40 and 2.02
(channel) and 2.30 and 2.28 (pipe). The pipe runs without a limiter on the
scalar: the limiter clips the smooth minimum of the temperature on the
axis and costs the order on the perturbed meshes (1.4 on the fourth
mesh). The driving force of the pipe converges to 8 μU/R² = 3.2 Pa/m.

| Check | Result |
|---|---|
| Pairing of perturbed hex, poly and tet boxes in x and z | every face matched, geometry unchanged to round-off |
| Periodic Poiseuille flow against the open channel, same meshes | velocity errors within a factor of two, order 1.98, 1.94 |
| Driving force at a fixed mass flow | 1.2011 Pa/m against 1.2 (+0.09 %) |
| Mass flow at a fixed pressure gradient, uniform meshes | errors 0.78 %, 0.20 %, 0.05 % (order 2.0) |
| Heat-flux jump and level | Δ = Q/ṁ = 4e-3 K and the bulk value to 1e-12 |
| Transient run with the heat-flux mode, restart after 6 of 12 steps | bit-identical, β in the restart file |
| Cases without periodic boundaries | results unchanged bit for bit |

| | |
|---|---|
| Code | `foehncraft_mesh::periodic` (`pair`), `Mesh::periodic`, kernels `fv::neighbour_centroid`, `fv::Boundary::across`, `foehncraft_flow::bc` (`BoundaryCondition::Periodic`, `periodic_pairs`, `Drive`), `simplec::Flow` (`body_force`, the mass-flow driving), `scalar` (`PeriodicMode`, the heat-flux jump and the level) |
| Tests | `tests/periodic.rs`, `periodic.rs` (`pairs_the_box_in_x_and_z`, `refuses_bad_pairs`) |
| Case files | `flow/poiseuille_periodic.toml`, `channel_heated_periodic.toml`, `pipe_heated_periodic.toml` |
| Recipe | [67](../recipes/boundary-conditions.md#periodic) |

```
target/release/foehn flow examples/cases/flow/poiseuille_periodic.toml --out output/pois_periodic
target/release/foehn flow examples/cases/flow/channel_heated_periodic.toml --out output/heated_periodic
target/release/foehn flow examples/cases/flow/pipe_heated_periodic.toml --out output/pipe_periodic
target/release/foehn flow examples/cases/flow/pipe_heated_periodic.toml --out output/pipe_mid \
    --set 'verification.levels=[1]' \
    --set 'output.lines=[{name="radius",from=[0.5,0.0078125,0.0],to=[0.5,0.4921875,0.0],points=32}]'
python3 -m tools.plot.periodic_chart docs/figures/gallery/periodic.svg \
    output/pois_periodic output/heated_periodic output/pipe_periodic output/pipe_mid
```

[↑ back to top](#boundary-conditions-of-the-flow-solver)

---

<a id="absolute-pressure"></a>
## 220. Absolute pressures: two fluids at 1 bar and 1.5 bar

<img src="../figures/gallery/pressure_values_two_parts.png" alt="Top: the absolute pressure in a parallel-plate counterflow exchanger 0.5 m long, the height stretched 25 times: the water channel below the steel plate is dark at 1.0 bar, the oil channel above it light green at 1.5 bar. Middle: the gauge pressure of each part relative to its own operating pressure, from 0 at the outlets to 3 Pa at the water inlet on the left and 45 Pa at the oil inlet on the right, white to red. Bottom: the gauge pressure along the middle of each channel, the water falling linearly from 3.1 Pa to 0 at x = 0.5 m, the oil rising from 0 at x = 0 to 45 Pa" width="1210" height="990">

**What the values mean.** From case format 0.7 every pressure value of
a case file is an absolute static pressure (`[fluid] pressure_values =
"absolute"`, the default, D239): an outlet to the atmosphere is
`pressure = 1.01325e5`, the water of a heat exchanger leaves at
`1.0e5` and the air at `1.5e5`. `pressure_values = "gauge"` keeps the
older meaning, values relative to the operating pressure,
p = p_op + p_gauge. Files of format 0.6 and older mean gauge: they read
with a note and keep their results bit for bit; `foehn case upgrade`
writes `pressure_values = "gauge"` into them.

**The operating pressure of every part.** The solver keeps working with the
gauge pressure p − p_op. With absolute values p_op is not a setting:
each fluid part (fluid cells connected without a solid between them,
D194) takes it from its own pressure boundaries, the static-pressure
outlets first (the smallest value), then the total-pressure outlets, then
the pressure inlets; a closed part from `[pressure] reference_value` or
`[initial] pressure`, else 101 325 Pa. The reader then writes every value
as its difference to the part's p_op (boundary values, expressions,
`[initial] pressure`, `[pressure] reference_value`) and sets
`[fluid] operating_pressure` and `[regions.<name>] operating_pressure`,
so that each gas's density follows the absolute pressure of its own part.
`case.resolved.toml` is that gauge twin. For a value within a factor 2 of
p_op the difference is exact (Sterbenz's lemma), so the absolute case and
the gauge case run the same bits.

**The case.** `cht_two_fluids_absolute.toml` is the water / oil
counterflow exchanger of gallery 206 with its two outlets at their
absolute pressures, the water at 1 bar and the oil at 1.5 bar. The log
names the convention and the two operating pressures:

```
[PRS005] pressure values are absolute (D239); the solver works with p - p_op, the outputs are absolute
[PRS005]   part 1 (regions water): operating pressure 100000 Pa (boundary xmax)
[PRS005]   part 2 (regions oil): operating pressure 150000 Pa (boundary oil_out)
```

and `flow.json` holds them as `pressure_values`. The `pressure` field of
`fields.vtkhdf` (attributes `pressure_values`, `operating_pressure`), the
probes and the mean pressures of `boundaries.csv` are absolute; forces
come from p − p_op as before. Every part starts at its own outlet
pressure (gauge 0 in each of its cells), so the jump of half a bar across
the plate never reaches the pressure correction. The residual history and
the mass flows are those of the gauge case to the last bit.

**Checks.** An absolute value of 0 or below is refused (PRS001, the
message suggests `pressure_values = "gauge"`); a value below 1 kPa in a
gas looks like a gauge value (PRS002, warning); a part whose pressure
boundaries are all expressions needs an explicit `operating_pressure`
(PRS003); one case-wide absolute `[initial] pressure` for parts of
different operating pressures is refused (PRS004).

| | |
|---|---|
| Code | `foehncraft_flow::pressure_values` (`resolve`, `Frame`), the changed-default rule of `foehncraft_mesh::case_format`, D239 |
| Case file | `flow/cht_two_fluids_absolute.toml` |
| Tests | `crates/foehncraft-flow/tests/pressure_values.rs` (`an_absolute_case_and_its_gauge_twin_run_the_same_bits`, `two_parts_take_their_own_outlet_pressures`, `older_files_keep_meaning_gauge`, `an_absolute_zero_is_refused`, `values_become_gauge_relative_to_the_outlet`), `case_format::tests::a_changed_default_keeps_the_meaning_of_older_files` |
| Recipe | [220](../recipes/boundary-conditions.md#absolute-pressure) |

```
target/release/foehn flow examples/cases/flow/cht_two_fluids_absolute.toml --out output/two_fluids_absolute
python3 -m tools.plot.pressure_values docs/figures/gallery output/two_fluids_absolute
```

[↑ back to top](#boundary-conditions-of-the-flow-solver)
