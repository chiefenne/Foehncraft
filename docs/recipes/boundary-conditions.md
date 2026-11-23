# Boundary conditions of the flow solver

Part of the [Recipes](../recipes.md): slip walls, wedges for axisymmetric flow, inlet profiles from tables, backflow at pressure outlets, convective outlets, periodic modules, absolute or gauge pressure values (D239). The basic set of conditions is in [recipe 20](flow.md#flow-boundaries); every key is on the [key page](../cases/flow-keys.md#boundary-conditions-d29).

| # | Recipe |
|---|---|
| 29 | [Use slip walls, wedges, table profiles, backflow treatment and convective outlets](#boundary-conditions) |
| 67 | [Compute a fully developed flow on one periodic module](#periodic) |
| 220 | [Give pressures as absolute values, or keep gauge values](#absolute-pressure) |

---

<a id="boundary-conditions"></a>
## 29. Use slip walls, wedges, table profiles, backflow treatment and convective outlets

**Slip walls** stand for a far field or a frictionless surface: no normal
flow and no shear stress, on any surface (`symmetry` is the same on a
plane of symmetry). The uniform flow between two slip walls is exact on
every mesh type:

```
target/release/foehn flow examples/cases/flow/channel_slip.toml
grep "velocity error" output/channel_slip/flow.log
```

```toml
[boundary.ymax]
type = "slip-wall"
```

**Axisymmetric flow on a wedge.** Mesh the half plane y ≥ 0 in 2D (y is
the radius) and revolve it; the wedge planes take `type = "wedge"`:

```
target/release/foehn flow examples/cases/flow/pipe_wedge.toml
# report: output/pipe_wedge/flow.md, plot: convergence.svg
```

```toml
[mesh.generate]
type = "hex"               # any type: hex, tet, poly, mixed
dimension = 2
upper = [2.0, 0.5, 0.1]    # x from 0 to 2, radius from 0 to 0.5
cells = [32, 8, 1]
wedge = 2.0                # degrees about the x axis
[boundary.zmin]
type = "wedge"
[boundary.zmax]
type = "wedge"
```

The axis needs no condition: the cells next to it have an edge there.
With `lower[1] > 0` the mesh is an annulus and `ymin` stays a boundary. A
wedge of θ is a slice of a polygonal pipe of 360°/θ sides; the difference
from the axisymmetric result is about θ²/8, so 2–5° suit applications and
a grid-convergence study takes 0.5–1°.

**Inlet profile from a table.** A measured or computed profile, linear
between the rows:

```
target/release/foehn flow examples/cases/flow/poiseuille_table.toml
head -5 examples/cases/flow/poiseuille_profile.csv
```

```toml
[boundary.inlet]
type = "velocity-inlet"
profile = "table"
table = "profile.csv"           # relative to the case file
table_coordinate = "y"          # or "x", "z", "radius"
# for "radius": the axis
# table_origin = [0.0, 0.0, 0.0]
# table_direction = [1.0, 0.0, 0.0]
```

Two columns are the coordinate and the normal speed into the domain, four
the coordinate and u_x, u_y, u_z. Faces beyond the first or last row take
its value (warning BC020). Changing the table changes the case hash.

**Backflow at a pressure outlet.** When a vortex or a recirculation crosses
the outlet, fluid enters through it. The default `backflow =
"total-pressure"` lets it enter normal to the face with the outlet
pressure as its total pressure; `"suppress"` closes the faces where it
would enter; `"zero-gradient"` is the treatment before stage 5:

```
target/release/foehn flow examples/cases/flow/cylinder_re100_short.toml
target/release/foehn flow examples/cases/flow/cylinder_re100_short.toml --set boundary.outlet.backflow=suppress --out output/cyl_suppress
grep -A6 "C_D,max" output/cylinder_re100_short/flow.md output/cyl_suppress/flow.md
```

**Convective outlet** (transient cases). The velocity leaves by
∂u/∂t + U_c ∂u/∂n = 0, so vortices pass with little reflection and the
outlet can sit a few diameters behind a body:

```
target/release/foehn flow examples/cases/flow/cylinder_re100_convective.toml
cat output/cylinder_re100_convective/flow.md
```

```toml
[boundary.outlet]
type = "convective-outlet"
# convection_velocity = 1.0   # U_c; default the mean outflow speed
```

The outlet flow is balanced with the inflow in every iteration; with no
pressure boundary the reference cell fixes the pressure level
(`[pressure] reference_point`, `reference_value`). In a steady case the
check refuses it (BC021). A restart and a rejected step continue bitwise:
the old boundary levels are in the restart file.

[↑ back to top](#boundary-conditions-of-the-flow-solver)

---

<a id="periodic"></a>
## 67. Compute a fully developed flow on one periodic module

Make the inlet and outlet of the module a periodic pair and give the
driving on one of them, the mass flow (kg/s) or the mean pressure
gradient (Pa/m, the pressure falling towards the partner):

```toml
[boundary.inlet]
type = "periodic"
partner = "outlet"
mass_flow = 0.1            # or: pressure_gradient = 1.2
# translation = [2.0, 0.0, 0.0]   # default: from the face centroids
[boundary.outlet]
type = "periodic"
```

The two boundaries must carry translated copies of the same faces (a
generated box keeps them so when it is perturbed); the mesh needs at least
two cells along the period. Several pairs (a channel periodic in x and z)
are fine; the pressure level is held in the reference cell. Check the
case, run it and read the driving force:

```
target/release/foehn flow examples/cases/flow/poiseuille_periodic.toml --check
target/release/foehn flow examples/cases/flow/poiseuille_periodic.toml
grep "periodic xmin" output/poiseuille_periodic/flow.log
```

The report (`result.json`, `periodic`) gives per pair the mass flow and
the driving force β (N/m³, the mean pressure gradient it stands for with
the opposite sign); the friction factor is f = β D_h/(½ ρ U_b²).

**A heated module.** A scalar heated through the walls by a given flux
crosses the pair with the jump Δ = Q/ṁ when one side of the pair says so;
`bulk` sets its mixing-cup value there (default 0). Fixed wall or inlet
values are refused with this mode, because the developed field is then
plainly periodic or decays (the wall-temperature mode comes with the
energy equation):

```toml
[boundary.inlet.scalar.theta]
periodic = "heat-flux"
bulk = 0.0
[boundary.wall.scalar.theta]
flux = 1e-3
```

```
target/release/foehn flow examples/cases/flow/channel_heated_periodic.toml
target/release/foehn flow examples/cases/flow/pipe_heated_periodic.toml
```

The Nusselt number follows from the report: per side of the pair
`mixed_mean` is the bulk value (side `b` is Δ above side `a`), per wall
`mean` the area mean; with the hydraulic diameter D_h,
Nu = q D_h / (Γ (T̄_w − T_b)), with the bulk value in the middle of the
module, (T_b,a + T_b,b)/2, against the mean wall value. Without the
heat-flux mode (`periodic = "periodic"`, the default) the scalar crosses
the pair unchanged.

[↑ back to top](#boundary-conditions-of-the-flow-solver)

---

<a id="absolute-pressure"></a>
## 220. Give pressures as absolute values, or keep gauge values

**When.** Always, from case format 0.7: the pressure values of a case file
are absolute static pressures unless the file says otherwise (D239).

```toml
case_format = "0.7"

[boundary.outlet]
type = "static-pressure-outlet"
pressure = 1.01325e5          # Pa, absolute: an outlet to the atmosphere

[boundary.inlet]
type = "total-pressure-inlet"
total_pressure = 1.02325e5    # 1 kPa above the outlet
```

The operating pressure p_op of each fluid part comes from its own outlet
(here 101 325 Pa), and the solver works with p − p_op. A case of several
fluid parts (a heat exchanger) needs nothing more: each part takes its own
outlet pressure, and its gas density follows it. To fix p_op yourself (an
outlet given only by an expression, PRS003):

```toml
[fluid]
operating_pressure = 1.5e5    # every part
[regions.water]
operating_pressure = 1.0e5    # one region's part
```

**Gauge values.** Values relative to the operating pressure, as in files
of format 0.6 and older:

```toml
[fluid]
pressure_values = "gauge"
operating_pressure = 1.0e5    # default 101 325 Pa
[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0
```

**Older files.** They read as gauge values with a note in the log, and keep
their results. `foehn case upgrade` writes the meaning into them:

```
target/release/foehn case upgrade examples/cases/flow/my_case.toml --write
```

**What to check.** The log lines `[PRS005]` and `flow.json` →
`pressure_values` state the convention and p_op of every part with its
source; `case.resolved.toml` shows the gauge values the solver used. The
`pressure` field, probes and `boundaries.csv` are in the file's
convention; pressure drops and forces are the same in both. An absolute
value of 0 Pa is refused (PRS001): it was meant as a gauge value.

**The two-fluid exchanger** (gallery 220):

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cht_two_fluids_absolute.toml --out output/two_fluids_absolute
grep PRS005 output/two_fluids_absolute/flow.log
python3 -c "import json; print(json.load(open('output/two_fluids_absolute/flow.json'))['pressure_values'])"
```

[↑ back to top](#boundary-conditions-of-the-flow-solver)
