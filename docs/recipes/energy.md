# Energy equation and heat transfer

Part of the [Recipes](../recipes.md): the energy equation of the fluid (stage 7, D150), its properties, thermal wall conditions and heat balance. The keys are on the [key page](../keys.md).

| # | Recipe |
|---|---|
| 154 | [Solve the temperature of a flow with temperature-dependent properties](#energy-fluid) |
| 156 | [Couple the temperature of a fluid and a solid (conjugate heat transfer)](#cht-solids) |
| 158 | [Couple a turbulent flow with a solid (turbulent conjugate heat transfer)](#cht-turbulent) |
| 159 | [Follow how a solid and a fluid heat up together (transient conjugate heat transfer)](#cht-transient) |
| 160 | [Put a contact resistance between two parts, orient an orthotropic solid, heat a component](#cht-contact) |
| 180 | [Bring the solids of a transient conjugate run to their state sooner (thermal acceleration)](#cht-acceleration) |
| 186 | [Switch on gravity and buoyancy; which pressure do my boundary values mean](#buoyancy) |
| 195 | [Let the density of a gas follow its temperature (low-Mach variable density)](#variable-density) |

---

<a id="energy-fluid"></a>
## 154. Solve the temperature of a flow with temperature-dependent properties

**The equation.** An `[energy]` table adds the temperature T of the fluid,
solved after the pressure correction of every outer iteration:

ρ c_p ∂T/∂t + ∇·(ṁ c_p T) − T ∇·(ṁ c_p) − ∇·(k ∇T) = q''' + Φ.

The matrix is in T, and a correction on the right-hand side makes the
converged result conservative in the sensible enthalpy h_s = ∫c_p dT
(zero at 298.15 K), so a c_p that changes with T does not lose heat. The
equation uses the convection schemes, limiters, turbulent Prandtl number
and thermal wall functions of the passive scalars (recipe 63).

```toml
[fluid]
material = "air"               # c_p(T), k(T) and μ(T) from the database
# or the values themselves (numbers or expressions of T):
# density = 1.0
# viscosity = 0.01
# specific_heat = 1000.0       # J/(kg K)
# conductivity = "0.02 + 7e-5 * (T - 300)"   # W/(m K)

[energy]
initial = 300.0                # K, a number or an expression of x, y, z
source = 0.0                   # q''' (W/m³), expression of x, y, z (and t)
viscous_heating = false        # Φ = μ (∇u + ∇uᵀ):∇u
temperature_dependent = true   # refresh c_p, k, μ from T every outer iteration
exact = "..."                  # optional: L2 error and observed order of a sweep

[solver.energy]
relaxation = 1.0               # default 0.9 steady
tolerance = 1e-8
max_iter = 200
```

`[fluid] specific_heat` and `conductivity` replace those of the material.
With `temperature_dependent = true` (the default) μ(T) of the material
also reaches the momentum equation; the density stays constant until the
variable-density step of stage 7 (7c).

**Thermal boundary conditions.** Walls are adiabatic unless a
`[boundary.<name>.energy]` table gives exactly one condition:

```toml
[boundary.hot.energy]
temperature = 320.0            # K, a number or an expression

[boundary.heated.energy]
heat_flux = 500.0              # W/m² into the fluid

[boundary.heater.energy]
heat_flow = 0.5                # W over the whole boundary

[boundary.outer.energy]
heat_transfer_coefficient = 10.0   # W/(m² K)
ambient_temperature = 293.15       # K: q = h (T∞ − T_w)

[boundary.inlet.energy]
temperature = 300.0            # inlets need one (THM011)
```

The condition follows from the keys, so there is no `type` key. Walls and
slip walls take any of them; walls without a table are adiabatic and the
check lists them once (THM005). Outlets keep a zero gradient, symmetry
planes and wedges pass no heat, and periodic pairs carry T across: an
energy table there is refused (THM003), as a flux at an inlet is.

The check before the run (`foehn flow --check`) also refuses:
- two conditions or none in one table, or a table of an unknown boundary
  (THM010), and an inlet without a temperature (THM011);
- a temperature, ambient temperature or start temperature ≤ 0 K, and a
  heat-transfer coefficient ≤ 0 (THM004);
- a steady case with no inflow and no temperature or convection boundary.
  That is THM001 when the temperature level is left open, or THM002 when a
  heat input has no way out;
- a boundary table without `[energy]` (THM009);
- a passive scalar named `temperature`, a fluid region with a material of
  its own, or an orthotropic fluid conductivity (ENE006), and a fluid
  without c_p or k (MAT007). Solid regions: recipe 156.

**The heat balance.** The run report (`flow.md`, `flow.json`,
`result.json`) and the log give, for the final state, the conduction
through the boundaries, the enthalpy carried in and out by the flow, the
sources (q''' and Φ) and the imbalance, and per boundary the mean
temperature, the heat flux (W/m²) and the heat flow (W). A Nusselt number
follows as Nu = q L/(k ΔT) from the wall flux; in the heated cavity of
gallery 154 it is 2.05.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cavity_heated.toml --out output/cavity_heated
target/release/foehn flow examples/cases/flow/couette_viscous_heating.toml --out output/couette
target/release/foehn flow examples/cases/flow/conduction_mms.toml --out output/conduction
```

**Limits of this step.** The property refresh and the enthalpy correction
run on the host: on a GPU every outer iteration copies the temperature
down and the properties up. With a turbulence model the wall functions and
the turbulence sources still use the reference viscosity. The
streamwise-periodic heat-flux mode of the passive scalars is not yet
available for the energy equation.

---

<a id="cht-solids"></a>
## 156. Couple the temperature of a fluid and a solid (conjugate heat transfer)

**One mesh, two regions.** The fluid and the solid share one mesh; a cell
zone marks the solid cells. With a mesh file the zones come from the file;
the box generator takes boxes, and a cell whose centroid lies in a box
goes to its zone (the last box that holds it wins, all other cells stay in
`fluid`):

```toml
[mesh.generate]
type = "hex"
dimension = 2
upper = [4.0, 1.2, 0.04]
cells = [80, 30, 1]
zones = [{ name = "wall", lower = [0.0, 0.0, 0.0], upper = [4.0, 0.2, 0.04] }]

[materials.wall-solid]
kind = "solid"
density = 2000.0               # kg/m³
specific_heat = 800.0          # J/(kg K)
conductivity = 140.8           # W/(m K), a number or an expression of T

[regions.wall]
kind = "solid"
zones = ["wall"]
material = "wall-solid"

[energy]
initial = 300.0
```

The flow runs on the fluid cells only. Every face between a fluid and a
solid cell becomes a no-slip wall of the flow named `interface-<region>`
(here `interface-wall`); it needs no `[boundary]` table and takes no
thermal condition, because the temperature is solved on all cells as one
matrix. At a face between two materials the conductance is the series
value of the two half cells, 1/((1 − w)/k_P + w/k_N), so a piecewise
linear temperature is exact and the heat flux is the same on both sides
to round-off. Outer boundaries of the solid are walls and take the
thermal conditions of recipe 154 (`temperature`, `heat_flux`,
`heat_flow`, `heat_transfer_coefficient`).

**The report.** The heat balance covers the whole domain. The boundary
table lists every interface with its mean temperature (the face value of
the series), heat flux and heat flow into the fluid, so the heat a wall
passes to the flow is read off directly; in the heated wall of gallery
156 it is all of the heater's 160 W. `fields.vtkhdf` and the lines are
written on the whole mesh (velocity and pressure zero in the
solids).

**The checks.** Besides those of recipe 154, `foehn flow --check`
refuses:
- a solid region without a material, or a solid material without
  density, specific heat or conductivity (MAT006);
- a zone the mesh does not have (REG002);
- solid regions without `[energy]` (REG001);
- solids together with a sliding mesh, a restart or passive scalars, an
  orthotropic solid conductivity, a periodic mesh, or a mesh with no fluid
  cell (ENE006). A turbulence model is allowed since PR 158 (recipe 158),
  a time loop since PR 159 (recipe 159).

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cht_composite_wall.toml --out output/cht_wall
target/release/foehn flow examples/cases/flow/cht_heated_wall.toml --out output/cht_heated
target/release/foehn flow examples/cases/flow/cht_heated_wall.toml --out output/cht_heated_100 \
    --set 'materials.wall-solid.conductivity=1408.45'
```

**Limits of this step.** Steady conjugate cases only (turbulent ones
since PR 158, recipe 158; transient ones since PR 159, recipe 159). The
gradients next to an interface still use the cell across it, so the
interface is exact on orthogonal meshes; the one-sided gradients of D157
follow. No contact resistance and no orthotropic k yet (7d2), no
non-matching interfaces (7d3), and the flux, gradient and viscosity
fields are mapped between the fluid and the whole mesh on the host in
every outer iteration. The material database has no solids yet, so cases
define their own.

<a id="cht-turbulent"></a>
## 158. Couple a turbulent flow with a solid (turbulent conjugate heat transfer)

**The case.** A conjugate case of recipe 156 with a `[turbulence]`
model. Nothing else changes: the solid region, its material and the
thermal conditions of the outer walls stay as they are, and the
interfaces `interface-<region>` still need no `[boundary]` table. The
thermal wall law and Pr_t come from `[energy]`, as for every wall of a
turbulent flow (recipe 154):

```toml
[turbulence]
model = "k-omega-sst"          # or k-epsilon, ..., with their wall treatments
wall_treatment = "automatic"

[energy]
initial = 300.0
turbulent_prandtl = 0.85       # c_p μ_t/Pr_t in the fluid
wall_function = "kader"        # or "jayatilleke", or "none" (laminar k/y)
```

**What happens at an interface.** The fluid side of a face between a
fluid and a solid cell is the thermal wall function of the turbulent
flow, h_w = Γ_w/y_P with Γ_w = μ c_p y*/T⁺(ỹ*) (D135; ỹ* = max(y*, 11.06)
with the scalable treatment), in series with the half cell of the solid:
q = (T_P − T_N)/(1/h_w + δ_N/k_s). y* and μ of the first fluid cell come
from the flow in every outer iteration, so the conductance follows the
boundary layer as it develops. With a resolved wall (y⁺ ≈ 1, `low-re`
or the automatic treatment near the wall) the law gives back the laminar
k/y_P. Inside the fluid the conductivity is k + c_p μ_t/Pr_t; inside the
solid it stays k_s.

**The report.** As in recipe 156: the heat flow and mean temperature of
every interface, the heat balance of the whole domain, and the first-cell
y⁺ of `interface-<region>` in the log (the wall-function range applies
there as at any wall of the flow).

**Starting the flow.** A case where the flow meets a solid head on (fins,
ribs, a block in the channel) starts from rest, `[initial] velocity =
[0, 0, 0]`: a uniform inlet velocity pushed into the solid's front face
upsets the first SIMPLEC iterations (the pressure-correction matrix stops
being positive definite). A lower momentum relaxation for the first
iterations (`[[schedule]]`, [recipe 39](runs.md#schedule)) helps further.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cht_turbulent_channel.toml --out output/cht_channel
target/release/foehn flow examples/cases/flow/cht_heat_sink_turbulent.toml --out output/cht_sink
```

**Limits of this step.** Steady only (transient conjugate heat transfer
since PR 159, recipe 159). The interface takes the wall function of the momentum
treatment of the flow; a wall-resolved mesh on the fluid side needs
y⁺ ≲ 1 like any wall. The turbulence fields are not written on the whole
mesh yet (velocity, pressure and temperature are).

<a id="cht-transient"></a>
## 159. Follow how a solid and a fluid heat up together (transient conjugate heat transfer)

**The case.** A conjugate case of recipe 156 or 158 with a `[time]`
table. The solid regions need nothing new: their density and specific
heat, which a steady run does not use, now set how much heat they store.

```toml
[time]
scheme = "bdf2"                # or "crank-nicolson", "bdf1"
time_step = 5.0
end_time = 3000.0

[solver]
outer_iterations = 3           # the flow at Courant numbers far above 1

[solver.momentum]
relaxation = 0.9
```

The time step follows the heat-up of the solid, not the flow: a wall of
heat capacity ρ c_p t per area that loses heat with h to the flow warms
with the time constant ρ c_p t/h, often minutes, while the flow passes in
seconds. At Courant numbers of 10 to 100 the flow needs two or three outer
iterations per step and a momentum relaxation below 1 (with one outer
iteration the relaxation must be 1, CASE005). A flow that does not
change, as here, settles in the first steps; the energy equation of all
regions is solved in every step to `[solver.energy] tolerance`. The
temperature takes BDF2 whatever the scheme of the flow, as every scalar
([recipe 63](scalars.md#passive-scalars)), so a Crank–Nicolson run gives the same second
order in time.

**What happens.** The fluid and the solid cells share one matrix with the
time terms of both, ρ c_p V/Δt of each cell's own material, so the heat
flow through an interface and the heat stored on its two sides come from
the same solve; there is no coupling iteration between a fluid and a
solid solver and no time lag at the interface. The time levels of the
temperature follow those of the flow: an adaptive step, a rejected step
(D65) and the return to the last valid state (D63) rewind the solid with
the fluid.

**The output.**
- `heat.csv`: the heat balance of every step, what the boundaries bring
  in by conduction and with the flow, the sources, the heat stored in the
  step per region (`fluid_storage`, `<region>_storage`), the heat flow
  through every interface (`interface-<region>_heat`, into the fluid),
  the imbalance and its ratio to the largest flow. The run record
  (`run.jsonl`) carries the same per step under `heat`.
- `probes.csv`: `<probe>_T`, the temperature of the probe's cell, also in
  a solid (velocity and pressure are zero there, as in the fields).
- `fields_<step>.vtkhdf` (`fields_every`) on the whole mesh.
- the report: the largest relative imbalance of a step
  (`transient.heat_relative_max`) and the heat stored per region in the
  last step.

**An exact solution.** `[energy] exact` takes `erf` and `erfc`, so a
conduction problem with an exact solution reports the L2 error of the
temperature and, with `[verification] levels` or `time_steps`, the
observed orders (`scalar_order`); `cht_contact_transient.toml` (two bodies
brought into contact) does this.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cht_heated_wall_transient.toml --out output/cht_wall_t
target/release/foehn flow examples/cases/flow/cht_contact_transient.toml --out output/cht_contact
target/release/foehn flow examples/cases/flow/cht_contact_transient.toml --out output/cht_contact_time \
    --set 'verification.time_steps=[20.0, 10.0, 5.0, 2.5]'
```

**Limits of this step.** No sliding meshes with solids, no restart of a
conjugate run yet, and the turbulence fields are not written on the whole
mesh. A periodic flow whose solid only follows the cycle-average (an
exhaust manifold over an engine cycle) is the cycle-averaged coupling of
D204, which comes later; until then such a case runs the full transient,
with the thermal acceleration of recipe 180 to bring the solids to their
cycle-mean state sooner.

<a id="cht-contact"></a>
## 160. Put a contact resistance between two parts, orient an orthotropic solid, heat a component

**Contact resistance and thin walls.** Every face between two regions of
a conjugate case is found from the cell zones; an `[interfaces.<name>]`
table only adds a resistance to the faces between two of them:

```toml
[interfaces.pad]
regions = ["chip", "board"]         # solid regions, a fluid region, or "fluid"
contact_resistance = 2.0e-4         # m²K/W

[interfaces.coating]
regions = ["fluid", "housing"]
thin_wall = { thickness = 5.0e-4, conductivity = 0.2 }   # t/k; or material = "<name>"
```

The face conductance is then the series
|S|/(R_P + R″_c + R_N) of the two half cells and the resistance (D157
§9.1); the heat flow is the same on both sides of the face, and the two
face temperatures differ by q R″_c. A `thin_wall` here is a resistance
without heat capacity; a `material` gives k(T) at the face temperature.
Two fluids cannot have an interface (a wall between two fluids is the thin
wall of D189, later), and an entry whose regions do not touch is refused
(CHT004). Typical values: a thermal pad or grease 10⁻⁵ to 10⁻⁴ m²K/W, two
machined metal faces pressed together 10⁻⁴ to 10⁻³, an oxide or paint
layer t/k.

**Orthotropic solids.** A solid material may have three principal
conductivities along two axes (the third is normal to both):

```toml
[materials.board-laminate]
kind = "solid"
density = 1900.0
specific_heat = 1100.0
conductivity = { form = "orthotropic", values = [25.0, 0.35, 25.0], axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]] }

[regions.board_2]                   # the same material, turned
kind = "solid"
zones = ["board_2"]
material = "board-laminate"
axes = [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]]
```

`axes` of a region replace those of its material, so one winding or
laminate material can be oriented differently in every part (MAT008
without axes). Each value can be a form of T. The matrix takes the
conductivity of the face normal, n̂·K·n̂; the rest of the flux, the part
that the tangential gradient drives through an inclined K, is a deferred
correction, so that a linear temperature in an orthotropic solid is exact
on any mesh (`cht_orthotropic_patch.toml`). Frames that follow the
geometry (cylindrical frames of e-motor windings, frames from fields) come
later (L48).

**Heat sources of a region.** `[[energy.sources]]` heats a single
region, in addition to `[energy] source`:

```toml
[[energy.sources]]
region = "chip"
power = 1.0                         # W over the region, spread by volume (a number or an expression of t)

[[energy.sources]]
region = "winding"
power_density = "2.0e5 * (1 + 0.004 * (t - 10))"   # W/m³, of x, y, z, t
```

The power enters the heat balance as a source (`sources` in the report
and in `heat.csv`); ENE008 for an unknown region or without exactly one
of the two keys. In a 2D case `power` is the power of the slice, so give
the depth of the mesh the value it has.

**Interface faces and their gradients.** Next to an interface the
gradient of T jumps; a cell there takes the face temperature of its own
side in place of the cell across in its least-squares gradient (one-sided
gradients, §9.2), and the flux through the face uses the face-normal
temperatures of both sides, so that skewed cells at an interface give no
error for a piecewise-linear temperature. This is automatic in every
conjugate case and needs no key. The report gives each pair of regions:
`interface-<solid>` between the fluid and a solid (heat into the fluid,
temperature of the fluid side) and `interface-<a>-<b>` between two solids
(heat into `a`, the first in alphabetical order, temperature of its side).

**Converge the energy solve.** A conjugate matrix with large conductivity
ratios needs more BiCGStab iterations than a fluid alone; if the run
converges with a heat balance that does not close (ENE007), raise
`[solver.energy] max_iter` (300 to 500) and tighten its `tolerance`
(10⁻¹⁰): the circuit board of entry 160 closes to 10⁻¹² with 500.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cht_circuit_board.toml --out output/cht_board
target/release/foehn flow examples/cases/flow/cht_orthotropic_patch.toml --out output/cht_patch
```

**Limits of this step.** The cross terms of an orthotropic solid at a
convection boundary are taken with the weight of the face value (exact
at a fixed temperature, a flux boundary carries its given flux); thin walls
with heat capacity, layers and in-plane conduction (D189), fluid–fluid
walls, and the entropy production of a contact resistance in the
second-law check come later. Periodic meshes with solids (the thick-walled
channel of V7.10) are not yet possible (ENE006).

<a id="cht-acceleration"></a>
## 180. Bring the solids of a transient conjugate run to their state sooner (thermal acceleration)

A solid warms up with its time constant, often minutes, while the flow
settles in seconds or follows a cycle of milliseconds; a transient run
that only needs the operating point (a steady state, or the cycle mean
of a periodic flow) spends most of its steps waiting for the solid. The
steady temperatures of a solid do not depend on its ρ c_p, so the run can
scale ρ c_p of the solids down until they have settled, and then give it
back (D213):

```toml
[energy.acceleration]
inertia_scale = "auto"     # or a factor 0 < f <= 1
# regions = ["wall"]       # default: every solid region
# restore = "settled"      # or a time in s at which the restore starts
# period = 0.04            # s, of a periodic flow: a window is one period
# window_steps = 10        # without a period: steps per window
# settle_tolerance = 1e-3  # change still to come, of the solids' rise
# restore_ramp = 3         # windows over which f rises back to 1 (3 with a period, 0 without)
# target_windows = 2       # auto: scaled time constant (2 with a period, 1 without)
# settle_windows = 2       # windows after the restore before the history counts as physical
```

**What happens.** With a factor, the solids run with f ρ c_p from the
first step. With `"auto"`, the first window runs with the true ρ c_p;
at its end the run measures each solid region's conductance G to the rest
(the change of its interface heat flow when it is 1 K warmer, from two
steady solves of the energy equation with every solid held and the flow as
it is, plus its outer boundaries with a temperature or a convection
condition), its true heat capacity C and its time constant τ = C/G, and
sets f = (target windows × window)/τ, at least 10⁻⁴. τ leaves out the
conduction inside the solid (a small Biot number); where it matters, the
acceleration is milder than asked for, never stronger. Every window the
run checks whether the solids' window-mean temperatures have settled; then
the true ρ c_p comes back: at once without a period (exact, the flow is
not periodic), raised geometrically over `restore_ramp` periods with a
period (a jump would leave the part of the accelerated swing that the
solid had at that moment in its deep temperature, where it decays only
with the true time constant; gallery 180).

**What is exact and what is not.** The steady state and the cycle mean of
a periodic flow are those of the true solid (to the settle tolerance;
the cycle mean to second order in the swing through h(T_w), μ(T) and
radiation). The history before the restore is not physical: the solids'
time runs 1/f times faster, and in a periodic flow their skin depth and
surface swing are 1/√f too large while accelerated (ACC003 when the skin
reaches through the region). The result gives `physical_from`, the restore
time plus `settle_windows` windows; use only what comes after it for
cycle averages, swings and histories. A run that ends before the restore
warns (ACC002).

**The output.**
- `heat.csv`: the column `inertia_scaling`, the heat the true ρ c_p would
  have stored in addition (W): storage + inertia_scaling is the storage of
  the true solids; the balance of the solved equation (`imbalance`) closes
  as without acceleration.
- the report, `transient.acceleration`: per region the factor, C, the
  measured G and τ (auto), the skin depth and thickness (with a period) and
  the heat booked as inertia scaling over the run (J); the times of the
  measurement, of the settling, of the start and end of the restore,
  `physical_from` and `history_physical`.
- the log and the run record: one ACC004 line for every step of it.

```
target/release/foehn flow examples/cases/flow/cht_heated_wall_accelerated.toml --out output/acc
target/release/foehn flow examples/cases/flow/cht_heated_wall_periodic_accelerated.toml --out output/per
```

**When not to use it.** A steady conjugate run has no ρ c_p of a solid in
its equations: nothing to accelerate (ACC001). When the heat-up history
itself is wanted (a component warming in a steady cooling flow, a drive
cycle), the solids must keep their true ρ c_p; the frozen-flow march of
D118, which keeps the true history with large steps, comes later. For an
exhaust manifold over many engine cycles the cycle-averaged loop of D204
(also later) solves the solid steady with the cycle-averaged heat transfer:
the limit f → 0 of this acceleration, without the swing.

<a id="buoyancy"></a>
## 186. Switch on gravity and buoyancy; which pressure do my boundary values mean

**Switch on gravity and the buoyancy.** Two tables, next to `[energy]`:

```toml
[gravity]
vector = [0.0, -9.81, 0.0]          # m/s²
reference_point = [0.0, 0.0, 0.0]   # x_ref of the reduced pressure (default)

[buoyancy]
model = "boussinesq"
expansion_coefficient = 3.333e-3    # β (1/K); default from the fluid's density
# reference_temperature = 300.0     # T_ref (K); default [fluid] temperature
```

The Boussinesq model keeps the density ρ₀ constant everywhere except in
the body force, which becomes −ρ₀ β (T − T_ref) g. It holds while β ΔT is
small; above 0.1, BUO002 warns, and the variable density of [recipe
195](#variable-density) is the model to use. Without `expansion_coefficient`, β comes from the
fluid's density at T_ref (1/T_ref for an ideal gas); a constant density
without β is refused (BUO003), and so is a buoyancy model without
`[gravity]` or `[energy]` (BUO001, BUO003). `[gravity]` with the energy
equation but no model warns (BUO004): the temperature then drives no
flow. The force enters the faces exactly like the pressure gradient (the
balanced force of D154), so a fluid at rest in a stratification stays at
rest on any mesh ([gallery 186](../gallery/energy.md#buoyancy)).

**Steady runs.** SIMPLEC and the coupled algorithm take a pseudo-time
step of 1/N, with N = √(|g| β ΔT / L) the buoyancy frequency, so that the
internal waves of a stratification settle instead of being overshot.
`[solver.pseudo_time] factor` scales the step; a factor above 1 rarely
helps, because the slow part is the plume or the roll finding its place.
A laminar plume takes about 2000 iterations on 6000 cells.

**Which pressure do my boundary values mean.** With `[gravity]` the
solver works with the **reduced pressure**

p′ = p − ρ_ref g·(x − x_ref),

the static pressure minus the hydrostatic pressure of a fluid at ρ_ref.
p′ is what drives the flow: a fluid at rest at ρ_ref has p′ = 0
everywhere, at every height. The result has both, `pressure` (p′) and
`static_pressure` (p).

The values of the pressure boundaries (`static-pressure-inlet`,
`total-pressure-inlet`, `static-pressure-outlet`,
`total-pressure-outlet`) are **p′ by default**:

| `pressure_kind` | The value means | 0 means | Use it for |
|---|---|---|---|
| `"reduced"` (default) | p′ = p − ρ_ref g·(x − x_ref) | an opening to an ambient at rest at ρ_ref, at every height | openings, windows, tall outlets, anything open to a still surrounding |
| `"static"` | the static pressure p; the solver uses p′_b = p_b − ρ_ref g·(x_f − x_ref) per face | 0 Pa static at every face: a uniform static value | a measured or specified static pressure at one height, a gauge reading |

<img src="../figures/gallery/buoyancy_pressure_kinds.png" alt="A liquid rising at 0.05 m/s in a channel 0.1 m wide and 2 m tall with gravity downwards. Left: the reduced pressure p prime, which falls by only 6.5 Pa from bottom to top, the friction. Middle: the static pressure p, which falls by 19.6 kPa, mostly the weight of the liquid. Right: p prime and p along the channel axis for an outlet given as p prime = 0 and for one given as the static pressure minus 19 620 Pa; the curves of the two runs lie on top of each other" width="800" height="431">

A liquid (ρ = 1000 kg/m³) rising at 0.05 m/s in a channel 2 m tall shows
both. `buoyancy_channel_reduced.toml` gives the outlet as p′ = 0, and
`buoyancy_channel_static.toml` gives it as the static pressure there,
p = −ρ g y = −19 620 Pa (relative to the inlet), with
`pressure_kind = "static"`. Both runs give the same flow to 10⁻¹¹. p′
falls by 6.5 Pa, the friction of the channel; p falls by 19.6 kPa, almost
all of it the weight of the liquid.

**The trap: a uniform static value over a tall opening.**
`buoyancy_open_box.toml` is a box of fluid at rest, open over its full
height of 1 m on one side. With the default kind, `pressure = 0.0` is the
ambient's own hydrostatic profile, and the fluid stays at rest (largest
speed 0). The same 0 Pa as a static pressure,

```
target/release/foehn flow examples/cases/flow/buoyancy_open_box.toml --out output/open_static \
    --set 'boundary.xmax.pressure_kind="static"'
```

is a uniform static pressure from bottom to top, while the fluid inside
carries 9.3 Pa of hydrostatic difference across the opening. Fluid is
pushed out at the bottom and drawn in at the top, at up to 3.5 m/s. The
physics has no such flow, and BUO006 warns when the hydrostatic difference
over a static-kind opening exceeds its dynamic pressure. Use the default
kind, or give the static profile as an expression of the height.

**Interpreting the results.** `static_pressure` is p = p′ + ρ_ref g·(x −
x_ref), the pressure a gauge would read relative to the reference point.
The wall forces of the report are integrated from the reduced pressure
p′. They leave out the uniform hydrostatic part, so a submerged body has
no Archimedes force ρ_ref g V in them; add it by hand where it matters.
The kinetic-energy ledger
(`balance.csv`) counts the work of the buoyancy force, so it still
closes.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/buoyancy_channel_reduced.toml --out output/chan_reduced
target/release/foehn flow examples/cases/flow/buoyancy_channel_static.toml --out output/chan_static
target/release/foehn flow examples/cases/flow/buoyancy_open_box.toml --out output/open_reduced
```

---

<a id="variable-density"></a>
## 195. Let the density of a gas follow its temperature (low-Mach variable density)

**When.** A gas whose temperature changes by more than about a tenth of
its absolute temperature: hot air in an oven or an enclosure with a
heater, an exhaust gas cooled in a duct, air that a heater expands in a
closed volume, a liquid heated by tens of kelvins on its way through a
pipe. The flow stays slow (Mach below
about 0.3), so the pressure waves do not matter and the density depends
on the temperature and on one thermodynamic pressure p₀, not on the local
pressure (the low-Mach model, D151). For a few kelvins the Boussinesq
model of [recipe 186](#buoyancy) is enough and cheaper.

**Switch it on.** One key in `[energy]`:

```toml
[fluid]
material = "air"                    # or a material with a density form
temperature = 600.0                 # the start state and T_ref

[energy]
initial = 600.0
density = "variable"                # ρ(T, p₀) in every cell and face
# thermodynamic_pressure = "constant"   # p₀ = [fluid] operating_pressure (default)
```

The density comes from the fluid's density form: an ideal gas ρ = p₀ M/(R
T) (`density = { form = "ideal-gas" }` and a `molar_mass`, as air of the
material library has), or a table or law of the temperature for a liquid.
The momentum equation takes the full stress with −⅔ μ ∇·u, the pressure
equation the storage V ∂ρ/∂t of a transient run, and the energy
equation ρ per cell. Cases without the key keep the constant density of
`[fluid] temperature` and are bitwise unchanged. With gravity, `[buoyancy]
model = "variable-density"` switches the density on as well and gives the
body force (ρ − ρ_ref) g; the Boussinesq model and a constant density
form are refused (ENE009).

**Open or closed.** In a domain open to an ambient (an inlet and an
outlet, a pressure boundary) p₀ is `[fluid] operating_pressure` and
constant. In a closed domain (walls, symmetry planes and nothing else)
the gas cannot leave, so p₀ follows from its fixed mass:

```toml
[energy]
density = "variable"
thermodynamic_pressure = "closed-volume"   # p₀(t) from the mass of the start state
# pressure_work = "thermodynamic"          # dp₀/dt as a heat source (the default here)
```

A closed domain with a constant p₀ is refused (BUO005): heating would
create or destroy mass. A closed volume needs an ideal gas, and an open
domain cannot be one (ENE009). The pressure work dp₀/dt is on by default
in a closed volume, so that a heated gas warms with c_v; `pressure_work =
"none"` leaves it out and the gas warms with c_p, which is the right
answer only for an open domain. `vd_closed_box_heating.toml` is the
check: a box at rest heated with 10⁴ W/m³ for 5 s reaches T and p₀
exactly, at any time step ([gallery 195](../gallery/energy.md#variable-density)).

**Steady runs.** SIMPLEC and the coupled algorithm both work and give
the same discrete solution (checked on the manufactured solution). A
closed cavity starts from rest; give a flow through a duct its inflow
velocity as `[initial] velocity`: from rest a heat source heats the still
gas far beyond its range before the flow arrives, and the run diverges. The large-ΔT cavity of gallery 195 takes
266, 808 and 2719 coupled iterations on 32², 64² and 128² cells clustered
to the walls, which the cold layer needs: on uniform meshes up to 96² the
heat flow is not yet in the asymptotic range.

**Interpreting the results.** The output has a `density` field. The
report has p₀, whether the volume is closed, the mass at the start and
at the end, the density range and the largest Mach number (ENE004 warns
above 0.3, where the low-Mach model no longer holds); `history.csv` of a
transient run has the columns `p0`, `mass`, `rho_min` and `rho_max`. A
mesh sweep gives the observed order and the GCI of p₀ and of the wall
heat flows. The ledgers book the new terms: the mass account the storage
of the density (it closes per step, `vd_closed_box_side_heated.toml`), the
heat account the storage ρ ∂h/∂t + h ∂ρ/∂t, and the kinetic-energy account
the pressure dilatation p ∇·u, the work of the expanding gas, as a
physical term, not numerical dissipation.

**Not yet.** Passive scalars, a restart, periodic meshes, convective
outlets, mass-flow inlets, wall mass fluxes and rotating zones are
refused with ENE006. The outlet value ⅔ρk of a turbulence model uses the
reference density.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/vd_cavity_large_dt.toml --out output/cavity
target/release/foehn flow examples/cases/flow/vd_heated_duct.toml --out output/duct
target/release/foehn flow examples/cases/flow/vd_closed_box_heating.toml --out output/box
target/release/foehn flow examples/cases/flow/vd_closed_box_side_heated.toml --out output/side
target/release/foehn flow examples/cases/flow/vd_mms_steady.toml --out output/vd_mms
```
