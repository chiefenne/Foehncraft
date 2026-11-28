# Two fluids (volume of fluid)

Part of the [Recipes](../recipes.md): two immiscible fluids with a sharp interface (L3, volume of fluid). The keys are on the [key page](../cases/flow-keys.md) (`vof.*`).

| # | Recipe |
|---|---|
| 317 | [Put water and air into one case](#vof-two-fluids) |
| 358 | [Let a liquid evaporate or boil](#vof-phase-change) |
| 368 | [Put dispersed bubbles into a liquid (two-fluid model)](#two-fluid) |

---

<a id="vof-two-fluids"></a>
## 317. Put water and air into one case

**1. Give the two phases and where the first one is at the start.** Phase 1 is α = 1.

```toml
[vof]
scheme = "thinc"          # or "cicsam", "upwind"
surface_tension = 0.07    # N/m (default 0)

[[vof.phase]]             # α = 1
name = "water"
density = 998.2
viscosity = 1.0e-3

[[vof.phase]]             # α = 0
name = "air"
density = 1.205
viscosity = 1.82e-5

[[vof.region]]            # in order; later regions overwrite earlier ones
shape = "box"             # box, sphere, cylinder or field
lower = [-1.0, -1.0, -1.0]
upper = [0.05715, 0.1143, 1.0]
```

A region with `alpha = 0.0` cuts fluid out again, for example the slot of Zalesak's disk. A `field` region takes an expression, `fraction = "y < 0.02 + 0.005*sin(100*x)"`.

**2. Set up a transient run with gravity.**

```toml
[gravity]
vector = [0.0, -9.81, 0.0]
[time]
time_step = 2.5e-4
steps = 1200
```

The pressure written is the reduced one, with the reference density `[fluid] density`. α takes sub-steps of its own (`max_courant`, default 0.25), so the time step is set by the momentum, not by α.

**3. Read the results.**
- `vof.csv` has, per step:
  - the volume of phase 1 and what crossed the boundaries;
  - the clipped volume;
  - the sub-steps;
  - the range of α;
  - the extent of the cells with α > 0.5 (the front of a dam break is `x_max`).
- The field files carry `volume_fraction` and `density`.
- VOF003 warns when more than 10⁻⁶ of the volume had to be clipped.

**Inflow of a phase.** `[vof.inflow] inlet = 1.0` lets phase 1 in through the boundary `inlet`. Other inflows bring α = 0.

**Not yet (VOF002):**
- steady runs, turbulence, a variable density of `[energy]`, scalars;
- regions, sprays, immersed bodies;
- the coupled algorithm, rotating zones, moving meshes;
- periodic boundaries, ranks and restarts.

Energy and phase change: [recipe 358](#vof-phase-change).

---

<a id="vof-phase-change"></a>
## 358. Let a liquid evaporate or boil

**1. Give both phases their heat capacity and conductivity, and switch on the energy equation.** The start temperature is never uniform in a boiling case: give the profile the interface sits in (the cases take the exact one).

```toml
[energy]
initial = "if(x < 1.0e-4, 383.1243 - 10.0 * erf(x / 4.67e-5) / erf(0.0677), 373.1243)"

[[vof.phase]]             # α = 1, the liquid
name = "water"
density = 958.372729
viscosity = 2.81660968e-4
specific_heat = 4216.61269
conductivity = 0.677207143

[[vof.phase]]             # α = 0, the vapour
name = "vapour"
density = 0.597623116
viscosity = 1.22312654e-5
specific_heat = 2077.39017
conductivity = 2.45677073e-2
```

The values are water and steam at 1 atm from `foehn steam --pressure 101325`.

**2. Say which phase is the liquid and where it boils.**

```toml
[vof.phase_change]
liquid = "water"                    # default: the denser phase
saturation_temperature = 373.1243   # K
latent_heat = 2.25654075e6          # J/kg
# or, for water: pressure = 101325.0 (IAPWS-IF97 gives both)
# relaxation = 100.0                # how hard the interface is held at T_sat
```

The surface tension defaults to IF97 at T_sat when `[vof] surface_tension` is not given; set it to 0 for the planar and spherical checks.

**3. Give the vapour a way out.** The vapour needs 1604 times the volume of the water it comes from, so the case needs an outlet (with `[vof.inflow]` for what may flow back in) or a free surface. A closed box would compress nothing: the pressure rows keep the volume.

```toml
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 101325.0
[vof.inflow]
xmax = 1.0
[solver]
algorithm = "piso"
outer_iterations = 10      # the volume source and the density jump need them
```

**4. Read the results.**
- `vof.csv` adds the evaporated mass (net, kg), the rate of the step (kg/s) and the number of interface cells.
- The field files carry `phase_change_rate` (kg/(m³ s), evaporation positive) and `temperature`.
- The check of the run is the volume account: the net boundary volume flow against the volume the phase change made.

**For spherical problems** `examples/meshes/sphere_sector.py` writes a thin sector of the sphere (planar side faces through the centre, `symmetry`), as `vof_scriven.toml` uses it.

**Not yet:** T_sat is fixed (no shift with the pressure or the curvature, no interface resistance), and the heat account of `heat.csv` does not close in a phase-change case (the volume and mass accounts do).

---

<a id="two-fluid"></a>
## 368. Put dispersed bubbles into a liquid (two-fluid model)

The two-fluid (Euler–Euler) model solves a momentum equation for the liquid and one for the gas, with one pressure (L3, MP3, D352 to D356). It is for bubbly flow, where the bubbles are smaller than the cells; for a resolved interface use `[vof]` (recipe 317).

**1. The liquid is the `[fluid]` of the case; add the gas and its bubble size.** Gravity comes from `[gravity]`.

```toml
[fluid]
density = 998.2
viscosity = 1.0e-3

[gravity]
vector = [-9.81, 0.0, 0.0]   # the pipe axis points up

[two_fluid]
surface_tension = 0.072      # N/m: Eötvös number of drag and lift
[two_fluid.gas]
name = "air"
density = 1.2
viscosity = 1.8e-5
diameter = 3.0e-3            # m, one bubble size
```

**2. Say where the gas enters.** The gas velocity defaults to the liquid's plus the terminal speed of the bubbles against gravity; give `velocity` to set it. The liquid inlet velocity is that of the liquid phase, so a superficial velocity j_L becomes j_L/(1 − α).

```toml
[two_fluid.inflow.xmin]
alpha = 0.0717               # gas fraction at the inlet

[boundary.xmin]
type = "velocity-inlet"
velocity = [1.1635, 0.0, 0.0]   # j_L/(1 − α) for j_L = 1.08 m/s
```

**3. Choose the forces (defaults shown).** Drag is always Ishii & Zuber. With a turbulence model the dispersion and the bubble-induced viscosity are on.

```toml
lift = "tomiyama"            # or "constant" with lift_coefficient, or "none"
lift_wall_damping = true     # no lift within a bubble radius of the wall
wall_lubrication = "hosokawa"   # or "none"
turbulent_dispersion = 1.0   # C_TD of Burns et al. (0: off)
bubble_viscosity = 0.6       # C_mu,b of Sato (0: off)
```

**4. Run it transient, with BDF1 and PISO.** Steady runs, BDF2 and an adaptive step are refused (TWF002), as are VOF, energy, scalars, sprays, regions, moving meshes and periodic boundaries in this version.

```toml
[solver]
algorithm = "piso"
outer_iterations = 10

[time]
scheme = "bdf1"
time_step = 2.0e-3
steps = 2000
```

**5. Read the results.** The field files carry `volume_fraction` (the gas fraction α_g) and `gas_velocity_x/y/z` next to the liquid `velocity`. A line in `[[output.lines]]` adds the columns `alpha_gas`, `u_gas`, `v_gas`, `w_gas`.


**Mesh and inlet.** The model averages over the bubbles, so keep the cells at about half a bubble or larger: 1.2 mm cells for 3 mm bubbles in the pipe of [gallery 368](../gallery/multiphase.md#two-fluid-pipe). A mesh graded to cells ten times smaller than the bubbles piles the gas into the wall cells. Give a pipe inlet the turbulence of a developed flow (`viscosity_ratio = 100` at Re ≈ 60 000), or the core stays laminar over tens of diameters and the dispersion is too weak.
