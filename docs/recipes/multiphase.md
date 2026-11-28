# Two fluids (volume of fluid)

Part of the [Recipes](../recipes.md): two immiscible fluids with a sharp interface (L3, volume of fluid). The keys are on the [key page](../cases/flow-keys.md) (`vof.*`).

| # | Recipe |
|---|---|
| 317 | [Put water and air into one case](#vof-two-fluids) |
| 358 | [Let a liquid evaporate or boil](#vof-phase-change) |

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
