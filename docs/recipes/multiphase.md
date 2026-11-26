# Two fluids (volume of fluid)

Part of the [Recipes](../recipes.md): two immiscible fluids with a sharp interface (L3, volume of fluid). The keys are on the [key page](../cases/flow-keys.md) (`vof.*`).

| # | Recipe |
|---|---|
| 317 | [Put water and air into one case](#vof-two-fluids) |

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
- steady runs, turbulence, energy, scalars;
- regions, sprays, immersed bodies;
- the coupled algorithm, rotating zones, moving meshes;
- periodic boundaries, ranks and restarts.

Energy and phase change follow in MP2.
