# Moving meshes and solids

Part of the [Recipes](../recipes.md): move the mesh in time, and later move solids with the fluid (L22, L52, L53). The keys are on the [key page](../keys.md) (`mesh_motion.*`).

| # | Recipe |
|---|---|
| 296 | [Move the mesh with a prescribed displacement (a piston)](#deforming-mesh) |

---

<a id="deforming-mesh"></a>
## 296. Move the mesh with a prescribed displacement (a piston)

**Add a `[mesh_motion]` table to a transient case.** The displacement of
every node is an expression of its place in the mesh file (x, y, z) and of
the time t, with the constants of the case:

```toml
[constants]
compression = 4.0     # V₀/V at the top dead centre
period = 1.0

[mesh_motion]
# every node moves down in proportion to its height: an evenly compressed
# cylinder whose top wall (the piston) follows a crank
displacement = [
  0.0,
  "-y * (1 - 1 / compression) * (1 - cos(2 * pi * t / period)) / 2",
  0.0,
]
```

The displacement at the start time may be non-zero; the mesh starts there.
Every step the nodes move to the end of the step, the geometry is computed
again and every face gets its swept volume, so that the volumes and the
fluxes agree exactly (the space conservation law, D91).

**Boundaries.** Walls move with their nodes: the wall velocity of the case
gives the tangential part, the mesh motion the normal part. Every other
boundary (inlets, outlets, symmetry planes, wedges) must stay in its plane;
its nodes may slide along it. A boundary that sweeps a volume stops the run
with MOV003 and its name.

**A closed volume needs a density that follows the pressure.** Compressing a
domain without a pressure boundary works only with the low-Mach variable
density of a closed volume (stage 7c): p₀ comes from the fixed mass and the
volume of the moved mesh, and the pressure work V dp₀/dt heats the gas.

```toml
[energy]
initial = 300.0
density = "variable"
thermodynamic_pressure = "closed-volume"
```

With a constant density the run stops with MOV004 at the first step that
changes the volume. A domain with a pressure boundary may change its volume
with any density: the flow through the boundary takes up the difference.

**Run and read.**

```sh
target/release/foehn flow examples/cases/flow/motion_piston_adiabatic.toml --out output/piston
```

The log says that the mesh deforms. `history.csv` has a column `volume`
(m³) for every step, and `result.json` has, under `transient.mesh_motion`,
the volume at the end, the smallest and largest ratio of a cell volume to
its volume in the mesh file over the run, and `closure_max`, the worst
closure of the space conservation law of a cell (round-off, about 10⁻¹⁵).
The fields of `fields_<step>.vtkhdf` are written on the moved mesh. Mass
flows in `boundaries.csv` and the balances are taken relative to the
moving faces, so a moving wall carries no mass.

**Choose the step.** The error of a slow compression falls with Δt² with
BDF2 (6.6·10⁻⁴ in p₀ at the top dead centre with 100 steps per crank
period, [gallery 296](../gallery/moving-solids.md#deforming-mesh)). A
displacement that folds a cell stops the run with MOV002, the time and the
place.

**Limits (stage MS0).** The motion is prescribed; a mesh-motion solver,
cell layering and bodies moved by the fluid come with MS1 and MS2. The
solver setup is rebuilt every step. Restarts, several ranks, turning zones,
regions, sprays, immersed bodies, periodic boundaries, an error-controlled
time step and level sweeps are refused (MOV001).
