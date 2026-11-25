# Moving meshes and solids

Part of the [Example gallery](../gallery.md): meshes whose nodes move in time, and later solids moved by the fluid (L22, L52, L53; [moving solids plan](../assessments/moving-solids-reed-valves-plan.md), D286 to D292).

| # | Example | Stage | Feature |
|---|---|---|---|
| 296 | [A piston compresses a closed gas volume on a deforming mesh](#deforming-mesh) | MS0 (D91, L22) | `[mesh_motion]`, swept volumes, moving walls, closed volume |

---

<a id="deforming-mesh"></a>
## 296. A piston compresses a closed gas volume on a deforming mesh

<img src="../figures/gallery/motion_piston.png" alt="Three panels of the piston. Left: the pressure p₀ against the volume over one crank period on logarithmic axes; the adiabatic run lies on the grey isentrope p₀V^γ from 101 kPa at full volume to 705 kPa at a quarter of it and back; the run with walls at 300 K stays below it, reaches 572 kPa and ends at 90 kPa, an open loop. Middle: the temperature of the cooled run at the top dead centre on the compressed mesh, a warm core near 500 K and cool layers along the four walls, with white velocity arrows pointing towards the piston at the top. Right: the error of p₀ against the time step on logarithmic axes, at the top dead centre and after the whole period, both parallel to a dashed Δt² line." width="800" height="235">

The mesh deforms with a node displacement given as an expression of the
node's place in the mesh file and of time (`[mesh_motion] displacement`,
D53). Every step the solver moves the nodes, recomputes the geometry and
gives every face its swept volume δV_f, the volume the face sweeps on its
way from the old to the new nodes (D91). The swept volumes of a cell close
exactly: Σ ±δV_f = V^(n+1) − V^n to round-off, on hexahedra, tetrahedra,
polyhedra and mixed meshes (the space conservation law, Demirdžić & Perić
1988). Convection uses the flux relative to the moving face, the time terms
use the volume of each time level (V^n and V^(n−1) for BDF2), and walls move
with their nodes.

`motion_piston_adiabatic.toml` is the evenly compressed cylinder of L22: a
gas column 2 cm wide and 4 cm high, closed by adiabatic walls; the top wall
is a piston that moves with a crank down to a quarter of the height and back
in 1 s, and every node moves in proportion to its height. The gas is the
closed volume of stage 7c (low-Mach variable density, p₀ from its fixed mass
and the volume of the moved mesh, heated by the pressure work V dp₀/dt). A
slow adiabatic compression is isentropic, p₀ V^γ = const with γ = 1.39983:

| | p₀ (kPa) | exact (kPa) | error |
|---|---|---|---|
| top dead centre (t = 0.5 s, V = V₀/4), 100 steps per period | 705.03 | 705.50 | 6.6·10⁻⁴ |
| after the period (t = 1 s) | 101.323 | 101.325 | 1.8·10⁻⁵ |

The error falls with Δt² (BDF2): halving the step from 0.02 s down to
0.0025 s divides it by 4.13, 4.06 and 4.02 at the top dead centre, and
after the period it falls from 1.5·10⁻⁴ to 2.9·10⁻⁷. The mass stays to
10⁻¹⁵ and the heat book of every step closes to 10⁻¹⁴.

`motion_piston_cooled_walls.toml` holds all walls and the piston at 300 K
(40 × 80 cells, 30 s). Thermal boundary layers grow along the walls; at the
top dead centre the core is near 500 K and p₀ is 572 kPa instead of 706 kPa.
The heat the walls take out opens the p–V loop: the expansion ends at
90.4 kPa, 11 % below the start.

<img src="../figures/gallery/motion_free_stream.png" alt="Two panels of the free-stream check. Left: a channel mesh of 20 by 10 cells whose inner nodes swing in a wavy pattern at the largest swing, coloured by the logarithm of the deviation of the velocity from the uniform 1 m/s; the cells show a random mix of values between 10⁻¹⁷ and 10⁻¹⁵ without a pattern. Right: the momentum and continuity residuals of every step over 0.5 s, flat between 10⁻¹⁶ and 3·10⁻¹⁶, with the worst closure of the space conservation law as a dashed line near 10⁻¹⁵." width="800" height="258">

`motion_free_stream.toml` is check MS0.1: a uniform flow of 1 m/s through a
channel on a perturbed mesh whose inner nodes swing by up to 3 cm (cell
volumes between 0.77 and 1.25 of their size in the mesh file). On the
moving mesh the flow stays uniform to round-off: the velocity deviates by at
most 1.3·10⁻¹⁵ m/s, the residuals stay below 3·10⁻¹⁶. The test does the same
on tetrahedral and polyhedral meshes and with BDF1, BDF2 and Crank–Nicolson.

**Limits of stage MS0.** The motion is prescribed; a mesh-motion solver and
bodies moved by the fluid are MS1 and MS2. The geometry and the solver setup
are rebuilt every step (an AMG setup per step). A restart, several ranks,
turning zones, regions, sprays, immersed bodies, periodic boundaries, an
error-controlled time step and level sweeps are refused (MOV001); only walls
move with the mesh, other boundaries may only slide in their plane (MOV003);
a closed volume of constant density cannot be compressed (MOV004). Turbulence
models and passive scalars run on the moving mesh but are not yet checked
against a reference.

Reproduce (the adiabatic piston and the free stream take seconds, the cooled
piston half a minute, the sweep about a minute):

```sh
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/motion_piston_adiabatic.toml --out output/adiabatic
target/release/foehn flow examples/cases/flow/motion_piston_cooled_walls.toml --out output/cooled
for dt in 0.02 0.01 0.005 0.0025; do for end in 0.5 1.0; do
  target/release/foehn flow examples/cases/flow/motion_piston_adiabatic.toml \
    --out output/sweep/dt${dt}_end${end} --set time.time_step=$dt --set time.end_time=$end
done; done
target/release/foehn flow examples/cases/flow/motion_free_stream.toml --out output/free \
  --set 'output.fields_times=[0.125]' --set output.history_every=1
pip install matplotlib h5py numpy
python3 -m tools.plot.mesh_motion docs/figures/gallery output/adiabatic output/cooled output/sweep output/free
```

Recipe: [296](../recipes/moving-solids.md#deforming-mesh). Case files: `examples/cases/flow/motion_piston_adiabatic.toml`, `motion_piston_cooled_walls.toml`, `motion_free_stream.toml`.
