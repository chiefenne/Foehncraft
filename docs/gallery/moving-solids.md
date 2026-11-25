# Moving meshes and solids

Part of the [Example gallery](../gallery.md): meshes whose nodes move in time, and later solids moved by the fluid (L22, L52, L53; [moving solids plan](../assessments/moving-solids-reed-valves-plan.md), D286 to D292).

| # | Example | Stage | Feature |
|---|---|---|---|
| 296 | [A piston compresses a closed gas volume on a deforming mesh](#deforming-mesh) | MS0 (D91, L22) | `[mesh_motion]`, swept volumes, moving walls, closed volume |
| 307 | [A breathing channel on three ranks, continued from a restart file](#deforming-mesh-ranks) | MS0 part 2 (D91, D93) | `[run] rank_count`, restart files of a moving mesh, GPU check |

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
are rebuilt every step (an AMG setup per step). Restart files and ranks came
with part 2 ([gallery 307](#deforming-mesh-ranks)); turning zones, regions, sprays, immersed bodies, periodic boundaries, an
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

---

<a id="deforming-mesh-ranks"></a>
## 307. A breathing channel on three ranks, continued from a restart file

<img src="../figures/gallery/motion_channel_ranks.png" alt="Three panels of the breathing channel. Left: the speed in a channel 2 m long and 0.5 m high on a mesh whose inner nodes are swung into gentle waves, slow blue and green layers along the two walls and a fast red core near 1.4 m/s; two white wavy lines cross the channel at about x = 0.67 m and 1.33 m, the cuts between the three ranks. Middle: a yellow stripe of dye entering at the middle of the inlet and spreading slightly as it reaches the outlet. Right: the axial velocity at a probe in the middle against time, rising from 1.0 to 1.38 m/s; a thick blue line for one rank, a thin orange line for three ranks on top of it, and black dots from a dashed restart line at 0.75 s onwards, all on the same curve." width="800" height="185">

Stage MS0 part 2 makes the deforming mesh of [296](#deforming-mesh) run on
ranks and restart, with the same bitwise rules as a fixed mesh (D42, D93).

`motion_channel_ranks.toml` is a channel 2 m × 0.5 m between two walls at
Re = 100 with a plug inflow of 1 m/s and a stripe of dye. Its inner nodes
swing back and forth by up to 3 cm in both directions (6 periods in 1.5 s,
cell volumes between 0.80 and 1.21 of their size in the mesh file); the
walls, the inlet and the outlet stay in place. The case runs on three
ranks (`[run] rank_count = 3`):

- **Ranks.** Every rank moves the whole mesh, as with the sliding mesh by
  imprinting ([gallery 305](meshes.md#run-ranks-imprint)). The partition of the start stays
  for the whole run (the cells keep their numbers, the white cuts move with
  the nodes); the geometry and the swept volumes of each rank's part are
  taken from those of the whole mesh. Every rank therefore solves bitwise
  the cells of the one-rank run, and rank 0 writes every file (fields,
  history, probes, balances, walls, restart files) byte for byte as one
  rank does.
- **Restart.** The restart file of a deforming mesh holds the nodes of the
  current time level, the cell volumes Vⁿ and Vⁿ⁻¹, the volume each face
  swept in the last step (BDF2 needs it) and the range of the volume
  ratios, as the optional geometry datasets of D91. A run continued from
  `restart_000075.h5` (t = 0.75 s, written on three ranks) starts on the
  mesh of that time and ends bitwise as the run without the interruption,
  on two ranks and on one. A restart whose nodes the displacement does not
  put there (a changed `[mesh_motion]`) is refused (MOV001).
- **GPU.** The grid flux, the relative flux and the volumes of the old
  levels are operations of the flow device, which runs the same kernel
  bodies on the CPU and the GPU (D31). CI has no GPU: the notebook
  `tools/notebooks/motion_colab.ipynb` (`python3 -m tools.bench.motion_cloud`)
  runs the free stream, this channel with its restart and the adiabatic
  piston on the CPU and on a Colab GPU and compares them (to rounding,
  since the CPU sums exactly and the GPU keeps its block sums, MPI step W2).

**Limits.** Every rank moves the whole mesh and rebuilds its part every
step: on this small mesh three ranks take about five times as long as one
(a scalable motion of the parts comes with the mesh-motion solver, MS1).
The variable density of stage 7c has no restart yet, so the closed-volume
piston does not restart and does not run on ranks before that work (the
variable density on ranks is PR #303). On the GPU the solver is built again
every step, with its device setup.

Reproduce (seconds):

```sh
cargo build --release -p foehncraft-cli
F=examples/cases/flow/motion_channel_ranks.toml
target/release/foehn flow $F --ranks 1 --out output/one --set 'output.fields_times=[1.0625]'
target/release/foehn flow $F --out output/three --set 'output.fields_times=[1.0625]'
target/release/foehn flow $F --ranks 2 --restart output/three/restart_000075.h5 --out output/cont
cmp output/one/fields.vtkhdf output/cont/fields.vtkhdf && echo bitwise
pip install matplotlib h5py numpy
python3 -m tools.plot.mesh_motion_ranks docs/figures/gallery output/one output/three output/cont
```

Recipe: [307](../recipes/moving-solids.md#deforming-mesh-ranks). Case file: `examples/cases/flow/motion_channel_ranks.toml`.
