# Flow cases (`examples/cases/flow/`)

Part of the [Case catalogue](../cases.md): incompressible laminar flow, steady and transient, grouped by physics; what a flow run writes.

On this page:

- [Steady laminar flow](#steady-laminar-flow)
- [Transient flow](#transient-flow)
- [Runs driven by programs, balances](#runs-driven-by-programs-balances)
- [The flow loop on a GPU](#the-flow-loop-on-a-gpu)
- [Boundary conditions](#boundary-conditions)
- [Derived fields](#derived-fields)
- [Expressions and materials](#expressions-and-materials)
- [Passive scalars](#passive-scalars)
- [Energy equation](#energy-equation)
- [Gravity and buoyancy](#gravity-and-buoyancy)
- [Lagrangian sprays](#lagrangian-sprays)
- [Porous regions](#porous-regions)
- [Immersed bodies](#immersed-bodies)
- [Deforming meshes](#deforming-meshes)
- [Species](#species)

Stage 4: steady incompressible laminar flow, SIMPLEC with the Rhie–Chow
interpolation ([explanation](../verification/stage4.md)). Stage 5: transient
flow with PISO, BDF2 or Crank–Nicolson, when the case has a `[time]` table
([explanation](../verification/stage5.md)). Run with
`foehn flow FILE`; `foehn flow FILE --check` only checks the case file
and the boundary conditions (D30); command-line options override the file
(`foehn flow --help`). Each run writes `flow.log` (the header says what
is written and how large it is), `residuals.csv` (one row per outer
iteration; `residuals_level<N>.csv` per level of a sweep), `boundaries.csv`
(per boundary: area, mass flow, mean pressure, force), `wall_<boundary>.csv` (per
wall face: position, pressure, wall shear stress), `line_<name>.csv` for
each sampling line, `fields.vtkhdf` (velocity, speed, pressure and, with an
exact solution, the velocity error; in blocks of 65,536 cells that ParaView
opens as one data set, MurmelCraft reads in parts and scripts read in one
piece with `foehncraft.vtkhdf.internal`), `case.resolved.toml` (the case with
every default filled in), `flow.md` and `flow.json` (the report),
`result.json` (the result summary for programs, without timings) and
`run.jsonl` (the run record, one JSON object per row; both in
[interfaces.md](../interfaces.md)), and the
charts `residuals.svg`, `convergence.svg` (sweeps) and `centreline_u.svg`,
`centreline_v.svg` (cavity against the reference data). The run fails when
a check of the case fails (observed order, error, reference deviation,
convergence, mass balance). CI runs every case on Linux and the quick ones
on macOS and Windows (R5); the sets are in `tools/check/ci_cases.py`, and the
long cylinder at Re 100 runs on main, nightly and with the label `full-ci`.

A transient case also writes `history.csv` (one row per time step: time,
Δt, Courant number, residuals, outer and CG iterations; the residuals of a
steady case stay in `residuals.csv`), `forces.csv` (force per wall and the
coefficients, every `output.monitor_every` steps or `monitor_interval`
seconds), `probes.csv`, `fields_<step>.vtkhdf` on the schedule of
`output.fields_every`, `fields_at`, `fields_interval` and `fields_times`
(the last two in simulated time, landed on exactly), `restart.h5` (D42), `coefficients.svg` (force coefficients;
for the periodic cylinder the evaluation window) and, for a time-order
sweep, `time_convergence.svg`. Control files in the output directory steer
a running case (D40): `.stop` ends it with its restart file and results,
`.result_N` writes the fields at step (or iteration) N, `.read_toml_N`
re-reads the case file at N and applies what may change during a run.

## Steady laminar flow

Couette and Poiseuille flow, the lid-driven cavity and the steady cylinder.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/couette_hex.toml` | generated hex 2D, 128 cells | plane Couette flow with convection: error at the solver tolerance (≤ 1e-9 required) | [recipe 18](../recipes/flow.md#flow-verification) |
| `flow/couette_tet.toml`, `couette_poly.toml`, `couette_mixed.toml` | generated 2D, perturbed by 15 %, 128 to 256 cells | the same Stokes flow on perturbed tetrahedra (prisms), polyhedra and mixed cells | [recipe 18](../recipes/flow.md#flow-verification) |
| `flow/poiseuille_hex.toml` | generated hex 2D, perturbed, 512 to 8 192 cells (3 levels) | plane Poiseuille flow, observed order of the velocity error 1.8–2.2 required | [gallery 13](../gallery/flow.md#poiseuille), [recipe 18](../recipes/flow.md#flow-verification) |
| `flow/poiseuille_tet.toml`, `poiseuille_poly.toml`, `poiseuille_mixed.toml` | generated 2D, 3 levels | the same on prisms, perturbed polyhedra and mixed cells | [gallery 13](../gallery/flow.md#poiseuille), [recipe 18](../recipes/flow.md#flow-verification) |
| `flow/cavity_re100.toml` | generated hex 2D, 64 × 64 | lid-driven cavity, Re = 100, against Ghia et al. (1982), deviation ≤ 2 % of the lid speed | [gallery 12](../gallery/flow.md#cavity), [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re1000.toml` | generated hex 2D, 128 × 128 | the same at Re = 1000 | [gallery 12](../gallery/flow.md#cavity), [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re1000_coupled.toml` | generated hex 2D, 128 × 128 | the cavity at Re = 1000 with the coupled algorithm and the automatic pseudo-time step: about 450 iterations instead of 1 397, same result | [gallery 78](../gallery/flow.md#coupled), [gallery 80](../gallery/flow.md#momentum-multigrid), [recipe 78](../recipes/flow.md#coupled) |
| `flow/adjoint_residual_cavity.toml` | generated polygons 2D, 24 × 24 | the steady residual R(q) of the adjoint (L17 stage A0) along a SIMPLEC run and a coupled run: both drive it to round-off | [gallery 271](../gallery/verification.md#adjoint-a0), [recipe 271](../recipes/verification.md#adjoint-a0) |
| `flow/cavity_re100_tet3d.toml` | generated tet 3D, 24 576 tetrahedra in a slab with symmetry planes | the cavity on fully three-dimensional tetrahedra, deviation ≤ 3 % | [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re100_poly3d.toml` | generated poly 3D, 5 445 polyhedra | the same on polyhedra | [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/cavity_re1000_poly3d.toml` | generated poly 3D, 28 227 polyhedra | Re = 1000 on polyhedra, deviation ≤ 3 % | [gallery 12](../gallery/flow.md#cavity), [recipe 19](../recipes/flow.md#flow-cavity) |
| `flow/taylor_couette_mrf.toml` | built-in annulus mesh, 32 × 8, 64 × 16 and 128 × 32 cells (3 levels), the inner half of the radial layers is the zone `rotor` | laminar Taylor–Couette flow at Re = 50 against the analytic solution with the inner half in a rotating frame (D55): L2 error 3.5, 0.69 and 0.13 %, observed order 2.3 and 2.4 | [gallery 94](../gallery/flow.md#mrf), [recipe 94](../recipes/flow.md#mrf) |
| `flow/taylor_couette_sliding.toml` | built-in annulus mesh, 16 × 8 to 128 × 64 cells (4 levels; the case file holds 32 × 16), the inner half of the radial layers is the zone `rotor` | laminar Taylor–Couette flow at Re = 50 from rest to the steady state, the rotor advancing by one of 32 sectors per time step (D56, D182): L2 error 6.8, 1.35, 0.28 and 0.06 %, observed order 2.2 to 2.3 | [gallery 95](../gallery/flow.md#sliding-mesh), [recipe 95](../recipes/flow.md#sliding-mesh) |
| `flow/taylor_couette_imprint.toml` | built-in annulus mesh with a split interface, 32 rotor and 48 stator sectors, 640 cells (the case file; 160 to 2 560 cells in the test) | laminar Taylor–Couette flow at Re = 50 from rest to the steady state with a time step of 0.02 s that no sector decides, the interface imprinted in every step (D182, variant B): L2 error 6.7, 1.29 and 0.25 %, observed order 2.4 | [gallery 96](../gallery/flow.md#sliding-mesh-imprint), [recipe 96](../recipes/flow.md#sliding-mesh-imprint) |
| `flow/couette_two_rotors.toml` | built-in annulus mesh with two split interfaces (`outer_rotor_layers`): an inner rotor ring of 32 sectors, a stator ring of 48 and an outer rotor ring of 40, 960 cells | laminar Couette flow between two cylinders turning at +10 and −5 rad/s, each cylinder a sliding rotor by imprinting with its own speed (D185), from rest to the steady state (500 steps of 0.04 s, 11 s): L2 error 1.8 %, the two torques 0.8 % from the exact value, the shaft powers add up to the dissipation | [gallery 102](../gallery/flow.md#several-rotors), [recipe 102](../recipes/flow.md#several-rotors) |
| `flow/axial_stage_uniform.toml` | built-in mesh of revolution (`[mesh.revolved]`): an annular duct with a rotor block on the hub, interface disc–cylinder–disc, 24 rotor and 34 stator intervals around the axis, 2 090 cells | uniform axial flow between slip walls through a sliding rotor at 40 rad/s (D186), 20 steps: L2 error 0.23 % from the step between the hub polygons of rotor and stator (second order in the pitch); with `r = [0.0, 0.04, 0.10]` a pot through the axis, 1.5·10⁻⁶ | [gallery 103](../gallery/flow.md#revolved-interfaces), [recipe 103](../recipes/flow.md#revolved-interfaces) |
| `flow/pot_solid_body.toml` | built-in mesh of revolution through the axis: a vessel with a rotor pot whose two discs cross the axis, 16 rotor and 24 stator intervals | solid-body rotation with every wall turning at the rotor's speed, kept by the imprinted meshes through half a revolution (D186): L2 error 0.8 % | [gallery 103](../gallery/flow.md#revolved-interfaces), [recipe 103](../recipes/flow.md#revolved-interfaces) |
| `flow/pot_solid_body_viscous.toml` | as `pot_solid_body.toml` | the same solid-body rotation at ten times the viscosity (ν = 10⁻² m²/s): the regression case of the PISO corrector on the wedge cells at the axis (PR #106), whose coupling around the axis is about 90 times their time term: L2 error 0.11 % (2.0 with Issa's corrector) | [gallery 103](../gallery/flow.md#revolved-interfaces), [recipe 103](../recipes/flow.md#revolved-interfaces) |
| `studies/mixer/mrf.toml`, `sliding.toml`, `mrf_turbulent.toml`, `sliding_turbulent.toml` (in `examples/studies/`) | Gmsh mesh of `examples/meshes/mixer.py` (`pip install gmsh`), 16 326 cells, conformal for the frame, split with 64 rotor and 96 stator intervals for the sliding mesh | a baffled stirred vessel in 2D, four blades turning at 10 rad/s, four baffles: steady in a rotating frame at Re = 100 (972 iterations, 22 s), and from rest through three revolutions with the sliding mesh by imprinting at Re = 400 (1 508 steps, 7 min); field output for the figures of `studies/mixer.py`, which runs the same case files on the vessel with five blades and three baffles (`mixer.py --blades 5 --baffles 3`, 16 848 and 16 897 cells; the figures show no four-blade rotor); the frame case also writes the fields of rotating machinery (relative and frame velocity, cylindrical components, flow angles, relative total pressure, rothalpy, D183); the `_turbulent` cases are the same vessel with water (Re = 20 000) and k-ω SST, the frame with the coupled solver (337 iterations, 27 s) and the sliding mesh through three revolutions (1 508 steps, 10 min) | [gallery 98](../gallery/flow.md#stirred-vessel), [gallery 99](../gallery/flow.md#rotor-fields), [gallery 100](../gallery/flow.md#turbulent-rotor), [gallery 101](../gallery/flow.md#torque-power), [recipe 98](../recipes/flow.md#stirred-vessel), [recipe 99](../recipes/flow.md#rotor-fields), [recipe 100](../recipes/flow.md#turbulent-rotor), [recipe 101](../recipes/flow.md#torque-power) |
| `studies/mixer/mrf_transient.toml`, `studies/twin_mixer/sliding.toml` (in `examples/studies/`) | Gmsh meshes of `examples/meshes/mixer.py` (conformal, 16 304 cells) and `examples/meshes/twin_mixer.py` (a stadium vessel with two rotors, each with its own interface circle of 64 rotor and 96 stator intervals, 16 953 cells) | the stirred vessel of `sliding.toml` as a transient frozen rotor (D185: `[[mrf]]` with `[time]`, 480 steps, 0.6 s) for the comparison with the sliding mesh, and a vessel with two stirrers turning at +10 and −15 rad/s, two sliding rotors by imprinting (800 steps of 1 ms); the figures run both on the five-blade meshes (`twin_mixer.py --right-blades 5`, three and five blades, 17 775 cells) | [gallery 102](../gallery/flow.md#several-rotors), [recipe 102](../recipes/flow.md#several-rotors) |
| `studies/three_mixer/sliding.toml` (in `examples/studies/`) | Gmsh mesh of `examples/meshes/three_mixer.py` (`pip install gmsh`): a stadium vessel with three rotors, each with its own interface circle of 64 rotor and 96 stator intervals; rotor C has seven circular-arc blades of constant thickness, 26 060 cells | three sliding rotors by imprinting at +10, −15 and +20 rad/s from rest (800 steps of 0.5 ms) | [gallery 225](../gallery/flow.md#curved-blades), [recipe 225](../recipes/flow.md#curved-blades) |
| `flow/cylinder_re20.toml` | built-in cylinder mesh, refine 1: 24 576 hexahedra | Schäfer & Turek 2D-1, steady, Re = 20: C_D, C_L, Δp within 0.5 % of the reference ranges | [gallery 15](../gallery/flow.md#cylinder), [recipe 23](../recipes/flow.md#cylinder) |
| `flow/naca0012_su2_re1000.toml` | SU2 from PyAero (`meshes/naca0012.su2`), 6 420 hexahedra in one layer | laminar NACA 0012, Re = 1000, 4°, coupled solver, force coefficients: a flow on an SU2 mesh | [gallery 127](../gallery/meshes.md#su2-airfoil), [recipe 127](../recipes/meshes.md#mesh-su2) |
| `flow/cylinder_potential.toml` | built-in cylinder mesh, refine 1: 24 576 hexahedra | the potential-flow start (D52) around the cylinder, written without outer iterations | [gallery 20](../gallery/flow.md#time-control), [recipe 28](../recipes/flow.md#time-control) |
| `flow/cylinder_re20_gci.toml` | built-in cylinder mesh, refine 0 to 2: 6 144 to 98 304 hexahedra (3 levels) | the benchmark 2D-1 on three meshes: observed order, extrapolated value and grid convergence index of C_D, C_L and Δp (D67); every mesh within 8 % of the reference ranges; coupled algorithm, 51 to 225 outer iterations (about a minute, in the `long` set) | [gallery 25](../gallery/verification.md#gci), [recipe 34](../recipes/verification.md#gci) |
| `flow/tjunction_mixing.toml` | Gmsh tetrahedra of `meshes/tjunction_gci.geo` (h = 0.0126 m, 4 852 polyhedra; the study adds 16 665 and 44 738) converted into polyhedra | laminar T-junction (Re 200 main pipe, 100 branch) with a passive scalar from the branch; three meshes with GCI: inlet pressure of the main pipe order 1.77, GCI 2.4 %; the scalar with the MLP limiter (D237), its difference across the pipe order 2.1, GCI 11 %; validation against published measurements pending | [gallery 196](../gallery/flow.md#tjunction-mixing), [recipe 196](../recipes/flow.md#tjunction-mixing) |
| `flow/diffuser_channel_profile.toml` | generated hex 2D, 2 × 96 cells, graded towards both walls (ratio 60), periodic in x | the periodic SST channel at centreline Re_H 20 000 (U_c = 1 m/s, U_b = 0.893 m/s) that gives the developed inlet table `diffuser_inlet_profile.csv` of the diffuser; `examples/studies/diffuser_inlet_profile.py` runs it and writes the table | [gallery 199](../gallery/turbulence.md#diffuser) |
| `flow/diffuser_sst.toml` | structured quadrilaterals from `meshes/diffuser.py` (SU2, level 0: 6 960 cells, Δx = 0.25 H along the diffuser, graded to y⁺ about 1 at both walls) | the asymmetric plane diffuser of Buice & Eaton (10°, area ratio 4.7, both corners rounded with 4.3 H as in the experiment) at Re_H 20 000 with k-ω SST `low-re` and the developed channel-profile inlet (`diffuser_channel_profile.toml`): potential start, upwind for 200 iterations then linear upwind (`[[schedule]]`), relaxation 0.8; converges in 580 iterations with separation at x = 1.65 H, reattachment at 29.25 H and C_p 0.54 at 45 H; the grid study (`examples/studies/diffuser_gci.py`) runs three levels; no measured value is used yet | [gallery 199](../gallery/turbulence.md#diffuser), [recipe 199](../recipes/turbulence.md#diffuser) |
| `flow/diffuser_sst_sharp.toml` | structured quadrilaterals from `meshes/diffuser.py --sharp` (SU2, level 0: 8 400 cells, diffuser cells shortened to 0.0175 H at the corners, aspect ratio about 11 at the kinks) | the same diffuser with sharp corners and the same settings: converges in 1 087 iterations; the wall shear stays near zero for 1.5 H behind the corner, the separated region runs from 1.82 H to 29.29 H; grid study with `diffuser_gci.py --sharp` on three meshes (reattachment GCI 0.51 %, separation 9 %) | [gallery 199](../gallery/turbulence.md#diffuser), [recipe 199](../recipes/turbulence.md#diffuser) |
| `flow/ubend_channel_profile.toml` | generated hex 2D, 2 × 128 cells, graded towards both walls (ratio 150), periodic in x | the periodic SST channel at Re_b 10⁵ (U_b = 1 m/s, U_c = 1.090 m/s) that gives the developed inlet table `ubend_inlet_profile.csv` of the U-bend; `examples/studies/ubend_inlet_profile.py` runs it and writes the table | [gallery 204](../gallery/turbulence.md#ubend) |
| `flow/ubend_sst.toml` | structured quadrilaterals from `meshes/ubend.py` (SU2, level 0: 8 832 cells, 96 cells over the turn, graded to first-cell y⁺ about 1.9 at both walls) | the two-dimensional U-duct of Monson et al. (1990), Rc/H = 1, at Re_b 10⁵ with k-ω SST `low-re` and the developed channel-profile inlet (`ubend_channel_profile.toml`): potential start, upwind for 200 iterations then linear upwind, relaxation 0.8; converges in 499 iterations with separation on the inner wall at θ = 138° and reattachment 2.3 H behind the turn; the grid study (`examples/studies/ubend_gci.py`) runs three levels; no measured value is used yet | [gallery 204](../gallery/turbulence.md#ubend), [recipe 204](../recipes/turbulence.md#ubend) |
| `flow/ubend_sst_sharp.toml` | Cartesian quadrilaterals from `meshes/ubend.py --sharp` (SU2, level 0: 12 624 cells, square divider end and end wall, square cells at the kinks) | the same duct with sharp edges: relaxation 0.5, linear upwind from iteration 601, relaxation 0.8 from 1 201; converges in 2 479 iterations, the separated region behind the divider end closes at 3.05 H; grid study with `ubend_gci.py --sharp` | [gallery 204](../gallery/turbulence.md#ubend), [recipe 204](../recipes/turbulence.md#ubend) |
| `flow/pipe_gmsh_raw.toml` | Gmsh tetrahedra without optimisation (`examples/meshes/pipe_gmsh_raw.msh`), 1 572 cells, 28 faces with a skewness above 0.9 | developing pipe flow at Re = 100 on a poor mesh: the cell-quality remediation (D50) switches to the Gauss pressure gradient and the limiter, reports it and converges; with `--set solver.remediation=off` the run diverges with a report that names the skewed faces | [gallery 26](../gallery/meshes.md#skewed-tets), [recipe 36](../recipes/meshes.md#skewed-tets) |
| `flow/ranks_channel_poly.toml` | generated poly 2D, 2 821 cells, perturbed by 10 % | a channel whose upper wall moves against the inflow, steady SIMPLEC from potential flow, solved again on 2, 3 and 5 ranks and bitwise the one-rank solution (MPI step W5) | [gallery 234](../gallery/meshes.md#flow-ranks), [recipe 234](../recipes/meshes.md#flow-ranks) |
| `flow/ranks_heated_channel_sst.toml` | generated poly 2D, 2 821 cells, perturbed by 10 % | a turbulent channel (k-ω SST, Re = 10 000) with a heated lower wall and the energy equation, solved again on 2, 3 and 5 ranks and bitwise the one-rank solution in every field (MPI step W5) | [gallery 244](../gallery/meshes.md#flow-models-ranks), [recipe 244](../recipes/meshes.md#flow-models-ranks) |
| `flow/ranks_restart_channel.toml` | generated poly 2D, 2 821 cells, perturbed by 10 % | the channel of `ranks_channel_poly` marched 40 time steps with BDF2, solved again on 2, 3 and 5 ranks; the restart gathered from five ranks must have the bytes of the one-rank restart, and a run continued from it with `--restart` reads it on any rank count (MPI step W6) | [gallery 262](../gallery/meshes.md#restart-ranks), [recipe 262](../recipes/meshes.md#restart-ranks) |
| `flow/ranks_run_channel.toml` | generated poly 2D, 2 821 cells, perturbed by 10 % | the channel of `ranks_channel_poly` solved on 5 ranks as the run itself (`[run] rank_count`, `--ranks N`); rank 0 writes every output file byte for byte as one rank does, plus `ranks.vtkhdf` with the rank of every cell (MPI step W7) | [gallery 279](../gallery/meshes.md#run-ranks), [recipe 279](../recipes/meshes.md#run-ranks) |
| `flow/ranks_run_transient_sst.toml` | generated poly 2D, 2 821 cells, perturbed by 10 % | a heated turbulent channel (k-ω SST, energy, a passive dye) with the upper wall moving against the inflow, marched 180 BDF2 steps on 5 ranks as the run itself; every output file, the heat balance of every step and the turbulence summary byte for byte one rank (MPI step W7 part 2) | [gallery 290](../gallery/meshes.md#run-ranks-physics), [recipe 290](../recipes/meshes.md#run-ranks-physics) |

## Transient flow

Time-stepped cavity, Taylor–Green vortex, periodic cylinder, error-based time step.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cavity_re100_transient.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 time-stepped (PISO, BDF2, Courant number 6) to its steady state: the SIMPLEC solution, deviation from Ghia ≤ 2 % | [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_bdf2.toml` | generated hex 2D, 64 × 64 | Taylor–Green vortex, plain PISO with BDF2 over five time steps: time order 1.8–2.3 required | [gallery 14](../gallery/flow.md#taylor-green), [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_cn.toml` | generated hex 2D, 64 × 64 | the same with Crank–Nicolson | [gallery 14](../gallery/flow.md#taylor-green), [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_space.toml` | generated hex 2D, 16² to 128² (4 levels) | Taylor–Green vortex, space order 1.8–2.2 required | [gallery 14](../gallery/flow.md#taylor-green), [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_cfl10.toml` | generated hex 2D, 64 × 64 | Taylor–Green vortex at Courant numbers 1.3 to 9 with outer iterations: time order 2 | [gallery 14](../gallery/flow.md#taylor-green), [recipe 21](../recipes/flow.md#transient-flow) |
| `flow/taylor_green_error.toml` | generated hex 2D, 64 × 64 | the time step from the BDF2 error estimate with rejected steps (D65), probe rows every 0.05 s and fields every 0.1 s landed on exactly (D66) | [gallery 20](../gallery/flow.md#time-control), [recipe 28](../recipes/flow.md#time-control) |
| `flow/cylinder_re100.toml` | built-in cylinder mesh, refine 0: 6 144 hexahedra | Schäfer & Turek 2D-2, periodic, Re = 100: C_D,max, C_L,max, St, Δp within 2.5 % (about 8 minutes; Linux, on main, nightly and with the label `full-ci`) | [gallery 15](../gallery/flow.md#cylinder), [recipe 23](../recipes/flow.md#cylinder) |

## Runs driven by programs, balances

Limits, exit codes, restart and the energy ledger.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cavity_re100_budget.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 with `[run.limits]` (wall time, iterations, output volume), run the way a script or an agent runs it; a tighter limit stops it with exit code 3 and `--restart latest` continues bitwise | [gallery 18](../gallery/runs.md#run-contract), [recipe 26](../recipes/runs.md#run-contract) |
| `flow/cylinder_steady_diagnosis.toml` | built-in cylinder mesh, refine 1, 24 576 hexahedra | the cylinder at Re = 400 run steady on purpose: the residuals stall, the diagnosis at iteration 2000 says so and stops the run (exit code 3, `RUN017`; Linux, set `long`) | [gallery 166](../gallery/runs.md#budgets), [recipe 166](../recipes/runs.md#budgets) |
| `flow/cavity_re100_columns.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 with `[log] columns`: residuals, a probe in the lower-left vortex, one cell, the largest velocity and lowest pressure with their cells, the lid mass flow; bitwise the results of `cavity_re100.toml` | [gallery 28](../gallery/runs.md#log-columns), [recipe 38](../recipes/runs.md#log-columns) |
| `flow/cavity_re100_log_steps.toml` | generated hex 2D, 32 × 32 | the transient cavity at Re = 100 for 100 steps with `[log] detail = "steps"` (box and row of each logged step, no outer-iteration table); a `[[schedule]]` entry changes `log.columns` at step 50 (header printed again); the Limitations and Events blocks (D214) | [gallery 27](../gallery/runs.md#log), [recipe 38](../recipes/runs.md#log-columns) |
| `flow/cavity_re100_schedule.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 started with momentum under-relaxation 0.5; `[[schedule]]` raises it to 0.9 from iteration 100 and thins the log from iteration 200; every change is in the log and in `run.jsonl` | [gallery 29](../gallery/runs.md#schedule), [recipe 39](../recipes/runs.md#schedule) |
| `flow/cavity_re100_energy.toml` | generated hex 2D, 64 × 64 | the cavity at Re = 100 with the balance table every 100 iterations: the lid power goes into viscous and numerical dissipation, the numerical part positive (D90 part 4) | [gallery 19](../gallery/runs.md#health), [recipe 27](../recipes/runs.md#health) |

## The flow loop on a GPU

CPU = GPU bitwise up to v0.7; since MPI step W2 (#219) the CPU sums exactly and the GPU in blocks, so they agree to rounding until v1.1. Not run by CI.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cavity_poly3d_gpu.toml` | generated poly 3D, 970 225 polyhedra | the steady flow loop on the GPU: 30 SIMPLEC iterations of the cavity at Re = 100 on the CPU and the GPU, bitwise equal, time per cell and iteration (GPU; not run by CI) | [gallery 16](../gallery/benchmarks.md#gpu-flow), [recipe 24](../recipes/gpu.md#flow-gpu) |
| `flow/cylinder_re100_gpu.toml` | built-in cylinder mesh, refine 2: 98 304 hexahedra | the transient flow loop on the GPU: 200 adaptive BDF2 steps of the Re = 100 cylinder on the CPU and the GPU, bitwise equal (GPU; not run by CI) | [gallery 16](../gallery/benchmarks.md#gpu-flow), [recipe 24](../recipes/gpu.md#flow-gpu) |

The adjoint-ready kernels of stage 5 (D59) have no case file: their
checks are tests on generated meshes (gallery 17, recipe 25).

## Boundary conditions

Slip walls, axisymmetric wedges, inlet profiles from tables, backflow
treatment at pressure outlets and the convective outlet (D29, stage 5
part 5b), periodic pairs with the heat-flux jump of the scalars (D29,
D109 (1), stage 6 part 6p); the keys are on the [key page](flow-keys.md#boundary-conditions-d29).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cylinder_re100.toml` | built-in cylinder mesh, refine 0: 6 144 hexahedra | Schäfer & Turek 2D-2, periodic, Re = 100: C_D,max, C_L,max, St, Δp within 2.5 % (about 8 minutes; Linux, on main, nightly and with the label `full-ci`) | gallery 15, recipe 23 |
| `flow/channel_slip.toml` | generated mixed 2D, perturbed by 20 %, 432 cells | uniform flow between slip walls, exact (error ≤ 1e-9 required) | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/poiseuille_table.toml` | generated hex 2D, 512 cells | the inlet profile read from `poiseuille_profile.csv`: the run with the exact parabola | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/pipe_wedge.toml` | generated hex 2D revolved into a 0.5° wedge, perturbed, 256 to 4 096 cells (3 levels) | Hagen–Poiseuille flow on the axisymmetric wedge, observed order 1.8–2.2 required | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/cylinder_re100_short.toml` | built-in cylinder mesh with the outlet 5.5 D behind the cylinder, 4 928 hexahedra | the Re = 100 vortex street through a static pressure outlet with backflow treatment (Linux, on main, nightly and with the label `full-ci`) | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/cylinder_re100_convective.toml` | as `cylinder_re100_short.toml` | the same with a convective outlet (Linux, on main, nightly and with the label `full-ci`) | [gallery 21](../gallery/boundary-conditions.md#boundary-conditions), [recipe 29](../recipes/boundary-conditions.md#boundary-conditions) |
| `flow/poiseuille_periodic.toml` | generated hex 2D, perturbed by 15 %, 512 to 8 192 cells (3 levels) | plane Poiseuille flow on one periodic pair driven by its mass flow 0.1 kg/s: observed order 1.8–2.2 required (1.98, 1.94), the driving force converges to 12 μU/h² = 1.2 Pa/m (+0.09 % on the finest mesh) | [gallery 67](../gallery/boundary-conditions.md#periodic), [recipe 67](../recipes/boundary-conditions.md#periodic) |
| `flow/channel_heated_periodic.toml` | as `poiseuille_periodic.toml` | the laminar heated channel on one periodic module with the heat-flux jump: observed order of the temperature error 1.8–2.2 required, Nu = 8.2460 against 140/17 = 8.2353 on the finest mesh | [gallery 67](../gallery/boundary-conditions.md#periodic), [recipe 67](../recipes/boundary-conditions.md#periodic) |
| `flow/pipe_heated_periodic.toml` | generated hex 2D revolved into a 0.5° wedge, perturbed, 64 to 1 024 cells (3 levels) | the short periodic heated pipe (one diameter, D111): Hagen–Poiseuille flow at a fixed mass flow and the heat-flux jump; orders 1.8–2.2 (velocity) and 1.8–2.6 (temperature) required, Nu = 4.3691 against 48/11 = 4.3636 on the finest mesh | [gallery 67](../gallery/boundary-conditions.md#periodic), [recipe 67](../recipes/boundary-conditions.md#periodic) |

## Turbulence models

RANS turbulence models with their wall treatments (stage 6, D129–D139);
the keys are on the [key page](flow-keys.md#turbulence-d129d139).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/tmr_plate_sst.toml` | generated hex 2D on the TMR grid family (8 + 24 cells along x, plate split from the symmetry plane at x = 0, wall grading 38 000), 3 072 cells (`verification.levels = [0, 1, 2, 3]`: 768 to 49 152) | the TMR zero-pressure-gradient flat plate at Re = 5·10⁶/m with k-ω SST and `low-re` (SIMPLEC on level 1; the four-level sweep with `solver.algorithm = "coupled"`): C_f(0.97) against the Coles–Fernholz fit (+1.45 % on the finest mesh, −6.1 % from White), the wall shear along the plate, the profile at x = 0.97, the D166 features | [gallery 68](../gallery/turbulence.md#sst), [recipe 68](../recipes/turbulence.md#sst) |
| `flow/channel_sst.toml` | generated hex 2D, 2 × 96 cells, graded towards both walls (ratio 60), periodic in x | fully developed channel flow at Re_τ = 395 driven by the DNS bulk velocity: U_b⁺ = 17.43 against 17.55 (Moser, Kim & Mansour), −0.7 % | [gallery 68](../gallery/turbulence.md#sst), [recipe 68](../recipes/turbulence.md#sst) |
| `flow/pipe_turbulent_sst.toml` | generated hex 2D on a 0.5° wedge, 2 × 48 cells graded towards the wall (levels 96 and 192 in the study), periodic module of 0.1 m | fully developed turbulent pipe flow at Re_b = 5 300 with k-ω SST (`low-re`): u_τ/U_b = 0.0714 against 0.0683 of the DNS (El Khoury et al. 2013), +4.5 %; GCI_fine 0.2 % on three meshes (`examples/studies/pipe_gci.py`) | [gallery 193](../gallery/turbulence.md#pipe-turbulent), [recipe 193](../recipes/turbulence.md#pipe-turbulent) |
| `flow/channel_rotating_sst.toml` | generated hex 2D, 2 × 96 cells graded 30 towards both walls, periodic in x | the channel at Re_b = 5 800 in spanwise rotation (Ro = 0.15, `flow.frame_rotation`) with the curvature correction of SST: wall shear stress 1.63 times higher on the unstable side (1.00 without the correction) | [gallery 113](../gallery/turbulence.md#curvature-correction), [recipe 113](../recipes/turbulence.md#curvature-correction) |
| `flow/backstep_outlet_backflow.toml` | generated hex 2D, the step mesh of gallery 74 cut at x = 4 H, 4 600 cells | the backward-facing step with the outlet inside the recirculation, SST automatic, k and ω of the backflow at the outlet | [gallery 140](../gallery/turbulence.md#outlet-backflow), [recipe 140](../recipes/turbulence.md#outlet-backflow) |
| `flow/square_duct_qcr.toml` | generated hex 3D, 2 × 72 × 72 cells on a quarter duct graded 12 towards the two walls, periodic in x, symmetry on the two bisecting planes | the square duct at Re_b = 3 500 with SST and QCR2000: eight corner vortices, largest secondary velocity 1.56 % of U_b (none with linear SST) | [gallery 115](../gallery/turbulence.md#qcr), [recipe 115](../recipes/turbulence.md#qcr) |
| `flow/cavity_sst_operator_check.toml` | generated tet 2D, 64 × 32 cells (4 096 triangular prisms), vertices perturbed by 15 % | the operator exactness check (D236) of a shallow lid-driven cavity at Re = 10 000 with SST, the curvature correction and QCR2000: the QCR stress takes the quadratic least-squares velocity gradient and the curvature correction a quadratic fit of S_ij on all 4 096 cells (D236 steps 2 and 3), so no `MESH005` or `MESH006` remains; the check goes to `flow.json` and `operators.csv`; converges in about 1 030 iterations | [gallery 189](../gallery/meshes.md#operator-check), [recipe 189](../recipes/meshes.md#operator-check) |
| `flow/tmr_plate_keps.toml` | generated hex 2D, 16 + 48 cells along x, 48 across with wall grading 361 (first-cell y⁺ 30 at x = 0.97), 3 072 cells | the TMR flat plate with the standard k-ε model, `scalable` wall functions and the sustaining terms: C_f(0.97) = 2.829e-3, −1.6 % from White's correlation | [gallery 71](../gallery/turbulence.md#k-epsilon), [recipe 71](../recipes/turbulence.md#k-epsilon) |
| `flow/channel_keps.toml` | generated hex 2D, 2 × 12 uniform cells (first-cell y⁺ 32), periodic in x | channel flow at Re_τ = 395 with the standard k-ε and `scalable`: U_b⁺ = 18.31, +4.3 % from the DNS | [gallery 71](../gallery/turbulence.md#k-epsilon), [recipe 71](../recipes/turbulence.md#k-epsilon) |
| `flow/pipe_rough_keps.toml` | generated hex 2D revolved into a 0.5° wedge, 2 × 40 uniform cells, periodic in x | fully developed rough-pipe flow at Re = 10⁶, k_s/D = 10⁻³, Cebeci–Bradshaw shift: f = 0.02081, +4.4 % from Colebrook; base case of the Moody study | [gallery 71](../gallery/turbulence.md#k-epsilon), [recipe 71](../recipes/turbulence.md#k-epsilon) |
| `flow/channel_zeta_f.toml` | generated hex 2D, 2 × 96 cells graded 60 towards both walls (first-cell y⁺ 0.27), periodic in x | channel flow at Re_τ = 395 with k-ε-ζ-f and the `compound` treatment: U_b⁺ = 17.89, +1.9 % from the DNS | [gallery 73](../gallery/turbulence.md#zeta-f), [recipe 73](../recipes/turbulence.md#zeta-f) |
| `flow/tmr_plate_zeta_f.toml` | generated hex 2D, 16 + 48 cells along x, 48 across with wall grading 17 555 (first-cell y⁺ 1 at x = 0.97), 3 072 cells | the TMR flat plate with k-ε-ζ-f, `compound`, a turbulent start and the sustaining terms: C_f(0.97) = 3.023e-3, +5.1 % from White's correlation (finding) | [gallery 73](../gallery/turbulence.md#zeta-f), [recipe 73](../recipes/turbulence.md#zeta-f) |
| `flow/channel_zeta_f0.toml` | generated hex 2D, 2 × 48 cells graded 27.3 towards both walls (first-cell y⁺ 1), periodic in x | channel flow at Re_τ = 395 with k-ε-ζ-f₀ (Popovac 2023) and `compound-2023`: U_b⁺ = 17.79, +1.4 % from the DNS | [gallery 74](../gallery/turbulence.md#zeta-f0), [recipe 74](../recipes/turbulence.md#zeta-f0) |
| `flow/tmr_plate_zeta_f0.toml` | the mesh of `tmr_plate_zeta_f.toml` (first-cell y⁺ 1 at x = 0.97), 3 072 cells | the TMR flat plate with k-ε-ζ-f₀ from the TMR free stream: turns turbulent on its own, every field converged, C_f(0.97) = 2.802e-3, −2.5 % from White's correlation | [gallery 74](../gallery/turbulence.md#zeta-f0), [recipe 74](../recipes/turbulence.md#zeta-f0) |
| `flow/backstep_zeta_f0.toml` | generated hex 2D, 125 × 80 cells with the step cut out (`remove`), 8 200 cells, first-cell y⁺ about 1.5 | the backward-facing step of Driver & Seegmiller (Re_H = 36 000) with k-ε-ζ-f₀: reattachment x_r/H = 6.38, +2 % from the experiment | [gallery 74](../gallery/turbulence.md#zeta-f0), [recipe 74](../recipes/turbulence.md#zeta-f0) |
| `flow/backstep_zeta_f0_coupled.toml` | the mesh of `backstep_zeta_f0.toml`, 8 200 cells | the same step with the coupled algorithm and pseudo-time steps, at most 150 iterations: x_r within 1 % after about 40 iterations instead of 550 to 600, the same x_r/H = 6.385 | [gallery 78](../gallery/flow.md#coupled), [gallery 80](../gallery/flow.md#momentum-multigrid), [gallery 81](../gallery/flow.md#turbulence-multigrid), [recipe 78](../recipes/flow.md#coupled), [recipe 80](../recipes/flow.md#momentum-multigrid), [recipe 81](../recipes/flow.md#turbulence-multigrid) |
| `studies/backstep/keps.toml`, `zf0.toml`, `sst.toml` (in `examples/studies/`) | the mesh of `backstep_zeta_f0.toml` with each interval halved, 250 × 160 cells with the step cut out, 32 800 cells | the same step with k-ε and scalable wall functions, k-ε-ζ-f₀ and k-ω SST, with field output: x_r/H = 5.36, 6.35 and 6.80 against the experiment's 6.26 after 3 000 iterations; the velocity and k figures of the gallery | [gallery 74](../gallery/turbulence.md#zeta-f0), [recipe 74](../recipes/turbulence.md#zeta-f0) |
| `flow/channel_eb.toml` | the mesh of `channel_zeta_f.toml`, 2 × 96 cells graded 60 towards both walls (first-cell y⁺ 0.27), periodic in x | channel flow at Re_τ = 395 with the elliptic-blending k-ε (Billard & Laurence 2012): C_f at 99.0 % of the DNS (the paper: 98 %), U_b⁺ = 17.64, +0.5 % | [gallery 88](../gallery/turbulence.md#elliptic-blending), [recipe 88](../recipes/turbulence.md#elliptic-blending) |
| `flow/tmr_plate_eb.toml` | the mesh of `tmr_plate_zeta_f.toml` (first-cell y⁺ 1 at x = 0.97), 3 072 cells | the TMR flat plate with the elliptic-blending k-ε from a turbulent start: C_f(0.97) = 2.889e-3, +0.5 % from White's correlation | [gallery 88](../gallery/turbulence.md#elliptic-blending), [recipe 88](../recipes/turbulence.md#elliptic-blending) |
| `flow/backstep_eb.toml` | the mesh of `backstep_zeta_f0.toml`, 8 200 cells | the backward-facing step with the elliptic-blending k-ε (SIMPLEC): x_r/H = 7.59, +21 % from the experiment (finding: the variable C*_ε2; 6.24 with it constant) | [gallery 88](../gallery/turbulence.md#elliptic-blending), [recipe 88](../recipes/turbulence.md#elliptic-blending) |
| `flow/backstep_validation.toml` | generated hex 2D, 136 × 96 cells with the step cut out, 11 040 cells, graded towards the walls, the step corner and the reattachment region (first-cell y⁺ about 1.5) | the step of Driver & Seegmiller in the layout of NASA's TMR (walls from x = −110 H, outlet at 50 H, no-slip top wall at 9 H) with k-ε-ζ-f₀ and the coupled algorithm: x_r/H = 6.33 (+1.1 %) as a wall crossing, C_p within an RMS of 0.018 of the experiment; the base case of the six-model study `studies/backstep_validation.py` | [gallery 89](../gallery/turbulence.md#step-validation), [recipe 89](../recipes/turbulence.md#step-validation) |
| `flow/jet_pipe_profile.toml` | generated hex 2D revolved into a 0.5° wedge, 2 × 48 cells graded towards the wall (first-cell y⁺ about 1) | the developed pipe flow at Re = 23 750 with k-ε-ζ-f₀ (periodic module driven by the mass flow, Darcy f 0.0253 against 0.0255 of Blasius); its line `profile` becomes the inlet table of the impinging jet with `python3 -m foehncraft.profiles` | [gallery 90](../gallery/turbulence.md#impinging-jet), [recipe 90](../recipes/turbulence.md#impinging-jet) |
| `flow/impinging_jet.toml` | generated hex 2D revolved into a 0.5° wedge, 50 × 52 intervals refined once (10 400 cells), first cell at the plate 5·10⁻⁴ D (y⁺ below 1), clustered at the pipe lip | the round jet from a pipe on a heated plate at H/D = 2, Re = 23 750, k-ε-ζ-f₀ with the coupled algorithm: the inlet table of the developed pipe flow with k, ε and ζ columns, a passive temperature with a uniform wall flux, Nu₀ = 157 at the wall probe `stag` (156 at `r01`, r/D = 0.1) and a secondary maximum of Nu at r/D ≈ 2; the base case of the study `studies/impinging_jet.py` | [gallery 90](../gallery/turbulence.md#impinging-jet), [recipe 90](../recipes/turbulence.md#impinging-jet) |

## Derived fields

Fields derived from the solution with the output or afterwards from a
restart file (D72, stage 5 part 5c); the keys are on the
[key page](flow-keys.md#derived-fields-d72).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/couette_derived.toml` | generated poly 2D, perturbed by 15 %, 153 cells | every derived field of plane Couette flow, exact (vorticity −1/s, Q = λ₂ = 0, wall shear stress 0.1 Pa); `foehn post derive` writes the same file from `restart.h5` | [gallery 22](../gallery/post-processing.md#derived-fields), [recipe 30](../recipes/post-processing.md#derived-fields) |
| `flow/cylinder_re100_short.toml` | built-in cylinder mesh, 4 928 hexahedra | vorticity and Q criterion of the vortex street at t = 8 s, from the restart file with `foehn post derive` | [gallery 22](../gallery/post-processing.md#derived-fields), [recipe 30](../recipes/post-processing.md#derived-fields) |

## Statistics of transient runs

Time means, RMS, covariances and phase averages of a transient run with a
confidence interval of each mean (D71); the keys are on the
[key page](flow-keys.md#statistics-of-transient-runs-d71).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cylinder_re100_statistics.toml` | built-in cylinder mesh, 6 144 hexahedra | mean, RMS and ⟨u'v'⟩ of the Re = 100 vortex street from t = 5 to 12 s, a wake probe and C_D, C_L with their phase averages over the shedding period (Linux, on main, nightly and with the label `full-ci`) | [gallery 110](../gallery/post-processing.md#transient-statistics), [recipe 110](../recipes/post-processing.md#transient-statistics) |

## Expressions and materials

Values as expressions of x, y, z and t, constants and tables (D53), the
material database and regions (D102, D112), walls with a mass flux
(D124 W3), stage 5 part 6a; the keys are on the
[key page](flow-keys.md#expressions-materials-and-regions-d53-d102-d112).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/couette_viscosity.toml` | generated hex 2D, perturbed by 10 %, 128 to 2 048 cells (3 levels) | Couette flow with the viscosity μ0 (1 + 4 y) from an expression: u = U ln(1 + a y)/ln(1 + a), observed order 1.8–2.2 required | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |
| `flow/couette_transpiration.toml` | generated hex 2D, perturbed by 10 %, 128 to 2 048 cells (3 levels) | Couette flow with blowing and suction through the walls (`mass_flux`), Re = 4, observed order 1.8–2.2 required | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |
| `flow/poiseuille_water.toml` | generated hex 2D, perturbed by 10 %, 2 048 cells | water at 313.15 K from the material database, in a region of the zone `fluid`; Poiseuille flow in a 1 mm channel, error ≤ 1e-3 required | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 32](../recipes/properties.md#materials) |
| `flow/poiseuille_expression.toml` | as `poiseuille_hex.toml`, 3 levels | the inlet parabola and half of the driving pressure gradient as expressions: the errors of `poiseuille_hex.toml` | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |
| `flow/taylor_green_expression.toml` | generated hex 2D, 4 096 cells | the decaying vortex prescribed on the sides as expressions of x, y and t, error ≤ 2e-3 at t = 0.4 s required | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |
| `flow/channel_ramp.toml` | generated hex 2D, 512 cells | start-up of a channel with the inflow ramped by a table in t, BDF2 with an adaptive step; steady Poiseuille flow at t = 4 s (error ≤ 2e-3) | [gallery 23](../gallery/properties.md#expressions-materials), [recipe 31](../recipes/properties.md#expressions) |

## Passive scalars

Scalars carried by the flow (D51) with the generic scalar transport
equation of stage 6 (D116), the bounded BDF2 (D28 (1)), the MLP vertex
limiter (D237) and the exact wall distance (D48), stage 6 part 6a; the keys are on the
[key page](flow-keys.md#passive-scalars-d51-d116).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/channel_heated.toml` | generated hex 2D, perturbed by 15 %, 512 to 8 192 cells (3 levels) | laminar heated channel: a passive temperature with uniform wall flux in Poiseuille flow, Pe = 1000; observed order of the temperature error 1.8–2.2 required, Nu = 8.2398 against 140/17 = 8.2353 on the finest mesh | [gallery 63](../gallery/scalars.md#passive-scalars), [recipe 63](../recipes/scalars.md#passive-scalars) |
| `flow/scalar_front_bdf2.toml` | generated hex 2D, 40 cells | a front of a scalar bounded by [0, 1] in a uniform flow, upwind, BDF2 at Courant number 0.8: the bounded BDF2 keeps it inside the bounds (plain BDF2 overshoots to 1.00033) | [gallery 63](../gallery/scalars.md#passive-scalars), [recipe 63](../recipes/scalars.md#passive-scalars) |
| `flow/scalar_oblique_step_mlp.toml` | generated polyhedra 2D, 24 × 24, perturbed by 30 % | an oblique step of a passive scalar in a uniform flow, nearly no diffusion: the MLP vertex limiter with Nishikawa's R3 (D237) keeps c within [0, 1.0021] (Venkatakrishnan: [−0.0175, 1.0178]); the run fails outside [−0.005, 1.005] | [gallery 211](../gallery/scalars.md#mlp-limiter), [recipe 211](../recipes/scalars.md#mlp-limiter) |

## Energy equation

The temperature of the fluid with c_p(T), k(T) and μ(T), viscous heating
and the thermal wall conditions (D150), stage 7 part 7a1; the keys are on
the [key page](../keys.md) and in [recipe 154](../recipes/energy.md#energy-fluid).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/cavity_heated.toml` | generated hex 2D, 4 096 cells | the Re = 100 cavity with the lid at 310 K and the bottom at 300 K, Pr = 0.71: forced convection, Nu = 2.05, heat balance closed to 10⁻¹¹ | [gallery 154](../gallery/energy.md#energy-fluid), [recipe 154](../recipes/energy.md#energy-fluid), [gallery 238](../gallery/energy.md#ledger), [recipe 238](../recipes/energy.md#ledger) |
| `flow/conduction_mms.toml` | generated hex 2D, perturbed by 15 %, 256 to 4 096 cells (3 levels) | check V7.1: conduction with k(T) = 1 + 0.01 (T − 300), manufactured solution; observed order of the temperature error 1.8–2.2 required | [gallery 154](../gallery/energy.md#energy-fluid), [recipe 154](../recipes/energy.md#energy-fluid) |
| `flow/couette_viscous_heating.toml` | generated hex 2D, 128 to 2 048 cells (3 levels), periodic | check V7.2: Couette flow heated by its own friction, Br = 1, the exact parabola; order 1.8–2.2 required, each wall takes half of the friction power | [gallery 154](../gallery/energy.md#energy-fluid), [recipe 154](../recipes/energy.md#energy-fluid) |
| `flow/cht_composite_wall.toml` | generated hex 2D, 120 cells, zones `a` and `b` | check V7.4: two solids and a fluid layer at rest, conductivity ratios 1 : 100 : 10⁴; the series heat flux 495.0005 W/m² at both walls and both interfaces, and the interface temperatures, to round-off | [gallery 156](../gallery/energy.md#cht-solids), [recipe 156](../recipes/energy.md#cht-solids) |
| `flow/cht_heated_wall.toml` | generated hex 2D, 80 × 30 cells, zone `wall` | conjugate heat transfer: a solid wall heated with 1000 W/m² from below under a laminar channel flow (Re = 100, Pr = 0.71, k_s/k_f = 10); all 160 W of the heater cross `interface-wall`, balance closed to 10⁻¹¹ | [gallery 156](../gallery/energy.md#cht-solids), [recipe 156](../recipes/energy.md#cht-solids) |
| `flow/cht_turbulent_channel.toml` | generated hex 2D, 40 × 33 cells, zone `wall` | check V7.36: a turbulent channel (Re = 2·10⁵, Pr = 0.71, k-ε) over a solid wall at 350 K with t/k_s = 0.5 m²K/W; the wall function in series with the solid gives 107.692 W through `interface-wall`, the fluid-only path 107.681 W | [gallery 158](../gallery/energy.md#cht-turbulent), [recipe 158](../recipes/energy.md#cht-turbulent) |
| `flow/cht_heat_sink_turbulent.toml` | generated hex 2D, 300 × 50 cells, zones `base`, `fin_1` to `fin_3` | turbulent conjugate heat transfer: an aluminium heat sink with three fins in a 5 m/s air channel (k-ω SST), 2000 W/m² from a chip below; all 2 W cross `interface-sink`, the metal stays within 0.7 K, balance closed to 4·10⁻⁷ | [gallery 158](../gallery/energy.md#cht-turbulent), [recipe 158](../recipes/energy.md#cht-turbulent) |
| `flow/cht_heated_wall_transient.toml` | generated hex 2D, 80 × 30 cells, zone `wall` | check V7.37: the heated wall of `cht_heated_wall` (lighter, ρ = 100) with the heater switched on at t = 0, 600 BDF2 steps of 5 s; the heat balance of `heat.csv` closes in every step to 2·10⁻¹¹, the end state is the steady one within 1.5·10⁻⁵ of the temperature range; probes with the temperature in the wall | [gallery 159](../gallery/energy.md#cht-transient), [recipe 159](../recipes/energy.md#cht-transient) |
| `flow/cht_heated_wall_accelerated.toml` | generated hex 2D, 80 × 30 cells, zone `wall` | check V7.41: the switched-on heated wall with `[energy.acceleration] inertia_scale = "auto"`: G = 8.96 W/K and τ = 286 s measured after the first window, f = 0.175, settled and restored at 450 s; at 600 s the heated face is within 0.0025 K of the true wall after 3000 s; `inertia_scaling` in `heat.csv` | [gallery 180](../gallery/energy.md#cht-acceleration), [recipe 180](../recipes/energy.md#cht-acceleration) |
| `flow/cht_heated_wall_periodic_accelerated.toml` | generated hex 2D, 80 × 30 cells, zone `wall` | the accelerated wall with a periodic inlet temperature (300 ± 50 K, 50 s): settled at 700 s, ρ c_p restored over three periods; the cycle mean at 1500 s within 0.0025 K of the true periodic state (a jump: 0.24 K too cold); ACC003 note | [gallery 180](../gallery/energy.md#cht-acceleration), [recipe 180](../recipes/energy.md#cht-acceleration) |
| `flow/cht_heated_wall_frozen.toml` | generated hex 2D, 80 × 30 cells, zone `wall` | check V7.43: the switched-on heated wall with a viscosity that halves over 17 K, warmed up with `[energy.frozen_flow]` (flow step 0.5 s, march step 10 s): held at 55 s, 7 refreshes, 557 steps in 17 s against 6000 steps in 155 s fully coupled (`enabled = false`); every probe within 0.058 K of the fully coupled history over a 21 K rise, balance closed to 3·10⁻¹² | [gallery 208](../gallery/energy.md#cht-frozen-flow), [recipe 208](../recipes/energy.md#cht-frozen-flow) |
| `flow/cht_pulsating_channel_cycle.toml` | generated hex 2D, 80 × 20 cells, zone `wall` | checks V7.39, V7.40: a pulsating laminar channel (inlet velocity ±50 %, temperature 300 ± 50 K, period 50 s) over a wall heated with 160 W, coupled by `[energy.cycle]`: converged after 4 updates and 10 cycles (500 s simulated, 5.5 s); 160.0002 W into the fluid, wall probes within 0.14 K of the fully coupled periodic state (`enabled = false`, 4000 s); ledgers to 10⁻¹³ | [gallery 190](../gallery/energy.md#cht-cycle), [recipe 190](../recipes/energy.md#cht-cycle) |
| `flow/cht_pulsating_channel_skin.toml` | generated hex 2D, 80 × 80 cells, zone `wall` | the pulsating channel over a wall with a thin skin (δ = 0.03 m in 0.2 m), coupled by `[energy.cycle]` with `skin = "effusivity"` and `resolved_cycles = 3`: 10 updates, ends at 1250 s (96 s) within 0.016 K of the fully coupled periodic state and within 0.1 % on the swing; the fully coupled run (`enabled = false`, 10 000 s) is within 0.1 K from 4400 to 5550 s | [gallery 202](../gallery/energy.md#cht-cycle-skin), [recipe 202](../recipes/energy.md#cht-cycle-skin) |
| `flow/cht_pulsating_channel_warmup.toml` | generated hex 2D, 80 × 20 cells, zone `wall` | check V7.44: the pulsating channel with a period of 5 s warmed up by `[energy.cycle] solid = "transient"` (`march_cycles = 20`, `march_steps = 5`, skin swing): 15 updates to 1500 s, 780 steps in 18 s; the wall at the end of every march within 0.21 K of the fully coupled cycle means over a 35 K rise, 0.01 K at the end; the fully coupled run (`enabled = false`) takes 6000 steps | [gallery 208](../gallery/energy.md#cht-frozen-flow), [recipe 208](../recipes/energy.md#cht-frozen-flow) |
| `flow/cht_counterflow_two_fluids.toml` | generated hex 2D, 200 × 18 cells, zones `water`, `plate`, `oil` | check V7.33 (planar form, D194): water and a light oil in counterflow across a steel plate, C_r = 1, two fluid parts each with its own outlet: ε = 0.7756 against ε–NTU 0.7751 with Nu = 5.385, the heat of the water, the plate and the oil agree, ledger to 3·10⁻⁵ (2670 iterations, 14 s) | [gallery 206](../gallery/energy.md#two-fluids), [recipe 206](../recipes/energy.md#two-fluids) |
| `flow/cht_two_fluids_absolute.toml` | generated hex 2D, 200 × 18 cells, zones `water`, `plate`, `oil` | absolute pressure values (D239, case format 0.7): the same exchanger with the water outlet at 1.0·10⁵ Pa and the oil outlet at 1.5·10⁵ Pa; each part takes its own outlet as its operating pressure (PRS005 in the log, `pressure_values` in `flow.json`), the `pressure` field and `boundaries.csv` are absolute, the solved fields are those of the gauge case bit for bit | [gallery 220](../gallery/boundary-conditions.md#absolute-pressure), [recipe 220](../recipes/boundary-conditions.md#absolute-pressure) |
| `flow/cht_two_fluids_start.toml` | generated hex 2D, 200 × 18 cells, zones `water`, `plate`, `oil` | start values per region and per fluid part (D238): the absolute exchanger with the oil driven by a total pressure 46 Pa above its outlet; water and oil start at their inlet temperatures (`"inlet"`), the plate at 320 K, the oil pressure on a Laplace solve between its inlet and outlet (REG010 in the log, `start_values` in `flow.json`); 2251 instead of 2616 iterations from a uniform start | [gallery 221](../gallery/energy.md#start-values), [recipe 221](../recipes/energy.md#start-values) |
| `flow/cht_contact_transient.toml` | generated hex 2D, 400 × 1 cells, zone `bar` | check V7.6: a solid at 400 K (k = 10) and a fluid at rest at 300 K (k = 1) brought into contact; the interface stays at T_i = 375.975 K, the erf profiles of `[energy] exact` to L2 0.0048 K at 500 s, second order in space and time | [gallery 159](../gallery/energy.md#cht-transient), [recipe 159](../recipes/energy.md#cht-transient) |
| `flow/cht_orthotropic_patch.toml` | generated hex 2D, 20 × 20 cells perturbed by 25 %, zones `laminate`, `block` | stage 7d2 patch test: an orthotropic laminate (k₁ = 20, k₂ = 2, axes turned by 30°), a fluid at rest and a block behind R″_c = 0.01 m²K/W (`[interfaces]`); piecewise-linear T exact to 3·10⁻¹¹ K, 200 W/m² through every wall and interface | [gallery 160](../gallery/energy.md#cht-contact), [recipe 160](../recipes/energy.md#cht-contact) |
| `flow/cht_circuit_board.toml` | generated hex 2D, 120 × 30 cells, zones `board`, `chip` | a 1 W chip (`[[energy.sources]]`) on a thermal pad (R″_c = 2·10⁻⁴ m²K/W) on an orthotropic board (k = 25 along, 0.35 through) in a 1 m/s air channel; chip 363.7 K (423.7 K with a plain laminate), balance closed to 7·10⁻¹³ | [gallery 160](../gallery/energy.md#cht-contact), [recipe 160](../recipes/energy.md#cht-contact) |

## Gravity and buoyancy

`[gravity]` and the Boussinesq buoyancy of stage 7b1 (D153) with the
balanced force of D154; the solved pressure is the reduced pressure
p′ = p − ρ_ref g·(x − x_ref), the output has the static pressure p as
`static_pressure` as well; the keys are on the
[key page](flow-keys.md#boundary-conditions-d29).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/buoyancy_enclosure_plume.toml` | generated hex 2D, 100 × 60 cells | air in a closed enclosure 10 cm × 6 cm over a floor strip up to 5 K warmer, Ra ≈ 1.1·10⁵: a laminar plume leaning to the nearer wall and one large roll, 1.9 cm/s; SIMPLEC with the buoyancy step, 2009 iterations, 30 s | [gallery 186](../gallery/energy.md#buoyancy), [recipe 186](../recipes/energy.md#buoyancy) |
| `flow/buoyancy_box_at_rest.toml` | generated hex 2D, 16 × 16 cells perturbed by 20 % | check V7.7 (a): a closed box with a stable stratification of 14 K along a gravity 37° from the y axis stays at rest, largest speed 1·10⁻¹⁵ m/s after 1000 iterations (tests/buoyancy.rs: also tet, poly, mixed 3D and the coupled algorithm) | [gallery 186](../gallery/energy.md#buoyancy), [recipe 186](../recipes/energy.md#buoyancy) |
| `flow/buoyancy_slot.toml` | generated hex 2D, 8 × 2 to 32 × 8 cells perturbed by 15 % (3 levels), periodic in y | the buoyant slot between a hot and a cold vertical wall against its exact cubic profile v = (g β ΔT D²/(48 ν))(ξ³ − ξ); observed orders 2.14 and 2.13 | [gallery 186](../gallery/energy.md#buoyancy), [recipe 186](../recipes/energy.md#buoyancy) |
| `flow/mixed_convection_channel.toml` | generated hex 2D, 8 × 2 to 32 × 8 cells perturbed by 15 % (3 levels), periodic in y with `mass_flow` | V7.9: fully developed mixed convection between a hot and a cold vertical wall (Aung & Worku 1986), Gr/Re = 200 with flow reversal at the cold wall, against the exact profile; observed orders 2.14 and 2.14 (2.16 and 2.08 with the net flow down) | [gallery 272](../gallery/energy.md#natural-convection), [recipe 272](../recipes/energy.md#natural-convection) |
| `flow/natural_convection_cavity_ra1e5.toml` | generated hex 2D, 32² to 128² cells clustered to the walls (3 levels) | V7.15: differentially heated square cavity, air, Ra 10⁵ (10³ to 10⁶ by `gravity.vector`), coupled: mean Nu 4.5213, order 2.13, GCI 0.019 %; Ra 10³, 10⁴, 10⁶: 1.1172, 2.2444, 8.8255 (GCI 0.004 %, 0.006 %, 0.124 %); de Vahl Davis (1983) and Le Quéré (1991) as context only (D201) | [gallery 272](../gallery/energy.md#natural-convection), [recipe 272](../recipes/energy.md#natural-convection) |
| `flow/natural_convection_plate.toml` | generated hex 2D, 32 × 40 to 128 × 160 cells graded to the plate and the leading edge (3 levels) | V7.12: isothermal vertical plate, Pr 0.71, Gr_L = 10⁸, open ambient (total pressure 0): extrapolated Nu_x at 0.25, 0.5, 0.75 m 0.37 %, 0.31 %, 0.26 % below the Ostrach similarity solution (`tools/verification/ostrach.py`), GCI ≤ 0.03 % | [gallery 272](../gallery/energy.md#natural-convection), [recipe 272](../recipes/energy.md#natural-convection) |
| `flow/buoyancy_channel_reduced.toml` | generated hex 2D, 10 × 100 cells | a liquid rising at 0.05 m/s in a channel 2 m tall with gravity, outlet p′ = 0 (the default kind): p′ carries the friction (6.5 Pa), p adds −ρ g y | [gallery 186](../gallery/energy.md#buoyancy), [recipe 186](../recipes/energy.md#buoyancy) |
| `flow/buoyancy_channel_static.toml` | as above | the same channel with the outlet as a static pressure, −19 620 Pa (`pressure_kind = "static"`): the same flow to 10⁻¹¹ | [gallery 186](../gallery/energy.md#buoyancy), [recipe 186](../recipes/energy.md#buoyancy) |
| `flow/buoyancy_open_box.toml` | generated hex 2D, 20 × 20 cells | a box at rest open over its full height on one side: the reduced value 0 keeps it at rest; with `--set 'boundary.xmax.pressure_kind="static"'` the same 0 Pa as a uniform static value drives a false flow of 3.5 m/s (BUO006) | [gallery 186](../gallery/energy.md#buoyancy), [recipe 186](../recipes/energy.md#buoyancy) |

## Variable density

The low-Mach variable density ρ(T, p₀) of stage 7c (D151): `[energy]
density = "variable"` or `[buoyancy] model = "variable-density"`; p₀ the
operating pressure of an open domain or, with `thermodynamic_pressure =
"closed-volume"`, that of the fixed mass of a closed one. The keys are
on the [key page](../keys.md) and in [recipe
195](../recipes/energy.md#variable-density).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/vd_cavity_large_dt.toml` | generated hex 2D, 32² to 128² cells clustered to the walls (3 levels) | air in a closed square cavity of 3.2 cm with walls at 960 K and 240 K (ε = 0.6, Ra ≈ 10⁵), the coupled algorithm: the heat flow through the hot wall 0.14858 W (per mm of depth), observed order 2.13, GCI 0.020 %; p₀ = 94 975 Pa, order 2.2, GCI 0.003 % | [gallery 195](../gallery/energy.md#variable-density), [recipe 195](../recipes/energy.md#variable-density) |
| `flow/vd_heated_duct.toml` | generated hex 2D, 100 × 4 cells | check V7c.1: gas heated by 100 K in a duct with slip walls; the mass flow is the same at both ends, ṁ c_p ΔT equals the source, the inlet pressure 0.392245 Pa is G²(1/ρ_out − 1/ρ_in) = 0.3922 Pa | [gallery 195](../gallery/energy.md#variable-density), [recipe 195](../recipes/energy.md#variable-density) |
| `flow/vd_closed_box_heating.toml` | generated hex 2D, 10 × 10 cells | check V7c.2: a closed box of gas at rest heated uniformly with 10⁴ W/m³ warms with c_v; T = 359.190 K and p₀ = 121 316 Pa after 5 s, both linear in t and exact to 10⁻¹⁰ with any Δt | [recipe 195](../recipes/energy.md#variable-density) |
| `flow/vd_closed_box_side_heated.toml` | generated hex 2D, 32 × 32 cells | check V7c.3: a closed box of air heated through one wall at 400 K for 2 s: p₀ rises to 117 630 Pa, the gas far from the wall is compressed to 321 K, the mass stays to round-off (10⁻¹⁵) | [gallery 195](../gallery/energy.md#variable-density), [recipe 195](../recipes/energy.md#variable-density) |
| `flow/vd_mms_steady.toml` | generated hex 2D, 16² to 64² cells (3 levels) | check V7c.4: manufactured solution with ρ = 1/T, a divergence-free mass flux, the full stress and a closed volume (`tools/verification/vd_mms.py`); observed orders 2.22 and 1.92 (velocity), 1.86 and 1.90 (temperature), 1.8–2.2 required | [recipe 195](../recipes/energy.md#variable-density) |

## Lagrangian sprays

`[spray]` with an included injector file (L12, SP1, D246): droplets as
parcels in a transient flow, one way or, with `[spray.coupling] mode =
"two-way"` (SP2, D247 to D252), with the drag as a momentum source of the
gas and, with `[energy]`, the heat as a source of the gas energy;
`[spray] dispersion = "eddy-lifetime"` (SP2c, D243) lets the turbulence
spread the parcels. In a steady case `[spray.steady]` (SP2s, D260) runs the
spray as streams through the steady gas, with relaxed sources. The keys are on the [key page](flow-keys.md); the run writes `spray.csv` and
`spray/parcels_<step>.vtkhdf`. `[[datum]]` tables (SP1b, D262) report the
liquid and the gas through planes, surfaces and boundaries and the droplets
in volumes, in `datums.csv` and `spray_hist_<datum>.csv`.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/spray_duct.toml` | generated hex 2D, 80 × 20 cells | 10 g/s of 150 µm water droplets sprayed from the upper wall into air at 4 m/s; about 70 % stick on the lower wall, the mass book closes; the injector is in `flow/spray_duct_injectors.toml` | [gallery 242](../gallery/sprays.md#spray-duct), [recipe 242](../recipes/sprays.md#spray-duct) |
| `flow/spray_datums.toml` | generated hex 2D, 80 × 20 cells | the duct of `spray_duct` with Rosin-Rammler sizes (X = 120 µm, n = 2.5; `flow/spray_datums_injectors.toml`) and seven datums: a cut plane, a plane on the faces, a disc, an STL surface (`flow/spray_datums_section.stl`), the outlet, the floor and a box; the outlet datum equals the escaped liquid to the bit, Dv50 91 µm at the outlet and 170 µm on the floor | [gallery 265](../gallery/sprays.md#spray-datums), [recipe 265](../recipes/sprays.md#spray-datums) |
| `flow/spray_two_way.toml` | generated hex 2D, 60 × 75 cells | two-way coupling: 50 g/s of water (50 to 400 µm, `flow/spray_two_way_injectors.toml`) sprayed down from the ceiling of a room of still air pulls the air into a jet of up to 6 m/s; the spray hands the air 0.77 N, 2.96 passes per step on average, the momentum pair cancels to the bit | [gallery 280](../gallery/sprays.md#spray-two-way), [recipe 280](../recipes/sprays.md#spray-two-way) |
| `flow/spray_cooling.toml` | generated hex 2D, 100 × 20 cells | two-way heat: 10 g/s of water at 20 °C (20 to 200 µm, `flow/spray_cooling_injectors.toml`) sprayed into air at 300 °C and 5 m/s; the droplets take about 6.6 kW from the air, which leaves about 50 K cooler; the enthalpy pair cancels to the bit, 3.7 passes per step with `max_passes = 4` | [gallery 283](../gallery/sprays.md#spray-cooling), [recipe 283](../recipes/sprays.md#spray-cooling) |
| `flow/spray_dispersion.toml` | generated hex 2D, 150 × 30 cells | eddy-lifetime dispersion: water droplets of 20 to 300 µm (`flow/spray_dispersion_injectors.toml`) shot along the axis of a duct with 10 % turbulence spread by the eddies alone (no gravity, no cone); after 1 m the small droplets spread 5.4 cm, the large ones 2.9 cm; no parcel lost | [gallery 285](../gallery/sprays.md#spray-dispersion), [recipe 285](../recipes/sprays.md#spray-dispersion) |
| `flow/spray_cooling_steady.toml` | generated hex 2D, 100 × 20 cells | the duct of `spray_cooling` as a steady state in the steady spray mode (`[spray.steady]`, D260): 400 water streams through the frozen gas, sources relaxed every 10 iterations; converges in 526 iterations and 58 passes, 6.39 kW from the air (the transient run: 6.58 kW), source residual 1e-4, all liquid leaves through the outlet | [gallery 292](../gallery/sprays.md#spray-steady), [recipe 292](../recipes/sprays.md#spray-steady) |

## Porous regions

`kind = "porous"` with `[regions.<name>.porous]` (D270, porous plan PM1):
a Darcy–Forchheimer resistance in the axes of the region, the velocity is
the superficial one; with `[energy]` and `solid` gas and matrix have
temperatures of their own (D271, PM2). The keys are on the [key page](flow-keys.md).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/porous_darcy_block.toml` | generated hex 2D, 100 × 4 cells | check V1: a porous block (K⁻¹ = 10⁶ 1/m², F = 10 1/m) in a duct with slip walls; the drop 4.780 Pa against the exact 4.8 Pa (0.41 %), second order in the cell size | [gallery 273](../gallery/porous.md#porous-flow), [recipe 273](../recipes/porous.md#porous-flow) |
| `flow/porous_honeycomb_channel.toml` | generated hex 2D, 240 × 80 cells | a ceramic honeycomb block (6.9 mm channels, ε = 0.256) in an air channel at 10 cm/s, channels along the flow; with `blocking = 1` the same resistance in every direction | [gallery 273](../gallery/porous.md#porous-flow), [recipe 273](../recipes/porous.md#porous-flow) |
| `flow/porous_schumann.toml` | generated hex 2D, 100 × 1 cells | check V4: air at 400 K into a ceramic bed at 300 K, gas and matrix against the exact solution of Schumann (1929) at 1000, 2000 and 3000 s; L1 error 5.2·10⁻⁵ of the step, second order in the cell size; heat balance to 10⁻¹³ | [gallery 275](../gallery/porous.md#porous-two-temperatures), [recipe 275](../recipes/porous.md#porous-two-temperatures) |
| `flow/porous_storage_charge.toml` | generated hex 2D, 40 × 100 cells | charging a honeycomb store: air at 573 K from the top right through a stack of ceramic bricks (6.9 mm channels, ε = 0.256) at 293 K; the thermocline half way after 6000 s, 200 steps in about 70 s | [gallery 275](../gallery/porous.md#porous-two-temperatures), [recipe 275](../recipes/porous.md#porous-two-temperatures) |

## Immersed bodies

`[immersed]` with `[[immersed.body]]` tables (L13 stage 1, D226): cylinders,
spheres and boxes put into any mesh by Brinkman penalisation, laminar and
isothermal only. The keys are on the [key page](../keys.md) (`immersed.*`);
`result.json` has the force on every body, `derived = ["solid-fraction"]`
writes χ.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/immersed_cylinder_re20.toml` | generated hex 2D, 440 × 82 cells (h = D/20), no cylinder in the mesh | Schäfer & Turek 2D-1 with the cylinder as a penalised body: C_D 5.581 (in range), Δp 0.9 % above, C_L 35 % above the range (smeared wall); fails beyond 40 % of the range midpoint | [gallery 270](../gallery/immersed.md#penalised-cylinder), [recipe 270](../recipes/immersed.md#penalised-body) |

## Deforming meshes

`[mesh_motion] displacement`: the nodes of the mesh move by an expression of
their place in the mesh file and of time (D91, L22, stage MS0); every face
gets its swept volume, walls move with their nodes, other boundaries may only
slide in their plane. The keys are on the [key page](../keys.md)
(`mesh_motion.*`); `history.csv` gets a column `volume`, `result.json` the
block `transient.mesh_motion`.

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/motion_free_stream.toml` | generated hex 2D, 20 × 10 cells perturbed by 20 %, inner nodes swinging by up to 3 cm | check MS0.1: a uniform flow of 1 m/s stays uniform on the moving mesh to 1.3·10⁻¹⁵ m/s, the space conservation law closes to 10⁻¹⁵ (tests/mesh_motion.rs: also tet and poly meshes, BDF1 and Crank–Nicolson) | [gallery 296](../gallery/moving-solids.md#deforming-mesh), [recipe 296](../recipes/moving-solids.md#deforming-mesh) |
| `flow/motion_piston_adiabatic.toml` | generated hex 2D, 10 × 20 cells, compressed evenly to a quarter | check MS0.2: a piston compresses a closed gas volume with adiabatic walls by a crank (ratio 4, 1 s): p₀ = 705.03 kPa at the top dead centre (isentrope 705.50 kPa), 101.323 kPa after the period; error second order in Δt, mass to 10⁻¹⁵ | [gallery 296](../gallery/moving-solids.md#deforming-mesh), [recipe 296](../recipes/moving-solids.md#deforming-mesh) |
| `flow/motion_piston_cooled_walls.toml` | generated hex 2D, 40 × 80 cells | the same piston with all walls at 300 K: thermal boundary layers, p₀ = 572 kPa at the top dead centre and 90.4 kPa after the period (an open p–V loop), heat book closed every step; 30 s | [gallery 296](../gallery/moving-solids.md#deforming-mesh) |

## Species

Mixtures of ideal gases (steam furnace stage SF1): `[species]` with the
GRI-Mech 3.0 species, N − 1 mass fractions transported and the closure
species by 1 − Σ; with the variable density ρ = p₀ M/(R T) of the mixture
and its c_p, λ and μ. The keys are on the [key page](flow-keys.md#species).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `flow/species_binary_diffusion.toml` | generated hex 2D, 100 × 2 cells | CH₄ diffusing in N₂ at rest: the exact decay of a cosine mode, BDF2 to t = 50 s; L2 error 9.5·10⁻⁶, species and element accounts closed to round-off | [gallery 278](../gallery/species.md#species-transport), [recipe 278](../recipes/species.md#species-transport) |
| `flow/species_coaxial_jet.toml` | generated hex 2D wedge, 90 × 40 cells | a laminar methane jet (r = 5 mm, 1 m/s) in a coflow of air at 0.2 m/s, isothermal, density from the mixture's molar mass (0.65 to 1.17 kg/m³), mixture properties from kinetic theory; 154 iterations, element accounts closed to 6.4·10⁻⁹ | [gallery 278](../gallery/species.md#species-transport), [recipe 278](../recipes/species.md#species-transport) |
