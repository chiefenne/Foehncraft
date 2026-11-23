# Changelog

## Unreleased

Nothing yet.

## 1.0.0 (2026-11-23) – milestone 1: a general purpose incompressible solver with heat transfer

**Foehncraft 1.0 in one paragraph.** A finite-volume solver for
incompressible and low-Mach flow on arbitrary polyhedral meshes in 2D and
3D, written in Rust, with one kernel source for the CPU and NVIDIA GPUs.
Steady and transient laminar and RANS flow (k-ω SST with curvature
correction and QCR, standard and realizable k-ε, k-ε-ζ-f, k-ε-ζ-f₀,
elliptic-blending k-ε, with every common wall treatment); SIMPLEC, PISO
and a coupled pressure–velocity algorithm with multigrid; the energy
equation with temperature-dependent properties, conjugate heat transfer
with any number of solids and fluids, gravity, Boussinesq buoyancy and the
low-Mach variable density; rotating machinery with the moving reference
frame and sliding meshes; readers for Gmsh, SU2, CGNS, OpenFOAM, Fluent
`.msh` and VTK meshes; restart, schedules, probes, statistics, ledgers and
divergence reports; a Python package, a typed tool registry and an MCP
server for scripts and assistants. Every feature has a verification grade
with its evidence on the [verification status](docs/verification/status.md)
page, and every run lists its limitations.

**Highlights since 0.7.**

- **Internal flows of D45 with grid convergence.** Turbulent pipe flow at
  Re_b 5 300 (u_τ/U_b +4.5 % from the DNS, GCI 0.2 %), the asymmetric
  plane diffuser of Buice & Eaton (reattachment 28.9 H, GCI 0.36 %), the
  U-bend of Monson et al. (separation at 139.8°, GCI 0.19 %), both also
  with sharp edges, and a laminar T-junction mixing a passive scalar on
  polyhedra; gallery 193, 196, 199, 204.
- **Several fluids in one conjugate case (D194).** Any number of fluid
  regions with their own liquids or gases beside any number of solids;
  each fluid part with its own pressure level. A shell and tube heat
  exchanger runs hot air, water and steel in one case on 0.48 M
  tetrahedra (application 203, a demonstration).
- **Gravity, buoyancy and variable density (stages 7b1 and 7c).**
  Boussinesq buoyancy with a balanced force (a stratified box at rest
  stays at rest to 10⁻¹⁵ m/s on any cell type) and the low-Mach density
  ρ(T, p₀) for open and closed domains, steady and transient, turbulent
  and conjugate.
- **Absolute pressures and start values (D239, D238; case format 0.7).**
  Pressure values in a case file are absolute by default (gauge optional,
  older files read as gauge); every region and fluid part starts from its
  own values and its own pressure boundaries. The case format is stable
  from 1.0 on.
- **Slow thermal transients made fast (D204, D118).** The cycle-averaged
  loop for periodic flows over slow solids, with a skin swing and fully
  coupled closing cycles, and the frozen-flow march: warm-ups in a tenth
  of the steps within a fraction of a kelvin of the fully coupled run.
- **Robust numerics on real meshes.** An operator exactness check of every
  mesh (D236) with quadratic fits that make the curvature correction and
  the QCR converge on tetrahedra; the vertex-based MLP limiter for scalars
  and temperature (D237); exact, order-independent sums, dot products and
  norms on the CPU, bitwise the same for any thread count (MPI step W2,
  with the canonical local mesh of W1 as groundwork for ranks).

**Known limits.** Runs use one process (threads, or one GPU): multi-rank
MPI is v1.1 (W1 and W2 are in, D235). Since the exact CPU sums, GPU
results agree with the CPU to rounding, not bitwise (GPU exact sums in
v1.1). No natural-convection benchmark validation yet (7b2 is v1.1), no
turbulent buoyancy closures (7e). The T-junction's branch-inlet pressure
does not converge with the mesh (GCI 29 %) and the D45 cases still lack
measurements, Re 10⁶ and a converged sharp-edged fine U-bend mesh; the
heat exchanger is a demonstration on one mesh. Several fluids do not yet
combine with variable density, gravity, rotating zones or passive scalars
(REG006); periodic meshes with solids are refused (ENE006). The GPU paths
of the stage 7 physics were not run on a device for this release, and the
generated verification manual was last run on 0.6.0.

**Changes.**

- **Gravity and the Boussinesq buoyancy (stage 7b1).** `[gravity]` with
  `vector` and `reference_point`, and `[buoyancy]` with
  `model = "boussinesq"`, `expansion_coefficient` (default from the
  fluid's density at T_ref), `reference_temperature` (default the
  `[fluid] temperature`) and `reference_density`. The solved pressure is
  the reduced pressure p′ = p − ρ_ref g·(x − x_ref); the output has the
  static pressure p as `static_pressure` as well. The values of the
  pressure boundaries are p′ by default; `pressure_kind = "static"` reads
  them as static pressures and converts them per face. Works with SIMPLEC
  and the coupled algorithm (the kernels are shared by CPU and GPU), with
  the energy equation and with solid regions (the fluid temperature drives
  the force).
- **Balanced force (D154).** The buoyancy enters every face as the
  pressure gradient does: the same face differences, the same
  least-squares operator in the cells, the matching Rhie–Chow term and
  boundary pressure. A closed box with a stable stratification along an
  oblique gravity stays at rest to 10⁻¹⁵ m/s on perturbed hex, tet, poly
  and mixed meshes, with SIMPLEC and the coupled algorithm (V7.7 (a)).
- **Buoyancy pseudo-time step.** Steady runs with buoyancy take a
  pseudo-time step 1/N (N the buoyancy frequency, scaled by
  `[solver.pseudo_time] factor`) in SIMPLEC and SIMPLE, and the coupled
  algorithm's automatic step takes the buoyancy velocity as a lower bound.
  Without it a stable stratification on a coarse mesh never settled.
- **Kinetic-energy ledger** counts the work of the buoyancy force.
- **Codes** BUO001 (no gravity), BUO002 (Boussinesq out of range), BUO003
  (no energy equation or no β), BUO004 (gravity with energy but no model),
  BUO006 (a static-kind pressure value over a tall opening); the
  variable-density model is refused with ENE006 until stage 7c.
- **Cases and docs.** `buoyancy_enclosure_plume`, `buoyancy_box_at_rest`,
  `buoyancy_slot` (analytic `buoyant-slot`, observed orders 2.14 and
  2.13), `buoyancy_channel_reduced`, `buoyancy_channel_static`,
  `buoyancy_open_box`; gallery and recipe 186 ("Switch on gravity and
  buoyancy; which pressure do my boundary values mean").
- **Infrastructure features can reach verified (D199).** Restart, logs and
  run record, Python package, command line, tool registry, journal and
  replay and the MCP server are graded through contract evidence: a row
  marked `infrastructure = true` lists its contracts (bitwise, schema,
  equivalence, replay, control, diagnostics, round-trip) and each needs a
  test of kind `contract`. Three new rows (tool registry, journal and
  replay, MCP); 48 of 97 features are verified, up from 41 of 93. New test:
  every row of `run.jsonl` follows its schema.
- **Operator exactness check of the mesh (D236).** Before the solve the
  velocity gradient and its gradient (the velocity Hessian of the
  curvature correction) are applied to linear and quadratic test fields on
  the mesh of the run. Cells where they are not exact are flagged:
  `MESH004` (a broken stencil), `MESH005` (with the curvature correction:
  it does not converge there), `MESH006` (with the QCR: the stress
  divergence does not converge there) and the limitation `LIM008`; the 50
  worst cells per operator go to `flow.json` and the `case.check` report,
  every cell to `operators.csv`. New key `solver.operator_tolerance`
  (0.01). No result changes. Gallery and recipe 189.
- **Cycle-averaged coupling of a periodic flow with its solids (D204).** `[energy.cycle]` runs the loop by itself on one mesh: a fluid phase holds the solids at their wall temperature and averages the heat transfer of every interface face (and the mass flux) over a cycle once the flow is periodic; a solid phase solves the energy equation steady with the fluid convected by the cycle-mean mass flux and the difference between the transient and the mean-flow heat flow added to the solid cells, so the wall answers its own temperature with the response of the flow; the new wall temperatures go back to the fluid phase, with a quasi-Newton (IQN-ILS) update after the first, until the wall temperature and the heat flow settle. On a pulsating channel over a heated wall the loop converges after 4 updates and 10 cycles; heat flow and wall temperatures are within 0.4 % of the rise of the fully coupled periodic state. The classic first-cell Robin solid phase stays as `solid_model = "robin"` for comparison: it converges no sooner than the fully coupled run there. `enabled = false` runs the same case fully coupled. Outputs `cycle_coupling.csv`, `cycle_faces.csv` and `transient.cycle_coupling` of the report; heat balances of every fluid step and solid phase to 10⁻¹³. Codes CYC001, CYC002, CYC004, CYC005 and LIM009; checks V7.39, V7.40; gallery and recipe 190.
- **Minimal path to v1.0.** Requirements section 7a and D235 record what v1.0 contains and what moves to v1.1.
- **Multi-rank MPI sized (spike).** Assessment `docs/assessments/mpi-multi-rank-sizing.md`: what exists, a test that shows two ranks do not yet give the matrix of one rank (local face orientation and face order), the canonical construction that does (bitwise on 1 to 4 ranks, hexahedra and polyhedra, new example `mpi_matrix` in the MPI job of CI), the measured cost of exact reductions and of MPI messages on one machine, and a sized plan of nine to ten pull requests. No solver change.
- **Turbulent pipe flow with grid convergence (D45).** New case
  `pipe_turbulent_sst.toml`: fully developed pipe flow at Re_b = 5 300 on a
  0.5° wedge with k-ω SST (`low-re`); friction velocity u_τ/U_b = 0.0714,
  +4.5 % from the DNS of El Khoury et al. (2013) (0.0683), GCI_fine 0.2 %
  on three meshes (`examples/studies/pipe_gci.py`). Gallery and recipe 193.
- **QCR stress divergence consistent on irregular meshes (D236 step 2).**
  On the cells where the operator check finds the face-neighbour velocity
  gradient first order, the explicit QCR stress takes a quadratic
  least-squares velocity gradient over the neighbours of the neighbours,
  and the face value of the stress is corrected for skewness with cell
  gradients of the stress. The stress force of the manufactured solution
  now converges on meshes with 15 % perturbed vertices (hexahedra 4.6 %,
  1.5 %, 0.65 %; tetrahedra first order) instead of staying at about 5 %.
  `MESH006` is now reported only for cells where the quadratic fit is rank
  deficient. QCR cases on meshes with flagged cells change slightly (the
  square duct, regular hexahedra with 286 flagged wall cells: largest secondary velocity 1.558 % of U_b, before 1.558 %). Gallery and recipe 194.
- **Low-Mach variable density ρ(T, p₀) (stage 7c, D151).** `[energy]
  density = "variable"` (or `[buoyancy] model = "variable-density"`, which
  switches it on) takes the density of the fluid's density form, an ideal
  gas ρ = p₀ M/(R T) or a table or law of the temperature, in every cell
  and face. The thermodynamic pressure p₀ is the operating pressure in a
  domain open to an ambient, or with `thermodynamic_pressure =
  "closed-volume"` follows from the fixed mass of a closed domain, which
  then stays to round-off. Works with SIMPLEC and the coupled algorithm,
  steady and transient, with the turbulence models and on the
  fluid side of conjugate heat transfer; the kernels are shared by CPU and
  GPU. Cases without the switch are bitwise unchanged.
- **The equations of a variable density.** The mass storage V ∂ρ/∂t in
  the pressure equation, the momentum time term with the density of each
  time level, the full viscous stress with −⅔ μ ∇·u, the energy equation
  with ρ per cell and, in a closed volume, the pressure work dp₀/dt
  (`[energy] pressure_work`, default `thermodynamic` there and `none` in
  an open domain), so that a closed gas heats with c_v.
- **Ledgers.** The mass account books the storage of the density, the
  kinetic-energy account the pressure dilatation p ∇·u as a physical term,
  the heat account the storage ρ ∂h/∂t + h ∂ρ/∂t; the report and `history
  .csv` carry p₀, the mass and the density range, and a mesh sweep gives
  the observed order and GCI of p₀ and of the wall heat flows.
- **Checks.** ENE004 warns above Mach 0.3; ENE009 refuses a constant
  density form, a closed volume in an open domain or without the switch,
  and the Boussinesq model with it; BUO005 refuses a closed domain with a
  constant p₀; ENE006 names what does not work with it yet (passive
  scalars, a restart, periodic meshes, convective outlets, mass-flow
  inlets, wall mass fluxes, rotating zones). BUO002 points to the new
  model.
- **Verification.** A heated duct against its exact one-dimensional
  balances (V7c.1), a closed box heated uniformly that warms with c_v and
  p₀(t) exact (V7c.2), a closed box heated on one side that keeps its mass
  (V7c.3), and a manufactured solution with ρ = 1/T in a closed volume
  (V7c.4, generated by `tools/verification/vd_mms.py`), second order for
  velocity and temperature on hexahedra and polygons, SIMPLEC and the
  coupled algorithm on the same discrete solution.
- **T-junction mixing on polyhedra with grid convergence (D45).** New case
  `tjunction_mixing.toml`: laminar T-junction (Re 200 and 100) with a passive
  scalar, Gmsh tetrahedra converted into polyhedra, three meshes of 4 852,
  16 665 and 44 738 cells (`examples/studies/tjunction_gci.py`). The main-inlet
  pressure converges at order 1.77 (GCI 2.4 %); the branch-inlet pressure does
  not (order 0.2, GCI 29 %), and the conversion fails with GEN001 on some
  meshes; both are written as gaps. Validation against measurements is
  pending. Gallery and recipe 196.
- **Curvature correction consistent on irregular meshes (D236 step 3).**
  DS_ij/Dt of the curvature correction of SST takes S_ij from a velocity
  gradient that is second order everywhere (the quadratic least-squares
  fit of the QCR on the cells where the face-neighbour gradient is first
  order), and on the cells where the least-squares gradient of S_ij is not
  consistent, the gradient of S_ij from a quadratic least-squares fit over
  interior cells. The manufactured solution converges at order 2 on
  regular and 15 % perturbed tetrahedra (error of k 0.54 %, 0.14 %,
  0.036 %) instead of staying at about 1.3 %; a channel in pure shear stays
  unchanged to round-off. The operator check reports
  `quadratic_hessian`; `MESH005` remains only where the fit fails, and a
  run whose cells take the fit writes `operators.csv` and the check in
  `flow.json`. Gallery and recipe 197.
- The face-based `.msh` reader is named after its format, Fluent `.msh` (`mesh.format = "fluent_msh"`; the old key `industry_msh` still works). Naming commercial codes is now allowed where appropriate (AGENTS.md).
- **Asymmetric plane diffuser (D45).** New cases `diffuser_sst.toml` and
  `diffuser_sst_sharp.toml`: the diffuser of Buice & Eaton (10°, area ratio
  4.7) at Re_H 20 000 (centreline velocity) with k-ω SST and a developed
  channel profile at the inlet (new case `diffuser_channel_profile.toml`,
  table by `examples/studies/diffuser_inlet_profile.py`), on structured
  meshes from `examples/meshes/diffuser.py` with the corners rounded as in
  the experiment (4.3 H) or sharp (`--sharp`, cells shortened towards the
  corners). With rounded corners the case converges on three meshes (GCI in
  gallery 199, `examples/studies/diffuser_gci.py`): the reattachment point
  (28.94 H on 111 360 cells) has a GCI of 0.36 %, the separation point
  (2.20 H) 4 % and the pressure recovery 0.03 % or less. With sharp corners
  (three meshes, 8 400 to 134 400 cells) the reattachment point (29.03 H)
  has a GCI of 0.51 % and the pressure recovery 0.15 % or less, within
  1.1 % of the rounded ones; the separation point converges at order 0.6
  only (GCI 9 %), since the wall shear stays near zero for 1.5 H behind
  the corner. No comparison with measurements yet (pending). Gallery and
  recipe 199.
- **Why the diffuser did not converge before.** A sharp wall corner makes
  the least-squares gradient of thin wall cells wrong by O(1), and linear
  upwind from the first iteration diverged on it; trapezoidal wall cells
  switched the D50 remediation on, on the coarse levels only; a uniform
  start broke continuity in the widening; and relaxation 0.5 only looked
  like a stall. The case now starts from potential flow with upwind for
  200 iterations, runs at relaxation 0.8, and judges convergence on probes
  in the log (`docs/pitfalls.md`).
- **The momentum scheme may change during a run.** `solver.momentum.scheme`
  may be scheduled (`[[schedule]]`) or re-read (`control.toml`) between
  `"upwind"` and `"linear-upwind"`: a first-order start, then second
  order. A change to or from `"central"` is refused.
- **A transient run may start from a steady restart file.** `--restart`
  with a steady restart file and a transient case starts the transient
  run from the steady fields at `time.start_time` (no older time levels);
  the log says so.
- **Multi-rank MPI, step W1: canonical local mesh.** The local mesh of every rank now has bitwise the geometry of the one-rank mesh, on any number of ranks and any partition (D93 (1), (3)): faces stored flipped are marked (`Mesh::face_flip`) and evaluated in the node order of the global mesh, every cell lists its faces in global face order, halo cells take their geometry from the owner. Checked in tests on 1 to 4 ranks, METIS and RCB, hexahedra, tetrahedra, polyhedra and mixed, 2D and 3D, after renumbering, and across MPI processes in CI (`mpi_matrix`). `foehn mesh` with a partition reports the flipped faces and the differing faces and cells per rank (always 0, the run fails otherwise) in `mesh_report.md` and `mesh_report.json` (`canonical`) and writes the cell field `flipped_faces`. Runs without a decomposition do not change. Gallery 200, recipe 200, case `mesh/canonical_ranks.toml`.
- **Fix: the bounded face velocity without a limiter crashed (D50).**
  `solver.face_velocity = "bounded"` set by hand with no momentum limiter
  read an empty limiter array and ended in a segmentation fault; the
  limiter's per-cell epsilon is now set up whenever the face velocity uses
  it, and an empty array is an error instead of a crash. The automatic
  remediation was not affected.
- **The skin swing, fully coupled cycles at the end and the restart of conjugate runs (D204).** `[energy.cycle] skin = "effusivity"` gives the held wall of the cycle-averaged loop the swing of a thin skin: the harmonics of the face heat flux over the averaged cycle times the surface response of a semi-infinite solid of effusivity √(k ρ c_p), with the swing amplitudes inside the quasi-Newton update; CYC003 warns where the skin depth reaches through the wall. `resolved_cycles` ends the loop with fully coupled cycles, started from the damped and delayed swing of the skin in depth, and reports the exact swing per face (`resolved_t_wall_min`, `resolved_t_wall_max`). On a thin-skin wall the loop ends at 1250 s within 0.016 K and 0.1 % on the swing of the periodic state, which the fully coupled run reaches only after 4400 to 5550 s. Conjugate runs now restart: `restart.h5` holds the temperatures of the solids with their old level and the state of the loop, and continuations are bit-identical in any phase of the loop; a restart with `[energy.acceleration]` is still refused. The heat ledger books what the swinging held cells store. New columns `swing_max`, `swing_change` in `cycle_coupling.csv` and `t_wall_amplitude`, `t_wall_min`, `t_wall_max` in `cycle_faces.csv`; gallery and recipe 202.
- **Engineering gallery: shell and tube heat exchanger mesh (application
  203).** `examples/applications/shell_tube_hx/mesh.py` turns a STEP
  assembly (housing, seven tubes, two tubesheets, three segmental baffles)
  into one conformal multi-region tet mesh with Gmsh: hot air on the shell
  side, water on the tube side and four metal zones, the four nozzles
  closed as named inlet and outlet faces, every interface shared node for
  node (983 249 cells at level 2). A new page of the engineering gallery,
  heat exchangers, holds the entry with a section and a cut-away of the
  mesh; recipe 203 gives the steps for one's own assembly. The flow run
  waits for two fluids with their own materials in one conjugate case
  (ENE006).
- **U-bend (D45).** New cases `ubend_sst.toml` and `ubend_sst_sharp.toml`:
  the two-dimensional U-duct of Monson et al. (1990), Rc/H = 1, at
  Re_b = 10⁵ with k-ω SST and a developed channel profile at the inlet (new
  case `ubend_channel_profile.toml`, table by
  `examples/studies/ubend_inlet_profile.py`), on structured meshes from
  `examples/meshes/ubend.py` with the rounded turn of the experiment or with
  sharp edges (`--sharp`: a square divider end and end wall). The rounded
  case converges on three meshes from a potential-flow start (GCI in
  gallery 204, `examples/studies/ubend_gci.py`): the inner wall separates at
  θ = 139.8° (GCI 0.19 %) and reattaches 2.19 H behind the turn (GCI 2.9 %);
  the wall pressures have a GCI of 0.27 % or less. The curvature correction
  of SST moves the separation to 127.5° and the reattachment to 2.57 H. The
  sharp-edged duct converges on the two coarser meshes with a gentler start
  (relaxation 0.5, linear upwind from iteration 601, relaxation 0.8 from
  1 201). No comparison with measurements yet (pending). Gallery and
  recipe 204.
- **Engineering gallery: flow cases of the shell and tube heat exchanger
  (application 203).** `case_coflow.toml` and `case_counterflow.toml` set
  the owner's conditions on the mesh of application 203: hot air at 300 °C
  and 15 m/s on the shell side against 1.5 bar, water at 20 °C and 5 m/s
  on the tube side against 1 bar, in coflow or counter flow, carbon steel,
  the housing held at 25 °C outside, k-ω SST, everything from rest. They
  run once several fluid materials share one conjugate case (D194).
- **Several fluids in one conjugate case (D194, stage MF-A).** Fluid regions take a material of their own (`[regions.<name>] material`, a liquid or a gas) and the temperature at which it is evaluated (`temperature`); density per cell and face, viscosity per face, c_p, k and μ(T) per fluid in the energy equation, ν and ρ per cell in the turbulence models and wall functions, and the density per face in the boundary conditions. No fixed limit on the number of regions (32-bit region index). Fluid cells that touch only through solids form separate fluid parts: each needs its own pressure boundary (REG003), two different fluids sharing a face are refused as a two-phase interface (REG004), the log lists the parts (REG005), and settings not yet built for several fluids say so (REG006: variable density, gravity, rotating zones, passive scalars, convective and periodic boundaries); the coupled solver takes the density per cell in its pseudo-time term and continuity scale. The heat ledger books each fluid's enthalpy with its own h_s. Check V7.33 (planar form): a water / oil counterflow exchanger across a steel plate gives ε = 0.7756 against ε–NTU 0.7751 (gallery and recipe 206).
- **The frozen-flow march and the warm-up history of the cycle loop (D118, D204).** `[energy.frozen_flow]` speeds up the thermal warm-up of a transient conjugate run: coupled steps until the flow has settled (a window test on the velocity with Aitken's extrapolation), then the temperatures of fluid and solids are marched alone with a longer `time_step` while velocity, pressure, mass flux and turbulence are held. A flow that follows the temperature (a viscosity form of T, buoyancy) is refreshed by coupled steps whenever a fluid temperature has moved by `refresh_change` (or every `refresh_every` frozen steps). A heated wall under a fluid whose viscosity halves over 17 K warms up in 557 steps instead of 6000 (17 s against 155 s), within 0.06 K of the fully coupled run; a flow that does not depend on the temperature is held once and matches to the time error of the long steps. `result.json` reports the phases, refreshes and solves per step; codes FRZ001 to FRZ004 and LIM010. `[energy.cycle] solid = "transient"` (with `march_cycles`, `march_steps`) turns the cycle-averaged loop into a warm-up history: each update marches fluid and solids under the cycle-mean flow and jumps the clock by whole periods; a pulsating channel warms up in 780 steps instead of 6000, within 0.6 % of the rise of the fully coupled cycle means. `cycle_coupling.csv` gains `march_time` and `march_jump`. A conjugate restart now also keeps μ(T) on the faces (a restart in a refresh was not bit-identical without it) and the state of the march; restarts are bit-identical in a frozen phase, in a refresh and after a jump of the clock. Gallery and recipe 208.
- **U-bend with sharp edges: third mesh (gallery 204).** The finest of the
  three sharp-edged meshes (201 600 cells) ran to its limit of 8 000
  iterations without converging (the separated region behind the divider
  end still grows slowly); with it, the reattachment has a GCI of 0.85 %
  and the end-wall pressure below 0.01 %, while the pressure 10 H behind
  the turn moved 4.6 % from the middle mesh. The level-2 values are marked
  provisional.
- **U-bend images drawn as in the literature (gallery 204).** The field
  images of the U-duct now show the flow entering at the upper left and
  leaving at the lower left, with inlet and outlet arrows; the case, the
  mesh and the results are unchanged.
- **Bounded convection of the scalars and the energy equation (D237).** `limiter = "mlp"` in `[solver.scalar]` and `[solver.energy]`: the multi-dimensional limiting process of Park & Kim (2012) bounds the linear-upwind reconstruction at every vertex by all cells that share it, with Nishikawa's R3 function (also R4, R5, Barth–Jespersen and Venkatakrishnan as `limiter_function`) and a threshold in the units of the field, ε = c (φ_max − φ_min) (`limiter_epsilon`, 0.01) or per cell (`limiter_epsilon_scaling = "cell"`); `limiter_relaxation` (0.5) damps the limit cycle of a steady run. Also new: the face limiter `barth-jespersen` (also for the momentum and the `foehn scalar` solver), the face clamp `bounded` and `verification.scalar_range`. On the oblique step on perturbed polyhedra c stays within [0, 1.0021] (Venkatakrishnan [−0.0175, 1.0178]); the air of the shell and tube exchanger stays within 0.6 K of its inlet temperature (Venkatakrishnan 12.6 K above). Venkatakrishnan stays the default; MLP is set in the T-junction (its scalar now converges with the mesh, order 2.1, GCI 11 %, gallery 196), the two heat-exchanger cases and the stirred-tank tracer. Refused for the momentum and turbulence (CASE004). Gallery and recipe 211.
- **Shell and tube heat exchanger with two fluids (application 203, D194).** The exchanger runs hot air (300 °C, 15 m/s, out against 1.5 bar) on the shell side, water (20 °C, 5 m/s, out against 1 bar) on the tube side and carbon steel in one steady conjugate case, coflow and counter flow, on the level 0 mesh (477 921 tetrahedra) with k-ω SST, the coupled solver and the vertex limiter of D237. The air gives up 1.72 kW in both cases and leaves at about 124 °C (ε = 0.63, NTU 0.99); with the water carrying 2 600 times the heat capacity flow of the air (C_r = 3.8·10⁻⁴) the flow arrangement drops out, as ε-NTU says. A demonstration: one mesh, and the water-side heat flows swing by a few percent with the residuals stalled near 10⁻³. New `flow_images.py` draws the temperature and speed sections; recipe 203 gives the run commands; a README teaser.
- **Steel temperatures of the shell and tube exchanger (application 203).** Two new images show the housing, the tubesheets, the tubes and the baffles coloured by temperature (20 to 30 °C), with a table of each part in coflow and counter flow: the steel follows the water, within about 4 K of its 20 °C, and is warmest near the air inlet (the housing at the inlet nozzle 27.7 °C, the rim of the first baffle 27.3 °C). The coflow images and the README teaser are redrawn from the state at iteration 3 250: the earlier ones were drawn one iteration after a restart, which reset the memory of the vertex limiter and let the air in the inlet nozzle overshoot to 360 °C for that iteration. The limits now say that the water in the end chambers falls to 15.7 °C in a few cells, below its 20 °C inlet.
- **Restart with the vertex limiter (D237).** With `limiter = "mlp"` the relaxation (`limiter_relaxation`) carries the limiter value ψ from one outer iteration to the next, so ψ is part of the state: restart.h5 now holds it (`scalar_psi_<name>` for the passive scalars, `cht_psi` for the energy of a conjugate case), and a run continued from it ends bit-identically to one without the interruption. Before, a restart started from ψ = 1 and its first iteration overshot: the air in the heat exchanger's inlet nozzle went up to 360 °C (inlet 300 °C) and settled again within 150 iterations. A file written by an older version still reads and starts from ψ = 1.
- **Start values per region and per fluid part (D238, requirement).** Each domain of a multi-region case gets its own start values (`[regions.<name>.initial]`); by default every fluid part starts at the pressure of its own pressure boundaries (a Laplace solve between several values), with the operating pressure per part from MF-B. Later: start values on volume selections of the mesh (cell zones, named cell sets, the regions of D73). Specified in stage 7 §9c.4 with codes REG007–REG010 and checks V7.45 and V7.46; not built yet.
- **Pressure values absolute by default, gauge optional (D239, requirement).** From case format 0.7 every pressure value in a case file is an absolute static pressure unless `[fluid] pressure_values = "gauge"`; the solver keeps solving the gauge part against an operating pressure that each fluid part takes from its own boundaries, and older files keep their gauge meaning (and their results) through a new changed-default rule of `foehn case upgrade`. Not built yet.
- **Roadmap: 7b2 moves to v1.1.** The laminar natural-convection benchmarks (square cavity Ra 1e3 to 1e6, vertical plate, mixed convection) and the heated D45 cases with SGDH are no longer part of v1.0 (maintainer decision 2026-11-23, D235, section 7a). The cavity will be graded on its own GCI and energy balance, because the published cavity values are numbers from other codes (D201).
- **Multi-rank MPI, step W2: exact reductions.** Sums, dot products, residual norms and min/max no longer depend on the order of the terms, the thread count or the number of ranks (D31 (1), D93 (4)). `foehncraft_core::exact::ExactSum` is a long accumulator over the whole binary64 range (Kulisch 1971; Collange et al. 2015), rounded once to nearest even, with NaN and infinities counted apart and a batched path (one integer bin per exponent). `Communicator::sum_exact`, `sum_exact_many` and `all_sum_i64` reduce it across ranks as one all-reduce of 71 integers; `min` and `max` are order-independent (`min_total`, `max_total`). The CPU backend and the flow CPU device compute dot products, norms and min/max exactly, the AMG eigenvalue estimate too, and about 30 sums over cells, faces and patches of the flow crate (balances, torques, mass flows, volumes, sources) are exact. One-rank results change once in the last bits against v0.7; no reference band is affected. **GPU:** the GPU keeps its block sums until v1.1 (owner decision of 2026-11-23), so CPU and GPU runs now agree to rounding, not bitwise. Partitioned `foehn mesh` runs check sums over the ranks against one rank (`reductions` in the report); the example `mpi_reductions` checks the same across 1 to 4 MPI processes in CI. Measured cost (D93, 4-core Xeon of the cloud session, two runs each against main): exact dot product 2.3 to 3.5 ns per value on one thread against 0.6 to 0.9 ns for the block sum, 0.7 to 0.9 against 0.3 to 0.6 ns on 4 threads (`exact_reduction_cost` example); per outer iteration +2 to +7 % on the cavity with 970 225 polyhedra (30 SIMPLEC iterations, 4 threads) and +18 to +25 % on the 64 × 64 cavity (4 096 cells, one thread, where the BiCGStab dot products weigh most); the same iteration counts and results to the printed digits. Gallery 219, recipe 219, case `mesh/exact_reductions.toml`.
- **Pressure values absolute by default, gauge optional (D239, case format 0.7).** Every pressure value of a case file (boundary `pressure` and `total_pressure`, numbers and expressions, `[initial] pressure`, `[pressure] reference_value`) is now an absolute static pressure; `[fluid] pressure_values = "gauge"` keeps values relative to the operating pressure. Each fluid part takes its operating pressure from its own pressure boundaries (static outlets first), or from `[fluid]` / the new `[regions.<name>] operating_pressure`; the reader turns the case into its gauge twin, so the solver is unchanged and an absolute case runs the same bits as its gauge twin, and each gas's density follows the pressure of its own part. The `pressure` field (with the VTKHDF attributes `pressure_values` and `operating_pressure`), probes and the mean pressures of `boundaries.csv` are written in the file's convention; the log (PRS005) and `flow.json` (`pressure_values`) state the convention and the operating pressure of every part. Checks PRS001 (absolute value ≤ 0), PRS002 (a gas below 1 kPa looks like a gauge value), PRS003 (operating pressure from expressions only), PRS004 (one absolute start pressure for parts of different operating pressures). The case-format migration has a second kind of rule, a changed default: files of format 0.6 and older read as gauge values (a note) and keep their results; `foehn case upgrade` writes `pressure_values = "gauge"` into them. New case `cht_two_fluids_absolute.toml` (water at 1 bar, oil at 1.5 bar). Gallery and recipe 220.
- **Start values per region and per fluid part (D238).** `[regions.<name>.initial]` sets `velocity`, `pressure`, `temperature` and `k` with `omega` or `epsilon` of a region (numbers or expressions of x, y, z; a solid takes `temperature` only), before `[initial]` and `[energy] initial`; `temperature = "inlet"` starts a fluid region at the inflow-weighted inlet temperature of its fluid part. Without a start pressure every fluid part now starts from its own pressure boundaries: at their value, or on a Laplace solve between two different values (a total-pressure inlet and an outlet), so a part driven by a pressure difference starts on its pressure line; not with `[initial] pressure`, a reference cell or gravity, so the cases of before keep their start. The log (REG010) and `flow.json` (`start_values`) list every start value and its source; checks REG007 (a key a region table cannot take), REG008 (a value not > 0), REG009 (an `"inlet"` temperature without an inlet). New case `cht_two_fluids_start.toml` (2251 instead of 2616 iterations from a uniform start). Gallery and recipe 221.
- **Docs: multi-rank MPI moves from v1.0 to v1.1 (D235, maintainer decision 2026-11-23).** Section 7a of the requirements, D235, the stage 8 plan and the MPI sizing note say so; W1 (#200) and W2 (#219) stay on `main` as groundwork. The README lists the internal-flow cases of D45 under "What works today", says that GPU results now agree with the CPU to rounding (since W2), and describes v1.0 as runs in one process.

## 0.7 (2026-11-22) – stage 7 (part): energy equation, conjugate heat transfer and agents


**Highlights.**

- **Energy equation of the fluid.** Conservative sensible-enthalpy temperature
  equation with c_p(T), k(T) and μ(T) from the material database (water and
  air tables, NASA-7/9, CHEMKIN reader, user database files), viscous
  heating, volumetric sources, thermal wall conditions (`temperature`,
  `heat_flux`, `heat_flow`, `heat_transfer_coefficient`) and a heat balance
  of every boundary.
- **Conjugate heat transfer.** Solid regions in steady and transient,
  laminar and turbulent flow: one energy matrix over fluid and solids with
  series conductances at material faces, the thermal wall function of the
  flow on the fluid side, contact resistances and thin walls, orthotropic
  solids, region heat sources, per-step `heat.csv` and probes of the
  temperature. Automatic thermal acceleration of solids in transient runs
  (D213) with a conservative restore. Exact checks V7.4, V7.6, V7.36,
  V7.37 and the acceleration checks V7.41 and V7.42.
- **Agents first.** The tool registry (D196): every command of `foehn` and
  every Python function runs through typed tools with schemas, a journal
  with replay (D197), the MCP server `foehn mcp` with read and write roots,
  run limits and study budgets (D198), the live control port with
  `run.status`, `run.stop` and `run.set` (L19), run limitations
  LIM001–LIM007 from the verification status page (D199), and an
  assistant-run heat sink study as gallery 181.
- **Robustness and clean room.** Bounded Rhie–Chow face velocity for
  skewed tetrahedra (D50), `unsafe` only in named modules and a
  broken-input corpus for every reader (D202), the ζ-f and ζ-f₀ wall sink
  that converges with the mesh, manufactured solutions of every
  turbulence model including ζ-f₀, curvature correction and QCR, and the
  clean-room rework of the reference numbers (D201, E1–E10).
- **Engineering gallery.** The stirred tank and the 2 L tank of Fyfe et al.
- **Plans and requirements.** Cycle-averaged CHT (D204), design
  optimisation (step 26), the immersed boundary method (L13, with the
  owner decisions), the adjoint solver (step 27, D220–D224) and the log
  redesign (D214). No code.

**Known limits.** Buoyancy (7b1) and variable density (7c) are not in this
release; periodic meshes with solids are refused (ENE006); the curvature
correction and the explicit QCR stress do not converge at full order on
tetrahedra (see the grades on the verification status page); the GPU path
of the energy equation is covered by host-side changes and was not run on
a device in this release; the generated verification manual was not
regenerated for 0.7.

**Changes.**

- **Engineering Gallery started (application 131).** A new track for industrial cases with complex geometry, next to the feature gallery: the index [docs/applications.md](docs/applications.md) with the check legend and trust levels, the rotating-machinery page, and the first entry, a baffled stirred tank with a three-blade hydrofoil (power number 0.29–0.30, pumping number 0.45 at 2 500 iterations). The entry is a demonstration: the steady frozen rotor does not settle, and the check card says so. The mesh script and the case are in `examples/applications/stirred_tank/`.
- **Mesh grading for the tank.** The tank mesh is graded from the walls, which keeps the blade–hub and sparger corners below skewness 1.1; without it they reach 1.4 and the pressure runs away.
- **Multi-zone SU2 and CGNS meshes for rotor/stator runs (D195, gallery 138,
  recipe 138).** `mesh.file` reads all zones of an SU2 file (`NZONE= n`) or a
  CGNS base into one mesh, each zone a cell zone, the nodes of different zones
  apart, so the interface is two boundaries for `[[sliding]]`. New key
  `mesh.merge_nodes` joins coincident nodes for touching stationary zones.
- **ζ-f and ζ-f₀ converge with the mesh at the wall (stage 6, maintainer
  decision 2026-11-20).** The sink of the k equation is now
  max(k/T_K, 2νk/d²), with d the wall distance and T_K the time scale
  with its Kolmogorov bound. It equals ε at the wall, so k falls like y²,
  and it is k/T_K away from walls, which keeps k positive at a laminar
  leading edge. The sink k/T_K of 6d2 let k fall like y only, and the
  wall value of ε then grew with every refinement: the converged channel
  of ζ-f₀ went from U_b⁺ 17.79 to 18.39 over four meshes. Now it gives
  17.34, 17.54, 17.58, 17.58 (+0.2 % from the DNS on 768 cells across).
  ζ-f₀ takes C_ε1 = 1.4 (1 + 0.009/ζ) in place of the printed
  1.4 (1 + 0.012/√ζ), which gives U_b⁺ −11 % with a consistent sink. The
  flat plate (C_f +1.8, +0.2, −0.2 % from White on three meshes, GCI 0.2 %) and the
  step (x_r/H 6.42) stay in their bands. The impinging jet converges with
  the mesh, but at Nu₀ = 196, +38 % from the experiment (the
  stagnation-point excess of k; finding). The ζ-f plate now converges
  (C_f +0.05 % from White on the middle mesh; before, it stopped at the
  iteration limit at +24 %).
- **New keys:** `turbulence.k_sink` (`"wall"`, the default; `"epsilon"`,
  the papers' −ε; `"kolmogorov"`, the old ζ-f₀ default) and
  `[turbulence.constants] c_epsilon1_power` (1, or ½ for the printed
  C_ε1″/√ζ of Popovac 2023).
- **Fix:** the ambient (sustaining) source of k in ζ-f and ζ-f₀ is the
  sink at the ambient state, so the free stream holds its k for every sink
  form.
- **The channel runs of the GCI study use the coupled algorithm**, at
  most 3 000 iterations each: ζ-f₀ converges to 10⁻⁹ in 436 to 883
  iterations instead of up to 35 000 with SIMPLEC, the k-ε models in
  under 700, ζ-f in under 2 600, with the same U_b⁺. SST on the finest
  mesh and EB on the two finest stop at the limit with a residual near
  10⁻⁷, within 0.3 % of their converged values (finding; the docs keep
  the SIMPLEC values of #128). The plate runs of ζ-f and ζ-f₀ use it too
  (178 to 496 iterations).
- **The channel case files of ζ-f and ζ-f₀ use the coupled algorithm**
  (about 230 and 370 iterations instead of 2 100 and 3 700).
- **The ζ-f plate converges** (`require_convergence = true`): C_f +2.0 %
  from White on 64 × 48 cells, +0.03 % on 128 × 96 (before +5.1 % with
  the k floor of the leading edge).
- **The inlet table of the impinging jet** (`impinging_jet_profile.csv`)
  is made again from the pipe flow with the new sink.
- **Gallery 141 and recipe 141**: the sink of k with two field images
  (the ζ-f plate near its leading edge, the wall layer of the jet).
- **Finding:** in the test of a turbulence patch carried by the rotor
  (`tests/sliding.rs`), the imprinted run of ζ-f₀ at a Courant number of
  1.7 drives ζ to its cap and diverges with any sink of k that is ε at the
  wall (the new default and the papers' −ε); k/T_K held there. The
  imprinted run of the test now takes half the time step (Courant 0.85).
- **Robust on skewed tetrahedra: bounded Rhie–Chow face velocity (D50).** On raw Gmsh tetrahedra, turbulent rotor cases (the stirred tank) still diverged with the cell-quality remediation. The cause was the face velocity of the mass flux, which is moved to the face centroid with the unlimited cell gradients. On skewed faces that feeds a velocity mode back into continuity with a gain above one. The remediation now also sets `solver.face_velocity = "bounded"`, the mean of the two limited extrapolations to the face centroid. It is exact for a linear field and stays between the values of the two cells. Meshes without skewed faces are unchanged. The new opt-in `solver.transposed_stress = "trace-free"` drops the grad-div part of the explicit transposed stress in turbulent flow.
- **Kernels: density per cell and per face (stage 7.0b).** The flux,
  momentum, scalar and turbulence kernels read ρ from cell and face fields
  instead of one number, on the CPU and the GPU. Results are bit for bit
  unchanged; the fields carry the variable density of stage 7.
- **Manufactured solutions of every turbulence model (stage 6
  specification 9.2).** `turbulence::mms` solves the turbulence equations
  of SST, the standard and the realizable k-ε, k-ε-ζ-f and the
  elliptic-blending k-ε alone in a frozen cellular flow with smooth
  manufactured k, ω/ε, ζ/φ, f/α and wall distance and every source term
  active; the L2 errors converge at order 1.8–2.3 on hexahedra, polyhedra
  and triangular prisms. New: `Flow::turbulence_iteration`,
  `Turbulence::set_added_zeta_sources` (with an added right-hand side of
  the f/α equation in the kernel) and `Turbulence::set_zeta_fields`.
- **Fix: realizable k-ε C_μ with the deviatoric strain.** W̃ of Shih et
  al. now comes from the trace-free strain; with the discrete trace it
  jumped across its range where the strain vanishes.
- **Requirements: oriented anisotropic material properties (L48).** Material
  frames for orthotropic conductivity and permeability in laminations,
  windings and magnets of e-motors; recorded with L30. No code.
- **Assessment: the ArtCraft development method.** `docs/assessments/artcraft-development-method.md` compares how the ArtCraft Crafting Apps are built with Foehncraft's practice and lists open questions for the maintainer (tool registry, mandatory roots for agents, `AGENTS.md`, cost gate, verification grades, clean room). No behaviour changes.
- **Materials: thermal data (stage 7.0c).** Every material may carry a
  specific heat c_p and a thermal conductivity k (any property form, or
  orthotropic with two axes). Solids need density, c_p and k (MAT006).
  User material database files (`[materials] include`, `FOEHN_MATERIALS`)
  extend the built-in database; `foehn materials --include FILE` lists them.
  Water and air have c_p and k tables from the IAPWS-95 and IAPWS 2011 (water)
  and Helmholtz-equation (air) formulations, evaluated with CoolProp 8.0 (MIT).
  NASA-7 and NASA-9 specific-heat forms and a CHEMKIN reader are in.
Record the ArtCraft challenge decisions D196-D203 and the L7 GUI rules; docs only.
- **Energy equation of the fluid (stage 7a1).** `[energy]` solves the
  temperature of the fluid, conservative in the sensible enthalpy, with
  c_p(T), k(T) and μ(T) from the material or from `[fluid] specific_heat`
  and `conductivity` (μ(T) reaches the momentum equation;
  `temperature_dependent = false` keeps the reference values). Viscous
  heating, volumetric sources and exact solutions. Walls take
  `temperature`, `heat_flux`, `heat_flow` (W) or
  `heat_transfer_coefficient` with `ambient_temperature` in
  `[boundary.<name>.energy]`; inlets take `temperature`. The report has a
  heat balance and the heat flow of every boundary. The condition follows
  from the keys (no `type` key). The thermal set-up check gives THM001 to THM005 and THM009 to
  THM011 (undetermined temperature level, heat with no way out,
  conditions on the wrong boundary kind, values out of range, adiabatic
  walls listed, tables without `[energy]`, malformed tables, inlets
  without a temperature); ENE006 for settings this version does not
  have. Checks V7.1 and V7.2 at second order (gallery and
  recipe 154).
State in D203 that Python stays the user front end.
- **Solid regions and conjugate heat transfer (stage 7d1).** A region of
  kind `solid` (`[regions.<name>] kind = "solid"`, its cell zones and a
  solid material with density, specific heat and conductivity) now joins
  the energy equation: the flow runs on the fluid cells, where the faces
  against a solid become the no-slip wall `interface-<region>`, and the
  temperature runs on all cells as one matrix, with the series conductance
  of the two half cells at every face between materials. The report gives
  the heat flow and mean temperature of each interface; the heat balance
  of the whole domain closes to round-off. The box generator takes
  `zones = [{ name, lower, upper }]` to mark cell zones by centroid.
  Steady laminar cases only for now: transient runs, turbulence, restart
  and passive scalars together with solids are refused (ENE006), as are a
  solid without a complete material (MAT006), solids without `[energy]`
  (REG001) and unknown zones (REG002). Check V7.4 (a composite wall with
  conductivity ratios 1 : 100 : 10⁴, exact to round-off) and a heated wall
  under a channel flow (gallery and recipe 156).
- **Requirement D204: cycle-averaged conjugate heat transfer for periodic flows** (engine exhaust manifolds), stage 7 page §9d: automated loop of transient fluid cycles with frozen solids, per-face cycle-averaged heat-transfer coefficient and reference temperature, steady solid solve, new wall temperatures until the tolerances are met. Plan: turbulent (7d5) and transient (7d6) CHT before buoyancy; checks V7.36–V7.40.
- **Turbulent conjugate heat transfer (stage 7d5).** Solid regions now
  work with every turbulence model. The fluid side of each fluid–solid
  interface is the thermal wall function of the flow (Kader or
  Jayatilleke, D135), h_w = ρ c_p u*/T⁺, in series with the solid's half
  cell, q = (T_P − T_N)/(1/h_w + δ_N/k_s); y* and μ of the first fluid
  cell are taken in every outer iteration, the fluid conducts with
  k + c_p μ_t/Pr_t, and μ(T) goes back to the laminar viscosity of the
  turbulence model. Check V7.36: a turbulent channel over a conducting
  wall gives the same heat flow as the fluid-only path with the wall
  function in series with a convection condition, to 10⁻⁴ for k-ε and
  SST (`cht_turbulent_channel.toml`). Gallery and recipe 158: a finned
  aluminium heat sink in a turbulent air channel
  (`cht_heat_sink_turbulent.toml`).
- **Fixes.** A `[[schedule]]` entry no longer drops the energy and scalar
  sub-tables of the boundaries (an energy case with an inlet and a
  schedule was refused). A temperature that turns NaN in a diverging run
  is reported as a divergence instead of stopping in the property
  interpolation.
- **Transient conjugate heat transfer (stage 7d6).** Solid regions now
  work with a `[time]` table. The fluid and the solids share one energy
  matrix with the time terms ρ c_p V/Δt of each cell's own material, so
  the interface needs no coupling iteration and has no time lag; the
  temperature's time levels follow those of the flow through adaptive,
  rejected and rewound steps. The heat balance gains the heat stored per
  region, and every transient run with the energy equation writes it per
  step to `heat.csv` (and to the run record), with the heat flow through
  each interface; the report gives the largest relative imbalance of a
  step. Probes give the temperature (`<probe>_T`), in a conjugate case on
  the whole mesh, solids included. Check V7.6: two bodies brought into
  contact keep the exact contact temperature (to 10⁻¹² K at 500 s) and
  follow the erf profiles at second order in space (2.02, 2.00) and time
  (2.05, 2.02) (`cht_contact_transient.toml`); check V7.37: the heated
  wall switched on closes its balance to 2·10⁻¹¹ in every step and ends
  at the steady solution (`cht_heated_wall_transient.toml`). Gallery and
  recipe 159.
- **`erf` and `erfc` in the expression language**, written in-house like
  the other functions (bitwise the same on every backend, within a few
  ulp of the C library), for exact solutions of transient conduction.
- **Sweeps.** A transient sweep over levels now checks the scalars with an
  exact solution (`scalar_max_error`, `scalar_order`) as a steady one
  does; a time sweep of a conjugate case takes the temperature on the
  whole mesh, and a fluid at rest no longer gives a NaN velocity order.
- **Contact resistance, orthotropic solids and region heat sources
  (stage 7d2).** `[interfaces.<name>]` puts a contact resistance
  (`contact_resistance`, m²K/W) or a thin wall (`thin_wall` with
  `thickness` and `conductivity` or `material`, the resistance t/k) into
  the faces between two regions, fluid–solid or solid–solid; the face
  conductance is the series of the two half cells and the resistance, and
  the report gives the heat flow and temperature of every pair of regions
  (`interface-<solid>`, and `interface-<a>-<b>` between two solids). Solid
  materials may now have an orthotropic conductivity (three principal
  values along two axes, each a form of T), oriented per region with
  `[regions.<name>] axes`: the matrix takes the conductivity of the face
  normal, the tangential rest of the flux is a deferred correction.
  `[[energy.sources]]` heats one region with `power` (W, spread by volume)
  or `power_density` (W/m³). Codes CHT004 and ENE008.
- **Exact interface flux on skewed meshes.** The cells next to a region
  interface take one-sided least-squares gradients (the face temperature
  of their own side in place of the cell across), and the interface flux
  uses the face-normal temperatures of both sides (D157 §9.2): a
  piecewise-linear temperature is exact on any mesh and at any
  conductivity ratio. Check: the patch test `cht_orthotropic_patch.toml`
  (a laminate with axes turned by 30°, a fluid at rest and a block behind
  a contact resistance, on a mesh perturbed by 25 %) is exact to
  3·10⁻¹¹ K; V7.4 with a contact resistance and a thin wall matches the
  series formula to 10⁻⁹. Gallery and recipe 160: a 1 W chip on a thermal
  pad on a circuit board with copper planes (`cht_circuit_board.toml`).
- **Heat balance on skewed meshes.** The reported conduction through a
  boundary with a fixed temperature now includes the non-orthogonal part
  of the flux, as the matrix has it, so the balance of the energy equation
  closes on skewed meshes too.
- **Generated meshes.** With `perturb`, nodes on a face plane of a `zones`
  box now move only along the plane, so the interfaces between the zones
  stay planes.
- **Stirred-tank images and the Fyfe 2 L tank.** The engineering-gallery tank (entry 131) now has velocity sections and an impeller picture. New entry 161 builds the 2 L lab tank of Fyfe et al. (OPRD 2024) from the paper's supporting information, runs it with the moving reference frame at 200 rpm and compares velocity planes and profiles with the paper's published CFD (power number about 2.0 against 1.8 from the paper's correlation). `tools.docs.images` also sizes images on the applications pages.
- **Rules and lessons in the repository.** New `AGENTS.md` holds the standing rules for agents and contributors as MUST and SHOULD, with the exact checks to run before pushing; new `docs/pitfalls.md` keeps the lessons of earlier work and the experiments not to repeat (steady-run budget and stopping, coupled-solver variants, skewed tetrahedra, observers, timing, iteration counts in tests) (D200). The pull request template asks what ran and what did not, and why (D203).
- **Stage 7 plan in the repository.** The working plan of stage 7 (pull request cut, lanes, decisions and status per pull request) is now `docs/requirements/stage7-plan.md`, linked from the stage 7 specification (D200).
- **One registry of tools (D196, part 1).** New crate `foehncraft-tools`: every engine function a front end calls becomes a typed tool with a stable id, a contract version, an action kind (read, write files, start a run) and input and output schemas from the same structs. `foehn commands` prints the catalogue as JSON (`docs/tools.md` is generated from it), `foehn call <id> --input file.json` runs any tool from JSON and prints a result envelope with status, exit code and diagnostics. In: the system, case-schema, case-check, case-resolve, diagnostic-code, material and expression tools; `foehn info`, `bench backends`, `schema`, `materials` and `flow --check` now go through them with unchanged text, JSON and exit codes. New commands `foehn case resolve`, `foehn diag codes` and `foehn expression`; new codes `TOOL001` and `TOOL002`.
- **Verification status page (D199, part 1).** Every feature a case file can switch on (91 rows: meshes, numerics, solvers, boundary conditions, turbulence models and wall treatments, heat transfer, energy, conjugate heat transfer, rotating machinery, runs, GPU and Python) has a grade, verified, partial or unverified, on the generated page `docs/verification/status.md`, with the evidence that checks it (cases with the bands of the verification manual, tests, recorded runs) and the known gaps. The grades are computed by `tools/release/verification_status.py` from the hand-written table `docs/verification/features.toml` and the manual, never typed in; CI fails when the table names a case, test or page that does not exist, leaves a catalogue case without a feature, or disagrees with the page. `docs/verification/grades.json` carries the same for the limitations of a run (part 2). First grades: 39 verified, 48 partial, 4 unverified (mass-flow inlet, static-pressure inlet, total-pressure outlet, buoyancy).
- **Iteration budget and diagnosis of steady runs (D198 (d)).** `foehn flow --check`, the log and the run record show the iteration budget and where it comes from. At `[run.limits] diagnose_at` (default 2000) a steady run that has not converged is diagnosed from its residuals: the trend of each equation over the last 500 iterations, the verdict (on track, slow, stalled, rising, oscillating) and advice, as a `diagnosis` row of `run.jsonl` and in the log; `diagnosis = "stop"` stops a run off track (exit code 3, `RUN017`). A budget above 10 000 iterations draws `RUN021` in the check. The diagnosis only reads the residuals: fields and residuals are bitwise those of a run without it.
- **Budgets per study (D198 (c)).** `[run.limits] study_cell_hours`, `study_gpu_hours` and `study_output` bound what all runs of a study spend together, kept in `study_budget.jsonl` next to the variant directories. The `budget` row of `run.jsonl` shows the D81 estimate (now with cell-hours and GPU-hours) beside each budget, the spend so far and what is left; a run that spends the rest stops (`RUN018`), the next one does not start (`RUN019`), an estimate that does not fit warns (`RUN020`). Gallery and recipe 166, case `cylinder_steady_diagnosis.toml`.
- **`unsafe` only in named modules (D202).** The workspace now denies `unsafe_code` and clippy's `undocumented_unsafe_blocks`. Twelve named modules allow unsafe code (the kernel bodies, the CUDA crate, the CPU and SELL-C backends, the host executor, the hypre interface and three tests of the kernel bodies); `python3 -m tools.check.unsafe_modules` lists them with their reasons and fails on `unsafe` anywhere else (CI). Six unsafe blocks got their missing `// SAFETY:` comment; the code is unchanged.
- **Broken-input corpus for every reader (D202).** Truncated files, damaged numbers, dropped and repeated lines, flipped bytes, garbage, the files of other formats and damaged CGNS nodes are fed to the Gmsh, SU2, industry .msh, VTK/VTU, OpenFOAM and CGNS readers, and damaged case files to the flow, mesh, scalar, Poisson and benchmark case readers. Each must end in an error with a diagnostic code, never a panic, an abort or a hang. The corpus found allocation aborts from counts a file declares (a mesh file declaring 10^14 elements aborted the program) and a read past the end of an OpenFOAM file; the readers now check declared counts against the file size and use checked arithmetic on ranges. Valid files read exactly as before.
- **New diagnostic codes `MESH003` and `CASE010`.** A mesh file that cannot be read (missing, unknown format, truncated, damaged) now ends in `MESH003` instead of `GEN001`; a case file that a reader refuses for a reason without a more specific code yet ends in `CASE010`. TOML and read errors of the mesh, scalar, Poisson and benchmark case files get `CASE001` to `CASE009` like the flow case files.
- **Runs and files as tools; the crate layering check (D196 part 2, D203).** New tools `case.estimate`, `case.upgrade`, `mesh.process`, `flow.run` (streams the run record and stops on request), `result.read` and `post.derive`; `foehn mesh`, `flow`, `case upgrade` and `post derive` go through them with unchanged text, JSON, exit codes and bitwise-equal results. Every function of the Python native module is now a wrapper of its tool, and `foehncraft.call(tool, **fields)` and `foehncraft.commands()` reach any tool from Python. `cargo xtask layering` runs in CI: crate dependencies must follow the layers, and the command line and the Python module call the engine only through the registry.
- **Requirements step 26: design optimisation.** Decisions D205–D212, L49 and the immersed-boundary method of L13 (cut cells as a mesh operation, Brinkman topology pre-stage): a differentiable chain of typed tools for outside optimisers, inexact one-shot steps, multi-fidelity trust regions with PyAeroFlow, gradient-enhanced Bayesian optimisation, surrogates from our own runs, optimisation as an agent tool with budgets, and zonal optimisation with adjoint-corrected boundaries ([step 26 page](docs/requirements/step26-design-optimization.md)).
- **Limitations of a run (D199, part 2).** Every result summary (`result.json`, `foehn flow --json`, Python `Result.limitations`) carries a `limitations` list with codes LIM001 to LIM006: a feature in use graded partial or unverified on the verification status page (the grades are compiled into the binary), a feature outside its validated range (first-cell y⁺, Reynolds, Prandtl and Courant numbers, skewness), a steady solve that stopped before its tolerance, one mesh only (no discretisation uncertainty), a GPU that is not verified (Tesla T4 so far), and automatic mesh remediation. The last frame of the log reads "Run finished with N limitations" and lists them, and `flow.md` starts with them. Status and exit code stay as they were (`finished`, 0); no flow result changes. Gallery 170 shows the grades per area and the raw Gmsh pipe with its three limitations; recipe 170 says what to do about each code.
- **Studies as tools; every command through the registry (D196 part 3).** New tools `study.scalar`, `study.poisson` and `study.g1`; `foehn scalar`, `poisson` and `g1` go through them with unchanged text, JSON, exit codes and results, so every command of `foehn` and every Python function now runs through the 19 tools of `docs/tools.md`. The command line depends only on `foehncraft-tools`, `foehncraft-flow` and the build crates, and `cargo xtask layering` rejects direct paths into the mesh, fv, bench, parallel and kernel crates from the front ends. Gallery 171 (the catalogue an agent reads), recipe 171 and the case page `docs/cases/tools.md` with JSON inputs in `examples/cases/tools/`.
- **Journal of tool calls and replay (D197).** `--journal FILE` (or `FOEHN_JOURNAL`, also read by the Python package) appends every tool call with its typed input (defaults filled in), status, exit code and duration to a JSON-lines journal (`foehn.journal`). A run writes the call that started it into `run.jsonl` (`{"type": "call"}` after the header; results unchanged). `foehn run --script FILE` (tool `journal.replay`) replays a journal or a run record without a window, skips calls refused when recorded, and stops at the first call that ends otherwise (`TOOL005`, exit 1) unless `--keep-going`. New codes TOOL003 to TOOL005. Gallery 172, recipe 172, `examples/cases/tools/cavity_session.jsonl`.
- **Requirements: layered distance-field immersed boundary (L13).** L13 target: wall on a thin layer mesh grown from a signed distance field at program start, outermost iso-surface as the smooth cut, imprint join; result of a 2D prototype (finned wall, narrow gap) on the [step 26 page](docs/requirements/step26-design-optimization.md#9-prototype-of-the-layered-distance-field-boundary-2d-2026-11-21).
- **Limitations in the tool-result envelope (D199, D196).** `foehn call` and `foehncraft.call` now copy the `limitations` of a tool's output (`flow.run`, `result.read`: LIM001 to LIM006 of the result summary) into the envelope's `limitations`, so an agent sees them without opening the result. Status, exit codes and results are unchanged.
- **MCP server `foehn mcp` (L18, D198 (a), (b), (e)).** An assistant (an MCP host) reaches every tool of the registry by its id over the Model Context Protocol on standard input and output, through the same call as `foehn call`, with progress notifications and cancellation (a cancelled run stops with a restart file, exit code 3). The server refuses to start without `--read-root` and `--write-root` (`MCP001`). An action table per project (`foehn-mcp.toml`, outside the write roots) sets each kind of action, start a run, exceed the estimate caps, write outside the write roots, overwrite results, publish, to default, allow, ask or refuse; "ask" puts the question to the person through the host (MCP elicitation), answered once or always. A run started through MCP needs `[run.limits] max_wall_time` and `max_output` (`MCP002`, with the fix). Every call is journaled with `"via": "mcp"` and the approval, so `foehn run --script` replays the session (refused calls skipped). New tool `result.page` pages the histories of a run (residuals, monitors, probes, statistics, the run record) by cursor; fields are named by path only. Tools tell what a call would write (`effects`), the approvals live in `foehncraft_tools::approval`. New codes TOOL006, MCP001 to MCP006. Guides [agents.md](docs/agents.md), [troubleshooting](docs/troubleshooting.md) and [install and update](docs/install.md) (D203). Gallery 175, recipe 175, `examples/cases/tools/mcp/`. Flow results unchanged.
- **Clean-room rework of the reference numbers (D201, audit E1 to E10).** No result of another solver is a reference any more. The plate of SST is checked against the Coles–Fernholz fit to measured layers (+1.45 % at Re_θ 7 848, band ±8 %) and White's correlation (−6.1 %, band widened from ±5 % to ±8 % for SST by the maintainer), with our own four-mesh Richardson value (2.7385e-3, GCI 1.76 %); `plate_channel_gci.py` and `tools/plot/sst_chart.py` compute the fit from the run's own profile and the TMR values are gone from the studies, cases, docs and the chart. The wall constants of `automatic` (n = 1.2, y₀ = 5), first fitted to the value of another solver, were rechecked against our own result on a 10-point grid and stay (largest deviation 2.5 %). The step no longer quotes a result of another code; the channel of the elliptic-blending model is judged against the DNS only (the test is `eb_channel_against_the_dns`); the "model in the thesis" column of the planned high-Pr tests and the Manceau entry of `DATA_SOURCES.md` are gone; DNS stays an accepted class (rule 5 of `DATA_SOURCES.md`, D201 and AGENTS.md say that D201 targets code-to-code RANS/CFD numbers). Entry 161 (Fyfe 2 L tank) no longer compares with the paper's CFD planes (script and figures removed); the bottom is rebuilt from the chamfers of the Supporting Information (read as one cone, flat to r = 7.5 mm, 15 mm high), the power-number reference is the measured literature band 1.0 to 1.7 (the earlier run gave 2.0: a finding), and the new `tracer.toml` with `mixing_time.py` sets up the mixing-time check against Table 3 of the paper (18.7 s, 200 rpm baffled) and the Grenville correlation; the runs wait for the maintainer's own mesh. New figure `fyfe_geometry.png`. Gallery 161 and recipe 161 reworked.
- **Live control port (L19, D198 (f)).** A run started with `foehn flow --control-port`, or any `flow.run` through `foehn mcp`, listens on a port of `127.0.0.1` only, with a random token in `control.token` (mode 0600 on Unix) beside `control.json` in its output folder; both disappear when the run ends. New tools `run.status` (the last iteration or step with its residuals, the wall time, what the run did with earlier changes), `run.set` (run-time settings, the keys a re-read may change, never `run.limits.*`) and `run.stop` (restart file, final output, exit code 3); new action kind `control`. Changes are applied after an iteration or step, checked against the case then (`RUN024` for one that does not fit), and written to the run record with the equivalent `[[schedule]]` entry (`RUN022` to `RUN024`, `TOOL007`). `foehn mcp` answers calls of the kinds read and control while a run is going, so an assistant can watch and stop a run it started; a journal replay skips live calls. Gallery 177, recipe 177, `examples/cases/tools/control/`, guide section in [agents.md](docs/agents.md). Flow results unchanged.
- **Requirements: staged plan for the immersed boundary method (L13).** Five stages (geometry kernel, penalisation pre-stage, cut cells, layered distance-field boundary, automatic hand-off to a body-fitted mesh) with verification cases, gallery cases, adjoint and CHT interaction, tool interfaces and open owner decisions; see the [L13 plan](docs/requirements/l13-immersed-boundary-plan.md). Plan only, no code.
- **Requirements step 27: staged plan for the adjoint solver (L17).** Proposal page with discrete against continuous adjoint (literature and computer-graphics benchmarks), a scatter-free reverse sweep on the D59 kernels, stages A0 to A6 (spike, steady laminar, geometry chain and morphing, turbulence, CHT, transient with checkpointing, extras), a gradient verification ladder, the tool and gradient interface, memory and cost relations, proposed entries D213 to D217 and L50, and ten open questions for the maintainer ([step 27 page](docs/requirements/step27-adjoint-solver.md)). Nothing is decided yet.
- **Thermal acceleration of solids in transient conjugate heat transfer (D213).** `[energy.acceleration]` scales ρ c_p of the solid regions by a factor (given, or `"auto"`: measured after the first window from each solid's conductance to the flow and its heat capacity, so that its time constant becomes one window, two periods in a periodic flow) until their temperatures have settled, then restores the true ρ c_p: at once in a steady flow, raised geometrically over `restore_ramp` periods in a periodic flow, where a jump would leave part of the accelerated swing in the deep solid. The steady and cycle-mean states are those of the true solids; the heated wall of gallery 159 gets there five times sooner, and with a periodic inlet the ramp keeps the cycle mean within 0.0025 K where a jump leaves 0.24 K. `heat.csv` books the heat the true solids would have stored as `inertia_scaling`; the report gives the measured conductance, time constant and factor, the restore times and from when the history counts as physical. Codes ACC001–ACC004; checks V7.41, V7.42; gallery and recipe 180.
- **An assistant runs a heat sink study through `foehn mcp`.** The fin-pitch sweep of the finned heat sink (four pitches, 24 to 48 mm) was run by an assistant under run limits and a study budget, with the refusals of the action table shown (no run limits `MCP002`, a folder outside the write root and a mesh over the cell cap `MCP003`), live status of a running case, histories read by page and a closing report with its limitations. The recorded journal replays with `foehn run --script`; a replay gives byte-identical residual histories. New: `tools/agent/stdio_host.py`, a small host that keeps `foehn mcp` open behind a socket for scripts and for assistants without MCP registration; figure tools `tools/plot/heat_sink_study.py` (the fields and the numbers) and `tools/plot/heat_sink_agent.py` (the session timeline with its refusals and budget, the speed fields and the temperature the assistant read). Gallery, recipe, `examples/cases/tools/mcp/heat_sink_study/`, the transcript in [agents.md](docs/agents.md). No engine change.
- **Requirements: log redesign review (D214).** The console log redesign is built (D125 to D128); the maintainer added `[log] detail`, a schedulable `log.columns`, a Limitations block at the start and an Events block at the end; see the [log review](docs/requirements/log-redesign-review.md). Plan only, no code.
- **Requirements: adjoint plan decisions (L17, step 27).** The maintainer's answers to the ten questions: discrete adjoint on the converged residual (D220), scatter-free reverse sweep from local dual Jacobians with a time-boxed spike on a compiler reverse-mode tool (D221), a gradient verification ladder including a whole-solver complex step, which makes the flow solver generic over `Real` (D222, stage A0b), adjoint tools and gradient contract (D223), non-chaotic transient adjoint only (D224), chaotic-flow sensitivities as research item L51; L17 staged A0 to A6, first demonstration a laminar bend pressure drop ([step 27 page](docs/requirements/step27-adjoint-solver.md)).
- **Manufactured solutions of k-ε-ζ-f₀, the curvature correction and the
  QCR (stage 6 specification 9.2).** The harness of PR #148 now covers
  ζ-f₀ (orders 1.85–2.27 on hexahedra, polyhedra, tetrahedra), the
  curvature correction of SST in a second prescribed flow whose strain axes
  turn (1.9–2.1 on hexahedra and polyhedra) and the explicit stress of the
  QCR2000/QCR2013 (order 2 on regular hexahedra and polyhedra). Findings,
  not fixed here: the curvature correction does not converge on tetrahedra
  (DS_ij/Dt from the least-squares gradient of the gradient), and the
  divergence of the explicit QCR stress is first order on tetrahedra and
  does not converge on randomly perturbed meshes. New:
  `Flow::explicit_stress_force`. The grades of the three features list
  the gaps. Gallery and recipe 184.
- **Requirements: L13 plan decisions.** The maintainer's answers to the eleven open questions are recorded in the [L13 plan](docs/requirements/l13-immersed-boundary-plan.md#9-owner-decisions-2026-11-22): tolerant dirty STL with Crane's methods as candidates, merged small cells, own near-body refinement, a 3D Python prototype before Rust, and a new stage 5 for prescribed rigid motion.
- **Requirements: ids D225 to D234 for the immersed boundary stages.** Geometry kernel, penalisation pre-stage, cut cells, layered boundary, hand-off and rigid motion are entered in the requirements table.
Nothing yet.

## Summary

| Version | Stage | Highlights |
|---|---|---|
| 1.0.0 (2026-11-23) | 8 (part) Milestone 1 | The internal-flow cases of D45 (pipe, plane diffuser, U-bend, T-junction on polyhedra) with grid convergence, several fluids in one conjugate case (D194), gravity, Boussinesq buoyancy and the low-Mach variable density (7b1, 7c), absolute pressure values and per-region start values (D239, D238; case format 0.7 stable), the cycle-averaged loop and the frozen-flow march for slow thermal transients (D204, D118), the operator exactness check with quadratic fits (D236), the MLP vertex limiter (D237), exact CPU reductions (MPI W1, W2). **Pipe u_τ +4.5 % from DNS (GCI 0.2 %); diffuser reattachment GCI 0.36 %; U-bend separation GCI 0.19 %; stratified box at rest to 10⁻¹⁵ m/s; closed box p₀(t) exact; warm-ups in a tenth of the steps within 0.06 K.** Gallery 186–221, recipes up to 221 |
| 0.7 (2026-11-22) | 7 (part) Energy, conjugate heat transfer, agents | Energy equation of the fluid with temperature-dependent properties, conjugate heat transfer with solid regions (steady and transient, laminar and turbulent, contact resistances, orthotropic solids, heat sources, thermal acceleration of solids), the tool registry, journal and replay, the MCP server `foehn mcp`, run limits and budgets, the live control port, run limitations from the verification status page, manufactured solutions of every turbulence model, the clean-room rework of the reference numbers. **CHT checks exact to round-off (V7.4, V7.6 to 10⁻¹² K, V7.36 to 10⁻⁴); heat balance closes to 10⁻¹¹ in every step; thermal acceleration ends within 0.2 % of the steady solution.** Gallery 131–184, recipes up to 184 |
| 0.6 (2026-11-20) | 6 RANS turbulence, heat transfer, rotating machinery | Six RANS models: k-ω SST 2003 (with the curvature and rotation correction and QCR2000/QCR2013), standard and realizable k-ε, k-ε-ζ-f, k-ε-ζ-f₀ and the elliptic-blending k-ε; wall treatments `low-re`, `automatic`, `scalable`, `standard`, `compound`, `compound-2023` and rough walls; passive scalars with turbulent diffusivity and thermal wall functions (also for high Prandtl numbers); periodic boundaries; rotating machinery with the moving reference frame, sliding meshes (conformal stepping and imprinted interfaces of any shape of revolution, several rotors), torque, shaft power and machine performance; the coupled pressure–velocity algorithm with pseudo-time steps and multigrid for momentum and turbulence; SU2, OpenFOAM, CGNS, industry `.msh` and VTK mesh readers; derived turbulence fields, transient statistics, the grid convergence study of the models; the new log, case format 0.6 and `[[schedule]]`. **TMR flat plate C_f +0.37 % (SST) and within −2.6 to +2.4 % for first-cell y⁺ 0.24 to 85; channel at Re_τ 395 U_b⁺ −0.7 % (SST) and EB C_f 97.8–99.1 % of DNS; backward-facing step x_r/H +0.5 % (k-ε against the published k-ε) and +0.6 % (k-ε-ζ-f₀ against the experiment); impinging jet Nu₀ −0.2 % (SST); heated channel Nu within −2.6 to +1.0 % of DNS; GPU bitwise equal to the CPU for 14 turbulence cases on a Tesla T4; the coupled algorithm settles the step in about 70 instead of 600 iterations.** Gallery 26–134, recipes 36–134 |
| 0.5 (2026-11-16) | 5 Transient flow, run control, Python | Transient incompressible flow: PISO with variable-step BDF2 and Crank–Nicolson, outer correctors for large time steps, a time step from the integration error, output at exact times, restart files and control files; the whole flow loop on NVIDIA GPUs; adjoint-ready kernels (dual numbers, transposed matrices, checkpointing); the run contract for scripts and agents (`--json`, `--set`, `result.json`, run record, exit codes and signals, case schema, diagnostic codes, estimate, limits, provenance); the conservation ledger and the divergence report; the bounded central scheme and the boundary conditions of D29 (slip, wedge, tables, backflow, convective outlet); derived fields; expressions, the material database and regions; the Python package `foehncraft` with studies; the grid convergence index of every sweep and a generated verification manual; all scripts in Python. **Taylor–Green time order 2.0–2.4 and space order 2.1; Schäfer & Turek 2D-1 inside all reference ranges, C_D extrapolated to 5.5797 ± 0.0004 from three meshes; GPU flow loop bitwise equal to the CPU on a Tesla T4 (722 ns per cell and iteration on 970 225 polyhedra); a case built in Python is bitwise its TOML file.** Commands `foehn materials`, `foehn post derive`, `foehn schema`; gallery 14–25, recipes 21–35 |
| 0.4 (2026-11-15) | 4 Laminar steady flow | First flow solver: steady incompressible laminar flow with SIMPLEC (and SIMPLE) on polyhedral meshes, relaxation-independent Rhie–Chow interpolation (Majumdar 1988), CG + AMG for the pressure correction, boundary conditions of D29 (walls also moving, velocity and mass-flow inlets, static and total pressure inlets and outlets, symmetry) with a check of their combination and guidance (D30), TOML case file with "did you mean" (D39), log, CSV, wall and line output (D43). **Couette exact to the iteration tolerance on all mesh types; Poiseuille observed order 1.80–2.00; lid-driven cavity within 0.4–1.4 % of Ghia et al. at Re 100 and 1000 in 2D and within 2.5 % on 3D tetrahedra and polyhedra.** New crate `cfdcraft-flow`, command `cfdcraft flow`, thirteen case files, gallery 12–13, recipes 18–20. CPU only; the GPU flow loop follows with stage 5 |
| 0.3 (2026-11-14) | 3 Pressure solver | Pressure linear solver: smoothed-aggregation algebraic multigrid as the preconditioner of CG (Vaněk, Mandel, Brezina 1996), symmetric V-cycle with Chebyshev or Jacobi smoothing, float32 hierarchy with float64 Krylov iteration (mixed precision, D37), one kernel source for the CPU and the NVIDIA GPU, optional hypre 3.2.1 reference (BoomerAMG-PCG, built from source). **10–13 iterations from 512 to 274 625 cells on all mesh types; same solution as hypre to 1.6e-9; on a Tesla T4 6.5 ns per cell and iteration against 10.5 for hypre at 912 673 cells, all GPU results bitwise equal to the CPU.** New command `cfdcraft poisson`, four Poisson case files, gallery 10–11, recipes 15–17 |
| 0.2 (2026-11-14) | 2 Operators on a scalar equation | Finite-volume operators on a scalar transport equation, one source for the CPU and the NVIDIA GPU: least-squares gradients (exact for linear fields, error 6e-15 on all mesh types), linear upwind with Venkatakrishnan and Venkatakrishnan–Wang limiters (no new extrema at K = 0), over-relaxed non-orthogonal diffusion, reproducible gather-based assembly into sliced ELL, BiCGStab. Observed order 1.87–2.08 on hex, tet, poly and mixed meshes in 2D and 3D. **GPU and CPU results bitwise equal on a Tesla T4** (912 673 polyhedra, 311 ns per cell and iteration). New command `cfdcraft scalar`, seven case files, gallery 8–9, recipes 12–14 |
| 0.1 (2026-11-14) | 1 Mesh | Polyhedral face-based mesh with exact geometry (every cell closed to 1e-15, volumes exact), Gmsh 4.1 reader, built-in test mesher (hex, tet, poly, mixed; 2D and 3D; graded, refined, perturbed), quality report, VTKHDF output for ParaView, METIS and RCB partitioning with a solver-ready halo layout and multi-field halo exchange (tested over MPI on 2–4 ranks), particle point location on all mesh types, golden-value test matrix, `cfdcraft mesh` command. G1 rerun on polyhedral and shuffled meshes: bitwise equal; quiet-machine verdict open (O15) |
| 0.0 (2026-11-14) | 0 Foundations and gate G1 | Cargo workspace, CI on Linux, macOS and Windows with MPI runs on 4 ranks, the G1 kernel benchmark. **G1 GPU verdict passed on a Tesla T4: single-source Rust kernels within 1.0 % of hand-written CUDA C, bitwise equal results**; **CPU verdict passed on an i9-9900K: single-source bodies within 1.8 % of hand-written loops, in two runs** |

Release 1.0.0 is milestone 1: incompressible and low-Mach flow, steady
and transient, laminar and RANS, with the energy equation, conjugate heat
transfer with several fluids and solids, buoyancy and rotating machinery,
on CPUs and NVIDIA GPUs in one process, verified case by case on the
[verification status](docs/verification/status.md) page. Multi-rank MPI,
natural-convection benchmarks and the rest of stage 7 follow in 1.1
(requirements section 7a). The Cargo version of the crates and of the
Python package is now 1.0.0; the release is the Git tag `v1.0.0` (D46).

Release 0.7 is the first part of stage 7: the energy equation of the
fluid, conjugate heat transfer with solid regions and the tools for agents
(registry, MCP server, live control port). Buoyancy and variable density
follow in later releases. The Cargo version of the crates is now 0.7.0;
the release is the Git tag `v0.7` (D46).

Release 0.6 closes stage 6: all six turbulence models, the wall
treatments, passive scalars with thermal wall functions and rotating
machinery are verified on the cases of the stage plan (TMR flat plate,
channel against DNS, heated channel and pipe, backward-facing step of
Driver & Seegmiller, impinging jet of the ERCOFTAC database,
Taylor–Couette flow, stirred vessels), with the findings and the open
points listed under "Known limits" of the 0.6 section and in the
[stage 6 verification notes](docs/verification/stage6.md). The generated
[verification manual](docs/verification/manual.md) lists every case with
its pass bands. The Cargo version of the crates is now 0.6.0; the release
is the Git tag `v0.6` (D46).

Release 0.5 closes stage 5: the acceptance checks of the stage plan pass
(Taylor–Green time and space order 2, also at Courant number 9 with outer
correctors; the cylinder benchmark of Schäfer & Turek, at Re 20 inside all
reference ranges, at Re 100 with C_L,max 0.5 % below its range on 24 576
cells; restart bit-identical; expressions from TOML and from Python
give the same run; the energy budget of the cavity closes; the GPU flow
loop gives bitwise the CPU results). The generated
[verification manual](docs/verification/manual.md) lists every case with
its pass bands. The Cargo version of the crates is now 0.5.0; the release
is the Git tag `v0.5` (D46).

Release 0.4 closes stage 4: all acceptance tests of the stage plan pass
(Couette at round-off, that is at the iteration tolerance; Poiseuille order
1.8–2.2 on all mesh types; lid-driven cavity in 2D and 3D at Re 100 and 1000
against Ghia et al.). The Cargo version of the crates is now 0.4.0; the
release is the Git tag `v0.4` (D46).

Release 0.3 closes stage 3: all acceptance tests of the stage plan pass
(iteration count nearly independent of the mesh size, same converged result
as hypre, GPU time per cell compared with hypre and reported). The Cargo
version of the crates is now 0.3.0; the release is the Git tag `v0.3` (D46).

Release 0.2 closes stage 2: all acceptance tests of the stage plan pass
(exact gradients, no new extrema, observed order 1.8–2.2 on all mesh types,
bitwise reproducible results, GPU and CPU bitwise equal on a Tesla T4). The
Cargo version of the crates is now 0.2.0; the release is the Git tag `v0.2`
(D46).

Release 0.1 closes stage 1: all acceptance tests of the stage plan pass
(closed cells, exact volumes, Gmsh files from Gmsh load, VTKHDF files load
with VTK's reader, point location on all mesh types). The Cargo version of
the crates is now 0.1.0; the release is the Git tag `v0.1` (D46).

Release 0.0 closes stage 0: both verdicts of gate G1 are in (GPU: Tesla T4,
CPU: i9-9900K, 2026-11-14, both PASS; O13 closed). The Cargo version of the
crates stays 0.0.0; the release is the Git tag `v0.0` (D46).

---

## 0.6 (2026-11-20) – stage 6: RANS turbulence, heat transfer and rotating machinery

**Highlights.**

- **RANS turbulence, six models.** k-ω SST 2003 (Menter, Kuntz & Langtry)
  with production limiters, sustaining terms and rough walls; standard and
  realizable k-ε; k-ε-ζ-f (Hanjalić, Popovac & Hadžiabdić 2004);
  k-ε-ζ-f₀ (Popovac 2023); the elliptic-blending k-ε (Billard & Laurence
  2012). One kernel source for the CPU and the GPU, the exact or the
  Poisson wall distance, inlet turbulence from intensity and length scale
  or viscosity ratio.
- **Wall treatments.** `low-re`, the SST `automatic` treatment that holds
  C_f within −2.6 to +2.4 % for first-cell y⁺ from 0.24 to 85,
  `scalable`, `standard`, `compound` and `compound-2023` wall functions,
  rough walls (Aupoix for SST, Cebeci & Bradshaw for the k-ε family).
- **SST options.** The curvature and rotation correction (Smirnov & Menter
  2009), a whole domain in a turning frame (`flow.frame_rotation`), and the
  quadratic constitutive relation QCR2000/QCR2013 for the secondary flow of
  non-circular ducts.
- **Heat and mass transfer in turbulent flow.** Passive scalars with
  ν_t/Pr_t diffusion, the Kader and Jayatilleke thermal wall laws, the
  Kays–Crawford turbulent Prandtl number and a high-Prandtl variant with
  variable properties (P-function); periodic boundaries with a heat-flux
  jump for developed heated channels and pipes.
- **Rotating machinery.** The moving reference frame (steady and
  transient), the sliding mesh by conformal stepping and by imprinting the
  interface at every step, several rotors at their own speeds, interfaces
  of any shape of revolution (cylinders, discs, cones), every turbulence
  model on all of them, speeds as `omega` or `rpm`; torque, shaft power,
  the angular-momentum balance and machine performance (head,
  efficiencies, Euler work); relative and blade velocities as derived
  fields.
- **Faster steady runs.** The coupled pressure–velocity algorithm
  (matrix-free, flexible GMRES with a SIMPLE-type preconditioner) with
  pseudo-time steps and a reduction multigrid for momentum and the
  turbulence equations: the backward-facing step settles in about 70
  outer iterations instead of about 600 with SIMPLEC.
- **Meshes from other tools.** SU2, OpenFOAM polyMesh, CGNS (HDF5), the
  face-based industry `.msh` and VTK (`.vtk`, `.vtu`) readers behind one
  mesh-format interface, besides Gmsh, and meshes written as `.vtu`; removed boxes in
  the built-in mesher (a backward-facing step) and annuli for rotors.
- **Results and verification.** Derived turbulence fields (time and length
  scales, Reynolds stresses, wall fields with the viscous weight), time
  means, RMS and covariances of transient runs, the grid convergence study
  of all models on the plate and the channel, the CPU–GPU parity run of the
  turbulence models (14 cases bitwise equal on a Tesla T4), and the
  published sources of all reference data
  (`docs/verification/DATA_SOURCES.md`).
- **Run experience.** The rewritten log with framed blocks and a
  configurable table, case format 0.6 with `foehn case upgrade`, solver
  settings per equation, settings that change during a run
  (`[[schedule]]`), automatic remediation of skewed tetrahedra (D50).
- **Requirements for what comes next.** Stage 7 (energy, buoyancy,
  conjugate heat transfer, thin walls and the heat of the wall film),
  data-driven models (D165–D169), LES and high-Prandtl flows, the
  elliptic-blending Reynolds-stress model (L47) and an applications track.

**Validation in numbers.** TMR flat plate C_f(0.97) +0.37 % from the TMR
value with SST (GCI 1.87 %), −2.6 to −0.3 % from White with scalable k-ε
for first-cell y⁺ 1 to 100, +0.5 % with EB; channel at Re_τ 395 U_b⁺
−0.7 % from the DNS with SST, EB C_f 97.8 / 99.0 / 99.1 % of the DNS at
Re_τ 180 / 395 / 590; backward-facing step x_r/H 5.372 with k-ε (+0.5 %
from the published k-ε, GCI 1.94 %), 6.296 with k-ε-ζ-f₀ (+0.6 % from the
experiment), 6.747 with SST (+7.8 %); impinging jet Nu₀ −0.2 % (SST) and
+7.2 % (k-ε-ζ-f₀) from the experiment; heated channel Nu within −2.6 to
+1.0 % of the DNS at Pr 0.71 for Re_τ 180 to 1020; heated pipe Nu 182.6
against Gnielinski's 180.2; Taylor–Couette flow in a moving frame at
observed order 2.3–2.4.

**Known limits.**

- **k-ε-ζ-f₀ near the wall.** Below first-cell y⁺ 1 the channel is not
  grid-convergent: U_b⁺ drifts up with refinement (17.79, 17.84, 18.11,
  18.39 on levels 0 to 3 against the DNS 17.55), because the k sink of the
  model is linear in k at the wall. A consistent sink is being worked out
  and will come in a later release; meshes with first-cell y⁺ of about 1
  give the values above.
- **k-ε-ζ-f** does not converge on the finest channel mesh of the GCI study
  in 40 000 iterations, and its plate boundary layer stays laminar from the
  TMR free stream.
- **ζ-f and EB on the step and the jet.** The reattachment length of ζ-f
  (7.74) and EB (7.56) is outside its band on the middle mesh (their
  finest mesh is still to run); Nu₀ of the jet is +22.9 % (ζ-f) and +76 %
  (EB).
- **SST plate order.** The observed order of C_f on the TMR plate is 0.43,
  below the band of the specification (the GCI is 1.87 %).
- **Realizable k-ε** gives 3 to 7 % low friction; the compound treatment
  overpredicts the wall shear below y⁺ 30.
- **Open checks that need the maintainer or data:** the GPU run of the two
  SST options (curvature correction, QCR) on the Colab T4; the GCI levels 1
  to 3 of the corrected study on the maintainer's machine; two DNS values
  (the rotating channel of Kristoffersen & Andersson 1993, the square duct
  of Pinelli et al. 2010); the NASA wall-mounted hump data (host not
  reachable from the build machine); manufactured solutions per
  turbulence model.

The parts in detail, in the order of their pull requests (the first
entries, before the fragments of `changelog.d/`, are at the end):

- **Pull requests no longer conflict on bookkeeping.** Changelog entries
  are fragments in `changelog.d/` that the release collects
  (`tools/release/changelog.py`); new gallery and recipe entries take the
  number of their pull request; the number tables of `docs/gallery.md` and
  `docs/recipes.md` are regenerated after each merge by the workflow
  `docs-index.yml` (`tools/docs/indexes.py`); topic pages merge with git's
  union driver. CI checks the lock file of the Python crate first. Rules in
  CONTRIBUTING.md, "Working in parallel".
- **Stage 6 specified (RANS turbulence and wall treatments).**
  `docs/requirements/stage6-turbulence.md` with decisions D129–D139:
  k-ω SST 2003 with Kato–Launder and Durbin production options, standard
  and realizable k-ε, k-ε-ζ-f, all-y⁺ momentum wall treatments, the
  thermal law of the wall (Kader, Jayatilleke) with constant or
  Kays–Crawford turbulent Prandtl number, turbulence numerics, checks,
  case keys, and a verification suite with GCI pass bands (flat plate,
  channel at Re_τ 180/395/590, backward-facing step, impinging jet). The
  elliptic-blending k-ε is open (O17).
- **Stage 7 specified: energy equation, buoyancy, conjugate heat transfer.**
  Decisions D150–D164 with the full text in
  `docs/requirements/stage7-energy-buoyancy-cht.md`: enthalpy conserved and
  the temperature correction solved, low-Mach variable density, SGDH and
  GGDH, Boussinesq and variable-density buoyancy with a balanced-force
  interpolation, fluid–solid interfaces in one temperature matrix, 25
  verification and validation cases with pass bands. The maintainer
  confirmed the new stage order (O18): energy as stage 7, v1.0 after
  stage 8.
- **Data-driven models assessed (step 20, open step O19).**
  `docs/requirements/step20-data-driven-models.md` reviews the literature of
  2025–2026 on learned RANS closure corrections, learned linear-solver
  components, neural surrogates and flow foundation models, judges them for
  industrial accuracy and proposes, not yet decided, D165–D169 (acceptance
  rules for learned components, invariant features in expressions, ensemble
  inversion, data and model licences, learned components that change only the
  path), L37–L38 and P13–P15. Learned models stay outside the solver and come
  in only as fields or compiled expressions.
- **Passive scalars (stage 6 part 6a).** `[scalars.<name>]` adds a scalar
  carried by the flow (D51), solved with the generic scalar transport
  equation of stage 6 (D116) in its bounded form, steady and transient,
  with single-source kernels for the CPU and the GPU. Boundary conditions per scalar in
  `[boundary.<name>.scalar.<scalar>]` (`value` or `flux`, `SCAL001`),
  sources, an `exact` solution with errors and observed orders, `bounds`
  with the bounded BDF2 (D28 (1)). Scalars appear in `residuals.csv`,
  `history.csv`, the console table (`residual.scalars`), the fields, the
  lines and the restart files. The laminar heated channel gives
  Nu = 8.2398 against 140/17 = 8.2353 and second order (gallery 63,
  recipe 63).
- **Exact wall distance.** A bounding volume hierarchy over the wall faces
  finds the nearest wall of every cell (D48, D108 part 1), deterministic
  for any number of threads; the derived field `wall-distance` writes it.
- **Smooth limits.** Hard and C¹ min, max and abs in
  `foehncraft_core::math` for the turbulence limiters (spec 2.6).
- **Data-driven models decided (O19 closed).** Adopted D165 (acceptance
  rules for learned components: do no harm on the stage-6 suite, improvement
  beyond the validation uncertainty, held-out cases, grid convergence with
  the model on, model card), D166 (invariant features in expressions, stage 6
  pull request 6b), D167 (ensemble inversion as a study type), D168 (data and
  model licences), L37 (closure campaign for SST and ζ-f on the diffuser and
  the U-bend, after v1.0) and P13 (learned smoother as a watch item).
- **Stage 6 specification extended.** The curvature and rotation
  correction of SST (Smirnov & Menter 2009) and the quadratic constitutive
  relation QCR2000/QCR2013 as options (D170, pull request 6h), the explicit
  algebraic Reynolds-stress model of Wallin & Johansson with the k-ω of
  Hellsten as a later item (D171), and the k-ω of Wilcox (2006) as an open
  question.
- **Learned components may change only the path (D169).** A learned solver
  component or a surrogate start reaches the same converged answer with or
  without it; CI checks this in reproducible mode.
- **Periodic boundaries (stage 6 part 6p).** Two boundaries of
  `type = "periodic"` that name each other as `partner` form a
  translational pair (D29): the faces are matched by their centroids, each
  side keeps its own faces, and every face loop, gradient and limiter sees
  the cell behind the partner face shifted by the period. The flow through
  the pair is driven by its `mass_flow` (a uniform body force adjusted
  after every pressure correction) or by a fixed `pressure_gradient`; the
  driving force is in the log, the report and the restart file. Checks
  BC030 (pairs that do not match or name each other) and BC031 (a pair that
  nothing drives). Generated boxes keep their perturbed periodic sides
  translated copies. Results of cases without periodic boundaries are
  unchanged bit for bit.
- **Heat-flux jump of the passive scalars (D109 (1)).**
  `periodic = "heat-flux"` in a scalar table of a periodic boundary carries
  the scalar across the pair with the jump Δ = Q/ṁ of the wall fluxes and
  sources, `bulk` fixes its level; refused with fixed wall or inlet values
  (D109 (4)). The periodic heated channel gives Nu = 8.2460 against
  140/17 = 8.2353 and the short periodic heated pipe, one diameter long,
  Nu = 4.3691 against 48/11 = 4.3636 (D111), both at second order; the
  periodic Poiseuille flow equals the open channel (gallery 67, recipe 67).
- **k-ω SST 2003 (stage 6 part 6b, D130).** `[turbulence] model =
  "k-omega-sst"` adds the SST model of Menter, Kuntz & Langtry (2003) on
  the scalar path of 6a: point-implicit (Newton) destruction, the
  production limiter (`menter`, `kato-launder`, `durbin`, `none`), the
  sustaining terms (`ambient`), the ν_t cap, constant overrides by name
  (TURB007), smooth limits for the adjoint. Wall treatments `low-re`
  (k = 0 and ω = 60ν/(β₁ d₁²) on the wall face) and `automatic` (blended
  viscous and log-law u_τ, ω and production in the first cell, fixed with
  a fractional weight in the row assembly); rough walls after Aupoix
  (`roughness`, BC027). Inlet turbulence from an intensity with a length
  scale or a viscosity ratio, or k and ω directly (BC028, BC029, TURB009).
  The wall distance is exact (6a) or the Poisson distance
  (`turbulence.wall_distance = "pde"`). Log and report: the first-cell y⁺
  per wall, the largest ν_t/ν, floor and cap counts (TURB001, TURB002,
  TURB008); `run.jsonl` rows `turbulence` and `turbulence_floors`.
- **Derived field `turbulence-features` (D166).** The invariant features
  of data-driven corrections (Pope's five invariants of the normalised
  strain and rotation, wall Reynolds number, time-scale and viscosity
  ratios, F₁, Q, the streamwise pressure gradient, the curvature of the
  strain) as datasets of the field files.
- **Verification of SST (gallery 68, recipe 68).** TMR flat plate:
  C_f(0.97) = 2.7007e-3 on 256 × 192 cells, +0.37 % from the TMR value,
  C_D +0.41 %; observed order 0.43 (GCI 1.87 %), below the band of the
  specification, recorded as a finding. Channel at Re_τ = 395:
  U_b⁺ = 17.43, −0.7 % from the DNS of Moser, Kim & Mansour. The
  automatic treatment holds C_f within 2.3 % for first-cell y⁺ of 1.7, 27
  and 92 and overpredicts it by up to 16 % in the buffer layer
  (TURB001 warns there). Decay of k and ω at second order in time; restart
  and 1 against 4 threads bit-identical. Curvature correction and QCR
  (D170) follow in 6h.
- **LES and high-Prandtl assessment (step 21).** Review of the
  scale-resolving and high-Pr prompt against the code and the literature
  (`docs/requirements/step21-les-high-prandtl.md`). Adopted D172 (turbulent
  Prandtl number as a closure field), D173 (variable properties, semi-local
  scaling and an integrated thermal wall function, stage 7), D174 (per-face
  convection blending), D175 (equation-free turbulence-model slot), L39–L43
  (LES numerics, SGS models, wall-modelled LES, divergence-free synthetic
  inlet, LES statistics and a campaign test tier, after v1.0), P16–P17 and
  O20. PANS is dropped altogether (X2, was P7).
- **High-Prandtl wall treatment from the TU Graz P-function paper.** D172
  gains the `kays-crawford-hp` turbulent Prandtl number (modified
  Kays–Crawford model of Irrenfried & Steiner with a Pr-dependent wall
  value), and D173's `"integrated"` thermal wall function is specified as
  the revised P-function (full van Driest mixing length with that Pr_t),
  with reference values checked against the paper's figures.
- **Standard and realizable k-ε (stage 6 part 6c, D131).** `[turbulence]
  model = "k-epsilon"` (Launder & Spalding 1974) or
  `"realizable-k-epsilon"` (Shih et al. 1995, variable C_μ and the
  destruction C₂ε²/(k + √(νε))) in one kernel. The standard model is
  bounded by Durbin's realizability limit by default (`production_limiter`
  also `kato-launder`, `menter` as min(P_k, 10ε), `none`); the ε
  destruction is linearised in the Picard form by default. ε at inlets and
  in `[initial]`, or ω converted with ε = C_μ k ω; sustaining terms with
  `ambient = { k, epsilon }`; a restart file of SST continues with k-ε and
  back. Log, residual and output columns name ε.
- **Wall functions (D134).** `wall_treatment = "scalable"` (default of the
  k-ε family, any y⁺, accurate from y⁺ 30), `"standard"` (30 ≤ y⁺ ≤ 300)
  and `"compound"` (Popovac & Hanjalić blend; with the high-Re k-ε claimed
  from y⁺ 30, owner decision: resolved walls go to ζ-f and EB). Rough walls
  shift the log law by ΔB of Cebeci & Bradshaw (`roughness`,
  `roughness_constant`, default 0.5). TURB001 bands per treatment.
- **Verification (gallery 71, recipe 71).** TMR flat plate with the
  scalable treatment: C_f(0.97) within −2.6 to −0.3 % of White for every
  first-cell y⁺ from 1 to 100. Channel at Re_τ = 395: U_b⁺ +4.3 % from the
  DNS at y⁺ 32. Rough pipe and the Moody study
  (`examples/studies/moody_study.py`): smooth pipes −2 to −5 % from
  Colebrook, fully rough +8 to +11 % with C_s = 0.5 (+1 to +3 % with 0.4).
  Findings: the realizable model gives 3 to 7 % low friction; the compound
  treatment overpredicts the wall shear below y⁺ 30.
- **High-Prandtl wall treatment with variable properties.** D172's
  `kays-crawford-hp` and D173's `"integrated"` thermal wall function now
  default to the variable-property P-function model of Steiner & Irrenfried
  (2019). Our 1D pipe check is within 3.1 % of the published DNS Nusselt
  numbers. The CC BY 4.0 thesis of Irrenfried (2020) supplies the reference
  data for the high-Prandtl checks, which answers O20.
- **k-ε-ζ-f model (stage 6 part 6d, D132).** `[turbulence] model =
  "k-epsilon-zeta-f"` (Hanjalić, Popovac & Hadžiabdić 2004): ζ = v²/k and
  the elliptic relaxation f added to k and ε, ν_t = C_μ ζ k T with the
  realizability and Kolmogorov bounds of T and L. Wall treatments
  `compound` (default, any y⁺: k and ζ by zero gradient, ε, P_k and f_w
  blended between sublayer and log law) and `low-re` (k = ζ = 0,
  ε_w = 2νk₁/d₁², f_w = −2νζ₁/d₁²). `zeta` at inlets (BC028 with another
  model) and in `[initial]`; ζ, f, T and v² in the output, the restart
  file and the log columns; a k-ε restart file continues with ζ = 2/3,
  f = 0. A solve lowers k, ε or ζ by at most a factor 10 per outer
  iteration (positivity safeguard).
- **Multigrid of the turbulence equations (D139).** Hierarchy slots on
  both backends (pressure, f, k, ε or ω, ζ). The elliptic f equation runs
  CG with the multigrid by default (`[solver.elliptic]`); the transport
  equations take it with `solver.turbulence.preconditioner = "amg"`,
  built from the symmetric part ½(A + Aᵀ). All hierarchies rebuild on the
  schedule of the pressure (`solver.amg_rebuild`), and their matrices go
  into the restart file (bit-identical continuation). The run-time table
  shows the turbulence model and the passive scalars.
- **Verification (gallery 73, recipe 73).** Channel at Re_τ = 395 resolved
  to y⁺ 0.27: U_b⁺ +1.9 % (`compound`) and +1.4 % (`low-re`) from the DNS.
  TMR flat plate at y⁺ 1: C_f(0.97) +5.1 % from White (finding, outside
  ±3 %). Findings: from the TMR free stream the boundary layer stays
  laminar (the case starts turbulent); near the leading edge the floor of
  k acts (TURB002).
- **k-ε-ζ-f₀ model (stage 6 part 6d2, D132).** `[turbulence] model =
  "k-epsilon-zeta-f0"` after Popovac (2023): f₀ = f + (N_ζ − 1)ζ/T takes
  the place of f and is zero at the wall, ν_t = (C_μ/ζ_eq) ζ k T with the
  constants of the paper. No code of the author is available; the default
  deviates from the printed model where the validation asked for it
  (specification 3.3.1): the sign of the (N_ζ − 1) term of the f₀
  equation, ε/k in the (N_ζ − 1) terms, the k sink k/T_K with the
  Kolmogorov-bounded k/ε, and the realizability coefficient a = 0.6 √2
  of T (with 0.6 the bound halves ν_t in separated shear layers). New
  constants `n_zeta` (ζ-f₀), `c_epsilon1_zeta` and `a` (both ζ-f
  models).
- **Compound wall treatment of 2023.** `wall_treatment = "compound-2023"`
  (default of ζ-f₀, option of ζ-f): the blend
  u⁺ = [y*^{−10/3} + u⁺_log^{−10/3}]^{−3/10} with the non-equilibrium
  factor ψ of the log law from convection and pressure gradient, ε and
  P_k of the first cell from the slope of the blend, f₀,w = 0.
- **Removed boxes in the built-in mesher.** `mesh.generate.remove =
  [{ lower, upper, name }]` cuts boxes out of the hex, tet and mixed
  meshes; the exposed faces form a new boundary (a backward-facing step).
- **Verification (gallery 74, recipe 74).** Channel at Re_τ = 395:
  U_b⁺ +2.7 / +1.4 / −4.2 / −1.5 % at first-cell y⁺ 0.3 / 1 / 10 / 30.
  TMR flat plate: C_f(0.97) −2.5 / +0.9 / −2.0 % from White's correlation
  at y⁺ 1 / 10 / 30, every field converged without a floor, and the
  boundary layer turns turbulent from the TMR free stream (ζ-f stays
  laminar). Backward-facing step of Driver & Seegmiller (new case
  `backstep_zeta_f0`): x_r/H = 6.38, +2 % from the experiment (ζ-f +25 %,
  SST +9 %).
- **Sayings**: two new entries in `docs/sayings.md`, *Deemed prudent* and *Residuals got stuck*.
- **Backward-facing step: three models, two field pictures (gallery 74).**
  `examples/studies/backstep/{keps,zf0,sst}.toml` run the step of Driver &
  Seegmiller on 32 800 cells with k-ε and scalable wall functions,
  k-ε-ζ-f₀ and k-ω SST: x_r/H = 5.36, 6.35 and 6.80 against the
  experiment's 6.26. The gallery shows the velocity with streamlines and
  the turbulent kinetic energy, one panel per model.
- **Figure scripts for 2D fields.** `tools/plot/backstep_paraview.py` makes
  the two pictures in ParaView's Python (`pvpython`; run in ParaView 5.11
  and 6.2), `tools/plot/backstep_chart.py` with matplotlib and h5py
  (recipe 74).
- **Backward-facing step: settled results and half the iterations.** The
  step cases start with the momentum under-relaxation 0.5 and switch to
  0.8 at iteration 100 (`[[schedule]]`), which halves the iterations to a
  settled reattachment length: about 600 on 8 200 cells and 1 750 on
  32 800 cells for k-ε, k-ε-ζ-f₀ and SST. The fine study cases run 3 000
  iterations, `backstep_zeta_f0` 1 500 instead of 4 000. Run to the end,
  SST reattaches at x_r/H = 6.80 (+9 %), not 5.79; k-ε at 5.36 and
  ζ-f₀ at 6.35 on the fine mesh. The figures of gallery 74 are redrawn
  from the settled runs; recipe 74 shows how to check when a result
  has settled.
- **Coupled pressure–velocity algorithm (D176, L11).** `[solver] algorithm
  = "coupled"` solves the momentum equations and the continuity equation
  of a steady outer iteration as one linear system for the increments of
  velocity and pressure: matrix-free with the kernels of the segregated
  solver (no block matrix), flexible GMRES with a SIMPLE-type
  preconditioner, `[solver.coupled]`. The discrete equations and the
  converged result are those of SIMPLEC, and the face fluxes satisfy
  continuity after every outer iteration; the kinetic-energy ledger of
  the cavity at Re = 1000 is the same with both algorithms. CPU and GPU
  from one source (two new kernels).
- **Pseudo-time steps for steady runs (D177, D178).** With the coupled
  algorithm the momentum equations take one pseudo-time step for all
  cells (automatic: the volume of the domain over its boundary area, over
  a velocity scale; `[solver.pseudo_time]`) and the relaxation (default
  0.98) only caps the local Courant number; the turbulence equations take
  the step min(Δτ, 0.5 T) with the time scale T of the model. The
  backward-facing step settles in 65 (k-ε-ζ-f₀), 75 (k-ε) and 80 (SST)
  iterations instead of 525 to 675, on 32 800 cells in 175 instead
  of 1 600; the channel with k-ε-ζ-f₀ converges in 323 instead of about
  2 100, the cavity at Re = 1000 in 452 instead of 1 397 (counts good to
  about 10 %: they react to the last digit of the step). In run time the
  step gains a factor of 1.4 to 1.5 so far, because a coupled iteration
  costs five to seven SIMPLEC iterations; small laminar cases do not gain
  yet (L44, a multigrid for the momentum sweep, is the next step). The
  turbulence step is available with SIMPLEC too
  (`solver.pseudo_time.turbulence = 0.5`, off by default: 230 iterations
  on the step). New cases
  `backstep_zeta_f0_coupled.toml` and `cavity_re1000_coupled.toml`,
  gallery 78, recipe 78, requirements step 22.
- **Log.** The column `K` (`linear.coupled`) shows the Krylov iterations
  of the coupled solve; `residuals.csv` and the run record carry them as
  `krylov` for coupled runs, and the summary of a run shows the
  pseudo-time step of its last iteration (`pseudo_step` in the run
  record). The keys of `[solver.coupled]` and
  `[solver.pseudo_time]` may change during a run.
- **Chart tool.** `tools.plot.settling_chart` draws a monitored value
  against the outer iteration for several runs, with the ±1 % band and
  the iteration from which each run stays inside it.
- **Intermediate results and balance tables no longer change a steady
  run (D31).** A `.result_N` request and the balance table of an
  iteration (`[balance] every`) updated the boundary values and gradients
  for their output and left them in place, so the next iteration started
  from values that had seen one more update: the path of a steady run
  depended on its output schedule (a duct with 24 × 8 cells: 7 % in the
  momentum residual of the iteration after a result; the step case 10⁻⁶
  with SIMPLEC and 10⁻⁴ with the coupled algorithm), though not its
  converged result. Both now put the boundary values, gradients and
  limiters of the iteration back. A run with result files or balance
  tables repeats the plain run bit by bit, with SIMPLEC and with the
  coupled algorithm (test
  `intermediate_results_and_balances_do_not_change_the_run`). Transient
  runs were not affected. The figure of gallery 78 is redrawn from
  histories recorded this way: the SST step settles in 80 coupled
  iterations (85 before). The step on 32 800 cells, run to the end
  with both algorithms (x_r/H = 6.352), settles in 175 coupled iterations
  against 1 600 with SIMPLEC, in 235 s against 327 s.
- **Multigrid for the momentum sweep of the coupled algorithm (D179, the
  momentum part of L44).** The preconditioner of the coupled algorithm
  now applies one V-cycle of a reduction multigrid to each velocity
  component instead of BiCGStab: classical strength of connection, the
  first pass of the Ruge–Stüben splitting into coarse and fine points, a
  restriction that approximates the ideal one by a two-term Neumann
  series (Manteuffel, Ruge & Southworth 2018), direct interpolation, and
  Jacobi relaxation on the fine points only (`amg::reduction`). The
  cycle runs on the level kernels of the pressure multigrid (CPU and GPU
  from one source, float32 hierarchy); the setup is deterministic. On
  the backward-facing step the Krylov iterations per outer iteration
  halve (10 → 5.3 on 8 200 cells, 14 → 6.0 on 32 800), the fine mesh
  settles in 150 instead of 175 iterations, and the run time to the
  settled result on one thread is 7.0 s (SIMPLEC 18.7 s) and 74 s
  (SIMPLEC 243 s, the BiCGStab sweep 131 s); the cavity at Re = 1000 on
  128² converges in 12.1 s (BiCGStab sweep 18.8 s, SIMPLEC 28.3 s).
  Converged results do not change. New keys `[solver.coupled]
  momentum_sweep` (`"multigrid"`, the default, or `"bicgstab"`),
  `momentum_cycles`, `momentum_reuse` and the table
  `[solver.coupled.multigrid]`. Gallery 80, recipe 80, requirements
  step 22 sections 3.1 and 7.7, verification stage 6 part 6f.
- **Kept hierarchies.** A hierarchy of the momentum sweeps is kept while
  the momentum matrix has changed by at most 2 % since its setup
  (`momentum_reuse`) and while its cycle still halves the residual;
  components with the same matrix share one, a component without a
  residual (the third one of a planar case) needs none. The restart file
  holds the matrices of the kept hierarchies (`momentum_amg_use`,
  `momentum_amg_matrix_<n>`), so a continued run stays bit-identical;
  older restart files still load. If a setup fails, BiCGStab does the
  sweeps of that outer iteration and the log carries a warning. The
  report of a run counts the setups of the momentum sweeps and gives
  the time of the turbulence model and of the scalars among its
  seconds per phase.
- **One thread without a pool.** On a pool of one thread (meshes below
  20 000 cells, or `--threads 1`) every kernel call and every level
  operation of a multigrid cycle was handed to the worker thread and
  waited for, about 2 000 times per outer iteration of a cavity run.
  Such loops, and loops of a single task on any number of threads (the
  coarse levels of a multigrid cycle), now run on the calling thread.
  The values are the same bit by bit. Small and medium cases gain with
  both algorithms: the step on 8 200 cells takes 18.7 s instead of 31 s
  with SIMPLEC (525 iterations) and 8.3 s instead of 21 s with the
  coupled algorithm and the BiCGStab sweep (65 iterations); the step on
  32 800 cells with SIMPLEC 243 s instead of 327 s.
- **Chart tool.** `tools.plot.krylov_chart` draws the Krylov iterations
  of the coupled solve against the outer iteration for several runs.
- **Multigrid for the turbulence equations and larger steps with it
  (D180, the turbulence part of L44).** With the coupled algorithm the
  transport equations of the k-ε family (k-ε, realizable k-ε, k-ε-ζ-f,
  k-ε-ζ-f₀) take one cycle of the reduction multigrid of D179 as the
  preconditioner of BiCGStab, with kept hierarchies (`[solver.turbulence]
  preconditioner = "reduction"`, the default there, and `reuse`), and
  the relaxations of the momentum and turbulence equations default to
  0.995 instead of 0.98. The relaxation acts as a local step that
  shrinks with the cells: at 0.98 it made the outer iterations grow with
  the mesh. The backward-facing step settles in 40 iterations on 8 200
  cells and 65 on 32 800 (65 and 150 before), k-ε in 45 (75), the
  channel converges in 231 (323). One thread: the fine step takes 56 s
  instead of 74 s; the coarse steps about the same time as before (7.5 s
  against 7.0 s, k-ε 5.7 s against 4.5 s), since an iteration costs
  more. SST keeps Jacobi and 0.98 (its steady state on the step is not
  held with accurate solves of its equations), laminar cases keep 0.98.
  Converged results do not change. The restart file holds the matrices
  of the kept hierarchies; a failed setup falls back to Jacobi with a
  warning; the report counts the setups. Gallery 81, recipe 81,
  requirements step 22 sections 4.2 and 7.8, verification stage 6 part
  6g.
- **Coupled algorithm: looser pressure solves in the preconditioner, and
  the Schur complement measured.** `[solver.coupled] pressure_tolerance`
  defaults to 0.3 instead of 0.1: the Krylov iteration corrects what the
  pressure solve of the preconditioner leaves, so the Krylov count stays
  and the CG iterations fall by a fifth to more than a quarter. One thread, to the settled
  result: the coarse step takes 7.4 s instead of 8.0 s, k-ε 5.9 s instead
  of 6.4 s; the fine step is unchanged (57.1 s against 56.9 s), since
  its pressure solves stop at their iteration limit either way; the
  results do not change. Two published alternatives to the SIMPLE-type Schur complement
  were measured and do not pay (SIMPLER: a second pressure solve for a
  few Krylov iterations; the least-squares commutator: fine on the
  cavity, at the Krylov limit on the step). Requirements step 22
  section 7.9; the next step is a pressure hierarchy refreshed in every
  outer iteration.
- **Coupled algorithm: the pressure hierarchy of the preconditioner
  follows its matrix (D181, stage 6 part 6h).** The multigrid hierarchy
  of the pressure solve inside the coupled preconditioner is built anew
  when its matrix has changed by more than
  `[solver.coupled] pressure_reuse` (5 %) since its setup, besides every
  `amg_rebuild` iterations as before (`pressure_reuse = 1` restores the
  periodic rebuild alone). A hierarchy up to nine iterations old needed
  9 to 20 CG iterations per solve where a fresh one needs two. One
  thread, to the settled result: the fine backward-facing step (32 800
  cells) takes 39.8 s instead of 57.1 s, 6.1 times faster than SIMPLEC;
  the coarse step 6.2 s instead of 7.4 s, k-ε 4.7 s instead of 5.9 s.
  The settle counts and the converged results do not change; laminar
  cases do not trigger the rule. Requirements step 22 section 7.10,
  gallery and recipe 83.
- **Multigrid setup: flat row builders (L44, no change of results).** The
  reduction multigrid of the coupled algorithm (momentum and turbulence
  equations) builds its strong graph, interpolation, restriction and
  Galerkin products into flat arrays, blocks of rows in parallel joined in
  row order, with one dense accumulator per worker instead of a vector
  per row. The matrices are bit-identical on any number of threads, and
  runs give the same residual history bit for bit. One thread: the
  multigrid setups of the fine backward-facing step take 4.9 s instead of
  6.9 s, the run 33.3 s instead of 39.8 s (7.3 times faster than
  SIMPLEC); the coarse step 5.7 s instead of 6.2 s. Requirements step 22
  section 7.11.
- **Passive scalars of a turbulent flow: turbulent diffusivity and thermal
  wall function (D135, stage 6 part 6f).** A scalar now diffuses with
  ν/Pr + ν_t/Pr_t on every face, and a wall face carries the conductance
  h_w = ρ u*/T⁺ of the wall law `scalars.<name>.wall_function` (`kader`,
  the default, valid for all y⁺ and any Pr, `jayatilleke`, or `none`).
  `turbulent_prandtl` (0.85) and `turbulent_prandtl_model` (`constant`
  or `kays-crawford`) set Pr_t. Fixed-flux walls report the wall value
  from T⁺. The wall file gains `<name>_wall`, `<name>_flux`, `<name>_h`
  and `<name>_t_plus`. The heated pipe at Re = 10⁵, Pr = 0.71 gives Nu
  182.6 (Kays–Crawford 181.3) against Gnielinski's 180.2; the flat
  plate St is +1.7 % to +7.5 % from the correlation for 0.25 ≤ x ≤ 1.
  Cases without a turbulence model do not change. Requirements step 6
  section 5, gallery and recipe 85.
- **Heated channel against DNS (stage 6 part 6f, specification 9.4).**
  `channel_heated_keps.toml` is the Re_τ = 395 channel with a passive
  temperature (Pr 0.71, uniform wall heat flux, periodic heat-flux jump).
  Its Nusselt number on the width 2h is 36.01 with the constant Pr_t and
  35.65 with Kays–Crawford, against 36.39 of the DNS of Abe, Kawamura &
  Matsuo (−1.0 % and −2.0 %). `examples/studies/channel_heat_dns.py`
  fetches the DNS tables of the JAXA database into a git-ignored folder
  and recomputes their Nusselt numbers for Re_τ 180 to 1020 at Pr 0.71
  and 0.025; no DNS data enters the repository (JAXA reserves the
  rights). Gallery and recipe 86.
- **Heated channel at every DNS point (stage 6 part 6f).** The channel of
  `channel_heated_keps.toml` at Re_τ 180, 395, 640 and 1020 and Pr 0.71
  and 0.025 against the DNS of Abe, Kawamura & Matsuo, on uniform meshes
  with first-cell y⁺ ≈ 30: Nu within −2.6 … +1.0 % at Pr 0.71 with both
  Pr_t models; at Pr 0.025 +5.9 … +7.2 % with Kays–Crawford and
  +12 … +22 % with the constant Pr_t. New case `channel_heated_pr0025.toml`;
  `examples/studies/channel_heat_dns.py --run` runs all sixteen points.
  The thermal law now takes T⁺ at the scalable treatment's floor
  ỹ⁺ = max(y⁺, 11.06), as the momentum law does (channel Nu on a y⁺ ≈ 8
  mesh +10 % instead of +28 %). CG stops when the preconditioned residual
  vanishes below the float32 range of the multigrid (an exactly
  one-dimensional flow), as BiCGStab did. The PTX module now has the
  entry points of the two passive-scalar kernels of pull request 85
  (`scalar_turb_faces`, `scalar_wall_faces`): without them a GPU run did
  not load the module; a test checks that the module has every entry the
  CUDA backend loads. Gallery and recipe 87.
- **Elliptic-blending k-ε (stage 6 part 6e, D133).** New model
  `elliptic-blending-k-epsilon`, the BL-v²/k model of Billard & Laurence
  (2012) with the equations and constants of the paper: k, ε, φ = v²/k
  and the blending parameter α (α − L²∇²α = 1, α = 0 at walls), the E
  term in the k equation and the variable C*_ε2, ν/2 in the molecular
  diffusion, ε_w = νk/y²; wall treatments `compound` (default) and
  `low-re`. It runs on the fields, solvers and multigrid of k-ε-ζ-f (φ in
  the place of ζ, α in the place of f); new element-wise kernels
  `squares_add` and `scale_vector` (also in the PTX module). Channel C_f
  at 97.8 / 99.0 / 99.1 % of the DNS at Re_τ 180 / 395 / 590 (the paper:
  100 / 98 / 100 %); TMR plate C_f(0.97) +0.5 % from White; backward-facing
  step x_r/H = 7.59, +21 % (finding: the variable C*_ε2;
  `constants = { c_epsilon4 = 1.83 }` gives 6.24 and a plate 4.9 % high).
  Output fields `phi` and `alpha`, log columns `ResPhi` and `ResAlpha`;
  the key `zeta` sets φ. New cases `channel_eb.toml`, `tmr_plate_eb.toml`,
  `backstep_eb.toml`; specification 3.4 corrected after the paper.
  Gallery and recipe 88.
- **Backward-facing step against the experiment (stage 6 part 6g,
  specification 9.5).** New case `backstep_validation.toml`: the step of
  Driver & Seegmiller (1985) in the layout of NASA's Turbulence Modeling
  Resource (walls from x = −110 H, outlet at 50 H, no-slip top wall at
  9 H). New study `examples/studies/backstep_validation.py`: six models
  on three meshes, the reattachment length with its GCI, C_p and C_f at
  the experiment's stations and the profiles at x/H = 1, 4, 6, 10; the
  experimental data of the TMR (CC0) in `examples/studies/backstep/data`.
  New output `[[output.wall_crossings]]`: the first point along a
  direction where the wall shear changes sign (reattachment or
  separation), in the report and with its observed order and GCI in a
  level sweep; `skip` passes the corner vortices behind a step.
  Results on the finest mesh: SST x_r/H 6.747 (+7.8 %), k-ε 5.372
  (+0.5 % from the published k-ε, GCI 1.94 %), realizable k-ε 6.407,
  k-ε-ζ-f₀ 6.296 (+0.6 %); ζ-f and EB 7.74 and 7.56 on the middle mesh,
  outside the band (their finest mesh follows). The run banner names
  the turbulence model and wall treatment instead of "laminar". Gallery
  and recipe 89.
- **Impinging jet at H/D = 2 (stage 6 part 6g, specification 9.6).** New
  cases `jet_pipe_profile.toml` (the developed pipe flow at Re = 23 750
  on a periodic wedge) and `impinging_jet.toml` (the round jet from that
  pipe on a heated plate, axisymmetric wedge, a passive temperature with
  a uniform wall flux). New study `examples/studies/impinging_jet.py`:
  six models on three meshes, Nu₀ with its GCI, Nu(r/D) with its
  secondary maximum, the stagnation-point anomaly with and without the
  production limiters, and the wall-jet velocity profiles; it compares
  with the experiment once the data are in
  `examples/studies/impinging_jet/data`. On 41 600 cells (with the
  fix of #91) Nu₀ at the axis: ζ-f 173.8 (GCI 0.12 %), SST 141.2
  (2.9 %), EB 249.2 (0.19 %); Nu at r/D = 0.1, the second probe `r01`:
  ζ-f 173.9, EB 205.0, SST 141.1 (GCIs 0.03 to 0.40 %), ζ-f₀ 151.7
  (GCI 5.3 %). The k-ε without a limiter has the anomaly (Nu₀ 211), its
  Durbin default removes 71 % of the excess k; SST has none (its ν_t
  bound holds on the axis). Gallery and recipe 90.
- **Inlet tables with the turbulence of the profile (specification 7.1).**
  The header of a `profile = "table"` file names its columns: `u_n` or
  `u_x`, `u_y`, `u_z`, and `k` with `epsilon` or `omega`, and `zeta`.
  The turbulence columns set the inlet's k, ε or ω and ζ face by face; a
  boundary that gives turbulence values as well is refused (BC029), an
  unknown column too (BC019). `python3 -m foehncraft.profiles LINE.csv
  TABLE.csv` makes the table from a line of a periodic pipe or channel
  run (new module `foehncraft.profiles`).
- **Wall probes report the passive scalars.** Each `[[output.wall_probes]]`
  gives the wall value and the diffusive flux of every scalar there, in
  the log and the report, and as GCI quantities of a level sweep
  (`wall_probes.<name>.<scalar>_wall`, `_flux`).
- **Fix: the thermal wall law of the k-ε family and ζ-f/EB takes
  u* = C_μ^¼ k_P^½ (specification 5.2), not u_τ.** u_τ vanishes with the
  wall shear, so the heat transfer of a stagnation point or a
  reattachment line went to zero (the k-ε Nu₀ of the jet was 37 instead
  of 156). SST keeps its blended u_τ. The heated channel moves by −0.1 to
  −1.2 points (still within 2.9 % of the DNS at Pr 0.71), the heated pipe
  by +0.3; the flat plate is +8 to +9 % for 0.4 ≤ x ≤ 1, and +14 % at
  x = 0.29 just after the k-ε start of the turbulent layer, where k
  overshoots (the test window moves from 0.25 to 0.4).
- **Fix: the full viscous stress on symmetry planes, slip walls and wedge
  planes (stage 6 part 6g).** With a turbulence model the momentum rows
  take the transposed part of the stress, ∇·[μ_eff (∇u)ᵀ], on the
  interior faces; on symmetry faces it was missing. There u_n = 0, so
  that part is μ_eff |S| ∂u_n/∂n along the normal, and the normal stress
  on the face is 2μ_eff ∂u_n/∂n (implicit, as before with the factor 1).
  On the two planes of a wedge it is half of the hoop stress
  −2μ_eff v/r²: the cells next to the axis of the impinging jet had v/r
  too high and Nu too low, more so on finer meshes (the open item of
  #90). Laminar runs do not change. New test
  `tests/symmetry_stress.rs`: a linear flow with constant stress balances
  in every cell, on a wedge and on a symmetry plane (the residual was
  4 % and 65 % of a face flux).
- **Validation: the impinging jet against the ERCOFTAC data (stage 6 part
  6g, part 4).** The six models of the jet study are compared with
  the experiment (ERCOFTAC classic database, case 025, H/D = 2: Nu/Re^0.7
  of Baughn et al. scaled to Re 23 750, and the mean velocity of Cooper
  et al.). Nu₀ of SST (−0.2 %) and ζ-f₀ (+7.2 %) is inside its band, that
  of ζ-f (+22.9 %) and EB (+76 %) outside; the RMS of Nu(r/D) is met only
  by SST (13.7 %); standard k-ε without a limiter is +49 % (the band asks
  for more than +50 %). The data are not in the repository (the wiki of the
  database is licensed CC BY-NC-SA 4.0): `examples/studies/impinging_jet/ercoftac_to_csv.py` makes the
  CSV files from the database files, `.gitignore` keeps them out, and the
  study's new `--no-experiment` makes the gallery charts without the
  points. Verification page (finding 8), specification 9.6 and gallery 90
  updated.
- **Requirements: sliding mesh in both variants (D182).** The sliding mesh
  is planned together with the moving reference frame in stage 6 (6g parts 5
  to 7): conformal stepping (variant A, D56 moves from stage 8) and general
  sliding with the interface imprinted at every step (variant B, was L4).
  Specification, L4, D56, D58 and the stage rows updated; documentation only.
- **Moving reference frame (D55) with laminar Taylor–Couette flow (stage 6
  part 6g, part 5).** `[[mrf]]` tables turn a cell zone steadily about an
  axis: the solved velocity is the absolute one, momentum, turbulence
  variables and scalars are convected with the flux relative to the frame,
  the rotating cells get the Coriolis source, continuity and the pressure
  correction stay on the absolute flux. New: `[mesh.annulus]` (the gap
  between two cylinders with an optional `rotor` zone), `[analytic]
  type = "taylor-couette"`, the diagnostics `MRF001` to `MRF003`, the case
  `taylor_couette_mrf.toml` (L2 error 0.69 % on 1 024 cells, observed order
  2.3 and 2.4, the same in every frame), spec 9.7, verification page, gallery
  and recipe 94. Steady only; the sliding mesh follows (D182).
- **Sliding mesh with conformal stepping (D56, D182, stage 6 part 6g, part
  6).** `[[sliding]]` turns a cell zone of a ring mesh by whole sectors in
  every time step (Δt = shift · 2π/(N |Ω|), refused with the right value
  otherwise, `SLD003`), so interface cells meet new neighbours without
  interpolation and the mesh stays conformal. The solved velocity is the
  absolute one, momentum and scalars are convected with the flux relative to
  the moving mesh (the frame flux of D55, without the Coriolis source), and
  only the time term sees the motion: the old time levels of a rotor cell are
  those of the cell that stood at its place. Laminar flow and passive
  scalars; turbulence models, pressure boundaries on the rotor and steady
  cases are refused (`SLD001` to `SLD004`). New: `[initial] from_analytic`
  for transient runs, the case `taylor_couette_sliding.toml` (L2 error
  1.35 % on 512 cells, observed order 2.2 to 2.3 from 128 to 8 192 cells),
  spec 9.7, verification page, gallery and recipe 95. The variant with a free
  pitch ratio (imprinting) follows.
- **Sliding mesh by imprinting (D182 variant B, stage 6 part 6g, part 7).**
  A rotor and a stator with different numbers of sectors meet on an
  interface that is imprinted at the rotor angle of every time step, so the
  time step is free. `[mesh.annulus] n_theta_stator` makes the two rings
  with the boundaries `rotor_interface` and `stator_interface`;
  `[[sliding]] rotor_interface / stator_interface` names them. Every
  intersection of a rotor face with a stator face is an ordinary internal
  face (`foehncraft_mesh::imprint`: polygons clipped in the unrolled
  surface of revolution), the cells close, mass and momentum flux are
  conserved to round-off, and the solver is built on that mesh in every step
  with the state, the old time levels and the face fluxes carried over. The
  scheme is fully implicit in the step (`extrapolation = false`, at least
  two outer iterations, `SLD003`); `SLD005` and `SLD006` report an interface
  that cannot be imprinted. Taylor–Couette flow from rest with 32 rotor and
  48 stator sectors: L2 error 6.74, 1.29 and 0.25 % on 160, 640 and 2 560
  cells (order 2.4); equal pitch agrees with conformal stepping (1.350 % and
  1.354 %); solid-body rotation is kept; a passive scalar ends at 0.5 %.
  Laminar flow and passive scalars, one rotor, a fixed time step, no
  restart, CPU. Also: the VTKHDF writer writes cells whose faces do not
  close on the nodes of their shape as polyhedra, and probes in a moving
  rotor find their cell again in every step. New case
  `taylor_couette_imprint.toml`, gallery and recipe 96, spec 9.7.
- **Faster rebuild of the sliding mesh by imprinting (D182 variant B).**
  The solver and its pressure multigrid are still built for every step (the
  hierarchy of the step before was measured and does not fit: the rotor
  turns by about one cell per step, and the conjugate gradients needed 5 to
  11 times the iterations), but each build costs less: the dense inverse of
  the coarsest multigrid level takes a third of the work over contiguous
  rows (5 ms to 0.9 ms for 175 unknowns, in every setup of every run),
  multigrid setups below 20 000 rows run on one helper thread instead of
  handing each short parallel loop to the pool, the mesh geometry and the
  discretisation of a step are computed once instead of twice, and the
  state is moved instead of copied. The case file `taylor_couette_imprint.toml`
  runs in 5.6 s instead of 8.4 s, 2 560 cells in 45 s instead of 51 s, with
  the same results (L2 errors unchanged, histories bit-identical or to
  round-off). Verification page stage 6, finding 3 of the sliding mesh by
  imprinting.
- **Walls of a turning zone push the fluid (fix).** A wall in a zone of
  `[[mrf]]` or `[[sliding]]` kept only the tangential part of its `velocity`,
  the rule for walls that stand (D53, `BC010`), so a blade, which moves
  normal to itself, slid along the fluid instead of pushing it. Such a wall
  now moves normal to itself with the zone (Ω × r) and along itself with the
  given velocity: absolute flux ρ **S**·**V**, relative flux zero, as D55
  states. Cylinders that turn about their axis (gallery 94 to 96) are not
  affected. New test: a box turning with its frame carries its fluid (L2
  1.8·10⁻² before, 1.3·10⁻³ after).
- **A stirred vessel with blades and baffles, and field figures of rotating
  machinery (gallery 98, recipe 98).** `examples/meshes/mixer.py` meshes a
  2D baffled vessel with a four-blade rotor in Gmsh (conformal for the
  frame, split 64/96 for the sliding mesh); `examples/studies/mixer/` has the
  steady frame case (Re 100) and the sliding-mesh start-up over three
  revolutions (Re 400); `examples/studies/mixer.py` runs both and draws the
  figures. `tools/plot/rotor_fields.py` draws any 2D rotor run: speed with
  streamlines, the speed relative to the blades in the rotor, vorticity,
  pressure, time series and zooms with the cell edges. Gallery 94, 95 and 96
  have field figures now too.
- **The flow a blade sees: relative and frame velocities of rotating
  machinery (D183, gallery 99, recipe 99).** With an `[[mrf]]` or a
  `[[sliding]]` table, `output.derived` takes five new fields in the
  notation of turbomachinery (absolute **c**, blade speed **u** = Ω × r,
  relative **w** = **c** − **u**): `relative-velocity` (**w** in each turning
  zone, **c** elsewhere, so a pump shows the relative flow in the impeller
  and the absolute flow in the volute), `frame-velocity`, `turbomachinery`
  (c_r, c_u, c_a, c_m, w_u about the axis and the flow angles α and β),
  `relative-total-pressure` and `rothalpy`. `output.frame = "<zone>"` shows
  every cell in the frame of one zone. Sample lines get the columns of **w**
  and of the cylindrical components, point probes of a transient run those
  of **w**; on the sliding mesh at the current position of the cell.
  `tools/plot/rotor_fields.py` draws any dataset of the file, also in the
  rotor only (`NAME@rotor`).
- **Fix:** line samples and the reference comparison of the sliding mesh by
  imprinting located their points on the mesh of the start of the run, whose
  rotor cells had moved on; they use the mesh of the end of the run now.
- **Turbulence on rotating machinery: the moving reference frame and both
  sliding meshes (gallery 100, recipe 100).** Every turbulence model now runs
  with `[[sliding]]` (it was refused, `SLD004`) and is checked with
  `[[mrf]]`. With conformal stepping the old levels of k, ε, ω, ζ and f move
  with the rotor and each step starts from the moved fields; with imprinting
  the turbulence state is carried through the rebuild of every step. The wall
  distance of the turning rotor is kept per rigid part (one bounding volume
  hierarchy each) and only the cells nearer to the edge of their part than to
  its walls search the other part in each step, with a capped query that
  starts from the last nearest face. Several independent `[[mrf]]` zones are
  tested (two zones turning in opposite senses).
- **Fix:** the pressure multigrid built a hierarchy of NaN when a coarse row
  had a non-positive lumped diagonal in the prolongator smoothing (`GEN001` in
  the first iteration of the turbulent stirred vessel with the frame); such a
  row keeps its own diagonal now.
- **Torque, shaft power and machine performance of rotating machinery
  (D184, gallery 101, recipe 101).** For `[[mrf]]` and `[[sliding]]`: the
  torque of the fluid on every wall about the axis of its rotor (the walls of
  a rotor are its `walls`, by default the walls inside its zone), in the
  boundary table, `boundaries.csv`, `forces.csv` and the report; the shaft
  power −Ω·T of each rotor, booked in the step and the final energy balance
  (`balance.csv` column `shaft_power`, the rows of the ledger tables) and in
  the compact balance line; the angular momentum balance about the axis of
  each rotor. `[[output.machines]]` names a rotor, an inlet and an outlet and
  gives the volume flow, the total pressure rise, the head, the hydraulic and
  total-to-static efficiency and the Euler work ṁ Ω Δ(r c_u). New log
  columns `rotor.<zone>.torque|power|angle` (the angle of every sliding rotor,
  in degrees, each with its own speed) and `machine.<name>.flow|pressure_rise|head|efficiency`
  (all appear by default). The convected kinetic energy of the ledger uses the
  flux relative to a turning zone. The torque of a moving wall includes the
  transposed viscous stress: the Laplacian form of the momentum rows misses
  the −u_θ/r of the shear stress (a turning cylinder: 35 % too small before,
  0.5 % against the exact torque now). `tools/plot/rotor_movie.py` writes the
  frames of an animation of a field series in ParaView.
- **Several sliding rotors and the transient moving reference frame (D185,
  gallery 102, recipe 102).** One `[[sliding]]` table per rotor, each with
  its own zone, speed, axis, origin and interface to the stator: all
  interfaces are imprinted together in every step at the angles of their
  rotors (variant B of D182; conformal stepping stays one rotor, an
  interface between two rotors is refused). The exact wall distance of
  turbulence models follows every rotor. The log headings of rotor columns
  name the rotor. `[mesh.annulus] outer_rotor_layers` and `n_theta_outer`
  make a second rotor ring at the outer cylinder; `examples/meshes/twin_mixer.py`
  makes a vessel with two stirrers. `[[mrf]]` now runs in transient cases
  too (the transient frozen rotor), also next to sliding rotors in other
  zones; `MRF001` refuses a zone that a frame and a sliding mesh both turn
  (it refused every transient frame before). New cases
  `examples/cases/flow/couette_two_rotors.toml`,
  `examples/studies/twin_mixer/sliding.toml`,
  `examples/studies/mixer/mrf_transient.toml`.
- **Sliding interfaces of any rotationally symmetric shape (D186, gallery
  103, recipe 103).** One interface per rotor, of any surface of revolution
  about its axis whose meridian is made of straight pieces: a cylinder, a
  disc (an axial gap between blade rows), a cone, or a chain of them such as
  a pot of disc, cylinder and disc around an impeller; several separate
  pieces (discs before and after a rotor row) and discs through the axis
  too. The meridian is found from the mesh; the faces are clipped in the
  unrolled (angle, meridian length) coordinates with a polar chart at the
  axis, and the corners of the pieces are inserted into the neighbouring
  face edges, so every cell of an imprinted mesh closes to round-off on any
  shape. The sliding log line names the shape. New built-in mesh
  `[mesh.revolved]`: a block grid in (r, z) turned about the axis, with
  rotor blocks (own pitch), solid blocks (walls such as an impeller disc),
  cones and wedges at the axis. New cases
  `examples/cases/flow/axial_stage_uniform.toml`,
  `examples/cases/flow/pot_solid_body.toml`,
  `examples/studies/revolved/disc_pot.toml`,
  `examples/studies/revolved/swirl_stage.toml`; new plot tool
  `tools/plot/meridional_fields.py`.
- **README up to date, highlight images, gallery links that land.** The
  README lists what works on main ("What works today": turbulence models
  and wall treatments, heat transfer with passive scalars and thermal wall
  functions, the coupled algorithm, rotating machinery, …) under a banner
  and six highlight tiles linked to their gallery entries
  (`assets/readme/`). Gallery and README images carry their width and
  height (`tools/docs/images.py`, checked in CI), so a link to an entry far
  down a topic page no longer lands thousands of pixels too high while the
  images above it load. Recipe 33's anchor is `python-api` (it clashed with
  the page heading). New rule in CONTRIBUTING.md and the pull request
  template: every pull request keeps the README and the docs current.
- **Sayings**: new entry *Better to know what's not happening* in `docs/sayings.md`.
- **PISO corrector of SIMPLEC type with outer iterations: stable on thin
  wedge cells at an axis at high viscosity (gallery 103).** With
  `outer_iterations` > 1 the pressure corrections of PISO use the SIMPLEC
  coefficient D^C = V/(a_P/α − Σ|a_N|) with the consistent neighbour
  correction Σ|a_N| (u'_N − u'_P)/a_P between two correctors, instead of
  Issa's D = V/a_P with H(u')/a_P. Where a_P is almost all coupling to the
  neighbours (the wedge cells about an axis at a fine pitch and high
  viscosity, a_P about 50 times the time term), Issa's form let the
  pressure grow in every outer iteration: the disc impeller in a pot
  diverged at ν = 1e-4 m²/s. It now runs at that viscosity
  (`examples/studies/revolved/disc_pot.toml` is back at ν = 1e-4). Plain
  PISO keeps Issa's form, which damps more there (bitwise unchanged); new
  key `solver.corrector` (`auto`, `issa`, `simplec`). A converged step is
  unchanged. New case `examples/cases/flow/pot_solid_body_viscous.toml`:
  solid-body rotation at ν = 1e-2 m²/s, L2 error 0.11 % (2.0 before).
- **Rotor images with five blades (gallery 96, 98 to 100, 102, README).** The
  flow fields of the four-blade stirred vessel drew a four-armed pinwheel
  that can be misread as a symbol. Every flow image of a four-blade rotor
  is replaced by the same case on a vessel with five blades and three
  baffles, and the twin vessel (three and four blades) by three and five
  blades. The four-blade meshes, case files, tests and the numbers and
  torque plots from them stay. `examples/meshes/mixer.py --blades
  --baffles`, `examples/meshes/twin_mixer.py --right-blades`,
  `examples/studies/mixer.py --blades --rotors` (all gallery and README
  rotor images in one script). New rule in CONTRIBUTING.md: no images that
  could be read as a symbol; rotors in published images have three or five
  blades.
- **Derived fields of the turbulence models (D72, gallery and recipe 108).**
  `[output] derived` and `foehn post derive` take `turbulence-scales`
  (ε, ω, the time scale of the model, k^{3/2}/ε, ν_t/ν for every model),
  `reynolds-stress` (the Boussinesq estimate ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ, six
  components), `wall-normal` (∇α/|∇α| of the elliptic-blending k-ε, for the
  GGDH of stage 7) and `turbulence-wall` (u_τ of the wall treatment and y*
  on the walls). New plot `tools/plot/turbulence_fields.py`; tests
  `tests/turbulence_derived.rs`.
- **GPU parity check of the turbulence models (specification 6.6).**
  `tools/bench/stage6_cloud.py` and `tools/notebooks/stage6_colab.ipynb`
  run every turbulence model on the TMR plate and the channel, the heated
  channel and pipe and the coupled step on the CPU and the GPU and compare
  each pair bitwise; `--cpu-only` compares one thread with four.
  `tools/check/compare_runs.py` no longer reports the hash of
  `case.resolved.toml` in `result.json` as a difference when only the
  device keys of `[run]` differ (`--threads`, `--backend`).
- **Statistics of transient runs (D71, gallery and recipe 110).**
  `[output.statistics]` takes time means, RMS and covariances of the fields
  (`statistics.vtkhdf`: velocity, pressure, k, scalars, the resolved
  Reynolds shear stresses and scalar fluxes) and of the probes and force
  coefficients (`statistics.csv`), weighted with the time step, from
  `start` to the end. Each series gets the 95 % confidence interval of its
  mean from batch means and a verdict from the lag-1 correlation of the
  batch means; a window that is too short warns with STAT001.
  `phase_period` writes phase averages (`statistics_phase.csv`). Checked on
  the decaying Taylor–Green vortex against its exact time moments
  (`tests/statistics.rs`).
- **Grid convergence of the turbulence models (spec 9.3, 9.4).**
  `examples/studies/plate_channel_gci.py` runs all six models on three
  meshes of the TMR flat plate and the Re_τ 395 channel and reports the
  observed order, the extrapolated value and GCI_fine (Celik 2008) of
  C_f, C_D and U_b⁺ against White, the TMR and the DNS. Findings on the
  stage-6 verification page: no model is in the asymptotic range on the
  plate meshes (SST within 0.14 % of the TMR value, GCI 0.47 %), ζ-f does
  not converge on the finer plates, the realizable k-ε with scalable wall
  functions moves with every refinement of the channel.
- **SST `automatic` wall treatment holds C_f in the buffer layer (6g).**
  The near-wall ω is blended with the exponent 1.2 and imposed, and the
  first cell takes the equilibrium production of Reichardt's law, both with
  the weight tanh((y⁺/5)²). On the TMR plate C_f now lies within −2.6 to
  +2.4 % of the reference for first-cell y⁺ 0.24 to 85 (before: up to
  +15.6 % at y⁺ 4.6). TURB001 warns for `automatic` only above y⁺ 300.
  The sweep is `examples/studies/sst_wall_sweep.py`; gallery 68 shows it.
- **SST `low-re`: the wall value of ω checked.** A sweep of the factor in
  ω_w = F · 6ν/(β₁ d₁²) from 1 to 100 on the plate and the channel keeps
  F = 10, which has the smallest mesh dependence (verification manual).
- **Curvature and rotation correction of SST (6h, D170, gallery 113,
  recipe 113).** `turbulence.curvature_correction = true` multiplies the
  production of k and ω by the rotation function f_r1 of Spalart & Shur
  (1997) in the form of Smirnov & Menter (2009), from r* = S/Ω and the
  material derivative of the strain rate (constants `c_r1`, `c_r2`, `c_r3`).
  It sees the rotation of `[[mrf]]` zones; not with `[[sliding]]` yet.
  In pure shear the result is unchanged to the last digit.
- **A whole domain in a turning frame.** `flow.frame_rotation = [Ωx, Ωy, Ωz]`
  solves for the velocity relative to a frame that turns at Ω, with the
  Coriolis force (the centrifugal force stays in the pressure), for flows
  that are periodic in the rotating frame. The new case
  `channel_rotating_sst.toml` is the channel at Re_b 5 800 in spanwise
  rotation (Ro 0.15). There the wall shear stress of the unstable side is
  1.63 times that of the stable side with the correction (1.00 without).
- The derived field `curvature-correction` writes f_r1. TURB010 notes the
  correction together with Kato–Launder production.
- **Requirements step 24: elliptic-blending Reynolds-stress model as a later
  item (L47).** Review of a supplied EB-RSM specification (Manceau 2015) on
  the new page `docs/requirements/step24-eb-rsm.md`, with the rotation term
  per frame formulation, tensor boundary conditions, realizability by
  anisotropy scaling and the checks; the algebraic heat-flux companion is
  P18, open points O21–O23. No code.
- **Quadratic constitutive relation of SST (6h, D170, gallery 115,
  recipe 115).** `turbulence.qcr = "qcr2000"` adds the non-linear stress
  of Spalart (2000) to the momentum equations, and `"qcr2013"` adds the
  isotropic term of Mani, Babcock, Winkler & Spalart (2013) (constants
  `c_cr1`, `c_cr2`). The stress goes onto the momentum faces as an
  explicit stress, the path the later algebraic stress model reuses. The
  new case `square_duct_qcr.toml` (Re_b 3 500) shows the eight corner
  vortices with a largest secondary velocity of 1.56 % of U_b. Linear
  SST gives none.
- **Requirements: O22 answered.** The rotating channel of the EB-RSM (L47)
  and the EARSM (D171) runs on the whole-domain relative frame
  `flow.frame_rotation` of PR #113; the step 24 page points to it, and the
  frame term C_ij in the stress equations stays part of L47.
- **Validation reference data: sources and attribution.** New page
  `docs/verification/DATA_SOURCES.md` lists the published sources of the
  reference values of the validation cases (channel, square duct, rotating
  channel, EB-RSM benchmark) with their citations and the quantities used. The
  repository keeps the quoted numbers only, not the data sets; no fetch script
  ships.
- **Validation data sources: citations verified.** The DOIs of Manceau 2015 (IJHFF 51), Pinelli et al. 2010 (JFM 644) and Gavrilakis 1992 (JFM 244) are confirmed in `docs/verification/DATA_SOURCES.md`; no entry is left as "to verify".
- **Validation data sources: channel DNS citation.** `docs/verification/DATA_SOURCES.md` now cites Lee & Moser (2015, JFM 774) for the Re_τ = 180–5200 channel set, with its DOI. The Re_τ 395 and 590 entries keep the Moser, Kim & Mansour (1999) citation, since the 2015 set has no case at those values.
- **Validation data sources: MKM 1999 archive.** `docs/verification/DATA_SOURCES.md` now names the Moser, Kim & Mansour (1999) archive page and its Re_τ = 395 profile listing as the sources of the Re_τ 395 and 590 reference values. The Re_τ = 590 listing is still marked "to verify".
- **ζ-f on the finest meshes.** The multigrid setup of an equation with a
  wide range of coefficients (the f equation of k-ε-ζ-f on a fine mesh)
  refused a healthy coarsest level as singular, because its last pivot was
  compared with the largest diagonal. For the turbulence equations each pivot is
  now compared with its own diagonal (the pressure keeps the strict test). A transport or elliptic equation whose multigrid setup
  fails no longer stops the run: its solve takes the Jacobi preconditioner
  and the log notes it. The GCI study runs the level-3 ζ-f plate and channel with
  the relaxations 0.5, since both diverge with the defaults.
- **GCI levels 1 to 3 on the maintainer's machine.** The plate and channel
  results of the finest three meshes of the turbulence GCI study are in the
  stage 6 verification notes, with findings: SST plate GCI 1.87 %
  (+0.37 % from the TMR value), k-ε plate 2.39 %.
- **Level-3 channels rerun in the cloud.** SST, ζ-f₀ and EB reproduce the
  maintainer's values exactly; ζ-f does not converge on level 3 in 40 000
  iterations. Recorded in the stage 6 verification notes (the cause of the
  drop follows in #128).
- **Industrial applications track proposed (step 25).** A proposal page for an "Applications" gallery with its cases and recipes next to the validation and feature gallery: layout, numbering, an entry template, checks without reference data (balances, convergence, mesh sensitivity, wall resolution, plausibility) with three trust levels, 14 candidate cases that run with the features on main, and the CI plan. Nothing is implemented yet ([step 25 page](docs/requirements/step25-industrial-applications.md)).
- **SU2 meshes and one interface for every mesh format (D187, gallery 127,
  recipe 127).** `mesh.file` reads SU2 native files (`.su2`, ASCII) besides
  Gmsh: tetrahedra, hexahedra, prisms and pyramids in 3D; a 2D file
  (triangles, quadrilaterals) is extruded one cell layer in z with the planes
  `zmin` and `zmax` (`mesh.thickness`, default the median edge length).
  Every `MARKER_TAG` becomes a boundary; `[mesh.boundaries]` maps it to the
  name of a boundary condition. The readers sit behind one format layer
  (`foehncraft_mesh::format`): the format comes from the new key
  `mesh.format` (`gmsh`, `su2`), else from the extension, else from the
  content, and every reader hands cells and named faces to one shared mesh
  builder, so a further format is one parser. `foehn mesh` takes an `.su2`
  file directly. The Gmsh reader moved behind the layer unchanged (a test
  compares fingerprints of all Gmsh fixtures). Gallery 127: a NACA 0012
  C-mesh from PyAero, read from SU2, with a laminar flow at Re 1000.
- **GCI study: the level-3 channel drop resolved.** The drop of U_b⁺ on the finest channel meshes (SST 17.346 → 17.009, EB 17.866 → 16.720) came from runs that stopped at the scaled residual 10⁻⁶ while U_b⁺ was still moving (up to 2.9 % from the converged value). `examples/studies/plate_channel_gci.py` now converges the channel to 10⁻⁹ and refines it in y only (same y nodes, two cells in x); converged, SST (GCI 0.12 %, p 1.0) and EB (below 0.01 %, p 2.1) meet the 0.5 % band on levels 1 to 3. Finding: ζ-f₀ is not grid-convergent below first-cell y⁺ 1 (its wall-cell k does not vanish like y²); the plate results are unchanged.
- **OpenFOAM and CGNS meshes (D188, gallery 129, recipe 129).** `mesh.file`
  reads an OpenFOAM case, `constant` or `polyMesh` directory (ASCII or
  binary, polyhedral cells, patches as boundaries, `cellZones` as cell zones)
  and CGNS HDF5 files (hexahedra, tetrahedra, prisms, pyramids, NGON/NFACE
  polyhedra; BCs and Gmsh physical names as boundaries). `mesh.format` takes
  `cgns` and `openfoam`. Decomposed OpenFOAM cases, `.gz` files and ADF
  files are refused with the fix.
- **The viscous weight of the wall treatment as a wall field (spec 2.3,
  D72).** `turbulence-wall` of `[output] derived` and `foehn post derive`
  gains the boundary dataset `blending`: 1 where the wall is resolved
  (`low-re`, the sublayer of `standard`), 0 where the log law holds (the
  log layer, `scalable`), in between 1 − tanh((y⁺/5)²) of `automatic`,
  e^{−Γ} of `compound` and Γ = dU⁺/dy⁺ of `compound-2023`. It shows which
  walls of a mesh are resolved. The restart file stores it as
  `turb_wall_blending`; from an older restart file `turbulence-wall` is
  refused with CASE005 (continue the run for one iteration).
- **Requirements: thin walls (D189) and the heat exchange of the wall film
  (D190), stage 7.** Thin walls on wall boundaries, baffles and region
  interfaces with layers, heat capacity, sources and in-plane shell
  conduction, solved in the one energy matrix; the wall film exchanges heat
  with plain walls, thin walls and solids (evaporation cools the wall,
  condensation heats it), with ledger pairs and per-face slots built in
  stage 7. Checks V7.26–V7.32, PR 7t, open step O24. No code.
- **Requirements: O24 closed.** The maintainer accepted the thin-wall and
  film points TW1–TW5 as recommended (baffles from face zones, thin walls in
  the energy matrix, shell edges on solids refused first, radiation slot,
  film evaporation needs the vapour species of L2); TW6 (a film without
  sprays) waits until the wall film exists. No code.
- **Industry `.msh` and VTK meshes, VTU output (D191, D192, gallery 134,
  recipe 134).** `mesh.file` reads the ASCII face-based industry `.msh`
  (polyhedra, zones; it shares `.msh` with Gmsh, the content decides) and VTK
  grids, legacy `.vtk` and XML `.vtu` in every encoding including zlib.
  `mesh.format` takes `industry_msh` and `vtk`. `output.vtu` and
  `foehn mesh --vtu` write the mesh as a `.vtu` with polyhedra and the names
  of zones and boundaries.
- **Stage 7 plan decided.** The maintainer chose 19 small pull requests in
  four lanes (groundwork, energy, buoyancy, conjugate heat transfer, then
  turbulent heat transfer), GPU kernels in every pull request, water by
  IAPWS tables with all other property forms kept, a solids starter set and
  user material database files (stage 7 page §16.1, D164). No code.
- **Fix: smooth limits keep floors (spec 2.6).** With
  `turbulence.smooth_limits = true` the SST flat plate did not converge:
  the smooth width followed the larger operand, so the floor of CD_kω
  became a large, non-monotone value. The width of `smin`/`smax` now
  follows the smaller operand, m = |a||b|/(|a| + |b|), in a form free of
  cancellation. Runs without smooth limits are unchanged. New test: the
  flat-plate C_f with smooth limits is within 0.090 % of the hard form
  (band 0.1 %).
- **Ledger as data (D119, stage 7 part 7.0a).** The conservation ledger
  books every quantity as an account of named terms (quantity, unit,
  region, term kind, name, exchange-pair id) and closes it with one
  routine, so the energy module and later modules add their quantities
  without new structs. Mass and kinetic energy are booked this way; every
  report (`balance.csv`, log tables, `flow.md`, `result.json`) is bitwise
  as before.
- **Turbulence values of the backflow at pressure outlets (BC028, gallery
  and recipe 140).** A `static-pressure-outlet` or `total-pressure-outlet`
  that gives `k` with `omega`/`epsilon`, `turbulence_intensity` with
  `turbulence_length_scale` or `viscosity_ratio`, or `zeta` now uses them on
  the faces where the flow enters (fixed values), decided per face in every
  outer iteration; the faces where it leaves stay zero gradient. Before,
  the entering fluid always carried the values of the cell and BC028
  warned that the values were unused. New case
  `backstep_outlet_backflow.toml`, plot `tools/plot/outlet_backflow.py`.
- **Fix: restarts of turbulent runs with pressure boundaries.** The
  modified pressure p̂ = p + ⅔ρk prescribed at the pressure boundaries came
  from the start fields until the first turbulence solve after a restart,
  so a continued run was not bit-identical to the uninterrupted one. It is
  now set from the loaded k.
- **Requirements: heat exchangers with several fluids and materials (D194).**
  Any number of fluid and solid regions with their own materials; fluid
  parts joined only through solids or thin walls get their own pressure
  reference and ledger; heat-exchanger outputs (duty, effectiveness, NTU,
  LMTD, UA); thin walls allowed on every face set, also inside solids and
  between two fluids (D189 amended). Checks V7.33–V7.35. No code.
- **Requirements: numbering.** The heat-exchanger and thin-wall decision of
  #142 is D194 (D191 belongs to the `.msh` reader of #134). No code.
- **Iteration limits of the steady cases.** No case file allows 10 000
  iterations any more. `cylinder_re20_gci.toml` runs the coupled algorithm
  (51, 90 and 225 outer iterations instead of 960, 3 007 and 9 876 with
  SIMPLEC at the momentum relaxation 0.5; about a minute instead of 12 on
  four cores) with `max_iterations = 1000`; C_D and Δp agree to 10⁻⁵, the
  lift moves by up to 0.8 % with the pseudo-time step (gallery 25).
  `cavity_re1000.toml` and `cavity_re1000_poly3d.toml` keep SIMPLEC
  (1 397 and 1 240 iterations) with `max_iterations = 2000`;
  `tmr_plate_sst.toml` keeps SIMPLEC on level 1 (635) with
  `max_iterations = 2000`, and its four-level sweep runs with
  `solver.algorithm = "coupled"` (622 and 411 iterations on levels 2 and 3
  instead of 2 292 and more than 4 400).
- **Requirements: models per fluid region.** Turbulence model, buoyancy
  model and Pr_t may be set per fluid region (MF1 decided, D194). No code.
- **CI: shorter waits with parallel work.** Pull requests run on Linux
  only (about 5 minutes); macOS and Windows join when a pull request touches
  build files, the toolchain, the workflow, the command line, MPI, VTKHDF,
  restart, run control or the console. Pushes to main run only what the
  gate did not (macOS and Windows tests and quick cases, MPI on macOS, the
  long cases, debug tests), skip docs-only merges, and a newer merge
  cancels the older run. The full tier runs nightly, by hand and with the
  label `full-ci`. HDF5 file locking is off in CI (flaky on macOS).
- **Banner font.** The name of the code on top of the log is set in the
  "Alpha" figlet font by Lennert Stock ("Foehn" over "craft", 133 columns
  with the box), as in the owner's example; the glyphs are in
  `banner_art.rs` with the credit. Only text changes.
- **Outer iterations in the log.** Every logged time step shows the
  residuals and the linear iterations of each outer iteration under its
  box. New `docs/sayings.md` (*Residuen-gaffen*), linked from the README.
- **Log look (D125).** The log opens with "Foehncraft" in block letters
  (see the banner font entry above). Every block is a box of asterisks as wide as its
  content and centred on the page, which has no fixed width any more
  (the banner is 133 columns). Summaries have blank rows, a heading and
  dashed rules; a logged time step starts with a box (step, time, dt,
  Courant number, outer iterations). Messages read
  `%FOEHN-I-CASE, text`: severity and the capitalised name of the source
  file that writes them. Only text changes; results, CSV and JSON are
  bitwise the same. Gallery 27 updated.
- **Cell-quality remediation (D50).** Steady runs on tetrahedra with skewed
  faces (Gmsh meshes without `Mesh.OptimizeNetgen = 1`) diverged with the
  least-squares pressure gradient. The solver now measures the face
  skewness at start-up; when faces exceed `solver.skewness_limit` (0.9),
  `solver.remediation = "auto"` (the default) switches to the Gauss
  pressure gradient (new key `solver.pressure_gradient`) and the
  Venkatakrishnan limiter, and reports it (`MESH001`, `MESH002`,
  `flow.log`, `flow.md`, `flow.json`). Meshes without such faces run
  bitwise as before (all catalogue cases). The raw Gmsh T-junction that
  diverged in iteration 87 converges; gallery 26, recipe 36, case
  `flow/pipe_gmsh_raw.toml` with its mesh in the repository.
- **Divergence report.** A pressure multigrid that breaks down in the
  middle of a run (the pressure-correction matrix no longer positive
  definite) is reported as a divergence (`RUN014`, `divergence.md`)
  instead of "the pressure level is not fixed"; `DIV003` names skewed
  faces anywhere in the mesh when the remediation is off.
- The mesh quality report counts the faces with a skewness above 0.9.
- `examples/meshes/tjunction_netgen.geo`: the T-junction with the Netgen
  optimisation; `tjunction.geo` loses a refinement field that Gmsh 4.15
  rejects (it never applied, the mesh is unchanged).
- `tools/plot/residuals_chart.py`: residual histories of several runs in
  one chart.

**Log output (D125–D127).** The log of `foehn flow` is rewritten so that a log passed on alone shows the run: a banner with version and build, framed blocks (input and output files, run control, mesh information with quality, fluid properties, solution settings as a matrix per equation, boundary conditions), tables with repeating header rows for iterations and time steps, summaries (convergence, verification, results per boundary, ledger, run time, warnings) and a closing frame with the exit code. Frames use asterisk rules; `[log] style = "box"` draws them with box characters. Console and `flow.log` are identical. Gallery 27, recipe 37.

**Configurable log table (D125, part 2).** New `[log] columns` (a list) chooses the columns of the iteration and time-step table: the residuals, the linear iterations, values at probes (`probe.<name>.<quantity>`) or in one cell (`cell.<index>.<quantity>`), the largest and smallest value of a quantity with the cell where it occurs (`max.velocity_magnitude`, `min.pressure`), the mass flow of a boundary and the time quantities; `[log] digits` sets the significant digits. Without `[log] columns` the table is the one of before, and no result changes with the columns (CPU and GPU), so `[log]` is not part of the case hash. `run.jsonl` and the machine formats go to 1.1 (additive): a `columns` row describes the table, a `table_row` row carries the values of every table line, and a `section` row carries every framed block of the log with its label, value and unit, so a monitor draws what the terminal shows. Gallery 28, recipe 38, case `flow/cavity_re100_columns.toml`.

**Case-file format 0.6 (D126).** New key `case_format` (absent = 0.5); old keys are mapped on reading and listed in the log; `foehn case upgrade [--write|--check]` rewrites files and keeps their comments.

**Renamed: "patch" is now "boundary" (D126).** Keys `[boundary.<boundary>]` stay as they were; the names that change are the file `patches.csv` → `boundaries.csv`, the field `patches` → `boundaries` in `result.json`, `flow.json` and the ledger, the quantity names `patches.*`, the VTKHDF array `patch` → `boundary`, and the words in the log and the reports. Scripts that read the old names must follow.

**Solver settings per equation (D127).** `[solver.momentum]` and `[solver.pressure]` replace the flat keys: `relaxation_velocity` → `momentum.relaxation`, `relaxation_pressure` → `pressure.relaxation`, `momentum_tolerance` → `momentum.tolerance`, `momentum_max_iter` → `momentum.max_iter`, `pressure_tolerance` → `pressure.tolerance`, `pressure_max_iter` → `pressure.max_iter`, `pressure_preconditioner` → `pressure.preconditioner`, `amg` → `pressure.amg`, and `scheme`, `limiter`, `limiter_k`, `central_blending`, `central_bounding` → `momentum.*`. All 73 example case files are upgraded. Old files still run (format 0.5).

**Settings that change during a run (D128).** New table `[[schedule]]`: entries with `at_iteration` (steady), `at_step` or `at_time` (transient) and `set = { key = value }` change settings from a prepared point of the run, so a whole run is set up before it starts (under-relaxation, tolerances, iteration limits, time step, output and log cadences). Schedulable are the keys a re-read may change (D40) without `run.limits.*` and `run.reread_every`; everything else (mesh, boundaries, schemes, fluid) is refused when the case is read, with the reason (`SCHED001`). Every change is a line in `flow.log`, a row in the "Run control" frame and an `event` row `RUN015` in `run.jsonl` (`key`, `old`, `new`); entries that have fired stay in force over a re-read, a restart applies the entries before its start. Case format stays 0.6 (a new optional table). Without `[[schedule]]` results, files and logs are bitwise unchanged. Gallery 29, recipe 39, case `flow/cavity_re100_schedule.toml`.

CPU results are bitwise unchanged against 0.5 (checked with `tools/check/compare_runs.py --across-versions` on the steady and the transient cavity; only the renamed names differ).

## 0.5 (2026-11-16) – stage 5: transient laminar flow, run control and Python

**Highlights.**

- **Transient flow.** PISO with variable-step BDF2 and Crank–Nicolson,
  outer correctors (Courant number 9 at time order 2), an adaptive time
  step from the integration error with exact step rejection, output at
  exact simulated times, potential-flow start, restart files that
  continue a run bit for bit, control files and re-reading the case file
  during a run, probes and force histories. Verified on the Taylor–Green
  vortex (time order 2.0–2.4, space order 2.1) and the cylinder benchmark
  of Schäfer & Turek (Re 20 and Re 100).
- **GPU.** The whole flow loop runs on NVIDIA GPUs from the same kernel
  source as the CPU, bitwise equal to the CPU on a Tesla T4.
- **Adjoint-ready kernels.** Dual-number derivatives of every kernel,
  transposed matrices as gathers, binomial checkpointing (D59).
- **Run contract and health.** `--json`, `--set`, `result.json`, the run
  record, exit codes and signals, the case schema, diagnostic codes, the
  estimate before a run, limits, provenance; the mass and kinetic-energy
  ledger and a divergence report with ranked causes and fixes.
- **Numerics and boundary conditions.** The bounded central scheme, slip
  walls, wedges, inlet profiles from tables, backflow treatment and
  convective outlets; derived fields (vorticity, Q, λ₂, wall shear stress,
  y⁺) with the output and from restart files.
- **Expressions and materials.** Boundary values, body forces and the
  viscosity as expressions of x, y, z and t; a material database with
  sources; regions.
- **Python.** The package `foehncraft`: cases, expression objects, runs in
  process with callbacks, results, materials and parametric studies; all
  scripts of the repository are Python now.
- **Verification.** The grid convergence index of every sweep of three or
  more meshes or time steps, and the generated verification manual.

The parts in detail, newest first:

**Scripts in Python (stage 5 part 6d, D61 (2)).** The shell scripts are
ported to Python and `tools/` is split into groups that run as modules
from the repository root: `tools/check/` (`ci_cases`, `compare_runs`,
`backend`, `c_free`, `kernels`, `vtkhdf`), `tools/bench/` (`build_hypre`
and the cloud GPU drivers with their shared steps in `cloud`),
`tools/mesh/`, `tools/plot/` (the charts), `tools/release/`
(`verification_manual`) and `tools/notebooks/` (the Colab notebooks).
`tools/ci_cases.sh SET` is now `python3 -m tools.check.ci_cases SET`,
`gpu/ptx-kernels/build.sh` is `python3 gpu/ptx-kernels/build.py`; CI,
docs, examples and error messages are updated, and CI lints the scripts
(ruff, syntax and names). See `tools/README.md`.

**Discretisation uncertainty and the verification manual (stage 5 part
6c, D67, D68).** A sweep of three or more meshes or time steps reports
for every monitored quantity (patch mass flows, mean pressures and forces,
force coefficients, benchmark values, viscous and numerical dissipation)
the observed order by the fixed-point iteration of Celik et al. (2008),
the extrapolated value, the band φ1 ± U, the grid convergence index and a
verdict (monotonic, oscillatory, divergent, unchanged): a table in
`flow.md` and `uncertainty` in `flow.json` and `result.json` (an
addition; `Result.uncertainty` in Python). Level sweeps no longer need an
exact solution or a benchmark; `expected_order` without `[analytic]` is
now refused (CASE005) instead of being ignored. New case
`cylinder_re20_gci.toml` (the benchmark 2D-1 on 6 144 to 98 304 cells,
`long` set): C_D extrapolates to 5.5797 ± 0.0004 inside the reference
range. New script `tools/release/verification_manual.py` writes
`docs/verification/manual.md`, `results.json` and a chart of the observed
orders from the case files and a catalogue run. Gallery 25, recipes 34
and 35 on new Verification pages; `tools/gci_chart.py`. Runs are bitwise
unchanged.

**The Python package (stage 5 part 6b, D61, D62).** `python/` holds the
package `foehncraft` with a native module (`crates/foehncraft-py`, PyO3,
built with maturin, abi3 for Python 3.10+): `Case` (read, change with
dotted keys, write; the TOML round trip holds for every example case;
relative paths kept when saved elsewhere), expression objects with Python
operators that render the strings of the case file (`foehncraft.expr`,
evaluated with the solver's own compiler), `run` in process or through
`foehn` with callbacks on every iteration or time step that can stop the
run, `Result` (summary, CSV tables, VTKHDF fields), the material database
and material builders, and `foehncraft.study` (grid, list, Latin
hypercube, Sobol; one case file per variant, `study.csv`, resume). The
engine gains `engine::Observer` (sees every row of the run record) and
`RunOptions::observer`; runs are bitwise unchanged. New CI job "Python
package" (ruff, pytest). New scripts `examples/python/poiseuille_inflow.py`
and `examples/studies/couette_viscosity_study.py`; gallery 24, recipe 33
on the new Python pages, case catalogue page "Python scripts".

**Expressions, materials and regions (stage 5 part 6a, D53, D102, D112,
D121, D124 W3/W4).** Boundary values (velocity, normal velocity, pressure,
total pressure, mass flow), the body force and the viscosity may be
strings of a small expression language in x, y, z and t, with
`[constants]` and piecewise linear `[tables]` (inline or from a file); the
strings are checked when the case is read (EXPR001–EXPR005, with the
place and a "did you mean") and compiled to a stack-machine program that
the solver evaluates with new in-house elementary functions
(`foehncraft_core::math`, D93 (6)), the same on every platform. Fluid
properties from a built-in material database (air, water, engine oil,
ethylene glycol, glycerol, each with its source) with `[fluid] material`,
`temperature` and `operating_pressure`; `[materials.<name>]` with `base`
and the forms constant, polynomial, Sutherland, power law, table (linear
or monotone cubic) and expression; the new command `foehn materials`
lists and evaluates them (also `--json`). `[regions.<name>]` assigns cell
zones to materials; non-fluid regions and unknown zones are refused
(REG001, REG002); different viscosities per region give a viscosity per
face. Walls take a `mass_flux` (blowing and suction), checked for closed
volumes (BC026), and expose per-face exchange coefficients
(`derived::wall_exchange`). The resolved property program is written into
`case.resolved.toml`; `properties_hash` and `materials_database` join the
provenance. The random-number key gains a stream word (`rng::keyed`;
meshes unchanged). `foehn flow --check` also checks the regions against
the mesh. New exact solutions `couette-viscosity` and
`couette-transpiration`; new cases `flow/couette_viscosity.toml`,
`couette_transpiration.toml`, `poiseuille_water.toml`,
`poiseuille_expression.toml`, `taylor_green_expression.toml` and
`channel_ramp.toml`; gallery 23, recipes 31 and 32 on the new
"Expressions and materials" pages; `tools/property_chart.py`. Cases
without the new keys run bitwise as before.

**Derived fields of stage 5 (part 5c, D72).** `[output] derived` adds
fields computed from the solution to every `fields*.vtkhdf`: `vorticity`
(components and magnitude), `q-criterion`, `lambda2`, `enstrophy`,
`strain-rate`, `pressure-coefficient` (reference state in
`[output.reference]`: velocity, pressure, density), `total-pressure`, and
on the patch blocks of the file `wall-shear-stress` and `y-plus`. The new
command `foehn post derive CASE` writes the same file afterwards from a
restart file (`--restart`, default `latest`; `--fields`, `--file`), equal
byte for byte to the file of the run. The VTKHDF writer gains boundary
face fields on the patch blocks; the `units` attribute lists the derived
datasets; the estimate counts them; `output.derived` and
`output.reference` may change during a run. `tools/render_mesh.py` gains
`--diverging`. New case `flow/couette_derived.toml` (exact on any mesh);
gallery 22, recipe 30 on the new Post-processing pages. Without `derived`
every output is as before.

**Boundary conditions of stage 5 (part 5b, D29).** New condition types
`slip-wall`, `wedge` and `convective-outlet` (∂u/∂t + U_c ∂u/∂n = 0 in
transient runs, the outlet flow balanced with the inflow, U_c by default
the mean outflow speed); inlet profiles from a table (`profile = "table"`,
`table`, `table_coordinate` x, y, z or the radius about an axis, the
table's sha256 in the case hash); backflow treatment at static and total
pressure outlets (`backflow = "total-pressure"`, the new default,
`"zero-gradient"`, the treatment before, or `"suppress"`, which closes the
faces where the flow would enter). Axisymmetric wedge meshes from any 2D
mesh of the box generator (`mesh.generate.wedge`, the axis without a
patch), and the exact solutions `uniform` and `pipe-poiseuille`. Seven new
diagnostic codes BC019–BC025. The restart file of a run with a
convective outlet gains the datasets `ub_old_x|y|z`. New cases
`flow/channel_slip.toml`, `poiseuille_table.toml` (with
`poiseuille_profile.csv`), `pipe_wedge.toml`, `cylinder_re100_short.toml`
and `cylinder_re100_convective.toml`; gallery 21, recipe 29. Runs whose
outlets see no backflow are bitwise as before; the Poiseuille cases on
perturbed prisms, polygons and mixed cells, where a few outlet faces see
inflow in the first iterations, change by 2e-6 relative
(`backflow = "zero-gradient"` gives the earlier results bitwise).

**Wall film requirements (D123, D124).** Requirements for the wall film of the Lagrangian sprays (L12): an Eulerian thin film on the wall faces with parcels for splash, rebound, separation and stripping; impingement regimes (Bai–Gosman, Kuhnke), capillarity and dry patches in the first version, CHT coupling of the heated film, ledger and reproducibility rules, staging F1 to F5 and validation cases; Milestone 1 keeps a surface mesh and surface halo of the wall faces, per-face wall boundary values including a wall mass flux, and gravity as a case vector. Documentation only.

**Numerics of stage 5 (part 5a).** A bounded central convection scheme
(`solver.scheme = "central"`, D25) with the bound β and a blend with linear
upwind; the time step from the error estimate of BDF2 and Crank–Nicolson
(`time.error_tolerance`, D65), with steps above the tolerance taken back
exactly and repeated shorter (`solves[].transient.rejected` in
`result.json`, an `error` value in the step rows of `run.jsonl`); output in
simulated time that the steps land on (`output.fields_interval`,
`fields_times`, `monitor_interval`, `run.restart_interval`,
`restart_times`, D66); and a potential-flow start
(`initial.potential_flow`, D52). The restart file gains the attribute
`dt_prev2` and, with the error estimate, the datasets `u_old2_x|y|z`;
older files still read. New cases `flow/taylor_green_error.toml`,
`flow/cylinder_potential.toml` and the variant `bounded-central` of
`scalar/oblique_step_2d.toml`; gallery 20, recipe 28,
`tools/time_step_chart.py`. Runs without the new keys are bitwise as
before.

**Material property requirements (D94–D103, L23–L25).** Requirements for a
material database and the handling of properties: database, case materials
and a resolved, hashed property program; property forms (constant,
polynomial, NASA-7 and NASA-9, Sutherland, tables with linear or monotone
cubic interpolation, polynomial c_p expressions, mixing rules, kinetic
theory, equations of state); checks; GPU evaluation; SI only; no property
changes at run time; a built-in starter set with IAPWS water and steam;
importers, table generators and a generic interface for external property
sources. Implementation starts in stage 5 PR 6 (μ and ρ) and continues with L1.

**Renamed from CFDCraft (and the short-lived Ennecraft) to Foehncraft**, after
the foehn wind of Innsbruck (see the README), so that the name cannot be
confused with other projects. Crates are now `foehncraft-*`, the command is `foehn`, the
environment variables are `FOEHNCRAFT_*`, and the repository is
`chiefenne/Foehncraft`. Entries below keep the names of their time.

**Faster CI.** Pull requests get a fast gate (about 5 minutes instead of
about 25): a `ci` build profile (release optimisation without LTO, results
bitwise equal), the case catalogue split over parallel jobs
(`tools/ci_cases.sh`), and a VTKHDF check that no longer computes unused
volumes of field files (11.5 minutes to seconds). Push to main, a nightly run
at 02:07 Vienna time and pull requests with the label `full-ci` run the full
tier: `release` profile, debug tests, every case including the cylinder at
Re 100, and the slow tests. One summary check, "CI OK". Main warms the `ci`-profile cache that pull
requests read (without it every PR compiled its dependencies from scratch).

**Run contract for scripts and agents (stage 5, part 4a).** `foehn flow`
can be driven by a program: `--set key=value` overrides (TOML values),
global `--json` (one JSON object on standard output, human text on
standard error) and `--quiet`, exit codes 0 finished, 1 error, 2 failed
checks, 3 stopped; SIGUSR1 writes a restart file, SIGTERM, SIGINT and
SIGUSR2 stop cleanly, `--restart latest` continues bitwise; budgets in
`[run.limits]` (wall time, steps or iterations, output volume). Every run
writes `result.json` (status, overrides, provenance with commit and case
hash, diagnostics, results, SHA-256 of the outputs; byte-identical between
identical runs) and the run record `run.jsonl`. 41 stable diagnostic
codes with key paths, locations and `--set` fixes; `foehn schema` prints
the JSON Schema of the case file, `docs/keys.md` is generated from it;
`--check` estimates memory, output and wall time; `--read-root` and
`--write-root` confine the files; CSV, VTKHDF and restart files describe
themselves; HDF5 files have no timestamps. The CLI is a thin layer on
`foehncraft_flow::engine`; `foehncraft-flow` without default features is
checked C-free in CI. Formats and codes: `docs/interfaces.md`; recipe 26,
gallery 18, case `cavity_re100_budget.toml`. No result changes: the
catalogue gives bitwise the fields of main.

**Health checks: conservation ledger and divergence report (stage 5,
part 4b).** Every flow run books mass and kinetic energy (D88): storage,
the terms of each patch (viscous work, convected energy, pressure work,
mass flow), the body force, the viscous dissipation of the discrete
operator, and the remainder, the numerical dissipation, which the second
law wants positive (D90 part 4), split into its pressure part (the
momentum interpolation) and the part of the convection and time
schemes. Steady runs give the
table every `[balance] every` iterations and at the end; transient runs a
compact line with every progress line and a time-integrated ledger of the
run (trapezoidal rule). `balance.csv`, `balance.svg`, `balance` rows in
`run.jsonl`, `solves[].balance` in `result.json`, tolerances with the
codes `BAL001`–`BAL003` (the new key `[balance] tolerance_energy` for
the negative numerical dissipation) and an optional stop. Measured: on the lid-driven
cavity the numerical dissipation is positive, about 1.5 % of the lid
power and falling with the mesh; on the Taylor–Green vortex it falls with
order 2.5, its pressure part with order 3. A run that diverges (D63) now writes its last valid state and
restart file, a VTKHDF snapshot of the failing cells and their
neighbours, and a report (`divergence.md`, `divergence.json`) with the
likely causes ranked by their evidence and a `--set` fix for each, and
ends with `RUN014` (exit code 1). A linear-solver breakdown counts as a
divergence, and a NaN in the continuity residual no longer slips past the
check. The restart file can hold point coordinates and cell volumes
(D91). Recipe 27, gallery 19, case `cavity_re100_energy.toml`,
`tools/balance_chart.py`, `tools/energy_budget_chart.py`. No result
changes: the catalogue gives bitwise the fields of part 4a.

### Adjoint-ready kernels (part 3)

Preparation for a discrete adjoint solver (D59; the solver itself is L17).
No result changes: all flow cases give bitwise the files of part 2.

- **Dual numbers:** `foehncraft_core::Dual<N>` implements `Real`, so every
  kernel body also runs on dual numbers and returns the exact derivatives
  of its results, from the same source. Switches (upwind direction,
  limiter branch, clipping) follow the values, i.e. they are frozen at the
  linearisation state; the limiter array ψ can be held constant.
  **Jacobians of 11 face, cell and boundary kernels agree with central
  differences in every input array (320 580 derivatives, largest relative
  difference 5e-7).**
- **Shape derivatives:** the geometry of `foehncraft-mesh` (face areas and
  centroids, cell volumes and centroids) is generic over the coordinate
  type, so derivatives with respect to the node coordinates come from the
  same code (165 024 derivatives against central differences, 2.8e-7).
- **Transposed matrices:** a transpose map of the matrix storage gives Aᵀ
  in the storage of A and the product Aᵀx without storing it, both gathers
  without atomics (`foehncraft_kernels::adjoint`; flow-device calls
  `transpose`, `matrix_product`, `transposed_product`; PTX entries
  `flow_transpose`, `flow_matrix_product`, `flow_transposed_product`).
  ⟨Av, w⟩ = ⟨v, Aᵀw⟩ to 1e-14 on the matrices of a SIMPLEC solve; CG with
  the multigrid and BiCGStab solve the transposed systems.
- **Binomial checkpointing** (`foehncraft_flow::checkpoint`, Griewank &
  Walther 2000): an optimal schedule (checked against a dynamic program),
  a driver for any primal run and the in-memory snapshot of the flow
  solver; the reverse sweep of a BDF2 run meets every level bitwise.
- Gallery 17, recipe 25, `tools/checkpoint_chart.py`, report
  [stage 5, part 3](docs/verification/stage5.md#adjoint-ready-kernels-part-3).

### GPU flow loop (part 2)

- **The flow loop on the GPU:** `run.backend = "cuda"` or
  `foehn flow --backend cuda` runs every flow case, steady and transient,
  on an NVIDIA GPU (compute capability 7.0 or newer), with all fields,
  matrices and old levels in GPU memory; results bitwise equal to the CPU
  (D31). **Tesla T4 (Colab): all 8 cases bitwise equal to the CPU; 722 ns
  per cell and SIMPLEC iteration on 970 225 polyhedra (7.1 times the two
  CPU threads of the VM), 1 049 ns per cell and PISO step on 98 304 cells
  (9.2 times).** Gallery 16.
- **One source down to the argument lists:** the 27 operations of the flow
  loop are entries of `foehncraft_kernels::entry` (a `#[repr(C)]` argument
  block and a body per index). The CPU loops and the PTX entries `flow_*`
  run the same body with the same block; the flow device is written once
  over an executor (`device::kernel`).
- **CUDA in one crate (D60):** new crate `foehncraft-cuda` holds all CUDA
  code (flow executor, the scalar and multigrid backend, the G1 backends
  and the CUDA C reference, moved from `foehncraft-fv` and
  `foehncraft-bench`); the solver crates find it through registered
  factories. CI lint `tools/check_backend.py`. The feature is now
  `foehncraft-cuda/cuda` (the CLI's `cuda` feature is unchanged).
- New tools `tools/compare_runs.py` (bitwise comparison of two runs),
  `tools/stage5_cloud.py` and `tools/stage5_colab.ipynb` (CPU against GPU
  on eight cases up to one million polyhedra); new cases
  `cavity_poly3d_gpu.toml`, `cylinder_re100_gpu.toml`; recipe 24.
- The CPU results of all flow cases are bitwise those of part 1. A missing
  NVIDIA driver now gives a message instead of a panic.

### Transient laminar flow (part 1)

Stage 5 comes in several pull requests and is released as 0.5 when it is
complete. The first part, the transient core and its benchmarks
([stage 5](docs/verification/stage5.md)):

- **PISO (Issa 1986)** with variable-step BDF2, Crank–Nicolson (first step
  BDF1, the start of Rannacher 1984) and BDF1; two correctors with the
  neighbour correction between them; optional outer iterations per step
  (PIMPLE-type) with a residual tolerance; the convecting flux and the
  deferred corrections extrapolated to the new time, which keeps plain
  PISO second order. A case is transient when it has a `[time]` table.
- **Adaptive time step** from the largest Courant number, with a growth
  limit inside the zero-stability range of BDF2.
- **Rhie–Chow time-step term (D27)** whose steady state is exactly the
  steady solution, for any Δt (cavity: equal to SIMPLEC to 2e-12).
- **Restart files (D42, HDF5):** a continued run is bit-identical to the
  run without the interruption. **Control files and re-reading (D40):**
  `.stop`, `.result_N`, `.read_toml_N`, `run.reread_every`. **Monitors
  (D43):** `history.csv`, `forces.csv`, `probes.csv`, field series,
  force coefficients.
- **Built-in cylinder mesh** `[mesh.cylinder]` (O-grid and wake block) and
  the parabolic inlet profile; the Schäfer & Turek benchmarks 2D-1 and
  2D-2 with their evaluation (C_D, C_L, St, Δp).
- **Results:** Taylor–Green time order 2.21 (BDF2) and 2.00
  (Crank–Nicolson) on the finest pair, space order 2.09, order 2 at
  Courant numbers up to 9 with outer iterations; cylinder Re 20 inside all
  reference ranges on 98 304 cells; Re 100 on 6 144 cells with C_D,max,
  St and Δp inside the ranges, on 24 576 cells C_L,max 0.5 % below its
  range and the others inside.
- `cfdcraft flow --restart FILE --end-time T`; seven new case files,
  gallery 14–15, recipes 21–23. The Release workflow now refreshes the
  cached README badges (recipe 6).

## 0.4 (2026-11-15) – stage 4: laminar steady flow

New crate `cfdcraft-flow`, kernel bodies `cfdcraft-kernels::flow` (one
source for the CPU and the GPU) and the command `cfdcraft flow`:

- **SIMPLEC and SIMPLE (D19).** Segregated outer iterations on the
  collocated, cell-centred polyhedral grid: momentum with implicit upwind
  and deferred linear-upwind convection, orthogonal diffusion implicit and
  the non-orthogonal and skewness corrections deferred, implicit
  under-relaxation; pressure correction with the SIMPLEC factor
  D^C = V/(a_P/α − Σ|a_N|), solved with CG and the stage-3 multigrid
  (hierarchy rebuilt every 10 outer iterations); reference cell when no
  boundary fixes the pressure.
- **Rhie–Chow interpolation (D27)** with the relaxation term of Majumdar
  (1988): the converged result does not depend on the under-relaxation
  (lid force equal to 3e-11 for α_u = 0.7 and 0.9).
- **Boundary conditions (D29):** wall (also moving), velocity inlet (vector,
  normal speed or the exact profile), mass-flow inlet, static- and
  total-pressure inlet and outlet (freely combined), symmetry (implicit per
  component). Periodic boundaries follow with stage 6.
- **Boundary-condition check (D30)** before every run and with
  `cfdcraft flow --check`: missing and unknown patches with "did you mean",
  net mass flow without a pressure boundary, over-specified flow paths,
  wall normal velocity, unused turbulence values, non-planar symmetry,
  reference-cell note; the guidance table of sensible combinations in the
  case catalogue.
- **Case file (D39):** `[fluid]`, `[flow]`, `[initial]`, `[pressure]`,
  `[boundary.<patch>]`, `[solver]`, `[analytic]`, `[verification]`, `[run]`,
  `[output]`; unknown keys and values with "did you mean"; the resolved case
  with all defaults in `case.resolved.toml`.
- **Output (D43):** `flow.log` with a header of what is written and how
  large, `residuals.csv`, `patches.csv` (mass flow, mean pressure, force),
  `wall_<patch>.csv` (pressure and wall shear stress per face, in global
  face order), `line_<name>.csv`, `fields.vtkhdf`, `flow.md`, `flow.json`
  and SVG charts (residuals, convergence, centre lines against the
  reference).
- **Verification:** exact solutions (Couette, Poiseuille) and the
  centre-line data of Ghia et al. (1982) for Re 100 and 1000; thirteen case
  files (CI runs all of them on Linux and the quick ones on macOS and
  Windows); tests for Couette, Poiseuille order, the cavity, relaxation
  independence, quick convergence with symmetry planes, the check, "did you
  mean" and bitwise equality on 1 and 4 threads.
- `Backend::set_matrix` (new matrix values on the same pattern) on the CPU
  and CUDA backends; `tools/render_mesh.py --no-edges`.
- Documentation: [stage 4 verification](docs/verification/stage4.md),
  gallery 12–13, recipes 18–20, flow cases, keys, boundary conditions and
  guidance in the case catalogue.

---

## 0.3 (2026-11-14) – stage 3: pressure solver

New module `cfdcraft-fv::amg`, kernel bodies `cfdcraft-kernels::amg` (one
source for the CPU and the GPU, float32 and float64) and the command
`cfdcraft poisson`:

- **Smoothed-aggregation multigrid (D35).** Strength of connection,
  three-phase greedy aggregation in the cell order, filtered and
  Jacobi-smoothed prolongator, Galerkin coarse operators, explicit
  restriction (gather, no atomics), dense coarsest level; deterministic host
  setup in float64. Plain aggregation as an option.
- **Symmetric V-cycle** with a Chebyshev (default) or damped Jacobi smoother
  in D⁻¹A, λ_max from a fixed power iteration.
- **Mixed precision (D37).** The hierarchy (including a float32 copy of the
  fine matrix) in float32, CG in float64; `precision = "double"` for
  comparison. R10 allows f32 only in the multigrid kernels and rejects
  contracted, approximate and flush-to-zero f32 arithmetic.
- **Flexible CG** (Notay 2000). `[solver] linear_solver = "cg"` and
  `preconditioner = "amg"` in `cfdcraft scalar` (also AMG with BiCGStab);
  the report shows the linear iterations and the multigrid setup time.
- **GPU.** PTX entries `amg_*_f32`, `amg_*_f64`, `blas_axpy`, `blas_xpby`;
  the CUDA backend runs the V-cycle with all vectors resident.
- **hypre reference (D20).** Cargo feature `hypre`: a C shim over hypre's
  public interface (BoomerAMG with PCG), `tools/build_hypre.sh` builds hypre
  from source (CPU or CUDA); a CI job builds it and compares.
- **`cfdcraft poisson`.** Mesh sweeps with several solvers and backends;
  report of iterations, setup and solve times, ns per cell and iteration,
  hierarchy, differences to hypre, bitwise comparison across backends;
  checks of iteration growth and agreement.
- **Verification.** Stage-3 tests (iterations nearly independent of the
  mesh size on all mesh types, mixed precision equals float64,
  SPD V-cycle, bitwise reproducibility, agreement with hypre), five Poisson
  case files and one scalar case, gallery 10–11, recipes 15–17, the
  explanation page `docs/verification/stage3.md`, `tools/stage3_colab.ipynb`.

## 0.2 (2026-11-14) – stage 2: operators on a scalar equation

New crate `cfdcraft-fv` (features `hdf5`, default on, and `cuda`) and new
kernel bodies in `cfdcraft-kernels` (`fv`, `blas`), one source for the CPU
and the GPU:

- **Least-squares gradient (D24).** Weighted least squares over the face
  neighbours; weights per cell-face entry precomputed on the host;
  fixed-value, fixed-gradient and symmetry boundaries keep linear fields
  exact on every mesh type.
- **Convection (D25).** First-order upwind and linear upwind as a deferred
  correction; Venkatakrishnan and Venkatakrishnan–Wang limiters.
- **Diffusion (D26).** Over-relaxed non-orthogonal correction, optionally
  limited.
- **Assembly (D31, D37).** Face kernel for the flux coefficients, cell
  kernel that gathers them into a sliced-ELL row (C = 32, σ = 1), the
  right-hand side and the Jacobi diagonal; implicit under-relaxation. No
  atomics; results independent of the thread count.
- **BiCGStab (D36)** with Jacobi preconditioning; dot products in fixed
  blocks of 256, the same order on the CPU and the GPU.
- **Backends.** CPU (rayon) and NVIDIA GPU (`cuda`: the bodies as PTX entry
  points `fv_*`, `blas_*` in `gpu/ptx-kernels`, all fields resident on the
  GPU).
- **`cfdcraft scalar`.** Case file with `[problem]` (manufactured solutions,
  boundary profiles), `[solver]`, `[run]` (backends compared with each
  other), `[verification]` (refinement sweeps, expected order, bounds),
  `[[variants]]` and `[output]` (VTKHDF fields, sampling line). Reports in
  Markdown and JSON with time per phase and ns per cell and iteration (R6),
  convergence and profile plots as SVG.
- **Verification.** Tests for exact gradients, no new extrema, order of
  accuracy 1.8–2.2 on hex, tet, poly and mixed meshes in 2D and 3D
  (release builds, R7), bitwise reproducibility; seven scalar case files,
  gallery examples 8 and 9, recipes 12–14, the explanation page
  `docs/verification/stage2.md`.
- **CI.** Every scalar case except the GPU one runs as a test (R5);
  `tools/check_kernels.sh` refuses atomics, f32 and fused or approximate f64
  arithmetic in the kernels and the PTX (R10).
- **GPU check.** `tools/stage2_colab.ipynb` and `tools/stage2_cloud.sh`
  (CPU against two GPU runs on 0.9 million polyhedra). Result on a Tesla
  T4 (Colab): both GPU runs, the four oblique-step schemes and the sweep
  are bitwise equal to the CPU; 7.1 s against 158.8 s on the 2-thread
  Colab VM.

---

## 0.1 (2026-11-14) – stage 1: mesh

- Release workflow (`.github/workflows/release.yml`): creates the tag
  `v<version>` on `main` and a GitHub release with the changelog section of
  that version as notes; dry-run by default (recipe 6).

New crate `cfdcraft-mesh` (features `hdf5` and `metis`, both on by default;
`mpi` for the MPI halo test):

- **Mesh structure (D23).** Faces as node lists, `owner` and `neighbour`
  per face, the faces of each cell with their orientation; patches and cell
  zones. Order: inner faces, processor faces, boundary faces by patch;
  cells inner, owned cells touching the halo, halo cells by source rank,
  each group in reverse Cuthill-McKee order (R1). u64 global ids for every
  cell, face and node in input order (R2). `validate()` checks index ranges,
  ordering and topological closure (every edge used twice).
- **Geometry.** Face area vectors and centroids from a triangle fan about
  the node average, cell volumes and centroids from tetrahedra on the same
  fan: Σ S = 0 to round-off for every cell, volumes add up to the domain
  volume.
- **Quality report.** Closure, total against enclosed volume, patch areas,
  non-orthogonality, skewness, face warp, aspect ratio, volume ratio,
  inverted faces, non-positive volumes; Markdown and JSON.
- **Gmsh reader (D33).** MSH 4.1 ASCII and binary (either endianness),
  tetrahedra, hexahedra, prisms, pyramids; physical surfaces as patches,
  physical volumes as cell zones, periodic links; refuses higher-order
  elements and other versions with a clear message.
- **Test mesher (D34).** Box meshes of hexahedra, tetrahedra (Kuhn split),
  polyhedra (dual of the tetrahedra), mixed cells; 3D or 2D one cell thick
  (D3); exponential grading per direction, one- or two-sided; nested
  refinement levels; node perturbation by the counter-based generator,
  boundary nodes kept on their faces. Conversion of any tetrahedral mesh into
  polyhedra with the exact volume (also on curved walls).
- **Renumbering.** Reverse Cuthill-McKee and shuffle, with the ids kept.
- **Partitioning (D32).** METIS (vendored) and recursive coordinate
  bisection; imbalance, cut faces and neighbours in the report.
- **Decomposition and halo (D38, R3).** Local meshes with halo depth as a
  parameter; per-neighbour send lists and contiguous receive ranges;
  several fields per message; host-staged transport now, GPU-aware selected
  at run time later on the same plan. Halo geometry filled by exchange.
- **VTKHDF writer (D41).** VTKHDF 2.4 `PartitionedDataSetCollection`:
  block `internal` with cell fields and one block per patch; standard VTK
  cells and `VTK_POLYHEDRON`; global-id order. HDF5 built from source and
  linked statically (needs CMake).
- **Point location and tracking (D47).** Tetrahedral-decomposition walk
  from cell to cell through the faces, with a bin search where the walk
  stops at a non-convex boundary; reports where a track leaves the mesh or
  the rank.

Tests: closure, volume and patch areas of 16 generated variants; Gmsh ASCII
against binary, zones, prisms, pyramids, periodic links, refusals; VTKHDF of
all types; point location on all types; decomposition with METIS and RCB on
1, 2, 3 and 5 parts and 1 and 2 halo layers; **golden values (R4)** for 32
cases of mesh type × dimension × perturbation × refinement × 1, 2 and 4
ranks, checking that owned volumes sum to the total and that every halo
cell has exactly one owner; the MPI halo exchange on 2 and 4 ranks in CI.

CLI: `cfdcraft mesh CASE.toml | FILE.msh` builds or reads, renumbers,
checks, partitions and writes the report and `mesh.vtkhdf`. G1:
`mesh.type = "poly" | "tet" | "hex" | "mixed"` runs the benchmark on
stage-1 meshes.

Tools: `tools/check_vtkhdf.py` (checks files with VTK's reader),
`tools/render_mesh.py` (gallery pictures), `tools/gmsh_test_meshes.py`
(the Gmsh test files), `examples/meshes/tjunction.geo`.

Documentation: gallery examples 3–7, mesh cases and keys in the case
catalogue, recipes 7–11, the G1 stage-1 rerun, implementation notes in
D23, D32, D38, D41 and D47, new open item O15.

Findings:

- Polyhedral boundary faces merged across a curved wall lose volume
  (0.34 % on the T-junction); keeping one planar piece per triangle there
  makes the dual exactly as large as the tetrahedral mesh.
- VTK computes the volume of a polyhedron as that of the convex hull of its
  nodes, so the VTKHDF check computes it from VTK's faces instead.
- **G1 on polyhedra (shared cloud machine, 2 million cells):** bitwise
  equal; the single-source bodies are 25–42 % faster on cell gather, equal
  on SpMV-CSR, 4–30 % slower on face flux, as on a hexahedral box of the
  same size. Shuffled order costs a factor 3 to 4. Quiet-machine and GPU
  reruns are open (O15).

---

## 0.0 (2026-11-14) – stage 0: foundations and gate G1

Workspace (`Cargo.toml`, Rust 2024 edition, minimum Rust 1.89):

- `cfdcraft-core` (`no_std`): the generic floating-point type `Real`
  (f64 reference path, f32 possible, D7), 32-bit local and 64-bit global
  indices, pairwise summation in a fixed order (D31), counter-based random
  numbers keyed by identity, step and event (D47).
- `cfdcraft-kernels` (`no_std`): the kernel bodies, one source for CPU and
  GPU: face flux, cell gather, SpMV in CSR and SELL-C-σ. Gather-based, no
  atomics. Feature `checked-indexing` for bounds-checked loads.
- `cfdcraft-parallel`: `Communicator` trait with `SerialComm` and `MpiComm`
  (feature `mpi`, system MPI library, D9); reproducible all-rank sums;
  multi-threaded sums bitwise equal to the sequential sum for any number of
  threads.
- `cfdcraft-bench`: gate G1. Synthetic periodic face-based mesh (natural or
  shuffled numbering), CSR and SELL-C-σ matrices, backends `cpu-serial`,
  `cpu-threads`, `cpu-native`, and with feature `cuda` `cuda-reference`
  (hand-written CUDA C, compiled by NVRTC at run time) and `cuda-rust` (the
  kernel bodies compiled to PTX). TOML case files, JSON, Markdown and SVG
  reports, G1 verdict per candidate/reference pair.
- `cfdcraft-cli`: the `cfdcraft` command with `info`, `bench g1` and
  `bench backends`.
- `gpu/ptx-kernels`: the kernel bodies as PTX entry points, built with the
  nightly toolchain (`build.sh`).

Continuous integration (`.github/workflows/ci.yml`, D22): format, clippy and
docs; tests in debug and release on Linux, macOS and Windows; minimum Rust
version; MPI smoke test with 2 and 4 ranks on Linux and macOS (O12); PTX
build of the single-source kernels; G1 smoke run on every platform.

Cloud GPU runs: `tools/g1_cloud.sh` (installs Rust, builds PTX and CLI for
the GPU found, runs G1 twice) and `tools/g1_colab.ipynb` (the same on a free
Google Colab T4).

CLI: `cfdcraft info` reports the communicator and the rank count.

Documentation: [gallery](docs/gallery.md), [case catalogue](docs/cases.md),
[recipes](docs/recipes.md), [G1 procedure and results](docs/benchmarks/g1.md).

Findings:

- **G1 GPU verdict, Tesla T4 (Google Colab), 10 million cells: PASS.** The
  single-source Rust kernels, compiled to PTX, run within 1.0 % of the
  hand-written CUDA C kernels on all four kernels (220–255 GB/s of 320 GB/s),
  in two runs, with results bitwise equal to the CPU. D8: one kernel source;
  D11 stands. The single-source GPU route needs nightly Rust (O14).
- **G1 CPU verdict, Intel i9-9900K (macOS), 10 million cells: PASS** in two
  runs: largest time ratio single-source against hand-written Rust loops
  1.018, results bitwise equal, 20–26 GB/s.
- The rustc NVPTX target supports compute capability 7.0 (Volta) and newer
  only; G1 needs such a GPU for the single-source route (O13). The
  maintainer's GTX 960 (compute capability 5.2) runs the CPU verdict and the
  CUDA reference (`g1_gpu_reference_4m.toml`, recipe 2); the GPU verdict is
  taken on a free Colab T4 (recipe 3).
- On the CPU the single-source bodies match hand-written loops once the
  parallel loop uses large enough work packets (64 Ki outputs): with 4 Ki the
  face-flux kernel was about 30 % slower on 10 million cells.
- A branching max/min in the upwind flux cost bandwidth on random mass-flux
  signs; `Real::max_of`/`min_of` are now branch-free.
