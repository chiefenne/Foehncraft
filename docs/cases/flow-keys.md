# Flow case-file keys

Part of the [Case catalogue](../cases.md): every key of a flow case, the boundary conditions (D29) and their sensible combinations (D30).

`[mesh]`, `[mesh.generate]` as in the [mesh cases](meshes.md#mesh-case-file-keys).
Unknown keys and values are errors with a "did you mean" hint (D39). SI
units. The complete list, generated from the code with types, defaults,
units and whether a key may change during a run, is
[keys.md](../keys.md); `--set key=value` overrides any key from the command
line (recipe 26).

| Key | Default | Meaning |
|---|---|---|
| `fluid.density` | 1 | ρ (kg/m³), or from `fluid.material` |
| `fluid.viscosity` | 1e-3 | dynamic viscosity μ (Pa s), or from `fluid.material`; an expression of x, y, z gives a viscosity that varies in space |
| `fluid.material` | none | a material of `[materials]` or of the database (`foehn materials`), below |
| `fluid.temperature` | 293.15 | the temperature at which the material properties are evaluated (K; isothermal until the energy equation) |
| `fluid.operating_pressure` | 101325 | the pressure at which they are evaluated (Pa) |
| `flow.convection` | true | false: Stokes flow |
| `flow.body_force` | `[0, 0, 0]` | body force per unit volume (N/m³); each component a number or an expression of x, y, z, t |
| `gravity.vector` | none | the gravitational acceleration g (m/s²); with it the solved pressure is the reduced pressure p′ = p − ρ_ref g·(x − x_ref) and the output has `static_pressure` p as well (D153, [recipe 186](../recipes/energy.md#buoyancy)) |
| `gravity.reference_point` | `[0, 0, 0]` | x_ref of the reduced pressure (m) |
| `buoyancy.model` | `"none"` | `"boussinesq"`: the body force [(ρ₀ − ρ_ref) − ρ₀ β (T − T_ref)] g with the balanced force of D154; needs `[gravity]` and `[energy]`. `"variable-density"`: the body force (ρ − ρ_ref) g with ρ(T, p₀) in every cell, and the variable density of `[energy] density = "variable"` with it (stage 7c, D151, [recipe 195](../recipes/energy.md#variable-density)) |
| `buoyancy.expansion_coefficient` | from the material | β (1/K); default −(1/ρ) dρ/dT of the fluid's density at T_ref (1/T_ref for an ideal gas) |
| `buoyancy.reference_temperature` | `fluid.temperature` | T_ref (K) |
| `buoyancy.reference_density` | ρ₀ | ρ_ref (kg/m³) of the reduced pressure |
| `flow.frame_rotation` | `[0, 0, 0]` | angular velocity Ω (rad/s) of a frame in which the whole domain is solved: the velocity is the one relative to the frame, the momentum equations take the Coriolis force −2ρ Ω × u, the centrifugal force stays in the pressure (p − ½ρ\|Ω × r\|² is solved); for periodic flows in system rotation (the rotating channel of 9.8); not with `[[mrf]]` or `[[sliding]]` |
| `initial.velocity` | `[0, 0, 0]` | start velocity (m/s) |
| `initial.pressure` | 0 | start pressure (Pa) |
| `initial.potential_flow` | false | start from potential flow instead of `initial.velocity` (D52): the irrotational velocity that carries the inflow to the pressure boundaries, one Laplace solve with the pressure multigrid and four non-orthogonal correction passes |
| `pressure.reference_point` | lower corner of the bounding box | the cell containing it holds the pressure level when no boundary fixes it |
| `pressure.reference_value` | 0 | pressure in that cell (Pa) |
| `boundary.<boundary>` | required for every boundary | condition, table below |
| `solver.algorithm` | `"simplec"` | `simplec` (Van Doormaal & Raithby 1984), `simple`, or `coupled`: momentum and continuity of an outer iteration as one linear system, with pseudo-time steps (D176 to D178, [recipe 78](../recipes/flow.md#coupled)) |
| `solver.coupled.tolerance`, `solver.coupled.max_iter` | 0.1, 30 | `coupled`: relative residual of the coupled system per outer iteration and the largest number of Krylov iterations (flexible GMRES) |
| `solver.coupled.momentum_sweep` | `"multigrid"` | `coupled`: the momentum sweep of the preconditioner, per velocity component: `multigrid`, V-cycles of a reduction multigrid (D179, [recipe 80](../recipes/flow.md#momentum-multigrid)), or `bicgstab` with a diagonal preconditioner |
| `solver.coupled.momentum_cycles`, `solver.coupled.momentum_reuse` | 1, 0.02 | `multigrid` sweep: cycles per sweep; a hierarchy is kept while the momentum matrix has changed by at most this part since its setup (0: a setup in every outer iteration) |
| `solver.coupled.multigrid.strength`, `.restriction_strength`, `.drop_tolerance`, `.max_coarse`, `.max_levels`, `.precision` | 0.3, 0.05, 10⁻⁴, 64, 30, `"mixed"` | the multigrid of the momentum sweeps: strength thresholds of the splitting and of the restriction, the part of the diagonal below which coarse entries are dropped, the size of the coarsest level, the number of levels, float32 (`mixed`) or float64 (`double`) hierarchy |
| `solver.coupled.momentum_tolerance`, `solver.coupled.momentum_max_iter` | 0.01, 10 | `bicgstab` sweep: tolerance and largest number of iterations per component |
| `solver.coupled.pressure_tolerance`, `solver.coupled.pressure_max_iter` | 0.3, 20 | `coupled`: the pressure solve of the preconditioner (CG with the preconditioner of `[solver.pressure]`); the Krylov iteration corrects what it leaves, so 0.3 costs no Krylov iterations (0.1 before #82) |
| `solver.coupled.pressure_reuse` | 0.05 | `coupled`: the multigrid hierarchy of the pressure solve is built anew when its matrix has changed by more than this fraction since its setup (the measure of `momentum_reuse`), besides every `solver.amg_rebuild` iterations; 0: whenever it has changed; 1: only every `amg_rebuild` iterations (before #83) |
| `solver.pseudo_time.step` | automatic | pseudo-time step Δτ of a steady run (s), the same in every cell; automatic: `factor` × L/U, L the volume of the domain over the area of its boundaries without symmetry planes, U the larger of the largest prescribed boundary speed (inlets, moving walls) and the root mean square of the cell speeds |
| `solver.pseudo_time.factor` | 1 | factor on the automatic step |
| `solver.pseudo_time.momentum` | true | the momentum equations of `coupled` take the step; false: relaxation only |
| `solver.pseudo_time.turbulence` | 0.5 with `coupled`, else 0 | c_T: the turbulence equations take the step min(Δτ, c_T T), T the time scale of the model; 0: relaxation only |
| `solver.momentum.scheme` | `"linear-upwind"` | convection of momentum: `upwind`, `linear-upwind` or `central` (bounded central, D25) |
| `solver.momentum.central_blending`, `solver.momentum.central_bounding` | 1, 0.1 | with `central`: share of the central correction (rest linear upwind) and the bound β (as for the scalar cases) |
| `solver.momentum.limiter` | `"none"` | gradient limiter: `none`, `venkatakrishnan`, `venkatakrishnan-wang`, `barth-jespersen` (`mlp` is for the scalars and the energy equation only) |
| `solver.momentum.limiter_k` | 5 | constant K of the limiter |
| `solver.nonorth_correction` | true | deferred non-orthogonal and skewness correction of the viscous term |
| `solver.momentum.relaxation` | 0.7; 0.98 with `coupled`, 0.995 for the k-ε family with `solver.turbulence.preconditioner = "reduction"` (D180) | under-relaxation α_u of the momentum equations (implicit, Patankar 1980); with `coupled` it only caps the local Courant number of the pseudo-time step at α/(1 − α) |
| `solver.pressure.relaxation` | 1 (SIMPLEC), 0.3 (SIMPLE) | relaxation α_p of the pressure correction |
| `solver.pressure_gradient` | `"least-squares"` | gradient of the pressure in the momentum source, the Rhie–Chow flux and the wall pressure: `least-squares` (second order on every mesh) or `gauss` (the face sum of the face pressures; robust on skewed faces) |
| `solver.face_velocity` | `"interpolated"` | face velocity of the Rhie–Chow mass flux: `interpolated` (w u_P + (1 − w) u_N moved to the face centroid with the interpolated gradient) or `bounded` (the mean of the two limited extrapolations to the face centroid; exact for a linear field, and on a skewed face it stays between the values the two cells see) |
| `solver.transposed_stress` | `"full"` | turbulent flow: the explicit transposed part of the viscous stress, `full` or `trace-free` (without −μ (∇·u) I: the same term where ∇·u = 0, without the grad-div part that cell gradients on tetrahedra leave in it; a further safety net for poor meshes, not switched by the remediation) |
| `solver.remediation` | `"auto"` | cell-quality remediation (D50): `auto` switches to the `gauss` pressure gradient, the Venkatakrishnan limiter and the `bounded` face velocity when a face of the mesh is skewed beyond `skewness_limit`, and reports it (`MESH001`, `MESH002`, the report); `off` keeps the schemes as set. A mesh without such faces runs unchanged ([recipe 36](../recipes/meshes.md#skewed-tets)) |
| `solver.skewness_limit` | 0.9 | faces with a larger skewness (distance of the face centroid from the line through the two cell centroids, over the centroid distance) are poor |
| `solver.max_iterations` | 2000 | largest number of outer iterations |
| `solver.tolerance` | 1e-8 | converged when every scaled residual (momentum, continuity) is below it; on fine wall-resolved meshes an integral quantity can still move after 10⁻⁶ is reached (the channel on 768 cells across: U_b⁺ 2.9 % from its converged value at 10⁻⁶, 0.001 % at 10⁻⁹; [verification](../verification/stage6.md#the-level-3-channel-drop-pr-128)), so check the monitored quantity or tighten the tolerance |
| `solver.momentum.tolerance` | 1e-2 steady, 1e-8 transient | relative tolerance of each momentum solve (BiCGStab, Jacobi) |
| `solver.momentum.max_iter` | 50 | largest number of BiCGStab iterations per solve |
| `solver.pressure.tolerance` | 1e-2 steady, 1e-8 transient | relative tolerance of each pressure-correction solve |
| `solver.pressure.max_iter` | 100 | largest number of CG iterations per solve |
| `solver.pressure.preconditioner` | `"amg"` | `amg` (stage 3) or `jacobi` |
| `solver.pressure.amg.*` | as `linear.amg.*` of the [Poisson cases](poisson.md#poisson-case-file-keys) | multigrid settings |
| `solver.amg_rebuild` | 10 | rebuild the hierarchy every N outer iterations |
| `case_format` | current (`"0.6"`; absent = `"0.5"`) | version of the case-file format; keys renamed later are mapped on reading and logged, `foehn case upgrade` rewrites the file |
| `log.style` | `"asterisk"` | frames of the log: `asterisk` (horizontal rules only) or `box` (drawn with box characters); console and `flow.log` are identical, without colour |
| `solver.log_every` | 10 | a log line every N outer iterations |
| `schedule[].at_iteration` | none | steady: the iteration from whose start the values of `set` hold (one trigger per entry, entries in ascending order, D128) |
| `schedule[].at_step`, `schedule[].at_time` | none | transient: the step or the time (s) from whose start the values hold; the run lands a step on `at_time` |
| `schedule[].set` | none | the settings that change, dotted keys of the case file (`solver.momentum.relaxation = 0.9`); only keys a re-read may change, without `run.limits.*` and `run.reread_every` |
| `analytic.type` | none | exact solution for the checks: `uniform`, `couette`, `poiseuille` or `pipe-poiseuille` |
| `analytic.lower`, `analytic.upper` | 0, 1 | channel walls (y) |
| `analytic.wall_velocity` | 1 | Couette: speed of the upper wall |
| `analytic.mean_velocity` | 1 | Poiseuille: mean velocity |
| `analytic.x_outlet`, `analytic.pressure` / `pressure_outlet` | 1, 0 | Poiseuille: outlet position and pressure; Couette: the constant pressure |
| `verification.levels` | none | refinement levels of `[mesh.generate]` or `[mesh.cylinder]`: a sweep; with `[analytic]` the observed order of the errors, with three or more levels the discretisation uncertainty of the monitored quantities (D67, [recipe 34](../recipes/verification.md#gci)) |
| `verification.expected_order` | none | accepted range of the observed order on the finest pair (needs `[analytic]`) |
| `verification.max_error` | none | largest accepted relative L2 velocity error |
| `verification.reference` | none | `ghia-re100` or `ghia-re1000`: centre-line data of the unit cavity |
| `verification.reference_tolerance` | none | largest accepted deviation, in lid speeds |
| `verification.lid_velocity` | 1 | lid speed for the reference comparison and the Reynolds number |
| `verification.require_convergence` | none | fail unless the outer iterations converge |
| `run.threads` | all cores | CPU threads (meshes below 20 000 cells run on one; the results are bitwise the same) |
| `run.backend` | `cpu` | where the flow loop runs: `cpu`, or `cuda` in a build with the `cuda` feature (NVIDIA GPU of compute capability 7.0 or newer); the results are bitwise those of `cpu` (D31). `foehn flow --backend` overrides it |
| `run.ptx` | `target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx` | PTX module of the `cuda` backend (`python3 gpu/ptx-kernels/build.py`); `foehn flow --ptx` overrides it |
| `output.directory` | `output/<case name>` | output directory |
| `output.vtkhdf` | true | write `fields.vtkhdf` |
| `output.history_every` | 1 | a row of `residuals.csv` every N outer iterations |
| `output.wall_csv` | true | write `wall_<boundary>.csv` for every wall |
| `output.lines` | none | sampling lines `[[output.lines]]` with `name`, `from`, `to`, `points` (101) |
| `output.fields_every` | 0 | transient: write `fields_<step>.vtkhdf` every N steps (0: only the final `fields.vtkhdf`) |
| `output.fields_at` | none | transient: also write the fields at these steps |
| `output.fields_interval`, `output.fields_times` | none | transient: also write the fields at the multiples of this simulated time (s) and at these times; the steps land on them exactly (D66) |
| `output.monitor_every` | 1 | a row of `forces.csv` and `probes.csv` every N steps (`probes.csv` has `<probe>_T`, the temperature of the probe's cell, with the energy equation; in a conjugate case the probes are found on the whole mesh, solids included) |
| `output.monitor_interval` | none | transient: a row of `forces.csv` and `probes.csv` at the multiples of this simulated time (s) instead of every `monitor_every` steps, landed on exactly (equally spaced samples for spectra, D66) |
| `output.probes` | none | point probes `[[output.probes]]` with `name` and `point`: velocity and pressure in `probes.csv`; with a turning zone also the velocity relative to it, `<name>_u_rel`, `_v_rel`, `_w_rel` (D183) |
| `output.wall_probes` | none | wall probes `[[output.wall_probes]]` with `name`, `boundary`, `point`: wall shear stress, C_f (with `output.reference.velocity`), y⁺ and the wall value and flux of each passive scalar (the Nusselt number of a heated wall, 6g) at a point of a wall, in the report and the uncertainty table of a level sweep |
| `output.wall_crossings` | none | wall crossings `[[output.wall_crossings]]` with `name`, `boundary`, `from`, `direction` (x), `kind` (`reattachment`: the wall shear along `direction` turns from negative to positive; `separation`: from positive to negative), `width` (only the faces within this distance of the line; default all), `skip` (leave out the faces closer than this to `from`, e.g. the corner vortices behind a step; default 0): the first such point after `from`, linear between the two face centroids, as the distance from `from` (the reattachment length of a step), in the report and with its GCI in a level sweep (6g) |
| `output.coefficients` | none | force coefficients C = 2F/(ρU²A) of the sum over `boundaries`, with `velocity` U, `area` A, `drag_direction` (x), `lift_direction` (y) |
| `output.derived` | none | derived fields in every `fields*.vtkhdf` (D72): `vorticity`, `q-criterion`, `lambda2`, `enstrophy`, `strain-rate`, `pressure-coefficient`, `total-pressure`, `relative-velocity`, `frame-velocity`, `turbomachinery`, `relative-total-pressure`, `rothalpy` (with a turning zone, D183), `wall-shear-stress`, `y-plus`, `wall-distance`, `turbulence-features`, `turbulence-scales`, `reynolds-stress`, `wall-normal` (EB), `turbulence-wall` (with a turbulence model), `curvature-correction` (SST with `turbulence.curvature_correction`); see [Derived fields](#derived-fields-d72) |
| `output.statistics` | none | transient: time means, RMS and covariances of the fields, probes and force coefficients from `start` (s), the batch-means confidence interval with `batches` (20) and the phase averages over `phase_period` (s) in `phase_bins` (24) bins, `tolerance` (10⁻³) of the verdict; `fields` (true) writes `statistics.vtkhdf` (D71); see [Statistics](#statistics-of-transient-runs-d71) |
| `output.machines` | none | `[[output.machines]]` tables `name`, `rotor` (a turning zone; default: the only one), `inlet`, `outlet` (lists of boundaries) and `gravity` (9.80665 m/s²): volume flow, total pressure rise, head, hydraulic and total-to-static efficiency and Euler work between inlet and outlet (D184, [recipe 101](../recipes/flow.md#torque-power)) |
| `mrf.walls`, `sliding.walls` | the walls inside the zone | the walls of the rotor, whose torque gives the shaft power −Ω·T (D184) |
| `output.frame` | none | the `zone` of an `[[mrf]]` or `[[sliding]]` table whose frame the fields of rotating machinery take in every cell (D183); default: each turning zone in its own frame, the rest in the fixed one |
| `output.reference.velocity`, `.pressure`, `.density` | `[output.coefficients] velocity`, 0, `[fluid] density` | reference state U_ref (m/s), p_ref (Pa), ρ_ref (kg/m³) of `pressure-coefficient` |

Transient cases (stage 5):

| Key | Default | Meaning |
|---|---|---|
| `time.scheme` | `"bdf2"` | `bdf1` (backward Euler), `bdf2` (variable steps) or `crank-nicolson`; the first step of BDF2 and Crank–Nicolson is a BDF1 step |
| `time.time_step` | required | Δt (s); with `adaptive` the first step |
| `time.end_time`, `time.steps` | one of them required | end time (s), number of steps |
| `time.start_time` | 0 | time of the initial fields (s) |
| `time.adaptive` | false | adapt Δt to the largest Courant number `max_cfl` |
| `time.max_cfl` | 1 | target Courant number with `adaptive` |
| `time.max_growth` | 1.2 | largest ratio of successive steps, below 1 + √2 (zero stability of BDF2) |
| `time.min_time_step`, `time.max_time_step` | 0, none | bounds of the adaptive Δt (s) |
| `time.error_tolerance` | none | with `adaptive` and `bdf2` or `crank-nicolson`: relative tolerance of the local time-integration error per step (D65); the step follows from the error estimate, inside `max_cfl` and `max_growth`, and a step above the tolerance is repeated shorter (at most 8 times, not below `min_time_step`). Keep it at least ten times above the linear-solver tolerances |
| `time.error_absolute` | 1e-6 | absolute tolerance of the error estimate (m/s), added to `error_tolerance` × the RMS velocity |
| `solver.algorithm` | `"piso"` with `[time]` | PISO (Issa 1986) |
| `solver.correctors` | 2 | pressure corrections per outer iteration (one diverges on meshes with high diffusion numbers, such as the cylinder) |
| `solver.outer_iterations` | 1 | outer iterations per time step (1: plain PISO; more: PIMPLE-type, for large steps) |
| `solver.outer_tolerance` | 1e-6 | ends the outer iterations of a step when the residuals at the start of one are below it |
| `solver.extrapolation` | true | extrapolate the convecting flux, velocity gradients and limiters to the new time in the first outer iteration (keeps plain PISO second order) |
| `solver.momentum.relaxation` | 1 transient | below 1 needs `outer_iterations` > 1 |
| `run.restart` | none | restart file to continue from (D42) |
| `run.restart_every` | 0 | write `restart.h5` every N steps or iterations |
| `run.restart_at` | none | also write `restart_<N>.h5` at these steps |
| `run.restart_interval`, `run.restart_times` | none | transient: write `restart.h5` at the multiples of this simulated time (s), and a numbered `restart_<N>.h5` at these times; landed on exactly (D66) |
| `run.restart_final` | true | write `restart.h5` at the end |
| `run.reread_every` | 0 | re-read the case file every N steps or iterations (0: only on `.read_toml_N`), D40 |
| `run.limits.max_wall_time` | none | stop after this wall time (s) with the restart file and results, exit code 3 (D83) |
| `run.limits.max_steps`, `max_iterations` | none | stop after this many time steps (transient) or outer iterations (steady) |
| `run.limits.max_output` | none | stop when the output directory holds this many bytes |
| `balance.every` | 0 | the full balance table (mass, kinetic energy per boundary and term) in the log every N iterations or steps; 0: at the end only (D89) |
| `balance.tolerance_step` | 1e-8 | relative mass imbalance of a time step or of the final steady state before the warning `BAL001` |
| `balance.tolerance_run` | 1e-4 | relative mass imbalance over the run (`BAL002`) |
| `balance.tolerance_energy` | 1e-2 | negative numerical dissipation of kinetic energy, relative to the largest term (`BAL003`) |
| `balance.action` | `"warn"` | `"stop"`: a transient run stops at the step with exit code 3 |
| `verification.time_steps` | none | a time-order sweep on one mesh: three or more steps, largest first |
| `verification.expected_time_order` | none | accepted range of the observed time order on the smallest steps |
| `verification.benchmark` | none | `schafer-turek-2d1` (steady, Re 20) or `schafer-turek-2d2` (periodic, Re 100) |
| `verification.benchmark_boundary` | `"cylinder"` | boundary of the cylinder |
| `verification.benchmark_tolerance` | none | largest accepted distance from the reference ranges, relative to their midpoints |
| `verification.evaluation_start` | half the run | start of the evaluation window of 2D-2 (s) |
| `analytic.type = "taylor-green"` | | the decaying Taylor–Green vortex with `velocity` U (1), `length` L (1, the period) and `pressure` p₀ (0) |
| `analytic.type = "uniform"` | | uniform flow with `velocity` ([1, 0, 0]), `pressure` (0) and `length` (1, for the Reynolds number): exact between slip walls |
| `analytic.type = "taylor-couette"` | | laminar Taylor–Couette flow between coaxial cylinders about the z axis with `inner_radius` r₁ (0.05), `outer_radius` r₂ (0.10), `inner_omega` ω₁ (1) and `outer_omega` ω₂ (0): u_θ = A r + B/r; the pressure is known up to a constant, so the error compares the volume means |
| `analytic.type = "pipe-poiseuille"` | | Hagen–Poiseuille flow along the x axis with `radius` R (1), `mean_velocity` (1), `x_outlet` (1) and `pressure_outlet` (0) |

Built-in cylinder mesh `[mesh.cylinder]` (D34; boundaries `inlet`, `outlet`,
`walls`, `cylinder`, `zmin`, `zmax`): an O-grid around the cylinder that
fills the square 0 ≤ x, y ≤ height and a graded wake block to the outlet.

| Key | Default | Meaning |
|---|---|---|
| `length`, `height`, `depth` | 2.2, 0.41, 0.01 | channel length, height and the layer thickness in z (m) |
| `center`, `diameter` | [0.2, 0.2], 0.1 | cylinder centre (x, y) and diameter (m) |
| `n` | 32 | intervals per side of the O-grid (4n around the cylinder) |
| `radial`, `radial_ratio` | 32, 6 | intervals from the cylinder to the square, ratio of the last to the first |
| `wake`, `wake_ratio` | 64, 8 | intervals in x of the wake block, ratio of the last to the first |
| `refine` | 0 | halve every interval this many times (refine 0: 6 144 cells, each level four times more) |

Built-in annulus mesh `[mesh.annulus]` (boundaries `inner`, `outer`, `zmin`,
`zmax`): the gap between two coaxial cylinders about the z axis, a ring of
hexahedra; one layer in z is a 2D case.

| Key | Default | Meaning |
|---|---|---|
| `r_inner`, `r_outer`, `depth` | 0.05, 0.10, 0.01 | the radii and the extent in z (m) |
| `n_theta`, `n_r`, `n_z` | 64, 16, 1 | intervals around the axis, from cylinder to cylinder and layers in z (`refine` halves `n_z` only above 1) |
| `radial_ratio`, `symmetric` | 1, false | grading of the radial intervals: ratio of the last to the first (of the middle to the ends when symmetric) |
| `rotor_layers` | 0 | the innermost radial intervals that form the cell zone `rotor` (the rest is `stator`); 0: one zone `fluid` |
| `refine` | 0 | halve every interval this many times (the grading is nested) |

Rotating zones, the steady moving reference frame (D55, [recipe 94](../recipes/flow.md#mrf)):
one `[[mrf]]` table per zone.

| Key | Default | Meaning |
|---|---|---|
| `zone` | all cells | the cell zone of the mesh that turns; without it all cells turn, which needs a single table |
| `omega` or `rpm` | none | the angular velocity in rad/s or in revolutions per minute (exactly one), positive about the axis by the right-hand rule |
| `axis`, `origin` | [0, 0, 1], [0, 0, 0] | direction of the axis (any length) and a point on it (m) |

The case is steady, or transient with a `[time]` table (the transient
frozen rotor, D185: the rotor stays where it is; `MRF001` for a zone that a
`[[sliding]]` table turns as well); the solved velocity is the
absolute one, its convecting flux, and that of the turbulence variables and
passive scalars, is relative to the frame. A wall that turns with the frame
takes its absolute velocity in `velocity` (an expression, D53).

The transient form, the sliding mesh with conformal stepping (D56, D182,
[recipe 95](../recipes/flow.md#sliding-mesh)): one `[[sliding]]` table, the
rotor of a ring of identical sectors that advances by whole sectors in every
time step.

| Key | Default | Meaning |
|---|---|---|
| `zone` | all cells | the cell zone of the mesh that turns |
| `omega` or `rpm` | none | the angular velocity (exactly one, nonzero) |
| `axis`, `origin` | [0, 0, 1], [0, 0, 0] | direction of the axis and a point on it (m) |
| `shift` | 1 | sectors the rotor turns per time step; the time step must be `shift` · 2π / (N · \|Ω\|) for the N sectors found in the mesh (`SLD003` gives the value) |

The case is transient (`SLD001` for a steady one), the time step is fixed
(`SLD003` for the adaptive step, the error tolerance or several verification
levels), and the mesh must be a ring of identical sectors about the axis
(`SLD002`/`SLD004`: the zone, the sector count, the rotor boundaries).
Turbulence models run on both sliding variants (the rotor carries k, ε, ω,
ζ and f with their old time levels; the wall distance of the turning rotor
is updated in every step). The solved velocity is the absolute one; the convecting
flux is relative to the moving mesh, and there is no Coriolis source.

The general form, with a pitch ratio that is not one and a free time step
(D182 variant B, [recipe 96](../recipes/flow.md#sliding-mesh-imprint)): the
same table with the two boundaries of the interface, and a mesh in which the
rotor and the stator are separate meshes (`[mesh.annulus] n_theta_stator`).

| Key | Default | Meaning |
|---|---|---|
| `rotor_interface`, `stator_interface` | none | the boundaries of the mesh that hold the faces of the rotor and of the stator on the interface (both or neither; they get no `[boundary]` table); the interface is imprinted at the rotor angle of every time step |

With them `zone` is needed, `shift` is not allowed (the time step is free),
and the scheme must be fully implicit (`[solver] extrapolation = false`,
`outer_iterations` ≥ 2, `SLD003`); a restart is not available (`SLD005`).

Several rotors (D185, [recipe 102](../recipes/flow.md#several-rotors)): one
`[[sliding]]` table per rotor, each with its own `zone`, speed, `axis`,
`origin` and pair of interface boundaries to the stator; every table needs
the interface by imprinting (`SLD002` for conformal stepping with several
rotors, for a zone named twice), and a rotor may not border another rotor
(`SLD005`: each needs its interface with the stator). The rotors turn at
their own speeds and are imprinted together in every step; the log columns
`rotor.<zone>.angle`, `.torque` and `.power` name each rotor.
`[mesh.annulus] outer_rotor_layers` (with `n_theta_outer`) makes a second
rotor ring at the outer cylinder, the zone `rotor_outer` with the interface
boundaries `stator_interface_outer` and `rotor_interface_outer`.

# Derived fields (D72)

`[output] derived` adds fields computed from the solution to every
`fields*.vtkhdf`; `foehn post derive CASE` writes the same file afterwards
from a restart file ([recipe 30](../recipes/post-processing.md#derived-fields)).
G is the cell velocity gradient (Gᵢⱼ = ∂uᵢ/∂xⱼ, the weighted
least-squares gradient of the solver), S = ½ (G + Gᵀ), Ω = ½ (G − Gᵀ),
‖·‖ the Frobenius norm.

| Name | Datasets | Definition | Unit |
|---|---|---|---|
| `vorticity` | `vorticity_x`, `_y`, `_z`, `vorticity_magnitude` | ω = ∇ × u | 1/s |
| `q-criterion` | `q_criterion` | ½ (‖Ω‖² − ‖S‖²), vortex cores where Q > 0 | 1/s² |
| `lambda2` | `lambda2` | middle eigenvalue of S² + Ω², vortex cores where λ₂ < 0 | 1/s² |
| `enstrophy` | `enstrophy` | ½ ω·ω | 1/s² |
| `strain-rate` | `strain_rate` | √(2 S:S) | 1/s |
| `pressure-coefficient` | `pressure_coefficient` | (p − p_ref)/(½ ρ_ref U_ref²); needs a reference velocity | - |
| `total-pressure` | `total_pressure` | p + ½ ρ u·u | Pa |
| `relative-velocity` | `relative_velocity_x`, `_y`, `_z`, `relative_velocity_magnitude` | **w** = **c** − **u**, the velocity relative to the turning zone of the cell; **c** (the solved velocity) in the other cells (D183) | m/s |
| `frame-velocity` | `frame_velocity_x`, `_y`, `_z`, `frame_velocity_magnitude` | **u** = **Ω** × (**x** − **x**₀) of the turning zone of the cell (the blade speed), 0 in the other cells | m/s |
| `turbomachinery` | `velocity_radial`, `velocity_circumferential`, `velocity_axial`, `velocity_meridional`, `relative_velocity_circumferential`, `flow_angle`, `relative_flow_angle` | c_r, c_u, c_a, c_m = (c_r² + c_a²)^½ and w_u about the axis of the zone, which points along **Ω** (so the blade speed and a swirl with the rotation are positive); α = atan2(c_m, c_u) between **c** and the circumferential direction, β = atan2(c_m, −w_u) between **w** and the negative circumferential direction (90° without relative swirl) | m/s, ° |
| `relative-total-pressure` | `relative_total_pressure` | p + ½ ρ w·w | Pa |
| `rothalpy` | `rothalpy` | I = p/ρ + ½ w·w − ½ u·u, constant along a relative streamline of a steady flow without losses | m²/s² |
| `wall-shear-stress` | `wall_shear_stress_x`, `_y`, `_z`, `wall_shear_stress_magnitude`, on the boundary blocks | tangential viscous stress of the fluid on the wall (as `wall_<boundary>.csv`), zero off walls | Pa |
| `wall-distance` | `wall_distance` | exact distance of the cell centroid to the nearest face of the `wall` boundaries (BVH closest-point search, D48, D108); needs a wall | m |
| `y-plus` | `y_plus`, on the boundary blocks | d √(τ_w/ρ)/ν, d the wall distance of the cell centroid, zero off walls | - |
| `turbulence-features` | `feature_lambda1` … `feature_lambda5`, `feature_wall_reynolds`, `feature_time_ratio`, `feature_viscosity_ratio`, `feature_shielding`, `feature_q`, `feature_pressure_gradient`, `feature_curvature` | the invariant inputs of data-driven corrections (D106, D166), with a turbulence model only: the five invariants of Pope (1975) of the normalised Ŝ = S/ω and Ω̂ = Ω/ω (tr Ŝ², tr Ω̂², tr Ŝ³, tr Ω̂²Ŝ, tr Ω̂²Ŝ²), the wall Reynolds number min(√k d/(50 ν), 2), the time-scale ratio ‖S‖/ω, ν_t/ν, the shielding function F₁, (‖Ω‖² − ‖S‖²)/(‖Ω‖² + ‖S‖²), the pressure gradient along the streamline u·∇p/(‖∇p‖‖u‖ + ε) and the streamline curvature of the strain, each bounded to [−1, 1] or [0, …] as in D166 | - |
| `turbulence-scales` | `dissipation_rate`, `specific_dissipation_rate`, `time_scale`, `length_scale`, `viscosity_ratio` | ε (β* k ω for SST), ω (ε/(C_μ k) for the ε models), the time scale T of the model (1/(β* ω) for SST, k/ε for the k-ε family, the bounded T of ζ-f and EB), k^{3/2}/ε and ν_t/ν; with a turbulence model only (D72) | m²/s³, 1/s, s, m, - |
| `reynolds-stress` | `reynolds_stress_xx`, `_yy`, `_zz`, `_xy`, `_yz`, `_xz` | the Boussinesq estimate ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ of u'ᵢu'ⱼ; with a turbulence model only | m²/s² |
| `wall-normal` | `wall_normal_x`, `_y`, `_z` | n = ∇α/\|∇α\| of the elliptic-blending k-ε (zero where ∇α vanishes), the wall direction of the GGDH of stage 7; EB only | - |
| `turbulence-wall` | `friction_velocity`, `y_star`, `blending`, on the boundary blocks | u_τ of the wall treatment, y* = d C_μ^¼ k_P^½/ν of the wall cell and the viscous weight of the wall treatment (1 where the wall is resolved, `low-re` and the sublayer of `standard`; 1 − tanh((y⁺/5)²) of `automatic`; e^{−Γ} of `compound`; Γ = dU⁺/dy⁺ of `compound-2023`; 0 in the log layer and for `scalable`), zero off walls; with a turbulence model only. h_w and T⁺ of the scalars are the columns `<name>_h` and `<name>_t_plus` of `wall_<boundary>.csv` | m/s, -, - |
| `curvature-correction` | `rotation_function` | f_r1 of the curvature correction of SST (D170): 1 in pure shear, above 1 where curvature or system rotation amplify the turbulence, below where they damp it, 1.25 where the vorticity vanishes; with `turbulence.curvature_correction` only | - |

The datasets are written in the order of this table, whatever the order of
the list. `output.derived` and `output.reference` may change during a run.
The fields of rotating machinery need an `[[mrf]]` or `[[sliding]]` table; a
cell outside every turning zone takes its cylindrical components about the
axis of `output.frame`, else of the first table. On the sliding mesh **x** is
the current position of the cell. A sample line (`[[output.lines]]`) of a case
with a turning zone has the columns `u_rel`, `v_rel`, `w_rel`, `c_r`, `c_u`,
`c_a`, `w_u` at the end.

# Statistics of transient runs (D71)

`[output.statistics]` averages a transient run in simulated time, each step
weighted with its time step ([recipe 110](../recipes/post-processing.md#transient-statistics)).

| Output | Contents |
|---|---|
| `statistics.vtkhdf` | `mean_<f>` and `rms_<f>` of `velocity_x`, `velocity_y`, `velocity_z`, `pressure`, `k` (with a turbulence model) and every scalar; `covariance_velocity_xy`, `_yz`, `_xz`; `covariance_velocity_<c>_<scalar>` for c = x, y, z. Needs `[output] vtkhdf` and `fields = true` |
| `statistics.csv` | per series (`probes.<name>.u`, `.v`, `.w`, `.p`, `coefficients.cd`, `.cl`): `mean`, `rms`, `min`, `max`, `half_width_95` (Student t of the batch means), `batch_correlation` r1, `verdict` |
| `statistics_phase.csv` | with `phase_period`: `phase` (bin centre), `weight` (time in the bin, s) and the mean of every series in the bin |
| `result.json` | `solves[].statistics`: window, samples, batches, the series and the field datasets |

The verdict is `sufficient` when r1 passes the one-sided 95 % von Neumann
test, r1 ≤ 1.645 √((B − 2)/(B² − 1)) for B batches, else `too short` with
the warning STAT001: the batch means are correlated, so the half width
understates the uncertainty (a trend, or too short a window). With
`phase_period` the batches hold whole periods. A series also passes when its half width, widened by
√((1 + r1)/(1 − r1)), is below `tolerance` (10⁻³) of √var + |mean|. The
section is refused for a steady run (CASE005), and `batches` < 2 or a
`phase_period` ≤ 0 give CASE004.

# Boundary conditions (D29)

| `type` | Keys | Velocity | Pressure | Mass flux |
|---|---|---|---|---|
| `wall` | `velocity` (default 0; only the tangential part is kept), `mass_flux` (kg/(m² s) into the domain through the wall, D124 W3: blowing > 0, suction < 0; a number or an expression) | prescribed (with the transpiration velocity m''/ρ along −n) | extrapolated | 0, or −m''·\|S\| |
| `slip-wall` | none | no shear stress, zero normal component (as `symmetry`, on any surface) | extrapolated | 0 |
| `velocity-inlet` | one of `velocity`, `normal_velocity` (into the domain), `profile = "analytic"` (the exact solution, in a transient case at each time level), `profile = "parabolic"` with `max_velocity` and `profile_axes` (developed laminar inflow, for example `["y"]`), `profile = "table"` with `table` and `table_coordinate` (a measured or computed profile, below) | prescribed | extrapolated | ρ u_b·S |
| `mass-flow-inlet` | `mass_flow` (kg/s, spread by area) | normal, from the mass flow | extrapolated | prescribed |
| `static-pressure-inlet` | `pressure` | normal to the face, from the flux | prescribed | Rhie–Chow |
| `total-pressure-inlet` | `total_pressure` | normal to the face, from the flux | p₀ − ½ρu_n² | Rhie–Chow |
| `static-pressure-outlet` | `pressure`, `backflow` | extrapolated; on backflow by `backflow` | prescribed; on backflow by `backflow` | Rhie–Chow |
| `total-pressure-outlet` | `total_pressure`, `backflow` | extrapolated; on backflow by `backflow` | p₀ − ½ρu_n² | Rhie–Chow |
| `convective-outlet` | `convection_velocity` U_c (default: the mean outflow speed) | ∂u/∂t + U_c ∂u/∂n = 0 (transient cases) | extrapolated | ρ u_b·S, balanced with the inflow |
| `symmetry` | none | no shear stress, zero normal component | extrapolated | 0 |
| `wedge` | none | as `symmetry`: the two side planes of an axisymmetric wedge mesh ([`mesh.generate.wedge`](meshes.md)) | extrapolated | 0 |
| `periodic` | `partner`, one of `mass_flow` (kg/s) and `pressure_gradient` (Pa/m) on one side of the pair, `translation` (m) | across the pair: the partner's cells are the neighbours, translated by the period | periodic, with the driving as a uniform body force along the translation | through the pair |

**Pressure boundaries with gravity** (`pressure_kind`, D153). With
`[gravity]` the solver works with the reduced pressure
p′ = p − ρ_ref g·(x − x_ref). The values of the four pressure boundaries
are p′ by default (`pressure_kind = "reduced"`): 0 is an opening to an
ambient at rest at ρ_ref, at every height. `pressure_kind = "static"`
reads the value as the static pressure p and converts it per face,
p′_b = p_b − ρ_ref g·(x_f − x_ref). A uniform static value over a tall
opening therefore pushes fluid in at the bottom and out at the top
(BUO006); give the hydrostatic profile as an expression, or the reduced
value. Without gravity the key changes nothing ([recipe
186](../recipes/energy.md#buoyancy)).

**Backflow at pressure outlets** (`backflow`). Where the flow enters
through an outlet face (a vortex crossing the outlet, a recirculation
reaching it):

| `backflow` | Velocity of the entering fluid | Pressure | Suits |
|---|---|---|---|
| `"total-pressure"` (default) | normal to the face | the prescribed value as the total pressure, p − ½ρu_n² | outlets that vortices or a recirculation cross; converges like an outflow |
| `"zero-gradient"` | extrapolated | the prescribed value as the static pressure | the treatment before stage 5; converges slowly where the flow enters |
| `"suppress"` | the face is closed for the iteration: a slip wall with the extrapolated pressure | | outlets where nothing may enter; a closed face opens when the prescribed pressure drives the flow out |

**Convective outlet.** The velocity leaves the domain by the convection
equation ∂u/∂t + U_c ∂u/∂n = 0, discretised with the backward difference
of the run between the face and its cell, so vortices pass the outlet with
little reflection and the outlet can sit close behind a body. It needs a
transient case. The outlet flow is balanced with the prescribed inflow in
every iteration; without a pressure boundary the reference cell fixes the
pressure level. U_c defaults to the mean outflow speed (the prescribed
inflow over the outlet area); give `convection_velocity` when nothing
prescribes an inflow.

**Periodic pairs** (stage 6 part 6p, D29, D109 (1); [gallery 67](../gallery/boundary-conditions.md#periodic),
[recipe 67](../recipes/boundary-conditions.md#periodic)). Two boundaries of
`type = "periodic"` that name each other as `partner` (one of them may
leave it out) form a translational pair: the faces of one are translated
copies of the faces of the other, matched by their centroids (within 1e-6
of the face size), with the translation from the area-weighted centroids
or from `translation` (from this boundary to its partner). Each side keeps
its own faces; the neighbour of a face is the cell behind its partner,
seen shifted by the period, in every face loop, gradient and limiter. The
pressure is periodic and its level fixed in the reference cell; the flow
is driven by a uniform body force β along the translation, either fixed
(`pressure_gradient = G`: the mean pressure falls by G per metre from the
boundary that gives it towards its partner, β = G) or adjusted after every
pressure correction so that the mass flow through the pair is `mass_flow`
(kg/s, from the boundary that gives it towards its partner). The log, the
report (`periodic`: `mass_flow`, `driving_force`) and the restart file
carry β. Without a driving and without a body force the check warns
(BC031): the flow decays. The mesh needs at least two cells across the
periodic direction, is renumbered before the pairing and runs on one
partition for now; generated boxes keep the perturbed nodes of
`xmin`–`xmax` (and of y, z) translated copies when the case pairs them.
Errors: BC030.

**Inlet profile from a table** (`profile = "table"`). `table` is a file
(relative to the case file) of rows of a coordinate and either the normal
speed into the domain (two columns) or the three velocity components (four
columns); `#` starts a comment, columns are separated by commas,
semicolons or blanks, a first line that is not numeric is a header, and
the coordinate increases. Between the rows the values are linear, beyond
the ends constant (the check warns, BC020). `table_coordinate` is `"x"`,
`"y"`, `"z"` or `"radius"`, the distance from the axis through
`table_origin` (default the origin) along `table_direction` (needed). The
table's sha256 is part of the case hash (D85), and the table must lie
inside the read root of the run (D82).

A header can name the columns after the coordinate (stage 6, 6g): the
velocity as `u_n` or as `u_x`, `u_y`, `u_z`, and the turbulence of the
inlet as `k` with `epsilon` or `omega`, and `zeta` (ζ of k-ε-ζ-f, φ of
the elliptic-blending k-ε), in any order. The turbulence columns set k,
ε or ω and ζ face by face (ε and ω converted with the model's C_μ), so
that the developed profile of a periodic pipe or channel run feeds a jet
or a duct (`python3 -m foehncraft.profiles LINE.csv TABLE.csv` makes the
table from a line of `[[output.lines]]`). A boundary that gives
turbulence values as well (`k`, `turbulence_intensity`, …, or `zeta`
with a `zeta` column) is refused (BC029), and so is a header that names
an unknown column (BC019).

**Wedge** (axisymmetric flow). `mesh.generate.wedge = θ` (degrees, up to
90) revolves a 2D mesh of the half plane y ≥ 0 about the x axis (y is the
radius) into a wedge one cell thick, symmetric about the x–y plane; nodes
on the axis merge, so the cells next to it become prisms or polyhedra with
an edge on the axis and there is no axis boundary. `zmin` and `zmax` are the
wedge planes and take `type = "wedge"`. The faces are planar, so a wedge
of θ is a slice of a pipe whose section is a polygon of 360°/θ sides: the
result differs from the axisymmetric one by about θ²/8 (1e-3 at 5°); the
check warns above 10° (BC024).

Inlets also take the turbulence values of the
[turbulence models](#turbulence-d129d139): `turbulence_intensity` with
`turbulence_length_scale` or `viscosity_ratio`, or `k` with `omega`
directly (or `k` with `epsilon` for the k-ε family). Pressure outlets take them, and `zeta`, for the backflow:
on the faces where the flow enters (ṁ < 0) they are fixed values, where
it leaves zero gradient, decided face by face in every outer iteration
([gallery 140](../gallery/turbulence.md#outlet-backflow)); without them
the entering fluid carries the values of the cell. In laminar flow the
check reports them as unused (BC009). Walls take `roughness`.

# Sensible combinations (D30)

The check runs before every solve and with `--check`. Errors stop the run
and name the boundary and what to change (a missing or unknown boundary with
"did you mean", a net mass flow with no pressure boundary to take it, a
flow path whose rate is prescribed twice, a wall with a normal velocity, a
negative mass flow); warnings and notes say what the run will do (no
pressure boundary: the reference cell holds the level; nothing drives a
flow; unused turbulence values; a symmetry boundary that is not planar; a
table that is missing or malformed; a convective outlet in a steady case
or without a convection velocity; wedge planes that do not form one
wedge).

| Combination | Suits | Pressure level |
|---|---|---|
| velocity-inlet + static-pressure-outlet (+ walls, symmetry) | internal flows with a known flow rate: channels, pipes, ducts, the Poiseuille case | outlet pressure |
| mass-flow-inlet + static-pressure-outlet | as above when the mass flow is known rather than a velocity profile | outlet pressure |
| total-pressure-inlet + static-pressure-outlet | pressure-driven flows: the flow rate is a result (fans, valves, plenum to plenum) | both |
| static-pressure-inlet + static-pressure-outlet | pressure-driven flows with a known static pressure drop | both |
| static-pressure-inlet + total-pressure-outlet | flows that exhaust into a large volume at a known total pressure | both |
| walls only, one wall moving (+ symmetry) | enclosed flows: the lid-driven cavity | reference cell |
| velocity-inlet + velocity-inlet + static-pressure-outlet | mixing of two streams: the T-junction | outlet pressure |
| velocity-inlet + convective-outlet (transient) | wakes and vortex streets with the outlet close behind the body | reference cell |
| velocity-inlet + static-pressure-outlet + slip walls | external flow in a finite box: the slip walls stand for the far field | outlet pressure |
| velocity-inlet + static-pressure-outlet + wedges (+ walls) | axisymmetric flow: pipes, nozzles, jets | outlet pressure |
| periodic pair (`mass_flow` or `pressure_gradient`) + walls (+ symmetry, wedges, further periodic pairs) | fully developed flow in one module: channels, pipes, ribbed passages, heat transfer with the heat-flux jump | reference cell |

# Expressions, materials and regions (D53, D102, D112)

Stage 5 part 6a ([gallery 23](../gallery/properties.md#expressions-materials),
[recipe 31](../recipes/properties.md#expressions),
[recipe 32](../recipes/properties.md#materials)).

**Expressions.** Wherever the table below says so, a value may be a string
of the expression language instead of a number. The variables are `x`,
`y`, `z` (m, the face or cell centroid) and `t` (s, simulated time); the
constants `pi`, `e` and those of `[constants]`; the operators
`+ - * / % ^` (also `**`), `< <= > >= == !=`, `and or not` (true is 1);
the functions `sin cos tan asin acos atan sinh cosh tanh exp log log10
sqrt abs floor ceil round sign erf erfc` (one argument; `log` is the
natural logarithm, `erf` the error function and `erfc` = 1 − erf, for the
exact solutions of transient conduction), `atan2 pow min max hypot step`
(two; `step(a, b)` is 1 where
a ≥ b), `clamp(a, lo, hi)`, `if(condition, a, b)` and every
`[tables.<name>]` as a function of one argument. Each string is checked
and compiled when the case is read; the solver evaluates the compiled
program, on the CPU and the GPU alike, with the in-house elementary
functions of D93 (6), so a result does not depend on the platform's
mathematics library. A number stays a number: a case without expressions
runs bit for bit as before.

| Key | Takes | Meaning |
|---|---|---|
| `constants.<name>` | number or expression of other constants | a named constant; the order does not matter, cycles are errors |
| `tables.<name>.x`, `.y` | numbers | a piecewise linear table, clamped at the ends, called as `name(arg)` |
| `tables.<name>.file` | path | the same from a two-column file (CSV or blanks, `#` comments, header allowed) |
| `boundary.<boundary>.velocity`, `normal_velocity`, `mass_flux` | x, y, z, t | inlet and wall values per face |
| `boundary.<boundary>.pressure`, `total_pressure` | x, y, z, t | pressure boundaries per face |
| `boundary.<boundary>.mass_flow` | t | one value for the boundary |
| `flow.body_force` | x, y, z, t | per cell; an expression is added to the momentum source, a uniform number keeps the kernel path |
| `fluid.viscosity`, `materials.<name>.viscosity` | x, y, z | per face (properties do not change during a run, D100) |

`t` is allowed only in a transient case; time-dependent values are
evaluated at the new time level of every step.

**Materials.** `fluid.material` names a material of the case's
`[materials]` or of the built-in database; a `[materials.<name>]` table
with `base = "<database name>"` overrides single properties of a database
entry. Each property is a number, an expression or a table:

| Key | Meaning |
|---|---|
| `materials.<name>.base` | the material it starts from |
| `materials.<name>.kind` | `gas`, `liquid` or `solid` |
| `materials.<name>.temperature_range` | `[T_min, T_max]` (K): outside it the properties are clamped, with warning MAT003 |
| `materials.<name>.molar_mass` | kg/mol (the ideal-gas density) |
| `materials.<name>.density`, `.viscosity` | a number, an expression, or a table with `form` = `constant` (`value`), `polynomial` (`coefficients` a₀, a₁, … in T), `sutherland` (`mu_ref`, `t_ref`, `s`), `power-law` (`value_ref`, `t_ref`, `exponent`), `table` (`temperature` and `values`, or `file`; `interpolation` = `linear` or `monotone-cubic`), `ideal-gas` (density p M / (R T)), each with `source` and optionally `range`; the complete list is in [keys.md](../keys.md) |

The density must be constant in the domain (a density that varies comes
with the energy equation, MAT002). The run writes the resolved property
program into `case.resolved.toml`, and its hash into the provenance of
every output file (`properties_hash`, `materials_database`).

**Regions.** `[regions.<name>]` assigns cell zones to a material:

| Key | Meaning |
|---|---|
| `regions.<name>.kind` | `fluid` or `solid` (with `[energy]`, conjugate heat transfer); `porous` is refused (REG001) |
| `regions.<name>.zones` | cell zones of the mesh (each zone in one region; unknown zones are REG002) |
| `regions.<name>.material` | its material (default: that of `[fluid]`) |
| `regions.<name>.axes` | the first two principal axes of an orthotropic conductivity in this region, `[[..], [..]]` in the mesh frame (default: those of the material; one material oriented differently in every part) |

**Interfaces** (conjugate cases, stage 7d2, [recipe 160](../recipes/energy.md#cht-contact)).
`[interfaces.<name>]` puts a resistance into the faces between two regions:

| Key | Meaning |
|---|---|
| `interfaces.<name>.regions` | the two regions: solid regions, a fluid region of `[regions]`, or `fluid` (every fluid cell); not two fluids (CHT004) |
| `interfaces.<name>.contact_resistance` | R″_c (m² K/W, ≥ 0) |
| `interfaces.<name>.thin_wall` | `{ thickness, conductivity }` or `{ thickness, material }`: the resistance t/k of a wall that is not meshed (no heat capacity) |

**Region heat sources.** `[[energy.sources]]` with `region` (a solid or
fluid region, or `fluid`) and one of `power` (W over the region, spread by
volume; a number or an expression of t) and `power_density` (W/m³, a
number or an expression of x, y, z, t), added to `[energy] source`
(ENE008).

Regions with different viscosities give a viscosity per face, the harmonic
mean across a zone boundary.

# Passive scalars (D51, D116)

Stage 6 part 6a ([gallery 63](../gallery/scalars.md#passive-scalars),
[recipe 63](../recipes/scalars.md#passive-scalars)). Each
`[scalars.<name>]` table adds the equation
ρ ∂φ/∂t + ∇·(ṁ φ) − φ ∇·ṁ − ∇·(ρΓ ∇φ) = ρ S with the current mass fluxes;
the bounded form (the term −φ ∇·ṁ) keeps the upwind operator an M-matrix
while continuity is not converged.

| Key | Default | Meaning |
|---|---|---|
| `scalars.<name>.diffusivity` | required | Γ (m²/s), a number or an expression of x, y, z (per face) |
| `scalars.<name>.initial` | 0 | start value, a number or an expression of x, y, z |
| `scalars.<name>.source` | 0 | source per unit mass S, expression of x, y, z and t; explicit |
| `scalars.<name>.bounds` | none | [lower, upper]: the bounded BDF2 (D28 (1)) and the extremes against them in the log |
| `scalars.<name>.exact` | none | exact solution: L2 and largest error, observed order of a sweep |
| `boundary.<boundary>.scalar.<name>.value` | per kind | fixed value at inlets and walls (x, y, z, t) |
| `boundary.<boundary>.scalar.<name>.flux` | per kind | diffusive flux Γ ∂φ/∂n into the fluid at walls, slip walls, symmetry planes and wedges (x, y, z, t) |
| `boundary.<boundary>.scalar.<name>.periodic` | `"periodic"` | on one boundary of a periodic pair: `"periodic"` (unchanged across the pair) or `"heat-flux"` (the jump Δ = Q/ṁ of D109 (1): Q = ρ times the flux and source integrals, ṁ the mass flow through the pair; refused with fixed values on walls or inlets) |
| `boundary.<boundary>.scalar.<name>.bulk` | 0 | `"heat-flux"`: the bulk (mixing-cup) value on this boundary, the level of a steady run |
| `solver.scalar.scheme` | `linear-upwind` | `upwind`, `linear-upwind`, `central` (bounded, with `central_blending` and `central_bounding`) |
| `solver.scalar.limiter`, `limiter_k` | `venkatakrishnan`, 5 | limiter of the scalar gradients: `none`, `venkatakrishnan`, `venkatakrishnan-wang`, `barth-jespersen` (faces), `mlp` (vertices, D237, [recipe 211](../recipes/scalars.md#mlp-limiter)); the same keys in `[solver.energy]` |
| `solver.scalar.limiter_function` | `r3` | `mlp`: `r3`, `r4`, `r5` (Nishikawa), `barth-jespersen`, `venkatakrishnan` |
| `solver.scalar.limiter_epsilon` | 0.01 | `mlp`: c of the threshold ε = c (φ_max − φ_min), in the units of the field; 0 is strictly bounded |
| `solver.scalar.limiter_epsilon_scaling` | `domain` | `mlp`: `domain` (one ε) or `cell` (ε √(V_c/V_domain), shrinks with the mesh) |
| `solver.scalar.limiter_relaxation` | 0.5 | `mlp`: β of ψ ← ψ + β (ψ_new − ψ) between outer iterations |
| `solver.scalar.bounded` | false | linear-upwind face value clipped to the interval of its two cells |
| `solver.scalar.relaxation` | 0.9 steady, 1 transient | under-relaxation |
| `solver.scalar.tolerance`, `max_iter` | 1e-3 steady, 1e-8 transient; 50 | relative tolerance and iterations of each BiCGStab solve |
| `solver.scalar.sweeps` | 1 | solves per outer iteration, each with new deferred corrections |
| `verification.scalar_order` | none | accepted observed order of the scalar errors (mesh sweep) or of the differences of successive solutions (time-step sweep) |
| `verification.scalar_max_error` | none | largest accepted L2 error on the finest mesh |
| `verification.scalar_range` | none | [lower, upper]: every scalar (and the temperature) within it on every mesh |

The defaults per boundary kind: zero flux on walls, slip walls, symmetry
planes and wedges; the value 0 at inlets; zero gradient at outlets.
Periodic boundaries carry the scalar across the pair, plainly or with
the heat-flux jump. Scalar names may not be
those of output fields, derived fields or the expression language
(`T`, `k`, `omega`, `Y`, ... are reserved for the energy equation, the
turbulence models and the species). Errors in these tables are `SCAL001`.

# Turbulence (D129–D139)

Stage 6 part 6b ([gallery 68](../gallery/turbulence.md#sst),
[recipe 68](../recipes/turbulence.md#sst)): the k-ω SST model of Menter,
Kuntz & Langtry (2003) with its production limiter, the wall treatments
`low-re` (integration to the wall, y⁺ ≲ 1) and `automatic` (any y⁺),
rough walls after Aupoix (2015). Part 6c
([gallery 71](../gallery/turbulence.md#k-epsilon),
[recipe 71](../recipes/turbulence.md#k-epsilon)): the standard and
realizable k-ε models with the wall functions `scalable`, `standard` and
`compound`, rough walls after Cebeci & Bradshaw (1977). Part 6d
([gallery 73](../gallery/turbulence.md#zeta-f),
[recipe 73](../recipes/turbulence.md#zeta-f)): the k-ε-ζ-f model of
Hanjalić, Popovac & Hadžiabdić (2004) with the treatments `compound` and
`low-re`, its elliptic f equation and the multigrid of the turbulence
equations. Part 6d2 ([gallery 74](../gallery/turbulence.md#zeta-f0),
[recipe 74](../recipes/turbulence.md#zeta-f0)): the k-ε-ζ-f₀ model of
Popovac (2023) with his compound treatment of 2023, robust from y⁺ 0.3
to 30. Part 6e ([gallery 88](../gallery/turbulence.md#elliptic-blending),
[recipe 88](../recipes/turbulence.md#elliptic-blending)): the
elliptic-blending k-ε of Billard & Laurence (2012) with the treatments
`compound` and `low-re`. Specification:
[stage 6](../requirements/stage6-turbulence.md).

| Key | Default | Meaning |
|---|---|---|
| `turbulence.model` | `"laminar"` | `"k-omega-sst"`, `"k-epsilon"` (Launder & Spalding 1974), `"realizable-k-epsilon"` (Shih et al. 1995), `"k-epsilon-zeta-f"` (Hanjalić, Popovac & Hadžiabdić 2004), `"k-epsilon-zeta-f0"` (Popovac 2023, with the deviations of specification 3.3.1), `"elliptic-blending-k-epsilon"` (Billard & Laurence 2012, the BL-v²/k model: k, ε, φ = v²/k and the blending parameter α) |
| `turbulence.wall_treatment` | SST `"automatic"`, k-ε `"scalable"` | SST: `"automatic"` (Menter, Kuntz & Langtry 2003: blended viscous and log-law u_τ, ω and production in the first cell) or `"low-re"` (k = 0 and ω_w = 60ν/(β₁ d₁²) on the wall face, d₁ the wall-cell height). k-ε family: `"scalable"` (log law at max(y*, 11.06), any y⁺, accurate from y⁺ 30), `"standard"` (log law, 30 ≤ y⁺ ≤ 300), `"compound"` (blended profile of Popovac & Hanjalić; with the high-Re k-ε claimed from y⁺ 30 like the wall functions). k-ε-ζ-f: `"compound"` (default: k and ζ by zero gradient, ε and P_k of the first cell and the wall value of f blended between the viscous sublayer and the log law) or `"low-re"` (k = ζ = 0, ε_w = 2νk₁/d₁², f_w = −2νζ₁/d₁²); also `"compound-2023"`. k-ε-ζ-f₀: `"compound-2023"` (default: the blend u⁺ = [y*^{−10/3} + u⁺_log^{−10/3}]^{−3/10} of Popovac 2023 with the non-equilibrium factor ψ from the convection and the pressure gradient, ε and P_k of the first cell from its slope; f₀,w = 0), `"compound"` or `"low-re"`. Elliptic-blending k-ε: `"compound"` (default, as for ζ-f, with ε_w = νk₁/d₁² in the blend and α_w = 0) or `"low-re"` (k = φ = α = 0, ε_w = νk₁/d₁²) |
| `turbulence.production_limiter` | SST `"menter"`, k-ε `"durbin"` | `"menter"` (SST: P̃ = min(ν_t S², 10 β* k ω); k-ε: min(P_k, 10ε)), `"kato-launder"` (ν_t S Ω), `"durbin"` (the realisability bound with `durbin_coefficient`), `"none"`; the realizable k-ε takes `"none"` only |
| `turbulence.k_sink` | `"wall"` | k-ε-ζ-f and k-ε-ζ-f₀: the sink of k; `"wall"`: max(k/T_K, 2νk/d²), ε at the wall (k ~ y², converges with the mesh) and k/T_K away from it (k stays positive in laminar regions); `"epsilon"`: −ε of the papers; `"kolmogorov"`: k/T_K, T_K = max(k/ε, C_T (ν/ε)^½), the ζ-f₀ default before #141 (specification 3.3.2, [gallery 141](../gallery/turbulence.md#k-sink)) |
| `turbulence.ambient.k`, `turbulence.ambient.omega` or `epsilon` | none | sustaining terms of Spalart & Rumsey (2007): destruction at these free-stream values is balanced, so k and ω (or ε) do not decay ahead of a body; SST takes `omega`, the k-ε family, ζ-f and the elliptic-blending k-ε `epsilon` |
| `turbulence.max_viscosity_ratio` | 1e5 | cap of ν_t/ν; capped cells are counted (TURB008) |
| `turbulence.smooth_limits` | false | smooth min and max in the model (specification 2.6, for the adjoint) |
| `turbulence.curvature_correction` | false | SST: the curvature and rotation correction of Spalart & Shur (1997) in the form of Smirnov & Menter (2009): the production of k and ω times f_r1 = max[0, min(f_rotation, 1.25)], from r* = S/Ω and r̃ with the material derivative of the strain rate (D170, specification 3.6); f_r1 = 1 in pure shear; sees system rotation from `[[mrf]]` zones and `flow.frame_rotation`; not with `[[sliding]]` yet; constants `c_r1`, `c_r2`, `c_r3` |
| `turbulence.qcr` | "none" | SST: the quadratic constitutive relation, an explicit non-linear stress in the momentum equations (D170, specification 3.6): `"qcr2000"` (Spalart 2000), −C_cr1 (O_ik τ_jk + O_jk τ_ik) with O = 2W/|∇u|; `"qcr2013"` (Mani, Babcock, Winkler & Spalart 2013) adds −C_cr2 μ_t √(2 S*:S*) δ_ij (in incompressible flow a pressure shift only); drives the secondary flow of non-circular ducts; constants `c_cr1`, `c_cr2` |
| `turbulence.wall_distance` | `"exact"` | `"exact"` (BVH, D48) or `"pde"` (Poisson distance, one solve at the start) |
| `turbulence.constants` | literature values | overrides by name: SST `beta_star`, `a1`, `kappa`, `gamma1`, `gamma2`, `beta1`, `beta2`, `sigma_k1`, `sigma_k2`, `sigma_omega1`, `sigma_omega2`, and `c_r1` (1), `c_r2` (2), `c_r3` (1) of the curvature correction, `c_cr1` (0.3), `c_cr2` (2.5) of the QCR; k-ε `c_mu`, `c_epsilon1`, `c_epsilon2` (realizable: C₂), `sigma_k`, `sigma_epsilon`, `a0` (realizable), `kappa`; ζ-f `c_mu`, `c_epsilon1`, `c_epsilon2`, `c_1`, `c_2`, `sigma_k`, `sigma_epsilon`, `sigma_zeta`, `c_t`, `c_l`, `c_eta`, `kappa`, `c_epsilon1_zeta` (C″ε1 = 0.012; ζ-f₀ 0.009), `c_epsilon1_power` (1: C_ε1 = C′ε1(1 + C″ε1/ζ); ½: C″ε1/√ζ, as printed for ζ-f₀), `a` (0.6, the realizability bound of T and L; ζ-f₀ 0.849); ζ-f₀ also `n_zeta` (N_ζ = 2); the elliptic-blending k-ε `c_mu` (0.22), `c_epsilon1` (1.44), `c_epsilon2` (1.83), `c_epsilon3` (2.3, the E term), `c_epsilon4` (0.4, the lower end of C*_ε2; 1.83 makes C*_ε2 constant), `c_1` (1.7), `c_2` (0.9), `sigma_k` (1), `sigma_epsilon` (1.5), `sigma_phi` (1), `c_t` (4), `c_l` (0.164), `c_eta` (75), `a` (0.6, the realizability bound of T in ν_t), `kappa`; a warning names each changed constant (TURB007) |
| `initial.k`, `initial.omega` or `initial.epsilon` | from the first inlet | start values, numbers or expressions of x, y, z; ω converts to ε = C_μ k ω and back |
| `initial.zeta` | 2/3 | start ζ = v²/k of k-ε-ζ-f and ζ-f₀, or φ = v²/k of the elliptic-blending k-ε, in (0, 2]; f starts at 0, α at 1 |
| `initial.turbulence_intensity`, `initial.viscosity_ratio` | 0.01, 1 | start values from the intensity of `initial.velocity` and ν_t/ν |
| `boundary.<inlet>.turbulence_intensity` | 0.05 | I, k = 1.5 (I U)² with U the inlet speed |
| `boundary.<inlet>.turbulence_length_scale` or `viscosity_ratio` | ν_t/ν = 10 | ω = √k/(β*^¼ ℓ) or ω = k/(ν · ratio) |
| `boundary.<inlet>.k`, `omega` or `epsilon` | none | k and ω or ε directly, numbers or expressions (ω and ε together: BC029) |
| `boundary.<inlet>.zeta` | 2/3 | ζ = v²/k at an inlet (k-ε-ζ-f, ζ-f₀; φ of the elliptic-blending k-ε; in (0, 2]; BC028 with another model) |
| `boundary.<wall>.roughness` | 0 | equivalent sand-grain roughness k_s (m); SST: the wall values of k and ω of Aupoix (2015); k-ε: the log-law shift ΔB of Cebeci & Bradshaw (1977); BC027 when k_s exceeds the first-cell height on more than a tenth of the wall |
| `boundary.<wall>.roughness_constant` | 0.5 | C_s of Cebeci & Bradshaw (k-ε family); 0.4 brings the fully rough friction factor of pipes to Colebrook's within 3 % |
| `solver.turbulence.scheme` | `"linear-upwind"` | convection of k and ω: `upwind`, `linear-upwind`, `central` (transient runs only, TURB006) |
| `solver.turbulence.limiter`, `limiter_k` | `venkatakrishnan`, 0 | limiter of the gradients of k and ω; the threshold is relative, ε = K (φ_max − φ_min) |
| `solver.turbulence.relaxation` | 0.8 steady, 1 transient; 0.98 with a pseudo-time step (`solver.pseudo_time.turbulence`), 0.995 with the `"reduction"` preconditioner | under-relaxation of k and ω, ε, ζ |
| `solver.turbulence.viscosity_relaxation` | 1 | relaxation of ν_t |
| `solver.turbulence.tolerance`, `max_iter` | 1e-2 steady, 1e-8 transient; 50 | relative tolerance and iterations of each BiCGStab solve |
| `solver.turbulence.preconditioner` | `"reduction"` for the k-ε family with `coupled`, else `"jacobi"` | preconditioner of the BiCGStab solves of k, ω or ε and ζ: `"jacobi"`; `"reduction"`, the reduction multigrid of D179 on each equation's matrix with the settings of `[solver.coupled.multigrid]`, kept while the matrix changes by at most `reuse` (D180); or `"amg"` (the aggregation multigrid of `[solver.pressure.amg]`, one hierarchy per equation, built from the symmetric part ½(A + Aᵀ) and rebuilt every `solver.amg_rebuild` iterations; not with a pseudo-time step) |
| `solver.turbulence.reuse` | 0.02 | `"reduction"`: a hierarchy is kept while Σ\|Δa\| / Σ\|a\| of its equation's matrix stays at most this value, and rebuilt when BiCGStab needs more than 8 iterations with it; 0 builds one for every solve |
| `solver.elliptic.preconditioner` | `"amg"` | the CG solve of the elliptic f equation of ζ-f (α of the elliptic-blending k-ε): `"amg"` (own hierarchy) or `"jacobi"` |
| `solver.elliptic.tolerance`, `max_iter` | 1e-2 steady, 1e-8 transient; 100 | relative tolerance and iterations of the f (α) solve |
| `solver.turbulence.destruction` | ω `"newton"`, ε `"picard"` | linearisation of the ω or ε destruction: `newton` (2βω in the diagonal) or `picard` |

The checks: TURB001 (the first-cell y⁺ does not suit the wall treatment
on more than a tenth of a wall: above 2 with `low-re`, or
above 300 with `automatic` (validated on the flat plate for y⁺ 0.24 to 85,
C_f within 2.6 %); below 30 with `scalable`, outside 30 … 300 with `standard`
and `compound` of the k-ε family; above 2 with `low-re` of ζ-f, ζ-f₀ and the elliptic-blending k-ε,
whose `compound` and `compound-2023` cover every y⁺), TURB002 (the floor of k or ω acted in the final tenth of the
solves), TURB003 and TURB009 (inlet turbulence out of the usual range, or
defaulted), TURB005 (no wall), TURB006 (a scheme the model does not allow),
TURB007 (changed constants), TURB008 (ν_t capped).
The log shows the first-cell y⁺ per wall, the largest ν_t/ν and the floor
counts; `run.jsonl` adds the rows `turbulence` and `turbulence_floors`.
