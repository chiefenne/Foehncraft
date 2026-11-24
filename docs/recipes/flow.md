# Incompressible flow

Part of the [Recipes](../recipes.md): Couette, Poiseuille, cavity and cylinder runs, boundary conditions, transient runs, time-step control.

| # | Recipe |
|---|---|
| 18 | [Verify the flow solver on Couette and Poiseuille flow](#flow-verification) |
| 19 | [Run the lid-driven cavity and compare with Ghia et al.](#flow-cavity) |
| 20 | [Set up the boundary conditions of a flow case](#flow-boundaries) |
| 21 | [Run a transient flow case and check its time order](#transient-flow) |
| 23 | [Run the cylinder benchmark of Schäfer & Turek](#cylinder) |
| 28 | [Let the error choose the time step, write output at exact times, start from potential flow](#time-control) |
| 78 | [Converge a steady case in fewer iterations with the coupled algorithm](#coupled) |
| 80 | [Choose and tune the momentum sweep of the coupled algorithm](#momentum-multigrid) |
| 81 | [Use the multigrid of the turbulence equations and its larger steps](#turbulence-multigrid) |
| 83 | [See and tune the pressure hierarchy of the coupled preconditioner](#pressure-hierarchy) |
| 85 | [Add heat transfer to a turbulent flow: wall function and turbulent Prandtl number](#thermal-wall-function) |
| 86 | [Compare the heated channel with the DNS Nusselt numbers](#heated-channel-dns) |
| 94 | [Run a rotating zone in a fixed frame (pump, fan, mixer, rotor)](#mrf) |
| 95 | [Turn a rotor through the mesh in time (sliding mesh, conformal stepping)](#sliding-mesh) |
| 96 | [Turn a rotor against a stator with a different pitch (sliding mesh by imprinting)](#sliding-mesh-imprint) |
| 98 | [Mesh and run a rotor with blades from Gmsh, and draw its relative flow](#stirred-vessel) |
| 99 | [Write the relative flow of a rotor: relative and frame velocity, flow angles, rothalpy](#rotor-fields) |
| 100 | [Run a turbulent rotor: the frame or the sliding mesh with a turbulence model](#turbulent-rotor) |
| 101 | [Read the torque and power of a rotor, and the head and efficiency of a machine](#torque-power) |
| 102 | [Run several rotors at their own speeds, and a transient frozen rotor](#several-rotors) |
| 103 | [Give a rotor an interface of any rotationally symmetric shape: discs, cones, pots](#revolved-interfaces) |
| 225 | [Mesh a rotor with curved blades of constant thickness in Gmsh](#curved-blades) |

---

<a id="flow-verification"></a>
## 18. Verify the flow solver on Couette and Poiseuille flow

The exact solutions of stage 4 ([explanation](../verification/stage4.md)):
Couette flow must come out at the iteration error on every mesh type,
Poiseuille flow with an observed order of 1.8–2.2. Each case checks its own
result and fails otherwise.

```
cargo build --release -p foehncraft-cli
for m in hex tet poly mixed; do
  target/release/foehn flow examples/cases/flow/couette_${m}.toml
done
for m in hex tet poly mixed; do
  target/release/foehn flow examples/cases/flow/poiseuille_${m}.toml
done
cat output/poiseuille_poly/flow.md          # errors per level, observed order, boundaries
# plots: output/poiseuille_poly/convergence.svg, residuals.svg
```

Keys used: `[analytic]` (`type = "couette"` or `"poiseuille"`),
`profile = "analytic"` at the velocity inlet, `[verification]` `levels`,
`expected_order`, `max_error`. More levels or other settings from the
command line:

```
target/release/foehn flow examples/cases/flow/poiseuille_hex.toml --levels 0,1,2,3
target/release/foehn flow examples/cases/flow/poiseuille_hex.toml --threads 1 --max-iterations 8000
```

The checks as tests (the solving ones in release builds):

```
cargo test --release -p foehncraft-flow --test stage4
```

[↑ back to top](#incompressible-flow)

---

<a id="flow-cavity"></a>
## 19. Run the lid-driven cavity and compare with Ghia et al.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cavity_re100.toml       # about 5 s
target/release/foehn flow examples/cases/flow/cavity_re1000.toml      # about 1 min
target/release/foehn flow examples/cases/flow/cavity_re100_poly3d.toml
cat output/cavity_re1000/flow.md
# centre lines against Ghia et al.: output/cavity_re1000/centreline_u.svg, centreline_v.svg
# velocity and pressure along the lines: line_vertical.csv, line_horizontal.csv
# wall pressure and shear stress: wall_ymax.csv (lid) and the other walls
pip install vtk
python3 -m tools.mesh.render_mesh output/cavity_re1000/fields.vtkhdf cavity.png --field speed --top --no-edges
```

Keys used: `[verification]` `reference = "ghia-re100"` or `"ghia-re1000"`,
`reference_tolerance` (lid speeds), `lid_velocity`; `[[output.lines]]`
for the sampling lines. The Reynolds number is ρ U L / μ with the unit
cavity; the case is refused if it does not match the reference.

Convergence: the outer iterations grow with the mesh and the Reynolds
number. A larger `relaxation_velocity` (0.9 to 0.95) needs fewer; lower it
(0.7, or `scheme = "upwind"` for the first iterations) if a run diverges:

```toml
[solver]
relaxation_velocity = 0.9
tolerance = 1e-8
max_iterations = 10000
```

[↑ back to top](#incompressible-flow)

---

<a id="flow-boundaries"></a>
## 20. Set up the boundary conditions of a flow case

Every boundary of the mesh needs a `[boundary.<boundary>]` table. Check a case
before running it: the check (D30) names each boundary that is missing,
unknown (with "did you mean"), over-specified or unsuitable, and says what
to change; `--check` stops after it.

```
target/release/foehn mesh examples/cases/mesh/hex_2d.toml     # boundary names: output/hex_2d/mesh_report.md
target/release/foehn flow my_case.toml --check
```

A channel with a known flow rate:

```toml
[boundary.inlet]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]        # or normal_velocity = 1.0, or mass-flow-inlet
[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.walls]
type = "wall"
```

A pressure-driven flow (the flow rate is a result):

```toml
[boundary.inlet]
type = "total-pressure-inlet"
total_pressure = 10.0
[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0
```

An enclosed flow (no pressure boundary: the reference cell holds the level):

```toml
[boundary.lid]
type = "wall"
velocity = [1.0, 0.0, 0.0]
[pressure]
reference_point = [0.5, 0.5, 0.0]
reference_value = 0.0
```

The table of sensible combinations and the cases each one suits is in the
[case catalogue](../cases/flow-keys.md#sensible-combinations-d30); every condition and
key is listed in [boundary conditions](../cases/flow-keys.md#boundary-conditions-d29).

[↑ back to top](#incompressible-flow)

---

<a id="transient-flow"></a>
## 21. Run a transient flow case and check its time order

A case is transient when it has a `[time]` table; it then runs PISO with
BDF2 (or Crank–Nicolson) instead of the steady SIMPLEC iterations.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml    # about 40 s
target/release/foehn flow examples/cases/flow/taylor_green_cn.toml
target/release/foehn flow examples/cases/flow/taylor_green_space.toml   # space order, about 30 s
target/release/foehn flow examples/cases/flow/taylor_green_cfl10.toml   # Courant 1 to 9, outer iterations
cat output/taylor_green_bdf2/flow.md
# plot: output/taylor_green_bdf2/time_convergence.svg, per step: history.csv
```

The smallest transient case file adds this to a steady one:

```toml
[time]
scheme = "bdf2"          # or "crank-nicolson", "bdf1"
time_step = 0.01         # s
end_time = 2.0           # or steps = 200
adaptive = true          # optional: Δt from the largest Courant number
max_cfl = 0.8
```

Plain PISO (one outer iteration, two correctors) is second order in time
up to Courant numbers of about 1. For larger steps add outer iterations,
which end when the residuals of a step are below `outer_tolerance`:

```toml
[solver]
outer_iterations = 50
outer_tolerance = 1e-8
relaxation_velocity = 0.9   # below 1 only with outer_iterations > 1
```

A time-order sweep runs the case once per step of
`verification.time_steps` (largest first) on the same mesh and reports
the order from the differences of successive solutions; with
`expected_time_order` it fails outside that range:

```toml
[verification]
time_steps = [0.01, 0.005, 0.0025, 0.00125]
expected_time_order = [1.8, 2.3]
```

Keys: [transient cases](../cases/flow-keys.md#flow-case-file-keys); the method and the
results: [stage 5](../verification/stage5.md).

[↑ back to top](#incompressible-flow)

---

<a id="cylinder"></a>
## 23. Run the cylinder benchmark of Schäfer & Turek

```
target/release/foehn flow examples/cases/flow/cylinder_re20.toml     # steady, Re 20, about 1 min
target/release/foehn flow examples/cases/flow/cylinder_re100.toml    # periodic, Re 100, about 8 min
cat output/cylinder_re100/flow.md
# C_D, C_L per step: output/cylinder_re100/forces.csv, chart: coefficients.svg
# probes in front of and behind the cylinder: probes.csv
pip install vtk
python3 -m tools.mesh.render_mesh output/cylinder_re20/fields.vtkhdf cylinder.png --field speed --top --no-edges --size 1800x420 --zoom 2.9
```

The mesh is built in (`[mesh.cylinder]`); `refine = 1` and `2` give 24 576
and 98 304 cells. Re 100 with `refine = 1` takes about 40 minutes on two
cores and brings C_L,max to 0.5 % of its range; with a finer mesh lower the
first `time_step` (the adaptive step starts from it). On 98 304 cells the steady case needs
`relaxation_velocity = 0.5` and about 12 000 outer iterations:

```toml
[mesh.cylinder]
refine = 2
[solver]
relaxation_velocity = 0.5
max_iterations = 12000
```

Keys: `[verification]` `benchmark = "schafer-turek-2d1"` or
`"schafer-turek-2d2"`, `benchmark_tolerance` (relative to the midpoints of
the reference ranges), `evaluation_start` (2D-2). Force coefficients of
any case: `[output.coefficients]` with `boundaries`, `velocity`, `area`;
point histories: `[[output.probes]]` with `name` and `point`.

[↑ back to top](#incompressible-flow)

---

<a id="time-control"></a>
## 28. Let the error choose the time step, write output at exact times, start from potential flow

**The time step from the error** (D65). With `adaptive = true` and
`error_tolerance`, BDF2 and Crank–Nicolson estimate the local error of
each step from the difference between the step and the quadratic
extrapolation of the last three levels, and take the next step that meets
the tolerance; `max_cfl` and `max_growth` stay upper bounds. A step whose
error is more than the tolerance is taken back exactly and repeated
shorter; the log names each one, `history.csv` holds the steps taken, and
`result.json` counts the repeats under `solves[].transient.rejected`.

```
target/release/foehn flow examples/cases/flow/taylor_green_error.toml --out output/tg_error
grep -E "rejected|error" output/tg_error/flow.log | head
target/release/foehn flow examples/cases/flow/taylor_green_error.toml --set time.error_tolerance=1e-6 --out output/tg_error6
python3 -m tools.plot.time_step_chart output/time_step_error.svg
```

```toml
[time]
scheme = "bdf2"            # or "crank-nicolson"
time_step = 0.001          # the first step; the estimate starts with the fourth
adaptive = true
max_cfl = 5.0              # upper bound
max_growth = 1.2           # upper bound of the step ratio
error_tolerance = 1e-4     # relative to the RMS velocity, per step
error_absolute = 1e-6      # m/s
```

The tolerance is per step: the error at the end of a run is up to the
number of steps times larger. Keep it at least ten times above the
tolerances of the linear solvers (`solver.momentum.tolerance`,
`pressure_tolerance`), whose noise the estimate would otherwise chase. With
plain PISO at Courant numbers above about 2 the error grows faster than
the step, and the estimate cuts the step there (the steep drops in the
gallery chart); outer iterations (`solver.outer_iterations`) let it take
longer steps. A restart continues such a run bitwise.

**Output at exact times** (D66). In simulated time instead of steps:

```toml
[output]
fields_interval = 0.1      # fields_<step>.vtkhdf at t = 0.1, 0.2, ... s
fields_times = [0.25]      # and at these times
monitor_interval = 0.05    # forces.csv and probes.csv rows every 0.05 s
[run]
restart_interval = 1.0     # restart.h5 at t = 1, 2, ... s
restart_times = [2.5]      # restart_<step>.h5 at these times
```

The step before a scheduled time is shortened to land on it, and when the
time lies within two steps, the two steps share the distance, so no step
falls below half the planned one. The files keep the step number in their
names; the time is inside (the root attribute `time` of the VTKHDF file,
the `time` column of the CSV rows).

**Potential-flow start** (D52). `initial.potential_flow = true` replaces
the uniform start velocity by the irrotational flow that carries the
inflow to the pressure boundaries: one Laplace solve with the pressure
multigrid and four passes of the non-orthogonal correction; walls and
symmetry planes let it slip. To see the start field itself, run no outer
iteration:

```
target/release/foehn flow examples/cases/flow/cylinder_potential.toml --out output/potential
python3 -m tools.mesh.render_mesh output/potential/fields.vtkhdf potential.png --field speed --top --no-edges --size 1800x420 --zoom 2.9
target/release/foehn flow examples/cases/flow/cylinder_re20.toml --set initial.potential_flow=true --out output/cyl_pf
```

The last run converges in 1527 outer iterations instead of 1570 from rest.
A restart file, when given, wins over the start.

[↑ back to top](#incompressible-flow)

---

<a id="coupled"></a>
## 78. Converge a steady case in fewer iterations with the coupled algorithm

One key switches a steady case from SIMPLEC to the coupled
pressure–velocity algorithm with pseudo-time steps
([gallery 78](../gallery/flow.md#coupled)):

```toml
[solver]
algorithm = "coupled"
```

or on the command line, for any steady case file:

```
target/release/foehn flow examples/cases/flow/channel_zeta_f0.toml --set solver.algorithm=coupled
```

Remove a `[solver.momentum] relaxation` and a relaxation `[[schedule]]`
that the case set for SIMPLEC: with the coupled algorithm the relaxation
is only a cap on the local Courant number α/(1 − α), and its default is
0.98. The result is the one of SIMPLEC; only the number of iterations
changes.

The step case of gallery 78 and the cavity:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml
target/release/foehn flow examples/cases/flow/cavity_re1000_coupled.toml
```

**Reading the log.** The table has one more column of linear
iterations: `K`, the Krylov iterations of the coupled solve of that outer
iteration (`linear.coupled` in `[log] columns`; `krylov` in
`residuals.csv` and the run record). `U V W` are then the multigrid
cycles of all momentum sweeps of the preconditioner (the BiCGStab
iterations with `momentum_sweep = "bicgstab"`, recipe 80) and `P`
the CG iterations of its pressure solves and of the final pressure
correction. A `K` at its limit (`solver.coupled.max_iter`, 30) in every
iteration means that the step or the cap is too large for the
preconditioner: lower `solver.pseudo_time.factor` or the relaxation.

**The pseudo-time step.** Default: automatic, Δτ = L/U with L the volume
of the domain over the area of its boundaries (without symmetry planes)
and U the larger of the largest prescribed boundary speed (inlets, moving
walls) and the root mean square of the cell speeds. The summary of the
run shows the step of the last iteration (`Pseudo-time step`), the run
record has it per iteration (`pseudo_step`).

```toml
[solver.pseudo_time]
factor = 4.0        # a larger step: fewer iterations, dearer linear solves
# step = 5.0        # or the step itself, in s
# momentum = false  # relaxation only (the cap alone sets the step)
# turbulence = 0.5  # c_T of the turbulence step min(Δτ, c_T T); 0: off
```

A larger `factor` pays on diffusion-dominated laminar cases (cavity at
Re = 400: 292, 130 and 84 iterations with 1, 4 and 20). On the turbulent
step the count does not change between 5 and 29 H/U. These keys and those
of `[solver.coupled]` may change during a run (re-read, `[[schedule]]`).
Iteration counts repeat exactly for the same input, but they are good to
about 10 % only: a change in the last digit of the step can move them
(cavity at Re = 1000: 423 iterations with `step = 0.25`, 452 with the
automatic step, which is the same to 16 digits). The result is the same.

**The cap.** `solver.momentum.relaxation` (0.98) caps the local Courant
number of the momentum equations at 49, `solver.turbulence.relaxation`
(0.98 with a pseudo-time step) that of the turbulence equations. A higher
momentum cap (0.99, 0.999) saves outer iterations on fine meshes and
costs more in each; lower it (0.95) if a run diverges at the start.
Leave the turbulence value: 0.9 makes the turbulence lag the momentum
equations and costs several times the iterations, 1 can diverge from a
uniform start.

**With SIMPLEC.** The turbulence step alone also helps the segregated
algorithm; it is off by default there, so that existing results do not
change:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0.toml \
    --set solver.pseudo_time.turbulence=0.5
```

The coarse step then settles in about 230 iterations at α_u = 0.8
instead of 525.

**Has the result settled?** As with SIMPLEC the residual of the step
case stalls near 10⁻⁴, so judge by a monitored value: write fields at
chosen iterations (`touch output/<case>/.result_50 .result_100 …` before
the run, recipe 74) and compare, or watch a probe in the log
(`[log] columns`, recipe 38). The figure of gallery 78 was drawn from
such histories:

```
python3 -m tools.plot.settling_chart docs/figures/gallery/coupled_step.svg "x_r / H" \
    "k-ε-ζ-f₀:SIMPLEC=simplec_zf0.csv" "k-ε-ζ-f₀:coupled=coupled_zf0.csv" …
```

with one two-column file (iteration, value) per run.

[↑ back to top](#incompressible-flow)

---

<a id="momentum-multigrid"></a>
## 80. Choose and tune the momentum sweep of the coupled algorithm

The preconditioner of the coupled algorithm (recipe 78) applies one
V-cycle of a reduction multigrid to each velocity component
([gallery 80](../gallery/flow.md#momentum-multigrid)). That is the
default; a case needs no key for it:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml
target/release/foehn flow examples/cases/flow/cavity_re1000_coupled.toml
```

**Reading the log and the report.** The columns `U V W` count the
multigrid cycles of the momentum sweeps of an outer iteration: one per
Krylov iteration (`K`) and component with the default settings, none for
a component without a residual (the third one of a planar case). The
report of the run (`flow.md`) has the time of the cycles under
`momentum solve` and that of the setups under `multigrid setup`, with
the number of hierarchies built for the momentum sweeps:

```
multigrid setup 1.78 (45 setups, 28 of the momentum sweeps)
```

**Compare with the old sweep.** BiCGStab with a diagonal preconditioner,
to `momentum_tolerance` or `momentum_max_iter`:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml \
    --set solver.coupled.momentum_sweep=bicgstab
```

The converged result is the same; the path to it and the counts differ
a little (the fine step settles in 175 instead of 150 iterations).

**The keys.**

```toml
[solver.coupled]
momentum_sweep = "multigrid"   # or "bicgstab"
momentum_cycles = 1            # cycles per sweep
momentum_reuse = 0.02          # keep a hierarchy up to this change of the matrix

[solver.coupled.multigrid]
strength = 0.3                 # strong couplings: above this part of the largest of a row
restriction_strength = 0.05    # the same for the entries that enter the restriction
drop_tolerance = 1e-4          # coarse entries below this part of the diagonal are dropped
max_coarse = 64                # unknowns of the coarsest level (dense inverse)
max_levels = 30
precision = "mixed"            # hierarchy in float32; "double": float64
```

- `momentum_reuse`: a hierarchy is kept while the sum of the absolute
  changes of the matrix coefficients since its setup stays below this
  part of the sum of their absolute values. 0 builds one in every outer
  iteration, which costs more than it gains (a setup takes as long as
  twenty to thirty cycles). Above about 0.05 a kept hierarchy stops
  converging; the guard then rebuilds it at once (a kept hierarchy must
  halve the residual in the first sweep of an outer iteration), so a
  large value costs time but not the run.
- `momentum_cycles = 2`: for a run whose `K` stays high although the
  step and the cap are the defaults. On the cases of the gallery it does
  not pay.
- `[solver.coupled.multigrid]`: leave the defaults. A larger
  `drop_tolerance` (10⁻³) gives a smaller hierarchy and is enough with
  the cap 0.98; without the cap it loses the convergence on closed
  streamlines.

The three keys of the sweep may change during a run (re-read,
`[[schedule]]`); those of `[solver.coupled.multigrid]` are fixed for a
run.

**If a setup fails** (a zero on the diagonal, a coarsening that stalls),
the run goes on: BiCGStab does the momentum sweeps of that outer
iteration, and the log carries one warning,
`the multigrid setup of the momentum sweep failed (…); BiCGStab does the
sweeps of such iterations`. Send such a case in: it should not happen.

**Restart.** The restart file holds the matrices of the kept
hierarchies (`momentum_amg_matrix_0` to `_2`); the continued run
rebuilds them and stays bit-identical. A restart file from before this
change works too (the hierarchies are built anew).

**The figure of gallery 80.** The Krylov iterations per outer
iteration come from `residuals.csv`:

```
python3 -m tools.plot.krylov_chart docs/figures/gallery/momentum_multigrid.svg \
    "8 200 cells:BiCGStab sweep=output/step_bicgstab" "8 200 cells:multigrid sweep=output/step_multigrid" \
    "32 800 cells:BiCGStab sweep=output/fine_bicgstab" "32 800 cells:multigrid sweep=output/fine_multigrid"
```

[↑ back to top](#incompressible-flow)

---

<a id="turbulence-multigrid"></a>
## 81. Use the multigrid of the turbulence equations and its larger steps

With the coupled algorithm (recipe 78) in a steady run, the transport
equations of k, ε and ζ of the k-ε family (k-ε, realizable k-ε,
k-ε-ζ-f, k-ε-ζ-f₀) take one cycle of the reduction multigrid of
recipe 80 as the preconditioner of BiCGStab, and the relaxations of the
momentum and turbulence equations default to 0.995 instead of 0.98
([gallery 81](../gallery/flow.md#turbulence-multigrid)). A case needs no
key for it:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml
```

The step settles in about 40 outer iterations and converges to its
tolerance in about 115 (recipe 83).

**Reading the report.** `flow.md` counts the hierarchies built for the
turbulence equations after those of the momentum sweeps:

```
multigrid setup 1.17 (4 setups, 76 of the momentum sweeps, 114 of the turbulence equations)
```

While the recirculation zone forms, most solves rebuild their
hierarchy; once the flow has settled they are kept. The BiCGStab counts
of k, ε and ζ are in `residuals.csv` (columns `bicgstab_k`, `bicgstab_epsilon`, `bicgstab_zeta`): 1 to 2
per solve with the multigrid, up to the limit 50 with Jacobi at 0.995.

**Back to the defaults of part 6f** (Jacobi and 0.98), for a comparison
or a case that does not take the larger steps:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml \
    --set solver.turbulence.preconditioner=jacobi \
    --set solver.turbulence.relaxation=0.98 --set solver.momentum.relaxation=0.98
```

Setting `preconditioner = "jacobi"` alone keeps the relaxations at 0.98
too; setting a relaxation keeps the others at their defaults. A
relaxation above 0.99 with Jacobi runs the solves into their limit and
the result jitters: do not combine them.

**The keys.**

```toml
[solver.turbulence]
preconditioner = "reduction"   # default with the coupled algorithm and the k-ε family;
                               # "jacobi" elsewhere; "amg" is the aggregation multigrid
                               # of recipe 73 (not with the coupled algorithm: its
                               # coarsest level is not positive definite there)
reuse = 0.02                   # keep a hierarchy up to this change of the matrix
relaxation = 0.995             # default with "reduction", else 0.98

[solver.momentum]
relaxation = 0.995             # the cap on the local step; default with "reduction", else 0.98
```

- `reuse`: as `momentum_reuse` of recipe 80, the sum of the absolute
  changes of the matrix since the setup over the sum of its absolute
  values; a hierarchy whose solve needs more than 8 BiCGStab iterations
  is rebuilt in the next solve anyway. 0 builds one in every solve.
- The multigrid settings are those of `[solver.coupled.multigrid]`.

**k-ω SST** keeps Jacobi and 0.98: with accurate solves of k and ω its
steady state on the step is not held (gallery 81). `--set
solver.turbulence.preconditioner=reduction` on an SST case is accepted
but not recommended.

**If a setup fails**, the solves of that equation take Jacobi until a
setup succeeds, and the log carries one warning, `the multigrid setup
of a turbulence equation failed (…); its solves take the Jacobi
preconditioner until a setup succeeds`.

**Restart.** The restart file holds the matrices of the kept
hierarchies (`turb_amg_matrix_<slot>`); the continued run stays
bit-identical. A restart file from before this change works too (the
hierarchies are built anew).

**The figure of gallery 81.** Write x_r every 5 iterations (recipe 78)
into a two-column file (iteration, value) per run, then:

```
python3 -m tools.plot.settling_chart docs/figures/gallery/turbulence_multigrid.svg "x_r / H" \
    "8 200 cells:part 6f (relaxations 0.98)=coarse_6f.csv" "8 200 cells:part 6g (0.995, multigrid)=coarse_6g.csv" \
    "32 800 cells:part 6f (relaxations 0.98)=fine_6f.csv" "32 800 cells:part 6g (0.995, multigrid)=fine_6g.csv"
```

[↑ back to top](#incompressible-flow)

---

<a id="pressure-hierarchy"></a>
## 83. See and tune the pressure hierarchy of the coupled preconditioner

With the coupled algorithm (recipe 78) each Krylov iteration solves the
pressure matrix of the SIMPLE-type preconditioner with CG and the
aggregation multigrid of `[solver.pressure]`. Its hierarchy is built
anew when that matrix has changed by more than 5 % since the setup,
besides every `solver.amg_rebuild` iterations, and the solve stops at
the relative tolerance 0.3 ([gallery 83](../gallery/flow.md#pressure-hierarchy)).
A case needs no key:

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml --threads 1
```

The step settles in about 40 outer iterations and converges to its
tolerance in about 115, with every check passed.

**Reading the report.** The first count under `multigrid setup` in
`flow.md` is the number of pressure hierarchies built:

```
multigrid setup 1.29 (10 setups, 76 of the momentum sweeps, 113 of the turbulence equations)
```

With the periodic rebuild alone it is one per `amg_rebuild` iterations
(4 in 40). The CG iterations of the pressure solves of an outer
iteration are in `residuals.csv` (column `cg`), the Krylov iterations in
`krylov`; with the hierarchy following its matrix the fine step needs
about 8 CG iterations per Krylov iteration instead of 16.

**Compare with the periodic rebuild** of the runs before #83, or with
the pressure tolerance 0.1 of the runs before #82:

```
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml --threads 1 \
    --set solver.coupled.pressure_reuse=1 --out output/step_periodic
target/release/foehn flow examples/cases/flow/backstep_zeta_f0_coupled.toml --threads 1 \
    --set solver.coupled.pressure_reuse=1 --set solver.coupled.pressure_tolerance=0.1 \
    --out output/step_before_82
```

The converged results are the same. To its tolerance this coarse case
takes 12.4 s with the defaults (115 iterations), 12.5 s with the
periodic rebuild (106) and 14.6 s with both settings of before #82
(117), on one thread. The new rule pays where the pressure solves are a
larger part of the run: on finer meshes and to the settled result
(gallery 83).

**The keys.**

```toml
[solver.coupled]
pressure_tolerance = 0.3   # relative tolerance of the pressure solve of the preconditioner
pressure_max_iter = 20     # CG iterations at most
pressure_reuse = 0.05      # new hierarchy when the matrix has changed by more;
                           # 0: whenever it has changed; 1: only every amg_rebuild iterations

[solver]
amg_rebuild = 10           # the periodic rebuild, as before
```

- `pressure_reuse` uses the measure of `momentum_reuse` (recipe 80): the
  sum of the absolute changes of the matrix since the setup over the sum
  of its absolute values. Smaller values build more hierarchies; on the
  step 0.05 was the best of 0.02, 0.05 and 0.2, and laminar cases do not
  reach it.
- A tighter `pressure_tolerance` costs CG iterations without saving
  Krylov iterations; a much looser solve (at most two CG iterations)
  costs Krylov iterations (requirements step 22, section 7.9).

**Restart.** The restart file holds the matrix of the last pressure
setup, as before; the continued run decides the next rebuild from it
and stays bit-identical.

[↑ back to top](#incompressible-flow)

---

<a id="thermal-wall-function"></a>
## 85. Add heat transfer to a turbulent flow: wall function and turbulent Prandtl number

A temperature of constant properties is a passive scalar (D51). Give it
the molecular diffusivity ν/Pr, its boundary values and, with a
turbulence model, nothing else: the defaults are the wall function of
Kader and Pr_t = 0.85.

```toml
[scalars.theta]
diffusivity = 1.4084507042253522e-5   # ν/Pr, Pr = 0.71
turbulent_prandtl = 0.85              # constant model
# turbulent_prandtl_model = "kays-crawford"   # raises Pr_t near walls, for low Pr
wall_function = "kader"               # "jayatilleke" (log region) or "none"

[boundary.plate.scalar.theta]
value = 1.0                           # or flux = q (K m/s: q_w/(ρ c_p))
```

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/pipe_heated_keps.toml
target/release/foehn flow examples/cases/flow/plate_heated_keps.toml
```

**Reading the result.** The wall file `wall_<boundary>.csv` gains four
columns per scalar: `theta_wall` (the wall value, for a fixed flux
φ_P + q/h_w), `theta_flux` (the flux into the domain), `theta_h` (the
exchange velocity flux/(wall − cell value), h_w/(ρ c_p) of a
temperature) and `theta_t_plus` (u\*/theta_h). The Stanton number is
`theta_flux/(U (theta_wall − θ_∞))`; the Nusselt number of the pipe
needs the wall value at the periodic plane (the mean over the wall
minus G L/2, gallery 85).

**Choosing.** Leave `wall_function` at `kader` on any mesh: it reduces to
the molecular conduction where y⁺ is below 1. `none` is for comparison
(Nu 39 % low at y⁺ 32). Use `kays-crawford` for fluids of low Prandtl
number, where the molecular conduction reaches further into the layer.
Roughness has no thermal correction: the smooth-wall T⁺ is used.

---

<a id="heated-channel-dns"></a>
## 86. Compare the heated channel with the DNS Nusselt numbers

```
cargo build --release -p foehncraft-cli
python examples/studies/channel_heat_dns.py --run
```

The script downloads the tables of the
[JAXA DNS database](https://jaxa-dns-database.jaxa.jp/channelflow.html)
(Re_τ 180, 395, 640 and 1020 at Pr 0.71 and 0.025) into `output/dns_jaxa`,
which git ignores, prints U_b⁺, Re_b, Θ_b⁺ and Nu on 2h and 4h next to
the Nu of each file's header, and with `--run` runs
`examples/cases/flow/channel_heated_keps.toml` (Re_τ = 395, Pr 0.71) and
prints its Nu against the DNS (−1.4 % against the integrated table, −1.0 %
against the header). The data are JAXA's (all rights reserved, citation of
Abe et al. requested): do not commit them.

For the other Reynolds and Prandtl numbers see [recipe 87](#channel-dns-reynolds).

<a id="channel-dns-reynolds"></a>
## 87. Run the heated channel at another Re_τ or Pr

```
cargo build --release -p foehncraft-cli
python examples/studies/channel_heat_dns.py --run
```

runs all eight DNS points (Re_τ 180, 395, 640, 1020 at Pr 0.71 and
0.025) with both Pr_t models and prints Nu against the DNS. Each case is
`channel_heated_keps.toml` with

- `fluid.viscosity = 1/Re_τ` (u_τ = 1 m/s, h = 1 m),
- `scalars.theta.diffusivity = ν/Pr`,
- `mass_flow = ρ U_b × 0.2 m²` with U_b = Re_b ν/(2h) from the DNS bulk
  Reynolds number,
- `cells = [2, N, 1]` with N ≈ Re_τ/30, so that the first cell centre
  h/N sits at y⁺ ≈ 30.

For a liquid metal (Pr ≲ 0.1) set `turbulent_prandtl_model =
"kays-crawford"` (`channel_heated_pr0025.toml`): the constant Pr_t = 0.85
gives Nu 12–22 % high. With the k-ε family keep the first-cell y⁺ near 30
on heated walls: the scalable treatment takes the momentum and the
thermal law at y⁺ ≥ 11.06, but the model itself has no damping, and Nu
rises to +5 % at y⁺ 16 and +10 % at y⁺ 8. For wall-resolved meshes
(y⁺ ≲ 1) use ζ-f or EB with the compound treatment.

<a id="mrf"></a>
## 94. Run a rotating zone in a fixed frame (pump, fan, mixer, rotor)

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/taylor_couette_mrf.toml
```

A zone of the mesh that turns steadily is a `[[mrf]]` table (D55):

```toml
[[mrf]]
zone = "rotor"            # a cell zone of the mesh (Gmsh volume group, [mesh.annulus] rotor_layers)
axis = [0.0, 0.0, 1.0]    # direction of the axis, any length
origin = [0.0, 0.0, 0.0]  # a point on the axis (m)
rpm = 1500                # or omega = ... in rad/s
```

- The case is **steady**, or **transient** with a `[time]` table (the
  transient frozen rotor, D185, recipe 102): the rotor keeps its position.
  The rotor that moves is the sliding mesh (D182).
- The solved velocity is the **absolute** one. Walls that turn with the zone
  take their absolute velocity from `velocity`: `["-spin * y", "spin * x", 0.0]`
  for a rotation about the z axis through the origin (D53). Walls that stand
  keep the default.
- The interface of a zone, and every boundary of it that is not a wall, should
  be a surface of revolution about the axis (a cylinder, a disc, a cone). The
  method conserves mass in the relative flux only there; `MRF003` names the
  number of cells where it does not.
- The result does not depend on the frame for an axisymmetric flow, which is
  the check of `taylor_couette_mrf.toml`: change `omega` of the table (also to 0)
  and the error against `[analytic]` stays at the discretisation error.
- With the Kato–Launder production of SST the vorticity is the absolute one in
  the zone (note `TURB004`).
- `[mesh.annulus]` makes the gap between two cylinders for such tests:
  `r_inner`, `r_outer`, `n_theta`, `n_r`, `rotor_layers`, `refine`.

[Gallery 94](../gallery/flow.md#mrf); the keys are in
[cases/flow-keys.md](../cases/flow-keys.md) and [keys.md](../keys.md).

<a id="sliding-mesh"></a>
## 95. Turn a rotor through the mesh in time (sliding mesh, conformal stepping)

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/taylor_couette_sliding.toml
```

A zone of the mesh that turns, followed in time, is a `[[sliding]]` table
(D56, D182, variant A). The mesh is a ring of N identical sectors about the
axis (the zone and the rest of the ring alike), and the rotor moves by whole
sectors:

```toml
[time]
scheme = "bdf2"
time_step = 0.019634954084936207   # 2 pi / (N * omega), N = 32 sectors
steps = 510

[[sliding]]
zone = "rotor"            # a cell zone of the mesh
axis = [0.0, 0.0, 1.0]
origin = [0.0, 0.0, 0.0]
omega = 10.0              # or rpm = ... (nonzero)
# shift = 1               # sectors per step; the step is shift * 2 pi / (N omega)
```

- The time step is **set by the geometry**: `shift · 2π/(N |Ω|)`. A different
  value is refused with the right one (`SLD003`); so are `end_time` that is
  not a whole number of steps, the adaptive step and the error tolerance.
  Write `steps`, or an `end_time` that fits.
- The count N is found from the mesh (the number of images of the outermost
  rotor cell under the rotation); a mesh that is not a ring of identical
  sectors, or an axis that misses its centre, is refused (`SLD002`).
- The case is transient (`SLD001` for steady ones, see recipe 94 for those),
  with passive scalars and any turbulence model (recipe 100); pressure
  boundaries on the rotor are refused (`SLD004`). The boundaries of the
  rotor must be walls that turn with it (`velocity = ["-spin * y", "spin * x", 0.0]`),
  slip walls, symmetry or wedge planes.
- The solved velocity is the **absolute** one. There is no Coriolis source;
  the interface at the middle of the gap should be a surface of revolution so
  that the mesh flux through it is zero.
- Output times that do not fall on a step fire at the step that reaches them;
  `end_time` and the time schedules are in whole steps.
- A fluid that stands in space is convected relative to the rotor at a CFL
  number near one by an implicit upwind scheme, which smears scalars. Make
  the sector small, or use the frame (recipe 94) when the flow is steady in
  the frame.
- The case file starts from rest and spins up; `[verification] max_error`
  makes the run fail above the band (L2 error 1.35 % on 512 cells).
  `[initial] from_analytic = true` starts a transient run from the
  `[analytic]` velocity and pressure instead (the solid-body check of
  `tests/sliding.rs`).

[Gallery 95](../gallery/flow.md#sliding-mesh); the keys are in
[cases/flow-keys.md](../cases/flow-keys.md) and [keys.md](../keys.md).

<a id="sliding-mesh-imprint"></a>
## 96. Turn a rotor against a stator with a different pitch (sliding mesh by imprinting)

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/taylor_couette_imprint.toml
```

When the rotor and the stator have different numbers of sectors, or the time
step should be free, the interface between them is imprinted in every step
(D182, variant B). The mesh has the rotor and the stator as two meshes that
meet in **two boundaries** on a circular cylinder about the axis (the built-in
annulus mesh makes them with `n_theta_stator`; from your own mesh give a cell
zone for the rotor and two boundary patches):

```toml
[mesh.annulus]
n_theta = 16
n_theta_stator = 24       # 0: one conformal mesh; with a value the interface is split
rotor_layers = 4          # radial layers of the zone "rotor"

[time]
scheme = "bdf2"
time_step = 0.02          # free
steps = 500

[[sliding]]
zone = "rotor"
axis = [0.0, 0.0, 1.0]
omega = 10.0              # or rpm = ...
rotor_interface = "rotor_interface"     # the boundaries of the two sides
stator_interface = "stator_interface"

[solver]
extrapolation = false     # fully implicit in the step
outer_iterations = 4      # 3 to 6 is usual
```

- The two interface boundaries get **no** `[boundary.<name>]` table (`SLD005`
  if they do): their faces are replaced by the imprint at every step. A name
  that is not on the mesh is refused with the nearest one.
- The time step is **free**: `shift` does not apply (`SLD002`). The step must
  be fixed (no `adaptive`, no `error_tolerance`, no sweep over `[verification]`
  levels, `SLD003`), and the scheme must be fully implicit: `extrapolation = false`
  and at least two `outer_iterations` (`SLD003`, with the reason). Two outer
  iterations are enough to keep the error of the case file to 1.35 %, four give
  1.29 %. Courant numbers near 1 relative to the rotor are fine.
- The interface must lie on one circular cylinder about the axis, with the same
  total area on both sides, rotor and stator cells must not share a face or a
  node, and the mesh must be extruded along the axis so that the cells next to
  the interface close (a faces-do-not-close warning, `SLD006`, otherwise).
- Laminar or turbulent flow (recipe 100) and passive scalars; boundaries of
  the rotor are walls turning with it (`velocity = [...]`), slip walls or
  symmetry planes (`SLD004`); one rotor; no restart (`SLD005`, and no
  restart file is written).
- Probes find the cell that is at their point at the time of the sample, and
  the fields file (`vtkhdf`) is written on the mesh of the end of the run with
  the cells at the interface as polyhedra.
- The run builds the solver and the pressure multigrid in every step: on the
  case file the setups are about a quarter of the solver time, on larger meshes
  less (13 % on 2 560 cells, 5 % on 41 000), where the pressure and momentum
  solves take the rest. Use variant A (recipe 95) when the two pitches are equal
  and the step may be one sector.

[Gallery 96](../gallery/flow.md#sliding-mesh-imprint); the keys are in
[cases/flow-keys.md](../cases/flow-keys.md) and [keys.md](../keys.md).

<a id="stirred-vessel"></a>
## 98. Mesh and run a rotor with blades from Gmsh, and draw its relative flow

```
pip install gmsh matplotlib h5py numpy pillow
cargo build --release -p foehncraft-cli
python3 examples/meshes/mixer.py
target/release/foehn flow examples/studies/mixer/mrf.toml --out output/mixer/mrf
target/release/foehn flow examples/studies/mixer/sliding.toml --out output/mixer/sliding
python3 -m tools.plot.rotor_fields mixer.png output/mixer/mrf/fields.vtkhdf \
    --quantities speed,relative,vorticity --omega 10 --rotor-radius 0.06 --outer-radius 0.1
```

(`python3 examples/studies/mixer.py` does all of it and draws the figures of
gallery 98 on the vessel with five blades and three baffles, `python3
examples/meshes/mixer.py --blades 5 --baffles 3` and `--set
mesh.file=../../meshes/mixer5_split.msh`: rotors in published images have
three or five blades, see CONTRIBUTING.md.) What a mesh of your own needs, with `examples/meshes/mixer.py` as
the pattern:

- **A cell zone for the rotor**: a Gmsh physical volume (`rotor`) that holds
  every cell that turns, bounded by a surface of revolution about the axis (a
  cylinder here). The rest can be one or more zones (`stator`).
- **For `[[mrf]]`**: one conformal mesh, rotor and stator sharing the nodes
  of the interface (Gmsh: `fragment` the two surfaces, extrude them in one
  call, `removeAllDuplicates`). The interface faces are interior faces, they
  get no boundary group.
- **For `[[sliding]]` by imprinting**: two meshes with their own nodes on the
  interface (no `fragment`), and the two sides as boundary groups
  (`rotor_interface`, `stator_interface`) that the `[[sliding]]` table names and
  that get no `[boundary.<name>]` table. The node counts on the circle may differ
  (64 and 96 here). The mesh must be extruded along the axis.
- **Walls of the rotor** (blades, hub) are a group of their own with the
  velocity of the rigid turn, `velocity = ["-spin * y", "spin * x", 0.0]` for a
  turn about z through the origin. In a turning zone a wall moves normal to
  itself with the zone, so the normal part of that velocity is used, too (a wall
  that stands keeps only the tangential part of a given velocity, `BC010`).
  Walls of the stator (vessel, baffles) keep the default.
- Write the fields with `output.vtkhdf = true` and `output.derived = ["vorticity"]`;
  for a transient run `output.fields_times = [...]` adds `fields_<step>.vtkhdf`
  at those times, on the mesh of that rotor angle.

`tools/plot/rotor_fields.py` draws the lower face of every cell of a 2D run:
`speed` (|u| with streamlines), `relative` (|u − Ω × r| inside
`--rotor-radius`, with the streamlines of the relative velocity; the stator is
grey), `vorticity`, `pressure`; several fields files make a time series of
the first quantity; `--mesh --window x0,x1,y0,y1` zooms in with the cell
edges; `--palette` writes a small PNG for the docs. The relative velocity is
what a blade sees: judge the rotor by it, the stator by the absolute one.

[Gallery 98](../gallery/flow.md#stirred-vessel); the keys are in
[cases/flow-keys.md](../cases/flow-keys.md) and [keys.md](../keys.md).

<a id="rotor-fields"></a>
## 99. Write the relative flow of a rotor: relative and frame velocity, flow angles, rothalpy

```toml
[output]
vtkhdf = true
derived = ["relative-velocity", "frame-velocity", "turbomachinery",
           "relative-total-pressure", "rothalpy"]
# frame = "rotor"   # optional: every cell in the frame of this zone
```

A rotor is judged by the flow its blades see. With an `[[mrf]]` or a
`[[sliding]]` table these derived fields (D183) are written with the fields,
in the notation of turbomachinery: **c** the absolute (solved) velocity,
**u** = **Ω** × (**x** − **x**₀) the blade speed, **w** = **c** − **u** the
relative velocity.

- `relative-velocity`: **w** and |w|. A cell in a turning zone is seen in the
  frame of its zone, any other in the fixed frame (**w** = **c**): a pump
  shows the relative flow in the impeller and the absolute flow in the volute
  in one field. In ParaView, `relative_velocity_x`, `_y`, `_z` make the vector
  for streamlines (Calculator: `relative_velocity_x*iHat + relative_velocity_y*jHat + relative_velocity_z*kHat`).
- `frame-velocity`: **u** and |u| (zero outside the zones).
- `turbomachinery`: c_r, c_u, c_a, c_m and w_u about the axis of the zone; the
  axis points along **Ω**, so the blade speed and a swirl with the rotation
  are positive. The flow angles in degrees from the circumferential
  direction: α (`flow_angle`) between **c** and **u**, β (`relative_flow_angle`)
  between **w** and −**u**, the convention of pump design (β = 90° when the
  relative flow has no swirl; compare β with the blade angle at the leading
  and trailing edges for incidence and deviation).
- `relative-total-pressure` p + ½ρ|w|², and `rothalpy` I = p/ρ + ½|w|² − ½|u|².
  Rothalpy is constant along a relative streamline of a steady flow without
  losses: where it falls through a rotor, energy is lost there. In the fixed
  parts it is the total pressure over ρ.
- `output.frame = "<zone>"` takes the frame of that zone in every cell (the
  stator seen from the rotor); it must name the `zone` of an `[[mrf]]` or
  `[[sliding]]` table.
- Sample lines get the columns `u_rel`, `v_rel`, `w_rel`, `c_r`, `c_u`, `c_a`,
  `w_u`, and point probes of a transient run `<name>_u_rel`, `_v_rel`, `_w_rel`.
  On the sliding mesh the position is that of the cell at the time of the
  sample (the rotor has turned).
- `foehn post derive` writes the same fields from a restart file of a frame
  run, and the list may change during a run (re-read, D40).

`tools/plot/rotor_fields.py` draws these datasets of a 2D run by name
(`--quantities relative_velocity_magnitude,relative_flow_angle,rothalpy`;
`relative_velocity_magnitude` with the streamlines of **w**).

[Gallery 99](../gallery/flow.md#rotor-fields); the keys are in
[cases/flow-keys.md](../cases/flow-keys.md#derived-fields-d72) and [keys.md](../keys.md).

<a id="turbulent-rotor"></a>
## 100. Run a turbulent rotor: the frame or the sliding mesh with a turbulence model

```toml
[turbulence]
model = "k-omega-sst"            # any model of [turbulence]
wall_treatment = "automatic"

[initial]
turbulence_intensity = 0.05
viscosity_ratio = 10.0

[[mrf]]                          # steady; or [[sliding]] in a transient case
zone = "rotor"
axis = [0.0, 0.0, 1.0]
omega = 10.0

[solver]
algorithm = "coupled"            # steady and from rest at a high Reynolds number
```

- A turbulence model works with `[[mrf]]` and with both sliding meshes
  (conformal stepping and imprinting), and with any number of `[[mrf]]`
  zones, each about its own axis. Nothing else is set: the turbulence
  variables are convected with the flux relative to each zone and produced
  by the strain of the absolute velocity, which is the same in every frame.
- Sliding mesh: the rotor carries k, ε, ω, ζ and f with it, their old levels
  included, so BDF2 stays second order. The wall distance of the turning
  rotor is updated in every step; the log line `wall distance (D48) of the
  turning rotor` says how many cells search the other part (the cells nearer
  to the edge of their part than to its walls; on the stirred vessel 2 118 of
  16 326). With `turbulence.wall_distance = "pde"` the imprinted mesh solves
  the Poisson distance again in every step instead (slower; the exact
  distance is the default). Conformal stepping turns the rotor by whole
  sectors, so its distance does not change.
- From rest at a high Reynolds number, SIMPLEC may diverge in the first
  iterations of a steady frame run (the stirred vessel at Re 20 000 in
  iteration 6, with or without the frame); the coupled algorithm converges.
- Put the first cell of a turning wall at y⁺ ≲ 1 or y⁺ ≳ 25 (`TURB001` warns
  otherwise): in the buffer layer the k of the first cell of the automatic
  treatment depends on the mesh along the wall, also without a frame
  (verification page, section on turbulence on the frame).
- Limits: no restart of a sliding run by imprinting (`SLD005`); several
  sliding rotors need the interface by imprinting (recipe 102); torque and
  power: recipe 101.

`examples/studies/mixer/mrf_turbulent.toml` and `sliding_turbulent.toml` are
the turbulent stirred vessel of gallery 100; `examples/studies/mixer.py
--turbulent` runs both and draws the figures with `tools/plot/rotor_fields.py`
(`--quantities speed,k,viscosity_ratio --nu 1e-6`).

[Gallery 100](../gallery/flow.md#turbulent-rotor); the keys are in
[cases/flow-keys.md](../cases/flow-keys.md) and [keys.md](../keys.md).

<a id="torque-power"></a>
## 101. Read the torque and power of a rotor, and the head and efficiency of a machine

```toml
[[mrf]]                      # or [[sliding]]
zone = "rotor"
axis = [0.0, 0.0, 1.0]
omega = 10.0
walls = ["blades", "hub"]    # default: the walls all of whose faces lie in the zone

[[output.machines]]
name = "pump"
rotor = "rotor"              # default: the only turning zone
inlet = ["inlet"]
outlet = ["outlet"]
gravity = 9.80665            # of the head H = dp_t / (rho g)
```

- The torque (N m) of the fluid on every wall about the axis of its rotor is
  the column `torque` of `boundaries.csv`, the column of the results table and,
  in a transient run, `<wall>_torque` of `forces.csv`; with it
  `<rotor>_torque`, `<rotor>_power` (shaft power −Ω·T), `<rotor>_angle` (degrees,
  sliding mesh) and per machine `<name>_flow`, `_dpt`, `_head`, `_power`,
  `_efficiency`. The end of the log and `flow.md` have a block "Rotors and
  machines".
- The log table gets the columns `rotor.<zone>.torque` and, on the sliding
  mesh, `rotor.<zone>.angle` by default, `machine.<name>.head` and
  `.efficiency` for every machine; `[log] columns` takes
  `rotor.<zone>.torque|power|angle` and
  `machine.<name>.flow|pressure_rise|head|efficiency`.
- The shaft power is in the energy ledger: the rows "shaft work of rotor … in"
  of the tables, the column `shaft_power` of `balance.csv`, the compact line of
  every step. `ROT001` warns when a rotor wall does not turn with its zone.
- The efficiency is hydraulic (flux of total pressure over shaft power; the
  inverse for a turbine), also total to static; the Euler work is ṁ Ω Δ(r c_u).
  Enthalpy rise, isentropic efficiency and pressure ratio of compressors and
  turbochargers come with the energy equation (stage 7, section 13.1).
- An animation of a field series: `fields_every = 10` in `[output]`, then
  `pvpython tools/plot/rotor_movie.py output/<case>` (colour map, isolines,
  colour bar and colours by options, see its header) and the `ffmpeg` command it
  prints.

[Gallery 101](../gallery/flow.md#torque-power); the keys are in
[cases/flow-keys.md](../cases/flow-keys.md) and [keys.md](../keys.md).

<a id="several-rotors"></a>
## 102. Run several rotors at their own speeds, and a transient frozen rotor

```toml
[[sliding]]                  # one table per rotor
zone = "rotor_left"
origin = [-0.08, 0.0, 0.0]   # each rotor its own axis
axis = [0.0, 0.0, 1.0]
omega = 10.0
rotor_interface = "left_rotor_interface"
stator_interface = "left_stator_interface"

[[sliding]]
zone = "rotor_right"
origin = [0.08, 0.0, 0.0]
axis = [0.0, 0.0, 1.0]
omega = -15.0                # its own speed and sense
rotor_interface = "right_rotor_interface"
stator_interface = "right_stator_interface"

[log]
columns = ["time.simulated", "rotor.rotor_left.angle", "rotor.rotor_right.angle",
           "rotor.rotor_left.torque", "rotor.rotor_right.torque"]
```

- Every rotor needs its own cell zone and its own pair of interface
  boundaries **to the stator** (the interface by imprinting, recipe 96);
  conformal stepping takes one rotor only (`SLD002`), and a rotor that
  borders another rotor is refused (`SLD005`). All interfaces are imprinted
  together in every step, at the angle |Ω| t of each rotor.
- The walls of each rotor turn with it: give them their absolute velocity,
  `["-w * y", "w * (x - x0)", 0.0]` for a rotor at `w` rad/s about the z
  axis through `(x0, 0)`.
- The log, `forces.csv` and the report have the angle, torque and shaft
  power of each rotor; the torque is about the axis along the rotor's own
  Ω, so a rotor that drives the fluid has a negative torque and a positive
  shaft power whatever its sense.
- Meshes: `examples/meshes/twin_mixer.py` (a vessel with two stirrers,
  `pip install gmsh`), `[mesh.annulus] outer_rotor_layers = 4` with
  `n_theta_outer` (a second rotor ring at the outer cylinder: zone
  `rotor_outer`, interface `stator_interface_outer` / `rotor_interface_outer`).
- A **transient frozen rotor**: a `[[mrf]]` table in a case with `[time]`
  (D185). The rotor zone keeps its position; the time term is that of the
  velocity in the rotating basis. It is exact for a rotor that is a body of
  revolution and a quick start-up for a bladed one, but it never moves the
  blades past the stator: use the sliding mesh for rotor-stator
  interaction. `[[mrf]]` zones and sliding rotors can share a transient case
  in different zones (`MRF001` for one zone turned by both).

```sh
python3 examples/meshes/twin_mixer.py
cargo run --release -p foehncraft-cli -- flow examples/studies/twin_mixer/sliding.toml
cargo run --release -p foehncraft-cli -- flow examples/cases/flow/couette_two_rotors.toml
python3 examples/meshes/mixer.py
cargo run --release -p foehncraft-cli -- flow examples/studies/mixer/mrf_transient.toml
```

The images of gallery 102 show the twin vessel with three and five blades
and the single vessel with five blades and three baffles (no four-fold
rotor patterns in published images, CONTRIBUTING.md); one script runs and
draws them:

```sh
python3 examples/studies/mixer.py output/mixer --rotors
```

[Gallery 102](../gallery/flow.md#several-rotors); the keys are in
[cases/flow-keys.md](../cases/flow-keys.md) and [keys.md](../keys.md).

<a id="curved-blades"></a>
## 225. Mesh a rotor with curved blades of constant thickness in Gmsh

```sh
pip install gmsh
python3 examples/meshes/three_mixer.py                 # three_mixer.msh, 26 060 cells
cargo run --release -p foehncraft-cli -- mesh examples/meshes/three_mixer.msh
cargo run --release -p foehncraft-cli -- flow examples/studies/three_mixer/sliding.toml
```

`curved_blade(rc)` in `examples/meshes/three_mixer.py` builds one blade
about the origin with the Gmsh OpenCASCADE kernel; `rotor()` turns it to its
place, as for the straight blades. The rules that make it mesh well:

- Make the blade the band between two **concentric arcs** about the centre
  of the camber circle (radii rc ∓ T/2): constant thickness, four curves.
- Close the tip with a half circle built as **two quarter arcs** through the
  nose point; `addCircleArc(start, centre, end)` of exactly 180° is
  ambiguous.
- End the root **well inside the hub** (here at 0.4 of the hub radius) and
  **fuse hub and blades** into one tool before `occ.cut`: a root that ends
  on the hub circle leaves cusps that `foehn mesh` shows as faces above 70°
  non-orthogonality.
- Check the mesh with `foehn mesh` (verdict, the non-orthogonality and
  skewness rows of `mesh_report.md`) before the run.

Change `R_CAMBER` for more or less curvature (a smaller radius bends more;
the bisection in `camber_angle` finds the tip and root on the camber
circle), `T_BLADE` for the thickness; the case file needs no change. For a
forward-curved rotor turn the rotor the other way (`omega = -20.0` and the
blade velocity with `speed_c = -20.0`). The images:

```sh
python3 -m tools.plot.rotor_fields three_mixer_vorticity.png \
    output/three_mixer/fields_000200.vtkhdf output/three_mixer/fields_000400.vtkhdf \
    output/three_mixer/fields_000600.vtkhdf output/three_mixer/fields.vtkhdf \
    --quantities vorticity --limits=-60:60 --dpi 100 \
    --labels "t = 0.1 s,t = 0.2 s,t = 0.3 s,t = 0.4 s"
python3 -m tools.plot.rotor_fields three_mixer_curved_rotor.png \
    output/three_mixer/fields.vtkhdf --quantities vorticity --limits=-60:60 \
    --mesh --window=-0.065,0.065,-0.095,0.035 --dpi 120
```

[Gallery 225](../gallery/flow.md#curved-blades); several rotors in
[recipe 102](#several-rotors).

<a id="revolved-interfaces"></a>
## 103. Give a rotor an interface of any rotationally symmetric shape: discs, cones, pots

```toml
[mesh.revolved]              # a block grid in (r, z) turned about the z axis
r = [0.0, 0.03, 0.045, 0.10] # radial breakpoints; r[0] = 0 reaches the axis
z = [0.0, 0.03, 0.045, 0.055, 0.07, 0.10]
n_r = [10, 5, 12]            # intervals per block
n_z = [8, 6, 4, 6, 8]
rotor = [[0, 1], [0, 3], [1, 1], [1, 2], [1, 3]]   # [r block, z block]
solid = [[0, 2]]             # no cells: the boundary wall_rotor (a disc)
n_theta = 48                 # around the axis: rotor
n_theta_stator = 72          # and stator
cone = 0.0                   # degrees: tilts the inner z breakpoints (cones)

[[sliding]]
zone = "rotor"
axis = [0.0, 0.0, 1.0]
omega = 20.0
rotor_interface = "rotor_interface"
stator_interface = "stator_interface"
```

- One interface per rotor (one pair of boundaries), of any shape that is a
  surface of revolution about the rotor's axis with a meridian of straight
  pieces (D186): a cylinder (radial gap), a disc (an axial gap between two
  blade rows), a cone, or a chain of them. A pot of disc, cylinder and disc
  around an impeller is one interface; so are two discs before and after a
  rotor row that reaches hub and casing. Discs may pass through the axis.
- The meshes of the two sides meet on that surface with their own nodes and
  their own pitch, and need a ring of nodes at every kink of the meridian
  (each face lies on one straight piece). A meridian that branches or closes
  on itself (a rotor ring that touches nothing) is refused. Curved
  meridians are later work; fully general non-matching interfaces are D58
  (stage 7).
- The log line of the sliding mesh names the shape it found ("a disc at
  z = …", "a cone from …", "a surface of revolution of 3 straight pieces in
  1 part, 2 through the axis").
- `[mesh.revolved]` names the boundaries `inner` (when r[0] > 0), `outer`,
  `zmin`, `zmax` and `wall` (around `solid` blocks); on rotor cells they get
  the suffix `_rotor` (walls that turn with the rotor: give them the
  velocity `["-w * y", "w * x", 0.0]`). The interface is always split into
  `rotor_interface` and `stator_interface`.
- Where the interface meets a wall (a hub), the rotor's and the stator's
  polygons of that wall differ by a small step: an error of second order in
  the pitch (0.2 % of a uniform flow at 24 and 34 intervals). Away from
  walls the imprint keeps a uniform flow to round-off level.
- Meridional images: `tools/plot/meridional_fields.py` averages the rings of
  a mesh of revolution and draws (r, z), mirrored.
- Wedge cells at the axis are thin around it: at high viscosity their
  coupling to the neighbours around the axis is many times their time term.
  With outer iterations PISO handles this since PR #106 (SIMPLEC-type
  correctors, `solver.corrector`); with fewer
  intervals around the axis (`n_theta`) the cells are less thin and the
  solves cheaper.

```sh
cargo run --release -p foehncraft-cli -- flow examples/studies/revolved/disc_pot.toml
cargo run --release -p foehncraft-cli -- flow examples/studies/revolved/swirl_stage.toml
cargo run --release -p foehncraft-cli -- flow examples/cases/flow/axial_stage_uniform.toml
cargo run --release -p foehncraft-cli -- flow examples/cases/flow/pot_solid_body.toml
cargo run --release -p foehncraft-cli -- flow examples/cases/flow/pot_solid_body_viscous.toml
pip install matplotlib h5py numpy
python3 -m tools.plot.meridional_fields disc_pot.png output/disc_pot/fields.vtkhdf \
    --quantities swirl,meridional --palette \
    --labels "swirl u_θ|meridional flow (u_r, u_z) with streamlines" \
    --interface "0,0.03,0.045,0.03;0.045,0.03,0.045,0.07;0.045,0.07,0,0.07"
```

[Gallery 103](../gallery/flow.md#revolved-interfaces); the keys are in
[keys.md](../keys.md).

---

<a id="tjunction-mixing"></a>
## 196. Mix two streams in a T-junction and check the mesh

Mesh the junction with Gmsh at a size `h` (Netgen optimisation breaks the
conversion into polyhedra, so leave it off), convert on reading and run:

```
gmsh -3 examples/meshes/tjunction_gci.geo -setnumber h 0.0079 -format msh41 -o tj.msh
target/release/foehn flow examples/cases/flow/tjunction_mixing.toml --out output/tj \
    --set 'mesh.file="tj.msh"'
```

`polyhedral = true` makes one cell per node. The branch carries a
`[scalars.theta]` value of 1, the main inlet 0 (recipe on scalars). For
the grid study run `python3 examples/studies/tjunction_gci.py output/tj_gci`
(three meshes, about half a minute); check the report for the observed
order of each quantity before trusting a number. A run that stops with
GEN001 has a boundary cell whose centre lies behind its face: take
another mesh size, or run the tetrahedra (`polyhedral = false`).

[↑ back to top](#incompressible-flow)
