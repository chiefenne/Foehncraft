# Pitfalls: lessons and experiments not to repeat

This page keeps what earlier work found out the hard way (decision D200),
so that nobody spends hours of runs finding it again. The rules that follow
from it are in [AGENTS.md](../AGENTS.md). Add a lesson here in the pull
request that learns it: what was tried, what it measured, and what to do
instead.

## Contents

- [Steady runs: budget and stopping](#steady-runs)
- [Coupled solver: experiments not to repeat](#coupled-solver)
- [Multi-rank: local geometry must be canonical](#multi-rank)
- [Skewed Gmsh tetrahedra](#skewed-tets)
- [Observers must not change the run](#observers)
- [Timing small meshes](#timing)
- [Iteration counts in tests](#counts)
- [Counts a file declares](#declared-counts)
- [Memory of a mesh reader](#reader-memory)
- [Node order of VTK's wedge](#vtk-wedge)
- [Arrays of a blocked VTKHDF file](#blocked-vtkhdf)
- [Reference numbers of other solvers](#other-solver-numbers)
- [Coupling loops: the solid must see the flow's response](#coupling-loops)
- [Limiters on tetrahedra and polyhedra](#limiters)
- [Porous resistances and the Rhie–Chow factor](#porous-rhie-chow)

<a id="steady-runs"></a>
## Steady runs: budget and stopping

- **Check at 2 000 to 3 000 iterations.** A steady run that has not settled
  by then has a reason (relaxation, tolerance, mesh, a model switch);
  diagnose it there and change the settings or stop. More iterations do
  not fix it and waste CI time. Never let a run go to 10 000 or more.
- **The residual is no stop measure.** On the backward-facing step the
  scaled residual sits at 5e-5 to 2e-4 when the reattachment length has
  settled, so a tolerance of 1e-6 is never reached (gallery and recipe 74).
- **The other direction too: a small residual is no proof.** In the
  turbulent channel the residuals passed 1e-6 while the bulk velocity
  U_b+ still moved: SST stopped 2.9 % low on the fine mesh, EB k-epsilon
  1.5 % high. The SST omega residual sits at 1e-11 to 1e-13 (its scale is
  dominated by the fixed near-wall rows) and never decides anything.
  Grid-convergence studies of wall-resolved steady runs need tolerance
  1e-9, or a check of the monitored quantity. An earlier study at 1e-6 was
  wrong by up to 3 %.
- **Judge settling by a monitored value** (reattachment length, torque,
  bulk velocity, Nusselt number), read from monitors or intermediate
  results.
- **Relaxation ramp for separated SIMPLEC cases:** momentum relaxation 0.5,
  then 0.8 from iteration 100 (`[[schedule]]`). With 0.5 throughout the
  count doubles; 0.9 stalls the residual at 1e-4. The count roughly
  triples when the spacing halves.
- **Uniform start velocity into a solid face** (for example a heat sink
  started with the inflow speed everywhere) upsets the first SIMPLEC
  iterations. This is an open finding; until it is fixed, start such cases
  from rest, as the heat-sink case does.

<a id="coupled-solver"></a>
## Coupled solver: experiments not to repeat

The coupled pressure-velocity solver (step 22,
[docs/requirements/step22-coupled-solver.md](requirements/step22-coupled-solver.md),
section 7) was measured in detail. Its speed work is finished; ask the
maintainer before starting more.

- **Schur-complement variants: do not try again.** SIMPLER (Klaij and Vuik
  2013) saved a little Krylov work but needs a second pressure solve and
  was slower (21.2 s against 16.3 s on the coarse step). The LSC commutator
  (Elman et al. 2006) sat at the Krylov limit and did not converge in 200
  iterations (69.3 s), a known boundary weakness. The cause is spectral: no
  diagonal approximation of the momentum inverse fits both smooth and
  rough modes (section 7.9). Only a coupled block multigrid would remove
  this.
- **The cap alpha_u = 0.98 stays.** With 0.999 the fine step settles in
  fewer iterations but every Krylov solve hits its limit; run time is
  about four times longer.
- **Turbulence relaxation 0.98 with the coupled solver.** 1.0 diverges
  k-epsilon in iteration 3; 0.9 is several times slower.
- **`[solver.turbulence] preconditioner = "amg"` fails with the coupled
  algorithm** (coarsest level not positive definite). The reduction
  multigrid (D180) is the one to use.
- **Algebraic multigrid choices for the momentum sweep were settled in the
  lab** (D179, section 7.7): smoothed-aggregation AMG was useless for the
  momentum block, reduction (AIR) multigrid works. Do not re-run that
  comparison.
- **Pressure preconditioner tolerance 0.3 and hierarchy reuse at 5 %
  change** (D181) are the measured optimum of 0.02, 0.05 and 0.2.
- **Small laminar cavities do not gain** from the coupled solver; that is
  expected, not a bug.
- **Measure a phase before naming it the limit.** An early note blamed the
  turbulence BiCGStab iterations for 42 % of the run; measured, the
  transport solves were 13 %.

<a id="skewed-tets"></a>
## Skewed Gmsh tetrahedra

- **Raw Gmsh tetrahedra have 0.6 to 0.8 % faces with skewness above 0.9**
  (up to 1.3), mostly tets with two wall faces. The least-squares pressure
  gradient on those faces drives divergence. Optimise the mesh with Netgen
  (`Mesh.OptimizeNetgen = 1;` or `gmsh -optimize_netgen`); `Algorithm3D`
  alone does not help. Grading the mesh from the boundary
  (`Mesh.MeshSizeExtendFromBoundary = 1`) also lowers the worst skewness.
- **Local fixes do not converge.** First order, local Green-Gauss or
  damping at the skewed cells alone were tried and failed; the cure is
  global. `solver.remediation = "auto"` (D50, the default) switches to the
  Gauss pressure gradient with a limiter and to the bounded Rhie-Chow face
  velocity (`solver.face_velocity = "bounded"`) when a face is skewed
  beyond `solver.skewness_limit`.
- **Not the cause on the stirred tank** (ruled out by switching terms off):
  the turning shaft, the pseudo-time step, the Coriolis term, the Reynolds
  number, the non-closure of the MRF interface (MRF003), non-orthogonal
  correction, a stronger limiter. The cause was the unlimited gradient
  shift of the face velocity on skewed faces. Reproducer:
  `examples/studies/skewed_tets.py`.
- **A home-made mini mixer with thin baffles** had more than 300 faces
  above 70 degrees non-orthogonality and diverged for another reason; it is
  not a representative test.

<a id="observers"></a>
### Polyhedra from Netgen tetrahedra

Converting tetrahedra into polyhedra (`mesh.polyhedral = true`) can end the
run with GEN001 (the centre of a boundary cell lies behind its face, found
on the T-junction of gallery 196, 2026-11-22): every mesh of the geometry
optimised by Netgen failed, plain Gmsh meshes converted on three of four
sizes (h = 0.0050 m failed, 0.0055 m passed). Use plain Gmsh tetrahedra,
try another size, or run the tetrahedra; the cause in the dual of boundary
nodes at the concave junction edge is not fixed.

### Sharp wall corners and thin wall cells (the asymmetric diffuser)

The steady SST diffuser of gallery 199 looked as if it did not settle on
fine meshes (2026-11-22). Four causes, each found by switching one thing:

- **A sharp corner of the wall breaks the least-squares gradient of the
  thin wall cells.** Across the corner the wall cell next door has the same
  wall distance but lies lower by about Δx tan θ / 2 (θ the kink angle), that
  is A tan θ / 2 wall-cell heights for a wall cell of aspect ratio A: 13
  heights at A = 150 and θ = 10°. The large wall-normal gradient leaks into
  ∂u/∂x, which is wrong by O(1) on every mesh with the same cell shapes (the
  face extrapolation error falls only as Δx). Linear upwind without a
  limiter from the first iteration on extrapolates it: the velocity of the
  wall row grows downstream by a factor 2 to 10 per iteration and the fine
  meshes diverged in about 10 iterations. Three remedies, each run on the
  diffuser (2026-11-22):
  - start with upwind for 200 iterations, then linear upwind
    (`[[schedule]]`): the sharp diffuser with Δx = 0.25 H (A = 150 at the
    kink) then converges on levels 0 and 1 (494 and 1 378 iterations);
  - shorten the wall cells towards the corner: A about 11 at the kink
    (`diffuser.py --sharp`, case `diffuser_sst_sharp.toml`), one wall-cell
    height of offset instead of 13, converges on levels 0 and 1 (1 087
    and 2 699 iterations) and reaches residuals of 4·10⁻⁷ on level 2 at
    the cap of 8 000 iterations, with the probes steady from iteration
    4 000;
  - round the corner where the geometry has a radius (the experiment has
    4.3 H): the error then falls as Δx² (0.10 and 0.017 of the velocity on
    levels 0 and 1, against 0.55 and 0.27 with the sharp corner).

  Upwind or the Venkatakrishnan limiter alone only hide the error. Behind a
  sharp corner the wall shear stays near zero for about 1.5 H, and small
  bubbles there come and go with the mesh: take the longest stretch of
  reverse flow as the separated region, not the first sign change.
- **Trapezoidal wall cells look skewed.** Where a channel widens, the wall
  cells (aspect ratio 300 to 600) are trapezoids whose centroids shift
  along the wall by more than their height; the face skewness (over the
  centroid distance) reaches 1.8 on Δx = 0.5 H and switches the D50
  remediation on, on the coarse levels only, so the levels of a grid study
  run with different numerics. Make Δx small enough (0.25 H here) that the
  skewness stays below 0.9 on every level (now at most 0.42, 0.25 and 0.13
  on the three levels).
- **A uniform start breaks continuity in a widening.** u = 1 everywhere
  carries 4.7 times the inflow out of the diffuser; the finer meshes
  diverged in iteration 6. Start from potential flow, and from upwind for
  the first 200 iterations (`[[schedule]]`), then linear upwind: the
  potential flow slips along the walls, and linear upwind on that start
  diverges on the coarse mesh.
- **SIMPLEC at relaxation 0.5 converges a large, weak separation very
  slowly.** The reverse flow is about 2 % of the inflow; at 0.5 the bubble
  was still growing after 1 700 iterations (reattachment 21.7 H) and the
  residual stalled near 5·10⁻⁵, which looked like a missing steady state.
  At 0.8 the same mesh converges to 10⁻⁷ with the reattachment at 29 H. A
  transient run from the stalled state moves the same way, so check with
  one before concluding that a steady state does not exist.

### Sharp 180° turns need a gentler start (the U-bend)

The U-bend of gallery 204 with a square divider end and a square end wall
(`ubend_sst_sharp.toml`, square cells at the kinks) did not run with the
settings of the rounded bend (2026-11-22):

- relaxation 0.8 from the potential-flow start diverged in 10 iterations
  in the cells at the concave corners of the end wall, where the potential
  flow stagnates; start at 0.5;
- switching to linear upwind and to relaxation 0.8 in the same iteration
  (201) diverged behind the divider end; switch the scheme first (601),
  the relaxation later (1 201);
- relaxation 0.5 throughout converges, but the separated region behind the
  divider end grows very slowly: on 50 400 cells it was still growing at
  6 000 iterations (residual 2·10⁻⁵). With the raise to 0.8 the same mesh
  converges in 5 964. Watch the probe in the bubble, not the residual.

## Observers must not change the run

Writing an intermediate result or a balance table between two outer
iterations used to update boundary values and gradients and leave them in
place. Both are not idempotent, so the next iteration started from a
different lagged state (momentum residual 7 % off on a small duct). Any
new observer (monitors, output, fields of later stages) wraps its work in
`Flow::lagged` / `Flow::set_lagged`, and `Lagged` is extended when a new
equation brings its own boundary values or gradients. The test
`intermediate_results_and_balances_do_not_change_the_run` checks this bit
by bit.

<a id="timing"></a>
## Timing small meshes

A thread pool of one thread used to hand every kernel launch to the pool
(about 2 000 hand-offs per outer iteration). It is fixed, but every run
time written down before that fix (2026-11-18) is too long; compare only
timings taken after it, on a quiet machine. Profile with `strace -c -f` and
`valgrind --tool=callgrind`.

<a id="counts"></a>
## Iteration counts in tests

Iteration counts repeat bitwise for the same input but react to the last
digit of an input (one ulp in the pseudo-time step: 423 against 452
iterations on the cavity), because of loose linear tolerances, truncated
sweeps and single-precision multigrid cycles. Quote counts to about 10 %
and never assert exact counts in tests.

<a id="declared-counts"></a>
## Counts a file declares

A mesh reader that trusts a count of the file (`NELEM= 99999999999999`,
`(10 (0 1 ffffffffffff 0 3))`, a Gmsh node tag of 2^63) and calls
`Vec::with_capacity`, `resize` or `vec![x; n]` with it aborts the whole
program: a failed allocation is not a panic and cannot be caught. The
broken-input corpus (D202) found this in four of the six readers. Check
the count against the size of the file first (`format::declared`: every
item takes one byte at least), map sparse tags with a hash map instead of
a table up to the largest tag, and compute `last - first + 1` with checked
arithmetic (tests run in debug, where an overflow panics).

<a id="reader-memory"></a>
## Memory of a mesh reader

A reader that loads the file whole and parses it into general containers
needs many times the size of the mesh. The VTK reader before PR #245 read
a 1.2 GB ASCII `.vtu` of 1 M polyhedra with a peak of 5.6 GB: the file, a
copy without the appended data, the text of every element again, every
number as `f64`, one `Vec` per polyhedron face in a hash map, a clone of
those for the builder, then a sorted copy of every face to match them. A
10 M-cell mesh would have needed some 56 GB. Read through a buffer, put the
numbers straight into the arrays the mesh keeps (`u32` node ids, flat face
streams, see `format::Polyhedra`), move arrays instead of cloning them,
drop each input as soon as it is used, and match faces in groups of their
smallest node (`builder::match_faces`) instead of sorting copies. The
peak is then about twice the finished mesh (1 GB per million
polyhedra). Measure the peak resident memory, not only the time
(`/usr/bin/time -v` or `-l` on macOS).

<a id="vtk-wedge"></a>
## Node order of VTK's wedge

Older VTK documentation and face tables of `vtkWedge` had the base
`(0, 1, 2)` pointing away from `(3, 4, 5)`, the opposite of the
hexahedron. VTK 9.7 lists `(0, 2, 1)` as the outward base, so `(0, 1, 2)`
points towards `(3, 4, 5)`, as for the hexahedron, the tetrahedron and the
pyramid; its parametric coordinates, its Jacobian and the volumes of
`vtkCellSizeFilter` and `vtkMeshQuality` all agree. FoehnCraft's VTKHDF
writer has always written that order; a "fix" proposed from the old
documentation would have turned every prism inside out in ParaView. Check
a node order with a volume computed by VTK itself
(`tools/check/vtkhdf.py`) or against VTK's `GetFace` table, as
`vtkhdf::tests::standard_cells_in_vtk_order` does, not with a comment.

<a id="blocked-vtkhdf"></a>
## Arrays of a blocked VTKHDF file

`fields.vtkhdf` is written in blocks (layout 1, shared with MurmelCraft
and KoflCraft): partitions of 65,536 cells along a Hilbert curve, each with
its own points and a layer of ghost cells that repeat the cells of the
neighbouring blocks with their values. Thirteen plot tools and examples
read `Points`, `Connectivity` and `CellData` of `internal` whole; on a
file of several blocks that indexes the points of one partition with the
local ids of another, counts every ghost twice in a sum or a volume
integral, and the cells no longer come in the order of their ids even in
one block. Read the arrays through `foehncraft.vtkhdf.internal` (the
tools through `tools.plot.vtkhdf`), which drops the ghosts and returns the
cells and points in the order of their identities, the arrays of the
single-piece file before; check a reader on a file of several blocks
(`vtkhdf::write_blocks` with a small block size), not only on a small case
that fits in one.

<a id="other-solver-numbers"></a>
## Reference numbers of other solvers

A value computed by another code (a benchmark page's mean of two solvers,
the authors' own model result in a paper, the velocity planes of a paper)
is not a reference (D201). Two cases showed how it creeps in: a model
constant of the automatic wall treatment was fitted so that the flat plate
met such a value (6g), and a gallery entry drew its figures against another
code's planes. Use measurements, DNS (published scalar values with their
citation), correlations and analytic solutions, and our own mesh sweep
(Richardson value, GCI) for verification. When a constant is fitted, fit it
to one of those and say which. The audit of 2026-11-21 and its rework of
2026-11-22 are in `docs/verification/DATA_SOURCES.md` (rule 5) and the
changelog fragment 176.

<a id="multi-rank"></a>
## Multi-rank: local geometry must be canonical

Found by the MPI spike of 2026-11-22
([assessment](assessments/mpi-multi-rank-sizing.md)). The local mesh of a
rank (`decompose::extract`) re-orients faces so that the owner is the lower
local cell and sums every cell over its faces in local face order. Geometry
and matrix rows then differ in the last bits from the one-rank run for about
40 % of the entities, and the halo cannot repair that, because the face
quantities are recomputed on each rank. **Do instead:** evaluate every face
with the node order and orientation of the global mesh (sign flip only),
sum every cell over its faces in ascending global face id, copy halo values
from the owner. Done in MPI step W1 (#200): `Mesh::face_flip` marks the
flipped faces, `extract` lists cell faces in global order, and
`Geometry::compute` is canonical on any local mesh; `foehn mesh` and the
example `mpi_matrix` (1 to 4 ranks, four cell types, in CI) check it.
**Flipping a face is not enough for the solver** (W3, #230): exact for the
geometry (a sign), but the face formulas are asymmetric (the owner weight w
against 1 − w, the upwind choice at a zero mass flux), so a flipped face
gave 146 of 210 cells a different φ on one rank. Since W3 `extract` keeps
every face in its global orientation (a processor face may be owned by the
halo cell), and every face is computed on both of its ranks from the same
inputs; `Mesh::face_flip` is all false there. Do not flip faces of a
rank's part again; code that needs the own cell of a processor face takes
the one of owner and neighbour below `n_cells`.
Do not compare partitions by tolerance when a bitwise comparison is
possible: a tolerance hides the cause (D93).

**A decision taken on the local mesh** (W5, #234). The flow solver picks a
pressure reference cell when no boundary fixes the pressure; on ranks a
part without the outlet decided "no pressure boundary" on its own and
pinned a cell, and every field of every rank differed from one rank. The
same holds for the mesh skewness that switches gradient corrections, and
for any choice made from what a rank sees. **Do instead:** reduce such a
decision over all ranks (`any`, `comm.max`) before it is used, and agree
on anything that steers a collective (for example the kind of a device
buffer, which decides whether it is exchanged) with `comm.min`, so that
no rank skips an exchange the others wait for. Start fields from a global
construction (potential flow, smooth fields) are computed on the global
mesh and each rank keeps its cells (`local.cell_position`).

**Sums over ranks in rank order.** A pairwise or block sum per rank,
combined in rank order, changes in the last bits with the rank count (0 to
2 ulps on the 1 280-cell box of gallery 219, for each of volume, first
moment and boundary area). **Do instead (W2, #219):** accumulate into
`foehncraft_core::exact::ExactSum` (or `.exact_sum()` on an iterator) and
reduce with `Communicator::sum_exact`; never `f64` partial sums. Several
sums go in one all-reduce with `sum_exact_many`. For min and max use
`min_total` and `max_total` (NaN propagates, −0 < +0), not `f64::min`,
whose result for ±0 depends on the order.

**Collective restart maps (W6, #262).** The gather and scatter of a
restart communicate inside `Restart::map_fields`, so every rank must map
the same arrays in the same order; a length check that fails on one rank
only would leave the others waiting in the all-gather. Agree such checks
with `comm.min` before the collective. The momentum sweep state keeps
`matrix` per hierarchy and `slot` per velocity component: mapping them
under one index dropped a hierarchy and made the gathered coupled-solver
file differ from one rank in `momentum_amg_use`.

**Runs on ranks take every branch together (W7, #279).** Every rank runs
the loop of the case, so a decision made from rank-local state parts the
ranks and leaves them waiting in different collectives: the wall clock of
`[run.limits]`, control files (rank 0 reads the output directory, the
others a scratch directory), signals and the live control port. Rank 0
decides and `OnRanks::agree` hands its bits to all ranks before the
branch. Output that reads the solver of a rank directly (the scalar and
energy summaries, the turbulence wall summary and floor counters, the
remediation counts) sees only its part: read it from the gathered
`Solution` or from the global discretisation, never from `fl`.

**A count over all ranks needs a total over all ranks (W7 part 2, #290).**
The device reductions on ranks are global, so a count such as the floored
cells of the turbulence model is the count of the whole mesh; comparing it
with `n_cells` of the solver (the rank's own cells) made the floor
statistics (TURB001) differ from one rank. Measure a global count against
`comm.all_sum_i64` of the own cells. The same holds for any summary that
loops over device arrays in cell or face order (heat balance, y⁺ per
wall, ν_t/ν max and its cell): gather the arrays with `OnRanks::fields`
and run the one-rank loop on the whole mesh, since plain `+=` sums in
another order give other bits even where the set of values is the same.

<a id="coupling-loops"></a>
## Coupling loops: the solid must see the flow's response

The first build of the cycle-averaged loop (D204) gave the solid phase the
Robin pair of the first fluid cell, as the classic mapped loop does. It
converged no sooner than the fully coupled run (gallery 190): the
first-cell conductance (205 W/K there) is the answer of the heat flux with
the fluid held, while the flow answers with about ṁ c_p and h A (8 and
30 W/K), because the near-wall fluid follows the wall and the bulk warms.
Aitken, IQN-ILS, a per-face and a global secant on top did not fix it (the
secants overshot by hundreds of kelvin). What worked: the solid phase
solves fluid and solid together, steady, with the cycle-mean mass flux, and
a per-face defect on the solid cells brings the heat flows to the measured
cycle means (3 to 4 updates). The defect goes into the solid cells only;
put into the fluid cells it changed the fluid's answer (75 W instead of
160 W). Before building a partitioned loop, compare the interface
conductance the partner sees with the true response of the other side.
The classic loop with an h against the bulk temperature (about h A) is
much closer than the first-cell h and converges in a few passes where
ṁ c_p is large against h A; it is slow where the bulk follows the wall.

Three more from the skin swing (#202):
- **A lagged swing oscillates.** The swing of the held wall changes the
  flux it is computed from; updated after the fact it swung between 0.08
  and 0.49 K per update. Its amplitudes belong in the unknowns of the
  quasi-Newton update next to the wall temperatures.
- **Release into fully coupled cycles from the full skin profile.** Freed
  from the mean, the solids lack the energy of the swing in the skin, and
  that error leaves with the time constant of the whole wall, not of the
  skin (0.3 K after three cycles); started from the damped and delayed
  swing in depth it was 0.016 K.
- **Material values, not the face values, before the first step.** The
  conductivity and c_p on the faces of the energy equation are the
  fluid's until the first refresh; an effusivity taken from them at the
  start was 12 % off on a wall with k_s = k_f, and wrong by the ratio of
  conductivities on any other.

Three more from the frozen-flow march (#208):
- **A settling test with extrapolation never passes for a flow that
  follows the temperature.** The flow drifts with the warming fluid, so
  the ratio of two window changes stays near 1 (0.98) and the change
  still to come, c r/(1 − r), stays large. For such a flow the change of
  the last window alone decides; the refreshes catch the drift.
- **A restart of a conjugate run must carry μ(T).** The viscosity on the
  faces of a flow that follows the temperature is lagged state; without it
  in `restart.h5` a restart in a refresh was 5·10⁻⁵ K off instead of bit
  for bit.
- **A viscosity written as an expression is not seen as a function of
  T.** Only the property forms (power law, polynomial, table, ...) tell
  the march that the flow follows the temperature; an expression string
  of T made the march hold the flow once and never refresh it. Use a
  form.

<a id="limiters"></a>
## Limiters on tetrahedra and polyhedra

Found with the vertex limiter of D237 (#211):
- **A face limiter does not keep a scalar bounded on tetrahedra.** Linear
  upwind with Venkatakrishnan or Barth–Jespersen keeps a face value in
  the range of the face neighbours of the upwind cell, not of the cell it
  flows into. The air of the heat exchanger reached 12.6 K above its
  inlet temperature with Venkatakrishnan and 44 K with Barth–Jespersen;
  a bounded oblique step on polyhedra stays bounded with Barth–Jespersen,
  so one mesh type proves nothing. Look at the extremes of the report
  against the boundary values on the real mesh.
- **A strictly bounded limiter stalls a steady run.** With ε = 0, the
  Barth–Jespersen function at the vertices or a per-cell ε that shrinks
  with the mesh, the limiter value of a few cells in a mixing layer
  switches every outer iteration and the scalar residual stops between
  10⁻⁵ and 10⁻³ while the flow has converged (T-junction). Keep ε at
  about 1 % of the range and relax ψ (`limiter_relaxation`).
- **A limiter that wins on bounds can lose on a grid study.** The vertex
  limiter clips smooth extrema within ε: the T-junction's scalar
  difference on the coarse mesh drops by a third, and the two cases with
  an exact smooth solution move away from it. Do not switch a default on
  a sharp-front case alone; run the smooth ones too.
- **A limiter with memory belongs in the restart file.** The relaxation
  of ψ makes the limiter value state: before #214 a restart started from
  ψ = 1, and in its first iteration the exchanger's air inlet nozzle went
  60 K above the inlet temperature (it healed within 150 iterations).
  Every quantity that carries from one outer iteration to the next goes
  into restart.h5, with a bit-identical restart test; do not take images
  or results from the first iterations after a restart of an older file.

<a id="porous-rhie-chow"></a>
## Porous resistances and the Rhie–Chow factor

Found with the porous regions of D270 (#273):
- **A scalar Rhie–Chow factor cannot carry a turned resistance tensor.**
  The factor D of the face fluxes is one number per cell. A honeycomb with
  its channels turned 30° against the mesh (resistance 10⁴ times larger
  across than along) stalled at a continuity residual of 5·10⁻³ whatever
  the relaxation; with the smallest principal value in D the pressure
  correction does not see the blocked direction, with the largest it does
  not see the open one. Axes along the mesh axes converge. Turned axes need
  a tensor (or face-projected) factor; until then they are refused
  (POR002), do not try a scalar compromise again.
- **Defer the off-diagonal part of the tensor.** Put Σ|a_ij| of the
  off-diagonal resistance on the diagonal and the off-diagonal product
  on the right side (deferred correction); the plain explicit
  off-diagonal term with a 10⁴ blocking ratio broke BiCGStab (NaN).
- **A block that closes the channel converges slowly.** Channels across
  the main flow (the block an obstacle) left the x-momentum residual near
  0.8 after 2 500 iterations while the pressure drop was settled; judge
  such runs by the monitored drop, not the residual.

<a id="chemistry-mechanisms"></a>
## Reaction mechanisms and their units

- **Check the units of a global mechanism against a second source.** The
  rate constants of Jones & Lindstedt (1988) are in m, kmol, s; read as cm,
  mol, s they are 178 to 1000 times too slow and the mechanism does not
  ignite. A later paper that tabulates them in cm, mol, s (Si et al. 2020,
  table 1) settles it; the unit test `the_built_in_mechanisms_load` keeps
  both numbers side by side.
- **A global mechanism need not end at its own equilibrium.** The 2-step
  Westbrook–Dryer mechanism dissociates CO₂ by a reaction of its own, not
  as the reverse of the CO oxidation, so its reactor ends far above the
  equilibrium temperature (277 K from 1400 K). That is the published
  mechanism, not a bug of the integrator: check an integrator with
  GRI-Mech 3.0, whose reversible reactions stop at the Gibbs equilibrium.
- **Finite differences of the Jacobian need an absolute floor in the
  step.** A species at zero (AR in the global mechanisms) with a step
  relative to its value gives a step of zero or round-off size; use
  `1e-6·|z| + 1e-9`.
