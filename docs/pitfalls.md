# Pitfalls: lessons and experiments not to repeat

This page keeps what earlier work found out the hard way (decision D200),
so that nobody spends hours of runs finding it again. The rules that follow
from it are in [AGENTS.md](../AGENTS.md). Add a lesson here in the pull
request that learns it: what was tried, what it measured, and what to do
instead.

## Contents

- [Steady runs: budget and stopping](#steady-runs)
- [Coupled solver: experiments not to repeat](#coupled-solver)
- [Skewed Gmsh tetrahedra](#skewed-tets)
- [Observers must not change the run](#observers)
- [Timing small meshes](#timing)
- [Iteration counts in tests](#counts)
- [Counts a file declares](#declared-counts)
- [Reference numbers of other solvers](#other-solver-numbers)

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
