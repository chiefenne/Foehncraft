# Step 26 – Design optimisation, ML surrogates and immersed boundaries

Requirements step 26 (2026-11-21), maintainer request: a design-optimisation
strategy for expensive, agent-driven CFD. Decisions D205–D212, L49 and the
amended L13 are in [requirements.md](../requirements.md). No code. Clean room (D201): literature only, no copyleft
source, no commercial solver names. Estimates are agent time.

Next free ids at the time of writing: D205, L49, P19, O25, step 26.

## 0. Where we stand (checked in the repository)

- **Adjoint.** What exists is the *adjoint-ready assembly* of D59: kernels are
  pure functions, they run on forward dual numbers and agree with central
  differences (`crates/foehncraft-flow/tests/adjoint.rs`), matrices can be
  applied transposed without atomics (`foehncraft-kernels/src/adjoint.rs`),
  and transient runs can checkpoint. The **discrete adjoint solver itself
  (L17) is not built yet**: no adjoint equations are solved, no surface
  sensitivities, no mesh morphing. Everything gradient-based below depends on
  L17. This is the honest correction to "the existing adjoint is a given":
  the foundation is a given, the solver is the first job.
- **Studies.** D62 (`foehncraft.study`: grids, Latin hypercube, Sobol, warm
  start from the nearest finished run), D105 (stepping API, state in and
  out), D106 (training-data export), D167 (ensemble Kalman inversion).
- **Agents.** D196 tool registry (PR 1 merged, #164), D197 journal, D198
  safety and budgets, D199 limitations (PR 1 merged, #165).
- **Learning.** D165 acceptance of learned components, D168 licences, P14
  "surrogate-assisted studies" was a proposal on the step 20 page and never
  adopted.
- **Immersed boundaries.** Already in the requirements as **L13** (maintainer
  2026-11-14): "cut cells or ghost cells on the polyhedral face-based mesh,
  no change to the solver data structures; scheme, accuracy targets and test
  cases are decided when the item starts". D108 part 3 (winding numbers,
  signed distance to STL) is reserved for it. So the method is still open:
  section 3 recommends one, the cards ask the owner.
- **pyaeroflow.** 2D unstructured cell-centred FV in JAX, incompressible
  (artificial compressibility, SA, SST, transition) and compressible steady.
  Being JAX, it is differentiable end to end, which makes it an unusually
  good cheap level (section 2.3).

## 1. Principle: the agent orchestrates, the optimiser optimises

An LLM agent is a poor optimiser (no reliable arithmetic over hundreds of
iterations, expensive per step, non-deterministic) and a good orchestrator
(sets up cases, picks a strategy, reads limitations, decides when a result
is not trustworthy, explains). So:

- the optimiser is deterministic code (SQP, interior point, trust region,
  Bayesian optimisation) running inside one tool call or one study;
- the agent chooses the formulation, the budget and the fidelity ladder,
  watches the run record, stops on limitations, and writes the report;
- every evaluation is a D196 tool call, journalled (D197) and replayable
  with `foehn run --script`, under the D198 budget.

## 2. The approved strategy, as requirements

### 2.1 A differentiable chain of typed tools (proposed D205)

Idea from simultaneous analysis and design and from Sharpe's AeroSandbox
work (P. D. Sharpe, MIT SM thesis 2021, "AeroSandbox: a differentiable
framework for aircraft design optimization"): every stage of the chain is a
function with a value and derivatives, so an outside optimiser can call it
like a symbolic-framework function.

Chain: design parameters p -> geometry (parametric CAD or free-form
deformation box, or a level set for immersed bodies) -> mesh morphing x(p)
(L4) -> flow state u(x) (primal solve) -> objective J(u, x) and constraints
c(u, x) (monitors of D43, integrals, ledger terms of D88).

Each stage is a registered tool with typed input and output and three
entry points: `value`, `jvp` (forward, from D59 dual numbers) and `vjp`
(reverse, the adjoint). The composition gives dJ/dp for all parameters at the
cost of one primal and one adjoint per objective, independent of the number
of parameters. Interfaces follow the standard nonlinear-programming
contract (x, f, grad f, c, Jacobian of c, optionally Hessian-vector
products), so SQP and interior-point solvers from Python (D61) plug in.
Open-source optimisers may be *used* through their public interface; the
licence check of D168 applies to what we ship.

Why not a full simultaneous analysis and design (state as optimisation
variable)? For 1-10 M cells the KKT system is too large for a general NLP
solver. We keep the nested form (solve the state, then the adjoint) and get
most of the benefit through 2.2.

Effort: wrapping the stages as tools is minutes to an hour per stage once
L17 and morphing exist; L17 itself is **not sized yet** (several PRs:
steady laminar adjoint, turbulence differentiated with D59 rather than
frozen, surface sensitivities, verification against dual-number and
finite-difference gradients).

### 2.2 One-shot and partially converged optimisation with warm starts (D206)

Converging primal and adjoint fully at every design step wastes most of the
budget. Literature: one-shot methods (Ta'asan 1991; Hazra & Schulz 2004),
piggyback iterations (Griewank & Faure 2002), inexact-gradient trust
regions (Carter 1991; Heinkenschloss & Vicente 2002).

Requirement: the optimiser may ask for a design step after a fixed number of
primal and adjoint iterations; each tool returns an error estimate of J and
of dJ/dp (residual norms times the adjoint-based error estimate, D67 when
available); the trust region accepts a step only when the predicted
reduction exceeds the estimated error, otherwise it asks for more
iterations. Every new design starts from the previous state (D52, D105).
The steady diagnosis rule of D198 (d) applies: a primal that does not settle
in 2000-3000 iterations is reported, not pushed further.

Failure mode to state plainly: on flows with a weak instability (separated
diffusers, wakes) partial convergence gives noisy J and gradients; the
error estimate then grows, the step is rejected, and the study falls back to
full convergence. That is the intended behaviour, not a bug.

### 2.3 Multi-fidelity with trust regions (D207)

Approximation and model management (Alexandrov, Dennis, Lewis & Torczon
1998): optimise a cheap model inside a trust region, corrected so that its
value and gradient match the expensive model at the centre (first-order
consistency); accept steps by the expensive model; adapt the radius.

Fidelity ladder:
1. pyaeroflow 2D section (seconds to minutes, differentiable in JAX), only
   where a 2D section is a physical proxy (airfoil-like blades, 2D fins,
   channel cross-sections);
2. Foehncraft on a coarse level of the D67 sweep;
3. Foehncraft on the production mesh;
4. final check run, grid-convergence verdict, limitations.

The tool contract carries a `fidelity` field, and results from different
levels are never mixed in one table without it (D106 already has a
fidelity column). Failure mode: if the cheap model has the wrong trend (2D
section of a strongly 3D flow), the correction is large and the trust
region shrinks to nothing; the study then reports "cheap model not useful"
and continues on level 2.

### 2.4 Gradient-enhanced Bayesian optimisation (D208)

For a small design space (up to roughly 10-20 parameters), several local
optima, or noisy objectives, a Gaussian-process (kriging) model with
gradient observations (gradient-enhanced kriging, Han et al.; Wu,
Poloczek, Wilson & Frazier 2017) uses each adjoint run as d+1 data points.
Acquisition by expected improvement or knowledge gradient, with
constraints. Above roughly 20 parameters the gradient-based path of 2.1 is
preferred, because GP cost grows with (n(d+1))^3.

### 2.5 ML surrogates from our own runs, with uncertainty (D209)

Surrogates are trained only on Foehncraft (and pyaeroflow) runs of the
study, never on third-party data (D168, D201). They carry an uncertainty
(GP variance, deep ensembles, or conformal prediction intervals) and are
used for screening and active learning: the next CFD run goes where the
surrogate is uncertain or promising. A surrogate value is never a result:
any design that is reported was computed by the solver, and D165 governs a
surrogate that enters a validation claim. This adopts the intent of the old
P14 proposal. Effort: the study-layer surrogate is a Python module of a
few hundred lines on top of D62 and D106, minutes to an hour once runs
exist.

### 2.6 Agent orchestration and budgets (D210)

An optimisation study through MCP is one tool (`study.optimise`) with a
typed formulation (parameters with bounds, objective, constraints,
fidelity ladder, budget in cell-hours, stopping rules). It refuses to start
without a budget (D198 c), writes every evaluation to the journal (D197),
and ends with "finished with N limitations" (D199): e.g. not converged,
gradient error above tolerance, best design outside the validated range of
a model, final mesh not grid-converged. The agent can pause, change the
budget (Ask action kind) or switch fidelity; it cannot accept a design the
final check rejected.

## 3. Immersed boundaries for our architecture

### 3.1 Candidates

| Family | How it works | Fit to us |
|---|---|---|
| Diffuse forcing (Peskin 2002) | Body force spread by a smoothed delta | Poor: smeared wall, no wall functions, first order |
| Volume penalisation / Brinkman (Angot, Bruneau & Fabrie 1999) | Darcy drag term with a large coefficient in solid cells; a density field for topology optimisation (Borrvall & Petersson 2003) | Very easy, smoothly differentiable, good for topology search; wall is smeared over a cell, wall shear and heat flux are inaccurate, RANS wall laws do not fit; recent work in an industrial aeronautical FV code (arXiv 2404.16132) shows it works for inviscid and laminar parts, with known limits for turbulent walls |
| Sharp ghost-cell / discrete forcing (Tseng & Ferziger 2003; Mittal et al. 2008) | Values in cells just inside the body set by interpolation through the image point | Second order for laminar flow, but not conservative, needs special stencils; for RANS wall models it shows spurious oscillations of skin friction and heat flux (Capizzano 2011; Tamaki & Imamura 2017 needed a modified eddy viscosity); weak on CHT, which needs a conservative interface flux |
| Conservative cut cell (Ye et al. 1999; Meyer, Devesa, Hickel & Adams 2010; Berger & Aftosmis 2012 for wall models) | Background cells are cut by the surface; the fluid part becomes a smaller cell with a new wall face; tiny cut cells are merged with a neighbour | Best fit: cutting a polyhedron gives a polyhedron, so the result is just an ordinary D23 mesh with real boundary faces |

### 3.2 Recommendation: cut cells as a mesh operation

On a face-based polyhedral solver the cut-cell method is not a solver
feature at all, it is a **mesh operation**: intersect the background mesh
with the surface (D108 part 3 for inside/outside and signed distance,
robust predicates as in D58), keep the fluid parts as polyhedral cells,
merge cut cells below a volume fraction (say 0.2-0.3 of the neighbour) into
the neighbour across a shared face, and hand the result to the unchanged
solver. Consequences:

- exact conservation, Rhie-Chow, AMG, halos, MPI, GPU kernels and the
  reproducible mode (D93) work untouched (L13's own condition "no change to
  the solver data structures" holds);
- the immersed wall is an ordinary wall face, so D13 wall functions, D135
  thermal wall laws and the D29 boundary set apply directly;
- CHT is natural: cut both ways, the solid parts become cells of the solid
  region (D112), and the cut face is an ordinary region interface (D160);
- the same imprinting code as D58 does the intersections, so one robust
  geometry kernel serves both.

Weak points, stated honestly:

- **Boundary layers.** A cut background mesh has no prism layers; y+ of cut
  cells varies by the cut volume fraction. With wall functions (y+ 30-300)
  results are useful for ranking designs; absolute drag, pressure loss and
  heat transfer need a body-fitted mesh with layers. Mitigation: refine the
  background near the surface (D57, L14) and, as option, extrude a thin
  layer of prisms from the cut surface into the cut cells ("cut and
  extrude").
- **Moving bodies.** Cells that change from solid to fluid between steps
  need special treatment (fresh-cell problem); for optimisation the shape
  is static per evaluation, so this does not matter now (L4 later).
- **Differentiability.** The cut topology changes when a surface crosses a
  node, so the mesh is not smooth in p. The fix is to use the *continuous
  shape derivative* on the wall faces (dJ = integral of the adjoint
  sensitivity g times the normal displacement over the surface), which only
  needs the converged primal and adjoint on the current cut mesh, not
  derivatives of the cutting. Topology changes (holes appearing) are the
  domain of penalisation (3.3).

### 3.3 Penalisation as a second, narrow tool

For early *topology* search (where should material be at all) a Brinkman
density field on a fixed mesh is smoothly differentiable and the adjoint of
D104-style source terms is already in reach of D59. The result is
converted to a surface (iso-surface of the density), cut into the mesh
(3.2) and from then on the shape is optimised with the sharp method. So
penalisation would be a topology pre-stage, not a second immersed method.

### 3.4 Owner idea (b): start immersed, switch to body-fitted later

Fits 3.2 well, with an explicit hand-off rule instead of "later":

1. early rounds on the cut mesh (fast to remesh: a new shape is a new cut,
   seconds, no mesher in the loop);
2. switch to a body-fitted mesh with layers (Gmsh script of D62, or
   morphing L4) when either the trust region has shrunk below the
   difference between the two mesh types measured on the current design,
   or the objective improvement per round falls below the uncertainty
   of the cut-cell result;
3. the switch is checked once: same design on both meshes, the difference
   recorded as a limitation and used as the first-order correction of
   2.3 (the cut mesh becomes a fidelity level).

"Unless our immersed method is good enough": that is decided per objective
by step 3, not in advance. For ranking (which fin pitch is best) cut cells
with wall functions are often enough; for absolute heat transfer or drag
they are usually not.

## 4. Owner idea (a): local optimisation in a cut-out region, then imprint

Procedure: global run; cut out a region around the part to optimise;
freeze boundary values from the global field on the cut-out boundary;
optimise the shape inside; imprint the region back (D58 non-conformal
interface does exactly this joining); global run; repeat 2-3 times, as in
the D204 cycle-averaged CHT loop.

What makes it work:
- the region is cheap (often 1-10 % of the cells), so many design steps
  run there;
- it is a Schwarz-type iteration; it converges when the shape change
  affects the region boundary only weakly;
- imprinting is already decided (D58), so joining back costs nothing new.

Where it fails, honestly:
- **Elliptic flows feel the shape upstream.** Incompressible flow is
  elliptic: blockage and pressure changes travel to the cut boundary. With
  frozen Dirichlet values everywhere, the optimiser can "win" by changes
  that would be rejected by the real outer domain (e.g. less blockage that
  really raises the flow rate, or a lower pressure loss that is really a
  shifted operating point).
- **Objectives that live outside the region** (outlet uniformity,
  downstream heat transfer) cannot be measured inside it.

Remedies, in order of value:
1. Natural boundary choice: velocity profiles upstream, pressure
   downstream, and zero-gradient tangential parts, never all Dirichlet.
2. **Adjoint-corrected local objective (new idea, recommended).** The
   global adjoint at the cut boundary gives dJ_global/d(boundary values).
   Optimise J_local + lambda . (change of the boundary fluxes the region
   produces) instead of J_local alone. To first order this is the global
   objective, which removes most of the "fake win". The global adjoint is
   only recomputed at each imprint, so it is cheap.
3. Check at each imprint: the predicted gain (local) against the realised
   gain (global); if they disagree by more than the uncertainty, enlarge
   the region. Stop after 3 cycles or when the boundary-value change is
   below tolerance (as D204).

This becomes a study type "zonal optimisation" (proposed L49).

## 5. Owner idea (c): an ML surrogate from one CFD run?

Short answer: **one run cannot train a surrogate of how the objective
changes with the shape**, however many cells it has. But one run plus its
adjoint does give something valuable.

Why not: 1-10 M cells are a lot of numbers, but they describe *one point*
in design space. The cells are strongly correlated (one smooth field), so
the information about "what happens if the shape changes" is essentially
zero. A network trained on that field learns that field. Earlier
assessment (step 20): foundation models for general flows reach a zero-shot
normalised RMSE around 0.3 on industrial cases, i.e. not usable for design
decisions; they were trained on thousands of runs.

**Several shapes in one run** (a "wind tunnel with many models"): the
shapes interact through blockage and wakes, so each sees a different
inflow and the samples are not independent; if they are spaced far enough
apart not to interact, it is equivalent to N separate runs in one mesh. It
saves setup overhead (one mesh, one AMG setup, one queue slot), not
information. It can be useful for small, local features in a developed
flow (rows of fins or roughness elements in a periodic channel, separated
by enough length) as a cheap batch, but every sample counts once.

What does work from one run:
1. **Primal + adjoint = local linear surrogate in all parameters.** J and
   dJ/dp for hundreds of parameters at the cost of about two runs; with
   Hessian-vector products (forward over reverse, D59 dual numbers on the
   adjoint) a local quadratic model. This is the real "one run gives a
   model" and it feeds D208 directly.
2. **Local closures and maps**, trained on the cells of one run (wall heat
   transfer from local features, correction fields D104): data-rich per
   run, but they generalise only within the same flow class, which is
   exactly what D165 tests with its held-out cases.
3. **Reduced-order models from one transient run** (POD or dynamic mode
   decomposition of the snapshots): good for the same geometry and nearby
   operating points, they fail for shape changes and extrapolation.
4. **Field surrogates over a design family**, trained on tens to hundreds
   of runs of the study (mesh-based graph networks or simple PCA plus GP):
   these are what 2.5 means; they need the runs, but the runs can be cheap
   levels (2.3) and placed by active learning.

## 6. Proposed entries

| ID | Title | Content (short) | Owner status |
|---|---|---|---|
| D205 | Differentiable chain of typed tools | 2.1 | approved 2026-11-21 (earlier summary) |
| D206 | Inexact one-shot steps with warm starts | 2.2 | approved 2026-11-21 |
| D207 | Multi-fidelity trust region, pyaeroflow as cheap level | 2.3 | approved 2026-11-21 |
| D208 | Gradient-enhanced Bayesian optimisation | 2.4 | approved 2026-11-21 |
| D209 | Surrogates from our own runs, with uncertainty | 2.5, replaces P14 | approved 2026-11-21 |
| D210 | Optimisation as an agent tool with budget and limitations | 2.6 | approved 2026-11-21 |
| L13 (amended) | Immersed boundaries by cut cells as a mesh operation, Brinkman topology pre-stage | 3.2, 3.3 | **decided 2026-11-21** (card 1: cut cells plus penalisation pre-stage) |
| D211 | Immersed first, body-fitted hand-off rule | 3.4 | **decided 2026-11-21** (card 2: automatic hand-off rule) |
| L49 | Zonal optimisation with adjoint-corrected boundaries | 4 | **decided 2026-11-21** (card 3: adjoint-corrected) |
| D212 | One-run findings: primal + adjoint as local model; no shape surrogate from one run | 5 | **decided 2026-11-21** (card 4: record findings, adopt primal + adjoint; no multi-shape feature) |

Order of work (all after the D196-D199 agent work and stage 7): L17 steady
adjoint and mesh morphing first, then D205/D210 tools, then D206-D209 in
the study layer, then L13 cut cells, then L49. Sizes: tool wrappers and
study-layer code minutes to an hour each; L17 and L13 not sized yet.

## 7. Sources

- P. D. Sharpe, AeroSandbox, MIT SM thesis 2021.
- N. M. Alexandrov, J. E. Dennis, R. M. Lewis, V. Torczon, Struct. Optim. 15
  (1998) 16-23.
- S. Ta'asan, ICASE report 1991; S. B. Hazra, V. Schulz, SIAM J. Sci.
  Comput. 2004; A. Griewank, C. Faure, Optim. Methods Softw. 2002.
- R. G. Carter, SIAM J. Numer. Anal. 1991; M. Heinkenschloss, L. N.
  Vicente, SIAM J. Optim. 2002.
- J. Wu, M. Poloczek, A. G. Wilson, P. Frazier, NeurIPS 2017.
- C. S. Peskin, Acta Numer. 2002; R. Mittal, G. Iaccarino, Annu. Rev. Fluid
  Mech. 2005.
- P. Angot, C.-H. Bruneau, P. Fabrie, Numer. Math. 1999; T. Borrvall,
  J. Petersson, Int. J. Numer. Meth. Fluids 2003.
- Y.-H. Tseng, J. H. Ferziger, J. Comput. Phys. 2003; R. Mittal et al.,
  J. Comput. Phys. 2008.
- T. Ye, R. Mittal, H. S. Udaykumar, W. Shyy, J. Comput. Phys. 1999;
  M. Meyer, A. Devesa, S. Hickel, N. A. Adams, J. Comput. Phys. 2010;
  M. Berger, M. Aftosmis, AIAA 2012-1301 (wall models on cut cells).
- F. Capizzano, AIAA J. 2011; Y. Tamaki, T. Imamura, AIAA J. 2017/2018.
- Volume penalisation in an industrial aeronautical FV solver,
  [arXiv 2404.16132](https://arxiv.org/pdf/2404.16132).
- Conservative Cartesian cut cells for incompressible flow,
  [arXiv 2211.10698](https://arxiv.org/pdf/2211.10698).
- Second-order thermal IBM for CHT, [arXiv 1411.1004](https://arxiv.org/pdf/1411.1004).
