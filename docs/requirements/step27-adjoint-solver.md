# Step 27 – Adjoint solver (L17): staged plan

Plan of 2026-11-22, maintainer request "Plan first" for L17. **No code.**
The maintainer answered the ten questions of section 9 on 2026-11-22; the
answers are in section 9 and adopted as D220 to D224, L51 and an amended
L17 in [requirements.md](../requirements.md). Where a later section still
says "proposed", the answer of section 9 governs. The ids were first
written as D213 to D217 and L50 and renumbered because other work took
those ids the same day. Clean room (D201):
published literature only, no copyleft source, no commercial solver names.
Time estimates are agent time (code in minutes to a few hours per pull
request, long runs one to two hours, CI separately); where a size cannot
be known before a measurement the text says "not sized yet".

Every gradient-based item of [step 26](step26-design-optimization.md)
(D205 to D212, shape parametrisation, L49) depends on this solver. Next
free ids after this page: D225, L52, P19, O25, step 28 (check the table
before taking one; other pages also reserve D213 and L50).

Contents: [0](#0-where-we-stand-checked-in-the-repository) where we stand,
[1](#1-discrete-or-continuous) discrete or continuous,
[2](#2-how-the-reverse-sweep-is-built) how the reverse sweep is built,
[3](#3-the-adjoint-linear-system) the adjoint linear system,
[4](#4-stages) stages, [5](#5-verification-ladder) verification,
[6](#6-tool-and-gradient-interface) interface, [7](#7-memory-and-cost)
memory and cost, [8](#8-proposed-entries) entries,
[9](#9-open-questions-for-the-maintainer) questions,
[10](#10-sources) sources.

## 0. Where we stand (checked in the repository)

What exists (all from stage 5 part 3, `docs/verification/stage5.md`):

- **Forward derivatives from the kernel source.** `foehncraft_core::Dual<N>`
  implements `Real`; 11 face, cell and boundary kernels agree with central
  differences to 5e-7 on every input array
  (`crates/foehncraft-flow/tests/adjoint.rs`). Switches (upwind direction,
  limiter branch, `max`/`min`) are *frozen* at the linearisation state.
  Smooth `smin`/`smax` and a smooth positive part exist for the turbulence
  and energy kernels (stage 6 §2.6, D160 §12).
- **Transposed products.** `foehncraft_kernels::adjoint` applies a matrix
  transposed as a gather through a transpose map, bitwise reproducible
  (D31); CG with multigrid and BiCGStab run on Aᵀ.
- **Shape derivatives of the geometry.** Areas, centroids and volumes
  generic over `Coord`, 165 024 derivatives checked at 2.8e-7. *Not yet:*
  the derived weights (|S|²/(d·S), the vector k, interpolation and
  least-squares weights) are `f64` only.
- **Checkpoint schedule.** Binomial checkpointing (Griewank & Walther 2000)
  with the optimal recomputation count and bitwise recomputed steps
  (`crates/foehncraft-flow/src/checkpoint.rs`).
- **Primal.** The steady solvers are SIMPLEC and the coupled solver; both
  converge to the same discrete equations (`coupled.rs`: "a converged result
  is the one of SIMPLEC"). The coupled solver applies the residual operator
  matrix-free and solves for increments with flexible GMRES and a SIMPLE-type
  preconditioner (Klaij & Vuik 2013). Turbulence models: k-ω SST, k-ε,
  realizable k-ε, k-ε-ζ-f, k-ε-ζ-f₀, elliptic-blending k-ε. Energy and CHT
  (stage 7) are in.

What does **not** exist: a reverse-mode (transposed Jacobian) evaluation of
a residual, any adjoint equation, an objective with a derivative, node
coordinate sensitivities of the *derived* geometry, mesh morphing (L4), a
gradient tool in the registry. The honest summary of step 26 section 0
stands: the foundation exists, the solver is the first job.

Two consequences found while reading the code, which shape the plan:

1. **The adjoint can be built on the residual, not on the algorithm.**
   For a steady run the converged state satisfies one set of discrete
   equations R(q, x, p) = 0 (q the state, x the node coordinates, p the
   other parameters). Whether SIMPLEC or the coupled solver got there does
   not matter, so the adjoint needs a pure, callable residual
   R(q, x, p) and no reversal of any iteration. The first work item is to
   expose that residual as a function (today it exists as parts of the
   algorithms). The state must be checked in the spike (stage A0): the
   Rhie–Chow face flux depends on the momentum diagonal, which depends on
   the face flux through convection, so the face fluxes are probably part
   of q (extra unknowns ṁ − F(u, p, ṁ) = 0), not a function of (u, p).
2. **The reverse of a gather is a scatter.** D31 forbids floating-point
   atomics. A face kernel's adjoint adds into its owner *and* neighbour
   cell, and a stencil kernel (gradient, limiter) into all cells of its
   stencil. The existing transposed row product is the pattern: express
   every reverse kernel as a gather over the inverse stencil (cell → faces
   that read it) with a fixed order. This is the same problem as "adjoint
   stencil loops" (Hückelheim et al. 2019) and the reason a naive reverse
   mode with atomics (as common in differentiable-simulation languages) is
   not acceptable for D31 and D93.

## 1. Discrete or continuous

| | Discrete adjoint (differentiate the discrete residual) | Continuous adjoint (derive PDE, then discretise) |
|---|---|---|
| Consistency | The gradient is the exact gradient of the discrete J that the solver computes; checkable to round-off against dual numbers | Gradient differs from that of the discrete J by a discretisation error that does not vanish when the primal converges; can mislead an optimiser near the optimum (Nadarajah & Jameson 2000; Peter & Dwight 2010 survey) |
| New physics | Free by D59: a new kernel generic over `Real` is differentiated by instantiation | Every equation, boundary condition, wall function and model needs its adjoint derived, discretised and verified by hand; every D29 boundary kind, D13 wall treatment, MRF, sliding interface, CHT interface doubles |
| Turbulence | Differentiated exactly with smooth limits, or frozen (a labelled approximation) | Usually frozen ("frozen eddy viscosity"), or derived by hand per model (Zymaris et al. 2009; Papoutsis-Kiachagias & Giannakoglou 2016) |
| Memory | Needs the Jacobian action: stored local Jacobians or recomputation (section 7) | Smaller: a PDE of the same size as the primal |
| Surface sensitivity | Node sensitivities dJ/dx; wall-face form by projecting on the normal | Natural: the Hadamard form (a function on the wall) |
| Topology change under shape change (cut cells, L13) | Not differentiable through the cutting | Boundary-integral form needs only the converged fields |

**Recommendation (D220 proposed): discrete adjoint, built on the kernels of
D59.** It is the form the requirements already prepare, it keeps the
kernel source as the single definition of the physics, and its gradient
is checkable against the forward derivative to round-off (D59's test). The
well-known weakness of discrete adjoints is approximation, not principle:
frozen limiters and frozen turbulence give gradient errors that can be
large (Dwight & Brezillon 2006), so the default differentiates the smooth
variants exactly and any frozen option is a *named* choice that lands in
the `limitations` list (D199).

The continuous idea is kept where it is the right tool: for **cut cells**
(L13) the shape derivative has an interior term and a boundary term for the
moving discontinuity, and the cut topology is not smooth in p. The plan
therefore also produces a *wall-face sensitivity* (the discrete adjoint
fields evaluated in Hadamard form, as in the hybrid discrete-continuous
approaches of the literature) and checks it against the exact node
sensitivities on body-fitted meshes where both exist. This is
exactly the step 26 section 3.2 hand-off (continuous shape derivative on
the cut wall faces), made verifiable.

### 1.1 Computer-graphics and differentiable-simulation benchmarks

Owner rule (2026-11-21): for geometry and interface topics benchmark
computer-graphics algorithms too. What transfers, and what does not:

| Idea | Source | Verdict for L17 |
|---|---|---|
| Adjoint as a first-class tool of graphics control: optimise smoke/fluid keyframes by the adjoint of a (discrete) Navier–Stokes solver | Treuille et al. 2003; McNamara, Treuille, Popović & Stam 2004 | Confirms the discrete route: they differentiate the discrete step, which is the same decision. Little new numerics (low-order grid solvers, no turbulence, no wall treatment) |
| **Preconditioned (Sobolev) shape gradients**: raw vertex gradients are rough and ill-conditioned; reparametrise with (I + λL) so a step moves smooth regions, which allows large steps and keeps mesh quality | Nicolet, Jacobson & Jakob 2021 ("Large steps in inverse rendering of geometry"); in CFD the same idea is Sobolev smoothing of shape gradients (Jameson & Vassberg 2001; Schmidt et al. 2013; Stück & Rung 2011) | **Adopt** as the default smoothing of wall sensitivities (stage A3); one operator, two literatures, cheap (one Laplacian solve on the surface mesh) |
| **Boundary term of a differentiated integral with a moving discontinuity**: interior derivative plus an edge/boundary integral; sampled by edge sampling or removed by reparametrisation of the integrand | Li, Aittala, Durand & Lehtinen 2018; Loubet, Holzschuch & Jakob 2019; Zhang, Miller, Yan, Gkioulekas & Zhao 2020 (path-space) | Same mathematics as the cut-cell topology problem (section above). Informs the *wall-face form* and tells us which part is missing if a cut changes between evaluations. Research reading for L13, not for the first adjoint |
| **Differentiate the fixed point, not the iterations** (implicit function theorem): one linear solve with the Jacobian transpose instead of unrolling | Christianson 1994; Griewank & Walther 2008; in learning Bai, Kolter & Koltun 2019, Blondel et al. 2022 | **Adopt**: it is the steady adjoint. Memory independent of the number of primal iterations |
| **Constant-memory reverse by replaying the forward pass** (replay instead of taping) | Nimier-David et al. 2020 (radiative backpropagation); Vicini, Speierer & Jakob 2021 (path replay backpropagation) | Same trade as our binomial checkpointing and our "recompute local Jacobians" mode (section 7): memory for time, with a fixed-order replay that stays deterministic. Confirms the design, no new method |
| **Source-transformation reverse for parallel kernels with atomics** | Hu et al. 2020 (DiffTaichi) | Atomic-based reverse conflicts with D31; use the gather-form reverse (section 2), keep as a speed comparison only |
| **Compiler-level reverse mode** for GPU kernels (LLVM) | Moses & Churavy 2020 (Enzyme) | **Candidate for a spike** (stage A0): permissive licence; in Rust still experimental, state to be checked before any reliance (section 2, option B) |
| Differentiable grid-free Monte Carlo (walk on spheres) | Miller, Sawhney, Crane & Gkioulekas 2024 ("Differential walk on spheres", arXiv 2405.12964) | Relevant only to L15 (radiation) where derivatives of a Monte Carlo estimate are needed; out of scope here, noted as L51-adjacent |
| Differentiable cloth / rigid contact with implicit differentiation | Liang et al. 2019; Qiao et al. 2020 | Same implicit-function idea as above; contact handling has no counterpart here |
| Differentiable CFD in array languages (automatic differentiation of a whole solver, unrolled) | Kochkov et al. 2021; Bezgin, Buhendwa & Adams 2023; Holl, Koltun & Thuerey 2020 | Good for 2D research and for the cheap fidelity level of D207 (PyAeroFlow in JAX). Unrolling a converged steady solve costs memory proportional to the iterations: not the route for 1–10 M cell 3D runs |

## 2. How the reverse sweep is built

The adjoint needs products with the transposed Jacobian: for a state
residual R(q; x, p), v̄ ↦ (∂R/∂q)ᵀ v̄, (∂R/∂x)ᵀ v̄, (∂R/∂p)ᵀ v̄, and for
the objective (∂J/∂q)ᵀ, (∂J/∂x)ᵀ, (∂J/∂p)ᵀ. Four ways to get them:

| Option | How | For | Against |
|---|---|---|---|
| **A. Local Jacobians from dual numbers, chained by hand** (proposed first) | The residual is a composition of a few kernel stages: gradient (linear in q, weights from the mesh) → limiter (frozen or smooth) → face flux kernel → cell sum. Each face/cell kernel is evaluated on `Dual<N>` over *its own* inputs (owner and neighbour state, their gradients, geometry: N of order 30) to get the small local Jacobian; the reverse sweep multiplies the adjoint by that transposed block and gathers into cells through the inverse stencil | Uses D59 as designed; no new AD tool; kernels stay `no_std` and GPU capable; deterministic gather form; exact, same source | Cost per product is about N times one face-kernel evaluation (forward mode over many local inputs), a factor 5–10 above an ideal reverse mode; chaining stages by hand is a structure that must be kept in step with the primal (a test enforces it, section 5) |
| **B. Compiler reverse mode (Enzyme-type) on the generic kernels** | Reverse-mode synthesis on the optimised LLVM IR of a kernel | Ideal reverse cost; works on the same source; supports GPU kernels | Rust support is experimental and nightly-only as far as known; toolchain pin (1.89 stable today); atomics in the generated scatter must be avoided for D31 |
| **C. Operator-overloading tape** (CoDiPack-style, Sagebaum, Albring & Gauger 2017) | A reverse-mode `Real` type that records operations | Mature technique in large codes | A tape per thread is heap-based; incompatible with `no_std` GPU kernels; memory per product large; hidden state against D59 (1) |
| **D. Hand-derived adjoint kernels** | Write the transposed kernels by hand | Fastest | Violates D59 (2) "no hand-derived or copied code"; breaks with every kernel change; rejected |

**Proposed (D221): start with option A** because it works with what exists,
and run a time-boxed spike on B (stage A0) to decide whether to switch
the *kernel-local* reverse step to a compiler tool later. The surrounding
structure (inverse-stencil gathers, objective and parameter plumbing, the
Krylov solve) is the same either way, so the choice can be changed
without redoing the rest.

**State and parameters.** q = (u, p, ṁ, turbulence and energy variables as
the case needs); x are the node coordinates through the *derived* geometry
(so the missing generic weights of section 0 must be made generic first);
p are boundary values, properties (D99), D104 correction fields and mesh
handles. D59 (5) already made them explicit kernel parameters.

**Frozen versus smooth.** Default: the smooth variants of D59 (3)
(`smin`/`smax`, smooth positive part) in the linearisation, with the
primal switched to `turbulence.smooth_limits = true` when the gradient
must match the primal exactly. Upwind switches and limiters have a
frozen variant only; a limiter that switches between two evaluations
makes J non-smooth in p, and that is stated, not hidden
(`limitations`: "limiter frozen; gradient is a one-sided derivative").

## 3. The adjoint linear system

Steady: solve (∂R/∂q)ᵀ λ = −(∂J/∂q)ᵀ, then
dJ/dp = ∂J/∂p + λᵀ ∂R/∂p (and the same with x).

- **Method.** Flexible GMRES (Saad 1993) on the *exact* transposed
  Jacobian product (section 2), right-preconditioned by the transpose of
  the *approximate* Jacobian that the coupled solver already uses
  (a SIMPLE-type sweep with BiCGStab per momentum component and CG with
  the pressure multigrid). This is the approximate-Newton–Krylov pattern of
  the adjoint literature (Hicken & Zingg 2010-type preconditioning of the
  exact Jacobian by a lower-order one; Lyu, Kenway & Martins 2015 for
  RANS). The transpose of the sweep reverses the order of its three
  parts and uses the transposed momentum and pressure matrices, which
  exist (D59 (4)); Ã⁻ᵀ is BiCGStab on Aᵀ, tested already.
- **Why not reverse the primal iteration (reverse accumulation of the
  fixed-point iteration, Christianson 1994; used for the unsteady and
  coupled problems of large open codes by Albring et al.).** The segregated
  SIMPLEC iteration contains inner linear solves and under-relaxation;
  reversing it means differentiating through all of them. The residual
  form needs none of that and does not depend on which primal algorithm
  converged. The fixed-point form stays an option for coupled
  multiphysics where only the partitioned iteration exists (stage A5).
- **Convergence and honesty.** Adjoint solves converge worse than the
  primal for separated and RANS flows (known; Lyu et al. 2015), and a
  steady adjoint of a flow that does not settle (shedding) is meaningless
  (Krakos et al. 2012 on limit-cycle flows). Rules: the adjoint runs only
  from a primal accepted by the D198 (d) diagnosis (2000 to 3000 iterations,
  residual trend); the result carries the adjoint residual reduction and
  the dual-weighted estimate of the gradient error (D206); an unconverged
  adjoint is a limitation, never a number without a note.
- **Hessian-vector products** (second-order optimisation, D206 trust
  regions with curvature) need forward-over-reverse. `Dual` holds an
  `f64` value, so a nested dual is not available: either a generic-value
  `Dual<T, N>` or a second tangent direction in the reverse kernels. Not
  needed for the first stages; **not sized yet**.

## 4. Stages

Order follows the maintainer's wish "steady laminar first, then
turbulence, CHT, transient". Each stage ends with: verification (section
5) green, a gallery example with at least one image (every gallery
example, D18 and the 2026-11-21 rule; a wall-sensitivity map in colour is
the image), case files, a recipe, and the grade in `features.toml` (D199).
Stages are pull requests or short series of them; sizes are agent time.

### A0. Spike and design (no solver change)

- Expose the primal residual R(q, x, p) as a callable (SIMPLEC and coupled
  share it); decide the state vector (face fluxes in or out).
- Measure option A on one face kernel: cost per transposed product,
  local-Jacobian size, inverse-stencil gather bitwise reproducibility
  (1 and 4 threads).
- Time-boxed check of option B (compiler reverse mode) on the same kernel:
  builds on our toolchain or not, cost, determinism.
- Result: a short measurement note and the answer to question Q2.
- Size: code minutes to under two hours, runs minutes. **Not sized
  further** until measured.

### A0b. Solver generic over `Real` (decision Q7)

Make the flow solver (not only the kernels) generic over `Real` so a whole
run can be complex-stepped; details and risks in
[section 9.1](#91-consequence-of-q7-a0b-a-solver-generic-over-real). The
`f64` path stays the production path, bitwise unchanged on the catalogue.
**Not sized yet**; first slice measured in A0.

### A1. Steady laminar incompressible adjoint, parameters without shape

- Reverse sweep (option A), FGMRES with transposed preconditioner,
  objectives: integrals over boundaries and regions from the D43 monitors
  (pressure or total-pressure drop between patches, force and moment on a
  patch, mass-flow weighted quantities).
- Parameters: boundary values (inlet speed, outlet pressure), properties
  (D99), D104 correction fields. No mesh sensitivity yet. This is already
  useful: field inversion (L37, D167 complement), sensitivity of a design to
  an operating point, and the first gradient for the optimiser of D205.
- Cases: 2D cavity-type, laminar channel with step, a laminar bend (the
  demonstration, see Q5).
- Size: **not sized yet** (several pull requests: reverse kernels,
  inverse-stencil gathers, objectives, Krylov/preconditioner, tool,
  verification).

### A2. Complete geometry chain and mesh morphing (L4 part)

- Make the derived geometry generic (weights, k, interpolation, least
  squares) so ∂R/∂x is exact from node coordinates.
- Mesh morphing with an *adjoint-capable* operator: free-form deformation
  box (Sederberg & Parry 1986: explicit, trivially differentiable), and one
  smoother for the volume mesh: radial-basis interpolation (de Boer, van der
  Schoot & Bijl 2007) or an elasticity / Laplacian-type solve with a
  Jacobian-based stiffening (Stein, Tezduyar & Benney 2003): its adjoint is
  one transposed solve. Choice by the spike (cost, quality, parallel
  scalability).
- Parameters p enter only through the geometry map x(p): the flow adjoint
  does not depend on the number of parameters. Geometry derivatives
  ∂x/∂p come from the FFD/RBF handle (analytic), or by finite differences of
  the *surface mesh generator* for parametric CAD (cheap: no flow solve),
  which is the honest answer to non-differentiable CAD kernels.
- Output: node sensitivities and a wall-face sensitivity field with
  Sobolev (Laplacian-preconditioned) smoothing, written as VTKHDF plus CSV;
  the coloured wall map is the gallery image (3 or 5 blade / non-symbolic
  geometries only, per the image rules).
- Size: **not sized yet**.

### A3. Turbulence

Order by use: k-ω SST first, then the k-ε family, then k-ε-ζ-f and
elliptic-blending. For each: kernels already generic with smooth limits;
the wall treatments (Kader, Spalding with a fixed Newton count, D13) are
differentiated by unrolling the fixed count (iterative-process
differentiation; Gilbert 1992); the wall distance enters the geometry chain
(derivative of the nearest-face distance, piecewise smooth, so the wall
distance of D108 becomes a differentiated kernel); `frozen_turbulence` is
an *option* with a limitation entry, default exact.
Acceptance: gradient of a wall-resolved flat plate drag and an
asymmetric diffuser pressure recovery against full-solve central
differences, frozen against exact reported side by side (gradient error of
the frozen variant as a recorded result, in the spirit of Dwight &
Brezillon 2006). Known hard case: separated flow, where the adjoint can
fail to converge; the plan reports it rather than hiding it.
Size: **not sized yet**; the first model sets the cost of the next ones.

### A4. Energy and conjugate heat transfer

Kernels are generic and carry the Jᵀ identity test (D160, V7.22). One
monolithic adjoint over the multi-region state (D112 equation domains), with
the series-conductance interface as one formula (no branch). Objectives:
wall heat flux, maximum or mean solid temperature (smooth p-norm of the
maximum), pressure drop plus heat transfer constraints: this is the heat
sink study of the demo thread. Cycle-averaged coupling (D204) is a
partitioned multi-rate loop: its adjoint is the reverse accumulation of the
loop (Gauss–Seidel style adjoint of coupled solvers; see
[arXiv 1811.00068](https://arxiv.org/pdf/1811.00068) for the technique in
multiphysics) and is a *later* item. Buoyancy (7b1) follows the stage-7
decision: later. Size: **not sized yet**.

### A5. Transient adjoint with checkpointing

The checkpoint schedule and bitwise-recomputed steps exist. Needed: the
per-step adjoint of the BDF steps and of the time-dependent boundary
inputs; geometry conservation law for moving meshes (with L4); a
dynamic-checkpoint variant for adaptive time steps (Wang, Moin &
Iaccarino 2009); optional disk and compressed snapshots for large
meshes. Scope rule (proposed Q8): **non-chaotic** windows (laminar
unsteady, URANS with a periodic or decaying state, a time horizon
shorter than the Lyapunov time). For chaotic and LES flows the plain
adjoint grows exponentially (Lea, Allen & Haine 2000; Wang 2013) and is
wrong; least-squares-shadowing methods (Wang, Blonigan & Wang 2014; Ni &
Wang 2017) are the research item **L51**, not part of L17. L38 (learned
unsteady corrections trained through unrolled steps, short horizons)
uses this stage. Size: **not sized yet**.

### A6. Extras on the same adjoint

- Adjoint-weighted error indicators for AMR (L14) and for a posteriori
  estimates of functionals (Becker & Rannacher 2001; Giles & Süli 2002;
  Venditti & Darmofal 2002); the D206 gradient error estimate uses the
  same dual-weighted residual.
- The zonal optimisation of L49 needs the *global* adjoint at the cut
  boundary: the adjoint of the boundary values is a by-product of A1 and
  A2 (sensitivity to boundary data on the cut), so no extra solver.
- Surface-form (Hadamard) sensitivity on cut-cell walls for L13.

### Order against the rest of the roadmap

A0 needs no solver change and runs beside the stage 7 work. Decision Q10:
A1 starts right after the current stage 7 part, A2 and A3 before the
later stage 7 items if the heat sink study needs gradients; A0b sits
between A0 and the verification of A1. Tool wrappers of
D205/D210 need A1 to have something to wrap.

## 5. Verification ladder

Gradients are the one output that is easy to get silently wrong, so
verification is a stage gate, not an afterthought. Each rung runs in CI
on small meshes and is a feature row in `features.toml` (D199).

1. **Dot-product test** (⟨(∂R/∂q) v, w⟩ = ⟨v, (∂R/∂q)ᵀ w⟩ to 1e-14) on
   the full residual, per stage and boundary kind; bitwise equal on 1 and
   4 threads (D31).
2. **Local-Jacobian test**: every kernel of the residual on `Dual<N>`
   against central differences (exists for 11 kernels, extended to every
   new kernel, as D59 asks).
3. **Independent derivative type**: the same kernels on a **complex**
   `Real` (comparisons on the real part; complex-step derivative,
   Squire & Trapp 1998; Martins, Sturdza & Alonso 2003) at the residual
   level: derivative to round-off with h = 1e-30, no subtraction error.
   This also tests the `Dual` implementation itself, since the two types
   share no code. A complex-step run of the whole *solver* needs a solver
   generic over `Real`: decided yes (Q7), work item A0b (section 9.1).
4. **Tangent-linear solve**: the same linear system solved forward
   (dual direction) and in reverse; dJ/dp from both must agree to the
   linear tolerance. Cheap, exact, catches errors in the reverse
   plumbing independent of FD noise.
5. **Full-solve finite differences**: central differences of J through
   two fully converged primal runs, with a step sweep and Richardson
   extrapolation, primal converged to 1e-12 on small meshes (the rule
   against long steady runs, D198 (d), applies: small meshes).
6. **Manufactured solutions** (MMS) for the *primal* already exist
   (stage 6 `mms.rs`); for the adjoint, a manufactured source so a chosen
   field solves the primal exactly gives a state at which dJ/dp of an
   integral objective has a known value on a mesh sequence: observed order
   of the *gradient* from D67 (Celik et al.). Where an exact value is not
   available, the rung is skipped and the feature grade says so.
7. **Shape derivative**: dJ/dp against central differences of re-meshed
   or re-morphed runs (same topology), with the step sweep; plus the
   wall-face (Hadamard) form against the node form on a body-fitted mesh.
8. **Gradient tool self-check** (`gradient.check`, section 6): the same
   comparisons on any user case, returning the largest relative error and
   a verdict; used by agents before an optimisation (D210), and as a
   limitation if skipped.

Grades follow D199: a rung that cannot be run is "partial", not "verified".

## 6. Tool and gradient interface

All entry points are D196 tools: typed input (D77 schema), typed output
(D80 envelope with `limitations`, D199), journalled (D197), `effects`
declared (D198 (b)), no interactive paths (D86, MCP rule of 2026-11-20).
Python calls the same ids (D61); the NLP contract of D205 is
`fun`/`jac` callables from `foehncraft.study`.

Sketch (not final):

```toml
[adjoint]
enabled     = true
method      = "discrete"          # D220
objective   = "dp_total"          # a D43 monitor or integral name, or an expression (D53)
frozen      = []                  # e.g. ["turbulence"], recorded as a limitation
tolerance   = 1e-8                # adjoint residual reduction
max_iterations = 200

[[adjoint.parameter]]
name = "inlet_speed"
kind = "boundary-value"           # | "property" | "correction-field" | "shape"
target = "inlet"

[[adjoint.parameter]]
name = "bump"
kind = "shape"                    # handles of a free-form box (A2)
handles = "examples/cases/bend/ffd.toml"
```

| Tool id (proposed) | Input | Output |
|---|---|---|
| `adjoint.solve` | case, objective(s), parameter list, tolerances | J, gradient per parameter (names, units), adjoint residual history, error estimates of J and dJ/dp (D206), `limitations`, `fidelity`, time and memory used |
| `adjoint.sensitivity` | case, objective, which boundary | wall-face sensitivity field (VTKHDF) and CSV, smoothed and raw, with the smoothing radius recorded |
| `adjoint.jvp`, `adjoint.vjp` | state or run id, direction or cotangent, parameter set | derivative of outputs along a direction (forward) or cotangent pulled back (reverse): the `value`/`jvp`/`vjp` of D205 |
| `gradient.check` | case, parameter, step list | the comparisons of rungs 3, 4, 5 (those cheap enough), maximum relative error, verdict |
| `adjoint.estimate` | case, objective | memory and time estimate before the run (D81), refuses past `[run.limits]` (D83) |

Properties required of the contract: structured input and output only;
the run record (D69) keeps objective, parameter definitions, tolerances,
the primal run it differentiates, and every approximation switch;
reproducible bitwise in the D93 mode (the reverse sweep is
gather-based with fixed reduction order); `limitations` always filled
(frozen limiter, frozen turbulence, adjoint not converged, primal not
settled, morphing quality, cut-topology change since the last
linearisation).

## 7. Memory and cost

Quantities below are *planning relations*, not promises. Every number must
be measured in A0/A1 and then replaced here.

Symbols: n cells, f faces, v variables (4 for laminar incompressible
flow, 6 with SST, more with energy/ζ-f), s snapshots, m Krylov vectors.

| Item | Relation | Remark |
|---|---|---|
| Primal state | n·v doubles plus face fluxes f doubles | already in the run |
| Adjoint vectors | λ, right-hand side, residual: 3 n v doubles | small |
| Krylov basis (FGMRES, restart m) | 2 m n v doubles (basis and preconditioned vectors) | m = 20 to 40; the dominant vector memory, tunable |
| Preconditioner | as the coupled solver's own (matrices, AMG hierarchy) | no extra beyond transposed maps (a gather map of 4 bytes per slot) |
| Local Jacobians, **stored** mode | per face a (v_out × N) block, 4 × 30 doubles ≈ 1 kB for laminar flow, ≈ 3·f kB per case | fast, small meshes (rule of thumb below ≈ 1 M cells) |
| Local Jacobians, **recompute** mode | no storage; cost per product ≈ N face-kernel evaluations | default for large meshes; the same "replay for memory" trade as checkpointing (Vicini et al. 2021) |
| Steady adjoint time | about one to a few times a primal solve for laminar flow; clearly more for RANS and separation (Lyu et al. 2015) | **to be measured**; D198 (d) caps the budget, the study layer sees an estimate before the run |
| Cost per optimisation iteration | one primal plus one adjoint, **independent of the number of parameters** | the reason for the adjoint (Giles & Pierce 2000) |
| Transient | with s snapshots and up to t recomputations of a step, t·n_steps − β(s+2, t−1) primal steps plus one adjoint step per step (existing schedule); at s = 20, 10 000 steps need about 3.8 primal steps per step on average | a snapshot is the restart content (state, old levels, fluxes), about 12–20 doubles per cell for laminar to RANS flows: 1–2 GB per 10 M cell snapshot; disk and compressed snapshots as options |
| Forward (tangent) alternative | one primal-sized solve **per parameter** | used for rung 4 and when parameters are fewer than objectives |

Cost sanity check for the agent budget (D198 c): an optimisation with
n_it design iterations costs n_it × (primal + adjoint) cell-hours; the
tool reports this before it starts.

## 8. Proposed entries

Adopted 2026-11-22 after the answers of section 9, with the changes noted
there. Next free ids: D225, L52.

| ID | Title | Content |
|---|---|---|
| D220 | Discrete adjoint on the D59 kernels | Section 1: discrete adjoint of the converged residual; exact (smooth) by default; frozen limiter or turbulence only as named options with a limitation; independent of the primal algorithm |
| D221 | Reverse sweep from local dual Jacobians, scatter-free | Section 2 option A first, compiler reverse mode as a spike (A0); every reverse kernel is a gather over the inverse stencil, bitwise reproducible (D31, D93) |
| D222 | Gradient verification ladder | Section 5: dot-product, local Jacobian, complex-step (residual level and, after A0b, the whole solve), tangent-linear, full-solve differences, shape derivative; a stage is "done" only with its rungs green or graded partial |
| D223 | Adjoint tools and gradient contract | Section 6: `adjoint.solve`, `adjoint.sensitivity`, `adjoint.jvp`/`vjp`, `gradient.check`, `adjoint.estimate`; fills `limitations`, error estimates, effects and budget |
| D224 | Scope rule for the unsteady adjoint | Section 4 A5: non-chaotic windows only; chaotic and LES sensitivities are a separate research item |
| L51 | Sensitivities of chaotic flows | Least-squares shadowing (Wang et al.) for LES/DES; open research; also differentiable grid-free radiation (differential walk on spheres) once L15 exists |
| L17 (amended) | Staged as A0–A6 (section 4) | Replace the one-line L17 text by the staged plan; no change of intent |

## 9. Decisions of the maintainer (2026-11-22)

Asked as ten cards, one at a time, each with a recommendation. Coding still
starts only on the maintainer's go (AGENTS.md).

| # | Question | Answer | Consequence |
|---|---|---|---|
| Q1 | Method | **Discrete adjoint** (recommended) | D220 |
| Q2 | Compiler reverse-mode tool | **Time-boxed spike only** (recommended) | D221; spike in A0, a nightly toolchain is allowed for the spike only |
| Q3 | Frozen turbulence and limiters | **Named options, exact by default** (recommended) | Section 2 and A3; always a `limitations` entry; gradient error of the frozen variant is measured |
| Q4 | Primal algorithm | **Residual function** (recommended) | SIMPLEC and the coupled solver share one callable residual; the coupled solver's preconditioner is reused transposed; its speed work stays closed |
| Q5 | First demonstration | **Laminar bend pressure drop** (recommended) | A1 demo with a coloured wall-sensitivity map; the heat sink study after A4 |
| Q6 | Shape sensitivity form | **Exact node form first** (recommended) | A2: node sensitivities with morphing first, wall-face form checked against them |
| Q7 | Whole-solver complex step | **Yes, make the solver generic** (*not* the recommendation) | New work item A0b below; rung 3 of section 5 extends to full-solve complex-step |
| Q8 | Transient scope | **Non-chaotic windows only** (recommended) | D224, research item L51 |
| Q9 | GPU | **CPU first, GPU after A1** (recommended) | Kernel structure GPU-capable from the start; bitwise parity after A1 |
| Q10 | Priority against stage 7 | **Spike now, A1 soon** (recommended) | A0 beside stage 7, A1 right after the current stage 7 part; A2 and A3 before the later stage 7 items if the heat sink study needs gradients |

### 9.1 Consequence of Q7: A0b, a solver generic over `Real`

The maintainer wants a whole run to be complex-steppable, so the flow
solver (orchestration, linear solvers, turbulence and energy hosts) must
become generic over `Real`, not only its kernels. Findings and plan, to be
confirmed in A0:

- The kernels are already generic (D59). What is `f64`-only: the derived
  geometry weights (also needed for shape derivatives, A2), the solver
  orchestration (`simplec.rs`, `coupled.rs`, `unsteady.rs`, `energy.rs`,
  `cht.rs`), the Krylov and multigrid code, and the device layer.
- A complex `Real` has to run through all of them, including GMRES with
  complex arithmetic (a conjugation-free, "analytic" inner product, as the
  complex-step technique requires), so every dot product and norm must use
  the plain (non-conjugated) bilinear form and the real part for
  convergence decisions. Comparisons and `max`/`min` act on the real part,
  like `Dual`.
- Costs: compile time and code size grow with the instantiations;
  speed of the `f64` path must not change. Rule: the `f64` path stays the
  production path and is bitwise unchanged (the 20 flow cases of the
  catalogue compare with `tools/check/compare_runs.py`, as in stage 5 part
  3). The generic code is monomorphised; the GPU path stays `f64`.
- The pay-off: one run with h = 1e-30 gives the gradient of the full
  solved J to round-off, independent of finite-difference noise: the
  strongest check of the adjoint, and also of the primal's own
  differentiability. It also makes the tangent-linear solve (rung 4) a
  by-product.
- Order: A0b follows the spike and the residual function (Q4) and comes
  before A1's verification, because A1's acceptance includes the full-solve
  complex-step comparison. **Not sized yet**: it touches the largest files of
  the flow crate; the spike measures a first slice (one steady laminar
  case, SIMPLEC only) before the full effort is estimated. Risk to state
  plainly: this refactor is larger than the adjoint's reverse sweep and
  could dominate A1; if the spike shows that, the question returns to the
  maintainer with the measured size.

## 10. Sources

Checked against a web search result on 2026-11-22 (title and existence;
bibliographic details of the others are from memory and must be verified
when the entry is adopted):

- B. Miller, R. Sawhney, K. Crane, I. Gkioulekas, *Differential walk on
  spheres*, 2024, [arXiv 2405.12964](https://arxiv.org/abs/2405.12964).
- B. Nicolet, A. Jacobson, W. Jakob, *Large steps in inverse rendering of
  geometry*, ACM TOG (SIGGRAPH Asia) 2021,
  [project page](https://rgl.epfl.ch/publications/Nicolet2021Large).
- J. Sagebaum, T. Albring, N. R. Gauger, *High-performance derivative
  computations using CoDiPack*, 2017,
  [arXiv 1709.07229](https://arxiv.org/pdf/1709.07229).
- *Accurate gradient computations for shape optimization via discrete
  adjoints in CFD-related multiphysics problems*,
  [arXiv 1811.00068](https://arxiv.org/pdf/1811.00068).
- A. McNamara, A. Treuille, Z. Popović, J. Stam, *Fluid control using the
  adjoint method*, ACM TOG (SIGGRAPH) 2004,
  [PDF](https://grail.cs.washington.edu/projects/control/fluidAdjoint.pdf).
- M. Nimier-David, S. Speierer, B. Ruiz, W. Jakob, *Radiative
  backpropagation*, ACM TOG 2020,
  [project page](https://rgl.epfl.ch/publications/NimierDavid2020Radiative);
  D. Vicini, S. Speierer, W. Jakob, *Path replay backpropagation*, ACM TOG
  2021, [record](https://infoscience.epfl.ch/record/288010?ln=en).

From memory, to verify before adoption:

- M. B. Giles, N. A. Pierce, Flow, Turbul. Combust. 65 (2000) 393–415
  (an introduction to the adjoint approach to design); already cited in D59.
- S. K. Nadarajah, A. Jameson, AIAA 2000-0667 (continuous against
  discrete adjoint); J. Peter, R. P. Dwight, Comput. Fluids 39 (2010)
  (survey of numerical sensitivity approaches); R. P. Dwight, J. Brezillon,
  AIAA J. 44 (2006) (approximations of the discrete adjoint).
- A. Jameson, J. C. Vassberg, AIAA 2001-2530; S. Schmidt, C. Ilic,
  V. Schulz, N. R. Gauger, AIAA J. 51 (2013) (Sobolev / shape-calculus
  gradients); J. Stück, T. Rung, Comput. Fluids 47 (2011).
- Z. Zymaris et al., Comput. Fluids 38 (2009); D. I. Papoutsis-Kiachagias,
  K. C. Giannakoglou, Arch. Comput. Methods Eng. 23 (2016) (continuous
  adjoint for turbulent flows).
- Z. Lyu, G. K. W. Kenway, J. R. R. A. Martins, AIAA J. 53 (2015)
  (RANS adjoint, approximate Newton–Krylov); J. E. Hicken, D. W. Zingg,
  AIAA J. 48 (2010).
- A. Griewank, A. Walther, *Evaluating Derivatives*, SIAM 2008; A. Griewank,
  A. Walther, ACM TOMS 26 (2000) (binomial checkpointing, in D59);
  B. Christianson, Optim. Methods Softw. 3 (1994) (reverse accumulation of
  attractive fixed points); J. C. Gilbert, Optim. Methods Softw. 1 (1992).
- Q. Wang, P. Moin, G. Iaccarino, SIAM J. Sci. Comput. 31 (2009) (minimal
  repetition dynamic checkpointing).
- D. J. Lea, M. R. Allen, T. W. N. Haine, Tellus 52A (2000); Q. Wang,
  P. Blonigan, Q. Wang, J. Comput. Phys. (least-squares shadowing, 2014);
  A. Ni, Q. Wang, J. Comput. Phys. 347 (2017); J. Krakos, Q. Wang,
  S. Hall, D. Darmofal, AIAA J. 50 (2012).
- J. R. R. A. Martins, P. Sturdza, J. J. Alonso, ACM TOMS 29 (2003);
  W. Squire, G. Trapp, SIAM Rev. 40 (1998) (complex step).
- T. W. Sederberg, S. R. Parry, SIGGRAPH 1986 (free-form deformation);
  A. de Boer, M. S. van der Schoot, H. Bijl, Comput. Struct. 85 (2007)
  (RBF mesh deformation); K. Stein, T. Tezduyar, R. Benney, J. Appl. Mech.
  70 (2003).
- R. Becker, R. Rannacher, Acta Numer. 10 (2001); M. B. Giles, E. Süli,
  Acta Numer. 11 (2002); D. A. Venditti, D. L. Darmofal, J. Comput. Phys.
  176 (2002).
- J. Hückelheim et al., ICPP 2019 (automatic differentiation of adjoint
  stencil loops); Y. Hu et al., ICLR 2020 (DiffTaichi); W. S. Moses,
  V. Churavy, NeurIPS 2020 (Enzyme).
- A. Treuille, A. McNamara, Z. Popović, J. Stam, SIGGRAPH 2003; T.-M. Li,
  M. Aittala, F. Durand, J. Lehtinen, SIGGRAPH Asia 2018; G. Loubet,
  N. Holzschuch, W. Jakob, SIGGRAPH Asia 2019; Y. Liang et al., NeurIPS
  2019; Y.-L. Qiao et al., ICML 2020; S. Bai, J. Z. Kolter, V. Koltun,
  NeurIPS 2019; M. Blondel et al., NeurIPS 2022.
- D. Kochkov et al., PNAS 118 (2021); D. Bezgin, A. Buhendwa, N. Adams,
  Comput. Phys. Commun. 282 (2023); P. Holl, V. Koltun, N. Thuerey,
  NeurIPS 2020.
- Klaij & Vuik 2013 and Saad 1993 are cited as in `coupled.rs`.
