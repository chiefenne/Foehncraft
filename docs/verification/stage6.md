# Stage 6 – RANS turbulence and wall treatments

Stage 6 of the [stage plan](../requirements.md#7-stage-plan) brings the
RANS turbulence models and their wall treatments
([specification](../requirements/stage6-turbulence.md)). This page grows
with the pull requests of the stage. Part 6a gives every later model its
transport equation, its wall distance and its solver path
([below](#passive-scalars-and-the-wall-distance-part-6a)), part 6p the
periodic boundaries that the channel cases need
([below](#periodic-boundaries-part-6p)), part 6b the k-ω SST model
([below](#k-ω-sst-part-6b)), part 6c the k-ε family with its wall
functions ([below](#k-ε-family-part-6c)), part 6d the k-ε-ζ-f model and
the multigrid of the turbulence equations
([below](#k-ε-ζ-f-and-the-scalar-multigrid-part-6d)), part 6d2 the
k-ε-ζ-f₀ model of Popovac (2023) with its wall treatment
([below](#k-ε-ζ-f₀-and-the-compound-treatment-of-2023-part-6d2)), part 6e the
coupled pressure–velocity algorithm and the pseudo-time steps that cut the
iterations of the steady turbulent cases
([below](#coupled-algorithm-and-pseudo-time-steps-part-6e)), part 6f the
multigrid of its momentum sweep
([below](#multigrid-for-the-momentum-sweep-part-6f)), part 6g the
multigrid of the turbulence equations and the larger steps with it
([below](#multigrid-for-the-turbulence-equations-part-6g)), part 6h the
pressure solve of the coupled preconditioner
([below](#the-pressure-solve-of-the-coupled-preconditioner-part-6h)), and
the elliptic-blending k-ε of Billard & Laurence (2012)
([below](#elliptic-blending-k-ε-d133)); the first validation case of
part 6g, the backward-facing step against the experiment with six
models on three meshes
([below](#backward-facing-step-against-the-experiment-part-6g-95)), and
the second, the impinging jet with its heat transfer and the
stagnation-point anomaly
([below](#impinging-jet-at-hd--2-part-6g-96)). The
commands are in
[recipe 63](../recipes/scalars.md#passive-scalars),
[recipe 67](../recipes/boundary-conditions.md#periodic),
[recipe 68](../recipes/turbulence.md#sst),
[recipe 71](../recipes/turbulence.md#k-epsilon),
[recipe 73](../recipes/turbulence.md#zeta-f),
[recipe 74](../recipes/turbulence.md#zeta-f0),
[recipe 78](../recipes/flow.md#coupled),
[recipe 80](../recipes/flow.md#momentum-multigrid),
[recipe 81](../recipes/flow.md#turbulence-multigrid),
[recipe 83](../recipes/flow.md#pressure-hierarchy),
[recipe 88](../recipes/turbulence.md#elliptic-blending),
[recipe 89](../recipes/turbulence.md#step-validation) and
[recipe 90](../recipes/turbulence.md#impinging-jet), the pictures in
[gallery 63](../gallery/scalars.md#passive-scalars),
[gallery 67](../gallery/boundary-conditions.md#periodic),
[gallery 68](../gallery/turbulence.md#sst),
[gallery 71](../gallery/turbulence.md#k-epsilon),
[gallery 73](../gallery/turbulence.md#zeta-f),
[gallery 74](../gallery/turbulence.md#zeta-f0),
[gallery 78](../gallery/flow.md#coupled),
[gallery 80](../gallery/flow.md#momentum-multigrid),
[gallery 81](../gallery/flow.md#turbulence-multigrid),
[gallery 83](../gallery/flow.md#pressure-hierarchy),
[gallery 88](../gallery/turbulence.md#elliptic-blending),
[gallery 89](../gallery/turbulence.md#step-validation) and
[gallery 90](../gallery/turbulence.md#impinging-jet), the generated
tables of every case in the [verification manual](manual.md).

## Passive scalars and the wall distance (part 6a)

**The equation.** Every scalar of stage 6 and later (the passive scalars
of D51 now, k, ω, ε, ζ, f of the turbulence models, the temperature of
stage 7) is one instance of the generic transport equation of D116,

ρ ∂φ/∂t + ∇·(ṁ φ) − φ ∇·ṁ − ∇·(ρΓ ∇φ) = ρ (S_u + S_p φ),  S_p ≤ 0,

written once as single-source kernels (`foehncraft-kernels/src/scalar.rs`,
CPU loops and PTX entries `scalar_faces`, `scalar_rows`, `scalar_scale`,
gathers in the stored face order, so CPU and GPU give the same rows bit
for bit, D31). The face coefficients are those of the stage 2 convection
and diffusion operators: upwind with the linear-upwind or bounded central
deferred correction under the Venkatakrishnan limiter, over-relaxed
non-orthogonal diffusion with the least-squares gradient of the scalar.
Three parts are new:

- **Bounded form.** The term −φ_P Σ_f ṁ_f, the net mass outflow of the
  cell from the current fluxes, is implicit. The upwind operator then has
  zero row sums and stays an M-matrix while continuity is still
  converging, so a scalar between two bounds stays between them in the
  outer iterations of a steady run (Patankar 1980; the form of the
  industry standards).
- **Sources.** S_u explicit, S_p ≤ 0 implicit in the diagonal; a
  positive S_p is taken as zero, so the diagonal can only grow (the
  turbulence destruction terms use it in 6b).
- **Bounded BDF2** (D28 (1)). BDF2 is backward Euler from the old value
  φ* = −(a₁ φⁿ + a₂ φⁿ⁻¹)/a₀ with the step Δt/a₀. Where φ* leaves the
  `bounds` of the scalar, the cell takes the coefficients of BDF1 for this
  step: a select per cell, no branch beyond it, and the history counts the
  cells. Crank–Nicolson runs integrate the scalars with BDF2. At the start
  of a BDF2 step the scalar is extrapolated, φ ← φⁿ + ω (φⁿ − φⁿ⁻¹) with
  ω = Δt/Δt_prev, before its deferred correction is evaluated: with one
  outer iteration per step (PISO) the correction would otherwise lag a
  whole step and the scheme would be first order in time (observed 1.23
  before, 1.98 after).

**Boundary conditions.** The face value is φ_f = a φ_P + b per boundary
face: a fixed value (a = 0, b = φ_b) at inlets and walls, a fixed
diffusive flux (a = 1, b = q d_n/Γ, D116) at walls, slip walls, symmetry
planes and wedges, zero gradient at outlets. The tables
`[boundary.<name>.scalar.<scalar>]` are read before the flow conditions
are parsed, so the flow keys of every boundary kind stay checked as
before; `SCAL001` reports a table that does not fit its boundary.

**Exact wall distance** (D48, D108 part 1). A bounding volume hierarchy
over the triangles of the wall faces (the fans around the node average that the face geometry uses; binned
surface area heuristic with 16 bins, at most 4 triangles per leaf) answers
closest-point queries (Ericson 2005, section 5.1.5, for the point on a
triangle). The traversal visits the nearer child first and also every box
at the same distance as the best found, and a tie between faces goes to
the smallest global id, so the result does not depend on the threads, the
partition or the order of the faces (D93). The tests compare it with a
brute-force search on random surfaces and with min(y, 1 − y) in a channel
on hex, poly and tet cells (exact to 1e-12).

**Smooth limits** (spec 2.6). `foehncraft_core::math` has the hard and
smooth min, max and abs that the turbulence limiters of 6b to 6e use;
`Switches::HARD` gives the published models exactly, `Switches::SMOOTH`
the C¹ versions for the adjoint (D104).

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Laminar heated channel (`channel_heated.toml`), 512 / 2 048 / 8 192 cells | observed order of the temperature error | 1.8 – 2.2 | 2.01, 2.14 |
| | Nu on the finest mesh against 140/17 = 8.2353 | ± 0.5 % | 8.2398 (+0.055 %) |
| Manufactured solution 1 + cos(πx/2) sin(πy) in Poiseuille flow, Pe = 50 (`tests/scalars.rs`) | observed order on hex cells | 1.8 – 2.3 | 2.07, 2.15 |
| | on polyhedra | 1.6 – 2.6 | 2.50, 2.34 |
| | on triangular prisms | 1.6 – 2.6 | 2.25, 1.91 |
| Taylor–Green vortex carrying a scalar, Δt = 0.04 … 0.005 s | observed time order, BDF2 | 1.8 – 2.6 | 2.02, 1.98 |
| | observed time order, Crank–Nicolson | 1.8 – 2.6 | 2.02, 2.02 |
| Bounded front (`scalar_front_bdf2.toml`) | largest value with `bounds = [0, 1]` | ≤ 1 | 1 − 7e-11 (plain BDF2: 1.00033) |
| Restart of a transient run with a scalar after 8 of 16 steps | the scalar and its old level at the end | bit-identical | bit-identical |
| Wall distance in a channel, hex, poly, tet | against min(y, 1 − y) | 1e-12 | passed |

The heated channel is the open variant of the 9.2 case: the developed
profile at the inlet and a zero-gradient outlet, whose missing axial
diffusion changes the temperature by about Γ G / U = 2e-6 K in the last
cell. The periodic variant with the heat-flux jump of D109 follows with
the periodic boundaries (part 6p).

## Periodic boundaries (part 6p)

**The pairs.** Two boundaries of `type = "periodic"` that name each other
form a translational pair (D29). `foehncraft_mesh::periodic::pair`
matches each face of one side with the face of the other whose centroid
lies at the translated position (a hash grid of the face centroids, within
1e-6 of the face size; the translation from the area-weighted centroids or
from `translation`) and checks that the two face areas are opposite. The
paired faces move to a block between the inner and the boundary faces,
side `a` and side `b` of each pair one after the other. Each side keeps
its own face with its own geometry, owned by its cell and listed only in
that cell: the neighbour of the face is the owner of its partner, seen at
its centroid plus the shift (−T on side `a`, +T on side `b`). The same
upper-triangular face and matrix layout then serves every kernel: the
face loops take the shifted centroid for the distance, the interpolation
weight and the least-squares stencil (`fv::neighbour_centroid`), the
transposed matrix pairs the two sides through the partner. A pair carries
the flux twice, once from each side, with the same values from the same
inputs, so mass is conserved across it to round-off. The pairing is the
last mesh step, after the renumbering; the mesh needs two cells across the
period and runs on one partition for now. Generated boxes give the nodes
of the two sides the same perturbation when the case pairs them.

**The driving.** The pressure is periodic, its level held in the reference
cell; a uniform body force β t̂ along the translation stands for the mean
pressure gradient. With `pressure_gradient = G` it is fixed, β = G. With
`mass_flow` it is adjusted after every pressure correction: a force
increment δβ moves the cell velocities by δβ D^C t̂ (the response of the
SIMPLEC rows) and the fluxes by the flux of that increment, so
δβ = (ṁ_target − ṁ)/K with K the flux of a unit increment through side
`b`; after the correction the mass flow through the pair is the target to
round-off. A uniform start velocity carries the target flow from the first
iteration. β is in the log, the report (`periodic`) and the restart file,
and it is part of the body force of the kinetic-energy balance.

**The scalars** (D109 (1)). A scalar crosses a pair unchanged
(`periodic = "periodic"`, the default) or with the jump of the heat-flux
mode, φ_N seen = φ_N ± Δ (minus on side `a`, plus on side `b`), in the
face coefficients (the jump goes into the source with the implicit
neighbour coefficient), the gradients and the limiters. Δ = Q/ṁ with Q
ρ times the integrals of the wall fluxes and the source and ṁ the mass
flow through the pair of the current fluxes, so the balance of the module
closes in every iteration. In a steady run the level is fixed after each
solve: the bulk value of the plane, the mass-flux mean of the face values
reconstructed with the gradients from both sides, is shifted to `bulk`. A
transient run keeps the level of its start. Fixed values on walls or
inlets are refused with the heat-flux mode (D109 (4)).

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Periodic Poiseuille flow (`poiseuille_periodic.toml`), mass flow 0.1 kg/s, 512 / 2 048 / 8 192 cells | observed order of the velocity error | 1.8 – 2.2 | 1.98, 1.94 (open channel 1.98, 1.89) |
| | L2 error against the open channel on the same meshes | factor 2 | 2.7e-3, 6.7e-4, 1.8e-4 (open 1.6e-3, 4.1e-4, 1.1e-4) |
| | driving force against 12 μU/h² = 1.2 Pa/m | | −0.62 %, −0.06 %, +0.09 % |
| | the same driven by 1.2 Pa/m on uniform meshes: mass flow against 0.1 kg/s | order 1.8 – 2.2 | +0.78 %, +0.20 %, +0.05 % |
| Periodic heated channel (`channel_heated_periodic.toml`), heat-flux jump | observed order of the temperature error | 1.8 – 2.2 | 2.40, 2.02 |
| | Nu on the finest mesh against 140/17 = 8.2353 | ± 0.5 % | 8.2460 (+0.13 %; open channel +0.055 %) |
| | jump Δ = Q/ṁ and the bulk value at the plane | 1e-12 | 4e-3 K, −17/70 K |
| Short periodic heated pipe (`pipe_heated_periodic.toml`), one diameter, wedge of 0.5°, 64 / 256 / 1 024 cells | observed order of the velocity error | 1.8 – 2.2 | 2.02, 1.98 |
| | observed order of the temperature error | 1.8 – 2.6 | 2.30, 2.28 |
| | Nu on the finest mesh against 48/11 = 4.3636 | ± 0.5 % | 4.3691 (+0.13 %) |
| | driving force against 8 μU/R² = 3.2 Pa/m | | −1.06 %, −0.28 %, −0.09 % |
| Pairing of perturbed hex, poly and tet boxes in x and z | every face matched, cell geometry | round-off | passed |
| Transient periodic run with the heat-flux mode, restart after 6 of 12 steps | fields, scalar and β at the end | bit-identical | bit-identical |
| Every case without periodic boundaries | results | bit-identical | unchanged |

The temperature of the pipe runs without a limiter: the Venkatakrishnan
limiter clips the smooth minimum on the axis in the prisms next to it, and
on the perturbed meshes the observed order falls to 1.65 and 1.42 (on a
fourth mesh); without it the order is 2.3 and 2.2, and on unperturbed
meshes 2.0 either way. The first pair of the channel sweep converges
faster than second order (2.40); the finest pair, which the check uses,
gives 2.02.

## k-ω SST (part 6b)

**The model** (D130). SST 2003 of Menter, Kuntz & Langtry: k and ω on the
generic scalar path of 6a with the F₁ blend of the constants, the cross
diffusion (1 − F₁) 2σ_ω2/ω ∇k·∇ω split by its sign (the positive part
explicit, the negative part in the diagonal), the eddy viscosity
ν_t = a₁k/max(a₁ω, S F₂), the production P̃ = min(ν_t S², 10 β* k ω)
(or Kato–Launder, Durbin, none), the sustaining terms of Spalart & Rumsey
and a cap of ν_t/ν. The destruction is linearised point-implicitly
(β*ω in the diagonal of k, 2βω with the Newton option for ω), so both
stay positive in a converging run; the floors (10⁻¹⁰ of the largest start or
inlet k, and the ω that gives ν_t/ν = 100 with it) are counted (TURB002 when they act in the final tenth of
the solves). Without walls the wall distance is infinite and F₁ = 0
(TURB005). The gradients of k and ω use the Venkatakrishnan limiter with
a threshold relative to the range of the variable, ε = K (φ_max − φ_min),
K = 0 by default: the dimensional threshold of the momentum equations is
far above the differences of k on a flat plate and would switch the
limiter off.

**The walls** (D134). `low-re`: k = 0 and ω_w = 60ν/(β₁ d₁²) as Dirichlet
values on the wall face, d₁ = 2 y_P the height of the wall cell, the value
of the TMR codes (the first-cell value 6ν/(β₁ y_P²) of the draft converged
at first order to +1.1 %; specification 4.3). `automatic`:
u_τ = (u_vis⁴ + u_log⁴)^¼ with the log law floored at y⁺ = 11.06, the
first-cell ω_P = (ω_vis² + ω_log²)^½ imposed with the weight m = 1 − b,
b = u_vis⁴/u_τ⁴, through a fractional fix of the row,
(1 − m)·row + m a_P (φ − ω_P), and the production of the first cell
blended likewise, (1 − m) P̃ + m u_τ³/(κ y_P). Rough walls take the wall
values of Aupoix (2015). The first-cell height y_P is floored at
10⁻³ V/A of the wall cell.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Flat plate (`tmr_plate_sst.toml`), `low-re`, 768 / 3 072 / 12 288 / 49 152 cells, first-cell y⁺ 0.09 on the finest | C_f(0.970084) on the finest mesh against the Coles–Fernholz fit C_f = 2[ln(Re_θ)/0.384 + 4.127]⁻² (Nagib, Chauhan & Monkewitz 2007, a fit to measured zero-pressure-gradient layers) at the run's own Re_θ = 7 848 (θ and the edge velocity from the profile at x = 0.97) | ±8 % | 2.7005e-3 (+1.45 %) |
| | observed order and GCI_fine of C_f (three finest meshes) | p 1.5 – 2.5, GCI ≤ 1 % | p = 0.44, GCI 1.76 %, extrapolated 2.7385e-3: **finding** |
| | C_D (one side), observed order and GCI_fine | GCI ≤ 1 % | 2.8605e-3; p = 0.75, GCI 0.83 % |
| | C_f(Re_x) against White (2006), 10⁶ ≤ Re_x ≤ 10⁷ | ±8 % | −8.5 % to −5.4 % (−6.1 % at x = 0.97): within the band but for the first 0.5 point at the upstream end: **finding**, the model's difference from the correlation |
| | C_f(Re_x) against Coles–Fernholz at the Re_θ of the plate (momentum integral of C_f/2), 10⁶ ≤ Re_x ≤ 10⁷ | ±8 % | −1.4 % to +2.5 % |
| | u⁺(y⁺) at x = 0.97 against Spalding's law, 1 ≤ y⁺ ≤ 300 | ±3 % | 2.3 % largest |
| The plate with `automatic`, 3 072 cells, first-cell y⁺ 1.7 / 4.6 / 9.3 / 27 / 92 | C_f(0.97) against our own finest value 2.7005e-3 | ±2 % | +1.9 %, +15.2 %, +14.4 %, −0.7 %, −1.1 % (the runs of the first version, re-expressed): **finding** in the buffer layer, TURB001 warns there |
| Channel at Re_τ = 395 (`channel_sst.toml`), periodic module, mass flow of the DNS bulk velocity, `low-re`, 96 cells across, first-cell y⁺ 0.28 | U_b⁺ against the DNS of Moser, Kim & Mansour (17.55) | ±3 % | 17.43 (−0.70 %; Dean's correlation 17.24: +1.1 %), Re_τ = 397.8 |
| | U_c⁺ against the DNS (about 20.1) | ±4 % | 19.64 (about −2.4 %) |
| | observed order and GCI_fine of U_b⁺ on 192, 384 and 768 cells across, converged to 10⁻⁹ (U_b⁺ 17.367, 17.333, 17.317; 96 cells: 17.426) | p 1.8 – 2.2, GCI ≤ 0.5 % | p = 1.01, extrapolated 17.300, GCI_fine 0.12 %: **finding** (first order; PR #128 replaced the values of the 10⁻⁶ runs, 17.427, 17.381, 17.346, p 0.40, whose finest mesh had stopped short, see [the level-3 channel drop](#the-level-3-channel-drop-pr-128)) |
| Homogeneous decay in a periodic box, k₀ = 1, ω₀ = 10, t = 0.5, 10 / 20 / 40 steps | BDF2 time order of k and ω against the exact solution (F₁ = 0) | 1.8 – 2.2 | k 1.92, 1.96 and ω 2.07, 2.03 (two pairs of meshes) |
| Steady plate, restart after half the iterations | fields at the end | bit-identical | bit-identical |
| Plate on 1 and 4 threads | fields after 8 iterations | bit-identical | bit-identical |
| Poisson wall distance on a plane wall | d against the exact distance | round-off | passed |

**Findings.** (0) U_b⁺ of the channel approaches 17.30 (the extrapolated value; Dean's correlation 17.24) at first order (p = 1.0 on the converged runs of PR #128; 0.40 on the first, partly unconverged ones) with the first-cell y⁺ falling from 0.28 to 0.03; the plate converges at 0.43: the wall value of ω, which grows as 1/d₁² with the mesh, converges slowly in both; a sweep of the wall value is an item for 6c. (1) The C_f of the plate converges monotonically but at a
low observed order on these meshes (0.43); the finest value is within
1.5 % of the Coles–Fernholz fit, and the
differences between the three finest meshes are below 0.7 %. A fifth mesh
(512 × 384) will show whether the order rises;
the band stays as it is. (2) The automatic treatment overpredicts C_f by
up to 16 % with the first cell in the buffer layer: the blended u_τ of the
original formulation and the log-layer production meet badly there. The
check TURB001 warns for y⁺ between 3 and 20; an improved blending is an
item for the wall-function work of 6c. *(Fixed in 6g, PR #112: see
[the automatic treatment in the buffer layer](#the-sst-automatic-treatment-in-the-buffer-layer-6g-pr-112).)* (3) White's correlation is a fit
to experiments that lies 5 to 9 % above the SST result along the whole
plate, while the Coles–Fernholz fit to measured layers lies within
−1.4 to +2.5 % of it: the White band of the specification was set for the
model and was widened to ±8 % (2026-11-22), which SST meets except at the
upstream end.

**Not in this part.** The manufactured solution of SST (moved to 6g with
the other model MMS), the D166 features as named variables of the
expression language (with their first consumer, D104), turbulence values
of the backflow at outlets (the cell's values are used; BC028 says so),
and curvature correction and QCR (D170, part 6h).

## k-ε family (part 6c)

**The models** (D131). The standard k-ε of Launder & Spalding (C_μ 0.09,
C_ε1 1.44, C_ε2 1.92, σ_k 1.0, σ_ε 1.3) and the realizable k-ε of Shih
et al. (1995): variable C_μ = 1/(A₀ + A_s U* k/ε) with A₀ = 4.0, the
destruction C₂ ε²/(k + √(νε)), C₁ = max(0.43, η/(η + 5)), σ_ε 1.2. One
kernel (`turb::keps`) for both, generic over the number type, so the
dual-number Jacobian sees the closures; the arccos of A_s is a helper of
`foehncraft-core::math` with its derivative. The standard model takes the
Durbin bound by default (inactive in equilibrium, it keeps the impulsive
start of near-wall cells bounded); Kato–Launder, the clip
min(P_k, 10ε) and `none` are options. The ε destruction is linearised
with the Picard form by default: the Newton form removes at most half
of ε per iteration and let k collapse in the decaying free stream of the
plate. Inlets and initial values take ε directly or ω (converted with
ε = C_μ k ω); a restart file of the other family is converted the same
way.

**The walls** (D134). All treatments take u* = C_μ^¼ k_P^½ and
y* = u* y_P/ν, with no iteration. `scalable`: the log law at
ỹ* = max(y*, 11.06), ε_P = u*³/(κ ỹ), P_k,P = (τ_w/ρ) u*/(κ ỹ).
`standard`: the log law from y* = 11.06, the viscous sublayer below it.
`compound`: the blend of Popovac & Hanjalić with Γ = 0.01 y*⁴/(1 + 5y*)
for u⁺, ε_P and P_k,P. ε_P is fixed in the near-wall cell (weight 1) and
written into the field before each solve; k and ε have zero gradient at
the wall. Rough walls shift the log law by ΔB of Cebeci & Bradshaw with
C_s = 0.5.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| TMR flat plate (`tmr_plate_keps.toml`), standard k-ε, `scalable`, 3 072 cells, sustaining terms at the inlet values, first-cell y⁺ 1 / 2 / 5 / 10 / 30 / 60 / 100 at x = 0.97 | C_f(0.97) against White, 2.875e-3 | ±5 % | −2.6, −1.7, −1.0, −0.3, −1.6, −1.5, −0.8 % |
| | y⁺ insensitivity against the y⁺ ≈ 1 result | ±5 % (scalable from y⁺ 30) | +1.1 to +2.4 % (every y⁺ within 2.4 %) |
| The plate, standard k-ε, `compound` | C_f(0.97) for y⁺ 30 / 60 / 100 | ±5 % | −2.3, −3.6, −3.4 % |
| | C_f(0.97) for y⁺ 1.3 / 2.6 / 6 / 11 | not claimed | +58, +52, +40, +21 %: **finding**, the high-Re model has no damping in the buffer layer; TURB001 below y⁺ 30 |
| The plate, realizable k-ε, `scalable` and `compound`, y⁺ 30 / 60 / 100 | C_f(0.97) against White | ±5 % | scalable −6.8, −5.2, −3.3 %; compound −7.4, −7.1, −5.9 %: **finding** |
| Channel at Re_τ = 395 (`channel_keps.toml`), standard k-ε, `scalable`, 12 uniform cells, first-cell y⁺ 32 | U_b⁺ against the DNS (17.55) | ±5 % | 18.31 (+4.3 %) |
| The channel on 4 to 40 cells (y⁺ 93 to 9), `scalable` | U_b⁺ | ±5 % from y⁺ 30 | +6.0 % (y⁺ 93), +4.3 to +4.9 % (y⁺ 47 to 9) |
| The channel, `compound`, y⁺ 93 to 11; graded mesh at y⁺ 0.4 | U_b⁺ | ±3 % at y⁺ ≈ 1 | +5.9 % (93) to +1.6 % (24), −4.4 % (16), −12.9 % (11); −32 % at y⁺ 0.4 (τ_w 2.2 times too high): **finding**, resolved walls go to ζ-f and EB (owner decision 2026-11-17) |
| The channel, realizable, `scalable`, y⁺ 93 to 9 | U_b⁺ | ±5 % | +6.0 to +10.7 %: **finding** |
| Rough pipe (`pipe_rough_keps.toml`), Re = 10⁶, k_s/D = 10⁻³ | f against Colebrook (0.01994) | ±5 % | 0.02081 (+4.4 %) |
| Moody study (`examples/studies/moody_study.py`), Re 3·10⁴ to 10⁷, k_s/D 0 to 10⁻² | f against Colebrook | ±5 % | smooth −5.0 to −1.7 %; fully rough +8 to +11 % (C_s = 0.4: +1 to +3 %); transitional −17 to +8 %: **finding** |
| Homogeneous decay in a box, k₀ = ε₀ = 1, t = 0.5, 10 / 20 / 40 steps | BDF2 time order of k and ε against the exact solution | 1.8 – 2.2 | asserted in `tests/keps.rs` |
| Sustaining terms in a box at rest, both models | k and ε constant | round-off | asserted |
| Restart of an SST run with the k-ε model | ε = C_μ k ω on the way in, ε in the new file | – | asserted |
| Plate with `compound` on 1 and 4 threads | residuals, wall shear after 8 iterations | bit-identical | asserted |
| Kernels | realizable C_μ in simple shear, the realizability bound ν_t ≤ k/(3λ_max), the sustaining terms, the source signs, dual-number derivatives; the wall treatments in the log layer and the sublayer, ΔB continuous at k_s⁺ 2.25 and 90 | – | kernel tests of `turb::keps` and `turb::wall` |

**Findings.** (1) The compound treatment with the high-Re k-ε models
overpredicts the wall shear below y⁺ 30 (up to 2.2 times at y⁺ ≈ 1): the
models have no damping in the buffer layer. The owner chose to keep it
as a finding and leave resolved walls to ζ-f (6d) and EB (6e); TURB001
warns. (2) The realizable model gives 3 to 7 % less wall friction than
the correlations on the plate and 6 to 7 % more U_b⁺ in the channel, in
both cases in the same direction; the log-layer consistency of its
constants holds (κ² = σ_ε(C₂√C_μ − C₁) = 0.168), so the near-wall cell,
where its C_μ responds to the poorly resolved velocity gradient, is the
first suspect; open for 6g. (3) With the default C_s = 0.5 the fully
rough friction factor lies 8 to 11 % above Colebrook; C_s = 0.4 puts it
within 3 %. In the transitional range the Nikuradse-type ΔB cannot follow
Colebrook's monotonic curve (up to 17 % below). (4) With the TMR inlet
values and no sustaining terms the free stream of the k-ε family decays
(the realizable k to zero): the plate runs with the sustaining terms.

**Not in this part.** The improved blending of the SST `automatic`
treatment and a sweep of the SST wall value of ω, both left open by 6b,
move to 6g with the model comparison; the GCI of the k-ε plate and
channel also follows there.

## k-ε-ζ-f and the scalar multigrid (part 6d)

**The model** (D132). Hanjalić, Popovac & Hadžiabdić (2004):
ν_t = C_μ ζ k T, C_μ = 0.22, C_ε1 = 1.4(1 + 0.012/ζ), C_ε2 = 1.9,
C₁ = 1.4, C₂′ = 0.65, σ_k 1, σ_ε 1.3, σ_ζ 1.2, C_T 6, C_L 0.36, C_η 85,
the realizability coefficient 0.6 (constants checked against the page of
the turbulence modelling benchmark working group, TMBWG). The bounds of T
and L take the strain magnitude S = (2 S_ij S_ij)^½, as the authors'
implementation does (6d amendment of specification 2.1, which read
(S_ij S_ij)^½: that form gave U_b⁺ +2.9 % in the channel against +1.9 %).
C_ε1 P is evaluated as 1.4(P + 0.012 P/ζ) with P/ζ = C_μ k T S², finite
as ζ → 0. One kernel (`turb::zeta_f`) gives T, L², ν_t, the source split
of k, ε, ζ and the right-hand side of f, generic over the number type
(dual-number test). Solve order per outer iteration: k and ε; the cell
values again with the new k and ε, then f (CG); the cell values again,
then ζ with the new f, bounded to [10⁻¹⁰, 2]; ν_t at the end. A transport
solve may lower k, ε or ζ by at most a factor 10 (`ZETA_F_DROP`; 6d
amendment): an unconverged BiCGStab iterate of the loose steady tolerance
can go negative where the solution is small, and the absolute floor
alone put ε near zero, which blew up T, L and P/ε in the source of f
(the plate at y⁺ 0.5 diverged in iteration 36).

**The walls** (D134). `low-re`: k = ζ = 0 (Dirichlet), ε_w = 2νk₁/d₁²
and f_w = −2νζ₁/d₁² on the wall faces. `compound` (default): k and ζ by
zero gradient (dk/dy = dζ/dy = 0 at the wall), ε and P_k of the first
cell fixed by the blend of Popovac & Hanjalić (2007) with C_μ = 0.09 in
the log-law values, f_w = e^(−Γ)(−2νζ₁/d₁²). Roughness with `low-re` is
refused.

**The multigrid** (D139). The backends hold hierarchies in slots
(0 pressure, 1 f, 2 k, 3 ε or ω, 4 ζ). The symmetric positive definite f
equation runs CG with its own hierarchy by default; the transport
equations take one with `solver.turbulence.preconditioner = "amg"`,
built from the symmetric part ½(A + Aᵀ) and preconditioning BiCGStab.
All rebuild on the schedule of the pressure (`solver.amg_rebuild`), and
the matrices of their last setups go into the restart file
(`turb_amg_matrix_<slot>`), so a continued run builds the same
hierarchies. A rebuild in every solve (the first version) cost 21 s of
the 56 s of the channel: on a small mesh the coarsest level holds all
cells and its dense inverse dominates.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Channel at Re_τ = 395 (`channel_zeta_f.toml`), 96 cells graded 60, first-cell y⁺ 0.27, `compound` | U_b⁺ against the DNS (17.55) | ±3 % | 17.89 (+1.9 %), 3 700 iterations |
| The channel, `low-re` | U_b⁺ | ±3 % | 17.79 (+1.4 %) |
| TMR flat plate (`tmr_plate_zeta_f.toml`), `compound`, 64 × 48 cells, first-cell y⁺ 1, turbulent start | C_f(0.97) against White, 2.875e-3 | ±3 % | 3.023e-3 (+5.1 %): **finding** |
| The plate on 64 × 96 cells (y⁺ 1) | C_f(0.97) | ±3 % | +4.7 %: mesh-independent, **finding** |
| The plate with an inlet of Tu 0.8 % (k = 10⁻⁴, ν_t/ν = 10), 48 cells | C_f(0.97) | – | +5.0 % |
| The plate from the TMR free-stream values (ν_t/ν = 0.009) | transition | – | laminar over the whole plate, C_f = 3.0e-4: **finding** |
| The plate at y⁺ 0.5 | stability | – | diverged without the relative bound; converges with it |
| Run time of the channel (192 cells, 3 000 iterations) | multigrid rebuild per solve against the pressure schedule | – | 56 s against 26 s (SST: 16 s) |
| Homogeneous decay, k₀ = ε₀ = 1, ζ₀ = 0.3, t = 0.5, 10 / 20 / 40 steps | BDF2 order of k, ε and of the return of ζ to 2/3 | 1.8 – 2.2 | asserted in `tests/zeta_f.rs` |
| Sustaining terms in a box at rest | k, ε constant, ζ = 2/3, f = 0 | 10⁻¹⁰ | asserted |
| Transient restart with the multigrid of all four equations | fields at the end | bit-identical | asserted |
| The plate on 1 and 4 threads (20 480 cells, multigrid of all equations) | residuals, wall shear after 8 iterations | bit-identical | asserted |
| A second hierarchy slot | the V-cycle of the pressure slot | bit-identical | asserted (`tests/stage3.rs`) |
| Kernels | T, L and ν_t in simple shear, the decay and sustaining limits, the source signs, the low-Re wall values, dual-number derivatives | – | kernel tests of `turb::zeta_f` |

**Findings.** (1) On the TMR plate the model gives C_f 5 % above White's
correlation, independent of the mesh (48 or 96 cells across), outside
the ±3 % band; in the channel it lies within 2 %. (2) From the TMR free
stream (ν_t/ν = 0.009) the boundary layer stays laminar over the whole
plate: the length scale of f in the weak free stream,
C_L C_η (ν³/ε)^¼ ≈ 7 mm, exceeds the laminar layer, and f drives ζ, and
with it ν_t, to zero across it; the case starts from a turbulent field.
(3) In the laminar part near the leading edge the published sink −ε of
the k equation has no positive steady solution once T sits at its
Kolmogorov bound (ε then no longer falls with k): the floor holds k in
about 120 cells and its residual stays at 7·10⁻⁶ while every other field
converges to round-off (TURB002). The ζ-f₀ variant of Popovac (2023),
with f₀ = 0 at the wall, follows as the next part for robustness in y⁺.

## k-ε-ζ-f₀ and the compound treatment of 2023 (part 6d2)

**The model** (D132, specification 3.3.1). Popovac (2023): f₀ = f +
(N_ζ − 1)ζ/T in place of f, N_ζ = 2, so f₀,w = 0;
ν_t = (C_μ/ζ_eq) ζ k T with C_μ = 0.09, ζ_eq = 0.41;
C_ε1 = 1.4(1 + 0.012/√ζ); the other constants of 2004. The same kernel
as ζ-f (`turb::zeta_f`, `FORM_F0`). The author's code is not available
(maintainer, 2026-11-17), so the default was chosen on the validation:

| Variant | Channel U_b⁺, y⁺ 1 | Plate C_f, y⁺ 1 | Step x_r/H |
|---|---|---|---|
| As printed (sign of (N_ζ − 1) in eq. 7 corrected, k sink −ε) | −11 % | +21 % (TMR), floor holds k | 7.34 (Gmsh mesh) |
| C_ε1 with 1/ζ (2004) | +8.5 %, not converged | – | – |
| k sink k/T with the realizability bound | +30 % (relaminarizes) | – | – |
| k sink k/T_K, a = 0.6 | +1.3 % | −4.0 % (White) | 7.67 |
| Default of 6d2 to #140: k sink k/T_K, a = 0.6 √2 | +1.4 % (not grid-convergent, #128) | −2.5 % (White) | 6.38 |
| **Default since #141: k sink max(k/T_K, 2νk/d²), C_ε1 = 1.4 (1 + 0.009/ζ), a = 0.6 √2** | **−1.2 %; +0.2 % mesh-converged** | **+1.8 % (White)** | **6.44** |

**The wall treatment** (D134, `compound-2023`). u* = C_μ^¼ k_P^½,
y* = u* y_P/ν; u⁺ = [y*^(−n) + u⁺_log^(−n)]^(−1/n), n = 10/3,
u⁺_log = ln(8.9 max(y*, 1))/(κψ); ψ = 1 − C_U max(y_P − y_v/2, 0)/(U_t κ u*)
with y_v = 11ν/u* and C_U = t·[(U·∇)U + ∇p/ρ] of the wall cell (steady
part), bounded to [0.5, 2] (the paper gives no bounds); ν_w = ν y*/u⁺.
Γ = du⁺/dy* from the blend, ε_P = Γ 2νk_P/y_P² + (1 − Γ) u*³/(κ y_P)
and P_k,P = u*² (τ_w/ρν) Γ, fixed in the cell; k and ζ by zero flux.
It also serves ζ-f (option).

**The mesher.** `mesh.generate.remove` takes the cells whose centroid lies
in a box out of a hexahedral, tetrahedral or mixed mesh and names the
exposed faces; the step case uses it.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Channel at Re_τ = 395 (`channel_zeta_f0.toml`), 48 cells, first-cell y⁺ 1 | U_b⁺ against the DNS (17.55) | ±3 % | 17.34 (−1.2 %), converged (#141; 6d2: 17.79) |
| The channel on 96 / 192 / 384 cells (#141) | U_b⁺ | GCI 0.5 % | 17.54 / 17.58 / 17.58 (+0.2 %), GCI 0.02 % |
| The channel at y⁺ 10 / 30 | U_b⁺ | – | −3.4 / −1.2 % (#141; 6d2: −4.2 / −1.5 %; ζ-f, 2007 treatment: −7.1 / +4.8 %) |
| TMR plate (`tmr_plate_zeta_f0.toml`) from the TMR free stream, y⁺ 1 | C_f(0.97) against White, 2.875e-3; convergence | ±3 % | 2.926e-3 (+1.8 %), converged to 10⁻⁷, no floor (#141; 6d2: −2.5 %) |
| The plate at y⁺ 10 / 30 | C_f(0.97) against White | ±3 % | −0.3 / −2.7 %, converged (#141; 6d2: +0.9 / −2.0 %) |
| Backward-facing step (`backstep_zeta_f0.toml`), 8 200 cells | x_r/H against 6.26 | ±10 % | 6.44 (+3 %; 6d2: 6.38) |
| The step on 32 800 cells | x_r/H | ±10 % | 6.50 (+3.8 %, coupled, #141), within 1 % of the coarse mesh (6d2: 6.35) |
| The step with ζ-f (2004, `compound`) | x_r/H | ±10 % | 7.82 (+25 %; 6.80 with `a = 0.849`) |
| The step with SST (`automatic`) on 8 200 / 32 800 cells | x_r/H | ±10 % | 6.81 / 6.80 (+9 %); the 5.79 reported first came from 4 000 iterations at α_u = 0.5, before the bubble had settled |
| The step with k-ε (`scalable`), 32 800 cells | x_r/H | – | 5.36 (−14 %) |
| The plate on 1 and 4 threads (20 480 cells) | residuals, wall shear after 8 iterations | bit-identical | asserted |
| Kernels | f₀ form against the 2004 source, the limits of the blend (y⁺ 30 to 300), Γ as the derivative of the blend, ψ, dual numbers | – | kernel tests of `turb::zeta_f` and `turb::wall` |
| Mesher | step volume, exposed area, unused nodes, taken names | exact | `removed_box_builds_a_step` |

**Findings.** (1) The realizability bound of T with a = 0.6 and
S = (2 S_ij S_ij)^½ is active in the separated shear layer of the step
and halves ν_t there; the bubble then grows to 7.4 to 7.7 H with both
ζ-f models. (2) The published k sink −ε has no positive steady k in the
laminar part of a boundary layer (the 6d finding); k/T_K removes it,
and since #141 max(k/T_K, 2νk/d²) also converges with the mesh at the
wall (see the GCI section).
(3) C_p, C_f and the profiles of the step against the experiment and a
GCI follow with the stage 6 validation. (4) Iterations of the step
(measured after 6d2, fields every 50 to 250 iterations): with the
momentum under-relaxation 0.5 throughout, x_r is within 1 % of its final
value after 950 (ζ-f₀), 1 250 (k-ε) and 1 500 (SST) iterations on 8 200
cells, and after about 2 900 (ζ-f₀) and more than 4 000 (SST) on 32 800
cells. With 0.5 for the first 100 iterations and 0.8 after them it takes
550 to 600 and about 1 750; 0.7 from the start works as well, 0.9 stalls
the residual at 10⁻⁴, and a turbulence relaxation of 0.95 oscillates. The
count grows about three times when the mesh spacing halves: the
relaxation sets a local pseudo time step of a few cell transit times, and
the bubble needs a fixed time to form. The scaled residual is 5·10⁻⁵ to
2·10⁻⁴ when x_r has settled; a residual tolerance does not show when the
step result is done.

## Coupled algorithm and pseudo-time steps (part 6e)

Specification: [requirements step 22](../requirements/step22-coupled-solver.md)
(D176 to D178), which also holds the measurements behind the choices.

**The algorithm** (D176, `foehncraft-flow::coupled`). One outer
iteration solves

    | A  G | |δu|   |r_u|
    | D  L | |δp| = |r_c|

for the increments, with the momentum rows A, the residuals and the
Rhie–Chow fluxes of SIMPLEC; then u += δu, p += δp, the fluxes follow
from the same increments, and one pressure correction removes the
continuity defect that the loose linear solve (0.1) leaves. The operator
is applied with the operations of the device (no block matrix); flexible
GMRES, right-preconditioned with a SIMPLE-type sweep: BiCGStab on the
momentum rows to 10⁻², the pressure matrix with ρ(D̄ + D̂̄) by CG and the
multigrid of the pressure, the velocity correction with
D̂ = min(D^C, V/((1 − α) a_P/α + ρV/(αΔτ))) (kernel `schur_diagonal`).

**The pseudo-time steps** (D177, D178). Momentum: ρV/Δτ (u − uᵏ) on the
rows through the time term of the transient solver and the
time-step-independent Rhie–Chow flux (D27) with the current iterate as
the old level; Δτ = `factor` × ΣV/ΣA_b / U, one value for all cells;
α_u = 0.98 as the cap on the local Courant number. Turbulence:
ρV max(1/Δτ, 1/(c_T T)) (φ − φᵏ) on the scalar rows (kernel
`pseudo_coefficient`), c_T = 0.5, relaxation 0.98; the residual scale
does not contain the term.

### Results

Backward-facing step of part 6d2, uniform start, iterations until x_r
stays within 1 % of its final value (a result every 5 iterations for the
coupled runs, every 25 for SIMPLEC, every 50 on the fine mesh; the fine
runs went on to 400 and 3 000 iterations for the final value):

| Test | Check | SIMPLEC (0.5, then 0.8 from 100) | Coupled |
|---|---|---|---|
| Step, k-ε-ζ-f₀, 8 200 cells (`backstep_zeta_f0_coupled.toml`) | iterations to 1 %; final x_r/H | 525; 6.384 | 65; 6.385 |
| Step, k-ε (`scalable`), 8 200 cells | the same | 600; 5.287 | 75; 5.287 |
| Step, k-ω SST, 8 200 cells | the same | 675; 6.82 | 80; 6.81 |
| Step, k-ε-ζ-f₀, 32 800 cells | the same | 1 600; 6.352 | 175; 6.352 |
| Channel at Re_τ = 395, k-ε-ζ-f₀ (`channel_zeta_f0.toml`) | iterations to 10⁻⁶; U_b⁺ | about 2 100; 17.79 | 323; 17.79 |
| Cavity at Re = 1000, 128² (`cavity_re1000_coupled.toml`) | iterations to 10⁻⁸; deviation from Ghia (u, v) | 1 397 (α = 0.9); 0.0040, 0.0136 | 452; 0.0040, 0.0136 |
| The same cavity, kinetic-energy ledger (D88) | numerical dissipation over lid power; its pressure part (W) | 0.637 %; 3.996223 × 10⁻⁶ | 0.637 %; 3.996224 × 10⁻⁶ |

| Test | Check | Band | Result |
|---|---|---|---|
| Cavity 16², Re = 100, converged to 10⁻¹¹ | lid force against SIMPLEC, for the automatic step, Δτ = 0.02 s, `factor = 20` with α = 0.9, and relaxation only | 10⁻⁸ | asserted (`coupled_result_is_the_simplec_result`) |
| Cavity 48², Re = 400 | iterations: SIMPLEC at α = 0.9; coupled with `factor` 1, 4, 20 | factor 4 at most 0.7 of SIMPLEC, 20 fewer than 4 | 255; 292, 130, 84 |
| Cavity 16², 10 iterations, lid speed 1 and 2 m/s | the automatic step (L = 0.25 m over the lid speed; the symmetry planes of the one cell layer must not enter); `factor`, `step`, `momentum = false` | 0.25 s and 0.125 s to 10⁻¹² | asserted (`automatic_step_is_the_length_over_the_prescribed_speed`) |
| Cavity 16², 24 iterations, restart at 10 | velocity, pressure, fluxes, multigrid matrix of the continued run | bit-identical | asserted |
| Cavity 150², 8 iterations, 1 and 4 threads | residuals, Krylov count, boundary forces | bit-identical | asserted |
| Input | `[time]` with `coupled`, relaxation 1, ranges of `[solver.coupled]` and `[solver.pseudo_time]`, the defaults | refused with CASE004 / CASE005 | asserted |
| Kernels | pseudo-time term of the scalar row and its residual scale, `pseudo_coefficient`, `schur_diagonal` | – | kernel tests |

Run time of the step with k-ε-ζ-f₀ until x_r has settled, one thread:
21 s against 31 s on the coarse mesh (65 against 525 iterations), about
235 s against 327 s on the fine one (175 against 1 600). A coupled
iteration costs five to seven SIMPLEC iterations (10 and 14 Krylov
iterations, each with a momentum sweep of 10 iterations and a pressure
solve), so the run time gains a factor of 1.4 to 1.5 where the count
gains 8 to 9.

**Findings.** (1) The coupled solve alone, with relaxation as before,
halves the count (300 on the coarse step, 700 to 800 on the fine one):
the larger part of the gain is the pseudo-time step. (2) On the step the
count does not depend on Δτ between 5 and 29 H/U, nor on having a global
step at all once the cap is 0.98; the global step makes the linear
solves cheaper (8 instead of 18 Krylov iterations per outer iteration).
(3) The cap on the local step sets the count on fine meshes and on the
cavity (Δτ given; 0.95: 1 012, 0.98: 423, 0.999: 166 iterations); a tighter
coupled solve does not (167). At 0.999 the momentum sweeps cost four
times the run time, so 0.98 is the default. (4) The momentum sweep must
reach 10⁻² to 10⁻³; truncated sweeps make the Krylov count and x_r
erratic. On the momentum matrix of the step at Δτ = 5 H/U, BiCGStab
needs 113 to 137 iterations with a diagonal preconditioner, 21 to 24
with symmetric Gauss–Seidel, 90 to 194 with the aggregation multigrid of
the symmetric part and 1 with a multigrid of approximate ideal
restriction (L44). (5) Turbulence: a global step is limited to 0.5 H/U;
with the time-scale rule c_T = 0.25, 0.5, 1 gives 75, 55, 45 to 50
iterations (cap 0.999); relaxation 1 lets k-ε diverge at the start, 0.9
makes the turbulence lag the momentum equations (x_r overshoots to 11 H
and has not settled after 150 iterations). (6) With SIMPLEC the
time-scale rule alone (`solver.pseudo_time.turbulence = 0.5`) gives 230
iterations at α_u = 0.8; it is off by default there. (7) On small
laminar cavities the coupled algorithm is not faster than SIMPLEC at
α = 0.9: the automatic step (0.25 L/U) is a careful one, and larger
steps cost more per iteration than they save (cavity at Re = 1000 on
64²: 584 iterations in 8 s with SIMPLEC; 625 in 13 s, 276 in 24 s and
178 in 43 s with `factor` 1, 4 and 20). The momentum sweep is the limit
in every one of these findings. (8) A run repeats bit by bit for the
same input, but its count reacts to the last digit of an input: the
cavity needs 423 iterations with Δτ = 0.25 s given and 452 with the
automatic step, which differs from 0.25 s in the 17th digit; the fine
step needed 145 iterations with a slightly shorter step (a velocity
scale that counted the symmetry planes) and needs 175 at
7.3 H/U. The
coupled solve stops at 0.1, its momentum sweep is cut off after 10
iterations and the multigrid cycle runs in single precision, so two runs
separate in the fifth digit within one outer iteration. The converged
result does not change; the counts of this section are good to about
10 %. The histories of the figure of gallery 78 come from intermediate
results every 5 iterations (every 25 with SIMPLEC). An intermediate
result leaves a run unchanged (stage 5, control files; test
`intermediate_results_and_balances_do_not_change_the_run`), so they are
the histories of the runs without them: the fine step gives x_r/H =
6.27, 6.34, 6.37 after 170, 180, 190 iterations either way. Before that
correction a result file shifted the path within this scatter (the SST
count read 85 instead of 80).

## Multigrid for the momentum sweep (part 6f)

Specification: [requirements step 22](../requirements/step22-coupled-solver.md)
sections 3.1 and 7.7 (D179).

**The method** (`foehncraft-fv::amg::reduction`). The momentum sweep of
the coupled preconditioner is one V-cycle of a reduction multigrid per
velocity component instead of BiCGStab. Setup on the host, per level:
strong dependencies (|a_ij| ≥ 0.3 max_{k≠i} |a_ik|); the first pass of
the Ruge–Stüben splitting into coarse points C and fine points F; the
restriction R = [Z I] with Z = −A_CF (2 D_F⁻¹ − D_F⁻¹ A_FF D_F⁻¹), an
approximate ideal restriction; direct interpolation; A_c = R A P with
the entries below 10⁻⁴ of the diagonal dropped; a dense inverse at 64
unknowns or fewer. The cycle is the one of the scalar multigrid on the
same level operations (CPU and GPU, float32 hierarchy), with an inverse
diagonal that is zero on the C points: Jacobi on the F points, once
before and once after the coarse correction.

**In the coupled algorithm** (`foehncraft-flow::coupled`). One hierarchy
per velocity component, shared by components with the same matrix and
built when a sweep first needs it; kept while the momentum matrix has
changed by at most 2 % since its setup (`momentum_reuse`) and while its
cycle still halves the residual of the first sweep of an outer
iteration; stored with its matrix in the restart file. A component
whose part of a Krylov vector is zero up to rounding takes no sweep. A
failed setup leaves the sweeps of that outer iteration to BiCGStab, with
a warning in the log.

**One thread without a pool** (`backend::cpu`, `device::cpu`). A pool of
one thread and loops of one task run on the calling thread instead of
being handed to a worker (same values); found with the measurements
below.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Splitting of a 5-point Laplacian and of a convection matrix | C and F points; every F point has a strong C neighbour; coarsening along the flow | – | asserted (`splitting_of_a_laplacian`, `splitting_follows_the_flow`) |
| Convection–diffusion in a solid-body rotation, 32² to 256², cell Péclet number 0 to 300 | factor per cycle with the diagonal raised by 2 % (the cap 0.98); without; operator complexity | 0.35; 0.5 to 0.7; 4 | 0.07 to 0.22 (0.15 and 0.16 on 128² and 256²); 0.25 to 0.53; 2.4 to 3.5 |
| Scalar transport on hexahedra, tetrahedra, polyhedra and a mixed mesh in 3D, cell Péclet number 0.5 and 20 (`tests/reduction.rs`, on the backend, float32 hierarchy) | residual after 8 cycles; operator complexity | 10⁻³; 5.5 | passes; 2.2 to 5.4 |
| The same hierarchy in two multigrid slots, float32 and float64, repeated | cycles | bit-identical repeat; the precisions agree to 10⁻⁵ | asserted (`precision_and_slots`) |
| Setup on 1 and 4 threads | hierarchy | identical | asserted (`setup_does_not_depend_on_threads`) |
| Matrix without couplings; dense inverse of a non-symmetric matrix | cycle solves; A A⁻¹ = I | 10⁻¹² | asserted |
| Cavity 24², converged to 10⁻¹¹ | lid force with the multigrid sweep, the BiCGStab sweep, two cycles, `momentum_reuse = 0`, a float64 hierarchy | 10⁻⁸ | asserted (`momentum_sweeps_give_the_same_result`) |
| Cavity 12 × 12 × 6 with a symmetry plane (the components have different matrices) | result against the BiCGStab sweep; a hierarchy per component | 10⁻⁷ | asserted (`multigrid_sweep_with_symmetry_planes_in_3d`) |
| Cavity 16² (restart at 10 of 24) and 24² (restart at 70 of 90, with kept hierarchies) | velocity, pressure, fluxes, kept hierarchies of the continued run | bit-identical | asserted (`coupled_restart_continues_bit_identically`) |
| Cavity 150², 8 iterations, 1 and 4 threads | residuals, Krylov count, boundary forces | bit-identical | asserted |
| Cavity 70² with `max_levels = 1` (the coarsest level would exceed the dense limit) | the run goes on with BiCGStab sweeps; one warning in `flow.log` | – | asserted (`failed_multigrid_setup_falls_back_to_bicgstab`) |
| Input | ranges of `momentum_cycles`, `momentum_reuse`, `[solver.coupled.multigrid]` | refused with CASE005 | asserted |

Runs on one thread, to the iteration from which x_r stays within 1 % of
its final value (step) or to 10⁻⁸ (cavity); Krylov iterations per outer
iteration in brackets:

| Case | SIMPLEC | Coupled, BiCGStab sweep | Coupled, multigrid sweep |
|---|---|---|---|
| Step, k-ε-ζ-f₀, 8 200 cells (`backstep_zeta_f0_coupled.toml`) | 525 iterations, 18.7 s | 65 (10.0), 8.3 s | 65 (5.3), 7.0 s |
| Step, k-ε-ζ-f₀, 32 800 cells | 1 600 iterations, 243 s | 175 (14.2), 131 s | 150 (6.0), 74 s |
| Cavity at Re = 1000, 128² (`cavity_re1000_coupled.toml`) | 1 397 iterations, 28.3 s | 452 (2.2), 18.8 s | 449 (2.1), 12.1 s |
| Step, k-ε (`scalable`) and k-ω SST, 8 200 cells | 600 and 675 iterations | 75 and 80 | 75 and 80 (2.2 and 5.1) |
| Channel at Re_τ = 395, k-ε-ζ-f₀ | about 2 100 iterations | 323 | 324 |

Final values as in part 6e: x_r/H = 6.385 (k-ε-ζ-f₀), 5.287 (k-ε), 6.81
(SST), 6.35 on the fine mesh; the cavity deviates from Ghia et al. by
0.0040 and 0.0136 with both sweeps.

**Findings.** (1) One cycle per component replaces a BiCGStab sweep of
10 iterations at a third of its cost, and the Krylov count per outer
iteration halves on the step (so do the pressure solves). (2) The fine
step settles in 150 iterations instead of 175 and in 74 s instead of
131 s; against SIMPLEC the run time gains a factor of 2.7 on the coarse
step, 3.3 on the fine one and 2.3 on the cavity. (3) A setup costs as
much as twenty to thirty cycles, and while the flow forms the matrix
changes by more than 2 % in most iterations: 96 setups in the 65
iterations of the coarse step (1.4 s of 7.0 s, against 0.3 s for all
cycles), 140 in 150 on the fine one, 28 in 449 on the cavity. A setup
in every outer iteration made the cavity slower than with BiCGStab,
which is why hierarchies are kept. Measured on step matrices, a
hierarchy from a matrix 1.3 %, 2.7 % and 5.4 % away converges with 0.27,
0.34 and 0.6 per cycle; from 6.4 % on it does not. (4) Without the cap
on the local step (relaxation 0.999) the fine step settles in about 120
iterations, but every outer iteration runs into the limit of 30 Krylov
iterations: the cap stays at 0.98 until the Schur complement
approximation is better. (5) The BiCGStab sweeps took a quarter of the
run time on the fine mesh; the cycles and their setups take 16 %, and
the turbulence model 35 % (its linear solves 21 %) and the pressure
solves 25 %. (6) Thread hand-off: with a pool of one
thread every kernel call went to the worker thread and back, about
2 000 times per outer iteration of the cavity. Running such loops on
the calling thread took the coarse step with SIMPLEC from 31 s to
18.7 s and with the coupled algorithm and the BiCGStab sweep from 21 s
to 8.3 s, without changing a value. The times of part 6e above are from
before that change. (7) Larger pseudo-time steps on small laminar
cases, which cost two to three times the default with the BiCGStab
sweep (finding 7 of part 6e), are now as cheap as it: cavity at
Re = 1000 on 64², SIMPLEC 584 iterations in 3.0 s, coupled with
`factor` 1, 4, 20: 624, 275, 170 iterations in 3.7 s, 2.8 s, 2.6 s.

## Multigrid for the turbulence equations (part 6g)

Specification: [requirements step 22](../requirements/step22-coupled-solver.md)
sections 4.2 and 7.8 (D180).

**The method** (`foehncraft-flow::turbulence`). With the coupled
algorithm and pseudo-time steps (steady runs), the transport equations
of the k-ε family (k, ε and ζ of k-ε, realizable k-ε, k-ε-ζ-f and
k-ε-ζ-f₀; multigrid slots 2 to 4) take one V-cycle of the reduction
multigrid of part 6f as the preconditioner of BiCGStab instead of the
diagonal (`[solver.turbulence] preconditioner = "reduction"`, the
default there). A hierarchy is kept while the matrix of its equation
has changed by at most 2 % since its setup (`reuse`, the measure of
`momentum_reuse`), and dropped when a solve with it needs more than 8
BiCGStab iterations. A failed setup leaves the solves of that equation
to Jacobi until a setup succeeds, with one warning in the log. The
restart file holds the matrices of the kept hierarchies
(`turb_amg_matrix_<slot>`). The elliptic f equation keeps CG with the
aggregation multigrid.

**The steps.** With the multigrid the relaxations default to 0.995
(`solver.turbulence.relaxation`, `solver.momentum.relaxation`) instead
of 0.98. A relaxation α acts as a local pseudo step α/(1−α) ρV/a_P,
about 25 h²/ν_eff in a flat cell at the wall at 0.98: it shrinks with
the cells, the rule min(Δτ, 0.5 T) of D178 does not. k-ω SST keeps
Jacobi and 0.98, laminar cases 0.98, transient runs and SIMPLEC are
unchanged; a key set in the case wins.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Channel at Re_τ = 395, k-ε-ζ-f₀, coupled, to 10⁻¹⁰ | driving force with the multigrid at 0.995 against Jacobi at 0.98; setups counted and fewer than three per iteration; fewer iterations; SST keeps the old defaults | 10⁻⁶ | asserted (`turbulence_multigrid_gives_the_same_result`) |
| The same channel, restart at 30 of 60 iterations (hdf5) | velocity, pressure, turbulence fields, kept hierarchies in the restart file | bit-identical | asserted (`turbulence_multigrid_restart_continues_bit_identically`) |
| Input | `reuse` outside [0, 0.2] | refused with CASE005 | asserted |

Runs on one thread, quiet machine, to the iteration from which x_r stays
within 1 % of its final value (step) or to the tolerance (channel);
part 6f in brackets:

| Case | Iterations | Run time |
|---|---|---|
| Step, k-ε-ζ-f₀, 8 200 cells (`backstep_zeta_f0_coupled.toml`) | 40 (65) | 7.5 s (7.0 s) |
| Step, k-ε-ζ-f₀, 32 800 cells | 65 (150) | 56 s (74 s) |
| Step, k-ε (`scalable`), 8 200 cells | 45 (75) | 5.7 s (4.5 s) |
| Channel at Re_τ = 395, k-ε-ζ-f₀ | 231 (323) | – |
| Step, k-ω SST; cavity at Re = 1000 | unchanged: 80; 449 | unchanged |

Final values: x_r/H = 6.385 and 6.353 (k-ε-ζ-f₀, coarse and fine), 5.287
(k-ε), as before. The step case file converges to its tolerance in 117
iterations with every check passed.

**Findings.** (1) The turbulence relaxation 0.98 made the outer count
grow with the mesh (65 to 150 iterations per four times the cells): on
the fine step, momentum caps of 0.99 and 0.995 gave 150 and 145
iterations, `pseudo_time.factor = 2` 140, `pseudo_time.turbulence = 1`
125, the turbulence relaxation 0.99 95. (2) With Jacobi a relaxation of
0.99 or more runs the BiCGStab solves of k, ε and ζ into their limit of
50 in every outer iteration and x_r jitters or runs away; with the
multigrid they need 1 to 2 iterations. (3) Both relaxations at 0.995
together give 65 iterations on the fine step, the count of the coarse
step before: the count grows by 1.6 per four times the cells instead of
2.3. (4) An iteration costs more: the momentum cap 0.995 doubles the
Krylov iterations (9 instead of 5 per outer iteration) and with them
the pressure solves, and while the flow forms the hierarchies of the
turbulence equations are rebuilt in most solves (114 setups in 120
solves on the coarse step, 172 in 195 on the fine one; a setup costs
about as much as twenty cycles). The coarse steps take about the same
time as before; the fine step is 1.3 times faster, and the gain grows
with the mesh. Against SIMPLEC the fine step is now 4.3 times faster
(243 s). (5) SST on the step: with accurate solves of k and ω (the
multigrid, or Jacobi to 10⁻⁴) the residuals grow again after about 160
iterations and x_r wanders from 6.81 to 6.85 H at any relaxation from
0.95 to 0.995; the loose Jacobi solves of today damp this. A relaxation
of ν_t of 0.2 cures it at the caps 0.98 (206 iterations to the
tolerance, x_r = 6.809 H) but not at the larger ones, so SST keeps the
old defaults. (6) The aggregation multigrid of the scalar equations
(`preconditioner = "amg"`, part 6d) fails with pseudo steps: its
coarsest level is not positive definite.

## The pressure solve of the coupled preconditioner (part 6h)

Specification: [requirements step 22](../requirements/step22-coupled-solver.md)
sections 7.9 and 7.10 (D181).

**The method** (`foehncraft-flow::coupled`). The pressure solve inside
the SIMPLE-type preconditioner of the coupled algorithm stops at the
relative tolerance 0.3 instead of 0.1 (`pressure_tolerance`): the
Krylov iteration corrects what it leaves. The aggregation multigrid of
its matrix Ŝ is built anew when Ŝ has changed by more than 5 % since
the setup (`pressure_reuse`, the measure of `momentum_reuse`), besides
every `amg_rebuild` iterations as before. The restart file and the
state kept across an intermediate result hold the matrix of the last
setup, as before.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Cavity 24², Re = 100, coupled, to 10⁻¹¹ | a new hierarchy whenever Ŝ changes (`pressure_reuse = 0`) against the periodic rebuild alone (`1`): more setups, the periodic count exactly every 10 iterations, the same lid force | 10⁻⁸ | asserted (`pressure_hierarchy_follows_its_matrix`) |
| Coupled tests of parts 6e to 6g, restart included (hdf5) | results with the new defaults | as before | asserted |
| Input | `pressure_reuse` outside [0, 1] | refused | asserted |

Runs on one thread, quiet machine, to the iteration from which x_r stays
within 1 % of its final value; part 6g with the tolerance 0.1 in
brackets:

| Case | Iterations | Run time | Pressure solves |
|---|---|---|---|
| Step, k-ε-ζ-f₀, 8 200 cells (`backstep_zeta_f0_coupled.toml`) | 40 (40) | 6.2 s (8.0 s) | 0.69 s (2.02 s) |
| Step, k-ε-ζ-f₀, 32 800 cells | 65 (65) | 39.8 s (56.9 s) | 5.81 s (20.0 s) |
| Step, k-ε (`scalable`), 8 200 cells | 45 (45) | 4.7 s (6.4 s) | 0.85 s (2.03 s) |
| Step, k-ω SST, 8 200 cells | 80 (80) | 5.4 s | 0.82 s |
| Cavity at Re = 1000 | 442 to the tolerance (449) | – | – |

The step case file converges to its tolerance in 115 iterations with
every check passed.

**Findings.** (1) Two published Schur complements were measured against
the SIMPLE-type one and do not pay on the step: SIMPLER (a pressure
prediction first) saves a few Krylov iterations for a second pressure
solve, the least-squares commutator runs into the Krylov limit near the
inlet, outlet and walls (section 7.9). (2) The pressure solves, not the
Schur complement, were the cost: a hierarchy up to nine outer
iterations old needed 9 to 20 CG iterations per solve where a fresh one
needs two. (3) The tolerance 0.3 alone gains 8 % on the coarse meshes
and nothing on the fine one (its solves stop at the limit of 20 either
way); with the hierarchy following its matrix the fine step needs 71 CG
iterations per outer iteration instead of 160 and runs 1.4 times
faster, 6.1 times faster than SIMPLEC (243 s). (4) Laminar cases do not
trigger the rule: the pseudo-time term holds D̂, and the cavity builds
the same 45 hierarchies either way.

## Thermal wall function and turbulent Prandtl number (part 6f)

Specification: [requirements stage 6](../requirements/stage6-turbulence.md)
section 5 (D135), `foehncraft-kernels::turb::thermal`.

**The method.** A passive scalar of a turbulent flow diffuses with
ρ(Γ + ν_t/Pr_t) on every face: Pr_t is constant (`turbulent_prandtl`,
0.85) or Kays & Crawford's correlation of ν_t/ν and Pr. On a wall face
the diffusivity is Γ_w = μ y⁺/T⁺, so the face conductance is
h_w = ρ u\*/T⁺ with T⁺ of Kader (1981; all y⁺, any Pr; the default) or
Jayatilleke (1969; log region, blended to Pr y⁺ by Kader's Γ). Fixed-flux
walls take b = q d_n/Γ_w, so the wall value is φ_P + q/h_w. y⁺ is that of
the momentum wall treatment (u\* of the k-ε family, the blended u_τ of
SST); with the scalable treatment T⁺ is taken at ỹ⁺ = max(y⁺, 11.06),
the virtual wall distance of its log law (part 6f, pull request 87:
before, the thermal law saw the true y⁺ and the channel Nu was +28 % on a
y⁺ ≈ 8 mesh). The default is `kader` where the specification wrote `none` with
`low-re`: the two agree for y⁺ below 1, and `kader` is also right at
larger y⁺, so no switch is needed. Both kernels are single-source
(GPU and CPU). Cases without a turbulence model do not run them.

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| T⁺ of Kader for Pr 0.01, 0.71, 7, 100 at y⁺ 1, 10, 100; Jayatilleke; Kays–Crawford | formulas evaluated apart in Python | 1e-12 relative | asserted (`turb::thermal`) |
| Limits | T⁺ → Pr y⁺ at the wall, log slope 2.12, Pr_t → 2 Pr_t∞ at ν_t → 0 | 1e-3 | asserted |
| Heated pipe, Re = 10⁵, Pr = 0.71, k-ε scalable, 36 cells | Nu against Gnielinski (1976), 180.2 | ±10 % | +1.6 % (Kader), +0.8 % (Kays–Crawford), +5.0 % (Jayatilleke); `none` −39 % (`heated_pipe_meets_gnielinski`) |
| Isothermal flat plate, Re = 5·10⁶ per unit length, 0.4 ≤ x ≤ 1 | St against 0.0287 Re_x^−0.2 Pr^−0.4 | ±10 % | +8.3 % at 0.55, +9.0 % at 0.97; +14 % at 0.29 (`heated_plate_meets_the_stanton_correlation`; see below) |
| Input | `turbulent_prandtl` ≤ 0, unknown wall function | refused | asserted |
| Heated channel, Re_τ = 395, Pr = 0.71, 12 cells, y⁺ 32 (`channel_heated_keps.toml`) | Nu on 2h against 36.39 of the DNS of Abe, Kawamura & Matsuo (JAXA DNS database, Re_b = 13 849) | ±7 % | −1.8 % (constant Pr_t), −2.8 % (Kays–Crawford) |
| Heated channel, Re_τ 180 … 1020, Pr 0.71 and 0.025 (table below) | Nu on 2h against the DNS | Pr 0.71: ±7 % both models; Pr 0.025: ±8 % Kays–Crawford | asserted for all eight points (`heated_channel_meets_the_dns_nusselt_number`) |
| Heated channel, Re_τ = 395, Pr = 0.71, 48 cells, y⁺ 7.8 | thermal law at the scalable floor | ±13 % | +10.1 % (+28 % with T⁺ at the true y⁺) (`the_thermal_law_follows_the_scalable_floor`) |

**The channel at all DNS points** (pull request 87). The case of
`channel_heated_keps.toml` with ν = 1/Re_τ, Γ = ν/Pr, the mass flow of
each DNS bulk Reynolds number (the Re_b of the file header of each
thermal case) and N uniform cells across 2h, chosen for a first-cell
y⁺ ≈ 30, the claimed band of the scalable wall functions. Reference: the
Nu of each file header, except Re_τ = 1020, Pr = 0.025 (header 278.71,
see below), where the Nu of the database's own profiles is used (9.81).
y⁺ and Re_τ are those of the run (u_τ from the driving force; standard
k-ε predicts a friction 3–7 % low at the fixed mass flow).

| Re_τ | Pr | Re_b | N | y⁺ | Nu DNS | Nu, Pr_t 0.85 | Nu, Kays–Crawford |
|---|---|---|---|---|---|---|---|
| 180 | 0.71 | 5 546 | 6 | 28.1 | 18.25 | 18.13 (−0.7 %) | 17.93 (−1.8 %) |
| 395 | 0.71 | 13 849 | 12 | 31.5 | 36.39 | 35.73 (−1.8 %) | 35.38 (−2.8 %) |
| 640 | 0.71 | 24 095 | 20 | 30.9 | 55.57 | 54.49 (−1.9 %) | 53.97 (−2.9 %) |
| 1020 | 0.71 | 41 441 | 32 | 31.3 | 83.00 | 82.88 (−0.1 %) | 82.16 (−1.0 %) |
| 180 | 0.025 | 5 731 | 6 | 28.9 | 5.48 | 6.16 (+12.4 %) | 5.80 (+5.9 %) |
| 395 | 0.025 | 14 147 | 12 | 32.1 | 6.60 | 7.88 (+19.4 %) | 7.07 (+7.2 %) |
| 640 | 0.025 | 24 428 | 20 | 31.3 | 7.89 | 9.57 (+21.3 %) | 8.42 (+6.7 %) |
| 1020 | 0.025 | 41 441 | 32 | 31.3 | 9.81 | 11.98 (+22.1 %) | 10.49 (+6.9 %) |

At Pr 0.71 both Pr_t models are within 2.9 % at every Re_τ. At Pr 0.025
(a liquid metal) conduction dominates: Pe_t = (ν_t/ν) Pr stays below
about 2, where the real Pr_t rises well above 0.85. The constant Pr_t
then over-predicts the turbulent heat flux and Nu by 12–22 %, growing
with Re_τ; Kays & Crawford, whose Pr_t tends to 2 Pr_t∞ at small Pe_t,
brings it to +6 … +7 %. The rest is not the mesh: at Re_τ 395, Pr 0.025,
Nu is 6.94, 7.07, 7.13 and 7.08 on 8, 12, 24 and 48 cells. At Pr 0.71
the k-ε model without damping sets the limit below y⁺ 30: Nu 35.52,
35.62, 38.21, 40.20 on 8, 12, 24, 48 cells (y⁺ 47, 32, 16, 8; with the
velocity scale u_τ of 6f, see below). Meshes
finer than y⁺ ≈ 30 near heated walls take ζ-f or EB with the compound
treatment, not the k-ε family.

**The velocity scale of the thermal law** (6g, pull request 90). 6f took
h_w = ρ u_τ/T⁺(y⁺) with u_τ = (τ_w/ρ)^½ for every model. Specification
5.2 asks for u* = C_μ^¼ k_P^½ for the k-ε family and ζ-f/EB (the blended
u_τ for SST). The two agree in an equilibrium log layer, but u_τ vanishes
with the wall shear: at the stagnation point of the impinging jet the
k-ε Nu₀ was 37, against 156 with u*. The thermal law now takes u* and
y* = u* y_P/ν for these models (the wall kernel writes y* next to y⁺).
The numbers above are those of u*: the channel moves by −0.1 to −1.2
points (Pr 0.71; Pr 0.025 unchanged, conduction dominates), the pipe by
+0.2 to +0.3 points, the fine channel of 48 cells from +10.5 to +10.1 %.
On the flat plate St follows k: the k-ε start of the turbulent layer at
x ≈ 0.16 leaves k above its equilibrium for a while, so St is +43 % at
x = 0.17, +14 % at 0.29 and +8 … +9 % from x = 0.4 on (with u_τ +1.7 %
at 0.29, +7.5 % at 0.97). The test now covers 0.4 ≤ x ≤ 1. ζ-f₀ on the
jet (y⁺ below 1) is unchanged to 0.1 %: below y⁺ ≈ 1 both scales give
the molecular conduction.

The Re_τ = 180 channel on 6 cells is exactly one-dimensional; its pressure
residual (about 10⁻⁴⁰) lies below the float32 range of the multigrid
hierarchy, so the preconditioned residual vanished and CG reported a
breakdown. CG now stops there, as BiCGStab did already.

The DNS data belong to JAXA (all rights reserved; citation of Abe et al.
requested), so the repository holds only the published Nusselt number.
`examples/studies/channel_heat_dns.py` downloads the tables to a
git-ignored folder, recomputes U_b⁺, Θ_b⁺ and Nu for the eight cases
(Re_τ 180, 395, 640, 1020 at Pr 0.71 and 0.025; the integrated Nu
agrees with the file headers to 0.4 % except `ch1020pr0025.dat`, whose
header repeats the Re_b of the Pr = 0.71 case and gives Nu = 278.71
where the profiles give 9.81) and compares our run.

Not in this part: the heated pipe with the other models, and the plate for x below 0.2, where
the k-ε start of the turbulent layer, not the thermal law, sets St
(97 % low before it, 43 % high just after it).

## Elliptic-blending k-ε (D133)

Specification: [requirements stage 6](../requirements/stage6-turbulence.md)
section 3.4 (D133), `foehncraft-kernels::turb::eb`.

**The model.** Billard & Laurence (2012), the BL-v²/k model, checked
term by term against the paper itself (Appendix A, eqs. A.1 to A.7,
Table A.10 and eq. 18): ν/2 in the molecular diffusion of k, ε and φ;
the E term C_ε3 (1 − α)³ (k/ε) 2ν ν_t (∂²U_i/∂x_j∂x_l)² in the k
equation; C*_ε2 = C_ε2 + α³(C_ε4 − C_ε2) tanh(|∂_j(ν_t/σ_k ∂_j k)/ε|^{3/2})
in the ε equation; the φ source (1 − α³) f_w + α³ f_h − φP_k/k +
(2/k)(ν_t/σ_k)∇φ·∇k with f_w = −εφ/(2k) and f_h = −(C₁ − 1 +
C₂P_k/ε)(φ − 2/3)/T; α − L²∇²α = 1; T = (k²/ε² + C_T²ν/ε)^½,
L = C_L (k³/ε² + C_η²ν^{3/2}/ε^{1/2})^½, ν_t = C_μ φ k min(T, T_lim),
T_lim = 0.6/(√6 C_μ φ (S_ij S_ij)^½); wall values k = φ = α = 0,
ε_w = νk/y². Constants: C_μ 0.22, C_ε1 1.44, C_ε2 1.83, C_ε3 2.3,
C_ε4 0.4, σ_k 1, σ_ε 1.5, σ_φ 1, C_T 4, C_L 0.164, C_η 75, C₁ 1.7,
C₂ 0.9. A published restatement of the model (the theory guide of an
open-source code) prints C*_ε2 with C_ε2 − C_ε4; the paper and this
code have C_ε4 − C_ε2, so C*_ε2 falls from 1.83 towards 0.4.

**The implementation.** The kernel of k-ε-ζ-f (`turb::zeta_f`) with
the form `FORM_EB`: φ takes the slot of ζ, α the slot of f, so both use
the transport and elliptic solvers, the multigrid and the wall
treatments of ζ-f. The E term needs (∂²U_i/∂x_j∂x_l)²: the least-squares
gradient of each of the nine velocity-gradient components, squared and
summed (`squares_add`). The switch of C*_ε2 needs ∇·(ν_t/σ_k ∇k): the
gradient of the vector (ν_t/σ_k)∇k (`scale_vector`), its trace. Both use
zero-gradient boundary values; α³ removes them next to walls. α is
clipped to [0, 1] and starts at 1, also on a restart from a
two-equation model.

**Finding: the form of the α equation.** The paper writes
α − L²∇²α = 1 with L² outside the Laplacian, and that form is solved here
(∇²α + (1 − α)/L² = 0: diffusivity 1, source (1 − α)/L²). The
conservative form −∇·(L²∇α) + α = 1 of the f equation of ζ-f, used
first, adds ∇L²·∇α, which near the wall, where L changes fast, raises α
and switches the E term off too early: U_b⁺ at Re_τ 395 came out 3.5 %
low (C_f +7 %), at 590 4.5 % and at 2 000 5.6 % low, independent of
the mesh (16.53, 16.86, 16.93 on 48, 96, 192 cells at Re_τ 395). An
independent one-dimensional solver of the model reproduced both forms
(published form: C_f at 97 % of the DNS at Re_τ 395).

### Results

| Test | Check | Band | Result |
|---|---|---|---|
| Channel at Re_τ = 395 (`channel_eb.toml`), 96 / 192 cells, first-cell y⁺ 0.28 / 0.14 | C_f/C_f,DNS (the driving force) against the DNS | U_b⁺ ±3 % | 99.0 / 98.8 %, U_b⁺ +0.5 / +0.6 %, converged (`eb_channel_against_the_dns`) |
| Channel at Re_τ = 590, 128 cells, y⁺ 0.25 | C_f/C_f,DNS against the DNS | U_b⁺ ±3 % | 99.1 %, U_b⁺ +0.5 % |
| Channel at Re_τ = 180, 96 cells, y⁺ 0.13, start I = 0.15, ν_t/ν = 20 | C_f/C_f,DNS against the DNS | U_b⁺ ±3 % | 97.8 %, U_b⁺ +1.1 %; from I = 0.05, ν_t/ν = 100 laminar (F = 0.260, the laminar 0.258) |
| The channel at Re_τ 395 with `algorithm = "coupled"` | driving force | – | 0.98993 against 0.98982 (SIMPLEC), converged in 315 iterations |
| TMR plate (`tmr_plate_eb.toml`), turbulent start, y⁺ 1 | C_f(0.97) against White, 2.875e-3 | ±3 % | 2.889e-3 (+0.5 %); k and φ held by the floor near the leading edge (TURB002) |
| The plate from the TMR free stream | C_f(0.97) | ±3 % | turbulent from x ≈ 0.25, +2.5 % |
| Backward-facing step (`backstep_eb.toml`), 8 200 cells, SIMPLEC | x_r/H against 6.26 | ±10 % | 7.59 (+21 %; finding), converged to 10⁻⁶ in 2 480 iterations |
| The step on 32 800 cells (each interval halved), SIMPLEC | x_r/H | ±10 % | 7.68 (+23 %), converged in 7 330 iterations; 1.2 % from the coarse mesh |
| The step with `algorithm = "coupled"`, 8 200 / 32 800 cells | x_r/H | ±10 % | 7.74 to 7.77 / 7.84, not converged (finding) |
| The step, the plate and the channel with C*_ε2 = C_ε2 (`c_epsilon4 = 1.83`) | x_r/H, C_f(0.97), C_f/C_f,DNS | – | 6.24 (−0.3 %), +4.9 %, 100.9 % |
| Homogeneous decay (box at rest, ν = 10⁻¹², BDF2 with 10, 20, 40 steps) | k, ε with C_ε2 = 1.83 and φ → 2/3 with C₁ = 1.7 against the exact solution; α = 1 | order 1.8 to 2.2 | asserted (`eb_homogeneous_decay_is_second_order_in_time`) |
| k-ε restart continued with the elliptic blending | φ in (0, 2], α in [0, 1] | – | asserted |
| The plate on 1 and 4 threads (20 480 cells) | residuals, wall shear after 8 iterations | bit-identical | asserted |
| Bad input | `compound-2023`, `scalable`, `automatic`, a production limiter; `sigma_zeta`; `c_epsilon3 = 0` | refused | CASE005, CASE002, refused |
| Kernels | sources against the formulas of the paper and their signs, the Kolmogorov bounds of T and L at a wall, dual numbers | – | kernel tests of `turb::eb` |
| ζ-f and ζ-f₀ after the change | channel driving force | bit-identical | 0.9624991 before and after |

**Findings.** (1) The backward-facing step: the variable C*_ε2 lowers
C_ε2 towards 0.4 where the turbulent transport of k is large against ε,
ε grows, ν_t drops and the bubble reaches 7.59 H (+21 %; ζ-f of 2004:
7.82, ζ-f₀: 6.38). With a constant C_ε2 the step is right (6.24) and the
plate is 4.9 % high, so the switch helps the attached boundary layer and
hurts the separated flow. The model stays as
published; `c_epsilon4 = 1.83` is the option. (2) With the coupled
algorithm the step does not converge: the residual of ε stays at 10⁻³
and x_r moves between 7.74 and 7.77 (SIMPLEC converges to 7.59); the
channel converges with both. (3) At Re_τ 180 the channel needs a strong
start. (4) As with ζ-f, the positivity floor holds k and φ in the
laminar part of the plate near the leading edge, so their residuals
stay at 5·10⁻⁶ and 1.4·10⁻⁵. (5) The draft of specification 3.4 is
corrected after the paper (ν/2 in the diffusion, the forms of T and L,
ε_w = νk/y², the E term in the k equation, the α equation with L²
outside the Laplacian). (6) Not in this part: the export of the wall
normal ∇α/|∇α| (specification 2.3; with GGDH in stage 7), the near-wall
profiles against the DNS, the step against the experiment's C_p and C_f
and a GCI (stage 6 validation).

## Backward-facing step against the experiment (part 6g, 9.5)

Specification: [requirements stage 6](../requirements/stage6-turbulence.md)
section 9.5; `examples/cases/flow/backstep_validation.toml`,
`examples/studies/backstep_validation.py`, `tests/wall_crossing.rs`.

**The case.** The step of Driver & Seegmiller (1985), Re_H = 36 000, in
the layout of NASA's Turbulence Modeling Resource (TMR): the walls start
at x = −110 H (the TMR has a further 20 H of symmetry plane ahead of
them; here the uniform inflow, I = 1 %, ν_t/ν = 10, is at the leading
edge), the outlet is at x = 50 H, and the top wall at y = 9 H is a
no-slip wall like the lower one (expansion ratio 1.125). The earlier
step cases (parts 6d2 and D133) have a slip top at 9 H and walls from
−85 H. The built-in mesher cuts the step out of the box: 136 × 96 cells
less the 56 × 36 under the step, 11 040 cells, graded towards the walls,
the step corner and the reattachment region (first-cell y⁺ about 1.5);
`verification.levels = [0, 1, 2]` halves every interval once and twice
(44 160 and 176 640 cells). The TMR's own grids could not be
downloaded here, so the family is not theirs.
The inflow boundary layer at x = −4 H: C_f = 2.95·10⁻³ with ζ-f₀ on the
coarse mesh (the experiment 2.88·10⁻³).

**The quantities.** The reattachment length is a wall crossing
(`[[output.wall_crossings]]`, new in this part): the first point after
`from` + `skip` along `direction` where the wall shear along `direction`
turns from negative to positive (`reattachment`) or the other way
(`separation`), linear between the two face centroids, reported as the
distance from `from`. It is in the report of every solve and, in a level
sweep, a monitored quantity of D67 (`wall_crossings.xr.position`) with
its observed order and GCI. The study evaluates the finest mesh against
the data of the TMR (CC0, `examples/studies/backstep/data`): the
reference speed is u at x/H = −4, y/H = 5 (1.05 m/s); C_p is shifted to
zero at x/H = 37.5, the last station of the shifted data; C_f =
τ_x/(½ρU_ref²). The RMS differences are taken at the experiment's
stations: 33 of C_p on the lower wall (x/H from −8.5 to 37.5), 20 of C_f
(−3.96 to 36) and 28 of u/U_ref per profile. The models' shear stress
−u'v' = ν_t ∂u/∂y on the lines uses the mean ν_t and the least-squares
slope of u over ±0.05 H (the lines sample cell values).

**The models and their algorithms.** k-ω SST (`automatic`, the ambient
values of the inflow), k-ε and realizable k-ε (`scalable`) and
k-ε-ζ-f₀ (`compound-2023`) run with the coupled algorithm, at most 800
iterations on each mesh. k-ε-ζ-f and the elliptic-blending k-ε
(`compound`) do not converge with it (finding (2) of D133); they run
with SIMPLEC, relaxation 0.5 and 0.8 from iteration 100, and the
pressure CG limit at 200 iterations: with the default 100 both diverged
on the middle mesh within six iterations (RUN014, DIV006).

### Results

*Since PR #141* ζ-f and ζ-f₀ take the wall sink of k (specification
3.3.2) and ζ-f₀ C_ε1 = 1.4 (1 + 0.009/ζ). The step of ζ-f₀ on 8 200
cells then gives x_r/H = 6.44 (+3 %) with SIMPLEC and with the coupled
algorithm (`backstep_zeta_f0_coupled.toml` converged in 165 iterations,
x_r within 1 % after about 50). The tables below are those of 6g, with
the sink k/T_K of 6d2 for ζ-f₀ and −ε for ζ-f.

| Test | Check | Band | Result |
|---|---|---|---|
| Step (`backstep_validation.toml`), k-ω SST, coupled, 11 040 / 44 160 / 176 640 cells | x_r/H on the finest mesh against 6.26 | ±10 % | 6.747 (+7.8 %) |
| | x_r/H on the three meshes, GCI_fine, p | GCI ≤ 2 %, p in [1, 3] | 6.793, 6.747, 6.747: oscillatory, the fine pair 0.008 % apart, GCI below 0.01 %, p not defined (**finding**) |
| | iterations | – | 800 (residual 1.8·10⁻⁴), 621 (converged), 800 (6·10⁻⁶) |
| The step, standard k-ε, `scalable`, coupled | x_r/H against 5.345 (Mizukami, Georgiadis & Cannon 1993) | ±5 % | 5.372 (+0.5 %; −14.2 % from the experiment) |
| | x_r/H on the three meshes, GCI_fine | ≤ 2 % | 5.237, 5.321, 5.372: p 0.69, GCI 1.94 %, converged on each mesh |
| The step, realizable k-ε, `scalable`, coupled | x_r/H | reported | 6.274, 6.337, 6.407 (+2.4 %): the steps grow (0.063, 0.070), no order and no GCI (**finding**) |
| The step, k-ε-ζ-f₀, `compound-2023`, coupled | x_r/H against 6.26; GCI_fine, p | ±10 % (the band of ζ-f); GCI ≤ 2 %, p in [1, 3] | 6.329, 6.297, 6.296 (+0.6 % on the finest mesh): monotonic, GCI < 0.01 %, p 4.8, outside [1, 3] because the fine pair is only 0.02 % apart; iterations 800 (residual 5·10⁻⁵), 238 (converged), 800 (1.9·10⁻⁵) |
| The step, k-ε-ζ-f, `compound`, SIMPLEC, 11 040 / 44 160 cells | x_r/H against 6.26 | ±10 % | 7.686, 7.741 (+23.7 %; **finding**); k at its floor in 44 cells, its residual stays at 1.4·10⁻⁶ |
| The step, elliptic-blending k-ε, `compound`, SIMPLEC, 11 040 / 44 160 cells | x_r/H against 6.26 | ±10 % | 7.446, 7.556 (+20.7 %; **finding**), converged in 2 754 and 8 684 iterations |
| C_p on the lower wall, 33 stations, finest mesh (ζ-f and EB: middle mesh) | RMS difference from the experiment | ≤ 0.02 for SST, ζ-f, EB | SST 0.023, k-ε 0.018, realizable 0.026, ζ-f₀ 0.017, ζ-f 0.037, EB 0.034 (**finding** for SST, ζ-f, EB) |
| C_p on the upper wall, 19 stations | RMS difference | reported | SST 0.017, k-ε 0.012, realizable 0.016, ζ-f₀ 0.012, ζ-f 0.015, EB 0.016 |
| C_f on the lower wall, 20 stations | RMS difference ×10³ | reported | SST 0.375, k-ε 0.373, realizable 0.153, ζ-f₀ 0.463, ζ-f 0.593, EB 0.410 |
| C_f at x/H = −3.96 (the inflow boundary layer) | against the experiment's 2.88·10⁻³ | reported | SST 2.71, k-ε 2.93, realizable 2.67, ζ-f₀ 3.02, ζ-f 3.13, EB 3.04 (·10⁻³) |
| The layout: slip top at 9 H and walls from −85 H (6d2) against the TMR layout, ζ-f₀, coarse meshes | C_p RMS; x_r/H | – | 0.033 against 0.018; 6.38 against 6.33 |
| Wall crossing on a laminar step (`tests/wall_crossing.rs`), Re 133 | position against the sign change of the wall file, linear between the face centroids; searched upstream from the outlet; within `width`; after `skip` | 10⁻¹² | asserted |
| | a kind the shear does not show; `skip` past the reattachment | no position | asserted |
| | three levels | GCI < 5 % | asserted |
| Bad input | a zero direction, `width = 0`, `skip < 0`, an unknown kind | refused | CASE004 (three), the parse error naming the two kinds |

**Findings.** (1) The top wall. With the slip plane of the earlier
step cases the C_p of ζ-f₀ is 0.033 RMS from the experiment, with the
no-slip wall of the experiment's channel 0.018; the reattachment moves
by less than 1 %. (2) The corner vortices. On the finer meshes the
vortices in the corner behind the step change the sign of the wall
shear before x = 1.5 H (SST on the finest mesh: negative to 0.14 H,
positive to 1.44 H), and the first crossing after the step was theirs
(SST 0.10 and 0.14). `skip = 2.0` starts the search behind them; the
position is still measured from the step. (3) SST: the coarse mesh does
not converge with the coupled algorithm (its residual stalls at
1.8·10⁻⁴; x_r changed by 0.17 % from 400 to 1 200 iterations), the
middle mesh does, and the finest mesh reaches 6·10⁻⁶. The TMR reports
that SST does not converge to a steady state on its refined grid
either. The x_r of the two finer meshes agree to 0.008 %, so the
coarse mesh is outside the asymptotic range and the sweep is
oscillatory: no observed order. ζ-f₀ behaves the same way (residuals
5·10⁻⁵, converged, 1.9·10⁻⁵ on the three meshes); its sweep is
monotonic, but the fine pair is only 0.02 % apart, so p = 4.8 says no
more than that the two finer meshes are mesh-independent within the
iteration error. (4) The k-ε models run with scalable
wall functions on a mesh made for resolved walls: the first-cell y⁺ of
the lower wall is 1.5, 0.7 and 0.3 on the three meshes, and the solver
says that the wall functions are claimed for y⁺ ≥ 30. Halving the cells
moves the first cell deeper into the sublayer, which changes what the
wall function does as well as the discretisation. That may be why the
order of k-ε is 0.69 and the steps of realizable k-ε grow; a sweep that
keeps the first cell would show it (not done). (5) ζ-f and EB reattach
21 to 24 % late, as on the mesh of 6d2 (ζ-f 7.82, EB 7.59 there), and
their C_p is 0.034 to 0.037 RMS from the experiment; the causes are those of
6d2 (the realizability bound of T) and D133 (the variable C*_ε2).
Neither converges with the coupled algorithm, and with SIMPLEC the
pressure CG needs 200 iterations on the finer meshes. Their finest mesh
needs an estimated 30 000 SIMPLEC iterations, hours on four cores, and
follows in the next 6g PR. (6) SST misses the C_p band by 0.003 on the finest
mesh (0.019 on the coarse one); its recovery after reattachment is
slower than the experiment's (C_f 1.6·10⁻³ at x/H = 30 against
2.1·10⁻³). (7) Run time, one thread per run with four to six runs
sharing four cores: k-ε 23 min, realizable 34 min, SST 36 min for the
three meshes; ζ-f and EB on the two coarser meshes 72 and 39 min.

## Impinging jet at H/D = 2 (part 6g, 9.6)

Specification: [requirements stage
6](../requirements/stage6-turbulence.md) section 9.6;
`examples/cases/flow/jet_pipe_profile.toml`,
`examples/cases/flow/impinging_jet.toml`,
`examples/studies/impinging_jet.py`, `tests/impinging_jet.rs`,
`tests/inlet_table.rs`, `tests/symmetry_stress.rs`.

**The case.** A round jet from a long pipe of diameter D hits a heated
plate at H = 2 D, Re = U D/ν = 23 750: the heat transfer of Baughn &
Shimizu (1989) and Baughn, Hechanova & Yan (1991), the velocity of
Cooper, Jackson, Launder & Liao (1993, Re = 23 000). Axisymmetric on a
0.5° wedge about the jet axis x: the plate is x = 0, the pipe exit the
plane x = 2 D for r < D/2, the rest of that plane an inlet at the total
pressure 0 for the entrained fluid, and the wall jet leaves at r = 6 D
at the static pressure 0. The built-in mesher gives 50 × 52 intervals,
the first cell at the plate 5·10⁻⁴ D (first-cell y⁺ below 1), clustered
at the pipe lip; the case refines once (10 400 cells), the study sweeps
`verification.levels = [0, 1, 2]` (2 600, 10 400 and 41 600 cells).

**The inlet.** The jet's turbulence comes from the pipe. A periodic
wedge pipe at the same Re (`jet_pipe_profile.toml`, a module of 0.1 D
with 48 cells across the radius, driven by the mass flow of the bulk
speed) samples the radius with each model; `foehncraft.profiles` writes
that line as a table with one row per cell and a wall row, and the
header names the columns: `r, u_n, k, epsilon, zeta` (`omega` for SST,
`zeta` is φ for EB). The turbulence columns set the inlet face by face
(specification 7.1, new in this part).

**The quantities.** A passive temperature θ (Pr = 0.71) with a uniform
flux q through the plate and θ = 0 at both inlets gives Nu(r) = q D/(Γ
θ_w(r)), Γ = ν/Pr, from `wall_xmin.csv`. Nu₀ is the value of the wall
probe `stag` at the axis, that is of the face next to it (r/D = 0.06,
0.03 and 0.015 on the three meshes); the probe `r01` gives Nu at r/D =
0.1. Wall probes now report the wall value and the flux of every scalar,
and these are D67 quantities of a level sweep
(`wall_probes.stag.theta_wall`), so both get an observed order and a GCI
(that of Nu is the same to first order). The study reports both: before
the fix of finding 2 the faces next to the axis were low with every
turbulence model. The excess k of the stagnation-point
anomaly is the largest k on the axis within x/D ≤ 0.25 of the plate
above that of k-ε-ζ-f₀ on the same mesh. The secondary maximum is the
largest local maximum of Nu for 0.8 ≤ r/D ≤ 4 at least 1 % above the
lowest Nu before it.

**The models and their algorithms.** k-ε-ζ-f₀ (`compound-2023`), k-ω SST
(`automatic`, Kato–Launder), k-ε and realizable k-ε (`scalable`) run
with the coupled algorithm, at most 4 000 iterations; the realizable k-ε
with the tolerance 10⁻⁴ (its residuals stall, finding 5). k-ε-ζ-f and
the elliptic-blending k-ε (`compound`) run with SIMPLEC, relaxation 0.5
and 0.8 from iteration 100 (EB 0.3, with 0.5 its coarse mesh diverges at
iteration 10), the pressure CG limit at 200 iterations; EB with the
tolerance 2·10⁻⁵ (finding 5). The limiter cases run on 10 400 cells.

**The experiment.** The data are those of the ERCOFTAC classic database,
case 025, files `ij2lr-*` (H/D = 2, Re = 23 000), uploaded by the
maintainer: Nu/Re^0.7 along the plate of Baughn et al. (labelled 1992 in
the file) and the single-wire mean velocity U/U_bulk of Cooper et al. at
r/D = 0.5, 1, 2, 2.5 and 3. The study takes Nu = (Nu/Re^0.7)·23 750^0.7
(factor 1 156.0, the source's own normalisation), so the experiment has
Nu₀ = 141.4, a minimum of 89.6 at r/D = 1.26 and a secondary maximum of
108.8 at r/D = 1.94 (the wiki of the database, cfd.mace.manchester.ac.uk/ercoftac,
is licensed CC BY-NC-SA 4.0 except where otherwise noted). That licence
is non-commercial and share-alike, which does not fit this MIT repository, so
the files are not in it (`.gitignore` keeps them out; the README of the
data folder says how to get them and `ercoftac_to_csv.py` makes the CSV
files the study reads), the charts of the gallery show the models only,
and this page quotes a few scalar values with their source.

### Results

*Since PR #141* (the wall sink of k for both ζ-f models, specification
3.3.2, and C_ε1 = 1.4 (1 + 0.009/ζ) for ζ-f₀) the jet gives, on 2 600 /
10 400 / 41 600 cells:

| Model | Nu₀ | p, GCI_fine | Nu(r/D = 0.1) | secondary maximum r/D, Nu | iterations |
|---|---|---|---|---|---|
| k-ε-ζ-f₀, coupled | 198.3, 195.9, 195.6 (+38 % from 141.4) | 3.49, 0.01 % | 196.0, 194.0, 194.2 | 2.03, 141.6 | 557, 884, 2 969 |
| k-ε-ζ-f, SIMPLEC | 176.8, 173.0, 172.5 (+22 %) | 2.86, 0.06 % | 175.9, 172.8, 172.7 | none | 518, 1 607, 3 383 |

ζ-f₀ now converges with the mesh (before 166.1, 156.9, 151.6, GCI
7.2 %), but Nu₀ rises from +7 % to +38 % of the experiment
(**finding**). The rise comes with the sink, not with C_ε1: the 6d2 sink
with the new C_ε1 gives 159 on the middle mesh, and a = 0.6 lowers the
new default only to 176. How the near-wall sink raises the k of the
stagnation region is open; a stagnation-point limiter for ζ-f₀ is the
follow-up. ζ-f
moves by less than 1 %. The u_r,max/U of ζ-f₀ at r/D 0.5, 1, 2, 2.5, 3:
0.715, 1.078, 0.742, 0.575, 0.462 (−2, −1, −6, −3, 0 % from the
experiment, as before). The tables below are those of 6g, before #141.

| Test | Check | Band | Result |
|---|---|---|---|
| Jet (`impinging_jet.toml`), k-ε-ζ-f₀, `compound-2023`, coupled, 2 600 / 10 400 / 41 600 cells | Nu₀ (axis face) on the three meshes, GCI_fine | GCI ≤ 3 % | 166.1, 156.9, 151.6: monotonic, p 0.67, GCI 7.2 % (**finding** 4; before the fix of **finding** 2 167.2, 155.4, 141.5, divergent) |
| | Nu at r/D = 0.1, GCI_fine, p | GCI ≤ 3 % | 164.2, 156.1, 151.7: monotonic, p 0.74, GCI 5.3 % (**finding** 4; extrapolated 145.5) |
| | secondary maximum; Nu at r/D 1, 2, 3, 4 | present | 128.7 at r/D 1.95 after the minimum 115.1 at 1.45 (10 400 cells: 130.1 at 1.97); 129.8, 128.4, 95.8, 69.2 |
| | iterations | – | 613, 862, 2 959 (converged; θ the last) |
| The jet, k-ε-ζ-f, `compound`, SIMPLEC | Nu₀ on the three meshes | GCI ≤ 3 % | 178.6, 174.6, 173.8: p 2.45, GCI 0.12 % (before the fix 179.3, 173.0, 167.2: no usable GCI) |
| | Nu at r/D = 0.1 | GCI ≤ 3 % | 177.7, 174.3, 173.9: p 3.17, GCI 0.03 % |
| | secondary maximum; Nu at r/D 1, 2, 3, 4 | present | 118.4 at r/D 2.27; 139.3, 117.4, 106.6, 81.1 |
| | iterations | – | 450, 1 254, 4 179 (converged) |
| The jet, elliptic-blending k-ε, `compound`, SIMPLEC (relaxation 0.3, tolerance 2·10⁻⁵) | Nu₀ on the three meshes | GCI ≤ 3 % | 236.7, 247.3, 249.2: monotonic, p 2.58, GCI 0.19 %; Nu peaks at the axis (249 at the axis face, 238 at r/D 0.055, 205 at 0.1, 147 at 0.5) |
| | Nu at r/D = 0.1 | GCI ≤ 3 % | 206.8, 205.6, 205.0: p 0.95, GCI 0.40 % |
| | secondary maximum; Nu at r/D 1, 2, 3, 4 | present | 93.4 at r/D 2.31 (weak); 113.9, 87.3, 79.1, 58.3 |
| | iterations | – | 456, 943, 2 123 (converged at 2·10⁻⁵; finding 5) |
| The jet, k-ω SST with `kato-launder`, `automatic`, coupled | Nu₀ on the three meshes | GCI ≤ 3 % | 151.9, 144.6, 141.2: p 1.01, GCI 2.9 % (before the fix 152.2, 143.0, 134.8: no usable GCI) |
| | Nu at r/D = 0.1 | GCI ≤ 3 % | 149.5, 142.4, 141.1: p 2.28, GCI 0.31 % |
| | secondary maximum; Nu at r/D 1, 2, 3, 4 | reported | 111.6 at r/D 1.84; 121.6, 110.2, 82.8, 61.7 |
| | iterations | – | 756, 978, 2 947 (converged) |
| The jet, standard k-ε (Durbin bound), `scalable`, coupled | Nu₀; Nu at r/D = 0.1 | reported | 155.0, 154.7, 155.0; 154.7, 155.2, 155.6: within 0.6 % on the three meshes (Nu₀ counts as divergent, steps of −0.3 and +0.3; at r/D = 0.1 p 0.52, GCI 0.67 %) |
| | secondary maximum; Nu at r/D 1, 2, 3, 4 | reported | 143.2 at r/D 1.55; 143.3, 128.8, 87.4, 64.1 |
| | iterations | – | 424, 1 063, 3 699 (converged) |
| The jet, realizable k-ε, `scalable`, coupled (tolerance 10⁻⁴) | Nu₀; Nu at r/D = 0.1 | reported | 208.3, 205.1, 207.2; 208.6, 205.4, 207.8 (oscillatory: the two coarser meshes stall); the anomaly: k on the axis 0.0552, 4.6 times that of ζ-f₀ |
| | secondary maximum; Nu at r/D 1, 2, 3, 4 | reported | 141.1 at r/D 1.65; 162.3, 135.3, 96.4, 71.0 |
| | iterations | – | 4 000 and 4 000 (stalled, **finding** 5), 783 (converged) |
| Stagnation-point anomaly, 10 400 cells, standard k-ε | Nu₀ with `none` | > +50 % of the experiment | 210.7: +49.0 % from the experiment's 141.4, **1 point short of the band** (+34 % above ζ-f₀'s 156.9); k_max on the axis 0.0787 against 0.0120 of ζ-f₀ |
| | excess k removed by `durbin` (default), `kato-launder` | ≥ 70 % | 71 % (Nu₀ 154.7), 111 % (Nu₀ 74.5: below ζ-f₀, k on the axis 0.0045) |
| Stagnation-point anomaly, 10 400 cells, SST | Nu₀ with `none` | > +50 % | 160.5 (+13.5 % from the experiment, +2 % above ζ-f₀); k_max on the axis 0.0104, below the 0.0120 of ζ-f₀: no excess (**finding** 3) |
| | excess k removed by `menter`, `kato-launder`, `durbin` | ≥ 70 % | `menter` and `durbin` 160.5, the same field (k 0.0104); `kato-launder` 144.6 (k 0.0038); no excess to remove |
| Nu₀ against the experiment (141.4), 41 600 cells | Nu₀ of the axis face | ζ-f, EB ±15 %; SST with a limiter ±20 % (ζ-f₀ as ζ-f) | SST −0.2 % and ζ-f₀ +7.2 % inside; ζ-f +22.9 % and EB +76.2 % (+45 % at r/D 0.1) **outside** (**finding** 8); reported: k-ε +9.7 %, realizable k-ε +46.6 % |
| Nu(r/D), 0 ≤ r/D ≤ 4 | RMS of the relative difference from the experiment | ≤ 15 % | SST 13.7 % met; EB 21.2 %, ζ-f₀ 23.7 %, ζ-f 31.5 % **not met**; reported: k-ε 31.0 %, realizable k-ε 45.9 % |
| Secondary maximum | r/D of the maximum | within ±0.5 of the experiment's 1.94 | ζ-f 2.27, EB 2.31 and ζ-f₀ 1.95 inside; SST 1.84 |
| Wall jet, u_r,max/U at r/D = 0.5, 1, 2, 2.5, 3 | difference from the experiment's 0.732, 1.089, 0.789, 0.594, 0.463 | reported | ζ-f₀ −2, −1, −6, −3, 0 %; SST +4, +3, 0, +6, +11 %; ζ-f −3, 0, +10, +16, +18 %; EB +2, +1, +5, +14, +23 %; k-ε −4, −8, −19, −14, −9 %; realizable k-ε −3, −9, −14, −9, −4 % |
| Wall probe with scalars (`tests/thermal.rs`, heated plate) | θ_w at the probe against the wall faces; flux; GCI quantities | interpolated between the two nearest faces; `wall_probes.<p>.theta_wall`, `_flux` | asserted |
| Inlet table with turbulence columns (`tests/inlet_table.rs`, the ζ-f₀ plate) | a table of constant columns (`u_x`, `u_y`, `u_z`, `k`, `epsilon`, `zeta`, in any order) against the constant inlet; ten times its k; turbulence values (`turbulence_intensity`, `zeta`) next to such a table; an unknown column | bit for bit; another run; refused | asserted; BC029, BC019 |
| `foehncraft.profiles` (`python/tests/test_profiles.py`) | one row per cell, the wall row; a laminar line | as specified; refused | asserted |
| The coarse jet (`tests/impinging_jet.rs`, 2 600 cells) | ζ-f₀ converged, probes against the wall table, Nu₀, secondary maximum; k-ε Nu₀ against Nu(r/D = 0.5) | to 10⁻⁹ and 10⁻³; 185 to 205 (150 to 185 before #141); r/D 1.5 to 2.5, 3 % above the minimum; > 0.9 | 198.3, maximum 138.6 at 2.01 above 124.2 (before #141 166.1, 128.2 above 116.8); k-ε 155.4 against 164.2 |
| The stress on symmetry and wedge planes (`tests/symmetry_stress.rs`) | a linear flow with constant stress, the viscous operator of the momentum rows per cell, without and with the transposed part: a radial flow v = a r, u = −2a x on a 5° wedge of perturbed cells; a stagnation flow onto a symmetry plane | 10⁻¹⁰ of a face flux | 8·10⁻¹⁶ and 1·10⁻¹⁴ (without the fix of finding 2: 3.6 % and 65 %) |

**Findings.** (1) The thermal law of the k-ε family (6f) took u_τ, which
vanishes with the wall shear: the k-ε Nu₀ of this jet was 37. It now
takes u* = C_μ^¼ k_P^½ as section 5.2 asks (the k-ε Nu₀ is 157 on 10 400
cells; the effect on the 6f cases is in the [6f
section](#thermal-wall-function-and-turbulent-prandtl-number-part-6f)).
(2) The axis cells, and the stress on symmetry faces (fixed in PR #91).
In the first runs Nu of the faces next to the axis was low with every
turbulence model, and more so on finer meshes: ζ-f₀ 1 % below the
plateau on 10 400 cells, 7 % (first face) and 2 % (second face) on
41 600; ζ-f 1 % and 4 %; SST 4 % on 41 600. In the axis row v/r was up
to 10 % above that of the next two rows and ν_t 1 to 3 % below, while
the laminar jet (Re = 500, the same meshes) was flat to the axis within
0.1 %. The cause: with a turbulence model the momentum rows take the
transposed part ∇·[μ_eff (∇u)ᵀ] of the stress (section 2.4 of the
specification) on the interior faces, but the symmetry faces (symmetry
planes, slip walls, the planes of a wedge) did not have it. On such a
face u_n = 0, so the transposed traction μ_eff (∇u)ᵀ·S has only the
normal part μ_eff |S| ∂u_n/∂n, and the normal stress on the face is
2μ_eff ∂u_n/∂n. On the two planes of a wedge that part is half of the
hoop stress −2μ_eff v/r², largest in the cells next to the axis and
growing as the mesh is refined (v/r² ≈ (∂v/∂r)/r). Laminar runs leave
the transposed part out on every face, so their operator was consistent.
The coefficient of a symmetry face now doubles with the transposed part
(2μ_eff D, implicit as before); `tests/symmetry_stress.rs` checks that a
linear flow with constant stress balances in every cell, a radial flow
on a wedge and a stagnation flow on a symmetry plane (the residual was
3.6 % and 65 % of a face flux, now 10⁻¹⁵ and 10⁻¹⁴). With the fix Nu is
flat to the axis on every mesh: on 41 600 cells ζ-f₀ 151.6 at the face
next to the axis and 151.7 at r/D = 0.1 (before 141.5 and 152.4), ζ-f
173.8 and 173.9 (167.2 and 174.2), SST 141.2 and 141.1 (134.8 and
141.1); beyond r/D = 0.25 Nu moved by less than 0.2 % (k-ε 0.4 %,
realizable k-ε 0.7 %), the wall jet by less than 0.01 U. Nu₀ at the
axis now converges for ζ-f (GCI 0.12 %) and SST (2.9 %). The symmetry
planes of the other turbulent cases (the z planes of the 2D cases, the
upstream plane of the flat plate, the top of the step) carry little
∂u_n/∂n, and their tests did not move outside their bands.
(3) The
stagnation-point anomaly and SST. SST's own bound of the eddy viscosity,
ν_t = a₁k/(F₂S), holds along the whole jet axis near the plate (with
`none` on 10 400 cells ν_t ω/k equals a₁ω/S within 4 % for x/D ≤ 0.6,
within 1 % from x/D = 0.06 on), so the production is a₁ k S there: the
2003 limiter (P ≤ 10 β* k ω, the ratio is 2 to 6) and the Durbin bound
(k/(√6 S), looser than a₁ = 0.31) do not act, and `none`, `menter` and
`durbin` give the same field. Only Kato–Launder changes it. The anomaly
that 9.6 asks SST without a limiter to show is not there to remove. (4)
ζ-f₀ near the stagnation point does not reach its asymptotic range on
these meshes: Nu at r/D = 0.1 falls by 8.1 and 4.4 (p = 0.74, GCI 5.3 %,
above the band of 3 %), at the axis by 9.2 and 5.3 (p 0.67, GCI 7.2 %);
ζ-f, EB and SST converge at r/D = 0.1 with p of 0.95 to 3.2 and GCIs
below 0.4 %, the standard k-ε changes by less than 0.6 % over the three
meshes, the realizable k-ε by 1.6 % (its coarser meshes stall). (5) Steady states that stall. EB holds k at its
floor in a few cells: the residuals of k and φ stay at 9.7·10⁻⁶ and
1.1·10⁻⁵ on the coarse mesh, unchanged to eight digits from iteration 2
000 to 20 000, every other residual below 10⁻⁸; the study takes 2·10⁻⁵
for EB. The realizable k-ε stalls on the two coarser meshes with the
coupled algorithm: on the coarse mesh its residuals stay at 1.2·10⁻⁴
(momentum) and 3.6·10⁻⁴ (k) from iteration 1 000 to 4 000, changing by
less than 1 % (before the fix of finding 2 unchanged to eight digits),
so the iteration has reached a state it no longer changes but that does
not satisfy the equations to the tolerance; the finest mesh converges
to 10⁻⁴ in 783 iterations. Its variable C_μ does not
prevent the anomaly here (Nu₀ 207, the axis k 4.6 times that of ζ-f₀).
(6) The passive scalar sets the iteration count. On 41 600 cells the
flow and turbulence residuals of the coupled runs are below the
tolerance 10⁻⁶ after 472 (ζ-f₀), 546 (k-ε) and 720 (SST) iterations, θ
only after 2 959, 3 699 and 2 947: four fifths of the iterations of the
finest mesh go to the passive scalar (one solve per outer iteration with
the relaxation 0.9 and the limited linear-upwind correction, each solve
two to nine BiCGStab iterations; the laminar jet still had 3·10⁻⁵ after
3 000). `solver.scalar.sweeps` and `relaxation` are the case's means; a
faster scalar path is solver work for later. (7) Run time, two runs of
two threads sharing four cores, the three meshes in one run: EB 8 min,
ζ-f 13, SST 11, k-ε 12, realizable k-ε 15, ζ-f₀ 16; the eight limiter
cases on 10 400 cells about a minute each; the study 42 min with the
pipes. (8) Against the experiment. At r/D = 1 the models are 18 to 45 % above
the experiment (96.4): EB 113.9, SST 121.6, ζ-f₀ 129.8, ζ-f 139.3. The
dip of the experiment (89.6 at r/D 1.26) comes earlier and is deeper
than that of ζ-f₀ (115.1 at 1.45) and SST (109.3 at 1.48). From r/D = 2
on SST is within 5 % (110.2, 82.8, 61.7 against 107.7, 79.9, 59.2),
ζ-f₀ 17 to 20 % above, ζ-f 9, 33 and 37 % above, and EB 19 % below at
r/D = 2 and within 2 % at 3 and 4. So the Nu₀ band of ζ-f (+22.9 %) and
EB (+76 %) is not met, nor the RMS band of any model but SST. EB peaks
at the axis (249 at the face next to it, 205 at r/D 0.1, 181 at 0.15;
the other models are flat there), but even at r/D = 0.1 it is +45 %.
The experiment's Nu is scaled from Re 23 000 to 23 750 with Re^0.7
(2.3 %); the models were not run at 23 000. The anomaly test of 9.6
(more than +50 %) is missed by 1 point on 10 400 cells (+49.0 %,
standard k-ε with `none`); the finest mesh was not run for the limiter
cases.

## Moving reference frame (part 6g, 9.7)

`[[mrf]]` (D55): a cell zone that turns steadily about an axis, the
absolute velocity as the solved variable, the flux relative to the frame as
the convecting flux of momentum, turbulence variables and scalars, the
Coriolis source in the rotating cells ([specification 9.7](../requirements/stage6-turbulence.md)).
The case is laminar Taylor–Couette flow: a cylinder of radius
r₁ = 0.05 m turning at ω₁ = 10 rad/s inside a fixed cylinder of
r₂ = 0.10 m, ρ = 1 kg/m³, μ = 5·10⁻⁴ Pa s (Re = ω₁ r₁ (r₂ − r₁)/ν = 50),
2D on the built-in `[mesh.annulus]` (`taylor_couette_mrf.toml`), exact
solution u_θ = A r + B/r with A = −ω₁ r₁²/(r₂² − r₁²), B = ω₁ r₁² r₂²/(r₂² − r₁²)
and the pressure from dp/dr = ρ u_θ²/r. The inner half of the radial layers
(the zone `rotor`) turns with the inner cylinder; its interface at the
middle of the gap is a circle about the axis. Errors are relative to the
speed of the inner wall (0.5 m/s).

| Cells | L2 error | L∞ error | pressure L∞ (of the range) | Iterations (SIMPLEC) |
|---|---|---|---|---|
| 256 | 3.46·10⁻² | 1.01·10⁻¹ | 4.4·10⁻² | 125 |
| 1 024 | 6.94·10⁻³ | 2.95·10⁻² | 3.5·10⁻² | 366 |
| 4 096 | 1.31·10⁻³ | 7.63·10⁻³ | 2.3·10⁻² | 1 169 |

Observed order 2.32 and 2.41 (the band is 1.8 to 2.6). The coupled algorithm
gives the same errors (L2 6.94·10⁻³ on 1 024 cells, 582 outer iterations).

**Frame independence** (1 024 cells, `tests/mrf.rs`). The flow is steady in
every frame about the axis, so the angular velocity of the zone must not
matter beyond the discretisation error: no `[[mrf]]` 6.26·10⁻³, Ω = ω₁
6.94·10⁻³, Ω = −ω₁ 5.72·10⁻³, Ω = 2.5 ω₁ 5.77·10⁻³; a table with Ω = 0 gives
exactly the result without it (bitwise). The test requires all four
below 1.2 % and within a factor of two. A passive scalar (θ = 1 inside, 0
outside, exact ln(r₂/r)/ln(r₂/r₁)) has the L2 error 2.12·10⁻³ in the frames
tried, equal to 11 digits.

**Solid-body rotation** (both cylinders and the whole mesh at 10 rad/s, the
frame turning with the fluid): L2 error 1.05·10⁻² in the frame and 0.83·10⁻² in
the inertial frame. In the frame the relative flux is zero and only the
Coriolis source and the pressure gradient remain, but their balance on a
polygonal ring is only second order (the cell gradient of p ∝ r² against the
Coriolis force at the cell centre), so the frame does not resolve the case
better than the inertial one does.

Findings:

1. The closure of the relative flux is exact for surfaces of revolution (the
   test finds a worst cell below 10⁻¹² of ρ|V||S|) and fails at the 10⁻³
   level when the axis is moved off the centre of the zone, so MRF003 has a
   meaning; the scale of the check is ρ|V||S| summed over the faces of the
   cell, not the frame flux itself, which is rounding noise on a face that
   the frame velocity slides along.
2. The explicit Coriolis source needed no under-relaxation at Re = 50 and
   ω up to 25 rad/s on these meshes (convergence with SIMPLEC at every
   frame tried).
3. Not tested yet: a turbulent rotating zone, the channel in spanwise rotation of
   9.7 (a periodic case, with 6p), and the combination with the energy
   equation (stage 7).

## Sliding mesh with conformal stepping (part 6g, 9.7)

`[[sliding]]` (D56, variant A of D182): the rotor of a ring mesh turns by
whole sectors in every time step, Δt = k · 2π/(N |Ω|), so the interface
cells meet new neighbours without interpolation and the mesh stays
conformal ([specification 9.7](../requirements/stage6-turbulence.md)). The
case is the Taylor–Couette flow of the section above (r₁ = 0.05 m,
r₂ = 0.10 m, ω₁ = 10 rad/s, ρ = 1 kg/m³, μ = 5·10⁻⁴ Pa s, Re = 50), 2D on
`[mesh.annulus]` with the inner half of the radial layers as the zone
`rotor`, now **started from rest** and followed for 510 steps on the case
file's mesh (t = 10 s, two viscous times (r₂ − r₁)²/ν), `taylor_couette_sliding.toml`,
against the exact steady solution u_θ = A r + B/r. Errors are relative to the
speed of the inner wall. BDF2, one sector per step, each mesh with the
sector angle and the time step halved:

| Cells (sectors) | Δt (s) | L2 error | Observed order |
|---|---|---|---|
| 128 (16) | 0.0393 | 6.76·10⁻² | |
| 512 (32) | 0.0196 | 1.35·10⁻² | 2.32 |
| 2 048 (64) | 0.00982 | 2.77·10⁻³ | 2.28 |
| 8 192 (128) | 0.00491 | 6.15·10⁻⁴ | 2.17 |

The test (`tests/sliding.rs`, three meshes) requires orders between 1.9 and
2.7 and the errors below 2·10⁻² and 4·10⁻³ on the 512 and 2 048 cells. The
512-cell case, 510 steps, runs in 0.7 s.

**Two sectors per step** (Δt twice as long, same mesh): L2 error 2.7·10⁻²,
four times the error of the 1.35·10⁻² (second order in the time step); the
test requires below 8·10⁻².

**Solid-body rotation** (both cylinders and the whole mesh at 10 rad/s, the
fluid started from the exact solution, 512 cells): the L2 error is 1.73·10⁻²
at 0.6 s and 1.57·10⁻² at 3 s, it does not grow (the test requires below 3·10⁻² and the
long run below 1.5 times the short one). A sliding mesh in which the old
time levels are the only thing that moves with the rotor keeps this state;
the first design moved the fields, gradients and fluxes of the rotor cells
with the data and diverged within 130 steps.

**Passive scalar** (θ = 1 inside, 0 outside, conduction, diffusivity 10⁻³
m²/s, rotor turning through it): L2 error 7.3·10⁻³ against the exact
ln(r₂/r)/ln(r₂/r₁) at the end (the test requires below 1.5·10⁻²): the old
level of a moved cell is what that cell held.

Findings:

1. **Only the time term sees the motion.** All fields, gradients, the
   Rhie–Chow velocity and the extrapolations stay Eulerian on the fixed
   lattice; the old time levels of a rotor cell are gathered from the cell that
   stood at its place one and two sectors earlier (a permutation and its
   square), and the convecting flux is the relative one. Shifting
   gradients or the pressure gradient with the data, or recomputing gradients
   from shifted unrotated fields, is wrong (28 to 61 % error) or unstable in
   the solid-body test. The vectors of the gathered levels are not rotated:
   the pressure gradient turns them, and BDF2 leaves an O(Δθ⁴) error per step
   (BDF1 loses cos Δθ per step and a first-order error).
2. A fluid that stands in space and is convected through the rotor at a CFL
   number of one with the implicit upwind scheme is smeared (scalars); this
   is the usual behaviour of a sliding or ALE scheme and is accepted.
3. Not tested yet: restart of a sliding run, a rotor with blades, a
   sliding mesh with the energy equation (stage 7), turbulence models
   (refused, `SLD004`, until PR #100: see the section on turbulence below),
   and the time-dependent speed Ω(t) of D182.

## Sliding mesh by imprinting (part 6g, 9.7)

The same Taylor–Couette flow as above (r₁ = 0.05 m, r₂ = 0.10 m, ω₁ = 10 rad/s,
ρ = 1 kg/m³, μ = 5·10⁻⁴ Pa s, Re = 50, started from rest, t = 10 s), now with
a **pitch ratio that is not one and a free time step**: variant B of D182,
`[[sliding]] rotor_interface / stator_interface` on `[mesh.annulus] n_theta_stator`
(the rotor with 2n sectors and the radial layers inside r = 0.075 m, the
stator with 3n sectors, the interface at r = 0.075 m in two boundaries,
`taylor_couette_imprint.toml`). In every step the rotor nodes are turned to
the angle Ω t, every rotor face is clipped with every stator face in the
unrolled interface cylinder, and the intersection polygons are the internal
faces of the interface (`foehncraft_mesh::imprint`); the solver is built on
that mesh with the state, the old levels and the fluxes carried over, fully
implicit (`extrapolation = false`, four outer iterations). BDF2, Δt = 0.04 s
on the coarsest mesh and halved with the mesh:

| Cells (rotor + stator sectors) | Δt (s) | L2 error | Observed order |
|---|---|---|---|
| 160 (16 + 24) | 0.04 | 6.74·10⁻² | |
| 640 (32 + 48) | 0.02 | 1.29·10⁻² | 2.39 |
| 2 560 (64 + 96) | 0.01 | 2.53·10⁻³ | 2.35 |

The test (`tests/sliding.rs`) requires orders between 1.9 and 2.8 and errors
below 2·10⁻² and 4·10⁻³ on the 640 and 2 560 cells. The 640-cell case, 500 steps,
runs in 5.6 s, the 2 560-cell case, 1 000 steps, in 45 s (8.4 and 51 s before
PR #97, measured on the same machine; finding 3).

**Equal pitch against conformal stepping** (32 + 32 sectors, 512 cells, the time
step of variant A, 2π/(32 · 10) s): L2 error 1.350·10⁻² against 1.354·10⁻² of
conformal stepping; the test requires agreement within 10 %. **Outer iterations:**
two give 1.35·10⁻² on the 640-cell case, four 1.29·10⁻², eight the same.
**Solid-body rotation** (both cylinders at 10 rad/s, started from the exact
solution; the rotor mesh turns through the fixed stator mesh, the fluid flows
through the stator cells): 1.02·10⁻² at 0.6 s and 0.97·10⁻² at 3 s, without
growth (the test requires below 3·10⁻² and the long run below 1.5 times the
short one). **Passive scalar** (θ = 1 inside, 0 outside): 4.9·10⁻³ against the
exact ln(r₂/r)/ln(r₂/r₁) (the test requires below 1.5·10⁻²). The **mass
balance** closes to 10⁻¹⁸ in every step (the test requires 10⁻¹²). **Outputs:**
a probe at r = 0.06 m in the rotor gives u_θ = 0.3483 m/s against the exact
0.3556 m/s, one at r = 0.09 m in the stator 0.0690 against 0.0704 m/s (the
2 % of the 640-cell mesh); the fields are written with the interface cells as
polyhedra.

Findings:

1. **The imprint closes the cells.** The pieces of a face of an extruded mesh
   add up to its area vector, and the slivers at the top and bottom edge of
   the interface cancel in the sum; the cell volumes and centroids are
   overridden by the rigid-body values of the base mesh, so the geometric
   conservation law is exact (w·S = 0 on a surface of revolution) and the
   closure defect |ΣS|/Σ|S| is round-off (`SLD006` warns above 10⁻⁸).
2. **Fully implicit in the step.** With the extrapolation of the previous
   step an explicit term is made of the unrotated vectors of the previous angle,
   which is stale on the turned mesh and unstable; with `extrapolation = false`
   and two outer iterations the old time levels of the BDF2 scheme, the cell
   content of the moved cells, are all the scheme needs of the step before.
3. **Rebuilding costs (PR #97).** The solver, its matrices and multigrid
   hierarchy are built again in every step because the face list changes. On
   the case file (main before #97, 8.4 s) 48 % of the solver time was the
   pressure multigrid setup and 1.15 s of the run the rebuild of mesh and
   solver around it. Measured and dropped first: **keeping the hierarchy of
   the step before.** The rotor turns by about one cell per step (0.2 rad
   against a sector of 0.196 rad), so its coarse levels no longer fit the
   interface: kept whole, the conjugate gradients needed 784 iterations per
   step instead of 72 (run 14.7 s); with only its finest level refreshed (new
   matrix, diagonal and λ_max, old coarse levels) 353 (8.7 s). A hierarchy
   built for every step is the cheapest. What #97 changed instead, each with
   the same result (histories bit-identical, or to 1e-13 where the coarse
   inverse is rounded differently; L2 errors unchanged): the dense inverse of
   the coarsest level, 80 % of a setup on this mesh, now goes through
   L⁻¹ by rows and A⁻¹ = L⁻ᵀL⁻¹ on the lower triangle (5 ms to 0.9 ms for
   175 unknowns); setups below 20 000 rows run on one helper thread instead
   of handing each of their short parallel loops (20 power iterations per
   level among them) to the pool; the mesh geometry of a step is computed
   once (it was twice, with cell volumes that the rigid-body cells replace);
   the discretisation once (the second one served only the least-squares
   weights of the pressure, now computed alone); the state is moved, not
   copied. Run time on the same machine, before and after:

   | Mesh | Steps | Before | After | Multigrid setup before / after |
   |---|---|---|---|---|
   | 640 cells (case file) | 500 | 8.4 s | 5.6 s | 3.36 / 1.28 s |
   | 2 560 cells | 1 000 | 51.2 s | 44.5 s | 6.17 / 4.92 s |
   | 40 960 cells | 20 | 30.8 s | 28.4–30.8 s | 1.71 / 1.37–1.57 s |

   On 40 960 cells the pressure and momentum solves are 90 % of the run and the
   difference is within the noise of the machine. What is left of the rebuild:
   the new solver (buffers and kernel data, 0.35 ms per step on 640 cells,
   26 ms on 40 960, 2 % of that run) and the discretisation; a rebinding in
   place would save these and is solver-speed work for later.
4. **Found on the way:** the VTKHDF writer classified a cell with 6 quadrilateral
   faces as a hexahedron even when its interface face had corners of its own
   (it panicked); a cell whose faces do not close on the nodes of its shape is
   now a polyhedron. Probes in the rotor sampled the cell that had been at
   their point at the start; they are located again in every step.
5. Not tested yet: a restart of a sliding run (refused, `SLD005`), turbulence
   models (refused, `SLD004`, until PR #100: below), more than one rotor, the
   time-dependent speed of D182, and the combination with the energy equation
   (stage 7). A rotor with blades: below.

## A rotor with blades: the stirred vessel (part 6g, PR #98)

The cases above turn cylinders, whose walls move along themselves. The 2D
baffled stirred vessel of gallery 98 (`examples/studies/mixer/`, Gmsh, 16 326
cells; four blades to r = 0.045 m inside the interface r = 0.06 m, four baffles
on the wall r = 0.10 m, Ω = 10 rad/s) runs with `[[mrf]]` (steady, Re = 100,
972 iterations) and with the sliding mesh by imprinting (64 against 96
interface faces, from rest through three revolutions, Re = 400, 1 508 steps).
The mass imbalance of the frame run is 5·10⁻¹⁸ kg/s.

Findings:

1. **A wall of a turning zone did not push the fluid.** The wall velocity of
   D53 keeps only its tangential part (a velocity normal to a wall that stands
   would make it an inlet, `BC010`), and that rule was applied to the walls of
   a turning zone, too: a blade, which moves normal to itself, slid along the
   fluid. In the frame run the fluid between the blades stood nearly still
   (highest speed 0.18 m/s, at the blade tips). A wall whose cell is in a zone
   of `[[mrf]]` or `[[sliding]]` now takes its normal velocity from the zone,
   n·(Ω × (x_f − x₀)), and its tangential velocity from `velocity`; its
   absolute flux is ρ **S**·**V** and its relative flux zero, as D55 states.
   With the fix the core between the blades turns with them (u_θ = 1.05 Ωr
   near the hub, 0.96 Ωr near the tips) and the highest speed is 0.43 m/s.
   The cylinders of the cases above are not affected (their rotation is
   tangential to them).
2. A square box turning with its frame, all four walls moving normal to
   themselves (`a_turning_box_carries_its_fluid`, 16 × 16 cells): L2 error of
   the rigid rotation 1.8·10⁻² before the fix, 1.3·10⁻³ after. The rest is the
   wall pressure, extrapolated with zero normal gradient where the true one is
   ρΩ²r·n; the test bound is 3·10⁻³.
3. Seen in the fields, not measured against data: in the frame of the blades
   the fluid of a blade passage turns in a closed eddy against the rotation
   (the relative eddy of impeller theory), and the blades shed starting
   vortices that wrap around the rotor and merge into a shear layer at the
   interface radius. No reference solution exists for this geometry; the case
   is a demonstration and a regression of the wall fix, not a validation.

## Fields of rotating machinery (D183, PR #99)

Five derived fields for `[[mrf]]` and `[[sliding]]`: the relative velocity
**w** = **c** − Ω × (x − x₀), the frame velocity, the cylindrical components
with the flow angles, the relative total pressure and the rothalpy.

1. Taylor–Couette flow in the frame (1 024 cells,
   `relative_and_frame_velocities_of_taylor_couette`): the frame velocity is
   Ω × x in the rotor to 10⁻¹² and zero in the stator, Ω × x everywhere with
   `output.frame = "rotor"`; c_u − w_u = Ωr and ρI = p + ½ρ|w|² − ½ρ|u|² to
   round-off; the largest deviation of c_u from the exact u_θ is 0.7 % of the
   wall speed, of c_r 2.9 % (the L∞ error of the velocity is 2.95 %), α below
   1.8°. The fields come from the restart file (`foehn post derive`), so the
   path after the run is checked as well.
2. The identities on 200 random velocities, points, axes and both senses of
   rotation (`turbomachinery_identities`): c = w + u, c_r² + c_u² + c_a² =
   |c|², c_u − w_u = |u| = |Ω| r, the angles from their definitions in [0°,
   180°], the fixed frame, a fluid turning with its frame (w = 0, α = 0), the
   axis itself (no radial or circumferential direction).
3. **Found on the way:** the sample lines and the reference comparison of the
   sliding mesh by imprinting located their points with the locator of the
   start of the run; the rotor cells had moved on, so a line in the rotor
   read the cell that had stood there at t = 0 and extrapolated its value
   with the gradient (0.68 m/s instead of 0.47 m/s at r = 52 mm). They use the
   mesh of the end of the run now; `probes_and_fields_follow_the_rotor` checks
   a radial line against the exact solution and c_u − w_u = Ωr in the rotor
   and 0 in the stator. The log columns `probe.<name>.<q>` are still located
   once, at the start (open).
4. On a sliding run of the stirred vessel (100 steps) the frame velocity in
   the written fields matches Ω r at the turned positions of the rotor cells.

## Turbulence on the frame and the sliding meshes (part 6g, PR #100)

Turbulent Taylor–Couette flow (k-ω SST, automatic wall treatment, ρ = 1000
kg/m³, ν = 2·10⁻⁶ m²/s, Re = Ω r₁ (r₂ − r₁)/ν = 12 500) and a turbulence patch
carried by the rotor; the sample line across the gap compares u_θ, k and ν_t.
The tests are in `tests/mrf.rs` and `tests/sliding.rs`, the wall distance in
`crates/foehncraft-mesh/src/query/turning.rs`.

1. **Wall distance of a turning rotor**
   (`turned_parts_give_the_brute_force_distance`,
   `the_imprinted_ring_agrees_with_the_wall_distance_of_its_mesh`): one
   bounding volume hierarchy per rigid part, the rotor walls at angle 0 and
   the stator walls; a cell searches the other part only when it is nearer to
   the edge of its part (the interface and the other open faces, without the
   planes that bound the whole domain at every angle) than to the walls of
   its part, with the query capped at that distance. On random points between
   turning blades and baffles it equals the brute-force distance at every
   angle tried, and on the split annulus (12 rotor against 18 stator sectors)
   the full distance of the imprinted mesh to 10⁻¹²; at most half the cells
   search the other part.
2. **The frame** (`turbulent_taylor_couette_does_not_depend_on_the_frame`,
   steady, 2 048 cells, first cells at y⁺ 0.8): the inner half of the gap in
   the frame of the cylinder against the inertial run, largest deviation on
   the line 0.24 % for u_θ, 0.43 % for k, 0.34 % for ν_t (of the largest
   value), less with twice the sectors. **Found:** with the first cell in the
   buffer layer (y⁺ 5.6, `TURB001`) its k differed between the frames by
   13 %, 44 % and 67 % with 64, 128 and 256 sectors; but in the inertial run
   alone this k fell by a factor of 2.8 from 32 to 256 sectors, the interior agreed to 0.4 %, and a
   restart of either run from the other's solution went back to its own. It is
   the first cell that the automatic treatment does not resolve, not the
   frame; the test uses the wall-resolved mesh.
3. **Two independent frames** (`two_independent_frames`, laminar): the inner
   half of the gap turning at 10 rad/s with the inner cylinder, the outer
   half at −4 rad/s with the outer one; L2 error against A r + B/r 0.78 %,
   the inertial run 0.72 %.
4. **Spin-up on both sliding meshes** (`turbulent_spin_up_agrees_with_a_fixed_mesh`,
   2 s from rest, 102 steps, 768 cells, the imprinted mesh with 32 rotor
   against 48 stator sectors): the profiles against the fixed mesh with the
   turning wall, largest deviation 1.6 % for u_θ, 1.5 % for k and 4.0 % for
   ν_t, the same for both variants; with half the width and half the step
   about 0.5 %, 1 % and 1.5 %.
5. **A patch carried by the rotor** (`turbulence_is_carried_by_the_rotor`,
   k-ε-ζ-f, fluid in solid-body rotation with the walls, ν = 5·10⁻⁴ m²/s, half
   a turn in 16 steps): the patch arrives at the opposite side with the peak
   k of the same patch decayed at rest within 12 % (conformal) and 7 %
   (imprinted), 7 mm short of the place of the resting peak (the velocity
   error of the coarse mesh), and leaves 10 % and 7 % of its peak behind.
   Since PR #141 (the wall sink of k) the imprinted run takes 32 steps:
   with 16 (Courant 1.7) the rotor drives ζ to its cap and the run
   diverges in step 8, also with the papers' −ε; the sink k/T_K of 6d2
   held there (finding). The peaks now agree within 4 % (conformal) and
   4 % (imprinted), at the place of the resting peak.
   **Found:** with conformal stepping ζ-f blew up within ten steps (`GEN001`,
   Courant 598): the time terms used the moved old levels, but each step
   started from the fields at the places of the fixed lattice, so the first
   iteration's coefficients belonged to turbulence that had moved on. The
   start of a step now comes from the moved content (extrapolated in time
   like the old levels), f included.
6. **The stirred vessel, turbulent** (gallery 100): water, Re = Ω R_tip²/ν =
   20 000, k-ω SST, 16 326 cells. **Found:** the frame run failed in the first
   iteration (`GEN001`, also on main): on the second level of the pressure
   multigrid seven rows had a lumped diagonal of the filtered matrix at or
   below zero, the power iteration returned a non-positive λ, and the
   smoothed prolongator and the third level became NaN. Such a row keeps its
   own diagonal now and λ falls back to the Gershgorin bound
   (`filtered_diagonal_stays_positive`). SIMPLEC then diverges from rest at
   this Reynolds number (iteration 6, also without the frame); the coupled
   solver converges (337 iterations to 10⁻⁵; the ω residual stays at
   2·10⁻⁶ while the others fall to 10⁻¹⁰). The sliding mesh by imprinting
   runs three revolutions from rest (1 508 steps, 10 min on four cores
   against 7 min laminar); 2 118 of 16 326 cells search the other part for
   their wall distance in each step. In the fields k and ν_t step at the
   interface, where the rotor cells are 1.5 times the size of the stator
   cells (64 against 96 intervals); a passive scalar diffusing from the
   blades to the vessel wall (no convection, D = 10⁻³ m²/s, 8 s) on the
   imprinted mesh agrees with the steady solution on the conformal mesh to 1
   to 2 % in every 2 mm ring from r = 50 to 72 mm, with the same steps in the
   ring averages on both meshes, so the transport across the interface is
   not at fault.
7. Not tested yet: turbulence with a sliding restart (refused, `SLD005`), the
   EB k-ε model on the sliding mesh, and a comparison with measured data of a
   stirred vessel (3D).


## Torque, shaft power and machine values (D184, PR #101)

The tests are in `tests/mrf.rs`, the code in `crates/foehncraft-flow/src/rotor.rs`.

1. **Torque of laminar Taylor–Couette flow** (`torque_and_shaft_power_of_taylor_couette_flow`;
   r₁ = 0.05 m, r₂ = 0.10 m, μ = 5·10⁻⁴ Pa s, ω = 10 rad/s, 1 024 cells):
   the torque of the inner cylinder −4π μ ω r₁² r₂² L/(r₂² − r₁²) = −2.094·10⁻⁶ N m
   to 0.47 % with the frame at the cylinder speed, 0.94 % in the inertial
   frame and 1.26 % in a frame at 25 rad/s; the outer cylinder carries the
   opposite torque (angular momentum imbalance 0.2 %, 0.7 %, 1.0 %). Without
   the transposed viscous stress the torque was 1.37·10⁻⁶ N m (35 % low): the
   Laplacian form ∇·(μ∇u) of the momentum rows has no ∇uᵀ term, which gives the
   shear stress μ(du_θ/dr − u_θ/r) its second part. A wall at rest is not
   affected (the term vanishes at a no-slip wall).
2. **Shaft power** in the same test: −Ω·T equals the work of the wall on the
   fluid to 10⁻¹⁰ (they are the same sum by construction for a wall that
   turns with the zone), the exact dissipation Ω |T| to 3 %, and the ledger
   books it: shaft work = wall power of the Laplacian form + work of the
   transposed stress (5·10⁻³).
3. **Machine values of plane Poiseuille flow** (`machine_values_of_poiseuille_flow`):
   through a zone at Ω = 0, Δp_t = −12 μ U L/h² = −2.4 Pa to 0.4 %,
   Q = U h d to 0.2 %, hydraulic power Q Δp_t, no Euler work, an undefined
   efficiency at zero shaft power.
4. **Sliding mesh:** a run of 40 steps of the stirred vessel shows the angle
   column (|Ω| t in degrees: 21.49° at t = 0.0375 s), the shaft power in the
   compact balance line and the torque columns of `forces.csv`.

## Several rotors and the transient frozen rotor (D185, PR #102)

The tests are in `tests/sliding.rs`, `tests/mrf.rs`,
`crates/foehncraft-mesh/src/imprint.rs` and `query/turning.rs`.

1. **Two rotors in one mesh** (`two_rotors_turn_at_their_own_angles`): an
   annulus of three rings (rotor of 12, stator of 18, outer rotor of 14
   sectors, two layers in z) assembled at four pairs of angles: every cell
   closes to 10⁻¹², the rotor cells are turned by the angle of their own
   rotor and keep their volumes, the pieces are those of the two interfaces
   alone, and the same rotor given twice is refused.
2. **Wall distance of two turning rotors**
   (`two_rotors_agree_with_the_wall_distance_of_their_mesh`): one hierarchy
   per part, the point turned into the frame of each other part, equals the
   wall distance of the assembled mesh at four pairs of angles to 10⁻¹².
3. **Couette flow of two counter-rotating cylinders**
   (`two_rotors_end_in_the_couette_flow_of_their_speeds`,
   `couette_two_rotors.toml`; r₁ = 0.05 m at ω₁ = +10 rad/s, r₂ = 0.10 m at
   ω₂ = −5 rad/s, each cylinder a sliding rotor by imprinting with its own
   pitch, 960 cells, 500 steps of 0.04 s from rest): L2 error 1.76 % against
   u_θ = A r + B/r; the torques −3.167·10⁻⁶ and −3.163·10⁻⁶ N m (each about
   the axis along its own Ω) against 4π μ L r₁² r₂² (ω₂ − ω₁)/(r₂² − r₁²) =
   −3.142·10⁻⁶ N m (0.8 % and 0.7 %); the shaft powers 3.17·10⁻⁵ and
   1.58·10⁻⁵ W add up to the exact Ω-weighted work 4.71·10⁻⁵ W (0.8 %) and to
   the dissipation of the ledger.
4. **Transient frozen rotor** (`the_spin_up_does_not_depend_on_the_frame`,
   `taylor_couette_mrf.toml` with `[time]`, 25 steps of 0.02 s): at
   t = 0.5 s the L2 distance from the steady state is 0.1299 in the inertial
   frame, 0.1296 in the frame of the cylinder and 0.1308 in the opposite
   frame (0.7 %); after 1 000 steps the transient frame ends at the L2 error
   of the steady frame, 6.93966·10⁻³ against 6.93966·10⁻³.

## Interfaces of any rotationally symmetric shape (D186, PR #103)

The tests are in `crates/foehncraft-mesh/src/imprint.rs`,
`crates/foehncraft-mesh/src/revolved.rs` and `tests/sliding.rs`.

1. **The meshes of revolution** (`revolved.rs`): a stage on the hub, a pot
   with wedge cells at the axis, a disc impeller as a solid block and the
   cones of `cone = 20`: every cell closes to 10⁻¹³, the cell count is the
   estimate, the boundaries are named as documented, crossing breakpoints
   and a block both rotor and solid are refused.
2. **Closed cells at every angle**
   (`surfaces_of_revolution_close_at_every_angle`): a stage (disc, cylinder,
   disc on the hub), its cones (25°), a pot whose two discs pass through the
   axis and a rotor row between two separate discs (−15° cones, two parts of
   the meridian), 12 rotor and 17 stator intervals, assembled at seven
   angles: every cell closes below 10⁻¹², the shape is recognised ("3
   straight pieces in 1 part, 2 through the axis", "2 straight pieces in 2
   parts"), and the pieces of every rotor face add up to its area within
   the difference between the chords and the surface. The cylinder tests of
   PR #96 and #102 pass unchanged.
3. **Uniform axial flow** (`uniform_flow_crosses_discs_cones_and_the_axis`,
   `axial_stage_uniform.toml`, slip walls, rotor at 40 rad/s, 20 steps):
   through the pot, away from all walls, the uniform flow is kept to
   1.5·10⁻⁶ (relative L2; the same mesh with aligned pitches and the rotor
   at rest: 1.8·10⁻¹²). Through the stage on the hub: 2.30·10⁻³, 7.57·10⁻⁴
   and 2.22·10⁻⁴ at 24/34, 48/68 and 96/136 intervals (ratios 3.0 and 3.4),
   with cones 2.36·10⁻³ and 7.87·10⁻⁴: the step between the hub polygons of
   rotor and stator where the interface meets the hub, second order in the
   pitch; with aligned pitches and the rotor at rest 4.3·10⁻¹².
4. **Solid-body rotation in the pot**
   (`solid_body_rotation_is_kept_by_a_pot_through_the_axis`,
   `pot_solid_body.toml`, every wall at 10 rad/s, 30 steps, half a
   revolution): L2 error 0.82 % against u = Ω × r (the cylinder imprint of
   PR #96 on its ring: below 3 %).

The disc impeller of gallery 103 (`disc_pot.toml`) at first ran only at
ν = 10⁻⁵ m²/s and diverged in step 17 at ν = 10⁻⁴ m²/s; the cause was
the PISO corrector, fixed in PR #106 (next section).

## The PISO corrector on thin wedge cells at an axis (PR #106)

**Finding.** `disc_pot.toml` at ν = 10⁻⁴ m²/s diverged in step 17, also
with half the time step and with the cell-quality remediation. Narrowed
down on that case with 40 steps:

| Variant | Result |
|---|---|
| as set (PISO, 2 correctors, 3 outer iterations, BDF2) | Courant number ×30 per step from step 4, diverged in step 17 |
| without the non-orthogonal and skewness correction | diverged in step 18 |
| Stokes flow (no convection) | ×45 per step, unbounded |
| rotor at rest (Ω = 10⁻⁶ rad/s), Stokes | ×45 per step |
| every axis cell in the stator (rotor block away from the axis), Stokes | ×40 per step |
| without the axis (`r[0]` = 1.5 mm, symmetry), as set | grows later and slower (Courant number 126 in step 30) |
| cell-quality remediation off / BDF1, Stokes | ×40 / ×110 per step |
| 50 instead of 1 000 momentum iterations | the same to four digits |
| 8 outer iterations, Stokes | ×5·10⁴ per step: ×3.9 per outer iteration |
| 4 correctors, Stokes | stable |

So it is neither the sliding mesh nor the interface nor the explicit
parts of the viscous term: it is an instability of the outer iteration of
PISO on the axis cells, and more correctors suppress it. The growing mode
is axisymmetric: a radial velocity towards the axis in the wedge cells,
outwards in the ring next to them, and a pressure that alternates along
the axis. Printed per stage, the pressure grows by 2.4 in every outer
iteration: the first corrector moves it by ten times its value and the
second takes most of that back.

**Cause.** The wedge cells at the axis are thin around the axis (72 to a
turn in the stator), so at high viscosity their a_P is almost all coupling
to their two neighbours around the axis, about 50 times the time term
ρV/Δt. Issa's corrector takes D = V/a_P, which neglects the neighbours: it
underestimates the velocity response to the pressure there by
a_P/(a_P − Σ|a_N|). The correctors then leave a pressure that the momentum
predictor of the next outer iteration, which sees the full coupling,
answers too strongly. (On a uniform cavity at the same diffusion number,
Σ|a_N|/a_P = 0.97, Issa's form is stable; the wedges add the rotation of
the Cartesian components between neighbours and the anisotropy.)

**Fix.** With outer iterations the correctors now use the SIMPLEC
coefficient D^C = V/(a_P/α − Σ|a_N|) (Van Doormaal & Raithby 1984), and
between two correctors the consistent neighbour correction
δ = Σ_N |a_N| (u'_N − u'_P)/a_P (D^C already counts Σ|a_N| u'_P;
`foehncraft_kernels::flow::neighbour_difference`). The Rhie–Chow flux
keeps D, so a converged step is the same as before. Changing only δ while
keeping D made it worse (×1 000 per step): the two belong together.
Plain PISO (one outer iteration) keeps Issa's form: there it damps more.
The lid-driven cavity stepped to its steady state with Crank–Nicolson at
Courant number 4 (`steady_state_of_the_time_stepping_is_the_simplec_solution`)
settles to 10⁻¹² with Issa's corrector, while with the SIMPLEC form its
residual grew to 0.17 in 1 200 steps; with three outer iterations the
SIMPLEC form settles too (10⁻¹⁶ in 300 steps). Hence `solver.corrector =
"auto"` (default): `simplec` with `outer_iterations` > 1, `issa` for plain
PISO; `"issa"` and `"simplec"` force one. Every plain-PISO run is bitwise
unchanged.

**Checks.**

1. **Solid-body rotation at high viscosity**
   (`solid_body_rotation_is_kept_at_high_viscosity_on_the_axis_wedges`,
   `pot_solid_body_viscous.toml`: `pot_solid_body.toml` at ν = 10⁻²
   m²/s, the wedges' coupling about 90 times their time term): L2 error
   0.11 % after 30 steps; with Issa's corrector 2.0 (the mode grows from
   round-off). `pot_solid_body.toml` at ν = 10⁻³: L2 0.82 % with either corrector
   (the error there is spatial).
2. **The disc impeller at ν = 10⁻⁴ m²/s** (`disc_pot.toml`, 400 steps,
   eight revolutions, 8.5 minutes on four cores): runs through, Courant
   number about 1.04 throughout; torque on the disc −2.28·10⁻³ N m (the
   two faces of a free disc: 1.40·10⁻³ N m, C_M = 3.87/√Re at Re = 180;
   the rim and the vessel add to it), shaft power 46 mW, the energy ledger
   closes with 6.8 % numerical dissipation, no step with negative
   numerical dissipation. The vessel is still spinning up after 2.5 s (the
   vessel walls take 1.29·10⁻³ N m of the 2.28).
3. **The other transient checks** (`tests/transient.rs`, `sliding.rs`,
   `mrf.rs`, `stage5.rs`): pass. Plain PISO is bitwise unchanged; the
   cases with outer iterations (sliding meshes, the transient frozen rotor)
   change within their tolerances.

## CPU–GPU parity of the turbulence models (6.6, PR #109)

`tools/bench/stage6_cloud.py` (Colab: `tools/notebooks/stage6_colab.ipynb`)
runs 16 cases, 300 outer iterations each, on the CPU and on the GPU and
compares every pair bitwise with `tools/check/compare_runs.py`: the five
channels (SST, k-ε, ζ-f, ζ-f₀, EB), the TMR plate with all six models,
the heated channel and pipe (passive scalar with the thermal wall function),
the coupled algorithm on the backward-facing step (8 200 cells), and since
PR #118 the two SST options of D170 (the rotating channel with the
curvature correction, the square duct with QCR2000).

| Check | Result |
|---|---|
| `--cpu-only --iterations 60` in the cloud (one thread against four; these meshes are below the 20 000 cells under which one thread is used, so this tests the script and the comparison) | all 14 pairs bitwise equal (9 to 11 files each) |
| Tesla T4 (sm_75), CPU (2 threads) against GPU, the first 14 cases, run by the maintainer on main 2026-11-20 | all 14 pairs bitwise equal (9 to 11 files each: fields, residual histories, wall values, restart) |
| The same with the two SST options of D170 (added after that run) | pending: run by the maintainer |

Run times of the T4 run (300 outer iterations, wall time with mesh
generation and output; ns per cell and iteration):

| Case | Cells | CPU s | GPU s | CPU ns | GPU ns |
|---|---:|---:|---:|---:|---:|
| channel_sst | 192 | 0.5 | 6.0 | 6 707 | 26 709 |
| channel_keps | 24 | 0.1 | 1.1 | 3 518 | 105 122 |
| channel_zeta_f | 192 | 0.7 | 2.3 | 10 639 | 33 212 |
| channel_zeta_f0 | 96 | 0.5 | 2.1 | 14 607 | 59 495 |
| channel_eb | 192 | 0.8 | 2.7 | 12 756 | 40 564 |
| tmr_plate_sst | 3 072 | 4.4 | 3.6 | 4 587 | 3 206 |
| tmr_plate_keps | 3 072 | 2.9 | 2.7 | 3 690 | 2 916 |
| tmr_plate_keps (realizable) | 3 072 | 5.8 | 3.5 | 6 041 | 3 077 |
| tmr_plate_zeta_f | 3 072 | 9.7 | 4.9 | 10 348 | 4 675 |
| tmr_plate_zeta_f0 | 3 072 | 12.6 | 4.6 | 13 369 | 4 339 |
| tmr_plate_eb | 3 072 | 11.8 | 4.8 | 12 576 | 4 580 |
| channel_heated_keps | 24 | 0.1 | 1.5 | 5 707 | 147 207 |
| pipe_heated_keps | 72 | 0.1 | 1.7 | 3 943 | 58 619 |
| backstep_zeta_f0_coupled | 8 200 | 27.4 | 10.4 | 27 136 | 8 678 |

On these small meshes the times measure the launch cost, not the
throughput. Below about 200 cells the GPU is 3 to 17 times slower. From
3 072 cells on it is 1.1 to 2.7 times faster than two CPU threads, and
2.6 times on the step.

`compare_runs.py` now leaves out the hash and size of `case.resolved.toml`
in `result.json`: the file's lines are compared without the device keys of
`[run]`, and `--threads` or `--backend` change its hash by design.

## Derived fields of the turbulence models (D72, PR #108)

`[output] derived` and `foehn post derive` gain `turbulence-scales`,
`reynolds-stress`, `wall-normal` (EB) and `turbulence-wall` (spec 2.3),
computed from the turbulence state of the solution, so the run and the
restart file give the same numbers.

| Check | Quantity | Band | Result |
|---|---|---|---|
| Channels of SST, k-ε and EB, 24 cells, 1 500 iterations (`turbulence_fields_of_every_model`) | ε = C_μ k ω, L = k^{3/2}/ε, ν_t/ν; T = k/ε (SST, k-ε) or the T of the model (EB) | 1e-10 relative; T bitwise | asserted |
| | trace of the Reynolds stress against 2k; sign of u'v' away from the centre line; u'v' with z zero | 1e-6 k | asserted |
| | u_τ of `turbulence-wall` against the wall treatment; y* = d C_μ^¼ k_P^½/ν; zero off walls | bitwise; 1e-10 | asserted |
| | EB: \|n\| = 1, n_y = ±1 pointing away from the nearer wall | 1e-12; 0.999 | asserted |
| Laminar case, SST with `wall-normal` (`turbulence_fields_are_refused_without_their_model`) | refused | CASE005 | asserted |
| The step, k-ε-ζ-f₀, 8 200 cells (gallery 108) | largest −u'v'/U_ref² in the shear layer | reported | 0.011 at x/H 4 to 7 |

Not in this part: the blending weight of the compound and automatic
treatments as a wall field (added in PR #130, below); the
heat-transfer coefficient and T⁺ of the scalars stay in the wall CSV.

### The wall blending weight (spec 2.3, PR #130)

The wall kernel writes the viscous weight of each wall treatment to an
output buffer (not read by the solver), the restart file stores it as
`turb_wall_blending`, and `turbulence-wall` gains the boundary dataset
`blending`. A restart file written before has no such dataset; `foehn post
derive` then refuses `turbulence-wall` with CASE005 and asks for one more
iteration.

| Check | Quantity | Band | Result |
|---|---|---|---|
| Channels of `turbulence_fields_of_every_model`, plus SST with `automatic` | `blending` in [0, 1], zero off walls | exact | asserted |
| | SST `low-re` 1, k-ε `scalable` 0 | bitwise | asserted |
| | SST `automatic` against 1 − tanh((y⁺/5)²) with the y⁺ of the same call | 1e-12 | asserted |
| | EB `compound` against e^{−Γ(y*)}, Γ = 0.01 y*⁴/(1 + 5 y*), with the y* of the final k | 1e-3 relative (k of the last wall call is one iteration older) | asserted |

## The SST automatic treatment in the buffer layer (6g, PR #112)

The first version of `automatic` (6b) overpredicted C_f of the plate
by up to 17 % with the first cell at y⁺ 3 to 20: the first cell took ω_P
with the weight 1 − u_τ,vis⁴/u_τ⁴, which is still small there, and the
log-layer production u_τ³/(κ y_P) with the same weight, so the cell's own
production from U_t/y_P, far too large above the sublayer, set k. The
sweep `examples/studies/sst_wall_sweep.py` (the plate on 64 × 48 cells,
the wall grading set for each first-cell y⁺ at x = 0.97) measured each
change in turn:

| Variant | y⁺ 1.4 | 4 | 8 | 24 | 85 |
|---|---|---|---|---|---|
| 6b: ω_P blend n = 2, weight 1 − b, log-layer production | +0.6 % | +13.4 % | +16.6 % | −0.2 % | −1.4 % |
| equilibrium production (u_τ⁴/ν)(1 − g) g, g = du⁺/dy⁺ of Reichardt | +0.6 % | +10.8 % | +4.2 % | −0.9 % | −2.9 % |
| `low-re` on the same meshes, for comparison | +0.6 % | +11.4 % | +46.3 % (not converged) | | |
| ω_P imposed fully, n = 1.2, equilibrium production | −7.5 % | −7.4 % | −1.4 % | +0.9 % | −2.1 % |
| **6g: as above with the weight tanh((y⁺/5)²) for ω_P and the production** | −1.0 % | −2.4 % | −1.4 % | +0.9 % | −2.1 % |

| Test | Check | Band | Result |
|---|---|---|---|
| The plate on 64 × 48 cells, first-cell y⁺ 0.19, 0.67, 1.2, 2.1, 3.3, 4.8, 7.0, 11.8, 21.6, 40.6, 76.6 reached (targets 0.24 to 84) | C_f(0.97) against our own `low-re` result on the finest mesh (2.7005e-3, y⁺ 0.09; the constants were first fitted to the value of another solver and were rechecked against this one on 2026-11-22) | ±2 % (9.3) | −0.9, −1.3, −1.0, −1.4, **−2.5**, **−2.2**, −1.9, +1.0, +1.2, −0.9, **−2.2** %: within 2.5 % |
| The plate on 128 × 96 cells, y⁺ 0.58, 2.9, 20.8 reached | the same | ±2 % | −1.0, **−2.6**, +0.3 % |
| Gallery 68 sweep (y⁺ 1.7 to 92, the chart) | the same | ±2 % | −0.5, **−2.3**, −0.8, +0.7, −1.8 % |
| Step, SST, 11 040 / 44 160 cells (`backstep_validation.py --models sst`) | x_r/H; C_p RMS | ±10 %; ≤ 0.02 | 6.761 / 6.743 (+7.7 %; before 6.793, 6.747); 0.0191 (before 0.023: now inside the band) |
| Jet, SST `kato-launder`, `automatic`, 2 600 cells (`impinging_jet.py --models sst --levels 0`) | Nu₀ | unchanged | 151.6 (before 151.9: −0.2 %; the first cells lie in the sublayer) |
| Kernel tests (`automatic_treatment_tends_to_its_limits`, `automatic_buffer_layer_weights`) | the limits of ω_P, the weight and the production; the production maximum u_τ⁴/(4ν) near y⁺ 11 | 10⁻³ to 2 % | asserted |

Four points of the sweep are 0.2 to 0.6 % outside the ±2 % band; TURB001
now warns only above y⁺ 300, where the log layer ends. The choices n =
1.2 and y₀ = 5 are fitted to this plate (n = 1.1 gives −2.8 to +1.2 %,
against the TMR value as first measured);
the near-wall cell of a separating or impinging flow still takes the
equilibrium production of a constant-stress layer above y⁺ ≈ 5, as every
wall function does.

**The wall value of ω in `low-re`.** The sweep of the factor F of
ω_w = F · 6ν/(β₁ d₁²) (F = 10 by default, Menter 1994) on the plate
(`tmr_plate_sst.toml`, `low-re`, levels 0, 1, 2, y⁺ ≈ 0.2) and the channel
(`channel_sst.toml`, levels 0, 1, 2; one build per factor):

| F | plate C_f(0.97) against our own finest value 2.7005e-3 (the first table was against the TMR value; the same runs, re-expressed) | channel U_b⁺ (DNS 17.55) |
|---|---|---|
| 1 | +13.6, +4.2, +1.8 | 16.940, 17.149, 17.232 |
| 3 | +5.5, +1.1, +0.5 | 17.216, 17.279, 17.296 |
| 10 | −0.3, −1.2, −0.5 | 17.427, 17.381, 17.346 |
| 30 | −3.9, −2.8, −1.3 | 17.575, 17.452, 17.380 |
| 100 | −7.1, −4.2, −1.9 | 17.707, 17.515, 17.411 |

Every factor tends to the same values as the mesh is refined (plate about
−0.1 to −1.5 %, channel 17.23 to 17.41 on the finest mesh). F = 10 has the
smallest mesh dependence on the plate, within 0.8 % already on the
coarsest mesh, and the middle of the spread in the channel. **Finding:**
the default stays F = 10 (the channel runs of this sweep stopped at the
case tolerance 10⁻⁶; on these three meshes that is within 0.1 % of the
converged values, PR #128). The −1.2 % of U_b⁺ in the channel does not come
from the wall value of ω.


## GCI of the turbulence models on the plate and the channel (9.3, 9.4, PR #111)

`examples/studies/plate_channel_gci.py` runs the gallery case of every
model on three meshes (the plate with `verification.levels` 0, 1, 2, the
built-in mesher halving every interval in both directions per level; the
channel, whose solution does not depend on x, refined in y only to the
same y nodes, since PR #128) with one thread each and evaluates the observed order p, the extrapolated value and
GCI_fine after Celik et al. (2008): C_f(0.97) and C_D of the TMR plate,
U_b⁺ of the channel at Re_τ = 395 (from the driving force at the DNS mass
flow). Iteration limits 6 000, 12 000 and 24 000 (the channel, converged
to 10⁻⁹: 3 000 with the coupled algorithm since PR #141, before 20 000
to 160 000 with SIMPLEC; the plate of ζ-f and ζ-f₀: 3 000, coupled); a
run that stops there is marked. The finest k-ε plate (49 152 cells) diverged in outer iteration 11
with the default relaxations (RUN014) and ran with 0.5 for momentum and
turbulence (the study sets this for that level).

## TMR flat plate, Re_L = 5·10⁶

References: C_f(0.97) White 2.87455e-03; Coles–Fernholz at the Re_θ of the finest run of the model (the SST row of the second table).

| model | cells | y⁺ (finest) | C_f ·10³ | p | GCI_fine | C_f extrap. ·10³ | finest vs White | C_D ·10³ | p | GCI_fine | iterations, converged |
|---|---|---:|---|---:|---:|---:|---:|---|---:|---:|---|
| k-ω SST, low-re | 768 / 3072 / 12288 | 0.19 | 2.6937 / 2.6682 / 2.6868 | – (oscillatory) | 0.47 % | – | −6.53 % | 2.8549 / 2.8259 / 2.8477 | – | 0.51 % | 226, 635, 2292 |
| k-ε, scalable | 3072 / 12288 / 49152 | 7.31 | 2.8289 / 2.8270 / 2.8156 | 2.61 (divergent) | 0.10 % | 2.8134 | −2.05 % | 2.7193 / 2.6334 / 2.6768 | – | 1.60 % | 236, 724, 5733 |
| realizable k-ε, scalable | 3072 / 12288 / 49152 | 7.00 | 2.6781 / 2.6128 / 2.5805 | 1.01 (monotonic) | 1.53 % | 2.5488 | −10.23 % | 2.5769 / 2.4216 / 2.4042 | 3.15 | 0.12 % | 6000*, 12000*, 24000* |
| k-ε-ζ-f, compound | 3072 / 12288 / 49152 | 0.23 | 2.9314 / 2.8756 / 2.8634 | 2.20 (monotonic) | 0.15 % | 2.8600 | −0.39 % | 3.0509 / 2.9953 / 2.9832 | 2.20 | 0.14 % | 205, 216, 496 |
| k-ε-ζ-f₀, compound-2023 | 3072 / 12288 / 49152 | 0.23 | 2.9257 / 2.8810 / 2.8689 | 1.88 (monotonic) | 0.20 % | 2.8644 | −0.20 % | 3.0595 / 3.0100 / 2.9963 | 1.85 | 0.22 % | 204, 279, 178 |

The ζ-f and ζ-f₀ rows are those of PR #141 (the wall sink of k,
specification 3.3.2; the coupled algorithm). Before, with −ε (ζ-f) and
k/T_K (ζ-f₀): ζ-f 3.0231 / 3.5516 / 3.3467 (oscillatory, GCI 7.9 %, all
three at the iteration limit), ζ-f₀ 2.8021 / 2.7923 / 2.7400 (divergent,
GCI 0.55 %, the finest at the limit).
| elliptic-blending k-ε, compound | 3072 / 12288 / 49152 | 0.23 | 2.8893 / 2.8408 / 2.8516 | – (oscillatory) | 0.85 % | – | −0.80 % | 2.9012 / 2.6880 / 2.6007 | 1.29 | 2.91 % | 6000*, 2605*, 4726* |

**SST on four meshes, against the fit to measurements (2026-11-22, D201
rework).** The plate of SST on levels 0 to 3 (coupled algorithm, the three
finest meshes for the GCI), with the Coles–Fernholz fit (Nagib, Chauhan &
Monkewitz 2007) at the Re_θ of the finest run (θ and the edge velocity
u_e, 0.26 % above the inlet velocity, from the profile at x = 0.97;
C_f on u_e²):

| model | cells | y⁺ (finest) | C_f ·10³ | p | GCI_fine | C_f extrap. ·10³ | finest vs White | Re_θ (finest) | finest vs Coles–Fernholz | C_D ·10³ | p | GCI_fine | iterations, converged |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|
| k-ω SST, low-re | 768 / 3072 / 12288 / 49152 | 0.09 | 2.6937 / 2.6682 / 2.6868 / 2.7005 | 0.44 (monotonic) | 1.76 % | 2.7385 | −6.06 % | 7848 | +1.45 % | 2.8549 / 2.8259 / 2.8476 / 2.8605 | 0.75 | 0.83 % | 226, 635, 622, 411 |

**The wall constants against our own result (2026-11-22).** The constants
n = 1.2 and y₀ = 5 of `automatic` were first fitted to the value of
another solver. Rechecked against our own `low-re` result on the finest
plate mesh (2.7005e-3, first-cell y⁺ 0.09), with the sweep
`examples/studies/sst_wall_sweep.py` on 64 × 48 cells and first-cell y⁺
(reached) 0.67, 2.1, 3.3, 7.0, 11.8, 21.6, 40.6, 76.6, the deviation of
C_f(0.97) in %:

| n, y₀ | 0.67 | 2.1 | 3.3 | 7.0 | 11.8 | 21.6 | 40.6 | 76.6 | largest |
|---|---|---|---|---|---|---|---|---|---|
| 1.2, 5 (the default) | −1.3 | −1.4 | −2.5 | −2.0 | +1.4 | +1.2 | −0.9 | −2.2 | 2.5 |
| 1.0, 5 | −1.3 | −1.8 | −3.9 | −4.2 | −1.4 | −1.1 | −2.8 | −3.4 | 4.2 |
| 1.1, 5 | −1.3 | −1.5 | −3.1 | −3.0 | +0.1 | +0.2 | −1.8 | −2.7 | 3.1 |
| 1.3, 5 | −1.3 | −1.2 | −1.9 | −1.1 | +2.5 | +2.1 | −0.2 | −1.8 | 2.5 |
| 1.5, 5 | −1.3 | −1.0 | −1.1 | +0.4 | +4.3 | +3.5 | +0.7 | −1.3 | 4.3 |
| 2.0, 5 | −1.3 | −0.8 | 0.0 | +2.8 | +7.2 | +5.7 | +2.0 | −0.7 | 7.2 |
| 1.2, 3 | −1.6 | −5.1 | −7.8 | −2.2 | +1.4 | +1.2 | −0.9 | −2.2 | 7.8 |
| 1.2, 4 | −1.4 | −3.0 | −5.3 | −2.1 | +1.4 | +1.2 | −0.9 | −2.2 | 5.3 |
| 1.2, 6 | −1.2 | −0.2 | +0.1 | −1.6 | +1.4 | +1.2 | −0.9 | −2.2 | 2.2 |
| 1.2, 8 | −1.2 | +1.2 | +4.0 | −0.2 | +1.4 | +1.2 | −0.9 | −2.2 | 4.0 |

The optimum has not moved: n = 1.2 and y₀ = 5 (equally 1.3 and 5, or 1.2
and 6, at 2.2 to 2.5 %) stay the best of the grid. The −1.3 % at the
smallest y⁺ is the discretisation error of the 64 × 48 mesh against the
finest one, not the treatment. The constants are unchanged.

## Channel, Re_τ = 395

Reference: U_b⁺ = 17.55 (DNS). Converged to the scaled residual 10⁻⁹
(PR #128; the table of levels 0 to 2 first given here, from runs stopped
at 10⁻⁶, was off by up to 0.9 % on level 2; the first level-3 runs by
up to 5.6 %), levels 0 to 3; p, GCI and
the extrapolated value from levels 1 to 3:

| model | cells | U_b⁺ | p | GCI_fine | extrapolated | finest vs DNS | iterations, converged |
|---|---|---|---:|---:|---:|---:|---|
| k-ω SST, low-re | 192 / 384 / 768 / 1536 | 17.426 / 17.367 / 17.333 / 17.317 | 1.01 (monotonic) | 0.12 % | 17.300 | −1.33 % | 4038, 8061, 20819, 75801 (SIMPLEC) |
| k-ε, scalable | 24 / 48 / 96 / 192 | 18.307 / 18.330 / 18.446 / 18.542 | 0.28 (monotonic) | 2.98 % | 18.985 | +5.65 % | 398, 401, 402, 427 |
| realizable k-ε, scalable | 24 / 48 / 96 / 192 | 18.691 / 19.100 / 19.477 / 19.614 | 1.46 (monotonic) | 0.50 % | 19.692 | +11.76 % | 602, 676, 691, 692 |
| k-ε-ζ-f, compound | 192 / 384 / 768 / 1536 | 17.869 / 17.899 / 17.910 / 17.911 | 2.77 (monotonic) | < 0.01 % | 17.912 | +2.06 % | 658, 772, 1203, 2531 |
| k-ε-ζ-f₀, compound-2023 | 96 / 192 / 384 / 768 | 17.339 / 17.543 / 17.576 / 17.583 | 2.16 (monotonic) | 0.02 % | 17.586 | +0.19 % | 436, 457, 535, 883 |
| elliptic-blending k-ε, compound | 192 / 384 / 768 / 1536 | 17.639 / 17.695 / 17.710 / 17.713 | 2.05 (monotonic) | < 0.01 % | 17.714 | +0.93 % | 5506, 9215, 24771, 90873 (SIMPLEC) |

Since PR #141 the study converges the channel with the coupled algorithm
(D178), within 3 000 iterations; the k-ε rows are those of SIMPLEC to
the digits shown. SST and EB keep their SIMPLEC values of PR #128: with
the coupled algorithm SST on 1 536 cells across and EB on 768 and 1 536
stop at 3 000 iterations with a scaled residual near 10⁻⁷ (U_b⁺ 17.340,
17.756 and 17.731, within 0.3 % of the converged values; a pseudo-time
factor c_T of 0.2 instead of 0.5 does not help; finding). The ζ-f and
ζ-f₀ rows are those of the wall sink of k (specification 3.3.2); with
the sink of 6d2 to #140 they read 17.886 / 17.929 / 17.940 / 18.029
(ζ-f, level 3 unconverged at 160 000 iterations) and 17.789 / 17.839 /
18.108 / 18.388 (ζ-f₀, divergent, GCI 51 %).

\* the run stopped at its iteration limit before the tolerance (`converged` false).

**Levels 1 to 3 (PR #123).** The first run of the study on level 3 (the
maintainer's Mac, 2026-11-20) failed in two k-ε-ζ-f runs with GEN001 from
the multigrid setup of the f equation. On the finest channel the
diagonals of the coarsest level span 10¹⁴, and the check compared each
pivot with the largest diagonal (a pivot of 2.7·10⁻⁴, 99 % of its own
diagonal, was refused). For the turbulence equations each pivot is now compared with its own
diagonal (the pressure keeps the strict test, which finds a matrix that
stops being positive definite early; a singular matrix is refused at any
scale), and a transport or
elliptic equation whose multigrid setup fails takes the Jacobi
preconditioner, with a note in the log, instead of stopping the run. With
that, both level-3 runs diverge with the default relaxations: the plate in
outer iteration 51 and the channel in iteration 1 182 (after its residuals
had stalled near 3·10⁻⁸ for U and 3.5·10⁻⁴ for ε). With the relaxations 0.5
for momentum and turbulence neither run diverged in the checks made in the
cloud (the plate to iteration 800, the channel past iteration 5 000, where
its residuals stand near 10⁻⁵ for U and 7·10⁻⁴ for ε, unconverged as on
level 2); the study sets them for these two runs, as it does for the k-ε
plate from level 2.

**Levels 1 to 3 on the maintainer's Mac (2026-11-20).** The study ran
with `--levels 1,2,3`; the runs that finished are below (the ζ-f plate
and channel, the ζ-f₀ and EB plates and the realizable k-ε plate on
level 3 were stopped as too slow on that machine and are missing).

| plate | cells | y⁺ (finest) | C_f ·10³ | p | GCI_fine | extrap. ·10³ | finest vs White |
|---|---|---:|---|---:|---:|---:|---:|
| k-ω SST, low-re | 3072 / 12288 / 49152 | 0.09 | 2.6682 / 2.6868 / 2.7007 | 0.43 (monotonic) | 1.87 % | 2.7410 | −6.05 % |
| k-ε, scalable | 12288 / 49152 / 196608 | 3.62 | 2.8270 / 2.8156 / 2.8059 | 0.24 (monotonic) | 2.39 % | 2.7524 | −2.39 % |

C_D of SST: 2.8259 / 2.8477 / 2.8606 ·10⁻³, p 0.74, GCI 0.84 %; of k-ε:
2.6334 / 2.6768 / 2.7459 ·10⁻³, p 0.67, GCI 5.31 %.

The channel values of that run (SST 17.381 / 17.346 / 17.009, EB
17.675 / 17.866 / 16.720, ζ-f₀ 17.840 / 18.106 / 18.378 on levels 1 to 3,
reproduced exactly in the cloud) are replaced by the converged table
above; the cause of their level-3 drop is the next paragraph.

What the plate shows: the SST plate on the three finest meshes meets the
2 % GCI band (1.87 %) and stays within 1.5 % of the Coles–Fernholz fit to
measured layers, while it is 6 % below White's correlation (outside ±5 %,
inside the widened ±8 %). The k-ε plate is just
outside the band (2.39 %). The SST plate converged to 10⁻¹¹ on levels 1
and 2 gives C_f 2.66823 and 2.68680 ·10⁻³, the values of the study at
its case tolerance 10⁻⁷ to five digits: the plate results stand.

### The level-3 channel drop (PR #128)

On level 3 of the first run, U_b⁺ of SST fell from 17.346 to 17.009 and
of EB from 17.866 to 16.720. These runs had not converged: the scaled
residuals (the largest of momentum, continuity, k, ω or ε and ζ or φ)
pass 10⁻⁶ while U_b⁺ is still moving. On 768 cells across, SST passes
17.009 at iteration 28 000 with the momentum residual at 2·10⁻⁷ and k
near 10⁻⁶, and then rises to 17.316 by iteration 95 000; EB passes its
first 10⁻⁶ at 17.95 and settles at 17.713. The first-cell y⁺ of level 3
is 0.034 (the cell centre; not 0.005 as first noted). The
mesh, the wall treatment and the y⁺ regime are not the cause: converged,
both models refine monotonically, with SST at first order as on the
plate and EB at p = 2.

U_b⁺ when the largest scaled residual first falls below each tolerance,
against the value converged to 10⁻¹² on the same mesh (y-refined meshes,
two cells in x):

| run (converged U_b⁺) | 10⁻⁶ | 10⁻⁷ | 10⁻⁸ | 10⁻⁹ |
|---|---:|---:|---:|---:|
| SST, 384 cells across (17.333) | +0.10 % | −0.01 % | −0.02 % | 0.00 % |
| SST, 768 (17.316) | −2.94 % | +0.18 % | −0.01 % | 0.00 % |
| EB, 384 (17.710) | +1.52 % | +0.03 % | −0.01 % | 0.00 % |
| EB, 768 (17.713) | +1.36 % | −0.24 % | +0.04 % | 0.00 % |
| ζ-f₀, 384 (18.388) | +0.94 % | +0.04 % | −0.01 % | 0.00 % |

(the U_b⁺ there read at the nearest multiple of 1 000 iterations). The
ω residual of SST stays between 10⁻¹¹ and 10⁻¹³ throughout: its scale
‖a_P ω‖ is dominated by the fixed wall-cell values of ω (10⁷ on level
3), so it never decides the stop. The study now converges the channel to
10⁻⁹ (iteration limits 20 000 to 160 000; SST needs 76 000 iterations on
768 cells across, EB 91 000) and refines it in y only: `nested_ratio`
gives the grading ratio whose nodes on 2^L times the base intervals are
those of level L, and the two cells in x stay. The x refinement adds only
cost: the 12 288-cell SST run of level 3 (16 cells in x), continued from
its stop at 17.009, reaches 17.316 at iteration 94 000 (largest residual
5·10⁻⁹), the value of the run with two cells in x.

**ζ-f₀ was not grid-convergent below y⁺ 1** (finding of #128, resolved
by #141, see below). Converged, U_b⁺
rises with every level (17.79, 17.84, 18.11, 18.39), with the first-cell
y⁺ from 1 to 0.11. The k of the wall cell does not vanish like y²: k⁺
is 1.15, 0.39 and 0.09 in the first cell on levels 0, 1 and 3 (the DNS,
k⁺ ≈ 0.08 y⁺², has about 0.08, 0.02 and 0.001 there), so the wall value ε = 2νk_P/y_P²
grows like 1/y_P (ε⁺ 14 on level 3, the DNS 0.17). The sink of k of the
6d2 default, k/T_K with the Kolmogorov bound T_K ≥ C_T (ν/ε)^½, is far
below ε in the viscous sublayer, and the compound-2023 treatment gives k
zero flux through the wall, so the sublayer keeps its k. Two checks
(experimental builds, not merged): the sink ε in the wall cell only
gives 15.70, 15.90, 16.27, 16.68 (still rising); the paper's sink ε
everywhere gives 15.686, 15.621, 15.621, 15.623, grid-converged but
11 % below the DNS (the 6d2 finding). A formulation that is both
grid-converged and accurate was left to the maintainer.

### The wall sink of k for both ζ-f models (PR #141)

Maintainer decision 2026-11-20: the sink of k is max(k/T_K, 2νk/d²), d
the wall distance (specification 3.3.2), for ζ-f and ζ-f₀, and ζ-f₀
takes C_ε1 = 1.4 (1 + 0.009/ζ). The sink is ε at the wall, so k falls
like y² and ε_w converges; away from the wall it is k/T_K, which keeps k
positive at a laminar leading edge. The candidates, all converged (the
channel to 10⁻⁹ on levels 0 to 3, the plate on 64 × 48 and 128 × 96,
the step on 8 200 cells, the jet on the three meshes of 6g):

| sink of k, C_ε1 of ζ-f₀ | channel U_b⁺, levels 0 to 3 | plate C_f vs White | step x_r/H | jet Nu₀ |
|---|---|---|---|---|
| k/T_K, 1.4 (1 + 0.012/√ζ) (6d2) | 17.79 / 17.84 / 18.11 / 18.39 | −2.5 / −2.9 % | 6.385 | 166 / 157 / 152 (GCI 7.2 %) |
| −ε, as printed | 15.69 / 15.62 / 15.62 / 15.62 | floor holds k | – | – |
| −ε, or a blend of −ε and k/T_K by ζ | – | leading-edge k at the floor, the ζ residual cycles | – | – |
| max(k/T_K, min(ε, 2νk/d²)), 0.009/ζ | grid-converged | converged | – | limit cycle (period about 80 iterations) at r/D 5.4 to 6 |
| max(k/T_K, 2νk/d²), 1.4 (1 + 0.038/√ζ) | – | – | 7.18 | – |
| **max(k/T_K, 2νk/d²), 1.4 (1 + 0.009/ζ)** (default) | **17.34 / 17.54 / 17.58 / 17.58** | **+1.8 / +0.2 %** | **6.44** | **198 / 196 / 196 (GCI 0.01 %)** |

The new default of ζ-f₀: channel +0.2 % from the DNS on the finest
mesh, p = 2.2; at first-cell y⁺ 10 / 30 −3.4 / −1.2 % (the channel) and
−0.3 / −2.7 % (the plate, from White); step x_r/H 6.44 (+3 %, coupled
and SIMPLEC alike). The impinging jet converges with the mesh but at
Nu₀ = 196, +38 % from the experiment's 141 (finding): the k of the
stagnation point. It comes from the sink (the 6d2 sink with the new
C_ε1 gives 159); a = 0.6 lowers it to 176 only. ζ-f₀ has no production
limiter at stagnation points; one is a follow-up. ζ-f (2004) with the
new sink converges on the plate (C_f +2.0 / +0.03 % from White on
64 × 48 / 128 × 96; before, the k floor of the leading edge stopped it
at the iteration limit, +16 to +24 %) and in the channel on all four
levels (+2.1 % on the finest).

**Findings.**

1. **Few models reach the asymptotic range on all three plate meshes.**
   C_f is oscillatory for SST and EB, divergent for k-ε, and monotonic for
   the realizable k-ε (p = 1.0) and, since PR #141, ζ-f and ζ-f₀. The spread of SST is
   small, 2.668 to 2.694 (GCI 0.47 %), so the order test is dominated by differences below
   1 %. SST and EB meet the 2 % GCI band of 9.3 without an order. Since
   PR #141 ζ-f and ζ-f₀ refine monotonically (p 2.2 and 1.9, GCI 0.15 and
   0.20 %, within 0.4 % of White on the finest mesh); before, the k floor
   of the leading edge kept the ζ-f plate from converging (+16 to +24 %).
2. **The k-ε family with scalable wall functions** refines its first cell
   to y⁺ 7 on the finest plate, below the floor y* = 11.06 of the
   treatment; C_f changes by 0.4 % (k-ε) and 3.6 % (realizable) over the
   three meshes, the insensitivity the treatment promises. The realizable
   model stays 10 % below White on the finest mesh (the 3 to 7 % of 6c on
   one mesh, growing with refinement).
3. **Channel** (converged, PR #128). SST (0.12 %, p = 1.0) and EB
   (below 0.01 %, p = 2.1) meet the 0.5 % GCI band on levels 1 to 3, and
   ζ-f on levels 0 to 2 (0.03 %; its level 3 does not converge with the
   relaxations 0.5); all three are within ±3 % of the DNS on the finest
   mesh (SST −1.3 %, EB +0.9 %). Since PR #141 ζ-f (p = 2.8, GCI below
   0.01 %, +2.1 %) and ζ-f₀ (p = 2.2, GCI 0.02 %, +0.2 %) converge with
   the mesh on levels 1 to 3 as well.
   The k-ε family with scalable wall functions on meshes refined towards
   the wall (24 to 192 cells) moves with every refinement: realizable
   +11.8 %, k-ε +5.7 % on the finest mesh, outside their ±5 % band.
4. **Still open**: the coupled algorithm for SST and EB on the finest
   channel meshes (residual stalls near 10⁻⁷); a stagnation-point
   limiter for ζ-f₀ (the jet's Nu₀ +38 %); the order of the plate (the plate's values are converged); a
   stop test on the monitored quantities (U_b⁺, forces), which would have
   caught the early stops of the channel without a tolerance per case.

## Statistics of transient runs (D71, PR #110)

`[output.statistics]` takes the time moments of the fields, probes and
force coefficients with weighted Welford updates (West 1979; each step
weighted with its time step) and the confidence interval of each mean from
batch means (Schmeiser 1982), with the von Neumann test of the lag-1
correlation of the batch means (Fishman 1978) as the test of the window;
with `phase_period` the batches hold whole periods.

| Check | Quantity | Band | Result |
|---|---|---|---|
| Weighted samples with unequal weights and a mean of 10⁶ (`weighted_moments_match_the_two_pass_values`) | mean, RMS, covariance against the two-pass values | 1e-9 relative | asserted |
| White noise, 20 000 samples (`batch_means_of_noise`) | the 95 % interval holds the true mean; half width against 1.96 σ/√N; verdict | factor 2; `sufficient` | asserted |
| A sine over whole periods; a quarter period (`batch_means_of_a_sine`) | mean 0, RMS a/√2, phase averages at the bin centres; verdict | 1e-3; 7 % of a; `sufficient` / `too short` | asserted |
| A sine over 9.3 periods in 10 batches; 20 periods (`batches_hold_whole_periods`) | batches of whole periods: 9 batches of one, 10 of two; verdict | `too short` without the period, `sufficient` with it | asserted |
| A drag-like series, 3.2 + 0.02 sin + a drift of 10⁻⁵ (`tolerance_passes_a_drift_far_below_the_size`) | verdict | `too short` with `tolerance = 0`, `sufficient` with 10⁻³ | asserted |
| The w of a probe in a 2D case, round-off (`round_off_component_of_a_probe_is_not_too_short`) | verdict | `sufficient` | asserted |
| Decaying Taylor–Green vortex, 32² cells, Δt 0.005, t 0.1 to 0.4 (`taylor_green_statistics_meet_the_exact_time_moments`): u = u₀(x) f(t), so the mean is u₀⟨f⟩, the RMS \|u₀\|σ_f, the covariance u₀v₀σ_f² | L2 of the mean / RMS / covariance fields against the exact moments | 1 % / 2 % / 3 % | 0.60 % / 1.02 % / 1.71 % |
| | probes a (0.125, 0.375) and b (0.3, 0.1): mean and RMS of u and v | 1 % of \|u₀\| | inside |
| | verdict of the decay (a trend) | `too short`, STAT001 | asserted |
| Cylinder at Re = 100, t = 5 to 12 s, 10 batches of two periods (`cylinder_re100_statistics.toml`, gallery 110) | every series of the probes and coefficients; the benchmark | `sufficient`; inside its ranges | all 14 `sufficient` (C_D 3.20207 ± 1.0e-4, C_L −0.0104 ± 5.4e-5; the drag and three probe series by the tolerance, r1 0.51 to 0.72); C_D,max 3.232, St 0.2955 |
| Steady case; `batches = 1`; `phase_period = 0` (`statistics_keys_are_checked`) | refused | CASE005; CASE004 | asserted |

Found while making gallery 110: with ten batches over nine shedding periods
each batch held part of a period, and the batch means of the lift carried
the phase from batch to batch (r1 0.71, `too short`, although the mean was
converged). With `phase_period` the batches now hold whole periods. The
plain test of the batch means then still flagged the drag (it drifts by
10⁻⁴ of its value over the window), so `tolerance` passes a mean whose
widened interval is below 0.1 % of its size.


## Curvature and rotation correction of SST (6h, D170, PR #113)

`turbulence.curvature_correction` multiplies the production of k and ω by
f_r1 = max[0, min(f_rotation, 1.25)] (Spalart & Shur 1997; Smirnov &
Menter 2009). The kernel `turb::sst::rotation_function` takes the
velocity gradient, DS_ij/Dt, Ω^rot of the frame and the spin of a relative
velocity. The flow builds DS_ij/Dt from the least-squares gradients of the
six components of S_ij (zero normal gradient on every boundary), the
velocity relative to an `[[mrf]]` zone, and BDF1 from the last time level
in a transient run. Ω_ij is the rotation tensor of the absolute velocity.

| Check (9.8) | Quantity | Band | Result |
|---|---|---|---|
| Channel without rotation, with and without the correction (`tests/curvature.rs`, 48 cells; the gallery channel, 192 cells) | driving force | ≤ 10⁻⁶ relative | equal to the last digit on both (r* = 1, r̃ = 0 exactly in a periodic shear flow) |
| f_r1 of a rotating shear flow by hand (kernel test) | f_r1 | 10⁻¹² | met |
| The same flow in the inertial frame (absolute gradient) and in the turning frame (relative gradient, spin) | f_r1 | 10⁻¹⁴ | met: the sign of the frame term of Ω_ij is the consistent one |
| Channel in spanwise rotation, Re_b 5 800, Ro = 2Ωh/U_b = 0.15 (`channel_rotating_sst.toml`, 192 cells) | τ_w lower/upper | DNS of Kristoffersen & Andersson (1993) ±15 % | SST 1.000 (no rotation effect, as expected for a linear model), corrected 1.627, unstable side the lower wall. **Open:** the DNS value (the paper's data are not in the repository) |
| Backward-facing step, 44 160 cells (level 1 of `backstep_validation.toml`) | x_r/H | reported | 6.743 SST, 6.869 corrected (+1.9 %; experiment 6.26): the correction damps the turbulence in the curved shear layer by up to 10 % and the bubble grows |
| Manufactured solution with the correction (PR #184) | observed order | 1.8–2.2 | 1.92 (hex), 1.99 (poly); no convergence on tetrahedra ([details](#manufactured-solutions-of-k-ε-ζ-f₀-the-curvature-correction-and-the-qcr-pr-184)); since PR #197 2.0 (hex), 2.0–2.1 (poly), 1.9–2.1 on regular and perturbed tetrahedra ([details](#consistent-curvature-correction-on-irregular-meshes-d236-pr-197)) |

**Findings.** (1) In the irrotational free stream r* = S/Ω is large and
f_r1 sits at its upper limit 1.25. The production there is zero, so the
results do not change, but the derived field shows 1.25 outside the shear
layers. (2) On the step the correction moves the reattachment away from
the experiment (+7.7 % to +9.7 % of 6.26 on this mesh). This is the known
behaviour of SST with this correction behind steps. (3) Sources differ in
the sign of the frame term of Ω_ij. The one used here gives the absolute
rotation tensor. With it, r* and r̃ both mark the side of the rotating
channel whose mean vorticity opposes the rotation as unstable, and the
inertial and the turning frame give the same f_r1.

## Quadratic constitutive relation of SST (6h, D170, PR #115)

`turbulence.qcr` adds the non-linear stress of Spalart (2000),
−C_cr1 (O_ik τ_jk + O_jk τ_ik) with τ = 2 μ_t S* and O = 2W/|∇u|
(C_cr1 = 0.3), and for `"qcr2013"` the isotropic −C_cr2 μ_t √(2 S*:S*) δ_ij
(C_cr2 = 2.5) of Mani et al. (2013). The kernel `turb::sst::qcr_stress`
gives the cell stress. `ExplicitStressFaces` interpolates it to the faces
and subtracts τ_f·S from the deferred part b_f of the momentum faces.

| Check (9.8) | Quantity | Band | Result |
|---|---|---|---|
| Simple shear by hand (kernel test) | the six components | 10⁻¹² | met: only τ_xx and τ_yy, of equal size and opposite sign |
| Square duct, Re_b 3 500, linear SST (`tests/curvature.rs`, 12 × 12; the gallery duct, 72 × 72) | largest secondary velocity | 0 | below 10⁻¹⁰ U_b on both |
| The same with QCR2000 and QCR2013 (`tests/curvature.rs`) | largest secondary velocity | > 0 | met for both |
| Square duct, QCR2000, quarter duct 2 × 48 × 48 and 2 × 72 × 72 (`square_duct_qcr.toml`) | largest secondary velocity / U_b, its position | DNS of Pinelli et al. (2010) ±30 %, the vortex pair present | 1.55 % and 1.56 % on the diagonal at 0.35 h and 0.33 h from the corner, directed into the corner; eight vortices in the whole duct. **Open:** the DNS value (the data are not in the repository) |
| QCR2013 against QCR2000, 48 × 48 | velocity | the same | the same to the printed digits; the isotropic term only shifts the pressure |
| Driving force, 72 × 72 | relative to linear SST | reported | +1.9 % (1.1150 against 1.0941 × 10⁻² N/m³) |

**Findings.** (1) The isotropic 2013 term only changes the pressure in
incompressible flow. A model that carries k could use (2/3) ρk in its
place ("QCR2013-V"); with the modified pressure p̂ = p + (2/3) ρk that
form is also only a pressure shift here. (2) The secondary velocity is
small (about 1.5 % of U_b), but the in-plane residuals fall like the
others: the 72 × 72 duct converges in 2 624 iterations, against 2 664
for linear SST on the same mesh.

## Smooth limits against the hard form (spec 2.6, PR #136)

`turbulence.smooth_limits = true` replaces every min and max of the models
with the C¹ form of spec 2.6 (width w = 10⁻³). The first run of the check
found the SST flat plate with smooth limits not converging (k residual
cycling at 10⁻³ for 10 000 iterations): the smooth width followed the
larger operand, m = ½(|a| + |b|), so the floor max(CD_kω, 10⁻¹⁰) sat at
about w²|CD_kω|/16 where CD_kω < 0. The width now follows the smaller
operand, m = |a||b|/(|a| + |b|), in a form free of cancellation; the runs
without smooth limits are unchanged.

| Check | Quantity | Band | Result |
|---|---|---|---|
| TMR plate SST, level 1 (64 × 48), hard and smooth (`smooth_limits_change_the_plate_friction_below_a_tenth_percent`) | C_f at x = 0.97; drag of the plate | 0.1 % | 0.090 %; 0.093 % (both converge in 635 iterations) |
| The same on level 0 (32 × 24) | C_f | 0.1 % | 0.1002 % (reported, not asserted) |
| `smooth_switches_are_bounded_and_c1` (core) | error of smin, smax against min, max, incl. (±10³, 10⁻¹⁰) | w·m/2, m ≤ min(\|a\|, \|b\|) | asserted |

## Manufactured solutions of every model (spec 9.2, PR #148)

`crates/foehncraft-flow/src/turbulence/mms.rs` solves the turbulence
equations alone in a frozen cellular flow, u = U₀ (sin πx cos πy,
−cos πx sin πy) with U₀ = 0.25 m/s in the unit square, ν = 10⁻⁵ m²/s.
k, ω or ε, ζ (φ), f (α) and the wall distance d are smooth manufactured
fields (d = 0.06 + 0.6 y + 0.1 x lets F₁ of SST run from 1 to 0.01); the
sides are velocity inlets with the exact k and ω/ε, ζ uniform there and
f with the zero gradient of its exact field. The added sources are the
residuals of the continuous equations written out in the test from the
published equations, independently of the kernels, with derivatives from
fourth-order differences. No switch of a model changes branch in the
domain (`manufactured_fields_stay_on_one_branch`), so every term is
smooth and active: production, destruction, cross diffusion and the F₁
blending of SST, the Durbin bound of the standard k-ε (inactive), the
variable C_μ of the realizable model, the scales T and L and C_ε1(ζ) of
ζ-f, and E, C*_ε2, the cross diffusion and the α equation of the
elliptic-blending model. 16 × 16 to 64 × 64 cells, perturbed by 15 %;
the residuals to 10⁻¹¹ (150 to 660 iterations).

| Model | Mesh | p of k, ω/ε (, ζ/φ, f/α), finest pair | L2 error of k, finest (relative) |
|---|---|---|---|
| k-ω SST | hex / poly / tet | 2.06, 2.04 / 2.10, 2.13 / 2.11, 2.17 | 2.7e-5 / 4.4e-5 / 1.5e-5 |
| standard k-ε | hex / poly / tet | 2.04, 2.07 / 2.10, 2.11 / 2.15, 2.24 | 2.2e-5 / 2.9e-5 / 1.1e-5 |
| realizable k-ε | hex / poly / tet | 2.05, 2.10 / 2.06, 2.12 / 1.82, 2.02 | 2.2e-5 / 3.3e-5 / 2.6e-5 |
| k-ε-ζ-f (2004) | hex | 2.03, 2.07, 1.89, 2.00 | 2.3e-5 |
| | poly | 2.09, 2.11, 2.06, 2.07 | 3.0e-5 |
| | tet | 2.13, 2.22, 2.29, 2.13 | 1.2e-5 |
| elliptic-blending k-ε | hex | 2.03, 2.07, 1.94, 1.95 | 2.3e-5 |
| | poly | 2.09, 2.11, 2.04, 2.04 | 2.9e-5 |
| | tet | 2.18, 2.25, 2.20, 1.93 | 1.1e-5 |

Bands: 1.8–2.2 on hexahedra, 1.6–2.4 on polyhedra and triangular prisms;
all met. A change of 2 % in one constant of the written-out model (γ₂ of
SST, C_ε2 of k-ε) drops the order to 0.2–0.7, so the sources check the
implemented equations and not only the discretisation. The kernel
tests with dual numbers (the Jacobians of the cell functions) stay green.

**Findings.** (1) The realizable C_μ took W̃ = S_ij S_jk S_ki/S̃³ from
the strain with its discrete trace. W̃ and its range ±1/√6 belong to a
trace-free strain; where S → 0 (the centre lines of the cellular flow)
the discrete divergence, O(h²), dominated S̃ and W̃ jumped across its
range, so C_μ jumped in a band of cells: order 1.37 for k on the
triangular prisms. C_μ now takes the deviatoric strain, which leaves an
incompressible flow unchanged to the discretisation error and gives
order 1.82. (2) The f equation of ζ-f is solved in the divergence form
−∇·(L²∇f) + f; the paper writes L²∇²f. The two differ by ∇L²·∇f; the
manufactured source takes the solver's form (the EB α equation is
solved with L² outside the divergence, as its paper writes it). (3) With
no walls the model sets d = 10¹⁰; the manufactured d is set with
`Turbulence::set_wall_distance` after attaching. (4) The ζ-f₀ form, the
curvature correction and the QCR are covered since PR #184 (next
section).

```
cargo test --release -p foehncraft-flow --lib turbulence::mms -- --nocapture
```

## Manufactured solutions of k-ε-ζ-f₀, the curvature correction and the QCR (PR #184)

The harness of PR #148 (frozen flow, written-out sources, 16 × 16 to
64 × 64 cells, residuals to 10⁻¹¹) now covers the last three items of
stage 6.

**k-ε-ζ-f₀.** The source is written from Popovac (2023) with our choices
of 6d2 and PR #141: ν_t = (c_μ/ζ_eq) ζ k T, the realizability bound with
a = 0.6 √2, the sink of k max(k/T_K, 2νk/d²), C_ε1 = 1.4 (1 + 0.009/ζ), the
ζ equation f₀ − ζ [P_k/ε + (N_ζ − 1)]/T and the f₀ equation with
+ (N_ζ − 1) ζ/T on the right. The manufactured fields keep every switch on
one branch (the sink is k/T_K = ε, 2ν/d² is smaller; `manufactured_fields_stay_on_one_branch`
checks it for both ζ-f forms). Orders of k, ε, ζ, f₀, finest pair:
hex 2.03, 2.07, 1.85, 1.98; poly 2.09, 2.11, 2.06, 1.99; tet 2.13, 2.22,
2.27, 2.06 (bands 1.8–2.2 and 1.6–2.4, met; the tet values lie above 2.2
as for ζ-f in PR #148).

**Curvature correction.** The cellular flow of PR #148 has S₁₂ = 0 and
r̃ = 0 identically, so f_r1 would not exercise DS_ij/Dt. The test uses a
second prescribed flow, the stream function of the cellular flow with the
modes (2, 1) and (1, 2) added (β = 1), tangential to the sides, with
strain axes that turn along the streamlines. The source multiplies the
production of k and ω by f_r1 after the limiter, f_r1 written out from
Smirnov & Menter (2009) with DS_ij/Dt = u·∇S_ij from nested fourth-order
differences. f_r1 is clipped to [0, 1.25] on much of the domain (the
figure of gallery 184); the clip is continuous, so the source has kinks
but no jumps and second order is kept. |r̃| reaches 0.25 where f_r1 is
not clipped, so a wrong DS_ij/Dt changes the source by about 10 % and the
error would stop converging. Orders of k, ω, finest pair: hex 1.92, 1.95;
poly 1.99, 2.06.

**Finding: tetrahedra.** The error of k and ω is 1.3 % and 0.7 % at all
three levels (`curvature_correction_on_tetrahedra_stays_bounded` asserts it
stays below 3 %). DS_ij/Dt takes the least-squares gradient of the
least-squares gradient of the velocity; that second derivative is not
consistent on irregular tetrahedra, so the production of both equations
carries an error that does not fall with the mesh. Hexahedra and
polyhedra are consistent. Listed as a gap of
`turbulence.curvature-correction` in the grades; a consistent gradient of
S_ij on tetrahedra (a face-based or higher-order reconstruction) is a
follow-up and not done here. Fixed in PR #197 (D236):
[consistent curvature correction](#consistent-curvature-correction-on-irregular-meshes-d236-pr-197).

**QCR.** The QCR changes the momentum equations, not the turbulence
equations, so the test checks the stress and its divergence as the
momentum equations take them: on the manufactured k, ω (ν_t from the
solver) and the cellular flow, `Flow::explicit_stress_force` (new, the
face sums of `ExplicitStressFaces` per unit volume) is compared with
∂_j τ_ij of the stress written out from the definitions of Spalart (2000)
and Mani et al. (2013), C_cr1 = 0.3, C_cr2 = 2.5, over the cells away
from the sides and from the lines where S = 0 (|S| of the isotropic
stress of QCR2013 has a kink there). Order of the L2 force error, regular
meshes: QCR2000 hex 2.06, poly 2.06, tet 1.03; QCR2013 hex 2.06, poly 1.94,
tet 1.01 (bands 1.8–2.2; first order 0.8–1.3 on tetrahedra). The coupled
momentum solve with QCR is not in the manufactured solution.

**Finding: the explicit stress is only consistent on regular meshes.** The
divergence is Green-Gauss with the linear interpolation of the cell
stress to the face. On hexahedra with randomly perturbed vertices (15 %)
the error stays at 5.5 % of the largest force at all three levels
(`qcr_stress_on_perturbed_meshes_stays_bounded`), on tetrahedra it is
first order even without perturbation. A skewness correction of the face
value with the gradient of the stress would remove it; the secondary flow
of the square duct (gallery 115) is a structured-mesh result and is not
affected at the quoted accuracy. Flagged for the owner; not changed here. Fixed in PR #194 (D236), next section.

## Consistent QCR stress on irregular meshes (D236, PR #194)

The finding above is fixed in two parts. (1) On the cells where the
operator check (D236, PR #189) finds the face-neighbour velocity gradient
first order, the explicit stress takes the velocity gradient of a
quadratic least-squares fit over the neighbours of the neighbours (third
ring where fewer than 2m points are in reach, m the number of unknowns;
boundary faces as points with their values; normal equations scaled by
the cell size and solved by Cholesky with a pivot tolerance of 10⁻⁹). The
fit is exact for quadratic fields (`foehncraft-fv` exactness tests,
`operators::tests::quadratic_fit_is_exact_on_tetrahedra`). (2) The face
stress is the weighted mean of the two cells plus ∇τ̄_f·(x_f − x_w), with
least-squares cell gradients of the stress, x_w the point of the line
between the centroids on the face; at a boundary face τ_P + ∇τ_P·(x_f − x_P).

L2 error of the QCR2000 force, 15 % perturbed vertices, three levels
(`qcr_stress_on_perturbed_meshes_converges`, band: order 0.9 to 2.2 and a
falling error):

| Mesh | Before | After | Order, finest pair |
|---|---|---|---|
| hexahedra | 8.5 %, 5.7 %, 5.5 % | 4.6 %, 1.5 %, 0.65 % | 1.2 |
| tetrahedra | 27.8 %, 14.6 %, 10.3 % | 15.6 %, 7.3 %, 3.6 % | 1.0 |

Either part alone was not enough in the prototype: the skewness
correction with the face-neighbour gradient still left an error that did
not fall. Regular meshes keep their orders (`qcr_stress_manufactured_solution_orders`:
hexahedra and polyhedra 1.8 to 2.2, tetrahedra 0.8 to 1.3). The first
order that stays on tetrahedra is the midpoint rule of the face sum on
irregular cells, not the face stress. The coupled momentum solve with QCR
is still not in the manufactured solution.

## Consistent QCR stress on irregular meshes (D236, PR #194)

The finding above is fixed in two parts. (1) On the cells where the
operator check (D236, PR #189) finds the face-neighbour velocity gradient
first order, the explicit stress takes the velocity gradient of a
quadratic least-squares fit over the neighbours of the neighbours (third
ring where fewer than 2m points are in reach, m the number of unknowns;
boundary faces as points with their values; normal equations scaled by
the cell size and solved by Cholesky with a pivot tolerance of 10⁻⁹). The
fit is exact for quadratic fields (`foehncraft-fv` exactness tests,
`operators::tests::quadratic_fit_is_exact_on_tetrahedra`). (2) The face
stress is the weighted mean of the two cells plus ∇τ̄_f·(x_f − x_w), with
least-squares cell gradients of the stress, x_w the point of the line
between the centroids on the face; at a boundary face τ_P + ∇τ_P·(x_f − x_P).

L2 error of the QCR2000 force, 15 % perturbed vertices, three levels
(`qcr_stress_on_perturbed_meshes_converges`, band: order 0.9 to 2.2 and a
falling error):

| Mesh | Before | After | Order, finest pair |
|---|---|---|---|
| hexahedra | 8.5 %, 5.7 %, 5.5 % | 4.6 %, 1.5 %, 0.65 % | 1.2 |
| tetrahedra | 27.8 %, 14.6 %, 10.3 % | 15.6 %, 7.3 %, 3.6 % | 1.0 |

Either part alone was not enough in the prototype: the skewness
correction with the face-neighbour gradient still left an error that did
not fall. Regular meshes keep their orders (`qcr_stress_manufactured_solution_orders`:
hexahedra and polyhedra 1.8 to 2.2, tetrahedra 0.8 to 1.3). The first
order that stays on tetrahedra is the midpoint rule of the face sum on
irregular cells, not the face stress. The coupled momentum solve with QCR
is still not in the manufactured solution.

Figure: [mms_turbulence.png](../figures/gallery/mms_turbulence.png) (gallery 184).

```
cargo test --release -p foehncraft-flow --lib turbulence::mms -- --nocapture
```

## Consistent curvature correction on irregular meshes (D236, PR #197)

The finding on tetrahedra of PR #184 is fixed with two quadratic
least-squares fits. (1) S_ij comes from a velocity gradient that is second
order on every cell: the fit of PR #194 on the cells where the operator
check finds the face-neighbour gradient first order. (2) On the cells
where it finds the gradient of the gradient inexact for a quadratic field,
the gradient of S_ij in DS_ij/Dt = u·∇S_ij is a quadratic least-squares
fit of S_ij over interior cells only (neighbours of the neighbours, a third
and fourth ring at corners), so no boundary values of S_ij are needed;
elsewhere it stays the least-squares gradient with zero normal gradient
at the boundaries. (2) is exact for the linear S_ij of a quadratic
velocity (`foehncraft-fv` exactness tests,
`operators::tests::quadratic_fit_is_exact_on_tetrahedra`).

A first version took ∇S_ij from the Hessian of the velocity fit, ½ (∂_k∂_j u_i +
∂_k∂_i u_j). It converged in the manufactured solution, but in a channel
without rotation it moved the driving force by 2 % (`tests/curvature.rs::pure_shear_is_unchanged`,
band 10⁻⁶): the second-order coefficients of a fit of the steep wall
profile U(y) pick up its non-polynomial part, so ∂_x S_ij ≠ 0 in pure
shear. Fitting S_ij itself keeps ∂_x S_ij = 0 exactly wherever the mesh is
symmetric in x, as before.

L2 error of k and ω in the manufactured solution of PR #184
(`curvature_correction_on_tetrahedra_converges`, band: order 1.6 to 2.4 and
a falling error):

| Mesh | h | k | ω |
|---|---|---|---|
| tetrahedra | 9.88e-3, 4.94e-3, 2.47e-3 | 5.05e-3, 1.38e-3, 3.61e-4 | 2.34e-3, 5.94e-4, 1.49e-4 |
| tetrahedra, 15 % perturbed | 9.88e-3, 4.94e-3, 2.47e-3 | 5.42e-3, 1.41e-3, 3.59e-4 | 2.47e-3, 6.26e-4, 1.51e-4 |

Orders of the finest pair: 1.93 and 2.00 (regular), 1.98 and 2.05
(perturbed); before, the error of k stayed at about 1.3 % at every level.
Hexahedra and polyhedra (`curvature_correction_manufactured_solution_orders`)
give 2.02, 2.03 and 2.03, 2.07 (before 1.92, 1.95 and 1.99, 2.06): their
wall cells take the fits too.
