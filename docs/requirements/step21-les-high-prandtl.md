# Step 21 – Scale-resolving numerics (LES) and high-Prandtl heat transfer: assessment

Requirements step 21 (2026-11-17). The maintainer wrote a prompt elsewhere
for "scale-resolving numerics (LES) and high-Prandtl-number heat transfer"
after stage 6b and asked for a critical review against the state of the art,
checked against the current code and requirements, with the changes to our
path named. This page is that review. The maintainer approved its direction
on 2026-11-17 ("update the requirements with the assessment only where
necessary") and decided at the same time to **drop PANS altogether**.

Decided from this page (section 8): D172–D175 (interfaces now), L39–L43
(LES after v1.0), P16–P17 (research items), X2 (PANS), O20 (papers and data
from the maintainer). No solver code is written before the maintainer's go
(project rule). Clean-room rule (D21): methods from the cited papers only;
other codes are not named.

Contents

1. [Short answer](#1-short-answer)
2. [What Foehncraft has today](#2-what-foehncraft-has-today)
3. [Low-dissipation numerics](#3-low-dissipation-numerics)
4. [Subgrid-scale models](#4-subgrid-scale-models)
5. [What the prompt leaves out](#5-what-the-prompt-leaves-out)
6. [High-Prandtl heat transfer](#6-high-prandtl-heat-transfer)
7. [Verification suite and cost](#7-verification-suite-and-cost)
8. [Items and their stages](#8-items-and-their-stages)
9. [Questions for the maintainer](#9-questions-for-the-maintainer)
10. [References](#10-references)

## 1. Short answer

The prompt points the right way (energy-conserving central convection,
algebraic subgrid models that need no test filter, a pluggable turbulent
Prandtl number, semi-local wall scaling) but it should not be run as written:

- **Four technical errors.** The Ducros sensor never switches in
  incompressible flow; "Choi / Mahesh–Ham" names two different things, and
  the first is already in the code; Vreman's model does not vanish in
  solid-body rotation and scales as y, not y³, at a wall; the coherent
  structure model is not dynamic in the Germano sense. The Taylor–Green
  reference should be the 512³ pseudo-spectral solution, not Brachet 1983.
- **Gaps.** Time integration decides whether "kinetic-energy preserving"
  holds at all; WALE is missing as the common baseline; wall-modelled LES
  is named but not specified, and for industrial internal flows it matters
  more than wall-resolved LES; averaging over homogeneous directions and
  spectra are missing; the inlet should be the divergence-free SEM.
- **High Pr.** The TU Graz work is a good anchor. The P-function with a
  modified Kays–Crawford Pr_t, in its variable-property form (Steiner &
  Irrenfried 2019; papers and thesis provided by the maintainer), is now
  specified and checked against the published DNS Nusselt numbers (D172,
  D173); the thesis tables (CC BY 4.0) give the reference data, only the
  DNS profiles are not public (O20). A gradient (Clark)
  SGS heat flux on its own is unstable and becomes a research item (P16).
- **Split in three.** High-Pr closures go into the interfaces of stage 7
  (D172, D173). Two small interface items keep LES open without work now
  (D174, D175). LES numerics, SGS models, WMLES, inlet and statistics
  follow v1.0 under L5 (L39–L43).

## 2. What Foehncraft has today

Checked against `main` at 6c90e7c (stage 6b merged):

| Topic | Today | Source |
|---|---|---|
| Convection | Linear upwind as deferred correction (D25); `solver.scheme = "central"` is a **bounded** central scheme: distance-weighted interpolation plus a gradient correction along the face (`w u_P + (1−w) u_N + ∇u_f·r_w`), bounded in the normalised variable (β = 0.1) and scaled by the upwind limiter; one blending number for all faces | `foehncraft-kernels/src/flow.rs` `momentum_face` |
| Rhie–Chow | Time-step-independent form (D27) with the time-step term γ_f Y_f; PISO solves for the **pressure correction** p′ with the latest pressure in the predictor (incremental form) | `flow.rs` `face_mass_flux`, `foehncraft-flow/src/transient.rs` |
| Time | BDF1, BDF2 (variable steps), Crank–Nicolson with θ = ½ and the pressure at n+½, first step BDF1; outer correctors per step optional (D28 (2)); CN refused for turbulence variables (stage 6 spec 6.1) | `transient.rs` |
| Energy ledger | Kinetic-energy budget per step with a **pressure part** (collocated −∫p∇·u) and a **scheme part** (convection and time scheme), pairwise sums, bitwise reproducible | `foehncraft-flow/src/balance.rs` |
| Taylor–Green | 2D cases (BDF2, CN, CFL 10, space order); central scheme: order 2, a tenth of the numerical dissipation of linear upwind (D25, part 5a) | `examples/cases/flow/taylor_green_*.toml` |
| Turbulence | `enum TurbulenceModel { Laminar, KOmegaSst }`, host trait per model, kernels monomorphised; invariants module (D166) | `foehncraft-flow/src/turbulence/mod.rs`, `foehncraft-kernels/src/turb/` |
| Gradients | Weighted least squares, exact for linear fields on every cell type (D24) | `foehncraft-kernels/src/fv.rs` |
| Pr_t | Specified as `constant` (0.85) or `kays-crawford` per scalar (D135, stage 6 spec 5.1); SGDH k_eff = k + μ_t c_p/Pr_t (D152); not yet in code | stage 6 and stage 7 pages |
| Thermal wall | Specified: Kader (1981) default, Jayatilleke option (D135) | stage 6 spec 5.2 |
| Statistics | Running means, rms, covariances, batch-means confidence interval (D71), stage 6 PR 6g | requirements |
| RNG | Counter-based streams keyed by entity id (D122) | requirements |
| Periodic boundaries | Driven by mass flow or pressure gradient, heat-flux scalars (6p, PR #67) | `foehncraft-mesh/src/periodic.rs` |

Two consequences. First, the present central scheme is **not** kinetic-energy
preserving on any mesh: the distance weights, the gradient correction and the
bounding all break the skew-symmetry (section 3.1). Second, the ledger already
measures exactly the quantity the prompt cares about, split into its two
sources, so every LES numerics change can be judged by a number rather than
by eye.

## 3. Low-dissipation numerics

### 3.1 Kinetic-energy-preserving convection: keep, with the exact conditions

With face mass fluxes ṁ_f that satisfy discrete continuity exactly
(Σ_f ṁ_f = 0 per cell) and the **arithmetic** face value
u_f = ½(u_P + u_N), the convective term C(u)_P = Σ_f ṁ_f u_f is
skew-symmetric: Σ_P u_P·C(u)_P = 0 on any polyhedral mesh, whatever the cell
shapes (Verstappen & Veldman 2003; Mahesh, Constantinescu & Moin 2004;
Felten & Lund 2006). The prompt is right on this point. The conditions that
make it true in the code:

1. **Weights independent of the mesh.** Felten & Lund (2006) show that the
   distance-weighted ("second-order") interpolation gives a kinetic-energy
   error and severe oscillations on curvilinear meshes, while the
   mesh-independent ½ removes the interpolation error. Our gradient
   correction ∇u_f·r_w also breaks the symmetry. A `kep` scheme therefore
   uses ½(u_P + u_N) and no correction.
2. **Divergence-free fluxes at the time level that convects.** Exact only
   if the pressure equation is converged tightly in every corrector; with
   loose tolerances the residual continuity error enters the energy budget.
   The −φ_P∇·F term that makes the scalar operators M-matrices (stage 6
   spec 6.1) must not be added to the `kep` momentum operator.
3. **A time scheme that conserves quadratic invariants** (section 3.4).
4. **No limiter and no bounding.** Bounding and limiting are what make the
   present central scheme robust; `kep` gives that up and needs the local
   stabilisation of section 3.2.

Accuracy: the ½-weights lower the local truncation error on stretched
meshes; the symmetry-preserving literature reports that the global error
still converges at second order on smoothly varying meshes (Verstappen &
Veldman 2003). We measure it with the existing space-order cases on
perturbed polyhedra before `kep` is offered for production.

### 3.2 Local stabilisation: not the Ducros sensor

Ducros et al. (1999) blend by Φ = (∇·u)² / ((∇·u)² + |ω|² + ε). In
incompressible flow ∇·u = 0, so Φ = 0 everywhere and the sensor never
switches (on the collocated grid the cell divergence is a small discretisation
residue, which makes Φ noise, not a signal). The sensor was built for shocks.

What works for incompressible LES on unstructured meshes, filled into the
per-face blending field of D174:

- **Mesh-quality blending.** A fixed small upwind (linear-upwind) fraction
  on faces with large skewness or non-orthogonality, from the quality
  fields of D50. This targets the source the prompt names (polyhedral
  skewness) directly and costs nothing at run time.
- **Wiggle detector.** Upwinding switched on only where the velocity shows
  odd–even oscillations, detected from sign changes of successive
  differences along P–N (Mary & Sagaut 2002 used such a detector for LES of
  a stalled airfoil). On polyhedra the "successive differences" are the
  face difference and the two reconstructed differences ∇u_P·d and
  ∇u_N·d.
- **The hybrid blending of DES** (Travin et al. 2000) for L5, which needs the
  same per-face field.

Every blend enters the ledger's scheme part, so the numerical dissipation it
adds is reported per run.

### 3.3 Rhie–Chow: the audit, not a new algorithm

The prompt's sentence mixes two different things:

- **Choi (1999)** made the momentum interpolation independent of the time
  step in a steady state. That is D27; it is in the code (the γ_f Y_f term
  of `face_mass_flux`).
- **Mahesh, Constantinescu & Moin (2004)** (with the energy discussion of
  Ham & Iaccarino 2004) is a different algorithm: a fractional-step method on
  a collocated grid in which the face-normal velocities are carried as their
  own variables, projected with the compact face pressure gradient, and the
  cell velocities are corrected with a least-squares reconstruction of the
  face pressure gradients chosen to keep the energy error small. It is not
  a coefficient in the Rhie–Chow formula.

The prompt's worry also points the wrong way. With the time term in a_P,
D_f = V/a_P behaves like Δt/ρ as Δt → 0, so the pressure-smoothing term
**shrinks** at small LES time steps; it does not over-damp. The kinetic-energy
error that the collocated pressure coupling leaves is of order Δt h² when the
total pressure is solved for and of order Δt² h² in the pressure-increment
form (Felten & Lund 2006; Komen et al. 2021, who reach practically full
conservation with the increment form in a PISO-type method). Our PISO is in
increment form already.

Decision: no change to D27 now. With L39 the Taylor–Green runs of section 7
measure the pressure part of the ledger against Δt and h; only if it does
not fall as Δt² h² does the formulation change, and then the candidate is the
energy-consistent cell-velocity correction of Mahesh et al., not a new
Rhie–Chow coefficient.

### 3.4 Time integration: missing from the prompt, decides the outcome

A skew-symmetric space operator conserves kinetic energy only together with
a time scheme that conserves quadratic invariants:

- **Implicit midpoint** conserves them exactly; **Crank–Nicolson** does so
  when the convecting flux is the one of t^{n+½} and each step is converged
  (outer correctors of D28 (2)); with one PISO pass and the extrapolated
  flux the step is not converged and a small error of either sign stays.
- **BDF2** damps each resolved mode by a factor of the fourth order in ωΔt
  per step. At CFL below 0.5 this mostly hits the scales near the grid
  cut-off, but it is a dissipation the SGS model does not know about; the
  ledger's scheme part measures it.
- **Explicit projection with low-storage Runge–Kutta** (Kim & Moin 1985 for
  the fractional step; Sanderse & Koren 2012 and Komen et al. 2021 for
  Runge–Kutta on collocated grids) needs no momentum linear solve and one
  pressure solve per stage. At the CFL below 1 of LES it is probably the
  cheapest option on a GPU, with implicit wall-normal diffusion where the
  viscous limit of thin wall cells would bind (to be measured, L39).

The prompt's "PISO + BDF2 / Crank–Nicolson" therefore becomes: CN with
converged outer correctors as the conservative default for LES, BDF2 allowed
and reported, explicit Runge–Kutta projection as an option to be measured.

### 3.5 Synthetic inlet turbulence

The original SEM (Jarrin et al. 2006) gives a velocity field that is not
divergence-free, which on a pressure-coupled incompressible solver produces
pressure noise at the inlet. The **divergence-free SEM** of Poletto, Craft &
Revell (2013) builds the fluctuation from vorticity-like eddies and is
divergence-free by construction; it takes the Reynolds-stress tensor and a
length scale per inlet face, which a RANS precursor (stage 6, the
`reynolds_stress` field of D129) gives. Digital filtering (Klein, Sadiki &
Janicka 2003) needs a structured inlet plane; on polyhedral inlet faces it
needs an auxiliary grid and interpolation, so it is not the first choice.

Reproducibility (D31, D93): eddy positions, signs and intensities come from
the counter-based streams of D122, keyed by (eddy id, regeneration count),
so they are the same on any partition; each inlet face sums the eddies
within its support (a gather). For channel and pipe validation a recycled
or periodic precursor is the simpler inflow and remains available (6p).

## 4. Subgrid-scale models

### 4.1 Properties

From the definitions, consistent with the comparison of Nicoud et al.
(2011); ν_sgs = (C Δ)² D(∇u) with ∇u the resolved velocity gradient:

| Model | ν_sgs at a wall | Pure shear | Solid-body rotation | 2D flow | Notes |
|---|---|---|---|---|---|
| Smagorinsky | O(1) (needs van Driest damping) | no | **zero** | no | reference only |
| WALE (Nicoud & Ducros 1999) | O(y³) | zero | not zero | not zero | the baseline everyone compares against; missing in the prompt |
| Vreman (2004) | **O(y)** | zero | **not zero** | not zero | prompt says y³ and "vanishes in rotation": wrong; its directional filter widths Δ_m need mesh axes, which polyhedra lack (use isotropic Δ) |
| σ-model (Nicoud et al. 2011) | O(y³) | zero | zero | zero | D_σ = σ₃(σ₁−σ₂)(σ₂−σ₃)/σ₁², σ_i singular values of ∇u |
| Coherent structure model (Kobayashi 2005) | O(y³) | zero | zero | not zero | C = C′ |F_CS|^{3/2} F_Ω, F_CS = Q/E, F_Ω = 1 − F_CS, **fixed** C′ (1/22 ◇); local, no test filter, not dynamic in the Germano sense; the prompt drops F_Ω |

◇ constants checked against the papers at implementation.

The derivation for Vreman: in solid-body rotation α = ∇u is antisymmetric
of rank 2, so β = Δ²αᵀα has a non-zero 2×2 principal minor and B_β ≠ 0.
At a wall only ∂u/∂y and ∂w/∂y are O(1), β is of rank one to leading order,
B_β = O(y²), and ν_sgs = O(y).

### 4.2 Implementation notes

- **Gradient.** The least-squares gradient of D24 is exact for linear fields
  on every cell type; node-based Green–Gauss adds nothing here.
- **Singular values for σ.** The eigenvalues of G = gᵀg in closed form
  (invariants and the trigonometric solution of the cubic), clamped against
  round-off when two eigenvalues meet; branch-free, so one kernel body for
  CPU and GPU (D8).
- **Filter width.** Cube root of the cell volume as default; the largest
  cell extent as option; the aspect-ratio factor of Scotti, Meneveau &
  Lilly (1993) as option for prism layers.
- **Slot.** All four models are algebraic: no transport equation, ν_sgs
  per cell from ∇u, Δ and the wall distance. That needs a model slot without
  equations in the turbulence interface (D175).
- **SGS scalar flux.** ν_sgs/Pr_sgs with a constant Pr_sgs through the Pr_t
  closure of D172. Sufrà & Steiner (2020) found Pr_sgs ≈ 0.5 in the core of
  heated pipe flow at Pr 10 and 20, rising to about 1 and above near the
  wall, and that WALE and the coherent structure model give realistic
  near-wall heat fluxes on a refined grid.
- **Dynamic procedure.** Germano et al. (1991) with Lilly's (1992)
  least squares needs a test filter, which on polyhedra is a face-neighbour
  average with a mesh-dependent filter ratio. Not needed for the four
  algebraic models; P17.

## 5. What the prompt leaves out

- **Wall-modelled LES (WMLES).** Wall-resolved LES of an industrial internal
  flow costs nearly as much as DNS near the walls; the industrial route is
  an equilibrium wall-stress model. The ODE of the thin boundary layer
  d/dy[(μ + μ_t) du/dy] = 0 with a damped mixing length is integrated from
  the wall to an exchange point, and the exchange point lies **above** the
  first cell (second to fourth cell, or about 0.1 δ) to avoid the log-layer
  mismatch (Kawai & Larsson 2012; reviews Larsson et al. 2016, Bose & Park
  2018). The same ODE with the energy equation and the Pr_t closure gives the
  thermal wall model. That ODE kernel is the one the high-Pr RANS wall
  function needs (section 6.3), so one kernel serves both (D173, L41).
- **Statistics.** D71 has running means and covariances. LES validation
  also needs averaging over homogeneous directions (planes, lines, the
  periodic directions), restart of the accumulators, and spectra and
  two-point correlations to check the box size (L43).
- **The DES route.** L5 lists DES next to LES. With PANS dropped, the hybrid
  route is SST-based (delayed DES, improved DDES), with the same per-face
  blending field (D174) and the shielding function already in the features
  of D166.
- **A cheap SGS check.** Decaying isotropic turbulence (Comte-Bellot &
  Corrsin 1971) on a 32³–64³ periodic box is the classical a-posteriori test
  of an SGS model and small enough for the hosted CI.

## 6. High-Prandtl heat transfer

### 6.1 The physics that sets the requirements

At high Pr the conductive sublayer is thinner than the viscous one
(δ_θ/δ_ν ∼ Pr^{−1/3}), and the smallest scalar scales are the Batchelor scales
η_B = η/√Pr (Batchelor 1959). The TU Graz DNS of heated pipe flow at
Pr_w = 56 needed a first-cell spacing of Δr⁺ = 0.05 (Irrenfried & Steiner
2024). Industrial cases with oils therefore rely on thermal wall functions in
RANS and on wall-modelled LES, not on resolved conductive sublayers. Two
effects decide the accuracy of those wall functions:

1. **The turbulent Prandtl number rises near the wall** well above its bulk
   value at high Pr: DNS of heated pipe flow at Re_τ = 360 gives
   Pr_t ≈ 1.0, 1.6 and 2.2 at the wall for Pr = 1, 10 and 20, falling to
   about 0.9 beyond y⁺ ≈ 10 (Irrenfried & Steiner 2016, figure 5; Steiner,
   Irrenfried & Brenn 2020; Irrenfried 2020).
2. **Temperature-dependent viscosity** changes the near-wall turbulence:
   at Re_τ,w = 380 and Pr_w = 56 the Nusselt number drops from 98.4 with
   constant properties to 75.9 with variable properties, while the skin
   friction changes little (Irrenfried & Steiner 2024).

### 6.2 Pluggable turbulent Prandtl number (D172)

The prompt is right that nothing may hardcode Pr_t. Today the specification
has `constant` and `kays-crawford` (D135); stage 7 uses μ_t/Pr_t in the
conductivity (D152) and in the buoyancy production (D155). D172 makes Pr_t a
field computed by a closure, read by all three:

- `constant` (0.85), `kays-crawford` (Kays 1994), as specified;
- **`kays-crawford-hp`**, the modified Kays–Crawford model of the TU Graz
  group (papers provided by the maintainer on 2026-11-17). With
  Pe_t = (ν_t/ν) Pr and a = √(2 (1/Pr_t,∞ − γ)):

      Pr_t = { γ + C Pe_t a − (C Pe_t)² [1 − exp(−a/(C Pe_t))] }⁻¹.

  Two published parameter sets exist, and the later one is the default:

  | Set | γ (wall value 1/γ) | C | Pr_t,∞ | Calibrated on | Source |
  |---|---|---|---|---|---|
  | `2019` (default) | μ_b* / (Pr_t,∞ + 0.166 Pr_w^0.7) | 1.5 | 0.85 | constant and variable properties, Re_τ 360–520, Pr_w 10–56, DNS and pipe experiments | Steiner & Irrenfried (2019), eqs. 9–11 |
  | `2016` | 1 / (Pr_t,∞ + Pr^0.88/10) | 0.8 | 0.9 | constant properties, Re_τ 360, Pr ≤ 20 | Irrenfried & Steiner (2016, eq. 23; 2017) |

  Both fix the one flaw of the original Kays–Crawford form (C = 0.2,
  Pr_t = 2 Pr_t,∞ at the wall for every Pr): the wall value rises with Pr
  as in the DNS. With constant properties (μ_b* = 1) the `2019` set gives
  1.02, 1.68, 2.20 and 3.42 for Pr = 1, 10, 20 and 50, the `2016` set 1.00,
  1.66 and 2.30 for Pr = 1, 10, 20; the DNS gives about 1.0, 1.6 and 2.2.
  For large Pe_t both tend to Pr_t,∞ (checked by expanding the exponential).
  The variable-property effect enters through Pe_t (with the local ν_t/ν)
  and through μ_b* = μ_b/μ_w, the ratio of the cross-section mean viscosity
  to the wall viscosity. A pipe bulk value does not exist in a general
  geometry, so Foehncraft defines it locally (our choice, to be checked by
  H-2): inside the wall function of D173 μ* is the viscosity at the end of
  the integration (first cell centre or exchange point) over the wall
  value, and Pr_w is the wall-face Pr; in the cell field away from walls the
  closure uses the local Pr and μ* = 1, where Pr_t is close to Pr_t,∞ in any
  case. It needs no wall distance, only ν_t, so it works in RANS and in the
  equilibrium wall model of L41. The set is chosen in the case file (C,
  Pr_t,∞ as overrides);
- `wall-profile`: the generic DNS-based form
  Pr_t = A + (Pr_t,∞ − A)/(1 + e^{−B y⁺})² of Steiner, Irrenfried & Brenn
  (2020). Their A and B were fitted per measurement case (A 2.2–2.6 and
  B 0.8–1.6 at Pr_w ≈ 22; A 4.4–4.5 and B 0.8–1.3 at Pr_w ≈ 50;
  Pr_t,∞ = 0.85), so they are **not** universal constants: the option takes
  A, B and Pr_t,∞ from the case file and has no defaults. With
  `kays-crawford-hp` available it is a secondary option;
- `sgs` with a constant Pr_sgs for LES (L40).

### 6.3 Thermal wall function with variable properties (D173)

The prompt's "revised P-function" has two published stages. The
**constant-property stage** is now in hand: Irrenfried & Steiner (2016, the
paper provided on 2026-11-17). Its content, which D173 adopts:

- The wall-layer relation dT⁺/dy⁺ = Pr_tot dw⁺/dy⁺ (momentum and heat fluxes
  constant near the wall) with
  Pr_tot = (1 + ν_t/ν) / (1/Pr + (ν_t/ν)/Pr_t) and
  dw⁺/dy⁺ = 1/(1 + ν_t/ν). Integrated to a w⁺ in the log region
  w⁺ = ln(y⁺)/κ + β this gives T⁺ = Pr_t,∞ (w⁺ + P), with the P-function
  P = ∫₀^∞ (Pr_tot/Pr_t,∞ − 1) dw⁺, the thermal resistance of the diffusive
  sublayer. D135's Jayatilleke form is this P with an empirical fit.
- Spalding's closed form of P rests on two simplifications that fail at high
  Pr: the leading term of van Driest's form near the wall,
  ν_t/ν ≈ (κ/A⁺)² w⁺⁴, and a constant Pr_t. (The paper prints eq. 8 as
  (κ/A⁺) w⁺⁴ and the factor of eq. 9 as (π/4) sin(π/4); both are misprints:
  the series of l⁺ = κ y⁺ (1 − e^{−y⁺/A⁺}) gives (κ/A⁺)², and integrating
  ∫₀^∞ dw⁺/(b + c w⁺⁴) gives the (A⁺/κ)^{1/2} and (π/4)/sin(π/4) of eq. 9
  only with c = (κ/A⁺)²; the product (π/4) sin(π/4) would give half the
  value of about 91 that figure 4 shows for Spalding at Pr = 20. We use the
  derived forms.) The paper replaces them by the full van Driest mixing
  length,
  ν_t/ν = l⁺² dw⁺/dy⁺ with l⁺ = κ y⁺ (1 − e^{−y⁺/A⁺}), A⁺ = 26, so that
  dw⁺/dy⁺ = 2 / (1 + √(1 + 4 l⁺²)), and by the Pr_t of `kays-crawford-hp`
  (set `2016`).
  (The true near-wall order of ν_t/ν is y⁺³ in DNS; the y⁺⁴ of van Driest
  does not matter for the integral, as the paper shows.) With κ = 0.34 and
  β = 4.5 at Re_τ = 360 (low-Re values, as in Wu & Moin 2008 and Durst et al.
  1995) the profile follows the DNS to the log layer at Pr = 1, 10, 20.
- The integral is cheap: a fixed number of quadrature points per wall face,
  or a table in (y⁺, Pr) built on the host (D31).

**Checked by us (2026-11-17), not only taken from the paper:** the integral
with these formulas, evaluated numerically for κ = 0.34, gives
P = 1.2, 59.2 and 105.5 for Pr = 1, 10, 20 with `kays-crawford-hp` (set
`2016`), and 99.7
at Pr = 20 with constant Pr_t = 0.9. The paper's figures show, at Pr = 20,
about 100 (constant Pr_t, figure 4), about 105 (eq. 7 with eq. 23, figure 8)
and about 106 (DNS). The integrated T⁺ at y⁺ = 100 is 17.5, 69.6 and 111.3
for Pr = 1, 10, 20, against about 19, 71 and 111 read from figure 9. The
reconstruction of eq. 23 is therefore right to the accuracy of reading a
figure (at Pr = 1 our T⁺ is about 1.5 lower), and
D173's `"integrated"` option gets unit tests with these numbers.

The **variable-property stage** is Steiner & Irrenfried (2019, provided
with four more papers and the thesis of Irrenfried (2020) on 2026-11-17).
It keeps the P-function and changes three things, with the wall values as
reference (subscript w, starred ratios ρ* = ρ/ρ_w, μ* = μ/μ_w, all at the
local temperature):

- **Momentum** (eq. 12, eq. 5):
  dw⁺/dy⁺ = 2 (τ⁺/μ*) / (1 + √(1 + 4 l_m⁺² τ⁺ ρ*/μ*²)) and
  μ_T/μ = (ρ*/μ*) l_m⁺² |dw⁺/dy⁺|, where τ⁺ is the total shear stress over
  τ_w.
- **Damping** (eqs. 6, 8, after Grifoll & Giralt 2000):
  l_m⁺ = κ y⁺ [1 − exp(−y⁺/(A_m⁺ μ_b*))] with
  A_m⁺ = A₀ [1 − exp(−y⁺/C_A)]^{1/2}, C_A = b A₀²/κ². This gives the correct
  near-wall order μ_T/μ ∼ y⁺³ that van Driest's constant A⁺ misses (y⁺⁴).
  Constants (their Table 1): κ = 0.40, A₀ = 33, b = 0.001 for Re_τ ≥ 750,
  the default; the low-Re values κ = 0.34, b = 0.0006 at Re_τ = 360 and
  κ = 0.36, b = 0.0008 at Re_τ = 500 only for the DNS checks.
- **Heat** (eqs. 1, 2): the enthalpy difference
  χ⁺ = (h_w − h) ρ_w u_τ / q_w = ∫ (1 + μ_T/μ)(q⁺/τ⁺) / (1/Pr + (μ_T/μ)/Pr_T) dw⁺
  with the local Pr = μ c_p/λ and Pr_T from `kays-crawford-hp` (set
  `2019`); q⁺ is the total heat flux over q_w. Enthalpy, not temperature,
  matches D92, so c_p(T) needs no extra term.
- **RANS form** (their section 6): the same integral in k-based units,
  y* = ρ_w C_μ^{1/4} k_P^{1/2} y/μ_w, blended with the viscous and
  conductive limits as w* = f₁ y* + (1 − f₁) w*_HRN and
  χ* = f₂ Pr_w y* + (1 − f₂) χ*_HRN, f₁ = exp(−3 μ_T/μ),
  f₂ = exp(−3 Pr_w μ_T/μ). This is the y* of D134, so the compound and
  scalable treatments keep their robustness where u_τ → 0.
- **Limits stated by the authors:** attached, near-equilibrium boundary
  layers; not separated or strongly non-equilibrium flow. Validated against
  their pipe experiments (Nu and c_f within 7 %, Pr_w ≈ 20 and 50,
  Re_b 3 500–27 000) and the water data of Hollingsworth et al. (1989).

D173 adopts this form with three generalisations for a general geometry,
all ours and all checked by H-2: q⁺/τ⁺ = 1 inside the wall layer (the pipe
expressions eqs. 3, 4 are used only in the 1D verification); μ_b* replaced
by μ* at the end of the integration (section 6.2); and the integration
upper limit at the first cell centre (RANS) or at the exchange point (L41).

The rest of D173:

- The wall-face interface of D134/D135 carries ρ, μ, k and c_p at the wall
  face (at T_w) and in the first cell, from the property database (D94–D103).
- **Semi-local scaling** (Huang, Coleman & Bradshaw 1995; Patel et al. 2015):
  y⁺_sl = ρ_P u_τ,sl y_P/μ_P with u_τ,sl = √(τ_w/ρ_P), next to the y⁺ of the
  wall values. (Name chosen because y* in D134 is already the k-based
  coordinate C_μ^{1/4} k^{1/2} y/ν.)
- **Integrated wall function** `"integrated"`: the wall-layer equations
  above integrated from the wall to the first cell centre on a fixed number
  of points, with μ(T), λ(T), c_p(T), ρ(T) at the local temperature, solved
  by a fixed number of fixed-point sweeps (no data-dependent branches,
  deterministic, one body for CPU and GPU). Options: the damping
  (`variable`, the 2019 default, or `van-driest` with A⁺ = 26, the 2016
  form) and the Pr_t closure of D172 (default `kays-crawford-hp`, set
  `2019`).
- Kader stays the default; Jayatilleke stays the option.

**Checked by us (2026-11-17) with the 2019 form and constant properties:**

| Pr | P (κ 0.40, b 0.001) | T⁺(y⁺ 100) | P (κ 0.34, b 0.0006) |
|---|---|---|---|
| 1 | 1.7 | 16.4 | 1.9 |
| 10 | 54.5 | 61.3 | 62.0 |
| 20 | 93.4 | 94.3 | 106.7 |
| 50 | 182.5 | 170.1 | 209.9 |

A full 1D pipe integration (eqs. 3 and 4 with the bulk velocity, mixed-mean
temperature) at the DNS conditions gives:

| Case | Nu, ours | Nu, DNS | Re_b ours / DNS | c_f ours / DNS |
|---|---|---|---|---|
| Re_τ 360, Pr 10 | 52.7 | 51.40 | 5 328 / 5 275 | 9.13 / 9.68 ·10⁻³ |
| Re_τ 360, Pr 20 | 67.4 | 65.95 | | |
| Re_τ 360, Pr 50 | 92.0 | 89.69 | | |
| Re_τ 500, Pr 10 | 76.7 | 74.89 | 7 782 / 7 680 | 8.26 / 8.70 ·10⁻³ |
| Re_τ 500, Pr 20 | 99.0 | 96.03 | | |

(Re_τ here is based on the diameter, as in the papers; DNS and thesis model
values from Irrenfried 2020, Table 6.3.) Our reconstruction is within 1.5 %
of the authors' own model numbers and within 3.1 % of the DNS Nusselt
numbers; the mixing-length momentum part is 5–6 % low in c_f at these low
Reynolds numbers. These numbers become unit tests of `"integrated"` (the
1D pipe driver is a test helper, not a solver feature). With constant
Pr_t the kernel must also reproduce the Kader profile to within the scatter
of Kader's own fit.

The same kernel, with the exchange point moved off the wall, is the
equilibrium wall model of L41.

### 6.4 SGS heat flux at high Pr

The gradient (Clark-type) model q_sgs ∝ Δ² ∇u·∇T on its own is not
dissipative enough and gives backscatter; it is known to be unstable without
clipping or a dissipative partner (Vreman, Geurts & Kuerten 1997 for the
momentum analogue). For high Pr, Irrenfried & Steiner (2017, a priori at
Pr 1, 10, 20) found the scale-similarity mixed dynamic model the most
accurate and tensor-diffusivity models no clear gain over cheaper ones;
Sufrà & Steiner (2020) found the numerical error of the scalar convection
larger, in rms, than the SGS terms themselves at high Pr. Irrenfried &
Steiner (2026, and the ETMM15 version) repeated the a-priori test with
variable properties (Re_τ,w 380 and 667, Pr_w 56 and 50): variable
properties change the SGS errors little; the dynamic two-parameter mixed
model (scale-similarity part plus dynamic eddy diffusivity) is better than
the dynamic Smagorinsky-type model, mainly in the normal flux components,
but neither captures the back-scatter, and the coefficients still deviate
from the a-priori values in the core. The gap between the Kolmogorov and
Batchelor scales, not the property variation, remains the main difficulty.
So the order is: first numerics (a bounded, low-dissipation scalar scheme)
and the eddy diffusivity with Pr_sgs, then mixed or gradient models as
research (P16), where the two-parameter mixed model is the first
candidate.

### 6.5 Baselines independent of the TU Graz data

- Kader (1981) and the Kays (1994) Pr_t for the profile (already D135).
- Gnielinski (1976) for Nu, with the property-ratio correction for liquids
  (Pr_b/Pr_w)^{0.11} (◇), and Petukhov (1970) with a viscosity ratio (◇).
- Public DNS: Pirozzoli (2023), turbulent pipe, Re_τ ≈ 1140,
  Pr 0.00625–16, open database; Schwertfirm & Manhart (2007), channel, Sc up
  to 49 (◇ data availability to be checked).

## 7. Verification suite and cost

### 7.1 Cases

| # | Case | Reference | What is checked | Tier |
|---|---|---|---|---|
| L-1 | Inviscid Taylor–Green, 32³, short time | exact: K constant | `kep` + implicit midpoint or converged CN: kinetic-energy drift at the iteration tolerance; the pressure part of the ledger falls as Δt² h² | CI |
| L-2 | Exact SGS properties | analytic fields | each model zero where section 4.1 says zero, ν_sgs ≥ 0, near-wall order from a Taylor-expanded wall field | CI (unit) |
| L-3 | Decaying isotropic turbulence, 32³–64³ | Comte-Bellot & Corrsin (1971) | energy decay and spectra per SGS model | CI short, campaign full |
| L-4 | Taylor–Green, Re 1600, 64³–256³ | pseudo-spectral 512³ (van Rees et al. 2011; high-order workshop case) | −dK/dt and 2ν × enstrophy against the reference to t = 10–20; their difference is the numerical dissipation, read from the ledger | campaign |
| L-5 | Channel, Re_τ 180 and 395, wall-resolved | Moser, Kim & Mansour (1999) | U⁺(y⁺), resolved Reynolds stresses, C_f per SGS model; mass-flow forcing of 6p | campaign |
| L-6 | Channel, Re_τ 1000–5200, wall-modelled | Lee & Moser (2015) | C_f and U⁺ with the exchange point above the first cell | campaign |
| L-7 | Channel inflow by DFSEM | the periodic channel L-5 | development length until C_f recovers | campaign |
| H-0 | 1D wall-layer kernel, constant properties | section 6.3 tables (our reconstruction, within 3.1 % of the DNS) | P, T⁺(100), pipe Nu, Re_b and c_f of `"integrated"` to quadrature accuracy | CI (unit) |
| H-1 | Heated pipe RANS, Pr 10, 20, 50 constant properties | TU Graz DNS, Irrenfried (2020) Table 6.3: Nu 51.40, 65.95, 89.69 (Re_τ 360), 74.89, 96.03 (Re_τ 500); c_f 9.68 and 8.70 ·10⁻³; Pirozzoli 2023 for Pr ≤ 16 | Nu, c_f, T⁺(y⁺) per Pr_t closure and wall function | CI (RANS, coarse) |
| H-2 | Heated pipe RANS, variable properties of the TU Graz heat-transfer oil, Pr_w 20–56 | DNS with variable properties, Irrenfried (2020) Table 5.9: Nu 58.11 (Re_τ,w 370, Pr_w 21), 68.66 (360, 50), 89.77 (520, 20); measurement series M20–M60 (Nu and c_f, Re_b 4 200–27 000); Irrenfried & Steiner (2024) for Pr_w 56; Gnielinski with property ratio | Nu (including the constant-to-variable property drop), c_f | CI (RANS, coarse) |
| H-3 | Heated channel/pipe WMLES, Pr 10–20 | TU Graz LES/DNS (O20) | Nu, T⁺, Pr_sgs | campaign |

The prompt's Pr = 50 wall-resolved LES is DNS-like in cost (section 6.1) and
not part of the plan; Pr 50 is covered by RANS (H-2) and by WMLES (H-3)
once the data are there.

Pass bands are set after the first campaign runs, as in stage 6 (D139):
a band is the reference uncertainty plus the GCI of the run (D67), never a
number chosen before the run.

### 7.2 Cost tier

LES channels, the Taylor–Green vortex at Re 1600 and the WMLES cases run for
hours to days and cannot run on the free hosted runners (D22). They form a
**campaign** tier (marker `campaign`, next to `slow` and `gpu`), run on the
maintainer's iMac or on Colab with a recorded commit, hardware and wall
time, and their verdicts are written into the verification manual. The hosted
CI runs only the CI rows above.

## 8. Items and their stages

Decided on 2026-11-17 with the maintainer's approval of this assessment.
Rows in [`requirements.md`](../requirements.md).

| ID | Item | Stage |
|---|---|---|
| D172 | Turbulent Prandtl/Schmidt number as a closure field (`constant`, `kays-crawford`, `kays-crawford-hp`, `wall-profile`, later `sgs`), read by D135, D152, D155 | where D135 is implemented (stage 6) for the first two; the others with stage 7, pull request 7e |
| D173 | Wall-face interface with variable properties, semi-local y⁺_sl, integrated thermal wall function `"integrated"` | stage 7, pull request 7e |
| D174 | Convection blending as a per-face field; modes `constant` (today's results bitwise) and `mesh-quality`; no Ducros sensor | stage 8 with D50 at the latest, earlier if a pull request touches the momentum face kernel |
| D175 | Algebraic (equation-free) model slot in the turbulence interface; filter width as a derived mesh field; singular values of ∇u in the invariants module | interfaces kept open in stage 6; first model with L40 |
| L39 | LES numerics: `kep` convection, CN with converged correctors and implicit midpoint, explicit Runge–Kutta projection to be measured, wiggle sensor, Rhie–Chow audit by the ledger | after v1.0, with L5 |
| L40 | SGS models WALE, σ, Vreman, coherent structure; SGS scalar flux with Pr_sgs | after v1.0, with L5 |
| L41 | Wall-modelled LES: equilibrium wall-stress and thermal wall model on the D173 kernel, exchange point above the first cell | after v1.0, with L5 |
| L42 | Divergence-free synthetic eddy inlet with reproducible eddies (D122) | after v1.0, with L5 |
| L43 | LES statistics (homogeneous averaging, restartable accumulators, spectra, two-point correlations) and the campaign test tier | after v1.0, with L5 |
| P16 | Gradient (Clark-type) and mixed scale-similarity SGS heat flux for high Pr, with clipping | research |
| P17 | Dynamic procedure with a polyhedral test filter | possible later |
| X2 | PANS in every form (on ζ-f or other models) | dropped by the maintainer 2026-11-17 (was P7) |
| O20 | TU Graz papers and DNS data for the high-Pr checks | answered 2026-11-17 (papers and thesis tables); raw DNS profiles only on request, section 9 |

## 9. Questions for the maintainer

1. **Papers (O20), answered 2026-11-17.** The maintainer provided
   Irrenfried & Steiner (2016, ETMM11), then Steiner & Irrenfried (2019) on
   variable properties, Steiner, Irrenfried & Brenn (2020), Irrenfried &
   Steiner (2026) with its ETMM15 version, and the thesis of Irrenfried
   (2020). The 2019 model is now the default of D172 and D173 (sections 6.2
   and 6.3); the 2016 constants stay as the set `2016`.
2. **Data (O20), mostly answered.** The thesis is published under CC BY
   4.0. Its tables can therefore ship with Foehncraft as reference data
   with attribution (D168): the DNS Nusselt numbers and friction factors
   (Tables 5.9 and 6.3), the measurement series M20–M60 (Tables 4.3–4.8) and
   the property fits of the heat-transfer oil of the rig (section 2.5),
   which also go into the material database (D94–D103) under a neutral name.
   Two points we found in the thesis: c_p is printed as
   0.818 + 3.664·10⁻³ T "J/(kg K)" but is in kJ/(kg K) (the figure shows
   1 700–2 700 J/(kg K)), and the DV520/20 Nusselt number is 89.77 in
   Table 5.9 and 90.68 in Table 6.4; we use Table 5.9, whose value the
   measurement (89.2 ± 4.4) confirms. Still not public: the DNS profiles
   (T⁺, Pr_t, fluxes over y⁺). The 2026 paper says "available on request".
   If the maintainer wants profile checks in H-1/H-2, a request to the
   authors is the way; until then the checks use the integral numbers above
   and the public Pirozzoli (2023) profiles up to Pr 16.
3. **Order after v1.0.** L5 (now L39–L43) stays after L27, L2, L12, L3, L4
   in the order of section 3 of the requirements. Say if LES should move
   up.

## 10. References

- G. K. Batchelor, *Small-scale variation of convected quantities like
  temperature in turbulent fluid. Part 1*, J. Fluid Mech. 5, 113–133 (1959)
- S. T. Bose, G. I. Park, *Wall-modeled large-eddy simulation for complex
  turbulent flows*, Annu. Rev. Fluid Mech. 50, 535–561 (2018)
- S. Choi, *Note on the use of momentum interpolation method for unsteady
  flows*, Numer. Heat Transfer A 36, 545–550 (1999)
- G. Comte-Bellot, S. Corrsin, *Simple Eulerian time correlation of full-
  and narrow-band velocity signals in grid-generated, "isotropic"
  turbulence*, J. Fluid Mech. 48, 273–337 (1971)
- F. Ducros, V. Ferrand, F. Nicoud, C. Weber, D. Darracq, C. Gacherieu,
  T. Poinsot, *Large-eddy simulation of the shock/turbulence interaction*,
  J. Comput. Phys. 152, 517–549 (1999)
- F. Durst, J. Jovanović, J. Sender, *LDA measurements in the near-wall
  region of a turbulent pipe flow*, J. Fluid Mech. 295, 305–335 (1995)
- F. N. Felten, T. S. Lund, *Kinetic energy conservation issues associated
  with the collocated mesh scheme for incompressible flow*, J. Comput. Phys.
  215, 465–484 (2006)
- M. Germano, U. Piomelli, P. Moin, W. H. Cabot, *A dynamic subgrid-scale
  eddy viscosity model*, Phys. Fluids A 3, 1760–1765 (1991)
- V. Gnielinski, *New equations for heat and mass transfer in turbulent
  pipe and channel flow*, Int. Chem. Eng. 16, 359–368 (1976)
- J. Grifoll, F. Giralt, *The near wall mixing length formulation
  revisited*, Int. J. Heat Mass Transfer 43, 3743–3746 (2000) (◇ volume and
  pages)
- F. Ham, G. Iaccarino, *Energy conservation in collocated discretization
  schemes on unstructured meshes*, CTR Annual Research Briefs, 3–14 (2004)
- D. K. Hollingsworth, W. M. Kays, R. J. Moffat, *The measurement and
  prediction of heat transfer in a turbulent boundary layer in water*, 7th
  Symposium on Turbulent Shear Flows, Stanford, Vol. 2 (1989)
- P. G. Huang, G. N. Coleman, P. Bradshaw, *Compressible turbulent channel
  flows: DNS results and modelling*, J. Fluid Mech. 305, 185–218 (1995)
- C. Irrenfried, *Convective turbulent near wall heat transfer at high
  Prandtl numbers*, doctoral thesis, Graz University of Technology,
  Monographic Series TU Graz, Computation in Engineering and Science 37
  (2020), DOI 10.3217/978-3-85125-724-3, licence CC BY 4.0
- C. Irrenfried, H. Steiner, *DNS of a turbulent heated pipe flow at high
  Prandtl numbers revisiting the P-function model*, 11th ERCOFTAC
  Symposium on Engineering Turbulence Modelling and Measurements (ETMM11)
  (2016) (◇ year and proceedings pages)
- C. Irrenfried, H. Steiner, *Subgrid-scale modeling of turbulent heat
  transport in forced convection at high molecular Prandtl numbers*, HEFAT
  2017, Portorož (2017)
- C. Irrenfried, H. Steiner, *DNS based analytical P-function model for
  RANS with heat transfer at high Prandtl numbers*, Int. J. Heat Fluid Flow
  66, 217–225 (2017) (journal version of the ETMM11 paper)
- C. Irrenfried, H. Steiner, *Modulation of turbulence flux budgets by
  varying fluid properties in heated high Prandtl number flow*, Flow
  Turbul. Combust. 112, 147–164 (2024)
- C. Irrenfried, H. Steiner, *Effect of variable fluid properties on the
  subgrid-scale modelling of turbulent heat transport*, Int. J. Heat Fluid
  Flow 122, 110706 (2026); conference version at ETMM15 (◇ year)
- N. Jarrin, S. Benhamadouche, D. Laurence, R. Prosser, *A synthetic-eddy
  method for generating inflow conditions for large-eddy simulations*,
  Int. J. Heat Fluid Flow 27, 585–593 (2006)
- C. Jayatilleke, *The influence of Prandtl number and surface roughness on
  the resistance of the laminar sub-layer to momentum and heat transfer*,
  Prog. Heat Mass Transfer 1, 193–321 (1969)
- B. A. Kader, *Temperature and concentration profiles in fully turbulent
  boundary layers*, Int. J. Heat Mass Transfer 24, 1541–1544 (1981)
- S. Kawai, J. Larsson, *Wall-modeling in large eddy simulation: length
  scales, grid resolution, and accuracy*, Phys. Fluids 24, 015105 (2012)
- W. M. Kays, *Turbulent Prandtl number – where are we?*, J. Heat Transfer
  116, 284–295 (1994)
- W. M. Kays, M. E. Crawford, *Convective Heat and Mass Transfer*, McGraw-Hill
  (1980)
- J. Kim, P. Moin, *Application of a fractional-step method to
  incompressible Navier–Stokes equations*, J. Comput. Phys. 59, 308–323
  (1985)
- M. Klein, A. Sadiki, J. Janicka, *A digital filter based generation of
  inflow data for spatially developing direct numerical or large eddy
  simulations*, J. Comput. Phys. 186, 652–665 (2003)
- H. Kobayashi, *The subgrid-scale models based on coherent structures for
  rotating homogeneous turbulence and turbulent channel flow*, Phys. Fluids
  17, 045104 (2005)
- E. M. J. Komen, J. A. Hopman, E. M. A. Frederix, F. X. Trias,
  R. W. C. P. Verstappen, *A symmetry-preserving second-order time-accurate
  PISO-based method*, Comput. Fluids 225, 104979 (2021) (◇ volume)
- J. Larsson, S. Kawai, J. Bodart, I. Bermejo-Moreno, *Large eddy
  simulation with modeled wall-stress: recent progress and future
  directions*, Mech. Eng. Rev. 3, 15-00418 (2016)
- M. Lee, R. D. Moser, *Direct numerical simulation of turbulent channel
  flow up to Re_τ ≈ 5200*, J. Fluid Mech. 774, 395–415 (2015)
- D. K. Lilly, *A proposed modification of the Germano subgrid-scale
  closure method*, Phys. Fluids A 4, 633–635 (1992)
- K. Mahesh, G. Constantinescu, P. Moin, *A numerical method for
  large-eddy simulation in complex geometries*, J. Comput. Phys. 197,
  215–240 (2004)
- I. Mary, P. Sagaut, *Large eddy simulation of flow around an airfoil near
  stall*, AIAA J. 40, 1139–1145 (2002)
- R. D. Moser, J. Kim, N. N. Mansour, *Direct numerical simulation of
  turbulent channel flow up to Re_τ = 590*, Phys. Fluids 11, 943–945 (1999)
- F. Nicoud, F. Ducros, *Subgrid-scale stress modelling based on the square
  of the velocity gradient tensor*, Flow Turbul. Combust. 62, 183–200 (1999)
- F. Nicoud, H. B. Toda, O. Cabrit, S. Bose, J. Lee, *Using singular values
  to build a subgrid-scale model for large eddy simulations*, Phys. Fluids
  23, 085106 (2011)
- A. Patel, J. W. R. Peeters, B. J. Boersma, R. Pecnik, *Semi-local scaling
  and turbulence modulation in variable property turbulent channel flows*,
  Phys. Fluids 27, 095101 (2015)
- B. S. Petukhov, *Heat transfer and friction in turbulent pipe flow with
  variable physical properties*, Adv. Heat Transfer 6, 503–564 (1970)
- S. Pirozzoli, *Prandtl number effects on passive scalars in turbulent
  pipe flow*, J. Fluid Mech. 965, A7 (2023); data at the authors' open
  database
- R. Poletto, T. Craft, A. Revell, *A new divergence free synthetic eddy
  method for the reproduction of inlet flow conditions for LES*, Flow
  Turbul. Combust. 91, 519–539 (2013)
- B. Sanderse, B. Koren, *Accuracy analysis of explicit Runge–Kutta methods
  applied to the incompressible Navier–Stokes equations*, J. Comput. Phys.
  231, 3041–3063 (2012)
- F. Schwertfirm, M. Manhart, *DNS of passive scalar transport in turbulent
  channel flow at high Schmidt numbers*, Int. J. Heat Fluid Flow 28,
  1204–1214 (2007)
- A. Scotti, C. Meneveau, D. K. Lilly, *Generalized Smagorinsky model for
  anisotropic grids*, Phys. Fluids A 5, 2306–2308 (1993)
- D. B. Spalding, *Monograph on turbulent boundary layers*, Imperial College
  Mech. Eng. Dept. TWF/TN/33 (1967)
- H. Steiner, C. Irrenfried, G. Brenn, *Near-wall determination of the
  turbulent Prandtl number based on experiments, numerical simulation and
  analytical models*, Heat Transfer Eng. 41 (15–16), 1341–1353 (2020)
  (journal version of the HEFAT 2017 paper)
- H. Steiner, C. Irrenfried, *Modelling of thermal wall boundary conditions
  with temperature-dependent material properties for use in RANS*, Int. J.
  Heat Fluid Flow 80, 108495 (2019)
- L. Sufrà, H. Steiner, *A priori assessment of subgrid-scale models and
  numerical error in forced convective flow at high Prandtl numbers*, Flow
  Turbul. Combust. 105, 377–392 (2020)
- A. Travin, M. Shur, M. Strelets, P. Spalart, *Physical and numerical
  upgrades in the detached-eddy simulation of complex turbulent flows*, in
  Advances in LES of Complex Flows, Fluid Mech. Appl. 65, 239–254 (2002;
  colloquium of 2000, cited as Travin et al. 2000 in D25)
- E. R. van Driest, *On turbulent flow near a wall*, J. Aeronaut. Sci. 23,
  1007–1011 (1956)
- W. M. van Rees, A. Leonard, D. I. Pullin, P. Koumoutsakos, *A comparison
  of vortex and pseudo-spectral methods for the simulation of periodic
  vortical flows at high Reynolds numbers*, J. Comput. Phys. 230,
  2794–2805 (2011)
- R. W. C. P. Verstappen, A. E. P. Veldman, *Symmetry-preserving
  discretization of turbulent flow*, J. Comput. Phys. 187, 343–368 (2003)
- A. W. Vreman, *An eddy-viscosity subgrid-scale model for turbulent shear
  flow: algebraic theory and applications*, Phys. Fluids 16, 3670–3681
  (2004)
- B. Vreman, B. Geurts, H. Kuerten, *Large-eddy simulation of the turbulent
  mixing layer*, J. Fluid Mech. 339, 357–390 (1997)
- X. Wu, P. Moin, *A direct numerical simulation study on the mean velocity
  characteristics in turbulent pipe flow*, J. Fluid Mech. 608, 81–112 (2008)
