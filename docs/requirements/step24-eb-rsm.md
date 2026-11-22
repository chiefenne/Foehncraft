# Step 24 – Elliptic-blending Reynolds-stress model (EB-RSM): assessment and later requirement

Requirements step 24 (2026-11-19). The maintainer supplied a specification
of an elliptic-blending Reynolds-stress model (Manceau 2015) with an
algebraic heat-flux companion (EB-AFM, Dehoux et al. 2012) and asked to
track it as a **future requirement**, after a critical review against the
literature and against Foehncraft. This page is that review and the
resulting requirement. No solver code is written before the maintainer's go
(project rule). Clean-room rule (D21): methods from the cited papers only;
other codes are not named.

Decided from this page (section 7): L47 (EB-RSM, after D171), P18 (EB-AFM),
O21–O23 (open points).

Contents

1. [Short answer](#1-short-answer)
2. [What Foehncraft has today](#2-what-foehncraft-has-today)
3. [The model](#3-the-model)
4. [Findings on the supplied specification](#4-findings-on-the-supplied-specification)
5. [Numerics and boundary conditions](#5-numerics-and-boundary-conditions)
6. [Verification and cost](#6-verification-and-cost)
7. [Items and their place](#7-items-and-their-place)
8. [References](#8-references)

## 1. Short answer

The core of the supplied model is sound: seven transport equations (six
stresses R_ij and ε) and the elliptic blending equation for α, the
wall-normal direction n = ∇α/|∇α| (no geometric normal needed), the α³
blending of the near-wall and homogeneous pressure–strain and dissipation,
generalised gradient diffusion (Daly & Harlow 1970) and the enhanced
C′_ε1. It reuses most of what stage 6 built for ζ-f and the elliptic-blending
k-ε (D132, D133): the α equation with L² outside the Laplacian on the
CG-AMG path, the bounded T and L, the exact wall distance (D48) and the
turbulence framework (D129).

Twelve points of the supplied text are wrong, incomplete or overclaimed
(section 4). The most important: the Coriolis term depends on the frame
formulation (absolute-velocity MRF needs half of it, the sliding mesh none),
the rotating-channel check needs a relative-frame option that D55 does not
have, and the tensor boundary conditions on symmetry, wedge and rotational
periodic faces are missing. With these corrections EB-RSM is a sensible
**optional** closure for secondary flows, swirl, system rotation and strong
curvature, planned after the explicit algebraic model of D171. The default
models stay SST and ζ-f.

## 2. What Foehncraft has today

| Area | State | Relevance |
|------|-------|-----------|
| Turbulence framework (D129) | generic bounded transport, field registry incl. `reynolds_stress`, `alpha`, `wall_normal`; one clip helper with smooth form (D59) | EB-RSM fits as one more model; stresses become solved, not estimated |
| ζ-f, ζ-f₀, EB k-ε (D132, D133) | α / f equations on CG-AMG, bounded T and L, `low-re`, `compound`, `compound-2023` | α equation, T, L and the low-Re wall values are reused |
| Explicit stress in momentum (D170 QCR, D171 EARSM later) | explicit non-linear stress path in the momentum assembly | the deferred anisotropic stress of section 5.1 uses the same path |
| MRF (D55) | **absolute velocity** solved, relative velocity in the convective flux | fixes the form of the rotation terms (finding 1) |
| Sliding mesh (D182, D185, D186) | inertial frame, rotor zone turns | no Coriolis term in the stress equations |
| Symmetry, wedge, rotational periodic faces (D29) | vector rotation of u; the missing transposed stress on wedge faces was a bug (fixed 2026-11-18) | the stress tensor needs Q R Qᵀ on these faces (finding 3) |
| Coupled solver, turbulence multigrid (D176–D181) | k-ε family on reduction multigrid; the EB k-ε step needs SIMPLEC on the finest mesh | EB-RSM is stiffer; convergence is a risk (section 6.2) |
| Stage 7 heat flux (D152, D172) | SGDH, GGDH, Pr_t closures | EB-AFM must not replace them (finding 11) |
| GPU parity, reproducibility (D31, D60, D93) | bitwise CPU/GPU comparison of each backend | applies to the new kernels unchanged |

## 3. The model

Taken from Manceau (2015). **Constants are not written here from memory**
(as for D133): they are taken from the paper at implementation and checked
by reproducing the paper's channel result first.

- k = ½ R_kk; anisotropy a_ij = R_ij/k − (2/3) δ_ij.
- T = max(k/ε, C_T (ν/ε)^½), L = C_L max(k^{3/2}/ε, C_η ν^{3/4}/ε^{1/4})
  (Kolmogorov lower bounds; not a realizability bound).
- α − L² ∇²α = 1, α = 0 on walls; n_i = ∂_i α / |∇α| (with a small
  regularisation where |∇α| → 0).
- DR_ij/Dt = P_ij + F_ij + G_ij + φ*_ij − ε_ij + D^ν_ij + D^t_ij, with
  P_ij = −(R_ik ∂_k U_j + R_jk ∂_k U_i), F_ij the frame term of section 3.1,
  G_ij the buoyant production (with stage 7 and the heat flux of P18 or
  D152), D^t_ij = ∂_k (C_μ/σ_k T R_kl ∂_l R_ij).
- φ*_ij − ε_ij = (1 − α³)(φ^w_ij − ε^w_ij) + α³(φ^h_ij − ε^h_ij), with
  ε^w_ij = (R_ij/k) ε, ε^h_ij = (2/3) ε δ_ij,
  φ^w_ij = −5 (ε/k) [R_ik n_j n_k + R_jk n_i n_k − ½ R_kl n_k n_l (n_i n_j + δ_ij)],
  φ^h_ij the **linear** SSG model (the quadratic SSG term is switched off in
  EB-RSM).
- Dε/Dt = (C′_ε1 P − C_ε2 ε)/T + ∂_k[(ν δ_kl + C_μ/σ_ε T R_kl) ∂_l ε],
  C′_ε1 = C_ε1 (1 + A₁ (1 − α³) P/ε), P = ½ P_kk.
- Wall values (low-Re): R_ij = 0, α = 0, ε_w = 2ν k/y² at the first cell
  (as in the EB k-ε of D133).

### 3.1 Rotation terms by frame formulation

For a field that is steady in a frame turning with Ω, the R_ij equation
gains a term that depends on how the mean velocity and the tensor
components are expressed. With C_ij = −2 Ω_k (R_im ε_jkm + R_jm ε_ikm):

| Formulation | Production from | Frame term F_ij | SSG rotation tensor |
|-------------|-----------------|-----------------|---------------------|
| Relative velocity, rotating basis | relative-velocity gradients | C_ij | W_ij(rel) + ε_mji Ω_m |
| **Absolute velocity, inertial basis (D55 MRF)** | absolute-velocity gradients | **½ C_ij** (the basis-rotation term of a co-rotating tensor) | W_ij(abs), no correction |
| Inertial frame, moving mesh (sliding, D182) | absolute-velocity gradients | **none** | W_ij(abs) |

The first row is the frame of `flow.frame_rotation` (PR #113: relative velocity, Coriolis force −2ρΩ×u in momentum, centrifugal part in the pressure, whole domain, refuses `[[mrf]]` and `[[sliding]]`); with it the stress equations take the full C_ij. The first two are the same equation (P_abs + ½C = P_rel + C), the check in
section 6.1. The supplied text adds the full C_ij in every case, which
double-counts rotation in our MRF and on sliding meshes.

## 4. Findings on the supplied specification

| # | Finding | Decision on this page |
|---|---------|-----------------------|
| F1 | Coriolis production C_ij (factor 2) is correct only in a relative-velocity frame; D55 needs ½C_ij, sliding meshes none; the same holds for the ε_mji Ω_m term in W_ij | Section 3.1 |
| F2 | The spanwise rotating channel cannot run with D55: Ω × r varies along the channel, so the absolute velocity is not streamwise periodic. The D171 EARSM checks name the same case | O22, answered: the whole-domain relative frame `flow.frame_rotation` (PR #113) runs the periodic rotating channel; the stress equations add C_ij there |
| F3 | No tensor boundary conditions: symmetry planes (R_nt = 0, reflected tensor), wedge/axis faces and rotational periodic faces (Q R Qᵀ) cannot be six scalar zero-gradients | Section 5.4, required |
| F4 | The stabilising viscosity C_μ k²/ε uses C_μ = 0.21, the GGDH coefficient; that is about 2.3 times the usual eddy viscosity. "Cancels to machine precision" holds only if the implicit and deferred νS terms use the identical face discretisation | Section 5.1 |
| F5 | Clipping shear components by Cauchy–Schwarz does not make the tensor positive definite; an active clip at convergence means the model is not solved | Section 5.3 |
| F6 | "Quadratic SSG" with g₂ = 0 is the linear SSG. The constant table uses the a_ij scaling (half the b_ij values of SSG); consistent, but the paper's own form must be used. C_ε2 = 1.83 matches the measured decay exponent, it is not "exact" | Section 3, constants from the paper |
| F7 | The cited "Wald, Lardeau, Manceau (2016), A wall-function-free approach for the EB-RSM" does not exist under that title; the paper is *Adaptive wall treatment for the elliptic blending Reynolds stress model* (Wald et al.). Billard & Laurence (2012) is the EB k-ε (D133), not an RSM wall treatment | References corrected; O21 |
| F8 | "Mixed y⁺ works seamlessly" is not established: at y⁺ ≈ 30 the length scale is L ≈ 1.8 y, so α and ∇α are not resolved in the first cell; the compound treatment (Popovac & Hanjalić 2007) was built for ζ-f | `low-re` first; adaptive treatment is O21 |
| F9 | The verification tests have no tolerances; DNS data licence not checked | Section 6, O23 |
| F10 | "Flagship general-purpose closure" conflicts with the defaults; cost (2–3 times per iteration, slower steady convergence) is not stated | Optional closure; section 6.2 |
| F11 | EB-AFM is offered for stage 7 in place of Pr_t, but stage 7 decided SGDH/GGDH (D152) and Pr_t closures (D172), and an AFM needs a solved R_ij; the supplied EB-AFM gives ranges instead of constants and a simplified near-wall and rotation form | Stage 7 unchanged; P18 |
| F12 | Overlap with D171 (EARSM): same target flows and cases | L47 after D171, shared cases |

## 5. Numerics and boundary conditions

### 5.1 Momentum coupling

The stress enters the momentum equation as an implicit stabilising
diffusion with ν_stab plus the deferred difference
−∇·[R − (2/3) k I + 2 ν_stab S] on the explicit-stress path of D170/D171;
the isotropic (2/3)ρk stays in the modified pressure p̂ (D129). ν_stab is
the model's own eddy-viscosity estimate (C_μ^{ν} R_nn T with the
wall-normal stress, as ζ-f, or 0.09 k²/ε as fallback; chosen at
implementation), not the GGDH C_μ. The implicit and the deferred νS terms
use the same face gradient, non-orthogonal correction and wall-face
treatment, so they cancel to round-off at convergence (checked: section
6.1). The same ν_stab is the `nu_t` field for the Pr_t closures (D172), the
wall shear of D134 and output; Rhie–Chow (D27) sees it through the
momentum diagonal as today.

### 5.2 Stress equations

Six segregated component equations with the same matrix structure; GGDH
split into an implicit isotropic part (non-negative, keeps the M-matrix
for the multigrid of D180) and a deferred anisotropic part; the
negative-definite parts of the pressure–strain and dissipation (−ε/k terms)
on the diagonal, the inter-component coupling deferred (Patankar). A 6 × 6
block-implicit cell coupling is a later option, not part of L47.

### 5.3 Realizability

R_ij must stay positive semi-definite. After each update, if a cell
violates it (negative normal stress, negative principal minor or
determinant), the anisotropy is scaled toward isotropy,
R = (2/3) k I + λ k a with the largest λ ∈ [0, 1] that is realizable, which
keeps k; k itself is bounded by the framework (D129). Through the one clip
helper of D129 (hard form in the primal, smooth form for linearisation and
adjoint), branch-free, deterministic on CPU and GPU. The number of limited
cells is counted and logged; validation runs require **zero limited cells
at convergence**.

### 5.4 Boundary conditions

- Walls: `low-re` (section 3) first; wall-function and adaptive treatments
  only after O21.
- Inlets: R_ij from a given tensor, or from k (or intensity) with the
  isotropic or a given anisotropy; ε as today; α = 1 or zero gradient.
- Symmetry planes: the reflected tensor (R_nt = 0 for the tangential
  directions t, zero normal gradient of R_nn and R_tt).
- Wedge/axis faces and rotational periodic faces (D29): the neighbour
  tensor rotated, Q R Qᵀ, which couples components across the face
  (deferred part), as for the transposed stress of momentum.
- MRF zone interfaces and sliding interfaces: R_ij expressed in the same
  (inertial) basis on both sides, so no transformation is needed in the
  D55 absolute formulation.

### 5.5 Interfaces and outputs

`turbulence.model = "eb-rsm"`; fields `reynolds_stress` (solved), `epsilon`,
`alpha`, `wall_normal`, `nu_t` (ν_stab); derived (D72): anisotropy
invariants (II, III, Lumley flatness A = 1 − (9/8)(A₂ − A₃)), the
barycentric map coordinates, the production tensor; Python bindings as for
the other models (D61); docs, gallery with field images, case and recipe
(D18).

## 6. Verification and cost

### 6.1 Checks

| # | Case | Reference | Target |
|---|------|-----------|--------|
| V1 | Frame equivalence | section 3.1 | the same steady field (to the solver tolerance) from the D55 absolute formulation and the relative formulation of `flow.frame_rotation` on an annular case |
| V2 | Stabilisation cancels | section 5.1 | converged fields independent of the scaling of ν_stab to the solver tolerance |
| V3 | Channel Re_τ 395 and 590, `low-re` | DNS of Moser, Kim & Mansour (1999) (the paper's own model results are not a target, D201) | C_f within 3 % of the DNS; U⁺ and the four stress profiles plotted; peak R_11⁺ and R_12⁺ within numbers set at implementation from the paper |
| V4 | Square duct, periodic | DNS of the duct at low Re (Gavrilakis 1992; Pinelli et al. 2010) | eight corner vortices; peak secondary velocity within a band set at implementation (RSMs commonly underpredict it) |
| V5 | Spanwise rotating channel | DNS of Kristoffersen & Andersson (1993) | asymmetric profiles, stabilised suction side; on `flow.frame_rotation`, case as `examples/cases/flow/channel_rotating_sst.toml` |
| V6 | U-bend, shared with D171 | as D171 | as D171 |
| V7 | Realizability | all of the above | zero limited cells at convergence |
| V8 | GPU parity | D60 | bitwise equal to CPU on the channel and the duct |

### 6.2 Cost and convergence

Seven transport equations and one elliptic equation against three or four
for ζ-f: about 2–3 times the turbulence work per iteration, and steady
convergence usually slower because the stresses feed momentum through the
deferred path. The elliptic-blending k-ε step already needs SIMPLEC on the
finest mesh (D133 finding). Solver-speed work for EB-RSM is not part of L47;
it is raised with the maintainer when the model exists (project rule).

## 7. Items and their place

| ID | Item | Place |
|----|------|-------|
| L47 | EB-RSM as an optional closure, sections 3–6 | After D171 (EARSM), after v1.0 |
| P18 | EB-AFM heat flux (Dehoux et al. 2012) with EB-RSM, constants and near-wall form from the paper; stage 7 (D152, D172) unchanged | After L47, next to P12 |
| O21 | Wall treatment of EB-RSM beyond `low-re` (adaptive treatment of Wald et al.; compound blending) | Assessment before its implementation |
| O22 | Whole-domain relative rotating frame for the rotating channel of L47 and D171 | Answered 2026-11-20: `flow.frame_rotation` (PR #113); L47 still adds the frame term C_ij of section 3.1 in the stress equations |
| O23 | Reference data licences for V3–V5 (channel and duct DNS) | Before the cases enter the repository |

## 8. References

- R. Manceau, *Recent progress in the development of the Elliptic Blending
  Reynolds-stress model*, Int. J. Heat Fluid Flow 51, 195–220 (2015)
- R. Manceau, K. Hanjalić, *Elliptic blending model: a new near-wall
  Reynolds-stress turbulence closure*, Phys. Fluids 14, 744–754 (2002)
- C. G. Speziale, S. Sarkar, T. B. Gatski, *Modelling the pressure–strain
  correlation of turbulence: an invariant dynamical systems approach*,
  J. Fluid Mech. 227, 245–272 (1991)
- B. J. Daly, F. H. Harlow, *Transport equations in turbulence*, Phys.
  Fluids 13, 2634–2649 (1970)
- J.-F. Wald et al., *Adaptive wall treatment for the elliptic blending
  Reynolds stress model* (HAL hal-01051800; details to be confirmed with
  O21)
- F. Dehoux, Y. Lecocq, S. Benhamadouche, R. Manceau, L.-E. Brizzi,
  *Algebraic modeling of the turbulent heat fluxes using the elliptic
  blending approach – application to forced and mixed convection regimes*,
  Flow Turbul. Combust. 88, 77–100 (2012)
- M. Popovac, K. Hanjalić, *Compound wall treatment for RANS computation of
  complex turbulent flows and heat transfer*, Flow Turbul. Combust. 78,
  177–202 (2007)
- K. Hanjalić, B. Launder, *Modelling Turbulence in Engineering and the
  Environment: Second-Moment Routes to Closure*, Cambridge University Press
  (2011)
- R. D. Moser, J. Kim, N. N. Mansour, *Direct numerical simulation of
  turbulent channel flow up to Re_τ = 590*, Phys. Fluids 11, 943–945 (1999)
- S. Gavrilakis, *Numerical simulation of low-Reynolds-number turbulent flow
  through a straight square duct*, J. Fluid Mech. 244, 101–129 (1992)
- A. Pinelli, M. Uhlmann, A. Sekimoto, G. Kawahara, *Reynolds number
  dependence of mean flow structure in square duct turbulence*, J. Fluid
  Mech. 644, 107–122 (2010)
- R. Kristoffersen, H. I. Andersson, *Direct simulations of low-Reynolds-number
  turbulent flow in a rotating channel*, J. Fluid Mech. 256, 163–197 (1993)
