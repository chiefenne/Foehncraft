# Stage 6 specification: RANS turbulence and wall treatments

Specification of stage 6 of the [stage plan](../requirements.md#7-stage-plan),
drafted on 2026-11-17 from the maintainer's stage-6 brief, before any code
(rule: requirements before code). It details the decisions **D129–D139** of
[`docs/requirements.md`](../requirements.md#2-decisions) and builds on the
earlier stage-6 decisions D12, D13, D28 (1), D29, D36, D48, D49, D51, D55,
D104, D108, D109 (1), D111, D116 and D124 (W4). Where this page and a D row
differ, the D row wins until this page is amended. *Amended 2026-11-17:*
the optional curvature and rotation correction and the quadratic
constitutive relation of SST (**D170**, section 3.6) and the explicit
algebraic stress model as a later item (**D171**, section 3.7).

Rules that hold for every item below:

- **Clean room** (D21): equations, constants, limiters and wall treatments
  come from the cited papers. No code of other solvers is read for
  implementation or copied. Every constant in the code carries the paper,
  and where possible the equation or table number, in a doc comment. Constants
  marked ◇ below are written from the literature as remembered at drafting
  time. They must be checked against the paper itself when the model is
  implemented (checklist item of the pull request), and the check is recorded
  in the doc comment.
- float64 on the reference path (D7). Gather-based assembly, no
  floating-point atomics, fixed summation order in global face order
  (D31, D93). The CPU and GPU run the same single-source kernel bodies
  (D8), bitwise equal.
- Adjoint-ready kernels (D59): pure functions of their inputs, generic over
  the `Real` type (dual numbers included), with smooth or frozen variants of
  every non-differentiable switch (section 2.6).
- The project says *reacting flows*, never the other word, and names no
  other CFD code.

Contents:

1. [Scope and reconciliation with earlier decisions](#1-scope-and-reconciliation-with-earlier-decisions)
2. [Common framework (D129)](#2-common-framework-d129)
3. [Turbulence models (D130–D133)](#3-turbulence-models-d130d133)
4. [Momentum wall treatment (D134)](#4-momentum-wall-treatment-d134)
5. [Thermal wall treatment and turbulent Prandtl number (D135)](#5-thermal-wall-treatment-and-turbulent-prandtl-number-d135)
6. [Numerics, GPU kernels and stability (D136)](#6-numerics-gpu-kernels-and-stability-d136)
7. [Boundary conditions, initial values and checks (D137)](#7-boundary-conditions-initial-values-and-checks-d137)
8. [Case-file keys (D138)](#8-case-file-keys-d138)
9. [Verification and validation suite, GCI pass bands (D139)](#9-verification-and-validation-suite-gci-pass-bands-d139)
10. [Implementation plan](#10-implementation-plan)
11. [Questions for the maintainer](#11-questions-for-the-maintainer)
12. [References](#12-references)


## 1. Scope and reconciliation with earlier decisions

### 1.1 What stage 6 delivers

| Item | Decision | Source of the requirement |
|---|---|---|
| Generic transport equation for scalars, passive scalars | D116, D51, D129 | earlier stage 6 scope |
| Exact wall distance from the wall-face hierarchy, PDE form as option | D48, D108 (1) | earlier stage 6 scope |
| k-ω SST (Menter, Kuntz & Langtry 2003) with the production limiter of the paper, Kato–Launder and Durbin options | D12, D130 | M1 decision 2026-11-14; brief 2026-11-17 adds the limiters |
| Standard k-ε (Launder & Spalding 1974) and realizable k-ε (Shih et al. 1995) | D12, D131 | M1 decision |
| k-ε-ζ-f (Hanjalić, Popovac & Hadžiabdić 2004) | D12, D132 | M1 decision (in place of RNG k-ε) |
| Elliptic-blending k-ε (Billard & Laurence 2012) | D133 (section 3.4) | brief 2026-11-17; maintainer decision 2026-11-17: joins ζ-f |
| Lag elliptic-blending k-ε (Lardeau & Billard 2016) | P11 (**proposed**) | brief, "optionally" |
| All-y⁺ momentum wall treatments: automatic (SST), compound (ζ-f, EB, optional for k-ε), scalable and standard wall functions, rough walls | D13, D134 | M1 decision; brief |
| Thermal law of the wall for any Pr (Kader 1981; Jayatilleke P-function), constant and Kays–Crawford turbulent Prandtl number | D135 | brief; needed by D111 |
| Turbulence numerics: bounded schemes, point-implicit sources, positivity, coupling order, bounded BDF2, AMG for scalars | D28 (1), D36, D136 | earlier scope; brief |
| Boundary conditions, inlet turbulence, checks and diagnostic codes | D29, D30, D137 | earlier scope; brief |
| Case-file keys | D138 | brief |
| Verification suite and GCI pass bands | D139 | brief; replaces the acceptance text of the stage-6 row where they differ |
| Moving reference frame, driven and rotational periodic boundaries with the affine map, study layer, statistics, derived fields, robustness controls | D55, D29, D109 (1), D62, D71, D72, D49 | earlier scope, **kept** in stage 6 unchanged (section 10) |
| Sliding meshes in both variants (conformal stepping and non-conformal by imprinting) next to the moving reference frame | D182, D56, D58 | **added** 2026-11-19 at the maintainer's request (section 10, 6g parts 5 to 7) |

### 1.2 Reconciliation

| Topic | Earlier decision | Brief of 2026-11-17 | This specification |
|---|---|---|---|
| Near-wall model for separation and heat transfer | ζ-f (Hanjalić 2004) with the compound wall treatment, in M1 (2026-11-14); PANS on ζ-f (P7) dropped by the maintainer 2026-11-17 (X2); literature line Hanjalić, newer Basara | Elliptic-blending k-ε (Billard & Laurence 2012), optionally the lag variant; compare with ζ-f | ζ-f stays as decided. EB k-ε **joins** ζ-f, after it in stage 6 (maintainer decision 2026-11-17, D133; section 3.5 compares the two). Nothing of ζ-f is dropped. |
| SST variant | SST 2003 | SST 2003 with Kato–Launder production limiter or Durbin realizability constraint | Kept. The 2003 limiter is the default; Kato–Launder and Durbin are options (D130). The ω production is written with the limited P̃_k, as corrected for the 2003 paper (section 3.1). |
| Realizable k-ε | M1 | "fast-converging industrial baseline" | Kept, unchanged (D131). |
| Standard k-ε | M1, for literature data | not mentioned | Kept (D12); its known stagnation-point anomaly becomes a test of the limiters (section 9). |
| Wall treatment | D13: standard and scalable WF for k-ε, automatic for SST, compound for ζ-f, thermal field with L1, roughness | All-y⁺ momentum and thermal treatments, CWT for ε models | Kept and detailed (D134, D135). The compound treatment is also offered for the k-ε family as an option. |
| Thermal wall functions and turbulent heat flux | "once L1 exists" (D13); energy equation is L1 (after M1) | stage 6: thermal law of the wall and Pr_t; stage 7: energy equation, GGDH, buoyancy | The thermal wall function and Pr_t are specified **and verified in stage 6 on a temperature carried as a constant-property passive scalar** (forced convection with constant properties is one-way coupled, so the passive scalar is exact for it). The energy equation, GGDH, buoyancy production and CHT belong to the stage-7 specification (D150 on), which takes the interfaces of D129 and D135 (agreed with the stage-7 thread). |
| Linear solver of the elliptic equations | D12: f of ζ-f "with the scalar path of stage 2, no AMG" | "leverage the AMG + CG solver for α" | Amended: f (ζ-f) and α (EB) are symmetric positive definite Helmholtz problems; they are solved with CG and the aggregation AMG of stage 3 by default, Jacobi-CG as option (D136). This uses D36's AMG for scalars. |
| Stage numbering | Stage 7 = parallel performance, internal flows, v1.0 | Stage 7 = energy, buoyancy, CHT | Not decided here. The stage-7 specification thread handles it. Stage 6 is the same stage in both plans. |
| Keys file | `docs/keys.md` is generated from the code (`foehn schema --markdown`, D77) | "TOML schema additions (`docs/keys.md`)" | The new keys are specified in section 8. `docs/keys.md` is regenerated in each implementation pull request, as the generator requires. |
| Crate names | `foehncraft-kernels` (no_std), `foehncraft-fv`, `foehncraft-cuda`, `foehncraft-flow` | same | Names match; new modules listed in section 10. |
| Wall-film and periodic constraints | D124 W4 (per-face exchange coefficients), D109 (1) (affine map on periodic pairs) | – | The wall treatment fills the shear, heat and mass slots of `derived::wall_exchange` (D124 W4). All turbulence variables cross periodic pairs through the affine map with a = 1, b = 0 (D109 (1)). |


## 2. Common framework (D129)

### 2.1 Notation

Incompressible, constant density ρ in stage 6 (variable density with
stage 7). ν = μ/ρ, ν_t = μ_t/ρ. Velocity gradient ∇u with components
∂u_i/∂x_j, from the least-squares gradient of D24.

- Strain rate S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i), rotation
  W_ij = ½(∂u_i/∂x_j − ∂u_j/∂x_i).
- Invariants: S = (2 S_ij S_ij)^½, Ω = (2 W_ij W_ij)^½, and
  |S| = (S_ij S_ij)^½ = S/√2 where a paper uses that form (EB,
  realizable k-ε). Each model section states which one it uses. *6d:*
  ζ-f takes |S| = S in the bounds of T and L (section 3.3).
- Production P_k = ν_t S² per unit mass (incompressible, Boussinesq),
  with the model's limiter where it has one.
- y: wall distance from D48 (exact closest-point distance by default).
  y_P: normal distance of the centroid of a near-wall cell from its wall
  face (section 4.1), used by the wall treatment.

### 2.2 Generic transport equation

Every turbulence variable is solved with the generic equation of D116,

  ∂(φ)/∂t + ∇·(F φ) − φ ∇·F − ∇·(Γ_φ ∇φ) = S_u + S_p φ,  S_u ≥ 0, S_p ≤ 0,

with F the face volume flux (mass flux over ρ), Γ_φ = ν + ν_t/σ_φ, and the
term −φ∇·F making the convection bounded while continuity is not yet
converged (section 6.1). The elliptic equations of ζ-f and EB use the
same interface with no time derivative and no convection:

  −∇·(L² ∇ψ) + ψ = S_ψ   (ψ = f or α).

### 2.3 Fields

Registered in the field registry (D115), written to the VTKHDF output on
request and to the restart file always:

| Field | Unit | Models | Notes |
|---|---|---|---|
| `k` | m²/s² | all | turbulent kinetic energy |
| `epsilon` | m²/s³ | k-ε family, ζ-f, EB; derived ε = β* k ω for SST | |
| `omega` | 1/s | SST | specific dissipation rate |
| `nu_t` | m²/s | all | kinematic eddy viscosity |
| `zeta` | – | ζ-f | v²/k |
| `f` | 1/s | ζ-f | elliptic relaxation function |
| `phi` | – | EB | v²/k (same physical quantity as ζ) |
| `alpha` | – | EB | elliptic blending parameter, 0 at walls, 1 far away |
| `wall_distance` | m | all | D48 |
| `turbulent_time_scale` | s | all | T of the model (k/ε, or bounded T of ζ-f and EB; 1/(β*ω) for SST) |
| `v2` | m²/s² | ζ-f (ζk), EB (φk) | wall-normal variance, for GGDH in stage 7 |
| `wall_normal` | – | EB | n = ∇α/\|∇α\| (zero where \|∇α\| vanishes), for GGDH in stage 7 |
| `reynolds_stress` | m²/s² | all | estimate of u'_i u'_j, symmetric, 6 components: Boussinesq (2/3)k δ_ij − 2ν_t S_ij for every model (stage 7 may refine it for ζ-f and EB) |

Wall fields (per wall face, boundary blocks of the output, D124 W6):
`y_plus`, `y_star`, `u_tau`, `wall_shear_stress`, `blending` (the Γ of the
compound treatment, or the automatic weight), and per scalar
`heat_transfer_coefficient` (h_w) and `t_plus`.

### 2.4 Coupling to momentum

- Momentum diffusion uses μ_eff = μ + ρν_t in the full stress
  ∇·[μ_eff(∇u + (∇u)ᵀ)]; the transposed part is a deferred explicit term
  (it vanishes for constant μ_eff and incompressible flow, not for a
  varying ν_t). It acts on symmetry faces too (symmetry planes, slip
  walls, the planes of a wedge): there u_n = 0, so its traction has only
  the normal part μ_eff |S| ∂u_n/∂n, and the normal stress on the face is
  2μ_eff ∂u_n/∂n, implicit like the Laplacian part. On a wedge it is half
  of the hoop stress −2μ_eff v/r² (6g: without it the cells next to the
  axis of the impinging jet were wrong).
- The isotropic part (2/3)ρk of the Reynolds stress is absorbed into a
  modified pressure p̂ = p + (2/3)ρk, which the pressure equation solves.
  Pressure boundary values given by the user are static pressures p and
  are converted with the boundary k. The output `pressure` is the static
  pressure p = p̂ − (2/3)ρk; `modified_pressure` is written on request.
  Forces and the pressure coefficient use p.
- The momentum wall shear at wall faces comes from the wall treatment
  (section 4) as a per-face effective viscosity, so the momentum kernel of
  stage 4 is unchanged in structure.

### 2.5 Hooks for later stages and correction fields

- **Model-correction fields (D104):** each turbulence equation takes the
  optional per-cell multiplier β on its production term and an additive
  source. Default off.
- **Buoyancy production (for stage 7):** the k, ε and ω equations each take
  a per-cell explicit part S_u ≥ 0 and implicit part S_p ≤ 0 named
  `buoyancy_production`, zero in stage 6. Stage 7 specifies G_b and its
  coefficients (D150 on).
- **Sustaining terms** (Spalart & Rumsey 2007) for the free-stream decay:
  optional ambient values (section 3.1 and 3.2).

### 2.6 Smooth and published switches

The models contain min, max, |·| and clipping. The primal solution uses the
**published hard form** by default, so that results can be compared with
the literature and the NASA Turbulence Modeling Resource (TMR) values. All
such operations go through one helper in `foehncraft-core::math`:

- `pmin(a, b)`, `pmax(a, b)`, `pabs(a)`: hard form (select instructions on
  the GPU, no branch);
- the same with a smooth variant, `smin(a, b; w) = ½(a + b − √((a − b)² + (w·m)²))`
  with m = |a||b|/(|a| + |b|) (half the harmonic mean; amended in PR #136,
  first m = ½(|a| + |b|)) and relative width w (default 10⁻³), and
  `smax` likewise, evaluated as min(a, b) − ½e²/(√(d² + e²) + d), d = |a − b|,
  e = w·m, free of cancellation. It is C¹ away from a = 0 and b = 0, and
  its error is bounded by w·m/2 ≤ w·min(|a|, |b|)/2, so a floor such as
  max(CD_kω, 10⁻¹⁰) stays a floor. (With the arithmetic mean that smooth
  floor sat at about w²|CD_kω|/16, and the SST flat plate with smooth
  limits never converged; PR #136);
- the linearisation and the adjoint (D59 (3)) always use the smooth
  variant or the frozen switch; `turbulence.smooth_limits = true` makes
  the primal use it as well. A test checks that the flat-plate C_f with
  smooth limits differs from the hard form by less than 0.1 % (PR #136:
  0.090 % on level 1 of `tmr_plate_sst.toml`, 0.1002 % on level 0).

### 2.7 Software structure

- The interface also admits models without a transport equation (ν_t per
  cell from ∇u, wall distance and a filter width; D175, step 21), so that
  the algebraic SGS models of L40 need no change of the framework.
- The model is an `enum TurbulenceModel` chosen on the host. Kernels are
  monomorphised per model (const generic or one function per model); no
  trait objects and no runtime branch on the model inside a kernel.
- A host-side trait `TurbulenceModule` in `foehncraft-flow` gives each
  model's fields, equations and their solve order, the wall-treatment kind
  and the log columns; it plugs into the module and step schedule of D118.
- Kernels (no_std, single source): `foehncraft-kernels::turb`
  (`sst`, `keps`, `zeta_f`, `eb`, `wall`, `thermal`, `invariants`). The
  device trait `FlowDevice` of stage 4 gets the matching operations; the
  CUDA device runs the same bodies (D8).


## 3. Turbulence models (D130–D133)

Equations are written per unit mass for constant ρ, with D/Dt = ∂/∂t + u·∇.

### 3.1 k-ω SST, 2003 version (D130)

Menter, Kuntz & Langtry (2003); constants of Menter (1994) where the 2003
paper keeps them; the TMR's SST-2003 page records the differences and the
correction of the ω production.

  Dk/Dt = P̃_k − β* k ω + ∇·[(ν + σ_k ν_t)∇k]

  Dω/Dt = (γ/ν_t) P̃_k − β ω² + ∇·[(ν + σ_ω ν_t)∇ω] + 2(1 − F₁) σ_ω2 (1/ω) ∇k·∇ω

- ν_t = a₁ k / max(a₁ ω, S F₂), a₁ = 0.31, with the strain-rate invariant
  S (2003 change from the vorticity of 1994).
- P_k = ν_t S², P̃_k = min(P_k, 10 β* k ω) (2003: factor 10, applied in
  both equations because the ω production uses P̃_k).
- F₁ = tanh(arg₁⁴), arg₁ = min[max(√k/(β* ω y), 500ν/(y² ω)), 4 σ_ω2 k/(CD_kω y²)],
  CD_kω = max(2 σ_ω2 (1/ω) ∇k·∇ω, 10⁻¹⁰) (2003 floor; in units of the
  per-unit-mass form, documented as such in the code).
- F₂ = tanh(arg₂²), arg₂ = max(2√k/(β* ω y), 500ν/(y² ω)).
- Blending of each constant: c = F₁ c₁ + (1 − F₁) c₂.
- Constants: β* = 0.09, κ = 0.41; set 1: γ₁ = 5/9, β₁ = 0.075, σ_k1 = 0.85,
  σ_ω1 = 0.5; set 2: γ₂ = 0.44, β₂ = 0.0828, σ_k2 = 1.0, σ_ω2 = 0.856.

**Production options** (`turbulence.production_limiter`, D138):

- `"menter"` (default): P̃_k as above.
- `"kato-launder"` (Kato & Launder 1993): P_k = ν_t S Ω instead of ν_t S²,
  then the 2003 clip. In a pure shear flow S = Ω and nothing changes; at a
  stagnation point Ω ≈ 0 and the spurious production vanishes. Ω comes
  from the absolute-frame velocity gradient, which is the solved variable
  in a rotating zone (D55); a note (TURB004) says so when both are active.
- `"durbin"` (Durbin 1996): the eddy viscosity is bounded by the
  realizability limit ν_t ≤ C_lim k/(√6 |S|), applied as
  ν_t = min(ν_t,SST, C_lim k/(√6 |S|)), C_lim = 1 (Durbin: "≤ 1").
  Derivation recorded in the code: the normal stresses
  (2/3)k − 2ν_t λ_i of the Boussinesq form stay non-negative when
  ν_t ≤ k/(3 λ_max), and λ_max ≤ (2/3)^½ |S| for a trace-free S_ij.
  Followed by the 2003 clip.

**Sustaining terms** (Spalart & Rumsey 2007; TMR "SST-sust"), optional
(`turbulence.ambient`): + β* ω_amb k_amb in the k equation and
+ β ω_amb² in the ω equation, so that free-stream values do not decay
between inlet and body. Off by default. The flat plate of section 9 uses
the TMR settings with them off and records the decay.

### 3.2 Standard and realizable k-ε (D131)

One implementation with two closures (D12).

  Dk/Dt = P_k − ε + ∇·[(ν + ν_t/σ_k)∇k]

**Standard** (Launder & Spalding 1974): ν_t = C_μ k²/ε,

  Dε/Dt = (ε/k)(C_ε1 P_k − C_ε2 ε) + ∇·[(ν + ν_t/σ_ε)∇ε],

C_μ = 0.09, C_ε1 = 1.44, C_ε2 = 1.92, σ_k = 1.0, σ_ε = 1.3. The
Kato–Launder option of 3.1 applies (P_k = ν_t S Ω), which is the form
Kato & Launder proposed for this model.

**Realizable** (Shih, Liou, Shabbir, Yang & Zhu 1995):

  Dε/Dt = C₁ S ε − C₂ ε²/(k + √(νε)) + ∇·[(ν + ν_t/σ_ε)∇ε],

C₁ = max(0.43, η/(η + 5)), η = S k/ε, C₂ = 1.9, σ_k = 1.0, σ_ε = 1.2;
ν_t = C_μ k²/ε with the variable

  C_μ = 1/(A₀ + A_s U* k/ε), U* = (S_ij S_ij + W_ij W_ij)^½,
  A_s = √6 cos φ_s, φ_s = ⅓ arccos(√6 W̃), W̃ = S_ij S_jk S_ki / |S|³,

W̃ clipped to [−1/√6, 1/√6] before the arccos (smooth clip in the
linearisation). A₀ = 4.0 ◇ (a value 4.04 is also in use; the paper's value
is taken). In a rotating zone the rotation term uses the rotation relative to
the frame, as in the paper; stated in the doc comment.

**Sustaining terms** for both, the analogue of 3.1 (Spalart & Rumsey
2007), optional with ambient values k_amb and ε_amb: the k equation gains
+ ε_amb and the ε equation + C_ε2 ε_amb²/k_amb (realizable:
+ C₂ ε_amb²/(k_amb + √(ν ε_amb))), so that the ambient state is a steady
solution of the decay.

**Notes of 6c (PR #71).**

- *Production default.* The standard k-ε takes the Durbin bound
  ν_t ≤ C_lim k/(√3 S) by default (`durbin`, C_lim = 1; the draft had
  `none`). It is inactive in an equilibrium log layer (C_μ S k/ε = 0.30
  against the bound 0.58), changes the channel of 9.4 by less than 0.2 %,
  and keeps the impulsive start of a flat plate bounded: with the compound
  treatment at y⁺ 30 the first iteration, k from the free stream and the
  shear of the no-slip start, gave P_k/ε ≈ 2·10⁴ in the wall cells, ν_t
  reached its cap and the run diverged in iteration 10. `none` stays
  available. The realizable model takes `none` only (its C_μ is the bound).
- *Free stream.* With the TMR inlet values (ε = C_μ k ω, ω = 125 U/L) the
  free stream decays: k by a factor of about 20 over the plate with the
  standard model; the realizable model takes k to zero at a finite
  distance (√(νε) > k in its destruction). The flat plate of 9.3 runs the
  k-ε family with the sustaining terms at the inlet values.
- *ε in the wall cells.* The fixed ε of the wall treatment (4.4) is also
  written into the field before the solve: a start value far from it (the
  free-stream ε 500 times the wall value) left the first solve negative in
  the wall row, ε at its floor and ν_t at its cap.
- *u* of the wall functions.* All k-ε treatments take u* = C_μ^¼ k_P^½
  (the realizable model with the constant C_μ = 0.09), y* = u* y_P/ν and
  ν_w = ν y*/u⁺(y*): no iteration.

### 3.3 k-ε-ζ-f (D132)

Hanjalić, Popovac & Hadžiabdić (2004); constants of the paper ◇.

  ν_t = C_μ ζ k T

  Dk/Dt = P_k − ε + ∇·[(ν + ν_t/σ_k)∇k]

  Dε/Dt = (C_ε1 P_k − C_ε2 ε)/T + ∇·[(ν + ν_t/σ_ε)∇ε]

  Dζ/Dt = f − (ζ/k) P_k + ∇·[(ν + ν_t/σ_ζ)∇ζ]

  L²∇²f − f = (1/T)(c₁ + C₂′ P_k/ε)(ζ − 2/3)

- T = max[min(k/ε, a/(√6 C_μ |S| ζ)), C_T (ν/ε)^½],
  L = C_L max[min(k^{3/2}/ε, k^{1/2}/(√6 C_μ |S| ζ)), C_η ν^{3/4}/ε^{1/4}].
- C_μ = 0.22, C_ε1 = 1.4(1 + 0.012/ζ), C_ε2 = 1.9, c₁ = 0.4, C₂′ = 0.65,
  σ_k = 1, σ_ε = 1.3, σ_ζ = 1.2, a = 0.6, C_T = 6, C_L = 0.36, C_η = 85 ◇.
- Wall values (integration to the wall): k_w = 0, ζ_w = 0,
  ε_w = 2ν k_P/y_P², f_w = −2ν ζ_P/y_P²; with the compound treatment these
  become blended values (section 4.4).
- Bounds: 0 ≤ ζ ≤ 2 (v² ≤ 2k), enforced by the source split of section
  6.2, not by clipping; a clip remains as the last guard and is counted.
- Curvature correction for the U-bend (P8) stays later. PANS on ζ-f (P7)
  was dropped by the maintainer on 2026-11-17 (X2, step 21).
- *6d (PR #73), as implemented:* the constants ◇ are confirmed against the
  page of the turbulence modelling benchmark working group (TMBWG). |S| in
  the bounds of T and L is the strain magnitude S = (2 S_ij S_ij)^½, as in
  the authors' implementation and on the TMBWG page (the draft took
  S/√2, which gave U_b⁺ +2.9 % in the channel against +1.9 %). C_ε1 P_k
  is evaluated as 1.4(P_k + 0.012 P_k/ζ) with P_k/ζ = C_μ k T S², finite
  as ζ → 0. ζ is floored at 10⁻¹⁰ and capped at 2 after its solve. A
  transport solve may lower k, ε or ζ of a cell by at most a factor 10
  per outer iteration (`ZETA_F_DROP`): without it an unconverged
  iterate drove ε to the absolute floor and the source of f blew up (the
  plate at y⁺ 0.5 diverged). Inlet and initial ζ default to 2/3 (keys
  `zeta`, in (0, 2]); f starts at 0; a two-equation restart continues
  with ζ = 2/3, f = 0. Findings of the flat plate (section 9.3): +5 %
  C_f, laminar from the TMR free stream, the k floor near the leading
  edge. *Since #141* the sink of k is that of 3.3.2, which removes the
  k floor: the plate converges, C_f +2.0 % (64 × 48) and +0.03 %
  (128 × 96) from White.

#### 3.3.1 k-ε-ζ-f₀ (Popovac 2023), part 6d2

Popovac (2023) solves f₀ = f + (N_ζ − 1)ζ/T in place of f, so that the
wall value f₀,w = (N_ζ − 2)·2νζ_P/y_P² is zero for N_ζ = 2: the stiff wall
value −2νζ/y² of f is gone. Constants of the paper's Table 1: C_μ = 0.09
with ν_t = (C_μ/ζ_eq) ζ k T, ζ_eq = 0.41 (so that C_μ/ζ_eq = 0.2195 takes
the place of 0.22), C_ε1 = 1.4(1 + 0.012/√ζ), N_ζ = 2; the others as in
3.3. Model `k-epsilon-zeta-f0`, default wall treatment `compound-2023`
(4.3).

  L²∇²f₀ − f₀ = (1/T)(c₁ − 1 + C₂′ P_k/ε)(ζ − 2/3) − (N_ζ − 1) ζ ε/k

  Dζ/Dt = f₀ − ζ P_k/(ε T) − (N_ζ − 1) ζ ε/k + ∇·[(ν + ν_t/σ_ζ)∇ζ]

*6d2, as implemented (maintainer instruction 2026-11-17: no code of the
author is available; take a default that may deviate from the paper but
fits the validation).* The deviations from the paper:

1. *Sign of the (N_ζ − 1) term in the f₀ equation.* The paper prints
   +(N_ζ − 1)ζ/T on the right of its eq. 7; the substitution of f₀ into
   the ζ-f equations gives −. With + the equilibrium ζ of the log layer
   falls to about 0.17 instead of ζ_eq = 0.41.
2. *ε/k in place of 1/T in the (N_ζ − 1) terms.* The two agree away from
   the wall; with ε/k the wall value f₀,w = 0 is consistent with the
   balance of ζ at the wall, where ε/k = 2ν/y²; the Kolmogorov-bounded 1/T
   is small there and leaves ζ without its wall sink.
3. *Sink of k: max(k/T_K, 2νk/d²)* (3.3.2; until #141 k/T_K alone),
   T_K = max(k/ε, C_T (ν/ε)^½), in place of −ε. With −ε the laminar
   leading edge of the flat plate holds k at its floor (the finding of
   6d) and the k residual stalls; with k/T_K every field of the plate
   converges to 10⁻⁷ without a floor. The realizability bound of T is
   left out of T_K: with it the plate relaminarizes and the channel gives
   U_b⁺ +30 %.
4. *Realizability bounds of T and L* take C_μ/ζ_eq = 0.2195, consistent
   with ν_t, not C_μ = 0.09.
5. *Realizability coefficient of T: a = 0.6 √2 = 0.849* (key `a`), in
   place of 0.6. With S = (2 S_ij S_ij)^½ this is the bound of Durbin
   (1996) written with (S_ij S_ij)^½. With 0.6 the bound is active in
   the separated shear layer of the backward-facing step and halves ν_t
   there: x_r/H came out at 7.4 to 7.7 against 6.26 (ζ-f: 7.8). With
   0.849 the bound no longer acts in the shear layer (a = 1.2 gives the
   same x_r) and x_r/H is 6.38; the channel moves by 0.1 %, C_f of the
   plate by +1.5 %.
6. *C_ε1 = 1.4 (1 + 0.009/ζ)* (#141, `c_epsilon1_zeta = 0.009`,
   `c_epsilon1_power = 1`), in place of the printed 1.4 (1 + 0.012/√ζ).
   With the sink of 3.3.2 the printed form gives U_b⁺ −11 % in the
   channel; 1.4 (1 + 0.038/√ζ) meets the channel but moves x_r/H of the
   step to 7.18, since √ζ is not small in the shear layer. The 1/ζ form
   acts only near the wall. `c_epsilon1_power = 0.5` restores the
   printed form.

The model as printed (deviations 1 and 4 only) gave U_b⁺ −11 % in the
channel; the cause is C_ε1 with √ζ, which the paper keeps (the 1/ζ form
of 2004 with 0.012 gave +8.5 % and did not converge with the sink
k/T_K). Results of the default since #141 (sections 9.3 to 9.6): channel
U_b⁺ 17.34 / 17.54 / 17.58 / 17.58 on 48 / 96 / 192 / 384 cells across
(first-cell y⁺ 1 to 0.1; +0.2 % from the DNS, observed order 2.2), −3.4
/ −1.2 % at first-cell y⁺ 10 / 30; flat plate C_f(0.97) +1.8 / +0.2 /
−0.2 % from White's correlation on 64 × 48 / 128 × 96 / 256 × 192 cells,
−0.3 / −2.7 % at y⁺ 10 / 30, all converged without a floor, and the
boundary layer turns turbulent from the TMR free stream on its own (ζ-f
stays laminar); backward-facing step x_r/H = 6.42 (+3 %); impinging jet
Nu₀ 196, mesh-converged, +38 % from the experiment (finding: the
excess of k at the stagnation point, which ζ-f₀ has no limiter for).
Until #141 (the sink k/T_K and the printed C_ε1) the channel drifted
from 17.79 to 18.39 over four meshes (+1.4 % on 48 cells). Keys:
`constants.n_zeta`, `constants.c_epsilon1_zeta` (C''ε1) and
`constants.c_epsilon1_power`; both ζ-f models take `c_epsilon1_zeta` and
`turbulence.k_sink`.

#### 3.3.2 The sink of k near walls (#141, maintainer decision 2026-11-20)

The level-3 channel of the grid-convergence study (#128) showed that
ζ-f₀ with the sink k/T_K does not converge with the mesh: U_b⁺ rose
from 17.79 to 18.39 over 48 to 384 cells across. A sink linear in k
with a rate that stays finite at the wall (1/T_K → 1/(C_T √(ν/ε)))
lets k fall like y instead of y², and the wall value ε_w = 2νk_P/y_P²
then grows like 1/y_P with every refinement. A sink that tends to ε at
the wall restores k ~ y². Both ζ-f models therefore take

  sink of k = k · max(1/T_K, 2ν/d²),  d the wall distance,

which is ε at the wall (k = ε y²/(2ν) there) and k/T_K away from walls,
where 2ν/d² is negligible. Unlike −ε or a blend with ε, it has no
feedback from ε: the forms max(k/T_K, min(ε, 2νk/d²)) and the blend
with w = min(1, 2νk/(d²ε)) met the channel and the plate but drove a
limit cycle (period about 80 iterations) in the near-laminar wall jet of
the impinging jet. The ambient (sustaining) source of k is the sink at
the ambient state, so the free stream holds its k for every form.

`turbulence.k_sink`: `"wall"` (the default), `"epsilon"` (−ε, the
papers), `"kolmogorov"` (k/T_K, the default of 6d2 to #140). EB takes
−ε and refuses the key. ζ-f (2004) gains as well: with −ε its plate
stopped at the k floor near the leading edge, with the new sink it
converges (C_f +2.0 / +0.03 % from White on 64 × 48 / 128 × 96).

### 3.4 Elliptic-blending k-ε (D133)

Billard & Laurence (2012), the "BL-v²/k" model. Variables k, ε, φ = v²/k and
the blending parameter α. *6e (PR #88):* the forms below are those of the
paper (Appendix A, Table A.10, eq. 18), checked term by term; they replace
the draft of this section, which was written from memory (ν instead of
ν/2 in the diffusion, max forms of T and L, ε_w = 2νk/y², E in the ε
equation, a C′_ε1):

  ν_t = C_μ φ k min(T, T_lim),  T = (k²/ε² + C_T² ν/ε)^½,
  T_lim = 0.6/(√6 C_μ φ (S_ij S_ij)^½),
  L = C_L (k³/ε² + C_η² ν^{3/2}/ε^{1/2})^½

  α − L²∇²α = 1,  α_w = 0  (L² outside the Laplacian; the conservative
  form −∇·(L²∇α) + α = 1 adds ∇L²·∇α and gives a channel C_f 7 % high at
  Re_τ 395, see the verification page)

  Dk/Dt = P_k − ε − E + ∇·[(ν/2 + ν_t/σ_k)∇k],
  E = C_ε3 (1 − α)³ (k/ε) 2ν ν_t (∂²U_i/∂x_j∂x_l)²

  Dε/Dt = (C_ε1 P_k − C*_ε2 ε)/T + ∇·[(ν/2 + ν_t/σ_ε)∇ε],
  C*_ε2 = C_ε2 + α³ (C_ε4 − C_ε2) tanh(|∂_j(ν_t/σ_k ∂_j k)/ε|^{3/2})

  Dφ/Dt = (1 − α³) f_w + α³ f_h − (φ/k) P_k + (2/k)(ν_t/σ_k) ∇φ·∇k + ∇·[(ν/2 + ν_t/σ_φ)∇φ]

with the near-wall part f_w = −(ε/2)(φ/k) and the homogeneous part
f_h = −(1/T)(C₁ − 1 + C₂ P_k/ε)(φ − 2/3). Constants of the paper: C_μ
0.22, C_ε1 1.44, C_ε2 1.83, C_ε3 2.3, C_ε4 0.4, σ_k 1, σ_ε 1.5, σ_φ 1,
C_T 4, C_L 0.164, C_η 75, C₁ 1.7, C₂ 0.9 (keys of `[turbulence.constants]`:
`c_mu`, `c_epsilon1` … `c_epsilon4`, `sigma_k`, `sigma_epsilon`,
`sigma_phi`, `c_t`, `c_l`, `c_eta`, `c_1`, `c_2`, `a` for the 0.6 of
T_lim). Acceptance requires that the channel at Re_τ = 395 stays within the
band of the DNS (section 9) before any other case. *6e:* C_f at 99.0 % of
the DNS.

- Wall values: k_w = 0, φ_w = 0, α_w = 0, ε_w = ν k_P/y_P² (the paper;
  half the value of ζ-f), or the compound value of section 4.4 with the
  same ν k/y² in its viscous part.
- Bounds: 0 ≤ α ≤ 1 (the discrete maximum principle; α is clipped to it as
  well), 0 ≤ φ ≤ 2.
- The wall-normal direction n = ∇α/|∇α| is exported (section 2.3); it is
  what GGDH and algebraic flux models of stage 7 use near walls. *6e:* not
  yet exported; it comes with its first user (GGDH, stage 7).
- Implementation (6e): the kernel of ζ-f with φ in the slot of ζ and α in
  the slot of f (section 3.5); the second derivatives of E and the
  transport of k in C*_ε2 are least-squares gradients of gradients.

**Lag EB k-ε (P11, proposed as possible later).** Lardeau & Billard (2016)
add a transport equation for the lag between stress and strain to the EB
model (a fifth equation, conference paper, AIAA 2016-1600). Proposal: keep
it as P11 and decide after the EB verification of stage 6, because the
gain is in non-equilibrium flows that are not in the stage-6 cases.

### 3.5 EB k-ε compared with ζ-f, and recommendation

| Aspect | ζ-f (Hanjalić et al. 2004) | EB k-ε (Billard & Laurence 2012) |
|---|---|---|
| Elliptic operator | −∇·(L²∇f) + f | α − L²∇²α (L² outside the Laplacian, as published; *6e*) |
| Matrix | symmetric, positive definite; M-matrix on orthogonal meshes, non-orthogonal part deferred (D26) | same |
| Condition number | grows like 1 + (L/h)²_max; L spans C_η η (near walls) to k^{3/2}/ε (outer flow), so the coefficient jumps by orders of magnitude: AMG needed | same operator, same scaling, same need for AMG |
| Right-hand side | dimensional (1/s), depends on P_k/ε, T and ζ − 2/3; changes sign | the constant 1 |
| Wall boundary value | f_w = −2νζ_P/y_P²: depends on the first cell (ζ_P and y_P²), couples back to ζ, and changes from face to face on polyhedral and prismatic wall layers with uneven y_P | α_w = 0: independent of the first cell |
| Bounds of the solution | none for f; f enters the ζ source with either sign | 0 ≤ α ≤ 1 by the discrete maximum principle (with an M-matrix and deferred correction limited, section 6.2) |
| Lagging | f must follow ζ every outer iteration | α depends only on L (k, ε); it can be solved every N outer iterations (`solver.elliptic.every`) |
| Near-wall source in the v²/k equation | through f (elliptic, global) | (1 − α³) f_w + α³ f_h: local, algebraic |
| Literature line | Hanjalić (PANS-ζ-f of Basara and co-workers dropped 2026-11-17, X2) | Billard & Laurence; elliptic blending also in later Reynolds-stress and algebraic heat-flux models |

The two matrices are equally conditioned; the difference is in the wall
value and the right-hand side. On arbitrary polyhedral meshes with uneven
first-cell distances the α equation is the better behaved one: its
boundary value is exact and its solution is bounded, whereas the f wall
value scales with 1/y_P² and links neighbouring wall faces of different
height through ζ. The compound treatment of ζ-f (section 4.4) removes most
of that sensitivity by blending f_w, which is why the 2004 and 2007 papers
pair them.

**Decision (maintainer, 2026-11-17): EB k-ε joins ζ-f, delivered after it in stage 6**
(pull request 6e, section 10). It reuses the k and ε equations, the v²/k
transport equation and the elliptic solver of ζ-f, so the extra code is one
source kernel and its tests. ζ-f stays the decided near-wall model (PANS on
ζ-f, P7, dropped 2026-11-17: X2). After the stage-6 validation (channel, backward-facing step,
impinging jet; section 9) the maintainer decides which of the two is the
recommended near-wall model for heat transfer in the user documentation.
Alternatives: EB replaces ζ-f (loses the maintainer's
literature line), or EB follows later as a P item (stage 6 stays smaller).


### 3.6 Optional corrections of SST: curvature and rotation, QCR (D170)

Two options of the SST model of 3.1, both off by default, so the published
SST 2003 stays the reference path. Both work on the velocity gradient that
the solver already has, so they add no transport equation.

**Curvature and rotation correction** (`turbulence.curvature_correction`):
the correction of Spalart & Shur (1997) applied to SST by Smirnov & Menter
(2009). The production P̃_k in the k and the ω equation is multiplied by

  f_r1 = max[0, min(f_rotation, 1.25)] ◇,

  f_rotation = (1 + c_r1) · 2 r*/(1 + r*) · [1 − c_r3 atan(c_r2 r̃)] − c_r1,

  r* = S/Ω,  r̃ = 2 Ω_ik S_jk [DS_ij/Dt + (ε_imn S_jn + ε_jmn S_in) Ω_m^rot] / (Ω D³),

  D² = max(S², 0.09 ω²),  c_r1 = 1.0, c_r2 = 2.0, c_r3 = 1.0 ◇,

with S and Ω the strain-rate and vorticity invariants of 3.1, Ω_ik the
rotation tensor in the frame of the zone and Ω^rot the angular velocity of
a rotating zone (D55; zero elsewhere). In a pure shear flow r* = 1 and
r̃ = 0, so f_r1 = 1 and the model is SST exactly; this is the first check.
DS_ij/Dt is the material derivative of the strain-rate tensor: its
convective part u·∇S_ij from the least-squares gradient of the six
components (one extra gradient of a symmetric tensor per outer iteration),
its time part from the BDF levels in a transient run. The upper limit 1.25
and the clip at 0 are those of the paper; the smooth limits of 2.6 replace
them in `smooth_limits` runs. *6h, PR #113:* Ω_ij is the rotation tensor of the
**absolute** velocity (with a velocity relative to a frame turning at Ω,
Ω_ij = W_ij − ε_ijk Ω_k): with that sign the unstable side of a rotating
channel gets f_r1 > 1 from both r* and r̃ (sources differ in the sign of
the frame term; the sign was fixed by this consistency and by the kernel
test that the inertial and the turning frame give the same f_r1). The
time part of DS_ij/Dt is BDF1 from the last time level (none in the first
step of a run or after a restart); the gradients of S_ij take zero normal
gradient on every boundary. Not with `[[sliding]]` yet (the strain rate
of the rotor cells would have to move with them). `flow.frame_rotation`
solves the whole domain in a turning frame for the rotating channel of
9.8 (the absolute velocity of `[[mrf]]` is not periodic there). A note (TURB code) says when the correction
is combined with Kato–Launder production, since both act on stagnation
and swirl regions.

**Quadratic constitutive relation** (`turbulence.qcr`): the Boussinesq
stress τ_ij = 2 μ_t S*_ij − (2/3) ρ k δ_ij (S* trace-free) is extended by
the QCR of Spalart (2000),

  τ_ij,QCR = τ_ij − C_cr1 (O_ik τ_jk + O_jk τ_ik),
  O_ik = 2 W_ik / √(∂u_m/∂x_n ∂u_m/∂x_n),  C_cr1 = 0.3 ◇,

with W_ik the rotation tensor, and, for `"qcr2013"`, the isotropic term of
Mani, Babcock, Winkler & Spalart (2013),

  − C_cr2 μ_t √(2 S*_mn S*_mn) δ_ij,  C_cr2 = 2.5 ◇.

For a model that carries k, the TMR page on QCR records how the 2013 term
relates to the (2/3) ρ k δ_ij term; which form is used is decided at
implementation against the paper and the TMR page (◇). The non-linear part
enters the momentum equation as the divergence of an explicit stress, a
deferred correction in the stored face order (D26, D31); the implicit
diffusion stays μ + μ_t. Its effect is the anisotropy of the normal
stresses that drives the secondary flow in non-circular ducts and corners,
which a linear model cannot produce. The production P_k of the turbulence
equations keeps the linear form (as the papers do).

*Implementation (PR #115):* the 2013 term is used in the paper's form
−C_cr2 μ_t √(2 S*:S*) δ_ij. The solver already carries (2/3) ρk in the
modified pressure, and an isotropic stress in incompressible flow only
shifts the pressure, so QCR2013 and QCR2000 give the same velocity. The
cell stress uses the cell velocity gradient and ν_t. The face value is
the interpolation of the two cells (the owner's on boundary faces), and
b_f of the momentum faces takes −τ_f·S after the face kernel. Both the
SIMPLEC and the coupled solver use that assembly.

### 3.7 Explicit algebraic Reynolds-stress model, later item (D171)

The explicit algebraic Reynolds-stress model (EARSM) of Wallin & Johansson
(2000) on the k-ω platform recalibrated for it by Hellsten (2005) is
planned after stage 6 (not in v0.6): the Reynolds-stress anisotropy is an
explicit tensor polynomial in the normalised strain and rotation tensors,

  a_ij = β₁ S_ij + β₃ (Ω_ik Ω_kj − II_Ω δ_ij/3) + β₄ (S_ik Ω_kj − Ω_ik S_kj)
         + β₆ (…) + β₉ (…)  (the full basis of the paper; in two dimensions
         β₁, β₃, β₄ remain) ◇,

with the β coefficients from the self-consistent solution of the paper
(a cubic for the production-to-dissipation ratio, solved in closed form, so
the model stays explicit and has no iteration). The k-ω equations are those
of Hellsten (2005) with its constants and cross diffusion (not SST's) ◇.
It reuses the explicit-stress path of the QCR (3.6) and the invariant
features of D166, and gives the stress anisotropy that the internal flows
of the validation (U-bend, square-duct secondary flow, swirl) need beyond
QCR. Its pull request, constants check and verification (rotating channel,
square duct, U-bend) are planned with the stage after v0.6; until then
D171 fixes only the model and the interfaces it needs, which stage 6 keeps
open: an explicit stress term in the momentum assembly (used first by the
QCR) and the strain and rotation tensors of the invariants module.


## 4. Momentum wall treatment (D134)

### 4.1 Wall-face interface

For every wall face f with near-wall cell P (the owner cell of the boundary
face):

- y_P = (x_P − x_f)·n_f, the normal distance of the cell centroid from the
  face plane, with n_f the outward unit normal; a floor of 10⁻³ of the
  cell thickness V_P/A_f normal to the face guards degenerate cells
  (counted in health); not 10⁻³ V_P^⅓, which on the thin wall cells of a
  unit-depth 2D mesh exceeds y_P itself (*6b*).
- U_t = |u_P − u_w − ((u_P − u_w)·n_f) n_f|, the tangential velocity
  relative to the wall velocity (D29 moving walls), and its direction t_f.
- Inputs: k_P, ν (face value), roughness k_s (section 4.5), the model's
  constants, and for scalars Pr (or Sc) and Pr_t.
- Outputs per face: u_τ, y⁺ = u_τ y_P/ν, y* = C_μ^{1/4} k_P^{1/2} y_P/ν,
  the wall shear stress τ_w (vector along −t_f), the effective wall
  viscosity μ_w = τ_w y_P / U_t used by the momentum kernel, the blending
  weight, and the first-cell values the model needs (ω_P, ε_P, the
  production P_k,P, f_w or α-related values), plus the per-face exchange
  coefficients of D124 W4 (shear; heat and mass from section 5).

The interface is model-independent: one wall-face kernel per wall
treatment, selected on the host. The same face list serves the thermal
treatment (section 5) and later the wall film (D123, D124).

**Near-wall cells with several wall faces** (corner cells, prism layers in
concave corners, polyhedra cut by two walls): the cell values that the
treatment imposes (ω_P, ε_P, P_k,P) are the area-weighted means over the
cell's wall faces, gathered in global face order from a per-cell list of
wall faces (CSR, built once on the host, D93 (1)), so the result is
bitwise reproducible.

### 4.2 Laws of the wall

- Log law u⁺ = (1/κ) ln(E y⁺), κ = 0.41, E = 9.793 (B = 5.2 ◇).
- Spalding (1961), one formula for all y⁺:
  y⁺ = u⁺ + e^{−κB}[e^{κu⁺} − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6],
  solved for u⁺ by Newton iteration with a fixed number of steps (no
  data-dependent loop on the GPU: 8 steps from the log-law guess, then a
  residual check counted in health).
- Viscous sublayer u⁺ = y⁺.

### 4.3 Treatments per model

| Treatment | Models | Wall shear | First-cell turbulence values | Claimed y⁺ range |
|---|---|---|---|---|
| `low-re` (integration to the wall) | SST, ζ-f, EB | τ_w = μ U_t/y_P | SST: k_w = 0 and the wall-face value ω_w = 10·6ν/(β₁ d₁²) with d₁ = 2 y_P, the wall-cell height (Menter 1994), as Dirichlet values (*6b*: replaces the first-cell value ω_P = 6ν/(β₁ y_P²) of the draft, see the note below); ζ-f, EB: ε_w, f_w or α_w of section 3 | y⁺ ≲ 1 |
| `automatic` (Menter, Kuntz & Langtry 2003) | SST (default) | u_τ = (u_τ,vis⁴ + u_τ,log⁴)^{1/4}, u_τ,vis = (ν U_t/y_P)^½, u_τ,log = κ U_t/ln(E y⁺) (fixed number of Newton steps, as for Spalding's law) | ω_P = (ω_vis² + ω_log²)^½, ω_vis = 6ν/(β₁ y_P²), ω_log = u_τ/(√β* κ y_P), imposed with the weight 1 − b, b = u_τ,vis⁴/u_τ⁴ (4.4) (*6g, PR #112:* ω_P = (ω_vis^1.2 + ω_log^1.2)^{1/1.2}, imposed with the weight tanh((y⁺/5)²), which also weighs the equilibrium production (u_τ⁴/ν)(1 − g) g of the first cell, g = du⁺/dy⁺ of Reichardt's law), on top of the wall-face ω_w of `low-re`; the log law takes y⁺ ≥ 11.06 (its lower end with κ = 0.41, E = 9.793), so that u_τ,log has no spurious root in the sublayer and b → 1 there (*6b*) | any; validated 0.5 ≤ y⁺ ≤ 100 |
| `compound` (Popovac & Hanjalić 2007) | ζ-f, EB (default), k-ε family (option) | τ_w from the blended profile: u⁺ = e^{−Γ} y⁺ + e^{−1/Γ} u⁺_log, Γ = 0.01 (y⁺)⁴/(1 + 5 y⁺) ◇ | ε_P = e^{−Γ} 2νk_P/y_P² + e^{−1/Γ} C_μ^{3/4} k_P^{3/2}/(κ y_P); P_k,P blended likewise between the integrated and the log-law value τ_w u_τ/(κ y_P); f_w (ζ-f) blended between −2νζ_P/y_P² and 0 | any; validated 0.5 ≤ y⁺ ≤ 100 |
| `compound-2023` (Popovac 2023) | ζ-f₀ (default), ζ-f (option) | u⁺ = [y*^{−n} + u⁺_log^{−n}]^{−1/n}, n = 10/3, u⁺_log = ln(8.9 y*)/(κψ) with u* = C_μ^{1/4} k_P^{1/2}; ψ = 1 − C_U (y_P − y_v/2)/(U_t κ u*), y_v = 11ν/u*, C_U = t·[(U·∇)U + ∇p/ρ], ψ in [0.5, 2] (the paper gives no bounds) | Γ = du⁺/dy*, from the blend; ε_P = Γ 2νk_P/y_P² + (1 − Γ) u*³/(κ y_P), P_k,P = u*² (τ_w/ρν) Γ, both fixed in the cell; k and ζ: zero flux at the wall; f₀,w = 0 | any; validated 0.3 ≤ y⁺ ≤ 30 (6d2: channel, plate, step) |
| `scalable` (Grotjans & Menter 1998) | k-ε family (default) | log law evaluated at ỹ* = max(y*, 11.06) | ε_P = C_μ^{3/4} k_P^{3/2}/(κ y_P) (with ỹ*), P_k,P = τ_w u_τ/(κ y_P) | y⁺ ≥ 1 without failure; accuracy claimed for y⁺ ≥ 30 (*6c*: the flat plate stays within 3 % of White's correlation from y⁺ 1 to 100) |
| `standard` (Launder & Spalding 1974) | k-ε family | log law with u_τ = C_μ^{1/4} k_P^{1/2} | as scalable, without the floor; below y* = 11.06 the viscous sublayer (τ_w = μU_t/y_P, P_k,P = 0) | 30 ≤ y⁺ ≤ 300 |

**Compound treatment with the k-ε family (*6c* finding, owner decision
2026-11-17).** The high-Re k-ε models have no near-wall damping, so
below the log layer the compound treatment overpredicts the wall shear
stress: flat plate +21 % at y⁺ 11, +40 % at 6, +58 % at 1.3; channel
at Re_τ 395 τ_w 2.2 times too high at y⁺ 0.4. The claimed range of
`compound` with the k-ε family is therefore 30 ≤ y⁺ ≤ 300 (TURB001
below), where it agrees with the wall functions; walls resolved to
y⁺ ≈ 1 take ζ-f (6d) or EB (6e). No two-layer model is added.

**Wall value of ω in `low-re` (*6b*, from the flat-plate runs).** Fixing
ω_P = 6ν/(β₁ y_P²) in the first cell, as the draft had it, made C_f(0.97)
of the plate converge at first order (p ≈ 0.9: −3.1 %, −1.2 %, −0.1 %
on 64 × 48 to 256 × 192 cells, against the TMR value of the day, since
dropped as a reference) towards about +1.1 % of it: the
1/y² profile of ω is not resolved by the linear flux between the first
two cells. The wall-face value 10·6ν/(β₁ d₁²) with d₁ the wall-cell
height, which pyaeroflow uses, gives −0.8 % already on
64 × 48 cells; pyaeroflow found the same (d₁ = the centroid distance, four
times larger ω_w, gives 2 % less). d₁ = 2 y_P is the wall-cell height of
hexahedra and prisms. *6g, PR #112:* the sweep of the factor (1, 3, 10,
30, 100) on the plate and the channel on three meshes keeps 10: it has the
smallest mesh dependence, and all factors meet on the finest mesh
(verification manual, stage 6).

The exponentials e^{−Γ} and e^{−1/Γ} are evaluated for every face (no
branch on the regime); Γ is floored at 10⁻³⁰⁰ before 1/Γ, and the in-house
exp (D93 (6)) returns exactly 0 below its underflow limit. For SST the
compound blending of Popovac & Hanjalić is not used: the automatic
treatment is the SST paper's own.

Popovac (2023) is cited for the later refinements of the compound treatment;
which of them enter is decided with the ζ-f pull request (6d) after reading
the paper (◇). *6d:* the maintainer decided (2026-11-17) to add the ζ-f₀
model of Popovac (2023) with its own compound treatment as the next
part, 6d2, for robustness in y⁺: f₀ = f + (N_ζ − 1)ζ/T is zero at the
wall, which removes the stiff wall value −2νζ/y² of f. 6d implements the
2007 treatment: k and ζ by zero gradient (dk/dy = dζ/dy = 0 at the wall),
ε and P_k of the first cell fixed by the blend with C_μ = 0.09 in the
log-law values, f_w = e^{−Γ}(−2νζ_P/y_P²).

### 4.4 Imposing first-cell values

ω_P (SST automatic) and ε_P (wall-function treatments) are fixed values in
the near-wall cell, with a weight m ∈ [0, 1] per cell: the row becomes
(1 − m)·(transport row) + m a_P (φ_P − value), so the off-diagonals are
scaled by 1 − m and b becomes (1 − m) b + m a_P·value (m = 1: the value is
fixed; the row stays diagonally dominant for every m). m is the area mean
of the faces' weights over the cell's wall faces (SST automatic: 1 − b of
4.3; 1 for the fixed values of the wall functions), applied without a
branch in the row-assembly kernel; the neighbour rows keep their
coupling (*6b*: the draft had the hard replacement only). For P_k,P the cell's production is replaced by the area-weighted
wall-face value. The compound treatment of ζ-f also blends the f wall
value. These replacements are part of the assembly, so the dual-number
Jacobian (D59) sees them.

### 4.5 Rough walls

Equivalent sand-grain roughness k_s per wall boundary (`roughness`, m),
k_s⁺ = u_τ k_s/ν:

- **Log-law shift** (all wall-function and compound treatments): u⁺ is
  reduced by ΔB(k_s⁺) of Cebeci & Bradshaw (1977): ΔB = 0 for k_s⁺ < 2.25;
  transitional ΔB = (1/κ) ln[(k_s⁺ − 2.25)/87.75 + C_s k_s⁺]·sin[0.4258 (ln k_s⁺ − 0.811)]
  for 2.25 ≤ k_s⁺ < 90; fully rough ΔB = (1/κ) ln(1 + C_s k_s⁺) for
  k_s⁺ ≥ 90; C_s = 0.5 (`roughness_constant`) ◇. The pieces join
  continuously; the smooth helper of section 2.6 joins them C¹ for the
  linearisation. *6c finding (rough pipe, 9.7):* with κ = 0.41,
  E = 9.793 the default C_s = 0.5 gives the fully rough intercept
  (1/κ) ln(E/C_s) = 7.3 instead of Nikuradse's 8.5, and the friction
  factor lies 8 to 11 % above Colebrook's; C_s = 0.4 gives +1 to +3 %.
  The default stays 0.5, the industry value that published k_s tables
  are used with; the recipe names 0.4 for Colebrook-type pipes.
- **k-ε family:** a rough wall face fixes ε like a smooth one; k and ε
  have zero gradient at every wall (k_w is not a boundary value of the
  wall functions).
- **SST** (`low-re` and `automatic`): rough-wall ω boundary value of
  Aupoix (2015), ω_w = (u_τ²/ν) ω_w⁺(k_s⁺) with the correlation of that
  paper ◇; with the automatic treatment ω_vis is replaced by this value.
- **ζ-f, EB**: the same ΔB in the compound profile.
- A warning (BC027) when k_s exceeds y_P on more than 10 % of a wall's
  area (the roughness elements then reach beyond the first cell centre).


## 5. Thermal wall treatment and turbulent Prandtl number (D135)

Stage 6 applies this to passive scalars (D51) and to a temperature carried
as a constant-property passive scalar (molecular Prandtl number
Pr = ν/a with a the thermal diffusivity, given as the scalar's
`diffusivity`). Stage 7 applies the same interface to the energy equation
(D92): the wall heat flux is then q_w = h_w (T_w − T_P) with the h_w below.

### 5.1 Turbulent diffusivity

Γ_φ = ν/Pr + ν_t/Pr_t (Sc and Sc_t for a mass-like scalar). Pr_t is a cell
and face field computed by the chosen closure, never a constant inside a
kernel (D172, step 21; stage 7 adds `wall-profile`). Options per
scalar (`turbulent_prandtl_model`):

- `"constant"` (default): Pr_t = 0.85 (`turbulent_prandtl`).
- `"kays-crawford"` (Kays 1994, after Kays & Crawford):
  Pr_t = 1 / {1/(2 Pr_t∞) + C Pe_t/√Pr_t∞ − (C Pe_t)² [1 − exp(−1/(C Pe_t √Pr_t∞))]},
  Pe_t = (ν_t/ν) Pr, Pr_t∞ = 0.85, C = 0.3. It raises Pr_t near walls and
  at low Pr (liquid metals).

GGDH and other anisotropic heat-flux closures are stage 7 (D150 on); they
use the fields of D129.

### 5.2 Thermal law of the wall

Wall-face interface as in 4.1, with Pr and Pr_t added. Outputs: T⁺ (or φ⁺),
the face conductance and the heat-transfer coefficient
h_w = ρ c_p u*/T⁺ (for a passive scalar: the exchange coefficient
ρ u*/φ⁺), with u* = C_μ^{1/4} k_P^{1/2} for the k-ε family and ζ-f/EB,
and the blended u_τ of 4.3 for SST. *6f (pull request 87):* with the
scalable momentum treatment T⁺ is taken at ỹ⁺ = max(y⁺, 11.06), the
virtual wall distance of its log law. *6g (pull request 90):* 6f took
u_τ = (τ_w/ρ)^{1/2} for all models; the k-ε family and ζ-f/EB now take
u* as above (T⁺ at y* = u* y_P/ν), so that the heat transfer does not
vanish with the wall shear at a stagnation point or a reattachment line
(the k-ε Nu₀ of the impinging jet was 37 with u_τ, 156 with u*).

- **`"kader"`** (default, all y⁺, any Pr; Kader 1981, wall form with the
  outer-layer factor dropped, y/δ → 0):
  T⁺ = Pr y⁺ e^{−Γ_t} + [2.12 ln(1 + y⁺) + β(Pr)] e^{−1/Γ_t},
  β(Pr) = (3.85 Pr^{1/3} − 1.3)² + 2.12 ln Pr,
  Γ_t = 0.01 (Pr y⁺)⁴/(1 + 5 Pr³ y⁺).
  Kader's paper covers temperature and concentration, so the same form
  serves a mass-like scalar with Sc.
- **`"jayatilleke"`** (log region only): T⁺ = Pr_t [u⁺ + P(Pr/Pr_t)],
  P = 9.24 [(Pr/Pr_t)^{3/4} − 1][1 + 0.28 e^{−0.007 Pr/Pr_t}]
  (Jayatilleke 1969). Blended with the conduction limit T⁺ = Pr y⁺ by the
  Γ_t of the Kader form, so that it is also defined for small y⁺.
- **Integration to the wall** (`"none"`): with a `low-re` momentum
  treatment the scalar flux is the molecular conduction ν/Pr (φ_P − φ_w)/y_P.
- With rough walls the thermal roughness correction is not part of stage 6;
  a note says that the smooth-wall T⁺ is used.

Boundary kinds of a scalar at a wall: fixed value, fixed flux (the flux
enters the first cell directly; the wall value is reported from T⁺), and
zero flux (default). Robin conditions and conjugate interfaces are stage 7.


## 6. Numerics, GPU kernels and stability (D136)

### 6.1 Discretisation

- **Convection:** linear upwind as a deferred correction on implicit
  upwind (D25) with the Venkatakrishnan limiter on by default for the
  turbulence variables, its threshold scaled with the variable:
  ε = K (φ_max − φ_min) over the domain (Wang 2000) with K = 0 by default
  (taken as 10⁻¹², the bounded Venkatakrishnan function). *6b:* the
  dimensional threshold (K h)³ with K = 5 of the draft is far above the
  differences of k (10⁻⁸ to 10⁻³ m²/s² on the TMR plate), switches the
  limiter off at the boundary-layer edge and left k floored in about 100
  cells with the residual stalled at 10⁻⁴; K = 10⁻³ still stalls at 10⁻⁶;
  first-order upwind as an option and for
  the first N iterations by the ramp of D49 (2). The bounded central
  scheme is refused for turbulence variables in steady runs (TURB006) and
  allowed in transient runs with β ≥ 0.1. The deferred correction is
  limited so that it cannot make the right-hand side negative (6.2).
- **Bounded form:** the term −φ_P ∇·F (with ∇·F the net volume outflow
  of the cell from the current fluxes) is added implicitly, so that the
  convection operator is an M-matrix whatever the continuity residual.
- **Diffusion:** over-relaxed non-orthogonal correction, deferred (D26);
  face diffusivity Γ_f by linear interpolation of ν_t (harmonic mean as
  option for strongly varying ν_t). The non-orthogonal correction of the
  elliptic equations is limited like the convective deferred correction.
- **Gradients:** least squares (D24) of k and ω for the SST cross
  diffusion; the same gradients serve the limiter.
- **Time:** BDF2 with the bounded fall-back of D28 (1), first order in the
  first step; Crank–Nicolson is refused for turbulence variables (not
  bounded), with a note that BDF2 is used for them.

### 6.2 Point-implicit linearisation and positivity

Each source is split into S_u ≥ 0 and S_p ≤ 0 with values φ* of the
previous outer iteration (Patankar 1980), so that the matrix keeps a
positive diagonal and the solution stays positive for positive boundary
values (M-matrix property). Every split uses select operations
(`max(x, 0)`, `min(x, 0)`), evaluated in all cells: no data-dependent
branch, so no warp divergence on the GPU.

| Equation | Term | S_u | S_p |
|---|---|---|---|
| k (all) | production P̃_k | P̃_k | 0 |
| k (SST) | −β* k ω | 0 | −β* ω* |
| k (ε models) | −ε | 0 | −ε*/k* |
| ω | (γ/ν_t) P̃_k | γ P̃_k/ν_t* (ν_t* floored, 6.2 floors) | 0 |
| ω | −β ω² | β ω*² (Newton) | −2β ω* |
| ω | cross diffusion CD = 2(1 − F₁) σ_ω2 ∇k·∇ω/ω | max(CD, 0) | min(CD, 0)/ω* |
| ε (standard, ζ-f) | (C_ε1 P_k − C_ε2 ε)/T | C_ε1 P_k/T | −C_ε2/T (*6c*: the default; `newton` gives C_ε2 ε*/T and −2C_ε2/T) |
| ε (realizable) | C₁ S ε − C₂ ε²/(k + √(νε)) | C₁ S ε* | −C₂ ε*/(k* + √(νε*)) |
| ζ | f − (ζ/k) P_k | max(f, 0) | min(f, 0)/ζ* − P_k/k* |
| φ (EB) | sources of 3.4 | positive parts | negative parts divided by φ* |
| f, α | elliptic | right-hand side (α: 1) | the unit coefficient (+1 on the diagonal) |
| all | deferred corrections D | max(D, 0) | min(D, 0)/φ* |
| k, ε, ω | buoyancy hook (stage 7) | its S_u | its S_p |

The divisions by φ* use φ* ≥ φ_floor. The Newton form of the ω
destruction makes the diagonal twice as strong as the Picard form and
converges faster; the Picard form (S_p = −β ω*) is an option. *6c:* for ε
the Picard form is the default (`solver.turbulence.destruction` unset):
the Newton form removes at most half of ε per iteration, and in a decaying
free stream k collapsed to its floor in 1 600 of 3 072 cells of the plate
while ε lagged behind; with Picard the plate converges in about 240
iterations.

**Floors** after each solve: k ≥ k_min, ε ≥ ε_min, ω ≥ ω_min with
k_min = 10⁻¹⁰ k_ref and ω_min, ε_min from k_min and the reference
viscosity ratio 10⁻², where k_ref is the largest inlet k (or the initial
k). With the splits above the floors should never act on a converging
run. Each activation is counted per equation in the health section
(D69) and gives TURB002 when more than 0.1 % of the cells needed it in an
iteration of the final tenth of the run (*6b*: in the first iterations the
inexact linear solves, relative tolerance 10⁻², can undershoot; this is
counted, not reported). ζ, φ are bounded to [0, 2] and α to [0, 1] the same way.

**Eddy-viscosity floor and bound:** ν_t ≥ 0 by construction; an optional
maximum viscosity ratio (`turbulence.max_viscosity_ratio`, default 10⁵)
caps ν_t/ν, counted in health.

### 6.3 Linear solvers

- k, ε, ω, ζ, φ: right-preconditioned BiCGStab with Jacobi (D36);
  AMG on the equation's own matrix as an option (`preconditioner = "amg"`)
  for stiff high-aspect-ratio layers (D36 *2026-11-15*).
- f, α: CG with the aggregation AMG of stage 3 (the matrices are
  symmetric positive definite, which the stage-3 AMG requires); Jacobi-CG
  as option. The hierarchy is rebuilt with `solver.amg_rebuild` like the
  pressure, and L² from the lagged k, ε.
- *6d:* the backends hold the hierarchies in slots (0 pressure, 1 f,
  2 k, 3 ε or ω, 4 ζ). The AMG of a transport equation is built from the
  symmetric part ½(A + Aᵀ) of its matrix (the stage-3 AMG needs a
  symmetric matrix; as a preconditioner of BiCGStab the symmetric part
  suffices). The keys: `solver.turbulence.preconditioner` (default
  `jacobi`) and `[solver.elliptic]` with `preconditioner` (default `amg`),
  `tolerance`, `max_iter`; the multigrid settings are those of
  `[solver.pressure.amg]`. Every hierarchy is rebuilt on the schedule of
  the pressure, and the matrix of its last setup goes into the restart
  file (`turb_amg_matrix_<slot>`), so a continued run builds the same
  hierarchies (D42).
- Tolerances per outer iteration: 10⁻² steady, 10⁻⁸ transient (as
  momentum), maximum iterations 50 (BiCGStab) and 100 (CG).

### 6.4 Coupling order

**SIMPLEC / SIMPLE outer iteration** (steps 1–7 of `simplec.rs` unchanged,
with μ_eff from the current ν_t and p̂ for the pressure):

8. velocity gradient of the corrected velocity; invariants S, Ω, |S|, W̃
   per cell (one invariants kernel);
9. wall treatment per wall face (section 4), then the per-cell gather of
   first-cell values;
10. the model's equations in this order, each assembled with the newest
    values of the ones before it:
    - SST: k, then ω (F₁, F₂ and CD_kω from the fields at the start of
      step 10);
    - standard and realizable k-ε: k, then ε;
    - ζ-f: k, ε, then f, then ζ;
    - EB: α (every `solver.elliptic.every` iterations, default 1), then
      k, ε, φ;
11. ν_t from the new fields, with optional relaxation α_νt (default 1);
    then the passive scalars (D51) with the new ν_t.

Under-relaxation defaults (steady): k, ε, ω, ζ, φ 0.8; f, α 1 (no time
derivative; they follow the lagged fields); transient 1.

**PISO** (plain): steps 8–11 once per time step after the pressure
correctors. **Iterative transient** (D28 (2), outer iterations per step):
once per outer iteration, so that a converged step is converged in all
variables. The residual of each turbulence equation (scaled as the
momentum residual) is part of the convergence test (`solver.tolerance` and
`outer_tolerance`).

**Wall distance** (D48) is computed once before the first iteration and
again only when the mesh moves (D56, L4).

### 6.5 GPU kernel design

- One **invariants** cell kernel (S, Ω, |S|, W̃ from ∇u), one **sources**
  cell kernel per model that writes S_u, S_p of all its equations and ν_t
  candidates, one **wall-face** kernel per treatment over the wall-face
  range, one **first-cell gather** kernel over the near-wall cells (CSR of
  wall faces per cell), and the generic face and row kernels of stage 2
  for convection and diffusion.
- No runtime branch on the model or regime inside a kernel: the model is a
  monomorphised instantiation; regimes are blended (compound, automatic)
  or selected with min, max and select.
- Transcendental functions (exp, ln, tanh, pow, cos, arccos, cube root as
  pow) come from the in-house routines of `foehncraft-core::math` on CPU
  and GPU alike (D93 (6)); no FMA contraction left to the compiler.
- All reductions (residual norms, y⁺ statistics, floor counts) are exact
  and order-independent (D31, D93 (4)); counts are integers.
- Memory per cell for SST: k, ω, ν_t, two old time levels of k and ω,
  S_u and S_p reused across equations, gradients of k and ω: about 20
  float64 values per cell on top of the flow solver; ζ-f and EB about 30.
  Measured and reported in the pull requests (D81 estimate updated).

### 6.6 Reproducibility, parity and adjoint

- Bitwise equal CPU and GPU results for each model on the flat plate and
  the channel (Colab T4 notebook, as for stage 5), and bitwise equal
  results on 1, 2 and 4 threads (CI).
- Dual-number Jacobian of every new kernel (sources, wall faces, Kader,
  Spalding Newton steps with fixed count) against central finite
  differences, as the stage-5 tests do (D59).
- Restart continues bitwise (D42), including the turbulence fields, the
  old time levels and the wall-treatment state.


## 7. Boundary conditions, initial values and checks (D137)

### 7.1 Inlets and backflow

Turbulence at velocity, mass-flow, static- and total-pressure inlets and for
the backflow of pressure outlets (keys already reserved in stage 5), one of:

- `turbulence_intensity` I and `turbulence_length_scale` l:
  k = (3/2)(I |u|)², ε = C_μ^{3/4} k^{3/2}/l, ω = k^{1/2}/(C_μ^{1/4} l);
- `turbulence_intensity` and `viscosity_ratio` r = ν_t/ν:
  ε = C_μ k²/(r ν), ω = k/(r ν);
- `k` with `epsilon` or `omega` directly (numbers or expressions, D53),
  as the TMR cases specify them;
- defaults when nothing is given: I = 0.05 and r = 10, with a note
  (TURB009) naming the boundary.

`zeta` (v²/k, ζ-f and EB, default 2/3) at inlets (D12). With a table
profile (D29) k and ε or ω can be columns of the table, so that a profile
from a periodic pipe or channel run (D111) feeds a jet or a duct.
*6g (pull request 90):* the header names the columns (`u_n` or `u_x`,
`u_y`, `u_z`; `k` with `epsilon` or `omega`; `zeta`), a boundary that
gives turbulence values as well is refused (BC029), an unknown column
too (BC019); `python3 -m foehncraft.profiles` makes the table from a
line of the periodic run.

### 7.2 Walls, symmetry, periodic, wedge

- Walls: section 4 (k_w = 0 in the integrated limit; zero normal gradient of
  k with wall functions), `roughness`, moving walls through U_t.
- Symmetry, slip wall, wedge: zero normal gradient of all turbulence
  scalars; f, α likewise (α is only 0 on no-slip walls).
- Periodic and rotational periodic: every turbulence variable through the
  paired-face map with a = 1, b = 0 (D109 (1)).
- Convective outlet: the turbulence scalars take the same convective
  condition as the velocity.

### 7.3 Initial values

`initial.turbulence_intensity` and `initial.viscosity_ratio`, or
`initial.k` with `initial.epsilon`/`initial.omega` (numbers or expressions);
by default the values of the first inlet in the order of the case file,
else I = 0.01 and r = 1 with the initial speed (or 1 m/s when that is 0).
ζ, φ start at 2/3, α at 1 with the elliptic equation solved before the
first iteration, f at 0.

### 7.4 Checks and diagnostic codes

The boundary-condition check (D30) and the validation of D39 gain these
codes (registry of D78, documented in `docs/interfaces.md` when
implemented):

| Code | Severity | Title |
|---|---|---|
| `BC027` | warning | A wall roughness k_s larger than the first-cell distance on more than 10 % of the wall's area. |
| `BC028` | warning | A turbulence value that the run does not use: `zeta` with SST or k-ε, `omega` with an ε model, or a turbulence value on a boundary other than an inlet or a pressure outlet (*6b*: the backflow took k and ω of the cell; *PR #140*: a pressure outlet uses its values on the faces where the flow enters, fixed values there and zero gradient where it leaves, decided per face in every outer iteration); BC009 stays for laminar cases. |
| `BC029` | error | Two incompatible turbulence specifications on one inlet (for example `turbulence_length_scale` and `viscosity_ratio`, or `k` with intensity). |
| `TURB001` | warning | First-cell y⁺ outside the claimed range of the wall treatment (section 4.3) on more than 10 % of a wall's area (reported per boundary with the area-weighted min, mean and max). *6b:* above 2 with `low-re`; between 3 and 20 with `automatic`, where the blending overpredicts the wall shear stress (section 9.3; *6g, PR #112:* the new weights of the automatic treatment hold C_f within 2.6 % for y⁺ 0.24 to 85, and TURB001 warns above y⁺ 300 only). *6c:* k-ε family: below 30 with `scalable`, outside 30 … 300 with `standard` and `compound`. |
| `TURB002` | warning | A positivity floor or bound acted on more than 0.1 % of the cells in an iteration (equation and cell count). *6b:* reported only for iterations in the final tenth of the solves; the run record counts all of them. |
| `TURB003` | warning | Inlet turbulence outside the usual range: I > 0.5, or a viscosity ratio below 10⁻² or above 10⁵. |
| `TURB004` | note | Kato–Launder production in a rotating zone: Ω is the absolute-frame vorticity. |
| `TURB005` | note | No wall boundary: the wall distance is infinite (F₁ = 0, F₂ = 0; the ζ-f and EB near-wall terms vanish). |
| `TURB006` | error | A convection scheme or time scheme that is not bounded for a turbulence variable (central in a steady run, Crank–Nicolson). |
| `TURB007` | warning | Model constants changed in `[turbulence.constants]`. |
| `TURB008` | warning | The ν_t cap `max_viscosity_ratio` acted (cell count). |
| `TURB009` | note | Inlet turbulence defaulted to I = 0.05, r = 10 (boundary named). |
| `TURB010` | note | The curvature correction of SST together with Kato–Launder production: both act on stagnation and swirl regions (the paper uses the standard production). | *6h, PR #113.*

The guidance table of D30 gets the recommended wall treatment per model
and y⁺ band and the inlet settings of the canonical cases.


## 8. Case-file keys (D138)

Additions to the case file (format version bumped with the first
implementation pull request, D126). `docs/keys.md` is generated from the
code; the table below is the target.

| Key | Type | Default | Unit | Run time | Description |
|---|---|---|---|---|---|
| `turbulence.model` | "laminar" \| "k-omega-sst" \| "k-epsilon" \| "realizable-k-epsilon" \| "k-epsilon-zeta-f" \| "elliptic-blending-k-epsilon" | `"laminar"` | | | Turbulence model (D12, D130–D133). |
| `turbulence.wall_treatment` | "automatic" \| "compound" \| "low-re" \| "scalable" \| "standard" | per model: SST `automatic`, ζ-f and EB `compound`, k-ε family `scalable` | | | Momentum wall treatment (D134). |
| `turbulence.production_limiter` | "menter" \| "kato-launder" \| "durbin" \| "none" | SST `menter`, standard k-ε `durbin` (*6c*), realizable `none` (the only value) | | | Production form of SST and standard k-ε (3.1, 3.2); for the k-ε family `menter` is the clip min(P_k, 10ε). |
| `turbulence.durbin_coefficient` | number | `1.0` | | | C_lim of the Durbin bound. |
| `turbulence.ambient` | { k, omega \| epsilon } | none | SI | | Sustaining terms (Spalart & Rumsey 2007). |
| `turbulence.curvature_correction` | boolean | `false` | | | SST: the curvature and rotation correction of Smirnov & Menter (2009) (3.6, D170); constants `c_r1`, `c_r2`, `c_r3`. |
| `turbulence.qcr` | "none" \| "qcr2000" \| "qcr2013" | `"none"` | | | Quadratic constitutive relation (3.6, D170); constants `c_cr1`, `c_cr2`. |
| `turbulence.max_viscosity_ratio` | number | `1e5` | | yes | Cap of ν_t/ν. |
| `turbulence.smooth_limits` | boolean | `false` | | | Smooth min and max in the primal solution (2.6). |
| `turbulence.wall_distance` | "exact" \| "pde" | `"exact"` | | | D48. |
| `turbulence.constants.<name>` | number | the paper's | | | Constant override by its literature name (`beta_star`, `a1`, `c_mu`, …); TURB007. |
| `solver.turbulence.scheme` | "upwind" \| "linear-upwind" \| "central" | `"linear-upwind"` | | | Convection scheme of the turbulence equations. |
| `solver.turbulence.limiter` | "none" \| "venkatakrishnan" \| "venkatakrishnan-wang" | `"venkatakrishnan"` | | | Limiter. |
| `solver.turbulence.limiter_k` | number | `0.0` | | | Relative threshold K of the limiter, ε = K (φ_max − φ_min) (6.1; *6b*: 0 is taken as 10⁻¹², the bounded Venkatakrishnan function). |
| `solver.turbulence.relaxation` | number | 0.8 steady, 1 transient | | yes | Under-relaxation of k, ε, ω, ζ, φ. |
| `solver.turbulence.tolerance` | number | 1e-2 steady, 1e-8 transient | | yes | Relative tolerance per solve. |
| `solver.turbulence.max_iter` | integer | `50` | | yes | BiCGStab iterations per solve. |
| `solver.turbulence.preconditioner` | "jacobi" \| "amg" | `"jacobi"` | | | Preconditioner (D36). |
| `solver.turbulence.destruction` | "newton" \| "picard" | ω `"newton"`, ε `"picard"` (*6c*) | | | Linearisation of −βω² and similar terms (6.2). |
| `solver.turbulence.viscosity_relaxation` | number | `1.0` | | yes | Relaxation α_νt of ν_t. |
| `solver.elliptic.tolerance` | number | 1e-2 steady, 1e-8 transient | | yes | f and α solves. |
| `solver.elliptic.max_iter` | integer | `100` | | yes | CG iterations. |
| `solver.elliptic.preconditioner` | "amg" \| "jacobi" | `"amg"` | | | Preconditioner of CG. |
| `solver.elliptic.every` | integer | `1` | | yes | Solve α (EB) every N outer iterations. |
| `boundary.<name> (type = "wall").roughness` | number | `0.0` | m | | Equivalent sand-grain roughness k_s. |
| `boundary.<name> (type = "wall").roughness_constant` | number | `0.5` | | | C_s of Cebeci & Bradshaw. |
| `boundary.<name> (inlets).k`, `.epsilon`, `.omega` | number \| string | none | SI | | Direct turbulence values (7.1). |
| `boundary.<name> (inlets).zeta` | number | `0.6667` | | | v²/k for ζ-f and EB. |
| `initial.turbulence_intensity`, `initial.viscosity_ratio`, `initial.k`, `initial.epsilon`, `initial.omega` | number \| string | from the first inlet | SI | | Start values (7.3). |
| `scalars.<name>.diffusivity` | number \| string | required | m²/s | | Molecular diffusivity (ν/Pr for a temperature). |
| `scalars.<name>.turbulent_prandtl` | number | `0.85` | | | Pr_t or Sc_t (D135). |
| `scalars.<name>.turbulent_prandtl_model` | "constant" \| "kays-crawford" | `"constant"` | | | D135. |
| `scalars.<name>.wall_function` | "kader" \| "jayatilleke" \| "none" | `"kader"` (`"none"` with `low-re`) | | | Thermal or mass law of the wall. |
| `scalars.<name>.initial` | number \| string | `0.0` | | | Start value. |
| `boundary.<name>.scalar.<scalar>.value`, `.flux` | number \| string | zero flux at walls, 0 at inlets | | | Per-physics sub-table of D116. |
| `output.fields` | adds `k`, `epsilon`, `omega`, `nu_t`, `zeta`, `f`, `phi`, `alpha`, `wall_distance`, `turbulent_time_scale`, `v2`, `wall_normal`, `reynolds_stress`, `modified_pressure`, the scalars | | | | Section 2.3. |
| `output.wall` | adds `y_plus`, `y_star`, `u_tau`, `blending`, `heat_transfer_coefficient`, `t_plus` | | | | Wall fields. |

Log (D125): residual columns `k`, `epsilon`, `omega`, `zeta`, `f`, `phi`,
`alpha` and one per scalar; the solution-settings matrix gets the
turbulence and elliptic rows (D127); the extremes block gets ν_t/ν max with
its cell and the wall y⁺ min, mean and max per wall boundary.


## 9. Verification and validation suite, GCI pass bands (D139)

### 9.1 How a verdict is formed

- **Discretisation uncertainty** from D67 (Celik et al. 2008) on the three
  finest of at least four systematically refined levels (refinement ratio
  r ≥ 1.3 per level, the same in each direction): observed order p,
  extrapolated value φ_ext and U_num = GCI_fine = 1.25 |e₂₁|/(r^p − 1)
  (relative). Oscillatory or divergent convergence is reported by D67
  and counts as a failed verification.
- **Verification verdict:** p inside the stated band and GCI_fine below
  the stated limit.
- **Validation verdict:** |φ_fine − φ_ref| ≤ B_model + U_num + U_ref, with
  B_model the model band in the tables below, U_ref the uncertainty of the
  reference (experimental uncertainty, or the spread of the DNS or TMR
  values), all relative. Model bands are set before the runs and are not
  changed to make a model pass; a model that misses its band is reported
  as a finding with the published result of the same model next to it.
- Every case runs on the hexahedral mesh of the built-in mesher and on
  its polyhedral variant (`mesh.generate.type = "poly"`) where the case
  allows; both must pass. Every case ships with gallery entry, case
  files and recipe (D18), and its numbers enter the verification manual
  (D68).

### 9.2 Code verification (exact solutions)

| Case | Quantity | Reference | Pass |
|---|---|---|---|
| Manufactured solutions for each model (construction after Eça, Hoekstra, Hay & Pelletier 2007): smooth k, ε or ω, ζ/φ, f/α and ν_t fields, all source terms of section 3 active, Dirichlet boundaries | L2 and L∞ errors of every variable | exact | p of the L2 error in [1.8, 2.2] on hexahedra, [1.6, 2.4] on polyhedra and tetrahedra; the same with the dual-number Jacobian test green |
| Homogeneous decay in a periodic box (no walls, uniform start) | k(t), ε(t) or ω(t) | exact ODE solutions: k-ε k = k₀(1 + (C_ε2 − 1)ε₀t/k₀)^{−1/(C_ε2 − 1)}; k-ω (SST with F₁ fixed at 1 or 0) ω = ω₀/(1 + βω₀t), k = k₀(1 + βω₀t)^{−β*/β} | time order of BDF2 in [1.9, 2.1]; error at the finest Δt below 10⁻⁶ relative; bounded BDF2 never leaves positivity |
| Sustaining terms | k, ω constant in time at the ambient values | exact | to round-off |
| Wall distance (D48) | d on planes, cylinders, spheres, and the PDE form | exact | exact (BVH) to round-off; PDE form within 1 % at the first cell layer |
| Wall functions as functions | u⁺(y⁺) of Spalding and the blended profile, T⁺ of Kader for Pr 0.01, 0.71, 7, 100, P-function | published formulas evaluated in high precision | relative error below 10⁻¹² (Newton steps converge) |
| Laminar heated channel, passive scalar: in 6a open (developed profile at the inlet, zero-gradient outlet), periodic with the heat-flux jump of D109 (1) once 6p is in | Nu | 140/17 = 8.235 | p in [1.8, 2.2], Nu within 0.5 % (D111) |

**Manufactured solutions (row 1, PR #148; details on the
[verification page](../verification/stage6.md#manufactured-solutions-of-every-model-spec-92-pr-148)).**
SST, standard and realizable k-ε, k-ε-ζ-f (2004) and the elliptic-blending
k-ε in a frozen cellular flow with every source term active: the L2
orders of all variables lie in 1.82–2.29 on hexahedra, polyhedra and
triangular prisms (bands met). Finding: the realizable C_μ now takes the
deviatoric strain (W̃ jumped where S → 0). Not covered: ζ-f₀ (its sink
is changing), the curvature correction and the QCR.

### 9.3 Zero-pressure-gradient flat plate

Set-up of the TMR flat plate verification case (2D, incompressible form of
its M = 0.2 case): Re = 5·10⁶ per unit length, plate from x = 0 to
x = 2 with a symmetry stretch ahead of it, the TMR's free-stream turbulence
values. The built-in mesher reproduces the TMR family of grids (the same
point counts and stretching, not the TMR files). Temperature as a passive
scalar with Pr = 0.71, isothermal plate.

| Quantity | Reference | Verification | Validation band B_model |
|---|---|---|---|
| C_f at x = 0.97 (SST, low-re treatment) | the Coles–Fernholz fit C_f = 2[ln(Re_θ)/0.384 + 4.127]⁻² (Nagib, Chauhan & Monkewitz 2007, a fit to measured layers) at the Re_θ of the run, and the Richardson-extrapolated value of our own mesh sweep (verification) | p in [1.5, 2.5], GCI_fine ≤ 1 % | ±8 % (decision 2026-11-22, D201: no value of another solver; the model itself lies about 6 % below White) |
| C_f(Re_x), 10⁶ ≤ Re_x ≤ 10⁷, every model with its default treatment | White (2006): C_f = 0.455/ln²(0.06 Re_x), and the Coles–Fernholz fit above | GCI_fine ≤ 2 % at x = 0.97 | ±5 % (White, k-ε family and ζ-f); ±8 % for SST against White |
| u⁺(y⁺) at Re_θ ≈ 10⁴ (y⁺ ≈ 1 meshes) | Spalding's law for 1 ≤ y⁺ ≤ 30, log law for 30 ≤ y⁺ ≤ 300 | – | ±3 % (SST, ζ-f, EB) |
| St_x(Re_x), 5·10⁵ ≤ Re_x ≤ 5·10⁶ | St = 0.0287 Re_x^{−0.2} Pr^{−0.4} (Kays, Crawford & Weigand 2005) | GCI_fine ≤ 2 % | ±10 % |
| y⁺ insensitivity: first-cell y⁺ of 0.5, 1, 2, 5, 10, 30, 60, 100 at x = 0.97 | the model's own y⁺ ≈ 1 result (low-re) | – | C_f within ±2 % (SST automatic), ±3 % (ζ-f and EB compound), ±5 % (k-ε compound, and scalable for y⁺ ≥ 30); St within ±5 % (Kader) |
| Smooth limits against hard limits (2.6) | the hard-limit result | – | C_f within 0.1 % |

**Results of 6b (PR #68, SST; details on the
[verification page](../verification/stage6.md#k-ω-sst-part-6b)).**
C_f(0.97) with `low-re` on 256 × 192 cells is 2.7005e-3, +1.45 % from
the Coles–Fernholz fit at the run's Re_θ (the TMR value used first was
dropped, D201): the validation row passes. Its convergence over the three
finest meshes is monotonic with p = 0.43 and GCI_fine = 1.87 %, outside
the verification band (finding; a fifth mesh decides whether the
asymptotic range is reached). C_f(Re_x) lies 5 to 9 % below White's
correlation (−6.1 % at x = 0.97) and −1.4 to +2.5 % from the
Coles–Fernholz fit: the ±5 % band against White does not suit SST and was
widened to ±8 % for it.
u⁺ at x = 0.97 follows Spalding's law within 2.3 %. The y⁺ study of
`automatic` gives +1.9 % (y⁺ 1.7), +15.2 % (4.6), +14.4 % (9.3),
−0.7 % (27) and −1.1 % (92), against our own value: the buffer layer misses the ±2 % band
(finding); TURB001 warns there and the blending is revisited with the
wall functions of 6c. *6g, PR #112:* the automatic treatment now blends
ω_P with the exponent 1.2 and imposes it, and the equilibrium production of
Reichardt's law, with the weight tanh((y⁺/5)²): C_f lies within −2.6 to
+1.2 % of our own finest-mesh result for first-cell y⁺ 0.24 to 85 on two meshes (the
sweep `examples/studies/sst_wall_sweep.py`; the ±2 % band is met for
y⁺ ≥ 4 on level 2). The temperature rows wait for 6f, the other models
for 6c–6e.

**Results of 6c (PR #71, k-ε family; details on the
[verification page](../verification/stage6.md#k-ε-family-part-6c)).**
TMR plate on 64 × 48 cells with the sustaining terms at the inlet values
(3.2 notes), the wall grading changed for each first-cell y⁺ at
x = 0.97; C_f(0.97) against White's 2.875·10⁻³:

| First-cell y⁺ | 1 | 2 | 5 | 10 | 30 | 60 | 100 |
|---|---|---|---|---|---|---|---|
| standard k-ε, `scalable` | −2.6 % | −1.7 % | −1.0 % | −0.3 % | −1.6 % | −1.5 % | −0.8 % |
| standard k-ε, `compound` | +58 % | +52 % | +40 % | +21 % | −2.3 % | −3.6 % | −3.4 % |
| realizable, `scalable` | −10.6 %¹ | | | | −6.8 %¹ | −5.2 % | −3.3 % |
| realizable, `compound` | | | | | −7.4 % | −7.1 % | −5.9 % |

¹ not converged to the tolerance (residual stalls near 10⁻⁵). The
C_f(Re_x) row passes for the standard model with both treatments in their
claimed ranges; the y⁺ insensitivity of `scalable` (all y⁺ within 2.4 % of
its y⁺ ≈ 1 result) and of `compound` from y⁺ 30 passes the ±5 % band.
The realizable model lies 3 to 7 % below White (finding, as in the
channel). The GCI of the k-ε plate follows with the model comparison of
6g.

### 9.4 Fully developed channel flow

Periodic plane channel driven by a mean pressure gradient that fixes
Re_τ (D29 *2026-11-15*), one cell column (2D, symmetric or full height),
refined in y only (four levels, r = 2); for the low-re models the finest
y⁺ of the first cell is 0.25. Re_τ = 180, 395 and 590 against the DNS of
Moser, Kim & Mansour (1999). Passive temperature with uniform heat flux
(the jump of D109 (1)), Pr = 0.71, against the DNS of Kasagi, Tomita &
Kuroda (1992, Re_τ = 150) and Kawamura et al. (1998, Re_τ = 180; Re_τ = 395 from the later DNS of the same group, chosen at implementation).

| Quantity | Models | Verification | Validation band |
|---|---|---|---|
| Bulk velocity U_b⁺ (equivalently C_f) | all | GCI_fine ≤ 0.5 % | ±3 % (SST, ζ-f, EB, realizable and standard k-ε with compound at y⁺ ≈ 1); ±5 % (k-ε family with scalable wall functions at y⁺ ≈ 30, Re_τ 395 and 590 only) |
| Centreline velocity U_c⁺ | all | GCI_fine ≤ 0.5 % | ±4 % |
| u⁺(y⁺), 5 ≤ y⁺ ≤ Re_τ | SST, ζ-f, EB | – | maximum deviation ≤ 5 % |
| k⁺ peak value and position; v²⁺ = ζk⁺ or φk⁺ near-wall profile | ζ-f, EB | – | peak position within ±20 %; reported against the model papers' own channel curves (the paper result is the acceptance target for EB at Re_τ = 395, section 3.4) |
| Nu (heat-flux jump), constant Pr_t and Kays–Crawford | all, Re_τ 150/180/395 | GCI_fine ≤ 1 % | ±7 % |
| Wall-function insensitivity: first-cell y⁺ from 1 to 100 | compound and automatic | – | wall shear within ±3 % of the y⁺ ≈ 1 result (stage-6 row of the plan) |

**Results of 6b (PR #68, SST, Re_τ = 395 only).** With `low-re`, 96
cells across and first-cell y⁺ 0.28: Re_τ = 397.8 from the driving force,
U_b⁺ = 17.43 (−0.7 % from the DNS, +1.1 % from Dean's correlation),
U_c⁺ = 19.64 (about −2.4 %); both inside their bands. On three meshes (96, 192, 384 cells
across) U_b⁺ = 17.427, 17.381, 17.346: observed order 0.40, extrapolated
17.237, GCI_fine 0.78 %, outside the verification band (finding; the
fourth mesh did not reach its tolerance). Re_τ = 180 and 590, the profiles and the wall-function
insensitivity follow with the model comparison of 6g.

**Results of 6c (PR #71, k-ε family, Re_τ = 395; details on the
[verification page](../verification/stage6.md#k-ε-family-part-6c)).**
U_b⁺ from the driving force on uniform meshes of 4 to 40 cells across
(first-cell y⁺ 93 to 9):

| Model, treatment | y⁺ 93 | 47 | 32 | 24 | 16 | 11–9 |
|---|---|---|---|---|---|---|
| standard k-ε, `scalable` | +6.0 % | +4.6 % | +4.3 % | +4.3 % | +4.5 % | +4.9 % |
| standard k-ε, `compound` | +5.9 % | +4.2 % | +3.1 % | +1.6 % | −4.4 % | −12.9 % |
| standard k-ε, `standard` | as `scalable` from y⁺ 24 | | | | −7.0 % | −25 % (y⁺ 13) |
| realizable, `scalable` | +7.1 % | +6.0 % | +6.5 % | +7.3 % | +9.2 % | +10.7 % |

The ±5 % row of the wall functions at y⁺ ≈ 30 passes for the standard
model (+4.3 %, `channel_keps.toml`). The realizable model misses it by
1.5 % (finding; on the flat plate its C_f is 3 to 7 % low as well, the
same direction). The ±3 % row of the compound treatment at y⁺ ≈ 1 does
not hold for the high-Re k-ε models (−32 % on a graded mesh at y⁺ 0.4,
the finding of 4.3); it moves to ζ-f and EB.

**Results of 6e (PR #88, elliptic-blending k-ε, `compound`, first-cell
y⁺ 0.13 to 0.28; details on the
[verification page](../verification/stage6.md#elliptic-blending-k-ε-d133)).**
C_f/C_f,DNS from the driving force at the DNS mass flow: 97.8 % at
Re_τ 180, 99.0 % at 395 (98.8 % on 192 cells), 99.1 % at 590, against
the paper's 100, 98 and 100 % (its Table 6); U_b⁺ +1.1, +0.5, +0.5 %.
The acceptance target of section 3.4 is met. At Re_τ 180 the start
(I 0.05, ν_t/ν 100) relaminarizes the channel; I 0.15 and ν_t/ν 20
keep it turbulent. Flat plate (9.3): C_f(0.97) +0.5 % from White.

**Results of 6g (PR #111, GCI of every model on three meshes; details
on the [verification page](../verification/stage6.md#gci-of-the-turbulence-models-on-the-plate-and-the-channel-93-94-pr-111)).**
Plate: no model shows an asymptotic order on levels 0 to 2 (SST, ζ-f and
EB oscillatory, k-ε and ζ-f₀ divergent, the realizable k-ε p = 1.0); SST
has GCI 0.47 % (levels 0 to 2); ζ-f does not
converge on the finer plates (C_f +16 to +24 %, finding). Channel: GCI
0.13 to 0.79 % for the low-Re models, the realizable k-ε with scalable
wall functions +11 % on the finest mesh (finding). Levels 1 to 3 belong
on the maintainer's machine.

### 9.5 Backward-facing step (Driver & Seegmiller 1985)

Expansion ratio 1.125, Re_H ≈ 36 000 on the step height H, inflow boundary
layer matched to the experiment's (δ ≈ 1.5 H at x/H = −4 ◇), as in the TMR's
2D backward-facing step case (its grid family reproduced by the built-in
mesher, three to four levels). Reattachment length x_R/H = 6.26 ± 0.10 from
the experiment.

| Quantity | Models | Verification | Validation band |
|---|---|---|---|
| x_R/H | SST, ζ-f, EB | GCI_fine ≤ 2 %, p in [1, 3] | ±10 % of 6.26 |
| x_R/H | standard and realizable k-ε | GCI_fine ≤ 2 % | ±5 % of a published result of the same model on a comparable grid (selected and recorded in the case before the run; the standard k-ε is known to reattach early). *6g:* standard k-ε with wall functions 5.22 and 5.47 from two solvers (Mizukami, Georgiadis & Cannon 1993, NASA Lewis, table 1), the band ±5 % about their mean 5.345; no published result of the realizable model on this case was found, its x_R is reported |
| C_p and C_f along the lower wall at the experiment's stations | all | – | reported against the experiment; RMS difference of C_p ≤ 0.02 for SST, ζ-f, EB |
| Velocity profiles at x/H = 1, 4, 6, 10 | all | – | reported (plots in the gallery) |

**Results of 6d2 (PR #74).** Built-in mesh with the step cut out
(`mesh.generate.remove`), wall from x = −85 H (boundary layer about 1.5 H
thick at x = −4 H), slip top at 9 H, outlet at 40 H, Re_H = 36 000;
x_r/H from the sign change of the wall shear, one mesh each (no GCI yet):
ζ-f₀ 6.38 (+2 %) on 8 200 cells, 6.38 on 32 800; with the published
realizability coefficient 0.6 7.44 / 7.67; ζ-f (2004, `compound`)
7.82 (+25 %; 6.80 with `a = 0.849`); SST (`automatic`) 5.79 (−7 %, on a Gmsh mesh of the same
layout). C_p, C_f and the profiles against the experiment follow with the
stage 6 validation.

**Results of 6e (PR #88).** The elliptic-blending k-ε on the mesh of
6d2 with SIMPLEC: x_r/H = 7.59 (+21 %; 7.68 on 32 800 cells), outside the band (finding). The
variable C*_ε2 of the paper causes it: with C_ε4 = C_ε2 (`c_epsilon4 =
1.83`) x_r/H = 6.24, but the flat plate is then +4.9 % from White
instead of +0.5 %. The model stays as published. The coupled algorithm
does not converge this case (residual of ε at 10⁻³, x_r 7.74 to 7.77;
7.84 on 32 800 cells).

**Results of 6g (PR #89; details on the
[verification page](../verification/stage6.md#backward-facing-step-against-the-experiment-part-6g-95)).**
The layout of the TMR: walls from x = −110 H, outlet at 50 H, a no-slip
top wall at 9 H (the slip top of 6d2 put the C_p of ζ-f₀ 0.033 RMS from
the experiment, the wall 0.018; x_r within 1 %); 11 040, 44 160 and
176 640 cells from the built-in mesher (the TMR's grids could not
be downloaded); x_r/H as the wall crossing `xr` with `skip = 2.0`, the corner
vortices of the finer meshes left out. On the finest mesh: SST 6.747
(+7.8 %), GCI below 0.01 % but an
oscillatory sweep (6.793, 6.747, 6.747), so p is not defined; k-ε 5.372,
+0.5 % from the published 5.345, GCI 1.94 %, p 0.69; realizable k-ε
6.407 (+2.4 %), its steps grow (no GCI); ζ-f₀ 6.296 (+0.6 %), monotonic, GCI below 0.01 %, p 4.8 (the fine pair 0.02 % apart). ζ-f and EB with
SIMPLEC on the two coarser meshes: 7.686, 7.741 (+23.7 %) and 7.446,
7.556 (+20.7 %), outside the band, as in 6d2 and 6e; their finest mesh
is a long run for the next 6g PR. C_p RMS on the finest mesh: SST 0.023,
k-ε 0.018, realizable 0.026, ζ-f₀ 0.017; ζ-f 0.037 and EB 0.034 on the
middle mesh: the band 0.02 is missed by SST, ζ-f and EB (findings). The
k-ε models run at first-cell y⁺ 1.5 to 0.3 with wall functions claimed
for y⁺ ≥ 30; a sweep that keeps the first cell may show whether that
causes their slow or growing steps.

### 9.6 Axisymmetric impinging jet, H/D = 2

Fully developed pipe jet (profile from a periodic pipe run, D111, fed in as a
table inlet with k and ε or ω columns, section 7.1) on a flat plate at
H/D = 2, axisymmetric on a wedge mesh (D29 *Part 5b*). Heat transfer:
Re = 23 750 with uniform wall heat flux, Baughn & Shimizu (1989) and Baughn
et al. (1991); velocity: Re = 23 000, Cooper, Jackson, Launder & Liao (1993).
Temperature as a passive scalar with Pr = 0.71 and uniform flux. The
digitised reference data (Nu(r/D); mean velocity and turbulence profiles at
the published r/D stations) are stored in the case with their source,
figure or table, and retrieval note.

| Quantity | Models | Verification | Validation band |
|---|---|---|---|
| Stagnation Nusselt number Nu₀ | ζ-f, EB | GCI_fine ≤ 3 % | ±15 % |
| Nu₀ | SST with `kato-launder` or `durbin` | GCI_fine ≤ 3 % | ±20 % |
| Secondary maximum of Nu(r/D) | ζ-f, EB | – | present, located at r/D within ±0.5 of the experiment's |
| Nu(r/D), 0 ≤ r/D ≤ 4 | ζ-f, EB, SST with limiter | – | RMS difference ≤ 15 % |
| Stagnation-point anomaly | standard k-ε without and with Kato–Launder; SST `none` against `menter`, `kato-launder`, `durbin` | – | without limiter Nu₀ over-predicted by more than 50 % (the anomaly is reproduced); the limiter removes at least 70 % of the excess k at the stagnation point. This tests the limiters, not the models |
| Wall-jet velocity profiles at r/D = 0.5, 1, 2, 2.5, 3 | all | – | reported |

**Results of 6g (PR #90, with the fix of PR #91; details on the
[verification page](../verification/stage6.md#impinging-jet-at-hd--2-part-6g-96)).**
`impinging_jet.toml` on a 0.5° wedge, the inlet the table of a periodic
pipe of each model (k, ε or ω and ζ columns, 7.1), 2 600, 10 400 and
41 600 cells. Nu₀ is the wall probe `stag` at the axis, next to it `r01`
at r/D = 0.1. The first runs had the face next to the axis low with
every turbulence model and more so on finer meshes (ζ-f₀ −7 % on the
finest mesh), while the laminar jet was flat to the axis: the symmetry
faces lacked the transposed part of the stress, on a wedge half of the
hoop stress (2.4; fixed in PR #91). Now Nu is flat to the axis, and the
GCI of Nu₀ at the axis is 0.12 % for ζ-f (173.8), 2.9 % for SST with
`kato-launder` (141.2) and 0.19 % for EB (249.2); at r/D = 0.1: ζ-f
173.9 (GCI 0.03 %), EB 205.0 (0.40 %), SST 141.1 (0.31 %), ζ-f₀ 151.7
(5.3 %, p 0.74: outside 3 %, not yet asymptotic). Secondary maximum:
ζ-f 118.4 at r/D 2.27, EB 93.4 at 2.31 (weak), ζ-f₀ 128.7 at 1.95.
Anomaly on 10 400 cells: k-ε `none` Nu₀ 210.7 (+34 % above ζ-f₀) with
6.6 times the axis k; `durbin` removes 71 % of the excess,
`kato-launder` 111 %. SST `none`, `menter` and `durbin` give the same
field (Nu₀ 160.5, the axis k below that of ζ-f₀): its bound
ν_t = a₁k/(F₂S) holds on the axis, so there is no anomaly for the
limiters to remove; `kato-launder` 144.6. The thermal law of the k-ε
family now takes u* (5.2; Nu₀ 37 with u_τ). **Against the experiment** (ERCOFTAC classic database, case 025, files
`ij2lr-*`: Nu/Re^0.7 of Baughn et al. scaled to Re 23 750, so Nu₀ 141.4
and the secondary maximum 108.8 at r/D 1.94): Nu₀ of SST with
`kato-launder` (−0.2 %) and ζ-f₀ (+7.2 %) is inside its band, that of
ζ-f (+22.9 %) and EB (+76 %) outside ±15 %; the RMS difference of
Nu(r/D) for 0 ≤ r/D ≤ 4 is 13.7 % for SST (met) and 21 to 32 % for
ζ-f₀, ζ-f and EB (not met); the secondary maximum lies within ±0.5 D of
the experiment's for ζ-f, EB and ζ-f₀; the standard k-ε with `none` is
+49 % (the band asks for more than +50 %: missed by 1 point on 10 400
cells). The wall-jet maxima, the findings and the scaling are on the
verification page (finding 8). The data are not in the repository (the
wiki of the database is licensed CC BY-NC-SA 4.0: non-commercial,
share-alike; part 4).

### 9.7 Cases kept from the earlier stage-6 scope

Unchanged from the stage plan: rough-pipe friction factor against
Colebrook (periodic pipe, fully rough and transitional, ±5 %; *6c*
results below), the study
over Re and k_s/D that reproduces the Moody chart (D62), laminar
Taylor–Couette flow in a rotating frame (D55) against the analytic
solution, the turbulent heated pipe as a passive scalar against Gnielinski
(1976) (±10 %, 10⁴ ≤ Re ≤ 10⁶, D111) with each model and both Pr_t
options, and the CPU–GPU and thread-count parity of section 6.6.

**Results of 6c (PR #71): rough pipe and Moody study.** Standard k-ε,
`scalable`, wedge pipe of one periodic module, `pipe_rough_keps.toml` and
`examples/studies/moody_study.py`; Darcy f from the driving force against
Colebrook:

| Re | smooth | k_s/D = 10⁻⁴ | 10⁻³ | 10⁻² |
|---|---|---|---|---|
| 3·10⁴ | −5.0 % | −6.1 % | −14.1 % | −10.0 % |
| 10⁵ | −3.8 % | −6.5 % | −17.3 % | +9.6 % |
| 10⁶ | −2.7 % | −13.5 % | +4.4 % | +11.0 % |
| 10⁷ | −1.7 % | +2.0 % | +7.9 % | +11.1 % |

Smooth pipes pass from Re 10⁵ (the wall functions of a 6-cell mesh at
Re 3·10⁴ sit at the band's edge). Fully rough (k_s⁺ > 90): +8 to +11 %
with the default C_s = 0.5, +1 to +3 % with C_s = 0.4 (finding of 4.5).
Transitional (k_s⁺ 3 to 50): the Cebeci–Bradshaw shift follows
Nikuradse's inflectional sand-grain curve and lies up to 17 % below
Colebrook's monotonic interpolation for commercial pipes; the ±5 % band
against Colebrook cannot hold there by construction of ΔB (finding). The
chart is gallery 71.

**Results of 6g part 5 (PR #94): the moving reference frame (D55).**
`[[mrf]]` turns a cell zone steadily about an axis. The solved velocity
is the absolute one **u**; the convecting velocity of the zone is
**u**ᵣ = **u** − **Ω** × (**x** − **x**₀):

    ∇·(ρ **u**ᵣ ⊗ **u**) + ρ **Ω** × **u** = −∇p + ∇·τ(**u**)

(a field that is steady in the rotating frame has ∂**u**/∂t = **Ω** × **u** −
(**Ω** × **r**)·∇**u** in the inertial one). Turbulence variables and passive
scalars φ, which are scalars of the frame, are convected with ρ **u**ᵣ; the
strain rate is the same in both frames, the vorticity of the Kato–Launder
production is the absolute one (TURB004). Discretely: continuity, the
Rhie–Chow flux and the pressure correction work on the absolute mass flux ṁ;
the convecting flux of momentum, turbulence and scalars is ṁ − ρ **S**·**V**
with **V** = **Ω** × (**x**_f − **x**₀) in a zone, the mean of the two cells'
values on an interface, the value of the cell on a boundary face (a wall
that turns with the frame has the absolute flux ρ **S**·**V** and so zero
relative flux); the Coriolis source −ρ V_c **Ω** × **u**_c of the rotating
cells is explicit. Σ_f ρ **S**·**V** = 0 for every closed cell with planar
faces, so the relative flux conserves mass wherever all faces carry their
exact frame flux: inside a zone, on boundaries and on interfaces that are
surfaces of revolution about the axis, the limit of every frozen-rotor
method; MRF003 reports the cells where the closure fails. The method was
steady only at first; D185 adds the transient frozen rotor (the time term
of the velocity in the rotating basis), and the moving rotor is the
sliding mesh (D182).
Laminar Taylor–Couette flow (r₁ = 0.05, r₂ = 0.10 m, ω₁ = 10 rad/s, Re = 50,
the inner half in the frame): L2 error of the velocity 3.46, 0.69 and 0.13 %
on 256, 1 024 and 4 096 cells, observed order 2.32 and 2.41, and on 1 024
cells 0.63 % without a frame, 0.69 % in the frame of the cylinder, 0.57 %
in the opposite frame and 0.58 % at 2.5 times the speed. Details on the
[verification page](../verification/stage6.md#moving-reference-frame-part-6g-97).

**Results of 6g part 6 (PR #95): the sliding mesh with conformal stepping
(D56, variant A of D182).** `[[sliding]]` turns the rotor of a ring of N
identical sectors by k sectors in every time step, Δt = k · 2π/(N |Ω|)
(`SLD003` names the value; the adaptive step and the error tolerance are
refused). The solved velocity is the absolute one, momentum and scalars are
convected with the flux relative to the moving mesh, ṁ − ρ **S**·**V** as in
the frame of D55 but without the Coriolis source (the basis does not turn);
Σ_f ρ **S**·**V** = 0 for every closed cell and the interface at the middle of
the gap is a surface of revolution, so the grid flux of the interface
vanishes and the geometric conservation law holds exactly. Implementation:
the solver stays a fixed-mesh solver, and all fields, gradients, the Rhie–Chow
velocity and the extrapolations are Eulerian on the lattice; only the time
term sees the motion: the old levels n and n − 1 of a rotor cell are those of
the cell that stood at its place one and two steps earlier (the permutation π
of the cells that an advance by k sectors induces, and π∘π), and the vectors
of these levels are not rotated (the pressure gradient turns them; BDF2
leaves an O(Δθ⁴) error per step). Shifting gradients, fluxes or the
extrapolated pressure gradient with the data was tried first and is unstable.
Laminar Taylor–Couette flow from rest (same geometry as above, t = 10 s):
L2 error of the velocity 6.76, 1.35, 0.28 and 0.06 % on 128, 512, 2 048 and
8 192 cells with the time step halved each time, observed order 2.3, 2.3 and
2.2; two sectors per step 2.7 % on 512 cells; solid-body rotation of the
whole mesh kept at 1.7 % (0.6 s) and 1.6 % (3 s); a passive scalar through the
turning rotor 0.7 %. Laminar flow and passive scalars only (turbulence models
were refused, `SLD004`, until PR #100), one rotor, a fixed time step, CPU. Details on the
[verification page](../verification/stage6.md#sliding-mesh-with-conformal-stepping-part-6g-97).

**Results of 6g part 7 (PR #96): the sliding mesh by imprinting (variant B of
D182).** A rotor and a stator with different numbers of sectors, and a free
time step. `[mesh.annulus] n_theta_stator` makes the two rings with the boundaries
`rotor_interface` and `stator_interface`, `[[sliding]]` names them. In every
step the rotor nodes are turned to the angle Ω t and the interface is imprinted:
every rotor face is clipped with every stator face in the unrolled (θ, z) plane of
the interface cylinder (convex polygons, Sutherland–Hodgman, `foehncraft_mesh::imprint`),
every intersection polygon is an internal face between a rotor and a stator cell,
and the faces of a cell add up to a closed surface (Σ S = 0 for extruded meshes),
so mass and momentum are conserved across the interface to round-off without
interpolation weights. Volumes and centroids are the rigid-body values, the
interface is a surface of revolution (w·S = 0), so the geometric conservation law
holds exactly. Implementation: the solver is built on the mesh of the angle in
every step (`Flow::on`, scalars and frame attached, `Transient::restore`), the
state and the old levels are carried (cell arrays are the same, faces of the
base mesh keep their flux and the pieces start from the flux of the linearly
interpolated velocities), and the step is fully implicit (`extrapolation = false`,
at least two outer iterations: the vectors of the previous angle are stale on
the turned mesh). D182 foresaw rebuilding only the interface part of the face list,
the matrices and the pressure hierarchy; the whole solver is rebuilt in every step
(PR #97 made the rebuild and the setup cheaper: the case file 8.4 s to 5.6 s, the
setups a quarter of its solver time; keeping the hierarchy of the step before
was measured and does not pay, the rotor turns by about a cell per step), a
rebinding in place is solver-speed work for later. Laminar Taylor–Couette flow from
rest with 16 + 24, 32 + 48 and 64 + 96 sectors (Δt halved with the mesh): L2 error
6.74, 1.29 and 0.25 %, observed order 2.4 and 2.35; equal pitch 1.350 % against
1.354 % of conformal stepping; solid-body rotation 1.02 % (0.6 s) and 0.97 % (3 s);
a passive scalar 0.49 %; mass balance 10⁻¹⁸. Laminar flow and passive scalars only
(turbulence since PR #100), one rotor, a fixed time step, no restart, CPU. Details on the
[verification page](../verification/stage6.md#sliding-mesh-by-imprinting-part-6g-97).

**Results of 6g (PR #100): turbulence on the frame and both sliding meshes.**
Every turbulence model runs with `[[sliding]]` and `[[mrf]]`. The rotor
carries k, ε, ω, ζ and f with their old levels (conformal stepping: moved by
the sector permutations, each step starting from the moved fields; imprinting:
carried through the rebuild of every step), and the wall distance of the
turning rotor is updated in every step from one bounding volume hierarchy per
rigid part, searched only from the cells nearer to the edge of their part
than to its walls. Turbulent Taylor–Couette flow (SST, Re 12 500): the frame
against the inertial run 0.24 / 0.43 / 0.34 % (u_θ, k, ν_t, wall-resolved
mesh); the spin-up on both sliding meshes against the fixed mesh 1.6 / 1.5 /
4.0 %, about 0.5 / 1.0 / 1.5 % on the next mesh; a ζ-f patch carried half a turn by the
rotor. Two independent frames (10 and −4 rad/s) 0.78 % against 0.72 %
inertial. The turbulent stirred vessel (water, Re 20 000) found a NaN
hierarchy of the pressure multigrid (fixed). Details on the
[verification page](../verification/stage6.md#turbulence-on-the-frame-and-the-sliding-meshes-part-6g-pr-100).

### 9.8 The SST options of 3.6 (D170)

Bands fixed before the runs, as everywhere in this section; the ◇ bands
are proposals to be confirmed with the papers.

| Case | Quantity | Reference | Band |
|---|---|---|---|
| Flat plate and channel at Re_τ = 395 with `curvature_correction` | C_f, u⁺ profile | SST without the option (pure shear: f_r1 = 1) | relative difference ≤ 1e-6 |
| Channel in spanwise system rotation, Re_τ = 194 (rotating zone, D55; periodic, 6p) | wall shear stress on the unstable and the stable side, mean-velocity profile | DNS of Kristoffersen & Andersson (1993) ◇ | ratio of the two wall shear stresses within ±15 % ◇ |
| Fully developed square duct, Re_b ≈ 3 500 (periodic, 6p) with `qcr = "qcr2013"` | largest secondary velocity, its position | DNS of Pinelli, Uhlmann, Sekimoto & Kawahara (2010) ◇ | within ±30 % and the corner vortex pair present; linear SST gives none (a check that the option acts) |
| Manufactured solution with the corrections switched on | observed order | — | 1.8 – 2.2 (the corrections must not lower the order) |


## 10. Implementation plan

Pull requests in order, each with tests, gallery entry, case files, recipe
(D18), changelog fragment and the verification-manual update; each starts
only on the maintainer's go. Branch names from the task (project rule).

| PR | Content | Main files | Acceptance |
|---|---|---|---|
| 6a | Generic scalar transport (D116), passive scalars and their boundary sub-tables (D51), bounded BDF2 (D28 (1)), smooth-limit helpers (2.6), exact wall distance with the BVH (D48, D108 (1)) and the derived field `wall-distance` (done, PR #63) | `foehncraft-kernels/src/scalar.rs`, `foehncraft-mesh/src/query/bvh.rs`, `foehncraft-flow/src/scalar.rs`, `foehncraft-core/src/math.rs` | 9.2 laminar heated channel (developed profile at the inlet, uniform wall flux, zero-gradient outlet; Nu within 0.5 %), wall distance, scalar MMS on hex, poly and tet cells, scalar time order |
| 6p | Periodic boundaries for the flow and the scalars (translational pairs, mass-flow or pressure-gradient driving, the heat-flux jump of D109 (1) for passive scalars; the wall-temperature map of D109–D111 with the energy stage) (done, PR #67) | `foehncraft-mesh` periodic pairing, `foehncraft-flow/src/bc.rs` | periodic Poiseuille and heated channel equal to the open ones; needed by the 9.4 channel at Re_τ 180/395/590 |
| 6b | k-ω SST 2003 with the three production options, low-re and automatic treatments, rough walls (Aupoix), inlet turbulence, TURB codes, log columns, p̂; the optional PDE wall distance with the key `turbulence.wall_distance` (the exact distance of 6a is the default); the invariant features of D166 as named expression variables (same definitions as the D106 export) (done, PR #68: the features as the derived field `turbulence-features`; their names in the expression language come with the first consumer, D104; the SST MMS moves to 6g, done for every model in PR #148; turbulence values of the backflow at outlets later, BC028, done in PR #140) | `foehncraft-kernels/src/turb/{mod,invariants,sst,wall}.rs`, `foehncraft-flow/src/turbulence/{mod,wall,inlet}.rs` | 9.2 decay (SST), 9.3 flat plate, 9.4 channel (SST) |
| 6c | Standard and realizable k-ε, scalable and standard wall functions, compound treatment for the k-ε family, Cebeci–Bradshaw roughness (done, PR #71: compound claimed from y⁺ 30 for the k-ε family, owner decision; Durbin bound and Picard ε destruction as defaults; findings: realizable 3–7 % low friction, rough pipe against Colebrook in 4.5 and 9.7) | `turb/keps.rs`, `turb/wall.rs` | 9.2, 9.3, 9.4 for the k-ε family, rough pipe and Moody study |
| 6d | k-ε-ζ-f with the compound treatment (Popovac & Hanjalić 2007; Popovac 2023 refinements as decided then), elliptic solver path (CG + AMG), the AMG preconditioner for the scalar equations (D36; 6a solves them with Jacobi-preconditioned BiCGStab, enough for the laminar cases) (done, PR #73: channel +1.9 % / +1.4 %; findings on the plate: +5 % C_f, laminar from the TMR free stream, the k floor near the leading edge; the elliptic solve lives in `turbulence/mod.rs`) | `turb/zeta_f.rs`, `foehncraft-flow/src/turbulence/elliptic.rs` | 9.2, 9.3, 9.4 incl. near-wall profiles |
| 6d2 | ζ-f₀ of Popovac (2023) with its compound treatment (maintainer decision 2026-11-17: robustness in y⁺ for industrial use) (done, PR #74: no code of the author available, own default of 3.3.1; channel +2.7 to −4.2 % and plate within 2.5 % of White over y⁺ 0.3 to 30, all converged; backward-facing step x_r/H 6.38 (+2 %; ζ-f +25 %, SST −7 %); the built-in mesher removes boxes, `mesh.generate.remove`; #141: the wall sink of k of 3.3.2 and C_ε1 = 1.4 (1 + 0.009/ζ), channel mesh-converged at +0.2 % from the DNS) | `turb/zeta_f.rs`, `turb/wall.rs` | 9.3 and 9.4 over first-cell y⁺ 0.5 to 100, same or similar results as ζ-f |
| 6e | EB k-ε (D133) (done, PR #88: channel C_f at 97.8 / 99.0 / 99.1 % of the DNS at Re_τ 180 / 395 / 590, plate +0.5 % from White; step x_r/H 7.59, +21 %, finding: the variable C*_ε2, 6.24 with it constant; the coupled algorithm does not converge the step; wall normal not yet exported) | `turb/eb.rs` | 9.4 at Re_τ = 395 against the DNS first, then 9.2–9.4 |
| 6f | Thermal wall treatment and Pr_t (D135) on passive scalars (done, PR #85; the default wall law is `kader` also with `low-re`, equal to `none` there; the channel Nu against DNS is open) | `turb/thermal.rs` | 9.2 wall-function functions, St on the flat plate, channel Nu, heated pipe against Gnielinski |
| 6g | Validation cases: backward-facing step, impinging jet; MRF (D55) with Taylor–Couette, then the sliding mesh in both variants (D182); statistics (D71) and derived fields (D72) of the turbulence models; GPU parity notebook; GCI reports (in several PRs; part 1 done, PR #89: the backward-facing step in the TMR layout, six models on three meshes, the reattachment point as a wall crossing `[[output.wall_crossings]]` with its GCI; part 2 done, PR #90: the impinging jet, the inlet table with turbulence columns, wall probes with scalars, the thermal law in u*, the comparison with the experiment once its data are in, made in part 4; part 3 done, PR #91: the full stress on symmetry and wedge planes, the axis cells of the jet; part 4 done, PR #92: the jet against the ERCOFTAC data of case 025, the comparison with the bands of 9.6, `ercoftac_to_csv.py`; part 5 done, PR #94: the moving reference frame (D55), `[[mrf]]`, `[mesh.annulus]`, laminar Taylor–Couette flow against the analytic solution in any frame, `crates/foehncraft-flow/src/frame.rs`, `tests/mrf.rs`, `taylor_couette_mrf.toml`; part 6 done, PR #95: the sliding mesh with conformal stepping (D56, variant A of D182), `[[sliding]]`, `crates/foehncraft-flow/src/sliding.rs`, `tests/sliding.rs`, `taylor_couette_sliding.toml`; part 7 done, PR #96: the sliding mesh by imprinting (variant B of D182), `[[sliding]] rotor_interface`, `[mesh.annulus] n_theta_stator`, `crates/foehncraft-mesh/src/imprint.rs`, `taylor_couette_imprint.toml`; PR #97: the rebuild of every step and its multigrid setup made cheaper; PR #98: walls of a turning zone move normal to themselves with it, the stirred vessel with blades and baffles in both methods, field figures, `examples/studies/mixer.py`, `tools/plot/rotor_fields.py`; PR #99: the fields of rotating machinery of D183, relative and frame velocity, cylindrical components, flow angles, relative total pressure, rothalpy; PR #100: turbulence on the frame and both sliding meshes, the wall distance of a turning rotor `crates/foehncraft-mesh/src/query/turning.rs`, the turbulent stirred vessel; PR #108: the derived fields of D72 for the turbulence models, `turbulence-scales`, `reynolds-stress`, `wall-normal`, `turbulence-wall`, spec 2.3 except the blending weight, which PR #130 adds as the dataset `blending` of `turbulence-wall`; PR #109: the GPU parity notebook `tools/notebooks/stage6_colab.ipynb` with `tools/bench/stage6_cloud.py`, 14 cases, the T4 run by the maintainer; PR #110: the statistics of D71, `[output.statistics]`, time means, RMS, covariances and phase averages of fields, probes and coefficients, batch means of whole periods with the von Neumann test, `cylinder_re100_statistics.toml`, `tests/statistics.rs`) | cases, gallery, `tools/notebooks/stage6_colab.ipynb`; `examples/cases/flow/backstep_validation.toml`, `examples/studies/backstep_validation.py`; `examples/cases/flow/{jet_pipe_profile,impinging_jet}.toml`, `examples/studies/impinging_jet.py`, `python/foehncraft/profiles.py`; `foehncraft-kernels/src/flow.rs`, `tests/symmetry_stress.rs` | 9.5, 9.6, 9.7, 6.6 |
| 6h | SST options of D170: curvature and rotation correction (Smirnov & Menter 2009), QCR2000 and QCR2013 with the explicit-stress path in the momentum assembly that D171 reuses (part 1 done, PR #113: the curvature correction, `turbulence.curvature_correction`, `flow.frame_rotation` for the whole domain in a turning frame, the derived field `curvature-correction`, `channel_rotating_sst.toml`, `tests/curvature.rs`; pure shear unchanged to the last digit, rotating channel at Ro 0.15 wall shear ratio 1.63 (SST 1.00), step x_r/H +1.9 %; the DNS comparison of 9.8 is open; part 2 done, PR #115: `turbulence.qcr` as an explicit stress on the momentum faces, `square_duct_qcr.toml`; linear SST no secondary flow, QCR2000 largest secondary velocity 1.56 % of U_b with the corner vortex pair; the DNS comparison is open) | `turb/sst.rs`, `turb/invariants.rs`, `foehncraft-flow/src/turbulence/stress.rs` | 9.8 |
| release | v0.6: version bump, changelog summary, README, verification manual, Release workflow (dry run first) | | all of section 9 green or findings recorded |

Validation cases shared with pyaeroflow, the maintainer's earlier
compressible code (it resolves walls, no wall functions), so that the two
codes can be compared on the same set-ups: the TMR flat plate (6b, done),
the heated flat plate (6f), the NASA wall-mounted hump (6g, and the
data-driven work of L37).

Order rationale: 6a gives every later model its transport, wall distance
and solver path; 6p brings the periodic boundaries that the channel cases
of 9.4 need, so it comes before the first model is verified there; SST first because it is the general-purpose default and
the flat plate and channel can be checked against the TMR at once; ζ-f
before EB because EB reuses its elliptic path; the thermal wall function
after the models because it needs u* from each of them; the SST options
(6h) after the validation cases, because their checks need the rotating
zone and the periodic duct. The EARSM of D171 is not in v0.6.


## 11. Questions for the maintainer

- **Q1 (closed 2026-11-17):** EB k-ε joins ζ-f in stage 6 (D133).
- **Q2:** Lag EB k-ε as P11 (possible later), decided after the EB
  verification (recommended), or in stage 6?
- **Q3:** The model bands of section 9 are a proposal; change any before
  the first run (they are fixed once a run has been made).
- **Q4 (open, not confirmed):** Add the k-ω model of Wilcox (2006) as a
  cheap variant on the SST code path (its stress limiter and cross
  diffusion; one more closure of the same kernels)? Not part of the plan
  until the maintainer confirms it.


## 12. References

- B. Aupoix, *Roughness corrections for the k–ω shear stress transport
  model: status and proposals*, J. Fluids Eng. 137, 021202 (2015)
- J. W. Baughn, S. Shimizu, *Heat transfer measurements from a surface
  with uniform heat flux and an impinging jet*, J. Heat Transfer 111,
  1096–1098 (1989)
- J. W. Baughn, A. E. Hechanova, X. Yan, *An experimental study of
  entrainment effects on the heat transfer from a flat surface to a heated
  circular impinging jet*, J. Heat Transfer 113, 1023–1025 (1991)
- F. Billard, D. Laurence, *A robust k-ε-v²/k elliptic blending turbulence
  model applied to near-wall, separated and buoyant flows*, Int. J. Heat
  Fluid Flow 33, 45–58 (2012)
- T. Cebeci, P. Bradshaw, *Momentum transfer in boundary layers*,
  Hemisphere (1977)
- I. B. Celik et al., *Procedure for estimation and reporting of
  uncertainty due to discretization in CFD applications*, J. Fluids Eng.
  130, 078001 (2008)
- D. Cooper, D. C. Jackson, B. E. Launder, G. X. Liao, *Impinging jet
  studies for turbulence model assessment — I. Flow-field experiments*,
  Int. J. Heat Mass Transfer 36, 2675–2684 (1993)
- D. M. Driver, H. L. Seegmiller, *Features of a reattaching turbulent
  shear layer in divergent channel flow*, AIAA J. 23, 163–171 (1985)
- P. A. Durbin, *On the k-ε stagnation point anomaly*, Int. J. Heat Fluid
  Flow 17, 89–90 (1996)
- L. Eça, M. Hoekstra, A. Hay, D. Pelletier, *On the construction of
  manufactured solutions for one and two-equation eddy-viscosity models*,
  Int. J. Numer. Methods Fluids 54, 119–154 (2007)
- V. Gnielinski, *New equations for heat and mass transfer in turbulent
  pipe and channel flow*, Int. Chem. Eng. 16, 359–368 (1976)
- G. Grotjans, F. R. Menter, *Wall functions for general application CFD
  codes*, ECCOMAS 98, 1112–1117 (1998)
- K. Hanjalić, M. Popovac, M. Hadžiabdić, *A robust near-wall
  elliptic-relaxation eddy-viscosity turbulence model for CFD*, Int. J. Heat
  Fluid Flow 25, 1047–1051 (2004)
- A. Hellsten, *New advanced k-ω turbulence model for high-lift
  aerodynamics*, AIAA J. 43, 1857–1869 (2005)
- C. L. V. Jayatilleke, *The influence of Prandtl number and surface
  roughness on the resistance of the laminar sub-layer to momentum and
  heat transfer*, Prog. Heat Mass Transfer 1, 193–329 (1969)
- B. A. Kader, *Temperature and concentration profiles in fully turbulent
  boundary layers*, Int. J. Heat Mass Transfer 24, 1541–1544 (1981)
- N. Kasagi, Y. Tomita, A. Kuroda, *Direct numerical simulation of passive
  scalar field in a turbulent channel flow*, J. Heat Transfer 114, 598–606
  (1992)
- M. Kato, B. E. Launder, *The modelling of turbulent flow around stationary
  and vibrating square cylinders*, 9th Symp. Turbulent Shear Flows, Kyoto,
  10-4 (1993)
- H. Kawamura, K. Ohsaka, H. Abe, K. Yamamoto, *DNS of turbulent heat
  transfer in channel flow with low to medium-high Prandtl number fluid*,
  Int. J. Heat Fluid Flow 19, 482–491 (1998)
- W. M. Kays, *Turbulent Prandtl number — where are we?*, J. Heat Transfer
  116, 284–295 (1994)
- W. M. Kays, M. E. Crawford, B. Weigand, *Convective heat and mass
  transfer*, 4th ed., McGraw-Hill (2005)
- R. Kristoffersen, H. I. Andersson, *Direct simulations of low-Reynolds-number
  turbulent flow in a rotating channel*, J. Fluid Mech. 256, 163–197 (1993)
- S. Lardeau, F. Billard, *Development of an elliptic-blending lag model
  for industrial applications*, AIAA 2016-1600 (2016)
- B. E. Launder, D. B. Spalding, *The numerical computation of turbulent
  flows*, Comput. Methods Appl. Mech. Eng. 3, 269–289 (1974)
- M. Mani, D. A. Babcock, C. M. Winkler, P. R. Spalart, *Predictions of a
  supersonic turbulent flow in a square duct*, AIAA 2013-0860 (2013) ◇
- F. R. Menter, *Two-equation eddy-viscosity turbulence models for
  engineering applications*, AIAA J. 32, 1598–1605 (1994)
- F. R. Menter, M. Kuntz, R. Langtry, *Ten years of industrial experience
  with the SST turbulence model*, Turbulence, Heat and Mass Transfer 4,
  625–632 (2003)
- R. D. Moser, J. Kim, N. N. Mansour, *Direct numerical simulation of
  turbulent channel flow up to Re_τ = 590*, Phys. Fluids 11, 943–945 (1999)
- NASA Langley Turbulence Modeling Resource (TMR): SST and SST-2003 pages,
  2D zero-pressure-gradient flat plate and 2D backward-facing step cases
  (values retrieved at implementation, retrieval date stored with the data)
- S. V. Patankar, *Numerical heat transfer and fluid flow*, Hemisphere (1980)
- A. Pinelli, M. Uhlmann, A. Sekimoto, G. Kawahara, *Reynolds number
  dependence of mean flow structure in square duct turbulence*, J. Fluid
  Mech. 644, 107–122 (2010)
- M. Popovac, K. Hanjalić, *Compound wall treatment for RANS computation of
  complex turbulent flows and heat transfer*, Flow Turbul. Combust. 78,
  177–202 (2007)
- M. Popovac, *Robust eddy viscosity turbulence modeling with elliptic
  relaxation and compound wall treatment*, Energies 16, 3685 (2023),
  doi:10.3390/en16093685 (CC BY; the ζ-f₀ model and the 2023 compound
  treatment of 6d2)
- T.-H. Shih, W. W. Liou, A. Shabbir, Z. Yang, J. Zhu, *A new k-ε eddy
  viscosity model for high Reynolds number turbulent flows*, Comput. Fluids
  24, 227–238 (1995)
- P. E. Smirnov, F. R. Menter, *Sensitization of the SST turbulence model
  to rotation and curvature by applying the Spalart–Shur correction term*,
  J. Turbomach. 131, 041010 (2009)
- P. R. Spalart, *Strategies for turbulence modelling and simulations*,
  Int. J. Heat Fluid Flow 21, 252–263 (2000)
- P. R. Spalart, C. L. Rumsey, *Effective inflow conditions for turbulence
  models in aerodynamic calculations*, AIAA J. 45, 2544–2553 (2007)
- P. R. Spalart, M. Shur, *On the sensitization of turbulence models to
  rotation and curvature*, Aerosp. Sci. Technol. 1, 297–302 (1997)
- D. B. Spalding, *A single formula for the "law of the wall"*, J. Appl.
  Mech. 28, 455–458 (1961)
- S. Wallin, A. V. Johansson, *An explicit algebraic Reynolds stress model
  for incompressible and compressible turbulent flows*, J. Fluid Mech. 403,
  89–132 (2000)
- F. M. White, *Viscous fluid flow*, 3rd ed., McGraw-Hill (2006)
- D. C. Wilcox, *Turbulence modeling for CFD*, 3rd ed., DCW Industries
  (2006) (Q4 only)
