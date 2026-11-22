# Stage 7 – Energy equation, buoyancy and conjugate heat transfer

Requirements step 19 (2026-11-17), asked for by the maintainer. This page is
the full specification; `docs/requirements.md` carries the decisions D150 to
D164 as short rows that point here, the step O18 (stage order and the
scope of milestone 1) and the changed rows of the stage plan. Thin walls
(D189) and the heat exchange of the wall film (D190) were added on
2026-11-20 (§9a, §9b, open step O24). No code is
written before the maintainer's go (project rule).

The rules of the project hold throughout: requirements before code; every
stage ends with results that can be checked against analytic solutions or
published data, with the discretisation uncertainty of Celik et al. (2008,
D67); clean-room numerics (D21): every equation, constant and interface
treatment below comes from the cited literature, and code from other solvers
is never copied or ported. Kernels stay `no_std`, float64 (D7), bitwise
reproducible (D31, D93), the same on CPU and GPU (D8, D60) and
adjoint-ready (D59).

Contents

1. [Place in the plan](#1-place-in-the-plan)
2. [What stage 7 builds on](#2-what-stage-7-builds-on)
3. [Energy equation (D150)](#3-energy-equation-d150)
4. [Low-Mach variable density (D151)](#4-low-mach-variable-density-d151)
5. [Turbulent heat flux (D152)](#5-turbulent-heat-flux-d152)
6. [Buoyancy (D153) and the balanced-force interpolation (D154)](#6-buoyancy-d153-and-the-balanced-force-interpolation-d154)
7. [Buoyancy in the turbulence models (D155)](#7-buoyancy-in-the-turbulence-models-d155)
8. [Thermal boundary conditions and their check (D156)](#8-thermal-boundary-conditions-and-their-check-d156)
9. [Fluid–solid interfaces (D157)](#9-fluidsolid-interfaces-d157)
   - 9a. [Thin walls: wall conduction without meshing the wall (D189)](#9a-thin-walls-wall-conduction-without-meshing-the-wall-d189)
   - 9b. [Wall film heat exchange with walls, thin walls and solids (D190)](#9b-wall-film-heat-exchange-with-walls-thin-walls-and-solids-d190)
   - 9c. [Several fluids and materials: heat exchangers (D194)](#9c-several-fluids-and-materials-heat-exchangers-d194)
   - 9d. [Cycle-averaged coupling of periodic flows (D204)](#9d-cycle-averaged-coupling-of-periodic-flows-d204)
   - 9e. [Thermal acceleration of solids (D213)](#9e-thermal-acceleration-of-solids-d213)
10. [Coupling strategy, linear solvers and order of solution (D158)](#10-coupling-strategy-linear-solvers-and-order-of-solution-d158)
11. [Multi-region partitioning and halos (D159)](#11-multi-region-partitioning-and-halos-d159)
12. [Kernels: GPU, reproducibility, adjoint (D160)](#12-kernels-gpu-reproducibility-adjoint-d160)
13. [Energy ledger and second law (D161)](#13-energy-ledger-and-second-law-d161)
14. [Verification and validation (D162)](#14-verification-and-validation-d162)
15. [Case-file keys and diagnostics (D163)](#15-case-file-keys-and-diagnostics-d163)
16. [Implementation plan (D164)](#16-implementation-plan-d164)
17. [Findings in the code today](#17-findings-in-the-code-today)
18. [Points for the maintainer](#18-points-for-the-maintainer)
19. [References](#19-references)


## 1. Place in the plan

Until now the energy equation and conjugate heat transfer were L1, the first
item after milestone 1, and stage 7 was multi-rank MPI with overlapped halos,
scaling, the internal-flow validation (pipe, diffuser, U-bend, T-junction),
rotating and non-conformal meshes, the reproducible mode and the tag `v1.0`.
The maintainer's prompt of 2026-11-16 makes energy, buoyancy and conjugate
heat transfer stage 7. New order (O18, confirmed by the maintainer on 2026-11-17):

| Stage | Content | Tag |
|---|---|---|
| 6 | Turbulence and wall functions, including the thermal wall function and the turbulent Prandtl number with temperature as a passive scalar ([stage 6 page](stage6-turbulence.md), D129–D139) | v0.6 |
| **7** | **Energy, buoyancy, conjugate heat transfer (this page, D150–D164)**, with the static non-conformal interfaces of D58 moved here from the old stage 7, because fluid and solid meshes made separately (D113) seldom match | v0.7 |
| 8 | The old stage 7 without D58: multi-rank MPI with overlapped halos (D38), scaling, the internal-flow cases of D45 (now also with heat transfer: heated periodic pipe D111 with every turbulence model, T-junction with hot and cold inflow), cell-quality remediation (D50), tet-to-polyhedra conversion (D54), conformal rotating meshes (D56, moved to stage 6 with D182 on 2026-11-19), the reproducible mode (D93) and its cost, GCI in the validation reports, user documentation | **v1.0** |

What this means for milestone 1: v1.0 moves one stage later, and its scope
grows from "incompressible flow" to "incompressible and low-Mach
variable-density flow with heat transfer, buoyancy and conjugate heat
transfer". Species and reacting flows (L2) stay after milestone 1, and so do
electric conduction and Joule heating (L27), which L1 had taken along: they
are the first item after v1.0, because they need nothing beyond the
temperature domain of this stage. The order after milestone 1 becomes L27,
L2, L12, L3, L4, L5.

Alternatives considered: (a) keep v1.0 at the end of the old stage 7 and do
energy as stage 8 (v1.1): a first release without heat transfer, which the
industrial reference class always has; (b) merge both into one stage 7: too
large for one verifiable step.


## 2. What stage 7 builds on

Agreed decisions that this stage implements or extends, and does not redo:

- **D112** regions and equation domains: one mesh, the temperature solved on
  the union of fluid and solid regions, so conjugate heat transfer is one
  matrix. D99: harmonic mean of k across fluid–solid interfaces. D116:
  generic transport equation with tensor diffusion, thermal boundary
  conditions as an `energy` sub-table of `[boundary.<name>]`, periodic pairs
  and region interfaces in one paired-boundary structure. D117 (1): region to
  region coupling computed once and booked on both sides. D118: modules,
  step schedule, multi-rate clocks, frozen flow. D119: ledger as data.
- **D92** energy as specific sensible enthalpy h_s, c_p in integral form,
  temperature recovered by a safeguarded Newton inversion; **D88–D90** the
  conservation ledger, the balance reports and the second-law checks, whose
  parts 1–3 become active here.
- **D94–D103** material database, property forms and tables, the property
  kernel, `[materials]` keys; D102 puts the thermal part (equations of state,
  k, c_p, solids with orthotropic k, Sutherland and NASA data, IAPWS water)
  with L1, so into this stage.
- **D109–D111** streamwise-periodic heat transfer (affine map on periodic
  pairs, thermal modes, Nu and f evaluation, the heated periodic pipe as the
  main validation case).
- **D13** and the stage-6 decisions (D129 exported turbulence fields, D134
  momentum wall treatment, D135 thermal wall treatment with Kader 1981 as the
  default and the Jayatilleke P-function as an option, and the turbulent
  Prandtl number, constant 0.85 or the variable form of Kays 1994): stage 7
  cites them and does not specify them again. Stage 6 leaves a source hook
  `buoyancy_production` in the k, ε and ω equations that stage 7 fills (D155).
- **L15** radiation: the photon Monte Carlo module comes later; this stage
  gives the radiation boundary condition to a surrounding temperature
  (D156) and keeps the interface the radiation module will use (a volume
  source and a wall heat flux per face, both with their linearisation).
- **D124** W3, W4, W7: wall boundary values as per-face arrays, per-face
  heat-transfer coefficients from the wall treatment, gravity as a case
  vector for all modules.
- **D58** non-conformal interfaces by imprinting: after imprinting the
  interface is conformal, so the fluid–solid treatment below needs no special
  case for non-matching meshes.


## 3. Energy equation (D150)

### 3.1 Equation

The conserved quantity is the specific sensible enthalpy h = h_s(T) of D92.
On the energy domain (fluid and solid regions):

    ∂(ρh)/∂t + ∇·(ρ u h) = ∇·(k ∇T) − ∇·q_t + Φ + W_p + S_h

| Term | Meaning | Where |
|---|---|---|
| ∂(ρh)/∂t | storage, BDF2 or Crank–Nicolson (D28) | fluid and solid |
| ∇·(ρ u h) | convection with the mass flux of the flow | fluid (u = 0 in solids; moving solids with D122 (3) later) |
| ∇·(k ∇T) | Fourier conduction, k scalar or tensor (D116) | fluid and solid |
| q_t | turbulent heat flux, SGDH or GGDH (D152) | turbulent fluid |
| Φ | viscous heating, off by default: laminar τ:∇u with τ = μ(∇u + ∇uᵀ − ⅔(∇·u) I); in RANS the mean part τ̄:∇ū plus the dissipation ρε of the turbulent kinetic energy (ε from the model of D129; for SST ε = β* k ω) | fluid |
| W_p | pressure work, `none` (default), `thermodynamic` (dp₀/dt of a closed volume, D151, L22) or `full` (Dp/Dt = ∂p/∂t + u·∇p) | fluid |
| S_h | volume heat sources: number, expression (D53) or total power over a region (W, spread by volume), with an optional linear part S_h = S_c + S_T T (S_T ≤ 0) | any region |

For a constant-density liquid h = ∫c_p dT; the pressure part of the enthalpy
of a liquid is neglected (the usual incompressible form, Bejan 2013 ch. 1).
Radiation (L15) enters later through S_h and the wall heat flux.

### 3.2 Solved unknown: temperature correction, conservative residual

D92 (1) keeps h as the conserved variable and recovers T by inversion.
Stage 7 keeps the conservative residual and changes the unknown of the
linear system: each outer iteration solves for the temperature correction
δT, with the residual written in h:

    R(T) = conservative finite-volume residual of the h equation above
    A δT = −R(T*),   T = T* + δT,   h_P = h(T_P) from the D92 coefficients

The matrix A is the linearisation of R about T*: storage ρ_P c_p(T*_P) V_P a₀/Δt
(a₀ the BDF coefficient), convection with the upwind cell's c_p (first-order
upwind implicit part, the higher-order face value as deferred correction,
D25), conduction directly in T. At convergence R = 0, so energy is conserved
exactly with temperature-dependent c_p and the ledger of D88 closes; c_p
appears only in A, which changes the speed of convergence, never the result
(the source-based enthalpy linearisation of Voller & Swaminathan 1991, which
also carries the latent heat of L35 later).

Why the change: (1) temperature, not enthalpy, is continuous across a
fluid–solid interface, so one matrix across materials with different c_p is
only well posed in T; (2) Robin, radiation and wall-function conditions are
linear or linearised in T; (3) the conduction-dominated solid rows are a
symmetric positive M-matrix in T, which the AMG handles at conductivity
ratios of 10⁴ (§10.2); (4) the deferred correction of D92 (3) between
(k/c_p)∇h and k∇T is no longer needed, because conduction is implicit in T.
D92 (3) is replaced by this; D92 (1), (2), (4), (5) hold unchanged. The
Newton inversion of D92 (1) stays for every place where h is given and T is
wanted: inlet enthalpy, the affine periodic map of D109 (which is written in
h), restart files, mixing (L2).

Periodic pairs (D109) with the map h_partner = a h + b: in the linear
system the partner's δT is a c_p(T)/c_p(T_partner) δT, the residual uses the
mapped h exactly.

### 3.3 Discretisation

- Convection of h: the face value from the cell values of h and the limited
  least-squares gradient of h (D24, D25), linear upwind by default, bounded
  central as an option; boundedness of T between its inflow and wall extremes
  is checked as for the passive scalars (D51).
- Continuity correction of the implicit part: the implicit operator
  subtracts c_p,P δT_P times the cell's discrete mass imbalance
  (∂ρ/∂t + Σ ṁ_f) so that the matrix stays diagonally dominant while the
  continuity equation is not yet converged (Patankar 1980, §5.2); the
  residual keeps the conservative form, so the converged result is
  unchanged.
- Conduction: over-relaxed non-orthogonal decomposition (D26) with k_f the
  arithmetic mean inside a region (D99) and the interface conductance of
  §9 across regions; tensor k by Γ·S in place of Γ S (D116).
- Steady runs: implicit under-relaxation per row (Patankar 1980), default
  α_T = 0.95 in fluid rows and 1 in solid rows, or the pseudo-transient
  local time step of D49 with a diffusion-number time step in solid rows,
  because a solid has no convective Courant number.
- Positivity: a correction that would take T to ≤ 0 K or out of the
  property range is damped by a factor applied to the whole update
  (a scalar from an exact reduction, so the result stays reproducible) and
  counted as a health number (D69).

### 3.4 Variable properties

Properties are cell fields evaluated by the property kernel of D99 after
the temperature update, once per outer iteration: ρ(T, p₀) (D151), μ(T),
k(T) (scalar or orthotropic tensor with axes per region), c_p(T) with h
and s° from the same coefficients (D92 (2)). Forms from D95: constant,
polynomial, NASA-7 and NASA-9, Sutherland and power law for μ and k of
gases, tables (linear for c_p so that h is exact, monotone cubic for
transport properties, D96), expressions (D53, also of T now: the
restriction to x, y, z, t of stage 5 is lifted for properties). The thermal
part of the database (D102): air, water (IAPWS-95 for the tables of
c_p, k, μ, ρ at the operating pressure, generated off line under D103),
the oils and glycols of stage 5 with thermal data, and solids: copper,
aluminium, steel, stainless steel, glass, a polymer, with sources (D97).


## 4. Low-Mach variable density (D151)

Without species, D17's low-Mach form arrives here; L2 adds species to it.

- **Density:** ρ = ρ(T, p₀) from the material's equation of state (ideal gas
  ρ = p₀ M/(R T), polynomial or table in T for liquids); it does not depend
  on the dynamic pressure p (the low-Mach expansion, Rehm & Baum 1978;
  Paolucci 1982), so the pressure equation stays elliptic and the pressure
  matrix symmetric, solved by AMG-CG as today (D35).
- **Thermodynamic pressure p₀:** open domains (any pressure boundary): p₀ =
  `fluid.operating_pressure`, constant. Closed domains: p₀(t) from the
  conservation of the total mass, M_total = ∫ρ(T, p₀) dV, solved for p₀
  after each temperature update (for an ideal gas p₀ = M_total R / (M ∫dV/T),
  an exact reduction, D31); dp₀/dt then enters W_p of the energy equation
  and the mass balance. This is the p₀(t) of L22, without mesh motion.
- **Continuity:** ∂ρ/∂t + ∇·(ρu) = 0 with ∂ρ/∂t from the BDF formula and the
  latest ρ; it is a known source of the pressure-correction equation in each
  outer iteration.
- **Face density:** ρ_f = w ρ_P + (1 − w) ρ_N in the mass flux and in the
  Rhie–Chow coefficient, the same ρ_f in the body force of D154.
- **Momentum:** ∂(ρu)/∂t + ∇·(ṁ u) = −∇p + ∇·τ + (ρ − ρ_ref) g, with μ(T) a
  cell field.
- **Validity check:** the Mach number of the cell velocities against the
  speed of sound of the material is a health number; above Ma 0.3 a warning
  (ENE004) says that the low-Mach form loses validity (L10).

**As built (PR 7c, #195, 2026-11-22).** The points that differ from the
text above or make it concrete:

1. The switch is opt-in: `[energy] density = "variable"`, or `[buoyancy]
   model = "variable-density"`, which turns it on. A fluid whose material
   has an ideal-gas form keeps the constant density of `[fluid]
   temperature` without it, so every case before 7c stays bitwise
   unchanged. The density form of the fluid gives ρ_op(T) at the
   operating pressure; the cells take ρ = ρ_op(T) p₀/p_op for an ideal gas
   and ρ_op(T) for a liquid. A closed volume of a liquid is refused
   (ENE009): it has no low-Mach solution.
2. p₀ of a closed volume is the exact reduction p₀ = M p_op / Σ ρ_op(T) V
   after every temperature update, with M the mass of the start state, so
   the mass stays to round-off; the mass source of the pressure equation is
   V(a₀ρ − a₁ρⁿ − a₂ρⁿ⁻¹)/Δt with the BDF coefficients, and the momentum
   time term takes the density of each time level. The pressure work is
   W_p = dp₀/dt (`pressure_work`, default `thermodynamic` in a closed
   volume, `none` in an open one; the `full` form stays with L22).
3. The energy equation keeps its bounded, non-conservative form with ρ per
   cell, so its storage needs no correction of the density levels; the heat
   ledger books the storage as V(ρ ∂h/∂t + h ∂ρ/∂t). The momentum equation
   takes the full stress with −⅔ μ ∇·u (the kernel's trace term with the
   factor ⅔). The kinetic-energy ledger books the pressure dilatation
   Σ p_P Σ_f ṁ_f/ρ_f as a physical term, and the mass ledger the storage of
   the density.
4. Checks: V7c.1 a heated duct against its exact one-dimensional balances,
   V7c.2 a closed box heated uniformly (T and p₀ linear in t, exact at any
   Δt), V7c.3 a closed box heated on one side (mass to round-off), V7c.4 a
   manufactured solution with ρ = 1/T in a closed volume (second order on
   hexahedra and polygons, SIMPLEC and the coupled algorithm on the same
   solution). V7.17 compares with no published numbers: its references
   (Vierendeels et al. 2003, Le Quéré et al. 2005) are results of other
   codes and fall under D201; the large-ΔT cavity of gallery 195 carries
   our own GCI instead (heat flow 0.020 %, p₀ 0.003 % on 32² to 128²
   clustered meshes).
5. Not yet with the variable density (ENE006): passive scalars, a restart,
   periodic meshes, convective outlets, mass-flow inlets, wall mass fluxes,
   rotating zones. The turbulence outlet value ⅔ρk uses the reference
   density.


## 5. Turbulent heat flux (D152)

Stage 6 (D135) owns the turbulent Prandtl number (constant 0.85 default,
Kays 1994 variable form as an option) and the thermal wall function. Since
step 21 (2026-11-17) Pr_t is a field from a closure (D172), read by the
conductivity below, by G_b (D155) and by the wall function; stage 7 adds the
`kays-crawford-hp` (default parameter set `2019`, Steiner & Irrenfried 2019)
and `wall-profile` closures and the variable-property wall treatment of D173
([step 21 page](step21-les-high-prandtl.md) §6). Stage 7
owns the closure of the turbulent heat flux in the energy equation:

- **SGDH** (simple gradient diffusion): q_t = −(μ_t c_p / Pr_t) ∇T, added to
  the conductivity of the conduction term, k_eff = k + μ_t c_p / Pr_t; with
  a variable c_p the c_p is the cell's, so the term stays a flux in T.
- **GGDH** (generalised gradient diffusion, Daly & Harlow 1970):
  ρ \overline{u′_i h′} = −C_θ ρ c_p T_t \overline{u′_i u′_j} ∂T/∂x_j, i.e. a
  tensor conductivity K_ij = C_θ ρ c_p T_t \overline{u′_i u′_j}. Inputs from
  D129: the time scale T_t (k/ε bounded below by the Kolmogorov scale, as the
  stage-6 models bound it) and the cell field `reynolds_stress` (six
  components): for ζ-f and the elliptic-blending model the stress with the
  wall-normal variance v² = ζk (or φk) along the wall-normal direction
  n = ∇α/|∇α|, for SST and the k-ε models the Boussinesq estimate
  ⅔kδ_ij − 2ν_t S_ij. Default C_θ = 0.3 (Ince & Launder 1989, GGDH in
  enclosures), a model constant in the case file; the value and its
  pairing with each model are checked against the papers before the code
  is written.
- **Implicit part of GGDH:** the face-normal component of K (n·K·n |S|) is
  implicit when positive and its negative part and the cross terms (K·S
  minus the normal part) are deferred corrections (D26 style); because the
  Boussinesq stress can be non-realisable, the implicit coefficient is made
  non-negative with the smooth positive part of D160, so the matrix stays an
  M-matrix. On walls the wall function of D135 replaces the turbulent flux
  of the first cell.
- **Later (P12):** the algebraic flux model (Hanjalić, Kenjereš & Durst 1996;
  Kenjereš, Gunarjo & Hanjalić 2005) with the temperature variance, for
  strongly buoyant and mixed convection; not in stage 7.


## 6. Buoyancy (D153) and the balanced-force interpolation (D154)

### 6.1 Gravity and the two buoyancy models (D153)

Gravity is one case vector `[gravity] vector` with a reference point x_ref
(D124 W7), seen by every module. The solved pressure is the reduced
pressure

    p′ = p − ρ_ref g·(x − x_ref)

so the hydrostatic part of a fluid at the reference density is removed
exactly and the body force is (ρ − ρ_ref) g. Two models:

- **Boussinesq** (`buoyancy.model = "boussinesq"`): ρ = ρ₀ constant
  everywhere except in the body force, (ρ − ρ_ref) g = −ρ₀ β (T − T_ref) g,
  with β from the case or the material (ideal gas: β = 1/T_ref). A
  check warns when β (T_max − T_min) exceeds 0.1 or the property change
  over the expected temperature range exceeds the limits of Gray &
  Giorgini (1976) (BUO002).
- **Variable density** (`"variable-density"`): ρ(T, p₀) of D151 everywhere,
  body force (ρ − ρ_ref) g with ρ_ref = ρ(T_ref, p₀) by default.

The form p − ρ g·x with a variable ρ (which gives a force −(g·x)∇ρ) is not
used: with ρ_ref constant the reduced pressure is exact for the reference
state and needs no density gradient, which is noisy at sharp density
changes.

Boundary values of pressure boundaries are values of p′ when gravity is
set: `pressure = 0` at an opening to an ambient at rest at ρ_ref is then
the correct hydrostatic outlet. A boundary may instead give the static
pressure (`pressure_kind = "static"`); it is converted per face by
p′_b = p_b − ρ_ref g·(x_f − x_ref). The output writes p′ (`pressure`) and the
static pressure p (`static_pressure`, a derived field of D72). Total
pressure boundaries use p′ + ½ρ|u|².

### 6.2 Balanced-force face interpolation (D154)

With a large body force the present scheme has spurious velocities near
walls (documented in docs/verification/stage5.md): the cell rows use
(b_P − ∇p_P) V with the least-squares ∇p, which is not exact for a
hydrostatic pressure that is quadratic when ρ varies, and the Rhie–Chow
face velocity sees only the pressure. The cure is to treat the body force
exactly as the pressure gradient, from the same face differences
(Gu 1991; Choi et al. 2003; Mencinger & Žun 2007; Zhang, Zhao & Bayyuk 2014;
the same principle as the balanced-force surface tension of Francois et al.
2006, which VOF (L3) will reuse):

1. **Reduced face jump.** For each face, with d = x_N − x_P and the face body
   force b_f = (ρ_f − ρ_ref) g (ρ_f as in D151, for Boussinesq from the
   interpolated T):

       J_f = (p′_N − p′_P) − b_f·d

   On a boundary face J_b = (p′_b − p′_P) − b_P·(x_b − x_P).
2. **Cell net force.** The cell term −∇p + b of the momentum rows is
   replaced by −G̃_P, where G̃_P is the least-squares gradient operator of
   D24, with its weights unchanged, applied to the jumps J_f instead of the
   pressure differences. If the discrete pressure satisfies J_f = 0 on every
   face (discrete hydrostatics), G̃_P = 0 exactly on every cell, on every mesh.
   The pressure gradient itself (∇p′ for output, gradients for the limiter
   and the remediation of D50) is unchanged.
3. **Rhie–Chow.** In `face_mass_flux` the term (jump − ∇p_f·d) becomes
   (J_f − G̃_f·d), with G̃_f interpolated from the cells as ∇p_f is today;
   the time-step-independent form (D27) is unchanged. In discrete
   hydrostatics the face flux is then exactly zero.
4. **Wall pressure.** The zero-gradient wall pressure becomes zero gradient
   of the reduced jump: p′_b = p′_P + b_P·(x_b − x_P) + G̃_P·(x_b − x_P), so
   that a wall carries the hydrostatic pressure without a jump. The same
   rule for symmetry, slip and wedge faces.
5. **Pressure correction.** The correction equation is unchanged (the body
   force is explicit in each outer iteration); the face body force is
   updated with T and ρ in every outer iteration and corrector.

The scheme keeps second order (G̃_P = ∇p − b to the accuracy of the
least-squares gradient on smooth fields) and reduces to the present one
without gravity, bitwise (the extra term is zero when g = 0, and the code
path is not entered).

**As built (PR 7b1, #186, 2026-11-22).** Four points differ from the
text above or make it concrete:

1. The face body force of the offsets is the arithmetic mean,
   o_f = ½(b_P + b_N)·d, not the interpolation with the face weight. For a
   linear stratification the hydrostatic p′ is quadratic along d, and the
   mean of the two end slopes is its exact difference on any face position,
   so the discrete hydrostatic state is exact on skewed and off-centre
   faces as well (the weighted form is not). The Rhie–Chow term uses
   o_f − B̄_f·d with B̄_f weighted as ∇p_f; the boundary faces take
   o_b = b_P·r, the faces with a zero-gradient pressure are left out of
   B_P as they are out of ∇p.
2. T_ref defaults to `[fluid] temperature` (the temperature of ρ₀), and
   ρ_ref to ρ₀; `[buoyancy] reference_temperature` and
   `reference_density` override them.
3. The steady algorithms take a buoyancy pseudo-time step: SIMPLEC and
   SIMPLE add ρV/Δτ with Δτ = L/U_b = 1/N (N = √(|g| β ΔT / L)), with the
   current iterate as the old level and the time term of the Rhie–Chow
   flux (D27); the automatic step of the coupled algorithm (D177) takes
   U_b as a lower bound of its speed. Without it the relaxation alone takes
   steps of about α/(1 − α) ρ h²/μ in a fluid at rest, N Δτ ≫ 1 on coarse
   meshes, and the explicit force overshoots the internal waves of a stable
   stratification without end (a 5 × 5 box kept 2·10⁻² m/s). With it the box
   is at rest to 10⁻¹⁵ m/s on every mesh, and the heated cavity at Ra 10⁴
   converges in 235 iterations instead of 202, to the same Nu.
4. BUO004 warns for `[gravity]` with the energy equation and no buoyancy
   model, without the Richardson estimate; BUO005 came with 7c (#195). The wall
   forces of the report are from p′ (the Archimedes part ρ_ref g V is not in
   them).

### 6.3 Body force with a density step

For a sharp density change (a hot layer over cold fluid, two regions of
different fluids, later VOF), b_f from the interpolated ρ_f gives a discrete
hydrostatic state with J_f = 0 on every face, so a fluid at rest stays at
rest to round-off (verified in V7.7, §14). Mencinger & Žun (2007) show that
an interpolation of the cell forces instead of the face form produces the
spurious currents this avoids.


## 7. Buoyancy in the turbulence models (D155)

Stage 6 leaves the hook `buoyancy_production` in the k, ε and ω equations
with an explicit part S_u ≥ 0 and an implicit part S_p ≤ 0 per cell (D129).
Stage 7 fills it:

- **Production of k by buoyancy:** G_b = −ρ β g_i \overline{u′_i T′} for
  Boussinesq and G_b = g_i \overline{ρ′u′_i} for variable density. With
  SGDH: G_b = β (μ_t / Pr_t) g·∇T (Boussinesq) and
  G_b = −(μ_t / (ρ Pr_t)) g·∇ρ (variable density); both are negative in a
  stable stratification (destruction of k).
  With GGDH: G_b = ρ β C_θ T_t g_i \overline{u′_i u′_j} ∂T/∂x_j; GGDH is required
  for a vertical boundary layer, where g·∇T ≈ 0 and SGDH gives almost no
  buoyant production (Ince & Launder 1989; Hanjalić 1994 for the reasoning).
  The choice follows `energy.turbulent_heat_flux` unless
  `buoyancy.production` names another form.
- **Split into the hook:** S_u = P⁺(G_b), S_p = −P⁺(−G_b)/k, with P⁺ the
  smooth positive part of D160, so that a destructive G_b is implicit and
  k stays positive (D129 positivity), with no branch per cell.
- **ε equation:** C_ε1 C_ε3 (ε/k) G_b with C_ε3 = tanh(|v_∥| / |u_⊥|)
  (Henkes, van der Vlugt & Hoogendoorn 1991): v_∥ the velocity component
  along g, u_⊥ the component normal to it; written with the smooth absolute
  value and a small velocity floor so that it is C¹ (D59). C_ε3 = 1 for a
  vertical plume, 0 in a horizontal shear layer. Option: a constant C_ε3.
- **ω equation (SST):** from the ε form through ε = β* k ω (the
  transformation of Wilcox 2006, ch. 4), with the SST coefficient α in place
  of C_ε1 − 1: S_ω,b = (ω/k) G_b ((α + 1) C_ε3 − 1), split into S_u and S_p as
  above. Published choices differ (some drop the C_ε3 part); the case file
  can select `buoyancy.omega_source = "transformed" | "none"`, and the
  vertical plate and the turbulent cavity (V7.13, V7.14) decide the default.
- **ζ-f and elliptic blending:** the production in the ζ (or φ) equation
  uses P_k + G_b in place of P_k (Kenjereš, Gunarjo & Hanjalić 2005).
- **Constants** are case settings with their sources in the schema (D77).


## 8. Thermal boundary conditions and their check (D156)

Each boundary takes a `[boundary.<name>.energy]` table (D116), in the
style of the `[boundary.<name>.scalar.<scalar>]` tables of the passive
scalars; the flow `type` is unchanged. The condition follows from the keys
given, so there is no `type` key: exactly one condition per table (THM010),
and a key the boundary kind cannot take is refused (THM003). Values are
numbers or expressions of x, y, z, t and, for Robin and radiation
conditions, T (D53). *Decided by the maintainer on 2026-11-21 (PR 7a1,
#154); the planned `thermal = { type = ... }` form is dropped.*

| Boundary | Condition | Keys of `[boundary.<name>.energy]` | Discrete form |
|---|---|---|---|
| wall, slip wall | adiabatic (default, no table; THM005 lists these walls once) | – | zero flux |
| wall, slip wall | temperature | `temperature` (K) | Dirichlet on the face; with a wall function the conductance is h_w of D135 |
| wall, slip wall | heat flux | `heat_flux` (W/m², positive into the domain) | Neumann; the face temperature is reconstructed for output |
| wall, slip wall | heat flow | `heat_flow` (W over the boundary) | heat flux = heat flow / area |
| wall, slip wall | convection | `heat_transfer_coefficient`, `ambient_temperature` | Robin: q = U (T_amb − T_w), U the series of h_ext, the optional thin wall and the near-wall conductance; implicit in T |
| wall, slip wall | radiation (7a3) | `emissivity`, `ambient_temperature` | q = εσ(T_amb⁴ − T_w⁴), Newton-linearised about the current T_w: h_rad = 4εσT_w*³, implicit; converges to the quartic exactly |
| wall, slip wall | convection and radiation (7a3) | `heat_transfer_coefficient`, `emissivity`, `ambient_temperature` | the sum of both |
| wall (any) | optional thin wall (7t1) | `[boundary.<name>.energy.thin_wall]`: `thickness` and `conductivity` or `material`; or `layers`, `conduction`, `edges` (§9a) | simplest form: a 1D resistance t/k in series, no heat capacity; the full thin wall of §9a (D189) adds layers, heat capacity, sources and in-plane shell conduction |
| inlets | temperature (required, THM011) | `temperature` (K) | Dirichlet h(T_in) on inflow |
| outlets | zero gradient; backflow (7a3) | `backflow_temperature` (K) | used on faces with reverse flow (D29 backflow) |
| symmetry, wedge | zero normal flux | – (any key: THM003) | zero normal flux |
| periodic | T carried across; heat-flux mode (7a2) | `bulk_temperature` of D109 | affine map in h |
| region interface | (no boundary) | `[interfaces.<name>]` with `contact_resistance` or `thin_wall` (§9) | conductance of §9 |

An inlet needs its temperature: a silent default (the initial temperature,
as the first draft had it) hides a forgotten inlet value, and the inlet
temperature sets the result of every heated flow. The start temperature is
`[energy] initial`.

Every wall value (temperature, heat flux, h_ext, T_amb) is a per-face array
(D124 W3) that another module can set each step (radiation L15, film D123,
external system models D105); the same holds on thin-wall faces and
region-interface faces, with the linearisation of the flux (§9b, D190).

**Check of the thermal set-up (D30 extended, codes `THM`):** run before the
solver, with the guidance of D30 (what to change):

- THM001 error: steady run whose boundaries and sources leave the
  temperature level undetermined (only adiabatic and heat-flux walls, no
  inflow, no temperature, Robin or radiation boundary); "add a temperature
  or convection boundary, an inflow, or run transient".
- THM002 error: steady run with a net prescribed heat input and no path for
  heat to leave (as THM001 with non-zero total): the temperature drifts.
- THM003 error: a thermal value on a boundary that cannot carry it
  (`heat_flux` on an inlet, `temperature` on a symmetry plane).
- THM004 error: emissivity outside (0, 1], heat-transfer coefficient ≤ 0,
  temperature ≤ 0 K.
- THM005 warning: energy on, a wall without a `[boundary.<name>.energy]`
  table (adiabatic assumed), listed once.
- THM006 warning: inlet or initial temperature outside the validity range of
  a material (with MAT codes of D98).
- THM007 error: the D109 classification (wall-temperature mode with mixed
  thermal walls, heat-flux mode with isothermal walls), as specified there.
- THM008 warning: viscous heating off although the Brinkman number estimate
  μU²/(kΔT) exceeds 0.1, or on although it is below 10⁻⁴.
- THM009 error: a `[boundary.<name>.energy]` table without `[energy]` (the
  energy equation is on when the case has an `[energy]` table; there is no
  `enabled` key).
- THM010 error (PR 7a1): a `[boundary.<name>.energy]` table of an unknown
  boundary, without exactly one condition, or with
  `heat_transfer_coefficient` and `ambient_temperature` apart.
- THM011 error (PR 7a1): an inlet without a temperature.
- Regions: REG001 (non-fluid region refused) stays for flow-only runs and is
  lifted for runs with the energy equation; MAT006–MAT009: a solid without
  conductivity, specific heat or density, a fluid without thermal data when
  energy is on, an orthotropic conductivity without axes.

Buoyancy checks (codes `BUO`): BUO001 error, buoyancy without gravity or
gravity of zero length; BUO002 warning, Boussinesq outside its range (§6.1);
BUO003 error, Boussinesq with a material whose density form is constant and
no β given; BUO004 warning, gravity set and buoyancy `none` with a
temperature range that would give Richardson numbers above 0.1; BUO005
error, variable density in a closed domain without `energy.thermodynamic_pressure
= "closed-volume"` (the mass cannot be conserved, like BC026);
BUO006 warning, a pressure boundary whose height span gives
ρ_ref g Δz above the dynamic pressure, with `pressure_kind` explained.


## 9. Fluid–solid interfaces (D157)

Fluid and solid cells of one energy domain meet at interior faces that the
domain view of D112 classifies as region interfaces. Conformal interfaces
come from one mesh or from parts joined with matching nodes (D113);
non-conformal ones are imprinted by D58 into ordinary faces first. The
interface treatment below therefore sees only conformal faces.

### 9.1 Conductance and face temperature

For an interface face f between cell P on side a and cell N on side b, with
unit normal n̂ (from P to N), normal distances δ_P = (x_f − x_P)·n̂ and
δ_N = (x_N − x_f)·n̂, and the face-normal temperatures

    T̃_P = T_P + ∇T_P·(r_P − δ_P n̂),   r_P = x_f − x_P   (likewise T̃_N)

(the tangential offsets of the non-orthogonal correction, deferred):

    q_f = G_f (T̃_P − T̃_N),   G_f = |S| / (R_P + R″_c + R_N)
    R_P = δ_P / k_P (conduction side) or 1/h_w,P (fluid side with a wall function, D135)

R″_c is the contact resistance (m²K/W) or t/k of a thin wall, zero by default.
For R″_c = 0 and conduction on both sides this is the harmonic mean of D99
(Patankar 1980, §4.2), written with the normal distances to the face rather
than the line weight, so it stays exact when the face is not halfway. The two
face temperatures follow from the series: T_f,a = T̃_P − q_f R_P/|S|,
T_f,b = T̃_N + q_f R_N/|S| (equal when R″_c = 0). With orthotropic k the
normal conductivity n̂·k·n̂ is used in R and the cross terms are deferred
(D116).

**Orientation of an orthotropic k (L48).** The material gives k₁, k₂, k₃ and,
as the default, the axes `axis_1`, `axis_2` (the second one is made normal to
the first, the third is their cross product); a region (D112) may give its own
axes, so one winding material can be oriented differently in every slot. The
tensor K = R diag(kᵢ) Rᵀ is a symmetric six-number cell field. Stage 7 takes
fixed axes; cylindrical frames (laminated stacks and windings of e-motors),
frames from expressions and from per-cell fields, and the rotation of a frame
with a rotating zone follow with L30 (L48, `requirements.md`). A test of
stage 7: the linear-field patch test with a rotated K on a skewed mesh, exact to
round-off (V7.4 variant).

The flux is one number per face, computed once and added to P and
subtracted from N (D117 (1)), so heat flux continuity
−k_a∇T_a·n = −k_b∇T_b·n holds to round-off, at every iteration, not only at
convergence.

### 9.2 One-sided gradients

∇T is discontinuous across an interface. The least-squares gradient of T in
a cell next to an interface uses the face temperature T_f of its own side
at x_f in place of the neighbour cell across the interface (a boundary-type
point, as at Dirichlet boundaries, D24), so each region's gradient sees only
its own material. The same for ∇h used by the convection limiter. With
piecewise-linear T in each region the gradients, the face temperatures and
the fluxes are then exact on any mesh, at any conductivity ratio (V7.4).
Because T_f depends on the gradients, the pair is updated once per outer
iteration (deferred); in a steady run this costs no extra iterations at
convergence, in a transient run it is converged by the outer correctors.

### 9.3 Large conductivity ratios

At k_s/k_f ≈ 10⁴ (copper and air) the interface conductance is dominated by
the fluid side and the solid side is nearly isothermal. The scheme above
loses no accuracy there (R_N ≪ R_P is just a small term in a sum). The
linear solver is the concern, addressed in §10.2. The interface heat flow
and the interface temperature are reported per interface in the log and the
run record.

### 9.4 Turbulent fluid side

With a wall function on the fluid side, R_P = 1/h_w,P with
h_w = ρ c_p u_τ / T⁺(y⁺, Pr) from D135 (Kader 1981 by default); the solid
side keeps δ_N/k_N. With the compound wall treatment (D13) the blending is
inside h_w. y⁺ and T⁺ per interface face are derived fields (D72). Built
in PR 7d5 (#158, V7.36): the conductance of an interface face is
1/(1/(h_w d_n) + w_s/k_s), with h_w = Γ_w/y_P, Γ_w = μ c_p y*/T⁺(ỹ*) from
y* and μ of the fluid view's first cell, taken in every outer iteration
(ỹ* = max(y*, 11.06) with the scalable treatment, as at the fluid's other
walls); inside the fluid the conductivity is k + c_p μ_t/Pr_t. For periodic flows with a slow structure (engine exhaust
manifolds) the cycle-averaged coupling of §9d (D204) uses this conductance.

### 9.5 Ledger

Each interface is an exchange pair of the ledger (D119) between the two
regions' energy balances: one face sum, two entries of opposite sign, which
cancel to the bit. The entropy produced by a contact resistance,
q²R″_c/(T_f,a T_f,b) per area, is a part of the second-law check (D90 (1)).


## 9a. Thin walls: wall conduction without meshing the wall (D189)

Added 2026-11-20 at the maintainer's request. A thin wall is a simplified
conjugate heat transfer: a wall, a plate or a sheet whose thickness is not
meshed. It lives on boundary faces (or on a face pair, below) and carries
heat through its thickness (1D) and, if asked, along its surface (2D, shell
conduction). It replaces a solid region where the wall is thin against its
extent, the usual case for sheet-metal ducts, housings, heat-exchanger
plates, glass panes, printed circuit boards and tank walls, and where meshing
the wall would cost more cells than the fluid.

The `thin_wall` of §8 (a series resistance t/k without heat capacity) is the
simplest case of this model and stays valid with the same meaning.

### 9a.1 Where a thin wall can sit

| Placement | Sides | Typical use |
|---|---|---|
| **Wall boundary** | fluid cell on one side, the boundary's `energy` condition (adiabatic, temperature, heat flux, convection, radiation, §8) on the outer side | a duct or a housing whose outside is ambient |
| **Baffle** (a pair of coincident boundaries, D116 paired-boundary structure) | a fluid cell on each side | a plate inside the flow: a fin, a separating plate of a plate heat exchanger, a shield |
| **Region interface** (§9) | fluid and solid, or two solids | a coating, a contact layer, a gasket with heat capacity |
| **Inside a solid** (a baffle in a solid region) | the same solid, or two solids, on both sides | a glue or solder layer, a crack, an air gap or a thermal pad inside a part |
| **Between two fluid regions** (§9c) | two different fluids | the tube wall of a tube-and-fin exchanger, the plate of a plate exchanger, without meshing the wall |

Thin walls are allowed on every face set of the energy domain, whatever its
two sides are (maintainer, 2026-11-20): boundaries, inside fluids, inside
solids, at fluid–solid, solid–solid and fluid–fluid interfaces. The model is
the same everywhere; only the two outer nodes connect to different
neighbours (a fluid cell with its near-wall conductance, a solid cell with
δ/k, or a boundary condition).

A baffle needs two faces at the same place. Mesh readers and the built-in
mesher keep internal faces of a named face zone and split them into a pair
of boundaries on request (`mesh.baffles = ["plate"]`); the two faces of a
pair are matched by face global id (conformal) or imprinted (D58).

### 9a.2 Through the thickness (1D)

A thin wall is a stack of layers, each with a thickness t_l, a material
(D94, solid or any material with k and c_p) and a number of nodes n_l
(default 1). Every wall face carries n = Σ n_l temperature nodes along its
normal; between nodes the conductance is k A/Δ with k(T) of the material
(harmonic across layer boundaries, plus an optional contact resistance
between layers); each node has the storage ρ c_p A Δ in transient runs and a
volumetric source if the layer has a power (a heater foil, a circuit board).
The two outer nodes connect to the fluid cell (through the near-wall
conductance of §9.1, with the wall function h_w of D135 on a turbulent side)
and to the outer condition or the second fluid cell.

- One layer, one node, steady, constant k: the resistance t/k of §8,
  bitwise the same face flux.
- Several nodes: the transient conduction profile inside the wall (a hot
  gas pulse on a cold wall); the node count per layer is a case key.

### 9a.3 Along the surface (2D, shell conduction)

With `conduction = "shell"` each layer node also conducts along the wall:
a finite-area discretisation on the surface mesh of the thin-wall faces
(D124 W1: edges as sorted node-id pairs, edge lengths, in-surface edge
normals). Across an edge e between faces a and b of the same layer the flux
is

    Q_e = k_l t_l L_e (T_a − T_b) / d_ab

with d_ab the distance between the face centres measured along the surface
(the two half-distances to the edge midpoint, so the surface may bend at the
edge) and the non-orthogonal correction deferred, as for the cells (D24).
Orthotropic in-plane k (a circuit board, a composite) is a 2×2 tensor in the
local face frame. The layer's in-plane and through-thickness conductances
form one system per thin wall.

Edges of a thin-wall face set get an edge condition: `adiabatic` (default),
`temperature`, or `joined` where two thin-wall sets share an edge (a T-joint
of plates: all faces meeting at the edge exchange through it). Where a thin
wall meets a solid region along an edge (a fin root on a meshed base), the
edge exchanges with the solid cells of the adjacent faces; this is open
point TW3 (§18) and refused (TWL004) until decided.

The shell model assumes t much smaller than the size of the faces and the
radius of curvature of the surface; TWL002 warns where t exceeds the local
face size or a tenth of the radius of curvature.

### 9a.4 Discretisation, coupling, solver

- **Unknowns:** the thin-wall nodes are extra unknowns of the energy domain,
  appended to the energy matrix after the cell rows, ordered by (face global
  id, node index) (D93). Conjugate heat transfer stays monolithic (D158):
  fluid cells, solid cells and thin-wall nodes are one matrix and one
  AMG-preconditioned solve. Through-thickness rows are tridiagonal per face
  and strongly coupled; the AMG aggregation keeps each face's stack together.
- **Conserved variable:** h of each node (D92, D150), temperature by the
  same Newton inversion; radiation and convection on the outer face use the
  linearisation of §8.
- **Ledger (D161):** each thin-wall set is a region of the thermal ledger
  with storage, sources, the exchange pair fluid–thin wall per side, the
  outer boundary heat flow and the edge heat flows; the second law gets
  k|∇T|²/T² of the wall and the contact-resistance term of §9.5.
- **Kernels (D160):** one series-conductance kernel for the through-thickness
  links, one edge kernel for shell conduction, both `no_std`, float64,
  generic for dual numbers, CPU = GPU; edge sums in a fixed edge order.
- **Partitioning (D159):** shell conduction across ranks needs the layer of
  wall faces of D124 W2, so W1 and W2 move from L12 into stage 7.
- **Output:** node temperatures per face and layer in the boundary blocks of
  the VTKHDF output (D124 W6) and in the restart; log columns
  `thinwall.<name>.heat_flow`, `thinwall.<name>.max_temperature`,
  `thinwall.<name>.min_temperature`.

### 9a.5 Keys (planned)

```toml
[boundary.duct]
type = "wall"
[boundary.duct.energy]
heat_transfer_coefficient = 8.0
ambient_temperature = 293.15

[boundary.duct.energy.thin_wall]
conduction = "shell"                # "normal" (1D, default) | "shell" (1D + in-plane)
edges = { default = "adiabatic", inlet_flange = { temperature = 330.0 } }

[[boundary.duct.energy.thin_wall.layers]]
thickness = 1.5e-3                  # m
material = "steel"
nodes = 3                           # through-thickness nodes of this layer

[[boundary.duct.energy.thin_wall.layers]]
thickness = 20e-3
material = "mineral-wool"
contact_resistance = 1e-4           # m²K/W to the layer before
power_density = 0.0                 # W/m³, or power (W over the set)

[mesh]
baffles = ["plate"]                 # internal face zones split into boundary pairs
```

Codes: TWL001 error (a thin wall on a boundary type that cannot carry it, or
energy off), TWL002 warning (shell assumption, t against face size and
curvature), TWL003 error (layer without thickness or without k, ρ, c_p in a
transient run), TWL004 error (a shell edge touching a solid region, until
TW3 is decided), TWL005 warning (a baffle pair whose faces do not match
within the D58 tolerance).


## 9b. Wall film heat exchange with walls, thin walls and solids (D190)

Added 2026-11-20 at the maintainer's request. The wall film (D123, with the
sprays L12, after milestone 1) is not in the code yet; stage 7 fixes the
interface so that the film can exchange heat with every kind of wall the
energy equation offers, and so that evaporation cools the wall and
condensation heats it.

### 9b.1 Topology

The film lives on wall faces and, new, on the fluid side of region-interface
faces and of thin walls. On every film face there are three parties in
series: the gas cell, the film (one or two enthalpy nodes, D123 F2), and the
wall side, which is one of

- a plain wall with a prescribed temperature or heat flux (the simple cases
  of D123 (3));
- the first node of a thin wall (§9a), on a wall boundary or on one side of
  a baffle;
- the solid cell of a region interface (§9, full CHT).

The film's wall-side temperature T_fw is the wall-side face temperature of
the wall's energy rows. W1 (surface mesh) therefore covers region-interface
faces and baffle faces as well as boundary faces.

### 9b.2 Heat paths and the sign of the effect

Per film face, positive into the film:

    q_wall = G_fw (T_w − T_fw)                        wall → film (conduction through the film, D123 F2)
    q_gas  = h_g (T_g − T_fs)                         gas → film surface (wall treatment, W4)
    ṁ″ (h_v(T_fs) − h_l(T_fs))                       latent heat leaving with the vapour (ṁ″ > 0 evaporation)
    deposit of impinging parcels: ṁ″_imp h_l(T_d)   (spray cooling)

Evaporation draws the latent heat from the film enthalpy; the film cools,
T_w − T_fw grows and the wall loses heat: the wall is cooled. Condensation
(ṁ″ < 0) releases the latent heat into the film, which heats the wall. The
sign is not prescribed: it follows from the balance per face, so a hot film
on a cold wall heats the wall even while it evaporates. On hot walls the
boiling curve of D123 F2 replaces G_fw: nucleate boiling takes the latent
heat at the wall (q_wall from the boiling correlation), and above the
Leidenfrost temperature the vapour cushion lowers the wall heat flux to the
film-boiling level, so the wall cooling drops sharply there. A face that
dries out keeps its gas-side wall treatment and its ledger entries go to
zero continuously (no jump in the wall heat flux by more than the dry-wall
value).

### 9b.3 Coupling

- The film module solves its own equations (height, velocity, enthalpy
  along the wall); the wall side solves the energy matrix. They exchange
  per face in every outer iteration (steady) or every outer corrector
  (transient) through the linearised flux q_wall = q* + (∂q/∂T_w)(T_w − T_w*):
  the derivative enters the diagonal of the wall's energy row, the film
  enthalpy equation gets the partner term. This is the coupling of D123 (3),
  stated for thin walls and solids.
- With boiling the derivative is large (a steep boiling curve); the
  implicit wall-side term keeps the exchange stable for any ratio of film
  to wall conductance (Giles 1997 stability argument for implicit flux
  coupling). An explicit coupling is not offered.
- The latent heat is not a separate source term: the ledger books enthalpy
  flows of the coefficient sets of D92, so mass and energy close to
  round-off across evaporation and condensation. A derived column
  `film.<name>.latent_heat_flow` reports ṁ″(h_v − h_l) summed over the faces
  for readability.

### 9b.4 Ledger (D161, D123 (5))

New exchange pairs: film–thin wall and film–solid region (in addition to
film–wall, film–gas and spray–film of D123 (5)). Each is one face sum booked
with opposite signs in the two balances, so the wall's thermal ledger shows
the heat taken by the film per boundary, per thin wall and per interface.
The second law gets the entropy of the film–wall heat flow across the
temperature difference, q_wall (1/T_fw − 1/T_w) ≥ 0 per face, and of the
phase change at the film surface.

### 9b.5 What stage 7 builds now for the film

Nothing of the film itself. Stage 7 provides: wall-side face temperatures and
linearised per-face heat-flux slots on wall boundaries, thin-wall faces and
region-interface faces (an extension of D124 W3 and W4 to interfaces and
thin walls); the surface mesh W1 and the wall-face halo W2 on those faces;
the ledger exchange-pair slots. A test module (a prescribed per-face heat
sink with its derivative, standing in for the film) checks the slots in
stage 7 (V7.32).


## 9c. Several fluids and materials: heat exchangers (D194)

**Stage MF-A implemented 2026-11-22 (#206)**: fluid regions with their own
material and temperature, properties per cell and face (flow, energy,
turbulence, boundary conditions), fluid parts as connected components with
REG003 (a part without a pressure boundary in a multi-part case), REG004
(two fluids sharing a face), REG005 (the parts in the log) and REG006
(settings not yet built for several fluids). V7.33 is checked in planar
form (`cht_counterflow_two_fluids`). Build plan of the remaining stages:
MF-B variable density per part and closed parts with their own reference
and p₀ (V7.34), MF-C turbulence, buoyancy and Pr_t per fluid region (MF1),
MF-D thin walls between fluids and the exchanger outputs (V7.33 thin-wall
run, V7.35).

Added 2026-11-20 at the maintainer's request: conjugate heat transfer must
model tube-and-fin and other heat exchangers, with several fluid regions of
different fluids and several solid regions of different materials, and no
limit on their number.

### 9c.1 Regions and domains

- **Any number** of fluid regions and solid regions (D112), each with its own
  material (D94–D103, case or database). No fixed upper bound in the code:
  regions and materials are lists, the per-cell region index is a 32-bit
  integer, and kernels read properties from per-cell (and per-face) fields
  (7.0b), never from a fixed table of a few materials.
- **The flow domain** is the union of all fluid regions; the **energy
  domain** the union of all regions. Fluid regions that touch only through
  solids or thin walls (water in the tubes, air around the fins) are
  separate connected parts of the flow domain: the momentum and pressure
  rows of different parts never couple, the matrix is block diagonal over
  the parts and one AMG setup serves all of them.
- **Pressure level per part:** each connected fluid part needs its own
  pressure reference (an outlet, or a reference cell when closed), its own
  mass ledger and, with D151, its own thermodynamic pressure p₀(t) when it is
  closed. The boundary-condition check (D29, D30) runs per part; a part
  without inflow and outflow and without a reference is refused (new code
  REG003).
- **Two fluids touching directly** (two immiscible liquids at rest, or a gas
  over a liquid) need a multiphase model (L3) and are refused: fluid regions
  of different materials may meet only through a solid region or a thin
  wall (REG004). Fluid regions of the same material may touch (they are one
  part of the flow).
- **Models per fluid part:** the turbulence model, the buoyancy model and the
  turbulent Prandtl number are case-wide by default and may be set per fluid
  region (`[regions.<name>] turbulence = "laminar"` for laminar oil in the
  tubes beside turbulent air; MF1 decided 2026-11-20).

### 9c.2 Interfaces

Every interface between two regions is a region interface of §9 (fluid–solid,
solid–solid) or carries a thin wall (§9a, also fluid–fluid). A fluid–solid
interface face exchanges heat through the near-wall conductance of the fluid
side (D135 with that region's material and turbulence model). Region
interfaces are found automatically from the cell zones; named interfaces
(`[interfaces.<name>]`) only add contact resistances or thin walls.

### 9c.3 Ledger, outputs, checks

- One mass and one kinetic-energy account per fluid part, one thermal
  account per region and per thin-wall set (D119); interface heat flows as
  exchange pairs. The report gives the heat-exchanger numbers per pair of
  fluids: heat duty, inlet and outlet bulk temperatures, the effectiveness
  and NTU against the case's capacity rates, log-mean temperature
  difference, overall conductance UA (log columns `exchanger.<name>.*`,
  with `[[output.exchangers]] hot = "water", cold = "air"`).
- Checks V7.33–V7.35 (§14.3b).
- Codes REG003 (fluid part without pressure reference), REG004 (fluids of
  different materials touching), REG005 info (number of fluid parts, regions
  and materials, listed in the log).


## 9d. Cycle-averaged coupling of periodic flows (D204)

**Core implemented (#190)**: the loop with the mean-flow solid phase (§9d.7);
the skin swing, the fully coupled cycles at the end and the restart (#202,
§9d.8).
Added 2026-11-21 at the maintainer's request as a **requirement for later**:
it is built after transient conjugate heat transfer (PR 7d6) and is not part
of the 7d5 and 7d6 pull requests. Those two are designed so that it fits
(§9d.6).

### 9d.1 The problem

An engine exhaust manifold sees a gas flow that is periodic over the engine
cycle (720° of crank angle for a four-stroke engine, 40 ms at 3000 rpm).
The gas temperature near the wall swings by several hundred kelvin within a
cycle; the structure needs minutes to warm up. The temperature swing in the
solid decays with the depth y as e^(−y/δ), with the skin depth of a harmonic
of angular frequency ω

    δ = √(2a/ω),   a = k/(ρ c_p)

(Carslaw & Jaeger 1959, §2.6). For stainless steel (a ≈ 4·10⁻⁶ m²/s) at
3000 rpm the fundamental has δ ≈ 0.23 mm, the tenth harmonic of the sharp
blow-down pulses δ ≈ 70 µm; cast iron is about twice as deep. Below this
skin the structure carries an almost steady temperature field. A fully
coupled transient conjugate run would need thousands of cycles to bring the
structure to its thermal state: far too costly.

### 9d.2 The practice it replaces

The maintainer's earlier practice: (1) run the transient gas flow over one or
more cycles with fixed wall temperatures; (2) cycle-average the gas-side heat
transfer, keeping its distribution over the wall and the energy balance;
(3) map it onto a separate structure mesh and solve the structure steady;
(4) give the new wall temperatures back to the gas side; two or three
repetitions. Its weak points: two meshes and an interpolation that does not
conserve the heat exactly, files moved by hand, a fixed number of loops
without a measure of convergence, the averaging definition left to the user
(mean flux, or h and a reference temperature), the start-up cycles counted
in the average, and no estimate of the surface temperature swing (thermal
fatigue, deposits).

### 9d.3 Requirement

One case, one mesh, one command; the loop runs by itself:

1. **Fluid phase.** The flow and the fluid energy run transient over
   `warmup_cycles` plus `average_cycles` cycles. The solid regions are
   frozen: each interface face is a wall at the current solid surface
   temperature T_w,f, constant in time (or the swing of item 5).
2. **Averaging.** Per interface face the run accumulates over the averaging
   cycles ∫ h_f dt and ∫ h_f T_ref,f dt, with h_f the fluid-side conductance of
   the face (1/R_P of §9.1: k/δ_P laminar, h_w of the thermal wall function,
   D135) and T_ref,f the fluid temperature that the face flux acts on.
   h̄ = ⟨h⟩ and T̄_ref = ⟨h T_ref⟩/⟨h⟩. Because T_w is constant over the
   cycle, h̄ (T̄_ref − T_w) equals the cycle-mean heat flux of the face
   exactly: the distribution and the energy balance of the cycle are kept,
   and the Robin form carries the first-order response of the gas side to a
   change of T_w.
3. **Solid phase.** The energy equation on the solid regions (the energy view
   restricted to the solids, the same interface faces and matrix rows as the
   conjugate solve, D158) with the Robin conditions (h̄, T̄_ref) on the
   interface faces; steady by default, or a transient march with large steps
   for warm-up and cool-down histories (the multi-rate schedule of D118).
4. **Update and stop.** New T_w per interface face (optional under-relaxation),
   then back to the fluid phase, which continues from its last state (no new
   start of the flow). The loop stops when the largest change of T_w is
   below `tolerance` and the interface heat flow changes by less than
   `heat_tolerance`, or after `max_updates` (warning). Each fluid phase also
   checks that the flow is periodic (the mean interface heat flow of the last
   two cycles within `cycle_tolerance`) and adds cycles up to `max_cycles`
   when it is not.
5. **Surface temperature swing** (`skin = "effusivity"`, optional). The
   Fourier coefficients q_n of the interface heat flux over the averaging
   cycles give, through the periodic solution of a semi-infinite solid,
   the surface swing T_w(t) = T̄_w + Σ_n Re[q_n e^(iω_n t) / (e √(iω_n))],
   e = √(kρc_p) of the solid at the face. The next fluid phase uses it as the
   wall temperature history, and it is written per face (amplitude, minimum
   and maximum of T_w) for thermal-fatigue and deposit studies. Valid where δ
   is small against the wall thickness and the radii of curvature, checked
   per face.
6. **Final resolved cycles** (`resolved_cycles`, optional): after
   convergence, one or two cycles of the fully coupled transient conjugate
   solve (7d6), started from the converged solid field, give the exact swing
   also at thin walls and edges for the cost of a few cycles.
7. **Outputs** (agents first, D77–D87): `cycle_coupling.csv` with one row per
   update (cycles run, interface heat flow per interface, largest and mean
   change of T_w, mean T_w); per-face fields h̄, T̄_ref, q̄, T_w and the swing;
   a `cycle_coupling` section in `result.json`; the ledger per update (the
   solid's heat balance to the solver tolerance; the cycle-mean balance of the
   fluid with the change of stored heat over the cycle). Every setting is a
   case key; nothing needs a window.
8. **Restart** stores the phase, the update and cycle counters, the
   accumulators and T_w per face, so a run can stop and continue in any
   phase.
9. **Other periodic flows** use the same loop: rotor–stator heat transfer
   (blade-passing period), pulsating jets, pulsed burners, reciprocating
   compressors. Engine port conditions come as tables over the crank angle
   (from a 1D engine simulation), periodic in the cycle; a 1D code coupled at
   run time goes through D105.

Why this is better than the old loop: the Robin pair (h̄, T̄_ref) instead of a
mapped flux makes the structure respond correctly to its own temperature,
so the loop converges in fewer updates (the remaining coupling is only the
dependence of h on T_w through the gas properties); one mesh and one set of
interface faces conserve the heat per face exactly; convergence is measured
and the loop stops by itself; start-up cycles are excluded and the cycle
periodicity is checked; the surface swing comes as a result instead of
being guessed.

### 9d.4 Keys

As built (#190, #202, #208) the table is `[energy.cycle]` (beside
`[energy.acceleration]` and `[energy.frozen_flow]`, which it excludes).

```toml
[energy.cycle]
period = 0.04              # s, a whole number (>= 2) of time steps
# enabled = true           # false: the same case fully coupled (its reference)
warmup_cycles = 1          # cycles before the averaging, after every update
average_cycles = 1
max_cycles = 6             # per fluid phase, when the flow is not yet periodic
cycle_tolerance = 0.005    # relative change of the mean interface heat flow between cycles
max_updates = 5
tolerance = 0.5            # K, largest change of T_w per update
heat_tolerance = 0.005     # relative change of the interface heat flow per update
relaxation = 1.0           # of the first update
update = "quasi-newton"    # later updates: "quasi-newton" | "aitken" | "relaxation"
solid_model = "mean-flow"  # "mean-flow" | "robin" (the first-cell Robin pair)
solid = "steady"           # "transient": a history, marched under the cycle-mean flow (#208)
march_cycles = 20          # with "transient": periods of physical time per update
march_steps = 4            # with "transient": steps of the march per update
skin = "none"              # "effusivity": the surface swing of a thin skin (#202)
skin_harmonics = 10        # harmonics of the period for the swing
resolved_cycles = 0        # fully coupled cycles after the loop (#202)
```

Codes: CYC001 error (a setting that cannot be used: a steady case, no
solid, with `[energy.acceleration]`, an adaptive step, a period that is
not a whole number of steps, a value out of range, `march_cycles` below
the cycles of a fluid phase or without `solid = "transient"`, a
`relaxation` other than 1 with it), CYC002 warning (the
fluid phase is not periodic after `max_cycles`), CYC003 warning (skin depth
larger than a third of the wall thickness V/A), CYC004 warning (stopped at
`max_updates` before the tolerances; not for `solid = "transient"`, whose
end at `end_time` is the end of the history), CYC005 note (a step of the loop);
LIM009 in the limitations of every run with the loop.

### 9d.5 Checks (planned)

V7.38–V7.40 in §14.3c.

### 9d.7 As built: the mean-flow solid phase (#190)

The first build used the Robin pair of item 3 with h̄ the conductance of
the first fluid cell. It converges, but slowly: in the gallery channel,
with the quasi-Newton update, the wall change per update grows from 1.6 K
to 6.7 K over the first six updates and is still 0.14 K after twelve; its
cycle mean is within 0.1 K of the periodic state at 1450 s simulated,
later than the fully coupled run itself (1350 s). The first-cell conductance (205 W/K for the whole wall there) says
how the heat flux answers a change of T_w with the fluid held; the fluid
does not hold: the near-wall fluid follows the wall and the bulk warms
downstream, so the true response of the cycle-mean heat flow is that of
the flow (about ṁ c_p = 8 W/K and h A = 30 W/K there). A solid phase that
sees 25 times the true conductance moves the wall too little per update.
This is the known weak point of the mapped Robin loop of the industry
practice, and it is why it needs relaxation and many repetitions.

The default solid phase (`solid_model = "mean-flow"`) therefore solves the
whole conjugate energy equation, fluid and solid, steady, with the fluid
convected by the cycle-mean mass flux of the averaging cycles (∫ṁ dt per
face over the cycles): solve A holds the solids at their current T_w and
gives the mean-flow heat flow Q_m per interface face; solve B frees
everything and adds the defect Q̄ − Q_m of each face (the cycle-mean heat
flow of the transient phase minus that of the mean flow) as a source to
the solid cell of the face. B's wall then answers its own temperature with
the response of the mean flow, and at the converged T_w the defect makes
its heat flows those of the periodic flow. The fluid temperatures of A and
B are thrown away; the fluid phase continues from its own last state.

Results on the gallery case (`cht_pulsating_channel_cycle.toml`, 80 × 20
cells, 2.5 s steps, period 50 s): converged after 4 updates and 10 cycles
(500 s simulated, 5.5 s on 4 threads); interface heat flow 160.0002 W of
160 W; probe temperatures within 0.14 K (0.4 % of the 33.7 K rise) of the
fully coupled run's cycle means at 4000 s, which is within 0.1 K of its
periodic state from 1350 s (about 12 s of its 36 s). The small rest is
the correlation of h(t) with the swing of T_w, which the held wall does
not have (the skin of PR C). The ledger: the heat balance of every fluid
step with the solids held (heat into the held cells booked as `held`) and
of every solid phase (the Robin heat or the defect booked as `robin`)
closes to 10⁻¹³ relative.

Updates after the first: IQN-ILS quasi-Newton (Degroote, Bathe &
Vierendeels 2009) on the wall temperatures by default, Aitken (Irons &
Tuck 1969; Küttler & Wall 2008) or fixed relaxation on request. Checks:
V7.39 and V7.40 pass (`tests/cht_cycle.rs`); V7.38 is not built yet (a gap
in features.toml).

### 9d.8 As built: skin swing, fully coupled cycles and restart (#202)

**Skin swing (`skin = "effusivity"`).** Over the averaged cycle the loop
accumulates per interface face the harmonics c_n = (2/T)∫q e^{−iω_n t} dt
of the heat flux into the solid per area (n ≤ `skin_harmonics` and half
the steps of a cycle), and ∫q² and ∫(T_w − T_s*) q for the resistance of
the solid half cell r_s (least squares). A semi-infinite solid of
effusivity e = √(k ρ c_p) answers with the surface swing θ_n = Z_n c_n,
Z_n = (1 − i)/(√2 e √ω_n) (Carslaw and Jaeger, §2.6); the held cell takes
(Z_n − r_s) c_n, area-weighted over its faces, around its mean at every
step of the next fluid phase. The swing changes the flux it comes from,
so its amplitudes are unknowns of the quasi-Newton update next to the
wall temperatures (a lagged swing alone oscillated: 0.16, 0.49, 0.34,
0.08 K); the loop converges when the change of the swing (`swing_change`)
is also below `tolerance`. ρ c_p and k come from the materials at the
current temperature (the face values of the energy equation hold the
fluid's until the first refresh, which made a first version 12 % off).
CYC003 where δ = √(2k/(ρ c_p ω)) exceeds a third of V/A. The heat the
swinging held cells store is booked with the `held` heat, so the ledger
closes in every step.

**Fully coupled cycles (`resolved_cycles`).** When the loop ends the
solids are released and run that many cycles fully coupled; per face the
heat flow and the range of T_w over the last cycle are reported
(`resolved_t_wall_min`, `resolved_t_wall_max`). With the skin swing every
solid cell within six skin depths of a held cell (found by a
breadth-first walk through the solid cells, depth along the normal of the
held cell's interface face) starts from s_n e^{−(1+i) d/δ_n} around its
mean, on its temperature and both old time levels. Released from the mean
alone the energy of the missing swing left only with the time constant of
the whole wall: 0.3 K off after three cycles on the gallery wall, against
0.016 K with the profile.

**Restart.** `restart.h5` holds for a conjugate case the temperature of
the full mesh and its old level by cell identity (`cht_cell_id`,
`cht_temperature`, `cht_temperature_old`) and the state of the loop as
JSON bytes (`cht_cycle_state`: accumulators, IQN history, held means,
swing, resolved cycles), floats round-tripped exactly. Continuations are
bit-identical in a fluid phase, with the swing and in the resolved
cycles; a loop state with another cell or face order is refused, a
restart with `[energy.acceleration]` too (ENE006).

Results on `cht_pulsating_channel_skin.toml` (80 × 80 cells, δ = 0.030 m in
a 0.2 m wall, time constant about 700 s): 10 updates (1100 s), then three
fully coupled cycles; at 1250 s every probe is within 0.016 K of the
fully coupled run's periodic state (10 000 s) and the swing within 0.1 %;
the fully coupled run is within 0.1 K from 4400 s (skin) to 5550 s (heated
face). Checks: `tests/cht_cycle.rs` (skin swing and resolved cycles,
CYC003, bit-identical restarts, steady restart).

### 9d.9 As built: the transient solid phase (#208)

`solid = "transient"` turns the loop into a warm-up or cool-down history
(item 3, the multi-rate schedule of D118). Each update marches, instead
of the steady solve B, fluid and solids together under the cycle-mean
mass flux and the defect of every interface face (with `solid_model =
"robin"` the solids under the Robin pairs) over M = max(`march_cycles`,
cycles of the fluid phase) periods, in `march_steps` BDF steps (the first
BDF1, then BDF2) from the temperatures at the end of the fluid phase.
The march covers the fluid phase that has run; the clock then jumps by
the rest, J = M − (time of the phase), a whole number of periods, so the
next fluid phase starts at the same phase of the inlet history and the
boundary values stay periodic. A march that would pass `end_time` is cut
there. The update has no quasi-Newton or relaxation (the march is the
history, not a fixed-point iterate), so `relaxation` must be 1. Boundary
values and sources are held at the time of the update during the march.
`cycle_coupling.csv` gains `march_time` and `march_jump`; the wall
temperatures at time + `march_jump` form the history, and the loop stops
when they change by less than `tolerance` per update, or at `end_time`
(no CYC004). The restart keeps the jump (the clock of `restart.h5` is the
jumped one).

Results on `cht_pulsating_channel_warmup.toml` (the pulsating channel
of §9d.7 with a period of 5 s, 80 × 20 cells, 0.25 s steps, 20 periods
per update in 5 steps of 20 s, skin swing): 15 updates to 1500 s, 780
steps in 18 s on 4 threads, against 6000 steps in 95 s fully coupled.
The wall at the end of every march is within 0.21 K of the fully coupled
run's cycle mean there over a rise of 35 K (largest at the first march,
at 100 s), and within 0.01 K at 1500 s; the balance of every march
closes to 2·10⁻¹² relative. Checks V7.44 (`tests/cht_cycle.rs`: within
1 % of the whole rise on 40 × 10 cells, bit-identical restart in the
fluid phase after a jump).

### 9d.6 What 7d5 and 7d6 keep open for it

- The fluid-side conductance of every interface face (laminar k/δ_P or the
  wall-function h_w of 7d5) is one number per face computed in one place,
  so the averaging of item 2 reads it without a second formula.
- The transient conjugate solve of 7d6 keeps the solid regions' rows in the
  same matrix with their own time terms; freezing the solids (item 1) and
  solving the solids alone (item 3) are restrictions of that matrix, not new
  solvers.
- The transient heat balance of 7d6 has a storage term per region, which the
  cycle-mean fluid balance of item 7 needs.


## 9e. Thermal acceleration of solids (D213)

Added 2026-11-22 at the maintainer's request (built in #180). The earlier
practice of the maintainer: reduce ρ c_p of the structure in the first
cycles, so that its thermal transient is shorter, as an initial flow field
shortens the flow's.

**What it changes.** For a factor f on ρ c_p of a solid: the steady
temperatures and the steady interface heat flows do not change (no ρ c_p in
the steady equation); the time constants ρ c_p L²/k and ρ c_p V/(hA) are f
times shorter; the skin depth of a period √(2k/(f ρ c_p ω)) and the surface
swing q̂/(e√ω), e = √(k f ρ c_p), are 1/√f times larger; the cycle mean of a
periodic flow is unchanged for linear physics and biased to second order in
the swing through h(T_w), μ(T), c_p(T) and radiation; the history before
the restore is not physical. ρ c_p of a fluid is never scaled (it also
multiplies the convection).

**Restore.** T stays continuous. In a steady flow a jump is exact. In a
periodic flow a jump catches the accelerated solid in its swing and leaves
that part in the deep temperature, where it decays only with the true time
constant: about (swing)/(ω τ_f) in a lumped model (8 % of the swing at
τ_f = 2 periods); a geometric ramp of f over n periods lets it die out
(0.3 % for n = 3, 0.05 % for n = 10). Default: a ramp over 3 periods with a
`period`, a jump without.

**Automatic factor.** After the first window with the true ρ c_p, each
accelerated region's conductance G to the rest: the heat flow out of the
region through its interfaces with the region held 1 K warmer, minus that
with it held as it is (steady energy solves, every solid held, the flow as
it is), plus Σ A k/d_n at outer temperature faces and Σ A h k/d_n/(h + k/d_n)
at convection faces. τ = C/G with C = Σ ρ c_p V; f = target/τ in
[10⁻⁴, 1]. Internal conduction is not in τ (a small Biot number); where it
matters f comes out larger, the acceleration milder.

**Settling.** Per accelerated region the window-mean temperature (a window
is a period, or `window_steps` steps); settled when the change of the last
window and the change still to come, extrapolated geometrically from the
last two windows (Aitken's Δ²), are below `settle_tolerance` of the rise
since the start, for every region, at least three windows after the
factor was set.

**Ledger.** The equation stores f times the heat of the true solid; per
region and step `inertia_scaling` = (1 − f) × true storage (W) in
`heat.csv` and the energy summary, the total per region over the run (J)
in the report. The imbalance of the solved equation is unchanged.

**Keys.** `[energy.acceleration]`: `inertia_scale` (number or `"auto"`),
`regions`, `restore` (`"settled"` or a time), `period`, `window_steps`,
`settle_tolerance`, `restore_ramp`, `target_windows`, `settle_windows`.
Codes: ACC001 error (steady case, no solid, a fluid region, values out of
range), ACC002 warning (ended before the restore), ACC003 note (skin depth
of the period larger than V/A of a region while accelerated), ACC004 note
(each step of the acceleration in the log and the run record).

**Relation to D204.** The solid phase of the cycle-averaged loop (a steady
solid with the averaged Robin pair) is the limit f → 0 of this acceleration
with the swing averaged out; the loop needs no factor. Checks V7.41, V7.42
(§14.3c).

## 9f. Frozen-flow march of transient conjugate runs (D118)

Built in #208 at the maintainer's request (D118 item: multi-rate clocks,
frozen flow; the third way to shorten a thermal warm-up beside the
acceleration of §9e and the loop of §9d).

**The problem.** A transient conjugate run steps the flow and the
temperatures with one step, set by the flow (a Courant number of a few to
ten). The wall warms over minutes or hours; the flow settles within a few
flow-through times. Fully coupled, every step of the warm-up is a flow
step: the heated wall of V7.37 with a 0.5 s step needs 6000 steps to
3000 s.

**The march.** `[energy.frozen_flow]` steps fully coupled until the flow
has settled, then holds velocity, pressure, mass flux and turbulence and
marches only the energy equation of fluid and solids (the one matrix of
D112, with its time terms, convected by the held mass flux) with
`time_step`, a multiple of the `[time]` step. Each frozen step solves the
energy equation repeatedly (`iterations`, 20 at most) until its change is
below 10⁻⁸ of the temperature span, for k(T), c_p(T) and the interface
conductances. The first frozen step after coupled steps is BDF1 (the old
levels belong to the flow's step), later ones BDF2; the first coupled
step after a frozen phase takes the BDF2 weights of the unequal steps,
which make it close to BDF1.

**Settling.** Every `window_steps` coupled steps the run compares the
velocity field with that of the last window: c = rms|u − u₋|/rms|u|. The
flow has settled when c and the change still to come, c r/(1 − r) with
r = c/c₋ (Aitken's Δ², as in §9e), are below `settle_tolerance`. For a
flow that follows the temperature, which drifts with the warming fluid,
c alone counts. `start` may instead give a time for the first hold.

**Refresh.** A flow follows the temperature when the viscosity depends on
it (a property form of D94–D103) or with buoyancy. Such a flow is
refreshed, by coupled steps until it has settled again, when a fluid
temperature has moved by more than `refresh_change` since the hold, or
every `refresh_every` frozen steps (also for a flow that does not follow
the temperature, e.g. time-dependent boundary values, FRZ004). A flow
that does not depend on the temperature is held once and never refreshed:
the march is then exact up to the time error of the long steps.

**Ledger and restart.** The heat balance closes in every step, frozen or
coupled (`heat.csv`). `restart.h5` holds the state of the march
(`cht_frozen_state`: phase, windows, reference temperatures, counters)
and, for a viscosity that follows the temperature, μ per face of the
fluid (`cht_viscosity`), which a CHT restart lacked before: a restart in
a frozen phase and in the coupled steps of a refresh continues bit for
bit.

**Keys.** `[energy.frozen_flow]`: `enabled` (false: the same case fully
coupled, its reference), `time_step`, `start` (`"settled"` or a time),
`window_steps` (10), `settle_tolerance` (10⁻³), `refresh_change` (2 K),
`refresh_every` (0: off), `iterations` (20). Codes: FRZ001 error (a
steady case, no solid, with `[energy.acceleration]` or `[energy.cycle]`,
variable density, an adaptive step, `time_step` below the `[time]` step,
a value out of range), FRZ002 warning (the flow did not settle), FRZ003
note (each hold and refresh, in the log and the run record), FRZ004
warning (time-dependent flow values); LIM010 in the limitations.
`result.json` has `frozen_flow`: the phases (kind, start, end, steps,
reason), frozen and coupled steps and time, refreshes, solves per frozen
step.

Results on `cht_heated_wall_frozen.toml` (the heated wall of V7.37,
80 × 30 cells, with a viscosity that halves over 17 K; flow step 0.5 s,
march step 10 s, `refresh_change` 2 K): held at 55 s, 7 refreshes, 557
steps (287 frozen, 270 coupled) in 17 s against 6000 steps in 155 s
fully coupled; every probe within 0.058 K of the fully coupled history
over a rise of 21 K, 0.03 K at the end (0.002 K with `refresh_change =
0.5`, 794 steps); the balance closes to 3·10⁻¹² relative. Checks V7.43
(`tests/cht_frozen.rs`).

**Relation to §9d and §9e.** The acceleration (§9e) changes the solid's
ρ c_p and keeps the flow step; the frozen flow keeps the physics and
takes long steps where the flow does not change; the loop (§9d) handles
a flow that never settles (periodic), and with `solid = "transient"`
uses the same idea, a march under the held (cycle-mean) flow. The three
exclude each other in one run.

## 10. Coupling strategy, linear solvers and order of solution (D158)

### 10.1 One matrix or two solvers

D112 decided that conjugate heat transfer is one temperature matrix. The
prompt asks for a comparison with partitioned coupling; the comparison
confirms D112.

| Criterion | Monolithic (one matrix, D112) | Partitioned (Dirichlet–Neumann or Robin–Robin with Aitken or IQN-ILS relaxation) |
|---|---|---|
| Conservation at the interface | exact per face at every iteration (one flux per face) | only when the coupling iterations converge; a residual mismatch per step otherwise |
| Stability | unconditional (implicit) | Dirichlet–Neumann is stable only when the Neumann side is the one with the larger conductance, and loses stability as the ratio approaches the other way (Giles 1997); Robin–Robin needs a well-chosen coefficient |
| Cost per time step or outer iteration | one linear solve on the energy domain | 5–30 coupling iterations at high ratios even with Aitken (Küttler & Wall 2008) or IQN-ILS (Degroote, Bathe & Vierendeels 2009) |
| Steady convergence | solid rows converge with the fluid in the same solve | the solid with its long diffusion time slows the coupling unless relaxed carefully |
| GPU | one AMG-preconditioned solve, no host round trip (D76) | interface exchange and relaxation per coupling iteration, reductions on the host |
| Different time scales (long thermal transients) | multi-rate with frozen flow (D118): the flow is held, the same energy-domain matrix is marched with large steps; still one matrix | natural, but this case is covered by D118 |
| Separate codes or external models | not applicable | the only route (P10, D105) |

Decision: monolithic only in stage 7. Partitioned coupling is used only
where the partner is outside the solver (external system models, P10); no
partitioned fluid–solid coupling inside the solver is planned.

### 10.2 Linear solver for the energy domain

- **Matrix:** non-symmetric in fluid rows (convection), symmetric in solid
  rows. BiCGStab (D36) preconditioned by the aggregation AMG (stage 6's
  AMG-preconditioned BiCGStab for scalars, D36), built on the energy-domain
  view, rebuilt when the coefficients change by more than a set fraction
  (`solver.energy.amg_rebuild`), float32 hierarchy (D37).
- **Coefficient jumps:** the strength of connection is measured on the
  symmetrically scaled matrix, |a_ij| / √(a_ii a_jj) ≥ θ (Vaněk, Mandel &
  Brezina 1996). Across a copper–air interface this ratio is about
  √(k_f/k_s) ≈ 10⁻², below the default θ = 0.04, so aggregates do not cross
  the interface and the coarse spaces keep the nearly isothermal solid as
  its own aggregates. Requirement: the iteration count of the composite
  cylinder (V7.5) at k_s/k_f = 10⁴ is at most 1.5 times that at ratio 1, and
  grows by at most 30 % per mesh doubling.
- **Solid-only domains** (pure conduction) take CG with the same AMG.

### 10.3 Order of solution per outer iteration

Steady (SIMPLEC, D19) and transient (PISO with outer correctors, D28):

1. momentum predictor with the body force of the current ρ or T (D154);
2. pressure correction(s) and flux correction with the balanced-force
   Rhie–Chow;
3. turbulence (stage 6) with G_b from the current T and ∇T;
4. energy on the energy domain (δT form, §3.2), with the interface
   temperatures of §9 from the previous iteration;
5. properties (D99), p₀(t) for closed volumes (D151), face densities;
6. convergence test, ledger at the cadence of D89.

The buoyancy coupling is explicit within one outer iteration and converged
by the outer iterations. A transient run with strong buoyancy needs outer
iterations for second order in time (verified by V7.21); without them the
buoyancy is lagged by one step and the run log says so (ENE005 info).
Option `solver.energy.position = "after-flow"` (default, order above) or
`"before-momentum"` (energy first, for flows driven only by buoyancy where
the temperature leads).


## 11. Multi-region partitioning and halos (D159)

Stage 7 runs its acceptance tests on one rank with threads; multi-rank
acceptance is stage 8. The multi-region data layout is built for ranks now:

- **Partition weights (D122 (1)):** METIS with two constraints, flow work
  (fluid cells 1, solid cells 0) and energy work (all cells 1), so each rank
  holds a balanced share of both; results never depend on the partition
  (D93 (2)).
- **Views per rank (D112):** each equation domain has a local view per rank
  with its own halo plan: the flow halo holds fluid cells only, the energy
  halo fluid and solid cells. A halo cell carries its region id, so the
  one-sided gradients of §9.2 work at rank boundaries.
- **Interface faces on rank boundaries** are computed once from the same
  inputs on both ranks (D93 (3)), as every processor face.
- **AMG on views:** the aggregation of D93 (5) on the view's global graph;
  aggregates may cross rank boundaries, never region interfaces with weak
  connections (by the strength rule, not by a special case).
- **Non-conformal parts (D58):** imprinting runs before partitioning, so the
  partitioner sees one conformal mesh.
- **Tests now:** the decomposition into 1, 2 and 4 parts of a fluid–solid
  mesh gives the same views, halo plans and assembled energy matrix in global
  order (a unit test on the decomposition of `foehncraft-mesh`), so stage 8
  only adds the communication.


## 12. Kernels: GPU, reproducibility, adjoint (D160)

- **Kernels:** `no_std` bodies in `foehncraft-kernels` (D60): the energy face
  flux (convection, conduction, GGDH, interface conductance), the cell terms
  (storage, sources, viscous heating, pressure work), the thermal boundary
  faces, the reduced face jump and the cell net force of D154, the buoyancy
  production of D155. Each a pure function (D59 (1)), generic over the real
  type (dual numbers for the Jacobians), float64 (D7).
- **No divergence on the GPU:** interior faces of the same region and region
  interfaces use one formula (the series conductance; inside a region
  R″_c = 0 and the arithmetic k of D99 is selected by a per-face factor
  s ∈ {0, 1} in a product, not a branch); gradient stencils take the
  neighbour value or the face temperature by a select. Boundary faces carry
  a kind tag per face as today. No virtual dispatch (D99).
- **Reproducibility:** gather-based assembly in global face order (D31, D93
  (1)), no floating-point atomics; the damping factor of §3.3, p₀ of D151
  and the ledger sums are exact reductions; interface fluxes are computed
  once (D117 (1)). CPU and GPU bitwise equal on the verification cases
  (V7.20), 1 and 4 threads bitwise equal, restart bit-identical.
- **Smooth functions for the adjoint (D59 (3)):** positive part
  P⁺(x) = ½(x + √(x² + δ²)) with δ relative to the local scale, smooth
  absolute value √(x² + δ²) in C_ε3, analytic derivatives of the radiation
  linearisation and of all property forms (D99); the upwind switch and the
  limiters keep their frozen variants of D59. Dual-number Jacobians of every
  new kernel against central differences, and ⟨Av, w⟩ = ⟨v, Aᵀw⟩ for the
  energy matrix, are tests (V7.22).


## 13. Energy ledger and second law (D161)

- **Quantities (D88, D119):** thermal energy ∫ρh dV per region with
  convective boundary flows ṁh, conductive heat flows per boundary and per
  interface, named sources (S_h, pressure work, radiation later); kinetic
  energy as in stage 5. Viscous heating Φ (and ρε in RANS) is an exchange
  pair kinetic → thermal, so the two ledgers close together. Gravity work
  ∫(ρ − ρ_ref) g·u dV is a named source of the kinetic-energy ledger. With
  Boussinesq the buoyancy work has no counterpart in the thermal ledger:
  the approximation is not energy-consistent, so the term is booked as a
  source without a partner and the report says so (not hidden, D88).
- **Closure:** steady: heat in = heat out to the solver tolerance (the
  default tolerance of D89 relative to the largest heat flow); transient:
  per step to the tolerance of D89.
- **Second law (D90 parts 1–3 become active):** per cell
  σ = k|∇T|²/T² + Φ/T, in RANS with the turbulent parts of Kock & Herwig
  (2004): ρε/T for the dissipation and (μ_t c_p/Pr_t)|∇T|²/T² for the
  turbulent heat flux (GGDH: ∇T·K·∇T/T², non-negative when K is made
  positive semi-definite, else a health number); contact resistances as in
  §9.5. The specific entropy s(T, p) comes from the same coefficients (D92
  (2)). The numerical entropy production is reported as in D90 (2).
- **Thin walls and the film (§9a, §9b):** every thin-wall set is a region of
  the thermal ledger (storage, sources, side exchange pairs, outer and edge
  heat flows); the exchange pairs film–wall, film–thin wall and film–solid
  have reserved slots now and are filled with the film (D123, D190).
- **Reports:** `balance.csv` gets the thermal columns per boundary, interface,
  thin wall and region; the end-of-run plot shows the heat flows (D89).


### 13.1 Turbomachines: compressors, turbines, turbochargers (D184)

The performance of rotating machinery of stage 6 (D184: torque, shaft power,
total pressure rise, head, hydraulic efficiency, Euler work, the angle of
every sliding rotor) is built for incompressible flow. The energy equation
adds what a turbocharger needs, as outputs of `[[output.machines]]` with the
same inlet, outlet and rotor:

- **Enthalpy rise** Δh_t = h̄_t,out − h̄_t,in of the mass-flow weighted total
  enthalpy h_t = h + ½|c|² (absolute velocity), and the shaft power balance
  P = ṁ Δh_t + Q̇ (the heat flow through the walls, the ledger of section
  13); the rothalpy I = h + ½|w|² − ½|u|² of D183 is constant along a relative
  streamline without losses and heat transfer.
- **Pressure ratio** Π = p̄_t,out / p̄_t,in (total to total) and total to
  static, with the mass-flow weighted total pressure of the ideal gas.
- **Isentropic efficiency** η_is = (h_t,out,is − h_t,in)/(h_t,out − h_t,in) for
  a compressor (the inverse for a turbine), with the isentropic outlet state
  from Π and γ (or from the property tables for a real gas); **polytropic
  efficiency** from the same states.
- **Corrected quantities** for the map of a machine: corrected mass flow
  ṁ √(T_t,in/T_ref)/(p_t,in/p_ref) and corrected speed N/√(T_t,in/T_ref).
- **Thermal effects of a turbocharger:** the heat flow from the turbine side
  through the bearing housing to the compressor side (conjugate heat transfer
  of the solid regions, D157) enters the compressor enthalpy rise; the
  efficiencies are reported with and without it ("adiabatic" and "diabatic"
  efficiency).
- **Compressibility:** the compressor and turbine need the high-Mach
  compressible flow (L-series, after milestone 1); the low-Mach energy
  equation of stage 7 covers fans, pumps and the compressor at low pressure
  ratio only. The outputs are specified now so that the case file and the log
  columns (`machine.<name>.enthalpy_rise`, `pressure_ratio`, `efficiency_is`)
  do not change when the compressible solver arrives.
- **Operating line** of a machine map: a study (D62) over the speed and the
  outlet pressure or mass flow, one row per point with the quantities above.

## 14. Verification and validation (D162)

Every case runs on hexahedral, tetrahedral, polyhedral and mixed meshes
where it is a verification case (perturbed meshes included), closes the
energy ledger (§13) and the second-law check, and ships with a gallery
example, case files and a recipe (D18). The discretisation uncertainty is
computed as in D67 (Celik et al. 2008): three meshes with refinement ratio
r ≥ 1.3 (2 where the mesh generator doubles), the observed order p from the
three solutions, the extrapolated value, GCI_fine with the safety factor
1.25, and the asymptotic-range check GCI_32/(r^p GCI_21) = 1 ± 0.1.

Pass bands: "order" means the observed order of the L2 error of T (and of
the named quantity) on the finest pair; "GCI" means GCI_fine of the named
quantity; "reference" means the extrapolated value against the reference.
Reference values marked "from the paper" are entered into the case's
reference file from the publication when the case is built and checked by
the maintainer; values given here are those of the cited papers.

### 14.1 Verification (exact or analytic solutions)

| # | Case | Checks | Pass band |
|---|---|---|---|
| V7.1 | Manufactured solution, energy only, k(T) and c_p(T) polynomial, ρ constant, imposed velocity field, Dirichlet and Robin boundaries | T, wall heat flow | order 1.8–2.2 on all four types; ledger closes to 10⁻¹⁰ relative |
| V7.2 | Couette flow with viscous heating (one wall isothermal, the other adiabatic): parabolic T profile, Brinkman number 1 | T profile, wall heat flux = μU²/H | order 1.8–2.2; wall heat flux = dissipated power to the solver tolerance |
| V7.3 | Adiabatic stirred box (lid-driven cavity, all walls adiabatic, transient) | rise of ∫ρh dV against the lid work minus the change of kinetic energy; entropy rise against ∫Φ/T | ledger per step to the tolerance of D89; entropy production ≥ 0 per cell |
| V7.4 | 1D composite wall of three solids, k ratios 1, 10², 10⁴, with and without contact resistance, on meshes whose cells do not line up with the layers in the transverse direction | T, interface temperatures, heat flux | error ≤ 10⁻⁹ relative to ΔT on every mesh type, perturbed (exact for piecewise-linear T, §9.2) |
| V7.5 | Composite cylinder, radial conduction through two materials with a curved interface, ratios 1 to 10⁴, Robin outer surface | T, interface heat flow | order 1.8–2.2; heat flow equal on both sides to round-off; AMG iterations as in §10.2 |
| V7.6 | Two semi-infinite bodies brought into contact (copper and air; and two equal solids): constant interface temperature T_i = (e₁T₁ + e₂T₂)/(e₁ + e₂), e = √(kρc_p), erfc profiles (Carslaw & Jaeger 1959, §2.15) | T_i(t), profiles | T_i within 0.1 % of ΔT after the first 10 steps; time order 1.8–2.2 (BDF2, Crank–Nicolson) and space order 1.8–2.2. **Met 2026-11-21 (#159)** with a solid against a fluid at rest (`cht_contact_transient.toml`): T_i within 3·10⁻³ K after 10 steps, orders 2.02/2.00 in space and 2.05/2.02 in time; copper–air, and equal bodies at 350 K to round-off (`tests/cht.rs`) |
| V7.7 | Hydrostatics: closed box at rest, (a) Boussinesq with a stable linear stratification, (b) variable-density isothermal gas column, (c) two regions of different constant densities (a density step), gravity oblique to the mesh axes | velocity, reduced and static pressure | max |u| ≤ 10⁻¹⁰ √(g L) after 1000 iterations on all four types, perturbed; static pressure equals the hydrostatic profile to round-off for (b) and (c) |
| V7.8 | Rayleigh–Bénard onset between rigid plates: growth rates of a small disturbance at Ra slightly below and above onset, critical Ra by interpolation | Ra_c | reference 1707.76 (Chandrasekhar 1961); within 0.5 % on the finest mesh; GCI(Ra_c) ≤ 0.5 % |
| V7.9 | Fully developed mixed convection in a vertical channel with walls at different temperatures, aiding and opposing flow, with flow reversal (Aung & Worku 1986) | velocity and temperature profiles | order 1.8–2.2 against the analytic solution |
| V7.10 | Thick-walled plane channel, laminar, streamwise periodic with the heat-flux mode of D109, k_s/k_f = 1 to 10⁴ | Nu at the interface, outer wall temperature | Nu = 140/17 = 8.235 (D111), outer wall temperature T_i + q t/k_s, order 1.8–2.2 |
| V7.11 | Slab with a heat source, convection and radiation on the outer face: T_w from the quartic balance | T_w, number of outer iterations | T_w to 10⁻¹² relative; ≤ 10 iterations from a cold start |
| V7.12 | Laminar free convection along an isothermal vertical plate (Ostrach 1953), Pr 0.71, Gr_x up to 10⁸ | Nu_x/(Gr_x/4)^¼, velocity and temperature profiles in similarity form | reference −θ′(0) from the similarity equations, solved by a shooting script in `tools/` and checked against Ostrach's table (0.505 at Pr 0.72); within 1 % for x more than 5 % of the plate length from the leading edge; GCI ≤ 1 % |
| V7.13 | Transient Boussinesq manufactured solution with outer correctors | time order of u and T | 1.8–2.2 for BDF2 and Crank–Nicolson at CFL 1 and 10 |
| V7.14 | Laminar heated periodic pipe with the energy equation (D111): heat-flux and wall-temperature modes | Nu | 48/11 = 4.364 and 3.657; order 1.8–2.2 |

### 14.2 Benchmarks (accurate numerical reference solutions)

| # | Case | Checks | Pass band |
|---|---|---|---|
| V7.15 | Differentially heated square cavity, air Pr 0.71, Ra 10³, 10⁴, 10⁵, 10⁶, 2D (one cell layer, D3) | mean Nu, u_max on the vertical midline and v_max on the horizontal midline with positions, Nu_max and Nu_min on the hot wall | de Vahl Davis (1983): mean Nu 1.118, 2.243, 4.519, 8.800; at Ra 10⁶ the more accurate 8.825 of Le Quéré (1991); extrapolated mean Nu within 0.5 % of the reference, velocities within 1 %, GCI(Nu) ≤ 1 % at Ra ≤ 10⁵ and ≤ 2 % at 10⁶ |
| V7.16 | Differentially heated cubic cavity, Ra 10⁴ to 10⁶ | mean Nu on the hot wall, midplane profiles | Tric, Labrosse & Betrouni (2000) and Fusegi et al. (1991), values from the papers; mean Nu within 1 % of Tric et al., GCI ≤ 2 % |
| V7.17 | Square cavity with large temperature difference, variable density, closed (p₀(t)), ε = (T_h − T_c)/(T_h + T_c) = 0.6, Ra 10⁶, Sutherland viscosity and conductivity | mean Nu on both walls, p₀/p_initial | Vierendeels, Merci & Dick (2003) and Le Quéré et al. (2005), values from the papers; Nu within 1 %, p₀ ratio within 0.1 %; total mass conserved to round-off |

### 14.3 Validation (experiments, turbulent)

Model error is measured here, so the bands are wider and are given per model
group; the numerical uncertainty must be small against them (GCI ≤ 2 %
of each compared quantity before the comparison is made).

| # | Case | Checks | Pass band |
|---|---|---|---|
| V7.18 | Turbulent natural convection along a heated vertical plate in air (Tsuji & Nagano 1988), Ra_x up to about 10¹¹, laminar–turbulent transition region excluded | Nu_x(Ra_x) in the turbulent region, mean velocity and temperature profiles | ζ-f and elliptic blending with GGDH, SST at y⁺ ≈ 1 with GGDH: Nu_x within ±10 % of the measurements; SGDH results reported to show its known underprediction, no band |
| V7.19 | Turbulent differentially heated cavity, air, Ra 1.58 × 10⁹ (Ampofo & Karayiannis 2003), 2D at the midplane | local Nu on the hot wall, mean T and v at mid-height, wall shear | near-wall models with GGDH: mean Nu within ±15 % of the measurement, profiles inside the stated experimental uncertainty in the boundary layers; transition location reported |
| V7.20 | Backward-facing step with a heated wall at constant heat flux (Vogel & Eaton 1985), Re_H 28 000, expansion ratio 1.25 | Stanton number along the heated wall (peak value and location), C_f, reattachment length | reattachment within ±10 % of the measured value; St peak within ±15 % (ζ-f, elliptic blending with the compound treatment) and ±25 % (wall functions), peak location within ±1 step height |
| V7.21 | Turbulent heated periodic pipe (D111): Pr 0.71, Re 10⁴ to 10⁶ as a study (D62) | Nu, mean temperature, turbulent heat flux | Gnielinski (1976) ±10 %; DNS of Redjem-Saad, Ould-Rouiss & Lauriat (2007) for profiles at their Re |

### 14.3a Thin walls (D189)

| # | Case | Checks | Pass band |
|---|---|---|---|
| V7.26 | Steady multi-layer thin wall (three layers, k ratios 1 to 10³, contact resistance) between a fluid and a convection boundary | wall heat flux, node and outer temperatures | equal to the series-resistance formula to round-off; equal to the resolved solid of V7.4 to its tolerance |
| V7.27 | Transient slab as a thin wall: sudden change of the ambient temperature on one face, convection on both faces, Bi 0.01 to 10 | node temperatures against the Fourier series of Carslaw & Jaeger (1959, §3.11) | time order 1.8–2.2; space order 1.8–2.2 in the node count; lumped limit at Bi 0.01 |
| V7.28 | Shell conduction in a straight fin of constant thickness with convection on both faces and an adiabatic tip, on quadrilateral, triangular and polygonal surface meshes, perturbed | T(x) against θ/θ_b = cosh(m(L − x))/cosh(mL), m² = 2h/(k t) (Bejan 2013), fin heat flow | order 1.8–2.2; fin heat flow within 0.1 % on the finest mesh |
| V7.29 | Thin wall against resolved CHT: the thick-walled periodic channel of V7.10 with t/H = 0.01 and 0.1, and a cylinder shell heated on one side with circumferential conduction (t/R = 0.01) | Nu, outer wall temperature, circumferential T | within 0.5 % of the resolved solid at t/H = 0.01; the deviation at 0.1 reported |
| V7.30 | Baffle: two laminar channels in counterflow separated by a thin plate (`mesh.baffles`) | outlet temperatures, plate temperature | within 0.5 % of the resolved plate; ledger closes |
| V7.31 | Closed box with a heater layer in a thin wall, transient, natural convection | ledger per step, entropy production | ledger to the tolerance of D89; entropy production ≥ 0 per cell and per wall node |

Film cases with thin walls and CHT are built with the film (D123 F2, D190):
Nusselt (1916) film condensation on a cooled vertical plate (film thickness
and heat transfer coefficient, the wall heated by the condensate), an
evaporating film on a thin wall heated from outside (wall temperature below
the dry case, q_wall = latent plus sensible to the ledger tolerance), and the
closed hot box of D123 with a thin wall and with a solid region (mass and
energy to round-off).

### 14.3b Several fluids and materials (D194)

| # | Case | Checks | Pass band |
|---|---|---|---|
| V7.33 | Laminar double-pipe counterflow exchanger: water in the inner tube, an oil in the annulus, the tube wall meshed (copper) and, second run, as a thin wall | outlet temperatures, heat duty, effectiveness | thin wall within 0.5 % of the meshed wall; effectiveness against ε–NTU with the fully developed Nusselt numbers of Shah & London (1978), within 2 % for a long exchanger; ledger closes to round-off between the two fluids and the wall |
| V7.34 | Three fluid parts and four materials in one case (two closed, one open), one closed part without reference | pressure levels per part, REG003 | each part reaches its own reference; the refused set-up gives REG003; results of each part equal to the part run alone (bitwise) |
| V7.35 | Tube-and-fin exchanger cell: one row of three tubes (water) with plate fins (aluminium, as shell thin walls joined to the tube walls) in air, periodic in the fin pitch | heat duty, fin efficiency, air-side pressure drop | fin efficiency against the straight-fin formula where it applies (±2 %); duty equal from the water side, the air side and the ledger to the solver tolerance; GCI ≤ 2 % |

### 14.3c Turbulent, transient and cycle-averaged conjugate heat transfer (7d5, 7d6, D204)

Added 2026-11-21 (maintainer: turbulent and transient CHT before 7b1 and
7d2). V7.6 (§14.1) is the exact check of 7d6.

| # | Case | Checks | Pass band |
|---|---|---|---|
| V7.36 | Turbulent channel (Re 2·10⁵, Pr 0.71) over a meshed solid wall whose resistance t/k_s is of the order of the fluid side's 1/h_w, outer face of the wall at a fixed temperature, k-ε (scalable) and SST (automatic) with the thermal wall function at the interface (7d5, `cht_turbulent_channel.toml`); against the same channel without the solid, whose wall carries the wall function in series with a convection condition h = k_s/t (the fluid-only path of 6f) | heat flow through the interface and the outer face; interface temperature; outlet bulk temperature; temperature drop across the wall; heat balance | interface heat flow equal to that of the outer face to 10⁻⁶; heat flow within 10⁻³ and interface temperature within 0.02 K of the fluid-only run; T_o − T_i = q t/k_s to 10⁻³; balance below 10⁻⁸ (`tests/cht.rs`). Changed from the first draft (a pipe in a steel tube against Gnielinski ±10 %): the comparison with the independent path tests the interface itself to 10⁻³, the Nusselt numbers of the wall function are checked in 6f |
| V7.37 | Heated solid wall under a laminar channel flow, heater switched on at t = 0 (the 7d1 case `cht_heated_wall` run transient with a lighter wall, ρ = 100, 7d6, `cht_heated_wall_transient.toml`) | heat balance per step: input = outflow + storage in fluid and solid; approach to the steady solution | balance to the solver tolerance in every step; the end state within 0.1 % of ΔT of the steady 7d1 result. **Met 2026-11-21 (#159):** 1.9·10⁻¹¹ in every step, end state within 1.5·10⁻⁵ of ΔT (`tests/cht.rs`) |
| V7.38 | Slab with a periodic gas temperature and heat-transfer coefficient on one face, a fixed temperature on the other (D204) | mean temperature profile after the loop; surface swing with `skin = "effusivity"` | mean profile within 0.5 % of ΔT of the fully coupled transient run in its periodic state after ≤ 3 updates; swing amplitude within 2 % of the Fourier solution |
| V7.39 | Pulsating laminar channel flow (inlet velocity and temperature periodic) over a thick conducting wall (D204) | interface heat flow, wall temperature | within 1 % of the fully coupled transient conjugate run in its periodic state; wall-clock time of both reported |
| V7.40 | Ledger of the cycle-averaged loop (D204) | heat balance of the solid per update; cycle-mean fluid balance with the change of stored heat | solid to the solver tolerance; fluid within 10⁻⁶ relative |
| V7.41 | The switched-on heated wall of V7.37 with `inertia_scale = "auto"` (D213, `cht_heated_wall_accelerated.toml`) | measured τ against the case's estimate; scaled τ one window; settled and restored within the run; end state against the steady solution; balance per step; `inertia_scaling` = (1 − f)/f × solved wall storage | τ within 10 %; f τ = window to 10⁻⁹; end state within 0.2 % of ΔT; balance below 10⁻⁹; the ratio to 10⁻¹² (`tests/cht_acceleration.rs`). **Met 2026-11-22 (#180)** |
| V7.42 | Fixed factor and restore at a time (D213) | true storage of the first step against the unaccelerated run; restore at the set time; periodic inlet: cycle mean after a 3-period ramp and after a jump against the unaccelerated periodic state | ratio between 0.9/f and 1/f; restore at the time; ramp within 0.01 K and ten times closer than the jump (`tests/cht_acceleration.rs`). **Met 2026-11-22 (#180)** |
| V7.43 | The heated wall of V7.37 with `[energy.frozen_flow]` (D118, `cht_heated_wall_frozen.toml`), a viscosity constant and following T | probe histories against the fully coupled run; refreshes; heat balance; restart in a frozen phase and in a refresh | constant viscosity: within 0.1 K over the history and 0.01 K at the end, no refresh; μ(T): within 0.15 K, at least 3 refreshes for temperature; balance below 10⁻⁹; restarts bit-identical (`tests/cht_frozen.rs`, 40 × 15 cells). **Met 2026-11-22 (#208)** |
| V7.44 | The pulsating channel with a period of 5 s warmed up by the loop with `solid = "transient"` (D204, D118, `cht_pulsating_channel_warmup.toml`) | wall temperatures at the end of every march against the cycle means of the fully coupled run; march and jump by whole periods; balance of each march; restart after a jump | within 1 % of the whole rise; balance below 10⁻⁹; restart bit-identical (`tests/cht_cycle.rs`, 40 × 10 cells). **Met 2026-11-22 (#208)** |

### 14.4 Implementation checks

| # | Check | Pass |
|---|---|---|
| V7.22 | Dual-number Jacobians of every new kernel against central differences; ⟨Av, w⟩ = ⟨v, Aᵀw⟩ for the energy matrix | 5·10⁻⁷ and round-off, as in stage 5 |
| V7.23 | CPU against GPU (Colab T4) on V7.5, V7.15 at Ra 10⁵ and V7.20; 1 and 4 threads; restart in the middle of V7.6 | bitwise equal; restart bit-identical |
| V7.24 | Views, halo plans and energy matrix of a fluid–solid mesh decomposed into 1, 2 and 4 parts (§11) | identical in global order |
| V7.25 | Runs without energy and gravity | bitwise unchanged against main on every flow case of the catalogue |
| V7.32 | Film interface slots (§9b.5): a test module sets a per-face heat sink q(T_w) with its derivative on a wall, a thin wall and a region interface | ledger pairs cancel to the bit; Newton-like convergence of the coupled wall temperature in ≤ 5 outer iterations for a steep q(T_w) (slope 10⁴ W/m²K) |

Thermal wall functions and turbulent Prandtl numbers are verified in stage 6
with temperature as a passive scalar (flat plate St, channel Nu, impinging
jet Nu); stage 7 repeats the channel case with the energy equation as a
regression check (same result to the solver tolerance when properties are
constant).


## 15. Case-file keys and diagnostics (D163)

`docs/keys.md` is generated from the code (`foehn schema --markdown`, D77)
and is regenerated when the keys exist; the planned keys are listed here.
Names follow the existing style (`[fluid]`, `[materials]`, `[regions]`,
`[boundary.<name>]` with the per-physics `[boundary.<name>.energy]` table of
§8, `[solver.<equation>]` of D127).

```toml
[energy]                            # the table turns the energy equation on; without it: flow only, bitwise as before
initial = 300.0                     # K, number or expression of x, y, z (7a1)
source = 0.0                        # q''' (W/m³), number or expression (7a1)
temperature_dependent = true        # c_p, k, μ follow T (7a1)
reference_temperature = 298.15      # K; h_s = 0 here (fixed at 298.15 K in 7a1)
viscous_heating = false             # Φ (and ρε in RANS)
pressure_work = "none"              # "none" | "thermodynamic" | "full"
thermodynamic_pressure = "constant" # "constant" (open) | "closed-volume" (p0(t), D151)
turbulent_heat_flux = "sgdh"        # "sgdh" | "ggdh"
c_theta = 0.3                       # GGDH constant (D152)
# the turbulent Prandtl number and its variable form are stage-6 keys (D135)

[[energy.sources]]
region = "heater"                   # a region of [regions]
power = 150.0                       # W over the region; or power_density (W/m³, number or expression)

[gravity]
vector = [0.0, 0.0, -9.81]          # m/s²
reference_point = [0.0, 0.0, 0.0]   # x_ref of the reduced pressure

[buoyancy]
model = "boussinesq"                # "none" | "boussinesq" | "variable-density"
reference_temperature = 298.15      # T_ref; default energy.reference_temperature
expansion_coefficient = 3.4e-3      # β (1/K); default from the material (ideal gas 1/T_ref)
reference_density = "auto"          # ρ_ref; "auto" = ρ(T_ref, p0)
production = "auto"                 # G_b form: "auto" (as turbulent_heat_flux) | "sgdh" | "ggdh" | "none"
c_eps3 = "henkes"                   # "henkes" (tanh form) | a number
omega_source = "transformed"        # "transformed" | "none" (SST)

[regions.fin]
kind = "solid"
zones = ["fin"]
material = "aluminium"

[materials.my-steel]
base = "steel"
conductivity = { form = "orthotropic", values = [45.0, 45.0, 15.0], axes = [[1,0,0],[0,1,0]] }
specific_heat = { form = "table", file = "steel_cp.csv" }

[boundary.hot]
type = "wall"
[boundary.hot.energy]
temperature = 320.0

[boundary.outer]
type = "wall"
[boundary.outer.energy]                 # convection and radiation: the keys of both
heat_transfer_coefficient = 10.0
ambient_temperature = 293.15
emissivity = 0.8
thin_wall = { thickness = 0.002, material = "glass" }

[boundary.inlet]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.inlet.energy]
temperature = 300.0

[boundary.outlet]
type = "static-pressure-outlet"
pressure = 0.0                      # reduced pressure p′ when gravity is set
pressure_kind = "reduced"           # "reduced" | "static"
[boundary.outlet.energy]
backflow_temperature = 300.0

[interfaces.fin-air]                # optional: a face set between two regions
regions = ["fin", "air"]
contact_resistance = 2.0e-4         # m²K/W; or thin_wall = { ... }

[solver.energy]
scheme = "linear-upwind"
limiter = "venkatakrishnan"
relaxation = 0.95                   # fluid rows; solid rows use relaxation_solid
relaxation_solid = 1.0
tolerance = 1e-3
max_iter = 50
preconditioner = "amg"
amg_rebuild = 0.1
position = "after-flow"             # "after-flow" | "before-momentum"

[output]
derived = ["static_pressure", "wall_heat_flux", "nusselt", "y_plus_thermal",
           "heat_transfer_coefficient", "entropy_production", "buoyancy_production"]
reference = { temperature = 293.15, length = 0.1, conductivity = 0.026 }
```

Run-time changes (D40, D128): thermal boundary values, source powers and
the ambient temperature may change during a run; models and regions may not.
New log columns (D125 kinds): `boundary.<name>.heat_flow`,
`interface.<name>.heat_flow`, `interface.<name>.temperature`,
`cell.max.temperature`, `cell.min.temperature`, `run.p0`.
New field-registry entries (D115): `temperature`, `enthalpy`, `density`,
`conductivity`, `specific_heat` (cell), face heat flux (boundary face),
restart with h and T at every time level.

Diagnostic codes (D78): THM001–THM009 and BUO001–BUO006 (§8), CHT001 error
(interface between regions not in one energy domain), CHT002 warning (gaps or
overlaps of an imprinted interface above the tolerance of D58), CHT003 info
(per-interface ratio of conductivities and the expected AMG behaviour),
CHT004 error (contact resistance on a face that is not a region interface),
ENE001 error (temperature recovery failed, D92), ENE002 warning (damped
temperature update, count), ENE003 warning (temperature outside a property
range, with D96 clamp), ENE004 warning (Mach number above 0.3), ENE005 info
(transient buoyancy without outer iterations), ENE006 error (an energy
setting the version does not have yet, PR 7a1), MAT006–MAT009 (§8), TWL001–TWL005 (§9a.5), CYC001–CYC004 (§9d.4). The
registry and `foehn explain <code>` (D78) carry the guidance text.


## 16. Implementation plan (D164)

The working plan of stage 7 (finer PR cut, lanes, status per PR) is on
the [stage 7 plan page](stage7-plan.md).

Six pull requests, each with its tests, its verification page section
(`docs/verification/stage7.md`), gallery examples, case files and recipes
on a new topic page `energy.md` in `docs/gallery/`, `docs/recipes/` and
`docs/cases/`, and a changelog fragment:

| PR | Content | Crates and files | Cases |
|---|---|---|---|
| 7a | Energy core: density, μ, k, c_p as cell fields in the flux kernels (finding F1), thermal property forms and database records (D102), energy module on the fluid domain with the δT form, SGDH, thermal boundary conditions without CHT, THM codes, thermal ledger, D90 parts 1–3 | `foehncraft-kernels::energy` (new), `flow.rs` (ρ field); `foehncraft-flow::energy` (new), `bc.rs` (thermal sub-table), `materials.rs`, `balance.rs`, `diag.rs`, `case.rs`; `data/materials/` | V7.1–V7.3, V7.11, V7.14, V7.25 |
| 7t | Thin walls (D189, §9a): layers with heat capacity and sources (1D), shell conduction on the surface mesh (D124 W1, W2 moved here), baffles from face zones, edge conditions, TWL codes; the per-face film slots with their linearisation on walls, thin walls and interfaces (D190, §9b.5) | `foehncraft-mesh` (surface mesh, `baffles`), `foehncraft-kernels::energy` (series and edge kernels), `foehncraft-flow::energy`, `thinwall.rs` (new), `balance.rs` | V7.26–V7.32 |
| 7b | Gravity, Boussinesq, balanced-force interpolation (D154), wall pressure, reduced pressure at boundaries, BUO codes | `foehncraft-kernels::flow` (`face_mass_flux`, momentum cell rows), `foehncraft-flow::buoyancy` (new), `bc.rs`, `derived.rs` | V7.7, V7.8, V7.9, V7.12, V7.13, V7.15, V7.16 |
| 7c | Low-Mach variable density, p₀(t), pressure work | `foehncraft-flow::simplec`, `transient`, `energy` | V7.7 (b), V7.17 |
| 7d | Conjugate heat transfer: domain views (D112), region interfaces, one-sided gradients, contact resistance and thin walls, solids with orthotropic k (D116 tensor), energy sources, D58 imprinting, METIS weights and per-view halo plans (D122 (1)), CHT codes | `foehncraft-mesh` (`view.rs`, `imprint.rs` new, `partition.rs`, `decompose.rs`), `foehncraft-fv::amg` (views), `foehncraft-flow::energy`; `foehncraft-parallel` (halo plans per view) | V7.4–V7.6, V7.10, V7.24; D58 cases of the old stage 7 (Poiseuille and T-junction across a non-conformal interface) |
| 7e | Turbulent heat transfer: GGDH, G_b in k, ε, ω, ζ/φ, C_ε3, Robin and radiation linearisation with wall functions; the `kays-crawford-hp` and `wall-profile` Pr_t closures (D172) and the wall treatment with variable properties, semi-local y⁺_sl and the integrated thermal wall function (D173, step 21; checks H-0 to H-2 of the step 21 page) | `foehncraft-kernels::energy`, stage-6 turbulence modules (hook) | V7.18–V7.21 |
| 7f | GPU parity of all stage-7 kernels, Jacobian tests, verification manual, release v0.7 | `foehncraft-cuda`, `gpu/ptx-kernels`, `tools/notebooks/stage7_colab.ipynb`, `tools/` | V7.22, V7.23 |

7t added 2026-11-20 (D189, D190): its 1D part needs only 7a and can come before
full CHT (7d); its shell part and the interface slots come after 7d's domain
views. 7a can start once stage 6 has merged its field and equation registries (D115,
D116) and the source hook; 7b and 7c can run beside 7d. The Python package
(D61) gets the new keys with each PR (case round trip bitwise, as in stage 5).


### 16.1 The finer cut (decided 2026-11-20)

The maintainer chose on 2026-11-20 to build stage 7 in 19 small pull
requests instead of the six above (which remain the grouping). Each brings
its GPU kernels with a CPU = GPU check (the GPU part of 7f moves into every
PR; 7f keeps the Colab notebook, the Jacobian sweep and the release), its
verification section, and, for physics, a gallery entry with 1–3 field
images, case files and a recipe. After 7a1 the buoyancy lane (B) and the
CHT lane (C) run in parallel.

| Lane | PR | Content | Checks |
|---|---|---|---|
| 0 | 7.0a | Ledger as data (D119): named rows per region; mass and kinetic energy moved onto it (stage 6 debt, finding F6) | V7.25 bitwise, ledger numbers identical |
| 0 | 7.0b | ρ, μ as cell and face fields in the flux, momentum and scalar kernels (F1), CPU and GPU | V7.25 bitwise |
| 0 | 7.0c | Thermal material data: c_p and k of the database fluids, NASA-7/9 forms and importer, water by IAPWS (below), solids starter set, orthotropic k, glycol–water mixtures, user database files, MAT006–MAT009 | unit tests against the sources |
| A | 7a1 | Energy module on the fluid, constant ρ: δT form, T(h), k(T), c_p(T), SGDH, thermal BCs, viscous heating, ENE and THM codes, thermal ledger | V7.1–V7.3, channel DNS regression |
| A | 7a2 | Wall-temperature periodic mode (D109 (2)) | V7.14 |
| A | 7a3 | Radiation to ambient, convection–radiation, backflow temperature | V7.11 |
| A | 7t1 | Thin walls 1D (D189): layers, heat capacity, sources, contact resistance | V7.26, V7.27 |
| B | 7b1 | Gravity, Boussinesq, reduced pressure, balanced force in SIMPLEC **and the coupled algorithm** (D176), wall pressure, BUO codes | V7.7 (a) |
| B | 7b2 | Laminar natural convection | V7.9, V7.12, V7.15 |
| B | 7b3 | Transient buoyancy | V7.8, V7.13, V7.16 |
| B | 7c | Low-Mach variable density, p₀(t), pressure work | V7.7 (b), (c), V7.17 |
| C | 7d1 | Domain views, solid regions, interface faces, one-sided gradients, AMG over the energy view | V7.4–V7.6 |
| C | 7d2 | Contact resistance, orthotropic k, energy sources | V7.4 with resistance, V7.10 |
| C | 7d3 | General static non-conformal interfaces by imprinting (D58 beyond surfaces of revolution, which `imprint.rs` already handles) | D58 cases, V7.4 on non-matching layers |
| C | 7d4 | METIS with two constraints, per-view halo plans | V7.24 |
| C | 7t2 | Shell conduction and baffles (D189), surface mesh W1 and wall-face halo W2, film slots (D190) | V7.28–V7.32 |
| C | 7d5 | Turbulent CHT: thermal wall function (D135) on the fluid side of every interface in series with the solid half cell (§9.4), turbulent conductivity c_p μ_t/Pr_t in the energy view (done, #158) | V7.36 |
| C | 7d6 | Transient CHT: time terms of the fluid and the solids in the one energy matrix, heat storage in the balance, output and probes on the full mesh (done, #159) | V7.6, V7.37 |
| C | later | Cycle-averaged coupling of periodic flows (D204, §9d), after 7d6; core implemented (#190), skin swing, resolved cycles and restart (#202), transient solid phase and the frozen-flow march of D118 (#208) | V7.38–V7.40, V7.43, V7.44 |
| D | 7e1 | `buoyancy_production` hook in k, ε, ω, ζ/φ (not delivered by stage 6), G_b with SGDH, C_ε3 | V7.25, V7.18 (SGDH, reported) |
| D | 7e2 | GGDH, G_b with GGDH | V7.18, V7.19, V7.21 |
| D | 7e3 | Heated backward-facing step | V7.20 |
| D | 7e4 | D172 Pr_t closures, D173 variable-property wall treatment | H-0 to H-2 of step 21 |
| – | 7f | Jacobians, Colab notebook, restart check, verification manual, release v0.7 | V7.22, V7.23 |

**7a1 done 2026-11-21 (#154):** the energy equation of the fluid with
c_p(T), k(T), μ(T), viscous heating, the thermal walls (`[boundary.<name>.energy]`
tables in the style of the passive scalars, now the form of §8 and §15),
THM001–THM005, THM009–THM011, ENE006, and the heat balance; V7.1 on hexahedra with a
Dirichlet boundary (the Robin condition by the exact series resistance in
`tests/energy.rs`), V7.2 with both walls isothermal. Moved to a later PR:
THM006–THM008, V7.3, the thermal account of the ledger in `result.json`,
the channel DNS regression of the energy equation, and the property
refresh on the device.

**7d1 done 2026-11-21 (#156):** solid regions (`kind = "solid"` with a
solid material) in the energy equation; the flow on a fluid view (a
submesh whose fluid–solid faces are the walls `interface-<region>`), the
energy on all cells with the series conductance at faces between
materials, the interface heat flow and temperature in the report, and
`zones` boxes in the box generator. V7.4 to round-off at ratios
1 : 100 : 10⁴, and a heated wall under a channel flow. Moved to a later
PR: one-sided gradients (§9.2, exact on orthogonal meshes until then),
V7.5 and V7.6, the AMG with strength on the scaled matrix over the energy
view, transient, turbulent and restarted conjugate runs, passive scalars
with solids, periodic meshes with solids, and maps on the device.

**7d5 done 2026-11-21 (#158):** turbulent conjugate heat transfer. The
fluid side of every interface face is the thermal wall function of the
flow's turbulence model in series with the solid half cell (§9.4), the
fluid conducts with k + c_p μ_t/Pr_t, and the fluid's outer walls keep
their wall functions in the energy view; μ(T) goes back to the laminar
viscosity of the turbulence model. V7.36 against the fluid-only path for
k-ε and SST, and a finned heat sink in a turbulent air channel (gallery
158). Found on the way: a `[[schedule]]` entry dropped the energy and
scalar sub-tables of the boundaries (fixed), and a uniform start velocity
against a solid's front face can make the SIMPLEC pressure-correction
matrix indefinite in the first iterations (a flow-solver item, not of the
conjugate coupling; the gallery case starts from rest). Moved to a later
PR: the turbulence fields on the whole mesh in the output, y⁺ and T⁺ of
the interfaces as derived fields (D72).

**7d6 done 2026-11-21 (#159):** transient conjugate heat transfer. The
fluid and the solid cells share the one energy matrix with the time terms
ρ c_p V/Δt of their own materials (BDF2, also under Crank–Nicolson as for
every scalar); the temperature's time levels in the energy view follow
those of the flow (step start, rejected steps of D65, the rewind of
D63). The balance has the heat stored per region; `heat.csv` writes it
per step with the heat flow of each interface, the run record carries it,
and the report gives the largest relative imbalance of a step. Probes
give the temperature, on the whole mesh in a conjugate case, and
`fields_<step>` are written there. The expression language gained `erf`
and `erfc` (in-house, D93 (6)), so `[energy] exact` carries the erf
solutions. Results: V7.6 with the interface at T_i to 10⁻¹² K at 500 s
(within 3·10⁻³ K, 0.003 % of ΔT, after 10 steps), the copper–air and
equal-body variants, L2 error 0.0784/0.0193/0.0048 K on 100/200/400
cells (orders 2.02, 2.00) and time orders 2.05, 2.02 (BDF2 and
Crank–Nicolson alike); V7.37 with the balance closed to 1.9·10⁻¹¹ in
every one of 600 steps and the end state within 1.5·10⁻⁵ of ΔT of the
steady 7d1 result. Found on the way: transient sweeps did not check the
scalars against their exact solutions (now as in steady runs), and a time
sweep of a conjugate case would have read the temperature on the fluid
mesh only (fixed before it could be hit). Still open: restart of a
conjugate run, solids with a sliding mesh, the turbulence fields on the
whole mesh.

**7d2 done 2026-11-21 (#160):** contact resistances and thin walls (t/k,
no heat capacity) at `[interfaces.<name>]` between fluid and solid or two
solids, orthotropic solid conductivity with axes per material or region
(the matrix takes n̂·K·n̂ of the face, the tangential rest is deferred, as
are the cross terms c = −(K n̂ − k_n n̂)·∇T of an orthotropic side in the
interface series of §9.1), one-sided gradients of §9.2 iterated three
times per sweep from the last values, `[[energy.sources]]` per region,
CHT004 and ENE008. Results: the patch test with a laminate turned by 30°,
a fluid at rest and a block behind R″_c = 0.01 m²K/W on a mesh perturbed
by 25 % (the V7.4 variant of §9.1) exact to 3·10⁻¹¹ K and 200 W/m² to 10⁻⁸
through every wall and interface; V7.4 with a contact resistance and a
thin wall equal to the series formula to 10⁻⁹; a chip on a thermal pad on
a circuit board with copper planes (gallery 160). The aggregation AMG of
the pressure as the preconditioner of the conjugate energy matrix was
tried and left out: on the heated wall of 7d1 it took 206 BiCGStab
iterations against 233 with Jacobi at four times the cost per iteration,
and on the circuit board it did not reach the tolerance; it needs the
strength on the scaled matrix of D158 (§10.2) first. A loose energy solve
that converges with an open heat balance (ENE007) is a matter of
`max_iter`: the circuit board closes to 7·10⁻¹³ with 500 Jacobi
iterations. Still open: V7.10 (periodic meshes with solids, ENE006), the
entropy production of a contact resistance in the second-law check, the
cross terms at convection boundaries (taken with the weight of the face
value), frames that follow the geometry (L48).

**Decided 2026-11-21 (maintainer):** turbulent and transient conjugate heat
transfer (7d5, 7d6) come before 7b1 and 7d2, as the more industry-relevant
steps; buoyancy (lane B) comes much later. The cycle-averaged coupling of
periodic flows (D204, §9d) is recorded as a requirement for after 7d6.

Other decisions of 2026-11-20:

- v0.6 is tagged before stage 7 code merges.
- The two stage 6 items that were not delivered (ledger as data, the
  buoyancy hook) are done in stage 7 (7.0a, 7e1).
- Water: IAPWS-95 (with the IAPWS formulations of μ and k) is evaluated in
  the database code into monotone tables (D96) at the operating pressure
  when the case is resolved; the kernel reads tables. Every other form of
  D95 (polynomial, Sutherland, power law, NASA-7/9, tables, expressions)
  stays available for every property, water included.
- Solids starter set: aluminium (pure and 6061), copper, carbon steel,
  stainless steel 304, titanium, fused silica glass, silicon, PTFE, from open
  sources; glycol–water at any fraction. User-defined materials with
  thermal data and solids in the case file (`[materials.<name>]`) and in user
  database files (`[materials] include = ["my_materials.toml"]`, the TOML
  format of D103), so users extend the fluid and solid database without
  changing the built-in one.
- Reference data: only published numbers in the repository; tables that may
  not be redistributed are fetched by a script into a git-ignored folder
  (as for the channel DNS); the maintainer supplies what is not openly
  available.
- Turbulent bands (V7.18–V7.21): ζ-f₀ and elliptic-blending k-ε with GGDH
  must pass; SST and k-ε are reported without a band (maintainer,
  2026-11-20).
- §18 point 3 and 5 defaults (reduced pressure at pressure boundaries, C_θ
  = 0.3, transformed ω source) as recommended (maintainer, 2026-11-20); the
  validation cases may change C_θ and the ω source.
- Reduced and static pressure get dedicated 2D examples with recipes in
  7b1 (maintainer, 2026-11-20): a tall 2D channel with gravity, run once
  with `pressure_kind = "reduced"` and once with `"static"` at the outlet
  (same flow, the static pressure field shows the hydrostatic part), a 2D
  box at rest where a static outlet value of 0 Pa without the hydrostatic
  part drives a false flow (the reason for the default), and a recipe
  "which pressure do my boundary values mean" with the field images of p′
  and p side by side.

## 17. Findings in the code today

Read on main at b15a488:

- **F1:** the flux and momentum kernels take the density as one scalar
  (`FluxInputs::rho`, `MomentumInputs::rho` in `foehncraft-kernels/src/flow.rs`);
  D102 planned the density field with stage 5 part 6a, which shipped constant
  density only ("the density must be constant in the domain until L1"). PR 7a
  turns it into cell and face fields.  **Done 2026-11-20 in PR 7.0b:** the kernels read ρ per
  cell (momentum time terms, scalar and pseudo-time rows, QCR stress, wall
  pressure gradient) and per face (mass flux, pressure-correction
  coefficient, boundary velocity and total pressure, turbulent face
  viscosities, wall μ_eff) from the mesh block, on the CPU and the GPU; the
  device fills both with the constant density and takes new fields with
  `FlowDevice::set_density`. μ was already a face field (D94). Still scalar
  until 7a1: the kinematic ν of the wall functions and turbulence sources,
  the Courant factor and the pseudo-time Schur diagonal of the coupled
  solver.
- **F2:** the momentum cell rows use (b − ∇p_P) V with the least-squares ∇p
  (`flow.rs`, the source term of the momentum cell kernel) and `face_mass_flux`
  uses the pressure jump only; with gravity this gives the spurious wall
  velocities noted in docs/verification/stage5.md. PR 7b (D154).
- **F3:** wall pressure is zero-gradient; with gravity it must follow the
  reduced jump (D154 (4)).
- **F4:** the flow solves on every cell and refuses non-fluid regions
  (REG001); PR 7d needs the flow view on fluid cells and the energy view on
  all cells (D112). **Done 2026-11-21 in PR 7d1.**
- **F5:** METIS gets one constraint and no weights (finding F6 of the
  multi-physics plan, D122 (1)); PR 7d.
- **F6:** `balance.rs` has fixed mass and kinetic-energy structs; the ledger
  as data (D119) is planned for stage 6 and is a prerequisite of 7a. **Done 2026-11-20 in PR 7.0a.**


## 18. Points for the maintainer

1. **O18, stage order and milestone 1** (§1), **decided 2026-11-17: as proposed**: energy, buoyancy and CHT as
   stage 7, the old stage 7 as stage 8 with v1.0; D58 moves into stage 7.
   Alternative: v1.0 without heat transfer and energy after it.
2. **D92 (3) replaced** by the temperature-correction form (§3.2); the
   conserved variable stays h. My recommendation; the reasons are in §3.2.
3. **Reduced pressure at pressure boundaries** when gravity is set (§6.1),
   with `pressure_kind = "static"` as the alternative per boundary.
4. **L27** (electric conduction and Joule heating) stays after milestone 1,
   first in the order.
5. **GGDH constant C_θ = 0.3** and the ω-equation buoyancy source as
   defaults, both decided by the validation cases.


6. **Thin walls and the film (D189, D190), added 2026-11-20**, **decided
   2026-11-20: TW1–TW5 as recommended, TW6 later when the wall film exists**
   (O24 closed):
   - TW1 *Baffles from meshes:* split named internal face zones into boundary
     pairs in every mesh reader and the mesher (`mesh.baffles`). Recommended.
   - TW2 *Thin walls in the matrix or beside it:* monolithic (extra rows of
     the energy matrix, recommended, as D158) or a separate solve per thin
     wall coupled per outer iteration (cheaper per iteration, slower and
     less robust at high wall conductance).
   - TW3 *Shell edges on solid regions* (a thin fin on a meshed base): join
     the edge to the solid cells now, or refuse (TWL004) in the first version.
     Recommended: refuse first, join in a follow-up. **Changed 2026-11-20 by
     D194:** fins of a tube-and-fin exchanger are shell thin walls whose
     edges sit on the tube walls (meshed or thin), so the edge joint is part
     of 7t2 and TWL004 is dropped.
   - TW4 *Thin walls with radiation between surfaces* (L15): the thin-wall
     outer node takes the radiation flux slot like any wall face. Recommended,
     no extra work now.
   - TW5 *Film evaporation needs a vapour species in the gas.* Evaporation
     and condensation move mass into and out of the gas, so even a
     one-component film needs the vapour mass fraction in the gas (species
     transport, L2) and variable density (D151). D123 F2 lists L1 only;
     recommended: F2 needs L1 and the single-vapour part of L2.
   - TW6 *When the film comes:* it stays with the sprays (L12) after v1.0;
     stage 7 builds only the interface (§9b.5). Alternative: a film without
     sprays (wetted walls, falling films, condensation) earlier, with L2.

7. **Several fluids and materials (D194), added 2026-11-20**, decided:
   - MF1 *Models per fluid region:* turbulence model, buoyancy model and
     Pr_t set per fluid region. **Decided 2026-11-20 (maintainer): yes, per
     fluid region**, case-wide values stay the default.

## 19. References

- W. Aung, G. Worku, *Theory of fully developed, combined convection
  including flow reversal*, J. Heat Transfer 108, 485–488 (1986)
- F. Ampofo, T. G. Karayiannis, *Experimental benchmark data for turbulent
  natural convection in an air filled square cavity*, Int. J. Heat Mass
  Transfer 46, 3551–3572 (2003)
- A. Bejan, *Convection Heat Transfer*, 4th ed., Wiley (2013)
- H. S. Carslaw, J. C. Jaeger, *Conduction of Heat in Solids*, 2nd ed.,
  Oxford University Press (1959)
- I. B. Celik et al., *Procedure for estimation and reporting of uncertainty
  due to discretization in CFD applications*, J. Fluids Eng. 130, 078001
  (2008)
- S. Chandrasekhar, *Hydrodynamic and Hydromagnetic Stability*, Oxford
  University Press (1961)
- S. K. Choi, S. O. Kim, C. H. Lee, H. K. Choi, *Use of the momentum
  interpolation method for flows with a large body force*, Numer. Heat
  Transfer B 43, 267–287 (2003)
- B. J. Daly, F. H. Harlow, *Transport equations in turbulence*, Phys.
  Fluids 13, 2634–2649 (1970)
- J. Degroote, K.-J. Bathe, J. Vierendeels, *Performance of a new
  partitioned procedure versus a monolithic procedure in fluid–structure
  interaction*, Computers & Structures 87, 793–801 (2009)
- G. de Vahl Davis, *Natural convection of air in a square cavity: a bench
  mark numerical solution*, Int. J. Numer. Methods Fluids 3, 249–264 (1983)
- M. M. Francois, S. J. Cummins, E. D. Dendy, D. B. Kothe, J. M. Sicilian,
  M. W. Williams, *A balanced-force algorithm for continuous and sharp
  interfacial surface tension models within a volume tracking framework*,
  J. Comput. Phys. 213, 141–173 (2006)
- T. Fusegi, J. M. Hyun, K. Kuwahara, B. Farouk, *A numerical study of
  three-dimensional natural convection in a differentially heated cubical
  enclosure*, Int. J. Heat Mass Transfer 34, 1543–1557 (1991)
- M. B. Giles, *Stability analysis of numerical interface conditions in
  fluid–structure thermal analysis*, Int. J. Numer. Methods Fluids 25,
  421–436 (1997)
- V. Gnielinski, *New equations for heat and mass transfer in turbulent
  pipe and channel flow*, Int. Chem. Eng. 16, 359–368 (1976)
- D. D. Gray, A. Giorgini, *The validity of the Boussinesq approximation for
  liquids and gases*, Int. J. Heat Mass Transfer 19, 545–551 (1976)
- C. Y. Gu, *Computation of flows with large body forces*, in: Numerical
  Methods in Laminar and Turbulent Flow 7, Pineridge Press, 294–305 (1991)
- K. Hanjalić, *Achievements and limitations in modelling and computation
  of buoyant turbulent flows and heat transfer*, Proc. 10th Int. Heat
  Transfer Conf., Brighton, vol. 1, 1–18 (1994)
- K. Hanjalić, S. Kenjereš, F. Durst, *Natural convection in partitioned
  two-dimensional enclosures at higher Rayleigh numbers*, Int. J. Heat Mass
  Transfer 39, 1407–1427 (1996)
- R. A. W. M. Henkes, F. F. van der Vlugt, C. J. Hoogendoorn,
  *Natural-convection flow in a square cavity calculated with
  low-Reynolds-number turbulence models*, Int. J. Heat Mass Transfer 34,
  377–388 (1991)
- N. Z. Ince, B. E. Launder, *On the computation of buoyancy-driven
  turbulent flows in rectangular enclosures*, Int. J. Heat Fluid Flow 10,
  110–117 (1989)
- B. A. Kader, *Temperature and concentration profiles in fully turbulent
  boundary layers*, Int. J. Heat Mass Transfer 24, 1541–1544 (1981)
- W. M. Kays, *Turbulent Prandtl number – where are we?*, J. Heat Transfer
  116, 284–295 (1994)
- S. Kenjereš, S. B. Gunarjo, K. Hanjalić, *Contribution to elliptic
  relaxation modelling of turbulent natural and mixed convection*, Int. J.
  Heat Fluid Flow 26, 569–586 (2005)
- F. Kock, H. Herwig, *Local entropy production in turbulent shear flows: a
  high-Reynolds number model with wall functions*, Int. J. Heat Mass
  Transfer 47, 2205–2215 (2004)
- U. Küttler, W. A. Wall, *Fixed-point fluid–structure interaction solvers
  with dynamic relaxation*, Comput. Mech. 43, 61–72 (2008)
- P. Le Quéré, *Accurate solutions to the square thermally driven cavity at
  high Rayleigh number*, Computers & Fluids 20, 29–41 (1991)
- P. Le Quéré, C. Weisman, H. Paillère, J. Vierendeels, E. Dick, R. Becker,
  M. Braack, J. Locke, *Modelling of natural convection flows with large
  temperature differences: a benchmark problem for low Mach number solvers.
  Part 1. Reference solutions*, ESAIM: M2AN 39, 609–616 (2005)
- J. Mencinger, I. Žun, *On the finite volume discretization of
  discontinuous body force field on collocated grids: application to VOF
  method*, J. Comput. Phys. 221, 524–538 (2007)
- W. Nusselt, *Die Oberflächenkondensation des Wasserdampfes*, Z. VDI 60,
  541–546 and 569–575 (1916)
- S. Ostrach, *An analysis of laminar free-convection flow and heat transfer
  about a flat plate parallel to the direction of the generating body
  force*, NACA Report 1111 (1953)
- S. Paolucci, *On the filtering of sound from the Navier–Stokes equations*,
  Sandia Report SAND82-8257 (1982)
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere
  (1980)
- S. Redjem-Saad, M. Ould-Rouiss, G. Lauriat, *Direct numerical simulation
  of turbulent heat transfer in pipe flows: effect of Prandtl number*, Int.
  J. Heat Fluid Flow 28, 847–861 (2007)
- R. G. Rehm, H. R. Baum, *The equations of motion for thermally driven,
  buoyant flows*, J. Res. Natl. Bur. Stand. 83, 297–308 (1978)
- R. K. Shah, A. L. London, *Laminar Flow Forced Convection in Ducts*,
  Academic Press (1978)
- E. Tric, G. Labrosse, M. Betrouni, *A first incursion into the 3D
  structure of natural convection of air in a differentially heated cubic
  cavity, from accurate numerical solutions*, Int. J. Heat Mass Transfer 43,
  4043–4056 (2000)
- P. Vaněk, J. Mandel, M. Brezina, *Algebraic multigrid by smoothed
  aggregation for second and fourth order elliptic problems*, Computing 56,
  179–196 (1996)
- J. Vierendeels, B. Merci, E. Dick, *Benchmark solutions for the natural
  convective heat transfer problem in a square cavity with large horizontal
  temperature differences*, Int. J. Numer. Methods Heat Fluid Flow 13,
  1057–1078 (2003)
- J. C. Vogel, J. K. Eaton, *Combined heat transfer and fluid dynamic
  measurements downstream of a backward-facing step*, J. Heat Transfer 107,
  922–929 (1985)
- V. R. Voller, C. R. Swaminathan, *General source-based method for
  solidification phase change*, Numer. Heat Transfer B 19, 175–189 (1991)
- T. Tsuji, Y. Nagano, *Characteristics of a turbulent natural convection
  boundary layer along a vertical flat plate*, Int. J. Heat Mass Transfer
  31, 1723–1734 (1988)
- D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries
  (2006)
- S. Zhang, X. Zhao, S. Bayyuk, *Generalized formulations for the
  Rhie–Chow interpolation*, J. Comput. Phys. 258, 880–914 (2014)
