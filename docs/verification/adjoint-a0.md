# Adjoint stage A0: measurement note (L17, D220, D221)

Stage A0 of the [adjoint plan](../requirements/step27-adjoint-solver.md)
(section 4): no solver change; the steady residual as one callable, the
state vector, option A of the reverse sweep measured on one face kernel,
a time-boxed check of option B, and the answer to question Q2. Pull
request #271, 2026-11-24. Gallery
[271](../gallery/verification.md#adjoint-a0), recipe
[271](../recipes/verification.md#adjoint-a0).

## 1. The residual as one callable

`Flow::steady_residual` (`crates/foehncraft-flow/src/residual.rs`) returns
R(q) of the steady discrete equations at the current fields:

| Part | Per | Definition |
|---|---|---|
| R_u | cell and component | b − A u, the unrelaxed momentum rows with the current face fluxes, the pressure gradient on the right-hand side |
| R_c | cell | Σ_f ± ṁ_f, the net outflow of the current face fluxes |
| R_ṁ | face | (ṁ − ṁ*(u, p, ṁ))/α_u, ṁ* the Rhie–Chow flux with the relaxation term of the current fields |

**State vector: the face fluxes are in q** (the question of section 0 of
the plan). The Rhie–Chow flux reads D = α V/a_P, and a_P holds the
convection coefficients of ṁ, so ṁ is defined only implicitly; q = (u, p,
ṁ). Divided by α_u, R_ṁ = (ṁ − ρ ū·S) + ρ (V/a_P)‾ |S|²/(d·S) [jump] holds no
relaxation factor.

Checks (`cargo test --release -p foehncraft-flow --lib steady_residual`,
perturbed-polygon cavity at Re 100):

- after 5 SIMPLEC iterations ‖R_u‖ = 2.4 × 10⁻², ‖R_ṁ‖ = 5.2 × 10⁻²
  (‖R_c‖ = 6.5 × 10⁻¹², since the pressure correction makes the fluxes
  conservative in every iteration); SIMPLEC after 600 iterations has all
  three at 1.4 to 1.6 × 10⁻¹⁶, the coupled solver after 150 iterations
  ‖R_u‖ = 1.0 × 10⁻⁹, ‖R_c‖ = 1.1 × 10⁻¹⁶, ‖R_ṁ‖ = 9.4 × 10⁻¹⁰: the callable
  is the residual both algorithms solve (decision Q4);
- the same state gives the same R for α_u = 0.8 and 0.5 (relative change
  6.8 × 10⁻¹⁴ of R_ṁ, 1.9 × 10⁻¹³ of R_u);
- a call between two outer iterations leaves the run bitwise unchanged
  (both algorithms; boundary values, gradients and limiters are put
  back).

Scope: steady laminar flow of constant density on one rank; turbulence,
transported scalars, variable density, time steps and several ranks are
refused with CASE005.

## 2. Option A: local dual Jacobians and the inverse-stencil gather

`foehncraft_kernels::reverse`: the inputs of one face of
`flow::momentum_face` are copied into a two-cell local problem and the
same kernel body runs on `Dual<26>` (ṁ_f, μ_f, u_P, u_N or u_b, ∇u_P,
∇u_N). x̄ = Jᵀ ȳ per face, then the cell parts are summed over the faces
of each cell in stored order (`Topology::cell_faces`): no atomics, every
value written by one index.

Verification (`cargo test --release -p foehncraft-flow --test adjoint reverse`,
mixed mesh of hexahedra, prisms, tetrahedra and polyhedra with every
boundary kind):

| Check (D222 rung) | Result |
|---|---|
| local copy against the primal kernel, every face | values equal |
| rung 1, dot product ⟨ȳ, J ẋ⟩ = ⟨Jᵀ ȳ, ẋ⟩ over the mesh, forward side from the primal kernel on `Dual<1>` | 5.8 × 10⁻¹⁷ relative |
| rung 5 at kernel level, Σ ȳ·y against central differences, steps 10⁻³ to 10⁻⁷ | 1.1 × 10⁻¹² at 10⁻³, 9.7 × 10⁻⁹ at 10⁻⁷ (round-off) |
| D31, 1 against 4 threads | bitwise equal |

Cost on 141 600 mixed cells, 322 230 faces, one thread
(`FOEHN_A0_CELLS=40,40,30`, cloud container):

| Step | ns per face | × primal face kernel |
|---|---|---|
| primal face kernel | 49 | 1 |
| local Jacobian on `Dual<26>` | 1 626 | 33 |
| pullback with a stored block | 201 | 4.1 |
| cell gather | 27 | 0.5 |
| **stored mode, per product** | | **4.6** |
| **recompute mode, per product** | | **38** |

Stored blocks take 1 040 B per face (335 MB here, about 2.4 kB per cell:
some 2.4 GB at 1 M cells), as the plan's rule of thumb said: stored below
about 1 M cells, recompute above.

## 3. Option B: compiler reverse mode, time-boxed

Enzyme through rustc's `#[autodiff_reverse]` on rustc 1.102.0-nightly
(2026-11-23, LLVM 23.1.6), applied to the unchanged `momentum_face` body on
one interior face (`tools/spikes/l17-enzyme-momentum-face`, not in the
workspace):

- **builds**, with two workarounds: the rustup component `enzyme`, a
  library path to rustc's LLVM, and the flag through `cargo rustc` (cargo's
  target probe fails with `RUSTFLAGS=-Zautodiff=Enable`);
- **exact**: Jᵀ ȳ equal to option A to 4.5 × 10⁻¹⁶, values bitwise equal;
- **fast**: 77 ns per reverse evaluation against 22 ns for the primal of
  the same local copy (3.5 ×), against 1 244 ns for option A recomputing
  its block on the same face (16 × slower than B);
- **deterministic by construction**: B differentiates the face-local
  function only; the gather over the inverse stencil stays ours, so no
  atomics appear;
- not checked: the GPU target (nvptx) with Enzyme, the stable toolchain
  (none: nightly only, the workspace pins stable).

## 4. Answer to Q2 and recommendation

D221 stands: **option A is the production route for A1** (stable
toolchain, `no_std`, the GPU source unchanged, exact), with the stored mode
below about 1 M cells and the recompute mode above. Option B is exact and
about 8 to 16 times cheaper than A's recompute mode and fits the same
structure (only the face-local step changes), so the choice can switch
later without redoing the rest: re-check when the Enzyme backend reaches
the stable toolchain or the GPU target.

## 5. Gaps (what A1 needs)

- The chain through the other kernels: least-squares gradient, limiter
  (frozen or smooth), boundary values, Rhie–Chow flux, rows; the
  dot-product test then on the whole R.
- Periodic faces (the neighbour does not list them): `reverse` returns
  `None` for them now.
- The transposed solve (FGMRES with the transposed coupled preconditioner),
  objectives, parameters and the `adjoint.*` tools of D223.
- A0b (a solver generic over `Real`, decision Q7) is the next item of the
  plan and is not started; its first slice will be measured before it is
  sized.
